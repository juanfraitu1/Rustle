#!/usr/bin/env python3
"""A held-out DONOR simulation (2026-10-20; the user chose "Donor simulation first" for the second round's held-out proof):
reads simulated from a never-used donor's own phased assembly (HG00733, HPRC year 1), so the donor-vs-CHM13 differences of
every copy are real, and every read's truth is exact: the CHM13 copy whose donor ortholog produced it.

    donorsim.py queries FAMILIES.json SET FAMS.json GENOME.fa OUT.fa
        each copy's CHM13 span +- FLANK, named copy|chrom|lo
    donorsim.py orthologs FAMILIES.json SET FAMS.json QUERIES.fa OUT.tsv HAP=PAF ...
        per copy and haplotype (minimap2 -x asm20 -c of the queries to the haplotype): the longest primary alignment with
        MAPQ >= MIN_MAPQ through which every exon end of every isoform lifts (the ortholog, its isoforms lifted); two copies
        whose orthologs overlap on a haplotype are both dropped (ambiguous)
    donorsim.py build FAMILIES.json SET FAMS.json ORTHOLOGS.tsv PROFILE.json OUT_DIR HAP=FASTA ... [--per-copy N] [--seed S]
        OUT_DIR/reads.fq and truth.tsv (hg002truth's columns, class `unique`, truth = the CHM13 copy; + hap, isoform): per
        expressed copy PER_COPY x its expression multiplier reads (realworld's: lognormal, non-coding members at the profile's
        measured ratio), split evenly over its haplotypes' orthologs and over isoforms by weight (main 2.0, others 0.5); each
        read the donor transcript's own sequence with the profile's retained introns, truncation, end jitter and HiFi errors.
"""
import collections
import csv
import importlib.util
import json
import os
import random
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FLANK, MIN_MAPQ = 20000, 20
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _sibling(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}_ds", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


RW = _sibling("realworld")


def lift_exons(row, exons):
    """Exons (query coordinates, half-open) lifted onto the target through a PAF row's alignment, sorted; None when an exon
    end falls outside it or the exons do not stay ordered. A '-' row aligns the query's reverse complement."""
    pos = [x for a, b in exons for x in (a, b - 1)]
    if row["strand"] == "+":
        m = RW.lift_positions(row["cg"], row["qs"], row["ts"], pos)
        new = [(m[a], m[b - 1] + 1) if m[a] is not None and m[b - 1] is not None else None for a, b in exons]
    else:
        rc = {p: row["qlen"] - 1 - p for p in pos}
        m = RW.lift_positions(row["cg"], row["qlen"] - row["qe"], row["ts"], list(rc.values()))
        new = [(m[rc[b - 1]], m[rc[a]] + 1) if m[rc[a]] is not None and m[rc[b - 1]] is not None else None for a, b in exons]
    if any(e is None or e[0] >= e[1] for e in new):
        return None
    new.sort()
    return new if all(x[1] <= y[0] for x, y in zip(new, new[1:])) else None


def pick_ortholog(rows, isoforms):
    """The ortholog among a query's PAF rows: the longest primary (tp P) with MAPQ >= MIN_MAPQ, when every isoform
    ({t: exons in query coordinates}) lifts through it -> {tname, strand, start, end, isoforms}, else None."""
    prim = [r for r in rows if r["tp"] == "P"]
    if not prim:
        return None
    best = max(prim, key=lambda r: r["te"] - r["ts"])
    if best["mapq"] < MIN_MAPQ:
        return None
    lifted = {t: lift_exons(best, ex) for t, ex in isoforms.items()}
    if any(v is None for v in lifted.values()):
        return None
    lo = min(e[0][0] for e in lifted.values())
    hi = max(e[-1][1] for e in lifted.values())
    return dict(tname=best["tname"], strand=best["strand"], start=lo, end=hi, isoforms=lifted, mapq=best["mapq"])


def drop_ambiguous(orth):
    """{(copy, hap): ortholog} less every pair whose orthologs overlap on the same haplotype contig."""
    bad = set()
    by = collections.defaultdict(list)
    for (c, h), o in orth.items():
        by[(h, o["tname"])].append((o["start"], o["end"], c))
    for (h, _), lst in by.items():
        lst.sort()
        for i, (s, e, c) in enumerate(lst):
            for s2, e2, c2 in lst[i + 1:]:
                if s2 >= e:
                    break
                bad |= {(c, h), (c2, h)}
    return {k: v for k, v in orth.items() if k not in bad}


def _family_copies(fam_p, set_name, sel_p):
    sel = set(json.load(open(sel_p)))
    out = []
    for f in json.load(open(fam_p))["sets"][set_name]:
        if f["family"] in sel:
            for k, c in enumerate(f["copies"]):
                out.append((f"{f['family']}.{k}", f["family"], c))
    return out


def cmd_queries(fam_p, set_name, sel_p, genome_p, out_p):
    import pysam
    fa = pysam.FastaFile(genome_p)
    n = 0
    with open(out_p, "w") as fh:
        for cid, f, c in _family_copies(fam_p, set_name, sel_p):
            lo = max(0, c["start"] - FLANK)
            hi = min(fa.get_reference_length(c["chrom"]), c["end"] + FLANK)
            fh.write(f">{cid}|{c['chrom']}|{lo}\n{fa.fetch(c['chrom'], lo, hi).upper()}\n")
            n += 1
    print(f"[donorsim] {n} flanked copy queries -> {out_p}", flush=True)
    return 0


def read_paf(p):
    rows = collections.defaultdict(list)
    for ln in open(p):
        x = ln.rstrip("\n").split("\t")
        tags = {t[:5]: t[5:] for t in x[12:]}
        if "cg:Z:" not in tags:
            continue
        rows[x[0]].append(dict(qlen=int(x[1]), qs=int(x[2]), qe=int(x[3]), strand=x[4], tname=x[5], ts=int(x[7]),
                               te=int(x[8]), mapq=int(x[11]), tp=tags.get("tp:A:", "P"), cg=tags["cg:Z:"]))
    return rows


def cmd_orthologs(fam_p, set_name, sel_p, queries_p, out_p, haps):
    lo_of = {}
    for ln in open(queries_p):
        if ln.startswith(">"):
            cid, _, lo = ln[1:].strip().split("|")
            lo_of[cid] = int(lo)
    orth, n = {}, collections.Counter()
    copies = _family_copies(fam_p, set_name, sel_p)
    for h, paf_p in haps.items():
        rows = read_paf(paf_p)
        for cid, f, c in copies:
            q = f"{cid}|{c['chrom']}|{lo_of[cid]}"
            iso = {t: [(a - lo_of[cid], b - lo_of[cid]) for a, b in ex] for t, ex in c["transcripts"].items()}
            o = pick_ortholog(rows.get(q, []), iso)
            n[f"{h}: {'ortholog' if o else 'none'}"] += 1
            if o:
                orth[(cid, h)] = o
    kept = drop_ambiguous(orth)
    n["ambiguous (dropped)"] = len(orth) - len(kept)
    with open(out_p, "w") as fh:
        fh.write("copy\thap\tcontig\tstart\tend\tstrand\tmapq\tisoforms\n")
        for (cid, h), o in sorted(kept.items()):
            fh.write(f"{cid}\t{h}\t{o['tname']}\t{o['start']}\t{o['end']}\t{o['strand']}\t{o['mapq']}\t"
                     f"{json.dumps(o['isoforms'], separators=(',', ':'))}\n")
    print(f"[donorsim] orthologs -> {out_p}: " + ", ".join(f"{k} {v}" for k, v in sorted(n.items()))
          + f"; copies with >= 1 ortholog {len({c for c, _ in kept})} of {len(copies)}", flush=True)
    return 0


def cmd_build(fam_p, set_name, sel_p, orth_p, prof_p, out_dir, haps, per_copy=100, seed=2026):
    import pysam
    sys.path.insert(0, os.path.dirname(HERE))
    from sim import stable_seed, write_fastq
    prof = json.load(open(prof_p))
    fas = {h: pysam.FastaFile(p) for h, p in haps.items()}
    orth = collections.defaultdict(dict)
    for r in csv.DictReader(open(orth_p), delimiter="\t"):
        orth[r["copy"]][r["hap"]] = (r["contig"], r["strand"], {t: [tuple(e) for e in ex] for t, ex in json.loads(r["isoforms"]).items()})
    os.makedirs(out_dir, exist_ok=True)
    reads, n = [], collections.Counter()
    for cid, f, c in _family_copies(fam_p, set_name, sel_p):
        if cid not in orth:
            n["copies without an ortholog"] += 1
            continue
        mult = RW.expr_mult(seed, cid, prof.get("expr_log_sd", 0.0)) * RW.paralog_mult(seed, cid, c.get("biotype", "protein_coding"), prof)
        if mult == 0.0:
            n["silent copies"] += 1
            continue
        w = {t: (RW.MAIN_WEIGHT if t == c["main"] else RW.OTHER_WEIGHT) for t in c["transcripts"]}
        tot = sum(w.values())
        hs = sorted(orth[cid])
        for h in hs:
            contig, ostrand, iso = orth[cid][h]
            strand = c["strand"] if ostrand == "+" else ("-" if c["strand"] == "+" else "+")
            fetch = (lambda fa_, ch: (lambda a_, b_: fa_.fetch(ch, a_, b_)))(fas[h], contig)
            cache = {}
            for t, ex in iso.items():
                k = int(round(per_copy * mult * w[t] / tot / len(hs)))
                for i in range(k):
                    r = random.Random(stable_seed(str(seed), cid, h, t, str(i)))
                    e2 = [tuple(x) for x in ex]
                    if len(e2) > 1:
                        keep = {j for j in range(len(e2) - 1) if r.random() < prof["kept_intron_per_intron"]}
                        e2 = RW.retain(e2, keep) if keep else e2
                    key = tuple(e2)
                    if key not in cache:
                        cache[key] = RW.splice(fetch, e2, strand, {})
                    seq = RW.draw_truncation(cache[key], prof["end_pairs"], r)
                    lo, hi = r.randint(0, 30), len(seq) - r.randint(0, 30)
                    if hi - lo < 100:
                        lo, hi = 0, len(seq)
                    b = seq[lo:hi]
                    b = RW.apply_edits(b, RW.draw_edits(b, r, prof))
                    reads.append((cid, f, h, t, (b, "~" * len(b))))
                    n[f"reads from {h}"] += 1
    names = RW.anon_names([(x[0], f"{x[2]}|{x[3]}") for x in reads], seed)
    order = sorted(range(len(reads)), key=lambda i: names[i])
    with open(f"{out_dir}/reads.fq", "w") as fq, open(f"{out_dir}/truth.tsv", "w") as tr:
        tr.write("read\tfamily\ttied_copies\tclass\ttruth\tdonor_nm\tdonor_margin\thap\tisoform\n")
        for i in order:
            cid, f, h, t, rd = reads[i]
            nm = "d" + names[i][1:]
            write_fastq(fq, nm, rd)
            tr.write(f"{nm}\t{f}\t\tunique\t{cid}\t\t\t{h}\t{t}\n")
    json.dump(dict(n), open(f"{out_dir}/build.json", "w"), indent=1)
    print(f"[donorsim] {out_dir}: {len(reads)} reads; " + ", ".join(f"{k} {v}" for k, v in sorted(n.items())), flush=True)
    return 0


def main(argv):
    haps = dict(a.split("=", 1) for a in argv if "=" in a and not a.startswith("--"))
    rest = [a for a in argv if not ("=" in a and not a.startswith("--"))]
    opt = {}
    for k in ("--per-copy", "--seed"):
        if k in rest:
            i = rest.index(k)
            opt[k], rest = int(rest[i + 1]), rest[:i] + rest[i + 2:]
    if rest[:1] == ["queries"] and len(rest) == 6:
        return cmd_queries(*rest[1:])
    if rest[:1] == ["orthologs"] and len(rest) == 6:
        return cmd_orthologs(*rest[1:], haps)
    if rest[:1] == ["build"] and len(rest) == 7:
        return cmd_build(*rest[1:], haps, per_copy=opt.get("--per-copy", 100), seed=opt.get("--seed", 2026))
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
