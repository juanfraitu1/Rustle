#!/usr/bin/env python3
"""One erased copy hidden in multimapping reads, simulated (2026-10-23). The user: "can we somehow only focus on one erased copy? We
could even simulate it so it's clean"; "the base case is absent copy that is hidden in multi mapping reads ... multi mapping reads are
everywhere there is similarity so pseudogenes should be trackable"; both cases below, genome-missing first.

Gorilla KB3781, a family with >= 3 copies. Every copy's transcript is its dominant fibroblast isoform (the most frequent intron chain
among the fibroblast BAM's primary reads inside the copy, ends at their medians) on both of his haplotypes: A = `_pri`, B = the blocks
lifted to the other haplotype through the frozen truth's asm5 alignments (isocon_sim.lift_blocks). Reads come from those transcripts
with isocon_sim's model (5' loss Exp(mean 100), 3' loss Exp(mean 10); 0.5% substitutions, 0.25% insertions, 0.25% deletions), so every
read's origin is known. The reference is a mini reference: the family's loci +- FLANK.
  - case `missing`: the erased copy E is masked (N) in the reference and absent from the copy table: the sample carries a copy the
    reference lacks; its reads can only land on the survivors.
  - case `unannotated`: E stays in the reference but not in the copy table; with --e-reads-per-allele 0 E is silent (a pseudogene).
Then the O3 port (o3port.py run, v5 defaults; no IsoCon side, there is no stage run) on the survivors, and where the family's reads
land outside the copy table.
    o3sim_one.py pick --out DIR
        the simplest family (fewest copies, fresh before development, then name; >= 3 copies; refabsent's 34 excluded) and erased
        copy (Amendment 7's order: last by chrom, start first) with every copy's isoform on both haplotypes, whose transcript is
        >= 0.97 identical over >= 80% of its length to >= 2 survivors' and for which >= half of 40 simulated E reads are
        multimappers in case `missing` (a record within AS_TIE of their best on >= 2 surviving copies) -> DIR/case.json, DIR/truth.fa
    o3sim_one.py sim --out DIR --name NAME --case missing|unannotated [--s-alleles AB|A] [--e-alleles AB|A]
                     [--reads-per-allele 100] [--e-reads-per-allele 100] [--seed 1]
    o3sim_one.py report --out DIR --name NAME   -> DIR/NAME.report.json, stdout
"""
import argparse
import collections
import csv
import json
import os
import random
import re
import statistics
import subprocess
import sys

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "rna_allele"))
import isocon_sim  # noqa: E402  lift_blocks, mutate

RA = "/mnt/linuxdisk/tmp/rna_allele"
SRC = "/home/juanfra/winloci_scratch/o3_collapse/method/intervals/data/intervals.tsv"
PRI = "/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta"
FIBRO = "/mnt/linuxdisk/home/juanfraitu/fibroblasts/GCA_029281585.2_flnc_mm.bam"
FLANK = 50_000
AS_TIE = 10
PICK_READS = 40
MM_R = ["-ax", "splice:hq", "-uf", "--eqx", "-Y", "-N", "50", "-p", "0.1", "--secondary=yes", "-t", "4"]
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s):
    return s.translate(COMP)[::-1]


def tsv(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def fasta(path):
    out, name = {}, None
    for ln in open(path):
        if ln.startswith(">"):
            name = ln[1:].split()[0]
            out[name] = []
        elif name:
            out[name].append(ln.strip())
    return {k: "".join(v) for k, v in out.items()}


def dominant_isoform(bam, chrom, s, e):
    """(exon blocks, strand, reads with that chain, reads): the most frequent intron chain among primary reads inside [s, e)."""
    chains, ends, strand = collections.Counter(), collections.defaultdict(list), collections.Counter()
    n = 0
    for rd in bam.fetch(chrom, s, e):
        if rd.is_unmapped or rd.is_secondary or rd.is_supplementary or rd.reference_start < s or rd.reference_end > e:
            continue
        n += 1
        introns, pos = [], rd.reference_start
        for op, k in rd.cigartuples:
            if op == 3:
                introns.append((pos, pos + k))
            if op in (0, 2, 3, 7, 8):
                pos += k
        key = tuple(introns)
        chains[key] += 1
        ends[key].append((rd.reference_start, rd.reference_end))
        strand["-" if rd.is_reverse else "+"] += 1
    if not chains:
        return None
    key, k = max(chains.items(), key=lambda x: (x[1], len(x[0])))
    cuts = [int(statistics.median(x for x, _ in ends[key]))] + [x for iv in key for x in iv] + \
           [int(statistics.median(y for _, y in ends[key]))]
    return [(cuts[i], cuts[i + 1]) for i in range(0, len(cuts), 2)], strand.most_common(1)[0][0], k, n


def transcripts(pri, cm, chrom, blocks, strand):
    """{'A': seq, 'B': seq or None} in transcription orientation."""
    sA = "".join(pri.fetch(chrom, x, y).upper() for x, y in blocks)
    out = {"A": rc(sA) if strand == "-" else sA, "B": None}
    c = cm[chrom]
    lb = isocon_sim.lift_blocks(f"{RA}/out/chr{c['chrom']}.paf", chrom, blocks)
    if lb:
        bfa = pysam.FastaFile(f"{RA}/t/chr{c['chrom']}.fa")
        sB = "".join(bfa.fetch(c["B_name"], t0, t1).upper() if st == "+" else rc(bfa.fetch(c["B_name"], t0, t1).upper())
                     for t0, t1, st in lb)
        out["B"] = rc(sB) if strand == "-" else sB
    return out


def compare(a, b, tmp):
    """(identity, mismatches, indel events, coverage of b) of b against a (minimap2 -c --eqx -x map-ont: asm20 misses most
    transcript pairs of a few kb; the longest alignment), or None."""
    with open(f"{tmp}/t.fa", "w") as f:
        f.write(f">a\n{a}\n")
    with open(f"{tmp}/q.fa", "w") as f:
        f.write(f">b\n{b}\n")
    paf = subprocess.run(["minimap2", "-c", "--eqx", "-x", "map-ont", "-t", "1", f"{tmp}/t.fa", f"{tmp}/q.fa"],
                         capture_output=True, text=True).stdout.splitlines()
    if not paf:
        return None
    p = max((ln.split("\t") for ln in paf), key=lambda f: int(f[10]))
    ops = re.findall(r"(\d+)([=XID])", next(x[5:] for x in p[12:] if x.startswith("cg:Z:")))
    return (int(p[9]) / max(1, int(p[10])), sum(int(n) for n, o in ops if o == "X"), sum(1 for n, o in ops if o in "ID"),
            (int(p[3]) - int(p[2])) / max(1, int(p[1])))


def mini_reference(copies, erased, mask, out):
    """The family's loci +- FLANK as segments (overlapping ones merged); E's interval set to N when `mask`. Returns
    {copy: (segment, start, end)} in segment coordinates."""
    pri = pysam.FastaFile(PRI)
    iv = sorted((c["chrom"], max(0, c["start"] - FLANK), min(pri.get_reference_length(c["chrom"]), c["end"] + FLANK))
                for c in copies.values())
    segs = []
    for ch, lo, hi in iv:
        if segs and segs[-1][0] == ch and lo <= segs[-1][2]:
            segs[-1][2] = max(segs[-1][2], hi)
        else:
            segs.append([ch, lo, hi])
    where = {}
    with open(f"{out}/ref.fa", "w") as f:
        for ch, lo, hi in segs:
            name = f"{ch}_{lo}_{hi}"
            seq = list(pri.fetch(ch, lo, hi).upper())
            for cid, c in copies.items():
                if c["chrom"] == ch and lo <= c["start"] and c["end"] <= hi:
                    where[cid] = (name, c["start"] - lo, c["end"] - lo)
                    if mask and cid == erased:
                        seq[c["start"] - lo:c["end"] - lo] = "N" * (c["end"] - c["start"])
            f.write(f">{name}\n{''.join(seq)}\n")
    pysam.faidx(f"{out}/ref.fa")
    return where


def simulate(tx, plan, rng, path):
    """plan: [(transcript name, n reads)] -> FASTA of reads named <transcript>_r<i>."""
    with open(path, "w") as f:
        for t, k in plan:
            s = tx[t]
            for i in range(k):
                a5 = min(int(rng.expovariate(1 / 100)), len(s) // 2)
                a3 = min(int(rng.expovariate(1 / 10)), len(s) // 4)
                f.write(f">{t}_r{i:03d}\n{isocon_sim.mutate(s[a5:len(s) - a3], rng)}\n")


def align(out, reads, bam):
    with open(bam, "wb") as bo:
        mm = subprocess.Popen(["minimap2", *MM_R, f"{out}/ref.fa", reads], stdout=subprocess.PIPE,
                              stderr=open(bam + ".mm2.log", "w"))
        subprocess.run(["samtools", "sort", "-o", "-", "-"], stdin=mm.stdout, stdout=bo, check=True)
        mm.wait()
    pysam.index(bam)


def landing(bam_path, where, copies_in_table):
    """Per read: its best AS, and the copies of the table holding a record within AS_TIE of it; plus every record (primary or
    secondary) outside the table's copies, for the locus scan."""
    best, recs = collections.defaultdict(lambda: -10 ** 9), collections.defaultdict(list)
    for rd in pysam.AlignmentFile(bam_path).fetch(until_eof=True):
        if rd.is_unmapped or rd.is_supplementary or not rd.has_tag("AS"):
            continue
        a = rd.get_tag("AS")
        best[rd.query_name] = max(best[rd.query_name], a)
        recs[rd.query_name].append((rd.reference_name, rd.reference_start, rd.reference_end, a, rd.is_secondary))
    tied, outside = {}, []
    for n, rs in recs.items():
        on = set()
        for seg, s, e, a, sec in rs:
            hit = [c for c in copies_in_table if where[c][0] == seg and s < where[c][2] and where[c][1] < e]
            if hit and a >= best[n] - AS_TIE:
                on.update(hit)
            if not hit:
                outside.append((seg, s, e, n, sec, a >= best[n] - AS_TIE))
        tied[n] = on
    return tied, outside, best


def pick(a):
    os.makedirs(a.out, exist_ok=True)
    rows = tsv(SRC)
    byfam = collections.defaultdict(list)
    for r in rows:
        byfam[r["fam"]].append(r)
    dev = {r["family_id"] for r in tsv(f"{RA}/a12/A12.copies.tsv")}
    bonly = {r["family"] for r in tsv(f"{RA}/refabsent/bonly.tsv")}
    # the simplest family first: fewest copies, then fresh before development, then name
    order = sorted((f for f, v in byfam.items() if len(v) >= 3 and f not in bonly), key=lambda f: (len(byfam[f]), f in dev, f))
    cm = {r["pri"]: r for r in tsv(f"{RA}/chrmap.tsv")}
    pri, bam = pysam.FastaFile(PRI), pysam.AlignmentFile(FIBRO)
    for fam in order:
        cps = sorted(byfam[fam], key=lambda r: (r["chrom"], int(r["clean_start"])))
        copies, tx = {}, {}
        for r in cps:
            iso = dominant_isoform(bam, r["chrom"], int(r["clean_start"]), int(r["clean_end"]))
            if not iso:
                break
            blocks, strand, k, n = iso
            t = transcripts(pri, cm, r["chrom"], blocks, strand)
            if not t["B"]:
                break
            copies[r["gene"]] = dict(chrom=r["chrom"], start=int(r["clean_start"]), end=int(r["clean_end"]), blocks=blocks,
                                     strand=strand, chain_reads=k, reads=n)
            tx[f"{r['gene']}|A"], tx[f"{r['gene']}|B"] = t["A"], t["B"]
        if len(copies) < len(cps):
            print(f"{fam}: a copy has no fibroblast isoform or no B allele; skipped", flush=True)
            continue
        for erased in [r["gene"] for r in reversed(cps)]:
            surv = [c for c in copies if c != erased]
            ident = {c: compare(tx[f"{c}|A"], tx[f"{erased}|A"], a.out) for c in surv}
            close = [c for c, v in ident.items() if v and v[0] >= 0.97 and v[3] >= 0.8]
            if len(close) < 2:
                print(f"{fam}: erase {erased}: {len(close)} survivors at >= 0.97 identity over >= 80% of its transcript", flush=True)
                continue
            where = mini_reference(copies, erased, True, a.out)
            simulate(tx, [(f"{erased}|A", PICK_READS // 2), (f"{erased}|B", PICK_READS // 2)], random.Random(1),
                     f"{a.out}/pick.fa")
            align(a.out, f"{a.out}/pick.fa", f"{a.out}/pick.bam")
            tied, _, _ = landing(f"{a.out}/pick.bam", where, surv)
            frac = sum(1 for v in tied.values() if len(v) >= 2) / PICK_READS
            print(f"{fam}: erase {erased}: {frac:.2f} of its reads tie on >= 2 survivors", flush=True)
            if frac >= 0.5:
                alle = {c: compare(tx[f"{c}|A"], tx[f"{c}|B"], a.out) for c in copies}
                case = dict(family=fam, erased=erased, survivors=surv, tie_fraction=frac, copies=copies,
                            erased_vs_survivor={c: v and dict(identity=round(v[0], 5), mismatches=v[1], indels=v[2], coverage=round(v[3], 3))
                                                for c, v in ident.items()},
                            alleles={c: v and dict(identity=round(v[0], 5), mismatches=v[1], indels=v[2]) for c, v in alle.items()})
                json.dump(case, open(f"{a.out}/case.json", "w"), indent=1)
                with open(f"{a.out}/truth.fa", "w") as f:
                    for k, s in tx.items():
                        f.write(f">{k}\n{s}\n")
                print(json.dumps({k: case[k] for k in ("family", "erased", "survivors", "tie_fraction", "erased_vs_survivor", "alleles")}))
                return
    print("no family meets the rule")


def sim(a):
    W = a.out
    case = json.load(open(f"{W}/case.json"))
    tx = fasta(f"{W}/truth.fa")
    E, surv, copies, fam = case["erased"], case["survivors"], case["copies"], case["family"]
    if a.s_het:   # synthetic heterozygosity: a survivor whose two real alleles differ at < 2 positions gets a B allele = A with
        rng_h = random.Random(a.seed + 1)  # s_het random substitutions (the certificate needs >= 2 differences to see two alleles)
        for c in surv:
            al = case["alleles"].get(c)
            if al is None or al["mismatches"] + al["indels"] < 2:
                sq = list(tx[f"{c}|A"])
                for i in rng_h.sample(range(len(sq)), a.s_het):
                    sq[i] = rng_h.choice([b for b in "ACGT" if b != sq[i]])
                tx[f"{c}|B"] = "".join(sq)
    plan = []
    for c in copies:
        alleles, n = (a.e_alleles, a.e_reads_per_allele) if c == E else (a.s_alleles, a.reads_per_allele)
        if n == 0:
            continue
        plan += [(f"{c}|{h}", n if len(alleles) == 2 else 2 * n) for h in alleles]
    os.makedirs(f"{W}/{a.name}", exist_ok=True)
    D = f"{W}/{a.name}"
    simulate(tx, plan, random.Random(a.seed), f"{D}/reads.fa")
    where = mini_reference(copies, E, a.case == "missing", D)
    align(D, f"{D}/reads.fa", f"{D}/reads.bam")
    with open(f"{D}/copies.tsv", "w") as f:
        f.write("family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\n")
        for k, c in enumerate(surv):
            seg, s0, e0 = where[c]
            f.write(f"{fam}\t{k}\t{c}\t{seg}\t{s0}\t{e0}\t1\t+\t0\t{s0}-{e0}\n")
    names = [ln[1:].strip() for ln in open(f"{D}/reads.fa") if ln.startswith(">")]
    tied, outside, best = landing(f"{D}/reads.bam", where, surv)
    with open(f"{D}/net.tsv", "w") as f:
        f.write("family\tread\tused\trole\tvia\tclass\tin_reads_tsv\n")
        for n in names:
            f.write(f"{fam}\t{n}\t1\t{'D' if n.startswith(E + '|') else 'S'}\t{'A' if tied.get(n) else 'B'}\tNA\t0\n")
    subprocess.run([sys.executable, f"{HERE}/o3port.py", "run", "--bam", f"{D}/reads.bam", "--copies", f"{D}/copies.tsv",
                    "--genome", f"{D}/ref.fa", "--net", f"{D}/net.tsv", "--out", f"{D}/port", "--dump", f"{D}/port.dump.jsonl"],
                   check=True, stderr=open(f"{D}/port.log", "w"))
    # where the family's reads land outside the copy table: merged loci of records (primary or secondary)
    loci = []
    for seg, s, e, n, sec, near in sorted(outside):
        if loci and loci[-1]["seg"] == seg and s <= loci[-1]["end"]:
            L = loci[-1]
            L["end"] = max(L["end"], e)
        else:
            L = dict(seg=seg, start=s, end=e, reads=set(), secondary=set(), near_best=set())
            loci.append(L)
        L["reads"].add(n)
        if sec:
            L["secondary"].add(n)
        if near:
            L["near_best"].add(n)
    eseg, es, ee = where[E]
    meta = dict(case=a.case, s_alleles=a.s_alleles, e_alleles=a.e_alleles, s_het=a.s_het, reads_per_allele=a.reads_per_allele,
                e_reads_per_allele=a.e_reads_per_allele, seed=a.seed, n_reads=len(names),
                loci=[dict(seg=L["seg"], start=L["start"], end=L["end"], reads=len(L["reads"]), secondary_only=len(L["secondary"]),
                           near_best=len(L["near_best"]), on_E=L["seg"] == eseg and L["start"] < ee and es < L["end"],
                           origin=dict(collections.Counter(n.split("|")[0] for n in L["reads"])))
                      for L in loci if len(L["reads"]) >= 2])
    json.dump(meta, open(f"{D}/meta.json", "w"), indent=1)
    print(f"{a.name}: {a.case}; reads {len(names)}; " + open(f"{D}/port.log").read().strip().splitlines()[-1])


def report(a):
    W, D = a.out, f"{a.out}/{a.name}"
    case, meta = json.load(open(f"{W}/case.json")), json.load(open(f"{D}/meta.json"))
    E, surv = case["erased"], case["survivors"]
    d = json.loads(open(f"{D}/port.dump.jsonl").readline())
    origin = lambda r: r.rsplit("_r", 1)[0]
    names = [ln[1:].strip() for ln in open(f"{D}/reads.fa") if ln.startswith(">")]
    # E reads: multimappers? (records within AS_TIE on >= 2 survivors), and their piles
    piles = {c: set(v["pile"]) for c, v in d["copies"].items()}
    e_reads = [n for n in names if n.startswith(E + "|")]
    e_multi = sum(1 for n in e_reads if sum(1 for c in piles if n in piles[c]) >= 2)
    e_any = sum(1 for n in e_reads if any(n in p for p in piles.values()))
    rows = []
    for m in d["merged"]:
        c = collections.Counter(origin(r) for r in m["reads"])
        rows.append(dict(reads=len(m["reads"]), piles=m["piles"], origin=dict(sorted(c.items())), fate=m["fate"], slot=m["slot"],
                         certified=m["certified"], private=len(m.get("private", []))))
    tx = fasta(f"{W}/truth.fa")
    contigs = fasta(f"{D}/port.contigs.fa") if os.path.exists(f"{D}/port.contigs.fa") else {}
    best = {}
    for c, s in contigs.items():
        hits = {t: compare(seq, s, D) for t, seq in tx.items()}
        hits = {t: h for t, h in hits.items() if h}
        if hits:
            t = max(hits, key=lambda t: (hits[t][0], -hits[t][1]))
            best[c] = dict(truth=t, identity=round(hits[t][0], 5), mismatches=hits[t][1], indels=hits[t][2], length=len(s))
    cert_e = sum(v for r in rows if r["certified"] for k, v in r["origin"].items() if k.startswith(E + "|"))
    cert_s = sum(v for r in rows if r["certified"] for k, v in r["origin"].items() if not k.startswith(E + "|"))
    rep = dict(name=a.name, meta=meta, E=E, E_reads=len(e_reads), E_reads_in_a_pile=e_any, E_reads_in_2plus_piles=e_multi,
               haplotypes=rows, contigs=best, E_reads_certified=cert_e, other_reads_certified=cert_s)
    json.dump(rep, open(f"{D}/report.json", "w"), indent=1)
    print(f"== {a.name} ({meta['case']}; survivors {meta['s_alleles']} x {meta['reads_per_allele']}/allele"
          + (f" (synthetic B alleles: {meta['s_het']} substitutions where the real ones differ at < 2)" if meta.get("s_het") else "")
          + f", erased {E} "
          f"{meta['e_alleles']} x {meta['e_reads_per_allele']}/allele)")
    print(f"   E reads {len(e_reads)}: in a survivor pile {e_any}, in >= 2 piles (multimappers) {e_multi}")
    for r in rows:
        print(f"   haplotype {r['reads']:4d} reads {r['origin']} piles {r['piles']}: {r['fate']}, slot {r['slot']}, "
              f"certified {r['certified']}" + (f", private alleles {r['private']}" if r["certified"] else ""))
    for c, b in best.items():
        print(f"   contig {c} ({b['length']} bp) = {b['truth']} at identity {b['identity']} ({b['mismatches']} mm, {b['indels']} indels)")
    print(f"   certified reads: from E {cert_e}, from survivors {cert_s}")
    for L in meta["loci"]:
        print(f"   reads landing outside the copy table: {L['seg']}:{L['start']}-{L['end']} reads {L['reads']} (secondary "
              f"{L['secondary_only']}, within {AS_TIE} AS of best {L['near_best']}), on E {L['on_E']}, from {L['origin']}")


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("pick")
    p.add_argument("--out", required=True)
    s = sub.add_parser("sim")
    s.add_argument("--out", required=True)
    s.add_argument("--name", required=True)
    s.add_argument("--case", required=True, choices=["missing", "unannotated"])
    s.add_argument("--s-alleles", default="AB", choices=["AB", "A"])
    s.add_argument("--e-alleles", default="AB", choices=["AB", "A"])
    s.add_argument("--reads-per-allele", type=int, default=100)
    s.add_argument("--e-reads-per-allele", type=int, default=100)
    s.add_argument("--seed", type=int, default=1)
    s.add_argument("--s-het", type=int, default=0, help="synthetic: survivors with < 2 allele differences get a B allele with "
                   "this many random substitutions")
    r = sub.add_parser("report")
    r.add_argument("--out", required=True)
    r.add_argument("--name", required=True)
    a = ap.parse_args(argv)
    {"pick": pick, "sim": sim, "report": report}[a.cmd](a)


if __name__ == "__main__":
    main(sys.argv[1:])
