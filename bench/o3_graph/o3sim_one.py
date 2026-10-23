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


def mismatches(a, b, tmp):
    """[(position in b, position in a)] of b's substitutions against a (minimap2 map-ont --eqx, the longest + strand alignment)."""
    with open(f"{tmp}/t.fa", "w") as f:
        f.write(f">a\n{a}\n")
    with open(f"{tmp}/q.fa", "w") as f:
        f.write(f">b\n{b}\n")
    paf = [ln.split("\t") for ln in subprocess.run(["minimap2", "-c", "--eqx", "-x", "map-ont", "-t", "1", f"{tmp}/t.fa", f"{tmp}/q.fa"],
                                                    capture_output=True, text=True).stdout.splitlines()]
    paf = [p for p in paf if p[4] == "+"]
    if not paf:
        return []
    p = max(paf, key=lambda f: int(f[10]))
    q, t, out = int(p[2]), int(p[7]), []
    for n, o in re.findall(r"(\d+)([=XID])", next(x[5:] for x in p[12:] if x.startswith("cg:Z:"))):
        n = int(n)
        if o == "X":
            out += [(q + k, t + k) for k in range(n)]
        if o in "=X":
            q, t = q + n, t + n
        elif o == "I":
            q += n
        else:
            t += n
    return out


def configure(case, tx, cfg, tmp):
    """The sample's transcripts for one configuration (cfg keys: s_het, e_frac, indiv_hom, indiv_het, seed): synthetic survivor
    heterozygosity, the erased copy moved toward its closest survivor, and the sample's own variants (another individual than
    the reference)."""
    tx = dict(tx)
    E, surv = case["erased"], case["survivors"]
    if cfg["s_het"]:   # a survivor whose two real alleles differ at < 2 positions gets a B allele = A + s_het substitutions
        rng_h = random.Random(cfg["seed"] + 1)
        for c in surv:
            al = case["alleles"].get(c)
            if al is None or al["mismatches"] + al["indels"] < 2:
                sq = list(tx[f"{c}|A"])
                for i in rng_h.sample(range(len(sq)), cfg["s_het"]):
                    sq[i] = rng_h.choice([b for b in "ACGT" if b != sq[i]])
                tx[f"{c}|B"] = "".join(sq)
    if cfg["e_frac"] < 1:   # each substitution of E against its closest survivor T kept with probability e_frac (else T's base)
        T = max(surv, key=lambda c: (case["erased_vs_survivor"][c] or {}).get("identity", 0))
        rng_f = random.Random(cfg["seed"] + 2)
        keep = {}
        for h in ("A", "B"):
            sq = list(tx[f"{E}|{h}"])
            for qp, tp in mismatches(tx[f"{T}|A"], tx[f"{E}|{h}"], tmp):
                if keep.setdefault(tp, rng_f.random() < cfg["e_frac"]) is False:
                    sq[qp] = tx[f"{T}|A"][tp]
            tx[f"{E}|{h}"] = "".join(sq)
    if cfg["indiv_hom"] or cfg["indiv_het"]:   # the sample's own variants: per copy, sites at rate indiv_hom carry one alternative
        rng_i = random.Random(cfg["seed"] + 3)  # base on both alleles, sites at rate indiv_het on one allele (positions scaled
        for c in case["copies"]:                # between the two alleles' lengths)
            A, B = list(tx[f"{c}|A"]), list(tx[f"{c}|B"])
            for i in range(len(A)):
                r = rng_i.random()
                if r >= cfg["indiv_hom"] + cfg["indiv_het"]:
                    continue
                j = min(len(B) - 1, round(i * len(B) / len(A)))
                alt = rng_i.choice([b for b in "ACGT" if b != A[i]])
                if r < cfg["indiv_hom"]:
                    A[i] = alt
                    B[j] = alt if B[j] != alt else B[j]
                elif rng_i.random() < 0.5:
                    A[i] = alt
                else:
                    B[j] = rng_i.choice([b for b in "ACGT" if b != B[j]])
            tx[f"{c}|A"], tx[f"{c}|B"] = "".join(A), "".join(B)
    return tx


DEFAULT = dict(case="missing", s_alleles="AB", e_alleles="AB", reads_per_allele=100, e_reads_per_allele=100, seed=1, s_het=0,
               e_frac=1.0, indiv_hom=0.0, indiv_het=0.0)


def run_config(W, name, cfg):
    """Simulate, align and run the port for one configuration -> W/name/ (reads, reference, BAM, copy table, net, port products,
    meta.json)."""
    case = json.load(open(f"{W}/case.json"))
    E, surv, copies, fam = case["erased"], case["survivors"], case["copies"], case["family"]
    D = f"{W}/{name}"
    os.makedirs(D, exist_ok=True)
    tx = configure(case, fasta(f"{W}/truth.fa"), cfg, D)
    with open(f"{D}/sample.fa", "w") as f:           # the sample's transcripts as simulated (the truth for this run)
        for k, v in tx.items():
            f.write(f">{k}\n{v}\n")
    plan = []
    for c in copies:
        alleles, n = (cfg["e_alleles"], cfg["e_reads_per_allele"]) if c == E else (cfg["s_alleles"], cfg["reads_per_allele"])
        if n:
            plan += [(f"{c}|{h}", n if len(alleles) == 2 else 2 * n) for h in alleles]
    simulate(tx, plan, random.Random(cfg["seed"]), f"{D}/reads.fa")
    where = mini_reference(copies, E, cfg["case"] == "missing", D)
    align(D, f"{D}/reads.fa", f"{D}/reads.bam")
    table = surv if cfg["case"] != "control" else list(copies)        # control: nothing erased, every copy in the table
    with open(f"{D}/copies.tsv", "w") as f:
        f.write("family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\n")
        for k, c in enumerate(table):
            seg, s0, e0 = where[c]
            f.write(f"{fam}\t{k}\t{c}\t{seg}\t{s0}\t{e0}\t1\t+\t0\t{s0}-{e0}\n")
    names = [ln[1:].strip() for ln in open(f"{D}/reads.fa") if ln.startswith(">")]
    tied, _, _ = landing(f"{D}/reads.bam", where, table)
    with open(f"{D}/net.tsv", "w") as f:
        f.write("family\tread\tused\trole\tvia\tclass\tin_reads_tsv\n")
        for n in names:
            f.write(f"{fam}\t{n}\t1\t{'D' if n.startswith(E + '|') else 'S'}\t{'A' if tied.get(n) else 'B'}\tNA\t0\n")
    subprocess.run([sys.executable, f"{HERE}/o3port.py", "run", "--bam", f"{D}/reads.bam", "--copies", f"{D}/copies.tsv",
                    "--genome", f"{D}/ref.fa", "--net", f"{D}/net.tsv", "--out", f"{D}/port", "--dump", f"{D}/port.dump.jsonl"],
                   check=True, stderr=open(f"{D}/port.log", "w"))
    json.dump(dict(cfg, n_reads=len(names), E_where=where[E]), open(f"{D}/meta.json", "w"), indent=1)


def summary(W, name):
    """One row: the erased copy's reads (in a pile, multimappers), its print (the largest haplotype holding >= 90% E reads) and
    what became of it (certified / look harder / an allele slot / below the floor ...), false calls (certified or look-harder
    haplotypes with < half E reads), the certified or look-harder contig against the sample's own E transcripts, and the
    landing loci on E."""
    case, meta = json.load(open(f"{W}/case.json")), json.load(open(f"{W}/{name}/meta.json"))
    D, E = f"{W}/{name}", case["erased"]
    d = json.loads(open(f"{D}/port.dump.jsonl").readline())
    names = [ln[1:].strip() for ln in open(f"{D}/reads.fa") if ln.startswith(">")]
    e_reads = {n for n in names if n.startswith(E + "|")}
    piles = {c: set(v["pile"]) for c, v in d["copies"].items()}
    in_pile = sum(1 for n in e_reads if any(n in p for p in piles.values()))
    multi = sum(1 for n in e_reads if sum(1 for p in piles.values() if n in p) >= 2)
    is_e = lambda m: sum(1 for r in m["reads"] if r in e_reads)
    prints = [m for m in d["merged"] if is_e(m) >= 0.9 * len(m["reads"])]
    pr = max(prints, key=lambda m: len(m["reads"]), default=None)
    fate = "none"
    if pr:
        fate = ("certified" if pr["certified"] else "look harder" if pr.get("look_harder") else
                f"slot {pr['slot']}" if pr["slot"] else pr["fate"])
    false_c = sum(1 for m in d["merged"] if m["certified"] and 2 * is_e(m) < len(m["reads"]))
    false_l = sum(1 for m in d["merged"] if m.get("look_harder") and 2 * is_e(m) < len(m["reads"]))
    sample = fasta(f"{D}/sample.fa")
    ident = None
    for path in (f"{D}/port.contigs.fa", f"{D}/port.look_harder.fa"):
        for c, sq in (fasta(path) if os.path.exists(path) else {}).items():
            hits = [compare(sample[f"{E}|{h}"], sq, D) for h in ("A", "B")]
            hits = [h for h in hits if h]
            if hits:
                ident = max(ident or 0, max(h[0] for h in hits))
    seg, es, ee = meta["E_where"]
    loci_e = 0
    if os.path.exists(f"{D}/port.loci.tsv"):
        for r in tsv(f"{D}/port.loci.tsv"):
            ch, se = r["locus"].rsplit(":", 1)
            s0, e0 = (int(x) for x in se.split("-"))
            if ch == seg and s0 < ee and es < e0:
                loci_e = max(loci_e, int(r["n_reads"]))
    return dict(name=name, E_reads=len(e_reads), E_in_pile=in_pile, E_multi=multi, print_reads=len(pr["reads"]) if pr else 0,
                print_E=is_e(pr) if pr else 0, fate=fate, false_certified=false_c, false_look=false_l, contig_identity=ident,
                locus_on_E=loci_e)


def sim(a):
    cfg = dict(DEFAULT, case=a.case, s_alleles=a.s_alleles, e_alleles=a.e_alleles, reads_per_allele=a.reads_per_allele,
               e_reads_per_allele=a.e_reads_per_allele, seed=a.seed, s_het=a.s_het, e_frac=a.e_frac, indiv_hom=a.indiv_hom,
               indiv_het=a.indiv_het)
    run_config(a.out, a.name, cfg)
    print(json.dumps(summary(a.out, a.name)))


def sweep(a):
    """The registered grid (2026-10-23): expression of the erased copy, its divergence from its closest survivor, and a nothing-
    erased control, each with the reference's own individual (KB3781) and with another individual (the sample's own variants at
    a human-like rate: 0.001 on both alleles + 0.001 on one) -> OUT/sweep.tsv"""
    other = dict(indiv_hom=0.001, indiv_het=0.001)
    grid = []
    for ind_name, ind in (("same", {}), ("other", other)):
        for het_name, het in (("survivors as KB3781", 0), ("survivors heterozygous", 6)):
            for n in (100, 30, 10, 5, 3):
                grid.append((f"expr{n}_{het}_{ind_name}", dict(DEFAULT, e_reads_per_allele=n, s_het=het, **ind),
                             f"expression {n}/allele; {het_name}; {ind_name} individual"))
            for f in (0.5, 0.25, 0.1, 0.05):
                grid.append((f"div{f}_{het}_{ind_name}", dict(DEFAULT, e_frac=f, s_het=het, **ind),
                             f"E toward its closest survivor (frac {f}); {het_name}; {ind_name} individual"))
            grid.append((f"control_{het}_{ind_name}", dict(DEFAULT, case="control", s_het=het, **ind),
                         f"control (nothing erased); {het_name}; {ind_name} individual"))
    for name, cfg, label in grid:
        if not re.search(a.select, name):
            continue
        run_config(a.out, name, cfg)
        r = dict(summary(a.out, name), label=label)
        json.dump(r, open(f"{a.out}/{name}/summary.json", "w"))
        print("\t".join(str(r[k]) for k in ("name", "E_reads", "E_in_pile", "E_multi", "print_reads", "print_E", "fate",
                                             "false_certified", "false_look", "contig_identity", "locus_on_E")), flush=True)
    keys = ["name", "label", "E_reads", "E_in_pile", "E_multi", "print_reads", "print_E", "fate", "false_certified", "false_look",
            "contig_identity", "locus_on_E"]
    with open(f"{a.out}/sweep.tsv", "w") as f:     # every grid run done so far, in grid order
        f.write("\t".join(keys) + "\n")
        for name, _, _ in grid:
            if os.path.exists(f"{a.out}/{name}/summary.json"):
                r = json.load(open(f"{a.out}/{name}/summary.json"))
                f.write("\t".join(str(r[k]) for k in keys) + "\n")


def report(a):
    """The run's summary, then every haplotype with its reads by origin and its fate."""
    D = f"{a.out}/{a.name}"
    print(json.dumps(summary(a.out, a.name)))
    d = json.loads(open(f"{D}/port.dump.jsonl").readline())
    for m in d["merged"]:
        c = collections.Counter(r.rsplit("_r", 1)[0] for r in m["reads"])
        print(f"   haplotype {len(m['reads']):4d} reads {dict(sorted(c.items()))} piles {m['piles']}: {m['fate']}, slot {m['slot']}, "
              f"certified {m['certified']}, look harder {m.get('look_harder', False)}")
    for f in ("look_harder.tsv", "loci.tsv"):
        if os.path.exists(f"{D}/port.{f}"):
            for r in tsv(f"{D}/port.{f}"):
                print(f"   {f.split('.')[0]}: {dict(r)}")


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("pick")
    p.add_argument("--out", required=True)
    s = sub.add_parser("sim")
    s.add_argument("--out", required=True)
    s.add_argument("--name", required=True)
    s.add_argument("--case", required=True, choices=["missing", "unannotated", "control"])
    s.add_argument("--s-alleles", default="AB", choices=["AB", "A"])
    s.add_argument("--e-alleles", default="AB", choices=["AB", "A"])
    s.add_argument("--reads-per-allele", type=int, default=100)
    s.add_argument("--e-reads-per-allele", type=int, default=100)
    s.add_argument("--seed", type=int, default=1)
    s.add_argument("--s-het", type=int, default=0, help="synthetic: survivors with < 2 allele differences get a B allele with "
                   "this many random substitutions")
    s.add_argument("--e-frac", type=float, default=1.0, help="keep each substitution of E against its closest survivor with this "
                   "probability (else the survivor's base): E's divergence scaled down")
    s.add_argument("--indiv-hom", type=float, default=0.0, help="the sample's own variants on both alleles, per base")
    s.add_argument("--indiv-het", type=float, default=0.0, help="the sample's own variants on one allele, per base")
    w = sub.add_parser("sweep")
    w.add_argument("--out", required=True)
    w.add_argument("--select", default="", help="run only the grid configurations whose name matches this regex")
    r = sub.add_parser("report")
    r.add_argument("--out", required=True)
    r.add_argument("--name", required=True)
    a = ap.parse_args(argv)
    {"pick": pick, "sim": sim, "sweep": sweep, "report": report}[a.cmd](a)


if __name__ == "__main__":
    main(sys.argv[1:])
