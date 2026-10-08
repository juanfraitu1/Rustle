#!/usr/bin/env python3
"""Arm S of Amendment 1 (docs/archive/2026-10/PREREG_rna_allele_haplotype_count_2026-10-01.md): truth transcripts of the S_fam copies on both haplotypes
of KB3781 and simulated Iso-Seq-like reads from them.

Haplotype A transcript = the copy's exon union on `_pri` (identical to the haplotype `_pri` took that chromosome from). Haplotype B
transcript = the same exon blocks lifted block by block through the frozen primary `_pri` -> B alignments (cs walk); a copy whose blocks
do not all lift at >= 95% gets no B transcript. B-only loci (truth_fam.tsv) = the family copy's exon-sum hit at that locus, its aligned
target blocks concatenated. Reads: 20 per transcript; 5' shortened by Exp(mean 100), 3' by Exp(mean 10); 0.5% substitutions, 0.25%
insertions, 0.25% deletions; seed 1.

    python3 isocon_sim.py --work /mnt/linuxdisk/tmp/rna_allele --pri GGO.fasta --out-dir isocon
"""
import argparse
import csv
import math
import random
import re

import pysam

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CS = re.compile(r"(:\d+|\*[a-z][a-z]|[+-][a-z]+|~[a-z]{2}\d+[a-z]{2})")


def rc(s):
    return s.translate(COMP)[::-1]


def blocks(s):
    return [tuple(int(x) for x in b.split("-")) for b in s.split(",") if b]


def lift_blocks(paf, pri, bl):
    """Per exon block: (t_lo, t_hi, strand) on the B chromosome, or None if < 95% of the block's bases align."""
    recs = []
    for ln in open(paf):
        f = ln.rstrip("\n").split("\t")
        if "tp:A:P" not in f[12:]:
            continue
        off = int(f[0].rsplit(":", 1)[1])
        q0, q1 = off + int(f[2]), off + int(f[3])
        if q1 <= bl[0][0] or q0 >= bl[-1][1]:
            continue
        cs = next(x[5:] for x in f[12:] if x.startswith("cs:Z:"))
        recs.append((q0, q1, f[4], cs, int(f[7])))
    out = []
    for x, y in bl:
        best = None                      # Amendment 2: the single record covering most of the block
        for q0, q1, strand, cs, t in recs:
            if q1 <= x or q0 >= y:
                continue
            tpos = []
            q = q0 if strand == "+" else q1
            for op in CS.findall(cs):
                if op[0] == ":" or op[0] == "*":
                    n = int(op[1:]) if op[0] == ":" else 1
                    for k in range(n):
                        qq = q + k if strand == "+" else q - 1 - k
                        if x <= qq < y:
                            tpos.append(t + k)
                    q = q + n if strand == "+" else q - n
                    t += n
                elif op[0] == "+":
                    n = len(op) - 1
                    q = q + n if strand == "+" else q - n
                elif op[0] == "-":
                    t += len(op) - 1
            if tpos and (best is None or len(tpos) > len(best[0])):
                best = (tpos, strand)
        if best is None or len(best[0]) < 0.95 * (y - x):
            return None
        tpos, st = best
        if max(tpos) + 1 - min(tpos) > 1.2 * (y - x) + 100:
            return None
        out.append((min(tpos), max(tpos) + 1, st))
    return out


def mutate(seq, rng):
    out = []
    for b in seq:
        r = rng.random()
        if r < 0.005:
            out.append(rng.choice([c for c in "ACGT" if c != b]))
        elif r < 0.0075:
            out.append(b); out.append(rng.choice("ACGT"))
        elif r < 0.01:
            continue
        else:
            out.append(b)
    return "".join(out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("work", "pri", "out_dir"):
        ap.add_argument("--" + k.replace("_", "-"), required=True)
    a = ap.parse_args(argv)
    W = a.work
    cm = {r["pri"]: r for r in csv.DictReader(open(f"{W}/chrmap.tsv"), delimiter="\t")}
    genes = {r["gene_id"]: r for r in csv.DictReader(open(f"{W}/genes.tsv"), delimiter="\t") if r["set"] == "S_fam"}
    truth = {r["gene_id"]: r for r in csv.DictReader(open(f"{W}/truth.tsv"), delimiter="\t")}
    pri = pysam.FastaFile(a.pri)
    # Amendment 2: each copy's primary annotated transcript (most junctions, ties: longest) from the 2026-09-29 truth GTF
    txex = {}
    for ln in open("/mnt/linuxdisk/tmp/rustle_figures_dev/copy_recovery_tools/ann/truth.ggo.gtf"):
        f = ln.rstrip("\n").split("\t")
        if len(f) > 8 and f[2] == "exon":
            gid = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
            tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
            txex.setdefault(gid, {}).setdefault(tid, []).append((int(f[3]) - 1, int(f[4])))
    prim = {g: sorted(max(t.values(), key=lambda e: (len(e), sum(y - x for x, y in e)))) for g, t in txex.items()}
    tx = {}   # id -> (family, copy, hap, seq)
    for g, r in sorted(genes.items()):
        fam = r["biotype"].split(":", 1)[1]
        bl = prim.get(g) or blocks(r["exons"])
        sA = "".join(pri.fetch(r["chrom"], x, y).upper() for x, y in bl)
        tx[f"{g}_A"] = (fam, g, "A", rc(sA) if r["strand"] == "-" else sA)
        if truth[g]["class"] in ("T2d", "T2i"):
            c = cm[r["chrom"]]
            lb = lift_blocks(f"{W}/out/chr{c['chrom']}.paf", r["chrom"], bl)
            if lb:
                bfa = pysam.FastaFile(f"{W}/t/chr{c['chrom']}.fa")
                parts = []
                for t0, t1, st in lb:
                    s = bfa.fetch(c["B_name"], t0, t1).upper()
                    parts.append(s if st == "+" else rc(s))
                sB = "".join(parts)
                tx[f"{g}_B"] = (fam, g, "B", rc(sB) if r["strand"] == "-" else sB)
    # B-only loci
    alias = {}
    for h in ("mat", "pat"):
        for acc, num, _ in csv.reader(open(f"{W}/{h}.len.tsv"), delimiter="\t"):
            alias[(h, num)] = acc
    hfa = {h: pysam.FastaFile(f"/mnt/linuxdisk/home/juanfraitu/gorilla_haps/{h}.fa") for h in ("mat", "pat")}
    for fr in csv.DictReader(open(f"{W}/truth_fam.tsv"), delimiter="\t"):
        for k, loc in enumerate(x for x in fr["B_only_detail"].split(";") if x):
            acc, se = loc.rsplit(":", 1)
            s0, e0 = (int(v) for v in se.split("-"))
            best = None
            for h in ("mat", "pat"):
                for ln in open(f"{W}/t1fam.{h}.paf"):
                    f = ln.rstrip("\n").split("\t")
                    m = re.fullmatch(r"chr(\w+?)_(mat|pat)_hsa[^_]+", f[5])
                    if not m or alias.get((m.group(2), m.group(1))) != acc or f[0] not in genes:
                        continue
                    if int(f[7]) < e0 and s0 < int(f[8]):
                        ident = int(f[9]) / max(1, int(f[10]))
                        if best is None or ident > best[0]:
                            best = (ident, h, f)
            if best:
                _, h, f = best
                cg = next(x[5:] for x in f[12:] if x.startswith("cg:Z:"))
                t, segs = int(f[7]), []
                for n, op in re.findall(r"(\d+)([MIDNSHX=])", cg):
                    n = int(n)
                    if op in "M=X":
                        segs.append((t, t + n)); t += n
                    elif op in "DN":
                        t += n
                s = "".join(hfa[h].fetch(acc, x, y).upper() for x, y in segs)   # index name chrN_<hap>_hsa* = this accession
                tx[f"{fr['family']}_Bonly{k}"] = (fr["family"], f[0], "Bonly", s if f[4] == "+" else rc(s))
    rng = random.Random(1)
    for fam in ("NPIP", "TBC1D3"):
        with open(f"{a.out_dir}/{fam}.truth.fa", "w") as tf, open(f"{a.out_dir}/{fam}.sim.fa", "w") as rf:
            for tid, (f_, g, h, s) in sorted(tx.items()):
                if f_ != fam:
                    continue
                tf.write(f">{tid}\n{s}\n")
                for i in range(20):
                    a5 = min(int(rng.expovariate(1 / 100)), len(s) // 2)
                    a3 = min(int(rng.expovariate(1 / 10)), len(s) // 4)
                    rf.write(f">{tid}_r{i:02d}\n{mutate(s[a5:len(s) - a3], rng)}\n")
        n = sum(1 for v in tx.values() if v[0] == fam)
        print(fam, n, "truth transcripts;", sum(1 for v in tx.values() if v[0] == fam and v[2] == "B"), "with a B allele;",
              sum(1 for v in tx.values() if v[0] == fam and v[2] == "Bonly"), "B-only")


if __name__ == "__main__":
    main()
