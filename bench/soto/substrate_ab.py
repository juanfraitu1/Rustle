#!/usr/bin/env python3
"""EXON-SUM vs GENOMIC-SPAN as the E_r substrate, on one labelled pair set, one aligner.

Usage: substrate_ab.py <copies.tsv> <copies.fa> <genome.fa> [chr1,chr15]

THE PROBLEM THIS TESTS. The exon-sum substrate is not uniform. A single-exon representative has no
introns to drop, so its "exon-sum" IS its genomic span (measured ratio 1.00); a spliced representative's
exon-sum is ~22% of its span. Since ~75% of representatives are single-exon stubs, E_r spends most of its
pairs comparing GENOMIC sequence on one side to SPLICED sequence on the other. That is the mechanism
behind the class-(b) failures: median best coverage 0.636 when both reps are spliced, 0.099 when one is
single-exon, 0.000 when both are, and 70.6% of the 180 missing pairs have coverage exactly 0.00 -- no
alignment record at all, which no coverage floor can rescue.

Aligning genomic spans on BOTH sides makes the substrate uniform. The stated cost is that introns diverge
faster than exons, so at a fixed identity floor the genomic substrate is STRICTER and can drop older
paralogs. This script measures both directions on the same pairs: what is rescued, and what is lost.

Both substrates get the SAME aligner and the SAME floors, so the only variable is the sequence.
"""
import subprocess
import sys
from collections import defaultdict

COP, FA, GENOME = sys.argv[1], sys.argv[2], sys.argv[3]
CHROMS = set(sys.argv[4].split(",")) if len(sys.argv) > 4 else None
BED = "bench/soto/80_fams.gene_preferred.bed"
OUT = "/tmp/claude-1000/-mnt-c-Users-jfris-Desktop/161b66ca-e17e-4e70-b1cf-8b479d6c328d/scratchpad"

copies = {}
for i, ln in enumerate(open(COP)):
    if i == 0:
        continue
    f = ln.rstrip("\n").split("\t")
    if len(f) < 9:
        continue
    chrom, s, e = f[3], int(f[4]), int(f[5])
    if CHROMS and chrom not in CHROMS:
        continue
    copies[f"{chrom}:{s}-{e}"] = {"fam": f[0], "chrom": chrom, "s": s, "e": e, "nexon": int(f[6])}

members = []
for ln in open(BED):
    c, ms, me, name, *_ = ln.rstrip("\n").split("\t")
    if CHROMS and c not in CHROMS:
        continue
    members.append((c, int(ms), int(me), name.split("|")[1]))
for k, d in copies.items():
    hit = [(min(d["e"], me) - max(d["s"], ms), fam)
           for (c, ms, me, fam) in members if c == d["chrom"] and d["s"] < me and d["e"] > ms]
    d["truth"] = max(hit)[1] if hit else None

# ---- genomic spans, same names as the exon-sum keys -------------------------------------------------
span_fa = f"{OUT}/spans.fa"
regions = "\n".join(f"{d['chrom']}:{d['s']+1}-{d['e']}" for d in copies.values())
open(f"{OUT}/regions.txt", "w").write(regions + "\n")
subprocess.run(f"samtools faidx {GENOME} -r {OUT}/regions.txt > {span_fa}", shell=True, check=True)


def key_of(name):
    return name.split("|")[2] if "|" in name else name


def norm(name):
    """samtools writes chr:START+1-END; map back to the 0-based key used everywhere else."""
    c, rng = name.rsplit(":", 1)
    a, b = rng.split("-")
    return f"{c}:{int(a)-1}-{b}"


def align(path, keyfn):
    # BOTH tiers, as `homology_edges_all_reps` runs them: asm20 cannot seed a pair below ~0.82
    # identity, so asm20 alone silently makes every floor under 0.82 behave identically.
    lines = []
    for tier in (["-x", "asm20"], ["-k", "11", "-w", "5"]):
        lines += subprocess.run(
            ["minimap2", *tier, "-c", "--eqx", "-t", "2", "-X", "-N", "50", "-p", "0.1", path, path],
            capture_output=True, text=True, check=True).stdout.splitlines()
    best = {}
    for ln in lines:
        f = ln.split("\t")
        q, t = keyfn(f[0]), keyfn(f[5])
        if q not in copies or t not in copies or q == t:
            continue
        cov = (int(f[3]) - int(f[2])) / min(int(f[1]), int(f[6]))
        ident = 1.0 - float(next((x[5:] for x in f[12:] if x.startswith("de:f:")), "1"))
        k = tuple(sorted((q, t)))
        if k not in best or cov > best[k][0]:
            best[k] = (cov, ident)
    return best


exon = align(FA, key_of)
span = align(span_fa, norm)

labelled = [(a, b) for a in copies for b in copies if a < b
            and copies[a]["truth"] and copies[b]["truth"]]


def edges(best, min_id, min_cov):
    return {k for k, (cov, ident) in best.items() if cov >= min_cov and ident >= min_id}


def report(tag, best, min_id, min_cov):
    E = edges(best, min_id, min_cov)
    t = sum(1 for (a, b) in labelled if (a, b) in E and copies[a]["truth"] == copies[b]["truth"])
    f = sum(1 for (a, b) in labelled if (a, b) in E and copies[a]["truth"] != copies[b]["truth"])
    print(f"  {tag:28s} TRUE {t:4d}   FALSE {f:4d}   precision {t/max(t+f,1):.3f}")
    return E, t, f


print(f"labelled pairs: {len(labelled)}   copies: {len(copies)}\n")
for min_id in (0.70, 0.80, 0.90, 0.95, 0.98):
    print(f"identity >= {min_id:.2f}, coverage >= 0.50")
    Ee, te, fe = report("exon-sum (shipped)", exon, min_id, 0.50)
    Es, ts, fs = report("genomic span", span, min_id, 0.50)
    print(f"  {'delta':28s} TRUE {ts-te:+4d}   FALSE {fs-fe:+4d}\n")

# The class-(b) question: TRUE pairs with NO exon-sum alignment record at all.
noaln = [(a, b) for (a, b) in labelled
         if copies[a]["truth"] == copies[b]["truth"] and tuple(sorted((a, b))) not in exon]
resc = [p for p in noaln if tuple(sorted(p)) in edges(span, 0.70, 0.50)]
print(f"TRUE pairs with NO exon-sum alignment at all : {len(noaln)}")
print(f"  ...of which the genomic span links (id>=0.70, cov>=0.50): {len(resc)} "
      f"({100*len(resc)/max(len(noaln),1):.0f}%)")

lost = [(a, b) for (a, b) in labelled if copies[a]["truth"] == copies[b]["truth"]
        and tuple(sorted((a, b))) in edges(exon, 0.70, 0.50)
        and tuple(sorted((a, b))) not in edges(span, 0.70, 0.50)]
print(f"TRUE pairs the exon-sum links but the span LOSES (intron divergence): {len(lost)}")

# Does the substrate mismatch actually explain the missing alignments?
print("\nmedian best exon-sum coverage of TRUE pairs, by representative pairing:")
for lab, pred in (("both spliced", lambda a, b: copies[a]["nexon"] > 1 and copies[b]["nexon"] > 1),
                  ("one single-exon", lambda a, b: (copies[a]["nexon"] > 1) != (copies[b]["nexon"] > 1)),
                  ("both single-exon", lambda a, b: copies[a]["nexon"] == 1 and copies[b]["nexon"] == 1)):
    for sub, bst in (("exon-sum", exon), ("span", span)):
        v = sorted(bst.get(tuple(sorted((a, b))), (0.0, 0.0))[0]
                   for (a, b) in labelled if copies[a]["truth"] == copies[b]["truth"] and pred(a, b))
        if v:
            print(f"  {lab:17s} {sub:9s} n={len(v):4d}  median={v[len(v)//2]:.3f}")
