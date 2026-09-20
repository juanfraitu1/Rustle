#!/usr/bin/env python3
"""Is junction-position CO-THREADING a discriminator between TRUE and FALSE homology edges?

Usage: junction_cothread.py <copies.tsv> <copies.fa> [chr1,chr15]

The question this answers. The E_r edge rule accepts a pair on identity + coverage. Coverage is
doing the duplicon-bridge defence, and it is the criterion that excludes borderline TRUE members
(a short locus aligning over only 0.40 of the longer one is rejected even when it is a real
paralog). So: is there a signal ORTHOGONAL to coverage that says "these two loci are the same
gene" -- one that would let the coverage floor be relaxed without admitting the false edges the
floor is currently holding off?

The candidate is junction POSITION. Intron/exon LENGTHS drift after duplication, so a
length-agreement test fails; but the SPLICE SITES themselves are inherited, so if two loci are
copies of the same gene, each junction of one should land on a junction of the other when
projected through the alignment. A shared repeat or domain produces alignment WITHOUT shared
junction structure, which is exactly the false edge coverage is there to reject.

Measurement, per ordered pair with an alignment:
    co-thread = (# junctions of the query that project to within TOL bp of a target junction)
                / (# junctions of the query inside the aligned block)

Labels come from Soto: a pair of loci both overlapping annotated family members is TRUE if the
two members are in the same Soto family, FALSE otherwise.

STRAND. A minus-strand locus's exon-sum is stored reverse-complemented, so a junction at
cumulative exonic position p sits at (total - p) in the stored frame and the order reverses.
Getting this wrong makes every plus/minus pair score 0.00 and looks like a real signal.
"""
import subprocess
import sys
from collections import defaultdict
from itertools import combinations

COP, FA = sys.argv[1], sys.argv[2]
CHROMS = set(sys.argv[3].split(",")) if len(sys.argv) > 3 else None
BED = "bench/soto/80_fams.gene_preferred.bed"
TOL = 10  # bp; a junction must project this close to count as co-threaded

# ---- copies, with junction positions in exon-sum coordinates ------------------------------------
copies = {}   # key "chrom:start-end" -> dict
for i, ln in enumerate(open(COP)):
    if i == 0:
        continue
    f = ln.rstrip("\n").split("\t")
    if len(f) < 10:
        continue
    fid, chrom, s, e, strand, exons = f[0], f[3], int(f[4]), int(f[5]), f[7], f[9]
    if CHROMS and chrom not in CHROMS:
        continue
    lens = [int(b) - int(a) for a, b in (x.split("-") for x in exons.split(","))]
    if strand == "-":
        lens = lens[::-1]          # stored sequence is the reverse complement
    cum, acc = [], 0
    for L in lens[:-1]:
        acc += L
        cum.append(acc)
    copies[f"{chrom}:{s}-{e}"] = {"fam": fid, "junc": cum, "chrom": chrom, "s": s, "e": e}

# ---- truth labels --------------------------------------------------------------------------------
members = []
for ln in open(BED):
    c, ms, me, name, *_ = ln.rstrip("\n").split("\t")
    if CHROMS and c not in CHROMS:
        continue
    members.append((c, int(ms), int(me), name.split("|")[0], name.split("|")[1]))
for k, d in copies.items():
    hit = [(min(d["e"], me) - max(d["s"], ms), g, fam)
           for (c, ms, me, g, fam) in members if c == d["chrom"] and d["s"] < me and d["e"] > ms]
    d["truth"] = max(hit)[2] if hit else None
    d["gene"] = max(hit)[1] if hit else None

# ---- all-vs-all alignment of the exon-sums --------------------------------------------------------
paf = subprocess.run(
    ["minimap2", "-x", "asm20", "-c", "--eqx", "-t", "4", "-X", "-N", "50", "-p", "0.1", FA, FA],
    capture_output=True, text=True, check=True).stdout


def cigar_ops(cg):
    n = 0
    for ch in cg:
        if ch.isdigit():
            n = n * 10 + ord(ch) - 48
        else:
            yield n, ch
            n = 0


def project_all(cg, qs, ts, qjunc):
    """Project query junctions through the alignment; returns (n_in_block, [target positions])."""
    blocks, qp, tp = [], qs, ts
    for n, o in cigar_ops(cg):
        if o in "=XM":
            blocks.append((qp, qp + n, tp))
            qp += n
            tp += n
        elif o in "DN":
            tp += n
        elif o == "I":
            qp += n
    inside = [j for j in qjunc if qs <= j < qp]
    out = []
    for j in inside:
        for a, b, t in blocks:
            if a <= j < b:
                out.append(t + (j - a))
                break
    return len(inside), out


def key_of(name):
    return name.split("|")[2] if "|" in name else name


best = {}
for ln in paf.splitlines():
    f = ln.split("\t")
    cg = next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), None)
    if cg is None:
        continue
    q, t = key_of(f[0]), key_of(f[5])
    if q not in copies or t not in copies:
        continue
    qs, qe, ts = int(f[2]), int(f[3]), int(f[7])
    cov = (qe - qs) / min(int(f[1]), int(f[6]))
    ident = 1.0 - float(next((x[5:] for x in f[12:] if x.startswith("de:f:")), "1"))
    n_in, projected = project_all(cg, qs, ts, copies[q]["junc"])
    tj = copies[t]["junc"]
    hit = sum(1 for p in projected if any(abs(p - x) <= TOL for x in tj))
    frac = hit / n_in if n_in else None
    k = tuple(sorted((q, t)))
    rec = (cov, ident, frac, n_in)
    if k not in best or cov > best[k][0]:
        best[k] = rec

# ---- compare TRUE vs FALSE edges ------------------------------------------------------------------
rows = []
for (a, b), (cov, ident, frac, n_in) in best.items():
    ta, tb = copies[a]["truth"], copies[b]["truth"]
    if ta is None or tb is None:
        continue
    rows.append((ta == tb, cov, ident, frac, n_in, a, b))

both_spliced = [r for r in rows if r[3] is not None and r[4] > 0 and copies[r[5]]["junc"] and copies[r[6]]["junc"]]
print(f"aligned pairs with truth labels : {len(rows)}   "
      f"(TRUE {sum(1 for r in rows if r[0])} / FALSE {sum(1 for r in rows if not r[0])})")
print(f"pairs where BOTH loci are spliced: {len(both_spliced)}   "
      f"(TRUE {sum(1 for r in both_spliced if r[0])} / FALSE {sum(1 for r in both_spliced if not r[0])})")
print(f"=> junction co-threading is UNDEFINED for {len(rows)-len(both_spliced)}/{len(rows)} "
      f"= {100*(len(rows)-len(both_spliced))/max(len(rows),1):.0f}% of labelled pairs\n")


def summarise(tag, sel):
    for lab, want in (("TRUE ", True), ("FALSE", False)):
        v = sorted(r[3] for r in sel if r[0] is want)
        if not v:
            print(f"  {tag} {lab}: n=0")
            continue
        med = v[len(v) // 2]
        print(f"  {tag} {lab}: n={len(v):3d}  median={med:.2f}  "
              f"frac>=0.80: {100*sum(1 for x in v if x >= 0.80)/len(v):.0f}%  "
              f"frac<=0.20: {100*sum(1 for x in v if x <= 0.20)/len(v):.0f}%")


print("junction co-threading, pairs where it is defined:")
summarise("all      ", both_spliced)
lowcov = [r for r in both_spliced if r[1] < 0.50]
print(f"\nthe pairs the coverage floor currently REJECTS (cov < 0.50), n={len(lowcov)}:")
summarise("cov<0.50 ", lowcov)

def auc(sel, idx):
    pos = [r[idx] for r in sel if r[0]]
    neg = [r[idx] for r in sel if not r[0]]
    if not pos or not neg:
        return None
    return sum((p > n) + 0.5 * (p == n) for p in pos for n in neg) / (len(pos) * len(neg))


for tag, sel in (("all pairs where defined", both_spliced), ("only cov < 0.50", lowcov)):
    a_j, a_c, a_i = auc(sel, 3), auc(sel, 1), auc(sel, 2)
    n_t = sum(1 for r in sel if r[0])
    n_f = len(sel) - n_t
    print(f"\nAUC, {tag}  (TRUE {n_t} / FALSE {n_f})   [0.5 = no signal]")
    for nm, v in (("co-threading", a_j), ("coverage    ", a_c), ("identity    ", a_i)):
        print(f"  {nm} = {v:.3f}" if v is not None else f"  {nm} = n/a")

# What a relaxed rule would actually do to THIS pair set. Only the pairs where the two rules differ
# matter, so they are counted directly rather than inferred from the AUCs.
print("\nedge-level effect on the labelled, aligned pairs (NOT a pipeline result -- the >=2-loci")
print("gate means the pipeline outcome can differ; see EDGE_RULE_DETECTION_COUPLING.txt):")
defined = {(r[5], r[6]) for r in both_spliced}
for name, keep in (
    ("cov >= 0.50  (shipped)", lambda r: r[1] >= 0.50),
    ("cov >= 0.30", lambda r: r[1] >= 0.30),
    ("cov >= 0.50 OR co-thread >= 0.80", lambda r: r[1] >= 0.50 or (r[3] is not None and r[3] >= 0.80)),
):
    t = sum(1 for r in rows if r[0] and keep(r))
    f = sum(1 for r in rows if not r[0] and keep(r))
    print(f"  {name:34s} TRUE edges {t:3d}   FALSE edges {f:3d}   "
          f"precision {t/max(t+f,1):.3f}")
