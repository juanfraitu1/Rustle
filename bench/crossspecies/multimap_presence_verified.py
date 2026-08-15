#!/usr/bin/env python3
"""Same test as multimap_presence.py, but with DNA-VERIFIED labels.

The first run scored AUC 0.9413 against labels that assumed "not in a tested family" == "single copy".
Spot-checking the four top-ranked negatives found THREE genuine segmental duplications (98.9% id over
32 kb, 98.4% over 18 kb, 99.3% over 4.7 kb), so that AUC was a FLOOR on a mislabelled negative set.

Here every candidate negative is relabelled by an INDEPENDENT criterion -- genome self-alignment -- and
the AUC is recomputed. The relabelling uses DNA paralogy, never the sec_frac signal being tested, so it
is not circular.

label: duplicated  <=>  a paralogous block of >= MIN_BLOCK bp at >= MIN_ID identity exists elsewhere
                        in the genome (best self-hit excluded)

usage: multimap_presence_verified.py BAM GENE_BED GENOME FAM1,FAM2,... [n_candidates]
"""
import random
import subprocess
import sys
from collections import defaultdict

bam, gene_bed, genome, fams_arg = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
N_CAND = int(sys.argv[5]) if len(sys.argv) > 5 else 60
FAMS = fams_arg.split(",")
MIN_READS, MIN_BLOCK, MIN_ID = 20, 1000, 0.90
WORK = "/home/juanfra/winloci_scratch/seedfam"

genes = []
for line in open(gene_bed):
    f = line.rstrip("\n").split("\t")
    if len(f) >= 4:
        genes.append((f[0], int(f[1]), int(f[2]), f[3]))

pos = [g for g in genes if any(g[3].startswith(f) for f in FAMS)]
posn = {g[3] for g in pos}
chroms = {g[0] for g in pos}
cand = [g for g in genes if g[0] in chroms and g[3] not in posn
        and not any(g[3].startswith(f) for f in FAMS)
        and 1000 <= g[2] - g[1] <= 200000]
random.seed(0)
random.shuffle(cand)


def counts(c, s, e):
    pri = subprocess.run(["samtools", "view", "-c", "-F", "2308", bam, f"{c}:{s+1}-{e}"],
                         capture_output=True, text=True).stdout.strip()
    sec = subprocess.run(["samtools", "view", "-c", "-f", "256", "-F", "2052", bam,
                          f"{c}:{s+1}-{e}"], capture_output=True, text=True).stdout.strip()
    return int(pri or 0), int(sec or 0)


rows = []
for c, s, e, n in pos:
    p, sc = counts(c, s, e)
    if p >= MIN_READS:
        rows.append(["FAM", n, c, s, e, p, sc, sc / (p + sc) if p + sc else 0.0])
taken = 0
for c, s, e, n in cand:
    if taken >= N_CAND:
        break
    p, sc = counts(c, s, e)
    if p >= MIN_READS:
        rows.append(["CAND", n, c, s, e, p, sc, sc / (p + sc) if p + sc else 0.0])
        taken += 1

with open(f"{WORK}/presence_loci.fa", "w") as fh:
    for tag, n, c, s, e, p, sc, f in rows:
        out = subprocess.run(["samtools", "faidx", genome, f"{c}:{s+1}-{e}"],
                             capture_output=True, text=True).stdout.splitlines()
        fh.write(f">{n}\n" + "\n".join(out[1:]) + "\n")

subprocess.run(["minimap2", "-x", "asm20", "-c", "--eqx", "-N", "100", "-p", "0.02",
                "-I", "2G", "-t", "4", genome, f"{WORK}/presence_loci.fa"],
               stdout=open(f"{WORK}/presence_loci.paf", "w"),
               stderr=subprocess.DEVNULL)

hits = defaultdict(list)
for line in open(f"{WORK}/presence_loci.paf"):
    f = line.rstrip("\n").split("\t")
    if len(f) < 12:
        continue
    q, nm, bl = f[0], int(f[9]), int(f[10])
    if bl >= MIN_BLOCK and nm / bl >= MIN_ID:
        hits[q].append(bl)

dup_by_dna = {}
for tag, n, c, s, e, p, sc, f in rows:
    h = sorted(hits.get(n, []), reverse=True)
    dup_by_dna[n] = len(h[1:]) > 0     # drop the self hit

P = [r for r in rows if dup_by_dna[r[1]]]
N = [r for r in rows if not dup_by_dna[r[1]]]
flipped = [r[1] for r in rows if r[0] == "CAND" and dup_by_dna[r[1]]]
missed = [r[1] for r in rows if r[0] == "FAM" and not dup_by_dna[r[1]]]

print(f"DNA-VERIFIED LABELS  (paralogous block >= {MIN_BLOCK} bp at >= {MIN_ID} identity)")
print(f"  family members    : {sum(1 for r in rows if r[0]=='FAM')}")
print(f"  candidate others  : {sum(1 for r in rows if r[0]=='CAND')}")
print(f"  -> duplicated     : {len(P)}")
print(f"  -> single-copy    : {len(N)}")
print(f"  candidates RELABELLED as duplicated: {len(flipped)}"
      + (f"  ({', '.join(flipped[:10])}{' ...' if len(flipped)>10 else ''})" if flipped else ""))
if missed:
    print(f"  ⚠ family members with NO DNA paralog: {', '.join(missed)}")


def med(v):
    v = sorted(v)
    return v[len(v) // 2] if v else float("nan")


pf, nf = [r[7] for r in P], [r[7] for r in N]
print(f"\n{'':<14}{'n':>5}{'median':>10}{'min':>10}{'max':>10}")
print(f"{'duplicated':<14}{len(pf):>5}{med(pf):>10.4f}{min(pf):>10.4f}{max(pf):>10.4f}")
print(f"{'single-copy':<14}{len(nf):>5}{med(nf):>10.4f}{min(nf):>10.4f}{max(nf):>10.4f}")

pairs = wins = 0
for a in pf:
    for b in nf:
        pairs += 1
        wins += 1 if a > b else (0.5 if a == b else 0)
print(f"\n  AUC = {(wins/pairs if pairs else float('nan')):.4f}   over {pairs:,} pairs")
print(f"  (mislabelled-negative version scored 0.9413 -- that was a FLOOR)")
best = (0, None)
for t in sorted(set(pf + nf)):
    tp = sum(1 for x in pf if x >= t)
    fp = sum(1 for x in nf if x >= t)
    j = tp / max(len(pf), 1) - fp / max(len(nf), 1)
    if j > best[0]:
        best = (j, t, tp, fp)
if best[1] is not None:
    _, t, tp, fp = best
    print(f"  best single cut sec_frac >= {t:.4f}: {tp}/{len(pf)} duplicated kept, "
          f"{fp}/{len(nf)} single-copy admitted")
print("\n  ⚠ A CUT IS NOT PROPOSED. sec_frac is reported per locus as EVIDENCE; the value above is a\n"
      "    descriptive separation, not a shipped threshold.")
