#!/usr/bin/env python3
"""Does the PRESENCE of multimappers at a locus mark it as a duplicated (family-member) locus?

The user's actual proposal -- presence, not width. An earlier attempt measured the multimap FOOTPRINT
and failed (CV 0.482 vs annotation 0.412), and failed for a structural reason worth keeping: adjacent
copies return the IDENTICAL footprint (NPIPA6 == NPIPA7, NPIPB6 == NPIPB7, NPIPB8 == NPIPB9), because
the track marks the duplicated REGION, not a copy. Width cannot resolve a copy boundary. Presence can
still answer a different and more useful question: IS THIS LOCUS DUPLICATED AT ALL?

SIGNAL (per locus, no NH tag needed, no pairwise information):
    sec_frac = secondary records / (primary + secondary) records overlapping the locus
A secondary alignment here means "a read whose best placement is elsewhere also fits HERE". It never
names where, so nothing from E_c enters. This is a node-CONSTRUCTION signal, not an edge rule.

⚠ minimap2 emits no NH tag; sec_frac is the NH-free equivalent of "these reads are multimappers".
⚠ The negative set is drawn from expressed chr16 genes NOT in any tested family. Some will be genuinely
  duplicated by chance -- that makes the test CONSERVATIVE, not wrong, and is reported.

usage: multimap_presence.py BAM GENE_BED FAM1,FAM2,... [n_negatives]
"""
import random
import subprocess
import sys

bam, gene_bed, fams_arg = sys.argv[1], sys.argv[2], sys.argv[3]
N_NEG = int(sys.argv[4]) if len(sys.argv) > 4 else 60
FAMS = fams_arg.split(",")
MIN_READS = 20

genes = []
for line in open(gene_bed):
    f = line.rstrip("\n").split("\t")
    if len(f) >= 4:
        genes.append((f[0], int(f[1]), int(f[2]), f[3]))


def counts(c, s, e):
    """(primary, secondary) alignment records overlapping the interval."""
    pri = subprocess.run(["samtools", "view", "-c", "-F", "2308", bam, f"{c}:{s+1}-{e}"],
                         capture_output=True, text=True).stdout.strip()
    sec = subprocess.run(["samtools", "view", "-c", "-f", "256", "-F", "2052", bam,
                          f"{c}:{s+1}-{e}"], capture_output=True, text=True).stdout.strip()
    return int(pri or 0), int(sec or 0)


pos = [g for g in genes if any(g[3].startswith(f) for f in FAMS)]
posn = {g[3] for g in pos}
chroms = {g[0] for g in pos}
cand = [g for g in genes if g[0] in chroms and g[3] not in posn
        and not any(g[3].startswith(f) for f in FAMS)
        and 1000 <= g[2] - g[1] <= 200000]
random.seed(0)
random.shuffle(cand)


def rows_for(gs, label, want=None):
    out = []
    for c, s, e, n in gs:
        p, sc = counts(c, s, e)
        if p < MIN_READS:
            continue
        out.append((n, c, s, e, p, sc, sc / (p + sc) if p + sc else 0.0))
        if want and len(out) >= want:
            break
    return out


P = rows_for(pos, "dup")
N = rows_for(cand, "single", want=N_NEG)

print(f"POSITIVES ({len(P)}) -- known multi-copy family members")
print(f"{'gene':<12}{'primary':>9}{'secondary':>11}{'sec_frac':>10}")
for n, c, s, e, p, sc, f in sorted(P, key=lambda x: -x[6]):
    print(f"{n:<12}{p:>9,}{sc:>11,}{f:>10.4f}")

print(f"\nNEGATIVES ({len(N)}) -- expressed genes on the same chromosome(s), not in any tested family")
Ns = sorted(N, key=lambda x: -x[6])
print(f"{'gene':<12}{'primary':>9}{'secondary':>11}{'sec_frac':>10}")
for n, c, s, e, p, sc, f in Ns[:10]:
    print(f"{n:<12}{p:>9,}{sc:>11,}{f:>10.4f}")
print(f"  ... ({len(N)-10} more, lower)" if len(N) > 10 else "")


def med(v):
    v = sorted(v)
    return v[len(v) // 2] if v else float("nan")


pf = [r[6] for r in P]
nf = [r[6] for r in N]
print(f"\n{'':<12}{'n':>5}{'median':>10}{'min':>10}{'max':>10}")
print(f"{'duplicated':<12}{len(pf):>5}{med(pf):>10.4f}{min(pf):>10.4f}{max(pf):>10.4f}")
print(f"{'single-copy':<12}{len(nf):>5}{med(nf):>10.4f}{min(nf):>10.4f}{max(nf):>10.4f}")

# AUC by rank (no threshold chosen)
pairs = wins = 0
for a in pf:
    for b in nf:
        pairs += 1
        wins += 1 if a > b else (0.5 if a == b else 0)
auc = wins / pairs if pairs else float("nan")
print(f"\n  AUC (P(dup ranks above single-copy)) = {auc:.4f}   over {pairs:,} pairs")
sep = min(pf) > max(nf)
print(f"  perfectly separable by sec_frac? {'YES' if sep else 'NO'}"
      + (f"   (gap {min(pf):.4f} vs {max(nf):.4f})" if sep else
         f"   (dup min {min(pf):.4f} <= single max {max(nf):.4f})"))
over = [r[0] for r in N if r[6] >= med(pf)]
if over:
    print(f"  ⚠ negatives at or above the duplicated median (likely genuinely duplicated, so this is\n"
          f"    a conservative test, not a false-positive count): {', '.join(over[:12])}")
