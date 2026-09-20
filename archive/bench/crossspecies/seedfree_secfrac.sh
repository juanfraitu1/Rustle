#!/bin/bash
# CAN MULTIMAPPER PRESENCE REPLACE THE SEED?
#
# Pipeline under test (NO annotation, NO seed):
#     sec_frac scan over chr16  ->  candidate intervals  ->  E_r  ->  connected components
#
# sec_frac = secondary / (primary + secondary) alignment records at a locus. It is a PER-LOCUS SCALAR:
# it says "reads placed elsewhere also fit here" without ever naming where, so no pairwise read-conflict
# (E_c) information enters O1. Multimapping as MAPPABILITY is node construction; multimapping as a
# RELATION is O2 and stays out.
#
# WHY THIS MATTERS: blind node construction tops out at purity 0.237, and that ceiling is the reason the
# definition still needs a seed. This is the only annotation-free node signal measured to work
# (AUC 0.94-0.98 duplicated vs single-copy; all 26 chr16 NPIP loci recovered at FULL SIZE -- median
# interval 20,000 bp against a 20,314 bp target).
#
# WHAT IS NOT YET KNOWN: at top-1000 bins the scan returns 585 intervals of which 31 are on NPIP loci.
# So it is CANDIDATE GENERATION, not locus calling. The open question -- and the whole point of this
# script -- is whether E_r plus the quasi-clique requirement turns those candidates into clean families
# or into a hairball.
#
# PRE-REGISTERED READING OF THE OUTCOME:
#   * NPIP emerges as one dense component, junk falls into small/sparse ones -> seed-free O1 works.
#   * NPIP is inside a large low-density blob -> sec_frac stays a CERTIFICATE and the seed stays.
#   * NPIP fragments across components -> the candidate intervals are the wrong nodes, not the rule.
#
# ⚠ ONE HEAVY JOB AT A TIME (WSL2). Serial by construction; outputs to winloci_scratch.
# ⚠ chr16_sub.bam is the WHOLE-chr16 slice of A119b (verified: -F 2308 counts match the full-BAM audit
#   at all 8 genes with recorded reference values). Not a -M -L subset.
# ⚠ Bin size and the top-N cut are SWEPT, not chosen. An absolute threshold has failed six times here.
set -uo pipefail
OUT=/home/juanfra/winloci_scratch/seedfam/seedfree
BAM=/home/juanfra/winloci_scratch/chr16_sub.bam
HS=/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0.fa
HERE="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$OUT"

echo "=== one streaming pass over chr16: per-bin primary and secondary counts ==="
python3 - "$BAM" "$OUT/bins.tsv" <<'PY'
import subprocess, sys
from collections import defaultdict
bam, out = sys.argv[1], sys.argv[2]
BIN = 5000
prim, sec = defaultdict(int), defaultdict(int)
p = subprocess.Popen(["samtools", "view", "-F", "2052", bam],
                     stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
n = 0
for line in p.stdout:
    i = line.find('\t'); j = line.find('\t', i+1)
    k = line.find('\t', j+1); m = line.find('\t', k+1)
    flag = int(line[i+1:j]); pos = int(line[k+1:m])
    (sec if flag & 256 else prim)[pos // BIN] += 1
    n += 1
p.wait()
with open(out, "w") as fh:
    fh.write("bin\tprimary\tsecondary\tsec_frac\n")
    for b in sorted(set(prim) | set(sec)):
        pr, se = prim[b], sec[b]
        fh.write(f"{b}\t{pr}\t{se}\t{se/(pr+se) if pr+se else 0:.6f}\n")
print(f"  alignments scanned {n:,}; bins written {len(set(prim)|set(sec)):,}")
PY

echo
echo "=== candidates -> E_r -> components, sweeping the top-N cut ==="
python3 "$HERE/seedfree_secfrac.py" "$OUT" "$HS"
echo DONE
