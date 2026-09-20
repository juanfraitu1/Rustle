#!/bin/bash
# P1 -- SEED-INVARIANCE.  Is F(s) a FIXED POINT of the seeding operator?
#
#   CLAIM: for every s' in F(s),  F(s') = F(s).
#
# If it holds, the family does not depend on WHICH annotated member you happened to possess, which is
# the property that separates "annotation as seed" from "annotation as answer". This is the direct
# answer to "isn't the seed arbitrary?"
#
# All 19 seeds go through ONE minimap2 pass as 19 query sequences, so the genome index is built once
# (~90s) instead of 19 times (~27 min). Per-seed loci are then derived independently from that PAF.
#
# ⚠ E_r ONLY. No read data enters this test -- seed-invariance is a property of the SUPPORT relation
#   (O1). The shared-read certificate is measured separately and never feeds back into the definition.
set -uo pipefail
OUT=/home/juanfra/winloci_scratch/seedfam
HERE="$(cd "$(dirname "$0")" && pwd)"
HS=/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0.fa
T=4

echo "=== building 19 seeds, one per annotated NPIP gene ==="
: > "$OUT/allseeds.fa"
while IFS=$'\t' read -r c s e g rest; do
  samtools faidx "$HS" "$c:$((s+1))-$e" 2>/dev/null \
    | sed "1s/.*/>${g}/" >> "$OUT/allseeds.fa"
done < "$OUT/hs_npip.bed"
grep -c '^>' "$OUT/allseeds.fa"

echo "=== one genome pass, 19 queries ==="
minimap2 -x asm20 -c --eqx -N 200 -p 0.02 -I 2G -t $T \
         "$HS" "$OUT/allseeds.fa" > "$OUT/allseeds_hits.paf" 2>/dev/null
echo "  records: $(wc -l < "$OUT/allseeds_hits.paf")"

python3 "$HERE/seed_invariance_report.py" "$OUT"
echo "DONE"
