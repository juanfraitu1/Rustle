#!/bin/bash
# DOES THE 9/19 FAILURE BELONG TO RNA, OR TO THE ANNOTATED TRANSCRIPT?
#
# What actually failed earlier was RefSeq NPIPB11: 8 exons, 6,117 bp, used as a probe -> 9/19.
# That is ANNOTATION, not RNA. This compares three probes built at the SAME single seed locus:
#
#   A  RefSeq transcript          6,117 bp   (measured: 9/19)
#   B  READ-SUPPORTED exons       IsoSeq primary reads at NPIPB11, >=3 reads  <-- the thesis substrate
#   C  full genomic span         25,154 bp   (measured: 19/19, density 1.000)
#
# If B lands near C, the loss was the annotation's, and RNA is not the weak link.
#
# ⚠ INPUT PROVENANCE (subset-BAM trap, 3 prior retractions): chr16_sub.bam is the WHOLE-chr16
#   slice of A119b.t2t.bam and its -F 2308 primary counts match the full-BAM audit exactly at all
#   8 genes with recorded reference values. All 19 NPIP genes are on chr16. Safe.
# ⚠ -F 2308 (primary only) before any per-read CIGAR statistic.
# ⚠ The multimap census at the end is a CERTIFICATE (E_c = evidence the loci belong together),
#   NEVER a family definition. O1 is E_r homology; E_c is the ambiguity oracle and belongs to O2.
set -uo pipefail
OUT=/home/juanfra/winloci_scratch/seedfam
HERE="$(cd "$(dirname "$0")" && pwd)"
BAM=/home/juanfra/winloci_scratch/chr16_sub.bam
HS=/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0.fa
SEED="chr16:29663351-29688504"     # NPIPB11, the ONLY locus given
T=4

echo "=== building probe B: read-supported exons at the seed locus ==="
samtools view -F 2308 "$BAM" "$SEED" \
  | python3 "$HERE/read_exons.py" > "$OUT/B11_readexons.raw.bed"
sort -k1,1 -k2,2n "$OUT/B11_readexons.raw.bed" \
  | bedtools merge -c 4 -o count \
  | awk -F'\t' '$4>=3 && ($3-$2)>=20' > "$OUT/B11_readexons.bed"
echo "  read-supported exon blocks (>=3 reads): $(wc -l < "$OUT/B11_readexons.bed")"
echo "  exonic bp: $(awk '{s+=$3-$2}END{print s+0}' "$OUT/B11_readexons.bed")   (RefSeq transcript was 6117)"

bedtools getfasta -fi "$HS" -bed "$OUT/B11_readexons.bed" -fo "$OUT/B11_probe_parts.fa"
python3 -c "
s=''.join(l.strip() for l in open('$OUT/B11_probe_parts.fa') if l[0]!='>')
open('$OUT/B11_probeB.fa','w').write('>NPIPB11_readexons\n'+s+'\n'); print('  probe B bp',len(s))"

echo
echo "=== probe B -> genome ==="
minimap2 -x asm20 -c --eqx -N 200 -p 0.02 -I 2G -t $T "$HS" "$OUT/B11_probeB.fa" \
  > "$OUT/probeB_hits.paf" 2>/dev/null
QL=$(grep -v '^>' "$OUT/B11_probeB.fa" | tr -d '\n' | wc -c)
awk -F'\t' '$11>=200 && ($10/$11)>=0.80 {print $6"\t"$8"\t"$9"\t"$10}' "$OUT/probeB_hits.paf" \
  | sort -k1,1 -k2,2n > "$OUT/probeB_hits.bed"
bedtools merge -d 10000 -i "$OUT/probeB_hits.bed" -c 4 -o sum \
  | awk -v q="$QL" -F'\t' '$4 >= 0.30*q' > "$OUT/probeB_loci.bed"
echo "  loci recovered (>=30% of probe): $(wc -l < "$OUT/probeB_loci.bed")"

echo
echo "=== multimap census: reads at NPIP loci touching >1 NPIP gene (CERTIFICATE ONLY) ==="
python3 "$HERE/multimap_census.py" "$BAM" "$OUT/hs_npip.bed" > "$OUT/multimap_census.txt" 2>&1
cat "$OUT/multimap_census.txt"

echo
python3 "$HERE/isoseq_probe_report.py" "$OUT"
echo "DONE"
