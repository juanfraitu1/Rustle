#!/bin/bash
# SEED SUBSTRATE TEST: annotated GENOMIC SPAN vs the seed gene's mRNA.
#
# Motivation (measured, not assumed): seeding with the 124.6 kb annotated NPIPB11 genomic span
# returns 25 gorilla loci, but only 10 of them carry ANY NPIP CODING sequence. The other 15 are
# 83-98% covered by aligned SEED sequence while holding 0% of the mRNA -- i.e. the edges were
# bought by the non-coding duplicated flank (LCR16a-scale material) inside the annotated span.
#
# This is the NODE-EXTENT problem again, and the fix is a SUBSTRATE change, not a taxon prior:
# seed with the transcript. Species-agnostic, no oracle, no per-family threshold.
#
# ⚠ Do NOT filter on "species X should have N copies". Literature copy numbers for morpheus/NPIP
#   are human 15, chimp 25-30, gorilla 17, orangutan 9; 1-2 copies is OLD WORLD MONKEYS. A rule
#   keyed to an expected count would delete true gorilla copies.
# ⚠ SERIAL, one heavy job at a time (WSL2 crash rule). Outputs to winloci_scratch.
set -uo pipefail
OUT=/home/juanfra/winloci_scratch/seedfam
HERE="$(cd "$(dirname "$0")" && pwd)"
GGO=/home/juanfra/winloci_scratch/GGO.fasta
HS=/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0.fa
T=4

mrna_seed () {   # $1=tag  $2=genome  $3=gene symbol  $4=gff
  local tag=$1 gen=$2 sym=$3 gff=$4
  echo "=================================================================="
  echo "[$tag] mRNA seed = $sym"
  local t0=$SECONDS
  awk -F'\t' '$3=="exon"' "$gff" | grep "gene=$sym;" \
    | awk -F'\t' '{print $1"\t"($4-1)"\t"$5}' | sort -k1,1 -k2,2n \
    | bedtools merge > "$OUT/$tag.cds.bed"
  bedtools getfasta -fi "$gen" -bed "$OUT/$tag.cds.bed" -fo "$OUT/$tag.cds_parts.fa"
  python3 -c "
import sys
s=''.join(l.strip() for l in open('$OUT/$tag.cds_parts.fa') if l[0]!='>')
open('$OUT/$tag.mrna.fa','w').write('>${sym}_mrna\n'+s+'\n')
print('  exons %d  mRNA bp %d' % (sum(1 for _ in open('$OUT/$tag.cds.bed')), len(s)))"

  minimap2 -x asm20 -c --eqx -N 200 -p 0.02 -I 2G -t $T \
           "$gen" "$OUT/$tag.mrna.fa" > "$OUT/$tag.mrnahits.paf" 2>/dev/null
  # loci = merged targets carrying >=30% of the mRNA at id>=0.80
  awk -F'\t' '$11>=200 && ($10/$11)>=0.80 {print $6"\t"$8"\t"$9"\t"$10}' "$OUT/$tag.mrnahits.paf" \
    | sort -k1,1 -k2,2n > "$OUT/$tag.mrnahits.bed"
  local QL; QL=$(grep -v '^>' "$OUT/$tag.mrna.fa" | tr -d '\n' | wc -c)
  bedtools merge -d 10000 -i "$OUT/$tag.mrnahits.bed" -c 4 -o sum \
    | awk -v q="$QL" -F'\t' '$4 >= 0.30*q' > "$OUT/$tag.cdsloci.bed"
  echo "  mRNA-seeded loci (>=30% of mRNA): $(wc -l < "$OUT/$tag.cdsloci.bed")   ($((SECONDS-t0))s)"

  awk -F'\t' '{print $1":"($2+1)"-"$3}' "$OUT/$tag.cdsloci.bed" > "$OUT/$tag.cdsloci.regions"
  samtools faidx "$gen" -r "$OUT/$tag.cdsloci.regions" > "$OUT/$tag.cdsloci.fa" 2>/dev/null
  samtools faidx "$OUT/$tag.cdsloci.fa"
# ⚠ TIER CORRECTION 2026-08-10 (defect B1): the all-vs-all E_r pass below now runs the SHIPPED
# command -- `minimap2 -c -X --no-long-join -t N -k 11 -w 5`, ONE tier. It previously ran
# `-c --eqx -N 200 -p 0.02` with NO `-X` plus an `-x asm20` leg that the shipped default SKIPS
# (RUSTLE_ER_SENSITIVE_ONLY has defaulted true since 2026-08-07). `-N`/`-p` are inert here; `-X`
# is the operative difference (it implies --dual=no). See bench/er_tier.sh for the full rationale.
# GENOME passes in this script are deliberately UNCHANGED: there -X would be wrong.
  minimap2 -c -X --no-long-join -t $T -k 11 -w 5 "$OUT/$tag.cdsloci.fa" "$OUT/$tag.cdsloci.fa" \
           > "$OUT/$tag.cds_ava.paf" 2>/dev/null
  echo "  all-vs-all records = $(wc -l < "$OUT/$tag.cds_ava.paf")   total $((SECONDS-t0))s"
  echo
}

mrna_seed GGO "$GGO" NPIPB11 /home/juanfra/winloci_scratch/GGO_tx.gff
zcat /mnt/c/Users/jfris/Desktop/Reference/chm13v2.0_RefSeq_full.gff.gz > "$OUT/hs.gff" 2>/dev/null
mrna_seed HS  "$HS"  NPIPB11 "$OUT/hs.gff"

python3 "$HERE/seed_cds_report.py" "$OUT"
echo "DONE"
