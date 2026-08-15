#!/bin/bash
# SEEDED-DNA FAMILY DEFINITION, one seed gene -> family, run identically on two species.
#
# CLAIM UNDER TEST: annotation supplies the SEED, the method supplies MEMBERSHIP and the EDGES.
#   Non-circular because the members recovered were NOT given: we hand the pipeline exactly ONE
#   annotated locus and ask what else in the genome belongs with it.
#
# WHY GORILLA IS THE INTERESTING PANEL: the NCBI mGorGor1 annotation names exactly ONE NPIP gene
#   (NPIPB11); CHM13 names NINETEEN. Gorilla is the real "complete genome, incomplete annotation"
#   case, and its validation truth comes from HUMAN annotation -- a source independent of the
#   gorilla annotation we seeded from.
#
# ⚠ Each species is seeded with ITS OWN annotated NPIPB11 span. Those spans differ 5x
#   (GGO 124.6 kb vs CHM13 25.2 kb). That is annotation slop, not biology -- reported, not hidden.
# ⚠ Repeat guard is an ALIGNED-BASE FLOOR (>=5 kb of seed sequence per candidate locus), not a
#   softmask gate: both genomes are ~53% softmasked at NPIP, so a softmask gate would delete
#   the family. 5 kb sits far below the ~16 kb shared cassette and far above Alu-scale hits.
# ⚠ SERIAL. One heavy job at a time (WSL2 crash rule). Outputs to winloci_scratch.
set -uo pipefail
OUT=/home/juanfra/winloci_scratch/seedfam
HERE="$(cd "$(dirname "$0")" && pwd)"
GGO=/home/juanfra/winloci_scratch/GGO.fasta
HS=/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0.fa
T=4
mkdir -p "$OUT"

# ---------------------------------------------------------------- per-species discovery
discover () {         # $1=tag  $2=genome  $3=seed region (samtools faidx syntax)
  local tag=$1 gen=$2 region=$3
  echo "=================================================================="
  echo "[$tag] seed = $region"
  local t0=$SECONDS

  samtools faidx "$gen" "$region" > "$OUT/$tag.seed.fa" 2>/dev/null
  echo "  seed bp = $(grep -v '^>' "$OUT/$tag.seed.fa" | tr -d '\n' | wc -c)"

  # 1. WHERE ELSE DOES THE SEED LIVE?  genome-wide, high -N / low -p so paralogs are not suppressed.
  if [ ! -s "$OUT/$tag.seedhits.paf" ]; then
    minimap2 -x asm20 -c --eqx -N 200 -p 0.02 -I 2G -t $T \
             "$gen" "$OUT/$tag.seed.fa" > "$OUT/$tag.seedhits.paf" 2> "$OUT/$tag.mm2.log"
  fi
  echo "  seed hit records = $(wc -l < "$OUT/$tag.seedhits.paf")   ($((SECONDS-t0))s)"

  # 2. CANDIDATE LOCI = merged target intervals carrying >=5 kb of aligned seed.
  awk -F'\t' '$11>=1000 && ($10/$11)>=0.80 {print $6"\t"$8"\t"$9"\t"$10}' \
      "$OUT/$tag.seedhits.paf" | sort -k1,1 -k2,2n > "$OUT/$tag.hits.bed"
  bedtools merge -d 10000 -i "$OUT/$tag.hits.bed" -c 4 -o sum \
      | awk -F'\t' '$4>=5000 && ($3-$2)>=5000' > "$OUT/$tag.loci.bed"
  echo "  candidate loci = $(wc -l < "$OUT/$tag.loci.bed")"

  awk -F'\t' '{print $1":"($2+1)"-"$3}' "$OUT/$tag.loci.bed" > "$OUT/$tag.loci.regions"
  samtools faidx "$gen" -r "$OUT/$tag.loci.regions" > "$OUT/$tag.loci.fa" 2>/dev/null
  samtools faidx "$OUT/$tag.loci.fa"

  # 3. E_r = the SHIPPED tier: sensitive (-k 11 -w 5) ONLY. The comment that used to sit here said
  #    "UNION of the asm20 tier and the sensitive tier, as shipped" -- that was WRONG on both counts.
# ⚠ TIER CORRECTION 2026-08-10 (defect B1): the all-vs-all E_r pass below now runs the SHIPPED
# command -- `minimap2 -c -X --no-long-join -t N -k 11 -w 5`, ONE tier. It previously ran
# `-c --eqx -N 200 -p 0.02` with NO `-X` plus an `-x asm20` leg that the shipped default SKIPS
# (RUSTLE_ER_SENSITIVE_ONLY has defaulted true since 2026-08-07). `-N`/`-p` are inert here; `-X`
# is the operative difference (it implies --dual=no). See bench/er_tier.sh for the full rationale.
# GENOME passes in this script are deliberately UNCHANGED: there -X would be wrong.
  minimap2 -c -X --no-long-join -t $T -k 11 -w 5 "$OUT/$tag.loci.fa" "$OUT/$tag.loci.fa" \
           > "$OUT/$tag.ava.paf" 2>/dev/null
  echo "  all-vs-all records = $(wc -l < "$OUT/$tag.ava.paf")   total $((SECONDS-t0))s"
  echo
}

# ---------------------------------------------------------------- run (SERIAL)
discover GGO "$GGO" "NC_073242.2:31337557-31462161"   # gorilla NPIPB11, the ONLY named NPIP
discover HS  "$HS"  "chr16:29663351-29688504"         # human   NPIPB11, 1 of 19 named

# ---------------------------------------------------------------- validation substrate
# Human NPIP gene sequences: independent truth for the GORILLA panel.
awk -F'\t' '{print $1":"($2+1)"-"$3}' "$OUT/hs_npip.bed" > "$OUT/hs_npip.regions"
samtools faidx "$HS" -r "$OUT/hs_npip.regions" > "$OUT/hs_npip_genes.fa" 2>/dev/null
samtools faidx "$OUT/hs_npip_genes.fa"

# gorilla loci -> human NPIP genes (cross-species orthology check)
minimap2 -x asm20 -c --eqx -N 50 -p 0.05 -t $T "$OUT/hs_npip_genes.fa" "$OUT/GGO.loci.fa" \
         > "$OUT/GGO_vs_hsnpip.paf" 2>/dev/null
echo "GGO loci -> human NPIP records = $(wc -l < "$OUT/GGO_vs_hsnpip.paf")"

# what does the GORILLA annotation call these loci? (should be mostly LOC..., i.e. unnamed)
awk -F'\t' '$3=="gene"' /home/juanfra/winloci_scratch/GGO_tx.gff \
  | awk -F'\t' '{n=split($9,a,";");g=".";for(i=1;i<=n;i++)if(a[i]~/^gene=/&&length(a[i])>5)g=substr(a[i],6);
                 if(g=="")g="."; print $1"\t"($4-1)"\t"$5"\t"g}' | sort -k1,1 -k2,2n > "$OUT/ggo_genes.bed"
# ⚠ the "." default is load-bearing: an empty 4th column leaves a TRAILING TAB and bedtools
#   rejects the whole file ("wrong number of fields"), which silently reads as "0 overlaps".
bedtools intersect -a "$OUT/GGO.loci.bed" -b "$OUT/ggo_genes.bed" -wa -wb > "$OUT/GGO_loci_genes.tsv"

echo "=================================================================="
python3 "$HERE/seed_family_report.py" "$OUT"
echo "DONE"
