#!/bin/bash
# RESUMABLE genome-wide measurement of dna_family_fallback: run --homology-primary + RUSTLE_DNA_FAMILY_FALLBACK=1
# PER-CHROM (each output persists; completed chroms are skipped on re-run), then combine + score. Survives the
# terminal-closing churn: a killed run only loses the in-flight chrom. Caveat: per-chrom homology misses
# cross-chrom MAIN families, but the fallback's orphan-locus projection is genome-wide per locus, so recall via
# the fallback and the full dna_family precision are both captured.
set -u
CACHE=/home/juanfra/winloci_scratch/soto_cache; PC=$CACHE/perchrom
BIN=/mnt/c/Users/jfris/Desktop/Rustle/target/release/gw_family_catalog
FA=/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa
[ -f "$FA" ] || FA=/mnt/wsl/PHYSICALDRIVE0p2/home/juanfraitu/winloci_data/chm13v2.0.fa
OUT=$CACHE/dna_perchrom; mkdir -p "$OUT"
CHROMS=$(ls "$PC"/chr*.bam 2>/dev/null | sed 's#.*/##;s/.bam//')
for c in $CHROMS; do
  if [ -f "$OUT/$c.copies.tsv" ]; then echo "[$c] done (skip)"; continue; fi
  RUSTLE_DNA_FAMILY_FALLBACK=1 timeout 1200 "$BIN" --bam "$PC/$c.bam" --fasta "$FA" --homology-primary --out "$OUT/$c" > "$OUT/$c.log" 2>&1
  echo "[$c] copies=$(($(wc -l <"$OUT/$c.copies.tsv" 2>/dev/null||echo 1)-1))  dna_family=$(($(wc -l <"$OUT/$c.dna_family.tsv" 2>/dev/null||echo 1)-1))  rc=$?"
done
# combine
{ echo -e "family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads"
  for c in $CHROMS; do [ -f "$OUT/$c.copies.tsv" ] && tail -n +2 "$OUT/$c.copies.tsv"; done; } > "$CACHE/homprim_pc.copies.tsv"
{ echo -e "family_id\tchrom\tstart\tend\tfamCN\tmin_locus_reads\tstatus\tprojection_loci"
  for c in $CHROMS; do [ -f "$OUT/$c.dna_family.tsv" ] && tail -n +2 "$OUT/$c.dna_family.tsv"; done; } > "$CACHE/homprim_pc.dna_family.tsv"
echo "=== combined: copies=$(($(wc -l <"$CACHE/homprim_pc.copies.tsv")-1))  dna_family=$(($(wc -l <"$CACHE/homprim_pc.dna_family.tsv")-1)) ==="
python3 /mnt/c/Users/jfris/Desktop/Rustle/bench/soto/soto_dna_family_score.py "$CACHE/homprim_pc.copies.tsv" "$CACHE/homprim_pc.dna_family.tsv"
