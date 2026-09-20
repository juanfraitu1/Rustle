#!/usr/bin/env bash
# Annotation-seeded DNA mode vs the read-only catalog (ledger §6be).
#
# Establishes: the SEED SUBSTRATE, not the preset, decides recall.
#   spliced CDS seed  -> 14/31   |  genomic seed -> 25/31  |  all-gene multi-seed -> 30/31
#
# Preserved because §6bc's lesson is that an unpreserved instrument's numbers cannot be defended.
# Report (seed, PRESET, -p, -N, node unit) with any copy count.
set -euo pipefail
BASE=${BASE:-/mnt/linuxdisk/home/juanfraitu}
WORK=${WORK:-$BASE/seedmode}
REF=${REF:-$BASE/npip_cat/npip3_contigs.fa}          # the 3 in-scope contigs
GFF=${GFF:-$BASE/winloci_data/GGO_genomic.gff}       # gorilla-native RefSeq
THREADS=${THREADS:-4}
mkdir -p "$WORK"; cd "$WORK"

# ---- seeds: every annotated gene AND pseudogene on the target contigs -------
awk -F'\t' '($1=="NC_073241.2"||$1=="NC_073242.2"||$1=="NC_073244.2") &&
            ($3=="gene"||$3=="pseudogene"){print $1"\t"$4-1"\t"$5"\t"$3"\t"$9}' "$GFF" > genes3.bed
awk -F'\t' '{printf "%s:%d-%d\n",$1,$2+1,$3}' genes3.bed > allgenes.regions
samtools faidx -r allgenes.regions "$REF" > allgenes.fa

# ---- project every seed genome-wide; keep alignment blocks >= 1 kb ----------
# asm10 was measured to add nothing over asm20 here (same 82 families, same 30/31).
for p in asm20 asm10; do
  minimap2 -x "$p" -N 50 -t "$THREADS" "$REF" allgenes.fa 2>/dev/null \
    | awk -v P="$p" 'BEGIN{OFS="\t"} $11>=1000 {print P,$1,$6,$8,$9,$5,$10,$11}'
done > multiseed.paf

# ---- evaluate: self-hits dropped, union-find into families, NPIP recall -----
python3 "$(dirname "$0")/seeded_mode_eval.py" "${1:-0.30}"
