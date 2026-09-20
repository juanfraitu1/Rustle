#!/usr/bin/env bash
# Reproduces Soto et al. 2025's disclosed "redundant gene" detection command exactly
# (section_I&II/A_SD98_regions.md, github.com/mydennislab/HSD_brain_evolution, fetched 2026-09-11):
#
#   bedtools intersect -wao -s -f 0.9 -a CHM13.combined.v4.txs.sd98.bed \
#   -b CHM13.combined.v4.txs.sd98.bed | awk '{if($4!=$10){print}}' | cut -f4,10
#
# "focusing on transcripts... fully contained (90%) within another gene with a different gene ID" --
# candidates for their manually-curated 71-gene removal (docs/o1_ledger.md §6ig). Confirmed: run against
# this project's own independently-recomputed SD98 gene set, this candidate list contains 100% (71/71) of
# the genes that differ between that recomputed set and Soto's own published 1,793-gene table -- the
# DETECTION step is fully reproducible. The final SELECTION among candidates is NOT: tested a biotype-
# priority + transcript-length heuristic and got 67.6% recall at 9.7% precision (§6ig) -- Soto's own text
# says the deciding criterion is whether a candidate "encoded an alternative protein", which is a
# protein-level judgment no coordinate-only rule can supply. Use this script to get the CANDIDATE list;
# do not treat its output as a final removal list without further (manual, or protein-level) review.
#
# Usage: soto_redundant_gene_candidates.sh <cat_v4.bed> <geneset.tsv-with-gene_id-column> <out.tsv>
#   <cat_v4.bed>   BED12, CAT v4 transcript annotation (gene_id in column 19, 1-indexed)
#   <geneset.tsv>  TSV with a header row and a `gene_id` column -- restricts which genes' transcripts
#                  enter the self-intersect (matching their own ".sd98"-suffixed input file's scope)
#   <out.tsv>      two columns: gene_id_A, gene_id_B (A's transcript is >=90% contained in B's, same
#                  strand, different gene ID) -- a CANDIDATE list, not a final removal list
set -euo pipefail
CAT_BED="$1"; GENESET="$2"; OUT="$3"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

tail -n +2 "$GENESET" | cut -f1 | sort -u > "$TMP/geneids.txt"

# clean BED6 (chrom start end GENE_ID score strand) so bedtools' -wao output columns 4/10 land on the
# A/B gene IDs exactly -- a 12-column BED here would put B's name at a different offset, a real bug this
# project hit on the first attempt (docs/o1_ledger.md §6ig).
awk -F'\t' 'NR==FNR{keep[$1]=1; next} $19 in keep {print $1"\t"$2"\t"$3"\t"$19"\t"$5"\t"$6}' \
  "$TMP/geneids.txt" "$CAT_BED" | sort -k1,1 -k2,2n > "$TMP/txs.bed6"

bedtools intersect -wao -s -f 0.9 -a "$TMP/txs.bed6" -b "$TMP/txs.bed6" \
  | awk -F'\t' '{if($4!=$10){print $4"\t"$10}}' | sort -u > "$OUT"

echo "[done] $(wc -l < "$OUT") candidate pairs -> $OUT" >&2
