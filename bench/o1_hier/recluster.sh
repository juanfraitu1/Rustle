#!/bin/bash
# Re-cluster one default catalog from its cached all-vs-all (docs/archive/2026-10/merged/PREREG_mcl_two_node_2026-10-13.md F2).
#   recluster.sh BIN SAMPLE OUTDIR [GRAPH_DUMP]
# --fasta / --emit-units are not passed: neither changes clusters.tsv or loci.tsv (human_testis: byte-identical, 5.7 s).
# GRAPH_DUMP: --dump-graph writes the PRE-MCL edges (u<TAB>v<TAB>weight) and changes no other output.
set -euo pipefail
BIN=$1; S=$2; O=$3
R=/mnt/linuxdisk/tmp/rustle_figures/runs/$S
mkdir -p "$O"
"$BIN" --paf "$R/$S.fam.loci.paf" --gff "$R/$S.fam.loci.gff3" \
  --min-exonic-bp 1 --min-shared-exon-frac 0.60 ${4:+--dump-graph "$4"} --out "$O/$S.fam" > "$O/$S.log" 2>&1
