#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT=$(cd "$(dirname "$0")/../.." && pwd)
OUT_DIR="$REPO_ROOT/bench/o1_outgroup_rooting_poc"
WORK_DIR="$OUT_DIR/work"
DATA_ROOT=/mnt/linuxdisk/home/juanfraitu
HUMAN_DATA="$DATA_ROOT/winloci_data"
GORILLA_DATA="$DATA_ROOT/gorilla_haps"

mountpoint -q /mnt/linuxdisk
test -r "$GORILLA_DATA/mat.fa"
test -r "$GORILLA_DATA/pat.fa"
mkdir -p "$WORK_DIR"
cd "$REPO_ROOT"

python3 bench/o1_outgroup_rooting_poc.py prepare \
  --human-fasta "$HUMAN_DATA/chm13v2.0.fa" \
  --human-gff "$HUMAN_DATA/HSA_genomic.gff" \
  --family-loci bench/o1_provenance_witness_prototype/HSA.GOLGA6_8.loci.tsv \
  --out-dir "$WORK_DIR" --flank-bp 25000

minimap2 -x asm20 -c --eqx --secondary=yes -N 100 \
  "$WORK_DIR/human_family_loci.fa" "$WORK_DIR/human_sources.fa" \
  > "$WORK_DIR/human_source_to_family.paf"

# Run serially: each prebuilt index uses about 13.5 GB RSS when loaded.
minimap2 -x asm20 -c --eqx --secondary=yes -N 100 \
  "$HUMAN_DATA/mGorGor1.mat.asm20.rebuild.mmi" "$WORK_DIR/human_all_loci.fa" \
  > "$WORK_DIR/gorilla_mat.loci.asm20.paf"
minimap2 -x asm20 -c --eqx --secondary=yes -N 100 \
  "$HUMAN_DATA/mGorGor1.pat.asm20.mmi" "$WORK_DIR/human_all_loci.fa" \
  > "$WORK_DIR/gorilla_pat.loci.asm20.paf"

# Strict flank orthology uses the FASTAs rather than the sensitive asm20 indexes.
minimap2 -x asm5 -c --eqx --secondary=yes -N 100 \
  "$GORILLA_DATA/mat.fa" "$WORK_DIR/human_all_flanks.fa" \
  > "$WORK_DIR/gorilla_mat.flanks.asm5.paf"
minimap2 -x asm5 -c --eqx --secondary=yes -N 100 \
  "$GORILLA_DATA/pat.fa" "$WORK_DIR/human_all_flanks.fa" \
  > "$WORK_DIR/gorilla_pat.flanks.asm5.paf"

python3 bench/o1_outgroup_rooting_poc.py analyze \
  --probes "$WORK_DIR/probes.tsv" \
  --human-links "$WORK_DIR/human_source_to_family.paf" \
  --mat-loci "$WORK_DIR/gorilla_mat.loci.asm20.paf" \
  --pat-loci "$WORK_DIR/gorilla_pat.loci.asm20.paf" \
  --mat-left "$WORK_DIR/gorilla_mat.flanks.asm5.paf" \
  --mat-right "$WORK_DIR/gorilla_mat.flanks.asm5.paf" \
  --pat-left "$WORK_DIR/gorilla_pat.flanks.asm5.paf" \
  --pat-right "$WORK_DIR/gorilla_pat.flanks.asm5.paf" \
  --out-dir "$OUT_DIR"
