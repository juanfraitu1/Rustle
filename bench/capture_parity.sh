#!/bin/bash
# Capture parity events from Rustle + instrumented StringTie on GGO_19 (chr19),
# and gffcompare Rustle's output vs the reference.
# Usage: bash bench/capture_parity.sh
set -euo pipefail
cd /mnt/c/Users/jfris/Desktop/Rustle

CHROM=NC_073243.2
BAM=GGO_19.bam                       # cwd copy (also at ../GGO_19.bam)
REF=/mnt/c/Users/jfris/Desktop/GGO_19.gtf
STEPS=graphnode_list,junction_accept,junction_raw,good_junction,graph_edge

# 1) Rustle de novo + parity log
RUSTLE_PARITY_LOG=/tmp/ru_parity.jsonl \
RUSTLE_PARITY_FILTER_STEPS=$STEPS \
RUSTLE_PARITY_FILTER_CHROM=$CHROM \
  ./target/release/rustle -L "$BAM" -o /tmp/ru_final.gtf 2>/dev/null

# 2) Instrumented StringTie + parity log
STRINGTIE_PARITY_LOG=/tmp/st_parity.jsonl \
STRINGTIE_PARITY_FILTER_STEPS=$STEPS \
STRINGTIE_PARITY_FILTER_CHROM=$CHROM \
  ./tools/stringtie/stringtie -L "$BAM" -o /tmp/st_final.gtf 2>/dev/null

# 3) Split per-step for the diff tools (graphnode_list lines into _gn files)
grep '"step":"graphnode_list"' /tmp/ru_parity.jsonl > /tmp/ru_gn.jsonl || true
grep '"step":"graphnode_list"' /tmp/st_parity.jsonl > /tmp/st_gn.jsonl || true

# 4) gffcompare Rustle vs reference (headline)
gffcompare -r "$REF" -o /tmp/ru_cmp /tmp/ru_final.gtf 2>/dev/null
echo "=== Rustle vs reference ==="
grep -E "Intron chain level:|Transcript level:" /tmp/ru_cmp.stats
echo "(parity logs: /tmp/ru_parity.jsonl /tmp/st_parity.jsonl; finals: /tmp/ru_final.gtf /tmp/st_final.gtf)"
