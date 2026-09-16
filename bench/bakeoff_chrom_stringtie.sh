#!/bin/bash
# StringTie arm of a human ordinary-chromosome assembler bakeoff (parameterized sibling of
# bench/bakeoff_chr20_stringtie.sh -- DO NOT overwrite that one; see bench/CHR20_ASSEMBLER_COMPARISON.md and
# docs/PREREG_gtf_refine_chr17_2026-09-16.md for context).
# usage: bakeoff_chrom_stringtie.sh CHROM
# Long-read mode (-L) at defaults, no annotation guidance -- same information our tool gets.
set -euo pipefail
CHROM=${1:?CHROM}
W=/mnt/linuxdisk/home/juanfraitu/bakeoff/human_${CHROM}
BAM="$W/${CHROM}.bam"
ST=/mnt/c/Users/jfris/Desktop/Rustle/tools/stringtie/stringtie
mkdir -p "$W/stringtie"
[ -s "$W/stringtie/st.gtf" ] || "$ST" -L -p 4 -o "$W/stringtie/st.gtf" "$BAM" 2> "$W/stringtie/st.err"
grep -vc '^#' "$W/stringtie/st.gtf"
