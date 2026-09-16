#!/bin/bash
# StringTie arm of the human chr20 ordinary-chromosome assembler bakeoff (bench/CHR20_ASSEMBLER_COMPARISON.md).
# Sibling of bakeoff_stringtie.sh (gorilla NPIP substrate) -- DO NOT overwrite that one.
# Long-read mode (-L) at defaults, no annotation guidance -- same information our tool gets.
set -euo pipefail
W=/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20
BAM="$W/chr20.bam"
ST=/mnt/c/Users/jfris/Desktop/Rustle/tools/stringtie/stringtie
mkdir -p "$W/stringtie"
[ -s "$W/stringtie/st.gtf" ] || "$ST" -L -p 4 -o "$W/stringtie/st.gtf" "$BAM" 2> "$W/stringtie/st.err"
grep -vc '^#' "$W/stringtie/st.gtf"
