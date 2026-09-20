#!/bin/bash
# StringTie arm of the human chr20 ordinary-chromosome assembler bakeoff (bench/CHR20_ASSEMBLER_COMPARISON.md).
# Sibling of bakeoff_stringtie.sh (gorilla NPIP substrate) -- DO NOT overwrite that one.
# Long-read mode (-L) at defaults, no annotation guidance -- same information our tool gets.
set -euo pipefail
W=/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20
BAM="$W/chr20.bam"
# StringTie resolution. The vendored `tools/stringtie` submodule was removed from the repo (its URL
# `../stringtie` pointed at a sibling checkout that does not exist, so `git clone --recursive` failed
# for anyone else). Set $STRINGTIE, or install it (`conda install -c bioconda stringtie`) and let PATH
# resolve it. On the original machine the old checkout was moved to ~/Desktop/stringtie.
ST="${STRINGTIE:-$(command -v stringtie || true)}"
[ -x "$ST" ] || { echo "stringtie not found: set \$STRINGTIE or install it (conda install -c bioconda stringtie)" >&2; exit 127; }
mkdir -p "$W/stringtie"
[ -s "$W/stringtie/st.gtf" ] || "$ST" -L -p 4 -o "$W/stringtie/st.gtf" "$BAM" 2> "$W/stringtie/st.err"
grep -vc '^#' "$W/stringtie/st.gtf"
