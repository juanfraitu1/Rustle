#!/bin/bash
# StringTie arm of the tool bakeoff (docs/PREREG_tool_bakeoff_2026-09-08.md).
# Long-read mode (-L) at defaults, no annotation guidance -- the same information we get.
set -euo pipefail
W=/mnt/linuxdisk/home/juanfraitu/bakeoff
BAM=/mnt/linuxdisk/home/juanfraitu/npip_cat/npip3.bam
# StringTie resolution. The vendored `tools/stringtie` submodule was removed from the repo (its URL
# `../stringtie` pointed at a sibling checkout that does not exist, so `git clone --recursive` failed
# for anyone else). Set $STRINGTIE, or install it (`conda install -c bioconda stringtie`) and let PATH
# resolve it. On the original machine the old checkout was moved to ~/Desktop/stringtie.
ST="${STRINGTIE:-$(command -v stringtie || true)}"
[ -x "$ST" ] || { echo "stringtie not found: set \$STRINGTIE or install it (conda install -c bioconda stringtie)" >&2; exit 127; }
mkdir -p "$W/stringtie"
# StringTie needs a coordinate-sorted BAM; npip3.bam already is.
[ -s "$W/stringtie/st.gtf" ] || "$ST" -L -p 4 -o "$W/stringtie/st.gtf" "$BAM" 2> "$W/stringtie/st.err"
grep -vc '^#' "$W/stringtie/st.gtf"
