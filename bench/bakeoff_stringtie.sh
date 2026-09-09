#!/bin/bash
# StringTie arm of the tool bakeoff (docs/PREREG_tool_bakeoff_2026-09-08.md).
# Long-read mode (-L) at defaults, no annotation guidance -- the same information we get.
set -euo pipefail
W=/mnt/linuxdisk/home/juanfraitu/bakeoff
BAM=/mnt/linuxdisk/home/juanfraitu/npip_cat/npip3.bam
ST=/mnt/c/Users/jfris/Desktop/Rustle/tools/stringtie/stringtie
mkdir -p "$W/stringtie"
# StringTie needs a coordinate-sorted BAM; npip3.bam already is.
[ -s "$W/stringtie/st.gtf" ] || "$ST" -L -p 4 -o "$W/stringtie/st.gtf" "$BAM" 2> "$W/stringtie/st.err"
grep -vc '^#' "$W/stringtie/st.gtf"
