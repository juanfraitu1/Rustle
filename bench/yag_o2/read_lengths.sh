#!/bin/bash
# Read lengths (longest record with a sequence) of the molecules named in NAMES, from BAM over REGION.
#   read_lengths.sh BAM REGION NAMES OUT.tsv
set -euo pipefail
samtools view -@2 -N "$3" "$1" "$2" \
  | awk -F'\t' '$10 != "*" { l = length($10); if (l > m[$1]) m[$1] = l } END { for (n in m) print n "\t" m[n] }' > "$4"
