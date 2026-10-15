#!/bin/bash
# `utilities sd-scan` vs the sandbox sedef-rs on chimp chr1 (plan docs/superpowers/plans/2026-10-13-o1-layers-phase-1-2.md,
# Task 4): on this one-sequence genome, candidates, blocks, pairs (= the sandbox's refined) and segments byte-identical
# after the header-name change ("chr1 mPanTro3-v2.1" -> "chr1"); on several sequences the vendored chain and dedup fixes
# change the output, so byte equality is not expected there (src/sedef/mod.rs). families.tsv is compared by its size
# profile, since the vendored families join a segment's two copies (src/sedef/mod.rs).
#   equivalence_chr1.sh UTILITIES_BIN [SANDBOX_DATA_DIR] [OUT_DIR]
set -euo pipefail
BIN=$1; D=${2:-/mnt/c/Users/jfris/Desktop/sd/sedef-rs/data}; O=${3:-/mnt/linuxdisk/tmp/o1_layers/equiv_chr1}
mkdir -p "$O"
/usr/bin/time -f $'wall_s\t%e\npeak_rss_kb\t%M' -o "$O/time.tsv" "$BIN" sd-scan --genome "$D/chr1.fa" --out "$O/chr1" --threads 5 2> "$O/log"
rc=0
for pair in candidates:candidates blocks:blocks refined:pairs segments:segments; do
  ref=${pair%%:*}; ours=${pair#*:}
  if cmp -s <(sed 's/chr1 mPanTro3-v2.1/chr1/g' "$D/chr1.$ref.tsv") "$O/chr1.sd.$ours.tsv"; then echo "$ours identical"; else echo "$ours DIFFERS"; rc=1; fi
done
awk -F'\t' '!/^#/{n++; c[n]=$2; tot+=$2; sp+=$3; if($2==1)s++} END{asort(c); printf "families %d one-copy %d copies %d largest %d %d %d %d %d span_sum %d\n", n, s, tot, c[n], c[n-1], c[n-2], c[n-3], c[n-4], sp}' "$O/chr1.sd.families.tsv"
cat "$O/time.tsv"
exit $rc
