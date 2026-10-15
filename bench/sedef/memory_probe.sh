#!/bin/bash
# Peak memory of `utilities sd-scan` on chimp chr1 and chr1-3, and its power-law projection to the whole genomes
# (plan docs/superpowers/plans/2026-10-13-o1-layers-phase-1-2.md, Task 4 gate: <= 20 GB to go on).
#   memory_probe.sh UTILITIES_BIN [OUT_DIR]
set -euo pipefail
BIN=$1; O=${2:-/mnt/linuxdisk/tmp/o1_layers/memprobe}; D=/mnt/c/Users/jfris/Desktop/sd/sedef-rs/data
mkdir -p "$O"
if [ ! -s "$O/chimp3.fa" ]; then
  # awk stops after the third record, so zcat ends on SIGPIPE: pipefail is off for this one pipeline
  set +o pipefail
  zcat "$D/mPanTro3-v2.1.fa.gz" | awk '/^>/{n++} n>3{exit} {print}' > "$O/chimp3.fa.part"
  set -o pipefail
  mv "$O/chimp3.fa.part" "$O/chimp3.fa"
fi
[ "$(md5sum < "$O/chimp3.fa" | cut -c1-32)" = 7a1f8b2fb4ae4faf71b6aca2d836af10 ] || { echo "memory_probe: $O/chimp3.fa md5 differs from 7a1f8b2fb4ae4faf71b6aca2d836af10" >&2; exit 1; }
md5sum "$O/chimp3.fa"
for item in chr1:$D/chr1.fa chimp3:$O/chimp3.fa; do
  tag=${item%%:*}; fa=${item#*:}
  /usr/bin/time -f '%e %M' -o "$O/$tag.time" "$BIN" sd-scan --genome "$fa" --out "$O/$tag" --threads 5 2> "$O/$tag.log"
  echo "$tag wall_s peak_kb: $(cat "$O/$tag.time")"
done
# a second chr1-3 run must write the same tables (multi-sequence determinism, Fact 2)
"$BIN" sd-scan --genome "$O/chimp3.fa" --out "$O/chimp3_rerun" --threads 5 2> "$O/chimp3_rerun.log"
for ext in candidates.tsv blocks.tsv pairs.tsv segments.tsv families.tsv duplicons.tsv; do
  cmp -s "$O/chimp3.sd.$ext" "$O/chimp3_rerun.sd.$ext" && echo "chimp3 rerun $ext identical" || echo "chimp3 rerun $ext DIFFERS"
done
grep -E $'^(blocks|pairs|pairs_ok|segments|families|duplicon_rows)\t' "$O/chimp3.sd.params.tsv" | tr '\n' ' '; echo
python3 - "$O" <<'EOF'
import math, sys
o = sys.argv[1]
pts = []
for tag in ("chr1", "chimp3"):
    bases = int(next(l.split("\t")[1] for l in open(f"{o}/{tag}.sd.params.tsv") if l.startswith("bases\t")))
    pts.append((bases, int(open(f"{o}/{tag}.time").read().split()[1])))
(x1, y1), (x2, y2) = pts
e = math.log(y2 / y1) / math.log(x2 / x1)
for g, size in (("CHM13 v2.0", 3_117_292_070), ("gorilla", 3_545_850_636)):
    print(f"projected peak {g}: {y2 * (size / x2) ** e / 1048576:.1f} GB (exponent {e:.2f})")
EOF
