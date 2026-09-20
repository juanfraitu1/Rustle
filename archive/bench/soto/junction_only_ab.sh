#!/bin/bash
# A/B the LOCUS DEFINITION: phase-1-only (junction sharing) vs today's phase-1 + span-aware merge.
#
# The two arms differ in exactly one thing, RUSTLE_LOCUS_JUNCTION_ONLY, which stops `collapse_parent`
# after the junction-share union-find. Everything else -- flags, binary, BAM, rep construction -- is held
# fixed, so any delta is attributable to the locus definition alone.
#
# Both arms run SERIALLY in this one script (WSL2 crash rule: never two heavy jobs at once).
set -uo pipefail
CACHE=${SOTO_CACHE:-/home/juanfra/winloci_scratch/soto_cache}
BIN=/mnt/c/Users/jfris/Desktop/Rustle/target/release/gw_family_catalog
FA=${SOTO_FASTA:-/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa}
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT=/home/juanfra/winloci_scratch/junction_only_ab
mkdir -p "$OUT"
FLAGS="--homology-primary --refine"

run_arm () {           # $1=name  $2=extra env assignment ("" for none)
  local name=$1 extra=$2
  echo "=================================================================="
  echo "ARM $name   flags: $FLAGS   env: ${extra:-<none>} RUSTLE_COTHREAD_REP=1"
  echo "=================================================================="
  local t0=$SECONDS
  env RUSTLE_COTHREAD_REP=1 $extra \
    "$BIN" --bam "$CACHE/soto_regions.bam" --fasta "$FA" $FLAGS \
    --threads 8 --out "$OUT/$name" > "$OUT/$name.log" 2>&1
  local rc=$?
  echo "  exit=$rc  wall=$((SECONDS-t0))s"
  if [ $rc -ne 0 ]; then echo "  FAILED — tail of log:"; tail -20 "$OUT/$name.log"; return 1; fi
  echo "  copies=$(($(wc -l < "$OUT/$name.copies.tsv")-1))  families=$(($(wc -l < "$OUT/$name.families.tsv")-1))"
  grep -E "cothread-rep|locus reps" "$OUT/$name.log" | tail -2
  echo "--- score ---"
  python3 "$HERE/soto_cache_score.py" "$OUT/$name.copies.tsv" 2>&1 | tail -25
  echo
}

run_arm baseline      ""                          || exit 1
run_arm junction_only RUSTLE_LOCUS_JUNCTION_ONLY=1 || exit 1

echo "=================================================================="
echo "LOCUS SIZE COMPARISON (copy span distribution, both arms)"
echo "=================================================================="
python3 - "$OUT/baseline.copies.tsv" "$OUT/junction_only.copies.tsv" <<'PY'
import sys, statistics as st
def spans(p):
    out=[]
    for i,ln in enumerate(open(p)):
        if i==0: continue
        f=ln.rstrip("\n").split("\t")
        if len(f)<6: continue
        try: out.append(int(f[5])-int(f[4]))
        except ValueError: pass
    return out
for lab,p in (("baseline",sys.argv[1]),("junction_only",sys.argv[2])):
    s=spans(p)
    if not s: print(f"{lab}: no copies"); continue
    s.sort()
    print(f"{lab:<15} n={len(s):<6} median={st.median(s):>9,.0f}  p90={s[int(.9*len(s))]:>9,.0f}  "
          f"max={max(s):>10,}  >100kb={sum(1 for x in s if x>100_000)}")
PY
echo "DONE"
