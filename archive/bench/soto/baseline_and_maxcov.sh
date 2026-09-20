#!/bin/bash
# TWO jobs, strictly SERIAL in one script (WSL2 crash rule: never two heavy runs at once).
#
#   ARM A  re-establish the Soto baseline on the ACTUAL shipped recipe. The committed 76.2% in
#          soto_member_detection.tsv was measured on --cross-chrom, which is now deprecated; after that
#          retirement the default is the homology (E_r) catalog, so every delta scored against the old
#          number is meaningless. This run becomes the reference.
#
#   ARM B  the first end-to-end test of RUSTLE_ER_COVERAGE_LONGER — coverage divided by the LONGER
#          sequence (BLAST qcovs AND scovs) instead of the shorter. Three independent measurements motivate
#          it (0.985-precision certificate; the mechanism behind the only RNA-vs-DNA partition
#          contradiction; the thing that would make component merging monotone). It is expected to COST
#          RECALL, because the duplicated unit is size-invariant while annotated spans are not.
#
# Both arms use identical flags and the same binary; only the env knob differs.
set -uo pipefail
CACHE=${SOTO_CACHE:-/home/juanfra/winloci_scratch/soto_cache}
BIN=/mnt/c/Users/jfris/Desktop/Rustle/target/release/gw_family_catalog
FA=${SOTO_FASTA:-/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa}
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT=/home/juanfra/winloci_scratch/maxcov_ab
mkdir -p "$OUT"
FLAGS="--refine"          # homology/E_r is the default now; no mode flag needed

[ -f "$CACHE/soto_regions.bam" ] || { echo "cache missing — run build_soto_cache.sh"; exit 1; }

run () {   # $1=name  $2=extra env ("" for none)
  local name=$1 extra=$2
  echo "=================================================================="
  echo "ARM $name   flags: $FLAGS   env: ${extra:-<none>}"
  date '+  start %H:%M:%S'
  local t0=$SECONDS
  env $extra "$BIN" --bam "$CACHE/soto_regions.bam" --fasta "$FA" $FLAGS \
      --threads 8 --out "$OUT/$name" > "$OUT/$name.log" 2>&1
  local rc=$?
  echo "  exit=$rc  wall=$((SECONDS-t0))s"
  if [ $rc -ne 0 ]; then echo "  FAILED:"; tail -25 "$OUT/$name.log"; return 1; fi
  echo "  copies=$(($(wc -l < "$OUT/$name.copies.tsv")-1))  families=$(($(wc -l < "$OUT/$name.families.tsv")-1))"
  grep -E "provenance|E_r edges" "$OUT/$name.log" | tail -2
  echo "--- member-detection score ---"
  python3 "$HERE/soto_cache_score.py" "$OUT/$name.copies.tsv" 2>&1 | tail -20
  echo
}

run baseline_homology ""                            || exit 1
run maxcov            RUSTLE_ER_COVERAGE_LONGER=1    || exit 1

echo "=================================================================="
echo "SIDE BY SIDE"
echo "=================================================================="
python3 - "$OUT/baseline_homology.copies.tsv" "$OUT/maxcov.copies.tsv" <<'PY'
import sys
from collections import defaultdict
def load(p):
    fam=defaultdict(list); n=0
    for i,ln in enumerate(open(p)):
        if i==0: continue
        f=ln.rstrip("\n").split("\t")
        if len(f)<6: continue
        fam[f[0]].append((f[3],int(f[4]),int(f[5]))); n+=1
    return fam,n
for lab,p in (("baseline(min)",sys.argv[1]),("maxcov(max)",sys.argv[2])):
    fam,n=load(p)
    sizes=sorted(len(v) for v in fam.values())
    spans=sorted(e-s for v in fam.values() for _,s,e in v)
    print(f"{lab:<16} copies={n:<6} families={len(fam):<5} "
          f"median_family_size={sizes[len(sizes)//2]:<3} largest={max(sizes):<4} "
          f"median_copy_span={spans[len(spans)//2]:,}")
PY
echo "DONE"
