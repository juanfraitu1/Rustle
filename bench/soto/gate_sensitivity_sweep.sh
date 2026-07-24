#!/bin/bash
# Readthrough / mis-chain GATE sensitivity sweep. Re-runs the cached per-chrom Soto recompute under several
# gate settings and tabulates copies / on-member% / member-recall, so we can see whether the gates that stop
# "readthroughs connecting copies" are load-bearing (strong influence) or robust. The gate knobs are read
# from env by the current binary (default = compiled constants; see denovo_pipeline.rs env_num):
#   RUSTLE_KEEP_READTHROUGH=1        -> DISABLE both readthrough + mis-chain gates
#   RUSTLE_READTHROUGH_MIN_DISTINCT  -> junctions a single-exon span must engulf to be dropped (default 5)
#   RUSTLE_MISCHAIN_GIANT_BP         -> "giant intron" bp threshold (default 50000)
#   RUSTLE_MISCHAIN_MIN_READS        -> giant intron dropped if <this many reads (default 3)
# Each config is a full per-chrom recompute (~40 min). Runs sequentially. Outputs gate_sensitivity.tsv.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
CACHE=/home/juanfra/winloci_scratch/soto_cache; PC=$CACHE/perchrom
SUM=$CACHE/gate_sensitivity.tsv
BED=$HERE/80_fams.chr.bed
echo -e "config\tgate_env\tcopies\ton_member%\trecall\tdetected" > $SUM

score_copies() {  # $1 = catalog.copies.tsv -> "copies<TAB>on%"
  python3 - "$1" "$BED" <<'PY'
import sys
from collections import defaultdict
cat, bed = sys.argv[1], sys.argv[2]
mem=defaultdict(list)
for l in open(bed):
    c=l.split("\t"); mem[c[0]].append((int(c[1]),int(c[2])))
def om(ch,s,e): return any(not(a>e or b<s) for a,b in mem.get(ch,()))
n=on=0
for l in open(cat):
    c=l.rstrip("\n").split("\t")
    if c[0]=="family_id" or len(c)<6: continue
    n+=1; on+= 1 if om(c[3],int(c[4]),int(c[5])) else 0
print(f"{n}\t{(100*on//max(n,1))}%")
PY
}

run_config() {  # $1=label  $2=gate_env
  local label="$1" genv="$2"
  echo "=== [$label] $genv ==="
  rm -f "$PC"/chr*.copies.tsv "$PC"/xchrom.copies.tsv   # clear stale (keep cached mini-BAMs)
  env $genv bash "$HERE/recompute_perchrom.sh" 5 > "$CACHE/gate_${label}.log" 2>&1
  cp "$CACHE/perchrom_catalog.copies.tsv" "$CACHE/gate_${label}.copies.tsv"
  local cs; cs=$(score_copies "$CACHE/gate_${label}.copies.tsv")
  local rec; rec=$(python3 "$HERE/soto_cache_score.py" "$CACHE/gate_${label}.copies.tsv" 2>/dev/null \
                    | grep -oE "NEW: [0-9]+/362 = [0-9.]+%" | sed 's/NEW: //')
  echo -e "${label}\t${genv:--(default gates ON)}\t${cs}\t${rec}" >> $SUM
  echo "[$label] copies/on% = $cs   recall = $rec"
}

# configs: default (gates ON) + gates OFF + threshold variations
run_config default            ""
run_config gates_off          "RUSTLE_KEEP_READTHROUGH=1"
run_config readthrough_strict "RUSTLE_READTHROUGH_MIN_DISTINCT=3"
run_config mischain_aggressive "RUSTLE_MISCHAIN_GIANT_BP=20000 RUSTLE_MISCHAIN_MIN_READS=5"

echo "=== GATE SENSITIVITY SUMMARY ==="; cat $SUM
