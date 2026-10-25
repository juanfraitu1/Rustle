#!/bin/bash
# o3port_golden.sh — the Rust copy_phasing against the Python o3port.py on frozen inputs (plan docs/superpowers/plans/2026-10-24-o3port-rust.md).
#   o3port_golden.sh py <name> [seed]   the Python side (PYTHONHASHSEED = seed, default 0) -> $GW/<name>/py.*  (seed != 0: py<seed>.*)
#   o3port_golden.sh rs <name>          the Rust side ($RS_BIN, --threads $RS_THREADS) -> $GW/<name>/rs.*
#   o3port_golden.sh diff <name> [a b]  byte-compare every product, the dump and the [o3port] stderr lines of prefixes a, b (py, rs)
#   o3port_golden.sh all [filter]       rs + diff on every case of o3port_golden.tsv whose name matches the grep filter
# Cases: o3port_golden.tsv (name, bam, copies, genome, net, isocon or -, flags). Inputs are read-only.
set -euo pipefail
REPO=$(cd "$(dirname "$0")/../.." && pwd)
GW=${GW:-/mnt/linuxdisk/tmp/o3port_golden}
RS_BIN=${RS_BIN:-/mnt/linuxdisk/home/juanfraitu/rustle_target/release/copy_phasing}
RS_THREADS=${RS_THREADS:-1}
PRODUCTS="candidates.tsv contigs.fa clusters.tsv reads.tsv families.tsv putative.tsv putative.fa loci.tsv o2.tsv dump.jsonl"

case_args() {   # the run arguments of case $1, output prefix $2
  awk -F'\t' -v n="$1" 'NR > 1 && $1 == n { print; f = 1 } END { exit !f }' "$REPO/bench/o3_graph/o3port_golden.tsv" |
  while IFS=$'\t' read -r _ bam copies genome net iso flags; do
    a="--bam $bam --copies $copies --genome $genome --net $net --out $2 --dump $2.dump.jsonl"
    [ "$iso" != "-" ] && [ -n "$iso" ] && a="$a --isocon $iso"
    echo "$a $flags"
  done
}

run_side() {   # side (py|rs), name, prefix tag, seed
  local side=$1 n=$2 tag=$3 seed=${4:-0} d=$GW/$2
  mkdir -p "$d"
  local args; args=$(case_args "$n" "$d/$tag") || { echo "no case $n" >&2; exit 2; }
  local t0=$(date +%s.%N)
  if [ "$side" = py ]; then
    PYTHONHASHSEED=$seed python3 "$REPO/bench/o3_graph/o3port.py" run $args 2> "$d/$tag.log"
  else
    "$RS_BIN" run $args --threads "$RS_THREADS" 2> "$d/$tag.log"
  fi
  rm -rf "$d"/tmp*                                    # the Python leaves its homology temp dir beside --out
  printf "%s\t%s\t%.1f s\n" "$n" "$tag" "$(echo "$(date +%s.%N) - $t0" | bc)" >> "$GW/times.tsv"
}

stderr_lines() {   # the per-family progress lines, without the program's name
  grep -E '^\[(o3port|copy_phasing)\]' "$1" | sed -E 's/^\[[a-z0-9_]+\] //'
}

diff_case() {   # name, prefix a, prefix b
  local d=$GW/$1 a=${2:-py} b=${3:-rs} bad=0
  for p in $PRODUCTS; do
    cmp -s "$d/$a.$p" "$d/$b.$p" || { echo "DIFF $1 $p"; bad=1; }
  done
  cmp -s <(stderr_lines "$d/$a.log") <(stderr_lines "$d/$b.log") || { echo "DIFF $1 stderr"; bad=1; }
  [ $bad = 0 ] && echo "SAME $1"
  return $bad
}

cmd=${1:?py|rs|diff|all}; shift
case $cmd in
  py)   s=${2:-0}; run_side py "$1" "$([ "$s" = 0 ] && echo py || echo py$s)" "$s" ;;
  rs)   run_side rs "$1" rs ;;
  diff) diff_case "$1" "${2:-py}" "${3:-rs}" ;;
  all)  bad=0
        for n in $(awk -F'\t' 'NR > 1 { print $1 }' "$REPO/bench/o3_graph/o3port_golden.tsv" | grep -e "${1:-.}"); do
          run_side rs "$n" rs; diff_case "$n" || bad=$((bad + 1))
        done
        echo "cases differing: $bad"; [ $bad = 0 ] ;;
  *) echo "usage: o3port_golden.sh py|rs|diff|all" >&2; exit 2 ;;
esac
