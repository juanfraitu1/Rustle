#!/usr/bin/env bash
# Authoritative DIRECT sweep (no LLM agents) over all candidates: per-locus eval
# + (for wins) read-backing verify. Deterministic and reproducible.
# Usage: winloci_sweep.sh CANDIDATES.json OUTDIR [PARALLEL=3]
set -uo pipefail
CJ="${1:?candidates.json}"; OUT="${2:?outdir}"; P="${3:-3}"
SH=/mnt/c/Users/jfris/Desktop/Rustle/bench/multi_copy_eval
GI=/mnt/c/Users/jfris/Desktop/Rustle/bench/paralog_secondary_scan/scan_out/gene_introns.tsv
mkdir -p "$OUT"
# emit TSV: gene chrom start end strand sib(blank if intergenic)
python3 - "$CJ" > "$OUT/args.tsv" <<'PY'
import json,sys
for c in json.load(open(sys.argv[1])):
    sib=c['dom_sibling'] if ':' not in c['dom_sibling'] else ''
    print('\t'.join(map(str,[c['gene_id'],c['chrom'],c['start'],c['end'],c['strand'],sib])))
PY
n=$(wc -l < "$OUT/args.tsv"); echo "[sweep] $n candidates, parallel=$P"

run_one() {
  IFS=$'\t' read -r g c s e st sib <<< "$1"
  GI="$GI" bash "$SH/winloci_eval.sh" "$g" "$c" "$s" "$e" "$st" "$sib" 2>/dev/null > "$OUT/$g.eval.json"
  # if win_vs_st, run verify
  cls=$(python3 -c "import json;print(json.load(open('$OUT/$g.eval.json')).get('classification',''))" 2>/dev/null)
  if [ "$cls" = "win_vs_st" ]; then
    GI="$GI" bash "$SH/winloci_verify.sh" "$g" "$c" "$s" "$e" "$st" 2>/dev/null > "$OUT/$g.verify.json"
  fi
  echo "  done $g [$cls]"
}
export -f run_one; export OUT SH GI

# concurrency-gated loop
i=0
while IFS= read -r line; do
  run_one "$line" &
  i=$((i+1))
  if [ $((i % P)) -eq 0 ]; then wait; fi
done < "$OUT/args.tsv"
wait
echo "[sweep] complete: $(ls $OUT/*.eval.json 2>/dev/null|wc -l) evals, $(ls $OUT/*.verify.json 2>/dev/null|wc -l) verifies"
