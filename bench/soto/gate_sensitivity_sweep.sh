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
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
CACHE=/home/juanfra/winloci_scratch/soto_cache; PC=$CACHE/perchrom
SUM=$CACHE/gate_sensitivity.tsv
BED=$HERE/80_fams.chr.bed

# ---- HARD ABORT BLOCK (docs/o1_ledger.md §6am) ------------------------------------------------------
# §6am: this sweep TRUNCATED $SUM (the committed result) and rm -f'd the cached per-chrom outputs BEFORE
# any binary was ever attempted, then scored whatever stale file survived and appended a verdict row.
# Every input is therefore checked HERE, ABOVE the first destructive step; a missing one aborts nonzero.
abort() { echo "SWEEP ABORT: $*" >&2; exit 2; }

# Prevents the whole sweep running with a dead binary: _detect_unit.sh takes the catalog binary from
# $GWCAT and the repo-local target/ DOES NOT EXIST, so the mandated CARGO_TARGET_DIR build is the default.
# Exported so sweep, driver and every unit agree on the binary that is checked here. (§6am)
GWCAT=${GWCAT:-/mnt/linuxdisk/home/juanfraitu/rustle_target/release/gw_family_catalog}
export GWCAT
[ -x "$GWCAT" ] || abort "gw_family_catalog not executable at $GWCAT (build with CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target --release, or set GWCAT=)"

# Prevents "FASTA found too late": _detect_unit.sh resolves the CHM13 FASTA itself and exits 2 per unit,
# which happens only AFTER this script has already truncated its result file. Same candidates, same order,
# resolved before anything is destroyed. (§6am)
SOTO_FASTA=${SOTO_FASTA:-}
if [ -z "$SOTO_FASTA" ]; then
  for cand in /mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa \
              /mnt/c/Users/jfris/Desktop/Reference/chm13v2.0.fa; do
    if [ -f "$cand" ]; then SOTO_FASTA="$cand"; break; fi
  done
fi
[ -s "$SOTO_FASTA" ] || abort "no CHM13 FASTA found (set SOTO_FASTA=)"
export SOTO_FASTA

# Prevents a unit dying on mini-BAM extraction after the truncation. (§6am)
SAMTOOLS=${SAMTOOLS:-/home/juanfra/miniforge3/bin/samtools}
[ -x "$SAMTOOLS" ] || abort "samtools not executable at $SAMTOOLS (set SAMTOOLS=)"
export SAMTOOLS

# Prevents scoring against a missing/empty truth set or driving with a missing driver/scorer. (§6am)
[ -s "$BED" ]                        || abort "member BED missing or empty: $BED"
[ -s "$HERE/recompute_perchrom.sh" ] || abort "driver missing: $HERE/recompute_perchrom.sh"
[ -s "$HERE/soto_cache_score.py" ]   || abort "scorer missing: $HERE/soto_cache_score.py"
[ -d "$CACHE" ]                      || abort "soto cache dir missing: $CACHE"
[ -d "$PC/beds" ]                    || abort "per-chrom BEDs missing: $PC/beds"
NBED=$(find "$PC/beds" -maxdepth 1 -name '*.bed' | wc -l)
[ "$NBED" -gt 0 ]                    || abort "no per-chrom BEDs in $PC/beds"
[ -s "$PC/xchrom.bed" ]              || abort "cross-chrom BED missing or empty: $PC/xchrom.bed"

# Prevents a unit whose mini-BAM is NOT cached from failing on a missing source BAM; the source is only
# required when something still has to be extracted, so a cache-complete run is unaffected. (§6am)
if [ ! -s "$CACHE/soto_regions.bam" ]; then
  for b in "$PC"/beds/*.bed "$PC/xchrom.bed"; do
    n=$(basename "$b" .bed)
    [ -f "$PC/$n.bam.bai" ] || abort "$CACHE/soto_regions.bam missing and no cached mini-BAM for $n"
  done
fi

# Preserves the previous table and stamps the replacement with its provenance, so a legitimate re-sweep is
# distinguishable from a run that quietly overwrote a real measurement (§6am: "it deletes the evidence and
# writes a fabricated replacement in the same run"). Pass the justification in REASON=...
if [ -s "$SUM" ]; then cp -p "$SUM" "$SUM.prev.$(date +%Y%m%d%H%M%S)"; fi
{ echo "# swept $(date '+%F %T')  binary=$GWCAT  commit=$(git -C "$HERE" rev-parse --short HEAD 2>/dev/null || echo NA)"
  echo "# reason: ${REASON:-UNSTATED (a re-sweep without a stated reason is not auditable)}"
  echo -e "config\tgate_env\tcopies\ton_member%\trecall\tdetected"; } > $SUM

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
  # Freshness stamp, written BEFORE the rm + recompute so each unit output can be PROVED to be new. (§6am)
  local stamp="$CACHE/.gate_${label}.stamp"; : > "$stamp"
  rm -f "$PC"/chr*.copies.tsv "$PC"/xchrom.copies.tsv   # clear stale (keep cached mini-BAMs)
  env $genv bash "$HERE/recompute_perchrom.sh" 5 > "$CACHE/gate_${label}.log" 2>&1 \
    || { echo "SWEEP ABORT [$label]: recompute_perchrom.sh exited nonzero; see $CACHE/gate_${label}.log" >&2; exit 3; }
  # Backstop: recompute_perchrom.sh exits 0 even when EVERY unit dies (xargs swallows the unit rcs) and its
  # combine step globs *.copies.tsv, so a dead run can leave a header-only or stale catalog that still
  # scores. Nothing is scored unless every expected unit wrote a copies.tsv newer than the stamp. (§6am)
  local expected=$(( NBED + 1 )) fresh=0 missing="" b n out
  for b in "$PC"/beds/*.bed "$PC/xchrom.bed"; do
    n=$(basename "$b" .bed); out="$PC/$n.copies.tsv"
    if [ -f "$out" ] && [ "$out" -nt "$stamp" ]; then fresh=$((fresh+1)); else missing="$missing $n"; fi
  done
  [ "$fresh" -eq "$expected" ] \
    || { echo "SWEEP ABORT [$label]: only $fresh/$expected units rebuilt; missing:$missing (see $CACHE/gate_${label}.log)" >&2; exit 3; }
  local CAT="$CACHE/perchrom_catalog.copies.tsv"
  { [ -f "$CAT" ] && [ "$CAT" -nt "$stamp" ]; } \
    || { echo "SWEEP ABORT [$label]: $CAT was not rebuilt by this config" >&2; exit 3; }
  # An empty catalog scoring as a row is the same defect wearing a different hat. (§6am)
  local ncopy=$(( $(wc -l < "$CAT") - 1 ))
  [ "$ncopy" -gt 0 ] \
    || { echo "SWEEP ABORT [$label]: rebuilt catalog has 0 copies -- refusing to score an empty catalog" >&2; exit 3; }
  cp "$CACHE/perchrom_catalog.copies.tsv" "$CACHE/gate_${label}.copies.tsv"
  local cs; cs=$(score_copies "$CACHE/gate_${label}.copies.tsv") || true
  # No score / no recall = an empty evidence set; refuse to append a verdict row that measures nothing. (§6am)
  [ -n "$cs" ] || { echo "SWEEP ABORT [$label]: score_copies produced no output" >&2; exit 3; }
  local rec; rec=$(python3 "$HERE/soto_cache_score.py" "$CACHE/gate_${label}.copies.tsv" 2>/dev/null \
                    | grep -oE "NEW: [0-9]+/362 = [0-9.]+%" | sed 's/NEW: //') || true
  [ -n "$rec" ] || { echo "SWEEP ABORT [$label]: soto_cache_score.py emitted no 'NEW: n/362' recall line" >&2; exit 3; }
  echo -e "${label}\t${genv:--(default gates ON)}\t${cs}\t${rec}" >> $SUM
  echo "[$label] copies/on% = $cs   recall = $rec"
}

# configs: default (gates ON) + gates OFF + threshold variations
run_config default            ""
run_config gates_off          "RUSTLE_KEEP_READTHROUGH=1"
run_config readthrough_strict "RUSTLE_READTHROUGH_MIN_DISTINCT=3"
run_config mischain_aggressive "RUSTLE_MISCHAIN_GIANT_BP=20000 RUSTLE_MISCHAIN_MIN_READS=5"

echo "=== GATE SENSITIVITY SUMMARY ==="; cat $SUM
