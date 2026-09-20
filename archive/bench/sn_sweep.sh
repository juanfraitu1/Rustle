#!/usr/bin/env bash
# Sensitivity sweep: loosen filter levers to recover missed reference transcripts.
# Each experiment is a named rustle+gffcompare run on GGO_19 de novo (-L).
# Results are collected into a Pareto table at the end.
#
# Usage: bash bench/sn_sweep.sh [config_pattern]
#   config_pattern: optional grep filter to run only matching configs
#
# All runs use -p 8. Sequential, ~3-5 min each, ~1.5 hr total for all 20 experiments.
set -euo pipefail

RUSTLE=./target/release/rustle
BAM=GGO_19.bam
REF=/mnt/c/Users/jfris/Desktop/GGO_19.gtf
OUTDIR=bench/sn_sweep
FILTER="${1:-}"

mkdir -p "$OUTDIR"

run_experiment() {
    local name="$1"
    shift
    local -a args=("$@")

    # Apply filter if specified
    if [[ -n "$FILTER" && "$name" != *"$FILTER"* ]]; then
        return 0
    fi

    local gtf="$OUTDIR/${name}.gtf"
    local cmp="$OUTDIR/${name}_cmp"

    # gffcompare sometimes drops the .stats suffix when the prefix contains dots
    # (e.g. f_0.008_cmp → creates f_0.008_cmp instead of f_0.008_cmp.stats).
    # Accept either naming.
    find_stats() {
        local p="$1"
        if   [[ -f "${p}.stats" ]]; then echo "${p}.stats"
        elif [[ -f "${p}" ]] && grep -q "Intron chain" "${p}" 2>/dev/null; then echo "${p}"
        fi
    }

    if [[ -n "$(find_stats "$cmp")" ]]; then
        echo "[SKIP] $name (already done)"
        return 0
    fi

    echo ""
    echo "=== $name ==="
    echo "  args: ${args[*]}"
    echo "  start: $(date '+%H:%M:%S')"

    # Extract env vars from args (prefix VAR=val)
    local env_args=()
    local rustle_args=(-L -p 8 -o "$gtf")
    for arg in "${args[@]}"; do
        if [[ "$arg" == *=* && "$arg" != --* && "$arg" != -* ]]; then
            env_args+=("$arg")
        else
            rustle_args+=("$arg")
        fi
    done

    env "${env_args[@]}" "$RUSTLE" "$BAM" "${rustle_args[@]}" 2>"$OUTDIR/${name}.stderr"
    gffcompare -r "$REF" -o "$cmp" "$gtf"

    local stats
    stats="$(find_stats "$cmp")"
    echo "  done:  $(date '+%H:%M:%S')"
    grep "Intron chain level:" "$stats" | sed "s/^/  IC: /"
}

# ─── Phase 1: isofrac (-f) sweep ────────────────────────────────────────────
# The predcluster sweep-line threshold.  Default = 0.01 (1%).
# Lower = more minority isoforms survive.

run_experiment "baseline"          # all defaults: -f 0.01, INTRA=0.01

run_experiment "f_0.008"   -f 0.008
run_experiment "f_0.005"   -f 0.005
run_experiment "f_0.003"   -f 0.003
run_experiment "f_0.001"   -f 0.001
run_experiment "f_0.0003"  -f 0.0003
run_experiment "f_off"     -f 0

# ─── Phase 2: global intra-gene coverage filter sweep ───────────────────────
# Applied genome-wide after all predcluster filters.  Default = 0.01 (1%).

run_experiment "intra_0.005"  "RUSTLE_INTRA_GENE_COV_FRAC=0.005"
run_experiment "intra_0.003"  "RUSTLE_INTRA_GENE_COV_FRAC=0.003"
run_experiment "intra_0.001"  "RUSTLE_INTRA_GENE_COV_FRAC=0.001"
run_experiment "intra_off"    "RUSTLE_INTRA_GENE_COV_FRAC_OFF=1"

# ─── Phase 3: pairwise (retained-intron) filter ─────────────────────────────
# --max-sensitivity disables the entire pairwise containment/RI pass.

run_experiment "maxsens"  --max-sensitivity true

# ─── Phase 4: junction mm rescue ────────────────────────────────────────────
# RUSTLE_JCT_ISOFRAC_RESCUE=1 rescues minority paths that have well-supported junctions.
# Previously: +6 TPs, +12 FPs at default -f.

run_experiment "jct_rescue"        "RUSTLE_JCT_ISOFRAC_RESCUE=1"
run_experiment "f_0.005_jct"       "RUSTLE_JCT_ISOFRAC_RESCUE=1"  -f 0.005
run_experiment "f_0.003_jct"       "RUSTLE_JCT_ISOFRAC_RESCUE=1"  -f 0.003

# ─── Phase 5: combined loose configurations ──────────────────────────────────
# Combine multiple levers to see if the gains compound.

run_experiment "f_0.003_intra_0.003"  \
    "RUSTLE_INTRA_GENE_COV_FRAC=0.003"  -f 0.003

run_experiment "f_0.001_intra_0.001"  \
    "RUSTLE_INTRA_GENE_COV_FRAC=0.001"  -f 0.001

run_experiment "f_0.003_maxsens"  \
    "RUSTLE_INTRA_GENE_COV_FRAC=0.003"  -f 0.003  --max-sensitivity true

run_experiment "all_loose"  \
    "RUSTLE_JCT_ISOFRAC_RESCUE=1"  "RUSTLE_INTRA_GENE_COV_FRAC=0.001"  \
    -f 0.001  --max-sensitivity true

# ─── Summary table ───────────────────────────────────────────────────────────
echo ""
echo "=== RESULTS ==="
printf "%-30s  %6s  %6s  %6s  %7s  %s\n" "config" "IC_Sn" "IC_Pr" "IC_F1" "chains" "tx_count"
printf "%-30s  %6s  %6s  %6s  %7s  %s\n" "------" "-----" "-----" "-----" "------" "--------"

for cmp_prefix in "$OUTDIR"/*_cmp "$OUTDIR"/*_cmp.stats; do
    # Deduplicate: skip if this is a .stats variant we already processed as bare
    [[ -f "$cmp_prefix" ]] || continue
    grep -q "Intron chain" "$cmp_prefix" 2>/dev/null || continue

    local_base=$(basename "$cmp_prefix")
    # Strip .stats suffix if present, then strip _cmp
    local_base="${local_base%.stats}"
    name="${local_base%_cmp}"
    stats="$cmp_prefix"

    sn=$(  grep "Intron chain level:" "$stats" | awk '{print $4}')
    pr=$(  grep "Intron chain level:" "$stats" | awk '{print $6}')
    chains=$(grep "Matching intron chains:" "$stats" | awk '{print $NF}')
    tx=$(  grep "Matching transcripts:"   "$stats" | awk '{print $NF}')

    # F1
    f1=$(awk -v sn="$sn" -v pr="$pr" \
        'BEGIN { if (sn+pr>0) printf "%.1f", 2*sn*pr/(sn+pr); else print "0.0" }')

    printf "%-30s  %6s  %6s  %6s  %7s  %s\n" \
        "$name" "$sn" "$pr" "$f1" "$chains" "$tx"
done | sort -k1,1
