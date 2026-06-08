#!/usr/bin/env bash
# Pinned configuration for the VG copy-recovery evaluation protocol.
# Source this from every stage script: `source "$(dirname "$0")/config.sh"`.
set -euo pipefail

# ── Data ────────────────────────────────────────────────────────────────────
SCRATCH="${SCRATCH:-/home/juanfra/winloci_scratch}"
BAM="${BAM:-$SCRATCH/GGO.bam}"
GENOME_FASTA="${GENOME_FASTA:-$SCRATCH/GGO.fasta}"
REF_GFF="${REF_GFF:-$SCRATCH/GGO_tx.gff}"

# ── Tools ───────────────────────────────────────────────────────────────────
ROOT="${ROOT:-/mnt/c/Users/jfris/Desktop/Rustle}"
RUSTLE="${RUSTLE:-$ROOT/target/release/rustle}"
STRINGTIE="${STRINGTIE:-$ROOT/tools/stringtie/stringtie}"
SAMTOOLS="${SAMTOOLS:-samtools}"
MINIMAP2="${MINIMAP2:-minimap2}"
GFFREAD="${GFFREAD:-/home/juanfra/miniforge3/envs/sqanti3/bin/gffread}"  # not on base PATH; sqanti3 env has 0.12.7
SQANTI3_DIR="${SQANTI3_DIR:-$HOME/tools/SQANTI3}"
SQANTI3_ENV="${SQANTI3_ENV:-sqanti3}"
# Python interpreter with pysam/pandas available (base env lacks them).
PY="${PY:-mamba run -n phasing_eval python}"

# ── Output ──────────────────────────────────────────────────────────────────
OUTDIR="${OUTDIR:-$ROOT/bench/copy_recovery_eval/results}"

# ── Arms (1=on) ─────────────────────────────────────────────────────────────
ARM_STRINGTIE=1                 # baseline (headline competitor)
ARM_RUSTLE_VG=1                 # headline method (decisive gate ON)
ARM_RUSTLE_VG_RAW="${ARM_RUSTLE_VG_RAW:-0}"      # optional ablation (gate OFF)
ARM_RUSTLE_PRIMARY="${ARM_RUSTLE_PRIMARY:-0}"    # optional attribution arm (-L primary-only)

# ── rustle --vg headline flags (win stack + decisive gate ON) ───────────────
RUSTLE_VG_ENV="RUSTLE_VG_TANDEM=1 RUSTLE_VG_TANDEM_PRIMARY_JUNCTIONS=1 RUSTLE_VG_DECISIVE_GATE=1 RUSTLE_VG_DECISIVE_GATE_MIN_PRIM=4"
RUSTLE_VG_RAW_ENV="RUSTLE_VG_TANDEM=1 RUSTLE_VG_TANDEM_PRIMARY_JUNCTIONS=1"   # no gate

# ── Thresholds (pinned) ─────────────────────────────────────────────────────
MIN_IDENTITY="${MIN_IDENTITY:-0.90}"   # U3 paralog identity floor
MIN_COV_FRAC="${MIN_COV_FRAC:-0.50}"   # U3 alignment length-coverage floor
GUARD_K="${GUARD_K:-1}"                # U5 min primary reads for authentic
# Restrict expensive runs to one chrom for smoke tests (empty = whole genome):
CHROM_SUBSET="${CHROM_SUBSET:-}"

# ── Provenance ──────────────────────────────────────────────────────────────
config_hash() {
  { echo "$MIN_IDENTITY $MIN_COV_FRAC $GUARD_K $RUSTLE_VG_ENV $RUSTLE_VG_RAW_ENV"; } | sha1sum | cut -d' ' -f1
}
mkdir -p "$OUTDIR"
