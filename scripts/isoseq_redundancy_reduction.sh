#!/usr/bin/env bash
# Reduce Iso-Seq read redundancy using PacBio smrttools (same roles as isoseq cluster + collapse).
#
# Two stages (conceptually):
#   1) Cluster  — isoseq cluster2 on **FLNC** (unmapped): one polished transcript per cluster.
#   2) Collapse — isoseq collapse on **genome-aligned** BAM: merge identical / fuzzy splice isoforms.
#
# Your minimap2 BAM (e.g. GGO.bam) skips (1). You can still run (2) to collapse redundant alignments.
# For the full PacBio stack, run cluster2 on FLNC, align with pbmm2 --preset ISOSEQ, then collapse.
#
# Requires: isoseq (smrttools), and for cluster/align: pbmm2, samtools. Optional: pbindex for FLNC.
#
# Usage:
#   ./scripts/isoseq_redundancy_reduction.sh collapse-mapped [--threads N] <mapped.bam> <out_prefix>
#       Writes <out_prefix>.gff (+ .abundance.txt, .read_stat.txt, ...). Uses bulk-friendly flags.
#
#   ./scripts/isoseq_redundancy_reduction.sh cluster-then-collapse \\
#       [--threads N] <ref.fa> <flnc.bam> <out_prefix>
#       cluster2 FLNC → transcripts.bam → pbmm2 ISOSEQ → collapse (3-arg; indexes FLNC with pbindex if needed).
#
# Examples:
#   ./scripts/isoseq_redundancy_reduction.sh collapse-mapped GGO.bam GGO_isoseq_collapsed
#   ./scripts/isoseq_redundancy_reduction.sh cluster-then-collapse /path/to/mGorGor.fa movie.flnc.bam GGO_from_flnc

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

THREADS=0
CMD="${1:?usage: $0 collapse-mapped|cluster-then-collapse [options] ...}"
shift || true

while [[ $# -gt 0 ]]; do
  case "$1" in
    --threads)
      THREADS="${2:?}"
      shift 2
      ;;
    *)
      break
      ;;
  esac
done

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "error: '$1' not in PATH (install smrttools / minimap2 stack)" >&2
    exit 1
  }
}

collapse_mapped() {
  local mapped="$1"
  local prefix="$2"
  need_cmd isoseq
  [[ -f "$mapped" ]] || { echo "error: not a file: $mapped" >&2; exit 1; }

  local out_gff="${prefix}.gff"
  local thread_arg=()
  if [[ "$THREADS" -gt 0 ]]; then
    thread_arg=( -j "$THREADS" )
  fi

  echo "Running: isoseq collapse (mapped-only) → ${out_gff}"
  echo "  (Warns if im/is tags missing — each read treated as one molecule.)"
  # Two-argument form: alignments + output GFF (FLNC optional in some builds; works with minimap2 BAM here).
  isoseq collapse "${thread_arg[@]}" \
    --do-not-collapse-extra-5exons \
    "$mapped" \
    "$out_gff"

  echo "Done. Sidecars: ${prefix}.*"
}

cluster_then_collapse() {
  local ref_fa="$1"
  local flnc_bam="$2"
  local prefix="$3"
  need_cmd isoseq
  need_cmd pbmm2
  need_cmd pbindex

  [[ -f "$ref_fa" && -f "$flnc_bam" ]] || {
    echo "error: need ref FASTA and FLNC BAM" >&2
    exit 1
  }

  local transcripts="${prefix}.transcripts.bam"
  local mapped="${prefix}.pbmm2.mapped.bam"
  local out_gff="${prefix}.collapsed.gff"

  local thread_arg=()
  if [[ "$THREADS" -gt 0 ]]; then
    thread_arg=( -j "$THREADS" )
  fi

  echo "Step 1/3: isoseq cluster2 (FLNC → polished transcripts)"
  isoseq cluster2 "${thread_arg[@]}" "$flnc_bam" "$transcripts"

  echo "Step 2/3: pbmm2 align --preset ISOSEQ (recommended for collapse)"
  pbmm2 align --preset ISOSEQ --sort -N 1 "$ref_fa" "$transcripts" "$mapped"

  if [[ ! -f "${flnc_bam}.pbi" ]]; then
    if command -v pbindex >/dev/null 2>&1; then
      echo "Indexing FLNC for collapse (pbindex)…"
      pbindex "$flnc_bam"
    else
      echo "error: ${flnc_bam}.pbi missing and pbindex not found; install smrttools pbindex or index FLNC." >&2
      exit 1
    fi
  fi

  echo "Step 3/3: isoseq collapse (bulk-style, 3-arg)"
  isoseq collapse "${thread_arg[@]}" \
    --do-not-collapse-extra-5exons \
    "$mapped" \
    "$flnc_bam" \
    "$out_gff"

  echo "Done: $out_gff (and ${prefix}.collapsed.*)"
}

case "$CMD" in
  collapse-mapped)
    [[ $# -eq 2 ]] || { echo "usage: $0 collapse-mapped [--threads N] <mapped.bam> <out_prefix>" >&2; exit 1; }
    collapse_mapped "$1" "$2"
    ;;
  cluster-then-collapse)
    [[ $# -eq 3 ]] || { echo "usage: $0 cluster-then-collapse [--threads N] <ref.fa> <flnc.bam> <out_prefix>" >&2; exit 1; }
    cluster_then_collapse "$1" "$2" "$3"
    ;;
  -h|--help|help)
    sed -n '1,35p' "$0"
    ;;
  *)
    echo "unknown subcommand: $CMD" >&2
    exit 1
    ;;
esac
