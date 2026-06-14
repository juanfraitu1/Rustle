#!/usr/bin/env bash
# make_ghost_fixture.sh — GHOST-copy paralog fixture for M5 new-copy recovery tests.
#
# Produces bench/fixtures/sim_ghost.{fa,fa.fai,reads.fa,bam,bam.bai}
#
# Design (see bench/sim/gen_ghost_fixture.py for full detail):
#   chrSIM_GHOST, two 3-exon gene copies at [10k-12k] (copy A) and [30k-32k] (copy B).
#   ~99.3% exon identity (6 PSVs, 2 per exon).
#
#   Read mix producing the GHOST condition (copy B has ZERO primaries):
#     copyA_full_*   40 reads — primary at A
#     copyA_skip_*   10 reads — primary at A (exon1+exon3 skip isoform)
#     copyB_xmap_*    8 reads — copy-A sequence → primary at A, SECONDARY at B
#
#   Result:
#     copy A: 50 primaries → fully assembled
#     copy B:  0 primaries → NO bundle forms (Layer 1 builds nothing at B)
#              8 secondaries → locus=None; ONLY M5 (--vg-layer2-new-copies) recovers B
#
# Dependencies (must be on PATH):
#   minimap2 >= 2.28   samtools >= 1.19   python3 >= 3.8
#
# Usage:
#   cd <repo-root>
#   bash bench/sim/make_ghost_fixture.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
OUT_DIR="${REPO_ROOT}/bench/fixtures"

echo "[make_ghost_fixture] Generating GHOST-copy paralog fixture..."
echo "[make_ghost_fixture] Output directory: ${OUT_DIR}"
echo "[make_ghost_fixture] minimap2: $(minimap2 --version)"
echo "[make_ghost_fixture] samtools: $(samtools --version | head -1)"

python3 "${SCRIPT_DIR}/gen_ghost_fixture.py" --out-dir "${OUT_DIR}"

echo ""
echo "[make_ghost_fixture] Done.  Key files:"
ls -lh \
  "${OUT_DIR}/sim_ghost.fa" \
  "${OUT_DIR}/sim_ghost.fa.fai" \
  "${OUT_DIR}/sim_ghost.reads.fa" \
  "${OUT_DIR}/sim_ghost.bam" \
  "${OUT_DIR}/sim_ghost.bam.bai"
