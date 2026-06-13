#!/usr/bin/env bash
# make_starved_fixture.sh — Genuine-starvation paralog fixture for Layer-2 recovery tests.
#
# Produces bench/fixtures/sim_starved.{fa,fa.fai,reads.fa,bam,bam.bai}
#
# Design (see bench/sim/gen_starved_fixture.py for full detail):
#   chrSIM_STARVED, two 3-exon gene copies at [10k-12k] (copy A) and [30k-32k] (copy B).
#   ~99.3% exon identity (6 PSVs, 2 per exon).
#
#   Read mix producing genuine starvation:
#     copyA_full_*   40 reads — primary at A
#     copyA_skip_*   10 reads — primary at A (exon1+exon3 skip isoform)
#     copyB_primary_* 1 read  — exact copy-B sequence → primary at B  (bundle survives)
#     copyB_xmap_*   20 reads — copy-A sequence → primary at A        (secondary at B)
#
#   Result:
#     copy A: 70 primaries → fully assembled (2 isoforms)
#     copy B:  1 primary  → bundle forms but NO transcript in default --vg output
#             70 secondaries → side-index carries B's splice coordinates
#
# Dependencies (must be on PATH):
#   minimap2 >= 2.28   samtools >= 1.19   python3 >= 3.8
#
# Usage:
#   cd <repo-root>
#   bash bench/sim/make_starved_fixture.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
OUT_DIR="${REPO_ROOT}/bench/fixtures"

echo "[make_starved_fixture] Generating genuine-starvation paralog fixture..."
echo "[make_starved_fixture] Output directory: ${OUT_DIR}"
echo "[make_starved_fixture] minimap2: $(minimap2 --version)"
echo "[make_starved_fixture] samtools: $(samtools --version | head -1)"

python3 "${SCRIPT_DIR}/gen_starved_fixture.py" --out-dir "${OUT_DIR}"

echo ""
echo "[make_starved_fixture] Done.  Key files:"
ls -lh \
  "${OUT_DIR}/sim_starved.fa" \
  "${OUT_DIR}/sim_starved.fa.fai" \
  "${OUT_DIR}/sim_starved.reads.fa" \
  "${OUT_DIR}/sim_starved.bam" \
  "${OUT_DIR}/sim_starved.bam.bai"
