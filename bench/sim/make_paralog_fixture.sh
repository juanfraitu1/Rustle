#!/usr/bin/env bash
# make_paralog_fixture.sh — Deterministic 2-copy paralog fixture for Layer-2 recovery tests.
#
# Tool choice
# -----------
# We do NOT use pbsim3 / badread.  HiFi reads are ~99.9% accurate, so emitting
# exact spliced transcript sequences (0 errors) is a valid and simpler approach.
# A small Python script (gen_paralog_fixture.py, same directory) constructs the
# genome and reads directly, then calls minimap2 + samtools.
#
# Dependencies (must be on PATH):
#   minimap2   >= 2.28
#   samtools   >= 1.19
#   python3    >= 3.8
#
# Usage:
#   cd <repo-root>
#   bash bench/sim/make_paralog_fixture.sh
#
# All outputs are written to bench/fixtures/.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
OUT_DIR="${REPO_ROOT}/bench/fixtures"

echo "[make_paralog_fixture] Generating synthetic 2-copy paralog fixture..."
echo "[make_paralog_fixture] Output directory: ${OUT_DIR}"
echo "[make_paralog_fixture] minimap2: $(minimap2 --version)"
echo "[make_paralog_fixture] samtools: $(samtools --version | head -1)"

python3 "${SCRIPT_DIR}/gen_paralog_fixture.py" --out-dir "${OUT_DIR}"

echo ""
echo "[make_paralog_fixture] Done.  Key files:"
ls -lh \
  "${OUT_DIR}/sim_paralog.fa" \
  "${OUT_DIR}/sim_paralog.fa.fai" \
  "${OUT_DIR}/sim_paralog.reads.fa" \
  "${OUT_DIR}/sim_paralog.bam" \
  "${OUT_DIR}/sim_paralog.bam.bai"
