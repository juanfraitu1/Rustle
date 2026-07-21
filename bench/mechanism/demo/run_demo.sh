#!/usr/bin/env bash
# V4c identity-gradient frontier demo — runnable driver.
#
# Regenerates the whole ladder + recovery anchor from scratch: for each target pairwise
# identity, builds a synthetic 4-copy tandem gene family at that divergence, simulates
# HiFi-like reads, aligns them, and runs the real `copy_assign` binary. Everything is
# seeded (Python's random.Random only — no wall-clock), so re-running reproduces the same
# qualitative story every time (see README.md for the small amount of expected read-level
# jitter from minimap2's own threading, disclosed there).
#
# Crash-rule compliant (see repo CLAUDE.md / advisor notes): every step below runs in the
# FOREGROUND, one at a time, serially. No `&`, no `nohup`, no background jobs, no waiter
# loops, no `pkill`. Output goes to a scratch dir under $HOME, never /tmp.
#
# Usage:
#   ./run_demo.sh          # full 6-point ladder (~90/94/96/98/99.5/100% identity)
#   ./run_demo.sh --fast   # 4-point representative subset (94/98/99.5/100%) — faster,
#                           # same qualitative story
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$HERE/../../.." && pwd)"
SCRATCH="${RUSTLE_MECH_DEMO_SCRATCH:-$HOME/winloci_scratch/mech_demo}"

if [[ ! -x "$REPO_ROOT/target/release/copy_assign" ]]; then
    echo "error: $REPO_ROOT/target/release/copy_assign not built." >&2
    echo "  Build it first: (cd $REPO_ROOT && cargo build --release --bin copy_assign)" >&2
    exit 1
fi

if ! command -v minimap2 >/dev/null 2>&1 && [[ ! -x "$HOME/miniforge3/bin/minimap2" ]]; then
    echo "error: minimap2 not found on PATH or at ~/miniforge3/bin/minimap2." >&2
    exit 1
fi

mkdir -p "$SCRATCH"

echo "== V4c demo: running foreground, serially, output -> $SCRATCH ==" >&2
python3 "$REPO_ROOT/bench/mechanism/sim_ladder.py" \
    --scratch "$SCRATCH" \
    --out "$REPO_ROOT/bench/mechanism/verification_results.json" \
    --demo-dir "$HERE" \
    "$@"

echo "== Done. Results written to: ==" >&2
echo "  $REPO_ROOT/bench/mechanism/verification_results.json (key \"v4c\")" >&2
echo "  $HERE/frontier.svg" >&2
echo "  $HERE/ref.fa, $HERE/reads.fq (regenerated example inputs, most-divergent ladder point)" >&2
