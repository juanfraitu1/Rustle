#!/bin/bash
# usage: bakeoff_chrom_ours.sh CHROM LABEL [extra copy_assign args...]
# Runs copy_assign --gtf over the whole chromosome into $W/LABEL/ours.gtf. Env vars this script is
# explicitly told to pass through (RUSTLE_LEGACY_PLACEMENT_DEDUP) survive from the caller; every OTHER
# RUSTLE_* variable is force-unset below.
#
# Hardening for the held-out chr17 test (docs/PREREG_gtf_refine_chr17_2026-09-16.md, "Run protocol (frozen)"):
# the PREREG names five --gtf-refine components and RUSTLE_LEGACY_PLACEMENT_DEDUP only. Several OTHER
# RUSTLE_* env vars exist in this codebase that can reach copy_assign's --gtf path or its dependencies
# (RUSTLE_TSS_SNAP[/_WINDOW/_FRAC], RUSTLE_FOOTPRINT_NODES, RUSTLE_JUNCTION_FUZZ_BP, RUSTLE_JUNCTION_MAJORITY,
# RUSTLE_JUNCTION_NC_MAX_BP, and any added later) and every one of them defaults OFF/byte-identical when
# unset -- so unsetting everything not named by the PREREG, unconditionally, guarantees the run matches the
# frozen rules regardless of whatever happens to be exported in the calling shell.
set -euo pipefail
CHROM=${1:?CHROM}; LABEL=${2:?LABEL}; shift 2
W=/mnt/linuxdisk/home/juanfraitu/bakeoff/human_${CHROM}
BIN=/mnt/linuxdisk/home/juanfraitu/rustle_target/release/copy_assign
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BASE_COMMIT=3b096012  # src/Cargo.* freeze target (the commit that froze these rules); NOT necessarily the
                       # PREREG doc's own latest commit -- see docs/PREREG_gtf_refine_chr17_2026-09-16.md
                       # "Frozen file versions"

if [ "$LABEL" = legacy ]; then
  if [ "${RUSTLE_LEGACY_PLACEMENT_DEDUP:-}" != "1" ]; then
    echo "REFUSING to run: label 'legacy' requires RUSTLE_LEGACY_PLACEMENT_DEDUP=1 (got '${RUSTLE_LEGACY_PLACEMENT_DEDUP:-unset}')." >&2
    exit 3
  fi
elif [ -n "${RUSTLE_LEGACY_PLACEMENT_DEDUP:-}" ]; then
  echo "REFUSING to run: label '$LABEL' requires RUSTLE_LEGACY_PLACEMENT_DEDUP to be unset (got '$RUSTLE_LEGACY_PLACEMENT_DEDUP')." >&2
  exit 3
fi

for v in $(compgen -e | grep '^RUSTLE_' || true); do
  [ "$v" = RUSTLE_LEGACY_PLACEMENT_DEDUP ] || unset "$v"
done

LEN=$(awk -v c="$CHROM" '$1==c {print $2}' "$W/${CHROM}.fa.fai")
mkdir -p "$W/$LABEL"; cd "$W/$LABEL"

# Working-tree check against BASE_COMMIT (not HEAD): catches uncommitted AND untracked src/Cargo.* edits,
# not just committed ones. Either check failing, or a git error from either command, -> SRC-CHANGED.
DIFF_STATUS=ok
git -C "$REPO_DIR" diff --quiet "$BASE_COMMIT" -- src Cargo.toml Cargo.lock || DIFF_STATUS=changed
STATUS_OUT=$(git -C "$REPO_DIR" status --porcelain -- src Cargo.toml Cargo.lock 2>&1) || DIFF_STATUS=changed
[ -z "$STATUS_OUT" ] || DIFF_STATUS=changed
if [ "$DIFF_STATUS" = ok ]; then
  SRC_STATUS=src-identical
else
  SRC_STATUS=SRC-CHANGED
fi
{
  echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "repo_head: $(git -C "$REPO_DIR" rev-parse HEAD)"
  echo "prereg_base_commit: $BASE_COMMIT"
  echo "src_check (src, Cargo.toml, Cargo.lock vs prereg_base_commit): $SRC_STATUS"
  echo "binary: $BIN"
  echo "binary_sha256: $(sha256sum "$BIN" | awk '{print $1}')"
  echo "env_RUSTLE_star:"
  env | grep '^RUSTLE_' || echo "  (none set)"
  echo "command: $BIN --gtf --bam $W/${CHROM}.bam --fasta $W/${CHROM}.fa --region ${CHROM}:1-${LEN} $* --out ours"
} > run_provenance.txt
cat run_provenance.txt

if [ "$SRC_STATUS" = "SRC-CHANGED" ]; then
  echo "REFUSING to run: src/Cargo.toml/Cargo.lock differ from the PREREG base commit $BASE_COMMIT." >&2
  echo "Rebuild the binary from a tree at $BASE_COMMIT (or a later commit with NO src/Cargo.* changes) before running the held-out chr17 test." >&2
  exit 3
fi

"$BIN" --gtf --bam "$W/${CHROM}.bam" --fasta "$W/${CHROM}.fa" --region "${CHROM}:1-${LEN}" "$@" --out ours \
  > ours.stdout.log 2> ours.stderr.log
