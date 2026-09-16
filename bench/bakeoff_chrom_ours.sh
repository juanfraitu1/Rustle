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
#
# Fix round 3 (finding A): BIN is now a PRIVATE copy made by the PREREG's frozen build step
# (/mnt/linuxdisk/home/juanfraitu/bakeoff/chr17_bin/copy_assign), not the shared rustle_target/release path --
# the shared path is workspace-relative and can be silently overwritten by an unrelated build in another
# checkout between the freeze-build and this run. This script refuses to run unless the private copy's
# CURRENT sha256 matches the sha256 recorded right after that build.
set -euo pipefail
CHROM=${1:?CHROM}; LABEL=${2:?LABEL}; shift 2
W=/mnt/linuxdisk/home/juanfraitu/bakeoff/human_${CHROM}
BIN=/mnt/linuxdisk/home/juanfraitu/bakeoff/chr17_bin/copy_assign
BIN_SHA_FILE=/mnt/linuxdisk/home/juanfraitu/bakeoff/chr17_bin/copy_assign.sha256
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BASE_COMMIT=3b096012  # src/Cargo.*/.cargo/config.toml freeze target (the commit that froze these rules); NOT
                       # necessarily the PREREG doc's own latest commit -- see
                       # docs/PREREG_gtf_refine_chr17_2026-09-16.md "Frozen file versions"

if [ "$LABEL" = legacy ]; then
  if [ "${RUSTLE_LEGACY_PLACEMENT_DEDUP:-}" != "1" ]; then
    echo "REFUSING to run: label 'legacy' requires RUSTLE_LEGACY_PLACEMENT_DEDUP=1 (got '${RUSTLE_LEGACY_PLACEMENT_DEDUP:-unset}')." >&2
    exit 3
  fi
elif [ -n "${RUSTLE_LEGACY_PLACEMENT_DEDUP+x}" ]; then
  # +x (not :-/-n) so a SET-BUT-EMPTY value is also refused, matching the PREREG text "set to anything at
  # all" (fix round 3, legacy guard nit).
  echo "REFUSING to run: label '$LABEL' requires RUSTLE_LEGACY_PLACEMENT_DEDUP to be unset (got '${RUSTLE_LEGACY_PLACEMENT_DEDUP}')." >&2
  exit 3
fi

for v in $(compgen -e | grep '^RUSTLE_' || true); do
  [ "$v" = RUSTLE_LEGACY_PLACEMENT_DEDUP ] || unset "$v"
done

LEN=$(awk -v c="$CHROM" '$1==c {print $2}' "$W/${CHROM}.fa.fai")
mkdir -p "$W/$LABEL"; cd "$W/$LABEL"

# Check 1: the SOURCE working tree against BASE_COMMIT (not HEAD, so uncommitted AND untracked edits are
# caught). .cargo/config.toml is included because it can change compiled behavior (e.g. rustflags) without
# touching src/ (fix round 3, cheap hardening). Either check failing, or a git error from either command,
# -> SRC-CHANGED.
DIFF_STATUS=ok
git -C "$REPO_DIR" diff --quiet "$BASE_COMMIT" -- src Cargo.toml Cargo.lock .cargo/config.toml || DIFF_STATUS=changed
STATUS_OUT=$(git -C "$REPO_DIR" status --porcelain -- src Cargo.toml Cargo.lock .cargo/config.toml 2>&1) || DIFF_STATUS=changed
[ -z "$STATUS_OUT" ] || DIFF_STATUS=changed
if [ "$DIFF_STATUS" = ok ]; then
  SRC_STATUS=src-identical
else
  SRC_STATUS=SRC-CHANGED
fi

# Check 2 (fix round 3, cheap hardening): the RUN-PROTOCOL files (this PREREG, the five bench scripts, the
# verdict script) have no fixed base commit to diff against -- they are allowed to sit at a LATER commit than
# BASE_COMMIT by design (fix rounds edit them after 3b096012). There is no way to diff against "the commit
# you are about to make", so instead require NO UNCOMMITTED edits to them at run time; whatever IS committed
# is governed by the PREREG's failure policy (an edit made after chr17 data exists voids the test regardless
# of this check).
PROTOCOL_FILES=(docs/PREREG_gtf_refine_chr17_2026-09-16.md bench/prep_chrom_ref.sh bench/bakeoff_chrom_ours.sh bench/bakeoff_chrom_stringtie.sh bench/bakeoff_chrom_flair.sh bench/chrom_score.sh bench/gtf_refine_verdict.py)
PROTOCOL_STATUS_OUT=$(git -C "$REPO_DIR" status --porcelain -- "${PROTOCOL_FILES[@]}" 2>&1) || PROTOCOL_STATUS_OUT="GIT-ERROR"
if [ -n "$PROTOCOL_STATUS_OUT" ]; then
  PROTOCOL_STATUS=PROTOCOL-CHANGED
else
  PROTOCOL_STATUS=protocol-clean
fi

# Check 3 (fix round 3, finding A): the binary is the PINNED PRIVATE COPY, matched against the sha256
# recorded right after the frozen build step -- not merely "some file exists at this path".
if [ ! -f "$BIN_SHA_FILE" ]; then
  BIN_STATUS=BIN-UNPINNED
  PINNED_SHA=""
  CURRENT_SHA=""
elif [ ! -f "$BIN" ]; then
  BIN_STATUS=BIN-MISSING
  PINNED_SHA=$(cat "$BIN_SHA_FILE")
  CURRENT_SHA=""
else
  PINNED_SHA=$(cat "$BIN_SHA_FILE")
  CURRENT_SHA=$(sha256sum "$BIN" | awk '{print $1}')
  if [ "$CURRENT_SHA" = "$PINNED_SHA" ]; then
    BIN_STATUS=bin-pinned
  else
    BIN_STATUS=BIN-CHANGED
  fi
fi

{
  echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "repo_head: $(git -C "$REPO_DIR" rev-parse HEAD)"
  echo "prereg_base_commit: $BASE_COMMIT"
  echo "src_check (src, Cargo.toml, Cargo.lock, .cargo/config.toml vs prereg_base_commit): $SRC_STATUS"
  echo "protocol_check (PREREG + 5 scripts + verdict script, uncommitted-edit check): $PROTOCOL_STATUS"
  echo "binary: $BIN"
  echo "binary_sha256: $CURRENT_SHA"
  echo "pinned_sha256 ($BIN_SHA_FILE): $PINNED_SHA"
  echo "bin_check: $BIN_STATUS"
  echo "env_RUSTLE_star:"
  env | grep '^RUSTLE_' || echo "  (none set)"
  echo "command: $BIN --gtf --bam $W/${CHROM}.bam --fasta $W/${CHROM}.fa --region ${CHROM}:1-${LEN} $* --out ours"
} > run_provenance.txt
cat run_provenance.txt

if [ "$SRC_STATUS" = "SRC-CHANGED" ]; then
  echo "REFUSING to run: src/Cargo.toml/Cargo.lock/.cargo/config.toml differ from the PREREG base commit $BASE_COMMIT." >&2
  echo "Rebuild the binary from a tree at $BASE_COMMIT (or a later commit with NO src/Cargo.*/.cargo/config.toml changes) before running the held-out chr17 test." >&2
  exit 3
fi
if [ "$PROTOCOL_STATUS" = "PROTOCOL-CHANGED" ]; then
  echo "REFUSING to run: uncommitted edits to the PREREG or the run-protocol scripts (see $PROTOCOL_STATUS_OUT)." >&2
  echo "Commit (or revert) those edits before running the held-out chr17 test." >&2
  exit 3
fi
if [ "$BIN_STATUS" != "bin-pinned" ]; then
  echo "REFUSING to run: binary not pinned ($BIN_STATUS). Run the PREREG's frozen build step to produce $BIN and $BIN_SHA_FILE." >&2
  exit 3
fi

"$BIN" --gtf --bam "$W/${CHROM}.bam" --fasta "$W/${CHROM}.fa" --region "${CHROM}:1-${LEN}" "$@" --out ours \
  > ours.stdout.log 2> ours.stderr.log
