#!/usr/bin/env bash
# Genome-wide downstream on/off for the contiguous-core family-merge gate (the deferred
# validation from core_gate_pipeline.md). Per-chrom SERIAL (OOM-safe, --vg is heavy), gate
# OFF (default) vs ON @0.13, coordinate-keyed diff (bench/gate_keyed_diff.py) so rayon
# line-order noise is ignored. Confirms genome-wide whether flipping the gate ON loses any
# real transcript (byte-identical-safety prerequisite for a default flip).
# Output: /tmp/gw/onoff_summary.tsv (one row/contig) + /tmp/gw/onoff_$C.diff (details).
#
# 2026-09-01 REPAIR (docs/o1_ledger.md §6am — this instrument was a CONFIRMED FALSE PASS):
# $RUSTLE pointed at the repo-local target/release/rustle, which DOES NOT EXIST, and the return
# code was never captured. The GTFs are written by rustle's -o (not a shell redirect), so a silent
# exec failure left a PRIOR run's off_/on_ GTFs in place and the keyed diff reported lost=0
# gained=0 — the strongest claim this script can make ("flipping the gate ON loses no real
# transcript") emitted from ZERO execution — and the per-contig `touch` then latched it as done.
# The science below is UNCHANGED; only guards were added.
set -euo pipefail
BAM=/home/juanfra/winloci_scratch/GGO.bam
FASTA=/home/juanfra/winloci_scratch/GGO.fasta
FAI=/home/juanfra/winloci_scratch/GGO.fasta.fai
# The mandated build dir is CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target; the
# in-repo target/ does not exist. Overridable for A/B-ing two binaries.
RUSTLE=${RUSTLE:-/mnt/linuxdisk/home/juanfraitu/rustle_target/release/rustle}
SAM=/home/juanfra/miniforge3/bin/samtools
PY=/home/juanfra/miniforge3/bin/python
DIFF=/mnt/c/Users/jfris/Desktop/Rustle/bench/gate_keyed_diff.py
OUT=/tmp/gw
SUM=$OUT/onoff_summary.tsv

# --- §6am guards: nothing may run, and $SUM may not be truncated, until every binary and input is
# --- present. A missing binary here is exactly what let the diff score a prior run's leftovers.
for x in "$RUSTLE" "$SAM" "$PY"; do
  [[ -x "$x" ]] || { echo "ABORT: required executable missing or not executable: $x" >&2
                     echo "  (build with CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target --release; override with RUSTLE=/path/to/rustle)" >&2
                     exit 2; }
done
# §6am: a missing/empty input must abort BEFORE the summary is truncated, not degrade into a zero row.
for f in "$BAM" "$FASTA" "$FAI" "$DIFF"; do
  [[ -s "$f" ]] || { echo "ABORT: required input missing or empty: $f" >&2; exit 2; }
done
# §6am: without a BAM index every per-contig `samtools view` fails and every contig silently skips.
[[ -s "$BAM.bai" || -s "$BAM.csi" ]] || { echo "ABORT: no index for $BAM (.bai/.csi)" >&2; exit 2; }
# §6am: the per-contig guide GTFs ARE the evidence set — with none present every contig is skipped
# and the aggregate would still print lost=0/gained=0 over nothing.
[[ -d "$OUT" ]] || { echo "ABORT: work dir $OUT does not exist; it must already hold the per-contig guide GTFs $OUT/st_<contig>.gtf" >&2; exit 2; }
ls "$OUT"/st_*.gtf >/dev/null 2>&1 || { echo "ABORT: no guide GTFs $OUT/st_*.gtf in $OUT — nothing to score" >&2; exit 2; }

FAILED=0
# §6am: one exit path for a contig whose run failed — count it, clean its temporaries, write NO row
# and create NO .done latch, so the contig is re-attempted instead of being frozen as a pass.
fail_contig() { # $1=contig  $2=why
  echo "[$1] FAIL: $2 — no row written, contig NOT latched" >&2
  FAILED=$((FAILED + 1))
  rm -f "$OUT/$1.bam" "$OUT/$1.bam.bai" "$OUT/$1.fa" "$OUT/$1.fa.fai"
}

echo -e "chrom\ttx_off\ttx_on\tlost\tgained\treal_attr\tcosmetic\tgated" > "$SUM"
CONTIGS=$(awk '$2 > 1000000 {print $1}' "$FAI")
T0=$(date +%s)
for C in $CONTIGS; do
  [[ -f "$OUT/onoff_$C.done" ]] && { echo "[$C] done (skip)"; continue; }
  [[ -s "$OUT/st_$C.gtf" ]] || { echo "[$C] no guide, skip"; continue; }
  CS=$(date +%s)
  # A per-contig extraction failure keeps the original behaviour: it falls through to the RN==0 skip
  # below. A wholesale failure now surfaces at the empty-summary guard at the end (§6am).
  "$SAM" view -b "$BAM" "$C" -o "$OUT/$C.bam" 2>/dev/null || true
  RN=$("$SAM" view -c "$OUT/$C.bam" 2>/dev/null || echo 0)
  if [[ "$RN" -eq 0 ]]; then echo "[$C] 0 reads, skip"; rm -f "$OUT/$C.bam"; continue; fi
  "$SAM" index "$OUT/$C.bam"
  "$SAM" faidx "$FASTA" "$C" > "$OUT/$C.fa" 2>/dev/null && "$SAM" faidx "$OUT/$C.fa"
  # §6am: drop the previous run's GTFs FIRST — stale off_/on_ GTFs are what a silently failed exec
  # was diffed against, and they are what produced lost=0 gained=0 with nothing executed.
  rm -f "$OUT/off_$C.gtf" "$OUT/on_$C.gtf"
  RC=0
  RAYON_NUM_THREADS=4 "$RUSTLE" --vg --vg-layer2 --genome-fasta "$OUT/$C.fa" \
      -G "$OUT/st_$C.gtf" -L "$OUT/$C.bam" -o "$OUT/off_$C.gtf" 2>"$OUT/off_$C.log" || RC=$?
  # §6am: the return code was never checked; a failed arm must not be scored as a result.
  [[ $RC -eq 0 ]] || { fail_contig "$C" "gate-OFF rustle exited $RC (see $OUT/off_$C.log)"; continue; }
  RC=0
  RAYON_NUM_THREADS=4 RUSTLE_VG_CORE_GATE_TRACE=1 RUSTLE_VG_FAMILY_MIN_CORE_COVERAGE=0.13 \
      "$RUSTLE" --vg --vg-layer2 --genome-fasta "$OUT/$C.fa" \
      -G "$OUT/st_$C.gtf" -L "$OUT/$C.bam" -o "$OUT/on_$C.gtf" 2>"$OUT/on_$C.log" || RC=$?
  # §6am: same for the ON arm — an unrun ON arm is what makes "the gate loses nothing" unfalsifiable.
  [[ $RC -eq 0 ]] || { fail_contig "$C" "gate-ON rustle exited $RC (see $OUT/on_$C.log)"; continue; }
  # §6am: rustle can exit 0 without writing -o; a diff of files this run did not produce is stale evidence.
  [[ -f "$OUT/off_$C.gtf" && -f "$OUT/on_$C.gtf" ]] || { fail_contig "$C" "rustle exited 0 but did not write both GTFs"; continue; }
  # `grep -c` prints 0 and exits 1 on no match; the old inlined `|| echo 0` appended a SECOND line and
  # corrupted the summary row. Fall back on the assignment instead (same value when it does match).
  GATED=$(grep -c 'would_gate=true' "$OUT/on_$C.log" 2>/dev/null) || GATED=0
  RC=0
  ROW=$("$PY" "$DIFF" "$OUT/off_$C.gtf" "$OUT/on_$C.gtf" 2>"$OUT/onoff_$C.diff") || RC=$?
  # §6am: a failed or empty diff would append a blank row that sums as lost=0 gained=0.
  [[ $RC -eq 0 && -n "$ROW" ]] || { fail_contig "$C" "keyed diff exited $RC with row='$ROW'"; continue; }
  echo -e "${ROW}\t${GATED}" >> "$SUM"
  echo "[$C] reads=$RN gated=$GATED  ->  $ROW  (wall=$(($(date +%s)-CS))s)"
  rm -f "$OUT/$C.bam" "$OUT/$C.bam.bai" "$OUT/$C.fa" "$OUT/$C.fa.fai"
  touch "$OUT/onoff_$C.done"
done
echo "=== GENOME-WIDE on/off DONE wall=$(($(date +%s)-T0))s ==="
# §6am: a partial run must not reach the aggregate verdict at all — totals over contigs that never
# ran are the false PASS this instrument was confirmed to emit. Failed contigs are unlatched, so
# simply re-running resumes them.
if [[ "$FAILED" -gt 0 ]]; then
  echo "ABORT: $FAILED contig(s) failed to run; the summary is INCOMPLETE and its lost/gained totals are NOT a verdict" >&2
  exit 3
fi
"$PY" - <<'PY'
import csv
import sys
rows=list(csv.DictReader(open("/tmp/gw/onoff_summary.tsv"),delimiter="\t"))
# §6am guard: no rows, or no transcripts in the OFF arm, means the diff had nothing that COULD be
# lost — "LOST=0 GAINED=0" would then be an empty-evidence pass rather than a measurement.
if not rows:
    sys.exit("ABORT: onoff_summary.tsv has NO rows - no contig was scored (all latched by stale "
             "onoff_*.done files? delete them to re-measure). lost=0/gained=0 here is vacuous.")
S=lambda k:sum(int(r[k]) for r in rows)
if S('tx_off')==0:
    sys.exit("ABORT: 0 transcripts in the gate-OFF arm across all %d contig(s) - nothing could be "
             "lost, so this is not a verdict." % len(rows))
print(f"contigs={len(rows)} tx_off={S('tx_off')} tx_on={S('tx_on')} "
      f"LOST={S('lost')} GAINED={S('gained')} real_attr_changed={S('real_attr')} "
      f"cosmetic_renumber={S('cosmetic')} gated_merges={S('gated')}")
PY
echo "ONOFF_GW_DONE"
