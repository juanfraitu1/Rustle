#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# Tandem-VG benchmark harness.
#
# Scores whether intra-bundle tandem decomposition (RUSTLE_VG_TANDEM=1) makes the
# per-copy resolution of collapsed tandem arrays more ACCURATE — not merely more
# granular — by comparing `--vg` baseline vs `--vg`+tandem against a reference
# annotation with gffcompare.
#
# WHY THIS EXISTS: on a real multi-copy bam, tandem mode turns one collapsed
# mega-bundle (e.g. RBMY1 → 54 jumbled transcripts tagged copy_id 3 / cc 1.000)
# into per-copy sub-bundles (→ ~29 transcripts with real per-copy
# capacity_confidence). That is clearly more INFORMATIVE; this harness measures
# whether it is more CORRECT (Sn/Pr vs annotation) — the remaining gate before
# promoting RUSTLE_VG_TANDEM to default-on. (See
# docs/superpowers/specs/2026-06-02-intrabundle-tandem-vg-design.md.)
#
# Usage:
#   bench/tandem_benchmark.sh <reads.bam> <genome.fa> [annotation.gtf] [region.bed] [out_dir]
#
#   <reads.bam>        long-read alignments (a real multi-copy region/chromosome bam)
#   <genome.fa>        reference FASTA covering the bam's contig(s)
#   [annotation.gtf]   OPTIONAL ground-truth GTF. If given, gffcompare Sn/Pr is
#                      reported for both modes. If omitted (or absent), the harness
#                      still runs both modes and reports the STRUCTURAL comparison
#                      (transcript counts, per-array resolution, per-copy cc) so it
#                      is useful before annotation is available.
#   [region.bed]       OPTIONAL BED of the tandem-array regions. When given, Sn/Pr
#                      is ALSO reported restricted to these regions (where tandem
#                      acts), so genome-wide noise doesn't dilute the signal.
#   [out_dir]          OPTIONAL output dir (default: /tmp/tandem_bench).
#
# Env:
#   RUSTLE        path to the rustle binary (default: target/release/rustle)
#   GFFCOMPARE    path to gffcompare    (default: gffcompare, or ~/miniforge3/bin/gffcompare)
#
# Note: `--vg` assembly is non-deterministic in transcript ORDER (parallel), so
# compare COUNTS and Sn/Pr (order-independent), never raw line diffs.
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

BAM="${1:?usage: tandem_benchmark.sh <reads.bam> <genome.fa> [annotation.gtf] [region.bed] [out_dir]}"
GENOME="${2:?need <genome.fa>}"
ANNOT="${3:-}"
REGION_BED="${4:-}"
OUTDIR="${5:-/tmp/tandem_bench}"

RUSTLE="${RUSTLE:-$REPO_ROOT/target/release/rustle}"
GFFCOMPARE="${GFFCOMPARE:-gffcompare}"
command -v "$GFFCOMPARE" >/dev/null 2>&1 || GFFCOMPARE="$HOME/miniforge3/bin/gffcompare"

[[ -x "$RUSTLE" ]] || { echo "ERROR: rustle not found/executable at $RUSTLE (set RUSTLE=...)" >&2; exit 1; }
[[ -f "$BAM" ]]    || { echo "ERROR: bam not found: $BAM" >&2; exit 1; }
[[ -f "$GENOME" ]] || { echo "ERROR: genome not found: $GENOME" >&2; exit 1; }
mkdir -p "$OUTDIR"

BASE_GTF="$OUTDIR/baseline_vg.gtf"
TAND_GTF="$OUTDIR/tandem_vg.gtf"
BASE_LOG="$OUTDIR/baseline_vg.log"
TAND_LOG="$OUTDIR/tandem_vg.log"

echo "═══════════════════════════════════════════════════════════════════════"
echo " Tandem-VG benchmark"
echo "   bam:        $BAM"
echo "   genome:     $GENOME"
echo "   annotation: ${ANNOT:-<none — structural comparison only>}"
echo "   regions:    ${REGION_BED:-<whole input>}"
echo "   out:        $OUTDIR"
echo "═══════════════════════════════════════════════════════════════════════"

# ── 1. Run both modes ────────────────────────────────────────────────────────
echo "[1/3] running --vg baseline (tandem OFF) ..."
RUSTLE_VG_ANCHOR_PRIOR=1 "$RUSTLE" -L "$BAM" --vg --genome-fasta "$GENOME" \
    -o "$BASE_GTF" 2>"$BASE_LOG"

echo "[1/3] running --vg + RUSTLE_VG_TANDEM=1 (tandem ON) ..."
RUSTLE_VG_TANDEM=1 RUSTLE_VG_ANCHOR_PRIOR=1 "$RUSTLE" -L "$BAM" --vg --genome-fasta "$GENOME" \
    -o "$TAND_GTF" 2>"$TAND_LOG"

tx_count() { grep -c $'\ttranscript\t' "$1" 2>/dev/null || echo 0; }
BASE_N=$(tx_count "$BASE_GTF"); TAND_N=$(tx_count "$TAND_GTF")

# ── 2. Structural / tandem-specific metrics ──────────────────────────────────
echo
echo "[2/3] structural comparison"
echo "   transcripts:        baseline=$BASE_N   tandem=$TAND_N   (delta $((TAND_N - BASE_N)))"
echo "   arrays decomposed:  $(grep -c 'copy sub-bundles' "$TAND_LOG" 2>/dev/null || echo 0)"
grep 'copy sub-bundles' "$TAND_LOG" 2>/dev/null | sed 's/^/      /' || true
echo "   weight-floor:       $(grep 'weight floor' "$TAND_LOG" 2>/dev/null | sed 's/.*floor: //' || echo 'not fired')"
TC=$(grep -c 'tandem_copy "true"' "$TAND_GTF" 2>/dev/null || echo 0)
TC_LOW=$(awk '/tandem_copy "true"/{m=match($0,/capacity_confidence "([0-9.]+)"/,a); if(a[1]+0<0.5)n++} END{print n+0}' "$TAND_GTF" 2>/dev/null || echo 0)
echo "   tandem-copy tx:     $TC   (of which cc<0.5, honest low-confidence: $TC_LOW)"

# Per-array transcript delta (if a region BED was supplied).
if [[ -n "$REGION_BED" && -f "$REGION_BED" ]]; then
  echo "   per-array transcript counts (baseline → tandem):"
  while read -r chrom start end label _rest; do
    [[ -z "${chrom:-}" || "$chrom" == \#* ]] && continue
    b=$(awk -v c="$chrom" -v s="$start" -v e="$end" '$1==c && $3=="transcript" && $4>=s && $4<=e' "$BASE_GTF" | wc -l)
    t=$(awk -v c="$chrom" -v s="$start" -v e="$end" '$1==c && $3=="transcript" && $4>=s && $4<=e' "$TAND_GTF" | wc -l)
    printf "      %-20s %s:%s-%s   %d → %d\n" "${label:-region}" "$chrom" "$start" "$end" "$b" "$t"
  done < "$REGION_BED"
fi

# ── 3. Accuracy vs annotation (gffcompare) ───────────────────────────────────
echo
if [[ -z "$ANNOT" || ! -f "$ANNOT" ]]; then
  echo "[3/3] accuracy scoring SKIPPED — no annotation supplied."
  echo "      Re-run with a ground-truth GTF as arg 3 to get Sn/Pr:"
  echo "        bench/tandem_benchmark.sh \"$BAM\" \"$GENOME\" <annotation.gtf> [region.bed]"
  echo "      A per-copy-attribution score additionally needs per-read truth"
  echo "      (e.g. a badread/gen_reads synthetic where each read's source copy is known)."
  exit 0
fi

command -v "$GFFCOMPARE" >/dev/null 2>&1 || { echo "ERROR: gffcompare not found (set GFFCOMPARE=...)" >&2; exit 1; }

# Optionally restrict the reference to the array regions.
REF="$ANNOT"
if [[ -n "$REGION_BED" && -f "$REGION_BED" ]]; then
  if command -v bedtools >/dev/null 2>&1; then
    REF="$OUTDIR/annot_regions.gtf"
    bedtools intersect -u -a "$ANNOT" -b "$REGION_BED" > "$REF" 2>/dev/null || REF="$ANNOT"
    echo "[3/3] accuracy vs annotation (restricted to ${REGION_BED##*/})"
  else
    echo "[3/3] accuracy vs annotation (bedtools absent → whole annotation)"
  fi
else
  echo "[3/3] accuracy vs annotation (whole input)"
fi

score() { # <query.gtf> <prefix> → echoes "Sn Pr" from the transcript-level line
  local q="$1" pfx="$2"
  "$GFFCOMPARE" -R -Q -r "$REF" -o "$OUTDIR/$pfx" "$q" >/dev/null 2>&1 || true
  # gffcompare .stats: "Transcript level:    <Sn>     |    <Pr>    |"
  awk -F'[:|]' '/Transcript level/{gsub(/[^0-9.]/,"",$2); gsub(/[^0-9.]/,"",$3); print $2, $3}' \
      "$OUTDIR/$pfx.stats" 2>/dev/null | head -1
}

read -r BASE_SN BASE_PR <<<"$(score "$BASE_GTF" gffcmp_baseline)"
read -r TAND_SN TAND_PR <<<"$(score "$TAND_GTF" gffcmp_tandem)"
: "${BASE_SN:=NA}"; : "${BASE_PR:=NA}"; : "${TAND_SN:=NA}"; : "${TAND_PR:=NA}"

echo
echo "   ┌──────────────┬──────────┬──────────┬──────────────┐"
echo "   │ mode         │   Sn (%) │   Pr (%) │ transcripts  │"
echo "   ├──────────────┼──────────┼──────────┼──────────────┤"
printf  "   │ %-12s │ %8s │ %8s │ %12s │\n" "--vg base" "$BASE_SN" "$BASE_PR" "$BASE_N"
printf  "   │ %-12s │ %8s │ %8s │ %12s │\n" "--vg+tandem" "$TAND_SN" "$TAND_PR" "$TAND_N"
echo "   └──────────────┴──────────┴──────────┴──────────────┘"
echo
echo "   VERDICT: tandem is promotion-worthy if Sn is up (recovers true per-copy"
echo "   isoforms) AND Pr is flat-or-up (the extra granularity is correct, not"
echo "   spurious). Sn↑/Pr↓ ⇒ over-enumeration; investigate before promoting."
echo "   gffcompare outputs: $OUTDIR/gffcmp_{baseline,tandem}.stats"
