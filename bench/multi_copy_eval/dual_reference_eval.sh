#!/usr/bin/env bash
# DUAL-REFERENCE evaluation: "does VG find >= baseline AND more?"
#
# Compares VG output against TWO references at once:
#   (1) the BASELINE floor = real StringTie's GTF (rustle baseline ≡ StringTie, so this is "what
#       the canonical tool finds"); answers "does VG recover >= baseline?"
#   (2) the ANNOTATION truth = a RefSeq GFF; answers "does VG find MORE of the truth than baseline?"
#
# Why NOT -G: guiding assembly with a GTF then scoring against it is CIRCULAR (rustle's own CLI warns
# "-G ... not for headline sensitivity/precision vs truth"). Using StringTie's de-novo output AS a
# reference is the honest equivalent of the "ensure ~100% baseline" goal — and rustle baseline already
# reproduces StringTie at 100/100.
#
# The decisive step: a transcript VG "loses" vs StringTie is only a Tier-1 violation if it is a REAL
# (RefSeq-matching) transcript. A lost StringTie FP (class != '=' vs RefSeq) is fine — good, even.
#
# Usage: dual_reference_eval.sh <reads.bam> <refseq.gff> [genome.fa] [vg_extra_env...]
#   e.g. dual_reference_eval.sh /tmp/tspy.bam /tmp/rbmy_truth.gff ../GGO.fasta \
#          RUSTLE_VG_TANDEM=1 RUSTLE_VG_TANDEM_PRIMARY_JUNCTIONS=1 RUSTLE_VG_ANCHOR_PRIOR=0 RUSTLE_VG_JOINT_STRAND_EM=0
set -euo pipefail
BAM="$1"; REF="$2"; GENOME="${3:-}"; shift; shift; [ $# -gt 0 ] && shift || true
VG_ENV=("$@")
HERE="$(cd "$(dirname "$0")" && pwd)"; ROOT="$(cd "$HERE/../.." && pwd)"
RUSTLE="$ROOT/target/release/rustle"; ST="$ROOT/tools/stringtie/stringtie"
W="${OUT:-/tmp/dualref}"; mkdir -p "$W"

tx() { gffcompare -r "$1" "$2" -o "$W/c" 2>/dev/null; grep "Transcript level" "$W/c.stats"|grep -oE '[0-9]+\.[0-9]+'|paste -sd'/'; }

"$ST" -L "$BAM" -o "$W/stringtie.gtf" 2>/dev/null
"$RUSTLE" -L "$BAM" -o "$W/baseline.gtf" 2>/dev/null
gfa=(); [ -n "$GENOME" ] && gfa=(--genome-fasta "$GENOME")
env "${VG_ENV[@]}" "$RUSTLE" --vg --vg-snp "${gfa[@]}" -L "$BAM" -o "$W/vg.gtf" 2>/dev/null

echo "## Dual-reference evaluation ($(basename "$BAM"))"
printf "  %-18s | vs RefSeq(truth) | vs StringTie(floor)\n" query
printf "  %-18s | %-15s | %s\n" "StringTie(=base)" "$(tx "$REF" "$W/stringtie.gtf")" "(self)"
printf "  %-18s | %-15s | %s\n" "rustle baseline"  "$(tx "$REF" "$W/baseline.gtf")" "$(tx "$W/stringtie.gtf" "$W/baseline.gtf")"
printf "  %-18s | %-15s | %s\n" "VG"               "$(tx "$REF" "$W/vg.gtf")"       "$(tx "$W/stringtie.gtf" "$W/vg.gtf")"

echo "  -- VG GAINS (RefSeq tx VG matches that StringTie does NOT = genuine discovery) --"
gffcompare -r "$REF" "$W/vg.gtf" -o "$W/vr" 2>/dev/null
gffcompare -r "$REF" "$W/stringtie.gtf" -o "$W/sr" 2>/dev/null
comm -23 <(awk -F'\t' 'NR>1&&$3=="="{print $2}' "$W/vr."*.tmap|sort -u) \
         <(awk -F'\t' 'NR>1&&$3=="="{print $2}' "$W/sr."*.tmap|sort -u) | sed 's/^/    + /'

echo "  -- VG LOSSES (StringTie tx VG drops) classified vs RefSeq (REAL=Tier-1 violation, FP=fine) --"
gffcompare -r "$W/vg.gtf" "$W/stringtie.gtf" -o "$W/sv" 2>/dev/null
for t in $(awk -F'\t' 'NR>1&&$3!="="{print $5}' "$W/sv.stringtie.gtf.tmap"); do
  cls=$(awk -F'\t' -v t="$t" 'NR>1&&$5==t{print $3}' "$W/sr.stringtie.gtf.tmap")
  echo "    - $t: vs-RefSeq=$cls  $([ "$cls" = "=" ] && echo 'REAL (Tier-1 LOSS)' || echo 'StringTie FP (ok)')"
done
echo "  VERDICT: VG >= baseline iff every LOSS is a StringTie FP; MORE iff any GAIN is listed."
