#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 4 ]]; then
  cat <<'EOF'
usage: scripts/run_rustle_gffcompare.sh <rustle_bin> <input_bam> <reference_gtf> <output_prefix> [rustle args...]

  <reference_gtf> is ONLY for gffcompare -r (truth for sensitivity/precision). Do NOT pass it to rustle
  unless you intentionally want guided assembly (-G); guided runs inflate metrics vs de novo truth tests.

De novo evaluation (recommended; comparable to StringTie without -G):
  scripts/run_rustle_gffcompare.sh \
    ./target/release/rustle \
    sample.bam \
    truth_reference.gtf \
    out_prefix \
    -L -p "$(nproc)"

Regional mini-BAM / single-locus GTF (locus-aligned gffcompare metrics):
  GFFCOMPARE_EXTRA_FLAGS='-R -Q' scripts/run_rustle_gffcompare.sh ...

Guided assembly (debug / parity only — not for headline benchmark numbers):
  ... same ...  -L -p 16 -G /path/to/same/or/other.gtf

StringTie de novo baseline (same BAM, for comparison):
  stringtie -L -p "$(nproc)" -o out_prefix.stringtie.gtf sample.bam
  gffcompare -r truth_reference.gtf -o out_prefix.stringtie.cmp out_prefix.stringtie.gtf
EOF
  exit 1
fi

rustle_bin=$1
input_bam=$2
reference_gtf=$3
output_prefix=$4
shift 4

query_gtf="${output_prefix}.gtf"
stats_prefix="${output_prefix}.gffcmp"
stats_file="${stats_prefix}.stats"
summary_file="${stats_prefix}.summary.txt"
report_file="${stats_prefix}"
json_file="${stats_prefix}.json"

"${rustle_bin}" \
  --help >/tmp/rustle_help_check.$$ 2>&1 || true

if rg -q -- '--input <INPUT>' /tmp/rustle_help_check.$$; then
  "${rustle_bin}" \
    --input "${input_bam}" \
    --output "${query_gtf}" \
    --reference "${reference_gtf}" \
    "$@"
else
  "${rustle_bin}" \
    -o "${query_gtf}" \
    "${input_bam}" \
    "$@"
fi

rm -f /tmp/rustle_help_check.$$

# Optional: export GFFCOMPARE_EXTRA_FLAGS='-R -Q' for regional/mini-BAM query GTFs (locus-aligned Sn/Precision).
# See BLOCKER_PATCH_GUIDE.md ("Same-locus / regional query GTFs").
# shellcheck disable=SC2086
gffcompare -r "${reference_gtf}" ${GFFCOMPARE_EXTRA_FLAGS:-} -o "${stats_prefix}" "${query_gtf}" >"${summary_file}" 2>&1

printf 'query_gtf\t%s\n' "${query_gtf}"
printf 'gffcompare_summary\t%s\n' "${summary_file}"
printf 'gffcompare_report\t%s\n' "${report_file}"

# gffcompare versions differ: some emit <prefix>.stats, others emit the summary into <prefix>.
if [[ -f "${stats_file}" ]]; then
  printf 'gffcompare_stats\t%s\n' "${stats_file}"
else
  printf 'gffcompare_stats\t%s\n' "${report_file}"
fi

# Best-effort JSON summary (stable fields for diffs/automation).
# gffcompare emits .tmap as: <prefix>.<query_gtf_filename>.tmap
tmap_file="${stats_prefix}.$(basename "${query_gtf}").tmap"
if [[ -f "${stats_file}" && -f "${tmap_file}" ]]; then
  python3 scripts/gffcompare_to_json.py \
    --reference-gtf "${reference_gtf}" \
    --stats "${stats_file}" \
    --tmap "${tmap_file}" \
    --out-json "${json_file}" >/dev/null 2>&1 || true
  if [[ -f "${json_file}" ]]; then
    printf 'gffcompare_json\t%s\n' "${json_file}"
  fi
fi

if [[ -f "${stats_file}" ]]; then
  awk '/^# Query mRNAs|^# Reference mRNAs|^\\s+Transcript level|Matching transcripts|Missed exons|Missed introns/ {print}' "${stats_file}"
elif [[ -f "${report_file}" ]]; then
  awk '/^# Query mRNAs|^# Reference mRNAs|^\\s+Transcript level|Matching transcripts|Missed exons|Missed introns|Sensitivity|Precision/ {print}' "${report_file}"
else
  awk '/^# Query mRNAs|^# Reference mRNAs|^\\s+Transcript level|Matching transcripts|Missed exons|Missed introns|Sensitivity|Precision/ {print}' "${summary_file}"
fi
