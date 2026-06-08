#!/usr/bin/env bash
# U1: generate per-arm transcriptomes from the full BAM (no windows).
# Canonical-sorts each GTF (rustle GTF line order has rayon nondeterminism).
source "$(dirname "$0")/config.sh"
set -euo pipefail
GENDIR="$OUTDIR/gtf"; mkdir -p "$GENDIR"

# Optional chrom subset for smoke testing (CHROM_SUBSET in config / env).
INBAM="$BAM"
if [ -n "$CHROM_SUBSET" ]; then
  INBAM="$OUTDIR/subset.bam"
  # CHROM_SUBSET may list several space-separated regions -> intentionally unquoted.
  # shellcheck disable=SC2086
  "$SAMTOOLS" view -b "$BAM" $CHROM_SUBSET > "$INBAM"
  "$SAMTOOLS" index "$INBAM"
fi

# Deterministic GTF order: by chrom, start, feature.
canon_sort() { sort -t$'\t' -k1,1 -k4,4n -k3,3 "$1" > "$1.sorted" && mv "$1.sorted" "$1"; }

PROV="$OUTDIR/provenance.json"
{
  echo "{"
  echo "  \"config_hash\": \"$(config_hash)\","
  echo "  \"rustle_commit\": \"$(cd "$ROOT" && git rev-parse HEAD)\","
  echo "  \"stringtie_version\": \"$($STRINGTIE --version 2>&1 | head -1)\","
  echo "  \"chrom_subset\": \"${CHROM_SUBSET:-WHOLE_GENOME}\","
  echo "  \"bam\": \"$INBAM\""
  echo "}"
} > "$PROV"

if [ "$ARM_STRINGTIE" = 1 ]; then
  # De novo (no -G): fair "find it from reads alone" comparison.
  "$STRINGTIE" -L -o "$GENDIR/stringtie.gtf" "$INBAM"
  canon_sort "$GENDIR/stringtie.gtf"
fi
if [ "$ARM_RUSTLE_VG" = 1 ]; then
  env $RUSTLE_VG_ENV "$RUSTLE" --vg --vg-snp --genome-fasta "$GENOME_FASTA" \
    -L "$INBAM" -o "$GENDIR/rustle_vg.gtf"
  canon_sort "$GENDIR/rustle_vg.gtf"
fi
if [ "$ARM_RUSTLE_VG_RAW" = 1 ]; then
  env $RUSTLE_VG_RAW_ENV "$RUSTLE" --vg --vg-snp --genome-fasta "$GENOME_FASTA" \
    -L "$INBAM" -o "$GENDIR/rustle_vg_raw.gtf"
  canon_sort "$GENDIR/rustle_vg_raw.gtf"
fi
if [ "$ARM_RUSTLE_PRIMARY" = 1 ]; then
  "$RUSTLE" -L "$INBAM" -o "$GENDIR/rustle_primary.gtf"
  canon_sort "$GENDIR/rustle_primary.gtf"
fi
echo "[U1] GTFs in $GENDIR"; ls -la "$GENDIR"
