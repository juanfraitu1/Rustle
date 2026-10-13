#!/bin/bash
# Fetch the Hallast et al. 2023 (HGSVC) per-man chrY assemblies and verify md5 (spec 2026-10-12-population-fixed-psvs §1).
#   fetch.sh OUTDIR
# Downloads run in parallel (xargs -P 6): EBI throttles per connection. Resumable: files whose md5 already matches are skipped.
set -euo pipefail
OUT=${1:?OUTDIR}
BASE=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/working/20230412_sigY_assembly_data
mkdir -p "$OUT"
curl -s --fail "$BASE/MANIFEST_20230412_sigY_assembly_data.tsv" -o "$OUT/MANIFEST.tsv"
awk -F'\t' '$1 ~ /verkko_1.0_release\/chrY\// && $1 !~ /\.fai$/ {print $1"\t"$3}' "$OUT/MANIFEST.tsv" > "$OUT/want.tsv"

fetch_one() {
  local path="$1" md5="$2" f got
  f=$(basename "$path")
  if [ ! -s "$OUT/$f" ] || [ "$(md5sum < "$OUT/$f" | cut -d' ' -f1)" != "$md5" ]; then
    curl -s --fail --retry 3 "$BASE/${path#20230412_sigY_assembly_data/}" -o "$OUT/$f.part" && mv "$OUT/$f.part" "$OUT/$f" || { echo "download failed: $f" >&2; return 1; }
  fi
  got=$(md5sum < "$OUT/$f" | cut -d' ' -f1)
  [ "$got" = "$md5" ] || { echo "md5 mismatch: $f" >&2; return 1; }
  samtools faidx "$OUT/$f"
}
export -f fetch_one
export OUT BASE

# want.tsv fields (path, md5) contain no spaces, so xargs -n 2 yields one pair per call.
cut -f1,2 "$OUT/want.tsv" | tr '\t' ' ' | xargs -P 6 -n 2 bash -c 'fetch_one "$0" "$1"'

n=$(wc -l < "$OUT/want.tsv")
{
  echo "# Hallast et al. 2023 Y assemblies (HGSVC, Verkko 1.0), fetched $(date -I)"
  echo "source: $BASE/assemblies/verkko_1.0_release/chrY/"
  echo "n_assemblies: $n (md5-verified against MANIFEST.tsv)"
  echo "samples:"; cut -f1 "$OUT/want.tsv" | xargs -n1 basename | cut -d. -f1 | sed 's/^/  /'
} > "$OUT/SOURCE.md"
echo "$n assemblies OK in $OUT"
