#!/usr/bin/env bash
# Regional two-locus parity: StringTie vs Rustle vs gffcompare -R -Q.
# Requires: rustle (cargo build --release), stringtie, gffcompare, samtools, grep on out_full.gtf
#
# Usage (from repo root):
#   scripts/run_two_locus_parity_bench.sh
#
# Optional env:
#   STRINGTIE GFFCOMPARE  Override tool paths (default: from PATH)

set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
P=bench/parity_two_locus
mkdir -p "$P"

ST="${STRINGTIE:-$(command -v stringtie)}"
GC="${GFFCOMPARE:-$(command -v gffcompare)}"
RUST="${RUST:-$ROOT/target/release/rustle}"

if [[ ! -x "$RUST" ]]; then
  echo "Build rustle first: cargo build --release" >&2
  exit 1
fi
if [[ -z "$ST" || -z "$GC" ]]; then
  echo "Need stringtie and gffcompare on PATH or STRINGTIE / GFFCOMPARE" >&2
  exit 1
fi

REF_FULL="${REF_FULL:-$ROOT/out_full.gtf}"
if [[ ! -f "$REF_FULL" ]]; then
  echo "Need reference slice source at REF_FULL=$REF_FULL" >&2
  exit 1
fi

grep 'gene_id "STRG.171"' "$REF_FULL" >"$P/ref_STRG171.gtf"
grep 'gene_id "STRG.212"' "$REF_FULL" >"$P/ref_STRG212.gtf"

B212="$P/GGO_19_STRG212_reads.bam"
if [[ ! -f "$B212" ]]; then
  samtools view -b "$ROOT/GGO_19.bam" NC_073243.2:30720000-30850000 -o "$B212"
  samtools index "$B212"
fi

run_cmp() {
  local tag=$1 bam=$2 ref=$3 out=$4
  "$GC" -r "$ref" -R -Q -o "$P/cmp_${tag}" "$out" >"$P/cmp_${tag}.summary.txt" 2>&1 || true
}

echo "=== StringTie STRG171 ==="
"$ST" -L -p 4 -o "$P/stringtie_STRG171.gtf" "$ROOT/GGO_19_STRG171_reads.bam"
run_cmp st_STRG171 "$ROOT/GGO_19_STRG171_reads.bam" "$P/ref_STRG171.gtf" "$P/stringtie_STRG171.gtf"

echo "=== Rustle default STRG171 ==="
"$RUST" -o "$P/rustle_default_STRG171.gtf" "$ROOT/GGO_19_STRG171_reads.bam" -L -p 4
run_cmp rd_STRG171 "$ROOT/GGO_19_STRG171_reads.bam" "$P/ref_STRG171.gtf" "$P/rustle_default_STRG171.gtf"

echo "=== Rustle + ref-junction-witness STRG171 ==="
"$RUST" -o "$P/rustle_parity_STRG171.gtf" "$ROOT/GGO_19_STRG171_reads.bam" -L -p 4 --ref-junction-witness
run_cmp rp_STRG171 "$ROOT/GGO_19_STRG171_reads.bam" "$P/ref_STRG171.gtf" "$P/rustle_parity_STRG171.gtf"

echo "=== StringTie STRG212 ==="
"$ST" -L -p 4 -o "$P/stringtie_STRG212.gtf" "$B212"
run_cmp st_STRG212 "$B212" "$P/ref_STRG212.gtf" "$P/stringtie_STRG212.gtf"

echo "=== Rustle default STRG212 ==="
"$RUST" -o "$P/rustle_default_STRG212.gtf" "$B212" -L -p 4
run_cmp rd_STRG212 "$B212" "$P/ref_STRG212.gtf" "$P/rustle_default_STRG212.gtf"

echo "=== Rustle + ref-junction-witness STRG212 ==="
"$RUST" -o "$P/rustle_parity_STRG212.gtf" "$B212" -L -p 4 --ref-junction-witness
run_cmp rp_STRG212 "$B212" "$P/ref_STRG212.gtf" "$P/rustle_parity_STRG212.gtf"

echo "Done. Key lines:"
grep -h -E 'Query mRNAs|Reference mRNAs|Transcript level|Matching transcripts' "$P"/cmp_{st,rd,rp}_STRG*.stats 2>/dev/null | paste - - - - | head -20 || true
echo "See $P/RESULTS.txt for interpretation."
