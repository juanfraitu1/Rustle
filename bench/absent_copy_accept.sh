#!/usr/bin/env bash
# Acceptance runner for the --absent-copies flag.
#
# Guard A: flag-OFF byte-identical — runs the binary on a SMALL representative region
#           WITHOUT --absent-copies and diffs against the frozen o2_definitive.assignments.tsv
#           for that region; asserts diff=0, exactly 8 columns, NO *.dna_needs.tsv produced.
#
# Guard B: ON==OFF when nothing is admitted — runs the SAME region WITH --absent-copies;
#           asserts assignments diff=0 vs the OFF run AND that <out>.dna_needs.tsv exists
#           (even when empty, the flag produces the file).
#
# Usage:
#   bash bench/absent_copy_accept.sh
#
# Requires:
#   /home/juanfra/winloci_scratch/GGO_mm.bam (+ .bai)
#   /home/juanfra/winloci_scratch/GGO.fasta  (+ .fai)
#   /home/juanfra/winloci_scratch/o2_definitive.assignments.tsv  (frozen reference)
#   target/release/copy_assign

set -euo pipefail

BINARY="/mnt/c/Users/jfris/Desktop/Rustle/target/release/copy_assign"
BAM="/home/juanfra/winloci_scratch/GGO_mm.bam"
FASTA="/home/juanfra/winloci_scratch/GGO.fasta"
O2_DEF="/home/juanfra/winloci_scratch/o2_definitive.assignments.tsv"

# Small region: NC_073224.2:101576782-101609689 (GWFAM3/CAFAM0, ~749 reads, finishes in ~2s)
REGION="NC_073224.2:101576782-101609689"
FAMILY="CAFAM0"

export PATH="/home/juanfra/miniforge3/bin:$PATH"

SCRATCH=$(mktemp -d /tmp/absent_accept_XXXXXX)
trap 'rm -rf "$SCRATCH"' EXIT

PASS=0
FAIL=0

pass()  { echo "  PASS: $*"; PASS=$((PASS + 1)); }
fail()  { echo "  FAIL: $*"; FAIL=$((FAIL + 1)); }

echo "=== absent-copy acceptance (Guard A + Guard B) ==="
echo "Region: $REGION  |  scratch: $SCRATCH"
echo ""

# -----------------------------------------------------------------------
# Guard A: flag-OFF byte-identical
# -----------------------------------------------------------------------
echo "--- Guard A: flag-OFF byte-identical ---"

echo "$REGION" > "$SCRATCH/region.txt"

"$BINARY" \
  --bam "$BAM" \
  --fasta "$FASTA" \
  --regions "$SCRATCH/region.txt" \
  --min-copies 2 \
  --skip-poa-diagnostic \
  --out "$SCRATCH/off" 2>/dev/null

# A1: no dna_needs.tsv produced by the OFF path
if [ -f "$SCRATCH/off.dna_needs.tsv" ]; then
  fail "A1  OFF path produced off.dna_needs.tsv — it should NOT exist"
else
  pass "A1  OFF path: no dna_needs.tsv (correct)"
fi

# A2: exactly 8 columns in assignments header
N_COLS=$(head -1 "$SCRATCH/off.assignments.tsv" | tr '\t' '\n' | wc -l)
if [ "$N_COLS" -eq 8 ]; then
  pass "A2  assignments.tsv has exactly 8 columns"
else
  fail "A2  expected 8 columns, got $N_COLS"
fi

# A3: diff vs frozen o2_definitive (filter to the small region's family)
grep "$FAMILY" "$O2_DEF" | sort > "$SCRATCH/def_fam.tsv"
grep "$FAMILY" "$SCRATCH/off.assignments.tsv" | sort > "$SCRATCH/off_fam.tsv"
DIFF_LINES=$(diff "$SCRATCH/def_fam.tsv" "$SCRATCH/off_fam.tsv" | wc -l)
if [ "$DIFF_LINES" -eq 0 ]; then
  pass "A3  OFF assignments byte-identical to o2_definitive for $FAMILY (diff=0)"
else
  fail "A3  diff=$DIFF_LINES lines vs o2_definitive for $FAMILY"
  diff "$SCRATCH/def_fam.tsv" "$SCRATCH/off_fam.tsv" | head -10
fi

echo ""

# -----------------------------------------------------------------------
# Guard B: ON==OFF when no copy is admitted
# -----------------------------------------------------------------------
echo "--- Guard B: ON==OFF when nothing admitted ---"

"$BINARY" \
  --bam "$BAM" \
  --fasta "$FASTA" \
  --regions "$SCRATCH/region.txt" \
  --min-copies 2 \
  --skip-poa-diagnostic \
  --absent-copies \
  --out "$SCRATCH/on" 2>/dev/null

# B1: dna_needs.tsv EXISTS when --absent-copies is set (even if empty)
if [ -f "$SCRATCH/on.dna_needs.tsv" ]; then
  pass "B1  ON path: dna_needs.tsv exists"
else
  fail "B1  ON path: dna_needs.tsv MISSING — flag did not produce it"
fi

# B2: assignments diff=0 vs OFF run (ON==OFF when nothing admitted)
ON_DIFF=$(diff <(sort "$SCRATCH/off.assignments.tsv") <(sort "$SCRATCH/on.assignments.tsv") | wc -l)
if [ "$ON_DIFF" -eq 0 ]; then
  pass "B2  ON assignments == OFF assignments (diff=0, freeze guarantee holds)"
else
  fail "B2  ON vs OFF diff=$ON_DIFF lines — assignments changed unexpectedly"
  diff <(sort "$SCRATCH/off.assignments.tsv") <(sort "$SCRATCH/on.assignments.tsv") | head -10
fi

# B3: no AC_* copies admitted (this region has no collapsed copies)
if grep -q "^AC_" "$SCRATCH/on.quant.tsv" 2>/dev/null; then
  fail "B3  unexpected AC_* copy admitted for this reference region"
else
  pass "B3  no AC_* copies admitted (expected for this reference region)"
fi

echo ""

# -----------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------
TOTAL=$((PASS + FAIL))
echo "=== SUMMARY: $PASS/$TOTAL guards PASS, $FAIL FAIL ==="
if [ "$FAIL" -eq 0 ]; then
  echo "ACCEPTANCE: PASS"
  exit 0
else
  echo "ACCEPTANCE: FAIL"
  exit 1
fi
