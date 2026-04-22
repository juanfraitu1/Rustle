#!/bin/bash
# Compare IsoSeq Collapse output with reference annotation
# Usage: ./compare_with_reference.sh <our_output.gtf/gff> <reference.gtf/gff> [output_prefix]

set -e

if [ $# -lt 2 ]; then
    echo "Usage: $0 <our_output.gtf/gff> <reference.gtf/gff> [output_prefix]"
    echo ""
    echo "Example:"
    echo "  $0 output.gtf GGO_genomic.gtf my_comparison"
    echo ""
    echo "This will:"
    echo "  1. Run gffcompare to classify isoforms"
    echo "  2. Generate statistics and summary"
    echo "  3. Create high-confidence ground truth based on reference matches"
    exit 1
fi

OUR_OUTPUT="$1"
REF_ANNOTATION="$2"
PREFIX="${3:-comparison}"

# Detect file format
OUR_FORMAT="gtf"
if [[ "$OUR_OUTPUT" == *.gff ]] || [[ "$OUR_OUTPUT" == *.gff3 ]]; then
    OUR_FORMAT="gff"
fi

REF_FORMAT="gtf"
if [[ "$REF_ANNOTATION" == *.gff ]] || [[ "$REF_ANNOTATION" == *.gff3 ]]; then
    REF_FORMAT="gff"
fi

echo "========================================"
echo "  IsoSeq Collapse Comparison Tool"
echo "========================================"
echo ""
echo "Our output:    $OUR_OUTPUT ($OUR_FORMAT)"
echo "Reference:     $REF_ANNOTATION ($REF_FORMAT)"
echo "Output prefix: $PREFIX"
echo ""

# Step 1: Clean reference (keep only transcript and exon features)
echo "Step 1: Preparing reference annotation..."
REF_CLEAN="${PREFIX}_ref_clean.gtf"

# Extract transcript and exon lines, handling both GTF and GFF
if [ "$REF_FORMAT" == "gff" ]; then
    # GFF format - convert IDs
    awk -F'\t' '$3=="transcript" || $3=="mRNA" || $3=="exon"' "$REF_ANNOTATION" | \
        awk -F'\t' '{
            gsub(/ID=/, "transcript_id \"", $9)
            gsub(/;Parent=/, "\";gene_id \"", $9)
            gsub(/;$/, "\"", $9)
            print
        }' > "$REF_CLEAN" 2>/dev/null || \
    awk -F'\t' '$3=="transcript" || $3=="mRNA" || $3=="exon"' "$REF_ANNOTATION" > "$REF_CLEAN"
else
    # GTF format
    awk -F'\t' '$3=="transcript" || $3=="exon"' "$REF_ANNOTATION" > "$REF_CLEAN"
fi

REF_COUNT=$(grep -c 'transcript_id' "$REF_CLEAN" 2>/dev/null || echo 0)
echo "  Reference transcripts: $REF_COUNT"

# Step 2: Run gffcompare
echo ""
echo "Step 2: Running gffcompare..."
gffcompare -r "$REF_CLEAN" -o "${PREFIX}" "$OUR_OUTPUT" 2>&1 | grep -E "(reference|query|matching)" || true

# Step 3: Parse and display results
echo ""
echo "Step 3: Classification Results"
echo "================================"

if [ -f "${PREFIX}.stats" ]; then
    cat "${PREFIX}.stats"
else
    echo "Warning: No stats file generated"
fi

# Count class codes
echo ""
echo "Detailed Classification:"
echo "------------------------"

if [ -f "${PREFIX}.combined.gtf" ]; then
    echo "  = : Exact match to reference        $(grep -c 'class_code "="' "${PREFIX}.combined.gtf" 2>/dev/null || echo 0)"
    echo "  c : Contained in reference          $(grep -c 'class_code "c"' "${PREFIX}.combined.gtf" 2>/dev/null || echo 0)"
    echo "  j : Novel, some reference junctions $(grep -c 'class_code "j"' "${PREFIX}.combined.gtf" 2>/dev/null || echo 0)"
    echo "  e : Novel, reference exons          $(grep -c 'class_code "e"' "${PREFIX}.combined.gtf" 2>/dev/null || echo 0)"
    echo "  k : Novel isoform                   $(grep -c 'class_code "k"' "${PREFIX}.combined.gtf" 2>/dev/null || echo 0)"
    echo "  u : Intergenic                      $(grep -c 'class_code "u"' "${PREFIX}.combined.gtf" 2>/dev/null || echo 0)"
    echo "  x : Exonic antisense                $(grep -c 'class_code "x"' "${PREFIX}.combined.gtf" 2>/dev/null || echo 0)"
    echo "  s : Intronic antisense              $(grep -c 'class_code "s"' "${PREFIX}.combined.gtf" 2>/dev/null || echo 0)"
    echo "  o : Overlap with reference          $(grep -c 'class_code "o"' "${PREFIX}.combined.gtf" 2>/dev/null || echo 0)"
    echo "  p : Polycistronic                   $(grep -c 'class_code "p"' "${PREFIX}.combined.gtf" 2>/dev/null || echo 0)"
fi

# Step 4: Create ground truth
echo ""
echo "Step 4: Creating Ground Truth Sets"
echo "===================================="

# Exact matches (FSM - Full Splice Match)
if [ -f "${PREFIX}.combined.gtf" ]; then
    grep 'class_code "="' "${PREFIX}.combined.gtf" 2>/dev/null | \
        awk -F'\t' '$3=="transcript"' > "${PREFIX}_exact_matches.gtf" || \
        touch "${PREFIX}_exact_matches.gtf"
    
    EXACT_COUNT=$(wc -l < "${PREFIX}_exact_matches.gtf")
    echo "  Exact matches (FSM):     $EXACT_COUNT"
    
    # Contained matches (ISM - Incomplete Splice Match)
    grep -E 'class_code "[cj]"' "${PREFIX}.combined.gtf" 2>/dev/null | \
        awk -F'\t' '$3=="transcript"' > "${PREFIX}_contained.gtf" || \
        touch "${PREFIX}_contained.gtf"
    
    CONTAINED_COUNT=$(wc -l < "${PREFIX}_contained.gtf")
    echo "  Contained (ISM):         $CONTAINED_COUNT"
    
    # High confidence (= and c)
    cat "${PREFIX}_exact_matches.gtf" "${PREFIX}_contained.gtf" 2>/dev/null | \
        sort -u > "${PREFIX}_high_confidence.gtf" || \
        touch "${PREFIX}_high_confidence.gtf"
    
    HC_COUNT=$(wc -l < "${PREFIX}_high_confidence.gtf")
    echo "  High confidence total:   $HC_COUNT"
fi

# Step 5: Summary
echo ""
echo "========================================"
echo "  Summary"
echo "========================================"

TOTAL=$(grep -c 'transcript_id' "$OUR_OUTPUT" 2>/dev/null | head -1 || echo 0)
if [ -f "${PREFIX}_exact_matches.gtf" ]; then
    EXACT=$(wc -l < "${PREFIX}_exact_matches.gtf")
    PERCENTAGE=$(awk "BEGIN {printf \"%.1f\", ($EXACT/$TOTAL)*100}")
    echo "Total isoforms:          $TOTAL"
    echo "Exact reference matches: $EXACT ($PERCENTAGE%)"
fi

echo ""
echo "Output files:"
echo "  ${PREFIX}.stats              - Comparison statistics"
echo "  ${PREFIX}.combined.gtf       - Combined annotation with class codes"
echo "  ${PREFIX}_exact_matches.gtf  - Exact reference matches (FSM)"
echo "  ${PREFIX}_high_confidence.gtf - High-confidence isoforms"
echo "  ${PREFIX}_ref_clean.gtf      - Cleaned reference"
echo ""
echo "To view detailed stats: cat ${PREFIX}.stats"
echo "To view exact matches:  head ${PREFIX}_exact_matches.gtf"

# Cleanup
rm -f "$REF_CLEAN"

echo ""
echo "Done!"
