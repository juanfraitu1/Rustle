#!/bin/bash
# Create ground truth by comparing against reference annotation
# Usage: ./create_ground_truth.sh <our_output.gtf> <reference.gtf> <output_prefix>

set -e

OUR_GTF="$1"
REF_GTF="$2"
PREFIX="$3"

if [ -z "$OUR_GTF" ] || [ -z "$REF_GTF" ] || [ -z "$PREFIX" ]; then
    echo "Usage: $0 <our_output.gtf> <reference.gtf> <output_prefix>"
    echo "Example: $0 output.gtf GGO_genomic.gtf ground_truth"
    exit 1
fi

echo "=== Creating Ground Truth ==="
echo "Input: $OUR_GTF"
echo "Reference: $REF_GTF"
echo "Output prefix: $PREFIX"
echo ""

# Step 1: Run gffcompare to classify our isoforms
echo "Step 1: Running gffcompare..."
gffcompare -r "$REF_GTF" -o "${PREFIX}_cmp" "$OUR_GTF"

# Step 2: Show classification summary
echo ""
echo "Step 2: Classification summary:"
echo ""
echo "Class codes:"
echo "  = : Exact match to reference"
echo "  c : Contained in reference (subset)"
echo "  j : Novel isoform with some reference junctions"
echo "  e : Novel isoform with reference exons but new junctions"
echo "  k : Novel isoform"
echo ""

# Count each class
for code in "=" c j e k; do
    count=$(grep -c "class_code \"$code\"" "${PREFIX}_cmp.combined.gtf" 2>/dev/null || echo 0)
    echo "  $code : $count"
done

# Step 3: Create high-confidence ground truth
echo ""
echo "Step 3: Creating high-confidence ground truth..."

# Extract FSM (Full Splice Match) and ISM (Incomplete Splice Match)
# These are the most reliable isoforms
grep 'class_code "="' "${PREFIX}_cmp.combined.gtf" | \
    awk -F'\t' '$3=="transcript"' > "${PREFIX}_exact_matches.gtf" 2>/dev/null || touch "${PREFIX}_exact_matches.gtf"

grep -E 'class_code "[cj]"' "${PREFIX}_cmp.combined.gtf" | \
    awk -F'\t' '$3=="transcript"' > "${PREFIX}_contained.gtf" 2>/dev/null || touch "${PREFIX}_contained.gtf"

# Get transcript IDs
grep -o 'transcript_id "[^"]*"' "${PREFIX}_exact_matches.gtf" 2>/dev/null | \
    sed 's/transcript_id "//; s/"$//' | sort -u > "${PREFIX}_exact_match_ids.txt" || touch "${PREFIX}_exact_match_ids.txt"

# Step 4: Filter our output for high-confidence isoforms
echo "Step 4: Extracting high-confidence isoforms..."

# Create high-confidence set (exact matches only for strict ground truth)
if [ -s "${PREFIX}_exact_match_ids.txt" ]; then
    # Get transcripts
    awk -F'\t' 'NR==FNR{a[$1]; next} $3=="transcript"{match($0, /transcript_id "([^"]*)"/, m); if(m[1] in a) print}' \
        "${PREFIX}_exact_match_ids.txt" "$OUR_GTF" > "${PREFIX}_ground_truth.gtf"
    
    # Get exons
    awk -F'\t' 'NR==FNR{a[$1]; next} $3=="exon"{match($0, /transcript_id "([^"]*)"/, m); if(m[1] in a) print}' \
        "${PREFIX}_exact_match_ids.txt" "$OUR_GTF" >> "${PREFIX}_ground_truth.gtf"
    
    hc_count=$(wc -l < "${PREFIX}_exact_match_ids.txt")
else
    hc_count=0
    touch "${PREFIX}_ground_truth.gtf"
fi

# Step 5: Summary
echo ""
echo "=== Results ==="
total=$(grep -c 'transcript_id' "$OUR_GTF" 2>/dev/null | head -1 || echo 0)
echo "Total isoforms in output: $total"
echo "Exact matches to reference: $hc_count"
echo "Percentage validated: $(awk "BEGIN {printf \"%.1f\", ($hc_count/$total)*100}")%"
echo ""
echo "Files created:"
echo "  ${PREFIX}_cmp.stats - Comparison statistics"
echo "  ${PREFIX}_cmp.combined.gtf - Combined annotation with class codes"
echo "  ${PREFIX}_exact_matches.gtf - Exact reference matches"
echo "  ${PREFIX}_ground_truth.gtf - High-confidence ground truth"
echo ""
echo "View statistics: cat ${PREFIX}_cmp.stats"
echo "View exact matches: head ${PREFIX}_exact_matches.gtf"
