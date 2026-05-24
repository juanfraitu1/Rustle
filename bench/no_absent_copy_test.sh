#!/usr/bin/env bash
# No-Absent-Copy Ground Truth Test
#
# For a given multi-copy gene family BAM + reference GFF + assembled GTF files,
# computes per-paralog read evidence and recovery rates.
#
# Defines "expressed" = copies with >= MIN_READS total reads (primary + supplementary).
# Tests: what fraction of expressed copies does each tool recover?
#
# Usage:
#   bash bench/no_absent_copy_test.sh <BAM> <REF_GFF> <RUSTLE_GTF> <ST_GTF> [MIN_READS]
#
# Example (GOLGA8):
#   bash bench/no_absent_copy_test.sh \
#     bench/multi_copy_eval/golga8_region.bam \
#     bench/multi_copy_eval/ref_GOLGA8.gff \
#     bench/multi_copy_eval/golga8_fp.gtf \
#     bench/multi_copy_eval/stringtie_golga8/out.gtf \
#     5
set -euo pipefail

BAM=${1:?Usage: $0 <BAM> <REF_GFF> <RUSTLE_GTF> <ST_GTF> [MIN_READS]}
REF_GFF=${2:?}
RUSTLE_GTF=${3:?}
ST_GTF=${4:?}
MIN_READS=${5:-5}
OUTDIR=$(dirname "$RUSTLE_GTF")/no_absent_copy
mkdir -p "$OUTDIR"

echo "=== No-Absent-Copy Ground Truth Test ==="
echo "BAM:      $BAM"
echo "Reference: $REF_GFF"
echo "Rustle:   $RUSTLE_GTF"
echo "StringTie: $ST_GTF"
echo "Min reads (expressed threshold): $MIN_READS"
echo ""

# Step 1: Extract gene loci from reference GFF
LOCI_FILE="$OUTDIR/ref_loci.tsv"
awk -F'\t' '$3=="gene" {
    match($9, /ID=([^;]+)/, m)
    print $1 "\t" $4 "\t" $5 "\t" m[1] "\t" $7
}' "$REF_GFF" > "$LOCI_FILE"
N_LOCI=$(wc -l < "$LOCI_FILE")
echo "Reference paralogs: $N_LOCI"

# Step 2: Count primary and total reads per locus
READ_FILE="$OUTDIR/read_counts.tsv"
echo -e "gene_id\tchrom\tstart\tend\tstrand\tprimary_reads\tall_reads\tsupp_reads\tcategory" > "$READ_FILE"

while IFS=$'\t' read -r chrom start end gene_id strand; do
    prim=$(samtools view -F 2064 -c "$BAM" "${chrom}:${start}-${end}" 2>/dev/null || echo 0)
    all=$(samtools view -c "$BAM" "${chrom}:${start}-${end}" 2>/dev/null || echo 0)
    supp=$((all - prim))
    if [ "$all" -ge "$MIN_READS" ]; then
        if [ "$prim" -ge "$MIN_READS" ]; then
            cat="primary"
        else
            cat="multi_map_only"
        fi
    else
        cat="absent"
    fi
    echo -e "${gene_id}\t${chrom}\t${start}\t${end}\t${strand}\t${prim}\t${all}\t${supp}\t${cat}"
done < "$LOCI_FILE" >> "$READ_FILE"

N_EXPRESSED=$(awk -F'\t' 'NR>1 && $9!="absent" {n++} END{print n+0}' "$READ_FILE")
N_ABSENT=$(awk -F'\t' 'NR>1 && $9=="absent" {n++} END{print n+0}' "$READ_FILE")
echo "Expressed copies (>= $MIN_READS reads): $N_EXPRESSED"
echo "Absent copies (< $MIN_READS reads):     $N_ABSENT"
echo ""

# Step 3: Run gffcompare for each tool
RUSTLE_CMP="$OUTDIR/rustle_cmp"
ST_CMP="$OUTDIR/st_cmp"
gffcompare -r "$REF_GFF" "$RUSTLE_GTF" -o "$RUSTLE_CMP" 2>/dev/null
gffcompare -r "$REF_GFF" "$ST_GTF"     -o "$ST_CMP"     2>/dev/null

# Step 4: Parse per-locus recovery from tracking file (class "=" = exact intron-chain match)
# tracking file: col3=ref_tx_id (format "gene_id|tx_id" or bare gene_id), col4=class_code
# Gene ID extraction: strip trailing "|<tx>" suffix if present, then match any ID= value.
RUSTLE_RECOVERED="$OUTDIR/rustle_recovered.txt"
ST_RECOVERED="$OUTDIR/st_recovered.txt"

awk -F'\t' '$4=="=" {
    id=$3; sub(/\|.*/, "", id)          # drop "|tx_id" suffix
    sub(/^[^|]*\|/, "", id)             # or keep just the bare gene part
    # Re-extract: col3 is ref_gene_id as set by gffcompare (ID= value from GFF)
    id=$3; sub(/\|.*/, "", id); print id
}' "${RUSTLE_CMP}.tracking" 2>/dev/null | sort -u > "$RUSTLE_RECOVERED" || true
awk -F'\t' '$4=="=" {
    id=$3; sub(/\|.*/, "", id); print id
}' "${ST_CMP}.tracking" 2>/dev/null | sort -u > "$ST_RECOVERED" || true

# Locus overlap from loci file (any transcript at locus, col4 != "-")
RUSTLE_OVERLAP="$OUTDIR/rustle_overlap.txt"
ST_OVERLAP="$OUTDIR/st_overlap.txt"
awk -F'\t' '$4 != "-" {id=$3; sub(/\|.*/, "", id); print id}' \
    "${RUSTLE_CMP}.loci" 2>/dev/null | sort -u > "$RUSTLE_OVERLAP" || true
awk -F'\t' '$4 != "-" {id=$3; sub(/\|.*/, "", id); print id}' \
    "${ST_CMP}.loci" 2>/dev/null | sort -u > "$ST_OVERLAP" || true

# Step 5: Per-locus result table
RESULT_FILE="$OUTDIR/results.tsv"
echo -e "gene_id\tcategory\tprimary\tsupp\ttotal\trustle_exact\trustle_overlap\tst_exact\tst_overlap" > "$RESULT_FILE"

while IFS=$'\t' read -r gene_id chrom start end strand prim all supp cat; do
    r_exact=0; grep -Fxq "$gene_id" "$RUSTLE_RECOVERED" 2>/dev/null && r_exact=1
    r_over=0;  grep -Fxq "$gene_id" "$RUSTLE_OVERLAP"   2>/dev/null && r_over=1
    # Coordinate-overlap fallback: GFF3 refs with no exon records produce no gffcompare matches.
    # Fall back to: does any assembled transcript overlap this locus?
    if [ "$r_over" -eq 0 ]; then
        r_tx=$(awk -F'\t' -v c="$chrom" -v s="$start" -v e="$end" \
            '$1==c && $3=="transcript" && $4<=e && $5>=s' "$RUSTLE_GTF" | wc -l)
        [ "$r_tx" -gt 0 ] && r_over=1
    fi
    s_exact=0; grep -Fxq "$gene_id" "$ST_RECOVERED"     2>/dev/null && s_exact=1
    s_over=0;  grep -Fxq "$gene_id" "$ST_OVERLAP"       2>/dev/null && s_over=1
    if [ "$s_over" -eq 0 ]; then
        s_tx=$(awk -F'\t' -v c="$chrom" -v s="$start" -v e="$end" \
            '$1==c && $3=="transcript" && $4<=e && $5>=s' "$ST_GTF" | wc -l)
        [ "$s_tx" -gt 0 ] && s_over=1
    fi
    echo -e "${gene_id}\t${cat}\t${prim}\t${supp}\t${all}\t${r_exact}\t${r_over}\t${s_exact}\t${s_over}"
done < <(awk -F'\t' 'NR>1' "$READ_FILE") >> "$RESULT_FILE"

# Step 6: Summary table
# Determine if gffcompare produced any exact matches (requires exon records in ref GFF).
# If not, use coordinate-overlap as the primary recovery metric.
n_exact_r=$(awk -F'\t' 'NR>1 && $6>0 {n++} END{print n+0}' "$RESULT_FILE")
n_exact_s=$(awk -F'\t' 'NR>1 && $8>0 {n++} END{print n+0}' "$RESULT_FILE")
if [ "$n_exact_r" -gt 0 ] || [ "$n_exact_s" -gt 0 ]; then
    MATCH_METRIC="exact intron-chain match"
    R_COL=6; S_COL=8
else
    MATCH_METRIC="coordinate overlap (GFF3 ref — no exon records for exact matching)"
    R_COL=7; S_COL=9
fi

echo "=== Per-category recovery ($MATCH_METRIC) ==="
printf "%-16s  %8s  %12s  %12s\n" "Category" "Copies" "Rustle" "StringTie"
printf "%-16s  %8s  %12s  %12s\n" "--------" "------" "------" "---------"

for cat in primary multi_map_only absent; do
    n=$(awk -F'\t' -v c="$cat" 'NR>1 && $2==c {n++} END{print n+0}' "$RESULT_FILE")
    r=$(awk -F'\t' -v c="$cat" -v col="$R_COL" 'NR>1 && $2==c && $col>0 {n++} END{print n+0}' "$RESULT_FILE")
    s=$(awk -F'\t' -v c="$cat" -v col="$S_COL" 'NR>1 && $2==c && $col>0 {n++} END{print n+0}' "$RESULT_FILE")
    if [ "$n" -gt 0 ]; then
        printf "%-16s  %8d  %12s  %12s\n" \
            "$cat" "$n" \
            "$r/$n ($(( r * 100 / n ))%)" \
            "$s/$n ($(( s * 100 / n ))%)"
    else
        printf "%-16s  %8d  %12s  %12s\n" "$cat" 0 "n/a" "n/a"
    fi
done
echo ""
echo "=== Expressed copies only (>= $MIN_READS reads) ==="
n_exp=$(awk -F'\t' 'NR>1 && $2!="absent" {n++} END{print n+0}' "$RESULT_FILE")
r_exp=$(awk -F'\t' -v col="$R_COL" 'NR>1 && $2!="absent" && $col>0 {n++} END{print n+0}' "$RESULT_FILE")
s_exp=$(awk -F'\t' -v col="$S_COL" 'NR>1 && $2!="absent" && $col>0 {n++} END{print n+0}' "$RESULT_FILE")
echo "Rustle:    $r_exp / $n_exp expressed copies ($(( r_exp * 100 / n_exp ))%)"
echo "StringTie: $s_exp / $n_exp expressed copies ($(( s_exp * 100 / n_exp ))%)"
echo ""
echo "Key: expressed = copies with >= $MIN_READS total reads (primary + supplementary)"
echo "     primary   = >= $MIN_READS primary reads (mappable without multi-mapping)"
echo "     multi_map_only = >= $MIN_READS reads but <$MIN_READS primaries (only via multi-mapping EM)"
echo "     match     = $MATCH_METRIC"
echo ""
echo "Results saved to: $OUTDIR/"
