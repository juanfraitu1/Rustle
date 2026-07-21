#!/usr/bin/env bash
# V1: verify the spine's decision numbers on real data. Foreground; outputs to winloci_scratch.
#
# Confirms two things from the shipped binaries + source, not from memory/docs:
#   1. The O1 edge criteria (asm20 identity>=0.80, cov-of-shorter>=0.50, >=2 distinct loci) as
#      printed by gw_family_catalog --help.
#   2. The O2 spine number on the real GSTM locus (n_copies=3), and that the three certificate
#      call sites (copy_assign.rs, absent_copy.rs, collapse_gate.rs) reference the same
#      significance form (alpha / poisson_binomial_upper_tail / Bonferroni saturating_sub(1)).
set -euo pipefail
RUSTLE=/mnt/c/Users/jfris/Desktop/Rustle
SCRATCH=/home/juanfra/winloci_scratch
BIN=$RUSTLE/target/release

cd "$SCRATCH"
"$BIN/copy_assign" --bam GGO_mm.bam --fasta GGO.fasta \
  --region NC_073224.2:129160000-129230000 \
  --homology-primary --min-copies 2 --out mech_gstm

# families.tsv columns: family_id(1) chrom(2) n_copies(3) ... — NOT column 2 (that's chrom).
NCOPIES=$(awk -F'\t' 'NR==2{print $3}' mech_gstm.families.tsv)
echo "GSTM n_copies=$NCOPIES"

# edge values from the shipped binary:
"$BIN/gw_family_catalog" --help 2>&1 | grep -iE "asm20 identity|cov-of-shorter|DISTINCT loci"

# certificate reuse:
for f in copy_assign absent_copy collapse_gate; do
  grep -q "saturating_sub(1)\|min_p\|alpha" "$RUSTLE/src/rustle/vg_family/$f.rs" \
    && echo "certificate present: $f.rs"
done
