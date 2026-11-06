# PREREG row 1562: flag satellite families by copy SEQUENCE alone (bench/satseq.py)

Written before any number of the rule was computed. Row 1561 (geometry) FAILED: satellite families are
6.7 kb two-exon loci, often multi-chromosome. The one thing a satellite copy cannot hide is its sequence: a
tandem array is internally periodic, so most of its k-mers recur within the same copy.

## Rule
For each copy sequence in PREFIX.fam.copies.fa: dup15 = fraction of 15-mer positions whose 15-mer occurs at
least twice within that copy. A family is flagged iff the median dup15 over its copies >= 0.5.
GC content is reported alongside but is not in the rule (HSat1A is AT-rich; other satellites are not).

## Truth and bars (row 1561's)
Truth = censat v2.1 majority-in, 550 of 1,565 (tools/annotated_families.py --majority-in). Precision >= 0.9,
recall >= 0.8. The rule's threshold is fixed here; no re-tune after the numbers.

## Prediction
Recall PASS (0.85-0.95: HSat1A, HSat3, HOR and bsat copies are all tandem; the misses are `ct`-zone families
that are real genes near the satellite). Precision risk: genuine tandem-repeat-bearing genes (ZNF arrays, mucins
VNTR exons, PCDH) could exceed 0.5; predicted precision 0.85-0.95, so the precision bar is the one at risk.
## Use if it passes
RUSTLE_EXCLUDE_SELFREPEAT=1 at assign: the same skip-list hatch as row 1559, computed from the catalogue itself,
no BED. The BED path stays for annotated genomes.
