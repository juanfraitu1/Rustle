# PREREG row 1561: flag satellite-like families by copy geometry alone (bench/satgeom.py)

Written before any number of the rule was computed.

## Input and truth
- Copies: /anvme/workspace/v138da13-rustle/asgsweep/gw_off/a1.fam.copies.tsv (1,565 families; every row counts,
  no `member_status` filter, as tools/annotated_families.py counts them).
- Truth: `tools/annotated_families.py COPIES ~/data/denovo_eval/chm13v2.0_censat_v2.1.bed --majority-in` =
  families with strictly more than half their copies overlapping a censat v2.1 element: 550 of 1,565
  (reproduced 2026-11-06; register row 1559's 618 counted the majority class, its PREREG gives 550).
  Note: the censat BED includes `ct` transition-zone elements (multi-Mb), so truth includes families in
  pericentromeric transition zones, not only satellite arrays.

## Rule (fixed, nothing tuned)
A family is satellite-like iff ALL of:
- >= 6 copies (rows);
- all copies on one chromosome;
- span = max(end) - min(start) <= 5,000,000 bp;
- median copy length (end - start, genomic span of the copy) <= 2,000 bp.

## Bars
Precision >= 0.90 and recall >= 0.80 against the 550. Both must hold to PASS.

## Expected failure modes (named now)
- False positives: tandem gene arrays (PCDH, UGT1A, OR clusters) are one-chromosome, compact, multi-copy; most
  have spliced copies longer than 2 kb, but single-exon OR copies (~1 kb) can pass. Every FP is listed with its
  gene_id so the class is visible.
- False negatives: satellite families with copies spread over several chromosomes (shared HSat/rDNA arrays on
  the acrocentrics) fail the one-chromosome clause; families with < 6 copies fail the count clause. FN counts per
  failing clause are reported.

## Prediction
Precision passes (>= 0.9), recall FAILS (one-chromosome and >= 6 clauses cut multi-chromosome and small
satellite families): recall 0.4-0.7.
