# PREREG row 1569: do O1's sequence families also share exon architecture? (a sequence-free second view)

## Why
Advisor: O1 is circular. O1 families are built from nucleotide alignment, and the CAT truth that scores them is also
alignment-based. A family view that uses no bases tests O1 without that loop: the lengths of a copy's internal exons,
in transcript order. Validation only; no default changes either way.

## Substrate (fixed before any score was computed)
A119b base run `~/runs/human_A119b/a1.*` (2026-10-29). Members: `a1.fam.copies.tsv` (1,565 multi-copy families,
7,987 copies; exons column, 0-based half-open). Null universe: `a1.fam.loci.gff3` (76,161 loci) minus every locus
whose span overlaps a family copy's locus span. Pair identity: `a1.fam.loci.gated.paf`, max col10/col11 over records
between the two loci's spans (names `chrom:start1-end`).

## Fingerprint and score
- Fingerprint: the representative's exons ordered 5'->3' by strand; internal exons = all but the first and the last
  (terminal lengths follow read ends). Testable = >= 3 exons (1,492 of 7,987 copies; 22,442 of 76,161 loci).
- Two internal exons match when their lengths differ by <= 3 bp.
- m(a,b) = longest in-order chain of matched internal exons (LCS; skips allowed on both sides);
  S = m / min(n_int(a), n_int(b)). A pair is ALIKE when S >= 0.5.

## Pair sets
- F: every unordered pair of testable copies in the same O1 family, spans not overlapping.
- F_lo: pairs of F with a PAF record and identity < 0.95. Reads seeded from secondaries at AS >= 0.98 x best can
  copy one locus's exons onto a near-identical partner; below 0.95 they cannot, so F_lo is the read-sharing control.
- N (null): for each pair in F, 20 pairs (a', b') of testable null-universe loci with n_int(a') = n_int(a) and
  n_int(b') = n_int(b) (n_int capped at 12), spans not overlapping, seed 1569.

## Bars
1. B1: ALIKE rate in F >= 3 x the rate in its matched N, one-sided Fisher p < 1e-6.
2. B2: the same in F_lo with its own matched N.
3. B3 (could structure alone define families?): over family copies plus null-universe loci with n_int >= 3, pairs (non-overlapping) whose
   internal-exon length tuples are identical (exact, 0 bp); fraction in the same O1 family >= 0.5.
PASS (the second view corroborates O1) = B1 and B2. B3 is scored separately: it decides only whether structure
could be a family DEFINITION.

## Descriptive (no bar)
Recall of B3 on O1 pairs with n_int >= 3; the "all members alike" share of families with >= 3 testable members vs
matched null families of the same sizes; alike rate by identity bin (< 0.80, 0.80-0.90, 0.90-0.95, 0.95-0.99, >= 0.99).

## Disclosed limits
Exon boundaries still come from the same read alignments (not the bases, but the aligner's placement). Coverage is
about 19% of copies; single- and two-exon copies (processed pseudogenes, fragments) are untestable. The null universe
can hold paralogs O1 missed, which makes the null rate higher (conservative for B1/B2, deflates B3).

## Prediction
B1 PASS, B2 PASS with a smaller ratio, B3 FAIL (identical architecture is common among unrelated multi-exon genes).

## Result (2026-11-07, `bench/structfam/structfam.py`, output `~/runs/structfam_r1569/a1.structfam.tsv`)
**PASS (B1 + B2); B3 PASS as registered but carried by near-identical copies (post hoc).** Testable: 1,620 of 7,987
copies (the GFF count of 1,492 missed 211 renamed tids), 401 families with >= 2 testable copies; null 19,756 loci.
- B1: 4,745 within-family pairs ALIKE 0.898 vs matched null 0.082 (94,900 pairs), 11.0x, p ~ 0. PASS.
- B2 (identity < 0.95, read sharing impossible): 919 pairs 0.738 vs null 0.161, 4.6x, p 3e-306. PASS.
- B3: 1,902 exact-tuple pairs, 1,691 in the same O1 family, precision 0.889. PASS as registered. Post hoc split:
  1,461 of the 1,691 are same-chromosome copies at identity >= 0.99, where shared reads can copy the exons. Outside
  >= 0.95 (identity < 0.95 or no PAF record): 128 same-family of 327, precision 0.39. Structure alone does not define
  families; it corroborates them.
- Descriptive: ALIKE by identity < 0.80 0.56 (52 pairs), 0.80-0.90 0.76, 0.90-0.95 0.72, 0.95-0.99 0.88, >= 0.99
  0.998; MCL-transitive pairs with no PAF record 0.79 (803). Exact-tuple recall on O1 pairs with n_int >= 3: 0.641.
  Families with >= 3 testable members where every pair is ALIKE: 107/167 (0.64) vs matched null families 0.0012.
Prediction right on B1/B2 (B2 ratio smaller), wrong on B3 as registered, right in substance after the post hoc split.
