# PREREG: how much of the source transcript a member place reproduces (register 1981)

2026-11-26. User goal: ways multi-mapping reads can help discover multi-copy gene families.
Written before the coverage of any row was looked at.

## Question
Rows 1773 / 1915 / 1967 / 1968: the places reached by a locus's multimapping reads are a near-complete
candidate list for its relatives, at 0.66-0.72 of rows on a relative with today's `members` defaults.
Row 1773 left "a second criterion on the candidates" untested. The families stage accepts a pair of assembled
loci when half of the exons are shared. Same idea for a place no locus was assembled at:

    coverage(row) = max_bases of the row / exonic length of the source locus's most-reads transcript, capped at 1

Rule under test: keep a row when coverage >= 0.5. The 0.5 is the families edge rule's value, taken as is.

## Data
Cached tables, coverage never computed on them: A119b (`r1915_members_blocks/a1.members.tsv`), WTC11 r3, WTC11 r2
(`r1967`, `r1968`), scorer rows of `r1968_members_tandem/*.interspersed/rows.tsv`. Baseline = today's default
(score >= 0.3 of the best, interspersed-repeat filter, neighbour exemption). "At a relative" = classes 1-3 or
sequence-similar, as in rows 1967 / 1968. WTC11 r1 (never used for `members`) is run only if K1 and K2 hold.

## Claims
- K1: the share of rows at a relative rises by >= 0.10 on each of the three samples.
- K2: reach of recent paralogues and of pseudogenes each >= 0.9 x the baseline's on each sample.
- Description: the same at 0.25 and 0.75; rows on no assembled locus kept.

## Reading
K1 and K2: candidate discovery rule, to be confirmed on WTC11 r1, then the user's call.
K1 only: coverage is a confidence column, not a filter. Neither: closed.

## Prediction
K1 holds (about 0.85). K2 holds for paralogues and fails for pseudogenes (processed pseudogenes are fragments).

## Outcome (2026-11-26, `r1981_member_coverage/three_samples.txt`, `classes_by_cov.txt`)
| sample | at a relative, today | coverage >= 0.5 | recent paralogues | older paralogues | pseudogenes | copies of one source gene |
|---|---|---|---|---|---|---|
| A119b | 0.657 | 0.719 | 0.94 x | 0.77 x | 0.58 x | 0.92 x |
| WTC11 r3 | 0.715 | 0.806 | 0.92 x | 0.68 x | 0.59 x | 0.98 x |
| WTC11 r2 | 0.724 | 0.809 | 0.93 x | 0.76 x | 0.65 x | 1.00 x |

K1 FAILS (+0.06 to +0.09, bar +0.10; predicted 0.85: wrong). K2 FAILS for pseudogenes (predicted), holds for recent paralogues.
Under half coverage 15-20% of the rows are on nothing annotated, above it 2-3%; the wrong rows that remain sit on another
annotated gene (16-24%). WTC11 r1 not run. Closed: coverage is neither a discovery rule nor a useful filter.
