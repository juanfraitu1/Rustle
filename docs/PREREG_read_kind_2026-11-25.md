# PREREG r1950: can evidence across reads tell a one-origin read from a several-copy read? (2026-11-25)

Why: row 1949. For a near-tied read the right statement is the primary alone when the read has one origin copy (0.883)
and every record within 10 `ms` when it fits several (0.961); no single statement serves both (best 0.848), and knowing the
kind would give 0.952. User, 2026-11-25: look at whether anything across reads tells the kinds apart.
Related, closed: rules that MOVE a read off its primary (rows 1777, 1783, 1923). This row never moves a read: it decides
only how many of the read's own near-tied places to keep beside the primary.

## Data (HG002 on CHM13; both quarters have been looked at)
Truth, records, kinds as rows 1796 / 1949. Fitting quarter = residues 0-24 (`r1785_truth_repair/hg`, `r1796.../hg_seen/sel.bam`);
judging quarter = 25-49 (`r1796.../hg_held`). Across-read evidence from the whole sample: the tie table of row 1949
(`r1949_assign_two_questions/hg.ties.tsv.gz` + `.placed.tsv`, loci of `g.multicopy.gtf`).

## Evidence across reads: how often two loci tie
For loci A and B: `rate(A, B)` = reads tied between A and B (both within 10 `ms` of the read's best) divided by those reads
plus the reads placed alone at A plus the reads placed alone at B. Near 1: the sample's reads cannot tell the copies apart
(identical where reads cover them); near 0: most reads at either copy are placed without question.
A record's locus = the `g.multicopy.gtf` locus on its contig and strand it overlaps most.

## Statements (per molecule near-tied on `ms`, kinds `one` + `several`)
- baselines: primary; exact band; within 10 (row 1949).
- **G (no fitting):** primary + every record tied exactly with the best.  (= exact band; printed as the reference.)
- **T(tau):** primary + exact ties + every record within 10 whose locus has rate(primary's locus, its locus) >= tau.
  A record on no locus is kept only when it ties exactly. tau from {0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9}: the value
  with the highest pooled `equals` on the fitting quarter, then frozen.
- **X(tau):** as T but exact ties are also subject to the rate (a record is kept only at rate >= tau, exact or not);
  same grid, own tau.
- oracle (kind known): primary for `one`, within 10 for `several`.

## Readings, fixed now (judging quarter)
SOLVES IT: a fitted statement reaches pooled >= 0.90 with `one` >= 0.70 and `several` >= 0.90.
HELPS: pooled above within 10 (0.848) with disjoint Wilson intervals, short of the above.
NO: otherwise.
Printed, no bar: for records within 10 that do not tie exactly, the share truly in the origin set by rate bin (the
evidence's raw separating power); the same for exact ties.

Prediction (mine): HELPS, not SOLVES: T around 0.87-0.89; the rate separates the non-exact records (true share rising
with rate) but exact ties at low rate are mostly true copies anyway, so X is not better than T.

Workspace `r1950_read_kind`. Script `bench/famwise_eval/read_kind.py`.

## RESULT (2026-11-25, job 4504755; `r1950_read_kind/kind.tsv`)
A first submission (job 4504754) was cancelled before any output was read: the tie table's rows of one read are not
adjacent and the pair counter assumed they were; fixed, nothing else changed.

Judging quarter (fitting quarter in brackets), share of molecules whose statement equals the origin set:
| statement | `one` | `several` | pooled |
|---|---|---|---|
| exact band | 0.785 (0.816) | 0.821 (0.821) | 0.817 (0.820) |
| within 10 | 0.003 (0) | 0.961 (0.959) | 0.848 (0.834) |
| T, tau 0.7 (fitted) | 0.754 (0.793) | 0.825 (0.825) | 0.816 (0.821) |
| X, tau 0.05 (fitted) | 0.150 (0.136) | 0.839 (0.833) | 0.757 (0.742) |
| oracle, kind known | 0.883 | 0.961 | 0.952 |

**NO.** T equals the exact band (2,660 against 2,662 molecules) and is below within 10; X is worse. My prediction (HELPS,
0.87-0.89) was wrong.
Why: how often two loci tie does not order the records. Records 1-10 `ms` below the best, share truly in the origin set by
rate: under 0.02: 0.89 (453 records); 0.02-0.05: 0.22; 0.05-0.1: 0.72; 0.3-0.5: 0.39; 0.5-0.7: 0.007 (139 records);
0.9-1: 0.86 (246). The same pattern on both quarters, so it is not noise: each bin is a few families with their own answer
(a pair whose reads are mostly placed alone can still share a stretch where a read fits both; a pair that ties in most
reads can still differ under a given read). Exact ties are true copies at nearly every rate (2,064 of 2,064 at >= 0.9;
134 of 137 under 0.02).
Reading: whether a near-tied read has one origin or several is not a property of the pair of loci; it depends on which
stretch of the copies the read covers and on the individual's own sequence there. A per-pair statistic over the sample
cannot give it. Not tested: evidence per POSITION (the sample's own bases at the sites where the two copies differ under
the read); rows 1777 and 1783 used that to move reads and failed, its use for the width of the set is untested.
