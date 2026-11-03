# PREREG O2 CONTESTED CENSUS (register 1508, 2026-11-03)

Why: O2 judges exact AS ties only (`as_margin` 0). The user's scoping rule (2026-11-03): a read belongs to O2 only when
minimap2 itself was uncertain; reads it placed uniquely have no business in O2; exact ties on both scores with the same
alignment shape are a K = 0 wall. The population worth the certificate is the near ties: a rival close on AS **or on ms**
(minimap2 ranks primaries and secondaries on `ms`, the max-scoring segment, not on AS), including the processed-copy class
where a spliced primary faces an unspliced rival with a higher AS (r1477: 6-10% of family-scope reads). ms has no naive
rule like the 10-AS Eichler rule, and whether ms and AS share a scale was unknown, so the ms window must come from the
measured distribution, not be assumed.

Change: `as_table` (commit f9df72c3) also records best/second ms, the primary's ms, MAPQ and splice shape, and the best-AS
non-primary record's AS, ms, splice shape and chromosome (columns 7-15; columns 1-6 unchanged). One pass per BAM, no new
alignment.

Census (descriptive, no verdict): genome-wide tables for A119b, HG002 and OR6737 (job 4421419,
`/anvme/workspace/v138da13-rustle/o2contested/*.molecules.tsv`), molecules with a primary and a rival, cross-tabulated by
AS margin (primary AS - rival AS: < 0, 0, 1-10, 11-50, > 50) x ms margin (same bins) x splice-shape pair
(spliced/unspliced primary x spliced/unspliced rival) x MAPQ (< 60 / 60). Also the per-record ms - AS gap by splice shape
(the scale check).

Pre-declared reading of the census (fixed before the genome-wide tables exist):
1. Scale: ms and AS are the same DP score when the alignment is unspliced (slim A119b/HG002: unspliced primaries ms - AS
   mean 1.1 / 1.0); spliced records carry AS penalties that ms does not (gap mean 176 / 153). So margins are on one scale
   per record and a 10-point window means the same on both, but the ms margin and the AS margin of one molecule differ
   whenever primary and rival differ in splice shape. The contested set is therefore defined on BOTH scores:
   `AS margin <= 10 OR ms margin <= 10`, negative margins included (rival better on that score).
2. The O2 target set for the proof experiment (next register row) = the contested set minus the K = 0 wall (AS margin = 0
   and ms margin = 0 and identical splice shape, same chromosome). Molecules with no rival, or a rival > 10 on both scores,
   are out of scope.
3. Sub-populations to be reported with their sizes, in the order they will be tested: (a) AS near tie, ms tie; (b) spliced
   primary vs unspliced rival (either score in window); (c) both spliced, different junction sets; (d) AS tie, ms clearly
   separated (> 50): minimap2 was sure on ms; these are reported but are candidates for OUT of scope; (e) MAPQ < 60 with
   both margins > 10 (minimap2 uncertain for another reason): reported only.
No bar: this row sizes populations; the proof row that follows carries the bar (HG002 truth: assigned fraction and
agreement vs the minimap2 primary and vs the 10-AS rule, per sub-population).

## RESULT r1508 (2026-11-03, job 4421419; `~/o2contested/census_r1508.txt`, tables in the workspace `o2contested/`)
Scale confirmed genome-wide: unspliced primaries ms - AS 0.4 to 1.8; spliced 175 to 221. Contested (either margin <= 10),
molecules with a primary and a rival: A119b 484,805 (8.4% of with-rival, 2.3% of all), HG002 1,375,995 (8.7% / 3.7%),
OR6737 68,285 (5.1% / 1.6%). Both margins > 10 ("placed"): 91-95% of with-rival; MAPQ < 60 among them 0.7-3.5%.
K = 0 wall (AS 0, ms 0, same shape, same chromosome): A119b 16.0%, HG002 3.2%, OR6737 5.8% of contested, i.e. the wall is
small and most of the contested set is informative. Largest class everywhere: spliced primary vs unspliced rival whose AS is
HIGHER while ms favours the primary by > 10 (the processed-copy class): 33% / 72% / 68% of the informative set; the rival's
AS lead is 11-200 on ~88%, 1-10 on 9-17%; MAPQ < 60 on 11-27%. Next: near tie on both scores, same shape (A119b 31%, HG002
5%, OR6737 12%); both spliced, rival higher AS, ms > 10 for the primary (10% / 6% / 8%). "AS tie, ms > 50" is 838 / 8,151 /
706 molecules: out of scope as pre-declared. Reading rules 1-3 stand; the proof row defines the target set as the
contested set minus the wall and tests, per class, on HG002 truth.
