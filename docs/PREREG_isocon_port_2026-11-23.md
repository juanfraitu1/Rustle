# PREREG: IsoCon's decision logic ported into Rustle, tested against IsoCon itself (row 1926, 2026-11-23)

Row 1925: IsoCon run on the reads of Rustle's tied groups counts expressed copies better than IsoCon with truth pools.
User: "my advisor's idea is to port in isocon's logic that makes it possible, not just run isocon".

## What was read in IsoCon 0.3.3 and what is ported (`bench/o2_families/recurrence_test.py`)
Source read: `hypothesis_test_module.py` (`statistical_test`, `arrange_alignments_new_no_realign`,
`raghavan_upper_pvalue_bound`, `get_correction_factor`), `functions.py` (`get_support`, `get_read_errors`,
`get_empirical_error_probabilities`), `isocon_statistical_test.py` (`stat_filter_candidates`).
- Ported as it is: support = reads carrying EVERY variant that separates a candidate from its nearest candidate; the
  null that all reads come from the neighbour; per-read probability of showing those variants by error from the read's
  OWN error counts (substitutions / 3, insertions / 4, deletions; never under |D| / length; an indel in a homopolymer
  counts by its length); Raghavan's bound on the weighted sum; the correction for the number of ways the variants
  could be placed; test only when the candidate has fewer than 5 x the neighbour's reads; alpha 0.01; candidates
  failing are merged into the neighbour and the tests repeated until nothing changes.
- Adapted to reference coordinates (not IsoCon's): candidates come from the reads' alignments at a locus of the tied
  group instead of all-against-all alignment and consensus; a difference seen in one read is dropped (IsoCon corrects a
  read towards its neighbours); reads are split by intron chain first and a chain or pattern contained in exactly one
  longer one is a shorter read of it (IsoCon's end-invariance; this is what removes its truncated candidates); a
  chained group is read at as many loci as it takes to cover its reads; a sequence needs >= 2 reads.
- Not ported: its alignment and consensus machinery, quality-value probabilities (planted reads have none).

## Development (said so the test is honest)
Built on seeds 1930-1932 of row 1924 in three steps, each after looking at their scores: (1) logic as above, one locus
per group, no isoform split: 0.17-0.42 exact; (2) split by intron chain + two-read floor: 0.35-0.79; (3) several loci
per chained group (TSPY's units were left without reads): 98 of 128 cells (0.766), control 40 of 42. Those seeds are
spent. Nothing was run on seeds 1933-1935.

## Test (frozen tool, the six variant sets and the control of row 1925; their scores for the other arms are known)
Same cells and truth as row 1925. Known there: hybrid (IsoCon on Rustle's groups) 72 of 126 (0.571), IsoCon with truth
pools 60 (0.476), every copy of the unit 63 (0.500), one per cell 32; control: hybrid 37 of 46, none over.
- T1: the port's exact rate >= the hybrid's (72 of 126) and above both constants.
- T2: control >= 0.80 exact and at most 2 cells over-counted.
- T3: above IsoCon with truth pools (60 of 126).
All three = the port replaces the idea of calling IsoCon, and is worth wiring as an opt-in column of the `o2` report
(sequences expressed among copies the reference cannot separate); wiring and any default are the user's call.

## Predictions (mine)
All three hold, with an exact rate near 0.7 (the development rate was 0.77; some of that is fitting to three seeds).
Limits known beforehand, unchanged: uniform simulated errors, no heterozygous variants (two alleles of one copy would
count as two sequences), sample = reference outside the planted variants, read names used to credit a sequence to a cell.

## Result (workspace `r1926_isocon_port`; cells with the count of expressed sequences exact)

| | development (seeds 1930-1932) | TEST (seeds 1933-1935) | confirmation (seeds 1936-1938, job 4496835) |
|---|---|---|---|
| ported test, with variants | 98 of 128 (0.766) | **115 of 126 (0.913)** | **101 of 126 (0.802)** |
| ported test, control | 40 of 42 | 46 of 46 | 43 of 44, none over |
| hybrid (IsoCon on Rustle's groups) | | 72 (0.571) | not run |
| IsoCon with truth pools | | 60 (0.476) | not run |
| every copy of the unit | 64 (0.500) | 63 (0.500) | 62 (0.492) |
| one per cell | 38 | 32 | 37 |

- **T1 HOLDS** (115 against the hybrid's 72 and the constants' 63 / 32). **T2 HOLDS** (46 of 46, none over).
  **T3 HOLDS** (115 against 60).
- The test rate (0.91) was above the development rate, so three never-used seeds were added after the test (not
  registered beforehand, said here): 0.80, control 43 of 44. **Quote 0.80, not 0.91**; the test seeds were favourable.
- Where it is still wrong (confirmation: over 16, under 9 of 126): mostly one sequence too many or too few in the large
  units (TSPY, 5-8 expressed copies). About 30% of the sequences it keeps hold reads of more than one true sequence
  (partial reads folded into a pattern that contains them), so the count is better than the read grouping under it.
- Prediction: all three hold (right), "near 0.7" (the confirmation says 0.8).
- **Verdict by the registered rule: the port replaces calling IsoCon.** No IsoCon install is needed, it runs in seconds
  on a planted set, and it is ahead of IsoCon itself on this count because genome alignments give it isoforms and read
  ends for free. Wiring as an opt-in column of the `o2` report and any default are the user's call.
- Limits, unchanged and untested: real HiFi errors, heterozygous variants (two alleles of one copy would count as two
  sequences), a sample that lacks or gains copies, real libraries (no truth). The credit of a sequence to a cell uses
  the read names.

**Wired 2026-11-23 (user: "yes lets wire it as an opt in column").** Stage `o2` with `RUSTLE_O2_SEQUENCES=1` (needs `--fasta`)
runs `recurrence_test.py` on the stage's own assembly, BAM and score table -> `PREFIX.o2.sequences.tsv`;
`o2_family_report.py --sequences` adds `tied_sequences` to both report tables (sequences read at that copy; family = sum).
A failure never fails the stage. Off by default.
End-to-end on real gorilla OR6737 through the driver (job 4496845, 799 s for the stage, exit 0; no truth): 325 groups,
392 sequences kept at 194 loci (130 loci with one sequence, 64 with 2-13), 95 tested, 282 with a non-reference
pattern; `tied_sequences` > 0 in 116 of 379 families.
