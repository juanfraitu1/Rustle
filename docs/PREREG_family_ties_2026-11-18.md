# PREREG r1790: are the near-tied reads concentrated in a few families? (2026-11-18)

Written before any count of this row. Workspace `r1790_family_ties` (`family_ties.py`). No new alignment and no new truth
run: the minimap2 BAMs and the read tables of r1776 / r1784 (HG002 on CHM13) and r1775 (KB3781 on its own assembly).

## Why
User, 2026-11-18: is the multicopy assignment problem a big deal, when minimap2's primary is right on about 90% of the
decidable near ties and most of the rest cannot be decided? Over all family reads of HG002 on CHM13 the near ties are
5.7% (7,770 of 136,074) and the decidable primary errors 0.05%. That average can hide families in which most reads are
near-tied. This row measures the spread per family. It tunes nothing.

## Definitions, fixed now
- Family = a component of the truth pairs (`famtruth`, level simiiformes), as in `neutral_truth.py select`. A read belongs
  to the family of the member its minimap2 primary overlaps (same strand, exon overlap; the first family in component order
  when several match). Reads that match no member in the re-aligned BAM are counted as `unplaced` and reported.
- Near tie = stratum S2 of `neutral_truth.py` (best rival locus within 10 AS). Status of a near tie from `reads.tsv`:
  tie, decidable, no near-perfect hit, outside ("no home" = the last two).
- A family is counted when it has >= 50 selected reads in the sample.

## Readings, fixed now
1. CONCENTRATED when the smallest set of counted families that holds half of all near ties is at most 5% of the counted
   families. Otherwise DIFFUSE.
2. LOCALLY LARGE when at least 10 counted families have >= 20% of their reads near-tied. Otherwise LOCALLY SMALL.
3. Answer to the user's question: "a big deal for named families" needs both CONCENTRATED and LOCALLY LARGE on HG002;
   DIFFUSE and LOCALLY SMALL = the problem is small at family level as well. Mixed results are reported as mixed.
4. Reported, not judged: per family the share of near ties that are true ties / decidable / no home, the primary-wrong
   count, the 20 families with the most near ties, the same table on KB3781 (own reference, expected lower).

## Predictions
HG002: CONCENTRATED and LOCALLY LARGE (10 to 40 families at >= 20%). In those families true ties are the largest status.
KB3781: CONCENTRATED, fewer than 10 families at >= 20%.

## Limits stated beforehand
- HG002 reads are a 25% hash sample, so the >= 50 read floor means >= about 200 reads in the library; only expressed
  families of one library each are seen.
- Truth families come from the annotation; copies the annotation lacks are not families here.
- The read statuses carry the known defects of the neutral truth (`docs/NEUTRAL_TRUTH_PLACEMENT_2026-11-17.md`); the
  near-tie share itself does not depend on the truth.
