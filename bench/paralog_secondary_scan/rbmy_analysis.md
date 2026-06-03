# RBMY1 — copy-level assignment validation on a testis Y gene family

**Question:** does the copy-assignment method actually work on a multi-copy
testis gene family — i.e. does it recover multiple expressed copies, or collapse
them? RBMY1 (RNA-binding motif protein, Y chromosome, family 1) is an ideal,
biologically-known test: a testis-specific Y ampliconic family.

## The locus

A 6-copy tandem array on chrY, `NC_073248.2:19,602,754–19,730,926` (~14 kb per
copy, − strand, ~23 kb spacing), annotated LOC129530238–243 with descriptions
"RNA-binding motif protein, Y chromosome, family 1 member B-like / F-J-like".

| label | gene | position |
|---|---|---|
| c1 | LOC129530243 | 19,602,754–19,616,644 |
| c2 | LOC129530239 | 19,625,715–19,639,621 |
| c3 | LOC129530240 | 19,648,670–19,662,577 |
| c4 | LOC129530238 | 19,671,638–19,685,525 |
| c5 | LOC129530241 | 19,694,606–19,708,531 |
| c6 | LOC129530242 | 19,717,578–19,730,926 |

**Pairwise identity:** the core c2–c5 are ~99.8% identical; c1 ~97.3%; c6 the
outlier at ~93–96%. A genuine identity gradient, as expected for a young amplicon.

## Reads, multimapping, and resolution

- **87 primary reads → 599 total alignments** (512 secondary): each read places at
  ~7 positions across the 6 near-identical copies. This is exactly the
  multimapping regime that naive tools mishandle (drop secondaries → miss copies;
  or count them all → inflate).
- Despite that, the reads **resolve to individual copies**. Primary placements
  distribute across the array:

| copy | primary reads | transcripts | coverage | capacity_confidence |
|---|---|---|---|---|
| c1 LOC243 | 8 | 2 | 7.9 | 1.000 |
| c2 LOC239 | 10 | 4 | 12.1 | 1.000 |
| c3 LOC240 | 14 | 3 | 10.9 | 1.000 |
| c4 LOC238 | 7 | 1 | 6.8 | 1.000 |
| c5 LOC241 | 45 | 6 | 38.0 | 1.000 |
| c6 LOC242 | 2 | 0 | 0.0 | — (below threshold) |

**5 of 6 copies are assembled, each anchored (`capacity_confidence 1.000`).** The
6th (c6) has only 2 reads — below the assembly threshold, not a mis-assignment.

## What it validates

1. **No collapse.** A multi-copy testis family is recovered *as* multi-copy —
   reads spread across the array (8/10/14/7/45/2) rather than piling on one copy.
   So the method is not systematically mis-assigning everything to a single copy.

2. **Resolution tracks identity, not the tool.** The RBMY core is ~99.8% identical
   and still resolves (the copies have enough copy-specific positions). Only a
   *near-perfect* duplicate — DAZ1/DAZ3 at 99.97%, ~10× closer — becomes
   genuinely non-identifiable, and there the tool **flags** the unanchored copy
   (`low_confidence`) instead of guessing.

3. **Biologically consistent.** RBMY is testis-specific; recovering 5/6 copies as
   expressed in a testis sample is the expected result. Combined with DAZL
   (cov 442) and DAZ1 (cov ~36), the Y/germ-cell gene families are clearly
   expressed and correctly assigned here.

Figure: `rbmy_analysis.png`. Data: `/tmp/rbmy/` (per-copy stats, identity matrix,
VG attribution).
