# PREREG: two tiers of the multi-copy mode, by whether a primary alignment carries the transcript (row 1952, 2026-11-25)

Diagnosis (workspace `false_isoforms_2026-11-25`, seeds 1964 / 1965, a description): at 60 reads per copy the tied
secondary alignments add no exact isoform to the mode (197 = 197) and add 50-56 chains that are no planted isoform of
their copy; 53-78 of the mode's non-exact chains per set are absent from the primaries-only assembly. FLAIR avoids this by
construction: its isoforms come from primary alignments and tied reads are shared among isoforms that exist. The user
agreed (2026-11-25) to test a TIER, not a filter: nothing is dropped.

## The tier (`bench/o2_families/support_tier.py`, from `tools/read_provenance.py`'s transcript table)
`primary` = at least one primary alignment has the transcript's intron chain (exactly, or after the junction replacement
the provenance tool infers), or for a single-exon transcript lies inside it. `tied_only` = none.

## Data
Planted human named families, FRESH seeds 1967, 1968, 1969 (ragged ends, today's defaults): all expressed (E) and half
silent (S) at 60 reads per copy; all expressed at 3 reads per copy. Real: gorilla OR6737, the mode of row 1937.

## Claims
- **C1, precision:** in each of the six sets at 60 reads the `primary` tier's isoform precision is at least 0.10 above the
  whole mode's, with isoform recall not lower by more than 0.01.
- **C2, false units:** in each half-silent set the `primary` tier has fewer false units than the whole mode.
- **C3, the other tier earns its place:** at 3 reads per copy, in each seed at least 25% of the mode's exact isoforms are in
  the `tied_only` tier (so dropping that tier would cost the low-coverage recall).
- **R, real sample:** on OR6737 the `primary` tier's junction precision at the named families is not lower than the whole
  mode's. The annotated copies of each tier are reported, no claim.
- Reported for every set: transcripts, exact isoforms, isoform recall and precision, units and false units of both tiers.

## Predictions (mine)
C1 holds (I expect about 0.80 against 0.65). C2 holds. C3 holds, with 40-50% of the exact isoforms at 3 reads in the
`tied_only` tier. R holds. If C3 fails, the tier is simply a filter and `RUSTLE_SEED_SPLICED=support` would do.

## Result (job 4504820; seeds 1967 / 1968 / 1969)

| | whole mode | `primary` tier | `tied_only` tier | Rustle primaries only | FLAIR | StringTie |
|---|---|---|---|---|---|---|
| isoform precision, all expressed, 60 reads | 0.654, 0.654, 0.659 | 0.811, 0.824, 0.804 | 0 of 58, 62, 54 exact | 0.831, 0.828, 0.812 | 0.769, 0.736, 0.729 | 0.769, 0.777, 0.766 |
| isoform recall, all expressed, 60 reads | 0.824, 0.824, 0.824 | 0.824, 0.824, 0.824 | 0 | 0.820, 0.824, 0.816 | 0.808, 0.803, 0.820 | 0.669, 0.669, 0.657 |
| isoform precision, half silent | 0.452, 0.417, 0.488 | 0.609, 0.621, 0.630 | 0 of 59, 83, 48 exact | 0.646, 0.620, 0.620 | 0.585, 0.597, 0.636 | 0.621, 0.590, 0.601 |
| false units, half silent (of 47, 49, 46) | 29, 34, 23 | 13, 13, 9 | 22, 28, 19 | 8, 12, 7 | 14, 16, 9 | 8, 15, 8 |
| exact isoforms at 3 reads (of 239) | 166, 162, 168 | 142, 140, 146 | 24, 22, 22 | 83, 70, 72 | 37, 49, 41 | 79, 64, 69 |
| isoform precision at 3 reads | 0.755, 0.726, 0.740 | 0.835, 0.819, 0.811 | 0.480, 0.423, 0.468 | 0.954, 0.921, 0.935 | 0.841, 0.961, 0.872 | 0.760, 0.744, 0.750 |
| expressed units found at 3 reads (of 108) | 103, 103, 103 | 94, 97, 101 | 34, 25, 28 | 66, 65, 64 | 40, 42, 42 | 84, 75, 81 |

Gorilla OR6737 (named families: annotated copies with a spliced transcript / junction precision / spliced transcripts):
whole mode 37 / 0.67 / 176; `primary` tier 35 / 0.68 / 125; `tied_only` tier 19 / 0.64 / 51. Genome-wide the tier split is
70,718 transcripts with a primary alignment and 2,455 without (3.4%).

- **C1 HOLDS in all six sets:** precision +0.14 to +0.20 at identical recall. At 60 reads the `tied_only` tier holds NOT
  ONE exact isoform (0 of 364 transcripts over the six sets): everything the tied secondaries add there without a primary
  is a chain the copy does not express.
- **C2 HOLDS:** 13, 13, 9 false units against 29, 34, 23.
- **C3 FAILS:** at 3 reads the `tied_only` tier holds 24 of 166, 22 of 162 and 22 of 168 exact isoforms (13-15%, bar 25%),
  at precision 0.42-0.48. My prediction (40-50%) was wrong. The mode's low-coverage recall does NOT come mainly from
  transcripts built on secondaries alone: the `primary` tier by itself has 140-146 exact isoforms against 70-83 for the
  primaries-only assembly. Tied alignments help mostly by adding SUPPORT to chains a primary alignment carries (they lift
  them over the read floor), which is the way FLAIR's paralog option uses them.
- **R HOLDS:** 0.68 against 0.67. But this measure cannot see the problem: a paralog's isoform laid onto a copy has that
  copy's annotated junctions, so the `tied_only` tier scores 0.64 on real data while it is 0% exact on planted reads.
  Junction precision at annotated copies is not evidence that a tied-only transcript is expressed.
- **What it gives.** The `primary` tier alone: isoform precision 0.80-0.82 at full recall at 60 reads (FLAIR 0.73-0.77,
  StringTie 0.77), and at 3 reads 0.59-0.61 recall at 0.81-0.84 precision (a FLAIR with two supporting reads and its
  paralog option, row 1948: 0.46-0.51 at 0.83-0.86). The `tied_only` tier costs 13-15% of the low-coverage exact isoforms
  and 2-9 units if left out.
- Not decided, not wired: the tier is a tag made after the run by `support_tier.py` from the provenance table. Whether the
  tag becomes part of stage `o2`, and whether any product should show only the `primary` tier, is the user's call.

## Addendum, written before running it (2026-11-25): three layers, and what stays unexplained
User: "an option for primaries only, and another layer for near tied multimappers to show how much it improves or what is
still left unexplained". The `primary` tier is split by how many primary alignments carry the chain (the provenance
table's count; the assembler's read floor is 2):
- `layer "primary"`: at least 2 primary alignments. The transcript stands without any tied read.
- `layer "tied_support"`: exactly 1. A primary alignment carries the chain and tied reads lifted it over the floor.
- `layer "tied_only"`: none.
Products: `PREFIX.o2.primary.gtf` (the first layer), `PREFIX.o2.tied.gtf` (the other two, tagged), `PREFIX.o2.layers.tsv`
(transcripts, loci and labels per layer; and of the reads with a rival within 10 `ms`, how many a transcript of each
layer carries and how many no transcript carries, by reason).
A DESCRIPTION on the nine sets of this row (no decision rule; the split into three is new and its numbers were not
looked at). What I expect: at 3 reads per copy `tied_support` holds most of what the tied reads add, about 60 of the
140-146 exact isoforms of the two primary-carried layers, at a precision between the other two layers; at 60 reads it is
nearly empty.

### Result of the addendum (job 4504877; the nine sets of this row, and OR6737)
| | `primary` | `tied_support` | `tied_only` | primaries-only ASSEMBLY, for reference |
|---|---|---|---|---|
| exact isoforms at 3 reads (seeds 1967 / 1968 / 1969) | 82, 69, 72 | 60, 71, 74 | 24, 22, 22 | 83, 70, 72 |
| isoform precision at 3 reads | 0.95, 0.93, 0.94 | 0.71, 0.73, 0.72 | 0.48, 0.42, 0.47 | 0.95, 0.92, 0.94 |
| expressed units found at 3 reads, cumulative (of 108) | 65, 63, 65 | 94, 97, 101 | 103, 103, 103 | 66, 65, 64 |
| transcripts at 60 reads, all expressed (exact) | 239 (197), 236 (197), 243 (196) | 4 (0), 3 (0), 2 (1) | 58 (0), 62 (0), 54 (0) | 236 (196), 238 (197), 240 (195) |
| false units, half silent, cumulative | 8, 11, 6 | 13, 13, 9 | 29, 34, 23 | 8, 12, 7 |

- The expectation held: at 3 reads `tied_support` holds 60-74 exact isoforms (as many as the `primary` layer) at
  precision 0.71-0.73, between the two others; at 60 reads it is nearly empty (2-12 transcripts).
- **The `primary` layer reproduces the primaries-only assembly** (82 / 69 / 72 exact against 83 / 70 / 72; 0.93-0.95
  precision either way) without a second run: it is the primaries-only product.
- **What the near ties add, in one line per coverage:** at 3 reads per copy they double the exact isoforms through
  `tied_support` (+60 to +74, precision 0.72) and add 22-24 more through `tied_only` (precision 0.45); at 60 reads they
  add nothing correct.
- **What stays unexplained** (reads with a rival within 10 `ms`, by where their primary alignment sits). Planted, seed
  1967, all expressed, 60 reads: 5,078 near-tied reads, 4,488 in a `primary`-layer transcript, 50 in the two tied layers,
  287 unspliced inside an exon, and 253 carried by no transcript (108 non-canonical junction, 84 unspliced with no model,
  the rest under the floor or removed by a rule). Gorilla OR6737: 16,275 near-tied reads; 10,455 in a `primary`-layer
  transcript (64%), 581 in `tied_support`, 66 in `tied_only`, 1,030 unspliced inside an exon, and 4,143 carried by no
  transcript (25%: 2,058 unspliced with no model there, 1,049 under the read floor, 900 with a non-canonical junction, 136
  removed by a later rule). On OR6737 the layers hold 68,674 / 2,044 / 2,455 transcripts; `tied_support` adds 1,010 loci
  that have no `primary`-layer transcript and `tied_only` 1,240 more.
- Reading note: a near-tied read is counted where its PRIMARY alignment sits. A `tied_only` transcript is built from
  secondary alignments of reads whose primaries are at another copy, which is why few reads are "carried" by that layer.
- Wired into stage `o2` (same switch `RUSTLE_O2_SUPPORT`): `PREFIX.o2.primary.gtf`, `PREFIX.o2.tied.gtf`,
  `PREFIX.o2.layers.tsv`; checked by a real run of the stage and the full test suite.

### What the 25% unexplained near-tied reads are (description, gorilla OR6737; workspace `unexplained_unspliced_2026-11-25`, jobs 4504913, 4504915; no decision rule)
**More than half are explained through their tie.** Of the 4,143 near-tied reads that no transcript carries at their
primary alignment, 2,287 (55%) have a tied alignment on exons of an assembled transcript elsewhere: under the read floor
866 of 1,049, non-canonical junction 605 of 900, removed by a later rule 131 of 136, unspliced with no model 685 of 2,058.
Counted that way, 1,856 of 16,275 near-tied reads (11%) are carried by no transcript at any tied place.
**The 2,058 unspliced reads with no model**, by where the primary sits (RefSeq): inside a gene but not on its exons 794
(39%), nothing annotated 736 (36%), an annotated single-exon gene 294 (14%), an exon of a spliced gene 207 (10%), an
annotated pseudogene 27 (1%). They form 1,243 loci: 987 hold one read, 52 hold five or more (582 reads, 28%).
- The clear class: 207 reads sit on single-exon genes and tie with a SPLICED alignment at an annotated gene. The largest
  loci are intronless copies of ribosomal-protein genes annotated as single-exon protein-coding `LOC` genes (tied with
  RPL21: 82 reads, RPL12: 34 + 13, RPL11: 26, RPL9: 20; also MICU2: 20). This is the processed-copy case: a read of the
  spliced parent ties on `ms` with its intronless copy (introns cost nothing on `ms`) and the aligner's tie-break puts the
  primary on the copy. Whether those copies are expressed themselves cannot be told from these reads.
- The bulk, about 1,350 reads, ties between two unspliced places that are intronic or unannotated (550 intron to intron,
  450 unannotated to unannotated, 210 + 137 mixed): repeats and unspliced material, mostly singletons. No model there is
  what a read floor is for.
- Only 27 reads sit on annotated pseudogenes, so "expressed pseudogenes left on the table" is not what this group is on
  this sample (the RefSeq annotation calls the ribosomal-protein copies genes, not pseudogenes).

### `ms` against AS for those reads (description, job 4504918, `scores.py`; primary against its best-scoring secondary)
| group | reads | `ms`: exact tie | `ms`: primary ahead by 1-10 | `ms`: primary behind | AS: primary higher | AS equal | AS: primary lower | median AS difference |
|---|---|---|---|---|---|---|---|---|
| all near-tied reads no transcript carries | 4,143 | 32% | 54% | 14% | 56% | 28% | 16% | +3 |
| unspliced, no model there | 2,058 | 38% | 54% | 8% | 62% | 31% | 7% | +5 |
| under the read floor | 1,049 | 32% | 54% | 13% | 53% | 33% | 14% | +3 |
| non-canonical junction | 900 | 17% | 52% | 31% | 42% | 16% | 42% | 0 |
| removed by a later rule | 136 | 31% | 68% | 1% | 68% | 31% | 1% | +5 |
| intronless copy tied with its spliced parent | 207 | 69% | 28% | 4% | 97% | 3% | 0% | +120 |

- For most of these reads the two scores tell the same story: a small lead or none (median AS difference 3 to 5).
- **The intronless-copy class is where they disagree completely.** On `ms` 69% are exact ties and 97% are within 10. On AS
  the intronless copy is ahead in 97%, by exactly 24 per intron of the parent's alignment (1 intron +24, 3 introns +72,
  5 introns +120, 6 introns +147, 11 introns +270; median +120). AS charges for introns and `ms` does not, so AS makes the
  intronless copy look like a clear winner for a read that fits the spliced parent equally well base for base.
- In 14% of all these reads (31% of the non-canonical group) the primary has a LOWER `ms` than a secondary: minimap2
  orders hits on an adjusted segment score, which the `ms` tag (the unadjusted one) does not always reproduce.
- Not checked: how many of the parent alignments of that class fall outside the assembler's AS-based seeding pool (AS >=
  0.98 x best); with a lead of 24 per intron most of them should.
