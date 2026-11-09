# PREREG row 1591: does exon-level sharing separate true from false paralog pairs better than exon coverage?

Written 2026-11-08, before the scoring script exists or any exon-sharing number is seen. Edge-level test on the
families stage's candidate pairs (no MCL). Advisor context: Soto et al. 2025 define families as SD98 genes that SHARE
EXONS; the SD98 half was tested today (row 1580 amendment 1: halves recall). This row tests the exon-sharing half, on
de novo loci.

## Prior work

r1569: O1 families share internal-exon length chains 11x a matched null (structure corroborates). r1570: structure
alone defines families far worse than sequence (.406/.580 vs Compara sequence .885/.804). The shipped edge gate
(`--min-shared-exon-frac 0.60`) asks the best record to cover >= 60% of the SMALLER member's exonic length with
exon-to-exon aligned bases: coverage, not exon correspondence.

## Scores (`bench/denovo_eval/exonshare.py`; per gated-PAF record, CIGAR `cg:Z` walked; loci = `PREFIX.fam.loci.gff3`)

- **COV** (the shipped gate's measure, recomputed): aligned bases with an exon on both sides / exonic length of the
  member with the smaller exonic length.
- **SHARE**: each exon of the member with fewer exons (ties: smaller exonic length) is projected through the alignment;
  it is shared when its projected bases overlap one exon f of the other member by >= 50% of its own length and >= 50%
  of f's. SHARE = shared exons / that member's exons.
- **SHARE_SS**: as SHARE, but a shared exon must also land its internal splice sites (exon boundaries that are not the
  member's outer ends) within 5 bp of an internal splice site of the other member. Single-exon members: SHARE_SS = SHARE.
  (Implementation, fixed before any run: a site lies between two bases, so a boundary base's mapped position x counts
  when x or x + 1 is within 5 bp of a site of the other member.)

A gene pair (`famtruth.map_copy` of each locus; both protein-coding, in the Compara universe, expressed per the run's
families GTF; different symbols, as `PAIRSCORE_SAME_SYMBOL=ignore`) takes the maximum of each score over its records.
Truth: Compara pairs at simiiformes among expressed genes (as `pairscore.py`). Recall = truth pairs passing / all
truth pairs; precision = truth / all scored coding pairs passing.

## Decision rule

Curves at thresholds 0.00, 0.05, ..., 1.00 for each score. R0 = recall of COV >= 0.60 (the shipped gate). For SHARE and
SHARE_SS: precision at R0 by linear interpolation between the two bracketing points of that score's curve (points
sorted by recall); gain = that precision - precision of COV >= 0.60.

- Dev (A119b, `r1581dev/d25_I2.8.fam.loci.*`, genes CHM13 CAT): the score (SHARE or SHARE_SS) with the larger gain is
  fixed as the verdict score and appended here before HG002 is scored. If neither gains > 0, stop: no held-out claim.
- Held-out (HG002, `r1580/hg002/plain_I2.8.fam.loci.*`, genes `DENOVO_HG002_GFF`): **PASS iff the verdict score's gain
  >= +0.03**; R0 outside that score's recall range = FAIL.

Reported, no bar: the full curves; gain within SD-born locus pairs (one CHM13 / HG002 SEDEF pair links them, the
`SdPairs::holds` rule) and within the rest; score medians of true vs false scored pairs (the chance level of each score
among aligned unrelated genes).

## What follows

- PASS: an exon-sharing edge rule in the families stage, registered next as a family-level arm (pairscore,
  precision at matched recall against the shipped gate), on top of `--sd-damp`.
- FAIL: the coverage gate already carries what exon correspondence adds; reported to the advisor as the answer on the
  exon-sharing half of Soto's definition.

## Setting chosen (2026-11-08, A119b dev, job 4453696, before HG002 is scored)

40,629 records scored; 16,528 gene pairs (1,927 truth) of 2,180 truth pairs. R0 = recall of COV >= 0.60 = .653 at
precision .778. SHARE at R0: .886 (**+0.108**); SHARE_SS: .878 (+0.100). **Verdict score: SHARE.** Descriptive: SD-born
pairs (2,138, 1,750 truth) COV>=0.60 .928/.789, SHARE +0.018; not SD-born (14,390, 177 truth) .126/.243, SHARE +0.177.
Medians true / false: COV .98 / .00, SHARE 1.00 / .00. Curves in `/anvme/workspace/v138da13-rustle/r1591/dev.4453696.out`.

## Outcome (2026-11-08, HG002 held-out, job 4453697): PASS

14,565 gene pairs scored (304 truth) of 414 truth pairs. R0 = recall of COV >= 0.60 = .309 at precision .292. SHARE at
R0: **.660, gain +0.368** (bar +0.03). PASS. SHARE_SS +0.370. Not an artefact of the interpolation: SHARE lies above
COV at every recall (.368 vs .229 at recall .42; .771 vs ~.31 at recall .27); the SHARE curve has a step between
thresholds 0.50 and 0.55 (.368/.418 -> .771/.268). Descriptive: SD-born pairs (307, 182 truth) COV>=0.60 .747/.632,
SHARE +0.061; not SD-born (14,258, 122 truth) .046/.107, SHARE +0.429. Medians true / false: COV .51 / .00, SHARE .50 /
.00. As on A119b, nearly all of the gain is in the pairs no SD holds: exon correspondence removes the domain-only and
partial-overlap false pairs that exon coverage admits.

Scope: edge level, recomputed coverage, no MCL. Next (registered separately): an exon-sharing edge rule inside the
families stage, scored family-level with pairscore against the shipped gate, with and without `--sd-damp`.
