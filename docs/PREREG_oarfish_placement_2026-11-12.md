# PREREG ABUNDANCE-EM PLACEMENT, oarfish (register 1648, 2026-11-12)

Why: user request; the second "alignment-free quantifier" of the 2026-11-12 option table. What it is: oarfish 0.11.0
read mode maps reads to transcripts with its own mapper (rammap, not minimap2) and resolves multi-mapping with an EM
over transcript abundances. So it is alignment plus an abundance prior: the r1412 class (EM overconfident at 0.99, dev
only), here on real reads with truth.

Samples, truth, reference: as r1647 (KB3781 r1642 set, HG002 r1644 set; `r1647_kallisto/{kb,hg}.ref.fa` with OTHER).
Abundance input: every primary read (`-F 0x904`) in the sample's panel regions (KB3781 `o2mstie/all.bed`, HG002
`r1644_hg002_unspliced/panel.bed`), not the contested reads alone.
Run: `oarfish --reads REGION.fa --annotated REF.fa --seq-tech pac-bio-hifi --write-assignment-probs`, defaults otherwise.
RULE (`bench/o2_contested/oarfish_place.py`): sum a read's posterior over all sequences of one locus (OTHER is its own
class); assign to the top locus iff its sum >= 0.99 and it is not OTHER.

Scorer, registered up front (r1647's two post hoc fixes): pick right = its span overlaps a truth copy
(`EXON_PLACE_RIGHT=locus`); primary right = the primary RECORD overlaps a truth copy (`EXON_PLACE_PRIMARY`, spans from
the region extraction). Validity = mechanics only (pilot parsed; done).

Bars, per sample, both must pass: 1. reach >= 30% of decidable molecules; 2. tool right / primary wrong >= 2 x the
reverse, first count >= 10; 3. negatives assigned <= 10%; 4. wall <= 30 min on 64 threads.
EM isolated (readout, the r1412 question): right / wrong on molecules the r1647 minimap2 control left as a tie on the
same reference: those are the reads abundance decides.
Stated bias: OTHER transcripts outside the panel regions receive almost no reads, so the EM leans to catalog copies; a
negatives FAIL is not interpretable without a whole-sample run.
Expected: bar 2 fails; on the tie subset the wrong rate is well above the 1-2% of the k-mer rule.

Outputs: `/anvme/workspace/v138da13-rustle/r1648_oarfish/`.

## RESULT r1648 (2026-11-12, job 4466969; `r1648_oarfish/{kb,hg}.result.txt`)
Abundance input: KB3781 5,107,217 primary reads of the regions (4,961,285 mapped), HG002 632,562. Walls 506 s and 280 s.
Scorer as registered (locus overlap; primary by record overlap).

| | KB3781 | HG002 |
|---|---|---|
| Bar 1 reach (>= 0.30) | 588 of 1,570 (0.375) pass | 572 of 2,412 (0.237) FAIL |
| oarfish right / wrong | 523 / 65 | 521 / 51 |
| primary right / wrong, same molecules | 528 / 60 | 520 / 52 |
| Bar 2: tool right, primary wrong vs reverse | 24 vs 29 FAIL | 36 vs 35 FAIL |
| Bar 3 negatives assigned (<= 0.10) | 0.311 FAIL | 0.517 FAIL |
| Bar 4 wall | pass | pass |
| EM isolated: molecules the minimap2 control left tied | 168: assigns 20, 4 right / 16 wrong | 606: assigns 108, 95 right / 13 wrong |
| kallisto, same scorer, for reference | 919 / 18, 15 vs 2 | 678 / 6, 59 vs 3 |

**FAIL as registered on both samples**, as expected: at posterior >= 0.99 oarfish equals the primary, with an 11% and
9% wrong rate. Where abundance alone decides, it is wrong 16 of 20 on KB3781 and 13 of 108 on HG002: the r1412 finding
(EM overconfident) on real reads. Negatives carry the stated region bias and are not interpreted. Closed.
(First HG002 scoring returned 0 picks: the r1644 names file carries a `>R|` prefix; fixed in the name list, no rerun of the tool.)
