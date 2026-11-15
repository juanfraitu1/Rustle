# Row 1673 — which haplotype supports the real loci among those supported only by secondary alignments (HG002)

Registered 2026-11-13 before running, on the user's "lets proceed with the next step" after row 1672.
Assembly layer, descriptive; no read is assigned to a copy, nothing is dropped, no rule is fitted (HG002 is spent).
Not for `advisor_share` until the user says so.

**Correction to the row 1672 report.** I offered this check as "real loci vs a by-product of the one-haplotype
alignment". That was too strong. The truth of row 1672 matches a molecule against the locus's counterpart on BOTH
haplotype assemblies, so a molecule perfect only on the maternal assembly is classed L only when its maternal hit
sits on the locus's own maternal counterpart. That is a real origin (the other allele of the same locus), not a
by-product. The split below separates two kinds of real.

**Step 0 (no BAM pass)** from `r1672_shadowflag/hg/{mat.sam,pat.sam,loci.tsv}` and the same GTF, strict truth
(slack 0), `bench/famcmp/sec_only_hap.py`. For each of the 106 class-L loci: molecules whose strict near hits
overlap the locus's hits; each molecule's matching hits are paternal (the aligned haplotype), maternal, or both.

Per-locus label, frozen: with >= 3 matching molecules,
- **aligned-haplotype**: >= half have a paternal match (a single-reference pipeline could see it);
- **other-allele only**: otherwise (real for the individual; the reference copy is not the better match);
- fewer than 3 matching molecules: thin.

Stop rule: if thin > 35 of 106, add a per-molecule table to `sec_only_loci.py` and rerun the full job; else
step 0 is the result. Limit, disclosed: step 0 does not recompute a molecule's primary locus, so molecules sampled
for another locus that also match here are counted.

**Sex-chromosome count (descriptive, no bar).** `hg002.pat.fa` has chromosomes 1-22 and Y, no X; HG002 is male, so
X-linked molecules have no true locus on this reference. Count, per row-1672 class, loci whose judged molecules
have their strict best hit on maternal chrX (`CM039033.1`) and nowhere on the paternal assembly. This is the one
mechanism here that would be an alignment by-product, and a candidate reason HG002 flags 7.6% of loci against
1.3% on KB3781. The per-class count needs each locus's molecules; from the SAMs alone it is reported as: molecules
(all sampled) by class of best strict hit (paternal / maternal X / maternal autosome / both), and for the 106 L
loci the share of matching molecules whose match is on maternal X.

## Result (2026-11-13, `r1673_sechap/step0.out`, `r1672_shadowflag/hg/L_haplotype.tsv`)

Step 0 stands (thin 4 of 106). Of the 106 class-L loci: aligned-haplotype 82, other-allele only 20, thin 4.
Sampled molecules with a near-perfect hit, 4,330: both haplotypes 2,818, maternal autosome only 681, paternal only
528, maternal X only 303. Class-L loci with any matching molecule also near-perfect on maternal X: 32.

Not registered, found on reading the table (post hoc, descriptive): **all 82 aligned-haplotype loci, and 103 of the
106 class-L loci, sit on unplaced scaffolds** (`JAKCWS…`, 491 contigs, 60.7 Mb) of the paternal assembly, not on
its chromosomes. Of the 1,518 flagged loci 1,237 are on scaffolds. Row 1672's classes split as:

| | L | both | P | other | mixed | thin | total |
|---|---|---|---|---|---|---|---|
| chromosomes | 3 | 12 | 135 | 120 | 3 | 8 | 281 |
| scaffolds | 103 | 76 | 98 | 260 | 160 | 540 | 1,237 |

On chromosomes the picture is KB3781's: 3 real (all other-allele only), 135 shadows; the 120 "other" are
consistent with X-linked molecules (no X in the reference) but were not tied to X per locus. Row 1672's frozen
bars are not re-judged; on chromosomes alone bar 1 would hold (3, 0.011) and bar 2 would not (0.49).
Reading: the real loci hidden by the primary choice on HG002 are short single-exon loci on unplaced scaffolds of a
draft haplotype assembly, where the primary went to another copy of the same sequence. Not checked: what the
scaffolds are (rDNA, satellite, duplicated chromosome sequence), whether the molecules' primaries sit on a
near-identical chromosomal copy; molecule primaries were not recomputed.

## Follow-up: what the 106 are (2026-11-13, user "can we check for that?"; descriptive, post hoc; `r1673_sechap/{content.py,content.out,L_content.tsv}`)

No BAM pass; from the row 1672 SAMs and the locus sequences placed on CHM13 (minimap2 asm20, best hit).

- Other copies of the locus sequence on the paternal assembly: on scaffolds only 83; a chromosomal copy 23
  (identical 12, < 1% different 11); none unique.
- Matching molecules' best paternal hit away from the locus (1,121 molecule-locus pairs): another scaffold, 1-2
  mismatches worse, 1,071; a chromosome, 1-2 worse, 41; further away 9.
- CHM13 placement of the locus sequence: chr13 44, chr21 18, chr15 14, chr14 12, chr22 4, other 3, no hit 11.
  92 of 95 placed loci are on the five acrocentric chromosomes. Content: rRNA gene 30, Alu 43, non-repeat 20
  (8 on TEC genes), satellite 1, other 2.

Reading: the real loci are copies of the rDNA array and its surrounding acrocentric short-arm sequence, left as
unplaced scaffolds in the HG002 paternal assembly. Each molecule matches one scaffold copy exactly and the next
copy at 1-2 mismatches; which copy gets the primary is the aligner's tie-break. They are real by the haplotype
truth and are a tandem-array class, not gene-family members. Not checked: the Alu-classed loci were not shown to
lie inside rDNA spacer specifically; molecule primaries were not read from the BAM.

## Follow-up 2 (2026-11-13, user "can we check for that?"; descriptive; `r1673_sechap/{alu_rdna.out,primaries.py,primaries.tsv}`, job 4470533)

- Distance of each CHM13 placement to the nearest rRNA gene body (RepeatMasker LSU/SSU-rRNA_Hsa >= 500 bp):
  within 20 kb for all 30 rRNA loci, all 44 Alu-dominated loci and 18 of 21 others; 3 others > 200 kb away
  (not inspected). 92 of 95 placed loci are inside an rDNA array.
- Primaries read from the BAM for the 1,121 matching molecule-locus pairs: 577 have their primary on another
  scaffold and no BAM record at this locus (molecules sampled for a neighbouring locus; the disclosed limit of
  step 0); 224 have no primary at any re-alignment hit (not found, not chased). Of the 303 pairs with a BAM
  record at the locus and a primary found: exact tie (equal NM, equal AS) 249 (247 scaffold, 2 chromosome);
  record here worse in the BAM (more mismatches, lower AS) 42; fewer mismatches but lower AS 9; equal NM, lower
  AS 3. Two further pairs have the primary on the locus itself.
  So where both records exist, 82% are exact ties broken by the aligner between array copies.

## Follow-up 3 (2026-11-13, user "ok, lets check those"; descriptive; `r1673_sechap/unfound.primaries.tsv`, job 4470538)

- The 42 pairs with a worse BAM record at the locus: mostly maternal-allele matches (the molecule is exact on the
  locus's maternal counterpart; on the paternal reference the primary's copy is the closer one). 16 come from one
  locus, paternal chr1:120,526,678-120,532,817 (all 16 matching molecules maternal only; e.g. NM 2 at the primary,
  10 here).
- The 3 loci outside the rDNA arrays: paternal chr1:120.53 Mb = CHM13 chr1:120.28 Mb, transcribed unprocessed
  pseudogene AC244669.1 in the 1q21 duplicated region, 30 molecules, maternal allele only; paternal
  chr19:38.53 Mb = CHM13 chr19:38.83 Mb, centromeric-type satellite, no gene, maternal allele only; scaffold
  JAKCWS010000191.1 = CHM13 chrY:9.16 Mb, FAM197Y5 tandem-array pseudogene, thin.
- The 224 pairs without a found primary are 23 molecules; all 23 primaries are on other unplaced scaffolds
  (MAPQ 0-1, NM 0-2 for 22 of them), i.e. further array copies that the re-alignment's 10 secondary hits did not list.
Conclusion: of 106 real loci, 1 (AC244669.1) is a gene-family-type case, and it is visible only through the
maternal allele; the rest are rDNA-array or other tandem/satellite copies.

## Follow-up 4 (2026-11-13, user "lets check those, also ... soft masked?"; descriptive, post hoc; `r1673_sechap/{softmask.py,softmask.out}`, `*/softmask.tsv`)

- The 11 loci without a CHM13 hit under asm20: 10 place exactly (map-hifi, 100% identity) on chr13/14/15 short
  arms, the 11th is a 163 bp 5.8S rRNA fragment identical in sequence to two of them. All 11 are rDNA; 103 of 106.
- AC244669.1 (paternal chr1:120,526,678-120,532,817) is a truth artifact, not a real hidden locus. The locus
  sequence's nearest maternal sequence (mat 120.776 Mb) is 14 mismatches away over 6.1 kb, and equally close to a
  different paternal copy (pat 120.899 Mb). The 16 molecules fit mat 120.76-120.78 Mb at 1-2 mismatches and
  pat 120.86-120.90 Mb at 2 (where their primaries are), against 4 at the locus. They come from the copy the
  primary chose; the truth called the locus real because it takes a locus's best hit on the other haplotype as its
  counterpart however distant. So 0 of the 106 is a gene-family locus hidden by the primary choice.
  Limit exposed: `copy_loci` has no identity floor on the cross-haplotype counterpart; "other-allele only" labels
  (20 loci) can contain such cases.
- Soft-masking (both references are soft-masked; fraction of lowercase bases over the locus's exons):

  | sample, where | class | n | median | share >= 0.5 |
  |---|---|---|---|---|
  | HG002 chromosomes | shadow (P) | 135 | 0.10 | 0.02 |
  | | elsewhere (other) | 120 | 1.00 | 0.75 |
  | | real (L) / both | 3 / 12 | 0.30 / 0.14 | 0.33 / 0.00 |
  | | not flagged | 17,676 | 0.18 | 0.07 |
  | HG002 scaffolds | L / both / P / other | 103 / 76 / 98 / 260 | 0.38 / 0.00 / 0.34 / 0.13 | 0.24 / 0.00 / 0.47 / 0.30 |
  | KB3781 | shadow (P) | 155 | 0.10 | 0.04 |
  | | elsewhere (other) | 10 | 0.61 | 0.60 |
  | | not flagged | 14,666 | 0.13 | 0.02 |

  Soft-masking does not separate real from shadow (shadows are unmasked paralog sequence; rDNA is largely
  unmasked). It does mark the "elsewhere" class: repeat-sequence loci whose molecules originate at neither the
  locus nor the primary's locus. Not registered as a rule; HG002 and KB3781 are both seen.

## Follow-up 5 (2026-11-13, user "lets checl for those"; descriptive, post hoc; `r1673_sechap/{floor.out,allele.out,reciprocal.py,reciprocal.out}`)

- Identity floor (0.3%, the truth's own near-perfect criterion) on the locus's maternal counterpart: 104 of 106
  counterparts are near-perfect; only the chr19 satellite locus loses its matching molecules. AC244669.1 passes
  the floor (14 / 6,138 = 0.23%), so a floor does not catch it; its problem is that the same maternal sequence is
  equally close to another paternal copy.
- The 20 "other-allele only" loci, by whether the locus is the best paternal fit of their maternal-only
  molecules: another paternal locus fits better 1 (AC244669.1, 16 of 16 molecules); no paternal hit at the locus
  in the re-alignment 19 (18 rDNA + the chr19 satellite; arrays have more copies than the 10 hits listed, so the
  test cannot decide them, and "allele vs copy" is not well defined inside an array).
- Earlier rows with stored haplotype truths (`hits.tsv`), re-read with a reciprocal filter: a copy's locus on a
  haplotype is dropped when a different catalog copy fits the same place better (mismatch rate lower by > 0.002
  over >= half of it). The filter is NOT validated and may drop true loci; read the figures as a sensitivity
  bound, not a correction. Assigned rows counted right, before -> after:
  o2mstie/truth2_ms10 33 -> 33 of 35; o2mstie/truth_all10 109 -> 109 of 259; r1653 hgpat 146 -> 146 of 176;
  r1657 hgpat 92 -> 82 of 104; r1644 420 -> 420 of 476; r1651 hg13 165 -> 162 of 231; r1653 hg13 293 of 390
  (from 297); r1657 hg13 28 -> 28 of 54. Scored rows whose truth becomes "outside the catalog" under the filter:
  20/70, 3,488/8,763, 18/1,714, 0/796, 0/4,189, 181/4,707, 32/1,607, 55/522.
  The O2 truth already sets aside reads whose truth spans two catalog loci, which is why "right" barely moves;
  the exposed mode is a place that matches no catalog copy well being credited to the nearest copy.
  Only r1657 on HG002's own haplotype moves materially (92 -> 82); no earlier verdict is re-judged here.

## Follow-up 6 (2026-11-13, user "lets checl for that"; descriptive, post hoc; `r1673_sechap/{gated_filter.py,t1657/}`; row 1657 is another session's, read only)

Correction to follow-up 5: the "92 -> 82" figure for row 1657 was computed on `truth/per_read.tsv`, which uses
`haplotype_truth`'s loose locus rule. Row 1657's quoted figures come from `bench/o2_contested/strict_truth.py`,
which since r1510 already requires a copy's best hit on a haplotype to be near-perfect (the weakness was known).

Row 1657 arms re-scored with `strict_truth.py`, before -> after removing copy places that fail the reciprocal
test under the gated rule (5 of 665 on the own haplotype, 1 of 625 on CHM13); "before" reproduces the stored
figures exactly:
- own haplotype: flag-arm base 160 / 26 -> 161 / 26; union 170 / 14 -> 171 / 14;
- CHM13: base 81 / 75 -> 81 / 75; union 108 / 105 -> 108 / 105.
Row 1657's verdict is unchanged.

The six dropped places checked by hand against positions: 5 right (MCL171 copies 0 and 1 claiming the maternal
places of copies 4 and 2, which match those at NM 0-1 and are syntenic; MCL528 copy 0 on copy 1's maternal place;
MCL1313 chr19 copy on a chr9 place owned at NM 0 by a chr9 copy), 1 ambiguous (MCL442, 1q21: copy 0's syntenic
maternal place matches copy 1 exactly and copy 0 at 4 mismatches; sequence and position disagree, consistent
with exchange between the copies). Six cases is a small check, not a validation of the filter in general.
