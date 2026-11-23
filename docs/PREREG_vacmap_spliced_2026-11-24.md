# PREREG: VACmap on Iso-Seq reads, against the genome and against isoform sequences (row 1936, 2026-11-24)

## Current state
Run 2026-11-24 (job 4499072, 13 min): NO GAIN on every arm, the padded-read addendum included (job 4499085). Outcomes at the end.

User (2026-11-24): the advisor asked for VACmap as well; "we both know it does not work with spliced alignment but what if
we can force it, maybe an option in vacmap or what if we extracted the loci with all isoforms". Written before any VACmap
alignment of these reads exists. Known at this point (setup only): VACmap (github micahvista/VACmap, cloned 2026-11-24) has
no spliced option (`-mode H|L|S|R|asm`; penalties `-globalpenalty`, `-localpenalty`); it reports a read as a primary record
plus supplementary records, so an intron can only appear as a deletion or as a split.

## Setup
- Sample and truth: HG002 Kinnex on CHM13, the 42,770 selected reads and the repaired truth table of row 1785
  (`r1785_truth_repair/hg/{reads.fa,reads.tsv}`); reused, as in rows 1776-1791, 1844.
- Arm G ("forced"): `vacmap -ref chm13v2.0.fa -read reads.fa -mode L` (HiFi mode), other options default.
- Arm GR: the same with `-mode R` (fixed variation penalty, described as more sensitive to gene conversion).
- Arm T (the user's idea): reference = the spliced sequences of every transcript of `annot_related/chm13.related.gtf`
  that overlaps a candidate locus of the table on the same strand (gffread); `vacmap -mode L` of the reads against them.
- Control TM: minimap2 `-ax map-hifi` (no spliced model) against the same transcript sequences, so that the effect of
  the reference is not credited to the aligner.
- Pick of an arm = the candidate locus of the read's table row that holds the most aligned bases of the read's primary
  and supplementary records on the locus strand (arms T, TM: the transcript's exons in genome coordinates); none =
  "elsewhere", counted wrong when the read is decidable. Baseline = the table's `primary` pick (minimap2 `splice:hq`).

## Step 0: what the alignments look like (reported, no bar)
Per arm: reads aligned; records per read; share of reads whose records together cover >= 90% of the bases the minimap2
primary aligned; arm G / GR: share of the minimap2 introns that reappear as a deletion or a split with both ends within
10 bases; MAPQ distribution on the near ties.

## Step 1: placement (the decision)
Decidable near ties (stratum S2) of the repaired table; baseline primary 553 right / 47 wrong of 600.
- GAIN: arm wrong <= 31 (a third of the room recovered, the bar of rows 1777-1783) AND, on the decidable clear-lead reads
  (S1), arm wrong <= primary wrong + 3.
- NO GAIN: arm wrong >= 47.
- PARTIAL: between.
Reported beside: fixes and breaks against the primary; the 36 decidable TBC1D3 reads of row 1844 (primary 16 / 20).

## Predictions (mine)
G and GR: most reads come back as several records, short exons are lost; placement NO GAIN, more wrong than the primary
(the candidates of a near tie are equally good for any scorer of the reference sequence). T: NO GAIN, close to TM; worse
where the annotation lacks the read's isoform. I do not expect VACmap's MAPQ to separate right from wrong picks.

## Limits
One sample, one individual, the seen quarter of the reads; the truth is the catalogue-free test with its known limits.
The transcript reference depends on the CHM13 annotation (Liftoff / CAT models). A NO GAIN says VACmap as released does
not place these reads better; it says nothing about its use on DNA.

## Outcome (2026-11-24, job 4499072; workspace `r1936_vacmap_tbc1d3/{T,TM,G,GR}.score.txt`, `pieces.txt`)

| arm | decidable near ties, right / wrong (600) | fixes / breaks | clear-lead reads wrong (14,596) | verdict |
|---|---|---|---|---|
| minimap2 `splice:hq` primary | 553 / 47 | | 18 | baseline |
| G, VACmap `-mode L` on the genome | 352 / 248 | 27 / 228 | 1,760 | NO GAIN |
| GR, VACmap `-mode R` on the genome | 341 / 259 | 13 / 225 | 524 | NO GAIN |
| T, VACmap `-mode L` on isoform sequences | 396 / 204 | 18 / 175 | 499 | NO GAIN |
| TM, minimap2 `map-hifi` on isoform sequences | 497 / 103 | 13 / 69 | 684 | NO GAIN |

- Step 0: G aligns 42,735 of 42,770 reads, median 2 records per read, 93% cover >= 90% of the bases minimap2 aligned, and
  51.4% of minimap2's introns reappear as a deletion or split within 10 bases (GR 42.8%, median 3 records). T: 85.6% of
  reads cover >= 90% (TM 86.0%): the annotation lacks the isoform of about one read in seven.
- Mechanism (post hoc, `pieces.py`): on the 600 decidable near ties G puts records on two or more candidate loci for 244
  reads, GR for 554; over all decidable reads 5,413 and 8,656 of 18,297.
- MAPQ does not separate: GR gives MAPQ 60 to 240 of its 259 wrong near-tie picks (308 of 341 right); G >= 30 to 89 of
  248 wrong and 51 of 352 right.
- TBC1D3, the 36 reads of row 1844 (no bar): primary 16 / 20, G 29 / 7 (fixes 14, breaks 1), GR 2 / 34, T 23 / 13,
  TM 18 / 18. One mode above the primary and the other far below it on the same reads; not followed up, and the truth
  at this family is the weakest of the table (row 1791).
- Predictions: NO GAIN held for all arms. "T close to TM" was wrong: VACmap on isoforms is 100 reads worse than minimap2
  on the same isoforms. "MAPQ does not separate" held.
- Reading: the isoform reference removes the intron problem and still loses 56 reads against the spliced genome
  alignment with minimap2 itself (TM), so the reference costs more than it gives; VACmap adds its own loss by spreading
  a read over copies. Closes VACmap as released for these reads. Not tried: tuned penalties (`-globalpenalty`,
  `-localpenalty`); a pick rule other than "most aligned bases".

## Addendum, registered before running (2026-11-24): reads padded with introns
User, after the outcome above: "my idea was to extract the full locus that spans all isoforms and extract that including
introns, then align that to genome"; asked where the locus sequence comes from, chose: each read turned into a DNA-like
sequence, its exons with the reference intron sequence filled in between, aligned to the genome by VACmap. Told beforehand:
the filled-in introns come from the copy minimap2 chose and pull the read back to it.
- Arm P1: every decidable read of the table rebuilt along its minimap2 primary record (read bases at match, insertion and
  clip columns; CHM13 sequence at every `N`), `vacmap -mode L` against the genome. Reads whose padded length exceeds
  500 kb are left out and count as unaligned.
- Arm P2 (the control for the caveat): the same read rebuilt along its best-AS record at ANOTHER locus (not overlapping
  the primary's), for the reads that have one.
- Pick and bars as above (GAIN = wrong <= 31 of the 600 decidable near ties, with the guard). Read beside the bar: the
  share of near ties where P1 picks the primary's locus and where P2 picks the runner-up's locus. If both are high, the
  placement is decided by the padding, not by the read.
- Predictions (mine): P1 repeats the primary on >= 95% of near ties (about 553 / 47, NO GAIN); P2 picks the runner-up's
  locus on >= 80% of them, so its right count collapses except where the runner-up is the origin.

### Outcome of the addendum (job 4499085, 8 min; `r1936_vacmap_tbc1d3/{P1,P2}.score.txt`, `pad.py`)
- P1 (padded along the primary, 18,297 reads): near ties 546 right / 54 wrong (fixes 1, breaks 8), the primary's locus on
  591 of 600; clear-lead reads 26 wrong (primary 18). NO GAIN. TBC1D3 36 reads: 16 / 20, all 36 at the primary's locus.
- P2 (padded along the best record at another locus, 15,129 reads): near ties 84 / 516; it picks the runner-up's locus
  on 522 of 597 (87%) and the primary's on 45. Clear-lead reads: runner-up 8,179 of 14,428 (57%), primary 5,408.
- Both predictions held (>= 95% and >= 80%). Reading: the reference introns decide where a padded read goes; on a near
  tie the read's own bases overrule them for 9 reads (P1) to 75 reads (P2) of 600. The padding cannot add information
  about the read's origin. `-mode R` on padded reads not run.

## Addendum 2, registered before running (2026-11-24): padded reads in the gene-conversion mode
User: "lets try the padded reads in vacmaps gene conversion mode". Arms P1R, P2R = the padded reads of P1, P2
(`P1.fa`, `P2.fa`, unchanged) with `vacmap -mode R`. Pick, bars and the padding-source reading as in the addendum.
Predictions (mine): P1R picks the primary's locus on >= 90% of near ties and is NO GAIN (wrong >= 47); P2R picks the
runner-up's locus on >= 70%. This mode cut pieces out of unpadded reads (554 of 600 on two or more loci), so I expect
more records per read than in mode L and somewhat less obedience to the padding, not better placement.
