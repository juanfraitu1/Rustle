# PREREG: do unspliced reads carry information for spliced transcripts? (row 1989, 2026-11-27)

The user's question: StringTie keeps single-exon reads in the same graph as the spliced ones, so the two share information;
can Rustle do something similar, and would it raise the metrics?

## Facts first
- StringTie (`rlink.cpp`, clone of 2026-11-25): an unspliced read is a transfrag over one or more contiguous nodes of the
  same splice graph. It adds coverage to those nodes, so to every transcript through them; contiguous transfrags get LESS
  weight than spliced ones (`nosplice`, line 6964); single-exon predictions inside a multi-exon one, or in its introns, are
  removed at the end by coverage (lines 19518-19625). It is shared coverage, not shared structure.
- Row 1893: with `-L` StringTie's graph builds almost nothing beyond read chains (95.5-99.97% of its output are read chains).
- Rustle counts a transcript's reads by intron chain. An unspliced read inside an exon supports no spliced transcript. Unspliced
  reads make single-exon models, which the sub-chain step and the shadow rule then judge against spliced ones.
- Closed, not repeated: single-exon floors and shadow ratios (rows 1832, 1847), single-exon rescue by family (1919, 1920),
  `RUSTLE_SEED_UNSPLICED=support` for secondaries (1520).

## What is tested: unspliced reads as support for a one-read spliced locus
Row 1987 left 7,000 one-read loci per human library outside the default because no other read shares their introns. An
unspliced read wholly inside one of the chain's exons is the evidence StringTie would count.
`bench/o2_truth/unspliced_support.py`: per transcript, unspliced primary reads lying wholly inside one of its exons (4 bp
slack), end exons and internal exons apart. Also the size of the resource: how many unspliced primary reads lie inside an
exon of an assembled spliced transcript, across one of its introns, or elsewhere.
WTC11 r1 explored, r2 judged once (tables of row 1987: all introns seen again in the other replicates, all annotated).
- **N1:** among the one-read loci that are NOT supported by another spliced read, those with >= 2 unspliced reads inside
  their exons number >= 500 on r2 and have all introns seen again at >= 0.60 (the class as a whole: 0.33-0.48).
- **N2:** the lift exists at all: >= 2 unspliced reads against none differ by >= 0.10 in all-introns-seen-again on r2.
Prediction: N2 holds weakly, N1 fails (unspliced reads pile up in noise and in intronic pre-mRNA as well as in exons).

## Result (job 4515521): no usable information for spliced loci
Size of the resource, WTC11 r1 (2.56 million molecules): 251,383 unspliced primary reads (10%); 86,262 inside an exon of
an assembled spliced transcript, 64,922 inside a transcript's span but not inside an exon, 100,199 elsewhere. r2: 77,391 /
39,819 / 89,923.
One-read loci that no other spliced read supports, by unspliced reads inside their exons (loci; all introns seen again;
all introns annotated):
| unspliced reads | r1 (explored) | r2 (judged) |
|---|---|---|
| 0 | 6,357; 0.275; 0.238 | 6,192; 0.368; 0.276 |
| 1 | 446; 0.339; 0.233 | 599; 0.546; 0.316 |
| 2-4 | 181; 0.243; 0.193 | 311; 0.521; 0.318 |
| 5-19 | 65; 0.200; 0.092 | 121; 0.579; 0.438 |
- **N1 FAILS:** on r2 about 440 loci have >= 2 such reads (bar 500) and they are seen again at 0.52-0.58 (bar 0.60).
- **N2 holds on r2 as written (+0.16) and is CONTRADICTED on r1** (0.24 and 0.20 against 0.275 with none): not a stable signal.
Where spliced support exists it is already decisive and unspliced reads add nothing: supported one-read loci 0.81 with
none, 0.77 with one; two-read chains of the main assembly 0.893 / 0.904 / 0.889 with 0 / 1-4 / >= 5 unspliced reads.
Reading: sharing coverage between unspliced and spliced reads, as StringTie's graph does, would not raise these metrics on
PacBio cDNA; nine in ten molecules are spliced and the unspliced ones do not mark real loci. Not built.
Not tested: unspliced reads for transcript ENDS or abundance, and direct-RNA / ONT libraries where far more reads are unspliced fragments.

## Addendum before running (2026-11-27, user: test the uses left open)
Data: WTC11 on nanopore, direct RNA and cDNA (LRGASP, mapped to CHM13 by row 1976: `r1976_read_types_real/ont_drna/x.bam`,
`ont_cdna_sub/x.bam`; 31% and 22% of their primary reads are unspliced against 10% on PacBio), with the PacBio libraries
of the same cell line as an independent judge; simulated fly / mouse / zebrafish for counts (the read name is the transcript).
**A. Libraries with many unspliced reads.** One-read loci on the two nanopore libraries (`RUSTLE_ONE_READ_LOCI=1` logic on
a primaries-only assembly), judged by "every intron is in a PacBio read of WTC11 r1-r3" and "all introns annotated".
- N3: among the loci no other spliced read supports, on direct RNA, those with >= 2 unspliced reads inside their exons
  have all introns in PacBio at >= 0.10 above those with none.
- N4 (the row 1987 default on nanopore, reported with it): supported loci have all introns in PacBio at >= 0.80.
**B. Ends** (`bench/o2_truth/unspliced_ends.py`). Assembled spliced transcripts whose chain is annotated, with >= 3
unspliced reads lying in the last exon (not crossing its splice site; they may run past the assembled end). The 3' end as
assembled against the median 3' end of those reads, each against the nearest annotated end of that chain. PacBio r1 and
direct RNA.
- E1: the unspliced reads' end is within 50 bp of the annotated end at least 0.05 more often than the assembled end.
**C. Abundance** (`bench/o2_truth/unspliced_counts.py`). Simulated sets: true count of a transcript = its simulated reads.
The assembler's `reads` of the exactly matching transcript, against `reads` + the unspliced reads inside its exons (a read
inside exons of several transcripts is shared in proportion to their `reads`).
- A1: adding them lowers the median relative error of the count by >= 0.02 on at least 2 of 3 sets.
Reported beside it for the real libraries: the share of primary reads that are unspliced and inside an exon of a spliced
transcript (the most a count could gain).
Predictions: N3 fails, N4 holds, E1 fails (the assembled end already rests on the spliced reads), A1 fails on simulated
reads (few are unspliced) while the real nanopore share is large enough to matter.

## Result of the addendum (jobs 4515540, 4515541)
**A. Nanopore, judged by the PacBio reads of the same cell line** (loci; all introns in a PacBio read; all introns annotated):
| class | direct RNA | cDNA (subsample) |
|---|---|---|
| one-read loci, every intron in another spliced read (the row 1987 default) | 3,902; 0.972; 0.905 | 3,841; 0.955; 0.869 |
| other one-read loci, no unspliced read inside the exons | 5,420; 0.602; 0.418 | 10,258; 0.306; 0.229 |
| ... 1 | 901; 0.710; 0.488 | 1,669; 0.332; 0.238 |
| ... 2-4 | 813; 0.750; 0.582 | 1,738; 0.324; 0.238 |
| ... >= 5 | 383; 0.695; 0.525 | 2,109; 0.175; 0.132 |
| main assembly, two-read chains | 6,718; 0.932; 0.754 | 7,748; 0.791; 0.570 |
- **N3 HOLDS on direct RNA** (>= 2 unspliced reads: 1,196 loci at 0.732 against 0.602, +0.13) and is ABSENT on cDNA (0.24
  against 0.31) as on PacBio r1. A direct-RNA effect only, and the class stays far under the two-read chains (0.93).
- **N4 HOLDS:** the supported one-read loci of row 1987 carry over to nanopore unchanged, above the two-read chains.
Unspliced primary reads inside an exon of an assembled spliced transcript: direct RNA 136,334 of 745,729 primary reads (18%).
**B. Ends: E1 FAILS, the unspliced reads are worse.** Within 50 bp of the nearest annotated 3' end, assembled / median of
the unspliced reads in the last exon: PacBio r1 0.669 / 0.633 (10,969 transcripts), direct RNA 0.677 / 0.600 (5,125).
Differing by > 10 bp, better / worse: 2,057 / 2,510 and 419 / 1,199.
**C. Counts: A1 FAILS as registered, direction right.** Simulated counts are already exact (median relative error 0.000 /
0.003 / 0.000); 342 / 351 / 435 transcripts gain >= 1 read and 83 / 81 / 88% of them move closer to the truth; total
counted / true 1.002 -> 1.013 (fly, overshoot), 0.974 -> 0.975, 0.928 -> 0.935. No count truth for real nanopore reads.
Nothing built. The one use that looks worth a build is a COUNT that includes unspliced reads on direct RNA (18% of the
reads); it needs a truth first (spike-ins or a matched quantification).

## Second addendum before running (2026-11-27, user: use the SIRVs, or Anaquin, as the count truth)
- Anaquin analyses sequins, synthetic RNAs that must be IN the tube; it cannot give a truth to a library made without
  them, and none of ours has sequins. A public set with sequins and SIRVs on direct RNA exists (SG-NEx); not downloaded.
- SIRV-Set 4 is what the LRGASP libraries carry (the mouse ES BAM has SIRV1-7, long SIRVs, ERCC contigs): 69 isoforms of
  7 genes ALL AT THE SAME MOLARITY (mix E0), 15 unspliced long SIRVs equimolar with them, 92 ERCC over six orders.
  So the truth is "every SIRV isoform equally abundant". The WTC11 nanopore libraries are LRGASP too; row 1976 mapped
  them to CHM13 only, so step 0 = map the raw reads to the spike-in contigs and see whether SIRV reads are there.
Measure (`bench/o2_truth/sirv_counts.py`), models = the TRUE spliced SIRV isoforms (no assembly in the way): count by exact
intron chain, against that count + the unspliced primary reads inside the isoform's exons (shared in proportion to the
chain counts). Libraries: mouse ES PacBio, WTC11 direct RNA, cDNA, PacBio.
- **S1 (isoforms):** on direct RNA the median |log2(count / median count)| over the spliced isoforms falls by >= 0.05.
- **S2 (genes):** on direct RNA the coefficient of variation over the 7 genes of (gene count / number of isoforms) falls.
Prediction: S1 fails (sharing in proportion keeps the ratios inside a gene), S2 unclear.

## Result of the second addendum (job 4515577): the SIRVs are in every WTC11 library; mixed answer
61 spliced SIRV isoforms in 7 genes, all equimolar. Count by exact chain -> with the unspliced reads inside the exons:
| library | reads on the SIRV genes | exact chain | spliced, other chain | unspliced inside an exon | median abs log2(count / median) | CV over genes of count per isoform | isoforms with no read |
|---|---|---|---|---|---|---|---|
| mouse ES PacBio | 44,502 | 37,900 | 3,927 | 1,279 | 0.771 -> 0.796 | 0.421 -> 0.383 | 0 -> 0 |
| WTC11 PacBio | 13,839 | 11,910 | 1,097 | 346 | 1.056 -> 1.041 | 0.487 -> 0.462 | 1 -> 1 |
| WTC11 nanopore cDNA | 117,980 | 65,023 | 37,428 | 11,754 | 1.349 -> 1.148 | 0.373 -> 0.295 | 4 -> 2 |
| WTC11 nanopore direct RNA | 4,000 | 1,752 | 1,502 | 591 | 0.886 -> 0.615 | 0.203 -> 0.287 | 5 -> 5 |
- **S1 HOLDS** on direct RNA (-0.27; my prediction was wrong) and on cDNA (-0.20); nothing on PacBio (+0.03, -0.02).
- **S2 FAILS** on direct RNA (the spread between genes RISES, 0.20 -> 0.29); it falls on the three other libraries.
- Direct RNA rests on 1,752 exact reads over 61 isoforms, five of them with none: thin.
- The larger loss on nanopore is not the unspliced reads: 38% (direct RNA) and 32% (cDNA) of the reads on the SIRV genes
  are spliced with a chain that is no isoform's chain (cut short or with a wrong junction), against 8-9% on PacBio.
Not built. A count for nanopore would have to place partial chains first; unspliced reads are the smaller part.
