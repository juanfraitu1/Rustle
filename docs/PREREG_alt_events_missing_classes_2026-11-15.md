# PREREG row 1744 (2026-11-15): the event classes row 1714 did not measure, and whether alignment loses short exons

User question: are all alternative splicing events modelled as the annotation expects, and would a graph (vg) help with
the tricky ones. Row 1714 covered shifted exon ends, alternative first / last exons, nested starts / ends and a cassette
control. Not measured: retained introns, mutually exclusive exons, microexons. Measurement only; written and committed
before the job ran. No default changes.

## Classes added to `bench/famcmp/alt_events.py`

Pairs of annotated RefSeq intron chains, as in row 1714 (both expressed = >= 2 full-chain primary reads each; read at
minor >= 5 reads; share with both chains emitted):

- **RI**: one chain equals the other without one internal intron (a terminal one is the nested class of row 1714).
- **MXE**: the chains differ in one internal exon replaced by a disjoint exon between the same outer splice sites.
- **MICRO**: cassette pairs whose skipped exon is <= 30 bp (a subset of CASSETTE).

**Skip reads** (the alignment question). The pair measure cannot see an exon the aligner misses: reads that skip it do
not spell the chain, so the pair is never "expressed". Hence a second table: for each annotated chain and each internal
exon that no annotated chain skips, the full-chain reads and the reads spelling the same chain with that exon skipped,
units with >= 5 such reads, by exon length (<= 15, 16-30, 31-60, > 60 bp). Unit = chain x exon (an exon shared by
isoforms counts once per chain). A skip read is a missed alignment or a real unannotated skip; the > 60 bp bin is the
control for the latter.

## Inputs

Row 1714's: OR6737 (`r1628/orrun/or.gtf`, `.final.gtf`), A119b (`r1617/default/a.gtf`, `.final.gtf`), with the benchmark
StringTie and FLAIR GTFs; KB3781 (`r1622/kb/kb.gtf`) reported without tools. Workspace `r1744_altevents2`.

## Reading rules (minor >= 5 reads; OR6737 and A119b decide, KB3781 reported)

1. **Per pair class (RI, MXE, MICRO):** fewer than 30 expressed pairs on a sample = too few to read there. Otherwise
   both-emitted share within 0.05 of that sample's CASSETTE share = modelled like ordinary splicing; lower = a gap,
   reported with the split major-only / minor-only / neither.
2. **Alignment (skip reads):** an alignment problem a graph could address is flagged when the skip-read share of the
   <= 30 bp exons (both short bins pooled) exceeds the > 60 bp control by >= 0.05 absolute on both deciding samples.
   Flagged = a graph test is worth registering, with that share as the number to beat. Not flagged = no alignment loss
   on short exons for a graph to recover; vg is not pursued for splicing events on this evidence.

## Limits

RefSeq pairs only; an RI chain's reads can be unspliced precursor RNA; a skipped short exon can also be absorbed as
an insertion or mismatches beside the junction, which gives the same skip chain and is counted; exons skipped by reads
that also differ elsewhere are not counted. The assemblies are the row 1714 files (default assembly rules have not
changed since).

## OUTCOME (2026-11-15, job 4472241, `r1744_altevents2/{or,a,kb}.txt`)

Pairs with both chains emitted, minor >= 5 reads (expressed pairs in brackets):

| class | OR6737 default | StringTie | FLAIR | A119b default | StringTie | FLAIR | KB3781 default |
|---|---|---|---|---|---|---|---|
| CASSETTE (control) | 0.975 (3,531) | 0.799 | 0.916 | 0.912 (7,382) | 0.695 | 0.876 | 0.942 (5,045) |
| MICRO (skipped exon <= 30 bp) | 0.985 (131) | 0.809 | 0.908 | 0.917 (218) | 0.679 | 0.913 | 0.952 (167) |
| MXE | 0.953 (232) | 0.707 | 0.901 | 0.879 (471) | 0.573 | 0.919 | 0.858 (379) |
| RI | 0.744 (43) | 0.698 | 0.953 | 0.682 (192) | 0.557 | 0.906 | 0.584 (89) |

Rule 1: MICRO and MXE are within 0.05 of the control on both deciding samples: modelled like ordinary splicing (MXE on
KB3781, reported only, is 0.084 below its control). **RI is a gap**: 0.231 and 0.230 below the control; the lost form is
almost always the minor one (major only 10 of 43 and 52 of 192); FLAIR is ahead on this class (0.953 / 0.906), StringTie
behind. `final.gtf` equals the default on all four classes.

Skip reads (chain x exon units with >= 5 reads; skip share of reads; units where skip reads >= full-chain reads):

| exon length | OR6737 units | skip share | skip >= full | A119b units | skip share | skip >= full | KB3781 units | skip share | skip >= full |
|---|---|---|---|---|---|---|---|---|---|
| <= 15 bp | 44 | 0.109 | 16 | 160 | 0.627 | 94 | 60 | 0.055 | 18 |
| 16-30 bp | 932 | 0.026 | 29 | 1,425 | 0.024 | 66 | 1,219 | 0.009 | 31 |
| 31-60 bp | 12,886 | 0.005 | 114 | 18,862 | 0.008 | 267 | 15,525 | 0.001 | 125 |
| > 60 bp (control) | 163,204 | 0.004 | 689 | 228,967 | 0.004 | 1,789 | 185,435 | 0.002 | 858 |

Rule 2, <= 30 bp pooled against the control: OR6737 0.040 vs 0.004 (+0.036), A119b 0.070 vs 0.004 (+0.066), KB3781
0.036 vs 0.002. **NOT flagged as registered** (OR6737 is under the +0.05 bar).

Post hoc, not a verdict: the pooled bar hides a concentrated loss. Exons of <= 15 bp are skipped by as many reads as
spell them in 36% (OR6737), 59% (A119b) and 30% (KB3781) of units, against 0.4-0.8% for exons over 60 bp; the default
emits only the skip chain for 13 of 44, 79 of 160 and 18 of 60 of them, and StringTie and FLAIR do the same (same
alignments). It is a small class (44 / 160 / 60 chain x exon units). 16-30 bp exons are skipped 5-7 x more often than
the control but still rarely (2.5-4.6% of units). Not shown: that the skip reads are misalignments rather than real
skipping (the length gradient and the three-sample agreement suggest it; no read was realigned).
