# Row 1857: polyA tails on a second library, and a 3' end placed at the tailed reads (pre-registered 2026-11-21)

User (2026-11-21): get a second untrimmed library, even a small one, and try it out. Row 1856 (mouse ES, the only library
with tails among the twelve) found no filter rule in polyA tails but good evidence for the END position (3' end annotated
within 50 bp: 0.80 of the chains with a majority of tailed reads, 0.37 of those with none). Workspace `r1857_polya_ends`.
Committed before the new library was downloaded or looked at, and before any end was re-placed on either library.

## Second library
ENCODE ENCSR507JOF (LRGASP, human WTC11 iPS cells, PacBio cDNA; the sister experiment of the mouse ES library ENCSR172GXL),
replicate 2, file ENCFF370NFS (1.08 GB, the smallest of the three). Mapped untrimmed to CHM13 v2.0 with the command of
`genomes/map_reads.sbatch`; truth = RefSeq on CHM13 (`ptr_ppa_rules/gffcmp/hsa.ref.gtf`, as for A119b). Assembled with
today's defaults. Whether it carries tails is checked first; if under 10% of its reads do, another LRGASP file is taken and
that is reported.

## Part A: does row 1856's picture repeat? (predictions, written now)
With `tails.py`'s definitions on the new library: (A1) at least 75% of the multi-exon chains have a majority of tailed reads,
correct or not (ES: 85%); (A2) NOTAIL (>= 3 reads, none tailed) holds under 5% of the wrong chains (ES: 1.6%); (A3) among
chains with >= 5 reads the 3' end is annotated within 50 bp at least 0.20 more often for majority-tailed chains than for
chains with no tailed read (ES: 0.80 against 0.37). REPEATS = all three.

## Part B: the end rule, measured as a re-placement (no assembler code)
For a multi-exon transcript of `final.gtf` with at least 2 tailed reads of its exact chain: NEW 3' end = the most frequent
aligned end position of its tailed reads (ties: the most downstream). Compared with the transcript's 3' end as emitted
today (no end layer, row 1846), over the same transcripts:
- share with the end within 50 bp, and within 10 bp, of an annotated transcript end on the strand (any transcript), for all
  such transcripts and for the correct chains only (end of a transcript with that chain);
- median distance moved, and the share moved by more than 50 bp.
**HELPS** on a library when the share within 50 bp rises by at least 1 point (row 1831's bar for the 5' layer) and the share
within 10 bp does not fall. **Candidate** = helps on mouse ES (development) and on WTC11 (check). A candidate becomes a build
row: an opt-in end layer that switches on when the read check finds tails.

## Expectation
Part A repeats (same protocol, same consortium). Part B: today's end is the transcript's extreme or modal read end over all
reads, tailed or not; the tailed mode should agree with it for most transcripts and correct the ones dragged by untailed
reads. I expect a gain of a few points at 50 bp and a larger one at 10 bp.

## Limits
- Two libraries of one protocol and one sequencing centre; nothing here says how other untrimmed libraries behave.
- RefSeq on CHM13 is a lifted annotation; annotated ends are weak truth (row 1831).
- A re-placement outside the assembler does not show interactions with the filter, which reads end features.
