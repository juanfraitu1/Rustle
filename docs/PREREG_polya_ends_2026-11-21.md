# Row 1860: polyA tails on a second library, and a 3' end placed at the tailed reads (pre-registered 2026-11-21)

User (2026-11-21): get a second untrimmed library, even a small one, and try it out. Row 1856 (mouse ES, the only library
with tails among the twelve) found no filter rule in polyA tails but good evidence for the END position (3' end annotated
within 50 bp: 0.80 of the chains with a majority of tailed reads, 0.37 of those with none). Workspace `r1860_polya_ends`.
Committed before the new library was downloaded or looked at, and before any end was re-placed on either library. (Registered first as row 1857; renumbered at once: another session had claimed 1857-1859.)

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

## Result (2026-11-21): the picture repeats; the end rule helps on mouse ES and HURTS on the second library: not a candidate

Second library: ENCFF370NFS mapped and assembled (job in `r1860_polya_ends/map.jobid`, exit 0; `reads/HSA/WTC11_r2_mm.bam`,
1,425,275 primary spliced reads; assembly `r1860_polya_ends/wtc`, repo and binaries of `commit.txt`). The read check finds a
tail on 99.9% of its reads (mouse ES 60%). Today's defaults on it: 20,016 correct of 31,336 multi-exon chains (0.639, RefSeq
on CHM13). Tables `ends.es.txt`, `ends.wtc.txt` (`ends.py`).

**Part A (predictions): REPEATS, 3 of 3.**
| | mouse ES | WTC11 | prediction |
|---|---|---|---|
| A1 chains with a majority of tailed reads (correct / wrong) | 0.850 (0.899 / 0.815) | 0.994 (0.997 / 0.989) | >= 0.75: yes |
| A2 NOTAIL, share of the wrong chains | 0.016 | 0.002 | < 0.05: yes |
| A3 end annotated within 50 bp, majority tailed minus no tailed read | 0.802 - 0.365 | 0.620 - 0.000 (10 chains) | >= 0.20: yes, on 10 chains |
On WTC11 nearly every read has a tail, so a tail tells nothing apart there: 48 chains have no tailed read. A filter from
tails is even less possible than on mouse ES.

**Part B, 3' end at the tailed reads' most frequent end** (transcripts with >= 2 tailed reads; share with the end within 50 /
10 bp of an annotated end, emitted -> re-placed):
| | mouse ES (development) | WTC11 (check) |
|---|---|---|
| all transcripts, any annotated end | 34,962: 0.708 -> 0.773 (+0.065); 0.635 -> 0.721 (+0.086) | 31,090: 0.661 -> 0.612 (-0.048); 0.604 -> 0.576 (-0.028) |
| correct chains, ends of that chain's transcripts | 16,061: 0.562 -> 0.626 (+0.064); 0.463 -> 0.560 (+0.097) | 19,946: 0.688 -> 0.620 (-0.068); 0.619 -> 0.581 (-0.038) |
| moved by more than 50 bp | 28.5% | 18.2% |
**HELPS on mouse ES, fails on WTC11** (5 points worse at 50 bp, 3-4 at 10 bp): not a candidate.

- Expectation: wrong for the check library. The gain on mouse ES does not come from the tail as such: where every read is
  tailed (WTC11) the tailed mode is simply the most frequent read end, and the end emitted today matches the annotation
  better than that. On mouse ES the 40% untailed reads pull today's end away, and restricting to tailed reads undoes it.
- So the usable statement is narrower than a rule: on a library where a large share of reads LACK tails, untailed reads are
  truncated at the 3' side and should not vote on the end. That is one library here.
- The two libraries differ in species and annotation (GENCODE on GRCm39, RefSeq lifted to CHM13), so "worse against RefSeq"
  may be the annotation's choice of the distal end; not separable with these data.

**Decision:** no polyA rule is built. Tails neither filter chains (rows 1856, 1860 A) nor place ends reliably (1860 B). The
read check's warning stays. The WTC11 library is now available as a thirteenth dataset (human iPS cells, 0.639 at today's
defaults) and is unused for any other rule.
