# Row 1864: which step loses real introns that have reads? (pre-registered 2026-11-21)

User (2026-11-21): derive rules for precision and recall from the WTC11 comparison. Row 1862's largest recall-side class:
true chains StringTie reports and Rustle lacks that carry an intron found in NO Rustle transcript (WTC11 1,543, mouse ES
533; on mouse ES 63% of them have StringTie coverage of 2 or more, so they are not all under the 2-read floor). Never
examined. Measurement by ablation; nothing is built. Workspace `r1864_lost_introns` (`arms.sbatch`, `lost.py`). Committed
before any arm ran. Development data: WTC11 replicate 2 and mouse ES. Replicates 1 and 3 stay unopened.

## Target set, fixed
LOST INTRON = an intron of a truth chain that StringTie OR IsoQuant reports correctly and Rustle's `final.gtf` lacks, where
the intron is in no multi-exon transcript of Rustle's `x.gtf` (today's defaults). Per lost intron, from the BAM: primary
records carrying exactly that junction (N), splice motif (GT-AG / GC-AG / AT-AC or other), and whether a Rustle intron lies
within 10 bp at both ends.

## Arms (one switch each against `base` = today's defaults, stage `assemble`)
| arm | switch | step |
|---|---|---|
| floor1 | `RUSTLE_MIN_CHAIN_READS=1` | the 2-read floor |
| dc0 | `RUSTLE_DUP_COUNT=0` | read counting (row 1850) |
| ss0 | `RUSTLE_POLISH_SPLICED_SHARE=0` | spliced share (row 1848) |
| jr0 | `RUSTLE_POLISH_JUNCTION_RATIO=0` | junction replacement (row 1845) |
| snapoff | `RUSTLE_POLISH_JUNCTION_SNAP=off` | junction snap |
| frac0 | `--polish-isoform-fraction 0` (driver snapshot) | isoform fraction |
| ism0 | `RUSTLE_POLISH_ISM_RATIO=0` | sub-chain collapse |
| ret0 | `RUSTLE_POLISH_RETAINED_RATIO=0` | retained-intron rule |
An arm RECOVERS a lost intron when the intron is in a multi-exon transcript of the arm's `x.gtf`.

## Reported, per library
Lost introns by N (0-1, 2, 3-5, 6-20, > 20) and motif; per arm: lost introns recovered, in all and among those with N >= 2;
correct chains gained and wrong chains added by the arm in `x.gtf` (the price); lost introns no arm recovers, by N.

## Reading, fixed now
- **The floor explains it** when floor1 recovers at least 70% of the lost introns on both libraries: then the class is the
  one-read tier's (row 1808) and no new rule follows.
- **A step is the pointer** when its arm recovers at least 20% of the lost introns with N >= 2 on both libraries. The rule
  that follows (a narrower form of that step) is its own registered row, checked on replicates 1 and 3 and the twelve datasets.
- Lost introns with N >= 2 that no arm recovers (at least 20% on both) point at pass 1 itself (junction acceptance, strand,
  region handling) and are examined by read in a follow-up.

## Expectation
On WTC11 the floor explains most (two thirds of the class had coverage under 2). On mouse ES I expect the isoform fraction
under real read counts (frac0, dc0) to recover the larger share of the introns with reads: row 1849 showed the 8% loss of
correct chains under real counts comes from the fraction. I do not expect the junction steps to matter (they move introns
by at most 10 bp, which this class excludes only partly).

## Limits
- One switch at a time: an intron lost by two steps is recovered by neither arm.
- N counts exact junctions in primary records; Rustle's pass 1 has its own acceptance rules.
