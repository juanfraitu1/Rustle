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

## Result of the ablation (2026-11-21)
Jobs 4491343 / 4491344, 18 arms exit 0; `r1864_lost_introns/lost.{wtc,es}.txt`.
Lost introns: WTC11 6,183 (in 2,084 missed true chains), mouse ES 2,334 (in 697). **99.6% / 99.7% of them have 2 or more
primary reads carrying the exact junction** (WTC11: 2 reads 1,479, 3-5 2,858, 6-20 1,623, over 20 199), 99.1% / 97.8% are
canonical, and under 1% have a Rustle intron within 10 bp. So the junctions are well supported; what is under the floor is
the CHAIN: no two reads carry the same full chain.
| arm | WTC11: lost introns recovered | `x.gtf` correct / wrong added | mouse ES: recovered | correct / wrong added |
|---|---|---|---|---|
| floor1 | 3,961 (64.1%) | +4,270 / +26,166 | 1,406 (60.2%) | +1,041 / +16,536 |
| frac0 | 405 (6.6%) | +897 / +2,663 | 428 (18.3%) | +1,293 / +22,110 |
| dc0 | 341 (5.5%) | +611 / +1,590 | 326 (14.0%) | +688 / +5,800 |
| ss0, jr0, snapoff, ism0, ret0 | 0-20 each | | 0-31 each | |
| no arm | 1,797 (29.1%) | | 468 (20.1%) | |
- "The floor explains it" (70% on both): NOT met (64% and 60%), though the floor is by far the largest step.
- "A step is the pointer" (20% of the N >= 2 introns on both): none (the isoform fraction reaches 18% on mouse ES, 7% on WTC11).
- **Lost introns with reads that no arm recovers: 29.1% and 20.0%, at or over the 20% bar on both: pass 1 itself.**
- Expectation: the floor part held; the fraction part was too optimistic (18% on ES, not the larger share).

## Follow-up, registered now before it is computed: the introns no arm recovers, by read
Set U = lost introns with N >= 2 that are in no arm's `x.gtf` (WTC11 1,790, mouse ES 466); control set C = lost introns that
floor1 recovers. For every primary read carrying an intron of U or C, flags: (a) the read has another junction that is not
GT-AG / GC-AG / AT-AC; (b) its junctions' motifs point to both strands; (c) MAPQ under 20; (d) it has a supplementary
alignment (SA tag); (e) it has an intron shorter than 30 bp or longer than 200 kb. Per intron: ALL-FLAGGED(x) = every one of
its reads carries flag x. Reported: the share of U and of C that is ALL-FLAGGED for each flag and for any flag; intron length
of U against C (median, share over 100 kb).
**Reading:** a flag EXPLAINS the class when at least 40% of U and at most a third of that share of C is ALL-FLAGGED for it, on
both libraries. An explaining flag names the pass-1 rule to look at (a registered build row). No explaining flag: the cause
is not in these read properties and needs a trace of pass 1 on a handful of loci.
