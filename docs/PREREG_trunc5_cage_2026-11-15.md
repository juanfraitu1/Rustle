# PREREG row 1742 (2026-11-15): are the chains `trunc5` drops real transcript starts? CAGE support on the two human samples

User question: can an external yardstick say whether the opt-in `trunc5` stage (row 1677) is worth turning on. Written
and committed before any figure below was computed.

## Why this measure

RefSeq chain identity is already known for the dropped chains (0.03-0.08 correct; rows 1677, 1692) and cannot tell a
degraded read from a real unannotated internal start. gffcompare class codes and SQANTI3's structural category are fixed
by the construction of the class (a 5' truncation is contained in a longer model; row 1716 met the same limit). The
evidence that is not fixed by construction is a CAGE peak at the 5' end. Row 1716 found the own in-peak flag
(`sq_nested.py`, `Peaks.dist == 0`) to agree with SQANTI3's `within_CAGE_peak`, so the own flag is used and SQANTI3 is not
re-run.

## Inputs (existing files, nothing re-assembled)

- A119b: `r1677_chainrules/port/a.gtf` (default `PREFIX.gtf`) and `a.trunc5.gtf`.
- HG002 Kinnex on CHM13: `r1692_trunc5_hg/hg.gtf` and `hg.trunc5.gtf`.
- CAGE: refTSS v3.1 lifted to CHM13, `r1716_sqanti/refTSS.hs1.bed` (strand-aware).
- RefSeq on CHM13 for the annotated / unannotated split: `r1716_sqanti/ref.gtf`.

## Measure

Multi-exon transcripts of `PREFIX.gtf`. Dropped = transcript id absent from the `trunc5` GTF; kept = the rest. For each,
the emitted 5' end. Share in a CAGE peak on the transcript strand (primary) and within 50 bp of one (reported).
Background = the same share at the position 100 bp downstream in the transcript (`tss_and_shift`), per group. Each
group also split by RefSeq exact intron chain (reported). Script `bench/famcmp/trunc5_cage.py`.

Known bias, same for both groups: `PREFIX.gtf` emits the most extended 5' end, which sits in a peak less often than the
commonest read start (row 1718: annotated chains 0.29-0.42). Hence the gate.

## Gate

Kept RefSeq-correct chains: in-peak share >= 3 x their background on each sample. Otherwise the run is INVALID and no
bar is read.

## Bars (read on the dropped group, both samples)

- **SAFE:** in-peak share <= 2 x the dropped group's background on both samples. Reading: the dropped chains are not
  starts; `trunc5` is eligible as a drop layer of `PREFIX.final.gtf` (user decision; the row 1677 / 1692 precision
  bars stay FAIL as registered).
- **MIXED:** in-peak share > 0.10 on either sample. Reading: a real subset is dropped; the stage should tag, not drop,
  and exempt supported starts.
- **Otherwise:** stays opt-in as it is.

## Limits

Human only (no CAGE set for gorilla). Says nothing about 3' ends. refTSS is a multi-tissue atlas lifted from GRCh38:
a testis- or sample-specific start can be real and outside every peak. No default changes in this row.

## OUTCOME (2026-11-15, workspace `r1742_trunc5_cage`, `a.txt` / `hg.txt`): INVALID as registered

| sample | group | chains | in peak | within 50 bp | background | ratio |
|---|---|---|---|---|---|---|
| A119b | dropped | 36,762 | 0.048 | 0.095 | 0.016 | 3.1 |
| A119b | dropped, RefSeq-correct | 873 | 0.257 | 0.409 | 0.077 | 3.3 |
| A119b | kept | 182,755 | 0.272 | 0.466 | 0.094 | 2.9 |
| A119b | kept, RefSeq-correct (gate) | 37,842 | 0.344 | 0.648 | 0.189 | **1.83** |
| HG002 | dropped | 16,889 | 0.125 | 0.187 | 0.029 | 4.4 |
| HG002 | dropped, RefSeq-correct | 546 | 0.416 | 0.604 | 0.115 | 3.6 |
| HG002 | kept | 136,144 | 0.493 | 0.709 | 0.140 | 3.5 |
| HG002 | kept, RefSeq-correct (gate) | 35,199 | 0.463 | 0.796 | 0.220 | **2.10** |

**Gate FAILS on both samples (1.83 x and 2.10 x against 3 x), so no bar is read.** Cause, from the table: the background
position 100 bp downstream is in a peak for 0.19-0.22 of annotated chains (promoter peak clusters, and `PREFIX.gtf`'s most
extended 5' end lying upstream of the main peak); the background that worked for nested starts in row 1716 (0.016) is
not a background for full-length starts. The gate was mis-designed, not the samples.

Post hoc, not a verdict:
- The dropped chains are in a peak 5.7 x (A119b) and 3.9 x (HG002) less often than the kept ones, and 3.1 x / 4.4 x
  above their own background. Had the bars been read: not SAFE on either sample, MIXED on HG002 (0.125 > 0.10), neither
  on A119b (0.048).
- Excess over background: about 1,200 (A119b) and 1,600 (HG002) dropped chains with a CAGE-supported start, 3.5% and
  3.3% of all supported multi-exon starts, while the stage drops 16.7% and 11.0% of the multi-exon chains.
- Reading: mostly truncations, with a real minority that a CAGE-free rule cannot name. Consistent with tagging rather
  than dropping; not evidence for a default. A valid rerun needs a background that leaves the promoter (for example
  the same offset on the opposite strand, or >= 1 kb), registered as a new row.
