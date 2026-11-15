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
