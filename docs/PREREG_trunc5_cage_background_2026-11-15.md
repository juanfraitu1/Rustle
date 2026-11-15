# PREREG row 1743 (2026-11-15): row 1742 rerun with a background outside the promoter

Row 1742 was INVALID: its background (100 bp downstream in the transcript) is still inside promoter peak clusters for
annotated chains (0.19-0.22 in a peak). Same inputs, same measure, same bars; only the background changes. Written
and committed before the new background was computed.

## What is already seen (not blind)

Row 1742's table: in-peak share of the dropped chains 0.048 (A119b) and 0.125 (HG002), kept 0.272 / 0.493, kept
RefSeq-correct 0.344 / 0.463. So the MIXED bar (> 0.10 on either sample) is already known to hold on HG002 if the gate
passes. New in this row: the gate and the SAFE bar, which depend on the background.

## Background

For each chain, the two genomic positions 1,000 bp upstream and 1,000 bp downstream of the emitted 5' end, same strand;
background share = in-peak share over both positions pooled, per group. Chosen because refTSS peaks are narrow (p95
56 bp) and 1 kb leaves the promoter cluster; alternative promoters within 1 kb still count against the gate, which is
the conservative direction.

## Gate and bars (unchanged from row 1742)

- Gate: kept RefSeq-correct chains in a peak >= 3 x their background on each sample, else INVALID.
- SAFE: dropped in-peak share <= 2 x the dropped group's background on both samples.
- MIXED: dropped in-peak share > 0.10 on either sample.
- Otherwise: stays opt-in.

Script `bench/famcmp/trunc5_cage.py ... flank`; workspace `r1742_trunc5_cage` (`a.flank.txt`, `hg.flank.txt`). Limits as
row 1742. No default changes in this row.

## OUTCOME (2026-11-15): gate passes; MIXED as registered

| sample | group | chains | in peak | background (1 kb flanks) | ratio |
|---|---|---|---|---|---|
| A119b | dropped | 36,762 | 0.048 | 0.0034 | 14.1 |
| A119b | dropped, RefSeq-correct | 873 | 0.257 | 0.0143 | 17.9 |
| A119b | kept | 182,755 | 0.272 | 0.0056 | 48.6 |
| A119b | kept, RefSeq-correct (gate) | 37,842 | 0.344 | 0.0089 | **38.8** |
| HG002 | dropped | 16,889 | 0.125 | 0.0056 | 22.4 |
| HG002 | dropped, RefSeq-correct | 546 | 0.416 | 0.0266 | 15.7 |
| HG002 | kept | 136,144 | 0.493 | 0.0076 | 64.5 |
| HG002 | kept, RefSeq-correct (gate) | 35,199 | 0.463 | 0.0089 | **51.8** |

- Gate: 38.8 x and 51.8 x against 3 x, passes.
- SAFE: fails on both samples (14.1 x and 22.4 x against 2 x).
- MIXED: holds on HG002 (0.125 > 0.10), as row 1742's table had already shown; A119b alone (0.048) meets neither bar.
- Size of the real subset (excess over background): about 1,600 (A119b) and 2,000 (HG002) dropped chains with a
  CAGE-supported start, 3.3% and 3.0% of all supported multi-exon starts; the stage drops 16.7% and 11.0% of chains.

Reading as registered: a real subset is dropped; the stage should tag, not drop, and exempt supported starts. `trunc5`
is not a default. Not shown: which of the dropped chains are the real ones without CAGE; anything on gorilla.
