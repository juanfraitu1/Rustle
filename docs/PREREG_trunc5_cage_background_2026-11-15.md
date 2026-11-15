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
