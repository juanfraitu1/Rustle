# PREREG r1617: single-read chains with pair support and known ends, held out on A119b (2026-11-10)

Rule chosen on gorilla OR6737 after seeing its results (r1616 screen: 1,957 added, 354 correct, 0.181; gorilla 2-read
default chains 0.170). User approved freezing and testing it ("yes, lets do that"). Written before any A119b run of it.

## Frozen rule R
A spliced chain carried by one read (in the `RUSTLE_MIN_CHAIN_READS=1` run, absent from the default run) is accepted when
1. every pair of adjacent introns occurs together, consecutively, in >= 3 primary records (one-intron chain: the intron in >= 3), and
2. its start and its end each fall within one 50 bp bin of the start / end of a default-run multi-exon model on its strand.
Evaluated from output files (`bench/famcmp/single_read_screen.py`); no Rust until it passes.

## Held-out test: A119b on CHM13, RefSeq chains (A119b.chm13v2.0_RefSeq.gtf.gz), default seed pool
- Bar 1: correct chains gained by R >= 100.
- Bar 2: share correct among R chains >= share correct among A119b default chains carried by exactly 2 reads.
- PASS = both. PASS means: worth porting to Rust as an opt-in, then a default decision on the full gffcompare table.
- Descriptive, no bar: the same counts for the other screen rows; 2-read default chains split by condition 1 at >= 5
  records (their own 2 plus 3 others).

## Outcome (2026-11-10, `r1617/screen.tsv`): FAIL as registered
- Floor 1 on A119b adds 122,194 chains, 3,330 correct (0.027).
- Rule R: 4,431 chains, 276 correct, 0.062. Bar 1 PASS (276 >= 100). **Bar 2 FAIL**: 0.062 < 0.066 (70,626 default
  2-read chains, 4,641 correct).
- Descriptive: pair >= 3 alone 0.036; ends alone 0.043; no new intron 0.057; R + no new intron 1,865 / 154 / 0.083;
  pair >= 5 + ends + no new intron 1,224 / 107 / 0.087 (both above 0.066, not the registered rule).
- 2-read default chains split by pair support >= 5: 0.069 vs 0.064, no separation.
- Not ported to Rust. Defaults unchanged.
