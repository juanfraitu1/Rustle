# PREREG: a second model for a chain with a second supported end (row 1993, 2026-11-28)

User: try the next lever of row 1992. There, a 3' end short of the annotated one was at a polyA site in 0.55-0.58 of the
cases (HG002, three tools): many chains use two polyA sites, and Rustle writes one model per intron chain. FLAIR writes up
to two starts and two ends per chain. The pipeline's `altends` stage only adds attributes (the commonest read end when it
is not the emitted one); it writes no model and was never judged against polyA sites.

## The rule, fixed before any result
From the reads with exactly the chain, after the primary end is set by the supported-end rule of row 1992:
a SECONDARY 3' end = a cluster of read ends (within 10 bp of its commonest position) that lies at least 50 bp inside the
primary end, holds >= 3 reads and >= 20% of the chain's reads, and is not followed by an A-rich stretch of the genome
(>= 12 A in 20 bp or 6 in a row: internal priming, the test `alt_ends.py` already has). At most one per chain: the best
supported such cluster. The same rule mirrored for a secondary 5' end is MEASURED beside it (no A test there), not proposed:
a start inside the primary one is what a cut-short read looks like.
`bench/o2_truth/alt_end_models.py`.

## Judges
- human (WTC11 r1, r2, HG002): within 25 bp of a polyA site of the atlas (3') / CAGE peak (5'); the primary ends' rate and
  the background 500 bp away beside it;
- every library with an annotation: within 25 bp of ANY annotated transcript end on that strand (an annotated polyA site of
  some isoform), primary ends beside it;
- simulated fly, mouse, zebrafish: the truth has one end per chain in 98-99.8% of chains, so a secondary end there is a
  false one unless it is within 25 bp of a true end of that chain.
## Claims
- **A1:** secondary 3' ends are at a polyA site at least 0.60 of the time on each of the three human libraries, and at
  least 0.80 of the primary ends' rate.
- **A2:** on each simulated set secondary 3' ends are written for at most 1% of the chains.
- **A3:** on the mouse libraries (ES, cortex) secondary 3' ends are within 25 bp of an annotated transcript end at least
  0.80 as often as the primary ends.
- **A4 (the A test earns its place):** the candidates it removes are at a polyA site at most half as often as the ones it keeps.
- **A5 (5', measured):** secondary 5' ends are near a CAGE peak less often than secondary 3' ends are at a polyA site.
Prediction: A1 holds, A2 holds, A3 holds, A4 holds, A5 holds.

## Result of the registered rule (job 4517926)
Secondary 3' ends: chains that get one; at a polyA site (human) or within 25 bp of an annotated transcript end (mouse,
simulated); the primary ends' rate beside it.
| library | chains with a secondary 3' end | rate | primary ends | removed as A-rich: n, rate |
|---|---|---|---|---|
| WTC11 r1 | 5,630 (12.7%) | 0.606 | 0.788 | 53, 0.170 |
| WTC11 r2 | 3,715 (11.0%) | 0.577 | 0.795 | 33, 0.242 |
| HG002 | 12,906 (10.8%) | 0.514 | 0.459 | 9,390, 0.073 |
| mouse ES | 6,250 (14.1%) | 0.421 | 0.753 | 290, 0.448 |
| mouse cortex | 3,182 (10.1%) | 0.478 | 0.800 | 621, 0.407 |
| fly / mouse / zebrafish simulated | 49 / 83 / 33 (0.68 / 0.30 / 0.23%) | 0.98 / 1.00 / 0.88 | 0.99 / 1.00 / 0.97 | 14 / 13 / 4 |
- **A1 FAILS** (0.606, 0.577, 0.514; bar 0.60 on each and 0.80 of the primary rate: met on WTC11 r1 for the first only).
- **A2 HOLDS** (under 1% of the chains; on simulated sets these are true second ends of the truth: 44 of 46, 50 of 63, 19 of 30).
- **A3 FAILS** (0.56 and 0.60 of the primary rate).
- **A4 HOLDS** on the human libraries (0.07-0.24 against 0.51-0.61); HG002 has 9,390 A-rich candidates, WTC11 about 40.
  On mouse the removed ones are as often annotated as the kept ones.
- **A5 FAILS, the other way round from my prediction:** secondary 5' ends are near a CAGE peak 0.941 / 0.947 / 0.630
  against 0.847 / 0.836 / 0.594 for primary 5' ends: the inner cluster is the main promoter and the supported-extended
  start of row 1992 an upstream minor one. On simulated sets they are false (5% of the chains, 4-18% at a true start).
The rule as registered is not good enough to write a second model: four in ten of its ends are not at a known polyA site.

## Development (WTC11 r1, mouse ES only) and rule 2, written before the other libraries are re-read
What separates the real ones on the two development libraries is the SHARE of the chain's reads at the cluster: 0.2-0.3
0.47 / 0.31; 0.3-0.5 0.56 / 0.37; 0.5-0.8 0.68 / 0.51; >= 0.8 0.83 / 0.63. A polyA signal (AATAAA / ATTAAA 5-40 bp
upstream) adds 0.15 / 0.06 (0.640 against 0.489; 0.437 against 0.379); the number of reads and the distance add little.
At >= 50% the inner cluster is the chain's MAIN end, and the primary end of row 1992 the less used, more distal one.
**Rule 2:** the registered rule, with the cluster holding >= 50% of the chain's reads AND a polyA signal upstream.
Judged on WTC11 r2, HG002 and mouse cortex (their overall tables above were seen; no rule was fitted on them):
- **B1:** at a polyA site >= 0.65 on WTC11 r2 and HG002, and >= 0.80 of the primary ends' rate on WTC11 r2.
- **B2:** mouse cortex: within 25 bp of an annotated transcript end >= 0.70 of the primary ends' rate.
- **B3:** simulated sets: still under 1% of the chains.
Prediction: B1 holds, B2 fails narrowly, B3 holds.

## Rule 2 judged (`rule2.py` in the workspace) and built
| library | chains with a second 3' end | rate (polyA site / annotated end) | registered rule |
|---|---|---|---|
| WTC11 r1 (development) | 2,085 (4.7%) | 0.722 | 0.606 |
| mouse ES (development) | 2,009 (4.5%) | 0.546 | 0.421 |
| WTC11 r2 | 1,369 (4.1%) | 0.698 | 0.577 |
| HG002 | 4,095 (3.4%) | 0.656 | 0.514 |
| mouse cortex | 1,190 (3.8%) | 0.582 | 0.478 |
| fly / mouse / zebrafish simulated | 13 / 18 / 6 (0.18 / 0.07 / 0.04%) | 0.92 / 1.00 / 1.00 | |
- **B1 HOLDS** (0.698 and 0.656; 0.88 of the primary rate on WTC11 r2). **B2 HOLDS** (0.582 = 0.73 of 0.800; my
  prediction was wrong). **B3 HOLDS.**
Built: `tools/end_choice.py ... GENOME.fa`, opt-in `RUSTLE_ALT_ENDS=1` in `all` (id `.a3`, `alt_end "3"`, `alt_end_reads`).
Check (job 4517936): WTC11 r2 1,368 second models, cortex 1,190, fly 13; every second model shares its chain with the
first one and differs only in the last exon, exons valid. On WTC11 r2 the second models are at a polyA site 0.698 and
THEIR first models 0.711: for these chains both ends are real. Fly through the whole `all`: no warning, 13 second models.
Gorilla KB3781 with the option: 2,477 second models; families with >= 2 copies 260 = 260, the same sets of loci, 724
copies. cargo test 1,591 passed. Default unchanged (off); the user decides.
Not done: a second 5' end (real on human by CAGE, false on simulated reads); the label and family-mode stages were not
run with two transcripts per chain beyond the families stage.
