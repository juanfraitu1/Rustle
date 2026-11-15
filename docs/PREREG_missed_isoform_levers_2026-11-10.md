# PREREG r1616: which threshold costs us the isoforms StringTie or FLAIR find? (2026-11-10)

User question (2026-11-10): "the transcripts that we miss but stringtie or flair find, what is missing on our end?"
Gorilla OR6737 first (quick test bed), then A119b for any arm that passes. Written before any arm below was run.

## Already seen (r1614, descriptive, `bench/famcmp/missed_chains.py`)
4,867 RefSeq chains StringTie or FLAIR output and ours does not: 3,360 have one full-chain primary record, 398 none,
1,109 two or more (510 / 334 / 265 at 2 / 3-4 / 5+). `--polish-isoform-fraction` off: +947 correct chains, precision
0.352 -> 0.301 (marginal precision of the added chains about 0.06); misses left 186 / 73 / 154. ISM ratio relaxed: 5+ misses
265 -> 195. So the ceiling of arm A is known: N = 3 at most about +620, N = 5 about +110.

## Arms (assemble only, default seed pool, same BAM)
- A. `--polish-fraction-min-reads N`, N = 3, 5, 10 (a chain with >= N reads is never dropped as a minor fraction).
- B. `RUSTLE_MIN_CHAIN_READS=1` (single-read spliced chains admitted).
- C. polish off entirely (`--assembly-polish none`, `RUSTLE_POLISH_TSS=off`): descriptive, what all filters together remove.

## Measure
Multi-exon intron chains, exact match to RefSeq (GCF_029281585.2). correct = chains of ours in RefSeq; precision = correct /
distinct chains of ours. Marginal precision of an arm = (correct gained) / (chains added) against the default run.

## Bars (registered)
- PASS (candidate default): correct gained >= 100 and marginal precision >= the default run's precision (overall precision
  does not fall). An arm that passes on gorilla must pass the same bar on A119b before adoption.
- Otherwise FAIL as a default; the trade is reported and stays opt-in. No bar is tuned after the result.
- Descriptive, no bar: new columns per missed chain (identical chain of ours on the other strand; near-exact reads, every
  junction within 10 bp; minimum primary records over its introns; partial reads no chain of ours explains).

## Outcome (2026-11-10, gorilla OR6737, `r1616/`; default 73,348 chains, 25,794 correct, precision 0.352)
| arm | chains | correct | precision | gained | added | marginal |
|---|---|---|---|---|---|---|
| A, N = 3 | 78,944 | 26,198 | 0.332 | +404 | 5,596 | 0.072 |
| A, N = 5 | 74,832 | 25,898 | 0.346 | +104 | 1,484 | 0.070 |
| A, N = 10 | 73,483 | 25,801 | 0.351 | +7 | 135 | 0.052 |
| B, floor 1 | 127,014 | 29,006 | 0.228 | +3,212 | 53,666 | 0.060 |
| C, polish off | 155,147 | 27,392 | 0.177 | +1,598 | 81,799 | 0.020 |

- **All arms FAIL** the registered bar (marginal precision < 0.352). Not run on A119b. Defaults unchanged; all stay opt-in.
- Descriptive columns (`r1614/default.cols.tsv`, 4,866 misses): identical chain of ours on the other strand 0. One-full-read
  class 3,360: every intron in >= 2 primary records 2,453; some intron in one record only 584 (all FLAIR-only); near-exact
  reads (junctions within 10 bp) 323. Zero-full-read class 398: all StringTie, every intron in >= 2 records.
- Post hoc, not registered: floor-1 chains whose every intron is already in a default model: 2,443 added, 370 correct,
  marginal 0.151 (1 new intron 0.062; 2+ 0.049). Still below the bar.
- StringTie and FLAIR were both run annotation-free (`benchmark_collapse/run_*.sbatch`).

## Exploratory screen, post hoc (2026-11-10, `bench/famcmp/single_read_screen.py`, `r1616/screen.tsv`)
User: "is there any way to accept the single read but with another requirement?" The 53,773 chains floor 1 adds, by feature
(marginal precision): all 0.060; all are canonical; min intron reads >= 3: 0.102; every adjacent intron pair in >= 3 primary
records ("pair >= 3"): 0.100; both ends within 50 bp of a default model's ends: 0.083; no new intron: 0.151;
pair >= 3 + ends: 1,957 added, 354 correct, 0.181; pair >= 3 + ends + no new intron: 957, 204, 0.213; pair >= 5 + ends + no
new intron: 589, 137, 0.233. Calibration: default chains by read count are 0.170 (2 reads), 0.216 (3), 0.287 (4-5),
0.355 (6-10), 0.638 (11+). So the best single-read rules match the 2-3-read chains already accepted, not the average.
Chosen on gorilla after seeing the results: a frozen rule needs a held-out run (A119b) before any claim.
