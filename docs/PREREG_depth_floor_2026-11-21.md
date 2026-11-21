# Row 1852: a read floor that scales with library depth (pre-registered 2026-11-21)

User (2026-11-21): do it, after row 1851 (a chain of N reads is worth far less in a deep library: 2 reads are 13% annotated
on mouse ES against 31% on heart; ES cut to heart's reported count has heart's precision). Measurement only, as a post-cut on
today's `final.gtf`; no assembler code. A build follows only on the user's word. Workspace `r1852_depth_floor` (`floor.py`).
Committed before the depths were counted and before any cut was made.

## The rule measured
Depth M = millions of mapped alignment records in the BAM (`samtools idxstats`, all contigs). Floor F = max(2, ceil(c x M));
a multi-exon transcript of `final.gtf` is dropped when its `reads` (every record counted, the default since row 1850) is
below F. Arms c = 0.25, 0.5, 1 reads per million, chosen a priori (1 per million is the round figure; row 1851's cuts on ES,
reads >= 5 to 10 at 8.7 M reads, fall between 0.5 and 1). Single-exon transcripts are not cut and not scored.

## Data
`r1849_dup_count_rerun/NAME/dc/x.final.gtf` (today's default) on the twelve datasets, truths of row 1849. Reused, said.

## Reading, fixed now
- An arm **pays** on a dataset when the chains it removes are at most half as precise as the dataset's `final.gtf` AND it
  keeps at least 90% of the correct chains. On a dataset where F = 2 the arm changes nothing and counts as neutral (left out
  of the count, reported).
- **Candidate** = an arm pays on at least 6 of 7 non-primate and 4 of 5 primate datasets, neutral ones aside, with at least
  three datasets changed. A candidate is put to the user as a rule to build (in the polish, before the filter).
- Reported per arm and dataset: F, chains removed, their precision, correct kept, precision after; and the spread of
  precision across the four mouse libraries before and after (does the floor make libraries comparable).

## Expectation
The simulated sets decide it. Mouse simulated has the depth of mouse ES by construction and a complete truth in which 2-read
transcripts are real, so a floor of 5-9 reads there removes true chains at the assembly's own precision: I expect it to FAIL
on mouse simulated (and fly, if deep) and to pay on the deep real libraries (ES, A119b, the gorillas). If so, the floor is a
statement about real-library noise that simulated reads do not carry, and the verdict by the bar is FAIL or undecided with
the table.

## Limits
- A post-cut on `final.gtf` is not the same as a floor inside the assembly (a removed chain's reads are not re-assigned).
- `idxstats` counts supplementary and secondary records too; depth is approximate.
- Real-library truths are incomplete: row 1851 showed most low-count "wrong" chains are canonical and seen by other tools.

## Result (2026-11-21): FAIL; a depth-scaled floor pays nowhere

`r1852_depth_floor/floor.txt` (`floor.py`, `depth.txt`). Depths M (millions of mapped records): fly 7.9, mouse simulated 29.5,
zebrafish 12.5, mouse ES 34.1, heart 4.6, adrenal 5.6, cortex 3.8, OR6737 10.7, KB3781 34.9, A119b 68.0, chimpanzee 4.0,
bonobo 3.7. (Row 1851's "8.7 M reads" for ES was the count of simulated reads; the ES BAM holds 34.1 M records.)

Per arm: floor F, chains removed (their precision), share of the correct chains kept, precision after.
| dataset | base | c = 0.25 | c = 0.5 | c = 1 |
|---|---|---|---|---|
| fly simulated | 6,753 / 0.973 | F = 2, neutral | F = 4: -1,008 (0.936), kept 0.860, 0.979 | F = 8: -1,856 (0.943), kept 0.741, 0.984 |
| mouse simulated | 25,077 / 0.970 | F = 8: -5,476 (0.935), kept 0.796, 0.979 | F = 15: kept 0.689, 0.982 | F = 30: kept 0.575, 0.986 |
| zebrafish simulated | 12,129 / 0.912 | F = 4: -2,111 (0.850), kept 0.852, 0.924 | F = 7: kept 0.754, 0.931 | F = 13: kept 0.623, 0.936 |
| mouse ES | 17,073 / 0.410 | F = 9: -19,593 (0.184), kept 0.789, 0.611 | F = 18: -25,055 (0.215), kept 0.684, 0.705 | F = 35: kept 0.580, 0.785 |
| mouse heart | 12,922 / 0.611 | neutral | F = 3: -4,413 (0.313), kept 0.893, 0.689 | F = 5: kept 0.764, 0.767 |
| mouse adrenal | 14,162 / 0.615 | neutral | F = 3: -4,770 (0.305), kept 0.897, 0.696 | F = 6: kept 0.732, 0.796 |
| mouse cortex | 15,490 / 0.520 | neutral | neutral | F = 4: -10,099 (0.253), kept 0.835, 0.658 |
| OR6737 | 24,315 / 0.420 | F = 3: -15,051 (0.198), kept 0.878, 0.499 | F = 6: kept 0.689, 0.602 | F = 11: kept 0.547, 0.678 |
| KB3781 | 23,407 / 0.440 | F = 9: -26,742 (0.257), kept 0.707, 0.626 | F = 18: kept 0.577, 0.695 | F = 35: kept 0.460, 0.753 |
| A119b | 34,882 / 0.236 | F = 18: -115,921 (0.153), kept 0.492, 0.536 | F = 35: kept 0.364, 0.619 | F = 69: kept 0.254, 0.698 |
| chimpanzee | 20,309 / 0.432 | neutral | neutral | F = 4: -22,647 (0.279), kept 0.689, 0.575 |
| bonobo | 17,383 / 0.358 | neutral | neutral | F = 4: -22,454 (0.200), kept 0.742, 0.494 |

- **No arm pays on any dataset** (removed chains at most half as precise as the assembly AND at least 90% of the correct kept).
  Where the removed chains are poor enough (mouse ES and OR6737 at c = 0.25, heart and adrenal at c = 0.5: 0.18-0.31) the
  floor still takes 10-21% of the correct chains; everywhere else the removed chains are more than half as precise as the rest.
- **On complete truths the low-count chains are true:** the chains a floor removes are 0.94 correct on fly and mouse simulated
  and 0.85-0.88 on zebrafish; precision rises by 0.006-0.024 for 14-43% of the correct chains. Expectation held.
- A floor is a plain trade: mouse ES 0.410 -> 0.611 for 21% of the correct chains; A119b 0.236 -> 0.536 for half of them.
- It does not make libraries comparable either: the four mouse libraries span 0.410-0.615 before, 0.520-0.705 at c = 0.5 and
  0.658-0.796 at c = 1. Depth in records is a poor proxy: ES and KB3781 have the same depth and differ in what a floor of 9
  removes (0.184 against 0.257 precise).

**Decision:** not a candidate; nothing built. A read cut stays what it is, a threshold the user can apply to `final.gtf`
(`reads` attribute) when precision matters more than recall. The other lead of row 1851, a yardstick for real libraries that
does not count every unannotated chain as wrong, is untouched.
