# Against the lab's own isoseq / StringTie / FLAIR runs (§6q7, 2026-09-19)

The lab already had annotation-free transcript sets for both project substrates, produced on a cluster:
`~/Desktop/isoseq_upload/` (isoseq collapse, human A119b and gorilla GGO_OR6737) and
`~/Desktop/benchmark_collapse/` (StringTie and FLAIR for both). **Those runs were reused as-is — nothing
was re-run.** Only our arm was produced here.

⚠**Substrate correction (register row 867).** Their runs use `A119b.t2t.bam` (1,104,846 chr20 records)
and `GGO_mm.bam`. The §6p8-§6q6 six-chromosome bakeoff used `human_testis.t2t.bam` — **188,864 chr20
records, ~6× shallower**. The two are not comparable: our chr20 output is 658 transcripts on one library
and 5,844 on the other. Every number below is on *their* BAMs.

Scoring is gffcompare against the RefSeq annotation of each species (CHM13 RefSeq for A119b, the
gorilla-native `GGO_genomic.gff` for GGO), not the tool-vs-tool comparison their `summarize.sh` runs.

## Human — A119b chr20 (1,104,846 records / 414,711 primary)

| arm | mRNAs | chains | chain Sn/Pr | tx Sn/Pr | beats (all 5 metrics) |
|---|---|---|---|---|---|
| **ours, shipped polish** | 5,844 | 1,064 | **24.8 / 19.9** | **23.3 / 18.2** | **StringTie, FLAIR** |
| ours, raw | 16,707 | 1,184 | 27.6 / 10.4 | 26.0 / 7.1 | – |
| ours, k3+majority raw | 20,699 | **1,259** | **29.4** / 8.2 | 27.7 / 6.1 | FLAIR |
| isoseq collapse | 64,384 | 1,253 | 29.2 / 3.0 | **28.0** / 2.0 | |
| StringTie | 5,731 | 861 | 20.1 / 16.8 | 19.0 / 15.1 | |
| FLAIR | 20,917 | 1,026 | 23.9 / 7.6 | 22.8 / 5.0 | |

### Where our chains go on deep data

| stage | mRNAs | chains | chain Pr |
|---|---|---|---|
| raw | 16,707 | 1,184 | 10.4 |
| mono floor + shadow only | 11,873 | **1,184 (−0)** | 10.4 |
| + ISM (ratio 0.7) | 11,709 | 1,129 (−55) | 16.7 |
| + isoform fraction 0.02 | 5,844 | 1,064 (−120) | 19.9 |

⭐The mono floor and shadow rule stay **free at 6× the depth** (0 chains, transcript precision 7.1 → 10.0).
⚠**The fraction filter costs 65 chains here against ~2 on the shallow library** — it is the depth-sensitive
component, and it was tuned on the shallow one.

⚠**Our raw output has 1,184 chains, already 69 below isoseq's 1,253 before any filtering**, so loosening
the polish cannot close that gap — the chains have to be assembled. `--read-isoform-k 3` +
`RUSTLE_JUNCTION_MAJORITY=1` does it: **1,259 chains from 20,699 transcripts against isoseq's 1,253 from
64,384**, beating isoseq on 4 of 5 (it misses transcript Sn by 0.3, 27.7 vs 28.0).

## Gorilla — GGO NC_073244.2, 80 Mb (473,231 records / 158,328 primary)

The thesis substrate. Scored against `GGO_genomic.gff` (5,936 reference transcripts on this contig).

| arm | mRNAs | chains | chain Sn/Pr | tx Sn/Pr | beats (all 5 metrics) |
|---|---|---|---|---|---|
| **ours, shipped polish** | 4,063 | 1,575 | **28.4 / 38.9** | **26.6 / 38.8** | **StringTie, FLAIR** |
| ours, raw | 7,407 | 1,625 | 29.3 / 24.7 | 27.5 / 22.0 | – |
| **ours, k3+majority + polish** | 5,298 | **1,688** | **30.4** / 31.9 | **28.5** / 31.9 | **isoseq, FLAIR** |
| ours, k3+majority raw | 10,933 | **1,787** | **32.2** / 17.7 | **30.3** / 16.4 | isoseq |
| ours, k3 + F=0.03 | 4,777 | 1,646 | 29.6 / 34.5 | 27.8 / 34.5 | FLAIR |
| ours, k3 + F=0.05 | 4,063 | 1,560 | 28.1 / 38.5 | 26.3 / 38.4 | StringTie, FLAIR |
| isoseq collapse | 20,643 | 1,655 | 29.8 / 9.0 | 28.1 / 8.1 | |
| StringTie | 3,735 | 1,374 | 24.7 / 37.0 | 23.2 / 36.8 | |
| FLAIR | 6,169 | 1,393 | 25.1 / 23.6 | 23.5 / 22.6 | |

⭐**On gorilla one setting beats isoseq outright on every metric** — k3+majority+polish, 1,688 chains at
31.9% chain precision against isoseq's 1,655 at 9.0%, from a quarter of the transcripts.

## The frontier, and what is still missing

**FLAIR is dominated on both species by every one of our settings that beats anything.** isoseq collapse
and StringTie sit at opposite ends of the recall/precision frontier — isoseq buys chains with 3-9%
precision and tens of thousands of transcripts, StringTie buys precision with 20-35% fewer chains — and we
sit between them, able to dominate **either** by choice of setting but **not both at once** (row 868):

| | to beat isoseq | to beat StringTie | our best nearby |
|---|---|---|---|
| A119b chr20 | ≥1,253 chains | ≥16.8% chain Pr | 1,259 chains @ 8.2% — needs 1,253 chains from ≤7,458 transcripts, has 20,699 |
| GGO NC_073244.2 | ≥1,655 chains | ≥37.0% chain Pr | **1,647 @ 34.5 — 8 chains and 2.5 points short** |

**Gorilla is within reach and human is not.** Closing gorilla needs ~330 non-matching transcripts removed
at F=0.03 without losing a chain; closing human needs ~13,000 removed.

⛔`--polish-fuzzy-junction 5` re-tested on gorilla: costs 126 chains (1,646 → 1,520) at the same
precision, and 2 bp is neutral (+0.2) — the §6q6 canonical-GT-AG explanation holds on a second species
(register row 869).

## Ranked levers

1. **Make `--polish-isoform-fraction` depth-aware.** It is the one filter whose cost scales with coverage
   (65 chains on A119b chr20, ~2 on the shallow library), and it is the binding constraint on both species.
2. **Turn on `--read-isoform-k 3` + `RUSTLE_JUNCTION_MAJORITY=1` for deep input.** That is where the
   +75 (human) / +212 (gorilla) raw chains come from; without them isoseq's recall is unreachable.
3. **Transcript-Sn against isoseq on human is 0.3 points** (27.7 vs 28.0) — a terminal-boundary problem,
   not a chain problem, and the §6p4 5'-truncation finding is the place to look.
