# Row 1836: floor on distinct records, read count on every record (pre-registered 2026-11-20)

User (2026-11-20), after row 1833 and its diagnosis (`docs/PREREG_count_all_reads_2026-11-20.md`): build the split. Counting
every read (`--keep-coordinate-duplicates`) raised precision on mouse, zebrafish, SIRV and OR6737 but flooded A119b, because half
of A119b's 2-record chains are identical copies of one record and passed the 2-read floor.

**The rule.** Opt-in `RUSTLE_DUP_COUNT=1` (new, `src/family.rs`, unit test
`dup_count_keeps_the_floor_on_distinct_records_and_counts_every_record`): a spliced chain is admitted on its DISTINCT
(start, end, chain) records, as by default; its `reads` is the number of records, as under `--keep-coordinate-duplicates`, so the
isoform fraction and the sub-chain ratio compare real counts. Unspliced records stay de-duplicated. Unset = byte-identical.

**Arms.** `base` and `dc` (`RUSTLE_DUP_COUNT=1`), stage `assemble`, `RUSTLE_PTC=0`, binaries built into `rustle_target_b` from the
commit of this document plus the switch (the shared binaries are untouched). Workspace `r1836_dup_count`. Registered before
either arm ran. `kd` of row 1833 is quoted beside the result.

**Data.** The seven datasets of row 1833 were seen under `kd` and shaped this rule: they are development here.
Check, not run under either variant before: mouse adrenal (ENCSR852HRM), mouse cortex (ENCSR862DMK), chimpanzee testis.
Unit: exact intron chains, `PREFIX.gtf`; `PREFIX.final.gtf` shown; gained and lost chains separately.

**Bars, `dc` against `base` on `PREFIX.gtf`:**
- **D3** on each of the three check datasets: precision higher, and at most 5% of base's correct chains lost.
- **A1** on A119b (development, the library the rule is built to protect): the same D3.
- **S1** (development): on the SIRV contigs of mouse ES precision is at least 0.80 (`kd` 0.862, base 0.589).
**Decision.** D3 on all three, A1 and S1 = candidate default; the user rules. Anything else = reported with the table.

**Expectation, written before the run.** The flood should go (the floor no longer sees copies). A119b's lost correct chains
may stay: its abundant chains are inflated by the same copies, and the isoform fraction then removes rare real isoforms next to
them. If A1 fails on that, the remaining candidate is the library check of row 1833 (the 2-record statistic), not this rule.

**Limits, stated in advance.** Single-exon transcripts keep de-duplicated counts while multi-exon counts grow, so the
single-exon read floor (a quantile of multi-exon support, row 1832) rises: single-exon models are reported per arm. The chain
score filter and the 5' truncation rule were fitted on de-duplicated counts. Two of three check datasets are mouse.

## Result (2026-11-20): D3 passes on the three check datasets, S1 passes, A1 FAILS on lost correct chains; not a candidate default

Jobs 4484252-4484261, both arms exit 0; `r1836_dup_count/score.NAME.txt`, `mono.NAME.txt`. `base` of the new binary is
byte-identical to the installed binary's on the seven datasets with a reference file. `PREFIX.gtf`, `dc` against `base`
(`kd` of row 1833 in the last column: precision, correct net):

| dataset | precision | correct chains (net) | correct gained / lost | wrong gained / lost | bar | `kd` |
|---|---|---|---|---|---|---|
| mouse adrenal | 0.540 -> 0.566 | 14,596 -> 14,344 (-1.7%) | +3 / -255 | +46 / -1,487 | **D3 PASS** | |
| mouse cortex | 0.465 -> 0.486 | 15,915 -> 15,733 (-1.1%) | +14 / -196 | +102 / -1,810 | **D3 PASS** | |
| chimpanzee testis | 0.351 -> 0.365 | 20,940 -> 20,473 (-2.2%) | +4 / -471 | +145 / -3,311 | **D3 PASS** | |
| A119b | 0.158 -> 0.172 | 39,049 -> 35,897 (-8.1%) | +84 / -3,236 | +2,879 / -38,025 | **A1 FAIL** (loss) | 0.139, -7.1% |
| SIRV contigs of mouse ES | 0.589 -> 0.877 | 53 -> 50 | 0 / -3 | 0 / -30 | **S1 PASS** | 0.862 |
| mouse ES | 0.331 -> 0.374 | 18,058 -> 17,365 (-3.8%) | +10 / -703 | +84 / -7,479 | dev | 0.372, -3.8% |
| mouse heart | 0.546 -> 0.570 | 13,234 -> 13,061 (-1.3%) | +11 / -184 | +67 / -1,244 | dev | 0.568, -1.1% |
| OR6737 | 0.321 -> 0.352 | 25,792 -> 24,570 (-4.7%) | +7 / -1,229 | +255 / -9,605 | dev | 0.346, -4.3% |
| fly simulated | 0.926 -> 0.939 | 7,082 -> 6,840 (-3.4%) | 0 / -242 | +2 / -130 | dev | 0.939, -2.4% |
| mouse simulated | 0.903 -> 0.927 | 27,731 -> 25,679 (-7.4%) | +35 / -2,087 | +13 / -975 | dev | 0.925, -6.4% |
| zebrafish simulated | 0.813 -> 0.836 | 12,574 -> 12,320 (-2.0%) | +6 / -260 | +16 / -497 | dev | 0.834, -0.3% |

`PREFIX.final.gtf` moves the same way (A119b 0.219 -> 0.232, 38,091 -> 35,175; adrenal 0.578 -> 0.605, -1.4%; chimpanzee
0.416 -> 0.430, -2.2%; OR6737 0.377 -> 0.416, -4.6%).

- **The flood is gone.** A119b: 35,146 fewer wrong chains than base (17%), where `kd` had 17,105 more; precision above base.
- **The loss stays, as expected.** A119b loses 3,236 correct chains (`kd` 3,646), mouse simulated 2,087: with real counts the
  isoform fraction and the sub-chain ratio remove rare isoforms beside abundant ones. That loss is the same under `kd` and
  `dc`, so it comes from the count, not from the floor. `dc` gains fewer correct chains than `kd` (the chains `kd` admitted on
  copies), which is why its net is 1.0-1.7 points lower on the simulated sets.
- So the split does what it was built for and turns the lever into a plain precision / recall trade that holds on every
  library: precision +0.013 to +0.043, correct chains -1.1% to -8.1%.
- Single-exon side effect, as stated: on real libraries 4-15% fewer truth single-exon transcripts found (A119b 234 -> 218,
  chimpanzee 104 -> 95); on the simulated sets a collapse (fly 343 -> 96, mouse 85 -> 30, zebrafish 144 -> 51), because the
  single-exon floor is a quantile of multi-exon support that now counts copies while unspliced records do not.

**Verdict by the rule: not a candidate default** (A1). The switch is kept as an opt-in (`RUSTLE_DUP_COUNT=1`), the safe form
of "every read counted". Open, the user's call: whether the isoform fraction (0.02) and the sub-chain ratio (0.3), both set on
de-duplicated counts, should be re-set under real counts; and counting unspliced records the same way to remove the
single-exon asymmetry.
