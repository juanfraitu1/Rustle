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
