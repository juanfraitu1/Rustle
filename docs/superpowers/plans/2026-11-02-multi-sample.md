# Multi-sample Runs Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `rustle_pipeline.sh` accepts several same-genome BAMs, builds one assembly + family catalog on the pooled reads, and
writes per-sample assignments and counts.

**Architecture:** `rustle_pipeline.sh` gains one early branch: a comma in `--bam` execs the new `tools/rustle_samples.sh` with
the same arguments, so the single-BAM path is untouched. `rustle_samples.sh` validates and pools the BAMs (`samtools merge`),
runs `rustle_pipeline.sh all --to families` on the pooled BAM, then `rustle_pipeline.sh assign` once per sample under
`PREFIX.samples/<S>` (family products symlinked in), then the new `utilities sample-counts` and the copy table.

**Tech Stack:** bash, samtools (pixi rustle env), Rust (`src/utilities.rs`, new `src/sample_counts.rs`), noodles-bam.

**Spec:** docs/superpowers/specs/2026-11-02-multi-sample-design.md

## Global Constraints
- One BAM (no comma): every product byte-identical to today; no new files.
- `--samples` defaults to BAM basenames without `.bam`; unique, `[A-Za-z0-9._-]+`; count must equal the BAM count. Exit 2 otherwise.
- All BAMs: identical `@SQ` lines, an index (`.bai` or `.csi`), >= 1 primary mapped read; read names unique across samples. Exit 2
  naming the BAM (or the shared read name and both samples).
- Counts: primary mapped reads only; compatible = read intron chain is a contiguous block of the transcript chain and the read
  lies in the transcript span +/- 10 bp (single-exon read: inside one transcript exon +/- 10 bp); 1/k to each of k; 2 decimals.
- `PREFIX.copy_counts.tsv` columns: `family_id copy_index copy_tid` then per sample `<S>_reads <S>_abundance` from
  `PREFIX.samples/<S>.assign.quant.tsv` (`n_reads_hard`, `abundance`); missing = 0.
- v1 refuses `--candidates`, `--legacy-catalog`, and runs no `flag` stage in multi-sample mode (exit 2 with the reason for the
  first two; the log says flag is skipped). Ruling recorded in the spec's out-of-scope list.

## Review Focus
1. Two BAMs with the same basename in different directories -> duplicate default names -> exit 2 telling the user to pass `--samples`.
2. Paths with spaces in `--bam` -> every BAM path quoted; tested.
3. Secondary / supplementary / unmapped records -> never counted (test fixture contains each).
4. A read on the opposite strand of a transcript with the same chain -> compatible (reads are unstranded cDNA strands in the
   BAM; the transcript strand is not a filter) -> counted; tested.
5. Re-running with the same BAMs -> the pool is reused (cache key), not rebuilt; tested by mtime.

---

### Task 1: `utilities sample-counts`

**Files:**
- Create: `src/sample_counts.rs` (lib module, `pub mod sample_counts;` in `src/lib.rs`)
- Modify: `src/utilities.rs` (subcommand `SampleCounts { gtf, bam (comma list), samples (comma list), out }`)
- Test: unit tests in `src/sample_counts.rs`; `tests/sample_counts.rs` (fixture BAM built with samtools from an inline SAM)

**Interfaces:**
- Produces: `pub struct Tx { id: String, gene: String, chrom: String, exons: Vec<(i64, i64)> }` (0-based half-open, sorted);
  `pub fn compatible(read_exons: &[(i64, i64)], tx: &Tx) -> bool`; `pub fn count(txs: &[Tx], reads: impl Iterator<Item = (String, Vec<(i64, i64)>)>) -> (Vec<f64>, u64, u64)`
  (per-tx counts, reads seen, reads compatible); CLI writes `transcript_id\tgene_id\t<S...>`, log line
  `[sample-counts] <S>: seen N, compatible C, unassigned N-C`.
- Read exons come from the CIGAR: M/=/X/D extend an exon, N closes it (intron); I/S/H/P skip.

- [ ] Step 1: unit tests, assertions:
  `compatible(&[(100,200),(300,400)], tx [(50,200),(300,500)])` true; tx `[(50,200),(300,500),(600,700)]` true (block at start);
  tx `[(50,200),(310,500)]` false (acceptor differs); read `[(100,200),(300,400)]` vs tx `[(95,200),(300,405)]` true and vs
  tx `[(105,200),(300,400)]` true (read start 5 bp outside, within 10) but vs `[(115,200),(300,400)]` false; single-exon read
  `[(120,180)]` vs tx exon `(100,200)` true, vs `[(100,150),(300,400)]` false; `count` with one read compatible with 2 txs ->
  0.5 each, seen 1, compatible 1.
- [ ] Step 2: `cargo test --release --lib sample_counts` -> FAIL (module missing).
- [ ] Step 3: implement `compatible` and `count` (index txs by chrom, sorted by start; candidate txs overlap the read span).
- [ ] Step 4: integration test `tests/sample_counts.rs`: SAM with 5 records (2 primary compatible, 1 secondary, 1 supplementary,
  1 unmapped) + 1 reverse-strand primary with a matching chain; two samples; assert TSV rows and that per-sample column sums
  equal the log's compatible count (Review Focus 3, 4).
- [ ] Step 5: `cargo test --release --lib sample_counts && cargo test --release --test sample_counts` -> PASS; commit.

### Task 2: `tools/rustle_samples.sh` + the driver branch

**Files:**
- Create: `tools/rustle_samples.sh`
- Modify: `tools/rustle_pipeline.sh` (before the shared `--bam` check: `case "$BAM" in *,*) exec bash "$(dirname "$0")/rustle_samples.sh" "$STAGE" "$@_original";; esac`
  -- keep the original argv in `ORIG_ARGS=("$@")` before the parse loop)
- Test: `tests/pipeline_driver_samples.rs`

**Interfaces:**
- Consumes: Task 1's CLI.
- Produces: `PREFIX.pooled.bam(.bai)`, `PREFIX.pool.key` (one line per BAM: path, size, mtime), `PREFIX.samples/<S>.assign.*`,
  `PREFIX.counts.tsv`, `PREFIX.copy_counts.tsv`, `PREFIX.samples.log`. Stages accepted: `all` and `pool`; others exit 2
  ("multi-sample runs use all or pool").
- Per-sample assign: `mkdir PREFIX.samples`; symlink `PREFIX.gtf`, `PREFIX.families.gtf`, `PREFIX.fam.*`, `PREFIX.splice.mmi`
  (if present) to `PREFIX.samples/<S>.<same suffix>`; run `rustle_pipeline.sh assign --bam <S.bam> --out PREFIX.samples/<S>` with
  the user's remaining options.

- [ ] Step 1: tests (driver run with `--bin` pointing at a stub dir whose binaries record argv, as `pipeline_driver_units.rs`
  does): count mismatch -> exit 2 "samples"; name `a b` -> exit 2; duplicate basenames -> exit 2 mentioning `--samples` (RF1);
  `@SQ` differ -> exit 2 naming the 2nd BAM; shared read name -> exit 2 naming both samples; `--candidates` -> exit 2; space in a
  path accepted (RF2); second `pool` run leaves `PREFIX.pooled.bam` mtime unchanged (RF5); single `--bam` never execs
  rustle_samples.sh (stub log shows no `pool`).
- [ ] Step 2: `cargo test --release --test pipeline_driver_samples` -> FAIL.
- [ ] Step 3: implement both scripts.
- [ ] Step 4: tests PASS; `bash -n` both scripts; commit.

### Task 3: single-BAM identity + docs

- [ ] Step 1: chr20 run of `rustle_pipeline.sh all --to families` and `assign` with the new driver vs commit 9cd6313b's driver
  (same binaries): `cmp` every PREFIX.* product -> identical.
- [ ] Step 2: `docs/DATA.md` section "Multi-sample runs" (usage, products, v1 limits); `--help` text lines in both scripts; commit.

### Task 4: PREREG r1503 panel

- [ ] Step 1: write `docs/PREREG_multi_sample_2026-11-02.md` + register row (bars from the spec's Testing section) and commit
  BEFORE the run.
- [ ] Step 2: sbatch: A119b and HG002 Kinnex alone and pooled on chr17+chr19 (`--contig chr17,chr19`), 2 h wall; gffcompare
  `final.gtf` vs chains/ref.gtf; counts check; assignment agreement report.
- [ ] Step 3: record the verdict in the PREREG file and register; genome-wide run only after a pass.
