# Multi-sample runs: pool, then split (design, 2026-11-02)

Goal: one run over several samples aligned to the same genome gives one transcript set and family catalog, plus per-sample counts.
User decisions: pool-then-split (not per-sample assembly + merge); equal-split transcript counts first.

## Interface
`tools/rustle_pipeline.sh STAGE --bam a.bam,b.bam [--samples A,B] --fasta G --out PREFIX ...`
- `--samples` defaults to the BAM basenames without `.bam`; names must be unique, `[A-Za-z0-9._-]+`; count must match `--bam`.
- One BAM (no comma): every product byte-identical to today; no new files.
- All BAMs must have identical `@SQ` lines (same genome); else exit 2 naming the first differing BAM.

## Data flow
1. pool (new, stage `pool`, cached by the BAM paths + sizes + mtimes): read-name collision check across samples (exit 2 with
   the first shared name and the two samples), then `samtools merge -@T` -> `PREFIX.pooled.bam` + `.bai`. The as-table and every
   later read-level stage use the pooled BAM.
2. assemble, families: unchanged, on the pooled BAM -> `PREFIX.gtf`, `PREFIX.final.gtf`, `PREFIX.fam.*`.
3. assign: once per sample BAM on the shared catalog -> `PREFIX.samples/<S>.assign.*` (same files as a single-sample assign).
   The pooled assign is not run.
4. counts (new): `utilities sample-counts --gtf PREFIX.final.gtf --bam a.bam,b.bam --samples A,B --out PREFIX.counts.tsv`.
   A primary, mapped read is compatible with a transcript of the same contig and strand-agnostic placement when the read's intron
   chain (N operations) is a contiguous block of the transcript's chain and the read lies inside the transcript's span
   +/- 10 bp; a single-exon read is compatible with a transcript whose exon contains it (+/- 10 bp). Each compatible read adds
   1/k to each of its k compatible transcripts. Output: `transcript_id gene_id <S1> <S2> ...` (counts to 2 decimals) and a
   footer-free summary to the log: per sample reads seen, compatible, unassigned.
5. copy table (new, in the driver): `PREFIX.copy_counts.tsv` = `family_id copy_index copy_tid` then per sample
   `<S>_reads <S>_abundance` from each `<S>.assign.quant.tsv` (`n_reads_hard`, `abundance`); missing rows = 0.

## Errors
Unreadable BAM, missing index (`.bai`/`.csi`), or a sample with zero primary reads on the genome: exit 2 naming it.

## Testing
- Unit (Rust): chain compatibility (contiguous block, gap, end tolerance, single-exon) and the 1/k split.
- Driver: single BAM byte-identical on a chr20 slice (all stage products vs today's build).
- PREREG r1503, panel chr17+chr19, A119b + HG002 Kinnex (CHM13), wall limit 2 h: (1) per-sample counts sum to compatible reads
  (exact); (2) pooled `=` vs RefSeq >= each single-sample run's; (3) pooled chain precision >= the lower single-sample value;
  report: per-sample assignment agreement (same copy) between shared-catalog and single-sample assign. Genome-wide after a pass.

## Out of scope
Different genomes per sample; transcript EM; per-sample-then-merge for many samples.
