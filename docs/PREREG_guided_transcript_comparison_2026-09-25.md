# Pre-registration: annotation-guided transcript comparison (figures 1-3, guided path) — 2026-09-25

Written 2026-09-25 before any guided number exists: no guided StringTie or guided FLAIR GTF is on this machine (every
`stringtie_guided_gtf` / `flair_guided_gtf` cell of `figures/samples.tsv` is `-`), so no table below has ever been
computed. The user supplies the GTFs later; the build computes them only then.

## Why

The tool comparisons of figures 1-3 are annotation-free on every side (user rule 2026-09-25): Rustle `assemble`
(reads + genome), StringTie `-L` without `-G`, FLAIR `bam2bed` + `collapse` without `-f`/`--gtf` (`flair correct`
skipped), IsoSeq collapse (verified in `benchmark_collapse/run_stringtie.sbatch`, `run_flair.sbatch`). A guided tool
run (the annotation given to the tool) answers a different question and must never share a panel, a table or a
ranking with the annotation-free runs.

## Rustle's counterpart

**None.** Rustle has no annotation-guided transcript assembly: `copy_assign` reads no annotation during assembly (its
`--gff` only tags catalog copies as annotated or not in the `--phase` copy graph), no `RUSTLE_*` switch feeds an
annotation to the assembler, and `tools/rustle_pipeline.sh assemble` takes `--bam --fasta` only (checked 2026-09-25
in `src/bin/copy_assign.rs` and the driver). Rustle's guided mode exists for families (fig. 7) and loci (fig. 8) only.
So the guided path compares the guided tools with the annotation; it has no Rustle row. If a guided Rustle mode is
ever added, it enters as its own registry column and this document gets an amendment first.

## What is computed (only when a sample's guided GTF exists)

For every registry sample with `stringtie_guided_gtf` and/or `flair_guided_gtf`, on the sample's usual genome-wide
evaluation scope (every contig its annotation covers), with the same code and parameters as the annotation-free path:

- fig. 1 (`fig1_guided`, drawn as `fig1g_guided`): gffcompare 0.12.10 intron-chain sensitivity and precision against
  the sample's annotation, and the sensitivity on the reference chains carried exactly by >= 2 primary reads (the
  `fig1_samples` rule).
- fig. 2 (`fig2_guided_categories`, `fig2_guided_filter`, drawn as `fig2g_guided`): SQANTI3 5.5.4 structural
  categories and rules-filter PASS.
- fig. 3 (`fig3_guided_bins`, drawn as `fig3g_guided`): the exact-chain fraction per tie bin and primary-support
  stratum, with the same cluster intervals; only on samples whose per-transcript counts the annotation-free build
  computes (the benchmark samples).

Every row carries `mode = annotation-guided`. The guided tables are separate files and separate figures.

## What is claimed

Nothing is ranked against Rustle. The guided numbers are reported as a description of the guided tools, with the
caveat printed in every table note and figure: **the guided tools were given the annotation they are scored
against, so their sensitivity and precision against it are not comparable with any annotation-free number.** A
guided-versus-annotation-free difference of one tool is reported as that tool's gain from the annotation, never as a
ranking of methods.

## What would change this

A guided Rustle transcript mode (then a like-for-like guided panel is designed and pre-registered as an amendment
before its first number), or a user decision to score the guided tools against an annotation they were not given.
