# Publication figures

Reproducible figures for the Rustle paper. Every plotted number lives in a small tidy table under `data/`, with a
provenance header (generator, git commit, date, and each input file with its size and mtime). The figures render
from those tables alone, on any machine; the tables regenerate from the raw inputs with one command.

```bash
python3 figures/make.py list                  # what each figure claims and which tables it reads
python3 figures/make.py plot all              # render out/<figure>.{pdf,png,svg} from data/*.tsv (seconds)
python3 figures/make.py check                 # every table present with provenance; every figure renders
cp figures/inputs.example.tsv figures/inputs.local.tsv   # then edit the paths for this machine
python3 figures/make.py data figN             # regenerate one figure's tables from the raw inputs (see below)
python3 figures/make.py data fig2 --set fig2_arms=rustle,stringtie   # override one input key for a call
python3 figures/make.py data fig4 --recorded  # figures 4-7: tabulate the recorded development runs (light) into data_recorded/
python3 figures/make.py plot fig4 --data figures/data_recorded --out /tmp/recorded   # render those
python3 figures/make.py samples --verify      # the sample registry, resolved; files, indexes and contig names checked
python3 figures/make.py runs --dry-run        # the genome-wide run queue, every sample x stage, with estimated cost
python3 figures/make.py runs --sample chimp_PTR --stage assemble   # run one genome-wide stage (HEAVY; cached)
python3 figures/make.py runs --sample chimp_PTR --stage catalog    # bounded call: exit 75 = run the same command again
python3 figures/make.py runs --sample S --stage ST --restamp PROOF.tsv   # re-stamp a stale stage from a cmp proof (no run)
```

## Requirements

Python 3 (built with 3.14) with matplotlib (3.10), numpy, pysam and scipy; on PATH: samtools, minimap2, gffcompare;
mmseqs and blastp only for the supplementary protein material (the translated-search comparator of fig. 6s-protein,
the protein-homology secondary references when requested, and the manual protein step of fig. S-P); SQANTI3 5.5 in its
own conda env (`sqanti3_env`, `sqanti3_dir`); the Rust binaries (`bin`; the `families` stage needs a `mcl_families`
that writes the copy table, `--help` naming `<out>.copies.tsv`). Text is set in Arial when the font files are found
(`figures/fonts/`, the Windows or macOS font folders; `figlib.FONT_DIRS`), else DejaVu Sans.

## Layout

| file | role |
|---|---|
| `figlib.py` | style (7 pt sans, 89/183 mm widths, editable PDF text), the validated tool palette, tidy-table IO with provenance, the UpSet helper, `run`/`fresh`/`work_dir` for cached heavy steps |
| `samples.tsv`, `samples.py` | the sample registry (6 samples, 4 species; cells name keys of the inputs file) and the genome-wide run cache (`make.py samples`, `make.py runs`; see below) |
| `assembly.py` | shared provisioning for the assembly figures (and the transcript-mode wording, `MODE_DENOVO_METHODS`, `GUIDED_NA`): our assembler in both seeding configurations from the run cache (`${work}/runs/<sample>/<sample>.gtf`, `<sample>.primary.gtf`), every arm and the annotation restricted to the evaluation contigs, gffcompare runs and parsers (`${work}/assembly/<species>/eval_<contigs>/`) |
| `_sqanti.py`, `_o2.py`, `_o1.py`, `_o1_recovery.py`, `_liftoff.py` | helpers of figures 2, 4-5, 6, 7 and 8 (`_liftoff.copy_pairs` / `pair_families`: the Liftoff copy-pair family reference of figs 6s, 7 and 8) |
| `_lrc.py`, `test_lrc.py` | %LRC (LRGASP long-read coverage: the share of each model's exonic length under primary aligned read bases), a reporting metric only, no figure reads it: `python3 figures/_lrc.py union|score|tables --sample ID` (cached under `${work}/lrc/`; docs/archive/2026-09/PREREG_lrc_metric_2026-09-29.md); `python3 figures/test_lrc.py` |
| `fig_*.py` | one module per figure: `META`, `build(cfg, data_dir, force)`, `plot(data_dir, out_dir)` |
| `captions/<fig>.md` | caption, claim, provenance and caveats of each figure |
| `inputs.example.tsv` | every raw input the builders read (this machine's paths), and every optional key as a `# key<TAB>value` comment; copy to `inputs.local.tsv` |
| `data/` | the plotted tables (21 tables, all rebuilt 2026-09-25 from the current binaries; `data_recorded/` holds `--recorded` tables and is not tracked) |
| `out/` | rendered figures (PDF with editable text, PNG at 300 dpi, SVG) |

## Rules every figure follows

- **Species are never pooled.** Gorilla (OR6737 testis, GGO reference) and human (A119b, CHM13 v2.0) are separate
  panels; every read-derived table has a `species` column.
- **Apples-to-apples comparisons, one mode per panel.** The tool comparisons of figures 1-3 are annotation-free
  (de novo) on every side: Rustle assemble (reads + genome), StringTie 3.0.1 `-L` without `-G`, FLAIR 3.0.1 collapse
  without annotation (`flair correct` skipped), IsoSeq collapse 26.2.0, all run on the same BAMs (checked in
  `benchmark_collapse/run_stringtie.sbatch` and `run_flair.sbatch`). Every tool-comparison figure prints that mode on
  its first line and in its caption. Annotation-guided tool runs (StringTie `-G`, FLAIR with the annotation) are
  registered per sample (`samples.tsv` columns `stringtie_guided_gtf`, `flair_guided_gtf`; empty today) and are
  scored by the same code into separate tables and figures (`fig1g_guided`, `fig2g_guided`, `fig3g_guided`) that
  build only when such a GTF exists; until then the figures say "guided comparison: not available (guided
  StringTie/FLAIR GTFs not supplied)". Rustle has no annotation-guided transcript assembly, so the guided figures
  have no Rustle row (docs/archive/2026-09/PREREG_guided_transcript_comparison_2026-09-25.md). Figure 7 compares Rustle's two
  FAMILY modes with each other (Rustle-internal, not a tool comparison); figure 8 compares Rustle's guided locus
  search with Liftoff like for like and scores the de novo loci against Liftoff as a reference. "Mode" always
  names its level (transcript, family, locus): GLOSSARY *Modes*.
- **One family definition, external references.** Rustle's families are its ONE default de novo family definition
  (user decision 2026-09-25): reads → loci seeded with secondary alignments within 2% of the best score → one
  representative per locus (its "positional exon sum") → families (`tools/rustle_pipeline.sh families` =
  `mcl_families --from-gtf`, exon-sum ≥ 0.60, MCL 2.8). Copy assignment (figs 4–5) consumes the same families through
  the stage's copy table (`<id>.fam.copies.tsv`). The `gw_family_catalog` copy catalog (`catalog` stage) is legacy:
  kept runnable, scored only as a labelled comparison. Protein is not in the default: the manual extra-sensitive step
  (`tools/protein_attach.py`) has its own supplementary figure (S-P). Main family figures use external references only
  (Ensembl Compara, Soto 2025 labelled not independent, Liftoff copies); the translated protein search and the
  protein-homology families appear only in supplementary figures, labelled comparator / secondary.
- **Same contigs for every method.** Figures 1-3 score every method genome-wide on every contig the sample's
  annotation covers (human: all but chrM, which CHM13 RefSeq does not annotate); the current human tables are still
  chr20-22 until the rebuild (each caption's Status box). Figures 4-8 name their scope.
- **Rustle's arm is the shipped pipeline** (`tools/rustle_pipeline.sh assemble`, loci seeded with secondary
  alignments within 2% of the molecule's best alignment score). "Rustle (primaries only)" is the same run with `--no-seed-secondaries`.
- **Palette** (validated for colour-vision deficiency on all pairs): Rustle blue, StringTie orange, FLAIR aqua,
  IsoSeq violet; primaries-only Rustle is hatched blue. Series are always labelled, never identified by colour alone.
- **Heavy steps run in the foreground, one at a time**, cached under `work` (default `/mnt/linuxdisk/tmp/rustle_figures`).
  Every step is cached, so re-running an interrupted `make.py data figN` resumes where it stopped. Never run
  `data all --force` (`--force` recomputes a figure's own steps; the genome-wide assemblies are recomputed only by
  deleting them).

## Samples and genome-wide runs

`samples.tsv` lists every sample; its cells are `${key}` references into the inputs file, which holds every machine path.
Numbers from different samples or species are never pooled; `id` names the sample in every table that has one.

| id | alias | species | tissue | reference | BAM records | lab baselines |
|---|---|---|---|---|---|---|
| `human_A119b` | `human` | human | unknown (the notes disagree) | T2T-CHM13 v2.0 | 68.0 M | StringTie, FLAIR, IsoSeq |
| `human_testis` | | human | testis (ENA ERR13885926) | T2T-CHM13 v2.0 | 9.1 M | none |
| `gorilla_OR6737` | `gorilla` | gorilla | testis | mGorGor1 GCF_029281585.2 | 10.7 M | StringTie, FLAIR, IsoSeq |
| `gorilla_KB3781` | | gorilla | fibroblast cell line (the reference individual) | mGorGor1 GCF_029281585.2 | 34.9 M | none |
| `chimp_PTR` | | chimpanzee | unknown | mPanTro3 GCF_028858775.2 | 4.0 M | none |
| `orangutan_PPY` | | orangutan | unknown | mPonPyg2 GCF_028885625.2 | 16.9 M | none |

Registry columns beyond the paths: `stringtie_gtf` / `flair_gtf` / `isoseq_gff` = the lab's annotation-free runs
(A119b and OR6737 only); `stringtie_guided_gtf` / `flair_guided_gtf` = annotation-guided runs of the same tools,
`-` for every sample until the user supplies them (a path, or a `${key}` of the inputs file). There is no Rustle
guided column: Rustle has no annotation-guided transcript assembly.

`make.py runs` runs every stage of `tools/rustle_pipeline.sh` GENOME-WIDE on the whole BAM of each sample, into
`${work}/runs/<sample>/` with the driver PREFIX `<sample>` (`assemble_primary`, `families_primary`:
`<sample>.primary`); the driver's `PREFIX.cache/` stays on. Stages: `assemble`, `assemble_primary`
(`--no-seed-secondaries`), `families` (THE default de novo family definition on `<sample>.gtf`, with its copy table
`<sample>.fam.copies.tsv` / `.fa` that copy assignment consumes), `families_primary` (the same on the primaries-only
assembly; OPTIONAL: only the supplementary fig. 6s-seeding reads it), `catalog` (LEGACY copy catalog, `--piecewise`;
no main figure needs it), `assign` (it reads the families' copy table since 2026-10-02, so it needs `families`; the
legacy catalog only with the driver's `--legacy-catalog`, not passed here; the driver's opt-in `candidates` stage, ruling
R14, is not a run-cache stage), `index` (a splice minimap2 index, only where none exists: chimpanzee, orangutan)
and `flag` (`--gff` = the sample's annotation, `--confirm` = the gorilla haplotype assemblies). A stage is skipped while its stamp `<stage>.key` matches the current key (BAM/FASTA/index/
annotation fingerprints, sha1 of the binaries, the driver's code hash, minimap2 version, `RUSTLE_*` environment,
upstream key); every run appends wall time and peak RSS to `<stage>.time`. It refuses to start while another heavy
process runs (one heavy run at a time); it does not take the machine lock itself, so wrap every call in
`flock -w 900 /mnt/linuxdisk/tmp/rustle_heavy.lock`. `--dry-run` prints the queue with an estimated wall time and peak
RSS per stage and the recorded runs each estimate rests on. Figure code reads the products with
`samples.product(cfg, sample, stage, name)` (a path; runs nothing) and never runs a stage implicitly.

**Bounded calls (exit 75 = run the same command again).** `catalog` builds its representatives one piece per unit
(`--piecewise`; a contig, or read-free sub-ranges of about 1 M records) within `catalog_budget_s` (420 s) per call,
then merges them in a later call. `families`, `families_primary` and `catalog` run their genome-wide minimap2
all-vs-all through `tools/mm2_shard.sh` (set in the child's environment; `runs_mm2_wrapper`): one shared index, query
shards mapped one at a time, each call bounded by MM2_SHARD_BUDGET_S (`runs_mm2_budget_s`, 480 s) and a deadline
(`runs_call_budget_s`, 560 s after the call starts). The shards concatenate to the single-run PAF byte for byte (cmp
on chr16 for both binaries), so the wrapper is not part of the stage key. `make.py runs` exits 75 when a call made
progress and work remains (pieces, the merge, shards, or `--max-stages` reached), and fails when an all-vs-all made no
progress in a whole call (one shard or the index does not fit the budget; raise the budget for that call).

**Re-stamping without a run.** When the binaries or the driver change but a stage's outputs provably do not,
`make.py runs --sample S --stage ST --restamp PROOF.tsv` writes the new stamp without running the stage. PROOF.tsv lists
`product<TAB>NAME<TAB>CACHED_PATH<TAB>FRESH_PATH` for every product of the stage (FRESH = made by the current code,
e.g. in a scratch PREFIX) and optionally `command<TAB>how FRESH was made`. It is refused unless the stage is stale only
through code fields (binaries, driver code, minimap2 build, or an upstream key whose stages are all fresh), every
product is listed once, every FRESH file is newer than the current binaries, and every pair is byte-identical (size and
sha1 computed at re-stamp time). The proof (path, sha1, text, per-product sha1, the old digest and the fields bridged)
is recorded in the stamp (`restamped`) and appended to `<stage>.restamp.log`.

The two assemblies figures 1-3 used (`${work}/assembly/{gorilla,human}/rustle*.genome.*`) were moved into the run
cache on 2026-09-25 (same inodes; the old names are symlinks) and adopted without re-assembly.

## Figures

| fig | module | shows | tables |
|---|---|---|---|
| 1 | `fig_intron_chain.py` | annotation-free comparison: gffcompare intron-chain sensitivity / precision of Rustle (both configurations), StringTie, FLAIR, IsoSeq collapse; sensitivity by read support; paired test on the >= 2-read stratum. Supplementary `fig1s_samples` (Rustle on all six samples), `fig1g_guided` (annotation-guided tools only, when registered) | `fig1_gffcompare`, `fig1_support`, `fig1_paired`; `fig1_samples`, `fig1_guided` |
| 2 | `fig_sqanti.py` | SQANTI3 structural categories, rules-filter PASS and FSM counts, same GTFs, annotation and mode as figure 1. Supplementary `fig2s_samples`, `fig2g_guided` | `fig2_sqanti_categories`, `fig2_sqanti_filter`, `fig2_sqanti_subcategories`; `fig2_samples_*`, `fig2_guided_*` |
| 3 | `fig_secondary.py` | annotated chains recovered by tie fraction of their reads; transcripts only reachable through secondary alignments within 2% of the best score; one example locus (annotation-free methods). Supplementary `fig3g_guided` | `fig3_ref_tie`, `fig3_bins`, `fig3_gain`, `fig3_example`; `fig3_guided_bins` |
| 4 | `fig_assignability.py` | per-read fate of simulated reads (MAPQ > 0, tied, decisive site, assigned, correct), one bar per sample; supplementary UpSet grid `fig4s_assignability_upset` | `fig4_assignability_upset` |
| 5 | `fig_assign_accuracy.py` | copy-assignment fraction assigned and fraction correct by identity band vs the aligner; hard-locus benchmark on real reads vs the annotation-free tools; supplementary `fig5s_margin_rule` (the alignment-score margin rule) | `fig5_assign_accuracy_bands`, `fig5_hard_locus`, `fig5_hard_locus_transcripts`, `fig5s_margin_rule_*` |
| 6 | `fig_family_spectrum.py` | the default de novo families across the protein-identity spectrum of Ensembl Compara pairs (a direct nucleotide alignment vs the families; precision against Compara); human samples. Supplementary `fig6s_seeding` (loci seeded with secondary alignments vs primaries only, each through the families stage, vs Compara ≥ 90% and Liftoff copy pairs; development: gorilla chr20 vs protein-homology families, secondary) and `fig6s_protein` (the translated protein search as a comparator, chr16) | `fig6_gw_*` (genome-wide) or `fig6_chr16_*` (development); supplementary `fig6s_seeding`, `fig6s_protein_tiers`, `fig6_gorilla_*` |
| 7 | `fig_family_recovery.py` | family recovery of Rustle's two FAMILY modes (de novo = the default definition, guided = the same rule on the annotation; Rustle-internal, not a tool comparison) against Compara families (primates), Soto 2025 (not independent), the NPIP set and Liftoff copy pairs (every species; de novo only): pairwise and one-to-one bipartite sensitivity / precision / F, per-family outcomes. Supplementary `fig7s_protein_homology` (protein-homology families, secondary) | `fig7_gw_*` (genome-wide) or `fig7_summary`, `fig7_per_family` (development) |
| 8 | `fig_loci.py` | loci in the Liftoff framework: Liftoff's self-lift (`-copies`) as the annotation-guided locus baseline, Rustle's guided search like for like, de novo loci, the default families (legacy catalog rows kept as a comparison) and missing-copy flags scored against it | the tables `fig_loci.py` lists (`fig8_*`) |
| 9 | `fig_loop.py` | the closed loop: tied reads assigned among the default families' copies (union test), given to their copy, re-assembled (`docs/archive/2026-09/PREREG_tied_read_loop_2026-09-25.md`; captions/fig9.md) | the tables `fig_loop.py` lists |
| S-P | `fig_protein_supp.py` | the manual extra-sensitive protein step (not the default) scored against Compara and Liftoff (`docs/archive/2026-09/PREREG_protein_attach_2026-09-25.md`; captions/figS_protein.md) | `figS_protein_*` |

The per-figure captions (claims, panels, n, provenance, caveats) are in `captions/`.

## Build order and cost: the genome-wide queue (2026-09-25, 4 threads, this WSL2 box)

Every heavy call runs in the foreground under the machine lock, one at a time, <= 10 min and <= 20 GB:
`flock -w 900 /mnt/linuxdisk/tmp/rustle_heavy.lock <command>`; if the lock is not acquired, skip and try later. Every
step is cached, so a repeated call resumes. "Repeat" = run the same command again while it exits 75. Costs are
measured (wall / peak RSS from `<stage>.time` and the phase-2 reports) unless marked "est.".

### 1. Pipeline runs per sample (`make.py runs`; order chimp_PTR, human_testis, orangutan_PPY, gorilla_OR6737, gorilla_KB3781, human_A119b)

| step | command | cost per sample (measured unless est.) |
|---|---|---|
| assemble | `make.py runs --sample S --stage assemble` | first run adds the best-AS table (`as_table`, one pass over the BAM: A119b 544 s). Re-run with the table present: A119b 452 s / 2.6 GB, KB3781 169 s / 1.5 GB, PPY 96 s / 1.3 GB, OR6737 79 s / 1.3 GB, PTR 49 s / 0.9 GB, testis 30 s / 0.8 GB |
| assemble_primary | `--stage assemble_primary` | A119b 418 s / 2.0 GB, KB3781 156 s, PPY 91 s, OR6737 75 s, PTR 37 s, testis 22 s (all <= 1.2 GB except A119b) |
| index | `--stage index` (chimp, orangutan only) | done: PTR 182 s / 20.6 GB, PPY 438 s / 20.8 GB; nothing may run beside it |
| catalog | `--stage catalog`, repeat until exit 0 | pieces: PTR 26 pieces, 857 s over 4 calls; testis 25 pieces in 1 call (~110 s); PPY 31 pieces of 74-259 s / 2.8-5.9 GB each, 3-4 per call (~50 min); A119b about 80 pieces of 60-120 s (1.5-2.5 h est.; chr1 cut into 6 pieces, largest 83 s / 5.8 GB; the chr13 rDNA piece cannot be cut below 4.78 M records: 110 s / 11.1 GB). Then the merge + the k11 all-vs-all of the representatives through the wrapper: PTR 45,488 reps (118 Mb), 12 shards in 3 calls (~9 + 9 + 4 min; last call 250 s / 1.2 GB); testis 17,541 reps (17 Mb), 22 s / 0.26 GB. Larger samples (OR6737, PPY, KB3781, A119b) est. from the dry-run: 0.5-6 h each |
| families | `--stage families`, repeat until exit 0 | THE default family definition, needed by figs 4–8 (with its copy table: a `mcl_families` that writes it). NEVER RUN GENOME-WIDE (est.): the all-vs-all of the de novo locus bodies, locus FASTA A119b 1,980 MB, OR6737 1,010, PPY 1,004, KB3781 877, PTR 791, testis 324 MB; chr2 alone (163 MB) took 416 s, time roughly quadratic per contig (dry-run range 14 min-17 h per sample). The wrapper builds one asm20 index per sample first (2 Gb target: est. ~5 min / 5-6 GB; if it does not fit one call the stage fails as "no progress": raise `runs_mm2_budget_s` for that call) |
| families_primary | `--stage families_primary`, repeat | OPTIONAL (only the supplementary fig. 6s-seeding reads it; no main figure): the same on the primaries-only assembly; about the same cost again (est. 4-31 h over six samples); queue it last or skip it |
| flag | `--stage flag` | est. 10-28 min / 13-19 GB in ONE process (scan 8-25 min + 1 min per index): over the 10-min rule; needs a per-contig scan in `missing_copy_flag` or the user's approval for one longer foreground run |
| assign | `--stage assign` | not read by any figure (figs 4-5 run their own copy assignment); est. 6 min (PTR)-5 h (A119b) in one process: skip unless needed |

`make.py runs --dry-run` prints the live queue with these estimates; a stage's `.time` replaces its estimate once it
has run.

### 2. Figure tables (`make.py data figN`; after the runs each figure needs)

| step | command | cost (est. from the phase-2 reports unless measured) |
|---|---|---|
| fig1 | `make.py data fig1 --set figs_budget_s=360`, repeat | 30-55 min in 5-9 calls: gffcompare of every method (human IsoSeq collapse 107 s / 5.1 GB measured), exact-chain read support per contig (human 11-23 min) |
| fig2 | `make.py data fig2 --set fig2_part=samples --set figs_budget_s=300`, repeat; then without `fig2_part` | 573 SQANTI3 calls, 12.7-25.2 h in total (< 3 GB per call); the Rustle-on-every-sample part first (3-5 h) |
| fig3 | `make.py data fig3 --set figs_budget_s=360`, repeat | 20-45 min: human per-transcript counts 17-40 min (67.6 M records, I/O-bound) |
| fig4, fig5 | `make.py data fig4 --set o2_samples=S`, repeat, one sample at a time; then `data fig5` | needs the sample's catalog. Simulation ~5 min (human), then one mapping part per call (~8 min; human 35-57 parts at 15-16 GB; ape and gorilla parts 19-21 GB, run alone), then one copy_assign shard per call. Experiment B (fig. 5c, A119b and OR6737): A119b 4-10 h, OR6737 1-2 h in family shards. The margin-rule supplement (fig5s) reuses these runs |
| fig6 | `make.py data fig6`, repeat | Compara (done, release 116); the human samples' spectrum T1/T2 all-vs-all through the wrapper (0.3-5 h, T3 never run genome-wide); needs `families` (human samples). Supplement fig6s-seeding: `primary_counts_gw` (one BAM pass per human sample over every gene record, est. 5-40 min), `families_primary` if run, the Liftoff tables of fig8; protein-homology families only with `--set fig6s_protein_homology=1` (not built by default) |
| fig7 | `make.py data fig7`, repeat | de novo families from the `families` stage; guided families per species through the wrapper (human 3-45 h, each ape 1-20 h, est.); Compara shared with fig6; Liftoff rows after fig8's self-lift and read support; protein-homology families only with `--set fig7_protein_homology=1` (supplement) |
| fig8 | `make.py data fig8 --set fig8_budget_s=540`, repeat | Liftoff self-lift per genome, one call per block of records (human 36 calls first); reads the runs' `assemble`, `families` (copy table), `flag` and, when present, the legacy `catalog` products (captions/fig8.md). Run it before `data fig7` / `data fig6` so their Liftoff rows exist |
| render | `make.py plot all`, `make.py check`, then `python3 figures/fig_*.py summary` for the caption numbers | seconds to minutes |

When a guided StringTie/FLAIR GTF is registered, `make.py data fig1|fig2|fig3` also builds the guided tables
(`fig2_part=guided` builds only those); nothing else changes.
