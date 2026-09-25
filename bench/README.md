# bench/ — the analysis scripts and per-topic reports

Regenerated 2026-09-24 after wave 7. bench/ holds **98 files**: **9 Python files** (6 command-line scripts, 3
libraries), 4 data tables, 80 reports, this README and 4 orphan READMEs. Wave 7 folded the **42** Python files that
waves 5 and 6 had kept (21 top-level, 12 in `layer_order/`, 9 in `soto/`) into these 9: 39 were removed, 6 are new and
3 were kept. Each replaced command was checked byte-identical against the old script on its recorded inputs; the
exceptions are listed under "Where the output is not the old output".

**Where the old files are.** At git tag `notebook-2026-09-24` (`git show notebook-2026-09-24:bench/<old path>`; to rerun a
provenance command exactly, `git checkout notebook-2026-09-24 -- bench/<old path>`), and copied to
`~/Desktop/Rustle_attic/2026-09-24/` (its `MANIFEST.tsv` gives the replacement and the reason for each). Scripts retired
by earlier waves are at `notebook-2026-09-23b` (wave 5) and `notebook-2026-09-23c` (wave 6); older ones sit under
`archive/` at `notebook-2026-09-20` / `-19` (see "Finding a script or output a document names" in
`docs/ACTIVE_WORKING_SET.md`). The reports (`*.md`) are the per-topic records that the ledger and register cite; they
are not pruned.

Every script runs from any directory (`python3 bench/score.py pairs --help`). Each module's docstring repeats its part
of the old-name table. Only the standard library loads at import time. pysam, numpy, scipy, sklearn and pyBigWig load
inside the subcommands that use them.

## Python files (9)

| file | kind | subcommands / main exports | what it is for |
|---|---|---|---|
| `lib.py` | library | `rc`, `ov`, `merge`, `UF`, `pairwise`, `bipartite_items` (item level), `bipartite_families` (family macro), `pair_scores`, `translate_refseq` / `translate_phased`, `longest_cds`, `gene_biotypes`, `gene_key_names`, `gtf_attr`, `load_compara`, `read_referee`, `soto_gene_family`, `families_on`, `paf_identity` / `paf_coverage`, `cigar_introns`, `sam_lines` (streamed samtools), `mcl` (calls `mcl_port`); the read-count rule `aligned_blocks` / `reads_with_block_in` / `reads_overlapping_span` / `spanning_genes` | Shared helpers for `score`, `sim`, `truth`, `guided_pipeline` and `layer_order`. The read-count functions are the ones `docs/THESIS_OBJECTIVES.md` rules 12 and 13 prescribe. |
| `score.py` | CLI, 14 subcommands | `pairs`, `spectrum`, `heldout`, `referee`, `edge-gap`, `rna-ceiling`, `members`, `adjudicated`, `protein`, `eichler`, `reads`, `bakeoff-calls`, `bakeoff-compare`, `locus-reads` | Scorers: family, pair and copy-assignment scoring against the project's truths. |
| `sim.py` | CLI, 5 subcommands | `chromosome`, `missing-copy`, `tandem`, `copies`, `excise`; library `simulate_reads`, `write_fastq`, `stable_seed` | Read simulators, with the truth in the read names. |
| `truth.py` | CLI, 4 subcommands | `nodes`, `adjudicated`, `protein`, `protein-referee`; library `excluded`, `pair_hsps`, `edges_from`, `load_genes`, `protein_edges`, `blastp_all_vs_all` | Truth builders. |
| `guided_pipeline.py` | CLI, flags only | `--workdir --gff --genes-gff-gz --genome --mmi --iqtree [--expected-units --reps --threads]`; library `gene_body_chains` (plus `merge`, `ov`, `rc`, `pairwise`, `bipartite`, re-exported from `lib`) | The guided O1 pipeline with the Addendum T fixes, on a leave-out of an annotated truth table. |
| `mcl_port.py` | library | `mcl(edges, inflation, prune, max_iter)` | Python MCL comparator. It is a thin shim over the bit-faithful Rust bin `mcl_port` (§6z3, r1047; `RUSTLE_MCL_PORT_BIN`). Since wave 7 its only importer is `lib.mcl` (its docstring's "eight bench scripts" is the wave-6 count). |
| `layer_order/lattice_common.py` | library | paths (`--root` / `LO_ROOT`), THE LEVEL TESTS (`tests`), `UF` / `components` / `truss3`, E1 catalogs (`catalog_context`, `catalog_keys`, `membership`), Soto mapping (`soto_load`, `soto_map_gene`), the EXPR read count (`count_reads`, `IntervalIndex`), clause-5 `c_tree`, scorers (`score_counts`, `score_lo`, `bip_jaccard`) | Library of the 2026-09-16 NPIP/TBC1D3 layer-order and nested edge-test lattice study. It keeps its name and old exports because the off-repo `LAT/corrections_pass2/*.py` import it by path. |
| `layer_order/npip_tbc1d3.py` | CLI, 11 stages + `all` | `expr-recount`, `corrected-tables`, `layer-order`, `lattice-edges`, `lattice-expr`, `lattice-levels`, `lattice-truth`, `lattice-filtration`, `lattice-check-c2`, `lattice-report`, `all [--with-check-c2]`; global `--root DIR` | Reproduces `LAYER_ORDER_NPIP_TBC1D3.md` and `NESTED_LATTICE_NPIP_TBC1D3.md` (L3 at 0.98, not the shipped 0.985). ⚠ The default root is the frozen results tree and every stage overwrites its outputs there, so point `--root` at a copy. It re-executes itself with `PYTHONHASHSEED=0`. Substrate: `docs/DATA.md`. |
| `soto/soto_replication.py` | CLI, 6 subcommands | `genesets`, `edges`, `cluster`, `dennislab`, `famcn`, `score` | Soto 2025 family replication. This is CONCORDANCE with Soto, not independent evidence (register T15/858). The §6ip headline is ARI 0.6959, 49.1% exact (median MAD). Recipe: `REPRODUCE.md` §5a; substrate: `docs/DATA.md`. |

## Data tables (4)

| file | what it is |
|---|---|
| `soto/soto_famCN_S1C.tsv` | Soto 2025 Table S1C (gene → family, famCN), the Soto truth. Read by `soto_replication.py`, `score.py heldout --soto` / `referee` / `edge-gap`, the layer-order Soto mapping, and `family_score --soto` (REPRODUCE). |
| `soto/soto_parCN_S1E.tsv` | Soto 2025 Table S1E, with dual CHM13 v1.0/v2.0 coordinates: the liftover anchors for `edges` and `famcn`. Wave 7 restored it from `cd37ccb0^` (wave 3 had archived it). |
| `soto/acro_extra_anchors.tsv` | Extra acrocentric liftover anchors (`--extra-anchors`, §6il). Restored from `cd37ccb0^` in wave 7. |
| `soto/shared_exons_2334_finalhuman.tsv` | The frozen §6ip `edges` output (4,192 shared-exon edges over 2,223 genes), so that `cluster` and `score` run from a clone. |

## Subcommands: what each was, and what it reproduces

| command | was | what it does | record |
|---|---|---|---|
| `score.py pairs` | `referee_band_score.py`; `identity_spectrum.py --catalog` | Pair-level family scoring. Recall is by annotated-mRNA identity band (`--bands paf:FILE`) or by Compara band; precision is over judgeable pairs. | register 1097-1098, 1101-1102 |
| `score.py spectrum` | `identity_spectrum.py` (tier mode) | Which edge tier (asm20 / sensitive k11 / protein) recovers which Compara paralogue pairs, by identity band, with the miss diagnosis. `--estimator de --coverage shorter_axis` is the builder's rule: a new arm, off by default. | §6zh, register 1096; REPRODUCE |
| `score.py heldout` | `heldout_family_score.py` | `mcl_families` clusters vs symbol-root truth (or `--soto`), per-family bipartite. `--exact-only` is the B4 fix, off by default. | register 966; `PREREG_heldout_families_2026-09-20` |
| `score.py referee` | `soto_vs_us_referee.py` | Ours vs Soto, both scored against the neutral protein referee. | register 934 (§6u5) |
| `score.py edge-gap` | `protein_edge_gap.py` | Does a protein-level edge close §6o8's no-edge gap? | `PREREG_protein_edges_2026-09-20`, `PROTEIN_EDGES_RESULT_2026-09-20` |
| `score.py rna-ceiling` | `rna_truth_from_protein.py` | A non-circular RNA-level truth from protein families, and the alignability ceiling it implies. | `DOMINANT_GAP_RESOLVED_2026-09-20` (§6t1) |
| `score.py members` | `member_completeness.py` | gffcompare `=` per universe gene. | register 1102 |
| `score.py adjudicated` | `adjudicated_truth.py score` | Catalogs vs the AK adjudicated two-annotation truth. | `PREREG_core_definition_2026-09-12` Addendum AK |
| `score.py protein` | `protein_families.py score` | Cross-annotation protein-family scoring by CDS overlap. | Addendum AN; ledger §6ko |
| `score.py eichler` | `eichler_compare.py` | The Eichler-style AS-margin assignment, compared with ours. | register 935-936; `EICHLER_COMPARISON_2026-09-21` |
| `score.py reads` | `copy_assign_read_truth.py score` | Copy assignment against per-read truth: OWN / PRIMARY / ANY readings, per divergence bin. | §6zf; `PREREG_o2_read_truth_2026-09-23`; REPRODUCE |
| `score.py bakeoff-calls` | `copy_assign_tool_bakeoff.py calls` | Per-molecule copy calls derived from a tool's GTF. | §6hz; `PREREG_tool_bakeoff_2026-09-08` |
| `score.py bakeoff-compare` | `copy_assign_tool_bakeoff.py compare` | Compares the per-tool calls on the hard (AS-tied) vs easy molecules. | PREREG hard_locus_bakeoff (5ca5c7e4) |
| `score.py locus-reads` | `locus_reads.py` | The correct read count at a locus: reads with an aligned block inside it, not reads that splice over it. | THESIS_OBJECTIVES rules 12-13 (§6cm) |
| `sim.py chromosome` | `ideal_chromosome_sim.py` | Ideal-chromosome read simulation (arms `ideal` / `trunc` / `rt`). | `PREREG_ideal_chromosome_sim_2026-09-21`; register 950 |
| `sim.py missing-copy` | `missing_copy_sim.py` | Missing-copy simulations with truth; modes `transcript` / `genomic` / `shuffled`. | `PREREG_o3_reference_bias_2026-09-23` (r1089), `PREREG_o3_rna_only_2026-09-23`; REPRODUCE positive control |
| `sim.py tandem` | `tandem_copy_sim.py` | Tandem / interleaved copy simulation; `--pipeline` also runs assembler, catalog and assignment. | §6zg, register 1095; REPRODUCE |
| `sim.py copies` | `copy_assign_read_truth.py sim` | Read-truth simulation from every copy of a catalog, mapped genome-wide. | §6zf; REPRODUCE |
| `sim.py excise` | `copy_assign_excision.py` | Removes copy X from a family, reruns `copy_assign`, and looks for the missing-copy signature. | PREREG adj/excise |
| `truth.py nodes` | `annotation_nodes.py` | Gene-level node tables from RefSeq / CAT / Ensembl. | Addenda AI/AJ/AK |
| `truth.py adjudicated` | `adjudicated_truth.py build` | Builds the AK two-annotation truth. | Addendum AK |
| `truth.py protein` | `protein_families.py build` | Protein-space multi-copy families. | Addendum AN; ledger §6ko |
| `truth.py protein-referee` | (new; inline before) | The protein referee as a `Gene Name` / `Family ID` table. | used inline by `referee`, `rna-ceiling`, `edge-gap` |
| `npip_tbc1d3.py expr-recount` … `lattice-report` | the 11 `layer_order/` scripts (table below) | One stage each; `all` runs both reproduce blocks in dependency order. | `LAYER_ORDER_NPIP_TBC1D3.md` §11, `NESTED_LATTICE_NPIP_TBC1D3.md` §11; registers 886, 887, 893, 1104 |
| `soto_replication.py genesets` | (new; the genesets were hand-made) | Derives the 2,334-gene and 1,793-gene genesets from S1C. | §6ip |
| `soto_replication.py edges` | `soto_replicate_from_sedef.py` | Steps 1-4: SEDEF CIGAR walk → shared-exon edge TSV. | ledger §6ie-§6ip |
| `soto_replication.py cluster` | `soto_cluster_from_shared.py` | Steps 5-6: components → famCN MAD split → families. | §6if-§6ip |
| `soto_replication.py dennislab` | `soto_cluster_dennislab_algorithm.py` | The Dennis-lab notebook algorithm (a parked arm). | §6if |
| `soto_replication.py famcn` | `famcn_from_wssd.py` | WSSD read-depth famCN at arbitrary v2.0 intervals. | §6ie, §6il, §6io |
| `soto_replication.py score` | `soto_score_against_truth.py` + `soto_bipartite_match_score.py` | ARI / exact / pair P-R-F1 (`--only pairs`) and bipartite family matching (`--only bipartite`) vs S1C. | §6ih-§6ip |

## Old name → new command (every file wave 7 removed)

A document that cites one of these paths or names resolves here. The old file itself is at tag `notebook-2026-09-24`.
The arguments are unchanged unless the "note" column says otherwise.

| old path (as cited) | now | note |
|---|---|---|
| `adjudicated_truth.translate` (function) | `lib.translate_phased` | 0-based segments, honours the first segment's phase; `protein_edge_gap.translate` is `lib.translate_refseq` (1-based, phase ignored) |
| `soto_vs_us_referee.protein_referee` (function) | `truth.protein_referee` / `truth.py protein-referee` | now writes `PREFIX.families.tsv` (`Gene Name`/`Family ID`) instead of returning the dict only |
| `bench/adjudicated_truth.py build …` | `bench/truth.py adjudicated …` | `build` had raised `TypeError` in `mcl_port` since 21d6c5c9 (B13); now fixed, and equal to the numpy-era build byte for byte |
| `bench/adjudicated_truth.py score …` | `bench/score.py adjudicated …` | |
| `bench/annotation_nodes.py KIND GFF CONTIGS OUT` | `bench/truth.py nodes KIND GFF CONTIGS OUT` | |
| `bench/copy_assign_excision.py FAM X OUT` | `bench/sim.py excise FAM X OUT [--bam --fasta --bin --threads]` | the hard-coded BAM / FASTA / BIN became flags, with the old values as defaults |
| `bench/copy_assign_read_truth.py sim TSV FA IDX OUT SEED` | `bench/sim.py copies TSV FA IDX OUT SEED [--threads 4]` | ⚠ stable seeds (B2): the reads differ from every earlier run |
| `CATALOG_TSV=CAT bench/copy_assign_read_truth.py score P O` | `bench/score.py reads --catalog CAT P O` | `CATALOG_TSV` is still honoured |
| `bench/copy_assign_tool_bakeoff.py calls …` | `bench/score.py bakeoff-calls …` | |
| `bench/copy_assign_tool_bakeoff.py compare …` | `bench/score.py bakeoff-compare …` | samtools output is streamed (B6) |
| `bench/eichler_compare.py …` | `bench/score.py eichler …` | |
| `bench/heldout_family_score.py …` | `bench/score.py heldout …` | `--exact-only` (the B4 fix) is off by default |
| `bench/ideal_chromosome_sim.py …` | `bench/sim.py chromosome …` | |
| `bench/identity_spectrum.py …` (tier mode) | `bench/score.py spectrum …` | `--estimator nm_bl --coverage query_over_min` are the defaults (the old rule) |
| `bench/identity_spectrum.py --gtf x --ref REF.gtf --fasta x --chrom C --compara CMP --out O --catalog COPIES [--universe U] [--referee R] [--gff-genes G]` | `bench/score.py pairs --members COPIES --genes REF.gtf\|G --chrom C --truth compara:CMP\|families:R [--universe U]` | the unused `--gtf/--fasta/--out` were dropped (B8) |
| `bench/locus_reads.py BAM CHROM START END` | `bench/score.py locus-reads BAM CHROM START END` | |
| `bench/locus_reads.py::reads_with_block_in`, `::spanning_genes`, `::reads_overlapping_span`, `::aligned_blocks` | `bench/lib.py::` the same names | docstrings kept verbatim |
| `bench/member_completeness.py ARM REF UNI GFF CHROM LABEL` | `bench/score.py members --gtf ARM --ref REF --universe UNI --chrom CHROM --label LABEL` | the unused `GENES_GFF` was dropped |
| `bench/missing_copy_sim.py MODE GTF FA CHROM K DIV OUT SEED` | `bench/sim.py missing-copy MODE GTF FA CHROM K DIV OUT SEED [--threads 4]` | ⚠ stable seeds (B2) |
| `bench/protein_edge_gap.py …` | `bench/score.py edge-gap …` | library: `longest_cds` → `lib.longest_cds`, `translate` → `lib.translate_refseq`, `protein_edges` → `truth.protein_edges` |
| `bench/protein_families.py build …` | `bench/truth.py protein …` | library: `excluded`, `pair_hsps`, `edges_from`, `load_genes` → `truth.*` |
| `bench/protein_families.py score …` | `bench/score.py protein …` | |
| `bench/referee_band_score.py CLUSTERS LOCI.gff3 LABEL` (in `/mnt/linuxdisk/tmp/gw22/sec/` = `$S`) | `bench/score.py pairs --members $S/CLUSTERS --genes $S/ref/NC_073244.2.genes.gff --chrom NC_073244.2 --truth families:$S/ref/NC_073244.2.tsv --expressed $S/ref/NC_073244.2.expressed.tsv --bands paf:$S/ref/referee_mrna.paf --label LABEL` | the unused `LOCI.gff3` was dropped |
| `bench/rna_truth_from_protein.py …` | `bench/score.py rna-ceiling …` | the referee alone: `bench/truth.py protein-referee` |
| `bench/sim_reads.py` (`from sim_reads import simulate_reads`) | `from sim import simulate_reads` (`write_fastq` too) | |
| `bench/soto_vs_us_referee.py …` | `bench/score.py referee … [--threads 4]` | the old file raised `NameError` on every run from 8db314c7 on (B1); its `pair_scores` / `gene_names` → `lib.pair_scores` / `lib.gene_key_names` |
| `bench/tandem_copy_sim.py …` | `bench/sim.py tandem …` | |
| `bench/layer_order/lo_expr_recount.py OUT [--ignore IDS] EXTRA…` | `bench/layer_order/npip_tbc1d3.py expr-recount OUT [--ignore IDS] EXTRA…` | every stage takes `--root DIR` (default: the frozen tree) |
| `bench/layer_order/lo_corrected_tables.py` | `npip_tbc1d3.py corrected-tables` | `catalog_keys` / `membership` / `CATALOGS` / `LIT` → `lattice_common` |
| `bench/layer_order/lo_analysis.py` | `npip_tbc1d3.py layer-order` | `c_tree` / `bip_jaccard` / `f1` / `hgnc_all` / `soto_ok` → `lattice_common`; `score` → `lattice_common.score_lo` |
| `bench/layer_order/lattice_edges.py` | `npip_tbc1d3.py lattice-edges` | its head (exec'd by check_c2) → `lattice_common.catalog_context` |
| `bench/layer_order/lattice_expr.py` | `npip_tbc1d3.py lattice-expr` | |
| `bench/layer_order/lattice_levels.py` | `npip_tbc1d3.py lattice-levels` | `shortest_path` → `lattice_common.shortest_path` |
| `bench/layer_order/lattice_truth.py` | `npip_tbc1d3.py lattice-truth` | `c2` / `score` → `lattice_common.c2` / `score_counts` |
| `bench/layer_order/lattice_filtration.py [--l1 c2_loose]` | `npip_tbc1d3.py lattice-filtration [--l1 c2_loose]` | |
| `bench/layer_order/lattice_check_c2.py` | `npip_tbc1d3.py lattice-check-c2` | it had not run since 2026-09-19 (F2, F2b); `denovo_shared_def.ExonIndex` → `lattice_common.IntervalIndex` |
| `bench/layer_order/lattice_report_tables.py` | `npip_tbc1d3.py lattice-report` | `bfs_path` → `lattice_common.shortest_path` |
| `bench/layer_order/soto_map.py` (`load`, `map_gene`, `FIELDS`) | `lattice_common.soto_load`, `soto_map_gene`, `SOTO_FIELDS` | |
| `bench/soto/soto_replicate_from_sedef.py …` | `bench/soto/soto_replication.py edges …` | `--s1e` is optional (default `bench/soto/soto_parCN_S1E.tsv`) |
| `bench/soto/soto_cluster_from_shared.py …` | `soto_replication.py cluster …` | default `--mad-statistic mean` kept |
| `bench/soto/soto_cluster_dennislab_algorithm.py …` | `soto_replication.py dennislab …` | default `--mad-statistic median` kept |
| `bench/soto/famcn_from_wssd.py …` | `soto_replication.py famcn …` | the `--s1e` default is relative to the module, not the cwd |
| `bench/soto/soto_score_against_truth.py …` | `soto_replication.py score --only pairs …` | `--truth` is optional (default S1C) |
| `bench/soto/soto_bipartite_match_score.py …` | `soto_replication.py score --only bipartite …` | a plain `score` prints both |
| `bench/soto/soto_attach_noncoding_members.py …` | `soto_replication.py cluster --full-geneset …` | its `main()` was already superseded in §6ii (same partition); `attach()` kept |
| `bench/soto/soto_replicate_clustering.py …` | not replaced | its minimap2 map-back path was superseded by the SEDEF path (§6ie); `load_exons()` kept |
| `bench/soto/rustlib.py` | not replaced | 0 importers since wave 5; cited only as provenance of the frozen E_r mirror (`docs/seeded_family_definition.md`, register 1044, `denovo_pipeline.rs` doc comments) |

**Earlier names of the same code.** Older documents cite some of these scripts by names they had before wave 7. Three
of those names never sat at a notebook tag; read them with `git show 8db314c7:bench/<name>`. A run older than the tag
used the code at its own commit, so check out that commit to reproduce it exactly.

| cited as | became | now |
|---|---|---|
| `bench/o2_read_truth_sim.py`, `bench/o2_read_truth_score.py` (tag `notebook-2026-09-23c`) | `o2_read_truth.py sim` / `score` (wave 6), then `copy_assign_read_truth.py` (publishing names, 2026-09-24) | `sim.py copies` / `score.py reads` |
| `bench/o2_read_truth.py` (commit 8db314c7) | `copy_assign_read_truth.py` | `sim.py copies` / `score.py reads` |
| `bench/tool_bakeoff.py` (`notebook-2026-09-23c`) | `o2_tool_bakeoff.py calls` (wave 6), then `copy_assign_tool_bakeoff.py calls` | `score.py bakeoff-calls` |
| `bench/hard_locus_bakeoff.py` (`notebook-2026-09-23c`) | `o2_tool_bakeoff.py compare`, then `copy_assign_tool_bakeoff.py compare` | `score.py bakeoff-compare` |
| `bench/o2_tool_bakeoff.py` (commit 8db314c7) | `copy_assign_tool_bakeoff.py` | `score.py bakeoff-calls` / `bakeoff-compare` |
| `bench/o2_excision.py` (`notebook-2026-09-23c`) | `copy_assign_excision.py` (a pure rename) | `sim.py excise` |
| `bench/o3_sim_copies.py` (commit 8db314c7) | `missing_copy_sim.py` | `sim.py missing-copy` |
| `/mnt/linuxdisk/tmp/gw22/o3/simB.py GTF FA K DIV OUT SEED` (out of repo; chr20 hard-coded) | — | `sim.py missing-copy genomic GTF FA chr20 K DIV OUT SEED` |
| `bench/neighbourhood_jaccard.py` (`notebook-2026-09-23c`) | retired in wave 6; its `gene_names` was inlined into `soto_vs_us_referee.py` | `lib.gene_key_names` (the analysis itself is only at the tag) |
| `bench/guided_min.py` (`notebook-2026-09-23c`) | retired in wave 6; its `load_genes` was inlined into `annotation_nodes.py` | `truth.nodes_load_genes` (the pipeline itself is only at the tag) |

## Where the output is not the old output

- **`sim.py copies` and `sim.py missing-copy` draw different reads from every earlier run (B2).** The old scripts
  seeded their per-copy and per-gene RNGs with Python's `hash(str)`, which is salted per process, so no earlier run
  could be reproduced read for read (PYTHONHASHSEED 0 and 1 gave different reads). They now seed with
  `sim.stable_seed()`. Two runs give the same output, and with the old seeding restored under a fixed
  PYTHONHASHSEED the new code equals the old code byte for byte. The numbers measured on the old runs (REPRODUCE's
  O2 read-truth "Expected" line, the missing-copy 40/40 positive control) must be re-measured, not compared read for read.
- **`npip_tbc1d3.py` pins `PYTHONHASHSEED=0`.** Seven layer-order outputs break ties by set iteration order (N1). The
  pinned rerun equals the old code under seed 0; the frozen 2026-09-16 files carry an unrecorded seed's tie order (the
  same rows and numbers, some in a different order).
- **`score.py pairs`** breaks two ties that the old scripts left to string-hash order by a fixed rule. On every
  recorded input the old output was the same under PYTHONHASHSEED 0 to 3, and the new output equals it.
- **`score.py spectrum`**: the mmseqs `t3.m8` is identical as a set, not line for line (multi-threaded line order).
- **Scoring disagreements between the old scripts were not unified.** Each subcommand keeps its own script's rule,
  under a distinct name (`bipartite_items` vs `bipartite_families`, `translate_refseq` vs `translate_phased`,
  `protein_edges` vs `edges_from`, the per-script member→gene rules). A better rule exists only behind a flag
  (`spectrum --estimator/--coverage`, `heldout --exact-only`).

## Orphan READMEs (4)

They describe outputs, or exporter scripts, that earlier waves archived; the files they describe are not in the tree.

| file | title |
|---|---|
| `mechanism/demo/README.md` | V4c — identity-gradient frontier demo |
| `soto/dna_graphs/README.md` | Soto DNA variation graphs -- all families |
| `soto/family_definition_graphs/README.md` | Rustle family-definition graphs |
| `soto/member_similarity_graphs/README.md` | Compact family-member similarity graphs |

## Reports (80)

| report | title | last change |
|---|---|---|
| `ASSEMBLER_WIDENING.md` | Assembler read-isoform widening + single-exon strand — vs `docs/PREREG_assembler_widening_2026-09-18 | 2026-09-18 |
| `ASSEMBLE_ONLY_MODE.md` | `--assemble-only`: the assembler product, with none of the all-vs-all work | 2026-09-19 |
| `ASSEMBLY_POLISH.md` | Assembly polish: matching StringTie in `--assemble-only` mode (§6p8, 2026-09-19) | 2026-09-23 |
| `CHIMERA_POLICY.md` | `--chimera-policy` — measured against `docs/PREREG_chimera_policy_2026-09-18.md` (md5 57f578ff0c4e2e | 2026-09-18 |
| `CHR16_JUNCTION_MAJORITY_ARM.md` | The chr16 arm `build_spliced_seq_with` demands before `RUSTLE_JUNCTION_MAJORITY` can be a default | 2026-09-18 |
| `CHR20_ASSEMBLER_COMPARISON.md` | Chr20 assembler comparison: ours vs StringTie vs FLAIR (gffcompare + SQANTI3) | 2026-09-16 |
| `CLUSTERING_OPERATOR_BAKEOFF.md` | Clustering-operator bakeoff on the L2 copy graph — DESCRIPTIVE, development families only | 2026-09-18 |
| `COPY_ASSIGNMENT_AND_GATE.md` | Copy Assignment And Gate (consolidated) | 2026-09-01 |
| `CROSS_SPECIES_NPIP_CONJUNCT.md` | Agent 3 of 4 — the certificates | 2026-09-17 |
| `DEFINITIONS_FORMAL.md` | Five Concepts, Four Baselines: Paralog, Segmental Duplication, Multi-Copy Gene Family, Expansion, an | 2026-08-14 |
| `DENOVO_PIPELINE.md` | Denovo Pipeline (consolidated) | 2026-09-01 |
| `FALSE_NEGATIVES.md` | False negatives: what the pipeline misses, and why | 2026-08-26 |
| `FAMILY_CERTIFICATES_NPIP_TBC1D3.md` | Family certificates on NPIP and TBC1D3: where each family is an exact connected component, and with  | 2026-09-17 |
| `FAMILY_DEF.md` | Family Definition (consolidated) | 2026-09-01 |
| `FAMILY_LEVELS_AND_RELATED.md` | Family Levels (RNA/DNA/Protein) & Related Methods (consolidated) | 2026-09-01 |
| `FEX_SWEEP_LORO.md` | Raising `--min-shared-exon-frac` from 0.30 to 0.60 — leave-one-region-out validated | 2026-09-18 |
| `GAP_CLOSED_FRACTION.md` | Gap-closed fraction: read-isoform widening and evidence-based admission floors, measured relative to | 2026-09-17 |
| `GATE_CENSUS_NPIP.md` | Where pass-1 skeletons die, and why the "ceiling" was not a depth limit | 2026-09-18 |
| `GENOME_WIDE_BAKEOFF_2026-09-22.md` | Genome-wide `--assemble-only` bakeoff: ours vs StringTie vs FLAIR | 2026-09-23 |
| `GFFCOMPARE_CHR20_2026_09_19.md` | gffcompare on chr20: ours vs StringTie vs FLAIR, with today's assembly settings | 2026-09-19 |
| `GTF_SECONDARY_POOL.md` | `RUSTLE_GTF_SECONDARY` — admitting secondary alignments into the `--gtf` read pool | 2026-09-18 |
| `IDEAL_WIDENING_K5.md` | Read-isoform widening at k = 5 on the idealized simulated substrate (agent 1 of 2) | 2026-09-17 |
| `ISOSEQ_FLAIR_MECHANISMS.md` | isoseq collapse and FLAIR: as comparison arms, and as mechanisms (§6q6, 2026-09-19) | 2026-09-19 |
| `JUNCTION_AND_READTHROUGH_RULES.md` | Adopting the two proposed definition changes: what each form costs, and why none of them is adoptabl | 2026-09-17 |
| `JUNCTION_MAJORITY_CHR16.md` | Splitting the `build_spliced_seq` bucket — and the chr16 arm the code asked for | 2026-09-18 |
| `LAB_DATASET_BAKEOFF.md` | Against the lab's own isoseq / StringTie / FLAIR runs (§6q7, 2026-09-19) | 2026-09-19 |
| `LATTICE_RULE_STRENGTHENERS.md` | Strengthening the DNA levels of the nested edge-test lattice: five strengtheners, swept, on NPIP and | 2026-09-17 |
| `LAYER_ORDER_NPIP_TBC1D3.md` | Layer order for NPIP and TBC1D3 (human T2T-CHM13): protein P, DNA catalog E1 ("D"), subfamily clades | 2026-09-24 |
| `LOCUS_ASSEMBLY_NPIP.md` | Locus assembly at NPIP — measured against `docs/PREREG_locus_assembly_2026-09-18.md` (md5 02237f9d6d | 2026-09-18 |
| `LOCUS_WIDTH_GAP.md` | Locus formation (node width and composition) vs annotated gene records — agent 1 of 2 | 2026-09-17 |
| `MERGED_LOCI_LAYER.md` | The MERGED-LOCUS layer: real fusions recorded as dual membership, outside the partition | 2026-09-19 |
| `NESTED_LATTICE_NPIP_TBC1D3.md` | Nested edge-test lattice on NPIP and TBC1D3 (human T2T-CHM13) | 2026-09-24 |
| `NODE_CUT_RULE.md` | The NODE CUT rule — measured against `docs/PREREG_node_cut_2026-09-18.md` (md5 2af3393070d6c2ded8db3 | 2026-09-18 |
| `NO_READTHROUGH_COUNTERFACTUAL.md` | A world with no readthroughs: it does not rescue the family definition, and on real data it cannot b | 2026-09-18 |
| `NPIP_DISHUCK_TRUTH_CHECK.md` | NPIP subfamily truth checked against Dishuck et al. 2025 | 2026-09-16 |
| `NPIP_IDEAL_EXPRESSION.md` | Agent 2 of 3 — certificates on the idealized synthetic substrate, plus the RNA-structure conjunct | 2026-09-17 |
| `NPIP_MEMBERSHIP_RESCORE.md` | NPIP membership re-run of the shared-junction conjunct t_J (k = 2, δ = 0) | 2026-09-17 |
| `NPIP_NODE_LADDER.md` | NPIP node ladder: annotated DNA nodes → de novo shared-definition nodes (human CHM13, 2026-09-17) | 2026-09-17 |
| `NPIP_NODE_RESCORE.md` | NPIP node arms re-scored: SCORER v2 (one-to-one copy mapping) + GATE v2 (agreement with guided A1) | 2026-09-17 |
| `NPIP_PER_MEMBER_METRICS.md` | NPIP per-member precision / recall / bipartite matching — human CHM13 (2026-09-16) | 2026-09-16 |
| `NPIP_PER_MEMBER_WHOLECHR16.md` | NPIP per-member metrics on whole-chr16 catalogs, with single-exon mapping variants — human CHM13 (20 | 2026-09-17 |
| `NPIP_READ_GROUP_NODES.md` | NPIP read-group nodes vs shipped de novo nodes (human CHM13, 2026-09-17) | 2026-09-17 |
| `NPIP_SIM_CEILING.md` | NPIP algorithmic-ceiling simulation — vs `docs/PREREG_npip_sim_2026-09-18.md` (md5 17be6b031092b79b3 | 2026-09-18 |
| `PERFORMANCE_AND_IO.md` | Performance And Io (consolidated) | 2026-07-09 |
| `PRECISION_LEVERS_CHR20.md` | Getting more transcripts, or better precision — measured on chr20; and what `-R`/`-Q` actually do | 2026-09-19 |
| `READ_ISOFORM_LOCUS.md` | Read-isoform locus and node-admission floor, pre-registered — agent 1 of 2 | 2026-09-17 |
| `RNA_LEVEL_DEFINITION.md` | An RNA-level family definition scored against RNA-DERIVED truth | 2026-09-18 |
| `SOTO_AS_A_REFINEMENT.md` | Are Soto's families a REFINEMENT of ours? (each Soto family inside one of ours) | 2026-09-18 |
| `SOTO_BENCHMARK.md` | Soto benchmark — RNA copy recovery on the human segmental-duplication catalog | 2026-07-16 |
| `SOTO_FAMILY_LOCUS_FIDELITY_2026-09-22.md` | How faithfully does each Soto family's genes get reproduced by our de novo loci? | 2026-09-22 |
| `SQANTI3_POLISH.md` | SQANTI3 on the polished `--assemble-only` output (§6q5, 2026-09-19) | 2026-09-19 |
| `STRAND_ASYMMETRY_AND_REALISTIC_SIM.md` | Strand fix for the shared-definition node set: four variants on two substrates (agent 1 of 2) | 2026-09-18 |
| `STRAND_V4_SCORER_V3.md` | SCORER v3 + genome-wide family panel: V1 / V3 / V4 on S-IDEAL, U00, U20, U40 (agent 1 of 2) | 2026-09-18 |
| `TBC1D3_GUITART_TRUTH_CORRECTION.md` | TBC1D3 subfamily truth corrected from Guitart 2024 Fig 6B/6C | 2026-09-16 |
| `THEORY.md` | Theory (consolidated) | 2026-09-01 |
| `TPM_AND_LOCUS_BED.md` | TPM in the GTF, and a locus BED that shows how close our loci are to the annotation (§6r5, 2026-09-1 | 2026-09-19 |
| `V4_GORILLA_TIEBREAK.md` | V4 on real gorilla data — the tiebreak (agent 1 of 2) | 2026-09-18 |
| `V4_HELDOUT_NOISE_CALIBRATED.md` | V4 on held-back substrates, with the leader rule's noise floor measured and the family clause calibr | 2026-09-18 |
| `V5_RETIRE_PLACEHOLDER.md` | Verification (independent recompute) — agent 2 of 2 | 2026-09-18 |
| `VALIDATION_AND_STATUS.md` | Validation, Reviews & Objective Status (consolidated) | 2026-09-01 |
| `mechanism/consolidation_divergences.md` | Consolidation divergences (deliverable B) | 2026-07-21 |
| `o1_expanded_family_audit/README.md` | Expanded O1 known-family purity audit | 2026-08-19 |
| `o1_fresh_emission_validation/README.md` | O1 fresh-emission validation | 2026-08-19 |
| `o1_gene_family_audit/README.md` | O1 typed gene-family audit | 2026-08-19 |
| `o1_golga2_subfamily_audit.md` | GOLGA2 versus GOLGA6/8: family or false merge? | 2026-08-26 |
| `o1_outgroup_rooting_poc/README.md` | O1 single-outgroup rooting proof of concept | 2026-08-19 |
| `o1_provenance_witness_prototype/README.md` | O1 duplication-provenance witness prototype | 2026-08-19 |
| `soto/MEETING_2026_08_04.md` | O1 status: what is established, what was refuted | 2026-08-04 |
| `soto/bam_tie_signals.md` | BAM signals for multimapper ties: what AS can't do, what `de` fixes, and what's still unused | 2026-08-14 |
| `soto/candidate18_classification.md` | SD classification of the 18 "no-Soto-homology" RNA candidates | 2026-08-22 |
| `soto/coordinate_version_check.md` | Coordinate-version check: Soto BED vs the A119b BAM (CHM13 v1 vs v2) | 2026-07-23 |
| `soto/dna_ceiling_rna_subset.md` | The T2T DNA graph is the 100% ceiling; RNA is a faithful subset | 2026-07-23 |
| `soto/dna_vs_rna_mode.md` | DNA vs RNA, one engine: reproducing Soto's families from the genome vs from transcripts | 2026-07-27 |
| `soto/exon_homogenized_nonexon_signal.md` | Exon-homogenized copies carry recoverable non-exon signal — PPIAL4 (ID_431) | 2026-07-21 |
| `soto/gate_robustness.md` | Soto gate robustness — do the readthrough / mis-chain gates distort detection? (2026-07-23) | 2026-09-01 |
| `soto/merge_quality_analysis.md` | Over- and under-merging in the DNA family partition — diagnosis and levers | 2026-09-01 |
| `soto/nonexon_rescue.md` | Non-exon-signal rescue POC — 24 K=0-bearing Soto families | 2026-07-21 |
| `soto/precision_recall_audit.md` | Soto precision/recall audit — current binary (2026-07-23) | 2026-07-23 |
| `soto/soto_sensitivity_precision.md` | Soto per-family sensitivity / precision — updated 2026-07-21 | 2026-07-21 |
| `soto/tied_seed_eval.md` | Tied-secondary seeding — Phase-1 benchmark result (2026-07-22) | 2026-07-23 |
