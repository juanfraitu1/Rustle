# Script inventory (2026-11-02)

Every `*.py *.sh *.R *.pl` outside `target/` (no .R/.pl exist): **687 scripts**. Generated read-only; nothing moved. Status rules: PIPELINE = reached from a shipped tool; PORTED = Rust port named in src/docs; KEEP-PY = final home is Python (scoring/report/orchestration/tests); TO-PORT = pipeline hot path still in Python; TO-MERGE = near-duplicate (survivor named); ARCHIVE = one-off study script (doc/register cited); UNKNOWN = undecided.

Key fact: `tools/rustle_pipeline.sh` invokes **no Python at all** (only `target/release/{utilities,as_table,copy_assign,mcl_families,gw_family_catalog,candidate_copies,missing_copy_flag}`); `.py` names in it are comments. The only Python on a shipped path is `bench/loop_home.py`, via the opt-in `tools/rustle_reassemble.sh`.

## Rust binaries (Cargo.toml `[[bin]]`; there is no src/bin/)

| bin | source |
|---|---|
| `mcl_families` | src/mcl_families.rs |
| `copy_assign` | src/copy_assign_bin.rs |
| `gw_family_catalog` | src/gw_family_catalog.rs |
| `readthrough_filter` | src/readthrough_filter.rs |
| `family_score` | src/family_score.rs |
| `missing_copy_flag` | src/missing_copy_flag.rs |
| `candidate_copies` | src/candidate_copies.rs |
| `as_table` | src/as_table.rs |
| `copy_phasing` | src/copy_phasing_bin.rs |
| `utilities` | src/utilities.rs |
| `family_gfa` | src/family_gfa.rs |
| `semisec` | src/semisec.rs |
| `semifam` | src/semifam.rs |
| `slimbam` | src/slimbam.rs |

Driver uses: utilities (sd-scan, candidate-augment, window-depth), as_table, copy_assign, mcl_families, gw_family_catalog, candidate_copies, missing_copy_flag. Not called by the driver: family_score, readthrough_filter, copy_phasing, family_gfa, semisec, semifam, slimbam.

## Non-bench scripts (one row each)

| file | status | note |
|---|---|---|
| figures/_liftoff.py | KEEP-PY | thesis figures/report (thesis repo; not in species-agnostic package) |
| figures/_lrc.py | KEEP-PY | thesis figures/report (thesis repo; not in species-agnostic package) |
| figures/_o1.py | KEEP-PY | thesis figures/report (thesis repo; not in species-agnostic package) |
| figures/_o1_recovery.py | KEEP-PY | thesis figures/report (thesis repo; not in species-agnostic package) |
| figures/_o2.py | KEEP-PY | thesis figures/report (thesis repo; not in species-agnostic package) |
| figures/_sqanti.py | KEEP-PY | thesis figures/report (thesis repo; not in species-agnostic package) |
| figures/assembly.py | KEEP-PY | thesis figures/report (thesis repo; not in species-agnostic package) |
| figures/fig_assign_accuracy.py | KEEP-PY | thesis figures/report (thesis repo; not in species-agnostic package) |
| figures/fig_assignability.py | KEEP-PY | thesis figures/report (thesis repo; not in species-agnostic package) |
| figures/fig_family_recovery.py | KEEP-PY | thesis figures/report (thesis repo; not in species-agnostic package) |
| figures/fig_family_spectrum.py | KEEP-PY | thesis figures/report (thesis repo; not in species-agnostic package) |
| figures/fig_intron_chain.py | KEEP-PY | thesis figures/report (thesis repo; not in species-agnostic package) |
| figures/fig_loci.py | KEEP-PY | thesis figures/report (thesis repo; not in species-agnostic package) |
| figures/fig_loop.py | KEEP-PY | thesis figures/report (thesis repo; not in species-agnostic package) |
| figures/fig_protein_supp.py | KEEP-PY | thesis figures/report (thesis repo; not in species-agnostic package) |
| figures/fig_secondary.py | KEEP-PY | thesis figures/report (thesis repo; not in species-agnostic package) |
| figures/fig_sqanti.py | KEEP-PY | thesis figures/report (thesis repo; not in species-agnostic package) |
| figures/figlib.py | KEEP-PY | thesis figures/report (thesis repo; not in species-agnostic package) |
| figures/make.py | KEEP-PY | thesis figures/report (thesis repo; not in species-agnostic package) |
| figures/samples.py | KEEP-PY | thesis figures/report (thesis repo; not in species-agnostic package) |
| figures/test_lrc.py | KEEP-PY | thesis figures/report (thesis repo; not in species-agnostic package) |
| tests/fixtures/candidate_copies/driver/run_e2e.sh | KEEP-PY | cargo test fixture generator |
| tests/fixtures/candidate_copies/make_fixture.py | KEEP-PY | cargo test fixture generator |
| tests/fixtures/o2_origin_resolve/make_bam.py | TO-MERGE | -> shared tests/fixtures/bamlib.py (lcg/filler/pattern/rec duplicated) |
| tests/fixtures/same_chrom_supplement/make_bam.py | KEEP-PY | cargo test fixture generator |
| tests/fixtures/xfam_conflict/make_bam.py | TO-MERGE | -> shared tests/fixtures/bamlib.py (lcg/filler/pattern/rec duplicated) |
| tools/audit_cleanup_candidates.py | KEEP-PY | dev-only: cleanup audit (not packaged) |
| tools/docs_index.py | KEEP-PY | dev-only: docs/INDEX.md generator (not packaged) |
| tools/identity_check.sh | KEEP-PY | dev-only: golden identity harness (shell; not packaged) |
| tools/minimap2_gate/build.sh | PIPELINE | builds `minimap2-gate`, the exon-gated all-vs-all `mcl_families --from-gtf` uses by default (O1_FAST_2026-10-21) |
| tools/mm2_shard.sh | PIPELINE | RUSTLE_MINIMAP2 wrapper for mcl_families / gw_family_catalog all-vs-all |
| tools/protein_attach.py | KEEP-PY | manual protein-homology attach (r1480 NULL; driver stage named in docstring no longer exists) |
| tools/rlock.sh | ARCHIVE | dev-machine job lock (WSL2 5-core box, /mnt/linuxdisk); Fritz uses Slurm |
| tools/rustle_pipeline.sh | PIPELINE | the driver (stages sd, assemble, families, catalog, candidates, assign, flag); calls Rust bins only |
| tools/rustle_reassemble.sh | PIPELINE | closed-loop pass 2 (opt-in; union certificate FAILED r1481) |
| tools/test_protein_attach.py | KEEP-PY | tests of protein_attach.py |
| tools/candidate_augment.py.legacy (not .py) | PORTED | `utilities candidate-augment` (driver `candidates` stage); delete with archive |

## bench/ top level and individually decided bench files

| file | status | note |
|---|---|---|
| bench/assembly_ends/end_score.py | KEEP-PY | assembly 5'/3' end scorer (PREREG_assembly_ends_2026-11-01) |
| bench/denovo_eval/famtruth.py | KEEP-PY | truth loader of pairscore/seqscore (human CHM13 + gorilla only) |
| bench/denovo_eval/pairscore.py | KEEP-PY | family scoring rule: pair precision/recall (+ --noncoding) |
| bench/denovo_eval/semifam.py | PORTED | `semifam` bin (src/semifam.rs, byte-identical; r1472) |
| bench/denovo_eval/seqscore.py | KEEP-PY | what genes really share (sequence tiers) |
| bench/f1_bridge.py | PORTED | `copy_assign --bridge-regroup` (family::bridge_regroup) |
| bench/family_container.py | PORTED | `mcl_families --emit-container` (rustle::family_container) |
| bench/guided_eval/cdsseed.py | TO-MERGE | -> famtruth.py (only load_genes("hg002") + HG002_XREF are used by pair/seqscore) |
| bench/guided_eval/guided3.py | TO-MERGE | -> famtruth.py (only load_genes("hg002") + HG002_XREF are used by pair/seqscore) |
| bench/hapcn_mat_pat_robust.py | TO-MERGE | -> hapcn_mat_pat.py (exec()s it with one edit) |
| bench/lib.py | KEEP-PY | shared helpers (sequence, annotation, truth, reads) |
| bench/loop_home.py | PIPELINE | called by rustle_reassemble.sh `home` step (only Python on any shipped path) |
| bench/mcl_port.py | PORTED | `utilities mcl-port` (bit-faithful; r1047); file is now a shim |
| bench/o2_p18_features_cond.py | TO-MERGE | -> o2_p18_features.py (exec()s it with one edit) |
| bench/o2_p20_sim19.py | TO-MERGE | -> o2_p20_sim20.py (16-line diff) |
| bench/o3_graph/o3port.py | PORTED | `copy_phasing` bin (r1402, 122 golden cases) -- but psv_vcf.py still imports it |
| bench/o3_graph/o3port_golden.sh | PORTED | golden/parity harness of the o3port -> copy_phasing port |
| bench/o3_graph/o3port_p4.py | PORTED | golden/parity harness of the o3port -> copy_phasing port |
| bench/o3_graph/o3port_p4.sh | PORTED | golden/parity harness of the o3port -> copy_phasing port |
| bench/o3_graph/psv_vcf.py | KEEP-PY | PSV VCF report (r1484); imports o3port.py |
| bench/o3_graph/test_o3port.py | PORTED | golden/parity harness of the o3port -> copy_phasing port |
| bench/o3_graph/test_psv_vcf.py | KEEP-PY | tests of psv_vcf.py |
| bench/o3_structure/proto.py | PORTED | `candidate_copies` (family::structure, PREREG_o3_structure_projection) |
| bench/o3_unexplained/psvclust.py | PORTED | `copy_assign --psv-clusters` (family::psv_clusters, md5-pinned v2) |
| bench/rna_allele/link_test.py | PORTED | `candidate_copies` link/merge rules (family::detect) |
| bench/rna_allele/merge_test.py | PORTED | `candidate_copies` link/merge rules (family::detect) |
| bench/rna_allele/yag_test.py | PORTED | `candidate_copies` link/merge rules (family::detect) |
| bench/score.py | KEEP-PY | shared scorers (family/pair/copy) used by figures/_o1,_o2 |
| bench/sdcn_noise_ex.py | TO-MERGE | -> sdcn_noise.py (4-line diff) |
| bench/semisec.py | PORTED | `semisec` bin (src/semisec.rs) |
| bench/sim.py | KEEP-PY | read simulators with truth in names (figures/_o2) |
| bench/stgaps/polya_mono.py | TO-MERGE | -> polya_mono_genomic.py (genomic form superseded read form) |
| bench/stgaps/polya_mono_genomic.py | PORTED | `copy_assign --polish-mono-primed` (default on 2026-10-26) |
| bench/test_family_container.py | PORTED | `mcl_families --emit-container` (rustle::family_container) |
| bench/test_loop_home.py | KEEP-PY | tests of loop_home.py |
| bench/test_units_from_annotation.py | KEEP-PY | tests of units_from_annotation.py |
| bench/truth.py | KEEP-PY | truth builders (figures/_o1_recovery, score.py) |
| bench/units_from_annotation.py | KEEP-PY | GUIDED detector -> copy_assign --bridge-units-list / RUSTLE_READTHROUGH_JUNCTIONS=list: |
| bench/famsim/* (12 files) | UNKNOWN | simulator package; design doc only (archive/2026-10/FAMSIM_DESIGN_2026-10-02.md), no register row/run; keep as test simulator or archive? |
| bench/*.py,*.sh other 38 files | ARCHIVE | one-off studies cited in LEDGER_O2 (o2_p16..p22, as_tie_p15, mm2_primary_rank, o2_tie_window_real), LEDGER_O3 (exhap*, sdcn*, bpclust, exceil, wgschk*, triochk, mmdemo_flank, o3_posctl.sh), LEDGER_O1/o1_ledger (guided_pipeline, hapcn*, dna_sd_atoms, semitruth, copy_support, npip_chains); registers 1105,1235,1239,1403,1433,1434... |

## bench/ subdirectories (ARCHIVE unless listed above)

| dir | n | PIPELINE | PORTED | KEEP-PY | TO-PORT | TO-MERGE | ARCHIVE | UNKNOWN | evidence (docs/) | register rows |
|---|---|---|---|---|---|---|---|---|---|---|
| annotation | 6 |  |  |  |  |  | 6 |  | PREREG_sqanti_single_copy_2026-11-01 | 1248 |
| asj_o2 | 5 |  |  |  |  |  | 5 |  | LEDGER_O2 (archive/2026-11) | - |
| assembly_ends | 1 |  |  | 1 |  |  |  |  | - | - |
| default_rescore | 4 |  |  |  |  |  | 4 |  | archive/2026-10/merged/DEFAULT_RESCORE_NPIP_2026-10-06 | - |
| denovo_eval | 40 |  | 1 | 3 |  |  | 36 |  | LEDGER_O1; DATA.md | 1435,1436,1437,1440,1453,1455… |
| entangled | 15 |  |  |  |  |  | 15 |  | LEDGER_O1/O2 (archive/2026-11) | 765,1257,1259,1260 |
| famsim | 12 |  |  |  |  |  |  | 12 | archive/2026-10/FAMSIM_DESIGN_2026-10-02 | - |
| g1 | 2 |  |  |  |  |  | 2 |  | archive/2026-09/PREREG_liftoff_loci_2026-09-25 | - |
| guided_eval | 6 |  |  |  |  | 2 | 4 |  | LEDGER_O1 GUIDED-S* | 1442,1444,1445,1446,1448,1449… |
| hierarchy | 15 |  |  |  |  |  | 15 |  | superpowers/plans/2026-10-14-soto-subset-cn; LEDGER_O1 | - |
| ideal_expression | 15 |  |  |  |  |  | 15 |  | LEDGER_O1 | 1253 |
| layer_order | 2 |  |  |  |  |  | 2 |  | o1_ledger | 1104,1105 |
| liftoff_pairs | 4 |  |  |  |  |  | 4 |  | archive/2026-09/PREREG_liftoff_loci_2026-09-25 | - |
| lrpap1_seeds | 3 |  |  |  |  |  | 3 |  | archive/2026-10/merged/PREREG_lrpap1_seed_reach_2026-10-11 (docstrings) | - |
| mechanism | 5 |  |  |  |  |  | 5 |  | archive/2026-09/PREREG_readthrough_v2_2026-09-26 | - |
| npip_read_pool | 6 |  |  |  |  |  | 6 |  | archive/2026-10/NPIP_READ_POOL_2026-10-01 | - |
| o1_hier | 15 |  |  |  |  |  | 15 |  | LEDGER_O1 | 1276,1277 |
| o1_reassign | 4 |  |  |  |  |  | 4 |  | archive/2026-10/merged/PREREG_post_mcl_reassign_2026-10-15 | - |
| o1_sd | 4 |  |  |  |  |  | 4 |  | LEDGER_O1 | - |
| o2_anchor | 3 |  |  |  |  |  | 3 |  | archive/2026-10/merged/PREREG_o2_anchor_2026-10-13 | - |
| o2_default_roster | 2 |  |  |  |  |  | 2 |  | archive/2026-10/merged/O2_DEFAULT_ROSTER_2026-10-06 | - |
| o2_tie | 2 |  |  |  |  |  | 2 |  | PREREG_o2_exact_tie_2026-10-13 | - |
| o2_truth | 4 |  |  |  |  |  | 4 |  | LEDGER_O2 | - |
| o2levers | 1 |  |  |  |  |  | 1 |  | LEDGER_O2; ACTIVE_WORKING_SET | 1481 |
| o3_fewread | 2 |  |  |  |  |  | 2 |  | archive/2026-10/merged/PREREG_o3_fewread_flags_2026-10-14 | - |
| o3_graph | 108 |  | 5 | 2 |  |  | 101 |  | LEDGER_O2/O3 | 1304,1307,1314,1318,1320,1323… |
| o3_maternal | 23 |  |  |  |  |  | 23 |  | LEDGER_O3 | - |
| o3_realdata | 18 |  |  |  |  |  | 18 |  | archive/2026-10/merged/PREREG_o3_realdata_2026-10-15 | - |
| o3_structure | 5 |  | 1 |  |  |  | 4 |  | superpowers/specs/2026-10-14-o3-structure-projection-design | - |
| o3_unexplained | 6 |  | 1 |  |  |  | 5 |  | LEDGER_O3 | 1088 |
| paf_duplicons | 4 |  |  |  |  |  | 4 |  | archive/2026-10/merged/O1_FAST_2026-10-21 | - |
| phased_toy | 1 |  |  |  |  |  | 1 |  | LEDGER_O3 | - |
| popfix | 12 |  |  |  |  |  | 12 |  | archive/2026-10/merged/PREREG_popfix_2026-10-12 | - |
| psv_ceiling | 2 |  |  |  |  |  | 2 |  | archive/2026-10/PSV_CEILING_2026-10-03 | - |
| psv_sun | 5 |  |  |  |  |  | 5 |  | archive/2026-10/merged/PSV_SUN_PILEUP_2026-10-12 | - |
| rep_rule | 2 |  |  |  |  |  | 2 |  | LEDGER_O1 | - |
| rna_allele | 27 |  | 3 |  |  |  | 24 |  | LEDGER_O3 (rna_allele_haplotype_count) | - |
| sedef | 7 |  |  |  |  |  | 7 |  | superpowers/plans/2026-10-13-o1-layers-phase-1-2; register 1278-1280 | - |
| seed_pool | 10 |  |  |  |  |  | 10 |  | archive/2026-10/merged/PREREG_seed_pool_real_reads_2026-10-07 | 1262 |
| soto | 5 |  |  |  |  |  | 5 |  | o1_ledger | 1105,1430 |
| soto_m2 | 22 |  |  |  |  |  | 22 |  | archive/2026-09/PREREG_unit_cover_2026-09-30 | - |
| soto_subset | 8 |  |  |  |  |  | 8 |  | superpowers/plans/2026-10-14-soto-subset-cn | - |
| sqanti_single | 1 |  |  |  |  |  | 1 |  | PREREG_sqanti_single_copy_2026-11-01 (r1482) | - |
| stgaps | 17 |  | 1 |  |  | 1 | 15 |  | LEDGER_O1; stgaps/FAMILY_LAYER_SPEC.md | 1423,1428,1479 |
| template_copies | 14 |  |  |  |  |  | 14 |  | superpowers/plans/2026-10-14-template-copies | - |
| tool_difference | 12 |  |  |  |  |  | 12 |  | LEDGER_O1/O2 | - |
| unmapped_rescue | 103 |  |  |  |  |  | 103 |  | O3_STATUS; THESIS_OBJECTIVES (o3_export.py) | 1065 |
| yag_o2 | 6 |  |  |  |  |  | 6 |  | archive/2026-10/merged/PREREG_popfix_2026-10-12 | - |

## Portability (PIPELINE, KEEP-PY, tools/*.sh)

Shipped driver code is clean: `rustle_pipeline.sh`/`rustle_reassemble.sh` have no paths or contig names outside comments; THREADS default 4 is overridable.

| file:line | problem | fix |
|---|---|---|
| tools/mm2_shard.sh:64 (doc :33) | cache root default `/mnt/linuxdisk/tmp/mm2_shard_cache` | default `${TMPDIR:-/tmp}/mm2_shard_cache` or under PREFIX |
| tools/mm2_shard.sh:113 | `THREADS=3` fixed for the old 5-core box | take minimap2 `-t` / `${RUSTLE_THREADS}` |
| tools/identity_check.sh:5,13,17,70-71 | golden dir `/mnt/linuxdisk/home/juanfraitu/...` and gorilla BAM slice | require DIR / env vars; dev-only |
| tools/identity_check.sh:57,75-82 | `--threads 2/4` fixed | `${THREADS:-4}` |
| tools/rlock.sh:13 | `D=/mnt/linuxdisk/tmp` lock dir; WSL2 5-core assumptions | archive (Slurm on Fritz) |
| bench/denovo_eval/famtruth.py:2,28-30,83 | human CHM13 v2.0 GENCODE/CAT-Liftoff + gorilla GGO_genomic.gff defaults under /mnt/linuxdisk; "Human only" | require `--genes GFF --paralogs TSV` (generic gene GFF3 + Ensembl paralogue export); no defaults |
| bench/denovo_eval/famtruth.py:42-51 | MAINSTAYS / GGO_MAINSTAYS / YAGS human/ape gene-symbol regexes | move to an optional per-species TSV |
| bench/denovo_eval/pairscore.py:9,25-35 | `--species human\|gorilla\|hg002` switch selects hard-wired truth loaders | species-free: truth from files |
| bench/denovo_eval/pairscore.py:44 | BioMart column `Human paralogue gene stable ID` | accept any `* paralogue gene stable ID` column |
| bench/denovo_eval/pairscore.py:47 / seqscore.py:26 | HGNC groups (`DENOVO_HGNC`), HG002 xref via guided3 | HGNC optional (already env); hg002 loader -> famtruth |
| bench/denovo_eval/seqscore.py:9,16,29 | `--species human\|gorilla`, default gorilla | as pairscore |
| bench/guided_eval/guided3.py:23-25 | chm13v2.0.fa, GGO.fasta, `/mnt/linuxdisk/tmp/hgnc/hierarchy.csv` hard-coded | do not ship; merge loader into famtruth |
| bench/units_from_annotation.py:8-11,97 | readthrough = RefSeq `description` contains "readthrough" | also honour Ensembl/GENCODE `gene_biotype`/tag `readthrough_transcript` |
| bench/sim.py:9,1007-1009 | defaults `/mnt/linuxdisk/home/juanfraitu/npip_cat/*`, rustle_target bin; chr20 note | require args; `--bin` default target/release |
| bench/sim.py:608 | minimap2 forced `-t 1` | fine for determinism; document |
| bench/score.py:8-10 | usage example paths `/mnt/linuxdisk/tmp/gw22`, NC_073244.2 | docstring only |
| bench/truth.py:476,696-743 | compara_chr16 path note; `nodes_gencode_like` assumes GENCODE/CAT gene records | document as one of several annotation readers |
| tools/docs_index.py:11; tools/audit_cleanup_candidates.py:323-324 | memory dir `~/.claude/projects/-mnt-c-Users-jfris-Desktop/memory` (old WSL machine) | flag/env; dev-only |
| tools/protein_attach.py:1-8 | docstring cites driver stage `families-protein` that the driver no longer has | fix docstring or restore stage |
| figures/figlib.py:128,346,354; figures/_liftoff.py:177,221-224,577; figures/samples.py:182,387,767; figures/fig_loop.py:107 | work/TMPDIR/mmi defaults under `/mnt/linuxdisk`, Windows font dirs | take `work` from inputs.tsv only; `$TMPDIR` |
| figures/_o1.py:66-67; _o1_recovery.py:83-85; _o2.py:78-83; fig_*.py SAMPLE_ORDER/labels | human chr16/chr6, gorilla NC_073244.2/NC_073234.2, sample ids A119b/OR6737/KB3781 hard-coded | thesis-only; move to samples.tsv if figures are ever released |
| bench/loop_home.py:19,401 | HG002 named in `--truth` help only | none needed |

## Counts

| status | n |
|---|---|
| PIPELINE | 5 |
| PORTED | 17 |
| KEEP-PY | 42 |
| TO-PORT | 0 |
| TO-MERGE | 9 |
| ARCHIVE | 602 |
| UNKNOWN | 12 |
| total | 687 |

## Actions (ordered)

1. **TO-PORT: none required.** The default pipeline is all Rust. Conditional: if the closed loop (`rustle_reassemble.sh`, union certificate FAILED r1481) is kept, port `bench/loop_home.py union` into `copy_assign --union-certificate` (emit `PREFIX.loop.home.tsv` directly); otherwise archive `rustle_reassemble.sh` + `loop_home.py` + `test_loop_home.py` together.
2. **Detach psv_vcf.py from o3port.py** (it imports 18 symbols: pile_records, observations, informative, ...): vendor them into psv_vcf.py or emit the VCF from `copy_phasing`; then o3port.py and its golden harness can be archived.
3. **Merge**: guided3.py/cdsseed.py hg002 loader -> famtruth.py; tests/fixtures `make_bam.py` x2 -> shared helper; archive-bound pairs (hapcn_mat_pat_robust, o2_p18_features_cond, sdcn_noise_ex, o2_p20_sim19, polya_mono) fold into their survivor or simply archive both.
4. **Release package (`rustle-eval`)**: pairscore, seqscore, famtruth, end_score, psv_vcf, units_from_annotation, score/lib/sim/truth, protein_attach (+tests). Move them out of bench/ into e.g. `python/rustle_eval/`.
5. **Portability fixes before release** (table above): famtruth/pairscore/seqscore species-free truth inputs (biggest); mm2_shard cache root + THREADS; units_from_annotation readthrough rule for non-RefSeq GFFs; sim.py defaults; protein_attach docstring.
6. **Archive**: move every ARCHIVE file (602: 601 bench + tools/rlock.sh) to `archive/bench/<dir>/` keeping paths cited by the register; PORTED study scripts go with them once their goldens are recorded. Fix the src doc comments that cite `bench/...` paths (e.g. family.rs cites rna_allele/*, o3_structure/proto.py) to the new location.
7. **Decide UNKNOWN**: bench/famsim (12) - keep as a test simulator in the package, or archive with FAMSIM_DESIGN.
8. Dev-only tools (docs_index, audit_cleanup_candidates, identity_check, figures/) stay in the thesis repo, excluded from the release package.

