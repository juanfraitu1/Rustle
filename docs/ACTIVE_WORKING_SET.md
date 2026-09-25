# Active working set — what is actually in use (2026-09-19; wave 4 applied 2026-09-23)

> **Wave 7, Rust side (2026-09-24):** catalog builder **9.4× faster** (human chr16 3,495 → 373 s, 13.9 → 7.0 GB; gorilla
> contig 292 → 45-53 s), every stage `cmp`-identical to the committed binary (727326fc) incl. three opt-in collapse modes:
> the locus collapse's all-pairs serial POA → span sweep + rayon + an exact substring bound (`core_coverage_reaches`);
> a coordinates-only second BAM pass; three O(n²) loops and the families-stage gene scan → exact indexed versions; rayon
> pool bounded to `--threads`. New **`run_cache.rs`** (`RUSTLE_CACHE_DIR`, driver default `PREFIX.cache`): reps after the
> collapse + every all-vs-all PAF, for fast re-runs and inspection (`--inspect`, `cache-ls`). Dead code removed: StringTie
> residue in `bam.rs`/`types.rs`, the `FamilyGraph` half of `family_graph.rs`, `minimizers.rs`/`bridge_detector.rs`/
> `repeat_catalog.rs` (survivors moved), `util/`, 24 fixtures (16 MB → 52 KB), crates `noodles`/`fixedbitset`/`rand`/
> `serde`; refuted switches `RUSTLE_JUNCTION_FUZZ_BP`, `_TERMINAL_QUANTILE`, `_ER_SUM_COVERAGE`, `_ER_CORE_COVERAGE`,
> `_SPLICED_REP`, `_LOCUS_LINK_MIN_READS` (their params.tsv rows kept as constants); `tools/` 5 → 2 (3 to the attic).
> `src/**/*.rs` 63,572 → 56,845 lines (56 → 51 files); tests 857 passed / 0 failed.
>
> **Wave 7 (2026-09-24, bench consolidation; same brief as wave 6):** `bench/` Python **42 → 9** (removed files at tag
> `notebook-2026-09-24` and in `~/Desktop/Rustle_attic/2026-09-24/`, manifest there). The 21 top-level scripts became
> `lib.py` (shared helpers) + `score.py` (14 scorer subcommands) + `sim.py` (5 simulators) + `truth.py` (4 truth
> builders), with `guided_pipeline.py` and `mcl_port.py` kept. The 12 `layer_order/` files became `lattice_common.py` +
> `npip_tbc1d3.py` (one subcommand per stage). The 9 `soto/` files became `soto_replication.py` (6 subcommands), with its
> two archived inputs restored. Every replaced command was `cmp`-checked against its old script on recorded inputs.
> Defects found and fixed: `soto_vs_us_referee.py` had raised NameError on every run since 8db314c7; the O2/O3
> simulators seeded their RNGs with per-process `hash()`, so no run was reproducible (now `stable_seed()`, which means
> the new reads differ from every earlier run); 4 of the 10 `layer_order` entry points could not run at HEAD
> (repo-root path bug), and `lattice_check_c2` had been dead since 2026-09-19; `adjudicated_truth.py build` had crashed
> in `mcl_port` since 21d6c5c9. Register rows 1104 (layer_order) and 1105 (the wave). **`bench/README.md` is the
> inventory and the old-name → new-command table**; the three redirect tables below ("Finding a script") cover the
> same moves. `bench/mcl_port.py` is now imported only by `bench/lib.py`.
> **Wave 6 (2026-09-23/24, consolidation; user: "fewer scripts, more cohesive and reusable, faster"):**
> Rust binaries **20 → 10** (retired to tag `notebook-2026-09-23c` / `~/Desktop/Rustle_attic/2026-09-23c/`: `asj`,
> `asj_verify` (dropped objective), `debug_poa`, `bam_null_probe` (probes), `index_bam`, `bam_header`, `filter_bam_by_as`
> (samtools does these), `gamma_refine`, `mcl_refine` (measurement scaffolding), `family_define` (legacy parity fixture
> whose input table was already retired)), and with them the 12 library modules only they reached (see
> `docs/MODULE_STATUS.md`; `lgamma` inlined into `missing_copy_flag_pass`). `mcl_families --from-gtf` runs the de novo family
> stage in one command (loci from the assembled GTF → all-vs-all → MCL), replacing the scratch `loci_from_gtf.py` +
> a hand-run minimap2. **`tools/rustle_pipeline.sh`** runs every stage (assemble → families → catalog → assign → o3)
> with the shipped defaults; verified end to end on a one-contig testis slice (10.6 min, 15.7 GB peak from the index
> loads). Python 42 → 38 by merging the O2 sim/score pair (`copy_assign_read_truth.py sim|score`) and the bakeoff pair
> (`copy_assign_tool_bakeoff.py calls|compare`) and inlining two single-use helpers. Speed: `missing_copy_flag` decodes
> insertion-carrying reads once in pass 1 (a structural-only locus no longer re-queries the BAM): 42 → 28 s on an
> 80 Mb contig, identical output. Tests 799 + integration, 0 failures.
> **Wave 5 (2026-09-23, `tools/cleanup_wave5_scripts.sh`, user: "only O1, O2, O3, the Soto replication and the live
> parts"):** `bench/` reduced from 212 Python scripts to **42** (list = the script's KEEP block + its import closure;
> `bench/README.md` is the reviewable inventory). 180 entries moved to `~/Desktop/Rustle_attic/2026-09-23b/` after
> tagging HEAD `notebook-2026-09-23b`: every refuted-arm and one-off probe script, the 13 `gen_*_fixture.py`
> generators (their imports were archived in wave 2, the JSON fixtures they made are committed), the Python parity
> oracles of ported Rust bins (`assembly_polish.py`, `readthrough_secondary_filter.py`, `o1_eval.py`, `igv_tracks.py`),
> the bench shell scripts, and `family_rna_refine.tsv` (r1033). The O3 scratch simulations were consolidated into one
> kept script, `bench/missing_copy_sim.py` (modes transcript / genomic / shuffled). Verified: 42/42 compile, no kept
> script imports a moved one.
> **Wave 4 (2026-09-23, `tools/cleanup_wave4_attic.sh`):** 73 entries / 3.6 MB moved to
> `~/Desktop/Rustle_attic/2026-09-23/` (manifest there): 23 superseded documents (the pre-September method
> docs `ONE_METHOD`/`NUMBERS`/`OBJECTIVES_AND_VERIFICATION`/`o1_investigations`/`METHOD_PSEUDOCODE`, the old
> README, `OPEN_ITEMS_2026-09-09`, `o3_missing_copy_evidence`, the sdd specs/plans in `docs/superpowers/`,
> `docs/{archive,artifacts,experiments}/`, four stale `bench/*.md`), 27 earlier data tables from the 09-04/05
> sweeps, 4 superseded scripts (the applied cleanup waves, `family_vg_report.py`; ⚠ `bench/mcl_port.py` was moved and put back — 8 scripts import it as a library (5 by wave 6's end; since wave 7 only `bench/lib.py`, lazily), the Rust `mcl_port` bin ports only its CLI),
> and 19 caches/run outputs. Tracked files were `git rm --cached` after tagging HEAD `notebook-2026-09-23`.
> `docs/CLEANUP_CANDIDATES.md` and its TSV (the 09-22 audit) went with them; re-run
> `tools/audit_cleanup_candidates.py` to regenerate.

> Waves 1 and 2 APPLIED (§6q8, §6q9). Wave 1: 101 tracked files → `archive/`, 33 untracked source →
> `archive/untracked/`, 50 caches/logs deleted. Wave 2: **`bench/` inverted — 428 of 781 scripts archived,
> leaving 353** (the 331 an authoritative document names, plus a 22-script import closure). After both:
> tests 883/0 + 28/0, and `bench/assembly_polish.py` still reproduces the Rust polish byte-for-byte.

Companion to `docs/CLEANUP_CANDIDATES.md`, which marks what is *probably dead*. This file marks what is
**live**, so the audit stops treating it as a candidate: `tools/audit_cleanup_candidates.py` anchors any
file named by a top-level `docs/*.md`, so everything listed here is protected by being listed here.

Re-derive the tiers with `git log --name-only`, the ledger, and the audit TSV; re-run the audit after
editing this file.

## Tier 1 — the live pipeline

What the current §6q7 workflow actually invokes. Touch these with care; they are load-bearing.

### Rust entry points (`src/bin/`, PROTECTED by the audit)

| binary | role |
|---|---|
| `copy_assign` | ⭐ the main pipeline: loci → isoform assembly → GTF; `--assemble-only` is the assembler product, `--assembly-polish` the §6p8-§6q6 filters; O2 assignment lives here too |
| `mcl_families` | the DNA-level family definition (`--min-exonic-bp`, `--min-shared-exon-frac`) |
| `family_define`, `mcl_refine`, `gw_family_catalog` | family catalog construction and refinement |
| `gamma_refine`, `parcn`, `filter_bam_by_as`, `index_bam`, `bam_header` | supporting steps |
| `asj`, `asj_verify` | ⚠ ASJ is a DROPPED objective — kept for provenance, not in the live path |
| `debug_poa` | diagnostic only |

### Scripts

| script | role | last used |
|---|---|---|
| `bench/assembly_polish.py` | Python mirror of the Rust `--assembly-polish` passes; byte-identical parity oracle | §6q6, 2026-09-19 |
| `readthrough_filter` (Rust bin) | ⭐ secondary-dominated readthrough flagging (opt-in, `--max-secondary-frac`, default 1.0 = no-op). Needs no `samtools` and no off-repo `dna_cert`; `--nodes`/`--exonless` default to the family_cert substrate | §6s4, 2026-09-20 |
| `bench/readthrough_secondary_filter.py` | the Python it was ported from; kept as the byte-parity oracle | §6n9/§6o0 |
| `gff_to_gtf` (Rust bin) | RefSeq GFF3 → gffread-style GTF; validated at 4,574 = 4,574 vs `chr20_ref.gtf`. **Needed because `gffread` is not installed on this machine.** | §6p9-§6r9 |
| `locus_bed` (Rust bin) | loci as BED + one-to-one match against an annotation (`size_ratio`) | §6r5-§6r9 |
| `tools/audit_cleanup_candidates.py` | this audit; read-only, re-runnable | §6q8 |
| `bench/layer_order/npip_tbc1d3.py` + `lattice_common.py` | the 2026-09-16 NPIP/TBC1D3 nested edge-test lattice + layer-order study (L3 at 0.98, not the shipped 0.985); one subcommand per stage (wave 7 folded the 12 `lo_*`/`lattice_*`/`soto_map` files into these 2). §6p1–§6p5 used the off-repo `lattice_rules/engine.py`, which reads `light/members.corrected.tsv` from `corrected-tables` | 2026-09-16 (outputs frozen); rerun-verified 2026-09-24 |

### Off-repo, but part of the live path

- `/mnt/linuxdisk/home/juanfraitu/lattice_rules/engine.py` — `L3_CUT`/`L4_CUT`, `l4_refine()` (§6p1/§6p5).
- `/mnt/linuxdisk/tmp/flair_shims/` — shims that make FLAIR 3.0.0 runnable here (§6q6); a second, broken
  FLAIR install on `linuxdisk` shadows the working one.
- Benchmark substrates under `/mnt/linuxdisk/home/juanfraitu/bakeoff/` — `human_chr{20,11,7,14,5,9}`
  (shallow `human_testis.t2t.bam`), `a119b_chr20` and `ggo_NC_073244.2` (the lab's deep libraries).
  ⚠ Register row 867: the two human libraries are **not** interchangeable.

## Tier 2 — reproduces a recorded result

**152 scripts are named directly in `docs/o1_ledger.md`, `docs/NEGATIVE_RESULTS_REGISTER.md` or a
`docs/PREREG_*.md`.** They are the provenance of published numbers: not live, but deleting one makes a
ledger claim unreproducible. Archive, never delete. Enumerate them with:

```sh
grep -ohE '(bench|tools|scripts|analysis)/[A-Za-z0-9_./-]+\.(py|sh)' \
     docs/o1_ledger.md docs/NEGATIVE_RESULTS_REGISTER.md docs/PREREG_*.md | sort -u
```

## What is in the folder but not in the repository

The working folder is not the repository. After §6s7 it is **289 MB**, of which `.git` is 257 MB; a
clone is **40 MB / 587 files**. Untracked caches and old run output (210 MB — ~107 one-off GTFs, a
103 MB `ri_kmer_ckpt.npz`, 16 run logs, the retired R scripts' data) were moved to
**`~/Desktop/Rustle_attic`**, which has its own README. None of it was ever committed, so none of it
is recoverable from git; delete that folder when you no longer want the old outputs.

⚠ `.gitignore` is an **allowlist** (`*`, then negations). A new top-level file is invisible to git
until it is negated there — that is how `REPRODUCE.md` stayed out of every clone (register 900).
After adding a root file, check `git ls-files --error-unmatch <file>`.

## Finding a script or output a document names

The working tree holds **585 files: source, documentation and test fixtures only** — 189 `.py` and
8 `.sh`, down from 342/30 before §6s1-§6s5. Everything else lives at a git tag, not on disk:
**`notebook-2026-09-19`** (1,720 entries: the §6r9 sweep) and **`notebook-2026-09-20`** (251 entries:
the retired-era R/probes and 157 superseded provenance scripts). Check the newer tag first.

```sh
git ls-tree -r notebook-2026-09-20 archive/ | grep <name>      # find it (then try -09-19)
git checkout notebook-2026-09-20 -- archive/bench/<name>.py    # bring it back
git show notebook-2026-09-20:archive/bench/<name>.py           # just read it
```

**So: if a document names `bench/X` and it is not there, it is `archive/bench/X` at that tag.** All of it
went in as tracked `git mv`s, so `git log --follow <old path>` still works too.

Seven paths named in older documents exist nowhere and predate this cleanup: three were deleted in
earlier `chore: prune` commits (`667f2e5c`, `a7d003a3` — recoverable from history) and four were never
tracked (`bench/em_coverage_sweep.py`, `bench/sun_catalog_fast.py`,
`bench/soto/soto_segdup_cn_refine.py`, `bench/copy_recovery_eval/results_genomewide/gw_run.sh`).

⚠**Five Tier-2 scripts were moved by wave 1 (§6q8) and now live under `archive/`.** The ledger and
register still name their original paths, so follow them there:

| ledger path | now at |
|---|---|
| `bench/denovo_shared_def.py` | `archive/bench/denovo_shared_def.py` |
| `bench/missing_copy_flag_pass.py` | `archive/bench/missing_copy_flag_pass.py` |
| `bench/vg_repeat_catalog.py` | `archive/bench/vg_repeat_catalog.py` |
| `bench/gw_rebuild.sh` | `archive/bench/gw_rebuild.sh` |
| `bench/gw_rebuild_v2.sh` | `archive/bench/gw_rebuild_v2.sh` |

They are tracked moves, so `git log --follow` and `git show HEAD~1:<old path>` both still work.

⚠**Wave 7 (2026-09-24) folded `bench/layer_order/` (12 files) into 2.** The ledger, the register and
`bench/TBC1D3_GUITART_TRUTH_CORRECTION.md` cite the old files by path and line. Those citations resolve at tag
`notebook-2026-09-24` (`git show notebook-2026-09-24:bench/layer_order/<file>`). To rerun, use:

| old path | now |
|---|---|
| `bench/layer_order/lo_expr_recount.py OUT [--ignore IDS] EXTRA…` | `bench/layer_order/npip_tbc1d3.py expr-recount OUT [--ignore IDS] EXTRA…` |
| `bench/layer_order/lo_corrected_tables.py` | `npip_tbc1d3.py corrected-tables` |
| `bench/layer_order/lo_analysis.py` | `npip_tbc1d3.py layer-order` |
| `bench/layer_order/lattice_edges.py` | `npip_tbc1d3.py lattice-edges` |
| `bench/layer_order/lattice_expr.py` | `npip_tbc1d3.py lattice-expr` |
| `bench/layer_order/lattice_levels.py` | `npip_tbc1d3.py lattice-levels` |
| `bench/layer_order/lattice_truth.py` | `npip_tbc1d3.py lattice-truth` |
| `bench/layer_order/lattice_filtration.py [--l1 c2_loose]` | `npip_tbc1d3.py lattice-filtration [--l1 c2_loose]` |
| `bench/layer_order/lattice_check_c2.py` | `npip_tbc1d3.py lattice-check-c2` (had not run since 2026-09-19) |
| `bench/layer_order/lattice_report_tables.py` | `npip_tbc1d3.py lattice-report` |
| `bench/layer_order/soto_map.py` (`load`, `map_gene`, `FIELDS`) | `lattice_common.soto_load`, `soto_map_gene`, `SOTO_FIELDS` |
| `lo_analysis.c_tree` / `bip_jaccard` / `hgnc_all` | `lattice_common.c_tree` / `bip_jaccard` / `hgnc_all` |
| `lattice_common.py` | unchanged name and exports (the off-repo `LAT/corrections_pass2/*.py` import it by path) |

Every stage takes `--root DIR` (or `LO_ROOT`) and overwrites its outputs there. The default is the frozen
`/mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3`, so to rerun, point it at a copy.

⚠**Wave 7 (2026-09-24) folded the 8 `bench/soto/*.py` Soto-replication scripts into `bench/soto/soto_replication.py`,
and removed `bench/soto/rustlib.py` (0 importers).** Ledger §6ie–§6ip and the register cite the old names; they resolve
at tag `notebook-2026-09-24`. The recipe and the §6ip numbers are in `REPRODUCE.md` §5a. To rerun (same flags unless
noted):

| old path | now |
|---|---|
| `bench/soto/soto_replicate_from_sedef.py` | `bench/soto/soto_replication.py edges` (`--s1e` optional, default `bench/soto/soto_parCN_S1E.tsv`) |
| `bench/soto/soto_cluster_from_shared.py` | `soto_replication.py cluster` (default `--mad-statistic mean` kept) |
| `bench/soto/soto_cluster_dennislab_algorithm.py` | `soto_replication.py dennislab` (default `--mad-statistic median` kept) |
| `bench/soto/famcn_from_wssd.py` | `soto_replication.py famcn` (`--s1e` default is now relative to the module, same file) |
| `bench/soto/soto_score_against_truth.py` | `soto_replication.py score --only pairs` (`--truth` optional, default S1C) |
| `bench/soto/soto_bipartite_match_score.py` | `soto_replication.py score --only bipartite`; plain `score` prints both |
| `bench/soto/soto_attach_noncoding_members.py` | dropped: `cluster --full-geneset` (§6ii: same partition); `attach()` kept |
| `bench/soto/soto_replicate_clustering.py` | dropped (minimap2 map-back path, superseded by §6ie); `load_exons()` kept |
| hand-made `soto_{1793,2334}_geneset.tsv` | `soto_replication.py genesets` (new; derived from S1C) |
| `bench/soto/rustlib.py` | not replaced; `git show notebook-2026-09-24:bench/soto/rustlib.py` |

Restored next to the module: `soto_parCN_S1E.tsv` and `acro_extra_anchors.tsv` (from `cd37ccb0^`; wave 3 archived them),
plus the frozen `edges` output `shared_exons_2334_finalhuman.tsv`.

⚠**Wave 7 (2026-09-24) folded the 21 top-level `bench/*.py` scripts into `lib.py` + `score.py` / `sim.py` / `truth.py`**
(`guided_pipeline.py` and `mcl_port.py` stay). The ledger, the register, the PREREGs and the result docs cite the old
names; they resolve at tag `notebook-2026-09-24`. To rerun (same arguments unless `bench/README.md` notes otherwise; its
table also gives the names these scripts had before 2026-09-24, such as `o2_read_truth.py` and `o3_sim_copies.py`):

| old path | now |
|---|---|
| `bench/referee_band_score.py`, `bench/identity_spectrum.py --catalog` | `bench/score.py pairs` (flags differ: see `bench/README.md`) |
| `bench/identity_spectrum.py` (tier mode) | `bench/score.py spectrum` |
| `bench/heldout_family_score.py` | `bench/score.py heldout` |
| `bench/soto_vs_us_referee.py` | `bench/score.py referee` (the old file raised NameError from 8db314c7 on) |
| `bench/protein_edge_gap.py` | `bench/score.py edge-gap` |
| `bench/rna_truth_from_protein.py` | `bench/score.py rna-ceiling` |
| `bench/member_completeness.py` | `bench/score.py members` (positional → flags) |
| `bench/adjudicated_truth.py score` / `build` | `bench/score.py adjudicated` / `bench/truth.py adjudicated` |
| `bench/protein_families.py score` / `build` | `bench/score.py protein` / `bench/truth.py protein` |
| `bench/annotation_nodes.py` | `bench/truth.py nodes` |
| `bench/eichler_compare.py` | `bench/score.py eichler` |
| `bench/copy_assign_read_truth.py sim` / `score` | `bench/sim.py copies` / `bench/score.py reads` |
| `bench/copy_assign_tool_bakeoff.py calls` / `compare` | `bench/score.py bakeoff-calls` / `bakeoff-compare` |
| `bench/copy_assign_excision.py` | `bench/sim.py excise` |
| `bench/missing_copy_sim.py` | `bench/sim.py missing-copy` |
| `bench/tandem_copy_sim.py` | `bench/sim.py tandem` |
| `bench/ideal_chromosome_sim.py` | `bench/sim.py chromosome` |
| `bench/sim_reads.py` (`simulate_reads`) | `sim.simulate_reads` |
| `bench/locus_reads.py` (`reads_with_block_in`, `spanning_genes`) | `bench/lib.py` (same names); CLI `bench/score.py locus-reads` |

⚠ `sim.py copies` and `sim.py missing-copy` seed with `stable_seed()`, so they draw different reads from every earlier
run (the old ones used per-process `hash()` and could not be reproduced anyway); re-measure before comparing.

## Tier 3 — anchored only transitively (the review pool)

Of 973 scripts in the repo, 575 are KEEP-CITED but only 152 are named directly; **the other 423 are
anchored only through another un-verified file.** ⚠ The 2026-09-16 round-1 verification found the
KEEP-CITED control itself was **4/12 dead** against a ≤25% bar, so transitive anchoring over-keeps. Tier 3
is where a purge should look after the candidate classes are dealt with — but only with per-file checks,
not in bulk.

## Numbers behind the tiers

| | files | note |
|---|---|---|
| scripts in repo | 197 | `.py` + `.sh`, after §6s5 (was 973 at the §6q8 audit) |
| named directly in ledger/register/prereg | 152 | Tier 2 — the KEEP seed |
| archived at `notebook-2026-09-20` | 251 | R figures, VG-HMM/StringTie probes, superseded provenance |
| archived at `notebook-2026-09-19` | 1,720 | the §6r9 sweep |
| Rust | 80 `.rs` / 69,824 LOC | vs 189 `.py` / 35,915 LOC — Rust is now ~2x the Python |

## Third-party tools — NOT vendored

The repo vendors no third-party assemblers. Install them yourself:

| tool | install | used by |
|---|---|---|
| StringTie | `conda install -c bioconda stringtie` | `bench/bakeoff_*_stringtie.sh`, `bench/gw_threeway.sh` (they take `$STRINGTIE` or resolve it on `PATH`) |
| gffcompare | `conda install -c bioconda gffcompare` | every bakeoff |
| FLAIR 3.0.0 | `conda create -n flair -c bioconda flair` | §6q6. ⚠A second, broken FLAIR install can shadow the working one; see the shim note above |
| SQANTI3 | checkout + its conda env | §6q5 |
| isoseq | `conda install -c bioconda isoseq` | ⚠needs native PacBio read names AND a relaxed `--min-aln-coverage`; the lab's runs are in `~/Desktop/isoseq_upload/` |
| gffread | not installed here — `tools/refseq_gff_to_gtf.py` stands in for `gffread -T` | reference GTF construction |

⚠**`tools/stringtie` was a git submodule and was removed (§6q8).** Its URL was `../stringtie`, a relative
path to a sibling checkout that does not exist — so `git clone --recursive` failed for anyone, including
on the original machine. The 55 MB checkout was moved to `~/Desktop/stringtie` locally; nothing in the
repo needs it.

⚠**Tool builds differ across the six-chromosome panel**: chr20's StringTie arm is 3.0.1 (that vendored
build), chr11/7/14/5/9 are 3.0.3 (conda). Register row 870.

## Retired because Rust replaced them (§6r9)

Removed from the tree; recover from history with `git show <rev>:<path>` (they are ordinary deletions,
after the `notebook-2026-09-19` tag, so the archive redirect above does not cover them).

| removed | replaced by | evidence |
|---|---|---|
| `tools/refseq_gff_to_gtf.py` | **`gff_to_gtf`** binary | byte-identical on chr20/chr16/chr7, **7.4× faster** (2.7 s vs 20 s) |
| `bench/locus_bed.py` | **`locus_bed`** binary | all three outputs byte-identical, 4× faster (0.10 s vs 0.40 s) |
| `bench/ism_collapse.py` | `bench/assembly_polish.py --support-ratio 999 --mono-quantile 0`, itself mirrored by `--assembly-polish` | byte-identical on chr20 (784 transcripts both ways) |

⚠The Python `refseq_gff_to_gtf.py` had its body DUPLICATED (an earlier docstring edit appended instead of
replacing), so it did the whole conversion twice — which is most of why it measured 20 s.
