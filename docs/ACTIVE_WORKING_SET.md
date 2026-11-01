# Active working set — what is actually in use (2026-09-19; wave 4 applied 2026-09-23)

Compact version (2026-11-01); the full wave history and every table verbatim: `archive/2026-11/ACTIVE_WORKING_SET.md`.
⚠ `tools/audit_cleanup_candidates.py` anchors files named in top-level `docs/*.md`; names dropped here survive only in the archive copy, so re-run the audit with `docs/archive/2026-11/` in scope before any purge.

## Current state (2026-11-01)
- **Rust binaries today (Cargo.toml, 14):** `mcl_families`, `copy_assign`, `gw_family_catalog`, `readthrough_filter`, `family_score`, `missing_copy_flag`, `candidate_copies`, `as_table`, `copy_phasing`, `utilities`, `family_gfa`, `semisec`, `semifam`, `slimbam`. Source is flat in `src/` (`copy_phasing/` module = the Rust O2, byte-identical to `bench/o3_graph/o3port.py`, register 1402). Module states: `MODULE_STATUS.md`.
- **Driver:** `tools/rustle_pipeline.sh` (stages assemble -> families -> catalog -> assign -> o3; `all --to/--from STAGE`); also `tools/rustle_reassemble.sh`, `tools/mm2_shard.sh`, `tools/minimap2_gate`, `tools/protein_attach.py`, `tools/docs_index.py`. Next target: one `rustle` binary + run directory (memory "next targets").
- **Live bench (2026-10/11):** `bench/denovo_eval/` (`pairscore.py`, `famtruth.py`, `semifam.py`, `sketchfam.py`, `dispersion.py`, `liftfloor.py`, `edgerecall.py`), `bench/o3_graph/` (`o3port.py`, `real_net.py`, `tie_ms_vs_as.py`, `psv_vcf.py`), `bench/o2levers/` (`levers_score.py`, `panel.families.txt`), `bench/stgaps/famasm.py` (FAMASM prototype, r1479 FAIL).
- **Machine:** Fritz cluster (paths in `DATA.md`, memory `fritz-env`). Laptop-only items below (`/mnt/linuxdisk`, `~/Desktop/Rustle_attic`, FLAIR shims) are not on Fritz. Third-party tools come from pixi envs (minimap2 2.30, samtools/bcftools 1.25, SQANTI3 6.1, gffcompare/gffread, IsoCon 0.3.3), so "gffread is not installed" is laptop-only.
- **Cleanup waves (history):** waves 1-2 (§6q8, §6q9; bench 781 -> 353 scripts), wave 4 (73 entries to the attic, tag `notebook-2026-09-23`), wave 5 (bench 212 -> 42, tag `notebook-2026-09-23b`), wave 6 (Rust binaries 20 -> 10, tag `notebook-2026-09-23c`; `mcl_families --from-gtf`; `tools/rustle_pipeline.sh`), wave 7 (bench Python 42 -> 9, tag `notebook-2026-09-24`, register 1104-1105; catalog builder 9.4x faster, `run_cache.rs`). `bench/README.md` is the inventory and old-name -> new-command table.

## Tier 1 — the live pipeline
As of 2026-09-19 (superseded by the binary list above; wave 6 retired `family_define`, `mcl_refine`, `gamma_refine`, `filter_bam_by_as`, `index_bam`, `bam_header`, `asj`, `asj_verify`, `debug_poa`). Scripts then live: `bench/assembly_polish.py` (parity oracle of `--assembly-polish`), `readthrough_filter` (opt-in `--max-secondary-frac`), `bench/readthrough_secondary_filter.py`, `utilities gff-to-gtf` (4,574 = 4,574 vs `chr20_ref.gtf`), `utilities locus-bed`, `tools/audit_cleanup_candidates.py`, `bench/layer_order/npip_tbc1d3.py` + `lattice_common.py`. Off-repo (laptop): `/mnt/linuxdisk/home/juanfraitu/lattice_rules/engine.py`, `/mnt/linuxdisk/tmp/flair_shims/`, bakeoff substrates `/mnt/linuxdisk/home/juanfraitu/bakeoff/` (register 867: the two human libraries are not interchangeable). (full: archive/2026-11/ACTIVE_WORKING_SET.md:72)

## Tier 2 — reproduces a recorded result
152 scripts named in `docs/o1_ledger.md`, the register or a `docs/PREREG_*.md`; archive, never delete. Enumerate: `grep -ohE '(bench|tools|scripts|analysis)/[A-Za-z0-9_./-]+\.(py|sh)' docs/o1_ledger.md docs/NEGATIVE_RESULTS_REGISTER.md docs/PREREG_*.md | sort -u` (add `docs/LEDGER_*.md` and `docs/archive/` for the ledgers since 2026-10). (full: archive/2026-11/ACTIVE_WORKING_SET.md:108)

## What is in the folder but not in the repository
Laptop: untracked caches moved to `~/Desktop/Rustle_attic` (not in git). `.gitignore` is an allowlist: a new top-level file is invisible until negated (register 900); check `git ls-files --error-unmatch <file>`. (full: archive/2026-11/ACTIVE_WORKING_SET.md:119)

## Finding a script or output a document names
If a document names `bench/X` and it is not there, look at a tag: `notebook-2026-09-24` (wave 7 folds), `notebook-2026-09-23c/b`, `notebook-2026-09-20`, `notebook-2026-09-19` (`git ls-tree -r <tag> archive/ | grep <name>`, `git show <tag>:<path>`). Wave-1 moves: `bench/{denovo_shared_def.py, missing_copy_flag_pass.py, vg_repeat_catalog.py, gw_rebuild.sh, gw_rebuild_v2.sh}` -> `archive/bench/`. Seven old paths exist nowhere (e.g. `bench/em_coverage_sweep.py`, `bench/sun_catalog_fast.py`).
Wave-7 folds (old -> new):
- `bench/layer_order/lo_*`, `lattice_*` -> `bench/layer_order/npip_tbc1d3.py` subcommands `expr-recount`, `corrected-tables`, `layer-order`, `lattice-edges`, `lattice-expr`, `lattice-levels`, `lattice-truth`, `lattice-filtration`, `lattice-check-c2`, `lattice-report`; `soto_map.py` -> `lattice_common.soto_load` / `soto_map_gene` / `SOTO_FIELDS`. Default root `/mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3` is frozen; rerun on a copy.
- `bench/soto/*.py` -> `bench/soto/soto_replication.py` `edges` (was `soto_replicate_from_sedef.py`), `cluster`, `dennislab`, `famcn`, `score --only pairs|bipartite`, `genesets`; `rustlib.py` removed. Recipe: `REPRODUCE.md` §5a.
- top-level `bench/*.py` -> `bench/score.py` (`pairs`, `spectrum`, `heldout`, `referee`, `edge-gap`, `rna-ceiling`, `members`, `adjudicated`, `protein`, `eichler`, `reads`, `bakeoff-calls`, `bakeoff-compare`, `locus-reads`), `bench/sim.py` (`copies`, `excise`, `missing-copy`, `tandem`, `chromosome`), `bench/truth.py` (`adjudicated`, `protein`, `nodes`), `bench/lib.py`. `sim.py` now seeds with `stable_seed()`: re-measure before comparing. Full table: `bench/README.md`.
(full: archive/2026-11/ACTIVE_WORKING_SET.md:131)

## Tier 3 — anchored only transitively (the review pool)
423 scripts (2026-09-16 audit) anchored only through un-verified files; the KEEP-CITED control was 4/12 dead, so purge only with per-file checks. (full: archive/2026-11/ACTIVE_WORKING_SET.md:239)

## Numbers behind the tiers
2026-09-19 counts: 197 scripts, 152 Tier 2, 251 at `notebook-2026-09-20`, 1,720 at `notebook-2026-09-19`, Rust 80 `.rs` / 69,824 LOC. Stale: Rust is now flat `src/` with 14 binaries. (full: archive/2026-11/ACTIVE_WORKING_SET.md:247)

## Third-party tools — NOT vendored
StringTie, gffcompare, FLAIR 3.0.0 (laptop shims), SQANTI3, isoseq, gffread: install yourself (on Fritz: pixi envs, memory `fritz-env`). `tools/stringtie` submodule removed (§6q8). Panel tool builds differ: chr20 StringTie 3.0.1, others 3.0.3 (register 870). (full: archive/2026-11/ACTIVE_WORKING_SET.md:257)

## Retired because Rust replaced them (§6r9)
`tools/refseq_gff_to_gtf.py` -> `gff_to_gtf` (byte-identical, 7.4x faster); `bench/locus_bed.py` -> `locus_bed` (4x); `bench/ism_collapse.py` -> `bench/assembly_polish.py --support-ratio 999 --mono-quantile 0` (784 transcripts both ways). Recover with `git show <rev>:<path>`. (full: archive/2026-11/ACTIVE_WORKING_SET.md:278)
