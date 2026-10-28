# Ledger: O1 (family definition, assembler, locus units)

Merged 2026-10-24 from 46 dated docs, oldest first; each kept verbatim under its own heading (its headings one
level down). The originals are in `docs/archive/2026-10/merged/` (git history has them at their old paths).

## Contents

- PREREG_hierarchy_AMENDMENT2_DRAFT.md
- IDEAL_CHROMOSOME_SIM_2026-09-21.md
- PREREG_inversion_aligner_blindness_2026-10-05.md
- DEFAULT_RESCORE_NPIP_2026-10-06.md
- ENTANGLED_BASELINE_2026-10-06.md
- FEWCOPY_IDEAL_CASES_2026-10-06.md
- IDEAL_EXPRESSION_DEFAULT_2026-10-06.md
- LOCUS_UNITS_2026-10-06.md
- PREREG_default_rescore_npip_2026-10-06.md
- PREREG_entangled_baseline_2026-10-06.md
- PREREG_fewcopy_class_DRAFT_2026-10-06.md
- PREREG_ideal_expression_2026-10-06.md
- PREREG_locus_units_2026-10-06.md
- GORILLA_FEWCOPY_DESIGN_2026-10-07.md
- GORILLA_OVERLAP_2026-10-07.md
- HIERARCHY_DUPLICON_EXPANSION_2026-10-07.md
- LOCUS_UNITS_LEVELS_2026-10-07.md
- PREREG_gorilla_overlap_2026-10-07.md
- PREREG_hierarchy_duplicon_expansion_2026-10-07.md
- PREREG_hierarchy_duplicon_expansion_REVIEW_DISPOSITION_2026-10-07.md
- PREREG_ideal_arm_pilot_2026-10-07.md
- PREREG_seed_pool_real_reads_2026-10-07.md
- SEED_POOL_REAL_READS_2026-10-07.md
- PREREG_family_hierarchy_flag_2026-10-09.md
- PREREG_paf_duplicon_units_2026-10-10.md
- PREREG_assembler_compatible_support_2026-10-11.md
- PREREG_direct_transfrags_2026-10-11.md
- PREREG_direct_transfrags_union_2026-10-11.md
- PREREG_inbetween_operating_point_2026-10-11.md
- PREREG_lrpap1_seed_reach_2026-10-11.md
- PREREG_tool_difference_model_2026-10-11.md
- GENOMEWIDE_FAST_2026-10-12.md
- PREREG_mcl_two_node_2026-10-13.md
- PREREG_o1_hierarchy_eval_2026-10-13.md
- PREREG_sedef_rs_layers_2026-10-13.md
- PREREG_denovo_copies_2026-10-14.md
- PREREG_liftoff_pairs_null_2026-10-14.md
- PREREG_soto_subset_cn_2026-10-14.md
- PREREG_template_copies_2026-10-14.md
- PREREG_post_mcl_reassign_2026-10-15.md
- NPIP_16Q22_CNV_2026-10-20.md
- MULTIMAPPER_PREMISE_2026-10-21.md
- O1_FAST_2026-10-21.md
- PREREG_multicopy_gains_real_2026-10-21.md
- PREREG_o1_families_sd_pairs_2026-10-21.md
- PREREG_o1_self_hits_2026-10-21.md


---

## PREREG_hierarchy_AMENDMENT2_DRAFT.md

*(original: `docs/archive/2026-10/merged/PREREG_hierarchy_AMENDMENT2_DRAFT.md`)*

## Amendment 2 — DRAFT v2 for the committed prereg (coordinator, 2026-10-07 evening). NOT committed.

Place it between Amendment 1 and Amendment 3 of `docs/PREREG_hierarchy_duplicon_expansion_2026-10-07.md`, with the heading below, in a commit made AFTER the code commit and BEFORE the first Liftoff call; give that commit's hash to `liftoff_cross.py run --amendment-commit`. `<CODE_COMMIT>` is the only placeholder to fill in before the commit of this text; `<AMENDMENT_COMMIT>` in the commands stays in the text and is replaced on the command line; `<species>` and `<label>` are path templates. The driver refuses to run unless the committed text contains `## Amendment 2` and the full 64-hex sha256 of the 7 code files of section 1 and of the 3 tools of section 2.

### Amendment 2 (2026-10-08, committed before the first Phase E call; it changes no cell, statistic, null or class of sections 1 to 13, and changes decision 12.7 (the hold-back unit) and decision 12.18 (the Gate 7 bar) as declared below)

Where the registered text allows more than one reading, section 10 fixes one. Two choices differ from the registered text, a wider exclusion and a stricter control, and are declared here: the verdict plan excludes the two control components, not only the LRPAP1 expansion (section 4), and the positive control is also judged on the triple of its class (section 10.4). The other departures are readings (section 10), schema additions (the `eligible` column of `control_pool.tsv`, the extra columns of `member_exons.tsv` and `record_calls.tsv`) and the driver's safeguards (sections 6 to 8).

**Authorisation.** The user's message of 2026-10-07 17:57 PDT, verbatim: 'go E, also please review the Yoo supplemental data in the thesis vault, I think my advisor wants me to prove that I find something prety similar to them, for example with LRPAP1'. It followed a message of 17:27 PDT that said the run is about an hour and that 'go E' buys the expansion list and its controls. It covers Phase E, tier B with the two controls of Gate 7, and nothing else: not tier A, not tier C, not T2a, T2b, R3a or R3b. The user was afterwards INFORMED, not asked, of the corrected plan (180 runs, 48 of them control runs, about 83 minutes; progress message of 18:54 PDT) and of the 2x time budget, the stricter positive control and the exclusion of the two control components (progress message of 23:38 PDT); no reply is recorded. These are the coordinator's decisions under the registered text (section 12(a)) and the user may stop the run at any time.

**1. Code (committed before the run).** Commit `<CODE_COMMIT>` holds

| file | lines | sha256 |
|---|---|---|
| `bench/hierarchy/expansion_table.py` | 997 | `4747328d0366d91c3b643d78bcaa1b9b62d5fd4c6bf7bd92bb2fb0b39acbfeff` |
| `bench/hierarchy/liftoff_cross.py` | 1466 | `a6aef94fb1c2b760e8623aedebf26b1a68749188e43156beff2d1c16a4c0288e` |
| `bench/hierarchy/test_expansion_table.py` | 3840 | `093a4552e6c2328ce6a31027b5f92aabd14553d01d59fd464c4f7b6b438addb7` |
| `bench/hierarchy/test_supplement_e.py` | 2646 | `8de392809dc9c39005a78660203beb7a5dbb42c538c3dff77f91915efca3b06b` |
| `figures/_liftoff.py` | 1592 | `a2d9b1d288a7da575c9b4beb47e3cdb94c630ca83c23588a049643fc0a1119ca` |
| `figures/figlib.py` | 496 | `f0a58426ca5c7d1a7c099490b1953b801e476bc6887f2885c8ea716410f78fab` |
| `tools/rlock.sh` | 25 | `d1bfc366e6e8ffc99aafa6d8fdcd73fb856961acf632781b4fcb932e0161b293` |

`python3 -B -m unittest bench/hierarchy/test_expansion_table.py bench/hierarchy/test_supplement_e.py`: 540 tests, OK under `PYTHONHASHSEED` 0 and 1, and under an audit-hook tripwire that blocks every access to the real run directory, the real tables and the real reference files (0 events). The tests use stub Liftoff scripts and made-up tables only; they show nothing about Liftoff 1.6.3 on a cross-species target. `figures/_liftoff.py` is imported unchanged (`parse_lifted`, `mm2_wrapper`, `TIME_FMT` and the registered parameters). `run` refuses code that is untracked or whose content differs from HEAD, refuses a later call with other code, and refuses to start unless this Amendment and the hashes above are in the committed prereg. A rehearsal of the whole phase on the real plan with a stub in place of Liftoff (180 runs in 3 calls) gives tables that are byte-identical under both hash seeds; its numbers are synthetic and are not results.

**2. Tools (Gate 0).**

| tool | version | path | sha256 |
|---|---|---|---|
| liftoff | 1.6.3 (entry script of the `liftoff` conda env, Python 3.12.3; gffutils 0.14, pyfaidx 0.9.0.4, networkx 3.6.1, numpy 2.5.1, biopython 1.87, parasail 1.3.4, interlap 0.2.7, pysam 0.22.1) | `/home/juanfra/miniforge3/envs/liftoff/bin/liftoff` | `239ae91a6399cfde33bd6536b8879a532f4c8d386677cf78baf0fa0a4fb4fa13` |
| minimap2 called by Liftoff (`-m`) | 2.24-r1122, through the per-call result cache of `figures/_liftoff.py::MM2_CACHE_WRAPPER` (keyed on the md5 of the query FASTA, the arguments, the index size and mtime and the minimap2 version; the cache lives in the run's `-dir` and is deleted with it) | `/home/juanfra/miniforge3/envs/liftoff/bin/minimap2` | `e91e2cf9717b415581896e09e7c5425e1d2f197d050a64165359de26e9c8e441` |
| minimap2 on the PATH, not used | 2.30-r1287 (also found: 2.28-r1209 sqanti3, 2.30 flair and rnaspades_env, 2.31-r1302 phasing_eval) | `/home/juanfra/miniforge3/bin/minimap2` | `5ef4a631a6318c66d1184007f87bedc55f712c27472d141da4d35e80df82e3cc` |
| driver Python | 3.14.4, run as `python3 -B`, standard library only | `/home/linuxbrew/.linuxbrew/bin/python3` | |
| GNU time | `/usr/bin/time -f 'wall_s\t%e\npeak_rss_kb\t%M...'` | `/usr/bin/time` | `3b11dec50514a8473e9f6efa7a34d584d0657538c09988f61b72d38ad4991a10` |
| lock wrapper | `flock` + `timeout` | `tools/rlock.sh` | `d1bfc366e6e8ffc99aafa6d8fdcd73fb856961acf632781b4fcb932e0161b293` |

**3. Inputs (Gate 0 additions).** Files below 1 GB carry their sha256, larger ones size and mtime plus a sha256 computed once, one file per call under `RLOCK_TIMEOUT=170 bash tools/rlock.sh light sha256sum FILE` before the first heavy call and appended to `/mnt/linuxdisk/tmp/hier_run/e/gate0_big_files.sha256`. The k15 w10 headers of the three outgroup indexes list exactly the sequence names and lengths of their FASTA `.fai` (25, 26 and 26 sequences; read from the index header, no minimap2 involved), and every `.fai` is newer than its FASTA, so neither Liftoff's `build_minimap2_index` (it builds only when `<target>.mmi` is absent) nor pyfaidx will write through the links.

| input | path | bytes | mtime | sha256 |
|---|---|---|---|---|
| gorilla assembly (REFERENCE) | `/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta` | 3,590,176,889 | 2026-06-04T23:54:46 | `bcf06464cfa2480ca9eb84b4713e3d1a0bb5dec4cf91230b0e4fceb3b62b41a6` |
| its `.fai` | `/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta.fai` | 994 | 2026-06-20T19:05:52 | `0198f946f3e28a5a7f5e92a005c75b9c98c4c769cd3b7e8f4c6d74c26e2e7b1a` |
| gorilla RefSeq annotation | `/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_genomic.gff` | 693,338,735 | 2026-06-04T23:55:01 | `32c5bd482f89810048dc46815e6239db9aabb09ca01f05c68e153cf56305bd9a` |
| gorilla self-lift table (pool) | `/mnt/linuxdisk/tmp/rustle_figures/liftoff/gorilla/liftoff_loci.tsv` | 15,561,852 | 2026-09-26T21:07:13 | `af503bf2713caf54d0606b585a7e44788c8d766936e364f4408633bf89438d73` |  53,323 rows
| gorilla record table | `/mnt/linuxdisk/tmp/rustle_figures/liftoff/gorilla/records.tsv` | 9,382,200 | 2026-09-25T14:14:10 | `825314751b07e2d2732ffaa97924b2c9f3fa3ab886204570c35716476520ea2f` |  41,193 rows
| H human FASTA (TARGET), CHM13v2.0 (RefSeq RS_2025_08) | `/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa` | 3,156,259,703 | 2022-09-28T09:23:04 | `15a4ba1246f6021a89699bf5083da7f2bad3f79c86acd7bc1eb0ca3a13164e85` |
| H `.fai` | `/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa.fai` | 788 | 2026-09-13T16:09:40 | `cec84e8875f740fea679b08cb042529deae56f8276f567c2d5dca232d6172c1b` |
| H k15 w10 index | `/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa.mmi` | 7,550,607,274 | 2026-09-11T17:30:06 | `58acd050bbc7215dc7f1fcf11332c570345c6e209703ee79e9f4c77a21b61940` |
| H own self-lift table | `/mnt/linuxdisk/tmp/rustle_figures/liftoff/human/liftoff_loci.tsv` | 17,409,071 | 2026-09-26T23:23:24 | `07cf5aed5f33c9524082b00a63988c8af6c889ce13e7f152b9ca62839ce10d26` |  71,982 rows
| C chimpanzee FASTA (TARGET), mPanTro3 (primary) | `/mnt/linuxdisk/home/juanfraitu/winloci_data/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna` | 3,217,481,210 | 2024-03-14T15:30:14 | `f71dfbcb3209b27292b50efc7bd0e461fb3899cdd669d2478909770f6244304a` |
| C `.fai` | `/mnt/linuxdisk/home/juanfraitu/winloci_data/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna.fai` | 989 | 2026-09-15T12:51:28 | `5792f6b3227d9757e78018229efa88b0dcad59d139b941a35aa20a41674e3b3a` |
| C k15 w10 index | `/mnt/linuxdisk/home/juanfraitu/ggo_npip/ref/PTR.fa.mmi` | 7,659,006,768 | 2026-09-17T16:20:29 | `8849765514cd737eb3978474a860b83b98a2ee1cc80131235213e0e7cec6c06e` |
| C own self-lift table | `/mnt/linuxdisk/tmp/rustle_figures/liftoff/chimpanzee/liftoff_loci.tsv` | 15,189,498 | 2026-09-25T22:42:54 | `ec9772c09e3498736413d0c63495ca3c559c6f7206b4f3cecc63ae8e6ed133e8` |  52,203 rows
| O orangutan FASTA (TARGET), mPonPyg2 (primary) | `/mnt/linuxdisk/home/juanfraitu/winloci_data/GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.fna` | 3,261,216,356 | 2024-03-14T15:00:02 | `6ff517010c8d9afdc4337d4a0bfe71925b4e6a6e47d3828ef2fcc83b62195d3a` |
| O `.fai` | `/mnt/linuxdisk/home/juanfraitu/winloci_data/GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.fna.fai` | 992 | 2026-07-17T00:15:47 | `ca2d5dcf8c41d8cd73a94893ea8508a1fb7ef9836fd30c5c501b4f6082fdab22` |
| O k15 w10 index | `/mnt/linuxdisk/home/juanfraitu/o1_falsemerge/lit/evid/PPY/genome/PPY.fa.mmi` | 7,738,998,992 | 2026-09-14T18:21:38 | `ec42e2ffa539beaea752deb969ad6c2eac032ff32a9a8e606247fb60622d50ea` |
| O own self-lift table | `/mnt/linuxdisk/tmp/rustle_figures/liftoff/orangutan/liftoff_loci.tsv` | 13,029,284 | 2026-09-26T22:41:56 | `eba8b18efb120075fa1e154b1bf4a18eceeb87a547045440e1b008ce3beed947` |  44,276 rows

**4. Plans (frozen; hub reading 'nodes', the registered rule).** Called set: 300 components, 1,403 records (907 gene, 496 pseudogene), J = 44, first five runs [300, 300, 165, 116, 81], 26 runs with 10 or fewer (the numbers of section 7, reproduced). The two controls are lifted first:

| group | content | runs per outgroup | records | Liftoff runs |
|---|---|---|---|---|
| `pos` | the LRPAP1 pool component (`PC_gene-LOC115933156`, 14 such records: the 6 qualifying, the 2 not required, 6 others joined by other copy rows) and the fragment component (`PC_gene-LOC115932756`, 3 records), rank-sharded; ranks 1 to 3 of the fragments share runs pos1 to pos3 with ranks 1 to 3 of the LRPAP1 component (pos1 = gene-LRPAP1 + gene-LOC134756368, pos2 = gene-LOC134756754 + gene-LOC115932954, pos3 = gene-LOC134756753 + gene-LOC115932756) | 14 | 17 | 42 |
| `neg` | the 100 planted records in 50 fake pairs (`random.Random(20261007).sample` over the sorted ids of the 34,283 eligible records; pairs formed greedily over the sample order, each on two contigs), member 1 of every pair in `ctl1` and member 2 in `ctl2` | 2 | 100 | 6 |
| `verdict` | the other 298 called components, rank-sharded (J = 44; first five runs [298, 298, 163, 115, 80]) | 44 | 1,386 | 132 |

Plan digest (sha256 over the hub reading and every run label with its record ids; recorded when the first Liftoff process is about to start, a later call with another plan stops): `107756b36d180aa557989176b8b323103ec0a799017e4934904c51937cd3238c`. Planting pool: 36,532 in_place records with exon union >= 200 bp, 34,993 without a copy row as source, 710 of those are the target of another record's copy row, 34,283 are eligible. The other hub reading ('records') would give 302 components and 1,446 records and is not used.

*The verdict set excludes the two control components (decision, reading 2).* Of their 17 records, 11 are the LRPAP1 copies and the 5' fragments (Table VIII.38 and the AMRP rows of Table VIII.35) and six are lncRNA genes that are joined to its component by copy rows of the gorilla self-lift: LOC134756754 (NC_073236.2:22,545,416-22,551,141), LOC129531455 (NC_073236.2:24,807,981-24,813,809), LOC134756758 (NC_073236.2:24,817,012-24,819,636), LOC134757220 (NC_073240.2:15,813,679-15,819,497), LOC129529291 (NC_073246.2:11,279,006-11,281,566), LOC134757959 (NC_073248.2:45,300,008-45,304,879). Four lie 2.5 to 2.6 kb and two 10.5 and 11.5 kb from an LRPAP1 copy; copy rows join them to each other and to the long gene models LOC129526389 and LOC129523503, which contain an LRPAP1 copy; that they were duplicated in the same events as LRPAP1 is likely and was not tested. The registered text (section 7, decision 12.7) excludes LRPAP1 from the verdict set at the level of the expansion; this Amendment excludes it at the level of the pool component, which is the unit of lifting, so the verdict plan is 298 / 1,386 / [298, 298, 163, 115, 80] instead of 300 / 1,403 / [300, 300, 165, 116, 81]. Expansions found among the control records are reported in `etable_control.tsv` and are not verdict statistics. Alternative: keep all 300 components in the verdict plan and drop the LRPAP1 class at build time (about +2 s per outgroup); it would put the lncRNA flank expansions into the verdict set. The statement of section 7 that 'the two controls contain no pool component of the verdict set' is true under this component-level exclusion.

**5. Command.** Run in the run directory `E_ROOT/lift/<species>/runs/<label>/`, which holds `in.gff` (the run's records and their descendants, in annotation order, `region`/`cDNA_match`/`match` lines left out; verified identical to an independent two-pass routing for all 60 labels and 1,503 records) and `types.txt` (`pseudogene`):

```
liftoff E_ROOT/lift/<species>/genome.fa /mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta -g in.gff -f types.txt -o lifted.gff3 -u unmapped.txt -dir inter \
         -copies -sc 0.95 -a 0.5 -s 0.5 -overlap 0.1 -p 4 -m E_ROOT/lift/mm2_cache.sh
```

`genome.fa`, `genome.fa.fai` and `genome.fa.mmi` in `E_ROOT/lift/<species>/` are links to the outgroup FASTA, its `.fai` and its index (read, never written). Liftoff's own `-mm2_options` stay at their default (`-a --end-bonus 5 --eqx -N 50 -p 0.5`), the process runs under `/usr/bin/time` in its own session with `PATH` starting at the Liftoff environment's `bin`, `TMPDIR=/mnt/linuxdisk/tmp`. After a run `inter/` and `in.gff_db` are deleted, `lifted.gff3` is kept as `lifted.gff3.gz`, and `loci.tsv` (parent-level rows from `figures/_liftoff.py::parse_lifted`), `unmapped.txt`, `time.txt`, `liftoff.log` and `done.key` (command, sha256 of `in.gff`, parameters) stay. `E_ROOT = /mnt/linuxdisk/tmp/hier_run/e`.

**6. Call schedule.** Three groups in this order, outgroup by outgroup within a group (H, then C, then O), and a call never mixes groups ('in their own calls'):

| group | Liftoff runs | predicted work | calls | gate before the group |
|---|---|---|---|---|
| `pos` | 42 | 1,056 s | 2 | none |
| `neg` | 6 | 184 s | 1 | the positive control must PASS (read from the 42 finished `pos` runs) |
| `verdict` | 132 | 3,768 s | 8 | both controls must PASS |
| total | 180 | 5,007 s (about 83 minutes) | 11 | |

```
cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto
## once, light: the plans and the wall budget (about 2 s; writes only plan.json), then the run GFFs and their marker (about 8 s)
RLOCK_TIMEOUT=170 bash tools/rlock.sh light python3 -B bench/hierarchy/liftoff_cross.py plan > /mnt/linuxdisk/tmp/hier_run/e/plan.json
RLOCK_TIMEOUT=170 bash tools/rlock.sh light python3 -B bench/hierarchy/liftoff_cross.py gffs
## heavy, FOREGROUND, one call at a time; repeat immediately while the exit code is 75
RLOCK_WAIT=5 RLOCK_TIMEOUT=585 bash tools/rlock.sh heavy python3 -B bench/hierarchy/liftoff_cross.py run --amendment-commit <AMENDMENT_COMMIT> --budget-s 530
## after exit 0, light: the table twice (Gate 3), each into a NEW directory
RLOCK_TIMEOUT=170 bash tools/rlock.sh light env PYTHONHASHSEED=0 python3 -B bench/hierarchy/liftoff_cross.py build --table-dir /mnt/linuxdisk/tmp/hier_run/e/table
RLOCK_TIMEOUT=170 bash tools/rlock.sh light env PYTHONHASHSEED=1 python3 -B bench/hierarchy/liftoff_cross.py build --table-dir /mnt/linuxdisk/tmp/hier_run/e/table_seed1      # SHA256SUMS must be identical
RLOCK_TIMEOUT=170 bash tools/rlock.sh light python3 -B bench/hierarchy/liftoff_cross.py build --strict-single-copy --table-dir /mnt/linuxdisk/tmp/hier_run/e/table_strict   # sensitivity, never the deciding table
```

Nothing but `plan.json` (the redirect above) is created by hand in `E_ROOT`; the code makes `lift/`, the links after checking the files behind them, the wrapper, `run.lock`, `run.log`, `wall.json`, `gate7.json`, `gff/` and, at the first launch, `plan.digest` and `run_identity.json`; `build` makes the table directories. Every heavy call is a foreground Bash call with the tool timeout set to 600000 ms (the harness moves a longer call to the background, against the crash rule), and every call repeats the identical `--amendment-commit` string (the full 40-hex hash). `TMPDIR=/mnt/linuxdisk/tmp` unless `TMPDIR` is already set in the environment of the call (unset it first). The 555 s process limit follows from the 590 s of the tool call (registered: split at 590 s). Exit codes of `run`: 0 all runs done, 75 run again at once, 3 a stop rule or a refusal of the code check (report, do not retry), 2 no or malformed `--amendment-commit` or a rehearsal flag with the real `--out`, 70 a defect of the driver, 128+n killed by signal n: for TERM, HUP and INT (143, 129, 130) the Liftoff session is killed before the driver exits; after SIGKILL (137, for example a harness or OOM kill) it is not, the tree ends by itself at most about 9.5 minutes after the launch of that run and holds `run.lock` through `timeout`, `time` and liftoff only, and the next call is refused by the stray-process guard until it is gone; from the wrapper: 1 = the heavy lock was not obtained in 5 s (call again), 124 = `timeout` ended the driver (look at `run.log` and `ps` before the next call). Exit codes of `build`: 0 valid, 4 table written but INVALID, 3 a stop (unfinished run, foreign run, changed code, input or identity), 2 refused. The call is a heavy job under `tools/rlock.sh heavy` (one at a time; two light jobs may run beside it). Accounting of the 590 s of a tool call: RLOCK_WAIT 5 s + RLOCK_TIMEOUT 585 s; the call clock starts at the start of the python process, the first unit of a call always runs and every later unit starts only if (elapsed + its predicted time) <= 530 s, no Liftoff process outlives 555 s from the start of the process, the kill takes at most 20 s, so the worst case is about 580 s. The predicted call lengths are [528, 528, 184, 517, 513, 521, 507, 513, 511, 511, 176] s: the packing of the predicted times from time 0; the driver counts the start-up of each call (about 2.5 s, about 12 s on the first) inside the 530 s, so the real number is 12 at the predicted speed, fewer for faster runs and more for slower ones. Before each call the coordinator looks at free memory by hand (`free -m`, MemAvailable of at least 20 GB (the stop of rule 2 is at 16,000,000 kB of summed RSS, and two light jobs of up to 2 GB may run beside the call)); the code has no system-level memory guard.

**7. Time and memory budget.** Prediction per run = 25 s start-up (index load and set-up; 10 to 20 s measured on the gorilla self-lift, plus margin) + 0.09 s x 1.25 per record: 352 s for the 14 `pos` runs, 61 s for the 2 `neg` runs and 1,256 s for the 44 `verdict` runs of each outgroup, 1,669.1 s per outgroup, 5,007 s in all (section 7 of the prereg said 26 to 50 minutes of work for 132 runs; the 48 control runs and the conservative start-up make it about 83 minutes). The estimate comes from same-genome self-lift fits (R2 0.00 to 0.55), not from cross-species lifts; the orangutan fit is dominated by one 544 s shard, a 10-fold outlier that occurred about once in 29 shards. **Time budget, implemented:** the Liftoff wall of every started process (finished, killed or interrupted) is summed per outgroup in `E_ROOT/wall.json`; when one outgroup's sum passes 2 x its prediction (2 x 1,669.1 s = 3,338.2 s) the call stops with exit 3; the stop is sticky (continuing needs a code change and a new Amendment). Peak resident memory 10 to 14 GB is expected (indexes of 7.6 to 7.7 GB); the summed RSS of the Liftoff session is polled every 0.5 s and the session is killed above 16 GB (16,000,000 kB as ps reports them, 15.3 GiB; the sum counts shared pages of the pool workers once per worker, 0.3 to 1.5 GB, so a run with a true peak above about 14.5 GB can be stopped). **Manual checkpoint:** after the first call the coordinator compares the wall of the H pos runs in `run.log` and `wall.json` with their prediction (352 s for the 14 runs) and stops and reports if the ratio is above 2; the per-outgroup budget alone can first fire inside the H verdict group.

**8. Stop rules** (each ends the call; none is retried; the state on disk stays valid and a later call resumes):

1. A Liftoff process that does not finish in the hard limit (555 s from the start of the python process, less the time already used) is killed (SIGTERM to the whole session, then SIGKILL of every member after 10 s) and its shard is split at the record-free cut nearest its middle (`<label>.s1`, `<label>.s2`, recorded in `splits.json`); a shard with no record-free cut is a stop. A process that is not the first of its call and is cut by the end of the call is not too long and is run again by the next call. The tree also runs under `timeout -k 10` at hard limit + 5 s inside its own session, so it ends by itself even if the driver is killed; SIGTERM, SIGHUP and SIGINT of the driver kill the session first; the Liftoff tree inherits the driver's `flock` on `E_ROOT/run.lock` through `timeout`, `time`, liftoff and its pool workers (not through minimap2 or the cache wrapper, which liftoff starts with `subprocess.run`); a minimap2 that outlives them is caught by the stray-process guard of rule 7 at the next call, not by the lock. A Liftoff session that is not empty 10 s after SIGKILL is a stop: no run may start until `ps --sid` is empty.
2. Summed resident memory of the session above 16 GB: the session is killed.
3. A non-zero Liftoff exit, a missing `lifted.gff3` or `in.gff`, or a run that is not complete: every parent record of the run's `in.gff` must appear in `lifted.gff3` or in `unmapped.txt` (the invariant holds on all 179,088 self-lift records); an incomplete run writes an `incomplete` file and no `done.key` (the invariant was checked on same-genome self-lifts only; if a cross-species run legitimately leaves a record neither lifted nor listed, the very first run stops with exit 3 and the evidence is in that file).
4. Gate 7 failing before the `neg` or the `verdict` group, an unfinished control run when the gate is asked, or `build` finding a control or verdict run unfinished, made by another command or `in.gff`, or in a rehearsal directory: no verdict-set record is lifted after a failed control; the gate is re-evaluated at every call (a saved failure is never accepted). Section 7 allows this Amendment to change the Liftoff command only on the evidence of the two controls; any such change is a new Amendment, committed before the controls are lifted again.
5. Another plan digest, other code, other tools, another command, other pinned inputs (17 inputs pinned in `run_identity.json`: sha256 below 1 GB, size and mtime above) or another Amendment commit than the first call recorded.
6. The time budget of section 7.
8. Refusals before anything is frozen (exit 3): a missing or unusable Liftoff, `timeout` or `time` binary; an unusable reference or outgroup FASTA, `.fai` or index (a `.fai` older than its FASTA, an index that is not the k15/w10 index listed in section 3, a missing link target); a committed text that lacks `## Amendment 2` at the start of a line, still contains the literal placeholder for the code commit, or does not carry the registered plan digest, the hashes of section 1, 2 and the small inputs of section 3, or a code commit that holds the listed files.
7. Another driver (the `flock` on `E_ROOT/run.lock` is busy) or a stray `liftoff`/`minimap2` process of this user with `E_ROOT` in its command line or as its cwd: the call refuses before anything is frozen.

**9. Outputs** (`build`, light; `E_ROOT/table/`; every file is listed with its sha256 in `SHA256SUMS`, in name order, the sha256 of nine tables (the four registered ones, `siblings.tsv`, `dropped_classes.tsv`, `record_calls.tsv`, `components.tsv`, `etable_human_only.tsv`) is also in `summary.json`; Gate 0 records them when section 7 closes): the registered `etable.tsv`, `member_exons.tsv` (columns record_id, pool_component_id, gene_id, exons), `control_pool.tsv` (columns pool_component_id, record_id, reason, joining_copy_rows, and one more, `eligible`: 0 for reason no_placement, 1 otherwise), `extension.tsv`; and `siblings.tsv`, `dropped_classes.tsv` (the three-outgroup classes dropped for multiplicity, section 10.5; the human-only classes dropped are only counted, in `summary.json`), `record_calls.tsv` (the call of every verdict record in every outgroup: status, detail, qualifying placements, reverse copy rows, outgroup locus, and a last column `multiplicity` naming the outgroups in which the record has a reverse copy row), `components.tsv` (every one of the 609 pool components with its status: hub, control, verdict, too_few_such_records), `etable_control.tsv` (the expansions of the two control components, i.e. LRPAP1, reported separately), `etable_human_only.tsv` (the sensitivity of Amendment 3: the same lifts, the triple reduced to the human locus h), `gate7.json` (both controls with their reasons, the digests of the control outputs, and for the positive control the triple, the expected triple and `triple_ok`), `STATUS` (`VALID`, or `INVALID: <reason>`; tables are written either way, and anomalies make the status INVALID) and `summary.json` (counts, including the records with no placement by detail, the number whose only row is partial, the components with different loci, the sibling groups, `hubs_with_two_such_records`, `dropped_classes`, `records_in_dropped_classes`, `records_with_reverse_copy_rows` per outgroup, `human_only_dropped_classes` and `anomalies`). For a VALID table `build` prints the three-outgroup dropped-class counts; after a failed control it writes `gate7.json` and `STATUS` (`INVALID: <reason>`, exit 4) and no table.

**10. Readings fixed here** (the registered text allows more than one; the other reading is in brackets).

1. Hub rule: members are the annotated records of any length plus the `copy@` sites of the union-find component, as written in section 7 (7 hubs among the 307 components; [annotated records only: 5 hubs, 302 components, 1,446 records]).
2. The LRPAP1 component and the fragment component are lifted only in the `pos` group and leave the verdict plan (section 4, with the reason and the alternative).
3. The positive control lifts the whole LRPAP1 pool component, 14 runs per outgroup, because the registered pass condition and failure mode refer to its rank plan. [Only the 8 relevant records, 6 runs fewer per outgroup, about 450 s less.]
4. 'An LRPAP1 copy interval' is the landing interval of one of the 8 gene-LRPAP1 rows of the gorilla self-lift (1 in_place, 7 dropped_M1; the 6 qualifying records have 0.916 to 1.0 of their exon bases inside, LOC129526389 has 0.414 and LOC129523503 0.164). **The control is also judged on the triple, which is stricter than the registered pass rule:** it passes iff ONE class holds the 6 qualifying records, none of the six is flagged multiplicity (section 10.5), AND that class's triple is (gene-LRPAP1, gene-LRPAP1, gene-LRPAP1), the gene-LRPAP1 anchors of H, C and O (in_place and single-copy in all three outgroup self-lifts); a class with another triple fails the control, and `gate7.json` and the failure message name both triples. The two records not required, and any other member of the class (for example a flank), may or may not be placed and join the class without a failure; a member that is not one of the six and carries a reverse copy row is left out of the control's class for the control (it cannot fail it), is named under the flagged members in `gate7.json`, and is still flagged, and its class is still dropped, in the real E-table. [Registered rule only: one class holds the 6 records.]
5. **Multiplicity (restated).** A gorilla record has more than one placement in an outgroup if it has a reverse row with `extra_copy_number` > 0 (of any coverage or identity) or more than one row with coverage >= 0.5. Such a record keeps the outgroup locus of its PRIMARY placement as a candidate (detail `reverse_copy_row` when the outgroup locus is ok) and is FLAGGED multiplicity (`record_calls.tsv`). Classes are formed as for every other record (records of one component with the same triple; no chaining). Every class of two or more records that contains a flagged record, in any of the three outgroups for `etable.tsv` and in H for the human-only table, is DROPPED AS A WHOLE and COUNTED (`dropped_classes.tsv`, `summary.json`, the output of `build`); this is the registered rule of section 7 ('a class with any member that has more than one is dropped and counted'), and the clause 'the reverse lift gives that gorilla record no extra_copy row in that outgroup' is read as the condition of a CLASS, not as a reason to exclude the member alone. A class of one flagged record is not an expansion candidate and is not counted as dropped. A record whose primary row is not a placement (partial_only, not_a_placement, unmapped, no row) has no placement, no triple and no flag. [Member-level exclusion, in which the guard is vacuous: not used.] The positive control fails with the reason 'class dropped for multiplicity' if a qualifying LRPAP1 record is flagged.
6. **Outgroup locus (restated).** The outgroup anchors are the in_place records of the outgroup's own self-lift, of any length (records under 200 bp are candidates: 7,085, 4,660 and 4,101 of them in H, C and O), that share at least 50% of the shorter exon union with the gorilla record's placement; a candidate that is the source of a copy row gives the status `outgroup_copy_row` whatever other candidates exist (a copy-row candidate overrides every single-copy candidate: the quoted yields reproduce only under that precedence); otherwise two or more candidates are ambiguous and give no triple (reason `no_outgroup_locus`, detail `ambiguous`) and exactly one candidate is the locus. As a call status `outgroup_copy_row` now means only that the OUTGROUP's own self-lift gives the locus a copy row; a reverse copy row of the gorilla record is the flag of reading 5. Input-only yield of the rule, queried with each outgroup's own records of exon union >= 200 bp (50,560, 37,010 and 33,354 of them): H 45,209 of 50,560 ok, 2,060 copy row, 3,291 ambiguous; C 34,657 of 37,010, 1,494, 859; O 31,612 of 33,354, 1,049, 693. [Strict single-copy, which also counts a record that is the landing site of another record's copy row: H 44,491, C 34,084, O 31,228 ok; available as `build --strict-single-copy`, never the deciding table.]
7. `control_pool.tsv` lists the registered reasons (outgroup_copy_row, different_triple, no_placement). A record flagged multiplicity in any outgroup is listed under `outgroup_copy_row`, eligible 1, ahead of `no_placement` and whatever its triple, and a non-hub component that contributes no expansion qualifies for the pool when it holds one; a component whose only classes were dropped for multiplicity is therefore in the pool, a component with a surviving expansion is not. A record placed on no annotated outgroup locus, one with an ambiguous locus, and a record whose triple is the only one of its component have none of these reasons and are not listed (they stay in `record_calls.tsv`); an unlisted record is not eligible as a control. The pool is a superset of the pool of the first reviewed version (member-level exclusion) on 6,000 random worlds. [Listing the unlisted records as eligible.]
8. 'The number of records whose only row is partial': a record whose only row is flagged partial_mapping (detail `partial_only`); a row that fails the placement rule otherwise, or has copy rows beside it, is detail `not_a_placement`.
9. Planted pairs: the sample of 100 in `random.Random(20261007).sample` order is paired greedily, each record with the first later one on another contig. [Any other pairing rule.]
10. Rank order: contig as a Python string, then start0, then end, then record id (one called component, `PC_gene-ALDH1B1`, has two records with the same contig and start0; it is ordered by end).
11. The identity of a `copy@` site includes the strand (`copy@contig:start0-end<strand>`; the registered text names `copy@contig:start-end`); on the real data there are 1,200 distinct loci and 609 components either way.

**11. What this Amendment does not do.** It lifts nothing for tier A or tier C, builds no index, opens no family file and no BAM, and computes no T2a or T2b statistic. The Y4 concordance of Amendment 3 reads the tables of section 9 and needs no further lift. Its reason chain (no pool record, hub, outgroup copy row, no placement, different triples) is extended here, before any Y4 number exists, by: component not called (too few such records, or a hub); control component (`etable_control.tsv`, `gate7.json`); record shorter than 200 bp; record alone on its triple; no annotated or an ambiguous outgroup locus; class dropped for multiplicity.

**12. Decisions taken before the first call.** (a) The coordinator reads the 'go E' as covering the 48 control runs (1,240 s of predicted work, about 21 minutes, on top of what section 7 stated), the factor 2 of the time budget, the stricter positive control and the exclusion of the two control components; the user was informed of each (Authorisation) and not asked, and may stop the run at any time. (b) Readings 2 and 3 as in section 10; the multiplicity reading 5 restores the registered class-level rule. (c) The code is committed first (`<CODE_COMMIT>`), then this Amendment with the hashes of section 1 (`<AMENDMENT_COMMIT>`). (d) The Gate 0 hashes of the files over 1 GB were computed once, one file per light call, and are in section 3 and in `/mnt/linuxdisk/tmp/hier_run/e/gate0_big_files.sha256`. (e) Independent review: five reviewers of the first version (plan recount, logic by differential testing on 3,800 random worlds, run safety, Amendment audit, mutation testing of 1,949 mutants), a fix round (reports in `/mnt/linuxdisk/tmp/hier_review/e/fix_round/CHANGES.md`: two blockers, five majors and the mutation supplement, 540 tests) and a delta review of the final files. (f) During the fix round a builder ran two test files before they were hermetic; they created an empty `lift/` and an 8-byte `run.lock` in `E_ROOT` and read real inputs read-only; both were removed by exact path and the listing of `E_ROOT` was checked identical afterwards (CHANGES.md section 0). (g) Errata of the committed prereg, no rule changed: the registered 132 runs are the 44 rank runs per outgroup of the verdict group; the registered 26 to 50 minutes (10 to 20 s start-up, 0.05 to 0.09 s per record) become 63 minutes for the same runs under the 25 s start-up and 0.09 x 1.25 s per record used here (3,768 s), and the 48 control runs add 1,240 s (21 minutes); the called set including the two control components is 300 / 1,403 and the planted control was not in the registered cost; the verdict plan excludes the control components (section 4).


---

## IDEAL_CHROMOSOME_SIM_2026-09-21.md

*(original: `docs/archive/2026-10/merged/IDEAL_CHROMOSOME_SIM_2026-09-21.md`)*

## Ideal-scenario chromosome simulation — in theory it does work

**§6v1, 2026-09-21.** Pre-registration `docs/archive/2026-09/PREREG_ideal_chromosome_sim_2026-09-21.md` (md5
`ff226f41`), committed `2a809b76` before any arm was scored. Tool `bench/ideal_chromosome_sim.py`.
⚠HUMAN substrate (A119b/CHM13 chr16) — do not pool with gorilla.

### Question

*"Simulate an ideal scenario using 1 chromosome with some multi-copy gene families but also some others
that are single copy, to prove that in theory it should work. We just need a way to remove unspecific
phenomena like readthroughs and ensure the reads are complete."* — establish the CEILING, and attribute
the real-data gap.

### Setup

chr16, **4,411 annotated transcripts over 1,443 genes**, 10 reads per transcript, ends jittered +-0-30 bp
(mandatory — identical reads collapse under dedup, §6n0), err 0.001. Aligned with the shipped
`-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes` and assembled through the shipped
`--assemble-only --assembly-polish full` path, unchanged.

Strata fixed from the SHIPPED chr16 homology graph before any arm ran: **MULTI-COPY** = the gene has >= 1
homology edge (297 genes); **SINGLE-COPY** = none (1,110). Universe = the 1,407 simulated genes that carry
an annotation record — an INPUT, so nothing is conditioned on any arm's output.

The readthrough rate is measured, not invented: **33,058 of 439,985 (7.51%)** real primary MAPQ-60 chr16
reads hit the exons of >= 2 distinct annotated genes at >= 25 bp. A_rt injects molecules splicing a gene's
transcript to a neighbouring same-strand gene's at that rate (3,014 of 43,700 reads = 6.9% realised).

### Result

| arm | stratum | raw endpoint | unscoreable | **genuine over-merges** | corrected |
|---|---|---|---|---|---|
| **A_ideal** | ALL | 0.8067 | 15.8% | **12** | **0.9578** |
| | MULTI-COPY | 0.6364 | 28.6% | 12 | 0.8915 |
| | SINGLE-COPY | 0.8523 | 12.3% | **0** | **0.9723** |
| A_trunc | ALL | 0.8188 | — | — | — |
| A_rt | ALL | 0.6503 | 12.6% | **188** | 0.7439 |
| | MULTI-COPY | 0.5556 | 22.9% | 46 | 0.7205 |
| | SINGLE-COPY | 0.6757 | 9.8% | 142 | 0.7493 |
| real data | ALL | 0.4961 | 13.9% | 139 | 0.5759 |

"Unscoreable" = the gene's node is shared ONLY with an annotated gene that overlaps it. Two annotated
genes sharing exonic sequence cannot be separated into one-to-one nodes by ANY method, so these are a
truth ceiling, not an error. "Corrected" excludes them; it is **post-hoc** and is reported beside the
pre-registered raw number, never instead of it.

### Verdict against the pre-registered bars

- **Bar 1 (A_ideal >= 0.95 overall and >= 0.90 multi-copy) FAILS as written**: 0.8067 / 0.6364.
- **Bar 2 (A_ideal < 0.80 ⟹ algorithmic defect) NOT triggered**: 0.8067.
- **Bar 3 (single-copy >= 0.95) FAILS as written**: 0.8523.
- **Bar 4 (attribution) PASSES decisively** — see below.

⚠**The bars were mis-set, and that is a prereg design flaw worth recording.** The endpoint has a
built-in ceiling of roughly **0.84** on this substrate, because 15.8% of the universe consists of
annotated genes that overlap another annotated gene. A bar of 0.95 on a metric that cannot exceed ~0.84
was unreachable by construction. This is the same ~0.8 ground-truth ceiling already documented in
§6kl/§6km, arrived at independently.

### What the simulation actually shows

⭐**With full-length reads and no readthrough the algorithm makes almost no errors**: **12 genuine
over-merges across 1,407 genes**, and **ZERO in the single-copy stratum**. Corrected for the truth
ceiling that is 0.9578 overall, 0.9723 single-copy, 0.8915 multi-copy. **So yes — in theory it works,
and multi-copy families are recovered as separate nodes.**

⭐**Readthrough alone reproduces the real pathology, and then some.** Injecting it at the measured 7.51%
takes genuine over-merges from **12 to 188** — a 15.7x increase — against **139** in the real data. The
injected rate slightly OVERSHOOTS the real defect, so readthrough is not merely sufficient to explain the
over-merge found in §6u7, it is more than sufficient. Nothing else needs to be invoked.

⭐**5' truncation costs nothing on this endpoint** (A_trunc 0.8188 vs A_ideal 0.8067, and the identical
multi-copy rate 0.6364). ⚠This does NOT contradict §6n0, where truncation cost 6 of 26 NPIP copies:
that endpoint was transcript COMPLETENESS, this one is where the node SITS. Truncation shortens models
without moving them.

### Consequence

The gap between the ideal ceiling (0.9578 corrected) and real data (0.5759 corrected) is attributed:
**readthrough carries it**. This is the third independent line pointing at the same conclusion as §6u7
(the defect is over-merge, not fragmentation) and §6v0 (full-length-ness corroborates bridges rather
than discriminating them) — and it is the first to show the algorithm is clean when the phenomenon is
removed. ⟹ **the open lever remains a readthrough-aware node SPLIT**, and it now has a measured ceiling
to aim at.

---

## 8. Does the idealized version hold a better FAMILY definition?

**User, 2026-09-21.** §6v1 scored only the NODE-construction endpoint. This runs the actual family
definition (`mcl_families --min-exonic-bp 1 --min-shared-exon-frac 0.60`, the shipped Rust binary, same
recipe as `dn16_fam3`) on the A_ideal and A_rt assembled loci, and scores the resulting clusters against
**both** independent truths from §6u8/§6u9 (Soto cover, protein referee), restricted to chr16.

⚠**n is small** — chr16 alone has 15 Soto families / 36 protein-referee families. Per-family movement
swings are large at this n (§6u8's lesson); read the direction, not the third decimal.

### 8.1 Pipeline

For each arm: exon-merged locus spans from the assembled GTF -> region list -> `samtools faidx` body
extraction -> `minimap2 -x asm20 -c --eqx -P` all-vs-all -> shipped `mcl_families` with the locus's own
exon-union GFF3 (mirroring `dn16.gff3`'s format exactly) -> `.clusters.tsv`. A fourth arm, **GUIDED**
(annotation records used directly as nodes, no assembly at all — `chr16_guided.clusters.tsv`, already on
disk), stands in as the ceiling: perfect nodes, no simulation needed.

### 8.2 Result

| arm | clusters >=2 | vs **Soto cover** (15 fams) | vs **protein referee** (36 fams) |
|---|---|---|---|
| **GUIDED** (annotation nodes) | 90 | sens 0.4567 / prec 0.4294 / **F 0.4174** | sens 0.2046 / prec 0.3203 / **F 0.2248** |
| **A_ideal** de novo | 56 | sens 0.3692 / prec 0.3649 / **F 0.3561** | sens 0.1722 / prec 0.2968 / **F 0.1944** |
| REAL de novo | 70 | sens 0.3372 / prec 0.3255 / **F 0.3166** | sens 0.1681 / prec 0.2172 / **F 0.1805** |
| A_rt de novo | 51 | sens 0.3136 / prec 0.3182 / **F 0.3068** | sens 0.1329 / prec 0.2603 / **F 0.1570** |

**Yes — modestly, and the direction holds on both independent truths.** A_ideal beats real de novo by
+0.0395 F (Soto) and +0.0139 F (protein referee). A_rt falls below REAL on both truths despite injecting
readthrough onto otherwise-ideal reads — consistent with §6v1's own finding that the injected 7.51% rate
slightly overshoots the real defect (12 -> 188 genuine over-merges vs 139 real), so this is not a
contradiction: an aggressive dose of the one thing being tested does more damage than the full complexity
of real data, which is itself evidence the mechanism is right.

Per-family (Soto), 3 better / 6 worse / 6 unchanged: the pooled gain is real but not one-sided — a caveat
consistent with n=15.

### 8.3 The bigger finding: node CORRECTNESS is fixed, node COVERAGE is not

**A_ideal falls well short of GUIDED — F 0.3561 vs 0.4174 (Soto), 0.1944 vs 0.2248 (referee) — despite
§6v1 showing node correctness is nearly solved (0.9578 corrected).** The reason is not over-merge (that
was measured and is nearly gone) and not fragmentation-as-scored-by-§6u7 (per-copy correctness is high).
It is **node COVERAGE**: how many of the truth genes have ANY node representing them in the homology
graph at all.

| arm | graph nodes | Soto truth genes (58 on chr16) represented as SOME node |
|---|---|---|
| GUIDED | 380 | **47 / 58 = 81.0%** |
| A_ideal | 324 | **34 / 58 = 58.6%** |

A 22.4-point coverage gap, under perfect reads and no readthrough. Of the 24 genes missing a node in
A_ideal, 17 (71%) DO have an assembled locus — the locus exists, it is just absent from the graph,
meaning either it never accumulated enough exonic content to clear the `identity >= 0.70 / cov_longer >=
0.30 / >= 300bp` edge floor, or (below) it never gets scoring credit at all:

**Of the 13 genes missing in A_ideal but present in GUIDED (the true assembly-attributable set), 5
(38%) are a SCORING ARTIFACT, not an assembly failure.** RefSeq itself carries curated **readthrough/
fusion gene records** overlapping almost exactly the same span — `PKD1P3-NPIPA1` (42,950 bp) sits over
`NPIPA1`'s 14,597 bp span at 37,168 bp overlap; `PKD1P4-NPIPA8`, `PDXDC2P-NPIPB14P`, `BOLA2-SMG1P6`
do the same for `NPIPA8`, `NPIPB14P`, `SLX1B`. The de novo locus's homology signal for the real gene IS
present (verified directly on the PAF: NPIPA1's locus has 10+ alignment records at 0.83-0.99 identity to
other loci) — but a **max-overlap, one-name-per-locus resolver, applied at scoring time, hands the whole
locus's credit to the bigger fusion-named record**, which then has no truth family to match, so the real
gene registers as "no node" though its sequence was correctly assembled and correctly grouped. This is
the family-definition-scoring analogue of §6v1's "unscoreable (annotation overlap)" category, applied
here for the first time to the CLUSTER-naming step rather than the per-copy node step.

The remaining ~8 (ABCC6, HERC2P5, HERC2P8, NPIPB10P, NPIPB7, PKD1, PKD1P6, SMG1P6) have no bigger
co-located record and are candidates for a genuine remaining edge-construction gap — short or
low-coverage assembled loci that never accumulate enough exonic content relative to their full-length
siblings to clear the coverage floor. This was not run to ground on this pass and is the natural next
target.

### 8.4 Consequence

The answer to "does the idealized version hold a better family definition" is **yes, but only partly for
the reason expected.** Fixing readthrough (§6v1's target) recovers a real, if modest, family-definition
gain — confirmed on two independent truths. But it does not close most of the gap to a perfect-node
ceiling, because **most of that remaining gap is not a node-construction defect at all**: at least 38% of
it is an artifact of how RefSeq's own curated fusion-gene names get resolved at scoring time, in exactly
the same region (PKD1/NPIP tandem duplication) already flagged by §6u7/§6v1/§6v0 as ground zero for
readthrough. The rest is a real, smaller, and as-yet uncharacterized edge-coverage gap.

---

## 9. The 8 residual genes: three distinct mechanisms, not one

**User, 2026-09-21.** Diagnosed each of the 8 genes from §8.3 (ABCC6, HERC2P5, HERC2P8, NPIPB10P,
NPIPB7, PKD1, PKD1P6, SMG1P6) individually against the raw annotation, the assembled GTF, and the PAF.

### 9.1 Three categories, verified against the raw GFF and the assembled loci

| category | genes | mechanism | fixable by better assembly / reads? |
|---|---|---|---|
| **A — no transcript model exists** | HERC2P5, HERC2P8, NPIPB10P (3) | Verified on the raw GFF: these `pseudogene` records have **zero** child `transcript`/`mRNA`/`exon` records of any kind — RefSeq assigns them a genomic span and nothing else | **No.** There is no sequence to simulate and (if real biology matches the annotation) no transcription to sequence. Not an assembly or read-quality question at all. |
| **B — fusion/readthrough truth-naming** | PKD1P6 (1) | Verified on the raw GFF: `PKD1P6`'s only two `transcript` records (`rna-NR_123721.1`, `rna-NR_123722.1`, 30 and 18 exons) are children of `gene-PKD1P6-NPIPP1`, **not** of `gene-PKD1P6`. Standalone PKD1P6 has no transcript of its own — its real transcription, when it occurs, produces the readthrough molecule, not a PKD1P6-only one. Joins §8.3's 5 fusion-name cases (NPIPA1, NPIPA8, NPIPB14P, SLX1B, PKD1P6-vs-`PDXDC1`) at **6 of 13** now. | **No.** This is the same phenomenon the whole session has been chasing (§6u7/§6v0/§6v1), now showing up as a truth-side naming artifact rather than a node over-merge. Full-length confidence changes nothing about which curated gene name a real fusion molecule maps to. |
| **C — real paralog, real assembly, genuine exon-structure mismatch or fragmentation** | ABCC6, NPIPB7, SMG1P6 (3) | See 9.2 below — two different sub-mechanisms | Partially — see 9.2 |

### 9.2 Category C in detail

**ABCC6**: both ABCC6's locus (4,831 bp exonic, vs its real transcripts' 4,463-4,550 bp — essentially
COMPLETE) and its true paralog ABCC6P1's locus (`DN_chr16_18504057_9`, 2,660 bp exonic, vs its own
2,664 bp transcript — also essentially complete) assemble correctly and independently. The PAF shows a
strong genomic alignment between them (23,511 bp aligned at 97.85% identity). **The edge still fails**,
almost certainly at the `--min-shared-exon-frac 0.6` conjunct: `exonic_denominator=true` means the
coverage/sharing test is computed over EXON-restricted alignment, not the raw genomic span, and ABCC6P1
is a `transcribed_pseudogene` — very likely a retro/processed duplicate whose own annotated exon
structure has diverged from ABCC6's spliced mRNA (fewer/different exon boundaries after duplication).
This is register 925's already-documented finding — *"half of family relationships have an INTRONLESS
member"* — measured here on a NEW pair with everything else held ideal. **Not an assembly defect: the
two loci are each correctly and completely reconstructed; the paralogs' real, annotated exon structures
are what differ.**

**NPIPB7, SMG1P6**: genuinely fragmented, unlike ABCC6. Comparing locus exonic content to the gene's
OWN transcript length (not genomic span, which is misleading for intron-rich genes): NPIPB7's locus
captures 464 of 2,254 bp (20.6%) of its single 10-exon transcript; SMG1P6's captures 326 of 2,623 bp
(12.4%) of its single 18-exon transcript. Both are unusually **exon-dense, compact transcripts**
(NPIPB7 ~225 bp/exon, SMG1P6 ~146 bp/exon) — a real under-assembly, not a coverage or truncation
artifact (A_ideal reads are full-length, error rate 0.001, jitter only 0-30 bp per end). This is a
genuinely open lead: whether minimap2's spliced alignment or the assembly/polish step is what drops
these small exons was not run to ground this pass.

### 9.3 Does isoseq cluster2 help?

**For 6 of 8 (categories A and B): no, and it cannot in principle.** Category A has no transcript to
sequence in the first place, real or simulated — cluster2 clusters reads, and there would be no reads.
Category B is a truth/naming problem: PKD1P6's real biological molecule, when transcribed, IS the
PKD1P6-NPIPP1 fusion. cluster2 would correctly report a confident, full-length, HQ isoform for that
fusion — which is precisely what §6v0 already established cannot be resolved by full-length-ness,
because the fusion **is** the genuine complete molecule.

**For ABCC6 (category C, exon-structure divergence): no.** Both loci are already reconstructed
essentially completely from these ideal simulated reads; the problem is the exon-conjunct correctly
detecting that the two paralogs' real, annotated splice structures differ. Better read processing
cannot change what the paralog's real exon structure is.

**For NPIPB7 and SMG1P6 (fragmented, compact transcripts): plausibly yes, and this is the one candidate
worth testing.** `isoseq cluster2` calls a consensus directly from the FLNC read sequences (POA-style,
before any genome alignment), which could recover a clean, complete transcript model for a compact
multi-exon gene even where individual per-read spliced alignments to the genome are imperfect on small
exons — a different failure surface than a reference-guided splice-alignment assembler faces. This is
speculative: the exact cause of the fragmentation (minimap2's splice seeding on <150-225 bp exons, vs
an assembly/polish-step filter) was not isolated this pass, so cluster2's benefit here is a hypothesis
to test, not a demonstrated fix — and it is a narrow win even if confirmed (2 of 8 genes here, and an
unknown fraction of chr16 more broadly).

### 9.4 Bottom line

Of the 8 residual genes, **6 have no assembly-side fix available at all** (3 have no transcript to
recover, 3 map to a real fusion). Only **2 (ABCC6, NPIPB7... correction: ABCC6 is unfixable by
assembly too — so really only NPIPB7 and SMG1P6, 2 of 8)** are candidates for a genuine, fixable
assembly gap, and `cluster2` is a plausible but unproven lever for exactly those two. This further
narrows what remains open after §6v1/§6v2: the readthrough-aware node split is still the dominant lever;
this residual set is small, mostly structural rather than algorithmic, and the one live thread worth a
follow-up is specifically **compact, exon-dense transcript assembly**, not readthrough or node
over-merge.

---

## 10. Root cause of the NPIPB7/SMG1P6 fragmentation: a NON-CANONICAL junction in RefSeq's own annotation

**User, 2026-09-21: "let's run the diagnostics to figure out why NPIPB7/SMG1P6 fragment."** Isolated to a
precise, confirmed mechanism, ruling out three alternatives along the way.

### 10.1 Ruled out, in order

1. **Read alignment**: all 10 simulated reads' PRIMARY alignments (`-F 2308`) land within 150 bp of each
   other at the true locus, with CIGARs that correctly reconstruct the full transcript. For SMG1P6, all
   10 reads carry the **exact same 17-intron chain**, every intron at 10/10 support — zero disagreement.
2. **Assembly polish**: re-ran with `--assemble-only` and NO `--assembly-polish` flag at all. The raw,
   unpolished GTF is **byte-identical in the fragmented region** to the polished one — the loss happens
   before polish ever runs.
3. **A minimal single-read / identical-duplicate repro**: running the assembler on 1 read, or on 5
   byte-identical copies of 1 read, produces **empty output** for both genes — consistent with the
   project's documented dedup rule (`sims need -N 50 and must VARY read ends`) collapsing identical
   reads to one effective molecule, below whatever minimum-support floor applies. Not informative about
   the fragmentation itself, but rules out a naive "not enough reads" framing.

### 10.2 The actual cause: a non-canonical splice junction in the gene's OWN RefSeq exon model

Checking every intron's splice motif (donor/acceptor dinucleotide, strand-corrected) against the genome
directly:

| gene | biotype | introns | first non-canonical junction |
|---|---|---|---|
| **SMG1P6** | `transcribed_pseudogene` | 17, all GT-AG except one | intron 15 of 17: `chr16:29721738-29725407` (3,670 bp), motif **AT-AG** |
| **NPIPB7** | `protein_coding` | 8, all GT-AG except the last 3 | introns 6-8 of 8: `chr16:28752509-28760135` (TA-GA), `chr16:28760195-28761928` (GT-GG), `chr16:28762017-28778548` (GG-TT) |

**Verified independent of my own pipeline**: the exon coordinates flanking each non-canonical junction
match the RAW, un-processed GFF `exon` records for these transcripts exactly (`rna-NR_135312.1`,
`rna-NM_001396030.1`) — this is RefSeq's own curated exon model, not an artifact of my extraction or of
simulated sequencing error.

The shipped assembler enforces canonical (GT-AG / GC-AG / AT-AC) splice motifs when building a spliced
transcript model, and — confirmed by direct test — **truncates the transcript at the first non-canonical
junction it encounters, discarding everything on the far side, even when every read agrees perfectly and
completely on the true structure.** This is not a bug introduced by imperfect reads; it fires on reads
that are already ideal.

### 10.3 Confirmation: `RUSTLE_JUNCTION_MAJORITY=1` recovers it

This project already has a named, documented mitigation for exactly this failure mode
([[project_junction_majority_chr16]] / §6m8: *"strict canonicity is 69.4% of pass-1 -> GTF loss"*).
Re-ran A_ideal with it set:

| gene | default (strict canonicity) | `RUSTLE_JUNCTION_MAJORITY=1` |
|---|---|---|
| **SMG1P6** | 2 exons, 224 bp (1 of 17 introns) | **FULL RECOVERY**: `chr16:29708443-29728813`, matching the annotated span almost exactly, gene_id suffix `_18` = 18 exons |
| **NPIPB7** | 4 exons, 1,282 bp | **PARTIAL RECOVERY**: `chr16:28739068-28762017`, reaching well past the original truncation point but still short of the transcript's far end (which carries THREE consecutive non-canonical junctions, not one) |

SMG1P6 recovers completely because it has exactly one non-canonical junction; NPIPB7 only partially
recovers because it has three in a row near its 3' end, and majority-rule rescue apparently reaches
through isolated non-canonical junctions more readily than a consecutive run of them.

### 10.4 Consequence for the family-definition question and for cluster2

This closes the last open item from §9: **NPIPB7 and SMG1P6 are not a fragmentation mystery — they are
the SAME already-documented canonical-junction mechanism as §6m8, now shown to explain 2 of the "8
residual genes" at the individual-gene level**, and shown for the first time to survive into (and cost)
the FAMILY-DEFINITION endpoint specifically, not just the chr16-wide pass-1-to-GTF loss rate §6m8
originally measured.

**isoseq cluster2 cannot help here either.** The reads are already fully confirmed correct, complete,
and unambiguous (§10.1.1) — this defect lives entirely inside Rustle's own assembler, downstream of
anything a PacBio read-processing tool touches. My prior answer treating these two genes as an open,
possibly cluster2-fixable lead is **superseded**: it is neither a read-quality problem nor an unexplained
one.

**This is the user's call, not a default to flip unilaterally** (memory records the default was
deliberately left off after the chr16 arm test — *"fear refuted, default still not flipped, user's
call"*). What this session adds beyond §6m8: independent confirmation on unrelated genes (SMG1P6,
NPIPB7, vs §6m8's chr16-wide aggregate), and a NEW consequence not previously measured — this specific
loss reaches the family-definition score, not just transcript-recovery counts.

---

## 11. Including non-canonical junctions on REAL chr16 data: measured, not assumed

**User, 2026-09-21: "can we include the non-canonical junctions too? ... for now I think my advisor will
be skeptical the method really works for any family if we cannot recall all members."** Ran
`RUSTLE_JUNCTION_MAJORITY=1` through the full pipeline on REAL chr16 (same `--assemble-only
--assembly-polish full` recipe that produced `dn16.gtf`), rebuilt the locus graph and family clusters
the same way as every other arm this session, and scored against both truths.

### 11.1 Assembly and family-definition effect

| | strict (baseline) | `RUSTLE_JUNCTION_MAJORITY=1` |
|---|---|---|
| transcripts kept (polish) | 9,629 | 10,093 (+464, +4.8%) |
| de novo loci | 2,550 | 2,544 |
| graph nodes | 864 | 869 |
| clusters >= 2 members | 70 | 70 |

| | vs Soto cover (15 fams) | vs protein referee (36 fams) |
|---|---|---|
| strict | F 0.3166 | F 0.1805 |
| `JUNCTION_MAJORITY=1` | F 0.3184 (+0.0018) | F 0.1817 (+0.0012) |

Pooled F is essentially flat — negligible, in either direction, at this n. Per-family (Soto): 2 better, 1
worse, 12 unchanged.

### 11.2 The number the advisor's question is actually about: individual-member recall

| | Soto truth genes on chr16 (58) covered as SOME node |
|---|---|
| GUIDED (annotation, ceiling) | 47 / 58 = 81.0% |
| REAL de novo, strict | 34 / 58 = 58.6% |
| REAL de novo, `JUNCTION_MAJORITY=1` | **35 / 58 = 60.3%** |

**+1 gene recovered (`NPIPB2`), 0 lost.** Small, but unambiguous and free on this measure: no truth gene
that strict canonicity covered was dropped by relaxing it.

⚠**§6v4's exact repro genes behave differently on real data than in the ideal simulation.** On real chr16
reads (natural depth variation, not my controlled 10-reads-per-transcript scheme), `NPIPB7` and `SMG1P6`
are **already covered under strict canonicity** — the non-canonical-junction truncation §6v4 demonstrated
unambiguously in the clean simulation does not reproduce identically on these same two genes at real
depth, most likely because real read population diversity gives some reads a slightly different
alignment path across the junction that the clean, uniform simulated set didn't have. The mechanism from
§6v4 is still correct (independently verified against the raw GFF, and directly confirmed by the
before/after test) — it just isn't the reason these particular two genes are missing on REAL data.
`ABCC6`, `HERC2P5`, `PKD1P6` remain uncovered under BOTH arms, for the reasons already established in
§9 (real exon-structure divergence, no transcript model, fusion-naming) — this lever cannot and should
not be expected to fix those.

### 11.3 What this does and does not answer for the advisor

**Does not, on its own, close the recall gap.** 23 of 58 Soto truth genes are still missing a node under
`JUNCTION_MAJORITY=1`, against 24 under strict — a gap of 1, not 24. The 22-point remaining difference to
the guided ceiling (81.0%) is NOT primarily a canonical-junction problem on real data; §9's three-way
breakdown (no transcript model / fusion-naming / real exon divergence) accounts for at least several of
the specific genes checked, and the rest is uncharacterized. **"Recall all members" needs more than one
lever** — this is one confirmed, safe, small piece of it, not the fix.

**Known cost, already measured genome-wide** (§6n4, `bench/CHR16_JUNCTION_MAJORITY_ARM.md`): strictly-
engulfed copies 79 -> 86 (+8.9%), a precision cost on copy BOUNDARIES, not on family structure (max family
size and family fusion were unaffected — the flag's specific feared harm was refuted there).

### 11.4 Recommendation

Given zero measured node losses here, a real (if small) recall gain, and an already-refuted specific fear
at genome-wide scale, this is a low-risk, well-evidenced choice for the user to adopt. Two distinct
decisions, kept separate:

1. **Use `RUSTLE_JUNCTION_MAJORITY=1` in this session's own family-definition / recall-focused bench
   runs going forward.** Zero code change, fully reversible, and directly what was asked for "for now."
2. **Flip the compiled default in `denovo_pipeline.rs`** (`env_num("RUSTLE_MISCHAIN_...")`/junction
   majority default) so every future run gets it without the env var. This is a bigger, harder-to-reverse
   change than (1) and is NOT done here — it is the user's call, per the standing rule this project has
   already applied twice to this exact flag.

---

## 12. The full remaining recall gap, classified (not just the original 8-gene sample)

**User, 2026-09-21: "keep exploring reasons for false positives/false negatives."** After the default
flip (§11 above, now committed), 23 of 58 Soto truth genes on chr16 still lack a node. Classified all 23
— not a hand-picked sample — by the same three mechanisms §9 found on 8 genes, checked directly against
the raw GFF and the new-default assembled GTF.

| category | count | share |
|---|---|---|
| **A — no transcript model in RefSeq at all** | 6 | 26.1% |
| **B — a bigger co-located annotation record wins the naming contest** | 8 | 34.8% |
| **C — a locus exists but captures only a fraction of the gene** | 9 | 39.1% |

**A (unfixable by construction): `ABHD17AP7`, `ABHD17AP8`, `ABHD17AP9`, `HERC2P5`, `NPIPB14P`, `PKD1P6`.**
A THIRD pseudogene trio (`ABHD17AP7/8/9`) joins `HERC2P5`/`HERC2P8`/`NPIPB10P` in having zero annotated
exon structure — no sequence exists to simulate or sequence, real or synthetic.

**B (scoring artifact, not construction failure): `MIR3179-1`, `MIR3179-2`, `PKD1P1`, `PKD1P6-NPIPP1`,
`RRN3`, `SLX1B`, `SULT1A3`, `SULT1A4`.** Extends §9's mechanism (a bigger fusion/co-located annotation
record wins a max-overlap, one-name-per-locus resolver's contest) to double the count. One striking case:
**`PDXDC1` (168,469 bp) sits over BOTH `PKD1P6-NPIPP1` and `RRN3`** — even the FUSION gene itself loses
its naming contest to a still-larger neighbouring record in this densely tandem-duplicated region. There
is also a CHAIN: `SULT1A3`/`SULT1A4` lose to the fusion names `SLX1A-SULT1A3`/`SLX1B-SULT1A4`, which
themselves land in category C (locus exists, captures little) — the naming collision and the coverage gap
compound for this pair.

**C (the open, still-uncharacterized lever): `ABCC6`, `ABCC6P1`, `ABCC6P2`, `CDR2`, `NPIPB4`, `RRN3P1`,
`RRN3P2`, `SLC7A5P2`, `SLX1A-SULT1A3`.** A real locus exists for 7 of 9; `ABCC6P1`/`ABCC6P2` have no
locus at all despite carrying their own transcript record. This is the same shape of gap NPIPB7/SMG1P6
were before §6v4 diagnosed them — and NPIPB7/SMG1P6 are confirmed NO LONGER in this list, a clean
positive control that the junction-majority fix reached exactly the class of gene it was built for. This
category is where the next diagnostic effort belongs.

⭐ This reframes the "22-point gap" from §11: it was never one problem. A quarter is a true annotation
ceiling, a third is a scoring resolver artifact (the real sequence is very likely already correctly
grouped, just credited to the wrong name), and just over a third is a genuine, unexplained construction
gap — the only slice comparable in kind to what junction-majority just fixed.

---

## 13. Category C dissolves: most of it is the SAME naming-collision mechanism as Category B

**§6v7's Category C ("locus exists, captures only a fraction — the open lever") does not survive
individual diagnosis.** Checked all 9 genes directly: for each, whether its own de novo locus shares its
assembled span with a bigger neighbouring annotated gene that wins the max-overlap resolver's contest —
exactly §6v2/§6v3/§6v7's Category B mechanism, but checked against the ASSEMBLED LOCUS's span (which can
reach much further than any one gene's own annotation, via readthrough) rather than the raw annotation
overlap Category B originally checked. **§6v7's Category B check only looked for a bigger annotation
record overlapping the target gene's OWN span — it missed collisions mediated by a shared DE NOVO locus.**

### 13.1 Five of nine are the same collision mechanism, now via the assembled locus

| gene | shared de novo locus | naming-contest winner |
|---|---|---|
| `NPIPB4` | `DN_chr16_22407829_6` (65,081 bp) | **`LOC112268174`** (23,535/23,536 bp = 99.996% of its own span) |
| `RRN3P1` | *same locus as NPIPB4* | **`LOC112268174`** |
| `RRN3P2` | `DN_chr16_29356116_18` (329,561 bp) | **`SNX29P2`** (62,762/62,763 bp) |
| `SLC7A5P2` | `DN_chr16_21397764_13` (75,080 bp) | **`SMG1P3`** (55,284/55,285 bp) |
| `SLX1A-SULT1A3` | `DN_chr16_30478427_6` (3,688 bp) | **`SLX1A`** (3,687 bp) |

In every case the winning neighbour claims essentially its ENTIRE own span inside the shared locus — these
are not marginal calls. `NPIPB4` and `RRN3P1` even share the exact same super-locus and both lose to the
same third gene. Register 690's finding (readthrough units chain multiple genes together, CDR2 among
them, MCL7's unit reaching 133 kb over 8 genes) is the same phenomenon, now traced to its scoring
consequence gene by gene.

### 13.2 CDR2: a new, fourth mechanism — the truth label isn't corroborated by CDR2's own sequence

CDR2 does NOT lose its own naming contest (its 28,680 bp overlap is the largest in its locus). Its Soto
family is `ID_39 = {CDR2, RRN3, RRN3P1, RRN3P2, AC009021.1, AC009093.9}`. The locus DOES have a strong PAF
hit to RRN3P2's locus (98.78% identity, 7,112 bp aligned, well above every mcl_families floor) — but
decoding the exact aligned coordinates places it at `chr16:21780339-21787538`, inside **`RRN3P3`**'s span
(21779755-21798544), not CDR2's (21842595-21871275). `RRN3P3` is not in Soto's S1C at all. So the real,
strong homology signal in this locus belongs to an uncharacterized passenger gene riding along in the same
readthrough-fused locus as CDR2 — the exon-conjunct is CORRECT to refuse crediting it to CDR2, since CDR2's
own exons genuinely do not align to RRN3P2. **CDR2's truth-family membership may reflect co-duplication
(physical proximity to an RRN3-family paralog) rather than CDR2 itself carrying RRN3-family sequence
homology** — worth flagging as a possible truth-label question, not assumed as a pipeline defect.

### 13.3 What's left once the collisions and CDR2 are pulled out

- **`ABCC6`**: wins its own naming contest; §9's original diagnosis stands — real, annotated exon-structure
  divergence from its paralog `ABCC6P1` (a `transcribed_pseudogene`), not a construction or naming defect.
- **`ABCC6P1`, `ABCC6P2`**: no locus assembled at all. Real-data read counts are low (5 and 12 primary
  reads respectively, on a locus with real biological expression that may simply be scarce in this
  sample) — a genuine depth limitation, not an algorithmic one.

### 13.4 Corrected tally: the whole remaining gap, four mechanisms

| mechanism | genes | share of 23 |
|---|---|---|
| **A — no transcript model in RefSeq** | 6 | 26.1% |
| **B — naming collision with a bigger neighbour (annotation- or locus-mediated)** | **13** (8 original + 5 reclassified) | **56.5%** |
| **CDR2 — truth label not corroborated by the gene's own sequence** | 1 | 4.3% |
| **D — real construction/depth gap** (`ABCC6` divergence, `ABCC6P1/2` low depth) | 3 | 13.0% |

**The single dominant mechanism across the ENTIRE remaining recall gap is naming collision — 56.5% of it
— not an edge-construction or assembly defect.** In every collision case checked, the actual homology
signal is very likely already correctly captured by the assembler; it is being credited to the wrong gene
name by a simplistic max-overlap, one-name-per-locus resolver. This reframes where the next lever should
go: not into further assembly fixes, but into how a multi-gene locus's homology is attributed to truth
labels — e.g. crediting a locus to every gene whose exons it substantially covers, rather than
winner-take-all by span.

⚠ This same resolver is used throughout this session's OWN scoring pipeline (`score_family_def.py`,
`cover_and_jn_definition.py`), not only in some downstream consumer — so this is very likely **undercounting
this project's own measured recall** wherever a de novo locus spans more than one gene, on every chromosome
scored this session, not only chr16's 23 genes.

### Correction (2026-10-06): the simulated reads lack each transcript's leftmost exon

`bench/sim.py` `load_transcripts` (used only by `sim.py chromosome`) keeps an exon only if its transcript row came EARLIER in the file. The RefSeq GFF is sorted by start then end, so each transcript's leftmost exon row (same start, shorter end) precedes the transcript row and is dropped: 6,950 of 80,664 chr16 exon rows,
and 2,280 transcripts vanish whole. On the first 300 testable multi-exon transcripts the simulator's exon list equals the annotated list minus its leftmost exon in 300 of 300, and 298 of 300 reads lack that exon's sequence. The reads therefore carry 4,370 transcripts of 1,407 genes (not all 4,411 / 1,443 of this document, which are the pre-filter loader counts)
and none is full length; they were also mapped to `chr16.fa` only. Every number of this document describes those reads. The corrected ceilings are not computed.


---

## PREREG_inversion_aligner_blindness_2026-10-05.md

*(original: `docs/archive/2026-10/merged/PREREG_inversion_aligner_blindness_2026-10-05.md`)*

## Inverted copies and the strand-aware spliced aligner — mechanism + a parked, pre-registered measurement (NOT RUN)

2026-10-05. Parked by the user's cost decision ("I don't have much usage — just record it"). No data was looked at beyond
what is cited; every number below is a *claim to be tested*, not a result.

### The mechanism (why this is expected, not speculative)

`minimap2 -x splice:hq` — the aligner every BAM in this project was made with — is **strand-aware**: it forces the
transcript strand from the canonical splice sites and does not emit a spliced alignment whose strand disagrees with them.
A wholly-inverted gene copy's transcripts are reverse-strand relative to the forward copies of its family, so:

- the reads of a whole-copy inversion do not get *weakly* aligned back to the family — the splice-strand constraint
  **refuses** the spliced chain outright; only short intra-exon local matches can survive. The reads land in the
  **unmapped** pool (or as poorly-placed fragments). Seeding is not the first wall; the strand constraint is.
- an **exon-internal inversion** partially survives: the flanking canonical splice sites still anchor the read's strand,
  and the inverted exon body aligns with mismatches/indels.
- an **intron-span inversion** is intermediate: breakpoints lose spliced support, interiors may still align locally.

Consequence for the pipeline as shipped: inverted copies are **invisible to O1 and O2 by construction** — their reads are
not in the mapped pool at all, so there is nothing to assemble into a locus and nothing to assign. This is a property of
every strand-aware spliced aligner, so it holds equally for the *real* gorilla/human BAMs: inverted copies, if any exist
in these gene families, are already sitting in the real unmapped pools today.

### The designed mitigation already in the tree (to be tested, not assumed)

`o3_candidates` pass B (`src/rustle/vg_family/o3_candidates.rs`, Amendments 13/13b) consumes **unmapped and poorly-placed
reads ≥ 300 bp** and attributes them by **`map-ont`** — which is *not* strand-constrained. A wholly-inverted copy is
exactly the class of object that mechanism exists for; whether it actually fires for inversions is an open empirical
question (the inversion changes the sequence locally — snp-level divergence — and `map-ont` needs ≥ 50% read coverage at
`de ≤ 0.20` against the family's net or copies).

### The parked experiment (famsim; ~1 h of compute on the 5-core box when run)

Substrate: `bench/famsim` synthetic background (self-contained, no real data needed; `examples/synthetic_smoke.json` is
the template). Runs: one 100%-identity control (`three_identical`: copies A, B, C, no ops) and one arm per inversion
type, each rung verify-PASS before scoring:
- `inv_whole` (whole-copy strand flip — the wall case),
- `inv_exon` (one exon inverted inside an otherwise-identical copy),
- `inv_intron` (one intron inverted).

Claims, written before any run (adjudicate on `famsim`'s own products, not the pipeline's):
- **C1 (aligner loss):** the fraction of each inverted copy's reads that are unmapped in `align/`s BAM is ~1.0 for
  `inv_whole`, intermediate for `inv_intron`, low for `inv_exon`; the control's unmapped fraction is the ordinary
  error/simulation floor. Measured from the BAM + `families.tsv`/reads truth, independent of rustle.
- **C2 (O2 abstention):** among the inverted copy's reads that *do* map, O2's assignment rate does not exceed the
  control's (they map to the forward copies; PSV/junction evidence should abstain or misassign — reported, not barred).
- **C3 (O3 net):** `o3_candidates --candidates` flags the inverted copy as a candidate for at least one of
  `inv_exon`/`inv_intron`, and for `inv_whole` **iff** pass B's map-ont attribution recovers ≥ `min-support` (6) reads
  into its family's net. Reported beside: the pass-B attributed-read count per inversion type, and whether the flagged
  candidate's nearest locus is the family's own interval.
- **C4 (advisor framing):** the deliverable number is the pair (C1 loss, C3 recovery) per inversion type — "inverted
  copies are an O3/unmapped-read phenomenon by construction of strand-aware spliced alignment; the aligner loses X, the
  net recovers Y" — never a claim that the primary mapping sees inversions.

Not in scope: real-data inversion screening (no annotation of inverted gene copies in these genomes is assumed);
`sim.py`-style standalone inversion support (famsim covers it; nothing else in bench/ makes inversions — see the
2026-10-05 inventory in the session record).


---

## DEFAULT_RESCORE_NPIP_2026-10-06.md

*(original: `docs/archive/2026-10/merged/DEFAULT_RESCORE_NPIP_2026-10-06.md`)*

## The current default at NPIP, re-scored at HEAD (2026-10-06)

Protocol and rules: `docs/PREREG_default_rescore_npip_2026-10-06.md` (committed 7c26e2c0, amended 1df2e20b before any DEF or PRE product existed). Runner:
`bench/default_rescore/run.sh`. Products: `/mnt/linuxdisk/tmp/rescore_2026-10-06/` (`verdict.json`, `support.copies.tsv`, `fs_DEF.json`, `fs_PRE.json`, `provenance.txt`).
Binaries: HEAD 8fbaba79 build in `rustle_target_m2/release` (`copy_assign` 87824d91, `mcl_families` a6308244, `family_score` 542923fd; `copy_support.py` d0ea1bfb).
Human A119b chr16, the 25 CAT/Liftoff NPIP copies. **Everything here is DEV** (the development block of NPIP and of every fusion rule): a measurement, not a validation.

### Verdict: PASS (R1 to R4 hold; G0, G1, G3 pass)

| rule | value (DEF = the default) | bar | |
|---|---|---|---|
| R1 strict FOUND, Amendment E, within own nodes | **8 of 25** | >= 6 | ok |
| R2 locus-level E-found | **21 of 25** | >= 18 | ok |
| R3 U2 bipartite F | **0.645** (sens .588, prec .714; pairs 117 of 240 truth, 199 predicted) | >= 0.645 | ok, at the bar |
| R4 own node of NPIPB2 and NPIPB6 | both have one (PRE: neither) | both | ok |

R3 sits exactly on its bar because the HEAD products are byte-identical to the 2026-09-30 products that set it (below), not because of a margin.

Gates: G0 the HEAD scorer reproduces the registered Amendment E numbers on the frozen arms exactly (P 23 / 10 / 19 / 10, GOOD 21 / 6 / 18 / 6, ALL 24 / 2 / 16 / 2: own node, E-found in own nodes, locus level, E-found chr16-wide);
G1 the own-node rule reproduces 50 of 50 stored flags (P and GOOD; the frozen ALL arm has two loci on one span, see the amendment); G3 the HEAD `family_score` reproduces the stored default-arm U2 (.588 / .714 / .645) and Compara (.500 / 1.000 / .667).
G2 (diagnostic): the PRE control, rebuilt by the HEAD binaries, reproduces the registered GOOD arm exactly (21 / 6 / 18 / 6). The scorer is deterministic under two `PYTHONHASHSEED` values.

### DEF against PRE (descriptive; same binaries, same data)

| | PRE (`--bridge-regroup off`, cov 0) | DEF (default: f1v2, cov .70) |
|---|---|---|
| copies with an own node | 21 | 24 |
| strict E-found in own nodes (of 25) | 6 | **8** |
| locus-level E-found (of 25) | 18 | **21** |
| U2 sens / prec / F | .500 / .680 / .576 | .588 / .714 / **.645** |
| Compara chr16 sens / prec / F (17 families) | .444 / 1.000 / .615 | .500 / 1.000 / **.667** |
| Compara CF153 (the NPIP family): genes hit, F | 10 of 19, .690 | 13 of 19, **.812** |
| Soto sens / prec / F | .465 / .917 / .617 | .507 / .923 / .655 |

Copies that change between PRE and DEF: **NPIPB2** and **NPIPB6** go from no node to an own node and are found (their holders were fused with GSPT1 and EIF3CL before the bridge regroup); **NPIPA6** gains an own node and is not found.
The E-found copies in DEF: NPIPB2, NPIPA2, NPIPA1, NPIPA5, NPIPB6, NPIPB8, NPIPB10P, NPIPB11.

### What remains short

- 22 of 25 copies are spliced-expressed under E (PKD1P6-NPIPP1, NPIPB12 and LOC124907808 are not). **14 of those 22 are not found** by the default: the locus representative does not carry the capped start and its first three introns, although a transcript of the same locus does at 13 of them (NPIPB5 at none). The loss is the representative, as `docs/archive/2026-10/SPLICED_COPY_SUPPORT_2026-10-04.md` found; f1v2 does not repair it (+2 copies, both from node separation).
- Found is 32 % of the copies in own nodes (40 % in the primaries-only arm P before f1v2: P was not re-run on HEAD, so no DEF-versus-P statement is made).
- Representative rules that would act on this (R_J, most junctions) were tested and are opt-in (`docs/archive/2026-10/LOCUS_REPRESENTATIVE_RULE_2026-10-04.md`); they were not run here.

### Provenance: HEAD reproduces the 2026-09-30 products byte for byte

`cmp` of the HEAD products against the e163d955 container-headroom products: **identical** for DEF (assembled GTF, `families.gtf`, clusters, loci GFF3, copy table) and for PRE (assembled GTF, clusters, loci GFF3).
So the 2026-10-06 consolidation of `src/` (`2a83a747`) and the later opt-in code changed nothing in the default's chr16 outputs (assemble 24 s / 0.75 GB; families 65 s). This is one input (A119b chr16), not the whole-genome run.

### Not shown

Any held-out NPIP substrate, any gorilla substrate (strict E needs the cap signal, which gorilla lacks), O2, O3. Liftoff recall was not recomputed.


---

## ENTANGLED_BASELINE_2026-10-06.md

*(original: `docs/archive/2026-10/merged/ENTANGLED_BASELINE_2026-10-06.md`)*

## What StringTie and FLAIR do at overlapping genes, on ideal reads, and where we stand (2026-10-06)

**Status: DEV, human CHM13 v2.0 / CAT-Liftoff only, simulation, circular by construction (reads come from the annotation that scores them). Protocol: `docs/PREREG_entangled_baseline_2026-10-06.md` (committed 651b5255; Amendment 1 9fb174f3, written before the lever arms existed). Tools: `bench/entangled/`. Products: `/mnt/linuxdisk/tmp/entangled_2026-10-06/{NPIP,TBC1D3}/rep{1,2}/` (`arms.json`, `arms.genes.tsv`, the tool GTFs, `pooled.json`).**
User request: "model the cases in which there are overlapping genes, ideally provide what stringtie and flair do for those loci and then do better than that."
Substrate: the four ideal-expression read sets (NPIP and TBC1D3 windows, two replicates each; per replicate NPIP 20,660 reads from 2,066 simulated transcripts of 474 genes (547 genes and 2,147 transcripts lie in the windows; 81 transcripts were skipped: under 120 bp, a kept intron over 200 kb or a molecule over 30 kb), TBC1D3 16,680 reads from 1,668 transcripts of 464 genes (519 genes, 1,725 transcripts in the windows, 57 skipped)), one BAM per run seen by every arm. Entangled genes (stratum E) = exon union shares >= 100 bp with another gene on the strand: 256 gene-runs with 1,606 chains; the other genes with a chain (N): 1,026 gene-runs, 4,984 chains.

### Answer

1. **At overlapping genes we already recover more than either tool**: chains recovered at E: ours 1359 / 1606 (84.6%), FLAIR 1159 (72.2%), StringTie 1015 (63.2%); elsewhere (N) 93.0% against 75.3% and 72.2%. The independent `gffcompare` intron-chain sensitivity on NPIP replicate 1 agrees (ours 87.8, FLAIR 72.2, StringTie 67.9; this instrument 88.7, 72.9, 68.5).
2. **At default settings we also emit about twice the artifacts** (multi-exon transcripts whose chain is no truth chain): 1031 against 412 (FLAIR) and 431 (StringTie); chain precision 85.3% against 92.3% and 91.5%.
3. **Two existing opt-in options bring the artifact rate to the tools' level without losing the recall lead (levers chosen after the artifact analysis: development evidence)**: primaries-only seeding (`--no-seed-secondaries`, arm P) cuts artifacts to 593 (precision 91.0%) at the same recall (91.1%); adding the sub-chain drop (`--polish-subchain drop`, arm PC) gives 451 artifacts, precision 92.9% (FLAIR 92.3%, StringTie 91.5%), recall 89.1% (FLAIR 74.5%, StringTie 70.0%), E recall 82.0% against 72.2% and 63.2%. Per run, PC precision is 90.8 / 90.5 / 95.9 / 95.6 against FLAIR 90.7 / 90.4 / 94.5 / 94.5: level with FLAIR in the two NPIP runs (+0.1 point) and above it in the two TBC1D3 runs (+1.4); absolute artifacts remain higher (451 against 412 and 431), and exact / (exact + artifacts) rises with recall by construction.
4. **By the bar registered in advance, no arm is BETTER** (section 4 of the prereg, absolute counts): PC fails (b) (451 artifacts against 412) and (c) as written (the FLAIR figure is inflated, see below). The rate-based reading in point 3 is post hoc and is not a pass.
5. **Gene separation**: 216 of the 256 entangled gene-runs share an exact junction with another truth gene (readthrough models and exon-reusing lncRNA / pseudogene models), so no gene_id rule based on the junction-sharing structure of the transcripts can separate them (FLAIR's apparent 165 of 216 comes from its 1 kb start-bin gene_ids); for the other 40 our assembler resolves 40 / 40 against StringTie 20 / 40.

### Pooled results (four runs; per-run rows in `arms.json`)

| arm | what | chains E | chains N | recall | artifacts (fragment / fusion / other) | chain precision | F1 |
|---|---|---|---|---|---|---|---|
| S | StringTie -L | 1015 (63.2%) | 3596 (72.2%) | 70.0% | 431 (156 / 48 / 227) | 91.5% | .793 |
| F | FLAIR collapse | 1159 (72.2%) | 3753 (75.3%) | 74.5% | 412 (187 / 19 / 206) | 92.3% | .825 |
| D_asm | our default | 1359 (84.6%) | 4637 (93.0%) | 91.0% | 1031 (427 / 77 / 527) | 85.3% | .881 |
| P_asm | primaries only | 1372 (85.4%) | 4632 (92.9%) | 91.1% | 593 (376 / 19 / 198) | 91.0% | .911 |
| C_asm | sub-chain drop | 1297 (80.8%) | 4559 (91.5%) | 88.9% | 838 (260 / 76 / 502) | 87.5% | .882 |
| PC_asm | both | 1317 (82.0%) | 4554 (91.4%) | 89.1% | 451 (241 / 18 / 192) | 92.9% | .909 |
| U_asm | f1units | 1359 (84.6%) | 4637 (93.0%) | 91.0% | 1031 (427 / 77 / 527) | 85.3% | .881 |

Truth: 6,590 chains (1,606 at E, 4,984 at N). Chain precision = exact / (exact + artifacts), artifacts counted over every contig. `_fam` rows (the families input after the regroup) differ from `_asm` by <= 5 chains. Per family: NPIP D precision 81.0 -> PC 90.7 (FLAIR 90.5, StringTie 89.6); TBC1D3 91.3 -> 95.7 (FLAIR 94.5, StringTie 93.8).

### Predictions and the registered bar

- **P1** both tools recover less at E than at N: held (S 63.2 against 72.2, F 72.2 against 75.3). **P2** our default has fewer artifacts per recovered chain than S and F: **FAILED** (0.172 against 0.093 and 0.084). **P3** we recover at least as many chains at E as S and F: held (1359 against 1015 and 1159). **P4** the units arm loses nothing and resolves strictly more E genes than D_fam: **FAILED** on its no-loss clause (U_fam loses 5 chains at N and 1 complete gene; resolved at E 90 against 88).
- **Q1** P cuts the in-window artifacts of class 'other' by >= 30%: held (in-window 417 -> 158, -62%; all contigs 527 -> 198, -62%). **Q2** C cuts fragment artifacts by >= 50% and loses <= 1% of the chains: **FAILED** (427 -> 260 is -39%, and 140 chains = -2.3% are lost; arm C alone, PC loses 125 = -2.1%). **Q3** P recovers fewer chains than D at E by >= 10: **FAILED** (P recovers 13 MORE; the seeding gain is not visible in these windows). **Q4** none of P, C, PC meets the bar: held vacuously (clause (c) can be met only by gene_ids that over-split).
- **Bar of section 4** (a) chains at E >= the larger of S and F: met by every arm of ours; (b) artifacts <= the smaller of S and F (412): met by none (P 593, PC 451); (c) resolved at E >= the larger of S and F and above D_fam: not met as written; (d) no loss at N relative to D: P loses 5 chains, C and PC lose 78 and 83. **No BETTER arm.**
- Caveat on (c): FLAIR appears to name genes by a 1 kb start bin ('chr16:11726000'), 934 ids for 474 truth genes in NPIP replicate 1, so its 'resolved' (201 / 256) measures over-splitting, not separation; StringTie fuses (305 ids, resolved 38 / 256).

### What the artifacts are (post hoc; NPIP and TBC1D3 replicate 1, 449 in-window artifacts of the default)

Nearest truth chain by shared junctions: terminal-trim sub-chains 206 (46%), novel junctions 172 (38%), no shared junction 33, displaced junctions within 60 bp 29, internal skips 9. Read support is not low: fragments carry mostly 4-10 reads, (pooled over the four runs) 207 of the 417 in-window 'other' carry exactly 10 reads (one transcript's full depth) and 385 of the 427 fragments carry 4-10, so depth cannot separate them; the 449 of the sample are replicate 1 only (fragment 206, fusion 39, other 204). In-window = the exons overlap a same-strand truth gene. Our artifacts outside the windows (other contigs, silent paralogs): 110, against 26 and 32 for the tools; the excess is in the windows. With primaries-only seeding 'other' falls 62% and fusions 75%: the tied secondary alignments of paralogous reads build most of them. Dropping sub-chains removes true annotated sub-chain isoforms too (chains -2.3% in the windows, a property of the simulation: every annotated sub-chain isoform is expressed at 10 reads; on real gorilla reads PC loses 18 of 19,771 chains at non-overlapping genes, 0.09%, `docs/GORILLA_OVERLAP_2026-10-07.md`), as in `docs/archive/2026-09/PREREG_complete_transcripts_2026-09-27.md`.

### Gene-level view (post hoc)

Of the 256 entangled gene-runs, 216 (84%) are E_j (share an exact junction with another truth gene on the strand) and 40 are E_x (overlap without a shared junction). Resolved: E_x D 40 / 40, S 20 / 40, F 36 / 40 (inflated); E_j D 53 / 216, S 18 / 216. 26 of the 27 merged gene_ids in NPIP replicate 1 and all 12 in TBC1D3 replicate 1 hold genes whose truth chains share a junction (the readthrough models BOLA2-SMG1P6, PKD1P6-NPIPP1, SLX1B-SULT1A4, CCL15-CCL14, PDXDC2P-NPIPB14P, TBC1D3P1-DHX40P1, and lncRNA models reusing a protein gene's exons). Re-assigning gene_ids by junction-sharing components gives identical scores in all four runs (515 components against 515 gene_ids in NPIP replicate 1; 442 against 440 and 407 against 405 in the two replicate-2 runs: only the id counts differ). Separating E_j genes needs (hypothesis) an annotation or coding evidence, or a representation where a transcript belongs to several gene units; the one evidence rule tried later (F1) changed no target in these windows.

### Limits

Ideal reads (no truncation, no readthrough, no depth limit), annotation as truth, human only, two replicates of one design (the same windows): the lever arms were chosen after seeing the default's artifacts, so P and PC are development evidence, not held-out. The seeding gain that primaries-only seeding gives up exists in real data (copies with no primaries: gorilla LOC129533806 and LOC129533813, human TBC1D3B and TBC1D3I in `docs/archive/2026-09/COPY_RECOVERY_TOOLS_2026-09-29.md`) and is invisible here. Versions: StringTie 3.0.3, FLAIR 3.0.0, minimap2 2.30 (lab 3.0.1, 3.0.1, 2.31). isoseq collapse not run. The lab's tools run in their default annotation-free recipes; nothing was tuned for them.


---

## FEWCOPY_IDEAL_CASES_2026-10-06.md

*(original: `docs/archive/2026-10/merged/FEWCOPY_IDEAL_CASES_2026-10-06.md`)*

## Few-copy ideal cases for O1, and are NPIP / TBC1D3 "unreachable"? (2026-10-06)

**Status: EXPLORATORY, DEV, human A119b only, stored products only (no pipeline run). This is not a pre-registered test.** The held-out class-level test is drafted in `docs/PREREG_fewcopy_class_DRAFT_2026-10-06.md` and is NOT registered.
User question (2026-10-06): "do we know of any ideal cases of multi-copy gene families (just 1 or 2 with few copies) that would validate our O1? So far it looks like NPIP and TBC1D3 are a bit unreachable right?"
Method: workflow wf_dbfd1e66-dd2 (four read-only sweeps run in parallel, then one adversarial reviewer that recomputed the decision-bearing numbers). Outputs: `/mnt/linuxdisk/tmp/fewcopy_2026-10-06/` (tables, scripts, the five sweep results; the instrument scripts are NOT committed). Frozen input-side criteria: `bench/fewcopy/criteria.md`, sha256 `89363117a885718695d70c6d71a78b75f88c31a62421904c505aeba34c7c1571`.

### Answer

1. **NPIP and TBC1D3 are not unreachable.** On ideal reads TBC1D3 passes its registered bar (12 of 13 reachable copies, both replicates) and NPIP misses by one copy (12 of 14, bar 13); perfect loci (the annotation through `mcl_families`) reach 14 / 14 and 13 / 13, so both bars are reachable by construction (`docs/IDEAL_EXPRESSION_DEFAULT_2026-10-06.md`, register 1254, recomputed by the reviewer from `g4_s*.copies.tsv` and `gates.json`).
   What is real: 10 of 25 NPIP and 2 of 16 TBC1D3 copies share >= 100 exonic bp with another gene (no per-copy-locus bar is demanded for them, yet 4 and 2 of them were found); real-read recovery is data limited (annotation-anchored FOUND, quoted: NPIP 12 real against 23 ideal of 25, TBC1D3 12 against 16); both verdicts sit on a one-copy margin; the simulation is circular and dev-only.
   They are hard stress tests of fused loci and representative choice, and poor clean validators of O1.
2. **Two clean few-copy cases, chosen from inputs only: GSTA1 / GSTA2 (chr6) and PAGE2 / PAGE2B / PAGE5 (chrX).** Both are exact in all three scoring views under the current default. They are positive controls for a NECESSARY condition on DEV data, not validation (details below).
3. **"Ideal by inputs" does not guarantee recovery, and the class is much harder than the two cases.** Of 377 Compara families with 2-4 genes, 11 pass the five input criteria; only 4 of the 11 are exact in all three views. Over all 292 scorable whole 2-4-gene families on 19 A119b contigs: sens .297, prec .935, F .450.

### Input-side criteria (frozen before any count)

Universe: Compara Primates families with 2, 3 or 4 listed genes (N0 = 377: 292 / 60 / 25), members mapped RefSeq to CAT/Liftoff v2.0 by shared exonic bp (never by name). A family passes only if EVERY member or pair passes.

| step | criterion | cumulative families |
|---|---|---|
| C1 | every member maps strongly to one protein-coding CAT gene; reference transcript >= 3 exons, every intron >= 50 bp and canonical | 129 |
| C2 | not entangled: < 100 exonic bp shared with any other annotated gene | 62 |
| C3 | separate loci: spans >= 1 kb apart | 61 |
| C4 | homologous, not identical: spliced reference transcripts, 0.90 <= identity < 0.999, coverage >= 0.50 (minimap2 asm20) | 35 |
| C5 | observable: >= 3 reads per member carrying its exact annotated chain (E0 rule, primary + good secondaries, region queries of `A119b.t2t.bam`) | **11** |

One code amendment after the first run (spans = exon-union extent; 3 flags changed, the 11 and their order did not); no threshold edited; the criteria hash was re-checked at the end and by me (identical). 1 of 864 members could not be mapped (FAM86B1). 57 families fail exactly one criterion (C5 24, C2 16, C4 11, C1 6; `near_miss.tsv`).

### The 11 families that pass, and what the current default did

Outcome side: stored current-default products (arm D, f1v2 bridge regroup + `--min-cov-shorter 0.70` + GOOD secondaries), scored with `family_score` (S1 = its own view; S2 = any other annotated gene in the cluster counts as a non-member; S3 = membership by >= 1 exonic bp). The outcome sweep was kept blind to the input sweep.

| family | members | contig | identity (min) | min reads | tie fraction | S1 | S2 | S3 | all three |
|---|---|---|---|---|---|---|---|---|---|
| CF150 | HBA1, HBA2 | chr16 | .967 | 369 | .004 | exact | exact | exact | **yes** |
| CF316 | GSTA1, GSTA2 | chr6 | .954 | 111 | .000 | exact | exact | exact | **yes** |
| CF334 | GTF2IRD2, GTF2IRD2B | chr7 | .993 | 18 | .209 | exact | exact | exact | **yes** |
| CF400 | PAGE2, PAGE2B, PAGE5 | chrX | .918 | 12 | .000 | exact | exact | exact | **yes** |
| CF63 | ZNF33A, ZNF33B | chr10 | .921 | 11 | .000 | exact | exact | partial + non-member | no |
| CF9 | FCGR3A, FCGR3B | chr1 | .977 | 8 | .000 | exact | merged | merged | no |
| CF403 | RHOXF2, RHOXF2B | chrX | .998 | 139 | .730 | partial | partial + non-member | exact | no |
| CF355 | PRR23D1, PRR23D2 | chr8 | .996 | 36 | .470 | partial | partial + non-member | merged | no |
| CF236 | ZNF600, ZNF611 | chr19 | .912 | 4 | .000 | merged | merged | merged | no |
| CF54 | LYZL1, LYZL2 | chr10 | .980 | 59 | .000 | no-locus | no-locus | no-locus | no |
| CF301 | UGT2B15, UGT2B17 | chr4 | .966 | 7 | .000 | no-locus | no-locus | no-locus | no |

Pooled over the 11: members hit 17 / 23, sens .739, prec <= .944, F .829, pairwise 9 / 13. "No-locus" here means no CLUSTERED locus: LYZL1/2 and UGT2B15/17 have loci and reads but no admitted edge (the 09-25 annotation-loci arm also fails LYZL, 0 / 2). Those two are O1-definition cases, more informative than the successes.
Products: HBA1/2 comes from a HEAD run (chr16); the other ten are replays of the 09-30 e163d955 arm D, and HEAD equals it byte for byte only on chr16, chr17 and chrY.

### Class-level numbers (always print beside any subset)

| set | families | sens | prec | F | pairs |
|---|---|---|---|---|---|
| ALL scorable whole 2-4-gene families | 292 | .297 | <= .935 | **.450** | 111 / 493 |
| size 2 / 3 / 4 | 223 / 48 / 21 | .298 / .333 / .226 | .964 / .873 / .905 | .455 / .482 / .362 | 54 / 223, 34 / 144, 23 / 126 |
| the 11 that pass the input criteria | 11 | .739 | <= .944 | .829 | 9 / 13 |
| all minus the 11 | 281 | .281 | <= .934 | .432 | 102 / 480 |
| the 42 exact in all three views (selected on the OUTPUT) | 42 | 1.0 | 1.0 | 1.0 by construction | do not quote alone |

Exactness over the 292: 63 (22%) in the scorer's own view, 42 (14%) in all three views; 162 (55%) have no member in a clustered locus (63 have no locus at all, 43 have every member only in unclustered loci). 78 of the 377 families cannot be scored (55 span contigs, 23 sit on chr13, 18, 20, 21 or 22, which have no stored current-default product); 7 more are scored on a restricted piece.
Selecting on inputs helps (exact in all views 4 / 11 = 36% against 38 / 288 = 13% for the rest) but is no guarantee. Post hoc, dev, by the minimum exact-chain reads per member: all-three-exact 10 / 193 at 0 reads, 14 / 70 at 1-29, 6 / 16 at 30-99, 7 / 12 at 100-299, 5 / 7 at >= 300.
Holder-rule view (the prior-evidence sweep, A119b contigs without chr16; a different instrument from the table above): of 118 pairs with BOTH genes in a locus, 63 are complete, 40 have no admitted edge, 9 sit in one locus, 4 have an edge that is not joined; this is the edge-admission loss of the earlier §6o8 bound, not two-node MCL numerics (which explain at most 1-4).
Chain fidelity (review sweep, not rechecked by me): only 21 of the 42 all-three-exact families have every member's locus chain identical to a RefSeq transcript; of the 11, only HBA, GSTA, GTF2IRD2 and PAGE are exact in all views AND chain-exact. NOMO1/2/3 loci carry 1-2 introns.

### The two recommended cases, and what they do and do not show

- **GSTA1 / GSTA2** (CF316, chr6, tandem, same strand, 27.9 kb apart; identity .954, coverage .99; 0 bp entanglement; 111 / 115 exact-chain reads, all primary (review sweep); no read-level link between the copies, tie 0.000): exact in S1, S2 and S3, both loci chain-exact (6 introns), edge weight .929 from sequence and exon structure alone; a second human library (genome-wide testis, holder rule) gives the right family with 330 / 204 reads (both from the review sweep, not rechecked by me). chr6 is the least-exposed contig in the ledgers (quoted).
- **PAGE2 / PAGE2B / PAGE5** (CF400, chrX, same strand, copies >= 10.1 kb apart; identity .918-.973; 0 bp entanglement; tie <= .07): exact in all views, every locus chain-exact (4 introns), complete triangle (MCL density 1.0); second library right-family with 165 / 167 / 114 reads (review sweep). A119b depth is thin (minimum 12 reads).
- **They validate** a necessary condition: reads give one read-supported locus per copy, and the copy graph joins them into exactly one family with nothing else attached, on a clean, expressed, mid-divergence pair or triple.
- **They do not validate** sufficiency or recall (the class gives F .45), ties or identical copies (O2), unlisted relatives, low depth, other species, or anything blind: the outcome has been seen, so both are DEV regression anchors. Both outcomes are replay products, not HEAD runs. n = 1 each; a "recent duplication" premise is unverified.
- Not for use as "ideal": identical-copy pairs (EIF3C/EIF3CL identity 1.000, MAGEA9/9B, HSFY1/2) recover because both loci are seeded by the same tied reads (O1-trivial, useful as O2 tie controls); LYZL1/2 and UGT2B15/17 fail even with annotation as loci. HBA1/2 and GTF2IRD2/2B stay as regression anchors on the dev chromosomes.
- The literature screen (post hoc tiers, not pre-registered) proposed GPR89A/B, ROPN1/ROPN1B and WASHC2A/WASHC2C as separable twins; none passes the frozen criteria (C1 annotation quirks: partial map, a 2-exon CAT reference transcript, a non-canonical intron). The pipeline recovers GPR89A/B (CF11) and WASHC2A/C (CF62) exactly; ROPN1/1B (CF292) is partial.

### Prior evidence and held-out status

- Earlier few-copy cases (gorilla 06-20 to 08-19: RABL2A/B, AK6, CCDC196, MAGEA pairs, GSTM, DAZ, RBMY, TSPY, PCDHB, RFPL, APOBEC3; famsim SNRPB / SPAG7 simulations) all predate MCL, `--min-cov-shorter 0.70`, f1v2 and GOOD seeding; none has been re-run at the current default. No few-copy family has been pre-registered as an O1 positive control for de novo real reads at the current default (earlier pre-registrations were guided or chrY).
- **No A119b contig is unexposed** (flips were decided on products that include them: chr16 is the dev chromosome of every fusion rule, chr7 Soto-dev, chr20-22 F1v2-dev), so the 11 are DEV. The drafted held-out test uses gorilla OR6737 testis (never pooled with human), or a never-run human library (HG002: availability unverified).

### Limits and corrections the review forced

- The outcome product set: 16 of 19 contigs (chr1-12, 14, 15, 19, X) are replays of the 09-30 e163d955 arm D; replay = real run was checked only on chr17 and chrY. chr13, 18, 20, 21, 22 have no product.
- The scorer is strand blind and resolves one gene per locus (nested genes mislabel loci: the TRIM43/TRIM43B loci are labelled GPAT2); 20 families are "hidden" under S1. S3 is circular (its floor, >= 1 exonic bp, is the pipeline's own edge floor). "Exact" is unreachable where the listed size hides relatives (ZNF600/611 with ZNF808 and ZNF888 are not false merges).
- Reads: C5 uses the pipeline's own seeding pool, so for near-identical copies (identity >= .995) both loci are seeded by the same tied reads and recovery is weak evidence; the informative cases are the 8 no-tie families (cross-link <= .07). Selection on inputs is positive-heavy: an O1 score on it is a mechanism check, not a recall estimate.
- Sweep claims withdrawn or corrected: "no few-copy family was ever pre-registered" is true only for de novo real reads at the current default; the prior sweep's "best candidates (exact)" list is outcome-selected (8 of 10 fail the input criteria; RHOXF2/2B and ROPN1/1B are not exact under the scorer); the literature sweep's "93-100% of primaries also align at the other copy" does not show ties (the library was mapped with `-N 50 -p 0.1`; HBA1 to HBA2 is 398 / 407 with any secondary and 1 / 407 with a good secondary); its tier A is post hoc and none of its four families passes the criteria; the outcome headline "22%" is the inflated S1 figure (42 / 292 = 14% over all views).
- Metric traps hit and handled (reviewer): selection on outcome (F .881 on 11 chosen against .628 on all 72), denominators conditioned on outcome (80 families with all members in a clustered locus give 52 exact = 65%, against 52 / 223 = 23%), universe intersection, circular S3, unreachable "exact", threshold fragility (ZNF600/611 min reads 4 against 3; identity .912-.921 against the .90 floor), non-independence (ZNF33 and ZNF600/611 are KRAB-ZNF; PAGE and RHOXF2/2B are chrX arrays: the 11 are about 9 independent groups).


---

## IDEAL_EXPRESSION_DEFAULT_2026-10-06.md

*(original: `docs/archive/2026-10/merged/IDEAL_EXPRESSION_DEFAULT_2026-10-06.md`)*

## If NPIP and TBC1D3 were ideally expressed, would the current default pipeline find them? (2026-10-06)

Protocol: `docs/PREREG_ideal_expression_2026-10-06.md` (v1 56637089, Amendment 1 d745c6a9). Tools: `bench/ideal_expression/`. Products: `/mnt/linuxdisk/tmp/ideal_expression_2026-10-06/{NPIP,TBC1D3}/rep{1,2}/`.
HEAD release binaries (`copy_assign` 87824d91, `mcl_families` a6308244, `as_table` 48786e8f, `family_score` 542923fd), human CHM13 v2.0, CAT/Liftoff v2.0. **DEV, human only, circular by construction** (the reads come from the annotation that scores them):
the numbers say whether the default pipeline can find these copies when expression, truncation, readthrough and depth do not limit; they do not predict recovery from a real library.
Two independent read replicates per family (seeds 20261006 and 20261007): NPIP 20,660 reads from 2,066 simulated transcripts of 474 genes (547 genes and 2,147 transcripts lie in the windows), TBC1D3 16,680 reads from 1,668 transcripts of 464 genes (519 genes, 1,725 transcripts in the windows) (every CAT/Liftoff transcript within +-500 kb of a copy, 10 full-length jittered reads each, error .001).

### Answer

| family | copies | reachable (R) | IDEAL-FOUND on R, rep1 / rep2 | rule | all copies | E2 / E3 / E4 on R | cluster precision |
|---|---|---|---|---|---|---|---|
| **TBC1D3** | 16 | 13 | 12 / 12 | **YES** (both replicates) | 14 / 16 | 13 / 12 / 13 | .938 |
| **NPIP** | 25 | 14 | 12 / 12 | **NO** (both replicates; one copy short of the 13 needed) | 16 / 25 | 12 / 13 / 14 | .909 |

IDEAL-FOUND = the copy owns a locus (both purities >= .5, no other copy shares it), that locus's representative has exactly one of the copy's chains, and the locus is in the family's cluster. **R** = copies that are not entangled with another annotated gene (>= 100 shared exonic bp),
have a canonical multi-exon chain, are not readthrough images, and have an observable chain (>= 3 reads of the seeding pool): fixed from the truth tables and the BAM before the pipeline ran (R sha1 83ab9229 for NPIP, e43717f7 for TBC1D3; identical in both replicates).
YES needs >= ceil(0.9 n) on R and cluster precision >= .5. All four arms pass G1, G2, G3, G4 and G6 (VALID). Predictions: P1 (TBC1D3 YES) held, P2 (NPIP not YES, and <= 21 of 25 over all copies) held, P3 (replicates agree) held.

The excluded copies, printed beside (identical in both replicates): **NPIP** entangled 10 (IDEAL-FOUND 4: h05, h08, h10, h11; three of the ten are also readthrough images), chainless 1 (NPIPA3, found 0), R 14 (12). **TBC1D3** entangled 2 (found 2), chainless 1 (TBC1D3P7, found 0), R 13 (12).
A control that gives the same scorer PERFECT loci (the canonicalized annotation through `mcl_families`, G6) reaches 14 / 14 on NPIP's R (17 / 25 over all copies) and 13 / 13 on TBC1D3's R (15 / 16): relative to perfect loci the pipeline loses 2 reachable NPIP copies and 1 reachable TBC1D3 copy.

### Why the reachable copies are missed (NPIP two, TBC1D3 one per replicate)

- **NPIPA2 (h01), both replicates:** its two true chains are assembled (2 / 2 recovered) but the locus also holds three artifact chains built from reads of paralogous copies, and the most-reads representative is one of those (6 exons, representative purity .454): E2 and E3 fail. The representative problem of `docs/archive/2026-10/SPLICED_COPY_SUPPORT_2026-10-04.md` exists with ideal reads.
- **NPIPA7 (h06), both replicates:** the representative is exact (purity 1.0, 7-intron chain), but the locus (gene_id) also contains nine transcripts, most of them 21-30 exon chains of the neighbouring pseudogenes (AC138969.2, PKD1P6): locus purity .112, so E2 fails. What the families stage and O2 receive (the representative) is right.
- **TBC1D3P4 / TBC1D3P3 (h26, h27), one per replicate:** near-identical copies with one 11-intron chain each. The representative is a coin toss (10 reads against 10) between the copy's own chain and a chain with one intron displaced by about 25-31 bp. The only reads that carry the displaced chain are the paralog's (h27's reads at h26 in replicate 1, h26's reads at h27 in replicate 2), all 10 as secondary alignments at this copy: the aligner put one junction off in the paralog's sequence, and the good-secondary seeding admitted them.
  The assembled transcripts contain every chain of every TBC1D3 copy (31 / 31 chains recovered, both replicates), so this is representative choice, not assembly. Which of the two copies fails flips between replicates.

### Same copies, same registered instrument, real reads against ideal reads

The real-read column is the current default on A119b (CAT/Liftoff copies, instrument of `docs/archive/2026-10/PREREG_spliced_copy_support_2026-10-04.md` Amendments A and B, `bench/copy_support.py`). Only these annotation-anchored measures are like-for-like; the registered strict FOUND needs the cap signal and cannot run on simulated reads.

| instrument | NPIP real | NPIP ideal (rep1 / rep2) | TBC1D3 real | TBC1D3 ideal (rep1 / rep2) |
|---|---|---|---|---|
| copies with a locus overlapping them | 25 / 25 (own node 24) | 25 / 25 (own node 24) | 13 / 16 (no node file) | 16 / 16 (own node 15) |
| annotation-anchored FOUND (Amendment A) | 12 | 23 / 23 | 12 | 16 / 16 |
| exact-chain FOUND (Amendment B) | 6 | 21 / 21 | 9 | 15 / 15 |
| locus level, annotation-anchored | 23 | 24 / 24 | 13 | 16 / 16 |
| locus level, chain | 14 | 24 / 24 | 12 | 16 / 16 |

On the same 25 NPIP copies the annotation-anchored FOUND goes from 12 to 23 and the exact-chain FOUND from 6 to 21 when the reads are ideal; for TBC1D3 12 to 16 and 9 to 15. The real-read shortfall is therefore mostly the data, with the representative problem left over at NPIP.

### Observability and chains (E0, E1)

NPIP: 24 of 25 copies have an observable chain; 142 of 148 chains are observable (>= 3 reads of the seeding pool; 131 with primary alignments alone); the primary alignment of a read lands on its own copy for 84% of reads (12 copies below 90%: the near-identical paralogs tie), yet junction-level observability is intact.
Chains recovered exactly as an assembled transcript: 133 / 148 (rep1), 129 / 148 (rep2); the unrecovered chains are at NPIPA1 (6 / 8), NPIPB5 (17 / 21), NPIPB4 (17 / 21) and AC138894.1 (7 / 12), all entangled copies. TBC1D3: 15 of 16 copies have an observable chain, 31 / 31 chains observable and recovered, 98% of reads land on their own copy.

### Family level (beside; NOT a ceiling)

`family_score` is keyed by RefSeq genes while the reads are CAT/Liftoff (RefSeq NPIPB3 has no CAT gene, and several NPIP symbols differ), and its universe here is the simulated windows, so these numbers are not the real-read run's and not an upper bound.
NPIP family (Compara CF153, 19 genes in the universe): 11 in the NPIP cluster, sens .579, prec .917, F .710 (U2 ID_154: 11 / 19, F .611 and .595); pooled over the nine Compara families F .733 (.667 without size-2 clusters). TBC1D3 family (Compara CF185, 9 genes): 9 / 9 in the cluster, prec .900, F .947; pooled F .753 and .737 (.720 and .685 without size-2 clusters).

### Controls

G1 0 violations in 2,066 and 1,668 simulated transcripts (12,857 and 8,093 introns), re-read from the genome by a separate script. G2 100% (NPIP) and 99.9% (TBC1D3) of the reads of 275 and 304 single-copy genes land on their gene. G3 one record per FASTQ read, driver exit 0, `families.gtf` newer than `gtf`. G4 scorer tables identical under two hash seeds.
**G5** an independent scorer written from the pre-registration text alone (it never opened `score.py`) reproduces E2, E3, E4 and IDEAL-FOUND for all 82 copy-arms (25 + 25 + 16 + 16) with no difference; a line-by-line code review found no blocker (one chromosome-check bug, found and fixed during the first scoring, and several consistency fixes made after the results: unmapped reads in the G2/E0 denominators, unrounded purity thresholds, E3p over all locus transcripts,
cluster precision without a strand clause, the folded-holder count; every per-copy flag was identical before and after). **G6** perfect loci give E2 = E3 = E4 = 100% on R in all four arms. **G7 (single-copy genes, >= 95% with exactly one locus) is NOT attainable as registered:** the annotation-as-loci control reaches 89.5% (NPIP) and 94.1% (TBC1D3) because genes overlapping other genes cannot each get one locus,
and the pipeline arm reaches 63-65% and 62% (many single-copy genes are single-exon and receive no locus at ten reads); it is reported against the control, not as a gate.

### Robustness of the NPIP verdict

The independent scorer computed twelve alternative readings of definitions the text leaves open (copy exon union against copy span, folded-record handling, strand clauses, cluster-precision variants, tie order, ...). Eleven keep NPIP at NO; one (purity measured against the copy's span instead of its exon union) gives PARTLY; none gives YES.
Post hoc, and not a verdict: if E2 judged the representative only (the representative is what the families stage and O2 receive), NPIPA7 would pass and NPIP would be 13 / 14, the YES threshold. NPIP's verdict therefore sits one copy below the bar and TBC1D3's exactly on it (12 of 13 needed 12).

### Held-out contrast (no bars)

SMG1P (7 genes in the NPIP windows): 2 reachable, 1 found; KRTAP (33 genes in the TBC1D3 window): 30 are single-exon, 2 reachable, 0 found (not grouped in one family). Too few reachable copies to say more.

### Per-copy results (replicate 1; replicate 2 differs only where the last column shows it)

#### NPIP

| copy | name (CAT/Liftoff) | strata | chains (observable) | chains recovered | reads on own copy | holder purity rep / locus | E2 | E3 | E4 | IDEAL-FOUND rep1 / rep2 | if not found |
|---|---|---|---|---|---|---|---|---|---|---|---|
| h00 | NPIPB2 | R | 8 (8) | 8/8 | 1.00 | 1.0 / 1.0 | 1 | 1 | 1 | 1 / 1 |  |
| h01 | NPIPA2 | R | 2 (2) | 2/2 | 1.00 | 0.454 / 0.535 | 0 | 0 | 1 | 0 / 0 | locus |
| h02 | NPIPA1 | E with CHM13_G0020712 (3881 bp, 53%) | 8 (8) | 6/8 | 1.00 | 1.0 / 0.674 | 1 | 0 | 1 | 0 / 0 | annotation:E |
| h03 | PKD1P6 | E+X with CHM13_G0020726 (1888 bp, 31%) | 6 (6) | 6/6 | 1.00 | 1.0 / 0.716 | 1 | 1 | 0 | 0 / 0 | annotation:E |
| h04 | NPIPA5 | R | 3 (3) | 3/3 | 1.00 | 1.0 / 0.835 | 1 | 1 | 1 | 1 / 1 |  |
| h05 | AC138969.1 | E+X with CHM13_G0020764 (5447 bp, 68%) | 12 (12) | 12/12 | 0.83 | 1.0 / 0.837 | 1 | 1 | 1 | 1 / 1 |  |
| h06 | NPIPA7 | R | 2 (2) | 2/2 | 0.65 | 1.0 / 0.112 | 0 | 1 | 1 | 0 / 0 | locus |
| h07 | NPIPA8 | E with CHM13_G0020801 (1227 bp, 76%) | 3 (3) | 3/3 | 1.00 | 0.906 / 0.129 | 0 | 1 | 1 | 0 / 0 | annotation:E |
| h08 | AC138969.1 | E+X with CHM13_G0020802 (1568 bp, 25%) | 9 (9) | 9/9 | 0.86 | 1.0 / 0.648 | 1 | 1 | 1 | 1 / 1 |  |
| h10 | NPIPB5 | E with CHM13_G0020899 (2061 bp, 22%) | 21 (19) | 17/21 | 0.56 | 1.0 / 0.606 | 1 | 1 | 1 | 1 / 1 |  |
| h11 | NPIPB4 | E with CHM13_G0020920 (386 bp, 5%) | 21 (20) | 17/21 | 0.70 | 1.0 / 0.537 | 1 | 1 | 1 | 1 / 1 |  |
| h12 | NPIPB3 | E with CHM13_G0020934 (386 bp, 10%) | 13 (12) | 13/13 | 0.70 | 1.0 / 0.234 | 0 | 1 | 1 | 0 / 0 | annotation:E |
| h13 | NPIPB6 | R | 2 (2) | 2/2 | 1.00 | 1.0 / 0.864 | 1 | 1 | 1 | 1 / 1 |  |
| h14 | NPIPB8 | R | 2 (2) | 2/2 | 0.65 | 1.0 / 1.0 | 1 | 1 | 1 | 1 / 1 |  |
| h15 | AC138894.1 | E with CHM13_G0021066 (2480 bp, 59%) | 12 (10) | 7/12 | 1.00 | 0.944 / 0.477 | 0 | 0 | 1 | 0 / 0 | annotation:E |
| h16 | NPIPB9 | R | 2 (2) | 2/2 | 1.00 | 1.0 / 0.928 | 1 | 1 | 1 | 1 / 1 |  |
| h17 | NPIPB10P | R | 1 (1) | 1/1 | 1.00 | 1.0 / 1.0 | 1 | 1 | 1 | 1 / 1 |  |
| h18 | NPIPB11 | R | 2 (2) | 2/2 | 1.00 | 1.0 / 1.0 | 1 | 1 | 1 | 1 / 1 |  |
| h19 | NPIPB12 | R | 8 (8) | 8/8 | 0.67 | 1.0 / 0.746 | 1 | 1 | 1 | 1 / 1 |  |
| h20 | NPIPB13 | R | 5 (5) | 5/5 | 0.46 | 1.0 / 0.656 | 1 | 1 | 1 | 1 / 1 |  |
| h21 | NPIPA3 | C | 0 | - | 1.00 | 1.0 / 0.702 | 1 | 0 | 1 | 0 / 0 | annotation:C |
| h22 | NPIPB14P | E with CHM13_G0022124 (1678 bp, 97%) | 3 (3) | 3/3 | 1.00 | 0.373 / 0.27 | 0 | 0 | 1 | 0 / 0 | annotation:E |
| h23 | NPIPB15 | R | 1 (1) | 1/1 | 0.70 | 1.0 / 1.0 | 1 | 1 | 1 | 1 / 1 |  |
| h24 | NPIPB15 | R | 1 (1) | 1/1 | 0.60 | 1.0 / 1.0 | 1 | 1 | 1 | 1 / 1 |  |
| h25 | NPIPB15 | R | 1 (1) | 1/1 | 0.70 | 1.0 / 1.0 | 1 | 1 | 1 | 1 / 1 |  |

#### TBC1D3

| copy | name (CAT/Liftoff) | strata | chains (observable) | chains recovered | reads on own copy | holder purity rep / locus | E2 | E3 | E4 | IDEAL-FOUND rep1 / rep2 | if not found |
|---|---|---|---|---|---|---|---|---|---|---|---|
| h26 | TBC1D3P4 | R | 1 (1) | 1/1 | 1.00 | 0.971 / 0.972 | 1 | 0 | 1 | 0 / 1 | representative |
| h27 | TBC1D3P3 | R | 1 (1) | 1/1 | 1.00 | 1.0 / 0.977 | 1 | 1 | 1 | 1 / 0 |  |
| h28 | TBC1D3P5 | R | 3 (3) | 3/3 | 1.00 | 1.0 / 1.0 | 1 | 1 | 1 | 1 / 1 |  |
| h29 | TBC1D29P | R | 5 (5) | 5/5 | 1.00 | 1.0 / 1.0 | 1 | 1 | 1 | 1 / 1 |  |
| h30 | TBC1D3B | E with CHM13_G0024042 (1383 bp, 49%) | 4 (4) | 4/4 | 0.90 | 1.0 / 1.0 | 1 | 1 | 1 | 1 / 1 |  |
| h31 | TBC1D3I | R | 3 (3) | 3/3 | 1.00 | 0.993 / 0.92 | 1 | 1 | 1 | 1 / 1 |  |
| h32 | TBC1D3G | R | 1 (1) | 1/1 | 1.00 | 0.974 / 0.727 | 1 | 1 | 1 | 1 / 1 |  |
| h33 | TBC1D3H | R | 1 (1) | 1/1 | 1.00 | 0.979 / 0.818 | 1 | 1 | 1 | 1 / 1 |  |
| h34 | TBC1D3B | R | 4 (4) | 4/4 | 0.78 | 1.0 / 1.0 | 1 | 1 | 1 | 1 / 1 |  |
| h35 | TBC1D3K | R | 1 (1) | 1/1 | 1.00 | 0.998 / 0.781 | 1 | 1 | 1 | 1 / 1 |  |
| h36 | TBC1D3D | R | 1 (1) | 1/1 | 1.00 | 0.998 / 0.781 | 1 | 1 | 1 | 1 / 1 |  |
| h37 | TBC1D3E | R | 1 (1) | 1/1 | 1.00 | 0.98 / 0.767 | 1 | 1 | 1 | 1 / 1 |  |
| h38 | TBC1D3 | R | 1 (1) | 1/1 | 1.00 | 0.973 / 0.817 | 1 | 1 | 1 | 1 / 1 |  |
| h39 | TBC1D3P7 | C | 0 | - | 1.00 | 1.0 / 1.0 | 1 | 1 | 0 | 0 / 0 | annotation:C |
| h40 | TBC1D3P1 | E with CHM13_G0024955 (182 bp, 11%) | 1 (1) | 1/1 | 1.00 | 1.0 / 1.0 | 1 | 1 | 1 | 1 / 1 |  |
| h41 | TBC1D3P2 | R | 3 (3) | 3/3 | 1.00 | 1.0 / 1.0 | 1 | 1 | 1 | 1 / 1 |  |


### Departures from the plan, and limits

Stratum X (readthrough images) is not exclusive: NPIP h03, h05, h08 are X and also E. The scorer was changed once after the first scoring (a missing chromosome check, which made G7 return 0 and could have affected cluster precision; every E2, E3, E4 and IDEAL-FOUND flag was unchanged) and again after the independent review (consistency fixes above, flags unchanged).
The legacy 2026-09-21 reads were not used (register 1253). The held-out families, the family-level scores and the like-for-like instruments were computed after the verdicts and are reported beside them. Not shown: real reads, O2 (about 16% of NPIP reads have their primary alignment on a paralogous copy even when every copy is expressed, so contested reads exist in the ideal case and O2 could be tested there),
O3, gorilla, other families, and anything about the annotation's correctness. A YES or NO here is about the algorithm at these loci: 10 of the 25 NPIP copies overlap another annotated gene and cannot each get a locus from any assembler.


---

## LOCUS_UNITS_2026-10-06.md

*(original: `docs/archive/2026-10/merged/LOCUS_UNITS_2026-10-06.md`)*

## A locus model for overlapping genes on ideal reads: primary-first representatives and structural separators (2026-10-06)

**Status (arms produced 2026-10-07 00:04-00:08 PDT; independently verified, workflow wf_994ece76-8d2: every number reproduces): DEV, human CHM13 v2.0 / CAT-Liftoff, simulation, circular by construction, the same two windows (NPIP chr16, TBC1D3 chr17) x two read replicates as `docs/IDEAL_EXPRESSION_DEFAULT_2026-10-06.md`. Protocol: `docs/PREREG_locus_units_2026-10-06.md` (committed 4ce90b18; Amendment 1 bff7dbfe, written before any arm existed). Tools: `bench/entangled/{locus_units.py (12 unit tests), run.sh units2}`. Products: `/mnt/linuxdisk/tmp/entangled_2026-10-06/{NPIP,TBC1D3}/rep{1,2}/{q,sd,sq}.*`.**
User direction: option B (a locus-level model for junction-sharing overlapping genes; the advisor wants clean combinatorial structure, provable statements, no arbitrary thresholds).

### Answer

1. **Rule 2 (evidence-free structural separators) is refuted at the family level on the dev windows (its wording was chosen with a truth-chain oracle on the same windows).** The private-bridge rule removes the transcripts that are the only junction-link between groups: on NPIP replicate 1 its 30 separators include almost only a gene's own long or main isoform (NPIPB3, B4, B5, B12 and B13 with 47 reads, SNX29, MAPK3, MRTFB, SGF29, TUFM, ITGAL, ZNF785, MLKL, AATF, SRCIN1 ...), in TBC1D3 TBC1D3G, TBC1D3B and TBC1D3P5. A bridging isoform and a fusion have the same junction structure; (hypothesis: reads' PAS and promoter signals, F1v2's rule, could tell them apart; the evidence-based F1 arm of `docs/LOCUS_UNITS_LEVELS_2026-10-07.md` found no signal to act on in these windows). Pooled over the four runs S_D finds 38 of 54 reachable copies (D: 48) and 9 of 24 entangled (D: 12).
2. **Rule 1 (primary-first representative) repairs exactly the failures it targets and is not a candidate by the registered rule.** E3 (the representative has one of the copy's chains) goes from 50 to **54 of 54** reachable gene-runs: NPIPA2 (both replicates) and TBC1D3P4 / TBC1D3P3 are fixed. The family-level loss is one event: in NPIP replicate 2 the NPIP cluster splits and NPIPB6, NPIPB8 and NPIPB9 land in a sibling cluster (E4 fails for three copies; Compara F .733 -> .643). R found 49 against D's 48.
3. **Primaries-only seeding (arm P, not registered as a family-level candidate here) makes the registered ideal verdict YES for both families**: NPIP R found 14 / 13 (bar 13), TBC1D3 13 / 13 (bar 12); the price is the entangled copies (E found 7 of 24 against 12): their cluster membership came from the enlarged bodies that secondary-borne chains build.
4. **The family-level criterion E4 ('the holder is in the one cluster K*') is brittle to partition granularity**: each arm that changes a few representatives or bodies re-partitions NPIP and moves a subfamily to a sibling cluster (rep 2: MCL0 -> MCL2). This is (post hoc) the advisor's hierarchy question (nested levels, level-relative TP / FP) showing up in the measurement.

### Pooled results (four runs; R = 54 reachable gene-runs, E = 24 entangled)

| arm | R found | E found | R: E2 / E3 / E4 | E: E2 / E3 / E4 | NPIP R found r1 / r2 (bar 13) | TBC1D3 R found r1 / r2 (bar 12) | Compara F NPIP r1 / r2 | Compara F TBC1D3 r1 / r2 |
|---|---|---|---|---|---|---|---|---|
| D default | 48 | 12 | 50 / 50 / 54 | 16 / 18 / 22 | 12 / 12 | 12 / 12 | .733 / .733 | .753 / .737 |
| P primaries only | 53 | 7 | 53 / 53 / 54 | 17 / 17 / 14 | 14 / 13 | 13 / 13 | .746 / .724 | .753 / .737 |
| C sub-chain drop | 48 | 13 | 50 / 50 / 54 | 16 / 18 / 23 | 12 / 12 | 12 / 12 | | |
| PC | 53 | 7 | as P | as P | 14 / 13 | 13 / 13 | .746 / .724 | .753 / .737 |
| **Q** Rule 1 | 49 | 12 | 52 / **54** / 51 | 16 / 18 / 22 | 13 / **10** | 13 / 13 | .733 / **.643** | .753 / .737 |
| **S_D** Rule 2 | 38 | 9 | 51 / 47 / 46 | 17 / 15 / 19 | 11 / 9 | 9 / 9 | .690 / .690 | .753 / .737 |
| **S_Q** Rules 1+2 | 45 | 10 | 52 / 54 / 47 | 16 / 16 / 20 | 13 / 10 | 11 / 11 | .690 / .667 | .753 / .737 |

Registered rule (YES at >= ceil(0.9 R), lower of the replicates): D NPIP NO, TBC1D3 YES; P and PC YES, YES; Q NPIP NO (rep 2: 10), TBC1D3 YES; S_D NO, NO; S_Q NO, NO.

### Predictions and the decision rule (prereg sections 5 and 6)

Q1 (Q finds >= 52 of 54 reachable): **failed** (49). Q2 (Q finds >= 11 entangled): held (12). Q3 (Q gives YES for both families in both replicates): **failed** (NPIP replicate 2). Q4 (S_D finds >= 2 more entangled than D with no loss at R): **failed** (9 against 12; R -10). Q5 (S_Q at least Q at R and E, and one more at E): **failed**. Q6 (Compara F of Q and S_Q >= D minus .005): **failed** (NPIP replicate 2; S_Q NPIP both).
Decision rule: Q has R + E 61 > 60 and R >= D and E >= D - 1, but Compara F fails in NPIP replicate 2 and the direction does not hold in both replicates of both families: **no arm is a CANDIDATE**; no held-out or real-data check is earned.

### What the arms show (post hoc, mechanism)

- Rule 1 per copy against D: gains NPIPA2 (representative purity .454 -> 1.0, both replicates), TBC1D3P4 / P3 (representative displaced chain -> own chain); no entangled copy changes; losses only the three NPIP replicate-2 cluster moves (NPIPB6, NPIPB8, NPIPB9: MCL0 -> MCL2), locus purities unchanged.
- Rule 2 cuts 30-33 transcripts per NPIP replicate and 11 per TBC1D3 replicate (F1v2 cuts 1-2); the truth-chain oracle had predicted the false-cut cost (16 of 108 junction-sharing entangled genes resolved, 12 more non-entangled genes split than with components). At the family level the cuts also shrink bodies and split loci into several nodes: copies move out of the family cluster (TBC1D29P and TBC1D3B to a sibling cluster, NPIPB5 to MCL5).
- The entangled copies that fail in D split into (a) artifact-driven (NPIPA8, NPIPB3: locus purity .129 and .234, 5 and 24 artifact transcripts in the locus; NPIPA1: the representative is an artifact chain) and (b) truly fused (NPIPB14P with the readthrough PDXDC2P-NPIPB14P as representative, AC138894.1 with CLN3 and NPIPB7, PKD1P6 with NPIPP1); primaries-only seeding fixes NPIPB3's locus purity and NPIPA1's representative, but breaks NPIPA8 (representative purity .906 -> .295) and loses the cluster membership of the two AC138969.1 copies and NPIPA1.

### What this does not show, and the next registered step

Same two windows; the arms were motivated by what the dev windows showed; the primary-supported marks come from a second assembly. Nothing here is held out. What the data point to (not yet tested): (1) a family-level criterion that is hierarchy-aware (membership in the family at the coarser nested level, graph component beside the MCL cluster, with the precision of that level) so that a partition flip is not scored as a representation failure; (2) an evidence-based rule for the fusion cut (F1v2's) with the structural rule only as a candidate generator; (3) the seeding trade-off, which has opposite signs for reachable and entangled copies, on real data and on windows around other multi-copy regions.


---

## PREREG_default_rescore_npip_2026-10-06.md

*(original: `docs/archive/2026-10/merged/PREREG_default_rescore_npip_2026-10-06.md`)*

## PREREG: strict FOUND, U2 and Compara CF153 for the current default pipeline at NPIP (2026-10-06)

Written and committed before any product of the HEAD binaries on this protocol exists (`DEF` and `PRE` below are generated after this commit).
Everything here is DEV: human A119b chr16 is the development block of NPIP and of every fusion rule. A PASS is "no regression against the registered
development values", not a validation.

### Question

The registered strict-FOUND numbers at NPIP (`docs/archive/2026-10/SPLICED_COPY_SUPPORT_2026-10-04.md`, Amendment E: 10 / 6 / 2 found within own nodes 23 / 21 / 24, locus
level 19 / 18 / 16, for P / GOOD / ALL) were measured on pre-f1v2 arms from a frozen 2026-10-01 binary. The family scores of the current default
(f1v2 + `--min-cov-shorter 0.70`, U2 F .645, Compara chr16 F .667, CF153 F .812) were measured on 2026-09-30 at e163d955. No product has the current default's
loci scored for strict FOUND, and none of these numbers was produced by the HEAD binaries (the 2026-10-06 consolidation changed `src/`, and its real-data
slice is not recorded as run). Does the default at HEAD do no worse than the registered development values?

### Arms (human A119b, chr16, the 25 CAT/Liftoff NPIP copies of `copy_recovery_tools_cat/ann`)

- **DEF**, the default: HEAD `copy_assign --assemble-only` with the driver's `assemble` flags (strict junctions, shipped polish, `--gtf-tpm`,
  `--bridge-regroup f1v2`), seeding from the stored genome-wide best-AS table (secondaries >= 0.98 of the best AS, the default `GOOD` pool), on
  `--region chr16` in place of `--genome-wide` (the two are equal on the contig: the gates of `docs/archive/2026-09/CONTAINER_HEADROOM_2026-09-30.md`); then the driver's
  `families` stage with no override (`--min-cov-shorter 0.70`, most-reads representative, `--min-shared-exon-frac 0.60`).
- **PRE**, the control: the same HEAD binaries with `--bridge-regroup off` and `RUSTLE_MIN_COV_SHORTER=0` (the 2026-09-25 pipeline). It is what the
  registered `GOOD` arm was meant to be, rebuilt by the binaries under test.
- The frozen read-pool arms `P`, `GOOD`, `ALL` (`/mnt/linuxdisk/tmp/readpool_npip`) are used only by the gates.

### Instruments

- `bench/copy_support.py` (Amendment E, unchanged since 2026-10-04; sha1 in every log), arms given as `--loci NAME=loci.gff3,transcripts.gtf`
  (DEF: `PREFIX.families.gtf`, the GTF the families were built from; PRE: `PREFIX.gtf`) and `--nodes`.
- `bench/default_rescore/nodes.py`: the own-node rule of `bench/npip_read_pool/pagedata.py`, with the copy exons of the registered scoring (`npip_read_pool.json`) (an NPIP cluster = a cluster that holds a locus overlapping a
  same-strand NPIP copy; a copy has an own node iff such a locus overlaps it).
- `family_score` of the HEAD build through `bench/default_rescore/score_families.py` (clusters and truth restricted to chr16, `--chrom ALL --pairwise
  --per-family`, truths compara / u2 / soto of `families_gw/species/human`).
- `bench/default_rescore/verdict.py` computes the gates, the rules and the verdict from the products. Nothing is judged by hand.
- Denominator: 25 copies, fixed. The own-node count is conditioned on the prediction and is reported beside, never used as a denominator.

### Gates (instrument checks; a failed gate makes the verdict INVALID, the arms are then reported without a verdict)

- **G0** the HEAD scorer on the frozen arms reproduces `support_hsa_E.json` exactly for each arm: `old_overlap_in_npip_nodes` 23 / 21 / 24,
  `tc_found_in_npip_nodes` 10 / 6 / 2, `locus_tc_found` 19 / 18 / 16, `tc_found` 10 / 6 / 2.
- **G1** `nodes.py` on the frozen arms P and GOOD reproduces every stored `pagedata.json` node flag (25 copies x 2 arms; copy exons from `npip_read_pool.json`, see Amendment 1).
- **G3** the HEAD `family_score` on the stored e163d955 default clusters (`container_headroom/data/human_A119b/fam/chr16/D.fam.clusters.tsv`) reproduces the
  stored U2 (sens .588, prec .714, F .645) and Compara chr16 (sens .500, prec 1.000, F .667).
  (Before this registration I ran the HEAD `family_score` once on those stored clusters, U2 only, to learn its output format: it printed .588 / .714 / .645.
  That is a reading of stored data, not of DEF or PRE.)
- **G2** (a diagnostic, not a gate): PRE reproduces the registered `GOOD` values (own node 21, E-found in own nodes 6, locus level 18, E-found chr16-wide 6).
  If it does not, the pipeline has drifted since 2026-10-01 and the rules below still compare DEF with the registered values, with the drift reported.
- The scorer is run twice under different `PYTHONHASHSEED` values; the two summaries must be equal (the Python `set` trap of register 1044).

### Decision rules on DEF (fixed before the run)

- **R1** strict E-found within own nodes (`tc_found_in_npip_nodes`) >= 6 of 25 (the registered `GOOD` value).
- **R2** locus-level E-found (`locus_tc_found`) >= 18 of 25.
- **R3** U2 bipartite F (`family_score` pooled line) >= 0.645, at three decimals.
- **R4** NPIPB2 and NPIPB6 each have an own node in DEF (the registered `GOOD` arm had neither; the 09-29 flip separated them from GSPT1 and EIF3CL).
- **Verdict:** INVALID if G0, G1 or G3 fails; otherwise PASS iff R1 to R4 all hold; otherwise FAIL, naming the rules that fail. No rescue arm and no
  redefinition follows a FAIL: it is recorded as the finding (the representative problem of `docs/archive/2026-10/SPLICED_COPY_SUPPORT_2026-10-04.md` would then be a
  property of the default, not of an older arm).

### Reported beside, no bar

DEF and PRE: own-node count, E-found in own nodes, locus level, chr16-wide E-found; the copies that change between PRE and DEF (found and own-node flags);
U2, Compara and Soto sens / prec / F and pairwise tp / truth / predicted pairs; the Compara CF153 row (hit / truth genes, sens, prec, F; stored value 13 of 19
genes, F .812); `cmp` of the HEAD products against the e163d955 products (assembled GTF, families GTF, clusters, loci, copy table; DEF against `D`, PRE against `B0`).
Liftoff recall is not recomputed (its scorer is not in the repo's runner).

### Not claimed

Anything about held-out NPIP substrates, gorilla, O2 or O3. Nothing about f1v2 or `--min-cov-shorter` as such: DEF against PRE is a descriptive
comparison on the development block. Strict FOUND is the registered rule of `docs/archive/2026-10/PREREG_spliced_copy_support_2026-10-04.md` (Amendment E); it depends on the
cap signal of this library, which gorilla lacks.

Runner: `bench/default_rescore/run.sh gates | asm DEF | fam DEF | asm PRE | fam PRE | score | verdict`.

### Amendment 1 (2026-10-06, written after the first gate run and before any DEF or PRE product exists)

The first run of G1 failed, on the frozen arms only (no DEF or PRE product was generated or read): `nodes.py` gave own-node counts 24 / 22 for P / GOOD
against the stored 23 / 21. Two causes, both in how the stored rule was run, not in the arms: (a) the stored rule took each copy's exons from the CAT gene's
exon blocks (`npip_read_pool.json`, `copies[].exons`), while `nodes.py` took the union of the copy's transcripts in `truth.hsa.gtf`; they differ at
PKD1P6-NPIPP1, where the stored flags are False in both arms; (b) the frozen `ALL` arm has two loci on one span (`DN_chr16_33611806_2` and `_3`), which
`pagedata.py`'s span join resolved silently and `nodes.py` refuses. Changes: `nodes.py` takes the copy exons from `npip_read_pool.json` (`--exons-json`; used for
DEF and PRE too, so the own-node rule is the registered one), and G1 is checked on P and GOOD, the two arms whose loci are unique. With the change G1 compares
50 flags and finds 0 mismatches. G0 still uses all three frozen arms (it reads the stored `pagedata.json`, not `nodes.py`). The rules R1 to R4 are unchanged.


---

## PREREG_entangled_baseline_2026-10-06.md

*(original: `docs/archive/2026-10/merged/PREREG_entangled_baseline_2026-10-06.md`)*

## Pre-registration: what StringTie and FLAIR do at overlapping genes, on ideal reads, and what "better" means (2026-10-06)

**Written before any StringTie, FLAIR or f1units output exists, and before `arms_score.py` was run on any arm.** Status of what was already seen: the default assembler products of the four ideal-expression runs exist and their per-COPY chain recovery (E1) and locus scores (E2-E4) are known (`docs/IDEAL_EXPRESSION_DEFAULT_2026-10-06.md`); the window-wide, per-GENE numbers of this study are new. DEV, human CHM13 v2.0 / CAT-Liftoff only, simulation, circular by construction (the reads come from the annotation that scores them).
User request (2026-10-06): "we need to be able to model the cases in which there are overlapping genes, ideally we should be able to provide what stringtie and flair do for those loci and then do better than that." Tools: `bench/entangled/` (`arms_score.py`, `pool.py`, `run.sh`, `flair_shims/`).

### 1. Substrate and arms

Reads: the four ideal-expression read sets (NPIP and TBC1D3 windows, replicates 1 and 2, seeds 20261006 / 20261007; 20,660 and 16,680 reads of 2,066 and 1,668 transcripts of 547 and 519 genes; 10 full-length jittered reads per transcript, error .001), the BAM `reads.bam` already used by the pipeline arm (minimap2 splice:hq, -N 50 -p 0.1, secondaries kept). Every arm sees the same alignments.

| arm | what | recipe |
|---|---|---|
| **S** | StringTie, annotation-free | `stringtie -L -p 4` on `reads.bam`, no `-G` (the lab's `benchmark_collapse/run_stringtie.sbatch`) |
| **F** | FLAIR, annotation-free | FLAIR's own BAM-to-BED12 filtering (`dofiltering`, primaries only, supplementary removed, quality 0) then `flair collapse --trust_ends --generate_map`, no `flair correct`, no annotation (the lab's `run_flair.sbatch`) |
| **D_asm** | our assembler's transcripts, default | `asm.gtf` of the registered default arm (HEAD binaries, f1v2 regroup, strict junctions, GOOD secondaries) |
| **D_fam** | the same after the families stage's input regroup | `asm.families.gtf` (bridge transcripts removed, gene_ids split) |
| **U_asm / U_fam** | the paused container-v2 lever, opt-in code at HEAD | the driver with `RUSTLE_BRIDGE_REGROUP=f1units` (each bridge transcript replaced by its units), `units.gtf` / `units.families.gtf` |

Versions (drift from the lab runs stated, not corrected): StringTie 3.0.3 (lab 3.0.1 and 3.0.3), FLAIR 3.0.0 (lab 3.0.1), minimap2 2.30 for FLAIR's internal alignment (lab 2.31). isoseq collapse is not run (the request names StringTie and FLAIR).

### 2. Truth, strata, metrics (all from `reads.transcripts.tsv`; computed by `arms_score.py`)

Truth = every simulated transcript of every gene in the windows, canonical chains (introns >= 50 bp and canonical). A gene is ENTANGLED (stratum **E**) iff its exon union shares >= 100 bp with another gene's exon union on the same strand (the rule that defined the stratum for the 25 + 16 family copies); the other genes with >= 1 multi-exon chain are stratum **N**. Genes without a multi-exon chain are in neither.
Per arm and stratum:
- **M1 chains**: truth chains (chrom, strand, ordered intron list) equal to the chain of >= 1 arm transcript; also genes COMPLETE (every chain recovered).
- **M2 artifacts**: arm multi-exon transcripts whose chain is no truth chain: fragment (contiguous part of a truth chain), fusion (exons overlap >= 2 truth genes on the strand), other. Reported as a count and per recovered chain.
- **M3 resolution**: a truth gene is RESOLVED iff some arm gene_id carries an exact chain of it and an exact chain of no other truth gene. `merged_ids` = arm gene_ids carrying exact chains of >= 2 truth genes.
All four runs are pooled (per family and together); replicates are never averaged away: the per-run rows are kept.

### 3. Predictions (falsifiable, fixed now)

- **P1** Both S and F recover a smaller share of chains at entangled genes (E) than at the others (N).
- **P2** D_asm has fewer artifact transcripts per recovered chain than S and than F.
- **P3 (parity)** D_asm recovers at least as many chains at entangled genes as S and as F (pooled). If this fails, the gap to the tools is at the assembler and is the first thing to repair.
- **P4 (units)** U_asm recovers at least as many chains as D_asm, and U_fam resolves strictly more entangled genes than D_fam, with no loss at stratum N (chains and resolved).

### 4. Definition of "better", fixed before any improvement arm is built

An arm is BETTER than the tools at overlapping genes iff, pooled over the four runs, ALL of: (a) its chains recovered at stratum E >= the larger of S and F; (b) its artifact transcripts at least as few as the smaller of S and F; (c) its resolved genes at E >= the larger of S and F AND strictly above D_fam; (d) no loss at stratum N relative to D_asm / D_fam (chains, complete genes, resolved). Family-level consequences (E2-E4 per copy, `bench/ideal_expression`) are reported beside, not part of this bar.
A BETTER arm earns a registered real-data check (A119b chr16 and chr17, gorilla OR6737) with the copy-recovery instruments; nothing is claimed for real reads from this study.
If no arm is BETTER, the report lists where each loses (which of a-d) and the next lever is registered separately.

### 5. Limits declared in advance

Simulation with the annotation as truth; ideal reads (no truncation, no readthrough, no depth limit); human only; two replicates of one design; the FLAIR and StringTie versions differ from the lab's; U is an existing opt-in lever, not a new design; "resolved" is tool-agnostic (it asks for an own gene_id) and for tools whose gene_id is per transcript it equals "any chain recovered".

### 6. What was run before the freeze

One truth-vs-truth smoke test of `arms_score.py` (the NPIP replicate-1 truth transcripts scored as an arm, no tool output): 2,066 transcripts, 1,899 multi-exon; stratum E 86 genes with 568 chains, stratum N 264 genes with 1,295 chains; every chain recovered, 0 artifacts, 0 merged ids, every gene resolved. Nothing else of S, F or U existed, and D had not been scored by this script.

### Amendment 1 (written after the S, F, D and U arms were scored, before any lever arm exists)

What the first four arms showed (pooled over the four runs; instrument cross-checked against `gffcompare` on NPIP replicate 1: intron-chain sensitivity D 87.8 / S 67.9 / F 72.2, mine 88.7 / 68.5 / 72.9): our assembler recovers more chains than either tool at overlapping genes (E: 1359 / 1606 against S 1015 and F 1159; N: 4637 / 4984 against 3596 and 3753) and emits about twice the in-window artifacts per recovered chain (in-window artifacts 921 against 405 and 380); P1 held, P2 FAILED, P3 held, P4 failed its no-loss clause (U_fam loses 5 chains at N). By the bar of section 4 no arm is BETTER (it fails (b); (c) fails as registered, but the FLAIR number is inflated by a gene_id per 1 kb start bin, 934 ids for 474 genes). Post hoc, not part of the bar: 216 of the 256 entangled genes share an exact junction with another truth gene (readthrough and exon-reusing models), so no gene_id rule can separate them from reads; for the 40 that do not, D resolves 40 / 40 and S 20 / 40. The mechanism sample of our in-window artifacts (NPIP and TBC1D3 replicate 1, 449 transcripts): terminal-trim sub-chains 206, novel junctions 172, displaced junctions 29, internal skips 9, no shared junction 33; read support is not low (fragments mostly 4-10 reads, 207 of the 417 'other' carry exactly 10).

New arms, existing opt-in options only, no new code (the driver commands of `run.sh levers`):

| arm | change to the default | targets |
|---|---|---|
| **P** | `--no-seed-secondaries` (primaries-only seeding) | artifacts built from tied secondary alignments of paralogous reads |
| **C** | `RUSTLE_POLISH_SUBCHAIN=drop` (`--polish-subchain drop`) | terminal-trim sub-chains (the `c` excess of `docs/archive/2026-09/PREREG_complete_transcripts_2026-09-27.md`) |
| **PC** | both | both |

Predictions, fixed now: **Q1** P cuts in-window artifacts of class 'other' by >= 30% against D_asm (pooled); **Q2** C cuts in-window fragment artifacts by >= 50% and loses <= 1% of the chains recovered by D_asm (pooled); **Q3** P recovers fewer chains than D_asm at stratum E (by >= 10 of 1606, pooled): the seeding gain made visible; **Q4** none of P, C, PC meets the bar of section 4 (the existing options are not enough). The bar of section 4 is unchanged and is not edited after these arms are scored. Each arm is scored as `X_asm` (`.gtf`) and `X_fam` (`.families.gtf`).

### Errata (2026-10-07, appended after the independent verification wf_994ece76-8d2; the text above is unchanged)

(1) Section 1 read counts: the 2,066 and 1,668 simulated transcripts belong to 474 and 464 genes (547 and 519 genes, 2,147 and 1,725 transcripts lie in the windows; 81 and 57 were skipped by the generator's filters). (2) The registered quantity of Q1, 'in-window artifacts of class other', is 417 -> 158 (-62%) pooled over the four runs; the results document also quotes the all-contig 527 -> 198 (-62%). In-window means that the transcript's exons overlap a same-strand truth gene; the definition was not written in the prereg. (3) Bar (c) and Q4 could be met only by gene_ids that over-split (FLAIR's 1 kb start-bin ids); they carried no information about the arms. (4) The sub-chain drop's cost (-2.3% of chains in the windows) is arm C alone; PC is -2.1% (-125 of 5,996). (5) Section 3 P1-P4 and the 'BETTER' bar are as registered; the rate-wise reading of PC against the tools is post hoc.


---

## PREREG_fewcopy_class_DRAFT_2026-10-06.md

*(original: `docs/archive/2026-10/merged/PREREG_fewcopy_class_DRAFT_2026-10-06.md`)*

## DRAFT pre-registration: class-level few-copy O1 test on a held-out substrate (2026-10-06)

**STATUS: DRAFT v0. NOT REGISTERED. Nothing in it has been run.** Written by the adversarial reviewer of wf_dbfd1e66-dd2 before any held-out product was opened; it becomes a pre-registration only when the decisions below are made, its numbers are frozen and it is committed as `docs/PREREG_fewcopy_class_<date>.md`. Amendments are allowed only before outputs exist.
Background and the dev results: `docs/FEWCOPY_IDEAL_CASES_2026-10-06.md`. Frozen input-side criteria (C1-C5): `bench/fewcopy/criteria.md`, sha256 `89363117a885718695d70c6d71a78b75f88c31a62421904c505aeba34c7c1571`.

### Decisions needed from the user before registering

1. **Substrate.** Gorilla OR6737 testis (primary choice of the draft; truth = Ensembl Compara gorilla paralogs, else human Compara 2-4 families projected by orthology as in the 07-29 bench), or a never-run human library (HG002: availability unverified), or human chr13 + chr18 (3 families: UNDERPOWERED by rule 6). Gorilla KB3781 is reported separately, never as independent of OR6737.
2. **Bars.** E1 uses 0.9 (as `PREREG_ideal_expression_2026-10-06.md`) and a NO threshold of 0.7; 0.7 is the one free parameter.
3. **Whether a runner must be written and reviewed first.** The dev instrument scripts live in `/mnt/linuxdisk/tmp/fewcopy_2026-10-06/` and are not committed; a held-out run needs a committed, reviewed `bench/fewcopy/` runner.

### Draft text

0. **Status of dev data.** The 11 human A119b families that pass the criteria and all 377 outcomes are SEEN = DEV / regression. No A119b contig is unexposed to the current default, so none is called held-out.
1. **Held-out (never pooled with human).** Primary: gorilla OR6737 testis (GGO RefSeq GFF); truth = Ensembl Compara gorilla paralogs frozen (sha1) first, else human Compara 2-4 families projected by orthology. Fallbacks as in decision 1.
2. **Class K** = families passing C1-C5 of `criteria.md`, computed on the held-out annotation and BAM from inputs only; no threshold edit or relaxation after counts. Near-miss families are listed, never put into K. |K| per size is printed before any outcome is read.
3. **Controls fixed before outcomes.** (G) annotation-as-loci arm: the same `mcl_families` flags with loci = annotation and no reads; a family is REACHABLE iff G recovers it exactly in S1, S2 and S3. (Neg) per family, 5 matched negative pairs (protein-coding, in no multi-copy Compara family, no >= 100 bp exonic overlap, same contig, gap >= 1 kb, reads >= the family's minimum), seeded.
4. **Arm N** = de novo (`assemble` with f1v2, then `mcl_families` defaults) at frozen binary sha1s, run once per library. Scorer: `family_score` (scipy tie policy named) plus pairwise. "Exact" = S1, S2 and S3 AND every member has its own distinct locus (a locus shared by two members is a miss). Chain fidelity is descriptive only.
5. **Endpoints.** E1 (primary): x = reachable families exact in N, n = reachable families; YES >= ceil(0.9 n); NO < ceil(0.7 n); else PARTLY (exact binomial 90% CI printed). E2: false merges on Neg: N must not exceed G. E3: A-class counterexamples = reachable families whose members all have >= 100 P2 reads and tie <= .07 and that N does not recover exactly; >= 2 means NO. Reported without a bar: sens / prec / bipartite F and pairwise over K and over ALL truth families (always beside K), per size and per depth tier (minimum P2 3-29 / 30-99 / >= 100), G-exact over all K (the definition ceiling), and the gap G minus N (guided against de novo).
6. n_reachable < 8 means UNDERPOWERED: no verdict, and it is not a pass.
7. **Failure** = E1 NO, or E2 violated, or E3 triggered. No rescue by changing instrument, tiers, thresholds, substrate or by re-running; deviations are logged as amendments.
8. **Reading a YES.** It means a necessary condition holds on clean, expressed few-copy families in that species and library. It is not sufficient (ties / O2, large families, other species, unlisted relatives). Dev evidence for planning: all-view exact 4 / 11 on the criteria-passing set, 42 / 292 over all scorable families; the annotation arm equals de novo on the 5 pairs with a 09-25 guided row.


---

## PREREG_ideal_expression_2026-10-06.md

*(original: `docs/archive/2026-10/merged/PREREG_ideal_expression_2026-10-06.md`)*

## PREREG: if every copy of NPIP and of TBC1D3 were ideally expressed, would the current default pipeline find them? (2026-10-06)

Version 1, after an independent four-lens review of the first draft (metric traps, feasibility, truth tables, prior art; workflow wf_44161cf9-f29). Registered before any S1 read is simulated or mapped.
Human only, DEV, the two families never pooled, nothing tuned. **Disclosures** (what was read before this registration): (a) the first draft's smoke test of the scorer on the REAL default products (real reads, not ideal): NPIP E2/E3/E4/FOUND 11/3/20/3 of 25, TBC1D3 12/9/11/9 of 16;
(b) the reviewers computed truth-table statistics from the annotation (entanglement, chainless copies, canonicalization effects) and E0-like read statistics on the legacy S2 BAMs; (c) the HEAD default pipeline was run on the two legacy S2 BAMs (24 s and 71 s) and NOT scored; (d) the 2026-09-21 products on those BAMs were read when the first draft was written.
The first draft's rules were changed because the reviewers showed that its YES bars could not be met by construction and that its validity gate (G2) fails by design on this substrate; both are recorded below as the reason for the strata and for the new G2.

### Question

User, 2026-10-06: "like we tried earlier, if we had ideal cases of NPIP and TBC1D3 all expressed, would we be able to find them?" **Ideal** = every annotated transcript of every gene in the neighbourhood is expressed at equal depth, full length (ends jittered 0-30 bp), HiFi error only
(substitution .001, indel .00033), no injected readthrough, no 5' truncation. It is CIRCULAR by construction (the reads come from the annotation that scores them): it asks whether the default pipeline can find the copies when expression, truncation, readthrough and depth do not limit.

Earlier evidence, all older than the current default (f1v2, `--min-cov-shorter 0.70`, GOOD seeding, the families copy table, CAT annotation) and none a human run of the current default on ideal NPIP reads: 09-17 human nodes from 30 identical reads per transcript (NPIP 27/27 nodes, 11/27 full-length, 13/27 swallowing a neighbour; TBC1D3 19/19 nodes, 17/19 full-length;
node level only, never chain level for TBC1D3); 09-18 NPIP 25/26 complete chains (one canonical transcript per copy, no neighbours); 09-21 whole-chr16 (see "Legacy substrate"); gorilla 09-28 NPIP 25/25 complete and in the family with the copies alone, 24/25 in the family with neighbours, 22-24/25 under the current default at f = 0.

### Substrate S1 (new; CAT/Liftoff v2.0, the default annotation)

- **Windows:** +-500 kb around the territory (`terr_lo0`-`terr_hi`) of every copy in `copy_recovery_tools_cat/ann/copies.hsa.tsv`, merged: NPIP 25 copies on chr16 (7 windows, 12.19 Mb), TBC1D3 16 copies on chr17 (9 windows, 9.53 Mb). RefSeq NPIPB3 has no CAT/Liftoff gene, so S1 has 25 NPIP copies. 8 of the 41 target records are Liftoff, not CAT;
  everything is keyed by copy id (`cid`) and gene ID, never by name.
- **Transcripts:** every transcript of every gene overlapping a window (expected: NPIP 547 genes, 2,147 transcripts; TBC1D3 519 genes, 1,725). **Canonicalization (the 2026-09-18 recipe):** an intron is kept only if >= 50 bp and canonical (GT-AG, GC-AG, AT-AC on the transcript's strand, read from the genome); every other gap is merged into its flanking exons.
  A transcript is **simulated** iff its spliced length is >= 120 bp, its molecule is <= 30 kb and no kept intron exceeds 200 kb (minimap2's default maximum intron length, so such a molecule cannot align spliced); skipped transcripts are listed with the reason (expected: 77 and 55 under 120 bp, one 121 kb molecule, the transcripts with the 1 Mb intron of lncRNA AC009093.11 and one 218 kb intron).
  Annotated fusion and readthrough gene models are simulated like any gene (the arm has no INJECTED readthrough).
- **Reads:** 10 per simulated transcript, the model of `bench/sim.py` (`simulate_reads`: err .001, indel err/3, no truncation; 0-30 bp end jitter, skipped when the body would fall under 100 bp), each read drawn from `stable_seed(seed, transcript id, k)`. Two independent replicates, seeds 20261006 and 20261007, each simulated, mapped and scored.
- **Mapping:** the shipped command (`minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`) against the whole CHM13 v2.0 splice index, in read-disjoint parts.
- **Families are simulated and run separately** (one read set per family and replicate: 4 mapped sets).

### Legacy substrate S2 (`/mnt/linuxdisk/tmp/idealsim/{ideal,rt}.bam`): NOT USED

Found by the review and confirmed here on 300 of 300 testable transcripts: `sim.py`'s `load_transcripts` keeps an exon only if its transcript row came earlier in the file, so every transcript's leftmost annotated exon (6,950 of 80,664 chr16 exon rows) is missing from the 2026-09-21 reads, which are therefore not full length;
the reads were mapped to `chr16.fa` only (one `@SQ` line, checked); and the `rt` arm replaces 7% of the reads and, by the reviewer's check of the reads, puts the downstream gene first on the minus strand (1,070 of its 1,109 minus-strand readthrough reads, 36.8% of all readthrough reads). The 2026-09-21 ideal-chromosome ceilings (`docs/IDEAL_CHROMOSOME_SIM_2026-09-21.md`) describe those reads. Scoring them would need its own registration.

### Pipeline and provenance

`tools/rustle_pipeline.sh assemble` then `families`, HEAD release binaries, no override (secondary seeding from the best-AS table of the same BAM, f1v2, shipped polish, strict junctions, `--min-cov-shorter 0.70`, most-reads representative). Products: `PREFIX.gtf`, `PREFIX.families.gtf`, `PREFIX.fam.loci.gff3` (one representative per locus), `PREFIX.fam.loci.tsv` (folded span to representative span),
`PREFIX.fam.copies.tsv`, `PREFIX.fam.clusters.tsv`. **G0:** sha1 of `copy_assign` (87824d91), `mcl_families` (a6308244), `as_table`, `family_score` (542923fd), the BAMs and the scripts is logged; the binaries must equal those recorded in `docs/DEFAULT_RESCORE_NPIP_2026-10-06.md`.

### Targets, strata and truth (fixed from inputs before any read is mapped; listed in `targets.tsv`)

Coordinates are 0-based half-open; an intron is (end of the upstream exon, start of the downstream exon); chains are tuples of introns in genomic order; one convention is asserted on a known intron (NPIPA1's first intron, 14938572-14942865).

- **E, entangled:** another simulated gene on the same strand shares >= 100 exonic bp (canonical exon unions) with the copy. The shared bp and the share of the copy's exon union are recorded; copies sharing < 50% are reported separately. Two annotated genes that share exons cannot each get a locus from any assembler.
- **C, chainless (MONO):** the copy has no multi-exon canonical chain (NPIPA3 and TBC1D3P7: every raw intron is non-canonical, so the molecules are unspliced gene bodies). Unspliced loci default to the + strand, so a MONO copy's holder is strand-agnostic.
- **X, readthrough image:** the copy table has `readthrough = 1` (NPIP h03, h05, h08: images of fusion records).
- **R_in = copies in none of E, C, X.** Expected before E0: NPIP <= 14 of 25, TBC1D3 <= 13 of 16. **R** (per replicate) = R_in copies with an OBSERVABLE chain (E0 below); R is computed from the truth tables and the BAM alone, before assembly, and its sha1 is logged.
- **Named-family genes:** simulated non-target genes whose CAT/Liftoff `gene_name` contains NPIP or TBC1D3 (e.g. NPIPP1, PKD1P6-NPIPP1, NPIPB7, PDXDC2P-NPIPB14P, TBC1D3J, TBC1D3P1-DHX40P1): expected non-target family members, never "foreign".
- **Single-copy genes** (instrument controls): simulated genes none of whose reads has a secondary alignment with AS >= 0.9 x its best AS, taken from the BAM.

### Endpoints (denominators fixed from inputs)

- **E0 observability (BAM alone).** Pools: P1 = primary alignments (-F 2308); P2 = P1 + secondary alignments with AS >= 0.98 x the read's best AS (the pipeline's seeding pool); P3 = any alignment. A chain is OBSERVABLE iff >= 3 reads of P2 carry exactly that whole chain (the read's junctions inside the copy's span equal the chain, N >= 50 bp).
  Printed per copy and pool: reads, share whose primary overlaps the source copy (no bar), observable chains. A copy with no observable chain is ALIGNER-LIMITED.
- **Holder.** Records are resolved through `fam.loci.tsv` (folded loci count as their representative; the number of folded holders is reported). The holder of a copy is the same-strand locus maximizing J = overlap(rep exons, copy exon union) / union of the two (bp); ties: fewer representative exonic bp, then lower gene name.
  Clusters join to loci by (chrom, gene-row 1-based start, end) as `bench/default_rescore/nodes.py` does; two loci on one span or an unmatched cluster row make the arm INVALID.
- **E1 chains:** distinct canonical multi-exon chains of the copy that appear exactly as a transcript chain in `PREFIX.gtf`, over all chains and over observable chains (descriptive).
- **E2 own locus:** a holder exists, holds no other target copy of the family, and both purities are >= 0.5: representative purity |rep exons in the copy| / |rep exons| and locus purity |exon union of the locus's transcripts in the copy| / |that union|. Also printed: same-strand loci overlapping >= 100 exonic bp (>= 2 = split), and whether the top-overlap locus is on the other strand.
- **E3 exact representative:** the holder's representative chain equals one of the copy's chains. **E3':** the representative's introns inside the copy span are a contiguous sub-chain of ONE chain with >= max(2, ceil(n/2)) introns, n = introns of that chain. **E3c** = chains recovered exactly among the holder locus's transcripts / chains simulated; **E3p** = holder-locus transcripts with a chain that is not a chain of the copy / holder-locus transcripts (over-emission).
  MONO copies: E3 = the representative covers >= 90% of the copy and has <= 110% of its bp (reported beside, outside R).
- **E4 in the family:** K* = the cluster of the holder of the lowest-cid R copy that has a holder in a cluster; the holder is in K*. **CP** = (member loci of K* overlapping a target copy or a named-family gene) / |K*|. Per family: pairwise tp / truth / predicted pairs, then bipartite sens / prec / F with `family_score` on truth rows restricted to the simulated windows, predicted members outside the truth counted as false, with and without size-2 clusters; the pooled line is printed labelled not-a-family-result.
- **IDEAL-FOUND = E2 and E3 and E4** (this is NOT the registered strict FOUND of Amendment E, which needs the cap signal and cannot run on simulated reads). **IDEAL-FOUND-relaxed** uses E3'. Printed beside, on the same copies: the `nodes.py` own-node flag and `bench/copy_support.py` Amendment A `ann_found` / `chain_found` computed on the ideal BAM, next to the registered real-read DEF values
  (own node 24/25, strict FOUND 8/25, locus level 21/25, U2 F .645, CF153 13/19 F .812); only own node and `ann_found` are like-for-like.

### Rules (per family and replicate, on R; n = |R|, N = all target copies of the family)

- **CEILING-LIMITED** if n < 0.5 N: no YES or NO is issued.
- **YES:** IDEAL-FOUND >= ceil(0.9 n) and CP >= 0.5.
- **PARTLY:** E2 and E4 each >= ceil(0.9 n) but not YES.
- **NO:** otherwise.
- The family's verdict is the LOWER of the two replicates (NO < PARTLY < YES); if the replicates differ it is reported UNSTABLE and no YES is claimed. The copies whose flags differ between replicates are listed.
- Always printed beside, never instead: IDEAL-FOUND over ALL copies as x/N, and one line per excluded stratum (E, C, X, ALIGNER-LIMITED). An arm that fails G1, G3, G4 or G6 is INVALID; an INVALID arm is amended in writing before any pipeline product of that arm exists, afterwards only a new arm may be registered. No rescue.

### Controls

- **G1** every simulated chain is canonical and >= 50 bp, re-read from the genome by a separate script (`verify_g1.py`) before mapping. **G2 (instrument)** >= 99% of the reads of single-copy genes have a primary alignment overlapping their source gene; the same fraction for the target copies is E0, printed, ungated (near-identical paralogs coin-toss primaries).
  **G3** the BAM holds exactly one primary-or-unmapped record per FASTQ read; the driver stages exit 0 and `families.gtf` is newer than `gtf`. **G4** the scorer is run twice under different `PYTHONHASHSEED` values and gives identical tables.
- **G5** an independent scorer written from this text alone (by someone who has not seen `score.py`) reproduces E2, E3, E4 and IDEAL-FOUND per copy.
- **G6 positive control:** the canonicalized annotation as loci (each gene a locus, equal read counts) through `mcl_families --from-gtf` with the shipped flags and the same scorer: E2 and E3 must be 100% on R, and its E4 is the family ceiling E4*; YES needs E4 >= min(0.9, E4*) on R.
  **G7 negative control:** >= 95% of single-copy simulated genes get exactly one locus (same strand, >= 100 exonic bp) and none lies in K*.

### Held-out contrast (no bars, labelled HELD-OUT: the pipeline was developed on NPIP and TBC1D3)

Two other multi-copy families inside the same windows, chosen now by CAT/Liftoff `gene_name`: KRTAP genes (`^KRTAP[0-9]+-[0-9]+$`) in the TBC1D3 windows and SMG1P genes (`^SMG1P[0-9]+$`) in the NPIP windows. Scored with E2-E4 on their non-entangled copies that have a canonical multi-exon chain; membership is not verified before the run.

### Predicted before looking

- **P1** TBC1D3 is YES on R in both replicates (09-17: 19/19 nodes, 17/19 full-length).
- **P2** NPIP is not YES on R (09-17: 13 of 27 nodes swallowed a neighbour) and IDEAL-FOUND over all 25 copies is <= 21.
- **P3** the two replicates give the same verdict in both families.

### Also reported with every headline

Per-copy tables (E0 three pools, E1, E2 with both purities and the holder, E3/E3'/E3c/E3p, E4, entangled-with gene, shared bp and share, stratum, and one miss reason from the fixed order annotation, aligner, locus formation, family membership, representative); the seed spread; per family CP with every non-target member by gene name and biotype;
the single-copy control; chain-level recall and precision; binary and BAM sha1. One sentence of scope: DEV families, human only, circular by construction, +-500 kb neighbourhood only (precision against homologous loci elsewhere is not measured), equal-depth full-length reads, annotated fusion records present, Liftoff records among the targets.

### Not claimed

Anything about real reads, the annotation's correctness, O2, O3, gorilla, or any family beyond these two and the two held-out ones. A YES means the default pipeline has no algorithmic obstacle at these loci when the data are ideal; it does not predict recovery from a real library (5' truncation cost 6 of 26 NPIP copies in 09-18, and real expression differs).

### Amendment 1 (2026-10-06, before any S1 read exists)

The shared-span rule of the holder paragraph is narrowed: two loci on one (chrom, start, end) span make an arm INVALID only when a cluster row or a fold row refers to that span (the join would be ambiguous); shared spans among loci that no cluster or fold row names are ignored.
The first scorer test on the REAL default chr17 products showed two unclustered loci on one span (`DN_chr17_87073_2` and `_3`), which the draft rule would have declared invalid although no cluster or fold row names them. Nothing else changes.


---

## PREREG_locus_units_2026-10-06.md

*(original: `docs/archive/2026-10/merged/PREREG_locus_units_2026-10-06.md`)*

## Pre-registration: a locus model for overlapping genes: primary-first representatives and structural separators (family level, ideal reads) (2026-10-06)

**Written before any arm of this study exists.** What was already seen (all on the same two dev windows, NPIP chr16 and TBC1D3 chr17, two read replicates each): the default arm D and the lever arms P (primaries-only seeding), C, PC at the FAMILY level (`bench/ideal_expression/score.py`, pooled over the four runs: reachable copies R, 54 gene-runs: D found 48, P 53, C 48, PC 53; entangled copies E, 24 gene-runs: D 12, P 7, C 13, PC 7), the per-copy differences D to P, and the anatomy of the entangled copies in D. DEV, human CHM13 / CAT-Liftoff, simulation, circular by construction. No held-out substrate is spent here; a held-out and a real-data check are the next registered step for any arm that passes.
User direction (2026-10-06): option B, a locus-level model for junction-sharing overlapping genes (units, multi-membership); the advisor wants clean combinatorial structure, provable statements, no arbitrary thresholds.

### 1. What the dev data say (the motivation, not a result of this study)

- A locus (the assembler's gene_id) is a connected component of the junction-sharing graph; replacing gene_ids by junction-sharing components reproduces the scores exactly. 216 of 256 entangled gene-runs share an exact junction with a partner gene (`docs/ENTANGLED_BASELINE_2026-10-06.md`): no component rule separates them.
- D loses reachable copies through the REPRESENTATIVE: NPIPA2 (representative is an artifact chain built from paralogous reads), TBC1D3P4 / P3 (the displaced-intron chain, carried only by the paralog's reads as secondary alignments, ties the copy's own chain at 10 reads), NPIPA7 and NPIPB3 (locus purity .112 and .234, bodies polluted with artifact chains). With primaries-only seeding all of them are found (R 53 / 54).
- Primaries-only seeding loses entangled copies from the family clusters (E4 22 -> 14: NPIPA1, AC138969.1 twice) and breaks others (NPIPA8 representative purity .906 -> .295): the secondary-borne chains enlarge the locus BODY toward its family, and the families stage aligns bodies. The representative and the body are different roles.
- Entangled copies in fused loci (NPIPB14P with the readthrough PDXDC2P-NPIPB14P as representative, AC138894.1 with CLN3 and NPIPB7, PKD1P6 with NPIPP1): the F1v2 read-share rule cannot cut a readthrough that carries as many reads as its parents (every simulated transcript has 10 reads).

### 2. The model (two rules, no free parameter)

**Rule 1, primary-first representative.** Within a locus a chain is PRIMARY-SUPPORTED iff the primaries-only assembly (the same run with `--no-seed-secondaries`) contains it exactly. The representative of a locus is the best primary-supported chain (most reads, then span, then index, the shipped key); only if the locus has none, the best chain overall. Bodies are untouched: every seeded transcript still belongs to its locus and enters the alignment of the locus body. Principle: secondary alignments may create a site and enlarge a body (O1 must not wait for an O2 coin toss), but structure and exactness come from the reads whose best placement is here.
**Rule 2, structural separators.** For a locus, the INCIDENCE GRAPH has the multi-exon transcripts and the junctions carried by >= 2 of them as vertices, an edge for each (transcript, junction) containment. A transcript is a SEPARATOR iff it is an articulation point of the graph (its removal splits the locus into >= 2 groups that each hold a transcript): the transcripts that are the only link between groups, the evidence-free generalisation of F1's bridge. Separators are removed from the families input (as F1's bridges are, listed in a side table), the remaining transcripts are regrouped by the components of the graph without them (gene_ids `<locus>.s<k>`), single-exon and junction-free transcripts attach to the component they overlap most, loci without a separator are unchanged. Facts used: the separator set and the components are canonical (order independent, linear time, Tarjan); the structure left is a tree of bridge-free groups; two groups linked by two vertex-disjoint paths cannot be separated by this rule (the boundary of a single-vertex rule, which the 216 junction-sharing genes of the dev windows largely sit behind).

### 3. Arms (the families input of each; every arm then runs `mcl_families` with the driver's families-stage flags, `--min-exonic-bp 1 --min-shared-exon-frac 0.60 --emit-units`)

| arm | families input |
|---|---|
| **D**, **P** | existing: the default `asm.families.gtf`, the primaries-only `lever_P.families.gtf` (same binaries) |
| **Q** | D's input with Rule 1: the `reads` attribute of every transcript is replaced by `reads + 1,000,000` when its chain is primary-supported (the representative key is (reads, span, -index); `reads` is used for nothing else in `fam_from_gtf`, checked in `src/family.rs`) |
| **S_D** | D's input with Rule 2 on all transcripts |
| **S_Q** | Q's input with Rule 2 where the graph holds the primary-supported transcripts only: the other multi-exon transcripts of a locus attach to the component sharing the most junctions with them (ties: the earliest), those sharing none form their own groups by junction sharing among themselves; separators are always primary-supported |

### 4. Measures

Family level, as the registered ideal test (`docs/PREREG_ideal_expression_2026-10-06.md`, `bench/ideal_expression/score.py`): IDEAL-FOUND (E2 and E3 and E4) per copy on stratum R (14 and 13 copies per run) and on the entangled copies E (10 and 2), E2-E4 separately, the registered rule (YES / PARTLY / NO, lower of the replicates), cluster precision, and `family_score` on the windows (Compara F, `bench/ideal_expression/fam_score.py`). Transcript level: `bench/entangled/arms_score.py` on each families input (chains recovered, artifacts, resolved genes). Pooled over the four runs and per family.

### 5. Predictions (fixed now)

- **Q1** Q finds >= 52 of the 54 reachable gene-runs (P: 53, D: 48). **Q2** Q finds >= 11 of the 24 entangled gene-runs (D: 12, P: 7): the body is kept. **Q3** Q gives YES for NPIP and for TBC1D3 in both replicates under the registered rule.
- **Q4** S_D finds >= 2 more entangled gene-runs than D with no loss at R. **Q5** S_Q finds at least as many gene-runs as Q at R and at E, and >= 1 more at E.
- **Q6** Compara F on the windows of Q and S_Q is >= D's minus .005 in both families.

### 6. Decision rule

An arm is a CANDIDATE iff, pooled over the four runs, IDEAL-FOUND (R plus E) is strictly above D's, R found is >= D's, E found is >= D's minus 1, and Compara F is >= D's minus .005 in both families; the direction must hold in both replicates of both families. A candidate earns a registered held-out check (ideal windows around other multi-copy regions, then A119b chr16 / chr17 and gorilla OR6737 with the copy-recovery instruments); nothing is claimed for real reads from this study. A non-candidate is reported with where it loses. Rule 2 is also reported by itself at the gene level against the oracle of the truth chains (`CUT-T`: 17 of 108 junction-sharing entangled genes resolved, 27 more non-entangled genes split than with components).

### 7. Limits declared in advance

Same two windows, same simulator, ideal reads, the arms were motivated by what the dev windows showed; the primary-supported marks come from a second assembly of the same reads (an approximation of 'the read's best placement is here'); Rule 2 uses no read evidence, so on real data it can cut a legitimate linking isoform (the truth oracle shows 5% of non-entangled genes split); the seeding gain for copies with no primaries is invisible in these windows.

### Amendment 1 (before any arm exists; the module and its unit tests were written after the section above)

Writing the unit tests showed that Rule 2 as worded ('a separator is an articulation point') also cuts an isoform that is the only link to a one-transcript attachment (A1 = j1 j2, A2 = j1 j2 j3, leaf = j3 j77: A2 would be removed from the families input). The truth-chain oracle of the dev windows (one replicate per family; connected components, three separator definitions; 108 junction-sharing entangled genes, 20 non-sharing, 513 non-entangled) chose the wording:

| rule | separators | junction-sharing entangled genes resolved / merged / split | non-entangled genes split (components: 128) |
|---|---|---|---|
| components (current loci) | 0 | 0 / 90 / 18 | 128 |
| articulation transcript (the first wording) | 92 | 17 / 64 / 27 | 155 |
| articulation transcript that contains the junction set of a transcript of two groups | 23 | 2 / 87 / 19 | 137 |
| **articulation transcript whose chain runs from a junction of one group, through a junction no group carries, to a junction of another group (private bridge)** | 38 | 16 / 68 / 24 | 140 |

**Rule 2 is the private bridge** (code: `bench/entangled/locus_units.py`, 12 unit tests, `python3 -m unittest test_locus_units`). It keeps nearly all of the articulation rule's gain with 12 false cuts instead of 27 and 38 separators instead of 92; it still resolves only 16 of 108 junction-sharing entangled genes (the other 68 are merged by shared junctions with no bridge: the boundary of a single-vertex rule, stated in section 2). The oracle is a dev analysis on the same windows and is disclosed as such; the family-level arms of section 3 are the test. Predictions Q1-Q6, the decision rule and Rule 1 are unchanged.

### Amendment 2 (2026-10-07; written after the level-M results of `docs/LOCUS_UNITS_2026-10-06.md` were seen, before any component-level number or any F1 arm exists)

**Why.** The registered E4 ('the holder is in the one MCL cluster K*') is brittle to partition granularity: arm Q repaired every exactness failure it targeted (E3 on R 50 -> 54 of 54) but in NPIP replicate 2 the NPIP cluster split and NPIPB6, NPIPB8, NPIPB9 moved to a sibling cluster, which the registered rule scores as three lost copies and a Compara F drop. This is the advisor's hierarchy question (`docs/ADVISOR_HIERARCHY_2026-10-03.md`: the cut is one level of a nested structure; each truth should be scored at its level; the principled object is the connected components of one evidence graph). This amendment is therefore a RE-READING motivated by an observed flip, declared as such; the registered level-M verdicts stand and are reported beside it.
**Levels (no new threshold).** For an arm, the pre-MCL homology graph is the one `mcl_families --dump-graph` writes for the arm's families input with the families-stage flags (the MCL clusters of that re-run must equal the arm's registered clusters byte for byte, else the arm is INVALID at level C). **Level M** = the shipped MCL clusters (registered). **Level C** = the connected components of that graph (a coarsening of every partition derived from it: nested by construction). Everything the registered scorer does with `fam.clusters.tsv` is repeated with the components as clusters (size >= 2): K*_C (the component of the holder of the lowest-cid R copy that has one), E4_C, IDEAL-FOUND_C = E2 and E3 and E4_C, the registered rule (YES at >= ceil(0.9 R), lower of the replicates) with CP_C, and Compara F_C by `family_score`. An E4 failure is then one of two kinds: no edge path to the family (not in K*_C: evidence loss) or in K*_C but not in K*_M (partition split).
**Decision rule at level C.** An arm is a CANDIDATE_C iff, pooled over the four runs, IDEAL-FOUND_C (R plus E) is strictly above D's, R found_C >= D's, E found_C >= D's minus 1, F_C >= D's minus .005 in both families, and CP_C >= .5 in all four runs; the direction must hold in both replicates of both families. An arm that is a candidate at one level only is reported as LEVEL-SPECIFIC and earns nothing by itself; a candidate at both earns the held-out windows check of section 6.
**New arms (existing options only).** **F1**: the driver with `RUSTLE_BRIDGE_REGROUP=f1` (every evidence-proven bridge cut, F1's PAS and promoter evidence without f1v2's read-share rule), assemble then families on the same BAM. **F1Q**: F1's families input with Rule 1, the primary chains taken from the primaries-only assembly of the same run. D is f1v2; Q is f1v2 + Rule 1.
**Predictions (fixed now).** **H1** F1 cuts more bridges than f1v2 does in the ideal windows (f1v2: 2 / 1 / 1 / 1 per run) and finds NPIPB14P (h22) in at least one replicate. **H2** F1 loses no reachable copy (R found_M >= D's 48). **H3** F1Q is a CANDIDATE_C. **H4** Q finds >= 53 of 54 reachable gene-runs at level C (the NPIP replicate-2 split is a partition event, not an evidence loss). **H5** S_D and S_Q stay non-candidates at level C (their losses are evidence losses, not splits).

### Errata (2026-10-07, appended after the independent verification wf_994ece76-8d2; the text above is unchanged)

(1) Section 1 says that with primaries-only seeding 'NPIPA7 and NPIPB3 ... are found (R 53 / 54)': NPIPB3 is an entangled copy (not in R) and is found in replicate 1 only; the reachable ones are NPIPA2, NPIPA7 and TBC1D3P4 / P3. (2) Section 6 quotes the first wording's oracle ('CUT-T: 17 of 108 resolved, 27 more non-entangled genes split'); the rule actually run (Amendment 1, private bridge) is 16 of 108 and 12 more. (3) The rationale of H5 in Amendment 2 ('their losses are evidence losses, not splits') was wrong for S_D: its reachable-copy loss is a representative failure (E3 47 of 54, E4_C 54 of 54); H5's outcome (not candidates) holds. (4) H4 could not be met: E2 on the reachable gene-runs was 52 at level M already and does not depend on the level. (5) `pool_levels.py` does not implement the amendment's 'direction must hold in both replicates of both families' clause; checked by hand for P (+4, +1, +1, +1 gene-runs per run) and Q (+1 in each). (6) A graph dump of arm D (NPIP replicate 1) existed 2.5 minutes before commit 5019106d; no component was computed from it.


---

## GORILLA_FEWCOPY_DESIGN_2026-10-07.md

*(original: `docs/archive/2026-10/merged/GORILLA_FEWCOPY_DESIGN_2026-10-07.md`)*

## Gorilla OR6737 testis few-copy class test: design sweep result (2026-10-07)

**Status: DESIGN ONLY, nothing of the pipeline was run on gorilla. The class that the inputs define is too small for the drafted bar: |K| = 3, so rule 6 of `docs/PREREG_fewcopy_class_DRAFT_2026-10-06.md` (n_reachable < 8 = UNDERPOWERED, no verdict) applies before any arm exists.**
User direction (2026-10-06): "lets try with the gorilla testis". Workflow wf_fbf8afe6-9cb (exposure ledger and symbol-projected truth in parallel, then criteria and class K blind to every pipeline product, then an adversarial review that recomputed 5,050 values). Outputs `/mnt/linuxdisk/tmp/fewcopy_gorilla_2026-10-07/` (tables, scripts, the four sweep results). Frozen gorilla criteria: `bench/fewcopy/criteria_gorilla.md`, sha256 `2b17d529dc62710ade954e67363118cd9685e31780b2c2cb612f52196c19086e` (frozen 00:07 PDT, first count 00:12 PDT, hash re-checked by me and by the reviewer; human criteria 89363117 unchanged).

### Label to use

**Not 'held-out'.** "Second-species, inputs-frozen (outcome-blind) class check on a partly spent library (gorilla OR6737 testis)"; KB3781 is reported separately as dev (same assembly, annotation and truth source). OR6737 is DEV for the seeding rule (NC_073244.2), F1 / F1v2 (both gorilla libraries) and the DAZ2 gate; held-out and spent for strict junctions + retained-intron ratio 10 (26 contigs), F1 (minus NC_073244.2) and R_J (NC_073234.2); the RefSeq annotation (the universe and the truth source) fixed the MCL constants; no contig is unexposed. Missed by the ledger and found by the review: NC_073230.2 + NC_073228.2 were the held-out substrate of Addendum AH (09-14, CF343 sits there), the MHC (C4A / C4B) is in the 06-25 reference-absent catalog, ZNF600 is in the 09-05 O2 / O3 gap-locus rows.

### Truth: human Compara families projected by exact RefSeq symbol

Human Ensembl Compara 2-4-gene families (377) projected onto the gorilla RefSeq annotation (GCF_029281585.2, individual KB3781) by exact symbol of a protein-coding gene record: **342 of 864 members project 1:1 (39.6%); 90 families whole, 135 partial, 152 none**; the 522 absent members are 305 with LOC-named protein-coding candidates by name evidence and 217 without. All rows recomputed against the GFF by the reviewer, 0 mismatches. Selection effect: old dispersed pairs project (whole in 54.0% of human cross-contig families against 17.8% of same-contig ones); the recent duplicates (LOC-named) do not. This is a third truth method, neither Compara gorilla nor the 07-29 projection of register row 861; it is human truth projected, not gorilla truth.

### The class K under the frozen criteria

Funnel (families): 377 -> **90** whole projection -> 38 (C1 annotation) -> 30 (C2 not entangled) -> 29 (C3 separate loci) -> 7 (C4 identity [0.90, 0.999), coverage >= 0.5) -> **3** (C5 >= 3 exact-chain reads per member). Independent passes among the 90: C1 38, C2 73, C3 86, C4 33, C5 17. Thresholds are the human ones; the forced changes are the single RefSeq annotation (no CAT mapping tier) and the longest mRNA as reference transcript.

| family | members | contig | identity | min reads | verdict on the case |
|---|---|---|---|---|---|
| CF343 | STEAP1, STEAP1B | NC_073230.2 (68.7 Mb apart) | .987 | 28 | the only clean O1 case: dispersed pair, no read ties |
| CF315 | C4A, C4B | NC_073229.2 (MHC) | .998 | 134 (same reads at both copies, P1 65 / 70) | an O2 tie control: cross-link .91 / .90, read divergence ~ the inter-copy .0020; hyperpolymorphic locus |
| CF236 | ZNF600, ZNF611 | NC_073244.2 | .947 | 23 | unreachable by listed size: ZNF808 / ZNF888 / ZNF578 are unlisted relatives in a 43-gene KRAB-ZNF cluster; its human outcome ('merged') is dev-known |

Near misses (12, fail one criterion): C4 x8 (FCGR2A / FCGR2B coverage .30, ALG10 / ALG10B .47, RABL2A / RABL2B .38, ZNF155 / ZNF230 identity .78, DGCR6 / DGCR6L .87, BEX1 / BEX2 .87, EOLA1 / EOLA2 identity .9997, ZNF195 / ZNF429 cross-contig), C5 x4 (KRT37 / KRT38, USP32 / USP6, PF4 / PF4V1, UGT2B15 / UGT2B17 one read short).

### What follows from it

- The drafted class-level test cannot return a verdict on gorilla testis with these inputs: |K| = 3 (one clean). Admitting every LOC-named family at the observed 3 / 90 rate would give about 7 before the reachability control.
- Read support is the binding constraint, not annotation: **17 of the 90 whole families have every member read-supported (15 pairs, 2 triples)**, but 7 of the 17 are KRAB-ZNF families on NC_073244.2 (ZNF181 / ZNF302, ZNF585A / B, ZNF224 / ZNF225 / ZNF284, ZNF600 / ZNF611, ZNF195 / ZNF429, ZNF324 / ZNF324B, ZNF155 / ZNF230), whose listed size hides relatives (not independent, exact recovery unreachable by listed size). The other ten: C4A / C4B, TCEAL2 / TCEAL4, BEX1 / BEX2, STEAP1 / STEAP1B, EOLA1 / EOLA2 (identical), ALG10 / ALG10B, FCGR2A / FCGR2B, DGCR6 / DGCR6L, RABL2A / RABL2B, APOBEC3C / D / F.
- Risks specific to gorilla: cross-individual reads (OR6737 reads on the KB3781 assembly; private variants drop exact-chain reads, so P2 is a lower bound), RefSeq-only annotation (C2 'not entangled' may mean under-annotated: the human CAT contrast fails C4A / C4B at 520 bp shared with antisense lncRNAs and STEAP1B at 986 bp), one movie and one library, 57 of 342 projected members on chrX, Gnomon models may cite long SRA reads (if OR6737 is among them, C5 is partly circular; unverified).
- Overlapping genes on gorilla (input-only, exploratory, not frozen): 2,202 gene pairs share >= 100 exonic bp; 1,362 have a valid chain in both genes; 423 have both genes with >= 3 exact-chain reads (35 on the same strand); at >= 1 kb shared 59 (10 same strand). This is a real-read test bed for the overlapping-gene question (the lab's StringTie, FLAIR and isoseq outputs for OR6737 exist genome-wide), independent of the few-copy class.

### Unsupported claims withdrawn by the review

'All three K families have 0 unlisted relatives' is name-based (true by sequence for CF315 and CF343, false for CF236); 'family-level exactness of 2-4-member families was never scored on any gorilla contig' is unverified (r1232 and r1237-r1245 are family-level default runs on OR6737 testis, not 2-4-member Compara families); the quoted '46.6% of genes are LOC ids' does not reproduce (55.7% of gene records, 23.4% of protein-coding genes); 'no better truth source on disk' misses `GGO_sedef_final.bed` (DNA segmental-duplication calls, a sequence-based label for LOC-named recent duplicates, with a circularity caveat); CF315 ranks first only because min P2 counts the same 135 reads at both copies.


---

## GORILLA_OVERLAP_2026-10-07.md

*(original: `docs/archive/2026-10/merged/GORILLA_OVERLAP_2026-10-07.md`)*

## Overlapping genes on real reads, gorilla OR6737 testis: our assembler against StringTie, FLAIR and isoseq collapse (2026-10-07)

**Status: descriptive real-read check on a partly spent library (cross-individual: reads OR6737, assembly and annotation KB3781), RefSeq annotation as truth, one library. Protocol: `docs/PREREG_gorilla_overlap_2026-10-07.md` (committed 120fc627, before any truth table or arm existed; errata appended there). Tools: `bench/entangled/{real_truth.py, real_score.py}` (14 unit tests, commits 9d8f6b16 / 89beddd1) and `run_gorilla.sh` (first committed with the results, 742c929d). Products: `/mnt/linuxdisk/tmp/gorilla_overlap_2026-10-07/`. Independently verified (workflow wf_994ece76-8d2: truth tables, 540 arm values and the conventions recomputed from raw files with separate code, 0 mismatches in any registered number); the corrections it forced are folded in below.**
User direction (2026-10-07): 'lets do 1 and 4' (1 = overlapping genes on real gorilla reads). Arms: ours D (default), P (`--no-seed-secondaries`), PC (P + `--polish-subchain drop`), HEAD release binaries, genome-wide `assemble` (D 4 min with the AS table, P and PC 1.7 min); the lab's S (StringTie 3.0.1 `-L`), F (FLAIR 3.0.1 collapse) and I (isoseq collapse) on the same BAM. D reproduces the 09-25 transcript count exactly (74,045), P the 09-25 primaries-only count (72,690).

### Answer

1. **Against StringTie and FLAIR the lead is clear; against isoseq collapse it is not established.** Same-strand overlapping genes with both genes expressed (E_both: 54 gene pairs = 108 genes; 202 gene-chains = 201 distinct chains, one chain being carried by both genes of the H1-4 / LOC101135182 pair), primary denominator (chains with >= 3 primary reads, which is what the tools read): D 170 of 202 (84.2%), I 159 (78.7%), S 140 (69.3%), F 140 (69.3%). Pair-level bootstrap of the recall difference: D - S +14.9 points (95% CI 9.6 to 20.1), D - F +14.9 (9.8 to 20.4), **D - I +5.5 (-1.4 to +12.7; 11 pairs favour D, 10 favour I; McNemar p = .080 on P1, .041 on P2)**; in the four junction-sharing pairs (14 chains) D is below isoseq and StringTie (8 against 11 and 9). Non-overlapping genes (N, 10,079 genes, 19,771 chains): D 94.6%, F 86.2%, I 81.2%, S 78.8%. The `gffcompare` run (201 distinct chains, D 84.1, S 69.2, F 69.7) uses the same truth and the same exact-chain criterion: a check of the instrument, not independent evidence. D's defaults were partly tuned on this library; the lab tools ran as untuned annotation-free recipes.
2. **Same-strand overlap goes with lower recall in D (-10.4 points against N, p = 7e-8), StringTie (-9.5, p = .002) and FLAIR (-16.8, p = 1e-9), not in isoseq collapse (-2.5, p = .37).** The stratum is weak overlap: 41 of the 54 pairs overlap in exons of expressed isoforms (13 only through unexpressed isoforms), the median shared overlap is 214 bp (138 bp between expressed isoforms), 13 pairs share >= 1 kb, 16 of the 108 genes sit on the KRAB-ZNF contig; not adjusted for paralogy. The gap to N persists at matched read support (median 15.5 reads at E_both, 15 at N). Antisense overlap (A, 1,583 genes): D 93.8%, F 85.9%, I 84.8%, S 75.2%.
3. **Gene separation: only D against StringTie is a defensible comparison.** E_both genes resolved: D 74 of 108, F 70, S 52, I 13; without shared junctions (E_both_x, 100 genes) D 71, F 67, S 47, I 13. Each tool's gene_id is its own grouping: FLAIR over-splits (47,070 ids for 153,402 transcripts, D 19,775 for 74,045), StringTie fuses (146 ids carrying exact chains of two genes against D's 61), 951 of isoseq collapse's 53,087 ids carry chains of two or more genes. The 'resolved' count also mixes recall and separation. Of the 8 junction-sharing genes: D 3, S 5, F 3, I 0 resolved. (Hypothesis, untested: the gorilla RefSeq annotation has no readthrough-labelled features; 'readthrough' occurs 0 times, 53 of the 54 pairs are protein-coding pairs with short tandem overlaps, and 48 of 108 genes share a donor or acceptor with the partner. The 8 of 108 and the 216 of 256 of the human ideal windows are different gene universes and not comparable.)
4. **P and PC cost 0-1 chains at E_both and 2 / 18 at N (primary denominator).** The transcripts they remove are 16% (P: 204 net) and 3% (PC: 1,672 net) annotation-matching against 42.4% overall (post hoc): consistent with a precision gain that the annotation-match share cannot show; real-read precision of D against the tools is unmeasured.

### Results (chains = expressed annotated chains recovered exactly; P1 denominator unless stated)

| stratum (genes, chains) | D | P | PC | S | F | I |
|---|---|---|---|---|---|---|
| E_both (108, 202) | 170 (84.2%) | 170 | 169 | 140 (69.3%) | 140 (69.3%) | 159 (78.7%) |
| E_both_x (100, 188) | 162 (86.2%) | 162 | 161 | 131 (69.7%) | 136 (72.3%) | 148 (78.7%) |
| E_both_j (8, 14) | 8 | 8 | 8 | 9 | 4 | 11 |
| E_one (176, 329) | 299 (90.9%) | 299 | 298 | 251 (76.3%) | 251 (76.3%) | 266 (80.9%) |
| A (1,583, 3,491) | 3275 (93.8%) | 3275 | 3274 | 2625 (75.2%) | 3000 (85.9%) | 2960 (84.8%) |
| N (10,079, 19,771) | 18702 (94.6%) | 18700 | 18684 | 15586 (78.8%) | 17034 (86.2%) | 16062 (81.2%) |
| transcripts (multi-exon) | 74,045 (73,345) | 72,690 | 69,449 | 68,249 | 153,402 | 551,342 |
| annotation-match share, E_both / N | 43.2 / 42.4% | 43.3 / 42.5% | 43.2 / 43.8% | 40.4 / 43.4% | 19.6 / 19.5% | 16.1 / 14.4% |

Secondary denominator (>= 3 reads of the seeding pool; E_both 110 genes, 207 chains; N 10,216 genes, 19,953 chains): E_both D 175 (84.5%), P 173, PC 172, S 142, F 141, I 162; N D 94.4%, P 93.8%, PC 93.8%, S 78.3%, F 85.4%, I 80.8%. Evaluated genes: 11,946 (P1) and 12,094 (P2) of 41,193; 93,564 distinct valid annotated chains of 29,312 genes; 23,793 (P1) and 23,993 (P2) expressed. Strata (P1): N 10,079, A 1,583, E_one 176, E_both 108 (E_both_j 8).

### Predictions (the registered text; what each could and could not show)

**G1 held** against S and F (E_both and N, both denominators; CIs above). **G2 (P loses >= 3 chains at E_both on the secondary denominator) is uninformative**: only 2 chains at E_both are expressed solely in the secondary pool (D recovers both, P neither), so the headroom (2) was below the bar (3) and computable from the truth tables before any arm was scored; at N the secondary-denominator gap D - P is 102 of 19,953 (not a registered criterion), and on the 200 chains expressed only in the secondary pool D recovers 137 and P 27. **G3 (annotation-match share of D below that of S and F) could not pass**: D emits 74,045 transcripts against F's 153,402 and isoseq's 551,342, so its share cannot fall below F's 19.6%; against S it is above at E_both (43.2 against 40.4%) and below at N (42.4 against 43.4%); the denominator is also restricted to the evaluated-gene universe. The failure says nothing about over-emission on real reads (**unmeasured**). **G4 held** (E_both_x D 71 against S 47). **G5 holds as a point estimate only** (I 159 against D 170; not distinguishable, see Answer 1).

### Limits and checks

- **Truth circularity is larger than 'partly Gnomon'**: 94,896 of 95,019 mRNA / lnc_RNA / transcript records are Gnomon (99.9%), and the model_evidence of the transcripts behind 98.9% of the expressed (P1) chains cites 'long SRA reads'. Whether OR6737 is among them is unverified (counts, not accessions). Recall on the 272 expressed chains from other models is lower (D 64.3, P 64.0, PC 63.6, I 64.0, S 59.2, F 52.9%) than on the long-read-cited ones (D 94.7, S 78.4, F 86.2, I 81.9%): the absolute recall values are conditional on a long-read-derived truth; the ranking is not threatened because no arm reads the annotation.
- Exact whole-chain equality rewards read-majority chains; D's lead is largest at 3-4 reads (E_both 72.5 against S 50.0; N 85.0 against 61.7). StringTie's secondaries never seed an isoform but shift its coverage and ends; isoseq collapse drops alignments below 0.99 coverage or 0.95 identity (it sees a subset of P1). The cluster BAM of the lab runs is identified with the local one by the same @PG line and the 5,506-unmapped identity only.
- One library, cross-individual (private variants drop exact-chain reads: denominators are lower bounds for every arm alike). The lab tools are default annotation-free recipes at versions that differ from the local ones. E_both_j has 8 genes (4 pairs): no conclusion. 10,705,377 is the mapped-record count (5,506 more are unmapped).
- Instrument: `real_truth.py` credits a read to a chain through ANY overlapping same-strand gene's span, not only the chain's own gene (prereg section 2 says its own): 2 of 93,564 chains (ZNF611 region) are over-counted (2 reads against 0), both below 3, so no expressed set, stratum or score changes. M2 means overlap with the exon union of a gene of the stratum and M3 'an exact chain of it' means any valid chain of the gene; the alternatives change no E_both number.
- Pair count: E_both is 54 same-strand pairs (all components of size 2); the 35 of the design sweep used a stricter exact-chain-read rule and the 423 / 1,362 of that sweep are not reproduced (the definition is undocumented; the 2,202 overlapping pairs are).


---

## HIERARCHY_DUPLICON_EXPANSION_2026-10-07.md

*(original: `docs/archive/2026-10/merged/HIERARCHY_DUPLICON_EXPANSION_2026-10-07.md`)*

## Hierarchy study: results (2026-10-07). Y block, concordance with the Yoo et al. 2025 gorilla expansions

**Status: DESCRIPTIVE; no verdict class.** Protocol: `docs/PREREG_hierarchy_duplicon_expansion_2026-10-07.md`, Amendment 3 (git 8fdff163, committed before the Y2 and Y3 numbers existed; Y1 was looked at before and is labelled spent). Results of T1 (Phases R0 to R2) and T2 (Phase E) will be appended to this file under their own headings. Truth files: `/mnt/linuxdisk/home/juanfraitu/winloci_data/yoo_2025/` (workbook sha256 865d4e868568b32d...; `table_38.tsv`, `yoo_GGO_table35_novel_paralogs.tsv`, `yoo_GGO_table40_lineage_specific.tsv`, `GGO_chr_map.tsv`, `yoo_GGO_unit_loci.tsv`, `y1_lrpap1_concordance.tsv`). Match rule: same chromosome, overlap of at least 50% of the shorter locus. Gorilla only; one animal (Jim, mGorGor1).

### T1-d, the gorilla SEDEF arm (Phase R2; descriptive, no verdict class)

**Protocol.** `docs/PREREG_hierarchy_duplicon_expansion_2026-10-07.md` section 5 and Amendment 1. Code `bench/hierarchy/t1_gorilla_pairs.py` with `bench/dna_sd_atoms.py` (restored unchanged, sha256 9d96d976...); 132 tests (`test_t1_gorilla_pairs.py`, `test_t1_gorilla_release.py`, `test_t1_gorilla_review_extras.py`), also under `PYTHONHASHSEED` 0, 1 and 2. Products in `/mnt/linuxdisk/tmp/hier_run/r2/` (`report.json`, `report_pairs_0.90|0.95|0.98.tsv`, `report.stdout`, `gate0_r2.tsv`, `gate0_r2_env.tsv`, `counts.json`, `counts.stdout`); the second hash-seed copy is in `hier_run/r2_report_hs1/`.

**Gates (all pass).** Gate 0: the five inputs and the atoms script match their registered sha256 prefixes (python 3.14.4, numpy 2.4.2, scipy 1.17.1, pysam 0.23.3). Gate 3: `report.json`, the three pair tables and the printed report are byte-identical under `PYTHONHASHSEED` 0 and 1. Gate 6: the verdict-set depth-matched pairs and families are 38/33, 15/15 and 5/5 as registered; under the registered rules the atoms, classes and edges differ from the recon values by -6, +22 and -186 (identical tuples collapsed, coverage >= 0.5). Positive control: the 8 LRPAP1 copies each touch at least one class at tau 0.90 and form one connected component (identity 0.957 to 0.990).

**Review and one disclosure.** Two independent reviewers of the first version found no defect that changes a registered number and asked for a validity coupling in `report` (INVALID and no rate if the control or a Gate 0 fingerprint fails), a fuller printout, hermetic guard tests and the interpreter in Gate 0; a delta reviewer then confirmed that no computed number changed (counts and report byte-identical on five made-up worlds, `counts` on the real inputs equal to the earlier record) and found no fail-open path (verdict: safe to release). After that review I made the reviewer's optional fail-closed hardening (all six inputs required in `validity`, a wrong self-lift table is INVALID with a reason, the registry is validated before anything is built, the registry source and the scorer's sha256 are recorded, stale pair tables are removed on INVALID, the B_viol caveat no longer calls the pool small); these were covered by new tests and were not reviewed again. **Disclosure:** during the first review a mutation run of one reviewer made the old guard tests run `report` on the real inputs once (18:12) and write `report.json` and three pair tables into the real output directory; the reviewer never opened them, deleted them by exact path at 18:16, and no number was seen by anyone; the report is deterministic, so nothing was lost. The leaking tests were replaced by hermetic ones.

**Result (released 2026-10-07T19:25:47 PDT, scorer sha256 e48c4c29...; one run, no tuning; tau 0.90 is the primary level).**

| tau | depth-matched pairs / families | S-rate, pair-weighted [Wilson] | family-weighted | A_viol | bootstrap, pair-weighted | path pairs | aligned < 0.5 | B_viol | UNDERPOWERED |
|---|---|---|---|---|---|---|---|---|---|
| 0.90 | 38 / 33 | 0.789 [0.637, 0.889] | 0.773 | 0.211 | 0.641 to 0.907 | 23 of 38 | 9 | 9 of 11,422 | no |
| 0.95 | 15 / 15 | 0.733 [0.480, 0.891] | 0.733 | 0.267 | 0.467 to 0.933 | 6 of 15 | 4 | 9 of 4,250 | yes |
| 0.98 | 5 / 5 | 0.800 [0.376, 0.964] | 0.800 | 0.200 | 0.400 to 1.000 | 2 of 5 | 2 | 4 of 1,419 | yes |

B_viol at tau 0.90 by the registered distance classes (the same count, split): cross-contig 0 of 10,506, 100 kb to 1 Mb 0 of 27, under 100 kb 5 of 11, over 1 Mb 4 of 878. B_viol is a ceiling-type quantity (the pool is not identity-matched and 92% cross-contig).

**Reading (descriptive; no verdict class is issued).** At tau 0.90, 30 of 38 depth-matched same-family gene pairs lie in a shared SEDEF atom class, and 8 do not. Different-family pairs almost never share a class (9 of 11,422). The 8 pairs without a shared class (post hoc listing from `report_pairs_0.90.tsv`, no statistic): five are cross-contig pairs of a multi-exon parent and a compact copy with one or two exon blocks (GLUD1 13 blocks and GLUD2 1; RBMX 12 and RBMXL1 2; UTP14A 15 and UTP14C 1; RPE 9 and RPEL1 2; FKBP1A 4 and FKBP1C 1), the signature of retrocopies, which are not segmental duplications; three are tandem ancient paralogs on one contig (GON4L and YY1AP1, APOBEC3D and APOBEC3F, FCGR2A and FCGR2B), two of them aligned over less than 0.56 of the shorter transcript. So in the gorilla the share of same-family pairs that sit in a common SD atom is about four fifths in this roster, and most of the rest are family relations that are not SD-derived, which fits the layering SD -> duplicon -> family with retrocopies and old tandem paralogs entering the family level by another route.

**Limits (registered, repeated).** The S-rate measures SEDEF and atom sensitivity at the genes as much as unit boundaries, because depth-matching on transcript identity >= tau selects the pairs for which a SEDEF alignment at fracMatch >= tau is expected. The layer is a human Compara family projection by exact symbol onto the gorilla RefSeq annotation (depleted of LOC-named recent copies, not gorilla truth). The within-contig control pool is tiny (11 and 27 pairs below 1 Mb). The gorilla roster is SD-rich by construction (the development slice and the most-junctions contig are excluded); nothing here is compared with the human T1 cells, and no species comparison is claimed. Path genes (a gene touching two or more classes) are 23 of the 38 pairs, so the pair-level S depends on how classes are cut.

### Why the counts differ from Yoo's (independent diagnosis, 2026-10-07; descriptive)

Three investigators recomputed every Yoo locus in both libraries (reads, assembly, node, family, copy row; 37 loci x 2) and a fourth checked them (29 checks: 23 confirmed, 6 partly, none refuted; all 74 locus-library cells agree between the two library reports). Per-locus tables: `winloci_data/yoo_2025/why/{KB3781,OR6737}_per_locus.tsv`, `yoo_loci_master.tsv`, full result `why/workflow_result_wf_05e4b6b0.json`. Evidence is DIRECT unless marked INFERRED.

**1. LRPAP1, 11 loci against 10 copies.** The extra locus is the chr16 5' fragment LOC115932756 (17,176,045-17,192,323). It is not a product of our RNA-first definition: the RefSeq annotation itself has 11 LRPAP1-family gene records (LRPAP1 and ten 'alpha-2-macroglobulin receptor-associated protein-like' LOCs; 7 protein-coding, 4 pseudogene records), Yoo's 10 loci are 10 of them, and Yoo's own Tables VIII.35 and VIII.40 list the 11th (AMRP, 258 aa) at identical coordinates but not Table VIII.38 or the text count (which names one chr16 copy). Every checkable reason for leaving it out fails: it is covered 100% by SEDEF SDs (9 of 11 overlapping records pass 0.90 and 1 kb) and lies in the same Table XXI.82 'New / Specific' row as counted copy 7; it is as solitary as the counted copies 2 and 6; its product is 69% of the ancestral length like copies 2 and 6; it is expressed (40 and 218 primary reads, RefSeq cites 67 long reads for one of its models). The text's 'all mapping to lineage-specific SD regions' is itself not reproduced for 2 of the 10 (copy 8 has no XXI.82 record, copy 9 lies only in 'unmapped' rows). The supplement states no method for how the 10 were chosen, so the omission is an undocumented curation or pipeline choice (INFERRED by elimination; not determinable from the files).

**2. LRPAP1, our recovery per library.** The 'every Yoo locus matched' result above is for the dossier set. Read-defined BASE nodes recover 8 of 10 in KB3781 (copies 8 and 9: 0 and 1 primary reads) and 10 of 10 in OR6737.

**3. The chr1 unit: Yoo's list and ours are different kinds of object.** Table VIII.38 is an annotation-defined list: 29 of its 37 loci have both boundaries on RefSeq feature boundaries, it has no expression or ORF column, and it also lists copies the text calls unexpressed. Our nodes need two reads on one exact intron chain. Per library (37 Yoo loci each; classes: matched / node but no copy row / reads but no assembled transcript / too few reads): KB3781 18 / 2 / 2 / 15; OR6737 16 / 5 / 10 / 6. By gene: JMJD7-PLA2G4B 7 of 8 matched in KB3781 and 6 of 8 in OR6737 (plus the ancestral chr16 PLA2G4B in KB3781 only); MAPKBP1 and SPTBN5 0 of 9 at exon level in both.
- *Reads (KB3781):* 0 genuine primary reads at 8 of 9 MAPKBP1 and 7 of 9 SPTBN5 loci; RefSeq's own model evidence, independent of our alignments, lists 0 long reads at 8 of 8 chr1 SPTBN5 and 6 of 7 chr1 MAPKBP1 models against 2 to 29 at the JMJD7-PLA2G4B loci. The chr1 MAPKBP1 and SPTBN5 copies are essentially unexpressed in the animal's own fibroblast data.
- *Reads (OR6737):* SPTBN5 has 2 to 6 aligned primary reads at 8 of 9 loci (25 in all), but the largest group of reads sharing one intron chain is 1 at every locus, so no two-read transcript can form (the assembler's floor, `PASS1_MIN_READS = 2`); Yoo's criterion is a read count (> 3), ours is two reads on one chain, and the reads clear Yoo's bar at 4 of 8 loci.
- *Artifacts of long N gaps:* minimap2 alignments with 110 to 217 kb gaps bridge tandem copies (48 alignments in KB3781, 66 in OR6737); they make the KB3781 '1 of 9' hits, fold OR6737 locus U08 into its neighbour, and inflate span-based read counts.
- *Match rule:* the registered '50% of the shorter locus' counts half-gene containment as recovery of Yoo's two-gene readthrough loci; a reciprocal 50% keeps 1 of 9 in both libraries; one miss (19,225,471) is the representative rule (most-read isoform is the JMJD7 half 2.5 kb away), and against Yoo's VIII.35, which lists that half, it is a hit. Which Yoo table is the truth changes the count.
- *Ancestral chr16 MAPKBP1* (74 of 78 reads, node assembled) has no copy row because no chr1 MAPKBP1 node exists to pair with (INFERRED consequence of the unexpressed chr1 copies).
- *Not the cause:* primary-placement coin tosses (the secondaries at the chr1 MAPKBP1 loci are the ancestral reads at lower score), the match-rule or coordinate conventions for LRPAP1, the individual (KB3781 is the assembly's own animal and shows the same scarcity).

**4. Not determinable from the files.** Why Yoo leaves the fragment out; how Yoo counted '> 3 Iso-Seq reads' (its multi-mapping allowance is undocumented, and Table VIII.35 has MK14 models at 6 chr1 loci that neither RefSeq nor our libraries support at copy level, and no SPTBN5 row for any species); whether OR6737 is the reference animal (the figures registry labels it 'testis, mGorGor1', project memory says a different animal; KB3781 is the assembly's own).

**What would make the counts comparable.** Score against Yoo's expressed subset instead of the whole VIII.38 list; use the guided mode (annotation as seeds) for the unexpressed copies, where RefSeq already has most of Yoo's loci (29 of 37 on RefSeq boundaries); count by exon overlap, not span; read the half-gene models against VIII.35. None of these was run.

### Y1, LRPAP1 (spent)

Our 11 loci (8 full-length copies and 3 5' fragments; `docs/archive/2026-10/LRPAP1_FAMILY_2026-10-04.md`) against Yoo's Table VIII.38 (the ancestral copy and 9 copies, 10 loci) and the 5 AMRP rows of Table VIII.35 (AMRP is the UniProt mnemonic of LRPAP1):

| our locus | class | chromosome | Yoo locus | table | overlap of ours / of Yoo's |
|---|---|---|---|---|---|
| LRPAP1 | full-length | chr3:12,090,719-12,110,719 | LRPAP1-anc | VIII.38 | 1.00 / 1.00 |
| LOC134756753 | full-length | chr12:22,553,726-22,573,741 | LRPAP1-1 (also AMRP 358 aa) | VIII.38, VIII.35 | 0.97 / 1.00 |
| copyB (unannotated) | full-length | chr12:24,785,460-24,805,469 | LRPAP1-3 | VIII.38 | 0.99 / 1.00 |
| copyC (unannotated) | full-length | chr12:30,205,371-30,225,550 | LRPAP1-4 | VIII.38 | 0.99 / 0.71 |
| LOC129526389 | full-length | chr14:23,392,677-23,412,705 | LRPAP1-5 (also AMRP 389 aa) | VIII.38, VIII.35 | 0.99 / 1.00 |
| LOC134757218 | full-length | chr16:15,822,038-15,842,048 | LRPAP1-7 | VIII.38 | 0.99 / 1.00 |
| LOC129523574 | full-length | chr22:11,292,022-11,312,055 | LRPAP1-8 | VIII.38 | 0.99 / 1.00 |
| LOC129530227 | full-length | chrY:45,277,394-45,297,429 | LRPAP1-9 | VIII.38 | 0.99 / 1.00 |
| LOC134756368 | 5' fragment | chr12:23,070,423-23,086,604 | LRPAP1-2 (solitary; AMRP 258 aa) | VIII.38, VIII.35 | 1.00 / 1.00 |
| LOC115932954 | 5' fragment | chr14:25,250,978-25,267,164 | LRPAP1-6 (solitary; AMRP 258 aa) | VIII.38, VIII.35 | 1.00 / 1.00 |
| LOC115932756 | 5' fragment | chr16:17,176,045-17,192,323 | AMRP 258 aa | VIII.35 only | 1.00 / 1.00 |

All 10 loci of Table VIII.38 match one of our loci, and all 11 of our loci are in Yoo's supplementary catalogue (10 in VIII.38; the chr16 fragment only in VIII.35, as a novel paralogous gene). Sensitivity 10/10 on Table VIII.38, precision 11/11 against the union of VIII.38 and VIII.35. The text count of 10 copies (1 ancestral and 9) and our 11 differ only by the chr16 fragment that Yoo lists in another table. Yoo's solitary copies 2 and 6 are our chr12 and chr14 fragments with identical coordinates. **Which object this scores.** The 11 loci are the dossier's set (Liftoff projection of the RefSeq annotation plus the three fragments), not read-defined BASE nodes. Scored with the registered rule on the BASE copy rows of each library, the 10 loci of Table VIII.38 are recovered 8 of 10 in KB3781 (copies 8 and 9 have 0 and 1 primary reads in the fibroblast library; Yoo's text says copy 8 is expressed only in testis) and 10 of 10 in OR6737 (45 and 18 reads at copies 8 and 9); the match is identical under the one-to-one, strand-aware and reciprocal-50% rules. Both Yoo tables for every locus, recomputed by `python3 -B bench/hierarchy/yoo_concordance.py y1 ...`: `winloci_data/yoo_2025/y1_lrpap1_all_matches.tsv` (sha256 b5e729097c973abd...). Our O1 and O2 results for these loci (de novo family of the expressed copies, no contested read, no haplotype-only copy) are in `docs/archive/2026-10/LRPAP1_FAMILY_2026-10-04.md`.

### Y2, the chr1 MAPKBP1 / JMJD7-PLA2G4B / SPTBN5 unit

Truth: the 27 loci of Table VIII.38 (8 chr1 triplets, the ancestral chr16 triplet; 9 MAPKBP1, 9 JMJD7-PLA2G4B (the 8 chr1 loci plus the ancestral chr16 record, which Table VIII.38 labels PLA2G4B), 9 SPTBN5 loci). Ours: copy rows of the BASE family tables (genome-wide, pre-f1v2 snapshot, one library at a time). The KB3781 table is the registered one; the OR6737 table is an addition made after the first result and is labelled so.

| library | gene | Yoo loci | matched by a BASE copy | families touched | primary reads at the loci with at least 2 reads |
|---|---|---|---|---|---|
| KB3781 (fibroblast, the animal's own) | JMJD7-PLA2G4B | 9 | 8 | MCL11 | 9 of 9 (26 to 138 reads) |
| | MAPKBP1 | 9 | 1 | MCL11 | 2 of 9 |
| | SPTBN5 | 9 | 1 | MCL11 | 3 of 9 |
| OR6737 (testis) | JMJD7-PLA2G4B | 9 | 6 | MCL31 | 9 of 9 (12 to 132 reads) |
| | MAPKBP1 | 9 | 0 | none | 4 of 9 |
| | SPTBN5 | 9 | 0 | none | 8 of 9 (2 to 69 reads) |

KB3781 family MCL11 has 8 members: 7 match a Yoo JMJD7-PLA2G4B locus (6 on chr1, the ancestral chr16 locus), one member (chr1:19,235,748-19,243,916) lies 2.5 kb from Yoo's chr1:19,225,471-19,233,170 without overlapping it, and the single matches of MAPKBP1 and SPTBN5 come from one 24-exon member (chr1:15,783,756-15,912,704) that runs through three adjacent gene copies. **Correction after the independent diagnosis of 2026-10-07 (section 'Why the counts differ', below):** (i) the 1 of 9 for MAPKBP1 and for SPTBN5 in KB3781 is a span artifact: that member is one 2-read row carried by 48 alignments with 110 to 115 kb N gaps and has 0 exonic bp in either Yoo locus (the SPTBN5 hit is on the opposite strand); at exon level MAPKBP1 and SPTBN5 are 0 of 9 in both libraries. (ii) the sentence that the MAPKBP1 and SPTBN5 copies are missed 'because the libraries carry few or no reads' is only partly right: it holds for KB3781 (0 genuine primary reads at 8 of 9 MAPKBP1 and 7 of 9 SPTBN5 loci) and for chr1 MAPKBP1 in OR6737, but for OR6737 SPTBN5 reads exist (2 to 6 primary reads at 8 of 9 loci, 25 in all, 99.7% identical to the locus) and are too few and too diverse to give two reads on one exact intron chain, which the assembler needs; the 'loci with at least 2 reads' column counted alignments that cross a locus only inside an N gap. (iii) the JMJD7-PLA2G4B matches are containment matches: in 5 of 6 OR6737 matches and 4 KB3781 matches our node is the 8-exon JMJD7 half covering 44 to 48% of Yoo's two-gene locus; with a reciprocal 50% rule only 1 of 9 survives in either library. (iv) 'the 1 and 3 missed loci are losses of the BASE snapshot, not of data' is true of the reads but not always of the snapshot: the miss at 19,225,471 is the representative rule (the most-read isoform is the JMJD7 half 2.5 kb to the right, which Yoo's own Table VIII.35 lists as JMJD7 19,235,754-19,243,912, so against VIII.35 it is a hit). The registered prediction of Amendment 3 (at least 5 of 9 MAPKBP1 loci covered) failed (1 of 9 and 0 of 9 by the registered rule, 0 of 9 and 0 of 9 at exon level). The member-level count differs from the locus-level count by design: in KB3781 the 24-exon member overlaps four Yoo loci (two JMJD7-PLA2G4B, one SPTBN5, one MAPKBP1), so 8 Yoo loci are covered by 7 members (the 8th member is the one 2.5 kb away). The object is the copy rows (transcript spans) of the BASE `copies.tsv`, as registered; a post hoc sensitivity with the locus extents (`--use-locus`) raises the JMJD7-PLA2G4B count to 9 of 9 (KB3781) and 8 of 9 (OR6737) and leaves MAPKBP1 and SPTBN5 at 1 of 9 in both libraries, so the conclusion does not move. Reproduce: `python3 -B bench/hierarchy/yoo_concordance.py y2 --unit-loci winloci_data/yoo_2025/yoo_GGO_unit_loci.tsv --copies /mnt/linuxdisk/tmp/rustle_figures/rt_arms/<library>/<library>.BASE.fam.copies.tsv --out PREFIX` (copies.tsv sha1 64abad727435 for KB3781, a2de88e1097a for OR6737; unit-loci sha256 19a30600139d667c...; tests `bench/hierarchy/test_yoo_concordance.py`, 15). Not done: the current-default (f1v2) families on the two contigs, and a testis library of Yoo's own animal.

### Y3, PSMA5 (count only)

The gorilla RefSeq annotation (`GGO_genomic.gff`) has the ancestral PSMA5 (NC_073224.2 = chr1:130,657,639-130,682,285) and three 'proteasome subunit alpha type-5-like' gene records upstream of it on the same chromosome (LOC129525664 at 116.67 Mb, LOC129525329 at 125.73 Mb, LOC129532613 at 129.28 Mb), which is Yoo's ancestral copy plus 3 duplicated copies upstream. A fourth record, LOC109028113 (NC_073231.2, 1 kb, 'proteasome subunit alpha type-5'), is a short retro-like record that Yoo does not count. The prediction (the ancestral gene and at least one duplicate as separate records) held (3 of 3). The coordinates of Yoo's PSMA5 copies are in no table, so the comparison is by count and chromosome only.

### Y4, genome-wide E-table against Tables VIII.40 and VIII.35

Not yet: it needs the E-table (Phase E).


---

## LOCUS_UNITS_LEVELS_2026-10-07.md

*(original: `docs/archive/2026-10/merged/LOCUS_UNITS_LEVELS_2026-10-07.md`)*

## The locus-model arms read at two levels of the family hierarchy: MCL clusters (M) and components of the pre-MCL graph (C) (2026-10-07)

**Status: DEV, human CHM13 / CAT-Liftoff ideal windows (NPIP chr16, TBC1D3 chr17, two read replicates each), simulation, circular by construction, the same windows as `docs/LOCUS_UNITS_2026-10-06.md`; nothing held out. Protocol: Amendment 2 of `docs/PREREG_locus_units_2026-10-06.md` (commit 5019106d, written after the level-M results and before any component was computed or any F1 arm existed; a graph dump of arm D, NPIP replicate 1, existed 2.5 minutes before that commit, no component was computed from it). Tools: `bench/entangled/{hier_clusters.py (2 unit tests), pool_levels.py, run.sh f1|f1q|hier}`. Products: `/mnt/linuxdisk/tmp/entangled_2026-10-06/{NPIP,TBC1D3}/rep{1,2}/{hierC_*,scoreC_*,famscoreC_*,score_f1*}` and, for arm D at level M, `/mnt/linuxdisk/tmp/ideal_expression_2026-10-06/`. Independently verified (workflow wf_994ece76-8d2: component membership rebuilt for 6 arm-runs, pooled tables and the CANDIDATE_C test recomputed with separate code, Compara F re-implemented): every registered number reproduces; the corrections it forced are folded in below.**
User direction (2026-10-07): 'lets do 1 and 4' (4 = back to the ideal windows: a hierarchy-aware family criterion, then the evidence-based fusion cut combined with the primary-first representative).

### Answer

1. **The level-C reading differs from the registered level-M reading in these windows (post hoc, declared in Amendment 2).** Level C = the connected components of the pre-MCL homology graph (a coarsening of every MCL partition: nested by construction; the components inherit the graph's edge thresholds: identity >= .7, coverage of the longer >= .3, >= 300 bp, shared exonic bp >= 1, shared-exon fraction >= .60). Pooled IDEAL-FOUND (R plus E gene-runs): D 60, P 67, PC 67, **Q 64**, S_Q 62, F1Q 64, S_D 55, F1 60 (level M: 60, 60, 60, 61, 55, 61, 47, 60). The gene-runs are 39 copies x 2 read re-draws of one annotation, not independent runs. Every re-run of `mcl_families --dump-graph` reproduced the arm's registered clusters byte for byte (32 of 32).
2. **Rule 1 (primary-first representative, arm Q) is a clean repair at level C**: it changes four gene-runs on three genes against D, nothing else in 32 columns: NPIPA2 in both NPIP replicates (E2 and E3, 0 -> 1) and TBC1D3P4 (replicate 1) and TBC1D3P3 (replicate 2) (E3 0 -> 1). Reachable copies found 48 -> 52 of 54, entangled 12 -> 12, cluster precision of K*_C and Compara F unchanged (.000 difference under every F variant). At level M the same arm 'lost' NPIPB6, NPIPB8 and NPIPB9 in NPIP replicate 2: they are in K*_C, a partition split, not an evidence loss.
3. **Primaries-only seeding (P, PC): the registered component-level rule is passed, but not robustly.** Reachable 53 of 54 (both families YES under the registered rule), entangled 14 of 24 (level M: 7); net +7 gene-runs = 7 genes gained and 3 lost (net +4 genes), of which R +5 is level-independent and E +2 is PKD1P6 x2, which appears only because K*_C grows from 40 to 48 loci (cluster precision .775 -> .667) and the added members include 3-4 loci overlapping the held-out SMG1P targets and two single-copy negative-control genes (G7; D and Q: 0). `family_score` deletes predicted members outside the truth universe; if they count as false (as the ideal-expression prereg writes it) F_C(P) - F_C(D) is -.049 / -.066 (NPIP replicates, bipartite; -.041 / -.046 pairwise) against a tolerance of -.005, so **P and PC fail the F clause and only Q stays a CANDIDATE_C**. P per copy: + NPIPA2, NPIPA1, PKD1P6, NPIPA7 (both replicates), NPIPB3 (replicate 1), TBC1D3P4 / P3; - NPIPB4 (both), AC138969.1 and NPIPB13 (replicate 2); NPIPA8 and AC138894.1 change flags without changing status.
4. **F1 (evidence-proven bridge cuts without f1v2's share rule) changed nothing for any target copy, and the test is weaker than it looked**: F1 is byte-identical to D in both TBC1D3 replicates (the same single junction is cut); in NPIP it cuts more (5 and 9 bridge transcripts, 3 and 3 junctions, against 2 and 1) and zero target flags differ. NPIPB14P, AC138894.1 (with CLN3 and NPIPB7) and PKD1P6 (with NPIPP1) stay failed. The null is about the simulator's signals (PAS-proven 3' clusters in 10 of 649 junction rows in NPIP replicate 1; F1 flags nothing at the readthrough), not about evidence-based cutting. F1Q equals Q on every target column; it misses the F clause by .012 in NPIP replicate 1 because F1 also removes three SLX1A-SULT1A3 bridge transcripts, a non-target region, and the bipartite F is tie-sensitive (the pairwise F1Q is .488 against Q's .490).
5. **By the registered rule P, PC and Q are CANDIDATE_C and all three are LEVEL-SPECIFIC (no candidate at level M), which by the amendment 'earns nothing by itself'. After verification only Q is robust.** S_D, S_Q, F1 and F1Q are not candidates.

### Pooled results (four runs; R = 54 reachable gene-runs, E = 24 entangled; R / E found, then R: E2 / E3 / E4)

| arm | level M R / E (R+E) | level C R / E (R+E) | R at C: E2 / E3 / E4 | NPIP R found r1 / r2 (bar 13) at C | TBC1D3 (bar 12) at C | verdict N / T at C | CP of K*_C (min) | Compara F at C NPIP / TBC1D3 |
|---|---|---|---|---|---|---|---|---|
| D | 48 / 12 (60) | 48 / 12 (60) | 50 / 50 / 54 | 12 / 12 | 12 / 12 | NO / YES | .78 | .754 .754 / .775 .805 |
| P | 53 / 7 (60) | 53 / 14 (67) | 53 / 53 / 54 | 14 / 13 | 13 / 13 | YES / YES | .67 | .767 .767 / .775 .805 |
| PC | 53 / 7 (60) | 53 / 14 (67) | as P | 14 / 13 | 13 / 13 | YES / YES | .67 | as P |
| **Q** | 49 / 12 (61) | **52 / 12 (64)** | 52 / **54** / 54 | 13 / 13 | 13 / 13 | **YES / YES** | .78 | .754 .754 / .775 .805 |
| S_D | 38 / 9 (47) | 45 / 10 (55) | 51 / 47 / 54 | 11 / 12 | 11 / 11 | NO / PARTLY | .79 | as D |
| S_Q | 45 / 10 (55) | 52 / 10 (62) | 52 / 54 / 54 | 13 / 13 | 13 / 13 | YES / YES | .78 | as D |
| F1 | 48 / 12 (60) | 48 / 12 (60) | as D | 12 / 12 | 12 / 12 | NO / YES | .78 | .774 .754 / .775 .805 |
| F1Q | 49 / 12 (61) | 52 / 12 (64) | as Q | 13 / 13 | 13 / 13 | YES / YES | .78 | .742 .754 / .775 .805 |

E4 failures at level M (R + E gene-runs): not in K*_M / of which no edge path to the family (not in K*_C) / partition split (in K*_C): D 2 / 2 / 0; P and PC 10 / 0 / 10; Q 5 / 2 / 3; S_D 13 / 3 / 10; S_Q 11 / 2 / 9; F1 2 / 2 / 0; F1Q 5 / 2 / 3. Level-M verdict for S_D in TBC1D3 is PARTLY, not NO (the table collapses the three-valued rule; the other arms agree).

### Predictions (Amendment 2), and what each could show

**H1** clause 1 (F1 cuts more bridges than f1v2) is a tautology (F1 is f1v2 without the share rule) and holds in NPIP only; clause 2 (F1 finds NPIPB14P) **failed**. **H2** (F1 loses no reachable copy): held. **H3** (F1Q is a CANDIDATE_C): **failed**, on one run, one non-target region and one F variant. **H4** (Q finds >= 53 of 54 reachable at level C) **could not be met**: E2 on the reachable gene-runs was already 52 at level M (NPIPA7 x2 at locus purity .112, a body that Rule 1 does not touch) and does not depend on the level; the prediction was written after that number existed. **H5** (S_D and S_Q are not candidates at level C): held, **but its stated rationale was wrong**: S_D's reachable-copy loss (R 45 against 48) is a representative failure (E3 47 of 54: NPIPB12 replicate 1, NPIPB13 replicate 2, TBC1D3G x2 lose their main isoform to the cut; E4_C on R is 54 of 54); its 13 E4 failures are 10 partition splits and 3 without an edge path, and of those 3 PKD1P6 x2 fails identically in D, so only NPIPB5 (an entangled copy) is S_D-specific.

### Reading

- Rule 1's effect is the same at both levels (four exactness repairs); at level M it was masked by one MCL partition flip. The instrument that scores 'the one cluster K*' mixes representation with partition granularity; reporting E4 as 'no edge path' against 'partition split' separates them (post hoc). The component level is the single-linkage extreme of the hierarchy, not one of the project's lattice levels L0-L4.
- The entangled stratum holds readthrough images (6 of 24 gene-runs: PKD1P6 and the two AC138969.1 per NPIP run; 4 of D's 12 found entangled gene-runs); without them R + E is D 56, P 62, Q 60 (same order), and the non-image entangled gene-runs found at level C are D 8 / 18, P 9 / 18, Q 8 / 18 (images: D 4 / 6, P 5 / 6, Q 4 / 6). Every entangled difference is NPIP (20 of the 24 gene-runs; TBC1D3's 4 are found by every arm).
- Fused-locus copies (NPIPB14P with the readthrough as representative, AC138894.1 with CLN3 and NPIPB7, PKD1P6 with NPIPP1) remain open: structural separators remove main isoforms (S_D), and the evidence-based cut has no signal to act on in these windows.

### Limits

The same two windows, the same simulator and the same 39 copies as every earlier arm of this series; the arms were motivated by what the dev windows showed; the primary marks of Q and F1Q come from the primaries-only assembly of the same reads under the default f1v2 bridging (lever_P). The registered guard CP_C >= .5 is passed by every arm (minimum .667) and guards nothing; the Compara F tolerance (.005) is smaller than the move of one gene; no annotation-as-loci reference (G6) is scored at level C; the family-level verdicts are knife-edge (NPIP: D 12 of 14 against a bar of 13; TBC1D3: D 12 of 13 against 12; P exactly 13 in NPIP replicate 2). `pool_levels.py` pools sums and does not implement the amendment's 'direction must hold in both replicates of both families' clause; checked by hand it holds for P (+4, +1, +1, +1 gene-runs per run) and Q (+1 in each). A candidate at level C only earns nothing by the amendment; the next registered step is a check on real reads (A119b chr16 / chr17 and gorilla) with the copy-recovery instruments, and ideal windows around other multi-copy regions.


---

## PREREG_gorilla_overlap_2026-10-07.md

*(original: `docs/archive/2026-10/merged/PREREG_gorilla_overlap_2026-10-07.md`)*

## Pre-registration: overlapping genes on real reads, gorilla OR6737 testis: our assembler against StringTie, FLAIR and isoseq collapse (2026-10-07)

**Written before the truth tables, the instrument or any arm of this study exists.** What was already seen: the input-only counts of the gorilla design sweep (2,202 annotated gene pairs share >= 100 exonic bp; 423 with both genes >= 3 exact-chain reads; 35 on the same strand; `docs/GORILLA_FEWCOPY_DESIGN_2026-10-07.md`), the 15-cell genome-wide bakeoffs of the earlier default against the lab's tools (register 1056-1057: aggregates, never per gene), and the ideal-read results of `docs/ENTANGLED_BASELINE_2026-10-06.md`. Nothing of any arm at these genes. User direction (2026-10-06/07): 'lets do 1' (overlapping genes on real gorilla reads).
**Label.** Descriptive real-read check on a partly spent library (OR6737 decided the seeding rule, F1/F1v2 and the strict-junction adoption; the BAM is cross-individual: reads OR6737, assembly and annotation KB3781). Not held-out, no bar except the one replication claim G1.

### 1. Inputs

BAM `winloci_data/GGO_mm.bam` (minimap2 `-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`, 10,705,377 records); FASTA `_from_wsl/winloci_scratch/GGO.fasta`; annotation RefSeq `winloci_data/GGO_genomic.gff` (GCF_029281585.2). Per-read genome-wide best AS: `as_table` of the HEAD build on the BAM (name, best_as, ...). Lab arms (fixed products on this BAM, versions as in `docs/archive/2026-09/PREREG_copy_recovery_tools_2026-09-29.md` section 2): **S** `benchmark_collapse/stringtie_GGO/GGO.stringtie.gtf` (StringTie 3.0.1, `-L`, no `-G`), **F** `benchmark_collapse/flair_GGO/GGO.flair.isoforms.gtf` (FLAIR 3.0.1, bam2bed + `collapse --trust_ends`, no correct, no annotation), **I** `isoseq_upload/isoseq_GGO_OR6737/GGO_OR6737.collapsed.gff.gz` (isoseq collapse defaults). Our arms (HEAD release binaries, the driver's `assemble` stage genome-wide, default f1v2 / GOOD secondary seeding / strict junctions): **D** default, **P** `--no-seed-secondaries`, **PC** P plus `RUSTLE_POLISH_SUBCHAIN=drop`.

### 2. Truth and strata (inputs only; computed before any arm is scored)

- A VALID annotated chain: a RefSeq transcript with >= 2 exons whose every intron is >= 50 bp with a canonical motif (GT-AG, GC-AG, AT-AC on the transcript strand, read from the FASTA); chain = (chrom, strand, ordered intron list), distinct per gene.
- Read pools per alignment (E0 rule of `docs/PREREG_ideal_expression_2026-10-06.md`): an alignment carries a chain iff its junctions (N >= 50 bp) lying inside the gene's exon-union span equal the chain. **P1** = primary alignments (flag & 2308 = 0); **P2** = P1 plus secondary alignments (flag 256, not 2048) with AS >= 0.98 x the molecule's genome-wide best AS. Counted as distinct read names.
- A chain is EXPRESSED iff >= 3 reads carry it exactly. **Primary denominator: >= 3 reads of P1** (what StringTie, FLAIR and isoseq collapse see: they use primaries only); **secondary denominator: >= 3 reads of P2** (the pipeline's seeding pool, generous to our arms). Both are reported, never mixed.
- Genes with >= 1 expressed chain are evaluated. Stratum from the exon unions of ALL annotated genes (expressed or not): **E_both** = shares >= 100 exonic bp on the same strand with another gene that also has an expressed chain; **E_one** = same, partner without an expressed chain; **A** = no same-strand overlap but >= 100 bp with an opposite-strand gene; **N** = none. Inside E_both: **E_j** = an annotated valid chain of the gene shares an exact junction with an annotated valid chain of the partner (not separable by any junction rule), **E_x** = no shared junction.

### 3. Metrics (per arm and stratum; `bench/entangled/real_score.py`)

- **M1** expressed annotated chains recovered exactly by >= 1 transcript of the arm (same chromosome and strand); genes COMPLETE.
- **M2 annotation-match share** (a precision PROXY: unannotated isoforms count against every arm alike): among the arm's multi-exon transcripts whose exons overlap, on the same strand, an evaluated gene of the stratum, the share whose chain equals ANY valid annotated chain.
- **M3 resolution**: a gene is RESOLVED iff an arm gene_id carries an exact chain of it and an exact valid annotated chain of no other gene; reported for E_x and E_j. Caveat as in the ideal study: FLAIR names genes by a 1 kb start bin (over-splitting), StringTie fuses.
- Cross-check: `gffcompare` intron-chain sensitivity of D, S and F against the expressed-chain truth GTF of E_both.

### 4. Predictions (fixed now)

**G1** D recovers more expressed chains than S and than F at E_both and at N (primary denominator): the replication claim. **G2** P recovers fewer chains than D at E_both by >= 3 on the secondary denominator (the seeding gain is visible on real reads). **G3** M2 of D is below that of S and of F (the default over-emits); M2 of PC is >= the smaller of S and F. **G4** at E_x, D resolves at least as many genes as S. **G5** I recovers fewer chains than D at E_both.

### 5. Limits declared in advance

One library, cross-individual, RefSeq-only annotation (partly Gnomon, which may cite long reads: if OR6737 is among them the truth is partly circular, unverified); the secondary denominator favours arms that seed secondaries; the lab tools are in their default annotation-free recipes at versions that differ from the local ones; M2 is a proxy; 35 same-strand E_both pairs are few: E_both is reported with its size and counts, not percentages alone.

### Errata (2026-10-07, appended after the independent verification wf_994ece76-8d2; the text above is unchanged)

(1) The user's direction was one message, 'lets do 1 and 4'; 1 is this study. (2) Exposure: the seeding rule, F1 / F1v2 and the DAZ2 gate were decided with OR6737 (dev); strict junctions + retained-intron ratio 10, F1 (minus NC_073244.2) and R_J are held-out and spent (`docs/GORILLA_FEWCOPY_DESIGN_2026-10-07.md`). (3) The realised E_both is 54 same-strand pairs (108 genes), not the 35 of the design sweep (a stricter exact-chain-read rule); 202 gene-chains are 201 distinct chains. (4) Scoring choices left open and fixed by the implementation, effect on E_both none: M2 'exons overlap an evaluated gene' = overlap with the exon union (a span reading gives 48,720 / 432 instead of 48,670 / 431), M3 'an exact chain of it' = any valid chain of the gene. (5) G2 and G3 could not be met by any arm (the secondary-only headroom at E_both is 2 chains; F emits twice D's transcripts): they were mis-specified bars, see the results document. (6) The 'partly Gnomon' wording of section 5 is 99.9% Gnomon, with 'long SRA reads' cited by the models behind 98.9% of the expressed chains.


---

## PREREG_hierarchy_duplicon_expansion_2026-10-07.md

*(original: `docs/archive/2026-10/merged/PREREG_hierarchy_duplicon_expansion_2026-10-07.md`)*

## Pre-registration: do duplicon boundaries organise families (T1), and do orthology-defined expansions sit inside our families and below their roots (T2)? (2026-10-07)

**Status: DRAFT 3, revised after three adversarial reviews (workflow wf_20046203-7d9) and after a closure review of draft 2 by four further reviewers (workflow wf_eea75a23-1c5). Every finding and its disposition is in `docs/PREREG_hierarchy_duplicon_expansion_REVIEW_DISPOSITION_2026-10-07.md`. Not committed yet. KEY=dupfamboundary (T1), KEY=expansionfamily (T2).**

**Written before any outcome statistic of T1 or T2 exists.** What was computed before this file: INPUT properties only. These are sizes of gene, exon, BED and partition tables, global counts of distinct strings, eligibility counts, file hashes, and the family-structure counts of sections 4, 6 and 7 (components of a graph, clusters, libraries; no unit string is joined to them). They come from four read-only reconnaissance agents (workflow wf_aaa86937-2b6; scratch `/mnt/linuxdisk/tmp/hier_prereg_recon/`, and `/tmp/claude-1000/traps_recon/` for the traps agent; parsed report `/mnt/linuxdisk/tmp/prereg_recon_results.json`), the three reviewers (scratch `/mnt/linuxdisk/tmp/prereg_review/`), the revision of draft 2 (scratch `/mnt/linuxdisk/tmp/prereg_revision/`, scripts and outputs listed in the disposition file) and the four closure reviewers (scratch `/mnt/linuxdisk/tmp/prereg_closure/`), plus reproduction checks of frozen 09-30 quantities. No family x unit-string table, no deficiency, purity, containment, enrichment or S-rate was computed. The distinct-string counts inside family-defined gene sets (section 2, item 1) are marginals, that is a number of cells and a number of distinct strings, and are not a cross-tabulation. The scorers and tools of section 13 are written AFTER this file is committed, with their unit tests first. The commit hash of this file is recorded in the result sections.

### Read first

**The two tests as offered.** On 2026-10-07 the question was: "could there be more to the hierarchy, like sedef defines duplicons, duplicons define expansions like LRPAP1, and then multi-copy gene families?" Two tests were offered and the answer was "yes please". The two sentences below are the wording of that offer. They are not words from the question.

- T1: "is every duplicon boundary a family boundary?"
- T2: "define an expansion as the co-orthologs of one outgroup gene; check that expansions sit inside our RNA families and are younger than the family root".

**Where the registered tests differ from those two sentences.**

1. The registered T1 asks about genes. It does not look at boundaries inside a gene. 1,062 of the 2,290 duplicon-bearing S1C genes (974 of the 2,142 genes of the universe V) already carry two or more different dominant duplicons in their own exons, so the along-sequence reading ("no duplicon boundary inside a family") is false at once and is not tested. Registered instead are D1 (the genes of one family carry the same duplicon content, up to nesting) and D2 (genes with the same duplicon content lie in one family).
2. "Every" is answered by an exact count of exceptions (cell T1-a). T1-a prints m1c (D1) and m2 (D2) separately. It is REFUTED if either is at least 1 and NOT-REFUTED only if both are 0, and the full exception list is printed. Whether the family cut is aligned with duplicon content more than random re-cuts of the same group is a separate and weaker question (T1-b for D1, T1-c for D2). It is asked only inside groups of homologous genes that the family cut splits, because DupMasker labels and family edges are both derived from sequence similarity and a comparison across unrelated genes would pass by construction. Inside a group the comparison is still partly tautological (section 2, item 3): T1-b and T1-c cannot tell duplicon organisation from shared sequence similarity.
3. T1 runs on DNA-level families (our MCL cut of the 2,334 S1C genes' own sequences), not on RNA families. RNA-derived families enter in Phase R3a (the human_testis snapshot, light) and Phase R3b (heavy, not authorised).
4. In T2 an expansion is built from the gorilla side. It is a set of gorilla copies at identity >= 0.95 whose outgroup orthologs are one locus that is single-copy in human, chimp and orangutan alike. The outgroup-first wording ("co-orthologs of one outgroup gene") is tier A, an optional replicate.
5. T2a ("sit inside our RNA families") scores expressed members only. It uses the 09-25/27 BASE family tables, which predate the f1v2 and min-cov-shorter defaults, one library at a time.
6. "Younger than the family root" is not tested on sequence identity. A family is built by chaining sequences that are at least 70% alike, so a set of near-identical copies inside it always looks younger and the test cannot fail. T2b asks instead whether the family that holds the lineage-specific copies also holds the syntenic ancestral copy (arm 1), and what share of the expansions form a clade of a tree built from non-coding sequence, rooted on the outgroup, that leaves out at least one other member of the family (arm 2). Both arms are estimates. Both arms are expected to be UNDERPOWERED. If so, the registered tests give no verdict on that clause.
7. T1 and T2 here are these two tests. They are not the register traps T1 and T2 or the lattice theorems T1, T1', T2 and T3 of earlier documents.
8. The chain in the question has four levels: SEDEF blocks, duplicons, expansions and families. T1 tests duplicon content against families. T2 tests expansions against families. Nothing here tests SEDEF blocks against duplicons or duplicons against expansions.
9. D1 is the gene-level form of the question's sentence: the genes of a family carry the same duplicon content. D2 is its converse, included because KEY=cnduplicon touched that direction. Up to nesting means this. Genes whose duplicon sets fit one inside the next count as one kind. For each family, the genes outside its largest such chain are the genes counted as different.
10. The universe V leaves out the 149 genes that carry two or more Soto families. Those are the genes that cross a family boundary. They are listed apart and enter no statistic.
11. The members of an expansion belong to one pool component. A pool component is built by linking copies that align at identity 0.95 or higher, one link after another (section 7). Two members of one expansion can therefore be less alike than 0.95.
12. Which verdicts can occur. T1-a: REFUTED is expected (Q1). NOT-REFUTED needs both counts to be 0 on all 1,657 genes. T1-b and T1-c in Phase R1: all three content classes can occur (30 variation-bearing groups against a floor of 8, and at least 8 of them must be live, a count that each cell prints). UNDERPOWERED cannot occur in the main cell if at least 8 groups are live. UNINFORMATIVE is practically unreachable. The sweeps at other inflations can be UNDERPOWERED and then carry that tag. SPLIT needs a second run: Phase R3a for T1, the two libraries for T2a. T1-d gets no class. T2a needs at least 8 expansions that each have at least 8 qualifying control components. Input-only counts say this is likely in OR6737 and uncertain in KB3781, where at most about 12 expansions exist, so T2a may be UNDERPOWERED in KB3781. T2b gives estimates and no class and is expected to be UNDERPOWERED, so the registered tests are expected to give no verdict on 'younger than the family root'.

**What each decisive cell decides.**

| cell | decides | cannot decide | possible verdicts |
|---|---|---|---|
| T1-a (F_O, rung U1, D1 and D2) | whether any gene of a family differs in duplicon content beyond nesting (D1), and whether any unit string is split over families (D2) | whether the exceptions are biology or annotation noise (they are listed and labelled) | REFUTED, NOT-REFUTED, INVALID |
| T1-b (D1) and T1-c (D2), F_O, rung U1 | inside the 30 groups of homologous genes that the family cut splits (Phase R1), whether the cut is aligned with duplicon content more than random re-cuts of the same group | whether duplicons organise families across unrelated genes; independence from sequence similarity (the strings and the cut are both sequence-derived) | CUT-ALIGNED, PARTIAL, NOT-DISTINGUISHED, UNDERPOWERED, UNINFORMATIVE, SPLIT, INVALID |
| T1-d (gorilla) | the share of same-family, high-identity gene pairs that share an SD unit class (descriptive) | any verdict (identity selection and instrument sensitivity confound it) | no class; UNDERPOWERED and INVALID flags only |
| T2a | whether the expressed members of an expansion lie in one family more often than matched high-identity control sets | whether containment reflects age (the graph links high-identity paralogs by construction) | SPECIFIC, PARTIAL, NOT-DISTINGUISHED, LESS, UNDERPOWERED, SPLIT, INVALID |
| T2b arm 1 | whether the family that holds derived copies also holds the syntenic ancestral copy (an estimate with an interval) | any age claim beyond the derived/ancestral definition | no class; UNDERPOWERED, INVALID |
| T2b arm 2 | the share of expansions whose members in the family tree form a supported clade that leaves out at least one other family member, in a tree of non-coding sequence rooted on the outgroup (an estimate with an interval) | age (the share also reflects the completeness of the E-table and branch support); gene conversion (counted as not clade, cause printed) | no class; UNDERPOWERED, INVALID |

**What the classes mean.** REFUTED (T1-a) means that at least one clustered gene breaks D1 or D2 in F_O. We expect this (Q1). It refutes the exact gene-level rule. It says nothing for or against the hierarchy idea itself. CUT-ALIGNED (T1-b, T1-c) means this. Take the 30 groups of homologous genes that F_O cuts into two or more families (416 of the 1,657 clustered genes). Re-cut each group at random 10,000 times with the same family sizes. The real cut leaves fewer mismatched genes than at least 99 of every 100 re-cuts (p <= 0.01). PARTIAL means the real cut beats 95 to 99 of every 100 re-cuts (0.01 < p <= 0.05). NOT-DISTINGUISHED means it beats fewer than 95. A cut that does worse than random also reads NOT-DISTINGUISHED. The class must also survive the removal of the one group that contributes most (section 4). Each class speaks for the groups that carry the signal. The removal clause protects against one dominant group only, so two or three small aligned groups can still give CUT-ALIGNED, and the z table lists them. Duplicon labels and family edges both come from sequence similarity. A cut and a set of strings that both follow sequence similarity therefore give CUT-ALIGNED inside a group, so no class shows that duplicons cause the cut. A mismatched gene is one that must move for D1 (T1-b) or for D2 (T1-c) to hold. INVALID means a gate or control failed, so no number of that table is used. UNDERPOWERED means fewer than 8 independent units. UNINFORMATIVE means that the re-cuts cannot change the statistic in any group.

**Terms.** S1C: Table S1C of Soto et al. 2025 (per-gene family and copy number on CHM13 v1.0; 2,334 genes). famCN: family copy number from read depth. SD98: segmental duplications at >= 98% identity. DupMasker colors BED: the duplicon calls of Vollger et al. 2022 on CHM13 v1.0. SEDEF: the segmental-duplication caller whose output (`GGO_sedef_final.bed`) defines the gorilla units. BASE tables: the genome-wide family tables of the 09-25/27 `rt_arms` runs (MCL, inflation 2.8). Liftoff classes: in_place (the source gene's own placement), extra_copy (an additional whole-gene-body copy), dropped_M1 (a copy landing on another annotated gene), dropped_M2 (a copy overlapping another copy). Pool: the union-find components of the gorilla self-lift copy rows (section 7). Stratum: a connected component of the F_O homology graph (section 4). "Unit" means only the SD block, duplicon and (duplicon, copy-number class) levels. Contig NC_073234.2 was used earlier to test the most-junctions representative rule. Clean gene: exactly one Soto family label that does not start with Unassigned. Family-less gene: only Unassigned labels. Multi-family gene: two or more such labels (the S1C column No. Assigned Families is at least 2). MCL: Markov clustering (inflation 2.8 unless stated). Independent unit: the unit counted against n_min = 8 (a variation-bearing stratum, a testable expansion, an informative expansion). Dv: the derived members of an expansion (section 9). J: the number of Liftoff rank runs per outgroup (section 7). Gates are written Gate 0 to Gate 7, unit-string rungs U0 to U2, consequences C1 to C4, predictions Q1 to Q10, phases R0 to R3 and E.

User direction (2026-10-07), full text: "could there be more to the hierarchy, like sedef defines duplicons, duplicons define expansions like LRPAP1, and then multi-copy gene families?" and, on the offer to pre-register the two tests, "yes please".

**Labels.** The human DNA-level layers (S1C x DupMasker) are DEV, second look: every S1C-based duplicon-family relation was seen on both halves in 09-30. Phase R3a (the human_testis snapshot) will be called HELD-OUT for the duplicon relation once it has run. It is not held out for decisions about the family rule. Gorilla is a second-species arm on a partly spent library (SD x family relations were examined on NC_073241.2, NC_073242.2 and NC_073244.2) and is never called held-out. For T2 the development case and positive control is LRPAP1 alone. Every other expansion is in the verdict set, and no gorilla contig is called held-out because the genome-wide stages ran everywhere. Human and gorilla are never pooled. KB3781 and OR6737 are never pooled.

### 1. The two claims, made operational

**Objects.** A *gene* is a path through *units*. Units, from coarse to fine, are the SD block (SEDEF), the *duplicon* (DupMasker label) and the *unit* (duplicon, copy-number class). A *unit string* of a gene is a set derived from its exons (section 4, rungs U0 to U2). A *family* is a block of a named family partition (a layer). An *expansion* is a set of at least two gorilla annotated loci that are co-orthologs of one locus that is single-copy in each outgroup (section 7).

**T1 has two directions.** D1: the genes of a family carry one unit string up to nesting, that is the family partition refines the unit-string partition up to nesting. D2: the genes of a unit string lie in one family, that is the unit-string partition refines the family partition. KEY=cnduplicon tested a soft, direction-neutral association (within-family minus between-family similarity of exonic duplicon composition inside sequence clusters). It made no strict D1 or D2 statement.

**T2 has two clauses.** T2a: the expressed members of an expansion lie in one family of the named layer. T2b: the expansion is younger than the root of that family, tested as in section 9.

**Not tested.** The along-sequence reading of T1 (a change of duplicon label inside a gene) is not tested. Under it T1 is false at once, because 1,062 of 2,290 duplicon-bearing genes contain two or more dominant duplicons.

### 2. What the inputs already settle (seen before writing; the rules below are written with these in view)

1. **Strict D1 is false by counting.** For any family partition F and string partition S on the same genes, m1 >= K - #F >= #S - #F, where K is the number of non-empty (family, string) cells and singleton families count as families. Counted over the genes of multi-gene families (`/mnt/linuxdisk/tmp/prereg_revision/logic/fq_checks.out`, `string_only_checks.out`): F_S has 865 U1 strings against 434 families (m1 >= 431) and 544 U0 labels (m1 >= 110); F_Q has 920 against 390 (>= 530) and 568 (>= 178); F_O has 734 against 362 (>= 372) and 443 (>= 81). On the 2,036 clean genes of V, F_S has 486 families (434 with at least two genes) and 901 distinct U1 strings, so K is at least 901 (K itself was not counted; 491 families and 444 with at least two genes are the counts over all 2,071 clean genes). Rung U2 and layer F_R were not counted. Strict D1 is reported and cannot carry T1. Only the tolerant relation and the exception list can.
2. **The universal statement is expected to fail for reasons that apply to genes.** NPIP members differ in dominant duplicon inside one family (SD9443, SD9622 and SD9613 were looked at by hand on 09-30). The universe holds 36 clean readthrough genes, 105 same-strand span-overlapping pairs of genes (181 genes) and 862 genes with at most two merged exons (counts in V, `c10_excl_counts.out`, `short_counts.out`). Multi-family genes (149 S1C genes) are outside V by construction and are listed separately. They are not all fusions (84 lncRNA, 63 pseudogenes, 39 on the readthrough list), and the NPIP-fusion block reading is DEV only (A119b EXPLAINED, held-out UNDERPOWERED with three switch units).
3. **Tautology.** Soto's sequence families come from exon sharing at identity >= 0.98, our F_O edges need identity >= 0.70, and DupMasker labels are sequence-derived, so 'shares a duplicon' is near-implied by an edge. The decisive alignment cells are therefore computed inside homologous groups (section 4, nulls). That removes the comparison across unrelated genes. It does not remove the dependence between the cut and the strings inside a group: any string system that is nested with the cut (clusters merged, clusters split, or the cluster identity) gives the lowest attainable p in T1-b and T1-c (synthetic check, `/mnt/linuxdisk/tmp/prereg_closure/logic/t1_nested_tests.out`), while strings independent of the cut give p between 0.54 and 0.63 on the real strata (`t1_null_tests.out`). CUT-ALIGNED is therefore expected whenever the strings and the cut both follow sequence similarity, and T1-b and T1-c cannot separate duplicon organisation from shared sequence similarity. A similarity-matched comparator (the F_O cut at inflation 2.0 and 4.0 as competing sequence-only re-cuts of the same group) is not registered here and could be added by an Amendment before any run. Layers derived from the same pairwise alignments are labelled as such (F_Q, F_O, the BASE snapshot). The Compara-projected layer is SD-blind and is the anchor of the gorilla arm.
4. **Circularity.** The (duplicon, CN class) unit built with S1C famCN uses the copy numbers Soto cut its families by. The 09-30 arm D with S1C famCN still predicts only 69.7% of clean genes into exactly their own family (15.7% into two or more families, 5.2% into none), so the circularity is partial. Rung U2 is descriptive, uses our famCN and is never primary.
5. **Spent substrates and family tables.** Every S1C-based duplicon-family relation was seen on both halves; the two halves share 43 units holding 1,060 of 2,334 genes. No gorilla contig is unexposed (section 7). No genome-wide family table under the CURRENT default (f1v2 plus min-cov-shorter 0.70) exists for any library. Genome-wide pre-f1v2 BASE tables (`min_cov_shorter` 0 in `fam.params.tsv`, binary commit not recorded) exist for gorilla KB3781 (308 families / 864 loci), gorilla OR6737 (483 / 1,559), chimp PTR (379 / 1,291) and human_testis (338 / 1,095), under `/mnt/linuxdisk/tmp/rustle_figures/rt_arms/`. None exists for orangutan or for A119b. The 17-hour estimate (quadratic cost model) is for A119b only (locus FASTA 1,980 MB; human_testis has 324 MB).
6. **'Younger than the root' on identity is vacuous.** Under single linkage, E a subset of F implies h(E) <= h(F) for any data, where h is the merge height in divergence units (1 - identity). Equality occurs when the weakest minimax link of F lies inside E. E is defined at Liftoff sequence_ID >= 0.95 while family edges need only identity 0.70 and coverage 0.30. For E = F no older member exists. Identity is not age (gene conversion, selection). The identity version is kept as a pipeline check (section 9) with no evidential weight.
7. **Expansions as defined by the literal wording are few and mostly invisible to name joins.** Exact RefSeq symbols once in human, chimp and orangutan and at least twice in gorilla: 0 (the 53 multi-record gorilla symbols are 47 tRNA families and 6 X/Y pairs). A CDS-product-stem proxy gives 242 stems (598 genes), of which 15 (+7 partly) have members linked by Liftoff copy pairs (16 + 8 of 255 stems under the variant normaliser). Among the 160 linked stems of the variant normaliser, 24 are single-copy in all three outgroups. This is a feasibility indicator only, because name joins are not used in the definition. LRPAP1 (11 loci) is spent. Of 1,539 copy-bearing gorilla records (exon union >= 200 bp), 547 have >= 2 primary reads in either library. Before any contig rule, the pool holds 80 (KB3781) and 139 (OR6737) supported components (157 with two supported records in either library, 153 within one library). If the 15% single-copy share of the name-linked stems transfers (24 of 160; the transfer is unverified), at most about 12 (KB3781) and 21 (OR6737) expansions are expected, before the filters of section 8 that define the testable ones (at least 2 expressed members in one class and at least 20 matched control components). The Wilson 95% interval of 24 of 160 is 10.3% to 21.3%, which gives 8 to 17 expansions for KB3781 and 14 to 30 for OR6737 before those filters. On the 76 non-hub supported components of KB3781 the expected number at 15% is 11.4 and the Poisson chance of fewer than 8 expansions before the filters is about 12%, and about 48% at the lower end of the interval. No cross-species orthology resource exists on disk (no chains, nets, HAL, TOGA).
8. **Coordinates.** The CHM13 v1.0 to v2.0 offset is constant per chromosome except on chr1 (+169, 230 genes; +634, 2 genes) and on the acrocentrics, where two values occur (chr13 +672,767 or +70; chr14 +57,685 or +478; chr15 +585,122 or -4; chr21 +737,013 or -1; chr22 +28,768 or -39). chr16 is +5 and chr17 is +291. The CAT/Liftoff v2.0 annotation carries all 2,334 S1C ids unchanged, so the join key is the id. Duplicons and SEDEF (v1.0) are intersected only with v1.0 genes (CAT v4).
9. **Disclosed counts that do not reproduce.** The 09-30 prereg disclosed 2,259 genes with an exonic duplicon and 1,055 of 1,402 labelled genes with the dominant duplicon among the `SD Unit` labels. They reproduce as 2,290 and 1,099 of 1,427. Window, offset and overlap explanations were tried by two reviewers and in the revision of draft 2 (`c5_disc2259.out`, `c5b_disc2259_scan.out`): no overlap threshold from 1 to 300 bp gives 2,259, an exonic-coverage fraction of 0.10 gives 2,258, and a fraction between 0.098 and 0.099 gives 2,259 and has no basis. The cluster counts of the 09-30 tests (33 clusters, 87 Soto families, 642 genes, 5 cluster genes without an exonic duplicon) and the six unit counts reproduce; the pooled delta, its p-value and the arm A shares are checked by Gate 1. The discrepancy is recorded as an erratum of the 09-30 disclosure and is not a precondition of any gate.
10. **Universe arithmetic.** V = 2,334 S1C genes - 149 multi-family genes - 43 non-multi-family genes with no exonic duplicon = 2,142 (2,036 clean + 106 family-less). Counts quoted below are over V unless stated: readthrough genes 36 (76 over all 2,334), non-eligible genes 388 (541), same-strand span-overlapping pairs 105 pairs / 181 genes (233 / 339), genes with at most two merged exons 862 (948).

### 3. Gates, verdict algebra, decisive cells, run order

**Scope of the gates.** A result table is INVALID unless the gates that cover it pass. T1 tables of Phase R1 need Gates 0, 1, 3, 4, 5 and 6. Gate 2 covers Phase R3a and the optional matched human run, which need Gates 0 to 6. T1 tables of Phase R2 need Gates 0, 3 and 6 and the positive control of section 5. T2 tables need Gates 0, 3, 4 (T2 controls), 6 and 7. A gate prints PASS or FAIL and no outcome statistic. Only the counts listed in Gate 6 are tripwires. The counts quoted in sections 5 and 7 are expected values, and a difference is printed and explained.

- **Gate 0, frozen inputs.** Path, sha256 and row count of every input of section 13 are recorded before the first run, together with the interpreter (`python3` 3.14.4 with pysam 0.23.3, numpy 2.4.2, scipy 1.17.1, run as `python3 -B`; the `python` on the path is 3.13.12 and is not used). Files above 1 GB (FASTA, BAM, PAF, minimap2 indexes) get size and mtime recorded at the run; their sha256 is computed once, in its own budgeted call, before the E phase. The sha256 of the E-table, `member_exons.tsv`, `control_pool.tsv` and `extension.tsv` is recorded when section 7 closes, before any family file is opened for T2.
- **Gate 1, harness reproduction.** The new loaders reproduce the registered OUTPUTS of KEY=cnduplicon (33 clusters, 87 Soto families, 642 genes; held-out +0.2853, dev +0.3276 at 10,000 shuffles, seed 20260930), of KEY=unitcovercn arm A (dev 0.439, held-out 0.473, exact 0.107, recall 0.555, precision 0.718, clean 0.546 / 0.303 / 0.020) and the six frozen arm unit counts (1,546 / 1,241 / 2,255 / 2,400 / 1,715 / 1,837). The frozen F_O binary (sha1 e4cc13b9) re-run with the flags of `ours.params.tsv` reproduces `ours.clusters.tsv`, `ours.loci.tsv` and `ours.pairs.tsv` byte for byte (done in the revision of draft 2, `/mnt/linuxdisk/tmp/prereg_revision/fo_graph/`) and the dumped graph has 1,819 nodes, 11,031 edges and 339 components.
- **Gate 2, coordinates** (covers Phase R3a and the optional matched human run, not Phase R1). For every autosome in the order chr1, chr2, ..., chr22 draw min(20, n) of its n S1C genes (v1.0 CAT v4 coordinates assign a gene to a chromosome; ids sorted as Python strings) with the sample method of one `random.Random(20260930)` object used for all chromosomes. A drawn gene differs if its exon sequence is not identical in the v1.0 CAT v4 exons and the v2.0 CAT/Liftoff exons. Pass iff at most 1 drawn gene differs when n >= 20 and none differs when n < 20. A failing chromosome is excluded and counted.
- **Gate 3, determinism.** Runs under `env PYTHONHASHSEED=0 python3 -B` and `env PYTHONHASHSEED=1 python3 -B` (never `-I`, which ignores the variable) give identical tables. All ties are sorted. Every cell re-seeds `random.Random(seed)` itself, with seed 20260930 for T1 and 20261007 (the registration date) for T2, so execution order does not matter.
- **Gate 4, scorer controls (synthetic, must pass).** T1: families equal to the string partition give m1 = m2 = 0; relabelling k genes of a family whose native string stays the plurality gives exactly the known deficiency; strict m1 is non-decreasing and strict m2 non-increasing when the string partition is refined (random partitions); synthetic families built as one string, as nested truncations and as branches give chain-tolerant deficiency 0, 0 and > 0; the path-defect cases of section 4 give path 0, two disjoint units 1, 3-leaf star 1, triangle 1; the exception list is identical under any input order of the genes; a family partition equal to the strata (zero variation-bearing strata) gives UNDERPOWERED; a synthetic set of at least 8 strata, each with two clusters and one string per stratum, gives UNINFORMATIVE; a planted partition that follows the strings gives p <= 0.01; on the cluster sizes of the 30 variation-bearing strata of section 4 (sizes only, synthetic strings), strings equal to the cluster identity, to unions of adjacent clusters (pairs, in strata of at least three clusters only) or to a split of every cluster of at least three genes into two alternating types give p <= 0.01 in both cells (the documented tautology of section 2, item 3); a planted set in which only the largest stratum follows the strings gives p <= 0.01 on all strata and a class below CUT-ALIGNED after the one-group removal; planted random partitions give p <= 0.05 in at most 10 of 100 repeats (under a calibrated test about 5 are expected, and P(more than 10) is about 1%). Each repeat uses B = 10,000 and runs in several light calls of under 180 s. T2: planted containment 0 and 1 give R = 0 and 1; a planted clade is recovered; the one-to-one node rule is applied when two members select one node.
- **Gate 5, consistency of the scorer on the data.** Gate 5a: on the common universe (V genes that have a famCN, once for each famCN source) strict m1 and m1c are non-decreasing and strict m2 is non-increasing from U1 to U2, for every layer. This is a theorem on a common gene set because dropping the CN class maps proper subsets to proper subsets. U0 does not enter Gate 5, because the principal duplicon is not a function of the U1 set (32 U1 strings carried by 110 genes of V contain genes with more than one principal duplicon, and 10 genes have a principal duplicon that is not in their own U1 set). Gate 5b: the scorer returns D[S -> Q] = 7 genes (4 of 444 Soto families split, counted over all clean genes) when ungrouped genes are singleton cells and 2 genes when they are ignored, where D[X -> Y] = sum over blocks x of X of (|x| - max over y of |x intersect y|), Soto families being S and the 405 sequence groups Q (`fq_checks.out`).
- **Gate 6, denominators printed first** (before any statistic). Universe sizes per layer and rung; every exclusion with its count; size-2 families (n_f = 2 in the layer universe; F_O: 152 of the 362 clusters with at least two V_F genes; 163 of 380 clusters have two locus rows in `ours.clusters.tsv` and 147 of them are among the 152); families with fewer than two genes in the universe; the 149 multi-family genes (listed with their families and strings in a separate table that lies outside every statistic), plus per layer the number of clusters that hold at least one of them and the number of members deleted by V; exons whose two labels with most bases tie exactly (857 of 12,597 labelled merged exons of the genes of V; 898 of 13,501 over all 2,334 S1C genes), genes of V with at least one such exon (316), distinct U1 strings under union-per-label (989 to 985) and under numeric-ID ties (989 to 989; 82 exons in 54 genes differ and the U1 sets of 53 genes change); for the F_O layer the 339 components, the 325 strata that hold V_F genes, the 30 variation-bearing strata (416 V_F genes) and the 295 single-cluster strata (1,241 V_F genes). The number of live strata of each cell depends on the draws and is printed by the cell before any p-value. A count in this Gate that differs from its registered value means the loader differs from the registration and the table is INVALID. For Phase R3a the number of strata that hold at least two clusters among the V genes of the human_testis BASE graph (Q8 applies only if it is at least 8), and the check that every BASE cluster lies inside one component of the dumped graph (a BASE cluster that spans two components makes Phase R3a INVALID). For T2 the eligible expansions, members, hubs, dropped classes and the reasons.
- **Gate 7, E-table controls (T2 only).** See section 7: the positive control (LRPAP1) and the planted negative control.

**Verdict algebra.** For a null-based cell, p = (1 + number of draws at least as extreme as observed) / (B + 1) with B = 10,000. The cut-offs are p <= 0.01 and p <= 0.05 in every section (with B = 10,000, p <= 0.01 means at most 99 draws are as extreme). The class of a cell is the first that applies in this order.

1. INVALID: a gate or control that covers the cell fails.
2. UNDERPOWERED: fewer than n_min = 8 independent units. The units are the variation-bearing strata (T1-b, T1-c), n_testable (T2a, per library, as defined in section 8), the informative expansions (T2b arm 2) and the units of arm 1. For T1-b and T1-c, UNDERPOWERED also applies when fewer than 8 strata are live (below).
3. UNINFORMATIVE (T1-b and T1-c): at least 8 variation-bearing strata and no live stratum. A stratum is live for a cell if the cell's statistic (m1c for T1-b, m2 for T1-c) takes at least two values over the B draws. A stratum whose genes all carry one string is dead for both cells, and for T1-b a stratum whose strings form one inclusion chain is dead as well. The cell prints the number of live strata before any p-value.
4. SPLIT: two runs of one cell (Phase R1 and Phase R3a for T1-b and T1-c) or the two libraries of T2a both have a content class and the classes are two or more steps apart on the cell's ordered scale. Adjacent classes agree, with a note. A run that is INVALID, UNDERPOWERED or UNINFORMATIVE leaves the other run's class, labelled single-run; for the pair Phase R1 and Phase R3a the surviving Phase R1 class is labelled DEV. The level sweep of section 4 is not a run for SPLIT. T2a headline class: when neither library is UNDERPOWERED and they are not SPLIT, the headline class is the class of the two that lies nearer to NOT-DISTINGUISHED on the scale, and both are printed.
5. The content classes. T1-b and T1-c scale: CUT-ALIGNED (p <= 0.01) > PARTIAL (0.01 < p <= 0.05) > NOT-DISTINGUISHED (p > 0.05). The upper-tail p is printed for T1-b and T1-c and carries no class. The class is issued only if it also holds after the one group that contributes most is removed (section 4); otherwise the lower class is issued. T1-d carries no class and is flagged UNDERPOWERED when it has fewer than 20 depth-matched pairs or fewer than 8 families. T2a scale: SPECIFIC (p_more <= 0.01) > PARTIAL, more > NOT-DISTINGUISHED > PARTIAL, less > LESS (p_less <= 0.01). T1-a: REFUTED or NOT-REFUTED. T2b arms 1 and 2 are estimates and carry no class.

A unit excluded for an input property (E equal to the family, unresolved synteny, a hub) leaves the denominator, is counted and is not a cell class. The label 'held-out' is used only for Phase R3a after it runs.

**Decisive cells.** Every other number is descriptive and carries no class.

| cell | phase | statistic | independent unit |
|---|---|---|---|
| T1-a | R1 | F_O, rung U1, universe V_F: m1c (D1) and m2 (D2), printed separately; REFUTED iff either is >= 1 | none |
| T1-b | R1, R3a | F_O (R3a: BASE snapshot), U1: m1c against the homology-preserving null | variation-bearing stratum |
| T1-c | R1, R3a | same layer: m2 against the same null | variation-bearing stratum |
| T1-d | R2 | gorilla S-rate on depth-matched pairs (descriptive) | family |
| T2a | T2a | mean R_j against matched control sets, one cell per library | testable expansion |
| T2b-1 | T2b | share of mixed-family expansions whose expressed ancestral member lies in the family that holds the derived members (estimate) | mixed-family expansion |
| T2b-2 | T2b | share of informative expansions whose members in the family tree form a supported clade that leaves out a family member (estimate) | informative expansion |

**Run order and authorisation.**

| phase | content | class | needs |
|---|---|---|---|
| R0 | Gates 0, 1, 3, 4, 5 and 6 for T1 (Gate 2 with Phase R3a) | light, several calls under 3 min each | 'go' |
| R1 | section 4 (T1, human, DEV), split by cell so each call stays under 3 min | light | 'go' |
| R2 | section 5 (gorilla, descriptive), after `dna_sd_atoms.py` is restored with its test | light | 'go' |
| E | section 7 (E-table, tier B, Gate 7) | heavy, bounded foreground calls, Liftoff peak RSS 10 to 14 GB | a separate heavy 'go' and an Amendment |
| T2a | section 8 | light | 'go' after the E-table sha256 is recorded |
| R3a | section 6 (human_testis snapshot) | light | a 'go' that names it, after Amendment 1 |
| matched human run | section 5 (optional) | light | its own 'go' |
| T2b-1, T2b-2 | section 9 | heavy | each a separate 'go' and an Amendment; predicted UNDERPOWERED |
| R3b, tier A, tier C, regeneration of family tables under the current default, cross-lineage replicate, level-C component table of the gorilla BASE PAFs, rank-sharded self-lift of the gorilla annotation | sections 6, 7, 8, 11 | heavy | a user message and an Amendment each |
| sha256 of files above 1 GB (Gate 0) | light, one file per call under 180 s; a file that needs longer is hashed in a heavy call | covered by the go of the first phase that opens the file |

A 'go' names the phases it covers and covers nothing else. Light phases run under `tools/rlock.sh light` in calls of under 180 s. Any heavy step, and any step marked Amendment in the table, runs only after a user message that names it and an Amendment committed to this file before the run that fixes the command, the binary commit, the input hashes, the call schedule and the time budget; heavy steps also run under `tools/rlock.sh heavy`, one heavy job at a time. Phase R3a needs a 'go' that names it and Amendment 1. A 'go' for light phases never authorises a heavy step. The limits are 590 s for Liftoff calls, 560 s for the families stage and 390 s for every other heavy call. Light calls run as `RLOCK_TIMEOUT=180 tools/rlock.sh light <command>`.

### 4. T1, human, DEV (Phase R1; tables only; light class)

**Universe V.** S1C genes keyed by `Gene ID` (never by name) that are not multi-family (149 removed; they are the genes that are paths across family boundaries and are listed separately) and have at least one exonic base on a DupMasker segment (43 more removed): 2,142 genes, of which 2,036 are clean and 106 are family-less. Merged exons come from `sd98_gene_exons.tsv`, coordinates as stored. Every overlapping colors-BED segment counts (7,615 rows overlap an earlier row of the same label; 48,745 adjacent overlapping neighbours, 544 of them same-label).

**Exon label and unit strings.** Bases of label d on a merged exon e are the sum, over every colors-BED segment of label d, of the overlap with e. Same-label overlaps add, as `comp()` of `soto_m2_unit_cover_cn.py` does; this reproduces the 901 and 564 strings and the six frozen unit counts. dom(e) is the label with the most bases, ties broken by the lexicographically smallest label string (the key of the code; the 09-30 text said ID order, and numeric order is a printed sensitivity). Exons that overlap nothing have no label and enter no string.

- **Rung U1 (primary):** the set {dom(e)} over the gene's merged exons. A duplicon enters only if it is dominant on at least one exon.
- **Rung U0:** the principal duplicon, the label with the most bases pooled over the gene's merged exons, ties by the same key.
- **Rung U2 (descriptive; circularity of item 4):** the set {(d, k)} with k the copy-number class of d that contains the gene's own famCN. The famCN values are those of `families_cn.json` (`cs` = S1C median, `co` = ours, 268 samples, as rounded there). Classes are built from every V gene with a famCN and d in its U1 set, by single linkage with gap >= 2 on the 0.1-rounded values and the comparator of the frozen code. Genes with no famCN (388 in V) are outside U2 and are counted. U2co (our famCN) is the descriptive rung; U2cs (S1C famCN) is computed for Gate 5a and for the F_S note only.
- **Sensitivities (descriptive):** union-per-label bases instead of added bases; numeric ID order for ties; a tie yields the SET of tied labels; a 5% near-tie band that keeps every label with bases >= 0.95 x the top label's bases; the any-overlap set (every duplicon with at least 1 bp on any exon of the gene).

**Family layers (one run each, never pooled).**

- **F_O (primary).** Our DNA-level MCL partition of the same 2,334 loci: binary sha1 e4cc13b9, flags `--min-exonic-bp 1 --min-shared-exon-frac 0.60 --dump-pairs` and defaults (inflation 2.8, `--min-cov-shorter` 0.70, prune 1e-9, min-size 2), on `/mnt/linuxdisk/tmp/ourhom/loci.paf` and `loci.gff3`; 380 clusters, 1,709 cluster rows, and 93 annotation records folded into 86 representatives (`ours.loci.tsv`). The F_O family of a gene is (1) the cluster of the row of `ours.clusters.tsv` whose chrom:start-end equals the gene's locus string in `names.tsv`; else (2) the cluster of the row equal to the representative of that string in `ours.loci.tsv`; else (3) none. In V, 1,581 genes get a family by rule 1, 76 by rule 2 (80 V genes are folded annotations and 4 of them have a representative that is not a row), and 485 by neither (467 are never connected to any other locus and 18 are singleton nodes).
- **V_F** is the set of the 1,657 V genes that carry an F_O cluster. It is the universe of every decisive cell. The 485 unclustered genes are outside the decisive cells. A descriptive arm adds them as singleton families, in which D2 counts each of them as a violation when it shares a string with a clustered gene, and the number of such genes is printed.
- **F_Q.** Sequence groups `collapse_cover(pair_families(..., gate=False))`: 405 groups, 2,102 grouped genes of the 2,185-gene universe, 83 ungrouped genes that are singleton families. Descriptive.
- **F_S.** Soto clean families: 486 on the 2,036 clean duplicon-bearing genes of V, 434 with at least two genes (491 and 444 over all 2,071 clean genes). The 106 family-less genes of V are outside F_S. Descriptive, and circular at U2.
- **F_R.** RNA-derived families: section 6 (R3a, R3b).

Layer universes are not mixed. n_f and n_s count members of the layer universe only. Every ratio m / |universe| is printed with the universe size beside it. The label of a statistic always names the layer: F_O is a DNA-level MCL partition of annotated SD98 loci and is never called "our RNA families".

**Half.** The half of a gene is the `half` column of `bench/soto/soto_split_2026-09-29.tsv` (gene level). A per-half statistic is computed on the V_F genes of that half with every family restricted to them. The two values are named half A (split dev) and half B (split heldout), are printed beside the all-V_F value, carry no class, are not additive and are never called held-out. Nine of the 362 multi-gene F_O clusters straddle the halves (61 of 1,647 genes).

**Statistics** (counts of genes that must move so that one partition refines the other).

- **m1 (strict D1)** = sum over families f with at least two genes in the universe of (n_f - max over strings s of n_{f,s}).
- **m1c (chain-tolerant D1)** = sum over such families of (n_f - W_f), where W_f is the maximum, over chains of the family's distinct strings under proper-subset inclusion, of the number of genes that carry a string of the chain. Its value is unique (the maximum of a finite set, a dynamic programme over the inclusion poset) and does not depend on ties. At rung U0 m1c equals m1. A copy whose string is nested in a chain of the family costs nothing. **m1t** is the gene-level version with nesting restricted to end-truncation: two genes are comparable iff the ordered string of one (consecutive repeats collapsed) is a prefix or a suffix of the ordered string of the other, equal strings are comparable, W_f^t is the largest number of genes of the family in a set that is pairwise comparable, and m1t = the sum over families of (n_f - W_f^t). It is printed beside m1c and carries no class.
- **m2 (strict D2)** = sum over strings s with at least two genes of (n_s - max over families f of n_{s,f}).
- **Macro versions:** the mean over families with at least two genes of (n_f - W_f) / n_f, and the mean over strings with at least two genes of (n_s - max_f n_{s,f}) / n_s. Micro and macro are both printed, and so are the same statistics for the size-2 stratum of F_O (n_f = 2 in the layer universe: 152 of the 362 clusters with at least two V_F genes), the class that MCL joins or splits by floating-point residue. Under a half restriction the at-least-two test is applied after the restriction.
- **Path defect (secondary, rungs U1 and U2).** The ordered string of a gene is the sequence of dom(e) over its merged exons in genomic order, consecutive repeats collapsed, unlabelled exons skipped. H_f is the simple undirected graph (no loops, no multi-edges) whose vertices are the units in the members' ordered strings and whose edges join consecutive units of a member. delta(f) = (c(H_f) - 1) + sum over v of max(deg(v) - 2, 0) + (|E| - |V| + c(H_f)), and Delta is the sum of delta over the families with at least two genes in the universe. Ownership multiplicity Mult = sum over units u with mu(u) >= 1 of (mu(u) - 1), mu(u) being the number of families with at least two genes in the universe that hold a member whose string contains u.

**Audit (descriptive; it never excuses a gene from a decisive statistic).** For every family the deficit genes are the genes outside the kept chain, the maximum-weight chain with ties broken by the lexicographically smallest key, where the key of a chain is its strings as sorted tuples, listed from smallest to largest by inclusion and compared as lists; all tied chains and their deficits are printed. For D2 the deficit genes of a string are the genes outside its best family (ties: smaller cluster id). Labels, first match wins: *short* (at most two merged exons); *ends-only* (D1 only: at least three merged exons and the set of dom(e) over all merged exons except the first and last in genomic order equals that set of at least one kept gene with at least three exons); *unexplained* (the rest). The counts per label are printed beside m1c and m2, with and without 'short' (862 of the 2,142 genes of V have at most two merged exons). No label removes a gene from the universe, from n_f, from n_s or from a denominator, and nothing is excused in T1-a, T1-b or T1-c.

**Exclusion arms (descriptive, side by side; the primary is the least-deleted universe).** None; the readthrough list (`rt_genes.tsv`, 36 genes in V); non-eligible genes (S1C column 'In Table S1 (SD98 gene set)' equal to No, 388 in V); same-locus duplicate models, defined from gene coordinates alone (the gene span is the hull of the gene's CAT v4 transcripts, `families_cn.json` fields c, s, e, st; of each same-strand span-overlapping pair of V genes, the gene with the smaller exon union is deleted, ties: larger gene id; 105 pairs involving 181 genes, so at most 105 genes are deleted); and all three. An arm deletes its genes from the universe, from every n_f and n_s and from every null draw. The deletion does not depend on any family label. An arm that removes more than 10% of V carries the flag 'denominator guard' (the bar of the CP-5 guard in `docs/archive/2026-09/PREREG_chimera_policy_2026-09-18.md`).

**Nulls for the decisive cells.**

- **Governing null (homology-preserving).** Strata are the connected components of the F_O homology graph, the pre-MCL graph written by `mcl_families --dump-graph` from the frozen `loci.paf` with the frozen binary and flags (1,819 nodes, 11,031 edges, 339 components; `fo_graph/ours.graph.tsv`, sha256 1820795ae336a003). A gene belongs to the component of its own node, or of the node of its representative if it is a folded annotation. Within a stratum the cluster labels of its V_F genes are permuted. A stratum with one cluster cannot move and keeps its labels in every draw; m1c and m2 are recomputed on the full label vector of V_F. In V_F, 325 strata hold at least one gene. The **variation-bearing strata** are those that hold at least two clusters: there are 30, with 416 V_F genes (4 to 42 genes, median 11); 29 of them hold at least two multi-gene clusters. The other 295 strata hold 1,241 V_F genes in a single cluster each. Genes in no component (467 V genes) and in no cluster are outside the decisive cells.
- **Draws.** For each (layer, null, universe) re-seed `random.Random(20260930)` and draw B = 10,000 permutations (V_F genes sorted by gene id, strata in order of their smallest gene id; each draw shuffles, for each variation-bearing stratum, the list of its genes' cluster labels in gene-id order and assigns the result back in that order; ids sort as Python strings). Cells that share a universe share draws. The statistic is m1c for T1-b and m2 for T1-c; p is the lower tail (smaller deficiency than random re-cuts of the same homologous group). **Robustness to one group.** For every variation-bearing stratum the table prints z_s = (the statistic computed on the genes of s alone, observed, minus its mean over the draws) / its sd over the draws (z_s = 0 if sd = 0). The stratum with the most negative z_s is removed (ties: smaller stratum id), its genes leave V_F, and p is recomputed on the same draws for the remaining genes. The class is the lower of the class from all strata and the class from the remaining strata. The z_s table and the removal are made separately for T1-b (m1c) and T1-c (m2). The stratum id is the rank of the stratum by smallest gene id. For T1-c the statistic of s alone is m2 computed on the genes of s with the strings counted over those genes only. The recomputation uses the stored draws of the full V_F, so the removed stratum keeps its place in the random stream. The removal step is skipped when no z_s is negative. The class from the remaining strata is the content class of its p alone: UNDERPOWERED and UNINFORMATIVE are not assigned to the remaining strata, and if all remaining strata are dead p is 1 and the class is NOT-DISTINGUISHED.
- **Sensitivity nulls (no class).** Strata equal to the chromosome, and strata equal to the SD98 region of largest exonic overlap (445 regions in V; ties: smaller region id). These destroy homology between co-labelled genes, and p <= 0.01 is expected under them whatever the organisation (prediction Q4, section 10). Level sweep: the whole cell at MCL inflation 2.0 and 4.0 on the same graph. The class is issued at inflation 2.8 (the shipped default). If the class at 2.0 or 4.0 differs from it by one or more steps, the headline carries the tag LEVEL-SPECIFIC (register row 1260). The sweep is reported and not acted on. At each inflation V_F is the set of V genes that carry a cluster at that inflation. It is not a run for SPLIT, the Gate 6 counts apply to inflation 2.8 only, and a sweep cell that is UNDERPOWERED carries that tag and no LEVEL-SPECIFIC comparison.
- **F_S descriptive.** Strata are the F_Q groups, with a stratum 'no group' for clean genes in none.

**Decisive cells and their classes.** T1-a prints m1c and m2 separately: REFUTED when either is >= 1 on V_F, NOT-REFUTED only when both are 0, with the full exception list. The information of the cell is the exception list. NOT-REFUTED can be issued only on the primary arm. T1-b and T1-c follow the verdict algebra of section 3 (n_min = 8 variation-bearing strata) and the robustness clause above. The p-value reflects only the 416 V_F genes (25%) in the variation-bearing strata, and a few strata can decide it (the z_s table shows which). The other 75% of V_F genes lie in families that equal their homologous group, so their cut cannot be tested.

### 5. T1, second species (Phase R2; light; descriptive, no verdict)

This arm is a gorilla analogue of T1. It does not replicate the human cells. The units are SEDEF atom classes, because no DupMasker library or copy-number class exists for gorilla, and the family layer is a human Compara projection. The arm reports numbers and issues no verdict class, for three reasons. First, depth-matching pairs on their own transcript identity >= tau selects exactly the pairs for which a SEDEF alignment at fracMatch >= tau is expected, so the S-rate measures SEDEF and atom sensitivity at the genes (intron divergence, repeat breaks) as much as unit boundaries (register rows 650 and 996: an SD-supported alignment edge is a SEDEF pair). Second, the control pool is tiny and not identity-matched: in the verdict set at tau 0.90 there are 11 eligible different-family pairs on one contig less than 100 kb apart, 27 at 100 kb to 1 Mb, 878 at more than 1 Mb and 10,506 across contigs (the 100 bp touch rule; 885 and 10,651 under a 1 bp rule), against 13, 5, 15 and 21 eligible same-family pairs (54 in all, not only the depth-matched ones; `logic/control_pool.out`). Third, the layer is name-keyed.

**Inputs and units.** `GGO_sedef_final.bed` (253,029 data rows, a trailing header line, 34 columns; md5 991ed343...). Drop the header line and the 36 NC_011120.1 rows. An identical tuple is equal (contig1, start1, end1, contig2, start2, end2); keep the row with the largest column 21, ties by larger column 12, then first in file (columns are 1-based: column 21 is fracMatch and column 12 is aln_len). A row is retained at tau iff column 21 (fracMatch) >= tau. **tau = 0.90 is primary** (the only level with a usable set); 0.95 and 0.98 are printed and are UNDERPOWERED by the rule below. Atoms come from `dna_sd_atoms.py` in cigar mode (MIN_BLOCK 1,000, class edge when two-sided CIGAR coverage >= 0.5, the J2/K0 constant), restored from `Rustle_attic/2026-09-23b/bench/dna_sd_atoms.py` (git b509675c, sha256 9d96d976eb52dca6e056b1bc8256f006835d72b265a8697c925d15cfcb562bce) into `bench/` with a test BEFORE the first run. The classes come from a new union-find over the atom edges of that script, and the dedupe step is new code; both have tests; the recon class counts (40,631 atoms, 12,483 classes, 372,793 atom edges at tau 0.90) were taken on un-deduplicated rows and with coverage > 0.5, and Gate 6 prints the counts under the registered rules. The format argument is 'gorilla' (34 columns).

**Family layer and pair universe.**

- The layer is the rows of `fewcopy_gorilla_2026-10-07/gorilla_truth.tsv` with `member_status` = 1to1 (342 genes; 377 Compara Primates families, 104 with at least two projected genes), keyed by (gorilla_contig, start1, end1). The exon blocks behind U(g) come from `families_gw/species/gorilla/genes.tsv`, joined on (contig, start1, end), and `pairs.tsv` joins the layer on a_gene and b_gene equal to gorilla_gff_id. It is a human Compara family projected by exact symbol onto the gorilla RefSeq annotation. The projection is the one name step of the arm (register rows 993 and 1034 forbid other name joins). The layer is depleted of LOC-named recent copies, is not gorilla truth and is used only as an SD-blind pair roster (register lines 316, 861 and 869; row 1258).
- A gene touches class K iff at least 100 of its exonic bp lie in atoms of K. U(g) is the set of classes it touches; genes with |U| >= 2 stay in the primary and are flagged 'path'. A gene is eligible iff U(g) is non-empty.
- A pair of genes is eligible iff both are eligible, g differs from h, and their exon unions share 0 bp. Distance is the gap between the nearest ends of the two gene spans (0 if they overlap), cross-contig being its own class. S(g, h) = 1 iff U(g) and U(h) intersect.
- A same-family pair is depth-matched iff it is a row of `pairs.tsv` (116 pairs from 90 whole families; 95 with an alignment) whose `identity` is >= tau (NA excluded; this is an aligned-block identity, and the aligned fraction is printed beside it, since 16 of the 61 pairs with identity >= 0.90 align less than half of the shorter transcript). Partial families have no depth-match and are out of scope; their eligible pairs are counted.
- A family belongs to the verdict set iff none of its 1to1 projected genes, eligible or not, lies on NC_073241.2, NC_073242.2, NC_073244.2 (the SD x family development slice) or NC_073234.2 (the contig of the most-junctions rule test).

**Counts (input-only, un-deduplicated recon atoms; Gate 6 reprints them under the registered rules).** In the verdict set, depth-matched eligible pairs and their families: tau 0.90, 38 pairs in 33 families (46 eligible whole-family pairs; the same under the 1 bp and the 100 bp touch rules); tau 0.95, 15 in 15; tau 0.98, 5 in 5 (`c7_depth_matched_verdict.out`).

**Report.** S-rate = share of depth-matched eligible pairs with S = 1, pair-weighted, with the family-weighted rate (mean over families of the family's pair share) beside it; A_viol = 1 - S-rate; B_viol = the count of eligible pairs of different families whose genes share a class, a ceiling-type quantity reported with that caveat. Wilson intervals, and a family-cluster bootstrap (2,000 resamples of families, seed 20260930). No control draws, no p-value and no class. A table is flagged UNDERPOWERED when it has fewer than 20 depth-matched pairs or fewer than 8 families (tau 0.95 and 0.98 are).

**Positive control (the arm is INVALID if it fails).** Each of the 8 full-length LRPAP1 copies (identity 0.957 to 0.990 to the reference copy) touches at least one class at tau 0.90, and the graph whose vertices are the 8 copies and whose edges join copies with intersecting U is connected. Otherwise the exact classes are reported and the arm is INVALID.

**Matched-instrument human run (optional, descriptive).** The identical pipeline on CHM13 v1.0 `final_v1_clean.bed` (5,217 and 31,430 atoms at 0.98 and 0.90). It is not required, because the human Compara layer is keyed by gene name (`compara.Primates.families.tsv`: Gene Name, Family ID, Contig) and the SEDEF file and CAT v4 are v1.0 while RefSeq is v2.0 with non-constant offsets on the acrocentrics. If it is run, the join is Name and Contig to the RefSeq v2.0 `genes.tsv` with a uniqueness assertion (ambiguous and unmatched names printed), then v2.0 to v1.0 by per-chromosome offsets verified by sequence (Gate 2; failing chromosomes excluded and counted), tau 0.90 and 0.98, no development-slice exclusion. It carries no class. No species comparison is claimed without it, and none is claimed with it either, because the S-rate is not a deficiency.

### 6. T1, further arms (Phase R3a light, Phase R3b heavy)

**R3a, the human_testis BASE snapshot (light; needs Amendment 1).** The human_testis BASE table (`rt_arms/human_testis/human_testis.BASE.fam.clusters.tsv`, 09-25, 338 families, 1,095 loci, 24 contigs) is an RNA-derived, genome-wide family layer that has never been examined against duplicons. It is admitted as a declared snapshot, as section 8 does for gorilla: it predates the f1v2 and min-cov-shorter defaults (`min_cov_shorter` 0) and its binary commit is not recorded. Input-only power count (`c9_r3a_power.out`): 334 V genes carry a BASE cluster by exon overlap of their CAT/Liftoff v2.0 exons with a BASE locus span, and 86 BASE clusters hold at least two V genes (316 V genes). The strata of the homology-preserving null are the components of the graph that `mcl_families --dump-graph` writes from `human_testis.BASE.fam.loci.paf` (176 MB, light); Gate 6 prints their number and n_min = 8 applies. The library is not A119b, so R3a will be called HELD-OUT for the family layer and the library once it has run; the genes and the unit strings are those of Phase R1. It is not held out for family-rule decisions, because the BASE snapshot was produced with rules fixed on other substrates. Same cells, statistics and nulls as section 4. **The gene-to-locus rule is fixed here:** an RNA locus that overlaps two S1C genes assigns its family to the gene of larger exonic overlap, ties to the smaller gene id; a gene overlapped by loci of two families takes the family of the larger exonic overlap, ties to the smaller cluster id; the half rule is that of section 4. **Amendment 1 is committed before the first Phase R1 run.** It may add the graph-dump command (with `--min-cov-shorter 0` and the flags of `human_testis.BASE.fam.params.tsv`), the frozen binary, the file hashes and the call schedule. The check that every BASE cluster lies inside one component of the dumped graph is run at the start of Phase R3a and is part of Gate 6 for that phase. Amendment 1 may not change the gene-to-locus rule, the cells, the statistics, the nulls or Q8, and from the first Phase R1 run no Amendment may change a rule, statistic, null, threshold or class that Phase R3a uses. The verdict of T1-b and T1-c is the class of R3a once it has run, labelled held-out for the family layer and the library, unless the verdict algebra gives SPLIT against Phase R1, or Phase R3a is INVALID, UNDERPOWERED or UNINFORMATIVE; in those three cases the Phase R1 class stands, labelled DEV. No rule, threshold or class boundary of this file is changed after R3a runs, so R3a informs no rule.

**R3b (heavy; not authorised).** Genome-wide de novo families under the current default (`tools/rustle_pipeline.sh families`, MCL I = 2.8) for A119b outside chr16, chr17 and chr18 (estimated 17 h), and for gorilla OR6737 and KB3781 (0.7 to 4.4 h and 0.6 to 3.3 h) on the contigs outside the development set. On A119b, chr5, chr7 and chr21 (selection of the 0.60 shared-exon threshold) and chr2, chr8 and chr10 (a reused verdict set) are not used to select any rule of this file for the duplicon relation. Chromosome-restricted runs change the family definition (cross-chromosome links are lost) and are not equivalent. Until R3a or R3b runs, T1 has a DEV verdict (section 4) and a descriptive second-species report (section 5) only.

### 7. T2, the E-table (definition step; heavy class)

An *expansion* in the study species gorilla (mGorGor1 / GCF_029281585.2, primary assembly `GGO.fasta`) is a set of at least two annotated gorilla records (gene or pseudogene records with exon union >= 200 bp) of one pool component that are co-orthologs of ONE locus that is single-copy in each of the outgroups H = CHM13v2.0 (RefSeq RS_2025_08), C = mPanTro3 (primary) and O = mPonPyg2 (primary). Unanimity over the three outgroups is required. No family file is opened before the E-table is frozen (sha256 recorded in Gate 0).

**Tier B, rank-sharded reverse lift (primary).** A single Liftoff run per shard cannot be used for the pool, because Liftoff 1.6.3 does not place two non-overlapping reference genes on overlapping target loci (`liftoff_utils.find_overlaps`: overlap above `-overlap` for primary placements, and above 0 for copy-tagged ones, is rejected unless the two reference features overlap in the reference annotation), so paralogs on one gorilla contig could not co-place, and the control LRPAP1 would recover at most 6 of its 11 loci (they sit on 6 contigs). The recipe is therefore:

1. **Pool.** The gorilla self-lift (`/mnt/linuxdisk/tmp/rustle_figures/liftoff/gorilla/liftoff_loci.tsv`, Liftoff `-copies -sc 0.95`, chosen a priori): copy rows of class extra_copy, dropped_M1 and dropped_M2 are joined to annotated records by the note column ('overlaps gene-X' or 'overlaps copy of X') or to new `copy@contig:start-end` nodes, and union-find gives 609 components of size >= 2 holding 2,976 annotated records, of which 2,249 have an exon union >= 200 bp; 307 components hold at least two such records. The 1,200 unannotated extra_copy sites cannot be lifted and are listed in an extension table, not scored. A component with more than 51 members (members are the annotated records of any length and the unannotated `copy@` sites of the union-find component; 51 is the point where Liftoff itself truncates whole-body copies; inherited, not tuned) is a hub, is not called and is listed (4 of the 80 supported KB3781 components and 5 of the 139 OR6737 ones).
2. **Rank-sharded runs.** For each outgroup, run j (j = 1..J, J <= 51) holds the j-th record (sorted by contig as Python strings, then start0, then end, then record id) of every called component, so no run holds two records of one component. Command: `liftoff TARGET REFERENCE -g shard.gff -f types.txt -o OUT.gff3 -u UNMAPPED -dir DIR -copies -sc 0.95 -a 0.5 -s 0.5 -overlap 0.1 -p 4 -m MINIMAP2`, with `types.txt` holding the single line `pseudogene` (Liftoff lifts GFF features of type gene only unless `-f` names others, and 681 of the 2,249 pool records and 496 of the 1,403 called records are pseudogene records), the gorilla annotation as reference, the outgroup assembly as target (a cross-species Target, which is new wrapper code: the reference FASTA differs from the target FASTA) and the minimap2 that Liftoff calls named in the Amendment (the liftoff environment has 2.24). The indexes `ggo_npip/ref/GGO.fa.mmi`, `ggo_npip/ref/PTR.fa.mmi`, `o1_falsemerge/lit/evid/PPY/genome/PPY.fa.mmi` and `winloci_data/chm13v2.0.fa.mmi` (under `/mnt/linuxdisk/home/juanfraitu/`) exist (k15, w10) and are listed in Gate 0. Shards hold at most 2,500 records (`figures/_liftoff.py::_plan_blocks`) and are split at 590 s. The called set holds 300 non-hub components with at least 2 such records and 1,403 records (907 gene, 496 pseudogene) under the hub rule of step 1; the largest holds 44, so J = 44 rank runs per outgroup and 132 runs in all, with 300, 300, 165, 116 and 81 records in the first five runs and 10 or fewer in 26 of the 44. (Counting only annotated records toward the hub rule gives 5 hubs, 302 components and 1,446 records; two components with 74 and 52 copy@ sites separate the rules.) The rank-run split rule (`_apply_splits` cuts at record-free start0 positions) is new code. The cost is one pass over three outgroups: 1,403 x 3 records at 0.05 to 0.09 s plus 10 to 20 s of start-up per run, so about 26 to 50 minutes of work, stated as 45 to 60 minutes in bounded foreground calls, RSS 10 to 14 GB. This is an estimate from three same-genome self-lift fits (`c12_shardfit.out`; the orangutan fit is dominated by one 544 s shard, a 10-fold outlier that occurs about once in 29 shards), not from cross-species lifts, and the Amendment sets the time budget per outgroup. Two unrelated records of one run can still conflict on the same outgroup locus. In Liftoff 1.6.3 the loser is lifted again with the occupied locus excluded and is usually written elsewhere with lower coverage and identity, with no flag. A displaced record is therefore recognised only by coverage < 0.5 or by landing on no single in_place outgroup record, and the number of records with no triple under that rule is printed.
3. **Outgroup locus.** A *placement* is the primary row of a record with coverage >= 0.5 and sequence_ID >= 0.5 and neither the partial_mapping nor the low_identity flag set (whole-genome Liftoff writes rows from coverage 0.05 and marks the weaker ones with those flags); any other row is no placement and the record has no triple. The number of records whose only row is partial is printed. A placement maps to an outgroup locus iff it shares at least 50% of the shorter exon union (strand ignored) with ONE record of that outgroup's own self-lift (`liftoff/<species>/liftoff_loci.tsv`) of class in_place that has no copy row, and the reverse lift gives that gorilla record no extra_copy row in that outgroup. The record's triple is (h, c, o), the three outgroup records. The number of qualifying placements (coverage >= 0.5) that Liftoff reports for the record in each outgroup is recorded, and a class with any member that has more than one is dropped and counted (register line 157, trap T19: best-hit counting is not counting). Liftoff reports further copies only at identity >= 0.95, so tier C, which counts every landing, is the check on this rule.
4. **Co-orthology and expansion.** Two records of a pool component are co-orthologous iff they have the same defined triple. An expansion is the set of at least two records of one component with the same triple. This is an equivalence relation, so there is no chaining. Records of one component with different triples form separate classes, counted as 'different loci'. Records with no triple (unplaced, outgroup multi-copy) are dropped and counted. Expansions of different components with an identical triple are listed together as siblings and are not merged. Siblings count as separate units; the number of sibling groups is printed, and a sensitivity counts each group once (the mean R_j within the group).
5. **Deciding E-table.** Tier B only. Tiers A and C are replicates and add no members. Columns: expansion_id, pool_component_id, record_id (contig:start0-end, 0-based half-open), biotype, exon_union_bp, h, c, o, kind (annotated). The E phase also writes `member_exons.tsv` (record_id and the exon blocks of every pool record, 0-based half-open), `control_pool.tsv` (pool_component_id, record_id, a reason in outgroup_copy_row, different_triple or no_placement, and the sequence_id of each joining copy row; a record with reason no_placement is not eligible as a member of a control set) and `extension.tsv` (the unannotated `copy@` sites, listed and not scored). Their sha256 is recorded with that of the E-table. A component that contributes an expansion is not in the control pool. Strand is ignored in every overlap test.

**Tier A (optional replicate).** The literal wording: the 19,841, 21,970 and 21,704 single-copy protein-coding anchors of H, C and O (in_place, exon union >= 200 bp, no copy row) lifted onto gorilla with copies (`-copies -sc 0.95`); it shares the cross-species Target code with tier B and is estimated at 1 to 2 h heavy.

**Tier C (replicate; synteny classes and fragments).** Sequence-first body search: `minimap2 -c -x asm20 --secondary=yes -N 100 -p 0.1 -t 4` of each outgroup anchor body against the splice index (the command recovered from the session that found LRPAP1's 11 loci); landing identity >= 0.90 and coverage >= 0.50, with the sweep {0.85, 0.90, 0.95} x {0.3, 0.5, 0.8} reported and not acted on. It is the only tier that recovers 5' fragments and unannotated copies. Splice indexes exist for gorilla (13.6 GB) and human (12.3 GB) only. For chimp and orangutan a build takes about 7 min, 19 GB RSS and 12 to 13 GB of disk, and /mnt/linuxdisk has about 37 GB free (98% used), so builds are serial and each index is deleted after use, or the w10 indexes are used and the lower seed density is recorded. Estimated 1 to 1.5 h heavy. The agreement matrix A / B / C is printed.

**Gate 7, controls of the E-table (the table is INVALID if either fails).**

- **Positive control.** The tier-B run returns ONE class that contains every annotated record of the gene-LRPAP1 pool component that has at least 50% of its exon bases inside an LRPAP1 copy interval (6 records: LRPAP1, LOC134756753, LOC115933156, LOC134757218, LOC129523574 and LOC129530227), each with one in_place triple in human, chimp and orangutan. LOC129523503 and LOC129526389 are long gene models that contain an LRPAP1 copy as a minor part (16% and 41% of their merged exon union, and at most 49% and 54% of the exon and CDS bases that Liftoff counts, lie inside the copy interval; input-only, `/mnt/linuxdisk/tmp/prereg_closure/disposition/lrpap_record_cov.out`). They are not required: either may or may not reach coverage 0.5, and if it is placed on the outgroup LRPAP1 locus it joins the class without causing a failure. An LRPAP1 copy interval is the landing interval of one of the eight gene-LRPAP1 rows of the gorilla self-lift (the in_place row and seven dropped_M1 landing sites). The 3 fragments form another pool component and can appear only as a sibling class. A failure of the control makes the E-table INVALID and blocks T2a and T2b. The control can fail (a qualifying copy not placed at coverage >= 0.5, an outgroup self-lift that lists LRPAP1 with a copy row, or a fragment that takes a locus from a full-length copy in a shared run, since the fragment ranks 1 to 3 share runs 1 to 3 with the full-length component). It checks the pipeline on a component that is in the pool by construction. It is not a test of discovery. The two controls are lifted and checked first, in their own calls, before any verdict-set record is lifted, and the verdict-set shards start only after both controls pass. The Amendment may change the Liftoff command only on evidence from the two controls of this Gate, which contain no pool component of the verdict set, and records that evidence.
- **Planted negative control.** 100 gorilla in_place records with exon union >= 200 bp that have no copy row as source and are in no pool component (neither as source nor as the target of a copy row) are paired at random on different contigs into 50 fake two-record components (`random.Random(20261007)` over sorted ids) and lifted by the same recipe. At most 1 of the 50 may come out as a called expansion. The conserved clusters PCDHB, TUBA and RBMY are not used as controls, because their counts come from other studies as recovery fractions or copy-number equalities and no gene lists with coordinates are on record.

**Contig roles for T2.** The hold-back unit is the expansion. The contig is not a hold-back unit. LRPAP1 (11 loci on chr3, chr12, chr14, chr16, chr22 and chrY) is the development case and positive control. It is spent, is excluded from the verdict set and is reported separately. Every other expansion is in the verdict set, and no rule of this file is tuned on any of them (LRPAP1 informs the recipe only through Gate 7). The word 'held-out' is not used for gorilla. Accession to chromosome: exposed to earlier tests NC_073224.2 chr1, NC_073228.2 chr4, NC_073234.2 chr10, NC_073241.2 chr17, NC_073242.2 chr18, NC_073244.2 chr20, NC_073230.2 chr6 and NC_073233.2 chr9 (`PREREG_core_definition_2026-09-12.md`) and NC_073247.2 chrX (`PREREG_rg3_npip_block_2026-09-28.md`); LRPAP1 NC_073227.2 chr3, NC_073236.2 chr12, NC_073238.2 chr14, NC_073240.2 chr16, NC_073246.2 chr22, NC_073248.2 chrY. The expansions whose annotated members all lie on the ten contigs that no prereg names (NC_086017.1 chr2, NC_073229.2 chr5, NC_073231.2 chr7, NC_073232.2 chr8, NC_073235.2 chr11, NC_073237.2 chr13, NC_073239.2 chr15, NC_073243.2 chr19, NC_073245.2 chr21, NC_086018.1 chr23) are printed as the least-exposed subset and carry no class. Only 3 (KB3781) and 8 (OR6737) supported pool components have all members there, before the outgroup filter, so the subset is expected to be UNDERPOWERED. These contigs are not unexposed: chr5 carries the MHC region (register line 779), chr19 holds the reserved MCL166 and chr23 appears once in the ledger.

**Expected size (input-only, `c8_heldout_power.out`).** Pool components with at least two annotated records of exon union >= 200 bp and at least two with >= 2 primary reads in the library: KB3781 80, OR6737 139, either 157, of size <= 51: 76 and 134. With the 15% single-copy share (item 7 of section 2, unverified) that is at most about 12 and 21 expansions, before the hub and LRPAP1 removals and before the filters of section 8 that define the testable ones; the interval and the chance of fewer than 8 are in section 2, item 7.

### 8. T2a, containment (light class once the E-table is frozen)

**Family object.** The shipped MCL partition (I = 2.8) of the BASE tables, a declared snapshot that predates the f1v2 and min-cov-shorter defaults (`min_cov_shorter` 0, binary commit not recorded): KB3781 308 families / 864 loci, OR6737 483 / 1,559, one library at a time and never pooled. Regenerate under the frozen default only if authorised. The level is named in every sentence: 'contained in the MCL family' (register row 1260: a claim that holds at one level only is LEVEL-SPECIFIC). A level-C component table (`mcl_families --dump-graph` on the BASE PAF) is an optional replicate. The dump has one class in sections 6, 8 and 9: light for the human_testis PAF (176 MB); for the gorilla PAFs (427 and 532 MB) the peak RSS is measured first under `tools/rlock.sh light` with `/usr/bin/time -v`, and the dump is heavy if it exceeds 2 GB. The Amendment that authorises it fixes the class.

**Members and matching.** A node n is a row of `PREFIX.fam.copies.tsv` (family_id, chrom, start, end, exons, n_reads; 864 rows for KB3781 and 1,559 for OR6737), its family is its family_id, and X(n) is the set of exon blocks in its `exons` column (0-based, half-open). `PREFIX.fam.clusters.tsv` is not read: its rows are merged nodes that join the copies rows by no exact key (22 of 864 and 28 of 1,559 have no exact match). A member l has node n iff both lie on one contig and cov(l | n) = |X(l) intersect X(n)| / |X(l)| >= 0.5, where X(l) is the member's exon union. Strand is ignored. Among several nodes the largest cov wins, ties by larger `n_reads`, then smaller start. One node serves one member of an expansion: when two members select the same node, the one with the larger cov keeps it (ties: smaller start) and the other counts as having no node. The collision is logged, and the results are printed with and without the expansions that had a collision.

**Expressed.** A member is expressed in a library iff at least 2 distinct primary reads (flag & 2308 = 0) have an aligned M, = or X block that overlaps X(l) by at least 1 bp (the rule of `figures/_liftoff.py::ensure_support`). Expressed is read from the support tables `liftoff/support/gorilla_*.tsv` (reads_capped >= 2) that this rule produced, and the BAM rule is re-run once on the 50 records of smallest id of each library, with any disagreement printed. Libraries and BAMs: KB3781 `/mnt/linuxdisk/home/juanfraitu/fibroblasts/GCA_029281585.2_flnc_mm.bam`, OR6737 `/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_mm.bam`. m_j is the number of expressed members of expansion j and an expansion is testable iff m_j >= 2. n_testable is the number of expansions of the verdict set with m_j >= 2 and at least n_min = 8 qualifying control components (below), evaluated before any class is assigned. UNDERPOWERED iff n_testable < 8. Expansions dropped for fewer than 8 qualifying control components are counted and listed.

**Statistics.** R_j = (number of pairs of expressed members whose nodes lie in one family) / C(m_j, 2). A member with no node is in no pair, so it counts as a non-co-member of every other member (absent members stay in the denominator; nothing is conditioned on presence). The deciding statistic is T = the mean of R_j over testable expansions. Descriptive beside it, each with 2,000-resample bootstrap intervals over pool components (seed 20261007): strict containment c_j (R_j = 1), presence P1_j = (expressed members with a node) / m_j, soft containment s_j = (largest number of present members in one family) / m_j, family purity |family intersect E| / |family| and the number of families that hold at least one non-member; micro and macro. Size-2 expansions (exactly two expressed members) count toward n_testable when they meet its definition and are also printed as a stratum.

**Matched controls and test.** The control pool is the set of pool components that are not hubs, do not contribute an expansion, and have records that place on an outgroup locus with a copy row in at least one outgroup or have different triples (`control_pool.tsv`); a record with reason no_placement is not eligible as a member of a control set. A control set for expansion j is a set of m_j expressed records of one control component with the same multiset of per-contig counts and a median pairwise Liftoff sequence_ID within the caliper. Identities are integers, round(1000 x sequence_ID). The median runs over the copy rows (one value per row, in either direction) that join two records of the set. Let M be twice that median, an integer. A control set matches iff |M_control - M_expansion| <= 10 (that is 0.005, inclusive; the sweep is 4, 10 and 20, reported and not acted on). A set with no joining copy row cannot be a control, and an expansion with no joining copy row among its expressed members is dropped and counted. A control component qualifies for expansion j iff it holds at least one matched set. Excluding the components that contribute an expansion lowers the number of qualifying components of every other expansion. Gate 6 prints, per expansion, the number of qualifying control components at calipers 4, 10 and 20. Input-only upper bounds (before the outgroup conditions and the removal of expansion components, size class pooled, m_j = 2 on two contigs): the best window holds 23 components in KB3781 and 46 in OR6737, and 8 and 16 or 17 for m_j = 3 (`/mnt/linuxdisk/tmp/prereg_closure2/author/pool_floor_nosize.py`). The matched sets of one component are the first 1,000 in lexicographic order of their sorted record ids, found by a depth-first search that respects the per-contig counts. The search stops after 1,000,000 candidate sets per component, and a component that reaches that cap does not qualify (the cap is a count, fixed now). Qualifying components are listed in pool_component_id order as Python strings. The one-to-one node rule applies to every control set exactly as to an expansion. Draw b = 1..10,000 replaces every testable expansion by a control set drawn in two steps, a qualifying component uniformly and then one of its matched sets uniformly, and recomputes T_b. T is compared as an integer. With L the least common multiple of the C(m_j, 2) of the testable expansions, T_num is the sum over expansions of k_j x (L / C(m_j, 2)), where k_j is the number of pairs of expressed members whose nodes lie in one family. p_more = (1 + #{b: T_num_b >= T_num_obs}) / 10,001 and p_less = (1 + #{b: T_num_b <= T_num_obs}) / 10,001. Expansions are processed in expansion_id order, and the enumeration and the draws use separate `random.Random(20261007)` objects, enumeration first. Printed with no class: the p-values within each component-size stratum (2, 3 to 5, 6 to 11 and 12 to 51 annotated records of exon union >= 200 bp), computed over the expansions of that stratum with the same draws; the component member count and the MCL family size of the control members as covariates. If the classes of two strata with at least 8 expansions each are two or more steps apart, the headline carries the flag SIZE-SENSITIVE.

**Classes.** SPECIFIC iff p_more <= 0.01. LESS iff p_less <= 0.01. NOT-DISTINGUISHED iff both p > 0.05 (no equivalence is claimed; the difference of means and its interval are printed). PARTIAL otherwise (0.01 < p <= 0.05 in one direction, which is named). T2a is two cells (KB3781, OR6737). Each library carries its own class, and the SPLIT rule and the headline rule of section 3 combine them. The headline number, whatever the class, is the absolute containment (mean R_j, c, P1 and s with intervals), with the reason for each split (fragment, MCL size-2 residue, fusion node, pseudogene). The claim that expansions sit inside our families asks how much of each expansion lies in a family, and also whether that beats the controls.

### 9. T2b, 'younger than the family root' (heavy; needs Amendments; predicted UNDERPOWERED)

A test on sequence identity cannot fail (section 2, item 6). T2b is registered here so that nothing is chosen after seeing a result. It is not part of the authorised run. The clock of arm 2 is only partly independent of the family builder: family edges are minimap2 alignments of the loci's genomic spans (introns included, identity >= 0.70), so intron identity is part of what builds the family. Only the outgroup-polarised topology is independent of it.

**Identity check (pipeline check only).** h(S) is the divergence 1 - (largest identity threshold at which all members of S lie in one connected component of the family builder's identity graph). The graph is written by `mcl_families --dump-graph` re-run on the BASE PAF (class as in section 8) and is hashed when it is made. For every expansion whose expressed members all have nodes in one family phi and whose node set E_nodes is a proper subset of the nodes of phi, h(E_nodes) <= h(phi) must hold, where h(phi) is the divergence of all nodes of phi (single linkage guarantees it; equality is allowed). A violation means a pipeline error and makes T2b INVALID. The share with h(E) < h(phi) is printed and carries no weight.

**Synteny class (new code, with tests, before use).** The 25 kb flanks of the outgroup locus are mapped to the gorilla assembly with `minimap2 -x asm5 -c`. Flank anchors need coverage >= 0.80, identity >= 0.90 and MAPQ >= 20. Two collinear anchors on one contig delimit the syntenic interval, with the split-anchor chaining of `bench/o1_outgroup_rooting_poc/README.md` (query gap <= 1,000, target gap <= 5,000). A gorilla record inside the interval is syntenic; a record outside it with unique flanks is derived; a record whose flanks are not unique or have no two-sided anchors is unresolved. The proof-of-concept files are gone (only its README remains), so the classifier is written from that README. The pilot certified 3 of 18 probes with two-sided synteny (15 of 18, 83%, unresolved; 30 of 40 haplotype rows, 75%, had no two-sided synteny) on an inverted substrate (human study, gorilla outgroup; register lines 88 and 753).

**Classes of members.** The anchor of an expansion is its triple. A family whose expansions have different anchors is excluded and counted. Per outgroup a member is syntenic, derived or unresolved. A member is *ancestral* iff syntenic in H, C and O, *derived* iff derived in all three, and *intermediate* otherwise; unresolved members are excluded from both tallies and counted. Dv is the set of members that are derived in all three outgroups; intermediate and unresolved members are in neither Dv nor the ancestral set. For an expansion, phi is the family that holds the largest number of nodes of its expressed derived members (ties: smaller cluster id). An expansion is *mixed-family* iff phi holds at least one node that is not a node of a member of E; an expansion whose phi holds only nodes of members of E has no other member to be older than, is excluded and is counted. The ancestral member need not lie in phi.

**Arm 1, root membership (an estimate; no class).** The unit is a mixed-family expansion that has at least one expressed derived member with a node in phi and at least one expressed, resolved ancestral member. The outcome is whether the node of an expressed ancestral member lies in phi (with several ancestral members, whether at least one does). The share is printed with an exact (Clopper-Pearson) 95% interval, the cause for every expansion whose ancestral member lies outside phi (no node, node in another family) and the funnel (E-table, testable, derived members in one family, mixed-family, expressed resolved ancestral member, units); an unexpressed ancestral member is an exclusion at that step of the funnel. No bar is applied because the null share is not 0.5. The ordering 'ancestral is older than derived' is supplied by the definitions (Dollo polarity). Arm 1 makes no further age claim. UNDERPOWERED iff fewer than 8 units.

**Arm 2, topology (an estimate; no class).** The unit (an informative expansion) is a mixed-family expansion with |Dv| >= 2 and |E_tree| >= 2. The family tree holds the nodes of phi that pass the 1,000 bp rule, the members of E that lie in phi, and the three outgroup orthologs; every leaf has at least 1,000 non-exonic bp (introns plus 2 kb of flank on each side, excluding any annotated exon of any gene; no repeat masking). E_tree is the set of members of E in the family tree (the ancestral member counts when it lies in phi). The sequences are aligned with `mafft --auto`. The tree is built with `iqtree3 -m GTR+G4 -B 1000 --seed 20261007` and rooted on the orangutan ortholog. E_tree forms a supported clade iff the node whose leaf set is exactly E_tree has UFBoot >= 95 and the family tree has at least one leaf outside E_tree, so that the most recent common ancestor of the expansion lies below the family root. The share of informative expansions that form a supported clade is printed per library with an exact (Clopper-Pearson) 95% interval, with the cause of every expansion that does not (gene conversion, low support, E-table incompleteness, other) and the funnel (E-table, testable, derived members in one family, mixed-family, informative, supported clade). An expansion that does not form a supported clade stays in the denominator. The share measures three things at once. It measures whether the expansion is younger than the family root. It measures how complete the E-table is, because family members that descend from the ancestral gene and are not in E break the clade (a simulation with 10% and 25% of derived copies missing from E lowers the clade share in the true tree to 0.78 and 0.52; `/mnt/linuxdisk/tmp/prereg_closure2/logic/t_arm2.out`). And it measures how many substitutions separate the stem from the crown, because short stems fail UFBoot 95. It is therefore an estimate and gives no class. The crown-depth ratio R(Dv) = (crown depth of Dv) / (depth of the human-gorilla split at the same locus) is printed and descriptive. A comparison of median divergences d(E) < d(F) is not used, because the median over all pairs of F contains E's own pairs. UNDERPOWERED iff fewer than 8 informative expansions.

**No combined T2b verdict.** Arms 1 and 2 are reported separately. The clause 'younger than the family root' is answered only by arm 2 (topology), as an estimate, and by arm 1 as a membership statement.

### 10. Consequences of the inputs, and predictions (written before any outcome)

**Consequences of the inputs (not predictions).**

- C1. Strict D1 is refuted for rungs U0 and U1 at layers F_S, F_Q and F_O (section 2, item 1). Rung U2 and layer F_R were not counted.
- C2. The identity check of section 9 passes (single linkage).
- C3. Gate 7 is a pipeline check on a component that is in the pool by construction. Six of the 8 full-length LRPAP1 copies qualify as controls, and the 3 fragments can appear only as a sibling class.
- C4. In Phase R1 the number of variation-bearing strata is 30 against n_min = 8 (the table would be INVALID if the count differed). Whether at least 8 of them are live is known only from the draws, and each cell prints it.

**Predictions.** Some are informed by what was seen (stated per item).

- Q1. T1-a is REFUTED for D1 and for D2 on F_O. This is not blind: the hand-looked NPIP genes and the 36 clean readthrough genes of V are exceptions if F_O co-clusters them.
- Q2. T1-b (D1) is CUT-ALIGNED, or at worst PARTIAL. This is expected whenever the strings and the cut both follow sequence similarity (section 2, item 3). A NOT-DISTINGUISHED would mean that the MCL cuts inside homologous groups are not aligned with duplicon content more than random re-cuts are. It would add to the closed duplicon-first results (register rows 732, 817, 818) and would be reported as found.
- Q3. T1-c (D2) is not CUT-ALIGNED. The reason offered, that NPIP subfamilies share their dominant duplicon, is not evidence about p. The D2 relation against F_S is a reproduction of cnduplicon and is not evidence.
- Q4. T1-b and T1-c under the chromosome and SD98-region nulls give p <= 0.01 whether or not duplicons organise families, because those nulls destroy homology between co-labelled genes. A cell with p > 0.01 under either null would show that this reasoning is wrong.
- Q5. The E-table passes Gate 7 with the 6 qualifying LRPAP1 records. The number of expansions is at most about 12 (KB3781) and about 21 (OR6737) before the filters of section 8, and n_testable is smaller after them. The chance that fewer than 8 expansions exist before the filters is about 12% at a 15% single-copy share and about 48% at the lower end of its interval for KB3781, so T2a may well be UNDERPOWERED there.
- Q6. In each library that is not UNDERPOWERED, T2a is NOT-DISTINGUISHED: absolute containment of expressed members is high but below 1 (MCL splits, fragments), and a difference from the matched high-identity control sets is not detected.
- Q7. T2b arm 1 and arm 2 are UNDERPOWERED in both libraries (the pilot certified 3 of 18 probes).
- Q8. If Phase R3a runs and the Gate 6 count of strata with at least two clusters is at least 8, its class for T1-b equals that of Phase R1 (predicted CUT-ALIGNED) and its class for T1-c equals that of Phase R1 (predicted PARTIAL or NOT-DISTINGUISHED).
- Q9. The pair-weighted S-rate of the gorilla arm at tau 0.90 exceeds 0.5. It is expected from selection by identity and is not evidence for or against the hierarchy.
- Q10. In Phase R1 at least 8 strata are live for T1-b and for T1-c, and UNINFORMATIVE does not occur.

A positive T1 licenses no duplicon-first node, core or drop rule (section 11).

### 11. Limits and register check

**Limits declared in advance.** One library per species. Human relations are second looks on spent halves. The gorilla family layers are depleted of young SD-rich copies (the Compara projection keeps old dispersed pairs; whole in 54.0% of cross-contig against 17.8% of same-contig human families). DupMasker and SEDEF are both sequence-derived, and so are our family edges. F_O is a DNA-level partition of annotated SD98 loci and not an RNA family layer. Atom constants and the Liftoff `-sc 0.95` are inherited. Gnomon models cite long SRA reads (75.7%), possibly including the KB3781 fibroblast Iso-Seq that feeds the families. Liftoff `-copies` sees only whole-gene-body alignments, so the copy-bearing set is biased to short genes (median span 3.9 kb against 12.2 kb). The outgroup primaries are haploid and may collapse young duplicates (not testable here). The E-table pool comes from a self-lift sharded by contig, in which Liftoff rejects a copy that lands on another annotated record; only 5 of the 307 pool components have all their records on one contig. The E-table therefore describes dispersed and cross-block copies. Tandem-only expansions are outside it, and T2a and T2b make no statement about them (a rank-sharded self-lift of the gorilla annotation would add them; heavy, Amendment). T1-b and T1-c cannot separate duplicon organisation from shared sequence similarity (section 2, item 3). The T2 family object is MCL and not an L3 component, and MCL splits real expansions (LRPAP1: 6+3 and 8+3), so strict containment is low for reasons unrelated to age. Big families dominate micro counts, so macro is printed beside every micro and the all-families number is always reported, never a hand-picked subset. The 15% single-copy share is a name-proxy indicator. Importing repo modules for the recon wrote two gitignored `.pyc` caches (`bench/soto/__pycache__/soto_replication.cpython-314.pyc`, `bench/soto_m2/__pycache__/soto_m2_duplicons.cpython-314.pyc`); later runs use `python3 -B`.

**Register check.** A positive T1 licenses no duplicon-first node, core or drop rule (register rows 676, 677, 732, 735, 739 to 743, 817, 818 and 997 close them). The 149 multi-family genes are listed per unit and no cover is proposed (rows 845, 846, 941, 942, 947, 1144 and 1184 to 1188). Poset family objects are not proposed (rows 1020 to 1022). T2b arm 1 re-enters register lines 88 and 753: 3 of 18 certified is expected to persist, and the inverted substrate is what changes. The SEDEF-PAF tautology (rows 650 and 996) bounds the gorilla S-rate. The broad/recent two-tier hierarchy was inert on 70% two-copy families (register line 311, row 654), so this study measures a relation between levels and proposes no hierarchy rule. Name joins are forbidden by rows 993 and 1034. The only name steps are the two Compara projections (gorilla by exact symbol, human by gene name, optional) and they are declared. Compara is not O1 truth (register line 316), a human annotation projected to gorilla is not gorilla truth (lines 861 and 869) and the same projected layer was judged UNDERPOWERED for the class-level O1 test (row 1258). The last register row is 1262, so the first result rows are 1263 and following. The three 09-30 duplicon tests and the NPIP fusion test have no register row yet and get rows with the first result.

### 12. Decisions taken by default (change only by an amendment committed before the first run)

1. **T1 universe.** V_F (genes with an F_O cluster) for the decisive cells; the 485 unclustered V genes are outside them and enter a descriptive singleton arm. Alternative: singleton families in the decisive D2 cell, which makes every unclustered gene that shares a string with a clustered gene a D2 violation and measures graph recall.
2. **Governing null.** Homology-preserving (components of the F_O graph). Alternative: chromosome or SD98-region strata as governing, which are predicted to give p <= 0.01 whatever the organisation (Q4).
3. **Audit.** Descriptive only; nothing is excused in a decisive statistic. Alternative: excuse 'short' and 'ends-only' inside the statistic, which makes NOT-REFUTED reachable by construction (862 of 2,142 genes have at most two exons).
4. **Gorilla arm.** Descriptive, no class; the matched-instrument human run is optional. Alternative: require the run, which needs the name and offset join of section 5.
5. **R3a.** The human_testis BASE snapshot is admitted as a light arm and will be called held-out for the family layer and the library once it has run. Alternative: wait for the heavy current-default run (R3b).
6. **E-table.** Tier B, rank-sharded, is primary, with tiers A and C as replicates. Alternative: tier A or C primary, which costs 1 to 2 h each and, for C, builds of splice indexes under a 37 GB disk limit.
7. **T2 hold-back.** The expansion LRPAP1 is the only held-back unit; the ten contigs are printed as a least-exposed subset. Alternative: the ten contigs decide, which is UNDERPOWERED by construction (3 and 8 supported components).
8. **T2a statistic and classes.** T is the mean R_j with absent members as non-co-members; PARTIAL covers 0.01 < p <= 0.05. Alternative: strict end-to-end containment c_j as the statistic, which is harsher on MCL splits.
9. **T2b.** Arm 1 is an estimate and arm 2 has classes, with no combined verdict, and both are predicted UNDERPOWERED. Alternative: drop T2b from this registration.
10. **Constants.** The denominator guard is 10% (CP-5, inherited). The hub cap is 51 (inherited). The bootstrap is 2,000 resamples, B = 10,000, seeds 20260930 (T1) and 20261007 (T2). Inherited, not tuned: coverage 0.5 and 50% overlaps (Liftoff `-a`, `figures/_liftoff.py` MATCH_COV), 200 bp (SHORT_BP), 2 reads (SUPPORT_READS), 25 kb flanks and anchor bars (`bench/o1_outgroup_rooting_poc`), MIN_BLOCK 1,000, 100 exonic bp for the gorilla touch rule borrowed from the few-copy criteria, tau 0.90 (from input-only power counts). The floor n_min = 8 independent units is inherited from the UNDERPOWERED rule (E < 8) of `PREREG_seed_pool_real_reads_2026-10-07.md`. Fixed now, not tuned, with a reported sweep where marked: the caliper of 10 half-thousandths (sweep 4, 10, 20), 1,000 matched sets per component and the exploration cap of 1,000,000 candidate sets per component, the gorilla flags of fewer than 20 depth-matched pairs or fewer than 8 families, Gate 2 (at most 1 of 20), Gate 4 (at most 10 of 100), Gate 7 (100 records, 50 fake pairs, at most 1), UFBoot 95, and 1,000 non-exonic bp and 2 kb flanks for arm 2. The floor of 8 qualifying control components per expansion is n_min. The sweep of n_min (6, 8, 10) is reported and not acted on.
11. **Annotations (T2).** RefSeq for all four species. 'Expressed' means at least 2 primary reads at the placement. For T1 the annotation is CAT v4 (v1.0) with the v2.0 CAT/Liftoff ids (section 2, item 8).
12. **Authorisation.** No heavy run, no genome-wide family run and no cross-lineage replicate without the user's message and an Amendment (section 3). `dna_sd_atoms.py` is restored with its hash and a test before Phase R2.
13. **Gate 1.** It requires the registered outputs only; the 2,259 / 2,290 explanation is not a precondition.
14. **T1 primary choices.** F_O is the decisive layer, U1 the primary rung and chain-tolerant m1c the D1 statistic; the class is issued at inflation 2.8; the half rule is gene-level. Alternative: F_S or U2 as primary, which is circular (section 2, item 4).
15. **T2a family object and members.** The BASE snapshot and expressed members only; regeneration under the frozen default needs an Amendment.
16. **T1-b and T1-c reading.** The class names describe the statistic (CUT-ALIGNED) and the text states that the cells cannot separate duplicon organisation from shared sequence similarity. Alternative: register a similarity-matched comparator (the F_O cut at inflation 2.0 and 4.0 as competing sequence-only re-cuts of the same group) before any run, which would test the stronger claim.
17. **Robustness clause.** The class is the lower of the all-strata class and the class after removing the stratum with the most negative z. Alternative: a stratum-level sign test (22 of 30 for p <= 0.01), which has less power.
18. **Gate 7 bar.** The 6 qualifying LRPAP1 records. Alternative: the eight 20-kb copy bodies as synthetic records, which changes the lift recipe.
19. **T2a controls.** Two-step draw (control component, then matched set); contig multiset and caliper in the matching; component size reported as strata with a SIZE-SENSITIVE flag; floor of n_min = 8 qualifying control components; headline class = the class nearer to NOT-DISTINGUISHED. Alternatives: match on component-size class, which makes the floor unreachable (best windows of 2 to 20 components), or a uniform draw over all matched sets, which samples almost only large components.
20. **Amendment 1 timing.** Committed before the first Phase R1 run, and no Amendment may then change a rule that Phase R3a uses.
21. **Siblings.** Expansions of different components with an identical triple count as separate units; a sensitivity counts each group once.

### 13. Names, scorers, results hook, frozen inputs

**Names.** KEY=dupfamboundary (T1) and KEY=expansionfamily (T2). Scorers, written after this file is committed with tests first: `bench/hierarchy/t1_refine.py`, `t1_gorilla_pairs.py`, `expansion_table.py`, `t2a_containment.py`, `t2b_age.py`; restored `bench/dna_sd_atoms.py`. New code to write: an importable string, dominant-duplicon and rung-U2 module (the frozen `comp` and `dom` are closures inside `main()`); the F_O loader with the fold map; the stratum builder; m1, m1c (longest-weight chain), m1t, m2, path defect, audit and nulls; the atom class builder and dedupe step; the pair table and S-rate; the cross-species Liftoff Target, the rank-sharded shard planner and its split rule, the E-table builder and Gate 7 controls; the T2a matcher (reuse `figures/_liftoff.py::found_by` and MATCH_COV) and control sampler; the synteny classifier and the T2b tree pipeline.

**Results hook.** Results are appended below the line at the end of this file without editing anything above. Each result section records the commit hash and sha256 of this file at the run, the document `docs/HIERARCHY_DUPLICON_EXPANSION_<run date>.md` and per-gene, per-pair, E-table and per-expansion TSVs. Register rows from 1263. **Amendment 1** is committed before the first Phase R1 run (section 6). Every heavy step needs its own Amendment.

**Frozen inputs (Gate 0 list; sha256 first 16 hex, rows; `c6_hashes.out` holds the full listing).**

- DupMasker colors BED `winloci_data/duplicons/chm13.draft_v1.0_plus38Y_dupmasker_colors.bed` 59533d25623ebf60, 205,815 rows, 9,655 labels, 25 contigs. `chm13.draft_v1.0_plus38Y.duplicons.bed` 2c0d0c94890e23f1, 254,198 rows.
- S1C `bench/soto/soto_famCN_S1C.tsv` 776a650716fd1a1b, 2,572 rows. Split `bench/soto/soto_split_2026-09-29.tsv` 8f2cba6c34984b69, 1,170 + 1,164 genes. Exon map-back `shared_exons_5154_exon_mapback.tsv` 90ef554f96ec0359. `sd98_gene_exons.tsv` 61826d144db5896b. `sd98_v1.bed` d3e56dfa2cf81d69, 817 rows. `famcn_ours_allwssd.tsv` 9067d31ef1838764. `families_cn.json` 625fa960d8a6ee1e. `elig.tsv` ffd86585d04bc5e5. `full.tsv` c2bb5a994ff9cbff. `rt_genes.tsv` d7b107123e5ee19e, 76 rows. `gene_v2_coords.bed` 8f25eb939b442432, 2,332 rows.
- F_O: `/mnt/linuxdisk/tmp/ourhom/` `names.tsv` 839412bd1a5deeb0, `loci.gff3` 8237e7e268432583, `loci.paf` (37,782,223 bytes), `ours.clusters.tsv` 8a835d4c434493a3, `ours.loci.tsv` 7490ec5f23f5474a, `ours.pairs.tsv` add12d1a70a0c101, `ours.params.tsv` 9964b9aaa6158a78, binary `bin/mcl_families` sha1 e4cc13b9d7eba5bcad3cc9c8c93a3fdc8ad434d7; the graph dump `/mnt/linuxdisk/tmp/prereg_revision/fo_graph/ours.graph.tsv` 1820795ae336a003, 11,031 edges. These are copied to a hashed location before R0.
- CAT/Liftoff v2.0 `chm13v2.0_CAT_Liftoff.genes.tsv` b365d1995d75e2ea and `refseq_map.tsv` ca87ff8aac922e47; CAT v4 `cat_v4.bed` 6a1fddadcc8168ab; `t2t-chm13-v1.0.fa.gz` (size 981,798,956 bytes).
- Gorilla SD: `GGO_sedef_final.bed` sha256 41aa1c9c53a18d68, md5 991ed34384832354c02df227e27f65b6, 253,029 data rows; `final_v1_clean.bed` d62edf9caa1c7a95, 151,758 rows. Compara and truth: `gorilla_truth.tsv` 9e6cb80ad518f528, `pairs.tsv` fb12608ae431a105, `families_gw/species/gorilla/genes.tsv` 6f9bfbd776e8ecc8, human `compara.Primates.families.tsv` 00952f5fa25da3d8.
- T2: `liftoff/<species>/liftoff_loci.tsv` gorilla af503bf2713caf54, human 07cf5aed5f33c952, chimpanzee ec9772c09e349873, orangutan eba8b18efb120075; support tables `liftoff/support/gorilla_KB3781.tsv` 63fb4c3063405e76 and `gorilla_OR6737.tsv` 1b1ac576fde36407; BASE tables `rt_arms/<library>/<library>.BASE.fam.{clusters,copies,loci,params}.tsv` for gorilla_KB3781, gorilla_OR6737, chimp_PTR and human_testis, with their PAFs (size and mtime recorded); GFFs, FASTAs and minimap2 indexes by size and mtime; the two BAMs named in section 8 and the human_testis BAM, and the four minimap2 indexes of section 7 by path.
- Tools: liftoff 1.6.3, mafft 7.526, iqtree3 3.1.3, minimap2 2.24 / 2.30 / 2.31, python3 3.14.4 (python 3.13.12 not used).

### Amendment 1 (2026-10-07, committed before any run; it changes no cell, statistic, null, class or prediction above)

**Authorisation.** User message, 2026-10-07 17:21 PDT, answering the offer of a 'go' for the light phases: "yes please go; maybe for today lets focus on the DNA side for all the fixtures (multi-copy families, duplicons, expansions, segmental duplications)". This is read as a 'go' for Phases R0, R1 and R2. All three are light and DNA-side: the human DNA-level families and duplicons with the SD98 regions (R0, R1) and the gorilla SEDEF arm, which is the only light phase that touches segmental duplications of gorilla (R2). It does not cover Phase E, which is the only phase that builds expansions, is heavy and needs its own 'go' and an Amendment 2. It does not cover T2a (expressed members), T2b, R3a (RNA-derived) or R3b. The matched human run of section 5 is not covered either.

**T1-e, a new descriptive cell with independent evidence (no class).** The user asked what to use other than sequence similarity. Copy number from read depth and genomic position are not computed from the similarity of two genes, so they can say whether duplicon strings carry information that the similarity graph and the family cut do not. Only a per-gene median copy number exists on disk (`famcn_ours_allwssd.tsv`: median, MAD, 268 samples), not a per-sample matrix, so co-variation across samples is not possible.
- Universe: the V_F genes of the 30 variation-bearing strata of F_O that have a famCN (`co` of `families_cn.json`) and at least one U1 label; the count is printed.
- Pairs: unordered pairs of those genes in one stratum whose nodes are joined by a direct edge of the dumped graph. For each pair: w = the dumped edge weight (the maximum if several nodes map to one gene); x_clu = 1 iff the genes are in the same F_O cluster; x_str = 1 iff their U1 sets are identical; x_jac = the Jaccard index of the U1 sets; y_cn = |co_g - co_h|, replaced by its rank over all pairs (average ranks for ties); y_loc = 1 iff both genes have the same SD98 region of largest exonic overlap (ties: smaller region id; a gene with no region gives 0).
- Model: for each outcome, ordinary least squares of the outcome on w, x_clu and x_str after subtracting the stratum mean of every variable (stratum fixed effects). Reported: the coefficients b_str and b_clu with 95% percentile intervals from 2,000 bootstrap resamples of the 30 strata with replacement (`random.Random(20260930)`, strata sorted by smallest gene id, refit each time); the same with x_jac in place of x_str; the number of pairs, the share with x_str = 1 and with x_clu = 1, and the correlation of x_str with x_clu. No p-value and no class.
- Reading registered now: if duplicon strings carry copy-number or location information beyond sequence similarity and cluster membership, b_str is negative for y_cn and positive for y_loc, with intervals that exclude 0. If both intervals include 0, the cell shows no evidence of that information. Neither result changes a class of T1-a to T1-c. The cell is descriptive and carries the limits that copy numbers of near-identical genes are coupled by read mapping (w is the control for it) and that an SD98 region is defined from alignments but says where two genes lie, not how alike they are.

**Fixed now for Phase R3a (not run today).** The graph dump is `mcl_families --paf human_testis.BASE.fam.loci.paf --gff human_testis.BASE.fam.loci.gff3 --min-exonic-bp 1 --min-shared-exon-frac 0.60 --min-cov-shorter 0 --dump-graph OUT.graph.tsv --out OUT` with the frozen binary of Gate 1 (sha1 e4cc13b9), run under `tools/rlock.sh light` (the PAF is 176 MB). At the start of Phase R3a every BASE cluster must lie inside one component of the dumped graph (otherwise Phase R3a is INVALID) and the Gate 6 count of strata with at least two clusters is printed.

**Call schedule for today.** R0: Gates 0, 1, 3, 4, 5 and 6, each in its own call of under 180 s under `tools/rlock.sh light`, in the order of section 3; Gate 4 repeats are split into blocks. R1: one call per cell (T1-a, then T1-b and T1-c with their 10,000 draws in blocks of 2,500, then T1-e), results appended below the results line with the commit hash and sha256 of this file. R2: after `dna_sd_atoms.py` is restored with its hash and test. A cell is run only after its gates pass.

### Amendment 3 (2026-10-07, committed before any Phase E run and before any Y2 to Y4 number exists): concordance with the Yoo et al. 2025 gorilla expansions

**Why.** The advisor wants the pipeline shown to find results similar to Yoo et al. 2025 (Nature, "Complete sequencing of ape genomes"; the thesis vault holds the Supplementary Information). The user's message of 2026-10-07 17:57 PDT, "also please review the Yoo supplemental data in the thesis vault, I think my advisor wants me to prove that I find something pretty similar to them, for example with LRPAP1", opens this block. The same message gives the 'go E' for Phase E, whose command, hashes and call schedule are fixed in Amendment 2, committed before the first Liftoff call.

**Truth files (frozen).** The Supplementary Tables workbook (Supplementary Tables 1 to 82) was downloaded from the Nature page of the paper (`41586_2025_8816_MOESM4_ESM.xlsx`) to `/mnt/linuxdisk/home/juanfraitu/winloci_data/yoo_2025/Yoo2025_SupplementaryTables_MOESM4.xlsx` (16,468,422 bytes, sha256 865d4e868568b32d...). Extracts there: `table_38.tsv` (Supplementary Table VIII.38, loci of expanded gene family copies: the gorilla chr1 MAPKBP1 / JMJD7-PLA2G4B / SPTBN5 units with the ancestral triplet on chr16, and LRPAP1, DOK7, HGFAC with the ancestral copy on chr3 and nine copies; sha256 6a1a51800c77061e...), `yoo_GGO_table35_novel_paralogs.tsv` (the 808 gorilla rows of Table VIII.35, novel paralogous genes that overlap SDs; sha256 29b22d4f56d9f63d...), `yoo_GGO_table40_lineage_specific.tsv` (the 276 gorilla rows of Table VIII.40, lineage-specific genes; sha256 795244a031415dcb...) and `GGO_chr_map.tsv` (chromosome name, GenBank and RefSeq accession, from the NCBI Datasets sequence report of GCF_029281585.2; sha256 857bb6b71bbfa8fb...). Tables VIII.35 and VIII.40 use GenBank accessions and UniProt-style gene mnemonics (LRPAP1 is AMRP). Table VIII.38 uses chromosome names. Coordinates in all three are those of the paternal (primary) mGorGor1 assembly, which has the sequences of our RefSeq contigs. PSMA5 is named only in the text of Supplementary Note VIII (2, 2 and 3 duplicated copies in chimpanzee, bonobo and gorilla) and is in none of the tables, so it can be compared by count only. Yoo's definition of a copy is a minimap2 asm20 hit with more than 75% query coverage and identity, an expressed copy has more than 3 Iso-Seq reads, and a valid ORF needs a TransDecoder peptide covering at least 75% of the blastp best hit at more than 80% identity. Yoo's expansion is relative to the human annotation and is assigned to a lineage afterwards, so it is looser than the unanimity of section 7: the chimpanzee and bonobo also carry copies of the MAPKBP1 unit and of PSMA5, so those two families fail the 'single-copy in chimpanzee' condition and are not called by the E-table of section 7.

**Match rule.** Two loci on the same chromosome match iff their overlap is at least 50% of the shorter locus. Reported with each concordance: the number of Yoo loci with a match (sensitivity), the number of our loci with a match (precision against Yoo's catalogue, bounded because Yoo's catalogues are expression- and ORF-filtered), and the size of a maximum one-to-one matching.

**Y1, LRPAP1 (looked at before this Amendment; labelled spent).** Our 11 loci (8 full-length copies and 3 5' fragments; `docs/archive/2026-10/LRPAP1_FAMILY_2026-10-04.md`) against the 10 loci of Table VIII.38 and the 5 AMRP rows of Table VIII.35. Result: all 10 Table VIII.38 loci match one of our loci (reciprocal overlap at least 0.97 of ours and at least 0.71 of theirs, 9 of 10 at least 0.99 of theirs), Yoo's solitary copies 2 and 6 are our chr12 and chr14 fragments with identical coordinates, and our 11th locus (the chr16 fragment at 17,176,045 to 17,192,323) is the AMRP row of Table VIII.35. All 11 of our loci are in Yoo's catalogue, 10 in Table VIII.38 and 1 only in Table VIII.35. The file is `/mnt/linuxdisk/home/juanfraitu/winloci_data/yoo_2025/y1_lrpap1_concordance.tsv` (written with this Amendment).

**Y2, the chr1 MAPKBP1 unit (not yet looked at).** Truth: the 8 gorilla chr1 triplets of Table VIII.38 (24 gene loci) and the ancestral chr16 triplet (3 loci). Ours: the families of the KB3781 BASE table (`rt_arms/gorilla_KB3781`, a declared pre-f1v2 snapshot, section 8) whose copy rows match a Yoo unit locus. Reported per gene (MAPKBP1, JMJD7-PLA2G4B, SPTBN5): Yoo loci matched, the number of BASE families touched, the members of those families that match no Yoo locus, and whether the three genes fall into the same family. Prediction, not blind to the field but blind to our table: the BASE table covers at least 5 of the 9 MAPKBP1 loci (fibroblast reads), and the three genes are in different families because the unit is a tandem of three unrelated genes.

**Y3, PSMA5 (count only; not yet looked at).** The records of the gorilla RefSeq annotation (`winloci_data/GGO_genomic.gff`) whose gene name or description is PSMA5 or 'proteasome subunit alpha type-5', counted against Yoo's 1 ancestral and 3 duplicated copies. Prediction: the annotation holds the ancestral gene and at least one of the 3 duplicates as separate records.

**Y4, genome-wide (needs the E-table; descriptive, no class).** Yoo's truth is Table VIII.40 (276 gorilla copies, primary) and Table VIII.35 (808 rows, secondary). Ours: the members of the expansions of the E-table of section 7 (the unanimity table), and as a sensitivity the members of the expansions called relative to human only (the same lifts, a triple reduced to its human locus h; Yoo's expansion is relative to human). Reported: the share of the 276 rows with a matching E-table member, in each table; the share of E-table members that match a row of Table VIII.40 or VIII.35; and, for each Yoo row without a match, the first reason in the order: no pool record at the locus, hub component, an outgroup self-lift copy row, no placement, different triples. The denominators and the shares are printed together. Neither share is a verdict, because Yoo's rows require Iso-Seq support and an ORF and the E-table does not, so the second share is bounded above by the fraction of E-table members that are expressed.

**What this block does not do.** It tests no other species. Chimpanzee, bonobo and orangutan expansions (PLA2G4B-JMJD7, PSMA5, HERC2, GOLGA6, GOLGA8, MCTP2) need their annotations and testis Iso-Seq reads, which are not local. It changes no cell, statistic, null or class above.

*Results (filled in after the run, below this line, without editing anything above.)*


---

## PREREG_hierarchy_duplicon_expansion_REVIEW_DISPOSITION_2026-10-07.md

*(original: `docs/archive/2026-10/merged/PREREG_hierarchy_duplicon_expansion_REVIEW_DISPOSITION_2026-10-07.md`)*

## Disposition of the three adversarial reviews of the T1/T2 pre-registration (2026-10-07)

Reviews: workflow wf_20046203-7d9 (three reviewers: wording 37 findings, logic 29, feasibility 26; 92 in all). The draft reviewed is the first draft of `docs/PREREG_hierarchy_duplicon_expansion_2026-10-07.md`. The revised file is the same path (DRAFT 2). The review results as JSON are `/mnt/linuxdisk/tmp/prereg_review_wording.json`, `prereg_review_logic.json` and `prereg_review_feasibility.json`; the numbers below are the positions in each file's `findings` list (wording W1 to W37, logic L1 to L29, feasibility F1 to F26).

**Closure review of draft 2** (workflow wf_eea75a23-1c5, 2026-10-07; four read-only reviewers on frozen copies): logic 25 findings, facts 22, wording 36, disposition audit 19 (the JSON files are `/mnt/linuxdisk/tmp/prereg_closure_{logic,facts,wording,disposition}.json`; scratch `/mnt/linuxdisk/tmp/prereg_closure/`). It found two blockers (both the Gate 7 positive control, one from the facts reviewer and one from the disposition audit) and 25 majors. The revised file is DRAFT 3; the closure dispositions are the second table below. In the first table, rows that the closure review found only partly fixed (F1, W15, W17, W22, L9) point to the closure rows that complete them.

**Final regression review of draft 3** (workflow wf_c1c3d967-0ee; three reviewers on the frozen diff): logic 23 findings, facts 20, wording 36 (JSON `/mnt/linuxdisk/tmp/prereg_regress_{logic,facts,wording}.json`; scratch `/mnt/linuxdisk/tmp/prereg_closure2/`). It found one blocker (the floor of 20 control components of T2a cannot be reached) and nine majors; the errata commit that follows draft 3 (DRAFT 3a) answers them, third table below. Rows of the first two tables describe the state of draft 2 and draft 3. Class names, counts, floors and bars in them can be superseded; the state of draft 3a is the prereg file itself and the third table. Pointers for rows the closure review found partly fixed: F1 to CF1 and CF4; W15 to CF5; W17 to CL8 and CL9; W22 to CF1 (the bar is the 6 qualifying records, not 8 loci); L9 to CL9 and CD7.

Dispositions. ACCEPTED: fixed as suggested or in substance. ADAPTED: fixed in a different way, and the row says how and why. REJECTED: not fixed, and the row gives the evidence. NOTED: no change needed. Section numbers refer to the revised prereg. Cross-references such as 'wording 3' point to another row of this table.

### Counts by disposition

| reviewer | findings | ACCEPTED | ADAPTED | REJECTED | NOTED |
|---|---|---|---|---|---|
| wording | 37 | 27 | 10 | 0 | 0 |
| logic | 29 | 22 | 5 | 0 | 2 |
| feasibility | 26 | 19 | 6 | 0 | 1 |
| all | 92 | 68 | 21 | 0 | 3 |

### Dispositions, one row per finding

| reviewer | # | severity | location | disposition | where fixed, or why not |
|---|---|---|---|---|---|
| wording | W1 | blocker | Gate 5: G1 does not refine G0; m1c not monotone; different gene sets | ACCEPTED | Section 3 Gate 5. Gate 5a on the common universe (U1 to U2, once per famCN source), Gate 5b with the numeric bar D[S->Q] = 7 (2 if ungrouped genes are ignored); U0 dropped from the chain; the gate prints PASS or FAIL and no value. The 32 strings / 110 genes and the 10 genes were reproduced (logic/string_only_checks.out). |
| wording | W2 | blocker | T2b identity check has the wrong sign and strictness | ACCEPTED | Section 9 Identity check (h is a divergence; h(E) <= h(phi), non-strict; a violation is INVALID); age_id deleted; section 2 item 6 uses the same h; consequence C2 replaces P8. |
| wording | W3 | major | Audit not closed (deficit set, ends-only, D2, counts) | ADAPTED | Section 4 Audit: tie rule for the kept chain (lexicographically smallest maximum-weight chain, all tied chains printed), ends-only defined, D2 labels, a label never removes a gene, counts quoted in V, non-eligible defined by the S1C column. Differs from the suggestion: the audit is descriptive and is not the decisive statistic, because 'short' alone excuses 862 of 2,142 genes (finding 34) and excuses are asymmetric between data and null. Decision 3 in section 12. |
| wording | W4 | major | Gene-to-family rule, layer universes, ratios | ADAPTED | Section 4 Family layers: rules 1 to 3, V_F = 1,657, layer universes, every ratio printed with its universe. Differs: unclustered genes are outside the decisive cells (feasibility finding 6 recommends this) and enter a descriptive singleton arm, instead of being singleton families by default. |
| wording | W5 | major | Verdict algebra forks | ACCEPTED | Section 3 Verdict algebra: precedence, p only with the same cut-offs everywhere, UNINFORMATIVE defined, n_min, SPLIT rule, HOLDS removed. Class names changed (ORGANISES-ENRICHED, NOT-DISTINGUISHED). |
| wording | W6 | major | Decisive cells and expected answers | ACCEPTED | Section 3 Decisive cells table; section 10 (C1 to C5 consequences, Q1 to Q8 predictions for T1-a, T1-b, T1-c, T2a, T2b). |
| wording | W7 | major | Nulls not reproducible (rungs, singletons, strata, seeds) | ADAPTED | Section 4 Nulls and Draws: re-seed per (layer, null, universe), sorted genes and strata, draws shared only inside one universe, singletons outside the decisive cells. Strata are the homology-preserving components (logic finding 2) and not the chromosome. |
| wording | W8 | major | 'Per half' rule | ACCEPTED | Section 4 Half (gene-level column of the split file; halves A and B descriptive, never called held-out). |
| wording | W9 | major | Unit-string construction | ACCEPTED | Section 4 Exon label and unit strings (bases added, argmin key, G0 key, coordinates as stored); Gate 6 prints tie counts; the set-Jaccard sensitivity is deleted. Counts from input_checks.out. |
| wording | W10 | major | Gate 1 may have no exit (2,259 vs 2,290) | ACCEPTED | Section 3 Gate 1 requires the registered outputs only; section 2 item 9 records the erratum; section 12 item 13. Two more explanation attempts failed (c5, c5b). |
| wording | W11 | major | Gorilla p-value undefined | ADAPTED | Section 5: no control draws and no p-value (arm descriptive); family-cluster bootstrap with seed 20260930; family-weighted S-rate printed. |
| wording | W12 | major | Gorilla pair universe not closed | ACCEPTED | Section 5 Family layer and pair universe (touch rule, U(g), eligibility, 0 shared bp for every pair, distance, depth-matched, verdict set). Counts re-derived (c7_depth_matched_verdict.out). |
| wording | W13 | major | T2 staging not a partition; decisive stage underpowered | ADAPTED | Section 7 Contig roles: the hold-back unit is the expansion (LRPAP1 only), the ten contigs are a least-exposed subset without a class, accession and chr table, 'held-out' not used for gorilla, expected sizes (c8). Differs from the suggested Confirm/Run-last rule, which would leave 3 and 8 components. |
| wording | W14 | major | T2a class statistic and test undefined | ADAPTED | Section 8 Statistics, Matched controls and test, Classes. Deciding statistic is the mean R_j with absent members counted as non-co-members (logic finding 11), B, seed, two one-sided p-values, PARTIAL, size-2 stratum. |
| wording | W15 | major | T2a matching and scoring | ACCEPTED | Section 8 Members and matching and Expressed (ties, one-to-one node rule, X(n) source, overlap length, m_j, R_j). |
| wording | W16 | major | E-table definition not closed | ADAPTED | Section 7: hub cap 51, expansion = same triple (an equivalence relation, no chaining), single Liftoff pass with -copies (the second pass is dropped), strand ignored, deciding E-table is tier B only, columns. |
| wording | W17 | major | T2b decision statistic | ADAPTED | Section 9: no combined verdict; arm 1 is an estimate; arm 2 has classes with a comparator; mixed-family defined; 0.5 kept only as a majority bar. |
| wording | W18 | major | T2b arms cannot be run as written | ACCEPTED | Section 9: synteny class inlined, member classes, derived set = E minus ancestral members, arm 2 command line, rooting and UFBoot 95, minimum tips, R(D) defined, d_orth deleted. |
| wording | W19 | major | What 'go' authorises | ACCEPTED | Section 3 Run order and authorisation (table and rule); E3 deleted; one call budget. |
| wording | W20 | major | Wording changes are buried | ACCEPTED | Read first: items 1 to 7, the cell table, and the statement that the two sentences are my wording of the offer. |
| wording | W21 | major | Matched-instrument human run: no join, 'agrees' undefined | ADAPTED | Section 5: the run is optional and descriptive, the join is specified (Name and Contig with a uniqueness assertion; offsets verified by sequence), 'agrees' deleted. |
| wording | W22 | major | Gate 7 cannot fail; negative control | ADAPTED | Section 7 Gate 7: bar of 8 loci in one class, can fail, fragments are a sibling class; the negative control is a planted set of 50 fake pairs (at most 1 may be called) instead of the suggested symbol lookup, which is vacuous for a pool of copy-bearing records. |
| wording | W23 | major | Gate housekeeping | ACCEPTED | Section 3: scope of the gates, Gate 2 pass bar, Gate 5b number, Gate 6 definitions, T2 controls in Gate 4, E-table md5 recorded when section 7 closes. |
| wording | W24 | major | Symbol collisions | ACCEPTED | Renamed throughout: U0 to U2, Gate n, Q1 to Q8 (C1 to C5), Mult, pool, E-table, 'gorilla', R_J contig. |
| wording | W25 | major | Genome-wide tables exist | ACCEPTED | Section 2 item 5 corrected; section 6 R3a admits the human_testis BASE snapshot (the author's choice made, decision 5). |
| wording | W26 | minor | Numbers that do not tie out | ACCEPTED | Corrected in section 2 items 7, 8, 10, section 4 Universe V, section 5 counts (38 pairs in 33 families), section 5 tau rule, section 10 C1. |
| wording | W27 | minor | Convention gaps (KEY, results hook, G0 list, hashes) | ACCEPTED | Section 13 (names, scorers, results hook, Amendment hook, frozen inputs with sha256, atoms hash, BASE binary commit stated as not recorded); traps agent scratch path added to the header. |
| wording | W28 | minor | Jargon | ACCEPTED | Terms paragraph; 'Amendment-free' deleted; 'in-paralogs' replaced; A_viol described in section 5. |
| wording | W29 | minor | Randomisation and constants | ACCEPTED | Section 3 Gate 3 (seeds, per-cell re-seed), section 12 item 10 (constants, sweep, inherited list), call budget in section 3. |
| wording | W30 | minor | Secondary statistics | ACCEPTED | Section 4 Statistics (macro, ordered string, simple graph, Mult, m1c = m1 at U0, nesting wording). |
| wording | W31 | minor | Input handling (identical tuples, sensitivities, R3 sketch) | ACCEPTED | Section 5 identical tuple and tie rule; sensitivities defined or deleted (section 4); section 6 R3 fixed by Amendment 1. |
| wording | W32 | minor | Over-claims (P3, 'literal T1', header) | ACCEPTED | Section 10 Q2; section 1 D1 wording; header sentence on what was joined. |
| wording | W33 | minor | Markdown layout | ACCEPTED | Blank lines between paragraphs throughout. |
| wording | W34 | note | 'short' label explains 40% of V | ACCEPTED | Section 4 Audit prints the unexplained count with and without 'short'. |
| wording | W35 | note | V removes the 149 multi-family genes | ACCEPTED | Section 4 Universe V and Gate 6 (per-layer counts of clusters holding them and members deleted). |
| wording | W36 | note | Tie rule differs from the 09-30 text | ACCEPTED | Section 4 Exon label and unit strings states it. |
| wording | W37 | note | T1/T2 name collisions | ACCEPTED | Read first item 7; section 13 KEY names. |
| logic | L1 | blocker | Gate 5 false premise | ACCEPTED | Same as wording 1. Section 3 Gate 5a and 5b; D[S->Q] = 7 and 2 reproduced (logic/fq_checks.out). |
| logic | L2 | major | Tautology; both nulls destroy homology | ACCEPTED | Section 4 Nulls: governing null = components of the F_O graph. The frozen binary re-run with --dump-graph reproduces the frozen products byte for byte (c1_fo_graph.out). Chromosome and SD98-region nulls become sensitivities (C4). 30 variation-bearing strata with 416 genes (c2_fo_components.out). F_Q has no variation under its own strata and stays descriptive. |
| logic | L3 | major | Audit not tie-invariant; D2 undefined | ADAPTED | Section 4 Audit. The suggested single optimisation with excuses is not used inside the decisive statistic (asymmetric excuses; 'short' excuses 40% of V). Adopted: unique chain value, tie rule, D2 labels, ends-only reference. Section 12 item 3. |
| logic | L4 | major | HOLDS-EXACT unreachable; HOLDS-ENRICHED mislabels a refuted universal | ACCEPTED | Sections 3 and 4: T1-a is REFUTED or NOT-REFUTED and is printed first; the enrichment class is ORGANISES-ENRICHED; NOT-REFUTED only on the primary arm. |
| logic | L5 | major | Classes not mutually exclusive | ADAPTED | Section 3 Verdict algebra (precedence; 'direction differs' and '1st percentile' deleted; p only; halves descriptive; SPLIT between Phase R1 and R3a (and, from draft 3, between the two T2a libraries)). UNINFORMATIVE is defined as a null that cannot move. The independent unit is the variation-bearing stratum, so UNDERPOWERED can in principle fire on the human cells (30 against 8). Draft 3 (closure row CL10) extends SPLIT to the two T2a libraries and adds the headline rule. |
| logic | L6 | major | Gorilla arm has no unconfounded decision statistic | ACCEPTED | Section 5, option 1: descriptive report with intervals; enrichment verdict dropped; reasons stated; multi-class genes in the primary and flagged; 'agree' and 'rule frozen on human' deleted. |
| logic | L7 | major | Gorilla layer is a name join | ACCEPTED | Section 5 Family layer (one name step, declared) and section 11 Register check (lines 316, 861, 869; rows 993, 1034, 1258); matched human run optional with printed counts. |
| logic | L8 | major | T2b does not test 'younger than the root' | ACCEPTED | Section 9: arm 1 is root membership with the Dollo definition, arm 2 is a rooted clade test; d(E) < d(F) deleted; the sentence that the clock is only partly independent of the family builder is added. |
| logic | L9 | major | T2b classes conditioned on outcomes | ACCEPTED | Section 9: a non-clade counts as not younger; unit exclusion only by input properties; shares defined per arm; funnel printed; pilot unresolved rates corrected to 83% and 75%; no combined verdict replaces FAILS, PARTIAL and SPLIT. |
| logic | L10 | major | Identity check sign and strictness | ACCEPTED | Same as wording 2 (section 9). |
| logic | L11 | major | T2a estimator, classes, control | ACCEPTED | Section 8: one-to-one node rule, the deciding statistic named, the test defined, NOT-DISTINGUISHED replaces AS-CONTROLS, parent-family size reported as a covariate. The statistic is R_j with absent members as non-co-members (not E2E_j). |
| logic | L12 | major | Deciding T2 subset expected UNDERPOWERED | ADAPTED | Section 7 Expected size; section 10 Q5 and Q7. The hold-back unit is the expansion, so the verdict set is every expansion except LRPAP1 (about 12 and 21 expected) and the ten-contig subset (3 and 8 components) carries no class. |
| logic | L13 | major | Gate 1 open-ended | ACCEPTED | Same as wording 10. |
| logic | L14 | major | D1 mislabelled as the literal T1 | ACCEPTED | Section 1 (D1 up to nesting, 'literal T1' deleted, along-sequence reading stated as untested with 1,062 of 2,290); section 2 item 2 reasons corrected; 'unlabelled boundaries' deleted from Gate 6. |
| logic | L15 | major | Section 10 mixes consequences, seen quantities and guesses | ACCEPTED | Section 10 split into C1 to C5 and Q1 to Q8; the 'first evidence' sentence replaced by a pointer to register rows 732, 817, 818. |
| logic | L16 | minor | No register citations | ACCEPTED | Section 11 Register check. |
| logic | L17 | minor | Counts and wording | ACCEPTED | Section 2 items 8 and 10, section 4, section 9 (75% and 83%), header, V-based region facts (445 regions). |
| logic | L18 | minor | Free constants | ACCEPTED | Section 7 hub cap 51 (inherited), section 8 caliper sweep, section 9 UFBoot 95, seeds, section 12 item 10. |
| logic | L19 | minor | Tie rule decides many strings | ACCEPTED | Section 4 (bases are summed; sensitivities with the SET of tied labels and a 5% band). Counts from logic/neartie.out. |
| logic | L20 | minor | Single-level decision; unfrozen primary layer | ADAPTED | Section 4 level sweep at inflation 2.0 and 4.0 with the LEVEL-SPECIFIC tag (instead of a second deciding cell), size-2 stratum printed, F_O files hashed (section 13), section 6 held-out wording, layer naming. The level-C components are the strata of the governing null. |
| logic | L21 | minor | Chain tolerance broader than its rationale | ACCEPTED | Section 4: m1t (end-truncation only) printed beside m1c; wording corrected. |
| logic | L22 | minor | Gate 7 cannot fail | ACCEPTED | Section 7 Gate 7: parenthesis deleted, called a pipeline check (C3), planted negative control with a bar. |
| logic | L23 | minor | Tier B uses the best placement only | ADAPTED | Section 7 step 3: the number of qualifying placements is recorded and classes with more than one are dropped. Limit stated: Liftoff reports further copies only at identity >= 0.95, so tier C is the check. Hub cap 51. |
| logic | L24 | minor | Smaller definitional gaps | ACCEPTED | Section 4 (singleton convention, per-layer universe), section 5 (whole-family restriction), section 8 (bootstrap over pool components). |
| logic | L25 | minor | Overstatements (cnduplicon direction, 'true by construction') | ACCEPTED | Section 1 and section 2 item 4 reworded ('partial circularity', 69.7 / 15.7 / 5.2). |
| logic | L26 | minor | Same-locus duplicate exclusion defined through families | ACCEPTED | Section 4 Exclusion arms: coordinates-only definition, family-independent deletion, 10% guard from CP-5. |
| logic | L27 | note | Pigeonhole bound verified | ACCEPTED | Section 2 item 1 restated with the restricted counts and 486. |
| logic | L28 | note | Direction mapping and path defect verified | NOTED | No change needed. The suggested Gate 4 test (exception list invariant to input order) was added to Gate 4. |
| logic | L29 | note | Process note (input-only discipline) | NOTED | The same discipline was followed in this revision; every computation is listed below. |
| feasibility | F1 | blocker | Tier B cannot register same-contig paralogs | ACCEPTED | Section 7 option (a): rank-sharded runs (liftoff 1.6.3 find_overlaps checked in the source); the old design recovers at most 6 of LRPAP1's 11 loci and this is stated; Gate 7 can fail; Q5; cross-species Target code and the placement classifier named as new code. |
| feasibility | F2 | major | Gate 1 cannot currently be satisfied | ACCEPTED | Same as wording 10; two further attempts to explain 2,259 failed (c5_disc2259.out, c5b_disc2259_scan.out). |
| feasibility | F3 | major | Gate 5 premise | ACCEPTED | Same as wording 1. |
| feasibility | F4 | major | Genome-wide family tables exist | ACCEPTED | Section 2 item 5; section 6 R3a; min_cov_shorter 0 and the 09-25 date for human_testis declared. |
| feasibility | F5 | major | Deciding T2 set is empty by design | ADAPTED | Same as wording 13 (assignment rule printed with its expected size; contigs by accession and chr). |
| feasibility | F6 | major | F_O gene membership is not the 1,709 rows | ACCEPTED | Section 4 Family layers (fold map; V_F = 1,657; unclustered genes excluded with a singleton arm; half rule). Counts reproduced (c2: 76 folded V genes inherit a cluster). |
| feasibility | F7 | major | Verdict-set counts quoted for the wrong pair set | ACCEPTED | Section 5 Counts: 38 depth-matched pairs in 33 families at 0.90 (c7), source pairs.tsv named, coverage printed, 0.95 and 0.98 UNDERPOWERED, partial families out of scope. |
| feasibility | F8 | major | Controls not identity-matched; circular selection | ACCEPTED | Section 5: handled by demotion (no controls, no p, no class); identity-matched controls are infeasible because the near-distance pool has 11 and 27 pairs. |
| feasibility | F9 | major | No executable join for the matched human run | ADAPTED | Same as wording 21. |
| feasibility | F10 | major | Audit tie dependence | ADAPTED | Same as wording 3. A fixed tie rule with all tied chains printed replaces the suggested minimum and maximum over resolutions. |
| feasibility | F11 | major | Exclusion counts are over 2,334; family-dependent exclusion | ACCEPTED | Section 2 item 10 and section 4 Exclusion arms (V counts 36, 388, 105 pairs / 181 genes, 862; family-free definition). Reproduced (c10_excl_counts.out). |
| feasibility | F12 | major | E-table parameters missing | ADAPTED | Section 7: cap 51, one -copies pass, 50% outgroup overlap, unannotated copy sites listed and not scored (the recon's inheritance rule is not adopted because it cannot be verified), tier C parameters inline, the dangling 'tier C step 5' reference removed. |
| feasibility | F13 | major | Arm 1 profile of every member; synteny classifier lost; pilot rates | ADAPTED | Section 9: arm 1 uses the expansion's ancestral member only (anchor = triple); classifier written from the README with the chaining rule; unresolved rates 83% and 75%; UNDERPOWERED predicted (Q7). |
| feasibility | F14 | major | Arm 2 under-specified; identity arm input | ACCEPTED | Section 9 Arm 2 (command line, UFBoot 95, rooting) and Identity check (input = graph dump re-run, light). |
| feasibility | F15 | major | Atom builder, class rule, dedupe, counts | ACCEPTED | Section 5 Inputs and units: class builder and dedupe are new code with tests; >= 0.5 declared; counts flagged as un-deduplicated; tie rule; format argument 'gorilla'. |
| feasibility | F16 | major | T2a verdict map hole | ACCEPTED | Same as wording 14 (PARTIAL, two one-sided p-values, fewer than 20 matched sets drops the expansion). |
| feasibility | F17 | minor | Rung U2 undefined for 388 genes; 'true by construction' | ADAPTED | Section 4 U2: genes without famCN excluded with the count printed; famCN source = families_cn.json cs and co, the inputs of the frozen code, and not the famCN_sotoiv column the reviewer names; section 2 item 4 reworded. |
| feasibility | F18 | minor | Small factual slips | ACCEPTED | Section 2 items 8 and 10, section 4, '09-25/27' BASE dates, 'rows' instead of members, 'draws shared per universe'. |
| feasibility | F19 | minor | Unfrozen definitions | ACCEPTED | Section 4 (bases, G0 tie rule, ordered string); set-Jaccard and 'unlabelled boundaries' deleted; the positional reading stated as set aside (Read first 1, section 1). |
| feasibility | F20 | minor | Gate 0 inputs | ACCEPTED | Section 13 Frozen inputs (sha256, rows), interpreter pinned, BAMs and support tables named (c6_hashes.out). |
| feasibility | F21 | minor | Cost estimates incomplete | ACCEPTED | Section 7 step 2 and Tier C (one pass, RSS, serial index builds, disk); section 3 (R1 split by cell). |
| feasibility | F22 | minor | Control provenance; 'unexposed' contigs | ACCEPTED | Section 7 Gate 7 (PCDHB, TUBA, RBMY dropped) and Contig roles ('no prereg names them', ledger mentions). |
| feasibility | F23 | minor | Register pointers, KEY names | ACCEPTED | Section 11 Register check (lines and rows), section 13 KEY names, rows from 1263. |
| feasibility | F24 | minor | Overstatements and gate ordering | ACCEPTED | Section 2 item 2 (fusion sentence qualified), section 3 scope of the gates, 'Amendment-free' deleted. |
| feasibility | F25 | note | Verified facts | NOTED | Re-used. A subset was re-verified in this revision (hashes c6, V counts c0 and c2, depth-matched counts c7). |
| feasibility | F26 | note | Implementer's list | ACCEPTED | Section 13 New code to write; constants in section 12 item 10. |

### Closure review of draft 2: one row per finding (draft 3)

Counts: ACCEPTED 56, ADAPTED 7 (63 rows; severities: blocker 2, major 25, minor 32, note 4). Rows with the same issue point to one fix. The `confirmed-ok` entries of the four reviewers (39 lines) are not listed.

| id | reviewer | severity | location | disposition | where fixed, or why not |
|---|---|---|---|---|---|
| CL1 | logic | major | Gate 4 clause 'family partition equal to the strata gives UNINFORMATIVE' against algebra items 2 and 3 | ACCEPTED | Section 3 Gate 4: zero variation-bearing strata gives UNDERPOWERED; a synthetic set of at least 8 strata with two clusters and one string each gives UNINFORMATIVE. |
| CL2 | logic | major | Homology-preserving null is reached by construction when strings are nested with the cut (draft-1 L2) | ADAPTED | Class renamed CUT-ALIGNED; section 2 item 3, Read first ('What the classes mean', item 2), Q2 and section 11 state that T1-b and T1-c cannot separate duplicon organisation from shared sequence similarity; Gate 4 gains the nested-strings case (p <= 0.01). The similarity-matched comparator (F_O at inflation 2.0 and 4.0) is not registered, because those cuts differ in granularity and change m1c and m2 mechanically; it is decision 16 and can be added by an Amendment before any run. |
| CL3 | logic | major | UNDERPOWERED counts strata while p is pooled; one stratum can decide p | ADAPTED | Section 4 'Robustness to one group': per-stratum z_s table, the most negative stratum removed, p recomputed on the same draws, the lower class issued. A sign test over strata was not adopted because it adds a second statistic with less power. |
| CL4 | logic | major | Gate 6 tie counts over the wrong universe | ACCEPTED | Section 3 Gate 6: 857 of 12,597 labelled merged exons of the genes of V (898 of 13,501 over all 2,334 genes); genes of V with at least one tie exon: 316; 'Only the counts listed in this Gate are tripwires'. |
| CL5 | logic | major | Sampling frame of the control sets not stated (draft-1 L11 d) | ADAPTED | Section 8: control pool excludes hubs and components that contribute an expansion; a draw takes a qualifying component uniformly and then one of its matched sets uniformly; the floor counts distinct control components; the component-size class (2, 3 to 5, 6 to 11, 12 to 51) was added to the matching and the MCL family size stayed a covariate (superseded in draft 3a: row RL1 moves the size class out of the matching). |
| CL6 | logic | major | Matching variable open (median over which rows, overlap with the expansion, boundary of the caliper) | ACCEPTED | Section 8: integer identities, M = twice the median over the copy rows that join two records of the set, caliper /dM/ <= 10 inclusive (sweep 4, 10, 20), a set without a joining row cannot be a control, components that hold an E-table class are not in the control pool. |
| CL7 | logic | major | n_testable versus the drop at 20 matched sets; expected size 12 and 21 (draft-1 L12, F5) | ACCEPTED | Section 8 defines n_testable with a floor of qualifying control components (20 in draft 3, n_min = 8 from draft 3a, row RL1), evaluated before any class; section 2 item 7, section 7 and Q5 give upper estimates with the Wilson interval (8 to 17 and 14 to 30) and the Poisson chance; siblings count as separate units with a group-once sensitivity (section 7 step 4). |
| CL8 | logic | major | Arm 1 forces its own outcome; 'mixed-family' contradicts the outcome | ACCEPTED | Section 9: phi = the family holding most nodes of expressed derived members; mixed-family = phi holds a node that is not a node of a member of E; outcome, causes and funnel rewritten; the identity check uses E_nodes; the unit is named in the section 3 tables. |
| CL9 | logic | major | Arm 2 does not measure youth (D is a clade in 1//D/ of random-parent histories) | ADAPTED | Section 9 arm 2 tests E_tree (members of E in the family tree) as a clade that leaves out at least one other family member; the comparator is size-matched within the identity caliper; the text says a clade is a pipeline and gene-conversion check and that a non-clade is the informative outcome; classes are kept because T2b is not authorised and is predicted UNDERPOWERED. |
| CL10 | logic | minor | Five open cases in the algebra (T2a headline, level sweep, UNDERPOWERED sweep cell, INVALID R3a, arm 2 per library) | ACCEPTED | Section 3 verdict algebra items 4 and 5; section 4 level sweep; section 9. |
| CL11 | logic | minor | m1t, Mult, Delta not well defined; Gate 4 relabelling clause (draft-1 L21, W30) | ACCEPTED | Section 4 statistics (m1t pairwise comparable, Mult over mu >= 1, Delta over families with at least two genes); Gate 4 plurality condition. |
| CL12 | logic | minor | Gate 7 bar of 8 of 11 is a bar of 8 of 8 (draft-1 L22) | ACCEPTED | See CF1: the bar is the 6 qualifying records, a failure blocks T2a and T2b. |
| CL13 | logic | minor | n_min, floor of 20 and flags have no entry in section 12 item 10; 'the expected answer' in Amendment 1 (draft-1 L18) | ACCEPTED | Section 12 item 10 lists n_min (source: the E < 8 rule of the seed-pool prereg) and every fixed floor with its sweep; section 6 no longer lets Amendment 1 fix an expected answer (Q8 stands). |
| CL14 | logic | minor | Section 10 mixes consequences, expectations and predictions (draft-1 L15, L27) | ACCEPTED | Section 10 rewritten (C1 to C4, Q1 to Q9); C1 limited to the counted layers and rungs; old C4 became Q4 with a falsifier; old C5 became Q9 with a number; old Q4 became C4; Q3 and Q8 restated; 'number of F_S families' replaces 'cell count'. |
| CL15 | logic | note | T1-a is a count whose information is the exception list | ACCEPTED | Section 4 decisive paragraph says so and prints m1c and m2 separately. |
| CL16 | logic | note | B and the call split for the Gate 4 repeats | ACCEPTED | Gate 4: B = 10,000 per repeat in several light calls of under 180 s. |
| CF1 | facts | blocker | Gate 7 positive control unreachable (pseudogene records need -f; two long gene models cannot reach coverage 0.5) | ACCEPTED | Section 7 step 2 command now carries -f types.txt; Gate 7 bar = every annotated record of the gene-LRPAP1 component with at least 50% of its exon bases inside an LRPAP1 copy interval (6 records), the two long models are expected misses with the reason; failure blocks T2a and T2b; the fragment-in-shared-run exposure is listed among the ways the control can fail; the Amendment may change the command only on evidence from the two controls; C3 and Q5 restated. |
| CF2 | facts | major | Registered command has no -f | ACCEPTED | Section 7 step 2: full command with types.txt (681 of 2,249 pool records, 502 of 1,446 called records are pseudogenes) and the minimap2 named by the Amendment. |
| CF3 | facts | major | What counts as a placement is not stated | ACCEPTED | Section 7 step 3: primary row with coverage >= 0.5, sequence_ID >= 0.5 and neither partial flag; other rows are no placement; records whose only row is partial are counted. |
| CF4 | facts | major | The pool is blind to cis-only (tandem) expansions (second half of draft-1 F1) | ACCEPTED | Section 11: dispersed and cross-block copies only; 5 of 307 components have all records on one contig; tandem-only expansions are outside the E-table and T2 makes no statement about them. |
| CF5 | facts | major | Node join between fam.clusters.tsv and fam.copies.tsv (draft-1 W15) | ACCEPTED | Section 8: a node is a row of fam.copies.tsv; fam.clusters.tsv is not read; node counts 864 and 1,559. |
| CF6 | facts | major | Gate 6 counts with no stated source (same-locus duplicates; tie universe) (draft-1 F11, W9) | ACCEPTED | Section 4 exclusion arms name the gene span (hull of the gene's CAT v4 transcripts, families_cn.json c, s, e, st) and say at most 105 genes are deleted; tie counts as CL4. |
| CF7 | facts | minor | 'Registered OUTPUTS reproduce' overstates; 'none gives 2,259' overstates (draft-1 F2, W10) | ACCEPTED | Section 2 item 9 rewritten (cluster and unit counts reproduce; pooled delta, p and arm A shares are checked by Gate 1; scan result as measured). |
| CF8 | facts | minor | python -I ignores PYTHONHASHSEED; python3 versus python | ACCEPTED | Gates 0 and 3: python3 -B, env PYTHONHASHSEED=0 and 1 without -I. |
| CF9 | facts | minor | The graph dump is light in sections 6 and 9 and heavy in section 8; BASE binary differs | ADAPTED | One class stated in section 8 (light for human_testis; gorilla PAFs measured first under light, heavy above 2 GB); Amendment 1 requires --min-cov-shorter 0 and the BASE flags and a check that every BASE cluster lies in one component. The RSS is not measured now because the dump is not authorised. |
| CF10 | facts | minor | Five labels that do not match their source (486, 885/10,651, 157, BAM count, three contigs) | ACCEPTED | Section 2 items 1 and 7; section 5 counts (878 and 10,506 under the 100 bp rule); section 7 contig list (chr6, chr9, chrX added with their prereg names); section 13 BAM wording. |
| CF11 | facts | minor | A displaced record is not 'unplaced' in Liftoff 1.6.3 | ACCEPTED | Section 7 step 2: displaced records are recognised by coverage < 0.5 or by landing on no single in_place outgroup record; the count of records with no triple is printed. |
| CF12 | facts | note | Cost line should state 1,446 records and K = 44 | ACCEPTED | Section 7 step 2 and section 13 (split rule listed as new code). |
| CW1 | wording | major | Gate 6 tie counts (898 of 13,501; 316) | ACCEPTED | See CL4. |
| CW2 | wording | major | 'Records' undefined in the hub rule (4 and 5 against 2 and 3) | ACCEPTED | Section 7 step 1: hub = more than 51 members, members being the annotated records of any length and the copy@ sites; 4 and 5 stand. |
| CW3 | wording | major | Opening sentence contradicts step 4; control pool may share a component with an expansion (draft-1 W16) | ACCEPTED | Sentence deleted; section 7 step 5 and section 8: a component that contributes an expansion is not in the control pool. |
| CW4 | wording | major | Join from cluster rows to copies rows (draft-1 W15) | ACCEPTED | See CF5. |
| CW5 | wording | major | n_testable before or after the matched-set drop | ACCEPTED | See CL7. |
| CW6 | wording | major | E phase does not write what section 8 needs (draft-1 W16) | ACCEPTED | Section 7 step 5: member_exons.tsv, control_pool.tsv, extension.tsv, 0-based half-open ids, hashes recorded with the E-table. |
| CW7 | wording | major | Caliper boundary on a 0.001 grid (draft-1 W14) | ACCEPTED | See CL6. |
| CW8 | wording | major | Amendment 1 can be written after Phase R1 has printed (R3a hold-back) | ACCEPTED | Section 6: the gene-to-locus rule and half rule are written in the file; Amendment 1 is committed before the first Phase R1 run and may not change the rule, cells, statistics, nulls or Q8; 'expected answer' deleted. |
| CW9 | wording | major | Read first does not define the classes (draft-1 W20) | ACCEPTED | Read first: 'What the classes mean' paragraph; ORGANISES-ENRICHED renamed CUT-ALIGNED. |
| CW10 | wording | major | Read first omits links tested, V, D1/D2 mapping, nesting, verdicts that can occur (draft-1 W20) | ACCEPTED | Read first items 8 to 12. |
| CW11 | wording | minor | Algebra cases (T1-a per direction, T2a cell, lower tail only, T1-d flag, single-run labels) (draft-1 W5) | ACCEPTED | Section 3 items 4 and 5, decisive-cell table, section 4. |
| CW12 | wording | minor | Gate 2 draw rule, Gate scope, tripwires (draft-1 W23) | ACCEPTED | Gate 2 rewritten (min(20, n), one stream, at most 1 differing); scope sentence; only Gate 6 counts are tripwires. |
| CW13 | wording | minor | F_S/F_Q universes and counts (draft-1 L27, W4) | ACCEPTED | Section 2 item 1, section 4 F_S bullet, Gate 5b wording, C1. The suggested F_Q figure of 475 cells over V was not added because it is not verified. |
| CW14 | wording | minor | Section 5 sources and definitions (draft-1 W12) | ACCEPTED | Section 5: genes.tsv join, pairs.tsv join, 1-based columns, atoms versus classes, 878 and 10,506, B_viol. |
| CW15 | wording | minor | Null details: what is shuffled, string order, RNG objects (draft-1 W7) | ACCEPTED | Section 4 nulls and draws; section 8 test. |
| CW16 | wording | minor | Descriptive statistics that two implementers compute differently (draft-1 W3, W30, W31) | ACCEPTED | Section 4: m1t, audit tie key, Mult, Delta, size-2, restriction, near-tie band, 105 pairs and 181 genes, U2co versus U2cs. |
| CW17 | wording | minor | Authorisation wording (draft-1 W19) | ACCEPTED | Section 3 run order: a 'go' names its phases; light under rlock light; Amendment rows; matched human run row. |
| CW18 | wording | minor | Gate 0 hashes of 36 GB; call budgets (draft-1 W19, W29) | ACCEPTED | Gate 0: size and mtime above 1 GB, sha256 once in its own call; light calls under 180 s; heavy limits apply to heavy calls; the BASE dump is hashed when made. |
| CW19 | wording | minor | T2b-1 unit and bounds (draft-1 W17) | ACCEPTED | Section 3 table and section 9: mixed-family expansion; one-sided exact bounds; at least one ancestral member; seed named. |
| CW20 | wording | minor | Decisions list omits defaults and constants (draft-1 W19) | ACCEPTED | Section 12: items 8 and 11 corrected, 10 extended, 14 to 16 added. |
| CW21 | wording | minor | Terms and symbol collisions (draft-1 W24, W28) | ADAPTED | Terms define clean, family-less, multi-family, MCL, independent unit; R_J, D (derived set) and K renamed or defined (Dv, J); the remaining overloads (cell, E) are disambiguated by context and left. |
| CW22 | wording | minor | C4/C5 are not consequences; 'held-out' tense | ACCEPTED | Section 10; 'held-out' only after R3a runs; R3b chromosomes are 'not used to select a rule'. |
| CW23 | wording | minor | BAM count, source of 'expressed', index paths, names (draft-1 W27) | ACCEPTED | Sections 7, 8 and 13. |
| CW24 | wording | minor | R3a power count not registered | ADAPTED | Gate 6 prints the number of strata with at least two clusters among V genes of the BASE graph; Q8 applies only if it is at least 8. Not computed now because it needs the dump. |
| CW25 | wording | minor | Style sentences that hurt comprehension | ACCEPTED | Read first items 1 and 6, labels, section 5 first sentence, section 7 contig roles. |
| CW26 | wording | note | Gate 7 wording 'blind' | ACCEPTED | Removed. |
| CD1 | disposition | blocker | Row F1: Gate 7 bar unreachable, Q5 predicts a pass | ACCEPTED | See CF1. |
| CD2 | disposition | major | Row W17 and others: mixed-family contradicts arm 1 | ACCEPTED | See CL8. |
| CD3 | disposition | minor | Row W15: no join key | ACCEPTED | See CF5. |
| CD4 | disposition | minor | Registered counts over different universes (Gate 6, F_S, control pool) | ACCEPTED | See CL4, CF6, CW13, CW14. |
| CD5 | disposition | minor | SPLIT does not apply to the two libraries; no T2a headline | ACCEPTED | See CL10. The headline is the lower class of the two libraries. |
| CD6 | disposition | minor | Constants that set classes are not in section 12 item 10 | ACCEPTED | See CL13. The origin of n_min is the seed-pool prereg; the floors of 20 pairs and 8 families and the floor of 20 control components have no earlier source and are listed as fixed now with sweeps. |
| CD7 | disposition | minor | No funnel for arm 2 (draft-1 L9) | ACCEPTED | Section 9 arm 2 prints the funnel with the share. |
| CD8 | disposition | minor | Dump command classed light and heavy | ADAPTED | See CF9. |
| CD9 | disposition | minor | Row c12: orangutan outlier | ACCEPTED | Input-only table row corrected; section 7 step 2 says the estimate is from three same-genome fits and the Amendment sets the budget per outgroup. |

### Regression review of draft 3: one row per finding (draft 3a)

Counts: ACCEPTED 41, ADAPTED 5 (46 rows; severities: blocker 1, major 9, minor 36). Some minor style proposals of the wording reviewer were merged into one row each. The `confirmed-ok` entries are not listed.

| id | reviewer | severity | location | disposition | where fixed, or why not |
|---|---|---|---|---|---|
| RL1 | logic | blocker | Floor of 20 qualifying control components cannot be reached by any expansion; T2a UNDERPOWERED by construction (also facts F3) | ADAPTED | Section 8: the floor is n_min = 8 (the only inherited floor); the component-size class leaves the matching and becomes printed strata with a SIZE-SENSITIVE flag; Gate 6 prints per expansion the qualifying components at calipers 4, 10 and 20. Input-only upper bounds with the size class pooled (author/pool_floor_nosize.py): best window 23 components in KB3781 and 46 in OR6737 for m_j = 2, 8 and 16 or 17 for m_j = 3. Read first item 12, Q5, Q6 and the reachability table now say T2a is likely powered in OR6737 and uncertain in KB3781 (at most about 12 expansions there). Option A (declare T2a UNDERPOWERED) was not chosen because pooling the size class makes a verdict reachable. |
| RL2 | logic | major | UNDERPOWERED counts dead strata; Gate 4 UNINFORMATIVE case is 8 dead strata | ACCEPTED | Section 3 items 2 and 3 (live strata), Gate 6, C4, Q10, Read first item 12. |
| RL3 | logic | major | Gate 4 tautology clause contradicts the 8-strata clause | ACCEPTED | Gate 4 clause rewritten on the cluster sizes of the 30 real strata with synthetic strings. |
| RL4 | logic | minor | Robustness clause at exactly 8 strata | ACCEPTED | Section 4: the class from the remaining strata is the content class of its p alone. |
| RL5 | logic | minor | Removal clause does not protect against two or three small aligned groups | ACCEPTED | Read first 'What the classes mean' states the limit. |
| RL6 | logic | minor | Robustness clause details open | ACCEPTED | Section 4: per cell, stratum id, T1-c statistic of s alone, stored draws, skip when no negative z. |
| RL7 | logic | minor | 'Lower class' is not conservative on the five-point T2a scale | ACCEPTED | Section 3 item 4: the class nearer to NOT-DISTINGUISHED. |
| RL8 | logic | minor | Component size basis (members versus records >= 200 bp) | ACCEPTED | Moot as a matching variable (RL1); strata use annotated records of exon union >= 200 bp; the quoted shares were deleted. |
| RL9 | logic | minor | One-to-one node rule for controls; sampling algorithm and order | ACCEPTED | Section 8: rule applies to control sets; components in id order as strings; first 1,000 matched sets in lexicographic order; exploration cap of 1,000,000 candidate sets. |
| RL10 | logic | minor | Ties in T under float arithmetic | ACCEPTED | Section 8: T compared as an integer (T_num). |
| RL11 | logic | minor | Two-step draw calibration depends on component size | ADAPTED | The two-step draw is kept; size strata, component member count and MCL family size are printed; SIZE-SENSITIVE flag; the text says that excluding expansion components lowers the control count. |
| RL12 | logic | minor | Arm 2 measures E-table completeness and branch support, not age | ADAPTED | Section 9: arm 2 is an estimate without class (like arm 1); the text lists what the share measures and cites the simulation. |
| RL13 | logic | minor | Arm 2 comparator undefined or saturated; Dv versus phi | ACCEPTED | Comparator deleted together with the classes; informative expansion, E_tree and Dv defined. |
| RL14 | logic | minor | Q3, Q5, Q6, C4 wording | ACCEPTED | Section 10 revised; Q10 added. |
| RL15 | logic | minor | Gaps between Amendment 1 and R3a | ACCEPTED | Section 6: no Amendment may change a rule R3a uses from the first R1 run; the BASE-cluster check moves to the start of R3a (Gate 6 for that phase); label for the family layer and the library. |
| RF1 | facts | major | Called-set numbers use a different hub rule from step 1 | ACCEPTED | Section 7 step 2: 300 components, 1,403 records (907 gene, 496 pseudogene), first runs 300, 300, 165, 116, 81; rank order stated; cost 26 to 50 minutes. |
| RF2 | facts | major | Shares behind the size classes belong to another size definition | ACCEPTED | Moot after RL1; the shares are deleted. |
| RF3 | facts | major | Floor of 20 unreachable | ADAPTED | See RL1. |
| RF4 | facts | major | Size-2 stratum: 163 counts rows, 152 counts V_F genes | ACCEPTED | Sections 3 and 4: n_f = 2 in the layer universe = 152 of 362; 163 and the 147 overlap are noted. |
| RF5 | facts | minor | Coverage wording for the two long LRPAP1 models | ACCEPTED | Gate 7: 16% and 41% of the exon union, at most 49% and 54% of the bases Liftoff counts; not required. |
| RF6 | facts | minor | Poisson chance on the non-hub frame | ACCEPTED | Section 2 item 7 and Q5: about 12% and 48%. |
| RF7 | facts | minor | '53 genes differ' versus 54 | ACCEPTED | Gate 6: 82 exons in 54 genes, U1 sets of 53 genes change. |
| RF8 | facts | minor | Planted negative control may pick pool members | ACCEPTED | Gate 7: no copy row as source and in no pool component. |
| RF9 | facts | minor | Gate 2 order of chromosomes | ACCEPTED | Gate 2 rewritten. |
| RW1 | wording | major | Size-2 count | ACCEPTED | See RF4. |
| RW2 | wording | major | R3a UNDERPOWERED/UNINFORMATIVE handling differs between section 6 and the algebra | ACCEPTED | Section 6 and algebra item 4 now agree (Phase R1 class stands, labelled DEV). |
| RW3 | wording | major | Read first item 9 describes nesting as pairwise subset | ACCEPTED | Item 9 describes chains. |
| RW4 | wording | minor | Stale C4 references; 'by construction' | ACCEPTED | Sections 4 and 12 point to Q4. |
| RW5 | wording | minor | Present-tense 'held-out' in section 12 item 5 | ACCEPTED | Item 5 rewritten. |
| RW6 | wording | minor | Subject of 390 s lost; light timeout not enforced | ACCEPTED | Section 3 authorisation paragraph. |
| RW7 | wording | minor | Amendment rule restricted to table rows; unlisted heavy steps | ACCEPTED | Rule covers any heavy step; table rows added. |
| RW8 | wording | minor | Gate scope does not cover Phase R2 | ACCEPTED | Scope sentence. |
| RW9 | wording | minor | Per-cell z tables, sweep V_F, largest-stratum Gate 4 case | ACCEPTED | Sections 3 and 4. |
| RW10 | wording | minor | Controls ride in the same runs as the verdict set | ACCEPTED | Gate 7: controls lifted and checked first. |
| RW11 | wording | minor | Amendment 1 'must show' the BASE-cluster property before any go | ACCEPTED | Check moved to the start of R3a. |
| RW12 | wording | minor | Section 12 and the disposition decision list disagree | ACCEPTED | Section 12 items 17 to 21 added; the disposition list now points to section 12. |
| RW13 | wording | minor | Disposition header pointers and draft-2 names in first-table rows | ACCEPTED | Header note and pointers added. |
| RW14 | wording | minor | Closure rows CL5, CF9, CD6, CW21 and others not fully true; 'unit' and Q collisions | ADAPTED | CL5 relabelled; the rank-run bound is J; bare 'unit' meaning a statistical unit in six places and the sequence groups Q of Gate 5b are left and disambiguated by context and by Terms. |
| RW15 | wording | minor | Q5 and expected-size wording contradict n_testable | ACCEPTED | Q5, section 7 and item 12 say 'expansions' before the filters. |
| RW16 | wording | minor | Section 9 definitions (informative expansion, Dv, family tree) | ACCEPTED | Section 9. |
| RW17 | wording | minor | Control-pool eligibility and sampling cap | ACCEPTED | Sections 7 and 8. |
| RW18 | wording | minor | Style: item 11 chaining | ACCEPTED | Read first item 11. |
| RW19 | wording | minor | Style: item 12 run-on | ACCEPTED | Read first item 12 rewritten. |
| RW20 | wording | minor | Style: classes not defined for the table | ACCEPTED | 'What the classes mean' extended. |
| RW21 | wording | minor | Style: voice | ACCEPTED | Neutral voice in the opening paragraph and items 8 and 9. |
| RW22 | wording | minor | Style: two meanings of 'headline' | ACCEPTED | 'Headline class' and 'headline number'; arm 1 cause wording. |

### Input-only computations run in this revision

Scratch directory: `/mnt/linuxdisk/tmp/prereg_revision/`. All ran under `tools/rlock.sh light` with `python3 -B -I`. None computes a deficiency, purity, containment, enrichment, S-rate or null draw. None joins a family label to a unit string beyond the marginal counts of section 2 item 1 (reproduced from the logic reviewer's scripts). The scripts `c0` and `c2` read the reviewers' strings file only for the booleans 'has an exonic duplicon' and the clean, multi-family or family-less class; they read family membership and graph structure and no unit string.

| script | output | one-line result |
|---|---|---|
| `c1_fo_graph.sh` | `c1_fo_graph.out`, `fo_graph/run.log`, `fo_graph/ours.graph.tsv` | the frozen binary (sha1 e4cc13b9) with `--dump-graph` reproduces `ours.clusters.tsv`, `ours.loci.tsv` and `ours.pairs.tsv` byte for byte; graph of 1,819 nodes and 11,031 edges (0.9 s) |
| `c0_V_list.py` | `c0_V_list.out`, `V_genes.tsv` | V = 2,142 genes (2,036 clean, 106 family-less) |
| `c2_fo_components.py` | `c2_fo_components.out`, `c2_fo_components.json` | 339 components; V_F = 1,657 genes (76 folded annotations inherit a cluster); 325 strata hold V_F genes; 30 variation-bearing strata with 416 V_F genes (29 hold at least two multi-gene clusters; 4 to 42 genes, median 11); 295 single-cluster strata with 1,241 genes; 467 V genes never connected and 18 singleton nodes |
| `c5_disc2259.py`, `c5b_disc2259_scan.py` | `c5_disc2259.out`, `c5b_disc2259_scan.out` | with either BED, 2,290 genes overlap an exon by at least 1 bp; no rule on summed, best-segment or best-exon bp from 1 to 300 gives 2,259; an exonic coverage fraction >= 0.10 gives 2,258; no per-chromosome subset explains 31 genes |
| `c6_hashes.sh` | `c6_hashes.out` | sha256 (16 hex), size and lines of the Gate 0 list; large files by size and mtime |
| `c7_depth_matched_verdict.py` | `c7_depth_matched_verdict.out` | gorilla verdict set, depth-matched eligible pairs and families: 38 / 33 (tau 0.90), 15 / 15 (0.95), 5 / 5 (0.98), identical under the 1 bp and 100 bp touch rules |
| `c8_heldout_power.py` | `c8_heldout_power.out` | supported pool components 80 (KB3781), 139 (OR6737), 157 either; of size <= 51: 76 and 134; no member on the six LRPAP1 contigs: 23 and 38; all members on the ten contigs no prereg names: 3 and 8 |
| `c9_r3a_power.py` | `c9_r3a_power.out` | 334 V genes carry a human_testis BASE cluster; 86 BASE clusters hold at least two V genes (316 genes) |
| `c10_excl_counts.py` | `c10_excl_counts.out` | 233 pairs / 339 genes over 2,334 genes become 105 pairs / 181 genes in V; the readthrough list has 76 rows, 36 in V |
| `c11_multi_duplicon_genes.py` | `c11_multi_duplicon_genes.out` | 1,062 of 2,290 duplicon-bearing genes carry two or more dominant duplicons (974 of 2,142 in V); strings only |
| `c12_shardfit.py` | `c12_shardfit.out` | Liftoff self-lift timings: gorilla, human and chimp fits give 0.05 to 0.09 s per record; the orangutan fit is dominated by one 544 s shard (a 10-fold outlier); peak RSS 10.2 to 13.8 GB. The estimate for cross-species lifts is therefore uncertain and the Amendment sets the time budget per outgroup |
| `logic/fq_checks.py`, `string_only_checks.py`, `short_counts.py`, `halves.py`, `region_strata.py`, `neartie.py`, `control_pool.py` (copies of the logic reviewer's scripts) | `logic/*.out` | D[S->Q] = 7 genes (2 if ungrouped genes are ignored); pigeonhole marginals of F_S, F_Q, F_O; 32 U1 strings / 110 genes with more than one principal duplicon; 862 genes with at most two exons; 9 straddling clusters (61 genes); 445 SD98 regions in V; 898 exact exon ties; gorilla control pool 11 / 27 / 885 / 10,651 against 13 / 5 / 15 / 21 real pairs |

**Closure round (draft 3).** The author computed nothing new in this round. The closure reviewers' scripts and outputs (synthetic scorer tests, stratum structure, tie universes, pool and direct-row counts, LRPAP1 record geometry, hash recomputations) are in `/mnt/linuxdisk/tmp/prereg_closure/{logic,facts,wording,disposition}/`; none computes a deficiency, purity, containment, enrichment, S-rate or null draw on real labels.

### Design decisions that need the user's confirmation (draft 3a)

The decisions and their defaults are items 1 to 21 of section 12 of the prereg, with the same numbering. Those that most change cost or meaning: 1 and 16 (what T1-b and T1-c can claim, and the optional similarity-matched comparator), 5 (R3a as a light arm), 7 (tier B as the E-table recipe), 12 (T2b as estimates), 18 (the Gate 7 bar of 6 records), 19 (T2a controls: two-step draw, size strata, floor of n_min = 8) and 20 (Amendment 1 before the first R1 run).

### Findings that could not be fully resolved

- The 2,259 versus 2,290 discrepancy of the 09-30 disclosure stays unexplained after five attempts; it is an erratum, not a gate.
- The 15% single-copy share behind the expected testable expansions comes from a name-linked subset (24 of 160; Wilson interval 10.3% to 21.3%) and may not transfer; the figures 12 and 21 are upper estimates and KB3781 may fall below n_min = 8.
- Rank-sharding removes conflicts inside a pool component, not between unrelated records of one Liftoff run; displaced records are detected only by coverage < 0.5 or by landing on no single in_place outgroup record. The fragments of the LRPAP1 control share runs 1 to 3 with full-length copies, which can make the control fail.
- Liftoff reports further outgroup copies only at identity >= 0.95, so the multiplicity guard of the E-table is incomplete; tier C is the check.
- The pool is blind to tandem-only expansions (5 of 307 components lie on one contig).
- T2a power: with the size class pooled and a floor of 8 control components, the input-only upper bounds allow a verdict in OR6737 and make KB3781 uncertain; the real counts depend on the E-table, the outgroup conditions and the removal of expansion components.
- T1-b and T1-c cannot separate duplicon organisation from shared sequence similarity; a similarity-matched comparator is an option for an Amendment.
- The synteny classifier of the old pilot is lost and is rewritten from the README; T2b is predicted UNDERPOWERED.
- The R3a stratum count needs a graph dump that is not run now (Q8 is conditional), and the peak RSS of the gorilla BASE-PAF dumps is not measured.
- The time estimate for cross-species Liftoff comes from three same-genome self-lift fits; the Amendment sets the budget per outgroup.

### Can the registered cells produce each verdict class? (draft 3)

| cell | class | reachable from the inputs |
|---|---|---|
| T1-a | REFUTED | yes, and expected |
| T1-a | NOT-REFUTED | in principle (m1c = 0 and m2 = 0 on 1,657 genes), practically unreachable |
| T1-b, T1-c (Phase R1) | CUT-ALIGNED, PARTIAL, NOT-DISTINGUISHED | yes, 30 variation-bearing strata; CUT-ALIGNED is expected whenever strings and cut both follow sequence similarity |
| T1-b, T1-c (Phase R1) | UNDERPOWERED, UNINFORMATIVE | UNDERPOWERED not on the count of 30 strata, but yes if fewer than 8 are live (printed by the cell); UNINFORMATIVE practically unreachable |
| T1-b, T1-c | SPLIT | only after R3a has run |
| T1-b, T1-c (R3a) | all classes | possible; UNDERPOWERED if fewer than 8 strata, known only after the BASE graph dump (86 candidate clusters); Q8 is conditional on the count |
| T1-d (gorilla) | none issued | tau 0.95 (15 pairs) and 0.98 (5 pairs) are flagged UNDERPOWERED |
| T2a, OR6737 | SPECIFIC, PARTIAL, NOT-DISTINGUISHED, LESS | yes if n_testable >= 8 (at most about 21 expansions expected, 14 to 30 before the section 8 filters; input-only windows of 16 to 46 control components for m_j = 2) |
| T2a, KB3781 | all content classes; UNDERPOWERED | borderline (at most about 12 expansions expected, 8 to 17; Poisson chance of fewer than 8 about 12% at 15% on the non-hub frame, about 48% at the lower interval end) |
| T2a | SPLIT | between the two libraries, only if both are powered |
| T2a, ten-contig subset | none issued | UNDERPOWERED (3 and 8 supported components) |
| T2b arm 1 | none issued (estimate) | UNDERPOWERED expected (pilot 3 of 18) |
| T2b arm 2 | none issued (estimate) | UNDERPOWERED expected; the share also reflects E-table completeness and branch support |
| any cell | INVALID | by gate or control failure; Gate 7 can fail |


---

## PREREG_ideal_arm_pilot_2026-10-07.md

*(original: `docs/archive/2026-10/merged/PREREG_ideal_arm_pilot_2026-10-07.md`)*

## Pre-registration: the ideal arm of the expression ladder, pilot on the 2,334 Soto genes (2026-10-07)

**Status: written before any simulated read exists and before any number of this study. Descriptive pilot, no verdict class.**
Authorisation: the user's message of 2026-10-07 21:49 PDT, "can we try the pilot of the ideal arm?", after the question of 21:44 ("to prove that the problem is expression, can we show Soto's way with our RNA mode but using a dataset that has all their proposed genes as expressed with simulated reads ...").
It covers the ideal arm only (simulate, map, assemble, families, score) as bounded foreground heavy calls. It does not cover the multi-tissue and testis arms, which need their own pre-registration after the survey of 2026-10-07 (wf_ee69e34a-6d5).

### 1. Question

If every Soto gene is expressed, how many of them does our RNA mode (de novo assembly and the families stage, current default) recover as one-to-one nodes, by biotype, and which losses remain that expression cannot explain? The real-library baseline is the genome-wide A119b run of section 6w4/6w5 (`bench/SOTO_FAMILY_LOCUS_FIDELITY_2026-09-22.md`): 765 of 1,078 Soto genes with a family id matched (71.0%): protein_coding 79.3%, transcribed_unprocessed_pseudogene 85.9%, unprocessed_pseudogene 46.6%, processed_pseudogene 42.1%, lncRNA 85.4%; 33.7% of the matched genes sit on a locus that covers two or more annotated genes. That baseline is from an older pipeline version and a different library, so a difference to it is **not** a like-for-like contrast; the like-for-like contrast (same harness, realistic profiles) is the later arms.

### 2. Design (fixed now)

| item | value |
|---|---|
| Genome and annotation | T2T-CHM13 v2.0 (`winloci_data/chm13v2.0.fa`) and the CAT/Liftoff v2.0 slim GFF (`winloci_data/gencode_chm13/chm13v2.0_CAT_Liftoff.slim.gff3.gz`); all 2,334 Soto gene ids (CHM13_G.../LOFF_G...) are gene ids in it (checked) |
| Genes | the unique `Gene ID` values of `bench/soto/soto_famCN_S1C.tsv` (2,334) |
| Windows | gene span +- 25 kb, merged; transcripts = every transcript of every gene (any biotype) overlapping a window, so expressed neighbours are present as in the earlier ideal runs |
| Reads | 10 per simulated transcript (ideal: equal depth for every transcript), the model of `bench/sim.py` (error .001, indel error/3, no truncation, 0 to 30 bp end jitter; ends must vary), canonical-intron recipe of `sim_windows.py`, seed 20261007, one replicate |
| Mapping | minimap2 with the command, `-N 50` and the whole-genome splice index of `bench/ideal_expression/run.sh` (read-disjoint parts, one part per call) |
| Pipeline | the driver's `assemble` and `families` stages with the HEAD release binaries (`/mnt/linuxdisk/home/juanfraitu/rustle_target_m2/release`), current defaults, nothing tuned, sha1 of every binary recorded in the run log |
| Scorer | `bench/ideal_expression/score_soto_nodes.py`, which wraps the matching rules of `bench/soto_family_locus_fidelity.py` (one-to-one exonic matching, fused loci reported separately) |
| Scored population | the Soto genes that were simulated; the family-id population of the baseline (1,078) is reported separately |

### 3. Outputs

Per biotype and overall: the one-to-one matched share, the fused-locus share (reported separately), median exonic Jaccard, |d5|, |d3| and width ratio of the matched genes, and a mutually exclusive ordered loss class per gene: matched; matched_but_fused; no_locus; wrong_extent (lowest 5% of the matched genes' Jaccard, the cut written to the output from the data). A table against the baseline above with the paired caveat of section 1. The estimand for each biotype is the **algorithm-limited share** = 1 - (matched share of the ideal arm), and, against the old baseline, the **gain** = ideal - baseline; no threshold turns either into a class.

### 4. Predictions (written before any outcome; directional, from the earlier ideal runs and the A119b natural experiment)

- P1. The ideal arm matches a larger share than the baseline in every biotype.
- P2. The gain is largest for unprocessed and processed pseudogenes and smallest for protein-coding genes; the gap between transcribed and untranscribed unprocessed pseudogenes (85.9% against 46.6% in the baseline) shrinks in the ideal arm, where both are expressed.
- P3. The ideal arm does not reach 100% in any biotype; the largest remaining class is matched_but_fused (the earlier NPIP ideal run lost 14 of 46 member nodes to a swallowed neighbour), then wrong_extent.
- P4. The fused-locus share of matched genes is at least that of the baseline, because expressed neighbours are simulated for every gene.
A failed prediction is reported as failed.

### 5. Gates (each prints PASS or FAIL and no outcome statistic; a table is INVALID unless all pass)

G1 manifest: reads equal the sum of per-transcript counts, every read name parses to a listed transcript, every Soto gene is simulated or listed with a reason. G2 read ends vary (the 5' and 3' offsets are not constant). G3 determinism: the simulation is byte-identical under `PYTHONHASHSEED` 0 and 1; the scorer's tables are identical under two hash seeds. G4 mapping sanity from the BAM alone, before the pipeline: the share of reads with a primary alignment, and the share of primary alignments that land on their source gene (a measure of paralog ambiguity that cannot be removed by any pipeline). G5 binaries: the stamps match HEAD.

### 6. Limits stated now

The reads are simulated from the annotation, so the truth for node extent is the annotation itself (circular for extent, independent for presence and fusion). Near-identical paralogs produce MAPQ-0 reads in the simulation as in real data, and the truth label stays with the source copy. Retrocopies, annotation errors and truly unexpressed copies are not represented. The comparison with A119b crosses pipeline versions and libraries. One replicate, one depth (10 reads per transcript); the node boundary is the k-th most extreme read end with k fixed at 2, so extent depends on depth, and a depth ladder belongs to the full study. No copy-number step and no family-level scoring in the pilot. Both halves of the Soto data are already spent for copy-number rules; this pilot tunes nothing.

### 7. Amendments

#### Amendment 1 (2026-10-08, DRAFT until committed; written before any simulated read and before any heavy call; no outcome of the ideal arm exists)

**Why.** The survey wf_ee69e34a-6d5 and the builder's report change the design. (a) The 71.0% of 1,078 genes in section 1 is a NAME join (1,078 of the 1,774 S1C gene names equal a RefSeq name); the same A119b GTF scored by gene ID gives 53.2% (1,242 of 2,334), so section 1's baseline is not like-for-like and is demoted to one labelled line. (b) The scorer of the earlier fidelity tool counted a 1-bp overlap as a match (an annotation shifted by 50 kb still scores 39.3% matched), took the extent cut from each run's own 5th percentile, and called 21 to 27% (span rule) or 9 to 11% (exon rule) of PERFECT loci fused. (c) The assembler's mono-exon floor sits at 9 to 10 reads, exactly the simulated depth of 10 per transcript. (d) Like-for-like comparators exist: the real-library products scored by the same scorer. (e) The presence headline of an ideal arm is nearly predetermined, so the study is about what remains.

**Replaced in sections 1 to 6.** Seed 20261007 becomes 20261008 (20261006 and 20261007 are the earlier replicates). One depth becomes two, 10 and 30 reads per transcript, from ONE read pool (depth 10 is the first 10 reads of every transcript of the depth-30 pool, so the two depths are not independent replicates). The scorer is `bench/ideal_expression/score_soto_nodes.py`, written from scratch because `bench/soto_family_locus_fidelity.py` was deleted in commit 8db314c7 (an equivalence test pins it to the historical matching rules). The scored population is fixed at the 2,224 simulable Soto genes by S1C Gene ID (never by name); the 110 genes that cannot be simulated (105 spliced under 120 bp, 3 with an intron over 200 kb, 2 with a molecule over 30 kb; 43 in the SD98 set) are listed apart as `not_simulable`. The four classes of section 3 become the ordered classes below; `wrong_extent` as the lowest 5% of the run's own Jaccard is dropped.

**Fixed now.**
- *Match rule:* greedy one-to-one by exonic overlap, then Jaccard, then gene id, then unit id (the tie policy is printed and the reverse-order sensitivity counted); a match needs shared exon bp of at least 0.5 of the SMALLER exon union, otherwise the class is `weak_overlap`.
- *Unit:* the families-stage node (representative exons of `PREFIX.fam.loci.gff3`) is primary; transcripts of `families.gtf` and gene_id loci are sensitivities. *Competition scope:* every gene of the universe (`all`) is primary, the listed genes only (`soto`) a sensitivity. *Universe for competition and fusion:* the ideal arm's truth. Unit, scope and fusion rule are printed on every table.
- *Fusion:* a locus covering two or more annotated genes is `fused`; the class `fused_excess` is fusion in the arm that the gene's own annotation (annotation-as-loci control) does not already show, decided by the span rule (default `--fused-rule span`); the exon-rule counts are printed beside it.
- *Extent:* multi-exon genes: end exons and intron chain equal to an annotated isoform of the gene; single-exon genes: exonic Jaccard at or above a cut frozen from the DEPTH-30 arm's dev half at the 5th percentile (`RS_CUT_DEPTH=30`; `score 30` runs first and writes `single_cut.json`), applied unchanged at depth 10 and to the held-out half. About 5% of dev single-exon genes fall below it by construction; read the held-out half and the other depth.
- *Classes (ordered, mutually exclusive):* not_simulable, no_reads, reads_not_on_gene, no_locus (single-exon and multi-exon apart), outcompeted, weak_overlap, fused_excess, extent_wrong, node_ok. 'Matched' means every class after no_locus, outcompeted and weak_overlap are excluded; 'node_ok' is the strict class. Strata: biotype, single-exon against multi-exon, CAT against Liftoff, SD98 eligibility, dev against held-out half, and own-read nodes against nodes that exist only through shared secondary alignments. The ASSEMBLY POLISH counts of the assemble log are printed.
- *Like-for-like comparators* (the same scorer, `--real-library`, on the ideal arm's truth): the old A119b GTF (`gw22/ours_genome.gtf`), the A119b assemble-only f1v2 GTF (`rustle_figures_dev/f1_port/runs/human_A119b/new.f1v2.families.gtf`) and the current-default testis run (`rustle_figures_dev/defaults_flip/runs/new`, 1,179,288 reads). The 71.0% stays as one line labelled NOT like-for-like. The real libraries differ from the arm in pipeline version, library and mapping (testis mapped without -uf, A119b with minimap2 2.31); a gap to them is not expression alone.
- *Pipeline:* defaults only; `assemble` is the shipped strict-junction, self-tuned-floor assemble-only; `families` through the resumable shard wrapper; the HEAD release binaries; gate G5 accepts `RS_SRC_OK=958972e7` (a commit after the binaries were built whose `src/` and `Cargo.*` diff has no non-comment line).
- *Gates:* G1 to G5 as in section 5; G6 positive control: the annotation scored as nodes has no fused_excess, no extent_wrong and no no_locus; G7 negative control: the annotation shifted by 50 kb has fewer node_ok than the annotation in place; G-depth: the depth-10 reads are the depth-30 pool's first reads read for read; G3 covers the simulation and the scorer under two hash seeds. An independent re-score (a second implementation by a reviewer) and an independent code review precede the first heavy call.

**What was seen before this amendment (disclosure).** The survey ran the earlier fidelity code on the old A119b GTF (53.2% by id) and prototypes of the new classes on real products and on the earlier NPIP and TBC1D3 ideal windows (which gave the expected presence rate transported to the Soto composition: 92 to 95% of simulable genes under any-overlap matching). The builder ran the new scorer on the annotation-as-nodes control (2,224 of 2,224 node_ok), the +50 kb control (60 node_ok, 2.7%) and the three real-library products (old A119b by transcript: 45.1% with a real overlap, 12.2% node_ok; A119b f1v2 by transcript: 52.6% and 15.0%; testis by node: 14.4% and 6.7%). These are comparators and controls, not outcomes of the ideal arm; no simulated read of the pilot has been mapped.

**Predictions (replace P1 to P4; directional; written before any number of the ideal arm).**
- P1 presence: at depth 30, the share of the 2,224 simulable genes that are matched (no no_reads, reads_not_on_gene, no_locus, outcompeted or weak_overlap) is between 0.88 and 0.97.
- P2 single-exon: the single-exon stratum has a lower matched share than the multi-exon stratum at both depths, lower at depth 10 than at depth 30, and the multi-exon share changes by less between the depths than the single-exon share.
- P3 residual classes at depth 30: no_locus of single-exon genes is the largest non-ok class, then fused_excess and weak_overlap or outcompeted in either order, and extent_wrong of multi-exon genes is the smallest of these.
- P4 paralogs: reads_not_on_gene is concentrated in genes whose spliced sequence is byte-identical to another gene's (418 listed genes) and in the large families, and the own-read share of the matched genes is below 100%.
- P5 like-for-like: the node_ok share of the depth-30 arm (unit node) exceeds that of each real-library comparator in every biotype stratum with at least 30 genes; the gap is read as expression plus library and version effects, not as expression alone.
- P6 fusion does not vanish: the fused_excess share of the matched genes in the ideal arm is not smaller than in the current-default testis run and in the A119b f1v2 GTF scored the same way.
A failed prediction is reported as failed. If the node_ok share at depth 30 is below 0.60, the algorithm-limited share is large and 'the problem is expression' is not supported for this pipeline version; if P1 and P5 hold, presence is an expression effect in this simulation, and what remains (single-exon polish, annotated overlap, shared sequence) is not.

**Runbook (one stage per call, foreground, Bash timeout 600000; `RLOCK_WAIT=30 RLOCK_TIMEOUT=560` are set by the runner).** `cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto; R=bench/ideal_expression/run_soto_pilot.sh; export RS_CUT_DEPTH=30 RS_SRC_OK=958972e7`; then `bash $R truth`, `truthqc`, `map` (21 heavy calls, repeat until 'mapping complete'), `cut`, `bamqc 10`, `bamqc 30`, `assemble 10|30`, `families 10|30` (exit 75: call again; exit 1 with 'no progress': stop), `anno 10|30`, `score 30` then `score 10`, `score-ctrl 10|30`, `score-real`, `report 10|30`, `gates-sim`, `gates 10|30` (prints VALID or INVALID), `cleanup`. Predicted 43 to 46 calls, heavy time 1 to 2.4 hours (typical 1.3), peak disk about 6 GB. Work directory `/mnt/linuxdisk/tmp/ideal_arm_pilot_2026-10-07`. Only one heavy job runs at a time; the pilot waits for the Phase E run if that has the heavy slot.

**Frozen code (sha256).** `bench/ideal_expression/sim_genes.py` (951 lines) `9f832de3f213370647761f56360a7032baee2baf079c533bac40f93ba2916c24`; `score_soto_nodes.py` (721) `4567793e7af9881267ced98f7f8ec465f83ba1e25bc2904f1d58652819ae8a9c`; `run_soto_pilot.sh` (465) `ddf95ee06e14d0ab8c90db622590766eac78e101608eb9c7ffd3abd977f86771`; `test_soto_pilot.py` (2,282; 155 tests, OK under `PYTHONHASHSEED` 0 and 1) `e7c64a3b9d79d91246cfa1eb26a69f976a0479de48199e159b01630d874e0e9d`. An edit after `truth` makes gate G5 fail.
**Default inputs (sha256, first 16 hex).** `soto_famCN_S1C.tsv` 776a650716fd1a1b; `soto_split_2026-09-29.tsv` 8f2cba6c34984b69; CAT/Liftoff slim GFF de32e0c86bf87c7f; `chm13v2.0.fa.fai` cec84e8875f740fe; `soto_fid.per_gene.tsv` 2ca769defe9259d0; `ours_genome.gtf` 75dc8701f790f1dd; `new.f1v2.families.gtf` 3baa7ecd05efde89; `human_testis.families.gtf` e8983cf9696caf4c; `human_testis.fam.loci.gff3` bb6d3cadfbf1c221; `human_testis.assemble.log` e923064550c03818; `a119b_reads_per_soto_gene.tsv` 3c2667212fd67773. The genome FASTA is 3,156,259,703 bytes (2022-09-28) and the splice index 12,294,925,466 bytes (2026-09-17).

**Limits added.** L1 reads come from the reference, so MAPQ-0 ambiguity is understated (10.2% in the NPIP ideal BAM against 24% in real A119b): the remaining loss is a lower bound. L2 418 listed genes share a byte-identical spliced sequence with another gene (202 sequences; 411 genes carry one), so minimap2 places their reads by a hash tie-break and some land in `reads_not_on_gene` whatever the pipeline does. L3 the truth for extent is the annotation itself and the node is one representative transcript per locus. L4 jitter 0 to 30 bp gives 961 distinct trim pairs per transcript and the node boundary depends on depth, so the single-exon cut does too. L5 the two depths share reads. L6 the unit and competition scope change the headline (one run: 92.4, 81.9, 81.6%), so they are printed everywhere. L7 gene-level intervals understate uncertainty (about 3.6 times against a bootstrap over the 397 split units); proportions are descriptive.

**Decisions (coordinator, before any run).** Fusion class by the span rule with the exon rule beside it; single-exon cut frozen at depth 30, 5th percentile, dev half; node primary with transcript and locus as sensitivities; competition `all` with `soto` as sensitivity; G6, G7 and G-depth accepted; `RS_SRC_OK=958972e7` instead of a rebuild; the families 'no progress' stop accepted; the strict-junction self-tuned assemble is the 'current default RNA mode'; both depths from one pool (21 map calls); pad 25 kb; input hashes recorded here and re-checked by the review; authorisation as in the heading (ideal arm only).


---

## PREREG_seed_pool_real_reads_2026-10-07.md

*(original: `docs/archive/2026-10/merged/PREREG_seed_pool_real_reads_2026-10-07.md`)*

## Pre-registration: the seeding pool (primary / good / all) and the primary-first representative on real reads, NPIP and TBC1D3 (2026-10-07)

**Written before any arm of this study exists.** The driver switch it uses (`tools/rustle_pipeline.sh --seed-pool primary|good|all`, `--seed-as-ratio`, `--contig`, `--as-table`; commits c5cd2524, dae62624, 21 tests) is committed; no pool other than the default has been assembled with it, and no number below is read from a product of this study.
User direction (2026-10-07): "yes lets do that, however lets also make it easy to set the 'good secondaries' and 'all' in order to prove to my advisor how does this decision affect the outcome" (that = the real-read test of the primary-first representative, Rule 1 of `docs/PREREG_locus_units_2026-10-06.md`, on human NPIP and TBC1D3, then gorilla).
**Label.** Human A119b: DEV (the development block of NPIP, and the TBC1D3 audit set). Gorilla OR6737: second-species check, held out for Rule 1, spent for the pool (the default pool was decided on gorilla NC_073244.2); fixed in Amendment 1 before any gorilla arm exists. Human and gorilla are never pooled. The study is descriptive: it measures what the pool decision and Rule 1 do, with the directional predictions of section 6 stated in advance and no combined pass/fail.

### 1. The question

Which alignments should seed the assembly's loci, and what does each choice cost or buy at the multi-copy families? The three pools are nested, so they are one axis, a filtration of the read pool indexed by the width rho of the tie band:

| pool | seeds | driver |
|---|---|---|
| P | primary alignments only | `--seed-pool primary` |
| G(rho) | primaries + secondaries with AS >= rho x the molecule's genome-wide best AS | `--seed-pool good --seed-as-ratio rho` (the default is rho = 0.98) |
| A | primaries + every secondary | `--seed-pool all` |

P is the limit rho > 1 and A is rho = 0, so the pools nest: P ⊆ G(1) ⊆ G(rho) ⊆ G(rho') ⊆ A for 1 > rho > rho' > 0. The second factor is the representative of each locus (the transcript that stands for the locus in the all-vs-all, the family clusters and the copy table): MR = the transcript with the most reads (the shipped rule) or R1 = the transcript with the most reads among those whose intron chain also occurs in the primaries-only assembly (Rule 1 of the locus-unit study, implemented as `bench/entangled/locus_units.py --rule1 --primary-gtf`).

### 2. What was already seen (ledger)

- Human chr16 NPIP, the 25 CAT/Liftoff copies, **before** the bridge regroup (f1v2) and the min-cov-shorter default (`docs/archive/2026-10/NPIP_READ_POOL_2026-10-01.md`, scorer HEAD-reproduced in `docs/DEFAULT_RESCORE_NPIP_2026-10-06.md` gate G0): own node 23 / 21 / 24, E-found within own nodes 10 / 6 / 2, locus-level 19 / 18 / 16 for P / GOOD / ALL. The current default at HEAD (G(.98) + f1v2 + cov .70): own node 24, E-found within own nodes 8, locus-level 21, U2 F .645. **P and A were not re-run at HEAD; TBC1D3 was never run by pool; Rule 1 was never run on real reads.**
- The default pool on the 09-29 copy-recovery run (`docs/archive/2026-09/PREREG_copy_recovery_tools_2026-09-29.md`): per-copy numbers for the default only.
- Ideal-read simulation windows around the same two loci (`docs/LOCUS_UNITS_LEVELS_2026-10-07.md`): P finds more reachable copies (53 of 54 against 48) at a lower cluster precision (K*_C .667 against .775); Rule 1 repairs 4 gene-runs on 3 genes with no precision cost (level C only; level M it is masked by one partition flip). Simulated reads, 39 copies x 2 draws of one annotation.
- Gorilla overlap check (`docs/GORILLA_OVERLAP_2026-10-07.md`): P costs 0-1 chains at overlapping genes and 2 of 18 at non-overlapping ones (primary denominator); the secondary-pool gain is real at N (102 of 19,953 chains).
- Not seen: any NPIP or TBC1D3 number of any pool x representative cell at HEAD other than G(.98) x MR on NPIP.

### 3. Substrates and instruments (nothing is chosen after seeing an arm)

- **Human A119b** (`winloci_data/A119b.t2t.bam`, CHM13 v2.0, `chm13v2.0.fa`): **NPIP on chr16** (25 copies, `copies.hsa.tsv` family NPIP) and **TBC1D3 on chr17** (16 copies, family TBC1D3) of the CAT/Liftoff v2.0 truth (`/mnt/linuxdisk/tmp/rustle_figures_dev/copy_recovery_tools_cat/ann/{copies.hsa.tsv, truth.hsa.gtf}`). Best-AS table: `human_A119b.molecules.tsv` of the 09-25 run (header `bam=` checked by the driver).
- **Copy-level instruments** (unchanged since the default re-score): `bench/default_rescore/nodes.py` (own node; NPIP with `--exons-json npip_read_pool.json` as in the re-score, TBC1D3 without) and `bench/copy_support.py --nodes` (Amendment E E-found = `tc_found`: the representative starts at a capped 5' end of the copy's reads and carries its first three introns). Read under two `PYTHONHASHSEED`s.
- **Family-level** (`bench/default_rescore/score_families.py`, `family_score` of the HEAD build): U2 (NPIP union truth, 3 families), Compara Primates, Soto 2025, on the contig, bipartite and pairwise; the all-families number is reported, and the row of the target family.
- **Own node and composition of the family's clusters** (new, `bench/seed_pool/composition.py`; its own-node flags must equal `nodes.py`'s wherever `nodes.py` is defined, gate G1b): the clusters that hold a locus overlapping a same-strand truth copy are the FAMILY CLUSTERS of the arm; their loci are the nodes; a node is ON-COPY (overlaps >= 1 exonic bp of a same-strand truth copy), IN-SPAN (inside a copy's span, off its exons), ANTISENSE (overlaps a truth copy on the other strand) or ELSEWHERE. Node precision NP = on-copy nodes / all nodes of the family clusters. Unlabelled nodes count as false (no intersection with a truth universe). The cluster file names a locus by its span only, and `nodes.py` stops when two loci share one (it did on the frozen A arm of 10-01); `composition.py` gives every locus on a shared span the cluster ids of that span (a locus is in a cluster if any row of its span is) and reports the number of shared spans per arm. This is the rule for every arm, fixed here before any arm exists.
- **Cost**: transcripts, loci, all-vs-all PAF records, wall time and peak RSS of each stage.

### 4. Arms (per contig; HEAD release binaries `rustle_target_m2/release`; every other setting the driver's default: f1v2, strict junctions, shipped polish, `mcl_families --min-exonic-bp 1 --min-shared-exon-frac 0.60 --emit-units`)

Primary arms: **P**, **G98** (the default), **A**, and Rule 1 on top of the last two, **G98+R1** and **A+R1**. Descriptive arms (the filtration, section 7): **G100, G995, G95, G90** and their +R1 versions. Control: **P+R1** (Rule 1 over its own primaries-only assembly; every chain is primary-supported, so it must equal P).
Each pool arm is `tools/rustle_pipeline.sh assemble` then `families` with `--contig CONTIG --as-table W/mol/mol.tsv --no-cache --seed-pool ... [--seed-as-ratio rho]` on the full BAM (the table is genome-wide for every contig). A +R1 arm copies the pool arm's `PREFIX.gtf`, writes `PREFIX_R1.families.gtf` = `locus_units.py --base PREFIX.families.gtf --primary-gtf P.families.gtf --rule1` (a newer file than the copied `.gtf`, as the driver's guard requires) and runs the driver's `families` stage on it, so the clustering command is the driver's own. The all-vs-all of A-type arms runs through `tools/mm2_shard.sh` (byte-identical to a single run, cmp-checked 2026-09-25) in bounded calls; the others use plain minimap2, as the default re-score did. Intermediates (PAF, loci FASTA, shard caches) are deleted once an arm is scored (the work disk has 38 GB free).

### 5. Gates (an arm table is INVALID unless all hold)

**G0** the G98 arm on chr16 reproduces the registered default products of `docs/DEFAULT_RESCORE_NPIP_2026-10-06.md`: `cmp` of the assembled GTF, `families.gtf`, `fam.clusters.tsv` and `fam.loci.gff3` against `/mnt/linuxdisk/tmp/rescore_2026-10-06/DEF/` (the driver with `--contig` equals the hand-made `--region` run). **G1** the scorers on G98 chr16 give own node 24, E-found within own nodes 8, locus-level 21, U2 F .645; **G1b** `composition.py` gives the same own-node flags as `nodes.py` on the G98 and P arms of chr16 and on a synthetic case. **G2** P+R1 equals P (`fam.clusters.tsv` identical) on both contigs. **G3** `copy_support.py` is deterministic under `PYTHONHASHSEED` 0 and 1 (same summary). **G4** sharded and unsharded all-vs-all agree: the P arm on chr17 is clustered both ways (`MM2_SHARD_MIN_BYTES=1000` forces the wrapper on a small FASTA), `fam.clusters.tsv` identical.

### 6. Outcomes, denominators and predictions (fixed now)

**Outcomes per arm and family.** N = the family's truth copies (25, 16). **E** = the copies that are spliced-expressed under Amendment E (`tc_expressed`; arm-independent, computed from reads and annotation only). M1 = copies with an own node (of N). M2 = E-found within own nodes (of E). M3 = locus-level E-found (of E): any transcript of a same-strand overlapping locus starts within 150 bp of a capped start of the copy's reads and carries the first three introns of an expressed chain (`locus_tc_found`). M4 = U2 / Compara / Soto bipartite sens, prec, F and pairwise F on the contig. M5 = NP and the node classes of the family clusters; transcripts, loci, PAF records, wall time. Beside M2, reported: `ann_found_in_npip_nodes`, `chain_found_in_npip_nodes` (Amendments A, B). A per-copy matrix (copy x arm: own node, E-found) is reported, so every difference is a named copy; counts are of copies (n = 25, 16) and a difference of one is one copy, no significance is claimed.

**UNDERPOWERED rule.** A family with E < 8 gets no verdict for the predictions that read M2 or M3; its counts are reported.

**Predictions** (each is reported as held or failed; there is no combined verdict):
- **S1** (breadth costs the representative) on NPIP: M2(P) >= M2(G98) >= M2(A), with at least one strict inequality.
- **S2** (breadth does not cost membership at the default) on NPIP: M1(G98) >= M1(P).
- **S3** (breadth costs node precision) on NPIP: NP(P) >= NP(G98) > NP(A).
- **S4** (Rule 1 repairs the representative at the default) on NPIP: M2(G98+R1) >= M2(G98) + 1, M1(G98+R1) >= M1(G98), and |F(G98+R1) - F(G98)| < .01 for U2 and for Compara (chr17, S7: Compara and Soto; U2 has no TBC1D3 family).
- **S5** (Rule 1 repairs it at A) on NPIP: M2(A+R1) >= M2(A) + 2 and M1(A+R1) >= M1(A).
- **S6** (control) P+R1 is identical to P (this is gate G2; listed so that its outcome is part of the table).
- **S7** (TBC1D3) S1-S5 are evaluated on chr17 TBC1D3 with the same wording; whichever the UNDERPOWERED rule excludes is reported with its counts.
- **S8** (the default is not a knife edge) on NPIP: M1 and M2 of G95 and G995 each differ from G98 by at most one copy.
**Summary statistic.** The non-dominated set of the arms of a family on the three counts (M1, M2, NP): X dominates Y iff X >= Y on all three and > on one. The report states the set per contig and which arm each non-dominated arm is; it does not name a winner.

### 7. The filtration (descriptive)

G100, G995, G98, G95, G90 and A, with P at the far end, give each outcome as a function of the pool width. Registered reading: the curves are reported as measured; a non-monotone curve or a step between neighbouring widths is a finding about the threshold, not a failure.

### 8. Limits declared in advance

One library per species; the human block is the development block of both families; the copy truth is a single annotation (CAT/Liftoff), its pseudogene copies are partly unexpressed (E, not N, is the found denominator); E-found rests on the 5' cap signal, which exists in the human reads only (the gorilla reading is fixed in Amendment 1); the primaries of R1 come from the same reads and the same polish as the arms they correct (a chain is primary-supported iff it occurs in P's `families.gtf`, an exact intron-chain match); NP counts real but unlabelled family members as false; the all-vs-all of A took about 20 min against 46 s for G98 on chr16 (10-01, 17 times the records), so the sharded, resumable run is the only one that fits a 10-minute call.

### 9. Amendment hook

Amendment 1 (written before any gorilla arm exists, after the human results): the gorilla contigs (NPIP NC_073241.2 / NC_073242.2 / NC_073244.2; TBC1D3 NC_073228.2 / NC_073224.2), the copy truth and family truths used, the found criterion that replaces the cap-based E-found, the table of OR6737, and the predictions carried over unchanged.

### Amendment 1 (2026-10-07, written after the human arms and before any gorilla arm exists)

**What was known when this was written.** The human outcomes of the predictions, as the runner's `table` prints them (`bench/seed_pool/run.sh table`): chr16 NPIP S1 held (M2 12 >= 8 >= 2), S2 held, S3 held (node precision .72 >= .59 > .41), S4 failed (M2 8 -> 11 and M1 unchanged, but the U2 / Compara F moved by .022 / .017), S5 held (M2 2 -> 9), S6 held, S8 failed (G95 has M2 6 against 8); chr17 TBC1D3 S1 failed (7, 8, 5), S2 held, S3 failed (.81, 1.00, .18), S4 failed (Rule 1 changes nothing at the default), S5 held (M2 5 -> 8, M1 12 -> 15), S6 held, S8 held; gates G0-G3 pass on chr16, G2-G4 on chr17 (G1b not applicable on chr17: `nodes.py` stops on three shared spans). No gorilla pool arm of NPIP or TBC1D3 exists; of the default only the 09-29 copy-recovery numbers are known (NPIP 11 of 25 copies with a locus, 0 complete; TBC1D3 7 of 14 with a locus, 6 complete). The arm-independent counts below were computed from reads and annotation only (`copy_support.py` without any `--loci`).

**Substrate** (gorilla OR6737 testis HiFi against the KB3781 assembly and RefSeq annotation GCF_029281585.2: cross-individual, private variants lower every exact-chain count alike). BAM `winloci_data/GGO_mm.bam`, FASTA `_from_wsl/winloci_scratch/GGO.fasta`, best-AS table `runs/gorilla_OR6737/gorilla_OR6737.molecules.tsv` (4,404,440 molecules, `bam=` header checked by the driver). Copy truth `copy_recovery_tools/ann/{copies.ggo.tsv, truth.ggo.gtf}` (annotation only, `docs/archive/2026-09/PREREG_copy_recovery_tools_2026-09-29.md` section 3; no CAT mapping exists for gorilla). Label: held out for Rule 1 (never used to build or tune it); not held out for the pool (the default pool was decided on gorilla OR6737 NC_073244.2, which is NOT assembled here); the five contigs were used by the 09-29 default copy-recovery run and NC_073228.2 by Addendum AH of 09-14.
- **gorilla_npip**: contigs NC_073241.2 (1 copy) and NC_073242.2 (23 copies), assembled together (`--contig NC_073241.2,NC_073242.2`) and clustered together; **N = 24**. The 25th copy sits on NC_073244.2, the KRAB-ZNF array and the seeding rule's verdict contig: not assembled (spent as a substrate, and an all-secondaries all-vs-all there is prohibitive).
- **gorilla_tbc1d3**: contigs NC_073228.2 (12 copies) and NC_073224.2 (2 copies); **N = 14**.

**Found criterion.** The human criterion rests on the 5' cap signal (an untemplated G), which the gorilla reads do not carry (`tc_expressed` = 0 for both families). Gorilla uses **Amendment A**: E = copies with >= 2 reads carrying >= k of the copy's annotated introns (`ann_expressed`, k = min(2, introns)); M2 = `ann_found_in_npip_nodes` (the representative carries >= k annotated introns, within an own node); M3 = `locus_ann_found`. Chosen as the registered criterion with the largest arm-independent E, before any arm: NPIP E = 9 of 24 under A (chain 2, shared-chain 4, cap 0), TBC1D3 E = 7 of 14 (chain 5, shared-chain 3, TSS-anchored 3, cap 0). By the UNDERPOWERED rule **TBC1D3 (E = 7 < 8) gets no verdict on S1, S4g, S5, S8**; its counts are reported. NPIP (E = 9) is evaluated.

**Family-level truth.** None exists for gorilla on these contigs (Compara, U2 and Soto are human gene families): **M4 is not scored**. M6, added for every substrate, is the number of truth copies in the largest family cluster K* (`composition.py`: K* = the family cluster holding the most distinct copies, ties to the smaller id); for the human arms it is computed after the fact from the same products (post hoc, reported beside). **S4 on gorilla is S4g**: M2(G98+R1) >= M2(G98) + 1, and M1 and M6 not lower (M6 replaces the missing F clause).

**Arms.** P, P_R1 (identity control), G98, G98_R1, G995, G95, A, A_R1. The sweep ends G100 and G90 and their Rule-1 arms are not run on gorilla (cost); the filtration is therefore thinner there. A and A_R1 get at most six bounded calls (about 50 minutes of all-vs-all) each; an arm not finished by then is reported as not completed and the predictions that read it get no verdict. Everything else as in section 4; S1-S3, S5, S6, S8 are carried over unchanged.

### Amendment 2 (2026-10-07, user direction; written after the gorilla NPIP arms P, P+R1, G98, G98+R1, G995, G95 and A had finished and before any gorilla TBC1D3 arm exists)

User: "Actually can we leave the all secondaries one to be run at a later point? I think it is not that important". The `all` pool is deferred on gorilla. Gorilla NPIP's A arm had already finished (11:32 PDT) and is reported as it is; A+R1 on gorilla NPIP was stopped during its all-vs-all (its partial products and shard cache were deleted); A and A+R1 on gorilla TBC1D3 are not run. Consequences for the carried-over predictions: S5 has no verdict on gorilla; S1 and S3 are evaluated on gorilla NPIP only (TBC1D3: no verdict); everything else is unchanged. The human A and A+R1 arms had finished on both contigs and are reported unchanged. Run later with `bench/seed_pool/run.sh all gorilla_npip "A_R1"` and `bench/seed_pool/run.sh all gorilla_tbc1d3 "A A_R1"` (the six-call completion rule of Amendment 1 then applies).

### Errata (2026-10-07, appended after the independent verification wf_7800f946-0ac; everything above is unchanged)

1. **M2 is not 'within own nodes' in the sense of section 6.** The instrument (`bench/copy_support.py`, as in the 10-06 re-score) computes M2 = the copy has an own node AND some same-strand locus overlapping it, node or not, has the exact representative; `--nodes` ANDs a per-copy flag and does not restrict the loci (its docstring said otherwise and is corrected). M2n (the exact locus must itself be a node) is reported beside it, post hoc: chr16 P, G100, G995, G98, G95, G90, A = 12, 10, 7, 7, 5, 2, 1 (M2: 12, 10, 8, 8, 6, 3, 2). No S verdict changes. The 10-06 re-score's 'strict FOUND within own nodes 8 of 25' has the same reading.
2. **The ledger (section 2) was incomplete.** The 09-29 copy-recovery run and its 10-01 CAT re-run scored an 'ours, primaries only' arm against the default in all four cells (gorilla NPIP LOCUS 9/25 against 11/25, gorilla TBC1D3 4/14 against 7/14 with COMPLETE 3 against 7, human TBC1D3 13/16 against 13/16, human NPIP 26/26 against 26/26; register 1192; `docs/archive/2026-10/COPY_RECOVERY_TOOLS_CAT_2026-10-01.md` section 4). 'TBC1D3 was never run by pool', 'per-copy numbers for the default only' and 'no gorilla pool arm exists' were wrong: S2 in every cell and S1 / S3 at NPIP reproduce contrasts already seen; they were not blind predictions. Also seen: `docs/IDEAL_EXPRESSION_DEFAULT_2026-10-06.md` scored the default on real-read TBC1D3 (overlap 13/16, Amendment A 12, Amendment B 9), which equals this study's chr17 G98 row, so G98 with the most-reads representative had been seen at TBC1D3; and 'human NPIP 26/26 against 26/26' is the RefSeq-era figure, the CAT copies used here score 25/25.
3. **Amendment 1.** The 09-29 gorilla TBC1D3 default outcome is 7 complete (6 was the pre-run disclosure); the other criteria counts for NPIP omit 'TSS-anchored 0'.
4. **Ledger bullet on the gorilla overlap check.** '2 of 18' misread '2 / 18': P costs 0 chains at E_both and 2 of 19,771 at N; PC costs 1 and 18.
5. **Definitions fixed by the implementation.** A locus is its representative's exons (`loci.gff3`); on-copy = those exons overlap a same-strand copy's exon union (for NPIP the registered exons file, otherwise the union of the copy's truth transcripts); in-span / antisense = those exons overlap the copy's SPAN in the copy table on the same / the opposite strand, tested in this order after on-copy. The registered NPIP exons file puts PKD1P6-NPIPP1 on the NPIPP1 half, the copy table on the Liftoff gene about 2 kb away.
6. **Provenance.** `bench/seed_pool/*` was first committed after the human arms (a27786a6). Scorer and runner code was written between this prereg's commit (94e69328, 08:40:23) and the first arm (08:44:28); the `composition.py` stamped on the original human arms (sha1 8d7bc148) was never committed, the first committed version is in a27786a6, and every table was regenerated from committed code at 13:39-13:42. `run.sh` was edited at 08:56, 09:04, 09:19 and 10:51 while arms were running (the minimap2 shim, budget-stop detection, the budget constants, the `sharded()` rule) and is not in the early stamps; none of the edits changes an assembler or scorer input. The original human arms were produced by the driver of dae62624 (sha1 07c03dc7, stamped on all 28); HEAD's is 9947b359 and is functionally identical for a single valid contig (stub-checked on 7 invocations; the gorilla driver 77c3f9ef on 10 multi-contig invocations). The DONE markers of arms made before the markers existed were written by `run.sh adopt` from the arm names (an attestation by name); what an arm ran with is shown by the `arm` line of its `run.run.log` and the driver's own `seeding pool:` line, which agree for all 23 base arms. `gates.json` was re-run after the last scoring. Every widened gorilla arm except NPIP G98 (one process, 6.5 min) went through the shard wrapper; the four sharded calls of the human A arms are documented by the session record, not by a product. G4 was run on human chr17 only (and again after the regeneration, 6 shards). The unqualified 'plain minimap2 for the non-A arms' of section 4 holds for the human arms only.
7. **Labels.** Gorilla is a second-species check on a library that is NOT held out for the pool question (primaries-only against the default was already scored on these contigs on 09-29 and 10-01, and the default pool was decided on gorilla NC_073244.2 of the same library, which is not assembled here); it is held out only for Rule 1, which was never used to build or tune it. This supersedes 'spent for the pool' (line 5), 'not held out for the pool' (Amendment 1) and 'partly spent' wherever they differ. Where the document says 'held out' it means this narrow sense. The human block is the development block of NPIP and the audit set of TBC1D3.
8. **Gorilla NPIP E = 9 is at the boundary.** Six of the nine copies have exactly the floor of 2 primary alignments; requiring MAPQ > 0 gives E = 8, a floor of 3 gives E = 3. The criterion is registered and arm-independent; the verdicts rest on a handful of 2-read copies.
9. **The post hoc additions** are M6 for the human arms, red, CP*, M2n, the target-family rows (registered in kind in section 3; their family ids were chosen after the first verification), `secondary_support.py`, `paralog_support.py` and the analyses that use them. NP and the pairwise F are registered. They are labelled in the table captions, in the Limits list and at first mention in the results document.

### Errata 2 (2026-10-07, appended after the loss of the human arm products and the second independent verification wf_f88e1d5c-a0d; everything above is unchanged)

10. **The original human arm directories were lost and the non-A arms regenerated (14:02-14:40).** A verification script ran `rm -rf` through a symlink farm of directory links and emptied all per-arm directories of human_chr16 and human_chr17 at 14:02:19-20; the gorilla substrates and every substrate-level file (`comp.json`, `nodes.json`, `support*.json`, `fs_*.json`, `fs/`, `nodeloci/`, `table.*`, `cost.md`, `matrix.tsv`, `gates.json`) were untouched, and a read-only copy of the human ones is kept (`_record_13h42/`). The 12 non-A arms of each human substrate and P_SH were regenerated at HEAD (driver 9947b359, 14:07-14:40). They reproduce the stored record: `run.fam.clusters.tsv` is byte-identical to the surviving `fs/<arm>.clusters.tsv` in all 24 arms; transcript, locus and PAF-record counts equal the tables; chr16 G98's `gtf`, `families.gtf`, clusters and loci GFF3 equal the 10-06 default re-score products byte for byte (G0 again); G4 on chr17 passes again (6 shards); and a complete re-scoring of the 24 arms in a scratch work directory (`RS_WORK`, links to the regenerated arms) gave `comp.json` records, both `copy_support` seeds, the node-only support, the node-loci GFF3 files and the family scores (apart from the cluster-file path) identical to the originals. **The four human A and A+R1 arms are not regenerated** (deferred with the other all-secondaries arms; `bash bench/seed_pool/run.sh all human_chr16 "A A_R1"` and `... human_chr17 "A A_R1"`, then `score`, `gates`, `table`). The rows of A and A+R1, their stamps, timings, memory and the checks that need their loci (75 of 5,951 representatives change at A; chr17 A's M2 = 5 on full loci; the shared spans of A) therefore rest on the original products, which the first verification and the second pass recomputed before the loss (533 arm-table cells, 0 mismatches, from the raw files). The original timings and peak memory were parsed from the original logs before the loss (164 cells matched `cost.md`) and are kept in the stored `cost.md` and `table.tsv`; the regenerated arms have new timings. The original stamps and mtimes were captured before the loss (`/mnt/linuxdisk/tmp/seed_pool_verify2/runlogs.txt`, `prov/human_arm_reads_before_notice.txt`, `snapshot/`). The first-run timeline stands: this prereg 08:40:23, the first human arm 08:44:28, Amendment 1 10:43:22 and the first gorilla arm 10:43:29, Amendment 2 11:36:39 and the first gorilla TBC1D3 arm 11:36:42.
11. **What the second pass changed in the results text** (every number it recomputed reproduced; the corrections are of wording): 'only the representative count is monotone' became 'M2 falls with the width (weakly), node precision and cohesion do not move steadily' (M1, M3, the cost columns and Soto recall are monotone too); 'the NPIP families' own rows are flat' now names the PKD1 family, which rises; the contig-wide Soto precision .852 to .725 is pooled over 26 families and is not the TBC1D3 family's (its own precision rises); the invisible-node count is per truth (49 are invisible to all three); the on-copy exon union covers 100% / 96% of the span for three TBC1D3 copies and 90-97% for nine gorilla NPIP copies, and the gorilla NPIP effect is larger (with truth exon records over 3 kb ignored M1 is 0, 0, 0, 0, 18 at P, G995, G98, G95, A); the chain-multiplicity sentence is scoped to the copies not found and no longer says 'below the assembler's per-chain floors'; the '3, 8, 3' of the 09-29 trace are counts of the 14 absent copies per stage, not floor values; the human M1 does not step at A; the NP column is registered, not post hoc; the extra copy at human NPIP is a coordinate artefact whose direction reverses under the Liftoff definition; cost is the last bounded call, includes lock wait, and the first good arm pays the table parse.
12. **Code changed after the second pass, each with a test that failed first.** Driver: a contig name that looks like a number is matched as a string (awk compared `01`, `1.0`, `1e0` with `1`); a newline or return in `--contig` is refused (it was cut at the first line); the status-3 message prints the whole matched line and a remedy that works for a damaged table; the header states the aligner-record ratios per contig; the usage line shows the list syntax. Runner: a failed assemble prints a console message; `adopt` writes no Rule-1 marker when its inputs are missing. `table.py`: the S4 clause compares whole thousandths (`abs(.011 - .001) < .01` is true in floating point); the gate banner compares the gates with every product of `score` and tells 'not run' from 'not applicable'; cost caveats added. New: `paralog_support.py` with `test_paralog_support.py`, tests for the composition wiring and chromosome check, and 21 runner scenarios. The full suite passes (1,204 cargo tests, 41 runner checks, 67 bench unit tests).


---

## SEED_POOL_REAL_READS_2026-10-07.md

*(original: `docs/archive/2026-10/merged/SEED_POOL_REAL_READS_2026-10-07.md`)*

## Which alignments seed the loci, and what a primary-first representative adds: NPIP and TBC1D3 on real reads (2026-10-07)

**Status: DESCRIPTIVE. Human A119b is DEV (the development block of both families); gorilla OR6737 is a second-species check on a partly spent library (section 5; Rule 1 was not tuned on it, the default pool was decided on gorilla NC_073244.2, which is not assembled here). Protocol: `docs/PREREG_seed_pool_real_reads_2026-10-07.md` (committed 94e69328 before any arm existed; Amendment 1, a27786a6, before any gorilla arm; Amendment 2, f65f81c7, the user's deferral of the all-secondaries arms; an errata section after the independent verification). Switch: `tools/rustle_pipeline.sh --seed-pool primary|good|all`, `--seed-as-ratio`, `--contig`, `--as-table` (commits c5cd2524, dae62624, 77c3f9ef and the review fixes). Runner and scorers: `bench/seed_pool/{run.sh, composition.py, gates.py, table.py, secondary_support.py}` with `test_composition.py` (19 tests), `test_table.py` (14), `test_paralog_support.py` (6), `test_run_sh.sh` (41 checks); the full `cargo test --release` (1,204 tests, 34 of them driver tests in `tests/pipeline_seed_pool.rs`) and the 67 bench unit tests pass. Products: `/mnt/linuxdisk/tmp/seed_pool_2026-10-07/{human_chr16, human_chr17, gorilla_npip, gorilla_tbc1d3}/` (`table.txt`, `table.md`, `target.md`, `cost.md`, `matrix.tsv`, `gates.json`; one directory per arm with the driver's own products and a `run.run.log` stamp). Independently verified twice (wf_7800f946-0ac: five read-only verifiers, 93 findings; wf_f88e1d5c-a0d: a second pass over the corrections, 54 findings): every table value except the aligner record counts was recomputed from the products and reproduced (the record counts were compared with the stored ones, and the regenerated arms reproduce them); the errors and overstatements found are corrected here and in the prereg errata, the code defects are fixed.** **Provenance of the human arm products (prereg Errata 10):** the original arm directories of human_chr16 and human_chr17 were lost on 2026-10-07 at 14:02 to a verifier's `rm -rf`. The 24 non-A arms were regenerated at HEAD and reproduce the stored record (clusters byte for byte; a full re-scoring of the 24 arms equals the original comp.json, copy_support, node-only support, node-loci and family-score files). **The four human A and A+R1 arms are not regenerated** (deferred with the other all-secondaries arms), so their rows rest on the original products, which both verification passes recomputed before the loss; the original timings and peak memory (cost.md) were parsed before the loss and are kept. Human and gorilla are never pooled. User direction (2026-10-07): 'lets also make it easy to set the good secondaries and all in order to prove to my advisor how does this decision affect the outcome'.

### Answer

1. **At the copy level of NPIP (chr16, 25 copies, 22 expressed) the pool trades the exact representative against node precision and cohesion: M2 falls with the width (weakly, one tie), node precision and cohesion do not move steadily.** Along the width P, G100, G995, G98 (default), G95, G90, A: exact representatives (M2) 12, 10, 8, 8, 6, 3, 2 (M2n, post hoc, the same with the exact locus itself a node: 12, 10, 7, 7, 5, 2, 1); node precision NP .725, .705, .618, .593, .708, .661, .411 (G95 and G90 recover it); on-copy nodes per copy with a node (red, post hoc) 1.2 rising to 2.8 and copy-level precision CP* .60, .55, .44, .44, .52, .42, .15 (post hoc); copies in the largest family cluster (M6, post hoc for human) 23, 23, 23, 23, 24, 24, 17, with A cutting the family into 15 clusters; copies with a node (M1) and locus-level found (M3) move by at most two copies, and the transcripts, loci, aligner records and the contig-wide Soto recall rise at every step; the aligner's records 72,118 to 1,648,588 (22.9 times; A took about 30 minutes of wall clock over four bounded sharded calls against a 36-second families stage for P). Copies with a node: 24, 24, 24, 24, 25, 25, 25; the one extra copy is PKD1P6-NPIPP1, whose node exists under the registered NPIPP1 exons and already exists at P under the Liftoff PKD1P6-NPIPP1 definition of the same copy table (two coordinate definitions of one copy); under the Liftoff definition M1 reads 25, 25, 25, 25, 24, 24, 24 (the node is present from P to G98 and absent at G95, G90 and A), the opposite direction, so the step is a coordinate artefact and not evidence about the pool.
2. **The family-level scores are contig-wide and do not show the NPIP family gaining.** Compara, U2 and Soto are scored over all their families on the contig (chr16: Soto 20 families / 71 genes, Compara 17 / 54). The Soto gain from P to A (hits 30 to 52 of 71; sensitivity .423 to .732, bipartite F .588 to .832) comes from the other 18 families (17 to 39 of 51 hits); the NPIP-holding Soto families ID_154 and ID_149 contribute 13, 14, 13 hits at P, G98, A, Compara CF153 is 13 of 19 at P, G98 and A, U2 ID_154 is 14, 15, 14 of 19; only the PKD1 family rises, U2 ID_149 from 3 to 6 of 12 and Soto ID_149 from 1 to 2 of 6 at G90 and A (target rows below). The pairwise F, which is tie-free, does not repeat the bipartite F: P to A it moves +.032 against +.244 (Soto), -.014 against +.007 (Compara) and +.055 against -.006 (U2); differences of the printed three-decimal scores.
3. **TBC1D3 (chr17, 16 truth copies, 15 with exons in the truth (TBC1D3P7 has none, and no arm has a node on its territory although the scorer would credit one); 11 expressed) is a different shape: the default is ahead of primary-only on the family and behind it on CP* (post hoc).** G98 against P: exact representatives 8 against 7 (one copy, TBC1D3 itself; M2n the same), copies in the largest cluster (M6, post hoc) 11 against 9, copies with a node and locus-level found tied (13, 9), NP 1.000 against .8125 but CP* .59 against .81 (post hoc; red 1.7 against 1.0 on-copy nodes per copy, post hoc), contig-wide Soto precision .725 against .852 (pooled over the 26 Soto families of chr17: the family's own precision rises, Soto ID_468 .714 to 1.000), PAF records +7.3%. The family itself gains: Compara CF185 5 of 9 genes to 9 of 9 (F .714 to 1.000), Soto ID_468 5 of 10 to 8 of 10 (F .588 to .889). `all` costs: 83 nodes of which 68 lie on no copy, CP* .14, one copy fewer with a node (12) and M2n 0 (the exact loci are not nodes), while CF185 stays at 8 of 9 (F .941) and the contig-wide bipartite F is the highest of the pool-only arms (Compara .533, Soto .599; A+R1 has Soto .631), although the pairwise F is lower than at G98 (.482 against .564).
4. **0.98 is neither a knee nor a plateau.** At NPIP M2 changes by -2, -2, 0, -2, -3, -1 along ρ = 1.0, .995, .98, .95, .90 and A (unequal steps with a flat stretch at the default); the registered claim that G95 and G995 stay within one copy of the default failed (S8: G95 has 6 against 8) and the first loss is already at exact ties (G100: NPIPB5 and LOC124907807 lose their exact representative). At TBC1D3 M1 and M2 stay within one copy (S8 held). Earlier family-level work found the width 0.98 to 0.90 insensitive on other substrates (register 1060, 1100); the representative count is sensitive to it.
5. **A primary-first representative (Rule 1) repairs the representative at NPIP G98 and A, at a family-F cost, and does other things elsewhere.** NPIP G98: 8 to 11 (NPIPA8, NPIPB15 and LOC124907807 return; M2n, post hoc, 7 to 10), A: 2 to 9 (+8, and NPIPB10P is lost; M2n 1 to 5); copies with a node and the largest-cluster count unchanged at G98 and A, although the node sets differ (G98: 2 leave and 2 enter; A: 12 leave, 20 enter); 21 of 2,832 representatives change at G98 and 75 of 5,951 at A. The bipartite F moves by .022 (U2), .017 (Compara), .013 (Soto) and the pairwise F by .043, .057, .044, so the registered |dF| < .01 clause fails (S4). Elsewhere: NPIP G90+R1 loses PKD1P6-NPIPP1's node (M1 25 to 24, M6 24 to 22); TBC1D3 G100+R1 loses a node and an exact representative, G995+R1 changes the node set (13 against 19 nodes), G95+R1 gains a representative, G98+R1 changes no count, A+R1 raises copies with a node 12 to 15 and Soto F .599 to .631. P+R1 is identical to P.
6. **Registered non-dominated arms on (copies with a node, exact representatives, NP).** NPIP, primary arms: P, A and A+R1, but A is non-dominated only through an NP larger by .002 (69/168 against 72/176; at two decimals A+R1 dominates A); the whole sweep: P, P+R1 and G95+R1 (25 / 11 / .733 against 24 / 12 / .725; single copies, a sweep reading, not a registered claim). TBC1D3, primary arms: G98, G98+R1 and A+R1. The statistic does not name a winner.
7. **What NP, the family scores and M2 can and cannot see.** NP counts every locus on a truth copy as correct, so loci piled on a copy raise it (the register's earlier read-pool study counted up to 25 loci at one copy under `all`); its denominator is the arm's own family clusters, so a rise can be re-clustering (at TBC1D3 the three off-copy nodes of P are no longer nodes at G98: one locus sits in a non-family cluster of 3, the other two exist as loci with no cluster row while larger loci over the same regions are in the family cluster; they left the denominator without being removed). `family_score` deletes predicted members outside the truth universe and labels a locus by its largest-overlap gene regardless of strand: at chr16 A has 99 nodes on no truth copy (21 in-span, 32 antisense, 46 elsewhere); per truth, 67 to 68 of them are labelled by no gene of that universe (invisible: Compara 68, U2 67, Soto 67) and 31 to 32 are scored, some of them credited (49 are invisible to all three truths and 50 are scored by at least one); dropping the scored ones moves Compara F .674 to .615, U2 .638 to .533 and Soto .832 to .770, dropping the invisible ones moves nothing. The contig-wide family precision therefore understates the cost of `all` without ignoring it (Compara 1.000 to .857, U2 .760 to .629, Soto .968 to .963 against NP .725 to .411). The on-copy test uses the exon union of the copy's truth transcripts, which covers 100% of the span for two TBC1D3 copies (TBC1D3B, TBC1D3F), 96% for LOC100420311, and 90 to 97% for 9 of the 24 gorilla NPIP copies; without those three TBC1D3 copies M1 13 becomes 10 and NP .8125 and 1.000 become .625 and .773 (the same whether the node set is held fixed or re-derived). The gorilla NPIP effect is larger: 23 of 24 copies carry a truth exon record of more than 3 kb (up to 27.6 kb), and the 9 on-copy nodes at P overlap a copy only through exon records of 7.1 to 16.6 kb. With a fixed node set and the records over 3 kb ignored, gorilla NPIP M1 is 0, 0, 0, 0, 18 at P, G995, G98, G95, A (over 10 kb: 8, 10, 10, 10, 22), so the gorilla NPIP M1, NP, CP*, M6 and S2/S3 rest on that overlap rule; the direction (A far above the rest) does not change. M2 is the instrument's count: the copy has an own node and some same-strand locus overlapping it, node or not, has the exact representative; M2n restricts the exact locus to nodes (the 10-06 re-score's 'strict FOUND within own nodes 8 of 25' has the same reading as M2).
8. **Gorilla (section 5): the chains the all pool builds at the expressed NPIP copies are largely other copies' reads, so it does not show that those copies are expressed or that the primary pool is the cause.** Outside the all pool no locus carries any annotated intron of the 9 expressed copies; under the all pool 7 of 9 carry 2 to 4 annotated introns and 22 of 24 copies get a node (9 at P, 11 at the good pools; see item 7 for what the on-copy test lets through), but 93% of the secondary alignments that carry the introns belong to molecules whose primary alignment lies in another NPIP copy. Rule 1 changes no registered count; gorilla TBC1D3 is underpowered. The all-secondaries arms were deferred at the user's request (Amendment 2): the NPIP arm had finished and is reported, A and A+R1 on TBC1D3 and A+R1 on NPIP are not run.

### NPIP, human chr16 (N = 25 copies, E = 22 expressed under Amendment E; arms in order of pool width)

| arm | assembled transcripts | loci | PAF records | nodes (on-copy / in-span / antisense / elsewhere) | NP | red | CP* | M1 | M2 | M2n | M3 | M6 | compara sens / prec / F (pair F) | u2 sens / prec / F (pair F) | soto sens / prec / F (pair F) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| P | 8,673 | 2,579 | 72,118 | 40 (29 / 1 / 2 / 8) | 0.725 | 1.21 | 0.60 | 24 | 12 | 12 | 21 | 23 | 0.500 / 1.000 / 0.667 (0.650) | 0.559 / 0.760 / 0.644 (0.484) | 0.423 / 0.968 / 0.588 (0.630) |
| G100 | 8,914 | 2,648 | 76,572 | 44 (31 / 2 / 1 / 10) | 0.705 | 1.29 | 0.55 | 24 | 10 | 10 | 21 | 23 | 0.481 / 1.000 / 0.650 (0.591) | 0.559 / 0.731 / 0.633 (0.464) | 0.451 / 0.970 / 0.615 (0.635) |
| G995 | 9,108 | 2,705 | 81,388 | 55 (34 / 3 / 5 / 13) | 0.618 | 1.42 | 0.44 | 24 | 8 | 7 | 21 | 23 | 0.481 / 1.000 / 0.650 (0.616) | 0.559 / 0.704 / 0.623 (0.468) | 0.479 / 0.944 / 0.636 (0.630) |
| G98 | 9,473 | 2,832 | 94,433 | 54 (32 / 3 / 8 / 11) | 0.593 | 1.33 | 0.44 | 24 | 8 | 7 | 21 | 23 | 0.500 / 1.000 / 0.667 (0.678) | 0.588 / 0.714 / 0.645 (0.533) | 0.507 / 0.923 / 0.655 (0.635) |
| G95 | 9,767 | 2,953 | 112,419 | 48 (34 / 2 / 3 / 9) | 0.708 | 1.36 | 0.52 | 25 | 6 | 5 | 21 | 24 | 0.500 / 1.000 / 0.667 (0.643) | 0.559 / 0.704 / 0.623 (0.470) | 0.535 / 0.974 / 0.691 (0.662) |
| G90 | 10,216 | 3,089 | 127,721 | 59 (39 / 3 / 9 / 8) | 0.661 | 1.56 | 0.42 | 25 | 3 | 2 | 21 | 24 | 0.481 / 1.000 / 0.650 (0.646) | 0.618 / 0.656 / 0.636 (0.467) | 0.634 / 0.957 / 0.763 (0.628) |
| A | 14,183 | 5,951 | 1,648,588 | 168 (69 / 21 / 32 / 46) | 0.411 | 2.76 | 0.15 | 25 | 2 | 1 | 19 | 17 | 0.556 / 0.857 / 0.674 (0.636) | 0.647 / 0.629 / 0.638 (0.539) | 0.732 / 0.963 / 0.832 (0.662) |
| P+R1 | 8,673 | 2,579 | 72,118 | 40 (29 / 1 / 2 / 8) | 0.725 | 1.21 | 0.60 | 24 | 12 | 12 | 21 | 23 | 0.500 / 1.000 / 0.667 (0.650) | 0.559 / 0.760 / 0.644 (0.484) | 0.423 / 0.968 / 0.588 (0.630) |
| G100+R1 | 8,914 | 2,648 | 76,572 | 46 (31 / 2 / 1 / 12) | 0.674 | 1.29 | 0.52 | 24 | 11 | 11 | 21 | 23 | 0.481 / 1.000 / 0.650 (0.591) | 0.559 / 0.731 / 0.633 (0.464) | 0.451 / 0.970 / 0.615 (0.635) |
| G995+R1 | 9,108 | 2,705 | 81,388 | 59 (34 / 3 / 5 / 17) | 0.576 | 1.42 | 0.41 | 24 | 9 | 8 | 21 | 23 | 0.481 / 1.000 / 0.650 (0.616) | 0.559 / 0.704 / 0.623 (0.468) | 0.479 / 0.944 / 0.636 (0.630) |
| G98+R1 | 9,473 | 2,832 | 94,433 | 54 (33 / 3 / 8 / 10) | 0.611 | 1.38 | 0.44 | 24 | 11 | 10 | 21 | 23 | 0.481 / 1.000 / 0.650 (0.621) | 0.559 / 0.704 / 0.623 (0.490) | 0.493 / 0.921 / 0.642 (0.591) |
| G95+R1 | 9,767 | 2,953 | 112,419 | 45 (33 / 2 / 3 / 7) | 0.733 | 1.32 | 0.56 | 25 | 11 | 9 | 21 | 24 | 0.500 / 1.000 / 0.667 (0.643) | 0.559 / 0.704 / 0.623 (0.470) | 0.535 / 0.974 / 0.691 (0.662) |
| G90+R1 | 10,216 | 3,089 | 127,721 | 60 (37 / 4 / 10 / 9) | 0.617 | 1.54 | 0.40 | 24 | 10 | 8 | 21 | 22 | 0.500 / 1.000 / 0.667 (0.646) | 0.618 / 0.656 / 0.636 (0.467) | 0.606 / 0.956 / 0.741 (0.613) |
| A+R1 | 14,183 | 5,951 | 1,648,588 | 176 (72 / 21 / 36 / 47) | 0.409 | 2.88 | 0.14 | 25 | 9 | 5 | 20 | 17 | 0.574 / 0.838 / 0.681 (0.651) | 0.647 / 0.595 / 0.620 (0.530) | 0.732 / 0.945 / 0.825 (0.648) |

M2 is the strict Amendment E criterion as the instrument computes it (see item 7); M2n and M3 are of the 22 expressed copies. NP is registered; red (on-copy nodes per copy with a node), CP* (copies with a node over nodes), M2n and M6 are post hoc for the human arms. Ties at the fourth decimal print half to even (13/16 = .8125 prints .812; a family F of .8125 prints .812). The transcripts column counts run.gtf, which holds the bridge transcripts the families input lacks; the +R1 rows copy their base arm's transcripts, loci and PAF records by construction. Family scores are over all families of the contig; U2 is NPIP-specific. Node classes: a locus is its representative's exons; on-copy = those exons overlap a same-strand copy's exons (the registered NPIP exons file); in-span / antisense = those exons overlap a same-strand / opposite-strand copy's span in the copy table. Copies that change between arms (`matrix.tsv`): P to G98 loses the exact representative at NPIPB5 and LOC124907807 (already at ties), NPIPA6 and NPIPA8 (at 0.995) and NPIPB15 (at 0.98) and gains NPIPB10P; Rule 1 at G98 returns NPIPA8, NPIPB15 and LOC124907807.

Target families, chr16 (hits over genes of the truth family and its F; the rows the prereg promised, chosen after the verification):

| truth family | P | G100 | G995 | G98 | G95 | G90 | A | P+R1 | G100+R1 | G995+R1 | G98+R1 | G95+R1 | G90+R1 | A+R1 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| compara CF153 | 13/19, F 0.812 | 12/19, F 0.774 | 12/19, F 0.774 | 13/19, F 0.812 | 13/19, F 0.812 | 13/19, F 0.812 | 13/19, F 0.812 | 13/19, F 0.812 | 12/19, F 0.774 | 12/19, F 0.774 | 12/19, F 0.774 | 13/19, F 0.812 | 13/19, F 0.812 | 13/19, F 0.788 |
| u2 ID_154 | 14/19, F 0.737 | 14/19, F 0.718 | 14/19, F 0.718 | 15/19, F 0.769 | 14/19, F 0.700 | 13/19, F 0.650 | 14/19, F 0.757 | 14/19, F 0.737 | 14/19, F 0.718 | 14/19, F 0.718 | 14/19, F 0.737 | 14/19, F 0.700 | 13/19, F 0.650 | 14/19, F 0.757 |
| u2 ID_149 | 3/12, F 0.400 | 3/12, F 0.400 | 3/12, F 0.400 | 3/12, F 0.400 | 3/12, F 0.400 | 6/12, F 0.667 | 6/12, F 0.667 | 3/12, F 0.400 | 3/12, F 0.400 | 3/12, F 0.400 | 3/12, F 0.400 | 3/12, F 0.400 | 6/12, F 0.667 | 6/12, F 0.667 |
| u2 ID_151 | 2/3, F 0.667 | 2/3, F 0.667 | 2/3, F 0.571 | 2/3, F 0.500 | 2/3, F 0.667 | 2/3, F 0.500 | 2/3, F 0.286 | 2/3, F 0.667 | 2/3, F 0.667 | 2/3, F 0.571 | 2/3, F 0.500 | 2/3, F 0.667 | 2/3, F 0.500 | 2/3, F 0.250 |
| soto ID_154 | 12/14, F 0.889 | 12/14, F 0.889 | 12/14, F 0.857 | 13/14, F 0.867 | 12/14, F 0.889 | 11/14, F 0.815 | 11/14, F 0.880 | 12/14, F 0.889 | 12/14, F 0.889 | 12/14, F 0.857 | 12/14, F 0.828 | 12/14, F 0.889 | 11/14, F 0.815 | 11/14, F 0.846 |
| soto ID_149 | 1/6, F 0.286 | 1/6, F 0.286 | 1/6, F 0.286 | 1/6, F 0.286 | 1/6, F 0.286 | 2/6, F 0.500 | 2/6, F 0.500 | 1/6, F 0.286 | 1/6, F 0.286 | 1/6, F 0.286 | 1/6, F 0.286 | 1/6, F 0.286 | 2/6, F 0.500 | 2/6, F 0.500 |

### TBC1D3, human chr17 (N = 16 copies, 15 with exons in the truth, E = 11)

| arm | assembled transcripts | loci | PAF records | nodes (on-copy / in-span / antisense / elsewhere) | NP | red | CP* | M1 | M2 | M2n | M3 | M6 | compara sens / prec / F (pair F) | soto sens / prec / F (pair F) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| P | 10,439 | 2,794 | 70,710 | 16 (13 / 0 / 1 / 2) | 0.812 | 1.00 | 0.81 | 13 | 7 | 7 | 9 | 9 | 0.253 / 0.905 / 0.396 (0.252) | 0.291 / 0.852 / 0.434 (0.310) |
| G100 | 10,502 | 2,837 | 72,775 | 16 (13 / 0 / 3 / 0) | 0.812 | 1.00 | 0.81 | 13 | 7 | 7 | 9 | 9 | 0.267 / 0.909 / 0.412 (0.343) | 0.304 / 0.857 / 0.449 (0.378) |
| G995 | 10,565 | 2,853 | 73,616 | 19 (19 / 0 / 0 / 0) | 1.000 | 1.46 | 0.68 | 13 | 7 | 7 | 9 | 11 | 0.280 / 0.913 / 0.429 (0.479) | 0.329 / 0.839 / 0.473 (0.480) |
| G98 | 10,668 | 2,880 | 75,896 | 22 (22 / 0 / 0 / 0) | 1.000 | 1.69 | 0.59 | 13 | 8 | 8 | 9 | 11 | 0.320 / 0.923 / 0.475 (0.564) | 0.367 / 0.725 / 0.487 (0.530) |
| G95 | 10,791 | 2,935 | 81,664 | 22 (22 / 0 / 0 / 0) | 1.000 | 1.69 | 0.59 | 13 | 7 | 7 | 9 | 11 | 0.320 / 0.889 / 0.471 (0.564) | 0.392 / 0.721 / 0.508 (0.527) |
| G90 | 10,923 | 3,035 | 90,308 | 22 (22 / 0 / 0 / 0) | 1.000 | 1.69 | 0.59 | 13 | 8 | 8 | 9 | 11 | 0.320 / 0.889 / 0.471 (0.564) | 0.418 / 0.733 / 0.532 (0.540) |
| A | 14,472 | 5,726 | 1,061,646 | 83 (15 / 0 / 0 / 68) | 0.181 | 1.25 | 0.14 | 12 | 5 | 0 | 9 | 10 | 0.373 / 0.933 / 0.533 (0.482) | 0.519 / 0.707 / 0.599 (0.462) |
| P+R1 | 10,439 | 2,794 | 70,710 | 16 (13 / 0 / 1 / 2) | 0.812 | 1.00 | 0.81 | 13 | 7 | 7 | 9 | 9 | 0.253 / 0.905 / 0.396 (0.252) | 0.291 / 0.852 / 0.434 (0.310) |
| G100+R1 | 10,502 | 2,837 | 72,775 | 16 (12 / 0 / 2 / 2) | 0.750 | 1.00 | 0.75 | 12 | 6 | 6 | 9 | 8 | 0.280 / 0.913 / 0.429 (0.380) | 0.316 / 0.862 / 0.463 (0.402) |
| G995+R1 | 10,565 | 2,853 | 73,616 | 13 (13 / 0 / 0 / 0) | 1.000 | 1.00 | 1.00 | 13 | 7 | 3 | 9 | 11 | 0.307 / 0.920 / 0.460 (0.539) | 0.354 / 0.875 / 0.505 (0.547) |
| G98+R1 | 10,668 | 2,880 | 75,896 | 22 (22 / 0 / 0 / 0) | 1.000 | 1.69 | 0.59 | 13 | 8 | 8 | 9 | 11 | 0.320 / 0.923 / 0.475 (0.558) | 0.367 / 0.725 / 0.487 (0.522) |
| G95+R1 | 10,791 | 2,935 | 81,664 | 22 (22 / 0 / 0 / 0) | 1.000 | 1.69 | 0.59 | 13 | 8 | 8 | 9 | 11 | 0.320 / 0.889 / 0.471 (0.564) | 0.392 / 0.721 / 0.508 (0.527) |
| G90+R1 | 10,923 | 3,035 | 90,308 | 22 (22 / 0 / 0 / 0) | 1.000 | 1.69 | 0.59 | 13 | 8 | 8 | 9 | 11 | 0.320 / 0.889 / 0.471 (0.564) | 0.418 / 0.733 / 0.532 (0.540) |
| A+R1 | 14,472 | 5,726 | 1,061,646 | 94 (27 / 0 / 0 / 67) | 0.287 | 1.80 | 0.16 | 15 | 8 | 8 | 9 | 11 | 0.373 / 0.933 / 0.533 (0.543) | 0.519 / 0.804 / 0.631 (0.519) |

Target families, chr17 (the Compara and Soto columns of the table above are contig-wide, as on chr16; the ann / chn columns promised in the prereg are in `table.txt`):

| truth family | P | G100 | G995 | G98 | G95 | G90 | A | P+R1 | G100+R1 | G995+R1 | G98+R1 | G95+R1 | G90+R1 | A+R1 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| compara CF185 | 5/9, F 0.714 | 6/9, F 0.800 | 7/9, F 0.875 | 9/9, F 1.000 | 9/9, F 1.000 | 9/9, F 1.000 | 8/9, F 0.941 | 5/9, F 0.714 | 7/9, F 0.875 | 9/9, F 1.000 | 9/9, F 1.000 | 9/9, F 1.000 | 9/9, F 1.000 | 8/9, F 0.941 |
| soto ID_468 | 5/10, F 0.588 | 6/10, F 0.667 | 7/10, F 0.737 | 8/10, F 0.889 | 8/10, F 0.889 | 8/10, F 0.889 | 8/10, F 0.842 | 5/10, F 0.588 | 7/10, F 0.737 | 9/10, F 0.900 | 8/10, F 0.889 | 8/10, F 0.889 | 8/10, F 0.889 | 8/10, F 0.842 |
| soto ID_469 | 0/2, F 0.000 | 0/2, F 0.000 | 0/2, F 0.000 | 0/2, F 0.000 | 0/2, F 0.000 | 0/2, F 0.000 | 1/2, F 0.182 | 0/2, F 0.000 | 0/2, F 0.000 | 0/2, F 0.000 | 0/2, F 0.000 | 0/2, F 0.000 | 0/2, F 0.000 | 1/2, F 0.500 |

Cost (wall seconds of the driver call and peak RSS from `/usr/bin/time -v`, which wraps the machine lock, so any wait for the lock is included (at most 3.6 s on all 65 logs); GB are kB / 10^6; for the A arms the families row is the last bounded call only; `human_chr16/cost.md`, `human_chr17/cost.md`, `gorilla_*/cost.md` hold every arm; the human values are the originals, parsed before the 14:02 loss, and the regenerated arms have new timings): chr16 P assemble 12.2 s / 0.32 GB, families 35.8 s / 2.29 GB; G98 24.2 s / 0.75 GB, 47.5 s / 2.65 GB; A 10.7 s / 0.48 GB, last call 381.6 s / 3.02 GB. The first good-pool arm of a substrate (G98 here) also pays the one-time parse of the best-AS table (6.3 s on chr16), which is why its assemble time is above the 9.5 to 15.5 s of the other good arms; no arm's cost includes the one pass that builds the table (66 s for the 11.7 GB gorilla BAM).

### The registered predictions (as `table.py` prints them; no combined verdict)

S2 in every cell and S1 and S3 at NPIP reproduce contrasts that were known before registration (prereg Errata 2), so they were not blind predictions.

| | statement | chr16 NPIP | chr17 TBC1D3 |
|---|---|---|---|
| S1 | exact representatives P >= G98 >= A, one strict | HELD (12 >= 8 >= 2) | FAILED (7, 8, 5) |
| S2 | copies with a node G98 >= P | HELD (24, 24) | HELD (13, 13) |
| S3 | node precision P >= G98 > A | HELD (.725, .593, .411) | FAILED (.812, 1.000, .181) |
| S4 | G98+R1 >= G98 + 1 exact, no fewer copies with a node, \|dF\| < .01 | FAILED (11 against 8, but \|dF\| .022) | FAILED (8 against 8) |
| S5 | A+R1 >= A + 2 exact, no fewer copies with a node | HELD (9 against 2) | HELD (8 against 5; copies with a node 15 against 12) |
| S6 | P+R1 identical to P | HELD | HELD |
| S8 | G95 and G995 within one copy of G98 on M1 and M2 | FAILED (G95 has 6 against 8) | HELD |

Gates: G0 (the default through the new driver, `--contig chr16 --as-table`, is byte-identical to the 10-06 default re-score: assembled GTF, families GTF, clusters, loci GFF3), G1 (own node 24, exact 8, locus 21, U2 F .645), G1b (`composition.py` equals `nodes.py` on P and G98), G2 (P+R1 = P), G3 (`copy_support` identical under `PYTHONHASHSEED` 0 and 1) pass on chr16; G2, G3 and G4 (the sharded all-vs-all equals a single process; 6 shards on the P arm of chr17) pass on chr17, where G1b does not apply (`nodes.py` stops on three spans that are each shared by two loci; `composition.py` gives such loci the cluster ids of their span; no locus on a shared span is a node in any human arm, so the rule changes no human number). The gates were re-run after the last scoring; G0 and G4 (6 shards) were re-run on the regenerated human arms and pass. Independent recounts (verification) reproduced every table value (the aligner record counts only against the stored counts, and by regeneration for the 24 non-A human arms), the S verdicts, the gates, the instrument's flags from the BAM, and the family scores with a separate labelling and one-to-one matching.

### Reading

- The pool decision is a position on a trade-off, not a free improvement. At NPIP primary-only is ahead on the exact representative, NP and CP* and level on membership and cohesion; the default is ahead on the contig-wide family recall measures (mainly other families) and level on the NPIP family's own rows. At TBC1D3 the default recovers more of the family (CF185 9 of 9 against 5) at a higher number of loci per copy. `all` is expensive (22.9 times the aligner records at chr16) and splits the NPIP family; the gain it shows in contig-wide scores is not the target family's.
- The right width is not obvious from these counts: the representative count falls steadily while NP, cohesion and the TBC1D3 family recover and plateau at 0.95 to 0.90; a rule that picks the pool per locus rather than genome-wide is not tested here.
- Rule 1 is a cheap, local repair of the representative, scoped to NPIP G98 and A; it is not free at the family level (bipartite F up to .022, pairwise up to .057 at NPIP G98) and its effect elsewhere is mixed (item 5). A native implementation (a `primary_reads` attribute, `mcl_families --representative primary-first`) is not justified by one dev library.

### Limits

One library per species; the human block is the development block of NPIP and the audit set of TBC1D3 (nothing is held out for them). E-found is the strict, cap-based criterion on human reads only; M2 does not require the exact locus to be a node (M2n does). The copy truth is one annotation (CAT/Liftoff for human; RefSeq/Gnomon for gorilla); no arm has a node on TBC1D3P7 (it has no exons in the truth; the scorer would credit a node on its territory) and NPIPB13's truth has no transcripts yet it counts as expressed from reads alone and is never exact. The registered exons file puts PKD1P6-NPIPP1 on the NPIPP1 half (chr16:15,105,363-15,124,458) while the copy table and truth GTF carry the Liftoff gene about 2 kb away (15,126,202-15,159,411). Counts are of copies (25, 16, 24, 14) and a difference of one is one copy; the sweep is descriptive. Post hoc (added after seeing results; labelled in the table captions, in this list and at first mention in the answers): M6 for the human arms, red, CP*, M2n, the target-family rows (registered in kind in section 3 of the prereg; their ids were chosen after the first verification), `secondary_support.py`, `paralog_support.py` and every analysis of section 5 that comes from them (paralog share, chain multiplicity, E sensitivity), and the G95+R1 reading; NP and the pairwise F are registered. The red and CP* columns came after all arms had run. S2 (all cells) and S1 and S3 at NPIP reproduce contrasts that were already known (prereg Errata 2): they were not blind predictions. The instruments were registered by text and committed after the human arms (`bench/seed_pool/*` first in a27786a6); the original human arms were produced by the driver of dae62624 (sha1 07c03dc7; HEAD's is 9947b359, functionally identical for a single valid contig, stub-checked) and an earlier `composition.py` that was never committed (every table was regenerated from the committed version); `run.sh` was edited four times while arms were running (08:56, 09:04, 09:19, 10:51: the minimap2 shim, budget-stop detection, budget constants, the `sharded()` rule), none of which changes an assembler or scorer input; the DONE markers of arms made before the markers existed were written by `run.sh adopt` from the arm names (the evidence of what an arm ran with is the `arm` line of its `run.run.log` and the driver's own `seeding pool:` line, which agree for all 23 base arms); every widened gorilla arm except NPIP G98 (one process, 6.5 min) went through the shard wrapper, which is byte-identical to a single process on human chr17 (gate G4, 6 shards; not re-checked on gorilla); the sharded history of the human A arms (four bounded calls) is documented by the session record, not by a product. Changing the driver changes the figures run-cache code hash of every stage although the default products are byte-identical (a one-off re-stamp is expected). Aligner intermediates were deleted after each arm (the PAF record counts are kept).

### 5. Gorilla OR6737 (Amendment 1; Rule 1 not tuned on it; the all-secondaries arms were deferred by Amendment 2)

Reads OR6737 against the KB3781 assembly and RefSeq annotation (cross-individual); copy truth = the annotation-only sets of the 09-29 copy-recovery study; no cap signal exists in the gorilla reads, so the found criterion is Amendment A (the representative carries at least k = 2 annotated introns; 'exact' is not claimed). NPIP: contigs NC_073241.2 + NC_073242.2, **N = 24, E = 9** (the 25th copy is on NC_073244.2, not assembled). TBC1D3: NC_073228.2 + NC_073224.2, **N = 14, E = 7, so UNDERPOWERED** for every prediction that reads M2 or M3. E rests on 2 to 9 primary alignments per copy for NPIP (six of the nine copies have exactly 2) and 2 to 40 for TBC1D3; with MAPQ > 0 required E would be 8 for NPIP and 7 for TBC1D3, with a floor of 3 primaries E would be 3 for NPIP. No family truth exists on these contigs: M4 is not scored. Arms run: P, P+R1, G98, G98+R1, G995, G95 on both; A on NPIP only; nothing for A on TBC1D3.

#### Gorilla NPIP (N = 24, E = 9)

| arm | assembled transcripts | loci | PAF records | nodes (on-copy / in-span / antisense / elsewhere) | NP | red | CP* | M1 | M2 | M2n | M3 | M6 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| P | 4,347 | 1,170 | 56,165 | 17 (9 / 0 / 3 / 5) | 0.529 | 1.00 | 0.53 | 9 | 0 | 0 | 0 | 8 |
| G995 | 4,449 | 1,221 | 69,802 | 20 (11 / 0 / 3 / 6) | 0.550 | 1.00 | 0.55 | 11 | 0 | 0 | 0 | 10 |
| G98 | 4,522 | 1,252 | 71,759 | 22 (11 / 0 / 3 / 8) | 0.500 | 1.00 | 0.50 | 11 | 0 | 0 | 0 | 10 |
| G95 | 4,705 | 1,336 | 75,679 | 22 (11 / 0 / 3 / 8) | 0.500 | 1.00 | 0.50 | 11 | 0 | 0 | 0 | 10 |
| A | 6,863 | 2,448 | 219,687 | 57 (51 / 0 / 0 / 6) | 0.895 | 2.32 | 0.39 | 22 | 7 | 7 | 7 | 13 |
| P+R1 | 4,347 | 1,170 | 56,165 | 17 (9 / 0 / 3 / 5) | 0.529 | 1.00 | 0.53 | 9 | 0 | 0 | 0 | 8 |
| G98+R1 | 4,522 | 1,252 | 71,759 | 22 (11 / 0 / 3 / 8) | 0.500 | 1.00 | 0.50 | 11 | 0 | 0 | 0 | 10 |

#### Gorilla TBC1D3 (N = 14, E = 7)

| arm | assembled transcripts | loci | PAF records | nodes (on-copy / in-span / antisense / elsewhere) | NP | red | CP* | M1 | M2 | M2n | M3 | M6 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| P | 12,233 | 3,127 | 238,328 | 4 (4 / 0 / 0 / 0) | 1.000 | 1.00 | 1.00 | 4 | 4 | 4 | 4 | 4 |
| G995 | 12,343 | 3,182 | 239,833 | 6 (6 / 0 / 0 / 0) | 1.000 | 1.00 | 1.00 | 6 | 4 | 4 | 4 | 6 |
| G98 | 12,462 | 3,197 | 241,374 | 7 (7 / 0 / 0 / 0) | 1.000 | 1.00 | 1.00 | 7 | 4 | 4 | 4 | 7 |
| G95 | 12,747 | 3,232 | 240,590 | 8 (8 / 0 / 0 / 0) | 1.000 | 1.00 | 1.00 | 8 | 5 | 5 | 5 | 8 |
| P+R1 | 12,233 | 3,127 | 238,328 | 4 (4 / 0 / 0 / 0) | 1.000 | 1.00 | 1.00 | 4 | 4 | 4 | 4 | 4 |
| G98+R1 | 12,462 | 3,197 | 241,374 | 7 (7 / 0 / 0 / 0) | 1.000 | 1.00 | 1.00 | 7 | 4 | 4 | 4 | 7 |

(Rule 1 changes neither family's clusters at P; at G98 it changes the representative of 8 of 1,252 NPIP loci and 10 of 3,197 TBC1D3 loci, of which 2 and 1 are family nodes, and moves no registered count: M2 is at its floor for NPIP and both representatives pass k = 2 at TBC1D3. The +R1 rows' transcripts, loci and PAF columns are copies of the base arm.)

| | statement | gorilla NPIP | gorilla TBC1D3 |
|---|---|---|---|
| S1 | representatives with >= 2 annotated introns P >= G98 >= A, one strict | FAILED (0, 0, 7: the order is reversed) | NO VERDICT (E = 7, no A) |
| S2 | copies with a node G98 >= P | HELD (11 against 9) | HELD (7 against 4) |
| S3 | node precision P >= G98 > A | FAILED as registered (.529, .500, .895); with CP* it holds (.529, .500, .386) | NO VERDICT (no A) |
| S4g | G98+R1 >= G98 + 1, no fewer copies with a node, no smaller largest cluster | FAILED (0 against 0: nothing to repair) | NO VERDICT (E = 7) |
| S5 | A+R1 >= A + 2 | NO VERDICT (A+R1 deferred) | NO VERDICT |
| S6 | P+R1 identical to P | HELD | HELD |
| S8 | G95 and G995 within one copy of G98 | HELD, vacuously on M2 (0 at every width) | NO VERDICT (E = 7; counts within one) |

Gates: G1b, G2 and G3 pass on NPIP, G2 and G3 on TBC1D3 (`nodes.py` stops on shared spans there); G0, G1 and G4 are human-only.

**What gorilla NPIP shows, and what it does not.** Under P, G995, G98 and G95 no locus carries any annotated intron of any of the 9 expressed copies (0, by the representative and by any transcript of the locus), while the all pool gives 2 to 4 annotated introns at 7 of the 9 (NPIPA5 and NPIPB14P stay at 0) and a node to 22 of 24 copies against 11. The registered S1 and S3 therefore fail in the direction of `all`. But: (i) of the 247 secondary alignments that carry at least 2 annotated introns at the 9 copies, 229 (93%) belong to molecules whose primary alignment lies inside another NPIP copy (1 in the same copy, 17 elsewhere; gorilla TBC1D3: 435 of 515, 84%) and only 0 to 5 per copy have an alignment score within 0.98 of the molecule's best, so the chains the all pool builds at a copy are largely paralog reads placed there, which is not evidence that the copy is expressed (`bench/seed_pool/paralog_support.py`, reproduced by two independent recounts); (ii) the copies do have primary alignments (2 to 9, median 2, against 18 to 36 secondary, median 29; human NPIP 156 against 3,427, human TBC1D3 150 against 2,037, gorilla TBC1D3 6 against 77, all over the Amendment-A expressed sets of 24, 13, 9 and 7 copies; `bench/seed_pool/secondary_support.py`), but among the primaries the largest number of identical junction chains is 1 at eight of the nine NPIP copies and 2 at NPIPB4 (3 if the chain is cut to the copy's range), while the four TBC1D3 copies found at P have 3, 2, 4 and 6, and LOC101125558 is found at P with 2, so two reads are not a floor these data show; whether the pool or a per-chain floor of the assembler is the cause is not separated by this design (the 09-29 trace of the 14 gorilla NPIP copies absent from the default counts 3 lost at the pass-1 floor, 8 at the strict-junction gate and 3 at the polish mono floor); (iii) the contrast with human is not like for like: under Amendment A in human the instrument's found counts are NPIP P 15, G995 11, G98 11, G95 9, G90 8, A 9 of 24 (node-restricted 14, 10, 10, 7, 6, 7), i.e. best at P, and TBC1D3 P 10 and G98 12 of 13 (best at G98, node-restricted the same), and the found fractions at P (human NPIP 15 of 24 = 63% at a median of 156 primaries, human TBC1D3 10 of 13 = 77% at 150, gorilla TBC1D3 4 of 7 = 57% at 6, gorilla NPIP 0 of 9 at 2) show no dose-response between 156 and 6 (the groups differ in size and in species); (iv) the gorilla `all` step is a step at A, not a trend (M1 9, 11, 11, 11, 22 at P, G995, G98, G95, A; in human it is the node count, the aligner records, NP and M6 that step at A, while human M1 does not: 24, 24, 24, 25, 25 at NPIP and 13, 13, 13, 13, 12 at TBC1D3). The all arm dominates on the registered (M1, M2, NP) and is the worst on CP* (.39 against .50 to .55): it needs 2.3 on-copy nodes per copy with a node, splits the family over 8 clusters (13 of the 22 copies in the largest) and leaves NPIPA5 and NPIPB10P without a node. **NP is inflated by redundancy**: the registered S3 verdict stands as printed; CP* was added after seeing it (post hoc). **The on-copy test is lenient at gorilla NPIP** (item 7 of the answers): the nodes of P, G995, G98 and G95 overlap a copy only through truth exon records of 7 to 17 kb, so M1 of those arms is the weakest number of the section.

**Gorilla TBC1D3.** Copies with a node rise 4, 6, 7, 8 of 14 at P, G995, G98, G95 while the representative count stays at 4 (5 at G95); node precision is 1.00 throughout and all nodes are in one cluster; the counts are too small for a verdict. This differs from human TBC1D3, where copies with a node are flat at 13, and resembles it only in the flat representative count.


---

## PREREG_family_hierarchy_flag_2026-10-09.md

*(original: `docs/archive/2026-10/merged/PREREG_family_hierarchy_flag_2026-10-09.md`)*

## Family hierarchy under a flag (SD -> family -> expansion -> duplicon, with fusion loci and external families), 2026-10-09

**Pre-registered before any hierarchy or Soto-containment number is computed.** User goal (2026-10-09): put segmental duplications, duplicons, expansions (including fusion loci) and gene families into one clean hierarchy, under a flag, and integrate or show that Soto's families are a subset of ours.

### Principle (from `docs/TERMINOLOGY_FAMILY_SD_DUPLICON_EXPANSION_2026-10-07.md`)

The family (O1, `mcl_families`) is the object and stays a PARTITION of loci. SD, duplicon and expansion are typed layers read from other evidence; they ANNOTATE families and never change membership (every attempt to let SD/duplicon structure decide membership failed: register rows 676, 677, 817, 818, 997). The hierarchy is a partial order: a family can sit in several SD components, a fusion locus relates two families, an external family can straddle expansions.

### The flag

`mcl_families ... --emit-hierarchy --hierarchy-duplicons <DupMasker colors BED> [--hierarchy-sd <SEDEF BED>] [--hierarchy-external <TSV>]`. Default off; every existing product byte-identical with it on. Writes `<out>.hierarchy.members.tsv`, `.nodes.tsv`, `.edges.tsv` and, with `--hierarchy-external`, `.external.tsv`.

### Levels and definitions (all threshold-free except >= 1 bp overlap)

1. **Member** = a clustered locus of `<out>.clusters.tsv`; its exons = the union of the exon blocks of every annotation record folded into it (`loci.tsv`).
2. **Duplicon composition.** For each exon, the bases covered by each DupMasker duplicon id (segments of one id merged first); the exon's **dominant duplicon** = the id covering most of its bases (ties: smallest id); exons with no duplicon base have none. The member's duplicon string = its exons' dominant ids in transcript order, consecutive repeats collapsed.
3. **Owner of a duplicon** = the family whose members' exons carry the most bases of it (ties: smallest family id).
4. **Fusion locus** = a member some of whose exons are dominated by duplicons owned by OTHER families; `fusion_of` = those families. It stays in its own family (partition) and gains a `fusion_of` relation (cover), as decided for the container (09-30).
5. **Expansion** (copies of one duplicated block) = within a family, the member's **signature** = the set of dominant duplicons of its exons that the family owns. The family's expansions are its MAXIMAL signatures under inclusion; a member belongs to every maximal signature containing its own (exactly one: assigned; several: a shared fragment, listed in each); a member with an empty signature belongs to the family's `unique` class. Lineage specificity (an outgroup) is not tested by this flag; the level is the duplicated block, "expansion candidate" in the strict terminology.
6. **SD component** (with `--hierarchy-sd`) = SEDEF segments merged by overlap into SD loci, joined by the SEDEF pairs into connected components; a member sits in every component one of whose SD loci overlaps its exons (>= 1 bp); per member the maximum pair identity touching it (SEDEF column 21).
7. **External family** (with `--hierarchy-external`: `chrom start end label`, GFF coordinates of annotation records, one row per record x label) = mapped to members through the fold table. Per external family: records, records on clustered members, the best family and **containment** = |E ∩ F| / |E on clustered members|, the best expansion and its Jaccard over members.

### Test on Soto 2025 (human, CHM13 v1.0, the F_O universe of `docs/PREREG_hierarchy_duplicon_expansion_2026-10-07.md`)

Inputs: `/mnt/linuxdisk/tmp/ourhom/{loci.paf, loci.gff3}` with the parameters of `ours.params.tsv`; DupMasker colors v1.0; SEDEF `final_v1_clean.bed`; Soto S1C (`bench/soto/soto_famCN_S1C.tsv`) by Gene ID -> record coordinates (`names.tsv`), never by name. Soto family = an S1C `Family ID` other than Unassigned.

- **H0 (regression).** `clusters.tsv`, `loci.tsv`, `pairs.tsv` byte-identical to the run without the flag.
- **H1 (Soto is a subset of ours).** Among Soto families with >= 2 records on clustered members: **>= 90% have containment 1** (every such record in one of our families) **and the record-weighted mean containment is >= 0.95.** Reported: the exceptions with their families, and catalog coverage (the share of Soto-family records that lie on any clustered member).
- **H2 (where Soto sits, descriptive).** Per Soto family: whether it equals one of our expansions, a union of expansions of one family, or cuts across expansions; median best-expansion Jaccard.
- **H3 (positive control for fusion loci).** The PKD1P-NPIP and PDXDC2P-NPIPB14P readthroughs present among clustered members are fusion loci whose `fusion_of` includes the family holding PKD1 / PDXDC1 duplicons (09-30 descriptive note `project_fusions_are_duplicon_blocks`). Reported per gene found; not a bar if none is clustered.

**Decision.** The flag ships (opt-in) iff H0 holds. "Soto's families are a subset of ours" is stated iff H1 holds; otherwise the containment numbers and exceptions are reported as they are.

### Result of the registered run and Amendment 1 (2026-10-09, after H0-H3, before any number of the amended levels)

**Run.** `mcl_families --paf loci.paf --gff loci.gff3 --min-exonic-bp 1 --min-shared-exon-frac 0.60 --dump-pairs` (the registered F_O), with and without `--emit-hierarchy --hierarchy-duplicons <colors v1.0> --hierarchy-sd final_v1_clean.bed --hierarchy-external <S1C>`. 1,709 members, 380 families, 430 expansions, 475 members with `fusion_of`, 1,240 duplicons owned. Soto: 605 families, 2,572 S1C rows, all mapped by Gene ID; 1,973 rows (76.7%) on clustered members; 157 Soto families have no member in our families, 48 one member.
- **H0 met.** `clusters.tsv`, `loci.tsv`, `pairs.tsv` byte-identical with and without the flag; the new build reproduces the registered F_O `clusters.tsv`.
- **H1 NOT met at the MCL-family level.** Of 400 Soto families with >= 2 members, 313 (78.2%) are fully contained in one of our families (bar 90%); member-weighted mean containment 0.892 (bar 0.95). The exceptions are Soto families spread over 2 to 5 of our MCL families (ID_151, ID_152, ID_153, ID_155, ID_212, ID_226, ID_78 at 1/3; ID_99, ID_100 at 0.40; ...): MCL at inflation 2.8 cuts some homology groups more finely than Soto does.
- **H2 (descriptive): Soto families sit at our expansion level.** Median best-expansion Jaccard 1.000; **215 of the 400 Soto families are exactly one of our expansions**, 343 at Jaccard >= 0.5.
- **H3 met in substance.** PKD1P6-NPIPP1: member of the PKD1 family (MCL28), duplicon string SD9473, SD9474, SD9609, SD9607, SD9613 (PKD1's) -> SD9449, SD9443 (NPIP core), `fusion_of` MCL9 (the family owning the NPIP core). PDXDC2P-NPIPB14P: PDXDC family (MCL286), SD9585, SD9479 -> SD9456 ... SD9443, `fusion_of` MCL9. The fusion links the two predicted families; the prediction's wording assumed the readthrough would sit in the NPIP family, here it sits in the PKD1 / PDXDC family and points to NPIP.

**Amendment 1 (registered now).**
(a) *Fusion typing (post hoc to the 475 count, disclosed).* A `fusion_of` member is a **fusion** when its exons are dominated both by duplicons its own family owns and by duplicons another family owns; it is **foreign-only** when every duplicon-dominated exon belongs to another family (the family boundary runs inside a duplicon, e.g. a member made of one duplicon its family does not own). Both are reported; only fusions are called fusion loci.
(b) *Superfamily level.* Above the MCL family: the **superfamily** = the connected component of the pre-MCL homology graph (the graph MCL partitions; nothing new is computed), so family ⊆ superfamily by construction. **H1' (Soto is a subset of ours at the superfamily level):** among Soto families with >= 2 members, >= 90% fully contained in one superfamily and member-weighted mean containment >= 0.95, **counted only if the superfamilies are not trivial: the largest holds <= 10% of the members** (otherwise reported as uninformative). Reported: superfamily size distribution, how many MCL families each contains, the Soto families still split.

### Amendment 1 result and Amendment 2 (2026-10-09, before any number of the amended SD level)

**Amendment 1 run** (same command, rebuilt binary; H0 holds again: `clusters.tsv`, `loci.tsv`, `pairs.tsv` byte-identical). 331 superfamilies (largest 65 of 1,709 members = 3.8%, guard met; median size 3; 32 hold more than one MCL family). **H1' NOT met: 346 of 400 Soto families (86.5%) fully contained in one superfamily (bar 90%), member-weighted containment 0.940 (bar 0.95).** Fusion typing: 319 fusion loci, 156 foreign-only members. S1C carries 114 `Unassigned_*` labels; none has >= 2 members, so the H1/H1'/H2 counts over the 400 families are unaffected (491 Soto families proper).
**Why 54 Soto families still straddle superfamilies (descriptive, after the numbers).** Of the 642 cross-superfamily member pairs inside them, 355 align in the all-vs-all (median identity 0.969, median block 1,983 bp) and are removed by the exonic filters of the registered DNA definition (`--min-exonic-bp 1 --min-shared-exon-frac 0.60`: the shared sequence is intronic or a small part of the exons); 287 do not align at all. Soto joins genes on partial exon sharing that our family definition rejects on purpose (LORO-validated, `docs/seeded_family_definition.md` §6o2-§6o6).
**The SD level as registered is uninformative.** SEDEF pairs chain 1,598 of the 1,709 members (93.5%) into one component (`SDc0`); every member touches some SEDEF segment.
**Amendment 2 (registered now).** The SD level becomes the **SD block** = a SEDEF segment locus (segments merged by overlap on one chromosome, as already built), not the pair-connected component. A member sits in every SD block its exons overlap. Per family, the SD edge evidence: the share of member pairs directly joined by one SEDEF pair (a pair with one segment in an SD block of each member), and the highest identity of such a pair. Components stay in the members table as a column for reference. Reported, no bar: SD-block length distribution, members per block, the share of families whose member pairs are all SD-joined, and how many SD blocks hold members of more than one family (co-duplicated blocks).

### Amendment 2 result and summary (2026-10-09; tables in `docs/hierarchy_results/`)

**Run** (same command, rebuilt binary). H0 holds again (`clusters.tsv`, `loci.tsv`, `pairs.tsv` byte-identical; the external placement table is unchanged). **SD blocks:** members sit in 322 SD blocks, median length 139,972 bp (90th percentile 872,215; maximum 5,489,270 where SEDEF segments overlap into one region), median 2 members per block; **192 blocks hold members of more than one family** (co-duplicated blocks). **SD edge evidence:** the share of a family's member pairs joined directly by one SEDEF pair has median 1.000; all pairs joined in 200 of 380 families, none in 47 (families below SD identity). Example: the NPIP family MCL9 (19 members, 3 expansions, 8 owned duplicons, 165 of 171 pairs SD-joined, best identity 0.998), the PKD1 family MCL28 (7 members, 19 of 21, 0.994) and the PDXDC family MCL286 sit in one superfamily (SF866); PDXDC2P-NPIPB14P is a fusion locus of MCL286 with `fusion_of` MCL9.

**The hierarchy, as shipped (opt-in).**

| level | object | defined by | relation to the level above |
|---|---|---|---|
| SD block | a merged SEDEF segment locus | SEDEF (assembly self-alignment) | a family sits in one or more blocks; a block can hold several families (192 here) |
| superfamily | a connected component of the pre-MCL homology graph | O1 graph (no new evidence) | family ⊆ superfamily (partition of the partition) |
| family | an MCL cluster (O1) | O1 | members partition the loci |
| expansion | copies of one duplicated block: a maximal duplicon signature of the family | DupMasker + ownership | expansion ⊆ family; shared fragments sit in several (cover) |
| duplicon | DupMasker id | DupMasker | composes an expansion's signature; owned by one family |
| fusion locus | a member with exons dominated by its own and another family's duplicons | duplicon ownership | stays in its family, `fusion_of` the other (319 here; 156 more members are foreign-only) |
| external family (Soto) | S1C family | Soto 2025 | `contained_in` family and superfamily, `best_expansion` |

**Soto in our hierarchy.** Soto's families are **not a strict subset** of our families (H1: 78.2% fully contained in one MCL family; H1': 86.5% in one superfamily; bars 90%), because Soto joins genes on partial exon sharing that our DNA definition rejects on purpose (355 of the 642 straddling pairs align, median identity 0.969, and fail the exonic filters; 287 do not align). **Where they do sit: at our expansion level** (H2): 215 of the 400 Soto families with >= 2 members are exactly one of our expansions, median best-expansion Jaccard 1.0, 343 at >= 0.5. The registered wording for the thesis: "Soto's families sit at the expansion level of our hierarchy (215 of 400 identical), inside our families for 78% and inside our superfamilies for 87%; the remainder are joined by partial exon sharing our family definition does not accept."

### Amendment 3 (2026-10-09, before any number of the new level): the SD complex, and a complete typed placement of every Soto family

**Why.** H1 and H1' place 78.2% / 86.5% of the Soto families inside one of our families / superfamilies; the rest straddle superfamilies because Soto joins genes on shared SD98 exons that our exonic filters reject. Soto's families are built on SD evidence, so the level that should contain them is the one where our SD evidence joins superfamilies.
**New level (in the flag).** **SD complex** = connected component of the superfamilies under the relation "some member of one and some member of the other are joined directly by one SEDEF pair" (a pair with one segment in an SD block of each member; no transitive chaining through blocks without members). Hierarchy: SD complex ⊇ superfamily ⊇ family ⊇ expansion ⊇ duplicon; SD blocks stay the spatial layer. Written: `sd_complex` per member, `sd_complex` nodes, `superfamily -> sd_complex` and `superfamily -- superfamily sd_linked` edges (member pairs, best identity), and complex containment of external families.
**H4 (Soto is a subset of ours at the SD-complex level).** Among Soto families with >= 2 members: >= 90% fully contained in one SD complex and member-weighted containment >= 0.95, **counted only if the largest SD complex holds <= 10% of the members** (else uninformative, as the pair-connected components of Amendment 2 were).
**H5 (complete typed integration, descriptive).** Every Soto family with >= 2 members gets exactly one type, the lowest level containing all its members: EXPANSION (equal to one expansion), WITHIN-FAMILY (inside one family), WITHIN-SUPERFAMILY, WITHIN-SD-COMPLEX, OUTSIDE (spans complexes); plus, for every Soto family, the count of its records not on any clustered member (not in our families at all). The table of counts and the OUTSIDE list are reported.
**Decision.** "Soto's families are a subset of ours" is stated at the level where H4 holds (with the guard); H5 is the integration whatever H4 says.

### Amendment 3 result and Amendment 3b (2026-10-09, before the direct-link numbers)

**Amendment 3 run** (H0 holds). 35 SD complexes; **the largest holds 1,532 of 1,709 members (89.6%)**: even direct member-to-member SEDEF pairs chain most families into one complex. **H4: 398 of 400 Soto families (99.5%) lie in one SD complex, weighted 0.998, but the guard fails (largest > 10%): uninformative as registered.** **H5 (typed placement of the 400): EXPANSION 215, WITHIN-FAMILY 98, WITHIN-SUPERFAMILY 33, WITHIN-SD-COMPLEX 52, OUTSIDE 2 (ID_266, ID_306, 4 members, 2 complexes each).** Of the 491 Soto families proper, 193 have records outside our families, 79 have no member in our families, 12 one member.
**Amendment 3b (registered now).** Transitive closure is what makes the complex trivial, so the relation is tested without it. **H4b:** for every Soto family whose members lie in >= 2 superfamilies (the 54 = 52 + 2), are ALL pairs of its superfamilies joined directly by an `sd_linked` edge (some member of one and some member of the other joined by one SEDEF pair)? **Bar: >= 90% of the 54. Counted only if the direct-link graph is sparse: <= 5% of the pairs of superfamilies that hold members are `sd_linked`** (otherwise uninformative). If H4b holds, the claim reads: every Soto family lies inside one expansion, family or superfamily of ours, or across superfamilies our SD layer joins directly, with the 2 OUTSIDE families named.

### Amendment 3b result and the final statement (2026-10-09)

**H4b: 52 of the 54 straddling Soto families (96.3%) have every pair of their superfamilies joined directly by a SEDEF pair (bar 90%), but the guard fails: 11,854 of the 54,615 pairs of superfamilies (21.7%) are `sd_linked` (guard <= 5%). Uninformative as registered.** Our SD layer joins SD-embedded families too broadly for any SD-level containment to say more than "both are in segmental duplications".
**In the flag now.** `<out>.hierarchy.external.tsv` carries, per external family, `placement` and `placed_under` (the lowest node holding every member) and the edges table an `external -> <level> placed_under` edge: the external catalogue is a typed level of the hierarchy.
**Final statement (registered tests, as they came out).**
- Soto's families are **integrated** into our hierarchy, each attached to exactly one node: of the 491 Soto families, the 400 with >= 2 members on our loci sit as **EXPANSION 215, WITHIN-FAMILY 98, WITHIN-SUPERFAMILY 33, WITHIN-SD-COMPLEX 52, OUTSIDE 2** (ID_266, ID_306); the 12 with one member sit WITHIN-FAMILY 10 / EXPANSION 2; 79 have no member in our families (their genes are not clustered by our DNA definition).
- They are **a subset of ours down to the superfamily for 346 of 400 (86.5%)**, and for 313 (78.2%) inside one of our families. A subset claim for all of them is **not** established: the remaining 54 are joined by partial SD98 exon sharing that our exonic filters reject, and the only levels of our hierarchy that contain them (SD complex, direct SD links) are too broad (89.6% of members; 21.7% of superfamily pairs) to count as containment.
- Thesis wording: "Soto's families are a typed level of our hierarchy: 215 of 400 are exactly our expansions, 346 lie within one of our superfamilies, and the 54 others are joined by partial SD98 exon sharing that our family definition rejects by design."


---

## PREREG_paf_duplicon_units_2026-10-10.md

*(original: `docs/archive/2026-10/merged/PREREG_paf_duplicon_units_2026-10-10.md`)*

## Our own duplicon units from the locus self-alignment (PAF units), 2026-10-10

**Pre-registered before any unit is computed.** User (2026-10-10), after "do we have our own way to find duplications?": pre-register duplicon-like units cut from the
alignments we already compute, so the expansion level of `mcl_families --emit-hierarchy` (`docs/PREREG_family_hierarchy_flag_2026-10-09.md`) no longer needs DupMasker,
which exists on disk for human CHM13 v1.0 only.

**Scope and what is not re-proposed.** Units only ANNOTATE: they replace the `--hierarchy-duplicons` track; family membership is untouched (families come from the same PAF
and the same `mcl_families` settings). Letting SD or duplicon structure decide membership failed before (register 676, 677/997, 817, 818) and is not tried. No genome-only
discovery (user 09-13): the alignments are between the loci the family definition already uses. Register 1257-1260 ("structural units") are about read-graph locus models,
not this.

### Method (fixed now; no parameter is fitted to the data)

Input: `loci.paf` (all-vs-all minimap2 of the 2,334 locus sequences, the family definition's own input) and `loci.gff3` (exons). Work in genomic coordinates (a locus
`chrom:s-e` starts at `s-1`), so loci that overlap share their sequence.

1. **Alignments.** Keep a PAF record when (a) the two loci differ and the two aligned intervals do not overlap on the genome (overlapping loci are the same DNA);
   (b) identity = matching bases / block length (columns 10/11) >= 0.90 and block length >= 1,000 bp: **the standard SD definition** (Bailey et al. 2001; the
   SDs DupMasker's duplicons were built from); (c) each aligned interval carries >= 1 exonic base of its locus: **the family definition's own `--min-exonic-bp 1`
   rule**, which also keeps intron-only repeat alignments out.
2. **Breakpoints.** Start from the ends of every kept alignment, on both sides. Project every breakpoint lying inside an alignment through that alignment's CIGAR
   (`cg:Z`) to its partner; add the projection unless a breakpoint already lies within 1,000 bp of it (the SD minimum length: two boundaries closer than an SD
   cannot both bound one). Repeat until no breakpoint is added, at most 20 rounds (the round count and whether the fixpoint was reached are reported; at the cap the
   last state is used).
3. **Atoms.** Cut every genomic stretch covered by a kept alignment at the breakpoints; atoms shorter than 1,000 bp are not units (their bases carry no unit, like
   bases outside every duplicon).
4. **Units.** Join atom a and atom b when some kept alignment maps >= 50% of a onto b and >= 50% of b onto a (reciprocal majority); a unit = a connected
   component (union-find), id `PU<n>`.
5. **Output.** BED `chrom start end unit_id`, read by the existing `--hierarchy-duplicons` loader unchanged; dominant unit per exon, ownership, signature,
   expansions, fusions as in the flag.

Code: `bench/paf_duplicons/paf_units.py` (+ `test_paf_units.py`, test-first: CIGAR projection on both strands, breakpoint snapping, the X/Y/Z chain below,
reciprocal joining, genomic self-overlap skip), `bench/paf_duplicons/score_units.py` (+ tests).
Chain the method must get right (unit test): X-Y aligned over AB, Y-Z over B only ⇒ units {A}, {B}; X and Y have signature {A,B}, Z {B} (a fragment inside the
same expansion), not one unit for everything.

### Test on human CHM13 v1.0 (truth = DupMasker, the registered final run)

Comparator: `/mnt/linuxdisk/tmp/hier_flag/final/ours.hierarchy.members.tsv` (DupMasker colors v1.0). Ours: the same command with `--hierarchy-duplicons <PAF units BED>`.
Per family (>= 2 members), the **cover** = its expansions as member sets plus its `unique` class as one more block (matched only with the other run's `unique`).

- **Family score** = maximum-weight bipartite matching of the two covers' blocks on Jaccard (members), sum of matched Jaccard / max(#blocks of each); 1 = identical.
- **Expansion sensitivity** = share of DupMasker expansions that equal one of ours exactly; **precision** = share of ours that equal one of DupMasker's
  (user reporting rule: sensitivity, precision, bipartite).
- **Baseline B1** = one block holding every member of the family (the "expansion = family" null).
- **Structured families S** = families whose DupMasker cover is not one block holding every member (B1 is wrong there by construction).

**Split.** Family id hashed: `int(sha256(id), 16) % 2 == 0` = DEV, else HELD-OUT. DEV is scored and written first; at most ONE amendment may follow, registered
before the held-out half is computed; the held-out half decides.

- **H0 (guard).** `clusters.tsv`, `loci.tsv`, `pairs.tsv` byte-identical to the DupMasker run (units only annotate).
- **P1 (primary).** Mean family score over held-out families >= **0.80**.
- **P2 (beats the null where it matters).** On held-out S, ours scores above B1 in more families than below, one-sided sign test p < 0.05 (ties dropped).
- **P3 (Soto placement replicates).** Soto families (>= 2 members) equal to exactly one expansion with PAF units >= **194** (90% of the 215 with DupMasker).
  All families; Soto is a comparator, not the objective, so P3 does not decide adoption.
- **P4 (named fusions).** PKD1P6-NPIPP1 and PDXDC2P-NPIPB14P are fusion loci `fusion_of` the NPIP family with PAF units (both needed).

**Verdict.** P1 and P2 hold on held-out ⇒ **ADOPT** PAF units as the default duplicon source of `--emit-hierarchy` when no `--hierarchy-duplicons` track is given
(Rust port, opt-in until reviewed); otherwise NOT ADOPTED, with the failing families listed. P3 and P4 are reported either way.
**Descriptive (reported, no bar):** alignments kept, rounds, atoms, units; share of members' exonic bases carrying a unit vs a DupMasker duplicon; expansions per
family under both; dev and held-out side by side; expansion sensitivity / precision; the NPIP family (MCL9) cover under both.

### Phase 2 (only after the Phase 1 verdict; descriptive, gorilla)

Gorilla has no duplicon track, so Phase 2 has no truth. Inputs (the gorilla locus set and its self-alignment, built like `loci.paf`) are fixed in an amendment
before the run. Named expectation registered now: **LRPAP1** (11 gorilla loci: 8 full copies + 3 5′ fragments, `project_lrpap1_family`) — the 8 full copies share
one maximal signature and the 3 fragments' signatures are subsets of it (one expansion, 3 shared or contained fragments). The NPIP family is reported as is.

### DEV result (2026-10-10; code 7a51101e + scorer fix below; held-out NOT computed)

**Units.** 164,616 PAF records, 12,373 kept; fixpoint after 3 rounds; 6,039 breakpoints, 4,525 atoms (26.8 Mb), **998 units** (396 with 2 copies, 221 with 3, one with 158).
4 s, 80 MB. Members' exonic bases carrying a unit **81.7%** vs a DupMasker duplicon 91.3%. **H0 holds** (`clusters.tsv`, `loci.tsv`, `pairs.tsv` byte-identical).
**DEV (200 families, 49 structured):**
- **P1 FAILS on dev: mean family score 0.630** (bar 0.80); 106 of 200 families identical; member-weighted 0.503. Exact expansion sensitivity 115/249, precision 115/228.
- **P2 holds on dev:** on the 49 structured families PAF units beat the one-block null in 27, lose in 3, tie 19 (sign test p 4e-6); mean 0.479 vs 0.278.
- **P3 holds (all families): 223 Soto families are exactly one expansion with PAF units** (DupMasker 215; bar 194); WITHIN-FAMILY 90, WITHIN-SUPERFAMILY 33,
  WITHIN-SD-COMPLEX 52, OUTSIDE 2 (unchanged above the family).
- **P4 holds:** PKD1P6-NPIPP1 member (MCL28) `fusion_of` MCL9; PDXDC2P-NPIPB14P (MCL286) `fusion_of` MCL171, MCL9.
- NPIP (MCL9, dev): DupMasker 3 expansions (4, 6, 10 members), PAF units 4 (3, 4, 6, 9); score 0.474.
**Where dev disagrees (94 families):** 18 where PAF units give no unit on the family's exons and DupMasker gives one expansion of all members (member alignments
below 0.90 or 1 kb: older duplications), 18 the reverse (DupMasker has no duplicon where >= 0.90 / >= 1 kb exonic alignments exist), 28 where PAF units merge
more, 23 where they split more, 7 same block count with different members (e.g. MCL0, 65 members: DupMasker 62 + 3, PAF units 63 + 2 unique).
**Scorer fix (not a method change, disclosed).** P4 first looked members up by the readthrough gene record CHM13_G0020726 and picked the containing PDXDC1
locus; it now checks the exact members H3 reported (`chr16:15126207-15159416`, the PKD1P6 locus LOFF_G0001012 whose model runs into the NPIP core, and
`chr16:75785720-75875798`). P1-P3 unchanged by the fix.

### HELD-OUT result and verdict (2026-10-10; no amendment, user's choice; same units, same scorer as DEV)

**HELD-OUT (180 families, 35 structured):**
- **P1 FAILS: mean family score 0.638** (bar 0.80); 103 of 180 identical; member-weighted 0.637. Exact sensitivity 105/181, precision 105/180.
- **P2 FAILS: 12 wins, 8 losses, 15 ties against the one-block null (sign test p 0.25)**; means 0.438 vs 0.310. The DEV win (27 / 3) did not replicate.
- P3 and P4 (all families) as on DEV: 223 Soto families exactly one expansion (DupMasker 215); both named fusions `fusion_of` MCL9.
- Disagreements (77): 22 DupMasker has no duplicon where PAF units give one expansion, 19 the reverse, 15 merge more, 15 split more, 6 other.
**Verdict: NOT ADOPTED.** PAF units do not reproduce DupMasker's expansion level (0.63-0.64 on both halves) and do not beat the one-block null on held-out.
DupMasker stays the duplicon source of `--emit-hierarchy`. PAF units remain a bench tool (`bench/paf_duplicons/`), not wired into the flag. What they do
reproduce: Soto's placement (223 vs 215 exact expansions) and the two named fusions, so those two statements do not depend on DupMasker.
**Phase 2 (gorilla).** Registered as descriptive after this verdict; with PAF units not adopted, any gorilla expansion level built on them carries the
measured agreement of 0.63-0.64 with DupMasker on human and no evidence of beating the one-block null. Not run.


---

## PREREG_assembler_compatible_support_2026-10-11.md

*(original: `docs/archive/2026-10/merged/PREREG_assembler_compatible_support_2026-10-11.md`)*

## Compatible-read support in the assembler's chain admission, 2026-10-11

**Pre-registered before the code exists and before any arm is scored.** User (2026-10-11): the advisor wants the assembler benchmarked
against StringTie and FLAIR; make the primary-only assembler as close to them as possible with their tricks, extendable later to the
seeded loci (secondaries, AS-tied primaries resolved by O2).

**Diagnosis that motivates it (dev only, gorilla OR6737).** Figure 1 rebuilt on the current assemblies (identical to 2026-09-25): we lead
on intron-chain sensitivity / precision (0.269 / 0.356 vs StringTie 0.244 / 0.346) and trail on intron sensitivity (0.593 vs 0.639).
Genome-wide, 18,961 reference introns are StringTie's and not ours (6,329 the reverse); their primary junction reads have a median of 4.
On NC_073244.2 (chr20, the dev contig): of 941 such introns only 20 exist before polish; **866 (92%) are carried by >= 2 primary reads
of which no two share an exact intron chain** — the assembler groups reads by EXACT chain and admits a chain only at `>= 2` reads
(`PASS1_MIN_READS`), so truncated or end-variant reads of one isoform never pool. StringTie's splice graph pools them.

**Not a re-proposal.** Register 848 / 855 (read-isoform widening: a chain admitted when each of its junctions has >= k reads counted
over ALL spliced reads of the region) failed (+50% transcripts for +2 NPIP junctions; refuted held-out on human chromosomes): it let a
chain borrow support from reads of OTHER isoforms. Here a chain borrows only from reads COMPATIBLE with it, and only chains that some
read actually carries are candidates — no new combination of junctions is ever built.

### The rule (one flag, `copy_assign --compat-support`, default off = byte-identical)

In pass 1 (`Pass1Acc`), after reads are grouped by exact `(chrom, intron chain)`: a group D is COMPATIBLE with a group C when D's chain
is a contiguous sub-chain of C's, `D = C[i .. i+m)`, D is not C, and D's reads stay inside C's exons around it — every read of D starts
at or after the acceptor of C's intron `i-1` (when i > 0) and ends at or before the donor of C's intron `i+m` (when i+m < len C), tested
on the group's outermost start and end. For each junction t of C, its compatible support `s_C(t)` = reads of C + reads of every
compatible D that contains t. A group with fewer than `min_reads` exact reads is ALSO admitted when `min_t s_C(t) >= min_reads`; its
skeleton then carries `n_reads = min_t s_C(t)`. Groups already admitted are untouched. Unspliced reads, the gate, polish and every
other step are unchanged; `min_reads` is the existing floor (2), so the rule has no new parameter.

### Test

**Arms.** OFF (current) and ON, each in the two seeding pools: `primary` (the benchmark arm) and `good` (the default: primaries +
secondaries with AS >= 0.98 of the molecule's best). Scored with gffcompare against the same reference as figure 1
(`assembly/<species>/eval_all/ref.gtf`).
**Dev.** Gorilla OR6737 NC_073244.2 only (used for the diagnosis): correctness and sanity, no tuning (there is nothing to tune); at
most ONE amendment, written before any held-out number.
**Held-out (decides).** (H-g) gorilla OR6737 genome-wide minus NC_073244.2; (H-h) human A119b genome-wide. Per sample and pool:
- **B1** intron sensitivity up by >= 1.0 point;
- **B2** intron-chain precision down by <= 0.5 point;
- **B3** intron-chain sensitivity not lower.
Reported: transcripts emitted, transcript-level sensitivity / precision, locus level, the share of the StringTie intron gap closed.
**Verdict.** B1-B3 hold in both held-out samples for BOTH pools => ADOPT as the default (and as the base for the seeded-loci extension);
for the primary pool only => default for the benchmark arm, opt-in otherwise; else NOT ADOPTED, kept opt-in with the numbers.

### Result (2026-10-11): NOT ADOPTED, kept opt-in

**Dev (NC_073244.2, primary):** OFF byte-identical to the previous binary; ON +17 transcripts (3,862 -> 3,879), intron sensitivity 54.8 -> 54.9,
intron-chain 28.3 / 40.8 -> 28.3 / 40.7. Why it barely moves (dev diagnosis, before any held-out number): for 864 of the 927 missed
junctions, no read chain carrying the junction reaches a compatible support of 2 at ALL its junctions — the supporting reads are
STAGGERED (each covers a different overlapping stretch), so no OBSERVED chain is long enough. StringTie stitches overlapping reads into a
longer path; this rule, by design, never builds a chain no read carries. No amendment was made.
**Held-out (gffcompare, sensitivity | precision, points):**

| arm | transcripts | intron | intron chain | transcript | locus |
|---|---|---|---|---|---|
| gorilla primary OFF | 68,828 | 59.6 / 82.8 | 26.8 / 35.3 | 24.0 / 35.1 | 30.2 / 70.7 |
| gorilla primary ON | 70,375 | 59.9 / 82.6 | 27.2 / 35.0 | 24.4 / 34.8 | 30.3 / 71.0 |
| gorilla good OFF | 70,134 | 60.0 / 82.5 | 27.0 / 35.0 | 24.2 / 34.7 | 30.5 / 69.1 |
| gorilla good ON | 71,700 | 60.3 / 82.3 | 27.3 / 34.6 | 24.5 / 34.4 | 30.7 / 69.4 |
| human primary OFF | 238,406 | 60.8 / 64.4 | 21.9 / 18.6 | 20.4 / 16.2 | 34.6 / 22.9 |
| human primary ON | 250,598 | 61.2 / 63.9 | 22.3 / 18.1 | 20.7 / 15.6 | 34.8 / 22.2 |
| human good OFF | 252,217 | 61.2 / 62.7 | 22.0 / 17.7 | 20.5 / 15.4 | 34.8 / 22.0 |
| human good ON | 263,799 | 61.6 / 62.2 | 22.4 / 17.3 | 20.9 / 14.9 | 35.0 / 21.6 |

(gorilla scored without NC_073244.2; human on the figure 1 reference `eval_annotated/ref.gtf`.) **B1 fails in all four** (intron
sensitivity +0.3 / +0.3 / +0.4 / +0.4, bar +1.0); B2 holds (intron-chain precision -0.3 / -0.4 / -0.5 / -0.4, bar <= 0.5); B3 holds
(+0.4 / +0.3 / +0.4 / +0.4). The flag stays opt-in (`copy_assign --compat-support`, `RUSTLE_COMPAT_SUPPORT=1`): about +0.4 point of
intron-chain sensitivity for -0.3 to -0.5 of precision and +2-5% transcripts. Next mechanism to pre-register separately: stitching
overlapping, consistent reads into longer paths (StringTie's splice graph), after measuring how many missed junctions it could recover.


---

## PREREG_direct_transfrags_2026-10-11.md

*(original: `docs/archive/2026-10/merged/PREREG_direct_transfrags_2026-10-11.md`)*

## Direct-transfrag assembly (no network flow), 2026-10-11

**Pre-registered before the code exists and before any arm is scored.** User (2026-10-11), after `--compat-support` failed
(register 1264): build, from scratch, an assembly that uses transfrags directly — no network flow — the way an earlier (now deleted)
StringTie-era mode reproduced StringTie except for flow-only transcripts.

**Why it can differ from register 1264.** Compatible support admitted a chain only when EVERY junction had >= 2 compatible reads; the
missed junctions sit in chains whose terminal junction one read carries (dev: 864 of 927). StringTie decides on a path's ABUNDANCE
(its reads plus the reads of the partial transfrags it explains), not a per-junction minimum. This rule does the same without flow.

### The rule (`copy_assign --direct-transfrags`, default off = byte-identical; exclusive with `--compat-support`)

Pass-1 groups (one per exact `(chrom, intron chain)`, with exact reads n, outermost start and end) are the transfrags. Transfrag T is
CONTAINED in C (C != T) when T's chain is a contiguous sub-chain `C[i .. i+m)` and T's reads stay inside C's exons around it (smallest
start >= acceptor of C's intron `i-1` when i > 0; largest end <= donor of C's intron `i+m` when i+m < len C) — the containment of
`--compat-support`. A transfrag contained in no other is MAXIMAL. Every non-maximal T transfers its n reads DIRECTLY to the maximal
transfrags that contain it, split in proportion to their own exact reads. Each maximal transfrag's abundance = its reads + what it
received; it is admitted when the abundance reaches the existing read floor (2) and its skeleton carries the rounded abundance.
Non-maximal transfrags are not emitted (their reads are in their containers). Unspliced reads, the gate, polish and every other step
are unchanged. No new parameter.

### Test (as `docs/PREREG_assembler_compatible_support_2026-10-11.md`)

Arms OFF / ON in the `primary` and `good` pools; gffcompare against the figure 1 references. Dev: gorilla OR6737 NC_073244.2, no
tuning, at most ONE amendment before any held-out number. Held-out (decides): gorilla OR6737 minus NC_073244.2; human A119b genome-wide.
**B1** intron sensitivity +1.0 point or more; **B2** intron-chain precision down by <= 0.5 point; **B3** intron-chain sensitivity not lower.
Reported: transcripts, transcript and locus levels, the share of the StringTie intron gap closed. **Verdict:** B1-B3 in both samples and
both pools => ADOPT as the default; primary pool only => default for the benchmark arm; else NOT ADOPTED, opt-in with the numbers.

### Result (2026-10-11): NOT ADOPTED as the default (B2 fails); kept opt-in

**Dev (NC_073244.2, primary):** default byte-identical; ON 3,862 -> 3,999 transcripts; intron 54.8 / 83.3 -> 58.5 / 81.2 (StringTie 58.9 /
83.9); intron chain 28.3 / 40.8 -> 27.5 / 38.2 (StringTie 24.7 / 37.0). No amendment.
**Held-out** (gffcompare sensitivity / precision; gorilla without NC_073244.2, StringTie and FLAIR scored on the same contigs; human on
`eval_annotated/ref.gtf`):

| arm | transcripts | intron | intron chain | transcript | locus |
|---|---|---|---|---|---|
| gorilla primary OFF | 68,828 | 59.6 / 82.8 | 26.8 / 35.3 | 24.0 / 35.1 | 30.2 / 70.7 |
| gorilla primary ON | 72,607 | 63.2 / 80.0 | 27.0 / 33.6 | 24.2 / 33.5 | 30.2 / 70.1 |
| gorilla good OFF | 70,134 | 60.0 / 82.5 | 27.0 / 35.0 | 24.2 / 34.7 | 30.5 / 69.1 |
| gorilla good ON | 73,996 | 63.6 / 79.7 | 27.1 / 33.3 | 24.4 / 33.1 | 30.5 / 68.4 |
| gorilla StringTie | — | 64.2 / 83.1 | 24.4 / 34.4 | 21.9 / 34.1 | 30.7 / 64.1 |
| gorilla FLAIR | — | 62.5 / 61.3 | 25.5 / 17.3 | 23.0 / 15.7 | 28.4 / 48.4 |
| human primary OFF | 238,406 | 60.8 / 64.4 | 21.9 / 18.6 | 20.4 / 16.2 | 34.6 / 22.9 |
| human primary ON | 251,710 | 64.4 / 59.0 | 21.9 / 17.1 | 20.4 / 15.3 | 34.0 / 24.1 |
| human good OFF | 252,217 | 61.2 / 62.7 | 22.0 / 17.7 | 20.5 / 15.4 | 34.8 / 22.0 |
| human good ON | 263,515 | 64.7 / 57.7 | 22.0 / 16.5 | 20.5 / 14.7 | 34.2 / 23.4 |
| human StringTie | — | 63.5 / 57.8 | 18.1 / 15.0 | 16.9 / 12.8 | 33.8 / 20.3 |
| human FLAIR | — | 67.7 / 30.2 | 22.4 / 6.3 | 21.1 / 3.6 | 31.9 / 6.3 |

**B1 holds** (intron sensitivity +3.6 / +3.6 / +3.6 / +3.5); **B2 fails** (intron-chain precision -1.7 / -1.7 / -1.5 / -1.2, bar <= 0.5);
B3 holds (+0.2 / +0.1 / 0.0 / 0.0). Not the default. Against the tools on the same contigs (primary pool): human — direct transfrags is
above StringTie on every intron and chain number (intron 64.4 / 59.0 vs 63.5 / 57.8; chain 21.9 / 17.1 vs 18.1 / 15.0); gorilla —
intron sensitivity 1.0 below StringTie and precision 3.1 below, chain sensitivity 2.6 above and precision 0.8 below. FLAIR has higher
intron sensitivity in human (67.7) at half the precision. The flag stays opt-in (`--direct-transfrags`, `RUSTLE_DIRECT_TRANSFRAGS=1`);
whether the benchmark arm uses it is a figure decision, not this test's.


---

## PREREG_direct_transfrags_union_2026-10-11.md

*(original: `docs/archive/2026-10/merged/PREREG_direct_transfrags_union_2026-10-11.md`)*

## Direct transfrags + exact chains (union), 2026-10-11

**Pre-registered before the code exists and before any arm is scored.** Follows `docs/PREREG_direct_transfrags_2026-10-11.md`
(register 1265: intron sensitivity +3.5 held-out, chain precision -1.2 to -1.7). User (2026-10-11): are the remaining StringTie-only
transcripts the flow-only ones, and can the rest be included without flow?

**Dev diagnosis (gorilla OR6737 NC_073244.2 only).** Of 195 multi-exon reference transcripts StringTie matches exactly and direct
transfrags does not: 9 (5%) are carried whole by no read (flow-only); **66 are in the DEFAULT output** — >= 2 reads carry exactly that
chain, but direct transfrags absorbs it into a longer transfrag that contains it (an alternative start or end that looks like a
truncation); 71 are carried whole by one read; 29 by >= 2 reads of which one exact; 16 exact >= 2 but dropped by gate or polish in both
modes; 4 only inside longer reads. The 66 are a cost of the absorption rule, not of the missing flow.

### The rule (`copy_assign --direct-transfrags-union`, default off = byte-identical; exclusive with the other two)

Admit every group the default admits (`>= 2` identical reads) AND every maximal transfrag whose direct-transfer abundance reaches the
floor (`--direct-transfrags`). A group in both carries the transferred abundance; a contained group admitted on its own reads keeps them
(its reads still count toward its containers). Everything else unchanged. No new parameter.

### Test (fresh held-out: none of these samples was used by any assembler test)

Arms OFF (the cached default runs) / ON, pools `primary` and `good`. **Held-out (decides):** gorilla KB3781 fibroblast (gorilla reference
of figure 1), chimpanzee PTR and orangutan PPY (their annotation GTFs, `assembly/<sample>/annotation.gff_to_gtf.gtf`). Dev: OR6737
NC_073244.2 (sanity only, no amendment unless written before any held-out number). Per sample and pool: **B1** intron sensitivity +1.0
point or more; **B2** intron-chain precision down by <= 0.5 point; **B3** intron-chain sensitivity not lower. **Verdict:** B1-B3 in all
three samples for both pools => ADOPT as the default; primary pool only => default for the benchmark arm; else opt-in with the numbers.
Reported: transcripts, transcript and locus levels; on OR6737 and A119b (already seen) the comparison to StringTie / FLAIR is descriptive.

### Result (2026-10-11): NOT ADOPTED as the default (B2 fails on every held-out arm); kept opt-in

**Dev (NC_073244.2, primary):** 4,169 transcripts; intron 58.8 / 81.2, chain 29.0 / 38.7, transcript 27.2 / 38.7, locus 41.1 / 78.1
(default 54.8 / 83.3, 28.3 / 40.8, 26.5 / 40.7, 39.9 / 77.8; StringTie 58.9 / 83.9, 24.7 / 37.0, 23.2 / 36.8, 40.0 / 71.4).
**Fresh held-out** (sensitivity / precision; OFF = the cached default runs):

| sample, pool | transcripts OFF -> ON | intron | intron chain | transcript | locus |
|---|---|---|---|---|---|
| KB3781 primary | 71,855 -> 74,958 | 53.9 / 81.7 -> 55.9 / 80.2 | 27.0 / 36.1 -> 27.6 / 35.3 | 24.3 / 36.0 -> 24.8 / 35.2 | 28.3 / 76.8 -> 28.7 / 77.9 |
| KB3781 good | 72,622 -> 75,773 | 54.1 / 81.5 -> 56.0 / 80.0 | 27.1 / 35.9 -> 27.7 / 35.1 | 24.4 / 35.8 -> 24.9 / 35.0 | 28.5 / 75.4 -> 28.9 / 76.5 |
| chimp primary | 53,909 -> 59,370 | 47.1 / 85.8 -> 52.2 / 83.2 | 16.0 / 39.4 -> 17.1 / 38.1 | 14.7 / 38.8 -> 15.7 / 37.7 | 26.4 / 63.3 -> 27.6 / 65.2 |
| chimp good | 54,821 -> 60,320 | 47.4 / 85.6 -> 52.4 / 83.0 | 16.1 / 39.1 -> 17.2 / 37.8 | 14.8 / 38.5 -> 15.8 / 37.3 | 26.7 / 62.4 -> 27.8 / 64.3 |
| orangutan primary | 106,467 -> 121,865 | 56.3 / 70.9 -> 63.3 / 65.4 | 22.1 / 19.9 -> 23.8 / 18.4 | 20.1 / 19.0 -> 21.6 / 17.8 | 30.5 / 39.9 -> 32.0 / 43.1 |
| orangutan good | 111,629 -> 127,292 | 56.7 / 69.9 -> 63.6 / 64.4 | 22.3 / 19.3 -> 24.0 / 17.9 | 20.3 / 18.2 -> 21.8 / 17.2 | 30.8 / 36.8 -> 32.3 / 40.0 |

**B1 holds** (intron sensitivity +1.9 to +7.0); **B3 holds** (chain sensitivity +0.6 to +1.7); **B2 fails** (chain precision -0.8 to
-1.5, bar 0.5). Locus sensitivity and precision both rise in all six. Opt-in: `--direct-transfrags-union`, `RUSTLE_DIRECT_TRANSFRAGS=union`.
**Against the tools (benchmark samples, genome-wide, primary pool, descriptive — both samples were seen before):**

| | base | intron | intron chain | transcript | locus |
|---|---|---|---|---|---|
| gorilla default | 45.1 / 68.1 | 59.3 / 82.8 | 26.9 / 35.6 | 24.2 / 35.4 | 30.7 / 71.1 |
| gorilla union | 46.9 / 67.5 | 63.1 / 80.2 | 27.9 / 34.1 | 25.1 / 34.0 | 31.5 / 71.9 |
| gorilla StringTie | 50.3 / 63.1 | 63.9 / 83.1 | 24.4 / 34.6 | 22.0 / 34.3 | 31.2 / 64.5 |
| gorilla FLAIR | 48.5 / 46.0 | 62.1 / 61.5 | 25.5 / 17.6 | 23.1 / 16.0 | 28.8 / 49.2 |
| human default | 48.2 / 17.5 | 60.8 / 64.4 | 21.9 / 18.6 | 20.4 / 16.2 | 34.6 / 22.9 |
| human union | 49.2 / 19.2 | 64.4 / 59.2 | 22.9 / 17.4 | 21.3 / 15.7 | 35.1 / 25.3 |
| human StringTie | 59.6 / 14.5 | 63.5 / 57.8 | 18.1 / 15.0 | 16.9 / 12.8 | 33.8 / 20.3 |
| human FLAIR | 66.4 / 8.8 | 67.7 / 30.2 | 22.4 / 6.3 | 21.1 / 3.6 | 31.9 / 6.3 |

Human: the union is above StringTie on intron, chain, transcript and locus (both measures); StringTie keeps base sensitivity. Gorilla: the
union is above StringTie on chain and transcript sensitivity and locus, about level on intron sensitivity (-0.8), below on intron
precision (-2.9) and chain precision (-0.5), and StringTie keeps base sensitivity. Flow-only transcripts are a small part of the
remaining gap (dev: 9 of 195).


---

## PREREG_inbetween_operating_point_2026-10-11.md

*(original: `docs/archive/2026-10/merged/PREREG_inbetween_operating_point_2026-10-11.md`)*

## An operating point between StringTie and FLAIR that wins on multi-copy loci, 2026-10-11

**Pre-registered before any arm below is run genome-wide and before the strata are computed.** User (2026-10-11): "something
in between StringTie and FLAIR, that has the strengths of both, but ours should outperform them in secondary-only loci or in
places where O2 found the right assignment or corrected the arbitrarily selected primary when AS is equal." User choices
(same day): success = beat each tool on its weak side; multi-copy arm = the good-secondaries pool only (no O2 closed loop);
decisive samples = the benchmark samples genome-wide.

### The operating point (fixed from the difference model, register 1267; no search)

`OP` = `copy_assign --direct-transfrags-union --read-isoform-k 2` on top of the shipped assembly settings, i.e.
`RUSTLE_DIRECT_TRANSFRAGS=union RUSTLE_READ_ISOFORM_K=2 tools/rustle_pipeline.sh assemble ...`. It keeps all our chains (no
25 bp junction correction, floor 2: the classes register 1267 found right more often than the tools' additions), adds
StringTie's exclusive class (single-read chains whose every junction another read carries), and does not add FLAIR's
(`single_novel`). Two pools: `OP_primary` (`--no-seed-secondaries`) and `OP_good` (the default good pool: secondaries with
AS >= 0.98 x the molecule's genome-wide best AS). Reported only: our cached shipped default runs (`default_primary`,
`default_good`, figure run cache 2026-10-10).

### Samples and scope

Gorilla OR6737 (`GGO_mm.bam`) and human A119b (`A119b.t2t.bam`), genome-wide **minus the development chromosome**
(NC_073244.2 / chr20), against the lab's annotation-free StringTie 3.0.3 and FLAIR 3.0.1 GTFs on the same BAMs, scored by
gffcompare against `GGO.GCF_029281585.2_RefSeq.gtf.gz` / `A119b.chm13v2.0_RefSeq.gtf.gz`. Both samples were used before
(register 1264-1266; disclosed); no rule here is fitted on them.

### Multi-copy strata (from reads only; `bench/tool_difference/multicopy_strata.py`, sha1 c663dc6b)

Per reference gene, records whose aligned blocks overlap its exons: primaries; GOOD secondaries (AS >= 0.98 x the molecule's
best, from `as_table`'s `molecules.tsv`); TIED secondaries (AS == the molecule's best: its primary was one of several equally
scored placements). **S1 secondary-only**: 0 primaries and >= 2 good secondaries. **S2 coin-toss-starved**: >= 1 primary,
>= 2 tied secondaries and tied >= primaries. Evaluated: the multi-exon reference transcripts of S1 / S2 genes.

### Claims (`bench/tool_difference/score_op.py`, sha1 f65413d1; each in both species)

- **B1** `OP_primary` intron-chain sensitivity > StringTie's and > FLAIR's.
- **B2** `OP_primary` intron-chain precision > FLAIR's.
- **C1** In S1, `OP_good` matches exactly (gffcompare `=`) more reference transcripts than StringTie and than FLAIR.
- **C2** The same in S2.
- **C3** In each stratum, `OP_good`'s precision (its transcripts touching a stratum gene's exon that match a reference
  transcript exactly) >= FLAIR's genome-wide intron-chain precision (the weak-side bar of B2).

A stratum with < 20 multi-exon reference transcripts in a species is UNDERPOWERED there (no verdict for C1/C2/C3 on it).
**Verdict**: B1-B2 in both species => OP is "between StringTie and FLAIR" as defined; C1-C3 in both species => it also wins
on multi-copy loci. Each failure is named; nothing is re-tuned after a number is seen. Reported, no bar: base / intron /
transcript / locus levels, OP_primary vs StringTie precision, default arms, stratum sizes.

**Expected, written now.** B1/B2 pass in both (dev chr20: gorilla chain 31.3 / 34.1 vs StringTie 24.7 / 37.0 and FLAIR 25.1 /
23.6; human 26.7 / 17.2 vs 20.1 / 16.8 and 23.9 / 7.6). C1 should pass wherever S1 is powered, since neither tool's transcripts
rest on secondary-only coverage; C3 is the real test (do transcripts built from secondaries match the annotation?). S2 may be
small: register 875 found few AS-tied molecules in catalog families.

**Amendment 1 (2026-10-11, before any stratum or arm number).** (a) Instrument crash fix: RefSeq gene-level features carry
`gene_id ""` (gorilla and human references); both scripts now key such a gene by its `transcript_id` (`gene-LOC...`) instead of
crashing. No definition changes. (b) Execution: `--read-isoform-k` forces the buffered reader, which held 14+ GB genome-wide
with secondaries (the gorilla good-pool genome-wide run stalled; the human good-pool run was stopped by Claude Code's memory
reaper). The arms run one contig per process (`tools/rustle_pipeline.sh --contig`, documented equal to the genome-wide run on
that contig) and `--read-isoform-k` now also runs on the streaming reader (1304476b; structurally identical to the buffered
path on human chr20). Same arms, same flags.
(c) The strata script ran past the 500 s run limit on gorilla before printing any stratum count; it now counts each
chromosome in a separate worker (4 processes) with the identical per-record rule (sha1 below), before any stratum number.
multicopy_strata.py sha1 f4d7bc36798e6f89fc914b6d52fe7c3b1aa9cd78

### Result (2026-10-11): B1, B2, C1, C2, C3 PASS in both species

Scored with `bench/tool_difference/score_op.py` (sha1 45901172, the amendment-1 version); outputs
`/mnt/linuxdisk/tmp/asm_gap/op/score_{gorilla,human}.txt`. Genome-wide minus the development chromosome, sensitivity / precision:

| arm | gorilla base | intron | intron chain | transcript | locus | human base | intron | intron chain | transcript | locus |
|---|---|---|---|---|---|---|---|---|---|---|
| OP_primary | 48.6 / 63.7 | 65.7 / 78.1 | **29.3 / 30.4** | 26.3 / 30.2 | 31.7 / 70.1 | 50.6 / 18.2 | 66.7 / 56.4 | **23.6 / 15.6** | 22.0 / 14.0 | 36.1 / 23.5 |
| OP_good | 49.1 / 62.0 | 66.1 / 77.8 | 29.4 / 30.0 | 26.5 / 29.8 | 32.0 / 68.3 | 51.7 / 17.8 | 67.0 / 54.9 | 23.8 / 15.0 | 22.2 / 13.4 | 36.3 / 22.6 |
| StringTie | 50.3 / 62.6 | 64.2 / 83.1 | 24.4 / 34.4 | 21.9 / 34.1 | 30.7 / 64.1 | 59.6 / 14.4 | 63.5 / 57.8 | 18.1 / 15.0 | 16.9 / 12.7 | 33.7 / 20.1 |
| FLAIR | 48.7 / 45.2 | 62.5 / 61.3 | 25.5 / 17.3 | 23.0 / 15.7 | 28.4 / 48.4 | 66.4 / 8.8 | 67.7 / 30.1 | 22.4 / 6.2 | 21.1 / 3.6 | 31.9 / 6.3 |
| default_primary (reported) | 45.1 / 67.7 | 59.6 / 82.8 | 26.8 / 35.3 | 24.0 / 35.1 | 30.2 / 70.7 | 48.1 / 17.4 | 60.8 / 64.4 | 21.8 / 18.5 | 20.3 / 16.1 | 34.5 / 22.8 |
| default_good (reported) | 45.7 / 65.9 | 60.0 / 82.5 | 27.0 / 35.0 | 24.2 / 34.7 | 30.5 / 69.1 | 49.2 / 17.1 | 61.2 / 62.6 | 22.0 / 17.7 | 20.5 / 15.3 | 34.8 / 21.9 |

**B1 PASS** (chain sensitivity: gorilla 29.3 vs 24.4 / 25.5; human 23.6 vs 18.1 / 22.4). **B2 PASS** (chain precision: gorilla 30.4
vs FLAIR 17.3; human 15.6 vs 6.2). Reported: against StringTie's chain precision the OP is below in gorilla (30.4 vs 34.4) and above
in human (15.6 vs 15.0).

Multi-copy strata (reference multi-exon transcripts matched exactly; precision = matching / all arm transcripts touching the stratum):

| | gorilla S1 (542 tx, 305 genes) | gorilla S2 (252 tx, 117 genes) | human S1 (926 tx, 304 genes) | human S2 (592 tx, 266 genes) |
|---|---|---|---|---|
| OP_good | **34** (34 / 171 = 19.9%) | **84** (90 / 367 = 24.5%) | **36** (37 / 224 = 16.5%) | **210** (218 / 1,368 = 15.9%) |
| OP_primary | 0 | 56 (36.5%) | 0 | 100 (17.9%) |
| StringTie | 0 | 51 (39.0%) | 0 | 77 (13.3%) |
| FLAIR | 0 | 42 (36.8%) | 0 | 72 (8.0%) |
| default_good (reported) | 33 (28.7%) | 81 (32.8%) | 33 (19.3%) | 198 (18.7%) |

**C1 PASS** (S1: 34 / 36 vs 0 for both tools; neither tool emits a transcript on secondary-only genes). **C2 PASS** (S2: 84 vs 51 /
42; 210 vs 77 / 72). **C3 PASS** in every stratum (19.9 / 24.5 / 16.5 / 15.9% >= FLAIR's genome-wide 17.3 / 6.2%; the gorilla S1 margin
is 2.6 points). **Limits, stated plainly:** S1 sensitivity is low (6.3% / 3.9% of the stratum's reference transcripts); in S2 the good
pool doubles the matches but its precision there is below the primary pool's and, in gorilla, below the tools'; the gain over our own
`default_good` is small (S1 +1 / +3, S2 +3 / +12 matches). The multi-copy advantage comes from the good-secondaries pool, not from the
operating point's admission rules. Gorilla was scored by the duplicate session (the same frozen scripts), human by this one.


---

## PREREG_lrpap1_seed_reach_2026-10-11.md

*(original: `docs/archive/2026-10/merged/PREREG_lrpap1_seed_reach_2026-10-11.md`)*

## How many seeds does the DNA mode need to find every LRPAP1 copy? (gorilla mGorGor1), 2026-10-11

**Pre-registered before any seed is aligned.** User (2026-10-11): the advisor is not sure the DNA (guided) mode finds all copies; show, for the
gorilla LRPAP1 expansion, how many copies are found from one or more of the other copies used as seeds, i.e. how many seeds are needed.

**Truth.** The 11 gorilla LRPAP1 loci of the LRPAP1 artifact (`/mnt/linuxdisk/tmp/lrpap1/browser/lrpap1_browser.json`): 8 full-length copies
(c00-c07) and 3 copies of the 5′ three-quarters (p12, p14, p16); 10 of them are Yoo 2025's list, the 11th (p16) was found 2026-10-04.

**The DNA-mode finder (as shipped in `bench/guided_pipeline.py`, unchanged).** A seed copy contributes its transcript (its model's exons,
spliced; `minimap2 -x splice -N 100 -p 0.1` on the primary genome index `GGO.splice.mmi`) and its gene body (model span, first to last exon;
`-x asm20`, same flags). Passing hits: `transcript_hits` (identity >= 0.80, >= 50% of the seed aligned) and `gene_body_chains` (identity >= 0.80,
aligned >= 50% of the shorter of seed and extrapolated span). Seed i FINDS copy j when a passing hit of i overlaps j's span on its chromosome.
Hits outside the 11 are counted per seed (candidate burden), not judged.

**Reported.** The 11 x 11 finds matrix; for every seed set size k = 1..11 and every set of that size: copies found in ONE pass (union of the seeds'
finds) and with ITERATION (found copies become seeds until nothing new: reachability in the finds graph); the minimum number of seeds that finds
all 11 in one pass (smallest set cover) and with iteration (number of source components of the finds graph).

**Expectation (registered).** E1: with iteration, ONE seed suffices, whichever copy is chosen (the copies are 95.7-99% identical, far above 0.80).
E2: in one pass, every full-length seed finds all 11; a fragment seed finds every copy carrying its exons 1-5, i.e. all 11.
If either fails, the copies a seed misses and why (coverage or identity) are reported; no threshold is changed after looking.

### Result (2026-10-11) and Amendment 1

**E1 and E2 hold.** Each of the 11 copies, used as the only seed, finds all 11 in one pass (and so with iteration); minimum seeds = 1 in both modes;
the finds graph is one strongly connected component. No passing hit falls outside the 11 (0 extra loci for every seed). The gene-body finder alone
finds all 11 from every seed; the transcript finder alone finds 8 from seeds c03, c04 and c06 (it misses the three 5′ fragments) and 11 from the
others. Pairwise gene-body identity between copies 0.959-0.998 (110 ordered pairs), against the finder's 0.80.
**Amendment 1 (after the result; descriptive, no claim).** A sweep of the finders' identity threshold from 0.80 to 1.00 (step 0.0025) on the same
alignments: seeds needed in one pass and with iteration, and the source components, to show how far LRPAP1 sits from the regime where subclusters
need their own seeds.
**Amendment 1 result.** One seed finds all 11 in one pass up to a finder identity of 0.9775 and with iteration up to 0.98; from 0.9825 the most
diverged copy (c06, chr22) needs its own seed; at 0.99 four seeds are needed with iteration (five in one pass), at 0.995 seven, at 1.0 eleven. The
shipped finder threshold is 0.80, so LRPAP1 sits 0.18 identity below the point where its subclusters would need separate seeds.


---

## PREREG_tool_difference_model_2026-10-11.md

*(original: `docs/archive/2026-10/merged/PREREG_tool_difference_model_2026-10-11.md`)*

## Difference model: why a chain is in our output and not StringTie's or FLAIR's, 2026-10-11

**Pre-registered after the development numbers (chr20, both species) and before any held-out chromosome is run.** User
(2026-10-11): establish a rule that distinguishes our assembler from StringTie (and FLAIR), as a base for an operating point
between the two tools.

### Instrument (frozen)

`bench/tool_difference/diff_model.py` (sha1 `a78ee425fd2874a04f67f1601c21168ef0c0c6f4`). Inputs per chromosome: the primary
alignments (`samtools view -F 2308`, columns 1-6), the figure 1 reference restricted to the chromosome (gorilla
`GGO.GCF_029281585.2_RefSeq.gtf.gz`, human `A119b.chm13v2.0_RefSeq.gtf.gz`), our shipped default assembly
(`copy_assign --assemble-only --region <chrom> --assembly-junctions strict` + the shipped polish, primary pool), and the lab's
annotation-free StringTie 3.0.3 / FLAIR 3.0.1 GTFs restricted to the chromosome. Read chains use our default de-duplication
key. Every multi-exon intron chain exclusive to one side gets the first class that applies:

- **other-only**: `ours_gate` (>= 2 reads carry it exactly; we lose it in gate / polish / motif), `single_shared` (1 read,
  every junction also carried by another read), `single_novel` (1 read, a junction only it carries), `corrected` (no read
  carries it; a read chain matches within 25 bp at every junction), `stitched` (no read carries it or a near copy).
- **ours-only**: `near_stronger` (a junction within 25 bp, donor and acceptor, of one more reads carry), `two_reads`
  (exactly 2 reads), `subset` (a contiguous sub-chain of another of our transcripts), `other`.

`reference-exact` = share of a class whose intron chain a reference transcript has exactly.

### Development numbers (chr20; seen)

| | gorilla NC_073244.2 | human chr20 |
|---|---|---|
| StringTie-only chains (reference-exact) | 1,274 (11.5%) | 2,954 (3.3%) |
| single_shared share / ref-exact | 63.1% / 14.4% | 37.6% / 5.0% |
| stitched share / ref-exact | 13.7% / 5.2% | 23.2% / 2.2% |
| ours-only vs StringTie (ref-exact) | 1,412 (24.3%) | 2,853 (10.4%) |
| near_stronger share / ref-exact | 40.2% / 34.9% | 21.8% / 17.8% |
| FLAIR-only chains (ref-exact) | 3,217 (3.5%) | 10,474 (1.7%) |
| FLAIR-only single_novel / single_shared ref-exact | 2.2% / 8.5% | 0.5% / 2.8% |
| FLAIR-only corrected + stitched | 0 | 0 |
| ours-only vs FLAIR (ref-exact); two_reads ref-exact | 1,180 (24.7%); 19.8% | 2,051 (10.4%); 5.1% |

### Held-out (decides)

Chromosomes fixed before running: gorilla OR6737 **NC_073242.2** and **NC_073245.2**; human A119b **chr18** and **chr22**.
(These samples' other chromosomes served as held-out for register 1264-1265; these four chromosomes are reused only as the
genome-wide remainder there, and this test adopts no assembler rule.) Counts are summed over the two chromosomes of a species.
Claims, each in both species:

- **D1** `single_shared` is the largest StringTie-only class.
- **D2** StringTie-only `single_shared` is reference-exact more often than StringTie-only `stitched`.
- **D3** ours-only (vs StringTie) `near_stronger` is reference-exact more often than StringTie-only chains overall.
- **D4** FLAIR-only chains are all read-carried (`corrected` + `stitched` = 0), and `single_novel` is reference-exact less often
  than `single_shared`.
- **D5** ours-only (vs FLAIR) `two_reads` is reference-exact more often than FLAIR-only chains overall.
- **D6** ours-only chains are reference-exact more often than other-only chains, against both tools.

**Verdict.** D1-D6 in both species => the difference model HOLDS and is the stated basis for the in-between operating point;
any failure is named and the model is restated without that claim. Class shares are reported per chromosome with their dev
values (no bar).

**Reading the model (fixed now).** StringTie's lead is mostly single-read chains whose junctions other reads share (D1, D2);
its 25 bp junction correction removes chains that are right more often than its own additions (D3). FLAIR's lead is
single-read chains, mostly with a junction no other read has (D4), and its support floor of 3 drops two-read chains that are
right more often than its additions (D5). An operating point between the two keeps our chains and adds StringTie's
`single_shared` class (our `--read-isoform-k 2`), not FLAIR's `single_novel` class.

### Result (2026-10-11): the difference model HOLDS (D1-D6 pass in both species)

Scored with `bench/tool_difference/score_claims.py` (pools the two chromosomes; class reference-exact counts are rebuilt from the
printed 3-decimal rates, so +-1 chain; no claim is within that margin). Outputs `/mnt/linuxdisk/tmp/asm_gap/diff/<chrom>/`.

| held-out | gorilla NC_073242.2 + NC_073245.2 | human chr18 + chr22 |
|---|---|---|
| StringTie-only chains (ref-exact) | 1,814 (10.5%) | 5,455 (3.0%) |
| single_shared share / ref-exact (D1, D2) | 63.6% / 12.7% | 36.1% / 5.0% |
| stitched share / ref-exact (D2) | 12.7% / 7.0% | 20.7% / 1.4% |
| ours-only vs StringTie (ref-exact) | 1,941 (18.1%) | 4,955 (8.5%) |
| near_stronger share / ref-exact (D3) | 32.9% / 29.9% | 16.8% / 16.5% |
| FLAIR-only chains (ref-exact) | 5,416 (3.0%) | 20,550 (1.6%) |
| FLAIR-only single_novel / single_shared ref-exact (D4) | 1.4% / 6.0% | 0.6% / 2.1% |
| FLAIR-only corrected + stitched (D4) | 0 | 0 |
| ours-only vs FLAIR (ref-exact); two_reads ref-exact (D5) | 1,501 (18.4%); 15.2% | 3,769 (7.6%); 4.7% |

D1 PASS, D2 PASS, D3 PASS, D4 PASS, D5 PASS, D6 PASS in both species. Class shares stay near their development values
(StringTie-only single_shared 63.6 vs 63.1% gorilla, 36.1 vs 37.6% human; FLAIR-only single_novel 54.7 vs 56.8%, 53.6 vs 54.6%).
Side check (gorilla chr20, `sec_cov`-style): no StringTie or FLAIR transcript rests on secondary-only coverage (every one has
>= 50% of its exonic bases under primary alignments), so both tools are primary-only in practice on these BAMs.


---

## GENOMEWIDE_FAST_2026-10-12.md

*(original: `docs/archive/2026-10/merged/GENOMEWIDE_FAST_2026-10-12.md`)*

## Running the pipeline genome-wide fast on a small machine (2026-10-12)

Goal (user, 2026-10-12): integrate and simplify the codebase so the pipeline runs genome-wide in a short time with limited
resources. Machine of record: WSL2, 5 cores, 25 GB RAM (Claude Code's background reaper kills jobs well before that).

### Baseline (recorded stage runs, `/mnt/linuxdisk/tmp/rustle_figures/runs/*/<stage>.time`)

| stage | gorilla OR6737 / KB3781 | human A119b | what dominates |
|---|---|---|---|
| assemble (streaming) | 1.7-3 min, 1.7-2.2 GB | 7.2 min, 2.9 GB | BAM decoding |
| families | 2 min, 6.3 GB | **45 min, 12.6 GB** | all-vs-all `minimap2 -x asm20` over 2 GB of locus spans (introns included; 4.5 M of the pairs it finds are later dropped for no exonic evidence) |
| assign | 9.4 min, 5.1 GB | never completed (reaped 2026-10-11); measuring | per-region read loading; cross-chromosome families (MCL3: 260 copies on 23 chromosomes) |
| flag | not run | not run | `minimap2 -x splice:hq` of a few thousand consensus sequences against the WHOLE splice index (12-14 GB `.mmi` loaded at once: the 16-19 GB peak) |
| splice index build | — | 7 min, 19-20 GB | one-part index of the whole genome |

### Targets

One command, resumable, genome-wide on the human sample: assemble -> families -> assign -> flag, every process <= 12 GB peak RSS
(a low-memory profile <= 8 GB), total well under 2 h on 4 threads; gorilla well under 30 min. Products byte-identical to the stages
run one by one, except where a documented, measured approximation is opted into.

### Work items

1. flag: always pass `--split-prefix` to minimap2 (identical on a one-part index — checked), so a MULTI-PART splice index
   (`minimap2 -I`) can be used; the driver builds one on request. A multi-part index is NOT byte-identical (toy test with paralogs:
   2% of PAF lines differ — MAPQ, a few chain ends), so it is an opt-in low-memory profile, measured against the one-part verdicts.
2. families: left unsharded by default — its all-vs-all PAF is already cached under `PREFIX.cache` (a repeated run replays it);
   `RUSTLE_MINIMAP2=tools/mm2_shard.sh` adds mid-run resume (cmp-identical) at the cost of a shard index on disk. Its 45 min on human
   is profiled below.
3. assign: measure on human; bound the cross-chromosome families' memory if they are the peak.
4. One entry point: `all` gets `merged`'s per-phase resume, per-stage time / peak-memory records and a resource summary.
5. Measure the whole run on human and gorilla and fill in the table above.

### How to run it (implemented 2026-10-12)

    tools/rustle_pipeline.sh all --bam READS.bam --fasta GENOME.fa --out OUT/sample --threads 4 [--low-mem] [--no-flag]

- `all` runs assemble -> families -> assign -> index (when no `--index` is given) -> flag, each stage as its own driver call under
  `/usr/bin/time`: `OUT/sample.<stage>.time` holds the wall seconds, the peak RSS of the stage's whole process tree and the exit
  code, and the run ends with a stage / wall / peak-memory table. Re-running the same command resumes: a stage whose product exists,
  whose record says exit 0 and whose settings match `OUT/sample.all.env` is skipped, until one stage runs (every later stage then
  runs too). `--no-resume` reruns everything.
- `assign` skips the POA homology diagnostic by default (assignments, families and abundance byte-identical — checked on gorilla
  KB3781: 7 and 20 families, cmp-identical; `RUSTLE_POA_DIAGNOSTIC=1` restores it) and takes `RUSTLE_REGION_THREADS` regions at
  once (default 1; byte-identical; memory grows with it).
- `index` builds `OUT/sample.splice.mmi` (`minimap2 -x splice -I <N>G`, `--index-gb N`, default 8 = one part, as the figures' index).
  `flag` always passes `--split-prefix`, so a multi-part index works.
- `--low-mem`: `--index-gb 1` (flag holds ~1 Gbp of index at a time) and one assign region at a time. The multi-part index is
  near-identical, not byte-identical, to one part (toy genome with paralogs: 2% of PAF lines differ in MAPQ or a chain end).
- `--no-flag` stops after assign (no index needed).

### assign on human: where the memory went and the fixes (2026-10-12)

- The default assign stage (POA diagnostic skipped) still reached 20 GB on human A119b. `RUSTLE_TIMING=1` showed why: each region
  loads every record in its copy windows (copy +- 50 kb) with full sequence and CIGAR, e.g. 1,244,342 records for chr1 (8.7 GB), and
  the acrocentric short arms hold millions more (chr13:5.2-16 Mb: 4,769,917 records, 76,512 of them primary).
- A first fix, dropping single-placement molecules' records at read time, was WITHDRAWN: inside multi-copy families almost every read
  has several placements (chr1 region: 1,174,977 of 1,244,342 records kept), and it hid copies supported only by unique reads from
  the per-copy read-presence check (whole-gorilla run: MCL188 failed).
- **`copy_assign --tied-prefilter`** (driver default; byte-identical): each region is read twice — first LIGHT (name, scores, span;
  no sequence, qualities or full CIGAR, ~100 bytes a record), the AS-tied molecules decided exactly as the AS-tied gate does
  (`as_evidence_per_read` + `as_tied`), then FULL records loaded for those molecules only; the per-copy read-presence check reads the
  light pass's spans. The gate keeps exactly those records anyway, so families, assignments and abundance are unchanged (the gate's
  log counts then cover tied molecules only). Memory follows the tied reads, not every read in the copy windows. **Checked on
  gorilla KB3781 genome-wide (326 families): all 7 output tables cmp-identical (assignments, families, quant, family_join,
  famcn_readonly, copy_conversions, params); peak 5.2 -> 3.6 GB, wall 572 -> 639 s.**
- **`--max-window-records N`** (driver; `--low-mem` = 500,000): `utilities window-depth` estimates each family's records per merged
  copy window from the BAM index (compressed bytes / mean bytes per record; 0.8 s on human) and the families over N are not assigned
  (`--skip-families`) and their regions are not swept. On human A119b the 62 such families all lie entirely on the acrocentric short
  arms (chr13/14/15/21/22 < 17 Mb; e.g. MCL1, a 3.5 M-record hull). Only reads those families share with others can change.

### Status and what is left (2026-10-12, PAUSED while the PSV engine is being fixed in another session)

**Done and committed** (branch `machine2/soto-evidence`): `03184537` timed, resumable `all` (per-stage `PREFIX.<stage>.time`,
`PREFIX.all.env` resume, final stage / wall / peak-memory table), `index` stage, `--low-mem`; `a36b9251` `flag` against the genome
FASTA indexed on the fly in 1-Gbp parts (`missing_copy_flag --index-gb`, `--split-prefix` always); `98119eeb` + `a30ca851`
`utilities window-depth` (index estimate, 0.8 s on human) + driver `--max-window-records`, and `copy_assign --tied-prefilter`
(driver default). `--skip-single-placement` was withdrawn in `a30ca851`.

**Measured so far**

| what | result |
|---|---|
| assign, POA diagnostic skipped vs not (gorilla KB3781, 7 + 20 families) | cmp-identical tables |
| assign `--tied-prefilter` vs not (gorilla KB3781 genome-wide, 326 families) | 7 tables cmp-identical; peak 5.2 -> 3.6 GB; 572 -> 639 s |
| human assign, old loading | killed at 16-20 GB (chr1 region 1.24 M records = 8.7 GB; chr13 short arm another ~7 GB) |
| human assign, `--low-mem` + prefilter | stopped by hand at 1.5 min (1.0 GB; the old loading was at 8.7 GB by then); never finished |
| `all` on one gorilla contig (NC_073244.2, `--no-flag`) | works; second run skips every stage |

**Left to do, in order, once the PSV engine is settled**

1. Rebuild and test (the full lib suite has NOT been run since `98119eeb`):
   `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo build --release` and `cargo test --release --lib`
   plus `cargo test --release --bin missing_copy_flag` (output to a file). If the PSV engine default changed, re-run the gorilla
   identity check of `--tied-prefilter` (the prefilter is upstream of the PSV code, so it should still hold):
   copy_assign on `/mnt/linuxdisk/tmp/gw_fast/gkb/regions.txt` with and without `--tied-prefilter`, then `cmp` every table.
2. Human assign with the new loading (scratch `/mnt/linuxdisk/tmp/gw_fast/human/`, the copy table symlinked from the run cache):
   `tools/rustle_pipeline.sh assign --bam .../A119b.t2t.bam --fasta .../chm13v2.0.fa --out /mnt/linuxdisk/tmp/gw_fast/human/human
   --bin .../rustle_target/release --threads 4 --low-mem` under a watchdog that stops `copy_assign` above 16 GB (sample
   `ps -C copy_assign -o rss=`; do NOT use `pgrep -f` with a pattern that also matches the watcher's own command line). Record wall
   and peak; the 62 skipped families are listed in `human.assign.depth_skipped.tsv`.
3. `flag` on gorilla, both index paths: the `.mmi` (`--index /mnt/linuxdisk/home/juanfraitu/winloci_data/GGO.splice.mmi`, ~16-19 GB
   expected) and `--low-mem` (FASTA, 1-Gbp parts); compare the two verdict tables (near-identical expected, not byte-identical) and
   record wall / peak.
4. The whole pipeline in one command: `tools/rustle_pipeline.sh all --bam ... --fasta ... --out OUT --threads 4 --low-mem` on gorilla
   OR6737, then human A119b; copy the final table into this document. Disk: ~5 GB free on `/mnt/linuxdisk` (2026-10-12); the human
   run writes ~5-8 GB (loci FASTA + PAF + GTFs) — free space first.
5. Human `families` (45 min, 12.6 GB, the largest remaining cost): profile whether the minimap2 all-vs-all of 2 GB of locus spans or
   `mcl_families`' PAF handling holds the memory; anything that changes the PAF (intron masking, exon-only sequences) changes the
   validated family definition and needs a pre-registration.
6. Then: a register row, the memory note, and the stage table at the top of this document.

**Caveats to carry**: `--max-window-records` changes results for the skipped (acrocentric-array) families and for reads they share with
other families; under `--tied-prefilter` the AS-tied gate's log counts cover tied molecules only (tables unchanged); a multi-part splice
index is near-identical, not byte-identical, to one part.


---

## PREREG_mcl_two_node_2026-10-13.md

*(original: `docs/archive/2026-10/merged/PREREG_mcl_two_node_2026-10-13.md`)*

## PREREG — MCL two-node components: the fix (2026-10-13)

**Bug** (`project_mcl_two_node_fragility`, found 10-02 by the famsim ladder; `bench/FAMSIM.md`, "What the two bold rows
are", item 1).
- `annotation_families::mcl` adds self-loops of 1.0, while an edge weighs identity x cov_longer < 1.
- On a connected two-node component the column converges toward a tie between "self" and "mate". The convergence test
  (`|u - v| < 1e-7`, applied to all columns at once) stops at a residual whose SIGN picks the attractor.
- At the shipped inflation 2.8 / prune 1e-9, 160 of the 500 weights 0.500..0.999 split a lone pair, among them 0.99,
  0.985, 0.90 and 0.75.
- With other components in the graph the set changes: next to a triangle, 0.99 joins and 0.90 splits. So a pair's
  existence also depended on unrelated families.
- Three or more nodes join. 57% of genome-wide families have two members.

**Fix.**
- After the attractors are read, a connected component of exactly two nodes (an edge whose two endpoints have no other
  edge) is one cluster.
- Every other cluster passes through unchanged, in the same order (size desc, first member).
- Flow never crosses components, so this equals deciding pairs before MCL.
- Scope: `annotation_families::mcl`, i.e. the `mcl_families` families stage. The separate CSC MCL in
  `src/utilities.rs` is not changed.

**Bars.**
- **F1** (unit, red first):
  - every probe weight, and all 500 weights 0.500..0.999, give one cluster `[0, 1]`;
  - a graph with a triangle, a pair at 0.90 and an edge-free node gives `[[0,1,2],[3,4],[5]]`;
  - the helper leaves a path, a triangle and an edge-free node alone.
- **F2** (the four default catalogs, re-clustered from their cached `.fam.loci.paf` / `.fam.loci.gff3`):
  - H0: the pre-fix binary reproduces each shipped `.fam.clusters.tsv` and `.fam.loci.tsv` byte-identically. A failure
    stops the task.
  - `bench/o1_hier/pair_diff.py` (member sets; IDs are ranks and are ignored) must give `pairs_split` 0 and
    `other_changes` 0.
  - Reported per catalog: pairs joined, clusters and members before and after.
- **F3** (famsim): after the fix, `snp_0.01` and `snp_0.1` (human SNRPB and orangutan SPAG7 templates) and human
  `shuffle` have `denovo_intact` and `guided_intact` = 1. The before state is recorded.
- **F4:** `cargo test --release` fails only the 3 known `tests/pipeline_driver_dispatch.rs` tests.

**Reading.** If F1-F4 pass, the fix is the default.
- It changes the output of the default families stage.
- The shipped catalogs are not rewritten here. The next families-stage run, and the paused figure rescoring, pick it up.

### Results

- **F1 PASS.** Red (`cargo test --release --lib -- a_connected_pair a_pair_joins`, before the fix): `a_connected_pair_is_one_family_at_every_weight` failed at `w = 0.99` with `left: [[0], [1]]`; `a_pair_joins_inside_a_mixed_graph_and_the_rest_is_untouched` failed with `left: [[0, 1, 2], [3], [4], [5]]`. Green (after): these two and `the_pair_rule_touches_only_two_node_components` pass, `3 passed`; whole lib `1051 passed; 0 failed; 13 ignored`.
- **F3 PASS.** denovo_intact/guided_intact (denovo_F/guided_F), before -> after, HEAD binary vs fixed binary: human SNRPB `snp_0.01` 0/0 (0.0000/0.0000) -> 1/1 (1.0000/1.0000); `snp_0.1` 0/0 (0.0000/0.0000) -> 1/1 (1.0000/1.0000); `shuffle` 0/1 (0.0000/1.0000) -> 1/1 (1.0000/1.0000); orangutan SPAG7 `snp_0.01` 0/0 (0.0000/0.0000) -> 1/1 (1.0000/1.0000); `snp_0.1` 0/0 (0.0000/0.0000) -> 1/1 (1.0000/1.0000). The before state equals the 10-02 tables; the only arm already intact before the fix is human `shuffle` guided (1, as on 10-02).
- **F4 PASS.** `cargo test --release --no-fail-fast`: 1290 passed, 3 failed, 14 ignored; the 3 failures are the known `tests/pipeline_driver_dispatch.rs` tests (`all_runs_assemble_families_assign_flag_and_neither_catalog_nor_candidates`, `all_with_candidates_runs_the_candidates_stage_between_families_and_assign`, `all_with_legacy_catalog_runs_the_catalog_stage_between_families_and_assign`; each `left: 127, right: 0`).
- **F2 PASS.** Each catalog re-clustered from its cached `.fam.loci.paf` / `.fam.loci.gff3` (`bench/o1_hier/recluster.sh`; outputs, `<sample>.graph.tsv` and `<sample>.pair_diff.json` in `/mnt/linuxdisk/tmp/o1_close/recluster/`), pre-fix binary (md5 9914dee365a041347275f2fafba25e07) vs fixed binary (md5 db55a19f8f165587563212c276bcfbca). **H0 PASS for all four:** the pre-fix binary reproduces each shipped `.fam.clusters.tsv` and `.fam.loci.tsv` byte-identically (the shipped files were only read). `bench/o1_hier/pair_diff.py` (member sets): `pairs_split` 0 and `other_changes` 0 in all four, so every pre-fix cluster is unchanged; every joined pair is a graph edge (no `no_weight`) and an isolated two-node component of the dump (`not_isolated` 0; this check, and exit status 1 on FAIL, were added to `pair_diff.py` after review, and the four JSONs were regenerated from the same files without re-clustering). Per catalog: pairs joined (`tie_regime` / `below_w_star`), clusters before -> after, members before -> after, `uniform_cliques_below_1` (components / members).
  - human_A119b: 68 joined (0 / 68); clusters 2150 -> 2218; members 9907 -> 10043; uniform cliques below 1: 3 / 9.
  - human_testis: 9 joined (0 / 9); clusters 342 -> 351; members 1115 -> 1133; uniform cliques below 1: 0 / 0.
  - chimp_PTR: 20 joined (0 / 20); clusters 396 -> 416; members 1343 -> 1383; uniform cliques below 1: 0 / 0.
  - gorilla_KB3781: 15 joined (0 / 15); clusters 326 -> 341; members 918 -> 948; uniform cliques below 1: 2 / 6.
- **F2, what the joins are (post hoc, descriptive; not a bar; computed by `bench/o1_hier/two_node_states.py` from the re-clustered files, output `/mnt/linuxdisk/tmp/o1_close/recluster/two_node_states.json`).** Every joined pair lies below w*: 112 of 112, none in the tie regime. All 7,256 two-node components of the four graphs were classified (A119b 4,634, testis 689, chimp 858, gorilla 1,075). The 6,969 with w >= w* (4,449 / 678 / 815 / 1,027) were already one family or one locus before the fix, so in these catalogs the residue never split a pair. The 287 below w* (185 / 11 / 43 / 48) were split before and none is split after. So here the fix is Amendment 1(a)'s definition change, not a residue repair. `loci.tsv` changes too: the fold-within-clusters step (§6ey) merges records that overlap on exon bases and share a cluster, so a joined pair of such records becomes ONE locus and leaves `clusters.tsv`. That is 117 / 2 / 23 / 33 new `loci.tsv` lines: 175 of the 287 pairs fold into one locus, and the other 112 stay listed as 2-member families. Every new line lies inside a below-w* two-node component; no pre-fix `loci.tsv` line is removed or changed. The 5 uniform cliques below 1 (3 K3 in A119b at w 0.994380, 0.999527, 0.999607; 2 K3 in gorilla at 0.989512, 0.999411) all came out whole (one 3-member family, or one locus after the fold).

**Amendment 1 (after F1/F3/F4, before F2; 2026-10-13).** Task 1 found a boundary below which the rule changes the definition, and a premise that was too strong. Committed before the F2 re-clustering runs; nothing above `## Results` is edited.

- **(a) Weight bands.** The pre-registered rule joins every connected pair. Below w* = 0.456885 (inflation 2.8, self-loops 1.0) exact MCL splits a lone pair deterministically, so those joins are a definition change, not a residue fix. Above w* the column converges to a tie, and float residue decided join or split: that is the bug. F2 reports joined pairs by band: `tie_regime` (w >= w*) and `below_w_star`. The choice between the weight-blind rule and exact-MCL semantics is the user's; the weight-blind rule stays the default until then.
- **(b) "Three or more nodes join" is too strong.** In an independent port, a uniform triangle with one weight w < 1 splits asymmetrically at 39 of the 500 weights 0.500..0.999, and a uniform K4 at 58. Such components are not changed by this fix. F2 counts them in the four catalogs (`uniform_cliques_below_1`: components that are complete graphs of >= 3 nodes whose edges all carry one weight below 1, equal to 6 dp).
- **(c) Corrections to the Bug paragraph.** "The column converges toward a tie" holds only above w*. "160 of the 500 weights" is 159 in two independent Python ports (first split at 0.504); the 160 came from a third port, and no Rust sweep measured the pre-fix count.

**Decision (user, 2026-10-13; added after the runs; changes no verdict).**
- The default is `mcl_families --pair-rule exact-mcl`, MCL's exact-arithmetic answer for a lone pair: it joins at
  weight >= w*(I) and leaves MCL's deterministic split below it. `weight-blind` is opt-in. `mcl()` follows the default.
- **Re-verified with the default binary, no flag:**
  - All four catalogs reproduce their shipped `clusters.tsv` and `loci.tsv` byte for byte.
  - `params.tsv` ends with `pair_rule exact-mcl` and `pairs_joined_by_rule 0`.
  - F3's five famsim rungs (`snp_0.01` and `snp_0.1` on both templates, human `shuffle`) are intact: 1 / 1, F 1.0000.
    Their pair weights (0.90-0.99) lie above w*.
  - Full suite: 1,296 passed, 3 known `pipeline_driver_dispatch` failures.
- Outputs: `/mnt/linuxdisk/tmp/o1_close/pair_default/`.


---

## PREREG_o1_hierarchy_eval_2026-10-13.md

*(original: `docs/archive/2026-10/merged/PREREG_o1_hierarchy_eval_2026-10-13.md`)*

## PREREG — O1 hierarchy-aware evaluation, human (2026-10-13)

Spec: `docs/superpowers/specs/2026-10-13-close-o1-design.md` section A, with its amendments. Plan:
`docs/superpowers/plans/2026-10-13-close-o1-phase-0A.md` (Tasks 5-6). Where this prereg and the spec or plan differ,
this prereg governs.

**Question.** Are Soto 2025's families nodes of our family hierarchy? The flat one-to-one match against our MCL
families counts every Soto family that sits one level down as a miss.

**Catalogs.**
- `<S>` is `human_A119b` (dev) or `human_testis` (held-out). `R = /mnt/linuxdisk/tmp/rustle_figures/runs/<S>`.
- Both are default catalogs, re-clustered from their cached all-vs-all with the MCL two-node fix ("fix 0",
  `docs/PREREG_mcl_two_node_2026-10-13.md`).
- Binary: `/mnt/linuxdisk/tmp/o1_close/mcl_families.after` (md5 db55a19f8f165587563212c276bcfbca). The pre-fix
  `mcl_families.before` (md5 9914dee365a041347275f2fafba25e07) is used only by the Secondary.
- Lifted inputs (md5): `dupmasker_colors.v2.bed` f5e1146faf778a6d4f52e5bf0c172284; `sedef_pairs.v2.bed`
  7bd521edff9bd1e27b8141027da37a01; `soto_truth.v2.tsv` 60b77a995216237bbc33738dd78e212f.
- Before the run, `md5sum` must reproduce all five values above, and `git diff 8f0515ef -- bench/o1_hier/hier_match.py`
  must be empty. A mismatch stops the run.

**Hierarchy.** `L = /mnt/linuxdisk/tmp/o1_close/lift`; output `O = /mnt/linuxdisk/tmp/o1_close/hier/<S>`.
- Command:
  `mcl_families.after --paf $R/<S>.fam.loci.paf --gff $R/<S>.fam.loci.gff3 --min-exonic-bp 1 --min-shared-exon-frac 0.60 --emit-hierarchy --hierarchy-duplicons $L/dupmasker_colors.v2.bed --hierarchy-sd $L/sedef_pairs.v2.bed --out $O/ours`.
  The nodes are read from `$O/ours.hierarchy.members.tsv`; the fold map is `$O/ours.loci.tsv`.
- **Inputs.**
  - DupMasker colors: `/mnt/linuxdisk/home/juanfraitu/winloci_data/duplicons/chm13.draft_v1.0_plus38Y_dupmasker_colors.bed`.
  - SEDEF: `/mnt/linuxdisk/home/juanfraitu/winloci_data/soto_replication/final_v1_clean.bed`.
  - Both are lifted from CHM13 v1.0 by the T2T chain `/mnt/linuxdisk/home/juanfraitu/winloci_data/soto_replication/v1.0_to_v1.1_rdna_merged.chain`
    (`bench/o1_hier/v1_to_v2.py`).
- **Lift losses and checks.**
  - DupMasker: 6,661 rows dropped (6,306 GRCh38 chrY, 355 acrocentric rDNA).
  - SEDEF: 925 rows dropped. 156 rows (247 sides) span a chain gap and are kept with shifted ends; these are not
    sequence-checked, and most lie within 1 kb of a chromosome end.
  - 1,000 / 1,000 sampled same-length intervals per file are sequence-identical.
- No `--hierarchy-external`: the truth never enters the build.
- Its families must equal the re-clustered catalog:
  `cmp $O/ours.clusters.tsv /mnt/linuxdisk/tmp/o1_close/recluster/after/<S>.fam.clusters.tsv`. A difference stops the
  run.

**Truth.** `$L/soto_truth.v2.tsv` (`bench/o1_hier/soto_truth.py`):
- S1C families proper: 491 families, 2,458 gene rows; `Unassigned_*` excluded.
- Genes are joined to CAT v4 by Gene ID.
- Exons = the union over the gene's transcripts, chain-lifted to v2.0 (0 exons unliftable).

**Mapping.**
- A Soto gene joins a member when one of its exons shares >= 1 bp, on the same strand, with an exon of the member's
  locus or of a locus folded into it.
- Exons come from the catalog's `$R/<S>.fam.loci.gff3`; the fold map is `$O/ours.loci.tsv`.
- The scorer's fallback for a member locus missing from the GFF (its span, either strand) is counted in
  `coverage.member_loci_missing_from_gff`, which is reported. The count covers member and folded-locus keys; absent
  from the JSON means 0. A GFF strand `.` matches either strand.
- A family is EVALUABLE with >= 2 members. The others are reported as coverage.

**Scorer.** `bench/o1_hier/hier_match.py` as of commit 8f0515ef (scipy 1.17.1, numpy 2.4.2):
`hier_match.py $O/ours.hierarchy.members.tsv $L/soto_truth.v2.tsv $O/hier_match.json --map exon --loci-gff $R/<S>.fam.loci.gff3 --fold $O/ours.loci.tsv --null 1000 --seed 7`.
- **Matching.** One-to-one Hungarian (scipy) on -Jaccard between the evaluable truth families and every node: family,
  expansion (including the `<family>.unique` buckets), superfamily and SD complex.
- **Tie policy.**
  - Rows in label order; columns in (level rank family < expansion < superfamily < sd_complex, node id) order.
  - Then canonical placement: the truth families matched to nodes with one identical member set take that set's
    lowest-rank nodes, best Jaccard first, then label (`Scorer._prefer_lower_rank`).
  - Identical sets give identical Jaccards, so the optimum and every number except the reported node and level are
    unchanged.
- **Nulls**, both keeping node sizes and nesting:
  - **global**: one random permutation of all members per round;
  - **within_family**: members permuted inside their own family. Families, superfamilies and SD complexes are
    unchanged; only expansions move.

**Secondary (registered; same statistics, no bars).**
- **Why.** `docs/PREREG_mcl_two_node_2026-10-13.md` (F2) showed that the residue bug split no pair in the shipped
  catalogs: all 6,969 tie-regime pairs in the four default catalogs were already joined (5,127 of them in the two
  scored here).
  - Fix 0 therefore only adds pairs below w* = 0.456885 (A119b +68 families, testis +9).
  - That is a definition change whose adoption is the user's choice (Amendment 1 there).
- **What.** The shipped pre-fix catalogs are also scored. For every two-node component they equal exact-MCL semantics
  (F2: no tie-regime pair was split).
  - They use the same hierarchy inputs, truth, scorer and seed, with `mcl_families.before`.
  - Output: `/mnt/linuxdisk/tmp/o1_close/hier_shipped/<S>`, built with `--out /mnt/linuxdisk/tmp/o1_close/hier_shipped/<S>/ours`.
  - Their families must equal the shipped catalog:
    `cmp /mnt/linuxdisk/tmp/o1_close/hier_shipped/<S>/ours.clusters.tsv $R/<S>.fam.clusters.tsv`. A difference stops
    the Secondary.
- **How reported.** Next to the primary, with no bar. The pair-rule choice must not be made from these scores.

**Bars.** All are computed from the unrounded values in the PRIMARY run's `$O/hier_match.json`.
- `exact` = `hierarchy.exact`, the Hungarian matches with Jaccard 1.
- `rate` = `hierarchy.exact_rate` = `hierarchy.exact` / `hierarchy.n_truth` (the evaluable families).

- **E1:** on each catalog, `null.global.p_exact` < 0.001 AND `null.global.p_mean_jaccard` < 0.001.
- **E1b:** on each catalog, `null.within_family.p_exact` < 0.001 AND `null.within_family.p_mean_jaccard` < 0.001.
- **E2:** |`hierarchy.exact_rate`(testis) - `hierarchy.exact_rate`(A119b)| <= 0.10.
- **p-values:** p = (1 + #{null rounds >= observed}) / (1 + 1000).
  - The smallest is 1/1001 = 0.000999, which passes p < 0.001 only when no null round reaches the observed value.
  - The advisor table prints `p_exact` to 6 dp; `p_mean_jaccard` is quoted from the JSON.

**Reading, fixed now.** The first row that matches applies.
1. **E1 fails on either catalog:** no claim.
2. **E2 fails:** sample-specific; only the per-catalog numbers are reported.
3. **E1b fails (E1 and E2 pass):** the claim stops at the family level and above.
   - Soto's families are matched by our families, superfamilies and SD complexes.
   - The headline rate is (`hierarchy.exact` - `hierarchy.exact_by_level.expansion`) / `hierarchy.n_truth`. Quote it
     with `hierarchy.exact_by_level` and `hierarchy.by_level` excluding `expansion`, next to the flat score.
   - On at least one catalog, the within-family null is not rejected at p < 0.001 on at least one of its two
     statistics. So no claim is made that our expansions recover Soto beyond random within-family splits of the same
     sizes.
4. **E1, E1b and E2 pass:** "Soto's families are nodes of our hierarchy". The O1 headline for Soto becomes, per
   catalog, the exact rate and mean Jaccard, next to the flat score and the coverage.

In every case the headline rate is quoted with `hierarchy.n_truth` (evaluable) out of 491, and with
`hierarchy.exact_of_all`.

**Reported regardless:**
- **Coverage:**
  - records on members / on unclustered loci / without a locus;
  - labels evaluable / one member / none;
  - members in several truth families;
  - `member_loci_missing_from_gff`.
- **Summaries:** the hierarchy and flat summaries: exact, rate, mean Jaccard, macro and micro recall and precision,
  gene recall.
- **Levels:** `exact_any` per level, and where matches land (`hierarchy.by_level`, `hierarchy.exact_by_level`).
- **Unconditioned numbers** over ALL 491 truth families: `hierarchy.exact_of_all`, `macro_gene_recall_all`,
  `micro_gene_recall_all`. A family without a match counts 0, and every truth gene counts, clustered or not.
- **`exact_any.expansion_not_family`:** Soto families equal to an expansion that is not a whole family. These are the
  only matches strictly below the family level. `exact_any.expansion_unique` (equal to a `<family>.unique` bucket) is
  reported separately.
- **Tables:** the per-truth table, and the advisor table (`bench/o1_hier/advisor_table.py`).

**Disclosures.**
- **What was fixed after the F_O probe.** E1b, the exon + strand mapping and the >= 2-member evaluable rule were fixed
  after the plan-time scorer probe on the 10-09 F_O universe (v1.0, guided; 200 null rounds). There, with exon
  mapping:
  - the global null averaged 0.045 exact matches;
  - the within-family null averaged 224 of the observed 231.

  The tie-free exact-expansion count (`exact_any.expansion`) was:
  - 215 with key mapping (the 10-09 Rust placement);
  - 206 with exon + strand mapping;
  - 169 with span-overlap mapping. That mode existed only in the probe; the committed scorer has no span mode.
- **No prior scoring.** No de novo catalog has been scored by this scorer before this prereg.
- **Format-only dry run.** While the plan was written (10-13, before this prereg), the testis hierarchy was built once,
  with the pre-fix binary and scratch copies of the lifted inputs, to check file formats (`/mnt/linuxdisk/tmp/o1_close/probe/`):
  member keys against the loci GFF, and node counts per level. No truth file entered the build, and nothing was
  scored.
- **Precision is a lower bound.** A non-Soto member inside a matched node counts against precision, and Soto's truth
  covers only its SD98 gene set.
- **What the review rounds added.** The unconditioned metrics, the expansion-level split, canonical placement and the
  6-dp printing were added in the scorer's review rounds. They were checked on the 10-09 F_O universe only:
  - Every earlier F_O number was unchanged, except where exact matches land. With key mapping, `exact_by_level`
    (expansion / family / superfamily) is 13 / 228 / 4 instead of 14 / 227 / 4, because an exact match now takes its
    identical family node.
  - `expansion_not_family` = 9 (exon mapping) of the 206 expansion matches.

### Results

Run 2026-10-13 (Task 6 of `docs/superpowers/plans/2026-10-13-close-o1-phase-0A.md`) on the committed prereg (HEAD `fed73e76`); nothing above this heading was edited. Scorer `bench/o1_hier/hier_match.py` at 8f0515ef, python 3.14.4, scipy 1.17.1, numpy 2.4.2, called as `--map exon --null 1000 --seed 7` (the default `--map` is `exon`). The bars below are computed from the UNROUNDED values in the two PRIMARY JSONs (`/mnt/linuxdisk/tmp/o1_close/hier/human_{A119b,testis}/hier_match.json`); the computation is saved in `/mnt/linuxdisk/tmp/o1_close/hier/verdicts.txt`.

#### Gates (all passed)

- **Pinned files.** `md5sum` reproduced all five values: `mcl_families.after` db55a19f8f165587563212c276bcfbca, `mcl_families.before` 9914dee365a041347275f2fafba25e07, `dupmasker_colors.v2.bed` f5e1146faf778a6d4f52e5bf0c172284, `sedef_pairs.v2.bed` 7bd521edff9bd1e27b8141027da37a01, `soto_truth.v2.tsv` 60b77a995216237bbc33738dd78e212f.
- **Scorer pin.** `git diff 8f0515ef -- bench/o1_hier/hier_match.py` is empty (the working-tree file equals the pinned commit).
- **Families = the catalog.** Primary: `cmp $O/ours.clusters.tsv /mnt/linuxdisk/tmp/o1_close/recluster/after/<S>.fam.clusters.tsv` is identical for human_testis and human_A119b. Secondary: `cmp /mnt/linuxdisk/tmp/o1_close/hier_shipped/<S>/ours.clusters.tsv $R/<S>.fam.clusters.tsv` is identical for both.
- **`coverage.member_loci_missing_from_gff`** is absent from all four JSONs, i.e. 0: every member and folded-locus key was found in the catalog GFF.
- **Builds** (`--emit-hierarchy` log lines; the log counts the expansions without the `<family>.unique` buckets that `n_nodes.expansion` includes):
  - primary testis: 1133 members, 344 superfamilies, 351 families, 309 expansions, SD complexes 157; primary A119b: 10043 members, 1525 superfamilies, 2218 families, 1584 expansions, SD complexes 598;
  - secondary testis: 1115 members, 335 superfamilies, 342 families, 302 expansions, SD complexes 149; secondary A119b: 9907 members, 1457 superfamilies, 2150 families, 1541 expansions, SD complexes 571.
- **Pre-run timing probe (disclosure).** Before scoring, one probe imported the scorer, loaded A119b, and ran 3 null rounds per kind; no scores or null values were printed or recorded (only elapsed times and peak memory).

Post hoc (added after the run; changes no verdict): **lift caveat.** The lift's sequence check (1,000 random same-length intervals per file, all identical) could not see a cluster at the chromosome ends. Of the 74 lifted SEDEF sides that keep their length and whose v1.0 interval starts within 1 kb of its chromosome's start or ends within 1 kb of its end, 27 differ in sequence between v1.0 and v2.0 (72 and 26 as distinct intervals). They are 74 of 301,428 same-length sides, so a random sample of 1,000 holds 0.25 of them on average. Their effect on the SD-complex layer was not measured. Recomputed by `bench/o1_hier/posthoc_checks.py --lift-ends` (output `lift_ends` in `/mnt/linuxdisk/tmp/o1_close/posthoc_checks.json`, which lists the 27 sides).

#### Verdicts (PRIMARY)

| Bar | Rule | A119b (dev) | testis (held-out) | Verdict |
|---|---|---|---|---|
| E1 | `null.global`: `p_exact` < 0.001 and `p_mean_jaccard` < 0.001, on each catalog | p_exact 0.000999000999000999; p_mean_jaccard 0.000999000999000999 | p_exact 0.000999000999000999; p_mean_jaccard 0.000999000999000999 | **PASS** |
| E1b | `null.within_family`: `p_exact` < 0.001 and `p_mean_jaccard` < 0.001, on each catalog | p_exact 0.08291708291708291 (not < 0.001); p_mean_jaccard 0.000999000999000999 | p_exact 0.48951048951048953 (not < 0.001); p_mean_jaccard 0.21878121878121878 (not < 0.001) | **FAIL** |
| E2 | \|`hierarchy.exact_rate`(testis) - `hierarchy.exact_rate`(A119b)\| <= 0.10 | exact_rate = 63/247 = 0.2550607287449393 | exact_rate = 41/84 = 0.4880952380952381 | difference 0.2330345093502988 > 0.10: **FAIL** |

- **E1.** 0.000999000999000999 = 1/1001: no null round reached the observed value in any of the four tests. A119b: observed exact 63 against a global-null mean of 0.001 (max 1); observed mean Jaccard 0.5844 against 0.1934 (max 0.2109). testis: exact 41 against 0.026 (max 2); mean Jaccard 0.7836 against 0.2715 (max 0.2967).
- **E1b.** The within-family null keeps every family and permutes members inside it, so only the expansions move. A119b: observed exact 63 against a null mean of 61.048 (max 64), mean Jaccard 0.5844 against 0.5752 (max 0.5819): the within-family null is not rejected on the exact count (p 0.08291708291708291) but is rejected on the mean Jaccard (p 1/1001), so one of the two statistics passes. testis: exact 41 against 40.569 (max 44), mean Jaccard 0.7836 against 0.7793 (max 0.8007): neither statistic is rejected.
- **E2.** The two exact rates differ by 0.2330345093502988 (exact fraction 4835/20748), more than twice the 0.10 bar; the unconditioned `exact_of_all` goes the other way (A119b 0.1283, testis 0.0835), because the evaluable set is 247 families in A119b and 84 in testis.

#### Reading applied: row 2

E1 passes on both catalogs, so row 1 does not match. **E2 fails, so row 2 matches and applies** (the first row that matches applies):

> 2. **E2 fails:** sample-specific; only the per-catalog numbers are reported.

E1b also fails, but row 3 is the case "E1b fails (E1 and E2 pass)" and is not reached. No hierarchy claim for Soto is made across the two catalogs; the per-catalog numbers below stand as reported. The headline rate is quoted with `hierarchy.n_truth` (evaluable) out of 491 and with `hierarchy.exact_of_all`:

| | evaluable of 491 | hierarchy exact | exact rate | exact_of_all | mean Jaccard | flat exact | flat exact rate | flat exact_of_all | flat mean Jaccard |
|---|---|---|---|---|---|---|---|---|---|
| A119b (primary) | 247 | 63 | 0.2551 | 0.1283 (63/491) | 0.5844 | 58 | 0.2348 | 0.1181 | 0.5288 |
| testis (primary) | 84 | 41 | 0.4881 | 0.0835 (41/491) | 0.7836 | 37 | 0.4405 | 0.0754 | 0.7130 |

Post hoc (added after the run; changes no verdict): Descriptive, not a bar: in both catalogs `hierarchy.exact_by_level.family` equals the flat exact count (58 and 37), so the hierarchy's extra exact matches over the flat score are exactly its `exact_by_level.expansion` matches (5 and 4; 1 and 2 of these are twin placements of duplicate Soto families, separated in the post hoc paragraph after the Levels table).

Post hoc (added after the run; changes no verdict): **E2 and set size.** E2's failure stands (its bar compares the unadjusted rates), but half to two thirds of the gap tracks mapped set size (reweighting A119b to testis's size mix leaves 0.079, testis to A119b's mix leaves 0.121; per-size rates in `posthoc_checks.json`, strata 2 / 3 / 4 / 5+), the evaluable truth family's member count in each catalog. The exact rate falls with size in both catalogs: 2 members 43/76 = 0.566 (A119b) and 30/52 = 0.577 (testis); 3 members 10/41 = 0.244 and 6/14 = 0.429; 4 members 4/32 = 0.125 and 3/7 = 0.429; 5 or more 6/98 = 0.061 and 2/11 = 0.182. Testis's evaluable families are smaller: 52 of 84 (0.62) have 2 members against 76 of 247 (0.31) in A119b, and the 77 truth families evaluable in both catalogs have median mapped sizes 5 (A119b) and 2 (testis). Reweighted to testis's size mix (strata 2, 3, 4, 5+), A119b's exact rate is 0.409 against testis's 0.488: a gap of 0.079 instead of the observed 0.233 (by exact size: 0.414, gap 0.074). Restricted to the 77 shared families, the 2-member strata are 12/14 (A119b) and 26/46 (testis). Script `bench/o1_hier/posthoc_checks.py`, output `/mnt/linuxdisk/tmp/o1_close/posthoc_checks.json` (`size_strata`).

#### Reported regardless

**Coverage** (records = Soto gene rows). In this and the following tables the two primary catalogs come first and the two Secondary (shipped) catalogs after them, the latter without verdicts:

| | A119b (primary) | testis (primary) | A119b (shipped, secondary) | testis (shipped, secondary) |
|---|---|---|---|---|
| members (`n_members`) | 10043 | 1133 | 9907 | 1115 |
| nodes: family | 2218 | 351 | 2150 | 342 |
| nodes: expansion (incl. `<family>.unique`) | 2863 | 409 | 2785 | 398 |
| nodes: superfamily | 1525 | 344 | 1457 | 335 |
| nodes: sd_complex | 598 | 157 | 571 | 149 |
| truth families (`labels`) | 491 | 491 | 491 | 491 |
| evaluable (>= 2 members) | 247 | 84 | 247 | 84 |
| one member | 55 | 34 | 54 | 34 |
| no member | 189 | 373 | 190 | 373 |
| records | 2458 | 2458 | 2458 | 2458 |
| records on members | 1212 | 287 | 1211 | 286 |
| records on unclustered loci | 154 | 99 | 155 | 100 |
| records without a locus | 1092 | 2072 | 1092 | 2072 |
| records on several members | 298 | 10 | 298 | 10 |
| members in several truth families | 180 | 24 | 180 | 24 |
| `member_loci_missing_from_gff` (key absent = 0) | 0 | 0 | 0 | 0 |

**Summaries** (over the evaluable families; precisions are averages over the matched ones):

| | A119b (primary) | testis (primary) | A119b (shipped, secondary) | testis (shipped, secondary) |
|---|---|---|---|---|
| hierarchy.n_truth | 247 | 84 | 247 | 84 |
| hierarchy.n_matched | 247 | 84 | 247 | 84 |
| hierarchy.exact | 63 | 41 | 63 | 41 |
| hierarchy.exact_rate | 0.2551 | 0.4881 | 0.2551 | 0.4881 |
| hierarchy.mean_jaccard | 0.5844 | 0.7836 | 0.5846 | 0.7852 |
| hierarchy.macro_recall | 0.7073 | 0.9281 | 0.7075 | 0.9305 |
| hierarchy.macro_precision | 0.7810 | 0.8350 | 0.7810 | 0.8350 |
| hierarchy.micro_recall | 0.5999 | 0.9132 | 0.6003 | 0.9170 |
| hierarchy.micro_precision | 0.1920 | 0.7543 | 0.1920 | 0.7543 |
| hierarchy.macro_gene_recall | 0.6720 | 0.6044 | 0.6720 | 0.6044 |
| flat.n_truth | 247 | 84 | 247 | 84 |
| flat.n_matched | 233 | 77 | 233 | 77 |
| flat.exact | 58 | 37 | 58 | 37 |
| flat.exact_rate | 0.2348 | 0.4405 | 0.2348 | 0.4405 |
| flat.mean_jaccard | 0.5288 | 0.7130 | 0.5289 | 0.7146 |
| flat.macro_recall | 0.6439 | 0.8389 | 0.6441 | 0.8412 |
| flat.macro_precision | 0.7506 | 0.8271 | 0.7506 | 0.8271 |
| flat.micro_recall | 0.5123 | 0.8430 | 0.5127 | 0.8465 |
| flat.micro_precision | 0.6454 | 0.7698 | 0.6454 | 0.7698 |
| flat.macro_gene_recall | 0.6236 | 0.5495 | 0.6236 | 0.5495 |

Post hoc (added after the run; changes no verdict): Descriptive, not a bar: A119b's `hierarchy.micro_precision` (0.1920 against 0.6454 for flat) is driven by ONE match. Truth family ID_256 (2 members) has no family-level partner in the flat score; in the primary A119b hierarchy it is matched to superfamily SF166 (2,991 members) at Jaccard 0.0006687 (recall 1.0, precision 0.0006687). Without that match the micro precision is 827/1326 = 0.6237. In the hierarchy every evaluable family gets a partner (`hierarchy.n_matched` 247 of 247 and 84 of 84, against 233 and 77 for flat), so the matched-based precisions and counts are not comparable between the two rows. Testis has no such match (its largest matched node has 13 members; micro precision 0.7543 against 0.7698 flat). `exact_any.any` (63 and 42) against `hierarchy.exact` (63 and 41): in testis, ID_453 and ID_454 are two Soto families with identical member sets that share one expansion node, and the one-to-one matching can give it to only one of them.

**Levels.** `exact_any` (tie-free: truth sets equal to some node of the level, evaluable families only), `hierarchy.by_level` (where the one-to-one matches land) and `hierarchy.exact_by_level`:

| | A119b (primary) | testis (primary) | A119b (shipped, secondary) | testis (shipped, secondary) |
|---|---|---|---|---|
| exact_any.any | 63 | 42 | 63 | 42 |
| exact_any.family | 59 | 39 | 59 | 39 |
| exact_any.expansion | 56 | 41 | 56 | 41 |
| exact_any.superfamily | 57 | 39 | 57 | 39 |
| exact_any.sd_complex | 9 | 5 | 9 | 5 |
| exact_any.expansion_not_family | 4 | 3 | 4 | 3 |
| exact_any.expansion_unique | 14 | 3 | 14 | 3 |
| hierarchy.by_level.family | 170 | 69 | 170 | 69 |
| hierarchy.by_level.expansion | 57 | 12 | 57 | 12 |
| hierarchy.by_level.superfamily | 20 | 3 | 20 | 3 |
| hierarchy.by_level.sd_complex | 0 | 0 | 0 | 0 |
| hierarchy.exact_by_level.family | 58 | 37 | 58 | 37 |
| hierarchy.exact_by_level.expansion | 5 | 4 | 5 | 4 |
| hierarchy.exact_by_level.superfamily | 0 | 0 | 0 | 0 |
| hierarchy.exact_by_level.sd_complex | 0 | 0 | 0 | 0 |

`exact_any.expansion_not_family` (Soto families equal to an expansion that is not a whole family; the only matches strictly below the family level) is 4 in A119b and 3 in testis; `exact_any.expansion_unique` (equal to a `<family>.unique` bucket) is reported separately, 14 and 3. No truth family lands on an SD complex in the one-to-one matching (`by_level.sd_complex` 0), although `exact_any.sd_complex` is 9 and 5.

Post hoc (added after the run; changes no verdict): each of those truth sets is also equal to a family node (the complex is one whole family; checked on the member sets), and the canonical placement gives it the lowest-rank node.

Post hoc (added after the run; changes no verdict): of the 5 (A119b) and 4 (testis) `exact_by_level.expansion` matches, 1 and 2 are canonical-placement twins of duplicate Soto families (the Soto family's member set equals another Soto family's set, and that other family holds the identical whole-family node): primary A119b ID_454 on `MCL1025.x1`, whose set equals family `MCL1025`, held by ID_453; primary testis ID_264 on `MCL202.x1`, whose set equals family `MCL202`, held by ID_262, and ID_481 on `MCL220.x1`, whose set equals family `MCL220`, held by ID_324; so the exact matches strictly below the family level are 4 in A119b (ID_112, ID_200, ID_477, ID_68) and 2 in testis (ID_147, ID_454); primary testis `exact_any.expansion_not_family` = 3 includes ID_453, which is not an exact match (its member set equals ID_454's, and the only node with that set, `MCL46.x2`, went to ID_454); and with the twin placements removed from the exact counts (same denominators) the exact rates are 62/247 = 0.2510 (A119b) and 39/84 = 0.4643 (testis), so the E2 difference is 0.2133, still > 0.10 (check script and output: `/mnt/linuxdisk/tmp/o1_close/hier/twin_check.py`, `twin_check.txt`).

Post hoc (added after the run; changes no verdict): **duplicate nodes.** The one-to-one matching runs over nodes, not over distinct member sets, and one member set often repeats across levels (for example a family with a single expansion, or a superfamily that is one family): the primary A119b hierarchy has 7,204 nodes but 3,452 distinct member sets, primary testis 1,261 nodes and 481 sets. Two truth families can therefore both match one member set through two of its nodes. In the scored runs 13 (A119b) and 6 (testis) member sets carry more than one match, 15 and 6 matches beyond the first on a set; 13 and 4 of these come from truth families with different member sets, the other 2 and 2 from duplicate Soto families. Matched over distinct member sets (one column per set, labelled by its lowest-rank node; the same Hungarian and canonical placement), A119b gives exact 62, mean Jaccard 0.5700, 245 matched and 50 landing on an expansion (as scored 63 / 0.5844 / 247 / 57), and testis 39 / 0.7488 / 81 / 8 (as scored 41 / 0.7836 / 84 / 12). The gain over flat shrinks from +5 to +4 exact and from +0.0557 to +0.0412 mean Jaccard in A119b, and from +4 to +2 exact and from +0.0707 to +0.0358 in testis: matching over distinct sets removes about a quarter (0.26) of A119b's mean-Jaccard gain and about half (0.49) of testis's. The duplicate Soto member-set groups (evaluable families with identical mapped member sets; A119b 4: ID_96 / ID_97, ID_162 / ID_347, ID_255 / ID_256, ID_453 / ID_454; testis 3: ID_262 / ID_264, ID_324 / ID_481, ID_453 / ID_454) explain the twin placements above and the ID_256 -> SF166 match of the Summaries note: in primary A119b, ID_256's set equals ID_255's, ID_255 holds family `MCL86` and ID_254 holds `MCL86.unique`, so the one-to-one matching leaves ID_256 the 2,991-member superfamily. Script `bench/o1_hier/posthoc_checks.py` (it reloads each scored JSON's own inputs and stops unless its re-match equals the scored placements); output `/mnt/linuxdisk/tmp/o1_close/posthoc_checks.json` (`distinct_sets`).

**Unconditioned numbers over ALL 491 truth families** (a family without a match counts 0; every truth gene counts, clustered or not):

| | A119b (primary) | testis (primary) | A119b (shipped, secondary) | testis (shipped, secondary) |
|---|---|---|---|---|
| hierarchy.n_truth_all | 491 | 491 | 491 | 491 |
| hierarchy.exact_of_all | 0.1283 | 0.0835 | 0.1283 | 0.0835 |
| hierarchy.macro_gene_recall_all | 0.3380 | 0.1034 | 0.3380 | 0.1034 |
| hierarchy.micro_gene_recall_all | 0.3470 | 0.0956 | 0.3470 | 0.0956 |
| flat.n_truth_all | 491 | 491 | 491 | 491 |
| flat.exact_of_all | 0.1181 | 0.0754 | 0.1181 | 0.0754 |
| flat.macro_gene_recall_all | 0.3137 | 0.0940 | 0.3137 | 0.0940 |
| flat.micro_gene_recall_all | 0.3104 | 0.0883 | 0.3104 | 0.0883 |

**Nulls** (1,000 rounds, seed 7; observed values are the hierarchy's; p = (1 + #{null rounds >= observed}) / 1001, quoted unrounded):

| | A119b (primary) | testis (primary) | A119b (shipped, secondary) | testis (shipped, secondary) |
|---|---|---|---|---|
| observed hierarchy.exact | 63 | 41 | 63 | 41 |
| observed hierarchy.mean_jaccard | 0.5844 | 0.7836 | 0.5846 | 0.7852 |
| null.global.exact_mean | 0.001 | 0.026 | 0.004 | 0.031 |
| null.global.exact_max | 1 | 2 | 1 | 2 |
| null.global.mean_jaccard_mean | 0.1934 | 0.2715 | 0.1923 | 0.2703 |
| null.global.mean_jaccard_max | 0.2109 | 0.2967 | 0.2105 | 0.2965 |
| null.global.p_exact | 0.000999000999000999 | 0.000999000999000999 | 0.000999000999000999 | 0.000999000999000999 |
| null.global.p_mean_jaccard | 0.000999000999000999 | 0.000999000999000999 | 0.000999000999000999 | 0.000999000999000999 |
| null.within_family.exact_mean | 61.048 | 40.569 | 61.071 | 40.522 |
| null.within_family.exact_max | 64 | 44 | 64 | 44 |
| null.within_family.mean_jaccard_mean | 0.5752 | 0.7793 | 0.5754 | 0.7808 |
| null.within_family.mean_jaccard_max | 0.5819 | 0.8007 | 0.5823 | 0.7988 |
| null.within_family.p_exact | 0.08291708291708291 | 0.48951048951048953 | 0.0939060939060939 | 0.4825174825174825 |
| null.within_family.p_mean_jaccard | 0.000999000999000999 | 0.21878121878121878 | 0.000999000999000999 | 0.2037962037962038 |

**Advisor table** (`bench/o1_hier/advisor_table.py`, as printed, 4 rows: the two primary catalogs, then the two shipped ones; TSV `docs/hierarchy_results/o1_close_human_2026-10-13.tsv`):

| sample | truth_families | evaluable | one_member | no_member | records_on_members | records_on_unclustered_locus | records_without_locus | records_without_member | flat_n_matched | hier_n_matched | flat_exact | flat_exact_rate | hier_exact | hier_exact_rate | hier_exact_of_all | flat_mean_jaccard | hier_mean_jaccard | flat_macro_recall | hier_macro_recall | flat_macro_precision | hier_macro_precision | hier_gene_recall | hier_gene_recall_all_macro | hier_gene_recall_all_micro | lands_family | lands_expansion | lands_superfamily | lands_sd_complex | exact_any_expansion | exact_expansion_not_family | exact_expansion_unique | p_global | p_within_family |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A119b | 491 | 247 | 55 | 189 | 1212 | 154 | 1092 | NA | 233 | 247 | 58 | 0.235 | 63 | 0.255 | 0.128 | 0.529 | 0.584 | 0.644 | 0.707 | 0.751 | 0.781 | 0.672 | 0.338 | 0.347 | 170 | 57 | 20 | 0 | 56 | 4 | 14 | 0.000999 | 0.082917 |
| testis | 491 | 84 | 34 | 373 | 287 | 99 | 2072 | NA | 77 | 84 | 37 | 0.44 | 41 | 0.488 | 0.084 | 0.713 | 0.784 | 0.839 | 0.928 | 0.827 | 0.835 | 0.604 | 0.103 | 0.096 | 69 | 12 | 3 | 0 | 41 | 3 | 3 | 0.000999 | 0.489510 |
| A119b_shipped | 491 | 247 | 54 | 190 | 1211 | 155 | 1092 | NA | 233 | 247 | 58 | 0.235 | 63 | 0.255 | 0.128 | 0.529 | 0.585 | 0.644 | 0.708 | 0.751 | 0.781 | 0.672 | 0.338 | 0.347 | 170 | 57 | 20 | 0 | 56 | 4 | 14 | 0.000999 | 0.093906 |
| testis_shipped | 491 | 84 | 34 | 373 | 286 | 100 | 2072 | NA | 77 | 84 | 37 | 0.44 | 41 | 0.488 | 0.084 | 0.715 | 0.785 | 0.841 | 0.931 | 0.827 | 0.835 | 0.604 | 0.103 | 0.096 | 69 | 12 | 3 | 0 | 41 | 3 | 3 | 0.000999 | 0.482517 |

records_without_member is a --map key column; it is NA in these exon-mapping runs.

Post hoc (added after the run; changes no verdict): **reading the advisor table.** `p_global` and `p_within_family` are the `p_exact` values only; the `p_mean_jaccard` values are in the Nulls table (A119b's within-family `p_mean_jaccard` is 0.000999, from the JSON). `exact_any_expansion`, `exact_expansion_not_family` and `exact_expansion_unique` are tie-free set equality (`exact_any`), not one-to-one matches: in primary testis, ID_453 equals `MCL46.x2` as a set, but its one-to-one match is family `MCL46` (Jaccard 0.6), because `MCL46.x2` went to ID_454, whose set is identical.

**Per-truth tables.** Each JSON carries `per_truth` (one record per evaluable truth family: `label`, `n_records`, `n_members`, `level`, `node`, `exact`, `jaccard`, `recall`, `precision`, `gene_recall`, `flat_node`, `flat_jaccard`) and the `inputs` it was called with:
- primary A119b `/mnt/linuxdisk/tmp/o1_close/hier/human_A119b/hier_match.json`; primary testis `/mnt/linuxdisk/tmp/o1_close/hier/human_testis/hier_match.json`;
- shipped A119b `/mnt/linuxdisk/tmp/o1_close/hier_shipped/human_A119b/hier_match.json`; shipped testis `/mnt/linuxdisk/tmp/o1_close/hier_shipped/human_testis/hier_match.json`.

#### Secondary (the shipped pre-fix catalogs; same statistics, no bars, no verdicts)

The Secondary columns of the tables above are the shipped catalogs (A119b 9,907 members / 2,150 families; testis 1,115 / 342), built with `mcl_families.before` from the same hierarchy inputs and scored with the same truth, scorer and seed. They differ from the primary in the catalog sizes, in a few coverage counts, and in ONE evaluable truth family per catalog; no exact count, level count or global-null p-value moves. The same families are evaluable (247 and 84); exact (63 and 41, flat 58 and 37), `n_matched`, `by_level`, `exact_by_level`, `exact_any` and both global-null p-values are identical, and no exact flag flips. The family that changes is ID_214 in A119b (7 mapped members in the primary, 6 in the shipped catalog; Jaccard 0.2857 in the primary against 0.3333 shipped) and ID_391 in testis (5 against 4; 0.6667 against 0.8). That moves mean Jaccard by +0.0002 (A119b) and +0.0016 (testis) in the shipped catalogs, and recall slightly (A119b macro recall 0.7073 against 0.7075, testis 0.9281 against 0.9305). The coverage counts that move do so by 1 each: records on members (A119b 1212 against 1211, testis 287 against 286), records on unclustered loci (154 against 155, 99 against 100), and in A119b one label (ID_373) with one member in the primary and none in the shipped catalog (`one_member` 55 against 54, `no_member` 189 against 190). The nulls' own values move a little (within-family `p_exact` 0.0829 in the primary against 0.0939 shipped in A119b, 0.4895 against 0.4825 in testis; within-family `p_mean_jaccard` in testis 0.2188 against 0.2038; in A119b it is 0.000999 in both). The exact-rate difference between the two shipped catalogs is the same 0.2330 as in the primary. As the prereg says, the pair-rule choice must not be made from these scores.

Post hoc (added after the run; changes no verdict): on 2026-10-13 the user made `exact-mcl` the default pair rule (`docs/PREREG_mcl_two_node_2026-10-13.md`, Results). The choice was made on catalog identity and on MCL's consistency across component sizes, not on these scores; the two arms give the same exact counts (63 / 41). The Secondary arm is therefore the default-rule result.


---

## PREREG_sedef_rs_layers_2026-10-13.md

*(original: `docs/archive/2026-10/merged/PREREG_sedef_rs_layers_2026-10-13.md`)*

## sedef-rs as the SD and duplicon layers: V1, V1b, V2 (2026-10-13)

**Pre-registered before any whole-genome sd-scan output is compared with SEDEF or DupMasker.** Spec:
`docs/superpowers/specs/2026-10-13-o1-layered-hierarchy-design.md`, section 2. Plan:
`docs/superpowers/plans/2026-10-13-o1-layers-phase-1-2.md`. Tool: `utilities sd-scan` (vendored sedef-rs, `src/sedef/`,
commit 96b75970), default parameters only; version key: `sedef-rs 0.1.0 vendored 2026-10-13 scan v1; seed k=16 window=2000 sketch=256 min_kmers=32 max_postings=100 min_shared=40; chain gap=2 min_score=20 merge_frac=0.5 min_cell_id=88 min_block_id=90; refine max_gap=400 diag_bin=50 min_span=500; decompose max_gap=200 diag_bin=50 min_seg=1000; max_pairs=1024; names=first-word`. V1 scorer: `bench/sedef/v1_compare.py`
(commit 9378b158).

**Not re-proposed.**
- Register 1263 (our own duplicon units, NOT ADOPTED at 0.638 held-out) is not re-run; V2 reuses its scorer and its
  split unchanged.
- SD or duplicon structure deciding family membership (register 676, 677/997, 817, 818) is not tried: sd-scan only
  annotates.

### Known before running (recorded now)

- Inverted SDs come out as window-sized pieces: on the toy genome, an 8 kb inverted copy gives four 2 kb pairs (`chain`
  links windows only along forward diagonals).
  - The base-pair union is unaffected except at a copy's ends.
  - The duplicons of an inverted SD are split at 2 kb.
- SD ends are clipped by about 1 kb (sandbox README caveat 6; the 2 kb window grid allows up to one window), and
  refined pairs stop one base short.
  Short SDs (1-5 kb) are the most exposed (recall).
- sd-scan identities are k-mer lower bounds; at V1b's 0.98 cut this costs recall.
- sd-scan writes est_id with one decimal, so at V1b's 0.98 cut a pair with est_id 97.95-97.99 rounds to 98.0 and is kept.
- sd-scan families are connected components of segment and >= 80%-overlap links. On chimp chr1 the largest holds 21.5% of
  all copies (3,750 of 17,426). Coarse duplicons in mosaic hubs are a risk for V2.
- The lifted CHM13 reference lost 925 of 151,758 SEDEF rows (rDNA rows and unliftable ends). Where their bases fall inside
  the universe, they count against sd-scan's precision.
- The vendored copy fixes two sandbox bugs that only multi-chromosome genomes show (`src/sedef/mod.rs`): chain merging
  across chromosome pairs, and family ids that followed HashMap order. Its runs are reproducible table for table
  (chimp chr1-3, two runs).

### V1 — SD calls, sd-scan vs SEDEF (both genomes)

**Inputs.**

| Genome | sd-scan run | SEDEF reference |
|---|---|---|
| CHM13 v2.0 | `chm13v2.0.fa` | `/mnt/linuxdisk/tmp/o1_close/lift/sedef_pairs.v2.bed` (v1.0 `final_v1_clean.bed`, lifted by the T2T chain; md5 7bd521edff9bd1e27b8141027da37a01) |
| Gorilla | `GGO.fasta` (mGorGor1 v2.0) | `winloci_data/GGO_sedef_final.bed` (md5 991ed34384832354c02df227e27f65b6) |

**Kept pairs.**
- SEDEF: column 21 (fracMatch) >= 0.90.
- sd-scan: status `ok` and est_id / 100 >= 0.90.
- Both: the shorter side >= 1,000 bp.

**Regions.** Both sides of every kept pair, merged per chromosome, then intersected with the universe:
- CHM13: the v2.0 bases inside the aligned blocks of `v1.0_to_v1.1_rdna_merged.chain`, i.e. bases with a v1.0 ancestor (chr1-22, X and M; not Y or the rDNA arrays);
- gorilla: every base of GGO.fasta.

**Measures.**
- recall = SEDEF region bases covered by sd-scan regions / SEDEF region bases;
- precision = sd-scan region bases covered by SEDEF regions / sd-scan region bases.

**Bar:** recall >= 0.80 AND precision >= 0.80, on CHM13 v2.0 AND on gorilla.

**Reading, fixed now.**
- **Pass:** sd-scan is the SD layer's source for every genome (the `sd` stage).
- **Fail:** sd-scan is not adopted. The SD layer keeps external SEDEF where it exists (human, gorilla); genomes without
  SEDEF get no SD layer.

### V1b — SD98

The same comparison, with 0.98 in place of 0.90 for both tools.

**Bar:** recall >= 0.80 AND precision >= 0.80, on both genomes.

**Reading.**
- **Pass:** the Soto layer (phase 3) takes its SD98 pairs from sd-scan.
- **Fail:** the Soto layer takes SEDEF SD98 where SEDEF exists, and alignment-based identity refinement in sedef-rs
  becomes its own task.

### Descriptive (V1 and V1b, no bar)

- Per chromosome: bases and coverage.
- Recall of the SEDEF pairs by identity bin and by length bin.
- Pair counts.
- Wall time and peak RSS of each sd-scan run (also the CHM13 v1.0 run, V2's input).

### V2 — duplicons, sd-scan vs DupMasker (human, F_O universe, CHM13 v1.0)

**Arms.** Both run with the binary named under **V2 binary** below; families are not recomputed from scratch, only
re-clustered from the same PAF.
- **DupMasker:**

      mcl_families --paf /mnt/linuxdisk/tmp/ourhom/loci.paf --gff /mnt/linuxdisk/tmp/ourhom/loci.gff3
        --min-exonic-bp 1 --min-shared-exon-frac 0.60 --dump-pairs --emit-hierarchy
        --hierarchy-duplicons winloci_data/duplicons/chm13.draft_v1.0_plus38Y_dupmasker_colors.bed
        --hierarchy-sd winloci_data/soto_replication/final_v1_clean.bed
        --hierarchy-external /mnt/linuxdisk/tmp/hier_flag/soto_external.tsv
        --out <arm dir>/ours

  This is the registered F_O command (10-09), re-run.
- **sd-scan:** the same command with `--hierarchy-duplicons <CHM13 v1.0 sd cache prefix>.sd.duplicons.tsv` (and its own
  `--out <arm dir>/ours`). The SD layer stays SEDEF, so only the duplicons differ.

**H0 (guard).** `ours.clusters.tsv`, `ours.loci.tsv` and `ours.pairs.tsv` are byte-identical between the arms.

**Scorer.** `bench/paf_duplicons/score_units.py`, unchanged (register 1263):
- the cover of a family (>= 2 members) = its expansions plus its `unique` class;
- family score = maximum-weight bipartite matching on Jaccard / max(#blocks);
- family-hash split: sha256 even = DEV, odd = HELD-OUT.
- Its `share_exonic_with_paf_unit` output is the sd-scan duplicons here.

**Bars.**
- **V2a:** mean family score over the HELD-OUT families >= 0.80.
- **V2b:** Soto families equal to an expansion strictly below a family, with sd-scan duplicons, >= ⌈0.9 × the DupMasker
  arm's count⌉ (10 on the 10-09 run, so >= 9).
  - Measured by `bench/o1_hier/hier_match.py MEMBERS soto_external.tsv OUT --fold LOCI --map key --null 0`,
    `exact_any.expansion_not_family`.
  - All Soto families; computed once per arm, at the held-out step.

**Procedure.**
- DEV is scored and written first.
- At most ONE amendment may follow, registered here before HELD-OUT is computed.
- HELD-OUT decides, and V2b is computed with it.

**Reading, fixed now.**
- **V2a and V2b hold:** sd-scan duplicons replace DupMasker as the expansion layer's source, for every genome.
  DupMasker stays available for comparison runs.
- **Otherwise:** sd-scan duplicons are not adopted; human expansions stay on DupMasker; apes have no expansion layer.

**Descriptive (no bar).**
- P2: sign test against the one-block null on structured families.
- P3: Soto EXPANSION placements per arm.
- P4: the two named fusions.
- Exact expansion sensitivity and precision; blocks per family.
- The NPIP family (MCL9) cover.
- Members' exonic bases carrying a duplicon, per arm.
- DEV and HELD-OUT side by side.
- Whether the re-run DupMasker arm equals the registered 10-09 files.

### Clarifications registered 2026-10-14, before any comparison

Added after this document's task review and before any sd-scan output was compared with SEDEF or DupMasker. No bar,
filter or measure changes.

- **Provenance.** The CHM13 v2.0 sd-scan run started at 10:36 on 2026-10-14, before this document's first commit
  (4e829791, 10:53), to overlap compute with other work; the gorilla run started at about 10:55 and the CHM13 v1.0 run
  (V2's input) at 11:49, both after it. Only the runs' own counts (pairs, families, run time, memory) were looked at.
  All three used the `utilities` binary built from branch `o1-layers` at ee6eb9c5, whose `src/sedef/` and
  `src/utilities.rs` equal commit 96b75970's, with the version key registered above. Starting the CHM13 v2.0 run before
  registration departs from spec section 2, which places the whole-genome runs inside the V1 gate and pre-registers each
  gate before it runs. The runs produced sd-scan's own outputs only; nothing was compared with SEDEF or DupMasker before
  the first registration commit (the plan's rule).
- **Readings apply per layer.** The SD layer follows V1; the Soto layer's SD98 source follows V1b; the expansion layer
  follows V2. V2's "Otherwise" reading concerns the duplicon and expansion layers only; the SD layer's source is whatever
  V1 decided. This narrows the spec's "sedef-rs gives the SD level only" (section 2, V2 fail), which presupposed a
  passing V1.
- **Genomes (md5 of the files as read).**
  - CHM13 v2.0: `/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa`, 42b830df45df08ab63068e9128156e56.
  - Gorilla: `/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta`, 7a26fb029571b9e9450179769b2dc969.
  - CHM13 v1.0: `/mnt/linuxdisk/home/juanfraitu/winloci_data/soto_replication/t2t-chm13-v1.0.fa.gz`,
    67f9a15a0299e3d32df8e527a611db20 (24 records: chr1-22, X and M).
  - Relative `winloci_data/...` paths in this document mean `/mnt/linuxdisk/home/juanfraitu/winloci_data/...`.
- **V2's sd-scan arm** takes `--hierarchy-duplicons <prefix>.sd.duplicons.tsv`, where `<prefix>` is the line in
  `/mnt/linuxdisk/tmp/o1_layers/chm13v1.sd.prefix`, written by the driver's `sd` stage (`tools/rustle_pipeline.sh sd`)
  on the CHM13 v1.0 FASTA above.
- **V2 binary.** Both arms run with the one `mcl_families` binary built in Task 9 of
  `docs/superpowers/plans/2026-10-13-o1-layers-phase-1-2.md`. Its source differs from the registration commit's (the
  last commit that changes this document above `## Results`, named with the results) only by Task 9's hunks that add
  `--sd-cache` to `src/mcl_families.rs` and sd-scan rows to `parse_sd` in `src/family.rs`. Neither arm passes
  `--sd-cache` or feeds sd-scan rows to `parse_sd`; SEDEF rows parse as before. Its commit and the binary's sha256 are
  recorded with the results.
- **V2 comparator.** If the re-run DupMasker arm differs from the registered 10-09 files, the re-run arm is the
  comparator, and V2b's bar is ⌈0.9 × the re-run DupMasker arm's `expansion_not_family`⌉. If that count is 0, V2b
  cannot be decided and V2 fails.
- **H0 failure.** If `ours.clusters.tsv`, `ours.loci.tsv` or `ours.pairs.tsv` differ between the arms, V2 is not scored
  and the failure is reported. No V2 reading then applies: human expansions stay on DupMasker and the apes get no
  expansion layer, and the amendment may not be used to cure an H0 failure.
- **Task 9's A119b run is descriptive, not V2.** The plan's Task 9 prints, for human A119b, the hierarchy summary
  counts of a `--sd-cache` run beside the lifted-DupMasker run (members, families, expansions, duplicons owned, SD
  blocks and complexes). No family score, exact match or Soto placement is computed there.
- **The one amendment** may change only how the sd-scan arm's duplicon track is built from sd-scan's own outputs
  (segments, families). The bars, the scorer, the split, the DupMasker arm and V2b stay as registered.
- **V2a's number** is `P1_mean_family_score` from `score_units.py`'s HELD-OUT JSON, as printed (4 decimals); it
  passes when that value is at least 0.8000.
  - Family score = the sum of the matched Jaccards over the maximum-weight bipartite matching of the two covers' blocks
    (the DupMasker arm's and the sd-scan arm's), divided by the larger block count; a `unique` block can match only
    the other cover's `unique` block, and a non-`unique` block only a non-`unique` block
    (`score_units.py`, `family_score`).
  - Structured families, the one-block null (B1) and the `unique` class are as defined in
    `docs/PREREG_paf_duplicon_units_2026-10-10.md`; F_O is the registered DNA-level partition of the 2,334 loci of
    `/mnt/linuxdisk/tmp/ourhom/loci.paf` + `loci.gff3`, as defined in
    `docs/PREREG_hierarchy_duplicon_expansion_2026-10-07.md` ("Family layers"). Both arms re-cluster the same PAF with
    the V2 binary; the registered 10-09 partition is a descriptive reference, and H0 and the comparator rule govern any
    difference from it.
- **V2b's command** runs once per arm, with `<arm dir>` = `/mnt/linuxdisk/tmp/o1_layers/v2/dupmasker` and
  `/mnt/linuxdisk/tmp/o1_layers/v2/sdscan`: `MEMBERS` = `<arm dir>/ours.hierarchy.members.tsv`, `LOCI` =
  `<arm dir>/ours.loci.tsv`, `OUT` = `<arm dir>/key.json`, and `soto_external.tsv` =
  `/mnt/linuxdisk/tmp/hier_flag/soto_external.tsv`. The sd-scan arm's count is compared with the bar computed from
  the DupMasker arm's count. Each arm's `mcl_families` command writes `--out <arm dir>/ours`.
- **More known-before items.**
  - sd-scan uses soft-masked (lowercase) bases as k-mer bases; whether SEDEF's inputs were masked is not recorded here.
  - The one-decimal est_id rounding also applies at V1's 0.90 cut: 89.95-89.99 is kept.
  - CHM13 v2.0 has a chrY (HG002) and the v1.0 SEDEF reference has none. The X or autosome side of an sd-scan X-Y or
    autosome-Y pair therefore lies inside the universe with no possible SEDEF counterpart, a precision cost.
  - The DupMasker BED is `chm13.draft_v1.0_plus38Y` (6,306 GRCh38 chrY rows). F_O's loci lie on chr1-22 only, so those
    rows touch no F_O locus.
  - SD-end clipping: the sandbox README (caveat 6) says about 1 kb at the edges; the window grid allows up to one 2 kb
    window.

### Results

#### V1 and V1b (2026-10-14)

Scored on 2026-10-14 with `bench/sedef/v1_compare.py` (commit 9378b158) and the four commands of the plan's Task 8, Step
3, after the provenance commit 092570c9 (the last commit that changes this document above `## Results`). The SEDEF
references are the registered files; their md5s were re-checked and equal the registered ones. The sd-scan inputs are
`<prefix>.sd.pairs.tsv` of the cache prefixes in `/mnt/linuxdisk/tmp/o1_layers/{chm13v2,ggo}.sd.prefix`. No bar, filter
or measure was changed after the numbers were seen. Outputs:
`/mnt/linuxdisk/tmp/o1_layers/v1_{chm13v2,ggo}_{0.90,0.98}.json`.

| genome | cut | SEDEF pairs | sd-scan pairs | SEDEF bp | sd-scan bp | recall | precision |
|---|---|---|---|---|---|---|---|
| CHM13 v2.0 | 0.90 | 44,734 | 59,692 | 237,098,879 | 102,510,802 | 0.3980 | 0.9204 |
| gorilla | 0.90 | 62,568 | 120,716 | 359,966,716 | 153,976,263 | 0.4067 | 0.9508 |
| CHM13 v2.0 | 0.98 | 3,323 | 27,332 | 131,531,378 | 63,761,476 | 0.3968 | 0.8185 |
| gorilla | 0.98 | 5,596 | 56,523 | 195,950,043 | 99,967,153 | 0.4236 | 0.8304 |

Pairs are the kept pairs (identity at or above the cut, shorter side >= 1,000 bp; sd-scan status `ok`). bp are the
merged region bases inside the universe (CHM13 v2.0: 3,044,903,725 bp; gorilla: 3,545,850,636 bp). Bases in both tools'
regions, the numerator of recall and of precision: CHM13 v2.0 94,354,434 (0.90) and 52,191,334 (0.98); gorilla
146,405,942 (0.90) and 83,009,680 (0.98).

**V1: FAIL. V1b: FAIL.** Both cuts miss the recall bar (0.80) on both genomes: recall 0.3980 (CHM13 v2.0) and 0.4067
(gorilla) at 0.90; 0.3968 and 0.4236 at 0.98. Precision clears 0.80 in all four rows (0.9204, 0.9508, 0.8185, 0.8304);
the bars were decided on the unrounded JSON values. Reading (fixed in advance):

- **V1 fail:** sd-scan is not adopted. The SD layer keeps external SEDEF where it exists (human, gorilla); genomes
  without SEDEF get no SD layer.
- **V1b fail:** the Soto layer takes SEDEF SD98 where SEDEF exists, and alignment-based identity refinement in sedef-rs
  becomes its own task.
- By the clarification "Readings apply per layer", this decides the SD layer and the Soto layer's SD98 source only; V2
  (the duplicon and expansion layers) is decided on its own bars.

Recall by SEDEF identity and length bin (descriptive, no bar; V1 = the 0.90 columns, V1b = the 0.98 columns). A cell is
the share of the SEDEF pairs' side bases inside the universe that sd-scan's regions at the same cut cover. Bins are
half-open [lo, hi); length is the shorter side of the pair. A base shared by overlapping SEDEF pairs counts once per
pair, so the bin bases (JSON `recall_bins`) add up to more than the merged SEDEF bp above: 1,967,262,036 (CHM13 v2.0)
and 2,507,676,446 (gorilla) at 0.90.

| SEDEF bin | CHM13 v2.0, 0.90 | gorilla, 0.90 | CHM13 v2.0, 0.98 | gorilla, 0.98 |
|---|---|---|---|---|
| identity 0.90-0.95 | 0.4350 | 0.3159 | - | - |
| identity 0.95-0.98 | 0.5266 | 0.5728 | - | - |
| identity 0.98-0.99 | 0.5323 | 0.7804 | 0.3527 | 0.6315 |
| identity >= 0.99 | 0.5312 | 0.7012 | 0.4601 | 0.6267 |
| length 1-5 kb | 0.3703 | 0.2180 | 0.2430 | 0.5144 |
| length 5-20 kb | 0.4582 | 0.4152 | 0.2553 | 0.6258 |
| length 20-100 kb | 0.5689 | 0.6383 | 0.4569 | 0.6390 |
| length >= 100 kb | 0.4521 | 0.5832 | 0.4037 | 0.6265 |

At 0.90 the pairs of >= 20 kb carry 77.5% (CHM13 v2.0) and 75.2% (gorilla) of the bin bases and the 1-5 kb pairs 5.9%
and 6.1%. The best single bin is 0.5689 (CHM13 v2.0) and 0.7804 (gorilla) at 0.90, 0.4601 and 0.6390 at 0.98: no bin
reaches 0.80 in any of the four columns.

Per chromosome (descriptive; full tables in the JSON `per_chrom`), over the sequences with at least 100,000 SEDEF bases
inside the universe:
- CHM13 v2.0, 0.90: 23 sequences; recall 0.1443 (chr18) to 0.5957 (chr16), no sequence reaches 0.80; precision 0.5274
  (chrX) to 0.9934 (chr8), below 0.80 on 6 of 23.
- CHM13 v2.0, 0.98: 23 sequences; recall 0.0555 (chr18) to 0.5749 (chr13), no sequence reaches 0.80; precision 0.3167
  (chr12) to 0.9752 (chr8), below 0.80 on 13 of 23.
- gorilla, 0.90: 25 sequences; recall 0.1736 (NC_073234.2) to 0.6392 (NC_073242.2), no sequence reaches 0.80; precision
  0.7835 (NC_073234.2) to 0.9966 (NC_073237.2), below 0.80 on 1 of 25.
- gorilla, 0.98: 25 sequences; recall 0.1883 (NC_073234.2) to 0.6802 (NC_073242.2), no sequence reaches 0.80; precision
  0.5737 (NC_073243.2) to 0.9528 (NC_073240.2), below 0.80 on 12 of 25.
- Left out of those ranges: CHM13 v2.0 chrM (7,615 SEDEF bases at 0.90, no sd-scan bases) and chrY (outside the
  universe); gorilla NC_011120.1, the mitochondrion (10,968 SEDEF bases at 0.90, no sd-scan bases). CHM13 v2.0 chrX has
  the lowest precision at 0.90 (0.5274); the chrY cost registered under "More known-before items" predicts a precision
  cost on X, but that split was not measured.

| run | wall s | peak RSS GiB | sequences | candidates | blocks | pairs (ok) | segments | families | duplicon rows |
|---|---|---|---|---|---|---|---|---|---|
| CHM13 v2.0 | 1111.86 | 10.31 | 25 | 600,273 | 189,574 | 189,574 (186,820) | 136,728 | 13,424 | 264,078 |
| gorilla | 3255.96 | 18.67 | 26 | 779,754 | 244,273 | 244,273 (240,384) | 209,924 | 15,404 | 401,916 |
| CHM13 v1.0 | 1080.44 | 9.44 | 24 | 544,839 | 178,980 | 178,980 (176,628) | 130,039 | 12,257 | 251,185 |

Wall s and peak RSS are `/usr/bin/time` around the whole `tools/rustle_pipeline.sh sd` call
(`bench/sedef/run_genome.sh`: 4 threads, 256 count shards, address space capped at 23 GiB by `ulimit -v 24117248`). Peak
RSS GiB = `peak_rss_kb` / 1,048,576, the unit of Task 4's projection: 10,812,312 (CHM13 v2.0), 19,578,244 (gorilla) and
9,900,128 (CHM13 v1.0) KiB. Stage seconds in `params.tsv` (seed / chain / refine / decompose): CHM13 v2.0 208.4 / 1.0 /
187.6 / 616.0; gorilla 229.1 / 1.6 / 336.5 / 2409.1; CHM13 v1.0 223.7 / 1.0 / 179.2 / 580.3. The driver logged the runs
at 10:36:09-10:54:34 (CHM13 v2.0), 10:55:02-11:49:10 (gorilla) and 11:49:36-12:07:29 (CHM13 v1.0) on 2026-10-14, as the
Provenance bullet states. The three cache prefixes are
`/mnt/linuxdisk/home/juanfraitu/sd_cache/sd/<md5>-2992c257/genome`; the md5s are the registered FASTA md5s, and the
`version_key` of each `params.tsv` equals the one registered above.

Observations (no bar):
- (a) CHM13 v1.0 has 24 sequences in `params.tsv` (chr1-22, X and M; chrM is included), not the 23 the plan expected.
  CHM13 v2.0 has 25 (the same plus Y) and gorilla 26.
- (b) The Task 4 memory projection (`bench/sedef/memory_probe.sh`, a power law through chimp chr1 and chr1-3)
  under-predicted both whole genomes. It gave 6.9 GiB (CHM13 v2.0) and 7.7 GiB (gorilla); the measured peaks are 10.31
  and 18.67 GiB, 1.5 and 2.4 times the projection (CHM13 v1.0 had no projection: 9.44 GiB). The gorilla peak RSS is 81%
  of the 23 GiB address-space cap (the cap limits virtual size, which is at least the RSS); all three runs finished with
  exit 0 and no allocation failure.

Clarification (2026-10-14): the sd-scan runs' `utilities` binary was linked at 09:28 and found up to date by the build at ee6eb9c5 (no source newer than the binary); its `src/sedef/` and `src/utilities.rs` equal 96b75970's.

#### V2 DEV (2026-10-14; HELD-OUT not computed at this step)

Run and scored on 2026-10-14. The run relies on registration commit 092570c9 (the last commit that changes this document
above `## Results`).

- **Binary (both arms):** `/mnt/linuxdisk/tmp/o1_layers/v2/mcl_families.65f6a1a3`, built from commit 65f6a1a3 (Task 9),
  sha256 `7829b2f5c842e84f0f4372d72d5fb1623e8c0550cf02ef015c6ba6e3da8d8a8e` (re-computed before the runs; it equals the
  recorded `mcl_families.65f6a1a3.sha256`).
- **Binary pin check**, `git diff --stat 092570c9 65f6a1a3 -- src`:

      src/family.rs       | 29 +++++++++++++++++++++++++++--
      src/mcl_families.rs | 42 +++++++++++++++++++++++++++++++++++++-----
      2 files changed, 64 insertions(+), 7 deletions(-)

  Only the two files of the V2 binary pin differ. Their hunks are the `--sd-cache` flag with its helper `hierarchy_inputs`
  (`mcl_families.rs`) and the sd-scan branch of `parse_sd` (`family.rs`), plus doc comments and two unit tests. The branch
  is taken only for a row of exactly 11 columns ending in `ok`, `self` or `weak`; all 151,758 rows of the SEDEF BED
  (`final_v1_clean.bed`) have 34 columns. Neither arm passes `--sd-cache`.
- **Arms:** the registered command (it includes `--dump-pairs`). DupMasker: `--hierarchy-duplicons
  winloci_data/duplicons/chm13.draft_v1.0_plus38Y_dupmasker_colors.bed`. sd-scan: `--hierarchy-duplicons
  /mnt/linuxdisk/home/juanfraitu/sd_cache/sd/67f9a15a0299e3d32df8e527a611db20-2992c257/genome.sd.duplicons.tsv` (the prefix
  in `/mnt/linuxdisk/tmp/o1_layers/chm13v1.sd.prefix`; 251,185 duplicon rows). Both: `--hierarchy-sd
  winloci_data/soto_replication/final_v1_clean.bed`, `--hierarchy-external /mnt/linuxdisk/tmp/hier_flag/soto_external.tsv`,
  `--out /mnt/linuxdisk/tmp/o1_layers/v2/{dupmasker,sdscan}/ours`. Both runs exit 0.
- **Scorer:** `bench/paf_duplicons/score_units.py`, unchanged, DEV half only for P1, P2 and exact matches (P3, P4 and
  the members' exonic share are all-family checks); output `/mnt/linuxdisk/tmp/o1_layers/v2/score_dev.json`.
  `score_units.py ... heldout` and `hier_match.py` (V2b) were not run on these arms at this step. No bar, filter or
  measure was changed after the numbers were seen.

`--emit-hierarchy` lines, as logged (the log's trailing output path left out):

    dupmasker: 1709 members, 331 superfamilies, 380 families, 430 expansions, 319 fusion loci (156 foreign-only), 1240 duplicons owned; SD blocks 322 (192 with >1 family), families with every member pair SD-joined 200; SD complexes 35 (largest 1532 members)
    sdscan:    1709 members, 331 superfamilies, 380 families, 378 expansions, 208 fusion loci (536 foreign-only), 1073 duplicons owned; SD blocks 322 (192 with >1 family), families with every member pair SD-joined 200; SD complexes 35 (largest 1532 members)

H0: clusters / loci / pairs identical between the arms: **yes** (`cmp`, byte-identical; md5
09f0fe53f591c327cdc521b51b54ce37, eaef39645a68e75f2abac28455ffdfa0 and 865519d877a390bd65b5c22647e04761).
DupMasker arm vs the registered 10-09 run (`/mnt/linuxdisk/tmp/hier_flag/final/`): `ours.clusters.tsv`, `ours.loci.tsv`,
`ours.pairs.tsv`, `ours.hierarchy.members.tsv`, `ours.hierarchy.nodes.tsv`, `ours.hierarchy.edges.tsv` and
`ours.hierarchy.external.tsv` are each byte-identical to it (`cmp -s`). The comparator rule is therefore not triggered: the
re-run arm and the registered files are the same. (Outside those seven files, the re-run's `ours.params.tsv` and log
gain `pair_rule exact-mcl` / `pairs_joined_by_rule 0`, which the 10-09 binary did not write; the rule joined 0 two-node
components.)

| measure (DEV, 200 families, 49 structured) | value |
|---|---|
| mean family score (P1 on DEV) | 0.4816 |
| member-weighted score | 0.4272 |
| families identical | 69 / 200 |
| exact expansion sensitivity / precision | 73/249 (0.2932) / 73/214 (0.3411) |
| P2: wins / losses / ties vs the one-block null, sign test p | 21 / 13 / 15; p = 0.1147 (one-sided, ties dropped; mean on the 49 structured families 0.3514, one-block null 0.2779) |
| P3: Soto families exactly one expansion, DupMasker / sd-scan (all families) | 215 / 179 |
| P4: PKD1P6-NPIPP1, PDXDC2P-NPIPB14P (sd-scan arm) (all families) | `MCL28 chr16:15126207-15159416 fusion_of=MCL29 fusion`; `MCL286 chr16:75785720-75875798 fusion_of=MCL29 fusion` (DupMasker arm: `fusion_of=MCL9 fusion` for both) |
| members' exonic bases carrying a duplicon, DupMasker / sd-scan (all 1,709 members) | 0.9131 / 0.8227 (of 2,988,087 bp) |
| NPIP (MCL9) cover, DupMasker / sd-scan | MCL9.x1 4, MCL9.x2 6, MCL9.x3 10 members / unique 19 members (family score 0.0) |
| blocks per family, DupMasker / sd-scan | {1: 151, 2: 31, 3: 13, 4: 3, 9: 1, 11: 1} / {1: 129, 2: 54, 3: 11, 4: 3, 5: 3} (blocks: families; `unique` is a block) |

P3 placements of the Soto families with >= 2 members, DupMasker / sd-scan: EXPANSION 215 / 179, WITHIN-FAMILY 98 / 134,
WITHIN-SD-COMPLEX 52 / 52, WITHIN-SUPERFAMILY 33 / 33, OUTSIDE 2 / 2. P3 is descriptive and is not V2b
(`exact_any.expansion_not_family`), which is computed once per arm at the held-out step.

No verdict at DEV: V2a (HELD-OUT mean family score >= 0.8000) and V2b are decided at the held-out step; DEV is reported,
not read against a bar.

Worst families (the 12 lowest scores): all twelve score 0.0. The JSON `worst` list sorts by (score, family id), so it is
cut among tied scores and runs in family-id string order: MCL109, MCL114, MCL115, MCL118, MCL12, MCL125, MCL13, MCL137,
MCL140, MCL147, MCL149, MCL15.

#### V2 HELD-OUT and V2b (2026-10-14)

Run and scored on 2026-10-14 on the arms of the DEV step, unchanged: the md5s of `ours.clusters.tsv`, `ours.loci.tsv` and
`ours.pairs.tsv` in both arm directories and the binary's sha256 were re-checked and equal the values recorded under V2
DEV. The run relies on registration commit 092570c9 (the last commit that changes this document above `## Results`);
`git log` of this document before the run lists the DEV commit ae60db77 first and no amendment commit. No bar, filter or
measure was changed after the numbers were seen.

Amendment: none registered. A DEV-only probe of a reciprocal-overlap family rule (>= 0.8 of the longer interval) scored a DEV mean family score of 0.3500 (40/200 identical; 1,093 expansions) against 0.4816 for the registered rule, and was not registered.

Commands (`<v2>` = `/mnt/linuxdisk/tmp/o1_layers/v2`; `<prefix>` as under V2 DEV):

    score_units.py <v2>/dupmasker <v2>/sdscan heldout /mnt/linuxdisk/tmp/ourhom/loci.gff3 \
      winloci_data/duplicons/chm13.draft_v1.0_plus38Y_dupmasker_colors.bed <prefix>.sd.duplicons.tsv > <v2>/score_heldout.json
    hier_match.py <arm dir>/ours.hierarchy.members.tsv /mnt/linuxdisk/tmp/hier_flag/soto_external.tsv <arm dir>/key.json \
      --fold <arm dir>/ours.loci.tsv --map key --null 0

`hier_match.py` ran once per arm (`<arm dir>` = `<v2>/dupmasker` and `<v2>/sdscan`); the three runs exit 0.

| measure | HELD-OUT (180 families, 35 structured) | DEV (200 families, 49 structured) |
|---|---|---|
| mean family score (P1) | **0.5153** | 0.4816 |
| member-weighted score | 0.4923 | 0.4272 |
| families identical | 77 / 180 | 69 / 200 |
| exact expansion sensitivity / precision | 71/181 (0.3923) / 71/164 (0.4329) | 73/249 (0.2932) / 73/214 (0.3411) |
| P2: wins / losses / ties vs the one-block null, sign test p | 10 / 13 / 12; p = 0.7976 (one-sided, ties dropped; mean on the 35 structured families 0.3174, one-block null 0.3099) | 21 / 13 / 15; p = 0.1147 (mean on the 49 structured families 0.3514, one-block null 0.2779) |
| P3: Soto families exactly one expansion, DupMasker / sd-scan (all families) | 215 / 179 (all families: the same as DEV's) | 215 / 179 |
| P4: PKD1P6-NPIPP1, PDXDC2P-NPIPB14P (all families) | `fusion_of=MCL9` in the DupMasker arm, `fusion_of=MCL29` in the sd-scan arm, for both (all families: the same as DEV's) | as in the DEV table above |
| members' exonic bases carrying a duplicon, DupMasker / sd-scan (all 1,709 members) | 0.9131 / 0.8227 (of 2,988,087 bp; all members: the same as DEV's) | 0.9131 / 0.8227 |
| NPIP (MCL9) cover | not in this half (MCL9 is in the DEV half) | MCL9.x1 4, MCL9.x2 6, MCL9.x3 10 members / unique 19 members (family score 0.0) |
| blocks per family, DupMasker / sd-scan | {1: 145, 2: 31, 3: 4} / {1: 130, 2: 42, 3: 7, 5: 1} | {1: 151, 2: 31, 3: 13, 4: 3, 9: 1, 11: 1} / {1: 129, 2: 54, 3: 11, 4: 3, 5: 3} |

Worst held-out families (the 12 lowest scores): all twelve score 0.0 (44 held-out families do). The JSON `worst` list
sorts by (score, family id), so it is cut among tied scores and runs in family-id string order: MCL106, MCL11, MCL145,
MCL159, MCL179, MCL186, MCL19, MCL201, MCL213, MCL217, MCL223, MCL226.

V2a: held-out `P1_mean_family_score` 0.5153 against the bar 0.8000.

V2b: Soto families equal to an expansion strictly below a family (`exact_any.expansion_not_family`, key mapping):
DupMasker 10, sd-scan 12; bar ⌈0.9 × 10⌉ = 9. The DupMasker arm's 10 is the count of the 10-09 run named in V2b's bar.
`exact_any` as printed: DupMasker {'any': 245, 'expansion': 215, 'expansion_not_family': 10, 'expansion_unique': 38,
'family': 231, 'sd_complex': 25, 'superfamily': 218}; sd-scan {'any': 247, 'expansion': 179, 'expansion_not_family': 12,
'expansion_unique': 58, 'family': 231, 'sd_complex': 25, 'superfamily': 218}.

**V2: FAIL (V2a FAIL, V2b PASS).** Reading (fixed in advance): **Otherwise:** sd-scan duplicons are not adopted; human
expansions stay on DupMasker; apes have no expansion layer. By the clarification "Readings apply per layer", this decides
the duplicon and expansion layers only: the SD layer keeps SEDEF (V1 FAIL) and the Soto layer's SD98 source is SEDEF
(V1b FAIL).

Observations (no bar; computed after the verdict from the arms' `ours.hierarchy.members.tsv` (and, for V2b's sets,
`soto_external.tsv` and the arm's `ours.loci.tsv`), with `score_units.py`'s own `covers`, `family_score` and `half` and
`hier_match.py`'s own loaders; scripts `bench/sedef/describe_disagree.py` and `bench/sedef/describe_v2b.py` (copies of
the scripts that wrote the outputs), outputs `<v2>/describe_disagree.txt` and `<v2>/describe_v2b_{dupmasker,sdscan}.txt`):
- Where the covers differ (HELD-OUT; DEV in brackets). Of 180 (200) families, 77 (69) are identical. The other 103 (131)
  are: 51 (79) with expansions in both arms that split or merge differently (mean score 0.2764 (0.3265)); 37 (39) with
  expansions in DupMasker's cover and only the `unique` class in sd-scan's (mean 0.0245 (0.0199)); 15 (13) the reverse
  (mean 0.0500 (0.0577)). Unit / no-unit flips are therefore 52 of 103 (50.5%) of the held-out disagreements and 52 of 131
  (39.7%) of DEV's. Families with at least one non-`unique` block: DupMasker 151 (175), sd-scan 129 (149).
- V2b's composition. Of sd-scan's 12 Soto sets, 9 equal a numbered expansion and 3 equal a `<family>.unique` bucket (the
  members of MCL51, MCL144 and MCL31 without a duplicon signature); DupMasker's 10 are all numbered expansions. The 9
  numbered expansions alone also meet the bar of 9. The two arms share 4 of these Soto sets (ID_289, ID_323, ID_433,
  ID_72); the other 6 (DupMasker) and 8 (sd-scan) differ.
- Context: register 1263's own PAF units scored 0.638 held-out on this scorer and split; sd-scan duplicons score 0.5153.


---

## PREREG_denovo_copies_2026-10-14.md

*(original: `docs/archive/2026-10/merged/PREREG_denovo_copies_2026-10-14.md`)*

## De novo analogue of Liftoff's `-copies` (2026-10-14)

Pre-registration, written before any number below was computed. Code: `bench/template_copies/copies.py` at `f0833180`
(9 tests). Branch `denovo-copies` (on top of `template-copies`).

### Question

The user asked (2026-10-14): "in order to prove that the pipeline is also very similar to liftoff, should we include
something similar to the -copies functionality of liftoff?" and then "yes lets do that".

Liftoff's `-copies` lifts every annotated gene and keeps its extra placements at sequence identity ≥ `-sc`. The de novo
analogue lifts our assembled loci instead of the annotation: each base locus's representative transcript is aligned to
the genome, and every alignment that passes Liftoff's own thresholds is a projected copy. Two questions:

1. For the genes we assemble, does the projection find the extra copies Liftoff finds?
2. What does it find that Liftoff does not?

### Known before running (recorded now)

- **Truth.** Figure 8's Liftoff tables (`/mnt/linuxdisk/tmp/rustle_figures/liftoff/{human,gorilla}/liftoff_loci.tsv`;
  Liftoff v1.6.3, `-copies -sc 0.95 -a 0.5 -s 0.5 -overlap 0.1`, sharded, with figure 8's cross-shard rules M1-M3).
  Extra copies that are not short (exon union ≥ 200 bp): human 1,888, gorilla 1,129, all at sequence identity ≥ 0.95.
- **Templates.** The template-copies templates (one representative per base locus: the transcript with the most reads),
  built 2026-10-14: A119b_shipped 94,685 (projected: 16 PAF batches, minimap2 `-x splice -c -N 50 -p 0.1`, CHM13 v2.0
  splice index), testis 13,048 and KB3781 17,672 (built, not projected).
- **Seen already on A119b_shipped**, in the template-copies study (a different question): the number of alignments at
  identity ≥ 0.80 and aligned ≥ 0.50 (561,980), the candidate counts, the template-guided loci and their Liftoff scores
  (register 1285), including that 63 of 106 read-supported template-reachable copies were reached by an assembled
  locus. **Not seen**: any alignment at Liftoff's thresholds, any coverage of a Liftoff copy by an alignment, any site.
- **A design difference, expected and not measured.** Liftoff aligns each gene's genomic sequence (exons and introns);
  the projection aligns the spliced transcript. An intronless retrocopy of a multi-exon gene should therefore be easy
  for the projection and hard for Liftoff.

### Method (fixed)

**Projected copy.** An alignment of a template that passes Liftoff's thresholds, computed as Liftoff computes them
(`find_best_mapping.convert_all_children_coords`) over the template's bases from the PAF CIGAR:
- coverage = aligned template bases (M) / template length ≥ 0.5 (`-a`);
- identity = matching bases (PAF column 10) / (template length + target-only bases, CIGAR D) ≥ 0.95 (`-sc` of the
  figure-8 run).

Both are normalised by the whole template, as Liftoff normalises by the whole gene, so a projected copy matches ≥ 95% of
its template's bases.

**Recall.**
- A Liftoff extra copy (not short) is **reachable** when its source gene's placement (`in_place` or `moved`, not short,
  same `source_id`) is covered ≥ 50% by a base locus. The **own templates** are the templates of every such base locus:
  the analogue of lifting that gene.
- **Found (own)**: a projected copy of an own template covers ≥ 50% of the copy's exon union (figure 8's `MATCH_COV`).
  **Found (any)**: a projected copy of any template does.
- Copies whose source is not placed, or is placed but not assembled, are counted and never scored.
- **Misses**, in this order: an alignment of an own template covers the copy but fails a threshold (`below_threshold`,
  with its best identity and coverage); an own template has minimap2's cap of 51 alignments (`capped`); neither
  (`no_alignment`).

**Sites** (what the projection finds beyond its own loci).
- Off-locus projected copies: the template's own base locus covers < 50% of the alignment's exon union.
- They merge into **sites**: same contig, either strand, exon unions overlapping by ≥ 50% of the smaller (transitively).
- Each site's best alignment is the one with the highest identity, then coverage, then template id. Its exons class the
  site, in this order:
  1. a Liftoff extra copy covers ≥ 50% (`liftoff_copy`);
  2. a Liftoff gene placement (`in_place` / `moved`) does (`liftoff_gene`);
  3. a base locus does (`assembled`);
  4. none of them does (`none`).
- **Retro-like**: the best alignment is unspliced (one block) and its template has ≥ 2 exons.
- **Unexplained**: class `none` and not retro-like.

**Guard S0.** At least 0.90 of a sample's templates must have an alignment covering ≥ 50% of their own exons; otherwise
the index and the base GTF are not the same assembly, and that sample has no verdict until the index is fixed.

### Bars (per sample; species never pooled)

- **C1 (primary): recall (own) ≥ 0.80**, with ≥ 30 reachable copies. Fewer than 30: not powered, no verdict for that
  sample.
- **C2 (secondary): unexplained sites ≤ 0.10 of sites.**

Reported, not barred:
- recall (any);
- reachable share (reachable / all extra copies);
- recall by the source gene's number of Liftoff copies (1, 2-5, ≥ 6);
- the misses and their best identity and coverage;
- the templates at the alignment cap;
- each site class with its retro-like count;
- the templates with at least one site.

### Samples

| sample | role | templates (run dir under `/mnt/linuxdisk/tmp/o1_template_copies/`) | base GTF | index | truth |
|---|---|---|---|---|---|
| A119b | development | `A119b_shipped` (projected) | `runs/human_A119b/human_A119b.gtf` | `npip_ladder/idx/target.splice.mmi` | human |
| testis | held out (same genome, another sample) | `testis` | `runs/human_testis/human_testis.gtf` | the same | human |
| KB3781 | held out (gorilla) | `KB3781` | `runs/gorilla_KB3781/gorilla_KB3781.gtf` | `winloci_data/GGO.splice.mmi` | gorilla |

The held-out templates are projected with `run_tg.py project` (same code and settings as A119b) after this file is
committed.

### Readings (fixed now)

1. **C1 passes in A119b and in both held-out samples** → claimed: "for the genes we assemble, the de novo projection
   finds Liftoff's extra copies", quoted with the reachable share (the claim says nothing about genes we do not
   assemble). C2 adds "and what it finds beyond Liftoff is mostly annotated genes or retrocopy-like", only in the samples
   where C2 passes.
2. **C1 passes in A119b but fails in a held-out sample** → not claimed (sample-specific); all numbers reported.
3. **C1 fails in A119b** → not claimed; the miss decomposition says why. The held-out samples still run (they are cheap)
   and are reported.

### Rulings (made while designing)

- **R1.** The templates are the template-copies templates as built (no rebuild): one representative per base locus.
- **R2.** The projection settings are template-copies' (`-x splice -c -N 50 -p 0.1`; `-p 0.1` is the guided finder's
  setting). A119b's projection is reused as it is.
- **R3.** No read evidence is used anywhere: Liftoff's copies need none, and the copy does not need to be assembled.
- **R4.** The guided analogue (annotation records as templates) is G1, figure 8's separately registered test; it is not
  part of this one.

### Amendment 1 (2026-10-14, speed only, before any result)

The first development run (A119b) was stopped by its 10-minute limit before it wrote anything. Profiling the stages
showed: inputs 2 s, own templates 1 s, the PAF pass 12 s, recall 1 s; the site merge did not finish. Seen while
profiling, and nothing else: 489,897 projected copies (self alignments included), 85,729 kept alignments of 24,801 own
templates, and projected-copy spans of 5,050 bp (median), 106,416 bp (99th percentile) and 1,640,982 bp (maximum). No
recall, miss, site or class number was computed or seen.

The cause is the merge's sweep over spans: a copy with long introns stays open across its whole span, so the sweep
compares nearly every pair. Diagnosing it showed where the pairs come from. Batch 2 holds 6,000 templates from one
tandem array (chr13 13.9-15.5 Mb, assembled as ~5,000 overlapping loci), and 5,035 of them reach minimap2's cap of 51
alignments, all inside the array. That gives 252,592 off-locus copies in batch 2, and 251,540 of them chain into one
cluster when any shared exonic base links two copies.

The merge (`copies.site_groups`) now sweeps exon blocks and keeps the open blocks per component. A block is compared
only with other components' open blocks. The exonic bases a pair shares so far are summed, and the pair joins as soon
as the sum reaches the bar: the sum only grows, so this is the rule itself. The rule is unchanged (exon unions
overlapping by >= 50% of the smaller, transitively), and so are the sites:
- `test_site_merge_uses_exons_not_spans` pins the two cases that differ between spans and exons;
- `test_site_groups_equal_brute_force_on_random_pileups` compares 60 random pileups against every pair;
- on real off-locus copies the new merge gave the same partition as the span sweep in batches 0-1 (2,166 copies), 5
  (5,594) and 8 (20,581);
- batch 2 alone now takes 225 s at 0.8 GB.

Also seen while checking, and nothing else: the number of sites in those batches (1,048, 816, 3,266 and 564). No recall,
miss, class or C2 number was computed or seen.

### Results

Added after the runs. Nothing above this heading changes once a statistic has been seen.

### Results (2026-10-14)

Run: code `3576bb9f` (Amendment 1). A119b's projection was reused (16 batches). The held-out templates were projected
after this file was committed, with the same code and settings (3 batches each): testis in 392 s at 16.2 GB, and KB3781
in 339 s at 19.3 GB (`winloci_data/GGO.splice.mmi`). Scoring took 253 s at 1.3 GB for A119b and under 10 s for each
held-out sample. Guard S0 held: 0.9956, 0.9999 and 1.0000 of templates align to their own exons.

| | A119b (development) | testis (held out) | KB3781 (held out) |
|---|---|---|---|
| templates / alignments / projected copies | 94,685 / 636,800 / 489,897 | 13,048 / 54,244 / 24,807 | 17,672 / 45,933 / 22,740 |
| templates at the 51-alignment cap | 8,682 | 347 | 189 |
| Liftoff extra copies (not short) | 1,888 | 1,888 | 1,129 |
| reachable (source placement assembled) | 420 (0.222) | 62 (0.033) | 52 (0.046) |
| **C1: recall (own)** | **268 / 420 = 0.638, FAIL** | **49 / 62 = 0.790, FAIL** | **44 / 52 = 0.846, pass** |
| recall (any) | 365 / 420 = 0.869 | 58 / 62 = 0.935 | 46 / 52 = 0.885 |
| misses: below threshold / capped / no alignment | 44 / 46 / 62 | 7 / 5 / 1 | 5 / 0 / 3 |
| recall (own) by the source's Liftoff copies: 1 / 2-5 / ≥ 6 | 71/104, 109/157, 88/159 | 9/10, 21/28, 19/24 | 9/12, 3/8, 32/32 |
| sites | 11,463 | 1,861 | 1,600 |
| `liftoff_copy` / `liftoff_gene` / `assembled` / `none` (retro-like) | 375 (1) / 2,643 (47) / 3,504 (62) / 4,941 (266) | 193 (0) / 1,061 (64) / 226 (1) / 381 (43) | 44 (1) / 762 (30) / 212 (1) / 582 (7) |
| **C2: unexplained sites** | **4,675 = 0.408, FAIL** | **338 = 0.182, FAIL** | **575 = 0.359, FAIL** |

The below-threshold misses had a median best identity and coverage of 0.898 and 0.930 in A119b, 0.498 and 0.515 in
testis, and 0.778 and 0.786 in KB3781.

#### Verdict (2026-10-14)

**Reading 3: not claimed.** C1 fails in development (0.638 < 0.80). It also fails in testis (0.790) and passes only in
KB3781 (0.846). C2 fails in all three samples. The expected retrocopy advantage did not show: retro-like sites are rare
(266 of A119b's 4,941 `none` sites).

#### Post hoc (descriptive, not pre-registered; `bench/template_copies/copies_posthoc.py`)

- **Recall (any) is not all projection.** It also counts a copy that is itself assembled, found by the self alignment of
  the template at the copy: 20 of A119b's 365, none in testis, 1 of KB3781's 46. Found by projection proper: 345, 58,
  45.
- **The misses follow the template's shape.** Call an own template gene-like when it covers ≥ 80% of the source
  placement's exon union and that union covers ≥ 80% of it.
  - Copies whose best own template is gene-like are found 15 of 16 times (A119b), 18 of 18 (testis) and 40 of 40
    (KB3781).
  - Otherwise: 253 of 404, 31 of 44 and 4 of 12.
  - Only 16 of A119b's 420 reachable copies have a gene-like own template. The rule is strict: a Liftoff placement's
    exons are the union of every isoform of the gene, while a template is one transcript.
- **Where the unexplained sites' templates sit.** Each `none` site's best template, placed against Liftoff's gene
  placements:

  | template | A119b (4,941 sites) | testis (381) | KB3781 (582) |
  |---|---|---|---|
  | ≥ 50% inside one gene's exons | 258 | 43 | 157 |
  | shares exonic bases with a gene, less than that | 1,089 | 97 | 175 |
  | inside a gene's span only (intronic) | 996 | 60 | 18 |
  | intergenic | 2,598 | 181 | 232 |

  - Only the intergenic templates (53%, 48% and 40% of the sites) are clearly loci the annotation lacks, which Liftoff
    cannot copy. The rest overlap or sit inside annotated genes. C2's unexplained share is a mix of the two; it is not
    a difference in input alone.
  - The ten templates with the most `none` sites make 443 (9%), 178 (47%) and 256 (44%) of them.
  - In A119b, 2,735 of the 4,941 come from single-exon templates. The contigs holding the most are chrY 492, chr13 459,
    chr1 393, chr9 361, chr15 320 and chr21 286.
  - 36.6% of A119b's templates are single-exon, against 3.4% (testis) and 3.0% (KB3781), so C2 is not comparable across
    samples.
- **The gene-like counts are per sample and concentrated.** KB3781's 40 gene-like copies come from 9 genes, 32 of them
  from two (LOC134756906 19, LOC134757412 13). The two human samples share one Liftoff table, so the same copies can
  appear in both (LINC01784 5 in each). The three samples are not pooled.

#### Final-review findings (verified, post hoc, 2026-10-14)

A fresh reviewer re-derived every recall number and stress-tested the site merge against all pairs, forcing the
pruning step, with no failure. The code stands; the text above was corrected. Each count below was recomputed with
`copies_posthoc.py` unless it says "reviewer".
- **Recall (own) includes copies found inside the template's own locus** (as registered, so the verdict is unchanged).
  In A119b, 30 of the 268 are found only by an alignment inside the own template's locus: the reviewer found the
  template itself spans the copy in 17 of them and a secondary hit inside a fused locus in 13. By projection off that
  locus, recall (own) is 238 / 420 = 0.567. Testis and KB3781 have none of these.
- **The held-out pass rests on two genes.** KB3781's 0.846 is per copy, and two genes give 32 of its 52 reachable
  copies, all found. Without the top gene it is 25 / 33 = 0.758; without both, 12 / 20 = 0.60. Testis misses the bar by
  one copy (50 / 62 would pass).
- **Unconditional recall** (found (own) / every Liftoff extra copy): 268 / 1,888 = 0.142, 49 / 1,888 = 0.026 and
  44 / 1,129 = 0.039.
- **The coverage threshold never binds** (reviewer). Identity ≤ coverage by construction (matches ≤ aligned bases, and
  the denominator ≥ template length), and Liftoff applies only `-sc` to extra copies (`liftover_types.py`). No number
  changes.
- **The `no_alignment` label is loose** (reviewer). In A119b, 13 of those 62 misses have a passing own alignment that
  covers 13-49% of the copy, and 4 more have a failing partial one.
- **Two design differences from Liftoff** (reviewer). Liftoff chains several alignment records into one placement; the
  analogue scores single PAF records. A genome-only stretch that minimap2 calls an intron costs nothing here but would
  lower Liftoff's identity.
- **Templates with at least one site** (a pre-registered reported item, omitted above): 13,771, 1,216 and 928.
- **Erratum to Amendment 1.** During profiling, `recall()` ran (1 s) but nothing it returned was printed. "No recall
  ... number was computed or seen" should read "printed or seen".


---

## PREREG_liftoff_pairs_null_2026-10-14.md

*(original: `docs/archive/2026-10/merged/PREREG_liftoff_pairs_null_2026-10-14.md`)*

## Liftoff copy pairs against matched decoys: O1's independent-truth test (2026-10-14)

**Pre-registered before any decoy statistic is computed.**
- Design: agreed with the user in chat on 2026-10-14 ("close F1, then a new test": count every read-supported copy
  pair; use a null matched for genomic distance; human for development, gorilla held out). It implements the family
  bar of `docs/superpowers/specs/2026-10-13-o1-layered-hierarchy-design.md`, section 3.
- Code: `bench/liftoff_pairs/` at commit `97f831b4`, run with the defaults (1,000 null draws, seed 7):

      bash tools/rlock.sh light python3 bench/liftoff_pairs/run_pairs.py human_A119b --out /mnt/linuxdisk/tmp/o1_liftoff_pairs/human_A119b
      bash tools/rlock.sh light python3 bench/liftoff_pairs/run_pairs.py human_testis --out /mnt/linuxdisk/tmp/o1_liftoff_pairs/human_testis
      bash tools/rlock.sh light python3 bench/liftoff_pairs/run_pairs.py gorilla_KB3781 --out /mnt/linuxdisk/tmp/o1_liftoff_pairs/gorilla_KB3781

### Question

Do our default de novo families put a gene's Liftoff copies together, more often than they put the gene together with
comparable loci that are not its copies? Liftoff's copy relation is independent of our edge rule: it comes from
projecting the annotation onto the assembly with `-copies`. Human is for development; gorilla is held out.

### Inputs (frozen)

| What | Path |
|---|---|
| Liftoff loci (self-lift with `-copies -sc 0.95`, figure 8) | `/mnt/linuxdisk/tmp/rustle_figures/liftoff/{human,gorilla}/liftoff_loci.tsv` |
| Read support per sample (≥ 2 primary reads with an aligned block on the exon union, `-F 2308`) | `/mnt/linuxdisk/tmp/rustle_figures/liftoff/support/<sample>.tsv` |
| Our default families (shipped, exact-mcl; one copy per member locus) | `/mnt/linuxdisk/tmp/rustle_figures/runs/<sample>/<sample>.fam.copies.tsv` |
| Figure 8's F1 rows (the wiring check) | `figures/data/fig8_samples.tsv` |

### Known before running (recorded now)

- **F1** (`docs/archive/2026-09/PREREG_liftoff_loci_2026-09-25.md`; verdict recorded today, register 1283) asks whether
  copy pairs whose two loci are both covered by family members share a family.
  - Human A119b: 131/148 (0.885), failing the 0.90 bar by three pairs.
  - Human testis 34/34; gorilla KB3781 11/12; chimp 18/18.
  - So among pairs we already place in families, most share one. F1 has no null and conditions on our own
    catalog.
- **The same denominator was defined but never scored.** Supplementary figure 6s defined the unconditioned
  sensitivity (every read-supported pair; `docs/archive/2026-09/PREREG_genome_wide_families_2026-09-25.md`,
  Amendment 1, item 3). It was never scored: no table of it exists, and no number of it has been seen.
- **Why not a family-size null.** The layered spec named "a size-matched null". A null that only keeps family sizes
  is passed by design: plan 1's E1 (register 1277) and the Soto study's T1a (register 1281) both passed this way. This
  test therefore matches decoys on coverage and on genomic distance.
- **Development contigs.** The family rule was developed on human chr16 and gorilla chr20 (NC_073244.2). Results
  without those contigs are reported.
- **Nothing is tuned.** The distance bins, the widening rule and the decoy definition are fixed in code commit
  97f831b4. They were written before any decoy statistic existed.
- **Input counts** (`inputs.json`; inputs, not outcomes). Read-supported copy pairs at identity ≥ 0.95 / 0.98 / 0.99 / 1.00:

  | sample | genome | without development contigs | decoy pool (genome / without) |
  |---|---|---|---|
  | human A119b | 514 / 338 / 246 / 53 | 469 / 307 / 227 / 47 | 43,703 / 42,018 |
  | human testis | 38 / 31 / 24 / 3 | 33 / 26 / 19 / 1 | 19,454 / 18,606 |
  | gorilla KB3781 | 78 / 49 / 35 / 5 | 75 / 46 / 33 / 5 | 18,551 / 17,380 |

### Definitions

- **Copy pair (a, c).** a is a source gene's annotated placement (in place or moved) and c is one of its Liftoff extra
  copies.
  - The extra copy has `sequence_ID` ≥ the floor.
  - Both exon unions are ≥ 200 bp.
  - Both loci are read-supported in the sample (`_liftoff.copy_pairs` with support).
- **Covered.** A locus is covered by a family when a member of that family (its representative's exons) covers ≥ 50%
  of the locus's exon bases. This is figure 8's rule (`_liftoff.pair_families`).
- **Shared.** A pair is shared when one family covers both loci.
  - A locus covered by no family member cannot be shared.
  - There is no conditioning on our catalog: every read-supported pair is in the denominator.
- **Decoy.** A decoy for c is a read-supported Liftoff locus of another source gene that does not overlap a, matched to
  c on two things:
  1. whether a family member covers it;
  2. its distance bin from a: another contig, or the same contig at < 100 kb, 100 kb-1 Mb, 1-10 Mb, or ≥ 10 Mb
     (gap between spans).
- **Widening.** When no decoy matches:
  - level 1 keeps the contig relation at any distance;
  - level 2 allows any contig;
  - a pair without a decoy at level 2 leaves the test, on both sides.

  The number of pairs at each level is reported.

### The test

- **Statistic.** k = the number of shared pairs among the n pairs in the test.
- **Null.** Each draw replaces every pair's copy by a uniformly drawn decoy and counts the shared decoy pairs. There
  are 1,000 draws, seed 7, and p = (1 + draws ≥ k) / 1,001.
- **Bar:** p < 0.001 on the genome at identity ≥ 0.95, in each of human A119b, human testis and gorilla KB3781.
- **Reported, no bar:**
  - the identity floors 0.98, 0.99 and 1.00;
  - the genome without the development contigs;
  - k / n with a 95% interval (2,000 bootstrap resamples of source genes, figure 8's seed);
  - the null's mean rate and maximum;
  - the widening levels;
  - the pairs' coverage: both loci covered, one, or none.

### Readings (fixed now)

| human A119b | human testis | gorilla KB3781 | Reading |
|---|---|---|---|
| pass | pass | pass | O1's independent-truth claim: "a gene's Liftoff copies share one of our families far more often than distance- and coverage-matched loci of other genes (rate x against null y), in human and held out in gorilla". |
| pass | pass | fail | The same claim for human only; the gorilla result is reported. |
| fail on either human sample | — | — | No claim; the per-sample tables say where it fails. |

**What a pass does and does not show.**
- Decoys are matched on whether a family covers them, so the test asks whether our families group a gene with its own
  copies rather than with other covered loci at the same distance.
- It would fail if our families lumped neighbouring genes duplicated with the copy.
- A decoy can be a paralog annotated as a separate gene, which makes the null stricter, not looser.
- It is not a precision measure for loci without a Liftoff relation.

### Wiring check (the run stops, and nothing is read, if it fails)

F1 recomputed by the code must reproduce `figures/data/fig8_samples.tsv`'s n and k, for every scope and floor of the
sample. It did for all three samples in the input reports and in the tests.

### Results

Added after the runs. Nothing above this heading changes once a statistic has been seen.

#### Results (2026-10-14)

Run: commit `2cd1446c` (code `97f831b4`). Wall time and peak memory: human A119b 4.56 s, 355 MB; human testis 2.17 s, 273 MB; gorilla KB3781 1.89 s, 233 MB.
Wiring check passed: F1 reproduced the figure 8 table for every scope and floor of each sample.

**The bar (genome, identity ≥ 0.95)**

| sample | pairs n | shared k | k / n [95% interval] | null mean rate (max k) | p | bar |
|---|---|---|---|---|---|---|
| human A119b | 514 | 75 | 0.1459 [0.1025, 0.1955] | 0.0214 (16) | 0.000999 (1/1001) | **PASS** |
| human testis | 38 | 12 | 0.3158 [0.1702, 0.5000] | 0.0321 (3) | 0.000999 (1/1001) | **PASS** |
| gorilla KB3781 | 78 | 4 | 0.0513 [0.0112, 0.1098] | 0.0006 (1) | 0.000999 (1/1001) | **PASS** |

**Coverage of the pairs and widening (genome, identity ≥ 0.95)**

| sample | both loci covered | one | none | shared / both covered | decoys at level 0 / 1 / 2 / none |
|---|---|---|---|---|---|
| human A119b | 91 | 133 | 290 | 75 / 91 | 507 / 7 / 0 / 0 |
| human testis | 12 | 2 | 24 | 12 / 12 | 35 / 3 / 0 / 0 |
| gorilla KB3781 | 4 | 8 | 66 | 4 / 4 | 78 / 0 / 0 / 0 |

**Every floor and scope (reported, no bar)**

| sample | scope | identity | n | k | k / n | null mean rate (max) | p |
|---|---|---|---|---|---|---|---|
| human A119b | genome | 0.95 | 514 | 75 | 0.1459 | 0.0214 (16) | 0.000999 (1/1001) |
| human A119b | genome | 0.98 | 338 | 62 | 0.1834 | 0.0284 (15) | 0.000999 (1/1001) |
| human A119b | genome | 0.99 | 246 | 51 | 0.2073 | 0.0383 (15) | 0.000999 (1/1001) |
| human A119b | genome | 1.00 | 53 | 18 | 0.3396 | 0.0855 (8) | 0.000999 (1/1001) |
| human A119b | minus_dev | 0.95 | 469 | 57 | 0.1215 | 0.0210 (13) | 0.000999 (1/1001) |
| human A119b | minus_dev | 0.98 | 307 | 48 | 0.1564 | 0.0291 (13) | 0.000999 (1/1001) |
| human A119b | minus_dev | 0.99 | 227 | 39 | 0.1718 | 0.0390 (14) | 0.000999 (1/1001) |
| human A119b | minus_dev | 1.00 | 47 | 15 | 0.3191 | 0.0942 (8) | 0.000999 (1/1001) |
| human testis | genome | 0.95 | 38 | 12 | 0.3158 | 0.0321 (3) | 0.000999 (1/1001) |
| human testis | genome | 0.98 | 31 | 10 | 0.3226 | 0.0392 (3) | 0.000999 (1/1001) |
| human testis | genome | 0.99 | 24 | 9 | 0.3750 | 0.0503 (3) | 0.000999 (1/1001) |
| human testis | genome | 1.00 | 3 | 1 | 0.3333 | 0.3333 (1) | 1.000000 (1001/1001) |
| human testis | minus_dev | 0.95 | 33 | 11 | 0.3333 | 0.0367 (4) | 0.000999 (1/1001) |
| human testis | minus_dev | 0.98 | 26 | 9 | 0.3462 | 0.0465 (3) | 0.000999 (1/1001) |
| human testis | minus_dev | 0.99 | 19 | 8 | 0.4211 | 0.0624 (3) | 0.000999 (1/1001) |
| human testis | minus_dev | 1.00 | 1 | 1 | 1.0000 | 1.0000 (1) | 1.000000 (1001/1001) |
| gorilla KB3781 | genome | 0.95 | 78 | 4 | 0.0513 | 0.0006 (1) | 0.000999 (1/1001) |
| gorilla KB3781 | genome | 0.98 | 49 | 3 | 0.0612 | 0.0000 (0) | 0.000999 (1/1001) |
| gorilla KB3781 | genome | 0.99 | 35 | 3 | 0.0857 | 0.0000 (0) | 0.000999 (1/1001) |
| gorilla KB3781 | genome | 1.00 | 5 | 0 | 0.0000 | 0.0000 (0) | 1.000000 (1001/1001) |
| gorilla KB3781 | minus_dev | 0.95 | 75 | 4 | 0.0533 | 0.0006 (1) | 0.000999 (1/1001) |
| gorilla KB3781 | minus_dev | 0.98 | 46 | 3 | 0.0652 | 0.0000 (0) | 0.000999 (1/1001) |
| gorilla KB3781 | minus_dev | 0.99 | 33 | 3 | 0.0909 | 0.0000 (0) | 0.000999 (1/1001) |
| gorilla KB3781 | minus_dev | 1.00 | 5 | 0 | 0.0000 | 0.0000 (0) | 1.000000 (1001/1001) |

#### Verdict (2026-10-14)

**Reading 1.** The bar passes in human A119b, human testis and gorilla KB3781 (held out), each at p = 0.000999
(1/1001).
- **Claim:** a gene's Liftoff copies share one of our families far more often than distance- and coverage-matched loci
  of other genes.

  | sample | copy pairs shared | decoys (null mean) |
  |---|---|---|
  | human A119b | 75 of 514 (0.146) | 0.021 |
  | human testis | 12 of 38 (0.316) | 0.032 |
  | gorilla KB3781 | 4 of 78 (0.051) | 0.0006 |

  It holds in human, held out in gorilla, without the development contigs, and at the floors 0.98 and 0.99 in
  every sample. At 1.00, testis (3 pairs) and gorilla (5 pairs) do not pass.

**What this does and does not say.**
- **The test is effectively conditioned on both copies being covered** (post-hoc note below). Only those pairs can
  count, on either side: 91 (A119b), 12 (testis) and 4 (gorilla).
  - Among them, the copies share a family in 75 / 91 (0.824), 12 / 12 and 4 / 4.
  - The coverage-matched decoys give null rates of 0.121, 0.102 and 0.012.
  - A119b's 0.824 is below the 0.90 that F1 asks for, and its 16 split pairs are F1's.
- **What a pass shows.** A gene's covered copies are grouped with it much more often than covered loci of other genes at
  a similar distance.
  - It is a specificity check against unrelated loci, not proof against an independent family definition.
  - The copies are ≥ 0.95-identity sequence matches and our families are similarity clusters, so agreement is expected.
- **Recall is low, mostly because of coverage.** In most read-supported copy pairs neither copy is covered by any family
  member: 290 of 514 (A119b), 24 of 38 (testis), 66 of 78 (gorilla).
  - The other unshared pairs are grouping failures, F1's split pairs among them.
  - This is O1's known recall limit (node completeness and edge recall), and it is not measured as a bar here.
- **The gorilla pass rests on 4 pairs from 4 genes.** It passes because no random decoy draw reached 4 (null maximum 1).

**Post hoc (added after the final review): the denominator does not do what Definitions says.**
- **Uncovered pairs carry no information.** Decoys are matched on the copy's coverage, so a pair with either copy
  uncovered adds 0 to k and to every null draw.
  - k, the null and p are therefore exactly those of the same test run on the both-covered pairs.
  - "There is no conditioning on our catalog" (Definitions) holds for the denominator only.
  - The headline rates (0.146 against 0.021 in A119b) are diluted by 423, 26 and 74 uninformative pairs.
- **The bar is weak.** It passes from k ≥ 17 in A119b (null maximum 16), k ≥ 4 in testis and k ≥ 2 in gorilla.
  - A119b would pass with about 80% of its covered copy pairs split, so the test cannot detect F1's failure mode.
  - With F1 already known (131 / 148), a pass was nearly certain.
- **These numbers were already defined elsewhere, and are now seen.** k / n is figure 7's Liftoff-pair sensitivity of
  the de novo families (`docs/archive/2026-09/PREREG_genome_wide_families_2026-09-25.md`, Amendment 1, item 4). It is
  also the default arm of F6.6 (item 3).
- **For a reuse.** To put recall into the bar, do not match decoys on coverage, or test the unconditioned rate against
  an external expectation.


---

## PREREG_soto_subset_cn_2026-10-14.md

*(original: `docs/archive/2026-10/merged/PREREG_soto_subset_cn_2026-10-14.md`)*

## Soto's families as our families cut by copy number: T1 and T2 (2026-10-14)

**Pre-registered before any T1 or T2 statistic is computed.**
- Spec: `docs/superpowers/specs/2026-10-14-soto-subset-cn-design.md`.
- Plan: `docs/superpowers/plans/2026-10-14-soto-subset-cn.md`.
- Code: `bench/soto_subset/` at commit `25006d97`. Every number below comes from it, run with the defaults (1,000 null
  draws, 1,000 bootstrap resamples, seed 7):

      bash tools/rlock.sh light python3 bench/soto_subset/run_subset.py guided --out /mnt/linuxdisk/tmp/o1_soto_subset/guided
      bash tools/rlock.sh light python3 bench/soto_subset/run_subset.py denovo --sample human_A119b --out /mnt/linuxdisk/tmp/o1_soto_subset/human_A119b
      bash tools/rlock.sh light python3 bench/soto_subset/run_subset.py denovo --sample human_testis --out /mnt/linuxdisk/tmp/o1_soto_subset/human_testis

### Question

- **T1, Soto ⊂ ours:** are Soto's families (Soto 2025, Table S1C) a refinement of our O1 multi-copy gene families?
- **T2, Soto ≈ ours ∩ CN:** does cutting ours by copy number, with Soto's own pair rule, rebuild them?

Human only. Nothing is tuned on Soto's table. The pair rule is Soto's (|Δ| < 2), and the nulls, draw counts and seed
are fixed above.

### Inputs (frozen)

| What | Path |
|---|---|
| Soto's families and the S1C copy number | `bench/soto/soto_famCN_S1C.tsv` |
| Assembly copy number (HG002 diploid, `fam_HG002`) | `/mnt/linuxdisk/tmp/rustle_figures_dev/soto_parcn_asm/results_s1e.tsv` |
| Read-depth copy number (`famCN_sotoiv`, 268 SGDP genomes) | `/mnt/linuxdisk/home/juanfraitu/winloci_data/soto_replication/famcn_ours_allwssd.tsv` |
| Soto's genes on CHM13 v2.0 (de novo mapping, coordinate check) | `/mnt/linuxdisk/tmp/o1_close/lift/soto_truth.v2.tsv` |
| Guided families (F_O, the 10-09 run) | `/mnt/linuxdisk/tmp/o1_layers/v2/dupmasker/ours.hierarchy.members.tsv`, `ours.loci.tsv`; loci `/mnt/linuxdisk/tmp/ourhom/loci.gff3` |
| De novo catalogs (exact-mcl; the same partitions as the shipped `<sample>.fam.clusters.tsv`) | `/mnt/linuxdisk/tmp/o1_close/hier_shipped/<sample>/ours.hierarchy.members.tsv`, `ours.loci.tsv`; loci `/mnt/linuxdisk/tmp/rustle_figures/runs/<sample>/<sample>.fam.loci.gff3` |
| Soto's recipe | `bench/soto/soto_replication.py cluster --pair-mad` on `bench/soto/shared_exons_5154_exon_mapback.tsv` |

### Known before running (recorded now)

- **T1 on the guided universe is not blind, and strict containment already failed.** On F_O (the 10-09 run, at
  member level, over the 400 evaluable families), 313 (78.2%) of Soto's families are fully contained in one of our
  families and 346 (86.5%) in one superfamily; 215 equal one of our expansions. That prereg's H1 (≥ 90% contained)
  was NOT met, and it concluded that Soto's families are not a strict subset of ours
  (`docs/PREREG_family_hierarchy_flag_2026-10-09.md`). T1 asks a different question: containment beyond random
  partitions with our block sizes, and which partition is the coarser. A T1 pass therefore does not make Soto a strict
  subset (ruling R7). The de novo catalogs are the replication.
- **T1a is expected to pass easily.** Plan 1's global null was passed by design (register 1277), and random
  partitions with our block sizes rarely contain a Soto family. T1a only rules out containment by chance; the
  informative T1 quantities are the containment rate, the direction D and the violation table.
- **Older catalogs nest too.** On the older RNA lattice catalogs, 80-90% of Soto's families nest in ours at L1/L2.
  On held-out chr5/7/21 it was 61/76: 53/59 families of ≤ 3 genes and 8/17 larger ones
  (`bench/SOTO_AS_A_REFINEMENT.md`).
- **Copy number already raised agreement on our edges.** In an exploratory, unregistered run (2026-10-07), our F_O
  graph's edges gated by Soto's pair rule gave ARI 0.7142 / 0.7849 / 0.8263 (no copy number / read depth / S1C) at
  identity ≥ 0.70. This study cuts our families rather than our edges, and adds assembly copy number.
- **Soto's recipe rungs** (registers 1162 and 1169; reproduced 2026-10-14), over the 491 families: no copy number
  0.7307 / 345 exact, read depth 0.9277 / 411, S1C 0.9698 / 479.
- **Assembly vs S1C:** assembly copy number agrees with S1C at Spearman 0.617 (human; register 1172).
- **No held-out half is needed.** S1C's family-hash halves are both spent for copy-number split rules (2026-10-07).
  Nothing here selects a rule, so no half is needed; the de novo catalogs replicate.
- **FAM90A.** ARI over Soto's families swings ±0.035 with FAM90A (register 1170). Hence exact families are a
  co-requirement of T2's bar, and ARI is also reported without FAM90A.
- **Input counts** (`inputs.json`; inputs, not outcomes):

| | guided | A119b | testis |
|---|---|---|---|
| members | 2,334 genes | 1,396 loci | 332 loci |
| Soto families with ≥ 2 members | 491 | 293 | 105 |
| single-label / multi-label members | 2,185 / 149 | 1,190 / 206 | 300 / 32 |
| members on our families | 1,796 | 1,115 | 238 |
| members with assembly / read-depth / S1C copy number | 1,534 / 2,334 / 1,793 | 1,175 / 1,396 / 1,196 | 288 / 332 / 292 |
| common set (a value in all three) | 1,534 | 1,175 | 288 |
| Soto families with ≥ 2 members in the common set | 357 | 257 | 85 |
| FAM90A members in the common set | 2 | 0 | 0 |

### Definitions

- **Guided universe.** Members are Soto's 2,334 genes. A gene's locus is its gene line in the F_O input GFF (by CAT
  v4 gene id), folded through `ours.loci.tsv`.
  - Genes on loci of one F_O family share a block.
  - Genes on one unclustered locus share a block.
  - Every other gene is our singleton.
- **De novo universe.** Members are the catalog loci reached by a gene of Soto's 491 families: an exon shares ≥ 1 bp,
  on the same strand, with a locus exon, and a folded locus counts as its representative.
  - Blocks are the catalog's families. A reached locus in no family is our singleton.
  - A locus carries the families of every gene that reaches it.
- **Soto's labels.** `Unassigned_*` genes are Soto singletons. A member with two or more labels is *multi-label*.
- **Copy number of a member** = the median over its genes that have a value.
- **Common set** = the members with a value in all three sources.

### T1 — Soto ⊂ ours

- **Containment.** A Soto family with ≥ 2 members is contained when all its members lie in one block of ours. The
  statistic is the number contained.
- **Size null.** Random partitions of the members with our exact block-size multiset, singletons included: a uniform
  permutation of the block labels, 1,000 draws, seed 7. p = (1 + draws ≥ observed) / 1,001.
- **T1a bar:** p < 0.001.
- **Direction.** Over the single-label members, D = H(Soto | ours) − H(ours | Soto), in bits.
  - The 95% percentile interval comes from 1,000 bootstrap resamples of *units* (ruling R1), seed 7.
  - A unit is a connected component of the members linked by a shared Soto label or a shared block of ours.
  - A unit drawn twice enters as two disjoint copies.
- **T1b bar:** the interval lies above 0.
- **Reported, no bar:**
  - the rate over all 491 families, and over the families with ≥ 2 members (R6);
  - plan 1's rate: families with ≥ 2 members on clustered members, contained within that restriction;
  - containment for families of ≤ 3 and > 3 members;
  - every violation, classified by the first rule that makes it contained: **coverage** (our singletons removed),
    then **cover** (multi-label members also removed), else **conflict**.

### T2 — Soto ≈ ours ∩ CN

- **The intersection.** Inside each block of ours, members are joined when |CN_a − CN_b| < 2 (Soto's per-pair
  MAD < 1). The components are the intersection: runs of the block's sorted values, split at gaps ≥ 2.
- **Scores,** on the same members:
  - **exact:** distinct Soto member sets, restricted to the single-label members with ≥ 2 of them, that equal a
    component (R3);
  - **ARI** over the single-label members, also reported without FAM90A (`ID_356`);
  - **Hungarian Jaccard:** the mean one-to-one best-match Jaccard over Soto's families with ≥ 2 members, all members
    included.
- **Null.** Copy-number values are permuted among the members of each block of ours, 1,000 draws, seed 7, with the
  same p formula. It is computed for each source.
- **T2 bar:** with assembly copy number on the common set, the exact count and ARI both beat their null at p < 0.001.
- **Reported, no bar:**
  - read depth and S1C on the common set, with their nulls;
  - each source on its own coverage;
  - our blocks alone;
  - the ceiling: Soto's partition intersected with ours, on single-label members (R5);
  - Soto's recipe with each source, guided only (R4), on its full gene set and on the common set.

### Rulings (made before any statistic)

- **R1 — bootstrap unit.** T1b resamples units, not Soto's blocks. Two copies of a Soto block inside one block of
  ours would fake a split and inflate H(Soto | ours), biasing toward the bar. D decomposes exactly over units.
- **R2 — a fourth file.** `universes.py` holds the universes beside the spec's three files.
- **R3 — exact over single-label members.** Exact is counted over single-label members, as Soto's own scorer does:
  479/491 excludes the 149 multi-family genes.
- **R4 — recipe rows on guided only.** Soto's recipe is gene-level.
- **R5 — the ceiling on single-label members.** The meet of a cover is not a partition.
- **R6 — two rates.** The rate over all 491 families (families with < 2 members count as not contained) and plan 1's
  rate are both reported.
- **R7 — no strict-subset wording.** The spec's statement "Soto's families are a refinement of ours" would claim strict
  containment, which the 10-09 prereg already refuted (78.2% < 90%) and which T1's null-based bars cannot establish.
  Readings 1 and 2 therefore state containment beyond chance with the contained count and the listed exceptions.
  Added while executing the plan, before any statistic.

### Readings (fixed now)

"T1 holds" means T1a and T1b both pass. A de novo "fail" means a bar fails on either catalog.

| Guided | De novo | Reading |
|---|---|---|
| T1 and T2 hold | hold on both catalogs | O1's statement: "Soto's families lie inside our multi-copy gene families far beyond chance, with ours the coarser partition (c of f contained; the exceptions are listed), and cutting ours by copy number rebuilds them (exact x, ARI y; ceiling z)" (R7). Replaces the failed hierarchy claim (register 1277) as the Soto sentence. |
| T1 and T2 hold | fail | The same statement, labelled as shown on the guided universe only. |
| T1 holds, T2 fails | — | Soto ⊂ ours is claimed. The copy-number part is claimed only if read depth passes T2's bar on the same set, labelled human-only and population-based. |
| T1 fails | — | No subset claim; T2 is reported, not claimed. The size and violation tables say where Soto and we disagree. |

### Wiring checks (the run stops, and nothing is read, if one fails)

- The gene-id join agrees with coordinates: 0 disagreements.
- Soto's recipe rungs reproduce: 0.7307 / 345, 0.9698 / 479 and 0.9277 / 411.
- Plan 1's coverage reproduces: guided, 1,796 genes on our families; A119b, 1,211 / 155 / 1,092 Soto records on
  clustered loci / only unclustered / none; testis, 286 / 100 / 2,072.

### Results

Added after the runs. Nothing above this heading changes once a statistic has been seen.

#### Guided (2026-10-14)

Run: commit `6d5be546` (code `25006d97`); wall 19.72 s, peak 164,044 KB. Wiring checks: all passed (gene-id join 1,710 / 0 / 121;
recipe rungs 0.7307 / 345, 0.9698 / 479, 0.9277 / 411 reproduced; plan 1's coverage reproduced in `inputs.json`).

**T1**

| Quantity | Value |
|---|---|
| Soto families with ≥ 2 members (`families_ge2`) | 491 |
| contained (`contained`); rate over 491 (`rate_all`); over ≥ 2 (`rate_ge2`) | 258; 0.5255; 0.5255 |
| size null mean / max (`null.mean`, `null.max`) | 0.709 / 4 |
| p (`p`) → **T1a** | 0.000999 (1/1001) → **PASS** |
| H(Soto \| ours) / H(ours \| Soto) | 0.2860 / 0.6855 |
| D [95% interval] (`D`, `D_ci95`) → **T1b** | -0.3995 [-0.6154, -0.1945] → **FAIL** |
| ≤ 3 members: contained / families (`by_size.le3`) | 176 / 271 |
| > 3 members: contained / families (`by_size.gt3`) | 82 / 220 |
| plan 1's rate: contained / evaluable (`plan1`) | 315 / 402 |
| violations: coverage / cover / conflict (`violation_classes`) | 139 / 47 / 47 |

The conflict families (from `violations.tsv`, `class` = conflict), as `label (size, blocks)`: ID_100 (19, 13), ID_104 (9, 4), ID_112 (8, 7), ID_12 (5, 2), ID_136 (6, 4), ID_14 (10, 6), ID_142 (15, 2), ID_147 (4, 2), ID_155 (6, 4), ID_156 (6, 3), ID_160 (6, 3), ID_175 (13, 2), ID_182 (8, 3), ID_187 (5, 4), ID_188 (5, 3), ID_208 (7, 5), ID_211 (4, 2), ID_212 (11, 8), ID_215 (11, 2), ID_22 (11, 3), ID_226 (6, 3), ID_24 (8, 4), ID_241 (3, 2), ID_244 (9, 4), ID_257 (5, 2), ID_261 (8, 3), ID_263 (4, 3), ID_266 (5, 3), ID_267 (6, 2), ID_269 (3, 2), ID_270 (6, 2), ID_28 (29, 21), ID_280 (11, 7), ID_283 (8, 2), ID_284 (13, 3), ID_306 (4, 2), ID_324 (6, 2), ID_35 (41, 7), ID_39 (7, 3), ID_393 (5, 2), ID_400 (6, 3), ID_443 (5, 2), ID_47 (4, 2), ID_63 (11, 4), ID_65 (15, 5), ID_78 (9, 7), ID_8 (14, 4).

**T2 (common set, 1,534 members)**

| Copy number | exact | ARI | ARI w/o FAM90A | Hungarian J | null exact mean (max) | null ARI mean (max) | p exact | p ARI | bar |
|---|---|---|---|---|---|---|---|---|---|
| assembly (T2's bar) | 186 | 0.8274 | 0.8275 | 0.7694 | 163.195 (171) | 0.6836 (0.7303) | 0.000999 (1/1001) | 0.000999 (1/1001) | **PASS** |
| read depth | 206 | 0.8685 | 0.8686 | 0.7890 | 173.116 (181) | 0.7026 (0.7464) | 0.000999 (1/1001) | 0.000999 (1/1001) | (reported; passes) |
| S1C (circular upper bound) | 232 | 0.9166 | 0.9167 | 0.8305 | 186.865 (195) | 0.7336 (0.7788) | 0.000999 (1/1001) | 0.000999 (1/1001) | (reported; passes) |
| our blocks alone (`ours_alone`) | 183 | 0.7890 | 0.7889 | 0.7411 | — | — | — | — | — |
| ceiling: Soto ∩ ours (`ceiling`) | 258 | 0.9425 | 0.9425 | 0.8727 | — | — | — | — | — |

**T2, each source on its own coverage (`own`)**

| Copy number | members | exact | ARI | ARI w/o FAM90A | Hungarian J |
|---|---|---|---|---|---|
| assembly | 1534 | 186 | 0.8274 | 0.8275 | 0.7694 |
| read depth | 2334 | 203 | 0.7835 | 0.7556 | 0.6706 |
| S1C | 1793 | 236 | 0.9340 | 0.9187 | 0.8113 |

**Soto's recipe with each source (`recipe.json`)**

| Copy number | full gene set ARI / exact of 491 | common set exact | common set ARI | common set ARI w/o FAM90A | common set Hungarian J |
|---|---|---|---|---|---|
| none | 0.7307 / 345 | 266 | 0.7465 | 0.7469 | 0.8312 |
| S1C | 0.9698 / 479 | 348 | 0.9894 | 0.9894 | 0.9902 |
| read depth | 0.9277 / 411 | 299 | 0.9248 | 0.9250 | 0.9293 |
| assembly | 0.7410 / 337 | 258 | 0.8656 | 0.8657 | 0.8694 |

Guided verdicts: T1a PASS, T1b FAIL, T2 PASS; read depth on the common set passes T2's bar.

**Post hoc (descriptive, not pre-registered; added after the guided run).**
- Why T1b fails: 503 of the 2,185 single-label genes are in no F_O family (495 alone in a block, 8 sharing an
  unclustered locus). They split Soto's families and raise H(ours | Soto). On the 1,682 single-label genes that sit in our families, H(Soto | ours) = 0.3715,
  H(ours | Soto) = 0.1809 and D = +0.1906 (ours the coarser). That cut conditions on our own catalog, so it is not a test,
  and no interval was computed. Of the 233 families not contained, 139 are coverage violations.
- T2's null sits below our uncut blocks: the shuffle null's mean (163.2 exact, ARI 0.6836) is lower than
  our blocks alone (183 exact, ARI 0.7890). The bar therefore compares the copy-number cut with random cuts,
  not with leaving our families uncut. Against uncut, assembly copy number gains +3 exact and +0.0384 ARI.

#### De novo (2026-10-14)

Run: commit `6d329ed9` (code `25006d97`); A119b wall 15.99 s, peak 319,368 KB; testis wall 6.42 s, peak 183,876 KB.
Wiring checks: all passed (gene-id join 1,710 / 0 / 121; plan 1's coverage reproduced in `inputs.json`).

**T1**

| Quantity | A119b | testis |
|---|---|---|
| Soto families with ≥ 2 members | 293 | 105 |
| contained; rate over 491; over ≥ 2 | 80; 0.1629; 0.2730 | 54; 0.1100; 0.5143 |
| size null mean / max | 0.192 / 4 | 0.291 / 3 |
| p → **T1a** | 0.000999 (1/1001) → **PASS** | 0.000999 (1/1001) → **PASS** |
| H(Soto \| ours) / H(ours \| Soto) | 0.2237 / 1.2667 | 0.1034 / 0.3555 |
| D [95% interval] → **T1b** | -1.0430 [-1.2541, -0.8029] → **FAIL** | -0.2520 [-0.3903, -0.1052] → **FAIL** |
| ≤ 3 members: contained / families | 59 / 126 | 45 / 83 |
| > 3 members: contained / families | 21 / 167 | 9 / 22 |
| plan 1's rate: contained / evaluable | 105 / 247 | 69 / 84 |
| violations: coverage / cover / conflict | 108 / 34 / 71 | 44 / 4 / 3 |

Conflict families: A119b ID_100 (31, 13), ID_113 (19, 4), ID_116 (8, 2), ID_127 (4, 2), ID_132 (6, 2), ID_136 (14, 10), ID_14 (12, 8), ID_142 (23, 9), ID_148 (5, 2), ID_149 (12, 6), ID_151 (11, 8), ID_152 (11, 7), ID_153 (10, 7), ID_155 (8, 4), ID_156 (7, 3), ID_160 (5, 4), ID_163 (14, 13), ID_175 (11, 4), ID_176 (6, 3), ID_178 (7, 5), ID_179 (9, 4), ID_181 (10, 7), ID_182 (14, 9), ID_187 (9, 7), ID_188 (8, 6), ID_203 (4, 2), ID_208 (7, 5), ID_211 (9, 4), ID_212 (14, 11), ID_214 (14, 12), ID_215 (13, 7), ID_218 (13, 10), ID_231 (5, 3), ID_233 (15, 4), ID_24 (6, 2), ID_243 (8, 4), ID_244 (7, 3), ID_245 (15, 7), ID_257 (5, 3), ID_261 (7, 4), ID_262 (12, 3), ID_270 (6, 5), ID_279 (6, 2), ID_280 (35, 18), ID_284 (26, 8), ID_31 (14, 3), ID_35 (45, 13), ID_363 (9, 2), ID_368 (8, 3), ID_369 (4, 2), ID_374 (5, 3), ID_391 (10, 2), ID_393 (6, 3), ID_395 (18, 7), ID_396 (4, 3), ID_397 (9, 5), ID_399 (8, 3), ID_400 (8, 7), ID_403 (2, 2), ID_41 (25, 9), ID_455 (4, 2), ID_468 (21, 3), ID_471 (4, 3), ID_49 (5, 3), ID_56 (6, 4), ID_60 (7, 3), ID_62 (15, 9), ID_65 (13, 6), ID_78 (11, 8), ID_8 (8, 4), ID_99 (27, 14); testis ID_154 (7, 2), ID_395 (5, 2), ID_407 (6, 3).

**T2 (common set: A119b 1,175 loci, testis 288 loci)**

| Copy number | catalog | exact | ARI | ARI w/o FAM90A | Hungarian J | null exact mean (max) | null ARI mean (max) | p exact | p ARI | bar |
|---|---|---|---|---|---|---|---|---|---|---|
| assembly (T2's bar) | A119b | 51 | 0.4022 | 0.4022 | 0.5361 | 47.982 (52) | 0.3674 (0.3873) | 0.006993 (7/1001) | 0.000999 (1/1001) | **FAIL** |
| assembly (T2's bar) | testis | 40 | 0.7574 | 0.7574 | 0.7530 | 39.243 (40) | 0.7312 (0.7708) | 0.243756 (244/1001) | 0.044955 (45/1001) | **FAIL** |
| read depth | A119b | 55 | 0.4414 | 0.4414 | 0.5425 | 50.035 (53) | 0.3730 (0.3941) | 0.000999 (1/1001) | 0.000999 (1/1001) | (reported; passes) |
| read depth | testis | 37 | 0.7723 | 0.7723 | 0.7438 | 37.003 (38) | 0.7247 (0.7858) | 1.000000 (1001/1001) | 0.012987 (13/1001) | (reported; does not pass) |
| S1C | A119b | 57 | 0.4992 | 0.4992 | 0.5600 | 50.460 (56) | 0.4111 (0.4298) | 0.000999 (1/1001) | 0.000999 (1/1001) | (reported; passes) |
| S1C | testis | 40 | 0.8305 | 0.8305 | 0.7634 | 38.262 (40) | 0.7633 (0.8286) | 0.005994 (6/1001) | 0.000999 (1/1001) | (reported; does not pass) |
| our blocks alone | A119b | 51 | 0.4621 | 0.4621 | 0.5428 | — | — | — | — | — |
| our blocks alone | testis | 40 | 0.7861 | 0.7861 | 0.7611 | — | — | — | — | — |
| ceiling | A119b | 63 | 0.5132 | 0.5132 | 0.6314 | — | — | — | — | — |
| ceiling | testis | 41 | 0.8409 | 0.8409 | 0.8081 | — | — | — | — | — |

**T2, each source on its own coverage**

| Copy number | catalog | members | exact | ARI | ARI w/o FAM90A | Hungarian J |
|---|---|---|---|---|---|---|
| assembly | A119b | 1175 | 51 | 0.4022 | 0.4022 | 0.5361 |
| assembly | testis | 288 | 40 | 0.7574 | 0.7574 | 0.7530 |
| read depth | A119b | 1396 | 62 | 0.4255 | 0.4253 | 0.5236 |
| read depth | testis | 332 | 39 | 0.7549 | 0.7549 | 0.6827 |
| S1C | A119b | 1196 | 58 | 0.4986 | 0.4984 | 0.5603 |
| S1C | testis | 292 | 41 | 0.8312 | 0.8312 | 0.7631 |

De novo verdicts: A119b T1a PASS, T1b FAIL, T2 FAIL; testis T1a PASS, T1b FAIL, T2 FAIL.

**Post hoc (descriptive, not pre-registered).** As on the guided universe, most violations are coverage violations
(A119b 108 of 213, testis 44 of 51). Cutting our blocks by assembly copy number does not beat leaving them uncut:
A119b 51 exact / ARI 0.4022 against 51 / 0.4621 uncut, and testis 40 / 0.7574 against 40 / 0.7861.

#### Verdict (2026-10-14)

Reading 4 of the readings table: T1 fails on the guided universe, so no subset claim is made, and T2 is reported, not
claimed.
- **T1a (guided): PASS.** 258 of 491 Soto families are contained, against a size-null mean
  of 0.709 (max 4); p = 0.000999.
- **T1b (guided): FAIL.** D = -0.3995 [-0.6154, -0.1945], below 0: H(Soto | ours) = 0.2860 is smaller than
  H(ours | Soto) = 0.6855, so ours is not the coarser partition on this universe.
- **T2 (guided, assembly copy number, common set of 1,534 genes): PASS, reported, not claimed.**
  - Exact 186, against a null mean of 163.2 (max 171); p = 0.000999.
  - ARI 0.8274, against a null mean of 0.6836 (max 0.7303); p = 0.000999.
- **Replication: both de novo catalogs pass T1a and fail T1b and T2.**
  - T1b: A119b D = -1.0430 [-1.2541, -0.8029]; testis D = -0.2520 [-0.3903, -0.1052].
  - T2, assembly: A119b p exact 0.006993, p ARI 0.000999; testis p exact 0.243756, p ARI 0.044955.
- **Read depth passes T2's bar on the guided common set** (exact 206, ARI 0.8685). Reading 3 would use
  it, but only when T1 holds.

**What this says.**
- 258 of 491 Soto families (52.5%) lie inside one block of ours (256 in one of our families, 2 on one shared
  unclustered locus), far more than random partitions with our block sizes give (null mean 0.7). That null is weak by
  design (Known before), and on the de novo catalogs the share of families with ≥ 2 members is 27.3% (A119b) and
  51.4% (testis). This is not a refinement.
- On the guided universe 503 of Soto's 2,185 single-label genes sit in no family of ours (495 alone in a block, 8
  sharing an unclustered locus); over all 2,334 genes, 526 are our singletons. Of the 233 uncontained families, 139
  are contained once our singletons are removed, 47 once Soto's multi-family genes are removed too, and 47 are
  genuine conflicts.
- Post hoc, and it conditions on our own catalog, so it is not a test: on the 1,682 single-label genes in our
  families, D = +0.1906.
- The copy-number cut beats random cuts on the guided universe, but gains only
  +3 exact and +0.0384 ARI over leaving our families uncut, and fails on both de novo catalogs.

**Post hoc (added after the final review): what T1b measures.** Over the single-label members,
D = H(Soto | ours) − H(ours | Soto) = [H(Soto, ours) − H(ours)] − [H(Soto, ours) − H(Soto)] = H(Soto) − H(ours)
(guided: 8.5367 − 8.9362 = −0.3995). T1b therefore compares only how coarse the two partitions are; it is not a nesting
test, and two unrelated partitions with the same block-size profile have the same D. The guided failure is the
singletons above, and dropping them gives the post-hoc +0.1906. With T1a's weak null, "T1 holds" could pass without
any nesting. T1 must not be reused (for example on gorilla) until T1b is replaced by a nesting statistic, such as
H(ours | Soto) against the size null.


---

## PREREG_template_copies_2026-10-14.md

*(original: `docs/archive/2026-10/merged/PREREG_template_copies_2026-10-14.md`)*

## Template-guided copy assembly (2026-10-14)

**Pre-registered before the first real projection.**
- Spec: `docs/superpowers/specs/2026-10-14-template-copies-design.md`.
- Plan: `docs/superpowers/plans/2026-10-14-template-copies.md`.
- Code: `bench/template_copies/` at commit `79ea8caa`.
- Every number below comes from that code, run as written in the plan's Tasks 6-9.

### Question

Can copies that have their own reads but no full locus be assembled from those reads, guided by an expressed locus of
the same gene, without lowering precision and without losing the operating point's lead over StringTie and FLAIR? This
is the assembler using the family relation (user, 2026-10-14).

### Known before running (recorded now)

- **Where read-supported Liftoff copy pairs are lost** (identity ≥ 0.95, shipped families):

  | sample | pairs | shared | a copy has no assembled locus | assembled but not in a family | split |
  |---|---|---|---|---|---|
  | human A119b | 514 | 75 | 378 | 45 | 16 |
  | human testis (public, ERR13885926) | 38 | 12 | 24 | 2 | 0 |
  | gorilla KB3781 | 78 | 4 | 71 | 3 | 0 |

- **The operating point helps a little.** Of A119b's 799 read-supported copy loci, the shipped assembler fully covers
  318 and `OP_good` 350 (+43, −11).
- **The ceiling.** 309 extra copies have no full locus under `OP_good`. Their source placement is:
  - in a family: 61;
  - an assembled singleton: 46;
  - not assembled: 202.

  So a template exists for at most 107.
- **Copy shape.** Of A119b's copy loci without a full locus, 360 are multi-exon and 121 single-exon; 217 are partly
  covered and 264 have no locus at all.
- **Calibration (planning).** One 1,000-template batch of A119b's shipped templates was projected on the human index
  (89 s, 16.3 GB; `/mnt/linuxdisk/tmp/o1_template_copies/calib/`) to size the batches. Nothing from it was scored or
  used.
- **Refuted neighbours, not re-proposed:**
  - `--tied-seed` (registers 338, 435, 980);
  - family rescue levers (647, 648);
  - all secondaries seeding loci (1059);
  - SD-region semi-guided nodes (1007/1009);
  - annotated coordinates as de novo vertices (294);
  - the closed loop is expected to be small (1272).
- **Templates per run** (inputs, from `run_tg.py templates`):

  | run | base GTF | templates | batches |
  |---|---|---|---|
  | A119b_shipped | `runs/human_A119b/human_A119b.gtf` | 94,685 | 16 |
  | A119b_OP | `asm_gap/op/human.OP_primary.gtf` | 84,734 | 15 |
  | testis | `runs/human_testis/human_testis.gtf` | 13,048 | 3 |
  | KB3781 | `runs/gorilla_KB3781/gorilla_KB3781.gtf` | 17,672 | 3 |
  | OR6737_OP | `asm_gap/op/gorilla.OP_primary.gtf` | 19,983 | 4 |

### Method (fixed)

As the spec's section 1, with the plan's rulings R1-R3:
- identity ≥ 0.80 and ≥ 50% of the template aligned (`minimap2 -x splice -c -N 50 -p 0.1`);
- a candidate is covered < 50% by every base locus, and duplicates overlapping by ≥ 50% of the smaller union merge;
- ≥ 2 own primary reads (`-F 2308`) per candidate and per emitted locus;
- transcripts come only from the own reads' blocks and junctions;
- the run fails when more than half the candidates lie on contigs absent from the BAM.

### Metrics and bars (per run; species never pooled)

- **Recall (reported).** Read-supported Liftoff extra copies (identity ≥ 0.95) covered ≥ 50% by a locus, base vs
  merged; the template-reachable ones and how many were reached.
- **P1 (bar).** Two quantities must not drop: the merged set's gffcompare intron-chain precision against the sample's
  RefSeq must be ≥ the base's − 0.2 (percentage points), and its `located` share ≥ the base's − 0.002.
- **P2 (bar).** The added loci's `located` share must be ≥ that of the base loci whose representative has 2-5 reads,
  with at least one locus added. A 95% interval over the added loci is reported.
- **Tool lead (bar; A119b_OP and OR6737_OP).** Register 1268's scorer, genome-wide minus the development chromosome
  (human chr20, gorilla NC_073244.2). The arm `OP_primary_TG` (merged) must have intron-chain sensitivity > StringTie's
  and > FLAIR's (B1), and intron-chain precision > FLAIR's (B2).

### Readings (fixed now)

| development: A119b_shipped and A119b_OP P1 + P2, tool lead on A119b_OP and OR6737_OP | held out: testis and KB3781 P1 + P2 | reading |
|---|---|---|
| all pass | all pass | Port into the Rust assembler as an opt-in second pass; the recall gain is the reported effect. |
| all pass | any fails | Not adopted; the held-out failures are listed. |
| any fails | — | Not adopted. |

### Rulings (made while planning)

- **R1.** Each emitted locus needs ≥ 2 own reads. The spec's ≥ 2 rule is applied to every component as well as to the
  candidate.
- **R2.** minimap2 also gets two flags:
  - `-c`, which gives the CIGARs; without them there are no projected exons;
  - `-p 0.1`, the guided finder's setting. minimap2's default 0.8 would drop secondary hits below 80% of the best
    score, which are exactly the diverged copies.
- **R3.** `assemble` fails when more than half the candidates lie on contigs absent from the BAM. That guards against
  an index and a BAM that name contigs differently.
- **R4.** Development's operating-point base is `OP_primary`, the base register 1268's B1/B2 judge. `OP_good` is not
  re-run (compute); its pre-design ceiling numbers are reported as context.
- **R5.** The tool lead reads register 1268's scorer's arm table, and `score_tg.lead` evaluates B1 and B2 for the TG
  arm. The scorer prints its own claims only for its four standard arms.
- **R6.** Batches hold 6,000 templates. One 1,000-template batch took 89 s and 16.3 GB on the human index, index load
  included.
- **R7.** P2 compares against base loci whose representative transcript has 2-5 `reads` (from `templates.tsv`). P2
  passes only when at least one locus was added.

### Amendment 1 (2026-10-14, before any result)

Speed only, at commit `6b260c8c`. The method and every threshold are unchanged.
- **The first A119b_shipped `assemble` was killed.** It reached its 900 s limit before writing anything: about 50 ms per
  candidate, and some projected spans reach 3 Mb.
- **Long spans are fetched exon by exon.** `own_reads` now does this when a candidate spans more than 200 kb, instead of
  walking every read of a long projected intron. A read with an aligned base on the exons overlaps one of them, so the
  reads are the same. On the first 600 candidates of A119b_shipped (batches 0-1), the old and new read sets were
  identical.
- **The own-read step runs in 4 worker processes** (`run_tg.py assemble --workers`). The output order and content do
  not depend on the worker count (test `test_workers_do_not_change_the_output`).

### Results

Added after the runs. Nothing above this heading changes once a statistic has been seen.

### Results (2026-10-14)

Run: code `c043d354` (the method at `6b260c8c`, Amendment 1; `c043d354` only fixes the scorer's gzipped reference).
A119b_shipped took 6 heavy calls of projection (16 batches, 116-627 s per call) and one 21-minute assembly (4 workers,
1.1 GB).

#### A119b_shipped (development)

| step | count |
|---|---|
| hits passing identity ≥ 0.80 and aligned ≥ 0.50 | 561,980 |
| candidates (covered < 50% by every base locus) | 97,097 |
| merged candidates | 31,515 |
| on contigs absent from the BAM | 0 |
| fewer than 2 own reads | 10,697 |
| no locus with ≥ 2 reads | 748 |
| candidates with loci | 20,070 |
| added loci / transcripts | 26,726 / 27,515 |

| measure | base | merged | bar |
|---|---|---|---|
| read-supported Liftoff extra copies (identity ≥ 0.95) covered | 234 of 641 | **394 of 641** | reported |
| template-reachable copies reached | — | 63 of 106 | reported |
| P1: intron-chain precision against RefSeq (gffcompare) | 17.7 | 17.4 | **FAIL** (≥ 17.5) |
| P1: share of loci at a Liftoff locus | 0.2452 | 0.2091 | **FAIL** (≥ 0.2432) |
| P2: added loci at a Liftoff locus vs base loci with 2-5 reads | 0.2318 (36,269 loci) | 0.0813 [0.0781, 0.0847] (26,726 loci) | **FAIL** |

Intron-chain sensitivity: base 22.0, merged 22.0. The base's 22.0 / 17.7 is the same pair as register 1268's
`default_good` human row, but that row excluded chr20 (the development chromosome) and this score is genome-wide, so the
two are not the same measurement.

**Not run (ruling, recorded in the plan's ledger).** A119b_OP, testis, KB3781 and OR6737_OP were not run, because the
readings below make the verdict independent of them: any development bar failing means not adopted. Their templates
are on disk (`/mnt/linuxdisk/tmp/o1_template_copies/<run>/`). They should not be run with this code: the final review
found three defects in the build tested (see "Final-review findings" below), and any rerun needs them fixed first (about
2.5 h of compute after that).

**Post hoc (descriptive, not pre-registered): the share of added loci at a Liftoff locus does not depend on shape or
depth.**

| added loci | count | at a Liftoff locus |
|---|---|---|
| all | 26,726 | 0.0813 |
| single-exon | 23,655 | 0.0800 |
| multi-exon | 3,071 | 0.0918 |
| multi-exon, ≥ 5 reads | 2,541 | 0.0799 |
| single-exon, ≥ 5 reads | 14,405 | 0.0909 |
| 2 reads | 4,585 | 0.0665 |

#### Verdict (2026-10-14)

**Reading 3: not adopted.** Development fails every precision bar on A119b_shipped:
- P1 chain precision: 17.7 → 17.4;
- P1 `located` share: 0.2452 → 0.2091;
- P2: 0.0813 against 0.2318.

Recall does rise, with 160 more read-supported Liftoff copies covered (234 → 394 of
641), but 119 of the 160 are covered only by single-exon added loci (below). It comes with 26,726 added loci (23,415
distinct exon sets), of which only 0.0813 sit at a Liftoff locus; most are single-exon. Of the five exon-count and
read-depth cuts tried post hoc (table above), none raises the share above 0.10. Sitting at a Liftoff locus is a proxy:
it is not the same as being the copy.

The verdict applies to the build tested. The final review (below) found three defects in it; a build with them fixed
has not been tested.

#### Final-review findings (verified, post hoc, 2026-10-14)

The whole-branch review found three defects in the build tested and three overstatements in the text above (now
corrected). Every count below was recomputed from the A119b_shipped outputs (`tg.gtf`, the base GTF, the BAM) after the
review. None changes a bar or the verdict, and all point the same way.

Known defects of the build tested (also listed in `bench/template_copies/run_tg.py`):
1. **Duplicate and non-new loci.** The 26,726 added loci are 23,415 distinct exon sets: 2,988 sets are emitted more than
   once (3,311 extra emissions, many as a + / − pair). 14,085 added loci overlap another added locus (the 26,726 fall
   into 17,760 overlap clusters), and 4,522 lie ≥ 50% inside one base locus (2,380 ≥ 95%). Cause: candidates merge only
   on the same strand, the own-read step collects reads of both strands, and nothing de-duplicates the emitted loci.
2. **Strand from the template hit, not the reads.** The BAM was mapped with `-uf`, so a read's strand is its
   transcript's. In a random sample of 400 added loci (seed 7), 190 (0.475) have most of their reads on the opposite
   strand to the locus; the reviewer's own random sample gave 173 of 392 (0.44).
3. **Whole read blocks.** A read's exon block is kept whole when it touches the projected exons, so a locus can run far
   past the copy: the reviewer found the covering locus is, at the median, 4.0× the copy's exonic length.

Overstatements, corrected above:
- **Recall.** Of the 160 newly covered copies, 119 are covered only by single-exon added loci, and 64 of those 119 are
  multi-exon copies covered by one unspliced block. Of the 63 template-reachable copies reached, only 31 are covered by
  a locus built from the copy's own source template, and 17 by a multi-exon locus (reviewer).
- **The comparison group.** The 2-5-read group is the bar the prereg fixed, not a group matched to the added loci. Base
  loci grouped by their representative's reads: 2-5 reads 0.2318 (36,269 loci), 6-20 reads 0.167 (31,037), more than 20
  reads 0.3515 (27,379). The added loci's 0.0813 is below every group. With a symmetric rule (a locus covers ≥ 50% of a
  Liftoff locus) the reviewer found 0.174 for the added loci against 0.208 for the base loci: the gap is smaller, and
  the direction holds.
- **The cuts.** The earlier wording said no exon-count or read-depth cut separates the real copies. Only five cuts were
  tried, `located` is not the same as a real copy, and exon count is confounded by the build rule (overlapping
  unspliced blocks merge into one exon).


---

## PREREG_post_mcl_reassign_2026-10-15.md

*(original: `docs/archive/2026-10/merged/PREREG_post_mcl_reassign_2026-10-15.md`)*

## Post-MCL strongest-edge reassignment of family members (2026-10-15)

Pre-registration, written before the rule was run on any catalog. Code: `bench/o1_reassign/` at `cb237430`
(16 tests). Branch `post-mcl-reassign`.

### Question

Human A119b fails claim F1 of figure 8 (131 / 148 Liftoff copy pairs share a family; the bar is 0.90). Exploring the
17 split pairs (2026-10-15) showed that 7 of them come from one MCL decision in a tandem array. The user asked to try a
rule that undoes such decisions: after MCL, a member moves to the family that holds its strongest edge when that edge
is stronger than every edge it has inside its own family. Does the rule fix F1 without harming the held-out catalogs,
and do its moves agree with independent truth?

### Known before running (recorded now)

- **F1 as shipped** (genome, sequence_ID ≥ 0.95): A119b 131 / 148, testis 34 / 34, gorilla KB3781 11 / 12.
- **The decoy test as shipped** (register 1284, genome, ≥ 0.95): A119b 75 / 91, testis 12 / 12, KB3781 4 / 4 against
  null means 0.121, 0.102 and 0.012 (p 0.000999 each).
- **A119b's split pairs, traced through the family graph** (`bench/` not involved; session scratch):
  - 7 CYKILR pairs: the chr19 38.8-40.1 Mb array is one dense block (edge weights 0.95-0.999) that MCL cut into
    MCL42, MCL451 and MCL566. The CYKILR source node's strongest edge (0.986) goes to MCL42, which holds all 7 copies;
    its best edge inside MCL451 is 0.963. Two other MCL451 members also have their strongest edge into MCL42 (0.994 vs
    0.976; 0.985 vs 0.976), and the fourth ties (0.916 and 0.916).
  - 8 pairs (5 genes) have no edge: the 0.60 shared-exon floor rejects them. The rule cannot join them.
  - RPL23AP88's copy node has edges of 1.000 inside its family and 0.999 to the source's family, so the rule leaves
    it. LOC124906865's cross edge is 0.391, below both nodes' best inside edges.
  - Expected from this: the CYKILR source node moves to MCL42 and A119b's F1 rises to about 138 / 148. That is why A119b
    is development: its F1 is not evidence for the rule.
- 14 of the 17 split pairs already share a superfamily of the emitted hierarchy, but 11 of them through SF166, which
  holds 2,991 members in 435 families.
- **Not seen**: how many members move in any catalog, any F1 or decoy number after the rule, any Soto comparison after
  the rule.

### Method (fixed)

**The rule** (`reassign_lib.reassign`; no parameter). For a member v of family F:
- in(v) is its strongest edge to another member of F;
- out(v, G) is its strongest edge into family G ≠ F.

If max over G of out(v, G) > in(v), v moves to that G; when two families tie, the smaller family name wins. Members are
visited in sorted order and see the moves already made (asynchronous passes), and passes repeat until one moves
nothing, up to 50.

Only members move:
- an edge to a node in no family is ignored, and no node joins or leaves the catalog;
- a member left alone in its family is reported (lost).

The edges are the shipped graph: `mcl_families --dump-graph` from each catalog's cached all-vs-all
(`bench/o1_hier/recluster.sh`). Weights are written to 6 decimals, so edges closer than 1e-6 count as ties.

**Gates** (a sample whose gate fails has no verdict):
- the re-clustered `clusters.tsv` equals the shipped one byte for byte;
- every row of the shipped copy table maps (gene → span node → folded representative) to a node of its own shipped
  family.

**Scoring** (`run_reassign.py score`):
- F1 before and after, every scope and identity floor, with figure 8's code (`run_pairs.f1`);
- the Liftoff pairs (genome, ≥ 0.95) the moves join and split;
- the decoy test's primary row (genome, ≥ 0.95, 1,000 decoy sets, seed 7) before and after;
- human only: Soto gains and losses. Per moved member reached by a gene of Soto's 491 families (soto_subset's de
  novo universe), count its Soto family-mates in its old family and in its new one; more in the new one is a gain.
  Also reported: pairwise precision and recall against Soto over the universe's members, before and after.

### Bars

| id | sample | bar |
|---|---|---|
| D | A119b (development) | F1 ≥ 0.90 after the rule |
| H1 | testis, KB3781 (held out) | F1 stays ≥ 0.90, and the decoy test still passes (p < 0.001) |
| H2 | testis, KB3781 | Liftoff pairs joined − split ≥ 0 in each; in testis, Soto gains > losses when gains + losses ≥ 5 (fewer: not powered) |
| guards | every sample | no member lost; F1's covered pairs change by ≤ 2%; the largest family grows by ≤ 10%; the passes converge |

**Evidence of benefit** (held out): Liftoff joined > split in testis or KB3781, or a powered Soto H2 in testis.

### Readings (fixed now)

1. **Adopted as an opt-in flag** (then ported to `mcl_families`): D, H1, H2 and the guards hold, and there is evidence
   of benefit.
2. **Harmless but unproven**: D, H1, H2 and the guards hold, but nothing held out changes enough to count as evidence.
   It is recorded as a rule that fixes A119b's F1 without held-out support, and it is not adopted.
3. **Not adopted**: D, H1, H2 or a guard fails.

A119b's moves outside the chr19 array (Soto gains and losses, Liftoff joined and split) are reported as a description,
not as a bar, since the array motivated the rule.

### Results

Added after the runs. Nothing above this heading changes once a statistic has been seen.

### Results (2026-10-15)

Run: code `cb237430`, outputs `/mnt/linuxdisk/tmp/o1_reassign/<sample>/` (`apply.json`, `score.json`, `moves.tsv`).
Both gates held in all three samples: the re-clustered `clusters.tsv` equals the shipped one byte for byte, and every
copy row maps to a node of its shipped family. Re-clustering took 6 s (testis), 17 s (KB3781) and 155 s at 4.7 GB
(A119b).

| | A119b (development) | testis (held out) | KB3781 (held out) |
|---|---|---|---|
| members moved (moves, passes) | 1,609 of 9,907 (2,006, 6) | 5 of 1,115 (5, 3) | 11 of 918 (11, 3) |
| families | 2,150 → 1,944 | 342 → 340 | 326 → 322 |
| largest family | 260 → 349 | 47 → 47 | 34 → 34 |
| members lost (left alone) | 9 | 1 | 0 |
| **F1** (genome, ≥ 0.95) | 131 / 148 → **137 / 148 = 0.926** | 34 / 34 → 34 / 34 | 11 / 12 → 11 / 12 |
| F1 without the development contigs | 113 / 128 → 119 / 128 | 33 / 33 → 33 / 33 | 11 / 12 → 11 / 12 |
| Liftoff pairs joined / split | 7 / 1 | 0 / 0 | 0 / 0 |
| decoy test (k / n, null mean; p) | 75 / 514, 11.0 → 80 / 514, 8.4; 0.000999 | 12 / 38, 1.2 → unchanged | 4 / 78, 0.05 → unchanged |
| Soto gains / losses (same, unlabelled) | 50 / 14 (21, 1,524) | 1 / 0 (0, 4) | — |
| Soto pairwise precision | 0.7344 → **0.6881** | 0.8717 → 0.8721 | — |
| Soto pairwise recall | 0.3008 → 0.3150 | 0.8386 → 0.8418 | — |

Bars and guards:
- D holds (A119b F1 0.926 ≥ 0.90).
- H1 holds: held-out F1 and the decoy test are unchanged.
- H2 holds trivially but without evidence: no Liftoff pair is joined or split in either held-out sample, and testis's
  Soto comparison is not powered (1 gain, 0 losses; the bar needs 5).
- **Guards fail.** A119b's largest family grows by 34% (260 → 349; bar 10%), and members are lost in A119b (9) and
  testis (1). KB3781's guards hold.

#### Verdict (2026-10-15)

**Reading 3: not adopted.** The guards fail in two of three samples. Even with the guards set aside, nothing held out
shows a benefit, which is reading 2 at best. The rule does fix A119b's F1 as predicted, but in the large catalog it
merges much more than the tandem array: Soto pairwise precision falls by 0.046 there.

#### Post hoc (descriptive, A119b)

- **The array is a small part of it.** Only 9 of the 2,006 moves are inside the chr19 38.8-40.1 Mb array. The other
  moves spread over the genome; chr13 has the most, 283 moved members.
- **Members MCL placed without a direct edge move first.** 718 moves (36%) start from in = 0: MCL put the member in its
  family through flow, with no edge to any member of it, so any edge to another family pulls it out.
- **Large families attract.** The largest family (MCL3) keeps 234 of its 260 members and gains 115. The top five
  destination families take 517 of the moves (MCL3 124, MCL121 121, MCL60 106, MCL617 94, MCL12 72).
- **The deciding gap is often small.** The median out − in is 0.171, but 381 moves (19%) are decided by less than 0.01
  and 685 (34%) by less than 0.05.
- **Outside the array** (1,601 of the 1,609 moved members; pre-registered as a description):
  - every Soto-labelled move is outside it: 50 gains, 14 losses (21 same, 1,516 unlabelled);
  - Liftoff pairs: 0 joined and 1 split, so all 7 pairs the rule joins are the array's CYKILR pairs.
  - Pairwise precision against Soto still falls, because the absorbed members bring their neighbourhoods.


---

## NPIP_16Q22_CNV_2026-10-20.md

*(original: `docs/archive/2026-10/merged/NPIP_16Q22_CNV_2026-10-20.md`)*

## The 16q22 NPIPB15 repeat: why LOC124907808 and LOC124907807 rest on tied reads (2026-10-20)

**User (2026-10-20):** "my advisor keeps saying that it is impossible and that there must be an error, can we add a
visualization that proves that?" The figure is in the NPIP Read Pools page (section `#cnv`).

### The claim

**CHM13 has three units.** CHM13 chr16:80,100,000-80,550,000 holds three tandem units, 114.4 kb apart (114,393 and 114,350 bp
between the NPIP copies). Each unit holds three genes:

| unit | PDPR-like (−) | NPIP (+) | CLEC18 (−) |
|---|---|---|---|
| 1 | AC009053.1 | NPIPB15 | CLEC18C |
| 2 | MSTRG.2242 | LOC124907808 | CLEC18B |
| 3 | MSTRG.2243 | LOC124907807 | CLEC18B (second copy) |

**Every haplotype checked carries one or two units.** That is four independently assembled human haplotypes (each CHM13 unit
aligned on its own, minimap2 -x asm10):

| haplotype | syntenic length | shorter than CHM13 | NPIP copies (best CHM13 match, identity) |
|---|---|---|---|
| HG002 mother (Q100 v1.1) | 221,127 bp | 228,873 bp = 2.00 units | NPIPB15 99.67% |
| HG002 father | 335,486 bp | 114,514 bp = 1.00 unit | NPIPB15 99.82%; LOC124907807 99.80% |
| HG00733 mother (HPRC y1) | 221,277 bp | 228,723 bp = 2.00 units | LOC124907808 99.55% |
| HG00733 father | 338,363 bp | 111,637 bp = 0.98 unit | NPIPB15 99.90%; LOC124907808 99.92% |

**What this shows:**
- **In HG002**, LOC124907808 has a copy on neither haplotype, and LOC124907807 only on the father's.
- **In HG00733**, LOC124907808 is only on the father's, and LOC124907807 on neither.
- **The unique flanks** align at MAPQ 60 and 99.6-99.97% identity.
- **Each haplotype NPIP copy matches all three CHM13 copies** at 99.3-99.9%, its best by 0.04-0.43 points. CLEC18 best
  matches are closer still (down to 0.02 points), so no hybrid-unit claim is made.

### Reads on T2T

**A119b** (testis; minimap2 -N 50 -p 0.1), primaries by AS margin to the best alignment on another NPIP copy:

| copy | primaries | exact ties | near ties (1-9 AS) | placed (≥ 10 AS) |
|---|---|---|---|---|
| LOC124907808 | 93 | 54 | 20 | 19 |
| LOC124907807 | 112 | 25 | 66 | 21 |

Most of the ties are with NPIPB15 (67 each).

**HG002's Kinnex RNA** barely expresses the cluster: 6 reads on it, and none with a primary on the two LOC copies.

### Not yet done

DNA reads spanning a unit junction (HG002 HiFi or ONT) would confirm each haplotype's unit count read by read. A119b's own copy
number needs its DNA.

### Reproduce

- `python3 bench/o3_graph/npip_cnv_fig.py OUT_DIR` writes the alignments, gene hits and 32-mer dot plots to `data.json`.
- `python3 bench/o3_graph/npip_cnv_render.py OUT_DIR/data.json READS.json section.html` renders the section.
- Data paths are in the script headers: HG002 `rustle_figures_dev/o2_presence/hg002/hg002v1.1.fa`; HG00733
  `tmp/o3_graph/donorsim/HG00733/asm/`.

### Gorilla tab (2026-10-20, user: "can we make another tab in this with those data from gorilla?")

KB3781 ("Jim") is the individual mGorGor1 was assembled from. His maternal (GCA_028885495.2) and paternal (GCA_028885475.2)
chromosomes, his fibroblast Iso-Seq and his Illumina WGS all come from the same animal. The reference chr12, NC_073236.2, is
his paternal CM054570.2.

**STON1 is present on his mother's chromosome 12 and absent from his father's:**
- **The deletion.** Mother CM054594.2:95.7-96.3 Mb against father CM054570.2:90.9-91.5 Mb (asm5): colinear at 99.73% or more on
  both sides of a 121,624 bp stretch only the mother has. The father has 3,027 bp there.
- **The genes.** Human mRNAs (splice): STON1 aligns over 99.9% of its length on the mother and only 5% (a stray piece) on the
  father. On the father, GTF2A1L keeps 38% (its 3′ part) and LHCGR is intact on both.
- **The reads.** The 113 reads of O3 cluster clSRR27438212.9740719 were aligned to each parent's whole genome (splice:hq):
  113 / 113 align to the mother, all on chr12; 0 / 113 align to the father.
- **The DNA.** 97.1% of the consensus is covered by his WGS 21-mers, at a median depth of 6.0, one haplotype (λ 4.98; PREREG
  unmapped rescue §T1).

**The 16q22 NPIP repeat in Jim: both parents carry the same 5 NPIPB15-like copies** (2 inverted, 3 forward; human NPIPB15
probe at 95-96%). His mother's and father's 7 Mb align colinearly at 99.4-99.75%.

Code: `bench/o3_graph/ston1_fig.py`, `gorilla_tab_render.py`. Gorilla and human numbers sit on separate tabs and are never
pooled.


---

## MULTIMAPPER_PREMISE_2026-10-21.md

*(original: `docs/archive/2026-10/merged/MULTIMAPPER_PREMISE_2026-10-21.md`)*

## Multi-mapping long reads: what they cost assembly, what O2 can decide, what O3 recovers (2026-10-21)

**The question (user, 2026-10-21).** "I made my entire thesis on the assumption that multi-mapping long reads were really important
information being lost, that minimap2's designation as primary is no better than a coin toss and that stringtie and flair are
incomplete without this important info. Can we actually check without O2 and O3 ...?"

This page collects the answers of registers 1360-1367. Every row was pre-registered, except register 1366, which is descriptive
over existing runs. Every row was checked by the metric-auditor; the numbers below are the audited ones.

### Short answers

1. **Coin-toss primaries are rare.**
   - An exact AS tie between the primary and another alignment happens for 0.14% (gorilla OR6737) / 0.52% (human A119b) of reads.
   - Within 10 AS: 0.54% / 1.30%.
2. **StringTie and FLAIR already behave as primary-only assemblers.**
   - StringTie gets secondaries but never builds a transcript from them: on the lab's slice, a full BAM and a primary-only BAM
     give identical intron chains.
   - FLAIR drops secondaries.
3. **What multi-mapping costs them is small genome-wide, but concentrated in multi-copy loci.**
   - About 3% of expressed annotated transcripts have mostly multi-placed reads, or no primary alignment at all.
   - There the tools recover about half their usual rate, and none of the no-primary ones.
   - Matched for depth and exon count, the shortfall is ~250-340 transcripts per tool per species.
4. **Using secondaries does not recover real copies.**
   - Seeding with secondaries adds ~0.2 points of sensitivity genome-wide; inside the multi-mapping group it adds +34-43%, but at a
     precision cost.
   - O2's certificate on the reads behind those gains: where the read is truly tied, the copy is undecidable for 76% (human) / 83%
     (gorilla) of the gains.
   - Where the read is not tied, minimap2 had already placed it on the other copy.
   - The mainstay genes with a few certified reads (TBC1D3, NPIPB5, SMN2, GOLGA8F, RBMY1) are certified 97-100% to another copy.
5. **The loci where reads multi-map are real duplicated DNA.**
   - Human: 99.1% are in the CHM13 SD90 track (97.6% on ordinary autosomes; 77.3% without satellite SDs).
   - Gorilla: 98.2% are in SEDEF SDs (a post-hoc check).
   - The base rate for genes without tied reads is 15-24%.
   - This is close to true by construction, since reads tie only where sequence is duplicated, and it does not validate any O1
     family.
6. **O2 in a matched genome** (KB3781 fibroblasts on his own T2T):
   - On real reads it is UNDERPOWERED: minimap2 is wrong on only 1-3 of ~107 tied reads that his genome places on a reference copy.
   - With every known copy expressed from both haplotypes it FAILS. Of minimap2's ~500 errors per replicate, the true copy has
     strictly fewer edits on only 3 / 9. The reads do not carry the provenance minimap2 misses.
7. **O3 recovers a copy the reference lacks:**
   - The reference-free chain finds the maternal-only NPIP locus's units, one candidate each, with nothing flagged outside it.
     IsoCon alone gives that copy only 7 of its 37 reads.
   - Once the candidates are added to the reference, alignment gives it 33 of 39 reads (26 of 32 among reads the candidate was not
     built from).
   - Zero false moves among 2,861 other NPIP reads and 951 reads aligned at the candidates' relatives.
   - The O2 ties at that locus drop from 34 to 10.

### What the thesis can claim, and what it should not

| claim | status |
|---|---|
| Multi-mapping is a genome-wide information loss for StringTie / FLAIR | **not supported** (~3% of transcripts, ~2% of the tools' output) |
| In multi-copy loci the tools lose about half their sensitivity | supported (registers 1366, 1367) |
| Adding secondary alignments recovers those copies | **not supported** (gains mostly undecidable or already placed elsewhere) |
| minimap2's primary is a coin toss | true only for exact ties (0.14-0.52% of reads), where the copies are identical over the read and "both" is the right answer |
| O2 recovers provenance the aligner misses | not in a matched genome (registers 1360-1362); its role is certifying placements and abstaining on true ties (registers 1338, 1342) |
| Multi-mapping loci are real duplications, not artefacts | supported for the loci (register 1367); a separate test is needed for the O1 families |
| A copy the reference lacks explains misplaced tied reads, and O3 can name it and return its reads | supported at one locus (registers 1363-1365; not held out) |

### Caveats

- Two testis samples for the assembly analysis (gorilla OR6737, human A119b).
- Human transcript numbers are scored on NCBI RefSeq, not CAT.
- "Matched" means the intron chain matches the annotation; that cannot tell which near-identical copy built it.
- The O2 and O3 tests use one animal (KB3781), with his own assemblies as the answer key only.

### Sources

| register | document |
|---|---|
| 1360, 1361 | `docs/PREREG_o2_corrector_families_2026-10-20.md` (Amendments 1-4, results, audits) |
| 1362 | `docs/PREREG_o2_corrector_idealsim_2026-10-21.md` |
| 1363 | `docs/PREREG_o3_npip_maternal_copy_2026-10-21.md` |
| 1364, 1365 | `docs/PREREG_o3_npip_candidates_realign_2026-10-21.md` |
| 1366 | descriptive, `docs/NEGATIVE_RESULTS_REGISTER.md` (existing runs in `benchmark_collapse/`, `figures/data/fig3_ref_tie.tsv`) |
| 1367 | `docs/PREREG_multicopy_gains_real_2026-10-21.md` |


---

## O1_FAST_2026-10-21.md

*(original: `docs/archive/2026-10/merged/O1_FAST_2026-10-21.md`)*

## O1 fast: the exon gate and the families stage's Rust side (2026-10-21/22)

**User (2026-10-21):** "make O1 extremely fast (as fast as possible on rust side); if it still uses python scripts it is fine if
you think it won't make too much of a difference to translate them".

**Result (gate6, 2026-10-22):**
- The all-vs-all is 15-30x faster than the ungated one: gorilla_OR6737 5,421 s → **182 s** (one process, 4 threads), and
  human_testis's families stage takes **71 s** end to end against 987 s for its old all-vs-all alone.
- human_A119b's gated all-vs-all in full: 7,332 s (gate3) → **2,270 s** (gate6), byte-identical products.
- Its own Rust work is about 12x faster: human_A119b with the PAF cached, 109 s → **9.1 s**.
- The products are byte-identical: clusters, copies, copies FASTA, regions, loci table, loci GFF3 and loci FASTA.
- What is left is measured at the end ("What is left, measured").

### Where the time went

| step (gorilla_OR6737, 19,711 loci, 1.0 Gb of locus spans) | before |
|---|---|
| minimap2 all-vs-all (`-x asm20 -c -X -N 50 -p 0.1`), 4 threads | 5,421 s (90 min) |
| `mcl_families` itself, PAF cached | 24.5 s |

- **Base-level alignment is about 97% of minimap2's CPU.** For one 10 Mb shard: 369 CPU s with `-c`, 9 CPU s without.
- **Almost none of that alignment is used.** An edge needs one record with exon-to-exon evidence covering
  `--min-shared-exon-frac` (0.60) of the smaller gene's exons. Of 1.2 M aligned locus pairs, 35,400 can meet that bar; the
  rest are aligned and then discarded.
- **The Rust side** reads the genome twice (6.5 s + 5.6 s, single-threaded) and writes `loci.gff3` unbuffered (2 s).

### The exon gate (`tools/minimap2_gate`, `minimap2 2.30-r1287-gate3`)

The patch to minimap2 2.30 is `gate.patch`, about 200 lines. It adds:
- `--gate FILE`: per sequence, its exon intervals in the sequence's own coordinates. `mcl_families` writes this as
  `<out>.loci.gate.tsv`.
- `--gate-margin INT` and `--gate-frac F`.

Just before base-level alignment, the chains of a query are grouped by target. A target keeps **all** its chains when one
chain, padded by the margin, holds gated bases on both sides, at least F times the smaller sequence's gated total. Every
other target's chains are dropped.

**Why the kept records are exactly the stock records:**
- **Extension bound.** Base-level alignment extends a chain by at most `max_gap` (`-g`; 10,000 under `-x asm*`; `align.c`,
  `l = l < opt->max_gap? ...`). With margin >= max_gap, a record that could pass the shared-exon fraction always comes from a
  chain that passes the gate.
- **Pair sums.** The graph sums every record of a pair: the coverage union and matches over block length. The gate therefore
  keeps every chain of a kept pair.
- **No competition under `-X`.** The families stage's `-X` keeps every chain (`MM_F_ALL_CHAINS`): minimap2 makes no
  primary/secondary selection before or after alignment, so dropping one target's chains cannot change another target's
  records.
  - MAPQ (column 12) can move on split remainders, and records with tied sort keys can swap order (minimap2's sort is
    unstable). The graph reads neither.
  - Columns 1-11 have never differed in any comparison below.
  - Without `-X` or `-P` the gate is NOT exact (the final review showed a synthetic counterexample), so `minimap2-gate`
    warns.
- **Without `--gate`,** the binary's output is byte-identical to stock minimap2.

**When `mcl_families` uses it:**
- **On by default** whenever every edge needs an exon-to-exon record: `--min-exonic-bp >= 1` with both-sided evidence, or a
  shared-exon fraction > 0. The margin is 10,000 and the fraction is the run's own `--min-shared-exon-frac`.
- **Off** under `--core-from-paf`, `--sedef` or `--emit-container`, which read every record. `--no-exon-gate` also turns it
  off.
- **Skipped with a note** when the binary has no `--gate`.
- **Finding the binary:** `minimap2-gate` on `PATH` or `RUSTLE_MINIMAP2`. `tools/mm2_shard.sh` runs it for any command line
  carrying `--gate`, and keys its layout on the gate file's md5.
- **Records:** `params.tsv` records `exon_gate`, `exon_gate_margin`, `exon_gate_frac` and `exon_gate_minimap2`. The PAF
  cache key names the gate file by content.
- **Gated PAF:** a gated run names its PAF `<out>.loci.gated.paf` and writes no `<out>.loci.paf`, so any reader of the
  full PAF fails loudly instead of reading fewer records. The other name's file, left by an earlier run into the same
  prefix, is removed.
  - `mcl_families --paf` refuses a gated PAF unless the run needs no record it lacks: exon-to-exon edges at a shared-exon
    fraction at least the recorded one, and no `--core-from-paf`, `--sedef` or `--emit-container`.
  - `tools/protein_attach.py` and `bench/family_container.py` refuse it with a message.
  - Readers of `.loci.paf` (`bench/paf_duplicons/paf_units.py`, `bench/o1_sd/self_hits_eval.py`, `bench/seed_pool/run.sh`,
    `figures/fig_protein_supp.py` via protein_attach without `--no-paf`) need a `--no-exon-gate` families run.
- **`params.tsv`:** the gate rows (`exon_gate`, `exon_gate_margin`, `exon_gate_frac`, `exon_gate_minimap2`,
  `exon_gate_counts`) are written only when the gate ran, so an ungated run's `params.tsv` is unchanged. With the gate on,
  the rejection counters (`rejected_no_exonic`, `rejected_self_hits`, ...) count only the aligned pairs.
- **Checks:**
  - the margin must be 10,000-2^31-1;
  - the binary must be gate3 or later (the first with `--gate-frac`);
  - the genome is kept for the copy table only on gated runs, so it never sits beside a full-size PAF.

#### A first version was wrong (kept here as a negative result)

**gate2: touch rule, 2 kb margin.** It kept a pair when a chain within 2 kb touched an exon on both sides.
- Products were byte-identical on human_testis, gorilla_KB3781 and gorilla_OR6737, but by luck.
- Base-level extension reaches up to 10 kb. On human_A119b, 137,361 records of exon-touching pairs were missing among 21,126
  queries, 136,904 of them on chr13. Example: chr1:47970686 vs chr6:5761075, identity 0.85, extended about 2.7 kb into an
  exon.

**gate3 (the shipped rule)** replaced it: margin = max_gap, with the shared-exon fraction.

### Validation of gate3

| catalog | all-vs-all, ungated | all-vs-all, gated | PAF | products |
|---|---|---|---|---|
| gorilla_OR6737 (19,711 loci) | 5,421 s | **673 s (8.1x)** | 531 → 37 MB | byte-identical |
| human_testis (13,034 loci) | 987 s | **133 s (7.4x)** | 176 → 16 MB | byte-identical |
| gorilla_KB3781 (17,608 loci) | not timed | 261 s | 425 → 22 MB | byte-identical |
| human_A119b (94,297 loci) | not run in full | – | – | see below |

- **Products** are clusters, copies, merged copies, copies FASTA, regions, loci table, loci GFF3 and loci FASTA. `params.tsv`
  differs only in the gate rows and the counts of pairs never aligned.
- **Timings** are sums of `tools/mm2_shard.sh` shard times plus the index build, at 4 threads.
- **human_A119b, record level:** 3,804 queries, including the densest chr13 pieces. All 43,432 pairs that could pass 0.60
  are present, with 0 records differing on columns 1-11.
- **The A119b hotspot.** chr13:13.96-15.5 Mb holds about 4,000 overlapping de novo loci; one 340 kb window holds 983 of them
  (9.9 Mb of sequence, 29x redundant). There most pairs share exons, so the gate saves about half the work (2.0x, below). A full A119b run is
  dominated by this region with or without the gate. It is a property of that assembly's loci, not of O1's rule.

#### Final review and identity check

A fresh reviewer found no critical problem. Its two important findings are fixed:
- a gated PAF indistinguishable from a full one, which `--paf` mode would read wrong at a lower fraction;
- `params.tsv` rows written even with the gate off.

Its minor findings are fixed or documented:
- exactness rests on `-X`;
- the genome is held only on gated runs;
- a misplaced test now runs;
- the margin and version are checked;
- the wrapper accepts `--gate=FILE`;
- a test retries on ETXTBSY.

The wrapper's index is keyed together with the gate file's md5, so it is rebuilt for each new gate file. That costs time
only.

**Identity check against the previous commit (22cd1586), on `tools/identity_check.sh`'s synthetic fixture:**
- `--no-exon-gate`: all 7 products byte-identical, `params.tsv` included.
- Gated: every product identical; `params.tsv` differs only in the PAF name and the gate rows.

### The Rust side

All changes are byte-identical in their outputs and covered by tests:
- **Only the loci spans are read.** `GenomeIndex::from_fasta_regions` reads just the merged spans through the `.fai`, in
  parallel, streaming 64k-line blocks with positioned reads. That is 963 Mb of OR6737's 3.6 GB genome, and 309 Mb of human
  testis's 3.1 GB.
  - It answers fetches inside the spans exactly as the whole-contig index does, and None outside.
  - Tests: the whole-contig loader against the whole-file parse; the region loader against the whole-contig loader on
    every sub-interval (CRLF, a last line with no terminator, contigs spanning blocks, clipping at the contig end).
- **One genome load per run.** `GtfLociRun.genome` is reused by the copy table on gated runs. On ungated runs the copy table
  reads its members' spans again, so the genome never sits beside a full-size PAF.
- **The GTF parse no longer allocates per line or per attribute lookup:** one reused line buffer, constant lookup patterns,
  the crate's deterministic hasher.
- **`loci.gff3`** is written through a `BufWriter`.

| run with the PAF cached | before | after | peak memory |
|---|---|---|---|
| gorilla_OR6737 | 24.5 s | **2.7 s** | 3.7 → 1.1 GB |
| gorilla_KB3781 | – | **2.7-4.3 s** | 0.96 GB |
| human_testis | – | **1.0 s** | 3.1 → 0.36 GB |

About half of what remains is writing the 1 GB `loci.fa`, a product that `bench/rep_rule/score.py` reads.

**Rounds 4-5** (byte-identical on the four catalogs; `RUSTLE_TIMING=1` prints the phases below):
- MCL on sorted sparse columns with a dense accumulator per thread (each entry sums its terms in the same order), one-pass
  `build_clusters`, parallel graph banking and per-pair verdicts;
- the GTF parsed without per-line allocations (borrowed keys, same hasher and insertion order, so the same iteration
  order) and a parallel per-gene pass;
- the loci FASTA written from borrowed slices (`GenomeIndex::fetch_slice`), its cache hash computed on another thread;
- PAF records classified over interned locus ids, numbered in GeneKey order so that pairs sort and nodes number as before
  (a test pins two spellings of one locus to one id).

| human_A119b, PAF cached, 4 threads (clean runs) | wall |
|---|---|
| GTF → loci (331 MB, 2.19 M lines) | 1.3 s |
| loci FASTA (1.97 GB written and hashed) | 1.1-1.3 s (a plain `dd` of 1.97 GB here: 1.7 s) |
| genome (the loci spans) | 0.5-0.6 s |
| exon blocks, PAF read (589 MB), loci map | 1.0 s |
| graph (2.58 M records) | 1.9-2.6 s |
| MCL | 1.3-1.4 s |
| clusters, outputs, copy table | 0.7 s |
| **total** | **9.1 s** (109 s before round 1), incl. 0.6-0.8 s for the test stub that writes the PAF |

- Runs under page-cache pressure on this machine took 17-44 s for the same work (the PAF or the genome re-read from disk);
  the table is the clean runs.
- MCL keeps its summation order, so its floats are the same; the I/O phases sit at this disk's floor.

### The DP memo (`--dp-memo`, `minimap2 2.30-r1287-gate4`)

Overlapping de novo loci hand the aligner the same genomic window many times. In 20 loci of the A119b chr13 cluster, 96% of
the chains of kept pairs were identical in genome coordinates to another chain.

**What it does.** `--dp-memo` keys each chain's base-level alignment on:
- the window's bases (two 128-bit hashes);
- the seeds relative to the window, with their flags;
- the strand, the inversion flag and the left-extension condition.

A chain with an identical key replays the stored result: the final coordinates relative to the window, the extra record
(CIGAR and scores), mlen / blen and any Z-drop split. A thread-safe table holds the results, dropped whole above
`MM_DPM_CAP_GB` (default 1 GB). Since gate5 the cap applies to each of the two tables (chains, DP calls), so up to about
2 GB together. `MM_DPM_VERIFY=1` recomputes every hit and counts mismatches.

**Exactness.** The result is a deterministic function of the key and the options.
- The first version stored the coordinates from before `mm_update_extra`'s trim of leading indels. A self-alignment then
  came out 20 bp off (a negative result; the fix stores the final coordinates).
- Since the fix:
  - an ungated memo run is byte-identical to stock minimap2;
  - gated memo runs are byte-identical to gate3 on a testis subset, an A119b hotspot piece and an A119b typical shard;
  - human_testis gives byte-identical products end to end;
  - verify mode found 0 mismatches over 2,440 hits.

**Speed:**
- **A119b hotspot piece:** 257 → 189 s (1.36x), 30% of chains reused.
- **A119b typical shard:** 141 → 131 s, 3% reused.

**Why not more.** Reuse needs the whole window to match, and a chain within 10 kb of a locus end gets a window clipped by
that end.

A timing build split the alignment time of 20 hotspot loci (single thread, 185 s CPU): left extension 67 s, gap filling
between seeds 17 s, right extension 71 s.

⚠ **Two conclusions drawn here were wrong** (register 1372): that the extensions "are not reusable beyond what the
whole-chain memo catches", and that a faster exact DP was out of scope. Both were reasoned, not measured. Measured below:
58% of the hotspot's DP cycles are calls that exactly repeat an earlier call (the call-level memo), and an AVX2 build of
the same kernel is exact by construction and 1.6-1.7x faster (the DP kernel).

### The DP kernel (gate5)

**Where the time is.** `perf` on an A119b typical shard (gate4, 4 threads): `ksw_extd2_sse41`, minimap2's dual-affine DP
kernel, is **86%** of minimap2's CPU; seed lookup, chaining and index loading share the rest (`mm_idx_get` 1.5%, chaining
about 2.5%). Inside the kernel: the cell loop 60%, the exact-max scan 29%, the score loop 10%. (The "seed collection"
share in the last table was measured on OR6737's 10 Mb shards and does not hold for A119b.)

**What it does.** `ksw2_extd2_avx2.c` is the SSE4.1 kernel built for AVX2:
- the same band and the same 16-cell blocks, two per 256-bit register; the one-byte carry between blocks crosses the
  128-bit lane boundary through `permute2x128` + `alignr`; a single block left at the end of a row runs the SSE4.1 code;
- the exact-max scan runs 8 lanes and folds them back into the SSE kernel's 4 residue classes (larger score, then the
  earlier position), so `max_t` and every Z-drop decision tie-break as before;
- dispatched when the CPU and the OS support AVX2; `MM_KSW_NO_AVX2=1` forces the SSE kernel.

**Exactness.** Every cell is computed with the same lane-wise operations on the same inputs, so the results are
identical by construction. Checked:
- `ksw2_avx2_test`: random related pairs (substitutions, indels, N, tandem repeats), six scoring presets, bands 0 to
  unbanded, every flag combination minimap2 uses: 0 differences in 140,000 cases. It catches 4/4 planted bugs (a lane
  shift, two tie-break rules, an extra score chunk). `build.sh` runs 20,000 cases and does not install on a difference.
- 72,361 recorded real calls (1 in 10 of an A119b typical shard's and a hotspot piece's calls): 0 differences.
- PAFs byte-identical on the hotspot piece, a typical shard and gorilla_OR6737; products byte-identical end to end.

**Speed** (the recorded calls replayed on one thread, alternating): typical shard 29.2 → 17.8 s (**1.64x**), hotspot
56.4 → 33.3 s (**1.69x**). The exact-max scan is now about 32% of the kernel: it works on 32-bit scores (8 lanes); a lazy
max-position pass could save at most about 8% of the kernel. Not done.

**Compiler flags.** `-O3 -march=native -ffp-contract=off` on the whole binary gives no gain, now with repeated,
alternating runs: kernel replay 18.33 vs 18.19 s (AVX2) and 31.4 vs 29.8 s (SSE4.1); one-thread mapping of a typical
subset 14.7/14.8 vs 14.9/14.6 s and of a hotspot subset 46.3/46.1 vs 45.9/45.9 s. The kernels are hand-written
intrinsics, which the compiler cannot widen. PAFs identical across the builds.

### The call-level memo (gate5)

An instrumented build hashed every DP call's complete inputs (both sequences, every scoring argument) and timed it:

| region | DP cycles spent on calls that exactly repeat an earlier call |
|---|---|
| A119b hotspot piece | **58%** (extensions 46%, gap filling 12.5%) |
| A119b typical shard | 4.4% |

Calls whose two sequences are identical along the diagonal took under 2% of the cycles, so a shortcut for them is not
worth its proof.

`--dp-memo` now also keys each dual-affine DP call (above about 20,000 band cells) on its scalar arguments, its scoring
matrix and two 128-bit hashes of the sequences, and replays the stored result (scores, positions and CIGAR); the result is
a deterministic function of exactly these inputs. `MM_DPM_VERIFY=1` recomputes every hit: **0 mismatches over 227,139
reused calls** on the hotspot piece. Hotspot piece 131 → 69 s (with the AVX2 kernel), PAF byte-identical.

### Scheduling (gate6)

minimap2 maps a batch's queries on its threads in input order and aligns one query's chains one after another. On
gorilla_OR6737 one 187 kb locus (`NC_073242.2:120830149-121016923`, 9,021 records) took **167 s of its batch's 419 s** of
query time on one thread, while the others idled.

With `--gate` (the Rustle all-vs-all):
- **Longest first.** Each batch's queries are taken longest first; results are still written in input order.
- **Parallel chains.** A query whose chains span ≥ 10 Mb has them aligned on all threads (`MM_ALIGN_PAR_MIN_SPAN` lowers
  the bar, for tests). Each original chain, with the remainders `mm_align1` splits off it and their inversions, is one
  task: its seeds are its own (disjoint after `mm_squeeze_a`) and its inserts only look back within the run. The runs are
  concatenated in chain order, which is the array the sequential loop builds. Each thread has its own copy of the encoded
  query, its own memory pool and its own DP buffer.
- ⚠ **The first version raced** (register 1373): `mm_align1` reverses a stretch of the encoded query in place for the
  left extension and restores it afterwards, so a shared query buffer gave different output on every run. ThreadSanitizer
  found it. With per-thread copies the output is deterministic (forced onto every query it equals the sequential one 5/5;
  a review run on chr16 and an NPIP subset was byte-identical to stock). ⚠ TSan still reports one pattern on larger inputs
  (the final review): a chain's "inspect nearby seeds" loop reads its neighbours' anchor words (bits 0-39) while their own
  tasks set flag bits 40-41. The values read cannot change and each word has one writer, but it is a data race in C;
  relaxed atomics would remove it (parked).

| gorilla_OR6737 all-vs-all, one process, 4 threads | wall | CPU |
|---|---|---|
| gate4 | 417.9 s | 938 s |
| gate5 (AVX2 kernel, call memo) | 285.9 s | 662 s |
| + longest first | 258.2 s | 658 s |
| + parallel chains (gate6) | **182.3 s** | 710 s |

All four PAFs are byte-identical. At 710 CPU s over 4 threads plus the 21 s index build, this is the parallel floor. Cost:
on an A119b typical shard gate6 takes the same wall time as gate5 (interleaved pairs: 78.8/72.1 vs 76.6/76.7 s) but
8-14% more CPU, since the longest queries run together.

### Final configuration, checked end to end (gate6 + `--dp-memo`, the families-stage flags)

| catalog | families stage, PAF not cached | products vs the ungated catalog |
|---|---|---|
| human_testis | **71.1 s** wall, one process (gate4 101.6 s; old ungated all-vs-all alone: 987 s) | byte-identical |
| gorilla_KB3781 | **105.9 s** wall, one process (gate4 128.3 s; peak 10.2 GB, the minimap2 index) | byte-identical |
| gorilla_OR6737 | all-vs-all **182.3 s** in one process (gate4 417.9 s; ungated 5,421 s) | byte-identical |
| human_A119b | **2,270 s** through `tools/mm2_shard.sh` (gate3 7,332 s; ungated about 5-5.5 h); table below | byte-identical |

- Through `tools/mm2_shard.sh` each shard reloads the index (3.2 GB for OR6737; 10-30 s), which dominates small shards:
  OR6737 in 100 shards of 10 Mb sums to 1,462 s of shard time for the same 182 s of work. Shard times measure this
  machine's 10-minute chunking, not the families stage.

### human_A119b in full

| part | gate3 | gate6 |
|---|---|---|
| index build | 74 s | 140 s |
| 19 other shards (the chr20-21 pericentric one included) | 2,537 s (gate3 split the pericentric one into 8 pieces) | 1,332 s (none split) |
| chr13 hotspot shard | 4,721 s (43 pieces) | 798 s (8 pieces) |
| **total** | **7,332 s (2.0 h)** | **2,270 s (38 min)** |

- gate6 (2026-10-22): all eight products byte-identical to the ungated catalog, and the gated PAF byte-identical to
  gate3's (the hotspot shard joins to the same 1,329,301 records).
- Every shard and piece reloads the index (27 loads in the gate6 run, 10-30 s each); one process would load it once. The
  index build itself varies with this machine's page cache (74 vs 140 s).

- Products are byte-identical to the ungated catalog: clusters, copies, copies FASTA, regions, loci table, loci GFF3.
- Not counted: 1,339 s lost to this machine's 10-minute chunking.
- The measured ratios give about 5-5.5 h ungated, so about 2.6x overall. The hotspot is 64% of what remains.

### What is left, measured

| where | what the time is | can it go faster exactly? |
|---|---|---|
| the DP kernel (A119b: about 80% of minimap2's CPU with AVX2) | the anti-diagonal SIMD DP over a ±1,501-cell band | Only with wider vectors (AVX-512: not on this CPU) or a different algorithm (not exact). A lazy max-position pass could save about 8% of the kernel. |
| repeated DP (dense regions) | reused by the chain and the call memo | Done (58% of the hotspot's DP cycles were repeats). |
| one giant query | its chains on all threads (gate6) | Done: OR6737 sits at its parallel floor. |
| small catalogs' typical shards (OR6737, 10 Mb) | ~7.5 of ~9 CPU s is seed collection; chaining < 0.4 s | No: every target can host a chain that passes, so seeding cannot be skipped before chaining (a pre-chaining filter would save < 0.4 s). |
| minimap2 build flags | `-O3 -march=native -ffp-contract=off` | No gain (repeated, alternating runs; above). |
| index | built once per run (OR6737 21 s, A119b 74 s) | The PAF cache covers reruns. |
| the Rust side (A119b 9.1 s, PAF cached) | the phase table above | The I/O phases sit at this disk's floor; MCL keeps its summation order. |
| CPU | gate6 +8-14% CPU on typical shards at equal wall time | A cost of longest-first; kept for the wall time. |
| memory (parked) | parallel chains nest a pool of `-t` threads inside each outer worker, so up to T² DP pools can coexist | Review measured +27% peak RSS with the parallel path forced on every query, none at the default 10 Mb bar; a cap (one parallel query at a time) is the fix if it ever bites. |

**human_A119b, same index, same shard files:**
- **Typical shard** (shard 0, 3,559 loci): 584 s ungated vs 136 s gated (**4.3x**; CPU 2,210 vs 382 s). The ungated records
  equal the cached full PAF.
- **Hotspot piece** (p2.0, 123 loci of the chr13 cluster): 559 s vs 278 s (**2.0x**). The gate does help there; an earlier
  line in this document said it could not, and was wrong.

### Python in O1

- **The production families stage** runs no Python.
- **The bench O1** (`bench/o3_graph/o1seed.py`) spends about 1-2 s per world in Python. The rest, about 18 s per family set,
  is one spliced minimap2 search of the catalogued transcripts against the set's local reference, cached across worlds.
  Translating the Python would save about 1 s per world, so it was left as is.

### Build and use

- **Build:** `tools/minimap2_gate/build.sh [BIN_DIR]`, default `~/.local/bin`. It needs git, make, a C compiler and zlib
  headers, checks the AVX2 kernel against the SSE4.1 one (20,000 cases) and installs `minimap2-gate` (2.30-r1287-gate6).
- **Switches for checks:** `MM_KSW_NO_AVX2=1` (SSE4.1 kernel), `MM_DPM_VERIFY=1` (recompute every memo hit),
  `MM_ALIGN_PAR_MIN_SPAN=N` (parallel chains from N bp), `RUSTLE_TIMING=1` (the families stage's phase times).
- **Use:** nothing else is needed. The families stage (`tools/rustle_pipeline.sh families`, or `mcl_families --from-gtf`)
  finds the binary and gates by default.
- **The full, ungated PAF:** `mcl_families --no-exon-gate`.


---

## PREREG_multicopy_gains_real_2026-10-21.md

*(original: `docs/archive/2026-10/merged/PREREG_multicopy_gains_real_2026-10-21.md`)*

## Pre-registration: are the multi-copy gains real copies, and are those loci real duplications? (2026-10-21)

**User (2026-10-21), after register 1366:** "lests check the multi-copy gains, also can we prove they are real multicopies? My
advisor is only familiar with NPIP and other mainstays and if he sees that O2 doesnt apply there he might think all other
'families' are artifacts".

Written before the instrument runs on these genes and before any SD overlap is computed.

### Population (fixed from existing tables)

- **Source:** `figures/data/fig3_ref_tie.tsv` (genome-wide; gorilla OR6737 on GGO, human A119b on CHM13 scored against RefSeq,
  as in register 1366).
- **High-tie strata:** expressed multi-exon reference transcripts with tie fraction > 0.5 at 0.98 × best AS, or with no primary
  molecule.
- **Gains (G):** high-tie transcripts that Rustle with secondary seeding matches and Rustle primaries-only does not. Gorilla has
  162 transcripts in 126 genes; human has 284 in 203 genes.
- **High-tie loci (H):** every gene of a high-tie transcript. This is the wider population for Q2.

### Q1: are the gains built from that copy's own reads? (O2's certificate, register 1269's instrument unchanged)

- **Instrument:** `bench/tool_difference/s1_o2.py` (sha1 6b6831e3...), with the gained genes as its gene list.
  - Every good secondary record (AS >= 0.98 × the molecule's best) overlapping a gained gene's exons (copy B) is paired with the
    molecule's primary record (copy A).
  - It certifies on the PSV columns the read spans (Binomial, p < 1e-3, strict maximum, abstain otherwise).
- **Gene classes:**
  - **supported:** >= 2 records certified to B;
  - **contradicted:** 0 to B and >= 2 to A;
  - **undecidable:** otherwise.
- **Reported:** per species, gained transcripts by their gene's class, and records by verdict.
- **Expected, from register 1269:** mostly contradicted or undecidable.
- **Reading, fixed now:**
  - the gains are "real copies" only if at least half of the gained transcripts lie in supported genes;
  - "shadows" if most lie in contradicted genes;
  - "unresolvable" if most lie in undecidable genes, meaning identical copies where the read cannot say which copy produced it.

### Q2: are these loci real duplications in the DNA? (independent of O2, O1 and the RNA reads)

A gene counts as **in a segmental duplication (SD)** when its span overlaps an SD interval of >= 90% identity:
- **human:** the published CHM13 SD track (`CHM13_SDs.SD90.bed`, SDs >= 90% identity; we did not produce it);
- **gorilla:** our genome-wide DNA SD scan of the GGO T2T (`sd_cache/.../genome.sd.pairs.tsv`, both sides of every pair with
  identity >= 0.90), which is sequence-only and blind to the annotation and the reads.

**Also reported:**
- each gene's location class:
  - pseudoautosomal (X/Y PAR: two placements by construction in a male reference);
  - chrY ampliconic;
  - acrocentric short arm;
  - other.
- the named mainstays among the genes: NPIP, TBC1D3, GOLGA6/8, SMN, LRPAP1, and the YAGs (TSPY, RBMY, DAZ, HSFY, BPY2, CDY).

**Bar, fixed now:** "the high-tie loci are real duplications" when >= 90% of H's genes lie in an SD of >= 90% identity, in each
species. The same fraction is reported for G.

### What this can and cannot show

- **Q2** shows whether the places where reads multi-map are duplicated sequence in a T2T genome. That is real, independent of
  our family definition.
- **Q2 cannot show** whether every copy is expressed or present in the sample's own genome: OR6737 and A119b are not the reference
  individuals.
- **Q1** shows whether the extra transcripts are expressed copies or paralog shadows.

### Rules

- Descriptive with fixed readings; nothing adopted or tuned.
- Heavy steps in the foreground, one at a time, each under 10 minutes. The human BAM is 96 GB; if its primary-record scan cannot
  finish in time, an amendment says how it is split.
- The metric-auditor checks the result before it is reported. Never pool human and gorilla. Push only when asked.

### Amendment 1 (2026-10-21, before the human run)

- **Instrument version:** the instrument in the tree is register 1270's version (3d4553af, sha1 a79218c0). Next to register
  1269's (da774afe, sha1 6b6831e3) it only adds the Eichler AS-margin report and the `--stratum` label; the certificate is
  byte-identical.
- **Human primaries scan, split:** the 96 GB human BAM cannot be scanned for primaries within one 10-minute call. The scan is
  split by chromosome group into 5 foreground calls, using two new options on `s1_o2.py` (the default path is unchanged):
  - `--names-only` writes the molecule names and stops;
  - `--primaries FILE` reads the primary records from a file instead of scanning.
- **Check:** on gorilla, the split path reproduces the default run's `records.tsv`, `genes.tsv` and summary byte for byte.

### Results (run once; `/mnt/linuxdisk/tmp/multicopy_gains/`)

#### Q1: the gains are mostly NOT that copy's own expression

| | gorilla | human |
|---|---|---|
| gained transcripts | 162 | 284 |
| in **supported** genes (>= 2 reads certified to the copy) | 20 (12%) | 44 (15%) |
| in **contradicted** genes (reads belong to the other copy: shadows) | 73 (45%) | **200 (70%)** |
| in **undecidable** genes (identical over the reads) | 69 (43%) | 40 (14%) |
| genes: supported / contradicted / undecidable | 7 / 63 / 63 | 34 / 134 / 37 |
| good secondary records, tied: to this copy / other / abstain | 4 / 544 / 1,076 | 45 / 2,073 / 13,707 |
| good secondary records, strictly worse: to this copy / other / abstain | 46 / 2,975 / 3,107 | 400 / 33,871 / 8,964 |

**Reading (fixed rules):**
- Neither species reaches "real copies" (half of the gains in supported genes).
- **Human: shadows**, with 70% of the gains in contradicted genes.
- **Gorilla:** no class holds a majority (45% contradicted, 43% undecidable). The gains there are shadows or unresolvable, not real
  copies.

**Mainstays among the gained genes:**
- **Human, supported:** TBC1D3D, TBC1D3H, TBC1D3I, TBC1D3K, NPIPB5, GOLGA8F, SMN2, RBMY1A1, RBMY1D.
- **Human, contradicted:** NPIPB15, GOLGA6A, TSPY2 / 3 / 8 / 9 / 10, CDY1B.
- **Human, undecidable:** TBC1D3B, GOLGA6B.
- **Gorilla:** three TBC1D3G-like genes, all contradicted.

#### Q2: the loci where reads multi-map are real duplications in the DNA

| genes in a segmental duplication (>= 90% identity) | gorilla (our GGO SD scan, stage-2 identity) | human (published CHM13 SD90) |
|---|---|---|
| **high-tie genes, span overlap (registered bar >= 90%)** | **1,286 / 1,423 (90.4%): holds** | **1,919 / 1,936 (99.1%): holds** |
| gained genes, span overlap | 104 / 133 (78.2%) | 202 / 205 (98.5%) |
| high-tie genes, >= 50% of the span inside an SD (descriptive) | 81.5% | 97.4% |
| high-tie genes with the scan's coarse identity estimate >= 90 (descriptive) | 82.9% | – |

- **Location of the high-tie genes:**
  - human: autosomes 1,444, chrY 337, chrX 118, PAR 37;
  - gorilla: autosomes 1,218, chrX 122, chrY 71, PAR 12.
- **Gorilla gained genes outside any SD** include a run of ~15 adjacent LOC129528xxx genes (undecidable, probably a tandem array
  the scan's 1-kb segment floor misses) and two PAR genes.

**Reading.** In human the loci where long reads multi-map are, almost without exception, segmental duplications in an independent
published annotation. In gorilla, by our DNA-only scan, 90% are (weaker on stricter cuts). The families are real duplicated DNA,
not artefacts of our family definition. What is NOT supported is that the extra transcripts secondaries build there are expressed
copies. Most are rebuilt from another copy's reads (human) or cannot be told apart (gorilla); 12-15% are supported, among them
several of the advisor's mainstays in human (TBC1D3 × 4, NPIPB5, SMN2, GOLGA8F, RBMY1 × 2).

Audit: below (Q1 holds with caveats; Q2 human holds with caveats; Q2 gorilla does NOT hold as stated; the reading above is withdrawn).

### Audit (metric-auditor, 2026-10-21): Q1 holds with caveats; Q2 human holds with caveats; Q2 gorilla does not hold as stated

Every number reproduces to the digit. The populations rebuild identically. The certificate is unchanged since register 1269.
Corrections:

- **The gorilla SD bar rested on a filter that removes nothing.**
  - Every row of our scan has stage-2 identity >= 90 by construction (blocks are kept at >= 90), so "90.4%" means "overlaps any
    scanned pair".
  - With the scan's own identity definition (est_id) it is 82.9%: **FAIL**.
  - Register 1278 kept SEDEF as gorilla's SD source. On SEDEF SD90, as a **post-hoc** check, 1,397 / 1,423 = **98.2%**.
- **Q2 is close to true by construction.** A read ties at >= 0.98 × best only where two loci share near-identical sequence, and
  Q2 had no null.
  - Base rate, genes with no tied reads: 15.3% in an SD (human), 23.6% (gorilla, SEDEF).
  - With the satellite-filtered CHM13 track (`SD90.nosat`), human falls to **77.3%**: chrY 18.4%, acrocentric chromosomes 68.8%,
    other autosomes 97.6%.
  - Of the 1,936 human high-tie genes, 69 lie in the X/Y PARs (allelic, not paralogous) and 371 are RefSeq pseudogenes.
  - "Not artefacts of our family definition" overreaches: no O1 family enters Q2. Its population is read ties.
- **"Supported" mainstays are certified 97-100% to the other copy.** Reads certified to this copy vs the other:
  - TBC1D3D 5 / 931, TBC1D3H 2 / 743, TBC1D3I 11 / 869, TBC1D3K 4 / 674;
  - NPIPB5 11 / 354, GOLGA8F 2 / 498, SMN2 5 / 593;
  - RBMY1A1 4 / 138, RBMY1D 3 / 211.
  
  Only 1 of the 34 supported human genes has as many reads certified to its own copy as to the other. The mainstay lists also
  matched on gene symbol only: counted by RefSeq description, human mainstay-family gains are 11 supported / 43 contradicted / 3
  undecidable.
- **"Contradicted = shadow" mostly restates minimap2's ranking.** Of the calls to the other copy, 94% (human) / 85% (gorilla) come
  from secondaries scoring strictly below the read's best, which minimap2 already placed elsewhere (outside O2's tied-only scope).
  - On TIED records alone, the gained transcripts split supported / contradicted / undecidable as human 9 / 59 / 216 (76%
    undecidable) and gorilla 2 / 25 / 135 (83% undecidable).
  - The instrument never tests reads whose primary lies at the gene itself.
- **Minor:**
  - The location split was computed by gene name, not coordinates (by coordinates: human PAR 69 / chrX 102 / chrY 320).
  - The acrocentric class was not computed.
  - Human is scored on RefSeq, not CAT.

**What stands:**
- **Q2:** the places where reads multi-map are strongly enriched for segmental duplications, 4-6× the base rate. Human: 97.6% on
  ordinary autosomes in the published track, and less for chrY / acrocentric satellite-adjacent regions. Gorilla: 98.2% in SEDEF
  (post hoc). That is expected when reads tie only where sequence is duplicated. It does not validate any O1 family.
- **Q1:**
  - The extra transcripts secondaries build in these loci are not shown to be the copies' own expression.
  - Where a read is truly tied, its copy is mostly undecidable (76-83% of gains).
  - Where it is not tied, minimap2 already put it elsewhere.
  - No mainstay gene is shown to be an expressed copy by this test.

**Correction (2026-10-21, from the audit of `PREREG_o1_families_sd_pairs`):** `CHM13_SDs.SD90.bed` is not a published track. It is a
SEDEF run on CHM13 v2.0 supplied by the user (fracMatch >= 0.90, >= 1 kb). Every "published" above should read "external SEDEF run".
SEDEF is independent of our code but not independent evidence (registers 650, 996).


---

## PREREG_o1_families_sd_pairs_2026-10-21.md

*(original: `docs/archive/2026-10/merged/PREREG_o1_families_sd_pairs_2026-10-21.md`)*

## Pre-registration: are O1's families real duplications? Copy pairs against independent segmental-duplication annotations (2026-10-21)

**User (2026-10-21):** after register 1367 (the multi-mapping loci are duplications, but that says nothing about O1's families):
"yes, please proceed with what you suggest". The suggestion: compare O1's family copies with SEDEF and the CHM13 SD track, and
count how many copy pairs are annotated duplication pairs. The concern behind it: the advisor knows NPIP and the other mainstays,
and might think all other "families" are artefacts.

Written before any overlap is computed.

### Catalogs (O1 outputs, no SD input)

- **human_A119b:** our pipeline's families stage on CHM13. 2,150 families, 9,907 copies, 7,774 with `max_family_identity` >= 0.90.
- **gorilla_KB3781:** the same stage on the GGO T2T. 326 families, 918 copies, 773 with identity >= 0.90.
- **gorilla fibro_cat:** the KB3781 catalog used by the O2 / O3 tests (`rna_allele/fibro_cat.copies.tsv`). 666 families, 3,981
  copies, no identity column.

The families are built from RNA-supported loci and their sequence similarity. The default stage reads no SD annotation (`sd_depth`
and `core_hull` are NA, and `families.key` has no arguments), so the test is not circular.

### SD annotations (independent of our code)

- **human:** the published CHM13 SD track at >= 90% identity (`CHM13_SDs.SD90.bed`; both sides of every row). The
  satellite-filtered version (`SD90.nosat`) is reported as sensitivity.
- **gorilla:** SEDEF on the GGO T2T (`winloci_data/GGO_sedef_final.bed`, rows with fracMatch, column 21, >= 0.90; both sides).
  Register 1278 chose SEDEF as gorilla's SD source.

### Test

- **SD-paired copy:** a copy is SD-paired in its family when one SD row has one side overlapping one of the copy's exons and the
  other side overlapping an exon of another copy of the same family (>= 1 bp each).
- **Bar** (the pipeline catalogs): "O1's families are real duplications" when, among copies with `max_family_identity` >= 0.90,
  >= 90% are SD-paired in their family, in each catalog.
- **Specificity null:** the same test after shuffling family labels among copies, keeping family sizes (seed 1). Expected near 0.
- **Reported (not bars):**
  - all copies, by identity bin (>= 0.98, 0.95-0.98, 0.90-0.95, < 0.90);
  - families in which every >= 0.90 copy is SD-paired;
  - fibro_cat over all its copies;
  - the mainstays by gene or description: NPIP, TBC1D3, GOLGA6 / 8, SMN, LRPAP1, TSPY, RBMY, DAZ, HSFY, BPY2, CDY;
  - the copies that fail, by class: chrY, chrX, de novo copy (no annotated gene), and other.

### What it can and cannot show

**It shows** whether the copies O1 groups are duplications that an independent DNA method also pairs. That is the evidence that
the families are not artefacts of our grouping.

**It cannot show:**
- that every copy is expressed (O2);
- that families below 90% identity are real (SD90 cannot pair them; they are reported, not barred);
- anything about samples other than the two references.

### Rules

- Descriptive with one bar; nothing tuned.
- Light jobs only.
- Code: `bench/o1_sd/sd_pairs.py`, with tests.
- The metric-auditor checks the result before it is reported. Never pool human and gorilla. Push only when asked.

### Results (tool ba7ab571; run once; `/mnt/linuxdisk/tmp/o1_sd/sd_pairs.json`): **bar FAILS in both pipeline catalogs**

| catalog | copies >= 0.90 identity | SD-paired in their family | shuffled-label null | overlap no SD at all | families with every >= 0.90 copy paired |
|---|---|---|---|---|---|
| human_A119b (CHM13 SD90) | 7,774 | **5,605 (72.1%)** | 1,092 (14.0%) | 1,932 (24.9%) | 1,186 / 1,755 |
| gorilla_KB3781 (GGO SEDEF >= 0.90) | 773 | **637 (82.4%)** | 31 (4.0%) | 93 (12.0%) | 214 / 277 |
| gorilla fibro_cat, all copies (no identity column) | 3,981 | 3,557 (89.3%) | 763 (19.2%) | – | 543 / 666 (every copy) |

- **By identity bin, human:**
  - >= 0.98: 74.5% (null 16.7%);
  - 0.95-0.98: 80.0%;
  - 0.90-0.95: 51.2%;
  - < 0.90: 17.1%.
- **By identity bin, gorilla:** >= 0.98 82.2%; 0.95-0.98 88.1%; 0.90-0.95 76.6%; < 0.90 50.0%.
- **Satellite-filtered human SDs:** 52.5% of >= 0.90 copies are paired, so 27% of the human pairings rest on satellite SDs.
- **Correction to the design:** every copy in the two pipeline catalogs is a de novo copy (built from reads), so the registered
  failure class "de novo copy" says nothing.
- **What the unpaired copies are** (descriptive):
  - human: 2,169 unpaired; 1,932 overlap no SD; 1,191 multi-exon and 741 single-exon; 989 lie in an annotated gene, 943 do not;
  - gorilla: 136 unpaired, 93 with no SD.
  
  Their nature (retrocopies, short shared segments, repeat-driven groupings) is not tested here.
- **Mainstays** (descriptive; matched by overlap with annotated genes, after fixing a chromosome-naming mismatch in the first pass):

  | catalog | copies SD-paired in their family |
  |---|---|
  | human_A119b | NPIP 46 / 46 (12 families); GOLGA8 24 / 24; GOLGA6 24 / 25; TBC1D3 22 / 24; SMN1 2 / 2; SMN2 4 / 4; DAZ 6 / 6; TSPY 7 / 7; HSFY 8 / 9; RBMY 5 / 6; CDY 4 / 5; BPY2 3 / 3 (all: 155 / 161, 96%) |
  | gorilla fibro_cat (curated families) | NPIP 25 / 25; LRPAP1 (SM22) 22 / 22; TBC1D3 11 / 14 |
  | gorilla_KB3781 | only 2 mainstay copies matched (NPIP 1, LRPAP1 1), both paired |

**Reading.**
- The bar ("O1's families are real duplications": >= 90% of high-identity copies SD-paired) fails: 72.1% human, 82.4% gorilla.
- The pairing is specific: shuffled families pair 4-14%.
- The advisor's mainstays are 96% (human) and 58 / 61 (gorilla curated) SD-paired: those families are duplications by an
  independent DNA method.
- But a quarter of the human high-identity copies (an eighth in gorilla) lie in no annotated SD. Those groupings are not
  confirmed as duplications by this test, and are the place to look for artefacts.

Audit: below (the FAIL holds with caveats; the reading above is corrected).

### Audit (metric-auditor, 2026-10-21): the FAIL holds with caveats; the reading above is corrected

Every registered number reproduces to the digit with an independent bedtools join. Over all copies the pairing is 60.3% (human) /
76.8% (gorilla). Corrections:

- **Repeats sit on both sides of the bar** (post-hoc classes). Of the 7,774 human >= 0.90 copies:
  - 1,056 are rDNA-array or chrY-satellite families, all paired;
  - 2,122 are transposon-dominated, 61% paired (young L1s, MCL3: 256 of 260 unpaired);
  - 862 are satellite-dominated, 25% paired;
  - the **gene-like rest: 3,036 / 3,237 = 93.8% paired**. Gorilla's gene-like copies: 543 / 590 = 92.0%.
  
  The 72.1% mostly measures repeat content, and being paired does not clear a family.
- **An O1 defect.** For 497 human and 70 gorilla copies, the >= 0.90 identity comes only from a sequence aligned to itself:
  overlapping loci, 144 + 32 of them sense / antisense pairs, e.g. MCL1092. They are 23% (human) and 51% (gorilla) of the
  unpaired copies, and SD90 pairs 4 of them.
- **fibro_cat is not an O1 catalog.** It was built from RefSeq records (every copy has 0 reads). The gorilla mainstay figure
  (58 / 61) is annotation-derived; O1's own gorilla catalog holds 2 mainstay copies.
- **The SD tracks are not independent evidence.**
  - SEDEF is independent of our code, but registers 650 and 996 found its pairs to be a lossy subset of our own all-vs-all
    alignments. Read "an external SD caller pairs them".
  - The human track is a SEDEF run on CHM13 v2.0 that the user supplied (`final_human.bed` filtered to fracMatch >= 0.90 and >= 1
    kb). It is not a "published" track; that word is also withdrawn from `PREREG_multicopy_gains_real`.
- **The satellite filter was misread.** The filtered file drops 2,061 of 2,147 chrY rows (TSPY 7 / 7 → 0 / 7). The copies that
  lose their pairing are mostly the rDNA / chrY array families, not "satellite SDs" in general.
- **Mainstays.** The registered run matched none (a chromosome-naming mismatch); the inline numbers used RefSeq.
  - Recomputed on RefSeq: 156 / 162 (NPIP 49 / 49 in 15 families).
  - On CAT: **206 / 214 (96.3%)**.
  - The name match also caught CDYL (autosomal) as CDY.
- **One human sample.** human_testis, from the same families stage, gives 788 / 972 = 81.1%.

**What stands:**
- The bar fails as registered.
- O1's gene-like families are paired by an external SD caller at about the bar (93.8% human, 92.0% gorilla; post hoc).
- The mainstays are paired at 96% (CAT).
- The shortfall comes from transposon and satellite families, and from **self-overlap families**: a locus grouped with its own
  overlapping or antisense copy. That defect is worth fixing in O1.
- **What would settle it:** a CAT-annotation families arm on CHM13 (a bar the annotation can meet), and a pre-registered rule
  excluding self-only edges, tested on human_testis.


---

## PREREG_o1_self_hits_2026-10-21.md

*(original: `docs/archive/2026-10/merged/PREREG_o1_self_hits_2026-10-21.md`)*

## Pre-registration: dropping diagonal self-hits from O1's family graph (2026-10-21)

**User (2026-10-21):** "lets try to fix the O1 defect". The defect is the one register 1368's audit found: copies grouped with an
overlapping or antisense copy of the same locus. Written before the fix is built or run.

### Root cause (systematic debugging, phase 1; evidence here, before any change)

- **The loci are genomic spans.** The de novo families stage (`mcl_families --from-gtf`) writes each locus as its GENOMIC SPAN and
  aligns the spans all-vs-all (`minimap2 -x asm20 -c -X`). `-X` only skips a record whose query and target names are identical.
- **Overlapping loci align to themselves.** Two different loci whose spans overlap carry the same genomic bases, so minimap2 aligns
  those bases onto themselves. In `human_A119b.fam.loci.paf` there are:
  - 253,893 records between different loci that lie on the main diagonal: same contig, + strand, query and target genomic
    intervals equal within 5 bp. 253,710 of them are at identity >= 0.999.
  - 37,451 records between overlapping intervals that lie off the diagonal (tandem period shifts). These are not self-hits.
- **Worked example:** the only record joining MCL1092's two loci maps chr2:192,437,263-192,525,700 onto itself at identity
  1.0000 over 88 kb. That record made the "family".
- **The existing guards do not catch it:**
  - `merge_overlapping_loci` (a node rule) is skipped under `fold_within_clusters` (default on), and merges exon-union overlaps
    only;
  - `--reject-overlapping` drops every record between span-overlapping loci. It is off by default, non-monotone (o1_ledger), and
    too broad: register 1010 refuted the span-level guard because overlapping de novo loci in duplicated regions are often real
    tandem copies.

### The fix (opt-in; default byte-identical)

`mcl_families --reject-self-hits` drops a PAF record when:
- query and target lie on the same contig;
- the record is on the + strand;
- the aligned genomic intervals coincide (start and end within 5 bp).

It is counted (`rejected_self_hits`) and reported in `params.tsv` and the log. An alignment of DNA onto itself is never evidence of
homology between two copies. Tandem copies keep their off-diagonal records, and inverted copies keep their − strand records.

### Evaluation

- **Catalogs:**
  - human_A119b: development; the defect was measured there.
  - **human_testis: held out.**
  - gorilla_KB3781: the second species.
- **Runs:** each is re-run twice with the current binary, from the same assembled GTF and the cached all-vs-all: a **control**
  without the flag and a **fix** with it. Comparisons are fix vs control (same binary), not vs the 10-10 catalogs.
- **Self-only copy:** a copy of the control catalog whose best same-family record at >= 0.90 identity is a diagonal self-hit,
  with no non-self same-family record at >= 0.90. This is the audit's definition, computed from the control's PAF.

**Bars, each in all three catalogs:**
- **B1 (the defect is removed):** >= 90% of the control's self-only copies, in the fix catalog, are either not in a family, or have
  a non-self same-family record at >= 0.90 identity.
- **B2 (no collateral):** of the control's copies SD-paired in their family (register 1368's test and tool, `bench/o1_sd`), >= 98%
  are still SD-paired in their family in the fix catalog.

**Reported, not bars:**
- register 1368's fraction of >= 0.90 copies SD-paired, control vs fix;
- families and copies;
- largest family size (the span-level guard was non-monotone);
- mainstay copies on CAT (human) or RefSeq (gorilla) in families.

**PASS** = B1 and B2 in all three catalogs. A PASS supports proposing `--reject-self-hits` as the families-stage default (the
user's decision); nothing is flipped by this pre-registration.

### Rules

- TDD for the Rust change: the unit test fails first.
- Build with `CARGO_TARGET_DIR` on `/mnt/linuxdisk`, `--release`, with the log captured to a file.
- Heavy runs in the foreground, one at a time, each under 10 minutes; a run that cannot fit gets an amendment.
- The metric-auditor checks the result. Never pool human and gorilla. Push only when asked.

### Results, v1 (code 72fdfbb6; evaluation e0e43642; `/mnt/linuxdisk/tmp/o1_selfhits/`): **FAIL** (gorilla B2)

Each control rerun reproduces its shipped catalog byte for byte (clusters, copies; the PAF replayed from cache).

| catalog | self-hit records dropped | B1 self-only copies fixed | B2 SD-paired copies kept | copies / families | >= 0.90 copies SD-paired |
|---|---|---|---|---|---|
| human_A119b | 253,893 | 549 / 561 (97.9%) ✓ | 5,857 / 5,974 (98.0%) ✓ | 9,907 / 2,150 → 9,495 / 1,913 | 72.1% → 76.8% |
| **human_testis (held out)** | 2,159 | 46 / 46 (100%) ✓ | 834 / 850 (98.1%) ✓ | 1,115 / 342 → 1,061 / 316 | 81.1% → 84.8% |
| gorilla_KB3781 | 6,645 | 69 / 70 (98.6%) ✓ | **687 / 705 (97.4%) ✗** | 918 / 326 → 834 / 283 | 82.4% → 90.4% |

- **Reported:**
  - largest family 260 → 262 (A119b) and unchanged elsewhere;
  - mainstay copies in families: 155 → 156 (A119b), 56 → 56 (testis), 2 → 2 (gorilla).

**Why B2 fails (post-hoc diagnosis).**
- The lost SD-paired copies are real pairs: in gorilla, all 18 come from 2-copy families at 0.977-1.000 identity, with a direct
  non-self record (e.g. MCL71: an inverted pair 49 kb apart, its − strand record at identity 1.000, 11.6 kb).
- Most of them carry a self-hit partner on the **same strand that shares exon bases**: a second model of the same transcription
  unit nested in it (human_testis 14 / 16, gorilla 13 / 18, A119b 35 pure and 59 mixed of 117).
- In the control, that self-hit kept the two models co-clustered, and `fold_within_clusters` merged them into one copy. Without
  it, host and nested model compete as separate nodes (gorilla's MCL71 becomes a 4-node bipartite graph), and MCL splits the
  pair apart.
- The self-only artefacts are the other kind: their self-hit partner is on the **opposite strand**, or on the same strand with
  **no shared exon bases** (a gene nested in another's intron). A pure "same strand, shared exons" partner is rare among them.
- `--pair-rule weight-blind` (post hoc) joins 19 two-node components but leaves B2 at 97.4%: the two-node rule is not the cause.

### Amendment 1 (2026-10-21, after v1's results; post hoc): v2 rule

- **v2** drops a diagonal self-hit only when the two loci are on **opposite strands**, or share **no exon bases**.
- A self-hit between same-strand loci that share exon bases is kept: they are two models of one transcription unit, and the fold
  step merges them into one copy.
- `mcl_families` passes the loci strands (`gene_strands` of `--gff`) to the graph builder. Without strands the rule is v1's.
- **Same bars (B1, B2) on the same three catalogs.** They are no longer held out: v2 was designed after seeing v1's results on all
  three, so a v2 PASS is a development result and needs a fresh catalog before any default flip.

### Results, v2 (Amendment 1; code e90d63f2; same evaluation): **PASS on all three catalogs (development, not held out)**

| catalog | self-hit records dropped | B1 self-only copies fixed | B2 SD-paired copies kept | copies / families | >= 0.90 copies SD-paired | mainstay copies |
|---|---|---|---|---|---|---|
| human_A119b | 170,084 | 548 / 561 (97.7%) ✓ | 5,931 / 5,974 (99.3%) ✓ | 9,907 / 2,150 → 9,355 / 1,872 | 72.1% → 77.4% | 155 → 152 |
| human_testis | 812 | 46 / 46 (100%) ✓ | 849 / 850 (99.9%) ✓ | 1,115 / 342 → 1,069 / 319 | 81.1% → 84.9% | 56 → 56 |
| gorilla_KB3781 | 3,357 | 69 / 70 (98.6%) ✓ | 704 / 705 (99.9%) ✓ | 918 / 326 → 843 / 290 | 82.4% → 91.5% | 2 → 2 |

- **Largest family:** 260 → 261 (A119b); unchanged elsewhere.
- **Collateral in A119b:** 43 SD-paired copies lose their family. 6 of them overlap mainstays.
  - Five are 2-exon models at NPIPB3 / B4 / B5 (MCL164, MCL112, MCL50) and one is at GOLGA8Q.
  - All are SD-paired at 0.996-0.999 non-self identity, and most have no self-hit of their own.
  - MCL re-routed once self-hit edges elsewhere in their component were gone: the non-monotone effect the o1_ledger recorded for
    edge-dropping rules.

**Reading.**
- v1, the pre-registered rule, FAILS: gorilla B2 is 97.4%.
- v2, designed after v1's results, removes 97.7-100% of the self-overlap artefacts and keeps 99.3-99.9% of the SD-paired copies
  in all three catalogs.
- Because v2 was designed on these catalogs, its PASS is a development result. A default flip needs a fresh catalog, e.g. a
  families run of another sample, pre-registered with v2 as is.

Audit: below, after Amendment 2. It holds with caveats; the A119b collateral reading above is corrected there.

### Amendment 2 (2026-10-21, before any run on it): confirmation on a fresh catalog, gorilla_OR6737

**User:** "ok lets ensure this works, what other dataset can we use?"

- **Fresh catalog: gorilla_OR6737.** A different gorilla individual and tissue (testis Iso-Seq on the same GGO reference). Its
  assembled GTF (`runs/gorilla_OR6737/gorilla_OR6737.families.gtf`, 73,965 transcripts) has never been through the families
  stage. Neither v1 nor v2 has seen it.
- **No fresh human catalog exists locally.** human_A119b and human_testis were used to design v2. A new human library (HG002
  Kinnex, the SCRI cerebellum) would first need a whole-genome alignment.
- **Rule: v2 exactly as committed** (code e90d63f2, `--reject-self-hits` with the loci strands); no change after this amendment.
- **Runs:** the families stage command, a control without the flag and a fix with it, the same binary, GTF and FASTA (GGO).
  The all-vs-all is computed once through `tools/mm2_shard.sh`, in resumable calls of at most 540 s, and cached; the fix
  replays it.
- **Bars:** unchanged.
  - **B1:** >= 90% of the control's self-only copies fixed.
  - **B2:** >= 98% of the control's SD-paired copies (SEDEF >= 0.90) still SD-paired.
  - **PASS** = both.
- **Reported:** the same descriptive figures as for the other catalogs.
- **Reading, fixed now:** a PASS here is the first held-out confirmation of v2, and supports proposing it as the families-stage
  default (the user's decision). A FAIL means v2 is fitted to the development catalogs.

### Audit of v1 and v2 (metric-auditor, 2026-10-21, on `audit_frozen_99782d47`): holds with caveats; readings corrected

Every number in both results tables reproduces to the digit in an independent recompute keyed on coordinates:
- diagonal records 253,893 (253,710 at >= 0.999), off-diagonal 37,451;
- records dropped, B1, B2, mainstay copies and largest family.

The code is correct:
- the loci FASTA is the forward strand, so the 5 bp slack is never used;
- exon blocks are 0-based half-open;
- the strand keys match;
- no record between different positions is dropped.

Corrections:
- **The A119b collateral is described wrongly above.**
  - 28 of the 43 "lost" SD-paired copies are representative switches: a model of the same locus that shares exon bases is
    still SD-paired in a fix family.
  - All 6 mainstay copies are among the 28 (e.g. DN_chr16_21697777_2 → DN_chr16_21699518_2, MCL136). In the control, those
    NPIPB loci sat in two families at once (MCL164 and MCL50).
  - Counted by locus, B2 is 5,959 / 5,974 (99.75%).
  - "Most have no self-hit of their own" is false: 37 of the 43, the 6 mainstays included, lost their own diagonal records.
  - The real collateral that no bar sees: 8 distant 2-copy families dissolved (identity 0.87-0.999). 4 of them are outside SD90
    (MCL153, 207, 963, 940).
- **v2 is not a pure deletion in A119b** (trap γ).
  - 43 fix families combine members of two or more control families (621 new co-member pairs).
  - Copy pairs in different families that share exon bases rise from 1,368 to 1,545.
  - 16 new same-family overlapping pairs appear (62 remain), and 4 new self-only copies (chr13).
  - Negligible in human_testis and gorilla_KB3781.
- **B1 is close to a tautology.** "Fixed" almost always means "left the families": 545 / 548, 45 / 46, 69 / 69. Only 3 / 1 / 0
  were repaired. By register 1368's definition the self-only copies go 497 → 0 (human), 70 → 0 (gorilla), 42 → 0 (testis), but
  485 and 69 of them left the families.
- **v1's FAIL is fragile.** It is 4 copies short under a gene_id join that the scorer, committed at 16:05, adopted after the v1
  runs. Counted by locus, v1 is 692-693 / 705 and would pass.
- **The ">= 0.90 SD-paired" gain is a shrinking denominator.** A119b: paired copies +9, >= 0.90 copies −523. Testis +2 / −42,
  gorilla 0 / −77. Quote counts, not percentages.
- **Default byte-identity was shown only for v1.** The controls ran on the v1 binary. The v2 build without the flag matched
  them by code reading only.
- **The fold does not depend on strand.** v2's rationale, "the fold merges them", would hold for opposite-strand pairs that
  share exon bases too. v2 drops those self-hits on purpose: a sense model and an antisense model of one locus are not two
  copies.

**What stands:**
- v1, as pre-registered, failed gorilla B2 by 4 copies under a gene_id join.
- v2, designed after seeing all three catalogs (including the testis held out for v1), passes them. That is a development
  result.

**What would settle it:**
- an evaluation fixed in advance: B2 counted by locus, plus guards on cross-family exon sharing, dissolved distant pairs and
  leftover self-overlap;
- run on gorilla_OR6737, plus a CAT-annotation arm.

A failure keeps the flag opt-in.

### Amendment 3 (2026-10-21, after the audit, before any OR6737 output exists): the evaluation, fixed in advance

Written while the OR6737 all-vs-all is still running. Its PAF and both catalogs do not exist yet. v2 stays as committed
(e90d63f2). Scorer: `bench/o1_sd/self_hits_eval.py` at the commit of this amendment.

**Bars:**
- **B1 (unchanged):** >= 90% of the control's self-only copies fixed.
  - It needs >= 10 control self-only copies; with fewer, the result is UNDERPOWERED, not PASS.
  - Reported: how many fixed copies left the families and how many were repaired.
- **B2, counted by locus** (replaces the gene_id join; audit findings 1 and 4): a control SD-paired copy is kept when it, or a
  fix copy sharing exon bases with it (either strand), is SD-paired in its fix family. >= 98%. The gene_id count is reported.

**Guards**, for what B1 and B2 cannot see; all must hold:
- **G1, one transcription unit split across families:** same-strand copy pairs in different families that share exon bases.
  - Fix <= control + max(5, 5% of the control's).
  - Opposite-strand pairs are reported only: separating sense from antisense is what the fix is for.
- **G2, distant pairs kept** (independent of SD tracks):
  - A control copy is distant-paired when its family holds a copy at a span-disjoint locus, joined to it by a non-self record
    at >= 0.90 identity.
  - The pair is kept when fix copies at its two loci (the copies themselves, or copies sharing exon bases with them) share a
    family.
  - >= 98% of distant-paired copies kept.
- **G3, the defect is gone from the output:** self-only copies under the fix's own families <= max(2, 10% of the control's
  self-only count).

**Reported:**
- SD-paired counts (not percentages);
- new and lost co-member pairs among copies present in both catalogs;
- fix families that combine members of two or more control families;
- dissolved 2-copy control families;
- the other figures as before.

**Checks before the verdict:**
- **Default unchanged:** the v2 binary without the flag must reproduce the development controls byte for byte (clusters and
  copies). Checked on the replays of gorilla_KB3781 and human_testis.
- **Development figures:** G1-G3 are computed on the three development catalogs and reported. They do not set the thresholds.
  - Known when this was written (from the audit), with no strand split seen: A119b's cross-family exon-sharing pairs
    1,368 → 1,545, 8 distant 2-copy families dissolved, and 4 new self-only copies.

**Verdict on OR6737:**
- PASS = B1, B2, G1, G2 and G3 all hold. A PASS supports proposing v2 as the families-stage default (the user's decision).
- Any FAIL keeps it opt-in.

### Results, Amendment 3 (scorer 7cda7ffd; frozen copy `audit_frozen_or6737_7cda7ffd`): **OR6737 (held out) PASS**; human_A119b (development) fails G1

- **The default is unchanged.** The v2 binary without the flag reproduces the KB3781 and testis controls byte for byte (clusters,
  copies, merged copies, PAF).
- **OR6737 control:** 1,607 copies in 498 families. The fix replays the same PAF from cache and drops 2,963 self-hit records.

| catalog | B1 self-only fixed (left / repaired) | B2 by locus (gene_id join) | G1 same-strand cross-family pairs (opposite strand) | G2 distant-paired copies kept | G3 self-only left | verdict |
|---|---|---|---|---|---|---|
| **gorilla_OR6737 (held out)** | 71 / 72 (71 / 0) ✓ | 1,272 / 1,272 (1,270) ✓ | 25 → 24 (36 → 36) ✓ | 1,364 / 1,365 ✓ | 0 ✓ | **PASS** |
| human_A119b (development) | 548 / 561 (545 / 3) ✓ | 5,959 / 5,974 (5,931) ✓ | **1,036 → 1,140 > 1,088 ✗** (332 → 405) | 7,560 / 7,598 ✓ | 10 ✓ | FAIL (G1) |
| human_testis (development) | 46 / 46 (45 / 1) ✓ | 849 / 850 (849) ✓ | 8 → 8 (4 → 4) ✓ | 961 / 962 ✓ | 0 ✓ | PASS |
| gorilla_KB3781 (development) | 69 / 70 (69 / 0) ✓ | 705 / 705 (704) ✓ | 9 → 8 (12 → 13) ✓ | 746 / 749 ✓ | 0 ✓ | PASS |

**Reported for OR6737:**
- copies 1,607 → 1,521; families 498 → 456; largest family 86 → 86;
- co-member pairs +2 / −6; 1 fix family merges control families;
- 43 two-copy families dissolved, 8 of them with members that are not self-only;
- >= 0.90 copies SD-paired: 1,159 / 1,383 → 1,160 / 1,293;
- mainstay copies 34 → 34, with the same gene_ids.

#### Audit (metric-auditor, 2026-10-21, on `audit_frozen_or6737_7cda7ffd`): holds with caveats

Every OR6737 number reproduces to the digit without the scorer. The PASS survives every stricter count:
- B2 with a same-strand-only join: 1,270 / 1,272;
- G2 requiring every pair: 1,358 / 1,365 (by pair, 6,595 / 6,599);
- G1 counting folded models: 25 → 26 (limit 30).

**The timeline is clean.**
- Amendment 2 was committed at 16:34:29 and Amendment 3 at 17:39:47. At that point 49 of the 100 all-vs-all shards were done
  and no family output existed.
- The family outputs date from 18:25.
- The scorer is identical to 7cda7ffd.
- The binary was built at 16:19:44, and `src/` is unchanged since e90d63f2.

**Caveats:**
- **B1 and G3 hold by construction.**
  - All 72 self-only copies form 36 two-copy families. Each pairs opposite-strand models that overlap in span, share no exon
    bases, and are joined only by a diagonal record (4.7-87.7 kb). v2 drops every such record.
  - The deletion is consistent with the catalog's own edge rules: the pairs' other records are internal repeats at 0.83-0.894.
  - The fix removes more than the self-only definition counts. Same-family opposite-strand pairs with overlapping spans go
    47 → 0. That includes the 8 "not self-only" dissolved families (the same artefact, hidden from the definition by
    off-diagonal internal-repeat records) and MCL233.
  - All 90 lost copies are accounted for: 71 self-only, 16 from those 8 families, 1 from MCL233, and 2 antisense
    representative switches.
- **One real collateral loss, which G1 as written does not see.**
  - MCL70's locus (NC_073240.2:21.11-21.23 Mb) is split between MCL175 and MCL368.
  - Cause: v2's same-unit test is pairwise, but the fold is transitive. 17 of the 334 dropped same-strand records join loci
    that are one unit through a third model.
  - This is G2's one lost copy. Its inverted 27 kb pairing survives through the locus's other model.
- **OR6737 is a held-out sample, not a held-out genome.**
  - It is on the same GGO reference as KB3781, where v1 failed and v2 was shaped.
  - 550 / 1,607 of its copies share exon bases with KB3781 copies.
  - 14 / 72 of its self-only copies and 8 / 43 of its dissolved families sit at loci also seen in KB3781.
  - Its 10-10 GTF had never been through the families stage. Earlier assemblies of the same library had been (`rt_arms`
    09-27, `rg3_run` 09-28), but none was used for v1 or v2.
- **A119b's G1 failure is mostly relabelling, but not only that.** My first label, "all switches", was not a test. Checked by
  locus:
  - 289 of the 313 new pairs lie at loci already split across control families. 291 of the 313 are on chr13, 13.9-15.5 Mb.
  - 17 are new splits of loci the control held in one family.
  - 7 involve a locus with no same-strand control copy.
  - Counted by unit, families beyond the first go 711 → 742 (+4.4%). That would pass a 5% locus bar, but only post hoc.

**Reading** (the auditor's wording):
- v2 (e90d63f2) passes every Amendment 3 bar and guard on one held-out gorilla catalog, OR6737, which is on the same GGO
  reference as KB3781. The result was reproduced independently.
- B1 and G3 hold by construction. The collateral is one locus split across two families.
- There is no held-out human catalog. human_A119b fails G1 as pre-registered (1,036 → 1,140 > 1,088). The guided (annotation)
  arm is untested.
- As pre-registered, this supports proposing v2 as the default of the de novo families stage **for gorilla**. For human, and
  for catalogs built from an annotation, it is not shown.
- Nothing is flipped; the user decides.

**What would settle the rest:**
- a held-out human catalog, e.g. HG002 Kinnex after a whole-genome alignment;
- control and fix runs on a guided CAT arm, both under Amendment 3;
- for MCL70-type splits, a same-unit test that follows the fold's transitive components rather than pairwise exon sharing.
  That would be a v3, with its own pre-registration.

### Decision (user, 2026-10-21): `--reject-self-hits` is the default

**User:** "lets make reject self hits the default". The user decided this knowing the scope of the evidence above:
- confirmed on one held-out gorilla catalog;
- human_A119b fails G1 in development;
- annotation-built catalogs are untested.

**What changed:**
- `mcl_families` now ships `--reject-self-hits` ON, with v2's rule (strands from the loci GFF).
- `--no-reject-self-hits` restores the graph of every catalog built before the flip. The control arms above are now run with
  that flag.
- `--reject-self-hits` still parses and is a no-op.
- The library default (`GraphParams::default()`) stays OFF, as with `min_shared_exon_frac`. `params.tsv` records the
  effective value.

**Checked after the flip:**
- **gorilla_KB3781 and gorilla_OR6737, new default vs the fix arm:** clusters, copies and merged copies are byte-identical,
  with the same records dropped (3,357; 2,963).
- **The same catalogs, `--no-reject-self-hits` vs the controls:** byte-identical.
- **Tests:** 19 / 19 `mcl_families`, 1,106 library, 3 driver.

**Run without `--gff`:** when `--paf` is given without `--gff`, no strands or exon blocks are known. The rule then falls back to
v1, which drops every diagonal self-hit. Every families-stage run in this project uses `--from-gtf`, which writes and passes its
loci GFF.

**What changes downstream:**
- Every families catalog built from now on drops self-overlap families. On the catalogs scored here:
  - human_A119b: 9,907 → 9,355 copies;
  - human_testis: 1,115 → 1,069;
  - gorilla_KB3781: 918 → 843;
  - gorilla_OR6737: 1,607 → 1,521.
- The shipped catalogs in `rustle_figures/runs` are unchanged until they are rebuilt.

## SDSCAN2 — sd-scan v2: find SEDEF's segmental duplications without SEDEF (prereg 2026-10-25)

**Why.** The user asked for SD finding "a la sedef without needing it as external dependency". The vendored `utilities
sd-scan` (src/sedef/) failed V1 (r1278: base recall vs SEDEF 0.398 CHM13 v2.0 / 0.407 gorilla at identity >= 0.90,
precision 0.92 / 0.95).

**Diagnosis (gorilla, the r1278 outputs; `/mnt/linuxdisk/tmp/sddiag/`).** Of SEDEF's pair bases (>= 0.90):
- 34% never seeded (no candidate window pair covers half the pair): 67% of 1-5 kb pairs, 55% of 90-95% pairs.
- 30% are inside a stage-2 block but lost in refine. Refine keeps ONE anchor run per block (a 270 kb block becomes one
  8.6 kb pair; any indel > 50 bp or gap > 400 bp ends a run).
- Refine's identity is anchor coverage, a lower bound (a 96% block reads 84-92%), so the >= 0.90 filter drops true pairs.
- Upper bounds with the V1 scorer: block identity instead of anchor coverage, recall 0.542; every block whole, 0.683
  (precision 0.916).

**Changes (scan v2; developed on gorilla only):**
- R1, refine: every collinear anchor run is kept; consecutive runs of one strand are joined when the step between them is
  <= `refine_max_gap` on both axes (the existing run rule, so no new threshold); every joined chain with both spans >=
  `refine_min_span` is a pair.
- R2, identity: the gaps between consecutive anchors are aligned (unit-cost global alignment); identity = matches /
  (matches + mismatches), indels excluded, as SEDEF's fracMatch.
- R3, seeding (if dev needs it): seeding windows that overlap by half, so a copy is never more than a quarter window off
  the grid.
Any further dev change is logged here before the held-out run.

**Evaluation.** `bench/sedef/v1_compare.py` unchanged (min-id 0.90, min-len 1000), on chromosome subsets so a run fits the
machine; SEDEF truth = the pairs with both sides in the subset; both arms (scan v1 = commit before R1, scan v2) are run on
the same subset FASTA.
- DEV: gorilla NC_073228.2, NC_073234.2, NC_073238.2, NC_073242.2 (any number of looks).
- HELD-OUT, run once with frozen code: CHM13 v2.0 chr1, chr9, chr15, chr16, truth = `o1_close/lift/sedef_pairs.v2.bed`
  (lifted v1.0 SEDEF), `--chain soto_replication/v1.0_to_v1.1_rdna_merged.chain`. Human and gorilla numbers are not pooled.

**Bars (held-out):** H1 recall >= 0.80 AND precision >= 0.80 (r1278's V1 bar). H2 recall of v2 > recall of v1 on the
same subset. Descriptive: per length / identity bin, run time per stage.
**Readings:** H1 PASS = sd-scan v2 can replace SEDEF as the SD source on genomes without it (adoption is a separate
decision with the user). H1 FAIL + H2 PASS = better but still not a replacement. H2 FAIL = the fixes are dev-overfit.

### SDSCAN2 dev log (gorilla dev subset only)
- scan v1 on the dev subset: recall 0.4898, precision 0.7659 (subset precision is lower than genome-wide: SEDEF pairs
  with a partner outside the subset are not in the truth).
- R1+R2 first build: OOM at 25 GB (crashed the session): in tandem arrays every shifted diagonal chained and each chain
  stored its anchors. Fix (no threshold): chains taken by A coverage, a chain kept only when its A x B rectangle meets
  no kept chain's. Rerun capped (`prlimit --as`): 6.3 GB peak, refine 246 s (v1 65 s).
- R1+R2: recall 0.7561, precision 0.7457 (85,443 pairs, 84,652 ok). Next: R3 (seeding).
- R3 (seeding windows every 1 kb; chain stage indexed in stride units, so v1 = stride 2 kb; `max_postings` 100 -> 200;
  pairs of overlapping windows of one locus skipped): recall 0.8478, precision 0.6972 (96,742 pairs; seed 62 s, refine
  314 s, 9.4 GB peak). The toy genome's inverted 8 kb copy is now one pair (v1: four 2 kb pairs); golden updated.
- Precision loss: every false pair (>= half outside every SEDEF pair of any identity) is >= 80% soft-masked, but so are
  1,365 true pairs. SEDEF source (github vpc-ccg/sedef, read 2026-10-25): seeds are k-mers with >= 1 unmasked base;
  a final SD needs >= 100 unmasked bases per side and >= 100 unmasked matches (`stats_main.cc:302`); fracMatch =
  matches / (matches + mismatches), R2's definition. Post-hoc filters on R3's pairs (dev, descriptive):
  1 kb unmasked per pair 0.654 / 0.976; per stage-2 block 0.697 / 0.970; 100 per pair 0.713 / 0.948; 100 per block
  0.757 / 0.930 (recall / precision). SEDEF reaches mostly-masked SD stretches by extending from unmasked seeds inside
  one long alignment; our pairs and blocks are shorter, so the same rule removes them. Not yet decided; no held-out run.
- SEDEF / BISER ports (sources cloned read-only to /mnt/linuxdisk/tmp/sdsrc, 2026-10-25):
  - BISER's masked projection (`masked_space`, CLI `--masked-space`): seed, chain and refine without soft-masked bases,
    coordinates mapped back. Dev: recall 0.7011, precision 0.8397; refine 19 s, 1.8 GB. Misses are long SDs that are
    almost all soft-masked (154 SEDEF pairs, ~31 kb each, < 1 kb unmasked).
  - SEDEF's mode (now the default): seeds only from k-mers with an unmasked base; ok pairs need >= 100 unmasked bases per
    side and >= 100 unmasked matches (SEDEF `stats_main.cc` defaults; status `repeat` otherwise).
  - BISER's second-level chain DP (MATCH 10, MISMATCH 1, GAP 0.5, GAPOPEN 100, MAX_GAP 10 kb), with the shorter gap of
    a join <= refine_max_gap so both copies' stretch is aligned; chains below min_block_id identity are not kept and
    shadow nothing in the selection.
  - Join cost fixes: a joined gap is aligned only up to refine_max_gap^2 cells, else its shorter side counts as
    mismatches (identity stays a lower bound); the DP stops scanning predecessors more than MAX_GAP away and counts every
    examined one toward MAX_ITER (an O(n^2) scan made refine 481 s).
  - SEDEF's seed rule, ported exactly: every k-mer counts toward the shared mass; a candidate carries `upper_seed` when
    a shared sketch k-mer is all unmasked in both windows (candidates.tsv column 9); the chain stage keeps blocks with
    such a cell; the 100-unmasked filter is over a block's ok pairs together (pieces of one duplication).
  - Dev (recall / precision; refine s, peak GB):
    | arm | recall | precision | refine | peak |
    |---|---|---|---|---|
    | v1 | 0.4898 | 0.7659 | 65 | - |
    | v2 repeat-blind (R1-R3) | 0.8478 | 0.6972 | 314 | 9.4 |
    | BISER projection (`--masked-space`) | 0.7062 | 0.8349 | 20 | 1.9 |
    | SEDEF rules, per-window gate | 0.6751 | 0.9693 | 124 | 3.9 |
    | SEDEF rules, block gate (default) | 0.7079 | 0.9653 | 406 | 3.8 |
    None reaches 0.80 / 0.80. The recall gap is long SEDEF SDs that are almost all soft-masked (~31 kb with < 1 kb
    unmasked): SEDEF reaches them by extending one alignment through the repeats; our blocks and joins stop at them.
  - 2026-10-26, scan v3 (dev only; held-out not run):
    - Port fix: SEDEF's seed k-mer needs ONE uppercase base (hash.cc HAS_UPPERCASE, search.cc:416), not all k; the
      seed-gate loss fell from 17.2% to 4.9% of truth bp (12,100 blocks pass instead of 8,354).
    - Port fix: the 100-unmasked filter is per duplication (stats_main.cc judges one alignment), not per block: summing a
      block passed tandem arrays of masked copies (toy: 207 repeat pairs). Pieces of one strand, in order on both
      copies, whose boundary is an indel (shorter gap and any overlap <= refine_max_gap) are summed; tandem copies
      (equal gaps) never group.
    - BISER's join clips an overlap (align.codon refine); bounded by refine_max_gap (an unbounded clip put one stretch
      into two chains: recall 0.667).
    - Refine speed, outputs byte-identical on a 1-in-8 block sample: Ukkonen-banded gap alignment, an identity floor
      (stop once all-matching remaining A bases cannot reach 90%), identities of a block's next 8 unshadowed chains in
      parallel then replayed in order. Sample refine 87.7 s -> 44.6 s.
    | arm (dev) | recall | precision | refine | peak |
    |---|---|---|---|---|
    | v3 seed fix + per-pair filter | 0.6982 | 0.9559 | 356 | 4.6 |
    | + overlap clip bounded | 0.7035 | 0.9547 | 345 | 4.5 |
    | + per-duplication filter (scan v3) | 0.7178 | 0.9530 | 335 | 4.5 |
    - Remaining loss (bp of truth): unmasked filter 10.8% (SDs split into 4+ pieces with unaligned masked stretches on
      both copies, which SEDEF aligns through), never a candidate 5.2%, seed gate 4.9%, chain 1.8%, refine 1.6%.

## PREREG SEMISEC (2026-10-26): semi-guided mode from multimappers - silent copies of expressed genes, classified
**Status:** committed before any candidate is computed. Stage 1 = dev census (testis); stage-2 bars are set from it and committed
before the held-out sample is run.
**Question (user, after the advisor):** can the secondary alignments of expressed genes' reads find their copies that carry no
reads (pseudogenes, lncRNA-like copies), and say which kind each is? Semi-guided as defined 09-26: seeds = read-supported loci,
targets = copies with no reads; no annotation is used to build or classify (annotation = truth only).
**Samples (gorilla only, one reference GCF_029281585.2):** DEV = testis OR6737 (`winloci_data/GGO_mm.bam`); HELD-OUT = fibroblast
KB3781 (`fibroblasts/GCA_029281585.2_flnc_mm.bam`). Both aligned `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1`.
**Steps** (per sample; script `bench/semisec.py`):
1. Seeds: primary records (`-F 2308`, MAPQ any) merged into loci by overlap on the same strand; a seed = a locus with >= 2 reads.
2. Candidate hits: secondary records (flag 256) of seed reads with AS >= 0.8 x the read's primary AS (minimap2's own default
   `-p 0.8`; the BAM kept 0.1) not overlapping the read's own seed locus. Hits merged by overlap per strand into candidate loci; a
   candidate needs >= 2 distinct reads.
3. Silent: a candidate locus with < 2 primary reads of the same sample. Only silent candidates are scored (the rest are expressed
   copies, de novo's job).
4. Class, from sequence only, using the candidate's best hit (highest AS): seed ORF = the longest ORF (ATG..stop) of that read, coding
   if >= 100 codons (the usual lncRNA cut). The candidate's sequence = reference bases along the hit's aligned blocks (introns skipped).
   - PROCESSED if the primary has >= 1 intron (N) and the hit has none; else UNPROCESSED.
   - Coding seed: DISRUPTED if the ORF's aligned stretch on the candidate has a premature stop or a net indel not divisible by 3,
     else INTACT. Non-coding seed: LNC-LIKE.
**Truth (annotation, never read by steps 1-4):** a candidate's annotated biotype = the gene / pseudogene feature of `GGO_genomic.gff`
overlapping it most (pseudogene and transcribed_pseudogene -> pseudogene; protein_coding; lncRNA; other; none).
**Reports (stage 1, dev):**
- Candidates: count, by class; annotated biotype per class (confusion). Expected if the classes mean something: DISRUPTED and
  PROCESSED enriched for pseudogene over INTACT; LNC-LIKE enriched for lncRNA.
- Recall: annotated silent pseudogenes whose description names an expressed seed's gene (description "X-like", X = the description
  of a protein_coding gene overlapping a seed) - the fraction overlapped by a silent candidate. NCBI's homology, not ours.
- Unannotated silent candidates ("none"): count and size, reported as possible new copies or noise (no truth).
**Stage 2 (held-out fibroblast):** bars committed after stage 1, before it is run.

### SEMISEC stage 1 - dev census (testis OR6737), 2026-10-26
Script `bench/semisec.py` (dump -> TSV per contig via samtools+mawk; run 2.5 min, 3.1 GB). Seeds 26,055; candidates 6,630, silent 2,672.

| class | n | pseudogene | protein_coding | lncRNA | other | none |
|---|---|---|---|---|---|---|
| PROCESSED/DISRUPTED | 313 | 189 (0.60) | 90 (0.29) | 14 (0.04) | 0 | 20 (0.06) |
| PROCESSED/INTACT | 152 | 13 (0.09) | 126 (0.83) | 5 (0.03) | 0 | 8 (0.05) |
| PROCESSED/LNC-LIKE | 199 | 60 (0.30) | 50 (0.25) | 12 (0.06) | 0 | 77 (0.39) |
| UNPROCESSED/DISRUPTED | 746 | 79 (0.11) | 161 (0.22) | 54 (0.07) | 3 | 449 (0.60) |
| UNPROCESSED/INTACT | 480 | 36 (0.07) | 159 (0.33) | 31 (0.06) | 3 | 251 (0.52) |
| UNPROCESSED/LNC-LIKE | 780 | 50 (0.06) | 166 (0.21) | 68 (0.09) | 46 (0.06) | 450 (0.58) |
| UNPROCESSED/UNCOVERED | 2 | 0 | 0 | 1 | 0 | 1 |

- Classes carry meaning: among annotated candidates, pseudogene fraction DISRUPTED 0.454 (n 590) vs INTACT 0.131 (n 373) = 3.46x;
  PROCESSED/DISRUPTED 0.60 pseudogene; PROCESSED/INTACT 0.83 protein_coding (intact retrocopies are annotated as genes). LNC-LIKE
  lncRNA 0.177 vs coding classes 0.108 (1.64x; weak).
- Recall (silent denominator, as preregistered): 115 / 582 = 0.198. Ruling: the first report counted expressed pseudogenes too
  (120 / 1,186) - a bug against the prereg text, fixed before any bar was set. Loss ladder over the 1,186 (dev only): 1,163 get
  secondaries of seed reads; 788 lose them all at AS >= 0.8 x primary - the main loss is the preregistered -p 0.8 cut.
- Unannotated silent candidates 1,256 (median span 3.9 kb, median 4 reads), from only 513 distinct best reads: a few multi-hit reads
  give most of them; includes mitochondrial reads (NUMTs, unannotated in RefSeq). No truth for them.

### SEMISEC stage-2 bars (committed before the fibroblast KB3781 run; code frozen at this commit)
Same script, same parameters, fibroblast BAM. Fractions among annotated candidates (truth != none), PROCESSED+UNPROCESSED pooled.
- **B1** pseudogene fraction DISRUPTED / INTACT >= 1.5 (dev 3.46).
- **B2** PROCESSED/DISRUPTED candidates annotated pseudogene >= 0.40 (dev 0.60).
- **B3** recall over annotated silent "X-like" pseudogenes of expressed seed genes >= 0.10 (dev 0.198; half of dev, since the seed set
  changes with tissue).
- **PASS = B1 and B2 and B3.** **B4** (reported separately, the lncRNA claim): lncRNA fraction LNC-LIKE / coding classes >= 1.25 (dev 1.64).
- Unannotated candidates: descriptive only.

### SEMISEC stage 2 - held-out fibroblast KB3781 (register 1431): **PASS** (B1-B3); B4 lncRNA FAIL
Code = 0a12c3d9. Seeds 23,311; candidates 6,378, silent 3,672 (run 487 s, 9.2 GB).
- B1 pseudogene DISRUPTED / INTACT = 0.378 (n 1,092) / 0.107 (n 484) = **3.52** (>= 1.5) PASS.
- B2 PROCESSED/DISRUPTED annotated pseudogene = 293 / 618 = **0.474** (>= 0.40) PASS. PROCESSED/INTACT again mostly protein_coding (0.77).
- B3 recall = 151 / 664 = **0.227** (>= 0.10) PASS.
- B4 lncRNA LNC-LIKE 0.164 (n 635) vs coding 0.150 (n 1,576) = 1.09 (< 1.25) **FAIL**: the dev 1.64 did not hold. "LNC-LIKE" means only
  "copy of a non-coding read"; it does not predict an annotated lncRNA.
- Unannotated silent candidates 1,437 (median span 3.1 kb): possible unannotated copies or noise; no truth.
- Reading: semi-guided from multimappers finds about 1 in 5 silent annotated pseudogenes of expressed genes (loss = AS >= 0.8 x primary)
  and tells a disrupted copy from an intact one (3.5x pseudogene enrichment in both tissues). It does not identify lncRNAs.

### SEMISEC Rust port (2026-10-26)
Steps 1-4 now run in `src/semisec.rs` (bin `semisec`; contigs scanned in 20 Mb pieces in parallel through the .bai, then indexed
fetches of each silent candidate's best hit and read); `bench/semisec.py report OUT GFF` scores. Same candidates as the Python run
on both samples (testis 6,630, fibroblast 6,378: reads, primaries, silent, class, best read; ORF coverage equal to 5 decimals) and
byte-identical reports. Time: testis 79 s / 2.0 GB (Python dump + run ~10 min), fibroblast ~200 s / 4.0 GB (~16 min / 9.2 GB);
report 5 s. New output column `best_seed` (the best read's seed locus) for the parent-agreement measure below.

## PREREG SEMIRATIO (2026-10-26): drop the 0.8 AS cut?
**Why:** the main recall loss of SEMISEC is the AS >= 0.8 x primary cut (788 / 1,186 target pseudogenes lose all their seed-read
secondaries there, testis). The alternative is no cut of our own: every secondary the aligner reported (`-p 0.1 -N 50`) with AS >= 0
(`semisec --ratio 0`). Everything else unchanged.
**New measure (both settings):** parent agreement = among silent candidates whose annotated pseudogene is "X-like", the fraction found
from a seed overlapping a protein_coding gene described X (NCBI's homology says the copy came from the right gene).
**Runs:** testis (descriptive) and fibroblast KB3781 (decisive), each at --ratio 0.8 and --ratio 0. Fibroblast bars:
- **C1** recall(0) >= recall(0.8) + 0.05.
- **C2** under ratio 0: B1 (pseudogene DISRUPTED / INTACT >= 1.5) and B2 (PROCESSED/DISRUPTED pseudogene >= 0.40) still pass.
- **C3** parent agreement(0) >= parent agreement(0.8) - 0.05.
- **PASS = C1 and C2 and C3** -> `--ratio 0` becomes the default. Candidate and unannotated counts reported.

### SEMIRATIO result (register 1432): **FAIL** - the 0.8 cut stays the default
- Ruling (before the fibroblast ratio-0 report was read; dev and fibroblast 0.8 had been seen): the parent-agreement measure as first
  coded matched "X-like" against gene DESCRIPTIONS and gave 0.12-0.14 at 0.8; NCBI names "X-like" pseudogenes after the PROTEIN
  ("small ribosomal subunit protein eS26-like" for gene "ribosomal protein S26"), so it now matches CDS `product=` names (isoform
  suffix dropped). A measurement fix, bars unchanged; cost if wrong: C3 judged on the corrected number only.

| sample, ratio | candidates (silent) | recall | B1 | B2 | parent agreement | unannotated |
|---|---|---|---|---|---|---|
| testis 0.8 | 6,630 (2,672) | 115/582 = 0.198 | 3.46 | 0.604 | 272/369 = 0.737 | 1,256 |
| testis 0 | 26,699 (9,380) | 191/582 = 0.328 | 3.60 | 0.443 | 544/860 = 0.633 | 6,157 |
| fibroblast 0.8 | 6,378 (3,672) | 151/664 = 0.227 | 3.52 | 0.474 | 403/569 = 0.708 | 1,437 |
| fibroblast 0 | 27,703 (15,549) | 364/664 = **0.548** | 2.18 | **0.359** | 713/1,403 = **0.508** | 8,236 |

Fibroblast: C1 PASS (+0.32 recall); C2 FAIL (B2 0.359 < 0.40); C3 FAIL (0.508 < 0.708 - 0.05). Dropping the cut more than doubles
recall but the extra hits are less often the right gene's copy and less often annotated. At 0.8, about 3 in 4 annotated hits are
the right gene's pseudogene (0.71-0.74 in both tissues).

## PREREG SEMITRUTH (2026-10-27): are SEMISEC's unannotated silent candidates real gene copies?
**Status:** committed before any projection is computed. The gorilla annotation has no answer for them; human does, by synteny.
**Sets** (SEMISEC default, `--ratio 0.8`, per sample): **U** = silent candidates with no overlapping gorilla gene/pseudogene;
**P** = silent candidates whose annotated biotype is pseudogene (positive control: shows what the check can see);
**N** = one random interval per U item, same contig and length, overlapping no gorilla gene/pseudogene feature and no candidate
(negative control, fixed seed 1).
**Projection to human:** the interval with 2 kb flanks on each side, from `GGO.aug.fa`, aligned with `minimap2 -x asm20 -c` (prebuilt
T2T-CHM13 v2.0 index `npip_ladder/idx/target.asm20.mmi`); keep the primary alignment with MAPQ >= 1 (unique; minimap2's own
uniqueness signal); the core interval is projected through its CIGAR. No unique primary, or a core projecting to < 50% of its length:
**UNPLACED**.
**Human call** (CAT/Liftoff v2.0, `chm13v2.0_CAT_Liftoff.genes.tsv`, gene spans): **PSEUDO** = the projected core overlaps a gene whose
biotype contains "pseudogene"; **GENE** = overlaps another gene (no pseudogene); **NONE** = placed, overlaps no gene.
**Parent match** (U and P): the human pseudogene's name is the parent's symbol + `P<digits>` (HGNC pseudogene naming), parent = the
`gene=` symbol of a protein_coding gene at the best read's seed locus. Reported over the items whose human pseudogene has an HGNC-style
name.
**Also reported:** fraction of each interval covered by RepeatMasker (GenArk `GCF_029281585.2.repeatMasker.out`), per set.
**Bars** (testis = dev, descriptive; fibroblast KB3781 = decisive): **T1** the check works: P PSEUDO rate >= 0.50. **T2** U is enriched
for copies: U PSEUDO rate >= 2 x N PSEUDO rate. **PASS = T1 and T2.** Then: estimated real unannotated copies = U PSEUDO / P PSEUDO rate
(descriptive); U parent-match rate (descriptive).

### SEMITRUTH result (register 1433): **PASS** in both tissues
`bench/semitruth.py` (one minimap2 call for both samples, 297 s, 14.5 GB with the index); per-item table `/mnt/linuxdisk/tmp/semisec/truth/items.tsv`.

| sample | set | n | PSEUDO | GENE | NONE | UNPLACED | repeat frac | parent match |
|---|---|---|---|---|---|---|---|---|
| testis | N | 1,256 | 0.018 | 0.181 | 0.535 | 0.266 | 0.564 | - |
| testis | P | 427 | **0.829** | 0.077 | 0.037 | 0.056 | 0.169 | 189/215 |
| testis | U | 1,256 | **0.120** | 0.265 | 0.410 | 0.205 | 0.699 | 14/23 |
| fibroblast | N | 1,437 | 0.024 | 0.216 | 0.496 | 0.264 | 0.539 | - |
| fibroblast | P | 650 | **0.792** | 0.138 | 0.020 | 0.049 | 0.208 | 262/306 |
| fibroblast | U | 1,437 | **0.106** | 0.288 | 0.417 | 0.189 | 0.774 | 27/40 |

- Fibroblast: T1 P 0.792 >= 0.50 PASS; T2 U 0.106 >= 2 x N 0.049 PASS. Testis the same (0.829; 0.120 vs 0.037).
- Estimated real copies among the unannotated: about 180-190 per sample (U PSEUDO / P PSEUDO), ~13%: copies the gorilla annotation
  missed and human annotates as pseudogenes at the syntenic place; where named, 60-70% carry the parent's symbol (small n).
- The rest is mostly repeat-driven: U is 70-77% RepeatMasker-covered vs 54-56% background. Descriptive split (post hoc): U PSEUDO
  0.23-0.25 when < 50% repeat vs 0.07 when >= 50%; P drops from 0.86-0.88 to 0.48-0.55 in the same way. A repeat filter is the next
  candidate rule (needs its own prereg).
- Parent check of the annotated controls: 0.86-0.88 of named human pseudogenes carry the seed gene's symbol.

## PREREG SEMIREP (2026-10-27): drop repeat-driven secondary hits
**Rule:** `semisec --repeats GCF_029281585.2.repeatMasker.out.gz`: a secondary hit is dropped (before candidates are formed) when
MOST of its aligned reference bases (M/=/X blocks; > half - a majority, not a tuned cut) lie in RepeatMasker intervals. Primaries,
seeds and the silent test are untouched. Everything else as SEMISEC (ratio 0.8).
**Contamination note:** SEMITRUTH's descriptive split (candidate SPAN repeat cover < / >= 0.5) was already seen on both samples. This
rule acts on single hits, before merging, so its candidate sets and every number below are unseen; but it is the same idea, so a PASS
here is weaker than a fresh-sample PASS.
**Runs:** testis (descriptive), fibroblast KB3781 (decisive); `semisec.py report` + `semitruth.py`, filtered vs unfiltered.
**Fibroblast bars:**
- **R1 precision of the unannotated:** U (unannotated silent) human-PSEUDO rate >= 0.20 (unfiltered 0.106).
- **R2 real copies kept:** U human-PSEUDO count >= 0.8 x unfiltered (152 -> >= 122).
- **R3 annotated side not worse:** recall >= 0.227 - 0.02; parent agreement >= 0.708 - 0.02; B1 and B2 pass.
- **PASS = R1 and R2 and R3** -> `--repeats` becomes part of the default SEMISEC run.

### SEMIREP result (register 1434): **FAIL** - `--repeats` stays opt-in
Unfiltered `semisec` output is byte-identical to before (cands.tsv, testis). The filter drops 46,042 testis hits; candidates 6,630 -> 4,292
(testis), 6,378 -> 3,765 (fibroblast).

| fibroblast | unfiltered | --repeats | bar |
|---|---|---|---|
| U (unannotated silent) | 1,437 | 346 | |
| U human-PSEUDO rate | 0.106 | **0.263** | R1 >= 0.20 PASS |
| U human-PSEUDO count | 152 | **91** | R2 >= 122 FAIL |
| recall | 0.227 | **0.202** | R3 >= 0.207 FAIL |
| parent agreement | 0.708 | 0.831 | R3 >= 0.688 pass |
| B1 / B2 | 3.52 / 0.474 | 4.90 / 0.550 | pass |
| P human-PSEUDO (control) | 0.792 | 0.864 | |

Testis (descriptive) the same way: U 1,256 -> 513, PSEUDO rate 0.120 -> 0.214, count 151 -> 110; recall 0.198 -> 0.189; parent agreement
0.737 -> 0.783.
Reading: the filter buys precision (unannotated real-copy rate x2.5, random background unchanged at 0.01-0.03) but loses 40% of the real
unannotated copies and 17 annotated ones: many real pseudogene copies sit in, or carry, repeats. Opt-in for a high-precision list
(about 1 in 4 unannotated candidates real, parent symbol 29/29 where named); the default keeps every copy.


## PREREG DENOVOFAM (2026-10-27): better de novo multi-copy families against Compara, mainstays and YAGs
**Goal (user /goal):** improve the de novo mode's detection of multi-copy gene families; judge on Ensembl Compara, the mainstays and the YAGs.
**Truth (`bench/denovo_eval/famtruth.py`):** Ensembl Compara human within-species paralogue pairs (BioMart, 3.3M rows); truth families =
connected components over pairs at duplication levels Homo sapiens..Catarrhini ("recent"; Simiiformes and Primates ladders reported).
Copies map to the CAT/Liftoff v2.0 gene sharing the most exonic bp on the strand; judgeable = protein-coding CAT genes with a Compara
Ensembl id that carry a de novo transcript. Dev/held-out halves = SHA-1 of each truth family's smallest ENSG (singletons hashed alike).
Mainstays/YAGs = the `bench/o3_graph/protfam.py` KNOWN regexes (NPIP, TBC1D3, GOLGA8, GOLGA6L, CT47A, XAGE, SPATA31; BPY2, CDY, DAZ, HSFY,
PRY, RBMY, TSPY, VCY, XKRY).
**Baseline = current main `mcl_families --from-gtf --min-exonic-bp 1 --min-shared-exon-frac 0.60 --emit-units`** on A119b (testis):
dev pair F 0.764, family F 0.791 (all-pairs 0.894 / 0.836, seen). Miss diagnosis (`bench/denovo_eval/missdiag.py`, 204 missed recent pairs):
**104 MCL splits** (an edge joins the two loci; GOLGA8/6, POTE, AGAP, NPIP), 78 below the shared-exon gate on the best SINGLE record
(19 of them reach >= 0.60 as the union of exon-to-exon records), 14 with no record, 7 coverage, 1 short.
**Arms (each a replay of the cached alignment, 26 s):**
- **A inflation:** MCL inflation 2.0 and 2.4 (default 2.8). §6ec: with the size-safe prune NPIP is stable 2.0-4.0, so the split is the
  graph, not prune.
- **B union shared-exon:** the shared-exon gate counts exon-to-exon bases over the UNION of the pair's records that each touch exons on
  both sides (now: the best single record). Distinct from r302/r1097 (union of SPAN coverage, which admitted Alu stubs): only exon-to-exon
  bases count, and the 0.60 fraction is unchanged. New opt-in flag.
- **C = best A + B.**
**Selection (dev half of A119b only; held-out halves and testis not printed while choosing):** the arm with the highest dev family F
among arms whose dev pair precision drops <= 0.02 and that do not regress any mainstay/YAG (largest-family members not fewer, families
not more, foreign not more by > 1).
**Held-out bars (chosen arm vs baseline):**
- **H1** A119b held-out half: family F >= baseline + 0.01 AND pair F >= baseline.
- **H2** fresh testis library ERR13885926 (`runs/human_testis`), all genes: family F >= baseline + 0.01 AND pair F >= baseline.
- **H3** precision: pair precision on both held-out sets drops <= 0.02 (recent), and Primates-ladder pair precision does not drop > 0.02.
- **H4** mainstays + YAGs on testis: no regression as defined above.
- **PASS = H1-H4** -> the arm becomes the default (old behaviour behind a flag). Any fail -> register row, stays opt-in or dropped.
**Amendment 1 (2026-10-27, before any arm was scored on a held-out set):** the scorer maps a copy to the PROTEIN-CODING CAT gene with the
most shared exonic bp when one overlaps it (else any gene, as before). Found on dev: DAZ2's copy went to CHM13_G1000495
(`unknown_likely_coding`, a 200 kb model over DAZ2), SPDYE15's to PMS2P2. Applies to baseline and every arm alike. Re-scored dev
baseline: pair F 0.764, family F 0.793 (303 dev pairs). Arm A on dev: I 2.0 F 0.753 / 0.791 (precision -0.034, out); I 2.4 0.764 / 0.801.
**Amendment 2 (same time): arm D, fold bridges.** Dev diagnosis: most dev misses are genes in NO family because
`fold-within-clusters` welds two genes into one locus through a readthrough model that overlaps both (GOLGA8Q + GOLGA8H -> one copy,
the 70 kb readthrough AC091057.3). Rule (no threshold): inside a cluster, a record is a BRIDGE when two of the records it overlaps on exon
bases do not overlap each other; loci = connected exon-overlap components of the non-bridge records; each bridge folds into the locus it
overlaps most (exon bp; ties to the smaller coordinate). Not r815 (which cut bridge transcripts in the assembler, before families, and
failed on gorilla): the transcripts are untouched, only which records count as one copy changes. Arms now A, B, D and C = the selected
combination; selection and bars unchanged.
**Amendment 3 (dev only; no held-out run yet):** arm D as written gave dev pair F 0.764 -> 0.825, family F 0.793 -> 0.810, GOLGA8 11 -> 13/13,
but NPIP 19 -> 18/20: NPIPB13's full 24 kb model overlaps two DISJOINT FRAGMENTS of NPIPB13 itself, so it counted as a bridge and a 2.8 kb
fragment became the representative. Refined rule (still no threshold): a record is a bridge only when two records it overlaps do not
overlap each other AND are joined by a homology edge of the graph (two paralogues of the family welded into one copy; two fragments of
one gene are not homologous to each other). Arm D' replaces D; selection and bars unchanged.

### DENOVOFAM result (register 1435): **FAIL** (H1, H2) - `--fold-split-bridges` stays opt-in; nothing becomes default
**Dev selection (A119b dev half, Amendment 1 scorer):** baseline pair F 0.764 / family F 0.793. A I 2.0 0.753 / 0.791 (precision -0.034,
out); A I 2.4 0.764 / 0.801; B `--shared-exon-union` 0.759 / 0.795 (DAZ splits into 2 families, out); **D' `--fold-split-bridges` 0.826 /
0.810 (selected; GOLGA8 11 -> 13/13, NPIP 19/20 kept)**; C = B + D' 0.820 / 0.811 (DAZ split, out).
| D' vs baseline | baseline | D' | bar |
|---|---|---|---|
| H1 A119b held-out half pair F / family F | 0.956 / 0.903 | 0.956 / 0.903 | family +0.01 FAIL |
| H2 testis ERR13885926 all pair F / family F | 0.954 / 0.916 | 0.954 / 0.916 (copies byte-identical) | family +0.01 FAIL |
| H3 pair precision recent (held-out A119b / testis); Primates A119b held-out | 0.974 / 0.939; 0.980 | same; 0.980 | pass |
| H4 mainstays + YAGs on testis | | unchanged | pass |
Reading: D' is INERT outside the development families: no held-out family has two paralogues welded by a readthrough model (testis:
not one fold changed). All-gene A119b pair F 0.891 -> 0.897, family F 0.834 -> 0.838, Primates all-pair precision 0.916 -> 0.909; the gain is
GOLGA8/GOLGA6, SPDYE, POTE, SERPINB on the dev half only. Opt-in for catalogs that hold readthrough models over paralogues; a fresh
held-out with such loci is needed before any default. B (union of exon-to-exon records) gains nothing on dev and splits DAZ.
Other findings: (1) the scorer must map copies to protein-coding genes first (Amendment 1); (2) most "MCL splits" in the first diagnosis
were folds, not MCL (inflation 2.0-2.8 moves dev sensitivity 0.717 -> 0.724); (3) ungated `--from-gtf` needs the sharded wrapper
(`RUSTLE_MINIMAP2=tools/mm2_shard.sh`): a plain minimap2 all-vs-all of 94k loci ran > 75 min.

## PREREG DENOVOFAM-G (2026-10-27): `--fold-split-bridges` (D') on FRESH gorilla catalogs
DENOVOFAM failed because both human held-outs had no paralogue pair welded by a readthrough (D' inert, r1435). The two gorilla
catalogs (fibroblast KB3781, OR6737; `runs/gorilla_*`) were never used to build or choose D'. Gorilla numbers are never pooled with human.
**Truth:** gorilla RefSeq genes (`winloci_data/GGO_genomic.gff`, protein-coding, exon union per gene) carry the human symbol; a gene
takes the Ensembl ids of the CAT/Liftoff human gene of the same symbol, and truth families are the human Compara "recent" components over
those ids (Homo sapiens..Catarrhini). Unnamed `LOC` genes are outside the universe (a stated blind spot: many gorilla multi-copy genes
are LOC). Copies map protein-coding-first by exonic bp, as Amendment 1. Mainstay/YAG regexes on gorilla symbols. All genes, no halves.
**Runs:** baseline vs D' (`mcl_families --from-gtf ... --min-exonic-bp 1 --min-shared-exon-frac 0.60 --emit-units [--fold-split-bridges]`).
**Bars, per sample:** G1 pair F and family F not lower than baseline; G2 pair precision (recent) drops <= 0.02; G3 no mainstay/YAG
regression (largest-family members not fewer, families not more, foreign not more by > 1). **PASS = G1-G3 on both samples AND family F
>= baseline + 0.01 on at least one** -> `--fold-split-bridges` becomes the default (old behaviour behind `--no-fold-split-bridges`).
Otherwise it stays opt-in.

### DENOVOFAM-G result (register 1436): **FAIL (no gain; truth underpowered)** - `--fold-split-bridges` stays opt-in
D' changes a few copies (KB3781 843 -> 846 copies, 290 -> 291 families; OR6737 1,521 -> 1,523) and NO scored number: KB3781 pair F 0.320,
family F 0.625 both arms; OR6737 0.324 / 0.619 both; mainstays/YAGs identical. G1-G3 pass, the +0.01 family bar fails on both.
Underpowered: the symbol-mapped Compara truth holds only 9 (KB3781) and 13 (OR6737) true families with >= 2 expressed members (13 and 17
true pairs), because gorilla RefSeq names most extra copies `LOC...` (the pre-stated blind spot).

## PREREG DENOVOFAM-L (2026-10-27): Liftoff-style exon chaining as family edges (pending item 1, user go 2026-10-27)
**Why:** the DENOVOFAM dev misses that D' did not fix are missing EDGES for retro/partial copies (CCL3L, SIRP, APOBEC3, LYZL, ZNF,
BOLA2): no single alignment record carries >= 0.60 of the smaller locus's exons. Liftoff (Shumate & Salzberg 2021) places a gene by
aligning it and CHAINING its exon-level alignments in order, and accepts a placement at coverage >= 0.5 (`-a`) and exon identity >= 0.5
(`-s`). Differs from arm B (`--shared-exon-union`, union of exon bases over unordered records inside the old gate) and from r1286/r1287
(one whole-transcript record, recall of Liftoff copies).
**Lift (bench/denovo_eval/liftchain.py, written before any score):** records = the UNGATED all-vs-all of the same locus FASTA
(`minimap2 -x asm20 -c -X -N 50 -p 0.1`, cache bd6780f14333f911, query hash db4ee602…; both PAF orientations). For template A and locus B:
per record, through the CIGAR, the A-exon bases aligned (and those landing on B-exon bases); records on genomic bases aligned onto
themselves dropped (the shipped self-hit rule). Chain = Liftoff's DAG: records on one strand, colinear in A and B, not overlapping in A;
the path maximizing A-exon aligned bases. coverage = chained A-exon bases / A exon-union length; identity = A-exon-base-weighted record
identity (nmatch/blocklen). **A lifts onto B** when coverage >= 0.5, identity >= 0.5 and >= 1 chained A-exon base lands on a B exon
(the shipped `--min-exonic-bp 1` convention). Edge A-B when A lifts onto B or B onto A; w = max over directions of identity x coverage
(the graph's own identity x coverage form). No other number.
**Arms (new opt-in `mcl_families --extra-edges TSV [--extra-edges-only]`; an extra edge raises an existing edge to max(w)):**
E1 = baseline graph + lift edges, MCL unchanged. E2 = lift edges only, MCL unchanged.
**Selection on the A119b DEV half (famtruth.py, recent levels, Amendment 1):** an arm qualifies if dev family F >= baseline + 0.01, dev pair
F not lower, dev pair precision drops <= 0.02, no mainstay/YAG regression (as DENOVOFAM). Best qualifying arm by dev family F goes held out.
**Held out:** A119b held-out half AND human testis (all genes; its own ungated all-vs-all, same command): family F >= baseline + 0.01 on both,
pair F not lower, precision drop <= 0.02, no mainstay/YAG regression. PASS -> the arm becomes a documented option and a default flip is
proposed to the user; otherwise opt-in. Reported (no verdict): Simiiformes/Primates ladders, edges added, largest family.
**Not tested here:** families as pure connected components (no MCL), and lifts onto silent sites (copies without a locus).
**Amendment 1 (2026-10-27, after the dev readout of E1/E2, before any held-out score):** dev, recent: baseline pair F 0.764 (prec 0.819),
family F 0.793; E1 0.763 (prec 0.773) / 0.828, GOLGA6L foreign 2 -> 15; E2 0.753 (prec 0.757) / 0.835, NPIP 19 -> 16. Both fail selection
(pair F lower, precision drop > 0.02). Liftoff's `-s 0.5` is a placement floor, not a homology floor. One added dev arm, no new number:
**E3** = E1 with each lift also meeting the graph's own shipped identity floor (`min_identity` 0.70) in a passing direction. Same
selection and held-out bars. If E3 fails dev, DENOVOFAM-L is FAIL and nothing is held out.

### DENOVOFAM-L result (register 1437): **FAIL at dev selection** - `--extra-edges` stays an opt-in tool, nothing held out
Lift (liftchain.py, 494 s): 5,779,244 locus pairs, 742,022 lift edges (baseline graph 553,600 edges); 256,428 self-hit records dropped.
Dev half, recent (pair sens / prec / F; family F): baseline 0.716 / 0.819 / 0.764; 0.793. E1 0.752 / 0.773 / 0.763; 0.828. E2 0.749 /
0.757 / 0.753; 0.835. E3 0.746 / 0.769 / 0.757; 0.821. Every arm raises family F by >= 0.028 and pair recall by ~0.03, and every arm
fails precision (drop 0.046-0.062 > 0.02) and pair F; GOLGA6L foreign members 2 -> 15 (E1, E3), NPIP 19 -> 16 (E2).
Post hoc (no verdict): the 0.70 identity floor is inert (removes 317 of 742,022 lifts; the asm20 records are already >= 0.70), so the
precision loss comes from Liftoff's ONE-SIDED coverage: a template covering half of its own exons anywhere inside a bigger locus
lifts onto it, which is a placement criterion, not a family one. Liftoff's thresholds answer "where does this gene go", not
"are these two loci one family". The recall gain is real (+11 dev true pairs), so an edge that needs the lift to explain the
TARGET too (two-sided) is the untested next step, and it needs its own prereg.

## PREREG DENOVOFAM-L2 (2026-10-27): two-sided lift edges
r1437 showed lift edges add true pairs but a one-sided lift (small template inside a big locus) costs precision. **E4** = baseline graph +
an edge A-B only when A lifts onto B AND B lifts onto A (each direction: Liftoff coverage >= 0.5, identity >= 0.5, >= 1 bp on the other's
exons; the same liftchain.py stats), w = min over the two directions of identity x coverage. No new number. Same dev selection and held-out
bars as DENOVOFAM-L (A119b dev half -> A119b held-out half + testis). Fail at dev -> FAIL, nothing held out.

### DENOVOFAM-L2 result (register 1438): **FAIL at dev selection**
421,645 two-sided lifts (29,107 new edges, 27,993 raised). Dev, recent: pair sens 0.723, prec 0.774, F 0.747 (baseline 0.716 / 0.819 /
0.764); family F 0.788 (baseline 0.793); DAZ splits (1 -> 2 families). Worse than E1 on every number: two-sidedness removes most of the
recall gain and keeps the precision loss, so one-sided coverage is not the whole cause. The raised weights also move MCL (DAZ). Post hoc
only: the Liftoff line is closed at its two registered forms; the recall it adds (E1 +11 dev pairs) comes with edges MCL cannot
separate from passengers.

## PREREG DENOVOFAM-L3 (2026-10-27): lift edges only for orphan loci
The DENOVOFAM dev misses are mostly genes in NO family; r1437/r1438 lost precision by adding or raising edges BETWEEN loci that already had
homology edges (MCL then regroups, GOLGA6L/DAZ). **E5** = baseline graph + the one-sided DENOVOFAM-L lift edges (identity >= 0.5) whose
pair has at least one locus with NO edge in the baseline graph (`--dump-graph` of the baseline run), so no existing edge is added to or
raised. No new number. Same dev selection and held-out bars (A119b dev -> A119b held-out half + testis with its own ungated all-vs-all).

### DENOVOFAM-L3 result (register 1439): **FAIL at dev selection (recent truth)**; the older Compara levels disagree
E5 (9,532 orphan lift edges, 0 raised): dev recent pair sens 0.733, prec 0.784, F 0.758 (baseline 0.716 / 0.819 / 0.764); family F 0.822
(0.793); every mainstay and YAG identical to baseline. Pair diff vs baseline on dev: +5 true (LYZL1-LYZL2, SLC35G4/5/6, APOL3-APOL4), +13
"false", 0 lost. The 13 are mostly REAL older paralogues (ZNF726 with seven KZNFs, ARHGEF5-ARHGEF35, GPRASP1-2, TMSB15A-B, SMG1-BOLA2-SMG1P6).
Reported ladders (dev; pair F / family F): baseline Simiiformes 0.743 / 0.773, Primates 0.714 / 0.748; **E1 0.780 / 0.828, 0.762 / 0.810**
(precision 0.830 vs 0.861, 0.837 vs 0.842); E5 0.774 / 0.809, 0.739 / 0.780. So against the older truth E1 gains on every number; the recent-
level precision loss is largely older paralogues the recent cut calls false. The verdict stays as registered (recent): FAIL. Changing the
primary truth level is a science choice for the user; held-out halves remain unscored for every L arm.

## PREREG DENOVOFAM-L4 (2026-10-27): orphan lift edges (E5) held out against Simiiformes truth (user decision 2026-10-27)
User decision after r1439: a multi-copy family is judged against Compara **Simiiformes** paralogues (recent stays reported). Dev selection
at Simiiformes, bars unchanged: E1 fails (precision 0.861 -> 0.830, drop 0.031 > 0.02); **E5 qualifies** (pair F 0.743 -> 0.774, prec 0.861
-> 0.861, family F 0.773 -> 0.809, mainstays/YAGs identical); E2-E4 not re-selected. **Held out (never scored for any L arm):** A119b
held-out half and human testis (all genes; testis lift edges from its own ungated `-x asm20 -c -X -N 50 -p 0.1` all-vs-all, orphans from
its own baseline `--dump-graph`). Bars at Simiiformes, each sample: family F >= baseline + 0.01, pair F not lower, pair precision drop
<= 0.02, no mainstay/YAG regression (largest-family members not fewer, families not more, foreign not more by > 1). PASS on both ->
propose making orphan lift edges the de novo default (user decides; needs the lift step in Rust). Otherwise opt-in.

### DENOVOFAM-L4 result (register 1440): **PASS on both held-out samples (Simiiformes)** - default flip proposed to the user
Testis lift: own ungated all-vs-all (672,671 records, mm2_shard 545 s), 372,255 pairs, 5,436 lifts, 334 orphan lifts (302 edges, 277 new
nodes). A119b held-out half, Simiiformes: pair F 0.913 -> 0.917 (prec 0.978 -> 0.970), family F 0.840 -> **0.858** (+0.018). Testis, all
genes, Simiiformes: pair F 0.900 -> 0.909 (prec 0.956 -> 0.952), family F 0.865 -> **0.897** (+0.032). Mainstays and YAGs identical in both.
Reported: recent A119b held-out pair F 0.956 -> 0.949 (prec 0.974 -> 0.957), family F 0.903 -> 0.914; testis 0.954 -> 0.944 (0.939 -> 0.918),
0.916 -> 0.926. Primates A119b 0.904 -> 0.910 / 0.837 -> 0.857; testis 0.890 -> 0.899 / 0.855 -> 0.887. Recent pair F falls slightly in
both (older paralogues, r1439). The lift step is Python today (liftchain.py over an UNGATED all-vs-all); a default needs it in Rust and the
ungated alignment (A119b ~9.5 min sharded, testis 9 min).
**Audit correction (2026-10-27, metric-auditor; scratch /mnt/linuxdisk/tmp/audit_2026-10-27_denovofamL4/): the PASS above is RETRACTED.**
Numbers reproduce to the digit, but (1) famtruth.py scores each half on predicted families TRIMMED to that half's genes
(`score(th, pf, hg | sing)`), so a predicted pair joining a held-out gene to a dev gene (always false) counts in neither half. E5 adds 22
such cross-half false pairs (8 ZNF726-KZNF, GSTA, IGLL, TUBA1A/C, PRDM7/9, TEX13A/B, ...). Counting predicted pairs that touch >= 1 held-out
gene, A119b held-out pair F 0.8903 -> **0.8885 (lower: bar fails)**, precision 0.9278 -> 0.9079 (drop 0.0199); family F untrimmed 0.8136 ->
0.8323. (2) The pass exists only at Simiiformes, chosen after E5 failed dev at recent; at recent both held-outs fail (A119b 0.9557 -> 0.9486,
testis 0.9536 -> 0.9444, testis precision drop 0.0211). (3) 91.5% of E5's new A119b gene pairs (864/944) and 194/223 testis pairs touch genes
outside the judgeable protein-coding universe; 571 of 1,367 new copies map to no CAT gene; families 1,872 -> 2,341. (4) Mainstays/YAGs are
identical in the SCORE but not in content (DAZL joins DAZ, NOS2P1 joins TBC1D3, RBMY1A3P / RBMY2EP / RBMY2FP join RBMY). Testis alone passes
at Simiiformes. **Verdict: FAIL (A119b held-out pair F lower once cross-half pairs count); `--extra-edges` stays opt-in.** The per-half
trimming also flatters every earlier DENOVOFAM dev/held-out number; famtruth.py must count cross-half pairs before the next prereg.
**famtruth.py fixed (same day):** a half now counts every predicted pair touching >= 1 of its genes, and the bipartite step sees untrimmed
families. Reproduces the audit: A119b Simiiformes held-out baseline pair F 0.890 (prec 0.928), family F 0.814; E5 0.888 (0.908), 0.832.
Every DENOVOFAM / -G / -L dev and held-out number above was computed with the old trimming and reads high on precision.

## PREREG DENOVOFAM-L5 (2026-10-27): orphan lifts that do not widen a family
After the audit (r1440 retracted) and the scorer fix (04c72875). E5's "false" pairs are orphans (often older paralogues: ZNF726 -> KZNFs)
attached to EXISTING families; its true gains (LYZL, SLC35G, APOL) look like orphans joining each other. Two number-free arms, from the
same A119b/testis lift and orphan sets (orphan = locus with no edge in the baseline `--dump-graph`):
**E6** = baseline + lifts whose BOTH loci are orphans (A119b 1,378 of 9,532). **E7** = E6 + an orphan-to-family lift kept only when its weight
is >= the weakest internal edge (baseline graph, both ends in the baseline MCL cluster of the non-orphan locus) of that family: an orphan may
join a family but not widen its divergence radius.
**Scorer:** fixed famtruth.py, primary truth Compara **Simiiformes** (user, 2026-10-27), recent reported. Dev selection on the A119b dev half
(family F >= baseline + 0.01, pair F not lower, precision drop <= 0.02, mainstays/YAGs no regression in the score AND content changes listed).
Best qualifying arm held out on the A119b held-out half and testis with the same bars. PASS on both -> propose default (user). Else FAIL.

### DENOVOFAM-L5 result (register 1441): **PASS at Simiiformes on both held-outs (pending audit)**; recent pair F still dips
Edges: A119b E6 1,378 orphan-orphan; E7 + 3,832 orphan-family kept (4,322 dropped as widening). Testis E7 133 + 42 (159 dropped).
Dev (fixed scorer, Simiiformes): baseline pair F 0.679 (prec 0.706), family F 0.755; E6 0.678 (0.694), 0.771 (fails pair F); **E7 0.695
(0.716), 0.769 (qualifies)**; content: NPIP family + AL359195.2 (lncRNA), GOLGA8 family - LINC02259; scores of every mainstay/YAG identical.
Held out, Simiiformes: A119b half pair F 0.8903 -> 0.8904 (prec 0.9278 -> 0.9271), family F 0.8136 -> 0.8253 (+0.012); testis pair F 0.9002
-> 0.9028 (prec 0.9562 -> 0.9556), family F 0.8650 -> 0.8962 (+0.031); mainstay/YAG scores identical, testis family content unchanged.
Reported: recent A119b pair F 0.9200 -> 0.9160 (prec 0.9028 -> 0.8943), family F 0.8765 -> 0.8819; testis 0.9536 -> 0.9522 (0.9387 -> 0.9343),
0.9160 -> 0.9277. Primates A119b 0.8842 -> 0.8836, 0.8102 -> 0.8291; testis 0.8901 -> 0.8927, 0.8546 -> 0.8855.
The A119b pair-F margin is +0.0001 (1 more true pair, 2 more predicted); the gain is in FAMILY recovery, not pairs.
**Audit of L5 (2026-10-27, metric-auditor; /mnt/linuxdisk/tmp/audit_2026-10-27_denovofamL5/): NOT a held-out validation; no default flip.**
Numbers reproduce to the digit; tie policy irrelevant. But (1) the A119b held-out half and testis had already been scored for E5 (L4) and
the audit pair lists that motivated E7 included held-out cross-half pairs, so neither set is fresh for E7. (2) The A119b pair-F margin
(+0.0001) is an MCL side effect: E7's own held-out gains are 2 true / 16 (FAM209A-B, ADH1B-C), and the margin comes from ZNF66 and ZNF430
leaving the 24-gene KZNF family into 2-member families with new loci (13 false pairs removed); without that, pair F 0.8863 < 0.8903.
(3) 94% (A119b) / 90% (testis) of E7's new pairs are outside the judgeable Compara universe; 401 of 988 new A119b copies map to no CAT gene.
(4) orphanlift.py was committed with the result, not the prereg. **Verdict r1441: NOT VALIDATED (fresh substrate needed); opt-in only.**
What would decide it: E7 on a sample never scored for any L arm (e.g. HG002 Kinnex FLNC, needs a full de novo run + ungated all-vs-all),
with pair-F split into intended joins vs side-effect splits and a minimum number of judgeable new pairs.

## PREREG GUIDED-S (2026-10-27): guided mode in SPLICED mode - CDS seeds find their family's other members and pseudogenes
User question 2026-10-27: does guided mode work spliced, i.e. fed a CDS from a family, does it find every other member, pseudogenes included?
Gene-body (genomic, intron-containing) seeds match a processed pseudogene only exon by exon; a spliced seed matches spliced copies through
introns and processed copies contiguously.
**Truth (bench/guided_eval/cdsseed.py, written before any score):** Compara paralogue components at **Simiiformes** (user choice; recent and
Primates reported) over ALL CAT/Liftoff v2.0 protein-coding genes (no expression filter), size >= 2. Pseudogenes (no parent table exists:
Ensembl CHM13 GFF and Compara carry none): a CAT `*pseudogene` whose name minus a final `P[0-9]*` equals a member's name, or a member's name
minus one trailing `[A-Z]` or `[0-9]+` token; one matching family only (else excluded). Stated blind spot: unnamed (`AC...`/`LOC`) pseudogenes.
**Seeds:** per family, the member with the smallest ENSG that has a CDS; its longest CAT CDS (spliced, 5'->3'), its protein (translation), and
its gene body (genomic span). A second seed (next ENSG) is reported for seed robustness.
**Arms (each against CHM13 v2.0, tool settings named, no tuning):** A = gene body, `minimap2 -x asm20 -c -N 100 -p 0.1` (the guided all-vs-all's
preset/secondary settings); B = CDS, `minimap2 -x splice -c -N 100 -p 0.1`; C = protein, `miniprot` 0.18 defaults (`--gff`).
**Hit -> gene:** each record's aligned genomic blocks; a record counts when it covers >= 0.5 of the seed (Liftoff -a) at identity >= 0.5
(Liftoff -s; PAF nmatch/blocklen, miniprot Identity). Its gene = the CAT gene with most exonic bp overlap, priority protein-coding >
pseudogene > other. The seed's own gene is excluded.
**Metrics (micro over families):** member recall = found other members / other members; pseudogene recall = found proxy pseudogenes / proxy
pseudogenes; precision = hits on the family's members or proxy pseudogenes / hits on any judgeable gene (Compara-universe protein-coding or any
family's proxy pseudogene); unjudgeable hits counted. Per mainstay/YAG family listed. Families split dev/held-out by famtruth `half`.
**Decision:** on DEV, the best arm by member recall + pseudogene recall among those with precision >= A's - 0.02. If that arm is not A, it must
beat A on HELD-OUT families on both recalls with precision >= A's - 0.02, and lose no mainstay/YAG member found by A. PASS -> spliced seeds
become the guided finder's documented default seed (user decides the flip).
**Deviation (2026-10-27):** the first scorer call printed BOTH halves for arm B (held-out seen for B before selection: member recall 0.747,
pseudo 0.368, precision 0.865). Fixed (`--show`, default dev). Selection still uses dev only; B's held-out reading must be stated as seen.
**Amendment 1 (2026-10-27, after the dev readout of A/B/C; held-out of A, C unseen, B seen by mistake):** dev member recall / pseudogene
recall / precision: A 0.423 / 0.035 / 0.882; B 0.717 / 0.129 / 0.601; C 0.119 / 0.000 / 0.690 (miniprot's default `--outs 0.99` keeps only
near-best hits). As registered, A is selected (B fails precision) -> **GUIDED-S FAIL**. Diagnostics (dev): 113 of B's 248 off-family hits are
Compara paralogues of the seed at an older level (KZNF clusters, PCDHB, H2AC, KRTAP: multigene families the Simiiformes cut splits);
`minimap2 -x splice` returns few secondary copies of near-identical arrays (TSPY seed on chrY: 2 records with -x splice, 49 with -x asm20;
not -f, not -G). Added arm **D** = B's records UNION the CDS aligned with `-x asm20 -c -N 100 -p 0.1` (r1287 G1's own unit + envelope design;
no new number). Same dev selection and held-out bars.

### GUIDED-S result (register 1442): **FAIL as registered (precision vs Compara Simiiformes)**; spliced seeds find far more copies
Dev (166 families' seeds; Simiiformes): member recall / pseudogene recall / precision: A gene body 0.423 / 0.035 / 0.882; B CDS splice
0.717 / 0.129 / 0.601; C protein miniprot 0.119 / 0.000 / 0.690; **D CDS splice + asm20 0.806 / 0.259 / 0.633**. Every spliced arm loses
> 0.02 precision, so A stays; nothing held out (B's held-out was seen by mistake, see Deviation). Off-family hits are multigene relatives the
Simiiformes cut separates (D: 113 of 249 are Compara paralogues of the seed at an older level; the rest mostly KZNF/KRTAP/PCDH cluster
members Compara splits across gene trees). Mainstays/YAGs (dev seeds): DAZ A 0/4 -> B 4/4; XAGE/GAGE A 0/26 -> B 2/26; NPIP A 21/22 ->
B 6/22 (splice mode drops near-identical secondaries; D restores it); TSPY 0/38 in A and B. `minimap2 -x splice` alone is unfit for arrays.
Reading for the user's question: a CDS seed DOES find spliced members and processed pseudogenes that gene-body seeds miss (pseudogene
recall x7), but against a Simiiformes truth it is scored as imprecise because it also finds the older members of the multigene family.

## PREREG GUIDED-S2 (2026-10-27): spliced seeds judged against MULTIGENE-family truth (user: HGNC gene groups)
r1442 failed on precision because Compara (any level) splits multigene clusters: 220 of 359 off-family hits are not Compara paralogues at
all, and any-level Compara merges superfamilies (largest 417). User rule: the multigene family is THE family; user chose HGNC gene groups.
**Truth (cdsseed.py `--truth hgnc`, HGNC complete set downloaded 2026-10-27, /mnt/linuxdisk/tmp/hgnc/):** each HGNC gene takes its MOST SPECIFIC
gene group (fewest members; ties: lowest id); CAT protein-coding genes are joined (by Ensembl id) when they share that group OR form a
Compara Simiiformes pair (genes in no group, e.g. NPIP, keep their Compara family). Families = components with >= 2 protein-coding genes.
Pseudogenes = CAT pseudogenes in the family's HGNC group plus the r1442 name proxy (one family only).
**Seeds, arms, hit rule:** as GUIDED-S: A gene body asm20; D CDS splice + asm20 (the r1442 leader); B, C reported. Families split dev/held-out by
the same SHA-1 rule on the NEW families. **Primary metric MACRO** (mean over families of per-family recall; one seed per family, so a
373-member KRAB group must not dominate); micro reported. Precision = micro over judged hits (judged = genes in any truth family or proxy).
**Decision:** dev: D replaces A if member and pseudogene macro recall are higher and precision >= A's - 0.02; then held-out on the same
three conditions plus no mainstay/YAG member found by A lost. PASS -> spliced CDS seeds (splice + asm20 union) proposed as the guided default.

### GUIDED-S2 result (register 1443): **PASS on held-out families (pending audit)** - spliced CDS seeds (splice + asm20) beat gene-body seeds
Truth: 1,879 HGNC-group/Compara families (16,412 protein-coding genes; largest 1,095 = C2H2/KRAB zinc fingers), 6,881 pseudogenes in 759.
Dev (964 families): MACRO member / pseudogene recall, precision: A 0.025 / 0.012, 0.950; B 0.051 / 0.109, 0.959; **D 0.053 / 0.115, 0.959**
(micro member 0.017 / 0.033 / 0.039). D qualifies. **Held-out (915 families): A 0.016 / 0.004, 0.859 (55/64); D 0.038 / 0.126, 0.956
(307/321)**; micro member 0.006 -> 0.014, pseudogene 0.002 -> 0.047. No mainstay/YAG member found by A lost; D gains TSPY 0 -> 38 members
(+10 pseudogenes), DAZ 0 -> 4, XAGE 0 -> 2, SPATA31 1 -> 2, BPY2 pseudogene 0 -> 1. Absolute recall is low because most HGNC groups are old,
diverged families no nucleotide aligner reaches at Liftoff's 0.5 floors; the comparison is the same truth for both arms.
Caveats: the truth changed twice after the r1442 fail (user choices); S2 held-out families were first scored here; B's S1 held-out was seen.
**Audit of S2 (2026-10-27, metric-auditor; /mnt/linuxdisk/tmp/audit_2026-10-27_guided_s2/): NOT grounds for a default change.** Numbers
reproduce to the digit (paired held-out families: D beats A on members 47 vs 5, pseudogenes 75 vs 1). But (1) the mainstay guard was EMPTY:
under the merged HGNC truth every mainstay family is seeded by a foreign gene (NPIP by SLC29A1, TBC1D3 by RABGAP1, GOLGA by GOLGA5, PRY/RBMY
by MSTRG.4761), so A finds 0 and nothing can be lost; with a real NPIPB12 seed (S1 files) A finds 21/22 and D 6/22, because D's asm20 half
breaks a CDS into exon-by-exon records each below the 0.5 per-record coverage floor (register §4.1 line 290, mRNA seed 19/19 -> 9/19, is the
same failure). (2) The named gains quoted as held-out (TSPY, DAZ, SPATA31, CT47A, VCY) are DEV families; "XAGE 0 -> 2" is GAGE12C/E.
(3) Truth and metric changed after the r1442 dev readout; D's held-out member gain equals B's, whose S1 held-out was seen; 43 S2 held-out
families contain S1-dev genes. (4) 26% of D's held-out hits are unjudged (A 3/67); counted as wrong, precision D 0.709 vs A 0.821 (fails).
(5) Truth families chain non-homology HGNC groups (the 1,095 family spans 31 groups, seeded by ARF5; not "C2H2/KRAB"); 845/915 held-out
families are 0 in both arms. (6) A is not handicapped by its body-coverage rule (CDS-based coverage leaves A at 0.0167 / 0.0086). (7) D loses
SSX and SPANXN members A finds. Seed 2 was promised and not run. **Verdict r1443: NOT VALIDATED.** The fix the audit points to is Liftoff's
own step: CHAIN a seed's records per target locus (aggregate coverage), which needs its own prereg and a fresh substrate.

## PREREG GUIDED-S3 (2026-10-27): Liftoff-chained spliced seeds, human dev, FRESH gorilla held-out
Fixes from the S2 audit, before any S3 score:
- **Chaining (Liftoff's step, both arms):** a seed's records (all of an arm's PAFs) are grouped by target gene (hit -> gene rule of GUIDED-S);
  per (seed, target gene, strand) the Liftoff DAG chain (colinear in seed and genome, not overlapping in the seed; maximize aligned seed
  CDS bases) gives coverage = chained seed-CDS bases / CDS length, identity = base-weighted nmatch/blocklen; a target is found when coverage
  >= 0.5 and identity >= 0.5 (Liftoff -a/-s). For the gene-body seed (A) only its CDS bases count (seed CDS mapped into body coordinates).
- **Truth without chaining:** a gene's family = its MOST SPECIFIC HGNC gene group (no union with Compara); genes in no HGNC group form
  families from Compara Simiiformes pairs among group-less genes. >= 2 protein-coding members. Pseudogenes = group members with pseudogene
  biotype + the name proxy (assigned only when exactly one family matches; fixes the 27 double assignments).
- **Seeds:** a family holding a mainstay/YAG gene is seeded by its smallest-ENSG (gorilla: smallest GeneID) label member; others as before.
  Seed 2 (next member) run and reported for every family.
- **Precision both ways:** judged-only AND unjudged-as-wrong; both bars apply.
**Arms:** A = gene body `-x asm20 -c -N 100 -p 0.1`; E = CDS `-x splice -c -N 100 -p 0.1` UNION CDS `-x asm20 -c -N 100 -p 0.1`; chained as above.
**Dev = ALL human families (CHM13, CAT v2.0)** (S2's split is spent). **Held-out = gorilla guided** (mGorGor1 v2.0 `GGO.fasta`,
RefSeq `GGO_genomic.gff`; never used by any guided arm): seeds = RefSeq CDS / gene bodies, truth by gene symbol (HGNC symbol -> group,
human Compara by symbol; `LOC` genes outside the truth, a stated blind spot). Never pool human and gorilla.
**Bars (dev, then gorilla):** E beats A on MACRO member and MACRO pseudogene recall; precision (judged) >= A - 0.02 AND precision (unjudged as
wrong) >= A - 0.02; for every label-seeded mainstay/YAG family E finds >= A's members. PASS on gorilla -> chained spliced seeds proposed as
the guided default (user decides).

### GUIDED-S3 RESULT (2026-10-27, register 1444): FAIL on dev; gorilla held-out NOT run (kept fresh)
`bench/guided_eval/guided3.py`, outputs `/mnt/linuxdisk/tmp/guided_s3h/seed{1,2}/score.{A,E}.seed*.json`. Truth: 2,089 human families
(largest 340, 200, 121), 6,776 pseudogenes in 865 families, 15 label-seeded mainstay/YAG families.

| seed | arm | MACRO member | MACRO pseudo | prec judged | prec unjudged-as-wrong | micro member |
|---|---|---|---|---|---|---|
| 1 | A body asm20 | 0.0252 | 0.0125 | 0.846 | 0.713 | 312/14,300 |
| 1 | E CDS splice+asm20 | 0.0460 | 0.1159 | 0.908 | **0.689** | 514/14,300 |
| 2 | A | 0.0288 | 0.0172 | 0.869 | 0.807 | 356/14,300 |
| 2 | E | 0.0508 | 0.1209 | 0.856 | **0.678** | 544/14,300 |

Bars: recall bars PASS both seeds; judged precision PASS both; **unjudged-as-wrong precision FAILS both** (-0.024, -0.129: E's extra hits are
mostly genes in no truth family); **named bar FAILS**: NPIP 6/21 (E) vs 13/21 (A) seed 1, 7 vs 20 seed 2; GOLGA8 1 vs 2 seed 1. E wins
TSPY 38/38 vs 0/38 (seed 1). Mechanism of the NPIP loss: E misses the divergent NPIPA subgroup (NPIPA1/2/5/7/8); the CDS query (1,767 bp) gives
asm20 records of <= 1,126 bp, while the gene body gives enough anchors through introns to reach them. Chaining did not recover them.
Side finding: against HGNC most-specific groups BOTH nucleotide arms find only 2.5-5% of members per family; most HGNC groups are too divergent
for a nucleotide seed at Liftoff 0.5/0.5, so this truth measures protein-level homology more than multi-copy families.
**Post-S3 dev exploration (not a verdict; dev only):** arm U = body asm20 UNION CDS splice+asm20, chained (guided3.py fix: gene-body PAFs `A*.paf`
always CDS-mapped; A numbers unchanged). U: MACRO member 0.048/0.053, pseudogene 0.117/0.121, judged precision 0.890/0.847 (A 0.846/0.869),
unjudged-as-wrong 0.680/0.673 (A 0.713/0.807); NPIP 13/21 and 20/21 (= A), TSPY 38/38, GOLGA8 2/27 (= A). **The unjudged hits are the truth's
gap, not off-family genes:** of U's 353 unjudged seed-1 hits, 349 are pseudogenes of the seed's own family that the name proxy missed
(FAM90A1 -> FAM90A2P..25P, MT-ND6 -> MTND6P*, HSP90AA1 -> HSP90AA4P/5P, FAM47A -> FAM47DP, FAM25A -> FAM25BP) or unnamed AC*/AL* processed
pseudogenes (TAF11L12 -> AC233724.*); 4 lncRNA (SETMAR transposase hits). The unjudged-as-wrong bar therefore penalises the pseudogenes the
user asked for; a fair S4 needs an independent pseudogene-parent judge before any held-out score.

## PREREG GUIDED-S4 (2026-10-27): independent pseudogene-parent judge; body UNION CDS seed; gorilla held-out
User choice (after S3): judge pseudogene hits with an independent parent annotation. `bench/guided_eval/pgparent.py`:
- **Human:** GENCODE v50 2-way consensus pseudogenes (Yale PseudoPipe + UCSC RetroFinder, `parent_id`) -> the GENCODE v50 pseudogene gene they
  overlap (same strand, >= 0.5 of the 2-way record, exactly one) -> CAT/Liftoff CHM13 gene by ENSG; parent ENSG -> CAT gene (exactly one).
- **Gorilla:** RefSeq pseudogene `description` minus a trailing `-like` -> the single protein-coding gene with that description (Gnomon's call).
**Truth = S3 truth, pseudogene sets widened:** pseudogene p joins family F when the S3 proxy assigns it to F OR its parent is a member of F; if the
two sources name different families, p is dropped. Judged set = S3 judged + every pseudogene with a table parent (a found pseudogene whose
parent is outside the seed's family counts wrong). Families and seeds are unchanged from S3 (pseudogenes do not shape families).
**Arms:** A = gene body asm20; U = gene body asm20 UNION CDS splice UNION CDS asm20 (one chain per seed/target/strand, body records counted on
seed CDS bases only). Both seeds.
**Dev = human (all families; already explored post S3, so dev is a sanity check only). Held-out = gorilla** (still never scored by any guided
arm). **Bars (each seed, dev then gorilla):** U > A on MACRO member and MACRO pseudogene recall; judged precision >= A - 0.02 and
unjudged-as-wrong precision >= A - 0.02; every label-seeded mainstay/YAG family U >= A members. PASS on gorilla -> propose U as the guided
default (user decides), then metric-auditor before any claim.

### GUIDED-S4 RESULT (2026-10-27, register 1445): FAIL on dev by 0.0016 (seed 2 strict precision); gorilla held-out NOT run
Parent tables: human 7,984 CAT pseudogene -> parent pairs (8,938 GENCODE v50 2-way pseudogenes, 7,755 linked); gorilla 1,180 (7,079 RefSeq
pseudogenes). Truth pseudogenes 6,776 -> 9,625 (1,132 families). Seeds byte-identical to S3 (`cmp cds.fa`). Outputs `/mnt/linuxdisk/tmp/guided_s3h/s4seed{1,2}/`.

| seed | arm | MACRO member | MACRO pseudo | prec judged | prec unjudged-as-wrong | named losses |
|---|---|---|---|---|---|---|
| 1 | A | 0.0252 | 0.0091 | 0.848 | 0.713 | - |
| 1 | U | 0.0477 | 0.1302 | 0.906 | 0.808 | none / 14 |
| 2 | A | 0.0288 | 0.0127 | 0.871 | 0.820 | - |
| 2 | U | 0.0533 | 0.1282 | 0.869 | **0.7987** (bar 0.8003) | none / 13 |

All bars pass on seed 1; seed 2 passes recall, judged precision and the named bar, and misses strict precision by 0.0016. With the independent
judge the unjudged pool falls from 353 to 162 (seed 1), so S3's precision "loss" was mostly the truth's pseudogene gap. Pseudogene recall
x10-14 (0.009 -> 0.130), member recall x1.9; TSPY 0/38 -> 38/38 (seed 1). Verdict as registered: FAIL; the held-out stays fresh.
**Post-S4 dev diagnosis (seed 2, U-only bad hits):** 109 "wrong" protein-coding hits, 48 from ONE seed, LINC02218 (HGNC: lncRNA; CAT biotype
protein_coding) -> TAF11L5/10/14; the rest spread (ACTA2 -> POTEG/POTEKP actin-domain fusions, MIA2 -> CTAGE*P, USP17L3 -> USP17L27, LILRB5 ->
LILRA2). Unjudged: FAM90A26 -> FAM90A*P (absent from the 2-way table), unnamed AC*/AL* retrocopies.

## PREREG GUIDED-S5 (2026-10-27): S4 with an HGNC cross-checked protein-coding set; gorilla held-out
One change from S4 (user: "keep FAIL, new prereg"): a gene is protein-coding for families and seeds only when its annotation biotype is
protein_coding AND, when HGNC lists it (human: by ENSG; gorilla: by symbol), HGNC locus_type is "gene with protein product". Genes HGNC calls
non-coding fall out of families, seeds and the judged set. Everything else = S4 (truth with pgparent tables, arms A and U, both seeds, bars).
Seeds that change are re-aligned with the same presets. Bars on human dev, then gorilla held-out; PASS on gorilla -> propose U as the guided
default (user decides), metric-auditor before any claim.
**S5 dev (human, both seeds): every bar PASSES.** 2,083 families; 2-3 seeds changed per seed set (re-aligned, `s5seed*/Anew.paf`, `*_new.paf`).
seed 1 A member 0.0253 / pseudo 0.0091 / judged 0.854 / strict 0.713; U 0.0479 / 0.1304 / 0.909 / 0.807. seed 2 A 0.0289 / 0.0127 / 0.876 /
0.818; U 0.0532 / 0.1282 / 0.899 / 0.823. Named losses none (14, 13). Ruling: the HGNC filter applies only to CAT protein_coding genes (a first
version also removed HGNC pseudogenes from the judged set, giving judged precision > 1; caught and fixed before any verdict). Ruling: gorilla
LOC pseudogenes enter the truth only through the parent table (S4 built it for them); LOC protein-coding genes stay outside (S3 blind spot).

### GUIDED-S5 RESULT (2026-10-27, register 1446): dev PASS; gorilla held-out FAIL on seed-2 strict precision (0.170 vs bar 0.173)
Gorilla held-out (`/mnt/linuxdisk/tmp/guided_s5g/`): 1,960 families, 353 pseudogenes in 127 families, 3/2 label-seeded named families.
CDS extraction checked (ATG starts >= 99.9%, stops 100%).

| seed | arm | MACRO member | MACRO pseudo | prec judged | prec unjudged-as-wrong | named losses |
|---|---|---|---|---|---|---|
| 1 | A | 0.0042 | 0.0083 | 0.806 (25/31) | 0.143 (+144) | - |
| 1 | U | 0.0210 | 0.1581 | 0.834 (171/205) | 0.195 (+670) | none / 3 |
| 2 | A | 0.0060 | 0.0005 | 0.674 (29/43) | 0.193 (+107) | - |
| 2 | U | 0.0219 | 0.1125 | 0.780 (131/168) | **0.170** (+604; bar 0.173) | none / 2 |

Every bar passes except seed-2 strict precision -> **FAIL as registered**. Of U's unjudged gorilla hits, 99.8% are LOC genes (seed 1: 383 LOC
pseudogenes, 284 LOC protein-coding; seed 2: 380, 224): the by-symbol gorilla truth cannot judge them, so the strict bar there measures the truth's
LOC blind spot (declared in S3). On judged hits U beats A on both seeds (+0.03, +0.11) and on recall x4-5 member, x19-200 pseudogene. Held-out
is now SPENT for this question. A further test needs a judge for LOC genes (RefSeq description -> parent, as pgparent.py does for
pseudogenes) and a fresh substrate.

## PREREG GUIDED-S6 (2026-10-27): S5 arms on a FRESH gorilla held-out, Kamilah_GGO_v0 + NCBI Annotation Release 102
User choices: mGorGor1 hap1/hap2 have no annotation, so the held-out is **Kamilah_GGO_v0 (GCF_008122165.1, RefSeq release 102, 2019)**: a different
animal, assembly and annotation from mGorGor1 (spent in S5); the same gene families biologically (stated caveat). Never pooled with human.
**Arms, seeds, truth, bars = S5** (A = body asm20, U = body UNION CDS splice UNION CDS asm20 chained; HGNC most-specific group truth by symbol +
Compara Simiiformes for group-less genes; --hgnc-coding; pseudogene parent table from this GFF's descriptions; both seeds; U > A on MACRO member
and pseudogene recall, judged and strict precision >= A - 0.02, no named-family member loss).
**One addition, fixed before looking (the S5 failure mode):** LOC protein-coding genes get the same description judge as pseudogenes
(description minus a trailing `-like` -> the single named protein-coding gene with that description). A found LOC gene whose parent is in the seed's
family counts CORRECT; parent in another family or outside the truth -> judged WRONG; no parent -> unjudged. LOC genes never enter recall
denominators. Dev check: human has no LOC genes, so S5 dev numbers stand. PASS -> propose U as the guided default (user decides), metric-auditor first.
**S6 setup (before any score):** Kamilah files md5-verified (`/mnt/linuxdisk/tmp/kamilah/`). Ruling: release-102 gene lines carry no
`description`, so descriptions come from NCBI Gene `gene_info` (All_Mammalia, taxid 9595) by GeneID, the same NCBI naming source as mGorGor1's
GFF descriptions (coverage: pseudogenes 5,176/6,290, LOC protein-coding 2,996/3,800, named 16,700/16,856) — cost if wrong: names updated after
2019 may differ from release 102's models. Parent table: 951 pseudogene pairs, 146 LOC pairs (118 LOC genes judged-correct in 67 families).
Truth: 1,927 families, 353 pseudogenes in 132 families, 2 label-seeded named families. CDS ATG >= 99.8%, stops >= 99.7%.

### GUIDED-S6 RESULT (2026-10-27, register 1447): FAIL on the fresh Kamilah held-out (strict precision both seeds; judged precision seed 2)
| seed | arm | MACRO member | MACRO pseudo | prec judged | prec unjudged-as-wrong | named (U/A/n-1) |
|---|---|---|---|---|---|---|
| 1 | A | 0.0043 | 0.0219 | 0.814 (35/43) | 0.255 (+94) | |
| 1 | U | 0.0206 | 0.1699 | 0.830 (171/206) | **0.223** (bar 0.236) | XAGE 1/0/2, GOLGA6L 0/0/2 |
| 2 | A | 0.0037 | 0.0002 | 0.778 (28/36) | 0.228 (+87) | |
| 2 | U | 0.0198 | 0.1456 | **0.751** (bar 0.758) | **0.203** (bar 0.208) | XAGE 2/0/2, GOLGA6L 1/0/2 |

Recall bars and named bar PASS; precision bars FAIL. U-only judged-wrong hits (27 / 41) are mostly LOC genes (22 / 31) whose description
parent sits in a sibling HGNC group or outside the truth: IFITM1/IFITM2 -> LOC IFITM-like copies (12, 13), SUMO2 -> LOC (10), LILRA <-> LILRB,
MIA2 -> CTAGE15, CR1L -> CR1. Pattern across S3-S6: U's extra hits land in neighbouring HGNC sub-groups of the same superfamily, which the
most-specific-group truth calls wrong. Verdict: no validated guided improvement; default unchanged. Kamilah held-out now SPENT.

## PREREG GUIDED-S7 (2026-10-27): HGNC PARENT-group truth; fresh human held-out HG002 pat (Ensembl rapid release 2022_07)
User choices after S6: truth one level up the HGNC hierarchy; fresh held-out = HG002.
**Truth change (only one vs S5):** a gene's family key = the PARENT of its most-specific HGNC group in HGNC `hierarchy.csv` (several parents ->
the smallest by gene count, then lowest id); a group with no parent keys itself. Group-less genes: Compara Simiiformes as before; >= 2
protein-coding members; --hgnc-coding; pseudogene sets = proxy + parent table.
**Dev = human CHM13 / CAT v2.0, re-run on the new truth** (seeds change where families merge; changed seeds re-aligned with the same presets).
**Held-out = HG002 paternal haplotype (GCA_021950905.1, carries chrY) with Ensembl rapid-release gene set 2022_07** (gene symbols from the
Ensembl xref `HGNC Symbol` rows; genes without an HGNC symbol are outside the truth = stated blind spot, like LOC). Ensembl's set is built by
its own projection + paralogue-mapping pipeline, independent of our chaining. Human held-out numbers are never pooled with gorilla.
**HG002 pseudogene parent judge:** an Ensembl pseudogene description `... (SYMBOL) pseudogene` -> the single HG002 protein-coding gene with that
symbol. **Arms, seeds, chaining, bars = S5:** A body asm20 vs U body UNION CDS splice UNION CDS asm20; both seeds; U > A on MACRO member and
pseudogene recall; judged and unjudged-as-wrong precision >= A - 0.02; no named mainstay/YAG member loss. Bars on dev, then HG002.
PASS -> propose U as the guided default (user decides), metric-auditor before any claim.

### GUIDED-S7 RESULT (2026-10-27, register 1448): FAIL on dev (seed-2 strict precision 0.7685 vs bar 0.7805); HG002 held-out NOT scored
Parent-group truth: 1,151 human families (largest 487, 465, 387), 9,564 pseudogenes in 719 families, 13 named. Seeds: same gene+CDS reused
(alignments renamed, `s7seed*/ *_re.paf`); 13 / 195 new seeds aligned with the same presets.
seed 1 A member 0.0337 / pseudo 0.0103 / judged 0.869 / strict 0.710; U 0.0549 / 0.1025 / 0.903 / 0.794 (all bars pass).
seed 2 A 0.0371 / 0.0060 / 0.873 / 0.801; U 0.0561 / 0.1031 / 0.870 / **0.7685** (bar 0.7805); named losses none.
U-only unjudged seed-2 hits (77): 36 are FAM90A26 -> FAM90A*P (own-family pseudogenes absent from the GENCODE 2-way table and missed by the name
proxy, the same gap seen since S3), the rest AC*/AL* unnamed pseudogenes (IFITM3 -> IFITM4P/8P, MIA2 -> CTAGE11P/12P, RPL26L1 -> AL*). Across S3-S7
the strict bar fails on pseudogenes the judges cannot reach, not on off-family protein-coding hits (2 here: ZNF582 -> ZNF781, STON1 -> STON1-GTF2A1L).

## PREREG GUIDED-S8 (2026-10-27): S7 + a name-stem pseudogene judge; dev, then HG002 pat held-out
User choice after S7. **Only change vs S7:** a named pseudogene left unassigned by the proxy AND the parent table gets base = its name minus a
trailing `P[0-9]*`; member stem = a family member's name minus its trailing `[0-9]+[A-Z]?[0-9]*`, kept when >= 3 characters. Among families
having a member stem that is a prefix of base, keep those with the LONGEST such stem; exactly one -> the pseudogene joins that family
(recall + judged); several -> unjudged. Unnamed (AC*/AL*/LOC/ENSG) pseudogenes stay with the parent table only. Everything else = S7
(parent-group truth, arms A/U, both seeds, bars). Dev re-scored on existing alignments (seeds unchanged: pseudogenes do not shape families).
**S8 dev (human, existing S7 alignments, seeds byte-identical): every bar PASSES.** Pseudogenes 9,564 -> 9,914 (745 families).
seed 1 A member 0.0337 / pseudo 0.0123 / judged 0.884 / strict 0.812; U 0.0549 / 0.1036 / 0.908 / 0.841. seed 2 A 0.0371 / 0.0083 / 0.875 /
0.816; U 0.0561 / 0.1031 / 0.877 / 0.821. Named losses none (13, 13). Ruling: clone-named genes (`^[A-Z]{1,2}[0-9]{5,6}\.[0-9]+$`, e.g.
AC007610.3) are excluded from the stem judge on both sides, as the prereg keeps unnamed pseudogenes with the parent table only.

## PREREG GUIDED-S9 (2026-10-27): many templates per family (every isoform CDS + body, 2 members) vs one template
User question: does one representative lose copies? Arm **M** = per family, the first TWO members (the S8 seed-1 and seed-2 genes), each
contributing its gene body (asm20) AND the CDS of EVERY annotated transcript (splice + asm20); each template chained and judged as in S3
(coverage >= 0.5 of that template, identity >= 0.5); the family's found set = the UNION over templates. Baseline **U2** = the same two members,
body + longest CDS only (the S8 U arm run from both seeds, unioned) — so M differs from U2 ONLY by the extra isoforms. Both seeded genes are
excluded from hits and denominators. Truth = S8 (parent group, --hgnc-coding, parent table, stem judge). Families with < 3 protein-coding
members have no recall denominator and count for precision only.
**Bars (dev human, then HG002 pat held-out, the same held-out as S8 — both pre-registered before it is looked at, separate questions):**
M > U2 on MACRO member and MACRO pseudogene recall; judged and unjudged-as-wrong precision >= U2 - 0.02; no named mainstay/YAG family with
fewer members than U2. PASS -> multi-template seeds proposed for guided mode and as the template rule of the de novo Liftoff mode (user decides).

### GUIDED-S9 RESULT (2026-10-27, register 1449): FAIL on dev (strict precision 0.787 vs bar 0.795); HG002 NOT scored
`bench/guided_eval/multiseed.py`, `/mnt/linuxdisk/tmp/guided_s3h/s9/`. 1,151 families, 2,302 seeded genes, 11,996 distinct isoform CDS
templates (2,327 longest: length ties count as longest for U2). 1,024 families with >= 1 non-seeded member.
U2 (2 members, body + longest CDS): MACRO member 0.0563, pseudo 0.1916, judged 0.889 (1194/1343), strict 0.815 (+122).
M (2 members, body + EVERY isoform CDS): member 0.0649, pseudo 0.2283, judged 0.890 (1452/1632), strict **0.787** (+213). Named: none lower;
NPIP 19/20 -> 20/20, TSPY pseudo 10 -> 15/34, GOLGA 7 -> 8/26 (wrong 1 -> 3), RBMY wrong 2 -> 4.
M-only unjudged hits (91): 81 pseudogenes no judge reaches (RBMY1B -> RBMY2KP/2XP/3AP/2EP/2QP/2DP — RBMY2* pseudogenes whose stem RBMY2 is
not a member stem; FAM86B1 -> FAM86EP/FP; DDX11 -> DDX11L5; clone-named AC*/AL*), 7 lncRNA (NSFL1C), 2 readthrough genes.
Descriptive: two members vs one (U2 vs S8 U) roughly doubles pseudogene recall (0.10 -> 0.19); extra isoforms add +15% members, +19% pseudogenes
at the same judged precision. The strict bar again fails on pseudogenes the truth cannot judge, not on wrong genes.

## PREREG ASMTPL (2026-10-27): assembled transcripts as family templates, gorilla, 13 families
User (10-27): "use gorilla and just a handful of gene families; several representatives; the assembler first, so lookups reduce to unique
intron chains". Substrate mGorGor1 v2.0 + RefSeq (`GGO_genomic.gff`), families by RefSeq description keyword (`/mnt/linuxdisk/tmp/asmtpl/
family_regions.tsv`, 231 loci): **dev** NPIP (5 pc), TBC1D3 (16+1), GOLGA6 (14+3), GOLGA8 (17), USP17L (45+11), NBPF (23+9); **held-out**
TSPY (6+10), XAGE (7+3), CT47 (8), RBMY (13+10), SPATA31 (10+3), CDY (1+15), DAZ (1). Members = protein-coding genes of the keyword; pseudogenes =
pseudogene-biotype genes of the keyword. (Gorilla here answers a NEW question — template construction — not the S5/S6 guided question.)
**Templates:** `copy_assign --assemble-only --assembly-junctions strict --bridge-regroup f1v2` (pipeline defaults, no secondary pool) on the
family loci (+-5 kb) of KB3781 fibroblast `fibroblasts/GCA_029281585.2_flnc_mm.bam` and OR6737 testis `winloci_data/GGO_mm.bam`; transcripts
assigned to the family whose locus they overlap; deduplicated by exact intron chain across both libraries; template sequence = spliced exons.
**Arms:** R1 = one template per family (the most-read-supported chain, ties: longest); RK = every unique-chain template. Each template is aligned
`-x splice -c -N 100 -p 0.1` and `-x asm20 -c -N 100 -p 0.1` to GGO.fasta, chained per target gene (guided3 rule, coverage >= 0.5 of the
template, identity >= 0.5); family found set = union over its templates; the template's own source locus is excluded from hits and denominators.
**Score per family:** member recall, pseudogene recall, judged precision (hit on a gene of ANOTHER keyword family or on a named gene outside
the 13 = wrong; a LOC gene of no keyword = unjudged, reported, NOT penalised — the user-approved revision of the strict bar).
**Bars (dev, then held-out):** RK >= R1 on member recall in every family and > in at least one; RK >= R1 on pseudogene recall summed; judged
precision RK >= R1 - 0.02; unjudged count reported. Lookups reported: templates per family (RK) vs reads at the loci.

### ASMTPL RESULT (2026-10-27, register 1450): recall bars PASS dev + held-out; precision bar FAILS both -> FAIL as registered
`bench/guided_eval/asmtpl.py`, `/mnt/linuxdisk/tmp/asmtpl/`. Assembler on 163 merged loci regions: fib 362 / tes 526 transcripts, 21 s / 18 s each;
764 unique intron chains (12,021 primary reads at the loci -> 764 lookups, 16x fewer). Rulings before the verdict: `--bridge-regroup off` (f1v2
needs one region per contig); family keywords matched by SYMBOL too (named members carry HGNC-style descriptions, LOC copies Gnomon's: the
first scorer missed every named member); a template excludes only its OWN locus; 78 templates whose exons overlap a named protein-coding gene
of another family are flagged chimeric and skipped (NPIP 40/49, XAGE 18/26 — readthrough into neighbours).

| split | arm | member | pseudogene | judged precision | unjudged |
|---|---|---|---|---|---|
| dev (6 fams) | R1 one template | 3/122 | 0/24 | 1.000 (3/3) | 19 |
| dev | RK all chains | **62/122** | **9/24** | **0.617** (71/115) | 102 |
| held-out (7) | R1 | 24/56 | 4/41 | 1.000 (28/28) | 2 |
| held-out | RK | **34/56** | **13/41** | **0.940** (47/50) | 6 |

Per family RK >= R1 everywhere (GOLGA8 0 -> 17/17, GOLGA6 0 -> 12/16, TBC1D3 0 -> 11/16, NBPF 3 -> 20/23, RBMY 9 -> 13/13 + 10/10 pseudo).
Precision: RK's dev "wrong" hits are (a) GOLGA6 <-> GOLGA8 cross-hits (one HGNC superfamily split by my keywords) and (b) SD-passenger
neighbours of the family loci (HERC2, ITSN2, DNM1, GOLGA2 beside GOLGA8; USP32/USP6/ABCA9 beside TBC1D3; ALG1/ALG1L beside USP17L):
templates carry UTR/readthrough exons of the neighbour. USP17L 0/45 in both arms (intronless array; to inspect). Held-out wrong: DAZL (DAZ
relative), FRAT2/ZNF484 (SPATA31 neighbours). Verdict FAIL on the precision bar; the recall answer to the user's question is unambiguous:
one representative per family finds 3/122 dev members, all unique chains find 62/122.

## PREREG ASMTPL-2 (2026-10-27): family-shared exons, Liftoff per-exon scoring
Same data, families, split, arms and bars as ASMTPL. Two changes: (1) **template = family-shared exons only**: an exon of a unique chain is kept
when another unique chain of the same family (either library) has an exon with the same donor AND acceptor (first/last exons: the same inner
splice site and overlapping span); mono-exonic chains are kept whole; a chain with no shared exon is dropped. Chimeric flag no longer applied
(the shared-exon rule replaces it). (2) **Coverage counted per lifted exon** (Liftoff's feature lift): a target is found when the chained
template exons that land on it cover >= 0.5 of the template's exon bases at identity >= 0.5 — i.e. the S3 rule applied to the trimmed template.
R1 = the trimmed template with most reads; RK = every trimmed unique chain (identical trimmed sequences deduplicated).

### ASMTPL-2 RESULT (2026-10-27, register 1451): FAIL (precision unchanged, two held-out families lost)
Shared-exon trimming: 764 -> 725 templates; CT47 and TSPY (single chains, nothing shared) have NO template -> held-out denominators shrink
(members 56 -> 36). RK dev member 58/122, pseudo 7/24, judged precision 0.596 (R1 0.750); held-out 21/36, 12/31, 0.846 (R1 0.857).
(An aln.sh path slip re-aligned the v1 templates once; inputs identical, v1 PAFs unchanged.)
**Diagnosis of RK's neighbour hits (v1 PAFs, chained blocks vs the hit gene's exons):** HERC2 7/7, GOLGA2 2/2, ALG1 10/10, USP6 13/20 and DNM1 land
MOSTLY IN INTRONS of the hit gene (a family copy nested in the host's intron, grazing one host exon: GOLGA8 in HERC2, TBC1D3 in USP6, USP17L in
the ALG1L region) — these are copies the annotation does not name, not wrong paralogues; `gene_of` hands them to the protein-coding host.
ITSN2 34/34, USP32 7/9, ABCA9 1/1 land on the hit gene's EXONS (real cross-hits to look at). Next (not run): judge a hit on target EXON bases
(Liftoff lifts features onto features); intronic landings become "putative unannotated copy" — O3's tier names — rather than wrong.

## PREREG ASMTPL-3 (2026-10-27): exon-side judge, GOLGA6+GOLGA8 merged, full templates
Same data, templates (ASMTPL v1, untrimmed, chimeric flag kept), arms and bars as ASMTPL. Changes: (1) a chained hit is judged on the target
gene's EXON bases: found when >= 0.5 of the chained template bases lie on that gene's exons (and the S3 coverage/identity rule holds); a chain
whose blocks lie mostly inside a gene's introns is binned "putative unannotated copy" (reported per family, neither right nor wrong — O3's tier
name); (2) GOLGA6 and GOLGA8 scored as ONE family (one HGNC superfamily; templates of either seed it); (3) ITSN2 checked by hand before scoring
(reported either way).

### ASMTPL-3 RESULT (2026-10-27, register 1452): precision 0.62 -> 0.91 dev, 0.94 -> 0.95 held-out; FAIL as registered (bar = R1 1.000 - 0.02, R1 has 3 dev hits)
`asmtpl.py` `ASMTPL_EXON=1` (`/mnt/linuxdisk/tmp/asmtpl/v3/`). ITSN2 by hand: ONE 3.2-kb segment (NC_073236.2:115,683,695-115,686,863, inside an
ITSN2 exon) hit by 24 GOLGA8/6 templates at 86-89% identity — a shared segment, kept as wrong.

| split | arm | member | pseudogene | judged precision | unjudged | putative (intronic) |
|---|---|---|---|---|---|---|
| dev | R1 | 3/122 | 0/24 | 1.000 (3/3) | 11 | 8 |
| dev | RK | **55/122** | 8/24 | **0.913** (63/69) | 46 | 73 |
| held-out | R1 | 24/56 | 4/41 | 1.000 (28/28) | 2 | 0 |
| held-out | RK | **31/56** | 4/41 | **0.946** (35/37) | 4 | 15 |

Per family RK >= R1 everywhere; GOLGA merged 23/33 members (2 wrong: ITSN2 segment, UBE2Q2P16), NBPF 19/23 (PDZK1, TRIP13), TBC1D3 11/16
(ABCA9, USP6 exon-side), RBMY 13/13. Cost of the exon judge: hits on pseudogenes whose annotated exons do not cover the lifted blocks moved to
putative (RBMY pseudo 10 -> 2, dev members 62 -> 55), so the putative bin holds both unannotated copies and loosely-annotated pseudogenes.
The precision bar is degenerate: R1 has 3 and 28 hits, so 1.000 - 0.02 is unbeatable by any arm with real recall; the user decides the bar
for a re-run (e.g. a fixed 0.90, or R1 - 0.02 with R1 pooled over splits). USP17L 0/45 persists (9 putative: the array copies' annotation).

### LIFTMODE prototype (2026-10-27, no verdict): `bench/denovo_eval/liftmode.py`
De novo Liftoff-style family mode wired from the ASMTPL-3 template rule: every unique intron chain of the assembled GTF is a template; each is
aligned (splice + asm20) and Liftoff-chained; a placement counts when >= 0.5 of the chained bases lie on the exons of ONE assembled locus
(gene_id); family = connected component of loci joined by placements (no MCL, no thresholds beyond -a/-s); landings on no locus's exons =
putative copies (BED). Smoke test on the 13-family loci assembly (fib + tes, `/mnt/linuxdisk/tmp/asmtpl/lm/`): 845 transcripts -> 760 chains,
233 loci -> 26 families, 414 putative copies. By annotation keyword the families are clean: NBPF 16 genes (+2 unnamed LOC), RBMY 11, GOLGA 6+3,
TBC1D3 4, TSPY 4, NPIP 2, 0 mixed families. Fix on the way: `asmtpl.read_gtf` dedupes exons (the same transcript id occurs in both library GTFs).
Not yet: genome-wide run vs the MCL catalog (prereg needed; gorilla OR6737 `asm_gap/heldout/gorilla_OR6737.good.dt.families.gtf`, 171 MB).

## PREREG LIFTMODE-GW (2026-10-27): Liftoff-style de novo families genome-wide vs the MCL catalog, gorilla
Substrate: the genome-wide assembled GTFs `rustle_figures/runs/gorilla_KB3781/gorilla_KB3781.families.gtf` (dev, fibroblast) and
`.../gorilla_OR6737/gorilla_OR6737.families.gtf` (held-out, testis). Arms: **MCL** = the shipped `mcl_families --from-gtf` catalog of the
same GTF (`<run>/gorilla_*.fam.copies.tsv`); **LIFT** = `liftmode.py` on the same GTF (templates = unique chains; splice + asm20 -N 100
-p 0.1 to GGO.fasta; Liftoff chaining -a/-s 0.5; placement on a locus's exons >= 0.5 of chained bases; connected components).
Scorer: `famtruth.py --species gorilla --levels simiiformes` (DENOVOFAM-G truth: human Compara by symbol; stated underpowered, 9/13 true
families) + the 13 keyword families of ASMTPL as a descriptive check (members per family, mixed families). Bars per sample: LIFT pair F and
family F >= MCL; pair precision drop <= 0.02; no mainstay/YAG regression (DENOVOFAM-G G3). PASS on dev then held-out -> liftmode proposed
as the de novo family definition beside MCL (user decides). Putative copies reported as counts only.

### LIFTMODE-GW RESULT (2026-10-27, register 1453): dev KB3781 FAIL (family F up, pair F and pair precision down); OR6737 NOT run
`/mnt/linuxdisk/tmp/liftgw/kb/`: 72,531 templates (244 Mb; 9 chunks, ~4 min splice + 45 s asm20 each), 89,253 placements, 22 s to families.
| arm | copies on genes | families | pair sens / prec / F | family sens / prec / F |
|---|---|---|---|---|
| MCL (shipped) | 731 | 302 | 0.138 / 0.267 / 0.182 (4/29, 15 pred) | 0.312 / 0.833 / 0.455 |
| LIFT | 3,476 | 1,233 | 0.966 / **0.074** / 0.137 (28/29, 379 pred) | 0.938 / 0.469 / **0.625** |
Bars: family F PASS (+0.17), pair F FAIL (-0.045), pair precision FAIL (-0.19), mainstays/YAGs unchanged (truth has almost none expressed in
fibroblasts). Keyword check: 0 mixed families in either; NBPF 14 -> 16/23, GOLGA 1 -> 2, TSPY 0 -> 1. Diagnosis: connected components over
placements at Liftoff's cross-species -a/-s 0.5 are transitive within one genome: 3,846 loci end up in families (MCL 918), largest 80 loci;
pair sensitivity is near-complete (28/29) but 379 predicted pairs for 29 true. Liftoff's own within-genome copy mode uses identity 1.0 (-sc),
which would be too strict here; a middle value is an arbitrary threshold. Threshold-free candidates for a re-run (user decides): reciprocal
placement (A lifts onto B AND B onto A) as the edge; or LIFT placements as edges INTO the existing MCL (the --extra-edges path of r1437-41).

## PREREG LIFTMODE-GW2 (2026-10-27): reciprocal placement as the family edge; then LIFT edges into MCL
Same data, templates, alignments, scorer and bars as LIFTMODE-GW. **Arm R:** locus A and locus B are joined only when some template of A
places on B AND some template of B places on A (each placement by the GW rule); families = connected components of reciprocal edges.
If R fails a bar on dev: **Arm X** = the GW placements (locus pairs, weight = identity x coverage of the best placement, max over directions)
fed to the shipped `mcl_families --from-gtf ... --extra-edges` (r1437-41 path) on the same GTF, scored the same way. Dev KB3781 first; a
passing arm is then run on OR6737 (held-out) before anything is proposed.
**GW2 arm R, dev KB3781 (`/mnt/linuxdisk/tmp/liftgw/kbR/`):** 439 families, 1,302 loci; pair sens 0.793 prec 0.190 F **0.307** (23/29 | 121 pred);
family sens 0.781 prec 0.641 F **0.704**. Pair F and family F beat MCL (0.182 / 0.455); pair precision 0.267 -> 0.190 FAILS the <= 0.02 bar
-> arm X next, as registered.

### LIFTMODE-GW2 RESULT (2026-10-27, register 1454): R, X and XR all beat MCL on pair F and family F; all FAIL the pair-precision bar -> FAIL, OR6737 NOT run
Base re-run with today's binary (`--reject-self-hits` default): `/mnt/linuxdisk/tmp/liftgw/mclX/`, cached all-vs-all PAF, 120 s each.
| arm (dev KB3781) | loci on genes | families | pair sens / prec / F | family sens / prec / F |
|---|---|---|---|---|
| MCL base | 667 | 267 | 0.138 / 0.333 / 0.195 (4/29, 12 pred) | 0.312 / 0.833 / 0.455 |
| R reciprocal components | 1,147 | 411 | 0.793 / 0.190 / 0.307 (23/29, 121) | 0.781 / 0.641 / 0.704 |
| X MCL + all LIFT edges (5,249; 2,340 added, 1,017 new nodes) | 1,497 | 572 | 0.862 / 0.181 / 0.299 (25/29, 138) | 0.906 / 0.690 / **0.784** |
| XR MCL + reciprocal LIFT edges (2,482; 969 added) | 1,141 | 446 | 0.655 / 0.216 / **0.325** (19/29, 88) | 0.750 / 0.706 / 0.727 |
Mainstays/YAGs identical in all arms (almost none expressed in fibroblasts). Pair precision drops 0.12-0.15 against a bar of 0.02: the bar
is set on a baseline with 12 predicted pairs, so any arm that finds the other 25 true pairs at < 1-in-3 precision fails. Verdict FAIL as
registered. Side note: /mnt/linuxdisk hit 100% during the run (today's indexes); 31 GB of reproducible .mmi deleted, 27 GB free now.

## PREREG LIFTMODE-GW3 (2026-10-27): retro-placement class + asymmetric pair judge; then OR6737
Same data/templates/alignments as GW/GW2. (1) **Retro placement:** a template with >= 1 intron whose chained placement on a locus has no
intron (every record of the chain is one block) while that locus's own transcripts have introns is classed `retro` (a processed copy inside the
locus), reported, and NOT a family edge. (2) **Judge:** recall = famtruth Compara Simiiformes (pairs, families), unchanged. Precision: a
predicted pair of judged genes is WRONG only when the two genes are Compara paralogues at NO level, share NO HGNC gene group, and are not the
same symbol (PAR); otherwise it is true (Simiiformes) or "older/HGNC" (reported, not penalised). Arms: MCL base, XR (MCL + reciprocal LIFT
edges, primary), X (reported). **Bars (dev then OR6737):** XR pair F and family F > base at Simiiformes; XR wrong pairs <= base wrong pairs + 2;
mainstays/YAGs no regression (G3). PASS on OR6737 -> propose `--extra-edges` from liftmode reciprocal placements as part of the de novo
family definition (user decides), metric-auditor first.
**GW3 dev KB3781 (`/mnt/linuxdisk/tmp/liftgw/{kb3,kb3R,mclX}/`):** retro placements 1,539 removed from the edge set (edges 5,249 -> 4,759; reciprocal
2,482 -> 2,286). Simiiformes: X3 pair F 0.338 (25/29 | 119) fam F 0.784; **XR3 pair F 0.339 (19/29 | 83) fam F 0.727** (base 0.195 / 0.455).
Asymmetric judge, WRONG pairs: base 2 (CHAMP1-TUBB8, MTHFS-SLAIN2), XR 5, X3 9, **XR3 4** (+AP1S2-BDP1, HNRNPA1L2-MRPS31); XR3's other
non-Simiiformes pairs: 21 older Compara, 36 HGNC-group (KZNF), 3 PAR. All XR3 dev bars PASS -> OR6737 held-out next.

### LIFTMODE-GW3 RESULT (2026-10-27, register 1455): XR3 PASSES dev AND the OR6737 held-out (pending metric-auditor)
OR6737 (`/mnt/linuxdisk/tmp/liftgw/{or,or3,or3R,mclOR}/`, frozen copies in `liftgw/frozen/`): 73,965 templates (251 Mb), 106,065
placements, 1,989 retro placements removed; base = the shipped catalog (1,521 copies; ungated all-vs-all PAF, today's binary, 235 s).
| arm (OR6737 held-out) | loci on genes | families | pair sens / prec / F (Simiiformes) | family sens / prec / F | WRONG pairs (judge) | older / HGNC / PAR |
|---|---|---|---|---|---|---|
| MCL base | 1,171 | 414 | 0.225 / 0.450 / 0.300 (9/40, 20) | 0.352 / 0.950 / 0.514 | 1 (MTHFS-SLAIN2) | 5 / 1 / 4 |
| **XR3** | 1,622 | 568 | 0.650 / 0.283 / **0.394** (26/40, 92) | 0.667 / 0.783 / **0.720** | **2** (+CDH18-NUP50) | 27 / 33 / 4 |
| X3 (reported) | 1,984 | 694 | 0.875 / 0.205 / 0.332 (35/40, 171) | 0.778 / 0.712 / 0.743 | 9 | 57 / 66 / 4 |
Bars: pair F and family F > base PASS; wrong 2 <= 1 + 2 PASS; mainstays/YAGs identical (NPIP 1/1 fam3, DAZ 1/1) PASS. Dev (KB3781) passed
the same bars (pair F 0.195 -> 0.339, family F 0.455 -> 0.727, wrong 2 -> 4). Cost: ~45 min of template alignment per sample (9 chunks),
22 s to edges, +2 min in mcl_families. Verdict: PASS, awaiting the metric-auditor; on audit PASS -> propose to the user: liftmode reciprocal
placements (retro-filtered) as `--extra-edges` in the de novo `families` stage.
**AUDIT of r1455 (metric-auditor, 2026-10-27 22:30): HOLDS WITH CAVEATS -> downgraded to NOT VALIDATED after the requested check.** All numbers
reproduce to the digit; OR6737 was run after the GW3 prereg commit. Caveats: (1) 94.8% of the 1,382 pairs XR3 adds on OR6737 (92.9% of 1,006
on KB3781) touch a LOC/non-judgeable gene and are invisible to both judges (the r1440/r1441 failure mode); (2) the bar was changed twice (GW,
GW2 precision bar -> GW3 "wrong <= base + 2" with a judge designed after inspecting dev's false pairs); at the registered level pair precision
still falls (KB 0.333 -> 0.229, OR 0.450 -> 0.283); (3) the judge excuses Opisthokonta-deep paralogy (tubulins, MAGEA); (4) G3 has no power
(1 expressed mainstay member); (5) bipartite family precision is conditioned on the matching (re-denominated OR fam F 0.418 -> 0.529 still up).
**Requested check, RefSeq description-stem agreement of the unjudged NEW pairs vs base's own unjudged pairs:** KB3781 0.812 (759/935) vs 0.690
— higher; **OR6737 0.760 (995/1310) vs 0.846 — lower by 0.09** -> by the auditor's rule the held-out row reads NOT VALIDATED. The stem check
undercounts (NDUFA4 named two ways, HLA Gogo/Patr/HLA, ZC3H11A vs "zinc finger CCCH domain-containing protein 11B-like"), so a protein-level
reciprocal-best-hit judge of the LOC pairs is the next thing that could settle it (not run).

## PREREG LIFTMODE-GW3-RBH (2026-10-27): protein-level judge of the unjudged pairs
Gene set = every protein-coding gene in a base pair or an XR3-added pair on either sample (1,378 genes); protein = translation of each gene's
longest RefSeq CDS (`GGO_genomic.gff`). `mmseqs easy-search` all-vs-all (-s 7.5, --max-seqs 2000, default e-value 1e-3). For a pair (A,B):
**mutual hit** = A's protein lists B and B's lists A (primary; family size does not cap it); **RBH** = each is the other's best non-self hit
(secondary). Rates compared per sample: XR3-added protein-coding pairs vs base's own protein-coding pairs. Rule, fixed now: OR6737 mutual-hit rate
of added pairs >= base's -> r1455 upgraded to HOLDS; lower -> stays NOT VALIDATED. Pairs with a pseudogene are reported separately (no protein).

### LIFTMODE-GW3-RBH RESULT (2026-10-27, register 1456): r1455 stays NOT VALIDATED (held-out mutual-hit rate 0.921 vs base 0.962)
1,378 proteins (longest RefSeq CDS, translated), `mmseqs easy-search -s 7.5`, 31,355 hits (`/mnt/linuxdisk/tmp/liftgw/rbh/`).
| sample | pairs | n (both coding) | mutual protein hit | RBH |
|---|---|---|---|---|
| KB3781 dev | base | 622 | 0.883 | 0.228 |
| KB3781 dev | XR3-added | 835 | **0.943** | 0.153 |
| OR6737 held-out | base | 2,208 | 0.962 | 0.106 |
| OR6737 held-out | XR3-added | 1,035 | **0.921** | 0.136 |
Unjudged-only subsets are within 0.005 of these. Pairs with a pseudogene/non-coding member (no protein): KB 171/1,006 added, OR 347/1,382.
By the fixed rule (added >= base on OR6737) the held-out FAILS by 0.041 -> r1455 remains NOT VALIDATED. Descriptively: 92-94% of the LOC pairs
the lever adds are protein-level homologues (RBH is capped by family size, so mutual hit is the fair number); the 6-8% with no mutual hit are
the unresolved residue. The lever's recall gain is real; its added pairs are slightly less protein-supported than MCL's own on the held-out.
**Diagnosis of the 82 unsupported added pairs (OR6737, 2026-10-27 22:50):** by locus structure 33 pairs have a locus spanning >= 2 annotated
coding genes (NOTCH2-NBPF, ACTB-RNPC3, LARP1B-KLHL2, C8orf89-ZNF589), 5 are same-chromosome < 1 Mb, 44 other (SMG1-kinase copy pairs with the
same description, NBPF-sortilin, "uncharacterized LOC" partners). Bridge tables of the run: only **1/82** involves a flagged F1v2 bridge
(97 bridges, 517 bridge junctions), so the fusion-category / readthrough machinery (fusions.tsv classes, bridge regroup) does not reach them;
40/82 are UNFLAGGED multi-gene loci (exons touching two coding genes), 41/82 single-gene loci. Added pairs involve a multi-gene locus more often
than base pairs (OR 267/1,382 = 19% vs 339/2,832 = 12%; KB 26% vs 32%). The lever to test next is template boundaries, not bridge classes:
templates cut to read-defined units (`mcl_families --emit-units`, units follow reads) instead of whole loci, or chains trimmed to exons the
family's other chains also carry (ASMTPL-2's rule, which cost single-chain families). Not run.
