# PREREG — MCL two-node components: the fix (2026-10-13)

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

## Results

- **F1 PASS.** Red (`cargo test --release --lib -- a_connected_pair a_pair_joins`, before the fix): `a_connected_pair_is_one_family_at_every_weight` failed at `w = 0.99` with `left: [[0], [1]]`; `a_pair_joins_inside_a_mixed_graph_and_the_rest_is_untouched` failed with `left: [[0, 1, 2], [3], [4], [5]]`. Green (after): these two and `the_pair_rule_touches_only_two_node_components` pass, `3 passed`; whole lib `1051 passed; 0 failed; 13 ignored`.
- **F3 PASS.** denovo_intact/guided_intact (denovo_F/guided_F), before -> after, HEAD binary vs fixed binary: human SNRPB `snp_0.01` 0/0 (0.0000/0.0000) -> 1/1 (1.0000/1.0000); `snp_0.1` 0/0 (0.0000/0.0000) -> 1/1 (1.0000/1.0000); `shuffle` 0/1 (0.0000/1.0000) -> 1/1 (1.0000/1.0000); orangutan SPAG7 `snp_0.01` 0/0 (0.0000/0.0000) -> 1/1 (1.0000/1.0000); `snp_0.1` 0/0 (0.0000/0.0000) -> 1/1 (1.0000/1.0000). The before state equals the 10-02 tables; the only arm already intact before the fix is human `shuffle` guided (1, as on 10-02).
- **F4 PASS.** `cargo test --release --no-fail-fast`: 1290 passed, 3 failed, 14 ignored; the 3 failures are the known `tests/pipeline_driver_dispatch.rs` tests (`all_runs_assemble_families_assign_flag_and_neither_catalog_nor_candidates`, `all_with_candidates_runs_the_candidates_stage_between_families_and_assign`, `all_with_legacy_catalog_runs_the_catalog_stage_between_families_and_assign`; each `left: 127, right: 0`).
- **F2 PASS.** Each catalog re-clustered from its cached `.fam.loci.paf` / `.fam.loci.gff3` (`bench/o1_hier/recluster.sh`; outputs, `<sample>.graph.tsv` and `<sample>.pair_diff.json` in `/mnt/linuxdisk/tmp/o1_close/recluster/`), pre-fix binary (md5 9914dee365a041347275f2fafba25e07) vs fixed binary (md5 db55a19f8f165587563212c276bcfbca). **H0 PASS for all four:** the pre-fix binary reproduces each shipped `.fam.clusters.tsv` and `.fam.loci.tsv` byte-identically (the shipped files were only read). `bench/o1_hier/pair_diff.py` (member sets): `pairs_split` 0 and `other_changes` 0 in all four, so every pre-fix cluster is unchanged; every joined pair is a graph edge (no `no_weight`). Per catalog: pairs joined (`tie_regime` / `below_w_star`), clusters before -> after, members before -> after, `uniform_cliques_below_1` (components / members).
  - human_A119b: 68 joined (0 / 68); clusters 2150 -> 2218; members 9907 -> 10043; uniform cliques below 1: 3 / 9.
  - human_testis: 9 joined (0 / 9); clusters 342 -> 351; members 1115 -> 1133; uniform cliques below 1: 0 / 0.
  - chimp_PTR: 20 joined (0 / 20); clusters 396 -> 416; members 1343 -> 1383; uniform cliques below 1: 0 / 0.
  - gorilla_KB3781: 15 joined (0 / 15); clusters 326 -> 341; members 918 -> 948; uniform cliques below 1: 2 / 6.
- **F2, what the joins are (post hoc, descriptive; not a bar).** Every joined pair lies below w*: 112 of 112, none in the tie regime. All 7,256 two-node components of the four graphs were classified (A119b 4,634, testis 689, chimp 858, gorilla 1,075). The 6,969 with w >= w* (4,449 / 678 / 815 / 1,027) were already one family or one locus before the fix, so in these catalogs the residue never split a pair. The 287 below w* (185 / 11 / 43 / 48) were split before and none is split after. So here the fix is Amendment 1(a)'s definition change, not a residue repair. `loci.tsv` changes too: the fold-within-clusters step (§6ey) merges records that overlap on exon bases and share a cluster, so a joined pair of such records becomes ONE locus and leaves `clusters.tsv`. That is 117 / 2 / 23 / 33 new `loci.tsv` lines: 175 of the 287 pairs fold into one locus, and the other 112 stay listed as 2-member families. Every new line lies inside a below-w* two-node component; no pre-fix `loci.tsv` line is removed or changed. The 5 uniform cliques below 1 (3 K3 in A119b at w 0.994380, 0.999527, 0.999607; 2 K3 in gorilla at 0.989512, 0.999411) all came out whole (one 3-member family, or one locus after the fold).

**Amendment 1 (after F1/F3/F4, before F2; 2026-10-13).** Task 1 found a boundary below which the rule changes the definition, and a premise that was too strong. Committed before the F2 re-clustering runs; nothing above `## Results` is edited.

- **(a) Weight bands.** The pre-registered rule joins every connected pair. Below w* = 0.456885 (inflation 2.8, self-loops 1.0) exact MCL splits a lone pair deterministically, so those joins are a definition change, not a residue fix. Above w* the column converges to a tie, and float residue decided join or split: that is the bug. F2 reports joined pairs by band: `tie_regime` (w >= w*) and `below_w_star`. The choice between the weight-blind rule and exact-MCL semantics is the user's; the weight-blind rule stays the default until then.
- **(b) "Three or more nodes join" is too strong.** In an independent port, a uniform triangle with one weight w < 1 splits asymmetrically at 39 of the 500 weights 0.500..0.999, and a uniform K4 at 58. Such components are not changed by this fix. F2 counts them in the four catalogs (`uniform_cliques_below_1`: components that are complete graphs of >= 3 nodes whose edges all carry one weight below 1, equal to 6 dp).
- **(c) Corrections to the Bug paragraph.** "The column converges toward a tie" holds only above w*. "160 of the 500 weights" is 159 in two independent Python ports (first split at 0.504); the 160 came from a third port, and no Rust sweep measured the pre-fix count.
