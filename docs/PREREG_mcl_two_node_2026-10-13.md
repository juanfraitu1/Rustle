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
- **F2:** Task 2
