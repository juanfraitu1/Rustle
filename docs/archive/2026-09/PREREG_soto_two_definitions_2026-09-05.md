# PREREG — the two O1 definitions on ONE truth: Soto 2025 families, HUMAN A119b/CHM13 v2.0 (2026-09-05)

Old definition (E_r / γ-quasi-clique, `gw_family_catalog --homology-primary`, ledger §6bx, `soto_adj/arm_off`):
detection any-overlap 293/362; pair precision in the [0.90,1.00) band **149/153 = 0.9739**; pooled (≥50 %)
precision 0.7198, recall|both detected 0.8743, recall|all Soto pairs 0.1729 (any-overlap: 0.4006 / 0.6866 / 0.5238).
New definition (MCL on the annotation, `mcl_families` at the canonical defaults: identity ≥0.70, cov_longer
≥0.30, ≥300 bp, exonic floor 1 bp, inflation 2.8, prune 1e-9, locus rule ON).
Paired substrate = the SAME slice the old arm saw: the 362 Soto member regions (`soto_adj/regions.bed`); the
annotation's 748 genes/pseudogenes overlapping them; reads `soto.bam` (the old arm's `-M -L` subset) for units.
⚠ Both definitions on a slice: only Soto's neighbourhoods exist, so false merges with the rest of the genome
cannot occur in either arm — the numbers are paired, not absolute. Same scorer `bench/soto_adjudicate.py`,
same three units, same bands (max within-family identity of the member).
Predictions:
- P1 pair precision in [0.90,1.00): MCL ≥ 0.95 (old 0.9739). Below 0.90 Soto is silent; not scored.
- P2 recall | both detected (≥50 % floor): MCL ≥ 0.87 (old 0.8743).
- P3 detection ≥50 %: MCL > 0.80 (old 0.3536) — the annotation is the node, so this is expected and is NOT
  credited to the definition; the recall|ALL-Soto-pairs difference that follows is a node effect.
Decision rule: P1 ∧ P2 ⟹ the switch costs nothing on the advisor's benchmark and gains recall through the
node; report both and keep the old definition runnable. ¬P1 ⟹ the switch is not warranted on precision;
report it and stop. Family-level exact-set match is reported, not predicted.
