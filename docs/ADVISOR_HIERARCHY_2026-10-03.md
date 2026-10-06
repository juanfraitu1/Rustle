# "Where is my hierarchy in the MCL?" — what the advisor most likely means, and what exists (2026-10-03)

Written for the next meeting. Sources: `docs/seeded_family_definition.md` §0★/§0★★★, ledger §5k, §6l, §6n7, §6kl,
memory `project_family_hierarchy`, register rows 9, 311, 457, 654.

## The likely meaning
MCL returns ONE flat partition at one inflation value (shipped: 2.8). Multi-copy families are not flat: superfamily ⊃
family ⊃ subfamily (NPIPA / NPIPB, TBC1D3 groups) ⊃ copies ⊃ alleles. "Where is the hierarchy" = where in the output is
that nesting, since a single cut asserts one resolution and hides the others. It is the threshold objection restated for
clusters: a cut is a parameter; a tree with levels is a structure.

The "true positives / false positives" reading is a consequence of the same point: with a hierarchy, TP and FP are
level-relative. NPIPA+NPIPB merged is a false merge at the subfamily level and a true family one level up. Our data
already behaves this way: Soto's families sit between our L2 and L3 (§6n7); NPIP is one family at L1/L2 and splits at L3;
precision against any single-resolution truth caps near 0.8 (§6kl). A tree lets each truth be scored at its own level.

## What exists
| item | status | where |
|---|---|---|
| Hierarchy inside MCL | none; inflation is a resolution knob, clusterings at different inflations are not nested; ours stable I = 2.0–4.0, shatter at ≥ 6 | §6ec, §6ew, register 670 |
| Nested edge-test lattice L0–L4 | proposed design, nesting is a theorem (T1), not the default; L3 (0.985) / L4 (0.995) shipped as thresholds | `seeded_family_definition.md` §0★★★ |
| What breaks nesting | the leader rule (RNA ⊆ DNA), Louvain (T1/T1′) | memory `project_leader_rule_breaks_nesting`, §6n5 |
| Average-linkage tree on gated distance | recovers NPIPA/NPIPB exactly (label-pure k = 2..8) where NO threshold can (A/B boundary = 0.3 % identity window) | memory `project_family_hierarchy`, register 457 |
| What a tree does NOT add | merge heights = 1 − identity relabelled (partial r −0.05); lift over one cut = +0.035 ARI = two gene attachments; best_k oracle-selected | §6l |
| Hierarchy as catalog repair | NO-GO: no coarse level at which NPIP reunites without the blob; broad/recent two-level hierarchy reaches 19 % of families (70 % are 2-copy) | §5k, register 9, 311 |

## The answer to give
"The MCL cut is one level of a nested structure, not the structure. The principled object is the lattice: connected
components of one evidence graph under nested edge tests, so nesting is provable; L3/L4 are its shipped levels, and an
average-linkage dendrogram on the gated distance draws it. The tree is a presentation and a scoring frame (each truth at
its level), not extra information: merge heights carry nothing identity does not."

## Question to ask him first
"Do you mean a nested sequence of clusterings — family ⊃ subfamily ⊃ copy — with provable nesting, so a truth is scored
at its own level? Or a confidence ranking of clusters?" Both are served by the same levels, since each level is a
stricter edge test; confirm which before building anything new.

## Do not re-propose
A subfamily level as a repair for fragmentation (§5k), a broad/recent two-tier hierarchy (rows 9, 311), MCL inflation
sweeps as a hierarchy (not nested), superfamily components over the cluster graph (row 654).
