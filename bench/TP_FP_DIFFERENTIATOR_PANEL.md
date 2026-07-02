# TP-vs-FP Family-Call Differentiator Panel — DEFINITIVE

Consolidated, size-controlled, held-out sweep of every candidate differentiator between
**true** and **false** family calls, at two levels (E_r **edge** = pair-of-loci homology
call; family **block** = merged component), validated against the gold diploid CN oracle.
The sharp new hypothesis under test: over-merges are dense (scalar density fails) but may be
**bimodal two-clusters** detectable by sub-structure features (edge-weight bimodality,
Laplacian spectral gap, best-2-cut / Louvain modularity).

- Code: `bench/tp_fp_panel.py` (deterministic, `PYTHONHASHSEED=0`, seed=0; three reruns byte-identical).
- Data: `bench/tp_fp_panel_edge.tsv` (2189 edges × 45 feats), `bench/tp_fp_panel_block.tsv`
  (606 families × 63 feats), `bench/tp_fp_panel.json` (both stages + verdict).
- Oracle: `bench/diploid_cn_oracle.tsv` (6 DNA-confirmed FP blocks).

---

## ⭐ ONE-PARAGRAPH VERDICT

**No usable NEW differentiator surfaced. The only deployable winners remain the two already
shipped: `aln_frac` (homology-extent) at the edge level and VG shared-multiplicity
(`min_shared_mult` / repeat-hub gate) as the one orthogonal axis — plus `core_recip` as a
small partially-orthogonal edge adder.** On the hard TP-vs-genuine-FP edge contrast, no
RNA-deployable single feature or logistic combination beats `aln_frac` (held-out 0.830)
out-of-fold within-block; every feature above it is either a near-monotone transform of
`aln_frac` (ρ 0.85–0.98, lift within fold noise) or DNA-only (`id`/`min_cov`, sparse
n=468). At the block level, the sharp bimodal-over-merge hypothesis is **refuted as a
general rule but partially vindicated on one block**: the NEW spectral family
(`spec_ratio32` 0.722, `spec_best2cut_mod` 0.655) are the only genuinely size-independent
RNA differentiators (survive a quadratic size basis) yet they are confounded by real
paralog sub-clades and **do not beat the reused base-graph combo (perfold 0.723); adding
sub-structure HURTS (0.676)**. On the confirmed FPs: GSTM2 hub-block 9 *is* genuinely
bimodal (best-2-cut / Louvain top-1%/128, hypothesis TRUE) but is already caught by shipped
features; **MAGE-549 is a dense uniform blob (best-2-cut ≈ 0, hypothesis FALSE)**, flagged
only by `af_bimod` (identity-distribution bimodality, above the real-dense max) which does
NOT generalize (held-out AUC 0.42); GSTM2 triangle-block 13 is a 3-node clique where no
sub-structure can exist. **Audit-corrected concentration result:** the over-merge residual
is **NOT concentration-DNA-bound** — the RNA-deployable NEW feature `ew_bimod` reaches
**4/6 oracle FPs @K=50, tying DNA `id_bimod`** — but it **IS generalization-DNA-bound**:
`ew_bimod` fails size control (in-sample lin/quad AUC 0.348/0.318) and the clean
matched-real-dense test does not flag MAGE; **only DNA `id_bimod` both concentrates 4/6 and
generalizes held-out (0.76).** Bottom line: the missed over-merge residual is real,
partly a two-cluster phenomenon (best-2-cut sees GSTM2-9) and partly a dense-uniform
identity phenomenon (MAGE), and **no single size-controlled RNA feature catches all of it
as a deployable gate — it remains generalization-bound to DNA identity.**

---

## Setup & contrast reconciliation

- **Edge panel** = 416 TP + 1773 genuine-FP E_r edges (those carrying computed features).
  TP = `in_dna_loose` (DNA/cDNA genome truth), independent of the scored RNA
  `aln_frac`/`core_recip`. Held-out = 5-fold **GroupKFold by E_r connected component** (437
  components; edge-count-balanced folds [480,428,427,427,427]; the 480-edge giant block is
  isolated → no within-block leakage). Univariate metric = **train-oriented fold-mean AUC**
  (sign on train, evaluated on test, averaged) — immune to the pooled-reorient artifact.
- **Contrast note:** the known "aln_frac 0.915" is **TP vs ALL-FP**. On the harder
  **TP vs genuine-FP-only** contrast used here, aln_frac = **0.830 held-out / 0.858
  in-sample**. Truthbar-FP are easy (aln_frac 0.945); the genuine-FP are the residual.
- **Block panel** = 606 families; label `is_Ep_impure` (≥2 protein families by E_p, the
  over-merge PROXY, 80/606); gold `oracle_fp` (6 DNA-confirmed). Ranking = **OLS
  size-residual** on `1 + log(n_members) + log(distinct_loci)` fit on TRAIN, evaluated on
  TEST fold-mean; a **QUADRATIC size basis** (adds squares+cross term) flags nonlinear size
  artifacts (a linear-resid discriminator that collapses/sign-flips under quad = size).
  Blocks verified disjoint (0 DNs in >1 of 606) → random 5-fold is leakage-free.
- **Deployability:** DNA identity (`id`/`min_cov`, `id_*`) and soft-mask
  (`vg_softmask`/`mask_*`/`bridge_mask`, VALIDATION-flagged not RNA-pure),
  `frac_repeat_bridged`/`frac_ge_gate` are tagged **VALID-ONLY** and excluded from the
  deployable set. Everything else is RNA-deployable.

---

## (1) EDGE-LEVEL ranked held-out AUC — TP vs genuine-FP (hard contrast)

Ranked by held-out fold-mean AUC. `new` = introduced this session; `re-conf` = reproduced a
prior result on this contrast; `FALSIFIED-prior` = the earlier "falsified" verdict was on a
different (easier) contrast — see note.

| rank | feature | held-out | sd | in-samp | n | status | note |
|---|---|---|---|---|---|---|---|
| 1 | id | 0.869 | .08 | 0.892 | 468 | **VALID-ONLY (DNA)** | sparse; NOT deployable |
| 2 | min_cov | 0.861 | .05 | 0.852 | 468 | **VALID-ONLY (DNA)** | sparse |
| 3 | contain_x_alnfrac | 0.840 | .10 | 0.867 | 2189 | new | ρ=0.94 w/ aln_frac → redundant |
| 4 | core_x_alnfrac | 0.836 | .10 | 0.868 | 2189 | new | ρ=0.98 w/ aln_frac → redundant |
| 5–6 | contain_min / recip_contain | 0.833 | .11 | 0.854 | 2189 | re-conf | exon-containment (prior "FALSIFIED" = diff contrast) |
| **7** | **aln_frac** | **0.830** | **.09** | **0.858** | 2189 | **SHIPPED (reference)** | the workhorse |
| 8 | contain_max | 0.827 | .11 | 0.857 | 2189 | re-conf | |
| 9–10 | n_shared_nodes / log | 0.820 | .07 | 0.858 | 2189 | new | VG shared-node count; ties aln_frac, ρ=0.85 |
| 11–12 | aln_len / log_aln_len | 0.820 | .07 | 0.825 | 2189 | re-conf | |
| 13 | **core_recip** | 0.802 | .09 | 0.814 | 2189 | **SHIPPED** | only partially-orthogonal within-block adder |
| 14–16 | vg_softmask / mask_max / bridge_mask | 0.712 | .10 | 0.774 | 2189 | **VALID-ONLY** | re-confirms soft-mask ~0.69–0.77 |
| 17–18 | n_exon_ratio / n_exon_min | 0.70 / 0.68 | | | 2189 | re-conf | structure, weak |
| 19 | **min_shared_mult** | 0.675 | .08 | 0.685 | 2189 | **SHIPPED** | repeat gate; **ρ=−0.42 w/ aln_frac = the one orthogonal axis** |
| 20–26 | kmFrac100_* / n_exon_* / shares_only_repeat_m5 | 0.62–0.67 | | | | re-conf | weak |
| 27–41 | max_shared_mult, nKmers/kmMed, mmMax*, mult_ratio, both_single | 0.50–0.58 | | | | **FALSIFIED** | non-generalizing |
| 42–45 | mmMed_max/mean/min, nReads_mean/max | 0.38–0.48 | | | | **FALSIFIED** | read-degree/multimapping outright anti-generalize |

**Beats / ties / loses vs aln_frac (0.830):**
- **Beat (deployable):** *none meaningfully.* The only things above aln_frac are `id`/`min_cov`
  (DNA, VALID-ONLY, sparse) and three near-monotone transforms of aln_frac itself
  (`contain_x_alnfrac` ρ=0.98, `contain_min`, `core_x_alnfrac`); +0.003–0.010 is within the
  ±0.10 fold sd = noise.
- **Tie (redundant homology-extent cluster, ρ 0.85–0.98):** n_shared_nodes, aln_len,
  contain_max, core_recip — all 0.80–0.83.
- **Lose / dead:** soft-mask (VALID-ONLY), min_shared_mult (orthogonal but weak), kmer, and
  **all read-degree/multimapping (`mm*`) features — re-confirmed non-generalizing (fold-mean
  < 0.50).**

**Edge combinations (out-of-fold):** the pooled-OOF gain (full_deployable 0.878 vs
aln_frac 0.845) is a **cross-block calibration effect**, NOT better within-block separation —
under the deployment-relevant **per-fold-mean** it vanishes (full 0.831, +0.0002). The one
genuine within-fold adder is interpretable **`aln_frac + core_recip` = perfold 0.846 (+0.016)**;
adding `min_shared_mult` HURTS (perfold 0.816, −0.014).

---

## (2) BLOCK-LEVEL ranked size-residualized held-out AUC (is_Ep_impure)

`raw|lin|quad` = in-sample size-residual AUC. A feature is a **genuine** size-independent
differentiator only if it survives the **quadratic** basis (`quad` stays discriminative and
doesn't sign-flip vs `lin`). `n` = blocks where defined (spectral/bimodality need multi-node
/ enough edges).

| rank | feature | held-out | sd | raw/lin/quad | n | tag | verdict |
|---|---|---|---|---|---|---|---|
| 1 | id_bimod | 0.760 | .21 | .57/.59/.57 | 28 | **VALID-ONLY (DNA)** | only feature that generalizes |
| 2 | **spec_ratio32** (λ3/λ2) | 0.722 | .09 | .69/.73/.71 | 178 | **SUBSTRUCT-NEW** | ✓ survives quad (genuine, confounded) |
| 3 | spec_ratio23 (λ2/λ3) | 0.708 | .10 | .31/.30/.30 | 178 | SUBSTRUCT-NEW | ✓ |
| 4 | spec_nl2 | 0.690 | .08 | .27/.29/.27 | 178 | SUBSTRUCT-NEW | ✓ |
| 5 | spec_ngap23 | 0.686 | .07 | .67/.70/.71 | 178 | SUBSTRUCT-NEW | ✓ |
| 6 | spec_l2 (fiedler) | 0.684 | .06 | .27/.29/.32 | 178 | SUBSTRUCT-NEW | ✓ |
| 7–8 | n_edges / ew_n | 0.665 | .08 | .65/.34/.51 | 606 | **SIZE-ARTIFACT** | collapses under quad |
| 9 | mem_per_locus | 0.664 | .10 | .48/.67/.33 | 606 | **SIZE-ARTIFACT** | sign-flip |
| 10 | spec_l3 | 0.660 | .06 | .29/.33/.33 | 178 | SUBSTRUCT-NEW | ✓ |
| 12 | frac_ge_gate | 0.656 | .08 | | 606 | base-reused, **VALID-ONLY** | |
| 13 | **spec_best2cut_mod** | 0.655 | **.02** | .69/.68/.68 | 178 | **SUBSTRUCT-NEW** | ✓ most stable NEW; the flagship 2-cut stat |
| 14 | mean_degree | 0.649 | .09 | .64/.33/.44 | 606 | **SIZE-ARTIFACT** | |
| 15 | fiedler (=λ2) | 0.642 | .04 | .35/.35/.38 | 606 | base-reused | |
| 16 | mincut | 0.640 | .05 | | 606 | base-reused | |
| 19 | **ew_bimod** | 0.633 | .19 | .41/**.35/.32** | 93 | **SUBSTRUCT-NEW** | inverts under size-resid → not a gate |
| 22 | edge_density | 0.628 | .05 | .33/.36/.37 | 606 | base-reused | scalar density fails (as predicted) |
| — | max_mult / med_mult | 0.574 / 0.603 | | | 606 | base-reused (SHIPPED gate) | multiplicity |
| — | louvain_mod | 0.476 | .05 | .60/.50/.52 | 606 | SUBSTRUCT-NEW | null after size-control |
| — | af_bimod | 0.415 | .27 | .45/.46/.43 | 60 | SUBSTRUCT-NEW | flags MAGE but non-generalizing |
| — | cross_chrom_frac | 0.472 | .03 | | 606 | DISPERSION-NEW | null |

**Genuine size-independent RNA differentiators (survive quad, not DNA):** the **spectral
family** — `spec_ratio32`, `spec_ratio23`, `spec_nl2`, `spec_ngap23`, `spec_l2`, `spec_l3`,
`spec_best2cut_mod`, `fiedler` — all NEW, all defined only on the 178 multi-node blocks.
**Size-artifacts** (linear-resid discriminative, collapse/flip under quad): `n_edges`,
`ew_n`, `mem_per_locus`, `mean_degree`, `max_degree`, `frac_in_giant`, `id_skew`, `af_kurt`.

**Block combinations (5-fold OOF perfold):** base-graph (clustering+multiplicity+cardinality,
reused) = **0.723**; substruct-NEW alone = **0.554**; all-combined = **0.676**. →
**Sub-structure does NOT beat the base-graph combo — adding it HURTS.**

---

## (3) DIPLOID-ORACLE precision@K — concentrating the 6 confirmed FP blocks (AUDIT-CORRECTED)

In-sample-oriented size-residualized precision@K (sign + residualizer fit on all blocks;
**only `*_combo_oof` are true OOF**). Computed for **every** ranked feature (the earlier
top-8-only cut hid `ew_bimod`). n=6 oracle FPs → a 1-block (4/6 vs 3/6) gap is within noise;
read alongside the held-out AUC and the clean residual test below.

| scorer | K=10 | K=20 | K=30 | K=50 | deploy? |
|---|---|---|---|---|---|
| **id_bimod** | 1/6 | 3/6 | 3/6 | **4/6** | **VALID-ONLY (DNA)** — *and* generalizes (AUC 0.76) |
| **ew_bimod** | 0/6 | 1/6 | 3/6 | **4/6** | **RNA-DEPLOYABLE (NEW)** — ties DNA on raw concentration |
| spec_l2 | 1/6 | 2/6 | 2/6 | 3/6 | RNA |
| af_bimod | 1/6 | 1/6 | 2/6 | 3/6 | RNA |
| n_edges / ew_n | 1/6 | 2/6 | 3/6 | 3/6 | RNA (size-artifact) |
| spec_best2cut_mod | 0/6 | 2/6 | 2/6 | 2/6 | RNA |
| louvain_mod | 1/6 | 1/6 | 2/6 | 2/6 | RNA |
| SUBSTRUCT_combo_oof | 0/6 | 1/6 | 2/6 | 2/6 | combo (true OOF) |
| ALL_combo_oof | 0/6 | 0/6 | 0/6 | 0/6 | combo (true OOF) |
| spec_ratio32 / spec_ngap23 / spec_ratio23 | 0/6 | 0/6 | 0/6 | 1/6 | RNA |

**AUDIT-CORRECTED conclusion:** the over-merge residual is **NOT concentration-DNA-bound** —
the RNA-deployable NEW feature `ew_bimod` reaches **4/6@50, tying DNA `id_bimod`** (it catches
MAGE-549 at rank ~15/606). The earlier "only DNA concentrates / best RNA ≤3/6" claim was an
over-generalization from a top-8-only scan and is **retracted**. However the residual **IS
generalization-DNA-bound:** `ew_bimod`'s 4/6 is fragile — its in-sample size-resid AUC
**inverts** (lin 0.348 / quad 0.318), held-out is a noisy 0.633±0.188, and by the clean
matched-real-dense test it does **not** flag MAGE (pct 0.593). `af_bimod` cleanly flags MAGE
(pct 1.0) but held-out AUC 0.42. **Only `id_bimod` (DNA) both concentrates 4/6 and generalizes
held-out (0.76).** No RNA feature does both.

---

## (4) ⭐ THE SHARP RESIDUAL TEST — does any NEW sub-structure feature catch the GSTM2/MAGE residual, size-controlled?

Reference = 128 matched real-dense pure families; `pct` = fraction of real-dense scored below
the target (>0.9 = top decile = clean catch). All four genuinely-missed / confirmed FPs are now
covered (block 110 added per audit). Block-9 dominant gene is `LOC115930164` (oracle fam14);
the "GSTM2" naming is inherited from the domain-hub over-merge description — FP status is real.

| feature | b9 GSTM2-hub | b13 GSTM2-tri | b549 MAGE | b110 missed | outcome |
|---|---|---|---|---|---|
| **spec_best2cut_mod** | 0.317 (**.99**) | −0.093 (.39) | −0.003 (.79) | −0.001 (.80) | catches **b9 only** (top 1%) |
| **louvain_mod** | 0.320 (**.99**) | 0.000 (.05) | 0.001 (.77) | 0.244 (.98) | catches b9 (+ b110) — not MAGE |
| ew_bimod | 0.615 (.76) | na | 0.490 (.59) | 0.675 (.85) | b9/b110 mid; **misses MAGE** |
| **af_bimod** | 0.507 (.64) | na | **0.908 (1.00)** | na | catches **MAGE only** (> real-dense max) |
| spec_l2 / spec_l3 / spec_gap23 | .48/.70/.90 | .46 | **.92/.99/.98** | .02/.00/.03 | flag MAGE — but = density×weight proxy (corr 0.87) |
| spec_ratio32 | 2.53 (.73) | 2.09 (.65) | 1.70 (.53) | 1.13 (.23) | none clean |

**Per-target definitive outcome:**

- **GSTM2 hub-block 9 — hypothesis TRUE.** Genuinely bimodal two-clusters:
  `best_2cut_modularity` = 0.317, `louvain_mod` = 0.320 → **top 1% of 128 real-dense** (only
  1/128 real-dense reaches 0.30). The flagship "2 weakly-joined clusters" hypothesis holds
  **for this block.** Caveat: block 9 is **already flagged by shipped features** (is_impure=1,
  `max_mult`=426 repeat gate) → not a truly-missed residual.
- **MAGE-549 — hypothesis FALSE.** `best_2cut_modularity` ≈ 0, `louvain_mod` ≈ 0, edge_density
  = 1.0: a **dense uniform complete blob**, a direct counterexample to graph-topology
  bimodality. It is caught **only** by `af_bimod` = 0.908, **above the real-dense max (0.866,
  p90 0.716)** — i.e. the over-merge shows up as **identity-distribution bimodality, not graph
  topology.** But `af_bimod`'s held-out AUC is 0.42 (real paralog families are often bimodal
  too) → **flags MAGE as an outlier, does NOT generalize as a gate.** The spectral λ-magnitudes
  that also flag MAGE (`spec_l3` etc., pct .99) are a **density×weight×size proxy** (corr 0.87
  among real-dense), not sub-structure.
- **GSTM2 triangle-block 13 — no sub-structure exists.** 3-node dense clique (edge_density=1);
  `best2cut` inapplicable, `af_bimod` na → **provably graph-invisible / DNA-bound.**
- **b110 (LOC134758618) — partial.** is_impure=0 missed FP; caught on raw concentration by
  `ew_bimod` (pct .85) and `louvain_mod` (.98), missed by best-2-cut and af_bimod.

**No single RNA feature catches all four as a top-decile / deployable gate.** The weak
"all-available-targets-above-median" bar (a 50%-FP operating point that silently drops the na
b13/b110 cells) is cleared by `ew_cv, ew_range, ew_bimod, af_mean, af_bimod` — but this is not
a usable threshold and none generalize.

---

## Is the residual now graph-invisible / DNA-bound? — DEFINITIVE

- **GSTM2 hub-block 9:** graph-visible (bimodal, best-2-cut top 1%) but **already caught** by
  shipped multiplicity+impurity → not the missed residual.
- **MAGE-549:** graph-topology-**invisible** (dense uniform); visible only as **identity
  bimodality (`af_bimod`)**, which does not generalize → **effectively DNA/identity-bound.**
- **GSTM2 triangle-block 13:** **provably graph-invisible → DNA-bound** (a 3-node clique has no
  detectable sub-structure).
- **Overall:** the residual is **NOT concentration-DNA-bound** (a NEW RNA feature `ew_bimod`
  ties DNA at 4/6@50) but **IS generalization-DNA-bound** — the only feature that both
  concentrates 4/6 and generalizes held-out is the DNA `id_bimod` (0.76). The best-2-cut /
  spectral hypothesis is **refuted as a general rule** (MAGE is a dense-uniform counterexample)
  and **partially vindicated on one block** (GSTM2-9).

---

## Did any usable NEW differentiator surface? — NO

**Winners remain the two already shipped:**
1. **`aln_frac`** (homology extent) — edge workhorse, 0.830 held-out; nothing deployable beats it.
2. **VG shared-multiplicity** (`min_shared_mult` / repeat-hub gate) — the one orthogonal axis
   (ρ=−0.42), plus VG node-multiplicity as the block repeat-hub gate.
   *(+ `core_recip` as a small +0.016 partially-orthogonal edge adder; `aln_frac + core_recip`
   is the one interpretable within-fold improvement.)*

**NEW features evaluated, verdict:**
- Spectral family (`spec_ratio32`/`ratio23`/`nl2`/`ngap23`/`l2`/`l3`/`best2cut_mod`) — genuinely
  size-independent (survive quad) but **confounded by real paralog sub-clades; do not beat the
  base-graph combo (0.723) and hurt it when added.** Interesting, not deployable.
- Edge-weight & identity bimodality (`ew_bimod`, `af_bimod`) — flag specific FPs (b9/b110 and
  MAGE respectively) but **fail size-control / do not generalize** (held-out 0.63 / 0.42).
- Interaction transforms (`contain_x_alnfrac`, `core_x_alnfrac`), `n_shared_nodes` — redundant
  with `aln_frac` (ρ 0.85–0.98).
- Read-degree / multimapping / k-mer / soft-mask — **re-confirmed weak-to-dead**; `mm*` outright
  anti-generalize.

---

## Audit fixes applied to `bench/tp_fp_panel.py` this pass

1. **Precision@K now computed for every ranked feature** (was top-8 only), surfacing the
   `ew_bimod` 4/6@50 tie with DNA `id_bimod`. Concentration claim corrected: residual is
   **not concentration-DNA-bound** but **is generalization-DNA-bound**.
2. **Block 110 (LOC134758618)** added to the per-feature residual test → all 6 confirmed FPs
   (not 3) covered.
3. Misleading flag `catches_all3_above_median` renamed **`catches_all_avail_above_median`**
   with `n_targets_available`, and the prose qualified: na targets are silently dropped; this
   is a weak 50%-FP bar, **not a deployable gate**.
4. Univariate precision@K explicitly labeled **in-sample-oriented**; only `*_combo_oof` are
   true OOF.
5. Block-9 "GSTM2" naming annotated as inherited (dominant gene `LOC115930164`), FP status real.

Determinism re-verified: edge TSV, block TSV, JSON, and stdout all **byte-identical** across
reruns (`PYTHONHASHSEED=0`, seed=0).
