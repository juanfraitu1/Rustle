# Precision–Recall Frontier — RNA multi-copy family definition

**Goal.** A HIGH-PRECISION family definition that avoids the residual over-merge false
positives (the GSTM2 domain-hub, the *fam17* repeat hubs, the collapsed-array oversize blobs)
even at a recall cost, with the whole P/R frontier exposed so the operating point can be picked.
**DNA (the phased-diploid CN oracle `diploid_cn_oracle.tsv`, 57 multi-copy genes) is
VALIDATION ONLY — it never gates an RNA edge.**

Builder / scorer: `bench/precision_recall_frontier.py` → `precision_recall_frontier.tsv` + `.json`.
Reuses the shipped machinery (`family_rna_refine`, `genome_family_def.refine_families`,
`graph_def_refine_sweep.eval_partition`, `rna_only_edge_oracle.oracle_residuals`,
`family_level_pr_current`) so the numbers line up with the committed catalogs. Deterministic
(`PYTHONHASHSEED=0`, `seed=0`, spectral bisection via `numpy.linalg.eigh`); reruns are
byte-identical; the default point reproduces the committed 606-family catalog exactly.

---

## TL;DR — recommended high-precision point

**`gamma = 0.40`, default edge rule (core≥0.19, aln≥0.24, repeat-gate ON), NO best-2-cut.**

It is the frontier point with the **fewest distinct over-merge FP blocks (4, down from 6)**, the
best comparable precision (**fixed-denominator P = 0.917** vs 0.875 default), and it holds recall
at **48/57 recovered** with **zero on-oracle genes lost**. It removes the two avoidable
collapsed-array OVERSIZE blobs (**MPHOSPH8, LOC134758618**) and both *fam17* repeat-hubs. It
**dominates** every "combined" (γ+best-2-cut) point on the real metrics. Its honest costs: a mild
over-split increase (undersize 33→37; 15 divergent-paralog pairs cut) and an **off-oracle**
KRAB-ZNF cost (γ≥0.27 knife-edge). The MAGE-class dense-uniform blobs survive — **no RNA point
removes them; that is the DNA-only floor.**

> Two honest caveats are baked into every number below and stated up front, because they change
> how the headline is read:
>
> 1. **The precision denominator MOVES.** `P_dedup = 1 − distinct_FP / oracle_mapped`, and
>    `oracle_mapped` is *not* fixed: it rises with splitting (48 default → 50 γ0.40 → 51
>    best-2-cut/combined) because splitting mints more oracle-mapped families. Part of any P gain
>    across points is denominator inflation, not FP removal. We therefore also report
>    **`P_fixed48` = 1 − distinct_FP / 48** on the fixed default denominator; that is the
>    comparable number.
> 2. **Recall `R` is a "recovered-SOMEWHERE" metric, blind to over-splitting.** A gene counts as
>    recovered if it survives in *any* multi-copy family, so splitting a real family into two
>    fragments, or shaving copies off, does **not** drop `R`. We therefore report two over-split
>    surrogates that *do* move: **undersize** (DNA diploid-CN ≫ surviving RNA loci) and
>    **tbCut** (divergent-paralog / TRUTHBAR co-membership pairs cut vs default). Read them next
>    to `R`.

---

## Metric definitions

- **P_dedup (moving)** = 1 − distinct_FP_blocks / oracle_mapped  — the 0.85→0.92 headline; **denominator moves**.
- **P_fixed48** = 1 − distinct_FP_blocks / 48 — same numerator, **FIXED** default denominator; comparable across points.
- **P_task** = 1 − flag_instances / oracle_mapped — the task-formula variant (also `P_task_fixed48` in the TSV).
- **R** = oracle multi_copy genes recovered-somewhere / 57 — **blind to over-splitting**.
- **P_Ep** = 1 − E_p-impure_blocks / n_families — protein-purity (also rises with n_families, so partly denom-inflated).
- **undersize** = # multi-copy families with DNA diploid-CN > 1.5×(surviving RNA loci) — the **over-split** tail.
- **tbCut** = # divergent-paralog (TRUTHBAR reciprocal-homology) pairs no longer co-membered vs default — **over-split** cost.

---

## FRONTIER TABLE

| # | point (operator) | nFam | orac | P_dedup(mov) | **P_fixed48** | R | **dFP** | GSTM2 | FOX | MAGE | undsz | tbCut | on-oracle lost |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | **DEFAULT** core≥.19,aln≥.24,γ0.20 | 606 | 48 | .875 | **.875** | **.842** | 6 | 2 | 1 | survives | 33 | 0 | — |
| 2 | F1-alnfrac core≥.26,**aln≥.72** | 418 | 42 | .881 | .896 | **.772** | 5 | **2 (NOT removed)** | 1 | survives | 30 | **518** | **4 genes** |
| 3 | aggressive **γ0.30** | 614 | 49 | .898 | .896 | .842 | 5 | 2 | 1 | survives | 35 | 7 | 0 |
| 3 | aggressive **γ0.40** ★ | 623 | 50 | .920 | **.917** | .842 | **4** | 2 | 1 | survives | 37 | 15 | 0 |
| 4 | best-2-cut **t0.15** | 634 | 51 | .902 | .896 | .842 | 5 | 1 (relabel) | 1 | survives | 39 | 37 | 0 |
| 4 | best-2-cut t0.20 | 626 | 50 | .880 | **.875** | .842 | 6 | 1 (relabel) | 1 | survives | 38 | 37 | 0 |
| 4 | best-2-cut t0.25 | 618 | 50 | .880 | **.875** | .842 | 6 | 1 (relabel) | 1 | survives | 38 | 29 | 0 |
| 4 | best-2-cut t0.30 | 614 | 50 | .880 | **.875** | .842 | 6 | 1 (relabel) | 1 | survives | 36 | 23 | 0 |
| 4 | best-2-cut t0.35 | 608 | 48 | .875 | .875 | .842 | 6 | 2 | 1 | survives | 33 | 7 | 0 |
| 5 | combined γ0.30+b2c t0.25 | 623 | 51 | .902 | .896 | .842 | 5 | 1 (relabel) | 1 | survives | 39 | 29 | 0 |
| 5 | combined γ0.30+b2c t0.30 | 621 | 51 | .902 | .896 | .842 | 5 | 1 (relabel) | 1 | survives | 38 | 23 | 0 |
| 5 | combined γ0.40+b2c t0.15 (DOMINATED) | 641 | 51 | .9216 | **.917** | .842 | 4 | 1 (relabel) | 1 | survives | 39 | 37 | 0 |
| 5 | combined γ0.40+b2c t0.25 (DOMINATED) | 629 | 51 | .9216 | **.917** | .842 | 4 | 1 (relabel) | 1 | survives | 39 | 37 | 0 |

`dFP` = distinct over-merge FP blocks (the real precision numerator). `GSTM2` = # multifam FP
blocks carrying the GSTM2 *name*. `FOX` = # multifam blocks carrying FOXO1 (never 0 anywhere).
`orac` = the moving denominator.

**Reading the table honestly:** on the FIXED denominator, only `γ0.40` reaches `dFP=4` /
`P_fixed48=0.917`. The two γ0.40+best-2-cut "combined" rows tie it at `0.917`/`dFP=4` —
best-2-cut adds **nothing** on top of γ0.40 (see below). best-2-cut at `t=0.20..0.35` gives
`P_fixed48=0.875` — **identical to the default** — i.e. zero precision gain, only over-split
cost.

---

## Which FPs each point removes — and the best-2-cut Q dichotomy

The default catalog's 6 distinct FP blocks are:
`GSTM2+LOC115930164+LOC115930576` (domain-hub), `GSTM2+LOC101129940` (triangle),
`FOXO1+LOC115933254`, `LOC129529978+LOC129529986` (X-array), `MPHOSPH8` (oversize),
`LOC134758618` (oversize).

Best-2-cut spectral (Fiedler) modularity **Q** of each named block explains *what is splittable*:

| block | n | Q | verdict |
|---|---|---|---|
| **GSTM2 domain-hub** GSTM2+LOC115930164+LOC115930576 | 15 | **0.335** | genuine 2-cluster → splittable, **but see below** |
| **FOXO1** FOXO1+LOC115933254 | 15 | 0.319 | fragments but genes stay co-membered → **never removed** |
| LOC134758618 (oversize) | 10 | 0.167 | splits only at t≤0.167 (γ0.40 also removes it) |
| **GSTM2 triangle** GSTM2+LOC101129940 | 3 | — (n<4) | **DNA-only floor** |
| **X-array** LOC129529978+LOC129529986 | 11 | **−0.001** | dense-uniform, edge_density≈1 → **DNA-only floor** |
| **MAGEA9** | 2 | — (n<4) | recall-side floor; **not a counted FP** |

**Precision attribution — where the `dFP 6→4` actually comes from:**

- **`γ0.40` removes {MPHOSPH8, LOC134758618}** — the two collapsed-array OVERSIZE blobs. This is
  the entire real precision gain (`dFP 6→4`).
- **best-2-cut does NOT reduce dFP.** It splits the GSTM2 domain-hub label
  (`gstm2_multifam_count 2→1`, Q=0.335) but the hub **remnant `LOC115930164+LOC115930576`
  stays a 2-gene multifam FP** — the block is *relabelled*, not removed. Net dFP change of adding
  best-2-cut on top of γ0.40 = **0** (4→4). It buys a cosmetic GSTM2-name split, not an FP.
- **FOXO1 is never removed** at any of the 13 points (`foxo1_count=1` everywhere). Its Q=0.319
  fragments the block (dl10→dl2) but keeps FOXO1 and LOC115933254 co-membered.

**Corrected residual-FP floor at the recommended γ0.40 point** (4 surviving distinct FP blocks):
1. `GSTM2+LOC115930164+LOC115930576` — the domain-hub (γ0.40 keeps it whole; best-2-cut would only
   relabel it, leaving the remnant as an FP);
2. `GSTM2+LOC101129940` — the GSTM2 triangle (n=3, Q undefined, **DNA-only**);
3. `LOC129529978+LOC129529986` — the X-array (dense-uniform, Q=−0.001, **the actual MAGE-class
   floor FP present here**);
4. `FOXO1+LOC115933254` — never removed at any RNA point.

(The earlier draft's floor list "GSTM2 triangle, X-array, MAGEA9" was wrong: it omitted the
GSTM2-hub remnant and FOXO1+LOC115933254, and listed MAGEA9, which is a 2-node block and **not a
counted FP** at any operating point.)

---

## Why aln_frac-tightening is the WRONG lever (F1-alnfrac proof)

Point **#2** tightens homology to `core≥0.26 AND aln_frac≥0.72`. The data confirm it is the
anti-lever:

- **It does NOT remove GSTM2** (`gstm2_multifam_count = 2`, unchanged) — GSTM2 and its paralogs
  share extensive sequence and pass any aln_frac.
- **It does NOT remove MAGE** (`mage_recovered = True`, X-array still present) — same reason.
- **It costs 4 real oracle genes** outright: **LOC101130854, LOC101130894, LOC109023386,
  LOC109029264** (R 0.842 → 0.772).
- **It shatters 518 divergent-paralog (TRUTHBAR) co-membership pairs** (`tbCut = 518`, vs ≤37 for
  any splitting point) and shrinks the catalog 606→418.
- Its lone apparent `dFP 6→5` (loses LOC134758618) is an edge-cut *dissolution*, not a clean
  split, on a collapsed/degraded oracle_mapped of 42.

**Conclusion:** tightening homology only kills divergent paralogs = pure recall loss; it never
touches the 2-cluster/dense-uniform over-merges. The precision lever is **splitting (γ / best-2-cut),
not a threshold** — exactly as claimed.

---

## The honest MAGE cardinality floor (no RNA point removes it)

`mage_recovered = True` and `mage_oversize_count = 0` at **all 13 operating points**. The MAGE-A
array is a **dense-uniform blob**: the representative **MAGEA9** is a 2-node block (Q undefined,
never even flagged — it is the *recall-side* floor), and the co-located **X-array
LOC129529978+LOC129529986** has best-2-cut **Q = −0.001, edge_density ≈ 1** and is the
dense-uniform block that **actually survives as a counted over-merge** at every point (it is the
*precision-side* floor). No threshold, no γ, no best-2-cut can bisect it — the copies are RNA-
identical, so only **DNA copy-number** resolves the MAGE-A cardinality (DNA CN 2 for ~18 real
loci). Every RNA operating point above LEAVES it standing; that is the irreducible floor.

---

## Recall cost per point (on-oracle + the over-split axis R is blind to)

- **Homology tightening (F1-alnfrac):** on-oracle R 0.842→0.772, **4 named genes lost**
  (LOC101130854/…894/LOC109023386/LOC109029264), **518 paralog pairs cut**. Zero FP benefit. The
  anti-lever.
- **Splitting operators (γ, best-2-cut, combined):** **zero on-oracle recovered-somewhere loss**
  (R = 0.842, 48/57, at every point) — but this metric is blind to over-splitting, so the honest
  cost is the over-split surrogates:
  - **γ0.40 (recommended):** undersize 33→37 (+4), **15 TRUTHBAR pairs cut**. Newly over-split
    high-CN gene-sets: GSTM2, FOXO1+LOC115933254. Plus the **off-oracle** KRAB-ZNF cost (below).
  - **best-2-cut t0.15:** undersize 33→**39** (+6), **37 TRUTHBAR pairs cut** — the
    *most* over-splitting of the on-frontier points. **The claim "best-2-cut-only = zero recall
    cost of any kind" is retracted** (it has undersize 39 and 37 cut pairs, worse than default's
    33/0). Newly over-split gene-sets include SNRNP25, LOC109023386, LOC115933254,
    FOXO1+LOC115933254, LOC115930164+LOC115930576.
  - **combined γ0.40+best-2-cut t0.15:** undersize 39, **37 TRUTHBAR pairs cut** — strictly worse
    over-split than γ0.40-alone (37 vs 15) for **no distinct-FP gain**.
- **Off-oracle γ caveat:** γ≥0.27 trips the ZNF716/KRAB-ZNF knife-edge (family density 0.261;
  `genome_family_def` GAMMA note) → γ0.30/0.40 split divergent KRAB-ZNF paralog families the
  sparse high-CN diploid oracle cannot see. This is an off-oracle recall cost **not** captured by
  `R`. γ0.20 (default / best-2-cut-only points) preserves them.

---

## Dominance: the "combined max-precision" point is dominated

| metric | γ0.40-alone | combined γ0.40+best-2-cut t0.15 |
|---|---|---|
| distinct_fp | **4** | 4 (tie) |
| P_fixed48 | **0.917** | 0.917 (tie) |
| P_dedup (moving) | 0.920 | 0.9216 |
| oracle_mapped (denom) | 50 | 51 |
| undersize | **37** | 39 (worse) |
| kept_TRUTHBAR | **773** | 751 (22 more pairs cut) |

The combined point's only edge — P_dedup 0.920 → 0.9216 — is **pure denominator inflation**
(50→51 oracle-mapped, `p_dedup_edge_is_denominator_only = true`); on the fixed-48 denominator both
are **0.917**. best-2-cut adds **net-zero distinct FP** (`best2cut_adds_zero_distinct_fp = true`)
and **over-splits strictly more** (+2 undersize, +22 cut pairs). **γ0.40-alone dominates.**

---

## Recommendation

**Primary high-precision operating point: `γ = 0.40` (default edge rule, no best-2-cut).**
- `dFP 6→4` — removes the avoidable collapsed-array OVERSIZE blobs (MPHOSPH8, LOC134758618) and
  both *fam17* repeat-hubs.
- `P_fixed48 = 0.917` (vs 0.875 default) — the best comparable precision on the frontier; moving
  `P_dedup = 0.920`.
- `R = 0.842` (48/57), **zero on-oracle genes lost**.
- Costs (stated, not spun): undersize 33→37, 15 divergent-paralog pairs cut, and the off-oracle
  KRAB-ZNF split (γ≥0.27). Leaves the DNA-only floor standing: GSTM2 domain-hub + triangle,
  X-array (MAGE-class), FOXO1+LOC115933254.

**If the off-oracle KRAB-ZNF cost must be avoided *and* separating the GSTM2 *name* from the hub
is an explicit goal:** use **best-2-cut-only at t=0.15 with γ held at 0.20**. It keeps the
divergent KRAB-ZNF paralogs (γ0.20) and splits the GSTM2 name off the hub — but be explicit that
this is **cosmetic for precision** (`dFP 6→5` only, the hub remnant LOC115930164+LOC115930576
remains an FP), and it carries the *largest* on-oracle over-split cost of the frontier (undersize
39, 37 pairs cut). It is a name-separation choice, not a precision choice.

**What no RNA point can do:** remove the GSTM2 triangle (n=3), the X-array, or the MAGE-A
cardinality. Those are the DNA-only floor — honest and irreducible.

---

## Wiring: a `--high-precision` flag on `family_rna_refine.py`

The recommended point is a **one-parameter change** to the shipped default and is worth wiring as
an opt-in flag (default OFF, so the shipped catalog stays byte-identical):

- `family_rna_refine.build_catalog(gamma=GAMMA)` — thread `gamma` through instead of using the
  module constant, and pass it to the single existing
  `G.refine_families(comps, …, gamma, SEED)` call (line ~266). The edge rule
  (core≥0.19, aln≥0.24, repeat-gate) and the allele-demote stage are unchanged.
- `--high-precision` sets `gamma = 0.40` (from the shipped 0.20). Relax the `assert abs(GAMMA-0.20)<1e-9`
  guard to apply only in the default path so the constant-drift check still protects the default.
- Emit the disclosure in the run header: (i) the off-oracle KRAB-ZNF caveat (γ≥0.27), and (ii) that
  precision is reported on both the moving and fixed-48 denominators. `genome_family_def.py`
  already exposes `--gamma`, so the underlying refiner needs no change.

Best-2-cut is **not** worth wiring as a precision flag (net-zero distinct-FP; cosmetic GSTM2
relabel at real over-split cost). If ever wired, it should be labelled a *name-separation*
option, not a precision option.

---

## Reproduce

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/precision_recall_frontier.py
```

Deterministic; reruns are byte-identical (`determinism_combined_point_byte_identical = true`).
Outputs: `bench/precision_recall_frontier.tsv`, `bench/precision_recall_frontier.json`
(the JSON carries `disclosure_moving_precision_denominator`, `recall_metric_caveat`,
`precision_attribution`, `dominance_analysis`, `residual_fp_at_recommended`,
`best2cut_Q_of_named_FP_blocks`, `gamma_offoracle_caveat`, and per-point `P_dedup_fixed48`,
`undersize_delta_vs_default`, `truthbar_pairs_cut_vs_default`,
`truthbar_pairs_broken_named_vs_default`, `undersize_genesets_new_vs_default`).
