# Copy assignment (O2): definition, decision rule, abstention, certificate

Status 2026-08-10. **This document describes what the shipped `copy_assign` binary computes.** Where an
existing document describes something else, that other document is named and corrected here (§8). Every
number below is measured; the script that produced it is named inline. Anything not measured is marked
**OPEN**.

---

## ⚠⚠ NOTICE — READ BEFORE QUOTING O2 ANYWHERE (2026-08-10)

Until today O2 had **no specification**, and the three places it *was* defined
(`bench/DEFINITIONS_FORMAL.md` lines 141, 225, 343) all said

> copy-assignment (O2) = **max-weight facility location**, assign-or-abstain, no 1/k

`grep -rn facility --include=*.rs src/` returns **ZERO hits**. Facility location is a **theory object**
(`bench/THEORY.md` §5c, Theorems 5–7) machine-checked only on synthetically enumerated instances by
`bench/mwca_approximation_check.py` — a pure-Python script (`numpy`/`scipy.optimize.linprog`, no
subprocess, no BAM) that never touches the binary. **No shipped code path solves it.**

The real object is stated in §1. The precise relationship is stated in §8, and it is not "the doc is
wrong about a detail": **facility location and the shipped rule differ in whether the objective couples
reads at all.** Facility location is only equal to the shipped rule in the degenerate case where every
facility is open and opening is free — which is exactly the case the shipped code is in, because the
copy set is fixed by O1 before assignment begins.

⚠ **Do not quote "facility location" as a description of any measured Rustle result.** It describes an
approximation-guarantee programme that has not been implemented.

---

## 0. What the definition must survive

Same three forms of arbitrariness the O1 spec is written against, re-aimed at O2:

1. **Rule-dependence** — if the answer changes with which of several coexisting gates ran, the
   assignment is an artifact of a build flag. (There are **two** gates in the binary and **two more** in
   unwired modules; §5 and §10 name all four and say which one shipped.)
2. **Threshold-dependence** — if the answer changes with α, the assignment is an artifact of a knob.
   (§1's table gives α a measured plateau.)
3. **Circularity** — if the evaluation consumes what produced the answer, nothing was shown.
   (§9.1: the headline metric **fails this**, and is therefore not quoted anywhere in this document.)

---

## 1. Definition

Fix a family `F` = an ordered set of **copies** `c_1 … c_n` (§2), and a set of **reads** `R` overlapping
them (§2). Each copy carries a feature vector; each read carries an observation vector over the same
features (§3).

> **Objective.** For an assignment `j : R → {1…n}`,
>
> `L(j) = Σ_{r ∈ R} log P(obs(r) | copy_{j(r)})`
>
> maximized over `j`. The emission model `P(obs | copy)` is §3's PSV term + junction term.
>
> **The objective has no term coupling two reads and no term coupling two copies.** There is no
> copy-opening cost, no cardinality constraint on the open set, and no abundance prior in the hard
> path. `L` therefore **decomposes**, and its maximizer is obtained read-independently:
>
> **Decision rule (per read).** `best(r) = argmax_i logl_i(r)`, ties broken to the **lowest copy index**
> (`copy_assign.rs:363-368`, `copy_assign.rs:412-417`). This is the *whole* optimization. Everything
> else in the pipeline is a gate on whether that argmax is allowed to be reported.
>
> **Abstention rule (per read), the shipped default gate.** Let `n` be the number of copies and
> `thr = α / max(n−1, 1)` — Bonferroni over the `n−1` competitors, so α bounds the **family-wide**
> misassignment rate rather than the per-comparison rate (`copy_assign.rs:449-462`). Then
>
> | outcome | condition |
> |---|---|
> | **`tied`** | `min_p(r) ≥ thr` — the read is **unresolvable in principle**: even if it supported `best` at *every* position that distinguishes `best` from its hardest competitor, the certificate could not clear `thr` |
> | **`assigned`** | `min_p(r) < thr` **and** `p(r) < thr` **and** `margin(r) > 0` |
> | **`ambiguous`** | `min_p(r) < thr` but the read's actual evidence does not clear it (or the argmax is not strict) |
>
> where `margin(r) = logl_best − logl_secondbest` (`+∞` for `n = 1`), `p(r)` is the **worst** (largest)
> pairwise significance over the `n−1` competitors and `min_p(r)` the **worst attainable** one; both are
> exact Poisson-binomial upper tails (§4). `margin > 0` is a **strict-MLE** clause: an exact log-LR tie
> is never assigned even when its certificate clears α.
>
> **`tied` vs `ambiguous` is the load-bearing distinction.** `tied` is a statement about the *substrate*
> (the K = 0 floor — this read could not be resolved by any amount of evidence at these copies);
> `ambiguous` is a statement about *this molecule* (it could have been resolved, it wasn't). Both
> abstain. **Neither ever splits weight 1/k** — a read contributes to exactly one copy's hard count or
> to none.

**Emission ≠ decision.** `<out>.assignments.tsv` prints `assigned_copy = best_copy` on **every** row
including `tied` and `ambiguous` ones (`bin/copy_assign.rs:1297`). The argmax is always reported; the
`status` column is what says whether it was believed. Consumers that read `assigned_copy` without
filtering on `status` are reading an unabstained argmax.

### Shipped values, with the measured role of each

| symbol | shipped | where it lives | measured role |
|---|---|---|---|
| **α** | **1e-3** | `copy_assign.rs:213`; CLI `--alpha` | **OPERATIVE — the only constant in the gate.** Plateau **[1e-3, 1e-1]**: over those two decades the assigned count moves **0.384%** (1,052,310 → 1,048,266 of 1,648,460 rows). The slope then steepens **5×** at 1e-4 (−1.81%) and **14×** by 1e-5 (−5.30%). **Shipped value is at the plateau's lower (conservative) edge**, exactly like `c` in the O1 spec. |
| `margin > 0` | strict MLE | `copy_assign.rs:456` | **OPERATIVE, small.** It is the *sole* decider (certificate already cleared, so it alone separates `assigned` from `ambiguous`) on **7,709 / 1,648,460 rows = 0.468%**. ⚠ These are **not** evidence-poor reads — the largest carries **`n_decisive = 2,266`** and still ties two copies exactly. A copy pair that is identical over everything a read spans produces this regardless of read length, which is the K = 0 floor showing up inside the *resolvable* set. |
| **τ (`margin`)** | 2.0 in the library default (`copy_assign.rs:212`), **6.9** at the CLI (`bin/copy_assign.rs:186`) | `AssignParams::margin` | **INERT by default — it is not consulted.** `use_margin_gate` is `false` (`copy_assign.rs:215`), and the gate branches on it (`copy_assign.rs:437`). ⚠ The two defaults **disagree**, so the number a reader finds depends on which file they open. Reachable only via `--margin-gate`. |
| `error_rate` | 0.003 | `copy_assign.rs:209`; CLI `--error-rate` | **CONDITIONALLY LIVE**: used only for PSV columns where the read carries **no per-base QV** (`copy_assign.rs:334-337`, `copy_assign.rs:426`). On QV-bearing HiFi BAMs the per-base `10^(−q/10)` wins everywhere and this constant never enters. ~~**OPEN**: the fraction of columns falling back to it has not been measured.~~ **MEASURED 2026-08-11 (§11)**: across the 5 gorilla control-panel families, **0 of 25,214 alignment records lack QUAL** (primary, secondary *and* supplementary all carry it — minimap2 `-Y`), so on this substrate the fallback **never fires** and `error_rate` is dead. ⚠ It stays live for QUAL-stripped BAMs. ⚠⚠ **What actually binds instead is the clamp inside `phred_err` (`copy_split.rs:189`): `10^(−q/10)` is floored at 1e-4, so every base is scored as at most Q40 — 74% of this substrate's bases are Q93. See §11.2; this constant, not `error_rate`, is what sets the certificate's floor.** |
| `junction_weight` | 5.0 | `copy_assign.rs:210` | **OPERATIVE where junctions distinguish.** A copy-specific junction is worth ±5.0 nats ≈ 7 PSV matches at ε = 0.003, so **junction evidence dominates PSV evidence** wherever both exist. It is a **fixed log-odds, not a calibrated one** — there is no measured plateau. ⚠ This is an absolute-threshold-shaped constant of the kind that has failed 8 times in this project. |
| `boundary_tol` | 4 bp | `copy_assign.rs:211` | Splice-jitter tolerance for junction matching (`boundary_present`, `copy_assign.rs:236`). **OPEN** — never swept. |
| `junction_err` | 1e-4 | `copy_assign.rs:214` | ε for a distinguishing junction in the certificate. Makes ~1 junction worth as much certificate as ~3 PSVs. **OPEN** — never swept. |
| `edit_rate` | 0.2 | `copy_assign.rs:217`; CLI `--edit-rate` | RNA-editing filter: inflates εⱼ on an A↔G column flagged heterogeneous, so an edited base cannot fake copy support. Affects the **certificate only, never the ranking** (`copy_assign.rs`, `assign_read_editing` doc). |
| `rna_editing_filter` | **true** | `copy_assign.rs:216` | DEFAULT ON. Disable with `--no-editing-filter`. |
| `PSV_MIN_ALLELE_READS` | 2 | `copy_assign_pipeline.rs:640` | Read-support PSV validation: a candidate column survives unless the reads positively contradict it. |
| `PSV_MIN_JUDGE_COV` | 4 | `copy_assign_pipeline.rs:644` | A column with `< 4` covering reads is **kept unjudged** — "we cannot invalidate what we cannot see". |
| `iterative_prune` | **false** | `copy_assign.rs:218`; CLI `--iterative-prune` | OFF. Refused outright under `--families` (it would change the roster O1 fixed). |
| `GATE_MIN_READS` | 3 | `denovo_assemble.rs:27` | Not part of the decision — the **invariance certificate** floor in `<out>.quant.tsv` (§6). |
| `QUANT_ERROR` | 0.01 | `copy_assign_pipeline.rs:954` | Soft-EM abundance only (§6), never the hard assignment. |

**The one number to quote for O2 is α = 1e-3, at the lower edge of a measured [1e-3, 1e-1] plateau.**
Everything else in the table is either inert (τ), conditional (`error_rate`), certificate-only
(`edit_rate`, `junction_err`), or unswept (`junction_weight`, `boundary_tol`) — and the unswept ones
should be reported as unswept, not as calibrated.

---

## 2. Input objects

O2 consumes three things. Only the third is new.

**(a) Reads — RECORDS in, MOLECULES out.** The pool is **every mapped alignment record** overlapping
the swept region, **including secondary and supplementary** ones
(`denovo_assemble.rs::reads_in_region`), not just primaries. That is deliberate and it matters: the
shipped BAMs are aligned with minimap2 `-Y` (see the `@PG` line of `GGO_mm.bam`), which writes the
**full SEQ onto secondary records**, so each of a multimapping molecule's records is an independently
scorable observation. Each record's observation is taken **in the frame of its mapped copy** —
`mapped_copy = best_overlap_copy`, the copy whose genomic span the record overlaps most
(`copy_assign_pipeline.rs:791-805`) — and reverse-complemented when that copy is `−`. The pool is *not*
filtered by MAPQ: MAPQ is recorded (for the invariance certificate) and never gates.

⚠ A record is an input, **never an output unit**. A molecule comes from exactly one copy, so its
records are reduced to **one result per (molecule, family)** before anything is emitted or counted
(`assign_family_detailed_once`, `mol_names`):

* **contradiction ⇒ abstain.** If ≥ 2 of a molecule's records are `assigned` and they name ≥ 2 **distinct**
  copies, the molecule is `ambiguous`. This is the assign-or-abstain contract applied at the molecule,
  and it introduces no constant.
* otherwise the representative is the record that observed the most (max `n_decisive`, then max
  `margin`, then min `p_value`, then lowest record index). When the records agree the decision is
  theirs unanimously; this choice only selects which certificate is printed.

**Why not "keep the strongest record".** On the gorilla control panel, arbitrating a contradiction by
margin coincides with simply believing minimap2's **primary flag** on **323/323** SHARP molecules and
**1,629/2,912 = 55.9%** panel-wide — i.e. on some families it *is* the primary flag, which is the defect
that retired the `uniq_agree` metric. There is no planted-truth evidence either way: on the `-Y` sim5x
ladder (`o2_14/mk_y.py`, 3 arms × 1,000 molecules, secondary records carrying full SEQ) **zero**
molecules contradict, so the simulator cannot adjudicate the arbitration — only the real substrate
produces the case, and there truth is unavailable. Abstention is therefore the conservative choice, and
its cost is measured below rather than assumed away.

**(b) A copy set.** Historically this was **re-derived inside `copy_assign`** by its own family
construction (pass-1 skeletons → assemble gate → membership oracle → co-location → refine → rescue),
which is a *second* family definition running beside O1's. Measured consequence at defaults: on GSTM the
O1 catalog emits 4 copies while `copy_assign` emits **0 families** on the same 6,031 reads.

**(c) NEW — the O1 catalog, via `--families` / `--copies-fa`.** `copy_assign --families <out>.copies.tsv
[--copies-fa <out>.copies.fa]` parses a `gw_family_catalog` catalog back into the copy set
(`vg_family/catalog_input.rs`, columns located **by header name**), and **switches off every family
CONSTRUCTION leg** in `detect_and_assign` (skeletons, assemble gate, collapse, membership oracle, POA
diagnostic, co-location, refine, thin-locus rescue). The catalog's `family_id` becomes the emitted
`family_id` — no `CAFAM{i}` is minted — and `<out>.family_join.tsv` records the O1↔O2 join by catalog
`tid`. **This is the mode in which "O1 and O2 talk about the same object" is true by construction rather
than by coincidence**; it reproduces the O1 roster exactly on 6 of 7 control-panel families and refuses
the 7th (RABL2, 5 contigs) loudly, because `copy_assign` is region-scoped.

⚠ **Under `--families`, α's Bonferroni denominator is the catalog's `n`.** That is the point of refusing
a truncated cross-chrom family: honouring only the in-region copies would tighten the certificate over
the wrong `K` and mislabel reads whose true copy is off-region.

---

## 3. The feature space

Two feature families, one shared likelihood.

**PSV columns.** Candidate columns come from **copy-vs-copy alignment** (`discover_psvs`,
`copy_assign_pipeline.rs:385`; default engine = a banded Gotoh affine DP, exact-poasta fallback when the
band is exceeded). They are then **validated by a read pileup**: a column is dropped only if it has
`≥ PSV_MIN_JUDGE_COV` coverage yet fewer than two alleles reach `PSV_MIN_ALLELE_READS`
(`read_supported_columns`, `copy_assign_pipeline.rs:652`). So a PSV is an *assembly* claim confirmed by
*molecules*, and a column no read can judge is kept, not discarded.

A **bubble** is `decisive` iff the copies carry ≥ 2 distinct non-`None` alleles there
(`BubbleGraph::from_copies`, `copy_assign.rs:51`). Per column the read spans, every copy gets
`log(1−e)` or `log(e/3)` with `e` = the read's own base QV, or `error_rate` if it has none.

**Junctions.** A read's intron boundaries are mapped into the mapped copy's **spliced** coordinate frame
via `gen2off` (`copy_assign_pipeline.rs:781-790`). A copy either has a boundary within `boundary_tol` or
it does not; the copy gets `+junction_weight` or `−junction_weight`. A junction is decisive when some
copies have it and others lack it.

`n_decisive` = decisive PSV columns spanned + decisive junctions carried. **Two assignments are computed
per read** (`ReadResult`, `copy_assign_pipeline.rs:808`): `psv` (PSV columns only) and `combined`
(PSV + junctions). **`combined` is the decision**; `psv` exists so `<out>.families.tsv` can report how
much the junction axis adds (the `junction_only` column, §6).

---

## 4. The certificate

For the pair (`best`, competitor `c`), collect every position that **distinguishes** them and that the
read observes: PSV columns where their alleles differ, and junctions one has and the other lacks. Each
such position `j` carries an error probability `εⱼ` (`e/3` for a PSV, `junction_err` for a junction,
raised to `edit_rate` for an editing-flagged column). Let `k` = how many of them support `best`.

* `p_{bc}` = `P(Σ Bernoulli(εⱼ) ≥ k)` — the exact Poisson-binomial upper tail, an `O(n²)` DP with no
  normal approximation (`poisson_binomial_upper_tail`, `copy_assign.rs:146`).
* `attain_{bc}` = `Π εⱼ` — the value `p_{bc}` would take if the read supported `best` at *every*
  distinguishing position. This is the **identifiability bound**.

Then `p(r) = max_c p_{bc}` (the least-significant competitor governs) and
`min_p(r) = max_c attain_{bc}`. `min_p ≥ thr` ⟹ **no possible molecule of this length and coverage
could resolve this read** ⟹ `tied`. That is the K = 0 floor, computed per read rather than asserted.

⚠ `min_p` is computed **relative to the same `best`** as `p` (`read_copy_evidence`,
`copy_assign.rs:311-314`). It is the bound for *this* argmax, not a copy-set-wide invariant.

---

## 5. Every constant, and where it comes from

| constant | value | file:line | kind |
|---|---|---|---|
| `AssignParams::error_rate` | 0.003 | `src/rustle/vg_family/copy_assign.rs:209` | compile-time default; CLI `--error-rate` |
| `AssignParams::junction_weight` | 5.0 | `copy_assign.rs:210` | compile-time, **no CLI, no env** |
| `AssignParams::boundary_tol` | 4 | `copy_assign.rs:211` | compile-time, **no CLI, no env** |
| `AssignParams::margin` (τ) | 2.0 | `copy_assign.rs:212` | compile-time; **inert** unless `--margin-gate` |
| CLI `--margin` (τ) | 6.9 | `src/bin/copy_assign.rs:186` | CLI; **inert** unless `--margin-gate` |
| `AssignParams::alpha` | 1e-3 | `copy_assign.rs:213` | compile-time default; CLI `--alpha` |
| `AssignParams::junction_err` | 1e-4 | `copy_assign.rs:214` | compile-time, **no CLI, no env** |
| `AssignParams::use_margin_gate` | false | `copy_assign.rs:215` | compile-time; CLI `--margin-gate` |
| `AssignParams::rna_editing_filter` | true | `copy_assign.rs:216` | compile-time; CLI `--no-editing-filter` |
| `AssignParams::edit_rate` | 0.2 | `copy_assign.rs:217` | compile-time default; CLI `--edit-rate` |
| `AssignParams::iterative_prune` | false | `copy_assign.rs:218` | compile-time; CLI `--iterative-prune` |
| `PSV_MIN_ALLELE_READS` | 2 | `copy_assign_pipeline.rs:640` | compile-time, **no CLI**; whole filter off via `RUSTLE_PSV_READFILTER=0` |
| `PSV_MIN_JUDGE_COV` | 4 | `copy_assign_pipeline.rs:644` | compile-time, **no CLI** |
| `QUANT_ERROR` | 0.01 | `copy_assign_pipeline.rs:954` | compile-time, **no CLI** |
| `GATE_MIN_READS` | 3 | `denovo_assemble.rs:27` | compile-time; certificate only |
| `PSV_BAND_MARGIN` | 1024 | `copy_assign_pipeline.rs:405` | compile-time; band slack (exact-poasta fallback preserves correctness) |
| `MOSAIC_EPS` | 0.01 | `copy_assign_pipeline.rs:1407` | compile-time; mosaic detector |
| `MAX_MOSAIC_SITES` | 250 | `copy_assign_pipeline.rs:1408` | compile-time; O(sites²) cap |
| `MosaicParams` (11 fields) | see `mosaic.rs:77-93` | `mosaic.rs:77` | compile-time defaults; 6 of 11 overridable via `RUSTLE_VG_MOSAIC_*` (`mosaic.rs:99-110`) |
| `RefineParams::min_identity` | 0.80 | `denovo_pipeline.rs` `RefineParams::default` | ⚠ **not exposed by `copy_assign`** — it hardcodes `Default::default()`. Moot under `--families` (refine is off). ⚠⚠ **Since X.4 (2026-08-11) this is NOT the floor refine actually applies on the default path:** refine's primary tier now comes from the shared `er_primary_tier()`, which under the shipped `RUSTLE_ER_SENSITIVE_ONLY` default is `-k 11 -w 5` at **`sensitive_identity` 0.60**. `min_identity` 0.80 binds only under `RUSTLE_ER_SENSITIVE_ONLY=0`. Read the run's own `<prefix>.refine.rule.tsv`, not this row. |
| `RefineParams::min_coverage` | 0.50 | `denovo_pipeline.rs:3202` | same |
| `RefineParams::sensitive_identity` | 0.60 | `denovo_pipeline.rs:3211` | env `RUSTLE_SENSITIVE_IDENTITY` |
| `--min-copies` | 2 | `bin/copy_assign.rs:112` | CLI; **not applied** under `--families` |
| `--win` | 5,000,000 | `bin/copy_assign.rs:131` | CLI; **not applied** under `--families` |
| `--read-cap` | 6000 | `bin/copy_assign.rs:123` | CLI; **documented NO-OP** — see §10 |

**Environment variables that reach the O2 path:** `RUSTLE_PSV_READFILTER`, `RUSTLE_PSV_MINIMAP2`,
`RUSTLE_PSV_POASTA`, `RUSTLE_PSV_BAND`, `RUSTLE_POA_CAP`, `RUSTLE_INTRON_PSV`, `RUSTLE_MINIMAP2`,
`RUSTLE_TIMING`, `RUSTLE_POSTERIOR_PRIOR`, `RUSTLE_VG_MOSAIC_*` (6), `RUSTLE_SENSITIVE_IDENTITY`,
`RUSTLE_NO_RECOMBINANT_ABSTAIN` (⚠ **has no effect** — its module is unreachable, §10).

---

## 6. Outputs — what every column means

### `<out>.assignments.tsv` — one row per (read, family)

| column | meaning |
|---|---|
| `read_name` | BAM query name. ⚠ **Not unique across FAMILIES**: a molecule appears once per family it was pooled into (control panel: 570 of MAGEA's 2,415 molecules are in two families). It is unique **within** a family — one row per (molecule, family), never one row per alignment record. |
| `family_id` | `GWFAM{i}` under `--families` (the catalog's own id); `CAFAM{i}` otherwise. |
| `assigned_copy` | **`best_copy` — the argmax, printed unconditionally**, including on `tied`/`ambiguous` rows. Filter on `status`. |
| `status` | `assigned` \| `ambiguous` \| `tied` (§1). |
| `n_decisive` | decisive PSV columns spanned + decisive junctions carried. `0` ⟹ nothing could distinguish. |
| `margin` | `logl_best − logl_second`, **printed at 3 decimals** — see §9.2. |
| `p_value` | worst-competitor certificate `p(r)`, `{:.3e}`. |
| `min_p_value` | identifiability bound `min_p(r)`, `{:.3e}`. Compare against `α/(n−1)`, **not** α. |
| `as_best`, `as_second`, `as_margin` | raw minimap2 `AS` of the read's best and runner-up placements and their difference. **Reported, never decisive** — raw `AS` is length-confounded (`read_conflict.rs:145-148`). |
| `as_per_base_best`, `as_per_base_2nd` | the length-fair versions. Also non-decisive. |

### `<out>.families.tsv` — one row per family

`family_id · chrom · n_copies · rescued_copies · collapsed_copies · n_reads · psv_cols ·
resolvable_psv · resolvable_j · junction_only · assigned_j · uniq_agree · uniq`

* `psv_cols` — surviving PSV columns after read-support validation.
* `resolvable_psv` / `resolvable_j` — reads with `n_decisive ≥ 1` on the **PSV-only** and the
  **PSV+junction** feature set. Their difference is what the junction axis buys.
* `junction_only` — reads resolvable *only* because of junctions.
* `assigned_j` — reads with `status == assigned` under the combined decision.
* `uniq` / `uniq_agree` — the headline metric's denominator and numerator. ⚠ **See §9.1 before using
  them.** `uniq` counts only reads that are **both** `assigned` **and** MAPQ > 0
  (`denovo_pipeline.rs:2218`), so its denominator is conditioned on the prediction.

### `<out>.quant.tsv` — one row per copy

`family_id · copy_index · copy_tid · copy_chrom · copy_start · copy_end · abundance ·
ci95_halfwidth · n_reads_hard · anchored_reads · tie_invariant · junction_invariant`

* `abundance` — **soft** per-copy fraction from an EM over per-read PSV likelihoods
  (`soft_quantify_em`, `copy_assign_pipeline.rs:961`), summing to 1. At zero informative PSVs it returns
  the **uniform prior** — the honest identifiability floor, deliberately not a guess. This is the one
  place a coupled (θ-carrying) objective exists, and **it does not feed the hard assignment**.
* `n_reads_hard` — reads with `status == assigned` and `best_copy == this copy`.
* `anchored_reads` — of those, how many had MAPQ > 0 (`anchored_support`, `bin/copy_assign.rs:921`).
* `tie_invariant` — `anchored_reads ≥ GATE_MIN_READS`. minimap2 sets MAPQ 0 exactly when the primary is
  an arbitrary tie, so a MAPQ > 0 read cannot be relabeled: **this copy exists under every tie-break**.
* `junction_invariant` — the copy is pinned by ≥ 3 reads carrying a copy-specific junction, which
  identifies it by splice structure regardless of the primary label.
* A copy is invariant overall iff `tie_invariant || junction_invariant`. **FALSE means the copy's
  existence leans on the arbitrary primary flag** — that is the column to check before claiming a copy.

### `<out>.famcn_readonly.tsv`

`family_id · chrom · n_copies · n_reads · chi_H · depth_cn · regime · famcn_readonly`.
`chi_H` = the read-conflict-graph lower bound (`chi_h_with_junctions`); `depth_cn` = the read-depth leg,
`NA` unless `--lambda-global`/`--lambda-file`; `famcn_readonly = max(chi_H, depth_cn)`; `regime` =
`reference_collapsed` when any collapsed copy was found, else `reference_resolved`.

### `<out>.family_join.tsv` — **only under `--families`**

`family_id · copy_index · copy_tid · catalog_family_id · catalog_copy_idx · chrom · start · end ·
n_reads_hard`. One row per assigned copy, keyed by catalog **`tid`** (never by position). This is the
O1↔O2 join key. The run **aborts** if any supplied catalog copy fails to appear here.

### Opt-in outputs

`<out>.posterior.tsv` (`--posterior`) — per-read softmax posterior over copies + the **consistent zone**
(copies above a 0.01 floor, `bin/copy_assign.rs:1316`); localizes even unassignable reads. Prior is
uniform unless `RUSTLE_POSTERIOR_PRIOR=abundance`.
`<out>.em.tsv` / `.em_abundance.tsv` (`--em`) — the soft relaxation.
`<out>.mosaic.tsv` — recurrent copy-switch breakpoints. ⚠ **report-only**, see §9.3.
`<out>.psv_reads.tsv` / `.psv_copies.tsv` / `.psv_cols.tsv` (`--dump-psv`) — the raw genotype matrix.
`<out>.phase_blocks.tsv` / `.phased_haplotypes.tsv` / `.phased_reads.tsv` / `.phase.gfa` (`--phase`) —
⚠ **a relabeling of the assignment, not an independent phasing**; see §9.4.

---

## 7. Verification — does the binary do what §1 says?

**Reproduction experiment** (`/home/juanfra/winloci_scratch/s1_unblock/o2spec/verify_o2_rule.py`).
§1's decision rule was re-implemented **from this document's text alone** (in Python, without consulting
the Rust) and run against the columns of **171 independent historical `copy_assign` runs** in
`/home/juanfra/winloci_scratch` spanning 2026-06-26 → 2026-07-21, **1,648,460 read-rows**, recomputing
`status` from `(p_value, min_p_value, margin, n_copies)` at α = 1e-3.

| | rows |
|---|---|
| total | 1,648,460 |
| **undecidable from the file** (printed `margin == "0.000"`, sign unknown) | 9,139 (0.55%) |
| decidable | 1,639,321 |
| **reproduced exactly** | **1,639,305 = 99.9990%** |
| residual | 16 (2 files, `mval_ID_386_*`) — printed `min_p` rounds to exactly `α/(n−1)` at 3 sig figs; a printing tie at the boundary, not a rule difference |

**Every row this document's rule can decide from the emitted file, it decides correctly.** The rule in
§1 is the shipped rule.

**α plateau.** Same corpus, sweeping α (harness: recompute `thr = α/(n−1)` per row):

| α | assigned | Δ vs 1e-1 |
|---|---|---|
| 1e-1 | 1,052,310 | — |
| 1e-2 | 1,050,351 | −0.186% |
| **1e-3 (shipped)** | **1,048,266** | **−0.384%** |
| 1e-4 | 1,033,243 | −1.81% |
| 1e-5 | 996,547 | −5.30% |
| 1e-9 | 953,739 | −9.37% |
| 1e-50 | 725,470 | −31.06% |

Per-panel the plateau is exact on 3 of 6 large panels (GSTM, DAZ, MEP1AP4 give **identical** counts over
[1e-3, 1e-1]); PCDHB and MAGEA move by < 0.15%. ⚠ **α is not inert** — unlike O1's τ it does bind, it is
just flat where it is set. Below 1e-6 it bites hard and unevenly (DAZ loses 43% of its assignments by
1e-12 while MAGEA loses 0.7%), so **an α change is not a uniform tightening across families.**

**OPEN.** `junction_weight`, `boundary_tol`, `junction_err` have never been swept. Their roles above are
mechanism statements, not measurements, and should be reported that way.

---

## 8. Reconciliation with `bench/DEFINITIONS_FORMAL.md`

**The real object is the per-read decomposed argmax of §1, with a Poisson-binomial abstention
certificate. It is not facility location, and it never was.**

The precise relation, which is why the confusion was easy to make:

> Facility location (`bench/THEORY.md` §5c; `bench/mwca_approximation_check.py`) maximizes
> `f(S) = Σ_r max_{c ∈ S ∩ N(r)} w(r,c)` over open sets `S ⊆ C` with `|S| ≤ K`. All of its content —
> submodularity, the `(1−1/e)` greedy guarantee, the LP bound — lives in the **choice of `S`**.
>
> In the shipped pipeline **`S` is not chosen by O2.** It is the family, fixed by O1 before any read is
> scored, and under `--families` it is literally the O1 catalog file. There is no opening cost and no
> cardinality budget. With `S` fixed and `K = |S|`, `f(S)` degenerates to `Σ_r max_c w(r,c)` — the
> per-read argmax. **The shipped rule is the degenerate instance of facility location in which the
> optimization has already been done by a different objective.**
>
> And it is not *only* that degenerate instance: facility location has **no abstention**. Every read
> takes its best open facility. The shipped rule's `tied`/`ambiguous` outcomes have no counterpart in
> `f`, and they are 36.3% of all reads (§9.1). So the shipped object is neither a special case of, nor
> generalized by, the formulation the document names.

**Action taken.** All three sites in `bench/DEFINITIONS_FORMAL.md` (lines 141, 225, 343) have been
corrected in place to name the shipped rule and to mark facility location as an **unimplemented design
goal**, with a pointer to this document. `bench/THEORY.md` §5c is *not* wrong and is left alone — it is
explicitly a theory section — but it should never be cited as a description of a Rustle **run**.

**"assign-or-abstain, no 1/k" survives verbatim.** That half of the sentence is exactly what ships. Only
the optimization name was wrong.

---

## 9. Known defects (stated, not hidden)

### 9.1 The headline metric is anti-correlated with correctness — do not quote it

`uniq_agree / uniq` ("unique-mapper agreement", `denovo_pipeline.rs:5835`) asks *how often the
assignment agrees with minimap2's primary flag*. Three independent faults:

1. **Its denominator is conditioned on the prediction.** `uniq` increments only when `assigned_j` is
   true (`denovo_pipeline.rs:2218`), so **abstention can never lower it.** This is the
   truth-contains-the-prediction shape that has already killed six metrics in this project.
2. **It measures the wrong thing.** *(Inherited result, not re-run here.)* Rewriting the sim5x BAMs so that *every* primary flag sits on the
   **wrong** copy leaves accuracy against the planted label unchanged (K=2 200/200, K=4 200/200,
   K=8 198/200) while this metric goes **100.0% → 0.0%**. A provably correct run scores zero.
3. **It is empty exactly where O2 exists.** It is MAPQ > 0-only, and O2's reason to exist is the MAPQ-0
   regime. Measured (same harness, `verify_o2_rule.py`) over 227 output sets / 486 family rows: `uniq` covers **13.49% of all reads** and
   **21.16% of assigned reads**, and is **undefined (`uniq == 0`) on 6.8% of family rows**. The
   aggregate it reports is 99.28% — a number produced on the one-eighth of the data where the question
   is easiest.

**Replacement is OPEN.** Until one exists, quote **`tie_invariant || junction_invariant`** from
`quant.tsv` (a certificate, not an accuracy).

⚠⚠ **The "63.7% assigned" figure is WITHDRAWN and must not be requoted.** It counted 1,649,606
**alignment-record rows** on a corpus produced before the record→molecule reduction, where a molecule
contributed one row per placement and 73.8% of multi-record molecules were counted as assigned twice, to
different copies. On the gorilla control panel the same statistic recomputed per MOLECULE is
**2,097 / 7,121 = 29.4%** (Wilson [28.4%, 30.5%]), against 6,723 / 7,121 = 94.4% of molecules "assigned"
under the old record-level count. The genome-scale rate has **not** been recomputed — it requires
re-running the corpus.

### 9.2 `assignments.tsv` cannot fully re-derive its own decision

`margin` is printed at `{:.3}` (`bin/copy_assign.rs:1654`) while the rule tests `margin > 0`. **538,907
of 1,648,460 rows (32.7%) print `0.000`**; on most of them the certificate has already decided the
outcome, but on **9,139 (0.55%)** the printed `0.000` is the *only* thing standing between `assigned`
and `ambiguous` and its **sign is unrecoverable from the file**, so an auditor cannot reproduce those
statuses. Likewise `p_value`/`min_p_value` at `{:.3e}` produce boundary ties (§7's residual 16).
**Fix: print `margin` at full precision, or emit the sign.** Not done here — it changes an output format
and needs its own regression.

### 9.3 The abstention arm for recombinant reads is detected but never applied

`mosaic.rs` runs on **every** family (`copy_assign_pipeline.rs:1538-1541`) and detects reads whose
per-site copy pattern **switches** at a contiguous breakpoint — a read that carries copy A's alleles in
one arm and copy B's in another, i.e. a molecule that belongs to **no single copy**. The count reaches
`<out>.families.tsv`/`.mosaic.tsv`. **It never changes `status`.** Such a read is still force-assigned
to whichever copy wins the argmax. The module that would abstain it (`recombinant_abstain.rs`,
documented "DEFAULT-ON") is unreachable — §10. **This is a real gap in the abstention rule, not a
documentation gap.**

### 9.4 `--phase` is a relabeling of the assignment

`<out>.phased_reads.tsv`'s `haplotype` is `best_copy` when `status == assigned` and `−1` otherwise
(`bin/copy_assign.rs:1532-1537`); `n_psv_spanned` is `n_decisive` and `margin` is `log_lr_margin`. It is
a projection of `assignments.tsv`, by construction.

**Measured** (`/home/juanfra/winloci_scratch/s1_unblock/o2spec/verify_phase_is_a_relabeling.py`): over
**10 historical `--phase` run pairs** with both files present (6 non-empty: `cg_gstm2`, `cg_gstm3`,
`v2_gstm`, `vgo2_phase`, `reads_ID_131`, `reads_ID_131_AMYLASE`), **119,524 read-rows**, the multiset of
`(read_name, family, haplotype)` equals the multiset of `(read_name, family, assigned_copy-if-assigned-
else −1)` **exactly — symmetric difference 0 on every file and 0 in total.**

⚠ A per-read-name join gives a spurious 91.3% instead, because **read names repeat within a family**;
the identity is only visible as a multiset. Anyone who has "checked" this pairwise and found
disagreement has hit that, not a real difference.

The `--phase` help text previously claimed a *"min-path-cover over the PSV graph"*. No such computation
exists in the O2 path. **Corrected 2026-08-10** — the flag help now says what it does.

### 9.5 Structural limits carried forward

* `copy_assign` is **region-scoped**: a family spanning several contigs (RABL2, 5) cannot be assigned.
  Under `--families` this is a loud refusal; without it, a silent truncation.
* `copy_assign` **exposes none of the O1 definition parameters** (it hardcodes `RefineParams::default()`
  and γ = 0.20). Under `--families` this is moot because refine is off; **without
  `--families` it means the two objectives can be run at different definitions with no warning.**
  ⭐ **X.4 (2026-08-11) closed the half of this that mattered most — the EDGE RULE.** Refine no longer
  resolves its own seed: `er_primary_tier()` / `er_sensitive_only()` are shared with O1's
  `homology_edges_all_reps_pooled`, so both sides run `-c -X --no-long-join -k 11 -w 5` at 0.60 by
  default (previously O2 ran `-x asm20` at 0.80 on 10 of 13 panel calls). It is **checkable per run**,
  not asserted: both sites emit a data-free `rule.tsv` and
  `diff <p>.rule.tsv <p>.refine.rule.tsv` was **EMPTY on all 5 refine calls** of the gorilla MAGEA
  control region, while the same diff prints 5 lines under `RUSTLE_ER_SENSITIVE_ONLY=0`.
  ⚠ **What is still NOT shared is the CLUSTERING** — O1 uses γ-quasi-clique, refine uses connected
  components — and `--min-identity` / γ still cannot be passed to `copy_assign`.
* `--read-cap` is a documented no-op (§10).

---

## 10. Orphaned modules — verdicts, with evidence

Measured with `grep -rn "<mod>::" src/ tests/ benches/` excluding the module's own file.

| module | lines | `#[test]` | call sites outside itself | verdict |
|---|---|---|---|---|
| `phasing.rs` | 742 | 12 | **0** | **DELETED** (this change) |
| `o2_materialize.rs` | 1,874 | 9 | **0** (two doc-comment mentions in `bin/copy_assign.rs` only) | **DELETE — not executed here** |
| `o2_columns.rs` | 334 | 2 | 1 — `o2_materialize.rs:44` | **DELETE with the cluster** |
| `o2_margin_gate.rs` | 355 | 3 | 1 — `o2_materialize.rs:45` | **DELETE with the cluster** |
| `recombinant_abstain.rs` | 956 | 6 | 1 — `o2_materialize.rs:46` | **KEEP AND CONNECT** |

**Test-count context.** The O2 modules a shipped binary actually executes carry 95 `#[test]`s
(`copy_assign.rs` 31, `copy_assign_pipeline.rs` 38, `em_copy_assign.rs` 11, `bin/copy_assign.rs` 15).
The five modules above carry **32**. So **25.2% of the "O2 tests" guard unreachable code** (17.8% if
`mosaic.rs`'s 22 and `copy_split.rs`'s 31 are also counted as O2). The complaint is real at either
denominator.

### `phasing.rs` → **DELETE. Executed.**

Evidence: (i) zero call sites; (ii) the flag it names, `--vg-phase`, **does not exist** — `RustleConfig`
has a `vg_phase` field (`types.rs:945`) that is set to `false` at construction (`types.rs:1292`) and
**never read anywhere**; (iii) it is the **wrong object**: it is *diploid* MEC over a *binary* allele
matrix with `h_B = ¬h_A`, over one copy's reads, whereas O2 is *k*-copy over 4-letter alleles across a
family. It is not a disconnected piece of O2; it is a different problem.

⚠ **It was load-bearing for a document.** `bench/THEORY.md:348` read *"Definition (MEC, as implemented in
`vg_family/phasing.rs`)"*, and `bench/soto/bam_tie_signals.md:132` cited it for polyA fields. **Deleting
the file without those edits would have created a new documented-≠-shipped gap of exactly the kind this
spec exists to close**, so both were corrected in the same change. That citation is also *why* the file
had to go rather than stay: a theory document asserting "as implemented" about code no binary runs is
the same defect as the facility-location line.

### The `o2_materialize` cluster → **DELETE 3, CONNECT 1. Not executed in this change.**

`o2_materialize` is a byte-parity Rust port of the Python prototype
`bench/o2_vg_visualization.py::materialize_family`. Nothing imports it, and it is the **sole** importer
of the other three.

* `o2_columns` (pileup → per-column alleles) duplicates `fill_psv_obs`
  (`fill_psv_obs`, `copy_assign_pipeline.rs:725`). **Redundant.**
* `o2_margin_gate` duplicates the **legacy τ-margin gate** at `MARGIN = 2.0` (`o2_margin_gate.rs:41`) —
  i.e. a **third** default for τ, next to the library's 2.0 and the CLI's 6.9. **This module is a
  concrete instance of the "which rule is O2?" problem** and is the strongest argument for deleting it.
* `o2_materialize` duplicates catalog loading, which `catalog_input.rs` now does for real (§2c).
* `recombinant_abstain` is **the exception**: it implements an abstention outcome the shipped gate does
  **not** have (§9.3). Its module doc says "DEFAULT-ON, disable with `RUSTLE_NO_RECOMBINANT_ABSTAIN=1`",
  and that sentence is **false in the shipped binary** — the env var has no effect because nothing calls
  `apply_abstain_to_vg`. It is genuinely unwired functionality, not dead code.

**Why the cluster deletion was not executed here:** deleting `o2_materialize` removes
`recombinant_abstain`'s only consumer, and connecting `recombinant_abstain` to the shipped path
**changes the decision rule** (it adds a fourth status). That needs its own RED-before/GREEN-after and a
control-panel re-run — it is not a documentation change and must not ride along inside one. The verdict
and the ordering are recorded here so the next run does not have to re-derive them:

1. Connect `recombinant_abstain` behind an opt-in flag; make `RUSTLE_NO_RECOMBINANT_ABSTAIN` mean
   something; add a fourth `status`.
2. Then delete `o2_materialize`, `o2_columns`, `o2_margin_gate` (2,563 lines, 14 tests) and the
   `--read-cap` flag with them.

### Verification of this change

`cargo test --release --all-targets` (exit code captured to a file, never piped):
**767 passed / 0 failed / 11 ignored**, versus the **779 / 0 / 11** baseline — exactly **−12**, which is
`phasing.rs`'s own 12 tests and nothing else. The 25-region gorilla control panel was re-run with the
rebuilt `gw_family_catalog` and `diff -rq` against the post-M1 reference
(`/home/juanfra/winloci_scratch/o1fix_o2audit/fix/out_m1fix/`) reports **no differences across all 75
files**. Expected by construction — a deleted unreachable module and a doc comment cannot reach a
shipped code path — but measured rather than assumed.

---

## 11. O2 on REAL data — commitment rate + abstention certificate (2026-08-11, O2.3)

**⚠⚠ NOTHING IN THIS SECTION IS AN ACCURACY.** Real data has no planted truth: not one read below has
a known true copy. Every number here answers *how often does the method commit, and with what
certificate* — never *how often is it right*. **A high commit rate is not evidence of correctness and a
low one is not evidence of error.** Anyone quoting `assigned` as `correct` is misreading the table, and
this sentence must travel with the numbers. The retired `uniq_agree` metric (§9.1) is not reported here
in any form; its denominator was conditioned on the prediction, so abstention could never lower it.

**Substrate.** Gorilla, and gorilla only: mGorGor1 (`GGO.fasta`, Kamilah) × `GGO_mm.bam` IsoSeq
(OR6737), the 25-region control panel at `/home/juanfra/winloci_scratch/o1_famctl/bams/`. The copy set
is the post-M1-fix O1 catalog (`/home/juanfra/winloci_scratch/o1fix_o2audit/fix/out_m1fix/`) consumed
through `--families … --copies-fa …`, so O1 and O2 talk about the same object by construction (§2c).
Scripts: `/home/juanfra/winloci_scratch/close_o1o2/{o23_report.py, mech_check.py, cert_check.py}`;
full output `/home/juanfra/winloci_scratch/close_o1o2/O2_REAL_DATA.txt`.

**The denominator, and why it is outside the method's control.** Every distinct read with a **primary**
alignment (`-F 2308`) inside the swept intervals — fixed by the region geometry before a single read was
scored, and never touched by the assignment. It was counted twice, once from the whole-genome
`GGO_mm.bam` and once from the panel slice: **identical on all 5 panels, 0 reads differ**, so the
`-M -L` slice introduces no bias here (the subset-BAM trap does not bite for primaries). **N = 5,378.**

A read gets one of four outcomes. `unevaluated` is not an O2 abstention — it means the read's family
never existed, which is an O1 scope fact; it is kept in the denominator precisely so the abstention
fractions cannot be inflated by dropping it. Reads in more than one family (548/5,378) are reduced by
`assigned > ambiguous > tied`.

⚠ **Which denominator to quote.** The fractions in 11.1 use the **fixed N**, and those are the numbers
to quote. Anywhere a rate below is written `x / evaluated` (11.3), the denominator excludes
`unevaluated` and is therefore conditioned on **O1's** family extent — not on O2's prediction (an
abstained read is still evaluated, so abstention still lowers every rate here, which is the failure
shape that killed seven earlier metrics). Those conditioned rates are supporting detail, not headlines.

### 11.1 Outcomes, stratified on MAPQ (Wilson 95%)

| stratum | n | assigned | ambiguous | tied (K=0) | unevaluated |
|---|---|---|---|---|---|
| **ALL** | 5,378 | **0.3466** [0.3340, 0.3594] | 0.5534 [0.5400, 0.5666] | 0.0084 [0.0063, 0.0112] | 0.0917 [0.0842, 0.0997] |
| **MAPQ == 0** (O2's whole remit) | 497 | **0.0302** [0.0184, 0.0492] | 0.9256 [0.8991, 0.9455] | 0.0262 [0.0153, 0.0442] | 0.0181 [0.0096, 0.0341] |
| **MAPQ > 0** (easy stratum) | 4,881 | **0.3788** [0.3653, 0.3925] | 0.5155 [0.5014, 0.5295] | 0.0066 [0.0046, 0.0092] | 0.0992 [0.0911, 0.1079] |

**The pooled row is dominated by the easy stratum and must not be quoted alone** — MAPQ > 0 is 90.8% of
the denominator, and the commit rate falls **12.5×** (0.3788 → 0.0302) when you look only at the reads
O2 exists for.

Per panel (all reads):

| panel | copies/family | n | assigned | ambiguous | tied | uneval | MAPQ==0 reads |
|---|---|---|---|---|---|---|---|
| APOBEC3 | 2 | 304 | 0.9539 | 0.0296 | 0.0099 | 0.0066 | 0 |
| GSTM | 4 | 1,551 | 0.8594 | 0.0858 | 0.0026 | 0.0522 | 0 |
| HERC2 | 3 | 604 | 0.0248 | 0.5066 | 0.0182 | 0.4503 | 1 |
| MAGEA | 9 and 2 | 2,332 | 0.0304 | 0.9455 | 0.0111 | 0.0129 | 496 |
| SHARP | 2 | 587 | 0.2641 | 0.5503 | 0.0017 | 0.1840 | 0 |

⚠ **A limit of this panel, stated up front: the MAPQ-0 stratum is essentially ONE family.** 496 of the
497 MAPQ-0 reads are MAGEA's; APOBEC3, GSTM and SHARP contribute **zero**. So the MAPQ-0 row is a
MAGEA number with a Wilson interval, not a five-family average, and it should be quoted that way.
⚠ RABL2, the 6th catalog family, is **refused** by `copy_assign` (5 contigs, region-scoped — assigning
it would silently truncate the copy set). The refusal is loud (exit 1), not a partial answer.
⚠ The 19 single-copy/pseudogene control regions carry **0 catalog copies** post-M1-fix, so O2 has no
object there and evaluates 0 reads. That is the correct behaviour and it is why every rate above is
quoted over the 5 multi-copy panels only.

### 11.2 The certificate that every commitment carries

Gate: commit iff `p(r) < thr = α/(n−1)`, α = 1e-3, `n` = the **O1 catalog's** copy count (§2c).

| family | n | thr | assigned rows | max p achieved | median p | violations |
|---|---|---|---|---|---|---|
| APOBEC3/GWFAM0 | 2 | 1.000e-3 | 363 | 1.000e-4 | 0 | 0 |
| GSTM/GWFAM0 | 4 | 3.333e-4 | 1,337 | 3.333e-5 | 1.6e-302 | 0 |
| HERC2/GWFAM0 | 3 | 5.000e-4 | 35 | 3.333e-5 | 0 | 0 |
| MAGEA/GWFAM0 | 9 | 1.250e-4 | 48 | 3.333e-5 | 1.9e-230 | 0 |
| MAGEA/GWFAM1 | 2 | 1.000e-3 | 26 | 1.000e-4 | 0 | 0 |
| SHARP/GWFAM0 | 2 | 1.000e-3 | 288 | 6.667e-6 | 2.9e-43 | 0 |

**2,097 assigned rows, 0 certificate violations, worst achieved `p/thr` = 0.2666.** The medians sit
hundreds of orders of magnitude below the bound.

⭐ **The `max p` column is not noise — it is a QUANTUM, and this run is what exposed it.**
`phred_err` (`copy_split.rs:185-190`) **clamps the per-base error probability to `[1e-4, 0.25]`**, so a
Q93 HiFi base (74% of this substrate's bases) is scored as if it were **Q40**. Hence the smallest
`εⱼ` any single feature can contribute is **1e-4/3 = 3.333e-5** for a PSV column and **`junction_err`
= 1e-4** for a junction, and those two constants are *exactly* the `max p` values in the table above:
3.333e-5 wherever a commitment rests on one PSV, 1.000e-4 wherever it rests on one junction. The worst
ratio in the whole panel, **0.2666, is precisely `3.333e-5 / 1.25e-4`** — one PSV against MAGEA's
9-copy Bonferroni threshold. Three consequences, all mechanical:

* the certificate is **strictly conservative** — it never trusts a base more than Q40 however good the
  instrument's QV is, which is why an independent recomputation using the *raw* QVs lands far below the
  binary's own `p` (below);
* it is **granular**: with `n = 2` a single distinguishing feature certifies; the bound cannot be
  approached continuously from above;
* ⚠ it puts a **family-size ceiling on single-feature evidence**. A lone junction certifies only while
  `1e-4 < α/(n−1)`, i.e. **n ≤ 10 copies**; a lone PSV only while `3.333e-5 < α/(n−1)`, i.e. **n ≤ 30**.
  Beyond those sizes a read needs ≥ 2 distinguishing features no matter how clean it is. The largest
  family here is n = 9, so the panel never reaches the ceiling — **but a real NPIP- or TBC1D3-sized
  family would**, and no earlier document records this. It is a *mechanism* statement (read off the two
  constants), not a sweep: `junction_err` and the `phred_err` clamp both remain **unswept** (§5).

⚠ **That table alone is nearly self-fulfilling and is reported as a consistency check, not as evidence.**
The binary derives `status` *from* `p`, so it can only fail if the gate is miswired. The claim that the
certificate is real rests on the independent recomputation below, which re-derives `p` from raw evidence
and never reads the binary's value.

**Cross-check, independent of the binary's own arithmetic** (`cert_check.py`). For 160 sampled assigned
reads (SHARP + APOBEC3) the certificate was rebuilt from raw evidence — per-copy PSV alleles and the
read's observed base from `--dump-psv`, and the read's **own per-base QV read out of the BAM** at each
PSV genome position — through an independently written exact Poisson-binomial DP. The binary's
`p_value` is never an input. Result: **160/160 clear `thr` on PSV evidence alone** (junction terms
dropped, so this is the conservative half of the certificate), worst recomputed `p/thr` = **2.3e-26**.
On the 6 sampled reads whose emitted `p` carried no junction term either, the recomputation agrees with
the binary within one order of magnitude on **6/6**. Where they differ the binary is the **more
conservative** of the two — e.g. 4.2e-34 against a recomputed 2.9e-114 on the same 17 columns — and the
reason is now known rather than guessed: the recomputation uses the read's **actual** QVs (Q93 ⟹
εⱼ ≈ 1.7e-10) while the binary floors every base at Q40 through the `phred_err` clamp above. **The
claimed bound is real, it is not tight, and the gap is a deliberate conservatism, not an error.**

### 11.3 K = 0 — the part no method can resolve

`tied` = `min_p ≥ thr`: even perfect support at every distinguishing position could not clear the
threshold. **45 / 4,885 = 0.0092 [0.0069, 0.0123]** of evaluated reads; **13 / 488 = 0.0266
[0.0156, 0.0450]** inside MAPQ == 0 — roughly **3× denser** in the hard stratum, which is what a
substrate floor should do. At the row level 977 of 7,691 (read, family) rows are tied, **950 of them
with `n_decisive == 0`** (nothing distinguishing was observed at all) and 27 with evidence whose bound
is still unattainable.

### 11.4 ⭐ What actually drives abstention on real data — and it is NOT the certificate

Splitting `ambiguous` by mechanism (the row's own `p` and `margin` against its `thr`):

| mechanism | rows | fraction |
|---|---|---|
| **(E)** evidence-limited — the representative record's certificate misses `thr`, or `margin == 0` | **2** | 0.0004 |
| **(C)** cross-record contradiction — the record *cleared* the gate, but ≥2 of the molecule's alignment records were assigned to **different copies**, so the molecule abstains (§2a) | **4,615** | **0.9996** |

**99.96% of O2's abstentions on this substrate are the molecule-level contradiction rule, not weak
evidence.** The mechanism was verified directly rather than inferred: reconstructing each record's
`best_overlap_copy` from the BAM and the catalog spans, **4,615 / 4,615** of the (C) molecules do place
records on ≥2 distinct copies, and **0** of them have only one record in the region. The per-read
certificate is essentially never the binding constraint here.

**Measured cost of that rule** (§2a says its cost is measured, not assumed): if the contradiction rule
were dropped and the strongest record simply believed, commitment would rise by **2,975 / 5,378 =
0.5532 [0.5399, 0.5664]** overall (0.3466 → 0.8998) and by **459 / 497 = 0.9235 [0.8968, 0.9438]** on
MAPQ == 0 (0.0302 → 0.9537). ⚠ **That is an abstention cost, not an error rate.** There is no truth
here that says which side of a contradiction is right — which is exactly why the conservative choice
was made. It is also why "believe the strongest record" is not a free upgrade: §2a records that on this
same panel arbitrating by margin coincides with simply believing minimap2's primary flag on 55.9% of
molecules, i.e. it partly re-imports the defect that retired `uniq_agree`.

### 11.5 Why reads were never evaluated

**493 / 493 (100%)** of the `unevaluated` reads have a primary alignment that overlaps **no catalog copy
span** — they sit in the swept box but outside the family. Not one read whose primary lands on a copy
failed to receive a row. `unevaluated` is therefore purely an O1 extent statement and carries no O2
information; it is retained in the denominator only so that the abstention fractions cannot be inflated
by quietly removing it.

### 11.6 The one-sentence version, with its guard rail

> On gorilla control-panel families, O2 commits **34.7% [33.4, 35.9]** of the reads with a primary
> alignment in the region and **3.0% [1.8, 4.9]** of the MAPQ-0 reads that are its actual remit; every
> commitment carries a Bonferroni-corrected certificate at α = 1e-3 that is met on **2,097/2,097**
> assigned rows with a worst-case `p/thr` of 0.2666 and survives independent recomputation from the
> BAM's own quality strings on 160/160 sampled reads; **0.9% [0.7, 1.2]** of evaluated reads (2.7% of
> MAPQ-0 reads) are `tied`, the identifiability floor no method can cross; and **99.96%** of the
> abstentions come from a molecule's own alignment records contradicting each other, not from the
> certificate. **None of this is an accuracy — no read here has a known true copy.**

---

## 12. One-paragraph summary for a reader in a hurry

O2 assigns one long read to one paralog copy inside a family that O1 has already fixed. The copies'
differences are PSV columns (copy-vs-copy alignment, confirmed by a read pileup) and copy-specific
junctions. Each read gets a per-copy log-likelihood over the features it spans; the decision is the
**argmax**, and because the objective contains no term coupling reads or copies, that per-read argmax
**is** the global optimum — there is no combinatorial optimization at assignment time. A read is
reported only if an exact Poisson-binomial certificate clears `α/(n−1)` against its hardest competitor
(α = 1e-3, at the lower edge of a measured [1e-3, 1e-1] plateau) and its argmax is strict; otherwise it
**abstains**, distinguishing "no molecule of this length could have resolved this" (`tied`, the K = 0
floor, computed not asserted) from "this one didn't" (`ambiguous`). Weight is never split 1/k. The unit
is the **molecule**, not the alignment record. Measured over 1.65 M read-rows from 171 runs the rule
above reproduces the emitted status on 99.9990% of the rows a reader can decide from the file; the
assigned-fraction on that corpus is withdrawn (§9.1) because it was counted per record.
