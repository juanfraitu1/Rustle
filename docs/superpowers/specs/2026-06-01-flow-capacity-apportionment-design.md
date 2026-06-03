# Mass-conserving, anchor-aware multimapper apportionment for trustworthy max-flow capacities

**Date:** 2026-06-01
**Status:** Design APPROVED (open questions resolved 2026-06-01) — ready for implementation plan
**Scope:** VG / `--vg`-gated only. Default de-novo (95.6/90.5) must stay byte-identical.

## Problem

In VG mode, multimapping long reads (paralog copies, especially **inverted** pairs like DAZ1(−)/DAZ3(+)) inflate per-copy abundance because secondaries are retained at uncorrected weight and that weight feeds the max-flow capacities directly. The chain is **exact-given-capacities**:

```
read.weight (=1/NH, bam.rs:629)
  ├─ node.coverage  += weight*bp           (map_reads.rs:331; graph.rs:121)
  └─ tf.abundance   += weight              (map_reads.rs:1855)
        → nodecov = max(abundin,abundout)  (global_flow.rs:180)
        → capacity[i][end] += tf.abundance (max_flow.rs:1280)   ← tf.abundance IS the edge capacity
        → Edmonds-Karp flux                (max_flow.rs:134/1487)
        → ecov = nflux * noderate          (path_extract.rs:8732)
        → Transcript.coverage              (path_extract.rs:9119-9123)
```

The max-flow is exact, so **wrong capacities ⇒ biased (not merely noisy) abundances** — ~2× for a 2-copy double-count. Empirically on DAZ (215 reads on both copies, anchored 177:1): uniform 1/n → DAZ3=107.5; uncorrected full inclusion → 215; mass-conserving EM → ~1.2 (truth ~2–6).

**Two distinct defects:**

1. **Apportionment gap.** `run_fingerprint_em` (vg.rs:4467) already softmax-normalizes a read's placements to sum to 1.0 (vg.rs:4757-4772, *mass-conserving*), but the M-step prior is **pileup-depth** (vg.rs:4694-4697, `log_priors = ln(copy_total/total + 1e-3)`), which mixes anchored and ambiguous mass and can **self-reinforce an already-double-counted copy**. The 0-diagnostic-site case (vg.rs:4531, DAZ) relies entirely on this prior.
2. **Inverted-pair gap.** `partition_and_remap_family_by_strand` (vg.rs:2242, called unconditionally at pipeline.rs:10667) groups `bundle_indices` by strand and **drops any multimapper that retains <2 placements after the split** (vg.rs:2283-2287). DAZ1(−)/DAZ3(+) are different loci (not caught by `is_strand_mirror`'s ≥90% co-location test, vg.rs:247-257), so a read shared only between them keeps 1 placement per strand-partition → dropped from both → never reweighted → stays at 1/NH in **both** copies → double-counted into both capacities.

No per-transcript capacity confidence exists: `copy_assignment_confidence` (vg.rs:1434-1464) is per-copy, read-count-averaged, multimap-only, and returns 1.0 when `n_multimap==0` (conflates "no multimappers" with "fully resolved"); `copy_independent_support` defaults to 1.0 on missing data.

## Capacity-coupling analysis (the core concern)

`read.weight` **is** the capacity unit, and the coupling is the point. The same scalar feeds both the noderate numerator (`coverage = Σ weight*bp`) and the denominator/capacity (`tf.abundance → nodecov`). For an isolated single copy, `noderate = coverage/nodecov` is invariant to uniform weight scaling, but **absolute capacity and thus flux scale linearly**. A multimapper counted in two copies at uncorrected 1/NH inflates each copy's capacity independently ⇒ each copy's flux/coverage independently inflated = the double-count.

Because BFS gates on `residual > EPSILON` (max_flow.rs:155), not magnitude, the **isoform path set is invariant** to positive capacity scaling — so a correct upstream weight fix changes *abundance numbers* without changing *structure* in the common case. Two structural-corruption escape hatches exist but are **latent (off by default)**: (a) abundance-ordered seed depletion (path_extract.rs:5885/5895, env-gated OFF; default is insertion order 5904-5906), and (b) gate-crossing (`min_cov_gate` path_extract.rs:8914, checktrf `readthr` path_extract.rs:42) — an inflated copy can pass a gate a true-weight copy fails.

**Therefore the fix is upstream:** correct `read.weight` *before* graph build. Everything downstream is exact given capacities, so this propagates with **zero flow-engine changes**. The in-flow alternative (scaling `max_flow.rs:1280`) is rejected — it would force consistent rescaling of `noderate` (1304-1306) and the reconversion at path_extract.rs:8732.

## Chosen design

- **DEFAULT (ships ON in VG mode): anchor-first one-pass apportionment** (panel 6/10) as the weight corrector, **plus the capacity-confidence output layer** (panel 5/10) as honest annotation.
- **OPT-IN (env, OFF by default): coupled EM↔flow fixed-point** (panel 4/10) for weak-anchor families only.
- **Shared prerequisite for both: the inverted-pair joint-EM-input fix.**

**Rationale.** Anchor-first is the highest-scored corrective and is the spec's primary (upstream-weight) intervention. It removes the M-step self-reinforcement by grounding the prior in **unambiguous mass** (unique reads + dNM-decisive via `anchor_read`, vg.rs:1500). For the DAZ 177:1 regime a single softmax with the anchored prior is equivalent to the converged EM. Its one weakness — the cold-start where *all* copies have zero anchors collapsing toward uniform 1/n — is covered by graceful-degrade and the opt-in coupled fallback. Coupled EM↔flow is the most principled (capacities become a fixed point of the flow they feed) but is expensive (N_outer × full assembly), can limit-cycle, and is degenerate exactly on the 0-site DAZ case — so it is opt-in, bounded to ≤4-copy families. Uncertainty-range does not fix the bias by itself (point estimate unchanged), so it is an add-on; its capacity-confidence + abstain + jointly-feasible MIN band are high-value and low-risk.

## Components

1. **Inverted-pair joint EM input (shared prerequisite).** New `vg.rs::family_for_em_input(family, bundles) -> FamilyGroup` returning the **unsplit both-strands** family for the EM **placement list + normalization** only. Wiring (pipeline.rs:10665-10669): build `em_input_partitions` (unsplit) separately from `em_hmm_partitions` (strand-split). `fam_pos` preserved as index into `bundle_indices` so write-back `global_bi/ri` (vg.rs:4852-4857) stays valid. Gate `RUSTLE_VG_JOINT_STRAND_EM` (default ON in VG).
   - **O1-resolved constraint:** `build_family_graph` hard-bails on mixed strands (family_graph.rs:432-433) and `enumerate_diagnostic_sites` assumes single-strand co-linear sequence, so the family graph is **NOT** rebuilt on the unsplit family. Instead: keep per-strand `family_graphs` (so the `fp` diagnostic-site term stays valid *within* a strand sub-family); set the `fp` term **neutral (0) for cross-strand placements** (which is already the reality for inverted near-identical copies — `fp.n_sites==0`, vg.rs:4531). The joint normalization then apportions cross-strand mass using the **strand-agnostic** signals only: `junc_scores`, `nm_scores`, and the anchored prior (all per-placement, computed from each alignment independently). Net: the fingerprint keeps working where it is valid; the cross-strand double-count is fixed by junction+identity+prior + the joint log-sum-exp.
2. **Anchored-mass prior (anchor-first core).** New `vg.rs::anchored_mass_per_copy(family, bundles, t, extent_frac) -> Vec<f64>`: per copy, sum `read.weight` of reads that are unique (`read_name_hash ∉ multimap_reads`) OR `anchor_read(...)==Owns` (vg.rs:1500, raw-dNM margin). Replace the M-step prior source (vg.rs:4688-4697): `log_priors[k] = ln(anchored_mass[k]/total + 1e-3)` computed once; run `max_iter=1` unless coupled mode. Phase 1/2 and normalization reused verbatim. Gate `RUSTLE_VG_ANCHOR_PRIOR` (default ON in VG).
3. **Persist read anchoring.** Add to `BundleRead` (types.rs:~142): `em_weight_gap: f64` (default −1.0), `em_n_sites: u32` (0), `em_anchored: bool` (default `nh<=1`). Set in write-back (vg.rs:4823-4866) where gap (4833) and max_sites (4834) already exist: `em_anchored = (gap>0.8) || (max_sites>0 && gap>0.5) || was-unique`.
4. **Capacity-confidence channel.** Add `anchored_coverage: f64` to `GraphNode` (graph.rs:111-174); accumulate in parallel with `coverage` (`+= em_anchored ? weight*bp : 0`). Per-tx `capacity_confidence = Σanchored / Σtotal` over path exons. Add `capacity_confidence: Option<f64>` and `abundance_min: Option<f64>` to `Transcript` (beside `copy_independent_support`).
5. **Coupled fallback (opt-in).** `RUSTLE_VG_COUPLED_EM` (OFF). When set AND ≥2 copies have ~0 anchor, run ≤N=3 damped outer iterations feeding flow-derived per-copy abundance back as prior. Bounded `n_copies<=4`.
6. **GTF emission.** gtf.rs:177-185 add `capacity_confidence` and `abundance_min` attrs (mirror `copy_confidence`), emitted only when `Some`.

## Algorithm

**Phase A — joint-strand EM input (pipeline.rs ~10665):** if `JOINT_STRAND_EM`, `em_input_partitions = map(family_for_em_input)` (unsplit); keep `em_hmm_partitions` (split) for graph/borrow. Build `family_graphs` against the EM input's copies so `fam_pos` aligns (see O1).

**Phase B — anchor-first apportionment (run_fingerprint_em):** Phase 1 (4548-4578) and Phase 2 (4580-4585) unchanged. Compute `anchored = anchored_mass_per_copy(...)` once; `log_priors[k] = ln(anchored[k]/total + 1e-3)`. If `total_anchored < eps` (all-zero-anchor): coupled mode → Phase D, else **graceful-degrade to the existing pileup prior for that family only**. Run the E-step once (4700-4773) with fixed prior; keep the score-gap gate and log-sum-exp normalization verbatim (mass sums to 1.0). Write-back sets weight + `em_weight_gap`/`em_n_sites`/`em_anchored`.

**Phase C — capacity-confidence (capacity build):** wherever `node.coverage += weight*bp` (map_reads.rs:331/562; bundle cov adds), also `node.anchored_coverage += (em_anchored?weight:0)*bp`, gated on `vg_mode`. At Transcript push (path_extract.rs:9119-9123), `capacity_confidence = (Σanchored/Σtotal).clamp(0,1)`, `abundance_min = coverage*capacity_confidence`.

**Phase D — coupled fallback (opt-in):** outer loop ≤3: E-step → write weights → build capacities + run flow for the family's bundles → next prior `ln(A_k/Σ + eps)` with damping `0.5`; stop at `max|ΔA_k|<0.01` or N. Bounded `n_copies<=4`.

**Phase E — confidence attach + abstain (pipeline.rs:18104-18110):** set `tx.capacity_confidence`, `tx.abundance_min`. If `< RUSTLE_VG_ABSTAIN_FLOOR` (0.05): tag low-confidence via `family_verdict` (**do not drop** — keep coverage so TPM/benchmark work). GTF emit when `Some`.

## Conservation guarantee

Per-read mass conservation preserved by reusing the unchanged log-sum-exp normalization (vg.rs:4757-4772): every multimapped read's placements sum to 1.0. Anchor-first changes only the prior + iteration count. The joint-strand input makes conservation span **both** strands (the read's 1.0 splits across DAZ1+DAZ3 instead of being 1/NH in each).

**Conservation holes (inherited, documented):** reads with <2 placements keep 1/NH (but joint-strand removes the biggest source); score-gap-gated reads keep their prior weight — to close this, when the gate fires, **explicitly set the read's per-placement weights to the normalized anchored prior** (mass-conserving fallback) rather than leaving raw 1/NH (small addition at vg.rs:4752). Coupled mode is a damped contraction (bounded N=3, d=0.5) to avoid starving weak-but-real copies. `abundance_min = coverage*capacity_confidence` is sub-conservative (correct lower bound); **no per-copy MAX is emitted** (non-additive).

## Capacity-confidence output

`capacity_confidence = anchored_coverage_sum / total_coverage_sum` over the transcript's path nodes. Strictly better than `copy_assignment_confidence` because it is per-transcript (not per-copy), bp-weighted (capacity units), and includes unique reads in the denominator. **Band:** only the jointly-feasible MIN is shipped (`abundance_min`); the point estimate stays primary. **Abstain:** `< floor` ⇒ tag low-confidence, keep coverage.

## Gating

All VG-gated; default de-novo byte-identical (entire VG family pipeline only runs under `config.vg_mode`; new fields written only in VG mode, emitted only when `Some`).
- **Default-ON in VG:** `RUSTLE_VG_JOINT_STRAND_EM=1`, `RUSTLE_VG_ANCHOR_PRIOR=1` (`=0` restores current VG behavior).
- **Default-OFF:** `RUSTLE_VG_COUPLED_EM=0` (+ `_ITERS=3`, `_DAMP=0.5`), `RUSTLE_VG_ABSTAIN_FLOOR=0.05`.
- No edits to `max_flow.rs`; default seed sort (path_extract.rs:5904-5906) untouched. `RUSTLE_VG_DROP_SECONDARY` remains the full escape hatch.

## Test & validation plan

1. **Regression guard (first):** default de-novo on GGO_19 with no `--vg` → 95.6/90.5, GTF byte-identical to HEAD. Assert `anchored_coverage` stays 0 when `!vg_mode`.
2. **DAZ headline (primary success):** `rustle --vg --genome-fasta …` on DAZ → DAZ3 coverage **~163 → ~1–2**; DAZ1 not over-deflated. Via `RUSTLE_VG_FP_ATTR_TSV` confirm shared reads appear in ONE PreEntry with `weight_sum==1.0` across both strands (proves joint-strand fix). Compare against `JOINT_STRAND_EM=0`.
3. **Synthetic oracle (unchanged):** fingerprint-EM synthetic (commit 9436e7e) stays 100% — anchor-first 1-pass on decisive fingerprints must still be exact.
4. **Mass-conservation unit test:** per multimap read, Σ weights across placements == 1.0 ± 1e-6; dedicated gate-fired-read test.
5. **Genome-wide scan (no fabrication regression):** re-run `bench/paralog_secondary_scan`; 93 expressed_real_copy retained, 89 pure_spillover get LOW capacity_confidence (band/abstain fires).
6. **Confidence sanity:** `0≤cc≤1`, `abundance_min≤coverage`, all-unique copies `cc~1`, DAZ3 `cc` low.
7. **Coupled opt-in:** convergence within N, no limit cycle; DAZ (0-site) excluded from coupled mode.

## Risks

1. **Build churn** — ~20 `Transcript` literals + `GraphNode`/`BundleRead` literals need the new fields; mitigate with `Default` impls; compiler enforces exhaustiveness (safe failure).
2. **Joint-strand fam_pos / graph alignment (O1-resolved, now lower risk)** — since the graph is kept per-strand (not rebuilt on the unsplit family), the `fp` term must map each unsplit placement to its strand-local graph copy index, or fall back to neutral when no same-strand graph copy exists. Mitigate: index `fp` by `(strand-subfamily, local copy)` and default to 0 cross-strand; `junc`/`nm`/prior carry the rest. `JOINT_STRAND_EM=0` = instant rollback to current behavior.
3. **Over-deflation of weak-but-real copies** — eps floor (never exactly 0); coupled opt-in/damped/bounded; genome-wide-scan is the false-negative guard.
4. **Structural corruption via gate/seed crossing** — corrected lower capacities may drop a copy below a gate, but that is *correct* if it existed only via double-count; abundance seed-sort stays OFF; check in genome-wide scan.
5. **Confidence-channel desync** — `anchored_coverage` must update everywhere `coverage` does; clamp at attach + per-node assert `anchored_coverage<=coverage`.
6. **Band/abstain semantics** — only jointly-feasible MIN shipped; no MAX.

## Decisions (resolved 2026-06-01)

- **O6 → Default-ON in VG mode.** `RUSTLE_VG_JOINT_STRAND_EM` and `RUSTLE_VG_ANCHOR_PRIOR` ship default-ON in `--vg`; `=0` restores prior VG behavior. Default de-novo (95.6/90.5) unaffected.
- **O2 → Graceful-degrade + flag** for 0-anchor families: fall back to the existing pileup-depth prior (status-quo point estimate) AND set low `capacity_confidence` / low-confidence tag. Coupled EM stays opt-in.
- **O5 → Tag, keep coverage** for abstained transcripts: low-confidence GTF attribute, transcript SET unchanged (TPM/benchmark preserved).
- **O3 → Single E-step** (`max_iter=1`) with the fixed anchored prior (panel recommendation; sufficient for decisive anchors).
- **O4 → Reuse project copy-anchor calibration** (raw dNM, `t=2`, `extent_frac=0.8`); expose as env later only if needed.
- **O1 → Resolved in code** (see below): per-strand graphs kept; cross-strand apportionment via junction + identity + anchored prior.

## Original open questions (for the record)

- **O1 — RESOLVED (verified in code):** `build_family_graph` bails on mixed strands (family_graph.rs:432-433); `enumerate_diagnostic_sites` assumes single-strand. So we do **not** rebuild the graph on the unsplit family. Resolution: keep per-strand graphs for the `fp` term (valid within strand), neutralize `fp` across strands (already true for inverted copies, `fp.n_sites==0`), and let joint normalization apportion cross-strand mass via junction + identity + anchored prior. No user decision needed.
- **O2:** for the 0-anchor case, default = graceful-degrade to existing pileup prior (status quo), with abstain flag firing. Acceptable, or prefer full-abstain reporting?
- **O3:** anchor-first single E-step (`max_iter=1`, recommended) vs a few iterations with the fixed prior. Default: single-pass.
- **O4:** reuse project copy-anchor calibration (raw dNM, t=2, extent_frac=0.8) vs expose as env from day one.
- **O5:** abstain action = tag (keep coverage, recommended) vs drop (changes GTF transcript set, risky).
- **O6:** ship `JOINT_STRAND_EM` + `ANCHOR_PRIOR` default-ON in VG (changes current `--vg` output) vs default-OFF for one release to A/B first.
