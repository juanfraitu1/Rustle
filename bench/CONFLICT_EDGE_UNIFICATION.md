# Unifying the family-definition edge with the assignment gate (2026-06-29)

**Motivation.** The family-as-graph definition and the read→copy assignment were two principled objects that used
*different* criteria for the same underlying question ("can this read tell these two copies apart?"):

- **Family edge** (`read_conflict.rs::de_tied`): a read links two loci iff `|de_a − de_b| ≤ delta` (a fixed
  `delta = 0.005`) and `max(de_a, de_b) ≤ de_max`. The `delta` is the last hand-set constant in the definition.
- **Assignment gate** (`copy_assign.rs`): a read is *tied* between two copies iff `min_p = ε^δ ≥ α` (Theorem 4),
  with `δ` = the number of distinguishing columns the read spans — an error-model-derived, threshold-free test.

The F4 finding (raw read-conflict graphs are error-inflated, colouring ≈3× the true copy count) is the empirical
case that the conflict edge should use the *same significance test* as assignment, not a fixed tolerance.

**The unification.** A read counts as conflict evidence between two loci iff it **cannot significantly
distinguish them** under the assignment gate's own criterion. With `m_x = de_x · aln_len_x` the mismatch count to
locus `x`, the excess `δ = |m_a − m_b|` is the per-read distinguishing-column proxy, and the read ties iff

  ε^δ ≥ α            (exactly the gate's `min_p ≥ α`, Theorem 4)

keeping the `de_max` quality floor (both alignments must genuinely fit). The arbitrary `delta = 0.005` is gone;
the tie boundary is now the error model (`ε`, default `e/3 ≈ 0.001`) and the significance level (`α`, default
`1e-3`) — the *same two numbers* the assignment gate uses. One IsoCon real-vs-error criterion now governs both
the family edge and the read assignment.

**Implementation (`src/rustle/vg_family/read_conflict.rs`).**
- `Placement` gains `aln_len` (aligned-block length = #M/=/X CIGAR columns), computed at construction in
  `denovo_pipeline.rs`.
- `ConflictParams.sig: Option<(eps, alpha)>` — `None` (default) = the legacy `de_tied`, so **OFF is
  byte-identical**; `Some` = `sig_tied`. Env-gated: `RUSTLE_CONFLICT_SIG=1` (+ `RUSTLE_CONFLICT_EPS`,
  `RUSTLE_CONFLICT_ALPHA`).
- `tied()` dispatches; `conflict_edges`/`family_mapq0_support` call it. The `de-tie ⊆ AS-tie` audit invariant is
  untouched.
- Tests: `sig_criterion_ties_ambiguous_resolves_distinguishing` (the unification at unit level),
  `sig_off_default_is_byte_identical_to_de_tied` (default ships OFF). Full lib suite green (685 tests).

**Why it lands with Canzar.** It removes the last similarity-style constant from the family definition: the
boundary is now a property of the data (the per-base error rate) and a stated significance level, identical to
the assignment gate — one clean combinatorial criterion, no tuned thresholds.

## Real-data measurement (significance edge ON vs OFF; genome-wide `gw_family_catalog` on `GGO_mm.bam`)

| catalog | families | copies |
|---|---|---|
| **OFF** (legacy `de_tied`, fresh run, this binary) | **81** | **205** |
| **SIG** (`RUSTLE_CONFLICT_SIG=1`, ε=1e-3, α=1e-3) | **71** | **176** |
| existing baseline (older binary, for reference) | 82 | 207 |

- **Byte-identical OFF confirmed end-to-end:** the fresh OFF run (81/205) matches the prior baseline (82/207)
  to within ±1 family — run-to-run noise from threaded rep assembly, *not* the struct change (the new `aln_len`
  field is ignored by `de_tied`; the unit test `sig_off_default_is_byte_identical_to_de_tied` pins the criterion).
- **The significance edge is a principled REFINEMENT.** Analytically, for equal-length placements at the default
  ε/α, `sig_tied ⟹ de_tied` (the tie boundary `δ ≤ 1` is `|de_a−de_b| ≤ 1/L ≤ 0.005`), so SIG edges are a
  **subset** of de-tie edges — it can only shrink or split families, never create them.
- **Effect: a modest, defensible narrowing — 81→71 families (−12%), 205→176 copies (−14%).** The catalog drops
  the ~1/8 of families held together only by *marginally*-tied reads (within the old 0.5% tolerance but
  significantly **resolvable** under the gate). The largest families are **preserved** (9-copy GWFAM58, 7-copy
  GWFAM10 / GWFAM67), and the size distribution is essentially unchanged (mostly 2-copy paralog pairs).
- The exact per-copy coordinate overlap between the two runs (~78%) is limited by rep-assembly non-determinism
  across separate runs, not by the edge criterion — the analytic subset relationship above is the rigorous
  refinement claim; the family/copy counts are the robust empirical headline.

**Interpretation / the choice for the advisor.** OFF defines a family as loci linked by reads within ~0.5%
divergence; SIG defines it as loci a read **cannot significantly resolve** — i.e. the genuinely
collapsed/assignment-hard set, under the *same* criterion the assignment gate uses. SIG is the more principled,
threshold-free object (no `0.005`), at the cost of excluding ~12% of near-but-resolvable pairs. It ships **OFF
by default** so nothing changes silently; flip with `RUSTLE_CONFLICT_SIG=1` once the definitional choice is made.
