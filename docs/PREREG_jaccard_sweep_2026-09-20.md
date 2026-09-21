# Pre-registration — SIMPLE JACCARD SWEEP vs the shipped rule

**Written 2026-09-20 before running either arm on the scored substrate.** User: *"my advisor insists we
are over-complicating, a simple jaccard sweep should work."*

Taking the claim seriously and testing it head-to-head, not arguing it.

## What is already measured, on both sides

- ⭐ **r772 supports him at the EDGE level.** "Jaccard/MinHash cannot replace the pairwise alignment for
  O1's edges" was **REFUTED on measurement** (`bench/jaccard_vs_align.py`, NPIP panel, 3,128 spans): only
  3 of 1,258 O1 edges fall below Jaccard 0.01; **t = 0.05 recovers 94.1%**, t = 0.01 recovers 99.8%. The
  reason the objection failed: O1 already demands `cov_longer ≥ 0.30`, so both are global-coverage
  criteria and they agree. ⚠ What Jaccard cannot do is the **core step**, which needs per-base
  coordinates per partner, not a scalar.
- ⚠ **Against him:** r293 (domain-sharer CREB1~METTL21A J = 0.406 **above every true paralog**, max
  0.313 — the ranges overlap), r359 (Jaccard penalises short-copy-vs-long-parent; EEF1A1 retrocopies
  J ≈ 0.07 at core identity 0.92), and **§6o6** (plain identity-threshold components are consistently
  worse than MCL + the shared-exon conjunct, LORO over 3 regions).
- ⚠ **And a reason the test may be uninformative either way:** §6o8/§6o9 found grouping is **saturated** —
  recall is bounded by the truth's divergence, not the operator. If both arms land in the same place,
  that is the expected result, and it argues for the SIMPLER arm.

## The two arms

**A — shipped (the "over-complicated" one):** `mcl_families --min-exonic-bp 1 --min-shared-exon-frac 0.60`
over the same all-vs-all PAF. Already computed in §6s8; not re-run, not re-tuned.

**B — the advisor's:** connected components of the pair graph keeping every pair with
`J = nmatch / (len_a + len_b − nmatch) ≥ t`, from the same PAF. **No MCL, no exon conjunct, no core
step, no corroboration.** Swept over **t ∈ {0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.50}**.

Both arms are scored by the **same** scorer (`bench/heldout_family_score.py`, bipartite
sensitivity/precision/F) against the **same** external truth (Soto S1C families, ≥ 3 members on the
chromosome), on chr2, chr8, chr10 (held out) and chr16 (development).

## The bar — committed now

Let **F_B\*** be arm B's pooled held-out F at its BEST t, and **F_A** arm A's pooled held-out F (0.7016).

| outcome | verdict |
|---|---|
| F_B\* ≥ F_A − 0.02 | ⭐ **HE IS RIGHT** — the simple sweep matches; recommend simplifying |
| F_A − 0.10 ≤ F_B\* < F_A − 0.02 | ⚠ **CLOSE** — report the gap and what the extra machinery buys |
| F_B\* < F_A − 0.10 | ⛔ **the machinery earns its place** — report the number, not the opinion |

⚠ **Declared before the run:** picking B's best t *after* seeing the scores is fitting B to the test set,
which flatters B. I will therefore report **both** B-at-best-t and B at the t that r772 already fixed
independently (**0.05**), and treat the latter as the honest comparison.
⚠ I will not re-tune arm A, drop a chromosome, or change the truth.
