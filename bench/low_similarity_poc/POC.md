# POC — graph-Viterbi assignment for low-similarity paralogs

**Date:** 2026-05-07
**Status:** mechanism demonstrated in unit tests; not yet wired into the EM pipeline.

## Goal

Extend the per-copy HMM-EM machinery to the **low-similarity band** (paralog
jaccard < 0.30) without disturbing the medium/high-similarity path that
already works on AMY and NBPF. The user's intuition: HMM should still be
useful, and aligning to the **graph itself** (not a specific paralog's path)
is the natural relaxation.

## What this POC adds

A single new function in `src/rustle/vg_hmm/scorer.rs`:

```rust
pub fn assign_via_graph_viterbi(fg: &FamilyGraph, read: &[u8])
    -> Vec<(CopyId, f64, f64)>;   // (paralog, recall, precision)
```

It runs the existing `viterbi_path` (unconstrained Viterbi over the full
family graph — same M/I/D trellis from slide 5, no new math) and then
scores each known paralog by **node-overlap** between its path and the
Viterbi trace. Recall = `|viterbi ∩ paralog_path| / |paralog_path|`,
precision symmetric. Sorted by recall descending.

The point: the HMM trellis stays exactly the same; we just stop forcing
the read to commit to a specific paralog's exact path. At low similarity
that constraint is exactly what hurts — reads with a normal error rate
collect penalty against any one path. Viterbi over the graph picks the
best path through ANY combination of nodes; we read the assignment off
which paralog those nodes belong to.

## Tests (in `vg_hmm::scorer::tests`)

Synthetic 4-node, 2-paralog graph with paralog-specific middle nodes
(deliberate low DNA jaccard between A and B):

```
              ┌─→ N1 (paralog A's middle: A-rich) ─┐
N0 (shared) ──┤                                    ├─→ N3 (shared)
              └─→ N2 (paralog B's middle: GC-only) ┘
```

### Test 1 — clean read traces paralog A's path

```
[POC] graph-Viterbi assignment on a low-similarity 2-paralog family:
[POC]   paralog A  recall = 1.000  precision = 1.000   ← assigned
[POC]   paralog B  recall = 0.667  precision = 0.667
[POC]   gap (recall) = 0.333   ← decisive
```

Recall gap of 0.333 because B shares only the 2 of 3 nodes with A's path
(N0, N3 — both shared exon-classes; N2 ≠ N1).

### Test 2 — sanity: read traces paralog B → B wins

Mirror of test 1. Confirms the assignment isn't biased toward CopyId 0.

### Test 3 — noisy read (≈15 % substitution rate) traces paralog A

```
[POC noisy read]  read = paralog A's path with ~15% per-base errors:
[POC]   per-copy forward log P:  A =  -40.79   B = -102.21   gap = 61.42
[POC]   graph-Viterbi recall:    A =   1.000   B =   0.667   gap = 0.333
```

Both methods favor A here. But:

- **Forward log P gap is in nats** — its absolute value depends on read
  length, error rate, and the per-base emission penalty. Hard to choose a
  threshold that works across families.
- **Viterbi recall is bounded in [0, 1]** — the same threshold (e.g.
  `recall ≥ 0.6` and `recall − second ≥ 0.2`) works across families
  regardless of read length or divergence.

This is the property that makes graph-Viterbi a reasonable score-gap
analogue for the low-similarity band.

## Where this would slot into the existing EM

Not yet implemented. The dispatch shape:

```rust
match family_similarity_band(fg) {
    Band::HighOrMedium => forward_against_path_for_copy(fg, read, path_c, c),
    Band::Low          => assign_via_graph_viterbi(fg, read),
}
```

Then the M-step still aggregates per-copy weights and produces priors
exactly as on slide 6 — the only change is *how the per-(read, copy)
score is computed* in the E-step.

## Honest caveats

1. **Family discovery first.** This POC assumes the family graph has
   already been built. At jaccard < 0.30 the existing k-mer linker won't
   put the two paralogs into the same family in the first place, so this
   E-step extension is meaningful only after a complementary linking
   signal lands (splice-site Jaccard, conserved-domain anchors, etc.).
2. **Synthetic, not yet on real BAM.** A real demo would need a
   genuinely-low-similarity family in GGO data. We don't have one
   built end-to-end yet because of caveat 1.
3. **Recall ties are possible.** If a read takes a recombined path (uses
   N1 from A AND N4 from B, say), both paralogs may tie at the same
   recall. That's a feature — it's a signal that the read isn't a clean
   instance of either known paralog and probably shouldn't be redistributed
   confidently. The score-gap rule can be re-used here: abstain when
   `top.recall − second.recall < δ`.

## Reproducibility

```bash
cd /scratch/jxi21/Assembler/Rustle
cargo test --profile dev-opt --lib vg_hmm::scorer::tests::poc -- --nocapture
```

Three tests, all passing. Source: `src/rustle/vg_hmm/scorer.rs`,
search for `assign_via_graph_viterbi` and `poc_graph_viterbi_*`.
