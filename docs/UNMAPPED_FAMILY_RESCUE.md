# Why family-similar reads can still be "unmapped"

A common worry when proposing to rescue unmapped reads against a gene-family
model: *if a read has any resemblance to known family copies, shouldn't the
aligner have found it already?* The answer is no — long-read aligners
(minimap2 in particular) use **seed-chain-extend with hard thresholds**, not
pairwise similarity. A read can carry substantial family-shared sequence and
still fail at any of several stages.

## Failure modes

| Failure mode | What happens in minimap2 | Why a k-mer / family-HMM scan still finds it |
|---|---|---|
| **1. Below-threshold chain score** | Read has matching k-mer islands separated by indels or divergence. Each island is too short to anchor a primary chain (minimap2's `-s` minimum chainable score, default ≈ 40 for HiFi). | Family scan sums *all* k-mer hits ignoring chain geometry. Fifty scattered 15-mer hits add up to a strong family signal but don't chain into a single alignment. |
| **2. Repetitive-seed masking** | Minimap2's `-f` flag drops the top X% most-frequent seeds. In a multi-copy family, those are exactly the *family-conserved* seeds — the most informative ones for family attribution. | The family k-mer index keeps every family k-mer, abundant or not. Profile-HMM emissions weight conserved columns *up*, not down. |
| **3. Divergence above the affine-gap budget** | A novel copy at ~85% identity with several copy-specific indels exceeds minimap2's score-per-base budget at the alignment-extension stage even after seeding succeeds. The candidate is dropped. | HMM forward scoring sums log-probabilities over all paths; there is no hard "give up" threshold, and indels are first-class transitions, not penalties against an affine model. |
| **4. Copy-specific structural difference** | A read spanning a copy-specific exon insertion gets fragmented into supplementaries with no primary alignment, or rejected by primary-alignment fraction filters. | A per-exon profile-HMM family graph allows the read to traverse a copy-specific branch at a bubble — the structural difference is modelled, not penalised. |
| **5. Reference is wrong (genome is an average)** | The true source copy is not in the linear reference at all (the assembly is a consensus over individuals; some paralogs exist in the donor but not the reference). Minimap2 has no positive anchor. | Cross-family conserved emissions still light up at conserved positions; the HMM lands the read on the closest paralog branch, and the divergence shows up as low per-column posterior — a *signature* of a missing copy. |

## The two questions, side by side

- **Aligner asks:** does a high-scoring chain exist between this read and the linear reference?
- **Family HMM asks:** what is the probability this sequence was emitted by *some* path through the family graph?

The second question is strictly more permissive in exactly the regimes that
matter for multi-copy families: high paralog-internal divergence,
copy-specific structural variation, and reference bias.

## Validation: per-read failure-mode classification

Rescue is only a defensible scientific claim if we can show *which* of the
five failure modes each rescued read fell into. Rustle's rescue pipeline
classifies every rescued read, distinguishing:

- **Buckets 1–2 (below-threshold / seed-masked):** the aligner could have
  found these with relaxed parameters. The HMM is acting as automatic
  re-parameterisation. Useful but not novel.
- **Buckets 3–5 (divergent / structural / reference-absent):** no minimap2
  parameter setting recovers these. The family HMM is doing irreducible
  work; these are the real novel-copy candidates.

For buckets 3–5 specifically, Rustle re-runs minimap2 as an external
subprocess with stepped-down parameters as verification — turning the
"this could not have been found by the aligner" claim into a reproducible
fact rather than a model assumption.
