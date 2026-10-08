# DP chaining as the pair score — tested, and it cannot help at this level

Run 2026-09-20 against `docs/PREREG_chained_jaccard_2026-09-20.md` (committed `62196fc2` before any
chained arm was scored). Tool: `bench/chained_pair_score.py`. Same truth (Soto S1C ≥ 3 members), same
scorer, same held-out chromosomes (chr2, chr8, chr10) as §6t3.

## Result

| arm | pooled held-out F | sens | prec |
|---|---|---|---|
| A — shipped: MCL + exon conjunct | **0.7016** | 0.808 | 0.722 |
| A− — MCL alone | 0.6394 | 0.814 | 0.586 |
| **B — global Jaccard, components** (the thing being improved) | **0.6044** | 0.713 | 0.639 |
| **D — chained, containment normaliser (PRE-REGISTERED)** | **0.2132** | 0.500 | **0.164** |
| D2 — chained, symmetric normaliser (POST-HOC) | 0.6044 | 0.713 | 0.639 |

**Pre-registered verdict: ⛔ WORSE** (D vs B, −0.391).

## Two separate failures, and the second one is the interesting one

**1. My normaliser was wrong, and I own it.** Arm D normalised the chain by `min(len_a, len_b)` —
containment — deliberately, to fix r359's short-copy-vs-long-parent penalty. It overcorrected: a short
gene chaining inside a long one scores ~1.0, so long genes become hubs. Measured: **864 chr8 pairs
saturate at ≥ 1.0** and **3,190 at ≥ 0.9**, and components blow up — largest **171 vs global Jaccard's
23** on chr2, **166 vs 49** on chr8. The pre-registration forbade re-tuning after seeing scores, so D
stands as the result.

**2. With that fixed, chaining changes nothing at all — and this is the real finding.** The post-hoc
symmetric variant D2 reproduces arm B to **every digit**, across the whole sweep
(0.2460 / 0.4537 / 0.6059 / 0.6109 / … / 0.6044). The reason, measured on chr8:

| | |
|---|---|
| gene pairs with anchors | 40,645 |
| **best chain uses more than one anchor** | **34 (0.08%)** |
| chained score beats the single best record | 34 (0.08%) |

⭐ **minimap2 has already chained.** Each PAF record *is* a maximal collinear chain; the records it
emits for a pair are separated by gaps far too large for a second-level chain to bridge, so the DP
collapses to "pick the best record" — which is exactly what arm B did. **Chaining above the PAF is
structurally redundant with minimap2's own chaining.**

## What this means for the goal

- ⛔ **The Jaccard operator cannot be improved by chaining at the PAF-record level.** Not "did not help
  here" — it is redundant by construction, and the 0.08% says so directly.
- ⭐ **The place a chaining gain could exist is BELOW minimap2, at the anchor (minimizer-hit) level** —
  i.e. inside the aligner, with different chaining parameters from the ones minimap2 chose. That is a
  much larger undertaking than a rescoring pass, and it competes with a tool tuned for exactly this.
  ⚠ The advisor dislikes minimizers, which is where such work would have to live.
- ⭐ The §6t3 conclusion is unchanged and remains the actionable one: **the operator is not where the
  value is** (MCL buys +0.035 over B; the exon conjunct buys +0.062, and only survives inside MCL —
  register 911). Effort spent on scalars between gene pairs is effort spent on the saturated axis
  (§6o8/§6o9).

## Limits

- One gap function (minimap2's form, `0.01·ā·g + 0.5·log2 g`) at one parameterisation. A far more
  permissive γ would join more records — but joining records minimap2 deliberately split is a claim
  about the alignment, not about family definition.
- 27 Soto families over three chromosomes; the 0.08% redundancy figure is chr8 only, though the identity
  of D2's and B's full sweeps implies it holds on chr2 and chr10 too.
