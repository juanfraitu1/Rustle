# Pre-registration — DP CHAINING as the pair score, vs global Jaccard

**Written 2026-09-20 before any chained arm is scored.** User: *"improve the jaccard operator to see if
it can be better at finding families, try DP chaining similar to minimap2's."*

## Why there is room

§6t3's arm B scored a pair by **global Jaccard** `J = nmatch / (len_a + len_b − nmatch)` on the pair's
**single best PAF record**. Two measured problems with that:

- **It discards almost all the evidence.** 99.3% (chr2), 99.7% (chr8), 99.9% (chr10) of gene pairs have
  **more than one** PAF record — median 2–4, max 638–999. Arm B looked at one of them.
- **It is a global criterion on a local phenomenon**, which the register already records as Jaccard's
  failure mode: **r293** (domain-sharer CREB1~METTL21A J = 0.406 above *every* true paralog, max 0.313)
  and **r359** (EEF1A1 retrocopies J ≈ 0.07 at core identity 0.92 — short copy vs long parent).

A minimap2-style colinear chain addresses both: it aggregates *all* anchors between the pair and
rewards **collinear structure**, which a domain-sharer does not have and a real duplicate does.

## The operator — frozen

For each unordered gene pair, anchors = its PAF records `(q_start, q_end, t_start, t_end, nmatch)`,
split by strand and sorted by `q_start`. Standard chaining DP:

```
f[i] = w_i + max(0, max over j<i colinear of ( f[j] − γ(gap(j,i)) ))
  colinear:  q_start_i >= q_end_j  and  t_start_i >= t_end_j      (same strand)
  w_i      = nmatch_i
  gap(j,i) = | (q_start_i − q_end_j) − (t_start_i − t_end_j) |     (diagonal difference)
  γ(g)     = 0 for g = 0, else 0.01 * avg_anchor_len * g + 0.5 * log2(g)      (minimap2's form)
```

Pair score = `max_i f[i]`, normalised as **`chain_cov = max_i f[i] / min(len_a, len_b)`** — the chained
matching bases as a fraction of the shorter gene, which is deliberately a *containment* measure rather
than a symmetric Jaccard, since r359's failure was exactly length asymmetry.

Families = **connected components** of pairs with `chain_cov ≥ t`, no MCL and no exon conjunct, so it is
a like-for-like replacement of arm B's scalar. Swept over
**t ∈ {0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70}**.

## Comparators, already measured on this exact substrate and scorer

| arm | pooled held-out F |
|---|---|
| A — shipped: MCL + exon conjunct | **0.7016** |
| A− — MCL alone | 0.6394 |
| **B — global Jaccard, connected components** | **0.6044** ← the thing being improved |
| C — Jaccard + conjunct, no MCL | 0.5869 |

## The bar — committed now

Let **F_D** be the chained arm's pooled held-out F with **t selected on chr16**, not on the held-out set
(the same discipline §6t3 applied to arm B).

| outcome | verdict |
|---|---|
| F_D ≥ 0.6394 (beats MCL-alone) | ⭐⭐ **MAJOR** — chaining replaces the operator, report for adoption |
| 0.6244 ≤ F_D < 0.6394 (beats B by ≥ 0.02) | ⭐ **IMPROVED** — chaining is the better scalar |
| 0.5844 ≤ F_D < 0.6244 | ⚠ **NEUTRAL** — not worth the extra code |
| F_D < 0.5844 (worse than B by > 0.02) | ⛔ **WORSE** — report it and keep global Jaccard |

⚠ Declared before the run: I will report **both** t-selected-on-chr16 and best-on-held-out, and treat the
former as the result. I will not re-tune γ, change the normaliser, drop a chromosome, or change the truth
after seeing scores. If chaining wins only at a t that does not transfer from chr16, that is a **NEUTRAL**,
not a win.
