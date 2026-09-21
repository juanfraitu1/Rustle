# A pair metric that tolerates length asymmetry — bakeoff result

Run 2026-09-20 against `docs/PREREG_asymmetry_metric_2026-09-20.md` (committed `b788bbe5` before any
new metric was scored). Tool: `bench/pair_metric_sweep.py`. Same truth (Soto S1C ≥ 3 members), same
scorer, same held-out chromosomes (chr2, chr8, chr10), t always selected on chr16.

## The blind spot is real, and bigger than expected

Recall on **true** (same-Soto-family) pairs, split by length ratio, at each metric's chr16-selected t:

| metric | similar length (179 pairs) | **asymmetric, min/max ≤ 0.5 (59 pairs)** |
|---|---|---|
| Jaccard `m/(la+lb−m)` | 96.1% | **0 / 59 = 0.0%** |
| **Ochiai `m/√(la·lb)`** | 95.5% | **0 / 59 = 0.0%** |
| Dice `2m/(la+lb)` | 95.0% | **0 / 59 = 0.0%** |
| guarded containment (ratio ≥ 0.25) | 97.8% | **31 / 59 = 52.5%** |
| containment `m/min(la,lb)` | 97.8% | **54 / 59 = 91.5%** |

⭐ **59 of 238 true pairs — 25% — are length-asymmetric, and every symmetric metric finds none of them.**
That quantifies r359 far beyond its single EEF1A1 example: it is not an edge case, it is a quarter of
the true pair population, and the incumbent metric has exactly zero recall there.

## But the principled candidate failed

| arm | held-out F | sens | prec | verdict |
|---|---|---|---|---|
| A — shipped: MCL + exon conjunct | **0.7016** | 0.808 | 0.722 | — |
| A− — MCL alone | 0.6394 | 0.814 | 0.586 | — |
| B — Jaccard (incumbent) | 0.6044 | 0.713 | 0.639 | — |
| **E1 — Ochiai (the pre-registered candidate)** | **0.6044** | 0.713 | 0.639 | ⚠ NEUTRAL |
| E2 — Dice | 0.5952 | 0.703 | 0.632 | ⚠ NEUTRAL |
| **E3 — guarded containment** | **0.6177** | **0.809** | 0.562 | ⚠ NEUTRAL (best of the four) |
| E4 — union(Jaccard, guarded) — POST-HOC | 0.6177 | 0.809 | 0.562 | collapses to E3 |

⛔ **Ochiai reproduces Jaccard exactly** (0.6044 / 0.713 / 0.639) and scores **0%** on the very
population it was chosen for. The geometric mean is better *shaped* — it cannot saturate the way
containment does — but that does not matter, because **a single global threshold, selected on a
population dominated by similar-length pairs, sits above anything an asymmetric pair can reach.**
E4 collapses to E3 for the same structural reason: containment ≥ Jaccard always, so `max` is the guard.

## The real conclusion: it is not a metric-shape problem, it is a one-threshold problem

No scalar tested clears the pre-registered IMPROVED bar (0.6244). The best, guarded containment, buys
**+0.0133 F**, converts **sensitivity 0.713 → 0.809** and **52.5% of the asymmetric blind spot**, and
pays **precision 0.639 → 0.562**. That is a real trade, not noise, and it points at the structure:

⭐ **Two populations need two criteria.** Similar-length pairs are served well by any symmetric metric
(~96%); asymmetric pairs are reachable only by containment, which alone builds hubs (r913: precision
0.164). Guarded containment is a first cut that recovers half the blind spot without the hubs — and
where the lost precision would have to come back from is **exactly what §6t3 measured**: the exon
conjunct inside MCL, which buys +0.136 precision and only survives inside MCL (r911).

**The next test that follows from this** — untested here — is guarded containment as the edge metric
*feeding MCL and the conjunct*, rather than feeding connected components. This bakeoff varied the
scalar while holding the operator at its weakest setting.

## ⭐ And it DOES help — but only in the right operator (POST-HOC)

The pre-registered bakeoff held the operator at connected components, its weakest setting. Feeding the
same scalars to **MCL (I = 2.8, the shipped default)** instead, t still selected on chr16:

| arm | held-out F | sens | prec |
|---|---|---|---|
| Jaccard → components (B, incumbent) | 0.6044 | 0.713 | 0.639 |
| guarded containment → components (E3) | 0.6177 | 0.809 | 0.562 |
| Jaccard → **MCL** | 0.6134 | 0.730 | 0.623 |
| ⭐ **guarded containment → MCL** | **0.6648** | 0.804 | **0.650** |
| A− — MCL alone, **shipped** identity/cov metric | 0.6394 | 0.814 | 0.586 |
| A — shipped: MCL + exon conjunct | 0.7016 | 0.808 | 0.722 |

⭐ **Guarded containment inside MCL beats the shipped pipeline's own metric at the same operator
setting: 0.6648 vs 0.6394, +0.0254** — and it does so while *raising* precision (0.586 → 0.650), not
trading it away. The same metric in connected components was NEUTRAL (+0.013). **The metric is worth
+0.047 more inside MCL than outside it.**

⚠ This is **post-hoc** — the pre-registration fixed the operator at connected components, so this arm
was not committed in advance and is reported as exploratory. The adoption test is different again:
guarded containment inside the Rust `mcl_families`, *with* the exon conjunct, which is the remaining
0.037 to the shipped number.

⭐ It is the same lesson as r911 from the other direction: **a scalar cannot be judged at the operator's
weakest setting.** Jaccard gains only +0.009 from MCL; guarded containment gains +0.047, because MCL's
flow can exploit the extra asymmetric edges that connected components turn into hubs.

## Limits

- 238 true pairs over three chromosomes; the asymmetric stratum is 59 of them.
- `min/max ≤ 0.5` is one cut; the blind spot's size will move with it, its existence will not.
- E3's 0.25 guard is the one tuned parameter in the bakeoff and was fixed before scoring, not swept.
