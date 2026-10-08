# Neighbourhood Jaccard — a real connectivity metric, and exactly where it is blind

Run 2026-09-21 against `docs/PREREG_neighbourhood_jaccard_2026-09-21.md` (committed `8ab4e902` before
any value was computed). Tool: `bench/neighbourhood_jaccard.py`. Held-out chr2/chr8/chr10, Soto S1C
labels, edges scored only where both endpoints are Soto-labelled.

    J_N(u,v) = |N(u) ∩ N(v)| / |N(u) ∪ N(v)|      over graph neighbours, u and v themselves excluded

## It survives the control that killed every other topology statistic

| population | TRUE / FALSE | raw AUC | **size-residualised AUC** | corr with log comp size |
|---|---|---|---|---|
| **pre-registered**: pre-conjunct graph edges | 197 / 24 | 0.739 | **0.740** | +0.419 |
| post-hoc power check: all labelled aligned pairs | 238 / 206 | 0.861 | **0.845** | +0.489 |
| *r523's edge betweenness, for contrast* | — | *0.683* | ***0.531*** | *−0.71* |

⭐⭐ **The residualised AUC does not drop.** That is the whole point: r523 established that topology
statistics on this graph are largely component size in disguise, and betweenness lost 0.152 AUC under
stratification. Neighbourhood Jaccard loses **0.001** on the pre-registered population and **0.016** on
the larger one. **r459's standalone AUC 0.826 is confirmed (0.861 raw here) and it is not an artefact.**

At **0.845 residualised it is the strongest separator measured anywhere in this project** — against
graph multiplicity 0.681 (§6t8), exonic fraction 0.655 (§6u0), alignment identity ≈ 0.50 (§6t7), and a
flat containment curve.

## And r459's structural warning is confirmed exactly, to three decimals

| component size | TRUE | FALSE | AUC |
|---|---|---|---|
| **2** | 15 | 20 | **0.500** |
| 3–4 | 39 | 31 | 0.872 |
| 5–9 | 57 | 64 | 0.808 |
| ≥ 10 | 122 | 66 | **0.924** |

⛔ **On isolated pairs the metric is exactly chance.** r459 said it in words — *"a 2-copy family is an
isolated edge with no common neighbours"* — and the measurement lands on 0.500. Meanwhile **30 of 50
truth families on these chromosomes (60.0%) have exactly two members**, matching r459's cited 57%.

So the metric is excellent on large components (0.924 at ≥ 10) and worthless on the modal family. It
cannot *define* families; it can *score* edges inside components that are already large.

## What follows, and what does not

- ⭐ **This is the connectivity metric the goal asked for, and it is real**: parameter-free, one line,
  survives size-residualisation, best-in-project AUC.
- ⛔ **It is not a family definition.** 60% of families are pairs, where it is chance. Any rule built on
  it inherits that blindness, which is the same shape as §6n5's trap — every triangle-based operator
  scored well on F while dissolving all two-member groups.
- ⚠ **It was already refuted as a hard filter** (r459: *"cuts 3 bridges at 9 paralogs = base rate"*), and
  nothing here re-proposes that. What this run establishes is that the *signal* is genuine, not that a
  filter built from it works — r663's lesson stands: a connectivity statistic must be evaluated as the
  catalog operation it implies, not as a pair score.
- ⭐ The defensible use is **an edge-level confidence score reported alongside large-component edges**,
  where §6t3 showed precision is the axis that matters, and explicitly abstaining on components of size 2.

## Limits

- 3 chromosomes, 50 truth families, 24 FALSE edges in the pre-registered population — the post-hoc
  broadening to 206 FALSE is what gives the stratified table its power, and it is post-hoc.
- Neighbourhoods come from the pre-conjunct graph (identity ≥ 0.7, cov_longer ≥ 0.3, ≥ 300 bp); a
  different base graph changes N(u) and therefore the metric.
- Component size and family size are correlated with the truth here, so the stratification is a control,
  not a proof of independence.
