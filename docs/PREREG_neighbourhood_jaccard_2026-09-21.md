# Pre-registration — a Jaccard-like metric on NODE CONNECTIVITY, with the size trap controlled

**Written 2026-09-21 before any neighbourhood Jaccard is computed.** User: *"figure out a jaccard like
metric for node definition or connectivity."*

This moves off the edge-scalar axis, which is exhausted (§6t3–§6u2: Jaccard, Ochiai, Dice, containment,
guarded containment, DP chaining, protein edges, junction sharing — every one at or below the shipped
rule), and onto §6o8's own priority list, where grouping was ranked last and node/connectivity work
above it.

## The metric

For an edge (u,v) of the shipped graph, with N(x) the set of graph neighbours of x:

    J_N(u,v) = |N(u) ∩ N(v)| / |N(u) ∪ N(v)|            (neighbours, excluding u and v themselves)

A real family should be clique-like — members share neighbours — while a bridge between two families
should be star-like.

## Two prior results define exactly what must be controlled

- **r459** — this metric, standalone **AUC 0.826**, the highest of any separator in the register. But it
  was refuted as a hard filter (*"cuts 3 bridges at 9 paralogs = base rate"*) and carries a structural
  warning: ⚠ ***"a 2-copy family is an isolated edge with no common neighbours — and 57% of families are
  pairs."*** For the modal family the metric is identically zero.
- **r523** — edge betweenness, **AUC 0.683 → 0.531 after size-residualisation**, with correlation
  **−0.71** against log component size. ⭐ **Topology statistics on this graph are largely component-size
  in disguise, and r459's 0.826 has not been size-residualised.** That is the first thing this run does.

Also on record: **r661** (triangle support, articulation points, k-core degree — *"conditioned on
`direct`, none raises F1"*) and **r663** (*"family min-degree ≥ 2"* gives pair-level precision 0.9047 but
**0.7348** when implemented as the catalog operation it implies). Neither is re-proposed; both say a
connectivity statistic must be evaluated as the operation it implies, not as a pair score.

## Method — frozen

Graph: the **shipped** DNA gene-body graph (`mcl_families --dump-graph`, conjunct applied), held-out
chr2/chr8/chr10. Labels: an edge is TRUE if both endpoints are in the same Soto family, FALSE otherwise;
only edges with both endpoints Soto-labelled are scored.

Reported, in this order:

1. **the pair fraction** — how many Soto truth families have exactly 2 members, i.e. the population where
   J_N is structurally blind;
2. **raw AUC** of J_N, comparable to r459's 0.826;
3. ⭐ **size-residualised AUC** — AUC computed WITHIN strata of the endpoints' component size
   (bins: 2, 3–4, 5–9, ≥10), then pooled by stratum weight. This is r523's control and it is the number
   that decides the run;
4. the correlation of J_N with log component size, comparable to r523's −0.71.

## The bar — committed now

| outcome | verdict |
|---|---|
| **size-residualised** AUC ≥ 0.75 | ⭐⭐ **A REAL CONNECTIVITY METRIC** — survives the control that killed betweenness |
| 0.65–0.75 | ⭐ **USABLE SIGNAL** — better than every edge scalar tested (best was multiplicity 0.681) |
| 0.55–0.65 | ⚠ **SIZE IN DISGUISE** — the r523 outcome |
| < 0.55 | ⛔ **NO** |

⚠ Declared now: a high **raw** AUC is not a result. If raw is high and residualised is not, the finding
is that **r459's 0.826 was component size**, and that is what will be reported. I will also report the
metric's coverage — the fraction of TRUE edges for which J_N > 0 — since a metric that is zero on the
modal family cannot define families however well it ranks the rest.
