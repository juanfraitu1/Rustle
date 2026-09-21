# Cliques, trusses, communities — the operator bakeoff on held-out substrate

Run 2026-09-20 against `docs/PREREG_operator_heldout_2026-09-20.md` (committed `054030c2` before any
operator was scored). Tool: `bench/operator_bakeoff_heldout.py`. Graph: the **shipped** DNA gene-body
graph via `mcl_families --dump-graph` (edge set and exon conjunct exactly as shipped). Substrate:
chr2/chr8/chr10, zero ledger and register exposure. Truth: Soto S1C, same scorer as every other arm.

This is the confirmation `bench/CLUSTERING_OPERATOR_BAKEOFF.md` (§6n5) said was required and could not
run: it tested these operators on NPIP and TBC1D3 and warned that picking a winner there is the dev-set
selection trap.

## Result — nothing beats MCL, and the dev-set ranking INVERTS

| operator | F | sens | prec | node coverage | 2-member groups | verdict |
|---|---|---|---|---|---|---|
| **label propagation** | **0.7295** | 0.872 | **0.740** | **100.0%** | **370** | ⚠ TIED |
| **MCL I=2.8 (shipped)** | 0.7123 | 0.838 | 0.715 | 95.4% | 358 | ⚠ TIED |
| louvain | 0.6535 | 0.881 | 0.594 | 100.0% | 362 | ⛔ WORSE |
| greedy modularity | 0.6318 | 0.831 | 0.580 | 100.0% | 362 | ⛔ WORSE |
| 3-clique percolation | 0.6138 | 0.756 | 0.604 | 37.9% | **0** | ⛔ WORSE |
| connected components | 0.5811 | 0.781 | 0.536 | 100.0% | 362 | ⛔ WORSE |
| 3-truss | 0.5776 | 0.741 | 0.544 | 37.6% | **0** | ⛔ WORSE |
| 4-clique percolation | 0.4731 | 0.617 | 0.420 | 25.5% | **0** | ⛔ WORSE |
| 4-truss | 0.4065 | 0.575 | 0.354 | 24.8% | **0** | ⛔ WORSE |

⭐⭐ **The dev-set winner becomes the loser.** On NPIP, §6n5 measured 4-clique percolation **best in
table** (F 0.536 against connected components' 0.287). On held-out substrate it is **0.4731 — second
worst of nine**. Louvain, the other dev-set standout (precision 0.081 → 0.472 on NPIP), lands at
0.6535, below the shipped operator. **§6n5's refusal to adopt from that table was correct, and this is
the measurement that shows it.**

⭐ **§6n5's coverage warning reproduces exactly.** Every triangle-based operator — 3-truss, 4-truss,
3-clique, 4-clique — **dissolves all 362 two-member groups** and covers 25–38% of nodes. Two is the
modal family size, so these operators do not define families so much as discard most of them.

## The one arm worth a follow-up

**Label propagation is better than the shipped MCL on every axis reported**: F +0.0172, sensitivity
+0.034, precision +0.025, coverage 100% vs 95.4%, and 370 two-member groups vs 358. It is also
**deterministic** — 1 distinct partition over 5 repeats on each of the three chromosomes (the
*asynchronous* `asyn_lpa` variant is seed-dependent, 3 partitions over 5 seeds, and is not what was
used).

⚠ **But the pre-registered bar for BETTER was +0.02, and this is +0.0172, so it is a TIE, not a win** —
and after a day in which arms of exactly this size evaporated (§6t6's +0.0133 was two families), that
distinction is worth keeping. It is a candidate for a pre-registered confirmation on more chromosomes,
not a reason to change the operator.

## What this settles

⛔ **No graph-theoretic operator tested is better at defining a family on held-out data.** Cliques and
trusses are much worse and destroy the modal family size; modularity communities (Louvain, greedy) are
worse; connected components are worse. This is the fourth independent confirmation of §6o8's *"grouping
rules — already saturated"*, now on a substrate chosen for having no exposure to the tuning history.

⭐ The operators spread over F 0.41–0.73 while every scalar tested in §6t3–§6t5 spread over 0.58–0.73,
so the operator is not *irrelevant* — a bad one costs a lot. It is that **the shipped one is already at
the top of the range**, and the headroom above it is a rounding error.

## Limits

- 27 Soto families over three chromosomes; F differences under ~0.02 are not resolvable here.
- The graph is the DNA gene-body graph, not §6n5's RNA L2 graph, so this is the same question on a
  different substrate, not a replication of their numbers.
- MCL is scored through `bench/mcl_port.py`, which is not bit-identical to the Rust MCL; all operators
  in the table go through the same Python path, so the comparison is internally consistent.
