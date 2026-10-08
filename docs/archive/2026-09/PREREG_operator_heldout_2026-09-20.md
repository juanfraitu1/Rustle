# Pre-registration — the clustering-operator bakeoff, on HELD-OUT substrate

**Written 2026-09-20 before any operator is scored on held-out data.** User: *"determine if other graph
theory concepts such as cliques, dense components, communities etc would be better at defining a spliced
family."*

## This is the experiment §6n5 said was required and could not run

`bench/CLUSTERING_OPERATOR_BAKEOFF.md` already tested exactly these operators — connected components,
3/4-truss, 3/4-clique percolation, Louvain, greedy modularity, label propagation — and opens with its own
warning:

> ⚠⚠ *"NPIP and TBC1D3 are DEVELOPMENT families. Picking a winner from this table is the dev-set
> selection trap that has cost retractions before. Nothing here is adopted; a candidate must be confirmed
> on a held-out substrate first."*

On development families several operators beat components by large margins (NPIP F: components 0.287 →
4-clique percolation **0.536**; Louvain precision 0.081 → **0.472**). **This runs the same question on
chr2/chr8/chr10 — zero ledger and register exposure — with Soto's published families as truth**, which
§6n5 did not have.

⚠ **Not a replication.** §6n5's graph was the RNA **L2** copy graph (369 nodes, 1,106 edges). Mine is the
**DNA gene-body** graph that `mcl_families` actually ships, obtained with `--dump-graph` so the edge set
and the exon conjunct are exactly the shipped ones. Same operator question, different substrate, and the
substrate is the one with held-out truth.

## The operators — every one on the identical graph

connected components (baseline) · **3-truss** · **4-truss** · **3-clique percolation** ·
**4-clique percolation** · **Louvain** · **greedy modularity** · **label propagation** ·
**MCL I=2.8** (the shipped operator, via `bench/mcl_port.py`).

## What must be reported together

§6n5's central lesson is that the F column hides the cost: **every triangle-based operator dissolved all
13 two-copy components**, and 2 is the modal family size. So each operator reports, on the held-out set:

1. pooled bipartite **F / sensitivity / precision** (`bench/heldout_family_score.py`, Soto truth);
2. **node coverage** — how many graph nodes end up in any group of ≥ 2;
3. **two-member groups retained**, the statistic that killed the truss and clique options before.

## The bar — committed now

Shipped MCL through the same port scores **F 0.7123** on this substrate (§6t6). An operator is:

| outcome | verdict |
|---|---|
| F ≥ 0.7123 + 0.02 **and** node coverage ≥ 90% of MCL's | ⭐⭐ **BETTER** — a genuine candidate, report for adoption |
| F ≥ 0.7123 + 0.02 but coverage < 90% | ⚠ **BOUGHT WITH COVERAGE** — §6n5's trap, report as such, do not adopt |
| within ±0.02 of 0.7123 | ⚠ **TIED** — consistent with grouping being saturated (§6o8) |
| F < 0.7123 − 0.02 | ⛔ **WORSE** |

⚠ Declared now: §6o8 measured that **grouping is saturated** and its priority list says *"do not spend
further effort here"*. The expected outcome is TIED, and a TIED result is a real finding, not a null —
it would mean the operator genuinely does not matter on held-out data. I will not pick a winner on the
best-of column, and I will not drop the coverage columns if an operator wins on F.
