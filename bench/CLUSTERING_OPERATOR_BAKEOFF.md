# Clustering-operator bakeoff on the L2 copy graph — DESCRIPTIVE, development families only

⚠⚠ **NPIP and TBC1D3 are DEVELOPMENT families. Picking a winner from this table is the dev-set selection
trap that has cost retractions before. Nothing here is adopted; a candidate must be confirmed on a held-out
substrate first.** Recorded because several operators beat the shipped one by large margins and that is
worth knowing.

Graph: the §0★★★ **L2** level (`records.tsv`, 369 nodes, 1,106 edges, cut 0.30). Same graph for every
operator — only the grouping rule changes.

## Quality

| operator | NPIP F | sens | prec | parts | TBC1D3 F | sens | prec | parts |
|---|---|---|---|---|---|---|---|---|
| connected components (baseline) | 0.287 | 0.926 | 0.081 | 2 | 0.471 | 0.386 | 0.286 | 8 |
| 3-truss components | 0.308 | 0.926 | 0.093 | 2 | 0.511 | 0.386 | 0.347 | 8 |
| 4-truss components | 0.335 | 0.926 | 0.111 | 2 | 0.686 | 0.386 | 0.725 | 8 |
| 3-clique percolation | 0.515 | 0.926 | 0.265 | 2 | 0.686 | 0.386 | 0.725 | 8 |
| **4-clique percolation** | **0.536** | 0.926 | 0.288 | 2 | 0.686 | 0.386 | 0.725 | 8 |
| **louvain (modularity)** | 0.515 | 0.855 | **0.472** | 3 | 0.471 | 0.386 | 0.286 | 8 |
| greedy modularity | 0.515 | 0.926 | 0.265 | 2 | 0.471 | 0.386 | 0.286 | 8 |
| **label propagation** | 0.515 | 0.926 | 0.265 | 2 | **0.686–0.727** | 0.386 | 0.725 | 8 |

## Coverage — the counterweight the F column hides

| operator | groups >= 2 | nodes covered | % of component nodes | **2-copy components kept** |
|---|---|---|---|---|
| connected components | 43 | 369 | 100.0% | **13/13** |
| 3-truss components | 21 | 190 | 51.5% | **0/13** |
| **4-truss components** | 7 | 119 | **32.2%** | **0/13** |
| 3-clique percolation | 27 | 190 | 51.5% | **0/13** |
| **4-clique percolation** | 9 | 119 | **32.2%** | **0/13** |
| **louvain** | 47 | 369 | **100.0%** | **13/13** |
| **label propagation** | 59 | 369 | **100.0%** | **13/13** |

**Every triangle-based operator dissolves all 13 two-copy components**, and 2 is the modal family size
(§1★.5). 4-clique percolation's best-in-table F 0.536 is bought by discarding 68% of the nodes — the same
cost that ruled out the truss options in §0★★★.

## Determinism (5 seeds)

| operator | NPIP F | TBC1D3 F | groups |
|---|---|---|---|
| **louvain** | 0.515–0.515 (sd **0.000**) | 0.471–0.471 (sd **0.000**) | 47 every time |
| label propagation | 0.515–0.515 (sd 0.000) | 0.686–**0.727** (sd 0.019) | 57–59 |

**Label propagation is not deterministic** — it moves on TBC1D3 across seeds. A definition cannot be
seed-dependent, so it is out regardless of its score.

## Reading

**Louvain is the only operator that improves on the baseline while keeping full coverage and determinism:**
NPIP F 0.287 → 0.515 with precision 0.081 → **0.472 (5.8x)**, all 13 two-copy components kept, identical
across seeds.

Its costs are real and must be stated with it:
- NPIP **sensitivity drops 0.926 → 0.855** and the family splits into 3 parts, not 2.
- **TBC1D3 is unchanged** (0.471) — the gain is NPIP-only in this table.
- ⚠ **Modularity is not a nested/laminar family.** It would break the T1/T1′ nesting theorems that are
  §0★★★'s actual theoretical contribution, and it carries a resolution limit. For an advisor who wants
  clean combinatorial structure and provable theorems, that is a large price for +0.228 F on one family.

**Recommendation: do not adopt anything from this table yet.** The one candidate worth a pre-registered
held-out test is Louvain, and the test must report the nesting loss as a cost, not omit it.
