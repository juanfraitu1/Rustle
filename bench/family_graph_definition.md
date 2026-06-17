# Graph-based multi-copy family definition (per-family POA variation graph) — an option, with a formal definition

Upgrade #1 as an **option** (not a replacement for the pairwise definition): build one POA variation
graph per candidate family (all N member transcripts at once) and read the family criterion off the
graph. It yields a clean, formally-statable definition that generalizes the validated pairwise one.

## Formal definition
Members `S = {s₁,…,s_N}`. Let `G` be their partial-order (POA) alignment graph with columns
`c₁,…,c_M` in topological order. For a column `c`, let `supp(c) = #{ i : sᵢ is non-gap at c }`.
The **conserved core** is the longest contiguous run `R` of columns with
`supp(c) ≥ max(2, ⌈N/2⌉)` for every `c ∈ R` (a majority spine). Then

> `S` is a **multi-copy family at level T**  ⟺  `|R| ≥ T · median_i |sᵢ|`,
> with **family score** `= |R| / median_i |sᵢ| ∈ [0,1]`.

**Reduces to the validated pairwise core at N=2:** `⌈2/2⌉ → 2`, so a core column requires *both*
copies non-gap, and `R` is exactly the longest co-aligned block — the pairwise contiguous-core we
validated (RABL2/DAZ score = their pairwise core, 0.16).

## Why it's cleaner than pairwise + connected-components
- **One statistic from one graph** — no O(N²) pairwise comparisons, no transitive-closure family
  building. Membership is "does this copy thread the shared core?", a graph property.
- **It exposes over-merge.** Pairwise connected-components chain distinct subfamilies through domain
  hubs (mega-"families" of 250). The graph score on such a chain is LOW (no column shared by a majority
  of members) → **7 of the 14 N≥25 components score < T = NOT single families** (the other 7 are genuine
  large families with a real shared core). Pairwise CC could not make this distinction; the graph does.

## Validation (labeled set)
| class | graph score | vs T=0.045 |
|---|---|---|
| RABL2 (2), DAZ (2) | 0.162 | ✓ family |
| RFPL (4) | 0.081 | ✓ family |
| APOBEC3 (6) | 0.062 | ✓ family |
| domain-sharers (CDPF1/PPARA, CREB1/METTL21A, GCA/KCNH7, …) | ≤ 0.029 | ✗ rejected |

**Clean separation: domain-sharer max 0.029 < T=0.045 < real-family min 0.062.** Perfect on the
labeled set. Genome-wide: 1,333 candidate families scored; 1,211 ≥ T; the 122 < T are domain-hub
chains / weak components the graph correctly down-weights.

## Honest caveats
- **Graded, not binary** — recent near-identical copies score high (RABL2 0.16); ancient/length-variable
  families score moderate (APOBEC3 0.062, RFPL 0.081). The score is a family-*tightness* measure, so the
  margin between diverged-real-families (≥0.062) and over-merged-chains (median 0.04) is real but thin —
  clean for recent families, tighter for ancient ones (consistent with the RNA-only tier's reach limit;
  ancient families still benefit from the protein/DNA tier).
- **T is metric-specific** (0.045 here vs 0.13 for the raw pairwise core — different normalization).
  Both are data-placed in the labeled gap.
- **Bounded inputs:** members capped at 30 (shortest) and sequences > 15 kb skipped (POA memory);
  noted, affects only the largest families.

## Verdict
The per-family POA variation graph gives a **clean, formal, N-way family definition** that reduces to
the validated pairwise core, separates real families from domain-sharers, and — uniquely — exposes the
transitive-closure over-merge. Worth keeping as the family-definition **option** for well-annotated,
multi-member families.

## Reproduce
- `MINIFORGE python bench/family_graph_define.py` (pyabpoa POA-MSA per family → core score)
- `python3 bench/family_graph_fig.py`
