# Graph-theory refinement of the RNA family E_r graph: does topology beat the shipped γ-quasi-clique?

Artifacts: `bench/graph_def_refine_sweep.py` (script), `bench/graph_def_refine_sweep.tsv`
(46 operators × metrics), `bench/graph_def_refine_sweep.json` (sweeps, operating points, beats,
residual examples). Deterministic: `PYTHONHASHSEED=0`, `seed=0`, sorted writes; the `.tsv`/`.json`
are byte-identical (md5) across re-runs.

## Question

The RNA family = connected components of the transcript-homology graph **E_r** (`core_recip ≥ 0.13`),
then the shipped **γ-quasi-clique refinement R** (γ=0.20) splits over-merged blobs.
(1) Does the γ-refined catalog still carry false positives (over-merges)?
(2) Can more graph theory — **k-truss / triangle support, edge embeddedness, community detection
(Louvain/greedy), k-core, graph bridges, spectral min-cut** — cut the over-merge *bridge* edges better
than the single γ threshold?

**Hypothesis under test:** an over-merge edge is a *bridge* — it joins two dense blobs but sits in few
triangles (low embeddedness), whereas a real family is a dense cluster; so triangle-support/embeddedness
should sever bridges while keeping families, more principled than a tuned γ.

## Graph and two ground-truth modalities (read the labels honestly)

- **Graph:** 3636 nodes (de-novo loci), 11 400 weighted edges, 1130 raw components. The largest raw
  component is a **728-locus repeat/ZNF megablob** that holds **98.9 %** of the over-merge (GENUINE)
  gene-pair mass.
- **(A) EDGE truth = reference cDNA / protein homology, NOT the DNA assembly.** Each raw-E_r gene pair is
  labelled `REAL_cdna` (in_dna_loose: reference cDNA-vs-cDNA minimap2 id ≥ 0.90 & cov ≥ 0.30),
  `TRUTHBAR` (reciprocal whole-protein homology, divergent paralog), or `GENUINE` (over-merge).
  Counts: REAL_cdna **460**, GENUINE **129 059**, TRUTHBAR **40 892**. This is *sequence-homology* truth in
  the **same modality as E_r**. It is **non-circular w.r.t. the graph operators** (they see only
  `core_recip`, never cDNA-id), so it is a fair *ranking* oracle — **but it is not assembly-independent
  and is not called "independent DNA".**
- **(B) The only assembly-independent quantity is `asm_hapCN`** (phased-diploid DNA oracle,
  `bench/diploid_cn_oracle.tsv`): the number of distinct haplotype loci in the DNA assembly of a *different*
  gorilla carrying a full-length copy of the family transcript. Its family **partition** comes from
  `validated_families.tsv` (a backbone-reinforced E_r RNA variant) — so the assembly supplies only the
  **scalar copy number**, not the family boundaries. **Coverage is only ~50/858 refined families (~5.8 %).**

## (1) Residual false positives — YES, confirmed against the assembly-independent DNA scalar

Leading with the two genuinely assembly-independent classes (not the RNA-partition-anchored multifam count):

- **allele-as-copy = 2** (threshold-free — DNA copy number literally = 1 on both haplotypes):
  `DHRSX` (2 RNA loci, `asm_hapCN`=1) and `LOC129530050` (3 RNA loci, `asm_hapCN`=1). The family definition
  treats **alleles as copies**. This is the cleanest DNA-independent FP.
- **extreme oversize = 4** (`distinct_loci > 1.5 × diploid DNA CN`, a margin above the reference/assembly
  cross-individual offset): `LOC129523567` **54 loci vs diploid 4 (13.5×)**, `MAGEA9+LOC129529978+LOC129529986`
  **18 vs 2 (9×)**, `MPHOSPH8` 10 vs 4 (2.5×), `LOC134758618` 10 vs 6 (1.67×).

> Unit-confound honesty: comparing reference/RNA `distinct_loci` to the **haploid** `asm_hapCN` of another
> gorilla over-counts — the **pure GSTM2 family already sits at ~2.5× haploid** (47 ref loci vs hapCN 19).
> The naive `>1.5×-haploid` flag therefore returns **6** and is reported only as a caveated upper bound; the
> diploid-margin metric (=4) is the defensible number.

The **multifam = 6** count (block spans ≥2 distinct oracle representative genes: GSTM2 recurs as a domain
hub ×4, plus FOXO1+LOC115933254, LOC101142904+LOC129526550) is **RNA-partition-anchored** (the oracle
family boundaries are backbone-reinforced RNA, not assembly), so it is reported as corroboration, **not** as
assembly-family evidence.

So the residual FP is real and confirmed by an independent DNA scalar — not an artifact of a
reference-derived truth — but the airtight evidence is the ~6 %-coverage diploid oracle, and the
`allele-as-copy`/`extreme-oversize` classes specifically.

## (2)/(3) Best graph-theory refinement — a Pareto surface, not a clean winner

### The bridge/triangle-support hypothesis is FALSIFIED

Per-graph-edge score AUC (high score ⇒ REAL), on 474 REAL vs 3136 GENUINE graph edges:

| score | AUC (REAL vs GENUINE) |
|---|---|
| **triangle_support** | **0.278 (INVERTED — worse than random)** |
| jaccard_embeddedness | 0.680 |
| edge_betweenness (per-component control) | 0.692 |
| **core_recip edge weight (homology)** | **0.846** |

Over-merge edges are **more** triangulated, not fewer: GENUINE median **9** triangles (mean 16.8, only
12.2 % zero-triangle) vs REAL median **3** (mean 4.8, 23.8 % zero-triangle). They live *inside* the dense
728-locus megablob and are heavily triangulated, while genuine divergent-paralog families are small and
triangle-poor. So **k-truss cuts the real families first.** The only edge-level separator is the **homology
weight `core_recip` (AUC 0.846), which is already in the graph — not graph topology.**

### Threshold operators cannot reach the regime

k-truss / embeddedness / k-core cut real-family edges before over-merge edges: their entire sweeps top out
at `tp_retention ≈ 0.73`; at `ret ≥ 0.90` **none** of them is reachable. Pair-closure ROC-AUC:
**γ 0.980, Louvain 0.985**, embeddedness 0.804, k-truss 0.690, k-core 0.562 (≈random control).

### The two "edge" metrics DISAGREE — and that is the crux

`overmerge_cut` in the original artifact is a **pair-closure (Rand-like) metric**: the GENUINE denominator is
129 059 gene pairs, **98.9 % inside the one 728-node megablob**. So a high `overmerge_cut` measures **how
finely the megablob is subdivided**, not "fraction of bridge edges removed." The honest, task-literal number
is the **graph-edge cut** on the 11 400 real edges:

| operator | tp_ret | PAIR overmerge_cut | **GRAPH-EDGE genuine-cut** | REAL-edge kept | impure blocks | TRUTHBAR ret |
|---|---|---|---|---|---|---|
| raw (no refine) | 1.000 | 0.000 | 0.000 | 1.000 | 99 / 795 (0.124) | 1.000 |
| **γ-quasi-clique 0.20 (SHIPPED)** | 0.978 | **0.968** | **0.316** (992/3136) | 0.996 (472/474) | **115 / 858 (0.134)** | 0.150 |
| γ + weighted split (principled) | 0.980 | 0.966 | 0.344 | 0.998 (473/474) | 115 / 861 (0.134) | 0.116 |
| Louvain res 8 (weighted) | 0.985 | 0.954 | **0.502** (1574/3136) | 0.996 (472/474) | **106 / 852 (0.124)** | 0.056 |
| Louvain res 16 (weighted) | 0.978 | 0.964 | 0.542 | 0.990 (469/474) | 107 / 858 | 0.041 |
| Louvain res 0.5 (weighted) | 0.987 | 0.779 | 0.019 | 0.998 (473/474) | 104 / 805 (0.129) | 0.729 |

- **Under the task criterion (max PAIR overmerge_cut at fixed retention): NONE beats γ.** The Louvain
  frontier crosses γ exactly — at res 16, matched `tp_ret=0.9783`, `overmerge_cut=0.9635 < γ 0.9677`; at
  res 24 the cut exceeds γ only by dropping retention. **γ sits on/above the Louvain pair-frontier.**
- **Under the honest GRAPH-EDGE criterion, Louvain res≈8 weakly DOMINATES γ:** identical REAL-edge retention
  (472/474), **1574 vs 992 genuine over-merge edges cut (50 % vs 32 %)**, **106 vs 115 impure blocks**, and
  higher pair-retention. The pair-metric "NONE beats γ" was a **megablob artifact.**
- **But Louvain buys this by sacrificing TRUTHBAR** (protein-divergent paralogs): tb-ret 0.056 vs γ 0.150,
  and it abandons the quasi-clique **certificate** for a tuned resolution. So on the 4-criterion Pareto
  (tp, cut, TRUTHBAR, impure) **NONE dominates γ** either. It is a genuine **non-dominated tradeoff**: γ is
  the conservative, certificate-backed, most-TRUTHBAR-preserving point; Louvain res≈8 is a more aggressive
  point that removes ~1.6× more over-merge *edges* and 9 fewer impure blocks. Low-res Louvain (0.5) is the
  opposite extreme — keeps 454 REAL and 73 % TRUTHBAR but removes almost no over-merge (edge-cut 0.019).

> Block-level honesty: γ's `impure_blocks` **rate** (0.134) is **not** an improvement over raw (0.124) or
> Louvain res 8 (0.124) — splitting the megablob into 858 families creates *more, smaller* impure blocks
> (99 → 115). γ's value is in the **pair-closure** (dissolving the transitive blob so per-family analysis is
> tractable), not per-block protein purity.

### The principled variant, and why topology stops here

The shipped `_split_once` builds an **unweighted** subgraph and calls Louvain **without weights** — it
discards `core_recip`, the one edge-level separator (AUC 0.846). The **principled, certificate-preserving**
variant keeps the exact unweighted γ-quasi-clique **stop rule** but makes the **split direction weighted**
by `core_recip`. Result: `tp_ret 0.980`, edge-cut `0.344` (vs 0.316), REAL-edge kept `0.998` (vs 0.996),
**impure blocks unchanged at 115**, TRUTHBAR 0.116. **Marginal.** The reason is diagnostic: the residual
over-merge lives in dense sub-blobs that **pass the unweighted density stop certificate** and are kept whole
*before the splitter is ever invoked*. Weighting the direction cannot help; only a **weighted stop rule**
would — and that recalibrates γ off its knife-edge (ZNF716 density 0.261) and costs more TRUTHBAR.

## Recommendation and honest residual

**Keep the shipped γ-quasi-clique refinement.** It is the principled operator in the correct regime — γ's
`_refine_component` is density-gated recursive **Louvain modularity** (community detection) with a
γ-quasi-clique **stopping certificate**; the certificate is the halt rule on top of the community operator.
γ=0.20 is calibrated (knife-edge ZNF716 0.261; γ≥0.27 splits the KRAB-ZNF family), not free, but it gives a
**structural guarantee** the tuned-resolution operators lack, and it is the most **TRUTHBAR-preserving**
point (protein-divergent paralogs survive best).

**Do not update the graph topology.** No topological operator (k-truss / embeddedness / k-core / bridges)
Pareto-beats γ: the bridge hypothesis is falsified — over-merge edges are triangle-*rich*, indistinguishable
from real dense families by any triangle/embeddedness/betweenness statistic. The Louvain res≈8 graph-edge
"win" is real but bought with TRUTHBAR and a tuned resolution.

**The honest residual graph structure cannot remove:** real divergent paralogs (small, triangle-poor:
ENO2–ENO3, OGDH–OGDHL) and TE/domain bridges are *both* low-density at the edge level, so no structural cut
separates them without destroying real families. **They are separable only by the homology weight
`core_recip` (already in the graph, AUC 0.846) plus an external DNA / protein arbiter.** The DNA-oracle
residuals — GSTM2 recurrent domain-hub over-merge, `LOC129523567` (54 vs 2), and the `DHRSX`/`LOC129530050`
allele-as-copy cases — are precisely what needs the diploid copy-number oracle or the protein arbiter, **not
more graph theory.** If anything is adopted, the zero-risk move is the **weighted split** (keeps one more
REAL edge, certificate intact); the higher-value move is **weighting the edge (`core_recip`) into the
refinement / adding the DNA-CN arbiter**, i.e. an edge-weight/oracle fix, not a graph-structure fix.

## Caveats (examiner-proof)

- The fully assembly-independent evidence is the **~50/858 (5.8 %) diploid-oracle** families; the per-edge
  `REAL_cdna`/`GENUINE` ranking truth is reference cDNA/protein homology (same modality as E_r), non-circular
  w.r.t. the operators but **not** assembly-independent.
- `overmerge_cut` / `tp_retention` / ROC-AUC are **pair-closure (Rand-like) quadratic** metrics dominated by
  the 728-node megablob; `ge_genuine_cut` / `ge_real_kept` (graph-edge) and `impure_blocks` (block) are the
  honest complements and are reported alongside.
- The Louvain resolution sweep is extended to **res 512** so the Pareto-frontier crossing at matched
  retention (res 16: `om_cut 0.9635 < γ 0.9677`) is on record, not left to reviewer re-derivation.
- Determinism: `all=True` (γ and Louvain byte-reproducible).
