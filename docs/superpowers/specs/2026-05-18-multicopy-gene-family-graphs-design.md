# Multi-Copy Gene Family Graphs — Design Spec
Date: 2026-05-18

## Goal

Produce a strict formal definition of a multi-copy gene family grounded in variation graph theory, validated empirically on annotated gene families (AMY, GOLGA, NPIP) from the GGO gorilla genome, and communicated via three publication-quality figures for advisor review.

## Scientific Question

What property of the variation graph distinguishes a **copy** (paralog locus) from an **isoform** (splice variant)? Candidates: genomic locus distinctness, copy-specific bubble depth, or read-linkage topology. The data will decide.

---

## Formal Definition (to be refined by data analysis)

### Objects

**Exon equivalence class** `v`: the union of overlapping exon intervals contributed by any annotated transcript across all loci in the family. Each `v` carries:
- genomic span `(chrom, start, end, strand)`
- copy membership set `C(v) ⊆ {1…n}` — which loci contribute to this node

**Variation graph** `G = (V, E, λ)`:
- `V` — set of exon equivalence classes
- `E ⊆ V × V` — directed junction edges (consecutive exons in any transcript)
- `λ: V → {shared, copy-specific}` — shared if `|C(v)| ≥ 2`, copy-specific if `|C(v)| = 1`

**Copy path** `P_i`: the ordered sequence of nodes in `V` traversed by locus `i` (one path per annotated gene locus).

**Isoform** of locus `i`: any path through `G` that is a sub-path or skip-variant of `P_i` — it may omit nodes of `P_i` via skip edges but does not traverse copy-specific nodes belonging to a different locus `j ≠ i`.

### Definition: Multi-Copy Gene Family

A **multi-copy gene family** is a tuple `F = (G, {P_1, …, P_n})` where:
1. `n ≥ 2` (at least two distinct loci)
2. Each `P_i` traverses at least one copy-specific node `v` with `C(v) = {i}` (loci are not identical)
3. At least one node `v ∈ V` is shared: `|C(v)| ≥ 2` (loci share conserved exon structure)
4. The graph `G` is connected (all loci reachable from a common source node)

**Operational criterion**: the family is detected when long reads multi-map across ≥2 loci — this is what triggers family construction; conditions 1–4 are what the resulting graph must satisfy.

### What this excludes (counterexamples)

| Case | Why excluded |
|---|---|
| Single-locus gene with multiple isoforms | `n = 1`, condition 1 fails |
| Two unrelated genes with no shared exons | condition 3 fails (no shared node) |
| Two identical copies (no divergence) | condition 2 fails (no copy-specific node) — degenerate; treat as one locus |
| Read-through transcripts linking two genes | condition 4 may hold but condition 3 will fail if the shared "node" is an artifact of run-on transcription, not homology |

---

## Architecture

```
GGO_genomic.gff
      │
      ▼
analysis/family_graphs/
  00_install.sh           mamba install all R packages
  01_parse_gff.R          rtracklayer → exon tables per family (RDS)
  02_build_graphs.R       igraph DAGs + copy paths + node labels (RDS)
  03_figure_problem.R     Figure 1 — multi-mapping ambiguity (3 panes)
  04_figure_method.R      Figure 2 — variation graph + copy paths (3 panes)
  05_figure_examples.R    Figure 3 — examples & counterexamples (3 panes)
  data/                   intermediate RDS
  figures/                output PDF + PNG
```

---

## Packages (all via mamba)

| Package | Channel | Role |
|---|---|---|
| `bioconductor-rtracklayer` | bioconda | parse GFF3 |
| `bioconductor-genomicfeatures` | bioconda | transcript/exon DB |
| `r-igraph` | conda-forge | graph construction & analysis |
| `r-ggraph` | conda-forge | graph layout & visualization |
| `r-ggplot2` | conda-forge | base plotting |
| `r-dplyr` + `r-tidyr` | conda-forge | data wrangling |
| `r-patchwork` | conda-forge | multi-panel assembly |

---

## Gene Families

| Family | GFF description pattern | Why |
|---|---|---|
| AMY | `amylase` | Small (5 members), well-characterized, salivary vs pancreatic split tests copy vs isoform distinction |
| GOLGA | `golgin\|GOLGA` | Primate-specific tandem duplications, strong copy-specific bubbles expected |
| NPIP | `nuclear pore.*interact\|NPIP` | Segmental duplications, tests read-linkage detection |

**Note**: gorilla GFF uses `LOC*` gene IDs. Families must be filtered by the `description` or `product` attribute in column 9, not by gene name prefix. Filter with case-insensitive regex on description field.

Start with AMY for validation, then scale.

---

## Figures

### Figure 1 — The Problem (3 panes)
| Pane | Content |
|---|---|
| A | Genomic strip: 3 AMY loci as colored exon blocks on chromosome line, scale bar |
| B | Same strip + read arcs: each arc spans the two loci a multi-mapping read hits, colored by read |
| C | Naive tool output: all arcs collapsed to locus 1, others empty. "reads mis-assigned or discarded" |

### Figure 2 — The Method (3 panes)
| Pane | Content |
|---|---|
| A | AMY variation graph DAG: gray nodes = shared exons, colored nodes = copy-specific, node size ∝ exon length |
| B | Same DAG + copy-path ribbons: one color per annotated locus |
| C | Zoomed subgraph: one shared node branching into two copy-specific bubbles; dashed edge = isoform skip. Anchors the formal definition. |

### Figure 3 — Examples & Counterexamples (3 panes)
| Pane | Content |
|---|---|
| A | **Example** GOLGA: DAG with copy-specific bubbles, 2+ loci. "multi-copy family ✓" |
| B | **Counterexample 1**: single-locus gene, skip edges only, no copy-specific nodes. "isoforms only ✗" |
| C | **Counterexample 2**: two disconnected DAGs, no shared nodes. "no shared exon structure ✗" |

---

## Validation Criteria

- AMY1A/B/C paths should share nearly all nodes (copy-specific bubbles only at SNP-level divergence)
- AMY1 vs AMY2 should branch at an early shared node with structurally distinct downstream subgraphs
- GOLGA copies should show prominent copy-specific bubble chains
- Single-gene counterexample must have zero copy-specific nodes by construction

---

## Out of Scope (Phase 1)

- Running rustle on GGO.bam (Phase 2)
- HMM-EM read assignment (Phase 2)
- Novel copy discovery (Phase 2)
- GFA/VG export

---

## Next Step

Implementation plan via writing-plans skill.
