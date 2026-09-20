# Multi-Copy Gene Family: Formal Definition

This note formalises the framework for multi-copy gene family analysis used by
the Rustle / family-graphs pipeline. All concepts are validated empirically
across ten real gene loci in the gorilla (GGO) annotation.

---

## 1. The variation graph

Let `L = {c_1, ..., c_N}` be a set of `N ≥ 2` annotated gene copies at clustered
genomic positions (typically the same chromosome). For each copy `c`, let
`Tx(c)` be the set of annotated transcripts of `c`, and for each transcript
`τ ∈ Tx(c)` let `E(τ)` be its ordered set of exons.

The **variation graph** `G(L) = (V, E)` is constructed by:

1. Pooling all exons across `L`: `U = ⋃_{c ∈ L} ⋃_{τ ∈ Tx(c)} E(τ)`.
2. Clustering `U` by sequence similarity (default: minimap2 all-vs-all
   with thresholds `min_identity = 0.85`, `min_coverage = 0.80`,
   cross-gene matches only). Each cluster becomes a single node `v ∈ V`.
3. For every transcript `τ`, the consecutive exons of `τ` become directed
   edges `(v_i, v_{i+1}) ∈ E` between their cluster representatives.

A node is **shared** if it groups exons from `≥ 2` distinct copies, and
**copy-specific** otherwise.

Each copy `c` is mapped to a **path** `π(c) ⊆ V`, the sequence of nodes that
its canonical transcript visits. Multi-isoform information is preserved as
`T(c) ⊆ Paths(G)`, the set of distinct paths corresponding to all annotated
isoforms of `c`.

---

## 2. Read-equivalence of copies

**Definition (read-equivalence).** Two copies `c, c' ∈ L` are
read-equivalent, written `c ≈ c'`, iff `T(c) = T(c')` as sets of paths in `G`.

**Proposition.** `≈` is an equivalence relation on `L`:
- *Reflexivity*: `T(c) = T(c)` by definition.
- *Symmetry*: set equality is symmetric.
- *Transitivity*: set equality is transitive.

The **quotient** `L / ≈` is the set of **equivalence classes** of read-
indistinguishable copies.

### Theorem 1 (Resolvability bound)

Let `c ≈ c'`. Then no read whose evidence consists of (a) the set of nodes it
maps to and (b) the set of edges (junctions) it spans can distinguish `c`
from `c'`.

*Proof sketch.* A read `r` is consistent with copy `c` iff `r`'s evidence
embeds into some path of `T(c)`. Since `T(c) = T(c')` by hypothesis, the same
set of paths is available for both copies, so `r` is consistent with both. □

**Corollary** (collapse necessity). Copies in the same equivalence class must
be quantified jointly. The minimum unit of resolution is the equivalence class.

### Lemma 1 (Witness)

If `c ≁ c'`, there exists a node `v ∈ V` or an edge `e ∈ E` that appears in
some path of `T(c)` but no path of `T(c')` (or vice versa). Such a `v` or `e`
is a **witness** for the distinction; any read mapping to it (or spanning it)
discriminates `c` from `c'`.

*Proof.* If `T(c) ≠ T(c')`, the set-symmetric-difference of their unions of
nodes (or edges) is nonempty. Any element of the difference is a witness. □

The function `class_witnesses(eq)` (in `02_build_graphs.R`) computes the
witness sets for every pair of distinct classes.

---

## 3. Component decomposition: extending to divergent paralogs

When paralogs in `L` diverge past the sequence-similarity threshold,
`G(L)` decomposes into multiple weakly-connected components. Let

```
G(L) = G_1 ⊔ G_2 ⊔ ... ⊔ G_K     (K ≥ 1)
```

with `n_k` denoting the number of copies whose path `π(c)` lies entirely in
`G_k` (so `Σ_k n_k = N`).

**Observation.** Because each `π(c)` is a connected path in a DAG, `π(c)` is
contained in exactly one component `G_k`. Therefore the family `L` partitions
into `L = L_1 ⊔ ... ⊔ L_K` where `L_k = {c : π(c) ⊆ G_k}`.

### Theorem 2 (Component-local quantification)

Reads are assigned within a single component. A read `r` mapping to node
`v ∈ V(G_k)` can only be assigned to copies in `L_k`. Cross-component
assignment never occurs.

*Proof.* A read is consistent with a copy `c` only via path-embedding (cf.
Theorem 1). Paths in `T(c) ⊆ Paths(G_k)` for some unique `k`, and `r`'s
evidence is in `V(G_k) ∪ E(G_k)`. The two constraints intersect iff
`c ∈ L_k`. □

**Corollary.** The equivalence-class quotient distributes over component
decomposition:

```
L / ≈ = (L_1 / ≈_1) ⊔ ... ⊔ (L_K / ≈_K)
```

so per-component analysis is sufficient and complete.

---

## 4. Taxonomy of multi-copy gene loci

Based on `(K, n_1, ..., n_K)`, every multi-copy locus falls into one of three
types:

| Type | Condition | Example |
|------|-----------|---------|
| **Coherent family** | `K = 1` and `n_1 ≥ 2` | RBMY, AMY, GOLGA8I, GOLGA8NL, GOLGA6C, GAGE |
| **Mixed cluster** | `K > 1` and some `n_k ≥ 2` | MAGEA (4 comps, 1 has 3 copies), APOL (3 comps, 1 has 4 copies) |
| **Fully fragmented cluster** | `K = N` (every copy is its own component) | DEFB (10/10), PRAMEF (4/4) |

All three are valid quantification targets. The pipeline distinguishes them
algorithmically via `analyse_locus_cluster(vg)`.

### Strict predicate

```
is_multicopy_family(vg) :=
    N ≥ 2  ∧  has shared and copy-specific nodes  ∧  K = 1
  ∧ ≥ 20% of nodes are shared
  ∧ max copies sharing a node ≥ ⌈N/2⌉
  ∧ |L / ≈| ≥ 2
```

(Used when the goal is a single connected VG.)

### Broad predicate

```
is_multicopy_locus(vg) :=
    N ≥ 2  ∧  type(vg) ∈ {coherent_family,
                          mixed_cluster,
                          fully_fragmented_cluster}
```

(Used when paralogs may be divergent.)

---

## 5. Empirical validation (10 gorilla loci)

Result of running the full classification on ten loci:

| Locus | N | K | classes | resolved | type |
|-------|---|---|---------|----------|------|
| RBMY      | 13 | 1  | 10 | partial | coherent_family |
| AMY       | 3  | 1  | 3  | yes     | coherent_family |
| GOLGA8I   | 7  | 1  | 7  | yes     | coherent_family |
| GOLGA8NL  | 6  | 1  | 6  | yes     | coherent_family |
| GOLGA6C   | 3  | 1  | 3  | yes     | coherent_family |
| GAGE      | 12 | 1  | 3  | partial | coherent_family |
| MAGEA     | 5  | 4  | 5  | yes     | mixed_cluster |
| APOL      | 6  | 3  | 6  | yes     | mixed_cluster |
| DEFB      | 10 | 10 | 10 | yes     | fully_fragmented_cluster |
| PRAMEF    | 4  | 4  | 4  | yes     | fully_fragmented_cluster |

**Resolution dial.** Rebuilding RBMY at `min_identity = 0.95` (default `0.85`)
increases nodes from 19 → 31 and resolves all 13 copies into 13 distinct
classes. The threshold trades VG complexity for resolvability.

---

## 6. Quantification protocol

Given a locus `L` and its VG `G(L)`:

1. Run `analyse_locus_cluster(vg)` → `(type, components)`.
2. For each component `G_k`:
   - Compute `L_k / ≈_k` via `family_equivalence_classes(sub_vg_k)`.
   - Per-class HMMs are built from the union of node sequences spanned by
     any path in any member's `T(c)`.
3. For each read `r`:
   - Identify the component(s) `G_k` whose nodes `r` aligns to (exactly one,
     by Theorem 2).
   - Score `r` against all class HMMs in `L_k / ≈_k`.
   - Hard-assign by best Δlog-likelihood when divergent; soft-assign via EM
     when tied. (See decision tree in `figures/fig_hmm7_golga8i_real.png`.)
4. Report abundance per equivalence class; sum within class to get
   per-class abundance.

---

## 7. Files

- `02_build_graphs.R` — `family_equivalence_classes`, `class_witnesses`,
  `analyse_locus_cluster`, `is_multicopy_family`, `is_multicopy_locus`,
  `build_variation_graph_seq`
- `12_methodology_review.R` — validation across AMY/GOLGA/NPIP/RBMY
- `13_subfamily_vgs.R` — GOLGA sub-family decomposition
- `17_validate_path_definition.R` — pairwise mechanism analysis
- `18_expand_families.R` — DEFB / GAGE / MAGEA / APOL / PRAMEF construction
- `19_family_taxonomy.R` — taxonomy figure across all three types
- `figures/fig_family_taxonomy.png` — visual taxonomy
- `figures/fig_hmm7_golga8i_real.png` — read-assignment decision tree
- `figures/expanded_families_validation.tsv` — full validation table
