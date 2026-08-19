# O1 duplication-provenance witness prototype

This experiment instantiates the typed model on fresh HSA families already on disk.
Every block entity is one real passing minimap2 PAF witness. RNA and DNA relations remain
separate; all directions are `UNROOTED`. Pairwise witnesses deliberately remain overlapping
and unmerged, so this is not yet the final non-transitive block-class decomposition.

| case | loci | fresh families | RNA/DNA/both | layer Jaccard | RNA comps | DNA comps | direction |
|---|---:|---|---:|---:|---|---|---|
| GOLGA6_8 | 19 | GWFAM16,GWFAM17,GWFAM18 | 76/84/58 | 0.5686 | 17,1,1 | 19 | UNROOTED |
| MAGEA | 10 | GWFAM33,GWFAM34 | 22/17/17 | 0.7727 | 7,3 | 7,1,1,1 | UNROOTED |
| RABL2 | 2 | GWFAM28 | 1/1/1 | 1.0000 | 2 | 2 | UNROOTED |
| NBPF_core | 19 | GWFAM0,GWFAM1,GWFAM2 | 49/100/42 | 0.3925 | 17,1,1 | 19 | UNROOTED |
| NBPF_repeat_bridge | 15 | GWFAM1,GWFAM15,GWFAM30,GWFAM8 | 26/14/6 | 0.1765 | 11,2,1,1 | 11,1,1,1,1 | UNROOTED |

## Interpretation

A connected component in one layer is not promoted to a family in another layer. In
particular, a DNA block witness records shared duplicated sequence but does not create an
RNA family edge. The GFA files are visual projections; `entities.tsv`, `relations.tsv`, and
the per-case locus/path tables are the normative evidence. Component-size lists include
isolated loci, so `15,1,1` means one 15-locus component and two singletons.

The next implementation step is to consolidate overlapping pairwise witnesses into stable
multi-locus block classes using reciprocal overlap plus quasi-clique cohesion. Until then,
the path tables explicitly say `NOT_YET_CONSOLIDATED`, and no ancestry is inferred.
