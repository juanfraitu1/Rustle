# `complete_poa_core` baseline evaluation

## What was tested
The existing Rust flag `DenovoConfig.complete_poa_core` (exposed as `--complete-core` in
`gw_family_catalog`) was already run on the gorilla BAM (`/home/juanfra/winloci_scratch/GGO.bam`)
with the cross-chromosome read-conflict catalog. Output files:

- Baseline (no completion): `gw_xcbase.{families,copies}.tsv`
- With completion: `gw_comp2.{families,copies}.tsv`

`poa_core_completion_adds` attaches any genome-wide rep whose contiguous POA core to an existing
family member clears `t_core = 0.13`. It is bounded by the minimizer-LSH candidate-pairs prefilter
and only considers free-rep ↔ family-member pairs.

## Evaluation script
`bench/eval_complete_core_baseline.py` maps every copy to a gene via `annot_intervals.tsv`, then
compares the two catalogs on:

- overall family/copy counts
- biotype/annotation of added copies
- coverage of literature-known multi-copy family windows
- coverage of diploid-CN oracle multi-copy genes

## Results

| metric | baseline | complete-core | delta |
|---|---:|---:|---:|
| families | 265 | 265 | 0 |
| copies | 869 | 924 | +53 |
| families with added copies | — | 35 | — |
| protein-coding added copies | — | 46 | — |
| lncRNA added copies | — | 2 | — |
| pseudogene / V_segment / unannotated added copies | — | 5 | — |

### Known-family windows
Only MAGEA and ANKRD18 gained copies inside their windows (+1 each), and in both cases the added
copy does **not** cover a new named gene — it is an additional isoform/fragment inside an already-
represented locus. All other showcase families (RABL2, APOBEC3, RGPD8, ZNF92, GSTM, HERC2) are
unchanged.

### Diploid-CN oracle
Of the 57 named multi-copy oracle genes with annotation coordinates, **51 are covered in both
catalogs**; completion covers **0 new oracle genes**.

### Annotation of added copies
The 53 added copies span 42 distinct gene symbols. They are real, mostly well-supported transcripts,
but they land in loci that are already adjacent to a conflict-graph family rather than recovering
isolated missed paralogs.

## Interpretation
`complete_poa_core` does **not** solve the discovery problem identified in the plan:

> missed diploid-oracle genes and known families have **zero candidate edges** to the current
> family graph.

The flag only attaches copies that are already reachable through the minimizer-LSH prefilter and
already share a contiguous POA core with a family member. It cannot recover a locus that is
structurally invisible to that prefilter (e.g. divergent enough that minimizers don't seed, or a
thin locus below the assembly gate that never became a `DenovoTranscript`/`rep`).

## Recommendation
Proceed to **Approach 1** from the plan: a genome-wide family-aware rescue / edge-generation pass
that scans **thin loci** (support 1–2 reads, below the `≥3` assembly gate) and POA-confirms them
against existing families, feeding confirmed edges back into the family graph. This directly
targets the "no candidate edges" root cause rather than re-scoring edges that already exist.
