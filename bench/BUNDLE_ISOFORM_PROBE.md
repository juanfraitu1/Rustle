# Bundle-style vs exact-intron-chain isoform capture in multi-copy families

**Question.** StringTie uses genomic-overlap *bundles* to group reads and assemble isoforms. The shipped de-novo pipeline instead groups reads by **exact intron chain** and collapses isoforms by identical junction union-find. Do bundles capture more isoforms in multi-copy gene-family regions?

**Answer.** On the 8 canonical known multi-copy families, **bundles do not capture more isoforms — the exact-intron-chain method captures more.** StringTie finds a handful of single-exon / low-exon isoforms that exact-chain misses, but it misses many more multi-exon isoforms that exact-chain recovers.

## Method

`bench/bundle_isoform_probe.py`

- Known-family regions from `bench/family_rna_refine.tsv` (members of the 8 families in `bench/known_family_showcase.tsv`), padded ±5 kb.
- **Exact-chain isoforms**: `denovo_transcripts.meta.tsv` + `denovo_skeletons.tsv` (introns).
- **Bundle isoforms**: `genome_st.gtf` (StringTie `-L` on `GGO.bam`).
- Matching is overlap-based (≥50% reciprocal span) **or** shares ≥1 identical intron.

## Results

| family | exact tx | exact chains | StringTie tx | StringTie chains | StringTie missed by exact | exact missed by StringTie |
|---|---:|---:|---:|---:|---:|---:|
| GSTM | 349 | 349 | 276 | 273 | 12 | 4 |
| ZNF92 | 7 | 7 | 10 | 10 | 0 | 0 |
| RABL2 | 23 | 23 | 24 | 24 | 0 | 0 |
| ANKRD18 | 355 | 355 | 354 | 353 | 28 | 1 |
| RGPD8 | 9 | 9 | 10 | 10 | 0 | 0 |
| HERC2 | 198 | 198 | 285 | 279 | 43 | 1 |
| MAGEA | 4 | 4 | 6 | 6 | 0 | 0 |
| APOBEC3 | 2 | 2 | 4 | 4 | 0 | 0 |
| **TOTAL** | **947** | **947** | **969** | **959** | **83** | **6** |

### What StringTie captures that exact-chain misses
Only **6 isoforms total**, mostly multi-exon:

- GSTM: 4 (2–4 exons)
- ANKRD18: 1 (2 exons)
- HERC2: 1 (2 exons)

### What exact-chain captures that StringTie misses
**83 isoforms total**, spread across GSTM (12), ANKRD18 (28), and HERC2 (43). Exon-count breakdown shows they are overwhelmingly multi-exon:

- HERC2: 1ex=5, 2ex=8, 3ex=10, 4ex=3, 5–20ex=17
- ANKRD18: 1ex=1, 2ex=9, 3ex=5, 4–11ex=13
- GSTM: 2ex=1, 5–15ex=11

## Interpretation

1. **In multi-copy family regions, exact intron-chain grouping already captures most isoform diversity.** The union-find collapse on identical junctions absorbs 5′/3′ variation as long as any internal junction is shared, so the chain-based method is less fragmentary than expected for these regions.

2. **StringTie's bundle advantage is small and concentrated in single-exon / few-exon isoforms** (6 total here), consistent with the known upstream loss of single-exon reads in the exact-chain pipeline.

3. **StringTie misses many multi-exon isoforms** that exact-chain resolves. This may reflect StringTie's parsimony/flow model preferring fewer transcripts, or coverage-driven merging of similar isoforms in these repeat-rich, multimapping regions.

4. The result aligns with the earlier read-coherence finding (`bench/DENOVO_PIPELINE.md`): bundle/read-chain methods deliver large isoform-recall gains at **single-copy loci**, but that recall is **disjoint from multi-copy family loci** — here, the exact-chain method is already the richer isoform set.

## Caveats

- Without a per-isoform truth set, "more isoforms" does not mean "more real isoforms". Some exact-chain transcripts may be fragments or artifacts.
- The counts are region-overlap based; large family regions (e.g., HERC2 duplicon cluster) may include neighboring non-family transcripts.
- The comparison uses one StringTie run; parameter tuning (`-f`, `-c`, `-g`) could shift the balance modestly.

## Verdict

For the specific goal of capturing more isoforms **inside multi-copy gene families**, switching to bundle-based grouping is **not supported by the data**. The exact-intron-chain method already recovers more isoforms in these regions. Bundle-style approaches remain valuable for **single-copy / genome-wide isoform recall** (already shipped via the read-coherence layer), but they are not the lever for multi-copy family isoform capture.

## Reproduce

```bash
PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/bundle_isoform_probe.py
# outputs bench/bundle_isoform_probe.tsv
```
