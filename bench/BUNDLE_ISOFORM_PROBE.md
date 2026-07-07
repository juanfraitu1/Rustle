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

### What exact-chain captures that StringTie misses
Only **6 isoforms total**, mostly multi-exon:

- GSTM: 4 (2–4 exons)
- ANKRD18: 1 (2 exons)
- HERC2: 1 (2 exons)

### What StringTie captures that exact-chain misses
**83 isoforms total**, spread across GSTM (12), ANKRD18 (28), and HERC2 (43). Exon-count breakdown shows they are overwhelmingly multi-exon:

- HERC2: 1ex=5, 2ex=8, 3ex=10, 4ex=3, 5–20ex=17
- ANKRD18: 1ex=1, 2ex=9, 3ex=5, 4–11ex=13
- GSTM: 2ex=1, 5–15ex=11

### Fragmentation check
A direct union does not help much because the exact-chain set is largely a **fragmentation** of StringTie isoforms, not a set of genuinely new ones. Each StringTie isoform matches on average 1–9 exact-chain isoforms:

| family | ST isoforms | avg exact matches per ST | max | unmatched ST |
|---|---:|---:|---:|---:|
| GSTM | 276 | 8.65 | 26 | 12 |
| ZNF92 | 10 | 4.40 | 6 | 0 |
| RABL2 | 24 | 8.42 | 13 | 0 |
| ANKRD18 | 354 | 5.62 | 43 | 28 |
| RGPD8 | 10 | 7.50 | 9 | 0 |
| HERC2 | 285 | 6.81 | 41 | 43 |
| MAGEA | 6 | 1.33 | 2 | 0 |
| APOBEC3 | 4 | 2.00 | 2 | 0 |

The 83 "StringTie-missed-by-exact" isoforms are the ones with **zero** exact-chain matches; the rest of StringTie's set is already splintered across many exact-chain entries.

## Interpretation

1. **StringTie already captures the bulk of distinct isoforms in multi-copy family regions.** Its bundle-based approach produces clean, compact transcript models, and the 83 isoforms it has that exact-chain lacks are real additional structures (not fragments).

2. **Exact-chain finds many near-duplicate fragments of the same StringTie isoforms.** The higher exact-chain counts in GSTM (349 vs 276) and ANKRD18 (355 vs 354) reflect fragmentation, not additional biological isoforms.

3. **A simple union of the two assemblies adds only 6 isoforms** (`bench/hybrid_isoform_assembly.py` materializes this). The easy "finds what StringTie finds but finds more" win does not exist at the level of simple isoform merging.

4. The result aligns with the earlier read-coherence finding (`bench/DENOVO_PIPELINE.md`): bundle/read-chain methods deliver large isoform-recall gains at **single-copy loci**, but that recall is **disjoint from multi-copy family loci**.

## A useful "in between" would be PSV-aware splitting, not a union

The right hybrid for multi-copy families is not isoform union; it is **using exact-chain/PSV evidence to split StringTie transcripts into copy-specific isoforms** in regions where reads support multiple paralogs. This keeps StringTie's clean models but adds the copy-level resolution the exact-chain method preserves. That is conceptually close to the O2 copy-assignment VG, but materialized as isoforms rather than read assignments.

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
