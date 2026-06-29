# Avenue A — mappability recovery of "empty" paralog copies: SIZED, essentially DRY (2026-06-08)

**Hypothesis.** Some annotated paralog copies get zero reads because minimap2's `-f`
frequent-minimizer filter (or top-N truncation) drops them; re-mapping with Winnowmap /
tuned `-f` would surface reads → more transcripts.

**Result: no meaningful headroom.** The "empty" copies are unexpressed, not unmappable.

## Sizing (3,015 copies in 817 paralog families, band 0.90–0.99)
| pool | copies | % |
|---|---|---|
| has confident primary reads | 2,196 | 73% |
| MAPQ-0-only (reads, no confident primary) | 499 | 17% |
| EMPTY (0 reads) | 320 | 11% |
| └ empty in an *expressed* family | 273 | 9% |

## The 273 candidate pool dissolves under inspection
- **51% (140) are single-exon processed pseudogenes** of highly-expressed genes
  (YWHAZ, RBBP4, HNRNPC…). Untranscribed. **Decisive test:** all 6,471 YWHAZ reads align
  *strictly better to the parent, 0 to its 929 bp pseudogene* (0 tied). Putting reads
  there = fabrication.
- **49% (133) are multi-intron**, but only **4 families** genome-wide have such an empty
  copy beside a co-located, well-expressed (≥20-read) sister — and on full-sensitivity
  re-mapping (`-N 50 -p 0.5`) **0/4 attracted any decisive reads**. The reads decisively
  belong to the expressed sister.
- The 499 MAPQ-0 pool is dominated by identifiability ties (Winnowmap provably can't
  break — see `winnowmap_floor_test.md`), not mappability.

## Why (mechanism)
minimap2 already evaluates reads against **all** gene-copy placements and assigns each to
its best copy. For a read's true home to be an empty copy, that copy would have to be
filtered out of the index — but `-f` only drops the top ~0.02% most-frequent minimizers
(satellite/centromeric repeats), **not** paralog gene copies (even 24-copy arrays). So
there is no reservoir of mis-filtered reads whose true home is an empty copy. Winnowmap /
`-f` tuning help mappability only in satellite/centromeric repeats, which are not
transcript-bearing genes.

## Consequence
The multimapping / VG well for *more transcripts* is at the annotation/identifiability
floor (consistent with the multimapper-treasure and information-borrowing findings).
Genuine "beat StringTie on recall" levers are **non-multimapping**:
- strand-aware bundling → convergent/antisense cut-off genes (validated +30, shippable),
- read-chain (IsoSeq-style collapse) recall,
- sequence-aware structural / gene-conversion emission (a transcript class ST cannot
  produce).
