# Multi-Copy Gene Family Experiment (GGO)

## Three families tested

| Family | Chromosome | Copies tested | Primary reads | Secondary reads |
|--------|-----------|---------------|---------------|-----------------|
| AMY    | chr1 (NC_073224.2) | AMY1B, AMY2A, AMY2B | 2+31+120 = 153 | 5+7+2 = 14 |
| GOLGA6L7 | chr19 (NC_073243.2) | L7_#1, L7_#2, L7_#3 | 26+6+5 = 37 | 55+70+70 = 195 |
| GOLGA8 | chr15 (NC_073240.2) | 8B-like, 8N_...21 | 130+3 = 133 | 0+34 = 34 |

Observation: the tighter the paralogs (GOLGA6L7 is the tightest family),
the higher the fraction of SECONDARY alignments — i.e., reads the aligner
considered 'equally good' at multiple copies. AMY copies have diverged enough
that secondary alignments are rare.

## The family VG for GOLGA6L7 (identical topology, shifted coordinates)

| intron | GOLGA6L7_1 | GOLGA6L7_2 | GOLGA6L7_3 | length |
|--------|-----------|-----------|-----------|--------|
| i01 | 104789919-104791223 | 104830738-104832042 | 104871546-104872850 | 1305 |
| i02 | 104791353-104792166 | 104832172-104832985 | 104872980-104873793 | 814 |
| i03 | 104792197-104792465 | 104833016-104833284 | 104873824-104874092 | 269 |
| i04 | 104792517-104792604 | 104833336-104833423 | 104874144-104874231 | 88 |
| i05 | 104792699-104792780 | 104833518-104833599 | 104874326-104874407 | 82 |
| i06 | 104792888-104794129 | 104833707-104834952 | 104874515-104875758 | 1242 |
| i07 | 104794189-104794517 | 104835012-104835340 | 104875818-104876146 | 329 |
| i08 | 104794660-104794953 | 104835483-104835776 | 104876289-104876582 | 294 |

**Interpretation**: all 3 copies share the same 9-exon, 8-intron topology. The
length of each intron is identical across copies (modulo tiny indels at i06-i08).
The only difference is the absolute START coordinate, which shifts by -70bp
(L7_1 → L7_2) and -82bp (L7_1 → L7_3). The copies are a near-perfect tandem
triplication with minor drift.

## Assembler performance per copy

| Region | gt mRNAs | primary reads | Rustle predictions | StringTie predictions | = matches |
|--------|---------|---------------|--------------------|-----------------------|-----------|
| AMY1B | 1 | 2 | 0 | 0 | 0 |
| AMY2A | 2 | 31 | 7 | 7 | 0 |
| AMY2B | 3 | 120 | 10 | 13 | 0 |
| GOLGA6L7_1 | 1 | 26 | 4 | 4 | 1 |
| GOLGA6L7_2 | 1 | 6 | 0 | 0 | 0 |
| GOLGA6L7_3 | 1 | 5 | 1 | 1 | 0 |
| GOLGA8B-like | 1 | 130 | 7 | 11 | 0 |
| GOLGA8N_x21 | 2 | 3 | 0 | 0 | 0 |

## Why multi-copy assembly is hard (GOLGA6L7 case)

| Copy | Primary reads | Secondary reads | Fate |
|------|---------------|-----------------|------|
| L7_1 | 26 | 55 | 2 predictions on + strand, 1 perfect match (8/8 introns) |
| L7_2 | 6  | 70 | **No predictions** — aligner gave primaries to L7_1/L7_3 |
| L7_3 | 5  | 70 | Single prediction on wrong strand (neighbor gene dominates) |

The aligner arbitrarily picks ONE copy as primary per read. For tandem paralogs,
this concentrates reads on L7_1 and starves L7_2/L7_3. A read-to-graph aligner
(vg giraffe) or a post-hoc EM redistribution can balance these.

## Rustle VG mode result on GOLGA6L7 cluster

After extending the multi-map scan to include secondary alignments:

- **Cross-bundle links found**: 217 multi-map alignments → 12 linked pairs
- **Family groups discovered**: 1 (covering 2 stranded sub-bundles in the L7 region)
- **EM iterations to converge**: 2 (with delta=0 — convergence is trivial)
- **Reads reweighted**: 12
- **Final GTF change**: None — predictions unchanged after EM

**Why VG mode does not yet help**: when paralogs are close enough (GOLGA6L7
spans ~85 kb), Rustle's bundle builder merges them into ONE bundle. VG mode
operates at the bundle boundary, so same-bundle paralogs are invisible to it.
To handle tandem paralogs we need either (a) a tighter bundle splitter that
breaks at paralog boundaries, or (b) an intra-bundle sub-cluster step that runs
assembly per paralog.

On the full chr19 benchmark: VG mode finds 23 families and reweights 416 reads,
but yields -4 matches relative to baseline (1468 vs 1472). The EM converges
trivially (all EMs in 2 iterations with delta=0) which means there is no
informative signal — most multi-mappers have identical compatibility across copies.

## What an assembled transcript looks like, projected on the family VG

- **RSTL.1.1** (8 introns, 8/8 family-matched)
  - ✓ L7_1:i01
  - ✓ L7_1:i02
  - ✓ L7_1:i03
  - ✓ L7_1:i04
  - ✓ L7_1:i05
  - ✓ L7_1:i06
  - ✓ L7_1:i07
  - ✓ L7_1:i08
- **RSTL.1.2** (7 introns, 6/7 family-matched)
  - ✓ L7_1:i01
  - ✗ novel:104791353-104792465
  - ✓ L7_1:i04
  - ✓ L7_1:i05
  - ✓ L7_1:i06
  - ✓ L7_1:i07
  - ✓ L7_1:i08

## Conclusion: how to compare family members

1. **At the topology level**: the family VG nodes/edges are canonical intron-exon
   positions. Every family member shares this abstract graph. That's the 'variation
   graph' the family defines.
2. **At the coordinate level**: each VG edge is decorated with {start, end} tuples
   per copy. The differences between tuples are indel drift (~70-100 bp).
3. **At the assembly level**: the classical per-copy primary-read approach fails
   for tightly-clustered paralogs because the aligner's primary pick is essentially
   random across near-identical sequences. Either (a) we split bundles at paralog
   boundaries and use multi-mapping reads to populate each copy, or (b) we assemble
   on the family VG directly and project paths back onto copies.
