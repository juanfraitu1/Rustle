# POA contiguous-core coverage family criterion

> *Naming note:* the separating axis is **contiguous-core coverage = largest ungapped block / shorter gene**. Because it divides by the SHORTER gene it equals the **maximum** of the two per-gene block ratios, so it is NOT a reciprocal coverage (a reciprocal = the *min* over both genes, i.e. divide by the LONGER gene). The word *reciprocal* applies only to the all-column literal metric `min(aligned/len_a, aligned/len_b)`. The two metrics use opposite denominators; only the all-column one is reciprocal.

> **POA contiguous-core coverage (largest ungapped block / shorter gene) >= 0.13 accepts 5/5 confirmed tandem dups and rejects 7/7 domain-sharers.**

## Definition under test

Two loci are in a multi-copy *transcript* family iff, aligned, the homologous region covers >= T of BOTH transcripts:

```
reciprocal_coverage(a,b) = min( aligned_len/len(a), aligned_len/len(b) )
```

i.e. the two genes must CO-ALIGN over a MAJORITY of *each* gene, not just one shared domain/exon. Both quantities (reciprocal aligned coverage, identity over aligned columns) are read off the alignment and are POA-derived. No DNA/protein domain annotation, no BLAST, no k-mer/minimizer is used as the criterion (minimizer-Jaccard is recomputed only to CONTRAST).

### POA-faithfulness note

POA generalises pairwise alignment to >2 copies by embedding each new sequence into a DAG of the previous ones; for a PAIR there is no DAG to extend, so the two-sequence POA instance IS a single pairwise global alignment. We use `Bio.Align.PairwiseAligner` in **global** mode (Needleman-Wunsch) with match=+2, mismatch=-1, gap-open=-5, gap-extend=-1 on the gene-representative sequences.

## Two POA coverage axes (the key result)

We measured coverage two ways, BOTH purely alignment-column-derived:

1. **all-column** reciprocal coverage -- the LITERAL definition above, summing every ungapped aligned column; and
2. **contiguous-core** coverage -- the largest *single* ungapped co-aligned block divided by the **shorter** gene (`biggest_block / min(len)`). This is NOT a reciprocal coverage: dividing by the shorter gene makes it the *maximum* of the two per-gene block ratios (a reciprocal would divide by the longer gene = the *min*). It is the metric that separates; the name 'reciprocal' belongs only to axis 1.

The literal all-column metric does **NOT** cleanly separate the classes: confirmed range [0.471, 0.842], domain-sharer range [0.268, 0.843] -- they OVERLAP. The reason is a global-alignment artifact: when the Needleman-Wunsch aligner is handed two NON-homologous transcripts it still pays a little gap cost to harvest scattered chance matches, fragmenting the alignment into hundreds of tiny blocks (domain-sharers here split into 89-579 blocks) whose aligned columns *sum* to a high reciprocal coverage even though the alignment is non-homologous filler (note the low all-column identity of those pairs in the per-pair table).

The **contiguous-core** metric removes that filler and is the POA-only fix. It is **SEPARABLE**: confirmed range [0.174, 0.608], domain-sharer range [0.008, 0.055], control range [0.007, 0.083]. A real copy has ONE large homologous block; a domain-sharer's largest single block covers only a few percent of the shorter gene -- indistinguishable from random cross-family controls.

## Result: the separating threshold

The principled, POA-only threshold (on the contiguous-core axis) is **T = 0.13** (midpoint between the lowest confirmed and the highest domain-sharer/control coverage -- they are strictly separated).

| category | n | accept at T (cov>=T) | reject at T (cov<T) |
|---|---|---|---|
| confirmed | 5 | **5** | 0 |
| domain-sharer | 7 | 0 | **7** |
| control (cross-family) | 60 | 0 | 60 |

Is the distribution bimodal/separable? **YES -- clean bimodal split (~3.2x vs domain-sharers, ~2.1x vs the nearest control).** (contiguous-core axis: confirmed >= 0.174 vs domain-sharer <= 0.055 and control <= 0.083; the binding barrier is the nearest control at 0.083.)

## Contrast: minimizer-Jaccard does NOT separate them as well

Recomputing the OLD criterion (canonical minimizer-Jaccard, k=15/w=10, blake2b-64 -- identical to `family_detection_validation.py` and rustle production) for the SAME labeled pairs:

| criterion | threshold | confirmed accepted | domain-sharers rejected | clean separation? |
|---|---|---|---|---|
| **POA contiguous-core coverage** | T=0.13 | 5/5 | 7/7 | **YES** |
| POA all-column coverage (literal) | T=0.91 | 0/5 | 7/7 | no (overlap) |
| minimizer-Jaccard (best-fit) | J>=0.278 | 1/5 | 6/7 | no (overlap) |
| minimizer-Jaccard (production 0.30) | J>=0.30 | 1/5 | 6/7 | - |

Even at its *best-fit* threshold the minimizer-Jaccard cannot cleanly separate the two classes -- confirmed and domain-sharer Jaccard ranges overlap, which is exactly why the similarity-only grouping conflated them. Confirmed Jaccard range [0.126, 0.322], domain-sharer Jaccard range [0.098, 0.417]. (Note RFPL1/2/3 -- genuine recent dups -- have Jaccard ~0.13-0.17, LOWER than several domain-sharers, so no Jaccard cutoff can both keep them and drop the domain-sharers.)

*Display-rounding note:* GPR39<->LYPD1's Jaccard is 0.2779, which rounds to 0.278 -- the same 3dp as the best-fit threshold J>=0.278. It sits just BELOW the threshold, so it is rejected (x < t); the apparent tie is only a rounding artifact in the table.

## Per-pair table

| pair | label | contiguous-core cov | all-column cov | aligned identity | minimizer-Jaccard | n blocks | biggest block | len A | len B | same univ. family |
|---|---|---|---|---|---|---|---|---|---|---|
| RFPL1 <-> RFPL3 | confirmed | **0.608** | 0.842 | 0.833 | 0.156 | 37 | 898 | 1476 | 1722 | yes |
| RFPL1 <-> RFPL2 | confirmed | **0.606** | 0.611 | 0.896 | 0.126 | 72 | 895 | 1476 | 2414 | yes |
| RFPL2 <-> RFPL3 | confirmed | **0.520** | 0.711 | 0.879 | 0.165 | 65 | 895 | 2414 | 1722 | yes |
| APOBEC3D <-> APOBEC3F | confirmed | **0.302** | 0.471 | 0.907 | 0.140 | 20 | 860 | 2846 | 4843 | yes |
| RABL2A <-> RABL2B | confirmed | **0.174** | 0.814 | 0.652 | 0.322 | 130 | 416 | 2715 | 2395 | yes |
| CASP8 <-> FLACC1 | domain-sharer | **0.055** | 0.613 | 0.663 | 0.200 | 259 | 163 | 2958 | 4795 | yes |
| ASDURF <-> ASNSD1 | domain-sharer | **0.031** | 0.330 | 0.765 | 0.098 | 89 | 21 | 675 | 2037 | yes |
| GPR39 <-> LYPD1 | domain-sharer | **0.020** | 0.651 | 0.646 | 0.278 | 191 | 42 | 3259 | 2142 | yes |
| CCDC188 <-> ZDHHC8 | domain-sharer | **0.014** | 0.661 | 0.655 | 0.220 | 320 | 53 | 3751 | 5638 | yes |
| CDPF1 <-> PPARA | domain-sharer | **0.012** | 0.640 | 0.632 | 0.160 | 579 | 79 | 6706 | 10397 | yes |
| CREB1 <-> METTL21A | domain-sharer | **0.010** | 0.843 | 0.577 | 0.417 | 430 | 68 | 8050 | 7022 | yes |
| GCA <-> KCNH7 | domain-sharer | **0.008** | 0.268 | 0.822 | 0.175 | 572 | 33 | 4136 | 15456 | yes |

Control pairs (random cross-family, n=60, seed=1729): contiguous-core coverage median 0.033, max 0.083; all-column coverage median 0.427; minimizer-Jaccard median 0.000.

## Interpretation

On the POA **contiguous-core** axis the distribution is **bimodal and cleanly separable** (~3.2x vs domain-sharers, ~2.1x vs the nearest control): every confirmed tandem/recent-duplicate pair has ONE large homologous co-aligned block (>= 0.174 of the shorter gene), while every domain-sharer's largest single block covers only a few percent (<= 0.055) -- indistinguishable from random cross-family controls (<= 0.083). A single principled threshold T = 0.13 accepts ALL confirmed and rejects ALL domain-sharers AND all controls.

The LITERAL all-column reciprocal coverage does NOT achieve this, because a global alignment of two non-homologous transcripts inflates coverage with gappy chance-match filler (visible as low all-column identity and hundreds of tiny blocks). Requiring the co-alignment to be a contiguous homologous CORE -- still a pure alignment statistic, no domain/BLAST/k-mer input -- is the **POA-only fix for the domain-sharing false positives.** The minimizer-Jaccard cannot do this at all: it rewards any shared subsequence, so a shared domain inflates the intersection exactly like whole-gene homology, leaving confirmed and domain-sharer Jaccard values overlapping.

## Honest caveats

- **Operational RNA-structural definition, NOT the gene family.** This defines a *POA-coalignable multi-copy transcript family* (do the mature transcripts co-align over a majority of each?). It is deliberately NOT claimed equal to the DNA/protein gene family from a gene tree -- proving that equivalence would need DNA/protein/domain evidence, which the constraint forbids. It is an RNA-level, alignment-only operational criterion.
- **Pairwise-global as the POA pairwise instance.** For two sequences the POA instance reduces exactly to a single global alignment; for >2 copies a true POA DAG would be used and could differ slightly at multi-copy column assignment. The pair tests here are the exact two-sequence case.
- **All-column vs contiguous-core coverage.** The clean separation is on the contiguous-core axis, NOT the literal all-column reciprocal coverage (which overlaps because a global aligner pads non-homologous pairs with gappy chance-match filler). Both are alignment-column statistics and thus both POA-only/within-constraint, but the operational criterion that WORKS is specifically 'largest single UNGAPPED co-aligned block >= T of the shorter gene'. NOTE this is the largest *ungapped* block, not a local-alignment span: a Smith-Waterman local alignment of these domain-sharers ALSO stitches chance matches and reports a high aligned-length (verified: domain-sharer local aligned-len/minlen reaches 0.91-0.97, overlapping the confirmed pairs) -- so it does NOT separate. It is the longest CONTIGUOUS gap-free homologous run that distinguishes a real copy (one long block) from a domain-sharer (many short blocks patched together). That largest-ungapped-block read-off is what the criterion uses.
- **Small labeled N.** Only the Compara-checkable named pairs are labeled: 5 confirmed and 7 domain-sharers. The separation is an exact count on a small sample, a directional sanity check, not a population rate.
- **Gene-representative sequences.** Coverage is computed on one representative sequence per gene (`gene_rep.fa`), not all isoforms; a different representative isoform could shift coverage at the margin.
- **RABL2A<->RABL2B label (assumption is doing real work).** Treated as confirmed-real per task spec because it is the flagship 2-copy tandem pair. Compara returned a fetch error for RABL2B, but note that RABL2A's OWN Compara query SUCCEEDED (para_status=ok, 68 paralogues) and does NOT list RABL2B's ENSG (ENSG00000079974). So the 'confirmed' label here rests on the task-spec/flagship assumption, not merely on a one-sided fetch error -- RABL2A's own successful list also omits RABL2B. Reassuringly, RABL2A<->RABL2B is the LOWEST confirmed pair (contiguous-core 0.174), yet still cleanly above T=0.13, so it does not prop up the separation.
- **Labels are external corroboration only.** Compara is used purely to *check* the POA criterion's verdicts; it is never an input to the criterion. The criterion is alignment-derived end to end.
- **Determinism.** Controls use a fixed RNG seed (1729); alignment and minimizer machinery are deterministic. Output is byte-stable.

