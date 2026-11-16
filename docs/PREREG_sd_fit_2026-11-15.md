# PREREG row 1755: do segmental duplications fit multi-copy gene families, and does exon correspondence fit better than the span rule (2026-11-15)

User question (2026-11-15): multi-copy families should be more likely inside SDs; how to prove it and how to use it better.
Already measured (not repeated): `sdbio.py` 2026-11-08 (80.5% of genes with a young paralog in an SD vs 4.5% of the rest;
91% of young Compara pairs SD-linked), row 1580 damping (default), rows 1581-1687 other SD uses (closed).

Annotation level only: CHM13 CAT protein-coding genes, human Compara pairs (young = human..Simiiformes last common ancestor),
native CHM13 SEDEF `final_human.bed` (identity col 21, CIGAR col 33). No Rustle loci, so no circularity with SD damping.
Script `bench/denovo_eval/sdfit.py`, workspace `r1755_sdfit`. Written before any number below was seen.

## A. Is the enrichment more than "similar genes are SDs by definition"? (decides the wording of the claim)

- A1 flank. For each span-linked young pair take, over its holding SD pairs, the largest min(flank on copy 1, flank on copy 2),
  flank = SD copy length minus its overlap with the gene span. **PASS: median >= 5,000 bp** (the duplicated segment is larger
  than the gene). FAIL: the SD is the gene itself; the enrichment claim must then be worded as sequence similarity.
- A2 null (descriptive, with empirical p): SD intervals circularly shifted per chromosome, 100 shifts (seed 1), genes fixed.
  Reported: share of young-paralog genes inside an SD (span rule), observed vs shifted.

## B. What are the young pairs that are not SD-linked? (descriptive)

Classes, first match: partial (both genes have exon bases in the two copies of one SD pair, span rule fails); one-sided
(only one gene, or both in unrelated SDs); retrocopy-like (one single-exon, one multi-exon, no SD); no SD at all.

## C. Exon correspondence through the SD alignment vs the span rule (decides whether a families step 0 is worth running)

Candidates: protein-coding gene pairs with exon bases in the two copies of one SD pair. CORR = exon bases of gene a that the
SEDEF CIGAR projects onto exon bases of gene b, over the smaller of the two genes' exon bases; best over SD pairs.
Rule C: CORR >= 0.5 (no span condition). Compared with the span rule (as `--sd-damp`) on the same genes, truth = young pairs.

- **PASS: precision of rule C >= span-rule precision + 0.10, and rule C reaches at least as many young truth pairs.**
  Then: next row = rule C as the damping link in the families stage, held out as row 1580 was.
- FAIL: the span rule stays; exon correspondence through SDs is closed.

Also reported, deciding nothing: precision against Compara pairs of any age, CORR thresholds 0.2 / 0.8.
