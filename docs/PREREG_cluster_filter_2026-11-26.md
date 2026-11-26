# PREREG: the 10-point lead between CLUSTERS of copies, on one reference (row 1974, 2026-11-26)

The user expects the advisor to ask for the Eichler lab's filter to be adapted. What carries over is its UNIT: the paper
compares paralog clusters and discards the rest (row 1972, `docs/EICHLER_FILTER_WHAT_IT_DOES_2026-11-26.md`). Our
single-reference version asked for a lead between individual copies and fires on 8.7% of HG002's near-tied reads (row
1954). Here the lead is asked between clusters: a read whose tied alignments all lie in one cluster is decided for that
cluster, whichever of its copies it came from.

## Data and truth
HG002 family reads on CHM13 (`r1776_hg002_room/mm.sel.bam`, `-N 50 -p 0.1`, every record). Truth for a read = the paralog
the paper's procedure retains on HG002's own haplotypes (row 1972; reads it marks ambiguous have no truth).

## Clusters of the CHM13 family genes (no phylogeny is used)
- **pipeline** (the registered one): genes whose loci the pipeline's `families` stage puts in one family (HG002 run
  `r1942_ms_seed_pool/hg.A0`); a gene in no family is its own cluster. This is what a de novo run has.
- reported beside it: **family** (the annotation's whole family, the coarsest), **99%** and **95%** (genes joined when
  their genomic sequences align over half of the shorter at that identity or more; single linkage).

## Rule (`bench/o2_truth/cluster_filter.py`), on `ms` and on AS
Each alignment of a read is given the family gene under it (>= 50 exon bases) and that gene's cluster. Per cluster the
read's best score. The cluster with the best score is RETAINED when it leads every other cluster, and every alignment on
no family gene, by 10 or more; otherwise the read is AMBIGUOUS. With one cluster per gene this is the per-copy rule.

## Measures
Over the reads that are near ties per copy on `ms` (another alignment within 10), and over all family reads: share
retained; of the retained with a truth, the truth paralog lies in the retained cluster; and the size of the retained
cluster (copies), so that "decided" is not bought with a cluster as wide as the family. TBC1D3 and NPIP on their own.

## Claims (pipeline clusters, `ms`)
- **A1:** of the reads near-tied per copy, at least 50% are retained for a cluster (per copy: 8.7% on AS, row 1954).
- **A2:** of those retained reads with a truth, the truth paralog is in the retained cluster for at least 95%.
- **A3:** the median retained cluster of those reads holds at most 5 copies.
## Predictions (mine)
A1 holds (60-80%). A2 holds. A3 uncertain: the pipeline's families may be as wide as the annotation's for the large
arrays; for TBC1D3 I expect one cluster holding the whole array, which decides every read and says little.

## Result (job 4510290). Truth: 114,850 reads retained by the paper's procedure on the haplotypes (1,249 ambiguous, no truth).
Clusterings over 1,459 family genes: copy 1,459 clusters; 99% identity 1,017 (largest 52); 95% 819 (largest 59);
pipeline 1,293; annotation family as given.

### The 7,051 reads that are near ties per copy on `ms`
| clusters | retained (ms) | truth paralog in the retained cluster | copies in the retained cluster, median / mean | retained (AS) | correct (AS) |
|---|---|---|---|---|---|
| copy (the per-copy rule) | 14.7% | 79.3% (630 of 794) | 1 / 1.0 | 16.1% | 74.9% |
| 99% identity | 78.4% | 98.0% (4,611 of 4,704) | 2 / 2.3 | 79.2% | 97.4% |
| pipeline families | 81.0% | 99.2% (4,851 of 4,891) | 2 / 2.6 | 81.9% | 98.8% |
| 95% identity | 85.9% | 99.8% | 2 / 3.9 | 86.1% | 99.7% |
| annotation family | 86.9% | 99.9% | 2 / 4.5 | 87.1% | 99.8% |
### All 135,640 reads with the primary on a family gene (ms)
copy 95.5% retained, 93.4% correct; 99% identity 98.8%, 96.3%; pipeline 99.0%, 95.1%; 95% identity 99.2%, 96.9%;
annotation family 99.3%, 99.8%.
### The two arrays
| | TBC1D3 (151 reads) | NPIP (1,715 reads) |
|---|---|---|
| per copy: retained / correct | 41.7% / 0 of 5 | 71.4% / 54.6% |
| 99% identity: retained / correct / copies in the cluster | 100% / 52 of 52 / 9 | 77.5% / 98.5% / 6 |
| pipeline families | 100% / 52 of 52 / 9 | 73.6% / 70.8% / 7 |
| 95% identity and annotation family | 100% / 52 of 52 / 9 | 79.8% / 99.9% / 23 |

### Claims (pipeline clusters, `ms`)
- **A1 HOLDS:** 81.0% of the reads near-tied per copy are retained for a cluster (per copy 14.7%).
- **A2 HOLDS:** the truth paralog is in the retained cluster for 99.2%.
- **A3 HOLDS:** the median retained cluster has 2 copies (mean 2.6).

### Reading
- **This is the part of the paper's idea that carries to one reference.** Asked between clusters of near-identical
  copies, the 10-point lead decides four near-tied reads in five and is right in 98-99%; asked between copies it decides
  one in seven and is right in 79%. Most ties are between two copies that are the same sequence for the read.
- **The price is resolution, stated per read:** the answer is "one of these 2 copies" in the median, "one of these 9"
  for TBC1D3 (the whole array is one cluster under every definition: every read is decided and none is placed), "one of
  these 6" for NPIP at 99% identity.
- **Which clustering:** the pipeline's own families work overall (99.2%) and fail on NPIP (70.8%: they split copies the
  reads cannot tell apart); clusters by sequence identity of the copies at 99% hold everywhere measured (98.0% overall,
  98.5% NPIP, 100% TBC1D3) at about the same resolution. An identity clustering needs the copies' sequences, which a
  run has (its own loci), not an annotation; here it was computed on annotated genes.
- `ms` and AS give the same picture (within one point); `ms` is slightly more often right.
- Limits: HG002 only; truth = the paper's procedure with paralog = nearest reference copy (row 1972's limits apply);
  clusters built on annotated genes, not yet on a run's own loci; reads whose truth is ambiguous are not judged.
