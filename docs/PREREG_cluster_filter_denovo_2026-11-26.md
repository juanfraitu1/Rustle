# PREREG: the cluster-level 10-point lead with clusters from the run itself (row 1975, 2026-11-26)

Row 1974: asked between clusters of near-identical copies, the 10-point lead decides 78-81% of HG002's reads that are
near ties between copies and is right in 98-99%. That used clusters built on annotated genes, one sample, and a threshold
(99% identity) read off that sample. The user: can it be generalized; run the test.

## What changes
- **Clusters from the run's own copies:** the sequences the `families` stage writes (`PREFIX.fam.copies.fa`, one per
  copy), compared with each other (`minimap2 -cx asm20 -X`), joined at 99% identity over at least half of the shorter
  (single linkage), ACROSS the pipeline's family calls (row 1974: those split NPIP). A locus with no copy sequence is
  its own cluster. No annotation is used to build them.
- **The verdict from the run's own loci:** each alignment of a read gets the assembled locus under it (>= 50 exon bases)
  and that locus's cluster; the cluster with the best score is retained when it leads every other cluster, and every
  alignment on no locus, by 10 (`bench/o2_truth/cluster_filter_run.py`).

## Data and truths
1. **HG002** (run `r1942_ms_seed_pool/hg.A0`, reads `r1776_hg002_room/mm.sel.bam`): truth = the paralog of row 1972's
   haplotype procedure, carried to the run's locus under it. Same reads as row 1974, new clusters.
2. **Planted reads with private variants** (new seeds 1975 and 1976, `plant.py --private 0.001`, scenarios E and S, 60
   reads per isoform): the sample is NOT the reference; truth = the planted copy. Unseen data for the 99% threshold.
3. **KB3781** (run `r1937_shared_floor1/kb`, truth `r1785_truth_repair/kb/reads.tsv`, a sample of family reads with
   their candidate loci and AS; the rule on AS, the table has no `ms`): another individual, another species.

## Claims (clusters at 99% identity; `ms`, AS on KB3781), over the reads that are near ties per locus
- **G1, HG002:** retained at least 65%; of the retained with a truth, the truth's locus is in the retained cluster for
  at least 95%.
- **G2, planted with private variants:** in each of the four sets retained at least 50% and right at least 95%.
- **G3, KB3781:** right at least 95% of the retained with a truth (the retained share is reported).
Beside each: the per-locus rule and clusters at 95%, and the size of the retained cluster.
## Predictions (mine)
G1 holds with fewer reads decided than on annotated genes (65-75%: the run's loci are not the annotation's). G2 holds.
G3 holds. If G1-G3 hold, building the verdict into the stage is justified; the default stays the user's decision.

## Result (jobs 4510302-4510305). Reads that are near ties per locus; clusters from the run's own copy sequences.
| data | per locus: decided / right | 99% clusters: decided / right | loci in the retained cluster, median / mean | 95% clusters: decided / right |
|---|---|---|---|---|
| HG002 (7,119 reads; ms) | 15.9% / 80.1% | 80.0% / 96.6% (4,660 of 4,823) | 2 / 2.1 | 80.1% / 96.6% |
| KB3781 (1,321 reads; AS) | 7.6% / 97.9% (46 of 47) | 74.3% / 99.3% (734 of 739) | 2 / 2.1 | 75.1% / 98.7% |
| planted, private variants, seed 1975 E (5,086; ms) | 11.0% / 100% | 83.3% / 100% | 7 / 12.4 | 95.6% / 100% |
| planted, seed 1975 S (2,790) | 11.0% / 99.7% | 83.5% / 100% (2,330 of 2,331) | 6 / 11.4 | 95.4% / 100% |
| planted, seed 1976 E (5,067) | 11.0% / 99.8% | 84.6% / 100% | 7 / 12.5 | 96.6% / 100% |
| planted, seed 1976 S (2,418) | 11.1% / 100% | 82.4% / 100% | 7 / 13.1 | 95.0% / 100% |
Clusters at 99%: HG002 40,805 over 42,083 loci (largest 217), KB3781 17,802 over 18,037 (largest 34), planted 41-52
over 134-154 loci (largest 33-34). On HG002, 553 of the 1,459 truth paralogs have an assembled locus under them.
All reads with a family primary, KB3781: per locus 96.0% decided, 96.8% right; 99% clusters 98.5%, 98.0%.

### Claims
- **G1 HOLDS:** HG002 80.0% decided, 96.6% right (row 1974 on annotated genes: 78.4%, 98.0%).
- **G2 HOLDS:** the four planted sets 82.4-84.6% decided, 100% right (one read wrong of 12,850).
- **G3 HOLDS:** KB3781 99.3% right, 74.3% decided.
My prediction that fewer reads would be decided with the run's own loci was wrong for HG002 (80.0% against 78.4%).

### Reading
- **It generalizes.** With clusters made from the run's own copy sequences and no annotation, on two real individuals
  of two species and on planted reads whose sample is not the reference, the cluster-level lead decides 74-85% of the
  reads that are tied between loci and is right in 96.6-100%. The per-locus rule decides 8-16%.
- **The 99% threshold held on data it was not read from** (KB3781, the planted sets); 95% decides a little more on
  planted reads at the same accuracy and changes nothing on the real samples. No reason found to move it.
- **Resolution depends on the family.** Real samples: the median retained cluster is 2 loci. The planted families were
  chosen for being large and nearly identical: there the median is 6-7 loci and the largest cluster 34, so "decided"
  means "one of seven". A verdict must always carry the cluster's size.
- **Error on HG002 is 3.4%** (1.4 points more than on annotated genes): the run's loci do not line up one to one with
  the annotation's genes, and the truth has to be carried from a gene to the locus under it.
- Limits: the truth of HG002 is row 1972's (its limits apply); KB3781 on AS only; planted reads have no chimeras; one
  private-variant rate (0.001); clusters need the `families` stage to have found the copies (a copy it misses is its
  own cluster).
- Not built into the stage. What a build would add: one all-against-all comparison of the copy sequences (seconds),
  the cluster id per locus, and per cluster the reads decided, discarded and the number of loci.

## Built (user 2026-11-26: "let's build it"): `tools/cluster_verdict.py`, stage `o2` with `RUSTLE_O2_CLUSTERS=1`
Opt-in (unset = the stage as before; the default is the user's decision). `RUSTLE_O2_CLUSTER_IDENTITY` (0.99).
Writes `PREFIX.o2.clusters.tsv` (cluster, loci, reads placed alone, tied reads decided, ambiguous reads touching it,
members) and `PREFIX.o2.cluster_loci.tsv` (locus, cluster, loci in the cluster, the same counts). No GTF attribute.
Two differences from what the claims above tested, so the accuracy was measured again on the built step:
- the verdict is read from the stage's tie table (`PREFIX.o2.ties.tsv.gz`: every tied read with its loci and `ms` gaps),
  so there is no extra pass over the reads;
- clusters are made on the loci that tied reads touch, from the sequence of each locus's longest transcript in the
  mode's GTF, because the mode adds loci the `families` stage does not list.

| built step (job 4512159) | tied reads | decided | right | loci in the decided cluster, median / mean |
|---|---|---|---|---|
| HG002 (truth: row 1972) | 7,476 | 83.7% | 97.3% (5,117 of 5,257) | 2 / 3.2 |
| planted, private variants, 1975 E | 5,090 | 82.6% | 100% | 7 / 11.2 |
| planted, 1975 S | 2,796 | 70.3% | 100% | 4 / 6.0 |
| planted, 1976 E | 5,077 | 82.8% | 100% | 7 / 11.2 |
| planted, 1976 S | 2,429 | 70.9% | 99.9% (1 wrong) | 5 / 8.0 |
HG002: 169 clusters over 384 loci, 65 with several loci (largest 49: an array on chr13 with no decided read, and one on
chr16); 5,798 of the 6,257 decided reads are decided for a cluster of several loci, 459 for a single locus. The 140
wrong reads have their true locus in another cluster (NBPF, SRGAP2C: human-specific duplications).
A scoring point found on the way: two assembled loci can lie on one annotated gene; the first scoring took one of them
as the truth and called 202 right reads wrong (93.5%). The scorer now accepts any locus on the truth gene (97.3%).
Stage checks: end to end on the four planted sets (exit 0, tables written); switch unset removes both tables; a value
other than 0 or 1 stops the stage with exit 2; whole `cargo test` 1,591 passed, 0 failed; the tool is in
`tools/test_selftests.py`.
