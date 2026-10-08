# Rescuing reads that became unmapped: a pool-first, reference-free net (2026-10-08)

Goal (user): design a better net or clustering that saves unmapped reads, from a parental haplotype or after erasing a copy, realigning and rescuing the lost reads.
Pre-registration `docs/PREREG_unmapped_rescue_2026-10-08.md` (body + Amendments 1 to 4); code `bench/unmapped_rescue/` (unit-tested); work dir `/mnt/linuxdisk/tmp/o3_rescue/`.
This write-up was rewritten after an independent review that recomputed every number: the held-out numbers reproduced, the dev-bed numbers changed because a BLAST target cap truncated hits (Amendment 4), and several statements were corrected.

## 1. What this is, and what is O3 and what is not

- **O3 (detect and flag a reference-absent copy)** needs a cluster of unmapped reads and its consensus; it does not need a family. Results 4.1 (clustering) and 4.4 (consensus fidelity, real bed) are the O3 results.
- **Family attribution** (sections 4.2, 4.3) puts the reads of a flagged cluster into a known family's net, so that a family-scoped stage can run. It is an O1-style membership question, not O2 (read-to-copy assignment by PSVs), which is not done here.

## 2. The method

1. **Pool first.** Only the unmapped reads, no reference. Each round maps every unclustered read to 2,000 seed reads (`minimap2 -x map-hifi -c`); an edge needs gap-compressed divergence `de` <= 0.00958 (the registered allele cutoff) and an alignment block >= 50% of the shorter read;
   components of >= 3 reads are clusters; the rest go to a new round with new seeds (at most 6 rounds). This replaces the all-vs-all of the prereg body, which did not finish: one shard of 10,400 reads ran past 10 minutes because the 7,000-read families make it quadratic.
2. **Consensus.** One abPOA consensus per cluster (<= 100 longest reads).
3. **Attribute the cluster (optional for O3).** Discontiguous megablast of the consensus against the surviving copies of the family catalog (erased copies removed; `-max_target_seqs 5000`); a family's score is the consensus bases covered by the HSPs of its best copy;
   the cluster goes to the best family iff its score is >= 1.10 x the runner-up's, else it abstains.
4. **Rescue.** Every read of an attributed cluster joins that family's net.

## 3. Beds

| bed | what | pool |
|---|---|---|
| A (dev) | 2026-08-14 whole-genome excision, 162 two-copy families erased from the full genome and the reads realigned | 59,791 deleted-copy reads left unmapped + 959 background reads (plus 1,719 reads without a label, mostly deleted-copy reads the start-only truth rule misses) |
| H (held out for the frozen rule) | A13's 53 held-out multi-copy families, masked genome | 5,313 deleted-copy reads left unmapped + the same 959; no survivor or unlabelled reads at all |
| M (real) | the 959 reads unmapped on the combined primary assembly of KB3781 against its maternal haplotype (124 or 125 map in one 71-kb stretch of mat chr12) | in A's pool, so M's clustering was dev-visible |

## 4. Results

### 4.1 Clustering (O3; U1 on bed H)

| bed, delta | deleted-copy reads in a cluster of >= 3 | purity | clusters |
|---|---|---|---|
| A, delta/2 | 53,320 of 59,791 (89.2%) | 1.0000 | 240 |
| A, delta (registered) | 57,664 (96.4%) | 1.0000 | 221 (1 to 7 per family, median 2) |
| A, 2 delta | 59,087 (98.8%) | 1.0000 | 203 |
| **H, delta (registered)** | **5,083 of 5,313 (95.7%)** | **1.0000** | 47 (1 to 4 per family); 0 background reads in a family cluster |

Purity is the share of clustered deleted-copy reads whose cluster's majority family is their own. **U1 passes** (bar: coverage >= 0.90 and purity >= 0.99). A has 1 background read in a family cluster.

### 4.2 Rescue into a family's net with the frozen rule (U2; frozen before bed H was run)

| bed, rule | rescued-correct (share of unmapped deleted-copy reads) | wrong joins (share of joined) | clusters attributed / abstained | cluster accuracy | copies reached |
|---|---|---|---|---|---|
| **H, frozen: cover score, margin 1.10** | 4,255 (80.1%) | 3 (0.07%) | 27 / 20 | 0.963 | 20 of 32 |
| H, strict (margin 1.00) | 4,255 (80.1%) | 474 (10.02%) | 38 / 9 | 0.684 | 20 of 32 |
| H, margin 1.50 | 4,243 (79.9%) | 3 (0.07%) | 25 / 22 | 0.960 | 19 of 32 |
| H, translated mmseqs, strict | 4,257 (80.1%) | 498 (10.47%) | 38 / 9 | 0.658 | 21 of 32 |
| A (dev), frozen | 55,901 (93.5%) | 53 (0.09%) | 126 / 95 | 0.937 | 76 of 118 |
| A (dev), strict | 55,912 (93.5%) | 626 (1.11%) | 179 / 42 | 0.676 | 78 of 118 |
| A (dev), margin 1.50 | 55,861 (93.4%) | 10 (0.02%) | 113 / 108 | 0.973 | 75 of 118 |
| A (dev), translated mmseqs, strict | 46,545 (77.8%) | 9222 (16.54%) | 175 / 46 | 0.589 | 72 of 118 |

- A13's net on bed H joined 557 unmapped reads (10.5%, 99.6% right; register row 1223). Quote the held-out precision as **1 wrong cluster of 27 (cluster accuracy 0.963)**; "3 wrong joins" is that one small cluster, read-weighted.
- The unthresholded top-1 family of the cover score is the true family for 26 of 38 clusters on H (68%; 4,255 of 4,669 reads) and 121 of 185 on A (65%; 55,912 of 56,521 reads): the margin rule removes the 12 and 64 wrong top-1 clusters, which are small.
- **The margin is fragile at 1.10.** Rescued / wrong joins by margin (A then H): A: 1.0: 55,912 / 626 | 1.05: 55,912 / 237 | 1.09: 55,903 / 219 | 1.1: 55,901 / 53 | 1.2: 55,886 / 46 | 1.5: 55,861 / 10 | 2.0: 55,466 / 7; H: 1.0: 4,255 / 474 | 1.05: 4,255 / 138 | 1.09: 4,255 / 129 | 1.1: 4,255 / 3 | 1.2: 4,252 / 3 | 1.5: 4,243 / 3 | 2.0: 4,052 / 3.
  On A the drop from 1.09 to 1.10 is one 166-read cluster (GWFAM316 attributed to GWFAM317 at ratio 1.092); on H it is the 116-read M pile and a 10-read cluster (ratios 1.097 and 1.095). Above 1.10 the wrong joins stay low, but going from 1.10 to 1.50 on A costs 40 reads and 1 copy and removes 43 wrong joins: 1.10 sits just above a handful of wrong clusters, it is not a plateau of identical answers.
- With the cap fixed, the strict rule is wrong at 1.1% of joined reads on A but still 10.0% on H; Amendment 1's account of why a margin is needed was overstated for A.
- Reads of foreign origin joining an attributed cluster are not judged: A has 1,620 (reads without a truth label, mostly deleted-copy reads mislabelled by the start-only truth rule, so A's rescued count is understated by about that many); H has none, so wrong joins for foreign reads are untested there.

**Why 20 to 25% of H's reads stay unrescued** (abstentions by where the true family ranks, H; reads in brackets): no nucleotide hit 8 clusters (414), true family ranked below 2nd 6 (228), true family 2nd 2 (94), true family has no hit while another does 3 (89), 1 background cluster; plus 230 reads in no cluster.
On A: no hit 33 (1,143), ranked below 2nd 43 (283), true family 2nd 1 (157), no hit for the true family 13 (116), true family 1st but margin small 3 (11). The biggest dev cluster that abstained under the capped output (GWFAM133, 3,408 reads) is attributed correctly once the cap is fixed.
The earlier statement that abstentions are mostly near-ties between two close families is wrong: for 9 of the 12 margin or tie abstentions on H the true family is not in the top two.

### 4.3 Comparators on the same reads (500 seeded unmapped deleted-copy reads per bed; rescued-correct / wrong / abstained)

| bed | B0: map-ont of single reads to the surviving genomic copies, with A13's coverage and divergence floors | B1: single-read dc-megablast, same cover score and margin | pooled, frozen rule |
|---|---|---|---|
| H | 0 / 0 / 500 | 410 / 4 / 86 | 406 / 0 / 94 |
| A (dev) | 0 / 0 / 500 | 475 / 1 / 24 | 469 / 1 / 30 |

- B0 is NOT A13's rule: A13 also aligned the reads to the family's own net reads and joined 557 reads; B0 aligns single reads to the intron-containing genomic copies and finds almost nothing to align (the usable minimap2 records are 33 and 47 reads of 500), so it measures minimap2's inability to align 15 to 30% divergent spliced reads to genomic copies. B1 > B0 is therefore vacuous.
- **U2's first clause passes on H** (rescued-correct 4,255 >= 2,656, wrong joins 3 = 0.07% <= 5%); **its second clause fails**: the pooled method does not beat single-read dc-megablast (406 against 410 on H, 469 against 475 on A).
  So the rescue comes from the sensitive nucleotide score and the margin rule, not from pooling. Pooling is about six times cheaper in wall-clock here (13 minutes against an estimated 80 for per-read alignment, scaled from the sample).

### 4.4 Consensus fidelity and the real bed (O3)

- **Fidelity (registered metric: identity x coverage >= 0.999 against the unmasked genome).** Over all clusters with a truth family: 14 of 40 on H and 65 of 212 on A; over attributed clusters: 12 of 27 (H) and 39 of 125 (A).
  By identity alone: 20 of 40 and 95 of 212. That every such cluster's best hit overlaps its erased interval (40 of 40, 212 of 212) follows from how the truth family is defined and is not a result.
- **Bed M.** 115 of the 125 reads that map to mat chr12 95.96-96.03 Mb form one cluster (92%; the other 10 are in no cluster) with no other read in it (the held-out bed's pool gives 116); its consensus is 5,694 bp and aligns to the maternal region at identity 1.0000 (5,557 of 5,557 matches, query coverage 0.976).
  **Attribution of this pile is not stable:** against the full 915-copy catalog its best family is GWFAM267 at 407 covered bases against 371 for GWFAM384 (ratio 1.097, just below the margin, so the frozen rule abstains); with the BLAST cap of 500 the same pile was attributed to GWFAM267 (reviewer's recomputation), and against bed A's reduced target set the call was 371 against 369.
  The best hit comes from a short element that hits 831 of the 915 copies (an independent recomputation puts it near 320 bp, 7% of the consensus). There is no family-specific evidence; treat it as unattributed, which for O3 is a complete flag.
- **Run composition.** The 115 reads come from all four SRA runs (8, 16, 82 and 9): not run-exclusive. The dominant run supplies 71% of the pile, close to its 72% share of the 124 chr12 reads that map, and above its 48% share of the 959 background reads.
- **Cost.** Clustering and consensus of the 62,471-read pool took about 10 minutes and the 221 consensus sequences about 2 minutes of alignment; aligning each read alone took 40 s per 500 reads on H and 69 s on A.

## 5. Verdicts against the registered bars and predictions

| bar | verdict |
|---|---|
| U1 (H): coverage >= 0.90, purity >= 0.99 | **PASS** (0.957, 1.0000) |
| U2 (H): rescued-correct >= 2,656 with wrong joins <= 5%, AND pooled > B1 > B0 | **FAIL as a conjunction**: first clause PASS, pooled does not beat B1 (406 vs 410; 469 vs 475 on A), B0 is not A13's rule |
| U3 (M): >= 90% of the 125 reads in one cluster, consensus identity >= 0.99, no other read in it | **PASS** (92%, 1.0000, none), clustering part dev-visible; the family attribution of the pile is unstable (4.4) |
| prediction: U1 holds | held |
| prediction: pooled > B1 > B0, rescued 20 to 60% of N_D | **failed**: pooled is slightly below B1, rescued 80.1% on H |

## 6. Exploratory arms after the held-out run (Amendments 2 and 3; bed H is no longer held out for them)

- **Family union** (the bases covered by any surviving member, a profile proxy): identical to the frozen rule on H and within 3 reads on A.
- **Oracle homology check** of the 13 abstained H clusters with a truth family against only their own family's survivors (blastn word 7, tblastx, permissive e-values): clear homology for 11 (for example 475 bits, E 5e-135, for the 200-read GWFAM158 cluster) against 30 to 150 bits for random families; none for GWFAM244 and one GWFAM348 cluster. The signal exists and the cover score ranks it badly.
- **Repeat-robust specific cover** (a consensus base covered by n families credits 1/n to each; Amendment 2) and with an evidence floor (Amendment 3, f selected on A as the smallest with wrong <= 1%: f = 0.02): A 55,920 rescued, 270 wrong (0.48%); **H 4,325 rescued, 271 wrong (5.90% > 5%): fails on H.** Without the floor H is 457 wrong (9.6%). A floor of 0.05 (not the selected one) gives H 4,321 / 27 wrong and A 55,852 / 14, reported and not claimed.

## 7. What did not go as predicted, and caveats

1. Pooling does not rescue more reads than aligning each read with the same sensitive score; it is cheaper and gives the consensus.
2. The registered translated attribution (mmseqs) was wrong at 10 to 17% of joined reads; the nucleotide cover score with a margin replaced it on the dev bed (Amendment 1). The margin 1.10 is a dev-bed choice that, on the corrected data, sits just above a few wrong clusters on both beds.
3. Both erasure beds are synthetic (a copy erased from the same animal's genome) and use catalog families as the target set; the real bed M has no family truth, only the maternal alignment.
4. Bed H was run once with the frozen rule (apart from a reporting crash that was fixed and re-run with identical numbers); the code that implements the rule was committed after that run (Amendment 4).
5. The unmapped pile of a novel gene (bed M) is clustered and reproduced, not named.

## 8. Not covered

Poorly placed reads (absorbed on a paralog rather than unmapped); the full 23 GB library as one pool; IG/TR screening of flagged piles beyond the run-composition check; independent replication; the per-family HMM profile arm (HMMER is installed in the `hmmer` mamba environment; not run in this version); any change to the shipped binary.
