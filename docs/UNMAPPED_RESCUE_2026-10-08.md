# Rescuing reads that became unmapped: a pool-first, reference-free net (2026-10-08)

Goal (user): design a better net or clustering that saves unmapped reads, from a parental haplotype or after erasing a copy, realigning and rescuing the lost reads.
Pre-registration `docs/PREREG_unmapped_rescue_2026-10-08.md` (body + Amendment 1, committed before beds H and M were run); code `bench/unmapped_rescue/` (unit-tested); work dir `/mnt/linuxdisk/tmp/o3_rescue/`.
Every number below is copied from `W/<bed>/<tag>/{metrics,rescue,fidelity}.json` and `W/<bed>/baselines.json`.

## 1. The method

1. **Pool first.** Only the unmapped reads, no reference. Each round maps every unclustered read to 2,000 seed reads (`minimap2 -x map-hifi -c`); an edge needs gap-compressed divergence `de` ≤ 0.00958 (the registered allele cutoff) and an alignment block ≥ 50% of the shorter read;
   components of ≥ 3 reads are clusters; the rest go to a new round with new seeds. (An all-vs-all of the pool is quadratic in the 7,000-read families: a first all-vs-all shard did not finish in 10 minutes, so the pool is clustered with seeds.)
2. **Consensus.** One abPOA consensus per cluster (≤ 100 longest reads).
3. **Attribute the cluster.** Discontiguous megablast of the consensus against the surviving copies of the family catalog (the erased copies removed); a family's score is the number of consensus bases covered by the HSPs of its best copy;
   the cluster goes to the best family iff its score is ≥ 1.10 × the runner-up's (a single-family hit counts), else it abstains. 1.10 was chosen on the dev bed, where 1.10 to 1.50 give the same answer.
4. **Rescue.** Every read of an attributed cluster joins that family's net.

## 2. Beds

| bed | what | pool |
|---|---|---|
| A (dev) | 2026-08-14 whole-genome excision, 162 two-copy families erased from the full genome and the reads realigned | 59,791 deleted-copy reads left unmapped + 959 background reads |
| H (held out) | A13's 53 held-out multi-copy families, masked genome | 5,313 deleted-copy reads left unmapped + the same 959 |
| M (real) | the 959 reads unmapped on the combined primary assembly of KB3781 against its maternal haplotype (125 map in one 71-kb stretch of mat chr12) | in the A pool (so M's clustering was dev-visible; see Amendment 1) |

## 3. Results

### Clustering (U1, held-out bed H; dev sweep on bed A)

| bed, δ | deleted-copy reads in a cluster of ≥ 3 | purity | clusters |
|---|---|---|---|
| A, δ/2 | 53,320 of 59,791 (89.2%) | 1.0000 | 240 |
| A, δ (registered) | 57,664 (96.4%) | 1.0000 | 221 (1 to 7 per family, median 2) |
| A, 2δ | 59,087 (98.8%) | 1.0000 | 203 |
| **H, δ (registered)** | **5,083 of 5,313 (95.7%)** | **1.0000** | 47 (1 to 4 per family); 0 background reads in a family cluster |

Purity is the share of clustered deleted-copy reads whose cluster's majority family is their own. **U1 passes** (bar: coverage ≥ 0.90 and purity ≥ 0.99).

### Rescue with the frozen attribution rule (U2)

| bed, rule | rescued-correct (share of unmapped deleted-copy reads) | wrong joins (share of joined) | clusters attributed / abstained | cluster accuracy | copies reached |
|---|---|---|---|---|---|
| **H, frozen: cover score, margin 1.10** | 4255 (80.1%) | 3 (0.07%) | 27 / 20 | 0.963 | 20 of 32 |
| H, strict (margin 1.00) | 4255 (80.1%) | 474 (10.02%) | 38 / 9 | 0.684 | 20 of 32 |
| H, margin 1.50 | 4243 (79.9%) | 3 (0.07%) | 25 / 22 | 0.960 | 19 of 32 |
| H, translated mmseqs, strict | 4257 (80.1%) | 498 (10.47%) | 38 / 9 | 0.658 | 21 of 32 |
| A (dev), frozen | 52496 (87.8%) | 57 (0.11%) | 125 / 96 | 0.928 | 75 of 118 |
| A (dev), strict | 52505 (87.8%) | 4032 (7.13%) | 179 / 42 | 0.659 | 77 of 118 |
| A (dev), translated mmseqs, strict | 46545 (77.8%) | 9222 (16.54%) | 175 / 46 | 0.589 | 72 of 118 |

For comparison, A13's net on bed H joined 557 unmapped reads (10.5%, 99.6% right; register row 1223).
The abstentions are legitimate: on H 11 clusters (355 reads) fail the margin because two families cover the consensus almost equally (for example 1,667 against 1,665 bases for GWFAM158), 8 clusters (414 reads) have no nucleotide hit at all (GWFAM244's 409 reads), 1 ties; a further 230 reads sit in no cluster (< 3 reads).
The strict rule and the translated score show the failure the margin removes: near-ties between paralogous families give 7 to 17% wrong joins.

### Comparators on the same reads (500 seeded unmapped deleted-copy reads per bed; rescued-correct / wrong / abstained)

| bed | B0: A13's floors on single reads (map-ont, coverage ≥ 0.5, de ≤ 0.20) | B1: single-read dc-megablast, same cover score and margin | pooled, frozen rule |
|---|---|---|---|
| H | 0 / 0 / 500 | 410 / 4 / 86 | 406 / 0 / 94 |
| A (dev) | 0 / 0 / 500 | 464 / 1 / 35 | 445 / 1 / 54 |

**U2's first clause passes on H** (rescued-correct 4,255 ≥ 2,656, wrong joins 3 = 0.07% ≤ 5%); **its second clause fails**: the pooled method does not beat single-read dc-megablast (406 against 410 reads on H, 445 against 464 on A).
So the rescue comes from the sensitive nucleotide score and the margin rule, not from pooling. What pooling adds is cost and structure, below.

### Candidate transcripts and the real bed (U3)

- **Fidelity.** The consensus of every cluster with a family lands on its erased copy in the unmasked genome: 40 of 40 on H and 212 of 212 on A. 20 of 40 (H) and 95 of 212 (A) reproduce it at identity ≥ 0.999: the erased copy's transcript is rebuilt from the unmapped reads alone.
- **Bed M.** 115 of the 125 reads that map to mat chr12 95.96-96.03 Mb form one cluster (92%; the other 10 are in no cluster) with no other read in it; its consensus is 5,694 bp and aligns to the maternal region at identity 1.0000 (5,557 of 5,557 matches, query coverage 0.976).
  The cluster is not attributed to a family (best family 371 against 369 covered bases): a flagged candidate, as the register requires. Its reads come from all four SRA runs in the proportions of the library (run-exclusivity screen passes: 82 of 115 from the dominant run, which carries 82% of the library). **U3 passes**, but its clustering part was visible in the dev run.
- **Cost.** Clustering and consensus of the 62,471-read pool took about 10 minutes and the 221 consensus sequences one 2-minute alignment (about 13 minutes in all); aligning each read alone took 40 s per 500 reads, about 80 minutes for the 59,791 deleted-copy reads (an estimate scaled from the sample), so pooling is roughly six times cheaper here.

## 4. Verdicts against the registered bars and predictions

| bar | verdict |
|---|---|
| U1 (H): coverage ≥ 0.90, purity ≥ 0.99 | **PASS** (0.957, 1.0000) |
| U2 (H): rescued-correct ≥ 2,656 with wrong joins ≤ 5%, AND pooled > B1 > B0 | **FAIL as a conjunction**: first clause PASS, pooled ≯ B1 (406 vs 410; 445 vs 464 on A) |
| U3 (M): ≥ 90% of the 125 reads in one cluster, consensus identity ≥ 0.99, no other read in it | **PASS** (92%, 1.0000, none), clustering part dev-visible |
| prediction: U1 holds | held |
| prediction: pooled > B1 > B0, rescued 20-60% of N_D | **failed**: pooled ≈ B1 (single-read), rescued 80.1% on H (above the range) |

## 5. What did not go as predicted, and caveats

1. Pooling does not rescue more reads than aligning each read with the same sensitive score (it rescues about as many, slightly fewer, with fewer wrong joins, in about 13 minutes of wall-clock against an estimated 80 for per-read alignment, roughly six times less).
2. The registered translated attribution (mmseqs) was wrong at 10 to 17% of joined reads; the nucleotide cover score with a margin replaced it on the dev bed (Amendment 1). The margin 1.10 is a dev-bed choice; H and M did not influence it.
3. Both erasure beds are synthetic: a copy erased from the genome of the same animal. The real bed M has no family truth, only the maternal alignment.
4. About 20% of the unmapped deleted-copy reads stay unrescued on H: reads of copies whose consensus has no nucleotide relative among the surviving copies (no hit), reads of copies tied between two close families, and reads in no cluster. Translated search and per-read attribution of singletons are untested extensions.
5. M's clustering part and bed A's background reads were visible in the dev run (stated in Amendment 1). Bed H was run once with the frozen rule, apart from a reporting crash that was fixed and re-run with identical numbers.
6. Families in the targets are the catalog families; a real library also has genes outside the catalog, which this evaluation does not test (an unmapped pile of a novel gene is the bed M case: clustered and reproduced, not named).

## 6. Not covered

Poorly placed reads (absorbed on a paralog rather than unmapped); the full 23 GB library as one pool (the unmapped pool there is about 960 reads, the problem is the erasure beds' 62,000); IG/TR screening of flagged piles beyond the run-exclusivity check; independent replication; any change to the shipped binary.
