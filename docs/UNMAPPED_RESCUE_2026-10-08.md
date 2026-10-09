# Rescuing reads that became unmapped: a pool-first, reference-free net (2026-10-08)

> **Independent review (2026-10-08).** Four fresh reviewers audited Amendments 9 to 24. Corrections, limitations and what changed are in the Errata section at the end of `PREREG_unmapped_rescue_2026-10-08.md`; the sentences below that were wrong have been corrected in place.

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
- **Per-family HMM profile (Amendment 5; HMMER 3.4 from the `hmmer` mamba env; H only, exploratory).** Universe: the 53 H families; profile = MAFFT alignment of the survivors' transcript consensus (148 survivor consensus sequences), `nhmmer` against the 47 cluster consensus sequences, cover score with the frozen margin. In the SAME universe the dc-megablast cover score gives 4,279 rescued / 59 wrong (1.36%) / 21 of 24 copies;
  **the profile gives 4,522 rescued / 125 wrong (2.69%) / 22 of 24 copies** (bit-score variant 4,479 / 183 wrong, 3.9%). By the amendment's condition (more rescued than the restricted baseline with wrong joins <= 5%) this is a gain of 243 reads. Where it comes from: two clusters of families with many survivors, the 200-read GWFAM158 cluster (9 survivors) and the 70-read GWFAM173 cluster (5 survivors), both abstained under the cover score and attributed correctly by the profile.
  What it costs: the GWFAM163 profile attracts small clusters of other families (five clusters of 56, 30, 10, 5 and 3 reads were attributed to it wrongly), 122 wrong reads in nine changed clusters, and the 24-read GWFAM440 cluster that the cover score had right now abstains. Among the 19 families with >= 3 survivors the profile rescues 496 reads with 14.3% wrong joins (baseline 226 with 20.7%), small numbers.
  So a profile helps where a family has many surviving members and needs a specificity guard against hub profiles; two clusters carry the whole gain. It is hypothesis-generating (H is not held out); a new multi-copy erasure bed would be needed to confirm it.

## 6a. Hybrid net (Amendment 6): pool first, then the leftover reads one by one

On bed H, 2,014 of the 6,272 pool reads (32.1%) are not rescued by an attributed cluster; attributing them alone (same cover score, margin 1.10) gives 4,421 rescued (83.2%) against 4,255 for the pooled rule and 4,357 for the B1 estimate scaled to the bed, with 71 wrong joins (1.58% of joined). The extra work is 32.1% of per-read attribution of the whole pool (506 s of about 1,580 s), above the registered 25%, so the registered rule is **not met** and the pool-first net stays "better than A13 (557) and about 6 times cheaper than per-read attribution", not "rescues more than per-read attribution". Bed A was not run: the rule needs both beds, so the H failure already fixes the verdict.

## 6b. The reconstructed loci as extra reference (Amendment 7, bed H)

Consensus sequences from a seeded half of the unmapped deleted-copy reads (plus the 959 background reads) are added to the reference by aligning reads to them and keeping the higher score (the primary a concatenated reference would give). On the other, never-clustered half: **2,505 of 2,657 held-out unmapped reads (94.3%) land on the consensus of their own family at de 0.0011**, absorbed deleted-copy reads move only when the consensus is closer (2,250 of 11,973; 2,247 to their own family; median de 0.1237 -> 0.0031), and **0 of 41,727 survivor reads move**, 0 of 959 background reads move onto a family cluster. All three bars met. This is the augmented-reference test of the earlier NPIPA2 result (row 1223 neighbourhood), now with consensus sequences that were built without any reference.

## 6c. Worked example: LRPAP1 and its gorilla copies (Amendment 8)

With the mother as reference, 381 of 2,875 LRPAP1 reads are poorly placed (de > 0.00958). They form 4 pure clusters (283, 21, 4, 3 reads). The 283-read consensus (1,492 bp) is the paternal c01 (0.9987, the chr12 22.55 Mb copy the mother's assembly lacks within the cutoff); its nearest maternal node is c00 at 0.9886, the registered lift site c01|c03 only 0.9839. The 4-read cluster (2,680 bp) is c01 at 0.9993 with the mother's best 0.9888 and is the one flagged at the registered 0.999 (O3). The 21-read cluster is a c00 allele (0.5% from the mother's), the 3-read cluster matches neither haplotype (0.97). Realigning all reads to the 4 consensus sequences: 84.8% of the held-out net reads move to the consensus of their own copy, and 2.8% of the reads of copies without a consensus move. Details and caveats in the prereg, Amendment 8 result.

## 6d. An IsoCon-style correction step on the consensus (Amendment 9): refuted

Aligning each cluster's reads to its consensus, testing every non-consensus allele against the cluster's error rate (Bonferroni 0.05) and replacing the consensus where a significant allele has the plurality changes the consensus in 53 of 221 clusters on dev A and 10 of 47 on H. Held-out H: identity x coverage >= 0.999 in **15 of 40 clusters against 14 before** (bar: at least 19), and 3 clusters get worse (bar: at most 2). **The decision rule is not met; the consensus-error explanation of the fidelity shortfall is refuted.** Most corrections are exon-scale deletions that trim a mixed-isoform consensus to its plurality isoform, not base repair. The residual edit columns are not allelic either (1.2% and 1.9% are significant minority-variant columns), and the registered metric turns out to be dominated by unaligned consensus ends: on dev A identity alone reaches 0.999 in 95 of 212 clusters, identity x coverage in 65. In the LRPAP1 example the 283-read c01 consensus has identity 1.0000 to pat c01 (0 edits over 1,490 aligned bases); the 0.9987 is 2 unaligned terminal bases, so the flag it missed was a coverage artifact of the registered metric, and the polish changes nothing there.

## 6e. A supervised synthetic world: erase one copy, recover it (Amendment 10)

30 synthetic families (random genes, 5 to 10 exons), three copies each: A and B present, E erased from the reference (divergence of E from its nearest survivor 0.5%, 1%, 2%, 4% or 8%, six families each), 80 reads per copy, exact truth. Recovery works wherever the erased copy is outside the allele cutoff: every E copy from 1% up has a pure cluster, its consensus is at least 0.999 identical to the true transcript (core identity, end gaps excluded) in all 26 clusters that exist (4 of the 0.5% copies have none), attribution is 100% right with 0 wrong, and the erased copy's reads land on its own consensus in 100% of cases. At 0.5% (inside the cutoff) only 12% of the erased copy's reads are in the net, as designed. **One registered bar fails: reads of the surviving copies move onto the erased copy's consensus in 30 to 50% of cases when it is within 2% of a survivor** (raw alignment scores of a spliced genome alignment and an unspliced transcript alignment are not comparable); a divergence rule (move only if the divergence drops) removes it (0 of 960 in every class) and leaves bed H unchanged (Amendment 11, exploratory). The ends: the simulated consensus is never longer than the true transcript and the aligner clips nothing, so jitter, 5' truncation and a polyA tail (which kills the metric entirely, 0 of 26) do not reproduce the real beds; the real beds' unaligned consensus ends are an untemplated 5' run of G (1 to 2 bp in 100 of 212 dev-A consensus sequences, none matching the genome next to it).

## 6f. Divergence move rule adopted; the untemplated 5' G run (Amendments 12 to 14)

The divergence rule is the default of the augmented-reference step: on a fresh synthetic world (seed 20261009, 32 families, three divergence classes not seen before) survivor reads moved 0 of 640 in every class (score rule: up to 320), erased-copy reads landed on their own consensus 100% wherever a pure cluster existed, and bed H and LRPAP1 are unchanged or better (LRPAP1 held-out net reads moving 172 of 191, other-copy reads 3 of 532). The detection edge in the synthetic world sits at the cutoff: 0% of an erased copy's reads leave it at 0.5%, 55% at 0.75%, 75% at 1%, all from 1.5%.
The 5' G run: simulated from the measured distribution (95% of reads, length 1 to 6, mode 2), it reproduces the unaligned 5' clip of the real consensus sequences (22 of 25 against 61% on dev A), and a trimming rule that measures the prefix against the family's surviving copies (T2) restores the registered metric in 21 of 23 synthetic families with no over-trimming and no false trims. **It is not adopted:** on real families the survivors are too divergent to anchor the 5' end (any alignment for only 44 of 221 consensus sequences), so T2 trims 3 of 212 on dev A and none of the 60 attributed clusters that carry the artifact. The artifact will need a library-level signature (the 5' soft-clipped pure-G run on cleanly aligned reads of surviving copies) rather than a per-family reference.

## 6g. The 5' G run read from the library (Amendment 15): works on synthetic, not adopted on real

73% of the cleanly aligned survivor reads of the real library have a 5' soft clip that is a pure G run of 1 to 3 bases (A, T or C: 80 reads in all): the artifact is a property of the library. A gate that tests this (G over-represented among pure 5' clips, binomial p < 1e-6) is open on the real library and on a synthetic world with the artifact and closed on one without. Trimming the leading G run (1 to 3 bases) of every consensus when the gate is open restores the registered identity x coverage metric on the synthetic world (12 to 22 of 24 erased copies, no harm) and raises it on the real beds (dev A 65 to 82 of 212, H 14 to 18 of 40). **It is not adopted:** it fails two of the three registered real-data bars (recall 92% against 95%; templated sequence lost at most 1 base in 67% against 90%), because real 5' ends are GC-rich: 61 of 197 trimmed consensus sequences had no unaligned G at all (105 had no pure-G clip of 1 to 3; the independent review corrected the accounting: templated bases lost at most 1 in 157 of 197 = 80%, still under the 90% bar), so their leading Gs were aligned to the erased copy and the alignment cannot say which are templated and which are artifact Gs that matched the G-rich flank. Practical consequence: the registered 0.999 flag is penalised by this artifact in about half of the clusters, and either the flag should ignore a leading G run of at most 3 bases when the library gate is open (score side, consensus untouched) or the consensus is trimmed and 1 to 3 templated bases are risked.

## 6h. Per-cluster correction of the 5' G run (Amendment 16, rule T4): best of the options, not adopted

Each read's own leading G run is r = j + a: j the templated G's (shared by the cluster) and a the artifact length (a library distribution, read from the aligned survivors' soft clips). The maximum-likelihood j per cluster, then removing only the consensus' leading run beyond j. On a synthetic world with 0 to 3 templated leading G's: the consensus is within one base of the right 5' end in 24 of 24 clusters (T3, trimming everything: 16 of 24), removes templated bases in 1 of 24 (T3: 7), and restores identity x coverage >= 0.999 in 24 of 24 (T3: 20, untouched: 16); the control without the artifact is untouched. Scored against the MODAL read start it fails the registered exactness bars (58% exact, artifact left in 38%); the independent review showed that convention is wrong when the jitter is one base (half the reads start one base in, so a consensus that begins exactly at the true transcript scores +1; the untouched control scores 21/24 exact). With a start-tolerant measure (the consensus begins at ANY observed read start) T4 is exact in 23 of 24 clusters, removes templated bases in 0 and leaves an artifact G in 1 on this world and on a fresh one (seed 20261301: 23, 0, 1; T3 15, 1, 8), i.e. it meets the synthetic exactness bars. It is still not adopted because the dev-A bar fails (below). On real dev A: 65 to 80 of 212 at 0.999 (bar met), but the registered bar that it removes at least the unaligned clip fails (64 of 100): 23 of the 36 shortfalls are clusters of fewer than 5 reads, where the rule does nothing, and 13 are short by one base.

## 6i. The gate-aware flag metric and the IsoCon-style partition (Amendments 17 to 20)

**Gate-aware identity x coverage (adopted, score side).** When the library gate is open, an unaligned leading G run of at most 3 bases leaves the denominator. Synthetic: 12 to 23 and 16 to 24 of 24 at 0.999, controls identical, 0 false passes; real: dev A 65 to 80, H 14 to 17; LRPAP1: the 283-read c01 consensus is 1.0000 to the paternal c01 and is now flagged (mat best 0.9899), the c00 allele cluster is flagged too (the flag never separated an allele from an absent copy).
**Partition into competing candidates (adopted as an opt-in step).** For a cluster: significant variant columns (the Amendment 9 test), pairs of variants that travel together across reads (hypergeometric, Bonferroni), linked blocks (a lone hotspot column is not one), split on the block that divides the reads most evenly, recurse. On a fresh synthetic world (36 families: control, an erased copy with a second isoform skipping an exon, two sibling erased copies 0.5% apart) it splits 0 of 12 controls, separates all 12 isoform clusters into pure leaves (baseline recovers 12 of 24 transcripts, PART 24 of 24) and has no spurious or redundant candidate; the frozen clustering already separated the sibling pairs in that world, while in two earlier worlds (which used the rejected `map-hifi` and edlib alignments) where it joined 5 and 1 pairs PART recovered 24 of 24 (from 19 and 23); a forced-join test with the adopted alignment by the independent review separated 0.5% sibling pairs in 12 of 12 but only 10 of 12 at 0.2%. Two implementation lessons: `minimap2 map-hifi` soft-clips a long deletion (use `splice:hq -uf` and count N as a deletion), and a scorer that takes the best of two orientations by core identity can score 1.0 on garbage (fixed; the first Amendment 18 to 20 verdicts were partly an artifact of it, see the prereg correction). On the real beds PART splits 72 of 212 (dev A) and 19 of 40 (H) clusters and, taking the best candidate per cluster, raises the gate-aware 0.999 count from 80 to 106 and from 17 to 26; it does not split the 22 clusters with more than 200 edits, so they remain unexplained.

## 6j. The 22 dev-A clusters with more than 200 edits (Amendment 21)

They are 3 to 12 reads each, and the reads are fine one by one (baseline divergence 0.0013, MAPQ 60). The consensus is what is wrong: in 16 of the 22 it does not explain its own reads (median read-to-consensus divergence 0.016 against 0.0015 elsewhere). The reads of such a cluster differ 4 to 10 fold in length and start tens of kb apart in the gene (667, 3,657 and 4,508 bases in one 3-read cluster): the edge rule lets a short read bridge two long reads that do not overlap, the component is a chain of partial or alternative-structure reads, and a global POA consensus of a chain is a blend. A truth-free flag (median read divergence to the consensus above the allele cutoff 0.00958) finds 16 of 22 on dev A and 2 of 3 on H, with 4% and 3% of the other clusters flagged; it marks the consensus unsupported and changes nothing. Six clusters on dev A (one on H) are not explained: their consensus fits their reads but differs from the erased copy by 238 to 614 edits. In reads they are 0.2% of the pool; in clusters 10%, and they cap the cluster-level fidelity count (80 or 106 of 212).

## 6k. Chain-aware clustering (Amendments 22 to 24): adopted, opt-in

The chain was a short read CONTAINED in two long reads that disagree: each containment is a proper overlap, and an alignment of two isoforms that differ by an exon-scale deletion has a tiny gap-compressed divergence, so single linkage joins them. Two changes remove it: the edge must be a proper overlap (a suffix-prefix overlap or a containment, not a shared middle), and, inside each component of at most 60 reads, a longest-first star clustering in which a read joins a cluster only through its representative and the pair divergence counts every gap (NM / block). Synthetic chain world (fresh seed, 12 control and 12 two-isoform families at low depth with fragments): isoform transcripts recovered 24 of 24 against 12 of 24, purity 0.91, control untouched. Real beds (corrected after the independent review found that the first star step used minimap2's default all-vs-all limits, which drop edges in components of 20 or more reads; fixed, tested, and re-run): clusters with more than 200 edits 22 to 13 (dev A) and 3 to 0 (H), unsupported consensus 24 to 14 and 3 to 0, erased-copy clusters at 0.999 38% to 52% (152 of 294) and 43% to 56% (31 of 55), at a cost of 1.4 and 2.7 points of read coverage; 71 of the 126 reads of the 22 dev-A clusters (and 14 of 14 on H) leave the clustering, so part of the improvement is removal. The star step touches only components of at most 60 reads (4 to 8% of the reads); larger clusters, where the mixed-isoform consensus problem of Amendment 9 lives, are untouched. Fresh worlds with the fixed step (depth 16, 22 and 50 reads per copy): isoform transcripts recovered 24 of 24 each, controls untouched, purity 0.89 to 0.93. The proper-overlap edge alone is not enough (Amendment 22 failed its bars) and neither is the star step with the gap-compressed divergence (Amendment 23).

## 6l. Specificity controls (Amendment 25)

**PART guard.** A hotspot (a recurrent sequencing error at fixed columns in a quarter of the reads) made PART split 10 to 12 of 12 single-transcript controls; a block-purity test (reads must carry all or none of a block's variants, and the non-carriers must not carry them above the background error rate) brings that to 1, 0 and 0 of 12 on two fresh seeds, keeps exon-skip and 0.5% sibling splits at 12 of 12, and leaves the fresh partition world unchanged (control 0 of 12, isoforms 24 of 24). On the real beds the many-leaf splits disappear (at most 7 leaves against 14) and the split share falls from 33% to 29% (dev A) and 43% to 36% (H), which is still not evidence of distinct transcripts: alleles are not excluded.
**Allele vs copy in the O3 call.** The registered flag (reconstructed on the truth genome, not identical on the reference) flags an allele of a present copy. Requiring a COPY to be more than the allele cutoff (0.958%) from every reference locus separates them: in a synthetic world 12 NULL families give no cluster, alleles at 0.25% and 0.5% give no COPY call (one ALLELE), erased copies are COPY 12 of 12, and on LRPAP1 the c00 allele is ALLELE while both c01 clusters are COPY (the c01 margin is 0.0005). Alleles 1% or more apart are called COPY (4 of 6 at 1%, 6 of 6 at 2%): RNA alone cannot separate them from paralogs.

## 7. What did not go as predicted, and caveats

1. Pooling does not rescue more reads than aligning each read with the same sensitive score; it is cheaper and gives the consensus.
2. The registered translated attribution (mmseqs) was wrong at 10 to 17% of joined reads; the nucleotide cover score with a margin replaced it on the dev bed (Amendment 1). The margin 1.10 is a dev-bed choice that, on the corrected data, sits just above a few wrong clusters on both beds.
3. Both erasure beds are synthetic (a copy erased from the same animal's genome) and use catalog families as the target set; the real bed M has no family truth, only the maternal alignment.
4. Bed H was run once with the frozen rule (apart from a reporting crash that was fixed and re-run with identical numbers); the code that implements the rule was committed after that run (Amendment 4).
5. The unmapped pile of a novel gene (bed M) is clustered and reproduced, not named.

## 8. Not covered

Poorly placed reads (absorbed on a paralog rather than unmapped); the full 23 GB library as one pool; IG/TR screening of flagged piles beyond the run-composition check; independent replication; a multi-copy erasure bed to confirm the HMM profile result (only H has surviving-member profiles); any change to the shipped binary.
