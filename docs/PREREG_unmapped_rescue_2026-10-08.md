# Pre-registration — rescuing reads that became unmapped: a pool-first, reference-free net

**Written 2026-10-08, before any clustering or attribution result exists for this question.** Goal (user, 2026-10-08): design a better net or clustering that saves unmapped reads, either
those that are unmapped on one parental haplotype or those left unmapped after a copy is erased from the genome and the reads are realigned.
Standing context: `docs/archive/2026-10/O3_CANDIDATES_ACCEPTANCE_A13_2026-10-03.md` (rows 1217, 1221, 1223): after the A13 net (attribution of EACH unmapped read by alignment to the family's reads and copies, read coverage >= 0.5, de <= 0.20)
23 of 53 deleted copies still had no read in their family's net and 7,184 deleted-copy reads were out of reach; 31-mers placed 0 of 5,312 unmapped reads. `docs/NEGATIVE_RESULTS_REGISTER.md` row 1065: a pile of unmapped reads is emitted as a FLAGGED CANDIDATE, never a family;
IG/TR and run-exclusivity screens are compulsory for anything called a new copy. This study rescues READS into a net; it does not define families.

## 1. The idea (frozen)

The A13 net judges every unmapped read alone against the family. Reads of one deleted copy are near-identical to EACH OTHER (one transcript, HiFi errors) even when each is 15-30% from every surviving relative. So:
1. **Pool first.** All-vs-all overlap of the unmapped pool alone (no reference). A read-read edge exists iff the overlap's gap-compressed divergence `de` <= delta = 0.00958 (the registered allele cutoff) and the aligned block covers >= 50% of the shorter read (the registered merge rule of Amendment 8).
   Clusters = connected components; a cluster needs >= 3 reads (the stage's `--min-cluster 3`).
2. **Consensus.** One abPOA consensus per cluster (<= 100 reads, longest first, seed 1).
3. **Attribute the cluster, not the read.** Translated search of the consensus (mmseqs2, `--search-type 2`, both strands, e-value 1e-5 as the external BLAST convention) against the surviving copies of the family catalog (915 copies minus the erased ones, DNA of the unmasked `_pri` intervals).
   The cluster goes to the family of its best hit iff that family's best bit score is strictly above every other family's (no tie); otherwise the cluster ABSTAINS.
4. **Rescue.** Every read of an attributed cluster joins that family's net.
Comparators: B0 = A13's single-read rule re-run on the same pool (minimap2 map-ont to the surviving copies, read coverage >= 0.5, de <= 0.20); B1 = single-read translated attribution (the same search as step 3, one read at a time, same tie rule). B1 isolates what pooling and consensus add.

## 2. Beds and what is held back

| bed | substrate | pool | truth | role |
|---|---|---|---|---|
| **A** (dev) | 2026-08-14 whole-genome excision, 162 two-copy families, `o3_excise/masked.bam` (347,798 reads, masked genome) | the 61,512 unmapped records + the 959 unmapped primaries of the fibroblast BAM (background) | a read is D of family F iff its baseline primary (`baseline.tsv`) lies in the erased interval of F | all design checks and the sensitivity sweep are on A |
| **H** (held out) | A13's 53 held-out multi-copy families, `rna_allele/linktest/R.bam` (masked), `labels.tsv` role D | its 5,313 unmapped records + the same 959 background reads | `labels.tsv` (copy of origin) | **not looked at until the rules in section 3 are frozen by a commit**; run ONCE |
| **M** (real) | the 959 unmapped primaries of the KB3781 fibroblast BAM against the maternal haplotype (`o3_mat/map/R_unm.mat.bam`: 132 map, 125 of them in one 71-kb stretch of mat chr12 95.96-96.03 Mb) | those 959 reads | the maternal alignment | held out like H: run ONCE, no tuning |

## 3. Rules frozen before the first run on any bed

Decided on bed A only, written here, committed, then H and M are run once with no change. Constants in section 1 are the only ones; the sweep on A reports delta/2, delta and 2 delta for the edge rule and e-values 1e-3 / 1e-5 / 1e-10, and the registered setting is delta and 1e-5.

## 4. Metrics (per bed, per method; read level unless stated; never pooled across beds)

- N_D = unmapped D reads; coverage = share of N_D in a cluster of >= 3; purity = share of clustered D reads whose cluster's majority truth family equals their own; fragmentation = clusters per deleted copy.
- Attribution: attributed / abstained clusters; cluster accuracy = attributed clusters whose majority family is the attributed family.
- **Rescued-correct** = D reads of an attributed cluster attributed to their own family; **wrong joins** = reads of an attributed cluster attributed to another family, and background reads in any attributed cluster.
- Deleted copies with >= 1 rescued-correct read (copies reached) against the copies that have unmapped reads.
- Candidate fidelity: the consensus of each attributed cluster aligned to the erased copy's transcript sequence (asm20): identity x coverage >= 0.999 counts as the copy reproduced from the unmapped reads alone.
- Bed M: whether the chr12 reads form one cluster, the consensus' identity to the maternal region, and what the 834 other reads do.

## 5. Bars (committed now)

- **U1 (the pool is coherent)** on H: coverage >= 0.90 and purity >= 0.99.
- **U2 (the net is better)** on H: rescued-correct >= 2,656 (half of the 5,312 unmapped D reads; A13 rescued 557) with wrong joins <= 5% of rescued reads (A13's registered false-move bound is 5%; it measured 0.77%), AND rescued-correct of the pooled method > B1 > B0 on the same pool.
- **U3 (real)** on M: >= 90% of the 125 chr12 reads in one cluster, consensus identity to the maternal region >= 0.99, and no background read in that cluster.
- Reported, not judged: copies reached, candidate fidelity, behaviour of the sweep, runtime.

## 6. Predictions, before looking

U1 holds (components at delta over reads of one transcript are nearly pure; the risk is fragmentation, not mixing). U2 is uncertain: attribution decides it, and translated search should reach coding families but not pseudogenes, lncRNA or very divergent copies;
I expect the pooled method to beat B1 and B1 to beat B0, and rescued-correct between 20% and 60% of N_D, so U2's 50% bar may FAIL. U3 holds for the cluster (one expressed sequence, 125 reads, 3.6 kb) with no family to attribute it to.

## 7. What would change the claim

If clustering holds (U1) but attribution fails (U2), the supported statement is "the unmapped reads of an erased copy are recoverable as coherent candidate transcripts (flagged, not family-attributed)", which is the register's own allowed class; the family label then needs the annotation-free evidence of the other layers.
If U1 fails on H, the pooled net is refuted as designed.

## 8. Not covered

Reads that map poorly (absorbed on a paralog) rather than unmap; non-excision absences (a copy real in another individual); the full 23 GB library as one pool (cost); any change to the shipped binary; IG/TR and run-exclusivity screens beyond reporting them for flagged clusters.

## 9. Execution

Python on top of minimap2 (`-x ava-pb`), pyabpoa, mmseqs2 and networkx, in `bench/unmapped_rescue/` with unit tests for the graph, the consensus selection, the tie rule and the scoring; heavy jobs foreground and serial under `tools/rlock.sh heavy` (`bash tools/rlock.sh`); work dir `/mnt/linuxdisk/tmp/o3_rescue/`.

---

## Amendment 1 (2026-10-08, after the dev-bed A results below, before bed H or any bed-M-specific step is run)

**Seen on bed A before this amendment (dev):** the seed-round clustering at the registered delta puts 57,664 of the 59,791 deleted-copy reads (96.4%) in clusters of >= 3 at purity 1.0000 (221 clusters, 1 to 7 per family, median 2).
Mmseqs translated attribution of the 221 consensus sequences (section 1 step 3 as written) rescued 77.8% of the reads but 15% of the joined reads were wrong joins and cluster accuracy was 0.59; the wrong clusters were near-ties between paralogous families (median best / runner-up bit-score ratio 1.03).
Replacing the score by the NUCLEOTIDE evidence of the whole consensus changed this: with discontiguous megablast (`blastn -task dc-megablast`, e-value 1e-5) the family score = the number of consensus bases covered by the union of the HSPs of the family's best target, and the dev curve is flat:
strictly-above rule 87.8% rescued / 7.1% wrong joins; relative margin >= 0.10 87.8% / 0.1%; 0.20 87.8% / 0.1%; 0.50 87.7% / 0.0%; an absolute margin of 10 bp 87.8% / 0.6%; a floor on the covered fraction of the consensus costs reads without helping (0.20: 68.3% / 0.7%).

**Frozen for beds H and M (replaces section 1 step 3; everything else stands):** attribute a cluster iff its best family's covered-base score is >= 1.10 x the runner-up family's score (a family with no hit scores 0, so a single-family hit is attributed); otherwise abstain. The 1.10 comes from bed A where 0.10 to 0.50 give the same answer; H and M report the sweep 1.00 / 1.10 / 1.50 and the translated mmseqs score as a comparator.
Comparators B0 (A13's single-read rule) and B1 (single-read dc-megablast, same margin rule) are run on a seeded sample of 1,000 unmapped deleted-copy reads per bed (cost); the pooled method's numbers are on all reads.

**Bed M is not fully held out.** The 959 background reads are part of bed A's pool, so the clustering outcome for the 125 chr12 reads was visible in the dev run. U3's clustering part is therefore reported as dev-visible; its consensus-fidelity part (identity of the consensus to the maternal region) and the attribution outcome are first evaluated after this amendment.

**Additional readouts (reported, not judged):** deleted copies reached; per-family cluster fragmentation; the clusters that abstain (their reads stay unrescued, with size); whether the consensus reproduces the erased copy (identity x coverage >= 0.999 against the unmasked genome at the erased interval).

---

## Amendment 2 (2026-10-08, after the frozen rule had been run once on bed H; written BEFORE any rescue number of the new score exists)

**What was seen on H and why this is exploratory.** Asked whether a family profile would resolve the attribution failures, I (a) scored families by the bases covered by ANY surviving member (a profile proxy): identical to the frozen rule on H (4,255 rescued, 3 wrong, 20/32 copies) and on A (52,493 vs 52,496); and (b) searched each of the 13 abstained H clusters that has a truth family against ONLY its own family's survivors with permissive settings
(blastn word 7, e-value 10; tblastx): 11 of 13 show homology far above the random-family background (for example 475 bits, E = 5e-135, for the 200-read GWFAM158 cluster; random families <= ~150 bits), two do not (GWFAM244: 31 bits; one GWFAM348 cluster: 35 bits).
So the signal exists and the frozen score ranks it badly. This diagnosis used the truth of bed H, so **H is no longer held out for anything derived from it**; the rule below is evaluated first on dev bed A (whose abstained clusters were not inspected), then on H as an exploratory arm and reported as such.

**New score (parameter-free; the margin rule is unchanged).** A family covers a consensus base if any of its copies has a dc-megablast HSP over it; a base covered by n families is credited 1/n to each. The cluster goes to the best family iff its score >= 1.10 x the runner-up (the frozen margin), else it abstains.
Hypothesis: repeat-derived hits give unrelated families near-equal cover, which the frozen score cannot separate from the true family's exon homology; crediting shared bases less should raise the true family's rank.
**Decision (before looking):** the arm is judged by wrong joins <= 5% of rescued reads AND rescued-correct above the frozen rule on A (52,496) and on H (4,255). If it only helps on H it is reported as unsupported (tuned on the held-out bed). The profile-HMM arm (nhmmer on a per-family MSA) is NOT run unless this arm shows the signal is rank-limited rather than absent.

**Amendment 2, result on dev bed A (the arm as registered):** 52,516 rescued (87.8%, +20 reads) but 4,041 wrong joins (7.15% of joined, cluster accuracy 0.667): it FAILS the <= 5% condition. Cause (55 clusters changed, 4 newly right, 51 newly wrong): crediting shared bases 1/n shrinks every score, so noise-level top hits such as 13 against 7 credited bases pass the 1.10 margin. The margin is relative; the arm has no evidence floor.

---

## Amendment 3 (2026-10-08, before any floor number exists)

Add an evidence floor to the Amendment 2 score: the cluster is attributed iff its best family's specific cover score is >= f x the consensus length AND >= 1.10 x the runner-up. **f is selected on dev bed A only:** the smallest f in {0.02, 0.05, 0.10, 0.20} with wrong joins <= 1% of joined reads; if none qualifies the arm is dropped. The selected f is then applied once to H (exploratory: H is no longer held out, see Amendment 2) with the same conditions as the frozen rule (wrong joins <= 5% of rescued; rescued-correct above 4,255).
If the arm beats the frozen rule on both A (52,496) and H (4,255) the claim is "a repeat-robust cover score with an evidence floor attributes more clusters"; if it helps only on H it is unsupported. The per-family HMM profile arm stays unrun unless this arm leaves a rank-limited remainder on A.

**Amendment 2 and 3 results, corrected after the cap fix (see Amendment 4).** With BLAST targets uncapped the arms read: Amendment 2 (parameter-free specific cover, margin 1.10): dev A 55,934 rescued (frozen 55,901), 610 wrong (1.08% of joined); held-out H 4,325 rescued (frozen 4,255), 457 wrong (9.56%).
It passes the <= 5% condition on A and FAILS it on H. Amendment 3 (selection rule: smallest f with wrong <= 1% on A): f = 0.02 (A 270 wrong, 0.48%); on H it gives 4,325 rescued and 271 wrong (5.90% > 5%): FAILS on H. The earlier figures in the Amendment 2 and 3 text (7.15% wrong on A, the 51 newly wrong clusters) were produced with the capped BLAST output.
Post hoc, not selected by the rule: f = 0.05 gives A 55,852 / 14 wrong and H 4,321 / 27 wrong (0.62%), better than frozen on H and worse on A; it is reported, not claimed.

---

## Amendment 4 (2026-10-08, disclosure after the independent review; no rule is changed)

1. **BLAST cap.** The dc-megablast runs used `-max_target_seqs 500`, below the catalog size (753 / 862 / 915); BLAST truncates in a non-score order, so a family's own copy can be dropped (87 of 188 dev consensus sequences and 10 of 39 held-out ones reached the cap). Re-run with 5,000 and a guard (`attribute.capped_queries`):
   bed H is unchanged; dev bed A changes (frozen rule 52,496 -> 55,901 rescued, 87.8% -> 93.5%; strict rule 7.1% -> 1.1% wrong joins, so the margin mattered less on A than Amendment 1 stated). The capped outputs are kept as `*.cap500.*`.
2. **Deviations never amended:** the clustering is seed rounds (`map-hifi` against 2,000 seeds, <= 6 rounds) instead of the all-vs-all of section 1, because the all-vs-all did not finish; the comparator sample is 500 reads per bed, not 1,000; the code implementing the frozen rule was committed after bed H was run (the dev numbers of Amendment 1 reproduce with it); the delta sweep and bed A's rescue file were produced after bed H.
3. **Bed M target set.** The M verdict in the first write-up used bed A's target set (162 families erased). In M's world nothing is erased, so the target set is the full 915-copy catalog.
4. **B0** was specified in the body as the single-read comparator against the surviving genomic copies; it is not A13's rule (A13 also aligned to the family's net reads), so B1 > B0 is vacuous.
5. **Fidelity metric.** The registered metric is identity x coverage >= 0.999; the first write-up reported identity alone.

---

## Amendment 5 (2026-10-08, before any profile number exists): the per-family HMM profile arm (exploratory; bed H is not held out)

**Why now.** After the cap fix, the abstentions that a better model could win are the rank-limited ones: on H the true family ranks below the top two for 6 clusters (228 reads), second for 2 (94) and has no hit while another family does for 3 (89): 411 reads at most (7.7% of N_D); on the two-copy dev bed a profile is a single sequence, so only H can test it.
**Arm.** Universe = the 53 H families (the only catalog families whose surviving copies have mapped reads in the masked run). For each survivor, the transcript consensus = abPOA of its <= 100 longest primary-aligned reads; per family an MSA (MAFFT) of its survivors' consensus sequences and `hmmbuild --dna` (HMMER 3.4 from the `hmmer` mamba env); `nhmmer` of each family HMM against the cluster consensus sequences (e-value <= 1e-3).
Score = the consensus bases covered by the family's nhmmer hits (union of the alignment spans); attribute iff best >= 1.10 x runner-up (the frozen margin). Variant reported, not judged: the best bit score with the same margin.
**Fair comparison.** The frozen rule re-evaluated in the SAME universe (dc-megablast HSPs restricted to the survivors of the 53 H families). The HMM arm is a gain only if it beats that restricted baseline in rescued-correct with wrong joins <= 5% of rescued reads; if a family has a single survivor its profile is that sequence and no gain is expected, so the result is also reported for families with >= 3 survivors separately.
**Not a held-out claim.** The abstained clusters were inspected on H, so a gain here is hypothesis-generating; it would need a new multi-copy erasure bed (erase one copy of other multi-copy catalog families) to be confirmed.

**Amendment 5 result (H; exploratory).** HMM profile: 4,522 rescued, 125 wrong (2.69%), 22 of 24 copies; restricted dc-megablast baseline in the same universe: 4,279 / 59 (1.36%) / 21 of 24; bit-score variant 4,479 / 183 wrong. The condition (more rescued with wrong joins <= 5%) is met: +243 reads, from two clusters (families with 9 and 5 survivors).
Costs: the GWFAM163 profile attracts five wrong small clusters; one cluster lost. Families with >= 3 survivors: 496 rescued / 14.3% wrong vs 226 / 20.7%. Not a held-out result.

---

## Amendment 6 (2026-10-08, after a stop-hook review of the claim; written before any hybrid number exists)

**The claim, stated against the existing approach.** The existing net on bed H is A13's (register rows 1221 and 1223): 557 of the 5,312 unmapped deleted-copy reads joined, 99.6% right. The pool-first net with the frozen rule rescues 4,255 (80.1%) of the same reads, 1 wrong cluster of 27, with the same truth. B1 (single-read dc-megablast with the same score and margin) is a variant of the NEW design, not the existing approach;
it shows the gain over A13 comes from the sensitive nucleotide cover score with a margin, and that pooling adds cost (about 6x less wall-clock) and consensus transcripts but no extra rescued reads (406 vs 410 of 500).
**Hybrid (parameter-free).** Pool first (clusters, consensus, cluster-level attribution as frozen); then every read NOT rescued by an attributed cluster (reads of abstaining clusters, reads in no cluster, background reads) is attributed alone with the same cover score and margin 1.10 (single-read dc-megablast against the same targets, cap 5000).
**Reading, read level, per bed (A dev first, then H):** rescued-correct = deleted-copy reads attributed to their own family; wrong = deleted-copy reads attributed to another family plus any background read attributed to a family; cost = the wall-clock of the extra per-read alignments.
**Decision before looking:** the hybrid supports "pooling does not have to cost rescue" only if on BOTH beds its rescued-correct >= the pooled frozen rule's and >= B1's estimate on the 500-read samples scaled to the bed, with wrong joins <= 5% of joined reads, and its extra alignment work is <= 25% of per-read alignment of the whole pool.
If it does not meet this, the pool-first net is reported as better than A13 and cheaper than per-read attribution but not as a net that rescues more than per-read attribution.

---

## Amendment 7 (2026-10-08, at the user's request, before any number): the reconstructed loci as extra reference (augmented-reference realignment) on held-out bed H

**Idea (user).** Use the loci reconstructed from the reads (the cluster consensus sequences) as additional reference and realign the reads with minimap2, so each read lands on the reconstructed locus (and the consensus can be placed among the family's nodes). No family attribution is needed for this step.
**Design.** Reference = the masked genome (the existing baseline alignments `linktest/R.bam`, `splice:hq -uf --eqx -Y -N 50 -p 0.1`) PLUS the consensus sequences. Instead of a new 13 GB index, the reads are aligned with the SAME minimap2 command to the consensus sequences alone; per read the primary becomes the consensus record iff its alignment score `AS` is strictly higher than the genome primary's (ties stay on the genome; an unmapped read takes any consensus record at query coverage >= 0.80). This is the primary a concatenated reference would give; MAPQ is not reproduced.
**Non-circular test (primary).** Half of the unmapped deleted-copy reads (seeded split, seed 1) build the clusters and consensus sequences (background reads stay in the building half); the OTHER half, which never entered clustering, and the reads that were absorbed on a paralog in the masked run (deleted-copy reads with a genome primary), and the survivors' reads (S) are realigned.
A full-pool version (consensus built from all unmapped reads, mapping the same reads back) is reported as a self-consistency check, not as a result.
**Readings.** (a) held-out unmapped D reads: share that get a consensus primary at identity >= 0.98 whose cluster's majority family is the read's own family; (b) absorbed D reads: share that move to an own-family consensus, and median `de` before and after; (c) FALSE MOVES: S reads whose primary becomes any consensus, and background reads mapping to a family-cluster consensus.
**Bars (before looking):** (c) false moves of S reads <= 5% (A13's bound; A13 measured 0.77%) and background reads <= 5%; (a) >= 80% of held-out unmapped D reads that belong to a family with a consensus in the building half (the clustering coverage on H was 95.7%, so a miss is mostly a copy with no cluster in the half); (b) reported, not judged.
**Prediction:** (a) holds near 90%; (b) moves a minority of absorbed reads (those at the distant paralogs; reads on a 1% paralog stay unless the consensus is closer); (c) holds below 1%.

## Amendment 8 (2026-10-08, at the user's request, before any LRPAP1 number): LRPAP1 and its gorilla copies as the worked example

**Substrate.** KB3781 fibroblast reads of the 11 LRPAP1 loci (2,875 net reads, `o3_mat/map/R_LRP.{mat,pat}.bam`); reference = the maternal haplotype; truth = the paternal one. The LRPAP1 copies are 96 to 99% identical to each other, the regime where a missing copy hides as an absorbed pile (small `de`), not as unmapped reads.
**Pipeline (the same truth-free steps).** (1) the net = LRPAP1 reads whose primary on `mat` has `de` > 0.00958 (the registered allele cutoff) or is unmapped; (2) seed-round clustering of the net (same edge rule), consensus per cluster; (3) align each consensus to `mat` and `pat` (identity x coverage; which locus it lands on, in particular pat chr12 22.55 Mb, the copy IsoCon's transcripts matched); (4) augmented reference as in Amendment 7 (mat primaries from `R_LRP.mat.bam` plus the consensus sequences): which reads move, their `de` before and after, false moves of reads of other copies;
(5) the consensus placed among the family's nodes: its identity to each LRPAP1 copy transcript on `mat` (the copies' RefSeq/Liftoff intervals), nearest node first.
**Readings and bar.** Report each step. Bar for the example to count as a demonstration: a cluster whose consensus matches a pat-only locus at identity x coverage >= 0.999 and a mat locus at < 0.999 (the O3 flag), AND false moves <= 5% of the reads that were on the other copies.

### Amendment 6 result (bed H, 2026-10-08; bed A not run, see below)

| H, 5,313 deleted-copy reads | rescued correct | wrong | unattributed | wrong / joined |
|---|---|---|---|---|
| pooled frozen rule | 4,255 (80.1%) | 3 | 1,055 | 0.07% |
| hybrid (residual reads one by one) | 4,421 (83.2%) | 71 | 864 | 1.58% |
| B1 estimate (410/500 scaled) | 4,357 (82.0%) | | | |

Residual reads: 2,014 of 6,272 (32.1%); extra BLAST 506 s (about 0.25 s per read, so per-read attribution of the whole pool is about 1,580 s).
**Decision rule:** rescued-correct >= pooled (4,421 >= 4,255, met), >= B1 scaled (4,421 >= 4,357, met), wrong <= 5% of joined (1.58%, met), extra alignment work <= 25% of the whole pool (32.1%, **NOT met**).
**Verdict: the hybrid does not meet the registered rule on H.** The pool-first net is reported as better than A13 and cheaper than per-read attribution, not as a net that rescues more than per-read attribution. The hybrid's +166 reads over the pooled rule and +64 over B1 are a descriptive observation, not a claim.
Bed A was not run: the rule needs BOTH beds, so the H failure fixes the verdict and a dev run cannot change it (disclosed; an hour of BLAST saved).

### Amendment 7 result (bed H, 2026-10-08, `bench/unmapped_rescue/run_augment.py`)

Building half: 2,656 of the 5,313 unmapped deleted-copy reads (seed 1) + the 959 background reads, clustered with the frozen rule (coverage 95.1%, purity 1.0000, 44 clusters, 0 background reads in a family cluster), abPOA consensus per cluster. The other 2,657 unmapped D reads never entered clustering.
Reads are aligned to the consensus sequences with the baseline command; a read moves iff the consensus alignment score is strictly higher than its genome primary's (unmapped reads: identity >= 0.98, query cover >= 0.80).

| class (non-circular half run) | n | moved | to own family | median de before -> after |
|---|---|---|---|---|
| held-out unmapped D | 2,657 | 2,505 (94.3%; 94.6% of the 2,648 whose family has a consensus) | 2,505 | unmapped -> 0.0011 |
| absorbed D (genome primary on a paralog) | 11,973 | 2,250 (18.8%; 82.4% of the 2,728 whose family has a consensus) | 2,247 | 0.1237 -> 0.0031 |
| S (survivor copies' reads) | 41,727 | 0 (0.00%) | | |
| background | 959 | 142 (all onto background-only clusters built from them), 0 onto a family cluster | | |

Bars: (a) >= 80% of held-out unmapped D with a family consensus: 94.6%, **met**. (c) false moves of S <= 5%: 0.00%, **met**; background onto a family cluster <= 5%: 0.0%, **met**. (b) reported: absorbed reads move only when the consensus is closer (median de 12.4% -> 0.3%), 99.9% of the moves onto the own family.
Full-pool self-consistency (47 consensus from all unmapped D, same reads mapped back): unmapped D 4,973 of 5,313 (93.6%), absorbed D 2,248 (2,243 own family), S 2 of 41,727 (0.005%, both own-family consensus), background 0 onto a family cluster. Not a result; shown to confirm the half run loses little.
Prediction check: (a) predicted near 90%, observed 94.6%; (b) predicted a minority, observed 18.8% (82% in families with a consensus); (c) predicted below 1%, observed 0.00%.
Caveats: the score compares a spliced genome alignment with an unspliced transcript alignment; S moving 0 of 41,727 is the control that the comparison does not favor the consensus by construction. A concatenated reference was not built (no 13 GB index); the primary is the one it would give, MAPQ not reproduced.

### Amendment 8 result (LRPAP1, KB3781 mother as reference, 2026-10-08, `bench/unmapped_rescue/run_lrpap1.py`)

**Net.** 381 of the 2,875 LRPAP1 reads have maternal divergence `de` > 0.00958 (none unmapped; 138 on the paternal haplotype). Clustering with the frozen edge rule (all net reads were seeds, so this is one all-vs-all round): 4 clusters, 311 reads, purity 1.000 by the paternal primary's locus. Cluster sizes 283 (c01), 21 (c00), 4 (c01), 3 (c00). A seeded half (190 reads) gives 3 clusters (146 c01, 10 c00, 3 c01), same placements.
**Consensus placement** (identity x coverage; all 11 LRPAP1 loci on `pat`, the 9 mat sites on `mat`; `c01|c03` is one mat site because both lift to the same 20 kb):

| cluster (reads) | best pat locus | best mat site (nearest node) | next mat sites | O3 flag at 0.999 |
|---|---|---|---|---|
| 283, 1,492 bp | **c01 0.9987** | c00 0.9886 | p12 0.9880, c05 0.9873, c01\|c03 0.9839 | no (misses by 0.0003) |
| 4, 2,680 bp | **c01 0.9993** | c05 0.9888 | c04 0.9885, c02 0.9884, c01\|c03 0.9855 | **yes** |
| 21, 1,493 bp | c00 0.9980 | c00 0.9946 | | no (an allele of c00, 0.5% from the mother's) |
| 3, 1,496 bp | c00 0.9725 | c00 0.9732 | | no (matches neither) |

**Augmented reference** (reads aligned to the 4 consensus sequences; a read moves iff the consensus score is strictly higher than its maternal primary's). Half run, non-circular: 162 of the 191 held-out net reads move (84.8%), all 162 to the consensus of their own copy; median `de` 0.0115 -> 0.0027. The 2,494 reads already within the allele cutoff on the mother: 1,709 move (1,691 to the consensus of their own copy, 18 to another), median `de` 0.0040 -> 0.0007; 1,189 of 1,432 c00 reads and 505 of 530 c01 reads. Reads whose copy has no consensus (c02 to c05, p12/p14, p16, unplaced; 532 reads): 15 move (2.8%; c02 5, c05 10); p12/p14 (246), p16 (39), c03 (50), c04 (94) 0 moves.
**Bar.** A cluster flagged at the registered 0.999 with the mother below it: yes, the 4-read c01 cluster (0.9993 vs 0.9888). False moves <= 5% of the reads of the other copies: 2.8%, met. Verdict: **met, through the small cluster**; the 283-read cluster (96% of the c01 reads) has the same placement but 0.9987 and is not flagged at the registered bar.
**Not established / caveats.** (1) Truth is the paternal primary's locus (alignment-derived, not independent); p12 and p14 are tied and scored as one group. (2) That 1,691 placed reads move to a closer consensus is mostly the individual's own c00 and c01 sequences replacing the mother's: the move is toward the animal's own transcript, not a rescue of reads that were lost, so it is not counted as a gain. (3) CORRECTED after Amendment 9: the 283-read consensus aligns to pat c01 with identity 1.0000 (1,490 matches, 0 edits); the 0.9987 is coverage, 2 terminal consensus bases that do not align, not a consensus error. (4) The registered lift places c01 on the mat site it shares with c03 (222 mismatches in 20 kb); the consensus is 1.6% from that site (0.9839) and nearer to c00 (0.9886): the reads support c01 as the chr12 locus the mother's assembly lacks within the cutoff, which is evidence, not a test (single animal, one library).

---

## Amendment 9 (2026-10-08, before any polish number exists): an IsoCon-style correction step on the cluster consensus

**Why.** The consensus of a cluster reaches identity x coverage >= 0.999 against the erased copy in only 14 of 40 clusters on H and 65 of 212 on dev A (register drafts). In the LRPAP1 example the 283-read c01 consensus is at 0.9987, 2 bp short of the registered flag. IsoCon's second half is a correction step: align the cluster's reads to the candidate, test every non-candidate allele against the error rate, and change the candidate where the evidence is significant. The pool-first net so far stops at the abPOA consensus. This amendment tests whether that step closes the gap. It was thought of after seeing the LRPAP1 cluster (so LRPAP1 is not held out for it); the polish has not been run on A or H.
**Step (truth-free; reads and consensus only).** For each cluster: (1) align its reads to its consensus (`minimap2 -ax map-hifi --eqx`, primary records), count per consensus position the four bases and deletions, and per gap the inserted strings; (2) `e` = non-consensus observations / all observations of the cluster (a conservative estimate: it includes true variants); (3) an alternative allele with `k` of `n` reads is SIGNIFICANT iff `P(Binomial(n, e) >= k) < 0.05 / (3 L)` (Bonferroni over the three alternatives at each of the `L` positions; an insertion string is tested at its gap with the same correction); (4) a significant allele that is held by more reads than the consensus allele at that column REPLACES it (substitution, deletion of the base, insertion of the string); a significant allele held by fewer reads is only RECORDED (a minority variant, the IsoCon candidate for a split, not applied here); (5) repeat from (1) until a round changes nothing (at most 5 rounds). No other parameter.
**Scoring.** Unchanged (`fidelity.py`): the consensus (before and after) aligned with `minimap2 -c --cs -x splice:hq -uf -N 5` to the unmasked genome (`--cs` only adds a tag), best hit overlapping the erased interval of the cluster's majority family, identity x coverage >= 0.999. Dev A first, then held-out H; the code is committed after A and before H.
**Readings.** R1: clusters at >= 0.999 before and after (A, then H), clusters whose identity x coverage drops, total edit distance (NM) before and after. R2 (dev A only; the genome is used for this diagnostic, never by the polish): among the clusters still below 0.999 after the polish, the share of consensus-vs-genome edit columns that are significant minority-variant columns of the cluster's pileup (allele-like columns).
**Decision, before looking.** The correction step is supported iff on H the count at >= 0.999 rises from 14 to at least 19 of 40 (about a fifth of the 26 below the bar) AND at most 2 of the 40 clusters lose identity x coverage, AND the dev A count rises. If H rises by fewer than 5, the answer for the consensus-error explanation is NO; then, if R2 shows that at least half of the residual edit columns are significant minority-variant columns, the remaining shortfall is allelic (the reads carry two haplotypes, the assembly one) and allele-aware splitting is the next amendment (not run here); otherwise the shortfall is between the reads and the assembly and is not a consensus problem. The +5 bar is chosen here, with no polish data on A or H.
**Also reported (not judged):** the LRPAP1 consensus sequences after the polish, placed as in Amendment 8; run time.
**Not in scope:** splitting a cluster by its variants (needs a within-cluster truth); attribution (unchanged).

### Amendment 9 result (2026-10-08, `bench/unmapped_rescue/{polish,run_polish}.py`; code committed after A and before H)

| | dev A (212 clusters) | held-out H (40 clusters) |
|---|---|---|
| identity x coverage >= 0.999, before -> after | 65 -> 70 (reproduces the stored 65) | 14 -> 15 (reproduces the stored 14) |
| clusters crossing the bar up / down | 8 / 3 | 1 / 0 |
| clusters whose identity x coverage drops / rises | 16 / 37 | **3** / 7 |
| total edit distance (NM) | 17,222 -> 15,942 (-7.4%) | 2,262 -> 2,094 (-7.4%) |
| corrections applied | 4,474 in 53 of 221 clusters (394 sub, 4,040 del, 40 ins) | 176 in 10 of 47 (16 sub, 159 del, 1 ins) |
| R2: residual edit columns that are significant minority-variant columns | 186 of 15,356 (1.2%) | 39 of 2,024 (1.9%) |

**Decision rule: NOT met.** H needed >= 19 of 40 (got 15) and <= 2 clusters losing identity x coverage (got 3); only the dev-A clause (count rises) holds. The consensus-error explanation of the shortfall is **refuted**. R2 is far below half, so the residual is not allelic either; per the rule the shortfall lies between the reads and the assembly, not in the consensus.
**What the corrections were.** 90% are deletions, concentrated in a few clusters whose consensus carried an exon most reads skip (e.g. 4,912 -> 2,882 bp, 2,030 deleted positions): a mixed-isoform cluster trimmed to its plurality isoform, not base-level repair. Base-level substitutions were 394 of 4,474 on A and 16 of 176 on H.
**Exploratory, on dev A only (not a result).** The metric is dominated by something polishing cannot touch: identity alone reaches 0.999 in 95 of 212 clusters but identity x coverage in 65, so 30 clusters are lost to coverage (unaligned consensus ends) alone; the median identity x coverage is 0.9976; 22 clusters have more than 200 edits against the assembled copy (the cluster is not a single copy end to end); NM is 0 in more than a tenth of the clusters. The registered 0.999 bar measures end-to-end agreement with the annotated copy more than consensus accuracy.
**LRPAP1 (reported, not judged).** The polish applied 0 corrections to the 4 consensus sequences; placements are identical to Amendment 8. Checked against the paternal genome: the 283-read consensus aligns to pat c01 over 1,490 of its 1,492 bases with 0 edits (identity 1.0000), so the 0.9987 is 2 unaligned terminal bases. Its 21 "minority variants" are single-base deletions held by 8 to 11 of about 275 reads each (HiFi dropout), not alleles.
**Stage 2 (allele-aware splitting) is not warranted by this amendment's rule** (R2 1.2% and 1.9%).

---

## Amendment 10 (2026-10-08, before any synthetic number exists): a supervised synthetic world, erase one copy, recover it

**Why.** Two things the real beds cannot answer: (i) what the unaligned consensus ends (1 to 4 bp in 21 of the 30 coverage-limited dev-A clusters, 5' extensions of 18 to 1,040 bp in 4; no polyA) are, because the true transcript is unknown there; (ii) how recovery depends on how close the erased copy is to a surviving copy, because the real families mix every divergence. A synthetic world gives exact truth: every read's copy, every copy's transcript and genomic sequence.
**World** (`bench/unmapped_rescue/synth_world.py`, built with `bench/famsim`; seed 20261008; nothing tuned afterwards). 30 families, 6 for each divergence class D in {0.5%, 1%, 2%, 4%, 8%}. Per family a random synthetic gene (5 to 10 exons, internal exons 90 to 350 bp, first 100 to 250, last 250 to 800, introns 300 to 2,500 bp, GT...AG), three copies on independent random background: A (the template, present), B (present, SNP rate 1.5 D), E (ERASED: the contig is omitted from the reference, SNP rate D; its nearest survivor is A at about D). 80 reads per copy (HiFi substitutions 0.001, short indels 0.0003, end jitter 30 bp at both ends; famsim's read model, which has no polyA, adapter or 5' degradation). Read-end variants on the same world: E0 as above; E1 = E0 + a 3' polyA tail of 20 to 30 A on every read; E2 = E0 + 5' truncation (30% of reads lose up to 30% of their length).
**Pipeline (the frozen one, truth-free):** reads aligned to the reference with the study's command; the net = unmapped or divergence > 0.00958; seed-round clustering with the frozen edge rule, components >= 3; abPOA consensus; dc-megablast cover-score attribution against the surviving copies (margin 1.10); realignment of all reads to the consensus sequences (Amendment 7 rule).
**Readings and predictions (bars set now):**
- S1 (the simulation behaves): for erased copies with D >= 2%, >= 95% of their reads are in the net; for D = 0.5%, < 20% (designed-undetectable: the cutoff is about 1%).
- S2 (consensus): for D >= 2%, a cluster of >= 3 reads exists for >= 90% of the erased copies, purity 1.000 by true copy, and its consensus aligns to the TRUE transcript of the erased copy with identity >= 0.999 over the aligned part in >= 90% of those clusters.
- S3 (the registered metric): identity x coverage >= 0.999 of the consensus against the true genomic copy (the metric used on the real beds), same clusters. Reported next to S2. If S3 << S2 the ends are the cause.
- S4 (ends, descriptive): per end, the consensus overhang beyond the true transcript and the true bases it lacks, and the aligner's soft clip, per read-end variant. Prediction: with jitter only, the consensus ends fall inside the read-end envelope (within 30 bp of the true ends); polyA (E1) appears as a 3' overhang of about 25 bp; the 1 to 4 bp clips of the real beds are NOT reproduced by E0.
- S5 (attribution): net reads of D >= 2% erased copies attributed to the correct family >= 90%, wrong <= 1%.
- S6 (augmented reference): erased-copy reads in the net land on their own consensus >= 90%; reads of surviving copies move <= 1%.
- S7 (the divergence curve, descriptive): S1, S2, S5 and S6 by D.
**Decision.** The supervised check "the pipeline recovers an erased copy when it is recoverable" passes iff S1, S2, S5 and S6 hold on E0. The ends question is answered by S3 vs S2 and S4: if E0 reproduces the real-bed signature (S3 << S2, clips of 1 to 4 bp) the mechanism is in the pipeline (aligner or consensus end handling); if only E1 or E2 reproduces it, it is in the reads; if none does, the real-bed ends come from something the simulation lacks (the individual's variation against the assembly, UTR structure), and the next step is read-end inspection on real data.

### Amendment 10 result (synthetic world, 2026-10-08; `synth_world.py`, `run_synth.py`, `ends.py`; seed 20261008; famsim `verify` 28/28 PASS on the families checked)

30 families x 3 copies (A, B present; E erased), 80 reads per copy, 7,200 reads; E0 = jitter only. The net is every read unmapped or with divergence > 0.00958.

| D (E vs nearest survivor) | E reads in the net (E0) | E copies with a pure cluster | consensus vs TRUE transcript, core identity >= 0.999 | identity x coverage >= 0.999 vs true genome copy (E0 / E1 polyA / E2 trunc) | attributed to the right family, wrong | E reads landing on own consensus | survivor reads moved, AS rule / divergence rule |
|---|---|---|---|---|---|---|---|
| 0.5% | 12.1% (E2 14.8%) | 2 of 6 | 2 of 2 | 2/2, 0/2, 2/2 | 57 of 58, 0 | 57 of 58 | 306 of 960 / 0 |
| 1% | 99.2% | 6 of 6 | 6 of 6 | 6/6, 0/6, 6/6 | 476 of 476, 0 | 476 | 480 of 960 / 0 |
| 2% | 100% | 6 of 6 | 6 of 6 | 6/6, 0/6, 6/6 | 480 of 480, 0 | 480 | 289 of 960 / 0 |
| 4% | 100% | 6 of 6 | 6 of 6 | 6/6, 0/6, 6/6 | 480 of 480, 0 | 480 | 0 / 0 |
| 8% | 100% | 6 of 6 | 6 of 6 | 6/6, 0/6, 6/6 | 480 of 480, 0 | 480 | 0 / 0 |

**S1 holds** (D >= 2%: 100% in the net; D = 0.5%: 12.1%, designed-undetectable; at 1% the cutoff is crossed, 99.2%). **S2 holds** (every D >= 1% copy has a pure cluster and its consensus is >= 0.999 identical to the true transcript over the aligned part; the 0.5% class has clusters for only 2 copies because only 58 of its 480 reads are in the net). **S5 holds** (0 wrong). **S6 FAILS under the registered score rule**: reads of the SURVIVING copies move onto the erased copy's consensus in 30 to 50% of cases for D <= 2% (bar <= 1%); they are all A reads (de 0.001 to 0.003 on the genome, AS 1221 to 1248) that score higher on the unspliced consensus (AS 1279, de 0.0138): a spliced genome alignment pays about 8 points per junction that an unspliced transcript alignment does not, so raw scores are not comparable. **Registered decision (S1, S2, S5, S6 all hold on E0): NOT met, because of S6 only;** recovery itself (S1, S2, S5) holds. The failure is in the move rule of Amendment 7; real bed H (D mostly >= 2%) did not show it (0 of 41,727), but a family whose erased copy is within about 2% of a survivor would.
**S3/S4, the ends question:** E0 does not reproduce the real-bed signature: identity x coverage equals the core identity (6/6), the aligner clips 0 bases (median), and the consensus is SHORTER than the true transcript (5' lacks 0 to 23 bases, 3' lacks 0 to 7: the read-end jitter), never longer. E2 (5' truncation) is the same. E1 (polyA) fails the metric completely (0 of 30 families; the 3' tail is 13 to 30 bases of the consensus, clipped by the aligner at a median of 24 to 28) which the real beds do not show (no A-rich clips). **None of E0, E1, E2 reproduces the real beds, so by the registered rule the real-bed ends come from something the simulation lacks; real-data inspection is the next step** (Amendment 11).

---

## Amendment 11 (2026-10-08, AFTER the Amendment 10 numbers; exploratory, post hoc, labelled as such)

Two things the synthetic world exposed, followed up on the real beds. Neither is a registered test; both are descriptive, dev A first. H was already used by the earlier amendments, so H numbers here are not held-out.
**(a) A divergence move rule.** A read moves to a consensus iff the consensus alignment's `de` is strictly lower than its genome primary's (unmapped reads: as before). Synthetic world: survivors moved 0 of 960 in every class and variant (AS rule: up to 480), erased-copy reads in the net still all land on their own consensus, and in the 0.5% class (copy hidden inside the allele cutoff) 103 of 422 absorbed erased-copy reads move to their own consensus where one exists. Real bed H, half run: held-out unmapped reads 2,505 of 2,657 (identical to the score rule), absorbed 2,122 of 11,973 (score rule 2,250), survivors 0 of 41,727: H's result does not depend on the rule.
**(b) The unaligned consensus ends of the real beds are an untemplated 5' run of G.** In consensus (transcript) orientation, on dev A 129 of 212 consensus sequences have a 1 to 20 bp unaligned 5' end, 70% G, 102 of them pure G; lengths 1 bp (56), 2 bp (44), 3 or more (29). None of 127 clips of 4 bp or less equals the adjacent genomic sequence (a random same-size window matches 24 of 127), so the aligner is not clipping genomic bases: the consensus carries bases that are not in the genome. 3' clips are rare (20) and not composition-biased. This is the template-switching / cap artifact (a few untemplated G at the 5' end of the cDNA); the earlier read-ends note (`project_read_proven_ends`) recorded no such signal for gorilla, which this contradicts for these consensus sequences and needs reconciling. Excluding a pure 5' G run of at most 3 bp from the coverage denominator turns identity x coverage >= 0.999 from 65 to 80 of 212 clusters (dev A) and from 14 to 17 of 40 (H): a metric that ignores the artifact would score the consensus higher, not a change to any registered result. The remaining shortfall is internal identity (identity alone >= 0.999 in 95 of 212) and 5' extensions of 18 to 1,040 bp (4 of 30 coverage-limited clusters).
**Next, proposed (not run):** add the untemplated 5' G run to the read simulator so the supervised world reproduces the real signature, then test trimming it from the consensus with exact truth.

---

## Amendment 12 (2026-10-08, before any number from the fresh world): adopt the divergence move rule in the augmented reference

**Decision.** The augmented-reference step (Amendments 7, 8, 10) moves a read to a consensus iff the consensus alignment's divergence `de` is strictly lower than its genome primary's (an unmapped read: consensus at identity >= 0.98 and query cover >= 0.80, as before). It replaces the raw-score rule, which compares a spliced genome alignment with an unspliced transcript alignment (about 8 points per junction in favour of the transcript) and moved 30 to 50% of surviving-copy reads in the synthetic world when the erased copy was within 2% of a survivor. The score rule stays available as an explicit option (`rule="score"`) so earlier results can be reproduced. The rule was found on the synthetic world of seed 20261008 (Amendment 10), so that world is no longer held out for it.
**Test (fresh world, never looked at).** `synth_world.py build --seed 20261009` with eight divergence classes, four families each (32 families): D in {0.5%, 0.75%, 1%, 1.5%, 2%, 3%, 4%, 8%}, three of which (0.75%, 1.5%, 3%) were not in the first world. Everything else as in Amendment 10 (E0 read model, frozen pipeline).
**Bars (set now).** (1) Survivor reads moved onto an erased-copy consensus <= 1% in every class (the score rule gave up to 50%). (2) Erased-copy reads in the net land on their own family's consensus >= 90% in every class with a pure cluster for at least 3 copies. (3) Recovery is unchanged by the rule: copies with a pure cluster and attribution results identical to the score-rule run on the same world (the rule only touches the move step). Reported, not judged: erased-copy reads absorbed inside the allele cutoff (not in the net) that move to their own consensus, by class.
**Also re-run under the new default (reported):** bed H half run (previous: held-out unmapped 2,505 of 2,657, absorbed 2,122, survivors 0 with the divergence rule) and the LRPAP1 example (Amendment 8).
**Decision rule.** The adoption stands iff bars (1) to (3) hold. If (1) fails in any class, the rule goes back to explicit-option status and the failure is reported.

### Amendment 12 result (fresh world seed 20261009, 32 families; 2026-10-08)

Net 2,015 of 2,560 erased-copy reads, 0 survivor reads; 26 clusters, 26 attributed.

| D | E reads in net | E copies with a pure cluster | survivor reads moved, divergence rule (score rule) | E net reads on own consensus |
|---|---|---|---|---|
| 0.5% | 0% | 0 of 4 | 0 of 640 (0) | |
| 0.75% | 54.7% | 3 of 4 | **0** (320) | 175 of 175 |
| 1% | 75.0% | 3 of 4 | **0** (240) | 240 of 240 |
| 1.5% | 100% | 4 of 4 | **0** (241) | 320 of 320 |
| 2% | 100% | 4 of 4 | **0** (71) | 320 of 320 |
| 3%, 4%, 8% | 100% | 4 of 4 | 0 (0) | 320 of 320 |

**Bars: (1) survivors moved 0 of 640 in every class (<= 1%): met. (2) erased-copy net reads on own consensus 100% in every class with >= 3 pure clusters: met. (3) recovery unchanged: the move rule acts after clustering and attribution, which are identical to the score-rule run: met.** The adoption stands; `augment.new_primary` and `move_metrics` now default to the divergence rule (`rule="score"` reproduces Amendments 7, 8, 10).
The detection edge sits where the cutoff says: 0% of the erased copy's reads leave it at 0.5%, 55% at 0.75%, 75% at 1%, all from 1.5%. Reported: erased-copy reads inside the cutoff that move to their own consensus: 65 of 145 at 0.75% (where a consensus exists).
**Re-runs under the new default.** Bed H half run: held-out unmapped 2,505 of 2,657 (94.3%), absorbed 2,122 of 11,973, survivors 0 of 41,727, background onto a family cluster 0 of 959: the H conclusions do not depend on the rule. LRPAP1 (Amendment 8): held-out net reads 172 of 191 move (90.1%; 162 under the score rule), already-placed reads 1,121 move (1,709 under the score rule; 1,118 to the right copy), other-copy reads 3 of 532 (0.6%; 15 under the score rule).

---

## Amendment 13 (2026-10-08, before any number): the untemplated 5' G run, simulated and trimmed

**What is known (Amendment 11, dev A).** 129 of 212 consensus sequences have an unaligned 5' end; 102 are pure G, 1 to 2 bp, none matches the adjacent genome. From the reads of the 53 dev-A clusters (>= 10 reads) whose consensus carries a 1 to 3 bp 5' G run: the median cluster has 88% of its reads starting with that run (quartiles 74%, 94%); over all their reads the leading G run is 0 bp in 5%, 1 bp 21%, 2 bp 41%, 3 bp 23%, 4 bp 7%, 5 or more 3%. Reads are in transcript orientation (the consensus is, and its 3' end is not polyA).
**Simulation (read variants on a fresh world, seed 20261010, same 30-family design as Amendment 10).** E4 = reads with end jitter 3 bp (a dominant start site, unlike E0's 30; no other change). E3 = E4 plus an untemplated 5' G run on 95% of reads, run length drawn from the measured distribution above (1: 0.22, 2: 0.43, 3: 0.24, 4: 0.07, 5: 0.03, 6: 0.01 of the reads that carry it).
**Trimming rule T (truth-free; uses what the pipeline already has).** For a cluster attributed to a family: take the surviving copy of that family with the largest union of dc-megablast HSP query spans (the cover score's best copy); the 5' PREFIX of the consensus is its bases before the first HSP's query start on that copy. If the prefix is 1 to 3 bases and all G, remove it from the consensus; otherwise change nothing. Unattributed clusters are not touched. No other parameter.
**Readings and bars (set now), classes D >= 1% (30-D classes with a pure cluster):**
- (1) E3 reproduces the real signature: identity x coverage >= 0.999 against the TRUE genome copy in <= 50% of E3 clusters (real beds: 31% dev A, 35% H) while >= 90% in E4 (the control). If not, the simulation does not capture the artifact and the rest is not interpreted.
- (2) T repairs it: after T, identity x coverage >= 0.999 in >= 90% of E3 clusters, and the 5' offset of the consensus against the TRUE transcript (ends.py) within [-4, +1] in >= 90%.
- (3) T does no harm: in E4 (no artifact) T trims at most 1 cluster; in E3 T trims a base that belongs to the true transcript (5' offset after T below -4) in at most 1 cluster.
- (4) Reported: T's decisions on E3 by D class (trimmed, not trimmed because abstained, not trimmed because prefix too long or not G).
**Real data (descriptive, dev A first; H is not held out for diagnostics).** Apply T with the real attribution and the real HSPs: how many consensus sequences it trims, agreement with the genome-based finding (a pure-G 1 to 3 bp 5' clip against the erased copy): precision, recall; and identity x coverage >= 0.999 after T by excluding the trimmed prefix from the denominator (no re-alignment needed when the prefix is within the unaligned clip).
**Decision.** T is adopted as a consensus clean-up step iff (1) holds and (2), (3) hold. If (1) fails, no conclusion about T; if (2) or (3) fails, T is reported as not adopted with the failure.

### Amendment 13 result (fresh world seed 20261010; 24 erased copies with a pure cluster at D >= 1%; 2026-10-08)

| | E4 (jitter 3, no artifact) | E3 (E4 + 5' G run) untrimmed | E3 after T |
|---|---|---|---|
| identity x coverage >= 0.999 vs true genome copy | 24 of 24 | **14 of 24 (58%)** | **24 of 24** |
| 5' offset vs true transcript within [-4, +1] | 24 of 24 | 21 of 24 | 24 of 24 |
| consensus with an unaligned 5' clip | 0 of 24 | 20 of 24 (lengths 2: 12, 1: 7, 3: 1) | |
| clusters trimmed by T / over-trimmed (offset below -4) | 0 / 0 | | 21 / 0 |

**Bar (1) NOT met (58% > 50%):** E3 produces the 5' clip itself (20 of 24 consensus sequences, real dev A 129 of 212) but, with an error-free consensus, a 1 to 2 bp clip on a 2.5 kb transcript scores 0.9996 or 0.9992 and still passes 0.999; only runs of 3 or a short transcript fail. The real beds fall to 31 to 35% because the clip adds to other small deviations. Bars (2) and (3) hold. **By the registered decision (needs (1)), T is not adopted.**
**On the real beds T does almost nothing (descriptive):** dev A: trimmed 4 of 212 (3 agree with the genome-based pure-G clip), prefix longer than 3 in 115, abstained 87; H: trimmed 0 of 40 (prefix too long 25, abstained 13). The dc-megablast HSP does not start at the first consensus base on a real survivor (a short or divergent first exon starts the HSP tens of bases in), so the HSP prefix is the wrong measurement outside the synthetic world, where the first exons are identical. identity x coverage >= 0.999 unchanged (65 and 14).

---

## Amendment 14 (2026-10-08, after the Amendment 13 numbers, before any T2 number): rule T2, aligner-based prefix

**Change.** The 5' prefix is measured by aligning the consensus (minimap2 `-c -x splice:hq -uf -N 5`, the command used for the fidelity metric) to the genomic spans of the surviving copies of its attributed family and taking the consensus bases before the alignment start (PAF query start) of the best record by matches; the decision step is unchanged (1 to 3 bases, all G, else no change; abstained clusters untouched). Nothing else changes.
**Test.** (a) Fresh synthetic world, seed 20261011 (same design, E3 and E4): bars (2) and (3) of Amendment 13, and bar (1') replacing (1): E3 shows an unaligned 5' clip in >= 50% of the consensus sequences and E4 in 0 (the property the real beds have: 61% dev A). (b) Real dev A, descriptive but with a bar for the decision: among attributed clusters that have a genome-based pure-G 1 to 3 bp 5' clip (against the erased copy), T2 recovers >= 80% of them, and >= 95% of T2's trims coincide with such a clip. H and the identity x coverage count after T2 are reported, not judged.
**Decision.** T2 is adopted as a consensus clean-up iff (1'), (2), (3) hold on the synthetic world AND the dev-A recall and precision bars hold. Otherwise it is reported as not adopted with the failure.

### Amendment 14 result (2026-10-08; synthetic world seed 20261011, real dev A and H)

**Synthetic (23 erased copies with a pure cluster at D >= 1%).** E4 (no artifact): identity x coverage >= 0.999 in 23 of 23, T2 trims 0, no over-trim (bar (3)). E3: an unaligned 5' clip in 22 of 25 consensus sequences (E4: 0 of 25), bar (1') met. Untrimmed E3: 12 of 23 (52%) reach 0.999; **after T2: 21 of 23 (91%)**, 5' offset within [-4, +1] in 22 of 23 (96%), over-trim 0, 18 clusters trimmed (the 2 failures: one prefix longer than 3, one not pure G). Bars (1'), (2), (3) met.
**Real dev A (the decision bar): not met.** T2 acts on 3 of 212 consensus sequences: `no_hsp` 122 (no alignment to any surviving copy of the attributed family), `abstained` 87, `too_long` 3. Of the 60 attributed clusters that have a genome-based pure-G 5' clip, T2 trims 0 (recall 0%, bar 80%). H: 0 of 40 (13 attributed clusters with a clip, 0 trimmed). The cause is alignment sensitivity: `minimap2 -x splice:hq` against the survivors yields any record for only 44 of 221 consensus sequences (median block 280 bp), because real survivors are 3 to 15% diverged paralogs, not the near-identical copies of the synthetic world. identity x coverage >= 0.999 unchanged (65, 14).
**Decision: T2 is NOT adopted** (the synthetic bars hold, the real-data recall bar fails). What the synthetic world establishes: when the survivors align at the first exon, the pure-G prefix is identifiable and trimming it restores the registered metric with no harm. What it does not: a way to identify the artifact on real families whose survivors are too divergent to anchor the 5' end. The artifact is a property of the library (a 5' soft-clipped pure-G run on reads that align cleanly elsewhere), which is not used by T or T2.

---

## Amendment 15 (2026-10-08, before any T3 number): rule T3, the artifact read from the library

**What the library shows (design check on the 41,727 surviving-copy reads of the real library, `linktest/R.bam`; looked at before this amendment).** Of 39,351 cleanly aligned survivor reads (primary, MAPQ >= 10, divergence <= 0.0096), 28,581 (72.6%) have a 5' soft clip that is a pure G run of 1 to 3 bases (1 bp 12,770; 2 bp 14,843; 3 bp 968); the same clip made of A, T or C: 62, 12, 5. The artifact is a property of the library, readable from reads that align perfectly elsewhere.
**Rule T3.** (1) LIBRARY GATE: from the cleanly aligned reads of the surviving copies count the 5' soft clips (read orientation) that are a pure single base of length 1 to 3; the artifact is declared present iff the G count has P(Binomial(n, 0.5) >= G count) < 1e-6, n = the pure-A/C/G/T clip count (needs n >= 20), i.e. G is not just one of four equally likely bases. (2) If present, remove the maximal leading G run of every consensus sequence if its length is 1 to 3 (a longer run is left alone); no reference, no attribution. If absent, nothing is trimmed.
**Test.** (a) Fresh synthetic world, seed 20261012 (E3 and E4; the gate is computed from the survivor reads of that world's own reference alignment). Bars: gate TRUE on E3 and FALSE on E4; E4: 0 consensus sequences trimmed; E3: identity x coverage >= 0.999 vs the true genome copy in >= 90% of the D >= 1% erased-copy clusters and the 5' offset against the true transcript within [-4, +1] in >= 90%. (b) Real dev A (decision bars, genome-based finding as reference): gate TRUE; of the consensus sequences with a pure-G 1 to 3 bp 5' clip against the erased copy, T3 trims >= 95%; among T3's trims, the trimmed run exceeds the unaligned clip by at most 1 base in >= 90% (templated G's are not lost); identity x coverage >= 0.999 after T3 (excluding the trimmed bases from the denominator, trimmed bases that the aligner had covered subtracted from the aligned length) >= 75 of 212. H reported, not judged (descriptive prediction 16 or more of 40).
**Decision.** T3 is adopted as the consensus clean-up iff the gate behaves (TRUE on E3 and real, FALSE on E4) and every bar in (a) and (b) holds; otherwise not adopted, with the failure.

### Amendment 15 result (2026-10-08; synthetic world seed 20261012; real library; dev A and H)

**Gate.** Real library: OPEN (G 28,614 of 28,694 pure 5' clips; A 63, T 12, C 5; 39,403 clean survivor reads; p = 0). Synthetic E3: OPEN (G 3,829, others 0); E4: closed (A 3, C 5, G 2, T 1; p = 1). The gate behaves.
**Synthetic (D >= 1%, 24 erased copies with a pure cluster).** E4: 0 trimmed, 24 of 24 at 0.999 untouched. E3: identity x coverage >= 0.999 in **12 of 24 -> 22 of 24 (92%)**, 5' offset within [-4, +1] in 14 -> **23 of 24 (96%)**, over-trim 0; the two untrimmed have a leading run longer than 3. Bars (a) met.
**Real dev A (decision bars): two of three FAIL.** T3 trims 197 of 212 consensus sequences (too_long 13, no run 2). Recall among the 100 consensus sequences with a genome-based pure-G 1 to 3 bp clip: **92 of 100 (92%; bar 95%) FAILS** (the other 8 have a leading run of more than 3). Templated-sequence check: the trimmed run exceeds the unaligned genome clip by at most 1 base in **132 of 197 (67%; bar 90%) FAILS** (excess 0: 39, 1: 93, 2: 51, 3: 14). identity x coverage >= 0.999: 65 -> **82 of 212** (bar >= 75) met. H (descriptive): 14 -> 18 of 40; trims 37 of 40, 19 of the 20 with a genome clip.
**Decision: T3 is NOT adopted** (two of the three registered real-data bars fail). What this shows: the artifact is a library property that a gate can read perfectly and that a leading-G trim removes, restoring the registered metric on synthetic worlds (and on the real beds, 65 -> 82 and 14 -> 18). What it cannot do is separate artifact G from templated G: real 5' ends are GC-rich, and 105 of the 197 trimmed consensus sequences have no unaligned G clip at all, so their leading G run was aligned to the erased copy (templated, or artifact G that matched the G-rich upstream flank by chance, which the alignment cannot tell apart). The registered over-trim bar uses the unaligned clip as the artifact's length, which therefore under-counts it; on real data the bar was mis-specified for what can be observed, but it was registered and it failed.

---

## Amendment 16 (2026-10-08, before any T4 number): rule T4, per-cluster correction of the 5' G run from the reads' own leading-G lengths

**Why.** T3 trims every leading G run of 1 to 3 and cannot tell the artifact's G's from the transcript's own (real 5' ends are GC-rich). The reads can: every read of a cluster carries the same templated run of j G's, then an artifact run a_i whose length distribution is a property of the library. So a read's leading G run is r_i = j + a_i (read starts are within a base or two of each other; the model is deliberately simple).
**Rule T4 (gate as in Amendment 15, nothing else added).** (1) Library artifact length distribution a(l), l = 0..8, from the cleanly aligned survivor reads: a(l) = (reads whose 5' clip is a pure G run of l bases) / (clean reads), a(0) = reads with no clip / clean reads; a floor of 1e-4 for lengths never seen. (2) For each cluster with >= 5 reads: r_i = the leading G run of each read's own sequence (no alignment, no CIGAR); j_hat = argmax over j in 0..6 of the sum of log a(r_i - j), a term with r_i < j taking 1e-6; ties go to the larger j. (3) The consensus has a leading G run of k; remove the first max(0, k - j_hat) bases (all G). Gate closed, or a cluster of fewer than 5 reads: no change.
**Test world (fresh, seed 20261013; 30 families, same classes).** Templates whose first exon begins with j = 0, 1, 2 or 3 G's (forced, equally many families) to reproduce the GC-rich 5' ends; E5 = end jitter 1 + the measured artifact (95% of reads, mode 2); E6 = the same without the artifact (control). Truth: s* = the modal read start (chain offset) among the cluster's reads; the correct consensus is the true transcript from s* on, so the 5' offset against the true transcript should be -s*; err = (offset + s*), 0 exact, > 0 artifact left in, < 0 templated bases removed. T3 (trim all) and the untrimmed consensus are scored on the same clusters.
**Bars (set now), D >= 1% erased-copy clusters:** (1) gate open on E5, closed on E6 (E6: nothing trimmed). (2) T4 err = 0 in >= 80% and |err| <= 1 in >= 95%. (3) T4 removes templated bases (err < 0) in <= 10% and leaves artifact (err > 0) in <= 10%. (4) T4 has strictly fewer err < 0 clusters than T3. (5) identity x coverage >= 0.999 against the true genome copy in >= 90% of clusters after T4.
**Real dev A (one bar, rest descriptive):** among consensus sequences with a genome-based pure-G 5' clip of g bases, T4 trims at least g in >= 90% (it does not leave unaligned artifact G's); identity x coverage >= 0.999 after T4 >= 75 of 212; reported: distribution of j_hat and of the trimmed amount against T3's. H descriptive.
**Decision.** T4 is adopted as the consensus clean-up iff bars (1) to (5) and the dev-A bar hold; otherwise not adopted, with the failure.

### Amendment 16 result (2026-10-08; synthetic world seed 20261013 with forced leading G; real dev A and H)

**Scoring check first.** The control E6 (no artifact, same templates) is exact in 24 of 24 clusters under all three treatments, so the err measure (5' offset against the true transcript + the modal true read start) is sound. (A first scoring run used the E0 read-start table by mistake, gave nonsense errors up to 30, and was discarded before any reading; the per-variant truth tables fixed it, the world is byte-identical.)
**Synthetic E5 (artifact + 0 to 3 templated leading G), D >= 1%, 24 erased-copy clusters:**

| | untrimmed | T3 (trim all, 1 to 3) | **T4** |
|---|---|---|---|
| err = 0 (exact) | 0 | 7 | **14 (58%)** |
| \|err\| <= 1 | 4 | 16 | **24 (100%)** |
| templated bases removed (err < 0) | 0 | 7 | **1 (4%)** |
| artifact G left (err > 0) | 24 | 10 | **9 (38%)** |
| identity x coverage >= 0.999 | 16 | 20 | **24 (100%)** |

j_hat equals the forced templated count or one less in every cluster (6, 6, 6, 7 clusters for j = 0 to 3). The 9 clusters with a G left all have the modal read start at transcript position 1 (jitter 1: half the reads start at 0, half at 1), where the templated count at the consensus start is one lower than the global estimate: a +1 ambiguity of the start itself.
**Bars: (1) met** (gate open on E5, closed on E6, nothing trimmed there). **(2) FAILS:** exact 58% (bar 80%); |err| <= 1 in 100% (bar 95%) met. **(3) FAILS:** templated removed 4% (met), artifact left 38% (bar 10%). **(4) met** (1 against T3's 7). **(5) met** (100%).
**Real dev A (bar: T4 removes at least the unaligned genome clip in >= 90%): FAILS, 64 of 100 (36% under-trimmed).** Of the 100 consensus sequences with a genome G clip, 77 come from clusters of >= 5 reads (T4 applicable): 64 ok, 13 short by one base; the other 23 are clusters of fewer than 5 reads, where T4 does nothing by rule. identity x coverage >= 0.999: 65 -> **80** of 212 (bar >= 75, met). j_hat: 0 in 98 clusters, 1 in 51, 2 in 4 (no estimate for 59 small clusters). T3 would remove more than T4 in 104 clusters (1 base: 65, 2: 30, 3: 9). H (descriptive): 14 -> 17 of 40; 6 of 20 clipped consensus sequences short (4 of them clusters < 5 reads).
**Decision: T4 is NOT adopted by the registered rule** (bars (2), (3) and the dev-A bar fail). In its favour and descriptively: it is the only treatment with |err| <= 1 everywhere, it removes templated bases 6 times less often than T3, it restores the metric in all synthetic clusters, and on dev A its shortfall is the 1-base start ambiguity (13 of 77) plus clusters too small for the estimate (23 of 100). The bars of (2) and (3) asked for exactness that the 1-base jitter of the read starts rules out.

---

## Amendment 17 (2026-10-08, before any gate-aware flag number): the registered identity x coverage ignores a leading G run when the library gate is open

**Change (score side only; the consensus is never touched).** The registered metric (identity x coverage >= 0.999 of the consensus against a genome copy, which is also the O3 flag's recovery hit) is ident x (aligned query bases / consensus length). Gate-aware version: if the library gate of Amendment 15 is open AND the consensus' unaligned 5' end (the bases before the alignment start, consensus = transcript orientation) is 1 to 3 bases, all G, those bases are removed from the denominator; in every other case the metric is unchanged. A 5' end that is not pure G, longer than 3, or any gate-closed library is scored exactly as before. 3' ends are never touched.
**Tests.** (a) Synthetic, fresh data not used for this rule: the worlds of seeds 20261012 (E3 artifact, E4 control) and 20261013 (E5 artifact + 0 to 3 templated leading G, E6 control), D >= 1% erased-copy clusters. Bars: gate-aware >= 0.999 in >= 90% of the clusters of E3 and of E5 (registered metric: 12 of 24 and 16 of 24); control worlds (gate closed) identical to the registered metric in every cluster; **no new false passes**: every cluster that passes only under the gate-aware metric has core identity >= 0.999 against the TRUE transcript (0 exceptions). (b) Unit test: a non-G or longer-than-3 5' end is scored as before. (c) Real dev A and H, reported not judged (the 80 and 17 of Amendment 11 were computed with the same exclusion from the genome clip, so these are not blind): gate-aware count. (d) LRPAP1 (Amendment 8): the O3 flag of each of the four consensus sequences, registered and gate-aware, from the existing alignments to `pat` and `mat` (a cluster's flag = best pat locus >= 0.999 and best mat site < 0.999).
**Decision.** Adopted as the flag/fidelity metric iff every bar in (a) and (b) holds.

### Amendment 17 result (2026-10-08)

**(a) Synthetic, D >= 1% erased-copy clusters (24 per variant):** E3 (artifact): registered 12 -> **gate-aware 23 (96%)**; E5 (artifact + 0 to 3 templated leading G): 16 -> **24 (100%)**; controls with the gate closed (E4, E6): identical in 24 of 24 clusters (24 of 24 at 0.999 both ways); **new passes without core identity >= 0.999 against the true transcript: 0** in all four. **(b)** unit tests: a non-G or longer-than-3 5' end, a closed gate, and a consensus with a real mismatch (99.5% identity) are scored as before. Bars (a), (b) met: **the gate-aware metric is adopted** (`flagmetric.gate_aware`; Python-side scoring and reporting only; the Rust `o3_candidates` flag is not changed).
**(c) Real beds (not blind: the same exclusion gave 80 and 17 in Amendment 11):** dev A 65 -> 80 of 212, H 14 -> 17 of 40.
**(d) LRPAP1 O3 flag (best `pat` locus >= 0.999 and best `mat` site < 0.999), registered -> gate-aware:** 283-read c01 consensus: pat c01 0.9987 -> **1.0000**, mat best c00 0.9886 -> 0.9899: **not flagged -> flagged**. 4-read c01 cluster: 0.9993 -> 1.0000 vs mat 0.9895, flagged both ways. 21-read c00 cluster: pat c00 0.9980 -> 0.9993 vs mat c00 0.9946 -> 0.9959: **not flagged -> flagged**; this is the c00 allele of Amendment 8, 0.4% from the mother's, so the flag as defined does not separate an allele from an absent copy (it never did; the registered bar is identity, not the 0.958% allele cutoff). 3-read cluster: its 5' end (GAG) is not pure G, unchanged, not flagged.

---

## Amendment 18 (2026-10-08, before any partition number): IsoCon-style partition of a cluster into competing candidate transcripts

**Why.** The frozen pipeline reports ONE consensus per cluster. A cluster that mixes two isoforms (exon skipping) or two near-identical copies (< 0.958% apart, joined by the edge rule) is reported as the majority isoform / a column-wise blend, and the other member is lost. Amendment 9 tested base-level correction (refuted) and the allele-like columns (1 to 2%), not the separation of isoform- or copy-scale members. IsoCon's remaining idea is exactly that: candidates compete for reads and a candidate stands only if statistically supported.
**Algorithm PART (no reference, no attribution).** For a cluster of >= 6 reads (at most 400, seeded sample for discovery): (1) its consensus C (abPOA of the 100 longest reads, as frozen); (2) align the reads to C (`minimap2 -ax map-hifi --eqx`), take the significant variant columns V of the Amendment 9 test (substitution, deletion, insertion alleles at Bonferroni 0.05; both the alleles the consensus carries as minority and those it lacks), dropping the first and last 15 consensus columns (read ends carry the jitter and the 5' artifact); (3) for every pair of variants in different columns, a one-sided hypergeometric (Fisher) test of positive association over the reads covering both, edge iff p < 0.05 / (number of pairs) and both carried by >= 3 reads; (4) a BLOCK = a connected component of >= 2 variants (a lone significant column, like a sequencing hotspot, is not a block); a read carries a block if it carries at least half of the block's variants that it covers; (5) split on the block whose smaller side is largest, if both sides have >= 3 reads; each side gets its own abPOA consensus and the procedure repeats (depth <= 6) until no block remains. The leaves, with their reads, are the candidates. Parameters: alpha 0.05, 15-column end window, 3 reads (the frozen minimum cluster), 6 reads to split, depth 6.
**Synthetic world (seed 20261014; `synth_world.py partition_specs`; 12 families per scenario; reads: end jitter 3, no artifact, 100 per copy).** All scenarios: two surviving copies A (2% SNPs from the template) and B (3%) and the erased copy or copies; divergence of the erased copy from A about 2%. CTL: one erased copy E (the template), one isoform. ISO: E has a second isoform (internal exon 3 skipped), reads 67 : 33. SIB: two erased copies E1 (the template) and E2 (0.5% SNPs from it), both erased, 100 reads each (0.5% < the edge cutoff, so the frozen clustering joins them).
**Scoring (exact truth).** Expected transcripts: each (erased copy, isoform) with >= 5 reads in the net. A candidate RECOVERS a transcript T iff its core identity (end gaps excluded, edlib NW) to T is >= 0.999 and not lower than to the other expected transcripts of the family (a blend does not count). A candidate is SPURIOUS if no expected transcript reaches 0.999, REDUNDANT if it duplicates an already recovered one. Baseline = the frozen consensus per cluster.
**Bars (set now).** P1 CTL: PART splits at most 1 of the 12 control clusters and recovers every transcript the baseline recovers. P2 ISO: baseline recovers <= 60% of the expected transcripts (the minority isoform is lost, prediction), PART >= 80%. P3 SIB: baseline <= 60%, PART >= 80%. P4: spurious + redundant candidates <= 10% of PART's candidates over the three scenarios, and leaf purity (reads of a leaf that belong to its best-matching transcript) >= 95%.
**Real dev A, descriptive (not judged; H after, not blind):** clusters split, candidates added, and whether the best candidate against the erased copy raises the gate-aware identity x coverage >= 0.999 count (65 -> 80 with the metric alone) for the clusters with > 200 edits.
**Decision.** PART is adopted as an opt-in post-clustering step iff P1 to P4 hold; otherwise not adopted, with the failure.

### Amendment 18 result (partition world seed 20261014; 36 families; 43 frozen clusters; 14 split; 57 candidates)

| scenario | expected transcripts | baseline (frozen consensus) recovers | PART recovers | clusters split / clusters | spurious + redundant candidates | leaf purity |
|---|---|---|---|---|---|---|
| CTL | 12 | 12 | 12 | 0 / 12 | 0 + 0 of 12 | 1.000 |
| ISO | 24 | 12 (50%) | **20 (83%)** | 9 / 12 | 0 + 1 of 21 | **0.849** |
| SIB | 24 | **19 (79%)** | **23 (96%)** | 5 / 19 | 0 + 1 of 24 | 0.958 |

**P1 met** (no control split, nothing lost). **P2 met** (baseline 50%, PART 83%). **P3 NOT met as written:** PART reaches 96% (bar >= 80%) but the baseline clause fails, baseline 79% not <= 60%, because the frozen clustering already separated 7 of the 12 sibling pairs. **P4 spurious + redundant met** (2 of 57 = 3.5%); **leaf purity NOT met for ISO** (0.849; bar 95%). **Decision: PART is not adopted by the registered rule.**
Cause of the ISO shortfall (3 of 12 clusters not split, all with a skipped exon of 271 to 322 bp; the 9 split include skipped exons of 93 to 310 bp): `minimap2 -x map-hifi` does not bridge a long deletion, it soft-clips the first part of the read (an iso1 read starts at consensus position 604 with a 286-base soft clip) so the skipped exon shows up as a change of coverage, not as a deletion run, and no variant column is found. This is an alignment limitation of the implementation, not of the partition logic.

---

## Amendment 19 (2026-10-08, after the Amendment 18 numbers, before any new number): global alignment for PART, fresh world

**Change.** PART aligns each read to the candidate with a semi-global edit-distance alignment (edlib `HW`: the read aligned in full inside the consensus, `task=path`, no distance cap), which keeps a long deletion as a deletion run, instead of `minimap2 -x map-hifi`. Reads are in transcript orientation (as the frozen pipeline's consensus is). Nothing else in PART changes (same variants, test, blocks, split rule, parameters).
**Test.** A fresh partition world, seed 20261015 (never run), E4 reads, 12 families per scenario. The bars of Amendment 18 are kept except that the baseline clause of P3 (a prediction about the frozen pipeline, not about PART) is replaced by PART >= baseline in every scenario: P1 (CTL: at most 1 of 12 clusters split, nothing lost); P2 (ISO: PART recovers >= 80%, and more than the baseline); P3 (SIB: PART >= 80%, and not fewer than the baseline); P4 (spurious + redundant <= 10% of PART's candidates over the three scenarios; leaf purity >= 95% in each scenario).
**Decision.** Adopted as an opt-in post-clustering step iff P1 to P4 hold on the fresh world. Real dev A is reported, not judged.

### Amendment 19 result (fresh world seed 20261015; 47 frozen clusters, 13 split, 61 candidates; edlib alignment)

| scenario | expected | baseline recovers | PART recovers | split / clusters | spurious + redundant | purity |
|---|---|---|---|---|---|---|
| CTL | 12 | 12 | 12 | 0 / 12 | 0 + 0 | 1.000 |
| ISO | 24 | 12 | **21 (88%)** | 12 / 12 | 0 + **4** | **0.829** |
| SIB | 24 | 20 | **21 (88%)** | 1 / 23 | 0 + 3 | **0.875** |

P1, P2 and the recovery clauses of P3 met; **P4 NOT met** (spurious + redundant 7 of 61 = 11.5% against 10%; leaf purity 0.829 and 0.875 against 95%). **Decision: not adopted.** Diagnosis (unit test, not a result): a unit-cost edit-distance alignment prefers a messy alignment of 40 mismatches to a clean 100-base deletion, so the isoform shows up through noise correlated across the reads of one isoform, not as a deletion run; the split works by accident and assigns reads loosely, hence the purity.

---

## Amendment 20 (2026-10-08, after the Amendment 19 numbers, before any new number): spliced alignment, N counted as a deletion

**Change.** PART aligns with `minimap2 -ax splice:hq -uf --eqx` (checked on the failing cluster of Amendment 18: all 33 exon-skipping reads come out as ONE alignment with a > 100 bp N gap, where `map-hifi` clips) and counts an N gap as a deletion run in both the pileup and the read-by-variant table (for PART only; the Amendment 9 polish is unchanged). Nothing else changes. **Test:** a fresh world, seed 20261016, same bars as Amendment 19 (P1 to P4). **Decision:** adopted as an opt-in post-clustering step iff P1 to P4 hold.

### Correction to the Amendment 18, 19 and 20 numbers (2026-10-08): a scoring bug in `parteval`, found by cross-checking the leaf composition

After the Amendment 20 run the leaf composition (`partition.json`) was perfect (all 12 ISO clusters split into a pure iso0 leaf of 67 reads and a pure iso1 leaf of 33) while the scorer reported purity 0.749 and 5 redundant candidates. The scorer took, for each candidate and expected transcript, the MAXIMUM core identity over the two orientations. A core identity is computed between the first and last run of 8 matches, so a garbage alignment (the wrong orientation, or an isoform against the other isoform) with one lucky run scores 1.0 over a tiny core. Fixed: the orientation is chosen by the GLOBAL identity and a core that covers less than half of the shorter sequence is `None` (tests added). Re-scoring the stored outputs of the three worlds (the partition output itself is unchanged):

| world (alignment) | CTL split / recovers | ISO baseline -> PART | SIB baseline -> PART | spurious + redundant of PART's candidates | purity CTL / ISO / SIB |
|---|---|---|---|---|---|
| 20261014 (Amendment 18, `map-hifi`) | 0 / 12 of 12 | 12 -> 21 of 24 (9 of 12 split) | 19 -> 24 | 0 of 57 | 1.00 / 0.905 / 1.00 |
| 20261015 (Amendment 19, edlib) | 0 / 12 of 12 | 12 -> 24 | 23 -> 24 | 1 of 61 | 1.00 / 0.994 / 1.00 |
| **20261016 (Amendment 20, spliced, the registered test)** | **0 / 12 of 12** | **12 -> 24 of 24 (12 of 12 split)** | 24 -> 24 | **0 of 60** | **1.00 / 1.00 / 1.00** |

**Verdicts with the corrected scorer.** Amendment 18 stays NOT adopted (P3 baseline clause: baseline 79% not <= 60%; ISO purity 0.905; `map-hifi` clips long deletions). Amendment 19 would have met P1 to P4 (the original "not adopted" came from the scoring bug; edlib's unit-cost alignment is still a poor aligner for long gaps, and its ISO split worked through correlated noise). **Amendment 20 meets P1 to P4 on the fresh world: PART is adopted as an opt-in post-clustering step** (`partition.py`, `PART_ALIGN=splice`). In the Amendment 20 world the frozen clustering had already separated every sibling pair (24 of 24 at baseline), so that world does not test the sibling gain; the 20261014 and 20261015 worlds, where the frozen edge rule joined 5 and 4 sibling pairs, show it (19 -> 24 and 23 -> 24 with 5 and 1 splits).

### Amendment 18/20 real beds, descriptive (PART with `splice:hq -uf`, N as deletion; candidates aligned to the unmasked genome; the "best candidate" picks, per cluster, the candidate that matches the erased copy best, so it is an oracle upper bound, not a pipeline output)

| | dev A (212 clusters on the erased copy) | H (40; not blind) |
|---|---|---|
| clusters split | 72 (34%); leaves per split cluster 2: 33, 3: 9, 4: 8, 5: 9, 6: 6, 7: 5, 12: 1, 14: 1 | 19 (48%); up to 12 leaves |
| gate-aware identity x coverage >= 0.999: frozen consensus -> best candidate | 80 -> **106** | 17 -> **26** |
| the split clusters only | 40 -> 66 | 9 -> 18 |
| edit distance to the erased copy, split clusters | 579 -> 239 | 195 -> 10 |
| the clusters with > 200 edits (the original motivation) | 22 of 212, **split 0**, best candidate at 0.999: 0 | 3, split 0, 0 |

PART finds structure in a third to a half of the real clusters and, taking the best candidate, raises the count at the registered threshold by 26 (A) and 9 (H); the cluster-level evidence is a real mixture of members. It does NOT touch the 22 clusters with > 200 edits: their reads are not a mixture separable by linked variant columns (they are far from the assembled copy as a whole, probably reads that belong to another copy or isoform structure than the one the genome label points at). Unknown on real data: how many of the 417 candidates are real transcripts (no truth), and whether clusters with 5 to 14 leaves are over-split (the synthetic control splits 0 of 12, but real reads carry correlated errors and read-end heterogeneity the simulation lacks).

---

## Amendment 21 (2026-10-08, before the H number): the 22 dev-A clusters with more than 200 edits, and a truth-free flag for them

**What dev A shows (diagnosis, looked at before this amendment).** The 22 clusters are tiny (3 to 12 reads, median 5) and their reads are individually fine: baseline alignment to the unmasked genome has median divergence 0.0013 and MAPQ 60, as for all other clusters. What is broken is the consensus: the median read-to-consensus divergence (reads aligned to their own consensus, `splice:hq -uf`) is 0.016 (quartiles 0.009, 0.039) against 0.0015 for the other small clusters and 0.0014 for the large ones. The reads of such a cluster differ strongly in length and span (e.g. 667, 3,657 and 4,508 bases, genome starts up to 41 kb apart): the edge rule (divergence <= 0.00958 and a block covering half of the shorter read) lets a short read bridge two long reads that do not overlap each other, the connected component is a chain of partial or alternative-structure reads, and a global POA consensus of such a chain is a blend, not a transcript.
**Flag Q (truth-free).** A cluster's consensus is UNSUPPORTED if the median divergence of its reads aligned to it exceeds 0.00958 (the allele cutoff the edge rule already uses; no new constant). On dev A: 24 of 212 clusters are flagged, 16 of the 22 with > 200 edits and 8 others (sizes 3 to 75, edits 1 to 151); the 6 unflagged big-edit clusters have a consensus that explains its reads (divergence 0.001 to 0.009) but differs from the erased copy as a whole, which Q cannot see.
**Test on bed H (not used for the diagnosis).** Bars: Q flags >= 2 of the 3 H clusters with > 200 edits and at most 15% of the other H clusters (dev A: 4%). Reported: the flagged clusters' sizes.
**Decision.** Q is adopted as a quality flag on the cluster consensus ("consensus unsupported") iff both bars hold; it changes no consensus.

### Amendment 21 result (2026-10-08; `consensus_support.py`)

| | dev A (design) | H (test) |
|---|---|---|
| clusters on the erased copy | 212 | 40 |
| flagged unsupported | 24 (sizes 3 to 75; 17 of them have 3 to 6 reads) | 3 (sizes 3, 8, 14) |
| of the clusters with > 200 edits flagged | 16 of 22 | **2 of 3** (bar >= 2) |
| other clusters flagged | 8 of 190 (4.2%) | **1 of 37 (2.7%)** (bar <= 15%) |

**Both bars met: Q is adopted as a quality flag ("consensus unsupported"); it changes no consensus.** What is explained: 16 of the 22 big-edit clusters on dev A (and 2 of 3 on H) are chains of 3 to 12 reads of very different spans whose POA consensus does not explain its own reads (median read divergence 0.016 against 0.0015). What is not: 6 big-edit clusters on dev A (and 1 on H) whose consensus explains its reads (divergence 0.001 to 0.009, 3 to 12 reads) but differs from the erased copy by 238 to 614 edits; their reads align to the unmasked genome individually at median divergence 0.0013, so the cause (reads from another copy or structure than the one the genome label points at, or a consensus that blends two structures the reads cannot expose) is not identified. In read terms the 22 clusters are 126 of about 59,000 reads (0.2%); in copy terms 10% of the clusters.

---

## Amendment 22 (2026-10-08, before any number from the new rule): chain-aware clustering, the edge must be a proper overlap

**What the pairs show (dev A, looked at before this amendment).** In the 16 flagged big-edit clusters, aligning the reads against each other pairwise (`splice:hq -uf`, the frozen edge rule applied to each pair) confirms only 104 of 216 pairs; 5 of the 16 clusters are not even connected by their own pairwise graph (3 have three isolated reads). The clusters were built from `map-hifi` read-to-seed edges, and the frozen edge accepts an alignment whose block covers only half of the shorter read (`block >= 0.5 x shorter`), wherever it sits: two reads that share a middle segment but continue differently at both ends are joined, and three such reads form a chain.
**Rule.** The frozen edge (gap-compressed divergence <= 0.00958 and block >= half of the shorter read) PLUS a proper overlap: at each end of the alignment at least one of the two reads must be reached, i.e. min(unaligned remainder of the query at that end, unaligned remainder of the target at that end) <= 0.00958 x block length (about 29 bases for a 3 kb alignment, the size of read-end ragged edges and of the 5' G run). A suffix-prefix overlap or a containment passes; an alignment that stops in the interior of both reads (a shared middle) does not. For a reverse-strand alignment the query start is paired with the target end. Nothing else changes (seed rounds, union-find components, minimum cluster of 3).
**Test.** Re-cluster dev A (design) and H (test) from scratch with the rule (`run_bed.py --proper --tag chain`), abPOA consensus, the Amendment 21 flag and the registered gate-aware metric on the clusters that overlap an erased copy. Bars, on H (test): (B1) reads of the erased copies clustered lose at most 3 points of coverage (frozen 95.7%); (B2) cluster purity by family stays 1.000; (B3) clusters with an unsupported consensus (Amendment 21) <= 1 of the 3 frozen; (B4) the gate-aware identity x coverage >= 0.999 count is not lower than the frozen 17, and the clusters with > 200 edits are <= 1 of 3; (B5) clusters on the erased copy do not fall by more than 5% (frozen 40). Dev A (design, direction only): unsupported <= 12 of 24, > 200 edits <= 11 of 22, fidelity count >= 80, coverage within 3 points.
**Decision.** Adopted as the cluster edge rule (opt-in flag `--proper`) iff B1 to B5 hold on H; otherwise not adopted with the failure.

### Amendment 22 result (2026-10-08; `run_bed.py --proper --tag chain`, `eval_clusters.py`)

| | H frozen | **H chain (test)** | dev A frozen | dev A chain (design) |
|---|---|---|---|---|
| clusters | 47 | 54 | 221 | 312 |
| erased-copy reads clustered | 95.7% | **94.0%** | 96.4% | 95.6% |
| purity by family | 1.0000 | **1.0000** | 1.0000 | 1.0000 (1 background read in a family cluster) |
| clusters on the erased copy | 40 | **49** | 212 | 304 |
| unsupported consensus (Amendment 21) among them | 3 | **1** | 24 | 22 |
| gate-aware identity x coverage >= 0.999 | 17 (43%) | **26 (53%)** | 80 (38%) | 143 (47%) |
| clusters with > 200 edits | 3 | **2** | 22 | 27 |

B1 met (coverage -1.7 points), B2 met, B3 met (1 of 3), B5 met (+23%), **B4 NOT met: the fidelity count rises (17 to 26) but the clusters with > 200 edits are 2, bar <= 1.** Dev A direction bars fail too (unsupported 22 against <= 12; > 200 edits 27 against <= 11). **Decision: not adopted by the registered rule.** What the numbers say: the proper-overlap edge makes more, smaller clusters (+23% on H, +43% on A) and the share of them that reaches 0.999 on the erased copy rises about 10 points on both beds, but the chain survives where a short read is CONTAINED in two long reads that disagree with each other: each containment is a proper overlap, so A-B and B-C are edges although A and C are not compatible. A pairwise edge rule cannot see that; it needs the transitive check (Amendment 23).

---

## Amendment 23 (2026-10-08, after the Amendment 22 numbers, before any new number): star clustering inside each component

**Change.** Inside every component of the proper-overlap clustering with at most 60 reads: all-vs-all `minimap2 -x map-hifi -c`, a pair is COMPATIBLE iff it passes the Amendment 22 edge (divergence <= 0.00958, block >= half the shorter read, proper overlap); then LONGEST-FIRST STAR clustering: the longest unassigned read is the representative, its cluster is itself plus the unassigned reads compatible with it, repeat; clusters of fewer than 3 reads are dropped (their reads leave the clustering). A read joins only a cluster whose representative it is compatible with, never through another member, so a short read contained in two incompatible long reads goes to the longer one. Components with more than 60 reads keep the Amendment 22 result.
**Test (a, blind, exact truth): a synthetic CHAIN world, seed 20261017.** 24 families: 12 CHAIN (erased copy E with two equally expressed isoforms, iso0 full and iso1 skipping internal exon 4; 16 reads per copy, 40% of the reads 5'-truncated by up to 70%, so that fragments starting after exon 4 are contained in both isoforms) and 12 CTL (one isoform, same depth and truncation). Frozen clustering, Amendment 22 clustering and Amendment 23 clustering are scored the same way as Amendment 18: expected transcripts = (erased copy, isoform) with >= 5 reads in the net; a candidate recovers a transcript if its core identity (orientation by global identity) is >= 0.999 and not lower than to the other expected transcript; purity = reads of a cluster that belong to its best-matching transcript. **Bars:** CHAIN: Amendment 23 recovers >= 80% of the expected transcripts with purity >= 90%, and more than the frozen clustering; CTL: Amendment 23 recovers every transcript the frozen clustering recovers, and loses at most 5 points of read coverage; coverage over both scenarios >= 85%.
**(b) Real dev A and H, descriptive (both used for Amendment 22):** counts as in Amendment 22.
**Decision.** Adopted as the cluster rule (opt-in) iff the bars of (a) hold; (b) is reported.

### Amendment 23 result (chain world seed 20261017; 12 CTL + 12 CHAIN families; frozen clustering, proper-overlap clustering, star clustering)

| | CHAIN expected / recovered | CHAIN reads clustered, purity | CTL expected / recovered | CTL reads clustered |
|---|---|---|---|---|
| frozen | 24 / 12 | 191 of 191, 0.50 | 12 / 12 | 192 of 192 |
| Amendment 22 (proper-overlap edge) | 24 / 12 | 191 of 191, 0.50 | 12 / 12 | 192 of 192 |
| **Amendment 23 (star step)** | 24 / **12 (50%)** | 182 of 191 (95%), **0.52** | 12 / 12 | 192 of 192 (100%) |

**Bar NOT met (CHAIN 50% against >= 80%, purity 0.52 against 90%): not adopted.** The CTL family is untouched (nothing lost, no coverage lost), the star step removes only 9 reads. Diagnosis: a full-length iso0 read and a full-length iso1 read (internal exon 4 skipped) align to each other in ONE alignment with a single 256 to 261 base deletion (block 2,122 of 2,124 bases), whose gap-compressed divergence `de` is 0.003: the pair passes the frozen rule, the proper-overlap test and the star compatibility, so the two isoforms stay together. The gap-compressed divergence hides exon-scale deletions by design.

---

## Amendment 24 (2026-10-08, after the Amendment 23 numbers, before any new number): compatibility counts every gap

**Change.** In the star step the pair's divergence is the full edit distance over the alignment block, NM / block (every inserted or deleted base counts), instead of the gap-compressed `de`; the pair is compatible iff NM / block <= 0.00958 (the same constant), the block covers half of the shorter read and the overlap is proper. A 260-base deletion in a 2,100-base block is 12%, far above the cutoff; HiFi indel errors (1 to 2 bases) are not. Nothing else changes (the seed rounds and the Amendment 22 edge keep `de`; the star step is applied to components of at most 60 reads).
**Test.** A fresh chain world, seed 20261018, same design and the same bars as Amendment 23: CHAIN: star clustering recovers >= 80% of the expected transcripts with purity >= 90% and more than the frozen clustering; CTL: recovers every transcript the frozen clustering recovers and loses at most 5 points of read coverage; coverage over both scenarios >= 85%. Real dev A and H reported, not judged (both used before).
**Decision.** Adopted as the cluster rule (opt-in) iff the bars hold.

### Amendment 24 result (fresh chain world seed 20261018; real beds descriptive)

| chain world | CHAIN expected / recovered | CHAIN reads clustered, purity | CTL expected / recovered | CTL reads clustered |
|---|---|---|---|---|
| frozen | 24 / 12 | 192 of 192, 0.50 | 12 / 12 | 192 of 192 |
| **proper-overlap edge + star step (every gap counted)** | 24 / **24 (100%)** | 191 of 192 (99%), **0.91** | 12 / 12 | 192 of 192 (100%) |

**All bars met** (CHAIN 100% against >= 80% and more than the frozen 50%; purity 0.91 against >= 90%; CTL nothing lost; coverage 99 to 100% against >= 85%). **Adopted as the cluster rule, opt-in:** `run_bed.py --proper` then `run_chain.py` (star step, components of at most 60 reads). Purity 0.91 is the ambiguity of contained fragments, whose bases are identical in both isoforms and which go to the longer representative.
Real beds, descriptive (both used before; frozen -> proper edge -> + star step):

| | H frozen | H proper | **H + star** | dev A frozen | dev A proper | **dev A + star** |
|---|---|---|---|---|---|---|
| clusters on the erased copy | 40 | 49 | **54** | 212 | 304 | **298** |
| erased-copy reads clustered | 95.7% | 94.0% | **92.8%** | 96.4% | 95.6% | **94.9%** |
| purity by family | 1.0000 | 1.0000 | **1.0000** | 1.0000 | 1.0000 | **1.0000** |
| unsupported consensus (Amendment 21) | 3 | 1 | **0** | 24 | 22 | **8** |
| clusters with > 200 edits | 3 | 2 | **0** | 22 | 27 | **7** |
| gate-aware identity x coverage >= 0.999 | 17 (43%) | 26 (53%) | **29 (54%)** | 80 (38%) | 143 (47%) | **156 (52%)** |

The chain problem is largely resolved on the real beds: the clusters with more than 200 edits fall from 22 to 7 (dev A) and from 3 to 0 (H), the unsupported consensus sequences from 24 to 8 and from 3 to 0, and the share of erased-copy clusters at 0.999 rises 12 to 14 points; the cost is 1.5 to 2.9 points of read coverage (reads of unsupported chains leave the clustering) and +40% clusters. The 7 dev-A clusters still above 200 edits are the previously unexplained kind (consensus fits its reads, not the erased copy).

---

## Errata and limitations after the independent review (2026-10-08)

Four fresh reviewers (a claims-vs-data audit and three code/method reviews of Amendments 9 to 24) worked from a frozen read-only copy. What they found, what was fixed, and what was not. Result logs for the print-only analyses are in `docs/unmapped_rescue_results/`.
**Fixed in code (tests added):**
1. *Star step all-vs-all* (Amendment 23/24): `chain.py` ran minimap2 with `-N 20` and the default `-p 0.8`, which drops edges in components of 20 or more reads and fragmented ordinary clusters (a depth-50 control world: 31 clusters instead of 12, 92% of reads, registered control bar failed). Now `-N 500 -p 0.1` plus an integration test. The registered test (depth 16, components below 20 reads) was unaffected; on fresh worlds with the fix, depth 16, 22 and 50: isoform transcripts recovered 24 of 24 each, controls intact (100% of reads, 12 clusters), purity 0.89 to 0.93. Real beds re-run (`chain3`): H 58 clusters, 93.0% of erased-copy reads clustered, 55 on the erased copy, unsupported 0, clusters over 200 edits 0, gate-aware 0.999: 31; dev A 304, 95.0%, 294, unsupported 14, over 200 edits 13, 152. **The earlier dev-A figures (22 -> 7, 24 -> 8) came from the defective step; the corrected ones are 22 -> 13 and 24 -> 14.**
2. *Partition scorer* (`parteval`): the first fix still accepted fragments (a 3' half of a transcript scored 1.0 against it and, by tie order, against the other isoform). Now `ends.recovered_identity` requires the core to cover 90% of the transcript and of the candidate, ties are deterministic; the three partition worlds re-scored give the SAME tables (leaf composition had been checked independently), so verdicts stand.
3. *T4 search limit* (a templated run of 7 or more G was cut back to 6): the limit now follows the consensus run; test added.
4. *Unit tests*: `test_partition.py` built its fixture from a fresh `Random(7)` per character (a homopolymer); fixed. The edlib test asserted the opposite of the documented limitation (unit-cost edlib does not keep a 100-base deletion); it now documents it.
5. *Accounting of templated bases lost by T3* (Amendment 15): a clip that began with G but was not pure G counted as zero unaligned G. Corrected: at most 1 base lost in 157 of 197 trims (80%; H 28 of 37), still under the 90% bar; trims with no unaligned G at all are 61 (not 105). Verdict unchanged.
6. *Start-tolerant exactness* (Amendment 16): the modal-read-start convention scores a consensus that begins exactly at the true transcript as +1 when the jitter is one base (the untouched control is 21 of 24 exact). With a start-tolerant measure T4 is exact in 23 of 24 clusters, removes templated bases in 0, leaves an artifact G in 1 (this world and a fresh one, seed 20261301: 23, 0, 1; T3 11/4/9 and 15/1/8). T4 meets the synthetic bars (2) to (5); it remains not adopted because the dev-A bar fails (under-trim 36%; 23 of the 36 are clusters of fewer than 5 reads where T4 does nothing).
**Corrections to claims:**
7. *Amendment 17* was not tested on data "not used for this rule": the worlds had been scored for T3/T4, whose trim is the same operation, and the guard against false passes cannot fire on error-free consensus sequences. Adopted as a reporting CONVENTION; it changes the registered O3 flag (LRPAP1 c01 and the c00 allele cluster flip to flagged). A truly fresh world (seed 20261301) and the reviewer's (20261101) agree: 19 and 14 -> 24 of 24, controls identical, no false passes.
8. *Amendment 24 purity bar* (0.90): across eight chain worlds the purity is 0.85 to 0.93 (below 0.90 in four with the defective step; 0.89 to 0.93 with the fixed one). The rule is judged on the recovery and control bars; purity is limited by fragments that are identical in both isoforms. "All bars met" holds for the registered seed only.
9. *Remaining over-200-edit clusters* (Amendment 21/24): 71 of the 126 reads of the 22 clusters leave the clustering (dev A), 14 of 14 on H, so part of the improvement is removal; some remaining ones are unsupported (not "consensus fits reads"). Six of the 22 dev-A clusters hold reads of BOTH orientations (the pool graph ignores strand and the consensus is not oriented), an undisclosed mechanism for "consensus unsupported". The star step refines only components of at most 60 reads (4 to 8% of the reads); the mixed-isoform consensus in large clusters is untouched.
10. *Partition (Amendments 18 to 20)*: (a) validated only against independent errors; with errors correlated across a read subset (a hotspot of fixed columns in a quarter of the reads) PART split 10 of 12 control clusters in the reviewer's test, and two linked heterozygous SNPs split 5 of 12, so the real-bed split counts (34 to 48%, up to 14 leaves) are not evidence of distinct transcripts; (b) it is blind to an extra exon carried by a minority (an insertion is one variant, never a block) and its sensitivity floor depends on the variant length (a 280-base skip in 8 of 100 reads is missed); (c) the sibling gain (19 -> 24, 23 -> 24) comes from worlds that used the rejected alignments, and the adopted configuration separated 0.5% pairs in a forced-join test but only 10 of 12 at 0.2%; (d) the real-bed "80 -> 106" counts `best or frozen` (monotone by construction), 8 of the 26 gains are a shorter leaf scoring 1.0, and 28 of the 72 split dev-A clusters have two largest leaves of equal length differing by 1 to 15 substitutions (an allele signature in a diploid animal). "A real mixture of members" is an interpretation. PART stays opt-in and is not validated for correlated errors or alleles.
11. *Stale or wrong statements fixed in place:* "held-out H" (Amendment 2 already declared H not held out for the arms that followed; "not tuned on H" is what holds); "30 of 30 families" and "0 of 30" (26 clusters exist); "joined 5 and 4 pairs" (5 and 1); in Amendment 16 the j_hat tallies (6, 6, 6, 7 = 25) were not per erased copy; in Amendment 14 both failures were `too_long` (the `not_g` cluster passes); in Amendment 15 the design-check paragraph (39,351; 28,581) differs from the result file (39,403; 28,614) and "p = 1" for controls is the n < 20 sentinel; in Amendment 21 "17 of them" is 16 and the median divergence is 0.0176 (0.0094, 0.0419); in Amendments 18/20 "417 candidates" covers all 221 clusters (407 for the 212); in Amendment 5 GWFAM163 attracts six wrong clusters (125 wrong reads in 10 clusters) and two clusters are lost, the +243 holds; in section 4 "26 of 38" is 26 of 39. The real beds are real reads with SYNTHETIC erasure, and the library gate and artifact distribution are computed from H's own survivor reads.
**Known statistical looseness (no verdict changes):** the binomial gate tests p = 0.5, not 0.25; the polish test divides alpha by 3L but tests more alternatives per column (anti-conservative, so Amendment 9's refutation is unaffected); the artifact length distribution is measured from soft clips and undercounts artifact G's that match the genomic flank (a(0) 0.21 against 0.05 simulated), which biases the templated-count estimate upward.
**Not checkable from the frozen copy, now backed by logs or still open:** the real-bed trimming and metric analyses (logs saved), the margin sweeps on dev A, the O3 maternal study (no result files in the copy), Amendment 8's score-rule per-copy breakdown (stored for the divergence rule only).

---

## Amendment 25 (2026-10-08, before any number): specificity controls for PART and for the O3 flag

**Why.** The independent review showed that the synthetic controls so far have independent errors, haploid copies and a complete net, so they cannot expose the two ways the pipeline can call something that is not there: (i) PART splitting a cluster whose reads share a recurrent sequencing error (a HiFi homopolymer hotspot), and (ii) the O3 flag calling an allele of a present copy, or a family with nothing missing, an absent copy.
**Part A: a guard for PART (block purity).** A linked block is a real haplotype or isoform when reads carry (almost) all of its variants or (almost) none; a hotspot is carried partially and ALSO appears at a few percent in the reads that do not carry it. For a candidate block with m variants over n reads (carriers = at least half of the variants they cover): `Dev` = the number of (read, variant) pairs in which the read disagrees with its own group (a carrier lacking the variant, a non-carrier having it); `e_bg` = the non-consensus allele rate outside the variant columns (floor 1e-4); the block is REJECTED iff P(Binomial(number of pairs, e_bg) >= Dev) < 0.05 / (number of candidate blocks); the next block is tried. No other change to PART.
Bars: read-level stress harness (`stress_part.py`, real `splice:hq` alignment, 12 transcripts per setting, fresh seeds). (A1) Hotspot control (75 good reads with 3% error at F fixed columns, 25 low-quality reads with probability p): splits <= 1 of 12 for (F, p) = (20, 0.3), (20, 0.5), (40, 0.5); independent-error control 0 of 12. (A2) True structure retained: an exon skip in 33% of the reads, split in >= 11 of 12; two siblings 0.5% apart (6 or more differing columns), split in >= 11 of 12; two linked heterozygous SNPs, split (allelic heterogeneity is reported by PART, classified by Part B). (A3) The registered partition world is not harmed: a fresh world, seed 20261401 (`--partition-world`): CTL splits 0 of 12, ISO recovers >= 22 of 24.
**Part B: separating an allele from an absent copy in the O3 flag.** For a cluster consensus: R = identity x coverage (gate-aware) of its best hit on the reference genome, T = the same on the truth genome (the other haplotype assembly, or in a synthetic world the genome with every copy). Registered flag (unchanged): T >= 0.999 and R < 0.999. New classification: COPY = T >= 0.999 and R < 1 - 0.00958 (more than the allele cutoff from every reference locus); ALLELE = T >= 0.999 and 1 - 0.00958 <= R < 0.999; PRESENT = R >= 0.999. Only COPY is an O3 call.
Bars: synthetic specificity world (seed 20261402, `--spec-world`): 12 NULL families (three copies, all present), 30 ALLELE families (A and B present; the individual is heterozygous at A, the second allele A' differs from A by 0.25%, 0.5%, 1%, 2%, 4%, six families each, A' absent from the reference), 12 ERASED families (copy E at 2% from its nearest survivor, the positive control). (B1) NULL: no cluster is classified COPY and none is flagged by the registered rule. (B2) ALLELE: for A' at 0.25% and 0.5%, 0 COPY calls (the reads either do not enter the net or form clusters classified ALLELE); the table for 1%, 2% and 4% is REPORTED, not judged (an allele more than 1% from its other haplotype is indistinguishable by sequence from a paralog; it is what RNA alone cannot decide). (B3) ERASED: >= 95% of the erased-copy clusters are COPY. (B4) LRPAP1, a sanity check on known data: the c00 allele cluster is ALLELE, the c01 clusters are COPY.
**Decision.** The guard is adopted into PART iff A1 to A3 hold; the three-way classification is adopted for the O3 call iff B1 to B4 hold. Failures are reported as such.

### Amendment 25 result (2026-10-08; `stress_part.py`, `flagmetric.o3_class`, `run_synth.py o3class`; logs in `docs/unmapped_rescue_results/`)

**Part A, the PART block-purity guard (read-level stress, 12 transcripts per setting, real `splice:hq` alignment).**

| setting | guard OFF (seed 5001) | **guard ON, seed 5001** | **guard ON, seed 5002** |
|---|---|---|---|
| independent-error control | 0 of 12 split | **0** | **0** |
| hotspot, F = 20 columns, p = 0.3 | 10 of 12 | **1** | **0** |
| hotspot, F = 20, p = 0.5 | 12 of 12 | **0** | **0** |
| hotspot, F = 40, p = 0.5 | 12 of 12 | **0** | **0** |
| exon skip in 33% of the reads | 12 of 12 | **12** | **12** |
| two siblings 0.5% apart | 12 of 12 | **12** | **12** |
| two linked heterozygous SNPs | 6 of 12 | 6 | 11 |

A1 met (hotspots <= 1 of 12 in every setting, control 0). A2 met for the skip and the siblings (12 of 12; the bar was >= 11); the heterozygous case is not a numeric bar and PART finds two linked SNPs in 6 and 11 of 12 transcripts (a 2-column block is at the edge of its power). A3 met on a fresh partition world (seed 20261401): CTL 0 of 12 split, ISO 24 of 24 recovered (12 of 12 split), SIB 24 of 24, purity 1.00, no spurious or redundant candidate. **The guard is adopted into PART.** Effect on the real beds (descriptive, same reads): dev A 63 of 221 clusters split (73 before), 330 candidates (417), at most 7 leaves (14), best-candidate gate-aware 0.999 count 103 (106); H 17 of 47 (20), 88 candidates (118), at most 10 leaves (12), 24 (26). The guard removes the many-leaf splits; the remaining splits are still not validated against alleles (see the Errata). Reproducibility note: two H clusters' identical consensus sequences gave alignments that differ by one edit between two minimap2 batches, which moves the H frozen count 17 -> 16.
**Part B, the allele-vs-copy classification (specificity world seed 20261402; E4 reads; library gate closed in this world).**

| scenario | true kind | reads in net | clusters | COPY | ALLELE | registered flag |
|---|---|---|---|---|---|---|
| NULL (3 copies, nothing missing) | present | 0 | 0 | 0 | 0 | 0 |
| ALLELE 0.25% | allele | 0 / 480 | 0 | 0 | 0 | 0 |
| ALLELE 0.5% | allele | 7 / 480 | 1 | **0** | 1 | 1 |
| ALLELE 1% | allele | 408 / 480 | 6 | 4 | 2 | 6 |
| ALLELE 2% | allele | 480 / 480 | 6 | 6 | 0 | 6 |
| ALLELE 4% | allele | 480 / 480 | 6 | 6 | 0 | 6 |
| ERASED (positive control) | copy | 960 / 960 | 12 | **12** | 0 | 12 |

B1 met (no cluster at all from the 12 NULL families, so nothing is flagged), B2 met (0 COPY calls for alleles at 0.25% and 0.5%; at 0.5% the one cluster is ALLELE where the registered flag calls it a copy), B3 met (12 of 12 erased copies are COPY). **B4 met on LRPAP1:** the 21-read c00 cluster (mat best 0.9959) is ALLELE, the 283-read and 4-read c01 clusters (mat best 0.9899, 0.9895) are COPY; the 3-read cluster has no call. Reported, not judged: an allele 1% from its other haplotype is already called COPY in 4 of 6 families and 2% or more in all (6 of 6): by sequence alone RNA cannot tell an allele that far apart from a paralog, which is the limit of any RNA-only O3 call; the c01 consensus sits only 0.0005 below the allele cutoff on LRPAP1 (0.9899 against 0.99042). **The three-way classification is adopted for the O3 call; the registered flag (COPY or ALLELE) is kept as the screening statistic.**

---

## Amendment 26 (2026-10-08, before any classification): the three-way O3 call on real data with an independent truth

**Setting.** The maternal-reference study (`O3_MATERNAL_REFERENCE_2026-10-08.md`) has real KB3781 reads, the animal's own two haplotype assemblies and, for each judged copy, the identity of its nearest relative on the reference haplotype (section 2 tables). The reference for the call is one haplotype (R), the truth the other (T). Pool = every one of the 35,094 reads of the 34 families and LRPAP1 whose primary alignment on the reference haplotype is unmapped or has divergence > 0.00958; clusters by the adopted rule (`--proper` edge, star step with the fixed all-vs-all, minimum 3 reads); abPOA consensus; R and T = gate-aware identity x coverage of the consensus' best hit on each haplotype (`flagmetric.o3_class`, library gate computed from this dataset's own shared reads).
**Prediction, written from the study's tables before any consensus is classified.** With the MOTHER as reference the only judged copy with no relative within the allele cutoff is GWFAM175_B0 (nearest relative 0.9126); every other judged copy has one at 0.9911 to 1.0000 (GWFAM208_B0 1.0000, 214_B0 0.9995, 205_B0 0.9935, 64_B17 0.9999, 64_B11 0.9980, 227_B0 0.9911, LRPAP1_p12 0.9959). So COPY is expected for GWFAM175_B0 and, from the LRPAP1 example, for the LRPAP1 c01 cluster; ALLELE or PRESENT for the others. With the FATHER as reference all 13 judged copies have a relative at 0.9927 or closer (GWFAM382_B0 matches neither haplotype within 5%): no COPY call is expected.
**Bars (set now).** (R1) Mother as reference: at least one cluster of >= 20 reads whose majority label is GWFAM175_B0 is COPY. (R2) Mother as reference: no cluster whose majority label is one of GWFAM208_B0, 214_B0, 205_B0, 64_B17, 64_B11, 227_B0 or LRPAP1_p12 is COPY, and no cluster of 'shared' reads (reads that map to the same place on both haplotypes), LRPAP1 reads excluded, is COPY. (R3) Father as reference: zero COPY calls among the clusters of the judged-copy labels, and at most 1 COPY call in total. (R4) LRPAP1 in the same flow (mother as reference): the c01 clusters are COPY and the c00 cluster is ALLELE. Reported, not judged: every COPY call in each direction with its reads, labels and identities; the calls the registered flag would make; the number of ALLELE and PRESENT clusters; (for bed H, descriptive) the three-way classes of the erased-copy clusters against the masked reference.
**Decision.** The three-way call is reported as supported on real data iff R1 to R4 hold; each failure is reported as such.

### Amendment 26 result (2026-10-08; `o3_haplotype.py`, logs `docs/unmapped_rescue_results/o3_haplotype.{mat,pat}.eval.log`)

Pools: mother as reference 2,888 net reads (GWFAM175_B0 283 of its 304 reads, 'shared' 2,392, LRPAP1 included), 37 clusters of >= 3 reads; father as reference 2,824 net reads, 36 clusters. The library gate is open in both (G 20,247 against A, C, T 34 pure 5' clips).

| mother as reference, truth father | clusters | COPY | ALLELE | PRESENT | no call | registered flag |
|---|---|---|---|---|---|---|
| GWFAM175_B0 | 1 (269 reads) | **1** | 0 | 0 | 0 | 1 |
| LRPAP1 | 6 | 2 (the 282-read and the 4-read c01 clusters) | 0 | 0 | 4 | 2 |
| 'shared' reads | 28 | **5** (356, 63, 28, 8, 65 reads) | 6 | 0 | 17 | 11 |
| 'ambiguous' | 2 | 0 | 0 | 0 | 2 | 0 |

Father as reference: 36 clusters; COPY 10 (7 'shared' with 210, 79, 134, 480, 3, 4 and 101 reads; 2 'ambiguous'; 1 LRPAP1); the judged-copy labels GWFAM70_B4 ALLELE, GWFAM26_B2 and GWFAM382_B0 no call.
**Bars. R1 met** (the GWFAM175_B0 cluster is COPY, R = 0.9306). **R2 NOT met** (no judged-locus cluster is COPY, but 5 of the 28 'shared' clusters are). **R3 NOT met** (no COPY among the judged-copy labels, but 10 COPY calls in total against <= 1). **R4 partly:** both c01 clusters are COPY, but the c00 allele cluster of the earlier flow is not formed in this pool (the four small c00 clusters have T below 0.999: no call), so ALLELE is not shown. **The three-way call is NOT supported on real data by the registered bars.**
What the failures are: every one of the 'shared' COPY clusters has its reference hit and its truth hit on the SAME chromosome of the two haplotypes (chr5_mat 49.05 Mb with chr5_pat 41.47 Mb, chr8, chr11, chr2, chr1), i.e. it is a diverged counterpart at the corresponding locus: R = 0.973 to 0.990 (and 0.83 to 0.93 in the father's direction), T = 1.0. Real haplotypes of this animal differ by 1 to 7% at a minority of expressed loci (5 of 28 'shared' clusters with the mother as reference, 7 of 24 with the father), above the 0.958% allele cutoff, so a purely sequence-based call over-calls. The GWFAM175_B0 hit (chr5_mat 47.66 Mb) sits at the same offset to its pat locus (7.64 Mb) as the neighbouring 'shared' cluster (7.58 Mb): its nearest mother's relative is at the orthologous position too, 7% away. What separates the two kinds is not distance but whether the T locus is the reference locus' best match in return (one-to-one orthology), Amendment 27.

---

## Amendment 27 (2026-10-08, after the Amendment 26 numbers, before any reciprocal result): a COPY must have no reciprocal best match

**Rule.** For a COPY-called cluster, take the reference locus X of its best hit and build the X TRANSCRIPT, the reference segments the consensus alignment covers (exons joined, introns removed, from the `cg` CIGAR of the existing PAF). Map the X transcript to the TRUTH haplotype (`splice:hq -uf -N 5`, the same index). X is the ortholog of the truth locus Y0 iff the best record of the X transcript (by matches) lies on the same chromosome and overlaps Y0 (the truth hit of the consensus) by >= 50% of the shorter. Reciprocal COPY -> ALLELE ("diverged ortholog"); non-reciprocal -> COPY confirmed (the reference locus belongs to another truth locus, so Y0 has no counterpart). No new constant.
**Predictions (from the study's own tables, before the test).** Mother as reference: the GWFAM175_B0 cluster and the LRPAP1 c01 clusters are non-reciprocal (the study found that their nearest relatives are claimed by other copies); the 5 'shared' COPY clusters are reciprocal. Father as reference: the 7 'shared' COPY clusters are reciprocal.
**Bars.** (S1) mother as reference: GWFAM175_B0 and the 282-read LRPAP1 c01 cluster stay COPY. (S2) at least 4 of the 5 'shared' COPY clusters (mother) and at least 6 of the 7 (father) become ALLELE. (S3) after the rule, COPY calls from 'shared' clusters are <= 1 in each direction. Reported: every call, the Y1 locus found for each X transcript, the 'ambiguous' and LRPAP1 COPY clusters in the father's direction.
**Decision.** The reciprocal rule is adopted into the three-way call iff S1 to S3 hold.

### Amendment 27 result (2026-10-08; `reciprocal.py`, `o3_haplotype.py reciprocal`, logs `o3_haplotype.{mat,pat}.reciprocal.log`; `o3_bed.py`)

**Mother as reference** (8 COPY clusters): the GWFAM175_B0 cluster (269 reads) stays COPY (its mother's locus, chr5_mat 47.66 Mb, back-maps to a different father's locus, chr5_pat 39.97 Mb, at 0.932); the 282-read and the 4-read LRPAP1 c01 clusters stay COPY (the mother's locus is c00, whose father's counterpart is c00 on chr3, not c01 on chr12); **all 5 'shared' COPY clusters are reciprocal** (356, 63, 28, 8, 65 reads: each back-maps onto its own father's locus at 0.977 to 0.989) and become ALLELE (diverged orthologs). S1 met, S2 met (5 of 5), S3 met (0 'shared' COPY).
**Father as reference** (10 COPY clusters): 4 of the 7 'shared' COPY clusters are reciprocal (210, 79, 134, 101 reads) and become ALLELE, as does the 5-read LRPAP1 cluster; **3 'shared' clusters stay COPY** (the 480-read cluster, whose father's locus chr5_pat 39.97 Mb back-maps to a DIFFERENT mother's locus, chr5_mat 47.55 Mb, at identity 1.0000; two chr1 clusters of 3 and 4 reads) and the two tiny 'ambiguous' clusters stay COPY. **S2 NOT met (4 of 7 against >= 6), S3 NOT met (3 against <= 1). The reciprocal rule is NOT adopted by the registered bars.**
What the residual calls are: the 480-read cluster is the mirror image of the mother-direction GWFAM175_B0 call. In the mother's direction the father's copy at 40.03 Mb has no ortholog (its nearest mother's locus, 47.66 Mb, 93% identical, is claimed in return by a different father's locus); in the father's direction the mother's copy at 47.66 Mb has no ortholog (its nearest father's locus is 93% identical and claimed by the father's other copy). The two haplotypes of this animal carry two paralogs of the GWFAM175 tandem array that the other lacks: a real copy-number polymorphism, expressed in 269 and 480 reads. The study's synteny-based truth list (loci absent from the father/mother) lists the first and labels the second 'shared'; the rule found a true absent copy the registered truth misses, not a false call. The bars assumed the list was complete (the study had already said it is not). The two chr1 clusters (3 and 4 reads, R = 0.83) and the two 'ambiguous' clusters are too small to judge.
**Bed H (erasure, masked genome as reference; descriptive):** 55 clusters on an erased copy: 31 COPY (the same 31 that reach 0.999), 24 no call (not reconstructed at 0.999), 0 ALLELE; 3 background clusters, no call. By construction the net holds only reads more than the allele cutoff from the reference, so an erased copy whose nearest survivor is within the cutoff leaves no cluster to call.
**Verdicts.** The three-way sequence-only call (Amendment 25) is supported where truth is exact (synthetic) and on bed H, and over-calls on real haplotypes where 5 of 28 (mother) and 7 of 24 (father) expressed 'shared' clusters are diverged orthologs 1 to 7% apart. The reciprocal rule removes 9 of those 12 and leaves one genuine copy-number polymorphism and four tiny clusters; by the registered bars it is not adopted, and the result is best read as: the call needs one-to-one orthology (not a distance) to separate alleles from copies, and real haplotypes show both a high rate of diverged orthologs and copies the synteny truth list does not contain.

---

## Amendment 28 (2026-10-08, before any result): discovery of reference-absent sequence from the transcriptome, the 959 reads unmapped on the primary assembly

**Question.** The aim is sequence the reference genome does not contain but the transcriptome reaches, with no truth built in. Real case: KB3781 Iso-Seq, the 959 reads that do not align to the combined primary assembly (`o3_mat/reads/R_unm.fa`). Independent confirmation is used only to grade the result: the animal's maternal and paternal assemblies (132 of the 959 reads map on the mother's, 3 on the father's, 125 of them in one 71-kb stretch of mat chromosome 12; 824 map on neither).
**Pipeline (no truth, no haplotype).** Cluster the 959 reads with the adopted rule (`--proper` edge, star step with the fixed all-vs-all, minimum 3 reads), abPOA consensus per cluster (<= 100 longest reads), and align each consensus (`splice:hq -uf -N 5`, gate-aware identity x coverage, the library gate computed from this dataset's shared reads, open) to the primary (R), to the mother's (Tm) and to the father's (Tp) assembly. Classes: CONFIRMED = R below 1 - 0.00958 and max(Tm, Tp) >= 0.999 (absent from the primary, present in the animal's other assembly); NOVEL = R below 1 - 0.00958 and max(Tm, Tp) < 0.999 (a reconstructed sequence none of the three assemblies holds); ALLELE-LIKE = R in [1 - 0.00958, 0.999); PRESENT = R >= 0.999.
**Control (specificity).** 959 reads sampled (seed 1) from `o3_mat/reads/R_LRP.fa` (LRPAP1 reads; every one of them maps on the primary), same pipeline, same classes.
**Predictions.** The 125-read stretch forms one cluster that is CONFIRMED against the mother's assembly; most of the other 824 reads stay unclustered; a few small clusters are NOVEL.
**Bars (set now).** (D1) the largest cluster (>= 100 reads) is CONFIRMED. (D2) at least 80% of the 135 reads that map on the mother's or father's assembly lie in CONFIRMED clusters. (D3) control: no CONFIRMED cluster and at least 98% of its clusters PRESENT or ALLELE-LIKE (no NOVEL). Reported, not judged: the number, reads and lengths of the NOVEL clusters, the reads left unclustered, the clusters of 'unmapped' reads that are PRESENT or ALLELE-LIKE.
**Decision.** The discovery run is reported as working iff D1 to D3 hold; every NOVEL cluster is a candidate for DNA verification (consensus k-mers in the individual's genomic reads), not a call.

### Amendment 28 result (2026-10-08; `discover.py`; logs `docs/unmapped_rescue_results/discover.*.log`)

959 unmapped reads -> 3 clusters of >= 3 reads (113, 4 and 3 reads; 120 reads clustered, 839 unclustered). The 113-read consensus (5,694 bp) aligns to the MOTHER's assembly at chr12_mat 95.96-96.03 Mb with identity 1.0000 over 5,557 of 5,694 bases (no hit on the primary), gate-aware identity x coverage Tm = 0.9759: **CLASS NOVEL, so D1 NOT met and D2 NOT met (0 of 134 confirmable reads in CONFIRMED clusters; 113 are in a cluster).** Control (959 clean LRPAP1 reads): 14 clusters, 7 PRESENT (362 reads), 6 ALLELE-LIKE (34), **1 NOVEL (462 reads, 4,244 bp, R = 0.77)**, 0 CONFIRMED: **D3 NOT met** (13 of 14 = 93%).
Why. (1) The 137 unaligned 5' bases of the 113-read consensus are real sequence of the same locus, not an artifact: normal composition, and mapped on their own they fall in the mother's locus in two pieces (46 bases at 96,028,629 and 90 bases at 95,981,351): small terminal exons the spliced aligner did not join, so the registered best-record coverage undercounts. (2) The control NOVEL cluster is a 462-read component (above the 60-read limit of the star step) whose consensus is a blend of reads that overlap in different parts of a multi-copy locus; the Amendment 21 support test would reject it. So the classes missed their own aim through two defects of the definition, not of the pipeline. Registered as written, D1 to D3 fail.

---

## Amendment 29 (2026-10-08, after the Amendment 28 numbers, before any revised classification): terminal-exon rescue and a support requirement

**Changes to the discovery classification only (clustering, consensus and alignments are unchanged).** (a) TERMINAL-EXON RESCUE: the consensus bases left unaligned at either end of the best record (>= 20 bases) are mapped (`minimap2 -x sr -N 5`) to the genomic window [record start - 5 kb, record end + 5 kb] of the same assembly; pieces with identity >= 0.98 are added to the coverage (union of query intervals) and to the identity (matches over block); applied to all three assemblies; the 5' G rule of Amendment 17 is unchanged. (b) A cluster whose consensus is UNSUPPORTED (median divergence of its reads to it above 0.00958, Amendment 21) is class UNSUPPORTED and is never CONFIRMED or NOVEL.
**Predictions (from the diagnosis above).** The 113-read cluster becomes CONFIRMED (the two pieces account for the 137 bases); the control 462-read cluster becomes UNSUPPORTED; the 4- and 3-read clusters stay NOVEL or UNSUPPORTED.
**Bars.** D1: the largest cluster of the 959 reads is CONFIRMED. D2: >= 80% of the 134 confirmable reads are in CONFIRMED clusters (113 of 134 = 84% if D1 holds). D3: on the original control (seed 1, not blind now) and on a FRESH control (959 reads, seed 2, blind): no CONFIRMED cluster and no NOVEL cluster. Reported: the NOVEL and UNSUPPORTED clusters of the 959 reads with sizes and lengths.
**Decision.** The discovery run is reported as working iff D1 to D3 hold; every NOVEL cluster remains a candidate for DNA verification, not a call.
