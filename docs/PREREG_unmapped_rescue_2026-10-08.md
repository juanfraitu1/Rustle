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
