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
