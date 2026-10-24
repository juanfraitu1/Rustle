# Ledger: general (handoffs, pending lists, register drafts, overviews)

Merged 2026-10-24 from 14 dated docs, oldest first; each kept verbatim under its own heading (its headings one
level down). The originals are in `docs/archive/2026-10/merged/` (git history has them at their old paths).

## Contents

- REGISTER_DRAFTS_machine2.md
- REGISTER_DRAFTS_o3_maternal.md
- REGISTER_DRAFTS_unmapped_rescue.md
- PENDING_2026-09-23.md
- ADVISOR_HIERARCHY_2026-10-03.md
- HANDOFF_2026-10-04.md
- PENDING_2026-10-04.md
- TERMINOLOGY_FAMILY_SD_DUPLICON_EXPANSION_2026-10-07.md
- PENDING_2026-10-11.md
- OPEN_IMPLEMENTATIONS_PHASED_2026-10-12.md
- HANDOFF_laptop_to_main_2026-10-15.md
- O1O2O3_FAST_REALISTIC_2026-10-16.md
- O1O2O3_INTEGRATED_OVERVIEW_2026-10-16.md
- HANDOFF_cloud_2026-10-21.md


---

## REGISTER_DRAFTS_machine2.md

*(original: `docs/archive/2026-10/merged/REGISTER_DRAFTS_machine2.md`)*

## Register drafts from the machine2 clone (renumber on merge into docs/NEGATIVE_RESULTS_REGISTER.md)

| M2-1 | 2026-09-30 | Soto nesting | The 99.1% sequence-only nesting of Soto families does not depend on the edge set (A2) | ⚠ **It does: 440/444 (99.1%) on exon map-back edges, 346/444 (77.9%) on SEDEF-projected edges (`shared_exons_2334_finalv1_native.tsv`); the 98 extra splits are missing edges (e.g. AMY1/AMY2A ID_131), not conflicts.** Quote 99.1% with its edge set. `docs/archive/2026-09/SOTO_M2_MEETING_EVIDENCE_2026-09-30.md` §1 |
| M2-2 | 2026-09-30 | Soto table audit | 09-10 biotype audit (605 families, 217 all-pseudogene) reproduces from S1C (B1) | ⭐ **Reproduces exactly only when the 114 `Unassigned_*` singletons count as families; on the 491 real families: 183 (37.3%) all-pseudogene, 253 (51.5%) no protein-coding member, 57.8% pseudogene members; 0/491 mismatches vs S1C's own `No. Protein Coding`.** §3.5 |
| M2-3 | 2026-09-30 | Soto table audit | 09-10 fragment audit (147/420 families bundle a < 20% member with a ≥ 80% one) reproduces (B2) | ⛔ **Does not reproduce (script gone, 420 denominator not rebuildable): 491/491 families have ≥ 2 members with exons; 154 (31.4%) hold a < 20% member, 45 (9.2%) also a second ≥ 80% member.** Quote 154/491. §3.6 |
| M2-4 | 2026-09-30 | Soto finer than ours | Soto's finer splits are backed by stable copy-number gaps (B3) | ⛔ **Not always: 40/47 single-exclusive-member Soto families share ≥ 98% exons with another family's member (NPIPB3/B4/B5 = ID_151/152/153, 98.6% identical); their famCN gaps (26.5/29.1/22.2) are not reproduced by our 268-sample re-measurement on the same interval (31.7/36.9/33.4, a different grouping); at famCN 20-35, 21% of genes differ by ≥ 2.6 copies between the two measurements.** §3.1-3.2 |
| M2-5 | 2026-09-30 | Soto per-family | Every Soto family can be labelled match / Soto smaller / we miss (our rule vs Soto's table) against our 0.9277 run | ⭐ **match 411 (= the ladder's exact, asserted), Soto smaller 37, we miss · our rule 24 (all our copy numbers), we miss · Soto's table 5 (own rule + own CN also split: ID_28/113/144/270; hand merge ID_482), both 14. Of 43 we split, 36 ours, 7 Soto's.** `soto_m2_families.py`, `docs/SOTO_M2_FAMILIES_2026-09-30.tsv` |
| M2-6 | 2026-09-30 | Soto biotypes | Soto's families survive removing pseudogenes / lncRNAs | ⛔ **No pseudogenes: 974/2,334 genes, 313/491 families gone (fewer than 2 genes or no copy-number backbone left); ours 143 exact of 178 (ARI .950); 8 hold together only through pseudogenes. No lncRNAs: 38 gone. Both: 352 gone.** Soto's released code keeps exons of every biotype (grep commented out) and pairs if either gene is coding or an unprocessed pseudogene. |
| M2-7 | 2026-09-30 | Soto CN cut | Loosening our copy-number cut (MAD < 1) recovers the families our CN splits | ⛔ **Every looser cut loses on DEV and HELD-OUT: < 1.25 recovers 4 / breaks 14 (held-out exact 234 → 230); < 2: 10 / 41 (215); < 8: 15 / 73 (199).** `soto_m2_loosen.py`. Fix the copy numbers, not the cut. |
| M2-8 | 2026-09-30 | Soto bipartite | One-to-one (Hungarian) matching of our families to Soto's | ⭐ **All genes: sensitivity .935 / precision .975 (micro), .935 / .982 (macro); tie order does not move precision; same 411, ours bigger 20, ours smaller 36, different 5, no partner left 19.** |
| M2-9 | 2026-09-30 | July artifact | The July "every member" page's "RNA missed" filter shows the 49 missed members | ⛔ **It showed the 313 FOUND rows (`d.rna!=='yes'`); fixed only in the merged page's embedded copy (`build_page.py`); the original artifact EeN5aPefUpJNjjoqgDVt6a still has it.** The two July pages also label 6 members differently (genuine-miss vs seeding gap). |


---

## REGISTER_DRAFTS_o3_maternal.md

*(original: `docs/archive/2026-10/merged/REGISTER_DRAFTS_o3_maternal.md`)*

## Proposed register rows — maternal / paternal reference study (2026-10-08)

Drafts only (not in `NEGATIVE_RESULTS_REGISTER.md`; numbers are assigned on entry). Source: `docs/O3_MATERNAL_REFERENCE_2026-10-08.md` (rewritten after an independent review); prereg `docs/PREREG_o3_maternal_reference_2026-10-08.md` (body + Amendments 1, 2).
Format as the register: `| # | date | study | result | caveat |`. Per prereg Q1/S7 nothing is pooled across loci.

| draft | 2026-10-08 | Do the reads of a copy the reference haplotype lacks hide as unmapped or multi-mapped reads? Real KB3781 fibroblast Iso-Seq against its own maternal (`mat`) and paternal (`pat`) assemblies, 34 families + LRPAP1; per LARGE copy (≥ 20 reads); registered bar of the 09-23 prereg (`bench/o3_maternal/fate.py`) | ⛔ **Refuted at the one copy that tests it.** 19 of the 20 judged large copies (7 with `mat` as reference, 13 with `pat`) have a relative on the reference within the allele cutoff (identity 0.9911-1.0000), so their reads are absorbed at divergence near zero. GWFAM175_B0 (304 reads, relative at 0.9126): 285 absorbed on another locus, 17 on the nearest relative, 2 tied, 0 unmapped, de 0.0662 on `mat` against 0.0007 on `pat`, MAPQ 60 on both. SUPPORTED at GWFAM205_B0 (81 of 82 tied, two equidistant relatives) and GWFAM390_B1 (20 of 20 tied). Control of all other reads 6.4% tied on `mat`, 5.0% on `pat` | ⚠ UNMAPPED is 0 observed (every read mapped on `_pri` first); one animal, one tissue; the verdict function uses fractions only (the bar's "de ≈ divergence" clause is not implemented); smaller beyond-cutoff copies reported not judged: GWFAM26_B2/B3/B6 on `pat` (17/15/8 reads) are 9/5/7 tied |
| draft | 2026-10-08 | Where unmapped reads go when the reference lacks the sequence: the 959 reads unmapped on the combined primary assembly, mapped to `mat` and `pat` | 132 map on `mat`, 3 on `pat`; 125 of the 132 lie in one 71-kb stretch of mat chr12 (95.96-96.03 Mb, MAPQ 60, median length 3,607 bp), hit by no family copy: an expressed mother-only sequence whose reads stay unmapped without it. The unmapped fate exists, at one locus | ⚠ the one read set not pre-filtered by `_pri`; not in any chain net (prereg S6 not followed); no annotation attempted |
| draft | 2026-10-08 | Amendment 1 prediction S1: median `de` lower on the haplotype that has the copy | ⛔ **FAILS at 4 of 8 (`mat`) and 11 of 13 (`pat`) large copies**: where the reference has an identical or near-identical copy the medians are equal, and for `pat` the label itself is wrong in places (GWFAM70_B4: de 0.0013 on the reference's own copy against 0.0090 on the labelled "absent" locus). Second clause (≤ 0.00958 on the present haplotype) holds for 20 of 21; GWFAM382_B0 is 5-6% divergent from BOTH genomes (a copy in neither assembly, not tested here) | ⚠ the synteny-lift truth misstates absence both ways and `bonly` cannot list copies on chromosomes `_pri` took from the reference haplotype |
| draft | 2026-10-08 | LRPAP1 on the mother's genome: which copies are absent? | By the registered rule only the chrY copy (sex). p12 is not absent: `truth_lift` class `T?`, lift 0.914, counterpart at mat chr12 27.43 Mb at 0.972; chr12 has 4 loci on pat and 3 on mat with many-to-one lifts. All 83 p12 reads are TIED with p14 on `pat`; on `mat` 66 absorbed on the p14 ortholog (0.9959), 16 tied. The 9 flagged IsoCon transcripts of the family match pat c01 (chr12 22.55 Mb) at ≥ 0.998 with best mat match 0.987-0.990, so the extra chr12 locus may be c01 | ⚠ corrects the artifact's earlier "mat 7 / pat 8" reading of p12 (prereg Amendment 2); the c01 reading is post hoc and untested |
| draft | 2026-10-08 | Recovery chain from the mother's genome alone (R1-R3, registered bars of Amendments 8-10), IsoCon 0.3.3 vs in-house `o3_candidates` WITH the Amendment 15/15b fix (acceptance not run) | IsoCon: **R1 PASS** (GWFAM175_B0 recovered), R2 PASS (13 copies, none flagged), ⛔ **R3 FAIL 6/24 = 25%** by the registered rule (bar ≤ 20%); post hoc, 4 of the 6 families (GWFAM309, 318, 70, LRPAP1) are sequences at ≥ 0.998 on the paternal genome with best maternal match 0.980-0.990 and no truth locus, only GWFAM169 and GWFAM382 match nothing. 1 of 16 expressed copies recovered, 9 flagged. In-house: ⛔ R1 FAIL (1 flagged, nothing recovered), R2 PASS, R3 PASS 1/24 | ⚠ one beyond-cutoff copy decides R1; the in-house binary is not the registered b29afa55 one; the truth table is incomplete for pat-sourced chromosomes |
| draft | 2026-10-08 | The same chain with the father's genome as reference (held-out direction, run after the `mat` verdicts were written) | R1 NOT TESTABLE (no large copy beyond the cutoff), R2 PASS (26 copies, none flagged), R3 PASS (IsoCon 2/21, in-house 0/21); no candidate falls in a truth locus; flagged IsoCon candidates in GWFAM162, 256 and 70 are mother-genome sequences at ≥ 0.998 with best paternal match below 0.99 | ⚠ reverse truth = `bonly.tsv` rows `hap = mat`; T1-candidate copies (GWFAM175:2, GWFAM491:1) excluded |
| draft | 2026-10-08 | The same copy seen by both runs (Amendment 1 S2) | GWFAM175_B0: 1 IsoCon output on the locus in each run; kept as a new copy and recovered in the `mat` run, matches the reference (0 new) in the `pat` run. Other copies: outputs exist in both runs (e.g. GWFAM208_B0 9 vs 5), none kept as new (linked within delta) | ⚠ p12 line: IsoCon outputs ≥ 0.999 on p12 (pat) 1, p14 (pat) 1, p14 (mat) 1, p12 counterpart in mat 0; the p12 and p14 outputs differ by 4 bp and link by delta; the p14 (mat) output is a different isoform assigned on 2 bp |


---

## REGISTER_DRAFTS_unmapped_rescue.md

*(original: `docs/archive/2026-10/merged/REGISTER_DRAFTS_unmapped_rescue.md`)*

## Proposed register rows — pool-first net for unmapped reads (2026-10-08)

Drafts only; rows from Amendment 9 on are corrected in part by the Errata of the prereg (independent review 2026-10-08), e.g. the chain-aware row now reads dev A > 200 edits 22 -> 13 and unsupported 24 -> 14 with the fixed all-vs-all, and the T4 row is exact 23/24 under the start-tolerant scoring; source `docs/UNMAPPED_RESCUE_2026-10-08.md` (rewritten after an independent review); prereg `docs/PREREG_unmapped_rescue_2026-10-08.md` (body + Amendments 1 to 4). Format `| # | date | study | result | caveat |`.

| draft | 2026-10-08 | Do the unmapped reads of an erased copy form coherent piles without a reference? Seed-round clustering of the unmapped pool (`minimap2 -x map-hifi`, de ≤ 0.00958, block ≥ 50% of the shorter read, components ≥ 3) | ✅ **U1 PASSES on the held-out bed**: 5,083 of 5,313 deleted-copy reads (95.7%) in 47 clusters at purity 1.0000, 0 background reads in a family cluster. Dev bed (162 families, 59,791 reads): 89.2% / 96.4% / 98.8% at δ/2, δ, 2δ, purity 1.0000 at all three | ⚠ synthetic erasure of the same animal's copies; the registered all-vs-all did not finish (quadratic in 7,000-read families), seed rounds replace it; a sketch-only overlap (`dv`, no alignment) reproduced the clusters 4× faster on an 8-family dev subset (not frozen) |
| draft | 2026-10-08 | Attribution of a cluster consensus to a family: mmseqs2 translated bit score, strictly-above rule (prereg body) | ⛔ **FAILS**: held-out 10.5% of joined reads wrong, dev 16.5% (cluster accuracy 0.59-0.66) | ⚠ replaced under Amendment 1 before bed H was run |
| draft | 2026-10-08 | The replacement: dc-megablast cover score (consensus bases covered by the family's best copy), attribute iff ≥ 1.10 × the runner-up, BLAST targets uncapped | ✅ held-out bed H: **4,255 of 5,312 unmapped deleted-copy reads rescued (80.1%), 1 wrong cluster of 27 (3 wrong joins, cluster accuracy 0.963), 20 of 32 copies reached**; A13 rescued 557 (10.5%). Dev bed A: 55,901 of 59,791 (93.5%), 53 wrong joins, cluster accuracy 0.937 | ⚠ the 1.10 margin is a dev-bed choice that sits just above a few wrong clusters (A: 219 wrong joins at 1.09, 53 at 1.10; H: 129 and 3); with the cap at 500 the dev numbers were 87.8% and 7.1% for the strict rule; abstentions are mostly copies without a nucleotide relative (H: 414 reads) or with the true family ranked below the top two, not near-ties |
| draft | 2026-10-08 | Is pooling what rescues the reads? B1 (single-read dc-megablast, same score and margin), pooled, 500 seeded reads per bed | ⛔ **U2's second clause FAILS**: pooled 406 vs B1 410 (H), 469 vs 475 (dev); pooling is ~6× cheaper in wall-clock | ⚠ B0 (map-ont of single reads to genomic copies) is not A13's rule (A13 aligned to the family's net reads too) and rescues 0 of 500 on both beds, so B1 > B0 is vacuous |
| draft | 2026-10-08 | Do the consensus sequences reproduce the erased copy? Registered metric: identity × coverage ≥ 0.999 against the unmasked genome | 14 of 40 (H) and 65 of 212 (dev) clusters; over attributed clusters 12 of 27 and 39 of 125; by identity alone 20/40 and 95/212 | ⚠ "lands on the erased interval" holds for all clusters but follows from how truth is assigned |
| draft | 2026-10-08 | Real bed: the 125 reads unmapped on the primary assembly that map to mat chr12 95.96-96.03 Mb | ✅ 115 of 125 (92%) in one cluster with no other read, consensus 5,694 bp at identity 1.0000 to the maternal region; reads from all four runs (not run-exclusive) | ⚠ clustering part visible in the dev run; family attribution is NOT stable: full catalog best family GWFAM267 407 vs 371 (ratio 1.097, abstains), driven by a short element hitting 831 of 915 copies |
| draft | 2026-10-08 | Can a family-level signal beat the cover score? Family union, oracle own-family search, repeat-robust specific cover with an evidence floor (Amendments 2-3; exploratory, bed H no longer held out) | Union: no change. Oracle: 11 of 13 abstained held-out clusters have homology to their own family (e.g. 475 bits) against ≤ 150 bits for random families. Specific cover + floor 0.02 (selected on A): H 4,325 rescued but 271 wrong (5.90% > 5%): FAILS | ⚠ f = 0.05 (post hoc) gives H 4,321 / 27 wrong and A 55,852 / 14: reported, not claimed |
| draft | 2026-10-08 | Per-family HMM profile (MAFFT of the survivors' transcript consensus, `hmmbuild --dna`, `nhmmer`; Amendment 5; H only, 53 families; exploratory) vs the dc-megablast cover score in the SAME universe | HMM **4,522 rescued, 125 wrong (2.69%), 22/24 copies** vs baseline 4,279 / 59 (1.36%) / 21/24: +243 reads, from two clusters of families with 9 and 5 survivors (200 and 70 reads); ⛔ the GWFAM163 profile attracts five small clusters of other families (hub), the 24-read GWFAM440 cluster is lost; families with >= 3 survivors: 496 rescued / 14.3% wrong vs 226 / 20.7% | ⚠ bed H not held out for it; two clusters carry the gain; needs a new multi-copy erasure bed; a family with one survivor gets a single-sequence profile and no gain |
| draft | 2026-10-08 | Hybrid net (Amendment 6): pool first, then residual reads attributed one by one, bed H | ⛔ registered rule NOT met: 4,421 vs pooled 4,255 vs B1-scaled 4,357, wrong 1.58%, but extra work 32.1% of per-read attribution (bar 25%) | ⚠ bed A not run (rule needs both beds; the H failure fixes the verdict); +166 reads is descriptive |
| draft | 2026-10-08 | Reconstructed loci as extra reference (Amendment 7): reads aligned to consensus built from half the unmapped deleted-copy reads, bed H | ✅ held-out unmapped reads 2,505/2,657 (94.3%) on own-family consensus at de 0.0011; absorbed 2,250/11,973 move (de 0.1237 -> 0.0031, 99.9% own family); survivors 0/41,727 move; background 0/959 onto a family cluster | ⚠ primary emulated (consensus alignment score strictly higher than the genome primary), no concatenated index; survivors are the control for the spliced-vs-unspliced score |
| draft | 2026-10-08 | LRPAP1 worked example (Amendment 8), mother as reference: 381 poorly placed reads -> 4 pure clusters -> consensus placed on mat/pat | ✅ bar met through the 4-read c01 cluster (pat c01 0.9993, mat best 0.9888, flagged); 283-read c01 cluster 0.9987 (not flagged at 0.999; identity 1.0000, 2 terminal bases unaligned), nearest mat node c00 0.9886; 84.8% held-out net reads move to own copy; 2.8% false moves | ⚠ truth = paternal primary placement; moves of already-placed reads are the animal's own sequence replacing the mother's, not rescue; c01 vs p12 is evidence, not a test |
| draft | 2026-10-08 | IsoCon-style correction of the cluster consensus (Amendment 9): reads re-aligned to the consensus, significant plurality alleles replace it; identity x coverage >= 0.999 against the erased copy | ⛔ **NOT met on held-out H**: 14 -> 15 of 40 (bar >= 19), 3 clusters worse (bar <= 2); dev A 65 -> 70 of 212; NM -7.4% both beds; residual edits 1.2%/1.9% allele-like | ⚠ 90% of corrections are exon-scale deletions (mixed-isoform consensus trimmed to its plurality isoform); the metric is dominated by unaligned consensus ends (A: identity >= .999 in 95/212, identity x coverage in 65); LRPAP1 c01 consensus is identity 1.0000, the 0.9987 is 2 terminal bases |
| draft | 2026-10-08 | Supervised synthetic world (Amendment 10): 30 synthetic families x 3 copies, one erased, 80 reads per copy, exact truth, pipeline frozen | ✅ recovery holds: D >= 1% every erased copy has a pure cluster, consensus >= 0.999 to the true transcript 30/30, 0 wrong attributions, erased reads on own consensus 100%; D = 0.5% only 12% of reads leave the cutoff (designed). ⛔ registered S6 FAILS: survivors' reads move onto the consensus 30-50% when D <= 2% | ⚠ cause = spliced-vs-unspliced raw alignment score (about 8 points per junction); divergence rule removes it (0/960) and leaves H unchanged |
| draft | 2026-10-08 | Ends of the consensus: do jitter, 5' truncation or a polyA tail reproduce the real beds' 1-4 bp unaligned ends? (S3/S4) | ⛔ none does: consensus is shorter than the true transcript, aligner clips 0; polyA kills the metric (0/30), absent in real data | ⚠ real beds: 5' untemplated G run, 129/212 dev-A consensus sequences have a 5' clip, 102 pure G, none matches the adjacent genome; excluding it gives 80/212 (A) and 17/40 (H) instead of 65 and 14 (descriptive, Amendment 11); contradicts the 'no gorilla CAP signal' note |
| draft | 2026-10-08 | Divergence move rule adopted (Amendment 12): fresh synthetic world seed 20261009 (32 families, new classes 0.75/1.5/3%), bed H, LRPAP1 | ✅ survivors moved 0 of 640 in every class (score rule up to 320 of 640); erased reads on own consensus 100%; H unchanged (2,505 / 2,657 held-out; 0 of 41,727 survivors); LRPAP1 172 of 191 held-out net reads (162 before), other-copy 3 of 532 (15) | ⚠ detection edge 0% at 0.5%, 55% at 0.75%, 75% at 1%, 100% from 1.5% |
| draft | 2026-10-08 | Untemplated 5' G run: trim a 1-3 bp pure-G prefix not covered by the family's surviving copy (rules T by dc-megablast HSP start, T2 by splice alignment), synthetic worlds seeds 20261010/20261011 with the measured G-run model (95% of reads, mode 2 bp) | T synthetic: 14/24 -> 24/24 but registered bar (1) NOT met (58% > 50%); T2 synthetic: 12/23 -> 21/23, 0 harm, bars met. ⛔ real dev A: T trims 4 of 212, T2 3 of 212, recall 0 of 60 attributed clusters with a genome-based G clip (bar 80%) -> NOT adopted | ⚠ real survivors are 3-15% diverged paralogs: splice:hq gives any record for 44 of 221 consensus sequences; the HSP starts tens of bases in; artifact needs a library-level signature |
| draft | 2026-10-08 | 5' G artifact read from the library (Amendment 15, rule T3): gate on the 5' soft-clip signature of survivor reads, then trim a 1-3 bp leading G run of each consensus | Gate: real library OPEN (G 28,614 vs A/C/T 80), synthetic E3 open, E4 closed. Synthetic (seed 20261012): 12 -> 22 of 24 at 0.999, offset window 14 -> 23, 0 harm, E4 0 trimmed. Real dev A: 65 -> 82 of 212, H 14 -> 18 of 40. ⛔ registered real bars FAIL: recall 92/100 (bar 95%), templated bases lost <= 1 in 132/197 = 67% (bar 90%) -> NOT adopted | ⚠ real 5' ends are GC-rich: 105/197 trimmed consensus sequences have no unaligned G clip; artifact G vs templated G not separable by alignment; alternative = ignore a leading G run <= 3 in the 0.999 flag when the gate is open |
| draft | 2026-10-08 | Per-cluster correction of the 5' G run (Amendment 16, T4): j = templated G's by maximum likelihood from the reads' own leading-G runs (artifact length distribution from aligned survivors), trim the consensus' run beyond j | Synthetic seed 20261013 (0-3 templated leading G): within +-1 base in 24/24 (T3 16/24), templated bases removed 1/24 (T3 7/24), identity x coverage >= .999 24/24 (T3 20, untrimmed 16), control untouched. ⛔ registered bars FAIL: exact 14/24 = 58% (bar 80%), artifact G left 9/24 (bar 10%); real dev A: 65 -> 80/212 (bar 75 met) but under-trim vs genome clip 36% (bar 10%) -> NOT adopted | ⚠ the 9 synthetic shortfalls and 13 of the 36 real are the +1 start ambiguity (jitter 1); 23 of the 36 real are clusters < 5 reads (T4 not applied) |
| draft | 2026-10-08 | Gate-aware identity x coverage (Amendment 17): the registered flag metric ignores a leading G run <= 3 when the library gate is open | ✅ synthetic E3 12 -> 23/24, E5 16 -> 24/24, controls identical, 0 false passes; real dev A 65 -> 80, H 14 -> 17 (not blind); LRPAP1 283-read c01 consensus 0.9987 -> 1.0000, flagged (mat best 0.9899) | ⚠ the c00 allele cluster (0.4% from the mother's) is flagged too: the flag does not separate an allele from an absent copy |
| draft | 2026-10-08 | IsoCon-style partition of a cluster into competing candidates (Amendments 18-20): linked significant variant columns (hypergeometric, Bonferroni) split the reads, recursive | ✅ fresh world seed 20261016: control 0/12 split; ISO baseline 12/24 -> PART 24/24, 12/12 clusters split into pure leaves; 0 spurious/redundant of 60; purity 1.00. Earlier worlds: siblings 19 -> 24, 23 -> 24. Real dev A: 72/212 clusters split, gate-aware 0.999 count 80 -> 106 (oracle best candidate), H 17 -> 26 | ⚠ 22 clusters with > 200 edits are NOT split (0 of 22); real-data over-splitting unknown (up to 14 leaves); scorer bug (best-of-two-orientations core identity) corrected after the first numbers; `map-hifi` soft-clips long deletions |
| draft | 2026-10-08 | The 22 dev-A clusters with > 200 edits against the erased copy (Amendment 21): truth-free flag "consensus unsupported" = median divergence of the cluster's reads to its own consensus > 0.00958 | ✅ H (test): flags 2 of 3 big-edit clusters, 1 of 37 others (2.7%); dev A (design): 16 of 22, 8 of 190 others (4.2%). Mechanism: tiny clusters (3-12 reads) that are single-linkage chains of reads with very different spans; reads fit the genome individually (divergence 0.0013), their POA consensus does not (0.016) | ⚠ 6 dev-A clusters (1 on H) have a consensus that fits its reads but differs from the erased copy by 238-614 edits: cause not identified; 0.2% of reads, 10% of clusters |
| draft | 2026-10-08 | Chain-aware clustering (Amendments 22-24): proper-overlap edge + longest-first star clustering inside components <= 60 reads with every gap counted (NM / block) | ✅ fresh synthetic chain world (seed 20261018): isoform transcripts 24/24 vs frozen 12/24, purity 0.91, control untouched, coverage 99-100%. Real (descriptive; fixed all-vs-all): clusters with > 200 edits 22 -> 13 (A), 3 -> 0 (H); unsupported consensus 24 -> 14, 3 -> 0; erased-copy clusters at .999 38% -> 52% (A), 43% -> 56% (H); cost 1.4 / 2.7 points of read coverage; 71 of the 126 reads of the A clusters leave the clustering. ⛔ earlier steps failed their bars: proper-overlap edge alone (Amendment 22: H > 200 edits 2 vs <= 1), star step with gap-compressed divergence (Amendment 23: CHAIN 12/24) | ⚠ real-bed numbers are not blind (both beds used for earlier amendments); the 7 remaining dev-A big-edit clusters are the 'consensus fits reads, not the erased copy' kind; purity 0.91 = contained fragments identical in both isoforms |
| draft | 2026-10-08 | Specificity controls (Amendment 25): PART block-purity guard; COPY / ALLELE / PRESENT classification of the O3 call | ✅ guard: hotspot controls split 10-12 of 12 -> 0-1 of 12 (two fresh seeds), exon skip and 0.5% siblings 12 of 12, fresh partition world unchanged (ISO 24/24, control 0/12). O3 class (seed 20261402): NULL 12 families -> no cluster; alleles 0.25% / 0.5% -> 0 COPY (registered flag calls the 0.5% one); erased copy COPY 12/12; LRPAP1 c00 allele = ALLELE, c01 = COPY | ⚠ alleles >= 1% apart are called COPY (4/6 at 1%, 6/6 at 2%): RNA cannot separate them from paralogs; LRPAP1 c01 is 0.0005 below the allele cutoff; real-bed splits (29% A, 36% H after the guard) still unvalidated against alleles |
| draft | 2026-10-08 | Three-way O3 call on the maternal-reference study (Amendments 26-27): reads of 34 families + LRPAP1 with unmapped or divergent primary on the reference haplotype, adopted clustering, COPY / ALLELE / PRESENT; then a reciprocal-best-match (orthology) test | Predicted calls hold: mother as reference COPY for GWFAM175_B0 (269 reads, R 0.93) and LRPAP1 c01 (282 and 4 reads); father as reference no COPY among judged-copy labels. ⛔ registered specificity bars FAIL: 5 of 28 (mother) and 7 of 24 (father) 'shared' clusters are COPY (diverged orthologs 1-7% apart). Reciprocal test: 9 of those 12 -> ALLELE (5/5 mother, 4/7 father), GWFAM175_B0 + c01 stay COPY; 3 residual 'shared' COPY incl. a 480-read cluster = mirror paralog of GWFAM175_B0 (real copy-number polymorphism missing from the synteny truth list). Bed H: 31/55 erased-copy clusters COPY, 0 ALLELE | ⚠ the reciprocal rule is NOT adopted by its bars (S2 4/7, S3 3 vs <= 1); orthology, not distance, is what separates allele from copy; the registered synteny truth list is incomplete in both directions |
| draft | 2026-10-08 | Discovery of reference-absent sequence from the 959 reads unmapped on the primary assembly (Amendments 28-31): adopted clustering, consensus, classes against the primary and the animal's two haplotype assemblies (confirmation only), terminal-exon rescue, support requirement | ✅ 3 clusters (120 reads): a 113-read, 5,694 bp consensus absent from the primary and 1.0000 identical to the mother's chr12 95.96-96.03 Mb (84% of the 134 confirmable reads); 2 supported clusters (4 and 3 reads, 2.4 kb) with NO alignment on any of the three assemblies. Clean controls: seed 3 nothing called; ⛔ fresh seed 4: 1 CONFIRMED cluster (9 reads, 1.3% of clustered reads) = the father's haplotype of chr11 locus, primary 1.1% away | ⚠ registered specificity bar not met on the fresh control (allele at the 1% boundary); the first definitions failed D1 (137-base unaligned terminal exons) and D3 (blended 462-read consensus; control not clean) and were repaired in Amendments 29-31 on the same data; the 2 NOVEL clusters need DNA verification (not done) |
| draft | 2026-10-08 | DNA verification of the discovery candidates (Amendment 32): exact 21-mer lookup of each reconstructed sequence in all of KB3781's whole-genome reads (277.6 M reads, 41.9 Gb) | ✅ 113-read sequence: 97.1% of k-mers >= 2, 99.8% of bases covered, median count 6.0 (one haplotype, lambda 4.98); PRESENT control 99.9% covered; ✅ the two NOVEL clusters (4 and 3 reads) are DNA-ABSENT: 0.3% / 0.2% of k-mers present = shuffled negatives (0.3% / 0.1%), no alignment on any gorilla assembly or on human CHM13 -> RNA-only products, rejected | ⚠ V2 as written (shuffled negatives covered <= 2%) not met: 4.8% for one (chance k-mers each cover 21 bases); k-mer fraction separation 97% vs 0.3%; the nature of N1/N2 (culture contamination, chimeric or trans-spliced products) not determined |

| draft | 2026-10-08 | Discovery on four more libraries (Amendments 33-34): unmapped reads of a second gorilla (testis), human A119b, chimpanzee, orangutan; same pipeline, classes vs the animal's own reference and the other primates' assemblies; species check of the 106 clusters | ✅ E1/P1: the second gorilla's largest cluster (30 reads, 5,523 bp, absent from the primary) is T1 of KB3781 at 99.82% end to end; 47/51 ELSEWHERE reads gorilla >= 0.99. ⛔ E2 specificity control met in 1 of 4 (gorilla 96.7%; human 86.7%, chimp 88.7%, orangutan 94.1%); rule defect: ELSEWHERE should need E > R. ⛔ P2: human clusters (619 reads) are no CHM13 hit and no single panel species >= 0.99 (82% of reads in 0.9-0.99): not human discoveries; chimp/orangutan: 0 same-species clusters. Unresolved: species of the human-library clusters (needs nt BLAST) | `docs/PREREG_unmapped_rescue_2026-10-08.md` Amendments 33, 33b, 34; `bench/unmapped_rescue/{discover,species_check}.py` |


---

## PENDING_2026-09-23.md

*(original: `docs/archive/2026-10/merged/PENDING_2026-09-23.md`)*

## Pending work — registered 2026-09-23 (user, 15:33), to be pre-registered and run later

Both items are the user's, recorded so they are not lost. Neither has been started. `docs/archive/2026-09/PENDING_2026-09-17.md`
items remain open unless closed elsewhere.

### 1. Quantification parity with StringTie / FLAIR / isoseq collapse on the transcripts all of us find

**Goal (user):** "ensure the quantification is comparable to StringTie's or the other tools for the matching
transcripts that we also find; proving that we are close would convince my advisor — TPMs and RPKMs if needed."

**What exists.** `copy_assign --gtf-tpm` (§6r5) writes count-based `TPM` and `cov` per transcript from the
`reads "N"` attribute (a long read is one molecule, so counts are NOT length-normalised: measured on chr20,
count-based ρ 0.879 vs StringTie's TPM on 507 shared chains, length-normalised only 0.714; vs FLAIR's isoform
counts ρ 0.755 on 388 shared chains). The lab's genome-wide outputs carry their own abundances:
StringTie `cov`/`FPKM`/`TPM` attributes, FLAIR `*.isoform.counts.txt`, isoseq `*.collapsed.abundance.txt`
(`count_fl`). The strand-aware chain key of `consensus.py` (`/mnt/linuxdisk/tmp/gw22/ext/`) already joins
transcripts across tools.

**Design to pre-register.** On the SHARED chains (ours ∩ tool, per species, human genome-wide and gorilla
genome-wide — the new-default assembler output): Spearman and Pearson (log1p) of our read count vs the
tool's count (isoseq FL, FLAIR count, StringTie `cov`), and of our TPM vs StringTie's TPM; RPKM as
reads / (spliced kb × M reads) for the advisor's convention, with the caveat stated up front that a
long-read count is already a molecule count and RPKM double-normalises. Report per tool, per species, with
the number of shared chains, and a per-family view (NPIP, TBC1D3) because that is where the thesis lives.
⚠ Known from r1058: our counts are DE-DUPLICATED (coordinate-identical reads collapse), the tools' are not —
expect a systematic ratio below 1 at deep loci; report it with `--keep-coordinate-duplicates` too.
Bar suggestion: ρ ≥ 0.85 against StringTie and isoseq on ≥ 1,000 shared chains per species.

### 2. O3 — does reference bias hide copies as unmapped or multi-mapping reads? (advisor's claim)

**Claim (advisor):** copy numbers match only when the sample is the individual the reference was built from;
otherwise the extra copies show up as unmapped reads or are disguised as multi-mapping reads.
**User's experience:** SNPs and small variants, yes; entire missing copies, no. **Goal:** machinery to prove or
disprove it.

**What exists.** O3's flag pass (`copy_assign --flag-missing-copies`, `o3_flag_pass.rs`: per-copy divergence
rate, Bonferroni threshold, `o3_candidate_loci.tsv`) is the detector for "reads at a copy that fit no reference
copy". Cross-sample material on disk: gorilla `GGO_mm.bam` (testis, animal OR6737) on the GCF_029281585.2
reference (a DIFFERENT animal, Kamilah), `fibroblasts/GCA_029281585.2_flnc_mm.bam` (⚠ verified 09-23 from its `@PG`: the reference
individual's own fibroblast FLNC reads on the SAME `GGO_GCF_029281585.2` primary — the accession in the
file name is the read set's, not another assembly's — so it is the same-individual control on the same
placement; the parental haplotypes are `o3_hapcnv/{pat,mat}.chr.fa`), `PPY_mm.bam` / `PTR_mm.bam` (orangutan / chimp reads on their own references),
`*_ds.bam` down-sampled copies, and the simulators (`bench/sim.py chromosome`, `bench/sim.py missing-copy`, and the read model
`sim.simulate_reads`; before wave 7, `bench/ideal_chromosome_sim.py`, `missing_copy_sim.py` and `sim_reads.py`).

**Design to pre-register — three arms, cheapest first.**
- **A. Simulation with truth.** Take one contig, make a "sample genome" by (i) DUPLICATING k annotated
  multi-exon genes at a chosen divergence (0.5-3% from their template, plus a few PSV-like fixed
  differences) and (ii) DELETING k others; simulate reads from the sample genome (`sim.simulate_reads` in `bench/sim.py`, ± 30 bp
  ends, 1-30 reads per transcript), map to the UNMODIFIED reference with the shipped minimap2 settings.
  Measure for the reads of each extra copy: unmapped %, secondary/AS-tied %, absorbed into the template
  copy % (and their `de`), and whether O3's flag pass calls the template locus. This decides the MECHANISM
  question with ground truth: at what divergence does an extra copy vanish into "unmapped", into
  "multi-mapping", or into a paralogue with inflated divergence.
- **B. Same-individual control vs different-individual, real data.** Gorilla reads on the reference's own
  second haplotype assembly (GCA) vs on the primary (GCF), and the OR6737 testis library vs the reference
  animal: compare unmapped fraction, MAPQ-0 / AS-tied fraction, O3 candidate count and per-family copy
  counts (NPIP-like LOCs, TBC1D3-like) between the two placements. If the advisor is right, the different-
  individual placement shows more unmapped/tied reads AND O3 candidates concentrated at multi-copy families;
  if the user is right, the excess is SNP-level `de`, not copy-level.
- **C. Cross-species stress test** (the extreme of the claim): human A119b reads on the gorilla reference and
  vice versa, same measurements — where copies genuinely differ between species, O3 must fire.

**Report** per arm: unmapped %, tied %, absorbed %, O3 candidates (true / false against the simulation's
truth), copies per family. Bar: the claim is SUPPORTED if, at ≤ 2% divergence, ≥ 20% of an extra copy's
reads are unmapped or tied rather than absorbed; REFUTED if ≥ 80% are absorbed into the template with
elevated `de` (the O3 signature). Hold-out: design on the simulation and the gorilla haplotype pair, judge on
the cross-species arm.

**Status 2026-09-23 evening:** item 1 measured (`docs/archive/2026-09/PREREG_quant_parity_2026-09-23.md`, r1088: ⚠ close; use `--keep-coordinate-duplicates` for the comparable count). Item 2 arm A measured (`docs/archive/2026-09/PREREG_o3_reference_bias_2026-09-23.md`, r1089: claim refuted with truth to 5% divergence); arms B (gorilla haplotype pair, needs the fibroblast reads mapped to GCF too) and C (cross-species) remain open.

### 3. O3 — the DNA step for the 52 RNA-only candidates (registered 2026-09-23, 19:28, user)

**Where O3 stands.** `o3_rna_flag` (§6ze, r1090/r1091) leaves 52 candidates on KB3781 fibroblast that no RNA test
kills: consistent PSVs at 1.7% over most exons, absent from the primary, both parental haplotypes, human, chimp
and orangutan; immune-cell genes in a fibroblast line. Whether the sequence is in the animal's genome is the one
question left, and it is a DNA question. Files: `/mnt/linuxdisk/tmp/gw22/o3/fib.o3_rna.tsv`, `fib.consensus.fa`.

**Design to pre-register.** Download the KB3781 (mGorGor1, BioSample SAMN04003007) HiFi DNA reads from SRA — the
assembly's own reads, so any copy the assembler collapsed or dropped is still in them. For each candidate
consensus: the PSV k-mers (k = 31 spanning each PSV site, both orientations) versus the matching host k-mers
(same positions, reference allele). Count both sets in the HiFi reads (jellyfish/KMC or a plain hash over the
FASTQ). Truth classes, fixed before counting: **present as a copy** if the PSV k-mers sit at ≥ 0.5× the host
k-mer depth (the depth a second copy of a diploid locus gives is ~1×; 0.5× allows a hemizygous copy) at ≥ 80% of
the sites; **absent from the genome** if ≤ 0.05× at ≥ 80% of sites (then the RNA has a non-genomic source:
contamination from an unsampled species/individual, editing, cell-line drift); ambiguous otherwise. Positive
control: KIAA0513 (the confirmed maternal copy, PSV k-mers must be at ~1×); negative control: the `scattered`
piles (no PSVs → nothing to count) and 50 random host k-mer sets at 2×. Report: candidates present / absent /
ambiguous, and the same for the 3 chimp-like and the orangutan-borderline ones. Bar: if ≥ 50% of the 52 are
"present", O3 has real reference-absent expressed copies in the matched animal and the RNA flag's precision on
real data is that fraction; if ≥ 80% are "absent", the flag finds non-genomic RNA and the thesis says so.
Cost: the SRA download (tens of GB, to `/mnt/linuxdisk`) is the only heavy step; counting is minutes.

**Also register (same item): HG002 as the human matched-individual test.** HG002 has public IsoSeq and its own
diploid T2T assembly (HPRC/T2T), i.e. the same "reads from the assembly's individual" design as KB3781 but in
human, where population copy-number variation (NPIP, TBC1D3, NBPF) is documented. Run `o3_rna_flag` with the
HG002 assembly as `--confirm` and CHM13 as the primary: candidates confirmed by HG002's own assembly are
individual-specific copies absent from the reference — the population version of the advisor's claim, with
truth. Needs: HG002 IsoSeq BAM on CHM13, the HG002 assembly and its splice index (~20 GB build).

**O1 side note (not O3):** the one improvement O1 still admits is a better truth — Compara gene trees instead of a
single annotation (r1031: the protein edges are right, the truth is missing). A download, not an experiment.

### Added 2026-09-24 — decisions for the advisor (§6zm, §6zk)

4. **Copy assignment: strict tie or expression prior?** `copy_assign --union-certificate` (shipped opt-in, row 1103)
   scores one certificate per tied read over every placement it touches, catalog copies across families plus
   outside loci rebuilt from the genome. It removes every foreign claim (844 → 0) and every wrong row on the chr16
   truth simulation — and ties every MAPQ-0 read, because each one's tied partner is an NM-identical genomic twin
   (e.g. an inverted exact duplicate the RNA catalog does not hold). The baseline's 157 "exact" own-family
   assignments were made blind to that twin. Making them assignable again means privileging the expressed catalog
   copy over an unexpressed identical locus, i.e. an expression prior. Thesis rule so far: no priors, assign-or-abstain.
   The user's/advisor's call; the flag is ready either way (default OFF today).
5. **Seeding default flipped** (`tools/rustle_pipeline.sh assemble` now seeds with tied secondaries, `--no-seed-secondaries`
   opts out; rows 1060/1100/1101): +30 loci on gorilla NC_073244.2 (97% annotated), ≥ 90%-identity referee pairs
   joined 18 → 36 of 108 at pair precision 1.000, streaming assembler validated identical to the buffered path
   (4 s / 0.3 GB vs 42 s / 7.8 GB on the testis slice). Reversible by one flag if the advisor prefers primaries only.


---

## ADVISOR_HIERARCHY_2026-10-03.md

*(original: `docs/archive/2026-10/merged/ADVISOR_HIERARCHY_2026-10-03.md`)*

## "Where is my hierarchy in the MCL?" — what the advisor most likely means, and what exists (2026-10-03)

Written for the next meeting. Sources: `docs/seeded_family_definition.md` §0★/§0★★★, ledger §5k, §6l, §6n7, §6kl,
memory `project_family_hierarchy`, register rows 9, 311, 457, 654.

### The likely meaning
MCL returns ONE flat partition at one inflation value (shipped: 2.8). Multi-copy families are not flat: superfamily ⊃
family ⊃ subfamily (NPIPA / NPIPB, TBC1D3 groups) ⊃ copies ⊃ alleles. "Where is the hierarchy" = where in the output is
that nesting, since a single cut asserts one resolution and hides the others. It is the threshold objection restated for
clusters: a cut is a parameter; a tree with levels is a structure.

The "true positives / false positives" reading is a consequence of the same point: with a hierarchy, TP and FP are
level-relative. NPIPA+NPIPB merged is a false merge at the subfamily level and a true family one level up. Our data
already behaves this way: Soto's families sit between our L2 and L3 (§6n7); NPIP is one family at L1/L2 and splits at L3;
precision against any single-resolution truth caps near 0.8 (§6kl). A tree lets each truth be scored at its own level.

### What exists
| item | status | where |
|---|---|---|
| Hierarchy inside MCL | none; inflation is a resolution knob, clusterings at different inflations are not nested; ours stable I = 2.0–4.0, shatter at ≥ 6 | §6ec, §6ew, register 670 |
| Nested edge-test lattice L0–L4 | proposed design, nesting is a theorem (T1), not the default; L3 (0.985) / L4 (0.995) shipped as thresholds | `seeded_family_definition.md` §0★★★ |
| What breaks nesting | the leader rule (RNA ⊆ DNA), Louvain (T1/T1′) | memory `project_leader_rule_breaks_nesting`, §6n5 |
| Average-linkage tree on gated distance | recovers NPIPA/NPIPB exactly (label-pure k = 2..8) where NO threshold can (A/B boundary = 0.3 % identity window) | memory `project_family_hierarchy`, register 457 |
| What a tree does NOT add | merge heights = 1 − identity relabelled (partial r −0.05); lift over one cut = +0.035 ARI = two gene attachments; best_k oracle-selected | §6l |
| Hierarchy as catalog repair | NO-GO: no coarse level at which NPIP reunites without the blob; broad/recent two-level hierarchy reaches 19 % of families (70 % are 2-copy) | §5k, register 9, 311 |

### The answer to give
"The MCL cut is one level of a nested structure, not the structure. The principled object is the lattice: connected
components of one evidence graph under nested edge tests, so nesting is provable; L3/L4 are its shipped levels, and an
average-linkage dendrogram on the gated distance draws it. The tree is a presentation and a scoring frame (each truth at
its level), not extra information: merge heights carry nothing identity does not."

### Question to ask him first
"Do you mean a nested sequence of clusterings — family ⊃ subfamily ⊃ copy — with provable nesting, so a truth is scored
at its own level? Or a confidence ranking of clusters?" Both are served by the same levels, since each level is a
stricter edge test; confirm which before building anything new.

### Do not re-propose
A subfamily level as a repair for fragmentation (§5k), a broad/recent two-tier hierarchy (rows 9, 311), MCL inflation
sweeps as a hierarchy (not nested), superfamily components over the cluster graph (row 654).

### Terminology (2026-10-07)
When he says Soto's or Yoo's sets are "multi-copy gene families", agree: they are, specifically the young SD-embedded
subset. SD, duplicon and expansion are the layers below the family, not alternatives to it; the wording and the
sources' own usage are in `docs/TERMINOLOGY_FAMILY_SD_DUPLICON_EXPANSION_2026-10-07.md`.


---

## HANDOFF_2026-10-04.md

*(original: `docs/archive/2026-10/merged/HANDOFF_2026-10-04.md`)*

## Handoff — 2026-10-04 (session of 10-03/10-04), for the next agent

### Where things are

- **Repo:** canonical `/mnt/c/Users/jfris/Desktop/Rustle` (never edit from a side session). Work was done in the worktree
  `/mnt/linuxdisk/home/juanfraitu/rustle_m2_soto`, branch `machine2/soto-evidence`; **`github/main` = the branch = `fe44124b`** (everything
  below is pushed). Build only with `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target_m2 cargo … --release`, output to a file.
  Heavy commands foreground via `bash tools/rlock.sh heavy <cmd>` (< 10 min each), light via `rlock light`; no background jobs, no `pkill -f`.
  Commits: author from the canonical repo's git config (name has a space: quote it), trailers `Co-Authored-By: Claude Fable 5.1
  <noreply@anthropic.com>` and `Claude-Session: https://claude.ai/code/session_01DAyQQ6R8drUxY5GsM5kNkb`; push `machine2/soto-evidence HEAD:main`.
- **This week's priority (user, 10-04 00:10): prove O1 and O2 work right.** O3 engineering is parked.
- **Standing rules:** pre-register decision rules before looking (`docs/PREREG_*.md`, commit before the run); hold a substrate back; never tune
  a registered rule to pass; register every result in `docs/NEGATIVE_RESULTS_REGISTER.md` (last row **1246**); memory files under
  `~/.claude/projects/-mnt-c-Users-jfris-Desktop/memory/` (index `MEMORY.md`, kept < 17.1 KB; topic files hold the detail).
- **Pending list:** `docs/PENDING_2026-10-04.md` (items 1-7), updated through the day.

### What was done (chronological, with the files)

1. **O3 `o3_candidates` stage — acceptance A13 passed, control A14 failed, default flip reverted** (`docs/archive/2026-10/O3_CANDIDATES_ACCEPTANCE_A13_2026-10-03.md`,
   `docs/archive/2026-10/O3_CANDIDATES_CONTROL_A14_2026-10-03.md`, rows 1221-1231; ledger `docs/superpowers/plans/2026-10-03-o3-candidates-a13-ledger.md`).
   The stage stays OPT-IN (`tools/rustle_pipeline.sh --candidates`). Cause of the control failure reproduced byte for byte
   (`docs/archive/2026-10/O3_CANDIDATES_CONSENSUS_DEFECT_2026-10-03.md`): the long-insertion vote without a majority test duplicates the consensus' own exons.
   **Fix + re-runs pre-registered as Amendment 15/15b** of `docs/archive/2026-10/PREREG_rna_allele_haplotype_count_2026-10-01.md` — NOT run (user's cost
   decision; ~6 agent dispatches, ~2 h compute). Whole-BAM cost (R23): one 50-family batch > 10 min.
2. **Row-1220 O2 prereg** `docs/archive/2026-10/PREREG_o2_tie_outside_2026-10-03.md` (tie-outside registration at the molecule's best score) — committed, not run.
3. **Strict "found copy" standard** `docs/archive/2026-10/PREREG_spliced_copy_support_2026-10-04.md` with Amendments A-E, results in
   `docs/archive/2026-10/SPLICED_COPY_SUPPORT_2026-10-04.md` (rows 1232-1233, 1237-1244), scorer `bench/copy_support.py` (all rules as columns:
   `strict_*` read-defined junctions, `ann_*` annotated intron coordinates, `chain_*` annotated chains (ISM/FSM), `xc_*` expressed chains
   (>= 3 unique reads), `td_*`/`te_*` TSS-anchored (annotated / modal start), `tc_*` capped starts). **The counting rule the user accepted is
   Amendment E:** a transcription start = a 20-bp bin where >= 3 cap-signal reads (untemplated 5' G) begin; a read supports a copy only if it
   starts within 150 bp of one and carries the first three introns of a chain expressed there; a copy is FOUND when its node's representative
   does the same. NPIP page `https://claude.ai/artifact/F3gJty4egn598SCZ9RBiM1` (version 8) restated: own node 23/21/24 -> FOUND 10/6/2.
   Key findings: the reads are spliced at 24/25 NPIP copies; the most-read REPRESENTATIVE is a 3' fragment; CAT/RefSeq NPIP models' first
   exons are 12-53 kb so the annotated TSS is useless there (`docs/archive/2026-10/NPIP_CHAIN_COMPARISON_2026-10-04.md`, row 1239: only 27% of multi-junction
   reads are an annotated chain); the cap signal marks the real starts at 17/24 copies (rows 1243-1244). Gorilla OR6737 has NO cap signal.
4. **Labeller page for the user's own marks:** `https://claude.ai/artifact/HhEYh41nPuBmXdaAQcnFdL` (db capability; exclusion marking: unmarked
   = supporting; one doc per locus in collection `labels`; read back with `ArtifactData list labels`). Empty so far.
5. **Representative rule test** (`docs/archive/2026-10/PREREG_locus_representative_rule_2026-10-04.md`, `docs/archive/2026-10/LOCUS_REPRESENTATIVE_RULE_2026-10-04.md`, rows
   1234-1236 + Amendment sections; ledger `docs/superpowers/plans/2026-10-04-representative-rule-ledger.md`): flag
   `mcl_families --representative most-reads|most-junctions` (+ driver `RUSTLE_REPRESENTATIVE`), byte-identical when unset. Held-out verdict:
   R_J (most junctions) stays OPT-IN — family clause fails on chr6 (one readthrough-chain representative; both arms ran with
   `RUSTLE_BRIDGE_REGROUP=off` because the Figure-7 GTFs predate f1v2), while found genes rise under R_J on every contig under the read-defined,
   coordinate and TSS-anchored rules and fall under the annotated-chain rules. Dev chr16: full-length representatives split NPIP into
   subfamilies (the advisor's hierarchy, as a side effect). Runner `bench/rep_rule/`.
6. **LRPAP1 dossier for the advisor** (`docs/archive/2026-10/LRPAP1_FAMILY_2026-10-04.md`, rows 1245-1246, page `https://claude.ai/artifact/UvFQf5QXfqZEzvkzcNt23a`
   version 3; work dir `/mnt/linuxdisk/tmp/lrpap1/`): **11 LRPAP1-like loci in gorilla** (8 full-length on chr3/12/12/12/14/16/22/Y at 0.957-0.990
   identity + 3 "solitary" copies of the 5' three-quarters on chr12/14/16), all RefSeq-annotated ("alpha-2-macroglobulin receptor-associated
   protein-like"); the advisor's reference counts 10 (ours minus the chr16 fragment). O1 de novo: fibroblast MCL2 = the 6 expressed full
   copies + MCL11 = the 3 fragments; testis MCL9 = all 8 + MCL44 = 3; full-length representatives; O2: 0 contested reads (aligner separates
   1-4% copies); O3: maternal 7 / paternal 8 (+Y), nothing missing; strict FOUND 11 of 11 across the two libraries.
7. Housekeeping: `docs/O3_STATUS.md` 10-04 header; `docs/MEMORY_DIGEST.md` "Compacted 2026-10-03"; memory topic files
   `project_rna_haplotype_count_parked`, `project_locus_read_pool`, `project_locus_representatives_study`, `project_lrpap1_family`,
   `project_week_2026-10-04_prove_o1_o2`, `project_closed_threads_2026-09-22_24`, `project_assembler_index_lines`.

### What is left (in the order the user cares about)

A. **O1/O2 proof (this week).** The evidence that exists is listed in `docs/PENDING_2026-10-04.md` ("This week"); the gaps opened this session:
   - **Representative rule, next prereg:** representative = the most-read chain that starts at a capped start (human) / modal start (gorilla),
     judged with the family metrics (Figure 7 scorers via `bench/rep_rule/`) and the Amendment-E FOUND counts; the f1v2 re-test of the chr6
     failure can use the frozen `f1_bridge.py`/`f1v2.py` under `/mnt/linuxdisk/tmp/rustle_figures/f1v2_frozen` (no re-assembly).
   - **Figure 6d / 8 denominators** ("read-supported" = >= 2 overlapping reads) to be re-scored under Amendment E (human) / C (gorilla);
     Figure 7 recovery with the strict FOUND rule (`bench/copy_support.py` generalised: see `bench/rep_rule/score.py h3-inputs`).
   - **Hierarchy for the advisor:** NPIP subfamilies under full-length representatives and the LRPAP1 8 + 3 split are two concrete cases
     of MCL partition vs the reference's cover (`docs/ADVISOR_HIERARCHY_2026-10-03.md`).
   - Row-1220 O2 prereg run (small: implement `--tie-outside-at-best`, A/B on `/mnt/linuxdisk/tmp/rustle_figures/o2sim/{human,gorilla}`).
B. **User's marks:** when `labels` in the labeller page is non-empty, derive the rule from what was left in and compare with Amendment E.
C. **O3 (parked):** Amendment 15/15b run; HG002 held-out when the data arrives; O2/O3 chapter drafts (in this session's scratchpad only —
   `CHAPTER_O2_DRAFT_2026-10-03.md`, `CHAPTER_O3_DRAFT_2026-10-03.md`; the numbers are in the docs cited above, so they can be re-derived).
D. **LRPAP1 follow-ups** (if the advisor asks): date the duplications; the transcribed pseudogene copies (chrY, chr22); haplotype
   resolution of the chr12 fragment; a guided-mode family run on the annotation for the 11 (the annotation already names all 11).

### How to resume quickly

- Read `MEMORY.md` (index) and the topic files named above; then `docs/PENDING_2026-10-04.md`.
- NPIP numbers: `docs/archive/2026-10/SPLICED_COPY_SUPPORT_2026-10-04.md` (sections in amendment order; the per-copy TSV `docs/SPLICED_COPY_SUPPORT_hsa_npip.tsv`
  carries every rule's columns); work dir `/mnt/linuxdisk/tmp/readpool_npip/` (`support_hsa_{A,B,C,D,E}.*`, `chains_hsa.*`).
- Representative rule: `/mnt/linuxdisk/tmp/rep_rule/<species>_<contig>/` (arms `R_M.*`, `R_J.*`, `h3{,A,B,C,D,E}.support.*`).
- LRPAP1: `/mnt/linuxdisk/tmp/lrpap1/` (families, O2, haplotype PAFs, support tables).
- The canonical checkout's `figures/inputs.local.tsv` has the machine paths; the worktree has a copy pointing `bin` at the m2 build.


---

## PENDING_2026-10-04.md

*(original: `docs/archive/2026-10/merged/PENDING_2026-10-04.md`)*

## Pending — 2026-10-04 (user decision: this week is for proving O1 and O2 work right; O3's stage work is parked here)

### Parked, with everything needed to resume

1. **Amendment 15 / 15b — the `o3_candidates` consensus fix and its re-runs** — **DONE 2026-10-05** (commits a13b817f,
   209a5fd7, 7e9421bd; result `docs/O3_CANDIDATES_AMENDMENT15_2026-10-05.md`, register 1247): the fix is in (cs
   normalisation, one-copy-per-long-insertion, the majority test); A13-1/2/3 pass, C2'/H3 hold, but C1' (10/53) and H1
   (6/30) fail, so the stage STAYS OPT-IN (`--candidates`, R14). R23 repeated: whole-BAM cost unchanged. Remaining open
   O3 items: HG002 (item 3), the few-read false-flag phenomenon (a future amendment's subject, named in register 1247).
2. **Row-1220 prereg — O2 tie-outside registration at the molecule's best score** (`docs/archive/2026-10/PREREG_o2_tie_outside_2026-10-03.md`): implement
   `--tie-outside-at-best` in `copy_assign` (the §6gz block, `src/bin/copy_assign.rs` ~:5895-5945; carry the outside record's AS in the
   registry), A/B on the fig4/fig5 human O2 inputs (`/mnt/linuxdisk/tmp/rustle_figures/o2sim/human`), gorilla held out. Small: one
   implementer, one review, ~20 min of runs. It belongs to the O2 proof, so it may be pulled forward this week.
3. **HG002 held-out of the O3 chain** (second individual; the user is obtaining the data).
4. **O2 and O3 chapter drafts** — written 2026-10-03 with every number checked against its result doc; held outside the repo until the
   O2 proof settles the operating points to quote (`CHAPTER_O2_DRAFT_2026-10-03.md`, `CHAPTER_O3_DRAFT_2026-10-03.md` in the session
   scratchpad; the O3 draft already carries the A14 outcome and the defect).

5. **Representative rule and recovery denominators (from `docs/archive/2026-10/SPLICED_COPY_SUPPORT_2026-10-04.md`, rows 1232-1233) — belongs to the O1
   proof, so it may run this week.** (a) **ANSWERED 2026-10-04** — prereg 3bf6aca1, result `docs/archive/2026-10/LOCUS_REPRESENTATIVE_RULE_2026-10-04.md`
   (rows 1234-1236): most-junctions is NOT adopted (held-out chr6 clause (a) failed, Compara F 0.444 -> 0.372; run on pre-f1v2 GTFs with
   `RUSTLE_BRIDGE_REGROUP=off`, forced) and stays opt-in; it does carry the copy's structure (H3 found genes up on all 7 contigs, NPIP
   strict 10 -> 20). Two follow-ups it opened, each a NEW pre-registration, held-out first, never tuned on chr16:
   (i) the f1v2 re-test — the same comparison on f1v2-regrouped GTFs. The frozen `f1_bridge.py`
   (`/mnt/linuxdisk/tmp/rustle_figures/f1_frozen/`) + `f1v2.py` (`/mnt/linuxdisk/tmp/rustle_figures/f1v2_frozen/`) regroup an existing
   base GTF from the BAM without a re-assembly. (ii) A support-aware representative rule: 16 of R_J's 46 H3 losses stay on the gene but
   carry fewer >= 3-read-supported junctions, because the junction count is blind to per-junction support (assembler floor 2 reads).
   (The original item: pre-register a representative rule that keeps the structure, the transcript with the most read-supported
   junctions or the exon union of its >= 2-read transcripts, and re-score the family metrics and the NPIP copy table on the held-out
   chromosomes.) (b) Re-score Figure 6d / Figure 8's "read-supported" loci (>= 2 overlapping reads) as spliced-expressed (>= 2 reads with >= 2
   read-supported junctions; intronless loci by >= 50% coverage) and Figure 7's recovery with the strict "found" rule
   (`bench/copy_support.py` generalised to a copy table + BAM); report old and strict side by side.

### This week: prove O1 and O2 are working right (the user's instruction, 2026-10-04 00:10)

Not planned here beyond the pointers that exist: O1 — `bench/famsim` (the simulator that proves its own condition, 2026-10-02), the
held-out family metrics (sensitivity / precision / bipartite, `docs/THESIS_OBJECTIVES.md` O1 row, `figures/` 6-8), the de novo vs guided
gap; O2 — the simulated-truth accuracy (`bench/o2_read_truth_sim.py` / `_score.py`: 157/157 own-family, 0 wrong), the excision abstention
sweep (247/262), the held-back family MCL58, the hard-locus bakeoff, figures 4-5, and item 2 above.

### Parked 2026-10-04 afternoon (user: "let's park it here for now")

6. **NPIP read support, Amendments A-E done** (`docs/archive/2026-10/SPLICED_COPY_SUPPORT_2026-10-04.md`, rows 1232-1244; page F3gJty4e v8): the counting rule
   is Amendment E (capped starts); the user's own marks go into the labeller page (artifact HhEYh41nPuBmXdaAQcnFdL, exclusion marking, `labels/`
   collection) — read them back with `ArtifactData list labels` and derive the rule from what was left in. Next prereg when resumed: the
   representative = the most-read chain starting at a capped start (human) / modal start (gorilla), judged with the family metrics.
7. **LRPAP1 family dossier** delivered (`docs/archive/2026-10/LRPAP1_FAMILY_2026-10-04.md`, row 1245, page): open items listed in the doc (dating the
   duplications; the transcribed pseudogene copies on chrY/chr22).

6. **Inverted copies vs the strand-aware spliced aligner — parked measurement** (`docs/PREREG_inversion_aligner_blindness_2026-10-05.md`,
   written before any run, NOT RUN by the user's cost decision 2026-10-05). Mechanism: `splice:hq`'s strand constraint (not seeding)
   refuses inverted-copy transcripts, so inversions are unmapped-pool/O3 material by construction — for the real BAMs too. The parked
   famsim experiment (three_identical control + inv_whole / inv_exon / inv_intron; claims C1-C4) would measure the aligner's loss and
   whether `o3_candidates`' map-ont pass B recovers each inversion type. ~1 h of compute when wanted.


---

## TERMINOLOGY_FAMILY_SD_DUPLICON_EXPANSION_2026-10-07.md

*(original: `docs/archive/2026-10/merged/TERMINOLOGY_FAMILY_SD_DUPLICON_EXPANSION_2026-10-07.md`)*

## Terminology: multi-copy gene family vs. segmental duplication, duplicon and expansion (2026-10-07)

**Why this file exists.** Recurring disagreement with the advisor: he calls Soto 2025's catalogue and the Yoo 2025
expansions "multi-copy gene families"; the thesis author has been calling them "families within segmental
duplications" / "gene expansions". Resolution: **both are right, because the terms name different levels of one
object.** The advisor names the object; the thesis names the layers below it. The phrasing "they are expansions,
not families" is a category error (an expansion is an event, a family is its result), and that phrasing, not the
content, is what made the exchanges go badly. Adopt the wording below in the thesis, in talks and with the advisor.

### Definitions (one object, four levels)

| term | what it is | kind | evidence it is read from |
|---|---|---|---|
| **multi-copy gene family** | a set of loci that are paralogs (descended from one ancestral gene by duplication) | the OBJECT; O1 defines it | connected block of the copy graph: nodes = genomic intervals carrying reads, edges = contiguous high-coverage homology on assembly sequence (`docs/seeded_family_definition.md` §1★) |
| **segmental duplication (SD)** | a pair of genomic segments ≥ 1 kb at ≥ 90 % identity (SEDEF/BISER convention; Yoo Supp. Note XXI) | EDGE EVIDENCE of paralogy for *young* copies; no gene or read needed | assembly self-alignment |
| **duplicon** | the ancestral duplication unit, typically sub-gene (DupMasker; Vollger 2022) | a NODE CLASS below the gene | assembly + ancestral-unit library |
| **expansion** | a lineage-specific increase in a family's copy number: the in-paralogs of one outgroup gene, i.e. a lineage-restricted clade of a block | an EVENT / a CLADE; needs an outgroup and an age signal | orthology + identity/dating |

Relations that hold:

- **SD-embedded families ⊂ multi-copy gene families.** Ancient families (olfactory receptors, globins, HOX) are
  multi-copy but far below 90 % identity, so no SD caller sees them and no read is ambiguous there. Soto's and Yoo's
  sets are the **young, high-identity, SD-embedded subset**, which is exactly the subset where MAPQ-0 ambiguity exists
  and where O2/O3 are needed. "Families within SDs" is therefore a *restriction* of the advisor's term, not an
  alternative to it.
- **An expansion produces a family in the lineage where it happened.** LRPAP1: single-copy in human/chimp/orangutan,
  one expansion in gorilla, one gorilla family of 11 loci (8 full + 3 5′ fragments; Yoo counts 10). One expansion,
  one family, same thing seen as event and as result.
- **Expansion ↔ family is not 1:1, and that is the thesis content, not a terminological point.** The unit of duplication
  is the DNA block, not the gene: Yoo's gorilla chr1 unit is MAPKBP1 + JMJD7-PLA2G4B + SPTBN5 co-duplicated ×8 (one
  expansion, three families); PKD1P-NPIP fusions are one path across a duplicon boundary (one family, two blocks;
  `docs/HIERARCHY_DUPLICON_EXPANSION_2026-10-07.md`, [[project_fusions_are_duplicon_blocks]]); Soto's cut follows
  duplicon boundaries and is a cover, not a partition ([[project_cn_cut_follows_duplicons]], [[project_unit_cover]]).
  The pre-registered T1 (duplicon boundary vs family boundary) and T2 (expansions inside families) tests
  (`docs/PREREG_hierarchy_duplicon_expansion_2026-10-07.md`) are the measurements of this non-1:1 relation.
- **Family and subfamily are levels of one tree** (NPIPA/NPIPB: label-pure at every k = 2..8 under average linkage,
  no pairwise threshold separates them; [[project_family_hierarchy]]). The advisor accepted that framing; this is the
  same move one level down.

### What the sources themselves say

- Yoo et al. 2025, Supplementary Note VIII: *"There have been a number of gene family expansions in the NHPs, with
  between 1394-2056 novel gene copies found across the 184-258 families."* Both words in one phrase; "expansion" is
  the event, "family" the set being counted. The thesis author is a contributing author of this note.
- Soto et al. 2025: title *"Human-specific gene expansions contribute to brain evolution"*; body uses "gene family"
  63 times, "gene expansion(s)" 7 times, "paralog(s)" 113 times (count on the vault full text). Their objects are
  "213 human-specific gene families / 1,002 paralogs" defined on SD98 (SDs > 98 % identity) by shared exons + famCN.
- The thesis's own scoring already treats Soto's sets as families (ARI 0.7096 against them; "Soto is a refinement of
  our L2/L3"). One cannot use a catalogue as family truth and deny that its entries are families.

### The sentence to use

> "Yes, they are multi-copy gene families. Specifically they are the young ones whose paralogy is still visible as
> segmental duplications, and that is the subset where reads cannot be placed. The SD, the duplicon and the expansion
> are the layers below the family: the SD is the edge evidence, the duplicon is the sub-gene unit, the expansion is
> the lineage-specific clade. O1 defines the family; the other three are how it got there."

### How to apply

- Concede the noun, keep the structure. Say "same object, I am naming the layer below", never "no, they are not
  families".
- Open definitional discussions with one instance (LRPAP1: one expansion, one family, 11 loci, 3 fragments) and ask
  for his definition first; agreement on an instance ends the argument faster than agreement on definitions.
- In the thesis text: "multi-copy gene family" is the object O1 defines; "SD", "duplicon" and "expansion" are typed
  layers from other evidence (SEDEF; DupMasker; orthology), related to O1 families by T1/T2. Defining SD/duplicon from
  the genome alone remains out of scope (genome-only discovery, user 2026-09-13); letting SD/duplicon structure decide
  membership failed every time it was tried (register rows 676, 677, 817, 818, 997).
- Never pool human (Soto) and gorilla (Yoo) numbers when quoting either catalogue.

Related: `docs/ADVISOR_HIERARCHY_2026-10-03.md`, `docs/seeded_family_definition.md` §1★, `docs/HIERARCHY_DUPLICON_EXPANSION_2026-10-07.md`;
vault `Topics/Terminology.md`; memory [[project_hierarchy_sd_duplicon_expansion_family]], [[reference_advisor_canzar]].


---

## PENDING_2026-10-11.md

*(original: `docs/archive/2026-10/merged/PENDING_2026-10-11.md`)*

## Pending work, 2026-10-11

The genome-wide rescoring started 2026-10-10 is **paused** until the figure set is agreed with the advisor, whose focus is the
StringTie / FLAIR benchmark. The decisions to take with him (samples, annotation-free vs guided, headline framing, truth and
metrics, family figures main vs supplementary) and the cost of each option are in the Claude Doc **"Benchmark Figure Decisions"**
(https://claude.ai/code/artifact/c1333a2e-d5c7-4fbd-a167-0d252c18f2a9). Current objective status: `docs/THESIS_OBJECTIVES.md`,
block "STATUS 2026-10-11".

### Done (run cache `/mnt/linuxdisk/tmp/rustle_figures/runs/<sample>/`, stamps fresh)

| Stage | Samples | Notes |
|---|---|---|
| `assemble`, `assemble_primary` | all six | refreshed for the current binaries |
| `families` (default, genome-wide) | human A119b (2,150 families / 9,907 loci), human testis (342 / 1,115), chimp (396 / 1,343), gorilla KB3781 (326 / 918) | first genome-wide runs; A119b took 196 shards, the chr13 short-arm shards 12-35 min each |
| `families_primary` | human testis | |
| `assign` | chimp (6 min, 5.4 GB), gorilla KB3781 (9 min, 5.3 GB) | |
| figure 8 tables | the four samples above | `figures/data/fig8_*.tsv`, 2026-10-11 |

### Paused or not started

| Work | State | Resume |
|---|---|---|
| O3 `candidate_copies` false flags: are they other apes' transcripts? (2026-10-14) | pre-registered, not run: `docs/PREREG_o3_candidates_foreign_check_2026-10-14.md`. A clean simulation gave 0 few-read false flags in 96 no-deletion units (`docs/PREREG_o3_fewread_flags_2026-10-14.md`), so no rule is fitted on simulation | main machine: align the a14a15 / a15h `iso.contigs.fa` unions to the gorilla primary + mat / pat and to CHM13, mPanTro3, mPonPyg2, mSymSyn1 (light: ~30 sequences; one `minimap2 -c -x splice:hq -uf` per genome) |
| Figure 7 genome-wide family scoring | human guided all-vs-all 180 / 183 shards (shard 180 = chrY ampliconic gene bodies, > 2 h alone); gorilla, chimp, orangutan guided runs not started (about 4-6 h each) | `python3 figures/make.py data fig7 --set fig7_budget_s=<seconds>` with `MM2_SHARD_BUDGET_S` / `MM2_SHARD_DEADLINE` set to match |
| `assign` human A119b | stopped by Claude Code's low-memory reaper (background job, idle session) | rerun in a session started with `CLAUDE_CODE_DISABLE_BG_SHELL_PRESSURE_REAP=1`, or only when asked |
| `flag` (O3 scans) KB3781, A119b, chimp | not started; 16-19 GB each; chimp needs its splice index (21 GB peak) | `make.py runs --sample S --stage flag` |
| `families` gorilla OR6737, orangutan; `families_primary` for five samples | not started | `make.py runs --sample S --stage families` |
| Figures 1 and 3 genome-wide for human; figure 2 genome-wide (SQANTI3, 13-25 h) | not started | `make.py data fig1|fig3|fig2` |
| Generalization pre-registration (species x tissue x family class, one species held out) | offered, not written | after the figure set is agreed |

### Machine notes

- `tools/mm2_shard.sh` needed its executable bit (097d06ad); the families stage failed with EACCES before.
- `/mnt/linuxdisk` filled to 5 MB free on 2026-10-10. Freed: 40 GB of old shard caches (user approved: the 09-29 o1_cover held-out runs, the 10-01
  rna_allele per-chromosome and heldout_cat runs, the 09-29 defaults_flip runs) and each finished sample's cache. About 25 GB free on 2026-10-11; the
  human guided cache (`mm2_shard_cache/a9578e75...`) is kept so figure 7 resumes at shard 180.
- Shards longer than the 10-minute call limit forced 2-hour calls in the background, one heavy job at a time.
- 2026-10-13: `/mnt/linuxdisk` was down to 3.2 GB free. The user approved removing the O1 false-merge study's leftovers, and only the parts no code references were deleted from `home/juanfraitu/o1_falsemerge`: about 67 GB, including `lit/phap`, `lit/evid/GGO` and `PTR`, `lit/sharedef_*`, `lit/nodemcl_*`, `lit/oracle_*`, `rebuild*` and `gw_atoms`. Kept because code still reads them: `lit/evid/PPY` (orangutan index for bench/hierarchy), `lit/aj_ho`, `lit/aj_dev`, `lit/annot_gencode` (cat_layer_order), `human2` and `step4` (layer_order). About 55 GB free afterwards. Other candidates, not approved and still in place: unreferenced `tmp/rustle_figures_dev` runs (~55 GB), `o3_excise` plus old subset BAMs in `_from_wsl/winloci_scratch` (~32 GB), and parked `tmp/rna_allele` runs (~33 GB).

### Milestone (deferred, decide later): reaching StringTie's remaining transcripts without flow

Recorded 2026-10-11 after the assembler work against StringTie (register 1264-1266; `docs/PREREG_direct_transfrags_union_2026-10-11.md`).
**Where we are.** `--direct-transfrags-union` (opt-in) is above StringTie in human on intron, chain, transcript and locus (both measures)
and in gorilla on chain and transcript sensitivity and locus; it costs 0.8-1.5 points of intron-chain precision against our default
(why it is not the default). StringTie still leads on BASE sensitivity in both species (gorilla 50.3 vs 46.9, human 59.6 vs 49.2) and,
in gorilla, on intron precision (83.1 vs 80.2).
**What remains (dev, gorilla NC_073244.2, 195 reference transcripts StringTie matches and pure direct transfrags does not):** 66 restored
by the union; 71 carried whole by one read; 29 by one exact read among longer ones; 16 exact >= 2 but dropped by gate or polish; 4 only
inside longer reads; **9 (5%) flow-only** (no read carries the whole chain).
**Candidate levers, not started:**
1. **Non-flow extension of truncated maximal transfrags** to a source / sink along the best-supported compatible continuation (StringTie's
   long-read back / forward extension is itself greedy, not flow). Ceiling: the flow-only class, about 5% of the remaining gap.
2. **Transcript ends** (what base sensitivity measures): where our boundary sits against StringTie's on matched transcripts; the 5′ end is a
   known dispersion problem (an extension rule failed before, §6w3), the 3′ end has the opt-in `--polish-tes`.
**Worth it?** Decide with the advisor's figure choice (`Benchmark Figure Decisions`): only if the benchmark figure is to beat StringTie on
every level. Each lever: pre-register, dev NC_073244.2, fresh held-out (KB3781 / chimp / orangutan are now spent for admission rules, so
use human testis and another gorilla chromosome set, or a new sample).

#### Both levers built as opt-ins (user 2026-10-11, after the FLAIR opt-ins)

| flag (pipeline variable) | what it does |
|---|---|
| `--extend-transfrags` (`RUSTLE_EXTEND_TRANSFRAGS=1`) | every admitted chain is ALSO emitted extended one junction at a time along the junction pair the most reads link (>= the read floor), upstream and downstream, ending at an observed read start / end; the first step needs the chain's own reads to stay out of the new intron. StringTie's greedy back / forward extension, no flow; it builds chains no single read carries |
| `--transfer-ends` (`RUSTLE_TRANSFER_ENDS=1`) | a chain's start / end also count the reads that share its first / last junctions and fit inside it (a splice graph's first / last exon spans every read through it) |

Default byte-identical on both contigs; unit tests in `src/family.rs`. **Dev effect (descriptive, not pre-registered;** primary pool, default
polish; sensitivity / precision; "StringTie-only" = reference transcripts StringTie matches exactly and the union does not):

| arm | gorilla NC_073244.2: n | base | chain | StringTie-only (of 133) | human chr20: n | base | chain | StringTie-only (of 88) |
|---|---|---|---|---|---|---|---|---|
| default | 3,862 | 44.3 / 77.1 | 28.3 / 40.8 | — | 5,522 | 51.2 / 23.6 | 24.7 / 21.0 | — |
| union | 4,169 | 46.9 / 76.3 | 29.0 / 38.7 | 0 | 6,068 | 52.7 / 24.4 | 25.6 / 19.5 | 0 |
| union + extend | 4,287 | 46.9 / 76.3 | 29.3 / 38.0 | 6 | 6,681 | 52.9 / 23.7 | 26.2 / 18.2 | 12 |
| union + transfer ends | 4,169 | 47.2 / 76.0 | 29.0 / 38.7 | 0 | 6,018 | 53.4 / 24.4 | 25.5 / 19.7 | 0 |
| union + both | 4,280 | 47.2 / 76.1 | 29.3 / 38.1 | 6 | 6,549 | 53.4 / 23.7 | 26.0 / 18.5 | 12 |
| default + extend | 3,984 | 44.7 / 77.1 | 28.7 / 40.1 | 11 | 6,050 | 51.6 / 21.9 | 25.3 / 20.0 | 15 |
| default + transfer ends | 3,862 | 44.5 / 76.9 | 28.3 / 40.8 | 5 | 5,502 | 51.6 / 23.6 | 24.7 / 21.1 | 5 |
| StringTie | 3,735 | 50.8 / 73.8 | 24.7 / 37.0 | | 5,731 | 60.5 / 17.7 | 20.1 / 16.8 | |

(intron / transcript / locus in `/mnt/linuxdisk/tmp/asm_gap/st_{g20,h20}/g_*.stats`; script scratchpad `st_arms.sh`.)
**Reading.** (1) Extension reaches 6 / 12 of StringTie's exclusive transcripts for -0.7 / -1.3 chain precision, about what the 5% flow-only
ceiling predicted: StringTie's remaining chain lead is not flow. (2) StringTie's BASE lead is not transcript ends. On reference transcripts
both tools match exactly, our terminal exons already cover MORE of the reference terminal exons than StringTie's (gorilla 79.3% vs 78.4%,
human 74.0% vs 68.9%; `--transfer-ends` raises ours to 80.3% / 78.9%). Of the reference exonic bp StringTie covers and the union does not
(gorilla 355 kb), 37% lie inside our transcripts' spans (exons of isoforms we do not output), 33% in reference genes where we output
nothing, 30% outside our spans in genes where we output something (other isoforms' terminal exons): volume, i.e. more isoforms and
low-coverage loci, not ends (scratchpad `base_gap.py`, `ends_cmp.py`). Closing it is the single-read / no-polish trade already measured
for FLAIR.

**StringTie's own default values (user, same day).** `-E 25` and `-a 10` as our nearest flags (`--read-junction-snap 25`: our rule, a junction
moves onto one within 25 bp that more reads carry, not StringTie's exact rule; `--terminal-exon-anchor 10`: first / last exons only):

| arm | gorilla: n | base | intron | chain | StringTie-only (of 133) | human: n | base | intron | chain | StringTie-only (of 88) |
|---|---|---|---|---|---|---|---|---|---|---|
| union | 4,169 | 46.9 / 76.3 | 58.8 / 81.2 | 29.0 / 38.7 | 0 | 6,068 | 52.7 / 24.4 | 65.1 / 60.6 | 25.6 / 19.5 | 0 |
| + snap 25 | 3,706 | 47.2 / 75.9 | 57.8 / 82.6 | 25.8 / 38.7 | 28 | 5,461 | 52.6 / 25.4 | 63.8 / 61.3 | 23.0 / 19.5 | 8 |
| + anchor 10 | 4,169 | 46.9 / 76.3 | 58.8 / 81.2 | 29.0 / 38.7 | 0 | 6,061 | 52.7 / 24.4 | 65.1 / 60.7 | 25.6 / 19.6 | 0 |
| + snap 25 + anchor 10 + extend + ends | 3,793 | 47.5 / 75.7 | 57.8 / 82.6 | 25.9 / 38.0 | 30 | 5,860 | 53.5 / 24.8 | 63.9 / 61.4 | 23.3 / 18.4 | 19 |
| + snap 25 + anchor 10 + single-read chains | 6,171 | 53.7 / 67.5 | 65.9 / 69.8 | 29.3 / 26.5 | 85 | 8,270 | 58.1 / 21.9 | 69.6 / 46.8 | 24.8 / 13.9 | 30 |
| StringTie | 3,735 | 50.8 / 73.8 | 58.9 / 83.9 | 24.7 / 37.0 | | 5,731 | 60.5 / 17.7 | 63.3 / 57.9 | 20.1 / 16.8 | |

**Reading.** StringTie's values do not help us: the 25 bp snap costs 3.2 / 2.6 points of chain sensitivity with no precision gain, the 10 bp
anchor is inert. But the snap makes us LOOK like StringTie (gorilla 3,706 vs 3,735 transcripts, chain 25.8 vs 24.7, intron 57.8 / 82.6 vs
58.9 / 83.9) and reaches 28 of its exclusive gorilla matches while losing more of ours: a 25 bp correction merges real alternative splice
sites, which is plausibly part of why StringTie's chain sensitivity is below ours (our rule, not StringTie's exact one). No code change.

**How StringTie chooses single-read chains (user, same day; scratchpad `single_read_choice.py`, dev contigs).** A single-read chain = an
intron chain exactly one primary read carries after our de-duplication key. Gorilla NC_073244.2: 16,156 such chains; StringTie emits 826
of them (reference-exact 14.2%), our `--min-chain-reads 1` emits 3,544 after polish (8.8%). Human chr20: 37,488; StringTie 1,650 (3.7%),
ours 4,064 (3.3%). The feature that separates them: every junction of the chain is also carried by ANOTHER read — 97.3% of StringTie's
kept chains vs 22.6% of our extras in gorilla, 67.3% vs 43.5% in human. Applied to the chains we emit, that rule keeps 1,162 at 14.9%
(gorilla; 173 reference-exact vs StringTie's 117) and 2,004 at 4.5% (human; 91 vs 61): StringTie's precision on these chains, with more
of them right. The rule is our existing `--read-isoform-k 2` (register 848 / 855: refuted against OUR default, as a family-recovery lever
and as an operating point). End to end:

| arm | gorilla: n | base | intron | chain | StringTie-only (of 133) | human: n | base | intron | chain | locus |
|---|---|---|---|---|---|---|---|---|---|---|
| union | 4,169 | 46.9 / 76.3 | 58.8 / 81.2 | 29.0 / 38.7 | 0 | 6,068 | 52.7 / 24.4 | 65.1 / 60.6 | 25.6 / 19.5 | 38.6 / 31.1 |
| union + read-isoform-k 2 | 5,110 | 48.6 / 73.8 | 61.1 / 78.9 | 31.3 / 34.1 | 63 | 7,207 | 54.8 / 22.9 | 67.7 / 57.3 | 26.7 / 17.2 | 39.9 / 29.1 |
| default + read-isoform-k 2 | 5,046 | 47.2 / 74.4 | 59.0 / 80.5 | 31.0 / 34.2 | 71 | 7,010 | 54.1 / 21.3 | 65.5 / 60.3 | 26.0 / 17.6 | 39.1 / 26.4 |
| union + single-read chains | 7,265 | 53.8 / 66.9 | 67.7 / 66.9 | 33.7 / 25.9 | 75 | 9,293 | 58.5 / 21.9 | 71.6 / 45.5 | 27.8 / 13.8 | 41.3 / 27.3 |
| StringTie | 3,735 | 50.8 / 73.8 | 58.9 / 83.9 | 24.7 / 37.0 | | 5,731 | 60.5 / 17.7 | 63.3 / 57.9 | 20.1 / 16.8 | 36.3 / 25.7 |

**Reading.** A rule without flow reproduces StringTie's single-read choice: `union + --read-isoform-k 2` reaches about half of StringTie's
exclusive gorilla matches (63 / 133) and keeps most of our precision (gorilla chain 31.3 / 34.1 vs StringTie 24.7 / 37.0; human 26.7 / 17.2
vs 20.1 / 16.8, above StringTie on both). It is not a default: against our own union it costs 4.6 / 2.3 points of chain precision (the
B2 bar is 0.5), the same verdict as 848 / 855. It is the StringTie-matched operating point if the benchmark wants one; that choice needs its
own pre-registration on a fresh sample. Note: `--read-isoform-k` runs on the buffered reader (not streaming), so it is slower genome-wide.
**Necessary, not sufficient.** Of the single-read chains whose junctions are all shared, StringTie keeps 804 / 6,901 (11.7%, gorilla) and
1,110 / 16,964 (6.5%, human); restricted to chains no longer read contains and whose every junction is >= 10% of its strongest competitor,
25.4% / 17.5%. Which of the eligible chains StringTie keeps is decided by its flow-based path selection and coverage, which no simple read
feature reproduced; the shared-junction rule matches the PRECISION of its choice (above), not the choice itself.

### Milestone (deferred, decide later): what FLAIR reaches that the union does not

Recorded 2026-10-11 (dev contigs only: gorilla OR6737 NC_073244.2, human A119b chr20; `--direct-transfrags-union` vs FLAIR 3.0.1 as run by
`benchmark_collapse/run_flair.sbatch`: `flair collapse --trust_ends --generate_map`, support 3 by default, no `flair correct`, no read
de-duplication, support counted after re-aligning reads to its isoforms). FLAIR matches **89 / 147** reference multi-exon transcripts the
union does not; the union matches **309 / 218** FLAIR does not.

| class (gorilla chr20 / human chr20) | count | what would reach it | worth it? |
|---|---|---|---|
| we assemble it, polish removes it | 14 / 53 | a looser polish | no: the polish dials were fitted for precision and no polish rule had positive value (register 853-873, 1066) |
| >= 2 exact reads but all at ONE coordinate (PCR duplicates under our key) | 0 / 14 | `--keep-coordinate-duplicates` | no: identical-coordinate reads add no independent evidence (1058, 1065, 1066); FLAIR's support there is duplicate-inflated |
| a non-canonical junction | 2 / 9 | relax `--assembly-junctions strict` | no: non-canonical junctions are mostly not real (849) |
| 1 exact read + reads that fit inside it, starved by the proportional split of direct transfer | **30 / 21** | credit each contained read to every maximal container it fits (admission only), or assign it to one container as FLAIR's alignment does | **candidate**: the one lever with a plausible precision budget |
| 1 exact read + near-identical reads (all junctions within 10 bp) | **11 / 11** | read-level junction correction before grouping (StringTie `-E` / `flair correct`) | candidate, small; post-assembly snaps were lossless but tiny (1140) |
| truly one read in our alignments | 32 / 36 | admit single-read chains | no: that is FLAIR's precision cost (gorilla chain precision 17.6 vs our 34) |
| other | 0 / 3 | — | — |

**Reading.** FLAIR's lead is mostly volume — single reads, duplicates, no polish — bought at half our precision; at most about a third of its
exclusive matches (41 / 32) is reachable with a principled change. Same rule as the StringTie milestone: pre-register, dev contigs, a fresh
held-out sample. Scripts: scratchpad `flair_only.py`, `never_assembled.py` (to be moved into `bench/` if the milestone is taken up).

#### FLAIR's levers as opt-ins (user 2026-10-11: "include those options on our side too, as an opt in")

Three new `copy_assign` flags (each default = off, byte-identical; unit tests in `src/family.rs`), alongside the existing
`--keep-coordinate-duplicates`, `--assembly-junctions majority` and `--assembly-polish none`:

| flag (pipeline variable) | what it does | FLAIR class it targets |
|---|---|---|
| `--transfer-credit full` (`RUSTLE_TRANSFER_CREDIT=full`) | with `--direct-transfrags[-union]`: a contained transfrag's reads count fully toward EVERY maximal container (admission only) instead of a proportional split | starved by the split |
| `--read-junction-snap W` (`RUSTLE_READ_JUNCTION_SNAP=W`) | before grouping, a read junction within W bp (donor and acceptor) of one carried by more reads (>= 2) moves onto it, to a fixed point | near-identical reads |
| `--min-chain-reads N` (`RUSTLE_MIN_CHAIN_READS=N`) | read floor for spliced chains (1 = single-read chains); unspliced reads keep the floor 2 | truly one read |

**Dev effect (descriptive, NOT a pre-registered test; no adoption decision).** Dev contigs only, primary pool, all with the default
polish unless noted; gffcompare sensitivity / precision; "FLAIR-only" = reference transcripts FLAIR matches exactly and the union does not
(gorilla 92, human 163 on these runs, single-exon included), and how many each arm recovers.

| arm | gorilla NC_073244.2: n | intron | chain | FLAIR-only | human chr20: n | intron | chain | FLAIR-only |
|---|---|---|---|---|---|---|---|---|
| default | 3,862 | 54.8 / 83.3 | 28.3 / 40.8 | 3 | 5,522 | 61.8 / 65.1 | 24.7 / 21.0 | — |
| union | 4,169 | 58.8 / 81.2 | 29.0 / 38.7 | 0 | 6,068 | 65.1 / 60.6 | 25.6 / 19.5 | 0 |
| union + credit full | 7,267 | 60.8 / 65.5 | 31.8 / 24.3 | 51 | 10,469 | 69.1 / 44.1 | 28.3 / 11.7 | 70 |
| union + snap 10 | 3,847 | 58.1 / 82.1 | 26.7 / 38.6 | 4 | 5,716 | 64.5 / 60.9 | 24.2 / 19.5 | 3 |
| union + min chain reads 1 | 7,265 | 67.7 / 66.9 | 33.7 / 25.9 | 38 | 9,293 | 71.6 / 45.5 | 27.8 / 13.8 | 24 |
| union + keep duplicates | 3,743 | 58.6 / 83.6 | 27.8 / 41.4 | 3 | 6,230 | 65.0 / 56.3 | 24.1 / 18.1 | 12 |
| union + credit full + snap 10 | 6,400 | 60.0 / 68.0 | 28.9 / 25.1 | 46 | 9,617 | 68.1 / 45.4 | 26.5 / 11.9 | 57 |
| union + all four | 7,375 | 66.2 / 65.0 | 30.5 / 23.0 | 50 | 10,706 | 71.1 / 39.5 | 26.4 / 10.7 | 56 |
| union + all four, no polish | 15,720 | 67.7 / 57.4 | 34.1 / 12.7 | 75 | 31,431 | 75.2 / 32.6 | 32.0 / 5.3 | 108 |
| FLAIR | — | 55.6 / 67.3 | 25.1 / 23.6 | | — | 67.9 / 33.1 | 23.9 / 7.6 | |

**Reading.** Full credit and single-read chains reach about half of FLAIR's exclusive matches and land at FLAIR's precision (gorilla
chain 24-26 vs FLAIR 23.6; human 12-14 vs 7.6), i.e. they buy FLAIR's volume at FLAIR's price while staying above it on chain
sensitivity and precision. Snap at 10 bp LOWERS chain sensitivity in both species (-2.3 / -1.4 vs the union): within 10 bp it merges real
alternative splice sites; a smaller window or a motif-aware snap would have to be pre-registered. Keeping duplicates is neutral in
gorilla and costs 4 points of intron precision in human. With every lever on and no polish, intron and chain sensitivity exceed FLAIR's in
both species but precision falls below it (chain 12.7 vs 23.6 gorilla, 5.3 vs 7.6 human): the polish is what keeps our precision. None of these is a default; a FLAIR-matched benchmark arm, if the advisor wants
one, would be `union + credit full` or `union + min chain reads 1`, chosen by pre-registration on a fresh sample.
Script: scratchpad `opt_arms.sh` (outputs `/mnt/linuxdisk/tmp/asm_gap/opt_{g20,h20}/`).

**Second batch (same day, user: "please proceed" with FLAIR's two remaining mechanisms).** Three more opt-ins, default byte-identical
(checked on both contigs):

| flag (pipeline variable) | what it does |
|---|---|
| `--transfer-credit best` (`RUSTLE_TRANSFER_CREDIT=best`) | each contained transfrag's reads go wholly to ONE maximal container (FLAIR gives each read one isoform): the one with the fewest exonic bp, then the most own reads |
| `--check-splice` (`RUSTLE_CHECK_SPLICE=1`) | FLAIR `--check_splice`: a spliced read supports nothing unless every junction has >= 4 of the 6 exonic bases around it aligned and no insertion over 3 bp among them |
| `--terminal-exon-anchor N` (`RUSTLE_TERMINAL_EXON_ANCHOR=N`) | the end rule of FLAIR `--stringent` (N = 25 there): first and last exons each >= N bp. FLAIR's other `--stringent` condition (80% of the isoform covered) is not built: under it a contained read never supports a longer isoform, which is our default already |

The two read filters run in the streaming `--assemble-only` reader only (the buffered path refuses them). They drop 3,642 (check-splice) and
3,418 (anchor) spliced reads on gorilla NC_073244.2, 6,817 together; 2,414 / 8,942 / 10,917 on human chr20.

| arm | gorilla: n | intron | chain | FLAIR-only (of 92) | human: n | intron | chain | FLAIR-only (of 163) |
|---|---|---|---|---|---|---|---|---|
| union | 4,169 | 58.8 / 81.2 | 29.0 / 38.7 | 0 | 6,068 | 65.1 / 60.6 | 25.6 / 19.5 | 0 |
| union + credit best | 4,616 | 59.4 / 78.7 | 29.6 / 35.6 | 17 | 6,772 | 66.2 / 57.4 | 25.9 / 17.6 | 17 |
| union + check-splice | 4,124 | 58.7 / 81.4 | 29.0 / 39.1 | 0 | 6,026 | 65.1 / 60.9 | 25.5 / 19.7 | 0 |
| union + anchor 25 | 4,104 | 58.5 / 81.5 | 28.7 / 38.9 | 0 | 5,978 | 64.7 / 60.9 | 25.5 / 19.8 | 1 |
| union + min1 | 7,265 | 67.7 / 66.9 | 33.7 / 25.9 | 38 | 9,293 | 71.6 / 45.5 | 27.8 / 13.8 | 24 |
| union + min1 + check-splice | 7,144 | 67.6 / 67.4 | 33.6 / 26.2 | 38 | 9,208 | 71.6 / 45.9 | 27.7 / 13.9 | 24 |
| union + min1 + anchor 25 | 7,190 | 67.4 / 67.1 | 33.6 / 26.0 | 35 | 9,223 | 71.4 / 45.5 | 27.6 / 13.8 | 24 |
| union + min1 + both filters | 7,074 | 67.4 / 67.6 | 33.5 / 26.4 | 35 | 9,129 | 71.3 / 45.9 | 27.6 / 14.0 | 24 |
| union + credit full | 7,267 | 60.8 / 65.5 | 31.8 / 24.3 | 51 | 10,469 | 69.1 / 44.1 | 28.3 / 11.7 | 70 |
| union + credit full + both filters | 6,929 | 60.6 / 66.8 | 31.5 / 25.3 | 49 | 10,130 | 68.9 / 45.1 | 28.3 / 12.1 | 70 |
| union + credit best + min1 + both filters | 7,363 | 67.5 / 66.4 | 33.5 / 25.4 | 41 | 9,279 | 71.4 / 45.2 | 27.5 / 13.4 | 33 |

**Reading.** (1) Single assignment (`best`) is the middle ground it was meant to be (chain precision -3.1 / -1.9 against the union vs
-14.4 / -7.8 for `full`) but reaches only 17 of FLAIR's exclusive transcripts in each species (`full` 51 / 70): the shortest container is
usually not the isoform FLAIR emits. Sensitivity bought per point of precision is about the same as `full` (gorilla) or lower (human).
(2) FLAIR's two read filters are nearly free (chain sensitivity -0.0 to -0.3) and nearly inert (chain precision +0.1 to +1.0); they do
NOT make single-read chains cheap (gorilla 25.9 -> 26.4, human 13.8 -> 14.0). The precision cost of single-read chains is not in
badly aligned junctions or short terminal exons. (3) With this batch every FLAIR mechanism we know of has an opt-in on our side; the
decomposition above has no class left without a lever.

### Multi-copy loci vs StringTie / FLAIR: where we stand (2026-10-12, for later)

Registers 1267-1270; preregistrations `PREREG_tool_difference_model_2026-10-11.md`, `PREREG_inbetween_operating_point_2026-10-11.md`,
`PREREG_s1_o2_certificate_2026-10-12.md`, `PREREG_o2_pruning_and_eichler_2026-10-12.md`; scripts `bench/tool_difference/`.

- **Operating point** `--direct-transfrags-union --read-isoform-k 2` (opt-in): above both tools on intron-chain sensitivity, above
  FLAIR on precision, in both species (1268). StringTie keeps a higher chain precision in gorilla.
- **Neither tool uses secondary alignments** (no StringTie or FLAIR transcript rests on secondary-only coverage).
- **Secondary-only genes (S1) are NOT a win**: O2's certificate puts their reads on the other copy; exactly 1 S1 gene per species is
  O2-supported (1269). Most S1 transcripts are shadows, often of a spliced parent gene on its processed retrocopy.
- **Coin-toss genes (S2)**: after O2 pruning we still match more reference transcripts than both tools (gorilla 72 vs 51 / 42,
  human 119 vs 77 / 72); precision above both in human, a few points below in gorilla (1270).
- **O2 vs Eichler's AS-margin rule**: they agree on 14-54% of calls; AS rewards unspliced alignments to retrocopies (no splice
  penalty), the PSV bases do not — O2 is right on that population (1270).
- **Open**: (1) copy-of-origin truth from the family simulator (`bench/famsim`) to score pruning and assignment directly; (2) a
  pruning rule finer than "gene contradicted" where copies are identical (human S2: 64% of records have no distinguishing column);
  (3) DONE (registers 1271-1272): O2 moves 1.0 / 2.8% of multi-placement reads (519 gorilla, 14,945 human); they mostly land on genes
  the tools already assemble or on unannotated sequence; only 25 / 44 genes with >= 2 moved reads lack a StringTie / FLAIR transcript
  (mostly ncRNA / LOC). Next lead: cluster the unannotated landing sites (10,248 human reads) into candidate unannotated copies (O3).

### O2 follow-ups (added 2026-10-13)

Context: the real-read truth tests of 10-13 (`docs/PREREG_o2_anchor_2026-10-13.md`, arms 1-4; register 1274-1275) and
the two O2 fixes they forced: exact ties decided by rounding (7ac557e7), and the tie-outside gate's region order
(014126be + d5d982a3).

| Work | State | What it needs |
|---|---|---|
| **Arm 4 (haplotype truth) on a second individual: HG002** (user: "add that to the list", not yet). The best O2 truth source: a read must match its own locus on mat or pat. Gorilla gave 28/29 right but was UNDERPOWERED and MCL26-dominated; HG002 is the human held-out | not started | **Present** in `/mnt/linuxdisk/tmp/rustle_figures_dev/o2_presence/`: `hg002/hg002v1.1.fa` (Q100 v1.1, mat + pat, 6.1 GB); `chain/` CHM13 <-> HG002 mat/pat chains; `dl/chunk0-1.bam` = the first 8 GB of the 58.9 GB PacBio HG002 Kinnex FLNC BAM. **RNA: downloaded 2026-10-13** to `/mnt/linuxdisk/home/juanfraitu/hg002_kinnex/flnc.bam` (63.2 GB, byte count as expected after a resumed download, `samtools quickcheck` OK, full read-through 37,200,651 records = the poly-A FLNC of `flnc.filter_summary.json`; `.pbi` and `flnc.report.csv` alongside; `/mnt/linuxdisk` now has ~440 GB free). **CHM13 mapping: on the user's cluster (decided 2026-10-13)**, with the exact O2 settings: `samtools fastq` -> `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes` against CHM13 v2.0 (md5 42b830df45df08ab63068e9128156e56; minimap2 2.30) -> sorted, indexed BAM, ~130 GB. A local probe on 20,000 reads (`/mnt/linuxdisk/tmp/hg002_probe/`) gave 250 reads/s here (about 41 h for all reads), 99.9% mapped, and 0.29% AS-tied reads (about 100k genome-wide). **Missing:** the cluster BAM; mat / pat splice indexes (~13-15 GB each; build on `/mnt/linuxdisk`); the catalog decision (HG002's own vs A119b's); an HG002 catalog (run the default pipeline on the HG002 RNA, or project the CHM13 A119b catalog through the chains); then `bench/o2_truth/haplotype_truth.py` unchanged, with its own prereg (bars H1-H3 as for gorilla). |
| O2 per-family blind spot: the true copy sits in another family | known (o2_read_truth, register 1275 arm-4 error) | the shipped `--union-certificate` ties everything (r1103); needs a decision with the advisor |
| Sole-candidate promotion of a read whose real origin is uncatalogued (arm 4: 1 molecule assigned in 2 families, NM 0 at an uncatalogued locus) | found 10-13 | a candidate gate (e.g. no promotion without decisive evidence), with its own prereg |
| Tie-outside gate residual: GTF routing under `--families --gtf` still reads the in-worker status | documented in `docs/PREREG_o2_tie_outside_gate_2026-10-13.md` | not a pipeline combination; fix only if that mode is used |
| Re-run every O2 number and figure on a binary >= 7ac557e7 | human contested assigned was -19% after the tie fix | before any O2 figure is shown |
| Arm 2 done soundly (internal introns only; no canonical pair within +-3 bp in the copy's genome) | arm 2 truth invalid as built | a fresh substrate (chimp) |
| Catalog-PSV-only read-star columns as an O2 option (popfix post hoc: 125 vs 3 YAG assignments) | untested | simulation truth first |
| Copy-group assignment under a flag (population sub-project A) | deferred by the user | after sub-project B (done 10-13) |
| Branch `machine2/soto-evidence` | DONE 2026-10-13: merged into `main` (fast-forward), pushed, then deleted locally and on GitHub (user); GitHub `main` = 49e62be3 | nothing |

### O1 two-node fix (2026-10-13)

- The shipped `runs/*` catalogs (`/mnt/linuxdisk/tmp/rustle_figures/runs/<sample>/`) predate the two-node MCL fix (`docs/PREREG_mcl_two_node_2026-10-13.md`), but under the default pair rule `exact-mcl` (user decision 2026-10-13) a families-stage run writes the same `clusters.tsv` and `loci.tsv` (verified on all four catalogs); `params.tsv` only gains the `pair_rule` / `pairs_joined_by_rule` rows. The paused figure rescoring needs no re-clustering for this fix; only the opt-in `--pair-rule weight-blind` changes the families (re-clustered copies in `/mnt/linuxdisk/tmp/o1_close/recluster/after/`).
- DECIDED 2026-10-13 (user): the default pair rule is `mcl_families --pair-rule exact-mcl`; `weight-blind` is opt-in. `exact-mcl` reproduces the four shipped catalogs' `clusters.tsv` and `loci.tsv` byte for byte. `params.tsv` ends with `pair_rule` and `pairs_joined_by_rule`: weight-blind 185 / 11 / 43 / 48 (A119b / testis / chimp / gorilla; 68 / 9 / 20 / 15 new 2-member families plus 117 / 2 / 23 / 33 pairs folded into one locus), exact-mcl 0 in all four.
- Uniform cliques below weight 1 (complete graphs of >= 3 nodes whose edges all carry one weight < 1) can still split by float residue: an independent port splits a uniform triangle at 39 and a uniform K4 at 58 of the 500 weights 0.500..0.999 (`docs/PREREG_mcl_two_node_2026-10-13.md` Amendment 1(b)). The four default catalogs hold 5 such components / 15 members (A119b 3 / 9, gorilla 2 / 6), all of which came out whole. Not fixed.
- Phase-B prerequisites (spec amendment 7, `docs/superpowers/specs/2026-10-13-close-o1-design.md`): match over distinct member sets; use the within-family null for E3; base V2's bar on `expansion_not_family`; stratify cross-catalog bars by mapped set size.

### O1 layers, phases 1-2: sedef-rs SD and duplicon layers (2026-10-14)

- **sd-scan.** `utilities sd-scan` (vendored sedef-rs, `src/sedef/`) and the driver's `sd` stage. Caches:
  `/mnt/linuxdisk/home/juanfraitu/sd_cache` (CHM13 v2.0, gorilla, CHM13 v1.0; prefixes in
  `/mnt/linuxdisk/tmp/o1_layers/*.sd.prefix`). `mcl_families --emit-hierarchy --sd-cache PREFIX` reads it.
- **Verdicts:** V1 FAIL (row 1278); V1b FAIL (1279); V2 FAIL (1280: V2a FAIL, held-out P1 0.5153 against 0.80; V2b PASS,
  12 against a bar of 9). Prereg: `docs/PREREG_sedef_rs_layers_2026-10-13.md`. Readings applied: the SD layer keeps
  SEDEF where it exists (V1); the Soto layer takes SEDEF SD98 (V1b); human expansions stay on DupMasker and the apes have
  no expansion layer (V2).
- `mcl_families --emit-hierarchy --sd-cache` builds its SD layer from every sd-scan `ok` pair (no identity or length floor, like the SEDEF input, which also carries pairs below 0.90); since V1 failed, SEDEF stays the SD source where it exists, and the A119b SD-block and complex counts from Task 9 (8,488 / 920 vs 1,548 / 571) are not like-for-like with SEDEF.
- Observation (not re-evaluated here): plan 1's hierarchy evaluation arms (`/mnt/linuxdisk/tmp/o1_close/hier/`, register 1277) were clustered with the then-default weight-blind pair rule (human A119b: 10,043 members / 2,218 families); the exact-mcl default gives 9,907 / 2,150 (`hier_shipped/`).
- **Next:** phase 3 of `docs/superpowers/specs/2026-10-13-o1-layered-hierarchy-design.md` (the Rust Soto core, S1; assembly
  copy number, CN1; S2). Its SD98 source follows V1b's reading (SEDEF SD98; alignment-based identity refinement in
  sedef-rs is its own task). Its plan is written from these results.
  - Consequences for the spec's section 3: with V2 failing, the apes have no expansion layer, so the Yoo expansion bar
    (scored at the expansion layer on apes) cannot be scored until an expansion source passes; with V1/V1b failing, a
    genome without SEDEF has no SD layer and no SD98 source. On disk only human (HSA) and gorilla (GGO) SEDEF tables
    exist, so chimp (section 3's dev ape) has neither; the current scope is human and gorilla (user, 2026-10-13).
  - For the identity-refinement task: on CHM13 v2.0, 107,038 of 186,820 `ok` pairs (57%) have est_id < 90 while every `ok`
    pair has stage2_id >= 90; V1 filtered on est_id (the registered k-mer lower bound), so part of the recall gap may be
    the identity estimate rather than detection. Register 'recall without the est_id filter' as a known-before item of
    the identity-refinement prereg before computing it.

### Making our sd-scan replace SEDEF and DupMasker (user, 2026-10-14; for later)

Asked: "what can we do to improve our own rs-sedef to match SEDEF and dupmasker and use them instead of something
external". Parked by the user (meeting) on 2026-10-14; nothing below has been run. Each step gets its own
pre-registration, with CHM13 dev chromosomes for development and the other chromosomes plus gorilla held back; the
V1/V1b/V2 verdicts above stand.

**Why replace them.** The real SEDEF takes 17-40 h on this machine and WSL2 rebooted mid-run (register 991); BISER
segfaults under WSL2 (register 928). sd-scan runs a whole genome in 18-55 min.

**SEDEF gap (V1 recall about 0.40, precision 0.92-0.95).** Cheapest and most diagnostic first:

| # | Change | Why (evidence) |
|---|---|---|
| 0 | Diagnostic, registered first: recall with `stage2_id` or with no identity filter | 57% of `ok` pairs have est_id < 90, yet all have stage2_id >= 90; V1 filtered on est_id. Separates the identity estimate from detection. |
| 1 | Real alignment identity: banded alignment of each refined pair, giving matches over aligned bases (SEDEF's fracMatch, column 21); reuse Rustle's semiglobal engine | est_id is a k-mer lower bound; this is also the follow-up queued after V1b failed. |
| 2 | Keep high-copy seeds: down-sample k-mers in > 100 windows instead of dropping them | Recall is low even for long, >= 99%-identity SDs (chr9 0.36, chr15 0.32): high-copy cores lose their seeds (sandbox README caveat 5). Sharded counting now keeps memory bounded. |
| 3 | Chain inverted copies (anti-diagonals) | Inverted SDs come out as 2 kb pieces, each losing its edges. |
| 4 | Recover SD ends: overlapping windows (1 kb stride) or extending refined pairs past block edges | Recall is lowest for 1-5 kb SDs (0.37 CHM13, 0.22 gorilla). |
| 5 | Match SEDEF's call semantics where they differ: seeds only in unmasked (uppercase) sequence, plus its final filters | SEDEF behaviour recorded in register 659-662. |

**DupMasker gap (V2 held-out 0.52).** Two routes:
- **Route B (recommended): rebuild DupMasker's function with its own data.**
  - DupMasker annotates a genome against a published library of about 11k ancestral duplicons (`dupliconlib.fa`,
    shipped with RepeatMasker; not on disk here).
  - Align that library to any genome with our own code, resolve overlaps as DupMasker does, and validate against
    DupMasker's CHM13 track.
  - It should match closely by construction. It gives duplicon ids consistent across human and gorilla, so gorilla gets
    an expansion layer. The external piece becomes a data file, not a tool.
- **Route A: fully de novo duplicons.** Cut sd-scan's genome-wide alignments at transitively projected breakpoints and
  group the atoms. It needs no DupMasker data, but the closest attempt (PAF units, register 1263) reached 0.64 against
  DupMasker and today's family rule 0.52, so its ceiling is uncertain.

**Suggested order:**
1. Steps 0-1 (about a day): they tell whether recall is mostly an estimate problem.
2. Steps 2-4 if detection is the gap.
3. Route B in parallel, since it does not depend on sd-scan's recall.


---

## OPEN_IMPLEMENTATIONS_PHASED_2026-10-12.md

*(original: `docs/archive/2026-10/merged/OPEN_IMPLEMENTATIONS_PHASED_2026-10-12.md`)*

## Open implementations — phased data (2026-10-12)

**Why this exists.** A collaborator suggested working with phased data. The toy
(`docs/PHASED_TOY_2026-10-12.md`, prereg c5431f41, P0-P3 PASS) shows the principle in the EASY case: on a collapsed
reference a heterozygous allele and an extra copy give identical evidence; with two haplotype assemblies a truth-free
rule separates them 100% (dispersed duplication, unique flanks, K = 2, perfect assemblies, 0.43% divergence).
Nothing below is built. **Build an item only once its trigger is met.**

### Step 0 — ask before building (blocks everything)

> "When you say phased data, which form, and for which individual: (a) haplotype-resolved assemblies (two FASTAs,
> e.g. hifiasm trio/Hi-C), (b) a phased VCF (GT `0|1`, `PS` blocks) against which reference, or (c) a BAM with
> `HP`/`PS` haplotags? Is it the same animal as the IsoSeq (KB3781) or another one?"

| form | O2 (which copy / haplotype a read is from) | O3 (copy absent from one haplotype) | item |
|---|---|---|---|
| (a) two assemblies | yes | **yes** — the only form the toy validated | I1 |
| (b) phased VCF, collapsed reference | only in well-resolved regions | **no** — in a collapsed SD the caller reports PSVs as het sites and the phaser phases them: the toy's collapsed arm | I3 → I2 |
| (c) haplotagged BAM | as (b) — tags inherit the VCF's errors in SDs | no | I2 |

If the collaborator can provide (a), push for it; it is the form that makes O3 testable.

### I1 — two-assembly mode (`--hap-ref hap1.fa,hap2.fa`) · size M · trigger: form (a)

- **What:** promote the toy rule to a tool. Per read: map to each haplotype; NM1 == NM2 → SHARED; else home =
  lower-NM haplotype; ALLELE if the home locus has a counterpart between the same flanks on the other haplotype, else
  EXTRA_COPY (O3 candidate on the home haplotype).
- **Replaces:** the ad hoc `O3_REF=mat|pat` path in `bench/o3_maternal/` (`truth.py` already computes "locus on OTHER
  with no counterpart on REF" from lifts) and `bench/rna_allele/` chrmap logic. Reuse those, do not fork them.
- **Open design points:**
  - Counterpart test: the toy used a fixed ±2 kb flank, which works only because the flank clears the 500-bp
    duplicated pad. Real SDs are tens of kb → use the lift/synteny already in `bench/o3_maternal/truth.py`, not a
    fixed flank.
  - Tandem duplications: "which copy is the extra one" is undefined; report the copy-count difference per array, not
    a per-copy label.
  - Where it lives: new `copy_assign` mode vs `missing_copy_flag` input vs a `utilities` subcommand — decide with the
    O2/O3 owners of those binaries.
- **Tests first:** the toy's S0/SA/SB as fixtures (`bench/phased_toy/phased_toy.py` already generates them, 1 s);
  then a KB3781 family with a known mat/pat difference (GWFAM175_B0, the one copy beyond allele divergence).

### I2 — read haplotype from BAM tags (`copy_assign --hap-tag HP`) · size S · trigger: form (b) or (c)

- **What:** read `HP` (and `PS`) aux tags from the input BAM; carry them as a label into the assignments output
  (read → copy, haplotype, phase set); report per-copy haplotype counts.
- **Do NOT** feed the tag into the assignment decision until validated: in SDs the tag inherits the phaser's
  collapsed-duplication errors (form (b) row above).
- **Validation before use:** on KB3781, tags derived from its own phased variants vs the O2 certificate's mat/pat
  assignment (`bench/asj_o2/`, T1 discordance .0014 held-out) — agreement inside vs outside SDs.

### I3 — phased VCF input · size XS · trigger: form (b)

- **No VCF parser.** Convert with `whatshap haplotag --reference ref.fa phased.vcf.gz reads.bam` → I2. Document the
  command (and the RNA caveat: spliced reads are fine, but allele-specific expression leaves some haplotypes unseen).

### I4 — naming collision with the existing `copy_assign --phase` · size XS · do anytime

- `--phase` WRITES `PS`/`HP` where **HP = the assigned COPY**, not a parental haplotype (help text at
  `src/copy_assign_bin.rs:1076-1087` already says so; register line 640). An input flag about parental haplotypes next
  to it will confuse readers and downstream tools that consume `HP`.
- **Options:** rename the output flag (e.g. `--emit-copy-tags`, keeping `--phase` as a deprecated alias) or name the
  new input flags without "phase". Pick one before I1/I2 land.

### I5 — toy rungs (harder cases) · size S each · trigger: before claiming anything beyond the easy case

Each a pre-registered rung of `bench/phased_toy/` (or a famsim two-haplotype mode, `bench/famsim/`):
1. tandem duplication (no unique flank between copies);
2. gene conversion between haplotypes / between copies;
3. divergence sweep down to the read error budget (PSVs/het sites < errors per read);
4. K = 3 copies — identifiability by separation + linkage fails at K >= 3 (register line 561);
5. imperfect assemblies (collapsed copy in one haplotype assembly);
6. real KB3781 mat/pat sequences of one family instead of simulated variants.

### Limits that travel with any claim

- Copy vs allele cannot be settled on RNA alone (register line 767) — phased **assemblies** are the referee, not the reads.
- Genome-wide whole-gene CN difference between KB3781's haplotypes was uninformative (control 5.4x the signal,
  register line 751) — expect few real mat/pat copy differences; GWFAM175_B0 is the known one.
- Phased data must be of the SAME individual as the RNA for the O3 test; another individual's haplotypes test
  population CN variation, not reference absence in this sample.


---

## HANDOFF_laptop_to_main_2026-10-15.md

*(original: `docs/archive/2026-10/merged/HANDOFF_laptop_to_main_2026-10-15.md`)*

## Handoff — laptop → main machine (2026-10-15)

Work done on the second machine (WSL laptop, 6 cores / 9 GB, no `/mnt/linuxdisk`) lives on the branch
**`laptop/o3-psv-multimappers`**, on top of `main` = `49e62be3`. It is NOT merged. `main` was not touched from the laptop.

### What is on the branch

| commit | what |
|---|---|
| `797f55f8` | `pixi.toml` / `pixi.lock` (minimap2 2.30, samtools 1.22.1, bench and IsoCon envs); `.gitignore` exceptions for the two fixture sets that were never committed (below) |
| `5022350c` | O3 positive control re-measured (`docs/PREREG_missing_copy_posctl_2026-10-14.md`, PASS); few-read false flags diagnosed (`docs/PREREG_o3_fewread_flags_2026-10-14.md`: F = 0 on simulation); exon-structure projection in `candidate_copies` (`structure.tsv`, `contigs.gtf`; `docs/PREREG_o3_structure_projection_2026-10-14.md`, PASS) |
| `efc84916` | O3 per the advisor's definition (cluster the reads O2 places on a family by their PSVs; a cluster no copy explains = a copy absent from the genome). Prototype `bench/o3_unexplained/psvclust.py`; `docs/PREREG_o3_unexplained_2026-10-15.md`, PASS |
| `fdebaec0` | Rust port: **`copy_assign --psv-clusters`** (`rustle::missing_copy::psv_clusters`, OPT-IN), v2 with the linkage test against RNA editing and shared errors; `docs/PREREG_o3_psv_clusters_v2_2026-10-15.md`, all bars PASS (425/426 absent copies, 0/152 controls flagged, 16/16 runs in parity with the prototype) |
| this file | the steps below |

Design notes: `docs/superpowers/specs/2026-10-14-o3-structure-projection-design.md`,
`docs/superpowers/specs/2026-10-15-o3-unexplained-multimappers-design.md`. `REPRODUCE.md`, `bench/README.md` and the O3 row of
`docs/THESIS_OBJECTIVES.md` are updated on the branch.

### Steps on the main machine, in order

Use the canonical repo's usual build target (`CARGO_TARGET_DIR=...` on `/mnt/linuxdisk`) and the usual locks.

#### 1. Fetch the branch

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle        # the canonical clone
git fetch github laptop/o3-psv-multimappers   # or `origin`, whichever remote points at juanfraitu1/Rustle
git switch -c laptop/o3-psv-multimappers --track github/laptop/o3-psv-multimappers
```

If `main` moved on the main machine since `49e62be3`, rebase the branch on it (or merge `main` into it) before step 3.

#### 2. Commit the test fixtures that exist only on the main machine

`.gitignore` ignores `*.gtf` and `*.fa`, so these were never committed. A fresh clone fails 10 tests without them. The branch
already has the `.gitignore` exceptions (`!/tests/fixtures/bridge_units/*.gtf`, `!/tests/fixtures/o2_origin_resolve/*.fa`).

```bash
git add tests/fixtures/bridge_units/plain.gtf \
        tests/fixtures/bridge_units/expected.units.gtf \
        tests/fixtures/bridge_units/s_f0.5.plain.gtf \
        tests/fixtures/bridge_units/s_f0.5.expected.units.gtf \
        tests/fixtures/o2_origin_resolve/copies.fa \
        tests/fixtures/o2_origin_resolve/genome.fa
git status --short tests/fixtures      # anything else the tests read that is still untracked? add it too
```

Check sizes first (`du -sh`); all should be small synthetic files.

#### 3. Run the full suite

```bash
cargo build --release --bins
cargo test --release 2>&1 | tee /tmp/rustle_test_laptop_branch.log
```

Expected on the laptop: 1,299 passed, 13 failed. The 13 are known:
- 10 are the missing fixtures (step 2). They should pass here.
- 3 are stale tests in `tests/pipeline_driver_dispatch.rs`. The driver's `all)` branch now calls `run_all`, which the test does not
  stub, so it exits 127. Fix the test (stub `run_all`) or record it as known.

Anything else failing is new and must be looked at before merging. `copy_assign` without `--psv-clusters` is byte-identical to
before on the simulation runs, but that was not checked on real data.

#### 4. Merge to `main` and push

Only after step 3 is clean (or only the 3 dispatch tests remain and are recorded):

```bash
git switch main && git merge --no-ff laptop/o3-psv-multimappers
git push github main
```

Add the laptop results to `docs/NEGATIVE_RESULTS_REGISTER.md` (the posctl, few-read, structure-projection, v1 and v2 preregs above).
They were not given row numbers on the laptop to avoid clashing with rows added on the main machine.

#### 5. M3 foreign-species check (pre-registered, not run)

`docs/PREREG_o3_candidates_foreign_check_2026-10-14.md`. Inputs are already on `/mnt/linuxdisk/tmp/rna_allele/{a14a15,a15h}/`.
It asks whether the remaining `candidate_copies` false flags best match another ape. Treat any result as evidence about the reads,
not about how they entered the library. That is not established.

#### 6. O3 on real data: pre-register first, then run

The simulation level is done. The next level, from the v2 prereg's reading, is:

- Pick real multi-copy families from the gorilla (KB3781) catalog. **Remove copies from the GENOME** (not just the catalog, which
  is what the 2026-09 excision experiments did, so reads stayed on the removed locus), re-map the real reads with the exact O2
  settings, run `copy_assign --families ... --copies-fa ... --no-as-tied-only --flag-missing-copies --psv-clusters`. The removed
  copy's reads should form an `unexplained` cluster on its nearest paralog. Then add the copy back: the cluster must disappear.
- **Single-copy-gene control**: run the same on single-copy genes. Any `unexplained` cluster labelled `copy` there is the false
  rate the real heterozygous allele background gives. Simulation cannot measure this.
- Write the prereg (bars, families, held-out split) and commit it before running.

Run with `--no-as-tied-only`. With exact AS ties only (the shipped gate), O3 sees only identical-twin cases. The absorbed reads of
a missing copy are mostly unique mappers.

#### 7. Smaller open items

- `o3_export.py` does not yet read `candidate_copies`' `contigs.gtf` (secondary under the advisor's framing).
- HG002 plan: `docs/PENDING_2026-10-11.md`.

### Laptop notes (if work continues there)

- pixi envs: `pixi run <cmd>` (default), `pixi run -e bench`, `pixi run -e isocon`. Build with `-j 4`
  (`.cargo/config.toml`'s `jobs = 48` is for the main machine).
- Data: `~/rustle_data/` (CHM13 chr20 + RefSeq annotation, the simulation worlds). Real BAMs will come from the NAS.
- The laptop keeps committing to `laptop/*` branches, never to `main`.

### Done on the main machine (2026-10-15)

- **Steps 1-4.** `main` (`4b69ad74`, 64 commits past `49e62be3`) was merged into the branch (`integrate/o3-psv-multimappers`,
  merge `90df29df`) with no conflicts.
  - Committed the fixtures: the six files above, plus `o2_origin_resolve/reads.bam`, `reads.bam.bai` and `genome.fa.fai`.
    These three are also needed (the `copy_assign_popfix` and `o2_origin_resolve` tests read the BAM), and the BAM's records
    equal the tracked `reads.sam`.
  - Fixed the 3 dispatch tests: they now run the driver's real `run_all` with a stub per-stage driver, so the `all` plan and
    its order are checked again. A broken plan fails them.
- **Suite on the merged tree:** `--lib` 1,099 passed; `--bins --tests` 1,356 passed, 0 failed, 14 ignored; no doctests.
- **Real data, without `--psv-clusters`.** `copy_assign` from the merged tree and from `main` gave byte-identical outputs on human
  testis (40 families; default gate and `--no-as-tied-only`: 5,240 assigned reads) and gorilla KB3781 (40 families default,
  20 with `--no-as-tied-only`: 27,525 assigned reads).
- **Register rows 1289-1293** record the posctl, few-read, structure-projection, v1 and v2 preregs.
- Still open: step 5 (M3 foreign-species check), step 6 (the real-data O3 prereg), step 7.


---

## O1O2O3_FAST_REALISTIC_2026-10-16.md

*(original: `docs/archive/2026-10/merged/O1O2O3_FAST_REALISTIC_2026-10-16.md`)*

## Faster family tests, realistic reads, and what O2's evidence rests on (2026-10-16)

> Update 2026-10-19. **O2's domain (user, 2026-10-19): AS-tied multimapping reads only** — a read with a rival alignment within 10 AS of its best (Eichler's within-10 rule; under splice:hq, 10 AS is two substitutions). A read with a clear best placement is already placed by minimap2: applying O2 to it is tautological. O2's point is to go beyond minimap2's shortcoming, its coin toss among AS-tied alignments, and assign those ambiguous reads to very similar copies by their PSVs, junctions and divergence, or abstain. Implementation: `fastset --o2-scope tied --as-slack 9` (`o2graph.py run W intega`), on test in round 21 (`docs/PREREG_o3_graph_round21_2026-10-19.md`); memory `feedback_o2_scope_ambiguous_only`.

**User goal (2026-10-16):** "integrate O1, O2 and O3 more so it's faster and cheaper to test full families; also improve the
simulations to account for more realistic ones and determine if PSVs can also Include SUNs and probably some influence of
the cigar string, nm, de tags".

Everything below is development work on family sets 8 and 9 (rounds 7 and 8), both already seen. No claim is made here.
The pre-registered test on a fresh set is round 9 (`docs/PREREG_o3_graph_round9_2026-10-16.md`): FAIL under both
profiles. Section 9 holds the independent review's corrections; where they disagree with a section above, section 9 wins.

### 1. Faster, cheaper family tests

**Where the time went (one round-8 world, 8 families, 2,537 reads):**

| step | wall | peak memory |
|---|---|---|
| simulate | 0.6 s | 38 MB |
| map to the whole genome (index load) | 79 s | 14.5 GB |
| current method (copy_assign) | 16 s | 1.06 GB |
| O1 (o1seed) | 0.6 s | 31 MB |
| O2 (integ) | 5 s | 51 MB |
| O3 (detect) | 17 s | 98 MB |

REMOVED also needs a masked whole-genome index: 146 s and 19 GB.

**What changed:**
- **`localref.py`: a homology-closed local reference.** It holds every locus any isoform of the set aligns to (minimap2
  splice:hq, `-N 500 -p 0.01`, more permissive than the reads' `-N 50 -p 0.1`), plus each copy's span, padded 10 kb and merged.
  - Set 9's is 4.6 Mb in 75 regions, against the 3.1 Gb genome; its index is 22 MB.
  - Reads are mapped with the reads' exact minimap2 command and lifted back to genome coordinates (SA tags included), so
    every later step is unchanged.
  - It is built once per family set (38 s, the only whole-genome step). Its masked REMOVED variant takes 0.5 s.
- **Validation against whole-genome mapping (all six round-8 worlds):**
  - Unique placements (MAPQ >= 1) get 100% the same primary.
  - 91-93% of reads have the identical full alignment set.
  - Every MAPQ-0 primary that moved landed on an equal-score alignment that whole-genome mapping also reports (one exception
    per world). Those are coin flips, decided arbitrarily by either index.
- **`fastset.py`: one command per family set.**
  - It simulates every world.
  - It maps each distinct read sequence once (INTACT and OUTSIDE share their reads, CROSS and BRANCH contain INTACT's, SILENT
    is a subset): 3,284 distinct sequences for 12,712 world reads.
  - It runs the current method and O1 -> O2 -> O3 per world, two worlds at a time under the light lock, then scores the end
    to end and the bars, timing every step.
- **Round 8 rerun with `fastset.py`:**
  - 237 s at 1.08 GB peak, against about 655 s at 19.4 GB for round 8's chain, plus the one-time 38-s local reference.
  - Same verdict. Five worlds have identical bars; SILENT differs by 36 correct reads, a MAPQ-0 coin flip.

### 2. Realistic reads

**`calibrate.py`** measures real Iso-Seq alignments made with the same minimap2 command on the same CHM13 v2.0:
- **Errors and donor alleles:** primary MAPQ-60 reads inside single-isoform genes, so every difference is an error or the
  donor's allele.
  - A mismatch shared by >= 3 reads at >= 20% of coverage is the donor's allele (het below 80%, hom above).
  - Indels are classed by homopolymer context.
- **Ends:** 5' and 3' truncation in transcript bases.
- **Retained introns:** one aligned block spanning a whole annotated intron.
- **Expression:** the spread of per-copy primary counts within families whose copies can be told apart.

| per aligned base | A119b (lab testis; development) | public testis (ERR13885926; held back) | old simulator |
|---|---|---|---|
| substitution errors | 0.028% | 0.092% | 0.1% |
| insertion events | 0.087% | 0.046% | 0.03% |
| deletion events | 0.125% | 0.033% | 0.03% |
| indels in homopolymers | 64-65% | 60-64% | none by design |
| donor het / hom alleles | 0.027% / 0.033% | 0.014% / 0.036% | het 0.05%, hom none |
| reads reaching the 5' end | 41% | 56% | all (0-30 bp jitter) |
| retained intron (per intron) | 0.89% | 0.14% | none |
| per-copy expression log-sd | 0.86 | 1.16 | 0 (60 reads per copy) |

**`realworld.py --sim-profile`** draws:
- each read's 5'/3' truncation from the calibrated (5', 3') pairs;
- errors at the calibrated rates, with homopolymer-biased indels and calibrated lengths;
- the donor's hom alleles on both haplotypes and het alleles on one;
- retained introns;
- a lognormal per-copy expression.

Without a profile it is byte-identical to before (round 8's INTACT and CROSS reads reproduced exactly).

**Check against real reads.** Uniquely placed reads (MAPQ 60), `calibrate.py tags`:

| | de (median) | NM/kb | mismatches/kb | ins/kb | del/kb | indels in homopolymers |
|---|---|---|---|---|---|---|
| real A119b | 0.0024 | 2.65 | 1.42 | 0.89 | 1.22 | 64% |
| real public testis | 0.0017 | 1.78 | 1.39 | 0.47 | 0.38 | 61% |
| old simulator | 0.0019 | 1.93 | 1.67 | 0.30 | 0.34 | 18% |
| A119b profile | 0.0028 | 2.87 | 0.90 | 0.82 | 1.26 | 70% |

The two real libraries differ (indels 2-3x apart), so a profile is a library, not "real reads" in general. Results are
reported under both profiles.

### 3. What realistic reads change (integrated pipeline v1 vs the current method, `fastset.py`)

| set, profile | bars | where it fails |
|---|---|---|
| 9, old simulator | FAIL (round 8) | CROSS errors 1 vs 0 |
| 9, A119b | FAIL | B1: INTACT 926 vs 932, CROSS 1,011 vs 1,030; B2: REMOVED 8 vs 0 |
| 9, public | FAIL | B2: REMOVED 7 vs 0 |
| 8, old simulator | PASS | |
| 8, A119b | FAIL | B1: INTACT -93, REMOVED -70, SILENT -19; B2: CROSS 13 vs 8, REMOVED 17 vs 3, INTACT and SILENT 1 vs 0 |
| 8, public | FAIL | B2: CROSS 10 vs 8, REMOVED 8 vs 4 |

The integrated pipeline's advantage on catalogued copies was partly an artifact of idealized reads. It keeps its advantage
where the catalog is incomplete: OUTSIDE errors 24 vs 154 and BRANCH 1,194 vs 932 correct (set 9, A119b).

**Why:**
1. **O2 calls more reads `unexplained`:** 197 of 2,537 under A119b, against 8 with idealized reads. Of 139 examined:
   - 38 align to no path;
   - 51 have NM above delta x read length, because a homopolymer indel of 4+ bases or a retained intron counts per base;
   - 21 are clipped >= 50 bases;
   - 29 are short truncated reads where a few edits exceed delta x read length.
2. **REMOVED's catalogued errors are present-copy reads inside true absent-copy calls:** 7 of 7 on set 9 and 14 of 15 on
   set 8 under A119b. A119b's noisier reads widen the clusters.

### 4. Do PSVs include SUNs?

**Yes: a SUN is the PSV column where one copy's base differs and every other copy agrees.**
- Of set 9's 284 PSV columns, 196 are SUNs (69%), 77 shared and 11 two-copy pairs. The de novo families of 10-12 gave 78%.

From the current O2's read-star proofs (`copy_assign --dump-star`, which changes no call), each decisive column is classed
among the read's candidates:
- **About 70% of decisive columns are SUN-type:** the called copy's own SUN, or another copy's SUN the read rules out.
- **Complete catalogs (INTACT, CROSS): no call is wrong** under any profile, with or without an own-SUN column (0 of 827 to
  1,226 per world).
- **Missing copies are where every error is.**
  - OUTSIDE under idealized reads: 60 of 121 wrong calls rest on the called copy's own SUN. The uncatalogued copy carries that
    "private" base too.
  - **A SUN is private only relative to the copies considered.** Genome-wide SUNs (unique k-mers, Sudmant 2010) would guard
    against copies in the genome but missing from the catalog. Nothing in the reference guards against copies absent from it.

### 5. CIGAR, NM and de

- **de flags absent-copy reads.**
  - In REMOVED, the current O2's wrong calls sit in the reads' top de quartile: 85 of 254 against 26 of 354 in the bottom
    quartile (idealized), and 39/191 against 12/202 under A119b.
  - Within a complete catalog it predicts nothing, because there are no errors.
- **Soft clips (CIGAR S).**
  - O2's distance counts every clipped base as an edit. That made round 8's error: one sequencing error among 4 clipped
    bases scored 1 vs 5.
  - Aligning clipped ends against the path (semi-global) removes that error. But it also removes the inflated margins of a
    few correct single-difference calls: about -3 correct per world on set 9, +3 to +8 on set 8.
- **NM vs de: count indel events, not indel bases.**
  - Real errors are mostly homopolymer indels, and a retained intron is one event of many bases.
  - A gap-compressed distance (mismatches + indel events, like minimap2's `de`, with semi-global ends) recovers O2 on
    realistic reads:
    - set 9, A119b: INTACT 926 -> 953 correct, unexplained 197 -> 106;
    - set 8, A119b: INTACT 2,141 -> 2,215 correct, wrong 1 -> 0; CROSS 2,628 -> 2,719 correct, wrong 9 -> 1.
  - Under idealized reads it costs some calls on set 8 (INTACT 2,476 -> 2,407), with fewer errors (CROSS wrong 11 -> 3).
- **The read's own de cleans absent-copy calls.**
  - Present-copy reads inside true absent-copy calls sit within delta of their own locus: 51 of 55 at de <= delta. True
    absent-copy members do not: 45 of 1,154 at de <= delta.
  - A per-read genome check (`detect.py --member-genome-check`, arm I) drops them.

### 6. Integrated v2 (opt-in; defaults byte-identical)

**What it runs:** `o2graph.py run W integ 2 --distance gapcomp` + `detect.py W --member-genome-check`, or in one command
`fastset.py --o2-distance gapcomp --member-genome-check`.

**Development results** (bars per world, as round 8; the per-read column also classes a REMOVED copy's read as detectable by
its own `de`, `summary.py e2e --per-read-detectability`):

| set, profile | v1 | v2, copy-level | v2, per-read |
|---|---|---|---|
| 9, idealized | FAIL: CROSS errors 1 vs 0 | PASS | |
| 9, A119b | FAIL: INTACT, CROSS, REMOVED | FAIL: REMOVED errors 1 vs 0 | PASS |
| 9, public | FAIL: REMOVED errors 7 vs 0 | FAIL: REMOVED errors 1 vs 0 | PASS |
| 8, idealized | PASS | PASS | |
| 8, A119b | FAIL: 4 worlds | FAIL: INTACT 2,215 vs 2,234 correct; REMOVED errors 5 vs 3 | FAIL: INTACT 2,215 vs 2,234 correct |
| 8, public | FAIL: CROSS, REMOVED | PASS | PASS |

- **Set 9's one REMOVED error under both profiles** is a read of the removed USP17L.8 (copy median `de` 0.018) that is
  near-identical to USP17L.9 over its own truncated span. Per-read detectability classes it unknowable.
- **Set 8 under A119b reads, INTACT's shortfall (-19)** is O2's `unexplained` reads (SPRR2 37, AQP7 12), which the current
  method assigns.
  - A MAPQ-60 genome fallback for unexplained reads was measured and rejected.
  - Outside BRANCH it places 17 of 329 reads wrong, in worlds where the current method makes almost no errors.
  - In BRANCH it places 159 of 220 wrong: novel-exon reads land uniquely on their parent copy.
- **v2 against v1 under idealized reads:** fewer correct reads on set 8 (INTACT 2,407 vs 2,476), fewer errors (CROSS 3 vs 11,
  REMOVED 0 vs 13).

**Round 9** tests v2 on a fresh set: `docs/PREREG_o3_graph_round9_2026-10-16.md`.

### 7. Full families

Uncapped (`families.py build --cap 1000`; `/mnt/linuxdisk/tmp/o3_graph/families_full.json`):
- set 2: NPIPA 6, NPIPB 14, TBC1D3 10, GOLGA8 18 and four others; 69 copies;
- set 9: USP17L 59, TSPY 35 and six others; 116 copies.

Runs with `fastset.py` (integrated v2, A119b profile; post hoc, both sets already seen):

| set | copies | wall time | peak memory | heaviest steps (summed over worlds) | bars |
|---|---|---|---|---|---|
| 2 (NPIP, TBC1D3, GOLGA8, ...) | 69 | 1,508 s | 14.4 GB (local reference build only) | O3 1,215 s, current 457 s, O2 434 s | PASS |
| 9 (USP17L 59, TSPY 35, ...) | 116 | 1,203 s | 14.3 GB (local reference build only) | current 579 s, O2 543 s, O3 494 s | PASS |

- **Full families are feasible but no longer cheap.**
  - O3's `vg msga` family graphs grow steeply with copies and locus length.
  - copy_assign peaks at about 5 GB per world on set 9.
- **Since then:**
  - `fastset.py` shares family graphs between worlds with one catalog (`O3_GRAPH_CACHE`). INTACT, CROSS, BRANCH and SILENT
    share the catalog, so their graph is built once.
  - `--heavy-current` runs copy_assign under the heavy lock.
- These sets used the old family rule (section 8).

### 8. Domain sharers in the simulation's families (user, 2026-10-16)

The question: do the families hold genes that share only a domain with the others (over-representation)? Would requiring
links over >= 80-90% of the shorter transcript remove them?

**Families today** (`families.py`):
- genes are taken by name pattern;
- two are linked when their canonical transcripts align at >= 0.90 identity over >= half the shorter one;
- a family is the largest single-linkage component.

**Yardstick (`domains.py`).** Each gene is graded against its family by spliced alignment of its isoforms to the other genes'
genomic LOCI, in both directions:
- `whole` (each locus holds >= 80% of the other's transcript);
- `fragment` (one direction only);
- `partial` (50-80%);
- `segment` (< 50% either way: a domain sharer).

It uses the partner's genome, not its annotation. An earlier locus-to-locus version failed on old duplicates whose introns
diverged (SSX, AMY scored 0) and was replaced.

**Results over all ten sets' candidate families, uncapped (105 families with >= 2 genes):**

| link rule | genes | whole | fragment | partial | segment | genes dropped (of which whole + fragment) |
|---|---|---|---|---|---|---|
| canonical transcripts, >= 0.5 of the shorter (today) | 644 | 591 | 41 | 10 | 2 | |
| canonical transcripts, >= 0.8 | 588 | 545 | 33 | 4 | 0 | 56 (41) |
| canonical transcripts, >= 0.9 | 548 | 506 | 26 | 4 | 1 | 96 (73) |
| **locus: canonical transcript >= 0.8 inside the other gene's locus, one direction** | 670 | 617 | 53 | **0** | **0** | |

- **Domain sharers were rare.**
  - Only 2 of 644 genes (KIR2DL4, KRTAP9-7) were segments.
  - 10 were partial: NBPF7, HERC2P7, GUSBP3, SPANXN1/5, KRTAP5-4/-10, KRTAP4-2/-12/-16.
  - Over-representation was 1.9% of genes, not more.
- **Transcript links over >= 80% do remove them, at about five genuine copies per non-copy removed.**
  - The losses are copies whose canonical transcripts carry different exon sets; NPIPB's 14 whole copies fall to 7.
- **The locus rule removes all 12 and adds 38 copies the transcript rule missed:**
  - pseudogene copies such as LRRC37A and GSTA7P;
  - diverged-annotation copies such as NPIPB 14 -> 17, SSX 5 -> 8, PSG 7 -> 9.
- **Linking on any isoform instead of the canonical one inflates families** through short single-domain isoforms (LILR 2 -> 9,
  MT1 2 -> 9). The rule uses the canonical transcript for that reason.
- **The rebuilt simulation sets** (`/mnt/linuxdisk/tmp/o3_graph/families_locus.json`, each set's original cap and selection)
  change in every set.
  - Membership changes, for example GOLGA6L 4 -> 7, SPATA31 7 -> 5, MAGEA 5 -> 3.
  - Four first-8 selections shift: set 5 VCX -> CTAG, set 6 ANKRD20A -> CYP2D, set 8 FAM138 -> PRSS, set 10 KRTAP4 -> KRTAP9.

**Runs on the rebuilt sets** (integrated v2, v1 profiles; post hoc; per-read scoring):

| set | old families | rebuilt families |
|---|---|---|
| 10, A119b | FAIL: B1 INTACT -40, CROSS -37, REMOVED -25, SILENT -12 | FAIL: B1 INTACT -54, CROSS -51, SILENT -20; REMOVED +20 |
| 10, public | FAIL: B1 INTACT -18, CROSS -14, SILENT -10; B2 INTACT 5 vs 1, CROSS 5 vs 2, SILENT 2 vs 1 | FAIL: B1 INTACT -43, CROSS -38, SILENT -31; B2 INTACT 5 vs 1, CROSS 5 vs 2, SILENT 2 vs 1 (OUTSIDE errors 8 vs 547) |
| 9, A119b | PASS | FAIL: B2 CROSS 3 vs 0 |
| 8, A119b | FAIL: B1 INTACT -19 | FAIL: B1 INTACT -41, REMOVED -74, SILENT -50 |

- **Rebuilding the families changes no conclusion.** In every rebuilt set the integrated pipeline still has fewer errors in
  every world and more correct reads wherever the catalog is incomplete (OUTSIDE, BRANCH), and still loses correct reads on
  catalogued copies.
- **Single-read bars move.** For example, rebuilt set 9's 3 CROSS errors are reads of the newly added OR4F17 that carry OR4F5's
  annotated structure; they are near-identical to OR4F5.

### 9. Independent review of this branch and round 9 (2026-10-16)

No critical finding: the round-9 FAIL reproduces exactly. The corrections:

- **Public-profile errors (round 9):** they come from O1 dropping expressed copies with 4-5 reads, under its fixed 5-read
  support rule (PGA.2 had 2 uniquely placed reads). They were not an O2 effect. Recorded in the round-9 document and register
  row 1307.
- **Cost:** the per-world memory figures (about 1 GB) leave out the one-time local reference build: 38-75 s at about 14 GB, once
  per family set. A new family set still peaks near 14 GB (19.4 GB before); every later run of that set stays under 2 GB per
  world.
- **Read profiles, v2 (`profiles/*_v2.json`; v1 kept, and reproduces its reads byte for byte):**
  - **errors and donor alleles from every MAPQ-60 read** (`calibrate.py errors`). Single-isoform genes had under-counted the
    donor's het alleles:
    - A119b: substitution errors 0.41, het 0.58, hom 0.39 per kb, 1.38 mismatches per kb in all, against 0.89 in v1;
    - the simulator gives 1.33 per kb against 1.42 in real reads (v1 gave 0.90);
  - **expression spread** over copies with reads, less Poisson counting noise: A119b 0.86, public 1.03 (v1 1.16);
  - **truncation** redrawn while it would leave < 200 bases. v1 fell back to full length for 17-61% of short transcripts' reads.
- **Not simulated: the public library's soft clips (12-14% of bases)** are the Tn5 mosaic-end adapter (`AGATGTGTATAAGAGACAG`),
  left untrimmed at both ends of 99% of reads.
  - That is a preprocessing gap, not read noise. Such reads must be adapter-trimmed before any O2 that compares a whole read
    with paths.
  - A119b is trimmed: 0.3% clipped.
- **Local vs whole-genome mapping on realistic reads** (rebuilt set 10, A119b INTACT):
  - 99.87% of unique placements are identical; the 2 that differ are the same copy with a different first exon at the read's
    5' end;
  - 96.3% of reads have identical alignment sets;
  - every moved MAPQ-0 primary moved to an equal-score alignment.
- **Wording:**
  - the public library is called "held back" above, but it shaped v2's development;
  - its BAM was mapped without `-uf`;
  - the real-reads table mixes medians (de, NM) with pooled rates;
  - "complete catalogs: no wrong call" holds for sets 8-9, not round 9 (the current method misplaces FCGR1.0 reads);
  - a read's MAPQ-0 tie choice also depends on its name, so mapping by sequence moves tie choices (bars unchanged).
- **Numbers that came from scratch scripts:**
  - "51 of 55 / 45 of 1,154" (per-read genome check) was counted on the v1 REMOVED worlds: set 9 A119b, public and idealized,
    and set 8 A119b. The reviewer's 46/47 and 45/786 come from other worlds.
  - The MAPQ-60 fallback numbers (17/329, 159/220) were counted on the v2 development worlds; no output was saved.
- **Code fixes, each tested:**
  - `realworld.py` loads `o3_fewread/world.py` by path; a same-named `world.py` shadowed it in mixed test runs;
  - `fastset.py` refuses a local reference built for another family set;
  - `localref.py` lifts a same-contig mate.

**Post hoc runs added after the review** (integrated v2, per-read scoring):

| run | catalogued copies, correct vs current | errors | verdict |
|---|---|---|---|
| set 10, idealized reads (review item 6) | INTACT +8, CROSS +22, SILENT -1 | REMOVED 3 vs 1 | FAIL by single reads |
| set 10, A119b v1 (round 9) | INTACT -40, CROSS -37, SILENT -12 | fewer everywhere | FAIL |
| rebuilt set 10, A119b v2 | INTACT -63, CROSS -61, SILENT -22 | fewer everywhere | FAIL |
| rebuilt set 10, public v2 | INTACT -41, CROSS -41, SILENT -27 | INTACT 3 vs 2, SILENT 3 vs 2 | FAIL |

- **The effect of realistic reads is now measured on the held-out set itself.** Idealized reads give the integrated pipeline a
  small gain on catalogued copies; both calibrated libraries turn it into a 1-4% loss.
- **The corrected expression spread does what the review predicted.** Under public v2, O1 drops 0-2 expressed copies per world,
  against 2-4 under v1.
- **The pattern does not change.** In every realistic run, the integrated pipeline has far fewer errors wherever the catalog is
  incomplete (OUTSIDE: 3 vs 407-547) and loses correct reads on catalogued copies.

**Open:**
- O1's fixed 5-read rule for keeping a copy loses low-expressed copies under realistic expression. A rule on uniquely placed
  reads would keep them; it is the next O1 candidate.
- O2's divergence gate on realistic reads remains the catalogued-copy weakness.

### 10. Closing the remaining gaps (after the goal check, 2026-10-16)

**O1 keeps a copy with one uniquely placed read of its own** (`o1seed.py`). The 5-read primary count stays, so an identical
copy whose reads all tie is still kept. Over 138 worlds:

| profile | expressed copies dropped: 5-read rule / new rule | silent copies kept: 5-read rule / new rule |
|---|---|---|
| public | 68 / 14 | 35 / 35 |
| A119b | 11 / 11 | 82 / 84 |
| idealized | 2 / 2 | 22 / 22 |

The 14 expressed copies still dropped have only MAPQ-0 (tied) reads, so no read can single them out. End-to-end effect, on the same reads (rebuilt set 10, public v2 profile, integrated v2):
- **Expressed copies dropped:** 1-2 per world -> 0 in every world. The 8 silent copies are still dropped.
- **Integrated errors:** 3 -> 0 in each of INTACT, CROSS and SILENT, so B2 (errors <= current) now holds in all six worlds.
- **Correct reads:** +2 to +4 per world.
- **B1 still fails** on catalogued copies (INTACT -39, CROSS -39, SILENT -23): O2's divergence gate, not O1.

**Full families: faster and cheaper.** All of the following are exact: outputs are identical with and without them.
- **Family graphs shared** between worlds with one catalog (`O3_GRAPH_CACHE`).
- **O3 arm I only** (`detect.py --arms I`): the legacy arms' alignment and clustering are skipped; arm I's calls are
  identical.
- **A world's families in parallel** (`O3_PROCS`): identical with 1 or 2 processes, and to the previous code.
- **O2 scoring pruned** (gap-compressed distance): byte-identical.
- **4 threads for steps that run alone:** the shared mapping, REMOVED's mapping, and copy_assign under the heavy lock.
  copy_assign's output is identical at 2 and 4 threads.

| uncapped set | before (no cache) | cache + arm I only | + parallel families, 4-thread lone steps |
|---|---|---|---|
| 2 (69 copies) | 1,461 s, O3 1,215 s | 1,316 s, O3 770 s, peak 1.5 GB | **967 s**, O3 529 s, shared mapping 105 s (was 204 s), peak 2.2 GB, PASS |
| 9 (116 copies) | 1,145 s, O3 494 s | 1,027 s, O3 214 s, peak 5.0 GB (copy_assign alone) | |

Times exclude the one-time local reference build. All runs PASS their bars post hoc (A119b v1 profile).
- **Set 2: 1.5x faster than the first full-family run, O3 2.3x faster.**
- **What remains:**
  - copy_assign, the comparator rather than the integrated pipeline: about 60 s per world;
  - O2's path alignment: minimap2 of every hosted read against every path, about 80 s per world;
  - arm I's pooled per-host clustering: about 90 s per world.
- **Peak memory** stays under 2.2 GB per set once the local reference exists; set 9's largest copy_assign reaches 5 GB, alone
  under the heavy lock.

**Untrimmed-library adapters are now simulated and trimmed:**
- `calibrate.py adapters` measures them: the public library's 5' end is Tn5 mosaic end + G run + SMART primer (48 bp), its 3'
  end polyA + `GCATCG` + mosaic end (40 bp), on 98.8% of reads.
- Profile v3 (`public_testis_v3.json`) attaches them at the calibrated clip lengths, with sequencing errors.
- `trim.py` / `fastset.py --trim-adapters` removes them before mapping.

**The check (rebuilt set 10, public v3, integrated v2, per-read scoring):**

| reads | integrated, INTACT: correct / errors | current, INTACT: correct / errors | integrated O2 calls `unexplained` |
|---|---|---|---|
| adapters left on | 0 / 1,874 | 1,145 / 0 | 2,094 of 2,094 reads |
| adapters trimmed (`--trim-adapters`) | 1,334 / 0 | 1,368 / 2 | as without adapters |
| no adapters (v2, same O1 rule) | 1,327 / 0 | 1,366 / 2 | |

- **Left on, the integrated pipeline collapses.**
  - O2 calls every read unexplained.
  - O3 makes 8 false absent-copy calls in INTACT, a world with no absent copy, and sweeps the present copies' reads into them.
  - The current method aligns to the genome and clips the adapters, so it is barely affected.
- **Trimming restores it to the adapter-free level.** Every read is trimmed, removing 85.5 bases on average, and no adapter
  core is left at either end.
- **Rule:** an untrimmed library needs `--trim-adapters` (or lima/refine) before the integrated pipeline.

### 11. Integrated v3: the catalogued-read loss and the false calls (after the second goal check, 2026-10-16)

Development sets 8, 9 and 10 (all spent), simulated fresh with the v2 profiles and the new O1 rule. Every run is scored per
read (`--per-read-detectability`); "dev runs" are the 6 set x profile runs.
- Runs: `/mnt/linuxdisk/tmp/o3_graph/fast/fbdev/` and `fbdev_exact/` (`--union-exact`).
- Error lists: scratch scripts `errs.py`, `vs.py`, `why.py`, `crosserr.py` (session scratchpad), summarised here.

**O2: the genome fallback (`integl`).** It applies to the reads O2's paths leave `unexplained`.
- **Scoring:** every primary and secondary genome alignment, as mismatches + indel events.
- **What refuses:** an insertion or clip >= 30 bp (that is O3's evidence), and the read's best alignment lying off the
  family's seeded copies.
- **Placement:** within delta, at margin 2.
- **It removes the round-9 loss.** Correct reads, fallback minus no fallback: +29 to +90 per world (set 8 A119b INTACT
  2,132 -> 2,196), with no new error outside BRANCH. In BRANCH it adds +5 to +16 errors, reads of near-identical absent
  copies; current's BRANCH errors are 46-225.
- **Three other modes were measured on the same runs and rejected:**
  - `combined`: the genome on every read, with the guard. Errors drop (OUTSIDE 76 -> 0), but correct reads collapse
    where arrays have identical uncatalogued units (set 9 INTACT 935 -> 534).
  - `seeded`: the genome on explained reads, without the guard. It removes CROSS's structure-only errors, but INTACT loses
    13-47 reads.
  - `resolve`: the genome on O2's tied and thin reads. It adds 0-1 read per world.

**O3: three false-call mechanisms, found read by read.**
1. **The mosaic rule's nearest copy took the host's base wherever a single read covered a position** (a position covered by
   one read is no variable site). The cluster then disagreed with its own copy:
   - in set 9 LCE1, a cluster of LCE1.1's reads (one carrying a donor allele) was 0.025 from LCE1.1, so it became a false
     mosaic;
   - with the reads' own bases (`cluster_bases`, f080a30f), the distance is 0.0008.
2. **Mixed clusters:** reads of several present copies clustered together, a false mosaic (set 8 FAM138, once the
   consensus is honest).
3. **Present reads swept into a true call:** 2 USP17L reads in set 9 REMOVED.
- **Member attribution (`--member-attribution`, 12f8f2af) fixes all three.** An absent-copy PSV call keeps only the reads its
  consensus explains better than every single catalogued copy, by O2's margin of 2.
  - O3 errors fall to 0 in INTACT, CROSS and SILENT;
  - REMOVED's correct reads change by -2 to +13 per run.

**O2: the CROSS errors are failed projections.** All 17 CROSS errors of set 8 (both profiles) are reads of two structures
that O2's union never gave to the true copy:
- **WASH T0075974's last exon is 6 bp.** The spliced aligner cannot place it, so the projection had 4 exons, not 5, and was
  rejected.
- **DDX11L T0075982 projected,** but was merged with a structure within 10 bp of the copy's own.
- **`--union-exact` (36a430f8) fixes both.** A projection may lack a terminal exon < 30 bp, and structures merge only when
  their junctions match exactly. CROSS errors: 14 -> 1 and 10 -> 1.
- The simulator's projection (liftover) and its defaults are unchanged, so the test stays independent of O2's projection.

**v3 (fallback + attribution + `--union-exact`) against the current method, correct / errors:**

| world | set 8 A119b | set 8 public | set 9 A119b | set 9 public | set 10 A119b | set 10 public |
|---|---|---|---|---|---|---|
| INTACT | 2,191/0 vs 2,167/0 | 2,250/0 vs 2,153/1 | 935/0 vs 907/0 | 873/0 vs 828/0 | **1,136/0 vs 1,138/2** | 1,057/0 vs 1,051/2 |
| CROSS | 2,708/1 vs 2,551/8 | 2,788/1 vs 2,558/9 | 1,024/0 vs 1,001/0 | 963/0 vs 917/0 | 1,256/0 vs 1,248/3 | 1,159/0 vs 1,146/3 |
| SILENT | 1,917/0 vs 1,860/0 | 2,085/0 vs 1,909/1 | 678/0 vs 644/0 | 628/0 vs 587/0 | 702/0 vs 689/2 | 613/0 vs 602/2 |
| REMOVED | 2,137/1 vs 2,059/12 | 2,158/0 vs 2,072/7 | 798/0 vs 782/0 | 713/0 vs 698/0 | 1,049/2 vs 1,031/10 | 1,021/2 vs 963/4 |
| OUTSIDE | 2,232/76 vs 1,979/150 | 2,273/62 vs 2,000/120 | 760/23 vs 674/152 | 681/19 vs 605/123 | 1,095/3 vs 863/169 | 1,028/3 vs 789/167 |
| BRANCH | 2,505/104 vs 2,203/137 | 2,509/166 vs 2,176/191 | 1,196/36 vs 907/60 | 1,080/30 vs 828/46 | 1,320/66 vs 1,138/225 | 1,192/81 vs 1,051/200 |
| verdict | PASS | PASS | PASS | PASS | FAIL (B1 INTACT, -2) | PASS |

- **The remaining shortfall is structural.** Set 10 A119b INTACT loses 15 short reads (184-484 bp) whose 2-5 sequencing edits
  exceed delta's share of their length.
  - A looser gate would place reads that per-read scoring counts as detectable in REMOVED.
  - The current method makes 0 errors in REMOVED in set 9.
- **Faster and cheaper, every change exact (outputs byte-identical):**
  - **O2 on a full set 2 world: 83 -> 49 s.**
    - The clipped-end edit distance uses Myers/Hyyro bit vectors; it equals the DP on 3,000+ random cases.
    - CIGAR events are parsed by regex.
    - Families run in parallel (`O2_PROCS`, `--o2-procs 2`).
  - **O3 arm I on full set 2 worlds: 55-59 -> 32 s.**
    - Isoform maps are built lazily.
    - PSV clustering keeps its consensus current instead of recomputing it per comparison, and computes linkage on read bit
      sets; it equals the reference implementation on 60 random hosts each.
  - **O3's graph:** building is cached (section 10). The remaining O3 time is GraphAligner on the hosted reads (about a
    third) and the per-host pool alignments.
- **Round 10 pre-registers v3** on an eleventh set: `docs/PREREG_o3_graph_round10_2026-10-16.md`.

### 12. After round 10: v4 (2026-10-16, evening)

Round 10 (`docs/PREREG_o3_graph_round10_2026-10-16.md`, register 1308) failed. Its independent review corrected the first
failure's mechanism. Development since then uses round 10's worlds, now spent, plus sets 8-10, all with the v2 profiles; that
is nine set x profile runs including round 10's untrimmed-library run.

**1. Attribution no longer keeps a read that no catalogued copy competes for** (88c445d9). Round 10's false TAF11L call came
from the host window:
- TAF11L.44's +-2 kb window holds 225 bp of TAF11L.43, and TAF11L.33 reads aligned to that fragment;
- no copy's map covers it, so attribution kept them by default.

038890fd (`--member-attribution-delta`) was built on the wrong mechanism and is not used.

**2. Copies within delta of each other are undecidable** (`--twins`; c3113f61, 7b8092d1):
- O1 marks a new copy within delta of a kept curated copy as its twin (genomic spans, asm20 `de`). O2 then lets the twin
  take a read only by the margin.
- O1 keeps a curated copy that its support rule would drop when it lies within delta of a kept one. Round 10's PRR20.1 had
  all its reads tied, a MAPQ-0 tie-break dropped it, and its reads went to the twin as errors.

**3. Delta is error-free, reads are not** (`--error-floor`; f2789fdd, 2afdbc3d).
- Delta (0.00958) is the p99 exonic divergence between KB3781's haplotype assemblies (register 1199).
- Every read-level comparison with it compared an error-inclusive `de`. The world's sequencing-error floor is now subtracted
  first, as the median `de` of its MAPQ-60 primaries (round 10: A119b 0.0027-0.0031, public 0.0015-0.0018). It applies to:
  - O2's `unexplained` gate (paths and genome fallback): explained up to (delta + floor) x read length;
  - O3's genome check, call-level and per member;
  - the scorer's per-read detectability: a REMOVED read is detectable when its `de` minus the floor exceeds delta, for both
    systems.
- It removes both remaining failure mechanisms:
  - set 10's short reads whose sequencing edits exceeded delta's share of their length;
  - array clusters whose members' `de` (0.0103-0.0118) is sequencing noise above a near-identical copy.

**4. O1-hosted O3 (`--o1-hosts`, f45e0fe5) is measured and not adopted.** It removes the integrated pipeline's last use of
copy_assign's run, but loses REMOVED and BRANCH recall in set 9's arrays (REMOVED -5 and -20, BRANCH -55 and -53). It
passed 6 of the 9 runs.

**v4 = v3 + 1 + 2 + 3, against the current method.** Correct / errors, floor-corrected per-read scoring:

| world | r10 A119b | r10 public | set 8 A119b | set 8 public | set 9 A119b | set 9 public | set 10 A119b | set 10 public |
|---|---|---|---|---|---|---|---|---|
| INTACT | 737/0 vs 645/0 | 654/0 vs 576/0 | 2,221/0 vs 2,167/0 | 2,254/0 vs 2,153/1 | 944/0 vs 907/0 | 876/0 vs 828/0 | 1,146/0 vs 1,138/2 | 1,059/0 vs 1,051/2 |
| CROSS | 779/0 vs 689/0 | 694/0 vs 616/0 | 2,742/1 vs 2,551/8 | 2,792/1 vs 2,558/9 | 1,037/0 vs 1,001/0 | 966/0 vs 917/0 | 1,265/1 vs 1,248/3 | 1,161/0 vs 1,146/3 |
| OUTSIDE | 982/83 vs 828/142 | 919/72 vs 802/115 | 2,264/82 vs 1,979/150 | 2,276/77 vs 2,000/120 | 777/43 vs 674/152 | 690/34 vs 605/123 | 1,105/8 vs 863/169 | 1,030/5 vs 789/167 |
| REMOVED | 969/0 vs 879/2 | 919/0 vs 843/3 | 2,126/0 vs 2,059/7 | 2,159/0 vs 2,072/4 | 802/0 vs 782/0 | 710/0 vs 698/0 | 1,032/1 vs 1,031/4 | 1,018/0 vs 963/4 |
| BRANCH | 1,043/59 vs 682/154 | 946/49 vs 576/139 | 2,536/110 vs 2,203/137 | 2,513/168 vs 2,176/191 | 1,205/36 vs 907/60 | 1,083/31 vs 828/46 | 1,329/74 vs 1,138/225 | 1,194/84 vs 1,051/200 |
| SILENT | 476/0 vs 464/0 | 439/0 vs 424/0 | 1,929/0 vs 1,860/0 | 2,023/0 vs 1,909/1 | 676/0 vs 644/0 | 629/0 vs 587/0 | 700/0 vs 689/2 | 609/0 vs 602/2 |

- **Floor-corrected scoring: PASS in all 9 runs,** including the untrimmed-library one (not tabulated).
- **Round 10's scoring: PASS in 5.** The 4 failures are all REMOVED B2: the wider gate places reads whose `de` lies
  between delta and delta + floor, which that rule counts as errors.
  - set 8 public 9 vs 7;
  - set 9 A119b 1 vs 0;
  - set 10 A119b 33 vs 10;
  - set 10 public 13 vs 4.
- Runs: `/mnt/linuxdisk/tmp/o3_graph/fast/r11dev6/` (`bars_floor.txt`, `bars_r10rule.txt`), from r11dev3's O1 output.

### 13. Rounds 11-13: v4 passes; v5 and v5' fail (2026-10-16/17)

| round | pipeline | families | verdict | what decided it |
|---|---|---|---|---|
| 11 (r1309) | v4 | 8, 37 copies | **PASS, 12/12** | the INTACT/CROSS/SILENT margins are the current method's own O3 artefact (review) |
| 12 (r1310) | v5 = v4 + O1-hosted O3 (`--o1-hosts`, tiered) | 8, 104 copies (80-copy DUX4) | FAIL, both profiles | O2: 2 margin-2 coincidences in DUX4, 5 structural errors (retained introns absorbed by a paralog's path as one gap-compressed edit; a read DNA-identical at two copies). O3 cost: 791-1,355 s per set |
| 13 (r1311) | v5' = v5 + `--branch-genome-novel` | 8, 89 copies (43-copy OR7E, 3 retrocopy families) | FAIL (A119b); public 6/6 | one false O3 call on three haplotype-2 reads of MRGPRX.1 in four worlds: a mosaic call whose divergence is mostly one read's sequencing indels, kept by attribution's margin of 2 that two donor alleles meet (review). No O2 error on catalogued reads |

- **v5 makes the integrated pipeline independent of the current method.** O3's hosts come from O1's seeded copies, tiered:
  primary alignment on a seeded copy, then a secondary, then a +-1 kb window.
- **v5' makes it cheap on large arrays.** O3's branch arm clusters only reads whose own genome alignment leaves >= 50 bases
  unexplained. On round 12's worlds this cut O3's longest world from 699 s to <= 78 s. In round 13 the integrated steps took
  223 and 216 s per set of six worlds; the current method took 188 and 183 s.
- **v6 (`--distance gapstruct`) is not adopted.** Gap compression applied only below 30 bp still fails round 12, and it costs
  errors and correct reads in sets 8 and 10.
- **Open after round 13:**
  - round 12's structural O2 errors, which round 13's families did not provoke;
  - the allele-versus-copy confound in O3's genome check (section 14).

### 14. After round 13: v8, then v9 under a diploid donor (2026-10-17)

Round 13's review traced the false call to three things together:
- a mosaic distance that counts one long read's sequencing indels at positions only it covers;
- attribution's fixed margin of 2, which two donor alleles meet exactly;
- the per-read genome check, the last filter.

Development uses 16 spent runs: sets 8-10 and rounds 10-13, both profiles, plus the two untrimmed-library runs. v5' passes 13
of them. The rule, fixed before any result was read: a lever is adopted only if it removes round 13's false call or reduces
errors, and turns none of the 13 passing runs into a failure (`fast/r14dev/<lever>/`).

| lever | round 13 A119b | runs passing | turned a passing run into a failure | cost |
|---|---|---|---|---|
| O3 `--allele-discount` v7f: donor alleles leave the genome check | fixed | 13 | set 10 A119b (REMOVED 1,032 -> 1,024 vs 1,031) | low-divergence removed copies whose reads O2 partly assigns |
| O3 `--member-attribution-delta`: attribution margin max(2, delta x length) | fixed | 13 | set 10 A119b (1,027) | REMOVED -1 to -34 in 11 runs |
| **O3 `--mosaic-shared`: the mosaic distance only over positions >= 2 members cover** | **fixed** | **14** | **none** | **none: identical or better in every run** (set 10 A119b REMOVED +49, round 11 A119b +11) |
| **O2 `--o2-genome veto`: a path call overruled where the genome places the read confidently** | (no O2 error) | **15** with `--mosaic-shared` | **none** | CROSS -1 and -2 in set 8 |

**The allele discount took four versions,** each lost on a different removed copy:
- v7 learned alleles from every explained read and lost CTSLP.1. CTSLP.1's own reads, which O2 explains on its paralog, carry
  its PSVs.
- v7b counted candidates per seeded locus, and CTSLP.1's PSVs were diluted over its paralog's length.
- v7c grouped co-carried candidates and lost CXADRP.1 and AGAP.5.
- v7f learned only from path-assigned reads, taking a variant as an allele when its assigned carriers outnumber the calls'
  members. It still lost two small calls of within-delta copies that are detectable only through high-error reads.

It stays opt-in, not adopted. A population cannot always tell a present copy's allele from a low-divergence missing copy's
PSV.

**The genome veto took three versions:**
- v1 scored alignments with `locus_distance`, which leaves a clip under 30 bp free. A paralog's alignment that clipped a read's
  differing end then beat the read's own locus, and true KRTAP2.1 calls were moved.
- v2 charged clipped ends and long insertions in full (`veto_distance`). Now a short terminal exon the aligner missed at the
  true locus moved a correct DDX11L call (primary MAPQ 1).
- v3 adds a guard: the read's primary alignment must lie on the genome's winner at MAPQ >= 30.
  - Audit of v2 over 13 runs: every wrong veto had MAPQ 0-1, and every veto that fixed an error had MAPQ 53-60.
  - It fixes round 12's four CDRT15 retained-intron errors (round 12 public now passes), cuts OUTSIDE errors in 8 runs, and
    adds no error anywhere.

**v8 = v5' + `--mosaic-shared` + `--o2-genome veto`: 15 of 16 development runs pass.**
- The one failure is round 12 A119b, with 2 DUX4 margin-2 coincidences and 1 KRTAP2 read that is DNA-identical at two copies,
  in INTACT, CROSS and SILENT. The current method abstains on all of them, and no sequence rule separates them.

**The simulator** (`--sim-het-draw`, round 13 review I1): each copy's het and hom site counts are drawn (Binomial over its
exonic bases) instead of rounded. With rounding, the public donor was effectively haploid: 86 of 89 copies had no het site.

**Under drawn het sites, v8 meets the allele confound again.** v8 was run on round 13's families re-simulated with
`--sim-het-draw` (`fast/r14dev/hetdraw13/`, both profiles).
- **Public passes.**
- **A119b fails INTACT 5 vs 0 and CROSS 6 vs 2.**
  - Every error is an O3 call on FTH1P.6's haplotype-2 reads.
  - FTH1P.6 is a 553-bp retrocopy that drew 3 het sites within 131 bp. Its short haplotype-2 reads (203-258 bp, MAPQ 60)
    carry all three, and at that length 3 alleles plus their sequencing errors clear delta plus the floor.
  - On FTH1P.4's window, the 5 shortest form a PSV call with 3 linked unexplained alts.
  - A copy-level view puts the alleles at 3 / 553 = 0.0054, under delta: a haplotype, not a copy.
- **v8's rules do not reach it.**
  - The mosaic fix does not apply, since this is a linked-PSV call.
  - The genome check reads each read's local divergence.
  - The genome veto is O2.

**O3 `--haplotype-check`.** A PSV call is measured as the haplotype it is, not as its members alone. Development took three
versions:
- **v1, genome coordinates.** The substitutions >= 2 members share, then every primary carrying one, over the positions >= 2
  carriers cover.
  - It removes MRGPRX.1 and FTH1P.6, and passes 15 of 16 development runs.
  - It misses an array. On round 10's families re-simulated with drawn het sites, 3 haplotype-2 reads of TAF11L.10 (a 45-copy
    array) are called on TAF11L.44's window. Their MAPQ-0 primaries lie on three different units, so in genome coordinates they
    share nothing.
- **v2, the host window.** The family's explained reads are aligned to the call's window, those carrying the call's alts are
  added, and the extended cluster's nearest-copy distance is taken over positions >= 2 of them cover.
  - Selecting carriers by any alt pulled in 861 reads of 12 FTH1P copies, since a mosaic cluster's alts are other copies'
    bases. The mixture looked as distant as the call, and true REMOVED calls were lost.
- **v3 = v2 with carriers selected by the call's unexplained alts only** (the alts no copy's map carries).
  - FTH1P.6: 61 FTH1P.6 carriers, distance 0.0019.
  - TAF11L.10's call is removed too.
  - Every true call in those worlds' REMOVED is unchanged.

**v9 = v8 + `--haplotype-check` (v3).**
- **On the 16 spent runs: 15 pass**, as with v8. The failure is again round 12 A119b. v9's O3 rules leave 13 runs identical to
  v5' and improve 3: round 13 A119b's false call removed; REMOVED +2, +11 (round 11 A119b) and +45 (set 10 A119b).
- **On 8 family sets re-simulated with the diploid donor: 8 of 8 pass.** These are rounds 9, 10, 11 and 13's families, both
  profiles (`fast/r14dev/hetdraw*`, `hetfresh*`, `hetc_v9c`).
- **Cost:** the check aligns the family's explained reads to each host window with a candidate call. O3's longest world on
  round 13 A119b went from 26 to 58 s.
- Round 14 pre-registers v9 with the diploid donor (`docs/PREREG_o3_graph_round14_2026-10-17.md`).

### 15. After round 14: v10, the O2 proof, and the cost of full families (2026-10-17)

**v10 = v9 + O3 `--haplotype-check-genome`.** This brings back the genome-coordinate haplotype check beside the window one
(f867ec7c).
- **On the 16 spent runs: 15 pass,** as with v9. The failure is again round 12 A119b. Per world, the results are identical to
  v9 or have fewer OUTSIDE errors (round 12 public: 130 -> 106).
- **On round 14's 5 seeds and 8 diploid-donor sets: all 18 runs pass.** Under v9, A119b passed only 2 of the 4 extra seeds.
- **These sets are now spent.** Round 15 is the held-out test.

**NPIP as one family** (families.py group `mainstay`): NPIPA and NPIPB copies link into one component under the locus rule
(23 of 23), and split in two they placed NPIPB reads on NPIPA copies. On the merged mainstay set (seed 306), both profiles
pass every bar.

**The O2 proof** (`docs/PREREG_o2_proof_2026-10-17.md`, register 1314) passed in every cell.
- The integrated O2 made 0 wrong assignments in 22,631, on 8 held-out near-identical families.
- MAPQ >= 1 is the strong baseline: it is more powerful on weak reads and makes 55 errors.
- The same runs' end-to-end OUTSIDE worlds failed: O1 did not seed near-identical uncatalogued twins (C4, HSFY), and O2's twin
  rule gave twin reads to the curated copy (VCY).

**Where a full family's time goes** (mainstay CROSS world, O3_TIMING, one world):
- **vg graph construction was 121 of O3's 158 s with a cold cache.** The cache is keyed by the graph's inputs, and each world's
  O1 catalog differs, so most worlds miss.
- vg msga is not deterministic either: two builds of NPIP's graph differ in copy17's paths.
- **O2:** 32 s, of which minimap2 of 1,829 NPIP reads against 307 paths is 23 s.

**`--branch-linear` (O3, 44798db2): the branch step without the family graph.**
- Candidates are spliced-aligned (minimap2 splice:hq) to the family's loci as linear sequences.
- A read base is novel only where NO locus alignment explains it (`novel_union`). That is how the graph explains NPIP.14's
  125-253-bp insertions: other copies' alignments hold them.
- **Mainstay CROSS world:** calls identical to the graph version; O3 158 s -> 34 s.
- **On all 34 spent runs** (16 rounding runs, round 14's 10 seed runs, 8 diploid-donor sets):
  - every verdict is v10's;
  - 14 runs have identical bars;
  - 19 gain 1-9 correct BRANCH reads (two REMOVED worlds +1), with errors equal or one fewer.
- **Measured and not adopted:** `--branch-genome-segments`, which clusters the genome alignment's own unexplained pieces. It
  called 32 present NPIP reads on the mainstay CROSS world, whose pieces another NPIP copy explains.

**Exact caches for a set's worlds (`fastset --aln-cache`):**
- **O2's path alignments** are cached per read sequence (`O2_ALN_CACHE`, f2e59107). minimap2 aligns each read on its own;
  equal-score ties are broken by a hash of the read's name, which changes only the tp / s1 / s2 tags, never NM, the CIGAR or
  the coordinates O2's distances read.
  - **Mainstay, against the original tables:** INTACT cold, then CROSS warmed by INTACT, then CROSS warm. All three are
    byte-identical.
- **O3's catalog-isoform-to-host alignments** (`O3_ALN_CACHE`, 7839138c) are keyed by window, query names and sequences.
- **Not cached:** O3's read-to-host alignments. They keep only minimap2's primary, so a name-hash tie could change which one
  is kept.

### 16. After round 15: the O2 margin (2026-10-17, evening)

**Round 15's errors, read by read.** Every one is a FAM90A read, and the current method makes none of them.
1. **Path-set asymmetry (most of them).**
   - FAM90A.45's read-derived paths (an exon skip) fit reads of its siblings, and the siblings' own paths lack that
     structure.
   - Each such read lay 16 edits from copy 45's paths and 64-378 from every sibling's. On the genome, the siblings hold it at
     5 and copy 45 at 16.
   - Union projection adds a paralog's isoform to a copy only for an intron chain new to that copy, so read-derived
     structures and longer first exons can stay on one copy.
   - The fixed margin of 2 edits rewarded a copy for having the read's path, not for its sequence.
2. **A locus-fallback call on clipped bases** (FAM90A.51's read): copy 45's genome alignment clipped the 11 bases that tell
   the copies apart. locus_distance leaves a clip under 30 bases free, so copy 45 won by 2.
3. **A call no rival's paths could answer** (FAM90A.12's 201-bp read): only copy 11's paths aligned it, and eight loci hold
   it exactly.

**Round 12 A119b's standing failure** (the one development run v10 and v11 never passed) adds two more mechanisms:
4. **Indel margins against many rivals** (DUX4, 80 copies).
   - Two path calls won by 2 edits that were two 1-bp homopolymer indels, against 63 rivals.
   - copy_assign ties these reads: its PSVs leave indels out by default (`indel_psv: false`).
5. **A retained intron** (KRTAP2.0).
   - The simulator keeps an intron in 5.8% of reads (the A119b profile, from real reads).
   - Copy 1's long annotated first exon absorbs the intron as one gap-compressed insertion (2), and copy 0's paths clip it
     (85). Both loci hold the read exactly.

**The error model behind mechanism 4** (the v2 profiles, from real FLNC reads):
- substitutions 4.1e-4 per base;
- indel events 2.1e-3 per base, 63% of them inside homopolymer runs of 3 or more;
- donor alleles 1.0e-3 per base.

A homopolymer indel is about 20 times likelier than any one substitution, yet the gap-compressed distance counted both as one
edit.

**Four rules, each in O2's genome step** (`o2.integv`; O3's input, `o2.integ`, is unchanged). Each was measured on all 46
spent development runs (276 worlds: rounds 10-15, sets 8-10, round 14's seeds, the diploid-donor sets and the O2 proof's
cells), against v11's table, with the production code (`o2graph.py run W integ|integv`; the integ step with
`--o2-significance` reproduced every original table exactly). Adoption rule: no world that passed may fail.

| rule | mechanism | correct | errors | verdicts |
|---|---|---|---|---|
| `--veto-abstain` | 1 | 0 | -95 (27 worlds) | round 15 A119b 503 FAIL -> PASS |
| `--locus-pairwise` (with veto-abstain) | 2 | -5 | -6 | round 15 A119b 502 FAIL -> PASS |
| `--o2-significance` | 4 | -421 | -14 | none alone |
| `--structure-check` (with veto-abstain) | 3, 5 | -868 | -24 | fixes round 15 public 503 CROSS |
| **all four (v12)** | 1-5 | **-1,294 (-0.21%)** | **-139 (-2.1%)** | **36 -> 40 of 46 runs pass; no world newly fails** |

- **`--veto-abstain` (1bd6693a): an interval rule.**
  - A copy's genome cost lies between its locus_distance (a clip under 30 bp free) and its veto_distance (every clipped base
    paid).
  - A path call is abstained (`genome_contradicts`) when a rival's worst case beats the call's best case by the margin, and
    no confident single winner takes the read. The veto's MAPQ >= 30 rule still moves reads with one confident winner.
  - **The first form counted the clip-paid cost alone.** It abstained correct calls whose own alignment missed a short exon:
    DDX11L.6's 15-bp terminal exon, and DAZ.2's read end.
- **`--locus-clipped` (not adopted):** the same interval rule on locus-fallback calls.
  - It abstained correct calls wherever every copy's alignment clips the read's ends (HNRNPCL, TAF11L, FAM236).
  - It failed r14 A119b s205's B1.
- **`--locus-pairwise` (10208f8c).**
  - Against each rival, a copy pays its alignment's edits plus the read bases the rival aligns and it leaves out. Bases every
    copy clips are nobody's evidence.
  - The one read it loses (ANP.1's) is the ambiguity no alignment settles: the true copy's alignment clipped 14 bases that
    the rivals align with 6-9 edits.
- **`--o2-significance` (1bd6693a).**
  - O2's path distances, with each edit priced at -ln of its type's rate in the world's uniquely placed reads: a mismatch to
    one given base, a homopolymer indel per homopolymer run, another indel per base. The CIGAR's context is read against the
    path and genome sequences.
  - Round 12 A119b: 7.8 for a mismatch, 3.9 for a homopolymer indel, 7.3 for another indel.
  - A path call stands only when it beats every rival by ln(R / alpha), with R the rivals and alpha = 1e-3, as copy_assign's
    certificate divides alpha by its rivals.
  - **Calibration, round 12 A119b INTACT:** the gate took 39 calls, all DUX4 on two homopolymer indels (7.7 against a
    threshold of 11.0); 2 were wrong. That is 5% wrong, against the model's R e^-7.7 = 3%.
  - It adds 1-3 CPU seconds per family world (NPIP 1.4 s, GOLGA8 1.2 s).
- **`--structure-check` (10208f8c).**
  - A rival whose own locus (an alignment inside its padded span) explains the read better than its paths do, by the margin,
    lacks the read's structure among its paths.
  - The call is then compared with that rival on the genome, and abstained (`structure_only`) unless its best explanation
    wins there by the margin.
  - DAZ.2's rivals align 40-58 kb past their loci, so they are no evidence.
  - It also abstains correct calls that rest on annotation structure the genome supports by one edit or not at all. The
    simulator rewards those calls, since its reads come only from annotated isoforms. The rule keeps the pipeline's own
    margin of 2 rather than a tuned exception.

**The PSV question, as far as O2 goes.**
- A SUN-type difference (a substitution) now weighs 7.6-7.8 nats and a homopolymer indel 3.9-4.6. The CIGAR's edit types
  and the read's own error profile set the margin, where the fixed 2 edits treated every edit alike.
- `de` still flags absent-copy reads (section 5). Path alignments' mismatches come from NM less the indel bases (their
  CIGARs carry M); the genome's from =/X.

**v12 = v11 + all four.** Over the 276 worlds:
- v12 makes 602,563 correct calls and 6,478 errors, against v11's 603,857 and 6,617.
- The current method makes 556,431 and 24,141.
- **Four failing runs pass:** round 15 A119b 502 and 503, round 15 public 503, and round 12 A119b. Round 12 A119b is the one
  development run that v10 and v11 never passed.
- **The six still failing are the O2 proof's cells,** on their OUTSIDE twins (register 1315).
- **The thinnest B1 margins left are 3-6 reads** (integrated correct less current correct):
  - the O2 proof's public 402 INTACT and CROSS: 13 -> 3 and 15 -> 5;
  - het9 INTACT and SILENT: 6-8 -> 4-6;
  - set 10 public SILENT: 7 -> 5.
  B1 is where v12's cost could fail a held-out world.
- The structure check and the significance gate cost about 30-36 correct reads per error they remove. That is the price of
  holding errors at zero where the current method abstains.
- **Round 16 tests v12 on families it has not seen** (`docs/PREREG_o3_graph_round16_2026-10-17.md`).

**Round 16 (register 1318): FAIL on correctness (4 of 6 cells pass; A119b 601 and 602 fail); PASS on cost.**
- **B2 held in all 36 worlds.** The O2 margin removed every error where the current method makes none. Over the 36 worlds:
  1,833 errors against 8,796, and 147,583 correct calls against 137,426.
- **The three failing bars are B1, by 4-7 reads,** in A119b INTACT and SILENT of small families (3-7 copies).
- **The integrated pipeline's fixed thresholds are the cause, not v12's rules** (v11's genome step fails the cell as well):
  - `unexplained` short reads, whose edits exceed (delta + floor) x length;
  - `thin` margin-1 calls among few rivals, which copy_assign's certificate assigns.
- **The cost bar passed in CPU** (-23% and -30% on the mainstay set). On round 16's families the integrated stages are also
  faster on wall time (68-100 s against 118-138 s per cell).

**After round 16 (register 1319): v13 = v12 + `--struct-rescue --short-path-align`.**
- **What round 16's missing reads were:** reads the current method placed and v12 did not.
  - Short reads over the fixed divergence gate.
  - Thin margin-1 calls. Gap compression let the read's best path be an isoform lacking one of its exons, so the paralogous
    bases inside that exon went uncompared.
  - 150-200-bp reads with no path alignment at all: map-hifi's minimum DP score of 200 is out of their reach.
- **The adopted fixes, on 52 spent runs:**
  - Thin and tied reads are decided again on structure-aware distances (gapstruct: a long indel counts its length).
  - Reads with no path alignment are realigned with -s 40.
  - The result: 44 -> 46 runs pass, round 16's three failing worlds pass, no world newly fails; +4,289 correct, +95 errors in
    OUTSIDE and BRANCH.
- **Not adopted:**
  - A length-calibrated gate: it places reads of absent copies.
  - A significance-based margin: about one error per 1,000 rescued reads, where the current method makes none.
- **Round 17's families are name-free** (`seqfamilies.py`, user, 2026-10-17). They are every protein-coding gene, LOC genes
  included, linked by sequence identity (the locus rule) and by structure (the alignment carries one gene's junctions onto
  the other's introns).

**Round 17 (register 1320): FAIL on correctness (0 of 6 cells), PASS on cost.**
- **The setup:** name-free families (40, mostly pairs).
- **What held:** B1 in all 36 worlds; v13 made 124,264 correct calls against 82,008.
- **What failed:** OUTSIDE's B2 in every cell. When a pair's twin is uncatalogued and O1 does not seed it, O2 places the twin's
  reads on the one remaining copy (the O2 proof's open failure). The current method abstains on them.
- **Cost:** the integrated stages used 20-29% less CPU on the NPIP-sized set, and 43% less on round 17's cells.

**After round 17 (register 1321): v14 = v13 + O1 `--ref-copies`.** O1 searches the genome for copies the catalogue lacks (the
family rule on the world's homology-closed local reference); round 17's spent cells went 0 -> 5 passing.

**Round 18 (register 1322): FAIL on correctness (4 of 6 cells pass; both seed-803 cells fail), PASS on cost.**
- **What held:** B2 in all 36 worlds, OUTSIDE included (0-11 errors against 69-198). v14 made 130,722 correct calls against
  86,384, and 3,782 errors against 10,431.
- **What failed:** B1 in five worlds of seed 803, all from one family (OR4F21, OR4F16).
  - The genome search adds five OR4F unprocessed pseudogenes, 98.8-99.4% identical, which the simulator never expresses.
  - v14 abstains the reads tied with one, or one substitution from one. The current method, which compares catalogued copies
    only (`--no-as-tied-only`), places them on OR4F21, and in this simulation it is always right.
- **The same tie in another form failed `--twin-guard` (register 1321):** a read tied with an expressed uncatalogued twin and a
  read tied with a silent paralog are the same exact AS tie.
  - The current method wins the second kind because the simulation never expresses paralogs outside the catalogue.
  - Real A119b testis expresses these five pseudogenes: 41 primary reads, against 59 on the two coding copies.
- **Cost:** -3% and -13% CPU on the NPIP-sized set. O1's genome search now costs 141-147 CPU seconds there.

### 17. The goal's PSV question, answered by objective (2026-10-18)

The goal asks whether PSVs can include SUNs, and how the CIGAR, NM and de tags should count. The answers, from sections 4, 5
and 16 and rounds 16-18:

**O1, which copies exist.**
- **The CIGAR's introns are O1's structure evidence.** The genome search (`--ref-copies`, register 1321) takes a locus as a copy
  only when a catalogued transcript aligns to it at identity >= 0.90 over >= 0.8 with its exon chain carried: every N of the
  spliced alignment at one of the transcript's junctions.
- **A SUN is private only relative to the copies considered** (section 4). Completing the catalogue from the genome makes
  "private" mean private among the loci the family rule reaches. Round 18 showed the cost (register 1322): five OR4F
  pseudogenes made OR4F21's catalogue-private differences shared, and its reads tied.

**O2, which copy a read came from.**
- **SUNs are most of the evidence:** about 70% of the decisive PSV columns are SUN-type, the called copy's own SUN or a rival's
  SUN the read rules out (section 4). With a complete catalogue no call was wrong, with or without an own-SUN column.
- **The CIGAR's edit types set the margin** (`--o2-significance`, section 16). Each edit is priced at -ln of its rate in the
  world's uniquely placed reads:
  - a substitution 7.6-7.8 nats;
  - a homopolymer indel 3.9-4.6;
  - any other indel 7.3.
  - A call must beat every rival by ln(R / alpha). One SUN decides against one or two rivals; one homopolymer indel never.
- **Indels count as events, not bases** (section 5). NM counts a retained intron or a 4-bp homopolymer indel as many edits. A
  gap-compressed distance (mismatches + indel events, like de) recovered the calls NM lost.
  - A long indel (>= 30 bp) counts its length again (`--struct-rescue`), because it can hide the paralogous bases inside a
    skipped exon (SPAG.1, register 1319).
- **Soft clips bound the genome cost from both sides** (`--veto-abstain`): a clip under 30 bp is free at the lower bound and
  paid in full at the upper.

**O3, reads of copies the catalogue lacks.**
- **de flags them.** A read whose de, less the world's error floor, exceeds delta is detectable as an absent copy's read (the
  scorer's per-read detectability, section 5).
- **The read's own de cleans the calls:** present-copy reads inside a true absent-copy cluster sit within delta of their own
  locus (51 of 55; `--member-genome-check`).

**What it bought, measured on held-out families:** B2 (errors no more than the current method's) held in every world of rounds
16, 17 (except OUTSIDE, fixed by `--ref-copies`) and 18. Over round 18's 36 worlds: 3,782 errors against 10,431. The open bar
was B1 in worlds with silent paralogs (register 1322). Round 19 passed it under the measured paralog expression (register
1324).

**Round 19 (register 1324): PASS on correctness (6 of 6 cells, all 72 bar checks) and on cost.**
- **The setup:** v14 unchanged; 35 fresh name-free families whose non-coding members are copies like any other (90 copies, 36
  non-coding), expressed at the level each testis library shows (profiles v4).
- **Correctness:** 85,440 correct calls against 62,070 and 2,355 errors against 8,028. The thinnest B1 margin is +29 reads; B2's
  thinnest are ties at zero errors in complete catalogues.
- **Cost:** -4% and -13% CPU on the NPIP-sized set. On round 19's small families the integrated stages use more CPU (849 against
  750 s), because O1's genome search costs about 20 CPU seconds a world.
- **What it rests on:** round 18's failure came from the simulation's silent paralogs. With them catalogued and expressed as real
  testis expresses them, the current method's errors rose (seed 803 INTACT: 34 -> 121 and 29 -> 134) while v14 stayed at 0
  (register 1323; not diagnosed read by read).


### 18. The variation graph is the model; linear alignment is the implementation (2026-10-18)

**The model (user, 2026-10-18).** A family is a variation graph of its copies' loci, each curated isoform a path through it.
- O2 assigns a read to the copy whose paths explain it, by the PSVs where the copies' paths differ.
- O3 calls a read of an absent copy where no path explains it: a novel sequence, or a switch between copies' paths.

**The implementation evaluates the same paths without building the graph,** for the goal's cost bar:
- **O2:** every isoform path of every copy (the paths `vg rna` put in the graph, plus paralog structures projected onto each
  locus) is a linear sequence; minimap2 aligns the reads to all of them, and the call is made on those path distances.
- **O3 (`--branch-linear`, 44798db2):** candidate reads are spliced-aligned to the family's loci as linear sequences. A base is
  novel only where no locus alignment explains it, which is the linear form of a path switch.

**The measured equivalence:**
- On the mainstay CROSS world the linear branch step's calls are identical to the graph's, and O3 falls from 158 s to 34 s.
- On 34 spent runs every verdict is the graph version's, and 19 runs gain 1-9 correct BRANCH reads.
- vg (`vg msga`, `vg rna`, GraphAligner) was 77% of O3's time on a full family's first world.

The graph code remains (`detect.py` without `--branch-linear`; `o2graph.py run WORLD nodes`), so the model can be run as stated
when a result must be shown on the graph itself.


---

## O1O2O3_INTEGRATED_OVERVIEW_2026-10-16.md

*(original: `docs/archive/2026-10/merged/O1O2O3_INTEGRATED_OVERVIEW_2026-10-16.md`)*

## O1 -> O2 -> O3 integrated: the overall picture (2026-10-16)

> Update 2026-10-19. **O2's domain (user, 2026-10-19): AS-tied multimapping reads only** — a read with a rival alignment within 10 AS of its best (Eichler's within-10 rule; under splice:hq, 10 AS is two substitutions). A read with a clear best placement is already placed by minimap2: applying O2 to it is tautological. O2's point is to go beyond minimap2's shortcoming, its coin toss among AS-tied alignments, and assign those ambiguous reads to very similar copies by their PSVs, junctions and divergence, or abstain. Implementation: `fastset --o2-scope tied --as-slack 9` (`o2graph.py run W intega`), on test in round 21 (`docs/PREREG_o3_graph_round21_2026-10-19.md`); memory `feedback_o2_scope_ambiguous_only`.

User (2026-10-16): "lets look at the overall picture first, probably it is better to look at O1, O2 and O3 integrated since
that is the end goal".

This is a post hoc summary across every simulated world built in rounds 1-7 (`docs/PREREG_o3_graph_*`, register 1297-1305).
No new claim is made here; the next step proposed below is a pre-registered end-to-end test.

### The integrated pipeline (frozen at round 7, `1e513942`)

1. **O1, RNA-supported seeding** (`o1seed.py`):
   - keep each curated copy that has primary-alignment support, and drop the silent ones;
   - add off-catalog loci with >= 5 uniquely placed reads;
   - give each copy the intron chains of near-unique reads (mapq >= 50) that it owns.
2. **O2, padded union** (`o2graph.py run W integ`): every family structure is threaded through every seeded copy, with path
   ends padded with the copy's own genome. Each read is assigned at margin 2, or abstains, or is `unexplained` beyond delta.
3. **O3, arm I** (`detect.py`):
   - clusters what O2 leaves unexplained, pooled over the family, with the mosaic rule;
   - adds branch calls from the family's vg graph;
   - keeps a call only when the genome does not explain its reads.

The comparison is the current method: copy_assign's O2 plus its `--psv-clusters` O3 on the curated catalog.

### End to end: every read's fate (`summary.py e2e`)

Each read is classed by the copy it came from:
- catalogued;
- uncatalogued but in the genome;
- absent and detectable;
- absent and within delta of the genome, so unknowable.

Its fate under each system is one of: assigned to its own copy (correct), to another copy (error), inside an absent-copy call
(correct only for an absent copy), or left unassigned. Pooled over 43 worlds, 8 family sets:

| world type | current: correct / error / unassigned | integrated: correct / error / unassigned | correct | errors |
|---|---|---|---|---|
| INTACT | 17,990 / 483 / 5,727 | 19,377 / 10 / 4,813 | x1.077 | -98% |
| CROSS (isoforms annotated on a paralog) | 18,470 / 508 / 6,954 | 20,002 / 22 / 5,908 | x1.083 | -96% |
| OUTSIDE (uncatalogued copy in the genome) | 15,084 / 2,762 / 6,354 | 19,550 / 304 / 4,346 | x1.296 | -89% |
| REMOVED (copy absent from the reference) | 16,116 / 623 / 6,128 | 18,065 / 26 / 4,457 | x1.121 | -96% |
| BRANCH (absent copy with a novel exon) | 11,418 / 1,955 / 5,064 | 15,008 / 83 / 3,346 | x1.314 | -96% |
| SILENT (a third of the copies unexpressed) | 3,085 / 46 / 956 | 3,403 / 0 / 684 | x1.103 | -100% |
| **all** | **82,163 / 6,377 / 31,183** | **95,405 / 445 / 23,554** | **x1.161** | **-93%** |

Unknowable placements (reads of absent copies within delta of the genome) were 1,333 for the current method and 1,652 for the
integrated pipeline.

- **Per world.** The integrated pipeline has at least as many correct fates and at most as many errors in 38 of the 43
  worlds. The five exceptions:
  - round 3's INTACT and CROSS: 1 and 2 fewer correct;
  - round 4's INTACT and CROSS: 5 and 12 fewer correct;
  - round 6's REMOVED: 227 more correct, 1 more error.
- **Not tuned on.** The 30 worlds never used to develop the integrated pipeline (the original development set, rounds 1, 2,
  3 and 5 without alleles, and round 7's held-out set) give x1.178 correct and 92% fewer errors. The 13 worlds it was
  developed on give x1.122 and 98% fewer.

### Where the errors come from

All 43 worlds, by the same rules as `summary.py e2e` (each error counted once):

| error source | current | integrated |
|---|---|---|
| a present copy's read inside an absent-copy call (O3) | 3,412 | 61 |
| a detectable absent copy's read placed on a present copy | 1,695 | 82 (all novel-exon, BRANCH) |
| an uncatalogued copy's read placed on another copy | 1,213 | 291 |
| a catalogued copy's read placed on another copy (O2) | 57 | 11 |
| **all** | **6,377** | **445** |

There were 109,493 reads from catalogued copies and 4,335 from uncatalogued ones.

### What each objective contributes

- **O3 makes the largest change.** More than half of the current method's errors (3,412 of 6,377) are absent-copy calls
  on present copies' reads.
  - Arm I made 117 true and 6 false calls over the 43 worlds.
  - All 6 false calls are one cluster: NPIPA.5's reads hosted on the sister family NPIPB in round 1's held-out set, just
    above delta (median `de` 0.0097, known since round 4).
  - The earlier arms made the same number: H made 5 false calls (the same NPIPB cluster); K made 5 (one WASH cluster in each
    round-7 world, 5 reads each), which arm I does not make.
  - Branch calls, the genome check, the mosaic rule and family pooling catch detectable absent copies' reads: those placed
    on a present copy fall from 1,695 to 82.
- **O1 fixes the catalog.** Uncatalogued copies' reads placed on another copy fall from 1,213 to 291. With seeding they go to
  their own copy: 72% of them on round 7's held-out set (348 of 481).
  - O1 dropped 22 of 26 silent copies and none of the 1,740 expressed ones.
  - Its observed structures supply isoforms annotated on a paralog (CROSS).
- **O2 was already precise on catalogued copies** (57 wrong of 109,493 reads); the graph O2 brings that to 11 and abstains
  less.

### What no step fixes

- **Near-identical copies whose reads tie at mapq 0.** They cannot be seeded, told apart, or dropped (CTAGE.2 and DDX11L.9 on
  round 7's set). Reads alone cannot resolve them.
- **Alleles versus copies below delta** (round 6). The same signature needs the individual's own haplotypes. HG002 has them;
  its Kinnex reads are being aligned on the cluster.
- **Removed copies near delta** (WASH at 0.0124, round 7), at the boundary of detectability.
- **A sister family's copy hosted at delta** (NPIPA.5's reads on NPIPB, median `de` 0.0097): its reads look like an absent
  copy. This is the only source of the integrated pipeline's false calls.
- **One bug, known:** an uncatalogued locus can be seeded into two adjacent families (WASH and DDX11L). A seeded locus must
  belong to one family. *Fixed in 1e3960fe; see the update below.*

### Why the per-objective bars failed although the whole works

The pre-registered O2 bars compared O2 alone with the current O2 alone:
- they did not count the current method's O3 calls on present copies, its largest error source;
- they counted a near-identical uncatalogued copy's reads placed on a paralog against O2, although no method can separate
  them.

End to end, those reads move from "error" (current: called absent, or placed on the wrong copy) to "correct" (seeded) or to
the same unresolvable placement.

### Next

1. **A pre-registered end-to-end test on a fresh family set.** Bars on read fates per world type: correct >= current,
   errors <= current. Include the one-family seeding fix and BRANCH as well as the round-7 worlds.
2. **Real data: HG002 Kinnex.** Truth for each read comes from its unique alignment to HG002's own maternal and paternal
   haplotypes (`/mnt/linuxdisk/tmp/rustle_figures_dev/o2_presence/hg002/hg002v1.1.fa`). Copies present in HG002 but absent
   from CHM13 give real O3 targets. It needs the Kinnex reads aligned to CHM13 (the pipeline) and to the HG002 assembly
   (the truth).

### Update after round 8 (2026-10-16)

**Two fixes since the tables above were written:**
- the one-family seeding rule (1e3960fe). Set 1's NPIPA/NPIPB locus had also been seeded twice, in all five of its worlds;
- the e2e scorer: a call on a seeded copy that matches no truth copy is now a placement, not an abstention (round 8 review I3).

**The same 43 worlds, re-scored:**
- integrated: 95,418 correct and 432 errors (was 95,405 and 445):
  - set 1's OUTSIDE: 3 more correct and 13 fewer errors;
  - set 1's CROSS: 10 more correct;
- REMOVED: 76 reads moved from abstain to unknowable (unknowable 1,728, abstain 4,381). No correct or error count changed;
- the integrated pipeline's error sources: uncatalogued 281 (was 291), absent-copy calls on present copies 58 (was 61),
  novel-exon 82, catalogued 11;
- arm I's false calls: 5, all of them the NPIPB cluster. The fix removed one.

**The worlds of one set share most of their reads:**
- INTACT, OUTSIDE and REMOVED simulate the same reads;
- CROSS and BRANCH contain all of INTACT's reads;
- SILENT's reads are a subset.

Pooled totals are therefore not independent counts; each world is the unit of comparison.

**Round 8, the pre-registered end-to-end test on a ninth set: FAIL by one read** (CROSS: 1 error against 0). B1 held in all six
worlds. See `docs/PREREG_o3_graph_round8_2026-10-16.md`.

### Update after rounds 9-11 (2026-10-16, evening): the integrated pipeline v4 passes on held-out full families

| round | pipeline | reads | verdict | what decided it |
|---|---|---|---|---|
| 9 (r1307) | v2 | calibrated v1 profiles | FAIL | 1-3.5% fewer correct catalogued reads (O2 `unexplained` gate on realistic reads); O1 dropping expressed copies |
| 10 (r1308) | v3 | calibrated v2 profiles | FAIL | one false O3 call (a host-window artefact in a 45-copy array); a new copy within delta of a catalogued one |
| 11 (r1309) | v4 | calibrated v2 profiles | **PASS, 12/12 worlds** | also under round 10's scoring, at copy level, and without H2AC/H2BC. The review found the INTACT/CROSS/SILENT margins are the current method's O3 artefact; with O3 removed from both, O2 is at parity under A119b (+3 to +5) |

**v4, end to end** (`fastset.py ... --o2-distance gapcomp --member-genome-check --member-attribution --union-exact
--o2-genome fallback --per-read-detectability --twins --error-floor`):
1. **O1** (`o1seed.py --twins`):
   - **Seeding:** the RNA-supported catalog, keeping a curated copy with >= 5 primaries or one owned uniquely placed read.
   - **Undecidable copies:** a curated copy within delta of a kept one is never dropped.
   - **Twins:** a new locus within delta of a curated copy is its twin.
2. **O2** (`o2graph.py run W integ 2 --distance gapcomp --union-exact --error-floor`, then `integl 2 --error-floor`):
   - **Paths:** every isoform structure is threaded through every seeded copy. Structures still project when the aligner
     drops a terminal micro-exon, and stay apart unless their junctions match exactly.
   - **Distance:** gap-compressed (indel events, clipped ends aligned against the path).
   - **Decision:** at margin 2. A twin wins a read from its curated copy only by the margin.
   - **Genome fallback:** a read the paths leave unexplained is placed on its own genome alignments (mismatches + indel
     events). Novel sequence >= 30 bp refuses, and the best alignment off the seeded copies keeps the read unexplained.
3. **O3** (`detect.py --member-genome-check --member-attribution --error-floor --arms I`):
   - **Clustering:** O2's unexplained reads, pooled over the family, with PSV clustering and the mosaic rule.
   - **Branch calls:** novel segments on the family's vg graph.
   - **Attribution:** a call keeps only the reads its consensus explains better than every catalogued copy that competes
     for them.
   - **Genome check:** per read, less the error floor.
4. **Delta and the error floor.** Delta (0.00958) is an error-free divergence: p99 between KB3781's haplotype assemblies.
   Every read-level comparison with it subtracts the world's sequencing-error floor first (median `de` of MAPQ-60 primaries).

**What remains (round 11's review and its own report):**
- **O3's read hosting still comes from copy_assign's run.** `--o1-hosts` removes that dependence, but in development it lost
  O3 recall in arrays; a tiered version is being measured.
- **Errors left in round 11:**
  - uncatalogued reads placed on a paralog (OUTSIDE 44 and 37, against the current method's 367 and 297);
  - BRANCH reads lacking the novel exon, which count as detectable by construction (118 and 119, against 313 and 297).
- **Not yet tested:** real reads end to end. Rounds 9-11 are simulations calibrated on two real libraries.

### Update after rounds 12-13 (2026-10-17): v5 and v5' fail

| round | pipeline | families | verdict | what decided it |
|---|---|---|---|---|
| 12 (r1310) | v5 = v4 + O1-hosted O3 | 8, 104 copies (80-copy DUX4) | FAIL, both profiles | O2: retained introns absorbed by a paralog's path (gap compression), a read DNA-identical at two copies, 2 margin-2 coincidences in DUX4 |
| 13 (r1311) | v5' = v5 + branch genome-novel prefilter | 8, 89 copies (43-copy OR7E) | FAIL (A119b); public 6/6, but the public donor is effectively haploid | one false O3 mosaic call on a present copy's haplotype-2 reads, on a paralog's window: one read's sequencing indels make its divergence, two donor alleles meet attribution's margin of 2 (review) |

- **The integrated pipeline is now independent of the current method:** O3's hosts come from O1.
- **It is cheap on full families:** 216-223 s for the integrated steps of a six-world set of 89 copies, with O3's longest
  world 72 s at 0.17 GB.
- **What remains is two confounds a sequence-only rule cannot always settle:**
  - an allele of a present copy against a low-divergence missing copy (O3);
  - a read's structure against a paralog's annotation (O2).
- **Development since round 13** (`docs/O1O2O3_FAST_REALISTIC_2026-10-16.md` section 14):
  - O3 `--allele-discount`: donor alleles, learned from the reads O2 assigns, leave the genome check;
  - O2 `--o2-genome veto`: a path call that the read's own genome alignments contradict by the margin is overruled.

### Update after round 14 (2026-10-17): v9 passes, with a diploid donor

| round | pipeline | families | verdict | what decided it |
|---|---|---|---|---|
| 14 (r1312) | v9 = v5' + `--mosaic-shared` + `--haplotype-check` + `--o2-genome veto` | 8, 74 copies (SD pseudogene families, a protocadherin cluster, olfactory subfamilies, 25 retrocopies) | **PASS, 12/12 worlds, but one fragile draw** (review) | 0 errors in INTACT, CROSS, SILENT and REMOVED; but at 4 further seeds A119b fails 2 (the allele confound again), public CROSS rests on the haplotype check, and the veto and mosaic-shared were not exercised |

**v9, end to end** (`fastset.py ... --o2-distance gapcomp --member-genome-check --member-attribution --union-exact
--o2-genome veto --per-read-detectability --twins --error-floor --o1-hosts --branch-genome-novel --mosaic-shared
--haplotype-check`). v4's description above still holds, plus:
- **O1** seeds the RNA-supported catalog; nothing from the current method.
- **O2: the genome veto.** After the paths and the genome fallback, a read the paths assign is overruled when its own genome
  alignments place it on another seeded copy by the margin, with its primary there at MAPQ >= 30. Retained introns are
  contiguous on the read's own locus.
- **O3**, hosted by O1:
  - the mosaic distance counts only positions >= 2 members cover;
  - a PSV call extended by the family reads O2 explains that carry its unexplained alts must stay > delta from every copy;
  - otherwise it is a present copy's haplotype.

**What remains:**
- **The allele confound (round 14 review).** It is the deciding open problem. A present copy's donor-allele reads, called on
  another locus's window, beat the genome check by 0.0005-0.0025 and the haplotype check by 0.00024 at further seeds of
  round 14's families under A119b. Two seeds of four fail. One arose on a copy O1 seeded at a retrocopy family's parent gene.
- **Round 12's A119b failure class:**
  - margin-2 coincidences among 80 near-identical units;
  - a read DNA-identical at two copies.
  The current method abstains on both, and no sequence rule separates them.
- **The secondary scorings:** under round 10's rule and at copy level, round 14 A119b fails REMOVED B2 (1 and 3 vs 0).
- **Not yet tested:** real reads end to end. Every round is a simulation calibrated on two real libraries.


---

## HANDOFF_cloud_2026-10-21.md

*(original: `docs/archive/2026-10/merged/HANDOFF_cloud_2026-10-21.md`)*

## Handoff — local WSL machine → cloud (2026-10-21)

### Branch

Clone the repository and work on **`o3-fast`**: `git clone git@github.com:juanfraitu1/Rustle.git && git checkout o3-fast`.
- It contains `main` (`github/main` is an ancestor) plus 283 commits: the O2 / O3 work of 10-16 to 10-21.
- It is not merged into `main`.
- Registers 1347-1369 and their pre-registrations are on it; the latest summary is
  [`MULTIMAPPER_PREMISE_2026-10-21.md`](MULTIMAPPER_PREMISE_2026-10-21.md).

### Environment

- **Tools:** `pixi install` (see `pixi.toml`) gives minimap2 2.30, samtools 1.22.1, pysam, numpy and scipy.
  - `pixi install -e bench` adds gffcompare, stringtie and liftoff.
  - `pixi install -e isocon` adds IsoCon 0.3.3, with its own Python 3.8.
- **FLAIR 3.0.0** is not in pixi. It is not needed to continue: its outputs already exist (below).
- **minimap2-gate:** `pixi run build-minimap2-gate` (or `tools/minimap2_gate/build.sh`; git, make, a C compiler, zlib headers).
  The families stage uses it to align only the locus pairs that can become edges: about 8x faster, byte-identical products
  (`docs/O1_FAST_2026-10-21.md`). Without it the stage runs the full all-vs-all, with a note.
- **Rust:** `pixi run build` (`CARGO_TARGET_DIR` defaults to `$HOME/rustle_target`; always `--release`; capture the log to a file).
- **Python tests:**
  - `python3 -m pytest -q bench/o3_graph` (481 pass at 3290e8a7);
  - `bench/rna_allele/test_npip_o3.py` (8);
  - `bench/tool_difference/test_s1_o2.py` (10);
  - `bench/o1_sd/test_sd_pairs.py` (5).
- **Paths:** scripts carry the local paths below (`/mnt/linuxdisk/...`, `/mnt/c/Users/jfris/Desktop/...`). On the cloud, either
  recreate those paths with symlinks, or pass the new paths where the scripts take arguments. The few hard-coded constants are:
  - `bench/rna_allele/npip_o3.py`: `BAM`, `GWFAM`, `TRUTH`, `GGO`;
  - `bench/o1_sd/sd_pairs.py`: catalogs and SD tracks;
  - `bench/o3_graph/donorsim.py` / `o2idealsim.py`: none (arguments only).

### State at handoff

| register | result | document |
|---|---|---|
| 1360, 1361 | the sample-aware corrector on NPIP / LRPAP1 tied reads in KB3781: UNDERPOWERED (final key) | `PREREG_o2_corrector_families_2026-10-20.md` |
| 1362 | the same corrector on ideal-expression reads from both haplotypes: FAIL (reads carry no positional evidence) | `PREREG_o2_corrector_idealsim_2026-10-21.md` |
| 1363 | reference-free O3 chain on KB3781 NPIP: finds the maternal-only units; FAIL on read recovery | `PREREG_o3_npip_maternal_copy_2026-10-21.md` |
| 1364, 1365 | candidates added to the reference: reads return (33 / 39), 0 false moves | `PREREG_o3_npip_candidates_realign_2026-10-21.md` |
| 1366 | premise check on existing StringTie / FLAIR runs (descriptive) | `MULTIMAPPER_PREMISE_2026-10-21.md` |
| 1367 | multi-copy gains are not shown to be real copies; the multi-mapping loci are SDs | `PREREG_multicopy_gains_real_2026-10-21.md` |
| 1368 | O1 family copies vs external SD calls: bar FAILS (72.1% / 82.4%); gene-like copies 93.8% / 92.0% paired; shortfall = transposon / satellite families and self-overlap families (an O1 defect) | `PREREG_o1_families_sd_pairs_2026-10-21.md` |
| 1369 | the self-overlap fix `--reject-self-hits`: v1 FAIL; v2 PASSES held-out gorilla_OR6737 (audit: holds with caveats); human_A119b fails guard G1 (development); now the DEFAULT (user decision) | `PREREG_o1_self_hits_2026-10-21.md` |

### Open items

1. **O1 self-overlap families** (registers 1368, 1369): fixed by `mcl_families --reject-self-hits`, **the default since
   2026-10-21** (user decision; `--no-reject-self-hits` restores the old graph) (`docs/PREREG_o1_self_hits_2026-10-21.md`).
   - The pre-registered v1 FAILED (gorilla B2, by 4 copies under a gene_id join).
   - v2 passed the three development catalogs.
   - Amendment 3 fixed the evaluation in advance: B2 by locus, plus guards G1-G3.
   - **Held-out gorilla_OR6737 PASSES every bar and guard** (audit: holds with caveats).
   - Under Amendment 3, human_A119b (development) **fails G1**: same-strand copy pairs split across families go
     1,036 → 1,140 > 1,088, mostly at loci already split on chr13.
   - **Supported:** v2 as the default for gorilla de novo catalogs. The user made it the default for all catalogs.
   - **Not shown** for human or for annotation-built catalogs.
   - **Next:**
     - a held-out human catalog (HG002 Kinnex after a whole-genome alignment) and a guided CAT arm, both under Amendment 3;
     - a transitive same-unit test (v3, its own pre-registration) for the one collateral split (MCL70).
2. **Classify the O1 copies that no SD pairs:** a quarter of the human high-identity copies, an eighth in gorilla. Check
   whether they are retrocopies (a single-exon copy with a multi-exon partner), groupings on short shared segments, tandem arrays,
   or repeat-driven. Pre-register first.
3. **Rerun the human Figure 3 table on CAT/Liftoff v2.0.** Register 1366's human numbers are scored on NCBI RefSeq.
4. **Generalize O3 beyond one locus:** a held-out reference-absent copy, then candidates plus realignment as in registers
   1363-1365.
5. **Older open threads:** `docs/PENDING_2026-10-11.md`, `docs/PENDING_2026-10-04.md`.

### Working rules this project follows

- **Pre-register first.** Write the decision rule (`docs/PREREG_*.md`, committed) before looking at results; change it only by
  a dated amendment.
- **Audit before reporting.** A fresh auditor (the `metric-auditor` agent) checks each verdict on a frozen, read-only copy with
  SHA1SUMS. Record the outcome in `docs/NEGATIVE_RESULTS_REGISTER.md`, and check that register before proposing anything.
- **Data invariants:**
  - never pool human and gorilla numbers;
  - human annotation = CAT/Liftoff v2.0;
  - CHM13 v1.0 ≠ v2.0;
  - primaries via `-F 2308`;
  - SAM is 1-based, exon tables 0-based;
  - gorilla SDs come from SEDEF, not our own `sd-scan` (register 1278: its stage-2 identity cut is a no-op).
- **O2 scope:** O2 = AS-tied multimappers only (a rival within 10 AS), allele-blind. A donor's own assembly may serve only as the
  answer key.
- **Families:** the multigene family (default `--rule multigene`) is the unit.
- **Push** only when asked.
- **Local-machine rules, not needed on a normal cloud machine:** the WSL crash rules (foreground, one heavy job at a time, no
  `pkill -f`) and the sparse-vhdx disk warnings.
- **Agent memory is not in this public repository.** To carry the local agent's notes over, copy
  `~/.claude/projects/-mnt-c-Users-jfris-Desktop/memory/` (2.8 MB) privately. The repository's `docs/MEMORY_DIGEST.md` holds the
  older compacted parts.

### Data to copy (local path → size → use)

#### A. Needed for almost everything (~14 GB)

| file | size | use |
|---|---|---|
| `/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta` (+ `.fai`) | 3.4 GB | gorilla reference (GGO T2T, NC_ names) |
| `/mnt/linuxdisk/home/juanfraitu/gorilla_haps/{mat,pat}.fa` (+ `.fai`) | 3.4 + 3.2 GB | KB3781's two haplotypes: the answer keys |
| `/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa` (+ `.fai`) | 3.0 GB | human reference (human work only) |
| `/mnt/linuxdisk/tmp/rna_allele/{chrmap.tsv,mat.len.tsv,pat.len.tsv,fibro_cat.copies.tsv}` | < 1 MB | which haplotype contig equals the reference; the KB3781 catalog |
| `/mnt/linuxdisk/tmp/rustle_figures/runs/*/*.fam.copies.tsv` | < 10 MB | O1 catalogs (human_A119b, gorilla_KB3781, ...) |
| `/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_sedef_final.bed` | 103 MB | gorilla SDs (SEDEF) |
| `/mnt/c/Users/jfris/Desktop/CHM13_SDs.SD90.bed`, `CHM13_SDs.SD90.nosat.bed` | 30 + 19 MB | human SDs |
| `/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_genomic.gff` | 662 MB (gzip it) | gorilla RefSeq annotation |

#### B. Reads: copy by task, or a subset

| file | size | advice |
|---|---|---|
| `/mnt/linuxdisk/home/juanfraitu/fibroblasts/GCA_029281585.2_flnc_mm.bam` (KB3781 fibroblasts, all records) | 22 GB | copy a **family-region subset** (recipe below) unless whole-genome work is planned |
| `/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_mm.bam` (gorilla OR6737 testis; the assembly benchmark) | 11 GB | only for the Figure 3 / StringTie-FLAIR work; subset possible |
| `/mnt/linuxdisk/home/juanfraitu/winloci_data/A119b.t2t.bam` (human testis) | 90 GB | **do not copy**; subset to family regions if human work is planned |
| `/mnt/linuxdisk/tmp/rustle_figures/runs/{gorilla_OR6737,human_A119b}/*.molecules.tsv` | 260 MB / 1.4 GB | per-read best / second AS (Figure 3, register 1366-style analyses) |
| `/mnt/c/Users/jfris/Desktop/benchmark_collapse/` | 13 GB | copy only the GTFs (`stringtie_*/`, `flair_*/*.isoforms.gtf`, `assembler_*/`, `ref/`); skip FLAIR's BED / FASTA (7 GB) |
| `/mnt/linuxdisk/tmp/rustle_figures/assembly/*/eval_*/` (reference GTFs, tool GTFs, tmaps) | 7.1 GB | only for rescoring Figure 3; the tables in `figures/data/` are in the repo |

**Family-region subset recipe.** It keeps every record of every read that touches a catalog copy (±50 kb). Secondaries elsewhere
in the genome stay, so the O2 analyses still work:

```bash
awk -F'\t' 'NR>1 {s=$5-50000; if (s<0) s=0; print $4"\t"s"\t"$6+50000}' fibro_cat.copies.tsv | sort -k1,1 -k2,2n > fam.bed
samtools view -M -L fam.bed -F 4 BAM | cut -f1 | sort -u > fam.names
samtools view -@ 4 -N fam.names -b -o subset.bam BAM && samtools index subset.bam
```

#### C. Do not copy (rebuild or regenerate)

- **Splice indexes** (`*.splice.mmi`, 13 GB each): rebuild with `minimap2 -x splice -d X.mmi X.fa`. That took 3.5 min and 23 GB RAM
  here; on a smaller machine, align against the FASTA directly.
- **Work directories:**
  - `/mnt/linuxdisk/tmp/o3_npip/` (17 GB, mostly an index copy);
  - `o3_graph/idealsim/` (2.5 GB);
  - `multicopy_gains/`, `o1_sd/`.
  
  Their results are in the documents above. Each has a frozen `audit_frozen_*` copy if a number needs rechecking.
