# Ledger: O3 (copies the reference lacks)

Compacted 2026-11-01. Each entry keeps the claim, the rule, the verdict as recorded and the deciding numbers; the full text of every
section is in `docs/archive/2026-11/LEDGER_O3.md` (line given per entry). Merged 2026-10-24 from 51 dated docs; originals in
`docs/archive/2026-10/merged/`.

## Current state (2026-11-03)
- **Framing:** O3 = detection of a copy the reference lacks from RNA; placement belongs to DNA (Amendment 48). An absent copy cannot be certified from RNA alone; an extra copy can (haplotypes > allele slots), subject to Theorem A. IsoCon's reference-free strength is its significance test on a pool of reads from one locus; nothing in it discovers a locus. Our port of that test is `candidate_copies` (clusters -> consensus -> significance merge -> flag / link). A silent pseudogene emits no reads: nothing to assemble; secondaries give membership of reference-present copies only, and a reference-absent copy shows as a mismatch pattern fitting no catalog locus (`candidate_copies`' job).
- **`candidate_copies` (shipped stage): OPT-IN** (`--candidates`, ruling R14). Phase 1 over a family pool (r1506: 2.0x, products identical). `--molecules P.molecules.tsv --min-ms-ratio R` (r1513): nets keep only secondaries at ms >= R x the read's primary ms; at R = 0.5 on the 40-family KB3781 panel it found 6 reference-absent contigs on Jim's haplotypes vs 3 for the default pool with fewer junk flags (19 vs 21) but no wall gain (447 vs 416 s): FAIL on the registered wall bar, opt-in; re-register without the speed bar. 14 of 38 junk contigs sit on two recurring repeat loci (CM054591.2:75.61 Mb, CM054595.2:100.75 Mb): a filter candidate.
- **Non-T2T genomes (r1504):** flag + candidates on A119b -> GRCh38 and GGO -> Kamilah find 1 of 32 / 0 of 37 missing copies (precision 0.012 / 0): not usable as a non-T2T missing-copy detector.
- **PSV-cluster O3 (`copy_assign --psv-clusters`):** PASS in simulation; on real reads NOT claimed (r1294). **O3 port (`copy_phasing`, r1402):** v6 PASSES dev and held-out P5; putative tier FAIL P6 / PASS P10.
- **Exon / SD level (Jim):** mat vs pat differ for 24% of families at exon level (r1416), WGS-confirmed core 3.0% SD / 6.3% families (r1419); BPCLUST (r1417) and per-exon certified tier (r1418) FAIL; reads could see 41.9% of gains (r1420): the detector is the limit; trio WGS confirms the NPIP and STON1 haplotype-only copies (r1422).
- **Defensible sentence (from O3_STATUS, 2026-11-01):** O3's RNA detectors are refuted, vacuous or unvalidated; O3 itself is untested in the regime where it would occur. Two DNA-level positive instances exist (GOLGA6L7 = the positive control; `LOC134757045` chr14 mat 5 / pat 3). No shipped binary emits a reference-absent copy at defaults.
- **Second individual on Jim's T2T (r1578/1579, 2026-11-08):** OR6737 testis vs depth-matched Jim: piles +5.7 points as registered (borderline; +2.2 on shared loci), three-haplotype test +2.5 on 7,056 shared loci (claim FAILS); 12 strong loci, 7 in families (NBPF, IGHV-like lead). Tissues not matched; testis-only loci untested.
- **Open:** detector sensitivity on exon gains; allele vs copy on short reads; SPDYE1's hidden copy (HG002); r1513 re-registration; the repeat-locus filter.

## Contents

Merged 2026-10-24 from 51 dated docs, oldest first (originals in `docs/archive/2026-10/merged/`). Entries below follow the original order; each `## <file>.md` heading is the source doc, followed by its own `## ` title. (full: archive/2026-11/LEDGER_O3.md:6)

## O3_CANDIDATES_AMENDMENT15_2026-10-05.md

## O3 candidates, Amendment 15/15b — the consensus-defect fix and its re-runs, 2026-10-05
Consensus-defect fix in `o3_candidates` (commit `a13b817f`: cs normalisation, one copy per long insertion, majority test for >= 20 bp insertions). Default-on flip needed A13-1/2/3 AND C1'/C2' AND H1/H3.
- A13 (53-family deletion held-out): A13-1/2/3 PASS (49 flagged vs 82 uncorrected; D right 7,882; false moves 0.77 %; unions keep 99.33 %).
- A14 no-deletion control: C2' HOLDS (0.54 %), **C1' FAILS** (10/53 = 18.9 %, bar <= 15.1 %). Held-out H: H3 HOLDS (0.26 %), **H1 FAILS** (6/30, bar <= 4/30).
- Verdict: **stage STAYS OPT-IN** (`--candidates`, ruling R14); remaining false flags = few-read consensus. Whole-BAM batch STOPPED at 570 s limit (R23 stands).
- Corrections 2026-10-06: uncorrected D right was 7,898 (fix cost 16 placements); A13-3 time 33.4 min vs 25.9 (+29 %, not ~15 %); GWFAM175_0/331_0 are not few-read flags. Verdicts unchanged. (full: archive/2026-11/LEDGER_O3.md:63)

## O3_MATERNAL_REFERENCE_2026-10-08.md

## What happens to the reads of a copy the reference haplotype lacks, and what a reference-free step recovers (2026-10-08)
KB3781 fibroblast Iso-Seq vs its own `mat` / `pat` assemblies (study of the prereg below).
- Q1: 19 of 20 large "missing" copies have a reference relative within the allele cutoff (id 0.9911-1.0000); only GWFAM175_B0 (0.9126) is a real test: its 304 reads are ABSORBED (285 elsewhere) at de 0.0662 vs 0.0007, a divergence pile, not a hole. mat: REFUTED at 4, SUPPORTED GWFAM205_B0, between bars GWFAM64_B11; pat: REFUTED 11/13, PARTLY GWFAM382_B0, SUPPORTED GWFAM390_B1. No judged copy unmapped.
- Unmapped: 125 of 132 primary-unmapped reads that map on mat sit in one 71-kb mat chr12 locus.
- Q2 (mat ref): IsoCon R1 PASS, R2 PASS, R3 FAIL (6/24, bar 20 %; 4 of 6 are pat sequences the truth table misses); in-house R1 FAIL, R2/R3 PASS. pat ref: R1 not testable, R2/R3 pass both arms.
- LRPAP1 p12 not mat-absent (Amendment 2); extra chr12 locus may be c01, untested. GWFAM382_B0 reads 5-6 % from both genomes (copy in neither). (full: archive/2026-11/LEDGER_O3.md:169)

## PREREG_o3_maternal_reference_2026-10-08.md

## Pre-registration — O3: what happens to the reads of copies the maternal haplotype lacks, and can a truth-free chain recover them?
Q1 fate (UNMAPPED / TIED / ABSORBED) of pat-only copies' reads on `mat`; Q2 can net -> cluster -> flag (IsoCon, in-house `o3_candidates`) recover them. Bars per LARGE locus: >= 80 % ABSORBED, <= 10 % UNMAPPED/TIED => REFUTED; 20-50 % => PARTLY; >= 50 % => SUPPORTED. R1 recover every beyond-delta LARGE locus; R2 no within-delta flag; R3 false-flag families <= 20 %. Prediction: REFUTED everywhere. Amendment 1 adds the symmetric `pat` run (held-out) and side-by-side readouts; Amendment 2: p12 is descriptive (`lrpap1_desc`), arm H includes the Amendment 15/15b fix (acceptance not run). Results: entry above. (full: archive/2026-11/LEDGER_O3.md:319)

## PREREG_unmapped_rescue_2026-10-08.md

## Pre-registration — rescuing reads that became unmapped: a pool-first, reference-free net
Pool unmapped reads first (all-vs-all, de <= 0.00958, >= 50 % block, clusters >= 3), one abPOA consensus per cluster, attribute the cluster (not the read) to a family, rescue its reads. Beds A (dev excision), H (held-out A13 53 families, run once), M (959 real unmapped KB3781 reads). Bars: U1 coverage >= 0.90, purity >= 0.99 on H; U2 rescued-correct >= 2,656 with wrong <= 5 % AND pooled > B1 > B0; U3 >= 90 % of the 125 chr12 reads in one cluster, id >= 0.99. 48 amendments follow (2026-10-08/09): margin cover score (A1), HMM profile (A5), hybrid net (A6, not met), augmented reference (A7, bars met), LRPAP1 (A8), consensus correction (A9, refuted), synthetic erase world (A10-12, divergence move rule adopted), 5' G run (A13-17, gate-aware metric adopted, trimming not), PART partition (A18-25, opt-in), chain-aware clustering (A22-24, opt-in), three-way allele/copy call (A26-27), discovery from the 959 reads + DNA k-mer check (A28-33), assembly-truth recall (A41-47: corrected-read support recall 0.566 on paternal truth, refused on precision in A46). **Amendment 48 (user):** O3 = detection, placement belongs to DNA; `o3_export.py` writes 2,559 KB3781 candidates (DIVERGED 2,357, NOVEL 202), `placement "requires_WGS"`. (full: archive/2026-11/LEDGER_O3.md:487)

## UNMAPPED_RESCUE_2026-10-08.md

## Rescuing reads that became unmapped: a pool-first, reference-free net (2026-10-08)
- U1 **PASS** on H (95.7 % clustered, purity 1.0000). U2 **FAIL as a conjunction**: first clause passes (4,255 = 80.1 % rescued, 3 wrong = 0.07 %, margin 1.10), but pooled does not beat single-read dc-megablast B1 (406 vs 410 on H). U3 **PASS** (92 %, consensus id 1.0000; attribution unstable).
- Margin 1.10 is fragile (sits just above a few wrong clusters). Fidelity id x cov >= 0.999: 14 of 40 on H.
- Exploratory: HMM profile +243 reads (hypothesis only); hybrid net not met; augmented reference 94.3 % land on own consensus, 0 of 41,727 survivors move; consensus correction refuted (15 vs 14 of 40, bar 19); divergence move rule adopted; G-run trimming not adopted, gate-aware metric adopted; PART and chain-aware clustering adopted opt-in; DNA k-mers confirm the 5.7 kb chr12 sequence, the 4- and 3-read clusters are RNA-only. (full: archive/2026-11/LEDGER_O3.md:1576)

## PHASED_TOY_2026-10-12.md

## Phased-data toy — results (2026-10-12)
30-kb CHM13 chr20 SNRPB toy: a collapsed reference sees the SAME 5 columns for an allele (SA, alt 0.50) and an extra copy (SB, 0.53); mapping to both haplotypes separately separates them. P0-P3 all **PASS** (SA 1.00 ALLELE, SB 1.00 EXTRA_COPY, 150/150 reads on true haplotype). Easy case by design (dispersed, 0.43 %, perfect assemblies). (full: archive/2026-11/LEDGER_O3.md:1766)

## PREREG_phased_toy_2026-10-12.md

## PREREG — phased-data toy (2026-10-12)
Collaborator question "would phased data help?": scenarios S0 homozygous, SA allele, SB extra copy on hap2; arms collapsed vs phased (truth-free NM rule + flank alignment). Bars P0 <= 2 % non-SHARED; P1 identical columns, alt within 0.10; P2 >= 95 % correct class; P3 >= 98 % true haplotype. Results above. (full: archive/2026-11/LEDGER_O3.md:1799)

## PREREG_missing_copy_posctl_2026-10-14.md

## PREREG — `missing_copy_flag` positive control, re-measured with stable seeds (2026-10-14)
Re-measure the "40/40 at 2 %" positive control (old `hash()`-seeded sim) with `sim.py missing-copy genomic`, 5 seeds x 40 genes, CHM13 chr20, plus a template-only negative control. Bars: P1a 2 % fired >= 190/200, P1b/P1c >= 0.95, P2a 1 % fired in [50,100], P2b >= 0.95, N <= 2/200. (full: archive/2026-11/LEDGER_O3.md:1851)

## OUTCOME (2026-10-14) — every bar passes; the control reproduces on stable seeds
All **PASS**: 2 % 199/200 fired, 198 RAC, 197 confirmed; 1 % 91/200 fired, 88/90 confirmed; N 0/200. Misses all at lncRNA / LOC genes. `REPRODUCE.md` quotes the pooled numbers. (full: archive/2026-11/LEDGER_O3.md:1947)

## PREREG_o3_candidates_foreign_check_2026-10-14.md

## PREREG — are `candidate_copies`' remaining false flags other apes' transcripts? (2026-10-14; to run on the main machine)
Follows the F < 5 branch below: align the 26 A14/H false-flag unions to own gorilla assemblies and to CHM13, chimp, orangutan, siamang; FOREIGN iff foreign >= 0.995 and foreign - own >= 0.00479. Reading: >= half FOREIGN => M3 (other apes' transcripts) is the main cause; <= 5 => not. Self-check: class-a flags must be OWN >= 0.999. (full: archive/2026-11/LEDGER_O3.md:1988)

## OUTCOME
Recorded as "(to be appended after the run)"; no result in this ledger. (full: archive/2026-11/LEDGER_O3.md:2049)

## PREREG_o3_fewread_flags_2026-10-14.md

## PREREG — `candidate_copies` few-read false flags: does a clean simulation produce them, and by which mechanism? (2026-10-14)
Do M1 (single-linkage mixing) or M2 (few-read consensus error) produce the A14/H false flags in a clean famsim world (48 families, seeds 1-2, E0/E3)? Branch: F >= 5 of 96 NULL units => develop a rule (Phase 2, seeds 3-5); F < 5 => go to the M3 foreign check. (full: archive/2026-11/LEDGER_O3.md:2056)

## OUTCOME (2026-10-14) — F = 0: the clean simulation does not reproduce the few-read false flags
**F = 0 of 96.** ERASED recall 48/48 at 6-80 reads (0/8 at 4 reads, below `--min-support 6`); true flags at <= 20 reads have id x cov >= 0.9978, so a read floor would only cost real copies. One M1 flag (131-read allele+copy chain). Next: the M3 foreign check. (full: archive/2026-11/LEDGER_O3.md:2140)

## PREREG_o3_structure_projection_2026-10-14.md

## PREREG — O3 structure projection: held-out simulation test (2026-10-14)
Claim: `candidate_copies` gives a reference-absent copy its exon structure within +-5 bp. P0 parity gate (Rust vs prototype) **PASS** 148/148 (one naming defect fixed). Bars on seeds 21-25: S1 >= 0.90, S2 >= 0.80 per structural scenario, S3 boundary recall/precision >= 0.97. (full: archive/2026-11/LEDGER_O3.md:2173)

## OUTCOME (2026-10-14) — S1, S2, S3 all PASS on the held-out worlds
S1 210/222 = 0.946; S2 all >= 0.80 (SHUFFLE 16/20 exactly at bar); S3 0.995 / 0.995. Limit: shuffled exons on a copy opposite its home paralog. Reconstruction: candidate_copies 242/270 vs IsoCon 256/270 (difference = detection at <= 1 %). `<out>.structure.tsv` / `.contigs.gtf` written by default. (full: archive/2026-11/LEDGER_O3.md:2253)

## PREREG_o3_family_best_2026-10-15.md

## O3: a read counts only for the family holding its best alignment (2026-10-15)
`copy_assign --psv-family-best` (opt-in, commit `7c9d2e8b`): a read enters family F's PSV clustering only if F holds one of its best alignments. Fresh held-out heldout3 (40 families). Bars S1 <= 0.10 families with a `copy` cluster; S1b <= 0.5 x old; S2 recovered >= old - 1; G guard.
- Result: **Reading 3: not claimed.** S1 FAIL (6/40 = 0.15), S1b FAIL (10 -> 6), S2 FAIL (7 -> 5 of 24; MCL117, MCL31 lost: their reads' best alignment is an uncatalogued locus), G pass.
- Post hoc: adding isoform content explains 1 of 136 clusters; a tolerant rule (median best-alignment `de` > delta) cuts families with a `copy` cluster 10 -> 2 on heldout3 (hypothesis only, seen substrates). (full: archive/2026-11/LEDGER_O3.md:2290)

## PREREG_o3_graph_real_2026-10-15.md

## Graph-based O2 and O3 on handpicked real families (simulation v2, 2026-10-15)
Simulated reads of handpicked CHM13 mainstay families (NPIPA/B, TBC1D3, GOLGA8/6L, CT47A, XAGE, SPATA31). Path O2 = reads to every curated isoform path; S4 = O3 on O2's `unexplained` reads + genome check. Bars O2a-d, O3a-d.
- First verdict (superseded): O2 reading 1 claimed (+14.3 %, precision 1.0); O3 BRANCH 8/8 (current 0), 0/32 false calls.
- Amendment 1 CROSS world (paralog-annotated structures): X1-X3 **fail** (path O2 precision 0.946, 115/273 correct), X4 pass.
- **Verdict, revised after the review:** O2 claimed only for a complete, copy-specific catalog, not real annotation (genome-separated gain +8.4 % < 10 %). O3 = branch calls + genome check, not O2-informed O3; REMOVED not measurable (6/8 within delta). GOLGA8 mosaic = next step. (full: archive/2026-11/LEDGER_O3.md:2562)

## PREREG_o3_graph_round4_2026-10-15.md

## Graph O2 as assign-or-abstain, and family-pooled mosaic-aware O3, on a fifth handpicked family set (2026-10-15)
Round 4: padded-union graph O2 (assign or abstain) and family-pooled O3 with the mosaic rule, held-out 8 families, seed 51.
- **Verdict, revised after the review (2026-10-16): O3 claimed (reading 2; B1-B4 pass):** novel-exon absent copies 7 of 8, detectable removed 6 of 6, false calls 0 of 24. Gains attributed to branch calls (novel exons), genome check (false calls), mosaic rule (PRAMEF); pooling added nothing. Not for sister-family hosting (dev set 2: false calls in 4 of 5 worlds).
- **O2 not claimed** (A1 fails, 0.969 in CROSS); genome-separable gain over current O2 at most +8.4 % across five sets, about 0 since. (full: archive/2026-11/LEDGER_O3.md:2950)

## PREREG_o3_graph_sim_2026-10-15.md

## O3 with a variation graph: three steps, tested on curated ideal scenarios (2026-10-15)
Three steps fed into `--psv-clusters` (S1 PSVs from every isoform, S2 branch reads off the loci graph, S3 genome check), famsim worlds from curated CAT isoforms (commit `883e54bc`).
- **Verdict: Reading 1 claimed, in simulation** (H1-H4 pass): BRANCH 12/12 vs 0/12 current; OUTSIDE false calls 20 -> 0, TWIN_ALT 9 -> 0; divergent absent copies kept (12/12 at 2 %, 8/12 at 1 %). Step 1 redundant with step 3 in these worlds; step 3 rejected a correct 1 % call on both sets. Register 1297. (full: archive/2026-11/LEDGER_O3.md:3253)

## PREREG_o3_graph_union2_2026-10-15.md

## Union-path graph O2 with padded path ends, and family-pooled O3, on a fourth handpicked family set (2026-10-15)
Round 3: padded path ends for union O2, family-pooled O3.
- **Verdict, revised after the review: O2 not claimed (V2, V3 fail)** (no gain, precision 1.0). **O3 not claimed (H3 fails, 2 of 4)**; branch calls 8/8, no false call after genome check; the skipped mosaic rule would have given 3 of 4 (post hoc, round 4's test). (full: archive/2026-11/LEDGER_O3.md:3448)

## PREREG_o3_graph_union_2026-10-15.md

## Union-path graph O2 and union-O2-informed O3 on a fresh set of handpicked families (2026-10-15)
Round 2: union-path O2 (every family structure threaded through every copy).
- **Verdict, revised after the review: O2 not claimed (U4 fails)**; +7.4 % at precision 0.9993 carried by FCGR's one non-canonical intron (x1.025 without). **O3 not claimed (G3 fails, 2 of 3 detectable)**; branch calls 8/8, no false call over 24 families. (full: archive/2026-11/LEDGER_O3.md:3733)

## PREREG_o3_psv_clusters_v2_2026-10-15.md

## PREREG — O3 PSV clustering v2 (linkage test) in Rust (`copy_assign --psv-clusters`), held-out on base and hard worlds (2026-10-15)
v2 adds a linkage test (a missing sequence's alts travel together; editing and errors do not) after v1 broke on A-to-I editing and 1 % error; hard worlds (allele background, conversion, 8-copy, 5x error, imbalance). Bars P0' parity, H1 >= 0.90, H2/H2b <= 3, H3 >= 0.80, H4 >= 0.90. (full: archive/2026-11/LEDGER_O3.md:4012)

## OUTCOME (2026-10-15) — all bars PASS; the Rust product scored
All **PASS**: parity 16/16 runs; H1 425/426 = 0.998 (existing detector 414/426); H2 0/152; H2b 0; H3 41/44 = 0.932; H4 135/135, 262/262. Copies at 0.3 % 24/24. `copy_assign --families --no-as-tied-only --psv-clusters` is O3's method at simulation level; next = real families. (full: archive/2026-11/LEDGER_O3.md:4078)

## PREREG_o3_realdata_2026-10-15.md

## O3 on real reads: copies removed from the genome, then added back (2026-10-15)
KB3781 real reads, copy masked from the genome and re-mapped, then added back; 60 held-out families + A13 dev. Bars R1 recovered >= 0.70; R2 cleared >= 0.90; R3 false `copy` <= 0.02 of 300 single-copy genes.
- **Verdict: Reading 3, not claimed.** R1 FAIL 4/12 = 0.333 (dev 15/25); R2 FAIL 3/4; R3 not measurable (one-copy families get no placed reads, so the 48 two-copy families and the control are untested). Register 1294.
- Follow-up (register 1295): of 57 held-out `copy` clusters in intact families, 35 are reads of another family's copy placed by `--no-as-tied-only`; other apes' transcripts remain; only GWFAM175 (dev) is a real reference-absent copy. (full: archive/2026-11/LEDGER_O3.md:4113)

## PREREG_o3_unexplained_2026-10-15.md

## PREREG — O3 (advisor's definition): PSV clustering of the reads O2 places on a family, held-out simulation (2026-10-15)
Cluster each copy's O2-placed reads by read-derived variant sites (`psvclust.py`); a cluster no copy explains = missing sequence. Bars D1 >= 0.90, D2 <= 2/120, D3 >= 0.80, D4 >= 0.90. (full: archive/2026-11/LEDGER_O3.md:4369)

## OUTCOME (2026-10-15) — D1-D4 all PASS on the held-out worlds
All **PASS**: D1 310/310 (existing 292/310); D2 0/120; D3 40/40; D4 74/78 = 0.949, 170/170. Arm T (exact ties only) recovers only twins (D1 80/310): O3 must cluster every placed read, not only ties. (full: archive/2026-11/LEDGER_O3.md:4435)

## PREREG_o3_graph_round10_2026-10-16.md

## Round 10: the integrated pipeline v3 under calibrated v2 reads, end to end, on an eleventh set of full families (2026-10-16)
Integrated O1 -> O2 -> O3 v3 vs the current method, 6 world types x 2 calibrated read profiles (A119b, public). Bars per world: B1 integrated correct >= current; B2 errors <= current; PASS needs all 12.

## Results (2026-10-16; run with the code of 58bca61d, started 19:30)
**Verdict: FAIL under both profiles. Not claimed** (register 1308). B1 holds in 11 of 12 (+2 to +370); A119b B2 fails in INTACT/CROSS/SILENT (8 errors vs 0: one false O3 call in the 45-copy TAF11L array, a mixed-cluster host-window artefact per review); public OUTSIDE B1 fails (-51: O1's new FAM236.1 within delta of FAM236.0). Catalogued-read loss of round 9 is gone. (full: archive/2026-11/LEDGER_O3.md:4466)

## PREREG_o3_graph_round11_2026-10-16.md

## Round 11: the integrated pipeline v4 under calibrated v2 reads, end to end, on a twelfth set of full families (2026-10-16)
v4 removes round 10's failure mechanisms; same B1/B2 bars over 12 worlds.

## Results (2026-10-16; run with the code of fe73bc8c, started 20:57)
**Verdict: PASS under both profiles. Claimed** (register 1309); B1 and B2 hold in all 12 worlds, also under round 10's rule and copy level. Review: PASS holds, but margins are the current method's own O3 artefact in chr6 histone clusters, not an O2 gain (O3 removed: A119b at parity); `--o3-procs 2` non-deterministic (fixed before round 12); zero errors rest on the call-level genome check (~0.005 margin). (full: archive/2026-11/LEDGER_O3.md:4765)

## PREREG_o3_graph_round12_2026-10-16.md

## Round 12: the integrated pipeline v5, independent of the current method, on a thirteenth set of full families (2026-10-16)
v5 hosts O3's reads from O1 (no dependence on copy_assign); fresh families incl. an 80-copy DUX4 array.

## Results (2026-10-16; run with the code of 00c241db, started 22:29)
**Verdict: FAIL under both profiles. Not claimed.** B1 holds in all 12 (+191 to +713); B2 fails in 4 worlds (A119b INTACT 3, CROSS 5, SILENT 2; public CROSS 2 vs 0): O2 margin-2 assignments to a near-identical wrong copy (DUX4, KRTAP2, CDRT15), never O3. Review: the FAIL stands; mechanism was wrong for 5 of its 7 distinct reads. Not cheaper on an 80-copy array (1,095 s vs 383 s). (full: archive/2026-11/LEDGER_O3.md:5074)

## PREREG_o3_graph_round5_2026-10-16.md

## Graph O2 with a structure-free fallback, and cross-family-guarded O3, on a sixth handpicked family set (2026-10-16)
Round 5: graph O2 with a structure-free (locus) fallback; O3 with a cross-family guard (K).
- **Verdict, revised after the review: O3 claimed without the removed-copy clause (reading 2):** novel-exon absent copies 8 of 8, false calls 0 of 24 (from branch calls + genome check; K's own parts untested). **O2 not claimed (C3, C4 fail):** +7.6-9.4 % correct but precision 0.834 counting uncatalogued-copy reads vs current 0.90-0.91. Post hoc: no per-read error cutoff separates near-identical uncatalogued copies. (full: archive/2026-11/LEDGER_O3.md:5348)

## PREREG_o3_graph_round6_2026-10-16.md

## O2 + an unannotated-copy flag, with heterozygous alleles, on a seventh handpicked family set (2026-10-16)
User question: is O2 misplacing reads because it detects unannotated copies? Flag tested on a fresh set with alleles.
- **Verdict, revised after the review: Not claimed (reading 2: H2 and H3 fail; H1 also fails).** On the fresh set O2's misplacements have no PSV signature and all 13 intact flags are heterozygous-allele clusters; the flag needs allele information (DNA or phasing). (full: archive/2026-11/LEDGER_O3.md:5651)

## PREREG_o3_graph_round7_2026-10-16.md

## Integrated O1 -> O2 -> O3 on an eighth handpicked family set, with alleles (2026-10-16)
Does O1 seeding (only RNA-supported copies) shrink the universe and help O2/O3?
- **Verdict: O2 not claimed (J2, J4 fail)**: +5-22 % correct, 72 % of uncatalogued reads recovered, but more wrong calls where uncatalogued near-identical copies were not seeded or seeded into two families (WASH/DDX11L), and only 71 % of silent copies dropped. **O3 not claimed (J7 fails, 2 of 3)**; no false call over 32 families. Bug: a seeded locus must belong to one family. Post hoc end-to-end fates favour the integrated pipeline (basis for round 8). (full: archive/2026-11/LEDGER_O3.md:5873)

## PREREG_o3_graph_round8_2026-10-16.md

## Round 8: the integrated O1 -> O2 -> O3, end to end, on a ninth handpicked family set (2026-10-16)
End-to-end bars per world: B1 correct >= current, B2 errors <= current.

## Results (2026-10-16; worlds built 13:22-13:31 with the code of e4d0703b)
**Verdict: FAIL. B2 fails in CROSS by one read. Not claimed.** B1 holds in all 6 worlds (x1.16-x1.51); the one error is TSPY.2 read s0001190 called TSPY.1 (a sequencing error at the only separating base + soft-clipped bases counted as edits). O3 no false call. (full: archive/2026-11/LEDGER_O3.md:6125)

## PREREG_o3_graph_round9_2026-10-16.md

## Round 9: the integrated pipeline v2 under realistic reads, end to end, on a tenth set of full families (2026-10-16)
v2 with realistic (A119b, public) read profiles; B1/B2 per world, 12 worlds.

## Results (2026-10-16; run after amendment 1 with the code of 8235f733)
**Verdict: FAIL under both profiles. Not claimed.** Wins where the catalog is incomplete (OUTSIDE +201/+204, BRANCH +120/+114); loses 1-3.5 % correct on catalogued copies (O2's divergence gate leaves realistic reads `unexplained`); public PGA errors. Review: verdict holds. (full: archive/2026-11/LEDGER_O3.md:6402)

## PREREG_o3_graph_round13_2026-10-17.md

## Round 13: the integrated pipeline v5', independent of the current method, on a fourteenth set of full families (2026-10-17)
v5' fixes round 12's structural (retained-intron) errors and the branch-arm cost; B1/B2 over 12 worlds.

## Results (2026-10-17)
**Verdict: FAIL (A119b). Not claimed** (register 1311). public passes all 6; A119b B2 fails INTACT 3 vs 0, CROSS 3 vs 2, SILENT 3 vs 0: one false O3 mosaic call on three short MRGPRX.1 haplotype-2 reads (review: driven by one 506-bp read's indels, not the alleles alone). Without OR7E, A119b's B1 margins turn negative. Branch prefilter holds cost (223 s vs 1,095 s). (full: archive/2026-11/LEDGER_O3.md:6660)

## PREREG_o3_graph_round14_2026-10-17.md

## Round 14: the integrated pipeline v9 on a fifteenth set of full families, with a diploid donor (2026-10-17)
v9 adds one rule per round 12/13 mechanism (incl. `--haplotype-check`), diploid donor under both profiles.

## Results (2026-10-17)
**Verdict: PASS under both profiles (12 of 12 worlds). Claimed** (B1 +52 to +298; B2 holds, 0 errors in INTACT/CROSS/SILENT/REMOVED). Review qualifiers: one draw (at 4 further seeds A119b fails 2); secondary scorings fail A119b REMOVED; public CROSS rests on the haplotype check; veto and shared-position mosaic not exercised; **the allele confound is not solved**; family set fixed 7.5 h before the method was frozen (undisclosed). (full: archive/2026-11/LEDGER_O3.md:6948)

## PREREG_o3_graph_round15_2026-10-17.md

## Round 15: integrated v11 on fresh families, three seeds per profile, and the cost of full families (2026-10-17)
6 cells (2 profiles x 3 seeds) x 6 worlds, B1/B2 per world, PASS needs all 72 checks; C1 integrated wall <= current on the mainstay set. Results: **Verdict: FAIL** (3 of 6 cells; every failing bar a few FAM90A errors from the fixed 2-edit margin among 57 rivals). **C1 FAILS** (416/418 s vs 320/333 s), though 2.4-2.9x faster than v9. (full: archive/2026-11/LEDGER_O3.md:7223)

## PREREG_o3_graph_round16_2026-10-17.md

## Round 16: integrated v12 (the O2 margin) on every unused family of round 15's pool, and the cost of full families in CPU seconds (2026-10-17)
Results: **Verdict: FAIL** (4 of 6 cells; A119b B1 short by 4-7 reads in INTACT/SILENT, small families where the current certificate wins). B2 holds in all 36. **C1 PASSES** (CPU 810/660 vs 1,046/949). (full: archive/2026-11/LEDGER_O3.md:7404)

## PREREG_o3_graph_round17_2026-10-17.md

## Round 17: integrated v13 on name-free families (linked by sequence identity and gene structure), three seeds per profile, and the cost of full families in CPU seconds (2026-10-17)
Results: **Verdict: FAIL, 0 of 6 cells. C1 PASSES.** B1 holds everywhere (124,264 vs 82,008 correct); OUTSIDE B2 fails 3-5x: an uncatalogued twin O1 did not seed has no rival in O2 (register 1314-1315's open failure). (full: archive/2026-11/LEDGER_O3.md:7589)

## PREREG_o3_graph_round18_2026-10-18.md

## Round 18: integrated v14 (copies the catalogue lacks, found in the genome) on fresh name-free families, three seeds per profile, and the cost in CPU seconds (2026-10-18)
Results: **Verdict: FAIL, 4 of 6 cells pass** (seed 803 B1: OR4F21 family plus five silent pseudogene copies O1's genome search adds; the current method's `--no-as-tied-only` ignores them). B2 holds in all 36; round 17's OUTSIDE failure gone. **C1 PASSES** (-3 % / -13 %). (full: archive/2026-11/LEDGER_O3.md:7780)

## PREREG_o3_graph_round19_2026-10-18.md

## Round 19: integrated v14 (unchanged) on fresh families whose non-coding members are copies like any other, expressed as real testis expresses them (2026-10-18)
Results: **Verdict: PASS, 6 of 6 cells, all 72 bar checks. C1 PASSES.** 85,440 vs 62,070 correct, 2,355 vs 8,028 errors; thinnest B1 margin +29. Claims nothing on real reads, FAM90A/DUX4-size arrays, or unannotated paralogs. (full: archive/2026-11/LEDGER_O3.md:7983)

## PREREG_o3_graph_round20_2026-10-18.md

## Round 20: integrated v15 (O1's genome search under the carried intron chain, every chain kept) on fresh families built by the same rule (2026-10-18)
Results: **Verdict: FAIL, 5 of 6 cells** (A119b 1001 B2: three ZNF479 reads moved to ZNF734P by the genome veto, which trusts a MAPQ >= 30 primary despite an AS-tied secondary). **C1 PASSES.** (full: archive/2026-11/LEDGER_O3.md:8173)

## PREREG_o3_graph_round21_2026-10-19.md

## Round 21: O2 only on ambiguous reads, by Eichler's rule (within 10 AS), on fresh multigene families (2026-10-19)
User rule: O2 only on reads with a rival within 10 AS; others are minimap2's. Results (register 1331): **Verdict: FAIL**, 5 of 6 cells, 69 of 72 checks (one NCOR1P4 read: the wrong copy's alignment clips 38 bases for free and outscores by 16 AS). C1 PASS. On AS-tied reads O2 placed 349, all right, abstained 6,651; minimap2 wrong on 4 of 49,914. Not adopted as registered. (full: archive/2026-11/LEDGER_O3.md:8355)

## PREREG_o3_graph_round22_2026-10-19.md

## Round 22: O2 only on ambiguous reads, with sequence ties, on fresh multigene families (2026-10-19)
Ambiguous = within 10 AS OR another locus explains the read with no more edits counting unaligned bases (`veto_distance`). Amendment 1: disk-full crash, five cells rerun once. Results: **Correctness: PASS**, 6 of 6 cells, 72 checks (78,421 vs 53,763 correct; 3,255 vs 12,600 errors). **C1: incomplete** (A119b passes; public cost run stopped by the user, UNDECIDED). Default adoption needs the public cost bar. (full: archive/2026-11/LEDGER_O3.md:8596)

## PREREG_o3_npip_candidates_realign_2026-10-21.md

## Pre-registration: do the maternal-only NPIP locus's reads find their copy once the O3 candidates join the reference? (2026-10-21)
Add register 1363's two candidates to the reference and realign. Bars H1 >= half of unit C's reads on C's candidate; H2 same for B; H3 <= 1 % false moves.
- Results: **PASS** (H1 33/39, H2 21/23, H3 0/2,861); AS-tied target reads 34 -> 10. Audit: PASS holds with caveats (H1 robust, 26/32 on non-IsoCon reads; H2 mostly self-fulfilling).
- Amendment 1 (candidates' relatives): **H3′ holds**, 1 of 8,339 moved and it is a target read; false moves 0. (full: archive/2026-11/LEDGER_O3.md:8777)

## PREREG_o3_npip_maternal_copy_2026-10-21.md

## Pre-registration: does the reference-free chain name KB3781's maternal-only NPIP locus as exactly one new copy, blind to his assemblies? (2026-10-21)
Positive control of the frozen chain (IsoCon, link at 0.96 %, flag). Bars E1 found, E2 exactly one, E3 holds >= half of R_T, E4 >= 90 % pure, E5 no false flag.
- Results: **FAIL (E3, E4, E5)**: one haplotype-only candidate in T (E1, E2 hold by one base), holding 7 of 94 target reads; IsoCon left 66 unassigned. Audit: FAIL holds; T holds three transcription units (C, B, A), each unit with a candidate got exactly one; E5's false flag is poly(A) consensus noise; E3 was out of reach by design. Detection works; read recovery depends on IsoCon converging. (full: archive/2026-11/LEDGER_O3.md:8985)

## O3_FINGERPRINT_SURVEY_2026-10-22.md

## O3 as an RNA-only "look harder" screen: synthesis of the vault survey (2026-10-22)
Design only (no held-out score). Advisor framing: O3 = a screening rank ("suggest when to look harder"), RNA vs one reference. Two bounds: an **extra** copy can be certified from RNA (haplotype excess over 2 allele slots per tied copy, bipartite matching); an **absent** copy cannot (statement 2), except by anchored ownership. Four tiers proposed (1 certified extra copy, 2 novel recurrent haplotype, 2b unexplained indel, 3 consider absent). Dead list: MEC as count, four-gamete, depth, soft clips, raw NM/`de`, reference-PSV-only O3. Quick checks: tied-read-only haplotype excess 0/131 HG00733, 0/50 HG002 (specific but blind; needs teacher alleles). (full: archive/2026-11/LEDGER_O3.md:9143)

## O3_SCREEN_SPIKE_2026-10-22.md

## O3 RNA-only "look harder" screen, tiers 1-2: a 45-minute spike on HG002 and HG00733 (2026-10-22)
`o3screen.py` on frozen O2-VC runs. v1 as specified: neither tier beats today's flag; tier 1 catches 2/5 HG002 positives but false-alarms 9/73 HG002 and 3/174 HG00733 (fragmentation + teacher leakage); tier 2 0/5. Post hoc v4 (merge + leak closure): specific (0/174, 2/73) but blind (0/5). AUC-like scores near chance (0.55/0.61 vs today's 0.56/0.64). USP17L11, POTEM, BMS1P1 have no taught copy. Absent copies invisible (0/46). Nothing claimable. (full: archive/2026-11/LEDGER_O3.md:9242)

## LANDING_LOCI_NEARBEST_2026-10-23.md

## Landing loci with near-best support, P2 control (2026-10-23)
P2 control (53 KB3781 families): 746 landing loci, 722 secondary-only (silent paralogs/pseudogenes, O1 members, not O3 copies); 13 near-best: 7 other families' copies (O1 boundary), 2 strong expressed loci missing from copy tables (GWFAM123, GWFAM99), 4 weak, 0 unannotated. Nothing is a copy absent from the reference. (full: archive/2026-11/LEDGER_O3.md:9370)

## O3_FIRST_PILOT_2026-10-23.md

## O1 -> O3 -> O2 pilot on HG002, with Strainy and the advisor's silent-copy check (2026-10-22/23)
HG002, 79 families. Family level: Strainy flags 1/5 (SPDYE1) with 7/73 false alarms; `candidate_copies` flagged 1/5, 14/73. Advisor's claim: 179 of 236 silent copies (76 %) are reached by family multimappers; divergence, not expression, limits it. Label check: only SPDYE1 is a solid extra-copy positive (6.4 kb donor insertion, 99 reads); USP17L11/BMS1P1 labels were liftover artifacts, so the family tables cannot rank arms. Read level: no arm finds SPDYE1's hidden copy (9 % of reads on SPDYE1.18; separable only by linkage); needs IsoCon-style significance + Strainy-style linkage + consensus. (full: archive/2026-11/LEDGER_O3.md:9410)

## O3_ONE_COPY_SIM_2026-10-23.md

## O3 on one erased copy, simulated (2026-10-23)
GWFAM248:3 erased (97.4 % to survivors), 100 reads per allele. The port pulls all erased-copy reads into one pure 184-read haplotype; **certified** only when survivors are heterozygous (contig 99.92 %), else **putative** (Theorem A: a homozygous survivor has a free slot). Silent pseudogenes collect secondary-only piles (516 reads). Sweep (v6, 40 runs): 0 false calls in every control incl. another individual; limits ~5-10 reads per allele and ~0.3 % divergence; found the 6-read-floor bug (P1 5,762 -> 6,129). Correlation clustering (`--clusterer cc`) keeps another individual's two alleles together (certified at 5 per allele). (full: archive/2026-11/LEDGER_O3.md:9512)

## O3_PORT_RESULTS_2026-10-23.md

## O3 port (Strainy's phasing + IsoCon's statistics): results on gorilla KB3781 (2026-10-23)
Port (`o3port.py`, frozen v4 at 1b4535cc): piles within 10 AS, significant informative columns, zero-disagreement clusters, DADA2 chimera rule, allele-slot b-matching certificate (>= 2 private alleles, >= 6 reads), plus IsoCon side.
- **Verdict: not adopted ("detects but not sensitive").** P1 **FAILS** (D right 3,995 vs bar 5,842; K 7,882; false moves 0.02 %); P2 **PASSES** (0 vs K 10 families); P3 H1 **PASSES** (3 vs K 6; one real maternal-only GWFAM390 copy).
- Lessons: RNA piles full of 2-7-read artifact clusters; the certificate is honest but conservative (Theorem A: GWFAM175's paternal-only copy lost to a free slot); losses split between no pile (2,265 D reads) and free slots (1,343). SPDYE1's hidden copy not certified (looks like a chimera). (full: archive/2026-11/LEDGER_O3.md:9637)

## PREREG_o3_port_2026-10-23.md

## Pre-registration: O3 port (Strainy's phasing + IsoCon's test), gorilla KB3781 (2026-10-23)
Arms P (port) vs K (`candidate_copies`, Amendment 15) on the stage's own bars: P1 D right >= 5,842 and false moves <= 5 %; P2 false-flag families <= 8/53; P3 held-out <= 4/30 and <= 5 %. Adopt iff all pass. Verdict above (not adopted). Development after the verdict: v5 fails P1 by 80 reads (5,762); v5.1 `--absorb` rejected; **v6** (6-read floor before the clique count) P1 6,129 **PASSES**, P2 0/53; correlation clustering P1 5,949, P2 1/53 (made default 6e7ed568 after P5); putative tier added (tie, then non-canonical-splicing signal); O2 v2 (assign only when forced over all piles). (full: archive/2026-11/LEDGER_O3.md:9763)

## PREREG_o3_port_v6_heldout_2026-10-23.md

## Pre-registration: the O3 port v6 on a fresh gorilla fibroblast held-out (P5), 2026-10-23
Bars H5-1 D right among attainable >= 0.50; H5-2 false moves <= 5 %; H5-3 false-flag families <= 3/20. Result: **the O3 port v6 PASSES the held-out** (label propagation 0.660, 0 false moves, 2/20; correlation clustering 0.829, 0, 1/20). Main loss is "no pile" (6 of 20 copies), not Theorem A. Choosing cc as default was made after seeing P5. (full: archive/2026-11/LEDGER_O3.md:9910)

## PREREG_o3_putative_clip_p10_2026-10-23.md

## Pre-registration: O3 putative copies with the soft-clip signal on a fresh gorilla fibroblast held-out (P10), 2026-10-23
Adds `foreign_clipping` (>= half of a slot holder's near-best reads clipped >= 50 bp, >= 6 at one boundary). Bars as P6. Result: **the putative tier PASSES** (HP1 6/6, HP2 0 control families, HP3 5/8 = 0.625); every flag carries the clip signature. (full: archive/2026-11/LEDGER_O3.md:10009)

## PREREG_o3_putative_heldout_2026-10-23.md

## Pre-registration: O3 putative copies on a fresh gorilla fibroblast held-out (P6), 2026-10-23
Putative tier (tie or non-canonical splicing signal). Bars HP1 precision >= 0.50, HP2 <= 12/80, HP3 recall >= 0.25. Result: **DOES NOT PASS**: HP1 2/2, HP2 1, HP3 **FAIL** 2/13 = 0.15 (two-copy families: only splicing can fire, and 11 of 13 splice canonically). Superseded by the P10 clip signal above. (full: archive/2026-11/LEDGER_O3.md:10066)

## PREREG_o3_two_tier_p4_2026-10-23.md

## Pre-registration: O3's two-tier output on a fresh gorilla fibroblast held-out (P4), 2026-10-23
Tier 1 certified (port v5) + tier 2 look harder (stage's other flags); bars H4-1 <= 9/60, H4-2 <= 5 %, H4-3 >= 0.50. **Not run**: the user stopped it before scoring; work moved to the one-copy simulation. (full: archive/2026-11/LEDGER_O3.md:10144)

## Copy phasing in Rust (`copy_phasing`, the port of `o3port.py`), 2026-10-24
Rust port `copy_phasing`: 81/81 golden cases byte-identical at 1 and 4 threads (after the final-review copy-order fix); P10 control arm Python 166.4 s / 1.9 GB vs Rust 19.1 s (1 thread), 5.4 s (4 threads). `candidate_copies --net-reads` added. Register row 1402. (full: archive/2026-11/LEDGER_O3.md:10210)

## 2026-10-24 REPORT (register 1408): how much does Jim differ from his own reference, family by family?
1,415 copies of 494 families mapped to reference, mat, pat: every reference copy is present at 100 % in mat or pat (1,414/1,415); copy-count differences mat vs pat in 27-35/494 families (5-7 %; 23 under all settings), mostly by 1 copy; most within-Jim difference is allelic. Revises July's "zero autosomal inter-haplotype CNVs". (full: archive/2026-11/LEDGER_O3.md:10241)

## PREREG SDCN (2026-10-25, REPORT, no pass/fail): copy number of exons and of SD segments, Jim's mat vs pat
REPORT, register 1416. Exon and SD-window copy counts in Jim's mat vs pat (exons remapped with `-x sr` after `asm20` failed its own QC; disclosed). mat != pat: exons 7.4 % (noise 0.2 %), 111/463 families (24.0 %) with >= 1 differing exon vs 5-7 % at whole-copy level; SD 5-kb windows 44.8 % (noise 6.6 %). Advisor's claim holds at SD level. Later qualified by WGSCHK (r1419). (full: archive/2026-11/LEDGER_O3.md:10261)

## PREREG BPCLUST (2026-10-25): do breakpoint groups in real RNA mark families whose exon copy number differs from the reference?
Ref-absent SPLIT/CLIP breakpoint groups (>= 3 reads, not chrM) flag a family; truth = 1416 exon DIFF. Bars H1 precision >= 0.40, H2 Fisher p < 0.05. Stage 2 held-out: **FAIL H1 and H2** (register 1417): precision 0.389, sensitivity 0.093, p = 0.10. Most exon copy differences leave no new junction in transcripts. (full: archive/2026-11/LEDGER_O3.md:10306)

## PREREG EXHAP (2026-10-25): does the certified tier, run per exon, find exon copies Jim has and his reference lacks?
`copy_phasing` per exon unit; family flagged if >= 1 certified copy. Held-out: **FAIL H1 and H2** (register 1418): precision 0.286, sensitivity 0.039, p = 0.31. Putative copies enriched (3/9 vs 0.03 base, p 0.002) but catch 3 of 150. Its reading ("RNA cannot see most gains") is **corrected by EXCEIL (r1420)**: the detector is the limit. (full: archive/2026-11/LEDGER_O3.md:10344)

## PREREG WGSCHK + EXCEIL (2026-10-25, reports with one reading rule each), committed before either is computed
- WGSCHK (register 1419): Jim's 12x WGS sides with the assemblies in 64.8 % of decided exons and 72.5 % of SD windows, below the 80 % rule: differences only **partly confirmed**; large (>= 25 %) differences are real (77 % / 94 %); WGS-supported core 29/463 families (6.3 %).
- EXCEIL (register 1420): level B (>= 2 reads closer to a haplotype-only exon sequence) reaches 41.9 % of gaining units (control 14.0 %) > 25 %: RNA's reach is NOT limited by the copies; the per-exon detector (2-4 %) is the gap (likely Theorem A free slots). (full: archive/2026-11/LEDGER_O3.md:10380)

## MMDEMO (2026-10-25, register 1421, descriptive demonstration on one gene, not a test)
PGM5: its 139 reads' 540 secondaries hit exactly the 6 loci a direct genome search finds (silent copies), so multimappers enumerate the reference's copies; but they do not show the copy is in the individual (the chr12 copy is absent from Jim's maternal chr12 yet still collects 92 secondaries). (full: archive/2026-11/LEDGER_O3.md:10431)

## PREREG TRIOCHK (2026-10-26): are Jim's haplotype-only copies real, and do his parents differ?
Trio WGS 21-mer dosage at NPIP maternal-only locus (N), STON1 121.6-kb insertion (S), father's deletion junctions (P). Result (register 1422): **all three CONFIRMED**, not assembly errors; parents differ (Dolly 2 copies, Trib 1, Jim 1 at both loci): a segregating insertion/deletion. (full: archive/2026-11/LEDGER_O3.md:10445)

## Rows 1504, 1506, 1513 (2026-11-02..03)
`archive/2026-11/PREREG_missing_copy_nont2t_2026-11-02.md` (r1504 FAIL: union recall 0.031 human, 0 gorilla; calls are fragments / mixtures of present loci), `PREREG_candidates_parallel_families_2026-11-02.md` (r1506 KEPT: clusters 385 -> 189 s, total 846 -> 646 s, 9 products identical), `PREREG_secondary_pool_2026-11-03.md` (r1513, numbers in the current state; scorer `bench/o2_contested/pool_score.py`, outputs workspace `o2pool/score{A,B}.contigs.tsv`).

## Do-not-quote list and what would settle O3 (folded from O3_STATUS.md, 2026-11-03)
Full text: `archive/2026-11/O3_STATUS.md`. Do not say: "289 candidates, 0 admitted" (12-family MEC panel only; 8/328 on 73 families); "TPR 0.7965 / FPR 0 / AUC 0.8924" (retracted oracle statistic); "8/101 = 7.92% haplotype CNV rate" (defensible <= 4/101, tandem-confirmed 1/101); "8/9/9 vs 5/6/8 haplotype deficit" (direction survives, magnitude not); "2.78% mat-vs-pat difference rate" (null 0.1512); "the DNA-side 0/817 is an independent second zero"; "the NPIP chr16:80,195,301 novel gene" (it is NPIPB15); "multi-copy families do not differ between individuals" (the zero is on 99 housekeeping loci, §6dq); "GOLGA6L7 is a reference-absent finding" (positive control). Three-way decomposition: specific RNA statistics REFUTED (k+1 MEC, between-copy depth, near-tie proxy, WGS k-mer multiplicity, `collapse_gate`); headline zeros VACUOUS / BLOCKED (FARDIV ∧ FARCLIP 0/915 cannot fire; unmapped route power 1/35; gate 5 rejects 172/320); the phenomenon DEMONSTRATED at DNA level. What would settle it: matched-tissue Iso-Seq from a second gorilla; put `LOC134757045` in front of the RNA detector; fix gate 5 before quoting any admission rate.

## Rows 1578-1579 (2026-11-07/08): a second gorilla's Iso-Seq on Jim's T2T
Advisor: a non-Jim sample misses many copies. Divergence piles (r1578) OR6737 vs depth-matched Jim: multi-copy 7.5% vs 1.8% (+5.7, NOT refuted as registered, borderline); shared loci 2.5% vs 0.3%. Three-haplotype test (r1579): 7.1% vs 4.6% (Jim = noise floor) on 7,056 shared loci, claim FAILS. Post hoc 12 strong loci (3 haplotypes AND a pile), 7 in de novo families: NBPF (MCL5, 60 independent pairs), IGHV-like (MCL67/MCL3), NACA-like, INTS4-like, SMPD4-like, ZNF717. Testis vs fibroblast: not the matched-tissue design asked for in the do-not-quote list; no third public HiFi Iso-Seq gorilla. Scripts `bench/crossind/{piles,threehap}.py`; outputs `/anvme/workspace/v138da13-rustle/crossind/`. `docs/PREREG_crossind_missing_copies_2026-11-07.md`.
