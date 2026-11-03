# Ledger: O2 (copy assignment of tied reads)

Compact version (2026-11-01). The verbatim ledger is `docs/archive/2026-11/LEDGER_O2.md`; each entry ends with a pointer to its line there. Original 53 dated docs: `docs/archive/2026-10/merged/`.

## Current state (2026-11-03)

- **Scope (user, 2026-11-03):** O2 is only for reads minimap2 could not place. The contested set = molecules whose runner-up is within 10 on `ms` (the score minimap2 ranks by), any AS, minus the K = 0 wall (same scores, same shape, same chromosome: nothing to decide). The AS-gap class (spliced primary vs unspliced higher-AS rival, ms margin > 10) is OUT (r1403/1404: primary right 99.6%). Sizes: A119b 283,924 molecules (wall 27%), HG002 197,872 (22%), OR6737 15,964 (25%), KB3781 52,061 (17%) (r1508 census, post hoc amendment).
- **Default `copy_assign`:** exact AS ties with the PSV certificate (`as_margin` 0), unchanged. **Opt-in:** `--tie-ms W` (ms gate; ms columns only under the flag; never fires on ms 0; needs an indexed BAM), `--primary-on-family` (r1514), `--tie-structure`, `--certify-all-records`, `--union-certificate` (r1481 FAIL), `RUSTLE_PSV_BLOCKS=1`, candidates stage.
- **Truth:** KB3781 fibroblast = Jim = the mGorGor1 individual; `bench/o2_truth/haplotype_truth.py` aligns reads and copies to his mat/pat haplotypes; `bench/o2_contested/strict_truth.py` reads it at slack 0 (unique best hit; copy loci need a near-perfect best per haplotype; `STRICT_LOCI=best2` = the old rule). Rows whose best hits sit on >= 2 loci are UNVERIFIABLE (the other-allele trap of r1413), never successes.
- **Proof status (r1509, r1510, r1514):** on the 40-family panel `--tie-ms 10` assigns 16x more than the default gate and at slack 0 is right on 284 of 286 vs the primary's 280 (r1509, post hoc truth; FAIL as registered). Held-out 454 families FAIL (r1510): O2 = primary (143/161 vs 147/157) because 77% of held-out assignments were molecules whose PRIMARY lies on no catalog copy (uncatalogued parent genes; MCL148 3,283, MCL13 148), pulled in by a secondary within 10 ms. `--primary-on-family` removes all of them (r1514: held-out assigned 4,529 -> 1,109; decidable 129/22 vs primary 125/26; panel keeps all 1,275 and gains 159; FAIL as registered on two side bars). On 151 decidable held-out molecules O2 beats the primary on 5 and loses on 1: the near-tie set's errors are shared with the aligner; this round's gain is correctness (no assignments to reads from outside the catalog), not precision.
- **Secondaries (r1511/1512):** useful fraction by ms ratio: <0.5 4.8%, 0.5-0.8 6.9%, 0.8-0.9 16%, 0.9-0.95 15%, 0.95-0.99 33%, >= 0.99 71%; the informative set is ratio >= 0.5 (minimap2's `-p 0.8` would drop 23% of useful secondaries; `-p 0.1 -N 50` was right). A gene-level graph from multimappers FAILS as family discovery (pair precision <= 0.032 vs O1 0.341; silent members not found).
- **Standing facts:** minimap2 ranks by `ms` (r1409); ms and AS are one unit per record, a molecule's margins differ by splice shape, so a near tie is read on both (r1508); assemblers (StringTie2/FLAIR) take the primary only, so on the near-tie set they commit to a coin toss; PSVs fitting no catalog copy = O3 input (r1484: O2 is a pileup + caller).
- **Next:** re-register `--primary-on-family` with bars on the decidable class (then default under `--tie-ms`); a shape guard (spliced read vs single-exon copy, the MCL148 class); the MCL183 allele-vs-paralog case (mat allele of copy 0 fits copy 1 at the PSVs) is the r1413 trap and needs haplotype-aware candidates.

## Contents

Merged 2026-10-24 from 53 dated docs (originals in `docs/archive/2026-10/merged/`), later sections appended; compacted 2026-11-01. Each entry below keeps its original `## ` headings; the dated-file heading and its title heading share one entry. (full: archive/2026-11/LEDGER_O2.md:6)

## O2_DEFAULT_ROSTER_2026-10-06.md
## O2 with read-level truth on the default families roster: NPIP, TBC1D3, Y (2026-10-06)
Simulation, DEV, A119b/CHM13. NPIP: B1 UNDERPOWERED (1 assigned); B2 PRIMARY UNDERPOWERED (12), B2 ANY UNDERPOWERED (29, 0 correct); B3 nothing assigned. TBC1D3: no MAPQ-0 read, all bars UNDERPOWERED. Contested simulated reads are rare (NPIP 33 of 1,099); real NPIP far more contested (24% MAPQ 0). Full chrY run (173-copy family) NOT done. (full: archive/2026-11/LEDGER_O2.md:65)

## PREREG_o2_default_roster_2026-10-06.md
## PREREG: O2 with read-level truth on the default families roster, at NPIP, TBC1D3 and the Y-ampliconic genes (2026-10-06)
Rules: B1 OWN 0 wrong/0 conflict with >= 30 assigned; B2 PRIMARY and ANY accuracy >= 0.95 (>= 0.90 per band); < 30 assigned = UNDERPOWERED; no rescue. Amendments 1-2 (skip-poa equivalence; chrY without the 173-copy family). Results above. (full: archive/2026-11/LEDGER_O2.md:136)

## PREREG_asj_o2_2026-10-10.md
## Allele-specific junctions as an O2 extension: haplotype assignment by the O2 certificate, allele vs copy by the haplotype assemblies, 2026-10-10
KB3781. T1 discordance dev 0.00099, held-out 872/883,715 -> 1,850/1,285,832 = 0.00144 (bar 0.005) PASS; T2 median concordance 1.000 ALLELIC (2,455 anchors) / 0.577 NON-ALLELIC (bars >= 0.95 / <= 0.80) PASS. Verdict: the claim holds. (full: archive/2026-11/LEDGER_O2.md:229)

## PREREG_o2_population_2026-10-12.md
## How many multi-placement reads can O2 assign, and how many does it move off the aligner's primary? 2026-10-12
Descriptive. Full population (Amendment 1): gorilla 49,444 multi-placement: moved 519 (1.0%), confirmed 76.4%, identical 7.9%, undecided 14.7%; human 532,876: moved 14,945 (2.8%), 62.5% / 10.6% / 24.2%. Many moved reads are the r1270 retrocopy pattern. (full: archive/2026-11/LEDGER_O2.md:339)

## PREREG_o2_pruning_and_eichler_2026-10-12.md
## O2 on the coin-toss genes, O2-pruned good pool, and Eichler's AS-margin rule, 2026-10-12
E1 holds (tied S2 secondaries certified to primary copy 41,569 vs 75 gorilla; 82,882 vs 386 human). E3 "agree mostly" FAILS: Eichler-O2 agreement S1 0.356/0.540, S2 0.141/0.300; disagreements are spliced parent reads aligning unspliced to retrocopies (no intron in 95%/92%). (full: archive/2026-11/LEDGER_O2.md:422)

## PREREG_popfix_2026-10-12.md
## PREREG — population-fixed PSVs for O2 on the human Y (2026-10-12 plan; written 2026-10-13, before any run below)
Follows YAG split-half 0.544 FAIL. Verdicts: P0-fixed PASS (9,670/9,670); (b) B1 UNDERPOWERED, B2 PASS (0 of 1,458 assigned); (c) C1 UNDERPOWERED; (c') C'1 UNDERPOWERED, so no registered verdict on donor-vs-artifact; (d) PASS; (e) VACUOUS. (full: archive/2026-11/LEDGER_O2.md:518)

## PREREG_psv_semiglobal_engine_2026-10-12.md
## PREREG — a log-odds, semi-global, anchored-band engine for O2's copy-vs-copy PSV alignment (2026-10-12)
Held-out 1 gorilla P1-P4 PASS; held-out 2 human P4 FAIL (10 families > 2 s), P1 0.997 PASS; held-out 3 chimp P1 0.993 vs 0.422, P2, P3 PASS, P4 per-family 2 s FAIL (MCL9 3.64 s). P5 PASS (assignments byte-identical). Verdict as registered: not all bars pass; adoption by user decision 2026-10-12 (see section). (full: archive/2026-11/LEDGER_O2.md:731)

## PREREG_psv_significant_blocks_2026-10-12.md
## PREREG — PSV columns only inside significant local blocks of the copy-vs-copy alignment (2026-10-12)
B1 PASS, B2 FAIL (aligned pairs retained 0.257, bar 0.99), B3 PASS (0.9987 non-homologous removed). Not default; OPT-IN `RUSTLE_PSV_BLOCKS=1`. (full: archive/2026-11/LEDGER_O2.md:867)

## PREREG_psv_sun_pileup_2026-10-12.md
## PREREG — SUNs vs our PSVs (M1) and a per-copy pileup check of PSVs (M2), gorilla KB3781 (2026-10-12)
Amendment 1: S1 PASS (99.6%), gate V FAIL (118,169/216,593 = 54.6% verify). Results in PSV_SUN_PILEUP below. (full: archive/2026-11/LEDGER_O2.md:932)

## PREREG_s1_o2_certificate_2026-10-12.md
## Are the secondary-only copies expressed? O2's certificate on their reads, 2026-10-12
Result: the secondary-only copies are NOT supported by O2; their reads belong to the other copy (qualifies register 1268 C1). (full: archive/2026-11/LEDGER_O2.md:1022)

## PREREG_yag_o2_2026-10-12.md
## PREREG — does O2 work on the Y ampliconic genes (YAGs)? (2026-10-12)
Arm A sim: 1,416 of 1,460 contested reads have an identical rival; A-P1 UNDERPOWERED, A-P2 FAIL (0/2), A-P3 PASS (0/1,458). Arm B A119b testis: B-P0 PASS, B-P1 split-half 31/57 = 0.544 FAIL, B-P2 UNDERPOWERED. (full: archive/2026-11/LEDGER_O2.md:1099)

## PSV_SUN_PILEUP_2026-10-12.md
## SUNs vs our PSVs (M1) + per-copy pileup check (M2) — results, gorilla KB3781 (2026-10-12)
Gate V FAIL: 54.6% of PSV columns verify (non-homologous artifacts); after fix 56.2%. M1: SUNs = 1,594 (30%) of homologous PSV columns. M2 per-copy pileup flags real deviations 61/62 vs 0.2% control, bar PASS. (full: archive/2026-11/LEDGER_O2.md:1182)

## PREREG_o2_anchor_2026-10-13.md
## PREREG — does O2 recover the right copy? Truth from the read's own unique anchor (2026-10-13)
Arm 1 gorilla dev: T0 truth gate 3,105/3,152 = 0.985 FAIL (no verdicts). Human held-out: A2 calibration PASS (0/42 undecidable assigned), A1 UNDERPOWERED (4 assigned). Arm 2 junction truth: J1 FAIL as registered (80/51) but not evidence against O2. Arm 3 split-half on reference animal: P0 PASS, S1 FAIL (0.224), S2 FAIL (0.097). Arm 4 haplotype truth: H1 28/29 right, UNDERPOWERED; wrong call = cross-family blind spot. (full: archive/2026-11/LEDGER_O2.md:1286)

## PREREG_o2_exact_tie_2026-10-13.md
## PREREG — O2 exact ties decided by rounding noise: the fix (2026-10-13)
Fix 7ac557e7. F1 PASS, F2 FAIL as registered, explained (60,189/60,214 rows identical), F3 PASS (3 MCL1623 wrong calls now ambiguous), F4 PASS. (full: archive/2026-11/LEDGER_O2.md:1646)

## PREREG_o2_tie_outside_gate_2026-10-13.md
## PREREG — the tie-outside gate depends on region sweep order: the fix (2026-10-13)
Fix 014126be + d5d982a3. G1-G4 PASS. (full: archive/2026-11/LEDGER_O2.md:1745)

## PREREG_o2_proof_2026-10-17.md
## Pre-registration: does O2 work? (held-out near-identical families, 2026-10-17)
Simulated, 8 held-out families, 2 profiles x 3 seeds x 6 worlds. Verdict: PASS (B1-B3 in all six cells, powered; mm2 errs on 53% of undecidable reads). (full: archive/2026-11/LEDGER_O2.md:1844)

## PREREG_o2_hg002_tied_2026-10-18.md
## Pre-registration: does O2 put AS-tied reads on the right copy? (HG002, the donor's own genome as truth, 2026-10-18)
Non-PAR world: UNDERPOWERED for all four systems (no system placed 100 of 584 truth reads; mm2 wrong on 58.6%); T3 holds; T5 holds for integrated, FAILS for current. Premise P1 FAILS: minimap2 right on 0.971 of untied reads (bar 0.99). PAR stratum reported, no bar. Corrections after review: no verdict change. (full: archive/2026-11/LEDGER_O2.md:2031)

## PREREG_o2_fair_heldout2_2026-10-19.md
## Pre-registration: O2 against minimap2's trusted calls on a deep new donor (2026-10-19)
Repeat of register 1341 (UNDERPOWERED) on a deep donor; bars F1-F4 unchanged. No result recorded in this section. (full: archive/2026-11/LEDGER_O2.md:2378)

## PREREG_o2_fair_heldout_2026-10-19.md
## Pre-registration: O2 against minimap2's own trusted calls, held out on a new donor (2026-10-19)
UNDERPOWERED: F1 fails (30 source reads placed, 100 needed); F2 0 wrong; F3 precision 1.000 vs MAPQ >= 1 baseline 0.894; F4 fails (ROPN1). (full: archive/2026-11/LEDGER_O2.md:2433)

## PREREG_o2_per_family_2026-10-19.md
## Pre-registration: O2 works family by family in the near-identical regime (2026-10-19)
Claimed 20 families work in sim, 26 on real reads; audit: the verdict does not hold (not held out; MAPQ >= 1 baseline passes in more families; in 46/48 cells O2 = minimap2 primary; floor of 5 too weak). (full: archive/2026-11/LEDGER_O2.md:2514)

## PREREG_o2_regimen_2026-10-19.md
## Pre-registration: O2's regimen on almost-tied reads, in simulation and on real reads (round 24 + A119b fragments, 2026-10-19)
Verdict FAIL: test A sim PASS (0 wrong of 349); test B A119b fragments B1 fails 7/312 = 2.24% (bar 2%), B2 precision 0.978 vs 0.656, B3 97.2%. (full: archive/2026-11/LEDGER_O2.md:2619)

## PREREG_o2_single_genome_2026-10-19.md
## Pre-registration: O2 in one genome -- the fewest-edit rule on CHM13-only worlds (round 23, 2026-10-19)
PASS under Amendment 1: F1 0 wrong of 2,045; F2 7,450/7,450 sets hold truth; F3 2,045 vs v17's 430; F4 PASS all cells. (full: archive/2026-11/LEDGER_O2.md:2767)

## PREREG_samplecopy_haploid_2026-10-19.md
## Pre-registration: the haploid rule of the O3 -> O2 sample-copy step, held out on HG002's chrX families (2026-10-19)
UNDERPOWERED: HS holds (0 placed), HB not reached (step without rule places none). (full: archive/2026-11/LEDGER_O2.md:2958)

## PREREG_samplecopy_hg002_2026-10-19.md
## Pre-registration: the sample's own copies place O2's ambiguous reads (HG002, held-out families, 2026-10-19)
Verdict FAIL: S1 4/86 wrong = 0.047 (bar 0.02), S3 fails (126 right vs 48 but 11 wrong vs 7); S2 holds (7 of 7,228 own_tie). (full: archive/2026-11/LEDGER_O2.md:3072)

## PREREG_samplecopy_par_2026-10-19.md
## Pre-registration: the O3 -> O2 sample-copy step, a second held-out test (HG002 pseudoautosomal families, 2026-10-19)
Verdict PASS: 508 right / 0 wrong / 114 abstain of 622, own_tie 0/6,681. (full: archive/2026-11/LEDGER_O2.md:3226)

## PREREG_samplecopy_yag_testis_2026-10-19.md
## Pre-registration: the O3 -> O2 sample-copy step on the Y-ampliconic families in testis, by split-half agreement (2026-10-19)
Y1 holds: 187/187 tied reads agree across half-models; reproducibility only, correctness not claimed (post-hoc doubt they are provenance). (full: archive/2026-11/LEDGER_O2.md:3330)

## O2_VAULT_METHODS_2026-10-20.md
## The ThesisVault's methods, assessed for O2 (2026-10-20)
Survey table: SDA/PSV lineage pertinent (reference PSVs fail fresh data, r1348/r1350; corrector the only held-out gain, r1343); set reporting no better than AS null; EM/MMR (row 1351) and facility location (row 1352) not applicable (exploratory); phasing outside O2. (full: archive/2026-11/LEDGER_O2.md:3467)

## O2_WINDOW_SWEEP_2026-10-20.md
## O2's window past Eichler's 10 AS: a development sweep on gorilla KB3781 (2026-10-20)
Dev only. Widening adds reads minimap2 already gets right (18 of 80,151 wrong in 30-49 band); nearly all minimap2 errors are within 30; the set holds truth on 99.6-100%. (full: archive/2026-11/LEDGER_O2.md:3537)

## PREREG_o2_corrector_2026-10-20.md
## Pre-registration: O2 as a corrector of minimap2 on held-out reads (2026-10-20)
Overall FAIL (S2 fails); H1 A119b passes (0.8171 vs 0.7147, 130 fixes / 3 breaks), S1 passes, H2 public testis UNDERPOWERED. Audit: holds with caveats. (full: archive/2026-11/LEDGER_O2.md:3616)

## PREREG_o2_corrector_families_2026-10-20.md
## Pre-registration: does the sample-aware corrector recover a share of minimap2's errors across many families? (2026-10-20)
Amendments 1-2: gorilla KB3781 mainstay families. UNDERPOWERED (minimap2 wrong on 2, bar 10); Amendment 3 rescore UNDERPOWERED (3); Amendment 4 positional key UNDERPOWERED, final for this experiment; audits: holds with caveats, direction not established. (full: archive/2026-11/LEDGER_O2.md:3811)

## PREREG_o2_fair_a119b_fresh_2026-10-20.md
## Pre-registration: O2 against minimap2's trusted calls on fresh A119b reads (2026-10-20)
PASS: F1 197 placed, F2 0 wrong (95% upper 1.9%), F3 precision 1.000 vs MAPQ >= 1 baseline 0.845, F4 18 families. Audit: holds with caveats; the public-donor test (register 1341) is still UNDERPOWERED. (full: archive/2026-11/LEDGER_O2.md:4121)

## PREREG_o2_known_family_smn_2026-10-20.md
## Pre-registration: a known family, SMN1/SMN2: does O2's pick of the real primary give better loci than minimap2's, StringTie's and FLAIR's? (2026-10-20)
Overall FAIL: L1 PASSES (each copy its own isoforms), L2 FAILS (StringTie/FLAIR do not keep the assignment; after Amendment 2 on FLAIR alone), L3 UNDERPOWERED (A119b barely expresses SMN2). (full: archive/2026-11/LEDGER_O2.md:4237)

## PREREG_o2_personal_copies_2026-10-20.md
## Pre-registration: personal copies against the corrector, on gorilla KB3781's own reads with exact truth (2026-10-20)
Verdict FAIL by the pre-registered bars: O2 changed 3 of minimap2's 887 calls (2 fixed, 1 broken); P holds only through SM392. Audit: holds with caveats. (full: archive/2026-11/LEDGER_O2.md:4441)

## PREREG_o2_second_round_2026-10-20.md
## Pre-registration: the second round -- loci rebuilt where O2's corrected multimappers have a say (2026-10-20)
HG002 Kinnex. FAIL: S1 and S2 fail, S3 holds (15 families, 29,965 truth reads); LOOP-O2 pooled TVD 0.112. Audit: FAIL stands. (full: archive/2026-11/LEDGER_O2.md:4878)

## PREREG_o2_second_round_donor_2026-10-20.md
## Pre-registration: the second round on a held-out donor (HG00733, simulated from its own assembly) (2026-10-20)
FAIL: gate holds (ORACLE 0.970); B1, B2 fail, B3, B4 hold; LOOP-O2 right 0.477 / wrong 0.079 / echo 0.428 vs primaries right 0.619. Audit: FAIL stands. (full: archive/2026-11/LEDGER_O2.md:5069)

## PREREG_o2_set_null_band_2026-10-20.md
## Pre-registration: the Pareto set against a null arm, on the near-tie band alone, fresh one-genome simulations (2026-10-20)
Round 26 band: PASS (8,955 reads, set holds truth 99.98%, one copy 56.3%, singletons right 99.96%). Substrate G gorilla: UNDERPOWERED (175 reads), would have been FAIL (H 98.86%, Q 97.6%). (full: archive/2026-11/LEDGER_O2.md:5274)

## PREREG_o2_smn_showcase_2026-10-20.md
## Pre-registration: the SMN1 / SMN2 showcase — PSV-decided reads, assign-or-abstain, and an assembler that honours them (2026-10-20)
FAIL (R2 on SMN2, R4 on SMN2): 3,623 tied reads, set decides 2,861; O2 moves 8 reads vs truth 616. Audit: FAIL stands. (full: archive/2026-11/LEDGER_O2.md:5432)

## PREREG_o2_subs_2026-10-20.md
## Pre-registration: O2 by substitutions only, on fresh substrates (2026-10-20)
S1 one genome FAIL (G, F; 0 fixes / 0 breaks), S2 A119b FAIL (G, S, F; 4 / 12), S3 gorilla no verdict (1 / 26). Audit: FAILs stand; S1 FAIL built in (no donor alleles). (full: archive/2026-11/LEDGER_O2.md:5572)

## PREREG_o2_window_set_2026-10-20.md
## Pre-registration: O2 as a pool reducer in a 30-point intron-free window, on fresh substrates (2026-10-20)
S1 one genome PASS; S2 HG002 FAIL (H, Q). Audit: S1 PASS certifies the 30-point window, not the Pareto set (null arm missing). Amendment 2 FCGR2A erratum, no verdict change. (full: archive/2026-11/LEDGER_O2.md:5697)

## PREREG_o2_corrector_idealsim_2026-10-21.md
## Pre-registration: the frozen corrector on ideal-expression reads simulated from both of KB3781's haplotypes (2026-10-21)
After register 1361 (real reads UNDERPOWERED). FAIL in both replicates (3781: 8 fixes / 7 breaks; 3782: 10 / 4; R and S fail). Audit: FAIL holds with caveats; TBC1D3 diagnosis withdrawn. (full: archive/2026-11/LEDGER_O2.md:5892)

## O2O3_EXPLORE_2026-10-22.md
## EXPLORATORY (post hoc, not pre-registered)
Can a truth-free O3 gate make O2-VC abstain on bad copies? No: on HG00733 no gate brings E1 under 2% (G3 best, 4.14%); PRAMEF7 (69 of 176 wrong) leaves no VCF trace; HG002 gates remove 93-98% of wrong but cost SMN2/PDPR. Recommendation: no VCF-level copy gate; RNA-only candidate G3+H5 (1.37% on HG00733 is a dev number). Registered verdicts unchanged. (full: archive/2026-11/LEDGER_O2.md:6083)

## O2VC_DEV_2026-10-22.md
## O2 variant calls: the development run on HG002 (S19SMN2, S19MSTRG3150), 2026-10-22
## P1 (keep): the one longcallR run
## P2 (unknown) and P3 (ref) reuse P1's VCFs (copied into dev_P2/vcf, dev_P3/vcf; byte-identical, 8 files)
## each policy scored twice: no panel, then the panel A view of SMN2
(The three lines above are command-block headings of the original.) Dev view on SMN2: P1 keep PASS (0/1221 wrong, 612 fixes); P2 FAIL (A3 612 > 612); P3 ref FAIL (K1 98/128). Frozen rule B8 chose P1 (P3 excluded; P1 = P2 tie). (full: archive/2026-11/LEDGER_O2.md:6427)

## O2_NAIVE_TIED_2026-10-22.md
## Naive tied-read policies (NAIVE-BEST, NAIVE-ALL) against minimap2's primary: read level and copy level, 2026-10-22
NAIVE-BEST no gain over primary on HG002 (wrong/placed 0.388 vs 0.388), modest on HG00733 sim (0.363 vs 0.382), far from VC (0.059/0.063). NAIVE-ALL recovers expressed copies (7/11, 17/18) but invents expression (prec 0.615 vs 0.784; 0.734 vs 0.836): recall-only lever. (full: archive/2026-11/LEDGER_O2.md:6609)

## PREREG_o2_variant_calls_2026-10-22.md
## Pre-registration: O2 from the sample's own variant calls (2026-10-22)
## A + C-human: HG002 real Kinnex reads, every pool family (chrM families out); truth = HG002's own assembly (answer key only)
## B: DAZ, A119b testis, no truth (descriptive)
## E: HG00733 simulated from its own haplotypes, held out (never run or scored during development)
(Panel headings of the original command block.) Results (code 07441e7c, audited): A GTF2IRD2 PASS (44 placed / 0 wrong / 24 fixes; caveats); PDPR PASS; NBPF1 FAIL (9/14 wrong); K1/K2 PASS (49/9,642; 12/15,725); C-human PASS by scorer but audit DOES NOT HOLD (relabelled C1 FAILS); B DAZ descriptive (VC places 0); E HG00733 FAIL (E1 176/2,801 = 6.3%; PRAMEF7+10 = 133). (full: archive/2026-11/LEDGER_O2.md:6726)

## PREREG_o2_variant_calls_expected_failures_2026-10-22.md
## Pre-registered expected failures: O2 from the sample's own variant calls (2026-10-22)
Committed before the dev run and panels; blind literature survey. (full: archive/2026-11/LEDGER_O2.md:6947)

## O2 from the sample's own variant calls: vault survey, pre-registered expected failures, novelty delta
Expected failures F1-F16 (teacher contamination, absorption, PRAMEF CNV, ASE, gene conversion, SNV-only, low expression, DAZ near-identity, label artifacts...) and guesses: GTF2IRD2 FAIL, NBPF1 FAIL, K1/K2 PASS, C1 FAIL, E lean PASS; plus novelty delta vs vault methods. (full: archive/2026-11/LEDGER_O2.md:6981)

## PREREG_o2_forced_heldout_2026-10-23.md
## Pre-registration: O2 as forced haplotype-to-copy matches, simulated gorilla held-out (P7), 2026-10-23
`o3port.py` 791c0bcc. Pooled over 7 runs: 53 placed, 24 wrong (45%; all GWFAM328); O2-1 FAILS (bar <= 2%). Verdict: DOES NOT PASS. (full: archive/2026-11/LEDGER_O2.md:7416)

## PREREG_o2_p9_abstain_2026-10-23.md
## Pre-registration: O2 v2 on fresh gorilla fibroblast families with a synthetic near-identical duplicate (P9), 2026-10-23
PASS, but only by abstaining: Q1 0 wrong of 0 placed, Q2 0 reads on the lacked copy X'. (full: archive/2026-11/LEDGER_O2.md:7474)

## PREREG_o2_semi_p11_2026-10-23.md
## Pre-registration: O2 `--semi` (two-site lead) on fresh gorilla fibroblast families (P11), 2026-10-23
FAILS S2 (2 reads on X' in the missing case); S1 0.09%, S3, power (4,586 placed) PASS. DOES NOT PASS; fixes 0, breaks 0. (full: archive/2026-11/LEDGER_O2.md:7527)

## PREREG_o2_v2_heldout_2026-10-23.md
## Pre-registration: O2 v2 (forced matches, weighed on every tied copy) on fresh erased-copy simulations (P8), 2026-10-23
Cannot be run: no fresh gorilla family meets the pick rule; nothing simulated or scored. (full: archive/2026-11/LEDGER_O2.md:7582)

## PREREG_o2_guard_p13_2026-10-24.md
## Pre-registration: O3 -> O2 with `--semi --guard` on fresh gorilla fibroblast families (P13), 2026-10-24
FAILS S2 (2 X|B reads on X'; indel beside a distinguishing site); S1 0.09%, S3, power, E1 PASS. DOES NOT PASS. (full: archive/2026-11/LEDGER_O2.md:7629)

## PREREG_o2_guard_p14_2026-10-24.md
## Pre-registration: O3 -> O2 with `--semi --guard` and the site homology check (P14), 2026-10-24
SITE_FLANK = 5 (register 1396). PASSES: S1 4/3,431 = 0.12%, S2 0 on X', S3, power, E1 (0 placed). (full: archive/2026-11/LEDGER_O2.md:7679)

## PREREG_o2_semi_sites_p12_2026-10-24.md
## Pre-registration: O2 `--semi` with supported distinguishing sites on fresh gorilla fibroblast families (P12), 2026-10-24
Register 1391. PASSES: S1 0.02%, S2 0 on X', S3, power (4,562). At d 1 every read abstains by design. (full: archive/2026-11/LEDGER_O2.md:7715)

## REPORT_o2_refabsent_2026-10-24.md
## Report-only: O3 -> O2 on KB3781's real reference-absent copies (refabsent panel), 2026-10-24
No verdict. O3 flags vs expressed B-only loci per family; addendum with the IsoCon side (`candidate_copies`, 34 families). Cited later as "157 placed, 0 wrong". (full: archive/2026-11/LEDGER_O2.md:7752)

## The AS gap is a splicing artefact, not a copy signal (2026-10-24, register 1403)
KB3781 fibroblast BAM: a secondary with higher AS than the primary is common and almost always a worse alignment (unspliced to a processed copy). The Eichler AS >= 10 move is not a copy signal. (full: archive/2026-11/LEDGER_O2.md:7829)

## PREREG: splice-aware tie window on simulated gorilla reads with known origin (P15), 2026-10-24
PASS on all three bars (register 1404): H1 R fires on 812, primary right on 809 (99.6%); H2 99.65%; H3 window C holds 270/270 with 38% fewer reads. `--tie-structure` wired in as opt-in (register 1405); default unchanged. (full: archive/2026-11/LEDGER_O2.md:7862)

## PREREG: O2 placement with and without `--tie-structure` on fresh simulated gorilla families (P16), 2026-10-24
FAIL by rule (H1 0/0, O2 assigns nothing); uninformative (register 1406). (full: archive/2026-11/LEDGER_O2.md:7916)

## REPORT (no verdict): `--tie-structure` vs the default window on real KB3781 reads with haplotype truth (P5, P10 control arms), 2026-10-24
Register 1407: the window changes nothing in O2 on real reads. Also REPORT register 1409: minimap2 ranks by `ms` (dp_max), not AS; exact ms ties broken by read-name hash. (full: archive/2026-11/LEDGER_O2.md:7953)

## P18 pre-registration (2026-10-25): does an error-model score beat minimap2's `ms` on near-tie reads?
FAIL (H1, H3), register 1410: 2 switches; H2 barely (40 vs 42 wrong of 222); H3 220 assigned, 40 wrong. `ms` stays the ranking. Reports: register 1411 no tag/CIGAR feature separates the true copy; register 1412 abundance tie-breaker gives no information under even expression (dev only). (full: archive/2026-11/LEDGER_O2.md:7990)

## PREREG P21 (2026-10-25): does copy_phasing (O3 -> O2, read clusters on alignment columns) fix the other-allele errors of P18?
FAIL (H2), register 1413: 1,190 tied rows, 9 assigned (all right), fixed 0 / broken 0; H1, H3 PASS. Register 1414: the two false certified copies are not recurrent errors (one is the r1403 AS intron artefact). (full: archive/2026-11/LEDGER_O2.md:8059)

## PREREG P22 (2026-10-25): `--certify-all-records` held out: precision and recall of certified copies
FAIL H2 (register 1415): flag keeps 8/11 masked copies (73% < 90%); H1 PASS (1 false copy, same in every arm). `--certify-all-records` stays opt-in. (full: archive/2026-11/LEDGER_O2.md:8108)

## PREREG O2-TIE-MS (register 1472, 2026-10-30): O2's tie window defined on `ms`, the score minimap2 ranks by
o3port.py `--tie-ms` on P12 + P14 substrates. FAIL as registered: M0 places 0 reads on every control/missing run (B1, B2 FAIL, B3 vacuous); R fails on P12 (4,377/2 vs 4,562/1, pre-SITE_FLANK), P14 reproduces. Fixes = 0 in all 450 runs; M10 = A within 3 reads. `--tie-ms` stays opt-in. Human A119b: within-family ms ties only 71.4% inside the AS window (gorilla 99.95%), so an ms-anchored window is the better human anchor; candidate for an HG002 prereg. (Register number here is as written; see report.) (full: archive/2026-11/LEDGER_O2.md:8147)

## O2 AMBIGUOUS DECOMPOSITION (2026-10-31, descriptive; HG002 `hg.assign.*` and A119b `a1.assign.*` from the FAMASM runs; debug `~/o2diag/`)
Distinct reads with a row: A119b 212,176 / HG002 141,265; assigned 4,331 / 17,805; tied wall 47,428 / 62,158; weak 12,673 / 3,878; origin rejected 164,293 / 67,571; orphan 21,438 / 16,349. Origin rejections: >= 5% substitution class is junk (correct abstention); < 5% class ties outside the catalog (83-94%), origin 1-5% divergent, recoverable only by widening the catalog (candidates / semisec, never run on A119b or HG002). Follow-up: 75% of < 5% primaries sit in another family's copy (`--union-certificate` off); bins are mixtures of CNV paralogs. Lever order: union certificate, candidates/semisec, correlation clustering (O3). (full: archive/2026-11/LEDGER_O2.md:8242)

## PREREG O2 LEVERS (2026-11-01, register 1481; dev A119b; runs `~/o2levers_A119b/`, scorer `bench/o2levers/levers_score.py`)
Arms on the shipped A119b base: U = `--union-certificate`; C = `candidates` stage + `assign --candidates`. Distinct-read scoring: base assigned 4,331, double-claimed 1,613, unique-mapper agreement 35.9%, lt5 class 20,899 with a row (986 assigned).
Bars: U1 assigned >= 4,331; U2 double claims 0; U3 agreement >= 35.9%; U4 lt5 assigned >= 986; C1 >= 1 flagged candidate; C2 lt5 rescued >= 25%; C3 assigned >= 4,331; C4 agreement >= 35.9%. PASS proposes an HG002 registration; FAIL closes the lever as a default.
Full runs gave nothing: U (job 4403494) cancelled at 19 h 20 (chr13:13.9-15.6 Mb hotspot, 89 families hold 86% of non-assigned rows from 30,155 reads, median 62 families per read); C (job 4403495) vanished with no exit line (sacct down).
Amendment (2026-11-01 20:45, before any panel run): fast panel `~/o2levers_fast/`, 57 hash-ordered non-hotspot families (14,966 non-assigned rows), BASE and U side by side, 1 h 30 wall; bars relative to BASE plus U5 cost (U wall <= 2 x BASE; wall-limit hit = FAIL). PASS licenses the hotspot cost question, then HG002; FAIL closes the lever. Arm C waits for its end reason, then same panel.
**Panel result: U FAIL** (see the RESULT subsection below). (full: archive/2026-11/LEDGER_O2.md:8306; amendment at :8431)

## PREREG TIE-MS PORT (2026-11-01, register 1483; follows register 1477 "worth porting"; runs `~/o2ties_real/`)
Port `copy_phasing run --tie-ms W` (pile window ms >= best ms - W; `--semi` ranks by ms when all records have it; default untouched). Real OR6737 (446 families) + KB3781 (280), arms A / M10 / M0 under `--semi --guard`. Bars: R byte-identical A; B1 >= 5% of A's tied reads untied under M10 on both; B2 M10 assigned >= 0.95 x A; B3 < 1% conflicting.
RESULT (Slurm 4403643, port 0311c388): **FAIL as registered**. R PASS (9 product files byte-identical). B1 PASS (16.1% / 14.1%). B2 FAIL (1,289 vs 1,758 = 73.3%; 3,092 vs 6,544 = 47.2%); post hoc, most lost assignments are reads whose rival the ms window excluded, so the bar counted the intended effect; verdict stands. B3 FAIL OR6737 (24/935 = 2.6%), PASS KB3781 (0.4%). M0 assigns 314 / 1,655.
Decision: `--tie-ms` stays OPT-IN. Next if revisited: HG002 truth with B2 on final placement. (full: archive/2026-11/LEDGER_O2.md:8330)

## PREREG PSV VCF REPORT (2026-11-01, register 1484; descriptive report, advisor's question "why no pileup + variant calling for PSVs?")
Report only, no default change. `bench/o3_graph/psv_vcf.py` writes O2's variable columns as a VCF (ID `family:copy:POS`), per-read decision paths, and deciding PSV sites, on OR6737 arm A of 1483; bcftools 1.25 mpileup | call -mv on the same pile reads; expectation >= 90% agreement.
Result: 136 families, 1,758 assigned (1,334 sites / 131 haplotype / 293 certified copy); 5,796 variable columns on 449 copies; 1,898 deciding PSV sites, 451 with the rival's base in >= 6 pile reads.
As registered bcftools calls 54.6% (3,162/5,796): it caps base quality at MAPQ and 24% of pile records are MAPQ-0. With `--mapq-reset 60` (post hoc, disclosed): 93.9% (5,445); the remaining 351 are mostly minor alleles < 10% that a diploid prior calls hom-ref. bcftools' 18,138 extra SNVs: one_allele 14,759, bonferroni 2,037, p 1,342; 504 indels.
Reading: O2 is pileup + caller without MAPQ capping and without ploidy. Outputs `~/o2ties_real/psv/or.psv60.*`. (full: archive/2026-11/LEDGER_O2.md:8378)

### O2 LEVERS fast-panel RESULT (2026-11-01, job 4405574, `~/o2levers_fast/`; bars frozen at 11f04e56)
Answer in 10 min (BASE 620 s, U 582 s). By distinct reads, BASE -> U: assigned 63 -> 59, weak 595 -> 884, origin_rejected 7,716 -> 7,087,
orphan 2,412 -> 2,362, tied 2,413 -> 2,807; reads assigned in >= 2 families 3 -> 0; unique-mapper agreement 12/22 -> 13/23; lt5 class
72 reads with a row, 0 assigned in both. Bars: **U1 FAIL** (59 < 63), U2 PASS (0), U3 PASS (56.5% >= 54.5%), U4 PASS (0 >= 0, vacuous),
U5 PASS (582 s <= 2 x 620 s). **FAIL as registered: the union-certificate lever closes as a default** (flag stays opt-in).
Mechanism: U loses 11 BASE assignments (10 now tied, 12 rows ambiguous) and gains 7, changes no copy; the 3 BASE double claims are
dropped, not resolved. Of 5,563 molecules in scope the certificate scored 1,952: 311 went to an OUTSIDE locus (the catalog lacks
their origin, as the 2026-10-31 decomposition predicted), 11 to a family, 795 tied, 835 ambiguous; it turns origin rejections into
ties and weak rows, not assignments. Caveat (disclosed, not a reason to re-run): 3,611 of 5,563 (65%) were left as today because a
tie partner lies in a region not loaded; the panel's region list holds only panel copies, and the genome-wide run had the same
class (e.g. 1,641 of 3,623 on ~xchrom~MCL33). The hotspot cost question is moot after a FAIL.

### O2 LEVERS arm C on the fast panel: frozen before the run (2026-11-01 ~21:00)
Same 57 families and BASE as the arm U panel (`~/o2levers_fast/base.*`). Pipeline steps restricted to the panel: `candidate_copies
--families <panel>` (delta 0.00958, max-reads 1000, the stage defaults) -> `utilities candidate-augment` -> the candidate families'
nets realigned to `aug.fa` (`minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`) -> copy_assign on `aug.*` for the
candidate families (`--only-families aug.families.txt`) and on the original inputs for the other panel families, the stage_assign
flags as in BASE; tables concatenated as `assign_with_candidates` does. Wall limit 1 h 30 on spr1tb; a limit hit = FAIL.
Bars (relative to BASE on the panel): C1 >= 1 flagged candidate; C2 of the lt5-class reads with a BASE row in the panel (72), >= 25%
have their realigned primary on a `cand_` contig; C3 assigned reads >= BASE (63); C4 unique-mapper agreement >= BASE (12/22).
PASS -> the hotspot and HG002 registration; FAIL -> arm C closes as a lever. If no candidate is flagged (C1 FAIL), the assign
steps are skipped and the arm fails on C1.

### O2 LEVERS arm C: PARKED (2026-11-01 22:00, user)
Job 4405613 (panel, bars at bdad5c49) cancelled at 1 h 04 by the user's decision to park, still inside `candidate_copies` pass A
(building the family nets): 430 GB read from the 96 GB BAM after 43 min. Cause measured on MCL331 (4 copies, ~13k records): samtools
counts the four copy regions in 0.4 s; `candidate_copies` pass A takes 8 s and reads 3.4 GB for the same regions, so each copy
query reads far more of the BAM than its region (mechanism not yet confirmed in the code). The genome-wide run (7,987 copies, job
4403495) would need days in this step alone. **The lever is untested; no verdict is recorded** (the limit was not reached, and the
cost is the tool's, not the idea's). To resume: collect the nets during pass B's existing sequential sweep against the merged copy
intervals instead of one indexed query per copy, test that the nets are identical on a few families, then rerun the panel with the
same bars and limit as a disclosed amendment.

## O2 on the ms near-tie set, KB3781 haplotype truth (2026-11-03, register 1508-1514)
Pre-registrations in `archive/2026-11/`: `PREREG_o2_contested_census_2026-11-03.md` (r1508), `PREREG_o2_ms_tie_proof_2026-11-03.md` (r1509), `PREREG_o2_ms_tie_heldout_2026-11-03.md` (r1510), `PREREG_secondary_usefulness_2026-11-03.md` (r1511/1512), `PREREG_secondary_pool_2026-11-03.md` (r1513, O3 side), `PREREG_o2_primary_on_family_2026-11-03.md` (r1514). Tools: `as_table` columns 7-17 (ms, MAPQ, splice shape, best-AS rival, positions) and `--records` (every primary/secondary record); `copy_assign --tie-ms W`, `--primary-on-family`; `bench/o2_contested/{panel,score,strict_truth,multimap_graph,pool_score}.py`; `bench/o2_truth/haplotype_truth.py`. Data: workspace `o2contested/*.molecules.tsv`, `o2contested/haps/KB3781.{mat,pat}.mmi`, runs `o2mstie/`, `o2contested/r1511/`, `o2pool/`.
- r1508 CENSUS: contested 2.3 / 3.7 / 1.6% of molecules (A119b / HG002 / OR6737); K = 0 wall 16 / 3 / 6% of it; the "either score within 10" set was too wide (the AS-gap class is 33-85% of it) -> amended post hoc to ms margin <= 10.
- r1509 PROOF (panel 40 families): identity PASS; MS10 contested 40,524 vs BASE 24,095 molecules, assigned 1,265 vs 78; registered truth (NM + 2 slack) FAIL (35 decidable, 2 wrong = primary); slack 0 post hoc: 284 right / 2 wrong of 286 vs primary 280 / 6; 648 unverifiable.
- r1510 HELD-OUT: FAIL (see current state); 131 of 161 wrong in MCL13 where both O2 and the primary sit on copy 0 and the reads come from an uncatalogued parent (NC_073230.2:52.6 Mb, MAPQ 60).
- r1511 REPORT bound = ms ratio 0.5; r1512 FAIL (graph); r1514 guard (see current state). Dumps: `o2mstie/*.strict_dump.tsv`, `~/o2contested/r1510_strict.txt`, `r1514_strict.txt`.
