# Ledger: O2 (copy assignment of tied reads)

Compact version (2026-11-01). The verbatim ledger is `docs/archive/2026-11/LEDGER_O2.md`; each entry ends with a pointer to its line there. Original 53 dated docs: `docs/archive/2026-10/MERGED_REMOVED.md` (directory removed 2026-11-08; verbatim in `docs/archive/2026-11/LEDGER_*.md`).

## Current state (2026-11-03)

- **2026-11-22, O2 REPLACED (user; register 1903; `O2_FAMILY_BENCHMARK.md`).** This ledger is the record of the LEGACY O2 (copy assignment of tied reads). O2 is now the benchmark of known gene families across tools, on planted reads of known origin and on real samples. Why: on every truth tried no assignment rule beat the primary alignment (rows 1775-1785), nothing reads the assignments (row 1699), and the family tables per tool with planted truth answer the question the objective was for. Nothing below is retracted. In `all` the stage `o2` replaced `assign` in every plan (later the same day also under `--candidates` and `--legacy-catalog`, and without a switch); the `assign` stage itself is unchanged and runs when called by name.
- **2026-11-17, placement on the catalogue-free truth (rows 1775-1778, 1781, 1783): `docs/NEUTRAL_TRUTH_PLACEMENT_2026-11-17.md`.** No rule beats the primary alignment on HG002's 710 decidable near ties; the default assign stage abstains safely. The truth has a defect (span overlap without strand or exons; unlocatable placements dropped), found on chr1: do not score another arm or quote the table before the repair.
- **2026-11-17 / 18, Winnowmap and the size of the problem (rows 1782, 1784, 1787, 1790, 1791): `docs/WINNOWMAP_AND_FAMILY_TIES_2026-11-18.md`.** Winnowmap, stock or with its copy-aware method patched on for spliced reads, is worse than minimap2's primary on HG002's decidable near ties (605 and 609 right against 639 / 638) and has no effect at tandem copies; it fixes seeding, families need ranking. Near ties are 5.7% of HG002's family reads and fixable primary errors 0.05%, but 14 of 155 expressed families have >= 20% of their reads near-tied (81% of all near ties). A `tie` of the neutral truth is SEVERAL COPIES (unplaceable; 44% of those families' near ties) or ONE PLACE (the truth cannot judge; 15%): do not call a tie undecidable without that split.
- **2026-11-17 / 18, truth repair, the 10-point filter's reach, and its name (rows 1785, 1786):** the truth is repaired under `NT_REPAIR=1` (opt-in): HG002 near ties 600 decidable, primary 553 right / 47 wrong; every rule re-scored, none gains, verdicts unchanged. The adapted Eichler filter: "take the best AS above 10" is wrong on 16.7% (HG002) / 40.0% (KB3781) of reads with a lead of 11-50, "keep the primary" on 0.9% / 0.08%; no threshold at 10; primary right 98% at a lead of 6-10. The source (Guitart 2024) was re-read: its filter compares primary alignments of separate haplotype mappings per paralog cluster, so everything here is an adaptation to one reference. Definition and wording in `copy_assignment_definition.md`; all docs, comments, tool output and Fig. 5s reworded 2026-11-18 (identifiers keep the name, user's decision). Summary: `NEUTRAL_TRUTH_PLACEMENT_2026-11-17.md`.

- **Scope (user, 2026-11-03):** O2 is only for reads minimap2 could not place. The contested set = molecules whose runner-up is within 10 on `ms` (the score minimap2 ranks by), any AS, minus the K = 0 wall (same scores, same shape, same chromosome: nothing to decide). The AS-gap class (spliced primary vs unspliced higher-AS rival, ms margin > 10) is OUT (r1403/1404: primary right 99.6%). Sizes: A119b 283,924 molecules (wall 27%), HG002 197,872 (22%), OR6737 15,964 (25%), KB3781 52,061 (17%) (r1508 census, post hoc amendment).
- **Default `copy_assign`:** exact AS ties with the PSV certificate (`as_margin` 0), unchanged. **Opt-in:** `--tie-ms W` (ms gate; ms columns only under the flag; never fires on ms 0; needs an indexed BAM), `--primary-on-family` (r1514), `--tie-structure`, `--certify-all-records`, `--union-certificate` (r1481 FAIL), `RUSTLE_PSV_BLOCKS=1`, candidates stage.
- **rDNA exclusion on HG002 (row 1741, 2026-11-15; `PREREG_rdna_sequence_rule_2026-11-13.md` Amendment 2):** confirming run through the driver, INCOMPLETE as registered: the default arm was OOM-killed (239.8 GB), the arm with `RUSTLE_EXCLUDE_FASTA=tools/data/rdna_human.fa` finished in 1,631 s at 31.5 GB with 2,898 assigned rows. Not a held-out test; row 1691 stays FAIL. On the older completed default run the listed families hold 0 of 17,960 assigned rows (descriptive). Still opt-in; open = a completed default arm on today's binaries (1 TB node), then the user's decision.
- **Truth:** KB3781 fibroblast = Jim = the mGorGor1 individual; `bench/o2_truth/haplotype_truth.py` aligns reads and copies to his mat/pat haplotypes; `bench/o2_contested/strict_truth.py` reads it at slack 0 (unique best hit; copy loci need a near-perfect best per haplotype; `STRICT_LOCI=best2` = the old rule). Rows whose best hits sit on >= 2 loci are UNVERIFIABLE (the other-allele trap of r1413), never successes.
- **Proof status (r1509, r1510, r1514):** on the 40-family panel `--tie-ms 10` assigns 16x more than the default gate and at slack 0 is right on 284 of 286 vs the primary's 280 (r1509, post hoc truth; FAIL as registered). Held-out 454 families FAIL (r1510): O2 = primary (143/161 vs 147/157) because 77% of held-out assignments were molecules whose PRIMARY lies on no catalog copy (uncatalogued parent genes; MCL148 3,283, MCL13 148), pulled in by a secondary within 10 ms. `--primary-on-family` removes all of them (r1514: held-out assigned 4,529 -> 1,109; decidable 129/22 vs primary 125/26; panel keeps all 1,275 and gains 159; FAIL as registered on two side bars). On 151 decidable held-out molecules O2 beats the primary on 5 and loses on 1: the near-tie set's errors are shared with the aligner; this round's gain is correctness (no assignments to reads from outside the catalog), not precision.
- **Secondaries (r1511/1512):** useful fraction by ms ratio: <0.5 4.8%, 0.5-0.8 6.9%, 0.8-0.9 16%, 0.9-0.95 15%, 0.95-0.99 33%, >= 0.99 71%; the informative set is ratio >= 0.5 (minimap2's `-p 0.8` would drop 23% of useful secondaries; `-p 0.1 -N 50` was right). A gene-level graph from multimappers FAILS as family discovery (pair precision <= 0.032 vs O1 0.341; silent members not found).
- **Standing facts:** minimap2 ranks by `ms` (r1409); ms and AS are one unit per record, a molecule's margins differ by splice shape, so a near tie is read on both (r1508); assemblers (StringTie2/FLAIR) take the primary only, so on the near-tie set they commit to a coin toss; PSVs fitting no catalog copy = O3 input (r1484: O2 is a pileup + caller).
- **r1542 (2026-11-05, PASS, adopted):** the primary-on-family guard is ON under `--tie-ms` (`--no-primary-on-family` = r1510's gate); byte-identical to r1514 / r1510 on KB3781; human A119b panel removes 8.7% of unguarded assignments, all guarded assignments have a primary in the region windows. `--tie-ms` itself stays opt-in. Closed 2026-11-05: shape guard vacuous after the primary-on-family guard (r1540, 0 rows); within-sample allele-site guard FAIL (r1541: fixes 4, costs 12 held-out); the MCL183 haplotype-aware version needs assemblies and is unregistered.

- **2026-11-10:** the assign router (r1600: one BAM read, schedule from exact counts) FAILED as registered (583 s vs 553 s; 58.5 vs 36.4 GB; tables identical) and is dropped; branch `r1600-router` kept. Molecules O2 cannot decide (AS within 0.98 of best on several loci) are now reported as tie groups: one model per locus, one count per group (register 1603-1611, see LEDGER_O1 update 2026-11-10).

- **2026-11-13, rows not carried above (r1685, r1696):** an O2 guard that abstains when the assigned copy shares an SD pair with a sibling copy STOP (at fracMatch 0.99 it removes 25 of 26 wrong rows and 237 of 393 right ones; every assigned copy shares an SD). Per-bubble vote on the family graphs of r1649 STOP at step 0 (565 assigned, 553 right / 12 wrong, but 8 vote-only molecules vs bar 30): containment, k-mers and the vote agree on these families; population-haplotype paths are the only untried graph use.

## Contents

Merged 2026-10-24 from 53 dated docs (originals in `docs/archive/2026-10/MERGED_REMOVED.md` (directory removed 2026-11-08; verbatim in `docs/archive/2026-11/LEDGER_*.md`)), later sections appended; compacted 2026-11-01. Each entry below keeps its original `## ` headings; the dated-file heading and its title heading share one entry. (full: archive/2026-11/LEDGER_O2.md:6)

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
## O2 on the coin-toss genes, O2-pruned good pool, and the AS-margin rule adapted from the Eichler lab's filter, 2026-10-12
E1 holds (tied S2 secondaries certified to primary copy 41,569 vs 75 gorilla; 82,882 vs 386 human). E3 "agree mostly" FAILS: margin-rule (adapted Eichler filter) vs O2 agreement S1 0.356/0.540, S2 0.141/0.300; disagreements are spliced parent reads aligning unspliced to retrocopies (no intron in 95%/92%). (full: archive/2026-11/LEDGER_O2.md:422)

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
## O2's window past the adapted Eichler filter's 10 AS: a development sweep on gorilla KB3781 (2026-10-20)
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
KB3781 fibroblast BAM: a secondary with higher AS than the primary is common and almost always a worse alignment (unspliced to a processed copy). The adapted-Eichler AS >= 10 move is not a copy signal. (full: archive/2026-11/LEDGER_O2.md:7829)

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
- 2026-11-10, correlation clustering in the assignment layer (opt-in measurements, no default change): r1624 allele vs paralog STOP without a run (the distinguishing quantity is the r1541 allele-site guard; class with truth = 12 reads, MCL183); r1625 pooling a family's unassigned reads with `copy_phasing run --semi --guard` (r1483 arm A, KB3781, strict haplotype truth): 120 calls on the 1,113 abstained decidable rows of the r1514 arm, 115 right / 5 wrong against the primary's 118 / 2 on the same rows, 7 calls in held-out families; FAIL on the aligner bar, closed. `PREREG_cc_readthrough_allele_2026-11-10.md`, `PREREG_pooled_family_reads_2026-11-10.md`.

## Placement without the aligner's primary: exon space, k-mers, graph (2026-11-12, register 1642-1650)
User question: is there a method other than minimap2's spliced genome alignment that places multi-mapping long reads
better? Set = the ms near-tie molecules of r1514 (KB3781) and, new, the same on HG002; truth = the individual's own two
haplotypes at slack 0, per molecule. Pre-registrations in `docs/`: `PREREG_exon_space_placement_2026-11-12.md` (1642),
`PREREG_exon_space_multi_isoform_2026-11-12.md` (1643), `PREREG_exon_space_unspliced_hg002_2026-11-12.md` (1644),
`PREREG_kallisto_long_placement_2026-11-12.md` (1647), `PREREG_oarfish_placement_2026-11-12.md` (1648),
`STEP0_graph_placement_2026-11-12.md` (1649), `PREREG_kallisto_placement_or6737_2026-11-12.md` (1650).
Tools: `bench/o2_contested/{exon_place,multi_iso_ref,kallisto_place,oarfish_place,graph_place,place_checks}.py`;
`strict_truth.py` `STRICT_DUMP_OUTSIDE=1` writes the outside-catalog rows (negative set).
- r1642 FAIL: unspliced `map-hifi` of the reads to one spliced sequence per copy: tool right / primary wrong 15, reverse
  44; assigns 0.456 of outside-catalog molecules.
- r1643 FAIL: several sequences per copy (catalog model, assembled transcripts, unspliced span): equals the primary (31 vs
  36), negatives 0.587.
- r1644 FAIL: the post hoc unspliced-primary split of r1643 on HG002: 11 vs 43 (4 vs 214 before the scoring fix).
  Exon-space re-alignment WITH minimap2 is closed, both references, both classes.
- r1647 INVALID as registered (validity gate 0.80 on pseudoaligned reads; got 0.41 / 0.70), OPEN at the time, since CLOSED by r1651 (HG002 on CHM13: reach 0.043, 8.3% wrong, FAIL; the rule works on the reference individual only): kallisto long mode
  (`-k 63`, `bus --long --threshold 0.8`), reference = per-copy sequences + every other assembly transcript as OTHER,
  assign iff the compatibility class is one locus. On its molecules 919 right / 18 wrong (KB3781) and 678 / 6 (HG002);
  tool right / primary wrong vs reverse 15 vs 2 and 59 vs 2 (after the primary-overlap correction); where O2 abstains it
  assigns 672 / 6 and 332 / 6; reach 0.60 / 0.28, four-fold lower on HG002 reads of the non-reference haplotype.
  minimap2 on the same reference stays level with the primary: the gain is the exact-k-mer rule, not the reference.
- r1649 STOP at step 0: reads aligned to a per-family seqwish graph (GraphAligner), assigned by path compatibility: 637
  assigned vs kallisto 638 on 1,026 molecules, the 613 shared all on the same locus, 18 graph-only (bar 30).
- r1648 FAIL: oarfish 0.11.0 (own mapper + abundance EM, all region reads), locus posterior >= 0.99: equals the primary
  (24 vs 29, 36 vs 35) at a 9-11% wrong rate; where abundance alone decides it is 4 right / 16 wrong (KB3781) and
  95 / 13 (HG002). Closed.
- r1650 PASS (reach and agreement only, no read truth): kallisto rule on OR6737 assigns 2,259 of 5,971 contested
  molecules (0.378; KB3781 like-for-like 0.252), 3 of 212 disagree with O2; 13.5% of picks leave the primary's locus
  (KB3781 3.4%), unjudgeable without truth.
Scoring lessons (apply to every placement row): (1) a pick whose span overlaps the truth copy is right
(`EXON_PLACE_RIGHT=locus`); HG002 has copies of different families on one span; (2) `strict_truth.py` judges the primary
by its leftmost position: recheck with the record's overlap (`EXON_PLACE_PRIMARY`) before quoting primary-wrong counts;
(3) a catalog-only reference force-places outside-catalog reads; (4) validity gates check mechanics, reach is a bar.
- r1651 FAIL: the k-mer rule on HG002 reads against CHM13 (A119b catalog, HG002 haplotype truth): reach 0.043, 166 right /
  15 wrong (8.3%) vs the primary 163 / 18, 175 of 181 picks on the primary's locus. Reference-individual only.
  **O2 across individuals, registered as an O2 error (user decision 2026-11-12):** `--tie-ms 10` on that set is 170
  right / 175 wrong. Moves off the primary's locus 2 right / 121 wrong; kept on it 168 / 54. Three families hold 171
  of the wrong rows (MCL189, MCL601, MCL223); the other nine are 138 / 4. Next O2 rows: a guard on cross-individual
  moves, a per-family fit check, the same breakdown on same-individual samples. (Not excluded: the truth placing a
  CHM13 copy on its paralog's HG002 locus.)
- r1652 (no new run): O2's moves off the primary's locus on same-individual samples: KB3781 2 right / 2 wrong (4 of 570
  assignments, underpowered), HG002-on-pat 1 / 22 (23 of 525, 20 in MCL338); pooled with CHM13 5 right / 145 wrong.
  By the rule the guard is for every sample. O2's contribution is assign-vs-abstain on the primary, not relocation.
  Next row: never-off-the-primary guard (identity, held-out, rows removed on A119b).
- r1653 FAIL: that guard on held-out families (60 each): removed rows 103 wrong / 39 right (bar 4 to 1); wrong rate
  6.1% -> 4.0% (HG002 own haplotype) and 12.3% -> 5.8% (CHM13); correct lost 1.6% and 2.3%. On the A119b default gate it
  would delete 83% of assignments (exact ties: the primary is arbitrary). Not adopted, not a general guard. Moved
  near-tie rows are wrong 59-78% vs kept 4-6%; a few families carry the wrong moves (MCL1256, MCL338, MCL586, MCL189).
- r1656 PASS (no new run, held-out sets of r1653): the families carrying O2's wrong near-tie assignments are those with
  SOLE-CANDIDATE rows (contested molecule, one candidate copy, no deciding site; register 787's class). Held out: 19
  right / 15 wrong and 13 / 71, against 953 / 48 and 997 / 70 for other rows; dropping them takes the wrong rate 6.1% ->
  4.8% and 12.3% -> 6.6%. They are 13% of the A119b default-gate assignments. Existing flag `--no-sole-candidate`.
  USER RULING 2026-11-13 ("Make it default"): under `--tie-ms` a single-candidate molecule is now `tied` unless
  `--sole-candidate` (`copy_assign`, `sole_on`). Replay of the r1657 HG002 own-haplotype panel (workspace
  `sole_default_2026-11-13`): new default byte-identical to `--no-sole-candidate`; 758 vs 812 assigned rows; slack 0
  140 right / 16 wrong (10.3%) vs 160 / 26 (14.0%) with `--sole-candidate`. Figures of rows 1509-1657 run under
  `--tie-ms` need `--sole-candidate` to reproduce. `--union-certificate` is affected too: same panel 754 assigned
  rows under the new default vs 803 with `--sole-candidate` (803 = the r1657 run; union table unchanged, 8,916 rows),
  so r1657's union figures (170 / 14) included sole rows; the union arm was not re-scored under the new default.
  SECOND RULING 2026-11-13 ("Extend to the default gate"): the promotion is now opt-in everywhere (`copy_assign
  --sole-candidate`; driver `RUSTLE_SOLE_CANDIDATE=1`); `--no-sole-candidate` is accepted and does nothing. Evidence:
  r1695 (HG002 own haplotype, default gate) sole rows 17 right / 4 wrong, other rows 5 / 0: thin, one sample; r1656 had
  found sole rows clean on the reference individual (KB3781 33 / 0, 8 / 0), so this ruling trades reach for safety away
  from the reference individual. Cost on A119b (soft mask + rDNA lists, same binary, `sole_default_2026-11-13/asg_mr`
  vs `asg_mrkeep`): assigned rows 5,784 -> 3,555; the only rows that change are the 2,229 sole rows (assigned -> tied);
  assigned molecules 4,109 -> 3,322. Plain default run: 5,399 held 2,159 sole rows (about 3,240 left; not rerun).
  A certificate-rejected single candidate keeps its `ambiguous` verdict (a first build relabelled 31,447 of them
  `tied`; fixed before this record). Not re-measured: `quant.tsv` (differs, sole rows no longer count 1) and the flag
  scan, which reads the assignments table. Earlier figures need the switch to reproduce.
  Kept parked by the user the same day: `family_gfa`, pending item 1b (E7 on HG002), the inverted-copies simulation.
- r1656 follow-up: KB3781 held-out sole rows 8 right / 0 wrong (clean and rare on the reference individual). 787 was
  shipped only as `tie_outside_catalog` (partner on no family's unit); a partner on ANOTHER family's unit still leaves a
  sole candidate, and that is 46 of 71 and 12 of 15 of the wrong sole rows. Lever: cross-family tied placements as
  candidates (787's rule in full), or `--no-sole-candidate`.
- r1657 cross-family candidates: step 0 STOP for a demote-only guard (a tied placement on another family's unit is on
  112 of 118 sole rows, right and wrong alike = `--no-sole-candidate`). `--union-certificate` under `--tie-ms 10`, fresh
  120-family draws: own haplotype 170 / 14 vs flag 140 / 16 (better), CHM13 108 / 105 vs 81 / 66 (FAIL, MCL720 adds
  17 / 31). FAIL as registered; flag stays the recommendation. The third CHM13 draw is 45% wrong WITHOUT sole rows
  (MCL598, MCL1339): cross-individual error in thin families is a wider class than sole candidates.
- r1659 truth audit (post hoc): across individuals `strict_truth`'s near-perfect gate on catalog copies drops a copy
  > 0.3% from the HG002 haplotype and a sibling inherits its locus. Artifacts: 81 of r1651's 175 wrong rows, 27 of 141
  and 35 of 75 on the CHM13 draws. Corrected r1651: 251 / 94 (wrong 0.27). r1656 PASS and r1657 FAIL both stand.
  Real residual example MCL1339: HG002 allele 3-4 edits from both CHM13 copies, sole-candidate pick goes to the wrong one.
Standing: no default changes. Every alternative to the primary is now closed or reference-individual only: exon-space
minimap2 (1642-1644), abundance EM (1648), graph path compatibility (1649), exact k-mers (1647 self only, 1651).
- Two measurements behind user questions the same day (2026-11-12, no register row, descriptive):
  - **Alignment length is not the placement problem.** KB3781, the 1,570 decidable contested molecules: the primary
    aligns a median 99.9% of the read whether it is on the right copy (1,379) or the wrong one (191); under 90% aligned
    on 1.1% vs 6.8% (about 13 of the 191). An aligner that returns longer alignments has nothing to add here. BLAST's
    longer hits (edge-rule step 0, 2026-11-07) were gene vs gene at ~85% identity, a different regime.
  - **Secondary and supplementary records.** OR6737 `GGO_mm.bam`: 10,710,883 records = 4.41 M primary, 6,266,238
    secondary (58%), 34,699 supplementary (0.3%). Supplementary records are skipped by `as_table`, `candidate_copies`,
    `semisec` and the assembler's loaders; dropping them saves nothing. Secondaries are set by `-N 50 -p 0.1`; r1511
    stands (useful share 71% at ms ratio >= 0.99, under 7% below 0.8; `-p 0.8` loses 23% of the useful ones). Decision:
    keep the alignment line, thin by ms ratio AFTER alignment when a smaller BAM is wanted; never re-align stricter.

## O2 AUDIT AGAINST THE ASSEMBLY GOAL (2026-11-13, register 1699, measurement only)

Goal audited (user, 2026-11-13): on AS or ms near ties the primary is a coin toss; O2 should re-pick it on PSVs or
shrink the pool of secondaries, and that should help assembly. State audited: working tree of 2026-11-13 (sole rows
tied under `--tie-ms`, uncommitted); A119b default products `r1690_softmask/asg_bed` (exact-AS gate, famreps, censat
skip list); KB3781 default `r1662_edgemiss/kb_all`; truth sets r1642 (KB3781), r1653 hgpat / hg13 (hg13 with the r1659
correction), all `--tie-ms 10`. Tool `bench/o2_contested/o2_funnel.py`; outputs
`/anvme/workspace/v138da13-rustle/o2_audit_2026-11-13/`. No default changed, no assemble run.

1. **Three near-tie definitions.** A119b molecules: seeding pool (second AS >= 0.98 x best) 532,876; goal scope
   (ms within 10) 306,295; shipped O2 gate (exact AS tie) 111,796. `copy_assignment_definition.md` states all three
   (header 10 AS, current state 10 ms, shipped exact). The default gate sees 21% of the pool's molecules.
2. **Scope.** Of the 111,796 exact ties 64,157 never reach the assign table: 51,129 in families on the skip list
   (satellite / rDNA, intended), 10,863 with primary and rival outside every family copy, 2,165 other.
3. **What O2 says.** Contested rows, A119b: origin-rejected 148,260 (72%), tied 44,793 (22%), other ambiguous 8,523,
   assigned 5,293 (2.6%, 2,159 of them sole candidates with no PSV; row 1695). KB3781: 49% / 50% / 0.6% / 138 rows.
   `ambiguous` is the origin certificate in 95-99% of rows.
4. **Premise, per sample (never pooled; the two HG002 sets are the same reads on two references).** Exact tie on AS
   and ms, decidable rows: primary right 0.42 on KB3781 (n 147); O2 assigns 1 of them on PSVs. HG002 gives 0.27
   (hgpat, 532) and 0.37 (hg13, 383), but with the leftmost-position rule that under-counts the primary there: not
   quoted as rates. Truth-free check on A119b exact ties: O2's PSV pick is the primary's locus in 953 of 2,903.
   AS tie with ms 1-10 (nearest to the shipped gate): KB3781 primary right 0.55 (n 22), O2 6 right / 7 wrong;
   hgpat 32 / 0; hg13 21 / 0. AS 1-10 and ms 1-10: KB3781 primary right 0.93, O2 right where the primary is wrong
   9, the reverse 2; hgpat 37 / 8; hg13 43 / 10. The primary is a coin toss only on exact ties, and there no site
   separates the copies (the K = 0 wall); where O2 has sites the primary is mostly right already.
5. **Output object.** One copy or abstention per (molecule, family). The pairwise certificates exist in the code and
   are not written, so no placement can be marked excluded while the molecule stays undecided (rows with >= 3
   candidates: tied 19,173, ambiguous with >= 2 decisive sites 24,184).
6. **No path into assembly.** `assemble` runs before `assign`; nothing reads the assignments. A pool secondary counts
   as a read where it lands (r1518: USP17L 3,567 "reads" over 130 primaries), so a corrected primary adds no support
   at its destination; it can only remove support elsewhere. Seeding new loci is not reachable through O2: its
   candidates are loci already assembled (r1523 0 of 3; r1671-1673).
7. **The shipped gate picks the stratum with the least evidence.** Through `copy_assign` as shipped, PSV-assigned
   molecules are 2,995 = 0.56% of pool molecules (96 families; 2,719 of 3,240 rows in 10 families; 428 of 771 copies
   could lose half their reads or fall under 2 in the worst case). The pairwise certificate over the whole pool
   (2026-10-12 population study, same 532,876 molecules, script `s1_o2.py` not `copy_assign`, no family catalog)
   gave confirmed 62.5%, moved 2.8% (14,945), identical 10.6%, undecided 24.2%; on true ties 3.2 / 1.6 / 48.1 /
   47.1%. So on about two thirds of the pool a certificate exists and it mostly confirms the primary: the reachable
   product is a smaller pool (the certified-losing secondaries of ~333,000 molecules), not many new primaries.
   KB3781 whole sample under `--tie-ms 10` without the primary-on-family guard (`o2mstie/all10`): 53,132 ms-near-tie
   molecules, 4,237 PSV-assigned (8.0%), tied 46%, ambiguous 25%, not in the table 21%; the guard later removed
   about three quarters of held-out assignments (r1514).
8. **Origin-rejected reads stay in the pool.** The origin certificate rejects 68% (KB3781) / 72% (hgpat) / 74% (hg13)
   of reads whose origin is outside the catalog and 10% / 41% / 34% of reads whose true copy is in it. On KB3781 the
   primary is wrong on 35% of the in-catalog rejected (57 of 164). These reads count at every near-tied placement.
Direct tests of the end claim so far: the 2026-10-12 pruned good-pool arm (gene-level rule: gorilla stratum precision
up, human lost every S1 match and 91 S2 matches) and the SMN showcase (FAIL, O2 moves 8 reads against 616 in truth).
Read-level removal of certified-losing records has not been run. The moved reads of the population study mostly land
where the tools already assemble the gene (25 / 44 genes with >= 2 moved reads and no StringTie or FLAIR model).
Conflict for the user: "choose a new primary" against the standing prohibition (never claim O2 places better than
minimap2) and the 2026-11-03 statement that nothing decides exact ties. Items 4 and 7 give the numbers.
Proposed, each needing a PREREG: (a) O2's gate = the pool (AS >= 0.98 x best, one definition for pool, gate and tie
groups), with per-placement verdicts (keep / excluded / undecided) as the output; (b) panel assemble arm on KB3781
with the excluded records removed from the pool, scored on new against lost chains; (c) the same for origin-rejected
reads, separately.

Follow-up rows (2026-11-14): pool reduction by placement verdicts passed the KB3781 panel (register 1701) and FAILED held out
on HG002 (register 1704: in the band AS >= 0.98 x best, 4.1% of exclusions remove the true placement, bar 2%); not
ported. Planted-case simulation (register 1705): the only planted class that failed is a converted allele.

## ASSIGN SKIP LISTS: SATELLITE-ONLY BED, SIDECAR, rDNA ON HG002 (2026-11-15..16, register 1750, 1758)

The censat BED was replaced by a satellite-only list in the assign route (register 1750; `docs/LEDGER_O1.md`, last section):
families in the transition regions (NBPF, POTE, GOLGA, TSPY) are assigned again, A119b 5,399 -> 6,707 assigned rows. `assign`
reads `FASTA.exclude.bed` when `RUSTLE_EXCLUDE_BED` is unset (empty = off); CHM13 has that sidecar. Older assign figures used
the full BED.

rDNA skip (`RUSTLE_EXCLUDE_FASTA=tools/data/rdna_human.fa`) as a default: NOT switched (register 1758, Amendment 3 of
`PREREG_rdna_sequence_rule_2026-11-13.md`). HG002 on its paternal assembly, both arms on one 1 TB node: list identical to the
tool's (20 families), other families 2,898 against 2,898 assigned rows, wall 1,811 s against 3,557 s, largest process 30.4 GB
against 417.7 GB. One bar failed: the listed families hold 51 of the default's 2,949 assigned rows (1.73%, bar 1%), 50 of them
in one 249-copy rDNA array family. Whether rows inside rDNA arrays are a product is the user's ruling.

