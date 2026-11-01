# Ledger: general (handoffs, pending lists, register drafts, overviews)

Compacted 2026-11-01. Merged 2026-10-24 from 14 dated docs plus later pending notes, oldest first. Every original `## ` heading is
kept below; each entry ends with a pointer to the verbatim text in `docs/archive/2026-11/LEDGER_GENERAL.md` (and the originals in
`docs/archive/2026-10/merged/`). Numbers, verdicts and rows are copied from the original; anything not here is in the archive.

## Current state (2026-11-01)

- This ledger holds no live decisions of its own: handoffs, pending lists and register drafts up to 2026-10-27. Decisions live in
  `docs/NEGATIVE_RESULTS_REGISTER.md`, `docs/MODULE_STATUS.md` and the per-objective ledgers (`LEDGER_O1/O2/O3.md`).
- Defaults named here that still hold: assembler seeds with tied secondaries (rows 1060/1100/1101); `--union-certificate` opt-in
  (row 1103); `o3_candidates` stage opt-in (rows 1221-1231, 1247); most-junctions representative opt-in (rows 1234-1236).
- Most recent rows (register 1469-1484, 2026-10): all family rules FAIL or stay opt-in (1469-1482); the multi-copy assembly mode
  (r1479) FAILED on HG002 and A119b; `--tie-ms` stays opt-in (r1483); r1484 adds the PSV VCF report (`bench/o3_graph/psv_vcf.py`).
  Family-stage defaults are unchanged.
- Open from this ledger: the 2026-10-27 pending list (`## PENDING_2026-10-27`) and the unfinished `family_gfa`
  (`## PENDING (2026-10-26): unfinished family_gfa`); quantification parity bar (item 1 of 09-23) was measured as r1088 "close".

## Contents
Section index (original order): REGISTER_DRAFTS_machine2.md, REGISTER_DRAFTS_o3_maternal.md, REGISTER_DRAFTS_unmapped_rescue.md,
PENDING_2026-09-23.md, ADVISOR_HIERARCHY_2026-10-03.md, HANDOFF_2026-10-04.md, PENDING_2026-10-04.md,
TERMINOLOGY_FAMILY_SD_DUPLICON_EXPANSION_2026-10-07.md, PENDING_2026-10-11.md, OPEN_IMPLEMENTATIONS_PHASED_2026-10-12.md,
HANDOFF_laptop_to_main_2026-10-15.md, O1O2O3_FAST_REALISTIC_2026-10-16.md, O1O2O3_INTEGRATED_OVERVIEW_2026-10-16.md,
HANDOFF_cloud_2026-10-21.md, then the 10-26 and 10-27 pending notes. (full: archive/2026-11/LEDGER_GENERAL.md:6)

## REGISTER_DRAFTS_machine2.md
Origin `docs/archive/2026-10/merged/REGISTER_DRAFTS_machine2.md`. (full: archive/2026-11/LEDGER_GENERAL.md:26)

## Register drafts from the machine2 clone (renumber on merge into docs/NEGATIVE_RESULTS_REGISTER.md)
Draft rows M2-1..M2-9 (2026-09-30, Soto audits):
- M2-1 Soto nesting depends on the edge set: 440/444 (99.1%) on exon map-back edges, 346/444 (77.9%) on SEDEF-projected edges; quote 99.1% with its edge set.
- M2-2 biotype audit reproduces only counting 114 `Unassigned_*` singletons; on 491 real families 183 (37.3%) all-pseudogene.
- M2-3 fragment audit (147/420) does not reproduce; quote 154/491 (31.4%) holding a < 20% member.
- M2-4 Soto's finer splits not always backed by stable CN gaps: 40/47 single-exclusive-member families share >= 98% exons with another.
- M2-5 per-family labels vs our 0.9277 run: match 411, Soto smaller 37, we miss (our rule) 24, we miss (Soto's table) 5.
- M2-6 Soto families do not survive removing pseudogenes (313/491 gone).
- M2-7 loosening our CN cut (MAD < 1) loses on DEV and HELD-OUT at every looser cut (held-out exact 234 -> 230 at < 1.25).
- M2-8 Hungarian matching: sensitivity .935 / precision .975 (micro).
- M2-9 July "RNA missed" filter showed the 313 FOUND rows; fixed only in the merged page. (full: archive/2026-11/LEDGER_GENERAL.md:30)

## REGISTER_DRAFTS_o3_maternal.md
Origin `docs/archive/2026-10/merged/REGISTER_DRAFTS_o3_maternal.md`. (full: archive/2026-11/LEDGER_GENERAL.md:45)

## Proposed register rows — maternal / paternal reference study (2026-10-08)
Drafts (KB3781 Iso-Seq on its own `mat`/`pat` assemblies, 34 families + LRPAP1; source `docs/O3_MATERNAL_REFERENCE_2026-10-08.md`).
Of 959 reads unmapped on the primary, 132 map on `mat`, 3 on `pat`; 125 of the 132 in one 71-kb stretch of mat chr12. Prediction S1
(lower median `de` on the haplotype with the copy) FAILS at 4 of 8 (`mat`) and 11 of 13 (`pat`) large copies. LRPAP1 on the mother:
only the chrY copy absent. Recovery chain from `mat`: IsoCon R1 PASS, R2 PASS, R3 FAIL; from `pat`: R1 NOT TESTABLE, R2 PASS, R3 PASS.
(full: archive/2026-11/LEDGER_GENERAL.md:49)

## REGISTER_DRAFTS_unmapped_rescue.md
Origin `docs/archive/2026-10/merged/REGISTER_DRAFTS_unmapped_rescue.md`. (full: archive/2026-11/LEDGER_GENERAL.md:65)

## Proposed register rows — pool-first net for unmapped reads (2026-10-08)
Drafts, partly corrected by the prereg Errata. Key verdicts: U1 PASSES (unmapped reads of an erased copy form piles, 5,083 of 5,313 on
held-out); mmseqs2 attribution FAILS (held-out 10.5% wrong), replaced by dc-megablast cover score (4,255 of 5,312 rescued, 80.1%);
U2 second clause FAILS (pooled 406 vs single-read 410); consensus identity x coverage >= 0.999 only 14 of 40 (H); real mat chr12 bed
115 of 125 in one cluster at identity 1.0000; HMM profile 4,522 rescued / 2.69% wrong (exploratory); hybrid net NOT met (extra work
32.1% vs bar 25%); reconstructed loci as extra reference 2,505/2,657 (94.3%); IsoCon-style correction NOT met (14 -> 15 of 40, bar >= 19);
synthetic world recovery holds 30/30; 5' G-run handling (T, T3, T4, gate-aware metric: E3 12 -> 23/24); partition PART 24/24;
"consensus unsupported" flag; chain-aware clustering 24/24 vs 12/24; specificity guard; COPY/ALLELE/PRESENT call; discovery from 959
reads, DNA k-mer verification (113-read sequence 97.1% k-mers >= 2); four more libraries. (full: archive/2026-11/LEDGER_GENERAL.md:69)

## PENDING_2026-09-23.md
Origin `docs/archive/2026-10/merged/PENDING_2026-09-23.md`. (full: archive/2026-11/LEDGER_GENERAL.md:105)

## Pending work — registered 2026-09-23 (user, 15:33), to be pre-registered and run later
1. Quantification parity with StringTie/FLAIR/isoseq on shared chains (bar suggestion rho >= 0.85 on >= 1,000 shared chains per
   species); measured as r1088: close, use `--keep-coordinate-duplicates`.
2. O3 reference bias (advisor's claim: missing copies hide as unmapped/multi-mapped). Arms A simulation, B same- vs different-individual,
   C cross-species; SUPPORTED if >= 20% of an extra copy's reads unmapped/tied at <= 2% divergence, REFUTED if >= 80% absorbed. Arm A
   measured r1089: claim refuted with truth.
3. DNA step for the 52 RNA-only KB3781 candidates (PSV 31-mer depth vs host; present if >= 0.5x at >= 80% sites); plus HG002 as the human
   matched-individual test.
4. Decision for advisor: `--union-certificate` (opt-in, row 1103) removes every foreign claim (844 -> 0) but ties every MAPQ-0 read;
   assigning them needs an expression prior (thesis rule: no priors).
5. Seeding default flipped to tied secondaries (rows 1060/1100/1101): +30 loci on gorilla NC_073244.2, referee pairs 18 -> 36 of 108.
(full: archive/2026-11/LEDGER_GENERAL.md:109)

## ADVISOR_HIERARCHY_2026-10-03.md
Origin `docs/archive/2026-10/merged/ADVISOR_HIERARCHY_2026-10-03.md`. (full: archive/2026-11/LEDGER_GENERAL.md:230)

## "Where is my hierarchy in the MCL?" — what the advisor most likely means, and what exists (2026-10-03)
MCL (inflation 2.8) gives one flat cut; families nest (family > subfamily > copies > alleles). Answer: the cut is one level of the
nested edge-test lattice L0-L4 (L3 0.985 / L4 0.995 shipped as thresholds); average linkage recovers NPIPA/NPIPB exactly (register 457)
but merge heights add nothing beyond identity. Do not re-propose: subfamily level as repair (§5k), broad/recent two-tier (rows 9, 311),
inflation sweeps, superfamily components (row 654). (full: archive/2026-11/LEDGER_GENERAL.md:234)

## HANDOFF_2026-10-04.md
Origin `docs/archive/2026-10/merged/HANDOFF_2026-10-04.md`. (full: archive/2026-11/LEDGER_GENERAL.md:283)

## Handoff — 2026-10-04 (session of 10-03/10-04), for the next agent
Week priority: prove O1 and O2. Done: `o3_candidates` A13 passed / A14 failed, stays OPT-IN (rows 1221-1231); strict "found copy" rule
= Amendment E (capped starts; NPIP own node 23/21/24 -> FOUND 10/6/2; rows 1232-1244); most-junctions representative stays OPT-IN
(rows 1234-1236); LRPAP1 dossier: 11 loci in gorilla (8 full + 3 fragments), O2 0 contested reads (rows 1245-1246). Last register row
then 1246. (full: archive/2026-11/LEDGER_GENERAL.md:287)

## PENDING_2026-10-04.md
Origin `docs/archive/2026-10/merged/PENDING_2026-10-04.md`. (full: archive/2026-11/LEDGER_GENERAL.md:370)

## Pending — 2026-10-04 (user decision: this week is for proving O1 and O2 work right; O3's stage work is parked here)
Amendment 15/15b DONE 2026-10-05 (register 1247): C1' (10/53) and H1 (6/30) fail, stage STAYS OPT-IN. Open then: row-1220 O2
`--tie-outside-at-best` prereg; HG002 held-out; chapter drafts; representative-rule follow-ups (f1v2 re-test, support-aware rule);
NPIP labeller marks; LRPAP1 follow-ups; inversion-vs-spliced-aligner prereg NOT RUN. (full: archive/2026-11/LEDGER_GENERAL.md:374)

## TERMINOLOGY_FAMILY_SD_DUPLICON_EXPANSION_2026-10-07.md
Origin `docs/archive/2026-10/merged/TERMINOLOGY_FAMILY_SD_DUPLICON_EXPANSION_2026-10-07.md`. (full: archive/2026-11/LEDGER_GENERAL.md:432)

## Terminology: multi-copy gene family vs. segmental duplication, duplicon and expansion (2026-10-07)
One object, four levels: multi-copy gene family (the object O1 defines); SD (>= 1 kb, >= 90% identity; edge evidence for young copies);
duplicon (sub-gene node class); expansion (lineage-specific event/clade). SD-embedded families are a subset of multi-copy families;
expansion and family are not 1:1. Concede the noun to the advisor; letting SD/duplicon decide membership failed every time (rows 676,
677, 817, 818, 997). Never pool Soto (human) and Yoo (gorilla) numbers. (full: archive/2026-11/LEDGER_GENERAL.md:436)

## PENDING_2026-10-11.md
Origin `docs/archive/2026-10/merged/PENDING_2026-10-11.md`. (full: archive/2026-11/LEDGER_GENERAL.md:511)

## Pending work, 2026-10-11
Genome-wide rescoring paused until the figure set is agreed with the advisor (Claude Doc "Benchmark Figure Decisions"). Done: first
genome-wide `families` runs (A119b 2,150 families / 9,907 loci; human testis 342 / 1,115; chimp 396 / 1,343; KB3781 326 / 918), assign
chimp + KB3781, figure 8 tables. Paused: figure 7 genome-wide, assign A119b, O3 flag scans.
- StringTie milestone (registers 1264-1266): `--direct-transfrags-union` opt-in beats StringTie on intron/chain/transcript/locus in human,
  costs 0.8-1.5 chain-precision points vs default; StringTie leads base sensitivity (gorilla 50.3 vs 46.9). Remaining gap only 9 (5%)
  flow-only. Opt-ins built (descriptive, not pre-registered): `--extend-transfrags`, `--transfer-ends`; StringTie's `-E 25` snap costs
  3.2 / 2.6 chain-sensitivity points; `union + --read-isoform-k 2` reaches 63 / 133 StringTie-exclusive gorilla matches but costs 4.6 / 2.3
  chain-precision points (B2 bar 0.5; same verdict as 848 / 855): not a default.
- FLAIR milestone: FLAIR matches 89 / 147 transcripts the union does not (union 309 / 218 the other way); at most ~41 / 32 reachable in
  principle. FLAIR-lever opt-ins (`--transfer-credit full|best`, `--read-junction-snap`, `--min-chain-reads`, `--check-splice`,
  `--terminal-exon-anchor`): full credit / single-read chains reach about half of FLAIR's exclusives at FLAIR's precision; read filters
  nearly inert. None is a default.
- Multi-copy loci vs StringTie/FLAIR (registers 1267-1272): `--direct-transfrags-union --read-isoform-k 2` above both on chain sensitivity
  (1268); secondary-only genes not a win (1269); after O2 pruning S2 genes 72 vs 51 / 42 (gorilla) (1270); O2 vs Eichler AS rule agree on
  14-54%; O2 moves 1.0 / 2.8% of multi-placement reads (1271-1272).
- O2 follow-ups 10-13 (register 1274-1275): HG002 haplotype-truth arm; per-family blind spot; re-run O2 numbers on a binary >= 7ac557e7.
- O1 two-node fix: default `mcl_families --pair-rule exact-mcl`, `weight-blind` opt-in (user, 2026-10-13).
- sedef-rs layers: V1 FAIL (1278), V1b FAIL (1279), V2 FAIL (1280; V2a held-out P1 0.5153 vs 0.80; V2b PASS). Plan to make sd-scan replace
  SEDEF (recall ~0.40, precision 0.92-0.95) and DupMasker (Route B with `dupliconlib.fa` recommended); parked, nothing run.
(full: archive/2026-11/LEDGER_GENERAL.md:515)

## OPEN_IMPLEMENTATIONS_PHASED_2026-10-12.md
Origin `docs/archive/2026-10/merged/OPEN_IMPLEMENTATIONS_PHASED_2026-10-12.md`. (full: archive/2026-11/LEDGER_GENERAL.md:849)

## Open implementations — phased data (2026-10-12)
Toy (`docs/PHASED_TOY_2026-10-12.md`, P0-P3 PASS): with two haplotype assemblies a truth-free rule separates allele from extra copy 100%
in the easy case. Nothing built; step 0 = ask which form (assemblies / phased VCF / haplotagged BAM). Items: I1 `--hap-ref` two-assembly
mode (M), I2 `--hap-tag HP` label only (S), I3 VCF via whatshap (XS), I4 rename `--phase` (HP = copy, not haplotype) (XS), I5 harder toy
rungs. Limits: copy vs allele not settled on RNA alone; phased data must be of the same individual. (full: archive/2026-11/LEDGER_GENERAL.md:853)

## HANDOFF_laptop_to_main_2026-10-15.md
Origin `docs/archive/2026-10/merged/HANDOFF_laptop_to_main_2026-10-15.md`. (full: archive/2026-11/LEDGER_GENERAL.md:936)

## Handoff — laptop → main machine (2026-10-15)
Branch `laptop/o3-psv-multimappers`: pixi envs, O3 positive control (PASS), few-read false flags diagnosed, O3 per the advisor's
definition (PSV clusters of O2-placed reads; a cluster no copy explains = absent copy), Rust port `copy_assign --psv-clusters` (OPT-IN).
Done on main 2026-10-15: merged (`integrate/o3-psv-multimappers`, merge `90df29df`), fixtures committed, suite 1,356 passed / 0 failed,
`copy_assign` byte-identical without `--psv-clusters` on real data; register rows 1289-1293. Still open then: M3 foreign-species check,
real-data O3 prereg (remove copies from the GENOME, single-copy control, run with `--no-as-tied-only`). (full: archive/2026-11/LEDGER_GENERAL.md:940)

## O1O2O3_FAST_REALISTIC_2026-10-16.md
Origin `docs/archive/2026-10/merged/O1O2O3_FAST_REALISTIC_2026-10-16.md`. (full: archive/2026-11/LEDGER_GENERAL.md:1069)

## Faster family tests, realistic reads, and what O2's evidence rests on (2026-10-16)
Update 2026-10-19: O2's domain is AS-tied multimapping reads only (rival within 10 AS, Eichler's rule). Development on family sets 8-10
toward the integrated O1->O2->O3 pipeline (`bench/o3_graph`), subsections 1-18 (faster tests, realistic reads, PSVs vs SUNs, CIGAR/NM/de,
v2..v14, independent review: section 9 wins where it disagrees). Round verdicts:
| round (row) | pipeline | verdict |
|---|---|---|
| 9 | v1 | FAIL both profiles |
| 11 (r1309) | v4 | PASS, 12/12 |
| 12 (r1310) | v5 (`--o1-hosts`) | FAIL, both profiles |
| 13 (r1311) | v5' (`--branch-genome-novel`) | FAIL (A119b); public 6/6 |
| 16 (1318) | v12 (O2 margin) | FAIL correctness (4 of 6), PASS cost |
| 17 (1320) | v13 (`--struct-rescue --short-path-align`) | FAIL correctness (0 of 6; OUTSIDE B2), PASS cost |
| 18 (1322) | v14 (O1 `--ref-copies`, row 1321) | FAIL correctness (4 of 6; seed 803 OR4F), PASS cost |
| 19 (1324) | v14 unchanged, paralogs expressed as in testis | PASS correctness (6 of 6, all 72 bar checks) and cost |
Round 19: 85,440 correct vs 62,070, 2,355 errors vs 8,028 (current method). PSV answer (section 17): ~70% of decisive O2 columns are
SUN-type; edits priced by rate (`--o2-significance`), indels as events; de flags absent-copy reads. Section 18: the variation graph is
the model, linear alignment the implementation (`--branch-linear`; O3 158 s -> 34 s, same verdicts on 34 spent runs).
(full: archive/2026-11/LEDGER_GENERAL.md:1073)

## O1O2O3_INTEGRATED_OVERVIEW_2026-10-16.md
Origin `docs/archive/2026-10/merged/O1O2O3_INTEGRATED_OVERVIEW_2026-10-16.md`. (full: archive/2026-11/LEDGER_GENERAL.md:1972)

## O1 -> O2 -> O3 integrated: the overall picture (2026-10-16)
Post hoc summary over rounds 1-7 (register 1297-1305), no new claim. Integrated pipeline (O1 RNA-supported seeding, O2 padded union at
margin 2, O3 arm I) vs current method over 43 worlds: correct 95,405 vs 82,163 (x1.161), errors 445 vs 6,377 (-93%); not-tuned 30 worlds
x1.178 / -92%. Largest change from O3 (current method's 3,412 absent-copy calls on present copies' reads vs 61). Unfixable: near-identical
mapq-0 twins, alleles vs copies below delta, NPIPA.5 hosted on NPIPB at `de` 0.0097. Updates: round 8 FAIL by one read; round 9 (r1307)
FAIL, 10 (r1308) FAIL, 11 (r1309) v4 PASS 12/12; 12 (r1310) and 13 (r1311) FAIL; 14 (r1312) v9 PASS 12/12 with a diploid donor; the allele
confound remains the deciding open problem; every round is a simulation (delta 0.00958). (full: archive/2026-11/LEDGER_GENERAL.md:1976)

## HANDOFF_cloud_2026-10-21.md
Origin `docs/archive/2026-10/merged/HANDOFF_cloud_2026-10-21.md`. (full: archive/2026-11/LEDGER_GENERAL.md:2209)

## Handoff — local WSL machine → cloud (2026-10-21)
Branch `o3-fast` (283 commits past main, registers 1347-1369). State: 1360-1361 corrector UNDERPOWERED; 1362 FAIL; 1363 reference-free O3
on KB3781 NPIP FAIL on read recovery; 1364-1365 candidates added to the reference, reads return (33 / 39), 0 false moves; 1366 descriptive;
1367 multi-copy gains not shown to be real copies; 1368 O1 vs SD calls FAILS (72.1% / 82.4%); 1369 `--reject-self-hits` v1 FAIL, v2 PASSES
held-out gorilla_OR6737, A119b fails G1 (1,036 -> 1,140 > 1,088), made the DEFAULT for all catalogs (user decision). Working rules
(pre-register, audit, pair precision/recall only, never pool species, O2 = AS-tied within 10 AS, `--rule multigene`), data copy lists and
the family-region subset recipe are in the archive. (full: archive/2026-11/LEDGER_GENERAL.md:2213)

## PENDING (2026-10-26): unfinished family_gfa
`family_gfa` (graph view of O1 families) is Task 1 of 6 only (`src/family_gfa.rs`, tests pass on main); tasks 2-6 not written. Spec
`docs/superpowers/specs/2026-09-30-family-gfa-design.md`, plan `docs/superpowers/plans/2026-09-30-family-gfa.md` (paths say `src/bin/`).
(full: archive/2026-11/LEDGER_GENERAL.md:2351)

## PENDING_2026-10-27
(full: archive/2026-11/LEDGER_GENERAL.md:2365)

## Pending work — registered 2026-10-27 (user), to be pre-registered after the DENOVOFAM arms (`docs/LEDGER_O1.md`)
1. Liftoff-style family mode: exon-level placement + chaining of read-supported templates with Liftoff defaults (near misses r1286 recall
   0.638 < 0.80, r1287 0.887 < 0.90 never chained exons). 1c. Exon-only k-mer sketch families with high-copy k-mer removal (user
   2026-10-29). 1b. PAUSED: fresh-sample test of DENOVOFAM-L5 arm E7 on HG002 Kinnex (per L5 audit r1441). 2. Does the de novo fix carry
   into semi-guided (SEMISEC)? 3. Strainy-style clustering (label propagation / recursion) for O1 families.
   (Later outcome per register 1469-1482: the sketch tree passed HG002 but was not shipped; semi-guided PASS not default.)
(full: archive/2026-11/LEDGER_GENERAL.md:2367)
