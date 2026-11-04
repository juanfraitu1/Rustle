# Ledger: O1 (family definition, assembler, locus units)

Compacted 2026-11-01. One short entry per section, in the original order; every `## ` heading of the full ledger is kept verbatim. The verbatim full text is `docs/archive/2026-11/LEDGER_O1.md` (each entry ends with its line number there); git history has every version.

## Current state (2026-11-03)

**Default families stage (de novo O1), unchanged through register 1514:** `mcl_families --from-gtf --min-exonic-bp 1 --min-shared-exon-frac 0.60 --emit-units` on the gated minimap2 all-vs-all, MCL 2.8, `--pair-rule exact-mcl`, `--reject-self-hits`, most-reads representative. Families are built from `PREFIX.gtf` (unlayered); since r1498 the assembly also writes `PREFIX.final.gtf` with the sub-chain drop (default on, `RUSTLE_ASSEMBLY_LAYERS=off`), and `PREFIX.gtf` stays byte-identical.

**Assembly defaults changed 2026-11-01/02:** `--polish-tss rescue` (r1485 PASS, `RUSTLE_POLISH_TSS=off` restores); sub-chain drop in `final.gtf` only (r1498 PASS: c share 10.1 -> 5.4%, = -0.44%, chain precision 17.6 -> 19.1). Opt-in assembly layers that PASSED on their own held-out but FAILED the families non-regression (r1490 E6: A119b pair precision 0.726 -> 0.700, recall 0.741 -> 0.714, single-exon copies 3,299 -> 217; r1491 E7: recall 0.741 -> 0.710): `--polish-tss-end median` (r1487), `--polish-tes-extend on` (r1488), `--mono-model cluster` (r1489). E2 sub-chain rescue STOPPED at development (r1486).

**Annotation-side rules (opt-in, none default):** readthrough detection: labels mode default, structural rule FAIL bar 1 (r1492, GENCODE recall 0.757), gene-level labels only (r1493); truth split of CAT collapsed genes below bar (r1494, recall -0.0049); annotation-guided fused-loci cut PASS panel and genome-wide (r1495: fused 547 -> 282, precision 0.726 -> 0.761, recall -0.0096) = recommended opt-in; shared-exon split FAIL (r1496 precision 0.698 -> 0.344; r1497 private units same); SQANTI3 rules filter FAIL (r1499, = -6.4%); own end-artifact rules E/J/K FAIL the 0.5% = bar (r1500-1502; `utilities end-artifacts` shipped as a tool, no layer).

**Opt-in family tools (all failed or unvalidated as defaults):** `--compat-support`, `--direct-transfrags[-union]`, `--fold-split-bridges` (r1435/1436), lift edges U1 and variants (r1437-1441, r1453-1457, r1467-1469), `--shared-exon-columns`, `--chain-gate[-union]` (r1475), `--representative reads-x-junctions` (r1473), dispersion (r1470), SKETCHFAM (r1471), semi-guided (`semisec` + `semifam`, PASS both held-outs r1472, auditor first), LOCUS-ABSORB (r1476 PASS by equality), FAMASM (r1479 FAIL).

**Comparisons:** IsoCon WORSE in both classes (r1478); SQANTI single-copy chr20: ours FSM+ISM 30.2% / artifact 5.0%, StringTie, FLAIR, IsoCon worse (r1482); multi-sample pooling A119b + HG002 more sensitive, not more precise (r1503).

**Open:** families GTF from pre-layer transcripts is done (the SPLIT); a fresh long-read RNA sample for any new family rule (KB3781, OR6737, A119b, HG002 are spent); NEXT TARGETS single `rustle` binary with stage hatches; FAMASM only with a stricter certificate on A119b.

## Contents

Merged 2026-10-24 from 46 dated docs (originals in `docs/archive/2026-10/merged/`); later sections appended directly. Each `## <file>.md` heading below is the wrapper for one merged doc; the section(s) after it are that doc. (full: archive/2026-11/LEDGER_O1.md:6)

---

## PREREG_hierarchy_AMENDMENT2_DRAFT.md

(full: archive/2026-11/LEDGER_O1.md:58)

## Amendment 2 — DRAFT v2 for the committed prereg (coordinator, 2026-10-07 evening). NOT committed.

Draft Amendment 2 to `PREREG_hierarchy_duplicon_expansion_2026-10-07.md`: Phase E (tier B) cross-species Liftoff 1.6.3 of gorilla pool components onto H (CHM13), C (mPanTro3), O (mPonPyg2) with `-copies -sc 0.95 -a 0.5 -s 0.5 -overlap 0.1`. Plan: 180 Liftoff runs (42 `pos` LRPAP1 control, 6 `neg` planted-pair control, 132 `verdict`), ~83 min predicted, 2x time budget, 16 GB RSS kill. Verdict plan excludes the two control components (298 / 1,386 records instead of 300 / 1,403); positive control judged also on triple (gene-LRPAP1 x3). Class with any multiplicity-flagged record dropped and counted. Authorised by user "go E" 2026-10-07 17:57 PDT. Code hashes, inputs, stop rules and readings 10.1-10.11 in full text. The three in-code-block comment headings of the command recipe are:
## once, light: the plans and the wall budget (about 2 s; writes only plan.json), then the run GFFs and their marker (about 8 s)
## heavy, FOREGROUND, one call at a time; repeat immediately while the exit code is 75
## after exit 0, light: the table twice (Gate 3), each into a NEW directory
(full: archive/2026-11/LEDGER_O1.md:62; command block lines 151-156)

---

## IDEAL_CHROMOSOME_SIM_2026-09-21.md

(full: archive/2026-11/LEDGER_O1.md:200)

## Ideal-scenario chromosome simulation — in theory it does work

§6v1 (HUMAN chr16, 1,407 genes, 10 ideal reads/transcript; prereg `2a809b76`). Bars 1 and 3 FAIL as written (A_ideal 0.8067 overall / 0.6364 multi-copy / 0.8523 single-copy; bars mis-set, endpoint ceiling ~0.84 from 15.8% overlapping-annotation genes); bar 2 not triggered; bar 4 (attribution) PASSES. A_ideal: 12 genuine over-merges, 0 single-copy, corrected 0.9578. Readthrough injected at measured 7.51% -> 188 over-merges (real data 139); truncation costs nothing on this endpoint. Lever: readthrough-aware node split. (full: archive/2026-11/LEDGER_O1.md:204)

## 8. Does the idealized version hold a better FAMILY definition?

Shipped `mcl_families` on ideal/rt/real loci, chr16, vs Soto cover (15 fams) and protein referee (36). F (Soto / referee): GUIDED 0.4174 / 0.2248, A_ideal 0.3561 / 0.1944, REAL 0.3166 / 0.1805, A_rt 0.3068 / 0.1570. Yes, modestly. Bigger gap is node COVERAGE (Soto genes with a node: GUIDED 47/58, A_ideal 34/58); 5 of 13 assembly-attributable misses are a scoring artifact (RefSeq fusion records like PKD1P3-NPIPA1 take the locus name). (full: archive/2026-11/LEDGER_O1.md:291)

## 9. The 8 residual genes: three distinct mechanisms, not one

8 genes split: A no transcript model (HERC2P5, HERC2P8, NPIPB10P); B fusion naming (PKD1P6); C real paralogs: ABCC6 (both loci complete, exon structure differs, edge fails at shared-exon 0.6), NPIPB7/SMG1P6 genuinely fragmented (20.6% / 12.4% of transcript). Only NPIPB7 and SMG1P6 (2 of 8) are fixable assembly gaps; isoseq cluster2 a plausible, untested lever. (full: archive/2026-11/LEDGER_O1.md:376)

## 10. Root cause of the NPIPB7/SMG1P6 fragmentation: a NON-CANONICAL junction in RefSeq's own annotation

Reads and polish ruled out; RefSeq's own exon models carry non-canonical junctions (SMG1P6 intron 15 AT-AG; NPIPB7 introns 6-8) and the assembler truncates at the first non-canonical junction. `RUSTLE_JUNCTION_MAJORITY=1`: SMG1P6 full recovery, NPIPB7 partial. Same mechanism as §6m8; supersedes the cluster2 lead of §9. Default flip left as user's call. (full: archive/2026-11/LEDGER_O1.md:450)

## 11. Including non-canonical junctions on REAL chr16 data: measured, not assumed

`RUSTLE_JUNCTION_MAJORITY=1` on real chr16: transcripts 9,629 -> 10,093; family F flat (Soto 0.3166 -> 0.3184, referee 0.1805 -> 0.1817); Soto genes with a node 34/58 -> 35/58 (+NPIPB2, 0 lost). NPIPB7/SMG1P6 already covered on real data. Known cost (§6n4): engulfed copies 79 -> 86. Recommended; default flip was later committed (see §12). (full: archive/2026-11/LEDGER_O1.md:527)

## 12. The full remaining recall gap, classified (not just the original 8-gene sample)

After the default flip, 23 of 58 Soto chr16 genes lack a node: A no RefSeq transcript model 6 (26.1%), B bigger co-located record wins naming 8 (34.8%), C locus captures a fraction 9 (39.1%). NPIPB7/SMG1P6 no longer in the list. (full: archive/2026-11/LEDGER_O1.md:603)

## 13. Category C dissolves: most of it is the SAME naming-collision mechanism as Category B

5 of 9 C genes lose naming via a shared de novo locus (e.g. NPIPB4/RRN3P1 -> LOC112268174). Final tally of 23: A 6 (26.1%), B naming collision 13 (56.5%), CDR2 truth label not corroborated 1 (4.3%), D real construction/depth gap 3 (13.0%: ABCC6, ABCC6P1/2). The max-overlap one-name resolver in `score_family_def.py` / `cover_and_jn_definition.py` likely undercounts recall. **Correction (2026-10-06):** `bench/sim.py chromosome` dropped each transcript's leftmost exon (4,370 transcripts of 1,407 genes, none full length); every number of this document describes those reads; corrected ceilings not computed. (full: archive/2026-11/LEDGER_O1.md:643; correction :712)

---

## PREREG_inversion_aligner_blindness_2026-10-05.md

(full: archive/2026-11/LEDGER_O1.md:721)

## Inverted copies and the strand-aware spliced aligner — mechanism + a parked, pre-registered measurement (NOT RUN)

Claim: `splice:hq` is strand-aware, so wholly-inverted copies' reads go unmapped and are invisible to O1/O2; `o3_candidates` pass B (map-ont on unmapped reads) is the designed mitigation. Parked famsim experiment (control + `inv_whole`/`inv_exon`/`inv_intron`), claims C1-C4 registered. NOT RUN (user cost decision). (full: archive/2026-11/LEDGER_O1.md:725)

---

## DEFAULT_RESCORE_NPIP_2026-10-06.md

(full: archive/2026-11/LEDGER_O1.md:786)

## The current default at NPIP, re-scored at HEAD (2026-10-06)

DEV, A119b chr16, 25 NPIP copies. **PASS**: R1 strict FOUND 8/25 (bar >= 6), R2 locus-level 21/25 (>= 18), R3 U2 F 0.645 (bar 0.645, at the bar because HEAD byte-identical to 2026-09-30), R4 NPIPB2+NPIPB6 own nodes. DEF vs PRE: Compara CF153 F .690 -> .812. 14 of 22 spliced-expressed copies not found: the locus representative lacks the capped start. (full: archive/2026-11/LEDGER_O1.md:790)

---

## ENTANGLED_BASELINE_2026-10-06.md

(full: archive/2026-11/LEDGER_O1.md:845)

## What StringTie and FLAIR do at overlapping genes, on ideal reads, and where we stand (2026-10-06)

DEV, CHM13 ideal reads (NPIP/TBC1D3 windows, 2 reps). Chains at entangled genes (E): ours 1359/1606 (84.6%), FLAIR 1159 (72.2%), StringTie 1015 (63.2%); but default makes ~2x artifacts (1031 vs 412/431; precision 85.3% vs 92.3/91.5%). Arm PC (`--no-seed-secondaries` + `--polish-subchain drop`): 451 artifacts, precision 92.9%, recall 89.1%. **Registered bar: no arm BETTER** (PC fails (b) 451 > 412). P2 and P4 FAILED; 216/256 E gene-runs share an exact junction with another gene, so no junction-based gene_id rule separates them. (full: archive/2026-11/LEDGER_O1.md:849)

---

## FEWCOPY_IDEAL_CASES_2026-10-06.md

(full: archive/2026-11/LEDGER_O1.md:899)

## Few-copy ideal cases for O1, and are NPIP / TBC1D3 "unreachable"? (2026-10-06)

EXPLORATORY, DEV, A119b, stored products. NPIP/TBC1D3 not unreachable (ideal: TBC1D3 12/13 passes, NPIP 12/14 misses by one). Input criteria C1-C5 over 377 Compara 2-4-gene families leave 11; 4 exact in all three views. Recommended DEV regression anchors: GSTA1/GSTA2 (CF316, chr6) and PAGE2/PAGE2B/PAGE5 (CF400, chrX). Class level (292 scorable 2-4-gene families): sens .297, prec .935, F .450. LYZL1/2 and UGT2B15/17 fail with no admitted edge. (full: archive/2026-11/LEDGER_O1.md:903)

---

## IDEAL_EXPRESSION_DEFAULT_2026-10-06.md

(full: archive/2026-11/LEDGER_O1.md:992)

## If NPIP and TBC1D3 were ideally expressed, would the current default pipeline find them? (2026-10-06)

DEV, circular, register 1254. TBC1D3 **YES** (12/13 reachable, both reps, cluster precision .938); NPIP **NO** (12/14, one short of 13; .909). Perfect-loci control 14/14 and 13/13. Misses are representative choice (NPIPA2 artifact representative, NPIPA7 locus purity .112, TBC1D3P4/P3 coin-toss with paralog secondaries). Real vs ideal annotation-anchored FOUND: NPIP 12 -> 23, TBC1D3 12 -> 16: real shortfall mostly data. P1-P3 held. (full: archive/2026-11/LEDGER_O1.md:996)

---

## LOCUS_UNITS_2026-10-06.md

(full: archive/2026-11/LEDGER_O1.md:1127)

## A locus model for overlapping genes on ideal reads: primary-first representatives and structural separators (2026-10-06)

DEV, ideal reads, NPIP/TBC1D3 x 2 reps. Rule 2 (private-bridge structural separators, S_D) **refuted** (R found 38 vs D 48; cuts a gene's own main isoform). Rule 1 (primary-first representative, Q): E3 50 -> 54/54 but NPIP rep 2 cluster splits (Compara F .733 -> .643); Q1, Q3-Q6 failed: **no arm is a CANDIDATE**. Primaries-only P gives YES for both families (NPIP 14/13, TBC1D3 13/13) at cost of entangled copies (7 vs 12). E4 brittle to partition granularity (motivates hierarchy-aware criterion). (full: archive/2026-11/LEDGER_O1.md:1131)

---

## PREREG_default_rescore_npip_2026-10-06.md

(full: archive/2026-11/LEDGER_O1.md:1175)

## PREREG: strict FOUND, U2 and Compara CF153 for the current default pipeline at NPIP (2026-10-06)

DEV, A119b chr16, 25 NPIP copies; arms DEF (default at HEAD) and PRE (`--bridge-regroup off`, cov 0). Gates G0/G1/G3; rules R1 strict E-found in own nodes >= 6/25, R2 locus-level >= 18/25, R3 U2 F >= 0.645, R4 NPIPB2+NPIPB6 own nodes. Amendment 1: copy exons from `npip_read_pool.json`; G1 on P and GOOD. Result: PASS (see `## The current default at NPIP, re-scored at HEAD (2026-10-06)`). (full: archive/2026-11/LEDGER_O1.md:1179)

---

## PREREG_entangled_baseline_2026-10-06.md

(full: archive/2026-11/LEDGER_O1.md:1265)

## Pre-registration: what StringTie and FLAIR do at overlapping genes, on ideal reads, and what "better" means (2026-10-06)

Arms S (StringTie -L), F (FLAIR collapse), D_asm/D_fam, U (f1units). Predictions P1-P4; BETTER iff (a) chains at E >= max(S,F), (b) artifacts <= min(S,F), (c) resolved at E >= max(S,F) and > D_fam, (d) no loss at N. Result: no BETTER arm; P1, P3 held, P2, P4 failed (see `## What StringTie and FLAIR do at overlapping genes, on ideal reads, and where we stand (2026-10-06)`). (full: archive/2026-11/LEDGER_O1.md:1269)
Amendment 1 added lever arms P / C / PC with predictions Q1-Q4 (Q1, Q4 held; Q2, Q3 FAILED). Errata 2026-10-07 after verification wf_994ece76-8d2. (amendment: archive/2026-11/LEDGER_O1.md:1318)

---

## PREREG_fewcopy_class_DRAFT_2026-10-06.md

(full: archive/2026-11/LEDGER_O1.md:1339)

## DRAFT pre-registration: class-level few-copy O1 test on a held-out substrate (2026-10-06)

DRAFT v0, NOT REGISTERED, nothing run. Held-out class K = families passing C1-C5 (`bench/fewcopy/criteria.md`), primary substrate gorilla OR6737; E1 YES >= ceil(0.9 n) reachable exact, NO < ceil(0.7 n); E2 false merges on negatives <= annotation arm; E3 counterexamples >= 2 = NO; n_reachable < 8 UNDERPOWERED. Needs user decisions (substrate, 0.7 bar, committed runner). Superseded in practice by `## Gorilla OR6737 testis few-copy class test: design sweep result (2026-10-07)`. (full: archive/2026-11/LEDGER_O1.md:1343)

---

## PREREG_ideal_expression_2026-10-06.md

(full: archive/2026-11/LEDGER_O1.md:1369)

## PREREG: if every copy of NPIP and of TBC1D3 were ideally expressed, would the current default pipeline find them? (2026-10-06)

S1 substrate: +-500 kb windows (NPIP 25 copies chr16, TBC1D3 16 copies chr17), every transcript at 10 ideal reads, 2 replicates; legacy 09-21 reads NOT USED (leftmost-exon bug). Strata E/C/X/R; IDEAL-FOUND = E2 own locus and E3 exact representative and E4 in family; YES >= ceil(0.9 n) and CP >= 0.5, verdict = lower replicate; controls G1-G7. Predictions P1 TBC1D3 YES, P2 NPIP not YES, P3 replicates agree. Amendment 1 narrows the shared-span INVALID rule. Result: P1-P3 held (see `## If NPIP and TBC1D3 were ideally expressed, would the current default pipeline find them? (2026-10-06)`). (full: archive/2026-11/LEDGER_O1.md:1373)

---

## PREREG_locus_units_2026-10-06.md

(full: archive/2026-11/LEDGER_O1.md:1478)

## Pre-registration: a locus model for overlapping genes: primary-first representatives and structural separators (family level, ideal reads) (2026-10-06)

Rule 1 primary-first representative (arm Q), Rule 2 structural separators (S_D, S_Q); predictions Q1-Q6; CANDIDATE iff IDEAL-FOUND (R+E) > D, R >= D, E >= D-1, Compara F >= D-.005, both replicates. Amendment 1: Rule 2 = private bridge (38 separators, 16/108 resolved). Amendment 2 (after level-M seen): level C = graph components; arms F1, F1Q; H1-H5. Errata 2026-10-07 (H4 unmeetable; H5 rationale wrong for S_D). Result: no CANDIDATE at level M (see `## A locus model for overlapping genes on ideal reads: primary-first representatives and structural separators (2026-10-06)`); level C in `## The locus-model arms read at two levels of the family hierarchy: MCL clusters (M) and components of the pre-MCL graph (C) (2026-10-07)`. (full: archive/2026-11/LEDGER_O1.md:1482)

---

## GORILLA_FEWCOPY_DESIGN_2026-10-07.md

(full: archive/2026-11/LEDGER_O1.md:1554)

## Gorilla OR6737 testis few-copy class test: design sweep result (2026-10-07)

DESIGN ONLY. Human Compara 2-4 families projected by exact symbol: 342/864 members, 90 families whole; frozen criteria (`criteria_gorilla.md` 2b17d529) leave **|K| = 3** (STEAP1/STEAP1B the only clean case; C4A/C4B an O2 tie control; ZNF600/ZNF611 unreachable) -> UNDERPOWERED by rule 6, no verdict. Label "second-species, inputs-frozen class check on a partly spent library", not held-out. Read support is the binding constraint. (full: archive/2026-11/LEDGER_O1.md:1558)

---

## GORILLA_OVERLAP_2026-10-07.md

(full: archive/2026-11/LEDGER_O1.md:1597)

## Overlapping genes on real reads, gorilla OR6737 testis: our assembler against StringTie, FLAIR and isoseq collapse (2026-10-07)

Descriptive, real reads, RefSeq truth (99.9% Gnomon, circular). E_both (54 pairs, 202 chains): D 170 (84.2%), I 159 (78.7%), S 140, F 140; D - S +14.9 (CI 9.6-20.1), D - I +5.5 (CI -1.4 to +12.7, not established). N: D 94.6%, F 86.2%, I 81.2%, S 78.8%. G1, G4 held; G2 uninformative; G3 could not pass; G5 point estimate only. P/PC cost 0-1 chains at E_both. Real-read precision unmeasured. (full: archive/2026-11/LEDGER_O1.md:1601)

---

## HIERARCHY_DUPLICON_EXPANSION_2026-10-07.md

(full: archive/2026-11/LEDGER_O1.md:1643)

## Hierarchy study: results (2026-10-07). Y block, concordance with the Yoo et al. 2025 gorilla expansions

DESCRIPTIVE, no verdict class. **T1-d gorilla SEDEF arm** (tau 0.90 primary): 30 of 38 depth-matched same-family pairs share a SEDEF atom class, S-rate 0.789 [0.637, 0.889]; different-family B_viol 9 of 11,422; misses are retrocopies (GLUD1/GLUD2 ...) and old tandem paralogs. **Y1 LRPAP1** (spent): Table VIII.38 10/10 matched, our 11 loci all in Yoo's catalogue (11th = chr16 fragment LOC115932756, in VIII.35 only); BASE nodes recover 8/10 KB3781, 10/10 OR6737. **Y2 chr1 MAPKBP1/JMJD7-PLA2G4B/SPTBN5**: JMJD7-PLA2G4B 8/9 KB3781, 6/9 OR6737 (containment matches); MAPKBP1, SPTBN5 0/9 at exon level; Amendment 3 prediction (>= 5/9 MAPKBP1) **failed**; cause = unexpressed copies / too few reads on one chain. **Y3 PSMA5**: ancestral + 3 duplicates held (3/3). **Y4**: not yet (needs Phase E). (full: archive/2026-11/LEDGER_O1.md:1647)

---

## LOCUS_UNITS_LEVELS_2026-10-07.md

(full: archive/2026-11/LEDGER_O1.md:1739)

## The locus-model arms read at two levels of the family hierarchy: MCL clusters (M) and components of the pre-MCL graph (C) (2026-10-07)

DEV, ideal windows, Amendment 2 (post hoc level C). Pooled IDEAL-FOUND R+E at level C: D 60, P 67, PC 67, Q 64, S_Q 62, F1Q 64, S_D 55, F1 60. Q (primary-first representative) is a clean repair (4 gene-runs, NPIPA2 x2, TBC1D3P4/P3; Compara F unchanged); P/PC fail the F clause once out-of-universe members count as false. **Only Q stays CANDIDATE_C, LEVEL-SPECIFIC (earns nothing by itself).** F1 changed nothing for any target. H1 clause 2, H3 failed; H2, H5 held; H4 unmeetable. (full: archive/2026-11/LEDGER_O1.md:1743)

---

## PREREG_gorilla_overlap_2026-10-07.md

(full: archive/2026-11/LEDGER_O1.md:1788)

## Pre-registration: overlapping genes on real reads, gorilla OR6737 testis: our assembler against StringTie, FLAIR and isoseq collapse (2026-10-07)

Descriptive check on a partly spent library (cross-individual). Truth = valid RefSeq chains expressed at >= 3 reads (P1 primary / P2 seeding pool); strata E_both, E_one, A, N (E_j / E_x). Predictions G1 D > S, F at E_both and N; G2 P loses >= 3 at E_both; G3 D's annotation-match share below S, F; G4 D resolves >= S at E_x; G5 I < D. Errata: E_both realised 54 pairs; G2, G3 mis-specified. Result in `## Overlapping genes on real reads, gorilla OR6737 testis: our assembler against StringTie, FLAIR and isoseq collapse (2026-10-07)`. (full: archive/2026-11/LEDGER_O1.md:1792)

---

## PREREG_hierarchy_duplicon_expansion_2026-10-07.md

(full: archive/2026-11/LEDGER_O1.md:1830)

## Pre-registration: do duplicon boundaries organise families (T1), and do orthology-defined expansions sit inside our families and below their roots (T2)? (2026-10-07)

DRAFT 3 (KEY=dupfamboundary T1, KEY=expansionfamily T2). T1 gene-level: D1 (family genes share duplicon content up to nesting) and D2 (converse); T1-a REFUTED iff any exception (Q1 expects REFUTED); T1-b/T1-c homology-preserving re-cut nulls (CUT-ALIGNED p <= 0.01 / PARTIAL / NOT-DISTINGUISHED), partly tautological; T1-d gorilla SEDEF descriptive; T1-e (Amendment 1) copy number / location. T2: E-table of gorilla expansions (single-copy in H, C, O) via cross-species Liftoff tier B; T2a containment; T2b younger-than-root, predicted UNDERPOWERED. Amendment 1: go for R0-R2 (DNA side). Amendment 2: Phase E (see `## Amendment 2 — DRAFT v2 for the committed prereg (coordinator, 2026-10-07 evening). NOT committed.`). Amendment 3: Yoo et al. 2025 concordance Y1-Y4. Results in this ledger: T1-d and Y1-Y3 only (`## Hierarchy study: results (2026-10-07). Y block, concordance with the Yoo et al. 2025 gorilla expansions`); no T1-a/b/c, T2 or Y4 result is recorded here. (full: archive/2026-11/LEDGER_O1.md:1834)

---

## PREREG_hierarchy_duplicon_expansion_REVIEW_DISPOSITION_2026-10-07.md

(full: archive/2026-11/LEDGER_O1.md:2190)

## Disposition of the three adversarial reviews of the T1/T2 pre-registration (2026-10-07)

Disposition tables for review wf_20046203-7d9 (92 findings), closure review of draft 2 wf_eea75a23-1c5 (two Gate 7 blockers, 25 majors -> DRAFT 3) and regression review of draft 3 wf_c1c3d967-0ee (one blocker: T2a floor of 20 control components unreachable -> DRAFT 3a). Ends with the reachability table of every verdict class (T1-a REFUTED expected; T2b UNDERPOWERED expected). (full: archive/2026-11/LEDGER_O1.md:2194)

---

## PREREG_ideal_arm_pilot_2026-10-07.md

(full: archive/2026-11/LEDGER_O1.md:2493)

## Pre-registration: the ideal arm of the expression ladder, pilot on the 2,334 Soto genes (2026-10-07)

Descriptive pilot, no verdict class: simulate every Soto gene (+-25 kb windows) at ideal depth, run the default, score one-to-one nodes by biotype. Amendment 1 (2026-10-08, DRAFT until committed): the 71.0% baseline is a name join (53.2% by gene id), new scorer `score_soto_nodes.py`, two depths 10/30 from one pool, 2,224 simulable genes, ordered loss classes, predictions P1-P6 (P1 matched share 0.88-0.97 at depth 30; node_ok < 0.60 would not support "the problem is expression"). No result of this pilot is recorded in this ledger. (full: archive/2026-11/LEDGER_O1.md:2497)

---

## PREREG_seed_pool_real_reads_2026-10-07.md

(full: archive/2026-11/LEDGER_O1.md:2581)

## Pre-registration: the seeding pool (primary / good / all) and the primary-first representative on real reads, NPIP and TBC1D3 (2026-10-07)

Descriptive; nested pools P ⊆ G(rho) ⊆ A (`--seed-pool`, `--seed-as-ratio`) x representative MR / R1 (Rule 1). Predictions S1-S8, UNDERPOWERED if E < 8. Human outcomes (Amendment 1): chr16 NPIP S1 held (M2 12 >= 8 >= 2), S2 held, S3 held (NP .72 >= .59 > .41), S4 failed (F moved .022/.017), S5 held (2 -> 9), S6 held, S8 failed; chr17 TBC1D3 S1 failed, S3 failed, S4 failed, S5 held, S2/S6/S8 held. Gorilla uses Amendment A found criterion (NPIP E = 9; TBC1D3 E = 7 UNDERPOWERED); Amendment 2 defers the `all` pool on gorilla. Errata: M2 is not strictly "within own nodes"; S2/S1/S3 at NPIP were not blind; human arm dirs lost by `rm -rf` and regenerated byte-identical. (full: archive/2026-11/LEDGER_O1.md:2585)

---

## SEED_POOL_REAL_READS_2026-10-07.md

(full: archive/2026-11/LEDGER_O1.md:2696)

## Which alignments seed the loci, and what a primary-first representative adds: NPIP and TBC1D3 on real reads (2026-10-07)

DESCRIPTIVE (human A119b DEV; gorilla OR6737 second-species). NPIP chr16 along P, G100, G995, G98 (default), G95, G90, A: M2 exact representatives 12, 10, 8, 8, 6, 3, 2; NP .725 ... .411; A costs 22.9x aligner records and splits NPIP into 15 clusters. TBC1D3: default beats P on family (CF185 9/9 vs 5/9) but CP* .59 vs .81. 0.98 is neither knee nor plateau (S8 failed at NPIP). Rule 1 repairs NPIP G98 8 -> 11, A 2 -> 9, but |dF| .022 (S4 FAILED). Gorilla NPIP: exact-intron chains under `all` are 93% other copies' reads; S1 FAILED, S2 HELD, S3 FAILED as registered. Native `--representative primary-first` not justified by one dev library. (full: archive/2026-11/LEDGER_O1.md:2700)

---

## PREREG_family_hierarchy_flag_2026-10-09.md

(full: archive/2026-11/LEDGER_O1.md:2850)

## Family hierarchy under a flag (SD -> family -> expansion -> duplicon, with fusion loci and external families), 2026-10-09

`mcl_families --emit-hierarchy` (opt-in): family stays a partition; SD block, superfamily (graph component), expansion (maximal duplicon signature), duplicon, fusion locus, external family are typed annotating layers. **H0 met** (byte-identical) -> flag ships opt-in. **H1 NOT met**: 313/400 Soto families (78.2%) in one MCL family (bar 90%), weighted 0.892. **H2**: 215/400 exactly one of our expansions. H3 met (PKD1P6-NPIPP1 fusion_of NPIP family). Amendment 1 **H1' NOT met**: 346/400 (86.5%) in one superfamily, 0.940. Amendment 3 H4 uninformative (largest SD complex 89.6%); H5 typing EXPANSION 215 / WITHIN-FAMILY 98 / WITHIN-SUPERFAMILY 33 / WITHIN-SD-COMPLEX 52 / OUTSIDE 2. Amendment 3b H4b 52/54 but guard fails (21.7% sd_linked): uninformative. Final: Soto is a typed level of our hierarchy, not a strict subset. (full: archive/2026-11/LEDGER_O1.md:2854)

---

## PREREG_paf_duplicon_units_2026-10-10.md

(full: archive/2026-11/LEDGER_O1.md:2949)

## Our own duplicon units from the locus self-alignment (PAF units), 2026-10-10

Units cut from `loci.paf` (identity >= 0.90, >= 1 kb, exonic, CIGAR-projected breakpoints, reciprocal-majority union-find) to replace DupMasker in `--emit-hierarchy`; DEV/HELD-OUT split by family-id hash. ADOPT iff P1 (mean family score >= 0.80) and P2 (beats one-block null, sign test) on held-out. DEV: P1 FAILS 0.630, P2 holds. **HELD-OUT: P1 FAILS 0.638, P2 FAILS (12 wins / 8 losses, p 0.25) -> NOT ADOPTED**; DupMasker stays. P3 holds (223 Soto families exactly one expansion vs 215), P4 holds (both named fusions). Gorilla Phase 2 not run. (full: archive/2026-11/LEDGER_O1.md:2953)

---

## PREREG_assembler_compatible_support_2026-10-11.md

(full: archive/2026-11/LEDGER_O1.md:3055)

## Compatible-read support in the assembler's chain admission, 2026-10-11

`copy_assign --compat-support`: admit a chain when every junction has >= 2 reads counting compatible sub-chain reads. Bars on held-out (gorilla minus NC_073244.2, human A119b) both pools: B1 intron sens +1.0, B2 chain precision -<= 0.5, B3 chain sens not lower. **NOT ADOPTED, kept opt-in** (register 1264): B1 fails in all four (+0.3/+0.3/+0.4/+0.4); B2, B3 hold. Supporting reads are staggered. (full: archive/2026-11/LEDGER_O1.md:3059)

---

## PREREG_direct_transfrags_2026-10-11.md

(full: archive/2026-11/LEDGER_O1.md:3131)

## Direct-transfrag assembly (no network flow), 2026-10-11

`copy_assign --direct-transfrags`: non-maximal transfrags transfer reads to their maximal containers; admit at abundance >= 2. **NOT ADOPTED as default (B2 fails), kept opt-in**: B1 holds (intron sens +3.6/+3.6/+3.6/+3.5), B2 fails (chain precision -1.7/-1.7/-1.5/-1.2), B3 holds. Human primary ON above StringTie on every intron and chain number (intron 64.4/59.0 vs 63.5/57.8). (full: archive/2026-11/LEDGER_O1.md:3135)

---

## PREREG_direct_transfrags_union_2026-10-11.md

(full: archive/2026-11/LEDGER_O1.md:3196)

## Direct transfrags + exact chains (union), 2026-10-11

`--direct-transfrags-union`: default admissions plus maximal transfrags by abundance. Fresh held-out KB3781, chimp PTR, orangutan PPY. **NOT ADOPTED as default (B2 fails on every held-out arm), kept opt-in**: B1 holds (+1.9 to +7.0 intron sens), B3 holds (+0.6 to +1.7), B2 fails (chain precision -0.8 to -1.5). Human union above StringTie on intron, chain, transcript, locus; flow-only transcripts a small part of the gap (dev 9 of 195). (full: archive/2026-11/LEDGER_O1.md:3200)

---

## PREREG_inbetween_operating_point_2026-10-11.md

(full: archive/2026-11/LEDGER_O1.md:3265)

## An operating point between StringTie and FLAIR that wins on multi-copy loci, 2026-10-11

OP = `--direct-transfrags-union --read-isoform-k 2` (from the difference model), pools primary/good, gorilla OR6737 + human A119b minus chr20. **Result: B1, B2, C1, C2, C3 PASS in both species.** Chain sens OP_primary gorilla 29.3 vs StringTie 24.4 / FLAIR 25.5, human 23.6 vs 18.1 / 22.4; chain precision > FLAIR (30.4 vs 17.3; 15.6 vs 6.2). Secondary-only stratum S1: OP_good 34 / 36 matches vs 0 for both tools. Limits: S1 sensitivity low (6.3% / 3.9%); the multi-copy advantage comes from the good-secondaries pool, not the admission rules; gain over `default_good` small. Amendment 1: crash fix, per-contig runs. (full: archive/2026-11/LEDGER_O1.md:3269)

---

## PREREG_lrpap1_seed_reach_2026-10-11.md

(full: archive/2026-11/LEDGER_O1.md:3369)

## How many seeds does the DNA mode need to find every LRPAP1 copy? (gorilla mGorGor1), 2026-10-11

Guided DNA-mode finder (`bench/guided_pipeline.py`), 11 LRPAP1 loci. **E1 and E2 hold**: any one copy as seed finds all 11 in one pass; 0 extra loci; transcript finder alone misses the 3 fragments from seeds c03/c04/c06. Amendment 1 sweep: one seed suffices up to finder identity 0.9775 (0.98 with iteration); shipped 0.80 sits 0.18 below the split point. (full: archive/2026-11/LEDGER_O1.md:3373)

---

## PREREG_tool_difference_model_2026-10-11.md

(full: archive/2026-11/LEDGER_O1.md:3411)

## Difference model: why a chain is in our output and not StringTie's or FLAIR's, 2026-10-11

Classes for chains exclusive to one tool (`diff_model.py`). Held-out gorilla NC_073242.2 + NC_073245.2, human chr18 + chr22. **Result: model HOLDS, D1-D6 PASS in both species** (register 1267). StringTie's lead = single-read chains whose junctions other reads share (single_shared 63.6% / 36.1%); FLAIR's = single_novel; our near_stronger and two_reads chains are reference-exact more often. Basis for OP (keep ours + StringTie's single_shared via `--read-isoform-k 2`). Neither tool uses secondary-only coverage. (full: archive/2026-11/LEDGER_O1.md:3415)

---

## GENOMEWIDE_FAST_2026-10-12.md

(full: archive/2026-11/LEDGER_O1.md:3502)

## Running the pipeline genome-wide fast on a small machine (2026-10-12)

Engineering note (machine: 5 cores, 25 GB). Baseline: human families 45 min / 12.6 GB (all-vs-all asm20 over 2 GB of locus spans); human assign reaped at 16-20 GB. Done: timed resumable `all` stage, `index`, `--low-mem`, `flag --index-gb`, `utilities window-depth` + `--max-window-records`, `copy_assign --tied-prefilter` (driver default; KB3781 7 tables cmp-identical, 5.2 -> 3.6 GB). `--skip-single-placement` withdrawn (MCL188 failed). PAUSED 2026-10-12 with a 6-item to-do list (human assign, flag both index paths, whole-pipeline run, families profiling). (full: archive/2026-11/LEDGER_O1.md:3506)

---

## PREREG_mcl_two_node_2026-10-13.md

(full: archive/2026-11/LEDGER_O1.md:3625)

## PREREG — MCL two-node components: the fix (2026-10-13)

Bug: lone two-node components split by float residue in `annotation_families::mcl`. **F1, F2, F3, F4 PASS** (F2: 0 pairs split, 0 other changes; joins A119b 68, testis 9, chimp 20, gorilla 15; famsim snp/shuffle rungs 0 -> 1 intact). Amendment 1: below w* = 0.456885 exact MCL splits deterministically, so those joins are a definition change; uniform triangles/K4 also split. **Decision (user): default `--pair-rule exact-mcl`** (join at w >= w*), `weight-blind` opt-in; all four catalogs byte-identical under the default. (full: archive/2026-11/LEDGER_O1.md:3629)

---

## PREREG_o1_hierarchy_eval_2026-10-13.md

(full: archive/2026-11/LEDGER_O1.md:3701)

## PREREG — O1 hierarchy-aware evaluation, human (2026-10-13)

Are Soto 2025's families nodes of our hierarchy? Hungarian one-to-one on Jaccard over family/expansion/superfamily/SD-complex nodes, A119b dev and testis held-out, global and within-family nulls. **E1 PASS (global nulls p 0.000999), E1b FAIL (within-family p_exact 0.083 / 0.490), E2 FAIL (exact rate 63/247 = 0.2551 vs 41/84 = 0.4881, diff 0.233 > 0.10).** Reading row 2: sample-specific, only per-catalog numbers; no hierarchy claim for Soto. Hierarchy exact over flat adds only expansion matches (5 and 4). Secondary (shipped / exact-mcl) gives same exact counts 63 / 41. (full: archive/2026-11/LEDGER_O1.md:3705)

---

## PREREG_sedef_rs_layers_2026-10-13.md

(full: archive/2026-11/LEDGER_O1.md:4038)

## sedef-rs as the SD and duplicon layers: V1, V1b, V2 (2026-10-13)

Vendored sedef-rs `utilities sd-scan` vs SEDEF (V1 at 0.90, V1b at 0.98; recall and precision >= 0.80 on CHM13 v2.0 and gorilla) and vs DupMasker (V2a held-out family score >= 0.80; V2b Soto expansion-not-family >= 9). **V1 FAIL, V1b FAIL** (recall 0.3980 / 0.4067 at 0.90; 0.3968 / 0.4236 at 0.98; precision 0.82-0.95). **V2 FAIL (V2a FAIL 0.5153, V2b PASS 12 vs bar 9).** Not adopted: SD layer keeps SEDEF, human expansions stay on DupMasker, apes have no expansion layer. (full: archive/2026-11/LEDGER_O1.md:4042)

---

## PREREG_denovo_copies_2026-10-14.md

(full: archive/2026-11/LEDGER_O1.md:4476)

## De novo analogue of Liftoff's `-copies` (2026-10-14)

Project each base locus's representative transcript to the genome at Liftoff thresholds (cov >= 0.5, id >= 0.95); C1 recall (own) of Liftoff extra copies >= 0.80 (>= 30 reachable), C2 unexplained sites <= 0.10. **Reading 3: not claimed.** C1: A119b 268/420 = 0.638 FAIL, testis 49/62 = 0.790 FAIL, KB3781 44/52 = 0.846 pass; C2 FAIL in all three (0.408 / 0.182 / 0.359). Retrocopy advantage did not show (266 of 4,941 `none` sites retro-like). (full: archive/2026-11/LEDGER_O1.md:4480)

---

## PREREG_liftoff_pairs_null_2026-10-14.md

(full: archive/2026-11/LEDGER_O1.md:4719)

## Liftoff copy pairs against matched decoys: O1's independent-truth test (2026-10-14)

Do our default families put a gene's Liftoff copies together more often than distance- and coverage-matched decoys? Bar p < 0.001 at identity >= 0.95 in A119b, testis, KB3781 (held out). **Reading 1, PASS in all three** (p 0.000999): shared pairs A119b 75/514 (0.146) vs null 0.021; testis 12/38 (0.316) vs 0.032; KB3781 4/78 (0.051) vs 0.0006. Not a precision measure. (full: archive/2026-11/LEDGER_O1.md:4723)

---

## PREREG_soto_subset_cn_2026-10-14.md

(full: archive/2026-11/LEDGER_O1.md:4938)

## Soto's families as our families cut by copy number: T1 and T2 (2026-10-14)

T1 Soto ⊂ ours (T1a containment vs size null p < 0.001; T1b D = H(ours|Soto) - H(Soto|ours) interval above 0); T2 ours cut by copy number (|Δ| < 2) rebuilds Soto (exact and ARI beat null). **Reading 4: T1 fails, no subset claim.** Guided: T1a PASS (258/491 contained), T1b FAIL (D -0.3995), T2 PASS (exact 186, ARI 0.8274) reported not claimed. De novo A119b and testis: T1a PASS, T1b FAIL, T2 FAIL. (full: archive/2026-11/LEDGER_O1.md:4942)

---

## PREREG_template_copies_2026-10-14.md

(full: archive/2026-11/LEDGER_O1.md:5260)

## Template-guided copy assembly (2026-10-14)

Assemble copies with own reads but no locus, guided by an expressed same-gene template. **Reading 3: not adopted** (register 1285). A119b_shipped dev fails every precision bar: P1 chain precision 17.7 -> 17.4 (bar >= 17.5), located share 0.2452 -> 0.2091, P2 0.0813 vs 0.2318. Recall rises 234 -> 394 of 641 read-supported Liftoff copies (119 of 160 by single-exon added loci); 26,726 added loci. Final review found three defects in the build; a fixed build untested. (full: archive/2026-11/LEDGER_O1.md:5264)

---

## PREREG_post_mcl_reassign_2026-10-15.md

(full: archive/2026-11/LEDGER_O1.md:5480)

## Post-MCL strongest-edge reassignment of family members (2026-10-15)

After MCL, move a member to the family holding its strongest edge (motivated by A119b figure-8 F1 131/148 < 0.90). Bars D (A119b F1 >= 0.90), H1/H2 held-out (testis, KB3781), guards. **Reading 3: not adopted**: guards fail in two of three samples; fixes A119b F1 but 2,006 moves (only 9 in the chr19 array), Soto pairwise precision -0.046; no held-out benefit. (full: archive/2026-11/LEDGER_O1.md:5484)

---

## NPIP_16Q22_CNV_2026-10-20.md

(full: archive/2026-11/LEDGER_O1.md:5626)

## The 16q22 NPIPB15 repeat: why LOC124907808 and LOC124907807 rest on tied reads (2026-10-20)

CHM13 chr16:80.10-80.55 Mb has three tandem units (NPIPB15, LOC124907808, LOC124907807; 114.4 kb apart); every haplotype checked (HG002 mat/pat, HG00733 mat/pat) carries 1 or 2 units, so the two LOC copies are absent from some haplotypes; A119b primaries at them are mostly exact or near ties (LOC124907808 54 exact ties of 93). Gorilla tab: KB3781 STON1 present on maternal chr12 only (121,624 bp deletion on paternal; 113/113 O3 reads align to mother); both parents carry 5 NPIPB15-like copies. Figure in the NPIP Read Pools page. (full: archive/2026-11/LEDGER_O1.md:5630)

---

## MULTIMAPPER_PREMISE_2026-10-21.md

(full: archive/2026-11/LEDGER_O1.md:5713)

## Multi-mapping long reads: what they cost assembly, what O2 can decide, what O3 recovers (2026-10-21)

Summary of registers 1360-1367. Exact AS ties 0.14% (OR6737) / 0.52% (A119b) of reads; StringTie/FLAIR are primary-only in practice; ~3% of expressed transcripts are multi-placed, where tools recover about half (shortfall ~250-340 per tool per species); seeding secondaries adds ~0.2 points genome-wide and gains are 76% / 83% O2-undecidable; multi-mapping loci are real SDs (99.1% human SD90, 98.2% gorilla SEDEF); O2 in matched genome UNDERPOWERED on real reads, FAILS on all-copies sim; O3 recovers the maternal-only NPIP copy (33/39 reads after adding candidates, 0 false moves). "Multi-mapping is a genome-wide loss" not supported; "tools lose about half in multi-copy loci" supported. (full: archive/2026-11/LEDGER_O1.md:5717)

---

## O1_FAST_2026-10-21.md

(full: archive/2026-11/LEDGER_O1.md:5797)

## O1 fast: the exon gate and the families stage's Rust side (2026-10-21/22)

Patched minimap2 (`tools/minimap2_gate`, `2.30-r1287-gate3`..`gate6`, `--gate FILE`, `--dp-memo`) aligns only locus pairs that can carry exon-to-exon evidence; on by default when `--min-exonic-bp >= 1`. Byte-identical products on every catalog. gorilla_OR6737 all-vs-all 5,421 s -> 182 s; human_testis families stage 71 s; human_A119b 7,332 s (gate3) -> 2,270 s (gate6), about 2.6x vs ungated ~5-5.5 h; Rust side 109 s -> 9.1 s. A first gate version was wrong (kept as a negative result). Remaining time: DP kernel and the chr13 A119b hotspot (64%). Build: `tools/minimap2_gate/build.sh`. (full: archive/2026-11/LEDGER_O1.md:5801)

---

## PREREG_multicopy_gains_real_2026-10-21.md

(full: archive/2026-11/LEDGER_O1.md:6166)

## Pre-registration: are the multi-copy gains real copies, and are those loci real duplications? (2026-10-21)

Gains = high-tie transcripts matched only with secondary seeding (gorilla 162, human 284). **Q1: gains mostly NOT the copy's own expression** (supported 12% / 15%; human 70% contradicted = shadows; gorilla 45% contradicted / 43% undecidable). **Q2**: high-tie genes in SD90: human 1,919/1,936 (99.1%) holds; gorilla 90.4% on our scan. **Audit**: Q1 holds with caveats; Q2 human holds with caveats (77.3% without satellites; close to true by construction); **Q2 gorilla does not hold as stated** (82.9% by est_id FAIL; SEDEF post hoc 98.2%); "supported" mainstays are certified 97-100% to the other copy; reading withdrawn. Register 1367. (full: archive/2026-11/LEDGER_O1.md:6170)

---

## PREREG_o1_families_sd_pairs_2026-10-21.md

(full: archive/2026-11/LEDGER_O1.md:6350)

## Pre-registration: are O1's families real duplications? Copy pairs against independent segmental-duplication annotations (2026-10-21)

Bar: >= 90% of O1 copies at max_family_identity >= 0.90 are SD-paired with another copy of their family (SEDEF / CHM13 SD track). **Bar FAILS in both catalogs: 72.1% human, 82.4% gorilla.** Audit: FAIL holds with caveats (repeats on both sides of the bar); post hoc, gene-like families are paired at about the bar (93.8% human, 92.0% gorilla). Settling it needs a CAT-annotation families arm on CHM13. (full: archive/2026-11/LEDGER_O1.md:6354)

---

## PREREG_o1_self_hits_2026-10-21.md

(full: archive/2026-11/LEDGER_O1.md:6493)

## Pre-registration: dropping diagonal self-hits from O1's family graph (2026-10-21)

Defect (register 1368 audit): overlapping de novo loci align to themselves (A119b: 253,893 diagonal records between different loci). `mcl_families --reject-self-hits` drops same-contig + strand coincident-interval records. v1 **FAIL** (gorilla B2 97.4%); Amendment 1 v2 rule (strands from loci GFF) PASS on the 3 dev catalogs (development); Amendment 3 fresh held-out **gorilla_OR6737 PASS**; human_A119b (dev) fails G1 (1,036 -> 1,140 > 1,088). **Decision (user, 2026-10-21): `--reject-self-hits` is the default** (`--no-reject-self-hits` restores old graphs). Copies: A119b 9,907 -> 9,355, testis 1,115 -> 1,069, KB3781 918 -> 843, OR6737 1,607 -> 1,521. (full: archive/2026-11/LEDGER_O1.md:6497)

## SDSCAN2 — sd-scan v2: find SEDEF's segmental duplications without SEDEF (prereg 2026-10-25)

After r1278 V1 FAIL (recall 0.398 / 0.407): scan v2 changes R1 (keep every collinear run), R2 (fracMatch identity), R3 (half-overlap seeding), then SEDEF/BISER ports. Held-out bar H1 recall and precision >= 0.80 on CHM13 chr1/9/15/16; H2 v2 recall > v1. Dev only (gorilla subset): v1 0.4898 / 0.7659; v2 repeat-blind 0.8478 / 0.6972; SEDEF rules block gate 0.7079 / 0.9653; scan v3 (2026-10-26) 0.7178 / 0.9530. None reaches 0.80 / 0.80 on dev; **held-out not run**. Remaining loss mostly long, mostly soft-masked SDs SEDEF aligns through. (full: archive/2026-11/LEDGER_O1.md:6849)

## PREREG SEMISEC (2026-10-26): semi-guided mode from multimappers - silent copies of expressed genes, classified

Seeds = read-supported loci; candidates = merged secondary hits (AS >= 0.8 x primary, >= 2 reads) with < 2 primaries (silent); classes PROCESSED/UNPROCESSED x DISRUPTED/INTACT/LNC-LIKE from sequence. Dev OR6737 census, then bars. **Held-out KB3781 (register 1431): PASS B1-B3**: B1 pseudogene DISRUPTED/INTACT 3.52 (>= 1.5), B2 PROCESSED/DISRUPTED pseudogene 0.474 (>= 0.40), B3 recall 151/664 = 0.227 (>= 0.10); **B4 lncRNA FAIL** (1.09 < 1.25). Rust port `src/semisec.rs` byte-identical reports. (full: archive/2026-11/LEDGER_O1.md:6947)

## PREREG SEMIRATIO (2026-10-26): drop the 0.8 AS cut?

`semisec --ratio 0` vs 0.8, fibroblast decisive. **FAIL (register 1432), 0.8 stays default**: C1 PASS (recall 0.227 -> 0.548), C2 FAIL (B2 0.359 < 0.40), C3 FAIL (parent agreement 0.508 < 0.708 - 0.05). (full: archive/2026-11/LEDGER_O1.md:7026)

## PREREG SEMITRUTH (2026-10-27): are SEMISEC's unannotated silent candidates real gene copies?

Project candidates to CHM13 (asm20, MAPQ >= 1) and call PSEUDO by CAT/Liftoff. **PASS in both tissues (register 1433)**: T1 P PSEUDO 0.829 / 0.792 (>= 0.50); T2 U 0.120 / 0.106 vs N 0.018 / 0.024 (>= 2x). ~180-190 real unannotated copies per sample (~13%); U mostly repeat-driven (70-77% RepeatMasker). (full: archive/2026-11/LEDGER_O1.md:7055)

## PREREG SEMIREP (2026-10-27): drop repeat-driven secondary hits

`semisec --repeats`: drop hits mostly in RepeatMasker. **FAIL (register 1434), opt-in**: R1 PASS (U PSEUDO rate 0.106 -> 0.263), R2 FAIL (count 152 -> 91 < 122), R3 FAIL (recall 0.202 < 0.207). Opt-in high-precision list. (full: archive/2026-11/LEDGER_O1.md:7095)

## PREREG DENOVOFAM (2026-10-27): better de novo multi-copy families against Compara, mainstays and YAGs

Compara-truth scorer `bench/denovo_eval/famtruth.py` (recent levels, dev/held-out halves by ENSG hash; mainstays/YAGs by regex). Arms A inflation, B `--shared-exon-union`, D/D' `--fold-split-bridges` (Amendments 1-3: protein-coding-first mapping; bridge rule refined to homology-joined pairs). Dev: D' selected (pair F 0.764 -> 0.826, family F 0.793 -> 0.810). **FAIL (register 1435): H1, H2 fail** (A119b held-out 0.956 / 0.903 unchanged; testis byte-identical); D' INERT outside dev families; stays opt-in. (full: archive/2026-11/LEDGER_O1.md:7130)

## PREREG DENOVOFAM-G (2026-10-27): `--fold-split-bridges` (D') on FRESH gorilla catalogs

**FAIL (register 1436): no gain, truth underpowered** (only 9 / 13 true families with >= 2 expressed members, LOC blind spot); G1-G3 pass, +0.01 family bar fails on both; stays opt-in. (full: archive/2026-11/LEDGER_O1.md:7193)

## PREREG DENOVOFAM-L (2026-10-27): Liftoff-style exon chaining as family edges (pending item 1, user go 2026-10-27)

Lift edges (`liftchain.py`, coverage >= 0.5, identity >= 0.5) via new opt-in `mcl_families --extra-edges [--extra-edges-only]`; arms E1/E2, Amendment 1 adds E3 (0.70 floor). **FAIL at dev selection (register 1437)**: every arm +0.028 family F but precision drop 0.046-0.062 > 0.02; GOLGA6L foreign 2 -> 15. One-sided coverage is a placement criterion, not a family one. (full: archive/2026-11/LEDGER_O1.md:7212)

## PREREG DENOVOFAM-L2 (2026-10-27): two-sided lift edges

E4 reciprocal lifts. **FAIL at dev selection (register 1438)**: pair F 0.747 (baseline 0.764), family F 0.788 (0.793), DAZ splits. (full: archive/2026-11/LEDGER_O1.md:7251)

## PREREG DENOVOFAM-L3 (2026-10-27): lift edges only for orphan loci

E5 = lift edges only where a locus has no baseline edge. **FAIL at dev selection, recent truth (register 1439)**: pair F 0.758 vs 0.764, family F 0.822 vs 0.793; 13 "false" pairs mostly real older paralogues; against Simiiformes/Primates E1 gains on every number. (full: archive/2026-11/LEDGER_O1.md:7264)

## PREREG DENOVOFAM-L4 (2026-10-27): orphan lift edges (E5) held out against Simiiformes truth (user decision 2026-10-27)

User decision: primary truth = Compara Simiiformes. E5 qualifies at dev; held-out first read PASS (A119b family F 0.840 -> 0.858; testis 0.865 -> 0.897). **Audit: PASS RETRACTED -> FAIL (register 1440)**: famtruth.py trimmed cross-half pairs; counted, A119b held-out pair F 0.8903 -> 0.8885 (lower); pass only at Simiiformes chosen after recent fail; 91.5% of new pairs outside the judgeable universe; mainstay content changed (DAZL joins DAZ). famtruth.py fixed same day; every earlier DENOVOFAM number reads high on precision. `--extra-edges` opt-in. (full: archive/2026-11/LEDGER_O1.md:7279)

## PREREG DENOVOFAM-L5 (2026-10-27): orphan lifts that do not widen a family

E6 orphan-orphan lifts; E7 = E6 + orphan-to-family lifts not widening the family's weakest edge. E7 qualifies at dev; held-out A119b pair F 0.8903 -> 0.8904, family F +0.012; testis family F +0.031. **Audit: NOT VALIDATED (register 1441), opt-in only**: held-out sets not fresh for E7; A119b margin +0.0001 is an MCL side effect (KZNF split); 94% / 90% of new pairs unjudgeable. Needs a fresh substrate (e.g. HG002). (full: archive/2026-11/LEDGER_O1.md:7311)

## PREREG GUIDED-S (2026-10-27): guided mode in SPLICED mode - CDS seeds find their family's other members and pseudogenes

Guided seeds: A gene body asm20, B CDS `-x splice`, C miniprot, D (Amendment 1) = B union CDS asm20; truth Compara Simiiformes + name-proxy pseudogenes. **FAIL as registered (register 1442), precision**: dev A 0.423 / 0.035 / 0.882; D 0.806 / 0.259 / 0.633 (member / pseudo recall / precision). CDS seeds find spliced members and processed pseudogenes (pseudo recall x7) but also older multigene relatives. B's held-out seen by mistake (deviation). (full: archive/2026-11/LEDGER_O1.md:7341)

## PREREG GUIDED-S2 (2026-10-27): spliced seeds judged against MULTIGENE-family truth (user: HGNC gene groups)

Truth = HGNC most-specific gene group OR Compara Simiiformes pair; macro recall. First read PASS on held-out (D 0.038 / 0.126, 0.956 vs A 0.016 / 0.004, 0.859). **Audit: NOT VALIDATED (register 1443)**: mainstay guard empty (foreign seeds); with real NPIPB12 seed A 21/22, D 6/22 (asm20 per-record coverage floor); named gains were dev families; unjudged hits 26%; truth changed after r1442. Fix: chain records per target locus (Liftoff step). (full: archive/2026-11/LEDGER_O1.md:7382)

## PREREG GUIDED-S3 (2026-10-27): Liftoff-chained spliced seeds, human dev, FRESH gorilla held-out

Chained (Liftoff DAG) seeds, truth = most-specific HGNC group; arms A body vs E CDS splice+asm20; both precisions barred. **FAIL on dev (register 1444); gorilla not run**: recall and judged precision pass, strict precision fails both seeds (0.689 vs 0.713; 0.678 vs 0.807), named bar fails (NPIP 6/21 vs 13/21: E misses divergent NPIPA). Post-S3 arm U = body UNION CDS; its unjudged hits are own-family pseudogenes the proxy missed. (full: archive/2026-11/LEDGER_O1.md:7416)

## PREREG GUIDED-S4 (2026-10-27): independent pseudogene-parent judge; body UNION CDS seed; gorilla held-out

`pgparent.py` (GENCODE v50 2-way parents; gorilla RefSeq description). **FAIL on dev by 0.0016 (register 1445)**: seed 2 strict precision 0.7987 vs bar 0.8003; all else passes; pseudogene recall x10-14 (0.009 -> 0.130); gorilla kept fresh. (full: archive/2026-11/LEDGER_O1.md:7461)

## PREREG GUIDED-S5 (2026-10-27): S4 with an HGNC cross-checked protein-coding set; gorilla held-out

`--hgnc-coding`. Dev PASS both seeds. **Gorilla held-out FAIL (register 1446)**: seed-2 strict precision 0.170 vs bar 0.173; 99.8% of unjudged hits are LOC genes (truth blind spot). Held-out spent. (full: archive/2026-11/LEDGER_O1.md:7494)

## PREREG GUIDED-S6 (2026-10-27): S5 arms on a FRESH gorilla held-out, Kamilah_GGO_v0 + NCBI Annotation Release 102

Kamilah_GGO_v0 + RefSeq 102, LOC protein-coding judged by description parent. **FAIL (register 1447)**: recall and named bars pass; strict precision fails both seeds (0.223 vs 0.236; 0.203 vs 0.208), judged fails seed 2 (0.751 vs 0.758). U's extra hits land in sibling HGNC sub-groups. Kamilah spent. (full: archive/2026-11/LEDGER_O1.md:7523)

## PREREG GUIDED-S7 (2026-10-27): HGNC PARENT-group truth; fresh human held-out HG002 pat (Ensembl rapid release 2022_07)

Truth = parent HGNC group. **FAIL on dev (register 1448)**: seed-2 strict precision 0.7685 vs 0.7805; seed 1 all pass; HG002 not scored. Failures are own-family pseudogenes no judge reaches (FAM90A*P). (full: archive/2026-11/LEDGER_O1.md:7552)

## PREREG GUIDED-S8 (2026-10-27): S7 + a name-stem pseudogene judge; dev, then HG002 pat held-out

Name-stem pseudogene judge added. **S8 dev: every bar PASSES** (seed 1 U 0.0549 / 0.1036 / 0.908 / 0.841; seed 2 0.0561 / 0.1031 / 0.877 / 0.821). HG002 held-out outcome: see register 1458-1460 and `## PREREG GUIDED-S10 (2026-10-27): the de novo lessons carried to guided mode — retro class + reciprocity; HG002 pat held-out` (no S8 HG002 score recorded in this section). (full: archive/2026-11/LEDGER_O1.md:7575)

## PREREG GUIDED-S9 (2026-10-27): many templates per family (every isoform CDS + body, 2 members) vs one template

Arm M (2 members, body + every isoform CDS) vs U2 (body + longest CDS). **FAIL on dev (register 1449)**: strict precision 0.787 vs bar 0.795; member 0.0563 -> 0.0649, pseudo 0.192 -> 0.228; HG002 not scored. Two members vs one roughly doubles pseudogene recall. (full: archive/2026-11/LEDGER_O1.md:7586)

## PREREG ASMTPL (2026-10-27): assembled transcripts as family templates, gorilla, 13 families

13 keyword families (dev NPIP, TBC1D3, GOLGA6, GOLGA8, USP17L, NBPF; held-out TSPY, XAGE, CT47, RBMY, SPATA31, CDY, DAZ); R1 one template vs RK every unique intron chain. **FAIL as registered (register 1450), precision**: members dev 3/122 -> 62/122, held-out 24/56 -> 34/56; precision dev 1.000 -> 0.617, held-out 1.000 -> 0.940. Wrong hits: GOLGA6/8 cross-hits and SD-passenger neighbours. 764 lookups vs 12,021 reads. (full: archive/2026-11/LEDGER_O1.md:7608)

## PREREG ASMTPL-2 (2026-10-27): family-shared exons, Liftoff per-exon scoring

Templates trimmed to family-shared exons. **FAIL (register 1451)**: dev precision 0.617 -> 0.596, members 58/122; CT47 and TSPY lose their only template. Most neighbour hits land in host-gene introns (unannotated nested copies). (full: archive/2026-11/LEDGER_O1.md:7646)

## PREREG ASMTPL-3 (2026-10-27): exon-side judge, GOLGA6+GOLGA8 merged, full templates

Judge on target exon bases; intronic landings binned "putative"; GOLGA6+8 merged. **FAIL as registered (register 1452), but the lever works**: precision dev 0.617 -> 0.913 (63/69), held-out 0.946 (35/37); members 55/122, 31/56 vs R1 3/122, 24/56. Bar R1 1.000 - 0.02 degenerate (R1 has 3 dev hits); retrospective bar in ASMTPL-4. LIFTMODE prototype (`bench/denovo_eval/liftmode.py`) wired from this rule: 26 clean families on the 13-family loci. (full: archive/2026-11/LEDGER_O1.md:7664)

## PREREG LIFTMODE-GW (2026-10-27): Liftoff-style de novo families genome-wide vs the MCL catalog, gorilla

LIFT (unique-chain templates, exon-side placements, connected components) vs shipped MCL, KB3781 dev. **FAIL on dev (register 1453), OR6737 not run**: family F 0.455 -> 0.625 but pair F 0.182 -> 0.137, pair precision 0.267 -> 0.074 (379 predicted pairs for 29 true): components over placements are transitive. (full: archive/2026-11/LEDGER_O1.md:7697)

## PREREG LIFTMODE-GW2 (2026-10-27): reciprocal placement as the family edge; then LIFT edges into MCL

Arms R (reciprocal components), X / XR (LIFT edges into MCL via `--extra-edges`). **FAIL on the precision bar only (register 1454), OR6737 not run**: pair F 0.195 -> 0.307 / 0.299 / 0.325, family F 0.455 -> 0.704 / 0.784 / 0.727; precision drops 0.12-0.15 vs bar 0.02 (base has 12 predicted pairs). (full: archive/2026-11/LEDGER_O1.md:7720)

## PREREG LIFTMODE-GW3 (2026-10-27): retro-placement class + asymmetric pair judge; then OR6737

Retro placements removed; asymmetric judge (WRONG = no Compara level, no HGNC group, not PAR). XR3 dev and OR6737 held-out bars pass (OR6737 pair F 0.300 -> 0.394, family F 0.514 -> 0.720, WRONG 1 -> 2). **Audit: NOT VALIDATED (register 1455)**: 94.8% of added OR6737 pairs unjudgeable; bar changed twice; registered pair precision still falls (0.450 -> 0.283); requested stem-agreement check OR6737 0.760 vs base 0.846. (full: archive/2026-11/LEDGER_O1.md:7742)

## PREREG LIFTMODE-GW3-RBH (2026-10-27): protein-level judge of the unjudged pairs

mmseqs all-vs-all of 1,378 longest-CDS proteins; rule: OR6737 added-pair mutual-hit rate >= base's. **FAIL by 0.041 (register 1456)**: OR6737 0.921 vs 0.962 (KB3781 dev 0.943 vs 0.883); r1455 stays NOT VALIDATED. 92-94% of added LOC pairs are protein homologues; 40/82 unsupported pairs come from unflagged multi-gene loci (only 1/82 a flagged F1v2 bridge) -> next lever = template boundaries. (full: archive/2026-11/LEDGER_O1.md:7779)

## PREREG LIFTMODE-GW4 (2026-10-27): template boundaries — representative-only (U1) vs locus-shared exons (U2) vs all chains (XR3)

U1 = one template per locus (most-reads chain), U2 = locus-shared exons. OR6737: U1 pair F 0.300 -> 0.440, family F 0.514 -> 0.710, WRONG 1, added-pair mutual-hit 0.968 (>= base 0.962); U2 fails the protein rule (0.922). **Audit: NOT VALIDATED (register 1457)**: OR6737 not held out for GW4 (arms designed from its failing pairs); U1's edges are a strict subset of XR3's; 81% of added pairs judged only post hoc; GW2 pair-precision bar still fails (0.450 -> 0.367). Needs a fresh gorilla sample (later U1-FRESH). (full: archive/2026-11/LEDGER_O1.md:7807)

## PREREG GUIDED-S10 (2026-10-27): the de novo lessons carried to guided mode — retro class + reciprocity; HG002 pat held-out

UR = retro class, URR = UR + reciprocity (back-lift onto seed). **URR FAILS dev on member recall (register 1458)**: 0.0549 -> 0.0446 / 0.0561 -> 0.0493 (Liftoff coverage is template-relative); UR small clean gain (+0.007 / +0.014 precision) but not the registered arm; HG002 prepared, not scored. (full: archive/2026-11/LEDGER_O1.md:7843)

## PREREG GUIDED-S11 (2026-10-27): length-symmetric reciprocity (URS); then Amendment: UR/URS vs U on HG002

URS (coverage of the shorter, targeted back-lift). **URS FAILS dev by 0.0015 / 0.0006 member recall (register 1459)** -> UR goes to HG002. **Amendment: UR vs U on HG002 pat FAIL (register 1460)**: seed 2 passes all; seed 1 NPIP 24 -> 17 (7 intron-variant copies binned retro), member 0.0440 vs 0.0442. Guided default unchanged; HG002 spent for the retro question. (full: archive/2026-11/LEDGER_O1.md:7874)

## PREREG SEMIFAM (2026-10-28): do SEMISEC's silent copies belong in the de novo families? (PENDING_2026-10-27 item 2)

SEMISEC silent candidates as extra family nodes. Ruling: `rdev` = OR6737, `rho` = KB3781 (labels swapped). **FAIL on dev KB3781 (register 1461), OR6737 not run**: S/S0 collapse (candidates are all-vs-all hubs; wrong 234 / 93); S2 (join seed's family) wrong pairs 2 -> 12 (bar <= 4), parent agreement 538 = 30% (bar 50%). (full: archive/2026-11/LEDGER_O1.md:7905)

## PREREG SEMIFAM-2 (2026-10-28): SEMISEC candidates admitted by the LIFT rule, not by the read's secondary hit

S3 = lift-confirmed candidates only (1,336 / 2,857). **FAIL on the parent bar only (register 1462)**: wrong pairs 2 -> 3, pair F 0.195 -> 0.208, family F 0.455 -> 0.471 pass; parent agreement 0.389 < 0.50 (nomenclature-dominated: ribosomal renaming, "uncharacterized LOC"). OR6737 not run. (full: archive/2026-11/LEDGER_O1.md:7938)

## PREREG STRAINY-A (2026-10-28): Strainy-style clustering on the de novo homology graph (PENDING_2026-10-27 item 3a)

LPA / Louvain / LPA-recursive / CC vs MCL on the dumped graph. First read: LPA passes dev and OR6737 (OR pair F 0.273 -> 0.364, wrong 8 -> 1). **Audit: NOT VALIDATED (register 1463)**; mechanism and "wrong 8 -> 1" headline RETRACTED: "singletons" are overlapping records folded together; MCL baseline rebuilt after prereg; whole gain is one KRAB-ZNF block; seed-dependent (KB seeds 2-5 fail). Closes as r922: TIED at the DNA level; do not re-propose. (full: archive/2026-11/LEDGER_O1.md:7956)

## PREREG U1-TSS (2026-10-28): representative chains that start or end inside a sibling's intron

Intronic-terminal chains are 26% of chains / 17% of reads on KB3781. **No effect (register 1464)**: U1f pair F 0.500 -> 0.462, family F 0.764 -> 0.717; U1p (polished re-assembly) pair F 0.500, family F 0.737. Does not skew the template result. (full: archive/2026-11/LEDGER_O1.md:7999)

## PREREG SEMIFAM-3 (2026-10-28): protein-level parent judge for the lift-confirmed SEMISEC copies

mmseqs translated search of family proteins vs candidate span, permuted-family control. **Mixed -> FAIL as registered (register 1465)**: judge passes (agreement 0.865 dev / 0.861 OR6737 vs control 0.006 / 0.003); famtruth bars fail on OR6737 by a denominator effect (pair F 0.300 -> 0.286, family F 0.514 -> 0.494, same hits, universe 40 -> 43 pairs). (full: archive/2026-11/LEDGER_O1.md:8028)

## ASMTPL-4 (2026-10-28): a non-degenerate precision bar for ASMTPL-3 — RETROSPECTIVE (numbers already seen in r1452)

Bar: 95% Wilson intervals of RK and R1 precision overlap. **Retrospective PASS (register 1466)**: dev R1 3/3 [0.438, 1.000] vs RK 63/69 = 0.913 [0.823, 0.960]; held-out 28/28 [0.879, 1.000] vs 35/37 = 0.946 [0.823, 0.985]; recall bars pass. Decides no default. (full: archive/2026-11/LEDGER_O1.md:8048)

## PREREG U1-FRESH (2026-10-28): the U1 lift-edge rule on a genome no GW/ASMTPL/DENOVOFAM row has touched — Kamilah_GGO_v0 + RefSeq 102

Full GW4 bar set on Kamilah with OR6737 reads. **WITHDRAWN (user, 2026-10-28 13:08)** before any alignment finished: Kamilah is not T2T, so the young copies are missing. No data looked at. (full: archive/2026-11/LEDGER_O1.md:8060)

## PREREG U1-FRESH-HG002 (2026-10-28): the U1 lift-edge rule on HG002 T2T Q100 paternal + Ensembl rapid release 2022_07

Replaces U1-FRESH: HG002 pat + Ensembl 2022_07 + HG002 Kinnex FLNC, same pipeline and full GW4 bars. Run as the intron-hit variant: see `## PREREG U1-FRESH-HG002-INTRON (2026-10-29, user "Sounds good, please proceed"): U1 lift edges + the intron-homology filter, HG002 T2T Q100 paternal held-out` (register 1469). (full: archive/2026-11/LEDGER_O1.md:8077)

## PREREG EDGE-HUMAN (2026-10-28, /goal "improve edge recall based on data we already have")

U1 lift edges on human A119b (dev) and testis ERR13885926 (held-out), Compara Simiiformes. Pre-MCL edge recall: A119b 0.807 -> 0.879, testis 0.824 -> 0.836. **Audit: NOT VALIDATED, PASS withdrawn (register 1467)**: edge-recall bar written after the result; registered MCL bars on testis fail (pair precision 0.956 -> 0.935, pair F 0.900 -> 0.895); dev gain mostly transitive closure; 91% / 80% of new pairs unjudged. Family F +0.03 on testis stands descriptively. (full: archive/2026-11/LEDGER_O1.md:8086)

## PREREG EDGE-HUMAN-2 (2026-10-28, user "yes" 16:05): family-calibrated identity floor for lift edges

Keep a U1 lift edge only if its identity >= the lowest within-family identity of its endpoints' shipped families. Retro rule fixed (NUP50 retrocopy in CDH18). Handed off to the cluster 2026-10-28 16:37 with a 5-step to-do. Dev result at KB3781 recorded in the register: **FAIL at dev on pair precision (register 1468)** (U1 0.297, floor 0.282 vs bar 0.321). (full: archive/2026-11/LEDGER_O1.md:8120)

## EXONSUPPORT (2026-10-28, user: "make the exons count better, based on the CIGAR"; improve the template pool before minimap2)

Descriptive (`bench/denovo_eval/exonsupport.py`, human slice chr20-21): 50.9% of transcripts carry a flagged exon; the "skip > cover" half flags legitimate cassette exons (withdrawn); 2,028 "cover < 2 primaries" flags suspected secondary-only. Secondary-inclusive rerun stopped at hand-off; no pool rule decided. (full: archive/2026-11/LEDGER_O1.md:8143)

## Column-exact shared exons: `mcl_families --shared-exon-columns` (2026-10-29, development on gorilla KB3781 fibroblast; opt-in)

New opt-in rule: shared-exon evidence = CIGAR columns joining exon base to exon base (default off, byte-identical). KB3781: 470 pairs pass the interval rule only, 0 the column rule only (352 overlapping loci of one transcription unit); recent pair F 0.310 -> 0.319, Simiiformes 0.257 -> 0.262, family F unchanged. Development, stays opt-in. Subsections in this section:
- **EDGE-HUMAN-2 RESULT, dev KB3781 (register 1468): FAIL** on pair precision (U1 0.297, floor 0.282 vs bar 0.321); the floor drops 12 of 24 true edges; U1 adds 0 wrong pairs, cost is Compara-unjudged/older pairs. Human testis held-out not spent.
- **Per-exon shared-exon rule** (descriptive): exon_frac >= 0.7 rescues 8 true pairs for 1 wrong but rides on 207 unjudged LOC pairs; not proposed.
- **Lift edges, development after EDGE-HUMAN-2** (post hoc, dev only): U1 + intron homology >= 500 gated-PAF columns is the only variant clearing every GW4 bar on dev (pair precision 0.358 >= 0.321, pair F 0.417, family F 0.585); registered as U1-FRESH-HG002-INTRON.
(full: archive/2026-11/LEDGER_O1.md:8156; EDGE-HUMAN-2 result :8212; per-exon :8247; lift edges dev :8269)

## PREREG U1-FRESH-HG002-INTRON (2026-10-29, user "Sounds good, please proceed"): U1 lift edges + the intron-homology filter, HG002 T2T Q100 paternal held-out

U1 lift edges kept iff >= 1 intron-intron aligned column in the base's gated all-vs-all (post hoc on KB3781, disclosed); HG002 Kinnex re-aligned to HG002 pat, Ensembl 2022_07 truth; GW4 bars. **RESULT (register 1469): FAIL on the held-out**: pair F 0.411 < base 0.413, precision 0.415 < 0.478; family F, WRONG (95 vs 94), mutual hit (0.993 vs 0.721) and judgeable bars pass. Gain is only on famtruth's dev split. Both rules stay opt-in. (full: archive/2026-11/LEDGER_O1.md:8318; result :8346)

## PREREG DISPERSION (2026-10-29, user "yes lets do that"): chromosome-dispersion filter on de novo families, held-out on human A119b (CHM13 v2.0) and gorilla OR6737 testis

Drop / recluster families spanning >= 4 chromosomes (`dispersion.py`), designed on HG002 (94 TE-glued WRONG pairs) + KB3781. Bars: pair F > base, precision >= base + 0.05, WRONG <= base/2, sens >= base - 0.02, family F >= base - 0.03. **RESULT (register 1470): FAIL on both held-outs**: A119b primary pair F 0.783 -> 0.782, WRONG 87 -> 67; OR6737 0.307 -> 0.299, WRONG 5 -> 3. TE-glued families are an HG002 phenomenon. Opt-in tool; +0.05 absolute precision bar unreachable from 0.81. (full: archive/2026-11/LEDGER_O1.md:8371; result :8406)

## PREREG SKETCHFAM (2026-10-29, user "lets try the sketch"): exon-only FracMinHash families instead of the minimap2 all-vs-all, held-out on A119b (CHM13) and OR6737 testis

`sketchfam.py`: exon-only k-mers (k 21, scale 20), hashes in > 20 loci dropped, >= 5 shared and containment >= 0.2, overlapping loci rejected, connected components; dev HG002 + KB3781. **RESULT (register 1471): FAIL on both held-outs**: A119b pair F 0.783 -> 0.561 (TSPY 34/34 -> 0/34: the 20-locus cap erases large genuine families), OR6737 0.307 -> 0.290. Opt-in tool. **SKETCHFAM-2** (dev note, all four samples now dev): chromosome-dispersion k-mer cap (`--max-chrom 3`) restores TSPY; MCL vs components trade precision/recall; needs a fresh held-out. (full: archive/2026-11/LEDGER_O1.md:8429; result :8456; SKETCHFAM-2 :8477)

## SCORING CHANGE (2026-10-29, user): pair precision and recall only; no family F, no bipartite family matching, no F; pseudogenes not a category

From here on the readout is `bench/denovo_eval/pairscore.py`: Compara pairs among expressed protein-coding genes, precision and recall only; famtruth family F, pair F and mainstay lines are diagnostics, not bars. The section then logs the sketch-engine development (all dev unless stated):
- **SCORING REFINEMENT**: `pairscore.py --noncoding` judges non-coding pairs by sequence (shorter covered >= 50% at >= 85%), else DOMAIN-ONLY. (:8506)
- **GORILLA TREE GRID** (`--cluster tree`, UPGMA cut on containment): k 17 scale 10 chosen; smaller k reaches only older paralogs; containment is one precision/recall axis because families are nested. (:8517)
- **MINIPROT MEMBER RECOVERY**: FAIL as an add-on on both gorilla sets (older relatives, retrocopies in hosts; nc precision falls). (:8577)
- **PER-FAMILY CUT, TRUSS, BRIDGES**: FAIL on gorilla dev; relative/gap cuts inert, triangle rules kill 2-gene families; recall ceiling 0.69; intron sharing does not separate Simiiformes from older layer. (:8602)
- **INTRON RESCUE** (`--weak-cont 0.1 --intron-cont 0.2`): first rule moving recall without precision (KB +0.06, OR +0.02 recall), nc precision -0.02. (:8652)
- **MAINSTAY DIAGNOSIS on A119b**: TBC1D3 duplicate families (two locus types), NPIP pieces (NPIPB3/B12/A3/B7 have no full-length locus), GOLGA6L fragments. (:8672)
- **MAINSTAY SCORE + FAMILY-OVERLAP MERGE**: `mainstay_score.py`; overlap merge changes only family counts, costs precision; not adopted. (:8694)
- **HYBRID sketch + base alignment edges**: FAIL on gorilla dev (unspliced locus PAF admits intronic repeats). (:8714)
- **GUIDED-MODE SKETCH**: Jaccard loses everywhere; guided containment is a sound mainstay catalogue (107/123, 0 unrelated). (:8724)
- **SPEED**: exact vectorised pair counting, byte-identical edges; KB3781 106 s -> 45 s; A119b > 1.5 h -> 3.6-5 min. (:8743)
- **MMSEQS CLUSTERING**: FAIL on both gorilla sets (unrelated false pairs); mmseqs stays a judge. (:8771)
- **SEQUENCE-DEFINED FAMILY SCORE**: `seqscore.py`, truth = expressed gene pairs aligning >= 50% at identity tier T (no Compara); real gap is recall at >= 98% (KB half missed), caused mostly by the chromosome cap. (:8790)
- **CHROMOSOME CAP REMOVED**: structural RNA excluded; no cap raises gorilla recall (KB @0.98 0.529 -> 0.765); + intron rescue best (0.843 / 0.910). (:8813)
- **A119b WITHOUT the cap**: human needs a cap (no cap precision 0.62 -> 0.47) but cap 20 = no cap on gorilla and beats base on human; ONE arm for both: `--cluster tree --k 17 --scale 10 --root-cont 0.15 --min-cont 0.2 --max-chrom 20 --weak-cont 0.1 --intron-cont 0.2`. (:8842)
- **PREREG SKETCH-SEQ-HG002: PASS on all five bars** (one run, audited; B1-B5 recall @0.85/@0.98, precision @0.85, precision @0.98 >= base - 0.02, gene-averaged recall @0.98); secondary Compara P/R 0.426 / 0.345 -> 0.505 / 0.377; gene-averaged precision goes down. Declared the family-engine candidate to ship (later not shipped; see Current state). (:8871)
- **GUIDED MODE with the final arm**: mainstays 110/123, 0 unrelated; sequence recall @0.98 0.973; GOLGA8 / GOLGA6L / SPATA31 split into real sequence sub-groups. (:8905)
- **STRONG-EDGE MERGE**: FAIL as a general rule (joins NPIP to PKD1, BOLA2/SMG1); alignment-based merge likewise; no protein extension for family definition. (:8917)
- **TWO-LEVEL idea** (smaller k upper level): would give a GOLGA superfamily; not built; decision left to user. (:8947)
- **GOLGA6L anatomy**: proximal vs distal groups share no nucleotide alignment, only 42-53% protein identity; GOLGA6L is a protein family name over two nucleotide families. (:8959)
- **SD CATALOGUE CHECK**: CHM13 SD90 draws the same GOLGA6L and SPATA31 sub-group lines as the sketch (no proximal-distal, no A-D pair); finds GOLGA8 13/13 (sketch 10/13). Sub-group splits are a property of the genome. (:8975)
- **SD-PROJECTION MERGE**: with the >= 50% fraction rule, opt-in (`sd_merge.py`), completes GOLGA8 13/13 in guided without touching NPIP; not part of the registered arm. (:8984)
- **FAMILY VARIATION GRAPHS**: the variation graph stays the thesis model; sketch only finds members; `family_graph.py` (seqwish) collapses each family into one connected graph; seed variants (spaced, translated) do not change containment. (:9012)
(full: archive/2026-11/LEDGER_O1.md:8497)

## Protein layer over the sketch families (2026-10-30, report; `bench/denovo_eval/protlayer.py`, `~/runs/protlayer/`)

Layer 2 = translated all-vs-all (mmseqs) of members; merge families by protein hits. Ungated: superfamily hairball (Compara precision 0.27 -> 0.09-0.10); NOT adopted as a merge. **Specific protein reconciler** (`--rbh --same-chrom`, 32-arm grid on guided A119b): id 0.7 rbh+chr merges 28 same-root local clusters at Compara precision -0.002 / recall +0.011; candidate second rule, unregistered (proposed arm id 0.7, cov 0.5, rbh + chr, frac 0.5). **On the de novo families (register 1480): NULL** — every arm leaves the 14 mainstay rows byte-identical; proposed arm merges 6 pairs (Compara +2 hits, +1 false); not pursued, no HG002 registration. (full: archive/2026-11/LEDGER_O1.md:9043; specific :9072; de novo :9103)

## PREREG SEMIGUIDED (2026-10-30): semi-guided families = de novo families + silent copies confirmed by the seed locus (user: "Lets focus on implementing the semi-guided mode, notice that in the artifact for CDY2B for example, if we used the secondaries we would be able to find at least part of the locus")

Every `semisec` candidate (Amendment 1: silent or not) overlapping no copy, confirmed by splicing the seed locus's exons onto it (cov >= 0.5, id >= 0.7), attached to the seed locus's family (`semifam.py --all-candidates`). Dev A119b (CDY2A/CDY2B attach); held-outs KB3781, OR6737 with Yoo mainstays. **RESULT (register 1472): PASS on both held-outs, all four bars** (B1 Yoo recalled_annot 14 -> 19, 15 -> 18; B2 unrelated 0; B3 precision 0.333 >= 0.321, 0.456 >= 0.409; B4 recall 0.609 -> 0.636, 0.677 -> 0.756). Honesty note: held-outs scored three times (first under a strand bug: FAIL B3); PASS after a FAIL under a bug is weaker. Rust port `src/semifam.rs` byte-identical. **Not a default; auditor pending.** (full: archive/2026-11/LEDGER_O1.md:9130; result :9161)

## PREREG REPRULE-2 (2026-10-30): the locus representative that agrees most with StringTie and FLAIR (user: "ensure that the assembly part is as similar as possible to stringtie and flair ... the chosen representatives are still different so it's time to measure and pick a new default")

Stage 1: reads x (junctions + 1) is the only rule above most-reads on StringTie exact / FLAIR member / RefSeq exact in both samples (A119b 0.593/0.670/0.503 vs 0.565/0.599/0.459). Stage 2 `mcl_families --representative reads-x-junctions`, held-outs KB3781 + HG002. **RESULT (register 1473): FAIL** — HG002 NPIP 11 -> 10 (B1); dev A119b recall 0.741 -> 0.679, NPIP 18 -> 15. Stays opt-in; `most-reads` remains the default. (full: archive/2026-11/LEDGER_O1.md:9197; result :9224)

## STRINGTIE NON-FLOW RULE REVIEW + THESISVAULT RULES (2026-10-30; reading only, no run yet; `~/stringtie_modded/rlink.cpp`)

Reading notes: six StringTie -L non-flow rules (greedy abundance-ordered grouping, measured as a pick rule: dead, StringTie exact 0.462 / 0.719 < most-reads; trustworthy ends; `compatible_long_core`; site scan; output-stage pairwise filters; terminal read exon cut) and ThesisVault rules (IsoCon nearest-neighbour test, Iso-Seq3 asymmetric end tolerance, Pal 2026 family-relative support, Yoo 2025 read-set Jaccard edge, FLAIR 10 bp junction correction). No run. (full: archive/2026-11/LEDGER_O1.md:9242)

## PREREG REPRULE-3 (2026-10-30): where the NPIP split of r1473 actually happens (diagnosis seen; one arm registered before any run)

Diagnosis: under reads-x-junctions MCL splits a near-clique of NPIP by an identity gradient (A vs B subfamily); `--inflation 2.0` not the lever. Arm S = identity saturation `w' = min(1, w/0.95)` via `--extra-edges-only`. **RESULT (register 1474): FAIL** — inert on base (0.725 / 0.741), rxj NPIP still 15 in 7; the split comes from MISSING cross edges, not weights. (full: archive/2026-11/LEDGER_O1.md:9259; result :9275)

## NEXT TARGETS (registered 2026-10-30, user decision): one binary, one object, stage escape hatches; prune as it comes

Plan: (1) single `rustle` binary with stage subcommands over one run directory; (2) escape hatches `--to` / `--from` / `--skip` (the `--to families` hatch answers the advisor's "loci before assignment" question); (3) acceptance = byte-identical products (`tools/identity_check.sh`); (4) prune OPT-IN modules no default or open prereg uses under an `archive` feature. Order 2, 1, 4. Subsection **REPRULE-3 chain gate, dev result**: the long NPIPA loci fail the interval gate (exons outside the homologous unit); `--chain-gate 3` on the base representative keeps every mainstay and gains pairs on dev (0.726 / 0.741 -> 0.734 / 0.759); NPIP not rescued under rxj. (full: archive/2026-11/LEDGER_O1.md:9283; chain-gate dev :9299)

## PREREG CHAIN-GATE (2026-10-30, registered before any held-out run): `--chain-gate 3` on the shipped representative

Base + `mcl_families --chain-gate 3`, held-outs KB3781 + HG002 with REPRULE-2 bars; Amendment 1 adds arm U `--chain-gate-union`. **Arm C RESULT (register 1475): FAIL on both held-outs** — KB3781 precision 0.341 -> 0.321 (B2), HG002 NPIP 11 -> 10 (B1); stays opt-in. Arm U dev: below arm C (0.727 / 0.750), not held out; opt-in. (full: archive/2026-11/LEDGER_O1.md:9321; arm C :9335; arm U :9349)

## PREREG LOCUS-ABSORB (2026-10-30, registered before any run): fold dominated overlapping gene_ids into one locus

`mcl_families --absorb-loci --absorb-ratio 0.1 --absorb-frac 0.5`: fold a same-strand overlapping gene_id with <= 0.1 x the reads and >= 0.5 of its exon bases inside a larger locus (112,404 overlapping gene_id pairs share no junction). Dev A119b: precision 0.726 -> 0.796, recall 0.741 -> 0.729, missing the registered dev gate (>= 0.731) by 0.002; the 26 lost true pairs are KRAB-ZNF pairs recalled only through absorbed hub loci. **Post-hoc amendment (user, 2026-11-01)** to run the held-out anyway. **Held-out OR6737 (register 1476): PASS B1 + B2, by equality** (mainstays and pairs identical: 0.429 / 0.239). Proposed as families-stage default (null effect on held-out; dev gate amended post hoc); auditor first, user decides; nothing changed in binaries. (full: archive/2026-11/LEDGER_O1.md:9356; dev :9376; lost pairs :9384; decision :9397; held-out :9405)

## TOOL COMPARISON ON MULTI-COPY LOCI (2026-10-31, descriptive, no prereg; `bench/denovo_eval/toolcompare.py`, runs `~/work/jfim/runs/toolcmp/`)

Per mainstay/YAG copy, modeled and junction precision. A119b (127 expressed copies): rustle 123 / 0.46, StringTie 107 / 0.47, FLAIR 98 / 0.35 (gap = TSPY 34/35 vs 19 vs 12, CT47A); OR6737 (30): 27 / 0.58, StringTie 28 / 0.59, FLAIR 23 / 0.55. Ahead on tandem arrays, tied elsewhere; junction quality at StringTie's level. (full: archive/2026-11/LEDGER_O1.md:9419)

## PREREG MS-VS-AS (2026-10-31, registered before any run): does minimap2's ms ranking disagree with O2's AS ranking often enough to matter?

Rule: port `--tie-ms` only if disagree AND distinct place AND outside the 10-AS window >= 5% of family-scope reads with >= 2 records on either sample. Amendment 1: the **Eichler 10-AS rule is a standing test** on every read-placement result. **RESULT (register 1477): bar met on both gorilla samples** (OR6737 10.2%, KB3781 5.8%) -> `--tie-ms` worth porting; 99% are spliced ms-best vs unspliced AS-best retrocopy rivals. Eichler rule decides 80.0% / 83.6%, agrees with primary 87.5% / 93.1%, retrocopy pick 11.9% / 6.8% of decisions. A119b appended later (in the ISOCON section): 1.9%, below the bar; verdict stands (rule was "either"). (full: archive/2026-11/LEDGER_O1.md:9434; result :9479; A119b addendum :9518)

## PREREG ISOCON-GORILLA (2026-10-31, registered before any run): IsoCon as the "truly de novo" comparator on single- and multi-copy targets

IsoCon 0.3.3 on alignment-defined pools, OR6737; 4 Yoo families (49 copies, 30 expressed) + 30 single-copy genes. **RESULT (register 1478): IsoCon WORSE in both classes** (multi-copy 20 vs ours 27 modeled; single-copy 27 vs 30, junction precision +0.05); a third to a half of its multi-copy candidates are `placed_ambiguous`. "Truly de novo" buys nothing on copy recovery here. (full: archive/2026-11/LEDGER_O1.md:9457; result :9502)

## PREREG FAMASM (2026-10-31, registered before any FAMASM number): multi-copy assembly mode on top of the normal assembly

`bench/stgaps/famasm.py`: per family pool reads (primary + secondaries within 2% of best AS), lift chains to the representative, family chains by the pass-1 floor (>= 2 molecules over all copies), admit at a copy only with reads certified to that copy (>= 2 whole or compatible). Amendment 1: single-record reads certified by placement; Amendment 2: rep junctions outside the locus are private. Bars M1 multi-copy recall +1.0, M2 precision >= BASE - 0.5, M3 single-copy |d| <= 0.1, M4 families. **HG002 held-out RESULT (register 1479): FAIL** (M1 PASS 11.3 -> 12.7, M2 FAIL 14.0 -> 8.5, M3 FAIL precision 19.2 -> 18.3, M4 PASS); not adopted, prototype stays in `bench/stgaps/`. Dev A119b (reported): same shape, M1 +0.7 FAIL, M2 12.7 -> 8.0. Next: stricter certificate on A119b only (HG002 spent). (full: archive/2026-11/LEDGER_O1.md:9525; result :9565; dev :9592)

## Assembly ends, families non-regression, annotation rules (2026-11-01..02, register 1485-1502)
Pre-registrations: `archive/2026-11/PREREG_assembly_ends_2026-11-01.md` (E1-E7), `PREREG_readthrough_detect_2026-11-02.md`, `PREREG_truth_split_collapsed_2026-11-02.md`, `PREREG_fused_loci_2026-11-02.md`, `PREREG_shared_exon_split_2026-11-02.md`, `PREREG_final_subchain_2026-11-02.md`, `PREREG_sqanti_filter_2026-11-02.md`, `PREREG_end_artifacts{,_shared,_pas}_2026-11-02.md`. Held-out = A119b chr11 + chr7 vs CAT (gffcompare) for the ends; A119b + OR6737 genome-wide pairscore for families.
- E1 `--polish-tss rescue` PASS -> default (chr11 chains 1,695 -> 1,702; precision 16.577 -> 16.505). E2 sub-chain rescue STOPPED (target class 2.1%, no PAS-proven 3'-short end). E3 median 5' ends PASS (starts within 50 bp 40.6 -> 64.6% chr11). E4 3' UTR extension PASS (base sensitivity +2.1 / +2.4). E5 single-exon cluster models PASS (base precision 25.8 -> 43.0 chr11). E6 families non-regression of E3+E4+E5 FAIL on A119b (pair precision 0.726 -> 0.700, recall 0.741 -> 0.714; single-exon copies 3,299 -> 217). E7 E3+E4 only FAIL (recall 0.741 -> 0.710; families split 1,565 -> 1,586). Resolution: assemble writes `PREFIX.final.gtf` with the layers; families keep `PREFIX.gtf`.
- Readthrough: structural rule recall 0.757 GENCODE / 0.828 RefSeq (bar 0.80) FAIL; labels = gene record only (r1492-1493). Truth split of 188 CAT collapsed genes into 444: recall -0.0049 (bar 0.005) opt-in (r1494). Fused-loci annotation-guided cut PASS: chr17+19 fused 106 -> 55, genome-wide 547 -> 282, precision 0.726 -> 0.761, recall -0.0096 (r1495). Shared-exon split FAIL (loci 60 -> 101, precision 0.698 -> 0.344; private units 60 -> 44, pairs unchanged; r1496-1497).
- Final sub-chain layer PASS and default in `final.gtf` (r1498). SQANTI3 rules filter: chain precision 19.1 -> 25.4 but = -6.4% FAIL (r1499). Own RT-switch / intrapriming detectors match SQANTI3 flags 100% in 33 s (vs 39 min); rules E (-2.3% =), J (-0.83%), K (-0.69%) all FAIL the 0.5% bar; `utilities end-artifacts` kept (r1500-1502).

## IsoCon salvage, tandem cuts, seeding pool, family panel, projection layer, base-mode family misses (2026-11-03..04, register 1515-1525)
Pre-registrations in `archive/2026-11/`: `PREREG_isocon_only_2026-11-03.md`, `PREREG_tandem_cut_2026-11-03.md`, `PREREG_tandem_cut_neighbour_2026-11-03.md`, `PREREG_seed_unspliced_primary_2026-11-03.md`, `PREREG_seed_unspliced_support_2026-11-03.md`, `PREREG_family_panel_tools_2026-11-03.md`, `PREREG_project_models_2026-11-03.md`, `PREREG_project_models_v2_2026-11-03.md`, `PREREG_family_mode_tied_seed_2026-11-04.md`, `PREREG_base_family_misses_2026-11-04.md`, `PREREG_strand_vote_dup_zmw_2026-11-04.md`. Tools: `bench/famcmp/family_panel.py`, `score34.py`, `cut_fused.py`, `bench/project/`.
- IsoCon-only transcripts (r1515): our candidate_copies chain with every consensus (`--all-clusters`) recovers 9.6% of the 15 isoforms only IsoCon finds at 58% artifacts: FAIL, closed. Panel table now carries `expr_recall` (reference isoforms with >= 2 exact primary reads, 936 on chr20: ours 67.1%, IsoCon 65.6%, FLAIR 59.2%, StringTie 53.3%) and `expr_recall_filt` after the SQANTI rules filter (60.7 / 61.9 / 54.5 / 48.5; the RTS flag costs us 60 true isoforms). Of the 308 expressed isoforms we miss, 71 with >= 10 exact reads are sub-chains the pass-2 collapse removes; rescuing contained chains by read floor has precision 0.06-0.08 vs 0.19 for what we emit, so no assembler lever remains there.
- Tandem cuts (r1516 Rule T, r1517 Rule T2 neighbour blocks): reach 5-6 of the 69 spliced tandem loci, 67% of cuts inside one gene: both FAIL. Class = 5-7 near-identical spliced pairs + 4 intronless blobs (USP17L 41 genes / 265 kb, TAF11L, OR10A2/4, SPINK7/9).
- Seeding pool (r1518, r1520): the good pool's near-tied UNSPLICED secondaries chain by span into array-wide single-exon blobs (USP17L: 3,567 "reads" over 130 primaries; 61 models >= 20 kb with a zero-primary run >= 3 kb genome-wide). `RUSTLE_SEED_UNSPLICED=primary` (drop them) removes every blob, pairs 0.728/0.742, fused 545, spliced models identical, but CAT single-exon genes hit 757 -> 342 (the r1060 pseudogene fill-in): FAIL. `=support` (they add support to a primary cluster, never seed or bridge) keeps 658 of 757 (76% of the loss recovered), blobs gone, pairs 0.726/0.735: FAIL on bar 4 (0.869 < 0.98) and FSM by one. Both opt-in; the default trade (99 pseudogene copies vs 61 blobs) is the user's call.
- Family panel across tools (r1521, descriptive, 420 primates-level truth families, A119b base mode): expression judged under P (>= 2 primaries) and T (>= 2 placements incl. tied secondaries). Ours hits 0.69 of the 87 T-only members (StringTie 0.00, FLAIR 0.05), 0.94 of MAPQ0-dominated members (0.83 / 0.84), fuses 49 families (59 / 69), but 0.82 of ordinary members vs 0.89 (3-9-molecule floor/polish), so StringTie/FLAIR complete more families (229/239 vs 204 under T). Decision (user): base mode need not win families; the family layer does (projection r1519, then tied-seed targets and double-supported splits as v2).
- Projection layer (r1519, FAIL as registered; chr20 own tier FSM+ISM 66.7% / artifact 11.1% vs base 16.3% / 39.9%, but 99% of projections fall in the unannotated chr13 hotspot with motif rate 0.35): each family's multi-exon models projected as spliced templates onto the other members (catalog unchanged, structure borrowed, reads never pooled), tiers own / tied / none at the target; `bench/project/project_models.py`. Four runs hit the 30-min wall (templates genome-wide; all-vs-all; repeat minimizers -> `-f 2000`; BAM pass per chain -> per target).
- Projection v2 (r1522, FAIL both arms): hotspot-out arm H shows the method's true shape: own tier 63 models, FSM+ISM 41% vs base 12%, but reach 0.41 and the none tier (1,004 models, 35% artifact) is unverifiable annotation transfer; the T >= 2 gate was vacuous (892/892). Closed: projection is a per-member support tool, not a model emitter; next family-layer row = tied-seed targets + fused-model cuts.
- Family-wide mode (r1523, FAIL as registered): tied-seed adds 0 of the 3 members reachable only by a tied chain (the good pool already seeds tied spliced secondaries; 65 base hits rest on them); cutting fused models at truth-member boundaries lowers fused families 49 -> 36 at zero hit cost but is annotation-guided. The misses are base-mode: of 222 missed expressed members only 34 carry a >= 2-read local chain; causes found = single-exon strand placeholder (RUSTLE_READ_STRAND opt-in) and a pass-1 drop of 2-read canonical chains (PRAMEF8). Tied-seed closed.
- Base-mode family misses (r1524, FAIL as registered on fused and chain-identity bars): the read-orientation strand vote alone is clean (every spliced chain identical, complete families 204 -> 209, ordinary hit 0.819 -> 0.829, 3,103 single-exon models now on the minus strand, 2,993 fewer single-exon models unexplained) and is the next default-on candidate; the pass-1 coordinate-duplicate admission rescues 24 of 34 family misses (+20 complete families) but admits 68,749 chains at 1% CAT match, so it must be narrowed; the catalog-family cut barely moves fused families (49 -> 47).
- Strand vote shipped (r1525 part 1 PASS): read-orientation strand for single-exon models is the assemble default; every spliced chain identical, +25 CAT single-exon gene hits, +5 complete truth families. Duplicate admission (r1525 part 2 FAIL): coordinate duplicates never share a ZMW in A119b, so the different-ZMW rule is vacuous and the 68,749 admitted chains stay at 1% CAT match; RUSTLE_DUP_ADMIT stays opt-in; the sequence-based duplicate test was pre-checked and closed: over a shared 1.5 kb window coordinate duplicates (391/400) and independent same-chain reads (397/400) are equally identical, so HiFi sequence carries no duplicate signal. Compatible-chain admission (r1526, FAIL as registered): the goal's rule (compatible reads pooled into the longest chain) is the existing `--direct-transfrags-union` (r1266); step 0 ceiling 40 of the 97 one-chain members; arm U reaches 19, complete 209 -> 219 and ordinary 0.829 -> 0.847, but fused 49 -> 54, base chains kept 0.856 (contained chains dropped by the sub-chain rule), new chains' CAT match 0.022, wall 2,881 s. Opt-in; closed. Of the 97, 24 have a single spliced read and 23 differ only by <= 10 bp junction offsets (`--read-junction-snap`). Junction snap (r1527, FAIL as registered): RUSTLE_READ_JUNCTION_SNAP=10 genome-wide moves 9.4% of base chains, complete 209 -> 208, fused 50, new chains' CAT match 0.032, reaches 1 of the 23 predicted members; opt-in, closed. The one-chain class is now closed on every read-level rule (duplicates, pooling, snapping); what remains for those members is single-read evidence, which no floor-2 assembler admits.
- Provenance route pre-check (2026-11-04, not registered; ~/r1528/step0_prov.py): assemble each member from reads O2 certifies to it. On the 28 missed T-only members (201 tied reads, A119b a1.assign of 2026-10-31): assigned to the member 4, to a sibling 15, tied 16, ambiguous 16, never in the assign table 146 (136 are secondaries at the member whose primary lies elsewhere, 119 on the same chromosome; 15 of the 28 members overlap no catalog copy at all). Ceiling today: no member reaches a 2-read certified floor. The bottleneck is upstream of O2 and of the assembler: the catalog misses the member (sketch / family discovery) or the read pair never enters the family's tie window. A certified-read assembly row is premature until the pool admits those reads.
- Shared-exon gate vs coverage (r1528, FAIL as registered): 631 of the panel's 1,136 expressed members are catalog copies; 299 uncovered members have a base model and are sketch singletons whose sibling pairs fail the 0.60 shared-exon gate (120 of 137 aligned pairs below it; 68 pairs do not align at all: PRAMEF/GAGE/PAGE/XAGE). Lowering the gate (new switch RUSTLE_MIN_SHARED_EXON_FRAC) to 0.45 / 0.30 raises coverage to 676 / 713 and pair recall to 0.753 / 0.779 but pair precision falls 0.726 -> 0.670 / 0.662 (new pairs mostly domain-level relatives). Gate stays 0.60. This is the family layer's coverage ceiling on the panel: the members it cannot reach are either DNA-divergent families or pairs the gate rejects for the same reason it rejects domain relatives.
- Transcript-level edge (r1529, STOPPED at step 0): pairscore's spliced-sequence rule reaches 30 of the 238 singleton expressed members with a sibling locus (bar 60) and passes 40% of the base's false coding pairs as registered (8% once the 205 same-symbol CAT-duplicate pairs, identical copies Compara cannot name, are excluded). The rejected sibling pairs are diverged at the transcript level too; no edge criterion on sequence reaches them. Family-layer ceiling on the panel recorded: 631 of 1,136 expressed members. Scorer note: about a third of pairscore's 'false' pairs are identical-copy duplicate ids, so the catalog's pair precision is understated.

