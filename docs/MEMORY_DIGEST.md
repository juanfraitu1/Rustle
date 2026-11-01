# Memory digest — the long form of the memory index (compacted 2026-11-01)

Compact, current version. The verbatim text (172,911 bytes) is `docs/archive/2026-11/MEMORY_DIGEST.md`; each entry ends
with a pointer to its line there. This file archives memory-index lines from 2026-07 to 2026-10-25; it is NOT a ledger
(measurements: `docs/LEDGER_O1.md`, `docs/LEDGER_O2.md`, `docs/NEGATIVE_RESULTS_REGISTER.md`). Older entries name
superseded catalogs (/494, /1415, 627/2,019 all superseded; see `NUMBERS.md`, `o1_catalog_provenance.md`) and the old
WSL machine paths (`/mnt/linuxdisk`, `/mnt/c/Users/jfris/Desktop/Rustle`); today's environment is Fritz (user memory `fritz-env.md`).

## Current state (2026-11-01)

Derived from register rows 1469-1484 and the user's auto-memory; nothing here is re-measured.
- **Defaults unchanged** through rows 1469-1484. Every family rule tested 2026-10-29..11-01 FAILED or stays opt-in.
- **Family-stage rules, all opt-in / not default:** U1-FRESH-HG002-INTRON r1469 FAIL (pair F 0.411 < 0.413, precision
  0.415 < 0.478); DISPERSION r1470 FAIL (A119b pair F 0.783 -> 0.782; OR6737 0.307 -> 0.299); SKETCHFAM r1471 FAIL
  (A119b pair F 0.783 -> 0.561, TSPY 34/34 -> 0/34; sketch tree PASSED HG002 per memory, not shipped);
  SEMIGUIDED r1472 PASS both held-outs (KB3781, OR6737) but "not a default; auditor first"; REPRULE-2 r1473 FAIL
  (HG002 B1, NPIP 11 -> 10); REPRULE-3 arm S r1474 FAIL (inert); CHAIN-GATE r1475 FAIL (opt-in only);
  LOCUS-ABSORB r1476 PASS B1+B2 by equality on OR6737 (null effect on gorilla).
- **Placement / O2:** MS-VS-AS r1477 bar met (OR6737 10.2%, KB3781 5.8%) -> `--tie-ms` worth porting; TIE-MS PORT
  r1483 FAIL (B2, B3 on OR6737 2.6%) -> `--tie-ms` stays opt-in, next = HG002 truth with B2 on final placement.
  O2 LEVERS r1481: full runs gave nothing (U cancelled at 19 h, C vanished), amended to a 57-family panel
  (`~/o2levers_fast/`): arm U (`--union-certificate`) FAILED the 57-family fast panel (assigned 59 < 63; rejections become ties, 311 molecules go to outside loci), stays opt-in; arm C (candidates) still to rerun on the panel. PSV VCF REPORT r1484 (report): bcftools calls 54.6% of O2's columns as registered,
  93.9% at `--mapq-reset 60` (post hoc).
- **Assembly / comparisons:** IsoCon r1478 WORSE in both classes; FAMASM r1479 FAIL on HG002 (multi recall
  11.3 -> 12.7 but precision 14.0 -> 8.5), not adopted; SQANTI single-copy r1482 (370 genes, report): ours FSM+ISM
  30.2% / artifact 5.0%, StringTie, FLAIR, IsoCon all WORSE by the rule. Protein reconciler r1480 NULL (no prereg).
- **Open / next (user memory):** one `rustle` binary + one run directory with stage hatches (`--to/--from` done,
  `--skip` to do), prune OPT-IN modules under an `archive` feature, acceptance = byte-identical products; advisor's
  "loci before assignment" = `families` stage product (`PREFIX.fam.copies.tsv`, `PREFIX.fam.loci.gff3`).
  FAMASM follow-up = stricter certificate on A119b. O2 lever = candidates/semisec on A119b, HG002 + 1-10 AS near ties.
- **Standing rules:** score families by pair precision/recall only (`pairscore.py`); the variation graph is the
  thesis model (sketch = member discovery only); Eichler 10-AS rule is a standing test on every placement result.
- SPENT samples (memory): KB3781, OR6737, A119b, HG002.

## Top-level index lines (2026-07 .. 2026-09-01, preamble, archive lines 1-120)

Unheaded block of ~50 index lines; highlights, newest first (full: archive/2026-11/MEMORY_DIGEST.md:1)
- Seeded mode 09-01: multi-seed 30/31 (CDS 14/31, genomic 25/31), a RECALL instrument not a definition. Cross-copy borrowing DEAD (borrow for assignment, never assembly).
- TSS/TES/UTR objection: premise FALSE (median end difference 1,190 bp; 3' flank 1,672 vs 5' 177 bp = 9.45x); concede coverage charged on shorter only (59.97% of edges would fail 0.50 on the longer).
- Advisor readiness audit: median edge identity 0.8287, only 5/121 families externally adjudicated; 47.7% tandem cases undetermined.
- minimap2 artifact telltales 3 real / 4 guarded / 13 refuted; `RUSTLE_JUNCTION_MAJORITY` measured (needs chr16 arm; later run, see 09-19 section).
- Engulfment = strand placeholder (58/60 are ('+','-')); neither fix ships. Spanning reads: O2 keep-current (proven); 71 copies engulfed, 34/121 families.
- Vacuous instruments are a class: 40 scripts audited, 10 confirmed FALSE PASSes, all 10 repaired; "show an instrument can FAIL before its PASS is evidence".
- PSV objections answered, do not build the VCF (genome_pos is not a locus); `RUSTLE_XFAM_RECONCILE` demotes 111/111. (Superseded in part: r1484 PSV VCF report, 2026-11-01.)
- k+1 MEC NO-GO; ML has no slot (labels); shared-read E_c edges REFUTED as a definition tier (fusion ~72% at any depth).
- O1 gap = primary-flag starvation: `RUSTLE_FLAGFREE_SITES` 14/31 -> 26/31 (opt-in). Node floor 2 SHIPPED (NPIP 12/31 -> 14/31). Flag audit: nothing qualifies to flip (124 flags).
- Haplotype CNV proven in one animal: GOLGA6L7 6 mat / 7 pat; `_pri` = 16 PAT / 9 MAT mosaic. O3 must read within-locus (1.27x vs 109.96x).
- Discover-then-define is CIRCULAR; never define membership by gene symbol; E_r identity clause DOES fire (nm/bl).
- Standing scope (08-25, user): TBC1D3 only + fibroblast-matching-T2T only; 627-family catalog is OR6737 testis (out of scope then).
- 08-23..08-25 ledger §4l-§4u: cross-substrate 87.06% of E_r edges recovered (ARI 0.7064); orientation guard filters an unmeasured strand field (5,928 single-exon reps all '+'); read-strand fix opt-in (break-even 627 -> 630); representative selection closed (0/215); stub defect no route; WGS gives no CN baseline.
- 08-11..08-21: 26.8% shipped copies single-exon stubs (541/2,019); SD-gap/primary-flag headlines refuted (C5 OFF); reads cannot fix the coverage hole; genome-anchored repeat gate (veto only); coverage hole not fixable; definition's named hole = scale-free coverage denominator; O2 reassignment truth UNDERPOWERED; O1 false-omission 9/162 = 0.0556; blind-spot audit: restate O2 target as AS near-ties (21.75%), defend O2 on abstention; O1 reach 22/40; O1 false-merge 2/150 = 1.33%; EM changes 0/3,081 decisions.
- Register then 631 entries (later 814, now >1484).

## O3 / detection

Haplotype deficit 8/9/9-vs-5/6/8 do NOT quote (record `-p`/`-N` with any copy count); unmapped-read route works only for UNIQUE sequence (STON1+GTF2A1L single detection, unconfirmed); whole-genome excision: ABSORBED 64.2% / ORPHANED 33.3%, S2 TPR 0.2703 / FPR 0.0200; DNA-side collapse zero at 816/817 (predicted, underpowered); literature: RNA niche empty; O3 candidates retracted (immunoglobulin/SHM; mGorGor1 = KB3781 male); gate 5 near-unpassable; O1 = missing from ANNOTATION, O3 = missing from GENOME. (full: archive/2026-11/MEMORY_DIGEST.md:121)

## O1 definition / mechanism

Family graph object: nodes = genomic intervals, edges = contiguous homology on assembly sequence, RNA sets nodes; families not one graph class; lambda per family SHIPPED. P1 seed-invariance a theorem; coverage does the discrimination; DNA and RNA give identical partition 7/7 on same loci; NPIP does NOT split on RNA (identity never fails 171/171); boundaries and edges coupled (recall 0.698 -> 0.788, precision 0.709 -> 0.269); primary-flag objection answered (95.9% invariant -> ABSTAIN). (full: archive/2026-11/MEMORY_DIGEST.md:141)

## Benchmarks / validation

Soto audit supersedes earlier Soto numbers (our merge was the bug; ARI 0.514 -> 0.757, precision 96.19%); A119b IsoSeq and Soto 2025 are HUMAN; O-2 closed (E_r union only at `refine`, never called by O1 catalog); joint DNA/RNA = property, not definition (length-matched DNA-only arm recovers 91/91). (full: archive/2026-11/MEMORY_DIGEST.md:163)

## O2 supporting

K=0 wall = READ not COPY; never quote 0/386 as a rate; PSV lever = the gate not scoring. (full: archive/2026-11/MEMORY_DIGEST.md:173)

## Infrastructure / defaults / dead-ends

Old WSL storage/mount notes (stale; Fritz now); defaults list (psv filter, tspy sim, MAPQ=0 = tie); dead ends in the register (FamilyGraph, k-mer Jaccard DEAD; vg TIES linear); many `project_*.md` not indexed. (full: archive/2026-11/MEMORY_DIGEST.md:177)

## Advisor / framing / prior art

Advisor Stefan Canzar (2016 multimapping = max facility location; wants clean graphs + theorems, dislikes 1/k); future avenues gated; prior art: Soto 2025, IsoCon closest (no allele-vs-copy), RNA niche empty, SDA same K=0 floor, Eichler TBC1D3 VG failed -> phylo. (full: archive/2026-11/MEMORY_DIGEST.md:187)

## User

Bioinformatics dev; gorilla great-ape pan-transcriptomic thesis under Canzar. (full: archive/2026-11/MEMORY_DIGEST.md:195)

## Compacted 2026-09-05 — full text of the index lines shortened on that date

09-01..09-05 O1/O2 arc (full: archive/2026-11/MEMORY_DIGEST.md:201):
- Excision §6ff rows 711-712: NPIP copy removed -> 100% abstain, 0 silent; wall = local identity. Read-star speed §6fe: 2x (1,299 -> 654 s).
- Genomic read-star = O2 default §6fd rows 709-710: paired 35 55.7% assigned @ 9,236/9,236; pairwise read-star default §6fc rows 704-708; read-star default + unit merge §6fb row 703 (canonical `rna_units_v9`, 349 units); D3 origin certificate §6fa rows 700-702.
- Five defaults flipped (user, §6ez): fold-within-clusters, exonless-span, exonic-both-sides, units-follow-reads ON, min_size 2; PSV read filter OFF; canonical `mcl_ann/rna_units_v8`. Row 694 root-caused §6ey rows 696-698 (fold-first locus rule). Q9+G5 closed, min_size 2 §6ew/§6ex rows 693-695.
- Two O1 definitions on one truth §6ev: RESTATEMENT ACCEPTED (user 09-05): O1 = SD-core family with read-supported units; E_r/gamma = earlier definition, opt-in `gw_family_catalog`. Framing audit 09-05.
- Worst-2 §6eu rows 689-691 (D1-D3); hygiene + sweep v2 §6et; S5 §6es (ALWAYS map through `<out>.family_join.tsv`); S2 §6er; S1 truth audited §6eq (55/117 anchors annotation gaps); roadmap `docs/archive/2026-09/ROADMAP_O1_O2.md`.
- O2-8c' shipped OFF (685/686); O2-8c not accepted (684); O1-10b shipped, §6el relabelled (682/683); O2 memory fixed §6em (25.4 GB -> 0.49 GB); O2-8a REFUTED (681); O2 on NPIP cores §6ej (678-680); core refinement OFF §6ei (677, DROP arm FAILS); duplicon-first §6eh; NPIP biology §6eg (676); locus node rule OFF §6ef (674/675); O2 on NPIP at scale §6ee (673).
- Default prune 1e-9 (user, §6ed, 672; canonical `mcl_ann/rna_bp1_p9`, `gw_bp1_p9`); MCL prune size-dependent §6ec (670/671); audit46 §6eb (668/669); "more RNA in defining families" NO §6ea (667); repeat library §6dz; disagreement cells §6dy; SEDEF = orthogonal instrument §6dx.
- O3 status 09-03 (`LOC134757045` first clean candidate; MAPQ filter RETRACTED); O3 positive control GOLGA6L7 +1 paternal; O1->O2 bridge 68/68; Soto sidelined, pivot to MCL-on-annotation + RNA (user 09-03); gate G: NPIP false merges are Alu; register-first rule; nothing from 09-01 precision-lever line shipped; coverage clause keep opt-in; bug-class sweep §6cn (§6be 50% -> 40.3%); Soto adjudication P 0.9739 / R 0.8743; human negative panel 5/5; dossier `docs/ADVISOR_QUESTIONS.md`.

## Compacted 2026-09-19 — full text of index lines moved out of MEMORY.md

09-07..09-19 (full: archive/2026-11/MEMORY_DIGEST.md:264):
- Thesis = 3 objectives (O1 family definition, O2 copy assignment assign-or-abstain, O3 reference-absent DETECT+FLAG; ASJ dropped); user 09-12 families may be defined in HUMAN (never pool). Scope 09-13: NO genome-only discovery.
- Old crash/disk rules (WSL2, stale). Invariants: `-F 2308`; never `HSA_genomic.gff`; gorilla-native `GGO_genomic.gff`. Hold a substrate back; pre-register; metric traps (21); never conflate O1 and O2.
- 09-16/17: leader rule breaks nesting; TBC1D3 Guitart groups = population paralog groups (16 claims retracted); nested edge-test lattice proposed; NPIP ladder (gap = node quality); family certificates (NPIP certified, TBC1D3 not); ideal expression NO CANDIDATE; cross-species conjunct not confirmed de novo; locus-formation gap SHIPPED 94ed994e read-isoform widening k=5 (§6m0).
- 09-18: no-readthrough counterfactual (§6m1); PKD1P6-NPIPP1 real fusion; chimera policy NOT ADOPTED (r845); node cut rule harmful (r846); RNA reconstructability 5/26; locus assembly NPIP (r847); assembler widening W-2 FAILED (r848); canonical junction ceiling 209/249 (r849); junction majority chr16 arm run §6n4 (families 282 -> 290, not flipped); sim ceiling 25/26 (§6n0/§6n1); `RUSTLE_GTF_SECONDARY` 10/26 -> 20/26 (O1 only; AS-tie echo collapse REFUTED r850); clustering bakeoff §6n5 (Louvain REFUTED at DNA level r851); Soto is a refinement of us §6n7; readthroughs §6n9/§6o0 (row 852); §6o1 RETRACTED; overall F band 0.57-0.63 (never quote 0.88 alone); goal run: best definition `--min-exonic-bp 1 --min-shared-exon-frac 0.60` LORO-validated (chr5/7/21 F 0.845).
- 09-19: L3 cut 0.98 -> 0.985 SHIPPED (§6p1/§6p2); self-tuning L3 FAILS (§6p3); CDS complete, TBC1D3 5'UTR 75.3% (§6p4); L4 subfamily `w_98 >= 0.995` selective SHIPPED (§6p5); `--assemble-only` SHIPPED (§6p6); chr20 levers (§6p7, ISM collapse Pr 44.6 -> 51.9); assembly polish 28/30 cells, 30/30 unreachable (row 863).
- Soto fix default `--min-shared-exon-frac 0.30` (§6kr-§6kt); current family definition then `docs/seeded_family_definition.md` §0; edge core LCS default (register 814); O2 AS-tied gate (§6hz 0.846 hard-locus bakeoff); resume pointer `docs/archive/2026-09/PENDING_2026-09-17.md`. Do NOT re-propose rows 1222, 294, 1072.

## Compacted 2026-09-21 — full text of index lines moved out of MEMORY.md

Four 09-20 lines superseded by §6t7 "ASYMMETRY THREAD CLOSED" (full: archive/2026-11/MEMORY_DIGEST.md:389):
- §6t3: advisor's Jaccard sweep F 0.6044 vs shipped 0.7016 (86%); r911 "drop MCL keep conjunct" FAILED (0.5869) — MCL makes the conjunct affordable.
- §6t4: DP chaining REFUTED, redundant by construction (34/40,645 = 0.08% multi-anchor); r913 `min(la,lb)` makes hubs.
- §6t5: r914 Ochiai = Jaccard; r915 asymmetry blind spot 25% of true pairs (59/238); r916 guarded containment -> MCL F .6648 vs .6394.
- §6t6 r917: guarded-containment weights in full rule F 0.7256 vs 0.7123, but two families only, post hoc in `mcl_port.py` — NOT a basis to change the definition.

## Compacted 2026-09-22

Node definition + edge construction, register rows 969-1031, report `docs/archive/2026-09/NODE_DEFINITION_2026-09-22.md` (full: archive/2026-11/MEMORY_DIGEST.md:405):
- §6x9 containment is a strict partial order (AUC 0.943); r1021 partition-vs-poset limitation; §6y0 r1022 held-out down-sets never beat their NULL; r1023 margin changed sign (always run a null arm); §6y1 r1024 intron-length trie 99.3% redundant.
- §6y2 r1026: 81.1% of chr16 referee same-family pairs (1,195/1,474) have NO alignment; scope not defeat (NPIP 90.6%, FAM90A 100%).
- §6y3-§6y5 r1027-r1031: protein sensitive mode NO; full protein homology CLOSED (FM 0.79-0.98 at every floor); edges right, truth missing (AKR1C1~AKR1C4) — needs Compara/gene trees.
- §6x1 union truth: guided best, lead widens (F .833/.762/.721); r991 semi-guided "prec 1.000" was universe intersection. §6x2 r994-r999 SEDEF NOT load-bearing (lossy copy of our PAF).
- §6x3-§6x5 `--min-cov-shorter` flag (then default OFF); r1009 harmful in semi-guided; r1011 FP side unsolved (3.0%); r1013 §6w2 premise WRONG (1.02x shrink); r1014 use C=0.70.
- §6x6/§6x7 node FP bounded (228 fused loci, 9.2%); r1017 oracle perfect split +0.021F, NPIP worse — do not build another trim/split/cut. §6x8 r1019: `--min-cov-shorter` is an EDGE-WEIGHT change, admission worth 0.000.
- §6w1 r969 colinearity guard NO (prec .247, AUC .490). §6w2 r970-972 binding gate `cov_longer` (99.8% of 1,026 evictions) — premise later REFUTED by r1013. §6u7 de novo defect is OVER-merge. §6w3 r973-975 5' boundary = dispersion not bias (carry 5' +-150-250 bp, 3' +-5-25 bp).

## Compacted 2026-09-23 — full text of index lines moved out of MEMORY.md

Assembler and bakeoff, 09-22/23 (full: archive/2026-11/MEMORY_DIGEST.md:444):
- §6za r1074-r1076 ASSEMBLER DEFAULTS CHANGED (user): `--assemble-only` = `--assembly-junctions strict` + `--polish-retained-ratio 10` (gorilla chain precision 33.1 -> 35.6); old = `majority` + ratio 0. r1077-r1085 every chain-level rule refuted; chain-rule search CLOSED.
- §6z9 r1067-r1073 single-read chains / non-canonical junctions: no learnable admission rule; FLOOR STAYS AT 2.
- §6z4 r1049-r1055 node formation: do NOT widen reps (fuller evicts, smaller makes hubs); FN remedy = edge escape §6x4.
- §6z2/§6z3 script reduction (r1045 bipartite precision not tie-invariant; r1047 `mcl_port` bin bit-identical).
- §6y6 r1032-r1034 `RUSTLE_ER_COVERAGE_LONGER_FLOOR` unreachable from `family_define`; `bench/family_rna_refine.tsv` parity claim FALSE; check contig names before choosing a truth.
- §6z7 r1059-r1061 seed read pool: GOOD secondaries (within 0.98 of genome-wide best AS) better held-out (needs `RUSTLE_GTF_SECONDARY_AS_TABLE`); coin-toss mis-seeding negligible.
- §6z5/§6z6 tool bakeoff: gorilla 10/15, human 12/15 cells; r1058 `copy_assign` dedupe defect (`--keep-coordinate-duplicates`, default OFF).

## Compacted 2026-09-25

(full: archive/2026-11/MEMORY_DIGEST.md:461)
- §6ze r1090/r1091 `o3_rna_flag` SHIPPED: KB3781 53 candidates, 1 haplotype-confirmed; r1094 structural detector ON BY DEFAULT: 64 candidates, 20 haplotype-confirmed.
- All structures closed (r1022-r1031, see 09-22).
- §6zc/§6zd r1088 quant parity CLOSE (Spearman .906/.786 gorilla); r1089 advisor's reference-bias claim REFUTED (extra copies 100% absorbed).
- §6z8 r1062-r1066 assembly-arm rule search CLOSED; dedupe key stays default.

## Compacted 2026-09-28

(full: archive/2026-11/MEMORY_DIGEST.md:468) §6zb r1086/r1087 streaming `--assemble-only --genome-wide` SHIPPED: whole human genome 482 s / 1.93 GB. r1013/r1017 repeated: over-merge does not cause evictions; do not build another trim/split/cut rule.

## Compacted 2026-09-29

(full: archive/2026-11/MEMORY_DIGEST.md:475) §6zh-§6zj r1096-r1099/r1115 identity spectrum vs Compara (union-coverage rule FAIL; chr16 catalog 85% Alu stubs; `RUSTLE_ER_COVERAGE_LONGER_FLOOR=0.30` stays opt-in). §6zg r1095 minimap2 does not join exons across tandem copies; interleaved dups invisible to `gw_family_catalog`. §6zf/§6zm r1092/r1093/r1103 `--union-certificate` SHIPPED opt-in (foreign rows 844 -> 0) but ties all 1,259 chr16 tied reads. §6za repeated.

## Compacted 2026-10-03

Lines removed/shortened from MEMORY.md 2026-10-03; contains the sub-sections below (full: archive/2026-11/MEMORY_DIGEST.md:482).

## Current work — 2026-09-22

09-22..10-02 working index (full: archive/2026-11/MEMORY_DIGEST.md:486):
- 10-02 famsim simulator (`bench/famsim`); MCL size-2 fragility (later fixed as `--pair-rule`, see 10-25 section).
- 10-01 RNA haplotype count (IsoCon chain) parked; Rust `o3_candidates` SHIPPED OPT-IN 10-03, A12 FAILS 34.7% vs 80% bar; rows 1194-1215.
- 09-30 unit cover held-out Jaccard .660; fusions = duplicon blocks (held-out UNDERPOWERED); Soto CN cut follows duplicons (p 1e-4); container v2 PAUSED.
- 10-01 CAT re-run: SQANTI3 ranking holds 27/30; copy recovery P3 FAILS on CAT NPIP. Soto meeting branch merged to main 10-01 (f38064fd).
- 09-29 copy recovery tools rows 1189-1193; cover growth REFUTED; O2 DNA k-mer presence gate REFUTED on HG002; O2 power (m=1 never assigns); gorilla NPIP loss no lever.
- 09-28 r1145-1147 F1 bridge loci; F1v2 DEFAULT ON 09-29. Soto ARI .7096 literal, .97 with their code.
- 09-25..27: readthrough ends filter r1117-r1120 opt-in; read-proven ends r1131-1135 opt-in; rep study r1121-r1127 PARKED; WGS trio CN 17/516; publication figures (register 1109-1115).
- 09-22 closed threads: `mcl_families --min-cov-shorter 0.70` SHIPPED (f2144faf), DEFAULT ON since 09-29 with F1v2 `--bridge-regroup f1v2`; r1026 alignment bound; SEDEF not load-bearing.
- Polish 28/30 cells; SQANTI3 FSM 46.6% / PASS 94.3%; resume pointer `docs/archive/2026-10/merged/PENDING_2026-09-23.md`; cleanup wave 7 (tag `notebook-2026-09-24`); §6s9 O1 emits a STRICT PARTITION; §6t1 divergence cliff at 0.60; §6t7 asymmetry closed; §6t8 VGs don't solve containment; §6t9 nothing beats MCL; §6u0-§6u4; §6u5/u6 Eichler AS-margin comparator (`--eichler-margin`); never `--region-threads`.

## O1 — family definition

(full: archive/2026-11/MEMORY_DIGEST.md:550) 10-03 hierarchy answer `docs/archive/2026-10/merged/ADVISOR_HIERARCHY_2026-10-03.md`; advisor `-G 50k` REFUTED; advisor Jaccard/anchor edge rule not a replacement; §6u8/§6u9 cover prediction and J_N REFUTED; definition then `docs/seeded_family_definition.md` §0; RNA level L3 `w_98 >= 0.985`, L4 0.995 selective; DNA best `--min-exonic-bp 1 --min-shared-exon-frac 0.60`; never quote F 0.88 alone (all 72 = 0.628); ceiling ~0.8; Soto between our L2 and L3; §6p3 self-tuning fails; Louvain not at DNA level (r851); LCS edge core (r814); §6o8 grouping SATURATED.

## Assembler / node construction

(full: archive/2026-11/MEMORY_DIGEST.md:569) §6u7 over-merge; §6w3 5' dispersion; one rep per locus; §6n0 sim ceiling 25/26; §6n2 `RUSTLE_GTF_SECONDARY` (O1 only); §6m8 JUNCTION MAJORITY recorded as default (`denovo_assemble.rs:1727`, verified 09-22; but see 09-23 §6za strict for `--assemble-only`); §6m7 canonical ceiling; §6p4 CDS complete, 5'UTR not.

## O2 / O3 / readthroughs

(full: archive/2026-11/MEMORY_DIGEST.md:581) O2 AS-tied gate (§6hz 0.846); PKD1P6-NPIPP1 real fusion; chimera policy r845 NOT ADOPTED; node cut r846 harmful; readthrough secondary filter opt-in (r852); FAM90A parked.

## Prior art & people

(full: archive/2026-11/MEMORY_DIGEST.md:590) Soto 2025 + full replication ARI 0.6959; 35.9% all-pseudogene; TBC1D3 old truth wrong (16 claims retracted); advisor/user/vault references.

## Older findings

(full: archive/2026-11/MEMORY_DIGEST.md:596) Many `project_*.md` not indexed; 34 older links dead (09-13 link check).

## Compacted 2026-10-13

(full: archive/2026-11/MEMORY_DIGEST.md:601) INVARIANTS corrected: CHM13 v1.0 and v2.0 differ on every chromosome (old constant-shift rule wrong); human annotation = CAT/Liftoff v2.0 since 10-01. 10-04 `o3_candidates` OPT-IN (b29afa55): A13 PASSES (45.7%), A14 control FAILS (35/53), Amendment 15/15b NOT run, rows 1194-1231. Closed-thread and assembler index lines moved to topic files.

## Compacted 2026-10-14

(full: archive/2026-11/MEMORY_DIGEST.md:610) 10-04 STRICT FOUND (Amendment E): NPIP loss = most-read REP; R_J OPT-IN. 10-03 PSV ceiling: 639 identical + 670 beyond-reach gorilla copies.

## Compacted 2026-10-17

10-04..10-09 index (full: archive/2026-11/MEMORY_DIGEST.md:618):
- 10-09 unmapped rescue: recall 15.3%, 82.5% wrong flags, O3 NOT claimable; species screen not adopted. `--emit-hierarchy` (fe5be654); PAF duplicon units NOT ADOPTED (r1263).
- 10-08 O3 maternal reference: missing-copy reads ABSORBED 92-97%, never unmapped; expression ladder pilot NOT run.
- 10-07 docs consolidated to `docs/archive/2026-09|10/` (958972e7); LRPAP1 browser; CN integration (edges cap ARI .83); Yoo tables local (LRPAP1 11/11); seed pool r1261-62; verifiers get copies not products; terminology settled (Soto/Yoo sets ARE multi-copy gene families); entangled baseline r1257-1260.
- 10-06 NPIP/TBC1D3/YAG audit (r1250); few-copy ideal cases r1255; ideal expression through default: TBC1D3 YES 12/13, NPIP NO 12/14 (r1254).
- 10-04 handoff `docs/archive/2026-10/merged/HANDOFF_2026-10-04.md`; LRPAP1 = 11 gorilla loci (rows 1245-1246); week goal = prove O1 and O2.
- Repeats of 09-25..10-02 lines (see Current work — 2026-09-22).

## Compacted 2026-10-25 (MEMORY.md "Current work" lines, verbatim)

10-11..10-14, most recent block (full: archive/2026-11/MEMORY_DIGEST.md:669):
- 10-14 Soto subset x CN: NO subset claim (T1b fails); main f93e9fa9.
- 10-14 Liftoff decoys PASS but weak (75/91, 12/12, 4/4); F1 A119b FAIL 131/148; G1 human FAIL .887 (r1287); post-MCL reassign NOT ADOPTED (r1288).
- 10-14 template copies NOT ADOPTED (r1285); de novo copies NOT CLAIMED (r1286).
- 10-14 O1 layers 1-2 in main; V1/V1b FAIL (.40), V2 FAIL (.52). 10-13 closing O1 plan 1: hierarchy claim FAILED prereg (E1b, E2); `--pair-rule` (exact-mcl = shipped default, user).
- 10-12 fast genome-wide pipeline PAUSED (`docs/GENOMEWIDE_FAST_2026-10-12.md`); YAG O2 split-half .544 FAIL; phased toy (cf4e8c0a).
- 10-11 assembler gap: `--direct-transfrags-union` (r1266, opt-in) beats StringTie in human; union+k2 PASSES (r1268). Rescoring PAUSED until figures agreed with advisor.
- 10-13 popfix: YAGs 0 portable (r1273). MCL two-node fragility fixed as `--pair-rule`; user chose default exact-mcl.
- After this block: rows 1469-1484 (see Current state).
