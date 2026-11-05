# Thesis objectives — the current three

Compact version (2026-11-01). The verbatim full text is `docs/archive/2026-11/THESIS_OBJECTIVES.md`; every `## ` heading below
keeps its exact text and ends with a pointer to the archive line. Thesis model: the **variation graph** is the proposal's family
model; linear alignment is the implementation, and any member-discovery step (minimap2 all-vs-all, sketch) feeds the graph,
it does not replace it.

## Current state (2026-11-03)

- **O1 default (unchanged):** `rustle_pipeline.sh families` = `mcl_families --from-gtf --min-exonic-bp 1 --min-shared-exon-frac 0.60`, minimap2 `-x asm20` all-vs-all of locus spans, MCL 2.8, pair rule `exact-mcl`, representative = most-reads, built from the unlayered `PREFIX.gtf`. Assembly defaults since 11-01/02: `--polish-tss rescue` (r1485) and the sub-chain drop in `PREFIX.final.gtf` (r1498); the end layers E3-E5 passed alone but failed the families non-regression (r1490/1491) and stay opt-in. Every family rule of 2026-10-29..11-02 FAILED or stays opt-in (r1469-1476, 1492-1497; fused-loci annotation-guided cut PASS as recommended opt-in, r1495).
- **O2 scope restated (user, 2026-11-03):** only reads minimap2 could not place: runner-up within 10 on `ms` (any AS) minus the K = 0 wall; exact ties on both scores are the aligner's coin toss and nothing decides them; the AS-gap class is out. Default binary still judges exact AS ties; `--tie-ms 10` (guard on by r1542, 2026-11-05) is the measured configuration: panel 284/2 vs primary 280/6 (r1509, slack-0 truth), held-out 129/22 vs 125/26 after the guard (r1514); on 151 decidable held-out molecules O2 beats the primary 5 to 1. Without the guard 77% of held-out assignments were reads from uncatalogued parent genes (r1510 FAIL). Prohibition stands: defend O2 on abstention and correctness, never "better than minimap2" at scale.
- **Secondaries:** informative down to ms ratio 0.5 (r1511); a gene graph from multimappers is not a family-discovery method (r1512 FAIL); the per-locus pool is (r1513: 6 vs 3 reference-absent contigs, opt-in).
- **Assembly side:** FAMASM FAIL (r1479); IsoCon WORSE both classes (r1478); SQANTI single-copy: StringTie, FLAIR, IsoCon worse (r1482); SQANTI3 filter and own end-artifact rules lose too many exact matches (r1499-1502).
- **O3:** detection only; `candidate_copies` opt-in; informative-secondary pool better on quality, no speed gain (r1513); non-T2T genomes not served (r1504); exon-level reach open (r1416-1420).
- **Spent samples (held-out use):** KB3781, OR6737, A119b, HG002. Experiments end ~early January 2027, then GitHub and species-agnostic scripts.

## Objectives table (2026-09-25 / 2026-10-19 restatements)

| # | objective |
|---|---|
| **O1** | One default de novo family definition at the RNA level (restated 2026-09-25, user decision): IsoSeq reads -> loci from primary + secondary alignments scoring >= 98% of the best (`rustle_pipeline.sh assemble`) -> one representative per locus (positional exon sum of the most-read transcript) -> `families` stage (rule above; edge needs >= 300 bp at >= 70% identity, >= 30% of longer locus exonic length, one exon-to-exon record over >= 60% of the smaller locus; families of >= 2 loci). Pair rule `exact-mcl` default (2026-10-13), `weight-blind` opt-in. O2 consumes the same families via `<id>.fam.copies.tsv`. Judged against external references only: Ensembl Compara, Soto 2025 (not independent), Liftoff. Closed loop (tied-read loop) is a pass-2 experiment, not the definition. Superseded 2026-09-05 SD-core definition kept in the archive for provenance. |
| **O2** | Only AS-tied multimapping reads (rival within 10 AS, Eichler, user 2026-10-19): assign to a copy on PSVs, junctions, divergence, or abstain. Contested set = AS near-ties, not MAPQ-0. No 1/k. 2026-10-10: same certificate assigns reads to haplotypes (ASJ as sub-result). |
| **O3** | Detect + flag expressed paths not explained by reference copies, stratified by whether orphaned reads have anywhere to go; FPR per stratum; completeness never claimed. Amendment 48 (2026-10-09): success = detection; placement is DNA's job. 2026-10-15 advisor: O3 extends O2 (PSV clusters, `copy_assign --families --no-as-tied-only --psv-clusters`, opt-in). 2026-10-19: O3 feeds O2 (`bench/o3_graph/samplecopy.py`); phased data never an input (register 1333). |

ASJ is DROPPED as an objective (cut for time), de-scoped not retracted; back as an O2 sub-result 2026-10-10. "O3" has meant
three things in this repo's history (EM, ASJ, reference-absent copies): resolve the scheme before quoting a number.
(full: archive/2026-11/THESIS_OBJECTIVES.md:1)

## ⭐ STATUS 2026-10-11 — where each objective stands (read before the 2026-08-23 scoreboard below)

Each objective attained in its restated, bounded form; none unconditionally. Genome-wide rescoring PAUSED 2026-10-11.
- **O1:** default families A119b 2,150 / 9,907 loci, testis 342 / 1,115, chimp 396 / 1,343, gorilla KB3781 326 / 918
  (unchanged under `exact-mcl`; `weight-blind` would give A119b 2,218 / 10,043). Fig. 8: 131 of 148 Liftoff copy pairs share a
  family (88.5%). Hierarchy: 215 of 400 Soto families exactly one expansion. Not claimed: hierarchy-aware Soto eval (r1277, E2
  failed 0.488 vs 0.255), Soto as refinement (r1281-1282, T1b failed), de novo Liftoff analogue (r1286), guided G1 (r1287,
  0.887 < 0.90). Liftoff pairs null PASS (r1284: 75/91, 12/12, 4/4 vs decoys 0.121, 0.102, 0.012), bar weak. Not adopted:
  template copies (r1285), strongest-edge reassignment (r1288, precision 0.734 -> 0.688). Limit: bipartite F about 0.63;
  RNA edge recall ceiling 0.052 human / 0.229 gorilla; own duplicon units not adopted (r1263).
- **O2:** haplotype test KB3781 discordance 0.0010 dev / 0.0014 held-out (bar 0.005), 636 ASJs in 286 genes; excision
  abstention AUC 0.80, MAPQ 0.49; hard-locus bakeoff 0.858 vs IsoSeq 0.732, FLAIR 0.473, StringTie 0.441, reverses on the
  contested stratum (0.692 vs 0.813, r1112). Limit: real check is on haplotypes, not paralogs.
- **O3:** erased copies 80.1% reads rescued, 0.07% wrong joins; assembly-truth recall 0.566 paternal (56 / 99), 0.646
  maternal; 2,559 candidates exported (`bench/unmapped_rescue/o3_export.py`); structure projection 210/222, boundary 0.995;
  PSV clusters 425/426 simulated, 0/152 intact flagged; real reads not claimed (r1294, 4 of 12 held out, bar 0.70).
  About 4 in 5 KB3781 flags are other apes' transcripts.
- **O3 reach 2026-10-25 (r1416-1420):** CN differs in 44.8% of 5 kb SD windows between the reference's haplotypes (upper bound;
  WGS-confirmed core 3.0%, exon differences 6.3% of families); breakpoint groups (1417) and certified tier (1418) fail
  held-out; 41.9% of exon gains leave >= 2 reads only the extra sequence explains (1420, control 14%). Open, not claimed.
- **O3 RESTATED (2026-08-19):** original target not achievable. Excision: ORPHANED 33.3% (92.7% reads unmapped), ABSORBED
  64.2% (1.75x depth). Strata: unique/unmapped works (M <= 6.4); paralogous/unmapped vacuous (π = 0.0286); paralogous/depth S2
  partial (TPR 0.2703 / FPR 0.0200; 45.78% of positives below 0.01 divergence). STON1+GTF2A1L n = 1, p = 0.055, UNCONFIRMED.
  Do not quote the 8/9/9-vs-5/6/8 haplotype deficit (does not reproduce).
(full: archive/2026-11/THESIS_OBJECTIVES.md:22)

# ⭐ SCOREBOARD — what is WON, what is DEAD, what is OPEN

Last updated 2026-08-23; superseded where it conflicts by the STATUS block above. Most O1 rows describe the legacy
`gw_family_catalog` (E_r, γ), not the default. (full: archive/2026-11/THESIS_OBJECTIVES.md:98)

## INDEX

Map of the 11 sections below; titles carry the verdicts, no derived tags (an auto-derived tag scored 11/22 = 50%).
(full: archive/2026-11/THESIS_OBJECTIVES.md:106)

## ⚠⚠ THE META-LEDGER — the failure modes that keep recurring

| # | failure mode | scar |
|---|---|---|
| 1 | Denominator conditioned on the prediction | killed 7 findings, 4 headlines in a week (0.8984 was base rate 0.9086) |
| 2 | A control can be broken too | gate control discarded 66% of evidence |
| 3 | Match null on SIZE distribution / COMPOSITION, not count | hapcnv floor 0.1512 vs signal 0.0278 |
| 4 | Judge node/edge changes on the emitted PARTITION | failed 3x end-to-end |
| 5 | T8: offline re-derivation is a hypothesis generator, never a test | 4 offline-proxy errors |
| 6 | argmax is not counting; record `-p` and `-N` | MAPKBP1 1/1 at -p 0.8, 9/8 at -p 0.1 |
| 7 | Never pool human and gorilla; state the unit | |
| 12 | A read that splices over a locus is not evidence for it; use `bench/lib.py::reads_with_block_in` | NPIPP1 inflated 3.4x (1,151 = 71.6% spliced across) |
| 13 | Ask what gene spans an oversized locus (`bench/lib.py::spanning_genes`) | SNX29 100.0%, PDXDC1 99.9% |
| 14 | Exemplars are not a measurement | 10-200x became 2.2x |
| 15 | Validated only where developed is not validated; name the held-out substrate | NPIP 26/26 -> 24/26 on human chr16+18 (r760) |
| 9 | Pre-declaring a confound protects only against that one | length vs exon structure |
| 10 | A parameter change makes an old number look unreproducible | reproduces at `RUSTLE_GATE_MIN_READS=3` |
| 11 | "Pure restriction" true of one layer, false of another | 29 copies created |
| 8 | Substrate provenance | shipped catalog from OR6737 testis, not KB3781 |
(full: archive/2026-11/THESIS_OBJECTIVES.md:126)

## O1 — define a multi-copy gene family topologically at the RNA level

Read first: default = `families` stage; most results here measured on the legacy E_r / γ catalog or opt-in definitions.
- **WON:** cross-species prereg test PASSED 5/5 (coverage clause; refused a default, §6bv); self-overlap defect 0/4,176 on
  the current path (3 contigs only, not genome-wide); shipped catalog 627 families / 2,019 copies (gorilla); E_r relation
  transports across animal and tissue (87.06% edges, ARI 0.9707 without largest component); orientation guard 28 -> 3 spurious
  edges; false-merge 2/150 = 1.33% at floor 3 (3/150 = 2.00% at floor 2; human specificity bound); false-omission 9/162 = 5.6%
  (gorilla, never pool); reach 22/40 = 0.5500; identity clause never binds 0/728; DNA = RNA partition 7/7; per-rep exon arrays.
- **DEAD:** relaxing canonical junctions (2 of 2,847 reps); fixing the coverage hole (no threshold exists); reads as third
  source (AUC 0.1259); R5 repeat gate (breaks P1); `min_shared_gmult` as admission (VETO only); density / Fiedler / min-cut
  (AUC 0.36 / 0.35 / 0.35); shared-exon pooling (8.29%); `RUSTLE_LOCUS_EXON_UNION` (-20 recall points);
  `--single-copy-baseline`; facility location / bipartite matching BANNED; `RUSTLE_SHARED_EXON` refuted §4w;
  `RUSTLE_ER_REPEAT_GATE` refuted §4x (NPIP 12/31 -> 7/31); lowering E_r coverage floor refuted §5c.
- **OPEN:** definition gaps 2026-09-05 (new definition accepted; old opt-in); positive stratum is the binding constraint
  (NPIP labels 21 of 678); seeded discovery 23/23 expressed NPIP vs 13/23 (n = 1 family, a demonstration); footprint nodes
  18/31 vs 12/31 (offline, T8); oracle ablation: definition owns 3%, node construction 58% (30/31 with perfect nodes);
  scale-free coverage denominator (41/494 = 8.30%, superseded catalog); stub defect real, not fixable by rep choice
  (E_r upgrade 0/215; `RUSTLE_SPLICED_REP` regressed chr7 F1 0.570 -> 0.411); orientation guard largest loss channel
  (167 = 0.1471); hierarchy NO-GO (§4m, legacy E_r); strand placeholder resolved, `RUSTLE_READ_STRAND` opt-in, line CLOSED;
  blob 480/2,019 = 0.2377; completeness deficit has no known repair.
(full: archive/2026-11/THESIS_OBJECTIVES.md:150)

## O2 — assign a read to a copy, or abstain

Domain: AS-tied reads only (user 2026-10-19), `fastset --o2-scope tied --as-slack 9`; one genome first (CHM13). Fewest-edit
rule round 23 PASS (r1338): of 7,450 reads tied within 10 AS, 2,045 placed right, 0 wrong (v17 430 / 0); the rest keep a set
(2.02 -> 1.73 copies) holding the true copy.
- **WON:** abstention TPR 0.5066 / FPR 0.0280, AUC 0.7995; MAPQ at chance (AUC 0.4944); contested set = AS near-ties
  (21.75%) not MAPQ-0 (0.0004); τ = ln((1-p)/p), P(misassign) <= 1e-3; per-read argmax optimal given copy set; O1 -> O2 is a
  conditioning (14/104,147 = 0.013%).
- **DEAD:** reassignment as headline (1.96%, never claim "better than minimap2"); EM changes 0/3,081; votes -> likelihood
  identical 16/16; winnowmap identical; λ / density / capacity predictors; structural-anchor truth (29/29 rejected); C5
  stays OFF (p = 0.2685).
- **OPEN:** `n_decisive = 0` floor is correct abstention; one excision panel only; whole-sample yield 2026-10-31 (HG002
  17,805 of 141,265; A119b 4,331 of 212,176) is a catalog-completeness + 1-10 AS near-tie problem. Later: r1481 O2 LEVERS
  arm U FAIL on the panel, arm C open, r1483 `--tie-ms` FAIL stays opt-in, r1484 PSV VCF report.
(full: archive/2026-11/THESIS_OBJECTIVES.md:292)

## O3 — detect + flag reference-absent copies (**completeness never claimed**)

- **WON:** excision positive control (162 two-copy families); two fates ABSORBED 64.2% / ORPHANED 33.3%; signature = unmapped
  reads, not clipping; S2 TPR 0.2703 / FPR 0.0200; sensitivity set by divergence (0.4500 vs 0.0588); DNA instrument FP 0/817;
  nobody has found a reference-absent copy from transcriptome data. Depth caller for ABSORBED (simulation) TPR 0.4248 / FPR
  0.0000, 0.7333 below 0.01 divergence (FPR not a real-data figure).
- **DEAD:** unmapped reads in paralogous stratum (π = 0.0286, vacuous); haplotype route (Jaccard 0.5263, CMH p = 0.5749);
  larger-compartment premise (~1.1x at best); O3 as O1+O2 residual (degenerates at |C|=1); depth residual; DNA-side collapse
  detection (0 at 816/817, underpowered by construction).
- **OPEN:** STON1+GTF2A1L n = 1, p = 0.055; 45.78% of positives below divergence line; n >= 3 excision panel (6/43 = 13.95%
  below de 0.01, wrong regime); larger compartment for the DNA side.
(full: archive/2026-11/THESIS_OBJECTIVES.md:338)

## ⚠ DROPPED — allele-specific junctions (ASJ)

Cut 2026-08-07 for time; the only objective the 2026-06-25 audit rated ATTAINED; de-scoped, not retracted.
(full: archive/2026-11/THESIS_OBJECTIVES.md:374)

---

**Older five-objective scheme (pre-2026-08-07, provenance only).** Claim: assemble isoforms of individual paralogous copies;
the variation graph (`--vg`) shares evidence across copies; an EM assigns ambiguous reads. Status key: ✅ done · 🔄 partial · ⛔ blocked.

## O1 — Formally define a multi-copy gene family in the VG  ✅

Family = connected component of the read-multi-mapping graph (splice-graph generalization) with M (>= 2 copies), H (shared
backbone), X (>= 2 expressed), I (graded identifiability). `classify_family` shipped (`src/rustle/vg.rs`), validated on 14 known
families. Fresh purity challenge 2026-08-16: 124/133 loci, 72/75 named targets re-emitted, 69/75 in the modal family; GOLGA2
typed as broad-family outgroup. Deferred: nested `RECENT_COPY` subfamily; provenance-model / outgroup rooting PoC
(`bench/o1_outgroup_rooting_poc/`, `UNROOTED`). Open: §9 split-range quantification.
(full: archive/2026-11/THESIS_OBJECTIVES.md:388)

## O2 — Recover copies that primary-only assembly misses  🔄 capability shown · ⛔ no validated real-data recovery

Secondary-covered copies StringTie drops; ~93 expressed_real_copy vs ~89 pure_spillover of 846 candidates (classification,
unreproducible from committed state). DAZ3 is a RETRACTED false positive (cov 4.04, abstains). Masked-copy recovery test not done.
(full: archive/2026-11/THESIS_OBJECTIVES.md:444)

## O3 — Assign ambiguous reads to the right copy (the EM)  ✅ synthetic

Fingerprint-EM on ΔNM (T=2 decisive; ΔNM=0 non-identifiable); synthetic read-to-copy 100% / 100%. Real-data validation thin.
(full: archive/2026-11/THESIS_OBJECTIVES.md:467)

## O4 — Assemble each copy's distinct isoforms / structural variants  ✅ synthetic (capability) · 🔄 real-data

Synthetic Sn/Pr 100/100 (baseline also 100/100: non-regression only). Benchmark 2026-06-03: gffcompare cannot discriminate;
copy attribution 100% (id 0.97, 5/5 seeds), abstains at identical copies; overconfident at 1-SNP (dec_acc 0.44 at dec_frac 0.75).
(full: archive/2026-11/THESIS_OBJECTIVES.md:477)

## O5 — Share evidence across copies via the graph  🔄 partial · ⛔ one blocker

Junction channel `lambda_j` inert in the EM. FAMASM (r1479) FAILED held-out: multi-copy recall +1.4 pt on HG002, precision
14.0 -> 8.5; sharing helps assignment, not assembly, at this certificate strength. Intra-bundle splitter blocker was a mirage
(GOLGA6L7 antisense-silent); real gap = paralog bundles rejected / partially assembled.
(full: archive/2026-11/THESIS_OBJECTIVES.md:495)

## Cross-cutting honest gaps

Most scores synthetic; 93/89 rate unreproducible; DAZ3 retracted; GOLGA6L7 splitter unsolved; default de novo headline held
at 95.6 / 90.5. Cited pruned docs are at git tag `notebook-2026-09-23` and `~/Desktop/Rustle_attic/2026-09-23/`.
(full: archive/2026-11/THESIS_OBJECTIVES.md:517)
