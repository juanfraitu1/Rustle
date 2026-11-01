# O1 LEDGER — what has and has not worked (compacted 2026-11-01)

> Compact version. The verbatim original is `archive/2026-11/o1_ledger.md`; every entry ends with a pointer to its line there.
> ⚠ **CATALOG PROVENANCE.** Figures quoted per **/494** (families) or **/1415** (copies) describe the superseded 2026-07-17 catalog. The current default emits **627 families / 2,019 copies**. See [`NUMBERS.md`](NUMBERS.md) and [`o1_catalog_provenance.md`](o1_catalog_provenance.md).

## Current state (2026-11-01)

- **This file is history.** New O1 work is recorded in [`LEDGER_O1.md`](LEDGER_O1.md). `o1_ledger.md` is the O1 history up to its last entry (§6r9, 2026-09-20). Later rows (to 1484) live in `LEDGER_O1.md` and [`NEGATIVE_RESULTS_REGISTER.md`](NEGATIVE_RESULTS_REGISTER.md).
- **Shipped defaults as of the last entry (2026-09-20):**
  - DNA family definition: MCL at `--min-shared-exon-frac 0.60`.
  - Lattice cuts: L3 0.985 plus the selective L4 subfamily split at 0.995 (`lattice_rules/engine.py`).
  - RNA side: the triangle-leader rule.
  - Assembler: `--assemble-only --assembly-polish full --polish-isoform-fraction 0.02 --polish-mono-shadow --polish-mono-quantile 0.82 --polish-ism-ratio 0.7`.
  - Reproduce path: Rust-only. `tools/genome_wide_sweep.sh` uses `--jobs 4`.
- **Opt-in as of 2026-09-20:**
  - Read-isoform widening k=5 (`RUSTLE_SHARED_DEFINITION`) and `copy_assign --read-isoform-k` (default 0).
  - `RUSTLE_JUNCTION_MAJORITY`, `RUSTLE_GTF_SECONDARY`, `RUSTLE_GTF_SECONDARY_AS_RATIO`.
  - `--chimera-policy` and `--max-secondary-frac` (default 1.0).
  - Merged-locus cover layer.
  - `--polish-ism-escape`, `--polish-fraction-min-reads`, `--polish-fuzzy-junction`, `--gtf-tpm`.
- **Not adopted:** strand variants V1-V5b, node cut, Louvain, identity-cut definitions, self-tuning L3 cuts.
- **Open at 2026-09-20:**
  - De novo family purity: recall 1.000, but NPIP is split across 9 clusters.
  - Non-circular RNA-level truth.
  - O2 decidability: YAG catalog 0.3% assigned, 94.0% tied.
  - Human recall gap against isoseq.
  - A family clause that separates merges from splits.
  - A paralog-divergence L3 cut.
- **Since then (register rows 1469-1484, recorded elsewhere):** defaults are unchanged.
  - FAIL or opt-in: r1469 U1-FRESH, r1470 DISPERSION, r1471 SKETCHFAM, r1473 representative rule (FAIL), r1474 identity weights (inert), r1475 chain gates (FAIL held-outs), r1478 ISOCON, r1479 FAMASM.
  - r1476 LOCUS-ABSORB: dev recall gate missed by 0.002, held-out OR6737 PASS by equality; PROPOSED default, not adopted (user decides).
  - Semi-guided families (r1472) are not the default.
  - r1480 protein reconciler: NULL.
  - r1481-1484 are O2 levers, SQANTI single-copy, the tie-ms port and the PSV VCF report.
  - See `LEDGER_O1.md` and the register for verdicts.

## INDEX
Section map of the original (titles carry the verdicts). The entries below keep every `## ` heading in its original order. (full: archive/2026-11/o1_ledger.md:14)

## 1. What SHIPPED and holds
`E_r` = identity >= 0.60 and coverage >= 0.50 of shorter, exon-sum reps: false-merge 2/150 = 1.33%. P1 seed-invariance is a theorem. γ-quasi-clique 0.20 (at γ->0 one component = 530 copies = 26.25%). λ certificates 627/627 (75 `cut_certified`). Orientation guard DEFAULT since 2026-08-20 (spurious edges 28 -> 3). `--refine` rejected (494-vs-627 discrepancy).
(full: archive/2026-11/o1_ledger.md:168)
## 2. What is TRUE and must be disclosed
Named hole: ~30 of 105 bad-family cases definitional (scale-free min-length denominator), exposure ceiling 41/494 = 8.30%. GFPT1xGFPT2 0.5353 < ATP1A1xATP4A 0.5689, no threshold separates. Identity never binds; coverage does 100%. Scope: expressed, well-delineated loci; dispersed families 0/8.
(full: archive/2026-11/o1_ledger.md:179)
## 3. What FAILED — by family, with the number that killed it
3a threshold/denominator/substrate CLOSED (HERC2 splits at c_long >= 0.034, first FP dies at 0.05; genomic span c=0.50: 0/14 FP but 30/150 TP; core denominator F1 0.704 -> 0.632). 3b other statistics fail (FLNC edge TP 130/144; block contrast AUC 0.4336; R5 breaks P1): FP class is the short-node class. 3c graph theory below chance (Fiedler AUC 0.35); false merges are dense hubs.
(full: archive/2026-11/o1_ledger.md:191)
## 4. What WORKS but only as a CERTIFICATE
Internal signals are sensitive but admit FPs; external ones are precise but low reach. Flank homology (portable, genome only) 0/14 FP, 30/150 TP; union of span/SD/flank 0/14 FP, 41/150 = 0.2733 TP, Wilson95 [0.0000, 0.2153]. Architecture = one high-sensitivity relation + certificates.
(full: archive/2026-11/o1_ledger.md:239)
## 4a. ⚠ The search for more FP rules has hit diminishing returns — structurally
Binding constraint is genome-wide TP cost (2.05% / 3.67% / 12.80% / 80%); reference space (nodes, graph, genome, SD catalog, strand, softmask) enumerated. What remains is a different KIND of reference (another individual/species, experiment).
(full: archive/2026-11/o1_ledger.md:274)
## 4b. ⭐ RECALL — what is missing so ALL members are found (2026-08-21)
Against SEDEF: ~1,374 expressed SD locations corrected to ~282, then RETRACTED by §4e (replacement: 275 loci across 120/627 families, 0.8067 node construction / 0.1919 `E_r` / 0.0014 γ). Tiers 2-3: 1,138 pairs, FOUND 810 = 0.7118, coverage 82, γ 246; "γ costs 3.0x coverage" RETRACTED, sign reversed (§4m).
(full: archive/2026-11/o1_ledger.md:296)
## 4c. The 114 questionable γ splits, adjudicated (2026-08-21)
Annotation adjudication: 75/114 uninformative (LOC ids), 11 γ correct, OVER-SPLIT 2/114; among 13 adjudicable 2/13 = 0.1538, Wilson95 [0.0432, 0.4229] (DHRSX). Counting LOC as correct had inflated exoneration ~7x. γ question unsettleable by annotation.
(full: archive/2026-11/o1_ledger.md:402)
## 4d. Can γ be SIMPLIFIED? (2026-08-21)
Edge veto cannot replace γ: dropping gmult >= 5 (63% of edges) still leaves a 206-node blob. Only one γ implementation (Louvain via `refine_families_exon_sum`) reachable from O1; CNM confined to legacy `driver.rs:386`. NP-hard: partition is a certified witness.
(full: archive/2026-11/o1_ledger.md:438)
## 4e. ⭐⭐⭐ The SD "node gap" AUDITED — it was the base rate (2026-08-21)
Never-a-node 0.8984 vs base rate 0.9086 / 0.9092: no gap; "γ costs 1/3,928" etc. RETRACTED. Survives: γ exonerated (107/17,924 = 0.60%); SD truth set cannot measure recall (53.39% below 0.50 floor). Primary-flag finding (45.28%) later RETRACTED by §4i. Defensible: 714 intervals -> 275 loci -> 120/627 = 19.14%; 0.8067 node / 0.1919 `E_r` / 0.0014 γ.
(full: archive/2026-11/o1_ledger.md:489)
## 4f. Design: permissive admission + PSV corroboration (2026-08-21)
Read presence cannot corroborate (phantom rep built from its reads); `E_r` cannot filter phantoms. Proposal: phantom retains ZERO reads after PSV assignment; slot nodes -> `E_r` -> PSV corroboration -> drop refuted -> γ -> λ; outcomes corroborated/refuted/ABSTAINED. Step 1 exon arrays done; three-arm test pre-registered (137 degree-0 intervals under exon-level coverage).
(full: archive/2026-11/o1_ledger.md:591)
## 4g. ⭐⭐ The exon arrays land: 45% of "a node is here" was an INTRON (2026-08-22)
Determinism gate passed. Of 344 intervals with a spliced rep, 45.06% INTRONIC (155) and 189 exonic; stub-only 55/3,928 = 1.40%. γ 1/189 (fourth exoneration). "E_r binding, 188/189 = 0.9947" RETRACTED (base rate 0.9945).
(full: archive/2026-11/o1_ledger.md:657)
## 4h. ⛔ SHARED-EXON / ISOFORM POOLING — CLOSED, NO-GO (2026-08-22)
Pool non-representative isoform exons (T8). Recovers 396/4,778 = 0.0829 edges (dead as replacement); genuine new pairs at `MIN_COUNT` >= 2: 9 vs pre-registered floor >= 100 (dead as addition). "MIN_COUNT >= 2 untested" RETRACTED (tested 2026-08-03: -40.73% copies). Two bugs (missing distinct-locus guard, 93.06% self-overlap, FIXED; MIN_COUNT counts isoform multiplicity, not landed). Gorilla RefSeq GFF exists (`GGO_genomic.gff`). γ non-monotone under edge addition.
(full: archive/2026-11/o1_ledger.md:699)
## 4i. ⛔⛔ §4e's PRIMARY-FLAG FINDING IS REFUTED BY ITS OWN CONTROL (2026-08-22)
First control tie (0.1236 vs 0.1180, p = 0.7701) withdrawn as `de <= 0.05` filter artefact (unfiltered p < 0.0001). Pre-registered volume-insensitive test `bench/o1_competitor_resolved.py`: 233/108,688 = 0.0021 vs 151/37,559 = 0.0040, Mann-Whitney p = 0.2685 TIE -> spillover. VERDICT: C5 stays OFF permanently, §4e REFUTED, O1 ⊥ O2 STANDS; "keep a site if >= 3 reads share an intron chain" REFUTED; §4f premise gone, do not build.
(full: archive/2026-11/o1_ledger.md:791)
## 4j. O1 → O2 IS A CONDITIONING, NOT A CORRELATION (2026-08-22)
O2 consumes only family cardinality (`copy_assign.rs:453`); collapsing all families to n = 2 moves 14/104,147 = 0.013% of verdicts. λ = 1 on 552/627 is mostly entailed (440 are 2-copy); informative 112/187. O2 defended on abstention: TPR 0.5066 / FPR 0.0280, AUC 0.7995 vs MAPQ 0.4944.
(full: archive/2026-11/o1_ledger.md:889)
## 4k. O3 AS A RESIDUAL OF O1+O2 — DEGENERATES ON ITS OWN POSITIVE CONTROL
Abstention residual undefined at 162/162 (two-copy panel); PSV residual = shipped S2; depth residual already refuted. O3 bar: TPR 0.2703 / FPR 0.0200; 45.78% of positives below de 0.01. Levers: n >= 3 excision panel (check divergence first), larger compartment.
(full: archive/2026-11/o1_ledger.md:925)
## 4l. ⭐⭐⭐ CROSS-SUBSTRATE REPLICATION — the RELATION transports, the PARTITION less so (2026-08-23)
OR6737 testis vs KB3781 fibroblast (depth-matched 1,627,629 primaries), real binary. `E_r` edges recovered 1,171/1,345 = 0.8706; ARI 0.7064 (perm p <= 0.005); pair Jaccard 0.5483; fibroblast merges 2.4x more than splits. Partition caveat later shown blob-driven (§4m: ARI 0.9707 without largest component).
(full: archive/2026-11/o1_ledger.md:943)
## 4m. ⛔⛔⛔ γ IS NOT THE PROBLEM — THE ORIENTATION GUARD IS (2026-08-23)
0/114 questionable γ splits share an `E_r` component (113-114/114 minus-strand only). Shipped-instrument loss: orientation guard 167 = 0.1471, coverage 147 = 0.1295, γ 11 = 0.0097 (reverses §4b). Hairball corrected to 480/2,019 = 0.2377. Hierarchy NO-GO (coarse level 78x less supported); MODIFIED-GO: `er_component_id`, `er_component_n_families`, `gamma_forced` columns as disclosure only.
(full: archive/2026-11/o1_ledger.md:1011)
## 4n. ⭐⭐⭐ THE ORIENTATION GUARD IS FILTERING ON AN UNMEASURED FIELD (2026-08-24)
`strand.unwrap_or('+')`: 5,928/5,928 single-exon reps '+' (33.07% unmeasured). Of 4,009 guard-blocked pairs 3,951 = 0.9855 involve a single-exon rep (2.50x vs kept); 58 both spliced. Read orientation recovers strand at 386/400 = 0.9650 vs 0.4867. Fix = measure strand; acceptance on partition, human panel, antisense rate, 58 pairs.
(full: archive/2026-11/o1_ledger.md:1086)
## 4o. ⭐⭐⭐ READ-STRAND MEASUREMENT — RUN, AND EVERY PRE-REGISTERED CRITERION PASSES (2026-08-24)
OFF inert (md5 35eb5c51c141cda25fc7c1866d310f1a). ON: reps 17,924 -> 15,626, families 627 -> 630, copies 2,019 -> 1,958; guard-blocked 4,009 -> 551; gained 972, lost 864, NET +108. All criteria PASS (0/58, 3/630, human 2/150). Node-construction fix; 61 copies left -> ships opt-in `RUSTLE_READ_STRAND`, default unchanged.
(full: archive/2026-11/o1_ledger.md:1156)
## 4p. ⚖️ WHERE THE 61 COPIES WENT — ROUGHLY HALF CORRECTION, HALF LOSS. **DO NOT FLIP THE DEFAULT.**
Churn 224 lost / 163 gained; single-exon losses 144/224 = 0.6429 vs 0.2680 (2.40x). "16 correct / 18 lost" of 35 dissolved RETRACTED IN PART by §4r (exact matching: 32 rep-re-picked / 3 opposite-strand / 0). DECISION: `RUSTLE_READ_STRAND` stays OPT-IN, default UNCHANGED; proposed abstention.
(full: archive/2026-11/o1_ledger.md:1221)
## 4q. ⛔ ABSTENTION DOES NOT RESCUE IT — THE 18 LOSSES ARE MARGIN-INVARIANT. **LINE CLOSED.** (2026-08-25)
Margin 0.90 real-binary run: 624 families / 1,948 copies, losses 18 in both arms (margin-invariant). CLOSED per pre-registered stop; `RUSTLE_READ_STRAND` and `RUSTLE_READ_STRAND_MARGIN` OPT-IN, default OFF. "18 losses" wording RETRACTED IN PART by §4r; conclusion stands (same 32 re-picked families).
(full: archive/2026-11/o1_ledger.md:1288)
## 4r. ⭐⭐⭐ WHY THOSE FAMILIES DISSOLVE — THE REPRESENTATIVE IS RE-PICKED, STUB → SPLICED (2026-08-25)
Exact matching: 32 dissolutions = rep re-picked (same 32 in both arms), 3/5 antisense, 0 same-strand-intact; 39/41 spans still have an ON rep; stub->spliced 27/39 = 0.6923 vs 0.8361 genome-wide (not enriched). Strand fix re-picks 2,306 reps, 1,928 stub->spliced (UNVALIDATED as improvement). Root cause: one-representative-per-locus. Merged sub-block "corroboration layer" (T8 certificates: flank 0/14 FP, 30/150; union 41/150; SD containment 136/627 = 21.7%; vetoes orientation 2.05%, gmult 3.67%).
(full: archive/2026-11/o1_ledger.md:1343)
## 4s. ⛔⛔ REPRESENTATIVE SELECTION IS CLOSED — `E_r` HAS NO UPGRADE PATH (2026-08-25)
Pre-registered kill (< 20% inherited): spliced container inherits stub's edge 0/215 = 0.0000 [upper 0.0170]. Spliced rep regressed at partition (`RUSTLE_SPLICED_REP` chr7 F1 0.570 -> 0.411, row 357). 1.49x projection: 2,891/4,778 = 0.6051 below floor (corrected by §4t). MODIFIED-GO: emit `frac_pure_intron` as disclosure; within-family 22/118 vs 16/118 = 1.37x.
(full: archive/2026-11/o1_ledger.md:1512)
## 4t. ⛔⛔ THE SCOPED ORIENTATION GUARD — SOUND IN PRINCIPLE, CATASTROPHIC IN EFFECT (2026-08-25)
Genuine upgrades blocked by guard 123/124 = 0.9919 (coverage only 1/215 = 0.0047, corrects §4s). Scoped exemption admits 2,224 pairs: 1,912/1,913 join different families, 77/627 touched, 17 blocks, new families 77. Human panel blind (0 qualifying pairs). Killed offline; `RUSTLE_ER_GUARD_SCOPED` ships off; stub defect has no remaining route, disclose.
(full: archive/2026-11/o1_ledger.md:1595)
## 4u. ⛔ WGS CANNOT GIVE A COPY-NUMBER BASELINE FOR O1 — AND THE ONE IT COULD GIVE NEEDS NO WGS (2026-08-25)
Killed at design: within-family identity median 0.8234 -> P(21-mer identical) 0.0169; 1,524/3,715 = 0.4102 share < 1% of 21-mers. CN_WGS = CN_assembly + collapse_deficit; latter is O3's. Not in definition (wrong thing, breaks portability, closed entry). Instead count copies in `GGO.fasta` (record `-p`/`-N`) as annotation.
(full: archive/2026-11/o1_ledger.md:1662)
## §4v — the completeness deficit is REAL, not a repeat artefact (2026-08-26)
3-contig fibroblast NPIP catalog (2,847 reps / 83 families). Edge rate falls with completeness: single-exon 227/1,234 = 0.1840 vs >= 15 exons 16/276 = 0.0580 (scale-free coverage denominator). Stub edges mostly not repeat-driven (median mult 1). "Veto is free precision gain" RETRACTED by §4x. Deficit SURVIVES (3.67x). Claim reproducibility and falsifiability, NOT adequacy.
(full: archive/2026-11/o1_ledger.md:1734)
## §4w — the shared-exon edge rule does NOT fix the completeness deficit (2026-08-26) — REFUTED
Absolute aligned-bp floor (`RUSTLE_SHARED_EXON`): >= 15-exon rate falls in every arm (0.0543 -> 0.0399 -> 0.0326); deficit ratio narrowing is invalid endpoint. OFF -> SE60 families 83 -> 59, NPIP 12/31 -> 10/31. Verdict REFUTED; isoform-pooled variant (2026-08-27) worse: 134 families, pure NPIP 3 -> 0.
(full: archive/2026-11/o1_ledger.md:1788)
## §4x — the repeat-hub veto is REFUTED for the catalog path (2026-08-26)
`RUSTLE_ER_REPEAT_GATE=20`: NPIP loci 12/31 -> 7/31, pure NPIP 3 -> 1; 7/11 lost copies at NPIP. Genome multiplicity cannot distinguish repeat from high-copy family (stub edge rate mult >= 20 -> 1.0000 tautological). REFUTED, flag retained default OFF, DO NOT ENABLE. First ON run was a silent no-op: make a no-op announce itself.
(full: archive/2026-11/o1_ledger.md:1851)
## §4y — reference bias is NOT the binding constraint; DEPTH is (2026-08-26)
NPIP: testis OR6737 3/31, fibroblast 41% subsample 1/31, fibroblast full depth 3 contigs 12/31: 12x from depth alone; premise false. Excision fates: ABSORBED 104/162 = 64.2% (1.75x depth ghost), ORPHANED 54/162 = 33.3% (COPY unit; 34.53% is READ unit). S2 TPR 0.2703 / FPR 0.0200, power 0.0588 below 0.01 divergence.
(full: archive/2026-11/o1_ledger.md:1895)
## §4z — a WGS depth caller for the ABSORBED stratum: sound design, BLOCKED ON DEPTH (2026-08-26)
Blocker 1 (0.65x) WITHDRAWN: truncated download; 102 WGS runs = 456.91x; HiFi SRR26152597 + SRR26152596 = 27.09x. Blocker 2: native compartment pre-registered < 1 collapse, so validate on the excision ablation (162 arms). Verdict REVISED: depth blocker gone.
(full: archive/2026-11/o1_ledger.md:1943)
## §5a — a DEPTH CALLER FOR THE ABSORBED STRATUM WORKS, and it is complementary to S2 (2026-08-26)
Simulation from complete `GGO.fasta`, aligned to masked. Sequence predicts 113/162 = 69.8% absorbable vs observed 64.2%. TPR flat in coverage; 15x knee (TPR 0.3982 / FPR 0.0177); refutes own ~25x requirement. Divergence < 0.01: TPR 0.7333 vs S2 0.0588; overall TPR 0.4248 / FPR 0.0000 / AUC 0.8034 (FPR optimistic, threshold not held out). Real-data arm UNRUN (download cancelled).
(full: archive/2026-11/o1_ledger.md:2027)
## §5b — simulation characterises O1's RECALL exactly, and CANNOT touch its precision (2026-08-26)
210 constructed pairs: hard cliff at 50% shared (0.50 -> 0.3667, 0.40 -> 0.0000), divergence irrelevant. Recall by floor: 0.40 -> 0.6095, 0.30 -> 0.7810, 0.20 -> 0.9286 (shipped 0.4810). Precision unanswerable: three negative constructions, three answers. Do NOT lower the floor on recall alone.
(full: archive/2026-11/o1_ledger.md:2105)
## §5c — the coverage floor CAN be lowered on edge evidence, and it FAILS end-to-end (2026-08-26)
Real negatives (87,990 genomic pairs): FPR 0.026287 at 0.50 vs 0.028094 at 0.20, TPR 0.4810 -> 0.9286. End-to-end: NPIP unchanged 12/31, pure NPIP 3 -> 1, largest 39 -> 104. VERDICT: do NOT lower the floor; target node construction (18 loci, corrected to 4 by §5d).
(full: archive/2026-11/o1_ledger.md:2167)
## §5d — ⚠ CORRECTION to §5c: node construction has 4 loci of headroom, not 18 (2026-08-26)
>= 50%-inside reads: NPIP loci with >= 3 reads 23/31 (not 28/31), reads 419 (not 5,544), both WITHDRAWN. 18 node-less: 7 no reads, 7 pooled chain < 3, 4 genuine downstream failures. Headroom: edge 1, node 4, fragmentation <= 7, untranscribed 7.
(full: archive/2026-11/o1_ledger.md:2229)
## §5e — ORACLE ABLATION: the definition costs ~3%, node construction costs ~58% (2026-08-26)
Oracle nodes + real rule: 30/31 NPIP loci grouped; oracle exon-sums 5/5; shipped 12/31. Loss: definition 1 = 3%, upstream 18 = 58% (7 untranscribed, 7 fragmentation, 4 no node). Open problem is node construction, not the definition.
(full: archive/2026-11/o1_ledger.md:2280)
## §5f — the canonical-junction rule is nearly INERT, and NPIP recall is at the EVIDENCE ceiling (2026-08-26)
`JUNCTION_MAJORITY=1` and `NC_MAX_BP=1e9` move 2 of 2,847 reps, NPIP unchanged 12/31, 0/4 rescued: do not change it. Hypotheses A-D all refuted; the 4 loci have 6-19 reads with no two agreeing on a structure. NPIP 12/31 at the EVIDENCE ceiling; node-construction line CLOSED for NPIP (§5e's 4 revised to ~0).
(full: archive/2026-11/o1_ledger.md:2329)
## §5g — the DNA substrate on RNA-found loci is WORSE, and it explains §5e (2026-08-27)
`--homology-genomic-span`: 91 families, NPIP 10/31 vs 12/31, pure 1 vs 3, 9.6x wall clock. Span of a wrong interval inherits RNA's boundary errors (vs §5e 30/31 on true spans). DNA-first discovery would change O1's claim (scope decision). Do not flip `--homology-genomic-span`.
(full: archive/2026-11/o1_ledger.md:2383)
## §5h — ⭐ REDEFINE THE NODE, NOT THE EDGE: footprint nodes recover +6 loci (2026-08-27)
Footprint node (union of bases covered by >= 2 reads, >= 300 bp, no chain requirement) prototype: NPIP 18/31 vs shipped 12/31; contamination 1/31 = 0.0323. OFFLINE (T8), node-level; superseded by end-to-end §5i (13/31).
(full: archive/2026-11/o1_ledger.md:2433)
## §5i — footprint nodes IMPLEMENTED and measured end-to-end: +1 locus, ships OPT-IN (2026-08-27)
`RUSTLE_FOOTPRINT_NODES` (OFF byte-identical): families 83 -> 87, NPIP in family 12/31 -> 13/31, pure 3 -> 3, 0 lost, 8 fusions. Prototype over-predicted (18 vs 13); `max_gap` re-chosen 5 kb (tuned on scoring substrate). SHIPS OPT-IN, default OFF. +1 later shown to be a lottery ticket (§5o).
(full: archive/2026-11/o1_ledger.md:2487)
## §5j — NPIP is ONE family fragmented into 5–6, not three real subfamilies (2026-08-27)
12 recovered NPIP loci across 5 families, split crosses NPIPA/NPIPB; 3 of 12 loci in two families; 17 NPIP reps alone give 6 components vs true spans 1 component of 30. Fragmentation caused by NODES. One graph per family is CIRCULAR; actionable item is double representation (locus collapse).
(full: archive/2026-11/o1_ledger.md:2532)
## §5k — ⛔ A FAMILY/SUBFAMILY HIERARCHY CANNOT RECONCILE THE FRAGMENTATION (2026-08-27)
γ cannot merge components; NPIP is one component only at floor 0.10, as [1359] with 1,342 non-NPIP reps; at shipped 0.50 the largest NPIP component holds 381 non-NPIP reps. No floor gives an NPIP-specific single component. Do not add a subfamily level; fix the NODES.
(full: archive/2026-11/o1_ledger.md:2572)
## §5l — graph-to-graph similarity: the implemented half is refuted, the other half targets the wrong axis (2026-08-27)
(a) pooled isoforms REFUTED: 134 families, 10/31 NPIP, 0 pure. (b) splice-graph structure unimplemented, can only buy precision, which is not binding (FPR +0.0018 vs TPR +0.4476; edge rule 3%). Best arm stays footprint node 13/31.
(full: archive/2026-11/o1_ledger.md:2621)
## §5m — ⭐⭐⭐ ONE SEED FINDS EVERY EXPRESSED MEMBER (2026-08-27)
Seed `NPIPB11` (asm20 -N 200 -p 0.1): round 1 16/31, converges round 3 at 25/31, precision 0.962, non-NPIP 1 every round. Expressed subset 23/23 = 1.000 vs RNA discovery 13/23 = 0.565; 6 unreached all zero reads. One seed suffices; closure transitive and terminates. n = 1 family (generalised in §5p).
(full: archive/2026-11/o1_ledger.md:2645)
## §5n — seeded windows through the real binary: 26/31, and the fragmentation dissolves (2026-08-27)
`--from-genome` on 26 seeded windows (DNA, no RNA nodes): 26/31 NPIP vs shipped 12/31; NPIP emitted as 2 families (GWFAM0 31 copies, GWFAM1 3). 5 copies outside truth set (purity 0.839, upper bound on FP). 22.6 GB peak. Not an RNA-level result.
(full: archive/2026-11/o1_ledger.md:2719)
## §5o — ⚠⚠ A FOURTH CONSTRUCTION ERROR, SHIPPED BEHIND A FALSE COMMENT; and seeded windows REJECTED (2026-08-27)
Footprints were not exempt from the canonical-motif test: 8/248 gaps canonical (p ≈ 0.24), guard killed 225 of 288 = 78%; FIXED (`build_footprint_seq`): 16/31 nodes, still 13/31 loci, families 100. Seeded windows -> RNA nodes: 12/31, pure 2: DO NOT IMPLEMENT. 22/31 ceiling is an oracle. §5m/§5n STAND.
(full: archive/2026-11/o1_ledger.md:2762)
## §5p — the one-seed closure GENERALISES: 65 human families, 65/65 converge (2026-08-27)
HUMAN CHM13 v2.0, Soto panel (never pool with gorilla NPIP). 65/65 converge, median recall 1.000 every stratum; SD-like 0.885 vs gene-like 0.895 (SD caveat does not manifest); non-member growth 18/65 = 0.277; size r = −0.043. NPIP vs Soto recall not commensurable.
(full: archive/2026-11/o1_ledger.md:2827)
## §5q — EDGE WEIGHTS: the partition discards a 0.9491-AUC signal, and using it changes nothing (2026-08-27)
Edges flattened to 1.0 before Louvain; identity AUC 0.9491 (NPIP-NPIP median 0.9878 vs NPIP-other 0.8281). `RUSTLE_ER_WEIGHTED_PARTITION=identity`: 88 families, NPIP 12/31, 5 NPIP families, pure 3, unchanged. γ test `induced_density` (`family_split.rs:511`) ignores weights. VERDICT: not the lever; flag default OFF. No graph concept promising (partition owns ~3%).
(full: archive/2026-11/o1_ledger.md:2870)
## §5r — ⭐⭐⭐ TIER-2 ADMISSION: expression and assemblability are different questions (2026-08-27)
Funnel: 31 -> 23 expressed -> 13 nodes (−10 expressed but no node) -> 12. Read-cluster span vs node-bearing locus clears `E_r` 10/10; non-NPIP controls 4/200 = 0.0200 (upper bound). Projected 12/31 -> 22/31; STILL A PROJECTION, no end-to-end run; output must carry evidence tier.
(full: archive/2026-11/o1_ledger.md:2924)
## §5s — ⭐⭐⭐ TOP-UP: coverage is NOT the limit — LOCUS COLLAPSE is (2026-08-27)
40 simulated reads per locus added (93.71% landed at source): NPIP in a family 12/31 -> 15/31 (+8 / −5), families 83 -> 76. Of 14 still-nodeless loci 13 pass the read gate and 13/14 have no overlapping rep: transcripts built then ABSORBED by locus collapse (`collapse_loci_span_aware`). Not covered by §4s closure. Spliced arm later invalidated (§5u); finding survives on unspliced arm.
(full: archive/2026-11/o1_ledger.md:2985)
## §5t — ⚠ CORRECTION: the ceiling is NOT expression, and the chain disagreement is an ALIGNER artefact (2026-08-27)
RETRACTED "31/31 unreachable": of 8 unexpressed loci after top-up 6/8 became nodes, 4/8 reached a family; 23-locus ceiling is a property of real-read substrate. Aligner claim (median 28 chains from one template; donor gap median 64 bp) RETRACTED by §5u.
(full: archive/2026-11/o1_ledger.md:3040)
## §5u — ⛔⛔ RETRACTION of §5t's aligner claim: the simulated templates were not transcripts (2026-08-27)
Projected templates were shattered blocks: 6,429 / 6,919 = 0.9292 under 20 bp. Retracted: 28-chain result, donor-gap distribution, §5t's qualification of §5f (§5f STANDS). Survives: valid unspliced arm 10/20 in a family; locus collapse finding; 6/8 top-up result.
(full: archive/2026-11/o1_ledger.md:3091)
## §6a — AMY: the family definition stated on DNA alone (2026-08-27)
SUPERSEDED IN PART BY §6c: `-p` seed-invariance result withdrawn (`-p` inert, `-X` implies `--dual=no`), human comparator wrong (CHM13 9-11 units, not 5), identity/coverage reading withdrawn. Stands: gorilla = 3 AMY copies (seed identities 0.9278 / 0.9303), complete triangle (γ = 1.0), MGAM/MGAM2 excluded by coverage (0.013-0.064). Addendum: `-k 11 -w 5` genome pass OOMs; SEDEF traps (71 Mb mtDNA record; repeat candidate).
(full: archive/2026-11/o1_ledger.md:3151)

## §6b — the DNA definition across four families: precision 1.0000, and orientation is the whole recall gap (2026-08-27)
DNA E_r on 30 CDS-envelope nodes (AMY, MAGEA, RFPL, GOLGA, HERC, SRGAP, NPIP) vs prereg key `amy/PREREGISTER.md`. Forward-only 9/13 = 0.6923 placed vs both strands 12/13 = 0.9231, 0 cross-family; misses were inverted duplications. Size-matched 600-gene negative panel: 0/205 (fwd) and 0/389 (both) false merges, Wilson upper 0.0184 / 0.0098.
SUPERSEDED IN PART BY §6c: "precision 1.0000" unsupportable, scale-free narrowing a zero-candidate result, gain is +23.08 POINTS from a node-extraction (orientation) bug, key wrong in 3 of 8 assertions.
(full: archive/2026-11/o1_ledger.md:3289)
## §6c — ADVERSARIAL AUDIT OF §6a/§6b: seven of eight claims broken or narrowed (2026-08-28)
25-agent audit: 7 of 8 REFUTED or OVERSTATED; counts reproduce, inferences do not. (1) identity clause "inert" (0 of 71,950 records < 0.60; later RETRACTED by §6i, then RESTORED by §6m). (2) scale-free hole real: 36/150 = 0.2400 spurious edges on repeat-rich 1 kb windows vs 0/150 (p = 2.7e-12). (3) recall gain = reference-orientation extraction bug; +23.08 points, CI [0.0, 37.5]. (4) `-p`/`-N` inert at E_r tier (byte-identical PAF). (5) human AMY is 9-11 units, not 5. (6) key wrong in 3 of 8 (GOLGA6L7/6L10 different subfamilies, RFPL4A 8-locus array, NPIPB11 not singleton). (7) correct precision: 0 cross-family among 209 eligible, Wilson upper 0.0180; decision gap empty (0.4086 vs 0.6339) so 0.50 floor UNTESTED.
(full: archive/2026-11/o1_ledger.md:3394)
## §6d — the node set rebuilt from SEQUENCE, not names (2026-08-28)
Name-blind node set (seed discovery, gene orientation, asm20). Claimed 7 clean components, forward-only guard costs nothing (58 vs 58), NPIP = SD blocks, preset changes copy count (MAGEA asm20 3 vs -k11 -w5 ~11-12).
LARGELY SUPERSEDED BY §6f: guard result a 0-of-0, components forced by stage-1 gate, NPIP intervals are 6 unannotated gene copies, hyphen mechanism from truncated index.
(full: archive/2026-11/o1_ledger.md:3539)
## §6e — the sensitivity test: AMY survives, MAGEA is 26 copies, and presets are NON-NESTED (2026-08-28)
Per-contig `-k 11 -w 5 -p 0 -N 500` re-run: AMY 3 copies under both presets (survives); RFPL 587,948 records -> 3 copies; MAGEA 26 copies (asm20 3), symbol grep sees 7/26 = 26.9%; claimed presets non-nested (GOLGA6L copy missed).
PARTLY CORRECTED BY §6f (gap is 3 vs 18, mechanism is scoring not seeding, naming from truncated index); non-nestedness downgraded to UNVERIFIED by §6j.
(full: archive/2026-11/o1_ledger.md:3626)
## §6f — SECOND AUDIT: §6d is largely an artefact of its own construction (2026-08-28)
15-agent audit: four claims REFUTED, one OVERSTATED. Two-stage design circular (58/58 within-seed pairs forced; 0/267 cross-seed). Guard result 0-of-0 (0/202 reverse candidates); reference-orientation relabel drops 58 -> 26, 32/32 rejections true edges, so §6d recommendation WRONG IN DIRECTION (shipped `--from-genome` already disables guard). NPIP: 0.8520-0.9935 of CDS union in 6 of 7 intervals = unannotated gene copies. Naming: symbol-only intervals 0/26 with full census. Preset gap 3 vs 18; gate is SCORING, seeding contributes 0.
(full: archive/2026-11/o1_ledger.md:3690)
## §6g — an honest precision/recall table, and the 0.50 coverage floor TESTED (2026-08-28)
Non-circular panel (annotation nodes, full-GFF product truth, 576 size-matched negatives). At shipped floor: pair recall 223/806 = 0.2767 [0.2469, 0.3086], component recall 48/72 = 0.6667, purity 1.0000, 0/195 cross-family; 3/184 family x negative (likely true NPIP paralogue). Coverage AUC = 0.9808. Recommended 0.50 -> 0.30 (component recall 0.8611).
Floor recommendation RETRACTED by §6h (register §4.1); P/R table stands.
(full: archive/2026-11/o1_ledger.md:3777)
## §6h — ⛔ RETRACT §6g's coverage-floor recommendation: the register already refuted it (2026-08-28)
Register §4.1 (REFUTED / RETRACTED, High redo risk) already had it: cov>=0.30 402 T / 132 F (P 0.753) vs shipped 340 T / 62 F (P 0.846). Adversarial genome-wide negatives (644,905 records): 0.50 -> 0.30 is a 1:1 trade (2 pos vs 2 neg). Shipped 0.50 stands; open problem = an adversarial AND annotation-independent negative set.
(full: archive/2026-11/o1_ledger.md:3863)
## §6i — ⛔ RETRACT §6c's "the identity clause is inert": I measured the WRONG ESTIMATOR (2026-08-28)
Claimed shipped rule uses `nm/bl` (denovo_pipeline.rs:3091), not `1 - de:f`; under nm/bl identity fires on 1.2-6.9% of records (min 0.1478), so E_r is two-clause.
Itself RETRACTED by §6m (wrong function read; §6c restored).
(full: archive/2026-11/o1_ledger.md:3921)
## §6j — the P/R table, corrected: recall depends on TRUTH GRANULARITY, and the protein tier is real (2026-08-28)
(1) nm/bl recompute changes nothing (223 edges, 0.2767, 48/72). (2) Prereg preset-union hypothesis REFUTED: asm20-only edges 0; §6e non-nestedness downgraded to UNVERIFIED. (3) Same edges: coarse pair recall 0.2767 (purity 1.0000) vs subfamily 0.8209 (component 57/62 = 0.9194, purity 0.9286). Protein tier (upper bound, in-frame CDS) coarse component recall 0.6667 -> 0.7361, 19 new edges 0 cross; worth end-to-end A/B with `protein_tail` ON.
(full: archive/2026-11/o1_ledger.md:3970)
## §6k — TBC1D3: subfamily as a LEVEL of one tree, selected by the widest merge gap (2026-08-28)
Gorilla TBC1D3 located by homology: 18 copies (11 full-length at 0.9713-0.9785, 5 partial, 2 distal) vs human 9; identity bimodal (gap 0.9054-0.9713). Average linkage on d = 1 - identity x coverage; widest merge gap 0.0427 selects k = 6, all clusters pure, the 11 full-length copies form exactly one cluster. No single global level transfers across families: read level per family. n = 1 family.
Hedge from §6l: the divergence ladder is `1 - identity` relabelled (partial r -0.0465, p = 0.8919); level selection unaffected.
(full: archive/2026-11/o1_ledger.md:4041)
## §6l — the hierarchy audit: what the tree does and does not add (2026-08-28)
(1) The 52% missing-record gap is cross-subfamily; within-subfamily only 13/201 = 0.0647 missing. (2) Coverage carries the signal: ARI 0.5849 (1 - identity) vs 0.8484 (1 - identity x coverage); edge-rule gating 0.8374. (3) Merge heights add nothing beyond identity (partial -0.0465, p = 0.8919). (4) Lift over single linkage +0.0351 ARI; best_k oracle-selected (upper bound). (5) Labels circular; genomic position is the non-circular validator (ARI 0.6333). No transferable height (intersection EMPTY).
(full: archive/2026-11/o1_ledger.md:4113)
## §6m — ⛔ RETRACT §6i, RESTORE §6c: I read the wrong function (2026-08-28)
§6i withdrawn: line 3091 is in `tier2_rescue`; main E_r (denovo_pipeline.rs:4665-4673) uses `1 - de` (fallback dead: 0 of 984,574 records lack `de:f:`). §6c RESTORED: 0 of 984,574 records below 0.60; identity clause cannot fire on E_r path. Defect exposed: tier2_rescue used stricter `nm/bl` (fixed in §6s).
(full: archive/2026-11/o1_ledger.md:4195)
## §6n — the precision/recall table with error DIAGNOSTICS (2026-08-28)
Main predicate `1 - de`; nm/bl gives identical 223-edge set. Subfamily: pair recall 165/201 = 0.8209, precision 165/223 = 0.7399, component 57/62; coarse: 223/806 = 0.2767, precision 1.0000. Misses: 23 fail coverage, 13 no record, 0 identity. 58 cross-subfamily edges are near-identical paralogues (annotation splits), 0 cross-coarse. 3 family x negative merges = NPIP SD (annotation disagreement). Only 3 unannotated copies, all MAGEA.
(full: archive/2026-11/o1_ledger.md:4249)
## §6o — "no alignment record" is mostly HARMLESS; the coverage floor is what actually breaks families (2026-08-28)
No-record pairs: 12/13 = 0.9231 still in same component; coverage-failed pairs only 12/23 = 0.5217. Of 12 separated pairs 11 fail coverage. Small-subfamily hypothesis refuted (1/12); connectivity matters (MAGE-B misses 23/66, loses nobody; NPIP 3/5, MAGE-D 2/4, MAGE-E 1/2). Only MAGEE1 ~ MAGEE2 (1 of 201, 0.50%) beyond every method.
(full: archive/2026-11/o1_ledger.md:4321)
## §6p — O2 preconditions on the fibroblast substrate: threading is possible, the population is 40 reads (2026-08-28)
PSV threading implemented (`vg_family/copy_assign.rs`; register's `psv_linkage.rs:292` cite stale). TBC1D3 ideal (0.9784 identity) but 9 primary FLNC reads: scope EMPTY for O2. Deep loci (11,563 reads) have 0 near-ties; only NPIPB11 has a target: 40 near-tie reads (0.3509). POWER statement, not refutation.
(full: archive/2026-11/o1_ledger.md:4377)
## §6q — O1↔O2 synergy: NOT a coupled definition, but O1's edge identity as an O2 PLANNER (2026-08-28)
Coupling definitions closed (all attempts zero); do not couple. Per-node max E_r identity predicts O2 near-tie rate: >= 0.95 0.7577 vs < 0.95 0.0000 (11,606 reads), Spearman 0.7173 (n = 20). 0.95 cut post hoc. Proposed reported certificate.
Strong form (zero below cut, "skip") REFUTED by §6r.
(full: archive/2026-11/o1_ledger.md:4442)
## §6r — HELD-OUT test of the §6q planner: the signal holds, the clean zero does NOT (2026-08-29)
Held-out shipped fibroblast catalog (356 families, 1,070 copies): >= 0.95 rate 0.5496 vs < 0.95 0.0867 (144,043 near-ties); Spearman 0.5804, permutation p = 0.0005. Signal replicates; strong form REFUTED; planner may RANK, must not SKIP. Emit identity as per-copy certificate; do NOT re-weight the partition.
(full: archive/2026-11/o1_ledger.md:4503)
## §6s — the three objectives end to end: the chain works, the substrate has no natural instance (2026-08-29)
Excision panel: excised median ratio 1.2744 vs control 1.0000; AUC = 0.8924, TPR 0.7965 at FPR 0.0000 (113 excisions) (reframed as ORACLE statistic by §6t). O2 near-tie tracks similarity not missing copy (depth skew 1.49 vs 1.56). RNA between-copy variance up to 109.96x buries 1.27x signal; O3 must be within-locus or DNA. Fixed: `tier2_rescue` now uses `1 - de` estimator and `ER_TIER_FLAGS`.
(full: archive/2026-11/o1_ledger.md:4546)
## §6t — ⚠ CORRECT §6s: the excision ratio is an ORACLE statistic, not a detection rate (2026-08-29)
Masked/complete ratio needs both references (conditioned on the answer); simulated 40x reads, `-ax map-pb --eqx -N 50 -p 0.1`. Survives as mechanism confirmation (1.2744 vs 1.0000, controls max 1.0015x). Shipped detector fires 0/915 whole-genome; `absent_copy.rs` within-locus design needs no change.
(full: archive/2026-11/o1_ledger.md:4614)
## §6u — ⭐⭐⭐ DIRECT PROOF of the advisor's phenomenon, inside ONE animal (2026-08-29)
KB3781 mat vs pat haplotypes, `-x asm20 -p 0 -N 500`. X-linked probes discarded (male, sex confound). 6/6 autosomal controls concordant; GOLGA6L7 MAT 6 / PAT 7: extra full-length unit (6,035 bp, id 0.9673, cov 1.000) in a regular ~40.8 kb tandem array. Primary assembly mosaic 16 PAT / 9 MAT chromosomes, so haplotype-specific copies invisible. n = 1; assembly-resolution alternative not excluded.
(full: archive/2026-11/o1_ledger.md:4666)
## §6v — scaling the haplotype scan: a rate, and why it must be read as an UPPER BOUND (2026-08-29)
101 autosomal probes (<= 20 kb): 8/101 = 0.0792 differ, Wilson [0.0407, 0.1486]. Direction 7 MAT>PAT / 1 PAT>MAT (sign p = 0.0703); copy ratio 1.0606 tracks assembly-size ratio 1.0464. 7.92% is an UPPER BOUND. GOLGA6L7 (PAT>MAT, against bias) remains the credible case.
(full: archive/2026-11/o1_ledger.md:4735)
## §6w — validating the 8: half are fragment artefacts, and the rate falls (2026-08-29)
Full-length-only (cov >= 0.90) test of §6v's 8: 4 ARTEFACT (fragments), 3 unresolved (dispersed), 1 REAL tandem (GWFAM327 24/19, bias-aligned). Rate 8/101 -> 4/101 = 0.0396 [0.0155, 0.0977]; tandem-confirmed 1/101 = 0.0099. Do not quote 7.92%; GOLGA6L7 (§6u) remains load-bearing.
(full: archive/2026-11/o1_ledger.md:4789)
## §6x — tier-2 admission: the mechanism finally fires, and the trade is real (2026-08-29)
Fixed genome-wide `footprint_skeletons` branch (blocked on ANY overlap; now coverage block). Tier-2: NPIP loci 12/31 -> 17/31, NPIP 5 -> 3 families, but purity 3 -> 0 (11/33 impurities distant). NOT shippable as default; ships OFF with `T2~` tid prefix.
(full: archive/2026-11/o1_ledger.md:4836)
## §6y — The remaining lever, closed; §5s RETRACTED (08-30)
`RUSTLE_COLLAPSE_EXONIC` measured: families 83 -> 94, NPIP purity 3/3, loci 12/31 (later REFUTED by §6aa/§6ab). §5s mechanism RETRACTED (locus collapse owns 1 of 15 residual loci); top-up substrate INVALID (18,467 bp injected vs 1.2-2.9 kb real). Residual is node construction: 13 of 15 missing loci have no overlapping node.
(full: archive/2026-11/o1_ledger.md:4879)
## §6z — Stage attribution: WHICH construction step drops the missing NPIP loci (08-30)
Instrumented shipped path: 11 Pass-1 / 4 gate / 2 collapse / 14 reach rep (12 in copies.tsv). Claimed dominant cause expression (median 4 reads) (RETRACTED by §6ad). All 26 read_floor rejections at pooled = 2 vs GATE_MIN_READS = 3. Reproduced exactly on canonical tree after provenance correction.
(full: archive/2026-11/o1_ledger.md:4912)
## §6aa — FLAG-DEFAULT AUDIT: nothing qualifies to flip (08-30)
124 `RUSTLE_*` flags; every measured O1 flag negative or break-even, rest unmeasured; nothing flipped. RETRACTION: §6y "purity HELD" too narrow; COLLAPSE_EXONIC creates 1,991 new co-membership pairs (2,778 -> 4,481), largest family 39 -> 60. Airtight panel blind to locus-collapse flags.
(full: archive/2026-11/o1_ledger.md:4954)
## §6ab — `COLLAPSE_EXONIC` REFUTED as a default: its merges are hub fusion (08-30)
New `bench/merge_precision_arms.py`: direct-E_r-edge backing of new pairs 14.92% vs baseline 49.53%; 1,244 of 1,991 in GWFAM8 (60 copies). REFUTED as default; stays OFF.
(full: archive/2026-11/o1_ledger.md:4989)
## §6ac — O1 NODE FLOOR: 3 → 2, SHIPPED (08-30)
Floor 2: NPIP 12/31 -> 14/31, purity 3 HELD, families 83 -> 121. First "do not ship" verdict was a size confound; size-robust 1.5 pairs per new edge (COLLAPSE_EXONIC 6.5). SHIPPED as `NODE_MIN_READS = 2`; `GATE_MIN_READS` STAYS 3. Pre-08-30 catalog numbers computed at 3. 815 pass / 0 fail.
(full: archive/2026-11/o1_ledger.md:5015)
## §6ad — THE GAP WAS NEVER EXPRESSION: flag-free site construction recovers 26/31 (08-30)
RETRACTS §6z "not expressed" and the "23/31 expressed ceiling": missing loci median 404 secondaries (primary-flag coin toss). `RUSTLE_FLAGFREE_SITES`: NPIP 14/31 -> 26/31, pure families 3 -> 6, 2.0 pairs/edge, denser in every band; recommended ON pending call.
CLOSED as opt-in by §6af (fails negative control).
(full: archive/2026-11/o1_ledger.md:5056)
## §6ae — SHARED-READ EDGES (E_c) AS AN ADDITIVE TIER: no depth threshold exists (08-30)
E_c not circular (54.5%), overlaps E_r only 11.6%, but fusion fraction ~72% at every shared-read floor (72.1% at 2, 71.6% at 100). E_c cannot be an additive definition tier; E_r defines, E_c scopes. Instrument note: first run 0 E_c edges was a `split` artefact.
(full: archive/2026-11/o1_ledger.md:5097)
## §6af — FLAGFREE_SITES CLOSED: recall real, but it fails the negative control (08-30)
Housekeeping negatives (n = 3 on contigs): ATP5F1A FAIL in GWFAM17 n=115 under FLAGFREE. Large-n uniquely-named genes in family: floor 3 12.78%, floor 2 16.63%, flagfree 33.08%. Density did not catch it. CLOSED: STAYS OPT-IN; reopenable as site-level admission bar.
(full: archive/2026-11/o1_ledger.md:5135)
## §6ag — Overlapping multimappers as site evidence: the cheap win (08-30)
Every placement proposes a site; guard: unspliced site needs a primary (`RUSTLE_FLAGFREE_UNSPLICED=1` disables). Guarded: NPIP 26/31, largest family 50, 3/3 controls PASS; residual 68.0% single-exon, large-n 23.29% vs 16.63%. Flag-free stays OPT-IN; needs full-contig run with all 30 housekeeping genes.
(full: archive/2026-11/o1_ledger.md:5172)
## §6ah — The genome-anchored repeat veto is a NULL on top of the site bar (08-30)
`RUSTLE_ER_REPEAT_GATE=50` on flagfree + bar: 92 of 4,507 reps flagged, 188 of 6,554 edges dropped, every endpoint unchanged (26/31, 50, 68.0%, 3/3). NULL; also undercuts the ML label plan (binding constraint is LABELS).
(full: archive/2026-11/o1_ledger.md:5214)
## §6ai — k+1 MEC for reference-absent copies: NO-GO, four independent kills (08-30)
NO-GO: (1) objective degenerate under missing data (GWFAM115 k=6 MEC = 0 at 15.3% scored); (2) clean negative GWFAM83 has largest drop 0.576, no elbow; (3) k = n_copies control fails (2 of 4 below 0.80); (4) substrate expects 0.20 true positives. Two own claims corrected (read sparsity, EM as soft MEC). Concept already ships (289 candidates, 0 admitted). MEC closed for O3 detection.
(full: archive/2026-11/o1_ledger.md:5241)
## §6aj — The advisor's two PSV objections, measured (08-31)
Obj 1 (VCF): PSVs are fixed differences between copy sequences; read gate default-on (36,414 -> 36,410); 0/36,410 columns read-unobserved; RNA-editing confound REFUTED (A->G/G->A ratio 0.453); VCF deliverable REFUTED (`genome_pos` is not a locus; `write_psv_vcf` broken, 0.5228 contradict genome). Obj 2: 0/490 pairs with zero PSVs (min 4); density = 1000 x (1 - identity) tautology; identity-regime claims REFUTED (rho = +0.109). Bug: 517 reads assigned in >= 2 families, 210 disjoint pairs.
(full: archive/2026-11/o1_ledger.md:5297)
## §6ak — Cross-family double assignment: fixed, and INDEPENDENTLY verified (08-31)
`RUSTLE_XFAM_RECONCILE` (off default / report / abstain), commit a8079cf. 519 contested pairs: 309 shared_locus, 99 readthrough_span, 111 cross_family_contradiction. Abstain: 111/111 contradictions demoted, rows invariant 79,175, OFF byte-identical. 824 passed / 0 failed. Second bug fixed (f7ed833, `byte_identity_gate.sh` vacuous). Default stays OFF.
(full: archive/2026-11/o1_ledger.md:5368)
## §6al — VACUOUS INSTRUMENTS: a defect class, one confirmed, scope unassessed (08-31)
Three instruments passed while unable to fail: `byte_identity_gate.sh` (fixed f7ed833), `RUSTLE_DEBUG_LOCUS` hook (fixed fef1317), `family_def_airtight_panel.py` (not fixed); plus a `split` analysis bug. Scope found, not assessed; stale-path claim CORRECTED BY §6am (symlink resolves). Rule: an instrument must be shown capable of FAILING before its PASS is evidence.
(full: archive/2026-11/o1_ledger.md:5415)
## §6am — The vacuous-gate audit: 10 confirmed false PASSes, and two of my premises corrected (08-31)
40 instruments classified statically: 10 FALSE_PASS, 19 FAILS_LOUD, 11 OK, 0 OVERTURNED. Premises corrected: winloci_scratch is a symlink (worse); `target/` absent, 28 unaudited dead-binary scripts. Worst: `soto/gate_sensitivity_sweep.sh` (destructive). §F2 81.7% probably real (control 25/82); scope ~6% of bench.
(full: archive/2026-11/o1_ledger.md:5460)
## §6an — The ten false-PASS gates repaired; every one retains residual holes (09-01)
12 files, 721 insertions, 89 deletions; all 10 can now fail; `gate_sensitivity.tsv` intact (md5 `18156526`). Residual holes remain in each (airtight panel 3 of 7 counterexamples pass on absence; salvage md5s never compared; etc.); 28 dead-binary scripts UNAUDITED.
(full: archive/2026-11/o1_ledger.md:5523)
## §6ao — Provenance of the numbers the false-PASS gates produced (09-01)
`salvage_validation.sh` / `gate_onoff_gw.sh`: no quoted numbers anywhere (o1_ledger.md:5289 was §6al's own sentence). `gate_robustness.md` sweep stands: 467 / 376 / 483 / 513 distinct copy counts; recall band 81.2-82.3% insensitive by construction.
(full: archive/2026-11/o1_ledger.md:5570)
## §6ap — Reads spanning two near-identical copies: keep current for O2, and a LIVE O1 defect found (09-01)
O2 keep current: spanning record interrogates disjoint read bases (1,806 distinct); `quant.tsv` md5 219f9007 across all four arms; population 99 rows = 0.128%, 93 one copy pair. LIVE O1 defect: 98 `SharedAcrossFamilies` overlaps, 71 copies engulfed by another family (34/121 families); `collapse_parent` phase 1 has no chimera test. Same-family spannable pairs invisible to xfam detector (`copy_assign.rs:659`).
(full: archive/2026-11/o1_ledger.md:5596)
## §6aq — Breakpoint concordance RUN: the gap is real splice, and the "copy" is a terminal exon (09-01)
REAL SPLICE: 329/329 bridging gaps canonical, 19/19 junctions GT..AG; 95/95 spanning reads share last acceptor 69,522,690. GWFAM111:3 is the TERMINAL EXON of GWFAM118:0, not a copy: root cause O1. Candidate rule untested at scale (n=1). Degeneracy (§6ax): discriminator invalid below ~3 bridging junctions.
(full: archive/2026-11/o1_ledger.md:5657)
## §6ar — The pre-registered engulfment test: PARTIALLY SUPPORTED, and it bounds the mechanism (09-01)
Prereg: engulfed copies split bimodally. First tally (17.5% / 68.4%) RETRACTED. On 27 informative loci: >= 0.90 10/27 = 37.0%, < 0.50 10/27 = 37.0%, middle 7/27; bimodal (17/27 = 63%), rule NOT refuted but explains ~a third. 23 no-junction loci (tandem mis-alignment reading later DEAD per §6as).
(full: archive/2026-11/o1_ledger.md:5700)
## §6as — Engulfment: the cause is the STRAND PLACEHOLDER, and NEITHER fix ships (09-01)
Cause: single-exon placeholder `'+'` (`strand.unwrap_or('+')`); 58 of 60 engulfments `('+','-')`; `strand_conflict` blocks collapse. A1 `RUSTLE_READ_STRAND=1` (engulfed 20, NPIP 15/31) and A2 `RUSTLE_COLLAPSE_UNSTRANDED=1` (engulfed 1, families 121 -> 87, net E_r -897) both fire prereg kill criteria: NOTHING SHIPS. Housekeeping endpoint dropped from scorer, fires. Split-only certificate vacuous (theorem of code). O2 mis-alignment reading DEAD (MAPQ median 60). 826 passed / 0 failed.
(full: archive/2026-11/o1_ledger.md:5740)
## §6at — Telltales of a minimap2 artifact: 38 enumerated, 12 measured, 3 real (09-01)
3 REAL-TELLTALE, 4 ALREADY-GUARDED, 13 REFUTED. Readthrough guard counts non-canonical noise: 55/157 = 35.03% drops flip under canonical-only map; RefSeq 330x separation; JUNB case. Disposition "filter support map to canonical" (WITHDRAWN by §6au). Guarded: mosaic reads (absorbed by O2), antisense introns, NC giant introns, fragmentary secondaries. Refuted incl. soft-clip asymmetry (1.044), s2/s1 (AUC 0.4840).
(full: archive/2026-11/o1_ledger.md:5819)
## §6au — The advisor is RIGHT on non-canonical junctions; §6at RESTATED (09-01)
NC share collapses with depth (28.34% at 1 read -> 0.20% at 500), but 3 of 7 deep NC sites have no canonical nearby (e.g. CT..AT, 2,407 reads). §6at canonical filter WITHDRAWN; replace with motif-free READ-MASS weighting. All-canonical rule in `build_spliced_seq_with` an unresolved cost; `RUSTLE_JUNCTION_MAJORITY` unmeasured.
(full: archive/2026-11/o1_ledger.md:5894)
## §6av — Cross-substrate recurrence CONFIRMS the non-canonical sites are real (09-01)
Testis OR6737 (different animal and tissue): all 3 sites recur at exact coordinates; CT..AT 87.6% and GT..GG 92.9% of covering reads (CT..CC 1.4%). Advisor CONFIRMED; §6at stays WITHDRAWN; `RUSTLE_JUNCTION_MAJORITY` highest-value untested flag.
(full: archive/2026-11/o1_ledger.md:5943)
## §6aw — `RUSTLE_JUNCTION_MAJORITY` measured: it works, and it is NOT yet a default (09-01)
Prereg K1-K4 ALL PASS: rescues 2 sites (0 -> 1 each), families 121 -> 117, copies 678 -> 700, NPIP 14/31 held, housekeeping 3/3, OFF byte-identical. Not a default: known harm (chr16 copies 66 -> 34) unreachable on these contigs; required chr16 arm. Stays OFF, MEASURED-POSITIVE.
(full: archive/2026-11/o1_ledger.md:5982)
## §6ax — READINESS AUDIT against the three questions the advisor will ask (09-01)
Q1: E_r relation reproduces across animal and tissue (87.06%, ARI 0.7064); but corroboration shares substrate, external evidence only at high identity (median edge 0.8287; 86.31% < 0.90), 5 + 1 families externally adjudicated, SD containment 33/121 = 27.27%, method NPIP-bound, zero experimental validation. Q2: 152 primaries span >= 2 same-family copies; 47.7% undetermined; §6aq discriminator degenerate (39/59 garbage). Q3: ape BAMs drop-in, but all four catalogs built with rejected `refine`; do not quote 149+84 (correct four-way histogram listed).
(full: archive/2026-11/o1_ledger.md:6035)
## §6ay — TSS/TES/UTR boundaries: the premise is FALSE, but the coverage statistic is ONE-SIDED (09-01)
Premise false: within-family median |delta exon-sum| 1,190 bp (ratio 1.43x); clipped sequence in PAF, not lost. Concession: coverage charged on shorter only; 1,035/1,726 = 59.97% of direct within-family pairs have cov_LONGER < 0.50. End difference is 3' UTR (3' flank 1,672 vs 5' 177 bp, 9.45x; placeholder-strand control 0.66x). `RUSTLE_ER_COVERAGE_LONGER` unmeasured.
(full: archive/2026-11/o1_ledger.md:6109)
## §6az — Should rep selection / exon-sum / seeding change? ONE adjustment, and it is DISCLOSURE (09-01)
Rep selection CLOSED (0/215 inheritance; SPLICED_REP chr7 F1 0.570 -> 0.411); exon-sum do not touch (EXON_UNION -20 points); seeding fine for discovery, circular as definition. `RUSTLE_ER_COVERAGE_LONGER` is the wrong fix (gate; up to 60% edges). Ship DISCLOSURE of `cov_longer` and flanks.
(full: archive/2026-11/o1_ledger.md:6176)
## §6ba — `cov_longer` and the flanks are now DISCLOSED, not gated (09-01)
Added `cov_longer`, `unaln_i`, `unaln_j` to E_r edge dump; additivity proven (3,141 edges, 50,256 cells 0 mismatches, projected file byte-identical md5 `5f921a83`); independent Python re-derivation 0 mismatches. 827 passed / 0 failed. cov_longer median 0.3169 overall (75.71% < 0.50); within-family 0.4444, 59.97% below 0.50.
(full: archive/2026-11/o1_ledger.md:6224)
## §6bb — Is a shared-multimapper graph as good as E_r? NO, and a 3-way gauge is worse (09-01)
Advisor's E_c replacement: NO (different, worse relation; loses most E_r edges, cliff at ~0.95, only 143/3,141 = 4.6% edges above 0.95); 3-way gauge worse; thresholds move to aligner, not vanish. Keep E_r defines / E_c scopes; multimapper count optional annotation.
Absolute counts (688 = 21.9%, 13,120) RETRACTED by §6bc; conclusions stand.
(full: archive/2026-11/o1_ledger.md:6262)
## §6bc — The E_c cliff RE-DERIVED: its location is DERIVABLE, and §6bb's counts are RETRACTED (09-01)
Re-derived (`bench/ec_identity_cliff.py`): 564/3,141 = 17.96%, 17,490 E_c-only pairs, top band 83.92%; at >= 50% overlap 0 of 1,203 sub-0.80 edges share a read. Cliff derivable: p >= (tauA + B)/(A + B), A = 1, B = 4 measured, p >= 0.99 at tau 0.95; 65/65 edges >= 0.99 have shared tied read. Artifact `docs/artifacts/ambiguity-cliff.html`.
(full: archive/2026-11/o1_ledger.md:6319)
## §6bd — CROSS-COPY BORROWING IS INERT OR DEAD — and what that costs the isoform claim (09-01)
`consensus.rs` dead (0 call sites); `rescue_thin_loci_iterative` wired but 0 `RC_` copies (0/678, 0/1,070); TOPO_BORROW 0 beneficiaries. Borrowing live only for ASSIGNMENT (O2 PSV matrix). 370/678 = 54.6% copies single-exon. Identical-chain pairs (n=6): above cliff unfalsifiable, below cliff real conservation; rule: report isoform sharing only where read-separable.
(full: archive/2026-11/o1_ledger.md:6380)

## §6be — SEEDED MODE: the seed SUBSTRATE dominates, and a MULTI-seed reaches 30/31 (09-01)
NPIP 31-locus truth, arm_f2. Spliced CDS seed 14/31 (retracted "seed cannot reach"); one genomic DNA seed 25/31; multi-seed (4,477 genes) 30/31, strict superset of reads. Only 10.8% of seeds fire. Weak as family finder: 43/121 = 35.5% of read families; n=74 repeat clique 0/74 corroborated. 0-corroborated loci: 123/305 = 40.3% expressed (CORRECTED §6cn, was 152/305). Seeded = recall instrument, annotation-free path stays the definition.
(full: archive/2026-11/o1_ledger.md:6428)

## §6bf — IS THE SEEDED CATALOG BETTER *AS A CATALOG*? NO — and the test cannot be neutral (09-01)
Coding-coherence (E_r rule on CDS), loci fixed, pair-scored. Genome-wide 118 shared loci: seeded F1 0.8115 vs read 0.8491; with gamma_refine on seeded 0.8257, recall equalised 0.8571, read wins on precision (0.8411 vs 0.7965). Verdict NO: seeded not better as catalog. Caveats: criterion not neutral; covers 4.9% of seeded catalog.
(full: archive/2026-11/o1_ledger.md:6512)

## §6bg — SOURCES OF FALSE MERGES IN BOTH MODES: transitive closure, and one oversized family (09-01)
CDS concordance = relative instrument only (0.8411 shared vs 0.1234 genome-wide; do not quote 0.1234 as FP rate). FP source 1: transitive-only pairs (fams >=11) 21-24% concordant vs 82-85% direct; direct-edge filter F1 0.7536 -> 0.7794. FP source 2 (GWFAM2 oversized) RETRACTED by §6bi (ZNF clusters). T8: offline proxy, no default moves.
(full: archive/2026-11/o1_ledger.md:6594)

## §6bh — CONNECTIVITY AS AN ORTHOGONAL LEVER: NO. And the direct-edge rule does NOT survive being a DEFAULT (09-01)
Connectivity levers (triangle, degree, articulation) all subsumed by direct (none add F1 over 0.7794). k-core as catalog form: seeded k-core 3 NPIP 30 -> 16/31, F1 falls. Pair-filter precision (0.9047) != k-core 2 (0.7348). Direct ships as DISCLOSURE only; cut_certified is the high-precision subset label.
(full: archive/2026-11/o1_ledger.md:6670)

## §6bi — ⭐⭐ E_c SPLITS THE LARGE FAMILIES, AND IT IS THE FIRST ORTHOGONAL LEVER (09-01)
E_c (shared multimapper) as within-family splitter, arm_f2: direct+E_c>=1 precision 0.2406 -> 0.5431; orthogonal to direct. Splits GWFAM79 (54) into 4, concordance 0.1624 -> 0.4612, loses 3/199 true pairs. GWFAM79/GWFAM2 are ZNF clusters -> §6bg FP source 2 RETRACTED; FP source 1 survives (+0.1538 non-ZNF). T8 hypothesis; E_c reaches 13% of pairs.
(full: archive/2026-11/o1_ledger.md:6725)

## §6bj — E_c REFINEMENT RUN AS A CATALOG CHANGE: safe, real, and NOT general (09-01)
E_c refinement as catalog op on arm_f2: guarded rule keeps NPIP 14/31, splits 4/121 families; unguarded refuted (64 families, 12/31). Not general: non-ZNF families median E_c coverage 0.000 (7/12 zero); non-ZNF precision 0.2800 -> 0.2917. Disposition: safe for near-identical arrays, ship guarded/as annotation only.
(full: archive/2026-11/o1_ledger.md:6786)

## §6bk — ⭐⭐ CHARGE COVERAGE ON BOTH SEQUENCES (BLAST qcovs/scovs): the strongest lever measured (09-01)
Longer-side coverage floor (shorter stays 0.50); offline reproduces 3,141 shipped edges. Median cov_longer 0.3175 (arm_f2). Floor 0.30: nZ precision 0.2593 -> 0.4419, F1 0.3182 -> 0.4043, NPIP 14/31 unchanged; floors >=0.40 drop NPIP to 12/31. Shared-domain mechanism is arm_f2-specific (see §6bl). Flag claim CORRECTED by §6bo.
(full: archive/2026-11/o1_ledger.md:6832)

## §6bl — TWO-SIDED COVERAGE ON A SECOND SUBSTRATE: the gain holds, the MECHANISM does not (09-01)
Genome-wide fibroblast (2,075 edges reproduced): floor 0.30 nZ F1 0.4306 -> 0.4595 (families -13%, copies -15%); symmetric 0.50 F1 0.5303 but -44% families. Gain reproduces; shared-domain mechanism does NOT (identity 0.771 vs 0.772). Disposition unchanged, justification weakened; still offline.
(full: archive/2026-11/o1_ledger.md:6909)

## §6bm — TWO-SIDED COVERAGE + CONNECTIVITY: the composite gain does NOT reproduce (09-01)
Global triangle/k-core disqualified (deletes all 2-member families, 255/356 shipped; arm_f2 NPIP 14 -> 13/31). Scoped composite (triangle in fams >=3): gw 0.5066 but arm_f2 0.4043 -> 0.4086 only, -22% copies. Recommendation: two-sided coverage 0.30 alone; connectivity not worth adding.
(full: archive/2026-11/o1_ledger.md:6953)

## §6bn — ⭐⭐ DENSE CLUSTERS INSIDE FAMILIES WORK; CLIQUES DO NOT (09-01)
Re-refine families >=3 at gamma' (scoped density). gamma'=0.70 best: arm_f2 nZ F1 0.4130, gw 0.4646, NPIP 14/31; within-ZNF control split 0.0126 vs kept 0.3866 (31x) arm_f2, 15x gw. Clique (1.00) too strict. Best within-family refinement measured; offline, T8.
(full: archive/2026-11/o1_ledger.md:7008)

## §6bo — NOTHING SHIPPED; and a CORRECTION to the flag claim (09-01)
NOTHING from §6bg-§6bn SHIPPED (only bench binary src/bin/gamma_refine.rs). CORRECTION: RUSTLE_ER_COVERAGE_LONGER replaces the denominator = symmetric 0.50 variant (rejected); asymmetric form needs a new constant. Tree (denovo_pipeline.rs:4900-4926) already argued it (monotone under improving evidence; 134/171 NPIP pairs reach 0.50).
(full: archive/2026-11/o1_ledger.md:7059)

## §6bp — ⭐⭐⭐ THE ASYMMETRIC TWO-SIDED COVERAGE CLAUSE, RUN END TO END: ALL THREE CRITERIA PASS (09-02)
New additive RUSTLE_ER_COVERAGE_LONGER_FLOOR, real binary on 3 contigs. K1 OFF byte-identical to arm_f2; K2 certificate distinguishes; K3 NPIP 14/31 -> 15/31, largest family 54 -> 32; edges 3,141 -> 1,652 (offline predicted 1,643). ALL THREE PASS; still NOT a default (needs human panel, non-ZNF substrate).
(full: archive/2026-11/o1_ledger.md:7096)

## §6bq — ⭐ PRE-REGISTRATION: what the coverage clause must do on the HUMAN negative panel (09-02)
Prereg on human 150-window negative panel (CHM13+A119b, standing 2/150). Score post-guard edges (3), not PAF (28). P1 W106 no family; P2 W063 survives (cov_longer 0.498); P3 2/150 -> 1/150; P4 no window gains; P5 edges 3 -> 1. Falsifiers stated; rate change not significant, mechanism is the evidence.
(full: archive/2026-11/o1_ledger.md:7142)

## §6br — ⛔⛔ THE ANNOTATION CANNOT CORROBORATE THE 0.80–0.90 BAND: a 39.9× enrichment that was ONE STRING (09-02)
Can gorilla RefSeq annotation corroborate the 0.80-0.90 identity band? NO. Two instrument defects fixed (large-gene attractor; proximity null). Band enrichment 39.90x then died: 776/788 agreements = string "zinc finger protein"; non-generic band [0.80,0.90) 0/8. §6ax concession STANDS; annotation cannot adjudicate (does not show edges wrong).
(full: archive/2026-11/o1_ledger.md:7215)

## §6bs — ⭐⭐ THE SELF-OVERLAP DEFECT IS REAL, CROSS-SPECIES, AND ABSENT FROM THE CURRENT PATH (09-02)
Blind-spot item 6 (copies overlapping own sister) MEASURED: refine-built GGO 7.09% [5.14, 9.69], PTR 7.10%; current path arm_f2 0/121 [0, 3.08], cf030 0/99. Substrates not matched (genome-wide vs 3-contig): do not state "fixed genome-wide".
(full: archive/2026-11/o1_ledger.md:7297)

## §6bt — ⭐⭐⭐ THE PANEL RUN: the mechanism transfers to HUMAN, and the BASELINE HAD MOVED (09-02)
Human panel on current binary: OFF 3/150 = 2.00% (W033, W034, W065; disjoint from documented W063/W106); ON 1/150 = 0.67%. 6/6 edge outcomes called by cov_longer < 0.30; P4 passes. 1.33% superseded by 2.00%. 6/6 new FP edges carry a 2-read endpoint (node floor 3->2, §6ac). Earlier "P1-P5 VOID" claim corrected by §6bt.1. Clause stays default OFF.
(full: archive/2026-11/o1_ledger.md:7343)

## §6bt.1 — ⭐⭐⭐ THE CONTROL: the frozen panel IS reproducible, and the node floor is the WHOLE cause (09-02)
RUSTLE_GATE_MIN_READS=3 reproduces documented 2/150 (W063, W106, 3 edges) exactly; cov_longer matches to six decimals. Panel reproducible; §6bq criteria LIVE. Node floor 3->2 is the entire 1.33% -> 2.00% difference (not significant; trade unpriced).
(full: archive/2026-11/o1_ledger.md:7429)

## §6bt.2 — ⭐⭐⭐ EVERY PRE-REGISTERED CRITERION PASSES, ON A SPECIES THE CLAUSE NEVER SAW (09-02)
Floor 3 + COVERAGE_LONGER_FLOOR=0.30 vs prereg §6bq (commit 5cbced4): P1-P5 all pass, 5/5; W106 died, W063 survived (0.497726); 1/150 = 0.67%; 9/9 edges across both floors. Do not report rate halved. Stays default OFF (28% copy cost unadjudicated).
(full: archive/2026-11/o1_ledger.md:7473)

## §6bu — ⚖️ ADJUDICATION OF `RUSTLE_ER_COVERAGE_LONGER_FLOOR`: the decision rule, fixed BEFORE the evidence (09-02)
Prereg decision rule for the clause (cost 678 -> 490 copies, -27.7%). D1 NPIP recall, D2 panel specificity (both met), D3 losses in defective strata (weak), D4 no penalty on corroborated copies, D5 trims not deletes families (required). Map: SHIP / KEEP OPT-IN / REJECT; most likely KEEP OPT-IN.
(full: archive/2026-11/o1_ledger.md:7523)

## §6bv — ⚖️⚖️ THE ADJUDICATION: **KEEP OPT-IN**, and four of the five criteria did not go as expected (09-02)
Verdict KEEP OPT-IN; SUPERSEDED by §6bx (REJECT). Catalog not monotone: 220 copies lost, 32 gained (29 new loci). D4 NULL (exon structure: single-exon 0.013 vs 0.009, multi-exon 0.791 vs 0.784); D5 fails (18 families deleted, 55 copies); D1 gain is one 274 bp 2-read node. Largest family 54 -> 32.
(full: archive/2026-11/o1_ledger.md:7569)

## §6bw — ⚖️ PRE-REGISTRATION: adjudicating the copy loss against SOTO's families (09-02)
Prereg Soto 2025 adjudication (human, 83 families / 362 members, CHM13). SD98 => only adjudicates high-identity stratum. Units: detection, pair P/R (primary), family. Rule: ON loses Soto-labelled copies => REJECT. Expression not limiting: 350/362 = 96.7% have >= 2 primaries, median depth 161.
(full: archive/2026-11/o1_ledger.md:7658)

## §6bx — ⭐⭐⭐ THE SOTO ADJUDICATION: high P/R in the stratum Soto covers, and the clause is REJECTED as a default (09-02)
Default path (arm_off): detection 293/362 = 0.8094 any-overlap; pair precision 0.7198 at >=50% (33.7x null, pooled); band [0.90,1.00) precision 0.9739 (149/153). Clause REJECTED as default: ON loses 7 Soto members (all SD98 = Yes), >=50% precision 0.7198 -> 0.7321 p ≈ 0.77. Supersedes §6bv. Disclosures: intervals unverified (§6cb), first o1-perp-o2 guard firing, -M -L subset BAM.
(full: archive/2026-11/o1_ledger.md:7740)

## §6by — ⭐⭐⭐ A FALSE-MERGE RULE VALIDATED ON EXTERNAL TRUTH: assert co-membership at GRAPH DISTANCE 1 (09-02)
On Soto pairs precision decays with graph distance: d=1 0.6308, d=2 0.2156, d=4 0.0118; holds in every size stratum (d=1 in fams 2-5 0.8168). Identity rule named circular. SHIPPED 2026-09-02 as certificate <out>.pairs.tsv (FamilyCertificate::pair_distance); catalog md5-identical, d=1 count 1,727 = sum n_edges. Gorilla 58.6% transitive. Replicates §6bg on external truth.
(full: archive/2026-11/o1_ledger.md:7830)

## §6bz — ⛔ LOCAL EDGE CONNECTIVITY DOES NOT BEAT THE DIRECT EDGE; and the precise statement O1 can make (09-02)
Local edge connectivity lambda (Menger) on Soto copy pairs (baseline P 0.4637): d=1 0.6998 vs lambda>=2 0.5677; d=1 + lambda>=2 0.7136 at -17% TP. Edge carries signal, not redundancy; size 16+ bad regardless. Statement: family = component, assertion = edge.
(full: archive/2026-11/o1_ledger.md:7950)

## §6ca — ⭐⭐ FAMILIES ARE HETEROGENEOUS, THE CATEGORIES COLLAPSE TO ONE AXIS, AND THE RULE'S VALUE IS CATEGORY-DEPENDENT (09-02)
Soto-metadata categories: outcome bimodal (14/76 families at pair recall 1.00, 15/76 at 0.00); Family MAD collapses to size. d=1 gain held out +0.23 to +0.38 in families <= 9, only +0.08 / +0.12 in 10+ (loses ~40% TP). Same rule, category-dependent trade; large families need a partition fix. Invalid parity split caught.
(full: archive/2026-11/o1_ledger.md:8021)

## §6cb — ⛔ THE SOTO INTERVALS CANNOT BE "CORRECTED" FROM THE GFF — the disagreement is NAMING, not coordinates (09-02)
PARTLY RETRACTED by §6cc (S1A is the right reference: median delta 1 bp). Name-matching to RefSeq: 260/362 = 71.8% matched; large deltas are paralog NAME SWAPS (SPDYE8/SPDYE17). Do not substitute RefSeq coords by name. Clean 64-gene subset underpowered (n=8 at >=50%). Verdict: keep CAT intervals. (§6cd later re-corrects S1A.)
(full: archive/2026-11/o1_ledger.md:8093)

## §6cc — ⭐⭐⭐ NPIP END TO END: all three objectives on the advisor's own family, against Soto's truth (09-02)
NPIP ID_154 on human default path: 13/13 detected, pair precision 0.8478 (first reported 1.0000, retracted); E_r gives 8 components, GWFAM126 density 0.6182 kept as one family (node construction is the limit). O2: MAPQ0 0.0599 vs 0.0004 genome-wide. O3 candidate RETRACTED (§6cd: it is NPIPB15; §6ci: O1 corollary not O3). Adjudicated: 14/14, precision 0.8558.
(full: archive/2026-11/o1_ledger.md:8171)

## §6cd — ⛔⛔ THREE CORRECTIONS TO §6cb/§6cc: THE LIFTOFF BED IS THE MORE ACCURATE SOURCE (09-02)
Three corrections: NPIPB15 BED position confirmed by RefSeq (S1A points at LOC124907807); O3 candidate was NPIPB15 (zero reference-absent candidates); S1A not adopted (RefSeq backs liftoff BED 8-3, s1a BED deleted). Adjudicated 16-member set: 14/14 detected, precision 0.8558 (89 TP / 15 FP). Gorilla 31-locus panel provenance sound but roster != ID_154.
(full: archive/2026-11/o1_ledger.md:8259)

## §6ce — ⚖️ PRE-REGISTRATION: `RUSTLE_ER_NO_STUB_EDGES=1` on the Soto substrate (09-02)
Prereg flag keeps stub loci but refuses their E_r edges (16 of 34 NPIP copies). P1 OFF byte-identical; P2 precision > 0.8193; P3 NPIP families < 10; P4 NPIPB10P the detection loss; P5 recall above delete-model 0.6044. Falsifier: no new family may appear. Fixed 91-pair denominator.
(full: archive/2026-11/o1_ledger.md:8326)

## §6cf — ⚖️ `RUSTLE_ER_NO_STUB_EDGES` RUN: it IS the fragmentation lever, precision 0.9655, and F1 is a wash (09-02)
P2 precision 0.9655, P3 families 3, P4 exactly NPIPB10P lost, P5 0.6154 (met by one pair); falsifier FIRED (3 pairs added, gamma non-monotonicity). F1 worse: NPIP 0.7816 -> 0.7517, held-out 0.4797 -> 0.4733; copies 879 -> 349. A high-precision mode, not a default.
(full: archive/2026-11/o1_ledger.md:8362)

## §6cg — ⭐⭐⭐ SEEDED SEARCH ON NPIP: 5–6 seeds recover 100% of the rest, with ZERO false positives (09-02)
DNA seeded search on chr16 (two disjoint seed sets): 10/10 and 11/11 non-seed ID_154 members recovered, incl. NPIPP1/NPIPB12/NPIPB10P. 10 "off-target" loci all RefSeq NPIP-family => 0 FP; ID_154 incomplete (>= 26 loci). Reachability only, not definition; the 10 lay outside the BAM slice.
(full: archive/2026-11/o1_ledger.md:8415)

## §6ch — ⭐ NPIP ON ITS COMPLETE ROSTER: 26/26 detected — and a VACUOUS precision metric, caught by its control (09-02)
26-locus roster (RefSeq rule UNION ID_154), all expressed; default arm 26/26 detected (no-stub 24/26). P = 1.0000 VACUOUS (degenerate one-family predictor scores 1.0000); precision stands at 0.8193 from multi-family run. Residual = fragmentation: 14 families for one.
(full: archive/2026-11/o1_ledger.md:8466)

## §6ci — ⛔ SCOPE CORRECTION: today's NPIP "novel locus" work is the O1 COROLLARY, not O3 (09-02)
Scope correction: unannotated in-assembly NPIP loci are the O1 corollary, not O3. Replicates 08-07 (8 then, 10 now RefSeq NPIP loci outside ID_154). No novel locus in today's run (all 7 unannotated copies single-exon). 7-step admission standard for an O1-corollary locus.
(full: archive/2026-11/o1_ledger.md:8522)

## §6cj — ⭐⭐⭐ THE FRAGMENTATION HAS A NAMED CAUSE AND AN UNUSED SIGNAL: reads link the pieces, the merge never looks (09-02)
PARTLY RETRACTED by §6cm (SNX29/PDXDC1 are real genes; 1,567 -> 457 reads). 80 copies over 26 genes = 3.08 loci/gene; 3,606 / 11,977 = 30.1% of reads touch >= 2 loci of one gene; 93.3% within-gene. Co-located merge only sees overlapping loci (26 of 34 fragment pairs disjoint). Offline K=10: 96 -> 42 loci, 3 cross-gene merges (T8).
(full: archive/2026-11/o1_ledger.md:8582)

## §6ck — ⚖️ THE LINKED-LOCUS MERGE, IMPLEMENTED AND RUN: right signal, WRONG SITE (09-02)
RUSTLE_LOCUS_LINK_MIN_READS shipped default OFF (byte-identical OFF, params row, o1-perp-o2 disclosure). K=10: 96 -> 93 copies (offline predicted 42); 76/79 linked pairs span different E_r blocks, so merge site is downstream of the split. Guard caught refactor blinding its scan. Suite 829 passed / 0 failed / 11 ignored.
(full: archive/2026-11/o1_ledger.md:8662)

## §6cl — ⭐⭐⭐ WHY THE BLOBS FORM: the canonical-junction filter SELECTS FOR the mis-chains it cannot see (09-02)
PARTLY RETRACTED by §6cm. Blob-node predictors weak (exon density 2.7x, giant-intron fraction 2.2x; max intron 0.96x); first framing circular. Canonical >10 kb filter keeps 6/6 canonical giant introns in models vs 98/234 = 41.9% in reads: filter selects for well-formed mis-chains. Intron-vs-exon-content bound = next test, not finding.
(full: archive/2026-11/o1_ledger.md:8743)

## §6cm — ⛔⛔⛔ RETRACTION: the "mis-chained giants" are REAL GENES, and my headline read count was 3.4× inflated (09-02)
RETRACTION: two of three "mis-chained giants" are real genes (SNX29 covers 100.0%, PDXDC1 99.9%); both candidate rules (intron/exon_sum, junction support) failed. NPIPP1 "1,567 spliced reads" RETRACTED -> 457 with aligned block (upper bound; 71.6% splice over). §6cj linkage, §6ck barrier, §6cl canonical mechanism stand.
(full: archive/2026-11/o1_ledger.md:8818)

## §6cn — ⭐⭐⭐ THE BUG-CLASS SWEEP: 8 classes, 28 candidates, 16 killed, 12 real — two of them mine (09-02)
Adversarial bug-class sweep: 28 candidates -> 16 REFUTED -> 12 survivors. Own errors: §6bx band stats (only precision 0.9739 is band; recall/33.7x pooled); seeded_mode_tier.py span bug (§6be 152/305 -> 123/305 = 40.3%). Self-matching guard fixed with fn_body(). 7 further survivors recorded unfixed; 5 candidates UNRESOLVED (agents died).
(full: archive/2026-11/o1_ledger.md:8880)

## §6co — ⛔⛔ THE EXON SUM CANNOT BE "APPROXIMATED" TO A GENOMIC COVERAGE — THE RE-NORMALISATION INVERTS THE DISCRIMINANT (09-02)
Can the genomic denominator be approximated from exon sum? NO: all four re-normalisations invert TP/FP ordering (TP median 0.7200 -> 0.3201, FP 0.5745 unchanged), because single-exon FP nodes have exon_sum/span 1.0000. Missing info is in the numerator. RUSTLE_ER_CORE_COVERAGE vacuous here (core == span 1533/1534). Closed, negative.
(full: archive/2026-11/o1_ledger.md:8954)

## §6cp — ⭐⭐⭐ THE TWO COVERAGES MUST BE USED TOGETHER: `span + coverage-longer 0.30` IS A PARETO WIN ON NPIP (09-02)
Five real-binary arms on soto_adj, NPIP one family (120 truth pairs). arm_spancl30 (genomic span + cov-longer 0.30): P 0.8558 -> 0.9722, R 0.5667 -> 0.6000, Pareto; span alone P 0.6364. recip50 rule near-vacuous (P = 1.0000 all arms). NOT shipped (coverage clause rejected on wider Soto adjudication). Qualified by §6cq.
(full: archive/2026-11/o1_ledger.md:9020)

## §6cq — ⚠⚠⚠ CORRECTION TO §6co/§6cp: THE 22 NPIP "FALSE MERGES" ARE A LABELLING ARTEFACT AT ANY-OVERLAP — AND THE ISOFORM UNION IS REFUTED BY CONSTRUCTION (09-02)
Corrections: FP count is any-overlap labelling-sensitive (P 0.7528 at floor 0 -> 0.9804 @0.10 -> 1.0000 @0.20); Gate G prereg (>= 10 of 22 must be genuine over-merges). §6cp Pareto = dominance, not significance; recall rise is gamma non-monotonicity. Isoform union REFUTED (disarms stub flag, inflates true members ~3x harder, 8/22 untouchable). Flag-masking interactions and source defects listed. (Floor reading inverted by §6cr.)
(full: archive/2026-11/o1_ledger.md:9079)

## §6cr — ⭐⭐⭐ GATE G: **PASS, 22/22**. THE NPIP FALSE MERGES ARE REAL, THE MECHANISM IS **Alu**, AND THE OVERLAP FLOOR IS DISQUALIFIED (09-02)
Gate G PASS: O1_ERROR 22/22, O3 0/22; precision 0.7528 unchanged. Fraction-of-gene floor DISQUALIFIED (removes labels, not errors). Mechanism Alu: edges only at -k11 -w5, identity 0.710-0.852 vs genuine 0.941-0.980; E_c 0 shared reads. Denominator caveat: 34 one-endpoint edges, P bounds 0.6634-0.7822. Truth set wrong on 11/22 labels; NPIPB12 BED 409 bp short (lead).
(full: archive/2026-11/o1_ledger.md:9175)

## §6cs — ⭐⭐⭐ THE FULL CAUSAL CHAIN: oligo-dT INTERNAL PRIMING ON Alu A-TAILS -> UNSPLICED PRE-mRNA "LOCI" -> Alu-Alu EDGES. AND THE PIPELINE UPPERCASES THE SOFT-MASK AWAY (09-03)
Chain: oligo-dT internal priming on Alu A-tails -> unspliced single-read stubs -> Alu-Alu edges. 51/67 true vs 0/22 false survive hard-masking (Fisher p = 5.5e-11); NPIP core stays connected. 72.1% unspliced primaries end on A-tract (null 1.1-5.0%). Stub node median span/read 1.11. GenomeIndex uppercases soft-mask away; repeat-hub veto uses asm20 (blind). >=0.90 identity cut post-hoc, not to ship.
(full: archive/2026-11/o1_ledger.md:9280)

## §6ct — ⚖️ PRE-REGISTRATION: THE REPEAT-JUSTIFIED EDGE RULE (`RUSTLE_ER_REPEAT_MASKED_EDGES`) — written 09-03 BEFORE any arm
Prereg RUSTLE_ER_REPEAT_MASKED_EDGES (zero fitted params: unchanged E_r rule on hard-masked reps). P1 NPIP precision 0.7528 -> >= 0.95; P2 13 genes one component; P3 held-out 82-family F1 >= OFF baseline (decides shipping); P4 survivors lower soft-mask. Ship iff all four; default OFF, params row mandatory.
(full: archive/2026-11/o1_ledger.md:9382)

## §6cu — ⚠⚠ CORRECTION TO §6cs's BAM PROVENANCE (the -M -L trap, caught 09-03). SUBSTANCE UNAFFECTED.
Correction: soto_adj BAM is soto.bam, an -M -L subset (362 regions, 14,882,755 bp) of A119b.t2t.bam, not the full BAM. Substance unaffected (node 891/893/1178 counts identical in both). Library stats were on full BAM, so priming mechanism is a library property.
(full: archive/2026-11/o1_ledger.md:9443)

## §6cv — ⚖️⚖️ THE REPEAT-JUSTIFIED EDGE RULE: IMPLEMENTED, RUN, **P1/P2 PASS, P3 SPLIT ⟹ STAYS DEFAULT-OFF** (09-03)
Implemented; 831 tests pass. Mask 51.0% of bases, 573/1,571 edges survive (as predicted). P1 PASS 51/0 = 1.0000; P2 PASS (genes 13, NPIP families 10 -> 3); P3 SPLIT (unconditioned F1 0.4540 -> 0.5668 any-overlap but 0.2788 -> 0.2559 at >=50%; prereg F1 underspecified); P4 vacuous. Verdict DEFAULT-OFF. 32 Soto members lost at >=50%. Superseded by §6cw.
(full: archive/2026-11/o1_ledger.md:9483)

## §6cw — ⛔⛔⛔ THE REPEAT-MASKING RULE IS **MECHANISTICALLY BACKWARDS**. IT PREFERENTIALLY KILLS RECENT DUPLICATIONS. DO NOT SHIP — AND THE USER'S CONTAINMENT RULE IS VINDICATED (09-03)
Masking is mechanistically backwards: killed edges incident to lost Soto members median identity 0.9728 vs survivors 0.8143 (5.4% >= 0.94). Recent SDs copy their repeats. Containment cov_longer >= 0.30 would spare 59/76 = 77.6%. Verdict: RUSTLE_ER_REPEAT_MASKED_EDGES must not ship (kept as labelling instrument); pursue containment; §6cv P1 not a method result.
(full: archive/2026-11/o1_ledger.md:9553)


## §6cx — ⛔⛔⛔ THE INTERNAL-PRIMING READ FILTER IS REFUTED. **EVERY STUB-TARGETING LEVER IS REALLY A LOW-DEPTH FILTER, AND RANDOM THINNING REPRODUCES IT** (09-03)
Internal-priming read filter: both refuters FATAL. Mechanism real (60.1% of unspliced primaries end on A-tract vs 3.8% strand-inverted) but in Soto SD regions 1.72x not 3.37x. Filter dies: negative controls fire (+300 bp null 27.4x > real 8.0x); depth confound (FALSE_NPIP median depth 2.5 vs TRUE 29; random thinning 10.41x vs observed 7.99x); global F1 0.7789 -> 0.7550.
Lesson: no stub-targeting rule may be reported without a depth-matched / random-thinning null. Containment remains the only surviving lever.
(full: archive/2026-11/o1_ledger.md:9612)

## §6cy — ⭐⭐ CONTAINMENT PASSES THE DEPTH NULL THAT KILLED EVERY OTHER LEVER — AND IT IS ALREADY IMPLEMENTED (09-03)
Containment (`RUSTLE_ER_COVERAGE_LONGER_FLOOR`, `arm_cl30`) separates within depth strata (FP median cov_longer 0.1655/0.2004 vs TP 0.3847/0.4461). But `min_depth >= 10` alone reaches precision 1.0000 (38 TP / 0 FP) vs cov_longer>=0.30 0.9545 (63/3): containment's advantage is recall.
Held-out 82 Soto: arm_cl30 F1 uncond 0.5091 vs arm_off 0.4540; flat at >=50% floor. Nothing to implement; §6bs rejection not overturned; precision/recall trade is the user's call.
(full: archive/2026-11/o1_ledger.md:9682)

## §6cz — ⭐⭐⭐ THE VG MINIMIZER-MULTIPLICITY ORACLE CLEARS **BOTH** NULLS AND IS THE FIRST CATEGORICAL SEPARATOR: FP median = **0.0** SHARED PRIVATE NODES IN EVERY STRATUM (09-03)
VG minimizer-multiplicity oracle (`bench/soto/vg_mult_oracle.py`, multiplicity over reps): keep edge iff >= 1 shared node of multiplicity < 5 (`M_OP = 5`, not fitted). TRUE 59/67, FALSE 1/22, precision 0.9833; FP median 0.0 low-mult shared nodes in every depth and length stratum.
Offline (T8) hypothesis only; does not overturn register's AUC 0.686 on GORILLA. NEXT: score held-out 82, then a binary arm.
(full: archive/2026-11/o1_ledger.md:9745)

## §6da — ⭐⭐⭐ MCL ON THE ANNOTATION + RNA CORROBORATION: NPIP RECOVERED **31/31** WITH **35/36 MEMBERS LOC-NAMED**, AND "DNA PROPOSES, RNA DISPOSES" IS MEASURED (09-03)
MCL on gorilla npip3 annotation (4,477 genes, weight identity x cov_longer; 903 nodes / 3,169 edges): NPIP 31/31 at I=1.4; principal cluster 24 members, 0 NPIP-symboled, 35/36 LOC-named across 3 NPIP clusters.
RNA corroboration bimodal: 50/98 clusters zero; NPIP 22/24 = 0.917; zero-group median density 1.000 vs 0.793 ("a repeat clique is a perfect clique"). Caveat: zero conflates clique with not-expressed. Counts later SUPERSEDED by §6dd.
(full: archive/2026-11/o1_ledger.md:9828)

## §6db — ⛔ THE EXTRA-LOCUS HALF: **0/2 SURVIVE.** BOTH ARE NEIGHBOURS' UNANNOTATED EXONS — AND THE FILTER'S DEFECT IS NAMED (09-03)
Extra-locus enumeration on 27 corroborated clusters: 2 candidates survived, both adjudicated (ARTEFACT = intron 1 of LOC101130877; ANNOTATED_ALREADY = read-through of LOC101151324). NO NOVEL LOCUS IS CLAIMED.
Filter defect named: discard should be "reads splice into a neighbour's annotated junctions", not span overlap. §6da definition half stands.
(full: archive/2026-11/o1_ledger.md:9891)

## §6dc — ⚠⚠⚠ MY MCL COVERAGE DENOMINATOR EXCLUDED 62.8% OF GENES BY CONSTRUCTION; AND THE AD HOC REFERENCE **CANNOT** MOVE THE O3 GATE (09-03)
Defect in §6da: cov_longer over gene SPAN removes 1,525/2,428 = 62.8% of eligible genes; exonic denominator takes graph 903 -> 1,850 nodes, no re-alignment needed. NPIP result stands.
Ad hoc reference cannot move gate 5 (query and target both genome-derived; GWFAM247 reconstruction 0.9961 vs true 0.9754); `linearize_certificate` already is stage 3 (opt-in). 73 families: 8/328 admitted. Corrections: "0 of 289" is MEC-panel-only; LIVE BUG synthetic tid keyed on host (`copy_split.rs:441`).
(full: archive/2026-11/o1_ledger.md:9952)

## §6dd — ⭐⭐ THE EXONIC DENOMINATOR, FIXED AND SWEPT: §6da SUPERSEDED. NPIP 31/31 WITH **26** TRUTH LOCI IN ONE CLUSTER, NO SUPERFAMILY (09-03)
Exonic denominator fixed (GFF 1-based join 4,460/4,477 = 99.6%): pilot 1,920 nodes / 24,286 edges. I = 2.8-3.2 plateau: NPIP 31/31, main cluster 43 with 26 truth loci, none >100; cliff at 3.6 (18/31).
§6da counts SUPERSEDED; qualitative findings stand pending re-measure. Lesson: a no-op result signals a failed join.
(full: archive/2026-11/o1_ledger.md:10007)

## §6de — ⭐⭐⭐ GENOME-WIDE: NPIP TRANSPORTS EXACTLY, THE SUPERFAMILY STAYS DEAD, AND "A REPEAT CLIQUE IS A PERFECT CLIQUE" REPRODUCES **SIZE-CONTROLLED** (09-03)
Genome-wide gorilla (41,193 genes, I = 2.8): 19,896 nodes / 1,264,135 edges, 3,310 clusters; largest 112, one >100. NPIP 31/31, principal cluster 44 members / 26 truth / 42/44 LOC.
Clique signature size-controlled: zero-corroboration 391/1,101 density 1.000 vs 0.700, both median size 4.0. Naming 73.0% / 44.8% (not the product-census 95.2% / 67.9%). Weakest: only 229/738 = 31.0% product families fully recovered.
(full: archive/2026-11/o1_ledger.md:10060)

## §6df — ⛔⛔⛔ `--vg-realign-correct` IS LENGTH-CONFOUNDED: A SHORT COPY IS EXCLUDED BY ARITHMETIC AND ITS READS ARE REASSIGNED **362/362 = 100%** UNCONDITIONALLY (09-03)
`--vg-realign-correct` on MCL `--families`: `aln_id` capped at len(t)/len(q) (721 bp copy cap 0.269) and `n_decisive` from edit distance -> 362/362 = 100.0% of linear_copy=1 rows reassigned unconditionally; 583/1,242 = 46.9% molecules affected.
VERDICT: DO NOT USE `--vg-realign-correct` until `aln_id` is length-normalised (default-off). Four own errors corrected (join many-to-one, tie-break rows).
(full: archive/2026-11/o1_ledger.md:10114)

## §6dg — ⭐⭐ §6df FIXED AND VERIFIED: THE UNCONDITIONAL-REASSIGNMENT PATHOLOGY COLLAPSES **362/362 -> 9/362** (09-03)
Fix: `aln_id_len_safe` + `psv_decisive_count` (PSV columns). Reassigned 626 -> 294; 721 bp copy reads reassigned 362/362 -> 9/362 = 2.5%; unique-mapper agreement 68/68 = 100% unchanged; 831 tests pass, fixtures clean. Fixtures now carry two PSV columns.
(full: archive/2026-11/o1_ledger.md:10177)

## §6dh — ⭐⭐ DISPERSED FAMILIES UNBLOCKED: THE 38-COPY / 89.5 Mb FAMILY NOW ASSIGNS, 292/293 UNIQUE-MAPPER AGREEMENT (09-03)
`--families` read loading from union of copy neighbourhoods (`COPY_READ_PAD = 50 kb`, 89.5 Mb hull -> ~3.2 Mb), dedup load-bearing. MCLFAM1 (38 copies) completes: unique-mapper agreement 292/293 = 99.7%, but 3,407/4,036 = 84.4% tied, 8 of 38 copies get decisive reads: "processable, not resolved".
`--vg-realign-correct` does not scale (needs minimizer prefilter).
(full: archive/2026-11/o1_ledger.md:10218)

## §6di — ⭐⭐⭐ THE OVER-MERGING WAS THE **PARTITIONER**, NOT THE READS: MCL ON THE SHIPPED E_r EDGES GIVES 63 TP / 4 FP WHERE γ-QC GIVES 67/22 (09-03)
Same shipped E_r edges (soto_adj/arm_off), only partitioner changed: γ-QC 67 TP / 22 FP (0.7528); MCL I=1.4 63/4 (0.9403); I=2.8 47/0 (1.0000). Offline (T8); inflation depends on graph density, re-sweep per substrate. NEXT: binary arm replacing γ-QC with weighted MCL.
(full: archive/2026-11/o1_ledger.md:10267)

## §6dk — ⭐⭐⭐ O3 POSITIVE CONTROL REPRODUCED FROM THE MCL ROSTER: GOLGA6L7 DISCORDANT BY ONE UNIT, **9/9 CONTROLS CONCORDANT**, TANDEM ARRAY MATCHING BASE-FOR-BASE (09-03)
Matched-individual O3 positive control with MCL-roster probes vs `o3_hapcnv/{mat,pat}.chr.fa` (KB3781): GOLGA6L7 5 vs 6 (discordant by one unit; recorded 6 vs 7); 9/9 autosomal controls concordant; tandem array coordinates identical to §6u.
Lesson: address a family by its cluster, never by a symbol. Open: absolute count probe-dependent.
(full: archive/2026-11/o1_ledger.md:10310)

## §6dl — ⭐ §6dk's OPEN QUESTION CLOSED: THE MISSING LOCUS IS THE 3,182 bp UNIT AT `cov_q` 0.485 vs 0.508 — THE **DIFFERENCE** IS PROBE-STABLE, THE **COUNT** IS NOT (09-03)
Missing GOLGA6L7 locus = 3,182 bp unit at cov_q 0.485 (full-span probe) vs 0.508 (§6u). pat - mat = +1 across cov_q 0.35-0.50; counts 6 vs 7 at 0.45-0.35. Rule: quote the DIFFERENCE; state the probe for absolute counts.
(full: archive/2026-11/o1_ledger.md:10372)

## §6dm — ⭐ TWO-SIDED COVERAGE IS INERT ON THE HAPLOTYPE SCAN, AND THE REASON REFRAMES THE COUNT: `cov_q` IS A COMPLETENESS CRITERION, NOT AN ARTEFACT FILTER (09-03)
Adding cov_locus >= 0.50 changes nothing at every floor; boundary unit cov_locus 0.986 vs cov_q 0.485, a complete match to a shorter unit: cov_q is a completeness criterion. Haplotype CNV in DNA is real; O3 (RNA detection) still unproven. Next: "can the RNA see THIS one?"
(full: archive/2026-11/o1_ledger.md:10404)

## §6dn — ⛔ THE HAPLOTYPE O3 SCREEN PROMOTES NOTHING: ALL 4 INSPECTED CANDIDATES FAIL, AND THE MISSING FILTER IS **MAPQ** (09-03)
Directional haplotype screen, 864 probes: 733/864 = 84.8% concordant; interaction test Fisher OR 1.82, p = 0.159 (not a result). Alignment contig map 23/24 clean. All 4 inspected candidates FAIL (LOC134758568 extras MAPQ 0, pat copies on unexamined contig). Recommended MAPQ > 0 (retracted in §6do). O3 status unchanged.
(full: archive/2026-11/o1_ledger.md:10442)

## §6do — ⭐⭐ A CLEAN REFERENCE-ABSENT CANDIDATE (`LOC134757045`, 2 maternal units), AND ⛔ §6dn's MAPQ FILTER RETRACTED BY THE KNOWN POSITIVE (09-03)
Retraction: MAPQ > 0 would delete GOLGA6L7's true extra unit (MAPQ 0). Discriminator is identity (>= 0.94). Candidate `LOC134757045` (lncRNA, CL493): mat 5 vs pat 3 units, 3 one-to-one, 2 maternal units cov 1.000 id 0.9961 absent from PAT-derived reference.
Haploid representation, not assembly collapse; no between-individual comparison. Next: RNA test, uninformative if silent. Standing doc `docs/O3_STATUS.md`.
(full: archive/2026-11/o1_ledger.md:10498)

## §6dp — ⛔ TWO-ANIMAL TEST: NO COLLAPSED EXTRA COPY AT THE TWO BEST-POWERED LOCI, AND THE READ-PRESENCE SCREEN IS **UNUSABLE** BETWEEN TISSUES (09-03)
Two animals (testis GGO_OR6737 vs fibroblast donor): read-presence screen UNUSABLE between tissues (33/52 differ, 16 vs 17, p = 1.0000; differential expression of existing copies). Allele-count test at CL100 (6,482 pos) and CL106 (2,678): 0 positions with >= 3 alleles in either animal. Two loci only; exclude mito contig.
(full: archive/2026-11/o1_ledger.md:10548)

## §6dq — ⛔⛔ THE TWO-ANIMAL ALLELE TEST IS A CLEAN ZERO ON 226,615 POSITIONS — AND THE FILTER THAT MAKES IT VALID EXCLUDES THE GENE CLASS THE HYPOTHESIS IS ABOUT (09-03)
Allele test scaled to 99 loci / 49 families, 226,615 positions: 0 in every 3-allele category. But expressed-in-both filter selects housekeeping genes: 49 of 3,310 families testable, the wrong 49. Fix is matched tissue, not depth. Not evidence against between-individual CNV.
(full: archive/2026-11/o1_ledger.md:10601)

## §6dr — ⭐⭐ O3 EVIDENCE INVENTORY: **NOT REFUTED — UNTESTED IN THE REGIME.** "WELL TESTED BUT LACKING SUPPORT" IS NOT DEFENSIBLE (09-03)
3-agent inventory verdict: O3 NOT REFUTED, untested in the regime; no shipped binary can emit a reference-absent copy at defaults (all four admission routes default-false). Corrected "289 candidates, 0 admitted" and "AN INDEPENDENT ZERO" in docs. Standing doc `docs/O3_STATUS.md`.
(full: archive/2026-11/o1_ledger.md:10643)

## §6ds — The MCL coverage floor had a UNIT MISMATCH; the remedy is an ADDITIVE exonic floor, not an exonic numerator (2026-09-04)
Defect: `annotation_families.rs` numerator in genomic-span bases vs exon-union denominator; 6,388/24,286 = 26.3% admitted edges have floor vacuous. MCL0 (84) = two cliques glued by one mosaic pseudogene; MCL2 (33) = repeat clique. `--exonic-overlap` shatters NPIP 43 -> 19/14/4/4 (wrong remedy).
`--min-exonic-bp 300` (additive) keeps NPIP 43/43, MCL3 22/27, splits MCL0, MCL2 1/33; 300 not pre-registered. Shipped all three clauses DEFAULT OFF; 775 passed / 0 failed.
(full: archive/2026-11/o1_ledger.md:10671)

## §6dt — There is no unbiased `min_exonic_bp` to choose, because the parameter DISSOLVES (2026-09-04)
Plateau 1-200 bp invariant; 300 costs MCL3 five members (min_bp anchor REFUTED). 53,305/63,361 = 84.1% pairs share exactly zero exonic bases: defensible clause is ">= 1 shared exonic base" (no free number); MCL3 27/27, MCL2 0/33. NEW101 (101-member cluster) unadjudicated. Default remains `min_exonic_bp: 0` (OFF); honest flip value is 1.
(full: archive/2026-11/o1_ledger.md:10732)

## §6du — NEW101 ADJUDICATED: not a family. And `frac_in` is the instrument `density` was missing (2026-09-04)
NEW101: frac_in 0.557, 94.6% intronic, 814 exonic bp -> NOT A FAMILY. Raising floor worsens cohesion. Shipped `frac_in` cohesion certificate column (reported, never a filter). Claimed 20/46 = 43.5% n>=5 clusters below 0.75 ("~40% hairball slices") - RETRACTED in §6dx. 776 passed. Defaults unchanged.
(full: archive/2026-11/o1_ledger.md:10783)

## §6dv — "DNA PROPOSES, RNA DISPOSES" IS NOT IMPLEMENTED: corroboration ANTI-correlates with cohesion (2026-09-04)
Pre-registered (`docs/archive/2026-09/PREREG_rna_corroboration_2026-09-04.md`), 4/4 predictions HELD: corroboration inverted, cohesive 0.8076 vs slice 0.8868 (Δ −0.0792), inversion grows with size; MCL0 real array 0.0000, NPIP 1.0000. "RNA disposes" is NOT implemented by corroboration; role = necessary condition only. Pair-level RNA = E_c, forbidden for O1.
(full: archive/2026-11/o1_ledger.md:10839)

## §6dw — STRUCTURAL TRAITS: the array signature is real but TYPE-SPECIFIC; readthrough is REFUTED (2026-09-04)
Annotation-only traits: readthrough REFUTED (real NPIP has 14 overlapping pairs). Array signature type-specific. Tandemness weak: r(tandem, frac_in) = +0.2211, OR 2.79; r(tandem, corroboration) = -0.1535. Core: 16/46 n>=5 clusters with frac_in >= 0.90 AND tandem >= 0.6 (post hoc cuts). Nothing ships.
(full: archive/2026-11/o1_ledger.md:10890)

## §6dx — RETRACTION of §6du's "43.5% hairball slices"; SEDEF is the instrument; λ and cross-species REFUTED (2026-09-04)
RETRACTED §6du "43.5% hairball slices" (circular). SEDEF: 17/20 low-frac_in clusters >= 0.64 SD-supported; artefacts are 3 SD-empty clusters (MCL28, MCL15, NEW101 23/661). MCL3 REAL label withdrawn (SEDEF 0.215). λ >= 2 REFUTED both directions; small-n chance P(≥5 of 14) = 0.91; cross-species replication REFUTED as discriminator.
§6dd plateau is NPIP-only. Register rows 300/301/286 already kill graph signatures. Next: adjudicate disagreement cells.
(full: archive/2026-11/o1_ledger.md:10937)

## §6dy — THE DISAGREEMENT CELLS ADJUDICATED, THEN ATTACKED: 39/46 SD-corroborated, 3 artefacts, BOUNDARIES certified for 3 (2026-09-04)
46 n≥5 clusters on (frac_in, frac_sd_in), disagreement cells pre-registered and attacked (wf_171430db-a4a). SEDEF floor corrected (fracMatch min 0.507, aln 900 bp). Cell A (17): 9 cuts by pre-registered rule (not 10/13); cuts uncertified; MCL23/44/33 separate families; do not lower inflation; tier-2 NEEDS-WORK, not shipped.
Cell B (3 SD-empty): member union-exonic-coverage calibrates (MCL28 0.000, MCL15 0.011, NEW101 0.000); MCL15 NEEDS-WORK. Cell C MCL4: REPEAT-WELDED pending a repeat library.
Tally: REAL-by-SEDEF 39 (3 anchored), ARTEFACT 3, UNDECIDABLE 1, not adjudicated 3; boundaries certified for 3 of 46. Hash every PREREG.
(full: archive/2026-11/o1_ledger.md:11017)

## §6dz — REPEAT LIBRARY ACQUIRED (UCSC GenArk, Dfam 3.8 + RepeatModeler): cell B CONFIRMED, MCL4 = a ≥673-copy UNCLASSIFIED young element, de novo library CANNOT arbitrate, gene-size rule = flag only (2026-09-04)
UCSC GenArk RepeatMasker (Dfam 3.8 curated) + RepeatModeler acquired. Cell B CONFIRMED: MCL28 L1 6/6, MCL15 ERV1 22/22, NEW101 L1 639/661 edges ≥0.80; FALSE clique SVA; REAL anchors tandem 0.000, NPIP 0/582. MCL4 0/172 on curated.
De novo prereg R1 FAILED (tandem array absorbed 1105/1105): de novo library is not an arbiter (register row 665). MCL4 = ≥673-copy young unclassified element (`rnd-1_family-7`), ARTEFACT-class for O1, flagged for O3. Gene-size rule P1 FAILS on all four statistics: ship as per-member flag only. Addendum: 7/8 giant NPIP models homologous end to end; 1/8 probable chimera.
(full: archive/2026-11/o1_ledger.md:11217)

## §6ea — "SHOULD RNA GET MORE PROTAGONISM IN DEFINING FAMILIES?" — TESTED ON NODE BOUNDARIES: NO. The giant NPIP models are transcription units (2026-09-04)
Pre-registered (`mcl_ann/adj/rnanode/PREREG.md`): RNA-FIXES REFUTED (1/8 giants shattered, needed 6); ANNOTATION-RIGHT holds 6/8 but test formally VOID (504 kb unique control scores 0.25). Read chains make boundaries less copy-consistent (NPIP controls 1.64, 5.68). Answer: no larger RNA role. Register row 667.
(full: archive/2026-11/o1_ledger.md:11341)

## §6eb — CATALOG AUDIT OF ALL 46 CLUSTERS: whole-gene closure measures INTRONS; the exonic form separates the anchors but is POST HOC; two SD-supported clusters are ERV-K (2026-09-04)
Pre-registered genome-wide closure audit (`audit46.tsv`): every prediction FAILED (P1 15/22 LOCAL; P2 10/17; P3 4 flagged) because whole-gene closure counts intronic repeats. Post-hoc exonic form separates all anchors/artefacts from MCL4 and MCL24 (must be pre-registered on `gw_cl_I2.8.tsv`). MCL32/MCL24 (U4) are ERV-K (R1 1.00/0.97). Nothing relabelled; register rows 668–669.
(full: archive/2026-11/o1_ledger.md:11408)

## §6ec — THE MCL PRUNE THRESHOLD IS SIZE-DEPENDENT: it dissolved the anchored tandem array genome-wide, and §6dd's "cliff at I=3.6" is a PRUNE artefact (2026-09-04)
Genome-wide exonic-closure prereg stopped at its fail conditions: tandem anchor unclustered genome-wide (84/84). Cause: absolute post-inflation prune 1e-5 empties near-uniform cliques with n+1 > 10^(5/I) (61 nodes at I=2.8). Prune 1e-9 prereg predictions all met (tandem 60+58; 138/142 3-contig clusters identical; largest 203, 3 clusters >100).
§6dd plateau/cliff is a PRUNE artefact; §6de "zero clusters over 100" was the prune. Shipped `--prune` (default 1e-5 unchanged). Register rows 670–671.
(full: archive/2026-11/o1_ledger.md:11467)

## §6ed — DEFAULT PRUNE FLIPPED TO 1e-9 (user decision); canonical catalogs regenerated; the genome-wide closure flag turns out to certify MCL CUTS, not elements (2026-09-04)
User decision: `--prune` default 1e-5 -> 1e-9. Canonical catalogs `mcl_ann/rna_bp1_p9.clusters.tsv`, `mcl_ann/gw_bp1_p9.clusters.tsv`; 138/142 clusters identical, no class change. Genome-wide exonic-closure P2/P3 FAIL: every genome-wide flag is an MCL cut, not an element -> candidate cut certificate (to pre-register). Register row 672. §6dd/§6de plateau statements superseded by §6ec.
(full: archive/2026-11/o1_ledger.md:11528)

## §6ee — O2 ON NPIP AT SCALE: it resolves the reads minimap2 already resolved, abstains on 99% of the multimappers, and half its ties were O1 NODE defects (2026-09-04)
O2 on NPIP (39-copy family, 222 s / 18.8 GB): 5,016 assigned, agreement 4,968/4,969, but 4,960 of 5,016 are MAPQ-60; multimappers 8 of 592 (1.4%). Copy 13 holds 4,743 assignments. Spearman −0.55 paralogue identity vs assigned share. 607 of 1,221 tied reads on overlapping copies (O1 node defect); keep-longer rule resolved 9/39 -> 16/27.
"O2 on NPIP is an abstention certificate, not a redistribution." Register row 673.
(full: archive/2026-11/o1_ledger.md:11575)

## §6ef — STEP O1-1: A LOCUS IS THE NODE. Shipped as `--merge-overlapping-loci` (OFF); measured; the NPIP boundary moves (2026-09-04)
Step O1-1 shipped: `mcl_families --merge-overlapping-loci` (OFF; byte-identical control). Exon-union-overlap rule folds 3,864 (9.5%) records into 3,154 loci; edges attributed to representative. Genome-wide 1,095 clusters / 8,893 members.
NPIP boundary moves: 48-member cluster = 33 NPIP loci + MCL7's 11 + 4. O2 locus-rule NPIP: resolved 17/43, MAPQ-0 assigned 1/78; 20.7 GB peak of 25. Default decision left to user.
(full: archive/2026-11/o1_ledger.md:11644)

## §6eg — WHAT NPIP IS, BIOLOGICALLY: the literature and a CHM13 cross-mapping say the locus-rule merge is the DUPLICATION BLOCK, not the gene family — and a quarter of the anchored "NPIP 44" are block chimeras (2026-09-04)
Literature: NPIP = LCR16a core duplicon; SMG1P = adjacent LCR16u (distinct origin). CHM13 cross-mapping: NPIP 44 -> 28 on NPIP loci, 7 ABCC1, 2 SORL1 (block chimeras); MCL7 -> 8 SMG1/SMG1P. Copy 13 (4,743 O2 assignments) is EIF3C/EIF3CL.
Attribution merges LCR16a with LCR16u (the duplication BLOCK): do NOT flip `--merge-overlapping-loci` with attribution (register row 676); re-adjudicate row 675.
(full: archive/2026-11/o1_ledger.md:11692)

## §6eh — SUGGESTIONS APPLIED: CHM13 landing for all 46 clusters; O2 re-read without EIF3C; and SEDEF recovers the ~23 kb LCR16a CORE inside every NPIP model — "define the duplicon first" WORKS (2026-09-04)
O1-9 CHM13 landing for all 46 clusters: real families one stem one chromosome; artefacts scattered; 16p blocks mixed stems. O2-7 without EIF3C: NPIP-proper MAPQ-0 0/68, MAPQ 1–59 8/524.
SEDEF core pilot: NPIP core ≥ half median 23,420 bp (LCR16a), EIF3C 808 bp, MCL7 0. "Define the duplicon first" = segmentation refinement over MCL clusters (O1-10), not replacement. Addendum: applies only to young SD families (old ZNF/OR clusters core = 0); core as fraction of exon-union to be pre-registered.
(full: archive/2026-11/o1_ledger.md:11758)

## §6ei — O1-10 IMPLEMENTED: `mcl_families --core-refine --sedef` (OFF). Acceptance: the TRIM arm holds on NPIP, the DROP arm FAILS on MCL5 and MCL25 (2026-09-04)
Pre-registered (`mcl_ann/adj/core/PREREG.md`): `--core-refine --sedef` (OFF), one constant "half". TRIM arm holds (NPIP chimeras cut to 23,228–25,456 bp hulls; median core 23,612 bp; MCL7 drops exactly PDXDC1, SLC7A5). DROP arm FAILS on MCL5 (2 real ZNF dropped) and MCL25; register row 677. Genome-wide 524 dropped. `dropped` is a flag, not a removal.
(full: archive/2026-11/o1_ledger.md:11843)

## §6ej — O2 ON NPIP's LCR16a CORES, AND A JUNCTION-ANCHORED TRUTH: O2 abstains on 98% including minimap2's unique reads, and its 11 confident calls are 4/11 WRONG by the reads' own splice structure (2026-09-04)
O2 on 35 NPIP cores: 19 assigned / 3,347 tied / 41 ambiguous; MAPQ-60 14/535 assigned. Junction-anchored truth (calibration 52/52): 117 anchored reads, O2 assigns 11, agreement 7/11. Causes: copies are annotation models; min_p counts star-projected columns. O2-8 defined. Register rows 678–679. (4 wrong calls later shown to be annotation-gap anchors, §6eq.)
(full: archive/2026-11/o1_ledger.md:11880)

## §6ek — O2-8a IMPLEMENTED AND REFUTED: the junction/PSV conflict abstain cannot see the discordance, because O2's junction term lives in spliced space and the truth lives on the genome (2026-09-04)
Pre-registered `--junction-conflict-abstain` (OFF): fires on 228 already-tied reads, 0 of 4 discordant; C1 FAILS. O2 junction term lives in spliced space, truth on genome. Register rows 680–681. O2-8 reduces to O2-8b (read-supported locus copies).
(full: archive/2026-11/o1_ledger.md:11934)

## §6el — O2-8b PASSES: copies built from the locus core + READ-SUPPORTED exon chain remove every confident wrong call (2026-09-04)
Pre-registered (`adj/o2rc/PREREG.md`): read-chain copies: NPIP anchored 5/5 (was 7/11), assigned 19 -> 96; 3-copy control 57/57; D1 ✓ D2 ✓. Peak 25.4 GB of 25. Later found produced by defective bench builder (§6en, register row 682).
(full: archive/2026-11/o1_ledger.md:11962)

## §6em — O2 MEMORY FIXED: the PSV-discovery DP kept three full score matrices; now rolling rows + one traceback byte per cell, byte-identical (2026-09-04)
`banded_msa_pair` rewritten to rolling rows + one traceback byte per cell, byte-identical (200-pair randomized oracle test). 3-copy control 25.4 GB -> 0.49 GB; NPIP 43 loci 20.7 -> 5.7 GB. Sweep unblocked on memory.
(full: archive/2026-11/o1_ledger.md:12001)

## §6en — O1-10b SHIPPED (`mcl_families --emit-units`), and it EXPOSED §6el: the bench chain builder was defective, and faithful read-chain units send O2's star projection to min_p 1e-270 (2026-09-05, 00:xx)
Shipped `mcl_families --emit-units --bam --fasta` (OFF). §6el bench builder defective (register row 682). Faithful units: majority-splice rule 123 / 2,185 / 666, anchored 3/4, wrong call min_p 3e-270 — star projection over units of unequal exon structure. O2-8c (genomic core columns) defined; register row 683.
(full: archive/2026-11/o1_ledger.md:12023)

## §6eo — O2-8c FIRST FORM IMPLEMENTED (`copy_assign --psv-genomic`, OFF) AND NOT ACCEPTED: genomic alignment over unit EXTENTS manufactures columns from non-homologous flanks; the span must be the SEDEF core hull (2026-09-05)
Pre-registered `copy_assign --psv-genomic` (OFF): global alignment over unit extents manufactures columns from non-homologous flanks (pair 28 spliced vs 3,055 genomic columns); GFF cores 2/6. Not accepted; register row 684. Next O2-8c′: span = SEDEF core hull.
(full: archive/2026-11/o1_ledger.md:12067)

## §6ep — O2-8c′ IMPLEMENTED: hull-span genomic PSV columns via minimap2 are CORRECT PAIRWISE and still not a family-level fix; the star projection and the low-MAPQ junction truth are the two remaining suspects (2026-09-05)
Shipped `core_hull` column; minimap2 hull alignment correct pairwise (240 / 597 / 238 columns vs global DP 3,210 / 1,211 / 2,228; register 685). Family level: 12 assigned, anchored 3/9, min_p 0 underflow. Not accepted (E1/E2 fail). Suspects: star projection (O2-8d) and low-MAPQ junction truth (register row 686). Spliced mode stays default.
(full: archive/2026-11/o1_ledger.md:12104)

## §6eq — S1 DONE: the junction truth audited. 55/117 anchors were annotation gaps; on the VALID truth the spliced mode makes NO wrong confident call on NPIP, the genomic modes do (2026-09-05)
Truth audit: 117 anchored -> 62 valid, 55 ambiguous; 9 of 13 disagreements sit on ambiguous anchors (register 686 confirmed). On valid truth, GFF-model spliced 5/5; genomic hulls 3/8; genomic modes OFF, off critical path.
Addendum CORRECTION: units scores compared mismatched copy indices; exact rescoring: faithful units 0/3 (nested pair copies 9/10). "Faithful units 3/3" withdrawn; S2 (units per LOCUS) must precede; register row 687.
(full: archive/2026-11/o1_ledger.md:12147)

## §6er — S2 DONE: THE UNIT IS THE CATALOG ROW. `mcl_families` now emits locus units by default; O2 consumes them with no script in between (2026-09-05)
Defaults flipped: `--merge-overlapping-loci` ON (representative-only), core refine with `--sedef`, units emitted with `--bam --fasta`; escapes byte-identical. Canonical `mcl_ann/rna_units_v1`: 111 clusters / 784 members, 751 units, 313/751 with ≥0.99 paralogue. NPIP 29/33 in one cluster; MCL7 separate.
O2 on NPIP locus units: 62 valid anchors -> 0 assigned, 60 ambiguous; MAPQ-60 placement 75/76. S1 ✓ S2 ✓.
(full: archive/2026-11/o1_ledger.md:12201)

## §6es — S5 DONE: O2 over the unit catalog. 87 families, 16,723 unit reads: 29% assigned, 7% tied, 64% ambiguous; MAPQ-60 placement agreement 95%; the nearest-paralogue FORECAST does not generalise (2026-09-05)
Scoring correction via `family_join.tsv`: faithful units spliced 3/3 (row 687 amended). Sweep `mcl_ann/sweep_v1`: 87 completed, 4 aborted; 16,723 unit reads: 29% assigned / 7% tied / 64% ambiguous; MAPQ-60 placement 4,451/4,685 = 95.0%. Forecast Spearman −0.04: nearest-paralogue forecast does not generalise (register row 688).
(full: archive/2026-11/o1_ledger.md:12236)

## §6et — O2 HYGIENE DONE and the sweep re-run on the expressed-unit catalog: 35 families, 0 aborts, same picture (2026-09-05)
Shipped `in_copy`, `catalog_copy_idx` columns; units counted inside emitted chain (suite 852/0/11). Catalog `mcl_ann/rna_units_v2`: 291 units. Sweep v2: 35 completed, 0 aborted; 16,650 reads, 29% / 6% / 64%; MAPQ-60 agreement 4,538/4,761 = 95.3%; forecast −0.05.
(full: archive/2026-11/o1_ledger.md:12276)

## §6eu — THE SWEEP'S TWO WORST FAMILIES ADJUDICATED: three defects, two fixed as OFF flags, one open (2026-09-05)
Pre-registered (`adj/worst2/PREREG.md`). D1 unit was a transcript fragment -> `--units-follow-reads` (OFF, clamped to annotated locus; unbounded form engulfs, row 690); MCL38 agreement 119/310 -> 424/426. D2 PSV read-support filter drops monomorphic true PSVs (row 689); filter off 95.3% -> 96.0%; default still ON. D3 records as independent observations OPEN (row 691). v3 + filter off MAPQ-60 4,674/4,704 = 99.4%.
(full: archive/2026-11/o1_ledger.md:12305)

## §6ev — THE TWO O1 DEFINITIONS ON ONE TRUTH: the switch costs no measurable precision on Soto's band and gains recall through the node; the annotation node loses members the reads found (2026-09-05)
Pre-registered (`soto_mcl/PREREG.md`): band [0.90,1.00) precision old 0.974 vs MCL units 0.954 (CIs overlap); recall|both 0.874 -> 0.940; detection 0.354 -> 0.658. P1 held, P2 held, P3 FAILED (detection 0.66–0.70 vs >0.80). Recall gain is the node; RNA admission (G5) unbuilt. Details `docs/archive/2026-09/O1_DEFINITION_SWITCH.md`.
(full: archive/2026-11/o1_ledger.md:12355)

## §6ew — Q9 (NPIPA vs NPIPB) UNDER THE CORE RULE: the gorilla loci land on B stems, three on both, and the graph never cuts A from B before it fragments B (2026-09-05)
Pre-registered (`adj/q9/PREREG.md`): 29 NPIP loci: 17 NPIPB-only, 3 AB, 9 chimeras, none NPIPA alone. P1 held (AB–B 0.9015 vs B–B 0.9007), P2 coverage does not separate, P3 held (20 NPIP-proper loci one cluster I 1.4–4.0). Human A/B split not carried by gorilla loci; qualitative (3 vs 17).
(full: archive/2026-11/o1_ledger.md:12373)

## §6ex — G5: the annotation node's misses decomposed on Soto; RNA admission of unannotated loci prototyped; min_size 3 → 2 (user decision) (2026-09-05)
54 missed Soto members: 19 in size-2 clusters, 9 scattered (MCL fragmentation, row 694), 11 annotation gaps. min_size 2: detection 254 -> 272/362, precision 0.949 unchanged; user flipped `--min-size` default to 2; catalog `rna_units_v4`. RNA admission prototype: 5 admitted, recovers 2 of 11 gaps (FAM95B1 658 reads); bench only.
(full: archive/2026-11/o1_ledger.md:12395)

## §6ey — ROW 694 ROOT-CAUSED AND FIXED: fold loci AFTER clustering, exon-less records are one exon, homology must be exon-to-exon (2026-09-05)
Row 694 root cause: fold-first locus rule discards overlapping records of different families; 160/747 records exon-less. Fix: three OFF flags `--fold-within-clusters`, `--exonless-span`, `--exonic-both-sides` (pre-registered). All three: detection 289 (0.798), band precision 548/581 = 0.943 vs reference 0.949; P1/P3 partly missed by one member; exon-to-exon alone on old order 0.751 -> 0.699.
Gorilla `rna_units_v6`: NPIP all 44 records one cluster (32 loci); L1 blob MCL0 DISSOLVES. Register rows 696–698; suite 855/0/11.
(full: archive/2026-11/o1_ledger.md:12421)

## §6ez — DEFAULTS FLIPPED (five, user decision), catalog `rna_units_v8`, O2 sweep on it, the three new NPIP loci landed, the threshold plateau measured (G3), D3 pre-registered, duplication blocks emitted (2026-09-05)
User flipped five defaults (three fold/exon flags + `--units-follow-reads` ON; PSV read filter OFF); canonical `rna_units_v8` (274 clusters, 161 pairs, 389 units). Sweep v7: 104 families, 69% assigned, MAPQ-60 99.76%; paired 35 families 28.2% -> 33.7%; NPIP 9.1% assigned.
Three new NPIP loci landed on CHM13. G3 pre-registered: P1, P3 held; plateau identity 0.60–0.80, cov 0.10–0.50, 100–500 bp; walls 0.85 / 0.60 / 1 kb (row 699). D3 pre-registered. `sd_blocks` emitted: 83 blocks.
(full: archive/2026-11/o1_ledger.md:12486)

## §6fa — D3 IMPLEMENTED: pooling the records REFUTED, "the read is the star" + an ORIGIN CERTIFICATE ships (OFF), and NPIP's abstention is finally attributed — to the node, not to identity (2026-09-05)
D3 form 1 (pool records' PSV readings) REFUTED (row 700): 25 / 44 / 8 wrong anchor calls; 34% of a read's columns disagree between placements on NPIP. Form 2 `--molecule-observations` (read is the star) + origin certificate (Binomial on error rate; row 701), OFF. NPIP 175 assigned, MAPQ-60 143/144, anchors 0 wrong; paired 35: 33.2% assigned @ 99.60%; 10–30-site bin 16.5% -> 40%.
NPIP: 285 of 433 low-MAPQ reads (66%) explained by NO unit — abstention attributed to the node, not identity. Suite 857/0/11.
(full: archive/2026-11/o1_ledger.md:12536)

## §6fb — READ-STAR IS THE DEFAULT (user), the unexplained reads counted and located, units that share exon bases merged, and the genome-wide rebuild launched (2026-09-05)
User flipped `--molecule-observations` ON (`origin_rejected` column). v8 sweep tied 76%: nested units (MCL108) -> v7's 69% assigned was inflated (row 703). Fix `--merge-overlapping-units` (ON); catalog `rna_units_v9` (349 units, 40 merged). 14,415 of 71,137 reads origin-rejected; 62 candidate loci (roster gaps). RNA admission: 4 loci (LCR16u 35.99 Mb, 123 reads). Genome-wide rebuild `gw_units_v1`: 2,296 clusters, 3,755 units.
(full: archive/2026-11/o1_ledger.md:12607)

## §6fc — ROSTER ADMISSION BY READ ORIGIN (O2 → O1), the read-star reporting floor tested, the dissolved genome-wide cluster classified (2026-09-05)
gw_p9 MCL8 = L1-welded artefact. `RUSTLE_STAR_P` 0.3 vs 0 (row 704): kept 0.3 as measured operating point. Roster admission pre-registered: 21 of 74 families get 61 admitted copies; defect row 705 -> pairwise certificates. Opt-in `--origin-substitutions-only` (1 of 20 anchors wrong, row 707), `--read-star-junctions` (row 708).
Final default: paired 35 32.1% assigned, MAPQ-60 5,024/5,029 = 99.90%; NPIP anchors 0 wrong. Admission P1 fails (rejections −3%, not 50%); P2–P4 hold; bench prototype.
(full: archive/2026-11/o1_ledger.md:12636)

## §6fd — GENOMIC READ-STAR: the read against each candidate's LOCUS, splice-aware, with a whole-read origin certificate (2026-09-05)
Pre-registered `--read-star-genomic`: splice-aware alignment to each candidate locus, whole-read origin certificate (unaligned bases counted), locus padded by longest read. NPIP 480 assigned, 301/301, anchors 12 right / 0 wrong; LCR16u 42 -> 252. Row 706 corrected: rejections were node EXTENT.
All four predictions held: paired 35 assigned 9,622 (55.7%), MAPQ-60 9,236/9,236 = 100.0% -> genomic form is the default (user decision); `--read-star-unit` restores §6fc.
(full: archive/2026-11/o1_ledger.md:12688)

## §6fe — SPEED: where the read-star time goes, what was taken, what was not (2026-09-05)
~125 of 131 s is external minimap2. Taken: threads 2 -> 4 (`RUSTLE_STAR_THREADS`), batch restricted to block-in-copy molecules: paired 35 1,299 s -> 654 s, assigned 9,622 -> 9,630. Not taken: aligning only to BAM-placed loci (makes single-candidate ties; row 700 defect).
(full: archive/2026-11/o1_ledger.md:12724)

## §6ff — EXCISION: remove a copy, follow its reads, look for "a copy should be here" (O3 in O2's vocabulary) (2026-09-05)
Pre-registered excision (`bench/o2_excision.py`): P1 (no silent misassignment) holds 4 of 5, fails on NPIP 13 (12 reads absorbed by copy 12 at 99.85% local identity; row 711). NPIP 2 74% absent = orphan reporting gap (row 712). Consistent-site signature 7–100× controls; ZNF copy reconstructed 0.994. Wall = local identity within sequencing error.
(full: archive/2026-11/o1_ledger.md:12744)

## §6fg — ORPHANS REPORTED (row 712 closed), `n_candidates` column (2026-09-05)
Orphans reported as `ambiguous`, `origin_rejected 1`, `n_candidates 0`; new `n_candidates` column. NPIP 3,074 -> 3,307 rows (21 orphans); prior verdicts identical. Row 712 closed. Suite 857/0/11.
(full: archive/2026-11/o1_ledger.md:12778)

## §6fh — L1 + L2: dropped members are O2 candidates; the locus extent is O1's (PREREG `docs/archive/2026-09/PREREG_loose_ends_L1_L2_2026-09-05.md`, md5 78181c4f) (2026-09-05)
L1: dropped members emitted as units (`member_status = dropped`; 17 on 3 contigs, incl. EIF3C 4,774 reads). L2: `locus_start`/`locus_end` read-supported extent clipped at family neighbours (rows 713–714 refuted alternatives). Catalog `rna_units_v10` (366 units); P1 gate held.
Paired 35: assigned 11,368 (46.5%) vs 10,607 (43.4%); NPIP anchors 16/16, 0 wrong; P3 held; P2 literal ≥ 100 missed (row 715). EIF3C reads single-candidate ties -> L3 pre-registered next.
(full: archive/2026-11/o1_ledger.md:12790)


## §6fi — L3: single-candidate ties under the pairwise rule (PREREG `docs/archive/2026-09/PREREG_loose_end_L3_2026-09-05.md`, md5 85ee0147) (2026-09-05)
L3 prereg: `RUSTLE_STAR_P=0` on NPIP/LCR16u. P1 held (NPIP anchors 15/15/0 wrong), P3 held (786/787, 915/915); P2 FAILED (row 716: dropped EIF3C member 77 → 433, not ≥ 4,000; 4,838 single-candidate ties remain); P4 FAILED (NPIP 97 s → 1,570 s). **Not adopted**; sole-candidate label left to user.
(full: archive/2026-11/o1_ledger.md:12839)

## §6fj — L3 label (user decision), L4 soft abundance, L6 the read-star proof dump (2026-09-06)
User decision: a certified sole candidate is `assigned`, with new column 18 `sole_candidate` (row 716). `tied` now = K = 0 pair only (233 of 24,462 unit reads on the paired 35); unexplained = `ambiguous` + `origin_rejected`. Also L4 soft abundance, L6 read-star proof dump.
(full: archive/2026-11/o1_ledger.md:12863)

## §6fk — L5: genome-wide catalog `gw_units_v2` and the genome-wide O2 run (2026-09-06, in progress)
L5: `gw_units_v2` catalog, **3,986 units** (v1 3,756): 3,143 kept_full, 611 kept_trimmed, 232 dropped members emitted; genome-wide O2 run (in progress at writing; harness memory kills, chunked driver).
(full: archive/2026-11/o1_ledger.md:12906)

## §6fl — THE O3 FLAG PASS (PREREG `docs/archive/2026-09/PREREG_o3_flag_pass_2026-09-06.md`, md5 400028a5; `bench/o3_flag_pass.py`) (2026-09-06)
O3 flag pass prereg: origin-rejected reads per best candidate vs certificate-accepted control, Poisson test. **50 flagged `missing_copy` (28.9 %, 25 families)**; P3's 25 % line missed (row 718); control < 1/kb only 456/763, P2 fails genome-wide (row 717). Most "missing" origins = annotated-no-unit or other-family loci, not reference absence.
(full: archive/2026-11/o1_ledger.md:12930)

## §6fm — A LOCUS NEVER CONTAINS ANOTHER CATALOG UNIT, AND A HIT MUST TOUCH THE UNIT (from the 15 suspect controls) (2026-09-06)
Diagnosis (row 719): 13 of 15 suspect controls = loci whose reference haplotype differs ~0.5 % from the reads; rule: a locus never contains another catalog unit, a hit must touch the unit. Reruns: 3 contigs 60 flagged (34.7 %), genome-wide 711 pairs, 99 flagged (13.9 %).
(full: archive/2026-11/o1_ledger.md:12973)

## §6fn — O3 items 2–4: reconstruction over every flag, the annotated-no-unit loci, the cut certificate (2026-09-06)
O3 items 2-4. Reconstruction (row 720): `reference_absent` is the verdict by construction, cannot separate absent copy from other causes. Cut certificate (row 721): 44 nested (4,042 reads, an O1 overlap not a cut), 27 `disjoint_no_homology` (4,179 reads).
(full: archive/2026-11/o1_ledger.md:13006)

## §6fo — THE SD CORE WITHOUT SEDEF: pairs derived from the catalog's own alignments (2026-09-06)
`mcl_families --core-from-paf`: SD core from the catalog's own asm20 alignments instead of SEDEF; self-contained fallback for substrates without an SD call. SEDEF-criterion depth rule on PAF pairs not run.
(full: archive/2026-11/o1_ledger.md:13043)

## §6fp — ONE SENSITIVITY, ONE SPECIFICITY FOR NPIP: reads of known origin, the two-form certificate, and what still abstains (PREREG `docs/archive/2026-09/PREREG_npip_known_origin_2026-09-06.md`, md5 c98df4cc + two amendments) (2026-09-06)
NPIP known-origin prereg (+2 amendments). Simulation: precision held (P1), sensitivity not. Real data P4 FAILED (row 722): two-form certificate changes nothing on real NPIP. The 71 flagged (row 723): 48 real read-throughs. Real contested NPIP: 472 of 1,000 below MAPQ 60, 182 assigned (38.6 %); anchors 5 assigned, 5 right, 0 wrong, 28 abstain.
(full: archive/2026-11/o1_ledger.md:13064)

## §6fq — THE MACHINERY IS FOR THE CONTESTED READS; ONE TABLE OVER EVERY READ (user, 2026-09-06 16:16)
Rule (user): certificate machinery decides contested (MAPQ < 60) molecules; uncontested assigned by certified call else placement. Certified-first justified by sim arm B (row 724). Result: 23,749 (97.1 %) vs 19,281 (78.8 %) assigned, MAPQ-60 agreement 23,248/23,248; all 76 families 48,012 of 49,449 (97.1 %).
(full: archive/2026-11/o1_ledger.md:13110)

## §6fr — TWO METRICS, NOT ONE: O1 over every locus (expressed or not), O2 over the reads the aligner could not place (user, 2026-09-06 16:47)
Two metrics: O1 locus metric over every truth locus (expressed or not; `adj/o1loci/npip_loci.tsv`), O2 over reads the aligner could not place. Posterior answers "most likely copy", certificate answers "may I claim it"; row 724 holds (MAPQ 0: minimap2 18 of 39 right).
(full: archive/2026-11/o1_ledger.md:13146)

## §6fs — BIPARTITE SIZE MATCHING ON THE NEW NODE (the advisor's 2026-07-28 instrument, rerun; user 2026-09-06 17:19)
`bench/soto/bipartite_size_match.py` rerun on new node; size is a measurement not a filter. Chain/model 0.60, locus extent vs gene span 1.02, vs core 1.42; quote the pair. 5 unmatched predictions = dropped members; 2 unmatched truths.
(full: archive/2026-11/o1_ledger.md:13189)

## §6ft — O1 EVALUATED BY SENSITIVITY, SPECIFICITY AND BIPARTITE COVERAGE; EVERY NPIP MEMBER ACCEPTED; PARTNERS AND THE READ-THROUGH CERTIFICATE (user, 2026-09-06 17:50; PREREG `docs/archive/2026-09/PREREG_npip_all_members_2026-09-06.md`, md5 167d233e)
`bench/o1_eval.py` (sensitivity, specificity, bipartite coverage); every NPIP member accepted. Polish 1 (`gff_fallback`) adopted; Polish 2 inclusive majority not adopted (row 725); partners adopted; read-through certificate adopted (escape `--no-readthrough-certificate`). 38 of 72 run-past reads now pass; NPIP-proper assigned 706 of 1,001 (70.5 %); cost 130 → 439 s.
(full: archive/2026-11/o1_ledger.md:13210)

## §6fu — THE ANNOTATION IS A PROPOSER, NOT THE DEFINITION: a degradation ablation and a core projection (user, 2026-09-07; PREREG `docs/archive/2026-09/PREREG_annotation_ablation_2026-09-07.md` md5 65efc5b2 and `docs/archive/2026-09/PREREG_core_projection_2026-09-07.md` md5 8f7849228fba, AMENDED)
Annotation ablation + core projection prereg (amended). A0 P0 HELD (byte-identical), A1 P1 HELD (±5 kb), A2 P2 HELD (±20 kb); survivor cores land on held-out locus 25/25 folds. Annotation is a proposer; only discovery of a family with no annotated member stays annotation-bound (G5). One family, one substrate.
(full: archive/2026-11/o1_ledger.md:13247)

## §6fv — THE SD-CORE DEFINITION ON HUMAN NPIP (user, 2026-09-07; PREREG `docs/archive/2026-09/PREREG_npip_human_2026-09-07.md`, md5 dfbe88b6)
Human NPIP (CHM13 chr16+18, `soto_mcl/npip_hsa`), never pooled with gorilla. P2 held 22/22, P4 held 0.92, P5 held (NPIPA+NPIPB one cluster MCL0); P1 held on best-overlap only. Human specificity 0.963. SEDEF sharpens core (0.89 → 1.00), `--core-from-paf` a fallback not an equivalent.
(full: archive/2026-11/o1_ledger.md:13300)

## §6fw — THE CONJOINED READ-THROUGH IS AN O1 OBJECT (user, 2026-09-07; PREREG `docs/archive/2026-09/PREREG_readthrough_object_2026-09-07.md`, md5 1a51fa3b)
Conjoined read-through as an O1 object (prereg 1a51fa3b). P0 HELD, P3 HELD (26/26, 25/25), P4 HELD (663 molecules assigned to read-through units, 551 previously `origin_rejected`); P5 FAILED as worded (242 move, 240 into a read-through unit). Suite 866/0/11.
(full: archive/2026-11/o1_ledger.md:13380)

## §6fx — Y AMPLICONIC FAMILIES: DAZ (user, 2026-09-07; PREREG `docs/archive/2026-09/PREREG_yags_daz_2026-09-07.md`, md5 f66fca74 after AMENDMENT 1)
DAZ (human chrY, prereg f66fca74). P1 HELD, P2 HELD 97.6 %, P3 FAILED (3.7 % assigned; DAZ 302 = 3.3 % clean), P4 HELD (0 read-through units), P6 HELD for XKRY/BPY/VCY, TSPY/RBMY/CDY/HSFY/PRY fragment. Correction: quote 39.6 % assigned / 28.3 % clean on the binary's copy-overlap definition.
(full: archive/2026-11/o1_ledger.md:13465)

## §6fy — THE FAMILY AS A CLOSURE: a fixed point, measured (user, 2026-09-07; PREREG `docs/archive/2026-09/PREREG_closure_definition_2026-09-07.md`, md5 8f144fd6)
Family as a closure (fixed point, no MCL; `bench/o1_closure.py`). P1 HELD (2/2/1 iterations), P3 HELD (gorilla 25/26, human 26/26); Q2 HELD (precision 1.000 from every start); Q3 holds gorilla, FAILS chrY; R1-R3 HELD (DAZ 4-copy level, shared segment 61.3 kb). Open: `ID_9` missed, human not rerun under AMENDMENT 3.
(full: archive/2026-11/o1_ledger.md:13559)

## §6fz — SOTO'S BENCHMARK IS A SEGMENTAL-DUPLICATION CATALOG, NOT A GENE-FAMILY ONE (user, 2026-09-07)
User claim that Soto's families are 98 % SD catalogs (pseudogenes, fragments, whole genes conflated) checked against CHM13 RefSeq biotypes: **the claim holds**. Exon-carrying lever unmeasured; 110 of 168 pseudogene members are `transcribed_pseudogene`.
(full: archive/2026-11/o1_ledger.md:13679)

## §6ga — THE SIZE REFERENCE MOVED TO THE ANNOTATION; THE EXONIC CORE IS NOT THE CORRECTIVE (user, 2026-09-07; PREREG `docs/archive/2026-09/PREREG_size_reference_2026-09-07.md`, md5 e5c01339)
Part A (adopted): size reference = annotation gene span; A1 HELD (0.75 → 0.82), A3 HELD (over-extended tail 37 → 1). Part B exonic core: B2 FAILS (protein-coding drop rate 0.067 → 0.117). Codon-level frame rule proposed, not built.
(full: archive/2026-11/o1_ledger.md:13741)

## §6gb — `--coding-core`: the frame test flips the class bias and DESTROYS NPIP (user, 2026-09-07)
`mcl_families --coding-core` shipped OFF (longest ORF ≥ half the family's best): flips the class bias and destroys NPIP. Cross-species data (`PTR_mm.bam`, `PPY_mm.bam`, `A119b.t2t.bam`) named as the open route.
(full: archive/2026-11/o1_ledger.md:13818)

## §6gc — THE READ-THROUGH JUNCTIONS REPLICATE IN CHIMPANZEE (user, 2026-09-07; PREREG `docs/archive/2026-09/PREREG_crossspecies_junctions_2026-09-07.md`, md5 b0794d03)
Read-through junctions in chimpanzee (prereg b0794d03): P2 HELD (negative 0.00), P3 HELD (read-throughs replicate 0.71), P1 FAILED (positive control 0.39 not ≥ 0.60). 12 non-replicating not shown to be artefacts.
(full: archive/2026-11/o1_ledger.md:13974)

## §6gd — NO minimap2 SETTING IS THE RIGHT PLACE FOR THIS: the aligner line closes (user, 2026-09-08; PREREG `docs/archive/2026-09/PREREG_minimap_readthrough_arms_2026-09-07.md`, md5 208648a5 + AMENDMENT 1)
minimap2 arms (prereg 208648a5 + AMENDMENT 1): B0 M1 HELD (8,719); M2 HELD (disabling long join destroys 20 of 31 conserved junctions); M3 HELD (not selective). Aligner line closes; only evidence for reality = §6gc (30 of 42, negative 0/42).
(full: archive/2026-11/o1_ledger.md:14007)

## §6ge — ⚠⚠⚠ MOST "READ-THROUGHS" ARE ORDINARY INTRONS OF ONE ANNOTATED GENE THAT OUR OWN UNITS SPLIT (user, 2026-09-08)
Most "read-throughs" are ordinary introns of one annotated gene split by our own overlapping units. Fix upstream: forbid or fold cross-family unit overlap.
(full: archive/2026-11/o1_ledger.md:14063)

## §6gf — THE TWO RULES THAT REMOVE THE MANUFACTURED READ-THROUGHS (user, 2026-09-08)
Two rules, both OFF by default, both needed (cross-family exon overlap + own-donor); read-through set shrinks 42 → 14. Suite 867/0/11.
(full: archive/2026-11/o1_ledger.md:14101)

## §6gg — THE 14 SURVIVORS, RE-TESTED IN TWO OUTGROUPS (user, 2026-09-08)
14 survivors re-tested in two outgroups: seven replicate in an outgroup, five in both, scrambled control in neither (annotated introns 32 → 1). Not "read-throughs are real biology" as a class.
(full: archive/2026-11/o1_ledger.md:14137)

## §6gh — THE ADVISOR'S FIVE COMMENTS: answers, and the framing they settle (user, 2026-09-08)
Advisor's five comments answered; revised pseudocode `docs/archive/2026-09/PSEUDOCODE_2026-09-08.md`. Deliverable = per-copy isoform GTF; flair/StringTie substitution test drops to secondary.
(full: archive/2026-11/o1_ledger.md:14166)

## §6gi — THE ISOFORM→COPY JOIN, BUILT: the sentence the advisor keeps asking for (user, 2026-09-08)
Isoform→copy join built (`--gtf`); a live emitter defect found on first run. Per-copy table now has position vs evidence columns (copy 49: 44 by position, 0 by evidence): overlapping a locus is not evidence of belonging.
(full: archive/2026-11/o1_ledger.md:14217)

## §6gj — StringTie ON THE SAME READS, INTERSECTED (user, 2026-09-08)
StringTie `-L` on the same gorilla BAM (7,223 transcripts): agree on 50 chains; ours explains 5 points more molecules exactly, StringTie recovers 16 points more junctions; only ours assigns a copy.
(full: archive/2026-11/o1_ledger.md:14282)

## §6gk — WHERE StringTie IS FOOLED ON THE Y AMPLICONIC FAMILIES (user, 2026-09-08)
StringTie on human chrY YAGs (930 transcripts). A false double-counting headline came from the empty-chain bug (register 757 reintroduced): true figures 17 % MAPQ-0, zero shared molecules. Any read-to-transcript matching must special-case the empty chain.
(full: archive/2026-11/o1_ledger.md:14358)

## §6gl — THE ISOFORM GTF NOW CARRIES ITS COPY, ITS EVIDENCE AND ITS STATE (user, 2026-09-08)
`--gtf` emitter now carries copy tag, evidence and state; position and evidence agree on 119 of 120 tagged isoforms. Suite 867/0/11.
(full: archive/2026-11/o1_ledger.md:14393)

## §6gm — THE DEFERRED TWO: a duplicate-id bug, and evidence for the rule defaults (user, 2026-09-08)
"Isoforms matching no read" = duplicate-id bug (`DenovoTranscript::tid` collides). Evidence for flipping the read-through rule defaults gathered; flip must be pre-registered (moves 60 unit rows across 8 families).
(full: archive/2026-11/o1_ledger.md:14422)

## §6gn — ⛔ THE DEFAULT FLIP IS REVERTED: the untested substrate failed (user, 2026-09-08; PREREG `docs/archive/2026-09/PREREG_cross_family_overlap_default_2026-09-08.md`, md5 92c2c1e0)
Cross-family exon-overlap default flip REVERTED (prereg 92c2c1e0): P1, P2, P3 (0.909/0.663), P5 (14, of which 1) HELD; P4 FAILED (human NPIP 24/26 and 0.960). **Stays an opt-in flag.**
(full: archive/2026-11/o1_ledger.md:14454)

## §6go — WOULD `flair correct` HAVE AVOIDED THE READ-THROUGHS? No — there is nothing to correct (user, 2026-09-08)
`flair correct` would not avoid the read-throughs: nothing to correct (14 survivors' ends sit at annotated-like sites); annotation-based correction rejected as in register 729; short-read route unavailable.
(full: archive/2026-11/o1_ledger.md:14481)

## §6gp — WHAT ELSE TO TAKE FROM flair: PRODUCTIVITY, PER COPY (user, 2026-09-08)
From flair take predictProductivity per copy; only 121 of 886 isoforms carry an assigned copy (22 copies). Next: `orf_aa`/`productive` in `--gtf`. Not built.
(full: archive/2026-11/o1_ledger.md:14510)

## §6gq — `copy_assign --productivity`: per-copy coding potential, shipped OFF (user, 2026-09-08)
`copy_assign --productivity` shipped **OFF**: `orf_aa` per transcript from genomic exon-sum; small denominators caveat. Suite 867/0/11.
(full: archive/2026-11/o1_ledger.md:14547)

## §6gr — "IS THIS AN ASSEMBLER?" — yes in shape, no in kind, and the difference is measurable (user, 2026-09-08)
"Is this an assembler?": yes in shape, no in kind; copy-of-origin is undefined for every existing assembler. Concession: StringTie recovers more read-supported junctions (58.6 % vs 42.1 %, §6gj).
(full: archive/2026-11/o1_ledger.md:14577)

## §6gs — What flair does with secondary alignments (2026-09-08, source read, not run)
flair 3.0.1 source read (not run) on secondary alignments: families > ~5 invisible (F2), alphabetical tie-break (F3), paralog mode inflates by k (F4), BAM overstates certainty (F5). Not a benchmark.
(full: archive/2026-11/o1_ledger.md:14609)

## §6gt — THE TOOL BAKEOFF: isoseq collapse vs ours on HUMAN NPIP (2026-09-08)
Tool bakeoff, human NPIP only (prereg 428a588b, amend acbe7b9c): isoseq never places 23 % of molecules, ambiguous ones disproportionately (30.2 % vs 16.3 %). "Ours is better" NOT supported: tie on placement and copy coverage; separation = conflation (57 vs 0) and hard molecules.
(full: archive/2026-11/o1_ledger.md:14676)

## §6gu — RUNTIME OF THE PAIRWISE DNA ALIGNMENT, AND JACCARD AS A REPLACEMENT (advisor Q3, 2026-09-08)
Advisor Q3: pairwise DNA alignment runtime vs Jaccard. 84 s for a panel; Jaccard prefilter saves nothing worth having since the core rule needs the alignment coordinates anyway.
(full: archive/2026-11/o1_ledger.md:14760)

## §6gv — O2's SCOPE CORRECTED: only AS-TIED MULTIMAPPERS (user, 2026-09-09; PREREG `docs/archive/2026-09/PREREG_as_tied_only_2026-09-09.md`, md5 `268fe69b`)
O2 scope = AS-tied multimappers only (prereg 268fe69b). Sole-candidate rulings conflict; table reported with sole candidates held OUT. Region-local caveat; held-back `fam_MCL2_073244` must be scored before restating (trap 15). Numbers superseded by §6hf/§6hk.
(full: archive/2026-11/o1_ledger.md:14813)

## §6gw — KNOWN SUBCLUSTERS: NPIP AND TBC1D3 vs THE PUBLISHED STRUCTURE (advisor, 2026-09-09; PREREG `docs/archive/2026-09/PREREG_known_subclusters_2026-09-09.md`, md5 `a443b12a`)
Known subclusters (prereg a443b12a): TBC1D3 not split = correct (P2, P3, P4 held); NPIP inflation never finds A/B (P1 FAILED, replicating §6ew). Gap test: an earlier p superseded by p = 0.072; P4 FAILED human NPIP Silverman p = 0.125. Nothing proposed as default.
(full: archive/2026-11/o1_ledger.md:14908)

## §6gx — THE GRAPH THAT REFLECTS WHAT O1 SHIPS: partition test on KEPT members only (user, 2026-09-09)
Partition test on KEPT members only (`gw_subfamily_scan.py --units`): 2 over-merge candidates (MCL38, MCL99), 2 borderline, 9 certified subfamily splits (biology). 67 families untestable at `MIN_COMP = 3`. Gorilla only.
(full: archive/2026-11/o1_ledger.md:15192)

## §6gy — O2 REFACTORED: the AS-tied gate runs BEFORE the certificate, default ON (user, 2026-09-09)
O2 refactor: AS-tied gate before the certificate, default ON. 4958 human `assigned` rows are tied with ONE candidate; L3 ruling vs 09-09 definition conflict.
(full: archive/2026-11/o1_ledger.md:15217)

## §6gz — ⛔⛔ THE SOLE-CANDIDATE CLASS WAS A CROSS-FAMILY LEAK: 4,706 EIF3C/EIF3CL reads "assigned" to an NPIP copy (user, 2026-09-09)
Sole-candidate class was a cross-family leak: 4,706 EIF3C/EIF3CL reads assigned to an NPIP copy. Tie-outside forced `Tied`, never `Assigned`; `tie_invariant` retired (row 786); held-back MCL2 still unscored.
(full: archive/2026-11/o1_ledger.md:15281)

## §6hb — REAL flair AND StringTie RUNS ARRIVE (user, 2026-09-09): §6gs's source-read predictions confirmed empirically, and a proper 4-way family-region bakeoff
Real flair 3.0.1 and StringTie 3.0.1 runs + new isoseq arm: §6gs source-read predictions confirmed empirically; 4-way family-region bakeoff with gffcompare cross-check. Framing unchanged (hard loci, copy attribute).
(full: archive/2026-11/o1_ledger.md:15356)

## §6hc — `--origin-drop-indels` + the excision control (2026-09-09; full record in `docs/archive/2026-09/PREREG_origin_drop_indels_2026-09-09.md` d4e85007 and `docs/archive/2026-09/PREREG_excision_odi_2026-09-09.md` f8b68981)
`--origin-drop-indels` + excision control (preregs d4e85007, f8b68981): rescues PSV-resolved assignments blocked by real indels; copy 2 excised 149/153 = 97.4 % abstain, copy 22 42/42. Default left to user (adopted ON in §6hf).
(full: archive/2026-11/o1_ledger.md:15426)

## §6hd — ⛔ `--admit-aligner-disagreement` REFUTED as a gate widening (2026-09-09; PREREG 58a5496b, outcome appended)
`--admit-aligner-disagreement` REFUTED (prereg 58a5496b): 2,278 human AS-unique≠primary molecules; 3× assignments at a quarter of the evidence, 31 % `other`; not recommended.
(full: archive/2026-11/o1_ledger.md:15437)

## §6he — ⛔ INDEL PSV COLUMNS ADD NOTHING; two defects exposed (2026-09-09; PREREG 021446fb, outcome appended; user's step-2 hold: "PSVs can also be indels")
Indel PSV columns add nothing (prereg 021446fb); two defects exposed (`psv_score` not pairwise). Register rows 793-794. Not applied.
(full: archive/2026-11/o1_ledger.md:15450)

## §6hf — ⭐⭐⭐ O2 DEFAULTS DECIDED (user, 2026-09-09, step 2) AND THE HEADLINE NUMBERS RE-BASED (step 3)
O2 defaults decided (user): `origin_drop_indels` ON, `admit_aligner_disagreement` OFF, `indel_psv` OFF. GTF vs O2 copy agree 140/153 = 0.915. Supersedes §6gv "759 / 69", §6gz "14 of 759", §6hb. Rows 791-794 stand. Itself superseded by §6hk.
(full: archive/2026-11/o1_ledger.md:15467)

## §6hg — THE EXCISION SWEEP: every human copy removed in turn (2026-09-09; PREREG 4c832450, outcome appended)
Excision sweep (prereg 4c832450): 219/230 = 95.2 % abstain when their copy is removed; 4 one-column assignments (10↔11) flip; 19/230 rest on ≤ 2 columns. Rows 795-796.
(full: archive/2026-11/o1_ledger.md:15500)

## §6hh — THE HARD-LOCUS BAKEOFF (advisor's metric, Part 0h; 2026-09-09; PREREG 5ca5c7e4, outcome appended)
Hard-locus bakeoff (prereg 5ca5c7e4): at hard loci we carry 2-3× flair/StringTie and, support-matched, 1.3× isoseq; unmatched, isoseq singletons out-carry us. Easy set is an unmatched control (row 797); singleton transcripts (row 798); copy attribution 140 of 153 at O2's copy (92 %).
(full: archive/2026-11/o1_ledger.md:15513)

## §6hi — THE HELD-BACK FAMILY: gorilla MCL58 (2026-09-09; PREREG 576f7f83, outcome appended)
Held-back gorilla MCL58 (prereg 576f7f83): 276 AS-tied, contested 40 = 2 assigned / 38 tied / 0 ambiguous; both assignments abstain when copy 2 excised. Refactored O2 transfers.
(full: archive/2026-11/o1_ledger.md:15535)

## §6hj — `--best-by-duel`: the read-star best by pairwise MAXIMIN (2026-09-09; PREREG 819c1615, outcome appended)
`--best-by-duel` (pairwise maximin bk, prereg 819c1615) behind a flag; 28 ambiguous → origin-rejected; certificate is bk-only (row 799). Default left to user (adopted ON in §6hk).
(full: archive/2026-11/o1_ledger.md:15543)

## §6hk — ⭐⭐⭐ `best_by_duel` DEFAULT ON (user, 2026-09-09) — THE SHIPPED O2 PIPELINE AND ITS NUMBERS, FINAL FOR THE DAY
`best_by_duel` default ON (user). Shipped O2: AS-tied gate → tie-outside never assigned → pairwise maximin. Own-state 263 assigned / 5,793 tied / 1,601 ambiguous over 7,657; GTF vs O2 166/182 = 0.912. Supersedes §6hf's 230 / 531 / 382 and "14.0 %".
(full: archive/2026-11/o1_ledger.md:15556)

## §6hl — ISOFORM PHANTOMS IN THE SHIPPED GTF, measured (2026-09-09; PREREG 62f34c64, outcome appended)
Isoform phantoms (prereg 62f34c64): 143/542 multi-intron isoforms lift to other copies; 53 of 542 in twin groups. Register row 800.
(full: archive/2026-11/o1_ledger.md:15578)

## §6hm — ⛔ THE POOLED (ISOFORM-LEVEL) CERTIFICATE IS NOT A PLACER (2026-09-09; PREREG 12576433, outcome appended)
Pooled isoform-level certificate is NOT a placer (prereg 12576433): only 27 abstaining reads would gain; origin certificate blind below ≈ 0.6 %. Row 801. Place by read-level evidence only, emit copy SETs.
(full: archive/2026-11/o1_ledger.md:15590)

## §6hn — THE GTF O2 BELIEVES: evidence-placed transcripts and copy SETS (2026-09-09; PREREG 95409846, outcome appended)
`bench/gtf_copy_set.py` (prereg 95409846): evidence-placed transcripts + copy sets; 684 → 676 family transcripts; tie set ≠ read-star candidate list (row 802).
(full: archive/2026-11/o1_ledger.md:15602)

## §6ho — THE LIFT-AWARE ISOFORM BAKEOFF, and the per-tool PHANTOM table (2026-09-09; PREREG 84e8ed50, outcome appended)
Lift-aware isoform bakeoff (prereg 84e8ed50), 3,633 hard spliced molecules: certified isoforms at O2's copy ours 100 %, flair 78 %; lift-aware carried rewards duplicates (row 803). Rule not yet in binary.
(full: archive/2026-11/o1_ledger.md:15612)

## §6hp — `--gtf-copy-set` PORTED INTO THE BINARY (2026-09-09 night; default OFF, the user decides)
`--gtf-copy-set` ported into `copy_assign` (default OFF, user decides); off = byte-for-byte (`a4d0f5cd`); assignments untouched (`8a057f68`, `c87c7f71`).
(full: archive/2026-11/o1_ledger.md:15628)

## §6hq — THE EXON-SUM WIDTH IS NOT SYSTEMATICALLY UNDER-REPRESENTED; a self-caught metric trap along the way (2026-09-09; PREREG 64f80dd7, outcome appended)
Exon-sum width not systematically under-represented (prereg 64f80dd7): median relative deficit ≈ 0 (−0.1 % both species); metric trap self-caught (row 805). Fix = containment guard on assembler boundaries, not implemented.
(full: archive/2026-11/o1_ledger.md:15650)

## §6hr — ⛔ `--min-isoform-fraction` DOES NOT SOLVE THE READTHROUGH PROBLEM (2026-09-10; PREREG fd894558, outcome appended)
`--min-isoform-fraction` (prereg fd894558, default off) does not solve the read-through problem; kept as relative-depth report. Rows 807-808.
(full: archive/2026-11/o1_ledger.md:15671)

## §6hs — ⭐⭐⭐ `--min-boundary-fraction` WORKS: the locus-boundary outlier test (2026-09-10; PREREG 35e290a9, outcome appended)
`--min-boundary-fraction` works (prereg 35e290a9; 0.10 with `--min-boundary-gap 5000`). Addendum: default ON at 0.10 (user, 2026-09-10); escape `--min-boundary-fraction 0` (`888794bf`); suite 871/0/11.
(full: archive/2026-11/o1_ledger.md:15688)

## §6ht — ⛔ IDENTITY-WEIGHTED DENSITY REFUTED GENOME-WIDE; `identity_gap.py` IS THE RIGHT INSTRUMENT (2026-09-10; PREREG b4b65ebd, outcome appended)
Identity-weighted density REFUTED genome-wide (prereg b4b65ebd): do not adopt; `identity_gap.py` is the validated answer to criticism 5. Rows 810-811.
(full: archive/2026-11/o1_ledger.md:15711)

## §6hu — CURRENT O1 STATUS: ZERO-CORROBORATION CLUSTERS ARE REAL BIOLOGY NOT A FALSE MERGE; NPIP'S GAP IS FRAGMENTATION NOT OMISSION (2026-09-10)
O1 status: 35/274 (12.8 %) gorilla clusters `corroborated=0.0` are real biology, not false merges; NPIP's gap is fragmentation (4 clusters: MCL1 25, MCL51 3, MCL183 2, MCL22 1) not omission. Rows 812-813.
(full: archive/2026-11/o1_ledger.md:15729)

## §6hv — SOTO'S OWN PUBLISHED FAMILIES AUDITED: 56.5% PSEUDOGENE, 35.9% OF FAMILIES ENTIRELY PSEUDOGENE, 35% CONTAIN A FRAGMENT BUNDLED WITH A FULL COPY (2026-09-10)
Soto's published families audited: 56.5 % pseudogene, 35.9 % of families entirely pseudogene, 35 % contain a fragment bundled with a full copy. Critique fully supported; `soto_families` rebuild still open.
(full: archive/2026-11/o1_ledger.md:15752)

## §6hw — A3 + B2 SHIPPED: sibling-identity reporting + evidence-backed singleton rescue (2026-09-10)
A3 shipped (sibling-identity reporting in `Assignment`) and B2 (evidence-backed singleton rescue, strand by majority vote). Base `--gtf` unaffected.
(full: archive/2026-11/o1_ledger.md:15785)

## §6hx — CORRECTION: B2's rescue mislabelled 11/14 transcripts; fixed + complete per-record provenance shipped (2026-09-10, same day, user-driven)
Correction: B2's rescue mislabelled 11/14 transcripts (grouped off-locus own alignments); fixed, plus per-record provenance (`--read-provenance`). 792/792 tests; flags off byte-identical.
(full: archive/2026-11/o1_ledger.md:15825)

## §6hy — HARD-LOCUS BAKEOFF RE-RUN, human chr16 MCL0 (26 copies, 65,341 AS-tied records): B2 is a strict superset (2026-09-10)
Hard-locus bakeoff re-run, human chr16 MCL0 (26 copies, 65,341 AS-tied records): B2 is a strict superset. Gorilla-side bakeoff and genome-wide timing not done.
(full: archive/2026-11/o1_ledger.md:15853)

## §6hz — QUICK CLOSES: D7, D4, A5, A7 (docs), B6 (`bench/isoform_bakeoff.py --summary`) (2026-09-10)
Quick closes D7, D4, A5, A7 (docs), B6 (`isoform_bakeoff.py --summary`, composite = carried AND address ⊆ own copies, row 803). Phantom rate: isoseq 85.3 %, ours lowest 66.4 %.
(full: archive/2026-11/o1_ledger.md:15886)

## §6i0 — A2 CLOSED: `origin_rejected` was tested against the duel winner only (2026-09-10)
A2 closed: `origin_rejected` now requires failing at every candidate, not just the duel winner (row 799). 39 reads move `origin_rejected → ambiguous` (1257→1218); off byte-identical; 792/792.
(full: archive/2026-11/o1_ledger.md:15921)

## §6i1 — B1 CLOSED: lift failures are now named, not a bare count (2026-09-10)
B1 closed: lift failures named per transcript/copy pair. Remaining open items C1/C2, D1, D5/D6, C6. Committed `d690715`.
(full: archive/2026-11/o1_ledger.md:15937)

## §6i2 — B8 CLOSED: `--gtf` run on the held-out family MCL58, both fixes hold on unseen data (2026-09-10)
B8 closed: `--gtf` on held-out gorilla MCL58 (9 copies); A2/B1/B2/A3/B6 fixes hold on unseen data. No code change.
(full: archive/2026-11/o1_ledger.md:15959)

## §6i3 — A6 CLOSED: the outside tie partner is now named, not a bareword (2026-09-10)
A6 closed: `--name-outside-tie` (default off) names the outside partner (row 790); merged spans ~29 on chr16, recurring EIF3C/EIF3CL cluster. 792/792; off byte-identical.
(full: archive/2026-11/o1_ledger.md:15975)

## §6i4 — B4 CLOSED: `isoform_copy_lift.py` now recognizes lifted placements as evidence (2026-09-10)
B4 closed: `isoform_copy_lift.py` reads `placed_by`; human MCL0 phantom groups 2 → 1, P5 = 0 (pass condition).
(full: archive/2026-11/o1_ledger.md:15990)


## §6i5 — A1 PARTIALLY TRACED: narrowed, not pinned; one of my own hypotheses caught and corrected (2026-09-10)
Traced molecule `5459_8000` through `no6` excision catalog. Own hypothesis (duplicate-locus bug) caught and corrected: harness deletes copy 6 and renumbers. Molecule has exact 2-way AS tie (2463), one placement inside the excised copy; drops out downstream of the AS-tie gate, before the certificate stage. Exact dropping line NOT found; left open.
(full: archive/2026-11/o1_ledger.md:16003)

## §6i6 — C6: NEW O1 HOLD-OUT FAMILY CHOSEN — FIRST PICK RETRACTED SAME HOUR, REPLACED (2026-09-10)
MCL2 uninformative (1 AS-tied read) so new O1 hold-out chosen by outcome-blind coordinate rule. First pick (repro MCL4) RETRACTED: MCL ids not stable across catalogs; it = audit46 MCL5. Exclusion must be by COORDINATES. 777 clusters survive; hold-out = `gw_units_v3` MCL166 (n=10, NC_073243.2:23,902,119-24,040,023, md5 `8b757279187bc51eb4dc6fea850361a8`); runners-up MCL185, MCL220. Use ONCE, not for tuning.
(full: archive/2026-11/o1_ledger.md:16029)

## §6i7 — REGISTER 757 RECURRED A SECOND TIME, IN THE SAME SESSION'S OWN B2/PROVENANCE CODE (2026-09-10)
Audit found register 757 recurred in `--rescue-singletons`/`--read-provenance` code (default OFF): unspliced records collided onto one gate-passed tid. Fixed (span containment): `contributed_gate_passed` 37,976 -> 32,646 (-5,330); rescue output made deterministic (BTreeMap); dead code removed. Default `--gtf` byte-identical; 792/792 lib tests.
(full: archive/2026-11/o1_ledger.md:16079)

## §6i8 — FOUR MORE AUDIT ITEMS CLOSED (2026-09-10, same evening)
Four audit items: `isoform_bakeoff.py` outside-token crash fixed; A6 outside-locus coordinate off-by-one fixed (`s0 + 1`); `DenovoConfig::from_env` doc fixed (6 env vars; `tied_seed` CLI-only, row 980); `family_define --rna-oracle` is a documented no-op, not a bug. No default behaviour change; 792/792.
(full: archive/2026-11/o1_ledger.md:16133)

## §6i9 — REMAINING FOUR: OVERLAPPING-REGIONS GUARD SHIPPED, TWO DOC-CURRENCY FIXES, ONE CLARIFICATION (2026-09-10, later same session)
`validate_no_overlapping_regions()` guard shipped (rejects overlapping `--regions`); register row 1059 stale line citations fixed; growth-extent PREREG md5 addendum (`8d5a63fb70875996d599c9f4218cd83c`; implemented, off, NOT measured). 792/792 is `--lib` only; 871/0/11 is full suite. Rows 980/1059 not bugs.
(full: archive/2026-11/o1_ledger.md:16164)

## §6ia — TASK 7: `--flag-missing-copies` REPRODUCTION GATE, 3 REAL BUGS FOUND AND FIXED (2026-09-11)
Rust `--flag-missing-copies` vs `bench/o3_flag_pass.py` on 76 families: first 272 vs 173 pairs tested. Port bugs fixed: missing truth gate, block-wise overlap, L2 locus-extent window, name-sorted 500 cap. After: flag matches 185/192 (96.4%), class 190/192 (99.0%); residual = O2 input drift. Citation error later corrected in §6ib.
(full: archive/2026-11/o1_ledger.md:16203)

## §6ib — TASK 7 CLOSEOUT: root-causing the 4 flag-flip pairs — 3 explained, 1 exposes a NEW unfixed defect (2026-09-11)
4 flag flips root-caused: all explained (1 by input drift, 1 by new defect, 2 compound). NEW unfixed defect: `copy_assign.rs:2525` control pool filters `status == Assigned` (no Python counterpart); 116 of 235 groups affected, 8 forced Untestable. Not fixed here (no re-run in scope).
(full: archive/2026-11/o1_ledger.md:16278)

## §6ic — TASK 7 FINAL CLOSEOUT: the §6ib defect fixed, tested, and re-verified against Python (2026-09-11)
§6ib defect fixed (status clause dropped), tests added (lib 816/0/11, bin 23). Re-run: flag matches 185/192 -> 190/192 (98.96%); class 190/192 unchanged; p diff 38 -> 26. Remaining 2 mismatches (`MCL1_073242:11`, `MCL6_073244:2`) = input drift. Task 7 closed.
(full: archive/2026-11/o1_ledger.md:16382)

## §6id — TASK 8: CROSS-INDIVIDUAL O3 DIFFERENTIAL — RUN, 2/2268 candidates, 90.9% no_signal, EXPLORATORY (2026-09-11)
Cross-individual O3 differential, testis (OR6737) vs fibroblast (KB3781): two axes differ. Original "2 candidates, 90.9% no_signal" RETRACTED (SUPERSEDED in continuation): `none` conflated untested with clean. Corrected final: 0 candidate differentials, 2 inconclusive (`MCL2209:0`, `MCL2288:0`), 2072 both-arms-testable, 0.000 no_signal; sanity gate VIOLATED. Round-2 fix: headline unchanged. EXPLORATORY; no statistical power.
(full: archive/2026-11/o1_ledger.md:16474)

## §6ie — SOTO REPLICATION RESUMED: fresh unmerged CHM13 v2.0 SEDEF input, 2 real bugs found+fixed, ARI 0.751 / exact 47.7% (2026-09-11)
Soto replication with unmerged CHM13 v2.0 SEDEF pairs (CIGAR projection, no minimap2). Bug 1: over-broad gene universe gave 149-gene mega-family; fixed by Soto's 1,793-gene set. Bug 2: mis-ported self-mapping guard (ARI 0.084); removed. Result: ARI 0.751 (vs 0.757), exact 234/491 (47.7%, vs 31.6%), pair P/R 0.925/0.634.
(full: archive/2026-11/o1_ledger.md:16677)

## §6if — SOTO'S ACTUAL RELEASED CODE READ DIRECTLY: real MAD-statistic mismatch found and fixed, new best ARI 0.767 (2026-09-11)
Read Soto's released code (`mydennislab/HSD_brain_evolution`). Finding 1: MAD is MEDIAN absolute deviation (scipy), not mean; median-MAD gives ARI 0.7665 (was 0.7510), F1 0.768, but exact 42.4% (was 47.7%); `--mad-statistic` default kept `mean`. Finding 2: their algorithm discards high-dispersion clusters then merges via coding genes; port `soto_cluster_dennislab_algorithm.py` WORSE (ARI 0.541/0.564); parked. Exact reproduction not achievable (manual curation).
(full: archive/2026-11/o1_ledger.md:16788)

## §6ig — Attempting to mimic Soto's manual "redundant gene" curation step (2026-09-11)
Soto step-2 manual redundant-gene curation. Our 1864-gene set = Soto's 1,793 + exactly 71 extra (the 71 removed). Their bedtools detection rediscovers 71/71 among 572 flagged genes; mechanical selection rule REFUTED (67.6% recall, 9.7% precision). No change to reported ARI numbers.
(full: archive/2026-11/o1_ledger.md:16886)

## §6ih — including pseudogenes/lncRNAs as family members: a real improvement, AND a self-caught scoring flaw (2026-09-11)
541 non-eligible member genes (lncRNA/pseudogenes) were absent from all runs. Naive add to graph REJECTED (precision 0.925->0.730). Attach-to-formed-family adopted (275/541 attached). Self-caught shrinking-denominator trap: §6if 0.7665 was 1,793-only; full 2,185-gene universe no-attach ARI 0.6394/0.6190. New headline median-MAD + attach ARI 0.6820 (mean-MAD+attach 0.6716). §6if re-scoped, not retracted.
(full: archive/2026-11/o1_ledger.md:16948)

## §6ii — `--full-geneset` lever added to `soto_cluster_from_shared.py` (2026-09-11)
Opt-in `--full-geneset` flag folds attachment into `soto_cluster_from_shared.py`. Bug caught pre-ship (non-eligible genes entered backbone: 434 families vs 390); fixed by eligible-only edge filter. Reproduces §6ih exactly: 390 families, 1,817 placed, ARI 0.6820, P/R/F1 0.837/0.577/0.683.
(full: archive/2026-11/o1_ledger.md:17022)

## §6ij — bipartite-matched precision/sensitivity against Soto's families (2026-09-11)
Hungarian bipartite precision/sensitivity vs Soto (scoring only; standing no-bipartite rule unaffected). 2,185 universe: median+attach MICRO 0.768/0.630, undetected 154/491; mean+attach MICRO 0.799/0.643, MACRO 0.658/0.634, undetected 143/491 (29.1%). Metric-dependent split vs ARI; report both. Largest missed: ID_176, ID_62, ID_175, ID_1, ID_328.
(full: archive/2026-11/o1_ledger.md:17054)

## §6ik — root-causing the 5 largest fully-missed Soto families (2026-09-11)
5 largest missed Soto families, three root causes: ID_176/ID_62/ID_1 `no_eligible_seed_isolated` (accepted cost of §6ih trade-off); ID_175 `liftover_regime_switch_drop` (fixable); ID_328 `sedef_coverage_gap_below_98pct` (gene-body identity ~0.913-0.916). No single fix.
(full: archive/2026-11/o1_ledger.md:17097)

## §6il — acrocentric liftover fix for ID_175, implemented and measured (2026-09-11)
`--extra-anchors` for acrocentric liftover. Narrow anchors ineffective (GAP_GUARD 50 kb); wide 60-80 kb anchors work. ID_175 recovered (`SEDEFFAM72`). Median+attach ARI 0.6820 -> 0.6827, undetected 154 -> 153; mean 0.6716 -> 0.6723, 143 -> 142. Noise-level otherwise, zero regressions.
(full: archive/2026-11/o1_ledger.md:17160)

## §6im — single-eligible-seed islands: a real, sizeable win (2026-09-11)
Single-eligible-seed islands: component with exactly one eligible gene cannot bridge, so admit wholesale (under `--full-geneset`). 53 components / 155 genes. Median: ARI 0.6827 -> 0.6945, exact 236/491 (48.1%), undetected 103/491; mean: 0.6849, 259/491 (52.7%), 92/491. Precision held (0.837->0.840). ID_62 not recovered (210-gene promiscuous component; corrects §6ik). Lesson for O1 flagged, not checked.
(full: archive/2026-11/o1_ledger.md:17225)

## §6in — ID_328's real root cause, and a weak-signal lever that DOESN'T ship (2026-09-11)
ID_328 real cause: SEDEF identity is wide-unit average (gene body 99.99%/99.94%) plus liftover-guard straddle; 5 more anchors added; standard pipeline unchanged (2017 lifted, 4153 edges). Opt-in weak-edge + own-famCN merge lever: 4/7 components pure false merges, ARI 0.6945->0.6962. NOT adopted. Principle: weak signal needs structural no-bridge guarantee.
(full: archive/2026-11/o1_ledger.md:17283)

## §6io — "what if we tried with the topped-up dataset?" (2026-09-11)
271-sample WSSD famCN via pyBigWig (new verified infrastructure). More samples do not fix it: ID_328/ID_321 CN ranges still overlap; false-merge rate ~50-57% and unstable. §6in ruling stands: not adopted.
(full: archive/2026-11/o1_ledger.md:17352)

## §6ip — "was the file I gave you too conservative?" (2026-09-11)
Unfiltered SEDEF `final_human.bed` (88,756 rows; >=0.98: 4,061 vs 2,156). Modest gain: median+attach+islands ARI 0.6959, 241/491 (49.1%), undetected 99/491; mean 0.6862, 264/491 (53.8%). Liftover drop 6.8% -> 45% (chrY/chrX unanchored; acrocentric arms 44-66%). New best headline.
(full: archive/2026-11/o1_ledger.md:17394)

## §6iq — attempted a systematic dense acrocentric anchor set; found a bigger wall instead (2026-09-11)
Dense acrocentric anchor attempt: not an anchor-density problem; v2.0 tandem arrays many-to-one onto small v1.0 windows on chr13/14/15/21/22. Genuine methodological ceiling; acrocentric-liftover thread closed.
(full: archive/2026-11/o1_ledger.md:17439)

## §6ir — the official, published v1.0 SD dataset (Vollger et al. 2021), replicated natively — worse, disclosed why (2026-09-11)
Vollger 2021 native v1.0 SD calls (5,268 rows >=98%, no CIGAR; linear-fraction projection). WORSE: median ARI 0.6316, 147/491 (29.9%), undetected 168/491; mean 0.6375. Precision drop 0.841->0.760. Union also worse (0.6549/0.6664). `final_human.bed` remains best.
(full: archive/2026-11/o1_ledger.md:17484)

## §6is — "can we somehow get the CIGAR info?" — yes, by deriving it ourselves (2026-09-11)
Derived real CIGARs for Vollger v1.0 pairs via mappy realignment (5,144/5,144, 0 mismatches). Improves over linear-fraction: ARI 0.6618/0.6612, exact 33.6%/37.7%, median P 0.815; union 0.6736/0.6769. Still below `final_human.bed` (0.6959/0.6862), which stays headline.
(full: archive/2026-11/o1_ledger.md:17540)

## §6it — checked Soto's own supplement for a coordinate error; found independent validation instead (2026-09-11)
Soto supplement checked: S1E byte-accurate (0 differences, 1,833 rows); coordinates are not off. S1A per-gene anchors match our derived offsets to the base pair; merged anchors give byte-identical edges (2,223 genes, 4,192 edges), confirming §6iq wall. S1D: hardest misses have 1 SD98 paralog (supports `no_eligible_seed_isolated`).
(full: archive/2026-11/o1_ledger.md:17583)

## §6iu — "do we already have all the info from Zenodo, or should we re-run SEDEF ourselves?" (2026-09-11)
Both Zenodo files merged (5,723 rows >=98%, +8.6%): essentially unchanged (median ARI 0.6618->0.6619; mean 0.6612). Fresh SEDEF re-run not recommended; `final_human.bed` headline stands (0.6959/0.6862).
(full: archive/2026-11/o1_ledger.md:17643)

## §6iv — actually ran LiftOff on our own gene set, to answer "why isn't this like LiftOff" with data (2026-09-11)
Ran LiftOff (`-copies`) v1.0->v2.0 on 2,334 genes: 2,332 lifted; TEKT4P2/ID_328 paralogs lifted individually, none as copies. LiftOff relocates known genes, does not define families; its 68 new loci map onto O3 (unannotated copies), not O1.
(full: archive/2026-11/o1_ledger.md:17670)

## §6iw — adapting LiftOff's primitive into a family-definition signal (2026-09-11)
LiftOff-style all-vs-all own-sequence realignment as family edges (`soto_liftoff_style_families.py`; `::coord` header bug fixed). Raw edge precision 89.3%; ARI 0.6752/0.6485, P 0.891, undetected 160/491. ID_328 split across 3 families. Not adopted as default.
(full: archive/2026-11/o1_ledger.md:17702)

## §6ix — O1's T_CORE has a provable, length-dependent max-divergence formula, validated 20/20 (2026-09-11)
T_CORE=0.13 implies d_max(L) = 1 - exp(-ln(L)/(T_CORE*L)) (derived, not fit). Median 913 bp -> min identity 94.4%; 3700 bp -> 98.3%. Validated 20/20 on SEDEF pairs realigned at gene-body scale (minimap2 proxy; re-validated vs real mechanism in §6j0).
(full: archive/2026-11/o1_ledger.md:17757)

## §6iy — a real minimap2/threshold parameter sweep, for "have we exhausted the parameter space" (2026-09-11)
Sweep of LiftOff-style thresholds (identity plateaus below 0.90) and minimap2 presets: asm20 ARI 0.7063, map-ont 0.7064, asm5 0.6698, ava-pb 0.6959. Higher ARI than `final_human.bed` but bipartite worse (undetected 99 -> 158/491). Not adopted as headline; result not fragile to parameters.
(full: archive/2026-11/o1_ledger.md:17805)

## §6iz — "can we improve de novo mode so it's as similar as possible to annotated mode?" — diagnosis + 5 ranked proposals (2026-09-11)
De novo vs annotated mode gap: ~58% from node construction. Already-killed ideas: rows 1222, 294, 1072. 5 ranked proposals (#1 SD-core projection for discovery, #2 `from_genome` DNA front end default, #3 length-dependent d_max for T_CORE, #4 controlled partitioner swap, #5 opt-in `--seed-from-gff`). User ruling: run 3, then 4.
(full: archive/2026-11/o1_ledger.md:17853)

## §6j0 — Proposal #3 re-validated against the REAL mechanism: agreement holds, but two real problems block shipping it (2026-09-11)
Proposal #3 vs real `confirm_edge` on 31 NPIP oracle loci: formula agrees 417/437 = 95.42%, but real T_CORE admits only 3/437 = 0.69% (genomic-span length dilutes match). New finding: `confirm_edge` severe per-pair cost variance (10 of 13 pairs timed out at 5 s). Ruling: do NOT ship; parked, not refuted. Tests 818/0/13.
(full: archive/2026-11/o1_ledger.md:17925)

## §6j1 — Proposal #4 executed: a controlled partition-algorithm swap, holding the node set fixed (2026-09-11)
Proposal #4: fixed graph, swap MCL vs gamma-quasi-clique. De novo nodes: both 3/31 NPIP; annotation nodes: both 31/31 (7 vs 9 blocks). Node construction CONFIRMED as the bottleneck; MCL leaves more singletons (73.75% vs 100% placed). Absolute numbers not comparable to historical 14/31. Ruling: proposals #1 and #2 worth building; #4 closed. New `mcl_refine` bin, `--dump-graph`.
(full: archive/2026-11/o1_ledger.md:17995)

## §6j2 — Proposal #6 added: Rustle already has an unwired bundle-construction primitive (2026-09-11)
Proposal #6: Rustle has an unwired bundle primitive (`types.rs::Bundle`; `bam.rs:569` `record_ref_span_capped` has zero callers). Use bundle construction as de novo node step; queued as cheaper first step before #1/#2.
(full: archive/2026-11/o1_ledger.md:18082)

## §6j3 — Proposal #6 executed: real evidence redirected the plan mid-flight, real gain found, real new risk found (2026-09-12)
Proposal #6 redirected: giant-intron guard already shipped (`retain_non_mischain`, default on); real failure is undersized/absent nodes. Measured existing `RUSTLE_FOOTPRINT_NODES`: NPIP 7/31 -> 12/31, but one node 16.09x OVERSIZED (not strand-mixing: 22/22 same strand). Ruling: PARKED, not shipped; CBundlenode depth-aware split flagged as next.
(full: archive/2026-11/o1_ledger.md:18115)

## §6j4 — Proposal #1 built: a real, cross-verified span-refinement primitive, with an honestly narrow ceiling on this substrate (2026-09-12)
Proposal #1: `SdPairs::single_span_core` (absolute depth floor min_partners=2, 15 kb slop) on 5 wrong-sized NPIP nodes (26/31 absent): 4/5 improve, 2 near correct ratio (1.07-1.18); Rust matches Python prototype. Ruling: PARTIALLY SHIPPED, not wired; edge effect not measured. Tests 822/0/14. Orphaned `paf_records_become_pairs_in_genomic_coordinates` test noted.
(full: archive/2026-11/o1_ledger.md:18205)

## §6j5 — Proposal #2 built: DNA self-alignment recovers ALL 26 RNA-absent NPIP loci at safe scale; genome-wide window generation hits a real memory wall (2026-09-12)
Proposal #2: Stage A ceiling 31/31 oracle loci linked by DNA self-alignment; Stage B `windows_from_sd_bed` + `--from-genome-sd` with region-restricted 62 windows: 99 reps, 31/31 covered incl. all 26 RNA-absent. 526-window run hit memory wall (>18.5GB), killed. Ruling: PARTIALLY SHIPPED (opt-in, region-targeted). Tests 823/0/16.
(full: archive/2026-11/o1_ledger.md:18286)

## §6j6 — Wiring #1+#2 together: BLOCKED by confirm_edge's own performance wall, reconfirmed independently (2026-09-12)
Wiring #1+#2 through unchanged `detect_edges`/`decompose_families`: full 391-rep run and 42-rep slice (8-min timeout) did not complete. Ruling: BLOCKED, not measured; `confirm_edge` performance wall confirmed independently. Next: per-pair TIME cap with skip-and-report.
(full: archive/2026-11/o1_ledger.md:18348)

## §6j7 — "could `CBundlenode` help the confirm_edge wall?" — profiled it directly: NO, and here's why (2026-09-12)
Profiled `confirm_edge` per pair: no single pathological pair; every pair 2.51s-19.75s (mean 9.26s) on 791-4,684 bp reps (full poasta path). Reps internally non-repetitive; cost is cross-pair similarity. Answer: NO, `CBundlenode` does not help; need per-pair TIME budget. Tests 823/0/18.
(full: archive/2026-11/o1_ledger.md:18409)

## §6j8 — the `confirm_edge` time budget, implemented and verified: real speedup, one new bottleneck found (2026-09-12)
`contiguous_core_coverage_bounded_budgeted` + `DetectParams::time_budget` (default None, byte-identical) with LCS fallback; 500ms chosen. Baseline and #1 cells complete (<100s): oracle 4/31; "false merges" 2/3 and 2/2. Combined cell blocked in `genome_reps`. Ruling then: fix SHIPS (later corrected by §6j9). Tests 827/0/19.
(full: archive/2026-11/o1_ledger.md:18466)

## §6j9 — §6j8's "false merges" adjudicated: 0/9 are false merges; the time budget leaks threads (OOM); the LCS fallback is not faithful (2026-09-12)
False merges adjudicated: 0/9 foreign reps are false merges (3' ends of same NPIP genes; oracle intervals truncate). Time budget leaks threads (OOM 25.7 GB); must not be enabled; §6j8 "SHIPS" wrong. LCS fallback NOT faithful: poasta rejects 8/18 LCS-admitted edges. NPIP problem is FRAGMENTATION. Related register row 813.
(full: archive/2026-11/o1_ledger.md:18568)

## §6ja — Which contiguous-core definition should `confirm_edge` use? Pre-registered test: LCS beats poasta on human (2026-09-12)
Pre-registered (PREREG_core_definition_2026-09-12.md) core-definition test at T_CORE 0.13. Human arm (deciding, 540 scored): POA P 0.927 R 0.547 F1 0.688 AUC 0.796; LCS P 0.914 R 0.813 F1 0.861 AUC 0.948; LCS-masked F1 0.577. Prereg decision: LCS meets the bar. Gorilla arm underpowered; 60/300 DIFF poasta timeouts. Default change needs family-level + held-out first.
(full: archive/2026-11/o1_ledger.md:18648)

## §6jb — LCS edge core: opt-in mode, family level, threshold, held-out (2026-09-13)
`EdgeCore {Poa default, Lcs}` opt-in (`RUSTLE_EDGE_CORE=lcs`). Family level human: LCS ARI 0.681 vs POA 0.525. Step 3 picked LCS >= 50 bp (F1 0.975). Held-out gorilla: flip criterion for 50 bp NOT MET (P 0.921); at natural prevalence 50 bp P 0.459 (register row 814); LCS @ 0.13 P 1.000 R 0.143 vs POA 0.063. 50 bp rejected; LCS @ 0.13 needs one more prereg confirmation.
(full: archive/2026-11/o1_ledger.md:18703)

## §6jc — Pre-registered confirmation: LCS @ 0.13 beats POA @ 0.13 on fresh held-out pairs — criterion MET (2026-09-13)
Addendum C fresh held-out gorilla sample: POA @ 0.13 P 1.000 R 0.060 F1 0.114 AUC 0.694; LCS @ 0.13 P 1.000 R 0.140 F1 0.246 AUC 0.921. Criterion MET; recommend flipping `edge_core` to `Lcs` (user decision). Open: 13% rules recall only ~6-14% of SEDEF paralog pairs.
(full: archive/2026-11/o1_ledger.md:18761)

## §6jd — Default flipped to the LCS edge core; poasta time budget deleted; what changes and what does not (2026-09-13)
User flipped default to `EdgeCore::Lcs` (`RUSTLE_EDGE_CORE=poa` escape); poasta time budget and thread code deleted (commit c97ebd7c). Tests 830/0/19. Affects only de novo `confirm_edge` callers; default homology `gw_family_catalog` (minimap2 E_r) and guided mode unchanged, so gap not yet moved.
(full: archive/2026-11/o1_ledger.md:18792)

## §6je — De novo <-> guided gap: node presence dominates; adding LCS edges to the homology E_r narrows the rest (2026-09-13)
Addendum D: de novo vs guided (`gw_units_v3`) with node set fixed. V3 union (E_r + LCS) NARROWS gap: R_G 0.0140 -> 0.0156, P_G 0.1957 -> 0.2392; V1/V2 do not. Only 3,073/67,844 = 4.5% guided pairs have de novo reps at both loci: node presence dominates (71.6% guided loci lack a rep). Next: V3 opt-in + rebuild + second substrate.
(full: archive/2026-11/o1_ledger.md:18818)

## §6jf — Opt-in LCS union in `homology_blocks`; rebuilt 3-contig catalogs confirm it narrows the de novo <-> guided gap (2026-09-13)
`RUSTLE_ER_UNION_LCS=1` opt-in in `homology_blocks` (commit 50736042; tests 834/0/19). Addendum E rebuilt 3-contig catalogs: default R_G 0.0736 P_G 0.3210 ARI 0.422; union R_G 0.0911 P_G 0.3772 ARI 0.506. Union NARROWS the gap. Not a default yet (NPIP-dense substrate; needs genome-wide + human).
(full: archive/2026-11/o1_ledger.md:18857)

## §6jg — NPIP and TBC1D3: do the modes recover the literature subclusters? Family yes, subfamilies only in pieces (2026-09-13)
Addendum G (descriptive): NPIP/TBC1D3 literature subclusters. Guided holds each family whole (1 family each); de novo NPIP in 5 (default) / 4 (union) families. No mode separates NPIPA/B or TBC1D3 cluster 1/2. Subfamilies only in pieces (NPIP A6+A9, A7+A8, union B4+B5; ARI 0.210 -> 0.300). TBC1D3 group claims CORRECTED by §6l7 (Fig 6C truth: M = B+H, CDKL = K+TBC1D3).
(full: archive/2026-11/o1_ledger.md:18885)

## §6jh — NPIP and TBC1D3 subclusters by gene-structure coverage: no arm separates the major boundaries; guided identity on transcripts comes closest to NPIPA/B (2026-09-13)
Addendum H: gene-structure coverage arms. Level 1 boundary recovered by no arm, no mode (exact 0/2 everywhere); closest guided identity NPIP (3 of 21 wrong side). Hypothesis not supported. Level 2 NPIP guided 5/12 exact, de novo 4/13. TBC1D3 results CORRECTED by §6l7 (by no arm; best 5/7 de novo coverage). Follow-up: TBC1D3 de novo vs guided is not a mode effect.
(full: archive/2026-11/o1_ledger.md:18935)

## §6ji — Missing and fragmented members: RNA locus switches, a readthrough bridge cut, DNA mode and guided on NPIP/TBC1D3; the bridge cut fails the gorilla hold-out (2026-09-13)
Addendum I: TBC1D3 lost to readthrough locus collapse. Opt-in R7 `RUSTLE_LOCUS_BRIDGE_CUT=1` selected (31/31 present, TBC1D3 rep recovered); no arm reduces NPIP fragmentation. DNA mode D1 groups NPIPA vs B via region context but collapses 21/31. Gorilla hold-out FAILS: R_G 0.0734 vs bar > 0.0736; stays opt-in, register row 815. Tests 838/0/19.
(full: archive/2026-11/o1_ledger.md:19018)

## §6jj — RNA fragmentation levers fail the gap hold-out; copy-level DNA nodes (SD atoms) put 4.9x more guided loci in families and are far closer to guided (2026-09-13)
Addendum J. J1: RNA fragmentation levers; F4 `LOCUS_LINK_MIN_READS=10` hold-out FAILS (R_G 0.0729 vs > 0.0736; row 816). J2: SD atoms (`bench/dna_sd_atoms.py`) match guided on NPIP/TBC1D3; NPIPA/B separation does NOT hold at copy level (L1 prec 0.545 < 0.775; row 817). Gorilla hold-out: DNA atoms CLOSER to guided (R_G 0.5585 vs 0.0736, P_G 0.4672). Tests 840/0/19.
(full: archive/2026-11/o1_ledger.md:19086)

## §6jk — Atom edges fixed (ATP5F1A was a retrocopy scored by gap-blind projection); DNA+RNA hybrid and genome-wide gorilla hold-out both pass (2026-09-13)
Addendum K. K0: ATP5F1A retrocopy edge fixed by exact CIGAR two-sided coverage (1.000 -> 0.210), housekeeping 0/3 (post-hoc, rebuild3 becomes development). K1 hybrid (expressed atoms) hold-out NARROWS gap: R_G 0.5174 > 0.1823, P_G 0.6312. K2 genome-wide hold-out: DNA atoms CLOSER to guided (R_G 0.5667 vs 0.0139, P_G 0.5596; 6,529 of 8,295 loci). Housekeeping panel not clean (guided 8/27). No Rust change.
(full: archive/2026-11/o1_ledger.md:19170)


## §6jl — Scope redefined by the user; Addendum L stopped (2026-09-13)
User scope: no genome-only mode. De novo = cluster read-supported loci (fair truth = guided restricted to read-supported loci, R_G 0.1823 / P_G 0.6018); guided = from (minimal) annotation, judged on width + breadth, leave-out simulates minimal annotation. Addendum L genome-wide RD rebuild and L1/L2 stopped, no decision.
(full: archive/2026-11/o1_ledger.md:19239)
## §6jm — Guided mode, first look: leave-out expansion on NPIP and TBC1D3 finds every hidden member at 50%, no over-merge, and single-seed expansion stops at subfamily boundaries (2026-09-13)
Addendum M, descriptive. Keep 50%: sensitivity 1.000 both families; keep 1: NPIP 0.257 ± 0.026, TBC1D3 1.000. Family-named precision 1.000 everywhere (no over-merge). Single NPIP seed reach follows NPIPA/NPIPB lines; width error = NPIP 5' truncation.
(full: archive/2026-11/o1_ledger.md:19256)
## §6jn — Guided rules on NPIP/TBC1D3: iteration crosses NPIPB -> NPIPA from some seeds but not on average; isoform and gene-span width rules trade truncation for overextension; zero over-merges (2026-09-13)
Addendum N. Iteration does NOT cross the subfamily barrier by the bar (keep-1 NPIP I 0.324, I+W2 0.476 <= 0.50); coverage direction is the lever. Width rules W1/W2 lower truncation 0.4-1.0 but raise overextension 1.8-3.6 (bar not met). 0 over-merges.
(full: archive/2026-11/o1_ledger.md:19313)
## §6jo — The literature's two levels on NPIP/TBC1D3: gene-body membership removes the under-merge (keep-1 NPIP 0.257 -> 1.000, no over-merge); one identity cut does not recover NPIPA/NPIPB (2026-09-13)
Addendum O. G1 gene body ADDRESSES under-merge (keep-1 NPIP 1.000, named precision 1.000, 0 cross-family); G1 width does not reduce truncation. G2 single identity cut: NPIPA/NPIPB recovered 0/10 runs (within-B divergence overlaps A vs B). TBC1D3 {H,K} split later corrected by §6l7.
(full: archive/2026-11/o1_ledger.md:19376)
## §6jp — Subfamilies as tree clades: a reference-projected alignment recovers NPIPA|NPIPB and every Dishuck NPIPB group in 10/10 leave-out runs, ~15x faster than MAFFT; hybrid width fixes G1's overextension (2026-09-13)
Addenda P (MAFFT, hybrid width) and Q (reference projection). Q: NPIPA|NPIPB 10/10 + reference; B3-5, B6-9, B12/13 10/10; A6-9 5/10. P2 hybrid width H IMPROVES on G1 (trunc+overext 9.8 vs 13.4 at 50%, 19.8 vs 21.6 at keep 1). TBC1D3 rows corrected by §6l7 (M {B,H} 9/10 projection).
(full: archive/2026-11/o1_ledger.md:19438)
## §6jq — The guided pipeline on the amylase family (CHM13): the unit that holds the family together is family-specific — AMY is held by coding sequence, NPIP by gene body (2026-09-13)
Addendum R. Gene-body membership UNDER-merges AMY at keep 1 (0.400 vs transcript 0.764), opposite of NPIP; no over-merge. H does not improve width on AMY at both levels. AMY1 clade recovered on reference (100/100), AMY2 never. Next: either-finder candidates, exon vs intron trees.
(full: archive/2026-11/o1_ledger.md:19504)
## §6jr — Either-finder candidates pass on all three families; exon and intron trees carry different subfamilies (NPIP introns, AMY exons); LOC124905662 is a partial AMY2A copy with a likely mis-joined first exon (2026-09-13)
Addendum S. Union addresses under-merge on all three families at keep 1, no over-merge; but at keep 50% loses on AMY (v1 U 0.800 vs G1 0.933) from single-linkage clustering. NPIP subfamilies from introns (10/10), AMY types from exons. LOC124905662 exons 2-8 = AMY2A exons 4-10 (AMY2Ap); exon 1 likely mis-join (0 junction reads).
(full: archive/2026-11/o1_ledger.md:19556)
## §6js — Guided fixes consolidated into one tool (`bench/guided_pipeline.py`): members from both finders without duplicates, CDS-envelope gene body, gene-span introns; 2 of 3 families pass every bar, AMY fails one by a mis-joined annotation (2026-09-13)
Addenda T/U (U1''). Addendum T B2 FAIL; final U1'': NPIP/TBC1D3 B1 PASS, B4 PASS; NPIP B2 FAIL by one run (A6-9 9 vs 10); AMY v2 B1 FAIL (0.933 vs 0.967) from LOC124905662's mis-joined envelope; no over-merge. Development families only, not validation. TBC1D3 clade cell RETRACTED by §6l7.
(full: archive/2026-11/o1_ledger.md:19629)
## §6jt — "Did some multi-copy families arise from duplicons?" First pre-registered test: NOT SUPPORTED (6/10 core-duplicon families vs >= 7), classifier discriminates retrocopies (0/3) (2026-09-13)
Addendum W. NOT SUPPORTED (6/10 < 7); retro 0/3. Misses are partial copies inside SDs scored GENE-ONLY. Register row 818. Followed by Addendum X (§6ju).
(full: archive/2026-11/o1_ledger.md:19672)
## §6ju — Genomic copy vs RNA-mediated copy on HELD-OUT families: SUPPORTED (segmental-duplication families 9/9 duplicon-derived, retrotransposition families 0/5) (2026-09-13)
Addendum X: GENOMIC iff >= 1,000 bp shared non-exonic; family DUPLICON-DERIVED iff >= 50% pairs GENOMIC. SUPPORTED (9/9 >= 7, 0/5 <= 1); SD pairs share 39-338 kb, retro <= 243 bp. Development 9/10 (LRRC37A fails: divergence/age limit). Not a core-duplicon origin test; held-out SD families chosen knowingly.
(full: archive/2026-11/o1_ledger.md:19707)
## §6jv — Mechanism layer (genomic copy / retrocopy / TE-derived gene): GENOMIC SUPPORTED (13/13 vs 0/16); RETROCOPY NOT SUPPORTED (5/16; hallmark requirement too strict); TE-DERIVED NOT SUPPORTED (ERV-env 5/8 vs >= 6) (2026-09-14)
Addendum Y (`bench/copy_mechanism.py`): GENOMIC / RETROCOPY (junction loss + poly(A) or TSD) / TE-DERIVED (>= 50% CDS in RepeatMasker TE). Verdicts as headed. Post hoc: intron loss alone separates (166/172 retro vs 0/54 SD); hallmarks found in only 34/172 (poly(A)), 32/172 (TSD). Register rows 819 (R1), 820 (R3).
(full: archive/2026-11/o1_ledger.md:19764)
## §6jw — Retrocopy = parent-intron loss, on FRESH families: NOT SUPPORTED at family level (10/16 vs >= 12) because 6/8 old retrogenes are unassessable; among assessed copies 213/222 RETROCOPY and SD 0/10 families, 0/19 members (2026-09-14)
Addendum Z (hallmarks report-only). NOT SUPPORTED (10/16 < 12; SD 0/10). The call made no error; old retrogenes below the 0.80 / splice-alignability floor need protein alignment. Register row 821. Fixed by §6kc.
(full: archive/2026-11/o1_ledger.md:19818)
## §6jt — De novo: the genomic-span E_r union helps on the development windows and fails the gorilla hold-out's precision guard; exon/intron clades carried over (2026-09-14)
Addendum V; D1 = `RUSTLE_ER_UNION_GENOMIC_SPAN=1` (opt-in). DN0 default byte-identical to §6jg. Hold-out: DN1 does NOT narrow the gap (R_G 0.4116 > 0.3887, P_G 0.2760 < guard 0.2845); 3 h 34 min vs 6 min 24 s. Register row 818 (as written; same row number as the duplicon §6jt). Flag stays opt-in.
(full: archive/2026-11/o1_ledger.md:19858)
## §6jx — Subfamilies as a threshold-free split-dominance tree from variation-graph bubbles: NOT SUPPORTED (literature groups 73 vs IQ-TREE 85), but fewer literature-conflicting splits (49 vs 72); AMY identical to IQ-TREE (2026-09-14)
Addendum AA (`bench/split_tree.py`, four-gamete dominance, no thresholds). Gate PASS. NOT SUPPORTED (73 < 85; conflicts 49 <= 72). Misses = ties/small counts on gap-free columns; NPIP/TBC1D3 conflict index ~0.36, AMY 0.015. §6l7 correction: Fig 6C truth 100 vs 84, still NOT SUPPORTED. Register row 822.
(full: archive/2026-11/o1_ledger.md:19905)
## §6jy — De novo D1 (E_r ∪ genomic-span edges) on the gorilla hold-out: does NOT narrow the gap (P_G falls 0.058 > 0.05 guard), 3h34m vs minutes (2026-09-14)
Addendum V decision: R_G 0.4116 > 0.3887 but P_G 0.2760 < 0.2845 guard. LCS union (§6jf) better on both (0.4710 / 0.3886). D1 stays opt-in. Register row 823.
(full: archive/2026-11/o1_ledger.md:19953)
## §6jz — The shared definition instantiated de novo (read loci + guided edge rule + components): highest de novo recall (R_G 0.509) but FAILS the precision guard (P_G 0.310 < 0.339); gap decomposition: missing nodes 219/741 dominate (2026-09-14)
Addendum AB (`bench/denovo_shared_def.py`). AB2 NOT SUPPORTED: R_G 0.5088 > 0.4710, P_G 0.3100 < 0.3386; loss from transitive components (largest family spans 16 guided clusters). Missing edges 126 -> 60; missing nodes 219 (29.6%) now dominant; class (c) only an upper bound. Register row 824.
(full: archive/2026-11/o1_ledger.md:19973)
## §6ka — The two levers from §6jz: read-locus nodes lift recall to 0.695, leader neighbourhoods lift precision to 0.581, but no pre-registered arm clears both bars (AC3 R_G 0.417 < 0.471) (2026-09-14)
Addendum AC. AC3 NOT SUPPORTED (R_G 0.4170 < 0.4710; P_G 0.5365 passes). Lever N cuts missing-node pairs 219 -> 61; lever G trades recall for precision. Frontier: AC2 0.695/0.300 ... AC1 0.340/0.581; F1 optimum not claimed. Register row 825.
(full: archive/2026-11/o1_ledger.md:20030)
## §6kb — Codon divergence does not explain unreachability uniformly; bridge-split grouping fails precision, but triangle-supported leaders clear the bar as the non-deciding arm (R_G 0.605 / P_G 0.419); a coarser truth explains little of the over-merge (2026-09-14)
Addendum AD. AD-1: synonymous saturation 2/7, protein rescue 3/7, both NOT SUPPORTED (three distinct causes: metric, saturation, diverged protein). AD2a bridge split NOT SUPPORTED (P_G 0.3023 < 0.3386); AD2b triangle (non-deciding) clears: R_G 0.6046 / P_G 0.4191, not adopted pending fresh substrate. AD-3 coarse truth +0.02-0.05 P_G only. Register rows 826, 827.
(full: archive/2026-11/o1_ledger.md:20083)
## §6kc — Retrocopy rule with a coding-sequence, gap-excluded identity and a protein-level fallback: SUPPORTED on fresh families (retro 13/13, SD 0/8); 6 of 7 previously unreachable retrogenes now assessed (2026-09-14)
Addendum AE-B: CDS query, gap-excluded identity >= 0.80, miniprot fallback. SUPPORTED (retro 13/13, SD 0/8). Development retro 15/16 (was 10/16). Of 7 Z-unassessed: 5 RETROCOPY, NAP1L3 UNRESOLVED, UBL4B UNASSESSED.
(full: archive/2026-11/o1_ledger.md:20158)
## §6kd — Triangle-supported leaders CONFIRMED on a fresh gorilla substrate: recall 0.198 vs LCS union 0.120, precision 0.676 vs 0.537 (2026-09-14)
Addendum AE-A (fresh contigs NC_073233.2 / NC_073240.2 / NC_073238.2). CONFIRMED: R_G 0.1983 > 0.1198, P_G 0.6761 >= 0.5370 - 0.05. Grouping clause (§0★★ clause 4) becomes triangle-supported leaders; not yet ported, no default changed.
(full: archive/2026-11/o1_ledger.md:20198)
## §6ke — Human third substrate: triangle-supported leaders CONFIRMED (bipartite F 0.558 vs union 0.515); split read groups (AF-3) NOT SUPPORTED by 0.0005 recall (2026-09-14)
Addendum AF-2/AF-3 (A119b chr15/17/22). AF-2 CONFIRMED: R_G 0.0723 > 0.0670, P_G 0.6822. AF-3 NOT SUPPORTED: R_G 0.0723 -> 0.0718. Bipartite F triangle vs union: 0.682/0.631, 0.691/0.644, 0.558/0.515. Register row 828 (AF-3).
(full: archive/2026-11/o1_ledger.md:20242)
## §6kf — Rust port parity EXACT on both gorilla substrates; the RNA-level truth of clause 6 equals the expressed DNA-level truth once the guided graph is exact; where the remaining de novo error sits (2026-09-14)
AF-1 port (`src/rustle/vg_family/shared_definition.rs`, `RUSTLE_SHARED_DEFINITION=1`): EXACT 76/76 and 72/72 families. Addendum AG RNA truth: expressed loci almost always connected through expressed loci (s1 741 -> 741); "unexpressed bridge" class was an artifact. Triangle s1 error: FN 354, FP 619 (499 from one KZFP family).
(full: archive/2026-11/o1_ledger.md:20289)
## §6kg — AF-1 default byte-identity PASS; the ceiling: annotation-quality nodes give bipartite F 0.955 under the guided construction, de novo nodes 0.726 — node construction is the bottleneck (2026-09-14)
Default path byte-identical to `catF_default` (PASS; full test run pending). `bench/node_graph_mcl.py` oracle nodes F 0.955 vs de novo 0.726; node span Jaccard 0.64. HGNC gene groups rejected as truth (guided MCL F 0.374).
(full: archive/2026-11/o1_ledger.md:20326)
## §6kh — Decomposing the de novo ceiling: exon/UTR completeness costs ~0.09 bipartite F, gene partition ~0.12; end-site gene splitting is the best read-only node rule so far (2026-09-14, development, substrate 1)
Development, post hoc. Read-covered exons lose 0.085 F; read gene partition 0.10-0.12. Edge rule not the lever (guided 0.3 + exon-to-exon best: oracle F 0.873, read components 0.751). Direct truth edges: oracle 557/561, AC 286/561. End-site splitting best read-only node rule.
(full: archive/2026-11/o1_ledger.md:20369)
## §6ju — Simulated IsoSeq: de novo finds all 12 amylase loci in one family at 40 and 10 reads per locus; the real-data misses were expression (2026-09-14)
Addendum W (`bench/amy_sim.py`; note: a second "Addendum W", md5 207f253e). DN0 and DN1 FIND ALL AMYLASE LOCI at both depths although every AMY1-type read is MAPQ 0; exon tree separates AMY1 (100/100) from pancreatic. Real testis misses (DN0 9/12, DN1 11/12) = expression limit.
(full: archive/2026-11/o1_ledger.md:20404)
## §6ki — Guided mode with half the annotation hidden, genome slice: bipartite F 0.844 against the full guided truth; the definition on full annotation 0.968 (2026-09-14, development, substrate 1)
`bench/guided_min.py`, random 50% seeds, DNA truth `gw_units_v3` (1,048 loci). 149 hidden loci never found: 87 in seedless families (unavoidable), 62 genuine misses; slice costs 0.03 F.
(full: archive/2026-11/o1_ledger.md:20442)
## §6kj — Node-rule development closed; held-out Addendum AH launched (2026-09-14)
No development arm adopted; 54/62 discovery misses get no seed alignment (path-connected). Addendum AH (md5 8da9e111; code aa869a08 / 8dd034a8): D (end-split read nodes) vs triangle on gorilla s3 + human chr15/17/22; G (50% annotation) vs seeds-only / all.
(full: archive/2026-11/o1_ledger.md:20465)
## §6kk — Addendum AH held-out: the read-gene-node definition raises recall and bipartite F over triangle on both substrates but fails the precision guard; guided minimal annotation passes on gorilla s3, misses by 0.013 on human; no de novo or minimal-annotation arm reaches the 0.90 goal bars (2026-09-14)
Verdicts as headed. GOAL bars (pair sens/prec, bipartite F >= 0.90) met only by full annotation on human (tautological); gorilla s3 full annotation 0.895 / 0.862. RNA-only nodes capped at F 0.873-0.893 on development. Register rows 829 (D), 830 (G human).
(full: archive/2026-11/o1_ledger.md:20491)
## §6kl — Addenda AI/AJ: the human guided truth is not reproducible from an independent annotation (E0 or E1); E1 (exon-to-exon edges) removes the repeat hub clusters and makes the partition deletion-stable, but fails its pre-registered HGNC guard on the hold-out (2026-09-14)
AI: GENCODE CAT through RefSeq construction NOT annotation-robust (71% false pairs in one 266-gene family; truth built at `--min-exonic-bp 0`). AJ (E1 = `--min-exonic-bp 1`): NOT SUPPORTED on hold-out chr16/19/20 (guard fails; HGNC groups lost were nomenclature/KZFP intronic-repeat edges). E1 partition stable at 50% deletion. Register rows 831 (AI), 832 (AJ).
(full: archive/2026-11/o1_ledger.md:20534)
## §6kp — Addendum X: provenance-typed de novo catalog — restricted DNA admission PASSES the gorilla hold-out (R_G 0.4089, P_G 0.3487 vs DN0 0.3887/0.3345); DNA-only edges are the least precise type (2026-09-14; was headed §6km in bd0d53c4, renumbered after a concurrent §6km)
Addendum X (commit eee52868), opt-in `RUSTLE_ER_UNION_GENOMIC_SPAN=restricted`. DN1r NARROWS the gap: R_G 0.4089 > 0.3887, P_G 0.3487 >= 0.2845. Expectation failed at edge level (rejected span-only 0.085 vs admitted 0.048; TX+DNA edges 0.776). Hold-out DN1r+B NOT RUN (killed). §6l7 correction on CDKL clade (suspended).
(full: archive/2026-11/o1_ledger.md:20617)
## §6km — Post-AJ exploration (looked-at substrates, no verdicts): the ~0.8 ceiling between two annotations is gene-model content, not repeats, MCL settings or partition instability (2026-09-14)
No verdicts. Under E1 definition stable to 50% gene deletion (F 0.97) and edge noise (0.93-0.94); remaining ~0.2 is how annotations model the same genes (exon-union lengths differ 2-3x), so single-annotation truth caps any method at ~0.8.
(full: archive/2026-11/o1_ledger.md:20720)
## §6kn — Adjudicated two-annotation truth (user choice, Addenda AK/AL): not achieved on development; chr5/7/21 hold-out untouched; the ~0.8 cross-annotation ceiling survives every variant (2026-09-14)
Development only. AK gate NOT VALID (evidence sensitivity 0.362 < 0.80; unscored 0.746). AL gate VALID (0.861 / 0.006 / 0.198) but scorer and union-find truth defects; AM matched-loci gate NOT VALID (unscored 0.566). Arbitrated truth is a third MCL partition at ~0.8. No hold-out number. Register rows 835 (AK), 836 (AL/AM).
(full: archive/2026-11/o1_ledger.md:20757)
## §6ko — Protein-space families (Addendum AN): reproducible across independent annotations on a fresh hold-out (GENCODE vs RefSeq 0.935 / 0.999 / F 0.977); how protein, RNA and DNA families relate (2026-09-14)
`bench/protein_families.py` (longest CDS, blastp, cover >= 0.30, MCL 2.8; r1/r2 exclude pseudogenes / V(D)J). DNA E1 comparison 0.617 / 0.529 / F 0.789. Protein and DNA spaces not nested (0.521 / 0.253). Qualification (post hoc, in §6kq): plain protein families are superfamilies; "goal met" is reproducibility only (aa >= 0.50 cross-group 0.4%, sensitivity 0.867). Register row 837.
(full: archive/2026-11/o1_ledger.md:20825)
## §6kq — Spike: pseudohaplotype panels for NPIP/TBC1D3 subfamilies — extra haplotypes make allelic divergence measurable and recover TBC1D3-CDKL; AE still not (2026-09-14, throwaway, not pre-registered) [§6l7 CORRECTED 2026-09-16: TBC1D3 subgroup claims in this section corrected — see §6l7]
Throwaway, not pre-registered. HG002 + 23-haplotype HPRC extension: allelic divergence measurable (0.03-0.06%); NPIP B3-5 and B12/13 one group each across 23 haps. TBC1D3 CDKL / AE claims RETRACTED by §6l7 (CDKL not recovered, M unsupported by last-match, AE untestable). Contains the §6ko superfamily qualification.
(full: archive/2026-11/o1_ledger.md:20893)
## §6kr — Soto 2025 families by name (descriptive): the DNA definition merges co-duplicated SD neighbours (precision 0.02-0.22 on 5/8 families); protein families are precise but miss partial and non-coding copies; only GPR89 and FAM72 are clean in every method (2026-09-14)
Descriptive. Projection: chr1 Soto families nest inside ours, 0.98 cut F 0.742 (universe); chr15/17 neither nests (0.681). Audit: we join genes sharing no exons (TP shared-exon median 0.89 / 0.52 vs FP 0.11 / 0.25); GOLGA8 split by MCL; 103/195 FN lack records. Next: shared-exon rule, -N completeness.
(full: archive/2026-11/o1_ledger.md:20989)
## §6ks — Addendum AP: a fraction-based shared-exon floor (`min_shared_exon_frac`) SUPPORTED on fresh, held-out Soto families — same TBC1D3 family, cleaner GTF2H2/TCAF/GTF2I precision, no recall loss (2026-09-14)
Addendum AP, T = 0.30 fixed on development; hold-out chr5/7/21, 11 Soto families. AP-1 SUPPORTED: bipartite F 0.831 -> 0.881, pairwise precision +0.185, zero recall loss (GTF2I 0.286 -> 1.000). TBC1D3 9 coding copies stay one family. Shipped opt-in; superseded by §6kt (default).
(full: archive/2026-11/o1_ledger.md:21112)
## §6kt — `min_shared_exon_frac` promoted from opt-in to the shipped default (user decision, 2026-09-14)
User decision: `mcl_families --min-shared-exon-frac` CLI default 0.0 -> 0.30 (`GraphParams::default()` stays 0.0). 852 lib tests pass before/after; no-flag run byte-identical to §6ks e1s. Not retroactive: all earlier truths (gorilla `gw_units_v3`, `human2/guided.clusters.tsv`, `lit/**/e1*`) are 0.0; `--min-shared-exon-frac 0.0` reproduces them.
(full: archive/2026-11/o1_ledger.md:21167)
## §6ku — Two detectability-limit spikes: minimap2's identity floor is k=19 exact-seed survival plus chaining, not an arbitrary asm20 choice; the DNA-blind/protein-sighted gap runs ~80% down to ~35-40% identity, driven by dS/dN decoupling, NOT single-codon degeneracy rescue (2026-09-14, descriptive, not pre-registered)
Descriptive. A: asm5/10/20 all k=19; empirical p50 0.786 (theory 0.701); k=12 moves it to 0.709 at ~30,000x spurious hits. B: blastp detects down to ~34-41% nt; gap ~80% to ~35-40%; mechanism dS/dN decoupling (dN/dS 0.10-0.17); single-codon synonymous share median 1.7% (matches AD-1).
(full: archive/2026-11/o1_ledger.md:21195)
## §6kv — Excision extension (PREREG `docs/archive/2026-09/PREREG_excision_extension_2026-09-15.md`): the missing-copy signature does NOT generalize past NPIP/ZNF875 — 5/5 new targets orphan instead of origin-rejecting (2026-09-15)
Fail rule triggers: 5 tandem-SD targets (identity 0.908-0.969), X's reads 99-100% O3 orphan; 0/5 flagged `missing_copy` (0 pairs tested). Signature is contingent on reads staying placeable-but-wrong; "reads vanish" case needs an orphan-count signal. Substrate `sweep_v13`. Register row 838.
(full: archive/2026-11/o1_ledger.md:21282)
## §6kw — The first non-synthetic, real-ground-truth positive: a DNA-confirmed haplotype-CNV family shows the
(heading continues: missing-copy consistent-site signature blind, without excising anything (2026-09-15)). `LOC134757045` silent (0 reads). GWFAM195 (RFPL4B-like; ref 2 units, PAT 4): `copy_assign` AS-tied 33/34, 0 assigned; raw consistent sites fire: copy0 5 sites 2.28/kb at 84-89% (load-bearing), copy1 4 sites 2.01/kb at 56% (vs <1/kb background).
(full: archive/2026-11/o1_ledger.md:21337)
## §6kx — `copy_assign` gains cross-chromosome families (user decision, 2026-09-15); the native O3 Poisson
(heading continues: test runs on both real haplotype-CNV candidates and reports `not_tested` on both, for two different, principled reasons). `copy_assign --families` accepts cross-chrom families (`~xchrom~` sweep key; dedup key gains `chrom`). Known limitation: PSV/mosaic layer has no chromosome field. 854/856 lib tests. `--flag-missing-copies`: GWFAM248 `not_tested` (3 molecules secondary-only), GWFAM195 `not_tested` (33/34 tied).
(full: archive/2026-11/o1_ledger.md:21407)
## §6ky — `docs/ADVISOR_QUESTIONS.md` refreshed end to end against everything since the 09-09 pass (2026-09-15)
5-agent audit §6hq-§6kx; Q1 superseded (triangle leaders, opt-in), Q2 gains the annotation ceiling (F 0.67-0.81, §6kl), Q9 rewritten, Part 3 +5 concessions; register 836->838, test baseline 824->854. Ledger not retro-tagged. TBC1D3 CDKL/AE claims later corrected by §6l7.
(full: archive/2026-11/o1_ledger.md:21476)
## §6kz — Q9's guided subfamily pipeline run on gorilla for the first time: it does not close the gap
(heading continues: (2026-09-15)). Truth proxy = CHM13 landing (id >= 0.90; 8 chimeric). B1/B4 bars fail (named precision 0.756-0.923 < 0.95; duplicates 0.200 > 0.05); ZERO recoverable clades in any run. Cannot separate noisy label vs no gorilla subfamily structure. Q9 still not closed on gorilla.
(full: archive/2026-11/o1_ledger.md:21508)
## §6l0 — Q6 (apes): the current-default `gw_family_catalog` DOES port to chimp/orangutan — and its known
(heading continues: NPIP-fragmentation defect turns out to be general, not gorilla/ape-specific (2026-09-15)). chr16 PTR/PPY/HSA, no code/flag change: ports (27 and 75 NPIP-like copies). Fragmentation universal: human 26/56 = 46.4%, chimp 35.7%, orangutan 15.7%; largest family <= ~19% of NPIP-like copies.
(full: archive/2026-11/o1_ledger.md:21554)
## §6l1 — Q5: a real minority of the "we do not know" bin is now cross-species-confirmed biology (2026-09-15)
Reconstruction not exact (60 pairs vs §6ax 153). Canonical-unannotated 33/60; of those 10/33 = 30.3% replicate in both apes (raw), 7/33 = 21.2% strict; negative controls 0%. Remaining 70-79% undetermined, not artifact.
(full: archive/2026-11/o1_ledger.md:21606)
## §6l2 — Feasibility spike: scoring a read-through overhang against candidate copies with O2's PSV certificate — architecturally clean, empirically too thin and catalog-dependent to build (2026-09-15)
No architectural blocker; empirical yield thin: population A 5/14 = 36% (zero-PSV failures at identity 0.984-0.999), population B 3/20 (15%). Verdict: not worth building as a general feature; at most opt-in pre-checked use.
(full: archive/2026-11/o1_ledger.md:21665)
## §6l3 — IGV output: two real bugs fixed on real data, plus an origin_rejected colour, and a real worked example for the advisor (2026-09-15)
`bench/igv_tracks.py` fixes: CHROM from family span (409/409 correct) and minus-strand REF/ALT (0/409 mismatch vs 686/818 = 83.9% before); `origin_rejected` colour (14 reads). Worked example TRIM `CAFAM0`: 898 (92.3%) assigned, 897/897 unique-mapper agreement. Addenda: TRIM reads not ambiguous; `--dump-star` vs `.psv_copies.tsv` col_index mismatch (do not cross-reference). RETRACTION: "788 = 788" was two positions inside copy 3 (copy 1 AS 727, NM 20); certificate result (copy 3, p = 1.000e-24) stands.
(full: archive/2026-11/o1_ledger.md:21724)
## §6l4 — Identity-vs-AS-tie table: real ties concentrate at LOW identity, not high (2026-09-15)
§6aj dump (490 pairs, 79,175 reads): 0/490 zero-PSV pairs; 96/79,175 = 0.121% exact AS ties; 72/96 (75%) in families whose max pair identity < 90%, GWFAM48 only 24/96. Fork's first summary overstated (caught).
(full: archive/2026-11/o1_ledger.md:21901)
## §6l5 — Why most low-identity AS-ties happen: an incomplete candidate set, not two similar copies (2026-09-15)
Of the 72 low-identity ties, 90% (categories 1+2) have the tie partner outside every catalog copy (e.g. `SRR27178663.557589` 31.7 kb before copy 2); 6 genuine two-copy ties. Reading: catalog-completeness symptom.
(full: archive/2026-11/o1_ledger.md:21948)
## §6l6 — `--discover-copies` real-data acceptance test: correctly computed as an out-of-catalog site, but never
(heading continues: reported — both §6l5 hand-verified examples are singleton ties with no second corroborating read (2026-09-15)). Catalog is `mec/batch.copies.tsv` (no `psv.copies.tsv`). Pass criterion NOT met: site n=1 vs `TIE_PARTNER_MIN_SUPPORT = 2` (left unchanged). "Works elsewhere" 189-row claim superseded by the addendum.
(full: archive/2026-11/o1_ledger.md:21984)
## §6l6 addendum (post-final-review fix) — the 189-row report was inflated by a cross-family pooling bug and
(heading continues: intron-chained cluster spans; corrected re-run gives 28 rows, and the core negative finding SURVIVES (2026-09-15)). Fixes `f48c8543` / `62334a0c` / `b2bdbf3a`: 189 rows were 42 sites x up to 8 family_ids; max width 657,105 -> 1,681 bp. Corrected 28 rows / 5 families, 0 attribution violations; three §6l6 example rows RETRACTED. Core negative stands (n=1). Limitation: a row is a block cluster, not a copy.
(full: archive/2026-11/o1_ledger.md:22109)
## §6l7 — TBC1D3 Guitart truth correction: the name-mapped subfamily truth was wrong for 7/9 copies (2026-09-16)
Source `bench/TBC1D3_GUITART_TRUTH_CORRECTION.md`. Name map AE={T,E}, CDKL={D,K} replaced by M = {TBC1D3B, TBC1D3H} (robust) and CDKL = {TBC1D3K, TBC1D3} (figure-dependent); TBC1D3D unresolved; variants F / S / V3; ARI(old, F) = −0.059. 33 inventory items (RELABEL / RETRACT / SUSPEND); verdicts of §6jx (register row 822), §6js bars unaffected; register row 778 claim retracted/suspended.
(full: archive/2026-11/o1_ledger.md:22209)

## §6l8 — NPIP de novo node construction (development): the gene-body edge already gives family recall 1.000; three node fixes do not close the strict-precision gap — NO CANDIDATE (2026-09-17)
Can de novo nodes match annotated nodes under the shipped shared-definition comparison (NPIP dev, CHM13, 27 copies)? Gene-body edge gives FAMILY R 1.000 in every arm; loss is node quality (A4 shipped: P strict 0.290, 5/27 full-length vs A1 0.675). Fold/read-group fixes; best gate-passing G3 gains +0.036 F strict, below 0.05 bar.
Decision (pre-declared): NO CANDIDATE. Register row 843.
(full: archive/2026-11/o1_ledger.md:22616)

## §6l9 — The shared-junction conjunct is an ANNOTATION-side rule: it replicates on human, gorilla and chimp wherever junctions come from annotation, and fails wherever they come from reads or de novo nodes (2026-09-17)
Shared-junction conjunct t_J (k=2, δ=0) replicates where junctions come from annotation (human h_join 1.000000 -> 0.999363, gorilla -> 0.988399, chimp -> 0.998233) and fails from reads/de novo (gorilla de novo h_join 0.997868 -> 0.000000 by annihilating members). No level certifies any family. Useful for guided mode/reporting only. Register row 844. Readthrough question parked by user.
(full: archive/2026-11/o1_ledger.md:22690)

## §6m0 — The locus-formation gap to annotated loci: cause found (one rep per locus), one lever that works (read-isoform widening), and a measured observability ceiling on the rest (2026-09-17)
Cause of locus-formation gap: catalog emits ONE rep per locus before consolidate. Lever: read-isoform widening (k=8 closes 30.5% of family F-strict gap; k=5 F strict 0.515). Rest bounded by observability (all chains close 34.4%; oracle 52.1%).
Addendum 1: widening ADOPTED at k = 5 inside opt-in `RUSTLE_SHARED_DEFINITION` (`ISOFORM_MIN_READS = 5`, off with `RUSTLE_SD_READ_ISOFORM=0`); ideal sim containment 0.579 -> 0.934, NPIP full-length 11/27 -> 22/27, family/locus FLAT.
Addenda 2-6: strand asymmetry real (wrong-strand node blocks); V1-V3 not adopted (realistic sim flips verdict); V4 passes held-out except copy clause on U00; gorilla: opposite-strand loci 9 -> 3 but real family splits (FAMILY R 0.68 -> 0.64) -> V4 NOT adopted; V5 no-op by construction, V5b REFUTED on chimp. Next: a clause separating merges from splits.
(full: archive/2026-11/o1_ledger.md:22746)

## §6m1 — The no-readthrough counterfactual: removing readthroughs does not rescue the family definition, and the last obstruction is internal (2026-09-18)
Removing all 209 readthrough records leaves NPIP h_join 1.000000 at all four levels; member-preserving + t_J gives first non-empty L2 interval (0.746986, 0.989040], shipped 0.30 below it. Residual boundary is INTERNAL via member PKD1P6-NPIPP1: end of the line for edge rules. Addendum: PKD1P6-NPIPP1 is a REAL fused transcript (110 reads, all MAPQ 60, GT-AG).
(full: archive/2026-11/o1_ledger.md:23016)

## §6m2 — `--chimera-policy`: a truth-side lever for known readthrough records. Well-defined, clean label, NOT ADOPTED — every gain is denominator shrinkage (2026-09-18)
Prereg PREREG_chimera_policy_2026-09-18.md (strict/exclude/multi; evaluation lever only). CP-1 ROBUST (E1S beats E1 under all policies); CP-2 NO-OP held-out; CP-3/5: all four ΔF >= 0.02 cells disqualified as denominator shrinkage (e.g. NPIP L2 exclude +0.040, universe −15.4%); CP-4 NO. NOT ADOPTED; kept as reporting option, OFF by default.
(full: archive/2026-11/o1_ledger.md:23090)

## §6m3 — the NODE CUT rule: split a chimeric record at its parent boundary. SAFE (loses no member) but HARMFUL on every other axis — NOT ADOPTED (2026-09-18)
Prereg PREREG_node_cut_2026-09-18.md. 195 of 209 readthroughs cut. NC-1 safety PASSED (NPIP sensitivity 0.926 unchanged); NC-2 FRESH FAILED at every level (L2 precision 0.074 -> 0.005, ΔF −0.243); NC-3 NO. Cutting doubles (7 of 8 chimeras both halves in NPIP L2) and shorter nodes add edges. NOT ADOPTED.
Addendum: RNA reconstruction: junctions yes, complete copy-assigned chain only 3/26 NPIP members. §6m4 (sub-section): ceiling with all alignments (perfect O2) 10/26 complete chains; 40 of 249 junctions (16.1%) have no read.
(full: archive/2026-11/o1_ledger.md:23147)

## §6m5 — LOCUS ASSEMBLY at NPIP: the bottleneck is EXACT-CHAIN COLLAPSE, not ambiguity and not depth. O2 contributed 0 junctions (2026-09-18)
Prereg PREREG_locus_assembly_2026-09-18.md; `copy_assign` as shipped. AS-1 integrity PASSED; AS-2 FAILS (60.4% < 75.5%; 0/9). O2 contributed 0 junctions (661 contested: 206 assigned). Mechanism: exact-chain collapse + per-chain floor (dropped junctions median largest chain 2 reads vs 32). Pre-declared ceiling from §6m4: <= 10/26 complete, <= 83.9% junctions.
(full: archive/2026-11/o1_ledger.md:23285)

## §6m6 — assembler read-isoform widening PORTED (opt-in) and single-exon strand fixed; W-2 FAILED, and §6m5's "+21 junctions" projection is RETRACTED (2026-09-18)
Prereg PREREG_assembler_widening_2026-09-18.md. Shipped opt-in `copy_assign --read-isoform-k` (default 0). W-1 PASSED (byte-identical), W-3 PASSED (+50% tx for +2 junctions), S-1 PASSED (`RUSTLE_READ_STRAND=1` single-exon '+' 1.000 -> 0.753), S-2 PASSED. W-2 FAILED: measured +2 / +10 junctions; §6m5's "+21 junctions" projection RETRACTED. Loss sits between pass-1 and GTF emission.
(full: archive/2026-11/o1_ledger.md:23339)

## §6m7 — the NPIP "ceiling" was NON-CANONICAL JUNCTIONS, not missing reads. §6m4's depth framing is WITHDRAWN (2026-09-18)
Gate census: only real filter is `build_spliced_seq` (rej seq 343, 31.4% at k=0). Of 249 NPIP junctions, P(observed|canonical) = 0.981, P(observed|non-canonical) = 0.100; canonical = 209. §6m4's "16.1% absent information" WITHDRAWN: those junctions are likely annotation artifacts. Canonical-only truth: complete 3/26 -> 8/26 (k=0).
(full: archive/2026-11/o1_ledger.md:23386)

## §6m8 — the `build_spliced_seq` bucket split: 69% is ONE bad junction in an otherwise canonical chain, and `RUSTLE_JUNCTION_MAJORITY` recovers it (2026-09-18)
Default strict rejects a chain for ONE bad junction; 238 of 343 rejected (69.4%) are recoverable. `RUSTLE_JUNCTION_MAJORITY=1` + k=3: canonical 170/209, complete 10/26; +146 junctions, 95.9% read-carried. Not the genome-wide chr16 arm: default must NOT be flipped. Addendum: of 39 still-missing, 35 (90%) read-observed (algorithmic); idealized BAM run INCONCLUSIVE.
(full: archive/2026-11/o1_ledger.md:23438)

## §6m9 — the emission path localised without new instrumentation, and NOT enforcing canonical junctions gains exactly ZERO (2026-09-18)
`--gtf-copy-set` and unbounded `nc_max` RULED OUT (+0). Not enforcing canonicity recovers 0 of 40 non-canonical junctions. `RUSTLE_GATE_MIN_READS=1`: 170 -> 177, 11/26, but tx 1,456 -> 2,815. Recommended arm: `--read-isoform-k 3` + `RUSTLE_JUNCTION_MAJORITY=1`. Addendum: non-canonical "introns" include 1-2 bp annotation artifacts; ideal BAM fails because 30 byte-identical reads dedup to n_reads = 1.
(full: archive/2026-11/o1_ledger.md:23528)

## §6n0 — NPIP ceiling simulation: 25/26 complete with full-length reads. The assembler has NO algorithmic defect; 5' TRUNCATION costs 6 of 26 copies (2026-09-18)
Prereg PREREG_npip_sim_2026-09-18.md. SIM-0, SIM-4 PASSED; SIM-1 (bar >= 24/26) PASSED: FL 25/26 complete, 208/209; 5'-truncated 19/26; real 10/26. 5' truncation costs 6 copies. SIM-3: upper bound only. (Closing "stop optimising the assembler" later withdrawn by §6n1.)
(full: archive/2026-11/o1_ledger.md:23609)

## §6n1 — the 5'-recovery ceiling on REAL data, and a bigger defect it exposes: 18/26 copies already have ONE read carrying the whole chain, yet only 10 are emitted (2026-09-18)
Real data: 18/26 copies have one read carrying the whole canonical chain, union oracle 22/26, assembler 10/26. Selection defect named (+8); qualifies §6n0. (Reinterpreted by §6n2 as the primary-only read pool.)
(full: archive/2026-11/o1_ledger.md:23654)

## §6n2 — 5' RECOVERY DELIVERED: `RUSTLE_GTF_SECONDARY` takes complete NPIP chains 10/26 → 20/26. §6n1's "selection defect" was really the primary-only read pool (2026-09-18)
Primary-only pool caps at 10/26 (what the assembler emits); primary + secondary 18/26 available. `RUSTLE_GTF_SECONDARY` (default OFF, byte-identical) gives canonical 197/209, 20/26 complete; cost 9x transcripts (1,456 -> 13,132). Use for O1 "which copies exist", not quantification. Corrects §6n1's "selection defect" framing.
(full: archive/2026-11/o1_ledger.md:23692)

## §6n3 — echo collapse by AS-tie: IMPLEMENTED and REFUTED. The multimapper echoes are not clear losers WITHIN SCOPE, so a region-local AS filter cannot see them (2026-09-18)
`RUSTLE_GTF_SECONDARY_AS_RATIO` (default 0.0, inert) shipped; REFUTED: 13,132 tx at 0.0 / 0.98 / 1.0. Region-local scope: 83.4% of molecules contribute one record. Real fix would be a genome-wide best-AS pre-pass; do not re-propose per-region AS collapse.
(full: archive/2026-11/o1_ledger.md:23735)

## §6n4 — the chr16 arm for `RUSTLE_JUNCTION_MAJORITY` is RUN: the feared family FUSION does not occur; the real cost is +8.9% engulfment (2026-09-18)
Full chr16 arm: families 282 -> 290, max family 71 -> 71 (fusion DID NOT HAPPEN); engulfed 79 -> 86 (+8.9%). Default NOT flipped (user's call; second substrate wanted). Addendum: TBC1D3 + secondary 147/150 (98.0%), 10/12 complete; family-level lattice metrics unchanged by the --gtf flags.
(full: archive/2026-11/o1_ledger.md:23771)

## §6n5 — clustering-operator bakeoff on the L2 graph: several beat components, but every triangle-based one dissolves the modal family size. Louvain is the only survivor, and it costs the nesting theorems (2026-09-18, DESCRIPTIVE)
DESCRIPTIVE, development families, nothing adopted. Clique/truss operators raise F (4-clique NPIP 0.536) but cover 32.2% and dissolve 0/13 two-copy comps; Louvain NPIP F 0.287 -> 0.515 at full coverage, but breaks T1/T1′ nesting. Label propagation seed-dependent. (Louvain lead refuted by §6n6/§6o2.)
(full: archive/2026-11/o1_ledger.md:23832)

## §6n6 — Louvain does NOT replace MCL at the DNA level: worse on both catalogs, both views, with pairwise precision collapsing 0.815 → 0.144 (2026-09-18)
Held-out chr5/7/21 Soto (11 families): E1 MCL universe F 0.831 vs Louvain 0.547; pairwise precision 0.815 -> 0.144. MCL's leaving genes unclustered (1,867 of 2,715) does real work. Answer: NO, keep MCL at DNA level.
(full: archive/2026-11/o1_ledger.md:23864)

## §6n7 — Soto's families ARE a refinement of ours at L1/L2 (~80-90%), EXACT for NPIP; L3 is where we overtake them (2026-09-18)
Soto families (a cover: 148 genes in >= 2) contained in one of ours: L1a/L1b 90.0%, L2 80.0%, L3 20.0%; NPIP six Soto families all inside our one at L1/L2. Soto sits between our L2 and L3. Held-out: E1 61/76 = 80.3%, E1S 55/72 = 76.4%; failures concentrated in large families.
(full: archive/2026-11/o1_ledger.md:23897)

## §6n8 — readthroughs joining paralogs are PRIMARY, not split alignments: supplementary support is 0.1%. Secondary matters only for a handful of outliers (2026-09-18)
Records carrying 195 readthrough fusion junctions: primary 57.2%, secondary 42.7%, supplementary 0.1%. Median per-record secondary fraction 0.000; pooled share driven by outliers (ERV3-1-ZNF117, PKD1P4-NPIPA8, PKD1P3-NPIPA1).
(full: archive/2026-11/o1_ledger.md:23934)

## §6n9 — secondary-dominated readthrough filter SHIPPED (opt-in): the evidence-based discrimination register 844 said was missing (2026-09-18)
SHIPPED opt-in `bench/readthrough_secondary_filter.py`, `--max-secondary-frac` default 1.0 (no-op). At 0.50 flags 9 of 189, keeps PKD1P6-NPIPP1 (0.0265) and PDXDC2P-NPIPB14P; answers register 844's objection. Node-rule effect not measured.
(full: archive/2026-11/o1_ledger.md:23971)

## §6o0 — would a GRAPH reference dissolve these readthroughs? Evidence says NO: the supporting reads are RESIDENT, not reference-forced (2026-09-18)
Primary fusion-junction reads show no reference-forced signature (unflagged MAPQ 60 97.2%, median de 0.0020, soft-clip 0.000; flagged 80.8% MAPQ 60). A graph reference unlikely to dissolve THESE readthroughs; blind to unmapped/elsewhere reads (O3). Graph test not run (no assets).
(full: archive/2026-11/o1_ledger.md:24008)

## §6o1 — OVERALL family metrics on ALL Soto families, not the 11 pre-chosen: F drops 0.893 → 0.628, and the loss is RECALL (2026-09-18)
All Soto families: E1S F 0.893 -> 0.628, loss is recall. RETRACTED by §6o2 (truth set conditioned on prediction).
(full: archive/2026-11/o1_ledger.md:24046)

## §6o2 — ⚠CORRECTION to §6o1, and the operator question settled on held-out data: MCL wins, F 0.807 (2026-09-18)
§6o1 RETRACTED. Node-independent truth: chr5/7/21 E1S F 0.807 (82 fams), chr15/17 0.693, chr1 0.706; band F 0.69-0.81. Held-out operator test: MCL 0.807 vs Louvain 0.598, best identity components 0.759: MCL stays. Gap is precision on chr15/17 (0.367) and chr1 (0.475).
(full: archive/2026-11/o1_ledger.md:24083)

## §6o3 — `--min-shared-exon-frac` 0.30 → 0.60: mean bipartite F 0.764 → 0.786, leave-one-region-out validated (2026-09-18)
f_ex sweep, 3 regions: 0.30 mean F 0.764 -> 0.60 0.786; LORO never worse (+0.039 / +0.009 / +0.012), 0.60 picked in 2 of 3. chr5/7/21 0.845 at/above annotation ceiling. Candidate for default; NOT changed (user's decision).
(full: archive/2026-11/o1_ledger.md:24124)

## §6o4 — the RNA-side of the goal: f_ex SATURATES there, and the §6n9 readthrough filter gives a small real gain while "drop all readthroughs" is denominator shrinkage (2026-09-18)
RNA-side L2: f_ex saturates (+0.011 F, 0.563 -> 0.574). §6n9 filter drop: F 0.563 -> 0.586 at −3.8% denominator (real gain); drop-all-readthroughs F 0.652 but −13.5% truth genes (shrinkage). RNA-side precision 0.302 binding.
(full: archive/2026-11/o1_ledger.md:24160)

## §6o5 — the RNA-side level that matches Soto is w_98 ≈ 0.990, one notch above the shipped L3: F 0.586 → 0.659 and pairwise precision 0.302 → 0.656 (2026-09-18, DEVELOPMENT ONLY)
DEVELOPMENT ONLY, NOT A CANDIDATE: L3 sweep on NPIP/TBC1D3 interior optimum w_98 0.990, F 0.619 (shipped 0.98) -> 0.659, pairwise precision 0.302 -> 0.656; pairwise sensitivity 0.885 -> 0.411. Needs a second lattice region.
(full: archive/2026-11/o1_ledger.md:24197)

## §6o6 — plain identity-threshold components are consistently WORSE than MCL + the shared-exon conjunct (LORO, 3 regions) (2026-09-18)
Identity-threshold components: LORO picks 0.970 in all folds (mean F 0.729) but worse than MCL + f_ex 0.60 in every region (by 0.086 / 0.030 / 0.056). RNA-side 0.990 not corroborated (different quantities). MCL + shared-exon conjunct stays.
(full: archive/2026-11/o1_ledger.md:24228)

## §6o7 — RNA-side definition search on TWO substrates with leave-one-substrate-out: nothing beats the shipped triangle-leader rule (2026-09-18)
Human chr15/17/22 and gorilla RNA graphs, leave-one-substrate-out: identity components and leaders + identity cut do not beat shipped triangle leaders (human F 0.285; gorilla pick transfers as 0.250, WORSE). Goal NOT met at RNA level: high P (0.775 / 0.863), low R (0.175 / 0.428).
(full: archive/2026-11/o1_ledger.md:24253)

## §6o8 — ⭐THE GOAL IS UNATTAINABLE WITH THE CURRENT RNA EDGE GRAPH: pairwise recall CEILING is 0.052 (human) / 0.229 (gorilla), because >50% of truth families have NO EDGE at all (2026-09-18)
Pairwise recall ceiling on full edge graph: human 0.052, gorilla 0.229; >50% of truth families (104/206, 400/739) have no edge. Grouping rules saturated; work must go to edge construction. (Reinterpreted by §6o9 as a truth artifact.)
(full: archive/2026-11/o1_ledger.md:24300)

## §6o9 — ⭐THE RNA RECALL CEILING IS A TRUTH ARTIFACT, NOT AN EDGE DEFECT: only 6.7% of guided-truth family pairs align at all as spliced RNA (2026-09-18)
Mapping artifact and node fragmentation refuted; only 96 of 1,431 within-family pairs (6.7%) align as spliced RNA, 71 (5.0%) pass the edge gate: the ceiling IS the alignable fraction. Alignability-derived RNA truth (F 0.736) is CIRCULAR, not to be quoted. Non-circular RNA truth needed.
(full: archive/2026-11/o1_ledger.md:24336)

## §6p0 — ⭐AN RNA-LEVEL DEFINITION THAT SCORES HIGH ON ALL THREE AXES, against RNA-derived truth, with a held-out family check (2026-09-18)
Family = components of L3 graph at w_98 >= 0.995. Dishuck Iso-Seq truth NPIP F 0.222 (0.980) -> 0.833 (0.995); mean of two views 0.451 -> 0.758; held-out TBC1D3 Guitart pairs survive 0.985-0.995, 0.999 refuted. Recommendation 0.995 (range 0.985-0.995). (Superseded: shipped 0.995 in §6p1, moved to 0.985 in §6p2.)
(full: archive/2026-11/o1_ledger.md:24389)

## §6p1 — L3 RAISED to 0.995 (shipped); the MERGED-LOCUS layer; and node fragmentation measured (2026-09-19)
L3 raised 0.98 -> 0.995 (commit 3b2b3279, `L3_CUT = 0.995`; earlier lattice numbers at 0.98). Merged-locus cover layer (curated readthrough, >= 3 primary fusion reads, not secondary-dominated): PDXDC2P-NPIPB14P and PKD1P6-NPIPP1 become singletons. Proximity merge gap <= 5 kb consolidates 23% of nodes at 96% label consistency; not a recall lever.
(full: archive/2026-11/o1_ledger.md:24436)

## §6p2 — ⚠TESTING THE NEW L3 CUT ON 268 MORE FAMILIES: 0.995 does NOT generalise; the cross-substrate optimum is 0.985 (2026-09-19)
268 Soto families, L3-equivalent graphs: mean F 0.980 0.658, 0.985 0.664 (optimum), 0.995 0.630. Action: L3 moved 0.995 -> 0.985 (survives both analyses: NPIP 0.720, Soto 0.664).
(full: archive/2026-11/o1_ledger.md:24483)

## §6p3 — SELF-TUNING THE L3 CUT: three parameter-free rules tried, ALL FAIL. The threshold encodes biology the graph does not contain (2026-09-19)
Self-tuning L3 cut: largest gap and max persistence 0.451, quantile rule 0.758 on NPIP but mean F 0.575 on guided catalogs (vs fixed 0.985 0.664): ALL FAIL. Cut encodes paralog divergence (NPIP identity ~0.981); next design = estimate divergence from an independent source (untested).
(full: archive/2026-11/o1_ledger.md:24523)

## §6p4 — WHERE THE INFORMATION IS LOST: not CDS, not 3'UTR — it is the 5'UTR, and only where truncation bites (2026-09-19)
Per-base coverage by feature: CDS 99.8% in NPIP and TBC1D3; TBC1D3 5'UTR 75.3% vs 3'UTR 100.0% (5' truncation signature); NPIP 5'UTR 98.0%. CDS-based rules safe; 5'-end rules unreliable for TBC1D3-like families.
(full: archive/2026-11/o1_ledger.md:24559)

## §6p5 — ⭐AN EMPIRICAL RULE THAT DOES HELP: stop forcing one cut to serve two granularities — ADD L4. Both levels then beat any single cut (2026-09-19)
id_gx best single feature (accuracy 0.978 at 0.9971). One cut cannot serve two granularities: add L4 subfamily (id_gx >= 0.995) under L3 (>= 0.985): L3 lit F 0.773 (P 1.000), L4 iso F 0.833; nesting 30/30. Addendum: L4 IMPLEMENTED with selective trigger (>= 2 parts of >= 2 nodes; 4 of 38 groups refined, mean 0.774); certificate trigger vacuous.
(full: archive/2026-11/o1_ledger.md:24591)

## §6p6 — `--assemble-only`: an assembler-only mode, byte-identical GTF at 3.1x the speed (2026-09-19)
`--assemble-only` skips detection/assignment: 17.3 s vs 53.3 s on NPIP windows, GTF byte-identical there. Addendum CORRECTION: not identical in general (chr16 region 2,791 vs 2,777, strict superset; 77x faster). chr20 gffcompare: matching chains 358 (recall flags) / 345 baseline vs StringTie 331, FLAIR 264; StringTie leads precision (47.4 vs 33.4).
(full: archive/2026-11/o1_ledger.md:24666)

## §6p7 — precision and recall levers on chr20, and what `-R`/`-Q` really do (2026-09-19)
`bench/ism_collapse.py` (drop contiguous sub-chains) is the precision lever: baseline intron-chain Pr 44.6 -> 51.9 (vs StringTie 47.4), −21 chains. Transcript precision still trails (41.6 vs 47.1). `-R -Q` makes denominators tool-dependent: characterise one tool only.
(full: archive/2026-11/o1_ledger.md:24733)

## §6p8 — assembly polish: `--assemble-only` reaches StringTie's precision band (2026-09-19)
Prereg PREREG_assembly_polish_2026-09-19.md. Shipped `--assembly-polish <none|mono|full>` (none byte-identical). chr20 `full` beats StringTie on all four axes (337 vs 331 chains). Held-out chr11: H2/H3/H4 pass, H1 fails by 0.5 points (−3.50% chains vs ≤3% bar). Mono floor free (0 chains lost). Precision only MATCHED on chr11.
(full: archive/2026-11/o1_ledger.md:24772)

## §6p9 — locus isoform fraction: `--assemble-only` matches or outperforms StringTie in 19/20 cells (2026-09-19)
`--polish-isoform-fraction` (default 0.0), F = 0.02 fixed on chr20: matches/beats StringTie in 19 of 20 cells over chr20/11/7/14 (miss: chr11 transcript precision 49.2 vs 50.0). Addendum B grid rule (k=3, F=0.05) REFUTED held out, register row 855.
(full: archive/2026-11/o1_ledger.md:24807)

## §6q0/§6q1 — the shadow rule: 28/30 cells against StringTie over six chromosomes (2026-09-19)
Shadow rule `--polish-mono-shadow` (+ quantile 0.82): 28 of 30 cells vs StringTie over six chromosomes, held-out chr9 5/5; chr5 misses two recall cells. Anti-strand span not a criterion (register 857). Closing attempts worse: rows 858, 859, 860, 855. Quantile fitted on four, validated on two.
(full: archive/2026-11/o1_ledger.md:24837)

## §6q2 — chr5's two cells are a measured limit, not a tuning gap (2026-09-19)
`--polish-ism-3p` 19-21/30 REFUTED (register 861); chr5 single-exon recall not recoverable (register 862). Pooled six chromosomes: 2,800 vs 2,761 chains (+1.4%), transcript precision 46.8% vs 45.2%. Ceiling stays 28/30.
(full: archive/2026-11/o1_ledger.md:24867)

## §6q3/§6q4 — ISM support ratio 0.7 shipped; chr5 and chr11 proven disjoint (2026-09-19)
`--polish-ism-ratio` SHIP 0.7 (default 1.0). Final shipped setting `--assemble-only --assembly-polish full --polish-isoform-fraction 0.02 --polish-mono-shadow --polish-mono-quantile 0.82 --polish-ism-ratio 0.7`: 28/30, pooled 2,809 vs 2,761 chains (+1.7%). chr5 and chr11 feasible regions DISJOINT (register 863); `--polish-ism-escape` a documented dial; `--polish-fraction-exempt` no help (register 864).
(full: archive/2026-11/o1_ledger.md:24888)

## §6q5 — SQANTI3 on the polished assembler: pooled win on every category (2026-09-19)
SQANTI3: chr20 polish FSM share 36.1% -> 51.2%, rules PASS 75.9% -> 94.8%. Pooled six chromosomes beats StringTie on every category (PASS 94.3% vs 93.6%, FSM 2,819 vs 2,794).
(full: archive/2026-11/o1_ledger.md:24917)

## §6q6 — isoseq collapse and FLAIR: arms and mechanisms (2026-09-19)
FLAIR on six chromosomes: pooled chains ours 2,809 / StringTie 2,761 / FLAIR 2,309. isoseq collapse "not runnable" (register 865, later RETRACTED in §6q7). `--polish-fuzzy-junction` (default 0) REFUTED at 5 bp: 15/30 vs 28/30 (register 866); sub-5 bp differences are real NAGNAG sites.
(full: archive/2026-11/o1_ledger.md:24939)

## §6q7 — against the lab's own isoseq/StringTie/FLAIR runs, human AND gorilla (2026-09-19)
SUBSTRATE CORRECTION (register 867): six-chromosome bakeoff used `human_testis.t2t.bam`, lab uses `A119b.t2t.bam`. A119b chr20: ours 1,064 chains beats StringTie 861 and FLAIR 1,026 on all five; gorilla NC_073244.2 ours 1,575 beats both, k3+majority+polish 1,688 @ 31.9% beats isoseq 1,655 @ 9.0%. No single setting beats all three (register 868); fuzzy 5 costs 126 chains on gorilla (register 869); register 865 RETRACTED.
(full: archive/2026-11/o1_ledger.md:24973)

## §6q8 — script/file stocktake: what is live is now marked (2026-09-19)
`docs/ACTIVE_WORKING_SET.md` created (three tiers; 152 scripts named directly). Audit at `174ae65a`: 2,927 files, 1,480 candidates, 343 MB. `tools/cleanup_wave1.sh` dry-run by default; not yet applied.
(full: archive/2026-11/o1_ledger.md:25008)

## §6q9 — REPRODUCE.md and docs/DATA.md: the repo is now runnable by someone else (2026-09-19)
`REPRODUCE.md` and `docs/DATA.md` written. Libraries aligned with different minimap2 invocations (`human_testis.t2t.bam` 2.30-r1287 without `-uf`; A119b and GGO_mm 2.31-r1302 with `-uf`). A119b chr20 `recall_raw` 1,259 chains vs `recall_polished` 1,101. A119b / GGO_OR6737 accessions still to be filled by user.
(full: archive/2026-11/o1_ledger.md:25033)

## §6q9b — wave 2: `bench/` inverted (2026-09-19)
`tools/cleanup_wave2_bench.py`: keep what authoritative docs name + import closure: 353 kept, 428 archived (`bench/` 781 -> 353); tests 883/0 + 28/0. Redirect rule: missing `bench/X` -> `archive/bench/X`.
(full: archive/2026-11/o1_ledger.md:25056)

## §6r0 — a reduced file tree was BUILT AND RUN before trusting it (2026-09-19)
Reduced tree (905 files) built and run: 883/0, chr20 REPRODUCE GTF byte-identical (658 mRNAs / 336 chains). Defect 1 (register 871) `concat!` path missed; Defect 2 (register 872) fresh clone failed 10 integration tests (gitignored fixtures), fixed.
(full: archive/2026-11/o1_ledger.md:25080)

## §6r1 — waves 3 applied and the notebook moved behind a tag (2026-09-19)
Wave 3 applied: working tree 906 files; `archive/` tagged `notebook-2026-09-19` and removed. Fresh clone validated (883/0, 18 integration suites, chr20 GTF byte-identical). History rewrite refused (would invalidate cited SHAs).
(full: archive/2026-11/o1_ledger.md:25107)

## §6r2 — the depth-aware isoform fraction: premise refuted (2026-09-19)
Depth-aware isoform fraction premise refuted (register 873): the 65 lost chains are sole covers with median 2 reads. `--polish-fraction-min-reads` (default 0) a recall dial (floor 3: 1,092 chains). Human recall gap is assembly, not polish. Method note register 874.
(full: archive/2026-11/o1_ledger.md:25128)

## §6r3 — should O2 emit transcripts from secondary alignments? Measured: no (2026-09-19)
Gorilla catalog (78 fam / 509 copies): 24,543 AS-tied, 21 CONTESTED, 2 assigned: an O2 secondary emit path adds 2 molecules, not built (register 875). O2's decision set bounded by catalog coverage (register 876).
(full: archive/2026-11/o1_ledger.md:25160)

## §6r4 — the Y ampliconic genes are where O2 finally has a population (2026-09-19)
YAG catalog (8 fam / 30 copies) on A119b chrY: 53.1% AS-tied, 3,641 CONTESTED, assigned 12 (0.3%), tied 3,423 (94.0%). Only DAZ and RBMY assign. Binding constraint is decidability, not emission.
(full: archive/2026-11/o1_ledger.md:25188)

## §6r5 — TPM in the GTF and a locus BED (2026-09-19)
`--gtf-tpm` (default off) count-based TPM: ρ 0.879 vs StringTie TPM (length-normalised 0.714). `bench/locus_bed.py` (evaluation only): chr20 median size ratio 0.994, ±10% 55.6%; StringTie matches more loci (315 vs 288); our q25 0.676.
(full: archive/2026-11/o1_ledger.md:25222)

## §6r6 — pipeline state test against O1, O2 and the advisor's questions (2026-09-19)
O1 guided PASSES end to end on chr16 (f_ex 0.60): 90 families, all 21 NPIP genes in ONE cluster; chr20 0 families correctly. Register 878: symbol-named PAF gives silent 0-node graph. O2 on YAGs abstains (94.0% tied, 0.3% assigned); never present it as solving YAG copy assignment. Q9 untested.
(full: archive/2026-11/o1_ledger.md:25253)

## §6r7 — O1 DE NOVO tested end to end: recall 1.000, the defect is purity (2026-09-19)
De novo A119b chr16: 2,550 loci, 70 families; NPIP covered 21/21, dominant cluster 20/21; defect is purity (37 loci in 9 clusters). Register 879: 5 kb merge does not close it (dominant 20/21 -> 17/21).
(full: archive/2026-11/o1_ledger.md:25282)

## §6r8 — making genome-wide validation fast: what actually works (2026-09-20)
Crate levers already present (register 882); `--region-threads` 1.02x and NOT byte-identical (registers 880, 881); BGZF threads nothing (883); streaming claim noise (884). Addendum: dropping seq/qual under `--assemble-only` cuts RSS −59% / −62%; 4 chromosomes 426 s -> 114 s at `--jobs 4` (3.74x), byte-identical (885). `tools/genome_wide_sweep.sh` defaults `--jobs 4`.
(full: archive/2026-11/o1_ledger.md:25311)

## §6r9 — Python → Rust: the reproduce path no longer needs Python (2026-09-20)
`gff_to_gtf` and `locus_bed` bins replace Python, byte-identical; reproduce path needs no Python. Addenda: 58 dead scripts removed (register 886 false positives; 887 lattice tools not clone-runnable); absolute paths removed; wave 7 folded `bench/layer_order/` 12 -> 2 (register 1104: earlier fix one dirname short; 61/61 identical) and `bench/soto/` into `soto_replication.py` (reproduces §6ip, e.g. median ARI 0.6959).
(full: archive/2026-11/o1_ledger.md:25371)

