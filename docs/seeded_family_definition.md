# The seeded multi-copy family: definition, properties, certificate

Compacted 2026-11-01. The verbatim pre-compaction text is `docs/archive/2026-11/seeded_family_definition.md`
(every entry ends with its line there); git history has every version. Numbers measured on CHM13 (human) and
NHGRI_mGorGor1-v2.0 (gorilla); scripts in `bench/crossspecies/`. Anything not measured is **OPEN**.

## Current state (2026-11-03)

- **Shipped default = §0★ unchanged**: reads → loci (`copy_assign --assemble-only`, secondaries at AS ≥ 0.98 x best,
  `--bridge-regroup f1v2`) → most-reads representative, positional exon sum → `mcl_families --from-gtf
  --min-exonic-bp 1 --min-shared-exon-frac 0.60` (asm20 all-vs-all, ≥ 300 bp at ≥ 70%, cov ≥ 0.30 longer or ≥ 0.70
  shorter, MCL I 2.8, nucleotide only) → copy table `PREFIX.fam.copies.tsv/.fa` for O2. MODULE_STATUS:
  `denovo_pipeline.rs`, `family_detect.rs` SHIPPED-DEFAULT; `fam_from_gtf.rs` OTHER-BINARY (`mcl_families --from-gtf`);
  `family_hierarchy.rs` (`--emit-hierarchy`) OPT-IN, annotation only.
- **Opt-in / proposed, not default**: legacy `gw_family_catalog`; protein attach (`tools/protein_attach.py`, manual);
  `RUSTLE_SHARED_DEFINITION` (§0★★); nested lattice and certificates (§0★★★, PROPOSED); representative variants.
- **Family-rule experiments 2026-10-29..11-01 (register 1469-1480), defaults unchanged:**
  - 1469 U1 intron-column lift edges: FAIL on HG002 (pair F 0.411 < 0.413, precision 0.415 < 0.478); opt-in.
  - 1470 DISPERSION (drop families on >= 4 chromosomes): FAIL on A119b and OR6737.
  - 1471 SKETCHFAM (FracMinHash families replacing minimap2): FAIL (A119b pair F 0.783 -> 0.561, TSPY 34/34 -> 0/34).
  - 1472 SEMIGUIDED (de novo + confirmed `semisec` candidates): PASS on both held-outs; not the default.
  - 1473 REPRULE-2 (representative = reads x (junctions+1)): FAIL on HG002 B1 (NPIP 11 -> 10).
  - 1474 REPRULE-3 arm S (edge weight saturated at 0.95): FAIL; split is from missing cross records, not weights.
  - 1475 CHAIN-GATE `--chain-gate 3`: FAIL (KB3781 precision 0.341 -> 0.321; HG002 NPIP 11 -> 10); opt-in only.
  - 1476 LOCUS-ABSORB: dev gate recall missed by 0.002 (0.729 vs >= 0.731); not adopted.
  - 1479 FAMASM multi-copy assembly mode: FAIL on HG002 (multi recall 11.3 -> 12.7 but multi precision 14.0 -> 8.5); not adopted, prototype `bench/stgaps/`.
  - 1480 protein reconciler on de novo families (A119b, descriptive): NULL; mainstay losses are assembly/nucleotide
    losses, not protein drift.
- **Families input (2026-11-02, r1490/1491/1498):** families are built from the unlayered `PREFIX.gtf`; the assembly's end layers (median 5' ends, 3' UTR extension, single-exon cluster models) passed their own held-out but cost pair precision 0.726 -> 0.700 / recall 0.741 -> 0.714 on A119b through the families stage, so they live in `PREFIX.final.gtf` only. Annotation-guided fused-loci cut PASS (r1495, recommended opt-in); shared-exon split FAIL (r1496/1497); readthrough structural rule FAIL (r1492); CAT truth split below bar (r1494).
- **Open**: held-out precision/recall of lattice levels L1-L3 / P1 (§0★★★.6); de novo node extent (H1) remains the
  dominant error (§4, §4a); reach bound §4a.

## INDEX

Original map of 14 sections; titles carry the verdicts (auto-derived tags scored 11/22 and were removed). Family vs
SD / duplicon / expansion terminology: `docs/archive/2026-10/merged/TERMINOLOGY_FAMILY_SD_DUPLICON_EXPANSION_2026-10-07.md`.
(full: archive/2026-11/seeded_family_definition.md:9)

## 0★. THE SHIPPED DEFAULT (2026-09-25, user decision) — ONE DE NOVO FAMILY DEFINITION AT THE RNA LEVEL

**What ships, stated once.** Everything below this section (§0★★★, §0★★, §1★ and the older sections) is design
history, proposals or opt-in machinery; where it says "current definition" it means the design reference of its date.
The default is:

1. **Reads → loci.** `tools/rustle_pipeline.sh assemble`: `copy_assign --assemble-only --genome-wide
   --assembly-junctions strict` + the polish flags, from primary alignments plus secondary alignments scoring ≥ 98% of
   the read's genome-wide best alignment score (`RUSTLE_GTF_SECONDARY_AS_RATIO=0.98` with the `as_table` best-AS
   table; register 1060). Since 2026-09-29 the assembly ends with the bridge-aware regrouping `--bridge-regroup f1v2`
   (the user's decision; `bench/ASSEMBLY_POLISH.md` addendum 3): a readthrough transcript that is the only link
   between two read-proven pieces of its `gene_id` and carries fewer reads than each becomes a `fusion_of` relation,
   and the pieces split. A locus is one `gene_id` of `PREFIX.families.gtf`, the assembled GTF without those bridges
   (`RUSTLE_BRIDGE_REGROUP=off`: one `gene_id` of `PREFIX.gtf`, the 2026-09-25 rule).
2. **Locus → representative: the positional exon sum.** The locus's representative is the transcript with the most
   reads (ties: the longer span). Its exon coordinates come from the reads; its bases from the genome ("exon-sum:
   the reference bases at the candidate's exon ranges", `docs/archive/2026-09/PSEUDOCODE_2026-09-08.md`). Loci whose exon unions
   overlap are folded into one node.
3. **Representatives → families.** `tools/rustle_pipeline.sh families` = `mcl_families --from-gtf --min-exonic-bp 1
   --min-shared-exon-frac 0.60`: all-vs-all of the loci's genomic spans (`minimap2 -x asm20 -c -X -N 50 -p 0.1
   --secondary=yes`); an edge needs ≥ 300 bp at ≥ 70% identity, coverage ≥ 0.30 of the longer locus's exonic length
   — or, since 2026-09-29, ≥ 0.70 of the SHORTER locus's exonic length, which is then the edge's coverage weight
   (`--min-cov-shorter 0.70`, the §6x4 containment escape, register 1006/1014; the user's decision;
   `--min-cov-shorter 0` = the 2026-09-25 rule; known regressions: NPIP-guided Soto F .833 → .800, semi-guided
   precision .973 → .833, register 1007/1009) —
   exonic bases on both sides, and one record joining exons to exons over ≥ 60% of the smaller locus's exonic length
   (`annotation_families.rs`); Markov clustering, inflation 2.8; families of ≥ 2 loci. Nucleotide only.
4. **The same families feed copy assignment.** With `--emit-units` the stage writes the copy table
   `PREFIX.fam.copies.tsv/.fa` (one copy per member locus: its representative's exons and spliced exon sum, in the
   contract `copy_assign --families` reads; `docs/archive/2026-09/PREREG_families_copy_table_2026-09-25.md`). O2 assigns tied reads
   among the copies of these families.
5. **Guided mode** (Fig. 7): the same family rule on the annotated gene and pseudogene bodies. The thesis goal is to
   reduce the difference between the two modes.

**Not the default.**
- `gw_family_catalog` (the `catalog` stage: primaries only, span-aware POA collapse, exon-sum k11 E_r edges,
  γ-quasi-clique) is **LEGACY**: kept runnable and scored only as a labelled comparison.
- **Protein** is an **optional, manually invoked extra-sensitive step** (`tools/protein_attach.py`;
  `docs/archive/2026-09/PREREG_protein_attach_2026-09-25.md`): for families with disconnected or strongly diverged members it
  attaches missing loci at the protein level to ONE existing family, never merges families, never runs by default;
  it appears only in a supplementary figure. `gw_family_catalog --protein-tail` stays an opt-in of the legacy catalog.
  The §0★★ clause-3 "miniprot fallback" is a label idea from the ledger, not code.
- `RUSTLE_SHARED_DEFINITION` (§0★★ clause 4 leaders), the SD-core certification, the nested lattice (§0★★★) and the
  representative variants `RUSTLE_LOCUS_EXON_UNION` / `RUSTLE_LOCUS_UNION_MIN_TX` / `RUSTLE_COTHREAD_REP` are opt-in
  or proposed (register rows 303, 1041, 1055, 393/484/841 for the representative variants).
- **The closed loop** (families → union-certificate assignment of tied reads → each assigned read given to its copy
  → re-assembly) is a pre-registered pass-2 experiment (`docs/archive/2026-09/PREREG_tied_read_loop_2026-09-25.md`), not part of the
  definition; families are frozen from pass 1 in it (register 395).

**How it is judged** (`docs/archive/2026-09/PREREG_genome_wide_families_2026-09-25.md`, Amendment 1; `figures/README.md`). External
references only in the main figures: Ensembl Compara (human; pairs by protein-identity band in Fig. 6, families at
Primates in Fig. 7), Soto 2025 (labelled not independent: the 0.60 threshold was chosen against it, register 903),
Liftoff copies (the self-lift's record / extra-copy pairs, every species; Figs 7 and 8). The annotation's
protein-homology families are a secondary reference of the supplementary figures only (fold-level families, mostly
outside the nucleotide regime: register 1101).

(full: archive/2026-11/seeded_family_definition.md:34)

## 0★★★. NESTED EDGE-TEST LATTICE (2026-09-16) — PROPOSED: nesting theorems proven, T3 conditional, precision OPEN; a design, not the default (§0★)

User-approved design (2026-09-16). Every level = connected components of ONE evidence graph under nested edge tests,
so nesting is a theorem at every level including the expression view. Motivation: (1) "every RNA family lies inside one
DNA family" holds for components, not for clause-4 leaders (6-node counterexample: reads A 100, B 50, C 10, D 5, E 4,
U 0; edges A–B, A–U, U–C, B–C, C–D, D–E; expressed-only leaders {C, D} cross two DNA families) or MCL; (2) independently
built layers are not nested (0.521 of E1 DNA-family pairs are protein-family pairs; 0.253 the other way, §6ko AN-2).
Cost of L0 = JOIN: on Soto TBC1D3 REFINE F 1.000 vs JOIN 0.818 / 0.643. Spec
`docs/superpowers/specs/2026-09-16-family-layer-order-design.md`; `bench/LAYER_ORDER_NPIP_TBC1D3.md`;
`LO/` = `/mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3/`; random-instance checks
`LO/lattice/proofcheck_math/lattice_bruteforce{,2}.py`. (full: archive/2026-11/seeded_family_definition.md:90)

### 0★★★.1 Objects

- **Nodes V = copies** (guided: annotated gene bodies + candidate loci; de novo: read-supported loci at gene level).
  Fixed within a run; NOT fixed across runs (de novo extents are read-derived; guided widths depend on seeds), so any
  cross-node-set statement (T3(b), T2 across samples in de novo) needs a **frozen** node set. Readthrough exclusion
  rule proposed, OPEN (would drop NPIP member PKD1P6-NPIPP1).
- **Evidence A**: DNA genome-alignment records (spliced transcript `-x splice -uf`, gene body `-x asm20`), blastp HSPs
  (guided only), SEDEF calls (carried, unused by tests). A DNA record links u→v when its target is off u's locus and
  touches v (exon block on v's exons, or gene body). sx(r) = min exon bases of u, v in r. R_D(e) = DNA records of e.
- **E = every pair with ev(e) ≠ ∅**; one graph for all levels.
- **Tests** (monotone form):
  - **t_D**: an exon record (identity ≥ 0.80 over ≥ 0.50 of the transcript) or a gene-body chain (pooled identity
    ≥ 0.80, ≥ 0.50 of the shorter body, target span on ≥ 1 exonic base); sense condition when both spliced. Shipped
    `gene_body_chains` is greedy with pooled identity and NOT monotone (203 of 5,000 random instances lose the edge).
  - **t_P**: both nodes carry a protein (longest annotated CDS, ≥ 10 aa, r2 biotype rule); union of HSPs at e ≤ 1e-5
    covers ≥ 0.30 of the longer protein. Shipped `protein_families.edges_from` (greedy by bitscore) NOT monotone (399 of
    5,000). L0 is guided-only.
  - **f_ex(e)** = max over R_D of sx(r) / min exonic length (adapted from `min_shared_exon_frac`, which uses E1 records only).
  - **w_98(e)** = max identity over one record with sx(r) ≥ 0.30 x min exonic length (single-record SD98 reading;
    gap-excluded identity chosen in the NPIP/TBC1D3 run).

| level | name | test | note |
|---|---|---|---|
| L0 | superfamily (guided-only) | t_P ∨ t_D | withdrawn as a level by §0★★★.7 |
| L1 | family | t_D | clause 2 |
| L2 | shared-exon unit | t_1 ∧ f_ex ≥ 0.30 | adapted from §6ks, §6kt |
| L3 | ≥ 0.985 identity unit (not the subfamily) | t_2 ∧ w_98 ≥ 0.985 | raised from 0.98 on 2026-09-19 (§6p0): NPIP vs Dishuck Iso-Seq F 0.833 (vs 0.222 at 0.98); 0.999 refuted |
| L4 | subfamily (Iso-Seq group), applied selectively | t_3 ∧ w_98 ≥ 0.995 where the L3 group subdivides (≥ 2 parts of ≥ 2 nodes) | new 2026-09-19 (§6p5); TBC1D3 specificity 0.000 → 0.750 at no recall cost |

L2/L3 names are not "SD"/"subfamily": L2 also holds for retrocopy–parent pairs; clause 5 stays the subfamily
definition. Constants: t_D 0.80/0.50; t_P 1e-5, 0.30, 10 aa; t_2 0.30; t_3 cut and its 0.30; expression reads ≥ 3;
de novo read-locus nodes ≥ 3 unique reads, exon depth ≥ 2. Levels G_k = (V, E_k), groups π_k = components.
Expression view X = expressed copies, groups = comp(G_k[X]); counts: any-overlap (lattice) or u ≥ 3 (MAPQ ≥ 1,
`bench/interval_expression.py`). (full: archive/2026-11/seeded_family_definition.md:128)

### 0★★★.2 Theorems

- **Lemma 0**: F ⊆ F′ ⇒ comp(F) ≤ comp(F′); comp(F ∪ F′) = join; comp(F ∩ F′) ≤ meet (can be strict). So L0 = join of
  protein and DNA components.
- **T1** levels nest; **T1′** all levels form one laminar tree; **T2** expression view nests at every level
  unconditionally (comp(G_{k+1}[X]) ≤ comp(G_k[X]) ≤ π_k|X); pooled views (full vs expressed across levels, nested
  samples at different levels, two incomparable samples) are NOT laminar (718 / 3,764 of 14,410 random draws). Zero
  violations on data are a code check, not evidence.
- **T3** (more information only merges or only splits; no copy moves). Hypotheses H1 fixed node content, H2 pairwise
  evidence, H3 monotone tests (holds for t_D exon disjunct, f_ex, w_98; fails for shipped body chains and t_P).
  (a) stronger tests ⇒ splits only (no hypotheses). (b) more copies (H1, H2) or more homology (H1-H3) ⇒ merges only.
  (c) more expression on one graph ⇒ merges only.
- **Counter-hypotheses**: 1 non-fixed node content; 2 evidence depending on other nodes (blastp E-value vs db size,
  `-N`/`-p`, best-hit tests, tie order); 3 non-monotone tests (greedy HSP, greedy chains, pooled identities);
  4 group-dependent grouping (leaders, MCL, split systems, seed neighbourhoods); 5 unique-read expression
  (TBC1D3P1-DHX40P1: 72 any-overlap, 0 unique reads); 6 monotone third-node filters (triangle support = 3-truss,
  4-truss, bridge removal) keep T1-T3. (full: archive/2026-11/seeded_family_definition.md:265)

### 0★★★.3 Filtration view: thresholds are cuts, not definitions

Sweeping one field w over base test b gives nested partitions = the single-linkage dendrogram; bottleneck h(u, v) from
a maximum spanning forest; 1 − h is a pseudo-ultrametric (ties at 1.0 common: 2,946 of 82,491 E1 identities). T3(b)
and T2 extend to the whole dendrogram. Each constant is a named cut of a sweep; report h per pair and merge heights.
Two-field sweeps have no single dendrogram. Carlsson & Mémoli (JMLR 2010) uniqueness applies only to a derived metric
space. (full: archive/2026-11/seeded_family_definition.md:434)

### 0★★★.4 Where the existing groupings sit (none is a lattice level; all stay report-only)

§6ko protein MCL families, MCL guided catalogs/truths (0.70/0.30), guided seed neighbourhoods, clause-4 leaders,
clause-5 split system, spec JOIN/REFINE: none is components of a monotone edge test; all report-only. Each leader family
lies inside one level-k group of its graph. Layer-order study: P vs D order undetermined (of 9 variants 6 D > P,
2 P > D, 1 tie); T2 precondition held for all 30 refinement pairs. (full: archive/2026-11/seeded_family_definition.md:472)

### 0★★★.5 The cost, stated first: components chain

Components (closest measured L1) have the best recall and the lowest precision (expressed-guided R_G / P_G):

| substrate | nodes | components | bridge-split | triangle leaders | one-hop leaders |
|---|---|---|---|---|---|
| gorilla hold-out, AB2 nodes | 2,431 | AB2 0.5088 / 0.3100 | — | — | AC1 0.3401 / 0.5806 |
| gorilla hold-out, + read-locus nodes | 2,664 | AC2 0.6950 / 0.3001 | 0.6100 / 0.3023 | AD2b 0.6046 / 0.4191 | AC3 0.4170 / 0.5365 |
| gorilla fresh (§6kd) | 2,853 | AC2 0.2273 / 0.3767 | 0.1240 / 0.5172 | 0.1983 / 0.6761 | AC3 0.1942 / 0.7015 |
| human chr15/17/22 (§6ke) | 6,220 | AC2 0.1606 / 0.1718 (F 0.522) | — | 0.0723 / 0.6822 (F 0.558) | AC3 0.0425 / 0.7248 |

Chaining is real (AB2 largest component 85 loci over 16 guided clusters); bridge removal failed (P_G 0.3023); L0 on
plain §6ko is degenerate (TBC1D3 component ≥ 4,522 genes); approximate (t_P ∧ aa ≥ 0.50) ∨ t_D still joins NPIPA1 and
TBC1D3 (1,017 genes, 7-hop path); guided mode loses "never chains" (approximate L3: 14-node NPIPA + PKD1P + readthrough
component; readthroughs removed → NPIPA 7-node component). L3's cut sits near paralog divergence (median p-distance
NPIP 0.0190, TBC1D3 0.0103). Precision options: (i) stricter pairwise tests / L2 (§6ks Soto F 0.831 → 0.881, never on
de novo: OPEN); (ii) triangle filter (dissolves every 2-copy family); (iii) leaders report-only; (iv) readthrough node
rule (not scored). (full: archive/2026-11/seeded_family_definition.md:490)

### 0★★★.6 Status

PROVEN: Lemma 0, T1, T1′, T2, T3(a), T3(b) under stated hypotheses, T3(c). Not met by shipped code: H3 for t_P and
body chains, H2 for t_P, H1 in de novo. L0 withdrawn (§0★★★.7). **OPEN: precision/recall of L1-L3.** NPIP/TBC1D3
descriptive run (`bench/NESTED_LATTICE_NPIP_TBC1D3.md`): L0 group 5,764 genes; L1 341 / 97; L2 142 / 24; L3 55 (25) /
11 (11); Soto in-group pair precision L2 0.187 / 0.392, L3 0.259 / 1.000; run used approximations (greedy t_P, catalog
all-vs-all t_D). Held-out not run: reserved human chr8/9/10/11 and a new gorilla contig set; name each level's truth
first, never score a level against a truth built with its own conjunct. Until then §0★★ clause 4 (triangle leaders)
stays the confirmed de novo grouping; no code changes. Audit notes reconcile 2,426 / 2,300 / 1,017 component sizes to
denominator and v-exon/strand variants. (full: archive/2026-11/seeded_family_definition.md:568)

### 0★★★.7 AMENDMENT (2026-09-17, user-approved direction): TWO PARALLEL STACKS REPLACE THE L0 JOIN; CERTIFICATES ARE THE REPORTING UNIT

PROPOSED; no code changes. (a) Join dropped: it puts NPIP in a 358-node component, TBC1D3 in 5,450. (b) DNA stack
L1/L2/L3 (both modes) and protein stack P1 = t_P (guided only; P2 = t_P ∧ aa ≥ θ not adopted); a family claim names its
stack. (c) Meet M = t_P ∧ t_D is a floor: no non-member for either family, but NPIP in 6 parts, TBC1D3 in 4.
(d) **Certificate**: h_join(S) = max boundary edge weight, h_split(S) = bottleneck of S's maximum spanning tree; S is an
exact component for cuts in (h_join, h_split]; a family claim needs a non-empty interval containing the level's cut,
reported as the margin. (e) t_P qualifies only with fixed `-dbsize` and union cover; then **the 21 NPIP proteins are an
exact P1 component for cover cuts in (0.128250, 0.840580]** (first certified family statement); TBC1D3 fails (exact only
in (0.688525, 0.900000]; 5,275 proteins at the shipped cut). (f) No DNA level certifies either family; a 320-row
strengthener sweep never made NPIP's margin positive (boundary held by co-duplicated neighbours at identity 1.000:
CLN3, EIF3CL, LOC100190986, LOC124907830/845, LOC128966632); F rose while no certificate appeared, so F is not an
isolation guard. (g) Open: held-out precision/recall, monotone body disjunct, readthrough node rule, no de novo protein
stack. Evidence `bench/FAMILY_CERTIFICATES_NPIP_TBC1D3.md`, `bench/LATTICE_RULE_STRENGTHENERS.md`.
(full: archive/2026-11/seeded_family_definition.md:652)

#### 0★★★.7c Shared-junction diagnostic (opt-in, never a default, not part of any shipped cut)

Witness record must carry ≥ k of the query copy's junctions onto the target's (δ bp); monotone (0 violations in 63,408
cells); only k = 2, δ = 0 has evidence. Cost is arm-dependent: annotated junctions cost nothing (NPIP outside 96 → 10),
read junctions break NPIP in two, de novo structure breaks NPIP into 5-6 and TBC1D3 into 8-10. Never certifies NPIP
(h_join 0.999363 vs h_split 0.971665). "Guarded" variant REFUTED. `bench/JUNCTION_AND_READTHROUGH_RULES.md`.
(full: archive/2026-11/seeded_family_definition.md:704)

#### 0★★★.7d Cut windows are reported, not acted on

Non-empty certificate intervals that miss the shipped cut are CUT WINDOWS, not grounds for a rule change. On record at
L3 under §0★★★.7c: four Dishuck groups exact in (0.997129, 0.997680] (the conjunct adds three; B15 already exact in the
baseline: never quote four gains); NPIPB exact in (0.980559, 0.983454], missing 0.98 by 0.000559. Needs a held-out
prereg. (full: archive/2026-11/seeded_family_definition.md:727)

#### 0★★★.7e Readthrough node rule (opt-in, measured, not adopted)

A readthrough record is dropped when it contains ≥ f of another node's exonic bases: f = 0.50 / 0.80 / 0.90 / 1.00 drops
178 / 128 / 92 / 20 of 209; no NPIP/TBC1D3 member dropped (PKD1P6-NPIPP1 survives). Invisible alone (NPIP h_join stays
1.000000); a containment rule only sees nodes (PDXDC2P-NPIPB14P survives). (full: archive/2026-11/seeded_family_definition.md:742)

## 0★★. DESIGN DEFINITION — ONE COPY GRAPH FOR BOTH MODES AND BOTH LEVELS (2026-09-14; opt-in `RUSTLE_SHARED_DEFINITION`, NOT the shipped default: §0★)

Design reference of 2026-09-14 (formerly titled "CURRENT DEFINITION"); clause 4 ships only as opt-in
`RUSTLE_SHARED_DEFINITION`; clause-3 miniprot fallback is not code. Nodes are genomic intervals, edges are homology on
assembly sequence, RNA decides which intervals are nodes and how far they extend, RNA never decides an edge.
1. **Node = a copy** at gene-body width with exon structure. Guided: annotated genes + chain-first candidate loci (§6js).
   De novo: read loci consolidated to gene level before edges; read-locus nodes (≥ 3 unique reads, exon depth ≥ 2)
   recover 43/58 node-less expressed loci, missing-node pairs 219 -> 61 (§6ka).
2. **Edge**: EXON edge (spliced transcript on v, identity ≥ 0.80 over ≥ 0.50) OR GENE-BODY edge (identity ≥ 0.80, ≥ 0.50
   of the shorter body); the union is required (NPIP held by gene bodies, AMY by coding sequence; §6jo, §6jq, §6jr).
3. **Edge labels = mechanism** vs the parent: GENOMIC (≥ 1 kb shared non-exonic non-repeat sequence; SD 9/9 vs retro
   0/5), RETROCOPY (introns lost at junctions; 8/8 fresh, 0/10 SD false calls), UNRESOLVED. Poly(A)/TSD reported only.
4. **Family = triangle-supported leader neighbourhood** (§6kd): leader by read support, + unassigned neighbours +
   copies adjacent to ≥ 2 of that star. Fresh gorilla R_G 0.198 vs 0.120, P_G 0.676 vs 0.537 over LCS union. Not a
   lattice level (read-order dependent); stays the confirmed de novo grouping until §0★★★.6 is measured. Superseded
   wording: family = connected component (γ = 0.20 dropped; components chain, P_G 0.310).
5. **Subfamily** = cluster of a supported pairwise-compatible split system (ML tree, SH-aLRT > 75, exon and intron;
   NPIPA|B 10/10); split-dominance estimator inadequate (73 vs 85); conflict index reported (NPIP/TBC1D3 ~0.36, AMY 0.015).
6. **DNA vs RNA**: RNA graph = DNA graph induced on expressed copies (reads ≥ 3), so every RNA family lies in one DNA
   family (components only; not leaders). Gap classes: (a) MISSING NODE, (b) MISSING EDGE (both de novo errors),
   (c) UNEXPRESSED BRIDGE (level difference). §6jz: missing edges 126 -> 60; missing nodes 219 pairs (29.6%) dominant.

Known exposures: 0.80 / 0.50 / 1 kb constants remain; subfamilies not in the Rust catalog; de novo instantiation
untested on a held-out substrate. (full: archive/2026-11/seeded_family_definition.md:755)

## ⚠⚠ TIER NOTICE — READ BEFORE QUOTING ANY NUMBER IN THIS DOCUMENT (2026-08-10)

Every `E_r` number is a statement about ONE alignment tier. SHIPPED = `minimap2 -c -X --no-long-join -t N -k 11 -w 5`
(`ER_TIER_FLAGS`, `denovo_pipeline.rs`); PANEL = `minimap2 {-x asm20 | -k 11 -w 5} -c --eqx -N 200 -p 0.02` (no `-X`,
old `bench/crossspecies/` scripts). On byte-identical FASTA the partition differed on 4/14 panels, edge count on 10/14.
FIXED 2026-08-10: one `ER_TIER_FLAGS` / `er_tier_command()` mirrored in `bench/soto/rustlib.py`; X.4 added
`er_sensitive_only()` / `er_primary_tier()` for both call sites and `<prefix>.rule.tsv`; parity =
`diff <prefix>.rule.tsv <prefix>.refine.rule.tsv`. Superseded 2026-08-13 (O-3): the shipped-default diff is one line on
the substrate axis (O1 swaps substrate, O2 refine unions one in), empty on every other knob; never "the same rule".
Scope: all-vs-all only; seed→genome pass unchanged. Resolved 2026-08-11: §2 numbers are SHIPPED-tier; P2 moved by one
edge (gorilla NPIP c = 0.4989); P4b RNA 338 → 260 edges; `E_RNA ⊂ E_DNA` is panel-specific (fails on 80-node
`o1_closure`: 33 RNA-only of 481; 1 of 33 joins two DNA components). §3b's N0 table is at a third tier, OPEN.
(full: archive/2026-11/seeded_family_definition.md:846)

### Coverage form (defect M1), fixed the same day

`c` was `(qend-qstart) / min(qlen,tlen)` (exceeded 1.0 on 21/255 records, max 1.142). Now the numerator's axis follows
the denominator: `qlen ≤ tlen ? (qend-qstart)/qlen : (tend-tstart)/tlen` (bounded, symmetric). Cost on the 25-region
gorilla control panel: 1 region changes (HERC2 gains a true third copy, identity 0.9940), 0 edges lost, 2 gained,
single-copy controls 0/7. (full: archive/2026-11/seeded_family_definition.md:915)

## 0. What the definition must survive

Three forms of arbitrariness to close, each by a measured property: (1) annotation-dependence (seed choice),
(2) threshold-dependence (γ), (3) circularity (evaluation consuming the labels that produced the answer).
(full: archive/2026-11/seeded_family_definition.md:936)

## 1★. THE OBJECT, RESTATED — ONE OBJECT AT TWO LEVELS (2026-08-14)

> **A multi-copy gene family is a connected block of the graph whose nodes are genomic intervals and whose edges are
> contiguous, high-coverage homology computed on assembly sequence. RNA determines which intervals are nodes and how far
> they extend. RNA does not determine an edge.**

Restatement only; §1 stays normative. (full: archive/2026-11/seeded_family_definition.md:949)

### 1★.1 Why this is one object and not two

RNA node bases come from the assembly, so DNA vs RNA is one relation over two node extents (S_g vs S_x). A contiguous
genomic window cut to the RNA node's length reproduces the DNA edge set (351 edges, symmetric difference 0, 91/91):
contiguity, not length, separates levels; the single-record requirement IS the contiguity requirement.
(full: archive/2026-11/seeded_family_definition.md:961)

### 1★.2 What graph class a family actually is — measured, and it is NOT one class

`o1_closure` panel (7 families, 80 nodes, CHM13, shipped tier): density 0.348-1.000 (median 0.864), λ 1-6; NPIP λ 5,
TBC1D3 2, MAGEA 6, GSTM 2 (complete), RABL2 1, APOBEC3 1, HERC2 1 (density 0.348). Retracts "families are complete
graphs" as the reason γ is inert; γ = 0.20 is still inert (margin 0.148); HERC2 would split at γ = 0.40. What ships is
the connected component; γ is a hairball floor load-bearing only blind (`ER_CC` 0.3000 below null 0.3408).
(full: archive/2026-11/seeded_family_definition.md:978)

### 1★.3 The DNA↔RNA relation, measured on the same nodes

|E_dna| / |E_rna| / shared 526 / 481 / 448; RNA-only 33 (6.9%); 1 of 33 joins two DNA components (APOBEC3); partition
equality 6/7; DNA edges retained 448/491 = 0.9124. Restates "RNA-only edges = 0" as panel-specific. Numbers on the
post-2026-08-11 rule (NPIP RNA 339, not stored 346). (full: archive/2026-11/seeded_family_definition.md:1006)

### 1★.4 ⚠ TWO CANDIDATE FRAMINGS TESTED AND REJECTED — record these, they look right

(a) "DNA quasi-clique, RNA components": REFUTED, inverted (RNA denser in 4 of 7 because RNA sheds nodes). (b) "RNA family
= k-core of DNA family": NOT ESTABLISHED; pooled AUC 0.7668 was a Simpson artefact; within-family sign test p = 0.8125.
Do not put it in the definition. (full: archive/2026-11/seeded_family_definition.md:1031)

### 1★.5 ⭐ WHAT THE GRAPH THEORY DOES BUY: λ AS A CERTIFICATE, NOT A DEFINITION

SHIPPED 2026-08-14: every emitted family carries edge connectivity λ (Stoer–Wagner, `family_split::edge_connectivity`):
at λ ≥ 2 no single record's loss splits the family. λ = 1 set {RABL2, APOBEC3, HERC2} = every historically fragile
family. Not a membership criterion: a 2-copy family has λ = 1 necessarily (`cut_certified = false` is not a defect flag;
test `certificate_two_copy_family_is_always_lambda_one`). λ says nothing about biological reality, missing copies (O3)
or assignability (O2). (full: archive/2026-11/seeded_family_definition.md:1049)

### 1★.6 What the restatement concedes — state these first

1. RNA cannot be shown to add edges in general (1 of 481); the joint object is a property (§P4c). 2. The node set is an
input (245/245 blind lost edges fail COVERAGE). 3. Reach bounded by divergence: 22/40 = 0.5500 [0.3983, 0.6929] on chr1
(Fisher p = 0.6090; §4a). 4. Families are not one graph class. (full: archive/2026-11/seeded_family_definition.md:1086)

## 1. Definition

(Normative design definition of 2026-08-09 for the O1 catalog; the shipped de novo default today is §0★.)
(full: archive/2026-11/seeded_family_definition.md:1099)

### ⭐⭐ 1.0 THE OPERATING ASSUMPTION (stated 2026-08-12, `/home/juanfra/winloci_scratch/o1_minann/O1_MINIMAL_ANNOTATION.md`)

> **O1 assumes the genome, a minimal annotation, and reads aligned to that genome.** "Minimal" is bounded on three
> channels: INTERVAL, SYMBOL, FAMILY. **O1's task is to supply the copy relation, for members the annotation never
> named.** It is not to discover loci ab initio.

- Granted: M1 coverage divides by the shorter sequence, so an edge needs ONE well-delineated endpoint (122/122; 4 exact
  nodes of 27 recover 79 of 80 pairs).
- CONCEDED permanently: the fully annotation-free claim (row O1.13). Detection was fine (0/27 loci absent; chr1 capture
  120/162 = 0.7407); delineation bound (blind nodes median 5.97x too long; 245/245 lost edges fail COVERAGE; 152/152
  reproduced).
- Annotation requirement (chr1, 40 families): interval completeness gain present at c = 1.00, gone by c = 0.50, +0.0000
  at c ≤ 0.25 (knee in (0.50, 1.00], not located, OPEN); boundary error ±2 kb costs nothing, cliff at 5 kb (false-merge
  0.0806 → 0.5000, Fisher p = 2.7e-11). The rule is STARVED, not broken (pairs 277 → 58 → 8 → 2, fire rate 0.12-0.18).
- Does not buy: τ = 0.50 headline gain +0.1000 [+0.0000, +0.2250], P = 0.0635, NOT significant; surviving claim m ≥ 3
  stratum 5/22 vs 0/22; dispersed families 0/8 (rule never fired); below ~0.75 identity = all-singletons floor; chr1 not
  genome-wide; robustness to a WRONG annotation untested, OPEN. (full: archive/2026-11/seeded_family_definition.md:1101)

**The definition (rewritten 2026-08-09).** The edge rule has one operative constant, `c`. Fix a genome `G` and a tier `T`.

> **Node set.** `V` = loci proposed from the data alone (read-supported loci from `pass1_skeletons` → `assemble_gate`,
> or under `--from-genome` duplicated genomic loci). `V` depends on no seed and no annotation.
>
> **Substrate `S`** (a term of the definition, added 2026-08-13): `S_x(v)` exon-sum or `S_g(v)` genomic span.
> `E_r(S_x)` and `E_r(S_g)` are different relations (244 human reps: |E_x| = 696, |E_g| = 1575, shared 630, ARI 0.4702);
> numbers under different substrates never sit side by side.
>
> **Support relation.** `E_r` = pairs `{u,v}` for which **some single alignment record** between `S(u)` and `S(v)` at
> tier `T` covers at least `c` of the shorter. (Code also checks identity ≥ τ; it behaves as `cov ≥ c`.)
>
> **Family.** Blocks of the γ-quasi-clique refinement of the connected components of `(V, E_r)`; the blocks PARTITION
> `V` (`family_split.rs:477`, only splits).
>
> **Seed — a QUERY, not a term.** `F(s)` = the block containing `s`; `F(s') = F(s)` for every `s' ∈ F(s)`. P1 is a
> theorem about a partition.
>
> **Certificate — emitted, never decided (2026-08-14).** Each family is emitted as `(B, λ(B))`, edge connectivity on
> the emitted copies; `λ ≥ 2` = no single record's loss splits it. Never a membership term (2-copy families have λ = 1).

Emission ≠ definition: only blocks with ≥ `--min-copies` (default 2) distinct loci are written. `<out>.families.tsv`
carries `n_edges`, `density`, `lambda`, `cut_certified` (λ computed after `distinct_locus_reps`; multi-record pairs are
one edge; `NA`, never `0`, on paths without an `E_r` graph and under `--refine`).

| symbol | shipped | measured role |
|---|---|---|
| `c` | 0.50 of the SHORTER locus | OPERATIVE, the only one (plateau claim later REFUTED, see below) |
| `S` | O1 catalog `S_x` (`homology_genomic_span` OFF); `refine` `S_x ∪ S_g` (gated, family-local) | OPERATIVE; the axis on which the two call sites differ |
| `T` | `minimap2 -c -X --no-long-join -k11 -w5` | fixed input, now defined once |
| `τ` | 0.60 (`sensitive_identity`) | INERT on [0, 0.65] |
| `γ` | 0.20 (`RUSTLE_GENOME_GAMMA`) | INERT on seeded real families (0.148 below sparsest, HERC2 0.348); load-bearing blind |
| `β`, `d` | 5 kb, 10 kb | not in the shipped pipeline (seeded probe's node leg); operative where they exist |

### Why each demoted constant is demoted (and what remains free)

- **τ = 0.60 INERT**: 634,502-record census min identity 0.6430 (0.6582 above the coverage floor); τ = 0 gives identical
  edge sets; first bite at 0.70; 99 adversarial pairs never below 0.6469. Divergence exits through coverage, not
  identity. Empirical, not a theorem (negative-scoring secondaries exist; every AS > 0 record has identity ≥ 0.6822).
- **γ = 0.20 INERT**: can only bind on |C| ≥ 3 below γ; β/d panel min density 0.4722 (panel-only; do not cross-quote
  with 0.348). Inert is not harmless: GOLGA8 blob (0.210) rejected, SPDYE+PMS2P (0.660) admitted. 0.20 ships (other docs
  said 0.40 / 0.30 / γ_core 0.13).
- **β, d NOT DEMOTED, operative node-construction choices**: β 0 → 30 kb gives 73 → 45 nodes (−38%), deletes read-
  supported members (MAGEA1, GSTM3 ...); APOBEC3A is a cut vertex (1 of 35 cells). Not in the shipped definition; the
  shipped node leg's constants (RNA `pass1_min_reads = 2`, `min_terminal_support = 2`, `GATE_MIN_READS = 3`; DNA
  `min_identity = 0.90`, `min_block = 1,000 bp`, `max_locus_span = 3 Mb`) must be quoted with any node set. WITHDRAWN
  2026-08-11: blind purity numbers 0.136 / 0.237 not quotable.
- **c = 0.50 OPERATIVE**: sweep of 27 values on the 25-region gorilla panel + 15/15 paralog panel; pre-M1 byte-identical
  interval [0.50, 0.60]; below it families are invented (0.10: 7/7 single-copy controls emit), above it members lost
  (0.65 GSTM). **CORRECTION 2026-08-11**: post-M1 the interval is [0.48, 0.60]; on 35 gene-tight windows the widest
  interval containing 0.50 is (0.4721, 0.5000], 0.50 on the upper boundary; row 1.7 records the plateau sweep as
  **REFUTED**. Defensible, not optimal (0.50 recall 6/15 with 1 false merge; 0.60 4/15 with 0). 0 of 1,383 cross-family
  edges at the operating point; absolute-bp floors refuted (5,000 / 2,000 / 500 bp per substrate). Scale-dependent (§1a).
  M1 coverage form fixed 2026-08-10 (+0.22 to +2.29 bite points, 1 region of 25 changes). Do not change `c`.
- **T** unified 2026-08-10 (`er_tier_flags_appear_as_a_literal_exactly_once_in_the_source`); panel minimum 0.4722 and
  the APOBEC3 RNA edge exist only under the panel invocation.
- Other code constants: panel node constructor `blocklen ≥ 1,000 bp`, `nmatch/blocklen ≥ 0.80`, 2 kb span floor;
  `RUSTLE_COVERAGE_SPLIT` default 0.0 (off). (full: archive/2026-11/seeded_family_definition.md:1244)

### The seed: what changed, and what `--seed` now is

The old seeded `V(s)` FAILS strict seed-invariance on human NPIP (19 seeds → 4 distinct `F(s)`, 64 of 171 pairs agree,
37.4%); `V(s)` inherits the seed's length; the seed never acts through the edge rule (contribution 0); seeding from
single-copy SDHA yields a false 14-copy family. Hence `V` is seed-free and P1 is a theorem. `--seed <chrom:start-end>`
on `gw_family_catalog` is a **projection over the emitted catalog** (byte-identical run plus `<out>.seed.tsv`; reports
the block by max overlap; `ABSTAIN_NO_OVERLAP`; use `--min-copies 1` to separate singleton from no locus;
`src/rustle/vg_family/seed_projection.rs`, `tests/gw_family_catalog_seed_projection.rs`). RABL2 demo: two seeds on
different members return the identical 5-member component. Withdrawn: the `bench/crossspecies` seeded probe computes
the same object (MAGEA4 1 vs 9; GSTM3 0 vs 4; HERC2 164 vs 2). asm20 tier retired 2026-08-07 (`RUSTLE_ER_SENSITIVE_ONLY`
default true; 0 unique edges on genomic, 1 on spliced). One definition since 2026-08-09: `refine_families_exon_sum` is
OPT-IN (`--refine`) on the homology catalog, default only for `--window-catalog` / `--cross-chrom`; flip cost zero at
shipped defaults (0 of 50 files), conditional on defaults. Annotation supplies `s` only.
(full: archive/2026-11/seeded_family_definition.md:1446)

### 1a. ⚠⚠ `c` IS SCALE-DEPENDENT — the edge SUBSTRATE is a free parameter of the definition

- **§1.1 What the binary emits (O-2, 2026-08-13)**: O1 catalog `homology_edges_all_reps_pooled` = `S_x` exon-sum,
  global swap (`homology_genomic_span` OFF); `refine_families_exon_sum` (`copy_assign`, conflict catalogs) = `S_x ∪ S_g`,
  gated family-local addition only for copysets the exon-sum core leaves disconnected. The shipped O1 `E_r` is the
  exon-sum ("RNA") relation. Genomic span stays OFF on unresolved sign (recall 0.1490 → 0.1780; precision 0.8501 →
  0.9876 not robust: 40.7% up / 39.4% down). Union fires only on co-location windows (`--cross-chrom` 288 edges added,
  11/26 partitions moved) and lands on the DNA-only partition (0/26 differ). Load-bearing case: gorilla SHARP `CAFAM0`
  (exon-sum coverage 0.0480 → 0 edges; genomic 1.0000 → 1 edge; 339 reads). Identity never fails; coverage is the only
  failure mode. `E_r` not invariant to node naming (8 of 260 = 3.1%).
- `c = 0.50 of min` is a completeness filter, not scale-invariant. 61-node human panel, pair-bite (SHIPPED, M1-fixed):
  genomic span 8.52%, pooled exons 28.83%, single rep transcript 31.33% (PANEL 4.47% / 7.57% / 34.36%); 0 cross-family
  edges everywhere. Median 0.50·min demand 6,338 / 3,548 / 876 bp. 28.31% vs 5.71% were the two tiers; quotable gap
  3.81x. Second 59-node panel M1: pooled exons 7.62% → 26.62% (3.49x, not 4.66x); quote levels, not ratios; "panel
  series not monotone" withdrawn.
- It is the LENGTH, not the splicing (truncated genomic bites harder, 61.70%). A coverage fix reaches at most 56% (81 of
  144 lost edges die on coverage, 63 have no record). Headroom 0.30 (genomic) vs 0.10 (rep). Only pair-level bite may be
  quoted. 5 of 7 families change partition with substrate.
- Recordable: `RUSTLE_ER_EDGE_DUMP=<prefix>` → `params.tsv` (`substrate_median_len_bp`,
  `coverage_floor_median_bp_demand`, `coverage_form`, derived `mm_args_sensitive`). MAGEA example: 875 vs 2,017 bp
  (2.31x) from the flag alone. (full: archive/2026-11/seeded_family_definition.md:1541)

#### Why this is documented rather than changed

- Genomic span as default: residual moves into locus boundaries (10 of 56 RNA loci span > 1 DNA copy), mixed gorilla
  result (HERC2 gain, GSTM loses 2, MAGEA fragments), ~450x cost. Not shipped.
- Absolute aligned-base floor: **REFUTED, do not re-propose** (safe N 5,000 / 2,000 / 500 bp, 10x spread).
- Union-of-records coverage `RUSTLE_ER_SUM_COVERAGE=1`: **REMOVED FROM THE CODE 2026-09-24** (r1097; ≥ 90% band
  0.903 → 0.806); recover from tag `notebook-2026-09-24`.
- Documenting first removed an undeclared free parameter without changing any number (D1 closed 2026-08-09).
(full: archive/2026-11/seeded_family_definition.md:1725)

## 2. Properties

(full: archive/2026-11/seeded_family_definition.md:1773)

### P1 — Seed-invariance (the anti-arbitrariness property)

> **Claim.** For all `s' ∈ F(s)`, `F(s') = F(s)`.

Since the §1 rewrite P1 is a THEOREM about the partition of a seed-free `V`. The panel below measures the retired
seeded probe (tier not applicable: seed→genome pass only). Seeded probe: **HOLDS ON MEMBERSHIP; boundaries are not
invariant**: 19/19 NPIP seeds recover all 19 genes and the identical gene set; |F(s)| 26 or 27; interval Jaccard median
0.963, min 0.256 (NPIPB8, shortest seed, boundary effect). Membership is seed-invariant; locus extent inherits the
seed's extent. (full: archive/2026-11/seeded_family_definition.md:1775)

### P2 — γ-independence

> **Claim.** If `F(s)` has induced density `d`, then `F(s)` is a γ-quasi-clique for every `γ ≤ d`.

Shipped tier 2026-08-11: human NPIP 27 nodes, 351 edges, d 1.000 (5.00x margin at γ 0.20; complete graph, so vacuous);
gorilla NPIP 25 nodes, 260 edges, d 0.867 (4.33x), one edge lost to a 0.0011 coverage near-miss (0.4989). Node sets are
tier-invariant by construction. (full: archive/2026-11/seeded_family_definition.md:1827)

#### Recomputed 2026-08-10 at BOTH tiers, on the §1a 61-node panel (`o1fix_o2audit/fix/s1a/`)

Genomic span: shipped d NPIP 0.926, MAGEA 0.939, HERC2 0.700 (min margin at γ 0.20 3.50x); partition differs by tier on
APOBEC3 and GSTM. Rep transcript: NPIP 0.569 (comps 24,1,1; γ 0.40 margin 1.42x shipped vs 1.09x panel), HERC2 no
family. (full: archive/2026-11/seeded_family_definition.md:1870)

#### The component-density DISTRIBUTION, both tiers, three substrates (2026-08-11, `close_o1o2/dens_both.py`)

20 family x substrate cells (`o1_closure`, human): **γ = 0.20 inert under both tiers and both coverage forms in all 20
cells.** γ = 0.40 binds on one component per tier, not the same one (PANEL HERC2 genomic d 0.3485; SHIPPED HERC2
pooled exons 0.3333 M1): a γ = 0.40 verdict is not tier-portable. Tier changes partition on 4 of 20 cells.
(full: archive/2026-11/seeded_family_definition.md:1898)

#### P2 on families OTHER than NPIP, at the shipped tier (2026-08-11, `multifamily_p1p4.py`)

Shipped d: PCDHB 0.6684, TUBA 1.0000, TBC1D3 0.6511, GOLGA 0.6106 (panel 0.3078: second γ = 0.40 flip, opposite
direction; edge Jaccard 0.4831); NBPF OPEN (do not fill from panel archive). TBC1D3/GOLGA node counts inflated by a
seed-set `startswith` defect (TBC1D30/31/32, antisense records), not a definition failure. Conclusions: γ 0.20 margin
≥ 2.84x everywhere; γ 0.40 margins tier-dependent; tier changes partition on 2 of 7 families.
(full: archive/2026-11/seeded_family_definition.md:1939)

### P3 — Non-circular recovery

> **Claim.** Members recovered were not supplied.

Tier-invariant by construction. Human: 1 seed in, 19/19 annotated NPIP genes out, one component of 27 (8 loci overlap no
NPIP gene). Gorilla: the annotation names exactly one NPIP gene; 25 loci (36 of 37 overlapping symbols are `LOC…`);
validated against human annotation: 23/25 loci with a qualifying edge (median identity 0.956), all 19 human NPIP genes
hit. (full: archive/2026-11/seeded_family_definition.md:1985)

### P4 — Substrate convergence

> **Claim.** `F(s)` computed from the genomic substrate and from the RNA substrate agree.

Probe→genome pass (tier n/a), human: RefSeq NPIPB11 transcript 9/19 genes, 21 loci; IsoSeq read-supported exons 19/19,
27 loci; genomic span 19/19, 27 loci. RNA and the T2T genome converge; annotation is the outlier.
(full: archive/2026-11/seeded_family_definition.md:2004)

### P4b — DNA graph vs spliced-RNA graph, same node set

Original table (RNA 338/351, Jaccard 0.963, 0 RNA-only) was PANEL tier with coverage aggregated over records, NOT the
shipped rule (relabelled 2026-08-10). (full: archive/2026-11/seeded_family_definition.md:2025)

#### The 27-node set AT THE SHIPPED TIER + SHIPPED single-record rule (2026-08-11) — this row is now CLOSED

DNA 351, RNA 260, shared 260, DNA-only 91, RNA-only 0, Jaccard 0.741; one component either way. 91/91 DNA-only are
coverage misses at identity 0.9709-1.0000 (max coverage 0.4989). Quote 0.741, not 0.963. Containment is not panel-
robust (80-node `o1_closure`: 33 RNA-only of 481); quote "only 1 of 33 joins two DNA components" (§1★.3). Diagnostic
defect fixed: earlier "no alignment record at all" lines from `graph_vs_graph_report.py` are void.
(full: archive/2026-11/seeded_family_definition.md:2054)

#### The same comparison at the SHIPPED tier + SHIPPED single-record rule (2026-08-10)

61-node panel totals: DNA 408, RNA 263, shared 250, DNA-only 158, RNA-only 13 (NPIP 7, MAGEA 4), Jaccard 0.594,
partition same 3/7. A qualification of P4b ("zero RNA-only under every rule tested" applied only to 2026-08 rules).
Aggregated coverage adds DNA edges about as fast as RNA edges at the shipped tier.
(full: archive/2026-11/seeded_family_definition.md:2091)

### P4c — IS THE DNA/RNA JOINT OBJECT DEFINITIONAL? **NO — IT IS A PROPERTY** (closed 2026-08-13)

Verdict `/home/juanfra/winloci_scratch/o1_joint/O1_JOINT_VERDICT.md`; tested definition in
`docs/o1_investigations.md#the-joint-dna-rna-family-definition-retracted` (§9). **Report the DNA/RNA relation as a
PROPERTY of this definition, never as a clause of it.** Joint partition = DNA-only on 0 of 5 disagreement families
(an entailment, re-quoted by O-2: under the shipped gated-union rule 2/5 and 2/7 at rep-transcript substrate); κ label
Fisher p = 0.40 (n = 6 cannot reach 0.05); RNA-only nodes 24 → 8 of 362 Soto members at the shipped tier (quote 8;
predicted direction of the tier correction REFUTED: DNA 351 → 349). Ceiling 13 of 362 (3.6%), all 13 with a DNA
paralog at 1 − de ≥ 0.9703: "loci RNA sees and DNA cannot" is EMPTY. RETRACTED, do not quote: the 2026-07-25 "DNA 90.3% /
RNA 84.5% / union 97%" (banned `nmatch/blocklen` metric). RNA column of the re-derivation is INHERITED. 27-node
redundancy null has power zero (complete graph); on 61-node NPIP RNA's losses are not random (p ≈ 2 × 10⁻⁴; NPIPB5 RNA
degree 0/25). The 91 lost edges are about CONTIGUITY: a contiguous genomic window of the RNA node's length recovers
91/91 (symdiff 0). An `E_r` edge needs at least ONE endpoint contiguous over the shared region.
Spec defect: the binary uses three substrates (refine `E_x ∪ E_g` additive tier default ON; catalog genomic span
default OFF); O-4 kept both defaults (catalog flip rejected on unresolved sign; contamination sign flips with truth).
O-3 fixed the certificate (`core_substrate`, `additive_genomic_tier`, refine DATA rows incl. `n_edges_genomic_tier_added`;
row X.2 re-quoted: diff is one line). `O1 ⊥ O2` is false at the node on the default catalog path
(`distinct_locus_reps` uses `reads_distinguish` and `conflict.min_reads`); `E_r` itself is clean.
(full: archive/2026-11/seeded_family_definition.md:2119)

### ⭐ RESOLVED 2026-08-13 — SCOPED, MEASURED AND MONITORED, NOT REMOVED

O1 decides membership by **SEQUENCE alone, EXCEPT** for a co-located same-strand pair of reps the junction-asymmetry
rule did not settle (true K = 0 case, junction-less overlapping reps), where the unique-mapper count decides. Not removed
(re-introduces advisor-flagged over-merges); not to be "fixed" by intron-chain comparison (isoforms). Measured inert
(`RUSTLE_LOCUS_AUDIT=1`): read leg decided 0 of 15 / 0 of 94 pairs, 0 firings on chr1. Monitored by an unconditional
`[o1-perp-o2] WARNING` line. Guard `homology_catalog_never_touches_the_conflict_graph` now derived from imports,
transitive, `DISCLOSED_EC_USES` ledger (adding to it is a thesis edit). "The RNA graph" is not one object (rep transcript
263 vs pooled exons 316 edges; 4 of 7 partitions differ); only "intersection damages HERC2" survives both. Joint object:
necessary to state O3, inert for O2 (K ± 1 moves 0 of 7,691 reads). What to adopt is a SUBSTRATE decision (+176 true /
+13 false edges at id ≥ 0.90, precision 0.908 → 0.916 on 540 reps). Do not quote RNA density 0.547 / 0.698.
(full: archive/2026-11/seeded_family_definition.md:2288)

### Variation graph: RETRACTED (2026-08-08)

The "single ~16.1 kb minigraph path" was the reference echoed back (`inserted 0 events`); the path/sub-path account of
`E_RNA ⊂ E_DNA` is withdrawn. The cassette result stands on the pairwise block analysis. Intron bubbles need a
base-level `vg construct`. (full: archive/2026-11/seeded_family_definition.md:2363)

### On lowering the identity floor below 0.80

The sensitive tier is already 0.60. Within-family density sweeps (0.80 → 0.698 ... 0.50 → 0.997) cannot justify a lower
global floor (tautology: no negatives); genome-wide id 0.70 gave one 776-copy hairball.
(full: archive/2026-11/seeded_family_definition.md:2381)

### None of P1–P4 uses read conflict — on purpose

> **`E_c` is a held-out test of `E_r`, not a term in it.**

Dense `E_c` inside a family confirms; dense `E_c` between families warns of a split (NBPF–NOTCH2NL, NPIP–SULT1A).
Feeding `E_c` into `E_r` would be training on the test set (O1 ⟂ O2 as measurement discipline).
(full: archive/2026-11/seeded_family_definition.md:2392)

## 3. The certificate (and why it is not the definition)

`E_c`, the read-conflict relation (pairs of loci a read cannot be placed between), is reported **alongside** the
family, never inside it.

> ⚠ **`E_c` must never enter `E_r`.** O1 (identify, support) and O2 (assign, ambiguity) are independent by
> construction; using read conflict to build families destroys the K=0 argument.

Its honest form is **abstention**. Human NPIP, 21,775 reads: 64.6% on > 1 NPIP gene, all 171/171 gene pairs linked.
Specificity on `human_testis.t2t.bam` (10 families, 135 genes, 1,312 within-family pairs):

| shared-read threshold | within-family linked | between-family linked | enrichment |
|---|---|---|---|
| ≥1 | 72.2% | 1.1% | **67×** |
| ≥3 | 69.6% | 0.4% | **179×** |
| ≥10 | 63.2% | 0.2% | **326×** |

Per family at ≥ 3: NPIP 171/171 (replicated in A119b and testis), NBPF 91/91, PCDHB 80%, TBC1D3 42%. Every
between-family link is a known co-duplication (NBPF–NOTCH2NL, NPIP–SULT1A). Within-family coverage is a lower bound
set by expression. (full: archive/2026-11/seeded_family_definition.md:2410)

### Certificate claims C1–C3 (properties of the EVIDENCE, not of the definition)

- **C1 completeness**: fraction of member pairs linked by a shared read (NPIP 171/171); a lower bound, never a failure
  claim alone.
- **C2 specificity**: within-family linkage must exceed between-family by a margin that grows with the threshold (67x /
  179x / 326x).
- **C3 abstention semantics**: a link records that the aligner could not separate two loci; an abstention, not a merge
  instruction. (full: archive/2026-11/seeded_family_definition.md:2457)

## 3a. THE THREE USES OF MULTIMAPPING — keep them apart

| # | use | object | side | role |
|---|---|---|---|---|
| **N0** | **presence**, `sec_frac(locus)` | per-locus **scalar** | **O1** | proposes candidate NODES |
| **C** | **conflict**, shared reads between loci | pairwise **relation** on a built family | held out | certifies the family |
| **K** | **copy number**, χ(H) of the conflict graph | pairwise relation, then a colouring | **O2** | lower-bounds the number of copies |

> **The one-line rule: a scalar over loci may build the graph; a relation between loci may only test it.**

K (`max(n_loci, χ(H)) ≤ true copy number`) and N0 are never one quantity. (full: archive/2026-11/seeded_family_definition.md:2471)

## 3b. N0 — Multimapper PRESENCE as a node signal (OFFICIAL STAGE)

Promoted to a named stage 2026-08-07: the only annotation-free node signal measured to work. Purity 0.136 / 0.237
WITHDRAWN 2026-08-11.

> **`sec_frac(locus)` = secondary records / (primary + secondary) records overlapping the locus.**

Stage contract: INPUT an aligned read set (no annotation, seed, catalog); OUTPUT ranked candidate intervals (candidates,
not loci; `E_r` filters next); GUARANTEE nothing pairwise is emitted. (full: archive/2026-11/seeded_family_definition.md:2492)

### ⭐⭐⭐ 3b-0. THE BLIND / SEED-FREE ROUTE, MEASURED END TO END (2026-08-12) — read this before the N0 table below

Human, `A119b.t2t.bam`, shipped tier; verdict `/home/juanfra/winloci_scratch/o1_blind/O1_BLIND_VERDICT.md`.
**Status: OPEN, with a measured ceiling and a named mechanism.**

> **Without a seed the method still FINDS the loci but cannot DELINEATE them, and the failure goes through exactly one
> clause of `E_r`.**

- Finding: 0/27 loci absent (one gene family; no TBC1D3: the "NPIP + TBC1D3" scope line is withdrawn); genome-wide
  38,130 nodes; chr1 gene_capture 120/162 = 0.7407. Delineation: median 5.97x too long, 14/27 FUSED.
- Mechanism: 6/251 seeded edges survive; 245/245 losses = identity passes, coverage < 0.50 (third substrate with that
  signature). TRIM 250/251 is an oracle ceiling; FILL is a tautology check.
- chr1 score (scope achieved chr1 only; shipped `-k 11` tier is NOT SCOPE-INVARIANT, `mid_occ` from the whole index):
  SEEDED 22/40 = 0.5500 (m ≥ 3: 5/22); BLIND ER_QC 18/40 = 0.4500 (m ≥ 3: 1/22, ARI 0.0467); BLIND ER_CC 12/40 = 0.3000
  below its null 0.3408; singletons 17/40. γ = 0.20 binds blind (contamination 0.7473 → 0.0976), inert seeded.
- Tried and failed: read-through filter REFUTED (0 of 6,824 links removed; read min-cut median 12, max 273); mis-chain
  gate does not clear null (p = 0.1534); depth split 178/251 NOT QUOTABLE (constant picked against truth; random-cut
  control 152/251; nulls must be size-matched).
- Works: junction-chain split (no new constant, `MIN_READS = 3`) 119/251 = 0.4741, p = 0.0011, FUSED 14/27 → 0/27,
  extent 5.97x → 1.89x. Oracle ladder: perfect split 29/251 (fails null), boundary only 170/251: **SPLITTING IS NOT THE
  LEVER; THE BOUNDARY IS, to about ±1 kb.** One well-delineated endpoint suffices (79 of 80 pairs).
- Exposures 1-10: panel truth is the seeded pipeline's own output (concordance, not accuracy); scope windows truncate
  reads (14.68%); blind survival significantly BELOW chance (p(≥ 6) = 0.9644); denominator conditioned on output; no
  precision term; chr1 clustering caveat (later withdrawn by §4a.1); min()-coverage forgives fragmentation REFUTED but
  survival non-monotone in extent; "27 loci in one component" useless (39 nodes); human only; 0.136 / 0.237 stay withdrawn.
- To close: K-ladder chr1..chr3 all-vs-all (~1.2 h + ~2.5 h); re-score panel against GFF; a seed-free BOUNDARY rule.
  Will not close it: node filters, read-dropping filters, `c` sweeps. (full: archive/2026-11/seeded_family_definition.md:2502)

### N0 measured, seed-free and annotation-free (5 chromosomes, fixed parameters)

**Third tier, human, OPEN**: scripts ran `-x asm20 -k11 -w5 -c --eqx -N 200 -p 0.02` with `nmatch/blocklen ≥ 0.60`.
Scripts fixed 2026-08-11 (`rustlib.py`, `*.er.paf`, `*.partition.tsv`); table superseded, not refreshed (cells are
oracle-selected `max(comps)`). On byte-identical `sec_frac` node sets the aligner moves ~13x more edges than the metric
(+446 / −129), partition changes on 6/6 chromosomes; all 129 tier-lost edges lost on COVERAGE. Panel-tier table: purity
0.279 (GOLGA8) to 0.857 (GOLGA6L); SPDYE 15/17 at 0.833; NPIP 18/19 at 0.556; families come out as components (chr16
157 candidates → 78 components). (full: archive/2026-11/seeded_family_definition.md:2664)

### N0's two failure modes — one fixable, one not

(1) Fixed candidate budget (engineering): chr2 got 97 candidates, ID_65 0/7 despite `sec_frac` 0.84-0.93; scale N to
duplicated content. (2) Co-duplicated families, not fixable by a threshold: chr7 SPDYE + PMS2P merged at density 0.660;
chr15 GOLGA8 blob at 0.210 rejected.

> ⟹ **N0 solves node DISCOVERY, not family SEPARATION.**

(full: archive/2026-11/seeded_family_definition.md:2732)

### What N0 does not do

Does not replace the seed (purity 0.56 seed-free vs 1.000 seeded on NPIP); does not measure copy number (K, O2); does
not define edges, and must never. (full: archive/2026-11/seeded_family_definition.md:2750)

## 3c. Multimapper PRESENCE — supporting measurements

`sec_frac` is a per-locus scalar ("is this locus duplicated"), a node question. chr16 A119b, 81 loci: duplicated median
0.94 (n 24), single-copy 0.05 (n 57); AUC plateau 0.94-0.98 over 11 of 15 labellings (0.9786 peak post hoc; the plateau is
the claim). SULT1A1 lowest positive (0.417) and no DNA paralog. Evidence for node selection, not a node builder; no cut
proposed. Multimapper WIDTH does not work (adjacent copies identical footprints; do not retry). MAPQ=0 is the wrong
proxy (evidence is in secondary records; minimap2 emits no NH tag). (full: archive/2026-11/seeded_family_definition.md:2756)

## 4. Known exposures (state these first, do not wait to be asked)

- PDXDC1 is a genuine false positive; human precision 19/27 = 0.70 to 25/27 = 0.93.
- Gorilla returns 25 NPIP loci vs published ~17. Unresolved.
- **Node extent, not homology, is the recurring failure** (coverage only; gorilla 207 kb tandem merge).
- Two-definition defect survives in `copy_assign --homology-primary` (refine runs unconditionally); not measured.
- `c` scale-dependent, substrate a free parameter (§1a; 6,338 vs 876 bp); unchanged by M1. M1 denominator defect CLOSED
  2026-08-10.
- τ and γ cannot bind; β/d are node constants (§1).
- `T` two invocations CLOSED 2026-08-10 (B1); backlog cleared 2026-08-11 (P2, P4b recomputed); §3b N0 table still third
  tier, OPEN.
- **Over-merge mode still UNTESTED (D4)**; next panel must be paralog pairs.
(full: archive/2026-11/seeded_family_definition.md:2800)

## 4a. ⭐⭐⭐ REACH — THE BOUND, AND WHAT SETS IT (2026-08-14)

Verdict `/home/juanfra/winloci_scratch/o1_dispersed/O1_DISPERSED_VERDICT.md`; human CHM13 chr1, RefSeq, shipped tier;
denominators W = 73 / WC = 460 / X = 12,508 fixed before any method.
- **THE BOUND**: O1 recovers **22 of 40 = 0.5500 [0.3983, 0.6929]** chr1 multi-copy families; chr1 representative by
  composition (16/40 vs 157/432, Fisher p = 0.6090).
- **Set by DIVERGENCE, not dispersal**: median protein identity ≥ 0.50 → 9/14 edge-bearing; below → 1/26. At matched
  divergence the coordinate stratum has no effect (p = 0.3714).
- Failure: dispersed families fire 0/73 [0.0000, 0.0500]; not aligner tuning (seed census 3/73 vs clustered 207/460;
  57/73 never get a record at any of 20 settings); aligned pairs span 126-312 bp, max coverage 0.0856. Homology exists
  (protein SW 70/73 at E < 1e-3; median protein identity 0.2809).
- Thesis sentence: O1 reaches recently duplicated families (~55%), not anciently duplicated ones (1 of 26 below 0.50
  protein identity; 0 of 73 dispersed pairs); boundary set by divergence, not moveable by any nucleotide aligner.
(full: archive/2026-11/seeded_family_definition.md:2851)

### 4a.1 What is WITHDRAWN or RESTATED by this pass (retractions travel together)

1. WITHDRAWN: "chr1 80.0% clustered vs 36.3%" (two different rules; matched 40.0% vs 36.3%, p = 0.6090); delete the
   "0.4500 biased upward" caveat. 2. RESTATE "0 of 8 dispersed" as "1 of 24 genome-wide-dispersed" (clustered 9/16 vs
   dispersed 1/24, p = 0.000333; makes O1 look worse, publish with 1). 3. "dispersed" → "anciently diverged" (window
   sweep 14 → 0). 4. No effective identity floor: 71 of 73 pairs produce no record. 5. Identity never the failure mode
   (seventh substrate); protein identity anti-informative (cross 0.3111 > within 0.2809). 6. Coarse truth key over-merges
   this stratum (fine dispersed denominator 27). 7. `div_strata` medians are genome-wide. (full: archive/2026-11/seeded_family_definition.md:2915)

### 4a.2 The three routes that were tried, and why the definition does not change

- **More sensitive nucleotide engine: DEAD.** Best `S09_k7w1_max` 14/73 with 1,172 cross-family, recovery 0.6500 →
  0.2500; `N_splice_k9w3` 48/73 with 40.7% of all cross pairs. Every gain is a gap-spanning artefact (M1′ takes W to 0).
- **Protein: DIAGNOSTIC, not definition.** 34/73 with 634 cross-family (95.6% OR↔OR), largest component 39 > 20, P3
  FAILS; hand-deleting OR↔OR may NOT re-score P3; coding-only; 29/73 domain-level; tautology exposure.
- **Shared coding-intron position + phase: clears all four bars but DEFERRED.** 5/73, ΔX_new 0, recovery 0.5500 →
  0.6250, dispersed-only ablation +0.0250; deferred (unswept representative policy fails P1; needs retained introns;
  annotation-dependent). Its +0.1000 gain WITHDRAWN. (full: archive/2026-11/seeded_family_definition.md:2954)

### 4a.3 Limits of §4a, before being asked

n = 8 families / 73 pairs (family instrument nearly powerless); no orthogonal truth (RefSeq description strings);
superfamilies, not duplication classes; human chr1 + RefSeq only; divergence explanation partly circular (non-circular
form p = 0.3714); chr1 divergence representativeness OPEN; "no nucleotide engine" bounded to this seed shape and 20
settings. (full: archive/2026-11/seeded_family_definition.md:3004)

## 5. Three false-positive filters that were tried and failed

1. Expected-copy-number prior per species: premise false, and an oracle.
2. Seeding with the mRNA instead of the genomic span: human recall 19/19 → 9/19.
3. Coding fraction as a per-locus certificate: runs backwards (annotated NPIP median 0.080 vs 0.443); "no coding
   sequence = flank" withdrawn. (full: archive/2026-11/seeded_family_definition.md:3022)

---

> **Moved documents (wave 4, 2026-09-23).** Files this document cites that were pruned from the working tree — `docs/o1_investigations.md`, `docs/OBJECTIVES_AND_VERIFICATION.md`, `docs/o3_missing_copy_evidence.md`, `docs/NUMBERS.md`, `docs/ONE_METHOD.md`, `docs/METHOD_PSEUDOCODE.md`, `docs/OPEN_ITEMS_2026-09-09.md`, `docs/superpowers/` — are at git tag `notebook-2026-09-23` (`git checkout notebook-2026-09-23 -- <path>`) and in `~/Desktop/Rustle_attic/2026-09-23/` (see its `MANIFEST.tsv`). The citations above are provenance and were left as written.
