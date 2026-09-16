# Layer order for NPIP and TBC1D3 (human T2T-CHM13): protein P, DNA catalog E1 ("D"), subfamily clades C, RNA operator EXPR

2026-09-16, **revised after audit** (see §12 Audit notes). Spec: `docs/superpowers/specs/2026-09-16-family-layer-order-design.md`
(SCOPE AMENDMENT; Definitions, incl. the AUDIT AMENDMENT). Results: `/mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3/integrate_slim/`
(abbreviated `IS/` below). The pre-audit report, tables and scripts are kept in `IS/pre_audit_snapshot/`. Scripts:
`bench/layer_order/lo_*.py`. Not committed.

"D" in this report is the **§6kl RefSeq E1 guided catalog** (identity ≥ 0.70, coverage ≥ 0.30, MCL), **not** §0★★ clause
2/4. On chr16/19/20, the exact catalog used for NPIP, E1 **failed its pre-registered HGNC guard** (§6kl: 935 → 476,
0.509). The NPIP/PKD1 boundary is annotation-dependent: GENCODE E1, built the same way, merges PKD1 into the NPIP group (§6).

## 0. Summary

| question | result (NPIP + TBC1D3 only; descriptive) | file |
|---|---|---|
| order P vs D (E1) | **Not determined.** Under the spec's P-universe rule the two tie: pooled 246/255 both ways (group level 1/2 vs 1/2); NPIP 210/210 both ways, which is **one group each** (1/1); TBC1D3 36/45 both ways (0/1). On the 28 genes that lie in every layer's universe, P and D are the **same partition** (207/207 both ways). The tie moves with the conventions: of 9 variants beyond the spec reading, 6 give D > P, 2 give P > D, 1 ties (§6). | `IS/tournament.tsv`, `IS/pd_variants.tsv` |
| C vs P and D | **By construction**, C lies below P and D: clades sit inside one family, and P and D each keep every family in one group. The measured content is narrower: **P and D split 0 clade pairs**, for the literature clades and for the §0★★ clause-5 split system (C_tree). | `IS/disagree_clades.tsv` |
| TBC1D3 and C | Literature C has **0** TBC1D3 pairs, so every literature-C comparison is **NA**, not "above". The clause-5 split system does have TBC1D3 subfamilies, {B,F,G,H} and {TBC1D3,D,E,K}; both lie inside P's and D's TBC1D3 group (12 of 12 pairs). | `IS/tournament.tsv` |
| order inside C | Literature levels C_L1 ⊇ C_mid ⊇ C_fine hold **by construction**, since C_mid := C_mid_ab ∨ C_fine. **As built, C_fine is above C_mid_ab** (0.267 [15] vs 0.250 [16]). C_tree_top ⊇ C_tree_min holds by construction. | `IS/tournament.tsv` |
| cost of forcing P and D to nest (TBC1D3, HGNC, bipartite F) | JOIN, U-restricted: 0.917. JOIN, whole catalog groups: 0.857. REFINE with USP6NL attached: 0.957. REFINE with USP6NL as a singleton: 0.909. **No JOIN-vs-REFINE ranking holds under every convention.** Each side turns on one gene: DHX40 for JOIN, USP6NL/TBC1D26 for REFINE. On Soto (TBC1D3), REFINE scores 1.000 under both conventions and JOIN 0.818 or 0.643. | `IS/truth_agreement.tsv` |
| truth recall, layer-independent | HGNC group 2227 ("TBC1 domain family members") is **superfamily-level**: 57 symbols, 54 RefSeq genes. P recovers 54 of 423 member-anchored pairs (**R 0.128**, P 1.000). D recovers 45 of 198 (R 0.227, P 0.714). For Soto, every RefSeq gene that maps into a member-holding family is already in U, so the U-based Soto scores are not inflated by the closure. | `IS/truth_member_anchored.tsv` |
| EXPR, any-overlap reads ≥ 3 (spec rule, primary) | NPIP: **13 of 27** members expressed. TBC1D3: **0 of 19**. Unique reads give 9; unique ignoring the 6 readthroughs that overlap a member gives 11. Sweep for t = 1..5 in §8. | `IS/member_expression.tsv`, `IS/expr_sweep.tsv` |
| EXPR nesting and T2 | EXPR(L) ⊆ L holds **by construction**. T2 needs E_M ⊆ E_L inside each L group; on these data that precondition **holds for all 30 refinement pairs**, so the 0 T2 violations are **guaranteed**, not evidence. The only EXPR splits occur at reads ≥ 1 in the unique modes, where D MCL4 splits into {DHX40, RNFT1-DT} and TBC1D3 copies. | `IS/expr_T2.tsv`, `IS/analysis.out` |

## 1. Scope

- **Descriptive study of two development families.** Both were used to build the DNA rules (§6js). This is not a claim about layer order in general, and it was not pre-registered.
- **Layers:**
  - P.
  - D (= E1).
  - C: the clause-5 split system C_tree is the ordered C layer. The literature-anchored clades are a CIRCULAR reference.
  - EXPR, applied to each layer.
- **Deferred:** the SD layers S1, S2 and S3 (`light/S1.*`, `heavy/S2.*`; user scope, 2026-09-16 15:03). They were not read.
- **Substrate:** RefSeq `chm13v2.0_RefSeq_full.gff.gz` and the testis Iso-Seq BAM `human_testis.t2t.bam`.

## 2. Layers as built

| layer | rule as run | provenance | universe inside U | groups holding members |
|---|---|---|---|---|
| **P** | §6ko (`bench/protein_families.py`): longest CDS, r2 biotype filter, ≥ 10 aa, blastp e ≤ 1e-5 against a 20,088-protein genome-wide database (§6ko used chromosome-trio databases), coverage ≥ 0.30 of the longer protein, weight = identity × coverage, MCL I = 2.8. MCL runs on the neighbourhood N_k. The builder's **stopping rule** accepted k = 2: the first k in {2, 3} whose member clusters equal those at k + 1. | `light/scripts/README.md`, reconstructed. Rounds 1-5 ran with a one-shot `layer_protein.py round N` that was later removed; its argv was not saved. Stability: `lo_p_variants.py` → `IS/P_stability_plain.tsv` | **34**: every U gene with a protein, incl. PKD1 and DHX40 as `P\|other` | PC1 = 21 NPIP. PC2 = 9 TBC1D3 + TBC1D26 + USP6NL |
| **D = E1** | `mcl_families --min-exonic-bp 1 --min-shared-exon-frac 0.0 --dump-graph`: identity ≥ 0.70, aligned ≥ 0.30 of the longer gene's exon union, ≥ 300 bp, ≥ 1 exonic bp on both sides, overlapping loci folded (a folded record **inherits** the cluster of its locus), MCL 2.8 / 1e-9 | c15_17_22: argv in `light/work/D/c15_17_22.e1.err`; `clusters.tsv` md5 aef97aa0 = `lit/aj_dev/refseq_e1`. c16_19_20: md5 8e753e3c = `lit/aj_ho/refseq/e1`, the catalog that failed the HGNC guard (§6kl) | **63** (chr15/17/22 and chr16/19/20 only) | MCL0 (30 genes, NPIP), MCL7 (10, the PKD1 group), MCL4 (22, TBC1D3) |
| **C_tree** (clause 5, ordered) | §0★★ clause 5. Take all SH-aLRT > 75 splits of the §6js exon **or** intron IQ-TREE, restricted to leaves common to both trees; cluster = the smaller side; keep splits compatible with every other kept split (Buneman). C_tree_top = maximal clusters; C_tree_min = minimal clusters | `light/C.supported_clades.tsv`; `lo_analysis.c_tree` → `IS/ctree_clusters.tsv` | **30** leaves (21 NPIP + 9 TBC1D3; NPIPB14P is missing from the exon tree) | top: NPIPA (7), {B10P,B15,B6-9}, {B12,B13,B3-5}, TBC1D3 {B,F,G,H}, {TBC1D3,D,E,K}. min: {A6,A9}, {A7,A8}, {B3,B5}, {B8,B9}, {B,H}, {TBC1D3,D,E,K} |
| literature C (CIRCULAR reference) | `light/scripts/layer_c.py`: a literature group is kept when some split of either tree has SH-aLRT > 75. TBC1D3 cluster1/2 are positional labels and were cleared (§3 #9) | `light/work/C.out` | 31 leaves | C_L1: NPIPA, NPIPB. C_mid_ab (as built): named NPIPB {B3,B4,B5,B11,B12,B13}. C_mid = C_mid_ab ∨ C_fine. C_fine: A6-9, B3-5, B6-9, B12/13 |
| **EXPR** | Primary reads (`-F 2308`). **any** (primary, spec): a CIGAR block overlaps an exon by ≥ 1 bp, strand ignored. *unique*: the read hits exactly one RefSeq record genome-wide. *unique_mr*: same, but ignoring the 6 readthrough records that overlap a member on the same strand (PKD1P3-NPIPA1, LOC131696449, PKD1P4-NPIPA8, PKD1P5-LOC105376752, PDXDC2P-NPIPB14P, TBC1D3P1-DHX40P1) | `heavy/scripts/expr_counts.py`. Recount: `lo_expr_recount.py --ignore …`, 0 of 321 genes differ on any/unique; unique_mr rises for NPIPA1 0→9, NPIPA9 0→3, NPIPB14P 0→2, TBC1D3P1 0→1 | all 68 U genes | — |

Checks:
- **P stability.** Member clusters are identical from k = 2 to k = 8 (|N_k| = 59, 138, 466, 727, 1,230, 1,972, 2,760; at k = 8, 2 unsearched genes were dropped). The stopping rule is therefore not binding here (`IS/p_variants.out`).
- **P at aa identity ≥ 0.50** (the §6ko qualification's non-superfamily stratum):
  - The member neighbourhood closes at |N_k| = 32 for k = 2..5.
  - MCL member clusters: 21 NPIP and exactly the 9 TBC1D3.
  - Exact components: {21 NPIP} and {9 TBC1D3 + USP6 + USP32}. USP6 is joined to all 9 TBC1D3 copies at aa 0.806-0.816.
  - TBC1D26's plain-rule P edges have aa 0.426-0.487; USP6NL's have 0.363-0.426.
- **PC1** is an exact connected component of the §6ko graph built from all 4,430 searches.
- **PC2** is an MCL cut inside a component. The component size depends on how equal-bitscore HSPs are ordered:
  - Shipped `bench/protein_families.edges_from` (greedy HSP order): component ≥ 4,522 genes, 965 never searched (`IS/p_components.out`).
  - Stable file order, as in the P verifier's reimplementation: 4,528 genes, 968 never searched (`IS/p_tie_order_component.out`).
  - Under both orders PC2 has direct edges to 12 outside genes (EVI5, EVI5L, GRTP1, SGSM3, TBC1D10A/B/C, TBC1D12, TBC1D14, TBC1D2, TBC1D28, USP6). The verifier's figure of 10 was not reproduced.
- **D groups are MCL cuts too** (`IS/d_groups_cut.out`):
  - MCL4 has edges to 4 distinct outside genes (RNFT1, USP6, RNFT1P3, NPEPPS), 18 gene-level edges in total.
  - MCL0 has edges to 25 outside genes (267 gene-level edges); MCL7 to 32 (153).
  - Between MCL0 and MCL7 there are weight-1.000 edges through the PKD1P*-NPIP readthroughs.

## 3. Build corrections (pre-audit, kept; corrected tables live next to the originals as `*.corrected.tsv`)

| # | issue (verifier, severity) | fix applied | effect |
|---|---|---|---|
| 1 | The member rule was applied differently to the two families (P and D verifiers) | One rule for both families. The RefSeq description must match `nuclear pore complex[- ]interacting protein` or `^TBC1 domain family member 3( \|[A-Z]\|$)`. A readthrough that names NPIP/TBC1D3 counts only if it overlaps no same-strand description-rule record. This restriction applies **to readthroughs only**; it is not a general "one record per copy" rule (see §4) | Members 44 → **46** (NPIP 27, TBC1D3 19). **Removed TBC1D3P1-DHX40P1**. P is unchanged |
| 2 | Non-coding members were written as P singletons | Non-coding members: `P\|NA`. **Audit:** `P.groups.corrected.tsv` now also carries PKD1 and DHX40 as `P\|other`, and `D.groups.corrected.tsv` carries TBC1D26 → `D\|c15_17_22\|MCL24`. Both tables are asserted equal to the universe labels | P universe 34; D universe 63 |
| 3 | PC2 is an MCL cut; tie order of equal-bitscore HSPs | §2 checks (both tie orders) | — |
| 4 | P provenance incomplete | `light/scripts/README.md` written (reconstructed) | — |
| 5 | D is not built with clause-2 edges (important) | **Relabel route.** "E1" in every headline row, plus the failed-guard and GENCODE sentence (top of report). No clause-2/4 rebuild | §6 GENCODE variant |
| 6 | D catalogs cover chromosome trios only | Genes outside both catalogs are *not in D*. No 6-chromosome catalog was built | NPIPB1P, TBC1D3P6, LOC124905656, LOC100420289 and USP6NL are outside D (probe in §9) |
| 7 | `truth_guided` E1 is D itself | E1 is not a truth. D vs E0 is reported only as construction sensitivity | §7 |
| 8 | C is literature-circular | **Audit:** the ordered C layer is now the clause-5 split system (C_tree). Literature C is a circular reference column | §5 |
| 9 | TBC1D3 cluster1/cluster2 are positional labels | Cleared in *literature* C_L1. In the clause-5 split system, {TBC1D3,D,E,K} vs {B,F,G,H,I} is a supported compatible split (intron 83/71), so it stays in C_tree | Literature C: 0 TBC1D3 pairs (NA). C_tree: 12 top pairs |
| 10 | The universe still carried the deferred S1 layer | U = members ∪ P ∪ D ∪ C | 79 → **68 genes** |
| 11-19 | `truth_guided` clusters; duplicate node keys; fold flags; weak truths; EXPR rows; seed_family; RNFT1 reads; unique definition | As before. **Audit:** any-overlap is now the primary count, and unique_mr has been added | §8 |
| 20 | As-built `C_mid` is not a hierarchy level | `C_mid` := C_mid_ab ∨ C_fine. **Audit:** C_mid_ab stays in every main table | §5 |

## 4. Universe U (68 genes) and members

- **Composition.**
  - 46 members.
  - Layers placing a gene with a member: P 32 genes, D 62, C_tree 30 leaves, literature C 31 leaves.
  - Layer universes: P 34, D 63, C_tree 30, C_lit 31.
  - Chromosomes: chr16 40, chr17 23, chr1 2, chr4 1, chr10 1, chr18 1.
  - Table: `light/universe.corrected.tsv`.
- **Member counts depend on a record convention.** LOC124905656 (chr1, 8.9 kb, "TBC1D3K-like") encloses member TBC1D3P6 on the same strand, and both records are kept. If the enclosed record is folded into its encloser, there are **45 members (TBC1D3 18)** (`IS/corrected_tables.out`). Neither record is in P, D or C, so no containment changes.
- **Fold inheritance (D).** A record that `mcl_families` folds into another locus is given that locus' cluster. This is how TBC1D3P1-DHX40P1, RNFT1-DT and DHX40P1 sit in MCL4 (folded into TBC1D3P1) and TBC1D26 sits in MCL24 (folded into ZNF286A). §6ks reported TBC1D3P1-DHX40P1 as "not in the main family" on the same catalog, because it scored the locus, not the folded record. 11 U genes are fold-inherited; §6 removes them as a variant.

**Member reconciliation with §6jg** (22 NPIP records, 9 TBC1D3 protein-coding copies; `IS/member_reconciliation.tsv`):

| family | members here | in §6jg truth | extra records and reason |
|---|---|---|---|
| NPIP | 27 | 22 | LOC124907807, LOC124907808: coding LOC records described as "NPIP family member B15". LOC124907834, LOC128966608: transcribed pseudogenes, "B13-like". PKD1P6-NPIPP1: a readthrough whose NPIP part has no gene record of its own |
| TBC1D3 | 19 | 9 | TBC1D3P1-P7 (7 pseudogenes) and LOC100420311 (chr17), LOC124905656 (chr1), LOC100420289 (chr4): §6jg's TBC1D3 truth is only the 9 protein-coding copies. TBC1D3P6 lies inside LOC124905656 |

## 5. Containment

- **Cell format.** c(row ⊇ col) [col pairs on genes in both universes] followed by g/G. G is the number of col groups with ≥ 2 genes; g is how many of them lie inside one row group. "vac" means 0 pairs, so the value is undefined.
- **Ordered layers:** P, D, C_tree_top, C_tree_min.
- **Reference columns** (literature, CIRCULAR): C_L1, C_mid, C_mid_ab, C_fine.

**Pooled** (genes in both: P/D 33, P/C_tree 28, D/C_tree 29, P/C_lit 28, D/C_lit 30):

| X ⊇ Y | P | D | C_tree_top | C_tree_min | C_L1 | C_mid | C_mid_ab | C_fine |
|---|---|---|---|---|---|---|---|---|
| **P** | — | 0.965 [255] 1/2 | 1.000 [53] 5/5 | 1.000 [11] 6/6 | 1.000 [87] 2/2 | 1.000 [27] 3/3 | 1.000 [15] 1/1 | 1.000 [16] 4/4 |
| **D** | 0.965 [255] 1/2 | — | 1.000 [58] 5/5 | 1.000 [11] 6/6 | 1.000 [112] 2/2 | 1.000 [27] 3/3 | 1.000 [15] 1/1 | 1.000 [16] 4/4 |
| **C_tree_top** | 0.256 [207] 0/2 | 0.257 [226] 0/2 | — | 1.000 [11] 6/6 | 0.411 [112] 1/2 | 0.815 [27] 2/3 | 0.667 [15] 0/1 | 1.000 [16] 4/4 |
| **C_tree_min** | 0.053 [207] 0/2 | 0.049 [226] 0/2 | 0.190 [58] 1/5 | — | 0.036 [112] 0/2 | 0.148 [27] 0/3 | 0.067 [15] 0/1 | 0.250 [16] 0/4 |
| **C_L1** | 0.420 [207] 0/2 | 0.455 [246] 0/2 | 0.793 [58] 3/5 | 0.364 [11] 4/6 | — | 1.000 [27] 3/3 | 1.000 [15] 1/1 | 1.000 [16] 4/4 |
| **C_mid** | 0.130 [207] 0/2 | 0.110 [246] 0/2 | 0.379 [58] 1/5 | 0.364 [11] 4/6 | 0.214 [126] 0/2 | — | 1.000 [15] 1/1 | 1.000 [16] 4/4 |
| **C_mid_ab** | 0.072 [207] 0/2 | 0.061 [246] 0/2 | 0.172 [58] 1/5 | 0.091 [11] 1/6 | 0.119 [126] 0/2 | 0.556 [27] 1/3 | — | 0.250 [16] 2/4 |
| **C_fine** | 0.077 [207] 0/2 | 0.065 [246] 0/2 | 0.276 [58] 0/5 | 0.364 [11] 4/6 | 0.127 [126] 0/2 | 0.593 [27] 2/3 | 0.267 [15] 0/1 | — |

**NPIP** (P/D genes 22):

| X ⊇ Y | P | D | C_tree_top | C_tree_min | C_L1 | C_mid | C_mid_ab | C_fine |
|---|---|---|---|---|---|---|---|---|
| **P** | — | 1.000 [210] 1/1 | 1.000 [41] 3/3 | 1.000 [4] 4/4 | 1.000 [87] 2/2 | 1.000 [27] 3/3 | 1.000 [15] 1/1 | 1.000 [16] 4/4 |
| **D** | 1.000 [210] 1/1 | — | 1.000 [46] 3/3 | 1.000 [4] 4/4 | 1.000 [112] 2/2 | 1.000 [27] 3/3 | 1.000 [15] 1/1 | 1.000 [16] 4/4 |
| **C_tree_top** | 0.240 [171] 0/1 | 0.242 [190] 0/1 | — | 1.000 [4] 4/4 | 0.411 [112] 1/2 | 0.815 [27] 2/3 | 0.667 [15] 0/1 | 1.000 [16] 4/4 |
| **C_tree_min** | 0.023 [171] 0/1 | 0.021 [190] 0/1 | 0.087 [46] 0/3 | — | 0.036 [112] 0/2 | 0.148 [27] 0/3 | 0.067 [15] 0/1 | 0.250 [16] 0/4 |
| **C_L1** | 0.509 [171] 0/1 | 0.533 [210] 0/1 | 1.000 [46] 3/3 | 1.000 [4] 4/4 | — | 1.000 [27] 3/3 | 1.000 [15] 1/1 | 1.000 [16] 4/4 |

**TBC1D3** (P/D genes 11; C_tree 9):

| X ⊇ Y | P | D | C_tree_top | C_tree_min | C_L1 / C_mid / C_mid_ab / C_fine |
|---|---|---|---|---|---|
| **P** | — | 0.800 [45] 0/1 | 1.000 [12] 2/2 | 1.000 [7] 2/2 | vac [0] |
| **D** | 0.800 [45] 0/1 | — | 1.000 [12] 2/2 | 1.000 [7] 2/2 | vac [0] |
| **C_tree_top** | 0.333 [36] 0/1 | 0.333 [36] 0/1 | — | 1.000 [7] 2/2 | vac [0] |
| **C_tree_min** | 0.194 [36] 0/1 | 0.194 [36] 0/1 | 0.583 [12] 1/2 | — | vac [0] |
| **C_L1 / C_mid / C_mid_ab / C_fine** | 0.000 [36] 0/1 | 0.000 [36] 0/1 | 0.000 [12] 0/2 | 0.000 [7] 0/2 | vac [0] |

**One common gene set** (`IS/tournament_common_genes.tsv`):
- The set is the 28 genes in every universe (P ∩ D ∩ C_tree ∩ C_lit): NPIP 19, TBC1D3 9.
- **P ⊇ D and D ⊇ P are both 1.000 [207], 2/2**, so on these genes P and D are the same partition.
- Every pooled verdict is identical to the verdict on the full universes.
- Across all 8 layers there is no intransitive triple. Most of these relations hold by construction (§6).

The clause-5 hierarchy as a whole (every compatible supported cluster, not only the two partition levels) sits inside P and D:
- NPIP: 12 clusters. 11 of 11 lie inside one P group (the clusters whose leaves are all in P's universe) and 12 of 12 inside one D group.
- TBC1D3: 4 of 4 inside one P group and 4 of 4 inside one D group (`IS/analysis.out`).

## 6. Order: what is measured, what is by construction

X is above Y iff c(X ⊇ Y) > c(Y ⊇ X). A comparison is NA when either side has 0 pairs.

| relation | status | evidence |
|---|---|---|
| P vs D | **undetermined**. Tie under the spec reading; the verdict changes with conventions (table below) | NPIP: one P group vs one D group (1/1). TBC1D3: 36/45 both ways (0/1). Pooled: 246/255 both ways (1/2) |
| P, D above C_tree / literature C | **by construction** that c(C ⊇ P) < 1 (subfamilies are finer than the family). **Measured:** c(P ⊇ C) = c(D ⊇ C) = 1, i.e. 0 clade pairs split | §5; `IS/disagree_clades.tsv` |
| C_L1 ⊇ C_mid ⊇ C_fine | **by construction** (literature hierarchy; C_mid := C_mid_ab ∨ C_fine) | as built, C_fine above C_mid_ab: 0.267 [15] vs 0.250 [16] |
| C_tree_top ⊇ C_tree_min | **by construction** (maximal vs minimal clusters of one laminar system) | 1.000 [11] 6/6 |
| TBC1D3: P, D vs literature C | **NA** (0 literature pairs) | — |
| TBC1D3: P, D vs C_tree | P, D above C_tree_top: 1.000 [12] vs 0.333 [36] | not vacuous |

The pre-audit report called the order "transitive". That is dropped: its pairwise values were computed on different gene intersections (P/D 33, P/C 28, D/C 30). On the one common gene set of 28 genes (§5) the pooled verdicts have no intransitive triple, but there P and D are the same partition, and the C relations hold by construction.

**How robust is the P~D tie?** (`lo_pd_variants.py` → `IS/pd_variants.tsv`, `IS/pd_variants.out`). Each variant re-closes U. Cells show pairs followed by group-level nesting.

| variant | U | TBC1D3 c(P⊇D) · c(D⊇P) | pooled c(P⊇D) · c(D⊇P) | verdict (pooled) |
|---|---|---|---|---|
| spec reading (P universe = every U gene with a protein; D as built) | 68 | 36/45 0/1 · 36/45 0/1 | 246/255 1/2 · 246/255 1/2 | tie |
| P verifier's universe rule (only the 32 genes P clusters with a member) | 68 | 36/36 1/1 · 36/45 0/1 | 246/246 2/2 · 246/255 1/2 | **P > D** |
| pre-audit group tables on disk (P 32 rows; no TBC1D26 row in D) | 68 | 36/36 · 36/36 | 246/246 · 246/246 | tie at 1.000 |
| P at aa ≥ 0.50, k = 2 MCL member clusters (TBC1D26 and USP6NL leave PC2) | 66 | 36/45 0/1 · 36/36 1/1 | 246/255 1/2 · 246/246 2/2 | **D > P** |
| P at aa ≥ 0.50, exact components (adds USP6 and USP32 to the TBC1D3 group) | 68 | 36/45 0/1 · 36/55 0/1 | 246/255 1/2 · 246/265 1/2 | **P > D** |
| D without fold-inherited membership (11 folded U records leave D) | 52 | 36/45 0/1 · 36/36 1/1 | 246/255 1/2 · 246/246 2/2 | **D > P** |
| P aa ≥ 0.50 MCL + D without folds | 50 | 36/45 · 36/36 | 246/255 · 246/246 | **D > P** |
| D = GENCODE E1 (same construction; RefSeq gene → GENCODE node by largest shared exonic bp) | 70 | 36/55 0/1 · 36/45 0/1 | 226/286 0/2 · 226/235 1/2 | **D > P** (NPIP too: 190/231 vs 190/190) |
| D = GENCODE E1, mapping kept only at ≥ 0.5 shared of the shorter exon union | 66 | 36/55 · 36/45 | 207/265 · 207/216 | **D > P** |
| P aa ≥ 0.50 MCL + D = GENCODE E1 | 68 | 36/55 · 36/36 | 226/286 · 226/226 | **D > P** |

Leave-one-out and leave-two-out (spec reading, `IS/pd_leave_out.tsv`):
- Only removing DHX40 (→ P > D) or TBC1D26 (→ D > P) changes a verdict; 2 of 68 genes.
- Over all 2,278 gene pairs removed, the pooled verdict is P > D 66, D > P 66, tie 2,146.

What drives each disagreement:
- **GENCODE E1** puts PKD1 (and PKD1P6/PKD1P3 records) and CLN3 in the NPIP group (`IS/gencode_map.tsv`). This agrees with §6kl, where 41% of the GENCODE-vs-RefSeq false pairs on this hold-out come from one PKD1 + PKD1P + NPIPA/B family.
- **In GENCODE**, USP6 joins the TBC1D3 group and TBC1D26 falls in a different group (MCL18).
- **DHX40 in RefSeq E1** is the §6ko AN-2 example of D-only co-membership through readthrough paths ("TBC1D3–DHX40").
- **The plain P rule** gives superfamilies (§6ko qualification: 49.5% of pairs in disjoint HGNC groups). At aa ≥ 0.50 that share is 0.4%, and TBC1D26 (aa 0.454-0.487 to members) leaves PC2.

## 7. Enforcement cost and truth agreement

**Conventions** (spec AUDIT AMENDMENT):
- JOIN is run two ways: on D groups cut to U, and on **whole catalog groups** with U re-closed.
- A gene outside the coarser layer's universe is a free choice under REFINE. Every variant is reported: *attach*, *single*, and for C *together*.

| enforced layer | pairs on the layer's own universe | change | notes |
|---|---|---|---|
| P JOIN over D, U-restricted | 265 → 276 | +11: DHX40 with PC2 | +30 genes (non-coding D genes) join the universe |
| P JOIN over D, whole groups | 265 → 276 | +11: DHX40 with PC2 | +37 genes. Whole MCL24 is pulled in: FOXO3B, TBC1D26-AS1, TBC1D27P, TBC1D28, ZNF286A, ZNF286A-TBC1D26, ZNF286B. No further chaining through the k = 2 P clusters (TBC1D28 is a P singleton; FOXO3B and ZNF286A lie outside N_2) |
| P REFINE in D, USP6NL attached | 265 → 255 | −10: TBC1D26 leaves PC2 | USP6NL stays (largest summed P weight) |
| P REFINE in D, USP6NL single | 265 → 246 | −19: TBC1D26 and USP6NL leave | — |
| C_L1 REFINE in D: attach / single / together | 126 → 126 / 112 / 112 | 0 / −14 / −14 | NPIPB1P is outside D |
| C_L1 REFINE in P: attach / single / together | 126 → 126 / 87 / 90 | 0 / −39 / −36 | NPIPB1P, B10P and B14P are outside P |
| C_tree_top REFINE in P: attach / single / together | 58 → 58 / 53 / 53 | 0 / −5 / −5 | NPIPB10P is outside P |
| C_mid, C_fine, C_tree_min REFINE in D or in P (every variant); C_tree_top REFINE in D | unchanged | 0 | — |

Every C change comes from the outside-universe convention. **No in-universe clade pair is split** by REFINE in D or in P. After enforcement:
- c(P_join ⊇ D) = 1.000 [711] 3/3, same for the whole-group JOIN.
- c(D ⊇ P_ref) = 1.000 [246] 2/2 under both conventions (`IS/analysis.out`).

**(A) Truth agreement on U.** Recall here is **conditioned on the prediction**: U is the closure of the scored layers, so a gene no layer found cannot be missed. Cells show pairwise precision / recall / bipartite F (§6ks count matching; Jaccard matching gives the same F on all 59 rows where both are defined). HGNC 2227 is **superfamily-level**. F is NA when either side has 0 pairs (`IS/truth_agreement.tsv`).

| layer / variant | gene set | HGNC, TBC1D3 | Soto, TBC1D3 | Soto, NPIP | Soto, pooled |
|---|---|---|---|---|---|
| P as built | P universe | 1.000 / 1.000 / 1.000 (n 12) | 0.800 / 1.000 / 0.909 (11) | 0.693 / 0.981 / 0.842 (19) | 0.717 / 0.986 / 0.867 (30) |
| P JOIN, U-restricted | P universe | 0.833 / 1.000 / 0.917 (12) | 0.655 / 1.000 / 0.818 (11) | = as built | 0.683 / 0.986 / 0.833 (30) |
| P JOIN, U-restricted | own universe | 0.725 / 1.000 / 0.857 (14) | 0.475 / 1.000 / 0.688 (16) | 0.588 / 0.897 / 0.806 (31) | 0.553 / 0.922 / 0.766 (47) |
| P JOIN, whole groups | P universe, re-closed | 0.725 / 1.000 / 0.857 (14: + TBC1D28, ZNF286A) | 0.418 / 1.000 / 0.643 (14) | = as built | 0.590 / 0.986 / 0.758 (33) |
| P JOIN, whole groups | own universe, re-closed | 0.532 / 1.000 / 0.737 (19) | 0.321 / 1.000 / 0.550 (20) | 0.588 / 0.897 / 0.806 (31) | 0.477 / 0.924 / 0.706 (51) |
| P REFINE in D, USP6NL attached | P universe | 1.000 / 0.818 / 0.957 (12) | 1.000 / 1.000 / 1.000 (11) | = as built | 0.751 / 0.986 / 0.900 (30) |
| P REFINE in D, USP6NL single | P universe | 1.000 / 0.655 / 0.909 (12) | 1.000 / 1.000 / 1.000 (11) | = as built | 0.751 / 0.986 / 0.900 (30) |
| D (E1; unchanged by either operator) | D universe | 0.682 / 0.818 / 0.800 (13) | 0.543 / 1.000 / 0.750 (16) | 0.588 / 0.897 / 0.806 (31) | 0.575 / 0.922 / 0.787 (47) |

JOIN vs REFINE, per variant:
- **HGNC (TBC1D3), F:** REFINE-attach 0.957 > JOIN-U 0.917 > REFINE-single 0.909 > JOIN-whole 0.857.
- **Soto (TBC1D3), F:** REFINE 1.000 (both conventions) > JOIN 0.818 / 0.643.
- **Neither statement holds under all variants on HGNC.** Each operator changes 1 gene per side: JOIN adds DHX40; REFINE removes TBC1D26, and USP6NL depending on the convention.
- DHX40 enters only through the readthrough fold. TBC1D26 and USP6NL are TBC-domain paralogs.

**(B) Member-anchored, layer-independent truth** (`IS/truth_member_anchored.tsv`). Pairs with ≥ 1 member endpoint, scored over every gene of each member-holding truth group, restricted to the layer's universe genome-wide (P: genes with a protein; D: catalog nodes).

| truth | layer / variant | genes (members) | truth pairs | P / R |
|---|---|---|---|---|
| HGNC 2227 (57 symbols, 54 RefSeq genes; 52 with a protein; 27 catalog nodes) | P as built | 54 (9) | 423 | 1.000 / **0.128** |
| | P JOIN U / JOIN whole | 54 (9) | 423 | 0.857 / 0.128 · 0.875 / 0.149 |
| | P REFINE attach / single | 54 (9) | 423 | 1.000 / 0.106 · 1.000 / 0.085 |
| | D (E1) | 30 (9) | 198 | 0.714 / **0.227** |
| Soto families holding a member (ID_149, 152-155, 468, 469; 66 Soto genes; 43 RefSeq genes map to them with flag ok, **all already in U**) | P as built | 30 (27) | 144 | 0.717 / 0.986 |
| | P JOIN U = JOIN whole | 30 (27) | 144 | 0.686 / 0.986 |
| | P REFINE attach = single | 30 (27) | 144 | 0.751 / 0.986 |
| | D (E1) | 47 (34) | 203 | 0.541 / 0.911 |

Reading (B):
- On HGNC, P's precision of 1.000 comes with a recall of 0.128: P recovers the TBC1D3 family, not the "TBC1 domain" superfamily group.
- On Soto, closing U over the layers inflates nothing: no RefSeq gene outside U maps (flag ok) into a member-holding Soto family.

**C against literature groups.**

| layer | truth | NPIP | TBC1D3 | pooled |
|---|---|---|---|---|
| C_L1 (CIRCULAR) | literature L1 | 1.000 / 1.000 / 1.000 (22) | — | — |
| C_mid_ab (CIRCULAR) | named NPIPB subfamily | 1.000 / 1.000 / 1.000 (22) | — | — |
| C_mid (= C_mid_ab ∨ C_fine) | named ∨ L2 | **not scored: identical by construction** | | |
| C_fine (CIRCULAR) | literature L2 | 1.000 / 1.000 / 1.000 (22) | P NA / R 0.000 / F NA (0 predicted pairs) | 1.000 / 0.889 / 0.967 (31) |
| C_tree_top (clause 5) | literature L1 | 1.000 / 0.411 / 0.765 (21) | — | — |
| C_tree_top (clause 5) | literature L2 | 0.348 / 1.000 / 0.667 (21) | 0.167 / 1.000 / 0.444 (9) | 0.310 / 1.000 / 0.600 (30) |
| C_tree_min (clause 5) | literature L2 | 1.000 / 0.250 / 0.833 (21) | 0.286 / 1.000 / 0.667 (9) | 0.545 / 0.333 / 0.778 (30) |
| C_tree rooted at its most balanced split (variant) | literature L1 | 1.000 / 1.000 / 1.000 (21) | — | — |

Literature groups vs the compatible split system (`IS/ctree_literature_groups.tsv`):
- NPIPA, NPIPB, A6-9, B3-5 and B6-9 are in it.
- **Two literature groups conflict with a supported split:**
  - B12/13 (intron 98.9/100) conflicts with the exon split {B12,B3,B4,B5} (90.3/89).
  - **The named NPIPB subfamily** {B3,B4,B5,B11,B12,B13} (intron 99.6/100) conflicts with {B10P,B11,B15,B6-9} (exon 86.6/75).
- So C_fine (B12/13) and C_mid_ab (named) each contain one group that fails the compatibility rule.
- TBC1D3 AE and CDKL have no supported split in either tree.

## 8. EXPR operator (testis; primary rule any-overlap, reads ≥ 3)

**Edges used:**
- P: `P.edges`, complete inside PC1 and PC2: 265 = C(21,2) + C(11,2).
- D: `D.edges.corrected`.
- C layers: co-membership.
- P_join: P ∪ D.

Member-holding L groups with ≥ 2 genes (`IS/expr_groups.tsv`, layer edges):

| layer group (genes) | any ≥ 3 → EXPR groups | unique ≥ 3 | unique_mr ≥ 3 |
|---|---|---|---|
| P PC1 (21) = D MCL0 (30) = P_ref PC1 | {A1, A2, A7, A9, B2, B4, B5, B6, B7, B9, B11} | {A2, A7, B2, B4, B5, B6, B9, B11} | + A1, A9 (10 genes) |
| P PC2 (11) | – (TBC1D26 expressed, no expressed partner: dropped) | – | – |
| D MCL7, PKD1 group (10) | {PKD1, PKD1P3-NPIPA1, PKD1P5-LOC105376752, PKD1P6, PKD1P6-NPIPP1} | – | – |
| D MCL4, TBC1D3 (22) | **{DHX40, DHX40P1, RNFT1-DT, TBC1D3P1-DHX40P1}** | **{DHX40, RNFT1-DT}** | {DHX40, RNFT1-DT} |
| C_L1 NPIPA (7) = C_tree_top NPIPA | {A1, A2, A7, A9} | {A2, A7} | {A1, A2, A7, A9} |
| C_L1 NPIPB (15) | {B1P, B2, B4, B5, B6, B7, B9, B11} | without B7 | without B7 |
| C_tree_top {B10P,B15,B6-9} / {B12,B13,B3-5} | {B6, B7, B9} / {B4, B5} | {B6, B9} / {B4, B5} | {B6, B9} / {B4, B5} |
| C_fine A6-9 / B3-5 / B6-9 / B12/13 | {A7, A9} / {B4, B5} / {B6, B7, B9} / – | – (A7 dropped) / {B4, B5} / {B6, B9} / – | {A7, A9} / {B4, B5} / {B6, B9} / – |
| C_tree_min (every cluster), C_tree_top TBC1D3 clusters | – | – | – |

**No EXPR(D) group in MCL4 is TBC1D3 expression.** DHX40 has 8/8 reads (any/unique) and RNFT1-DT 52/47.

**Threshold sweep** (`IS/expr_sweep.tsv`, layer edges). Cells show members expressed, NPIP of 27 / TBC1D3 of 19.

| mode | t ≥ 1 | t ≥ 2 | t ≥ 3 | t ≥ 4 | t ≥ 5 | L groups split |
|---|---|---|---|---|---|---|
| any (primary) | 24 / 4 | 21 / 0 | **13 / 0** | 9 / 0 | 6 / 0 | none at any t |
| unique | 20 / 2 | 14 / 0 | 9 / 0 | 5 / 0 | 5 / 0 | t ≥ 1: D MCL4 and P_join → {TBC1D3, TBC1D3E} \| {DHX40, RNFT1-DT} |
| unique_mr | 23 / 3 | 17 / 0 | 11 / 0 | 6 / 0 | 6 / 0 | t ≥ 1: D MCL4 and P_join → {DHX40, RNFT1-DT} \| {TBC1D3, TBC1D3E, TBC1D3P1} |

**Distinct EXPR groups** at any ≥ 3: 9 with layer edges, 10 with co-membership. The extra group is P_join with TBC1D26. Across all layers × modes × edge settings at t = 3 the count is 102, but those are the same few groups repeated.

Checks:
- **EXPR(L) ⊆ L:** 0 violations. This holds **by construction**, since components are computed inside L groups.
- **Co-membership vs layer edges:** they differ on 1 of 90 member-holding group × mode rows (P_join, any: TBC1D26 has no P or D edge to the expressed genes).
- **T2 precondition E_M ⊆ E_L inside each L group:** holds for **all 30** pairs (L, M) where M refines L on these data (`IS/analysis.out`, `[T2 precondition]`). By the corrected T2 (spec AUDIT AMENDMENT, definition doc §0★★ clause 6 note), EXPR(M) therefore refines EXPR(L) automatically.
- **T2 checks:** at t = 3 with layer edges there were 63 (any), 50 (unique) and 62 (unique_mr), all with 0 violations. They are **guaranteed**, not evidence. With co-membership edges T2 is guaranteed regardless. Counterexample for when the precondition fails: spec, EXPR.
- **EXPR(L) truth agreement not reported.** At t = 3 no L group splits in any mode, so EXPR(L) is L restricted to expressed genes minus dropped singletons, and its truth scores add nothing to §7.

Unexpressed members (`IS/member_expression.tsv`):
- **TBC1D3, all 19** in every mode. The protein-coding copies have 0-1 reads.
- **NPIP, any < 3 (14):** LOC124907807, LOC124907808, LOC124907834, LOC128966608, A5, A6, A8, B3, B8, B10P, B12, B13, B14P, B15.
- **NPIP, unique < 3 (18):** the 14 plus A1, A9, B7, PKD1P6-NPIPP1. **unique_mr < 3 (16):** the 14 plus B7 and PKD1P6-NPIPP1.
- **NPIPB1P** (6/6) is expressed but outside both P and D.

Why any and unique differ (co-hit reads; `lo_expr_cohits.py` → `IS/expr_cohits.out`):

| gene | any / unique / unique_mr | explanation |
|---|---|---|
| TBC1D3P1-DHX40P1 | 72 / 0 / 0 | 66 of 72 reads also hit RNFT1. Alignment strand same/opposite is 65/7 over **all 72 reads** |
| DHX40P1 | 5 / 0 / 0 | all 5 also hit RNFT1-DT and the readthrough; all 5 on the opposite strand |
| NPIPA1 | 10 / 0 / 9 | all 10 also hit the PKD1P3-NPIPA1 readthrough; 1 also hits NPIPA2 |
| NPIPA9 | 3 / 0 / 3 | all 3 also hit the PKD1P5-LOC105376752 readthrough |
| NPIPB14P | 2 / 0 / 2 | both also hit PDXDC2P-NPIPB14P |
| NPIPB7 | 3 / 0 / 0 | all 3 also hit CLN3 (not a readthrough) |
| NPIPB9 | 16 / 9 / 9 | 6 also hit EIF3C |
| TBC1D26 | 7 / 0 / 0 | all 7 also hit LOC105371559; 5 also hit ZNF286A-TBC1D26 |

## 9. Gene-level disagreements

**P groups these genes with members; D does not** (`IS/disagree_P_not_D.tsv`):

| gene | P group | D status | blastp to members: identity / coverage of the longer protein (9 of 9 edges) |
|---|---|---|---|
| TBC1D26 (chr17) | PC2 | D MCL24. **D has no TBC1D3 edge for TBC1D26:** its record's only catalog edges are ZNF286A-TBC1D26 1.000, TBC1D26-AS1 1.000, TBC1D28 0.985 and TBC1D27P 0.813. It is labelled MCL24 because the record is folded into the ZNF286A locus, together with ZNF286A-TBC1D26 and TBC1D26-AS1 (`IS/d_groups_cut.out`) | 0.454-0.487 / 0.600-0.689 (best TBC1D3F 0.487 / 0.689) |
| USP6NL (chr10) | PC2 | outside D (no catalog) | 0.363-0.380 / 0.364-0.412 (best TBC1D3I 0.372 / 0.412) |

**D groups these genes with members; P does not** (31 genes, `IS/disagree_D_not_P.tsv`):
- Identity and coverage are those of the **best edge to a member in the same D group**, from the catalog PAF under the D edge rule. Identity × coverage reproduces the D weight on 28 of 28 edges.
- Flags: f = folded into another locus; rt = readthrough record.

| group | genes (identity / coverage of the longer gene's exon union) |
|---|---|
| MCL0, non-coding NPIP members (outside P) | NPIPB10P 0.990/1.000 · NPIPB14P 0.963/0.991 · LOC124907834 0.981/1.000 · LOC128966608 0.982/1.000 |
| MCL0, non-member non-coding | LOC100505915 0.982/1.000 · LOC112268174 (f) 0.988/1.000 · LOC124907811 (f) 0.986/1.000 · LOC124907844 (f) 0.995/1.000 · PDXDC2P-NPIPB14P (f, rt) 1.000/1.000 |
| MCL7, the PKD1 group | Its only member is PKD1P6-NPIPP1 (rt, folded into PKD1P6). PKD1 (coding; P other) 0.889/0.960 · PKD1P1 (f) 0.914/1.000 · PKD1P2 0.896/1.000 · PKD1P3 (f) 0.988/0.580 · PKD1P6 1.000/0.467 · PKD1P3-NPIPA1 (rt) 0.939/1.000 · PKD1P4-NPIPA8 (rt) 0.935/1.000 · PKD1P5-LOC105376752 (rt) 0.931/1.000 · LOC131696449 (rt) 0.948/1.000. **The four readthroughs also have weight-1.000 D edges to NPIP members in MCL0** (NPIPA1, NPIPA8, NPIPA9, NPIPA6) |
| MCL4, non-coding TBC1D3 members (outside P) | TBC1D3P1 0.987/0.984 · P2 0.976/1.000 · P3 0.920/1.000 · P4 0.920/1.000 · P5 0.846/1.000 · P7 0.810/0.693 · LOC100420311 0.910/0.994 |
| MCL4, non-members | TBC1D29P 0.944/1.000 · LOC100420408 0.869/0.534 · TBC1D3P1-DHX40P1 (f, rt) 1.000/1.000 |
| MCL4 genes with **no member edge** | RNFT1-DT (f; 1.000 to TBC1D3P1-DHX40P1) · DHX40P1 (f; 1.000 to TBC1D3P1-DHX40P1) · **DHX40** (coding; P other; 0.991 to RNFT1-DT) |

**Members that P and D place differently** (`IS/disagree_members_P_vs_D.tsv`):
- The 9 TBC1D3 coding copies differ by one non-member each: TBC1D26 in P, DHX40 in D.
- PKD1P6-NPIPP1 lies outside the main NPIP D group (MCL7).
- **0** of the 21 NPIP coding members are placed differently.

**Clade conflicts** (`IS/disagree_clades.tsv`):
- 0 clade pairs are split by P or D at any level: C_L1 87 (P) / 112 (D); C_mid 27; C_fine 16; C_tree_top 53 / 58; C_tree_min 11 / 11.
- Clade genes outside P: NPIPB10P, NPIPB14P, NPIPB1P (C_tree: NPIPB10P, NPIPB1P). Outside D: NPIPB1P.

**Off-catalog members** (probe rerun identical; `IS/probe/probe.out`, `short_probe_opts.out`; minimap2 of gene spans against the 52 genes of MCL0 and MCL4):

| gene | asm20 defaults | sensitive seeds (-k13 -w5) |
|---|---|---|
| NPIPB1P (chr18) | 27 records to MCL0; best 0.971 identity over 0.991 of its body | — |
| LOC124905656 (chr1) | 10 records to MCL4; best TBC1D3P7 0.795 identity over 0.491 of its body | — |
| TBC1D3P6 (chr1) | nothing | TBC1D3P5 0.849 over 465 bp |
| LOC100420289 (chr4) | nothing | TBC1D3P5 0.809 over 388 bp |
| USP6NL (chr10) | 5 records, each ≤ 0.7% of its 151 kb body (repeat-like) | — |

## 10. Known limits

- **Scope:** two development families; descriptive; not pre-registered; S1/S2/S3 deferred.
- **D is E1**, not §0★★ clause 2/4:
  - It failed its HGNC guard on chr16/19/20 (§6kl).
  - GENCODE E1 gives D > P.
  - It is per chromosome trio: chr1/4/10/18 genes lie outside it, and no chr16↔chr17 edge can exist.
  - Folded records inherit a cluster (§4).
- **P is coding-only and bounded:**
  - 16 of 46 members have no CDS.
  - Stable for k = 2..8, but PC2 is an MCL cut inside a TBC-domain component of ≥ 4,522 genes.
  - The plain rule builds superfamilies (§6ko qualification).
  - E-values come from a genome-wide database.
- **The P~D verdict** depends on conventions (§6). No convention is privileged by the spec except the P-universe rule.
- **C_tree:**
  - The estimator is clause 5 as in the definition doc.
  - Smaller side = cluster is a rooting convention; the rooted variant is in §7.
  - One reference alignment and one haplotype.
- **Truths:**
  - HGNC gives no group to any NPIP gene and is superfamily-level for TBC1D3.
  - Soto needs an exon-overlap mapping (20 of 68 U genes are not flag ok).
  - Literature C is circular.
  - TBC1D3 L2 groups are an assumed name mapping (§6jg).
- **EXPR:** one sparse tissue; 0 expressed TBC1D3 members; no split at t ≥ 2; T2 guaranteed on these data.

## 11. Reproduce (foreground, in order)

```
python3 bench/layer_order/lo_expr_recount.py IS/expr_recount.tsv --ignore PKD1P3-NPIPA1,LOC131696449,PKD1P4-NPIPA8,PKD1P5-LOC105376752,PDXDC2P-NPIPB14P,TBC1D3P1-DHX40P1 USP6NL LOC100420408 LOC100420311 TBC1D29P LOC124905656 LOC100420289   # 25 s
python3 bench/layer_order/lo_corrected_tables.py > IS/corrected_tables.out   # light/*.corrected.tsv (6 s)
python3 bench/layer_order/lo_p_variants.py > IS/p_variants.out               # P stability k=1..8, aa>=0.50, N_2 clusters (10 s)
python3 bench/layer_order/lo_analysis.py > IS/analysis.stdout                 # containment, tournaments, enforcement, truths A/B, EXPR, disagreements (4 s)
python3 bench/layer_order/lo_pd_variants.py > IS/pd_variants.out              # P~D variants, GENCODE E1 mapping, leave-out (5 s)
python3 bench/layer_order/lo_crosscheck.py > IS/crosscheck.out                # independent P/D/C_L1 re-derivation from the original tables
python3 bench/layer_order/lo_p_components.py > IS/p_components.out            # PC1/PC2 component vs MCL cut
python3 bench/layer_order/lo_d_groups_cut.py > IS/d_groups_cut.out            # D group cuts; TBC1D26 raw catalog edges
python3 bench/layer_order/lo_expr_cohits.py IS/expr_cohits.tsv TBC1D3P1-DHX40P1 DHX40P1 DHX40 RNFT1-DT NPIPA1 NPIPA9 NPIPB7 TBC1D26 PKD1P6-NPIPP1 NPIPB9 NPIPB14P TBC1D3P1 TBC1D3 TBC1D3E > IS/expr_cohits.out   # 35 s
bash    bench/layer_order/lo_offcatalog_probe.sh                             # IS/probe/* (rerun identical)
```

```
python3 bench/layer_order/lo_p_tie_order.py > IS/p_tie_order_component.out    # PC2 component under shipped vs stable HSP tie order (8 s)
```

## 12. Audit notes (each finding → resolution)

| # | finding (severity) | resolution |
|---|---|---|
| 1 | REFINE's outside-universe convention was unstated; "REFINE costs less than JOIN on both references" (important) | **Fixed.** Every variant is reported (§7): P attach/single, C attach/single/together. The directional sentence is replaced by per-variant F: HGNC REFINE 0.957 / 0.909 vs JOIN 0.917 / 0.857; Soto REFINE 1.000 / 1.000 vs JOIN 0.818 / 0.643. Stated: "no ranking holds under every convention; n = 1 gene each side". C_L1 in D 126 → 112 (single, together); in P 87 (single), 90 (together). Spec REFINE definition amended |
| 2 | P universe rule; group tables missing P\|other and TBC1D26 rows (important) | **Fixed.** `lo_corrected_tables.py` writes PKD1/DHX40 `P\|other` rows and TBC1D26 → MCL24 and asserts both tables equal the universe labels (`IS/corrected_tables.out`). P-universe rule stated in §2 and in the spec. The verifier's 32-gene rule gives **P > D** (246/246 vs 246/255); the pre-audit tables gave a tie at 1.000. Both are in the §6 variant table, not framed as "drop one gene" |
| 3 | The order "D ~ P > C_L1 > C_mid > C_fine, transitive" was mostly built into the construction (critical) | **Fixed.** §0/§6 rewritten: P~D undetermined; C below P/D by construction, with the measured content "0 clade pairs split"; within-C orders by construction; TBC1D3 literature-C comparisons NA (the tournament verdict is now NA for 0-pair sides, `lo_analysis.verdict`); group-level n in every cell; "transitive" dropped and a one-common-gene-set tournament added (28 genes: P = D, 207/207; verdicts identical to the full universes; no intransitive triple); the as-built C_mid_ab is in the main tables (C_fine above C_mid_ab 0.267 vs 0.250) |
| 4 | Principled alternatives break the tie (important) | **Fixed.** §6 variant table with 9 variants plus leave-out. §6ko AN-2 ("TBC1D3–DHX40") and the §6ko qualification (superfamilies) are cited. aa ≥ 0.50 components are reported but **not adopted as the P layer**: it would be a post hoc rule change, and it does not resolve P~D. Components give **P > D** (USP6 and USP32 join the TBC1D3 group; D pairs neither), while the aa ≥ 0.50 MCL clusters give D > P. The audit's D > P came from the MCL reading only. Reproduced: leave-two-out 66 / 66 / 2,146 |
| 5 | D is not §0★★ clause 2/4; failed guard; GENCODE behaviour unmentioned (important) | **Relabel route** (the finding's second option): "E1" in the headline rows, plus the failed-guard and GENCODE sentence at the top. In addition, the real GENCODE E1 catalogs were mapped onto U (`IS/gencode_map.tsv`). GENCODE puts PKD1 in the NPIP group (confirmed), and also CLN3 and USP6. Verdict **D > P** for NPIP (190/231 vs 190/190) and pooled, robust to a ≥ 0.5 mapping filter. No clause-2/4 rebuild on chr1/4/10/16/17/18 (not attempted in this correction pass) |
| 6 | Truth recall conditioned on the prediction; HGNC 2227 is superfamily-level; F rewards singletons (important) | **Fixed.** New (B) member-anchored, layer-independent scores (§7): HGNC 2227 has 57 symbols, 54 RefSeq genes. P R 0.128 (P 1.000); D R 0.227 (P 0.714). For Soto, every mapped family gene is already in U, so (A) is not inflated. "Superfamily-level" is stated wherever HGNC precision is quoted. Bipartite F is NA when either side has 0 pairs (C_fine TBC1D3 now P NA / R 0.000 / F NA) |
| 7 | JOIN/REFINE do not follow the spec (U-cut JOIN; attach rule) (important) | **Fixed.** Whole-group JOIN with U re-closed (+37 genes: MCL24's FOXO3B, TBC1D26-AS1, TBC1D27P, TBC1D28, ZNF286A, ZNF286A-TBC1D26, ZNF286B). Scored on the P universe (HGNC F 0.857, reproducing the audit) and on its own universe (0.737). REFINE attach and single are both reported. Spec JOIN definition amended |
| 8 | EXPR nesting and T2 guaranteed by construction; T2 needs E_M ⊆ E_L (important) | **Fixed.** EXPR(L) ⊆ L and co-membership T2 are labelled by construction. The precondition was measured: it **holds for all 30 refinement pairs**, so the layer-edge T2 checks are also guaranteed. The precondition and the counterexample were added to the spec and to `docs/seeded_family_definition.md` §0★★ clause 6. Threshold sweep t = 1..5 for 3 read modes; the only splits are at t ≥ 1 (unique, unique_mr). Distinct EXPR groups reported (9 at any ≥ 3) |
| 9 | "TBC1D3 has no supported clade" contradicts clause 5 (important) | **Fixed.** The ordered C layer is now the clause-5 split system (C_tree_top, C_tree_min; the whole hierarchy nests in P and D, §5). The sentence is corrected: C_tree has TBC1D3 subfamilies {B,F,G,H}, {TBC1D3,D,E,K}, {B,G,H}, {B,H}. Clause 5 is cited as the source of the rule. Literature C is kept as the circular column; TBC1D3 cluster1/2 stay cleared there because the literature defines them positionally |
| m1 | JOIN scored on P's coding universe while the text said "own universe" | **Fixed.** Both universes are reported (U-restricted own: HGNC 0.857, Soto 0.688 on TBC1D3). U closure and the MCL24 genes of a whole-group JOIN are named |
| m2 | TBC1D3 order through vacuous values; within-C by construction | **Fixed** (§0, §6) |
| m3 | Named NPIPB subfamily also incompatible | **Fixed** (§7): conflicts with {B10P,B11,B15,B6-9} (exon 86.6/75); C_fine and C_mid_ab each contain one conflicting group |
| m4 | Wording: outside genes vs edges; "best member edge"; TBC1D26 fold; strand 65/7 | **Fixed** (§2: 4/25/32 outside genes, 18/267/153 gene-level edges; §9: "best edge to a member in the same D group", plus the MCL7 readthroughs' 1.000 edges to MCL0; TBC1D26's raw catalog edges; §8: "65/7 over all 72 reads") |
| m5 | LOC124905656 encloses TBC1D3P6 | **Fixed** (§3 #1, §4): the rationale is limited to readthroughs; the member count is 45 (TBC1D3 18) if enclosed records are folded; flagged per member (`span_inside_member_record`) |
| m6 | EXPR(L) truth agreement missing | **Fixed** (§8): one line explaining why it is not reported |
| m7 | PC2 component figures disagree with the P verifier | **Fixed** (§2): shipped greedy order gives 4,522 / 965; stable file order gives 4,528 / 968 (reproduced). Outside neighbours are 12 under both orders; the verifier's 10 was **not reproduced** |
| m8 | k acceptance is a stopping rule; stability beyond k = 4; round 1-5 provenance | **Fixed**: called a stopping rule; stability k = 2..8 recomputed (`IS/P_stability_plain.tsv`); `light/scripts/README.md` records the removed `round N` command (argv not saved) |
| m9 | EXPR primary should be any-overlap; unique inconsistent with the member rule; no sweep | **Fixed**: any is primary; unique_mr added (NPIP 11 at t ≥ 3); sweep table |
| m10 | Member reconciliation with §6jg; fold inheritance vs §6ks; C_mid truth circular twice | **Fixed** (§4 table; fold-inheritance convention stated and contrasted with §6ks; C_mid vs named ∨ L2 marked "identical by construction", not scored) |
