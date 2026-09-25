# Nested edge-test lattice on NPIP and TBC1D3 (human T2T-CHM13)

2026-09-16. This report tests the user-approved "nested edge-test lattice" design (16:39) on the data of the NPIP/TBC1D3
layer-order study.

- **Background:** `docs/seeded_family_definition.md` §0★★ (clauses 2, 4, 5, 6 and the 2026-09-16 clause-6 scope note) and
  §0★★★ (the lattice definition this run instantiates); `docs/superpowers/specs/2026-09-16-family-layer-order-design.md`
  (Definitions, AUDIT AMENDMENT); `bench/LAYER_ORDER_NPIP_TBC1D3.md`.
- **Results:** `/mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3/lattice/`, written `LAT/` below.
- **Scripts:** `bench/layer_order/npip_tbc1d3.py lattice-*`, with the library `bench/layer_order/lattice_common.py`
  (wave 7, 2026-09-24; the old `lattice_*.py` files are at git tag `notebook-2026-09-24`; §11 maps old to new).
- **Traceability:** every table is regenerated from the result files by `npip_tbc1d3.py lattice-report`, which writes
  `LAT/report_tables.md`. Sections there are named T-…, and each table below cites its section. Numbers are formatted from
  unrounded values.
- **Status:** descriptive. Not pre-registered. Both families were used to develop the DNA rules (§6js). Not committed.
- **Corrected version (correction pass after the audit of the 17:12 version).** The 17:12 report, its result files and
  its scripts are kept in `LAT/pre_correction_1712/`. §12 lists every audit finding and what changed. Three changes
  bring the tests closer to the definition (`docs/seeded_family_definition.md` §0★★★.1):
  1. **Unrounded thresholds.** The 17:03 edge table stored identities at 4 decimals, which let 2 edges below 0.98 into L3.
  2. **L1 checks the target side, as the shipped clause 2 does.** The hit must overlap v's exons, and spliced pairs
     must pass the shipped strand check (`bench/denovo_shared_def.py`, archived: `notebook-2026-09-19:archive/bench/denovo_shared_def.py`). The 17:12 L1 is kept as the variant
     "no v-exon/strand".
  3. **L3 uses the single-record w_98 (gap-excluded identity), not a pooled identity.** Gap exclusion is a choice. It
     is the convention the 17:12 run fixed before any result. Pooled identities are kept as variants and are not
     monotone.

## 0. Summary

**Levels.** Each level is the connected components of the edges passing its test. Tests, level names, theorems (T1,
T2, T3) and hypotheses (H1–H3) are those of `docs/seeded_family_definition.md` §0★★★.1–.2. L1–L3 each add a conjunct to
the level above, and L0 is L1 OR the protein test, so t_3 ⇒ t_2 ⇒ t_1 ⇒ t_0. **Where this run approximates a test or
makes a choice the definition leaves open, the bullet says so** (details in §1).

- **L0 superfamily (guided-only):** protein OR DNA. *Approximation:* the protein test is the shipped §6ko greedy HSP
  cover, not the definition's union cover (fails H3), and only 4,430 of 20,088 proteins were searched.
- **L1 family:** DNA, clause 2. *Approximation:* the records are all-vs-all gene-body alignments inside two
  chromosome-trio catalogs, not genome alignments. The gene-body disjunct uses the shipped greedy chains (fails H3). The
  exon disjunct is a proxy, because no spliced-transcript alignment exists. As in the definition, the hit must overlap
  v's exons, and spliced pairs must pass the shipped sense (strand) check.
- **L2 shared-exon unit:** L1 AND f_ex ≥ 0.30. *Choice:* f_ex is taken over E1 records only (≥ 300 bp, identity ≥ 0.70),
  as in the shipped code; the definition takes all DNA records.
- **L3 ≥ 0.98 identity unit:** L2 AND w_98 ≥ 0.98, from a single record, as defined. *Choice:* gap-excluded identity. It
  is not a subfamily call; clause 5 stays the subfamily definition.
- **All levels, choice:** pairs whose annotated intervals intersect (same-locus) are excluded from every test, the
  protein test included. The definition only drops DNA records that land on the query's own locus.

"Members" are the 46 annotated NPIP/TBC1D3 records of the prior study (`light/members.corrected.tsv`): NPIP 27,
TBC1D3 19. The uncorrected `light/members.tsv` has 44.

| question | result | source |
|---|---|---|
| Edge table | **432,211 gene pairs** over V (8,070 genes): 344,800 with blastp HSPs, 86,836 with gene-body PAF records, 2,637 with S2 edges. **Re-implementations reproduce the dumps and the shipped code exactly:** E1 graphs 12,629 of 12,629 edges (0 weight mismatches); S1 graphs 9,780 of 9,780; the §6ko edge set 216,024 (symmetric difference 0); `light/P.edges.tsv` 265 of 265; the shipped gene-body chain finder on 78,423 of 78,423 pairs; the shipped v-exon-overlap and strand loop of `denovo_shared_def.py`, run on those chains, on 78,423 of 78,423. The 17:03 build is reproduced on all 432,211 rows: 0 mismatches over 55 columns once the new L1 requirements are dropped | `LAT/edges_build.out`, `LAT/check_c2.out`, `LAT/corrections_pass2/cmp_build.out`, T-edges |
| Can clause 2 (0.80 / 0.50) be evaluated from the E1 dump? | **No; only approximately, from the PAF behind it.** The dump holds one weight = identity × coverage (longer exon-union denominator, records ≥ 300 bp at ≥ 0.70, exon gates), which is not clause 2. **Gene-body clause:** the shipped chains, identity ≥ 0.80, aligned ≥ 0.50 of the shorter body, and the chain's target span overlaps v's exons. **Exon clause:** a proxy, because no spliced-transcript alignment exists: records at identity ≥ 0.80 whose target overlaps v's exons must cover ≥ 0.50 of u's exon union. Both clauses apply the shipped strand check when both genes are spliced | §1 |
| How far the L0 closure goes | **Primary closure of U: 5,769 genes in 16 hops.** This is a lower bound: 1,288 of its proteins were never searched, and 3,455 of its genes lie outside both E1 catalogs. Protein edges alone reach 4,576 genes; primary clause-2 edges alone reach 458. **With the union-cover protein test of the definition: 5,802, and 33 of those lie outside V**, so V is not closed under that form. The 17:12 closure (no v-exon/strand requirement) was 6,416 | T-closure |
| Level sizes (primary) | **L0:** NPIP and TBC1D3 share **one group of 5,764 genes** (42 members). **L1:** NPIP 341 (26 members); TBC1D3 97 (16). **L2:** NPIP 142 (26); TBC1D3 24 (16). **L3:** NPIP 55 (25) plus {NPIPB14P, LOC124907811}; TBC1D3 11 (11) plus {TBC1D3P3, TBC1D3P4}; LOC100420311, TBC1D3P5 and TBC1D3P7 are singletons. 4 members are singletons at every level: NPIPB1P, TBC1D3P6, LOC124905656, LOC100420289 (no protein, outside both catalogs). The 17:12 L1 gave 6,412; 604 · 897; 159 · 31; 73 · 16 | T-levels |
| T1 and T2 on real data | **0 violations in 302 checks.** The checks cover 17 variants, 3 expression sets and both 3-truss expression views. They hold by construction, so a nonzero count would have meant a bug. **Non-vacuity (primary):** of the multi-gene blocks, L1 splits 22 of 115 L0 blocks (1 holds members), L2 splits 84 of 192 L1 blocks (2), and L3 splits 207 of 270 L2 blocks (2). The one member-level L0 → L1 split, NPIP vs TBC1D3, is forced by the data. DNA evidence exists only inside one chromosome trio (0 cross-trio L1 edges), so L1–L3 cannot join the two families, and L1–L3 sizes are lower bounds for copies in other trios | T-sanity, T-crosstrio |
| Chaining | **L0 joins NPIP to TBC1D3** by a 7-edge path: BNIP3P16 (a 2,051-bp span-exon pseudogene aligned inside NPIPB9) → ZNF28 → ZNF232 (protein edge, aa 0.553) → KRT17P4 → LGALS9B (249,716 bp) → TBC1D3. **No protein-side strengthening tested separates the families.** (P ∧ aa ≥ 0.50) ∨ D still joins them (1,017 genes), as does the union-cover protein test. Only replacing D by E1 at 0.80 / 0.50 separates them (94 / 5,003). **L1 hubs** come through the gene-body disjunct: short span-exon pseudogenes (VN1R91P, 889 bp, degree 95; BNIP3P16, 63) and a 250-kb record that contains copies (LGALS9B, 41). The intron-landing hubs of the 17:12 report came from the missing v-exon requirement: TM2D3, MIR4713HG, LOC101926889 and LOC105372586 keep 6, 2, 3 and 6 of their 108, 60, 104 and 69 edges | T-L0path, T-L0variants, T-hubs, §5 |
| What L2 removes from TBC1D3, judged against Soto | L2 removes DHX40, DHX40P1, TBC1D3P1-DHX40P1 and CA4; TBC1D26/28 and RNFT1-DT are already apart at L1. It keeps USP6, USP32, USP32P1-4, TBC1D29P and LOC100420408. **Soto families:** USP32P1-3 are ID_60 and USP32/USP32P4 are ID_91; DHX40/DHX40P1 are ID_341, separate from TBC1D3's ID_468. **Soto pair precision inside the L2 TBC1D3 group is 0.392 (60 of 153 pairs)** | T-chaining, T-soto-neighbours, T-truth-ingroup |
| What stays with NPIP, judged against Soto | PKD1 and PKD1P1/P2/P3/P6 stay at L1–L3, linked through PKD1P–NPIP readthrough records. **Soto puts PKD1, PKD1P1-3, P6 and the readthroughs LOC131696449, PKD1P3-NPIPA1 and PKD1P5-LOC105376752 in the same family as NPIPA1, NPIPA9 and the member PKD1P6-NPIPP1 (ID_149), so by Soto PKD1 joining NPIPA1/A9 is not a false merge.** By Soto the L3 NPIP group still merges six families: ID_149 (11 labelled genes) with ID_154 (16, most NPIP members, including the anchor NPIPB2), ID_41 (SMG1P, 10), ID_302 (BOLA2, 3), ID_152 (1) and ID_153 (1). The L2 group also holds ID_28 (PLA2G10*P, 11). Soto pair precision inside the NPIP group is 0.187 (L2) and 0.259 (L3) | T-truth-ingroup, T-soto-members |
| Truth agreement: bipartite F (one-to-one Jaccard) on the prior study's gene sets; in-group precision next to it | **Soto, TBC1D3 (16 member genes):** 0.688 → 0.750 → 0.875 → **1.000** (L0 → L3). **Soto pair precision of the whole anchor group:** 0.015 → 0.078 → 0.392 → 1.000. **On V:** 0.015 → 0.115 → 0.244 → 0.440 (all 550 Soto-labelled genes of V). **Soto families are SD98 duplications with a shared-exon map-back, the same conventions L2 and L3 test, so agreement at L2/L3 is partly by construction.** **Soto, NPIP (32):** 0.531 at L0–L2, 0.594 at L3; D/E1 MCL scores 0.806 on its 31 genes (lattice L3 0.581). **Literature NPIPA\|NPIPB (21 C_tree leaves):** 0.634 at every primary level; C_tree_top 0.765, and L3 with gap-inclusive w_98 also 0.765. **Literature L2:** primary 0.238, gap-inclusive w_98 0.667, C_tree_min 0.833. **Clause 5 was developed on NPIP (§6jp–§6js), so C_tree is not an independent comparator** | T-truth-U, T-truth-ingroup, T-truth-V |
| Threshold filtration (single linkage on fixed evidence) | On the primary L3 field every literature or clause-5 boundary lies between 0.978354 and 0.999014. **On that field** (w_98 gap-excluded, shared-exon ON): NPIPA and NPIPB stop sharing a group above 0.980559 (first grid value 0.99). NPIPA is an exact group in (0.980559, 0.989897]. **NPIPB never is:** NPIPB14P joins the rest of NPIP at 0.978354. B3-5 appears in (0.997341, 0.998666], B6-9 in (0.994729, 0.996797] and B12/13 in (0.997341, 0.999014]. A6-9 and TBC1D3 {B,F,G,H} never appear. **The gap convention moves the boundaries:** with gap-inclusive w_98, NPIPA is exact at the grid value 0.98 | T-appearance, `LAT/filtration.txt` |
| Triangle-support (3-truss) variant | **The theorems still hold** (0 violations). **Sizes, primary → 3-truss:** L1 NPIP 341 → 230 (33% removed) and TBC1D3 97 → 63 (35%); L2 142 → 64 and 24 → 23; L3 55 → 37 and 11 → 11. It dissolves every 2-gene component: 87, 120 and 76 at L1–L3. Two of them held members: {TBC1D3P3, TBC1D3P4}, **a true 2-copy family (both Soto ID_469)**, and {NPIPB14P, LOC124907811}. PKD1 and PKD1P1/P2/P3/P6 stay with NPIP at every triangle level. The two large, dense families here cannot measure the cost for 2-copy families | T-sizes, T-triangle, §8 |
| Expression view (testis) | **Any-overlap reads ≥ 3:** 13 NPIP and 0 TBC1D3 members. There is one member-holding component per level: 16, 16, 16 and 15 genes, each with 12 members. **Any-overlap is not copy-resolved** (reads of any MAPQ). **Unique reads ≥ 3:** 9 members; components of 8, 8 and 8 genes, then 6 + 2 at L3. **At reads ≥ 1 the L0 view separates NPIP (48 genes) from TBC1D3 (3,022).** On the observed graph, every L0 path between them passes through a copy with 0 any-overlap reads. The protein search and the DNA catalogs are truncated, so an expressed bridge through unobserved edges is not ruled out | T-expr |

**What this says about the design, on these two families (descriptive):**

- **Nesting holds; T3 is not established here.**
  - T1 and T2 hold by construction. Each finer test implies the coarser one, so E_{k+1} ⊆ E_k, which supplies the
    E_M ⊆ E_L condition the audit found missing.
  - T3(b) is not established for this operationalisation (§4). *More copies* needs H1 and H2, and H2 fails: E-values
    depend on database size, and the minimap2 `-N`/`-p` caps depend on the target set, which here is the node set
    itself. *More homology* also needs H3, and H3 fails: the shipped t_P and gene-body chains are greedy (the pooled L3
    identity variants are not maxima over records either).
- **The cost is as large as feared.**
  - Components chain into superfamily-scale L0 (one 5,764-gene group) and into large L1 groups: NPIP 341 genes,
    including 74 ZNF-named and 20 BNIP3P records.
  - The shared-exon test (L2) removes most L1 non-members (NPIP 315 → 116, TBC1D3 81 → 8), but not all co-duplications.
    TBC1D3 keeps USP6 and the USP32 family (Soto in-group precision 0.392). NPIP keeps SMG1P, BOLA2, PLA2G10*P and 40
    ZNF-named records (0.187).
  - The 3-truss removes 33% of the NPIP and 35% of the TBC1D3 L1 group and keeps the PKD1 co-duplication, which Soto
    calls one family with NPIPA1/A9. It dissolves the true 2-copy family {TBC1D3P3, TBC1D3P4}.
- **Whether a fixed identity cut gives subfamilies depends on the gap convention.**
  - With gap-excluded w_98 (primary), L3 keeps 25 of the 26 connected NPIP members in one group; only NPIPB14P leaves.
  - With gap-inclusive w_98, L3 splits NPIP into three groups, with B2, B11 and B14P as singletons:
    - A (+PKD1P6-NPIPP1), 20 genes: 8 members and 12 non-members (PKD1P1/P2/P3/P6, the 4 PKD1P readthroughs and 4 LOC
      records; PKD1 itself is not in it);
    - B3-5/B12/B13 (+2 LOC members), 13 genes: 7 members and 6 non-members (SMG1P1/P3/P4 and 3 LOC records);
    - B6-9/B10P/B15 (+2 LOC members), 8 genes, no non-members.
  - The gap-inclusive variant was seen after the primary was fixed, so it is not adopted. Fix the gap convention on
    held-out data.
  - Under either convention several literature groups never appear as an exact single-linkage group at any cut: NPIPB,
    TBC1D3 {B,F,G,H}, and A6-9 under the primary.

## 1. What was built (operationalisation)

**Nodes.** RefSeq gene and pseudogene records (`chm13v2.0_RefSeq_full.gff.gz` via `light/work/refseq/genes.tsv`), at
their annotated span. Nothing is read-extended; the node set is frozen for the whole run.

- Readthrough records are nodes: 56 in V, one of them a member (PKD1P6-NPIPP1). A variant removes them (§3, §5).
- Overlapping records are not merged.

**Graph.** V is U closed under the union of the L0 edges of every operationalisation below, including the 17:12 L1
without the v-exon/strand requirements: 8,070 genes, identical to the 17:03 V. V contains the primary closure (5,769
genes). **One variant is cut at the edge of V: the union-cover protein test.** Its closure of U has 5,802 genes, 33 of
them outside V (§2). The 5 U genes outside its joint group on V have no edges, so its joint group has 5,802 − 5 = 5,797
genes, of which the 5,764 in V are reported. Every other variant's components lie inside V.

**Tests** (`lattice_common.tests`). A missing attribute makes its clause unsatisfiable; nothing is imputed. All
comparisons are on unrounded values.

| level (§0★★★.1 name) | test | operationalisation (approximations and choices in bold) | monotone in the evidence (H3)? |
|---|---|---|---|
| L0 superfamily (guided-only) | not same-locus ∧ (t_P ∨ t_1) | t_P = §6ko edge as shipped (`bench/protein_families.edges_from`): blastp e ≤ 1e-5 (search cutoff), coverage ≥ 0.30 of the longer protein, any identity. **Greedy HSP selection, not the definition's union cover** (genome-wide 216,024 vs 229,831 qualifying pairs; union-cover closure row in §2). **Evidence exists only for pairs where at least one protein was searched (4,430 of 20,088)** | **No:** the greedy HSP selection fails H3, and E-values depend on the database size (H2) |
| L1 family | not same-locus ∧ (GB ∨ EX), shipped strand check | **Records: all-vs-all gene-body alignments inside one catalog (chr15/17/22 or chr16/19/20), not the definition's genome alignments; no spliced-transcript record exists.** **GB:** chains built exactly as `bench/guided_pipeline.gene_body_chains` (gap ≤ query body, span ≤ 2 × query body, query order) on the pair's gene-body PAF records, both directions. A chain passes at identity ≥ 0.80 with aligned query bp ≥ 0.50 × min(body u, body v), **and its target span must overlap ≥ 1 exonic base of v**. **EX (proxy):** records with identity ≥ 0.80 **whose target overlaps v's exons** cover ≥ 0.50 of u's exon union. **The exon union stands in for u's spliced transcript, and record identity includes introns.** **Strand check** (`denovo_shared_def.py cmd_families`): when both genes are spliced (exon union ≥ 2 blocks; 75,346 of the 86,836 PAF pairs in V), v's strand must equal u's for a '+' hit and differ for a '−' hit. Gene bodies are extracted on the genomic + strand | GB: **no** (greedy chains, pooled identity). EX proxy: yes (a union over records). Both: the minimap2 `-N 50 -p 0.1` caps depend on the target set (H2) |
| L2 shared-exon unit | t_1 ∧ f_ex ≥ 0.30 | The `mcl_families --min-shared-exon-frac` quantity, re-implemented (reproduces the S1 dumps exactly). f_ex = max over single records of min(exonic bp of u inside the record, exonic bp of v inside it) / min(exon-union length u, v). A record without exon features counts as one span exon. **Taken over E1 records only (≥ 300 bp, identity ≥ 0.70): the shipped form, a choice.** The definition takes all records: 8,129 L2 edges instead of 8,034 | yes (a maximum over records) |
| L3 ≥ 0.98 identity unit | t_2 ∧ w_98 ≥ 0.98 | w_98 = max identity over single PAF records of the pair whose shared exonic bases are ≥ 0.30 × min(exon-union length u, v): one record witnesses both the identity and the shared exons (§0★★★.1). **Identity gap-excluded (Σ matches / Σ CIGAR M, the SEDEF `fracMatch` analogue). This is a choice, fixed before any result: it was already the 17:12 primary convention.** Restricting w_98 to E1 records changes 139 of 7,294 passing key pairs genome-wide. Variants: w_98 gap-inclusive; identity pooled over the E1 records (gap-excluded = the 17:12 L3; gap-inclusive); S2 map-back `max_identity` | w_98: yes. **Pooled variants: no** (a new record can dilute the pooled identity) |

**Same-locus.** Two records whose intervals intersect on one chromosome form a same-locus pair: a readthrough and its
component gene, or nested records. 2,218 of the table's pairs are same-locus.

- Their "alignment" is the identity map of shared bases, not homology between two copies.
- They are excluded from every test, the protein test included. The `with_same_locus` variant re-admits them.
- **This is a choice beyond §0★★★.1.** The definition drops only DNA records whose target overlaps the query's own locus,
  and it says nothing about same-locus HSPs.

**Variants reported** (T-levels):

- `with_same_locus`.
- `L0=(P_and_aa>=0.50)_or_D`: the §0★★★.2 T3(a) "strengthen one disjunct" stack.
- Three L1 forms that relax the new requirements: v-exon overlap without the strand check; v-exon overlap on the exon
  proxy only; neither (the 17:12 L1).
- `L1=c2x_extrapolated`: the guided finder's own denominator, min(query body, extrapolated target span), with the v-exon
  and strand requirements.
  - On a gene-body PAF the extrapolation is clipped at the target body's end, so a 600-bp overlap at two body ends passes
    the denominator: NPIPA8–PKD1P1, 600 of 600 bp identical.
- `L1=E1_as_built`: the D layer's graph.
- `L1=E1_at_0.80/0.50`: E1 edges with identity ≥ 0.80 and E1 coverage ≥ 0.50.
- Four L3 identity forms: w_98 gap-inclusive; pooled gap-excluded; pooled gap-inclusive; the S2 map-back identity.
- `17:03_tests_exact`: the 17:12 report's tests (no v-exon/strand requirement, pooled gap-excluded L3), unrounded.
- `triangle`: the 3-truss of each primary level (and of `17:03_tests_exact`).
- `no_readthrough_nodes`.

## 2. Edge table provenance (`LAT/edges.tsv`, one row per pair inside V with any evidence; floats at full precision)

| attribute group (columns) | source | exact or approximate; missing ⇒ |
|---|---|---|
| protein: `p_aa_identity`, `p_cov_longer`, `p_weight`, `p_max_bitscore`, `p_qualifies_6ko`, `p_aa50`, `p_directions`, `p_searched_a/b`; `p_cov_union_longer`, `p_qualifies_union` | `light/work/P/blastp.tsv` (1,436,646 HSP lines; outfmt `qseqid sseqid nident length qstart qend sstart send bitscore`); `proteins.index.tsv`; `searched.txt` | Exact: the shipped greedy HSP order reproduces all 216,024 §6ko edges. The union-cover columns merge every HSP interval on the longer protein (229,831 qualifying pairs genome-wide, 13,807 union-only, 0 greedy-only). They feed one closure row and no level. **`p_evalue` was not saved per HSP;** the table records "≤ 1e-5 (search cutoff)". A pair with no HSP cannot satisfy the protein clause |
| DNA, E1 as recorded: `d_e1_identity`, `d_e1_identity_gapexcl`, `d_e1_cov_longer` and its denominator, `d_e1_gate_exonic`, `d_e1_edge`, `d_e1_weight`, `d_e1_records` | `o1_falsemerge/human2/genes.asm20.paf` (chr15/17/22) and `lit/aj_ho/refseq/all.paf` (chr16/19/20). Both are `minimap2 -x asm20 -c -X -N 50 -p 0.1` all-vs-all runs on gene-body spans (`human2/run_chunks.sh`, `bench/node_graph_mcl.py`). Exon unions come from `light/work/refseq/exons.tsv` and `lit/aj_ho/refseq/nodes.tsv` | **Exact** for E1: the weight identity × cov_longer reproduces all 12,629 dump edges. The dump itself stores only the weight. Pairs across the two catalogs, or with a gene on another chromosome, have no PAF, so t_1–t_3 are unsatisfiable |
| DNA, clause 2: `d_c2_genebody` (+ `d_c2_gb_best_frac`, `d_c2_gb_chain_identity`), `d_c2_exon` (+ `d_c2_exon_best_frac`), `d_c2_approx`, `d_both_spliced`; variants `d_c2nostrand_approx`, `d_c2exontgt_approx`, `d_c2loose_*`, `d_c2x_*` | the same PAFs | **Approximate** (§1). Two checks against shipped code (`LAT/check_c2.out`): (1) the raw chains agree with `gene_body_chains` on 78,423 of 78,423 pairs; (2) the `cmd_families` edge loop of `denovo_shared_def.py` (ExonIndex hit on v's exons, strand check), run on those chains, agrees with `d_c2x_genebody` on 78,423 of 78,423. The primary uses the same target/strand code with the shorter-body denominator. `-X` keeps one direction per pair; the other is rebuilt by swapping coordinates. Pairs beyond minimap2's caps are absent. The `d_c2loose_*` columns equal the 17:03 primary columns on every row |
| shared exon: `d_shared_exon_bp`, `d_shared_exon_frac`, `d_shared_exon_denominator`; `d_shared_exon_frac_allrec` | the same PAF records (primary: ≥ 300 bp, ≥ 0.70; `allrec`: all records) | **Exact** mcl_families quantity (S1 dumps 9,780 of 9,780). **It is taken from the single best record**, so a copy pair whose alignment is split into many records can score low: TBC1D3–TBC1D3P1-DHX40P1 scores 0.151 at pooled identity 0.983. 2,358 PAF pairs have no qualifying record: their fraction is 0 and their E1 identity is NA |
| L3 identity: `d_w98_gapexcl` (primary), `d_w98_gapincl`; pooled `d_e1_identity_gapexcl`, `d_e1_identity` | the same records | w_98 is NA when no single record witnesses shared-exon ≥ 0.30 (9,962 pairs have a gap-excluded w_98). Short records can reach 1.000: LOC100190986 (2,453-bp lncRNA) vs NPIPB5 at w_98 1.000 |
| S2: `s2_max_identity`, `s2_n_shared_exons`, `s2_n_projected_exon_pairs`, `s2_same_locus` | `heavy/S2.edges.tsv` (3,971 edges; 321 genes of the SD98 closure, all mapped) | Only 2,637 pairs in V have one. Missing ⇒ the S2 identity variant is unsatisfiable |
| clade support: `ctree_family`, `ctree_smallest_common_cluster`, `ctree_support`, `ctree_same_top_cluster` | `lattice_common.c_tree` on `light/C.supported_clades.tsv` (clause-5 split system) | Annotation only; used by no test. 246 pairs annotated (30 leaves) |

**Closure** (T-closure; `LAT/closure.tsv`).

| closure of U (68 genes) | genes | new genes per hop | never-searched proteins | genes without a §6ko protein | genes outside both E1 catalogs | outside V |
|---|---|---|---|---|---|---|
| **primary L0** | **5,769** | 38, 156, 1,152, 807, 803, 711, 798, 833, 208, 86, 54, 21, 20, 8, 5, 1 | 1,288 | 832 | 3,455 | 0 |
| union of all reported operationalisations = V | 8,070 | 89, 945, 2,090, 2,321, 1,248, 691, 476, 104, 28, 7, 3 | 1,909 | 2,433 | 3,606 | 0 |
| protein edges only (greedy cover) | 4,576 | 14, 16, 79, 328, 260, 503, 742, 787, 1,311, 400, 41, 27 | 965 | 34 | 3,316 | 0 |
| primary clause-2 DNA edges only | 458 | 26, 141, 82, 80, 40, 13, 8 | 135 | 268 | 5 | 0 |
| union-cover protein ∨ primary clause 2 | 5,802 | 38, 156, 1,164, 880, 851, 746, 1,187, 366, 169, 76, 48, 23, 19, 7, 3, 1 | 1,292 | 832 | 3,483 | **33** |
| 17:12 L0 (no v-exon/strand requirement) | 6,416 | 69, 695, 1,753, 1,688, 971, 582, 270, 91, 61, 159, 7, 2 | 1,464 | 1,276 | 3,487 | 0 |
| 17:12 clause-2 DNA edges only | 1,506 | 57, 253, 346, 405, 230, 93, 30, 22, 2 | 390 | 949 | 5 | 0 |

- **The closure is truncated in two ways**, so every L0 size below is a lower bound:
  - **Protein:** an edge between two never-searched proteins is unobservable. The TBC-domain component was already
    ≥ 4,522 genes in the prior study.
  - **DNA:** evidence exists only inside one chromosome trio. There is no chr16↔chr17 DNA edge, and the 3,606 genes of V
    outside both catalogs have none.
- **DHX40 is an example.** Its protein was never searched, so its L0 edges are DNA edges only.

## 3. Levels (primary), groups holding members, pulled-in non-members

Primary edge counts per level: 224,107, 10,276, 8,034 and 5,754 (T-edges). Group sizes are written as size (members /
non-members) (T-levels; `LAT/member_groups.tsv` lists every group; `LAT/groups.tsv` gives the label of every gene per
level and variant).

| level | NPIP group(s) | TBC1D3 group(s) | member singletons |
|---|---|---|---|
| L0 | **5,764 (42 / 5,722), one group for both families** | — | NPIPB1P, TBC1D3P6, LOC124905656, LOC100420289 |
| L1 | 341 (26 / 315) | 97 (16 / 81) | the same 4 |
| L2 | 142 (26 / 116) | 24 (16 / 8) | the same 4 |
| L3 | 55 (25 / 30); {NPIPB14P, LOC124907811} | 11 (11 / 0); {TBC1D3P3, TBC1D3P4} | the same 4 + LOC100420311, TBC1D3P5, TBC1D3P7 |

**Non-members pulled in** (T-nonmembers)

- **L0:** 5,722 non-members. The group has 4,937 protein-coding genes; 714 are on chr19, 498 on chr1, 467 on chr17 and
  428 on chr15. 1,288 of its genes have never-searched proteins, and 3,451 lie outside the catalogs.
- **L1, NPIP (315 non-members):** chr19 159, chr16 124, chr20 58 (all 341 genes); 74 ZNF-named and 20 BNIP3P records.
- **L1, TBC1D3 (81 non-members):** 96 of the 97 genes are on chr17. They include USP6, USP32, USP32P1-4, DHX40, DHX40P1,
  TBC1D3P1-DHX40P1, CA4, KRT14/16/17 and their pseudogenes, LGALS9/9B/9C/9DP, CDRT15*, FAM106*, SRP68*, MIR4713HG and
  YWHAE.
- **L2, NPIP (116):** the group of 142 has 82 genes on chr16, 57 on chr19 and 3 on chr20.
  - 40 ZNF-named records and 1 BNIP3P (BNIP3P16).
  - PKD1, PKD1P1/P2/P3/P6, and the readthroughs PKD1P3-NPIPA1, PKD1P4-NPIPA8, PKD1P5-LOC105376752, LOC131696449 and
    PDXDC2P-NPIPB14P.
  - SMG1 and SMG1P1-7; PLA2G10CP/EP-KP; BOLA2, BOLA2B, BOLA2-SMG1P6; PDXDC1, PDXDC2P; SLC7A5, SLC7A5P1/P2;
    VN1R81P/82P/91P; MTDHP2/P5; PABPN1P2; SIRPD; LINC01859.
  - 33 further LOC records.
- **L2, TBC1D3 (8):** USP6, USP32, USP32P1-4, TBC1D29P, LOC100420408.
- **L3, NPIP (30):**
  - PKD1, PKD1P1/P2/P3/P6 and the 4 PKD1P readthroughs (PKD1P3-NPIPA1, PKD1P4-NPIPA8, PKD1P5-LOC105376752,
    LOC131696449).
  - SMG1P1-7; BOLA2, BOLA2B, BOLA2-SMG1P6; SLC7A5P2.
  - 10 further LOC records.
- **L3, TBC1D3:** no non-members. The 17:12 L3 group held 5 lncRNA records whose exons align inside TBC1D3 copies. On the
  17:12 route (LOC105371853–TBC1D3D: shared exon 1.000, identity 0.982) the lncRNA and the copy are both annotated '+'
  but align on '−'. That makes the pair antisense, and the strand check removes it.

**Operationalisation variants** (group sizes; members in parentheses where not all connected members are in the group):

| variant | L0 | L1 NPIP · TBC1D3 | L2 | L3 |
|---|---|---|---|---|
| primary | one group, 5,764 | 341 · 97 | 142 · 24 | 55 (25) + {NPIPB14P}+1 · 11 + {P3,P4} |
| with same-locus links | one group, 6,260 | 435 · 121 | 148 · 24 | 56 (25) + {NPIPB14P}+5 · 11 + {P3,P4} |
| L0 = (P ∧ aa ≥ 0.50) ∨ D | one group, 1,017 | = primary | = primary | = primary |
| L1: v-exon overlap, no strand check | one group, 6,050 | 378 · 111 | 159 · 31 | 73 · 16 (11) + {P3,P4} |
| L1: v-exon overlap on the exon proxy only | one group, 6,162 | 474 · 725 | 159 · 31 | 73 · 16 (11) + {P3,P4} |
| L1: no v-exon overlap, no strand check (17:12 L1) | one group, 6,412 | 604 · 897 | 159 · 31 | 73 · 16 (11) + {P3,P4} |
| L1 = c2x (finder denominator, v-exon + strand) | one group, 5,825 | 373 · 105 | 143 · 25 | 57 · 11 + {P3,P4} |
| L1 = E1 as built (D graph) | one group, 5,295 | 117 · 39 | 85 · 25 | 74 (25) · 11 + {P3,P4} |
| L1 = E1 at 0.80 / 0.50 | **NPIP 94 · TBC1D3 5,003 (separate)** | 92 · 28 | 71 · 23 (15) | 61 (25) · 11 + {P3,P4} |
| L3 = w_98 gap-inclusive | = primary | = primary | = primary | NPIP split: {A1,A2,A5-9,PKD1P6-NPIPP1}+12 · {B3,B4,B5,B12,B13,LOC124907834,LOC128966608}+6 · {B6-9,B10P,B15,LOC124907807/808}; singletons B2, B11, B14P · TBC1D3 {9 coding} · {P1,P2} · {P3,P4} |
| L3 = pooled gap-excluded (17:12 identity, unrounded) | = primary | = primary | = primary | 55 (25) + {NPIPB14P}+1 · 11 + {P3,P4} (5,720 edges) |
| L3 = pooled gap-inclusive | = primary | = primary | = primary | NPIP split: {A1,A2,A6-9}+12 · {B6-9,B10P,B15,LOC124907807/808} · {B5}+6 · {B3,B12,LOC124907834,LOC128966608}; singletons A5, B2, B4, B11, B13, B14P, PKD1P6-NPIPP1 · TBC1D3 {9 coding} · {P1,P2} · {P3,P4} |
| L3 = S2 map-back | = primary | = primary | = primary | NPIP split: {A1,A2,A5-9}+7 · {B3,B4,B5,LOC128966608}+6 · {B6-9,B10P} · {B11,B12,B13,LOC124907834} · {B15,LOC124907807/808}; singletons B2, B14P, PKD1P6-NPIPP1 · TBC1D3 {9 coding} · {P1,P2}; P3, P4 singletons |
| 17:12 tests, unrounded | one group, 6,412 | 604 · 897 | 159 · 31 | 73 · 16 (11) + {P3,P4} (5,868 edges; the 17:12 report said 5,870) |
| no readthrough nodes | one group, 5,706 | 328 · 94 | 123 · 24 | 37 (24) + {NPIPB14P}+1 · 11 + {P3,P4} |

- **w_98 vs pooled (gap-excluded) at L3** (T-L3variants). Of the 8,034 primary L2 edges, 5,710 pass both, 44 only w_98
  and 10 only the pooled identity. The member groups are identical.
- **The gap-inclusive identity variants match the literature groups better (§6).** That was seen after the primary was
  fixed, so neither is adopted.

## 4. T1 / T2 sanity (`LAT/sanity.tsv`; T-sanity)

"Coarse blocks" are the multi-gene blocks of the coarser partition. "Split" means the finer partition actually divides
them, which shows the check is not vacuous.

| check | checks | groups checked (all, singletons included) | violations | coarse blocks ≥ 2 genes | of which split | member-holding split |
|---|---|---|---|---|---|---|
| T1: components of G_j refine G_i, i < j (17 variants) | 102 | 727,594 | **0** | 18,131 | 12,463 | 156 |
| 3-truss components refine plain components (primary, 17:12 tests) | 8 | 48,016 | **0** | 1,308 | 996 | 16 |
| T2b: G_k[X] refines G_k restricted to X (any ≥ 3 / any ≥ 1 / unique ≥ 3) | 16 / 16 / 16 | 37,677 / 56,096 / 27,719 | **0** | 917 / 1,574 / 546 | 212 / 285 / 150 | 16 / 21 / 16 |
| T2a: G_{k+1}[X] refines G_k[X] | 12 / 12 / 12 | 34,236 / 52,012 / 25,257 | **0** | 1,097 / 1,559 / 803 | 637 / 834 / 495 | 6 / 21 / 5 |
| 3-truss views: T2a and T2b for both views, and comp(Δ(G_k[X])) ≤ comp(Δ(G_k)[X]) | 36 / 36 / 36 | 95,795 / 144,478 / 69,657 | **0** | 1,474 / 2,245 / 1,216 | 837 / 1,070 / 761 | 26 / 54 / 28 |

**Primary, per level pair** (T1; multi-gene coarse blocks / split / member-holding split): L1 in L0 115 / 22 / 1; L2 in
L1 192 / 84 / 2; L3 in L2 270 / 207 / 2. **T2b at reads ≥ 3:** L0 43 / 13 / 1, L1 63 / 23 / 1, L2 69 / 14 / 1, L3 34 / 5 / 1.
The single member-holding L0 block that L1 splits is the joint NPIP–TBC1D3 group. **That split is forced by the data:
0 L1 edges cross chromosome trios**, because the DNA catalogs are per trio (T-crosstrio).

**Why each is guaranteed.**

- **T1 and T2a:** each finer test implies the coarser one (L1–L3 add a conjunct; L0 = L1 OR protein), so
  E_{k+1} ⊆ E_k and E_{k+1}[X] ⊆ E_k[X].
- **T2b:** an induced subgraph's edges are a subset of the full graph's edges.
- **3-truss:** truss(H) ⊆ truss(G) whenever H ⊆ G, and truss(G[X]) ⊆ truss(G)[X]. The truss is the largest subgraph in
  which every edge lies in a triangle.
- The E_M ⊆ E_L precondition of the prior study's layer operators is automatic here.

**T3 is not established for this operationalisation.**

- **T3(a)** (a stronger nested stack only splits) is a theorem that needs no hypothesis on nodes or evidence. The data
  hold one instance, the stack with L0 = (P ∧ aa ≥ 0.50) ∨ D: each of its 4,530 L0 groups lies inside one primary L0
  group (the NPIP–TBC1D3 group shrinks from 5,764 to 1,017 genes), and its L1–L3 labels equal the primary labels
  (`LAT/groups.tsv`).
- **T3(b)** (more copies or more homology only merges; §0★★★.2): *more copies* needs H1 and H2; *more homology* also
  needs H3. This run's tests do not meet H2 or H3.
  - **H3 fails** for t_P and for the gene-body chains. t_P selects HSPs greedily: 399 of 5,000 random instances lose an
    edge (§0★★★.2). The chains are greedy with a pooled identity: 203 of 5,000. The pooled-identity L3 variants fail H3
    as well. The primary L3 (single-record w_98), f_ex and the exon proxy are maxima or unions over records and meet H3.
  - **H2 fails.** blastp E-values depend on database size. The minimap2 `-N 50 -p 0.1` caps depend on the target set,
    and the target set is the node set itself (all-vs-all gene bodies).
  - **H1 holds** for this frozen annotated node set.
- **The 0-violation T2b count is not a T3 test.** It runs on one frozen evidence table, so it cannot show that adding
  copies or evidence only merges.
- **What a real T3 test needs:** the monotone test forms plus new evidence. One example is searching the 1,288
  never-searched proteins of the primary closure with a fixed `-dbsize`. No such test has been run.

**Which records become nodes** drives the PKD1 merge.

- The 56 readthrough records in V are annotated records with fixed extents. So this is a node-set question (the proposed
  readthrough rule of §0★★★.1), not a node-extent question.
- Without them, PKD1, PKD1P1/P2/P3/P6 and 4 lncRNA LOC records form their own 9-gene group at L2 and at L3 (§5).

## 5. Chaining evidence

**L0 joins the two families** (T-L0path, primary). One shortest path, from BFS with sorted neighbours:

1. NPIPB2 – NPIPB9: protein aa 0.604 (coverage 0.919); exon proxy 0.765; E1; shared exon 0.765.
2. NPIPB9 – **BNIP3P16**: a 2,051-bp pseudogene with a single span exon, aligned inside NPIPB9 (gene-body 0.682, chain
   identity 0.878, shared exon 0).
3. BNIP3P16 – ZNF28: gene-body 0.847 (chain identity 0.811).
4. ZNF28 – ZNF232: **protein only, aa 0.553** (coverage 0.418); the genes are on chr19 and chr17, in different catalogs.
5. ZNF232 – KRT17P4: gene-body 0.585 (chain identity 0.839).
6. KRT17P4 – **LGALS9B**: gene-body, chain identity 0.994. LGALS9B is a 249,716-bp record.
7. LGALS9B – TBC1D3: gene-body 0.850, **chain identity 0.907** (pooled E1 identity: 0.907 gap-inclusive, 0.919
   gap-excluded); shared exon 0.

The 17:12 path ran through USP31–USP6 (protein aa 0.312). Its links LOC105372335–LOC101926889 and LOC101926889–USP31
were exon-proxy-only edges with 0 shared exonic bp, and they fail the v-exon requirement.

**Which L0 forms join the families** (T-L0variants; components on V):

| L0 edge form | joined | NPIPB2 group | TBC1D3 group |
|---|---|---|---|
| primary: P ∨ D | yes | 5,764 | 5,764 |
| (P ∧ aa ≥ 0.50) ∨ D | **yes** | 1,017 | 1,017 |
| union-cover P ∨ D (cut at V: 5,797 genes with the 33 outside V, §1) | yes | 5,764 | 5,764 |
| P ∨ D, gene-body disjunct only | yes | 5,659 | 5,659 |
| P ∨ D, exon-proxy disjunct only | yes | 4,998 | 4,998 |
| P ∨ E1 as built (0.70 / 0.30, exon-to-exon gate) | yes | 5,295 | 5,295 |
| P ∨ E1 at 0.80 / 0.50 | **no** | 94 | 5,003 |
| P only | no | 21 | 4,517 |
| D only (= L1) | no | 341 | 97 |
| 17:12: (P ∧ aa ≥ 0.50) ∨ D without v-exon/strand | yes | 2,300 | 2,300 |

- **No protein-side strengthening tested separates the families.** (P ∧ aa ≥ 0.50) ∨ D keeps the same 7-edge path,
  because ZNF28–ZNF232 has aa 0.553.
- **Only the DNA side separates them, and the reason is E1's thresholds, not exon-to-exon evidence.** P ∨ E1 at
  0.80 / 0.50 separates them. P ∨ E1 as built joins them, although it has the same exon-to-exon gate. The difference
  is E1's 0.80 identity and 0.50 coverage on its pooled identity over the longer exon union.
- **The L0 size is a lower bound.** L0 unites protein superfamilies through DNA hubs. The TBC-domain component alone was
  ≥ 4,522 genes in the prior study, and 1,288 of the primary closure's proteins were never searched.

**L1 hubs** (T-hubs). Degree and body length inside the member-holding L1 groups:

| L1 form | group | internal L1 edges | E1 edges | edges with ≥ 1 shared exonic bp | exon-proxy-only edges (with 0 shared exonic bp) | top-degree nodes |
|---|---|---|---|---|---|---|
| primary | NPIP (341) | 1,285 | 860 | 988 | 130 (0) | VN1R91P (95; 889 bp, span exon) · BNIP3P16 (63; 2,051 bp, span exon) · NPIPB14P (40) · PDXDC2P-NPIPB14P (35; 90,079 bp) · PKD1P4-NPIPA8 (35) |
| primary | TBC1D3 (97) | 374 | 288 | 319 | 21 (0) | LGALS9B (41; 249,716 bp) · USP6 (26) · TBC1D3, B, D, E, F, G (20 each) |
| 17:12, no v-exon/strand | NPIP (604) | 2,362 | 1,239 | 1,406 | 487 (274) | LOC101926889 (104) · VN1R91P (95) · LOC105372586 (69) · BNIP3P16 (64) · PDXDC1 (52) |
| 17:12, no v-exon/strand | TBC1D3 (897) | 4,024 | 2,430 | 2,711 | 498 (250) | TM2D3 (108) · MIR4713HG (60) · ARHGAP11B-DT (59) · GOLGA8T (58) · GOLGA8J (57) |

**Where the hubs come from.** In the 17:12 report they were blamed on clause 2 itself, but the two kinds have
different causes.

- **Intron-landing hubs came from the missing v-exon requirement** (hub attribution table in T-hubs):
  - TM2D3 had 108 edges; 107 of them passed by the exon proxy alone, 100 of those with 0 shared exonic bp. It keeps 6.
  - MIR4713HG: 60, 60, 57; keeps 2.
  - LOC101926889: 104, 104, 100; keeps 3.
  - LOC105372586: 69, 69, 63; keeps 6.
- **Hubs that survive come through the shorter-body gene-body disjunct:**
  - A short span-exon pseudogene aligned inside other genes: "≥ 0.50 of the shorter body" admits VN1R91P (95 edges,
    all gene-body) and BNIP3P16 (63).
  - A very long record containing copies: TBC1D3 → LGALS9B, gene-body 0.850 at chain identity 0.907, shared exon 0.
  - The strand check also removes 32 of PDXDC1's 48 v-exon-overlapping edges inside its 17:12 L1 group.

**Neighbours of interest** (T-chaining). "with" = in the anchor's group (NPIPB2 or TBC1D3); apart(n) = in its own group of
n genes.

| gene | L0 | L1 | L2 | L3 | L1Δ | L2Δ | L3Δ | L2 / L3 without readthrough nodes |
|---|---|---|---|---|---|---|---|---|
| PKD1, PKD1P1, P2, P3, P6 | with | with | with | with | with | with | with | apart(9) / apart(9) |
| PKD1P3-NPIPA1, PKD1P4-NPIPA8, PKD1P5-LOC105376752, LOC131696449 (readthroughs) | with | with | with | with | with | with | with | removed |
| PDXDC2P-NPIPB14P (readthrough) | with | with | with | apart(3) | with | with | apart(1) | removed |
| NPIPB14P (member) | with | with | with | apart(2) | with | with | apart(1) | with / apart(2) |
| DHX40, DHX40P1 | with | with | **apart(3)** | apart(3) | apart(1) | apart(1) | apart(1) | apart(2) / apart(2) |
| TBC1D3P1-DHX40P1 (readthrough) | with | with | **apart(3)** | apart(3) | with | apart(1) | apart(1) | removed |
| RNFT1-DT | apart(1) | apart(1) | apart(1) | apart(1) | apart(1) | apart(1) | apart(1) | apart(1) / apart(1) |
| RNFT1, RNFT1P3 | apart(2) | apart(2) | apart(2) | apart(1) | apart(1) | apart(1) | apart(1) | apart(2) / apart(1) |
| TBC1D26, TBC1D28 | with | **apart(14)** | apart(5) | apart(2) | apart(1) | apart(1) | apart(1) | apart(5) / apart(2) |
| USP6, TBC1D29P, LOC100420408 | with | with | with | apart(1) | with | with | apart(1) | with / apart(1) |
| USP6NL | with | apart(1) | apart(1) | apart(1) | apart(1) | apart(1) | apart(1) | apart(1) / apart(1) |

**How PKD1 reaches NPIP** (T-chaining paths, T-L2paths). NPIPB2 → LOC131696449 (readthrough PKD1P1-NPIPA5L; shared exon
0.548, w_98 0.980) → PKD1 (shared exon 0.398, w_98 0.982).

- The same two edges carry the link at L1, L2 and L3. The pooled identity of the second edge is 0.972, so under the 17:12
  L3 the route ran through LOC100288162 instead.
- **Without readthrough nodes,** PKD1, PKD1P1/P2/P3/P6 and 4 lncRNA LOC records form their own 9-gene group at L2 and L3.

**Which "stays / leaves" verdicts a truth supports** (Soto flag-ok labels, T-soto-neighbours, T-truth-ingroup):

- **PKD1 joining NPIPA1/A9 is not a false merge by Soto.** Soto labels PKD1, PKD1P1, P2, P3, P6 and the readthroughs
  LOC131696449, PKD1P3-NPIPA1 and PKD1P5-LOC105376752 with NPIPA1, NPIPA9 and the member PKD1P6-NPIPP1 as ID_149 (NPIPA6
  is a weak match). Most other NPIP members are ID_154.
- **The Soto-false content of the NPIP group:**
  - at L2 and L3, ID_149 merged with ID_154 (16 labelled genes, including the anchor NPIPB2), and at L3 also the
    single-gene ID_152 and ID_153;
  - at L3, ID_41 (SMG1P1, P2, P4-P7 and 4 LOC records) and ID_302 (BOLA2, BOLA2B, LOC124907841);
  - at L2 also ID_28 (PLA2G10EP-KP and 4 LOC records, including LOC100505915).
  - The ZNF-named and BNIP3P records have no ok Soto label, so no truth used here scores them.
- **TBC1D3 side:** DHX40/DHX40P1 are Soto ID_341, separate from TBC1D3's ID_468, so their L2 split agrees with Soto.
  TBC1D3P1-DHX40P1 carries both labels (ambiguous). USP32P1-3 (ID_60) and USP32/USP32P4 (ID_91) stay at L2, against
  Soto.

**How the ZNF genes enter L2** (T-L2paths): NPIPB2 → NPIPB5 → LOC128966632 (SMG1-like) → SMG1P7 → **VN1R91P** (single
span exon; shared exon 0.619 and 0.711) → ZNF737 → ZNF728 → ZNF429.

**Diagnostic, not a level:** the shared-exon clause is made unsatisfiable on edges that touch pseudogene records whose
exon union equals the whole body.

- The NPIP L2 group shrinks to 53 genes, with 0 ZNF.
- It also drops 1 NPIP and 5 TBC1D3 members: the exon-less member pseudogenes. The TBC1D3 group becomes 18 genes with
  11 members.
- PKD1 stays.

## 6. Truth agreement vs the prior study's report-only groupings (`LAT/truth.tsv`, `LAT/truth_ingroup.tsv`; T-truth-U, T-truth-ingroup, T-truth-V)

**Circularity, stated first.**

- **Soto.** Soto families are SD98 (≥ 98% identity) duplications, kept only where the map-back covers the shared exons.
  Those are the conventions L2 (shared exon ≥ 0.30) and L3 (identity ≥ 0.98) test. Soto agreement at L2 and L3 is
  therefore partly by construction.
- **Clause 5 (C_tree).** It was developed on NPIP against the literature groups (§6jp–§6js), so it is not an independent
  comparator on these families.
- **C_L1 / C_fine** are the literature groups themselves.

**Metrics:**

- Pairwise precision and recall.
- Bipartite F with one-to-one Jaccard matching (`lattice_common.bip_jaccard`). F is NA when either side has 0 pairs.
- **In-group pair precision:** all labelled genes of the anchor's group, members and pulled-in non-members.

**Scorer check:** the scorer reproduces 12 of 12 rows of `integrate_slim/truth_agreement.tsv` (P and D as built; Soto and
HGNC) (`LAT/truth.out`).

**Gene sets:**

- **U:** each report-only grouping's own universe inside U, with the lattice scored on the same genes. These scores see
  members only; read the in-group table next to them.
- **V:** recall is conditioned on the lattice's own closure.

**Bipartite F on the same genes (U):**

| truth | family | genes | L0 | L1 | L2 | L3 | L1Δ | L2Δ | L3Δ | L3 w_98 gap-incl | L3 pooled gap-excl | L3 S2 | 17:12 L3 | report-only grouping on these genes: F (pair P / R) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Soto (flag ok) | TBC1D3 | 16 (D universe) | 0.688 | 0.750 | 0.875 | **1.000** | 0.839 | 0.839 | 0.933 | 0.933 | 1.000 | 0.897 | 1.000 | D/E1 MCL 0.750 (0.543 / 1.000) |
| Soto | TBC1D3 | 11 (P universe) | 0.818 | 0.909 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | P MCL 0.909 (0.800 / 1.000) |
| Soto | NPIP | 31 (D universe) | 0.516 | 0.516 | 0.516 | 0.581 | 0.516 | 0.548 | 0.581 | 0.712 | 0.581 | 0.615 | 0.516 | D/E1 MCL 0.806 (0.588 / 0.897) |
| Soto | NPIP | 19 (P universe) | 0.789 | 0.789 | 0.789 | 0.789 | 0.789 | 0.789 | 0.789 | 0.600 | 0.789 | 0.480 | 0.789 | P MCL 0.842 (0.693 / 0.981) |
| Soto | pooled | 47 (D universe) | 0.340 | 0.596 | 0.638 | 0.723 | 0.624 | 0.645 | 0.696 | 0.787 | 0.723 | 0.716 | 0.681 | D/E1 MCL 0.787 (0.575 / 0.922) |
| HGNC (superfamily-level) | TBC1D3 | 13 (D universe) | 0.923 | 0.880 | 0.960 | 0.917 | 0.960 | 0.960 | 0.917 | 0.917 | 0.917 | 0.917 | 0.833 | D/E1 MCL 0.800 (0.682 / 0.818) |
| HGNC | TBC1D3 | 12 (P universe) | 0.917 | 0.818 | 0.909 | 0.909 | 0.909 | 0.909 | 0.909 | 0.909 | 0.909 | 0.909 | 0.909 | P MCL 1.000 (1.000 / 1.000) |
| literature NPIPA\|NPIPB | NPIP | 21 (C_tree leaves) | 0.634 | 0.634 | 0.634 | 0.634 | 0.634 | 0.634 | 0.634 | 0.765 | 0.634 | 0.727 | 0.634 | C_tree_top 0.765 (1.000 / 0.411) |
| literature NPIPA\|NPIPB | NPIP | 22 (lit leaves) | 0.651 | 0.651 | 0.651 | 0.619 | 0.651 | 0.651 | 0.619 | 0.743 | 0.619 | 0.706 | 0.651 | C_L1 1.000 (circular) |
| literature L2 | NPIP | 21 (C_tree leaves) | 0.238 | 0.238 | 0.238 | 0.238 | 0.238 | 0.238 | 0.238 | 0.667 | 0.238 | 0.762 | 0.238 | C_tree_min 0.833 (1.000 / 0.250); C_tree_top 0.667 |
| literature L2 | NPIP | 22 (lit leaves) | 0.227 | 0.227 | 0.227 | 0.273 | 0.227 | 0.227 | 0.273 | 0.682 | 0.273 | 0.773 | 0.227 | C_fine 1.000 (circular) |
| literature L2 | TBC1D3 | 9 | 0.222 | 0.222 | 0.222 | 0.222 | 0.222 | 0.222 | 0.222 | 0.222 | 0.222 | 0.222 | 0.222 | C_tree_min 0.667; C_tree_top 0.444 [corrected 2026-09-16: **C_tree_min 0.778; C_tree_top 0.556** under Fig 6C — see `bench/TBC1D3_GUITART_TRUTH_CORRECTION.md`] |

**Pair precision / recall on all of U** (T-truth-U-pairs), L0 → L3:

- **Soto TBC1D3 (16):** 0.475 / 1.000 → 0.543 / 1.000 → 0.722 / 1.000 → 1.000 / 1.000.
- **HGNC TBC1D3 (14):** 0.846 / 1.000 → 0.818 / 0.682 → 1.000 / 0.682 → 1.000 / 0.545.
- **Soto NPIP (32):** 0.376 / 1.000 at L0–L2; 0.431 / 1.000 at L3.

**In-group pair precision of the anchor's whole group** (T-truth-ingroup):

| level | NPIP group: size (members); Soto-labelled genes; same-family pairs / all = precision; Soto families | TBC1D3 group: same columns |
|---|---|---|
| L0 | 5,764 (42); 474 labelled; 1,645 / 112,101 = 0.015; 170 families | the same group |
| L1 | 341 (26); 58; 279 / 1,653 = 0.169; 11 | 97 (16); 51; 100 / 1,275 = 0.078; 15 |
| L2 | 142 (26); 55; 278 / 1,485 = 0.187; 9 | 24 (16); 18; **60 / 153 = 0.392**; 4 (ID_468 11, ID_60 3, ID_469 2, ID_91 2) |
| L3 | 55 (25); 42; 223 / 861 = 0.259; 6 (ID_154 16, ID_149 11, ID_41 10, ID_302 3, ID_153 1, ID_152 1) | 11 (11); 11; 55 / 55 = 1.000; 1 (ID_468) |
| L3Δ | 37 (25); 29; 175 / 406 = 0.431; 4 | 11 (11); 1.000 |
| 17:12 L2 / L3 | 159 (26): 0.187 / 73 (26): 0.194 | 31 (16): 0.392 / 16 (11): 1.000 |

**On V** (T-truth-V; pooled; recall conditioned on V). Bipartite F and pair precision, L0 → L3:

| truth | genes | F, L0 → L3 | pair P, L0 → L3 | L1Δ | L1 = E1 at 0.80 / 0.50 | L3 w_98 gap-incl | L3 S2 | 17:12 L3 |
|---|---|---|---|---|---|---|---|---|
| Soto | 550 | 0.193 → 0.486 → 0.677 → 0.737 | 0.015 → 0.115 → 0.244 → 0.440 | 0.551 | 0.571 | 0.734 | 0.616 | 0.691 |
| HGNC | 5,666 | 0.137 → 0.435 → 0.455 → 0.440 | 0.009 → 0.165 → 0.555 → 0.822 | 0.442 | 0.443 | 0.438 | 0.432 | 0.439 |

**Reading:**

- **TBC1D3.** Soto F on the member genes rises from 0.688 to 1.000, and the L3 group equals Soto ID_468 (11 genes).
  - Soto shares the SD98 and shared-exon conventions, so L2/L3 agreement is partly by construction.
  - The whole L1 and L2 groups have Soto in-group precision 0.078 and 0.392. L2 is not "enough" for TBC1D3: it keeps
    the USP32 family (ID_60, ID_91).
- **NPIP.** One group holds all connected NPIP members down to L2, and all but NPIPB14P at L3.
  - Soto F moves only at L3 (0.531 → 0.594), and in-group precision stays low (0.187, 0.259).
  - Part of that is Soto splitting NPIP itself (ID_154 vs ID_149).
- **Against the report-only MCL groupings:** on NPIP they do better than the lattice's components (Soto: D 0.806 vs 0.581
  on 31 genes; P 0.842 vs 0.789 on 19 genes).
- **Against the literature, the result is specific to the identity choice.**
  - With the primary gap-excluded w_98, C_tree does better: NPIPA|NPIPB 0.765 vs 0.634; literature L2 C_tree_min 0.833
    vs 0.238.
  - With gap-inclusive w_98, L3 equals C_tree_top on NPIPA|NPIPB (0.765) and scores 0.667 on literature L2. That variant
    was seen after the primary was fixed.
  - Clause 5 was developed on NPIP, so neither comparison is independent.

## 7. Threshold filtration: single linkage on fixed evidence

**What this is, and what it is not.** One field is swept over a frozen evidence table. Sweeping a threshold gives nested
partitions by construction: a single-linkage dendrogram. It shows how the groups at one level refine as a cut rises. It
is **not** the "more information" picture: adding copies or evidence is T3(b), which is not established here (§4).

**Construction.** Inside the member-holding primary L1 groups, keep edges whose field w is ≥ t, for t in {0.70, 0.80, 0.90,
0.95, 0.98, 0.99, 1.00}. An edge with a missing field passes no t.

- **Shared-exon OFF:** L1 edges, with w = the pooled gap-excluded identity. No single-record field is defined without the
  shared-exon witness. This sweep is not one of the §0★★★.3 level sweeps, and its field is not a maximum over records.
- **Shared-exon ON:** L2 edges, with w = w_98 gap-excluded, the primary L3 field. At t = 0.98 the result is exactly L3.
- Full listing with every non-member, including ON with the pooled identity: `LAT/filtration.txt`. Group table:
  `LAT/filtration_groups.tsv`.
- Notation: `+n` = non-members in the same group; A/B = NPIPA/NPIPB; T3 = TBC1D3.

**NPIP (L1 group: 341 genes, 26 members), OFF (pooled) → ON (w_98)**

```
L1 / L2          {all 26 members +315}                         | {all 26 +116}
w >= 0.70-0.80   {all 26 +313}                                 | {all 26 +116}
w >= 0.90        {all 26 +218}                                 | {all 26 +66}
w >= 0.95        {all 26 +86}                                  | {all 26 +61}
w >= 0.98        {all 26 +51}                                  | {all but B14P +30} {B14P +1}                      (= L3)
w >= 0.99        {B6 B7 B8 B9 B10P B15 LOC124907807/808}       | {B3 B4 B5 B11 B12 B13 LOC124907834 LOC128966608 +6}
                 {B3 B4 B5 B12 B13 LOC124907834 LOC128966608   | {B6 B7 B8 B9 B10P B15 LOC124907807/808}
                  +8}                                          | {A1 A5 A6 A7 A8 A9 PKD1P6-NPIPP1 +7}
                 {A1 A5 A6 A7 A8 A9 +7}                        | {A2} {B2} {B14P}
                 {A2} {B2} {B11} {B14P} {PKD1P6-NPIPP1}        |
w >= 1.00        26 member groups; {B5 +1} {B13 +1}            | 26 member groups; {B5 +1}
```

**TBC1D3 (L1 group: 97 genes, 16 members), OFF (pooled) → ON (w_98)**

```
L1 / L2          {all 16 +81}                                  | {all 16 +8}
w >= 0.70-0.80   {all 16 +81}                                  | {all 16 +8}
w >= 0.90        {all but P7 +76} {P7}                         | {14 members +6} {P5 +2} {P7}
w >= 0.95        {T3 B D E F G H I K P1 P2 +3} {P3 P4 +60}     | {T3 B D E F G H I K P1 P2} {P3 P4} {P5 +1}
                 {P5 +1} {LOC100420311} {P7}                   |  {LOC100420311} {P7}
w >= 0.98        {T3 B D E F G H I K P1 P2 +3} {P3 P4 +29}     | {T3 B D E F G H I K P1 P2} {P3 P4}                (= L3)
                 {LOC100420311} {P5} {P7}                      |  {LOC100420311} {P5} {P7}
w >= 0.99        {9 coding} {P1 P2 +3} {P3 P4 +14}             | {9 coding} {P1 P2} {P3 P4} {LOC100420311} {P5} {P7}
                 {LOC100420311} {P5} {P7}                      |
w >= 1.00        16 singletons                                 | 16 singletons
```

**Appearance thresholds** (T-appearance). A reference group G appears at t when two conditions hold:

- **Whole:** t ≤ the minimum bottleneck inside G.
- **Separated:** t > the maximum bottleneck from G to the rest of its reference set.

Reference sets: NPIP = the 21 literature records that are catalog nodes (NPIPB1P is outside both catalogs);
TBC1D3 = the 9 clause-5 leaves. Intervals are printed at 6 decimals from unrounded fields.

| group | pooled gap-excl, OFF | pooled gap-excl, ON | pooled gap-incl, OFF | **w_98 gap-excl, ON (L3 field)** | w_98 gap-incl, ON |
|---|---|---|---|---|---|
| NPIPA \| NPIPB (no A–B pair together) | t > 0.981757 (grid 0.99) | t > 0.980559 (0.99) | t > 0.974041 (0.98) | **t > 0.980559 (0.99)** | t > 0.971665 (0.98) |
| NPIPA exact | (0.981757, 0.989326] | (0.980559, 0.989326] | (0.974041, 0.978889] | **(0.980559, 0.989897]** | (0.971665, 0.981585], grid 0.98 |
| NPIPB exact | never | never | never | **never** | never |
| A6-9 | never | never | (0.988425, 0.988689] | **never** | (0.988425, 0.988689] |
| B3-5 | never | never | never | **(0.997341, 0.998666]** | never |
| B6-9 | (0.994729, 0.996785] | (0.994729, 0.996785] | never | **(0.994729, 0.996797]** | never |
| B12/13 | (0.998402, 0.998932] | (0.996872, 0.996899] | never | **(0.997341, 0.999014]** | (0.988086, 0.990873], grid 0.99 |
| named NPIPB {B3,B4,B5,B11,B12,B13} | (0.987108, 0.987703] | (0.987108, 0.987703] | never | **(0.987810, 0.992067], grid 0.99** | never |
| TBC1D3 {B,F,G,H} (clause 5) | never | never | never | **never** | never |
| TBC1D3 {TBC1D3,D,E,K} (clause 5) | (0.997303, 0.997982] | (0.997303, 0.997982] | (0.996591, 0.997159] | **(0.997303, 0.997982]** | (0.996591, 0.997159] |
| {B,F,G,H} \| {TBC1D3,D,E,K} | t > 0.997303 (1.00) | t > 0.997303 (1.00) | t > 0.996591 (1.00) | **t > 0.997303 (1.00)** | t > 0.996591 (1.00) |

The table is restricted to the primary L1 groups. On the 17:12 L1 groups with unrounded identities
(`LAT/filtration_appearance.17_03_tests_exact.tsv`, T-appearance.17_03_tests_exact), the pooled gap-excluded ON window
for B12/13 is (0.996872, 0.996905]. The 17:12 table said "never" because it rounded to 4 decimals.

**Why most groups never appear** (merge heights, `LAT/filtration.txt`, end section).

- **NPIPB:** on the L3 field, NPIPB14P joins the rest of NPIP at 0.978354, below the A|B split at 0.980559. The
  whole-NPIPB condition therefore fails before the separation condition holds.
- **A6-9:** A1 joins {A6, A9} at 0.995679, above the height at which {A7, A8} joins them (0.994447).
- **{B,F,G,H}:** {B,H} joins the group already holding {TBC1D3,D,E,F,G,I,K} at 0.996053. F (0.996240) and G (0.997303)
  merge with the {TBC1D3,D,E,K} side first.
- **The reason:** single linkage merges through whichever member pair is closest, and paths may run through non-member
  records.

## 8. Triangle-support variant (3-truss; T-sizes, T-triangle)

**Construction.** An edge is kept iff it lies in ≥ 1 triangle of the same level graph. The fixed point is reached after
one pass: deleting an edge that is in no triangle cannot remove any other edge's triangle (peel depth 1 at every level,
`LAT/levels.out`). Edges kept per level: 223,059, 9,716, 7,680 and 5,602 (dropped 1,048, 560, 354 and 152).

| level | primary: components ≥ 2 genes / 2-gene components / singletons | 3-truss: components ≥ 2 / 2-gene / singletons | NPIP group | TBC1D3 group | member-holding groups dissolved |
|---|---|---|---|---|---|
| L0 | 115 / 55 / 1,587 | 95 / 0 / 2,359 | 5,764 → 4,657 (one group for both) | — | none |
| L1 | 192 / 87 / 6,242 | 85 / 0 / 6,824 | 341 → **230 (33% removed)** | 97 → **63 (35%)** | none |
| L2 | 270 / 120 / 6,594 | 129 / 0 / 7,058 | 142 → 64 (55%) | 24 → 23 (USP32P4 drops) | none |
| L3 | 153 / 76 / 7,353 | 63 / 0 / 7,591 | 55 → 37 (33%) | 11 → 11 | **{TBC1D3P3, TBC1D3P4}: both Soto ID_469, a true 2-copy family**; {NPIPB14P, LOC124907811} |

Under the 17:12 tests, unrounded: L1 604 → 340 (44% removed) and 897 → 78 (91%); L3 73 → 67 and 16 → 13 (T-sizes).

**What the 3-truss changes:**

- **Drops, by level:**
  - L2: 78 NPIP genes, including 39 of the 40 ZNF-named records.
  - L3 (NPIP): 18 non-members: BOLA2, BOLA2B, BOLA2-SMG1P6, SMG1P1-7, SLC7A5P2 and 7 LOC records.
- **Separation from the anchor:**
  - DHX40, DHX40P1, TBC1D26 and TBC1D28 are apart at L1Δ.
  - **PKD1, PKD1P1/P2/P3/P6 and the four PKD1P readthroughs stay with NPIP at every triangle level.**
- **Truth scores** (§6):
  - Soto TBC1D3: L1Δ 0.839, L2Δ 0.839, L3Δ 0.933, vs L3 1.000. The loss is the dissolved {TBC1D3P3, TBC1D3P4}.
  - HGNC TBC1D3 on D's genes: L1Δ 0.960, vs L1 0.880.
- **The 2-copy cost cannot be measured on these families.** They are large and dense, and the one member-holding 2-copy
  family they contain is dissolved. §0★★★.5 calls 2-copy families the modal family size.

**Expression view.** Two views exist and both nest (Δ(G_k[X]) ⊆ Δ(G_k)[X]); T-expr gives both.

- **The view reported here is comp(Δ(G_k[X])):** the truss of the induced expressed subgraph.
  - At any-overlap reads ≥ 3, the member component has 16, 16 and 16 genes (12 members) at L0–L2.
  - **At L3 it is two groups:** {NPIPA1, A2, A7, A9, PKD1P6-NPIPP1} + PKD1, PKD1P3-NPIPA1, PKD1P5-LOC105376752
    (8 genes), and {NPIPB2, B4, B5, B6, B7, B9, B11} (7 genes).
- **The other view, comp(Δ(G_k)[X]):** 16, 16, 16 and 15 genes, one group with 12 members at every level.
- **At reads ≥ 1** (NPIP / TBC1D3 components):
  - Δ(G_k[X]): L0 46 (23 members) / 2,500 (4); L1 46 / 16; L2 30 / 9; L3 15 (15) + 10 (7) / 4 (4).
  - Δ(G_k)[X]: L0 48 / 2,715; L1 48 / 19; L2 32 / 9; L3 25 (22) / 4.
- **Correction of the 17:12 §8.** Under the 17:12 tests with unrounded identities, Δ(G_3[X]) at reads ≥ 3 is also two
  groups: {NPIPB2, B4, B5, B6, B7, B9, B11} + LOC124907832 (8 genes) and {NPIPA1, A2, A7, A9, PKD1P6-NPIPP1} + 2
  (7 genes). The 17:12 report's single group of 15 came from two sub-0.98 edges that rounding admitted: NPIPA2–NPIPB2
  (0.9799886) and LOC112268174–NPIPB9 (0.9799584).

## 9. Expression views (testis, `human_testis.t2t.bam`; T-expr)

**Counts.** `npip_tbc1d3.py lattice-expr` counts all 8,070 V genes with the `expr-recount` rule (one function,
`lattice_common.count_reads`), using primary reads (`-F 2308`).

- **Any-overlap:** a read counts for a gene if one of its blocks overlaps an exon of the gene, at any MAPQ.
- **Unique:** the read's blocks hit exons of exactly one RefSeq record genome-wide.
- **Conventions.** Any-overlap (strand ignored, ≥ 1 bp of an exon) is the §0★★★.1 any-overlap convention, the lattice's
  choice. Unique is neither §0★★★.1 convention; in particular it is not u ≥ 3 (MAPQ ≥ 1, ≥ 1 aligned block inside the
  node interval, the convention of §6jz–§6kf). Any-overlap counts are not copy-resolved: TBC1D3P1-DHX40P1 has 72
  any-overlap reads and 0 unique reads (`LAT/expr_counts.tsv`).
- **Check:** the 241 genes also in `integrate_slim/expr_recount.tsv` are identical, and `expr_counts.tsv` is
  byte-identical to the 17:04 run (`LAT/expr_counts.out`).

| expression set | expressed V genes | expressed members (NPIP / TBC1D3) | L0 | L1 | L2 | L3 |
|---|---|---|---|---|---|---|
| any-overlap reads ≥ 3 | 3,063 | 13 / 0 | 16 genes (12 members) | 16 (12) | 16 (12) | 15 (12) |
| any-overlap reads ≥ 1 | 4,776 | 24 / 4 | **3,022 (4 TBC1D3) and 48 (23 NPIP), separate** | 48 (23) and 27 (4) | 42 (23) and 9 (4) | 27 (22) and 4 (TBC1D3, E, G, P1) |
| unique reads ≥ 3 | 2,207 | 9 / 0 | 8 (8) | 8 (8) | 8 (8) | 6 (6) and 2 (2) |

**Members in the expressed components:**

- **Any-overlap ≥ 3:** NPIPA1, A2, A7, A9, NPIPB2, B4, B5, B6, B7, B9, B11 and PKD1P6-NPIPP1.
  - At L0 the component also holds PKD1, PKD1P3-NPIPA1, PKD1P5-LOC105376752 and PKD1P6.
  - NPIPB1P (6 reads) is expressed but has no edges.
- **Unique ≥ 3:** NPIPA2, A7, NPIPB2, B4, B5, B6, B9 and B11, with no non-members. L3 separates {NPIPA2, A7} from the six
  NPIPB records.

**Clause-6 class (c), on the observed graph.** At reads ≥ 1 the L0 view splits NPIP from TBC1D3, whereas L0 itself holds
them in one group. By T2b the split can only refine the L0 group.

- **On the observed graph,** every L0 path between them passes through a copy with 0 any-overlap reads: an unexpressed
  bridge.
- **The observed graph is truncated.** 1,288 proteins of the primary closure were never searched, and no DNA evidence
  crosses chromosome trios. An expressed bridge through unobserved edges is not ruled out.

## 10. Known limits

- **Scope.**
  - Two development families: descriptive, not pre-registered.
  - One haplotype (CHM13) and one annotation (RefSeq).
  - One tissue (testis); 0 TBC1D3 members are expressed at reads ≥ 3.
- **Clause 2 is approximate** (§1, §2):
  - The exon clause is a proxy on gene-body alignments. Record identity includes introns, and the exon union stands in
    for the transcript.
  - The gene-body clause uses the shorter-body denominator on gene-body targets.
  - minimap2 `-X -N 50 -p 0.1` limits which pairs exist, and those caps depend on the target set (H2).
  - The v-exon and strand requirements are checked against the shipped loop only for the gene-body disjunct. No
    spliced-transcript alignment exists to check the exon disjunct.
- **Monotonicity (H3)** fails for the shipped t_P, for the gene-body chains and for the pooled-identity variants, and
  **pairwise evidence (H2)** fails for blastp E-values and the minimap2 caps (§4). T3(b) is established neither for more
  copies (needs H1, H2) nor for more homology (also needs H3).
- **The E1 catalogs are per chromosome trio:**
  - No chr16↔chr17 DNA edge exists, so L1–L3 cannot join NPIP and TBC1D3.
  - Members on chr1/4/18 (NPIPB1P, TBC1D3P6, LOC124905656, LOC100420289) are DNA singletons.
  - The two catalogs cover 12,048 genes.
- **Protein evidence is genome-wide but the closure is truncated:**
  - 4,430 of 20,088 proteins were searched; 1,288 of the primary closure's proteins never were (e.g. DHX40).
  - E-values were not saved.
  - V is not closed under the union-cover protein test (33 genes outside V).
- **Node conventions drive chaining:**
  - Readthrough and nested records are nodes (same-locus links excluded; variants reported).
  - Exon-less records count as one span exon, which inflates the shared-exon fraction for pseudogenes (the ZNF route).
  - f_ex and w_98 are taken from single records.
- **L3 identity:**
  - The gap convention decides whether L3 splits NPIP: gap-excluded w_98 keeps every connected member except NPIPB14P in
    one group; gap-inclusive w_98 splits it into three groups.
  - The better literature scores of the gap-inclusive variants were seen after the primary was fixed.
- **Truths:**
  - HGNC 2227 is superfamily-level, and HGNC gives NPIP no group.
  - Soto needs the exon-overlap mapping (flag ok only), and it shares the SD98 and shared-exon conventions with L2/L3.
  - Literature C is circular for C_L1 / C_fine, and clause 5 (C_tree) was developed on NPIP.
  - TBC1D3's literature L1 (cluster1/cluster2) is positional; its literature L2 has only two non-trivial groups, AE
    (TBC1D3, E) and CDKL (D, K). [corrected 2026-09-16: its literature L2 has two non-trivial groups read post hoc
    from Guitart Fig 6B/6C, M (B, H) and CDKL (K, TBC1D3); TBC1D3D's group is unresolved. See
    `bench/TBC1D3_GUITART_TRUTH_CORRECTION.md`.]
  - V-scope recall is conditioned on the lattice's own closure.
- **Not done:**
  - Leaders as a report-only grouping.
  - A split-support L3.
  - A clause-2 rebuild with real spliced-transcript alignments.
  - A T3 test with added annotation or evidence, using the monotone test forms.

## 11. Reproduce (foreground, in order; `LAT` = the results directory)

Wave 7 (2026-09-24) replaced the 8 `lattice_*.py` scripts with one CLI, `bench/layer_order/npip_tbc1d3.py`, with one
subcommand per old script. The old scripts are at git tag `notebook-2026-09-24`. Every stage **overwrites** its outputs
under the root. The default root is the frozen `/mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3`, so to rerun,
copy `light/ heavy/ integrate_slim/ lattice/` and pass `--root COPY` (or set `LO_ROOT`). `npip_tbc1d3.py all
--with-check-c2 --root COPY` runs both reproduce blocks: this one and `LAYER_ORDER_NPIP_TBC1D3.md` §11.

```
python3 bench/layer_order/npip_tbc1d3.py --root COPY lattice-edges      # LAT/{edges,nodes,closure,edges_all_c2}.tsv, edges_build.out (64 s, 2.5 GB)   (was lattice_edges.py)
python3 bench/layer_order/npip_tbc1d3.py --root COPY lattice-expr       # LAT/expr_counts.tsv, expr_counts.out                     (43 s)   (was lattice_expr.py)
python3 bench/layer_order/npip_tbc1d3.py --root COPY lattice-levels     # LAT/{levels,member_groups,groups,sanity,expr_views,chaining,triangle_drops}.tsv, levels.out (44 s)   (was lattice_levels.py)
python3 bench/layer_order/npip_tbc1d3.py --root COPY lattice-truth      # LAT/truth.tsv, truth_ingroup.tsv, truth.out               (8 s)   (was lattice_truth.py)
python3 bench/layer_order/npip_tbc1d3.py --root COPY lattice-filtration # LAT/filtration.txt, filtration_groups.tsv, filtration_appearance.tsv (5 s)   (was lattice_filtration.py)
python3 bench/layer_order/npip_tbc1d3.py --root COPY lattice-filtration --l1 c2_loose   # the same with suffix .17_03_tests_exact (17:12 L1 groups) (5 s)
python3 bench/layer_order/npip_tbc1d3.py --root COPY lattice-check-c2   # LAT/check_c2.out (shipped chains and shipped v-exon/strand loop) (31 s)   (was lattice_check_c2.py)
python3 bench/layer_order/npip_tbc1d3.py --root COPY lattice-report     # LAT/report_tables.md (every table in this report)          (26 s)   (was lattice_report_tables.py)
```

Rerun check (wave 7, on a copy of the frozen directory). The old scripts, with only the repo-root path fixed, and the
new CLI give byte-identical files for all 61 outputs of both reproduce blocks. The exceptions are the timing tokens of
`edges_build.out`, `expr_counts.out` and `levels.out`. `check_c2.out` also equals the frozen 2026-09-16 file. Four
caveats:
- `lattice_check_c2.py` had not run since 2026-09-19. It imported the archived `denovo_shared_def`, and it `exec`'d the
  head of `lattice_edges.py`, which has needed `__file__` since §6r9. The baseline needed both patched.
- Seven outputs depend on Python's string-hash seed: ties are broken by set iteration order. They are
  `LAT/{member_groups,groups,triangle_drops}.tsv`, `LAT/report_tables.md`, `IS/expr_groups.tsv`, `IS/expr_sweep.tsv`
  and `IS/analysis.out`.
- The frozen files hold one random seed's tie order: the same rows and numbers, some in another order. The CLI pins
  `PYTHONHASHSEED=0`.
- `LAT/corrections_pass2/{doc_numbers,strand_diag}.py` import the new `lattice_common` by path. Rerun on the new
  outputs, they reproduce their frozen `.out` files byte for byte.

Shared definitions (paths, the four tests and their variants, union-find, 3-truss, split counts) live in
`bench/layer_order/lattice_common.py`. The rebuild check against the 17:03 table is
`LAT/corrections_pass2/cmp_build.py`.

## 12. Correction log and audit notes

**Every audit finding was applied; none was judged wrong.** Findings are listed in audit order, with the section where
each correction lives.

- **(1) L3 thresholded 4-dp-rounded identities (important).**
  - `edges.tsv` now stores floats at full precision, and every table is formatted from unrounded values.
  - Under the 17:12 tests, unrounded: 5,868 L3 edges (was 5,870) and 5,711 3-truss edges (was 5,713). The two admitted
    edges were NPIPA2–NPIPB2 (0.9799886) and LOC112268174–NPIPB9 (0.9799584).
  - Δ(G_3[X]) at reads ≥ 3 is two groups (§8), and the B12/13 ON window is (0.996872, 0.996905] (§7). All reproduce
    the audit (T-edges, T-expr, T-appearance.17_03_tests_exact).
- **(2) The exon proxy ignored target exons (important).** The variant row "v-exon overlap on the exon proxy only" gives
  L0 6,162 and L1 474 / 725, as the audit found. It is superseded as primary by (4).
- **(3) Pooled L3 identity is neither w_98 nor monotone; T3 coverage was overstated (important).**
  - w_98 is computed as specified. The gap-excluded form is now primary; the gap-inclusive form is a variant (§1, §3).
  - The non-monotonicity of pooled identity, greedy t_P and greedy chains, and the failure of H2, are stated (§1, §4).
  - The "C_tree beats every level" sentence is now tied to the identity choice (§6).
  - Under the 17:12 L1, w_98 gives 5,905 L3 edges with the audit's member groups (T-levels).
- **(4) L1 omitted "the hit must overlap v's exons" (important).**
  - The requirement is now primary, together with the shipped strand check.
  - The strand part is validated against the shipped `cmd_families` loop on 78,423 of 78,423 pairs.
  - The 17:12 L1 (604 / 897) is kept as a variant. Hubs are re-attributed (§5), and the definition doc's §0★★★.5 numbers
    are updated.
  - **Extension, not a disagreement.** The audit's corrected numbers (L0 6,050; L1 378 / 111; L2 and L3 "unchanged"
    159 / 31 and 73 / 16) were computed without the strand check, and they reproduce exactly in the "v-exon overlap, no
    strand check" variant. With the strand check, which the fix also asked for, L2 and L3 do change: L0 5,764; L1 341 /
    97; L2 142 / 24; L3 NPIP 55 (NPIPB14P leaves) and TBC1D3 11.
  - The strand check removes 995 L1 edges (247 touching a member), 422 L2 edges (87) and 151 L3 edges (65). All 87
    member-touching L2 edges it removes pair a member with a lncRNA record aligned antisense
    (`LAT/corrections_pass2/strand_diag.out`).
- **(5) "L2 is enough for TBC1D3" was an overclaim (important).** Replaced. In-group precision is reported next to every
  U-scope score (§0, §6).
- **(6) Soto and C_tree circularity was not stated (important).** Stated in §0 and §6. L3 is not called "exact" without
  it.
- **(7) The PKD1 / DHX40 verdicts had no named truth (important).** Tied to Soto families (§0, §5).
- **(8) L3 was not the specified test, and the summary depended on an unflagged gap convention (important).** w_98 is
  computed as specified, the gap convention is stated as a choice, and the summary bullet is conditional (§0).
- **(9) T3 overclaim; §7 framed as "more information"; readthrough wording (important).** §4 replaces the T3 paragraph.
  §7 is retitled as a threshold filtration. The readthrough example is now a node-set question (§4).
- **Minor findings, all applied:**
  - L0 uses the greedy t_P form; a union-cover closure row is added (§1, §2).
  - The E1-record restriction of f_ex is named as a choice, with the all-records count 8,129 vs 8,034 on the primary L1
    (§1). The audit's 8,562 vs 8,467 was on the 17:12 L1.
  - Both 3-truss expression views are reported (§8).
  - Double rounding is gone. The audit's cells reproduce: Soto pooled, 47 genes, 17:12 L1 = 0.574 in T-truth-U; NPIPB2–NPIPB9
    aa 0.604.
  - The L0 join is not pinned on a weak protein edge. (P ∧ aa ≥ 0.50) ∨ D is added, and the E1 0.80 / 0.50 separation
    is attributed to its thresholds (§5).
  - Split counts over non-singleton blocks, and the cross-trio statement, are added (§0, §4).
  - {TBC1D3P3, P4} = Soto ID_469; percentages replace "most" (§0, §8).
  - Named-edge slips are fixed:
    - (a) the 17:12 lncRNA L3 route was LOC105371853–TBC1D3D, and those lncRNAs now fail the strand check;
    - (b) LGALS9B chain identity is 0.907;
    - (c) the TBC1D3 literature labels are described correctly (§10);
    - (d) level names follow the definition doc.
  - Class (c) is qualified to the observed graph, and a unique-read view is added (§9).
- **One number in the definition doc was checked under the rounding fix and holds.** In §0★★★.3, ties at exactly 1.0
  (E1 identity 2,946 of 82,491; shared-exon fraction 5,138 of 84,679) give the same counts at full precision
  (`LAT/corrections_pass2/doc_numbers.out`).

## 13. Open issues

**Open issues (reconciliation pass, 2026-09-16):** no blocking finding from either track's re-check is open. Both re-checks (definition track and data track) returned ok, and every finding was minor. This pass applied the minor findings still open after the fixes. Definition track: body-record linkage, both directions and the sense condition written into t_D; "a t_D sweep moves L0–L3"; approximate L0 sizes called indicative; chr8–11 shared with the deferred layer-order study; f_ex's missing identity floor stated. Data track: union-cover L0 cut at V (5,797 genes with the 33 outside V, §1); the stale definition-doc audit note; the Soto ID_149 + ID_154 merge inside the NPIP group; non-member counts of the gap-inclusive L3 split (PKD1 itself is not in the A group, contrary to the re-check's wording); the LGALS9B identity wording; option row (ii) of §0★★★.5 now reports this run's 3-truss result. **Still open, not blocking:** T3(b) is untested, because H2 (E-values, minimap2 caps) and H3 (greedy t_P, greedy gene-body chains) are unmet. Clause 2 is approximated: gene-body records only, no genome or spliced-transcript alignment, and an exon-disjunct proxy checked against no shipped code. The w_98 gap convention must be fixed on held-out data. The readthrough node rule is proposed, not adopted. The run is not pre-registered, and the held-out run (human chr8–11 or new gorilla contigs) is not done. Leaders as a report-only grouping and a split-support L3 are not done. `corrections_rerun/rerun_checks.py` and `recheck_definition/recheck_data.py` now read the strict columns, so their saved outputs match only `LAT/pre_correction_1712/edges.tsv`. The `audit_skeptic` scripts were not re-run after the correction pass. Nothing is committed.
