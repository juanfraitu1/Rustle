# Family hierarchy under a flag (SD -> family -> expansion -> duplicon, with fusion loci and external families), 2026-10-09

**Pre-registered before any hierarchy or Soto-containment number is computed.** User goal (2026-10-09): put segmental duplications, duplicons, expansions (including fusion loci) and gene families into one clean hierarchy, under a flag, and integrate or show that Soto's families are a subset of ours.

## Principle (from `docs/TERMINOLOGY_FAMILY_SD_DUPLICON_EXPANSION_2026-10-07.md`)

The family (O1, `mcl_families`) is the object and stays a PARTITION of loci. SD, duplicon and expansion are typed layers read from other evidence; they ANNOTATE families and never change membership (every attempt to let SD/duplicon structure decide membership failed: register rows 676, 677, 817, 818, 997). The hierarchy is a partial order: a family can sit in several SD components, a fusion locus relates two families, an external family can straddle expansions.

## The flag

`mcl_families ... --emit-hierarchy --hierarchy-duplicons <DupMasker colors BED> [--hierarchy-sd <SEDEF BED>] [--hierarchy-external <TSV>]`. Default off; every existing product byte-identical with it on. Writes `<out>.hierarchy.members.tsv`, `.nodes.tsv`, `.edges.tsv` and, with `--hierarchy-external`, `.external.tsv`.

## Levels and definitions (all threshold-free except >= 1 bp overlap)

1. **Member** = a clustered locus of `<out>.clusters.tsv`; its exons = the union of the exon blocks of every annotation record folded into it (`loci.tsv`).
2. **Duplicon composition.** For each exon, the bases covered by each DupMasker duplicon id (segments of one id merged first); the exon's **dominant duplicon** = the id covering most of its bases (ties: smallest id); exons with no duplicon base have none. The member's duplicon string = its exons' dominant ids in transcript order, consecutive repeats collapsed.
3. **Owner of a duplicon** = the family whose members' exons carry the most bases of it (ties: smallest family id).
4. **Fusion locus** = a member some of whose exons are dominated by duplicons owned by OTHER families; `fusion_of` = those families. It stays in its own family (partition) and gains a `fusion_of` relation (cover), as decided for the container (09-30).
5. **Expansion** (copies of one duplicated block) = within a family, the member's **signature** = the set of dominant duplicons of its exons that the family owns. The family's expansions are its MAXIMAL signatures under inclusion; a member belongs to every maximal signature containing its own (exactly one: assigned; several: a shared fragment, listed in each); a member with an empty signature belongs to the family's `unique` class. Lineage specificity (an outgroup) is not tested by this flag; the level is the duplicated block, "expansion candidate" in the strict terminology.
6. **SD component** (with `--hierarchy-sd`) = SEDEF segments merged by overlap into SD loci, joined by the SEDEF pairs into connected components; a member sits in every component one of whose SD loci overlaps its exons (>= 1 bp); per member the maximum pair identity touching it (SEDEF column 21).
7. **External family** (with `--hierarchy-external`: `chrom start end label`, GFF coordinates of annotation records, one row per record x label) = mapped to members through the fold table. Per external family: records, records on clustered members, the best family and **containment** = |E ∩ F| / |E on clustered members|, the best expansion and its Jaccard over members.

## Test on Soto 2025 (human, CHM13 v1.0, the F_O universe of `docs/PREREG_hierarchy_duplicon_expansion_2026-10-07.md`)

Inputs: `/mnt/linuxdisk/tmp/ourhom/{loci.paf, loci.gff3}` with the parameters of `ours.params.tsv`; DupMasker colors v1.0; SEDEF `final_v1_clean.bed`; Soto S1C (`bench/soto/soto_famCN_S1C.tsv`) by Gene ID -> record coordinates (`names.tsv`), never by name. Soto family = an S1C `Family ID` other than Unassigned.

- **H0 (regression).** `clusters.tsv`, `loci.tsv`, `pairs.tsv` byte-identical to the run without the flag.
- **H1 (Soto is a subset of ours).** Among Soto families with >= 2 records on clustered members: **>= 90% have containment 1** (every such record in one of our families) **and the record-weighted mean containment is >= 0.95.** Reported: the exceptions with their families, and catalog coverage (the share of Soto-family records that lie on any clustered member).
- **H2 (where Soto sits, descriptive).** Per Soto family: whether it equals one of our expansions, a union of expansions of one family, or cuts across expansions; median best-expansion Jaccard.
- **H3 (positive control for fusion loci).** The PKD1P-NPIP and PDXDC2P-NPIPB14P readthroughs present among clustered members are fusion loci whose `fusion_of` includes the family holding PKD1 / PDXDC1 duplicons (09-30 descriptive note `project_fusions_are_duplicon_blocks`). Reported per gene found; not a bar if none is clustered.

**Decision.** The flag ships (opt-in) iff H0 holds. "Soto's families are a subset of ours" is stated iff H1 holds; otherwise the containment numbers and exceptions are reported as they are.

## Result of the registered run and Amendment 1 (2026-10-09, after H0-H3, before any number of the amended levels)

**Run.** `mcl_families --paf loci.paf --gff loci.gff3 --min-exonic-bp 1 --min-shared-exon-frac 0.60 --dump-pairs` (the registered F_O), with and without `--emit-hierarchy --hierarchy-duplicons <colors v1.0> --hierarchy-sd final_v1_clean.bed --hierarchy-external <S1C>`. 1,709 members, 380 families, 430 expansions, 475 members with `fusion_of`, 1,240 duplicons owned. Soto: 605 families, 2,572 S1C rows, all mapped by Gene ID; 1,973 rows (76.7%) on clustered members; 157 Soto families have no member in our families, 48 one member.
- **H0 met.** `clusters.tsv`, `loci.tsv`, `pairs.tsv` byte-identical with and without the flag; the new build reproduces the registered F_O `clusters.tsv`.
- **H1 NOT met at the MCL-family level.** Of 400 Soto families with >= 2 members, 313 (78.2%) are fully contained in one of our families (bar 90%); member-weighted mean containment 0.892 (bar 0.95). The exceptions are Soto families spread over 2 to 5 of our MCL families (ID_151, ID_152, ID_153, ID_155, ID_212, ID_226, ID_78 at 1/3; ID_99, ID_100 at 0.40; ...): MCL at inflation 2.8 cuts some homology groups more finely than Soto does.
- **H2 (descriptive): Soto families sit at our expansion level.** Median best-expansion Jaccard 1.000; **215 of the 400 Soto families are exactly one of our expansions**, 343 at Jaccard >= 0.5.
- **H3 met in substance.** PKD1P6-NPIPP1: member of the PKD1 family (MCL28), duplicon string SD9473, SD9474, SD9609, SD9607, SD9613 (PKD1's) -> SD9449, SD9443 (NPIP core), `fusion_of` MCL9 (the family owning the NPIP core). PDXDC2P-NPIPB14P: PDXDC family (MCL286), SD9585, SD9479 -> SD9456 ... SD9443, `fusion_of` MCL9. The fusion links the two predicted families; the prediction's wording assumed the readthrough would sit in the NPIP family, here it sits in the PKD1 / PDXDC family and points to NPIP.

**Amendment 1 (registered now).**
(a) *Fusion typing (post hoc to the 475 count, disclosed).* A `fusion_of` member is a **fusion** when its exons are dominated both by duplicons its own family owns and by duplicons another family owns; it is **foreign-only** when every duplicon-dominated exon belongs to another family (the family boundary runs inside a duplicon, e.g. a member made of one duplicon its family does not own). Both are reported; only fusions are called fusion loci.
(b) *Superfamily level.* Above the MCL family: the **superfamily** = the connected component of the pre-MCL homology graph (the graph MCL partitions; nothing new is computed), so family ⊆ superfamily by construction. **H1' (Soto is a subset of ours at the superfamily level):** among Soto families with >= 2 members, >= 90% fully contained in one superfamily and member-weighted mean containment >= 0.95, **counted only if the superfamilies are not trivial: the largest holds <= 10% of the members** (otherwise reported as uninformative). Reported: superfamily size distribution, how many MCL families each contains, the Soto families still split.

## Amendment 1 result and Amendment 2 (2026-10-09, before any number of the amended SD level)

**Amendment 1 run** (same command, rebuilt binary; H0 holds again: `clusters.tsv`, `loci.tsv`, `pairs.tsv` byte-identical). 331 superfamilies (largest 65 of 1,709 members = 3.8%, guard met; median size 3; 32 hold more than one MCL family). **H1' NOT met: 346 of 400 Soto families (86.5%) fully contained in one superfamily (bar 90%), member-weighted containment 0.940 (bar 0.95).** Fusion typing: 319 fusion loci, 156 foreign-only members. S1C carries 114 `Unassigned_*` labels; none has >= 2 members, so the H1/H1'/H2 counts over the 400 families are unaffected (491 Soto families proper).
**Why 54 Soto families still straddle superfamilies (descriptive, after the numbers).** Of the 642 cross-superfamily member pairs inside them, 355 align in the all-vs-all (median identity 0.969, median block 1,983 bp) and are removed by the exonic filters of the registered DNA definition (`--min-exonic-bp 1 --min-shared-exon-frac 0.60`: the shared sequence is intronic or a small part of the exons); 287 do not align at all. Soto joins genes on partial exon sharing that our family definition rejects on purpose (LORO-validated, `docs/seeded_family_definition.md` §6o2-§6o6).
**The SD level as registered is uninformative.** SEDEF pairs chain 1,598 of the 1,709 members (93.5%) into one component (`SDc0`); every member touches some SEDEF segment.
**Amendment 2 (registered now).** The SD level becomes the **SD block** = a SEDEF segment locus (segments merged by overlap on one chromosome, as already built), not the pair-connected component. A member sits in every SD block its exons overlap. Per family, the SD edge evidence: the share of member pairs directly joined by one SEDEF pair (a pair with one segment in an SD block of each member), and the highest identity of such a pair. Components stay in the members table as a column for reference. Reported, no bar: SD-block length distribution, members per block, the share of families whose member pairs are all SD-joined, and how many SD blocks hold members of more than one family (co-duplicated blocks).

## Amendment 2 result and summary (2026-10-09; tables in `docs/hierarchy_results/`)

**Run** (same command, rebuilt binary). H0 holds again (`clusters.tsv`, `loci.tsv`, `pairs.tsv` byte-identical; the external placement table is unchanged). **SD blocks:** members sit in 322 SD blocks, median length 139,972 bp (90th percentile 872,215; maximum 5,489,270 where SEDEF segments overlap into one region), median 2 members per block; **192 blocks hold members of more than one family** (co-duplicated blocks). **SD edge evidence:** the share of a family's member pairs joined directly by one SEDEF pair has median 1.000; all pairs joined in 200 of 380 families, none in 47 (families below SD identity). Example: the NPIP family MCL9 (19 members, 3 expansions, 8 owned duplicons, 165 of 171 pairs SD-joined, best identity 0.998), the PKD1 family MCL28 (7 members, 19 of 21, 0.994) and the PDXDC family MCL286 sit in one superfamily (SF866); PDXDC2P-NPIPB14P is a fusion locus of MCL286 with `fusion_of` MCL9.

**The hierarchy, as shipped (opt-in).**

| level | object | defined by | relation to the level above |
|---|---|---|---|
| SD block | a merged SEDEF segment locus | SEDEF (assembly self-alignment) | a family sits in one or more blocks; a block can hold several families (192 here) |
| superfamily | a connected component of the pre-MCL homology graph | O1 graph (no new evidence) | family ⊆ superfamily (partition of the partition) |
| family | an MCL cluster (O1) | O1 | members partition the loci |
| expansion | copies of one duplicated block: a maximal duplicon signature of the family | DupMasker + ownership | expansion ⊆ family; shared fragments sit in several (cover) |
| duplicon | DupMasker id | DupMasker | composes an expansion's signature; owned by one family |
| fusion locus | a member with exons dominated by its own and another family's duplicons | duplicon ownership | stays in its family, `fusion_of` the other (319 here; 156 more members are foreign-only) |
| external family (Soto) | S1C family | Soto 2025 | `contained_in` family and superfamily, `best_expansion` |

**Soto in our hierarchy.** Soto's families are **not a strict subset** of our families (H1: 78.2% fully contained in one MCL family; H1': 86.5% in one superfamily; bars 90%), because Soto joins genes on partial exon sharing that our DNA definition rejects on purpose (355 of the 642 straddling pairs align, median identity 0.969, and fail the exonic filters; 287 do not align). **Where they do sit: at our expansion level** (H2): 215 of the 400 Soto families with >= 2 members are exactly one of our expansions, median best-expansion Jaccard 1.0, 343 at >= 0.5. The registered wording for the thesis: "Soto's families sit at the expansion level of our hierarchy (215 of 400 identical), inside our families for 78% and inside our superfamilies for 87%; the remainder are joined by partial exon sharing our family definition does not accept."

## Amendment 3 (2026-10-09, before any number of the new level): the SD complex, and a complete typed placement of every Soto family

**Why.** H1 and H1' place 78.2% / 86.5% of the Soto families inside one of our families / superfamilies; the rest straddle superfamilies because Soto joins genes on shared SD98 exons that our exonic filters reject. Soto's families are built on SD evidence, so the level that should contain them is the one where our SD evidence joins superfamilies.
**New level (in the flag).** **SD complex** = connected component of the superfamilies under the relation "some member of one and some member of the other are joined directly by one SEDEF pair" (a pair with one segment in an SD block of each member; no transitive chaining through blocks without members). Hierarchy: SD complex ⊇ superfamily ⊇ family ⊇ expansion ⊇ duplicon; SD blocks stay the spatial layer. Written: `sd_complex` per member, `sd_complex` nodes, `superfamily -> sd_complex` and `superfamily -- superfamily sd_linked` edges (member pairs, best identity), and complex containment of external families.
**H4 (Soto is a subset of ours at the SD-complex level).** Among Soto families with >= 2 members: >= 90% fully contained in one SD complex and member-weighted containment >= 0.95, **counted only if the largest SD complex holds <= 10% of the members** (else uninformative, as the pair-connected components of Amendment 2 were).
**H5 (complete typed integration, descriptive).** Every Soto family with >= 2 members gets exactly one type, the lowest level containing all its members: EXPANSION (equal to one expansion), WITHIN-FAMILY (inside one family), WITHIN-SUPERFAMILY, WITHIN-SD-COMPLEX, OUTSIDE (spans complexes); plus, for every Soto family, the count of its records not on any clustered member (not in our families at all). The table of counts and the OUTSIDE list are reported.
**Decision.** "Soto's families are a subset of ours" is stated at the level where H4 holds (with the guard); H5 is the integration whatever H4 says.
