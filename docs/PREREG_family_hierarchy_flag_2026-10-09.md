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
