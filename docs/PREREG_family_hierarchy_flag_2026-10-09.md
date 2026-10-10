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
