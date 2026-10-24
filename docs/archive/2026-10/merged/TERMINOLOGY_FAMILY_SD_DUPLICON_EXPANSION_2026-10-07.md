# Terminology: multi-copy gene family vs. segmental duplication, duplicon and expansion (2026-10-07)

**Why this file exists.** Recurring disagreement with the advisor: he calls Soto 2025's catalogue and the Yoo 2025
expansions "multi-copy gene families"; the thesis author has been calling them "families within segmental
duplications" / "gene expansions". Resolution: **both are right, because the terms name different levels of one
object.** The advisor names the object; the thesis names the layers below it. The phrasing "they are expansions,
not families" is a category error (an expansion is an event, a family is its result), and that phrasing, not the
content, is what made the exchanges go badly. Adopt the wording below in the thesis, in talks and with the advisor.

## Definitions (one object, four levels)

| term | what it is | kind | evidence it is read from |
|---|---|---|---|
| **multi-copy gene family** | a set of loci that are paralogs (descended from one ancestral gene by duplication) | the OBJECT; O1 defines it | connected block of the copy graph: nodes = genomic intervals carrying reads, edges = contiguous high-coverage homology on assembly sequence (`docs/seeded_family_definition.md` §1★) |
| **segmental duplication (SD)** | a pair of genomic segments ≥ 1 kb at ≥ 90 % identity (SEDEF/BISER convention; Yoo Supp. Note XXI) | EDGE EVIDENCE of paralogy for *young* copies; no gene or read needed | assembly self-alignment |
| **duplicon** | the ancestral duplication unit, typically sub-gene (DupMasker; Vollger 2022) | a NODE CLASS below the gene | assembly + ancestral-unit library |
| **expansion** | a lineage-specific increase in a family's copy number: the in-paralogs of one outgroup gene, i.e. a lineage-restricted clade of a block | an EVENT / a CLADE; needs an outgroup and an age signal | orthology + identity/dating |

Relations that hold:

- **SD-embedded families ⊂ multi-copy gene families.** Ancient families (olfactory receptors, globins, HOX) are
  multi-copy but far below 90 % identity, so no SD caller sees them and no read is ambiguous there. Soto's and Yoo's
  sets are the **young, high-identity, SD-embedded subset**, which is exactly the subset where MAPQ-0 ambiguity exists
  and where O2/O3 are needed. "Families within SDs" is therefore a *restriction* of the advisor's term, not an
  alternative to it.
- **An expansion produces a family in the lineage where it happened.** LRPAP1: single-copy in human/chimp/orangutan,
  one expansion in gorilla, one gorilla family of 11 loci (8 full + 3 5′ fragments; Yoo counts 10). One expansion,
  one family, same thing seen as event and as result.
- **Expansion ↔ family is not 1:1, and that is the thesis content, not a terminological point.** The unit of duplication
  is the DNA block, not the gene: Yoo's gorilla chr1 unit is MAPKBP1 + JMJD7-PLA2G4B + SPTBN5 co-duplicated ×8 (one
  expansion, three families); PKD1P-NPIP fusions are one path across a duplicon boundary (one family, two blocks;
  `docs/HIERARCHY_DUPLICON_EXPANSION_2026-10-07.md`, [[project_fusions_are_duplicon_blocks]]); Soto's cut follows
  duplicon boundaries and is a cover, not a partition ([[project_cn_cut_follows_duplicons]], [[project_unit_cover]]).
  The pre-registered T1 (duplicon boundary vs family boundary) and T2 (expansions inside families) tests
  (`docs/PREREG_hierarchy_duplicon_expansion_2026-10-07.md`) are the measurements of this non-1:1 relation.
- **Family and subfamily are levels of one tree** (NPIPA/NPIPB: label-pure at every k = 2..8 under average linkage,
  no pairwise threshold separates them; [[project_family_hierarchy]]). The advisor accepted that framing; this is the
  same move one level down.

## What the sources themselves say

- Yoo et al. 2025, Supplementary Note VIII: *"There have been a number of gene family expansions in the NHPs, with
  between 1394-2056 novel gene copies found across the 184-258 families."* Both words in one phrase; "expansion" is
  the event, "family" the set being counted. The thesis author is a contributing author of this note.
- Soto et al. 2025: title *"Human-specific gene expansions contribute to brain evolution"*; body uses "gene family"
  63 times, "gene expansion(s)" 7 times, "paralog(s)" 113 times (count on the vault full text). Their objects are
  "213 human-specific gene families / 1,002 paralogs" defined on SD98 (SDs > 98 % identity) by shared exons + famCN.
- The thesis's own scoring already treats Soto's sets as families (ARI 0.7096 against them; "Soto is a refinement of
  our L2/L3"). One cannot use a catalogue as family truth and deny that its entries are families.

## The sentence to use

> "Yes, they are multi-copy gene families. Specifically they are the young ones whose paralogy is still visible as
> segmental duplications, and that is the subset where reads cannot be placed. The SD, the duplicon and the expansion
> are the layers below the family: the SD is the edge evidence, the duplicon is the sub-gene unit, the expansion is
> the lineage-specific clade. O1 defines the family; the other three are how it got there."

## How to apply

- Concede the noun, keep the structure. Say "same object, I am naming the layer below", never "no, they are not
  families".
- Open definitional discussions with one instance (LRPAP1: one expansion, one family, 11 loci, 3 fragments) and ask
  for his definition first; agreement on an instance ends the argument faster than agreement on definitions.
- In the thesis text: "multi-copy gene family" is the object O1 defines; "SD", "duplicon" and "expansion" are typed
  layers from other evidence (SEDEF; DupMasker; orthology), related to O1 families by T1/T2. Defining SD/duplicon from
  the genome alone remains out of scope (genome-only discovery, user 2026-09-13); letting SD/duplicon structure decide
  membership failed every time it was tried (register rows 676, 677, 817, 818, 997).
- Never pool human (Soto) and gorilla (Yoo) numbers when quoting either catalogue.

Related: `docs/ADVISOR_HIERARCHY_2026-10-03.md`, `docs/seeded_family_definition.md` §1★, `docs/HIERARCHY_DUPLICON_EXPANSION_2026-10-07.md`;
vault `Topics/Terminology.md`; memory [[project_hierarchy_sd_duplicon_expansion_family]], [[reference_advisor_canzar]].
