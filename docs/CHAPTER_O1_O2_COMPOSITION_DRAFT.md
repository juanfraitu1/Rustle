# O1 → O2 by composition: families of transcribed loci, and assignment among exactly those loci (draft, 2026-09-05)

## 1. The two objects, and the one that is not the thesis's
A duplicated region of a genome is two things at once. SEDEF, or any segmental-duplication caller, sees a
**duplication block**: a mosaic of modules copied together. A gene family is a **set of loci sharing a
transcribed core**. On 16p the two disagree by construction: LCR16a carries NPIP, the adjacent module LCR16u
carries the SMG1 pseudogenes, and Johnson et al. (2001) showed they have divergent origins and concerted
expansion. The catalog must define the second object and must be able to say why it did not build the
first. Every instrument in this chapter is judged by whether it can tell the two apart on NPIP.

## 2. O1: the family as a set of units
**Node.** A unit is a locus, not an annotation record: annotation records whose exon unions overlap on a
contig are one locus (representative = longest exon union; ledger §6ef, §6er). The locus's core is the
segment linked by SD pairs to at least half of its family (§6ei, one constant, "half", used three times);
chimeric models are trimmed to it (NPIP: seven 125–308 kb models → 23–25 kb LCR16a hulls; EIF3C, an
808 bp sliver of the core, and the ABCC1-region records fall out). The unit's exon chain is read-supported:
a base is exonic when ≥ 3 reads cover it and more reads cover it than splice over it (§6en). The GFF is the
seed; the reads are the node.
**Edge.** Homology between annotation records' sequences (asm20 alignment; identity ≥ 0.70, coverage of the
longer ≥ 0.30, ≥ 300 bp; an additive exonic floor of 1 bp, which is where the distribution's wall is, §6dt), and
the alignment must map exon bases of one record onto exon bases of the other (§6ey: records are genomic spans, a
nested pseudogene carries its host's bases). The three thresholds are points inside a measured plateau — identity
0.60–0.80, coverage 0.10–0.50, 100–500 bp leave every anchored family and the Soto scores unchanged on both
substrates; the walls are at identity 0.85, coverage 0.60 and 1 kb (§6ez). Records are clustered first and
folded into loci afterwards, inside their cluster, so two records are one locus only if they overlap on exon
bases and share a cluster (§6ey; the fold-first order lost every pseudogene nested in another family's gene).
**Family.** A set of loci that share a duplicated core: each member's core is the segment of its locus linked
by SD pairs to at least half of the set. That is the object; MCL over the homology graph (inflation 2.8, prune
1e-9, smallest family 2 — a pair sharing a core is the minimum object, §6ex) is the pre-clustering that proposes candidate sets, and the core rule, SD corroboration, the repeat
library and the CHM13 landing certify them. The inflation is not a tuned constant: with a size-safe prune the
anchored families are unchanged from 2.0 to 4.0 (§6ec); the earlier "cliff at 3.6" was the prune emptying the
columns of any near-uniform clique above 61 nodes. The earlier definition (γ-quasi-clique partition of the
read-transcript graph E_r) stays runnable and is scored against this one on Soto in `O1_DEFINITION_SWITCH.md`:
band-[0.90,1) precision 0.974 vs 0.954 (CIs overlap), recall among detected pairs 0.874 → 0.940.
**Certificates carried on every unit:** SD depth and core length; curated-repeat fraction of the chain
(Dfam 3.8); nearest-paralogue identity; the duplication block each core hull belongs to and the families it is
directly SD-linked to (`blocks.tsv`, §6ez: NPIP's hulls are directly linked to LCR16u's — one block, two
families); and, per family, the CHM13 landing of its members (one human gene stem on one chromosome for every
anchored real family; scattered for every artefact; mixed 16p stems for every duplication-block "family", §6eh).

## 3. What the certificates decided on the 3-contig substrate (46 clusters with n ≥ 5)
39 SD-corroborated (3 anchored: NPIP, the CGB/NTF tandem array split at a cut vertex); 3 artefacts
(L1, ERV1, L1 — every edge a curated interspersed repeat, §6dz); 1 element-welded cluster (MCL4: a ≥ 673-
copy unclassified young element with an embedded MER1A, absent from the curated library, §6dz/§6eg);
boundaries certified for 3. The genome-wide catalog at prune 1e-9 recovers the tandem array (dissolved at
1e-5) and its three clusters above 100 members are the two repeat blobs and one 119-member protein-coding
family (§6ec).

## 4. O2: assignment among exactly O1's units
`copy_assign` consumes the unit table as its copy set — no detection, no re-derivation of the roster
(§6da–§6er). PSV columns are the star projection of the units' spliced sequences; a read is assigned when
its evidence clears `min_p < α/(n−1)` against every competitor, tied when no column distinguishes, ambiguous
when the margin is insufficient. The truth for what O2 assigns is junction-anchored reads, **audited**: of
117 anchors on NPIP, 55 were annotation gaps (the read carries the same splice on another copy at equal or
better edit distance); on the 62 valid anchors the spliced path makes no wrong confident call (GFF cores
5/5, read-chain units 3/3, §6eq–§6es), while the genomic-hull column modes do (3/8) and stay off.

## 5. The result on NPIP, and catalog-wide
**NPIP (29 loci after the locus rule; LCR16u separate):** O2 abstains on every read whose copy of origin is
at issue — 0 of 62 valid anchors assigned, all ambiguous — and places MAPQ-60 reads where minimap2 placed
them 75 of 76 times. The "4,743 assigned reads" of the annotation-level family were EIF3C, the block's
housekeeping gene, not NPIP (§6eg). The multimapping pool on the LCR16a cores is 592 of 6,545 reads and
98.6% of it is unassignable with 2–5 kb molecules: **on NPIP, O2 is an abstention certificate.**
NPIPA versus NPIPB (the advisor's Q9): the 29 gorilla loci land on CHM13 as 17 NPIPB-only, 3 on both NPIPA2 and
NPIPB13, and 9 chimeric models; identity and coverage do not separate the A-landing loci and MCL keeps all 20 in
one cluster to inflation 4.0, its first cut being exactly the four core-0 chimeras (§6ew). The A/B split is a
human 16p-position label, not a cut in this catalog.
**Catalog-wide (35 SD-evidenced families on the 3-contig substrate, paired across unit catalogs, §6eu,
Figs. 1–2):** with the hull-clipped units and the PSV read-filter, 29% assigned / 6% tied / 64% ambiguous
and MAPQ-60 placement agreement 95.3%; with units that follow the reads inside the annotated locus and the
filter off, 28% assigned / 8% tied / 63% ambiguous and agreement 99.4% (4,674/4,704): the same number of
assignments, with the confident wrong calls gone. Nearest-paralogue identity does not predict abstention across families
(ρ −0.05; row 688).

## 6. What is a limitation and what is a defect
- Node boundaries are copy-inconsistent in the annotation and RNA cannot repair them (§6ea): the read
  chain is the transcribed unit, not the gene.
- A library-silent element (MCL4) is caught only by de novo family identity; the curated column is blind.
- Cohesion is slice-conditioned; the reported object is the genome-wide catalog.
- SD corroborates the duplication, not the gene: MCL32/MCL24 are ERV-K inside SDs (row 669).
- Two O2 defects found by adjudicating the sweep's worst families (§6eu, rows 689–691), both fixed or fixable
  without a new constant: (i) the read-support PSV filter deleted every column of an unexpressed paralogue and
  left the expressed copy's own polymorphisms as the "PSVs" (ZSCAN5: 4 of 216 columns, 11 confident wrong
  calls) — off is strictly better on the sweep; (ii) the unit's chain was clipped to the SD hull, so a gene's
  main exon outside the shared core was not in its unit and O2 decided on a secondary record over the
  neighbouring copy (ZNF569-like: 190 wrong calls) — the chain now follows the reads within the annotated
  locus (Fig. 3). Still open: a molecule's secondary records are scored as observations (row 691).
- O2's hygiene closed (§6et): `in_copy` and `catalog_copy_idx` columns; unit reads counted inside the chain.

## 7. The sentence
DNA defines the edge, RNA defines the node where the annotation is absent, SEDEF and the repeat library audit
both, and O2 assigns — or, on NPIP, declines to assign — among exactly the units O1 defined.

## 8. Figures
- Fig. 1 `docs/figures/fig_sweep_status.png` — assigned / tied / ambiguous share per sweep family, three unit catalogs.
- Fig. 2 `docs/figures/fig_sweep_agreement.png` — MAPQ-60 placement agreement per family, hull-clipped units vs units that follow the reads.
- Fig. 3 `docs/figures/fig_znf569_locus.png` — the ZNF569-like locus: GFF exons, SD core hull, the clipped unit and the unit that follows the reads, over MAPQ-60 read coverage.
- Fig. 4 `docs/figures/fig_npip_anchors.png` — the 62 valid NPIP junction anchors by arm (assigned-agrees / wrong / tied / ambiguous).
Rendered by `bench/chapter_figures.py`.
