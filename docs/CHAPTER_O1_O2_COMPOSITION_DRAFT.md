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
**Edge.** Homology between units' sequences (asm20 alignment; identity ≥ 0.70, coverage of the longer
≥ 0.30, ≥ 300 bp; an additive exonic floor of 1 bp, which is where the distribution's wall is, §6dt).
**Family.** MCL over that graph, inflation 2.8, prune 1e-9. The inflation is not a tuned constant: with a
size-safe prune the anchored families are unchanged from 2.0 to 4.0 (§6ec); the earlier "cliff at 3.6" was
the prune emptying the columns of any near-uniform clique above 61 nodes.
**Certificates carried on every unit:** SD depth and core length; curated-repeat fraction of the chain
(Dfam 3.8); nearest-paralogue identity; and, per family, the CHM13 landing of its members (one human gene
stem on one chromosome for every anchored real family; scattered for every artefact; mixed 16p stems for
every duplication-block "family", §6eh).

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
**Catalog-wide (87 SD-evidenced families, 16,723 unit reads, §6es):** 29% assigned, 7% tied, 64%
ambiguous; MAPQ-60 placement agreement 95%; the MAPQ<60 pool is 7% of reads, 12% of it assigned.
Nearest-paralogue identity does not predict abstention across families (ρ −0.04; row 688).

## 6. What is a limitation and what is a defect
- Node boundaries are copy-inconsistent in the annotation and RNA cannot repair them (§6ea): the read
  chain is the transcribed unit, not the gene.
- A library-silent element (MCL4) is caught only by de novo family identity; the curated column is blind.
- Cohesion is slice-conditioned; the reported object is the genome-wide catalog.
- SD corroborates the duplication, not the gene: MCL32/MCL24 are ERV-K inside SDs (row 669).
- O2's remaining hygiene: statuses are reported for reads that overlap no unit; the catalog index must be
  read through `family_join.tsv`; 4 families abort on units without reads inside the chain.

## 7. The sentence
DNA defines the edge, RNA defines the node where the annotation is absent, SEDEF and the repeat library audit
both, and O2 assigns — or, on NPIP, declines to assign — among exactly the units O1 defined.
