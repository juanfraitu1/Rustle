# Loci, copies, and isoforms — the definitional hierarchy

A recurring source of confusion in multi-copy gene-family work is conflating three things that this framework
keeps strictly separate. This note fixes the hierarchy, the orthogonality of copies and isoforms, and the two
regimes where "copy = locus" holds or bends.

## The hierarchy

```
FAMILY   = a set of paralogous LOCI                         (read-conflict component / homology)
  └─ COPY    = one LOCUS — a genomic position (chrom:start–end)   ← the unit the DEFINITION is built on
        └─ ISOFORM  = a splice path (intron chain) of that copy    ← a WITHIN-copy axis
```

- **A family is a set of loci.** The read-conflict graph's nodes are loci (the reps produced by
  `collapse_loci`, which folds *all* the transcripts/isoforms at one genomic position into one locus). A family
  is a connected component of **loci**; a "copy" in the catalog *is* a locus.
- **`collapse_loci` collapses by genomic position, not by splice pattern.** So two loci that happen to express
  the same isoform remain two separate copies — isoforms are never what defines a copy.

## Copies and isoforms are ORTHOGONAL axes

In the variation-graph view, the family is one PSV-aware variation graph with two kinds of variation:
- **PSV bubbles** = sequence differences between copies → the **copy** axis.
- **Splice bubbles** = junction (intron-chain) choices → the **isoform** axis.

A single sequenced molecule is a **path through both**: it picks a copy (a PSV path) *and* an isoform (a splice
path). The two are independent — a copy can express many isoforms, and an isoform can be shared across copies.

| feature a read carries | what it distinguishes |
|---|---|
| a **PSV** (base differs between copies) | the **copy** (genetic) |
| a **copy-specific junction** (one copy has the exon, another lacks it) | the **copy** (genetic) — used in assignment as an extra "column" (§6 junction-as-column lift) |
| an **alternatively-spliced** junction (a copy *sometimes* uses it) | nothing about the copy — it's a *within-copy* isoform choice |
| a junction **shared** by all copies of the family | nothing — non-distinguishing |

**Rule of thumb:** a copy-specific junction distinguishes copies; a shared or alternatively-spliced isoform does
not. Isoforms are reported *per copy* (the assembled splice variants at that locus); they are output, not identity.

## "The same isoform at different loci" — three cases, all handled

This is the case that prompts the worry, and it is not a problem because the definition is locus-anchored, not
isoform-anchored:

1. **Different loci, same isoform, but the copies differ in SEQUENCE (PSVs).** The read is assigned by its PSVs;
   the shared isoform contributes nothing. Resolved — the locus is recovered from the sequence.
2. **Different loci, same isoform, AND sequence-identical over what the read observes.** This is the **K=0
   identifiability floor**: the read maps equally well to both loci (MAPQ 0), and the locus is *not recoverable
   from a single RNA read*. The method certifies it **Tied** (min_p = 1, a real impossibility certificate)
   rather than guessing. The two loci exist; one molecule simply cannot say which.
3. **Different loci, DIFFERENT isoforms (a copy-specific junction).** The copy-specific junction is a
   distinguishing feature, used alongside the PSVs to assign the read.

So "the same isoform at two loci" is, by construction, *non-confusing* — it is either irrelevant to assignment
(the PSVs do the work), or it is the honest K=0 floor (correctly abstained on), never a wrong merge or split.

## The two regimes: where "copy = locus" holds vs bends

- **Reference-present (the clean case):** each paralog occupies its own genomic position, so **copy = locus**
  exactly. The read-conflict family is a set of paralogous loci; assignment routes each read to a locus.
- **Collapsed / reference-absent:** several copies share *one* reference locus — their reads all pile there at
  MAPQ 0. Here **copy ≠ locus**: the copies are **PSV-haplotypes within a single locus**, separated by phasing,
  and the locus alone cannot distinguish them. This is the regime the whole copy-assignment apparatus exists for
  (and where reference-absent copies live).

## The one substantive boundary: copy-specific junction vs alternative isoform

There is a genuinely hard distinction the framework must respect: a junction present in some molecules and absent
in others can be **genetic** (a *copy-specific* junction — between copies) or **regulatory** (an *alternatively
spliced* isoform — within a copy). Telling them apart is the **copy-vs-allele problem**:
- The allele-specific-junctions (ASJ) machinery uses per-molecule **allele→junction linkage** (the long-read
  advantage) to separate them where a heterozygous anchor exists.
- The genuinely ambiguous case — a junction that varies *between paralog copies* that are otherwise
  indistinguishable in RNA — is the **irreducible RNA boundary** and is deferred to DNA (paralog-specific copy
  number, parCN). The framework *flags* these rather than guessing.

## Summary for the thesis

- Families and copies are defined on **loci** (genomic position + sequence), never on isoforms.
- Isoforms are an **orthogonal within-copy axis**, modelled as splice bubbles / junction-columns; they are
  reported per copy and used in assignment *only* where a junction is copy-specific.
- "Same isoform at different loci" is handled by the locus/PSV anchor; the irreducible case is the K=0 floor,
  honestly abstained on.
- The only deep question is copy-specific-junction vs alternative-isoform (copy-vs-allele), addressed by ASJ
  linkage and, for the irreducible case, DNA/parCN.
