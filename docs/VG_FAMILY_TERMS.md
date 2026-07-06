# `vg_family` terminology: locus, isoform, copy, family

The live thesis code lives in `src/rustle/vg_family/`. This note fixes the
meaning of a few overloaded words that are used in different senses across the
module and in the companion formal notes.

> **Rule of thumb.** When you see the word **"locus"** in `vg_family`, ask
> whether it means a **gene locus** (a splice-junction community of isoforms)
> or a **physical genomic span** (a `(chrom, start, end)` interval used for the
> ≥2-distinct-loci certificate). They are not the same object.

## Canonical definitions

### Isoform

A distinct spliced transcript observed in the RNA data. Operationally it is a
specific intron chain — a ordered list of exon intervals. In code:

- `DenovoTranscript` (`family_detect.rs`) is one assembled isoform.
- `FamilyPath` (`layer2.rs`) is one recovered isoform inside a multi-copy
  family, annotated with its copy and source (`Native`, `Transferred`,
  `PsvLinked`).

An isoform is **not** a gene; one gene can produce many isoforms.

### Gene locus

The set of isoforms that are alternative splice products of the **same gene**.
In the pipeline a gene locus is collapsed from the raw assembled isoforms by
**shared splice junctions** (`family_detect::collapse_loci`): if two isoforms
use the exact same intron `(donor, acceptor)` coordinates on the same
chromosome, they belong to the same gene locus. One representative isoform is
kept for downstream family detection.

A gene locus is therefore defined by splicing, not by span overlap. This is
intentional: dense genomes would chain many unrelated genes into one bogus
locus if we used span overlap alone.

### Copy / member

One gene locus that belongs to a multi-copy gene family. Copies are homologous
loci. In the copy-assignment stage, reads are assigned to a specific copy.

### Family

A set of copies (gene loci) that are homologous. In code:

- `detect_edges` (`family_detect.rs`) builds homology edges between locus
  representatives using POA contiguous-core coverage.
- `decompose_families` (`family_split.rs`) refines raw connected components
  into cohesive families (γ-quasi-clique refinement).

### Physical span / distinct locus

The `(chrom, start, end)` interval of a gene or copy. `family_definition.rs`
uses this notion for the multi-copy certificate: a family must contain at least
**two distinct physical loci** (`distinct_loci`), where "distinct" means the
spans do not reciprocally overlap by ≥ `LOCUS_OVERLAP` (50%).

This is a genomic-coordinate check, independent of splicing. It answers the
question "do these members sit at different places in the genome?" rather than
"are these isoforms of the same gene?".

## Why two notions of locus are necessary

- **Gene-locus collapse** (shared junctions) prevents a single gene with many
  alternative splice isoforms from being mistaken for many paralog copies.
- **Physical-span distinct-loci count** prevents two nested or overlapping
  annotations of the same genomic gene from being counted as two copies of a
  family (e.g., MAGEA9 nested inside a LOC entry).

## Known limitation: isoforms with no shared junction

The shared-junction rule can be too strict: two genuine isoforms of the same
gene may have **disjoint intron sets** (for example, alternative first/last
exons that use completely different splice sites). In that case they are
currently treated as separate gene loci and may enter family detection as false
paralogs. The span-aware recovery in `family_detect::collapse_loci_span_aware`
mitigates this by additionally merging locus representatives that are
span-overlapping and either strongly contained or highly homologous.

## Mapping to the formal notes

- `DNA_FAMILY_DEFINITION_FORMAL.md` and `PROTEIN_FAMILY_DEFINITION_FORMAL.md`
  use **locus** in the physical-span sense (NCBI RefSeq gene loci).
- `family_definition_formal.md` uses **expressed locus** in the gene-locus
  sense (de-novo assembled isoforms collapsed by shared junctions).
- This note aligns the code vocabulary with the formal notes: a **gene locus**
  in code corresponds to an expressed locus in the RNA formal note, and the
  **physical span** corresponds to the DNA/protein note's locus.
