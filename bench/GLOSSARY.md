# Glossary — the terms, in one line each

Crisp working definitions for the defense. The rigorous, baseline-organized versions (paralog, segmental
duplication, multi-copy family, expansion, reference-absent copy) are in `bench/DEFINITIONS_FORMAL.md`.

## The sequence world (what a molecule is)

| term | one-line definition |
|---|---|
| **Gene copy** (copy) | one of several near-identical genomic instances of a gene, made by duplication. *In our tool: a copy = one **path** through the family's variation graph (a choice of allele at each bubble).* |
| **Isoform** | one of several transcript variants of a **single** gene, made by alternative **splicing** (different exon/junction combinations). *A path through the junctions — a different axis from copy.* |
| **Exon / intron** | exon = a segment kept in the mature mRNA; intron = a segment spliced out. |
| **Splice junction** | the donor–acceptor boundary where two exons are joined after an intron is removed. *In the graph, a branch point.* |
| **Allele** | the specific base present at a variant position on one molecule/copy. |
| **PSV** (paralogous sequence variant) | a position where the **copies** of a family differ — a column carrying ≥2 different alleles across copies. The signal that distinguishes copies. *In the graph, a bubble.* |
| **SUN** (single unique nucleotide) | a **PSV whose allele is private to exactly one copy** (Sudmant 2010) — so a single read over it pins that copy. **SUN ⊆ PSV** — every SUN is a PSV, but a PSV that only splits copies into groups (e.g. 2 vs 2) is not a SUN. |

## The evolution world (how copies arise and change)

| term | one-line definition |
|---|---|
| **Paralog** | two genes whose last common ancestor is a **duplication** event, within a lineage. *(Fitch 1970.)* |
| **Ortholog** | two genes whose last common ancestor is a **speciation** event — the "same" gene in two species. *(Fitch 1970.)* |
| **Multi-copy gene family** | a set of ≥2 paralogs **present in one genome** — a descriptive state. *Our primary object (O1).* |
| **Segmental duplication (SD)** | a duplicated genomic **block** (≥1 kb, high identity, with flanking context) — a DNA/structural object; the main mechanism producing recent paralogs. **SD98** = ≥98% identity. *(Bailey 2002.)* |
| **Expansion** | a family that **gained** copies on a lineage vs its ancestor/outgroup — a cross-species, directional change. *Needs ≥2 species to polarize.* |
| **Contraction** | the opposite — a lineage **lost** copies. |
| **Pseudogene** | a gene copy that has **lost function** (disrupting mutations and/or no expression). *Unprocessed = a duplicated copy that decayed; **processed** = a retrotransposed, intronless copy.* |
| **Retrocopy / retrogene** | a copy made by reverse-transcription of an mRNA and reinsertion → **intronless** (a signature we use to flag retrocopies). |
| **Exonization** | a formerly **non-coding** sequence (an intron, or a transposable element such as an Alu) becoming incorporated as a **new exon** — a route by which new gene structure arises. |

## Our method objects

| term | one-line definition |
|---|---|
| **Variation graph** | the family's sequence graph — a shared **backbone** (spine) + **PSV bubbles** (copy axis) + **junction branches** (isoform axis). Copies and isoforms are **paths** through it. |
| **Bubble** | a place in the variation graph where paths diverge — a PSV (copy-distinguishing) or a junction (isoform-distinguishing). |
| **Family (O1)** | a **cohesive homology cluster** of loci — each homologous to ≥ γ of the others (a γ-quasi-clique). |
| **Copy number** (χ_H) | the **fewest copy-paths that cover all the reads** — a minimum path cover. Our per-genome copy count. A **lower bound**: exon-identical copies collapse to one (the K=0 floor). |
| **Copy assignment (O2)** | which copy-path each read lies on — a max-weight **facility location**, **assign-or-abstain**, **no 1/k**. |
| **Assign-or-abstain / no 1/k** | assign a read to a copy only when a calibrated significance test passes; otherwise **abstain** (certify it unresolvable) rather than split its weight 1/k across copies. |
| **Reference-absent copy (O4)** | a copy present in the sequenced **individual** but collapsed/absent/too-divergent in the reference **assembly** (the fewest copies needed exceeds the annotated copies). |
| **Allele-specific junction (O3)** | an **allele** linked to the splice **junction** it co-occurs with on the **same read** — allele-specific splicing off single molecules, no phasing. |
| **K=0 floor** | the identifiability limit: exon-identical co-located copies carry **no PSV**, so they are genuinely RNA-unresolvable — we certify them **TIED**, not guess. |

## The one-sentence relationships (the ones that trip people)

- **SUN ⊆ PSV** — a SUN is a PSV that is private to one copy.
- **paralog** = the relationship · **multi-copy family** = the per-genome state (what we measure) · **segmental duplication** = a recent DNA mechanism · **expansion** = the cross-species change in the state.
- **copy** (paralog, PSV/bubble axis) vs **isoform** (splice variant, junction axis) — two orthogonal path-structures in the same variation graph.
- **copy number** = fewest copy-paths covering the reads (a minimum path cover / facility-location count) — *not* framed as a conflict-graph colouring.
