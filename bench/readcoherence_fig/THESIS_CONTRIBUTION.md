# Copy-Resolved Isoform Assembly on a Spliced Family Variation Graph
*Proposed thesis contribution — rustle (long-read transcript assembly)*

## Gap
Per-molecule long-read isoform assembly is mature: "collapse" by intron chain (cDNA_Cupcake, IsoSeq, IsoQuant, FLAIR, LRAA). Graph/allele-aware methods are HAPLOTYPIC: per-molecule allelic isoforms (HapIso), or haplotype-aware quantification of KNOWN transcripts on spliced pangenome graphs (RPVG, pantas). PSV-based PARALOG resolution exists only at the DNA/gene level (Paraphase). The unoccupied cell: **de-novo, per-molecule isoform assembly resolved to paralog COPIES.**

## Setting
A gene family of unknown copy number; copies near-identical, separated by paralog-sequence variants (PSVs = fixed inter-copy base differences). Long reads multimap across copies; each read is a SPLICED PATH (intron chain) carrying an ALLELE VECTOR at PSV columns. Build a **spliced family variation graph** G: nodes = exonic segments + PSV columns (allele variants); a TRANSCRIPT = a path specifying junctions and allele choices; a COPY = a maximal allele-consistent haplotype with its isoform set.

## Problem (joint decomposition)
Given the multiset of read paths on G, find the minimum set of (copy, isoform) paths and a read assignment such that:
- (read-coherence) every emitted isoform is spanned by at least one molecule whose junctions match;
- (allele-linkage) each read's alleles are consistent along its single path;
- (copy-consistency) all isoforms of a copy share one allele-haplotype;
- (shared evidence) reads without decisive PSV evidence are apportioned across copies, not forced to one.
This is a constrained path-cover / flow-decomposition on G under linkage constraints. It FLIPS the resolve-one (facility-location) frame: multimappers are pooled shared evidence, disambiguated per-molecule only where PSVs are decisive.

## Identifiability theorem (the guarantee)
Copies c_i, c_j are DISTINGUISHABLE iff they differ at >= K PSV columns that are jointly spanned by >= m reads above the per-base error floor e. Under this condition the (copy, isoform) decomposition, restricted to spanned isoforms, is unique up to the unidentifiable residue; below it the copies are provably merged. This characterizes the recoverable regime (K independent error-agreements are improbable at long-read error rates) and turns the empirical "indistinguishability wall" into a theorem.

## What exists vs. what is new
rustle already provides the ingredients: family variation graph, PSV columns + linkage, EM read-reweighting, read-chain (per-molecule) extraction, copy certificates. NEW: the joint formulation — a single per-molecule key (intron chain (+) PSV-haplotype) decomposed on G — the identifiability theorem, and the regime characterization.

## Evaluation
Restrict to identifiable families (>= K PSVs); report copy-resolved exact-isoform (FSM) recovery vs StringTie / IsoQuant / per-copy truth; present the non-identifiable regime as a limit-of-the-data result with honest identifiability accounting.

## Positioning
- IsoQuant / Cupcake / FLAIR / LRAA : de-novo per-molecule isoform — but single-locus; paralogs out of scope.
- HapIso : de-novo per-molecule isoform split by variant — ALLELIC (2 haplotypes of one gene), not paralog copies.
- RPVG / pantas (spliced pangenome) : graph + variants — QUANTIFY known transcripts, allelic/population, not de-novo copy assembly.
- Paraphase : PSV-based paralog resolution — DNA/gene level, not isoforms.
- THIS THESIS : de-novo per-molecule isoform assembly x paralog-copy resolution, on one spliced family variation graph, with provable identifiability.

Refs: cDNA_Cupcake; IsoSeq; IsoQuant (Nat Biotech 2023); FLAIR2 (Genome Biol 2024); HapIso; RPVG / spliced pangenome (Nat Methods 2022); pantas (PLOS CB 2024); Paraphase (Nat Commun 2025).
