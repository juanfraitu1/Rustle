# Rustle — Algorithms & Theory

This document explains *why* Rustle works, not just *what* it does. It's aimed at readers who want to defend the method, extend it, or reason about its failure modes — in practice, advisors, reviewers, and yourself in three months.

Sections:

1. [The problem](#1-the-problem)
2. [Splice graphs](#2-splice-graphs)
3. [Network flow: formulation and why it works](#3-network-flow-formulation-and-why-it-works)
4. [Path extraction: from flow to transcripts](#4-path-extraction-from-flow-to-transcripts)
5. [Variation graphs for gene families](#5-variation-graphs-for-gene-families)
6. [Multi-mapping: the problem](#6-multi-mapping-the-problem)
7. [EM solver: derivation and convergence](#7-em-solver-derivation-and-convergence)
8. [Flow solver: two-pass redistribution](#8-flow-solver-two-pass-redistribution)
9. [Novel copy discovery: k-mer rescue of unmapped reads](#9-novel-copy-discovery-k-mer-rescue-of-unmapped-reads)
10. [SNP-based copy assignment](#10-snp-based-copy-assignment)
11. [Putting it together: the gene-family workflow](#11-putting-it-together-the-gene-family-workflow)
12. [Failure modes and honest limitations](#12-failure-modes-and-honest-limitations)
13. [Further reading](#13-further-reading)

---

## 1. The problem

**Input:** aligned long reads (PacBio IsoSeq, ONT) spanning multi-exon transcripts.

**Output:** a GTF of assembled transcripts with per-transcript coverage (reads / bp) and relative abundance.

**What makes it hard:** a gene produces multiple *isoforms* — different combinations of the same exons via alternative splicing. Reads provide noisy, partial evidence of each isoform. From a pile of reads, we must:

- Reconstruct which full-length isoforms exist
- Estimate how much each is expressed
- Avoid inventing isoforms that no real transcript supports

Long reads help (they're longer than most transcripts), but they still:
- Miss ends (incomplete 5'/3')
- Have alignment errors near splice sites
- Don't distinguish isoforms with shared 5' or 3' halves
- Cannot uniquely place reads in regions with paralogs (multi-mapping)

A **multi-copy gene family** is a cluster of paralogous genes (evolutionarily related copies) with near-identical sequence in some exons. In gorillas and other great apes, biomedically important families like olfactory receptors, amylases, TBC1D3 and the aquaporins can have 5-30 copies in the reference — and up to several additional copies not in any reference. Reads from these regions map equally well to multiple copies or fail to map at all, and standard assemblers either discard them or bias assignment to whichever copy was added first.

Rustle's two contributions:

1. A splice-graph + max-flow assembler that matches StringTie's accuracy on standard loci (faithful port of [Pertea et al. 2015](#13-further-reading), modernized in Rust).
2. A **variation graph (VG) mode** that links paralogous bundles, resolves multi-mappers probabilistically across the family, and rescues reads whose best copy isn't in the reference at all.

---

## 2. Splice graphs

A **splice graph** is a directed acyclic graph (DAG) where:

- **Nodes** are exonic segments (genomic intervals that appear in some transcript).
- **Edges** are either **junctions** (spliced intron jumps, donor→acceptor across a gap) or **contiguous** (two adjacent nodes from the same exon split at a boundary).
- Two special nodes are added: **source** (transcript starts) and **sink** (transcript ends).

```
                        ┌────────→ exon 3a ─────┐
                        │                        ↓
source → exon 1 → exon 2                         exon 4 → sink
                        │                        ↑
                        └────────→ exon 3b ─────┘
```

A **transcript** is then a path from source to sink. An isoform that skips exon 3a and uses 3b corresponds to one such path. The graph encodes *every transcript the data supports* but nothing more.

### Where the nodes come from

Two operations segment the genome into nodes:

1. **Bundle detection** finds a contiguous region of read coverage (a locus).
2. Within a bundle, every **junction boundary** splits one contiguous "bundle node" into two graph nodes. An exon used in isoform A but skipped in isoform B becomes its own node because its boundaries are junction endpoints.

Long-read mode adds a third segmentation:

3. **Longtrim** — when many reads start or end at the same position (a TSS or polyadenylation site) without a junction there, split a node at that position and add a source or sink edge. This makes the graph encode transcription boundaries, not just splicing.

### Why a *graph* and not a list of transcripts

Two isoforms that differ by one exon share every other exon. Storing them as two separate strings wastes space and forbids reasoning about their joint evidence. A graph stores each exon once and lets each transcript be a path — *combinatorial in paths, linear in data*. This is the central data-structure decision that underpins every modern transcript assembler (Cufflinks, StringTie, Scallop, IsoQuant).

---

## 3. Network flow: formulation and why it works

### 3.1 The assembly-as-flow analogy

Think of each read as a unit of *mass* (reads / bp) flowing along its splice pattern through the graph. A transcript is then a *path* from source to sink carrying some amount of flow. The sum of all transcript flows through any edge must equal the evidence observed at that edge (read coverage).

Formally: a splice graph `G = (V, E)` gains a capacity function `c : E → R≥0` where `c(u, v)` is the read-support passing from `u` to `v`. We want to find a set of paths `P₁, P₂, …` each with an abundance `a_i`, such that for every edge `(u, v)`:

    Σ_{i : (u,v) ∈ P_i} a_i ≤ c(u, v)

and the collection explains as much of the observed flow as possible. This is **flow decomposition** — writing the edge-capacity function as a sum of path-flow functions.

### 3.2 Why this works

Transcripts are, operationally, *paths through splicing*. The flow that arrives at a node equals the flow that leaves it (conservation) because every read that enters an exon must also leave it via exactly one of the outgoing splice or contiguous choices. Max-flow with source→sink conservation is therefore a **mechanistic model** of read production, not just a numerical trick.

The optimization we actually run — maximum flow from source to sink — computes the *total assembly mass the graph can support*. Decomposing that flow into paths yields a set of transcripts whose total flow is maximal (so we explain all the evidence we can) and whose individual flows give per-transcript abundance.

### 3.3 Edmonds-Karp in practice

Rustle uses the **Edmonds-Karp** variant of Ford-Fulkerson: each iteration finds a *shortest* augmenting path (BFS in the residual graph) and pushes its bottleneck capacity. Total runtime is `O(VE²)` which is ample for per-locus splice graphs with ≤ a few hundred nodes.

The twist for transcript assembly: we don't care about the max-flow *value*, we care about the **set of paths** that realize it. After running max-flow, Rustle reads off transcripts by **seeded path extraction** (next section).

### 3.4 The honest limitation

Max-flow decomposition is not unique. Two different path sets can both realize max-flow, and only one of them may correspond to real biology. Rustle (like StringTie) picks paths by a greedy heuristic — take the highest-abundance seed, extend source→sink, subtract that path's flow, repeat — which tends to favor the *dominant* isoform first. This is a reasonable prior but can miss rare isoforms that share edges with abundant ones. Section 12 spells out when this matters.

---

## 4. Path extraction: from flow to transcripts

After max-flow, Rustle extracts transcripts one by one:

1. Pick a **seed transfrag**: a high-abundance path fragment directly witnessed by a long read.
2. **Extend left** (back_to_source) following the edges with highest remaining flow until reaching a hardstart or source.
3. **Extend right** (fwd_to_sink) following highest-flow edges until reaching a hardend or sink.
4. **Subtract** this path's flow from all its edges.
5. **Repeat** with the next-best seed until no seed has enough flow.

### Why seed from long reads

A long read that spans 15 exons *directly observes* a 15-exon isoform. Using it as a seed anchors the extracted path on real evidence rather than a greedy BFS choice. This is the key reason long-read assembly beats short-read assembly at isoform-level resolution: each read is (usually) a full-length path through the graph, so there's no combinatorial ambiguity to resolve.

### Rescue and redistribution

Some transfrags don't survive back/forward extension (boundary mismatch, zero residual flow). These enter **checktrf**:

- **Redistribute**: if the transfrag's intron chain is a contiguous subsequence of an already-kept transcript, add its abundance there instead of emitting a new transcript.
- **Rescue**: if not redistributable and still has meaningful coverage vs the max-cov kept path, emit as an independent transcript.
- **Drop**: if abundance is below `readthr` or length < `min_transcript_length`, discard.

This three-way split handles the noisy residue of flow decomposition — it's where Rustle's current precision/sensitivity trade-offs live (see our tuning of `RUSTLE_RESCUE_FRAC` and `RUSTLE_NO_STRICT_REDIST`).

---

## 5. Variation graphs for gene families

### 5.1 What a variation graph is

A linear reference represents each genomic position once. A **variation graph** (VG) represents alternate sequences as parallel paths sharing nodes where they agree and diverging where they differ:

```
                ┌─→ A → T → C ─┐
 5' → → → → → → ┤               ├ → → → → → 3'
                └─→ A → G → C ─┘
```

The left path has `ATC`, the right has `AGC` — they share the `A` and `C` endpoints (edge-adjacent nodes) but diverge at the middle base (a SNP). Extend this to exon-scale and you get a graph that encodes all known alleles, haplotypes, or paralog variants in one structure. See [Paten et al. 2017](#13-further-reading), [Garrison et al. 2018](#13-further-reading).

### 5.2 Why paralogs share a VG

Paralogs arose by duplication of a common ancestor. They share exon *sequences* wherever the copies haven't diverged, and differ where they have. A variation graph of a gene family captures this directly:

- **Shared exons** → shared nodes (one node per exon, used by every copy)
- **Copy-specific exons** → parallel paths through the family-node region
- **Copy-specific SNPs within shared exons** → bubble-nodes (alternate bases)

This is useful *because reads from the family are hard to place on a linear reference*: a read covering a shared exon has equal alignment score at every copy. On the VG, that same read traverses the shared node exactly once, and its support is immediately attributable to the family — disambiguation then happens at the copy-specific nodes, where the read may or may not agree with each copy's path.

### 5.3 How Rustle discovers families

Rustle doesn't require a pre-built VG. Instead it *detects* families from the evidence:

1. **Multi-mapping links:** if many reads align equally well to bundles A and B (supplementary alignments), A and B are likely paralogs. Rustle builds a graph over bundles where an edge connects A↔B whenever they share ≥ `--vg-min-shared` multi-mappers. Connected components in this bundle-graph are candidate families.

2. **K-mer similarity:** for each family, Rustle computes the Jaccard index of exonic k-mers (k=25 by default) between bundles. High overlap confirms paralogy; spurious multi-mapping links get filtered out.

The family output (via `--vg-report`) is a TSV: one row per family with member bundles, their coverages, and the estimated copy count.

---

## 6. Multi-mapping: the problem

In a reference with *N* copies of gene `X`, a read from the shared exon of copy `k` aligns equally well to all *N* copies. Aligners typically report an `NH` tag (number of hits) and flag one as primary. Downstream tools then face a choice:

- **Use only primary:** throws away information; primary is essentially arbitrary.
- **Use all with weight `1/NH`:** spreads evidence uniformly, but the true biological distribution is almost certainly not uniform.
- **Discard all `NH > 1`:** throws away up to 50% of reads in tandem-duplicate regions.

Each of these under-uses the data. Rustle's VG mode offers three better strategies.

### The ground-truth question

When a read maps to `N` copies, where does it really come from? Three scenarios:

- **One copy is expressed.** Sequencer can't tell, but biology gives 100% to that copy.
- **Multiple copies are expressed.** The molecule genuinely came from one of them, but in shared regions we can never know which. The *probability-over-copies* is the honest answer: 60% at A, 40% at B.
- **Novel copy.** The read matches an allele that isn't in the reference at all. Assigning it to any reference copy is wrong; it needs its own new bundle.

EM and Flow handle scenarios 1 and 2 (fractional assignment). Novel copy discovery handles 3. These are the two `--vg-solver` modes plus the discovery pass.

---

## 7. EM solver: derivation and convergence

### 7.1 Setup

Let `r = 1, …, R` index multi-mapping reads, `k = 1, …, K` index copies in the family. Each read `r` has a set of candidate copies `C_r ⊆ {1,…,K}` (from alignment). We want weights `w_{r,k} ∈ [0, 1]` with `Σ_k w_{r,k} = 1` for each `r`, assigning read `r`'s unit of mass across its candidate copies.

Introduce a latent indicator `Z_{r,k} = 1` if read `r` came from copy `k`, else `0`. Let `θ_k` be the relative expression of copy `k`, with `Σ_k θ_k = 1`. Let `P(r | k)` be the likelihood of observing read `r` given it came from copy `k` — in practice, this is a function of the *junction compatibility*: how well do `r`'s splice junctions match those already assembled at copy `k`?

Complete-data log-likelihood:

    L(θ, Z) = Σ_r Σ_k Z_{r,k} [log θ_k + log P(r | k)]

We don't observe `Z`. EM alternates:

**E-step:** compute posterior `w_{r,k} = E[Z_{r,k} | r, θ] = (θ_k · P(r | k)) / Σ_j (θ_j · P(r | j))`.

**M-step:** update `θ_k = (Σ_r w_{r,k}) / R`.

Both steps have closed forms. Each iteration monotonically non-decreases the observed-data log-likelihood (Jensen's inequality — see [Dempster et al. 1977](#13-further-reading) or any graduate text). Standard EM convergence: reach a stationary point (local optimum).

### 7.2 What `P(r | k)` looks like in Rustle

For each read `r` and candidate copy `k`:

    P(r | k) ∝ exp( − α · disagree(r.junctions, k.junctions) )

where `disagree(.)` counts splice junctions in `r` that don't appear in the set of junctions assembled at copy `k`, and `α` tunes sharpness. When `α` is large, the read strongly prefers copies whose splicing matches — a read with junctions A-B-C has near-zero probability at a copy with junctions A-B-D.

This is what we earlier called *junction compatibility*: the same notion as read-to-transcript alignment in [RSEM](#13-further-reading) or [Salmon](#13-further-reading), specialized to the junction-set view that splice graphs naturally provide.

### 7.3 Convergence behavior

In our tests, EM converges within 10-20 iterations on most families. The `--vg-em-iter N` flag sets the cap (default 20). Families with extreme ambiguity (say, 8 copies with 99% shared sequence) take longer but generally stabilize within 50 iterations. We re-weight reads before reassembly, so EM runs *once per dataset* before the main assembly pass — this is why the overhead is negligible.

Key property we exploit: **EM naturally handles "belongs in both."** If copies A and B are both expressed and a read is equidistant, EM converges to weights `w ≈ θ_A / (θ_A + θ_B), θ_B / (θ_A + θ_B)`. The *honest probabilistic answer*. Winner-take-all methods (primary only, first-match) cannot produce this.

### 7.4 When EM fails

- **All-ambiguous reads**: if every read in the family is multi-mapping and no copy has *any* uniquely-mapping support, `θ` is unidentifiable — EM will converge to whatever the initialization prefers. In practice this means purely-identical paralogs with no SNPs stay indistinguishable. This is a limitation of the data, not the method; use `--vg-snp` with a genome FASTA to add per-base information, or accept that you can only assay the family in aggregate.
- **Local optima**: EM finds *a* stationary point, not necessarily the global MLE. We initialize θ from uniform and from the uniquely-mapped fraction; in practice both converge to the same answer for biological data, but adversarial constructions exist.

---

## 8. Flow solver: two-pass redistribution

The Flow solver (`--vg-solver flow`) uses the *assembly* itself as compatibility evidence.

### Algorithm

1. **First pass:** run core assembly on every family bundle with uniform multi-mapper weights (`1/NH` each). This gives a first draft of transcripts and their abundances at each copy.

2. **Re-weight multi-mappers:** for each multi-mapping read, look at the transcripts assembled at each of its candidate copies. Assign the read proportional to the coverage of the transcript(s) it supports at each copy:

        w_{r,k} ∝ Σ_{t ∈ transcripts(k) | r compatible with t} cov(t)

3. **Second pass:** reassemble with the new weights.

### Why this helps

The first pass uses only alignment geometry; the second pass uses assembled structure. A read that looks ambiguous at the alignment level may support only one isoform-structure once that structure is visible. Flow catches those cases.

### Trade-off vs EM

- Flow uses assembled *transcripts* as the compatibility signal; EM uses raw *junctions*. Assembled structure is richer (full path evidence) but requires a first-pass assembly (more compute).
- Flow is harder to debug (two passes, two sets of transcripts, multi-mapper weights that depend on the first-pass output).
- EM is theoretically cleaner (monotonic convergence); Flow is more pragmatic.

---

## 9. Novel copy discovery: k-mer rescue of unmapped reads

### The problem

A paralog that's absent from the reference — say, because the reference lineage lost it or the assembly failed there — has **no genomic coordinates in the BAM**. Reads from it either:

- Didn't map at all (flag `0x4`), or
- Got mis-assigned to a related copy that's in the reference.

Either way, no bundle is ever created for the missing copy. A linear-reference pipeline *cannot* discover it; there's literally no locus for the assembler to visit.

### The k-mer solution

From all assembled transcripts of the known family, extract k-mers (k=25 default). Build a hash set `F = {family k-mers}`. Then scan **unmapped reads** and **supplementary alignments to uncovered regions**:

- Count how many of each read's k-mers appear in `F`.
- If the fraction exceeds a threshold, mark it as a **family candidate**.

Cluster family candidates by approximate region (based on supplementary-alignment coordinates, or via k-mer overlap with other candidates). Each cluster becomes a **synthetic bundle** fed back into the standard assembly pipeline, tagged `copy_status: novel` in the output.

### Why it works when the aligner failed

The aligner is stuck in a linear-reference paradigm: a read is placed wherever its score is highest, with no awareness of "this read matches *a family*." K-mer matching ignores the linear reference entirely — it asks "does this read resemble the family, regardless of where it's supposed to land?" Matches rescue reads the aligner couldn't confidently place.

### The failure mode

K-mer matching gives no coordinate. If a novel copy is in a region that has *no* supplementary alignments anywhere (reads came back completely unmapped), the novel-copy bundle has no genomic anchor and the resulting transcript is emitted at a placeholder locus. This is a known open problem and why we currently gate novel-copy discovery behind `--vg-discover-novel` (off by default).

---

## 10. SNP-based copy assignment

When paralogs share exon sequences but differ by point mutations, alignment alone is ambiguous but **sequence variants** disambiguate.

### What we parse

For every aligned read, the BAM `MD` tag encodes the reference sequence at mismatched positions. Rustle parses `MD` into a list `mismatches = [(ref_pos, query_base), …]` per read.

### Building diagnostic variant profiles

For each position inside a family bundle, compute allele frequencies using *only uniquely-mapping* reads at each copy. A position is **copy-diagnostic** if:

- Copy A has allele `X` at frequency ≥ 0.8
- Copy B has allele `Y` at frequency ≥ 0.8
- `X ≠ Y`

Multi-mapping reads that carry allele `X` at such a position get upweighted at copy A (and downweighted at copy B) — this is a direct likelihood-ratio update on the compatibility score.

### Why this matters for EM/Flow

Without SNP information, identical-sequence copies are *indistinguishable to EM*. With SNPs, even a single diagnostic position per copy breaks the symmetry and EM converges to the correct assignment. In practice, biologically distinct paralogs almost always differ at at least a few positions — SNP mode is how we turn those differences into assignments.

---

## 11. Putting it together: the gene-family workflow

End-to-end, VG mode with all features engaged:

```
 1. bundle each locus normally
 2. detect families (multi-map links + k-mer Jaccard)
 3. build per-family variation graph:
      - shared exons = shared nodes
      - copy-specific exons = parallel paths
      - diagnostic SNPs = bubbles
 4. for each family:
      a. collect multi-mapping reads
      b. compute compatibility(r, k) from junctions (+ SNPs if --vg-snp)
      c. solve EM / Flow → weights w_{r,k}
      d. redistribute read mass across copies
 5. scan unmapped reads for family-kmer matches → novel bundles
 6. run core assembly on each (original + novel) bundle with re-weighted reads
 7. tag output with family_id, copy_id, copy_status
```

The output looks like any other GTF except that transcripts from gene families have extra attributes:

    gene_id "RSTL.1"; transcript_id "RSTL.1.1";
    cov "12.5"; FPKM "0.5"; TPM "150.0";
    family_id "FAM_3"; copy_id "2"; family_size "5";
    copy_status "reference";

A novel copy gets `copy_status "novel"`.

---

## 12. Failure modes and honest limitations

### 13.1 Flow decomposition is non-unique

Max-flow has one optimal value but potentially many path decompositions. Rustle's greedy seeded decomposition is a consistent choice but not necessarily the biologically correct one. Two known consequences:

- **Minor isoforms sharing edges with the dominant one can be missed.** They get their flow absorbed into the dominant path during seed extraction. Observable: StringTie on our benchmark emits 14 isoforms at STRG.247, Rustle emits 10 — the 4 missing are low-cov alternates whose flow was consumed by the top isoform.
- **Combinatorial over-enumeration at near-neighbor junction sites.** When two junctions share an acceptor but have donors 50 bp apart (aligner jitter vs real variant), flow can split into all combinations, producing more paths than biology.

### 13.2 EM is a local method

EM is guaranteed to non-decrease the likelihood each iteration but not to find the global MLE. Families where unique-mapping support is very uneven can converge to a local optimum that reflects the initialization rather than the data. In practice this has not been a dominant source of error on our benchmarks; on adversarial data it would require random restarts.

### 13.3 Novel-copy discovery is anchor-weak

See §9: if a novel copy has zero supplementary alignments anywhere in the reference, k-mer matching gives no genomic coordinate. The transcript still gets emitted but with placeholder coordinates. This is an open issue.

### 13.4 SNP mode requires alignment-level reliability

`MD` tag accuracy depends on the aligner's indel-correction behavior. Minimap2 with `-ax splice` on IsoSeq data is reliable; some other aligner/preset combinations produce MD tags that encode insertions as mismatches, which poisons the diagnostic variant detector. `--vg-snp` outputs a warning when diagnostic-site quality looks low.

### 13.5 Graph over-fragmentation

Our graph construction produces ~28% more nodes than StringTie's on the same reads, largely from zero-length placeholder nodes created by simultaneous End→Start events at the same coordinate. These are cosmetically clutter but do not cause transcript-level errors in our measurements. Optional collapse is available via `RUSTLE_COLLAPSE_ZEROLEN=1` (trades 1430 fewer nodes for ~5 match loss on GGO_19 — not enabled by default pending coverage-attribution investigation).

### 13.6 The parity honest statement

Rustle is within **1% sensitivity and 84% transcript precision** of StringTie on the GGO chr19 benchmark. The remaining gap is almost entirely in *flow-decomposition path enumeration* (99.2% of missed references have all their junctions in Rustle's KEEP set — the paths exist in the graph but aren't extracted). This is documented extensively in `docs/STRINGTIE_PARITY_SYSTEMATIC.md` and traced per-locus via the `RUSTLE_SEED_STATS` diagnostic.

---

## 13. Further reading

### Splice graphs and transcript assembly

- Pertea et al. (2015). **StringTie enables improved reconstruction of a transcriptome from RNA-seq reads.** *Nature Biotechnology* 33(3):290-295.
- Trapnell et al. (2010). **Transcript assembly and quantification by RNA-Seq reveals unannotated transcripts and isoform switching during cell differentiation.** (Cufflinks). *Nature Biotechnology* 28:511-515.
- Shao & Kingsford (2017). **Accurate assembly of transcripts through phase-preserving graph decomposition.** (Scallop). *Nature Biotechnology* 35:1167-1169.

### Max-flow and decomposition

- Ford & Fulkerson (1956). **Maximal flow through a network.** *Canadian Journal of Mathematics* 8:399-404.
- Edmonds & Karp (1972). **Theoretical improvements in algorithmic efficiency for network flow problems.** *Journal of the ACM* 19(2):248-264.
- Vatinlen et al. (2008). **Simple bounds and greedy algorithms for decomposing a flow into a minimal set of paths.** *European Journal of Operational Research*.

### Variation graphs

- Paten et al. (2017). **Genome graphs and the evolution of genome inference.** *Genome Research* 27:665-676.
- Garrison et al. (2018). **Variation graph toolkit improves read mapping by representing genetic variation in the reference.** (vg). *Nature Biotechnology* 36:875-879.
- Siren et al. (2021). **Pangenomics enables genotyping of known structural variants in 5202 diverse genomes.** *Science* 374:eabg8871.

### EM for multi-mapper resolution

- Dempster, Laird, Rubin (1977). **Maximum likelihood from incomplete data via the EM algorithm.** *Journal of the Royal Statistical Society B* 39:1-38.
- Li & Dewey (2011). **RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome.** *BMC Bioinformatics* 12:323.
- Patro et al. (2017). **Salmon provides fast and bias-aware quantification of transcript expression.** *Nature Methods* 14:417-419.

### Long-read transcript assembly

- Kovaka et al. (2019). **Transcriptome assembly from long-read RNA-seq alignments with StringTie2.** *Genome Biology* 20:278.
- Tang et al. (2020). **Full-length transcript characterization of SF3B1 mutation in chronic lymphocytic leukemia reveals downregulation of retained introns.** *Nature Communications* 11:1438.

### Gene family / paralog analysis

- Demuth & Hahn (2009). **The life and death of gene families.** *BioEssays* 31:29-39.
- Zhang (2003). **Evolution by gene duplication: an update.** *Trends in Ecology & Evolution* 18:292-298.

---

*This document reflects the theory behind Rustle as of April 2026. For the current empirical state of the benchmark and the open engineering issues, see `README.md` and `docs/STRINGTIE_PARITY_SYSTEMATIC.md`.*
