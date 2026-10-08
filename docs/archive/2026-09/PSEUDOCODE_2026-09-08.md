# How the pipeline goes from reads to gene families

> **Revision note for Stefan, 2026-09-08.** Five changes, each marked **[Q1]**–**[Q5]** where it answers one
> of your numbered comments. Nothing else in the text has moved: no stage was added, removed or reordered, and
> no threshold changed value. Four of the five changes only *report* something the pipeline already computed;
> the one behavioural change is in step 7, where the single γ cut becomes one level of a hierarchy that is now
> reported in full. Each change can be accepted or rejected on its own.

## What this is and is not  **[Q1 — new paragraph]**

Stages 1 to 5 build **loci**, not transcripts. For each locus exactly **one** intron chain is elected and its
exon-sum goes forward; every alternative chain at that locus is discarded. That is the opposite of what flair
and StringTie are for — their product is the isoform repertoire, and this pipeline throws that repertoire away
on purpose, because the object being defined is a genomic locus and its sequence.

Two consequences worth stating plainly. **We do not claim to assemble transcripts better than they do**, and a
comparison on isoform-level metrics would be the wrong test. What we claim is narrower: that the locus set and
the one sequence elected per locus are a sound input to family definition. **Where an assembler is a valid
substitute for stages 1–5, it should be used** — the family definition in steps 6 and 7 consumes a FASTA of one
sequence per locus and does not care how that FASTA was produced. The stages exist because the assemblers'
output needed post-processing into exactly this object anyway, not because they were judged inadequate.

## What each stage produces

  read            a nucleotide sequence with one alignment to the genome (chromosome, start position, CIGAR).
  exon            a genomic range (chromosome, start, end), read off the CIGAR.
  intron chain    an ordered list of (donor, acceptor) coordinate pairs. A key used to group reads, not a sequence.
  candidate       a group of reads sharing one intron chain, or for unspliced reads a group with overlapping spans, plus an extent (chromosome, start, end) and a read count.
  exon-sum        nucleotides: the reference bases at the candidate's exon ranges, concatenated. The reads supply the coordinates, the genome supplies the bases. Read sequences are not used after this point.
  locus           one elected candidate, carrying a range (chromosome, start, end, strand), its exon-sum sequence, and the pooled read support.
  homology edge   a relation between two loci, decided by aligning their exon-sums. Position is not used.
  graph           nodes are loci, edges are homology edges.
  component       a maximal set of loci reachable from one another by chains of edges.
  family          a component, split where it is too sparse, with co-located loci merged, kept only if at least 2 distinct loci remain.

Stages 1 to 3 work in coordinates; the exon-sum converts coordinates back into nucleotides. Position is not consulted again until the co-located merge in step 7.

## 1. Reads to exons

Primary alignments only (samtools -F 2308: no unmapped, secondary or supplementary records). An N in the CIGAR is a spliced-out intron, and the blocks between the N operations are the exons. A D is a deletion, not an intron, and does not split an exon.

**[Q2 — new paragraph] Why primary only, and what it costs.** Among equally good placements of a multi-mapping
read the primary flag is an arbitrary tie-break, so the first question is whether the result depends on it. It
does not, for 6 of 7 named families: relabelling the flag adversarially and re-running returns the same locus
set. The seventh, TSPY, drops from 5 loci to 2, and that is the honest exception — it marks the point where
the copies are too similar for any placement rule to separate them.

The restriction is nevertheless a real loss, and it is measured rather than assumed. A secondary-visible tier
recovers **26 of 31** NPIP loci against **14 of 31** for the primary-only path, so secondaries do carry copies
the primary flag hides. That tier is not in the shipped path because its **precision is unvalidated**: the same
secondaries also spill over from the anchor copy, and at the loci it adds, 21,770 of 23,807 secondary reads
carry the anchor's alleles rather than the new locus's. Primary-only is therefore a stated conservatism with a
known cost, not a claim that secondary alignments are uninformative.

## 2. Reads to candidate transcripts

Two rules, applied separately.

Spliced reads: two reads belong to the same candidate if their intron lists are identical, junction for junction, where a junction is the exact triple (chromosome, donor, acceptor).

Unspliced reads are excluded from that rule, because an empty intron list would pool them chromosome-wide. They are clustered separately, per chromosome, by single-linkage span overlap with no threshold: any overlap links. This rule produces the single-exon candidates.

A chain needs at least 2 reads to form a candidate at all. A candidate is then kept if at least 2 reads support it, and that count is pooled over the locus rather than counted per candidate: support is summed across the connected component of the junction-incidence graph, in which two candidates are adjacent if they share an exact junction. The two floors are equal, so the pooled test does not currently reject a candidate that formed. The pooling rule is what allows the floor to be raised without discarding fragmented loci, where the reads at one locus split across many chains and no single chain reaches the floor alone.

A candidate's start is the 2nd smallest read start and its end is the 2nd largest read end, so a single read cannot extend it. The trimming is per candidate, so two candidates that later join the same locus can have different extents. One locus yields many candidates.

## 3. Candidate transcript to exon-sum

Take the reference bases at the read-derived exon coordinates and concatenate them.

Every junction must be canonical and all junctions must agree on strand, or the candidate is dropped. Canonical means GT..AG, GC..AG or AT..AC on the plus strand, and their reverse complements CT..AC, CT..GC or GT..AT on the minus strand: six motifs.

Three length bounds apply: span at most 3 Mb, spliced length at least 100 bp, spliced length at most 300 kb.

A single-exon candidate has no junctions, so the canonical test is vacuous for it: it passes automatically, and its strand cannot be determined this way.

## 4. Two structural filters

Each is a conjunction, and neither drops on length alone.

Drop a single-exon candidate whose span entirely contains 5 or more distinct junctions witnessed in the primary reads, each carried by at least 2 reads. Distinct rather than total is what is counted. The rule is applied at read level, not over assembled transcripts.

Drop a candidate with an intron longer than 50 kb whose exact junction is supported by fewer than 3 primary reads.

## 5. Candidates to loci

Two candidates are placed in the same locus if any of three conditions holds:

  (a) they share at least one identical junction;
  (b) they are on the same chromosome and strand and their genomic spans overlap by at least half of the shorter span;
  (c) they are on the same chromosome and strand, their spans overlap, and the POA contiguous-core coverage of the shorter exon-sum is at least 0.50.

Condition (c) admits alternative-first-exon and alternative-last-exon isoforms whose junction sets are disjoint. Apply the three conditions repeatedly until nothing further merges. Conditions (b) and (c) form a fixed-point loop rather than a single pass, so the result is not order-independent.

The merge produces a set of candidates. The locus that goes forward is one elected member: the candidate with the most reads, ties broken by the larger span. The locus's chromosome, start, end, strand and sequence are that member's. The other members contribute their reads to the support count and are not carried forward.

A locus therefore begins and ends where its best-supported candidate begins and ends. It is not the union of its members' exons and not their outer envelope. The exon-sum passed to step 6 is one intron chain, not the locus's full exon repertoire.

The elected candidate can be single-exon even where spliced candidates exist at the same locus, because the election is on read count. For such a locus the exon-sum is the span.

The read floor in step 2 pools over the connected component of the junction-incidence graph. The locus here is built by the three conditions above, which include span overlap and POA core. These are different groupings over the same candidates and they do not coincide in either direction.

## 6. Loci to the homology graph

All representatives in the run are written to one FASTA, and that file is passed to minimap2 as both the target and the query:

    minimap2 -c -X --no-long-join -t <threads> -k 11 -w 5   reps.fa   reps.fa

The comparison is all-vs-all, not one-vs-all, and global rather than per family. There is no pre-grouping, and families are the output of this step rather than an input to it. The shipped path is single-tier, at identity 0.60 or above, with -k 11 -w 5 seeding.

An edge requires identity at least 0.60 and coverage at least 0.50 of the shorter sequence, both on a single alignment record, in the forward orientation.

The minimizer index acts as a prefilter, so most pairs produce no alignment record at all: most non-edges are "no alignment found" rather than "failed a threshold".

-X implies --dual=no, so each pair is emitted in one direction only, (a,b) or (b,a) but never both.
--dual=no keeps whichever sequence has the lexicographically smaller FASTA header, and the headers are the representative index, so the query is not necessarily the shorter sequence. Any per-record statistic must therefore state its axis.

**[Q3 — new paragraph] Runtime, and why not Jaccard.** The step is cheap because it aligns exon-sums, not
genomic spans: 961 gorilla loci totalling 2.2 Mb take **17.9 s wall, 64 s CPU, 1.1 GB**, yielding 4,786 edges.
For scale, the same command on 1,929 human chromosome-17 exon-sums (6.7 Mb) takes **11.7 s**; run on the
corresponding *genomic spans* (58 Mb) it does not finish in 38 minutes, which is the reason the object aligned
here is the exon-sum.

A MinHash Jaccard sketch over the identical sequences was measured as the alternative: **11.5 s** for all
461,280 pairs, so **1.6× faster** — not a saving that changes anything, because alignment is not the
bottleneck. It is also much less sensitive:

| Jaccard threshold | edges | alignment edges recovered | edges not in the alignment set |
|---|---|---|---|
| 0.05 | 3,116 | 64 % | 4 % |
| 0.10 | 2,828 | 59 % | 2 % |
| 0.20 | 2,469 | 52 % | 0 % |
| 0.30 | 2,175 | 46 % | 0 % |

At its most permissive it misses **a third** of the edges. The reason is structural, not tuning: Jaccard is a
*global containment* measure over the whole sequence, while the edge rule is *local* — 60 % identity over half
the shorter sequence on one alignment record. Two loci that share one conserved region and differ elsewhere
score low on Jaccard and are genuine edges. Jaccard is a reasonable prefilter and a poor edge criterion.

## 7. Homology graph to final families

Take the connected components of the edges.

A connected component is a maximal set of loci in which every locus is reachable from every other by a chain of edges; maximal means no further locus can be added without breaking that. Every locus lands in exactly one component, so the components partition the loci.

Membership is therefore transitive. If A-B and B-C are edges, then A, B and C form one component even when A and C have no edge between them, no measured identity and no alignment record. A component is held together by chains, not by pairwise resemblance. A single spurious edge merges two whole components rather than two loci.

Then refine the components.

The refinement in full. For a set of loci S, the induced subgraph is S together with those edges whose endpoints both lie in S. Its density is

    density(S) = 2m / (n(n-1))

where n is the number of loci in S and m the number of edges in the induced subgraph. The denominator n(n-1)/2 is the number of possible pairs, so density is the fraction of possible pairs that are edges: 1.0 if every locus is joined to every other, and near 0 for a long chain.

The rule is: a block of 3 or more loci whose density is below gamma = 0.20 is split, and each part is tested again, until every surviving block either has 2 or fewer loci or has density at least gamma. Blocks of 2 or fewer are exempt, because a single edge between 2 loci already gives density 1.0 and the test cannot fail.

The split itself is Louvain community detection on the induced subgraph, run at resolution 1, 2, 4 and 8, taking the first resolution that returns 2 or more communities. If none does, the block is split into its connected components instead.

**[Q5 — changed] The graph carries identity.** Previously the graph was unweighted: identity and coverage
decided whether an edge exists and were then discarded, so every edge counted the same in the density. Each
edge now keeps its alignment identity, and two numbers are reported for every block alongside its density: the
**mean** and **minimum** pairwise identity among its edges. Density says how *connected* a block is; identity
says how *similar* its members are, and a block can score well on one and badly on the other. The split rule
itself is unchanged — this adds a reported quantity, not a new criterion.

**[Q4 — changed] Report the hierarchy, not one cut.** The recursion above stops at a single value of gamma and
emits one partition. It now records **every level it passes through** and emits the whole nesting, so gamma is
a level that can be read off the hierarchy rather than a threshold the answer depends on. Because each edge
now carries identity, the natural axis for that hierarchy is identity itself: report the partition of each
component as the identity floor rises, which needs no resolution parameter and no arbitrary cut.

*Worked example, TBC1D3 (human, CHM13, chromosome 17 exon-sums, 1,929 loci, 805 edges, 135 components).*
The component containing TBC1D3 holds **32 loci: 11 named TBC1D3 records, 5 USP32 records, and 16 others**.
That is not an over-merge — TBC1D3 arose inside a USP32 duplication, so the 17q12 block containing both is
the expected object at a permissive identity floor.

| identity floor | blocks | largest | what the largest blocks contain |
|---|---|---|---|
| 0.60 | 1 | 32 | 11×TBC1D3 + 5×USP32 together — the duplication block |
| 0.70 | 3 | 23 | 10×TBC1D3 + 5×USP32 |
| 0.90 | 8 | 22 | 10×TBC1D3 + 5×USP32 |
| **0.95** | 11 | 13 | **TBC1D3 (10) separates from USP32 (6)** |
| **0.98** | 19 | 6 | TBC1D3, B, F, G, H, I \| TBC1D3D, K + LOC102724956 \| USP32P1/P2/P3 |
| 0.99 | 23 | 6 | the same 6-member TBC1D3 core persists |

Read down the column: the family is one object at 0.60, separates from its parent USP32 duplication at 0.95,
and splits into two TBC1D3 subgroups at 0.98 — which is the identity threshold Soto 2025 uses to define a
segmental duplication. The hierarchy makes that visible; a single cut would have reported only one of these
levels and hidden the rest.

The graph the refinement runs on is unweighted for the density test; identity is carried alongside (see **[Q5]**).

The refinement only ever splits. It starts from the connected components and never merges across them, and no locus moves between components. This is what makes the output invariant to the order in which loci are supplied and to the choice of starting locus.

A block that passes the test is a gamma-quasi-clique: a set of loci in which at least gamma of the possible pairs are edges. Finding a maximum such set is NP-hard, so the block produced here is a certified witness that a set meets the bar, not a proof that it is the largest set that does.

This step occupies the same position that Markov clustering occupies in TribeMCL, OrthoMCL and OrthoFinder: it partitions an all-vs-all homology graph into families. Taking the connected components alone, without it, chains subfamilies together through shared-domain hubs.

Within a block, two loci that overlap on the same strand are merged into one unless at least 3 reads map uniquely (MAPQ above zero) to one of them, which is the evidence that they are two loci rather than one counted twice.

Keep the block as a family only if at least 2 distinct loci remain.

---

## Summary of the five changes

| # | your comment | what changed | kind |
|---|---|---|---|
| Q1 | overlaps with flair / StringTie | new opening paragraph: stages 1–5 build loci, not transcripts, and an assembler may be substituted for them | wording only |
| Q2 | why primary alignments only | new paragraph in stage 1 with the invariance measurement (6/7 families), the exception (TSPY 5→2) and the measured cost (14/31 vs 26/31) | wording only |
| Q3 | runtime, and versus Jaccard | new paragraph in step 6 with wall/CPU/memory and the measured Jaccard comparison | wording only |
| Q4 | report the hierarchy, show it for TBC1D3 | step 7 emits every level of the recursion instead of one gamma cut, with the TBC1D3 table | **behaviour** |
| Q5 | similarity as well as density | each edge keeps its identity; mean and minimum identity reported per block | reporting; split rule unchanged |

**Still open, and not addressed here:** a direct comparison of stages 1–5 against flair or StringTie on the
same reads. Q1 above argues the object differs; it does not yet show a measurement.
