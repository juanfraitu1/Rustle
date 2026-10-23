# O3 as a port of Strainy's phasing plus IsoCon's statistics (design, 2026-10-23)

**Status:** the user approved the direction on 2026-10-23: "we are porting strainy's method right, not just calling it? ... O3 might
be a mixture of some of strainy's behavior plus some of isocon's"; then "yes please", with gorilla first. The advisor called Strainy
"very important" at the 2026-10-23 meeting. Decision rules: `docs/PREREG_o3_port_2026-10-23.md`, committed before any run.

## 1. Why
- **The O3 question.** O3 (RNA-only, one reference: `feedback_o3_rna_only`) has to show where this sample's reads need a copy that
  the reference lacks. It must then give that copy its reads back, so that O2 stops forcing them onto a reference copy.
- **Why the current stage fails.** The 2026-10-23 HG002 pilot (`docs/O3_FIRST_PILOT_2026-10-23.md` §4) shows how. SPDYE1's hidden
  copy is 99 reads at 9% of the pile on SPDYE1.18. It is separable only by alleles linked across columns.
  - The shipped stage `candidate_copies` (IsoCon-like, clusters cut at delta 0.96%) merges it.
  - Strainy run off the shelf merges it too: its allele-frequency floor is 0.2.
- **What each source contributes.** Strainy brings the structure (reference frame, read graph, refinement, no preset count) and
  IsoCon the statistics (a low-abundance variant is real by a test, not by a frequency floor; reference-free consensus).
- **Why gorilla first.** On KB3781 the reference is the same animal's genome, and its mother's and father's haplotypes give
  every copy's two alleles. A copy missing from the reference can only be one we mask, or one carried on one haplotype only.

## 2. The method (one family at a time; input = O1's family copies + the family's read net, as `candidate_copies` builds it)
1. **Piles (Strainy's collapsed unitig = a reference copy).** A read joins the pile of every copy holding one of its
   alignments within 10 AS of its best: O2's tie rule, so near-identical copies share reads. The pile is the read's aligned
   bases on the copy's exons, in genome coordinates. Bases come with their quality; a base below Q10 counts as missing.
2. **Real columns (IsoCon's test in place of an allele-frequency floor).** A pile column is a candidate when >= 2 reads carry
   the same non-reference base.
   - **The test.** The allele is real when its count is improbable under sequencing error alone: a Poisson-binomial tail on
     the carriers' error probabilities from their base qualities, Bonferroni over the pile's candidate columns, alpha 0.01.
   - **Editing.** A>G / T>C columns are kept but marked; the report says how many clusters depend on them.
3. **Read graph (Strainy, clusteringDivergence 0).**
   - Two reads of a pile are linked when they share >= 1 real column and disagree at none.
   - Reads that share no real column are neither linked nor in conflict. That removes the partial-coverage false alarms of the
     2026-10-22 spike's tier 1.
4. **Clusters.** Deterministic semi-synchronous label propagation (networkx `label_propagation_communities`) gives the
   clusters. A cluster needs >= 2 reads; the rest are unclustered.
5. **Refinement (Strainy).**
   - Inside each cluster, step 2 is recomputed on the cluster's reads alone.
   - If any column is real inside it, the cluster's graph is rebuilt on those columns and re-clustered.
   - This repeats until no cluster splits.
   - A cluster's haplotype is its majority base at every real column it covers.
   - Clusters whose haplotypes agree at every column they share are merged. That joins one haplotype that was split where the
     reads show no heterozygosity.
6. **Certificate (allele slots; Hall's condition).**
   - **Slots.** Each copy has 2 slots, or 1 for a non-PAR chrX/chrY copy (KB3781 and HG002 are male).
   - **Candidate edges.** A haplotype may take a slot of every copy whose pile holds >= half of its reads.
   - **Matching.** A maximum b-matching assigns haplotypes to slots.
   - **The flag.** Haplotypes left unmatched are certified extra copies: the reads demand more sequence than the reference
     copies they align to can carry. That is the copy the reference lacks.
7. **Candidate sequence (IsoCon's consensus, in the reference frame).**
   - The certified haplotype's copy (the pile holding most of its reads) is taken as its exon union.
   - The cluster's consensus is applied to it: the majority base at every position >= 2 members cover, and the majority indels.
   - That patched sequence is the candidate's contig.
8. **Products.** The same products as `candidate_copies` (`candidates.tsv` with `flagged` = certified, `contigs.fa`,
   `clusters.tsv`, `reads.tsv`, `families.tsv`). The existing acceptance harness (`bench/rna_allele/accept_o3_candidates.*`) then
   scores the port on the same panels with the same bars: realignment to genome + contigs, D right, false moves, the mat / pat
   classification.

**Constants.** All are project conventions already in use: >= 2 reads, Q10, 10 AS, alpha 0.01. There is no delta and no
allele-frequency floor.

**Not in v1** (deferred, written down here):
- Reads with no alignment to any copy (unmapped, or attributed by the net's pass B) take no part in step 1. A reference-free
  path for them (IsoCon proper) is v2. The 10-01 IsoCon chain and `candidate_copies` already cover that path.
- The port is not combined with `candidate_copies`' own candidates in v1. It is scored alone, so its effect is measured
  alone.

## 3. Test bed (gorilla KB3781 fibroblast Iso-Seq on its own `_pri` genome; truth from mat / pat), in order
1. **Development 1, the deletion panel** (Amendment 7's 53 families, work dir `rna_allele/a13a15`: `R.bam` on the masked
   genome, `net_reads.tsv` with role D/S). Measures: D right against the comparator C, and false moves.
2. **Development 2, the no-deletion control** (the same 53 families on the unmasked genome, `a14a15`). Measures: families with
   a false flag, each flag classified against mat / pat (a haplotype-only / b allele / c unmatched / pri).
3. **Held-out** (`a15h`: Amendment 10's 30 disjoint families). Run once, with code frozen after the development runs.
4. **After gorilla:** HG002 (CHM13) as a second substrate. Its numbers are never pooled with gorilla's. Report-only there:
   SPDYE1's hidden copy (does the port certify it, and do its 99 reads land on the candidate?).

## 4. What is left to the builder
- Python prototype first: `bench/o3_graph/o3port.py` with tests. A Rust port with parity tests comes only if it is adopted.
- Execution: native, by the session (the user asked to reduce approval gates; one whole-branch review at the end).
