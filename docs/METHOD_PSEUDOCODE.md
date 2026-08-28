# O1 — how a multi-copy gene family is built, step by step

Every constant is the shipped default, verified against the source or against a run's own
`params.tsv` / `rule.tsv` certificate. Each step names the function that implements it, so any claim can
be checked against one place in the code.

⚠ **Two runs are quoted below and they are different objects.** The **NPIP catalog** (3 contigs,
full-depth fibroblast: 2,847 reps → 83 families) and the **gorilla testis catalog** (genome-wide:
17,924 reps → 627 families / 2,019 copies). Every figure says which.

---

## 1. Reads → exons

**Primary alignments only** — `-F 2308` drops unmapped, secondary and supplementary records.

`N` in the CIGAR is an intron spliced out of the transcript; the aligned blocks *between* the `N`
operations are the exons. `D` (deletion) is **not** an intron and does not split an exon.

*Code:* `bam::exons_from_cigar`. SAM `POS` is 1-based; exon arrays here are 0-based half-open.

---

## 2. Reads → candidate transcripts

**Two separate rules run here — one for spliced reads, one for unspliced. Both produce candidates.**

### 2a. Spliced reads — grouped by exact intron chain

Two reads belong to the same candidate iff their intron lists are **identical, junction for junction**:
same count, and every `(chromosome, donor, acceptor)` triple equal to the base. Not "similar", not
"overlapping".

- **Extent.** The candidate's start is the **2nd smallest** read start and its end the **2nd largest**
  read end (`min_terminal_support = 2`), so one runaway read cannot stretch the model.
- **The read floor is applied to the LOCUS, not to the candidate.** `pool_locus_support` is **true** by
  default: support is summed over the connected component of the **junction-incidence graph** — two
  candidates adjacent iff they share an exact `(chrom, donor, acceptor)`. A 2-read candidate therefore
  survives when its locus totals **≥ 3** (`GATE_MIN_READS = 3`). *Why:* at DAZ2 twelve reads fragment
  into nine chains whose best support is 2 — every one would die at a per-candidate threshold.
  The earlier pass-1 grouping keeps a chain at `PASS1_MIN_READS = 2`.

One locus yields many candidates: NPIPB9's 206 reads give 54 distinct intron lists, the largest holding
only 21%.

### 2b. Unspliced reads — clustered by span overlap

Unspliced reads are **excluded from 2a entirely** (an empty intron chain would pool them
chromosome-wide) and are clustered separately: per chromosome, by **single-linkage span overlap**,
threshold-free — any overlap links.

> This second rule builds the single-exon representatives, and they are **not a rounding error**:
> **1,234/2,847 = 43.3%** of reps in the NPIP catalog, **5,928/17,924 = 33.1%** genome-wide in testis.
> A single-exon candidate is its own junction-incidence component, so for it the ≥ 3 floor is simply its
> own read count.

*Code:* `pass1_skeletons_robust_with` (2a), `cluster_unspliced` (2b), `locus_support` (the pooling).

---

## 3. Candidate → exon-sum sequence

Take the **reference** bases at the read-derived exon coordinates and concatenate. **The reads give the
coordinates; the genome gives the bases.** The result — introns removed — is the *exon-sum*, and it is
the sequence every later comparison uses.

**Every junction must be canonical *and* all junctions must agree on strand**, or the whole candidate is
dropped. Canonical means **`GT..AG`, `GC..AG`, `AT..AC`** on `+`, and their reverse complements
**`CT..AC`, `CT..GC`, `GT..AT`** on `−`. An intron shorter than 4 bp is non-canonical by construction.
NPIPB12 has 8 of 9 junctions canonical and one `CT..AT` — in neither list — so a 109-read model is
discarded.

Three hard length bounds also apply here: **span ≤ 3,000,000 bp**, and **spliced length in
[100, 300,000] bp**.

> A single-exon candidate has no junctions, so this test is **vacuous** for it — it passes automatically
> and **its strand cannot be determined this way**. Historically every such rep carried a `'+'`
> placeholder.

*Code:* `build_spliced_seq`, `genome.rs` junction motifs. *(Opt-in `RUSTLE_JUNCTION_MAJORITY` tolerates a
minority of non-canonical junctions under 10 kb; **off by default**.)*

---

## 4. Two structural filters, each a conjunction

- **Readthrough.** Drop a **single-exon** candidate whose span **entirely contains ≥ 5 DISTINCT
  junctions witnessed in the primary reads**, each carried by **≥ 2** reads on the same chromosome.
  *Distinct*, not total: TSPYL1 has 51 junction observations inside its span but only 4 distinct ones.
  The rule is **read-level**, not transcript-level — phrased over assembled transcripts it misses the
  RFPL readthrough entirely.
- **Mis-chain.** Drop a candidate with an **intron > 50 kb** whose **exact junction** is supported by
  **< 3 primary reads**.

Neither of these two drops on length alone. *(The three bounds in step 3 are genuine length cut-offs.)*

*Code:* `read_junction_support`, `denovo_pipeline.rs` readthrough/mis-chain gates.

---

## 5. Candidates → loci

Two candidates go in the same locus if **any** of:

- **(a)** they share **≥ 1 identical junction** — the exact `(chrom, donor, acceptor)` triple;
- **(b)** same chromosome **and** same strand **and** their **genomic spans** overlap by **≥ 50% of the
  shorter span**;
- **(c)** same chromosome **and** same strand **and** any non-zero span overlap **and** POA
  contiguous-core coverage of the shorter **exon-sum** **≥ 0.50** — this admits alternative-first/last-exon
  isoforms whose junction sets are entirely disjoint.

⚠ **(b) and (c) measure genomic spans, introns included.** A giant intron satisfies (b) vacuously: at
NPIPB9 the selected representative has 24 aligned blocks, **none inside the gene**, joined by one
104,410 bp intron spanning the whole of it — membership by absence. `RUSTLE_COLLAPSE_EXONIC=1` switches
(b) to exons; off by default.

⚠ **This is not order-independent.** Phase (a) is an order-independent union-find over all candidates.
Phases (b)/(c) are a **fixed-point loop**: each round re-picks one representative per locus and compares
only representatives. So `pick_locus_rep` is *part of* the merge, not a step after it, and membership
depends on which representative each round selects. The run is **deterministic for a fixed input order**;
it is **not order-invariant**, and A~B, B~C does not imply {A,B,C}. *(There is no k and no clustering
algorithm — but the earlier claim that the outcome cannot depend on ordering was wrong.)*

Each locus keeps **one representative**: most reads, ties broken by longer span.

> **Every other isoform at that locus is discarded for the rest of the pipeline.** NPIPA2 has 14
> candidates and 6 distinct first exons; one survives. This is the largest known information loss in the
> method — 46% of representatives covering a known family member are single-exon stubs, and 53% of those
> loci have a discarded spliced chain carried by ≥ 3 reads.

*Code:* `collapse_loci_span_aware` → `pick_locus_rep`.

---

## 6. Loci → the E_r homology graph

**This is the step the earlier draft described incorrectly.**

### 6a. One all-vs-all, over every representative in the run — not per family

All representatives are written to **one FASTA**, and that same path is passed to minimap2 as **both
target and query**. The shipped argv, pinned byte-for-byte by a test and reproduced in each run's own
`.args` file:

```
minimap2 -c -X --no-long-join -t <threads> -k 11 -w 5   reps.fa   reps.fa
```

- **All-vs-all, not one-vs-all.** The search is unrestricted — no pre-grouping.
- **Global, not within candidate families.** *Families are the OUTPUT of this step*, so there is no
  "candidate family" to run inside.
- **Once, not three times.** Single tier (`sensitive_only = true`; the certificate records
  `additive_sensitive_tier = same-as-primary (not re-run)`).
- ⚠ **minimap2's minimizer index is a prefilter.** Only **75,272 of C(2847,2) = 4,051,281 pairs = 1.86%**
  produced any record at all in the NPIP run. So 98.1% of non-edges are *"no alignment found"*, not
  *"failed a floor"*.

### 6b. What `-X` does, and the consequence

`-X` implies `--dual=no`: each pair is emitted in **one direction only** — `(a,b)` or `(b,a)`, never both.
A pair may still carry **many records** (mean 1.58, max 1,419 in the NPIP run) and on either strand
(1,216 pairs carry both a `+` and a `−` record).

⚠ **The query is not necessarily the shorter sequence.** `--dual=no` keeps the lexicographically smaller
FASTA header — and the headers are the rep **index**, which is arbitrary with respect to length. Measured
on the NPIP PAF the query is the longer sequence in **60,747/119,282 = 50.9%** of records: a coin flip.
*(A source comment says ~60%; that figure does not reproduce on this run — quote the mechanism, not
either number.)* **Any per-record statistic must therefore choose its axis explicitly.**

Bench scripts that ran without `-X` silently built a *different* graph on byte-identical input —
partition differed on 4/14 panels.

### 6c. An edge exists iff ONE SINGLE record clears BOTH floors

| test | value | parameter |
|---|---|---|
| identity | `1 − de` **≥ 0.60** *(fallback `nmatch/blocklen` when `de:f:` is absent)* | `sensitive_identity` |
| coverage | aligned span on the **shorter** sequence ÷ its length **≥ 0.50** | `min_coverage` |
| orientation | the record must be **forward (`+`)** | `alignment_orientation` |

- Coverage form, exactly: `ql ≤ tl ? (qe−qs)/ql : (te−ts)/tl` — **the numerator's axis follows the
  denominator.**
- **ANY single record** clearing both floors makes the edge. Records are **not** summed.
- The **forward-only guard** drops every minus-strand record. It is valid only because these are
  transcript-oriented sequences. ⚠ It filters on a field that is a **placeholder** for single-exon reps
  (step 3) — 98.55% of the pairs it blocks involve such a rep.
- ⚠ `params.tsv` also prints `min_identity_asm20 = 0.80`. That is **inert here** — the asm20 tier is not
  run on this path. Do not read it as a second threshold.

*Code:* `homology_edges_all_reps_pooled` → `nucleotide_edges_scored`.

---

## 7. E_r graph → families

The partition is **connected components, then refined**:

1. Start from the **raw connected components** of the edge graph (singletons included, so the output
   partitions every node).
2. Recursively **split** any block of **≥ 3 nodes** whose induced density is **< γ = 0.20**. Blocks of
   **≤ 2 nodes are exempt from the test entirely**. Density = the fraction of possible internal edges
   present; the graph γ runs on is **unweighted** — identity and coverage decide whether an edge exists
   and are then discarded.
   > ⚠ γ is **inert on 79.11%** of the catalog: for four families in five the block *is* the connected
   > component. It matters where it matters — with γ removed, one component holds **466 copies /
   > 38 families** in the testis catalog.
   > The operation is **split-only**, which is what keeps P1 seed-invariance intact.
3. Within a block, two loci whose **spans overlap on the same strand** are merged back into one,
   **unless ≥ 3 reads map uniquely** (MAPQ > 0) to one of them — unique mappers being the evidence that
   these are genuinely two loci rather than one counted twice.
4. Keep the block as a family only if **≥ 2 distinct loci** remain (`min_copies = 2`).

A run reports both numbers: the NPIP catalog logs
`2438 γ-quasi-clique blocks → 83 families (≥ 2 distinct loci)`.

⚠ γ-quasi-clique partitioning is NP-hard, and **nothing is asserted at runtime**. γ ≥ 0.20 holds by
construction on the *pre-merge* block; the certificate the run actually writes (`n_edges`, `density`,
`lambda`, `cut_certified` in `families.tsv`) is a λ/density report on the **post-merge** locus graph, and
is not a γ certificate.

*Code:* `family_split::gamma_quasi_clique_partition`, reached via `homology_blocks_pooled_with_edges`.
γ from `RUSTLE_GENOME_GAMMA` **in the catalog binary only** — two other call sites hardcode 0.20.

---

## What the earlier draft got wrong

| the draft said | actually |
|---|---|
| "within each candidate family we run minimap2" | **globally, over all representatives at once** — families are this step's output |
| "all-vs-all **three times**, union the results" | **once**; single sensitive tier, `-k 11 -w 5`, identity ≥ 0.60 |
| "take **connected components** of those edges" | connected components **then split** any ≥ 3-node block below γ = 0.20 — and γ is **inert on 79.11%**, so "connected components" is nearly right, not wrong |
| "the outcome does not depend on any ordering" | **false for the locus merge** (step 5): phases (b)/(c) are a fixed-point loop over representatives. Deterministic, not order-invariant |
| "at least 3 reads carry that exact list" | the ≥ 3 floor is **pooled over the locus**, not per candidate |
| "every junction must be canonical" | **six motifs**, not two — and they must also **agree on strand** |

Also absent and load-bearing: the **coverage ≥ 0.50** clause, the **forward-only orientation guard**, the
**separate unspliced clustering rule** that builds a third to a half of all nodes, and the fact that
step 5 **permanently discards every non-representative isoform**.

---

# Where the method actually loses copies — measured, 2026-08-27

The steps above describe what the pipeline does. This section says which of them is responsible for what
is missed, because that was not obvious and four repair attempts were aimed at the wrong stage.

## The oracle decomposition

Substituting each stage with a perfect version, on the gorilla NPIP family (31 loci located by homology
to the 19 human copies; fibroblast reads from the assembly's own animal):

| nodes | edges | NPIP loci grouped |
|---|---|---:|
| oracle — the 31 true locus spans | **the real shipped rule** | **30/31** |
| real — what the pipeline builds | the real shipped rule | **12/31** |

⭐⭐ **GIVEN A NODE, THE DEFINITION GROUPS IT.** The edge rule owns ~**3%** of the loss; everything
upstream of it owns ~**58%**. Of the 18 NPIP loci that never become nodes: **7 are not transcribed at
all**, 7 have 1–6 reads, and 4 clear the read gate yet still build nothing because **no two reads agree
on a splice structure** (14 reads → 14 mutually conflicting chains).

⟹ **the open problem in O1 is node construction, not the family definition.** Four routes that enrich the
EDGE were measured and refuted: shared-exon matching, pooled-isoform exons, the genomic-span substrate,
and lowering the coverage floor. The last is the sharpest — it has genuinely strong edge-level evidence
(FPR +0.0018 for TPR +0.4476 against 87,990 real genomic pairs) and still cannot move NPIP off 12/31,
because a perfect edge rule cannot connect a locus that has no node.

## Two things that are true of the definition, and worth stating plainly

- **Identity never binds.** Repeat-driven cross-homology between random genomic regions is universal —
  50.5% of random 30 kb pairs align and 49.0% clear identity ≥ 0.60. The COVERAGE clause is what separates
  "shares an Alu" from "is a paralogue".
- **Completeness is penalised.** Edge formation falls as a transcript model gets more complete:
  single-exon reps form an edge at 0.1840, ≥15-exon reps at 0.0580. The mechanism is the coverage
  denominator — a stub needs half of ~2 kb, a complete model half of ~30 kb across diverged UTRs. This
  survives a repeat attack (the deficit is 3.67× after vetoing repeat-driven stub edges) and **has no
  known repair**.

## Optional modes that exist, with their measured effect

| flag | what it does | measured |
|---|---|---|
| `RUSTLE_FOOTPRINT_NODES=1` | a node may be the union of bases covered by ≥ 2 reads, exonic only, with NO requirement that reads agree on a splice structure | +1 NPIP locus (12/31 → 13/31), purity held at 3 pure families, 0 lost. **Ships OFF.** |
| `RUSTLE_READ_STRAND=1` | measure a single-exon rep's strand from read orientation instead of the `'+'` placeholder | two-sided: 972 edges gained, 864 lost; 16 antisense families correctly dissolved, 18 genuinely lost. **Ships OFF.** |

## Finding a family's members from ONE annotated copy

`bench/family_closure.py` — align a seed to the genome, use its hits as the next round's seeds, repeat to
a fixed point.

| | |
|---|---|
| gorilla NPIP, seeded from `NPIPB11` alone | converges in 3 rounds at 25/31 loci; **23/23 = 1.000 of the EXPRESSED members** |
| the pipeline's own RNA discovery | 13/31 loci, 13/23 expressed |
| human, Soto's 65-family panel | **65/65 converge**, median recall **1.000** per stratum |
| segmental-duplication-like families | 0.885 vs 0.895 for gene families — **no drift** |

⚠ It finds **LOCI, not nodes**. Feeding these windows to RNA node construction was measured and REJECTED:
it loses to the unseeded footprint pass (12/31 vs 13/31). It narrows the annotation dependency from many
seeds to one; it does not remove it.
