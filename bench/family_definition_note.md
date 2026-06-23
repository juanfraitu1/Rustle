# A self-contained RNA-level definition of multi-copy gene families

**Two relations on expressed loci — read-confusability and backbone-homology — and the family is where both hold.**

This note states the definition in closed form. The empirical record (panel ledger, genome-wide run, survey of
rejected alternatives) is in `family_definition_formal.md`; this note is the formal object and is self-contained
given one named input (the locus set, §1).

---

## 0. The object, in one sentence

> A **multi-copy gene family** is a maximal set of expressed loci that is both **read-coupled** (its loci are mutually
> confusable to the aligner) and **backbone-connected** (its loci share a common gene backbone).

The two italicised conditions are two relations defined below; a family is a connected component of the second
relation restricted to a coupled group of the first. Neither relation alone is the right object — that is the point
(§4).

---

## 1. Primitives

Fix a set of spliced long reads aligned to a genome.

- **Vertices $V$.** The **de-novo loci** emitted by the upstream read-coherence assembly: each vertex is an expressed
  transcription unit — a maximal set of reads sharing splice structure at a genomic position. $V$ is therefore
  *read-derived and label-free* (no annotated gene is consulted), but it is an output of that named procedure, not a
  closed-form predicate in this note; the structural claims below that rely on tight per-locus vertices (P6) are
  stated as conditional on it.
- **Copy model $S(v)$.** The **all-isoform exon-union** of a locus $v$: every exonic base any read at $v$ aligns
  through (a retained intron appears where reads span it). A *copy* is a locus, and its isoforms **aggregate** into one
  $S(v)$ — so alternative splicing within a copy is one model, not several copies.
- **Divergence.** Each alignment of a read $r$ to a locus carries gap-compressed per-base divergence $d(r,v)$
  (minimap2 `de`). Each alignment **record** is attributed to its single best-overlap locus — ties broken toward the
  *smallest containing* locus, which makes the attribution single-valued and routes a record enclosed by a larger
  nesting locus to the most specific one — so a read's placement set has at most one entry per record.

---

## 2. The two relations

**Read-confusability $\sim_R$** — the aligner cannot tell the loci apart. With tie tolerance $\Delta$, ceiling
$\mathrm{DE_{max}}$, quorum $k$:
$$
i \sim_R j \iff \big|\{\, r : d(r,i),d(r,j)\le \mathrm{DE_{max}},\ |d(r,i)-d(r,j)|\le \Delta \,\}\big| \ge k .
$$
Because attribution is per-record, $r$ contributes only when it is a genuine *alternative* placement (two distinct
records), never when one locus merely nests in another. Let $R^\* $ be the **read-coupling** relation: $i\,R^\*\,j$
iff $i,j$ lie in the same connected component of $\sim_R$.

**Backbone-homology $\sim_B$** — the loci are the same gene over a real fraction of *both* copies:
$$
i \sim_B j \iff \min\!\big(\mathrm{cov}_i,\mathrm{cov}_j\big)\ge \tau ,
$$
with $\mathrm{cov}_i$ the aligned fraction of $S(i)$ (minimap2 `asm20`). Reciprocal coverage is what distinguishes a
shared backbone (both copies largely homologous) from a shared element (one copy aligns only over a short insert).
A separate identity floor is unnecessary: the aligned region is intrinsically high-identity — every pair clearing
$\tau$ has alignment identity $\ge 0.96$ (measured, $n=393$) — so identity never binds independently and is **not a
free parameter** (§6).

> **Granularity (resolves an apparent code/note discrepancy).** $\sim_B$ is a **copy-level** relation: one alignment
> of the whole exon-union models $S(i),S(j)$, scored by reciprocal coverage. The production family-graph builder
> *additionally* applies a **contiguous-core** test at the **exon level** (per shared exon node). The two are
> complementary, not competing — copy-level coverage certifies "these loci are copies"; exon-level contiguity governs
> which individual exons fuse into one shared graph node. Contiguous-core is **not** substitutable at the copy level:
> on whole exon-union models it fragments at exon boundaries and rejects real copies (measured: it drops ~20 % of clean
> two-copy families that are 99 % identical over 80 % of their length).

---

## 3. The definition

> **Definition.** A **multi-copy gene family** is a connected component, of size $\ge 2$, of the graph
> $G = \big(V,\ \sim_B \cap R^\*\big)$ — i.e. a $\sim_B$-connected component *within* a read-coupled group. Its
> **copies** are the locus models deduplicated by exon-union identity.

Equivalently and operationally: take the connected components of $\sim_R$ (read-coupled candidate families); within
each, keep the connected components of $\sim_B$ (the backbone-coherent cores). A locus that is read-coupled to a group
but shares a backbone with no member of it is in no family.

*Why these coincide (and why $R^\*$, not $\sim_R$).* $R^\*$ is an equivalence relation, so every $\sim_B\cap R^\*$ edge
lies within a single $R^\*$-class and connectivity cannot cross classes; within a class $\sim_B\cap R^\*=\sim_B$. Hence
$\mathrm{components}(V,\sim_B\cap R^\*)=\bigcup_\text{classes}\mathrm{components}(\sim_B|_\text{class})$. The closure is
load-bearing: the analogous statement with the *direct* relation $\sim_R$ is false — in a $>\!2$-copy array whose
distal copies never co-place within one read's placement set, $\sim_R$ would fragment a coupled group that $R^\*$
keeps whole so $\sim_B$ can re-knit it.

---

## 4. Why both relations — and why this is the right object

Each relation alone is a known, deficient definition; the contribution is that the family is their conjunction.

- **$\sim_R$ alone is the read-conflict graph.** Necessary — if two loci are not confusable, assigning a read between
  them is trivial and they are not one unit — but it *over-couples*: a shared transposable element or a processed
  retrocopy makes the reads of two unrelated loci cross-map, fusing them.
- **$\sim_B$ alone is sequence homology** (the "biological" definition). Necessary — copies must share gene structure
  — but it is annotation-style: it admits unexpressed copies (out of RNA scope) and says nothing about whether the
  reads are actually confusable. (It does *not*, on its own, exclude domain-sharers; those are removed by $\sim_R$, P1.)

A family is the **conjunction**: loci that are at once *coupled* ($R^\*$ — the copy-assignment problem is non-trivial)
and *the same gene* ($\sim_B$ — it is well-posed). $\sim_R$ supplies the coupling, $\sim_B$ supplies the identity; the
elegance is intrinsic — the "validation" is not a separate mechanism but **connectivity of $\sim_B$ inside a coupled
group**, so a backbone-less bridge locus is simply not in any component, and a genuine pair is a single
$\sim_B$ edge inside its coupled group (no density required).

---

## 5. Properties

**Structural — hold by construction** (P1 additionally requires tight vertices, P6).

- **P1 (domain-sharers excluded) — by construction *given* tight vertices (P6).** Single-valued per-record best-overlap
  gives nested/adjacent single-domain genes $0$ conflicting reads ⟹ not $\sim_R$ ⟹ no family. The *implication* is
  structural, but its premise — best-overlap is single-valued and faithful — holds only when vertices are tight
  transcription units; under coarse vertices best-overlap can mis-attribute and manufacture domain-sharer edges, so
  **P1 inherits P6's conditionality** (it is structural-given-P6, not unconditional). *Panel: 5/5 domain-sharers
  excluded, though they pass DNA homology.*
- **P3 (no isoform pollution).** Copies are exon-unions, so within-copy alternative splicing is an intra-copy bubble,
  never an extra copy. *100s–1000s of isoforms per locus collapse to one $S(v)$.*
- **P4 (pairs admitted).** A genuine pair is a single $\sim_B$ edge inside its coupled group — no triangle or density
  is required (which is why clique/$k$-truss objects, needing density, are wrong: 57 % of families are pairs).

**Empirical / conditional — hold on the measured substrate.**

- **P2 (repeat/retro bridges rejected).** Two loci whose reads cross-map through a shared element but which are *not
  copies of each other* fail $\sim_B$ (one side aligns only over the short shared insert) ⟹ rejected. *Decisive
  measurement: `OCLN~SEPTIN7` (3,369 confusing reads) has backbone coverage $0.05$, `BCAS4~CCDC30` $0.09$, vs real
  copies `RABL2A~RABL2B` $0.94$, `CCDC196~LOC` $0.97$.* (A retrocopy that genuinely **is** a copy of its parent shares
  a backbone and is correctly kept; $\sim_B$ separates "retro-mediated bridge between non-copies" from "retrocopy of
  parent" exactly by the backbone.)
- **P5 (resolvable paralogs excluded — by quorum).** A diverged copy is excluded iff *fewer than $k$* reads tie at the
  `de` floor; this is a quorum statement, not a per-read one — inclusion can be carried by a quorum even when the
  median per-read gap exceeds $\Delta$ (e.g. RABL2, median gap $0.006$), and `EEF1A1`'s retro-pseudogene is excluded
  because $0$ reads tie. The definition is thus *recent-paralogy ∩ read-confusability*, a strict subset of paralogy.
- **P6 (well-posed vertices).** Reciprocal coverage in $\sim_B$ is meaningful only when $S(v)$ is a tight
  transcription unit; this is conditional on the de-novo locus pipeline (§1) producing per-locus, not coarse, spans —
  the one external assumption.

---

## 6. Parameters

| symbol | role | value | basis |
|---|---|---|---|
| $\Delta$ | $\sim_R$ tie tolerance | 0.005 | single-read divergence resolution at HiFi error: per-read SE $\sqrt{\epsilon/L}\approx0.0009$, tie statistic $\sqrt{2\epsilon/L}\approx0.0013$, $\sim 4\sigma$ |
| $\mathrm{DE_{max}}$ | $\sim_R$ divergence ceiling | 0.05 | loose copy-vs-distinct-gene ceiling |
| $k$ | $\sim_R$ quorum | 3 | **quorum classifier, load-bearing** (not an inert floor): it admits true positives the per-read tie test alone misses — RABL2's median per-read $|\Delta d|=0.0061>\Delta$, so only the quorum of $k$ tied reads carries it (§5, P5) |
| $\tau$ | $\sim_B$ reciprocal coverage | 0.30 | the one data-calibrated threshold: set in the empty gap between repeat-bridges (one-sided coverage $\le 0.1$) and validated copies ($\ge 0.31$ genome-wide). $\tau=0.30$ sits at the gap's lower-copy edge — permissive (it admits partial-coverage copies); it could be centred lower in the gap at a small recall cost |
| $\mathrm{GUARD}$ | $\sim_B$ rejection guard | 20 reads | a backbone-isolated locus is rejected as a bridge only with $\ge\mathrm{GUARD}$ reads (enough to model confidently); below it the locus is held out as unmodelled, not rejected |

An identity floor $\iota$ is **omitted as measured-inert**: every pair clearing $\tau$ already has alignment identity
$\ge 0.96$ ($n=393$, min $0.963$), so identity never binds independently. $\Delta$ is a measurement constant; $k$ is a
load-bearing quorum classifier; $\tau$ is the single empirical threshold, placed in a measured gap. $\mathrm{GUARD}$
gates *rejections* only, so §3/§4's "a backbone-less locus is in no family" is exact for well-expressed loci and a
hold-out (not a rejection) for sparse ones.

---

## 7. Evidence (summary; full record and the rejected-alternatives survey in `family_definition_formal.md`)

- **Panel.** 17 hand-labelled IsoSeq candidates: **TP = 7, TN = 10, FP = 0, FN = 0**.
- **$\sim_B$ is the decisive lever.** Backbone coverage separates repeat/retro bridges (one-sided coverage $\le 0.1$)
  from validated copies — the clean panel cases sit high (RABL2A~RABL2B $0.94$), and genome-wide the lowest validated
  copies reach $\ge 0.31$, so the gap is $[0.1,\,0.31]$ with $\tau=0.30$ at its upper edge. This is the population no
  read-level predicate (coverage, intron, junction concordance) could separate, because that contamination is
  read-indistinguishable from real paralogy.
- **Robust to multimapper sampling (new $-N\,50\,-p\,0.1$ BAM).** The old alignment used minimap2's default secondary
  cap, undersampling $\sim_R$. Re-aligning with $-N\,50\,-p\,0.1$ surfaces $12$–$65\times$ more cross-mapping; $\sim_R$
  then *recovers* hidden dispersed paralogs (per-chrom $\sim_B$-validated copies grow, e.g. chrY $21\to35$) while
  $\sim_B$ prunes every extra bridge (e.g. $10\to59$ on one chromosome). Family identity tracks copy-model homology,
  not read counts (`family_def_newbam_validation.md`).
- **Genome-wide, OOM-safe** (memory > 15 GB free throughout; copy models built only for the ~1,300 family loci,
  reads capped, streamed). Over de-novo loci: **212 candidate → 196 validated families**; cross-chromosome bridges
  (components spanning $\ge 3$ chromosomes) **20 → 12**; **14** well-modelled backbone-less members rejected as
  confident bridges; size-2 families 80 → 73. The genome-wide refinement is *modest* — most de-novo candidate
  families are already coherent, so $\sim_B$ mainly removes the residual cross-chromosome bridges; the *decisive*
  evidence for $\sim_B$ is the per-candidate separation above, not a large genome-wide rejection count.
- **Honest cost.** $\sim_B$ is a precision filter with a small recall cost: among the 14 rejected members the
  rejected edges are $9$ no-homology vs $12$ DNA-paralog/sub-bar, so it is conservative — calibrated to discard
  confident bridges, not to maximise rejection.

---

## 8. Relation to the copy-assignment problem

A family is, by construction, the unit on which read-to-copy assignment is at once necessary ($R^\*$ couples its reads)
and well-defined ($\sim_B$ guarantees a shared backbone, hence comparable copies). The definition hands the
copy-assignment / identifiability machinery exactly its input — mutually confusable, genuinely homologous copies — and
nothing else.

---

## 9. What it does not claim

- It defines families that are **expressed and read-confusable**; unexpressed or read-resolvable paralogs are out of
  scope **by design** (they pose no copy-assignment problem), not missed.
- "Family" here is **read-indistinguishable homologous copies** — a strict subset of sequence paralogy; genuine
  repeats and bridges between non-copies are rejected by $\sim_B$, and read-resolvable diverged paralogs by $\sim_R$.
- **Copy number is deferred.** Counting copies requires PSV-haplotype resolution across the shared backbone and is not
  computed here; the backbone supplies its substrate, and exonically identical copies collapse to one backbone (the
  identifiability floor).
- The single external dependency is the de-novo locus set (§1); the closed-form claims are conditional on it producing
  tight transcription-unit vertices.
