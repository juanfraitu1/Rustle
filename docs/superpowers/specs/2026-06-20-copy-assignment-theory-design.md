# Copy assignment as constrained path-cover: a hardness/recovery dichotomy — design

*Advisor interest #2 (assign reads to copies under ambiguity), theory-first deliverable. The output is a
written theory note: a formal problem statement plus two theorems with full proofs, citing the existing
empirical results. No new solver or validation is built in this milestone.*

## Goal

State per-read copy assignment as a clean combinatorial optimization on the family variation graph and prove a
**dichotomy**: the problem is NP-hard in general, but exactly recoverable in polynomial time under an
identifiability condition whose threshold coincides with the measured K-bound. This unifies the three advisor
interests — family detection (#1), copy assignment (#2), and the identifiability theorem — under one condition.

## Non-goals (deferred to a follow-up)

- Building the explicit constrained-flow / LP path-cover solver, or any new code.
- Fresh empirical validation (we cite the existing sim5x K-ladder and GGO silver-standard results).
- The abundance/quantification coupling (per-copy expression as a constraint), polyploid/>2-allele extensions
  beyond what the proofs require, and read-error models richer than the single error-tolerant variant noted below.

---

## The formal problem — Constrained Minimum Copy-Path-Cover (CMCPC)

**Family variation graph.** `G = (V, E)` is a DAG; vertices are exonic segments, edges are splice junctions /
adjacencies. A distinguished set of **PSV columns** `Φ = {1,…,L}` are positions carrying alleles, each with an
alphabet `Σ_j` (`|Σ_j| ≥ 2`).

**Copy.** A copy is a pair `(π, a)`: a path `π` through `G` (its isoform structure) and an allele-vector
`a ∈ ∏_{j∈Φ(π)} Σ_j` over the PSV columns the path crosses — a haplotype *with* isoform structure.

**Read.** A read `r` is a subpath of `G` observing a subset `S_r ⊆ Φ` of columns with alleles `o_r : S_r → Σ`,
and a subset of junctions. A read originates from exactly one copy and is **consistent** with copy `(π,a)` iff
`o_r(j) = a(j)` for all `j ∈ S_r ∩ Φ(π)` and its junctions lie on `π`.

**CMCPC.** Find a minimum-cardinality set of copies `P = {(π_1,a_1),…,(π_k,a_k)}` and an assignment
`α : R → P` with every read consistent with `α(r)`. Equivalently, a **minimum flow-decomposition** of the
read-flow on `G` into paths, **constrained** by allele-linkage (reads sharing a PSV column route together).

> **The framing's spine.** Unconstrained minimum path-cover / flow-decomposition on a DAG is polynomial
> (Dilworth; max-flow). The **allele-linkage constraint** makes it NP-hard. **Identifiability** makes it
> polynomial again. The dichotomy is exactly this constraint switching on and off.

**The allele-vector core.** Per design decision (b), the theorems are proved on the core obtained by collapsing
isoform structure to a single shared backbone: copies are allele-vectors `c ∈ Σ^L`, reads are partial
observations `(S_r, o_r)`, and **Minimum Copy Cover (MCC)** asks for the minimum number of allele-vectors (and
a consistent assignment) explaining all reads. Path/isoform structure re-enters as a corollary (below).

---

## Theorem 1 — Hardness

> **MCC is NP-hard.** *(Hence CMCPC, which contains MCC as the single-backbone case, is NP-hard.)*

**Vehicle: reduction from graph coloring; framing: the K-copy generalization of MEC.** Define the read
**conflict graph** `H`: reads `r, r'` are adjacent iff `∃ j ∈ S_r ∩ S_{r'}` with `o_r(j) ≠ o_{r'}(j)`. Two
reads can share a copy iff they pairwise agree at every co-observed column (an independent set in `H`); since
consistency is per-column and per-column agreement is pairwise, **pairwise non-conflict ⇒ joint consistency**,
so a copy is exactly an independent set of `H` and `\mathrm{MCC} = χ(H)`, the chromatic number. Any graph `H`
is realizable as a conflict graph (one PSV column per edge `(u,v)`: `u` observes allele 0, `v` observes allele
1, others abstain), so **graph coloring ≤ₚ MCC**.

**MEC framing.** The error-tolerant variant — fixed `k` copies, minimize total allele mismatches (flips) — is
the **K-copy, variation-graph generalization of Minimum Error Correction** (haplotype assembly), which is
NP-hard already at `k = 2`. Thus hardness is robust to both the model (error-free min-#copies via coloring) and
the genomics-standard formulation (error-tolerant fixed-k via MEC).

**Proof obligations (for the plan):** (i) the conflict-graph realizability gadget, with `|Σ_j|` bounded; (ii)
the `MCC = χ(H)` equivalence including the pairwise⇒joint argument; (iii) a clean statement relating the
error-tolerant variant to MEC (citation + the generalization direction).

---

## Theorem 2 — Exact recovery under identifiability

> **Let `C* = {c_1,…,c_K}` be the true copies. Under identifiability condition `𝒞`, the minimum copy cover is
> unique and equals `C*`, and is computable in polynomial time.**

**Condition `𝒞`.**
- **(C1) Pairwise distinguishability.** Every pair `c_i, c_j` differs at ≥ 1 PSV column (the per-pair PSV count
  `K_{ij} ≥ 1`; the robust regime is `K_{ij} ≥ 2`).
- **(C2) Linkage / coverage.** The reads connect the distinguishing columns: for every pair of true copies, the
  reads drawn from them co-observe a distinguishing column densely enough that the **agreement graph** (complement
  of `H` restricted to true-copy membership) has the `c_i` as its connected components, and forcing two copies'
  reads into one color class necessarily creates a conflict. (Error-tolerant form: per-column per-copy coverage
  `≥ m` so the majority allele is correct w.h.p.)

**Proof strategy.** Under `𝒞`, `H` is a **cluster graph** — a disjoint union of `K` cliques, one per true copy:
same-copy reads pairwise agree (non-adjacent), and (C2) forces some cross-copy conflict for every pair, so no
independent set spans two copies. On a cluster graph, minimum coloring is polynomial and unique, with
`χ(H) = K`; recovering the colors = the connected components of the agreement graph = `C*`. This is the
K-copy, variation-graph generalization of the classical result that **error-free haplotype phasing is
polynomial when reads connect all heterozygous sites with distinguishing overlaps**.

**The K-bound, recovered as a theorem.** `𝒞`'s distinguishability axiom (C1) is exactly the K-bound:
- `K_{ij} = 0` (copies identical over the transcribed region) ⇒ (C1) fails ⇒ no cross-copy conflict ⇒ minimum
  cover *merges* `c_i, c_j` (fewer copies) ⇒ true copies **provably unrecoverable**. This is the MAGEA
  co-located regime (resolvable fraction 0/494) — RNA-level non-identifiability as a hardness corollary.
- `K_{ij} ≥ 2` with linkage ⇒ `𝒞` holds ⇒ exact recovery. This is the sim5x K-ladder regime (K≥2 → 100%).

**Proof obligations (for the plan):** (i) the precise, minimal form of (C2) that makes `H` a cluster graph —
the central rigor risk; (ii) uniqueness (not just optimality) of the minimum cover under `𝒞`; (iii) the
polynomial recovery procedure and its correctness; (iv) the error-tolerant restatement with coverage `m` and a
w.h.p. guarantee; (v) explicit statement that `K_{ij}=0 ⇒` unrecoverable.

---

## Corollary — paths / isoforms (joint copies + isoforms)

Re-attach isoform structure: a copy is `(π, a)`, and reads carry junction observations. The corollary lifts
the allele-vector recovery to **paths**: under `𝒞` augmented with junction-linkage (reads connect the
junctions distinguishing isoforms, the same condition applied to the path alphabet), the minimum **constrained
path-cover** is unique and equals the true copies-with-isoforms, in polynomial time — recovering #copies *and*
their isoform structure jointly. Plain min-path-cover on the DAG is the unconstrained relaxation; the corollary
shows allele+junction linkage promotes it from "a minimum cover" to "*the* true cover" under `𝒞`.

**Proof obligation:** show the allele-vector argument composes with the DAG path-cover (flow) so the joint
object inherits uniqueness/recovery — or identify the extra condition junctions require.

---

## Empirical evidence (cited, not produced here)

- **THM 2 / `K≥2` regime:** sim5x K-ladder (`bench/sim_reads.py`, `/home/juanfra/winloci_scratch/sim5x/`),
  K≥2 → 100% correct assignment; GGO silver-standard (resolution vs uniquely-mapped reads) = 100%
  (`detect_and_assign` smoke run, MAGEA region: 1026/1026).
- **THM 2 boundary / `K=0`:** MAGEA co-located arrays, resolvable fraction 0/494; copies sequence-identical over
  the transcribed exon (`NM_A == NM_B`). (`bench/resolution_improvement_bound.md`.)
- **The same condition governs detection (#1):** the family is the conflict-graph component; `bench/family_definition_formal.md`.

---

## Document structure (the theory note to be written)

1. Introduction — copy assignment as the unit downstream of family detection; the multimapping-resolution lineage.
2. The family variation graph and CMCPC (definitions).
3. The allele-vector core (MCC).
4. Theorem 1 (hardness) — coloring reduction; MEC framing; proof.
5. Theorem 2 (recovery) — condition `𝒞`; cluster-graph proof; the K-bound corollary.
6. Corollary — paths/isoforms (joint inference).
7. Empirical corroboration (cited).
8. Discussion — the dichotomy as the through-line across interests #1–#3; the K=0 frontier (where interest #3,
   allele-specific junctions, or longer reads become the lever).

## Risks / open rigor points

- **(C2) is the crux.** The minimal correct linkage/coverage condition making `H` a cluster graph must be stated
  precisely; too weak ⇒ spurious merges, too strong ⇒ vacuous. Mitigation: anchor to the established error-free
  phasing condition and generalize.
- **Min-#copies vs error-tolerant objective.** Coloring (min-#copies) and MEC (min-flips, fixed k) are different
  optimizations; the note must keep them distinct and state which theorem governs which.
- **Paths corollary** may need a condition beyond the allele-vector core (junction linkage); if it does not
  compose cleanly, fall back to stating the core result and treating paths as future work (per the user's
  "interject if a better option appears").
- **`|Σ_j|` in the reduction** — keep the alphabet small (binary suffices for coloring) so hardness is not an
  artifact of large alphabets.
