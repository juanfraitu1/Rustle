# Theorem 3 — polynomial recovery under Strong Separation: design

*Adds one theorem to the existing copy-assignment theory note (`bench/copy_assignment_theory.md`), closing the
explicitly-deferred polynomial-recovery remark and completing the hardness/recovery dichotomy into a full
algorithmic picture. Theory only; the "algorithm" is stated and proven, plus an exhaustive verification check —
no production solver.*

## Goal

State and prove **Theorem 3**: under Strong Separation the true copy set is recoverable in `O(n²·m)` time
(n = #reads, m = #columns), not merely *identifiable*. The proof is a direct corollary of the already-proven
structural fact (Strong Separation ⟹ conflict graph `H` is complete K-partite). Back it with an exhaustive
verification check that reuses the Theorem-2 enumeration harness.

## Non-goals

- No production / Rust solver, no integration with `detect_and_assign`. The algorithm is a mathematical object
  + an exhaustive check.
- No claim of optimality of the `O(n²·m)` bound (it is *a* polynomial bound; tighter is not pursued).
- No error-tolerant / noisy extension (that was a separate option, not chosen).

## The theorem

**Setup (unchanged from §5).** Reads `R` (|R| = n) over columns `[m]`; conflict graph `H` (reads adjacent iff
they co-observe a column with differing alleles); Strong Separation = for all i≠j, every read of copy i
conflicts with every read of copy j. Lemma 1: MCC = χ(H). Theorem 2: under Strong Separation, `C*` is the
unique minimum cover, and `H` is exactly the **complete K-partite** graph on the true read-classes.

**Compatibility graph.** `H̄` = the complement of `H` on the read set: `{r, r'} ∈ H̄` iff `r, r'` do **not**
conflict (compatible — they agree at every co-observed column, possibly observing disjoint columns).

> **Theorem 3 (polynomial recovery).** Under Strong Separation, the connected components of the compatibility
> graph `H̄` are exactly the true copy read-classes `C*`; computing them, and the allele-vector of each class
> as the union of its reads' observations, recovers `C*` in `O(n²·m)` time.

**Algorithm `RECOVER(R)`.**
1. Build `H̄`: for every read pair `(r, r')`, add an edge iff compatible (scan co-observed columns). `O(n²·m)`.
2. Connected components of `H̄` (union-find or BFS). `O(n²)`.
3. For each component, the allele-vector = the union of its reads' `(column → allele)` observations (well-defined:
   the component is internally compatible, so no column gets two alleles). `O(n·m)`.
Return the components (read-partition) and their allele-vectors.

**Proof.** By Theorem 2's structural fact, Strong Separation ⟹ `H` is the complete K-partite graph whose parts
are the true read-classes `P_1, …, P_K`. The complement of a complete K-partite graph is the disjoint union of
K cliques, one per part (within a part: no `H`-edges, so all `H̄`-edges present → a clique; across parts: all
`H`-edges, so no `H̄`-edges). Hence `H̄`'s connected components are exactly `P_1, …, P_K`. Each `P_i` is
internally compatible, so its reads' observations agree column-wise and union to a single allele-vector — the
true copy `c_i` (it agrees with every read of `P_i`, which all originate from `c_i`). Step 1 is `O(n²·m)`,
steps 2–3 are dominated by it. ∎

**Self-certifying addendum.** `RECOVER` can verify Strong Separation on its input in the same `O(n²)` by
checking `H̄` is a disjoint union of cliques (each component is complete in `H̄`). If the check fails, the input
does not satisfy Strong Separation — the instance is not certified identifiable — and `RECOVER` reports
*"not strong-separated"* rather than returning a partition that Theorem 3 does not guarantee. (At K=0 / the
recombination regime, `H̄` is not a disjoint-clique union, so the certificate correctly fails.)

## The check

`check_thm3_recovery_algorithm` (added to `bench/copy_assignment_theory_checks.py`):
- Implement `recover(reads)` → (partition, per-component allele-vectors) via compatibility-graph connected
  components, and `is_disjoint_clique_union(reads)` for the self-certifying check.
- Reuse the existing exhaustive enumeration harness (same `(copies, windows)` generator as
  `check_thm2_strong_exhaustive`, same K∈{2,3}, L=3 bound): for **every** instance satisfying Strong Separation
  (`cond_strong`), assert `recover(reads)` returns the **true** partition (0 mismatches over all such instances).
- Also assert: on the explicit recombination witness (K=3, sep+link holds but strong fails), the self-certifying
  check `is_disjoint_clique_union` returns False (RECOVER refuses rather than mis-recovers).
- Deterministic, exhaustive over the bounded space; runs within the suite's existing ~25s envelope.

## Note edits (ripples)

- **New theorem block** placed immediately after Theorem 2 / the Proposition in §5 (numbered Theorem 3). It must
  reference Theorem 2's complete-K-partite fact explicitly (the proof depends on it).
- **§5 deferred-recovery Remark**: change from "polynomial recovery is deferred" to "polynomial recovery is given
  by Theorem 3."
- **§8.4 discussion**: currently presents poly-recovery as an open/deferred question; update to state it is
  resolved by Theorem 3 under Strong Separation (general case remains NP-hard by Theorem 1 — the dichotomy is now
  hardness vs. poly-time-recovery-under-identifiability, fully closed on both axes).
- **§1 intro**: the one sentence about "the deferred algorithmic question" updates to note Theorem 3 closes it
  under Strong Separation.

## Risks / open points

- The proof leans entirely on Theorem 2's "H is complete K-partite under Strong Separation" fact. That fact was
  independently verified in the final review (0 deviations over 808 K=3 strong instances). The new check
  re-exercises it via `recover` agreeing with the true partition, so a regression in that fact would surface.
- `is_disjoint_clique_union` must be correct (a component being connected in `H̄` is not enough — it must be a
  *clique*). The check must test a non-strong instance where some component is connected-but-not-complete to
  confirm the certificate rejects it (the recombination witness provides this).
