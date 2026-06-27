# Theorem 3 (Polynomial Recovery) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add Theorem 3 (under Strong Separation the true copies are recoverable in O(n²·m) time) to the copy-assignment theory note, backed by an exhaustive verification check, closing the deferred-recovery remark.

**Architecture:** Task 1 adds the algorithm (`recover`) + self-certifying predicate (`is_disjoint_clique_union`) + an exhaustive check to `bench/copy_assignment_theory_checks.py`, TDD'd via the existing Theorem-2 enumeration harness. Task 2 writes the Theorem 3 block + proof into `bench/copy_assignment_theory.md` and updates the three "recovery deferred" spots. The proof is a corollary of an already-proven fact (Strong Separation ⟹ conflict graph H is complete K-partite), so the executable check is the load-bearing verification.

**Tech Stack:** Markdown + LaTeX math; Python 3 (stdlib `itertools` + the existing helpers in the checks file; `networkx` already imported there).

## Global Constraints

- Theory only — no Rust / production solver, no `detect_and_assign` integration.
- Reuse the existing checks-file helpers; do not duplicate them. Reads are `frozenset` of `(column, allele)` pairs throughout.
- The algorithm operates ONLY on `reads` (it does not see the true labels/copies); the check compares its output to the true partition.
- The exhaustive check reuses the SAME enumeration bound as `check_thm2_strong_exhaustive` (K∈{2,3}, L=3, exactly 2 windows/copy; K=2 over all windows, K=3 over size≥2 windows). Deterministic.
- `python3 bench/copy_assignment_theory_checks.py` must exit 0 with all checks OK.
- Theorem 3's proof MUST cite Theorem 2's "H is complete K-partite under Strong Separation" fact (the proof depends on it).

## Before you start (both tasks)

Read `bench/copy_assignment_theory_checks.py` and `bench/copy_assignment_theory.md` first. In the checks file, note the EXACT signatures of the existing helpers you will reuse — run:
```
grep -n "^def observed\|^def partition_of\|^def cond_strong\|^def conflict_graph\|def check_thm2_strong_exhaustive" bench/copy_assignment_theory_checks.py
```
`cond_strong`'s parameter order matters — use whatever the file actually defines (do not assume). `observed(read)` returns `{column: allele}`. `partition_of(labels)` returns the canonical `frozenset(frozenset(read-indices))` partition.

---

## Task 1: `recover` + `is_disjoint_clique_union` + exhaustive check

**Files:**
- Modify: `bench/copy_assignment_theory_checks.py`

**Interfaces:**
- Consumes: `observed`, `partition_of`, `cond_strong`, and the enumeration generator pattern from `check_thm2_strong_exhaustive`.
- Produces: `compatibility_edges(reads) -> list[(i,j)]`, `recover(reads) -> frozenset[frozenset[int]]`, `is_disjoint_clique_union(reads) -> bool`, `check_thm3_recovery_algorithm()`.

- [ ] **Step 1: Write the check first (TDD), plus the explicit recombination-witness assertion**

Add to `bench/copy_assignment_theory_checks.py` (after the Theorem-2 helpers/checks). Use the existing module helpers; do not redefine `observed`/`partition_of`/`cond_strong`.

```python
def compatibility_edges(reads):
    """Edges of the compatibility graph H-bar (complement of the conflict graph): the pairs of reads that
    do NOT conflict (they agree at every co-observed column)."""
    edges = []
    for i, j in itertools.combinations(range(len(reads)), 2):
        oi, oj = observed(reads[i]), observed(reads[j])
        if all(oi[c] == oj[c] for c in (oi.keys() & oj.keys())):
            edges.append((i, j))
    return edges


def recover(reads):
    """RECOVER (Theorem 3): connected components of the compatibility graph H-bar, via union-find.
    Returns the read-partition as a frozenset of frozensets of read indices. O(n^2 * m)."""
    n = len(reads)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, j in compatibility_edges(reads):
        parent[find(i)] = find(j)
    comp = {}
    for i in range(n):
        comp.setdefault(find(i), set()).add(i)
    return frozenset(frozenset(c) for c in comp.values())


def is_disjoint_clique_union(reads):
    """Self-certifying check (Theorem 3 addendum): is the compatibility graph a disjoint union of cliques?
    Equivalently, is every connected component of H-bar internally complete (all pairs compatible)? This
    holds iff the input satisfies Strong Separation. O(n^2)."""
    n = len(reads)
    compat = [[i == j for j in range(n)] for i in range(n)]
    for i, j in compatibility_edges(reads):
        compat[i][j] = compat[j][i] = True
    for cls in recover(reads):
        members = sorted(cls)
        for a in range(len(members)):
            for b in range(a + 1, len(members)):
                if not compat[members[a]][members[b]]:
                    return False  # same component but not directly compatible -> not a clique
    return True


def check_thm3_recovery_algorithm():
    """EXHAUSTIVE certificate of Theorem 3: over every Strong-Separation instance in the K in {2,3}, L=3
    enumeration, RECOVER returns exactly the TRUE partition, and is_disjoint_clique_union accepts it.
    Also: the explicit recombination witness (sep+link holds, strong fails) is REJECTED by the certificate."""
    L = 3
    cols = list(range(L))
    all_windows = [frozenset(s) for k in range(1, L + 1) for s in itertools.combinations(cols, k)]
    windows_for = {2: all_windows, 3: [w for w in all_windows if len(w) >= 2]}
    copyvecs = [tuple(c) for c in itertools.product((0, 1), repeat=L)]

    n_strong = 0
    recover_mismatch = 0
    cert_reject_on_strong = 0
    for K in (2, 3):
        windows = windows_for[K]
        for copies in itertools.combinations(copyvecs, K):
            copies = list(copies)
            for assign in itertools.product(itertools.combinations(windows, 2), repeat=K):
                reads, labels = [], []
                for ci, wins in enumerate(assign):
                    for w in wins:
                        reads.append(frozenset((c, copies[ci][c]) for c in w))
                        labels.append(ci)
                if not cond_strong(copies, reads, labels):   # use the file's actual cond_strong signature
                    continue
                n_strong += 1
                if recover(reads) != partition_of(labels):
                    recover_mismatch += 1
                if not is_disjoint_clique_union(reads):
                    cert_reject_on_strong += 1

    assert recover_mismatch == 0, f"RECOVER disagreed with the true partition on {recover_mismatch} strong instances"
    assert cert_reject_on_strong == 0, f"self-certify wrongly rejected {cert_reject_on_strong} strong instances"

    # Recombination witness: sep+link holds, strong fails -> certificate must REJECT (refuse to recover).
    wc = [(1, 1, 0), (0, 0, 1), (0, 1, 1)]
    wwins = [(0, {1, 2}), (0, {0, 1}), (1, {0, 1, 2}), (1, {0, 1}), (2, {0, 1}), (2, {1, 2})]
    wreads = [frozenset((c, wc[ci][c]) for c in w) for ci, w in wwins]
    assert not is_disjoint_clique_union(wreads), "certificate must reject the non-strong recombination witness"

    print(f"    [thm3] strong instances={n_strong}: RECOVER==true 100%, self-certify accepts all; witness rejected")
    return f"Thm 3: RECOVER == true partition on all {n_strong} strong instances (exhaustive K=2,3/L=3); witness rejected"


CHECKS.append(check_thm3_recovery_algorithm)
```

(If `cond_strong` in the file takes a different argument order, adjust the call accordingly — the behavior is "every cross-copy read pair conflicts.")

- [ ] **Step 2: Run the suite; the new check must fail first only if helpers are missing, else pass**

Run: `python3 bench/copy_assignment_theory_checks.py 2>&1 | tail -8`
Expected: all checks `OK` including `Thm 3: RECOVER == true partition on all N strong instances ...`. If `recover_mismatch > 0`, STOP — that would contradict Theorem 2's structural fact and means a real bug in `recover` or `compatibility_edges` (debug those, do not weaken the assertion). Confirm exit 0.

- [ ] **Step 3: Sanity-check the time**

Run: `time python3 bench/copy_assignment_theory_checks.py >/dev/null 2>&1`
Expected: completes in roughly the existing envelope (the Thm-3 enumeration skips the expensive `_min_covers`, so it adds little). If it exceeds ~60s, note it in the report but do not change the bound.

- [ ] **Step 4: Commit**

```bash
git add bench/copy_assignment_theory_checks.py
git commit -m "theory(copy-assign): Theorem 3 recovery algorithm + self-certify + exhaustive check"
```

---

## Task 2: Theorem 3 prose + proof + ripple updates

**Files:**
- Modify: `bench/copy_assignment_theory.md`

**Interfaces:**
- Consumes: the executable check from Task 1 (referenced by name in the prose).
- Produces: the §5 Theorem 3 block; updates to the §5 recovery Remark, §8.4 discussion, and §1 intro.

- [ ] **Step 1: Add the Theorem 3 block immediately after the Theorem 2 Proposition in §5**

Find the end of the §5 Proposition / before the "(polynomial-time recovery, deferred)" Remark. Insert:

```markdown
### Theorem 3 (polynomial-time recovery)

Theorem 2 identifies the optimum with the truth; Theorem 3 shows the truth is also *computable* in polynomial
time under the same condition — so the hardness of the general problem (Theorem 1) gives way, under Strong
Separation, not merely to a unique optimum but to an efficient algorithm.

Let the **compatibility graph** $\bar H$ be the complement of the conflict graph on the read set: $\{r, r'\}$
is an edge of $\bar H$ iff $r$ and $r'$ do **not** conflict (they agree at every co-observed column).

> **Theorem 3.** Under Strong Separation, the connected components of $\bar H$ are exactly the true copy
> read-classes $P_1, \ldots, P_K$. The algorithm `RECOVER` — build $\bar H$, take its connected components, and
> set each component's allele-vector to the union of its reads' observations — returns the true copy set $C^\*$
> in $O(n^2 m)$ time ($n = |R|$, $m$ the number of columns).

*Proof.* By Theorem 2 (and its structural core), Strong Separation makes $H$ the **complete $K$-partite** graph
whose parts are the true read-classes $P_1, \ldots, P_K$. The complement of a complete $K$-partite graph is the
**disjoint union of $K$ cliques**, one per part: within a part $P_i$ there are no $H$-edges, so every pair is a
$\bar H$-edge (a clique); across parts every pair is an $H$-edge, so no $\bar H$-edge. Hence the connected
components of $\bar H$ are exactly $P_1, \ldots, P_K$. Each $P_i$ is internally compatible, so its reads agree
at every column any of them observes; the union of their observations is therefore a single allele-vector,
consistent with every read of $P_i$ — i.e. the true copy $c_i$. Building $\bar H$ tests each of the
$\binom{n}{2}$ read pairs over at most $m$ co-observed columns ($O(n^2 m)$); connected components and the
per-component union are dominated by this. $\square$

**Self-certifying recovery.** `RECOVER` can verify Strong Separation on its own input in the same $O(n^2)$ by
checking that $\bar H$ is a disjoint union of cliques — every connected component is internally complete. When
this fails the input is *not* Strong-Separated (e.g. the $K = 0$ regime, or the $K \geq 3$ recombination witness,
where $\bar H$ has a connected-but-incomplete component), and `RECOVER` reports *not identifiable* rather than
returning a partition Theorem 3 does not guarantee. The algorithm thus knows exactly when it is entitled to its
answer.

**Verification.** `check_thm3_recovery_algorithm` enumerates every Strong-Separation instance over $K \in \{2,3\}$,
$L = 3$ and confirms `RECOVER` returns the true partition on **all** of them, the self-certifier accepts each, and
the recombination witness (Strong Separation fails) is **rejected** — exhaustive, not sampled.
```

- [ ] **Step 2: Update the §5 deferred-recovery Remark**

Find the Remark that defers polynomial recovery (it says recovery / the polynomial-time algorithm is out of scope / deferred). Replace its deferral wording so it points to Theorem 3, e.g. change the sentence asserting deferral to:

```markdown
**Remark (recovery).** Polynomial-time recovery under Strong Separation is given by **Theorem 3** below
(connected components of the compatibility graph). The general problem remains NP-hard (Theorem 1); the
dichotomy is therefore fully closed on both axes — hard in general, efficiently solvable under the
identifiability condition.
```
(If the Remark currently sits *before* the Theorem 3 block you inserted, either move this Remark to just after Theorem 3 or reword "Theorem 3 below" to "Theorem 3 above" so the reference direction is correct. Ensure the reference resolves.)

- [ ] **Step 3: Update §8.4 (discussion of the deferred algorithm)**

Find §8.4 (the discussion paragraph presenting polynomial recovery as open/deferred). Update it to state the question is resolved by Theorem 3: under Strong Separation, `RECOVER` (compatibility-graph components) computes the true copies in $O(n^2 m)$; the general case stays NP-hard (Theorem 1). Keep the honesty that Strong Separation is conservative-sufficient (the empirical resolver may succeed more broadly) — Theorem 3 is about the *certified* regime.

- [ ] **Step 4: Update the §1 intro sentence about the deferred algorithmic question**

Find the §1 sentence mentioning "the deferred algorithmic question (the polynomial recovery solver)" (or similar). Update it so the dichotomy reads: general MCC is NP-hard (Thm 1); under Strong Separation the truth is the unique optimum (Thm 2) **and** is recovered in polynomial time (Thm 3). Remove the "deferred" framing from the intro.

- [ ] **Step 5: Run the suite (markdown-only edits — safety check) and read the section**

Run: `python3 bench/copy_assignment_theory_checks.py 2>&1 | tail -3`
Expected: exit 0, all checks `OK` (prose edits don't touch checks). Read §5 end-to-end: Theorem 3's statement, the proof's dependence on the complete-K-partite fact, the self-certifying addendum, and confirm no remaining text calls polynomial recovery "deferred"/"out of scope" anywhere (§1, §5 Remark, §8.4 all updated).

- [ ] **Step 6: Commit**

```bash
git add bench/copy_assignment_theory.md
git commit -m "theory(copy-assign): Theorem 3 (poly recovery under Strong Separation) prose + close deferred-recovery in §1/§5/§8"
```

---

## Self-Review

**Spec coverage:** Theorem 3 statement + algorithm + O(n²·m) proof → Task 2 Step 1; self-certifying addendum → Task 1 (`is_disjoint_clique_union`) + Task 2 Step 1; exhaustive check (recover==true over all strong instances + witness rejected) → Task 1; the three ripple updates (§5 Remark, §8.4, §1) → Task 2 Steps 2–4. All spec sections covered.

**Placeholder scan:** none — full code and full proof text are inline. The only conditional instruction (`cond_strong` argument order) is explicitly delegated to "use the file's actual signature," which the implementer resolves by reading the file (not a placeholder in the deliverable).

**Type/name consistency:** `recover`, `compatibility_edges`, `is_disjoint_clique_union`, `check_thm3_recovery_algorithm` are defined in Task 1 and referenced by those exact names in Task 2's prose and in `CHECKS.append`. Reads are `frozenset((column, allele))` consistently. `partition_of(labels)` and `cond_strong(...)` are reused from the existing file with their actual signatures.
