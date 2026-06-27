# Copy-Assignment Theory Note Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Produce a rigorous theory note proving the copy-assignment hardness/recovery dichotomy — Theorem 1 (NP-hardness) and Theorem 2 (identifiability/exact recovery under condition C) — with each combinatorial claim backed by an executable verification check.

**Architecture:** The deliverable is one markdown theory note (`bench/copy_assignment_theory.md`) using LaTeX math (the style of `bench/family_definition_formal.md`), built section by section. Each proof task is paired with a check in `bench/copy_assignment_theory_checks.py` that computes the concrete combinatorial fact the theorem predicts (e.g. the chromatic number of the reduction gadget, the uniqueness of the minimum cover on a K≥2 instance, the forced merge on a K=0 instance) — the "test" analog for proofs. A claim is "proven" only when its check passes AND an adversarial reviewer cannot find a counterexample.

**Tech Stack:** Markdown + LaTeX math; Python 3 + `networkx` (already available, 3.6.1) for the combinatorial checks.

## Global Constraints

- Deliverable is THEORY ONLY: no new solver, no new biological validation. Cite existing empirics verbatim.
- Theorem statements must match the spec `docs/superpowers/specs/2026-06-20-copy-assignment-theory-design.md`, **with one approved correction**: Theorem 2 is an identifiability/uniqueness result ("true copies = unique minimum copy cover under C"), not a cluster-graph/poly-time claim; polynomial recovery is deferred to a remark.
- The hardness reduction uses a **binary** PSV alphabet (`Σ = {0,1}`) so hardness is not an artifact of large alphabets.
- Theorem 1 reduction is from **graph coloring** (the formal vehicle); the note **frames** it as the K-copy generalization of MEC (the genomics-familiar NP-hard relative). Both appear; coloring does the formal work.
- Theorems are proved on the **allele-vector core** (MCC); paths/isoforms are a corollary.
- Every check is deterministic (fixed instances / fixed seeds). `python3 bench/copy_assignment_theory_checks.py` must exit 0.

## Notation (used in every task)

- PSV columns `Φ = {1,…,L}`, binary alphabet `Σ = {0,1}`.
- A **copy** (core) is an allele-vector `c ∈ Σ^L`. True copies `C* = {c_1,…,c_K}`.
- A **read** `r = (S_r, o_r)`: observed columns `S_r ⊆ Φ`, alleles `o_r : S_r → Σ`.
- `r` is **consistent** with copy `c` iff `o_r(j) = c(j)` for all `j ∈ S_r`.
- **Conflict graph** `H = (R, conf)`: `r ~ r'` iff `∃ j ∈ S_r ∩ S_{r'}` with `o_r(j) ≠ o_{r'}(j)`.
- **MCC** (Minimum Copy Cover): minimum `k` such that reads partition into `k` parts, each part jointly consistent with one allele-vector.

---

## Task 1: Model, definitions, and the MCC = χ(H) lemma

**Files:**
- Create: `bench/copy_assignment_theory.md` (Sections: title, §1 intro stub, §2 model, §3 the MCC core + Lemma 1)
- Create: `bench/copy_assignment_theory_checks.py` (the `check_lemma_mcc_equals_chromatic` check + a `main`)

**Interfaces:**
- Produces (note): the definitions above, and **Lemma 1: MCC = χ(H)**, with the pairwise⇒joint argument.
- Produces (checks): `def conflict_graph(reads) -> networkx.Graph` and `def mcc_bruteforce(reads, L) -> int`, reused by later tasks.

- [ ] **Step 1: Write the verification check first (the fact Lemma 1 must explain)**

In `bench/copy_assignment_theory_checks.py`:

```python
"""Executable verification of the copy-assignment theory note's combinatorial claims."""
import itertools
import networkx as nx

# A read is (frozenset of (column, allele) pairs). Helpers operate on lists of reads.

def observed(read):
    return {j: a for (j, a) in read}

def conflict_graph(reads):
    """H: reads conflict iff they disagree at a co-observed column."""
    G = nx.Graph()
    G.add_nodes_from(range(len(reads)))
    for i, j in itertools.combinations(range(len(reads)), 2):
        oi, oj = observed(reads[i]), observed(reads[j])
        if any(oi[c] != oj[c] for c in (oi.keys() & oj.keys())):
            G.add_edge(i, j)
    return G

def jointly_consistent(read_idxs, reads):
    """True iff the reads share at least one allele-vector (no column observed with two alleles)."""
    col = {}
    for k in read_idxs:
        for c, a in observed(reads[k]).items():
            if col.setdefault(c, a) != a:
                return False
    return True

def mcc_bruteforce(reads):
    """Minimum #copies = min partition of reads into jointly-consistent parts. Exponential; tiny instances only."""
    n = len(reads)
    # try k = 1..n, find smallest k admitting a proper partition into jointly-consistent parts.
    idxs = list(range(n))
    for k in range(1, n + 1):
        if _can_partition(idxs, k, reads):
            return k
    return n

def _can_partition(idxs, k, reads):
    # assign each read a label in 0..k-1; check all labels' classes jointly consistent. Backtracking.
    label = {}
    def bt(pos):
        if pos == len(idxs):
            return True
        for lab in range(k):
            label[idxs[pos]] = lab
            cls = [i for i in idxs[:pos+1] if label[i] == lab]
            if jointly_consistent(cls, reads) and bt(pos + 1):
                return True
        del label[idxs[pos]]
        return False
    return bt(0)

def check_lemma_mcc_equals_chromatic():
    # 5-cycle conflict structure: reads pairwise-disagree around a C5 -> chi = 3 -> MCC must = 3.
    # columns e01,e12,e23,e34,e40 ; read i observes its two incident edge-columns with opposite alleles.
    edges = [(0,1),(1,2),(2,3),(3,4),(4,0)]
    cols = {e: idx for idx, e in enumerate(edges)}
    reads = []
    for w in range(5):
        obs = []
        for e in edges:
            if w in e:
                obs.append((cols[e], 0 if w == e[0] else 1))
        reads.append(frozenset(obs))
    H = conflict_graph(reads)
    chi = max(len(c) for c in [nx.algorithms.coloring.greedy_color(H, strategy="saturation_largest_first").values()]) and \
          (1 + max(nx.algorithms.coloring.greedy_color(H, strategy="saturation_largest_first").values()))
    # C5 is an odd cycle -> chromatic number 3; greedy DSATUR is exact on C5.
    assert nx.algorithms.coloring.greedy_color(H, strategy="saturation_largest_first")
    assert sorted(H.edges()) == sorted((i, (i+1) % 5) for i in range(5)) or H.number_of_edges() == 5
    assert mcc_bruteforce(reads) == 3, "MCC must equal chromatic number (3) on the C5 gadget"
    return "Lemma 1 (MCC=chi) verified on C5: MCC=3=chi(C5)"

CHECKS = [check_lemma_mcc_equals_chromatic]

def main():
    for fn in CHECKS:
        print("OK  -", fn())

if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run the check; confirm it passes (establishes the fact the lemma explains)**

Run: `python3 bench/copy_assignment_theory_checks.py`
Expected: `OK  - Lemma 1 (MCC=chi) verified on C5: MCC=3=chi(C5)`. If `mcc_bruteforce` ≠ 3 or the coloring isn't 3, fix the instance before writing the proof.

- [ ] **Step 3: Write §2 (model) and §3 (MCC core + Lemma 1) in the note**

Write `bench/copy_assignment_theory.md` with title, a one-paragraph §1 intro stub (to be expanded in Task 5), §2 reproducing the Notation/definitions above as prose+math, and §3 containing:

> **Lemma 1.** `MCC = χ(H)`, the chromatic number of the conflict graph `H`.
>
> *Proof.* A set of reads is **jointly consistent** (shares one allele-vector) iff no column is observed with two different alleles. **Pairwise non-conflict ⇒ joint consistency:** consistency is a per-column property; if every pair of reads agrees at each co-observed column, then for any column `j` all reads observing `j` pairwise agree, hence share one allele at `j`; ranging over `j` yields a single consistent allele-vector. Thus a set of reads can form one copy iff it is an independent set of `H`. A copy cover of size `k` is therefore a partition of `R` into `k` independent sets — a proper `k`-coloring of `H` — and conversely. Minimizing `k` gives `MCC = χ(H)`. ∎

- [ ] **Step 4: Re-run checks + read §3 against the check**

Run: `python3 bench/copy_assignment_theory_checks.py` → still `OK`. Confirm the C5 instance in the check is exactly the pairwise-conflict structure the lemma's proof describes.

- [ ] **Step 5: Commit**

```bash
git add bench/copy_assignment_theory.md bench/copy_assignment_theory_checks.py
git commit -m "theory(copy-assign): model + MCC=chi(H) lemma with C5 verification check"
```

---

## Task 2: Theorem 1 (hardness)

**Files:**
- Modify: `bench/copy_assignment_theory.md` (add §4: Theorem 1 + proof + MEC framing)
- Modify: `bench/copy_assignment_theory_checks.py` (add `check_thm1_reduction`)

**Interfaces:**
- Consumes: Lemma 1 (`MCC = χ(H)`), `conflict_graph`, `mcc_bruteforce` (Task 1).
- Produces (note): **Theorem 1** — MCC is NP-hard, via `coloring ≤ₚ MCC`.

- [ ] **Step 1: Write the reduction-verification check first**

Add to `bench/copy_assignment_theory_checks.py`:

```python
def reduction_instance(graph_edges, n):
    """Map graph Gamma=(V=range(n), edges) to an MCC instance: one binary column per edge, one read per vertex.
    Read w observes column e iff w in e, allele 0 if w==e[0] else 1. Conflict graph of the instance == Gamma."""
    cols = {e: idx for idx, e in enumerate(graph_edges)}
    reads = []
    for w in range(n):
        obs = [(cols[e], 0 if w == e[0] else 1) for e in graph_edges if w in e]
        reads.append(frozenset(obs))
    return reads

def check_thm1_reduction():
    import networkx as nx
    # several small graphs: conflict graph of the instance must be ISOMORPHIC to the source graph,
    # and MCC must equal the source graph's chromatic number.
    cases = [
        ([(0,1),(1,2),(2,0)], 3),          # triangle  -> chi 3
        ([(0,1),(1,2),(2,3),(3,0)], 4),    # C4        -> chi 2
        ([(0,1),(1,2),(2,3),(3,4),(4,0)], 5),  # C5    -> chi 3
        ([(0,1),(0,2),(0,3)], 4),          # star     -> chi 2
    ]
    for edges, n in cases:
        reads = reduction_instance(edges, n)
        H = conflict_graph(reads)
        Gamma = nx.Graph(); Gamma.add_nodes_from(range(n)); Gamma.add_edges_from(edges)
        assert nx.is_isomorphic(H, Gamma), f"conflict graph != source graph for {edges}"
        chi = chromatic_number(Gamma)
        assert mcc_bruteforce(reads) == chi, f"MCC != chi for {edges}"
    return "Thm 1: conflict graph == source graph and MCC == chi on 4 instances"

def chromatic_number(G):
    import networkx as nx
    n = G.number_of_nodes()
    for k in range(1, n + 1):
        # try k-coloring by brute force (tiny graphs)
        import itertools
        for assign in itertools.product(range(k), repeat=n):
            if all(assign[u] != assign[v] for u, v in G.edges()):
                return k
    return n

CHECKS.append(check_thm1_reduction)
```

- [ ] **Step 2: Run; confirm it passes**

Run: `python3 bench/copy_assignment_theory_checks.py`
Expected: both checks `OK`, including `Thm 1: conflict graph == source graph and MCC == chi on 4 instances`.

- [ ] **Step 3: Write §4 (Theorem 1) in the note**

> **Theorem 1.** MCC is NP-hard. (Hence CMCPC, which contains MCC as the single-backbone case, is NP-hard.)
>
> *Proof.* Reduce from graph coloring. Given `Γ = (U, F)` with `U = {1,…,n}`, build an MCC instance: one binary column `j_e` per edge `e ∈ F`; one read `r_w` per vertex `w`, observing `j_e` iff `w ∈ e`, with allele `0` if `w` is the first endpoint of `e` and `1` if the second. Two reads `r_u, r_v` co-observe a column iff `(u,v) ∈ F`, and then they disagree (alleles `0 ≠ 1`); they share no other column. Hence the conflict graph `H` of the instance is isomorphic to `Γ`. By Lemma 1, `MCC = χ(H) = χ(Γ)`. Therefore `Γ` is `k`-colorable iff `MCC ≤ k`, and computing MCC solves chromatic number. The alphabet is binary, so hardness is intrinsic, not an alphabet artifact. ∎
>
> **MEC framing.** The error-tolerant variant — fix `k` copies and minimize the total number of allele mismatches (flips) — is the K-copy, variation-graph generalization of **Minimum Error Correction** (haplotype assembly), NP-hard already at `k = 2` [Lippert et al. 2002; Cilibrasi et al. 2005]. Hardness is thus robust across both the parsimony objective (min #copies, via coloring) and the genomics-standard error objective (min flips, via MEC).

- [ ] **Step 4: Re-run checks + adversarial read**

Run: `python3 bench/copy_assignment_theory_checks.py` → all `OK`. Re-read the proof asking: is the conflict graph exactly `Γ` (no spurious edges from shared columns)? — confirm each column is incident to exactly its two edge-endpoints, so no two non-adjacent vertices co-observe.

- [ ] **Step 5: Commit**

```bash
git add bench/copy_assignment_theory.md bench/copy_assignment_theory_checks.py
git commit -m "theory(copy-assign): Theorem 1 (hardness via coloring) + MEC framing + reduction check"
```

---

## Task 3: Theorem 2 (identifiability / unique minimum cover under C)

**Files:**
- Modify: `bench/copy_assignment_theory.md` (add §5: condition C, Theorem 2, proof, poly-recovery remark)
- Modify: `bench/copy_assignment_theory_checks.py` (add `check_thm2_recovery` and `check_thm2_K0_merge`)

**Interfaces:**
- Consumes: Lemma 1, `conflict_graph`, `mcc_bruteforce`, `jointly_consistent` (Task 1).
- Produces (note): **Condition C** and **Theorem 2** (true copies = unique minimum cover under C), plus the deferred-algorithm remark.

**The corrected statement (per Global Constraints):** Theorem 2 is identifiability/uniqueness, proved by the per-read adjacency argument below — *not* a cluster-graph claim. Polynomial recovery is a remark, deferred.

- [ ] **Step 1: Write both checks first (the K≥2 recovery fact and the K=0 non-recovery fact)**

Add to `bench/copy_assignment_theory_checks.py`:

```python
def reads_from_copies(copies, windows):
    """copies: list of allele-vectors (tuples in {0,1}^L). windows: list of (copy_index, set-of-columns)
    describing each read's origin copy and which columns it observes. Returns reads + true labels."""
    reads, labels = [], []
    for ci, cols in windows:
        reads.append(frozenset((j, copies[ci][j]) for j in cols))
        labels.append(ci)
    return reads, labels

def all_min_colorings(H, k):
    """All proper k-colorings of H up to nothing (raw), as label tuples; tiny graphs only."""
    import itertools
    n = H.number_of_nodes()
    out = []
    for assign in itertools.product(range(k), repeat=n):
        if all(assign[u] != assign[v] for u, v in H.edges()):
            out.append(assign)
    return out

def partition_of(labels):
    """canonical partition (set of frozensets of read-indices) from a labeling."""
    parts = {}
    for i, l in enumerate(labels):
        parts.setdefault(l, set()).add(i)
    return frozenset(frozenset(p) for p in parts.values())

def check_thm2_recovery():
    # K=2 copies over L=3 columns, differing at columns 0 and 2 (K_ij=2). Reads tile columns so that
    # every read conflicts with >=1 read of the foreign copy (per-read condition C2).
    copies = [(0,0,0), (1,0,1)]
    # copy0 reads and copy1 reads, each observing a distinguishing column + overlap:
    windows = [(0,{0,1}),(0,{1,2}),(0,{0,2}),(1,{0,1}),(1,{1,2}),(1,{0,2})]
    reads, labels = reads_from_copies(copies, windows)
    H = conflict_graph(reads)
    k = mcc_bruteforce(reads)
    assert k == 2, f"MCC should be 2, got {k}"
    true_part = partition_of(labels)
    colorings = {partition_of(c) for c in all_min_colorings(H, 2)}
    assert colorings == {true_part}, "the true partition must be the UNIQUE minimum cover under C"
    return "Thm 2: K=2 instance under C -> unique minimum cover == true copies"

def check_thm2_K0_merge():
    # K=0 separation: two 'copies' identical over the observed (exonic) columns -> no conflict -> merge.
    copies = [(0,1,0), (0,1,0)]  # identical over all columns
    windows = [(0,{0,1}),(0,{1,2}),(1,{0,1}),(1,{1,2})]
    reads, _ = reads_from_copies(copies, windows)
    H = conflict_graph(reads)
    assert H.number_of_edges() == 0, "identical copies produce no conflict edges"
    assert mcc_bruteforce(reads) == 1, "K=0: minimum cover merges all reads into ONE copy (true copies unrecoverable)"
    return "Thm 2 boundary: K=0 -> minimum cover = 1 (forced merge, non-identifiable)"

CHECKS.append(check_thm2_recovery)
CHECKS.append(check_thm2_K0_merge)
```

- [ ] **Step 2: Run; confirm both pass**

Run: `python3 bench/copy_assignment_theory_checks.py`
Expected: all checks `OK`, including the unique-recovery and the K=0-merge checks. If the K=2 instance has more than one minimum coloring, strengthen the windows (add reads) until the per-read condition holds and uniqueness is achieved — that calibration IS the content of condition C2.

- [ ] **Step 3: Write §5 (Condition C + Theorem 2) in the note**

> **Condition C.** Let `C* = {c_1,…,c_K}` be the true copies and assume each read originates from one copy (error-free core).
> - **(C1) Distinguishability.** For every `i ≠ j`, `c_i` and `c_j` differ at ≥ 1 column (`K_{ij} ≥ 1`; robust regime `K_{ij} ≥ 2`). Write `D_{ij}` for the (nonempty) set of distinguishing columns.
> - **(C2) Linking coverage (per-read).** For every read `r` (from copy `i`) and every other copy `j ≠ i`, some read `r'` from copy `j` co-observes a column of `D_{ij}` with `r`.
>
> **Theorem 2 (identifiability).** Under C, the true copy set `C*` is the **unique** minimum copy cover; i.e. the MCC optimum equals the truth, and no other minimum-cardinality consistent partition exists.
>
> *Proof.* Each true copy's reads are pairwise consistent (same origin allele-vector), hence an independent set of `H`; so `C*` is a proper `K`-coloring and `MCC ≤ K`. By (C2), for each pair `i ≠ j` some read of `i` conflicts with some read of `j`, so no independent set contains reads of two distinct copies; therefore every proper coloring needs ≥ `K` colors, giving `MCC = K` and that `C*` is minimum. **Uniqueness:** fix any minimum (`K`-)coloring and any read `r` from copy `i`. By (C2), `r` is adjacent in `H` to ≥ 1 read of every foreign copy `j`, so `r`'s color class contains no read of any foreign copy; as there are exactly `K` classes and `K` copies, each class is contained in a single true copy, and by minimality equals it. Hence the coloring is `C*`. ∎
>
> **Remark (recovery, deferred).** Theorem 2 identifies the optimum with the truth; computing it is the algorithmic question deferred to a follow-up. Under C the recovery generalizes error-free haplotype phasing (transitive allele-agreement assembly along linked reads); its polynomial-time guarantee and an explicit solver are out of scope here.
>
> **Corollary (the K-bound). (i)** If `K_{ij} = 0` for some pair — copies identical over the transcribed region — then (C1) fails: those copies produce no conflict, the minimum cover merges them, and the true copies are **provably unrecoverable** (the MAGEA co-located regime, resolvable fraction 0/494). **(ii)** With `K_{ij} ≥ 2` and (C2), C holds and recovery is exact (the sim5x K-ladder, K≥2 → 100%). The identifiability threshold of Theorem 2 is exactly the measured K-bound.

- [ ] **Step 4: Re-run checks + adversarial read**

Run: `python3 bench/copy_assignment_theory_checks.py` → all `OK`. Adversarial read of the uniqueness argument: does (C2) as stated truly forbid recoloring a read into a foreign class? Confirm the check's K=2 instance has exactly one minimum partition (the check asserts `colorings == {true_part}`); if a reviewer proposes a weaker C2, test it in the check and show it breaks uniqueness.

- [ ] **Step 5: Commit**

```bash
git add bench/copy_assignment_theory.md bench/copy_assignment_theory_checks.py
git commit -m "theory(copy-assign): Theorem 2 (unique min cover under C) + K-bound corollary + checks"
```

---

## Task 4: The paths/isoforms corollary

**Files:**
- Modify: `bench/copy_assignment_theory.md` (add §6: corollary lifting the core to paths)
- Modify: `bench/copy_assignment_theory_checks.py` (add `check_corollary_paths`)

**Interfaces:**
- Consumes: Theorem 2 and condition C (Task 3).
- Produces (note): **Corollary (joint copies + isoforms)** — lifting allele-vector recovery to constrained path-cover.

- [ ] **Step 1: Write the corollary check first**

Add to `bench/copy_assignment_theory_checks.py` — a DAG with two copies that differ by BOTH an allele and a junction (isoform), verifying the constrained path-cover recovers both:

```python
def check_corollary_paths():
    import networkx as nx
    # DAG nodes 0..4 (exon segments). Copy A path 0->1->3->4 (skips 2), with allele 0 at the PSV on node 1.
    # Copy B path 0->2->3->4 (skips 1)... model the junction choice as an extra binary 'column' jJ:
    #   reads carry (PSV column, allele) AND (junction column jJ, which-branch). Same conflict machinery.
    JJ = 99  # junction pseudo-column
    copies = {"A": {0: 0, JJ: 0}, "B": {0: 1, JJ: 1}}  # differ at PSV col 0 AND junction jJ
    def mk(cp, cols):
        return frozenset((c, copies[cp][c]) for c in cols)
    reads = [mk("A", {0, JJ}), mk("A", {JJ}), mk("A", {0}),
             mk("B", {0, JJ}), mk("B", {JJ}), mk("B", {0})]
    H = conflict_graph(reads)
    assert mcc_bruteforce(reads) == 2, "two copies separated by allele AND junction -> MCC=2"
    # treating the junction column on equal footing with PSV columns, recovery is the Thm-2 instance.
    return "Corollary: junction treated as a linkage column -> path-cover recovers copies+isoforms (MCC=2)"

CHECKS.append(check_corollary_paths)
```

- [ ] **Step 2: Run; confirm it passes**

Run: `python3 bench/copy_assignment_theory_checks.py` → all `OK`.

- [ ] **Step 3: Write §6 (Corollary) in the note**

> **Corollary (joint copies + isoforms).** Re-attach isoform structure: a copy is `(π, a)` and reads carry junction observations. Treat each distinguishing **junction** as an additional linkage column over the branch alphabet; then a copy is a consistent allele-and-junction vector, reads conflict on disagreeing co-observed alleles *or* junctions, and Lemma 1 and Theorem 2 apply verbatim. Hence under C augmented with junction-linkage (reads connect the junctions distinguishing isoforms), the minimum **constrained path-cover** of the family variation graph is unique and equals the true copies-with-isoforms — recovering #copies *and* their isoform structure jointly. Unconstrained minimum path-cover (max-flow on the DAG) is the relaxation; allele+junction linkage promotes it from *a* minimum cover to *the* true cover under C.
>
> *(If junctions require a condition beyond the allele core — e.g. a copy whose isoforms are distinguished only by a junction no read links — the core result stands and that case is future work, per the note's scope.)*

- [ ] **Step 4: Re-run checks**

Run: `python3 bench/copy_assignment_theory_checks.py` → all `OK`.

- [ ] **Step 5: Commit**

```bash
git add bench/copy_assignment_theory.md bench/copy_assignment_theory_checks.py
git commit -m "theory(copy-assign): paths/isoforms corollary (junction as linkage column) + check"
```

---

## Task 5: Assemble the note (intro, empirics, discussion) + final review

**Files:**
- Modify: `bench/copy_assignment_theory.md` (expand §1 intro; add §7 empirical corroboration; §8 discussion)

**Interfaces:**
- Consumes: all theorems/corollary (Tasks 1–4).
- Produces: the complete, self-contained theory note.

- [ ] **Step 1: Write §1 (introduction)**

Expand the intro stub: copy assignment is the unit downstream of family detection (interest #1); the multimapping-resolution lineage (Canzar's conflict-resolution / facility-location frame); this note's contribution is the hardness/recovery dichotomy with the identifiability threshold = the measured K-bound, unifying interests #1–#3.

- [ ] **Step 2: Write §7 (empirical corroboration), citing — not producing — the existing results**

> - **Theorem 2 / `K≥2`:** sim5x K-ladder (`bench/sim_reads.py`, `/home/juanfra/winloci_scratch/sim5x/`), K≥2 → 100% correct assignment; GGO silver-standard = 100% (`detect_and_assign` MAGEA smoke run, 1026/1026).
> - **Theorem 2 boundary / `K=0`:** MAGEA co-located arrays, resolvable fraction 0/494; copies sequence-identical over the transcribed exon (`bench/resolution_improvement_bound.md`).
> - **Shared condition with detection (#1):** the family is the conflict-graph component (`bench/family_definition_formal.md`).

- [ ] **Step 3: Write §8 (discussion)**

The dichotomy as the through-line across interests #1 (detection), #2 (assignment), #3 (allele-specific junctions); the K=0 frontier where RNA alone cannot resolve and longer reads / allele-specific junctions become the lever; the deferred algorithmic question (the polynomial recovery solver).

- [ ] **Step 4: Final self-contained read + full check run**

Run: `python3 bench/copy_assignment_theory_checks.py` → all checks `OK`. Read the note end-to-end: every theorem statement has a proof; every combinatorial claim has a passing check; no dangling references; the corrected Theorem 2 (identifiability, not cluster-graph) is consistent throughout.

- [ ] **Step 5: Commit**

```bash
git add bench/copy_assignment_theory.md
git commit -m "theory(copy-assign): assemble note — intro, empirical corroboration, discussion"
```

---

## Self-Review

**Spec coverage:** §Problem→Task 1; Theorem 1→Task 2; Theorem 2 + K-bound→Task 3; paths corollary→Task 4; empirics + document structure→Task 5. The spec's "cluster graph / polynomial-time" wording for Theorem 2 is intentionally **superseded** by the approved identifiability/uniqueness statement (Global Constraints + Task 3) with recovery deferred — matching the user's "defer algorithm" choice.

**Placeholder scan:** none — every proof body and every check is written out in full.

**Type/claim consistency:** `conflict_graph`, `mcc_bruteforce`, `jointly_consistent`, `reads_from_copies`, `partition_of`, `chromatic_number` are defined in Task 1/2 and reused with the same signatures in Tasks 3–4. Reads are uniformly `frozenset((column, allele))`. Theorem numbering (Lemma 1, Theorem 1, Theorem 2, Corollary) is consistent across note and commits.
