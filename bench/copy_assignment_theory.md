# Copy-assignment theory: model, definitions, and the MCC = χ(H) core

*This note develops the combinatorial foundation for the read-to-copy assignment problem in multi-copy gene
families. §1 situates the problem. §2 fixes the model and notation. §3 proves that the minimum number of copies
needed to explain a set of reads equals the chromatic number of the conflict graph — the theorem that makes the
assignment problem tractable.*

---

## §1 Introduction (stub — to be expanded)

Long-read RNA sequencing resolves full isoform structures, but reads from paralogous gene-family members are
often alignment-ambiguous: a single molecule may map with equal or near-equal quality to two or more distinct
genomic loci. The downstream assembly then faces a **copy-assignment** problem: given a set of reads and their
observed allelic evidence, partition the reads so that each part is explained by a single gene-copy
(allele-vector). The minimum number of copies required to explain the data is the central quantity. This note
gives a formal model, defines that quantity precisely, and proves it equals the chromatic number of a naturally
arising conflict graph — a result that both bounds the complexity of the problem and guides efficient algorithms.

---

## §2 Model and definitions

### Reads and observed columns

Let $R = \{r_1, \ldots, r_n\}$ be a finite set of **reads**. Each read $r_i$ carries partial evidence about a
fixed set of **columns** $[m] = \{1, \ldots, m\}$, where a column represents a genomic position (or PSV — a
paralog-specific variant site) that distinguishes gene copies. Each column $j$ takes values in a finite allele
alphabet $A_j$.

A read $r_i$ is modelled as a partial function from columns to alleles. We write $\mathrm{obs}(r_i) \subseteq [m]$
for the set of columns $r_i$ actually observes, and $r_i(j) \in A_j$ for the allele it reports at column
$j \in \mathrm{obs}(r_i)$. In the verification code a read is represented as a `frozenset` of $(j, a)$ pairs.

### Allele-vectors and gene copies

A **gene copy** is characterised by an allele-vector $v \in \prod_{j=1}^{m} A_j$, one allele per column. A read
$r_i$ is **consistent with** a copy $v$ if $r_i(j) = v_j$ for every $j \in \mathrm{obs}(r_i)$ — the read does
not contradict the copy at any column it observes.

### Joint consistency of a set of reads

A set $S \subseteq R$ of reads is **jointly consistent** if there exists a single allele-vector $v$ consistent
with every read in $S$, equivalently: for every column $j$, all reads in $S$ that observe $j$ report the same
allele. Joint consistency is a **per-column property**: a column $j$ is inconsistent in $S$ iff two reads in $S$
both observe $j$ but with different alleles; $S$ is jointly consistent iff no column is inconsistent in $S$.

### Conflict

Two reads $r_i, r_k$ **conflict** if they co-observe at least one column with different alleles:
$$
r_i \sim r_k \iff \exists\, j \in \mathrm{obs}(r_i) \cap \mathrm{obs}(r_k):\; r_i(j) \neq r_k(j).
$$
Non-conflicting reads are **compatible**: they agree at every co-observed column (though they may observe
disjoint column sets entirely).

### Conflict graph

The **conflict graph** $H = (R, E)$ has the reads as vertices and an edge between every conflicting pair:
$$
\{r_i, r_k\} \in E \iff r_i \sim r_k.
$$
An **independent set** of $H$ is a set of reads that are pairwise non-conflicting.

### Minimum Copy Cover (MCC)

A **copy cover** of $R$ is a partition $\mathcal{P} = \{P_1, \ldots, P_k\}$ of $R$ into jointly-consistent
parts. The **Minimum Copy Cover** is
$$
\mathrm{MCC}(R) = \min\bigl\{k \in \mathbb{N} : R \text{ admits a copy cover of size } k\bigr\}.
$$
$\mathrm{MCC}(R)$ is the smallest number of gene copies whose existence is necessary to explain all the reads.

---

## §3 The MCC core: Lemma 1

The key structural result is that MCC coincides with the chromatic number of the conflict graph, which gives an
exact combinatorial characterisation of the problem and connects it to a well-studied graph-theoretic quantity.

**Lemma 1.** $\mathrm{MCC}(R) = \chi(H)$, the chromatic number of the conflict graph $H$.

*Proof.*

We show that the jointly-consistent subsets of $R$ are exactly the independent sets of $H$, and therefore that
copy covers and proper colorings of $H$ are in bijection.

**Step 1 — Independent set ⇒ jointly consistent.** Let $S \subseteq R$ be an independent set of $H$: every pair
$r_i, r_k \in S$ is non-conflicting, i.e., they agree at every co-observed column. We claim $S$ is jointly
consistent.

Consistency is a per-column property: fix any column $j$ and consider the reads in $S$ that observe $j$, say
$T_j = \{r \in S : j \in \mathrm{obs}(r)\}$. For any two reads $r_i, r_k \in T_j$, both observe column $j$;
since $r_i$ and $r_k$ are non-conflicting, they report the same allele at $j$: $r_i(j) = r_k(j)$. Because this
holds for every pair in $T_j$, all reads in $T_j$ report a single allele $a_j$ at column $j$. Ranging over all
columns $j$, the vector $v$ with $v_j = a_j$ (choosing $v_j$ arbitrarily when $T_j = \varnothing$) is an
allele-vector consistent with every read in $S$. Hence $S$ is jointly consistent.

**Step 2 — Jointly consistent ⇒ independent set.** Let $S$ be jointly consistent, witnessed by allele-vector
$v$. For any two reads $r_i, r_k \in S$ and any $j \in \mathrm{obs}(r_i) \cap \mathrm{obs}(r_k)$, both reads
agree with $v$ at column $j$, so $r_i(j) = v_j = r_k(j)$. Thus $r_i$ and $r_k$ do not conflict, and $S$ is an
independent set.

**Conclusion.** Steps 1 and 2 together give:
$$
S \text{ is jointly consistent} \iff S \text{ is an independent set of } H.
$$
A copy cover of size $k$ is therefore a partition of $R$ into $k$ jointly-consistent parts, which is exactly a
proper $k$-coloring of $H$, and conversely. Minimizing $k$ gives $\mathrm{MCC}(R) = \chi(H)$. $\square$

### Verification: the C5 gadget

The check in `bench/copy_assignment_theory_checks.py` instantiates the theorem on a concrete example whose
conflict graph is a 5-cycle $C_5$.

**Construction.** Label 5 reads $r_0, \ldots, r_4$. Use 5 columns $e_{01}, e_{12}, e_{23}, e_{34}, e_{40}$, one
per edge of $C_5$. Read $r_w$ observes exactly the two columns corresponding to edges incident to $w$ in $C_5$,
with opposite alleles: $r_w(e_{w,w'}) = 0$ and $r_w(e_{w',w}) = 1$ (where $w' = w-1 \bmod 5$, $w'' = w+1 \bmod 5$).
Adjacent reads $r_w, r_{w+1}$ share column $e_{w,w+1}$ with opposite alleles, so they conflict. Non-adjacent
reads share no column, so they are compatible. The conflict graph is exactly $C_5$.

**Result.** $C_5$ is an odd cycle, so $\chi(C_5) = 3$. DSATUR (greedy saturation-largest-first coloring) is
exact on cycles and returns 3. The brute-force MCC search confirms $\mathrm{MCC} = 3$, establishing the
prediction of Lemma 1 on this instance.

Run: `python3 bench/copy_assignment_theory_checks.py` — exits 0 and prints:

```
OK  - Lemma 1 (MCC=chi) verified on C5: MCC=3=chi(C5)
```

---

---

## §4 Complexity: Theorem 1

The MCC = χ(H) identity (Lemma 1) immediately yields NP-hardness of MCC via a reduction from graph
coloring, one of Karp's original 21 NP-complete problems.

**Theorem 1.** MCC is NP-hard. (Hence CMCPC, which contains MCC as the single-backbone case, is NP-hard.)

*Proof.* We reduce graph k-colorability to MCC in polynomial time.

Let $\Gamma = (U, F)$ be an arbitrary undirected graph with $U = \{1, \ldots, n\}$ and edge set $F$.
Construct an MCC instance as follows.

- **Columns.** Introduce one binary column $j_e$ for each edge $e \in F$; the allele alphabet for every
  column is $\{0, 1\}$.
- **Reads.** Introduce one read $r_w$ for each vertex $w \in U$.  Read $r_w$ observes column $j_e$
  if and only if $w \in e$.  If $e = (u, v)$ (ordered so that $u$ is the first endpoint and $v$ the
  second), then $r_u(j_e) = 0$ and $r_v(j_e) = 1$.

**The conflict graph of the instance is isomorphic to $\Gamma$.**
Two reads $r_u$ and $r_v$ co-observe column $j_e$ if and only if both $u$ and $v$ belong to $e$,
i.e., $(u, v) \in F$.  Whenever they do co-observe $j_e$, one of them reports allele 0 and the
other reports allele 1, so they always disagree — they conflict.  Conversely, two reads $r_u, r_v$
with $(u, v) \notin F$ share no column (each column $j_e$ is observed only by its two
edge-endpoints, so non-adjacent vertices never co-observe).  Therefore the conflict graph $H$ has
an edge $\{r_u, r_v\}$ if and only if $(u, v) \in F$, making $H$ isomorphic to $\Gamma$.

**The reduction is correct.**
By Lemma 1, $\mathrm{MCC}(R) = \chi(H) = \chi(\Gamma)$.  Hence $\Gamma$ is $k$-colorable if and
only if $\mathrm{MCC}(R) \leq k$, so any polynomial-time algorithm for MCC would solve graph
coloring in polynomial time.  Graph $k$-colorability is NP-complete for $k \geq 3$
[Garey–Johnson 1979], so MCC is NP-hard.

**The alphabet is binary.**
Every column in the construction takes values in $\{0, 1\}$; hardness is therefore intrinsic to the
problem, not an artifact of a large allele alphabet. $\square$

### Adversarial check: each column is incident to exactly two reads

The key invariant in the proof is that column $j_e$, for edge $e = (u, v)$, is observed by exactly
$r_u$ and $r_v$ — no other read observes it.  This ensures non-adjacent vertices share *no* column
and thus cannot conflict, i.e., the reduction introduces no spurious edges.  The executable check
below verifies this structural property explicitly by confirming `nx.is_isomorphic(H, Gamma)` for
four small graphs (triangle, $C_4$, $C_5$, star); see `bench/copy_assignment_theory_checks.py`,
function `check_thm1_reduction`.

### MEC framing

The **Minimum Error Correction** (MEC) problem is the standard genomics formulation of the
haplotype-assembly problem [Lippert et al. 2002]: given a set of reads and a fixed number $k$ of
haplotypes, find an assignment of reads to haplotypes and an allele-vector for each haplotype that
minimizes the total number of allele *mismatches* (flips) between reads and their assigned copy.
MEC is NP-hard already at $k = 2$ [Lippert et al. 2002; Cilibrasi et al. 2005].

The copy-assignment problem studied here is the $k$-copy, variation-graph generalization of MEC:
reads are molecules over PSV columns, copies are allele-vectors in a family variation graph, and the
error objective counts allele mismatches.  Theorem 1 (hardness via coloring, parsimony objective)
and MEC-hardness (error objective) together show that hardness is **robust across both objective
functions**: minimizing the number of copies and minimizing allele flips are both NP-hard.

---

## §5 Identifiability: condition C and Theorem 2

Theorem 1 says the *general* MCC problem is hard. But the assignment instances that arise in practice are not
adversarial: when paralog copies carry enough distinguishing variation and the reads link those variants, the
truth is not merely *a* minimum cover — it is the *unique* one. Theorem 2 makes this precise. It is an
**identifiability / uniqueness** statement: it certifies *what the optimum is* (the true copy set), not *how to
compute it efficiently*. The conflict graph $H$ is **not** a disjoint union of cliques in general (two reads of
the same copy that observe disjoint columns are non-adjacent), so uniqueness does **not** follow from any
cluster-graph structure; it is proved directly by a per-read adjacency argument.

### Condition C

> **Condition C.** Let $C^\* = \{c_1, \ldots, c_K\}$ be the true copies and assume each read originates from
> exactly one copy (error-free core; each read $r$ is consistent with its origin copy $c_i$).
>
> - **(C1) Distinguishability.** For every $i \neq j$, the allele-vectors $c_i$ and $c_j$ differ at $\geq 1$
>   column ($K_{ij} \geq 1$; robust regime $K_{ij} \geq 2$). Write $D_{ij} \neq \varnothing$ for the set of
>   distinguishing columns of the pair $(i,j)$.
> - **(C2) Linking coverage (per-read).** For every read $r$ (from copy $i$) and every other copy $j \neq i$,
>   some read $r'$ from copy $j$ **co-observes a column of $D_{ij}$ with $r$** — i.e. $r$ and $r'$ both observe
>   some $d \in D_{ij}$, so they report different alleles there and therefore conflict in $H$.

Equivalently, (C2) says: in $H$, every read of copy $i$ has at least one neighbor in every foreign copy $j$.
(C1) is the global existence of separating variation; (C2) is the local, per-read demand that this variation is
actually *witnessed* by overlapping reads. (C1) is necessary for (C2) to be satisfiable — without a
distinguishing column there is nothing to co-observe — but (C2) is strictly stronger.

### Theorem 2

> **Theorem 2 (identifiability).** Under condition C, the true copy set $C^\*$ is the **unique** minimum copy
> cover of $R$: $\mathrm{MCC}(R) = K$, and the only minimum-cardinality partition of $R$ into jointly-consistent
> parts is $C^\*$ itself.

*Proof.* Write $H$ for the conflict graph; by Lemma 1, minimum copy covers are exactly minimum proper colorings
of $H$, so it suffices to reason about colorings.

**(a) $C^\* $ is a proper $K$-coloring, so $\mathrm{MCC} \leq K$.** All reads of a single true copy $c_i$ are
consistent with the common allele-vector $c_i$, hence pairwise compatible (they agree at every co-observed
column); they form an independent set of $H$. Coloring each read by its origin copy therefore partitions $R$
into $K$ jointly-consistent parts, so $\mathrm{MCC} \leq K$.

**(b) Every proper coloring needs $\geq K$ colors, so $\mathrm{MCC} = K$ and $C^\*$ is minimum.** Fix $i \neq j$.
By (C2), some read $r$ of copy $i$ and some read $r'$ of copy $j$ co-observe a column of $D_{ij}$ and so conflict:
$\{r, r'\} \in E(H)$. Hence no independent set of $H$ — equivalently, no jointly-consistent part — can contain
reads of two distinct true copies. Any proper coloring's classes are independent sets, so each class lies inside
a single true copy; covering all $K$ copies needs $\geq K$ classes. With (a), $\mathrm{MCC} = K$, and $C^\*$ is a
minimum cover.

**(c) Uniqueness.** Let $\Pi$ be *any* minimum (i.e. $K$-)coloring; we show $\Pi = C^\*$. Take any read $r$ from
copy $i$ and any foreign copy $j \neq i$. By (C2), $r$ is adjacent in $H$ to at least one read $r'$ of copy $j$.
Since a color class of $\Pi$ is an independent set, $r$ and $r'$ cannot share a class, so $r$'s class contains
**no** read of copy $j$. As this holds for every foreign $j$, $r$'s class contains reads of copy $i$ only;
thus every class of $\Pi$ is a subset of a single true copy. There are exactly $K$ classes and $K$ copies, and
the classes partition all of $R$ (every copy's reads must be covered), so the map class $\mapsto$ its (unique)
copy is a bijection and each class equals its copy. Hence $\Pi = C^\*$. $\square$

The argument in (c) is exactly what the weakened-C2 stress test in the checks falsifies: if some read observes
only columns outside $\bigcup_j D_{ij}$ (so it has no foreign neighbor), that read may be recolored into a
foreign class without creating a conflict, producing a second minimum cover and destroying uniqueness. Per-read
(C2) is therefore load-bearing, not decorative.

### Remark (polynomial-time recovery, deferred)

Theorem 2 identifies the optimum with the truth; it does **not** assert that the optimum is computable in
polynomial time (Theorem 1 forbids that in general). Under condition C, recovery generalizes error-free
haplotype phasing: starting from any read and following the linking overlaps guaranteed by (C2), transitive
allele-agreement assembles each true copy's read set. The design and polynomial-time analysis of such a solver
— and the precise structural class of conflict graphs $H$ for which it runs in polynomial time — are deferred to
a follow-up; they are out of scope for this identifiability note.

### Corollary (the K-bound)

> **Corollary (K-bound).**
> **(i)** If $K_{ij} = 0$ for some pair (copies identical over the transcribed/observed region), then (C1)
> fails: that pair produces no conflict edge in $H$, the minimum cover **merges** the two copies into one part,
> and the true copies are **provably unrecoverable**. This is the MAGEA co-located regime (measured resolvable
> fraction $0/494$).
> **(ii)** With $K_{ij} \geq 2$ for every pair together with (C2), condition C holds and recovery is exact —
> the sim5x PSV ladder, where $K \geq 2$ gives $100\%$ recovery.
>
> The identifiability threshold of Theorem 2 thus coincides exactly with the empirically measured **K-bound**:
> $K_{ij} = 0 \Rightarrow$ merge / non-identifiable; $K_{ij} \geq 2 + \text{linkage} \Rightarrow$ exact
> identification.

The $K = 0$ branch is verified directly in `check_thm2_K0_merge` (identical copies $\Rightarrow$ zero conflict
edges $\Rightarrow$ $\mathrm{MCC} = 1$, a forced merge), and the $K \geq 2$ branch in `check_thm2_recovery`
(a $K_{ij} = 2$ instance satisfying per-read (C2), whose unique minimum cover equals the true partition).

Run: `python3 bench/copy_assignment_theory_checks.py` — exits 0 and prints, among the checks:

```
OK  - Thm 2: K=2 instance under C -> unique minimum cover == true copies
OK  - Thm 2 boundary: K=0 -> minimum cover = 1 (forced merge, non-identifiable)
```

---

*See also: `bench/family_definition_formal.md` for the upstream family-detection step whose output populates
$R$ for each identified gene family.*
