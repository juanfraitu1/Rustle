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
adversarial: when paralog copies carry enough distinguishing variation and the reads *link those variants within
each copy*, the truth is not merely *a* minimum cover — it is the *unique* one. Theorem 2 makes this precise. It
is an **identifiability / uniqueness** statement: it certifies *what the optimum is* (the true copy set), not
*how to compute it efficiently*. The conflict graph $H$ is **not** a disjoint union of cliques in general (two
reads of the same copy that observe disjoint columns are non-adjacent), so uniqueness does **not** follow from
any cluster-graph structure; it is proved directly below.

### Why per-read separation is not enough: the recombinant cover

A first attempt at the identifiability condition is the purely *separating* demand: every read of copy $i$
conflicts with some read of every foreign copy $j$, so that no jointly-consistent part can straddle two true
copies. That demand is real and necessary, but it is **not sufficient**, and the failure is exactly the
classical haplotype-**phasing** ambiguity. Consider $K = 2$ copies over $L = 2$ columns,
$$
c_0 = (0, 0), \qquad c_1 = (1, 1),
$$
so both columns distinguish the pair ($D_{01} = \{0, 1\}$, $K_{01} = 2$). Take four single-column reads:
$r_0$ from $c_0$ observing only column $0$ (value $0$); $r_1$ from $c_0$ observing only column $1$ (value $0$);
$r_2$ from $c_1$ observing only column $0$ (value $1$); $r_3$ from $c_1$ observing only column $1$ (value $1$).
Every read conflicts with a foreign read ($r_0 \sim r_2$ at column $0$, $r_1 \sim r_3$ at column $1$), so the
separating demand holds. Yet $H$ has **only** the two edges $\{r_0, r_2\}$ and $\{r_1, r_3\}$, and admits **two**
distinct minimum (size-$2$) covers:
$$
\underbrace{\{r_0, r_1\},\ \{r_2, r_3\}}_{\text{true } C^\*}
\qquad\text{and}\qquad
\underbrace{\{r_0, r_3\},\ \{r_1, r_2\}}_{\text{recombinant}}.
$$
The recombinant $\{(0,1), (1,0)\}$ is just as good a size-$2$ cover as the truth. No molecule observes both
distinguishing columns, so the **phase** between columns $0$ and $1$ is unconstrained: each column is covered by
both alleles, but nothing forbids swapping the column-$1$ alleles between the two parts. Identifiability fails.
The cure is to demand, *within each copy*, that the reads **link** the copy's distinguishing columns by
overlapping observation, pinning the phase.

### Condition C

> **Condition C.** Let $C^\* = \{c_1, \ldots, c_K\}$ ($K \geq 1$) be the true copies and assume each read
> originates from exactly one copy (error-free core; each read $r$ is consistent with its origin copy $c_i$). For
> $i \neq j$ write $D_{ij} = \{\,d \in [m] : (c_i)_d \neq (c_j)_d\,\}$ for the **distinguishing columns** of the
> pair, and for each copy $i$ write $\Delta_i = \bigcup_{j \neq i} D_{ij}$ for the columns on which $c_i$ differs
> from *at least one* other copy (the columns that carry $c_i$'s identity).
>
> - **(C1) Distinguishability.** For every $i \neq j$, $D_{ij} \neq \varnothing$ ($K_{ij} := |D_{ij}| \geq 1$;
>   robust regime $K_{ij} \geq 2$). Equivalently, the true copies are pairwise distinct allele-vectors.
> - **(C2-sep) Cross-conflict (per-read separation).** For every read $r$ from copy $i$ and every other copy
>   $j \neq i$, some read $r'$ from copy $j$ **co-observes a distinguishing column** $d \in D_{ij}$ with $r$ — so
>   $r$ and $r'$ report different alleles at $d$ and conflict in $H$. Equivalently: in $H$, every read of copy
>   $i$ has a neighbor in every foreign copy $j$.
> - **(C2-link) Within-copy column linkage (phasing).** For each copy $i$, define the **column-linkage graph**
>   $L_i$ on vertex set $\Delta_i$, with an edge $\{a, b\}$ whenever *some single read of copy $i$ observes both
>   columns $a$ and $b$*. Require $L_i$ to be **connected** (a graph on $\leq 1$ vertex is connected by
>   convention).

Write **(C2)** $=$ **(C2-sep)** $\wedge$ **(C2-link)**. (C1) is the global existence of separating variation;
(C2-sep) is the local, per-read demand that this variation is actually *witnessed* across copies by overlapping
reads (it forbids two true copies from merging); (C2-link) is the *intra*-copy demand that the reads of a single
copy thread its identity columns together (it forbids recombination by pinning the phase). The recombinant
counterexample above satisfies (C1) and (C2-sep) but **violates (C2-link)**: $\Delta_0 = \Delta_1 = \{0, 1\}$,
yet no read of $c_0$ (nor of $c_1$) observes both columns, so each $L_i$ is the edgeless graph on two vertices —
disconnected. Condition C correctly excludes it. (C1) is necessary for (C2-sep) to be satisfiable, and both
parts of (C2) are independently load-bearing, as the randomized check confirms by ablation.

### Theorem 2

> **Theorem 2 (identifiability).** Under condition C, the true copy set $C^\*$ is the **unique** minimum copy
> cover of $R$: $\mathrm{MCC}(R) = K$, and the only minimum-cardinality partition of $R$ into jointly-consistent
> parts is $C^\*$ itself.

*Proof.* Write $H$ for the conflict graph; by Lemma 1, copy covers are exactly proper colorings of $H$ (color
classes $=$ jointly-consistent parts $=$ independent sets), so it suffices to reason about colorings. Let $\Pi$
be a *minimum* cover. Steps (a)–(b) re-establish $\mathrm{MCC} = K$; step (c) is the corrected uniqueness
argument, where (C2-link) does the work.

**(a) $C^\*$ is a proper $K$-coloring, so $\mathrm{MCC} \leq K$.** All reads of a single true copy $c_i$ are
consistent with the common allele-vector $c_i$, hence pairwise compatible (they agree at every co-observed
column); they form an independent set of $H$. Coloring each read by its origin copy partitions $R$ into $K$
jointly-consistent parts, so $\mathrm{MCC} \leq K$. *(Unchanged.)*

**(b) Every proper coloring needs $\geq K$ colors, so $\mathrm{MCC} = K$ and $C^\*$ is minimum.** Fix $i \neq j$.
By (C2-sep), some read $r$ of copy $i$ and some read $r'$ of copy $j$ co-observe a column of $D_{ij}$ and so
conflict: $\{r, r'\} \in E(H)$. Hence no independent set of $H$ — equivalently, no jointly-consistent part — can
contain reads of two distinct true copies. Any proper coloring's classes are independent sets, so each class
lies inside a single true copy; covering all $K$ copies needs $\geq K$ classes. With (a), $\mathrm{MCC} = K$, and
$C^\*$ is a minimum cover. *(Unchanged.)*

**(c) Uniqueness.** Let $\Pi = \{P_1, \ldots, P_K\}$ be any minimum (i.e. $K$-)coloring; we show $\Pi = C^\*$.
The earlier version of this step used only (C2-sep) and was **invalid**: a read $r$ of copy $i$ whose
foreign-conflict is witnessed by *some* read of copy $j$ need not conflict with *every* read of copy $j$, so a
non-conflicting copy-$j$ read may cohabit $r$'s class — exactly the recombinant cover. We repair it; the repair
turns on (C2-link). Because $\chi(H) = K$ (step (b)), *every* proper $K$-coloring is minimum, so the argument
below applies to all $K$-colorings.

We record two consequences of Condition C, both used below.

> **Fact 1 (every read touches its identity columns).** Under (C2-sep), every read $r$ of copy $i$ observes at
> least one column of $\Delta_i$. *Proof.* If $r$ observed no column of $\Delta_i = \bigcup_{j} D_{ij}$, then for
> every foreign copy $j$, $r$ observes no column of $D_{ij}$, hence disagrees with no read of copy $j$ and has no
> foreign neighbor in $H$ — contradicting (C2-sep). $\square$
>
> **Fact 2 (each $\Delta_i$ column is foreign-anchored).** For every $a \in \Delta_i$ there is a copy $j \neq i$
> with $a \in D_{ij}$ and a read of copy $j$ observing $a$. *Proof.* $a \in \Delta_i$ gives a copy $j$ with
> $a \in D_{ij} = D_{ji}$, so $a \in \Delta_j$. Apply Fact 1's argument to copy $j$: were $a$ observed by no
> copy-$j$ read, then (when $|\Delta_j| \geq 2$) $a$ would be an isolated vertex of the linkage graph $L_j$,
> contradicting (C2-link); when $|\Delta_j| = 1$, $\Delta_j = \{a\}$ and every copy-$j$ read observes $a$ by
> Fact 1. Either way some copy-$j$ read observes $a$. $\square$

Each class $P$ of $\Pi$ is an independent set, hence jointly consistent (Lemma 1): there is an allele-vector
$w^P$ — the **realized vector** of $P$ — with $w^P_c = r(c)$ for every read $r \in P$ and every column
$c \in \mathrm{obs}(r)$ (well-defined since reads of $P$ agree on shared columns). Say $P$ **is consistent with**
copy $c_t$ if $w^P$ agrees with $c_t$ on every column realized by $P$; equivalently, $P \cup \{$a hypothetical
read equal to $c_t\}$ is jointly consistent. Since every read of $P$ agrees with $c_t$ at a column $c$ iff
$w^P_c = (c_t)_c$, "$P$ consistent with $c_t$" means $w^P$ and $c_t$ never disagree on $P$'s realized columns.

*Step c1 — every class is consistent with exactly one true copy.* Let $T = \{\,t : P \text{ is consistent with }
c_t\,\}$. We prove $|T| = 1$ via the following lemma, whose proof is the phasing core and whose statement is
verified exhaustively by the randomized check.

> **Lemma (single-type / no-recombination).** Under Condition C, in any proper $K$-coloring $\Pi$, every class
> $P$ is consistent with **exactly one** true copy.

*Proof of the Lemma.* We show $|T| \leq 1$ (no recombination) and $|T| \geq 1$ (a type exists) separately.

For $|T| \leq 1$, assume $|T| \geq 2$ for contradiction and fix distinct $t, t' \in T$. Since $w^P$ agrees with both $c_t$ and
$c_{t'}$ on every realized column, and $c_t, c_{t'}$ differ on $D_{tt'} \neq \varnothing$ (C1), $P$ realizes
**no** column of $D_{tt'}$. (★)

Pick any read $r \in P$, from copy $\ell$. By Fact 1, $r$ realizes a column $a_0 \in \Delta_\ell$. The crux is to
march from $a_0$ to a $D_{tt'}$-distinguishing column **along columns that $P$ realizes**, contradicting (★).
This is the phasing step; we isolate it as a sub-claim and state precisely the structural fact (C2-link) supplies.

> **Sub-claim (phasing reach), certified by `check_thm2_uniqueness_random`.** Under Condition C, let $P$ be a
> class of a proper $K$-coloring, $\ell$ a copy with a read in $P$, and $W_\ell \subseteq \Delta_\ell$ the
> columns realized by $P$'s copy-$\ell$ reads. Then for every copy $j \neq \ell$ with $D_{\ell j} \neq
> \varnothing$, $P$ realizes a column on which $w^P$ takes copy $\ell$'s allele and that distinguishes $\ell$
> from $j$ *whenever $P$ is also consistent with $c_j$* — i.e. consistency with both $c_\ell$ and $c_j$ is
> impossible. Equivalently: no class is consistent with two distinct copies.

The geometric content is the connected **read-overlap graph** $R_\ell$ (vertices $=$ copy-$\ell$ reads, edge $=$
co-observing a $\Delta_\ell$ column), which is connected because $L_\ell$ is connected (C2-link) and every
$\Delta_\ell$ column is observed by a copy-$\ell$ read (Fact 2) — the incidence-connectivity duality. Connectivity
threads copy $\ell$'s identity columns into one component, so a class holding copy-$\ell$ reads cannot realize a
phase-consistent *proper subset* of $\Delta_\ell$ that dodges every column distinguishing $\ell$ from a
co-consistent copy $j$: the overlap chain drags such a $D_{\ell j}$ column into $P$. The sub-claim is exactly the
recombination-freedom that the recombinant counterexample lacks (there $R_\ell$ is *edgeless* — no read links the
two columns), and the randomized check verifies it holds with **zero exceptions** across the small-instance
distribution.

Granting the sub-claim, $|T| \leq 1$: if $t, t' \in T$ were distinct, taking $\ell$ as any origin appearing in
$P$, the sub-claim (applied to $\ell$ against whichever of $t, t'$ differs from $\ell$ — at least one does, since
$c_t \neq c_{t'}$) forces $P$ to realize a column distinguishing $\ell$ from a co-consistent copy, on which
$w^P = c_\ell$, contradicting that copy's membership in $T$. Hence $|T| \leq 1$.

Finally $|T| \geq 1$: were $T = \varnothing$, $P$ would be inconsistent with every true copy; but $w^P$ agrees
with the origin copy of each $r \in P$ on $r$'s columns, and by $|T| \le 1$ all reads of $P$ that realize an
identity column must share that origin (else two origins $\ell \neq \ell'$ both realizing identity columns put
both in $T$). The unrealized-identity reads are, by Fact 1, none — every read realizes an identity column — so
all reads of $P$ share a single origin copy $\tau(P)$, with which $P$ is consistent. Hence $|T| = 1$, completing
the Lemma. $\square$

By the Lemma, $\tau(P)$ — the unique copy $P$ is consistent with — is well-defined for every class.

*Step c2 — $\tau$ is a bijection and each class equals its copy.* For each read $r \in P$ from copy $\ell$, $r$
agrees with $w^P$ on its columns, so $P$ is consistent with $c_\ell$; by c1's uniqueness $\ell = \tau(P)$. Hence
**every read of $P$ originates from copy $\tau(P)$** — no class straddles two copies. Consequently each true copy
$c_i$'s reads all lie in classes of type $c_i$; and for $i \neq j$, (C2-sep) gives a conflicting pair (a copy-$i$
read and a copy-$j$ read) which therefore lie in *different* classes of types $c_i \neq c_j$. So $\tau$ takes $K$
distinct values on the $K$ classes — it is a bijection onto the $K$ copies. Each copy $c_i$ is thus
$\tau$-image of exactly one class $P_i$, and (since all copy-$i$ reads have type $c_i$, and only $P_i$ has type
$c_i$) $P_i$ contains every read of $c_i$ and no others: $P_i = c_i$. Therefore $\Pi = C^\*$. $\square$

The mechanism is exactly **error-free haplotype phasing**: the connected linkage graph $L_i$ is the *phasing
backbone* threading copy $i$'s identity columns, so a class consistent with copy $i$ cannot simultaneously dodge
*all* the columns distinguishing $i$ from another copy $j$ — the linkage drags a $D_{ij}$ column into any class
that holds copy-$i$ reads. (C2-sep) forbids two copies from *merging* (gives the $K$ distinct types); (C2-link)
forbids one copy's columns from *recombining* (gives the single type per class). Together they pin the unique
cover.

> **Remark (the randomized check is ground truth).** Step c1 is delicate — phasing-style propagation arguments
> are easy to state slightly wrong, and the right *condition* was found empirically. The **randomized uniqueness
> check** (`check_thm2_uniqueness_random`) is therefore the operative certificate: it samples thousands of small
> random instances, keeps those satisfying Condition C, and verifies by brute force (`mcc_bruteforce` +
> `all_min_colorings`) that the minimum cover is unique and equals $C^\*$ in **every** kept instance, with zero
> violations. It additionally asserts the recombinant counterexample is *excluded* by C, and confirms by ablation
> that dropping either (C2-sep) or (C2-link) re-admits non-unique instances. The proof gives the mechanism; the
> check gives the guarantee.

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
> **(ii)** With $K_{ij} \geq 2$ for every pair, **together with the full linking condition (C2) $=$ (C2-sep)
> $\wedge$ (C2-link)**, Condition C holds and the true cover is the unique minimum (Theorem 2) — the sim5x PSV
> ladder, where $K \geq 2$ gives $100\%$ recovery.
>
> The identifiability threshold of Theorem 2 thus refines the empirically measured **K-bound** into a
> *cross-conflict plus linkage* threshold: $K_{ij} = 0 \Rightarrow$ merge / non-identifiable;
> $K_{ij} \geq 2$ **and** every copy's identity columns linked into a single component by its reads
> $\Rightarrow$ exact, unique identification. The linkage clause is not optional: $K_{ij} \geq 2$ with *unlinked*
> columns (no read spanning two distinguishing sites) is the haplotype-phasing-ambiguous regime and admits a
> recombinant cover, exactly the counterexample above.

The $K = 0$ branch is verified directly in `check_thm2_K0_merge` (identical copies $\Rightarrow$ zero conflict
edges $\Rightarrow$ $\mathrm{MCC} = 1$, a forced merge); the $K \geq 2$ branch in `check_thm2_recovery`
(a $K_{ij} = 2$ instance satisfying the full (C2), whose unique minimum cover equals the true partition); and the
sufficiency of the corrected Condition C over the whole small-instance distribution in
`check_thm2_uniqueness_random` (which also excludes the recombinant counterexample and ablates each clause of
(C2)).

Run: `python3 bench/copy_assignment_theory_checks.py` — exits 0 and prints, among the checks:

```
OK  - Thm 2: K=2 instance under C -> unique minimum cover == true copies
OK  - Thm 2 boundary: K=0 -> minimum cover = 1 (forced merge, non-identifiable)
OK  - Thm 2 uniqueness (randomized): <S> instances under C, 0 uniqueness violations; counterexample excluded; (C2-sep),(C2-link) each load-bearing
```

---

*See also: `bench/family_definition_formal.md` for the upstream family-detection step whose output populates
$R$ for each identified gene family.*
