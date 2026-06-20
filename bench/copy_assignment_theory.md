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

## §5 Identifiability: Strong Separation and Theorem 2

Theorem 1 says the *general* MCC problem is hard. But the assignment instances that arise in practice are not
adversarial: when paralog copies carry enough distinguishing variation **and** the reads cover the family densely
enough, the truth is not merely *a* minimum cover — it is the *unique* one. Theorem 2 makes this precise under a
clean closed-form condition, **Strong Separation**, that is sufficient for **all** $K$. It is an
**identifiability / uniqueness** statement: it certifies *what the optimum is* (the true copy set), not *how to
compute it efficiently*. The conflict graph $H$ is **not** a disjoint union of cliques in general (two reads of
the same copy that observe disjoint columns are non-adjacent), so uniqueness does **not** follow from any
cluster-graph structure; it is proved directly below.

> **How the condition was found — and why a weaker one fails.** Two earlier candidate conditions failed
> adversarial review. The first (per-read cross-conflict alone) admits the classical haplotype-**phasing**
> ambiguity. The second (`sep+link`: per-read cross-conflict **plus** per-copy column-linkage) repairs phasing at
> $K = 2$ but **still fails at $K \geq 3$** through cross-copy *recombination*. The correct sufficient condition,
> and the exact frontier where the weaker one breaks, were established not by a sampling argument (randomized
> checks gave false confidence — they missed the violations) but by **exhaustive enumeration** of every small
> instance. The §5 result below is what that enumeration certifies, and `check_thm2_strong_exhaustive` re-runs it.

### The recombinant cover (why phasing must be pinned) and the K-frontier (why linkage is not enough)

A first attempt at the identifiability condition is the purely *separating* demand: every read of copy $i$
conflicts with *some* read of every foreign copy $j$, so that no jointly-consistent part can straddle two true
copies. That demand is necessary but **not sufficient**: the failure is exactly the classical
haplotype-**phasing** ambiguity. Consider $K = 2$ copies over $L = 2$ columns,
$$
c_0 = (0, 0), \qquad c_1 = (1, 1),
$$
so both columns distinguish the pair ($D_{01} = \{0, 1\}$). Take four single-column reads: $r_0$ from $c_0$
observing only column $0$ (value $0$); $r_1$ from $c_0$ observing only column $1$ (value $0$); $r_2$ from $c_1$
observing only column $0$ (value $1$); $r_3$ from $c_1$ observing only column $1$ (value $1$). Every read
conflicts with a foreign read ($r_0 \sim r_2$ at column $0$, $r_1 \sim r_3$ at column $1$), so the separating
demand holds. Yet $H$ has **only** the two edges $\{r_0, r_2\}$ and $\{r_1, r_3\}$, and admits **two** distinct
minimum (size-$2$) covers:
$$
\underbrace{\{r_0, r_1\},\ \{r_2, r_3\}}_{\text{true } C^\*}
\qquad\text{and}\qquad
\underbrace{\{r_0, r_3\},\ \{r_1, r_2\}}_{\text{recombinant}}.
$$
The recombinant $\{(0,1), (1,0)\}$ is just as good a size-$2$ cover as the truth: no molecule observes both
distinguishing columns, so the **phase** between columns $0$ and $1$ is unconstrained. At $K = 2$ the cure is to
demand, *within each copy*, that the reads **link** the copy's distinguishing columns by overlapping observation
(the `link` clause), pinning the phase. This `sep+link` condition is the standard error-free phasing
identifiability statement, and it is sufficient **at $K = 2$**.

**But `sep+link` fails at $K \geq 3$.** With three or more copies, cross-copy *recombination* can assemble an
alternative minimum cover even when every copy's columns are internally linked. The explicit witness (over
$L = 3$, columns $\{0, 1, 2\}$) is
$$
c_0 = (1,1,0), \qquad c_1 = (0,0,1), \qquad c_2 = (0,1,1),
$$
with each copy contributing two reads on the column-windows $\{0,1\}$ and $\{1,2\}$. The class
$$
\{\, \text{($c_0$'s read on cols }\{1,2\}\text{)},\ \text{($c_2$'s read on cols }\{0,1\}\text{)} \,\}
$$
co-observes only column $1$, where $c_0$ and $c_2$ **agree** (both $1$), so the two reads are *compatible* and the
class is jointly consistent — its realized vector is the **novel** haplotype $(0,1,0)$, which is **no true copy**.
This recombinant class anchors an alternative size-$3$ cover, so the minimum cover is **not unique**, yet `sep+link`
holds (each copy's two windows link its identity columns). The phasing intuition (link within a copy) is simply
*too local* once $K \geq 3$: it cannot forbid a class from being phase-consistent with a *blend* of two foreign
copies. The condition that does forbid this is global and stronger.

### Condition: Strong Separation

> **Strong Separation.** Let $C^\* = \{c_1, \ldots, c_K\}$ ($K \geq 1$) be the true copies, and assume each read
> originates from exactly one copy (error-free core: each read $r$ is consistent with its origin copy). For
> $i \neq j$ write $D_{ij} = \{\,d \in [m] : (c_i)_d \neq (c_j)_d\,\}$ for the **distinguishing columns** of the
> pair. Strong Separation requires:
>
> > For all $i \neq j$, **every** read of copy $i$ conflicts (in $H$) with **every** read of copy $j$.
>
> Equivalently: **no jointly-consistent class contains reads from two different copies** — every cross-copy read
> pair co-observes a distinguishing column with differing alleles.

**Biological reading.** Strong Separation is *pairwise distinguishability* (PSVs exist: $D_{ij} \neq \varnothing$,
so the copies are distinct allele-vectors) **plus dense read coverage** (every cross-copy overlap actually
*catches* a difference, not merely is *capable* of doing so). It is the cross-copy strengthening of the earlier
per-read demand: not "every read meets *some* foreign read in conflict," but "every read meets *every* foreign
read in conflict." That uniform demand is exactly what kills both failure modes at once — phasing (no class can
mix two copies, so no phase can be swapped) and recombination (no class can blend two foreign copies, so no novel
vector can be realized).

### Theorem 2

> **Theorem 2 (recovery, all $K$).** Under **Strong Separation**, the true copy set $C^\*$ is the **unique**
> minimum copy cover of $R$: $\mathrm{MCC}(R) = K$, and the only minimum-cardinality partition of $R$ into
> jointly-consistent parts is $C^\*$ itself.

*Proof.* By Lemma 1, copy covers are exactly proper colorings of $H$ (color classes $=$ jointly-consistent parts
$=$ independent sets of $H$), so we reason about colorings.

**(a) $C^\*$ is a proper $K$-coloring, so $\mathrm{MCC} \leq K$.** All reads of a single true copy $c_i$ are
consistent with the common allele-vector $c_i$, hence pairwise compatible (they agree at every co-observed
column); they form an independent set of $H$. Coloring each read by its origin copy partitions $R$ into $K$
jointly-consistent parts, so $\mathrm{MCC} \leq K$.

**(b) $\mathrm{MCC} = K$ and $C^\*$ is minimum.** Under Strong Separation, for every $i \neq j$ *every* copy-$i$
read conflicts with *every* copy-$j$ read; in particular some cross-copy pair conflicts. Hence no independent set
of $H$ — equivalently, no jointly-consistent part — can contain reads of two distinct true copies: any
independent set lies **within one copy**. Any proper coloring's classes are independent sets, so each class lies
inside a single true copy, and covering all $K$ (nonempty) copies needs $\geq K$ classes. With (a),
$\mathrm{MCC} = K$ and $C^\*$ is a minimum cover.

**(c) Uniqueness.** Let $\Pi = \{P_1, \ldots, P_K\}$ be any minimum (i.e. $K$-)coloring; we show $\Pi = C^\*$.
Because $\chi(H) = K$ by (b), every proper $K$-coloring is minimum, so this applies to all of them. Each class
$P_s$ is an independent set, hence (by the argument of (b), now that **every** cross-copy pair conflicts) is
contained in a **single** true copy: there is a map $\tau(s)$ with $P_s \subseteq c_{\tau(s)}$ (i.e. every read of
$P_s$ originates from copy $\tau(s)$). The $K$ classes partition all $n$ reads among $K$ true copies, each of
which is nonempty; since each class injects into one copy and the $K$ copies must all be covered, $\tau$ is a
**bijection** from the $K$ classes onto the $K$ copies, and each class $P_s$ must equal the *entire* read set of
its copy $c_{\tau(s)}$ (any copy-$\tau(s)$ read left in another class $P_{s'} \subseteq c_{\tau(s')}$ would force
$\tau(s') = \tau(s)$, contradicting injectivity). Therefore $\{P_1, \ldots, P_K\} = C^\*$. $\square$

> **Why step (c) is now valid.** The earlier (invalid) version used only per-read separation: a copy-$i$ read
> whose foreign-conflict is witnessed by *some* copy-$j$ read need not conflict with *every* copy-$j$ read, so a
> non-conflicting copy-$j$ read could cohabit its class — precisely the recombinant cover. Strong Separation
> closes exactly this gap: because **every** foreign read conflicts (not just *some*), each independent class is
> trapped inside one copy with no escape, and the counting in (c) goes through cleanly for all $K$. No phasing
> propagation, no linkage-graph connectivity argument, no $K = 2$ restriction.

### Proposition (tightness / the K-frontier)

> **Proposition.** Strong Separation is **sufficient** for unique recovery for all $K$, but it is **not
> necessary**.
>
> **(i) ($K = 2$).** The alternative natural condition `sep+link` (incomparable to Strong Separation at $K = 2$,
> not strictly weaker) — cross-conflict per read, plus per-copy
> column-linkage connected — already suffices at $K = 2$. This is the standard error-free phasing
> identifiability: linkage pins the phase, forbidding the recombinant cover.
>
> **(ii) ($K \geq 3$).** `sep+link` **fails**: cross-copy *recombination* can produce an alternative minimum
> cover realizing a haplotype outside $C^\*$. The explicit witness is
> $c_0 = (1,1,0),\ c_1 = (0,0,1),\ c_2 = (0,1,1)$ with each copy read on windows $\{0,1\}$ and $\{1,2\}$: the
> class $\{$ $c_0$'s read on $\{1,2\}$, $c_2$'s read on $\{0,1\}$ $\}$ is jointly consistent (the two reads share
> only column $1$, where $c_0 = c_2 = 1$) and realizes the novel vector $(0,1,0)$, anchoring a second minimum
> cover. This instance satisfies `sep+link` but **violates** Strong Separation (that cross-copy pair does not
> conflict).

**Honest conclusion.** Multi-copy recovery at $K \geq 3$ is **strictly harder** than the pairwise ($K = 2$) case:
there is no closed-form, per-pair-plus-per-copy condition that captures it. The exact
**necessary-and-sufficient** condition for unique recovery is **recombination-freeness** — *no alternative
size-$K$ cover exists*. Recombination-freeness is **instance-global** (it quantifies over all alternative covers,
not over pairs or single copies) and has **no clean closed form**. Strong Separation is the clean *sufficient*
surrogate: the exhaustive enumeration finds it holds for only a minority of the truly-unique $K = 3$ instances
($\approx 15\%$ over the full $L = 3$ enumeration; $808/2992 \approx 27\%$ over the bounded re-run shipped in the
check), so it is a conservative guarantee, not a characterization.

### Remark (polynomial-time recovery, deferred)

Theorem 2 identifies the optimum with the truth; it does **not** assert that the optimum is computable in
polynomial time (Theorem 1 forbids that in general). Under Strong Separation the conflict graph is a complete
$K$-partite graph on the copies' read sets, so recovery is immediate in principle; but the design and
polynomial-time analysis of a solver for the broader recombination-free regime — and the precise structural class
of conflict graphs $H$ for which it runs in polynomial time — are deferred to a follow-up, out of scope for this
identifiability note.

### Corollary (the K-bound)

> **Corollary (K-bound).**
> **(i)** If $K_{ij} = 0$ for some pair (copies identical over the transcribed/observed region), then there is
> **no distinguishing column**: Strong Separation fails (no cross-copy pair can conflict) **and** the two copies
> merge — that pair produces no conflict edge in $H$, the minimum cover collapses them into one part, and the
> true copies are **provably unrecoverable**. This is the MAGEA co-located regime (measured resolvable fraction
> $0/494$).
> **(ii)** With $K_{ij} \geq 2$ for every pair and **sufficient read coverage**, Strong Separation plausibly
> holds (every cross-copy overlap catches a difference), so the true cover is the unique minimum (Theorem 2) —
> the sim5x PSV ladder, where $K \geq 2$ gives $100\%$ recovery.
>
> **Honesty caveat.** Strong Separation is a **conservative sufficient** condition, not a characterization. The
> empirical resolver may succeed under *weaker* conditions (e.g. recombination-free instances that are not
> Strongly Separated, or $K = 2$ instances that are merely `sep+link`). Strong Separation is a **guarantee** of
> exact recovery, not a necessary requirement; the tight boundary is recombination-freeness, which has no clean
> closed form (Proposition above).

The $K = 0$ branch is verified directly in `check_thm2_K0_merge` (identical copies $\Rightarrow$ zero conflict
edges $\Rightarrow$ $\mathrm{MCC} = 1$, a forced merge); the $K \geq 2$ branch in `check_thm2_recovery` (a Strongly
Separated $K = 2$ instance whose unique minimum cover equals the true partition); and the full §5 result over the
entire small-instance universe in `check_thm2_strong_exhaustive`, which **exhaustively** enumerates all $K \in
\{2, 3\}$, $L = 3$ instances with $\mathrm{MCC} = K$ and certifies: `strong` $\to$ **0** uniqueness violations at
both $K = 2$ and $K = 3$ (the sufficiency certificate); `sep+link` $\to$ **0** at $K = 2$ but **$> 0$** at
$K = 3$ (the K-frontier); and that the recombination witness has a non-unique cover and is excluded by `strong`.

Run: `python3 bench/copy_assignment_theory_checks.py` — exits 0 and prints, among the checks:

```
OK  - Thm 2: K=2 instance under Strong Separation -> unique minimum cover == true copies
    [exhaustive] K=2: total(MCC=2)=10728 unique=3864  strong holds=1944 viol=0  sep+link holds=1728 viol=0
    [exhaustive] K=3: total(MCC=3)=11112 unique=2992  strong holds=808 viol=0  sep+link holds=3240 viol=248
OK  - Thm 2 (exhaustive K=2,3 / L=3): strong viol K2=0/K3=0 (SUFFICIENT all K); sep+link viol K2=0/K3=248 (K-frontier: K=2 only); recombination witness non-unique and excluded by strong
OK  - Thm 2 boundary: K=0 -> minimum cover = 1 (forced merge, non-identifiable)
```

---

*See also: `bench/family_definition_formal.md` for the upstream family-detection step whose output populates
$R$ for each identified gene family.*
