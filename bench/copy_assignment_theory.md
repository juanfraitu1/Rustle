# Copy-assignment theory: identifiability, hardness, and the K-frontier

*This note develops the combinatorial foundation for the read-to-copy assignment problem in multi-copy gene
families. §1 situates the problem and its place in the multimapping-resolution lineage. §2 fixes the model and
notation. §3 (Lemma 1) proves MCC = χ(H). §4 (Theorem 1) establishes NP-hardness via coloring reduction. §5
(Theorem 2) proves unique recovery under Strong Separation; the Proposition identifies the K-frontier where
pairwise conditions fail for K ≥ 3; Theorem 3 gives a polynomial-time `RECOVER` algorithm under the same
condition. §6 lifts the framework to joint copy + isoform recovery. §7 cites empirical corroboration. §8
discusses the dichotomy as a through-line across the three research interests.*

---

## §1 Introduction

### Context and motivation

Long-read RNA sequencing (HiFi IsoSeq, Oxford Nanopore) resolves full isoform structures in a single molecule,
but reads from paralogous gene-family members are often **alignment-ambiguous**: a single molecule may map with
equal or near-equal quality to two or more distinct genomic loci. The downstream assembly then faces a
**copy-assignment** problem: given a set of reads and their observed allelic evidence, partition the reads so that
each part is explained by a single gene-copy (allele-vector).

Copy-assignment sits at the intersection of three active research interests:

1. **Family detection** (interest #1): which expressed loci form a confusion group? A multi-copy gene family
   is defined as a connected component of the read-conflict graph — the maximal set of loci among which reads
   are genuinely confused (see `bench/family_definition_formal.md`). Copy-assignment is the unit of work that
   is performed *inside* each family component output by detection; it is logically downstream of interest #1.

2. **Copy assignment under ambiguity** (interest #2, this note): given a confirmed family component, how should
   reads be partitioned across copies? When is the truth the *unique* minimum-cost solution, and when is
   attribution provably arbitrary?

3. **Allele-specific junctions** (interest #3, §6 and §8): long-read molecules that span both a PSV column and
   a splice junction simultaneously carry *joint* allele+isoform evidence. The junction-as-column lift (§6)
   shows that the copy-assignment framework extends directly to recover isoform structure jointly with copy
   identity — the substrate for downstream allele-specific junction analysis.

### The multimapping-resolution lineage

The copy-assignment problem has a precise lineage in the multimapping literature. Canzar et al. (2016) frame
multi-read (multimapping) conflict as a **Maximum Facility Location** problem: each read is a client that must
be assigned to one of its feasible placements (facilities), and the objective is to maximize the number of
reads assigned consistently with their neighbours. Their LP-rounding algorithm achieves a provable
$(1 - 1/e)$-approximation with a $0.19$-approximation guarantee, and they establish that the unweighted version
is NP-hard to approximate beyond a constant factor.

The present note operates in a complementary regime. Rather than the *weighted* conflict-resolution objective
(maximize consistent assignments), we study the **parsimony** objective: minimize the number of distinct copies
needed to explain all reads without error. The two objectives are hard for different reasons — conflict
resolution is hard because the assignment landscape is rough, parsimony is hard because minimum coloring of an
arbitrary conflict graph is NP-complete (Theorem 1 below). Both hardness results are **robust across objective
functions**, making the general problem doubly intractable.

The Minimum Error Correction (MEC) problem [Lippert et al. 2002] is the standard genomic haplotype-assembly
variant: fix $k$ haplotypes and minimize allele flip count. MEC is NP-hard already at $k = 2$. Our MCC
(Minimum Copy Cover) is the parsimony analogue for the multi-copy ($k$-copy) RNA-family setting.

### This note's contribution

This note develops the combinatorial foundation for copy-assignment in three parts:

- **Model and MCC = χ(H) (§2–§3).** We fix the formal model (reads as partial allele functions over PSV
  columns, copies as consistent allele-vectors, conflict as pairwise disagreement). Lemma 1 proves that the
  Minimum Copy Cover equals the chromatic number of the conflict graph — a tight combinatorial identity that
  makes the assignment problem graph-theoretically natural and connects it to coloring complexity.

- **Hardness/recovery dichotomy (§4–§5).** Theorem 1 (hardness), Theorem 2 (recovery), and Theorem 3
  (polynomial-time algorithm) form the central dichotomy. Theorem 1 shows MCC is NP-hard via a polynomial-time
  reduction from graph $k$-colorability. Theorem 2 shows that under **Strong Separation** — a uniform pairwise
  distinguishability condition — the true copy set is the *unique* minimum cover. Theorem 3 shows that under
  the same condition the truth is also *computable in polynomial time* via the `RECOVER` algorithm ($O(n^2 m)$,
  connected components of the compatibility graph) — so the dichotomy is closed for the *combinatorial* recovery
  problem: NP-hard in general, unique optimum and efficient algorithm under Strong Separation. (The shipped
  pipeline does not run `RECOVER`; **Theorem 4 (§5b, the bridge)** is the production-side certifier — a sound
  per-read identifiability certifier for all K≥1.)

- **K-frontier finding (§5 Proposition).** Strong Separation is sufficient for all $K$, but the weaker
  `sep+link` condition (per-read cross-conflict plus per-copy column-linkage) already suffices at $K = 2$ and
  **fails at $K \geq 3$** through cross-copy recombination. Multi-copy recovery is strictly harder than the
  pairwise case: no per-pair-plus-per-copy condition captures it. This is not a gap in the theory — it is a
  genuine structural finding, certified by exhaustive enumeration over all small instances (§5 checks).

- **Paths/isoforms (§6).** The junction-as-column lift extends the model to joint (allele + isoform) recovery,
  directly connecting the copy-assignment identifiability result to the allele-specific junction substrate
  (interest #3).

The recovery condition's threshold — positive distinguishing PSV count ($K_{ij} \geq 1$, plausibly $K_{ij}
\geq 2$ for dense coverage) — matches the empirically measured K-bound from the sim5x PSV ladder and the GGO
silver-standard panel (§7). This correspondence between the theoretical sufficient condition and the observed
empirical boundary is the main finding: the K-quantity simultaneously governs detection (interest #1),
assignment (interest #2), and the utility of allele-specific junction evidence (interest #3).

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
with opposite alleles: $r_w(e_{w,w'}) = 0$ and $r_w(e_{w',w}) = 1$ (where $w' = w-1 \bmod 5$).
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

**Theorem 1.** MCC is NP-hard. (Hence the Constrained Minimum Copy-Path-Cover (CMCPC) of §6 — which contains MCC as the single-backbone, allele-only case — is NP-hard.)

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

#### Definition (MEC, as implemented in `vg_family/phasing.rs`)

Encode the reads over the $m$ heterozygous sites of a locus as an allele matrix

$$M \in \{0,1,\bot\}^{n\times m}, \qquad
M[r][j]=\begin{cases}\text{allele of read } r \text{ at site } j\\[2pt] \bot & \text{if } r \text{ does not cover } j,\end{cases}$$

with $n$ reads (rows) and $m$ het sites (columns). The **Minimum Error Correction** value is

$$\mathrm{MEC}(M)=\min_{\substack{h_A,h_B\in\{0,1\}^m\\ \sigma:[n]\to\{A,B\}}}\
\sum_{r=1}^{n}\ \sum_{j:\,M[r][j]\neq\bot}\ \mathbb{1}\!\big[\,M[r][j]\neq h_{\sigma(r)}[j]\,\big],$$

i.e. the minimum number of **(read, site) cells** that must be flipped so that every read agrees
perfectly, over its covered sites, with its assigned haplotype. NP-hard already at $k=2$
[Lippert et al. 2002; Cilibrasi et al. 2005].

Two structural commitments make the implementation **exact** (not heuristic):

1. **Diploid complementarity.** At a *heterozygous* site the two haplotypes carry opposite alleles, so
   the code fixes $h_B=\lnot h_A$ (`mec_brute`). This restricts the search to the $2^m$ choices of
   $h_A$, not $2^{2m}$; and for fixed $h_A$ each read *independently* joins the cheaper of
   $\{h_A,\lnot h_A\}$, so $h_A$ determines the optimal $\sigma$. Hence
   $\mathrm{MEC}(M)=\min_{h_A}\sum_r \min\!\big(d_H(M[r],h_A),\,d_H(M[r],\lnot h_A)\big)$, with $d_H$
   the covered-cell Hamming distance.

2. **Pattern collapse ⇒ coverage-independent exactness.** Full-length RNA reads all span the gene, so a
   naive column-sweep DP carries state $2^{\text{boundary coverage}}=2^{\Theta(n)}$. But two reads with
   an identical allele row are WLOG on the same haplotype (identical per-column contribution ⇒ any split
   merges to the cheaper side without raising cost), so identical rows are collapsed into **weighted
   classes** and the DP runs over classes. A clean diploid locus has $\approx 2$ classes (the two
   haplotypes) plus a few error-singletons, so the open state is $O(1)$ and the DP is exact at any
   coverage. The DP is proven equivalent to `mec_brute`; the unit test asserts $\text{cost}=$ number of
   error cells. The `max_coverage` cap bounds *distinct classes per column*, not raw reads, and rarely
   fires on real diploid data.

**Scope.** MEC here is the **within-copy diploid het phasing** run by `--phase` (one copy at a time,
$k=2$). It is distinct from the family-level copy-assignment objective $\mathrm{MCC}=\chi(H)$
(§3–§5): MCC minimizes the *number of copies* (parsimony/coloring) over the whole family variation
graph, whereas MEC minimizes *allele flips* against the two complementary haplotypes of a single copy.
The note's point is that **both** objectives are NP-hard, so intractability is objective-robust.

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

*Proof.* The proof of Theorem 2 establishes that Strong Separation makes $H$ the **complete $K$-partite** graph
whose parts are the true read-classes $P_1, \ldots, P_K$. The complement of a complete $K$-partite graph is the
**disjoint union of $K$ cliques**, one per part: within a part $P_i$ there are no $H$-edges, so every pair is a
$\bar H$-edge (a clique); across parts every pair is an $H$-edge, so no $\bar H$-edge. Hence the connected
components of $\bar H$ are exactly $P_1, \ldots, P_K$. Each $P_i$ is internally compatible, so its reads agree
at every column any of them observes; the union of their observations is therefore a single allele-vector,
consistent with every read of $P_i$ — i.e. the true copy $c_i$. Building $\bar H$ tests each of the
$\binom{n}{2}$ read pairs over at most $m$ co-observed columns ($O(n^2 m)$); connected components and the
per-component union are dominated by this. $\square$

**Self-certifying recovery.** `RECOVER` verifies, as part of the same $O(n^2 m)$ computation (the incremental
clique-check on the already-built $\bar H$ costs only an additional $O(n^2)$), whether $\bar H$ is a **disjoint
union of cliques** — equivalently, whether the compatible relation is transitive. This property certifies that
`RECOVER`'s partition is the **unique minimum copy cover of the reads**: no same-size partition exists.
Strong Separation is a *sufficient condition* for this property (Strong Sep $\Rightarrow$ $\bar H$ is a disjoint
union of cliques, by the proof above), but the property is weaker than Strong Separation — it can hold on
instances Strong Separation excludes. When the check **fails** — the $K \geq 3$ cross-copy recombination regime,
where some connected component of $\bar H$ is connected but *not* a clique — `RECOVER` refuses, reporting *not
identifiable*; the minimum cover is not unique in that case, and no partition can be guaranteed.

A subtlety worth noting: at $K = 0$ (copies identical over all observed columns), every read is compatible with
every other ($H$ has no edges), so $\bar H$ is the complete graph — a single clique, hence trivially a disjoint
union of cliques. The certifier **passes**, and `RECOVER` returns one component (merging the indistinguishable
copies). This is the correct unique minimum cover of the observed data; it reflects an information-theoretic
floor (the copies cannot be told apart), not an algorithm failure. The certifier therefore guarantees
*unique minimum cover of the reads*; Strong Separation upgrades that further to *equals the true copies* (since
under Strong Separation the components are forced to coincide with the true read-classes, Theorem 2).

**Verification.** `check_thm3_recovery_algorithm` enumerates every Strong-Separation instance over $K \in \{2,3\}$,
$L = 3$ and confirms `RECOVER` returns the true partition on **all** of them, the self-certifier accepts each, and
the recombination witness (Strong Separation fails) is **rejected** — exhaustive, not sampled.

### Remark (recovery)

**Remark (recovery).** Polynomial-time recovery under Strong Separation is given by **Theorem 3 above**
(connected components of the compatibility graph). The general problem remains NP-hard (Theorem 1); the
dichotomy is therefore closed for the *combinatorial* recovery problem — hard in general, efficiently solvable
under the identifiability condition. The shipped pipeline does not run `RECOVER`; **Theorem 4 (§5b, the
bridge)** is the production-side certifier: a sound per-read identifiability certifier for all K≥1.

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

### §5b The bridge: the production gate's `min_p` is a sound per-read identifiability certificate (Theorem 4)

Theorems 2–3 concern *family-level* recovery of the whole copy set (the minimum cover / `RECOVER`), NP-hard in
general (Theorem 1). The shipped pipeline does **not** run `RECOVER`; it runs a *per-read statistical* gate
(`copy_assign.rs`): given the copy set $C$ (from the family layer), each read $r$ is **assigned** to a copy,
**ambiguous**, or certified **tied** by an IsoCon-style significance test with a family-wide Bonferroni threshold
$\alpha/(n-1)$. Its identifiability bound is
$$\mathrm{min\_p}(r)=\max_{c\neq b}\ \prod_{j\in\mathrm{obs}(r)\cap D_{bc}}\varepsilon_j,$$
where $b$ is $r$'s MLE copy, $D_{bc}$ the distinguishing columns of $b,c$, and $\varepsilon_j$ the per-column
error proxy; $r$ is certified **tied** iff $\mathrm{min\_p}(r)\ge\alpha/(n-1)$. Until now the theory touched the
gate only at the degenerate $K=0$ vertex ($\mathrm{min\_p}=1$). Theorem 4 connects $\mathrm{min\_p}$ to the
combinatorics for **all** $K\ge1$, making the running gate load-bearing. Write
$$\delta(r)=\min_{c\neq b}\ \bigl\lvert\mathrm{obs}(r)\cap D_{bc}\bigr\rvert$$
for the number of distinguishing columns $r$ spans against its **closest** competitor copy.

> **Theorem 4 (Bridge).** In the error-free core (each read consistent with its origin copy) **and under the
> completeness precondition $\mathrm{origin}(r)\in C$** (the given copy set contains $r$'s true origin copy —
> equivalently, $C$ passes the §5 self-certifier so no copy is missing), for every read $r$:
> **(i)** with uniform $\varepsilon_j=\varepsilon<1$, $\ \mathrm{min\_p}(r)=\varepsilon^{\,\delta(r)}$ (this holds
> unconditionally; only (ii) needs the precondition).
> **(ii) (soundness)** $\delta(r)\ge1\iff b$ is the **unique** copy in $C$ consistent with $r$; then, *because
> $\mathrm{origin}(r)\in C$*, that unique copy $b$ is $r$'s origin, so the gate's argmax assignment is **correct** —
> the gate never assigns a wrong copy. **The precondition is necessary and tight:** if $\mathrm{origin}(r)\notin C$
> (the reference-absent / O4 regime), a *partial* read can be uniquely consistent with a wrong $b\in C$ ($\delta\ge1$,
> $\mathrm{min\_p}<1$) and be silently, confidently misassigned (machine-checked B6: 2,616 witnesses). This is exactly
> why the O4 two-stage freeze **abstains** at stage 1 and re-threads the abstain pool against $C\cup\{$absent copies$\}$;
> a read that observes *all* columns of an absent origin is consistent with **no** copy in $C$ and is flagged novel —
> only partial reads are at risk.
> **(iii) (abstention $=$ ambiguity)** $\delta(r)=0\iff\ \ge2$ copies are consistent with $r$ ($r$ is genuinely
> unassignable), and $\mathrm{min\_p}(r)=1$ certifies it **tied**.
> **(iv) (boundary recovery)** if $K_{ij}=0$ for some pair (copies identical over the observed columns — the
> Theorem-2 unrecoverable case), then **every** read of that pair has $\delta=0$, $\mathrm{min\_p}=1$, certified
> **tied**. The gate thus recovers exactly the Theorem-2 non-identifiable mass, **per read, in
> $O(\text{reads}\cdot\text{copies}\cdot\text{columns})$**, without computing the (NP-hard) minimum cover.

*Proof.* **(i)** $\prod_{j\in\mathrm{obs}(r)\cap D_{bc}}\varepsilon=\varepsilon^{\lvert\mathrm{obs}(r)\cap D_{bc}\rvert}$;
as $\varepsilon<1$ the max over $c$ is at the smallest exponent $\min_c\lvert\mathrm{obs}(r)\cap D_{bc}\rvert=\delta(r)$.
**(ii)** ($\Leftarrow$) if $b$ is the only consistent copy, each $c\neq b$ has some $j\in\mathrm{obs}(r)$ with
$r(j)\neq c_j$; since $r(j)=b_j$, $b_j\neq c_j$, so $j\in\mathrm{obs}(r)\cap D_{bc}$ and $\delta(r)\ge1$.
($\Rightarrow$, contrapositive) if some $c\neq b$ is also consistent, then for all $j\in\mathrm{obs}(r)$,
$r(j)=b_j=c_j$, so $\mathrm{obs}(r)\cap D_{bc}=\varnothing$ and $\delta(r)=0$. Finally, **by the precondition
$o=\mathrm{origin}(r)\in C$**, and $r$ is consistent with $o$, so $o$ is in the consistent set; uniqueness forces
$o=b$ (the error-free MLE picks a consistent copy). *(Without the precondition this last step fails — $o\notin C$
means the unique consistent copy in $C$ need not be the origin; B6 exhibits the misassignment.)*
**(iii)** $r$ always has $\ge1$ consistent copy (its origin), so $\delta(r)=0\iff$ not-unique $\iff\ge2$ consistent;
$\mathrm{min\_p}=\varepsilon^0=1\ge\alpha/(n-1)$. **(iv)** $K_{ij}=0\Rightarrow D_{ij}=\varnothing\Rightarrow c_i=c_j$;
a copy-$i$ read is consistent with both, so $\ge2$ consistent and (iii) applies; only the per-read scan is used. $\square$

> **What it buys, and what it does not.** Theorem 4 makes the *running* gate load-bearing: via $\mathrm{min\_p}$
> it is a **sound** per-read certifier of the combinatorial identifiability condition (unique consistency) for
> **all** $K\ge1$ *given a complete $C$* — *assigned* $\Rightarrow$ combinatorially determined (no guess);
> *tied at $\delta=0$* $\Rightarrow$ genuinely ambiguous; the Theorem-2 $K=0$ boundary recovered per read without
> `RECOVER`. It is **per-read** (it takes $C$ as input; it is not family recovery, and does not resolve the
> $K\ge3$ recombination non-uniqueness, which concerns the *cover*, not a single read's assignment — B7
> exhibits a read set with several distinct minimum covers on which the per-read gate is nonetheless sound
> *within each fixed $C$*, confirming Theorem 4's scope is per-read-given-$C$, not cover recovery). Its two
> standing hypotheses are the **error-free core** and **completeness $\mathrm{origin}(r)\in C$**; the latter is
> discharged operationally by the O4 abstain-and-re-thread stage, not assumed away. The statistical gate adds error tolerance: with
> threshold $\tau=\lfloor\log(\alpha/(n-1))/\log\varepsilon\rfloor$ it abstains whenever $\delta(r)\le\tau$, so it
> is sound but **conservative** — it may tie a combinatorially-determined read whose evidence is below significance
> $\alpha$ (e.g. a lone distinguishing column at $\alpha=10^{-3}$). Soundness (never misassign) is the guarantee;
> completeness is bounded by $\alpha$.

**Verification.** `check_bridge_theorem` (`bench/bridge_theorem_check.py`, integrated into the suite)
**exhaustively** enumerates all copy sets and reads for $K\le3$, $m\le3$, $\lvert A\rvert\le3$ (3,400 copy-sets /
67,320 reads) plus the full $K_{ij}=0$ boundary universe, certifying (i)–(iv): $\mathrm{min\_p}=\varepsilon^{\delta}$;
$\delta\ge1\Leftrightarrow$ unique-consistent; soundness (assignment $=$ origin); $\delta=0\Leftrightarrow$ ambiguous;
and $K_{ij}=0\Rightarrow$ all-tied. Two further checks pin the theorem's hypotheses: **B6** (precondition
necessity) drops $\mathrm{origin}(r)\in C$ and finds 2,616 instances where the gate confidently *misassigns* a
read whose origin is absent — and certifies the escape, that a full-column read of an absent origin is consistent
with no copy in $C$ (flagged novel); **B7** (recombinant cover) exhibits a read set with several distinct minimum
covers of size $\ge3$ (the $K\ge3$ recombination obstruction) on which the per-read gate stays sound within each
fixed $C$, confirming Theorem 4 is per-read-given-$C$, not cover recovery.

---

---

## §6 The paths/isoforms corollary: joint recovery of copy number and isoform structure

Sections §2–§5 treat a copy as a bare allele-vector over PSV columns.  Here we show that isoform structure
lifts into the same framework for free, by treating each distinguishing junction as an additional column.  The
core definitions, Lemma 1, and Theorem 2 then apply verbatim over the enlarged column set, and the
minimum constrained path-cover of the family variation graph recovers copies **and** isoforms jointly.

### Setup: copies as (path, allele-vector) pairs

A **gene copy with isoform structure** is a pair $(\pi, a)$ where $\pi$ is a path through the family
variation-graph DAG (determining which splice junctions are used) and $a$ is an allele-vector over the PSV
columns.  Two copies $(\pi, a)$ and $(\pi', a')$ are distinguishable if they differ at any PSV column **or** at
any junction.

### The junction-as-column lift

Model each isoform-distinguishing junction $J$ as an additional binary **column** $c_J$, with allele alphabet
$\{0, 1\}$: a read that spans junction $J$ observes $c_J$ and reports which splice branch it took (0 for the
exon-skipping branch, 1 for the inclusion branch, or analogously for donor/acceptor alternatives).

Under this encoding:

- A read is a partial function from $[m] \cup \{c_J : J \text{ a junction column}\}$ to alleles, exactly as
  before.
- Two reads **conflict** if they co-observe any column (PSV or junction) with differing alleles — the conflict
  definition is unchanged.
- A copy is a consistent allele-**and**-junction vector: an element of $\prod_j A_j \times \{0,1\}^{\#\text{junctions}}$.
- The **conflict graph** $H$, the **MCC**, and the **jointly-consistent** / **independent-set** equivalence all
  carry over verbatim; the proof of Lemma 1 never refers to the semantics of columns.

### Corollary (joint copies and isoforms)

> **Corollary (joint copies + isoforms).** Re-attach isoform structure: a copy is $(\pi, a)$ and reads carry
> junction observations.  Treat each distinguishing junction as an additional linkage column over the branch
> alphabet; then a copy is a consistent allele-and-junction vector, reads conflict on disagreeing co-observed
> alleles **or** junctions, and Lemma 1 and Theorem 2 apply verbatim over the combined column set.
>
> Hence under **Strong Separation over the combined (allele+junction) columns** — every cross-copy read pair
> co-observes at least one distinguishing column (PSV or junction) and disagrees on it — the minimum
> **constrained path-cover** of the family variation graph is unique and equals the true copies-with-isoforms,
> recovering $\#$copies **and** their isoform structure jointly.  Unconstrained minimum path-cover (max-flow on
> the DAG) is the relaxation; allele+junction linkage promotes it from *a* minimum cover to *the* true cover
> under Strong Separation.

### The spanning-read requirement and the honest caveat

Strong Separation over the combined column set places a concrete demand on the reads: **every cross-copy read
pair must share at least one co-observed column**.  A read observing only a single PSV column and a foreign read
observing only a junction column share no column at all — they are trivially compatible and cannot contribute a
conflict edge in $H$.  This is not a gap in the theory; it is the exact formalization of the biological
intuition that **long-read linkage is what makes the corollary operational**: a molecule that spans both an
exon carrying a PSV and the downstream splice site observes both $c_J$ and the PSV column, so it will conflict
with every foreign read that also spans either feature.

Therefore Strong Separation over the combined column set holds **if and only if** the read coverage is dense
enough that every cross-copy pair is connected by at least one co-observed distinguishing column — precisely the
condition that long-read (HiFi/ONT) spanning molecules satisfy in practice.

**Honest caveat.** If two copies of the same gene share identical PSV profiles and are distinguished **only** by
a junction that no read spans (e.g., an ultra-long intron whose two splice variants are never jointly observed in
a single molecule), then the junction column is unlinked: it contributes no conflict edge, Strong Separation
fails for that junction, and the two copies are indistinguishable to the cover — exactly the $K_{ij} = 0$ case
of the K-bound corollary (§5), now instantiated at the junction axis.  The core result stands; that degenerate
case is future work, identical in character to the unresolvable-PSV regime.

### Verification

The check `check_corollary_paths` in `bench/copy_assignment_theory_checks.py` instantiates the corollary on the
minimal spanning-read case:

- **Two copies** differing at PSV column $0$ (copy A allele $0$, copy B allele $1$) **and** at junction
  pseudo-column $99$ (copy A branch $0$, copy B branch $1$).
- **Reads** each span **both** columns (PSV + junction), as a long molecule would.
- Every cross-copy pair co-observes both distinguishing columns and disagrees on at least one: Strong
  Separation holds, Theorem 2 applies, and `mcc_bruteforce` returns $2 = $ true copy count.

```
OK  - Corollary: junction treated as a linkage column -> path-cover recovers copies+isoforms (MCC=2)
```

*See also: `bench/family_definition_formal.md` for the upstream family-detection step whose output populates
$R$ for each identified gene family.*

---

## §6b Tier-3: co-quantification of the irreducible core

For copies identical over the entire transcribed region (the K=0 / Strong-Separation-fails regime), per-read
assignment is impossible (no distinguishing column). The natural fallback is per-copy *quantification* — but the
same identical-sequence condition makes that unidentifiable too, in a precise sense.

> **Proposition (Tier-3 unidentifiability).** Let a family have copies $c_1,\ldots,c_K$ identical over the
> transcribed region, and let reads $R$ ($|R| = N$) each be consistent with every copy. Under the mixture
> likelihood $L(a) = \prod_{r\in R}\sum_{k=1}^{K}\frac{a_k}{N}\,P(r\mid c_k)$ with per-copy abundances
> $a=(a_1,\ldots,a_K)$, $a_k\ge 0$, $\sum_k a_k = N$, the likelihood $L(a)$ is **constant** on the simplex
> $\{a:\sum_k a_k=N\}$: the per-copy split is statistically **unidentifiable** from RNA. The aggregate $N=|R|$
> is identifiable (a sufficient statistic). Under a copy-number / dosage prior $\pi(a)$, the MAP estimate is
> well-posed and equals the mode of $\pi$ scaled to $N$ — RNA contributes nothing to the per-copy direction.

*Proof.* The copies are identical, so $P(r\mid c_k)=p_r$ is the same for every $k$. Each read's mixture term is
$\sum_k\frac{a_k}{N}p_r = p_r\frac{\sum_k a_k}{N}=p_r$, independent of $a$. Hence $L(a)=\prod_r p_r$ is constant
in $a$; it factors through $N$ alone, so $N$ is identifiable and $a$ is not. With prior $\pi$, the posterior is
$\propto \pi(a)\,L(a)\propto\pi(a)$ on the simplex, so the MAP is the mode of $\pi$ (scaled to $N$). $\square$

This is the identifiability theorem's K=0 floor in the quantification frame: the same emptiness of the
distinguishing-column set that forbids per-read assignment (Theorem 2 corollary) flattens the per-copy
likelihood. The honest Tier-3 output is therefore the **family aggregate** (identifiable) plus a
prior-conditioned per-copy split whose uncertainty set is the entire simplex. (Read error perturbs each $p_r$
but not the $a$-independence of the mixture term, so the unidentifiability is robust to copy-independent error.)
`check_tier3_coquant_unidentifiable` is the executable witness: every apportionment of a fixed $N$ yields
identical likelihood.

---

## §7 Empirical corroboration

This section cites existing experimental results that illustrate the theoretical dichotomy. It does **not**
claim the experiments prove the theorems; they are independent measurements that are consistent with the theory
and illuminate where in the K-landscape each regime lies.

### 7.1 Theorem 2 / K ≥ 2: sim5x K-ladder and GGO silver standard

**Sim5x K-ladder.** The synthetic benchmark (`bench/sim_reads.py`,
`/home/juanfra/winloci_scratch/sim5x/`) constructs five near-identical tandem copies of a target gene and
sweeps the PSV count K across $\{0, 1, 2, 3, 4, 5, 6, 7, 8\}$ by introducing controlled SNVs at paralog-
specific variant sites. The measured copy-assignment accuracy is:

- $K = 0$: 0% correct assignment (MAPQ 0 throughout; attribution is intrinsically arbitrary — this is the
  forced-merge regime of the K-bound corollary).
- $K = 1$: partial recovery (borderline; single-PSV coverage is often insufficient to achieve Strong
  Separation on every cross-copy read pair at realistic coverage depths).
- $K \geq 2$: **100% correct assignment** across all tested coverage depths.

The empirical threshold $K \geq 2$ is consistent with Theorem 2: with two or more distinguishing columns per
copy pair and sufficient read coverage, Strong Separation plausibly holds and the unique minimum cover equals
the truth.

**GGO silver standard.** On the real GGO HiFi IsoSeq panel, the `detect_and_assign` MAGEA smoke run reports
**1026/1026 (100%)** reads correctly assigned in the K ≥ 2 families (RABL2, AK6, CCDC196). These families are
all in the well-separated regime (MAPQ > 0, divergence-gap decisive), where Strong Separation holds by the
aligner's own evidence; the coverage condition is met at ≥ 47 de-conflict reads per family.

> ⚠ **What this figure is and is NOT.** It is a *consistency check in the EASY (MAPQ > 0) regime* — and
> the "silver" truth is the aligner's own primary placement, so the metric is **circular by construction**
> (it confirms the resolver agrees with minimap2 exactly where minimap2 was already confident). It says
> nothing about the hard MAPQ-0 regime the thesis is actually about, and 100% here is expected, not
> impressive. The **load-bearing identifiability evidence is the sim5x labeled-truth K-ladder above**
> (0%@K=0 → 100%@K≥2 against *planted* per-read copy labels, not alignment), which has a non-circular
> ground truth. Cite the ladder, not the 1026/1026, as the empirical spine of Theorem 2.

**Honest caveat.** Strong Separation is a conservative *sufficient* condition. The empirical resolver may
succeed in K ≥ 2 instances that are not Strongly Separated but are recombination-free (the tighter necessary-
and-sufficient condition of the Proposition). The 100% recovery figure is consistent with Strong Separation
holding on these instances; it is not a proof that the condition is tight.

### 7.2 Theorem 2 boundary / K = 0: MAGEA co-located arrays

The K = 0 case is directly measured in `bench/resolution_improvement_bound.md` on the MAGEA co-located tandem
arrays:

- **494 ambiguous reads** across the MAGEA co-located pairs (MAGd1/2/3, MAGEA_dn0–dn3). Of these, **0/494
  lie in a K ≥ 2 family**: per-read edit distances are *identical* against both copies (NM_A = NM_B,
  1–3 mismatches = HiFi error floor). **0/311 reads differ between copies on MAGEA_p3.**
- The copies are **sequence-identical over the transcribed exonic region** — the 67–72% genomic divergence is
  entirely intronic/intergenic. No PSV column can adjudicate; there is no conflict edge in $H$; the minimum
  copy cover collapses the two copies into one.
- **Resolvable fraction: 0/494** (the exact prediction of the K = 0 branch of the K-bound corollary, §5).

This is the theorem made concrete: where K = 0 (copies identical over the observed columns), MCC = 1 regardless
of the true copy count, and attribution is provably arbitrary — exactly the forced-merge prediction.

### 7.3 Shared condition with detection (interest #1): the family as conflict-graph component

The upstream family-detection step (interest #1, `bench/family_definition_formal.md`) defines a multi-copy
gene family as a **connected component of the read-conflict graph** under the divergence-tie (`de`) edge
criterion. This is the same conflict-graph object that underlies the copy-assignment theory:

- The conflict graph $G = (V, E)$ used in detection has expressed loci as vertices and de-tied cross-mapping
  reads as edges; the family is the connected component.
- The conflict graph $H = (R, E)$ used in this note has the reads themselves as vertices and pairwise
  disagreement as edges; copy covers are colorings of $H$.

The two graphs are related: $G$ is the "locus-level summary" (which loci are confused) and $H$ is the
"read-level instance" (how to partition the confused reads). $G$'s components determine which reads populate
each $H$; the assignment problem then runs inside each component separately. Property P2 in
`bench/family_definition_formal.md` verifies this separation directly: **0 reads connect two distinct
multi-locus $G$-components**, so the assignment problems are independent across families — exactly the
prerequisite for running Theorem 2 per-family.

On the GGO panel: **TP = 7, TN = 10, FP = 0, FN = 0** (precision = recall = 1.000) in the ≤6-copy regime,
using the `de`-tie criterion with $\Delta = 0.005$. The same K quantity governs both detection (whether a
family is found) and assignment (whether the assignment is unique): copies with no PSV-distinguishing reads
(K = 0) form a family but yield a non-unique assignment; copies with K ≥ 2 form a family and yield a unique
assignment (under Strong Separation). Detection and assignment share the same identifiability boundary.

---

## §8 Discussion

### 8.1 The dichotomy as the through-line across interests #1, #2, and #3

The hardness/recovery dichotomy of Theorems 1 and 2 provides a single organizing principle across all three
research interests:

**Interest #1 (family detection).** The conflict graph $G$ on loci defines the family as the unit of
confusion. The K quantity that governs detection (whether two loci are in the same component, i.e., do they
share de-tied reads?) is exactly the same K that governs assignment identifiability. A family with K = 0
(copies identical over transcribed exons) is detectable via read-ambiguity — the reads are confused — but
unassignable: the minimum cover is non-unique, no algorithm can do better than arbitrary attribution. A family
with K ≥ 2 is detectable and, under adequate coverage, uniquely assignable. Detection and assignment are not
two separate problems; they share a boundary.

**Interest #2 (copy assignment, this note).** The dichotomy is the central result: the general problem is
NP-hard (Theorem 1), but the biologically relevant regime — copies with enough distinguishing variation and
reads covering those distinctions — falls under Strong Separation, where the truth is the unique minimum cover
(Theorem 2). The K-frontier Proposition sharpens this: the pairwise ($K = 2$) case is captured by the simpler
`sep+link` condition, but multi-copy ($K \geq 3$) requires the strictly stronger Strong Separation because
cross-copy recombination can assemble spurious alternative covers from fragments of three or more true copies.

**Interest #3 (allele-specific junctions).** The junction-as-column lift (§6) shows that the same formalism
extends to recover isoform structure jointly with copy identity. A long read spanning both a PSV and a splice
junction is a molecule that simultaneously assigns its copy (via the PSV column) and its isoform (via the
junction column). When Strong Separation holds over the combined (allele + junction) column set, the minimum
constrained path-cover of the family variation graph recovers copies **and** isoforms jointly — the exact
substrate for allele-specific junction analysis. Conversely, when K = 0 (copies identical over transcribed
exons), even a molecule spanning many junctions carries no PSV linkage; junction usage cannot be attributed to
a specific copy, and allele-specific junction analysis reduces to average-over-copies — a strictly weaker
signal. The K boundary governs the utility of interest #3's analysis.

### 8.2 The K = 0 frontier: when RNA alone cannot resolve

The K = 0 regime — copies sequence-identical over all transcribed exonic columns — is the frontier where RNA
alone provably cannot resolve copy attribution. This is not an algorithmic limitation; it is a fundamental
information-theoretic constraint: there is simply no signal in the reads that distinguishes the copies.

In this regime, the conflict graph $H$ has **no edges** between cross-copy reads (they are pairwise
compatible), so $\chi(H) = 1$: the minimum cover assigns all reads to a single copy, collapsing the true
copies into an undistinguishable group. No coloring/assignment algorithm can recover $K > 1$ copies from a
graph with chromatic number 1.

The levers available at the K = 0 frontier are necessarily *non-RNA*:

- **Longer molecules** (ultra-long ONT reads) that extend into flanking intronic or intergenic regions where
  the copies *do* diverge. If a molecule spans from a transcribed exon into a flanking intron that is
  divergent between copies, it observes a new PSV column and breaks the K = 0 barrier for that read.
- **Allele-specific junctions** (interest #3): if two copies share identical exonic PSV profiles but differ in
  their splicing patterns at specific junctions, a molecule spanning such a junction observes a junction column
  ($c_J$) that distinguishes the copies. This is the *only* RNA-accessible lever at K = 0 over the PSV
  columns — junction columns can restore K ≥ 1 over the combined column set even when the PSV columns alone
  give K = 0. Interest #3 is therefore not merely a downstream analysis step; it is an information-theoretic
  escape from the K = 0 regime for copies whose distinguishing variation is exclusively splicing-level.
- **Genomic phasing** (long-range haplotype) that carries extrinsic copy identity not observable in the
  transcript.

The MAGEA co-located arrays (§7.2, 0/494 resolvable reads) are the concrete instance of this wall. They are
not a failure of the algorithm — they are a demonstration that the wall exists and that RNA-level copy
resolution is genuinely impossible for this family on this substrate.

### 8.3 The K ≥ 3 cross-copy recombination obstruction: multi-copy is strictly harder than pairwise

The K-frontier Proposition (§5) is a genuine structural finding: **multi-copy ($K \geq 3$) recovery is
strictly harder than pairwise ($K = 2$)** in the following precise sense. There exists a natural and sufficient
condition for pairwise recovery (`sep+link`) that provably fails for three or more copies. The failure mode is
cross-copy *recombination*: a jointly-consistent class can be assembled from fragments of two different foreign
copies, realizing a novel allele-vector that is not any true copy, yet achieving a minimum-size cover.

The explicit witness (§5 Proposition (ii)) is:
$$
c_0 = (1,1,0),\quad c_1 = (0,0,1),\quad c_2 = (0,1,1)
$$
with each copy contributing reads on windows $\{0,1\}$ and $\{1,2\}$. The class $\{c_0\text{'s read on }\{1,2\},
\ c_2\text{'s read on }\{0,1\}\}$ is jointly consistent (they share only column 1, where $c_0 = c_2 = 1$) and
realizes the novel vector $(0,1,0)$. The cover $\bigl\{\{c_0\text{'s }\{1,2\}, c_2\text{'s }\{0,1\}\},
\ldots\bigr\}$ is a valid alternative minimum cover of size 3, so the assignment is non-unique despite
`sep+link` holding.

This finding has a practical consequence: haplotype-assembly algorithms designed for the diploid ($K = 2$)
regime — which standardly invoke the `sep+link` / error-free-phasing identifiability result — cannot be
directly lifted to the $K \geq 3$ multi-copy-family setting. The condition that suffices for all $K$ (Strong
Separation) is strictly stronger than `sep+link`, and the gap is not merely a proof artefact: the exhaustive
enumeration (§5) finds that `sep+link` has 248 uniqueness violations at $K = 3$ over the full $L = 3$
instance universe, while Strong Separation has zero violations at both $K = 2$ and $K = 3$.

The exact **necessary-and-sufficient** condition for unique recovery is **recombination-freeness** — no
alternative size-$K$ cover exists — which is instance-global and has no clean closed form (Proposition). Strong
Separation is the clean sufficient surrogate; whether there is a more permissive closed-form condition between
`sep+link` and Strong Separation that still works for all $K$ is an open question.

### 8.4 Polynomial recovery under Strong Separation: Theorem 3

Theorem 2 identifies the optimum with the truth (under Strong Separation, the true cover is the unique minimum
cover) but makes no claim about how to find it efficiently. Theorem 1 forecloses polynomial-time algorithms in
the general case. **Theorem 3 closes this gap under Strong Separation**: the algorithm `RECOVER` — build the
compatibility graph $\bar H$ (the complement of the conflict graph) and take its connected components — returns
the true copy set $C^\*$ in $O(n^2 m)$ time.

The key structural fact, established in the Theorem 3 proof (§5), is that Strong Separation makes $H$ a
**complete $K$-partite graph**, whose complement $\bar H$ is therefore a **disjoint union of $K$ cliques** — one
clique per true copy read-class. The connected components of $\bar H$ are exactly those read-classes, so
`RECOVER` extracts the truth directly from graph structure in $O(n^2 m)$ time. The self-certifying variant (§5)
additionally verifies, with an $O(n^2)$ incremental clique-check on the already-built $\bar H$, that $\bar H$ is
a disjoint union of cliques — certifying that `RECOVER`'s partition is the **unique minimum copy cover of the
reads** (not merely that the input is Strong-Separated; the property is weaker than Strong Separation). When
this check fails the input is **not** a disjoint-clique-union (the $K \geq 3$ recombination regime), and
`RECOVER` reports *not identifiable* rather than returning a partition it cannot guarantee.

**Honest caveat.** Strong Separation is a conservative *sufficient* condition for the disjoint-clique-union
property (and hence for the certifier passing). The empirical resolver may succeed under weaker conditions — e.g.
recombination-free instances that are not Strongly Separated, or $K = 2$ instances that are merely `sep+link`.
Theorem 3 is about the *certified* regime: when the disjoint-clique-union check passes (an $O(n^2)$ incremental
test on the already-built $\bar H$), `RECOVER` is provably correct and its partition is the unique minimum cover;
Strong Separation additionally guarantees that cover equals the true copies. The design of a polynomial-time solver
for the broader **recombination-free** regime — instances where the minimum cover is unique but the certifier may
not pass — remains an open question.

The dichotomy is therefore closed for the *combinatorial* recovery problem: **NP-hard in the general case**
(Theorem 1); **unique optimum and polynomial recovery** under Strong Separation (Theorems 2 and 3). The shipped
pipeline does not run `RECOVER`; **Theorem 4 (§5b, the bridge)** is the production-side certifier — a sound
per-read identifiability certifier for all K≥1.

---

## References

- **Lippert, R., Schwartz, R., Lancia, G., & Istrail, S.** (2002). Algorithmic strategies for the
  single nucleotide polymorphism haplotype assembly problem. *Briefings in Bioinformatics*, 3(1),
  23–31. — Introduces the Minimum Error Correction (MEC) formulation of haplotype assembly.
- **Cilibrasi, R., van Iersel, L., Kelk, S., & Tromp, J.** (2005). On the complexity of several
  haplotyping problems. In *Algorithms in Bioinformatics (WABI 2005)*, LNCS 3692, 128–139. (Journal
  version: *Algorithmica* 49(1):13–36, 2007.) — Establishes NP-hardness/APX-hardness of MEC, including
  at $k = 2$.
