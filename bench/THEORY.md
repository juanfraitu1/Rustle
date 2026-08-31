# Theory (consolidated)

> Merged from 3 source docs (verbatim; git keeps the originals' history).

**Contents:** copy_assignment_theory · F4_SCOPE · IDENTIFIABILITY_LIMITS · em_consistency


---

## copy_assignment_theory

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

1. **Family detection** (interest #1): which expressed loci are copies of one another? A multi-copy gene
   family is defined as an $R$-refined connected component of the **transcript-homology graph $E_r$** (the RNA
   homology oracle; see `bench/DEFINITIONS_FORMAL.md`). The **read-conflict (de-tie) graph is *not* the
   family boundary** — it is demoted to the **within-family O2 copy-assignment *decomposition*** (the confusion
   partition *inside* a fixed family). Its **cardinality** $\chi(H)=\mathrm{MCC}$, by contrast, is **re-promoted to
   an O1 family property** — the copy *count*, including the unassignable copies (§3 Remark;
   `bench/DEFINITIONS_FORMAL.md`). Copy-assignment is the unit of work performed *inside* each
   family; it is logically downstream of interest #1.

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
unique-mapper agreement panel (§7). This correspondence between the theoretical sufficient condition and the observed
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

### Remark (counting vs. assigning): χ(H) is a copy *count*, well-defined even when reads are unassignable

Lemma 1 separates two questions that Theorem 2 later ties together, and the separation is load-bearing for the
family-level (O1) reading of the conflict graph:

- **Counting copies** = $\mathrm{MCC} = \chi(H)$ — a function of the **conflict structure alone**. It is defined
  for *every* instance, with no separation hypothesis, and in particular **even when no read can be assigned to a
  specific copy**. It answers "how many copies must exist to explain these reads?"
- **Assigning reads** to named copies requires the **Strong-Separation** hypothesis of Theorem 2 (§5); that is the
  *resolvable* / SUN-identifiable core, a strictly stronger demand than counting.

> **Proposition (count needs only conflict; assignment needs Strong Separation).** $\chi(H)$ is a functional of the
> conflict relation $\sim$ **alone** — invariant under any relabelling of per-read origins, evaluable with no
> Strong-Separation hypothesis, no SUN, and no coverage condition; in particular it is defined even when **no** read
> of the family is assignable. *Single-read* assignability of a read $r$ ($\lvert N(r)\rvert=1$, Theorem 4(ii)) is a
> **strictly stronger** demand: it requires a single-position Strong-Separation witness (a SUN, §5·SUN), which
> Tier-2 copies lack by definition. Counting therefore consumes strictly less information than assigning — the
> count survives where the assignment abstains. This is the load-bearing separation: **the conflict says how many
> copies there are even when they are not resolvable; the resolvable subset is the nice core; the counted-but-
> unassignable copies are real and matter.**

> **Lemma (χ(H) = number of distinct observed hap-vectors = `psv_graph`'s $K$ *on the chosen H-substrate*).**
> ($\chi(H)$ is **substrate-dependent** but a **single invariant on one $H$**: the read-level `psv_graph` graph
> gives the *assignment-realized* count $\sum=354$, while the depth-independent **copy-consensus** graph
> (`copyonly_K`) gives the O1-authoritative count $\sum=361$; these are two substrates of the *same* Lemma, not two
> invariants — the O1 copy number uses copy-consensus, see §7.4 and the reconciliation below.) In the full-span regime
> (each read spans the family's PSV columns — the long-RNA case), group the reads by realized PSV hap-vector. Two
> reads with the **same** hap-vector agree at every co-observed column ⇒ non-adjacent; two reads with **distinct**
> hap-vectors differ at $\ge1$ co-observed column ⇒ adjacent. Hence $H$ — taken on the **de-tied / significance-
> gated** graph, **not** the raw allele-disagreement graph, whose error edges inflate the colouring to a median
> $\approx3\times K$ (`bench/conflict_graph_structure.py`) — is **complete multipartite** on the hap-vector groups,
> so $\chi(H)=n_{\text{groups}}$: precisely the per-family $K$ that `psv_graph_genomewide.json` already emits (its
> `group_sizes` are the multipartite parts; singleton parts $=$ Tier-1 $\uplus$ Tier-2, non-singleton parts $=$
> Tier-3 collapse). The O1 copy number is thus a field the pipeline already computes. (Partial reads only weaken
> this to $\chi(H)\le n_{\text{groups}}$ via recombinant covers — the $K\ge3$ non-uniqueness of the §5 Proposition;
> RNA long-reads sit near the full-span equality case.)

Hence a family can be **counted** as $\chi(H)$ copies while only a **subset** is single-read **assignable**: the
distinct-hap-vector copies that carry no single-position Strong-Separation witness (Tier-2, §5·SUN) are counted by
$\chi(H)$ yet are not single-read taggable. Because merging can only *reduce* the number of colors, $\chi(H)$ is a
**lower bound** on the true copy number: identical copies that share a hap-vector (Tier-3, the K-frontier) collapse
to one color, so $\chi(H)$ under-counts them by exactly $\sum_{\text{collapsed groups}}(\text{size}-1)$; recovering
those needs read-**depth** / DNA (parCN). The full honest chain is

$$
\underbrace{n_{\text{resolvable}}}_{\text{Tier-1, assignable}}\ \le\ \underbrace{\chi(H)}_{\text{conflict count (incl. Tier-2 unassignable)}}\ \le\ \underbrace{\chi(H)+\text{collapsed\_excess}}_{\text{= }n_{\text{loci}}\text{ (reference-resolved) = true\_copy\_lower\_bound}}\ \le\ \text{true copy number},
$$

with the last gap = reference-**absent** copies (O4). (The two middle terms swap order in the reference-*collapsed*
regime — a single reference locus hiding several hap-vectors, the SDA/Vollger case — where $n_{\text{loci}}\le\chi(H)$;
the general statement is $\max(n_{\text{loci}},\chi(H))\le\text{true}$.) This is why $\chi(H)$ is promoted to an **O1
family property** (the copy number) in `bench/DEFINITIONS_FORMAL.md`, while read→copy assignment stays
O2. Genome-wide on GGO (154 co-located families / 412 copies, `bench/copy_number_catalog.py`), on two instantiations of
$H$ that agree on **143/154** families. On the **read-observed** graph (`psv_graph_genomewide.json`, the Lemma's
directly-plumbed $\chi$): $\sum\chi(H)=354=322$ **read-level** singleton parts $+\,32$ collapsed parts, so
$354+58=412$ and **58** identical (Tier-3) copies go uncounted across **30/154** strictly-under-counting families
(the read-level singleton count **322** is a read-sampling quantity, **not** the tier count — the SUN tiers are
defined on the copy-consensus graph, where Tier-1 $\uplus$ Tier-2 $=339\ne322$, the gap being the 11 divergent families).
On the depth-independent **copy-consensus** graph (`copyonly_K`, read-sampling-invariant, authoritative for the O1
count): $\sum\chi(H)=361=338$ Tier-1 $+\,1$ Tier-2 $+\,22$ collapsed-group reps, missing **51** Tier-3 copies. The
two differ by exactly $+7$ over the 11 divergent families (10 consensus-over-splits at $+1$ each, minus fam 22's
read-over-split of $-3$); both obey $\chi(H)\le\sum n_{\text{loci}}=412$, so the choice is a graph-substrate choice,
not a change to the invariant.

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

#### Definition (MEC — a THEORY object; ⚠ NOT implemented in shipping Rustle)

> ⚠ **Corrected 2026-08-10.** This heading read *"as implemented in `vg_family/phasing.rs`"*. That
> module had **zero call sites** and named a CLI flag (`--vg-phase`) that never existed, and it was
> **deleted** on 2026-08-10. It also solved a *different* problem from the one this section is about:
> it was **diploid** MEC over a **binary** allele matrix at **one** locus (with `h_B = ¬h_A`), whereas
> O2 is *k*-copy over 4-letter alleles across a family. The `k = 2` construction below therefore
> stands as **mathematics with no shipped implementation**; the shipped O2 decision rule is a per-read
> argmax + Poisson-binomial certificate, specified in `docs/copy_assignment_definition.md`.


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

### §5c The capstone: copy assignment as facility location, with an LP-rounding approximation (Theorems 5–7)

Theorems 1–4 give *hardness* (Thm 1), *exact recovery under Strong Separation* (Thm 2–3), and a *per-read
certificate* (Thm 4) — but no **guarantee for the hard general instance**. This is exactly the gap the advisor's
own programme fills: Canzar et al. (2016) cast multimapping resolution as **maximum facility location** and solve
it by **LP relaxation + rounding** with provable approximation ratios. We adopt that paradigm for copy assignment.

**The problem (MWCA — Max-Weight Copy Assignment).** Clients are reads $R$; facilities are candidate copies $C$;
$N(r)\subseteq C$ is the set of copies *consistent* with read $r$ (no contradicted PSV); $w(r,c)\ge0$ for
$c\in N(r)$ is the read's evidence weight for copy $c$ — the production gate's per-read log-likelihood, the
principled replacement for the uniform $1/k$ split. Open a set $S\subseteq C$ with $\lvert S\rvert\le K$ and serve
each read by its best open consistent copy; the value is
$$f(S)=\sum_{r\in R}\ \max_{c\in S\cap N(r)} w(r,c)\qquad(\text{0 if }S\cap N(r)=\varnothing),\qquad \max_{|S|\le K} f(S).$$
By Lemma 1 a feasible $S$ of size $K$ is exactly a proper $K$-colouring of the read-conflict graph $H$ (each copy
= an independent set), so MWCA is the **weighted soft-assignment relaxation of the minimum copy cover**.

**The objective fork (stated up front, because the two halves differ sharply).** Minimising the cover size
$K=\chi(H)$ is graph colouring, hence **inapproximable** within $n^{1-\epsilon}$ for general $H$ — we keep that
only as the honest hardness boundary (Thm 1) and a *structured-case* result (below). The **max-weight** objective
above is the tractable one, and it is Canzar's facility-location problem verbatim.

> **Lemma (submodularity).** $f$ is monotone and submodular.
> *Proof.* Monotone: adding a copy can only raise a read's best option. Submodular: for $S\subseteq T$ and
> $e\notin T$, a read's marginal gain $\max\{0,\,w(r,e)-\max_{c\in S\cap N(r)}w(r,c)\}$ from opening $e$ is
> non-increasing in the incumbent best, which only grows from $S$ to $T$; summing over reads preserves it. $\square$

> **Theorem 5 (hardness).** MWCA is NP-hard.
> *Proof.* Weighted maximum coverage (NP-hard) is the special case $w\equiv1$, $N(r)$ arbitrary; equivalently it
> inherits hardness from Thm 1 via Lemma 1. $\square$

> **Theorem 6 (LP-rounding approximation — the capstone).** The natural LP relaxation
> $$\max \sum_{r,c} w(r,c)\,x_{r,c}\quad\text{s.t.}\ \textstyle\sum_{c} x_{r,c}\le1,\ \ x_{r,c}\le y_c,\ \
> \sum_c y_c\le K,\ \ x_{r,c}=0\ (c\notin N(r)),\ \ 0\le x,y\le1$$
> is an **upper bound** on the integer optimum, and the greedy / LP-pipage algorithm returns $S$ with
> $f(S)\ge\bigl(1-\tfrac1e\bigr)\,\mathrm{OPT}$. With pairwise copy-conflicts added (open $S$ must be conflict-free)
> the same LP rounds to a **constant factor** in Canzar's facility-location regime.
> *Proof.* Cardinality-constrained monotone-submodular maximisation: greedy is $(1-1/e)$-approximate (Nemhauser–
> Wolsey–Fisher), tight; the LP is a relaxation of the integer program, so $\mathrm{LP}\ge\mathrm{OPT}$. The
> conflict variant is submodular maximisation over an independence system, where dependent LP-rounding gives the
> Canzar constant. $\square$

> **Theorem 7 (integrality bridge — unifies §2–§5).** On a **Strong-Separation** instance the MWCA LP is
> **integral**, its optimum equals the true-cover value $\sum_r \max_{c\in N(r)} w(r,c)$, greedy is **exact**, and
> the per-read certificate $\mathrm{min\_p}$ (Thm 4) is the **complementary-slackness witness** of the tight
> assignment.
> *Proof.* Strong Separation $\iff$ every read is consistent with **exactly one** copy ($\lvert N(r)\rvert=1$,
> Thm 2). Then each $x_{r,\cdot}$ has a single free coordinate forced to its bound, the assignment polytope is
> integral, and opening the $K$ realised copies attains $\sum_r w(r,\text{its copy})$ — which greedy reaches by
> opening those copies. Unique consistency is precisely $\delta(r)\ge1$, i.e. $\mathrm{min\_p}(r)=\varepsilon^{\delta}<1$
> (Thm 4), the dual slack certifying read $r$'s assignment is determined, not guessed. So the approximation
> *collapses to exact recovery*, recovering Thm 2–3 as the integral special case. $\square$

This makes the theory one object: **NP-hard in general (Thm 1/5) → $(1-1/e)$ LP-rounding (Thm 6) → exact and
integral under Strong Separation (Thm 7) → per-read $\mathrm{min\_p}$ certificate (Thm 4)** as the dual witness.

**The flow-decomposition reading (the thesis framing).** With copies = paths through one PSV-aware variation
graph and reads = evidence on edges, MWCA is a **max-weight flow decomposition**: open $\le K$ paths carrying the
read-evidence flow. Path-flow LPs carry near-total-unimodular structure, which is the natural route to Theorem 7's
integrality and to *better-than-constant* ratios on the series-parallel / interval instances RNA actually
produces — and it lands the "flip Canzar: multimapping reads are **shared** evidence" message
([[project_thesis_framing_family_vg]]).

**Structured instances, and what the real graphs actually look like.** Where the read-conflict graph $H$ is
**chordal** (or perfect), $\chi(H)=\omega(H)$ is computable in polynomial time, so the *minimum cover itself* is
exact — strictly stronger than the worst-case bound. Empirically (`bench/conflict_graph_structure.py`, GGO
co-located families) this holds for **about one quarter** (7/30) of families. The dominant empirical fact, however,
is that the **raw** read-conflict graph is heavily **error-inflated**: its colouring exceeds the true copy
count $K$ by a **median $\approx3\times$** (far more in the worst cases — a $K=2$ family whose raw graph needs
tens of colours), because raw allele-disagreement counts
*sequencing-error* edges as conflicts. This is the empirical case that Lemma 1's $\mathrm{MCC}=\chi(H)$ is a
statement about the **error-free / de-tied** conflict graph — and hence that one should **not** solve $\chi$ on
the raw graph but use the per-read significance gate (Theorem 4) and MWCA's evidence-weighted relaxation
(Theorem 6), which absorb the error model instead of paying for it in spurious chromatic number.

**Verification.** `check_mwca` (`bench/mwca_approximation_check.py`, integrated into the suite) exhaustively
enumerates the small PSV universe (534 instance/$K$ pairs): the Lemma (diminishing returns on every $S\subseteq T$),
Theorem 6 ($f(\text{greedy})\ge(1-\tfrac1e)\,\mathrm{OPT}$ and $\mathrm{LP}\ge\mathrm{OPT}$ on every instance, plus a
genuine greedy-gap witness $5<6$), and Theorem 7 over the full Strong-Separation universe (166 instances: LP
integral, $=$ true-cover value, greedy exact, $\mathrm{min\_p}$ certificate tight).

---

---

## §5·SUN — the single-position private-allele witness (Sudmant 2010)

Strong Separation (§5) and the `min_p` bridge (§5b, Theorem 4) both turn on the per-read predicate
$\lvert N(r)\rvert=1$ (a read consistent with **exactly one** copy; Theorem 7). Sudmant et al. (2010, *Science*
330:641, doi:10.1126/science.1197005) named the concrete, single-position object that forces this predicate from
**one** read: the **Singly Unique Nucleotide (SUN)** — a position where one paralog's base is *unique* among all
copies of the family (a private/singleton allele). In their read-depth pipeline (→ QuicK-mer2 / fastCN,
Eichler lab) depth over SUNs yields **paralog-specific copy number** (parCN); they identified 4.1M SUN positions
genome-wide. We import the SUN as the concrete, checkable, single-read witness for **per-copy** identifiability,
sitting one level below the family-level Strong Separation / K-frontier objects.

**Definition (SUN / column-private allele).** Fix the family's copies $C^\*=\{c_1,\dots,c_K\}$ as allele-vectors
over the PSV columns $[m]$ (§2). Column $p\in[m]$ is a **SUN for copy $c$** if $c$'s allele is unique among all
copies at $p$:
$$ \mathrm{SUN}(c,p)\iff \forall c'\neq c:\ (c)_p\neq(c')_p \iff p\in\bigcap_{c'\neq c} D_{c,c'}. $$
i.e. at $p$ the base $(c)_p$ is a **singleton** among the $K$ copy-alleles (a PSV column, §2, has $\ge2$ distinct
alleles; a SUN column additionally has one copy *alone* against the rest). A copy $c$ is **SUN-identifiable** if
it has $\ge1$ SUN column. This depends **only on the copy sequences** (asm20 self-alignment) — no reads, no SEDEF.

*Scope of "unique".* Uniqueness is **within-family and conditional on catalog completeness**: unique among the
family's *enumerated* copies $C^\*$, not genome-wide. Sudmant's SUNs are genome-wide-unique; a called SUN here is
genome-wide-unique only if $C^\*$ is complete — a missing paralog carrying $(c)_p$ would demote it. The columns
come from an **all-to-longest-copy (star) projection**, not a true multiple-sequence alignment, and record only
aligned substitution columns (`matches_only`); this is **conservative** — it can only *under*-count SUNs (it drops
gap/ambiguous columns and copy-private indels, themselves excellent paralog markers).

**Lemma (SUN $\Rightarrow$ single-read determinism $=$ the $\lvert N(r)\rvert=1$ witness).** If $\mathrm{SUN}(c,p)$
then every read $r$ with $p\in\mathrm{obs}(r)$ and $r(p)=(c)_p$ is consistent with $c$ and **inconsistent with every
other copy**: $N(r)=\{c\}$, $\lvert N(r)\rvert=1$.
*Proof.* For $c'\neq c$, $(c')_p\neq(c)_p=r(p)$, so $r$ contradicts $c'$ at $p$; and $r$ agrees with $c$ everywhere
a $c$-read observes. $\square$
A SUN therefore discharges the Theorem 7 hypothesis and (Theorem 4(ii)) the **per-read gate** assigns $r$ to $c$
with $\mathrm{min\_p}=\varepsilon^{\delta(r)}<1$ **soundly from one position, one read** — no co-observation, no
phasing, no linkage — and $r$ can **never be misassigned to another true copy**. This is a **per-read** guarantee
(the gate, Thm 4); it is the *coverage-free* strengthening of Strong Separation's per-read demand. It is **not**
a cover-level guarantee about the *copy* (see "The recombination link" below): the copy can still be dissolved by
an alternative minimum cover in the NP-hard `RECOVER`/MCC problem — which is precisely why production runs the
per-read gate, not `RECOVER`.

**Honest relation to Strong Separation (NOT iff).** SUN is neither sufficient nor necessary for *family-level*
Strong Separation:
- **SUN $\not\Rightarrow$ family Strong Sep.** Strong Separation is an all-pairs, coverage-dependent property
  (every cross-copy read pair conflicts). A single copy's private column constrains only $c$-vs-others at $p$;
  another pair $(c',c'')$ may share every column a read observes and stay unseparated. (Witness
  $c_0{=}(1,0),c_1{=}(0,0),c_2{=}(0,1)$: $c_0$ has a SUN at col 0, but the $(c_1,c_2)$ pair is under-separated for
  reads that miss col 1.)
- **Strong Sep $\not\Rightarrow$ SUN.** A Strong-Separated copy need not carry a private column — its identity can
  be a **combination** of individually-shared alleles. Machine-checked witness $c_0{=}(0,0),c_1{=}(0,1),c_2{=}(1,0)$:
  all three hap-unique and Strong-Separated under full-length reads, yet **$c_0$** has **no** SUN (col 0 shared with
  $c_1$, col 1 with $c_2$) — so *some* copy is Strong-Separated without a SUN (SUN not necessary; here $c_1,c_2$ do
  have SUNs). A genuine **all-no-SUN** Strong-Separated family needs $K\ge4$: $\{(0,0),(0,1),(1,0),(1,1)\}$ — every
  column's two symbols each appear twice, so **no** copy has a SUN, yet all four vectors are distinct (S4).

  So SUN is a **one-directional, single-position, single-read sufficient witness** for the per-read
  identifiability predicate $\lvert N(r)\rvert=1$ that Strong Separation (Thm 2/7) and the gate (Thm 4) rely on —
  *locally stronger* (one read / one column) than Strong Separation, *globally weaker* (silent about other pairs and
  coverage). $\mathrm{SUN}(c)\Rightarrow$ [$c$ is single-position-separated from all others] $\wedge$ [reads over the
  SUN are deterministically $c$]; the converse fails in both directions.

**The 3-tier per-copy identifiability ladder (refines psv_graph's family $K$).** `psv_graph_genomewide.py`
reports, per family, $K=$ the number of distinct copy hap-vectors (`n_groups`); a copy is "resolvable" iff its full
hap-row is a **singleton** group. The SUN lens splits each copy:
1. **SUN-identifiable (Tier 1):** $\exists$ SUN column. Single-read **gate-deterministic** — a read over the SUN
   column certifies $\lvert N(r)\rvert=1$ (Thm 4/7) and is **never misassigned to another true copy**. (This is
   per-read immunity, **not** cover-level immunity of the copy — see "The recombination link" below.)
2. **hap-vector-unique-only (Tier 2):** full hap-row is a singleton group but the copy has **no** SUN. Uniqueness is
   *combination-based* $\Rightarrow$ **no single read pins it**; resolving it needs a read/class **co-observing
   $\ge2$ PSVs** (linkage/phasing). At the per-read gate a single-column read on this copy is always ambiguous
   ($\lvert N(r)\rvert\ge2$), and its multi-PSV identity is exactly the fragment a **recombinant** read can spoof.
3. **K-frontier / unresolvable (Tier 3):** full hap-row equals another copy's (non-singleton group) $\Rightarrow$
   PSV-identical, `NM:i:0` collapse; Strong Sep fails on that pair; Thm 4(iv) certifies $\mathrm{min\_p}=1$ (tied).

The ladder is a **strict refinement** of $K$: {resolvable copies} $=$ Tier 1 $\uplus$ Tier 2 (the singleton groups
`psv_graph` already counts), {frontier} $=$ Tier 3. **SUN-identifiable $\subsetneq$ hap-vector-unique** whenever
Tier 2 is nonempty (unique-but-no-private-column copies exist). $\mathrm{SUN}(c)\Rightarrow$ $c$'s hap-vector is
unique (Tier 1 $\subseteq$ singleton), machine-checked with **0** violations.

**The recombination link — per-read immunity, NOT cover-level immunity (Proposition).** A SUN gives $c$
**per-read gate immunity**: any read over the SUN column $p$ is pinned to $c$ ($\lvert N(r)\rvert=1$, Lemma) and
the gate (Thm 4) can never route it to another *true* copy. Equivalently, at the **read** level,
$$\{\text{reads ambiguous among }\ge2\text{ true copies}\}\ \subseteq\ \{\text{reads carrying no copy's SUN allele}\}$$
(if $r$ carries any copy's private allele it is single-read pinned; contrapositive of the Lemma, machine-checked
S2). This is a statement about **reads** and the **per-read gate**.

It does **not** upgrade to **cover-level** immunity of the *copy* in the NP-hard `RECOVER`/MCC problem (§5,
Theorem 2). The naive claim — "no alternative cover can realize $c$'s private allele, so a SUN copy is unspoofable"
— is **false**: a phantom class can carry $(c)_p$ by absorbing $c$'s **own** SUN-covering reads. On the theory's
canonical $K=3$ witness $c_0{=}(1,1,0),\,c_1{=}(0,0,1),\,c_2{=}(0,1,1)$ with its exact 6-read set (each copy on
windows $\{0,1\},\{1,2\}$), $\mathrm{MCC}=3$ admits a **second** minimum cover
$$\{\underbrace{(0,0,1)}_{=c_1}\},\ \{\underbrace{(0,1,0)}_{\text{phantom}}\},\ \{\underbrace{(1,1,1)}_{\text{phantom}}\}$$
that **dissolves $c_0$** — even though $c_0$ has SUNs at **both** cols 0 and 2. The phantom $(1,1,1)$ carries
$c_0$'s private $0{\to}1$ at col 0 and $(0,1,0)$ carries $c_0$'s private col-2 allele, precisely because $c_0$'s own
reads ($c_0$ on $\{0,1\}$ and on $\{1,2\}$) are absorbed into them. So there is **no** cover-level containment
$\{\text{spoofable}\}\subseteq\{\text{no-SUN}\}$, and no "hijacked copy is exactly $c_2$" reading: across the two
minimum covers, $c_0$ (SUN-rich) and $c_2$ (no SUN) **both** dissolve; $c_1$ happens to survive here but that is
incidental, not proven general.

**The right conclusion (and why it matters for the pipeline).** SUN determinism lives at the **per-read gate**
(Thm 4), which the shipped pipeline runs — not in the NP-hard `RECOVER`/MCC cover, which no SUN protects. This is a
positive argument *for* the gate over cover-recovery: the gate inherits SUN's single-read soundness ($\lvert
N(r)\rvert=1$), whereas the minimum-cover objective is genuinely ambiguous under the $K\ge3$ recombination
obstruction regardless of how many private columns a copy has. Tier 1 vs Tier 2 is therefore a **per-read**
distinction (single read pins it vs needs $\ge2$-PSV co-observation), not a cover-survival guarantee.

**Verification.** `bench/sun_theory_check.py` (exhaustive over distinct copy sets $K\in\{2,3,4\}$, $m\in\{1,2,3\}$,
$\lvert A\rvert=3$): **(S1)** SUN $\Rightarrow$ unique hap-vector — **0 violations / 1,252,380 SUN-copies**;
**(S2)** every read over a SUN column carrying the private allele has $N(r)=\{c\}$ — **0 violations / 6,675,294
reads** (per-read gate immunity); **(S3)** the canonical $K\ge3$ witness — factual report only: recombinant class
realizes the novel vector $(0,1,0)$, unique no-SUN copy $=c_2$ (reported as a *value*, not read as an immunity
claim); **(S3\_cover)** the **load-bearing counterexample** — enumerates all minimum covers of the 6-read set,
confirms $\mathrm{MCC}=3$, exhibits the alternative cover that **dissolves the SUN-rich $c_0$** with phantoms
carrying $c_0$'s private alleles, and simultaneously confirms per-read gate immunity still holds
(`cover_level_copy_immunity_FALSE=true`, `per_read_gate_immunity_holds=true`); **(S4)** NOT-iff both directions —
a Strong-Separated hap-unique copy with **no** SUN (SUN not necessary), plus a $K=4$ instance
$\{(0,0),(0,1),(1,0),(1,1)\}$ where **no** copy has a SUN yet the family is fully Strong-Separated. All green.

The real-family catalog `bench/sun_catalog_fast.py` (copy-only asm20 aligned-pairs, no reads/SEDEF) tiers every copy
of the **154** GGO validated multi-copy families (**412 copies**): **Tier 1 SUN-identifiable = 338 (82.0%)**,
**Tier 2 hap-vector-unique-only = 1**, **Tier 3 frontier/collapsed = 73 (17.7%)**. It re-checks
$\mathrm{SUN}\Rightarrow$ unique-hap with **0 violations** on real data, and the strict subset
**SUN-identifiable $\subsetneq$ hap-vector-unique** ($338<339$) is witnessed by exactly **1** Tier-2 copy. The
**empirical** message is sharp and honest: on this substrate the Tier-1/Tier-2 gap is nearly vacuous — essentially
every resolvable gorilla copy earns its uniqueness through at least one *single-position private allele*, so **82%
of copies are single-read gate-deterministic** (a read over the private column is never misassigned to another true
copy) and the $K\ge3$ Tier-2 regime — where no single read pins the copy and $\ge2$-PSV co-observation is needed —
is empirically rare (1/339). (This copy-level catalog drops `psv_graph`'s read-support gate, so it is a
*superset* of the observed catalog: 135 families carry a SUN-identifiable copy at the sequence level vs the 126
read-supported in `psv_graph_genomewide.json`; identifiability is a property of the copies, observability of the
reads.) See `bench/sun_identifiability_catalog.json`.

**Canonical deliverable + on-real-data Strong-Sep machine-check.** `bench/sun_identifiability.py` is the
consolidated per-copy catalog (copy coordinates + gene annotation), emitting `bench/sun_identifiability.tsv`
(per-copy: family, copy, chrom, start, end, gene, n_psv, n_sun, tier, hap_unique, group_size) and
`bench/sun_identifiability.json` (genome-wide tier breakdown, K-refinement, examples). Beyond the exhaustive
abstract check `sun_theory_check.py` (S1–S4, S3_cover), it runs a **self-consistency** machine-check on the real
families: for every copy it recomputes — from the raw allele dicts, *not* from the `has_sun` flag — the
single-position Strong-Separation witness set $W(c)=\{p:\forall c'\neq c,\ (c)_p\neq(c')_p\}=\bigcap_{c'\neq
c}D_{c,c'}$, and verifies (i) $W(c)\neq\varnothing\iff c$ is SUN-identifiable, (ii) every read over a $p\in W(c)$
carrying $(c)_p$ has $N(r)=\{c\}$, $\lvert N(r)\rvert=1$, (iii) SUN $\Rightarrow$ singleton hap-vector. Result on
the 412 copies: **338 SUN-identifiable $=$ 338 with a single-position Strong-Sep witness**
(`SUN_equals_witness=true`), **0** violations of any kind (`all_green=true`). **Honest scope:** the
witness recompute uses `all(hap[c]\neq hap[c'])` while `has_sun` uses `count(hap[c])==1` — these are the **same
predicate** (copy $c$ contributes exactly one occurrence to a column), so `SUN_equals_witness` is true **by
construction**: a no-coding-bug / internal-consistency check, **not** independent corroboration of SUN. Genuine
independence comes from the abstract checks (S1/S2/S3_cover) and from re-deriving the tier counts from raw
sequence. K-refinement: exactly **1** family (`family 42`, an 8-copy `LOC129529768`/`LOC129529*` block, 49 PSVs)
is `fully_resolvable` by `psv_graph`'s $K$ yet contains a Tier-2 copy (`copy4`, `LOC129529774`: unique hap-vector,
`n_sun=0`, `n_witness=0`) — the lone copy with no single-position private allele, whose resolvability needs
multi-PSV co-observation; $K$ over-counts the single-read gate-taggable copies by exactly that 1. Output is
byte-identical across re-runs (deterministic asm20).

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

*See also: `bench/DEFINITIONS_FORMAL.md` for the upstream family-detection step whose output populates
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

### 7.1 Theorem 2 / K ≥ 2: sim5x K-ladder and GGO unique-mapper agreement

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

**GGO unique-mapper agreement.** On the real GGO HiFi IsoSeq panel, the `detect_and_assign` MAGEA smoke run reports
**1026/1026 (100%)** reads correctly assigned in the K ≥ 2 families (RABL2, AK6, CCDC196). These families are
all in the well-separated regime (MAPQ > 0, divergence-gap decisive), where Strong Separation holds by the
aligner's own evidence; the coverage condition is met at ≥ 47 de-conflict reads per family.

> ⚠ **What this figure is and is NOT.** It is a *consistency check in the EASY (MAPQ > 0) regime* — and
> the "unique-mapper agreement" truth is the aligner's own primary placement, so the metric is **circular by construction**
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

### 7.3 Shared condition with detection (interest #1): the O2 conflict-graph decomposition inside a family

The upstream family-detection step (interest #1, `bench/DEFINITIONS_FORMAL.md`) defines a multi-copy
gene family as an $R$-refined **connected component of the transcript-homology graph $E_r$** (the O1 homology
oracle). The **read-conflict (de-tie) graph is *not* the O1 family boundary** — it is the **O2**
copy-assignment decomposition *inside* a fixed $E_r$ family (`DEFINITIONS_FORMAL.md`, with
$E_c^{\mathrm{sig}}\subseteq E_c\subseteq E_r^{\mathrm{asym}}$). That O2 conflict-graph object is what
underlies the copy-assignment theory here:

- The conflict graph $G = (V, E)$ (an **O2 object, one per $E_r$ family**) has that family's expressed loci as
  vertices and de-tied cross-mapping reads as edges; its components are the exact-decomposition units of
  assignment.
- The conflict graph $H = (R, E)$ used in this note has the reads themselves as vertices and pairwise
  disagreement as edges; copy covers are colorings of $H$.

The two graphs are related: $G$ is the "locus-level summary" (which loci are confused, **inside one $E_r$
family**) and $H$ is the "read-level instance" (how to partition the confused reads). $G$'s components
determine which reads populate each $H$; the assignment problem then runs inside each component separately.
Property P2 verifies this separation directly: **0 reads connect two distinct multi-locus $G$-components**, so
the O2 assignment sub-problems are independent across those components (and a fortiori across the
homology-disjoint $E_r$ families) — exactly the prerequisite for running Theorem 2 per O2-component.

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
- **Sudmant, P. H., Kitzman, J. O., Antonacci, F., Alkan, C., Malig, M., Tsalenko, A., Sampas, N.,
  Bruhn, L., Shendure, J., & Eichler, E. E.** (2010). Diversity of human copy number variation and
  multicopy genes. *Science*, 330(6004), 641–646. doi:[10.1126/science.1197005](https://doi.org/10.1126/science.1197005)
  (PMID 21030649). — Defines the **Singly Unique Nucleotide (SUN)**: a position where one paralog's base
  is unique among the copies of a highly-identical gene family. Identifies **4.1 million** SUN positions
  genome-wide and uses short-read depth over them to genotype **paralog-specific copy number** (parCN) —
  the marker later systematized by QuicK-mer2 / fastCN (Eichler lab). §5·SUN imports the SUN as the
  concrete single-position, single-read witness for the per-copy identifiability predicate $\lvert N(r)\rvert=1$.


---

## F4_SCOPE

# F4 — scope: the Canzar-shaped theory capstone (LP-rounding approximation for copy assignment)

**Goal.** Complete the identifiability theory's missing piece: an **approximation algorithm with a provable
ratio** for the NP-hard general case, cast in the advisor's own paradigm (Canzar 2016 = multimapping resolution
as **facility location**, solved by **LP relaxation + rounding** with approximation guarantees). Today we have
hardness (Thm 1), exact recovery under Strong Separation (Thm 2–3), and a per-read certificate (Thm 4) — but
**no guarantee for the hard general instance**, which is exactly the shape his 2016 paper fills.

This is a scope/plan, not a proof. It fixes the formulation, the target theorems, the algorithm, the machine-
check, and the honest boundaries, so the proof work is well-defined.

---

## 1. The formulation (faithful to Canzar's facility-location frame)

**MAX-WEIGHT COPY ASSIGNMENT (MWCA).**
- **Clients** = reads `R` (partial allele-functions over PSV columns).
- **Facilities** = candidate copies `C` (allele-vectors / haplotypes; from the family layer, or the maximal
  independent sets of the de-tie conflict graph `H̄`).
- **Compatibility** `N(r) ⊆ C` = the copies read `r` is consistent with (no contradicted PSV).
- **Weight** `w(r,c) = ` the read's log-likelihood for copy `c` (the production gate's per-column LLR), defined
  only for `c ∈ N(r)`. This is the principled replacement for `1/k` the advisor dislikes — a real evidence
  weight, not a uniform split.
- **Decision:** open a set `S ⊆ C` of copies and assign each read to one open compatible copy
  `σ(r) ∈ S ∩ N(r)` (or leave it unassigned), to **maximize `Σ_r w(r, σ(r))`** subject to a budget `|S| ≤ K`
  (the cover size, `= χ(H)` by Lemma 1) — or, in the facility-cost variant, minus `Σ_{c∈S} f_c`.

This is **uncapacitated facility location / max-coverage with a cardinality budget** — Canzar's exact paradigm.
Two objective variants, with very different approximability (the key honest fork, see §2).

**Bridge to the existing object.** `S` is a conflict-free cover (each copy = an independent set in the read-
conflict graph `H`; Lemma 1: a feasible `S` of size `K` ⟺ a proper `K`-coloring of `H`). So MWCA is the
**weighted, soft-assignment relaxation of the minimum cover** — the cover decides *which* copies, MWCA decides
the best evidence-weighted assignment given that choice.

---

## 2. The objective fork (decide first — it determines what's provable)

| Objective | What it is | Approximability | Verdict |
|---|---|---|---|
| **MIN-COVER** (minimize `K = χ(H)`) | fewest copies to explain all reads | graph coloring → **inapproximable** within `n^{1−ε}` for general `H` (unless P=NP) | **Do NOT target in general.** Only tractable on *structured* `H` (perfect/interval/bounded-degree). Keep as an honest hardness boundary + a structured-case result (see §5). |
| **MAX-ASSIGNMENT** (maximize weight, `≤K` copies) | best evidence-weighted assignment under a copy budget | submodular max-coverage → **(1 − 1/e)** via greedy/LP-pipage; with pairwise copy-conflicts → Canzar-style **constant-factor** LP-rounding (his 0.19–0.44 regime) | **TARGET THIS.** Clean constant-factor guarantee, Canzar's exact method, evidence-weighted (no 1/k). |

**Recommendation: target MAX-ASSIGNMENT.** It is the version that *has* a constant-factor guarantee, it is
literally Canzar's facility-location problem, and it uses the gate's LLR weights (his aesthetic). The MIN-COVER /
coloring side becomes the honest *inapproximability boundary* of the chapter (we already have its NP-hardness in
Thm 1), plus a tractable result on the structured conflict graphs our data actually produces.

---

## 3. Target theorems (numbered to extend the note)

- **Thm 5 (NP-hardness of MWCA).** MWCA is NP-hard (reduce from Thm 1 / max-coverage). *(Likely short.)*
- **Thm 6 (LP-rounding approximation — the capstone).** The natural LP relaxation of MWCA
  (`x_{r,c} ∈ [0,1]` assignment, `y_c ∈ [0,1]` opening; `Σ_c x_{r,c} ≤ 1`, `x_{r,c} ≤ y_c`, `Σ_c y_c ≤ K`,
  `x_{r,c}=0` for `c ∉ N(r)`) admits a rounding achieving a **constant factor** of the optimum:
  `(1 − 1/e)` for the cardinality-budget / no-conflict case (pipage / dependent rounding), degrading to a
  Canzar-style constant when pairwise copy-conflicts are added. *(The core proof effort.)*
- **Thm 7 (integrality bridge — ties F4 to Thm 2–4).** Under **Strong Separation**, the MWCA LP is **integral**
  and its optimum equals the true cover: the approximation collapses to **exact recovery**, recovering Thm 2/3
  as the integral special case, and the per-read certificate (Thm 4, `min_p`) is exactly the **LP dual /
  complementary-slackness witness** that certifies a read's assignment is tight. *(This is the elegant glue that
  makes the whole theory one object: hardness → approximation → exact-under-separation → per-read certificate.)*

---

## 4. The VG / flow-decomposition reading (the thesis framing, optional second lens)

Per `project_thesis_framing_family_vg`: copies = paths through one PSV-aware variation graph; reads = evidence on
edges. MWCA then reads as a **min-cost / max-weight flow decomposition**: open `≤K` paths (copies) carrying the
read-evidence flow, maximizing assigned weight. LP-rounding of path-flow = the same theorem in flow language, and
flow LPs are often integral (totally-unimodular structure) — a candidate *route to Thm 7* and to better-than-
constant ratios on the structured (interval/series-parallel) instances RNA produces. Worth exploring as the proof
vehicle for Thm 6/7; it also lands the "FLIP Canzar: multimapping = shared evidence" thesis message.

---

## 5. Structured-instance angle (where MIN-COVER becomes tractable — bonus)

Real read-conflict graphs are not arbitrary: PSV-column conflicts give them interval/comparability-like
structure. If `H` (or `H̄`) is shown to be **perfect** (or interval, or bounded-treewidth) on co-located
families, then `χ(H)` = MIN-COVER is **poly-time exact** there — a much stronger statement than the constant-
factor MAX-ASSIGNMENT bound, for the instances that matter. **Action:** characterize the conflict-graph class
empirically (are the GGO family conflict graphs perfect? chordal? interval?) and prove the structural lemma if so.
This is the highest-payoff side quest and may be easier than the general approximation.

---

## 6. Machine-check plan (mirror the existing exhaustive style)

`bench/mwca_approximation_check.py`, exhaustively over small universes (K≤3, m≤3 PSV columns, |R| small):
- enumerate instances; solve the **LP** (scipy.optimize.linprog or PuLP) and the **integer optimum** (brute force).
- verify the **rounded** solution achieves the claimed ratio vs the integer OPT on *every* instance (the bound is
  not just expected-case);
- verify **Thm 7**: on Strong-Separation instances the LP is integral and `= ` the true cover, and `min_p`
  matches the dual slack. (Reuse the Strong-Sep generator from `copy_assignment_theory_checks.py`.)
This makes Thm 6/7 machine-witnessed the way Thm 1–4 already are.

---

## 7. Honest boundaries (state up front, the way the note already does)

- MAX-ASSIGNMENT is approximable; **MIN-COVER (coloring) is not** in general — the chapter must lead with the
  objective choice, not blur them. (This is the single most important honesty point.)
- The constant factor (≈`1−1/e`, or Canzar's `~0.19` with conflicts) is a *worst-case* guarantee; real instances
  do far better (the structured-case §5 + integrality §3 explain why). Report both.
- `w(r,c)` inherits the gate's assumptions (error model, editing filter); the approximation is over the *given*
  weights, not a claim the weights are perfect.
- This is **theory**, not a production change — like Thm 1–4, `RECOVER`/MWCA need not run in the pipeline; the
  shipped gate is Thm 4. F4 guarantees the algorithm you *could* run, in his paradigm.

---

## 8. Work plan (phases, each independently checkable)

1. **Fix the formulation + objective** (this doc; recommend MAX-ASSIGNMENT). — *decision*
2. **Thm 5** NP-hardness of MWCA. — *short proof + note section*
3. **Thm 6** LP + rounding + ratio, MAX-ASSIGNMENT (start no-conflict `1−1/e`; then add conflicts → Canzar
   constant). — *core effort; try the flow-decomposition vehicle §4*
4. **Thm 7** integrality under Strong Separation + `min_p` = LP-dual witness. — *the glue to Thm 2–4*
5. **Structured-instance lemma** (§5): characterize + (if perfect) exact `χ(H)` on co-located families. — *bonus*
6. **`bench/mwca_approximation_check.py`** exhaustive machine-check (§6). — *witnesses 5–4 like the rest*
7. **Write §5c "MWCA / facility-location" into `bench/copy_assignment_theory.md`** + update the checks suite.

**Smallest valuable slice (if time-boxed):** steps 2 + 3-no-conflict + 6 give a clean, machine-checked
`(1−1/e)` LP-rounding theorem in Canzar's paradigm — already the missing capstone. Steps 4/5 are the elegance
multipliers.


---

## IDENTIFIABILITY_LIMITS

# The identifiability boundaries — what is principled-limit, not loose-end

A central thesis claim is that copy-resolution from RNA has **identifiability boundaries** that the theory
*predicts* and the methods *certify*, rather than gaps to be closed by more engineering. Two objectives sit
at those boundaries by their nature. Framing them as boundaries (and reaching-but-not-crossing them, with
certified abstention) is the result — crossing them requires orthogonal data (DNA), explicitly outside the
RNA-only scope.

## K=0 boundary — exonically-identical co-located copies

**Statement.** Two co-located copies whose spliced (exonic) sequences are byte-identical share *no*
distinguishing feature in RNA: a read carries the same bases/junctions whichever copy it came from. By
construction there is no PSV and no copy-specific junction → the read is **unidentifiable**.

**Theorem-predicted, not a failure.** This is the K=0 vertex of the identifiability spectrum (copy_assignment
theory: a copy is resolvable iff it carries a distinguishing feature; the Strong-Separation / K-frontier
results characterize exactly when a unique cover exists). K=0 violates the precondition.

**Empirically confirmed, and the substrate is pinned (re-derived 2026-08-11).** **GORILLA** — mGorGor1
(`GGO.fasta`), reads `GGO_mm.bam` (OR6737 IsoSeq; note `GGO.bam` is a *symlink* to it), **no annotation**:
the pair is de-novo family `DNFAM1`, coordinates in `bench/copy_resolution_census.tsv`. At the two
exonically-identical co-located loci the locus references are byte-identical — re-aligned today with
`minimap2 -x asm20 -c --cs`, `NC_073247.2:164381193-164384845` vs `164442446-164446047` gives
**nmatch 3599 / blocklen 3599**, and `164397062-164401094` vs `164424321-164430225` gives **4030 / 4030**,
identity 1.000000 both — and **0/386 reads** are resolvable. ⚠ State the denominator's definition when
quoting it: 386 = reads aligned at **both** loci carrying a MAPQ-0 alignment **and holding a primary
record at one of the two** (311 + 75); without the primary restriction the same query returns 1,692.
The **numerator is entailed, not measured**: identical references ⟹ the distinguishing-column set is
empty ⟹ 0 resolvable at any denominator. A hard floor, not a tuning issue.

**Now CERTIFIED by the significance gate.** The IsoCon gate routes such reads to `Tied` because every
competitor has an empty distinguishing set → `min_p_value = 1.0 ≥ α`. So the method does not guess (no 1/k);
it *proves* unresolvability per read. Validated on the sim5x **K=0 rung → 100% Tied**.

**Escapes (require leaving RNA-only or aggregating).** (i) reference-vs-reference NM gate (use the assembly,
not the reads); (ii) aggregate per-family quantification (the count is identifiable even when the partition
is not — copy_assignment theory §6b: famCN is a sufficient statistic, parCN is not); (iii) DNA.

## O4 boundary — copy vs allele (reference-absent confirmation)

**Statement.** Distinguishing a genuine extra **copy** from an extra **allele** (heterozygous site) is
information-theoretically impossible from RNA alone: both add one haplotype's worth of sequence to the pool.
RNA carries expression, not copy number.

**So O4 is DETECT-and-FLAG, bounded by design.** We *can* flag candidate reference-absent copies from
RNA-visible signals — collapsed-CNV (≥12 balanced co-segregating alt columns → 18 strong candidates) and
divergent-unmapped — and quantify the flag precision (**FP bound 7.39%**, `project_objectives_status`). What
RNA cannot do is **confirm** copy-vs-allele; that needs **DNA copy-number (parCN, Soto-2025-style)**. The
divergent-unmapped route was additionally **dry on T2T** (the reference already contains them), a second,
data-specific boundary.

**Framing.** O4 attained = *the RNA-attainable half* (detect + flag + FP-bound). The unattainable half
(confirmation) is an information boundary, and the catalog flags candidates *for* a DNA follow-up rather than
overclaiming them — the honest, theorem-consistent position.

## Why this is a contribution, not a shortfall

The methods **assign-or-abstain up to the boundary** and **certify** the boundary per read (Tied via
`min_p ≥ α`; flag-not-confirm for O4). That is exactly the Canzar-aesthetic result: a clean combinatorial
identifiability criterion, provable boundaries, no arbitrary 1/k apportionment across the unidentifiable.
Reaching the boundary with a certificate is the deliverable; the orthogonal-data escapes are named, not
pretended.


---

## em_consistency

# EM copy-assignment = soft SDA PSV-clustering: the MLE consistency theorem under the heuristic

*The EM copy-assignment engine (`src/rustle/vg_family/em_copy_assign.rs`) is the maximum-likelihood **soft
relaxation of SDA's PSV correlation-clustering** (Vollger et al., Nat Methods 2019 — the advisor's prior-art
paper), run on the thesis's PSV-aware variation graph, and it is **consistent** in the identifiable regime.
Every piece is forced: the per-read likelihood is SDA's attraction/repulsion made continuous, the E-step is
SDA's read partition made soft, and the identifiable regime is exactly the Strong-Separation ⟹ unique
minimum copy cover (MCC = χ(H)) result of the copy_assignment_theory section above. The consistency theorem
is the provable layer **under** SDA's NP-complete heuristic — and it explains SDA's 91–93% accuracy floor.
Design spec: `docs/superpowers/ARCHIVE_INDEX.md` → `2026-07-08-em-copy-assignment-design.md` (in git history).*

## Derivation — the EM is SDA's PSV graph made continuous

**SDA's PSV graph (Vollger 2019, `reference_sda_vollger`).** SDA (Segmental Duplication Assembler) resolves
**collapsed** segmental duplications from long DNA reads. After detecting a collapse by **read-depth excess**,
it calls **PSVs** (second-most-frequent base at positions whose total depth ≈ a multiple of single-copy
coverage — the coverage gate separates PSVs from heterozygous alleles) and builds a signed **PSV graph**:
nodes = PSVs; a read carrying two PSVs **on one molecule** = an **attraction** edge (same copy); two PSVs
**mutually exclusive** across reads = a **repulsion** edge (different copies). SDA runs **correlation
clustering** on this signed graph — *ab initio*, **no preset cluster number**, **NP-complete**, solved
heuristically with **15 random initializations** — then **WhatsHap** partitions reads by the PSV clusters.
SDA validates at **91–93%** of PSVs correctly assigned (SRGAP2 / NOTCH2NL, BAC ground truth) and states the
limit outright: *"virtually identical duplications cannot be distinguished and will require even longer reads
(>100 kb)."*

**The same object on the PSV-aware VG** (`project_psv_aware_vg`, `psv_linkage.rs`): bubbles `[m]` = PSV
columns (alphabets `A_j`), parallel paths = copies, copy `k` = allele-vector `θ_k ∈ ∏_j A_j` (built from
read-supported PSV bubbles, SDA-style, not a handed-down catalog); a read `r` is a partial path `obs(r) ⊆ [m]`
with allele `r(j)`. Two PSVs co-carried by one read = SDA's attraction edge. This is the identical
column/allele model of the §2 above.

**Per-read likelihood = SDA's ±1 edges made continuous.** Emission at a bubble
`q_j(a\mid b) = 1−e_j` if `a=b`, else `e_j/(|A_j|−1) ≈ e_j/3` (DNA), with `e_j` = per-column error. Then
`L_{rk} = ∏_{j∈obs(r)} q_j(r(j)\mid(θ_k)_j)` (this is the existing `ReadEvidence.logl`). At a bubble carrying
copy `k`'s private allele (a SUN column, `reference_sudmant_2010_sun`), the match term `log(1−e_j)` is a
**soft attraction to `k`** and the mismatch term `log(e_j/3)` for every other copy is a **soft repulsion** —
SDA's ±1 edge as a continuous log-likelihood. A read spanning **no distinguishing** bubble has equal `L_{rk}`
across the tied copies: SDA's *un-attractable* read = the K-frontier.

**E-step = soft WhatsHap; M-step = path abundance.** With copy abundances `π = (π_1,…,π_K)` on the simplex:
E-step `γ_{rk} = π_k L_{rk} / ∑_j π_j L_{rj} = softmax_k(\log L_{rk}+\log π_k)` (a fractional read→copy
assignment, SDA's hard partition made soft); M-step `π_k = \tfrac1N ∑_r γ_{rk}` (**copy-path abundance** —
the quantity SDA has no parameter for, since it assembles rather than quantifies). Convergence: the
observed-data log-likelihood `ℓ(π) = ∑_r \log ∑_k π_k L_{rk}` is **non-decreasing** each sweep (the tested
monotone invariant), stopping at `Δℓ < ε(1+|ℓ|)`. `θ` is **fixed** from the VG (copy-path refinement —
re-estimating `θ_k` from γ-weighted reads, the direct analog of SDA re-clustering PSVs — is deferred).
**In one line:** the EM is SDA's PSV correlation-clustering with signed edges → substitution log-likelihood,
hard WhatsHap → soft responsibilities, plus a copy-abundance parameter.

## The model and its assumptions

One PSV-aware VG per family: bubbles `[m]`, `K` copy-paths `θ_k`. Read `r` has hidden origin `z_r ∼
Categorical(π)`, draws `obs(r)` from a coverage law independent of `z_r`, emits `r(j) ∼ q_j(\cdot\mid
(θ_{z_r})_j)`. Marginal = the **finite mixture** `P(r) = ∑_k π_k L_{rk}`. **Assumptions.** (A1) error-free
core / completeness: each read originates from one of the `K` VG copy-paths (the standing hypothesis of
Theorem 4; the reference-absent/O4 escape is handled by abstain-and-re-thread, not here). (A2) bounded known
per-column error `e_j < (|A_j|−1)/|A_j|` (matched allele is modal). (A3) coverage law independent of origin.
These are the combinatorial theory's assumptions given a probability measure.

## Theorem (EM consistency in the identifiable regime)

With `D_{ij} = {d:(θ_i)_d≠(θ_j)_d}` the distinguishing bubbles and `δ(r) = \min_{k≠b}|obs(r)∩D_{bk}|` (the §5b
`δ`; `b` = MLE copy), the per-read certificate is `min_p(r) = ε^{δ(r)}`, `Certified` iff `min_p(r) <
α/(K−1)`, `SoftZone` otherwise (`label_read`).

**Denominator caveat (do not conflate).** The shipped `label_read` uses `α/(K−1)` with `K` = number of copies
(per-copy Bonferroni over a read's `K−1` competitors); Theorem 4 above writes `α/(n−1)` with `n` = reads.
These are **not** the same symbol — reconcile which population each Bonferroni family is over before equating.

> Suppose the family satisfies **Strong Separation** (§5): every copy pair is separated by a bubble reads
> actually span, so `δ(r) ≥ 1` (`min_p < 1`) is achievable for every cross-copy comparison. Then:
> **(a) Identifiability** — with `θ` fixed the mixture components are **pre-labelled** (each `L_{k·}` tied to a
> named `θ_k`), so **no relabelling ambiguity**; `π` and per-copy assignments are uniquely determined.
> **(b) MLE consistency** — as per-copy coverage `n→∞`, `π̂_n → π*` almost surely, and the EM from a generic
> start converges to it. **(c) Assignment consistency** — for every identifiable read (`δ≥1`), the MAP
> `ẑ_r = \arg\max_k π̂_k L_{rk} → z*_r` for `n` large, and `γ_{r,z*_r} → 1`. **(d) The non-identifiable class
> is honest** — for `δ(r)=0` (`min_p=1`, SoftZone) the posterior stays at the prior `γ_{rk}=π_k` for all `n`;
> the EM makes **no hard call and never imposes a 1/k split**.

This is exactly the regime where the theory proves **Strong Separation ⟹ `C*` is the unique minimum copy
cover, `MCC = χ(H) = K`** (Lemma 1 + Theorem 2); the EM converges to that cover.

**Proof ingredients.** Finite-mixture identifiability + MLE consistency (Redner & Walker, *SIAM Review*
26:195–239, 1984; Teicher, *Ann. Math. Statist.* 34:1265–1269, 1963), specialized to the discrete PSV
emission. **Linear independence** is *not* by "distinct product-multinomials are independent" (false for
`K≥3`): instead the `K×K` matrix `M = [f_k(θ_j)]` with `f_k(θ_j)=∏_d q_d((θ_j)_d\mid(θ_k)_d)` is **strictly
diagonally dominant** in the well-separated regime (each `f_k` peaks at its own centre `θ_k`, off-diagonal
entries lose a factor `(1−e_d)→e_d/(|A_d|−1)` at each distinguishing column), hence invertible
(Levy–Desplanques) ⟹ `{f_k}` independent ⟹ mixture identifiable. **Global (not merely local) EM
convergence** holds **because `θ` is fixed**: then `ℓ(π)=∑_r \log∑_k π_k L_{rk}` is a sum of logs of functions
**affine in `π`**, hence **concave on the simplex**, and the mixture-proportion EM is coordinate/MM ascent to
the **global** MLE (unique under strict concavity, all `π*_k>0`). If `θ` were **also** estimated (the deferred
refinement), `ℓ(π,θ)` is non-concave and only the *local* Redner–Walker guarantee returns — the shipped
engine estimates `π` with `θ` fixed, so it lives in the globally-convergent concave case.

**Coverage attribution (what the limit does and does not buy).** The likelihood-ratio favouring the truth is
`L_{r,z*_r}/L_{r,k} = ((1−e)/(e/3))^{≥δ(r)} > 1` — a **per-read `δ/e` quantity that does NOT grow with
coverage** (`δ(r)` is fixed by which bubbles the single molecule spans). Under A1 (`e→0`) per-read assignment
is **exact at any coverage** (correct with `n=1` copy-read as with `n=100`). What `n→∞` buys is the
**abundance** consistency `π̂→π*` of (b); it enters (c) only through the prior factor `π̂_k/π̂_{z*_r}`, which
matters solely in the noisy `e>0` regime.

## Consequence — explaining SDA's 91–93% floor

SDA runs *hard* correlation clustering and **forces** a label on every PSV/read, including the un-attractable
K-frontier ones. Two costs: (i) the NP-complete heuristic can stall at a non-global clustering; (ii)
hard-calling the genuinely unidentifiable mass **must** misassign a fraction. The theorem decomposes this: on
the **identifiable** set (`δ≥1`) the ML soft relaxation is **consistent** — per-read accuracy high at every
coverage, `→100%` under A1; it does **not** floor. On the **unidentifiable** set (`δ=0`, K-frontier/K=0) no
method resolves the reads; SDA hard-calls and eats a misassignment rate whereas the EM **abstains**
(SoftZone). **What the theorem supplies:** the mechanism / *kind* of floor (a hard-caller scores ≈100% on the
identifiable set and pays a forced error on the unidentifiable residue, pinning *overall* accuracy below 100%).
It does **not** derive the specific **91–93%** — that value is the *instance-specific unidentifiable fraction*,
which the theory does not supply. Keep the two curves **distinct**: the EM's **identifiable-set** accuracy
`→100%` is a different measurement from SDA's **overall** 91–93%.

**Prediction (coverage sweep, `bench/em_coverage_sweep.py` → `bench/THEORY.md`).** On a planted sim
with known `θ*, π*, z*`, over per-copy coverage `{1,2,5,10,20,50,100}×`, two distinct axes: *(per-read `δ/e`,
not coverage-driven)* assignment accuracy on identifiable reads high at every coverage `→100%` under A1;
*(coverage axis)* abundance error `‖π̂−π*‖₁ → 0`; and K=0 families **stay SoftZone** at every coverage (the
boundary is a boundary, not a coverage artifact). These are the theorem's falsifiable content.

**Column filters (which bubbles enter the graph).** The raw allele-disagreement graph is error-inflated
(§5c: colouring ≈3× the true `K`), so columns are gated by the same per-column filters, reused unchanged:
**IsoCon** (`reference_isocon_sahlin`, per-position real-variant-vs-error), **SUN** (`reference_sudmant_2010`,
single-position private markers — a read over a SUN column is single-read pinned, the `δ`-contributing
attraction part (c) relies on), **Clair3-RNA** (`reference_clair3_rna`, flags A→I editing).
`em_assign_family` computes `detect_editing_columns` and passes it into `read_copy_evidence`, so an editing
column is downweighted in the EM likelihood as in the hard gate. **Scope caveat:** the shipped `--em` path is
**PSV-only** — it threads no copy-specific junctions and no per-base quality
(`ReadFeatures::junctions`/`psv_qual`, `CopyProfile::junctions` left empty), so its per-read labels can differ
from the hard `.assignments.tsv` gate on reads whose call depends on junction or per-base-quality evidence.
This does not affect the abundance/consistency result: junctions/quality/editing change only `min_p` and
per-column weights, not the shape of the abundance fixed point.

## Ties to the combinatorial theory (above)

- **Lemma 1 (MCC = χ(H)):** the EM's hard limit (`γ` → indicators) is a `χ(H)=K`-colouring of the conflict
  graph `H`.
- **Theorem 2 (unique cover under Strong Separation):** the EM's identifiable regime **is** this hypothesis;
  its MAP partition converges to the provably-unique cover `C*`.
- **Theorem 3 (`RECOVER`, poly-time):** the EM's hardened responsibilities reproduce the connected-component
  partition of the compatibility graph — the statistical route to the same object.
- **Theorem 4 (`min_p` bridge):** the per-read certificate the EM attaches via `label_read` (Certified =
  part (c) applies; SoftZone = part (d) abstention).
- **Theorem 7 (integrality bridge):** under Strong Separation the MWCA LP is integral; the EM's soft
  assignment converges to that same integral optimum — the EM is the mixture-model face of the
  `NP-hard (Thm 1/5) → (1−1/e) LP-rounding (Thm 6) → integral under Strong Separation (Thm 7) → per-read
  min_p (Thm 4)` object, not a third ad-hoc solver.
- **The K=0 floor = SDA's ">100 kb reads" = SoftZone:** when a copy pair shares every spanned bubble the
  mixture likelihood is flat (§6b), every read is `δ=0`, `min_p=1`, and the EM returns a posterior-invariant
  class (`γ=π`) — **never a hard 1/k**. The EM reaches the identifiability boundary and certifies it.

