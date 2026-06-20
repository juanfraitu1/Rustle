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

*See also: `bench/family_definition_formal.md` for the upstream family-detection step whose output populates
$R$ for each identified gene family.*
