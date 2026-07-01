# A Formal Definition of DNA-Level Multi-Copy Gene Families for the Gorilla Pan-Transcriptomic Thesis

**Scope.** This note fixes, for the rest of the thesis, what a **multi-copy gene family at the DNA level** is, as **one canonical combinatorial object** rather than a raw single-linkage catalog. It sits deliberately parallel to two companion notes and must compose with them without contradiction:

- `bench/SEGDUP_DEFINITION_FORMAL.md` — the **segmental-duplication predicate** `SD(·)` and the segdup **edge oracle** $E_a$ (plus $E_b$ exonic homology, $E_c$ read-ambiguity, and the lattice among them). *That note flags the raw connected components of $E_a$ as "**pre-refinement**", explicitly deferring the fix to "the ≥ 2-distinct-loci / homology-component refinement" (its §3.3, §3.4). **This note is that refinement, made a well-defined operator.***
- `bench/family_definition_formal.md` — the **RNA (transcript) family** $E_r$: a $\gamma$-quasi-clique-refined connected component of the **transcript-homology** graph $G_R=(V_R,E_r)$ (the fourth homology oracle, same skeleton and operator $R$ as this note). **The read-conflict graph $E_c$ (with the copy-assignment count $\mathrm{MCC}=\chi(H)$ and the SUN 3-tier ladder) is demoted there to the within-family O2 copy-assignment substructure**, with $E_c^{\mathrm{sig}}\subseteq E_c\subseteq E_r^{\mathrm{asym}}$ the sole clean containment. (Supersedes the older reading in `bench/FAMILY_DEF.md` §`family_definition_formal`, where the RNA family was $E_c$ itself.)

The combinatorial **skeleton** is shared across all three (nodes = loci, edge = a homology tie, family = connected component with $|C|\ge 2$); **the only thing that differs is the edge oracle** (`SEGDUP_DEFINITION_FORMAL.md` §3.1). What this note adds is the missing **refinement operator $R$** that turns a raw $E_a$ component — which single-linkage **over-merges** — into the canonical family, together with its provable properties, the multi-copy predicate, and the honest naming of the one irreducible cohesion parameter.

All empirical statements were re-verified against `/home/juanfra/winloci_scratch/GGO_genomic.gff` (34,114 RefSeq gene loci) and `/mnt/c/Users/jfris/Desktop/final.bed` (SEDEF self-alignment) via `bench/genome_family_def.py` on 2026-06-30; the load-bearing counts and the machine-checked properties are inlined and reproducible (Appendix A).

> **One-line thesis claim.** A DNA-level multi-copy gene family is a **cohesive community** of the genome segdup graph $G_{\mathrm{DNA}}=(V,E_a)$ satisfying the **cohesion certificate**: a block $C$ produced by the cohesion-refinement operator $R$ from a connected component of $G_{\mathrm{DNA}}$, subject to $C$ being a $\gamma$-quasi-clique (internal edge density $\ge\gamma$) with **$\ge 2$ distinct genomic loci**. The **certificate is the canonical object** — a *property* every emitted family provably carries; exact max-$\gamma$-quasi-clique partition is NP-hard, so $R$ returns **one certificate-passing *witness* partition** (seed-fixed), not a unique family set (§3.2, §6). It is the **same skeleton, different edge oracle** as the RNA read-conflict O2 substructure ($E_a$ read-independent vs $E_c$ read-derived; the RNA *homology* family is $E_r$, the true $E_a$ parallel — §5.2); $R$ is the DNA-side cure for a **transitive-closure** over-merge that — unlike the RNA side's per-edge significance gate — **cannot** be fixed by bridge removal (only 0.3 % of the blob's edges are bridges), and so is **necessarily community/clique-based**.

---

## 1. The raw object, and why it is not yet the canonical family

### 1.1 The genome segdup graph $G_{\mathrm{DNA}}$

Fix the haploid T2T genome $G$ = GCF_029281585.2. From `SEGDUP_DEFINITION_FORMAL.md` §3.1 (Layer 1, Catalog A), define

$$
G_{\mathrm{DNA}}=(V,\,E_a),\qquad V=\{\text{NCBI RefSeq gene loci}\}\ (|V|=34{,}114),
$$
$$
\{u,v\}\in E_a \iff \exists\,(I,J)\in S:\ \mathrm{SD}(I,J)\ \wedge\ \mathrm{cov}(u,I)\ge \tfrac12\ \wedge\ \mathrm{cov}(v,J)\ge \tfrac12,\qquad \mathrm{cov}(g,X)=\tfrac{|g\cap X|}{|g|},
$$

where $S$ is the SEDEF self-alignment and $\mathrm{SD}(\cdot)$ is the segdup predicate of the companion note. **Two honestly-named oracle settings** (this note is agnostic to which; both are supported):

- **default (raw ~50 %-floor superset)** — $\mathrm{SD}(\cdot)$ relaxed to SEDEF's on-disk floor ($\rho\gtrsim 0.5$, $\ell\ge 1$ kb); chosen for **recall** of divergent paralog families (CEACAM/KRAB-ZNF/PRSS/IFITM/ULBP all sit **below** the Bailey 90 % cliff and are lost under strict $\mathrm{SD}(\cdot)$; the APOBEC3D/F = 88.4 % witness generalizes).
- **strict** (`--bailey-sedef`) — the Bailey `SD(·)`: $\rho\ge 0.90\wedge\ell\ge 1$ kb $\wedge\ \neg\mathrm{TE}$.

The **raw DNA family** is a connected component of $G_{\mathrm{DNA}}$ with $|C|\ge 2$, computed by union-find (`genome_family_def.py:build`). On the default oracle: **623 raw families / 5,162 member genes / 261 cross-chromosome** (91,247 gene–gene edges).

### 1.2 The over-merge is a dense transitive blob, not a bridge chain (load-bearing)

Single-linkage takes the **transitive closure** of $E_a$. In a segdup-rich region this collapses biologically unrelated loci into one giant component. Verified on the **largest** raw component:

| property | value |
|---|---|
| genes | **1,547** (885 protein-coding, 416 lncRNA, 81 snRNA, 70 V_segment, 56 rRNA, 24 snoRNA, 15 tRNA) |
| edges | 24,431 |
| internal density $\rho_{\mathrm{in}}$ | **0.0204** |
| **graph bridges** | **71 (0.3 % of edges)** |
| median degree | 23 |
| content | DGCR6, PRODH, DRD5, SLC6A8, OTOP1, BCAP31, ANKRD36C, TPTE2, … — the **22q11 / pericentromeric** segdup-rich region |

**Consequence (this reshapes the refinement).** The over-merge is a **dense heterogeneous blob**: unrelated loci (DGCR6 $\leftrightarrow$ DRD5 have **no direct edge**) are chained only *transitively* through dense real $\ge 50\%$-homology edges. Therefore:

- **Bridge removal / 2-edge-connectivity does NOT fix it** — only 71 of 24,431 edges are bridges; cutting them leaves the blob essentially intact.
- **`SD(·)` alone does NOT fix it either** — imposing the Bailey floors removes the *repeat-bridged* pseudo-edges (TRNAV-CAC ×173, rRNA ×70) but a **mosaic SD** sharing a valid $\ge 90\%$ block with family $X$ and another with $Y$ still chains $X$–$Y$ through two *individually valid* edges (max family still 317 under `--bailey-sedef`; `SEGDUP_DEFINITION_FORMAL.md` §3.3: "`SD(·)` necessary but not sufficient").

The cure must operate on **cohesion**: a real family (CEACAM5/6/7 = a size-10 near-clique; a KRAB-ZNF cluster) is **densely mutually homologous**, whereas the blob is a low-density union of such clusters glued at mosaic hubs. This is a **community/clique** criterion, and the cohesion cut is an **irreducible parameter** (named honestly in §6 — *not* threshold-free).

---

## 2. The canonical object: the cohesion certificate (and its seed-fixed witness catalog)

> **Definition (DNA-level multi-copy gene family).** Let $\gamma\in(0,1]$ be the cohesion parameter. A **DNA-level multi-copy gene family** is a vertex set $C\subseteq V$ such that
> $$
> C\in R(\kappa)\ \text{ for the connected component }\kappa\text{ of }G_{\mathrm{DNA}}\text{ that contains it,}\qquad \text{and}\qquad \boxed{\ \Lambda(C)\ge 2\ }
> $$
> where $R$ is the cohesion-refinement operator of §3 (so $C$ is a **$\gamma$-quasi-clique**: internal density $\rho_{\mathrm{in}}(C)\ge\gamma$, or $|C|\le 2$) and $\Lambda(C)$ is the **distinct-locus count** of §4 (the multi-copy predicate). The predicate "$C$ is a $\gamma$-quasi-clique with $\Lambda(C)\ge 2$ inside a raw $E_a$ component" is the **cohesion certificate**, and **it — not any one partition — is the canonical object**: exact maximum-$\gamma$-quasi-clique partition is NP-hard (§6), so $R$ is a heuristic that returns **one certificate-passing *witness* partition** (deterministic under the disclosed splitter seed $=0$; §3.2). A concrete witness is emitted to the **refined catalog** (`genome_families_refined.tsv`, prefix `GRFAM`) under the **opt-in `--refine` flag**; a default run (flag absent) emits only the raw single-linkage catalog and is **byte-identical** to the pre-refinement output (the default `*.tsv`/`.json` are untouched, and the refined catalog is written to its own file), exactly mirroring the RNA engine's byte-identical-when-OFF `--refine` discipline (§5.2). The raw catalog is retained as the airtight, refinement-free baseline.

$\rho_{\mathrm{in}}(C)=\dfrac{2\,|E_a(C)|}{|C|\,(|C|-1)}$ is the induced edge density of the subgraph $G_{\mathrm{DNA}}[C]$ ($=1$ by convention for $|C|\le 1$). A $\gamma$-quasi-clique is the standard cohesive-subgraph notion: at $\gamma=1$ it is a clique; the requirement encodes "densely mutually homologous," which the transitive blob fails and CEACAM5/6/7 passes.

**Result (default oracle, $\gamma=0.20$, seed $=0$).** $623$ raw components $\Rightarrow$ **708 refined families / 5,136 member genes / 316 cross-chromosome** in the seed-0 witness; the certificate holds **$708/708$**. Of these, **$614$ are splitter-invariant** — raw components already dense enough that the step-1 gate emits them **whole, before the splitter runs** (the five verified real families are all among them, §3.4) — and only **$9$ low-density mega-components are actually split by the seed-dependent heuristic**. The single 1,547-gene blob $\Rightarrow$ **48 cohesive families** *in the seed-0 witness* (sizes 207, 189, 104, 78, …); this blob-derived **count is a witness artifact, not an invariant** — it is $48/41/39$ at seeds $0/1/7$ (§3.2, §6). What every seed shares is the **certificate** (every emitted family a $\gamma$-quasi-clique with $\ge2$ loci) and the $614$-family invariant core. Every one of the five verified real families **survives whole under every seed** (§3.4), because it is gate-emitted before the splitter.

Each family is, invariantly, a **sub-object of a raw $E_a$ component** ((R1), seed-independent) — exactly the "post-refinement" reading the segdup note demands.

---

## 3. The refinement operator $R$

### 3.1 Definition

$R$ acts on a raw component $\kappa$ (a node set) and returns a **set of node-set blocks** — a partition of $\kappa$:

$$
R:\ 2^{V}\longrightarrow 2^{2^{V}},\qquad R(\kappa)=\text{leaves of the density-gated recursive split of }G_{\mathrm{DNA}}[\kappa].
$$

**Recursion (density-gated).** For a block $B$ (initially $B=\kappa$):

1. **Certificate gate.** If $|B|\le 2$ or $\rho_{\mathrm{in}}(B)\ge\gamma$, **stop** — emit $B$ as a leaf (it is a $\gamma$-quasi-clique). *This gate is evaluated on the whole block first, so a cohesive family is never split — the fix for modularity's over-splitting (§3.4).*
2. **Split.** Otherwise partition $B$ by a **guaranteed-progress splitter** $\sigma(B)$ (§3.2) into $\ge 2$ strictly smaller blocks, and recurse $R$ on each.

$R$ applied to a whole catalog splits every raw component and then keeps only leaves passing the multi-copy predicate $\Lambda(\cdot)\ge 2$ (`genome_family_def.py:refine_families`).

### 3.2 The splitter $\sigma$ (a heuristic witness-generator; the *certificate* is the clean object)

$\sigma$ is a **guaranteed-progress** splitter: adaptive-resolution modularity (Louvain, resolution $\in\{1,2,4,8\}$, seed $=0$), falling back to Kernighan–Lin bisection (seed $=0$), then deterministic halving — so on any $B$ with $|B|>2$ it returns $\ge 2$ blocks each strictly smaller than $B$.

**Honest separation of object from algorithm (Canzar-style).** The clean, well-defined object is the **certificate *property*** — a leaf is admitted iff it is a $\gamma$-quasi-clique with $\ge2$ loci — and *this property* is **$\sigma$-invariant**: every splitter that makes progress yields an $R$ whose every leaf passes the certificate. What is **not** $\sigma$-invariant is the **specific witness partition**. Because finding a *maximum* $\gamma$-quasi-clique partition is **NP-hard** (§6), $\sigma$ can only reach **one certificate-passing partition among exponentially many**, and different seeds carve a sub-$\gamma$ region differently. This is **measured, not hypothetical**: re-refining the 1,547-gene blob at seeds $0/1/7$ gives **$48/41/39$** families, and only **$22$ of the $48$** seed-0 blob families reproduce exactly at seed $1$ — *yet the certificate holds at every seed* ($41/41$, $39/39$). So:

- **$\sigma$-invariant** (canonical): the certificate; the $614$ **gate-emitted** components (step-1 density $\ge\gamma$ before $\sigma$ ever runs), which **include all five verified real families** (§3.4); (R1) each family $\subseteq$ one raw component; idempotence and $\gamma$-monotonicity.
- **$\sigma$-dependent** (witness artifact): the exact partition of the $\approx9$ **sub-$\gamma$ mega-components**, hence the blob-family **count** ($\sim39$–$48$), the total refined **count** ($708$ at seed $0$; $699$ at seeds $1,7$), and the `GRFAM` index labels. These are reported as the **seed-0 witness**, not as invariants.

We fix seed $=0$ to make the emitted witness **deterministic and byte-identical** (Appendix A), and disclose it as a load-bearing constant (§6). (Contrast with using **raw** modularity as the operator: at resolution 1 it over-splits CEACAM $12\to[8,4]$ and ZNF716 $18\to[7,7,4]$; the step-1 certificate gate is exactly what vetoes those splits — that veto is $\sigma$-invariant, the carving of the blob is not.)

### 3.3 Properties (all machine-checked green, `bench/genome_family_def.py` + Appendix A)

Let $R(\kappa)$ denote the leaf set for component $\kappa$; write $R$ also for its catalog-wide lift.

- **(R1) Partition / sub-components of the raw component.** $R(\kappa)$ is a partition of $\kappa$; in particular every family $C\in R(\kappa)$ satisfies $C\subseteq\kappa$. So a refined family is always **contained in exactly one raw $E_a$ component** — $R$ only ever *splits*, never merges or crosses component boundaries. *Machine-checked: each of the 708 blocks $\subseteq$ one raw component = True.*
- **(R2) Certificate.** Every emitted block is a $\gamma$-quasi-clique or has $|B|\le 2$: $\forall B\in R(\kappa),\ |B|\le 2 \ \vee\ \rho_{\mathrm{in}}(B)\ge\gamma$. Guaranteed by step 1 termination + $\sigma$'s progress. *Machine-checked: 708/708 blocks $\gamma$-cohesive-or-$\le 2$ = True.*
- **(R3) Idempotent.** $R\circ R=R$: re-refining an already-refined family returns it unchanged, because every leaf passes the certificate gate (step 1) on the first re-entry and is emitted immediately. *Machine-checked: $R(R(\text{catalog}))$ = $R(\text{catalog})$ as sets, $708\to708$ identical = True.*
- **(R4) Monotone (refining) in $\gamma$.** For $\gamma'\ge\gamma$, $R_{\gamma'}$ **refines** $R_{\gamma}$: every $\gamma'$-block is contained in some $\gamma$-block (a stricter cohesion cut can only subdivide). *Machine-checked: $R_{0.30}$ refines $R_{0.15}$ = True; $|R_{0.15}|=678 \le |R_{0.20}|=708 \le |R_{0.30}|=753$.* (This is a **finiteness / consistency** property, not a claim that $\gamma$ is free.)
- **(R5) Termination.** $\sigma$ strictly decreases block size and step 1 halts at $|B|\le 2$; the recursion tree has depth $\le|\kappa|$, so $R$ terminates. Worst case leaves are size-$\le 2$ pairs (trivially $\gamma$-cohesive), so (R2) holds unconditionally.

> **Proposition (canonical certificate).** For any raw component $\kappa$ and $\gamma\in(0,1]$, $R(\kappa)$ is a partition of $\kappa$ into blocks each of which is a $\gamma$-quasi-clique or has size $\le 2$; $R$ is idempotent and $\gamma$-monotone. *(R1)–(R5), machine-checked.*

### 3.4 The verified real families survive whole — and are $\sigma/$seed-invariant (default oracle, $\gamma=0.20$)

The five families the refinement was required to preserve — plus two bonus textbook families — are emitted as **single** cohesive communities. All seven are **their own raw $E_a$ component** with density $\ge\gamma$, so the **step-1 gate emits them whole *before* the splitter $\sigma$ runs**; hence their preservation is **seed-independent** (verified at seeds $0/1/7$), unlike the blob partition of §3.2. (`GRFAM` ids are seed-0 witness labels; the families sort ahead of the flagged mosaic hubs of §3.5.)

| family | id (seed-0) | genes | distinct loci | $\rho_{\mathrm{in}}$ | biotype purity | note |
|---|---|---:|---:|---:|---:|---|
| CEACAM5/6/7 | GRFAM75 | 12 | 12 | 0.470 | 0.833 | near-clique, whole |
| KRAB-ZNF (ZNF716…) | GRFAM54 | 18 | 18 | **0.261** | 0.611 | divergent, whole (**lowest** verified real density — sets the $\gamma$ ceiling) |
| IFITM1/2 | GRFAM138 | 7 | 7 | 0.524 | 1.000 | whole |
| ULBP1 | GRFAM133 | 7 | 7 | 0.667 | 1.000 | whole |
| MAGEA12 / CSAG1 | GRFAM91 | 10 | 10 | 0.733 | 0.800 | **kept together**, whole |
| PRSS2/3 (trypsinogen) | GRFAM57 | 17 | 17 | 0.316 | 0.529 | whole (see note) |
| LILRA2/5 | GRFAM172 | 6 | 6 | 0.467 | 1.000 | whole |

They sit in a density band **$[0.261,\,0.733]$**, while the transitive blob is $0.0204$ — a **$\ge 13\times$ gap** that $\gamma=0.20$ occupies (§6). None is part of the 1,547-gene blob, so the refinement is **$\approx 99\%$ inert (614 / 623 raw components pass through unchanged)**, acting only on the **9** low-density mega-components. The genuinely conflated loci are inside those 9: the 22q11 / pericentromeric anchors DGCR6, PRODH, DRD5, TPTE2, BCAP31, CEP170 (all verified blob members) are transitively glued into the blob and are exactly what $R$ separates into distinct cohesive communities.

> **Honest disclosure (calibration = the survival set; the $\gamma$ knife-edge).** $\gamma=0.20$ was **calibrated on exactly these seven families** — there is **no held-out family**, so "they survive whole" is guaranteed **by construction**, not an independent validation pass. It is genuinely **$\gamma$-sensitive**: the lowest real family, KRAB-ZNF at $\rho_{\mathrm{in}}=0.261$, is its own component and is gate-emitted whole **iff $\gamma\le 0.261$**, so $\gamma=0.20$ preserves it with a **margin of only $0.061$**, and **any $\gamma>0.261$ splits the KRAB-ZNF family** (empirically $\gamma=0.30\Rightarrow[7,7]$). The preservation guarantee is therefore a **calibrated, not free, choice**; what makes $\gamma=0.20$ defensible is the **$\ge13\times$ empirical gap** below (§1.2, §6), not a claim of uniqueness.

*Note on GRFAM57 (PRSS).* Its named members are **PRSS2, PRSS3** (trypsinogen serine proteases, both sub-90 % divergent) together with the two flanking genes **TMEM45B, UBE2R2** co-covered by the same trypsinogen segdup block, plus 13 LOC pseudogenes; **PRSS1 is not present** in this refined family under the current annotation. Its biotype purity $0.529$ sits just above the $0.50$ coherence line (§3.5); it is reported as a textbook trypsinogen segdup rather than a pure-trypsinogen orthogroup — an honest note, not a purity claim.

### 3.5 The certificate is necessary but **not sufficient**: the repeat-mosaic residual

A $\gamma$-quasi-clique with $\ge 2$ loci is **necessary** for a multi-copy family but **not sufficient** for biological coherence. In the seed-0 witness, **9 of the 708** families clear the cohesion certificate yet are **repeat-mosaic hubs** — dense enough to pass $\gamma$, but a mix of *unrelated biotypes* spread across many contigs, i.e. co-clustered by shared interspersed/segmental-repeat content rather than by being copies of one gene. The two largest:

| id (seed-0) | genes | contigs | $\rho_{\mathrm{in}}$ | biotype purity | biotype mix |
|---|---:|---:|---:|---:|---|
| GRFAM699 | 207 | 23 | 0.265 | **0.435** | 90 protein_coding / 85 lncRNA / 32 snRNA |
| GRFAM700 | 104 | 16 | 0.226 | **0.481** | 50 protein_coding / 38 lncRNA / 16 snRNA |

**We make the insufficiency auditable in the catalog itself**, not just in prose. Every refined family carries two coherence descriptors: `biotype_purity` (dominant-biotype fraction $\max_b|\{i\in C:\mathrm{biotype}(i)=b\}|/|C|$) and a boolean `repeat_mosaic` flag = $[\,\text{purity}<0.50\,]$ — "no biotype is a majority," so it cannot be copies of a single gene. **9/708** families are flagged; they are **sorted to the *end*** of `genome_families_refined.tsv` (so the incoherent hubs no longer sit at the top) but are **not filtered out** — the certificate object is left exactly as defined; the flag is a *descriptor*, not a gate. All seven verified real families pass (purity $0.529$–$1.000$; §3.4).

This residual is the **TE/repeat cut** the segdup note discloses, and it is **oracle-removable but not for free**: switching to the strict Bailey $\mathrm{SD}(\cdot)$ oracle (`--bailey-sedef`, which excludes high-copy interspersed repeats) shrinks the catalog to **489 families / max 124** but still leaves a large residual hub (12 `repeat_mosaic`-flagged) *and* drops the divergent real families CEACAM/KRAB-ZNF/PRSS/IFITM/ULBP (all sub-90 %). So **neither operating point fully cleans the catalog** — the honest statement is: the certificate is necessary; the `repeat_mosaic` column and the raw-vs-strict oracle are the two disclosed handles on the remaining insufficiency (`PURITY = 0.50` is itself a named parameter, §6).

---

## 4. The $\ge 2$-distinct-loci multi-copy predicate

A "family" must have **$\ge 2$ physical copies**, not one locus wearing two overlapping annotations (e.g. in GRFAM1, `SLC25A15` at `NC_073238.2:55,065,401–55,137,826` and the nested `LOC109029455` at `55,072,914–55,134,160` reciprocally overlap $\gg 50\%$ and count as **one** locus). Formalize the **distinct-locus count** $\Lambda(C)$ by collapsing co-located annotations:

$$
u \approx_{\mathrm{loc}} v \iff \mathrm{chrom}(u)=\mathrm{chrom}(v)\ \wedge\ \frac{|u\cap v|}{|u|}\ge \tau_{\ell}\ \wedge\ \frac{|u\cap v|}{|v|}\ge \tau_{\ell},\qquad \tau_{\ell}=0.50\ (\text{reciprocal overlap}),
$$
$$
\Lambda(C)=\big|\,C/\!\approx_{\mathrm{loc}}^{*}\,\big|\quad(\text{number of connected components of }\approx_{\mathrm{loc}}\text{ restricted to }C).
$$

**Multi-copy predicate:** $C$ is admitted as a family iff $\Lambda(C)\ge 2$ (`genome_family_def.py:distinct_loci`). This is the DNA analog of the RNA family's $|C|\ge 2$ **distinct-locus** requirement (the segdup note's "$\ge 2$-distinct-loci" clause). On the default run it removes the residual single-physical-locus blocks (5,162 raw member genes $\to$ 5,136 after refinement + $\Lambda\ge 2$). *Machine-checked: all 708 refined families have $\Lambda\ge 2$ = True.*

---

## 5. Composition with the segdup predicate $\mathrm{SD}(\cdot)$ and relation to the RNA O2 read-conflict substructure $E_c$ / $\mathrm{MCC}=\chi(H)$

### 5.1 Composition with $\mathrm{SD}(\cdot)$ (the segdup lattice is preserved)

$G_{\mathrm{DNA}}$ is built **on** the segdup edge $E_a$, which is built **on** the predicate $\mathrm{SD}(\cdot)$ (`SEGDUP_DEFINITION_FORMAL.md` §3.1). $R$ acts **inside** the $E_a$ circle: since $R$ only splits components (R1), a refined family is a **subset of a raw $E_a$ component**, so $R$ **does not disturb the segdup lattice**:

$$
E_c^{\mathrm{sig}}\subseteq E_c\subseteq E_r^{\mathrm{asym}}\subseteq_{\mathrm{proj}} E_b^{\mathrm{asym}}\qquad\text{(the containment chain; }E_c\subseteq E_r^{\mathrm{asym}}\text{ is now tightest, sharing the node set }V_R\text{)};\qquad E_a\ \text{incomparable to }E_b,E_r,E_c.
$$

The tightest clean containment moved down one level to the transcript oracle $E_r$ (`bench/family_definition_formal.md` §4, `bench/SEGDUP_DEFINITION_FORMAL.md` §3.3 Update); $E_c\subseteq E_b^{\mathrm{asym}}$ stays valid by projection. $R$ acting inside the $E_a$ circle does not disturb any of this.

$R$ is precisely the operator the segdup note **defers to** in §3.3/§3.4 ("`SD(·)` is necessary but not sufficient… the ≥ 2-distinct-loci / homology-component refinement is the rest"; "the formal components should be read as *post-refinement*, not raw union-find"). Composition is exact: **$R$ takes the pre-refinement $E_a$ components and returns the canonical families**, at either the raw-superset or the strict-$\mathrm{SD}(\cdot)$ oracle setting.

### 5.2 Relation to the RNA read-conflict O2 substructure — **same skeleton, different edge oracle**

> **Reframe alignment (2026-06-30).** Post-`bench/family_definition_formal.md`, the RNA *homology family* (O1) is the transcript-homology component $E_r$ — the true homology-oracle **parallel of $E_a$** (both are read-independent homology graphs refined by the same $R$). The read-conflict graph $E_c$ is **not** the RNA family; it is the within-family **O2** copy-assignment substructure ($E_c^{\mathrm{sig}}\subseteq E_c\subseteq E_r^{\mathrm{asym}}$). The comparison below therefore pairs $E_a\leftrightarrow E_c$ **only at the over-merge/refinement-failure-mode level** (both are graphs that over-merge and need refinement); the homology-oracle correspondence proper is $E_a\leftrightarrow E_r$. Read "RNA read-conflict substructure" for "RNA family" throughout this subsection.

The RNA **O2 read-conflict substructure** (inside a fixed $E_r$ family; `bench/family_definition_formal.md` §5) is the connected component of the read-conflict graph $G_c=(V,E_c)$, $|C|\ge 2$, where an edge is a read-**ambiguity** tie (de-tie / significance-tie). The DNA family is the connected component of $G_{\mathrm{DNA}}=(V,E_a)$, refined by $R$. The **skeleton is identical** (`SEGDUP_DEFINITION_FORMAL.md` §3.1: `conflict_families()` and `build` are the *same* union-find); the **only** difference is the edge oracle:

| | RNA read-conflict substructure (O2) | DNA family (this note) |
|---|---|---|
| edge oracle | $E_c$ — read-ambiguity (read-derived) | $E_a$ — genome self-alignment (read-independent) |
| raw object | component of $E_c$ | component of $E_a$ |
| **over-merge failure mode** | repeat-bridge / mosaic **reads** | **transitive closure** over dense real $E_a$ edges |
| **refinement** | **per-edge** significance gate $E_c^{\mathrm{sig}}\subseteq E_c$ ($\min_p=\varepsilon^\delta$, Thm 4) — removes weak *edges* | **whole-community** $R$ ($\gamma$-cohesion) — bridge/edge removal **can't** (0.3 % bridges) |
| multi-copy predicate | $|C|\ge 2$ distinct loci | $\Lambda(C)\ge 2$ distinct loci (§4) — **shared** |

The two refinements **differ because the failure modes differ**: the RNA over-merge is a few spurious *edges* (repeat-bridged reads), curable by a per-edge statistical gate; the DNA over-merge is a *dense transitive region* with almost no bridges, curable only by a *community* cut. Both then apply the **same** $\ge 2$-distinct-loci predicate. This is the honest DNA$\leftrightarrow$RNA correspondence: **one skeleton, two oracles, two failure modes, two refinements, one multi-copy predicate.**

Both refinements also share the same **engineering discipline** — the parallel is the **opt-in / byte-identical-when-OFF** contract, not a claim that the two operators are the same. On the RNA side there are in fact **two distinct** refinements, both opt-in: the **per-edge significance gate** $E_c^{\mathrm{sig}}\subseteq E_c$ of the comparison table (the $\min_p=\varepsilon^\delta$ tie-gate, toggled by the env var `RUSTLE_CONFLICT_SIG=1`), and the separate **exon-sum homology** refinement (`gw_family_catalog --refine`, `denovo_pipeline::refine_families_exon_sum`). The *edge-oracle* analogue of the DNA $R$ is the **first** (an edge-level cohesion gate); what the DNA `genome_family_def.py --refine` copies from **both** is the discipline: refinement is an opt-in step whose *absence* leaves the raw catalog **byte-identical** to the unrefined pipeline (verified in Appendix A). The refined catalog is the seed-fixed **witness** of the certificate-defined canonical object, while the raw single-linkage catalog is preserved unchanged as the airtight baseline, so every claim about $R$ is checked against a fixed reference.

### 5.3 Relation to $\mathrm{MCC}=\chi(H)$

On DNA there are no reads, hence no conflict hypergraph $H$ and no chromatic number — $\chi(H)$ is an **RNA-side** quantity (copy-assignment cover, `FAMILY_DEF.md` Lemma 1). The DNA family instead supplies the **reference copy number**
$$
K_{\mathrm{DNA}}(C)=\Lambda(C)\quad(\text{distinct loci in the refined family}).
$$
When a DNA family $C_{\mathrm{DNA}}$ and an RNA family $C_{\mathrm{RNA}}$ coincide on their shared loci, the two are linked by a **clean inequality**:
$$
\boxed{\ \mathrm{MCC}=\chi(H)\ \le\ K_{\mathrm{DNA}}(C)\ }
$$
— RNA can resolve **at most** as many copies as the genome contains, and the **K-frontier** is exactly the strict case $\chi(H)<K_{\mathrm{DNA}}$: read-identical copies (SEDEF $\rho\to 1$, PSV density $\to 0$; RFPL4A/DAZ) that the genome distinguishes but reads cannot, forcing abstention (`SEGDUP_DEFINITION_FORMAL.md` §4.2). So $R$ furnishes the $K$ against which $\chi(H)$ is measured; the DNA family is the **ground-truth copy count**, the RNA O2 resolution $\chi(H)$ the **read-resolvable** subset (inside the $E_r$ family), and their agreement region $E_a\cap E_c$ (RABL2A/B triple-core) is where both oracles see the same copies.

---

## 6. Honest irreducible parameters (do **not** claim threshold-freedom)

Four parameters are genuinely irreducible, plus one **disclosed splitter seed**; naming them all is the rigor Canzar rewards.

1. **COVER $=0.50$** — the gene-node **projection** fraction: a gene is "in" a segdup side iff $\ge 50\%$ of its span is covered. Inherited unchanged from $E_a$ (`SEGDUP_DEFINITION_FORMAL.md` §2.3 knob 5); orthogonal to cohesion, not dissolved here.
2. **The edge oracle: raw ~50 %-floor superset vs strict $\mathrm{SD}(\cdot)$** (`--bailey-sedef`) — a **recall/precision** choice on the SEDEF identity floor, *not* threshold-free. Default = raw (keeps sub-90 % divergent families CEACAM/KRAB-ZNF/PRSS/IFITM/ULBP); strict = Bailey $\rho\ge 0.90$ (cleaner repeat-wise: catalog $489$ / max $124$, but drops those families). This is the same honest cut the segdup note carries; $R$ runs on either.
3. **The cohesion parameter $\gamma$ (default $0.20$)** — the minimum internal density of a family (the $\gamma$-quasi-clique cut). It is **the irreducible community/clique cut**: there is **no threshold-free boundary** between "densely mutually homologous family" and "transitively glued blob." It is chosen — honestly, empirically — in the **gap** between the verified-real-family density band $[0.261,0.733]$ (§3.4) and the transitive-blob density $0.0204$ (§1.2), a $\ge 13\times$ gap. **Knife-edge (disclosed).** $\gamma=0.20$ was calibrated **on** the seven flagship families (its calibration set == its "survives-whole" set, §3.4 — no held-out family), and it sits just below the **lowest** real family, KRAB-ZNF at $0.261$: the preservation **margin is only $0.061$**, and **any $\gamma>0.261$ splits KRAB-ZNF** ($\gamma=0.30\Rightarrow[7,7]$). By (R4) the catalog varies **monotonically and predictably** with $\gamma$ ($678\!\to\!708\!\to\!753$ families over $\gamma\in\{0.15,0.20,0.30\}$), so $\gamma$ is a *calibratable* cohesion level — **not** a hidden knob, but a knob whose defensibility rests entirely on the empirical $\ge13\times$ gap, not on threshold-freedom.
4. **The coherence descriptor `PURITY` $=0.50$** (§3.5) — the dominant-biotype fraction below which a certificate-passing family is flagged `repeat_mosaic`. A **reporting/annotation** parameter (it never filters a family, only sorts flagged hubs to the end); it exposes the certificate's necessary-but-not-sufficient gap. Named, not hidden.

A minor nuisance parameter, $\tau_\ell=0.50$ (reciprocal-overlap for the $\approx_{\mathrm{loc}}$ locus-collapse, §4), is stable and orthogonal.

**The splitter seed $=0$ is a disclosed constant, not a free parameter of the object — but it *is* load-bearing for the witness.** Reconciling §3.2 with the NP-hardness of maximum $\gamma$-quasi-clique partition: the **certificate property** (every family a $\gamma$-quasi-clique with $\ge2$ loci) is **seed-invariant**, and so are the $614$ gate-emitted components (including all five real families) and properties (R1)–(R5). But the **specific witness partition** is **not** seed-invariant: the blob yields $48/41/39$ families at seeds $0/1/7$, only $22$ of the $48$ seed-0 blob families reproduce at seed $1$, and the total count is $708$ (seed 0) vs $699$ (seeds 1, 7). We therefore fix seed $=0$ to emit **one deterministic, byte-identical witness** and report **that** witness's counts (708, the 48-blob split, the `GRFAM` ids) as witness artifacts, never as invariants.

**What is NOT claimed.** $R$ is **not** threshold-free; the certificate *object* (a $\gamma$-quasi-clique with $\ge 2$ loci inside a raw component) is clean and provable, but $\gamma$ and the raw-vs-strict oracle are irreducible, `PURITY` is a named annotation cut, and the exact family *set* is a seed-fixed witness of an NP-hard partition. Every emitted family carries its density as a verifiable certificate (`genome_families_refined.tsv` column `density`) and its coherence as `biotype_purity`/`repeat_mosaic`. Whether a large certificate-passing community (e.g. the size-207 pericentromeric hub GRFAM699, $\rho_{\mathrm{in}}\ge\gamma$ but `repeat_mosaic`-flagged) is *one biological family* is a $\gamma$-relative, seed-relative statement, disclosed as such, not a monophyly claim.

---

## Appendix A — verification & reproduction log (2026-06-30)

Reproduce the **canonical refined** catalog (seed-0 witness): `/home/juanfra/miniforge3/bin/python bench/genome_family_def.py --threads 12 --refine` (add `--gamma G`, `--seed S`, `--bailey-sedef`). Reproduce the **raw baseline** (byte-identical to the pre-refinement output): drop `--refine` (the default). The refined output is byte-identical across `--threads` and across repeat runs at fixed `--seed`; the blob-derived family *count* changes with `--seed` (see the seed-sensitivity row).

| check | source / command | result |
|---|---|---|
| raw catalog (default oracle) | `genome_family_def.py:build` | 623 families / 5,162 genes / 261 cross-chrom (91,247 edges) |
| largest raw component | `networkx` on $E_a$ | 1,547 genes (885 PC, 416 lncRNA, 81 snRNA, 70 V_seg, 56 rRNA); 24,431 edges; $\rho_{\mathrm{in}}=0.0204$; **71 bridges (0.3 %)**; median deg 23 |
| bridge removal insufficient | 71 / 24,431 bridges | confirmed: cutting all bridges leaves blob intact |
| refined catalog ($\gamma=0.20$, **seed-0 witness**) | `refine_families` | **623 raw $\to$ 708 refined / 5,136 genes / 316 cross-chrom**; blob $\to$ 48 cohesive families (207,189,104,…) — the counts are a seed-0 witness of an NP-hard partition (see seed-sensitivity row) |
| **refinement is $\approx 99\%$ inert** | map $R$ leaves $\to$ raw components | **614 / 623** raw components emitted **unchanged** (already $\gamma$-quasi-cliques); only **9** low-density mega-components are split — the five verified real families are among the 614 (their own components, never in the blob) |
| **(R2) certificate** | `genome_family_def.py` in-run check | **708/708** blocks $\gamma$-cohesive-or-$\le 2$ = True |
| **(R1) partition-refinement** | `props` machine-check | each block $\subseteq$ exactly one raw component = True |
| **(R3) idempotent** $R\circ R=R$ | `props` machine-check | $708\to708$, identical as sets = True |
| **(R4) monotone in $\gamma$** | `props` machine-check | $R_{0.30}$ refines $R_{0.15}$ = True; $|R|=678/708/753$ at $\gamma=0.15/0.20/0.30$ |
| **seed-sensitivity** (witness vs invariant) | `refine_families` at seeds $0/1/7$ | blob $\to$ **$48/41/39$** families; total $708/699/699$; only **$22$ of $48$** seed-0 blob families reproduce at seed $1$; **certificate holds at every seed** ($41/41$, $39/39$). $\Rightarrow$ certificate = invariant; family *set/count* = seed-fixed witness |
| **$\sigma$-invariant core** | gate-emitted (density $\ge\gamma$ before $\sigma$) | **614/623** components + all 7 flagship families are seed-independent; only the $\approx9$ sub-$\gamma$ mega-components' partition is seed-dependent |
| **multi-copy** $\Lambda\ge 2$ | `distinct_loci` | all 708 refined families have $\ge 2$ distinct loci = True |
| **coherence: `repeat_mosaic` residual** | `biotype_purity < 0.50` | **9/708** certificate-passing families flagged (mixed-biotype hubs), sorted to end; certificate is necessary-not-sufficient (§3.5); strict oracle+refine still leaves 12 flagged / max 124 |
| verified families survive whole (seed-independent) | `genome_families_refined.tsv` | CEACAM5/6/7 (GRFAM75, 0.470, purity 0.833), ZNF716 (GRFAM54, 0.261, 0.611), IFITM1/2 (GRFAM138, 0.524, 1.000), ULBP1 (GRFAM133, 0.667, 1.000), MAGEA12/CSAG1 (GRFAM91, 0.733, 0.800), PRSS (GRFAM57, 0.316, 0.529), LILRA (GRFAM172, 0.467, 1.000) — all `repeat_mosaic`=0 |
| raw modularity over-splits (why the gate is needed) | Louvain res=1 | CEACAM 12$\to$[8,4], ZNF716 18$\to$[7,7,4] — vetoed by the §3.1 step-1 certificate gate |
| **default (no `--refine`) byte-identical raw** | `md5sum` / `diff` vs committed `genome_families_{annotated,protein_coding,blocks}.tsv` + `genome_family_def.json` | **all 4 IDENTICAL** (refinement is opt-in and additive; the default catalogs + stats json are untouched, refined output goes to its own `genome_families_refined.tsv`) |
| `--refine` leaves the default catalogs untouched | `diff` the 3 default TSVs before/after a `--refine` run | identical (only `genome_families_refined.tsv` is added; the stats `.json` gains `refined_*`/`gamma`/`seed`/`purity`/`locus_overlap` keys) |
| determinism (repeat + threads) | `md5sum genome_families_refined.tsv` across two `--refine` runs and `--threads 12` vs `4` | **byte-identical** (seeded Louvain + KL, integer-node sets, order-independent parallel union) |
| 1,547-gene blob split (seed-0 witness) | map refined families $\subseteq$ raw blob | **48 cohesive families** (sizes 207, 189, 104, 78, 66, 60, 53, 53, 49, 48, 46, 43, …); covers 1,521/1,547 blob genes (26 dropped by $\Lambda<2$); the two largest are the **`repeat_mosaic`-flagged** residual hubs GRFAM699 (207, 23 contigs, $\rho_{\mathrm{in}}=0.265$, purity 0.435) / GRFAM700 (104, 16 contigs, 0.226, purity 0.481) — sorted to the **end** of the catalog, not the top |

**Outputs.** `genome_families_refined.tsv` (**canonical seed-0 witness**, written **only under `--refine`**; columns `family_id, n_genes, n_loci, n_contigs, cross_chrom, density, biotype_purity, repeat_mosaic, dom_biotype, compara_paralog, genes` — `density` is the per-family $\rho_{\mathrm{in}}$ certificate, `biotype_purity`/`repeat_mosaic` are the §3.5 coherence descriptors; rows are sorted genuine-families-first then largest-first, with the 9 `repeat_mosaic` hubs at the end) · `genome_families_annotated.tsv` / `_protein_coding.tsv` (raw, always written, byte-identical with or without `--refine`) · `genome_families_blocks.tsv` (Catalog B) · `genome_family_def.json` (stats; the default run writes a **15-key** baseline, byte-identical to the committed HEAD json; `--refine` adds `refined_families`, `refined_member_genes`, `refined_cross_chrom`, `refined_repeat_mosaic`, `gamma`, `seed`, `locus_overlap`, `purity` for a **23-key** file).

**Definitional one-liner for the thesis body.** A *DNA-level multi-copy gene family* is a **cohesive community** of the genome segdup graph $G_{\mathrm{DNA}}=(V,E_a)$ meeting the **cohesion certificate**: a block $C$ that the density-gated operator $R$ carves from a connected component so that $C$ is a $\gamma$-quasi-clique ($\rho_{\mathrm{in}}\ge\gamma=0.20$) with $\ge 2$ distinct loci. **The certificate — not any single partition — is the canonical object**: exact max-$\gamma$-quasi-clique partition is NP-hard, so $R$ emits one **seed-fixed witness** (deterministic, byte-identical) whose $\gamma$-cohesive core (614 gate-emitted components + the flagship families) is seed-invariant while the blob-derived family count is a witness artifact. It is the **same skeleton, different edge oracle** as the RNA read-conflict O2 substructure ($E_a$ read-independent vs $E_c$ read-derived; the RNA *homology* family proper is $E_r$, the $E_a$ parallel); its refinement is **community-based** (not the RNA side's per-edge significance gate) because the DNA over-merge is a bridge-poor transitive blob; and it supplies the reference copy number $K_{\mathrm{DNA}}=\Lambda(C)\ge\chi(H)=\mathrm{MCC}$ against which RNA copy-assignment is measured. The certificate is necessary-not-sufficient (a `repeat_mosaic` purity descriptor flags the mixed-biotype residual), and the irreducible knobs — COVER $=0.5$, the raw-vs-$\mathrm{SD}(\cdot)$ oracle, the cohesion $\gamma$ (knife-edge $0.061$ above KRAB-ZNF), `PURITY`$=0.5$, and the disclosed splitter seed $=0$ — are named, not hidden.
