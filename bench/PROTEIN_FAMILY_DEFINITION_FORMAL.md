# A Formal Definition of Protein-Level Multi-Copy Gene Families for the Gorilla Pan-Transcriptomic Thesis

**Scope.** This note fixes, for the rest of the thesis, what a **multi-copy gene family at the protein level** is — the field-standard "gene family" notion (OrthoFinder / Pfam / MCL / mmseqs-cluster lineage: *two genes are in one family iff their protein products are homologous*). It is the **fourth edge oracle** $E_p$ in the family lattice, and it is built on the **same combinatorial skeleton and the same cohesion operator $R$** as the three companion objects, differing **only** in the edge oracle. It composes with — and must not contradict — three companion notes:

- `bench/SEGDUP_DEFINITION_FORMAL.md` — the segdup predicate $\mathrm{SD}(\cdot)$ and the genomic edge oracle $E_a$, the exonic-homology oracle $E_b$, the read-ambiguity oracle $E_c$, and the lattice $E_c^{\mathrm{sig}}\subseteq E_c\subseteq E_b$ with $E_a$ incomparable.
- `bench/DNA_FAMILY_DEFINITION_FORMAL.md` — the DNA-level family: the cohesion-refinement operator $R$ ($\gamma$-quasi-clique + $\ge 2$-distinct-loci) applied to a raw $E_a$ component. **This note reuses that $R$ verbatim** (`genome_family_def.refine_families`, imported).
- `bench/family_definition_formal` (in `bench/FAMILY_DEF.md`) — the RNA read-conflict family: component of $E_c$, refined by $E_c^{\mathrm{sig}}\subseteq E_c$, with copy count $\mathrm{MCC}=\chi(H)$.

The **skeleton is shared** across all four (nodes = gene loci, edge = a homology tie, family = refined connected component with $|C|\ge 2$); **the only thing that differs is the edge oracle**. What this note adds is the **protein oracle $E_p$** — the *biological* axis, orthogonal to the genomic ($E_a$), exonic ($E_b$) and read ($E_c$) *nucleotide* axes — the level that captures **ancient / WGD-era paralogs the segdup recency filter structurally misses** (globin witness) and **excludes the non-coding segdup over-merge blob** (DGCR6/DRD5), so the thesis's **operational** family — the read-conflict unit $E_c$, **restricted to its coding-both core** — is a **subset/refinement of the biological family** $E_p$ (a *conditional* containment: $E_c$ and $E_p$ are **incomparable** in general, §4.3).

All empirical statements were re-verified against `/home/juanfra/winloci_scratch/proteins.fa` (22,614 proteins, one per gene, header = gene symbol), `.../aln.m8` (mmseqs all-vs-all, 1,436,582 rows), `.../protclust_cluster.tsv` (mmseqs easy-cluster), `/home/juanfra/winloci_scratch/GGO_genomic.gff` (34,114 RefSeq gene loci) and `/mnt/c/Users/jfris/Desktop/final.bed` (SEDEF), via `bench/protein_family_def.py` on 2026-06-30. Load-bearing counts and machine-checked properties are inlined and reproducible (Appendix A).

> **One-line thesis claim.** A protein-level multi-copy gene family is a **cohesive community** of the protein-homology graph $G_p=(V,E_p)$ satisfying the **same cohesion certificate** as the DNA level: a block $C$ produced by the operator $R$ from a connected component of $G_p$, a $\gamma$-quasi-clique ($\rho_{\mathrm{in}}(C)\ge\gamma=0.20$) with **$\ge 2$ distinct genomic loci** ($\Lambda(C)\ge 2$). The edge oracle is **reciprocal whole-protein homology**: $\{u,v\}\in E_p$ iff the two genes' proteins align in **both** orientations, each hit with **E-value $\le\alpha_p$** and **reciprocal coverage $\min(q_{\mathrm{cov}},t_{\mathrm{cov}})\ge c_p$** — a significance floor plus a whole-protein (not single-domain) coverage floor, **with no identity floor** (the design choice that lets $E_p$ reach the twilight zone). This is the **same skeleton, same operator, different oracle** as $E_a$; $E_p$ is the biological family. $E_c$ and $E_p$ are **incomparable** in general, but every read-conflict edge with **coding-both** endpoints lies in $E_p$ (the *conditional* containment $E_c\cap\{\text{coding-both}\}\subseteq E_p$; both real $E_c$ edges are coding paralogs and confirm it, 2/2), placing the operational copy-assignment **coding core** **inside** it. The irreducible knobs are named, not hidden: $\alpha_p$ (a calibratable Type-I level), $c_p$ (the domain-vs-whole-protein cut, the coding analog of the segdup `COVER`), the scoring scheme $S$ (BLOSUM62 + gap penalties), and the cohesion $\gamma$ — **no threshold-freedom is claimed.**

---

## 1. The canonical object

### 1.1 The protein-homology graph $G_p$

The node set is **shared, with zero mapping.** `proteins.fa` carries exactly **one protein per gene**, and its FASTA header **is** the gene symbol (`>A1BG`, `>A2M`), equal to the GFF `gene=` name that labels the $E_a$ / $E_c$ nodes. So the four oracles live on **one** ground set $V$ = the 34,114 RefSeq gene loci; the protein oracle touches the **22,614** of them that are protein-coding with an annotated product, and **all 22,614 headers map to a gene node (0 loss)** — verified `proteins_mapped = 22614 / 22614`.

For genes $u,v$ with proteins $p_u,p_v$, define the **reciprocal whole-protein homology** oracle. Let the mmseqs all-vs-all report local alignments; a **directional** hit $(q\!\to\!t)$ is *retained* iff

$$
q\neq t,\qquad E(q\!\to\!t)\le \alpha_p,\qquad q_{\mathrm{cov}}\ge c_p\ \wedge\ t_{\mathrm{cov}}\ge c_p,
$$

i.e. it is significant (Karlin–Altschul E-value below the level $\alpha_p$) **and** covers at least a fraction $c_p$ of **both** proteins. Then

$$
\{u,v\}\in E_p \iff (u\!\to\!v)\ \text{retained}\ \wedge\ (v\!\to\!u)\ \text{retained}\qquad(u\neq v),
$$

symmetric and irreflexive ⇒ undirected graph $G_p=(V,E_p)$. **Defaults $\alpha_p=10^{-5}$, $c_p=0.50$** (§2). At the defaults: **1,413,976 non-self hits → 570,565 pass the (E, cov) filter → 19,611 are one-directional-only and DROPPED by reciprocity → 275,477 undirected $E_p$ edges.**

**Two clauses, both load-bearing (verified):**

- **Reciprocity** (both orientations) is not cosmetic: dropping it lets a one-directional *partial* hit (a small protein fully covered by a large multi-domain one, passing only in the short→long direction) glue unrelated families, inflating the largest component. Requiring both orientations makes a chance edge rarer still (§2).
- **Whole-protein reciprocal coverage** ($\ge c_p$ of **both**) demands **coding homology over most of both ORFs**, not one shared domain. This is exactly what **excludes single-domain sharers** (§4.4) and the **non-coding segdup blob** (§4.2) from $E_p$ — the segdup can span two loci genomically with no shared ORF, so those loci carry no $E_p$ edge.
- **No identity (`fident`) floor** — deliberate. Significance + reciprocal coverage alone define homology. This is precisely what lets $E_p$ reach **ancient/divergent** paralogs an identity gate would sever: the globin superfamily is one $E_p$ family down to **24.0% aa identity** (MB↔HBE1, $E=7.8\times10^{-9}$; §4.1), well inside the twilight zone where all three *nucleotide* oracles ($E_a,E_b,E_c$) are saturated.

The **raw protein family** is a connected component of $G_p$ with $|C|\ge 2$ (union-find — the *same* skeleton as `conflict_families()` / the DNA `build`). **3,030 raw families / 15,432 member genes; largest raw component 769.**

### 1.2 The family: same certificate, same operator $R$ as the DNA level

> **Definition (protein-level multi-copy gene family).** Let $\gamma\in(0,1]$. A **protein-level multi-copy gene family** is a vertex set $C\subseteq V$ such that
> $$
> C\in R(\kappa)\ \text{for the connected component }\kappa\text{ of }G_p\text{ containing it,}\qquad\text{and}\qquad\boxed{\ \Lambda(C)\ge 2\ }
> $$
> where $R$ is the **identical** cohesion-refinement operator of `DNA_FAMILY_DEFINITION_FORMAL.md` §3 (imported verbatim from `genome_family_def.refine_families`), so $C$ is a **$\gamma$-quasi-clique** ($\rho_{\mathrm{in}}(C)\ge\gamma$, or $|C|\le 2$), and $\Lambda(C)$ is the **distinct-locus count** (§below). The predicate "$\gamma$-quasi-clique with $\Lambda\ge2$ inside a raw $E_p$ component" is the **cohesion certificate**, and **it — not any one partition — is the canonical object** (exact max-$\gamma$-quasi-clique partition is NP-hard; $R$ returns one seed-fixed *witness*).

$\rho_{\mathrm{in}}(C)=2|E_p(C)|/(|C|(|C|-1))$ is the induced-density on the $E_p$-subgraph. $\Lambda(C)$ collapses reciprocally-overlapping ($\ge 50\%$ both ways) co-located annotations to one physical locus and counts the rest, so $\Lambda(C)\ge 2$ means **$\ge 2$ physical loci** — genuine multi-copy, not one gene wearing two overlapping labels. This is the **same** multi-copy predicate as the DNA and RNA levels, and it is meaningful here **because each gene node carries real genome coordinates** (one protein per gene = one locus).

**Result (defaults, $\gamma=0.20$, seed $=0$, `--refine`).** $3{,}030$ raw $\Rightarrow$ **$3{,}106$ refined families / $15{,}432$ member genes / $2{,}769$ cross-chromosome.** Certificate **$3{,}106/3{,}106$** $\gamma$-quasi-clique-or-size-$\le2$; **all $3{,}106$ have $\Lambda\ge 2$**; **$0$ `repeat_mosaic`** (§4.2). Largest family **769** (`PRFAM0`, the 7-TM GPCR/olfactory superfamily, $\rho_{\mathrm{in}}=0.463\ge\gamma$ ⇒ kept whole). Each family is a **sub-object of one raw $E_p$ component** (property (R1)), exactly the post-refinement reading the lattice demands.

**Top families** (verified in `protein_families_refined.tsv`):

| id | genes | loci | contigs | $\rho_{\mathrm{in}}$ | mmseqs-cluster purity | content |
|---|---:|---:|---:|---:|---:|---|
| PRFAM0 | 769 | 769 | 23 | 0.463 | 0.624 | 7-TM GPCR / olfactory receptor superfamily |
| PRFAM1 | 617 | 616 | 25 | 0.365 | 0.323 | KRAB / C2H2 zinc-finger (ZNF/ZBTB/ZSCAN) |
| PRFAM2 | 171 | 171 | 23 | 0.508 | 0.760 | RAS/RAB/ARF/RHO GTPases (**contains RABL2A/B**) |
| PRFAM3 | 124 | 123 | 21 | 0.411 | 0.629 | serine proteases (PRSS/TMPRSS/MASP) |
| PRFAM5 | 82 | 82 | 11 | 0.772 | 0.939 | keratins / intermediate filaments |
| PRFAM61 | 25 | 25 | 6 | 0.427 | 0.640 | Ig-domain CEACAM/CD superfamily |
| PRFAM135 | 15 | 15 | 5 | 0.829 | 0.933 | **globin superfamily** (§4.1 witness) |
| PRFAM665 | 4 | 4 | 1 | 1.000 | 1.000 | APOBEC3B/D/F/G (**contains the $E_c$ edge**) |

---

## 2. The mmseqs construction and the honestly-named irreducible parameters

### 2.1 The all-vs-all build

$E_p$ is built from **one** mmseqs `easy-search` of `proteins.fa` against itself:

```
mmseqs easy-search proteins.fa proteins.fa aln.m8 tmp \
    --format-output 'query,target,evalue,pident,qcov,tcov' -e 1e-5 -s 7.5 --max-seqs 300
```

(`mmseqs` v18 `18.8cc5c`; all tmp/big outputs under `/home/juanfra/winloci_scratch`, never `/mnt/c`.) `bench/protein_family_def.py` parses each row, keeps directional hits with `evalue≤`$\alpha_p$ ` AND min(qcov,tcov)≥`$c_p$, and emits the reciprocal, irreflexive edge set. `--max-seqs 300` caps hits/query; the 769-gene GPCR superhub exceeds 300, so a few of its *edges* may be truncated, but the **single-linkage component still forms** (component membership is robust to edge truncation), and its certificate density is computed on the retained edges. The mega-components can be made edge-complete with `--max-seqs 1000` if needed.

### 2.2 The irreducible parameters — named, calibratable, not threshold-free

The protein oracle has **three** irreducible knobs plus the scoring scheme. They are the exact analogs of the segdup note's $\alpha_{\mathrm{GW}}$ / `COVER` / scoring scheme, and of the DNA note's $\gamma$.

1. **$\alpha_p$ — the E-value significance level (default $10^{-5}$).** This is a **calibratable Type-I / FDR level**, not an ad-hoc biological constant. mmseqs's E-value is already database-size-corrected, so the expected number of **chance (false-homology) directional hits** across the whole all-vs-all is
   $$\mathbb{E}[\text{chance hits}]\approx \alpha_p\cdot n_{\mathrm{queries}} = 10^{-5}\times 22{,}614 \approx 0.226 < 1,$$
   i.e. **Bonferroni-safe**; requiring **both** orientations makes a chance *edge* far rarer still. It is the protein analog of the segdup Karlin–Altschul $\alpha_{\mathrm{GW}}$ and — by schema, not by number — of the read gate's `min_p`$=\varepsilon^\delta$. **Sensitivity (measured):** loosening $\alpha_p$ by $100\times$ to $10^{-3}$ (expected chance hits $\to 22.6$) changes the object by only $\sim1\%$ — $275{,}477\to278{,}369$ edges ($+1.05\%$), $3{,}106\to3{,}082$ families — so $\alpha_p$ sits in a **stable regime**, not a knife-edge; the choice is disclosed, not free.

2. **$c_p$ — the reciprocal whole-protein coverage floor (default $0.50$).** This is the **one genuinely arbitrary knob**: the domain-vs-whole-protein boundary. It is the coding analog of the segdup **`COVER`$=0.50$** gene-projection fraction. It is what makes $E_p$ *whole-protein* homology (excluding single-shared-domain edges and the non-coding blob, §4). **Sensitivity (measured):** tightening $c_p$ to $0.80$ drops the edge set by $37\%$ ($275{,}477\to173{,}552$) and shrinks the largest component $769\to712$ — so $c_p$ **is** load-bearing, and is declared arbitrary; it is an optional granularity handle (§2.3).

3. **$\gamma$ — the cohesion parameter (default $0.20$, transferred unchanged from the DNA level).** Discussed in §2.3.

**The scoring scheme $S$ (BLOSUM62 + gap penalties + X-drop)** is the irreducible substitution model that fixes the **length↔identity exchange rate** and from which $\lambda,K$ (hence the E-value) derive — the same "irreducible scoring scheme" caveat the segdup note names for its Karlin–Altschul E-value. It is not eliminated; it is *where* the two soft constants (length, identity) are relocated into one significance statement.

**No threshold-freedom is claimed.** $\alpha_p$, $c_p$, $\gamma$ and the scoring model are named parameters whose sensitivity is reported. What is bought is a *reduction*: an identity-and-length cliff becomes one significance level ($\alpha_p$) plus one coverage cut ($c_p$) plus the scoring model — and, unlike an identity floor, this predicate reaches the twilight zone.

### 2.3 $\gamma$ transfers unchanged — and its role flips (honest disclosure)

$\gamma=0.20$ is **transferred verbatim** from the DNA level, **not recalibrated** — because the whole point of the four-level lattice is *one operator + one constant across levels differing only in the edge oracle*; silently retuning $\gamma$ would break "same operator." But its role **flips**:

- On **DNA**, the raw $E_a$ graph is a density-$0.020$ single-linkage repeat blob, so $\gamma=0.20$ is a **load-bearing knife-edge** (just below the lowest real family, KRAB-ZNF at $0.261$).
- On **protein**, reciprocal whole-protein homology yields raw components that are **already near-cliques**: **$3{,}020/3{,}030 = 99.7\%$** of raw components are already $\gamma$-cohesive, so $\gamma=0.20$ is **near-inert** — it splits only the **10** raw components with density $<0.20$, a net **$+76$** families ($3{,}030\to3{,}106$). The median raw component is a clique.

So on the protein level the **real granularity knob is the mmseqs cutoff** ($\alpha_p$, $c_p$) — or the field-standard `--min-seq-id` — **not $\gamma$**. The load-bearing region for $\gamma$ is the **domain-superfamily mega-components** (GPCR $\rho_{\mathrm{in}}=0.463$, KRAB-ZNF $0.365$, GTPase $0.508$), which $\gamma=0.20$ keeps **merged**. These are the protein analog of the $E_a$ repeat-array over-merge: **one shared domain bridging many distinct sub-families**. Raising $\gamma$ toward and past their densities ($\approx 0.37$–$0.51$), or tightening $c_p$ / raising `--min-seq-id`, carves them into **orthogroups** — offered as an **optional sensitivity row**, never baked into the canonical object. (By property (R4) inherited from $R$, the catalog varies monotonically and predictably with $\gamma$.)

---

## 3. Cross-validation against the field-standard clustering

### 3.1 The field standard: mmseqs easy-cluster (OrthoFinder / MCL analog)

The canonical protein-family baseline is a **greedy set-cover / MCL-style clustering**. We run mmseqs's `easy-cluster` (the standard analog):

```
mmseqs easy-cluster proteins.fa protclust tmp --min-seq-id 0.30 -c 0.5 --cov-mode 0 -s 7.5
```

→ `protclust_cluster.tsv` (rep⟨TAB⟩member, one row per gene). `--cov-mode 0` = bidirectional coverage, matching the reciprocal-coverage $E_p$ edge; `--min-seq-id` is the orthogroup-granularity control. Result: **10,154 clusters over 22,614 proteins.** This is an **independent, non-circular** yardstick: it never sees genome coordinates, segdup calls, or reads — only the protein sequences.

### 3.2 Agreement (verified)

Comparing our $\gamma$-community (refined $E_p$) families against the mmseqs-cluster partition, over the 15,432 $E_p$ member genes:

| metric | value | reading |
|---|---:|---|
| **gene-weighted** dominant-cluster purity | **0.860** | the honest headline: 86% of $E_p$ member genes lie in their family's single most-common cluster (size-weighted, not inflated by trivial families) |
| per-family-mean dominant purity | 0.965 | mean over families, each counted once — **inflated** by the 1,562 size-2 families (mean purity 0.995) |
| **large-family** ($\ge30$ genes, $n=44$) mean purity | **0.719** | on the biologically-interesting multi-copy superfamilies the two catalogs agree least — this is the granularity residual (§below) |
| families fully inside one cluster | 2,769 / 3,106 = 89.2% | $E_p$ family ⊆ one mmseqs cluster |
| Adjusted Rand Index | **0.553** | moderate — dominated by the superfamily-granularity difference below |

**Honest reading.** The headline is the **gene-weighted 0.860**, *not* the per-family mean 0.965 (which is inflated by the 1,562 trivial size-2 families at 0.995): the two catalogs agree on **86% of member genes**, and $89.2\%$ of $E_p$ families nest inside one cluster. Agreement is strongest on **recent coding families** — CEACAM/Ig (PRFAM61, purity 0.640), keratins (PRFAM5, 0.939), APOBEC3 (PRFAM665, 1.000), and notably the **globins** (PRFAM135, **0.933** — the field standard *also* groups them; §4.1). It is weakest exactly where expected — the **44 large ($\ge30$-gene) families average only 0.719** — driven by the **~4 domain superfamilies** where $E_p$ keeps the superfamily **merged** while the identity-gated set-cover **splits** it into subfamilies: GPCR (PRFAM0, mmseqs-purity **0.624**), KRAB-ZNF (PRFAM1, **0.323**), GTPase (PRFAM2, **0.760**). This is the same **structured** residual the ARI (0.553) reports — the **same granularity knob** ($\gamma$ / $c_p$ / `--min-seq-id`, §2.3), disclosed as a sensitivity axis, not a disagreement of substance: raising $\gamma$ or `--min-seq-id` (or $c_p$, e.g. $0.80\Rightarrow$ ARI 0.609) moves $E_p$ toward the mmseqs orthogroup granularity and raises the ARI.

---

## 4. The four-level lattice: $E_p$ / $E_a$ / $E_b$ / $E_c$

**Four oracles on one shared node set $V$** (34,114 gene loci; 22,614 with a protein, all mapping with zero mapping). Same skeleton (nodes, edge = homology tie, family = refined component $|C|\ge2$), differ **only** in the edge oracle:

| oracle | edge $\{u,v\}$ holds iff | axis | read-dependence |
|---|---|---|---|
| $E_p$ (protein, **new**) | proteins reciprocally whole-protein homologous ($E\le\alpha_p$, $q_{\mathrm{cov}},t_{\mathrm{cov}}\ge c_p$) | **biological** | read-independent, **coding** |
| $E_a$ (DNA segdup) | a SEDEF segdup pair links the two loci | genomic | read-independent |
| $E_b$ (exonic) | shared segment exonic at $\ge1$ endpoint, genomically homologous at the other | exonic | read-independent |
| $E_c$ (RNA) | a read de-ties the two loci (read-conflict tie) | expressed | **read-derived** |

$E_p$ is the **biological axis, orthogonal to the three nucleotide axes** ($E_a$ genomic, $E_b$ exonic, $E_c$ read). The segdup note's clean chain $E_c^{\mathrm{sig}}\subseteq E_c\subseteq E_b^{\mathrm{asym}}$ (its **only** unconditional containment; $E_a$ incomparable) is preserved **unchanged**; $E_p$ adds **no new unconditional containment** — it is incomparable with $E_a$, $E_b$ and $E_c$ — only one **conditional** result: $E_c$ restricted to coding-both endpoints $\subseteq E_p$ (§4.3).

### 4.1 $E_p$ vs $E_a$ — INCOMPARABLE; $E_p\setminus E_a$ witness = the GLOBINS (ancient paralogs $E_a$ misses)

**Witness (verified).** The seven symbol-annotated globins —
**HBZ, HBM, HBQ1** (α-cluster, tandem on `NC_073242.2`), **HBE1** (β-cluster, `NC_073233.2`), **MB** myoglobin (`NC_086018.1`), **CYGB** cytoglobin (`NC_073228.2`), **NGB** neuroglobin (`NC_073239.2`) — **span 5 chromosomes** yet form **one $E_p$ family** (`PRFAM135`, together with 8 LOC globins, 15 genes, $\rho_{\mathrm{in}}=0.829$). All $\binom{7}{2}=\textbf{21/21}$ globin pairs are co-family (16 direct reciprocal $E_p$ edges + transitive closure).

**Ancient-divergence check (verified from `aln.m8`):** the deep pairs sit squarely in the **twilight zone** —

| pair | aa identity | E-value |
|---|---:|---:|
| HBZ–HBE1 | 39.5% | $5.6\times10^{-29}$ |
| MB–CYGB | 29.4% | $4.9\times10^{-17}$ |
| HBE1–NGB | 28.2% | $2.0\times10^{-6}$ |
| MB–HBZ | 27.0% | $2.6\times10^{-14}$ |
| HBM–CYGB | 26.1% | $6.2\times10^{-11}$ |
| MB–HBE1 | **24.0%** | $7.8\times10^{-9}$ |

(E-values are the **canonical minimum-E-value row per unordered pair** in `aln.m8`, reproducible via the Appendix-A globin extractor; 16 of the 21 pairs have a direct reciprocal $E_p$ edge, the other 5 are co-family by transitive closure.)

**None of the 21 globin pairs is an $E_a$ segdup edge ($\textbf{0/21}$) and none is a Bailey $\ge90\%$ segdup ($\textbf{0/21}$)** — verified by the same one-pass `final.bed` overlay that emits both edge sets. These paralogs (α/β split ~450 Mya; myo/cyto/neuro deeper) have diverged **far below** the segdup identity floor and are dispersed across five chromosomes, so they are exactly the **ancient / WGD-era paralogy the segdup recency filter is designed against** (`SEGDUP_DEFINITION_FORMAL.md` §1.1: "excludes ancient WGD ohnologs"). **$E_p\setminus E_a$ is non-empty, cross-chromosomal, and biologically real** — and the field-standard mmseqs clustering **corroborates** it (PRFAM135 mmseqs-purity 0.933). This is $E_p$'s distinct biological value: **protein homology reaches past the twilight zone that all three nucleotide oracles are saturated under.** (Note: HBA1/HBA2/HBB/HBD/HBG are not separately symbol-annotated in this gorilla T2T GFF — LOC-named or collapsed — but the seven present symbols already establish the claim.)

$E_a\setminus E_p$ is witnessed by the blob (§4.2). So $E_p$ and $E_a$ are **incomparable**.

### 4.2 $E_a\setminus E_p$ — the 22q11 non-coding-segdup over-merge blob does NOT form in $E_p$

**Witness (verified).** In $E_a$ the pericentromeric 22q11 region is a **1,547-gene over-merge blob** (the largest raw $E_a$ component; `SEGDUP`/`DNA` notes), in which six probe genes **DGCR6, DRD5, SLC6A8, PRODH, OTOP1, BCAP31** sit in **one** family — $\binom{6}{2}=\textbf{15/15}$ pairwise co-family in $E_a$. In $E_p$ they land in **SIX DISTINCT protein families** — $\textbf{0/15}$ pairwise co-family:

| gene | $E_p$ family | note |
|---|---|---|
| DGCR6 | PRFAM1057 | |
| DRD5 | **PRFAM0** | dopamine receptor is a 7-TM GPCR |
| SLC6A8 | PRFAM70 | SLC6 transporter |
| PRODH | PRFAM1375 | proline dehydrogenase |
| OTOP1 | PRFAM588 | otopetrin |
| BCAP31 | PRFAM1658 | |

**Mechanism.** The reciprocal **whole-protein coverage** clause ($q_{\mathrm{cov}},t_{\mathrm{cov}}\ge c_p$) requires **coding homology over most of both proteins**, which these genes do not share; the segdup spans the loci **genomically** through non-coding/pericentromeric homology **without any shared ORF**. So $E_p$ **excludes the blob by construction** — an **independent, non-circular** confirmation that the blob is a **segdup artifact, not a gene family** (the protein oracle never sees segdup coordinates). This is why $E_p$ is **cleaner than the raw $E_a$ superset**: it excludes the blob at the *oracle* level, before any community cut, and consequently has **0 `repeat_mosaic` families** vs the DNA level's 9 (the coding oracle is biotype-pure by construction). (Second instance: MAGEA1/MAGEA12 and CSAG1 are co-located by a cancer-testis segdup in $E_a$ but are **distinct** protein families in $E_p$.)

### 4.3 $E_c$ vs $E_p$ — INCOMPARABLE in general; a CONDITIONAL containment on the coding-both core

$E_c$ (read-conflict) and $E_p$ are **incomparable** — the same status as $E_a$ vs $E_p$ — and this is *forced* by the companion segdup note, which must not be contradicted:

- **$E_c\setminus E_p\neq\varnothing$.** The read-conflict predicate (`read_conflict.rs`: de-tie on per-placement `de` + `aln_len`, with **no coding/biotype gate**) lets a read from an expressed locus $u$ multimap onto an **intronic / UTR / unexpressed-paralog / retrocopy-pseudogene / lncRNA** stretch of $v$ — expressed-sequence-identical at $u$, but with **no protein at $v$** (or no whole-protein reciprocal coverage). Such an edge is in $E_c$ but **not** in $E_p$. This is exactly the asymmetry for which the segdup note makes $E_b$ deliberately one-sided ($E_c\subseteq E_b^{\mathrm{asym}}$ is that note's *only* unconditional containment).
- **$E_p\setminus E_c\neq\varnothing$.** Divergent-but-not-read-confusable paralogs (the globins: MB vs HBB never share a read) are in $E_p$ but not $E_c$.

So **$E_c\not\subseteq E_p$ and $E_p\not\subseteq E_c$** (a non-empty $E_c\setminus E_p$ *refutes* containment — it does not license "$E_c\subseteq E_p$, not equality"; that would require $E_c\setminus E_p=\varnothing$). What *is* true is a **conditional containment** on the coding-both core:

$$\boxed{\ E_c\ \cap\ \{\text{both endpoints protein-coding, shared stretch coding at both}\}\ \subseteq\ E_p\ }$$

i.e. **the operational family, restricted to its coding-both core, is a subset of the biological family.** Mechanism: a coding-both read-tie means the two expressed sequences are near-identical over a ~2.5 kb read $\Rightarrow$ their proteins are reciprocally whole-protein homologous $\Rightarrow\{u,v\}\in E_p$. It is a **strict** subset even there — $E_c$ keeps only the read-**confusable** members (the $\chi(H)$ / K-frontier core), not the whole evolutionary family: within globins only near-identical pairs ever read-conflict; MB vs HBB never do.

**Empirical support (weak by construction).** The only real, gene-symbol-labelled $E_c$ edges on this substrate are **2** with $n_{\mathrm{conflict}}>0$ — **APOBEC3D↔APOBEC3F** → **PRFAM665** and **RABL2A↔RABL2B** → **PRFAM2**, both direct $E_p$ edges and co-family (`ec_in_ep_family = 2/2`, `ec_is_ep_edge = 2/2`). But **both are annotated protein-coding paralogs**, a sample that *structurally cannot* contain the $E_c\setminus E_p$ counterexample class, so they confirm **only the coding-both conditional**, not any unconditional $E_c\subseteq E_p$. The counterexample class is *guaranteed* by the segdup note's own asymmetric-$E_b$ construction, but is **not yet a byte-verified row** (no larger labelled gene-pair $E_c$ edge list exists on disk). (APOBEC3D/F is *also* the $E_c\setminus E_a$ witness of the segdup note — their only covering SEDEF pair is 88.4% < 90% $\Rightarrow\notin E_a$ — yet they are one protein family, so the coding-both conditional holds precisely where $E_c\subseteq E_a$ fails.)

### 4.4 $E_p$ vs $E_b$ — INCOMPARABLE

- **$E_p\setminus E_b$:** the ancient globins. MB–HBE1 (24% aa identity) is protein-significant ($E=7.8\times10^{-9}$) but **DNA-saturated** — not nucleotide-alignable within $E_b$'s divergence bound `de_max`, so $\notin E_b$. Protein homology reaches past the twilight zone the three nucleotide oracles live under.
- **$E_b\setminus E_p$:** single-shared-exon / domain-sharers (below $c_p=0.5$ of the protein) — e.g. ASDURF/ASNSD1 — have partial exonic homology but no **whole-protein** reciprocal coverage ⇒ $\notin E_p$.

### 4.5 Lattice summary and edge-level overlaps

**Precise structure.** The segdup note's clean chain is preserved **unchanged**, and its **only unconditional containment stays the only one** — $E_c\subseteq E_b^{\mathrm{asym}}$. $E_p$ is a new circle **incomparable with all three** existing oracles:

$$
E_c^{\mathrm{sig}}\subseteq E_c\subseteq E_b^{\mathrm{asym}}\qquad(\text{the ONLY unconditional containment in the lattice}),
$$
$$
E_p\ \text{incomparable with}\ E_a,\ E_b,\ E_c;\qquad E_a\ \text{incomparable with}\ E_b,\ E_c,\ E_p.
$$

The one positive $E_p$ result is **conditional**: $E_c\cap\{\text{coding-both}\}\subseteq E_p$ (§4.3). Non-containments: 4 from the segdup note $\{E_a\text{ vs }E_b,\ E_a\text{ vs }E_c\}$ + **6 new** $\{E_p\not\subseteq E_a,\ E_a\not\subseteq E_p,\ E_p\not\subseteq E_b,\ E_b\not\subseteq E_p,\ E_p\not\subseteq E_c,\ E_c\not\subseteq E_p\}$, all witnessed above.
- **Agreement region $E_a\cap E_p$** = recent **coding** segdup families (CEACAM/Ig, PRFAM61).
- **Quad-core $E_a\cap E_b\cap E_c\cap E_p$** = **RABL2A/B** (PRFAM2) — all four oracles see the same two copies.

**Edge-level overlaps** (`protein_family_def.json`, on the shared node set). **Symbol note:** here $E_a$ = the **raw $\ge 50\%$-floor `final.bed` superset** ($|E_a|=46{,}206$ gene edges), *distinct* from `SEGDUP_DEFINITION_FORMAL.md`'s **canonical** $E_a=\mathrm{SD}(\cdot)\ge 90\%$ Bailey predicate (26,892 gene edges here); the witnesses additionally check the Bailey $\ge90\%$ set (globins: 0/21 under **both**).

| quantity | value |
|---|---:|
| $\|E_p\|$ | 275,477 |
| $\|E_a\|$ (segdup gene edges) | 46,206 |
| $\|E_p \setminus E_a\|$ | 264,654 |
| $\|E_a \setminus E_p\|$ | 35,383 |
| $\|E_p \cap E_a\|$ | 10,823 |

**Check:** $264{,}654 + 10{,}823 = 275{,}477 = \|E_p\|$. ✓ The **$E_a$-only 35,383 decomposes cleanly**: **25,242** edges have a **non-protein-coding endpoint** (structurally impossible in $E_p$: $46{,}206 - 20{,}964$ protein-coding-both) + **10,141** protein-coding pairs joined by **non-whole-protein** segdup homology (intronic / partial-domain: $20{,}964 - 10{,}823$). On the **20,964 protein-coding-both $E_a$ edges, 10,823 (51.6%) are also $E_p$ edges.**

### 4.6 The operational family's coding core is a subset/refinement of the biological family

Three precise levels (all computed in `bench/protein_family_def.py`):

1. **Read-conflict (edge level):** the *conditional* containment $E_c\cap\{\text{coding-both}\}\subseteq E_p$ (§4.3); the **2** real $E_c$ edges are both coding paralogs and confirm it (`ec_in_ep_family = 2/2`), so the copy-assignment **coding core** is contained in a genuine protein family. This is **not** an unconditional $E_c\subseteq E_p$ — non-coding-partner read-conflicts are in $E_c\setminus E_p$, so $E_c$ and $E_p$ are incomparable.
2. **DNA segdup family (family level):** **470 / 607 = 77.4%** of canonical refined $E_a$ families (GRFAM) with a protein member sit **fully inside one $E_p$ family**. The **137** that span $>1$ $E_p$ family are segdups **co-locating distinct protein families** (MAGEA+CSAG; the residual repeat-mosaic hubs).
3. **Coding segdup (edge level):** **10,823 / 20,964 = 51.6%** of protein-coding-both $E_a$ edges are $E_p$ edges (§4.5).

So the statement is **not** $E_a\subseteq E_p$ ($E_a$ is incomparable — the blob witness) **nor** an unconditional $E_c\subseteq E_p$ ($E_c$ is incomparable too — non-coding read-conflicts). It is precisely: **"the coding-both core of $E_c$ lies in $E_p$, and the coding, expressed part of $E_a$ lies in $E_p$."** Consequently $E_p$ acts as an **independent refiner of $E_a$** (splitting co-located distinct protein families — it cures the blob the DNA community-cut has to work to separate), while $E_a$ refines $E_p$ in the **recent-vs-ancient** direction (segdups are the recent, high-identity slice of the protein family). The thesis's operational family's **coding core** is a **subset/refinement of the biological family $E_p$** (the conditional above; $E_c$, $E_a$ each incomparable with $E_p$ in general).

---

## 5. Results and Limitations

### 5.1 Results

- **Canonical object:** protein-level multi-copy gene family = $\gamma$-quasi-clique community (density $\ge\gamma=0.20$, $\ge2$ distinct loci) of the reciprocal-whole-protein-homology graph $G_p=(V,E_p)$ — **same skeleton, same operator $R$, different (coding) oracle** as the DNA family. Certificate **3,106/3,106**; **0 `repeat_mosaic`**; deterministic downstream of the fixed `aln.m8` (seed 0; the Python catalog is byte-identical across reruns / thread counts / `PYTHONHASHSEED`, md5 `b3b3b9e3…`).
- **Reproduces the field standard (honestly weighted):** vs mmseqs easy-cluster (10,154 clusters, OrthoFinder/MCL analog), **gene-weighted** dominant purity **0.860** (per-family mean 0.965, inflated by the 1,562 size-2 families; large $\ge30$-gene families, $n=44$, mean **0.719**), **89.2%** of $E_p$ families nested in one cluster, ARI **0.553** — the residual is the disclosed superfamily-granularity axis.
- **Captures ancient paralogs $E_a$ misses:** the globin superfamily (PRFAM135), 5 chromosomes, down to 24% aa identity, **0/21 segdup pairs, 0/21 Bailey$\ge$90% pairs** — $E_p\setminus E_a$ non-empty and biologically real.
- **Excludes the segdup over-merge blob by construction:** the 22q11 six-gene probe → **6 distinct $E_p$ families (0/15 co-family)** vs **15/15 in $E_a$** — an independent, non-circular confirmation that the blob is a segdup artifact; hence $E_p$ is **cleaner than the raw $E_a$ superset** (0 vs 9 `repeat_mosaic`).
- **Positions the operational family:** the *conditional* containment $E_c\cap\{\text{coding-both}\}\subseteq E_p$ (both real $E_c$ edges confirm it, 2/2), and 77.4% of coding DNA families sit inside one $E_p$ family — the operational copy-assignment **coding core** is a **subset/refinement of the biological protein family**. Lattice: $E_c^{\mathrm{sig}}\subseteq E_c\subseteq E_b^{\mathrm{asym}}$ is the **only** unconditional containment; $E_p$ and $E_a$ are incomparable with all other oracles.

### 5.2 Limitations (honest)

1. **Irreducible parameters, not threshold-freedom.** $\alpha_p$ (calibratable Type-I level, Bonferroni-safe at $10^{-5}$), $c_p$ (the one genuinely arbitrary domain-vs-whole-protein cut, = the segdup `COVER` analog), the scoring scheme $S$ (BLOSUM62 + gaps — where the length↔identity trade-off is relocated, not removed), and $\gamma$ are all named knobs with reported sensitivity. **No cut vanishes.**
2. **$\gamma$ is near-inert here (role flip).** 99.7% of raw $E_p$ components are already $\gamma$-cohesive; the real granularity control on the protein level is the mmseqs cutoff ($\alpha_p,c_p$) / `--min-seq-id`, not $\gamma$. The domain-superfamily mega-components (GPCR/ZNF/GTPase) are kept **merged** at $\gamma=0.20$; whether they are "one family" or should be carved into orthogroups is a $\gamma$/`--min-seq-id`-relative statement, disclosed as a sensitivity row, not a monophyly claim.
3. **Superfamily merging drives the moderate ARI (0.553).** $E_p$ (no identity floor) is deliberately **more inclusive** than the identity-gated field standard on domain superfamilies — a feature for reaching ancient paralogs, but it means the two catalogs disagree on GPCR/ZNF granularity. This is disclosed, not error.
4. **Canonical relative to a fixed `aln.m8`; `--max-seqs 300` truncation.** Determinism is certified **downstream** of the archived `aln.m8`: the Python catalog is byte-identical across reruns / thread counts / `PYTHONHASHSEED` (md5 `b3b3b9e3…`), but the canonical $E_p$ object is defined **relative to that one archived search**, not diffed against a fresh mmseqs rebuild. Under `--max-seqs 300` a few edges of the 769-gene GPCR superhub are capped, so component **membership** is robust (single-linkage) but the certificate **density value** on the mega-components (e.g. GPCR $\rho_{\mathrm{in}}=0.463$) is a **truncation-dependent lower bound** — it decides keep-whole-vs-split, and here sits far from $\gamma=0.20$ (so the 10 splits are robust), but is not fully canonical until `--max-seqs` is raised (e.g. 1000).
5. **Coordinates are haploid-T2T.** $\Lambda(C)\ge2$ counts distinct **loci in the single haploid assembly**; a collapsed paralog mis-assembled as one locus is under-counted (the same residual assembly caveat as the segdup note). The protein oracle validates *homology*, not per-copy read-depth copy number.
6. **$E_c$ is thin (2 edges) on this substrate**, and **both are coding-both paralogs**, so they can confirm **only the coding-both conditional** $E_c\cap\{\text{coding-both}\}\subseteq E_p$ — **not** any unconditional $E_c\subseteq E_p$, which is **false**: non-coding / UTR / retrocopy read-conflicts lie in $E_c\setminus E_p$, so $E_c$ and $E_p$ are incomparable (§4.3). The conditional itself is *structural* (a coding-both read-tie requires near-identical expressed sequence ⇒ whole-protein homology). No larger labelled gene-pair $E_c$ edge list exists on disk to exhibit the $E_c\setminus E_p$ counterexample class as a byte-verified row.

---

## Appendix A — verification & reproduction log (2026-06-30)

Reproduce: `/home/juanfra/miniforge3/bin/python bench/protein_family_def.py --refine`. Inputs in `/home/juanfra/winloci_scratch`: `proteins.fa`, `aln.m8` (mmseqs `easy-search`, cols `query,target,evalue,pident,qcov,tcov`), `protclust_cluster.tsv` (mmseqs `easy-cluster`). Outputs (in `bench/`): `protein_families_refined.tsv` (canonical, `--refine` only), `protein_families.tsv` (raw baseline, always), `protein_family_def.json`.

| check | source | result |
|---|---|---|
| node set / mapping | `load_genes` + `proteins.fa` | 34,114 gene nodes; 22,614 proteins; **22,614/22,614 map** (0 loss) |
| $E_p$ build | `aln.m8`, $\alpha_p=10^{-5}$, $c_p=0.50$ | 1,413,976 non-self hits → 570,565 pass → 19,611 one-way dropped → **275,477 edges** |
| chance-hit calibration | $\alpha_p\cdot n_{\mathrm{queries}}$ | $10^{-5}\times22{,}614=\textbf{0.226}<1$ (Bonferroni-safe) |
| raw families | union-find on $E_p$ | **3,030 families / 15,432 genes; max 769**; 3,020/3,030 (99.7%) already $\gamma$-cohesive |
| refined ($\gamma=0.20$, seed 0) | `refine_families` (imported) | **3,030 raw → 3,106 refined / 15,432 genes / 2,769 cross-chrom** |
| certificate | in-run density check | **3,106/3,106** $\gamma$-quasi-clique-or-$\le2$ = True |
| multi-copy $\Lambda\ge2$ | `distinct_loci` | all 3,106 families $\ge2$ distinct loci = True |
| coherence residual | `biotype_purity < 0.50` | **0/3,106 `repeat_mosaic`** (coding oracle ⇒ biotype-pure) |
| determinism (downstream of a **fixed** `aln.m8`) | rerun + `PYTHONHASHSEED` | `protein_families_refined.tsv` byte-identical (md5 `b3b3b9e3…`); the canonical $E_p$ is defined **relative to the archived `aln.m8`** — the Python catalog is what is certified byte-identical, not a fresh mmseqs rebuild |
| $\alpha_p$ / $c_p$ sensitivity | rerun `--evalue` / `--cov` | $\alpha_p\,10^{-5}\!\to\!10^{-3}$: 275,477→**278,369** edges (+1.05%), 3,106→**3,082** fams (stable); $c_p\,0.50\!\to\!0.80$: 275,477→**173,552** edges (−37%), max comp 769→**712** (load-bearing) |
| xval vs mmseqs cluster | `protclust_cluster.tsv` | **10,154 clusters**; dominant purity **gene-weighted 0.860** / per-family-mean 0.965 / large($\ge30$,$n=44$) **0.719**; **2,769/3,106 (89.2%)** nested; **ARI 0.553** |
| superfamily granularity | per-family mmseqs purity | GPCR PRFAM0 **0.624**, KRAB-ZNF PRFAM1 **0.323**, GTPase PRFAM2 **0.760** ($E_p$ merges, mmseqs splits) |
| **WITNESS 1 — globins $E_p\setminus E_a$** | `protein_families_refined.tsv`, `final.bed` overlay | HBZ/HBM/HBQ1/HBE1/MB/CYGB/NGB → **PRFAM135**, 5 chroms; **21/21 co-family**; **0/21 $E_a$**, **0/21 Bailey$\ge$90%**; MB–HBE1 24.0% aa id $E=7.8\text{e-}9$ |
| **WITNESS 2 — blob exclusion $E_a\setminus E_p$** | `protein_families_refined.tsv`, $E_a$ overlay | DGCR6/DRD5/SLC6A8/PRODH/OTOP1/BCAP31 → **6 distinct $E_p$ families**; **0/15** co-family in $E_p$ vs **15/15** in $E_a$ |
| $E_c\cap\{\text{coding-both}\}\subseteq E_p$ (conditional) | `family_def_readconflict.tsv` | **2/2** read-conflict edges (both coding paralogs) co-$E_p$-family AND direct $E_p$ edges (APOBEC3D/F→PRFAM665; RABL2A/B→PRFAM2) — confirms the coding-both conditional **only**; $E_c$, $E_p$ incomparable in general |
| edge-level lattice | `protein_family_def.json` | $\|E_p\setminus E_a\|$=264,654, $\|E_a\setminus E_p\|$=35,383, $\|E_p\cap E_a\|$=10,823; 264,654+10,823=275,477 ✓; 51.6% of PC-both $E_a$ edges are $E_p$ |
| family-level operational-subset | GRFAM ⊆ $E_p$ | **470/607 = 77.4%** refined DNA families inside one $E_p$ family; 137 span $>1$ |

**Definitional one-liner for the thesis body.** A *protein-level multi-copy gene family* is a **cohesive community** of the reciprocal-whole-protein-homology graph $G_p=(V,E_p)$ meeting the **same certificate** as the DNA family: a $\gamma$-quasi-clique ($\rho_{\mathrm{in}}\ge\gamma=0.20$) with $\ge2$ distinct loci, carved by the **same operator $R$** from a raw $E_p$ component. Its edge is **significant reciprocal whole-protein homology** ($E\le\alpha_p=10^{-5}$, $q_{\mathrm{cov}},t_{\mathrm{cov}}\ge c_p=0.50$, **no identity floor**), which reaches the twilight zone: it captures the ancient globin paralogs the segdup recency filter misses ($E_p\setminus E_a$, 0/21 segdup pairs) and excludes the non-coding 22q11 segdup blob by construction ($E_a\setminus E_p$, 6 distinct families). It is the **biological** fourth oracle — orthogonal to the genomic ($E_a$), exonic ($E_b$) and read ($E_c$) nucleotide axes — reproducing the field-standard clustering to **gene-weighted 0.860** dominant purity / 89.2% nested (ARI 0.553; per-family mean 0.965 inflated by size-2 families, large families 0.719; superfamily-granularity residual disclosed), and the thesis's operational family's **coding core** is a **subset/refinement** of it (the *conditional* $E_c\cap\{\text{coding-both}\}\subseteq E_p$, both real $E_c$ edges confirm; 77.4% of coding DNA families inside one $E_p$ family; $E_c$ and $E_p$ incomparable in general). The irreducible knobs — $\alpha_p$, $c_p$, the scoring scheme, and $\gamma$ — are named, not hidden; **no threshold-freedom is claimed.**
