# A Formal Definition of Segmental Duplications for the Gorilla Pan-Transcriptomic Multi-Copy-Family Thesis

**Scope.** This note fixes, for the rest of the thesis, what a *segmental duplication* (SD) is, in three registers that must agree: (1) the **operational** field-standard definition (Bailey 2002 / Eichler lineage / SEDEF), including exactly what our `final.bed` encodes; (2) a **relational reframe** — *interpretive, not an operationalized predicate* — that **reinterprets** Bailey's two soft biological constants (1 kb, 90%) as one significance tail against an explicit null, while **naming the irreducible knobs honestly** (the scoring scheme, α, and the repeat-exclusion cut *R*); and (3) a **graph formalization** in which families are connected components, together with the precise — and deliberately *non-nested* — relationship between **genomic-segdup homology** ($E_a$), **exonic/transcript homology** ($E_b$), and **read-ambiguity** ($E_c$). Section (4) wires this back to our copy-assignment framework (`genome_family_def`, the RNA read-conflict graph, MCC = χ(H), the K-frontier, and the SEDEF-identity-vs-PSV break-test). **The formal lattice of §3 is proved on the Bailey-thresholded predicate `SD(·)`; §2's E-value is a reframe, computed nowhere in the pipeline.**

All empirical statements were re-verified against `/mnt/c/Users/jfris/Desktop/final.bed` and `bench/genome_family_def.py` on 2026-06-30; verification commands and counts are inlined where they are load-bearing.

> **One-line thesis claim.** A segmental duplication is a pair of genomic intervals whose local self-alignment is *significant above the genome's unique-alignment background*, excluding high-copy interspersed repeats — operationally a SEDEF pair with ℓ ≥ 1 kb, ρ ≥ 90%, ¬TE. This is **axis (a)** only. It is **not** the same object as our RNA read-conflict family, and it does **not** contain it: the only clean containment in the whole picture is *read-ambiguity ⟹ (asymmetric) exonic homology* ($E_c \subseteq E_b$), with genomic segdup ($E_a$) an **incomparable third circle** — $E_c \setminus E_a$ witnessed by APOBEC3D/F (∈ $E_c$, genomic SD only 88.4% < 90% ⇒ ∉ $E_a$), $E_a \setminus E_c$ by the unexpressed-SD mass. RABL2A/B is in all three (a verified triple-core SD, **not** a retrocopy).

---

## 1. Operational definition (Bailey / Eichler / SEDEF)

### 1.1 The field-standard predicate

Following Bailey et al. 2002 (*Science* 297:1003; verbatim-checked against `/mnt/c/Users/jfris/Desktop/bailey_famcn.md`), a **segmental duplication** is a genomic segment satisfying **all** of:

1. **length ≥ 1 kb**,
2. **sequence identity ≥ 90%** between copies,
3. present in **≥ 2 copies** in the (single-haplotype) genome,
4. **not** a high-copy interspersed repeat / transposable element (TEs are masked before detection and re-filtered after).

What each criterion is *for* — this matters because each one is exactly the boundary against one of the look-alikes the thesis must exclude:

- **≥ 1 kb** excludes short, exon-level, or spurious homology. A single shared exon or protein domain is sub-kilobase; an SD must be a *contiguous genomic block*, not merely shared coding sequence. **This is precisely what separates genomic SD homology (axis a) from mere exonic/transcript homology (axis b).**
- **≥ 90% identity** is a **recency filter**. Bailey frames SDs as recent (≲ 40 My); at the neutral primate substitution rate, 90% identity ≈ the last ~40 My ("On the basis of neutral expectation of divergence, this corresponds to duplications that have emerged over the last ~40 My"). This is the criterion that **excludes ancient whole-genome-duplication ohnologs** (the **2R** vertebrate WGD, ~450–500 Mya; the teleost-specific 3R never occurred in the primate lineage and is *not* cited here), which have diverged far below 90%.
- **≥ 2 copies** — duplication, not unique sequence.
- **TE exclusion** — high-copy interspersed repeats (Alu/LINE/SINE/HERV/DNA transposons) are genome-wide multi-copy and often ≥ 90% identical, but are *not* lineage-specific SDs. They are RepeatMasker-masked before detection and re-filtered after ("due to the presence of high-copy-number repeats, which were then filtered").

### 1.2 Detection: the WGAC ∩ WSSD intersection

The core 2002 contribution is that a single assembly self-alignment cannot, by itself, tell a true duplication from an unmerged **allelic** overlap ("the inability to distinguish between allelic and duplication sequence overlap"). The fix is the intersection of two orthogonal methods:

- **WGAC** (Whole-Genome Assembly Comparison): BLAST-based self-alignment of the assembly, keeping pairs ≥ 1 kb and ≥ 90% identity. WGAC *alone* flagged 16.5% of the genome as duplicated, but four of five alignments > 98% identity were **allelic overlaps / mis-assembly artifacts**, not paralogy.
- **WSSD** (Whole-Genome Shotgun Sequence Detection): WGS read depth in 5-kb windows; depth over-representation (> 3 SD above the unique-region mean) marks duplication. Read depth tracks **known copy number at R² = 0.96** — i.e. WSSD measures copy number directly from depth.
- **SD = WGAC ∩ WSSD.** Filtering WGAC with WSSD removed ~85% of the > 98%-identity alignments (the allelic artifacts) and dropped "duplicated" genome from 16.5% to **5.2%** (130.5 Mb).

A sequence difference between two SD copies is a **paralogous sequence variant (PSV)** — categorically distinct from an allelic SNP. Bailey estimated ~100,000 PSVs were contaminating dbSNP ("roughly one of two SNPs is, in fact, a paralogous sequence variant rather than an allele"). **The PSV/SNP distinction is the conceptual hinge on which the entire copy-assignment program turns.**

**Enrichment** (Bailey Table 2 — the families this thesis targets): Immunoglobulin/MHC (4.0×; class I 12.0×), somatotropin (31.8×), trypsin serine proteases PRSS (4.3×), KRAB-box ZNF (3.5×), cytochrome P450 (5.5×), galectins (7.8×), defensins (22.5×). These are the same families our genome family graph recovers (CEACAM, MAGE-A, KRAB-ZNF, PRSS1/2/3).

### 1.3 SEDEF / BISER = modern WGAC, and *only* WGAC

SEDEF (Numanagić et al. 2018) and its faster successor BISER **are** modern WGAC: a three-stage *seed → chain → align* pipeline with an explicit error/statistics filter, run on a single assembled genome. **It operationalizes axis (a) — genomic self-homology — and only (a).** It is read-independent: it never sees an RNA read and never sees WGS depth.

Our exact run (`/home/juanfra/winloci_scratch/sedef_build/sedef.sh`, default `stat_params=""`): `sedef search -k 12 -w 16` (minimizer seeds, both strands) → bucketing + `sedef align generate -k 11` → `sedef stats generate` → `final.bed`, then `sort | uniq`. **No `awk` post-filter is applied.** SEDEF's default error model (verified from the binary strings):

| parameter | default | role |
|---|---|---|
| maximum divergence rate between SDs | 0.30 | chain/seed model floor (~70% model identity) |
| maximum small mutation rate | 0.15 | edit-error component |
| gap frequency rate | 0.005 | |
| minimal allowed SD length after gap-splitting | 1000 | the 1 kb floor |
| **maximum WGAC-scaled error rate allowed for SD** | **0.5** | **the reporting floor — identity down to 50%** |

**Load-bearing consequence (must be stated honestly in the thesis).** `final.bed` is a **SEDEF candidate superset at a ~50%-identity / ~1 kb floor, not the textbook ≥ 90% SD set.** To recover the Bailey/Eichler operational SD you must post-filter (Check F below). On a *finished haploid* T2T this WGAC-alone object is nevertheless far closer to a clean SD oracle than WGAC-alone was in 2002 — see §1.5.

### 1.4 What `final.bed` encodes (verified)

`final.bed` = 253,029 data rows + 1 trailing `#`-prefixed header line (a SEDEF concat artifact — the header is the **last** line; a naïve `head -1` parser sees a data row, so any parser must skip `#`-prefixed lines). 34 tab-separated columns; the columns this thesis uses:

| col (1-based) | field | meaning |
|---|---|---|
| 1–3 | chrom1,start1,end1 | region A (0-based half-open) |
| 4–6 | chrom2,start2,end2 | region B |
| 9–10 | strand1,strand2 | col 10 = `-` ⇒ inverted/reverse-complement pair |
| 11 | max_len | length of the longer side |
| **12** | **aln_len** | **alignment length** |
| 17 / 18 | matchB / mismatchB | matched / mismatched aligned columns |
| 19–20 | transitions / transversions | (sum = mismatchB) |
| **21** | **fracMatch** | **fraction identity = matchB / alnB, gap-excluded** |
| 22 | fracMatchIndel | identity counting indels in the denominator |
| **23** | **jck** | **divergence = Jukes-Cantor distance** |
| 24 | k2K | Kimura 2-parameter distance |
| 26–28 | uppercaseA/B/uppercaseMatches | **non-RepeatMasked content — the TE/satellite handle** |
| 29 | aln_matches | |

Encoding cross-checks **pass** on row 0: fracMatch 0.948223 = 5256/5543; jck 0.053651 = −¾·ln(1 − 4⁄3·0.051777); transitions+transversions (259+28) = mismatch (287); matchB+mismatchB (5543) = alnB; max_len (5545) = endA − startA. So **cols 12 / 21 / 23 = aln_len / fracMatch / Jukes-Cantor divergence**, confirmed.

**Empirical distributions** (re-verified 2026-06-30, n = 253,029):

- **Identity (col 21):** min **0.5071**, p1 0.692, p5 0.770, **median 0.8696**, max 1.0. Rows < 0.90 = **183,226 (72.41%)**; the floor is SEDEF's `--max-error 0.5` default, **not** 90%. (Bailey's own 2002 set sat at ~91–100%; this raw file is a candidate superset.)
- **Length (col 12):** min **900**, p1 940, **median 4,301**, p99 123,849, max 3,261,496; < 1 kb = **8,671 (3.43%)** (soft edge from gap-splitting). Bins: [1k,5k) 52.3%, [5k,10k) 30.4%, [10k,50k) 10.9%, ≥ 50k 3.0%.
- **Divergence (col 23):** min ≈ 0, median 0.1433, p95 0.275, p99 0.397, max 0.803. (Identity and JC are the same axis reparameterized: the 0.803-JC tail *is* the 0.507-identity tail.)
- **Geography:** **84.28% interchromosomal**, 15.72% intrachromosomal; **50.75% inverted** (col 10 = `-`); **1,939 exact self-pairs** (regionA == regionB) that a downstream consumer should drop.

**Canonical-filter survival** (the operational SD recovery — Check F):

```
raw                                   253,029  (100.0%)
fracMatch >= 0.90  AND  aln_len >= 1000  66,142  ( 26.14%)   <- Bailey >=90% / >=1kb
   AND uppercase(non-repeat) match-frac >= 0.50  27,623  ( 10.92%)   <- + TE/satellite exclusion
```

(uppercaseMatches / aln_matches: median 0.456; **54.92%** of pairs are majority repeat-masked.)

**Consumer audit (Check G).** `bench/genome_family_def.py:79–93` (`load_sedef_pairs`) keeps **only** the six coordinate columns `f[0..5]` and applies **no** identity/length/TE filter. Our genome family graph is therefore currently built on **all 253,029 candidate pairs at the 50% floor**, so the documented repeat-array over-merge (TRNAV-CAC ×173, an rRNA cluster ×70, repeat-lncRNA arrays, DNFAM0 = 728 members) is a **direct, removable artifact** of not applying the canonical filter above — *not* a property of segdups. The thesis definition must state the predicate `SD(·)` explicitly and treat ≥ 1 kb / ≥ 90% / ¬TE as **part of the predicate**.

### 1.5 The allele/paralog problem, resolved structurally on a haploid T2T

In 2002 the deepest exclusion — **paralogy vs allelism** — was unresolvable from assembly alone and required WSSD depth (depth ~1× = allele, ~2× = duplication). Here it is largely resolved **by construction**: GCF_029281585.2 is a **single-haplotype T2T** assembly, so there is no second haplotype for an *allelic* overlap to arise from. Two distinct correctly-assembled loci are paralogous *by construction*, and (modulo assembly error) every SEDEF self-alignment column difference (col 21 < 1.0) is a **PSV, not an allelic SNP**. On a finished haploid assembly the haploidy does WSSD's specific allele-resolving job, so WGAC-alone (SEDEF) is far closer to a clean SD oracle than it was in 2002.

**Residual gap (honest).** Haploidy removes the *allelic-overlap* artifact WSSD originally targeted, but SEDEF-on-T2T carries **no read-depth term**, so two other assembly-side artifacts survive that only depth could catch: (i) a **collapsed** segdup mis-assembled as one locus is invisible to it (no per-copy copy number / parCN), and (ii) a **false duplication / haplotype-switch mis-assembly** can manufacture a spurious SEDEF self-pair. Neither is an allele, but both mean "every column difference is a PSV" holds only *modulo assembly correctness*. This is the residual allele/copy ambiguity, and it is exactly why we have no gorilla WGS parCN and why **O4 copy-vs-allele stays partly open**.

---

## 2. Relational reframe (interpretive, not an operationalized predicate)

Canzar's aesthetic — clean combinatorial structure, provable statements, no arbitrary thresholds, no `1/k` — motivates **reinterpreting** Bailey's two soft *biological* constants (1 kb, 90%) as **one significance tail against an explicit null**, and — equally important — **naming which knobs are genuinely irreducible**. This section is an **interpretive reframe, not an operationalized predicate**: **no code in this thesis computes the segdup E-value**, `final.bed` is produced by SEDEF's hard-floored error model (§1.3), and the formal lattice of §3 is proved on the **Bailey-thresholded** predicate `SD(·)`, not on the E-value object defined here. The reframe's claim is therefore the modest, defensible one — *two explicit constants become one significance statement (whose tuning is relocated, not removed)* — not "threshold-free."

### 2.1 The relational segdup edge

Fix the haploid genome *G* (GCF_029281585.2) and a per-base neutral / chance-match noise model. Two genomic intervals *A*, *B* would form a **segmental-duplication pair** iff there exists a **maximal** local alignment *a* : *A* ↔ *B* whose alignment statistic is **significant above *G*'s unique-alignment background** — formally its Karlin–Altschul E-value

$$E(a) = K\,|G|^2\,e^{-\lambda S(a)} < \alpha_{\text{GW}}$$

(genome-wide-corrected), equivalently the log-likelihood ratio of the duplication hypothesis vs the unique-background null exceeds the level-α critical value — **AND** neither *A* nor *B* is explained by membership in a high-copy interspersed-repeat family (genomic copy number < *R*).

**What this move does — and does not — buy.** "Maximal" + "significant" let length and identity enter through a *single* monotone score *S(a)* instead of two standalone cliffs: a short match carries a large E-value (chance), a long 86% gorilla-specific dup can be significant while a 300 bp 96% match is not. But the trade-off surface is **not parameter-free** — it is *set by the scoring scheme*. *S(a)* is computed under a substitution matrix, a gap penalty, and an X-drop/extension cutoff; **those weights ARE the length↔identity exchange rate** the reframe celebrates collapsing, and λ, *K* are themselves derived from them. So the two biological constants are **relocated into the shape of *S* (plus α)**, not eliminated. The honest statement is a *reduction* of two interpretable-but-arbitrary biological numbers to one significance level plus a scoring model — itself a parameter.

### 2.2 The background model — and the dual null that yields read-ambiguity

The null is *G*'s **unique-alignment background**: the distribution of best local-alignment scores a single-copy segment achieves against the rest of *G*. Three operationalizations of the **same significance schema** (each *analogous* to — not identical with — our family-edge gate `min_p = ε^δ`):

1. **Karlin–Altschul E-value.** Under an iid-composition null, P(best local score ≥ S) ≈ 1 − exp(−K|G|²e^{−λS}); a segdup is an upper-tail event (E < α). This is the *homology analog* of the per-read Poisson-binomial p-value — there δ PSV disagreements at error rate ε, here matched bases against a chance-match rate, with a |G|² genome-wide correction (vs the read gate's Bonferroni over n−1 loci). Same **schema**, different data + null + correction.
2. **Minimizer / k-mer uniqueness** (read-free; ties to our minimizer-identifiability predictor). A single-copy locus's minimizers are near-unique in *G* (multiplicity 1); a duplicated, *correctly-assembled* locus's minimizers occur ≥ 2 times. A segdup is a run significantly over-multiplic above the unique mode. This is the **read-free analog of WSSD's** depth-over-representation logic, but it is *not* WSSD: Bailey's WSSD measured independent WGS **read depth vs a reference** (depth ~ copy number, R² = 0.96), whereas this is self-multiplicity of *G* against itself, with no depth term. The two diverge exactly on the collapse case — a *collapsed* segdup is multiplicity-1 in a haploid assembly (§1.5), so R² = 0.96 does **not** transfer to this null. The right prior art is **PSV-multiplicity correlation-clustering (SDA, Vollger 2019)**, which resolves collapsed duplications from long-read self-multiplicity; we borrow its logic, not WSSD's read-depth statistic.
3. **Neutral-divergence / chance-match floor.** Identity significantly above the by-chance-match rate (and divergence, col 23, significantly below the neutral/saturation ceiling) ⇒ shared ancestry, not coincidence — the *dual* of the RNA error floor ε.

**The unification — same schema, two nulls (not one common level α).** One *schema* — *indistinguishability/over-similarity above a noise floor* — evaluated against **two different nulls** that share no common α:

- Against the **genome unique-alignment null** ⇒ **(a) segdup edges** ("significantly *more* similar than unique background").
- Against the **read error-floor null** ⇒ **(c) read-conflict/tie edges** via `min_p` ("*not* significantly more divergent than sequencing error"; p ≥ α fails to reject a tie).

These differ in data (assembly vs reads), noise model (chance-match vs sequencing error ε), and multiple-testing correction (|G|² vs Bonferroni n−1); they are **the same significance schema, not the same number**. Their **set-difference (a) \ (c)** — old/divergent segdups that are significant genomic homology yet carry enough PSVs that reads resolve them — is the regime the break-test probes (SEDEF %identity vs read-PSV-density; §4.3, bounded). So the Venn among (a) genomic, (b) exonic, (c) read-ambiguity is **not ad hoc**: it falls out of **which null** the shared schema is run against.

### 2.3 Residual thresholds — stated, not hidden

**Reinterpreted** (two explicit constants → one significance tail, tuning relocated):

- (i) the **1 kb length floor** and (ii) the **90% identity floor** → no longer two standalone biological cliffs but a single significance surface in *S(a)* (length↔identity traded continuously).

This is a genuine *reduction in the number of ad hoc biological constants* (two → ~one statistical level). It is **not** a removal of all tuning — see below.

**Irreducible** (do **not** claim threshold-freedom):

1. **The scoring scheme** — the substitution matrix, gap penalty, and X-drop/extension cutoff that define *S(a)*. This is the **knob that encodes the very length↔identity trade-off §2.1 collapses two constants into**; it cannot be omitted from this list. (It is the one knob the earlier draft of this note wrongly left out.)
2. **α_GW**, the genome-wide significance **level** — a chosen Type-I / FDR rate. It is *analogous* to (not the same number as) our family gate's Bonferroni α/(n−1), is statistically interpretable (expected false segdup pairs), and *calibratable* rather than an ad hoc biological constant.
3. **The repeat copy-number cutoff *R*** — the "not a high-copy interspersed repeat" clause is **irreducibly** a threshold on copy number / repeat-family size. There is **no threshold-free boundary** between a ~200-copy gorilla-specific segdup family and a low-copy transposon family; separating "locally-duplicated" from "interspersed-mobile" needs either *R* or an orthogonal RepeatMasker/Dfam class label (itself thresholded). **This is a genuine ontological cut, not a tunable knob** — and it is exactly where our catalogs over-merge (TRNAV-CAC ×173, rRNA ×70, repeat-lncRNA arrays, DNFAM0 = 728) and exactly what Bailey buries in "EXCLUDING high-copy interspersed repeats." **This is the one place threshold-freedom most honestly fails.**
4. **k / minimizer window w** — a seed-resolution scale for the uniqueness null (a nuisance parameter; stable across reasonable k but not parameter-free).
5. **The gene-node COVER fraction** (our COVER = 0.50) — a downstream annotation-**projection** choice for turning interval-pairs into gene edges; orthogonal to homology significance and **not** dissolved by this reframe.

**Reality check (do not overclaim).** No code computes the E-value; `final.bed` is hard-floored by SEDEF (`max-error 0.5`, `min-len 1000`, `max-div 0.30`; §1.3), and in *every* actual run we recover the operational SD by applying ≥ 1 kb / ≥ 90% / ¬TE (only **26.1%** of the 253,029 pairs survive; a significance/identity gate at *some* α still prunes ~4×). The reframe replaces a **hard cliff with a smooth tail in principle**, and reduces the *biological* constant count, but it relocates the trade-off into the scoring scheme + α and leaves *R* irreducible. It does **not** make a cut vanish, and §3 below is proved on the thresholded `SD(·)`.

### 2.4 Why this fits Canzar (and where the honest line is)

- **Clean combinatorial object — two oracles, one node set.** Both definitions are connected-components (|C| ≥ 2) over the **same** node set *L*; they differ **only** in the edge oracle (genomic self-alignment vs expressed-read ambiguity). The object of interest is therefore their **agreement region $E_a \cap E_c$ and the lattice relating them (§3.3)** — *not* a graph homomorphism. (An earlier draft asserted "O1 corroboration is a homomorphism between two instances of one predicate"; that line is **withdrawn** — it is refuted by the §3.3 incomparability proof, since $E_c \not\subseteq E_a$ via APOBEC3D/F means the identity map on *L* preserves no edge in that direction. The honest statement is "two distinct edge oracles on a shared node set, whose overlap and non-overlap we characterize.")
- **Interpretive parallel to our certifiers (not a proven equivalence).** The E-value plays, *by analogy*, the role our proven certificates play on the read side — Thm 4 (`min_p = ε^δ`, a sound per-read certifier ∀ K ≥ 1) and Thm 7 (Strong-Sep ⇒ LP integral + `min_p` = dual witness). This parallel is **conceptual**: because no code computes the segdup E-value, there is no machine-checked "exact recovery + dual witness" on the genome side, only on the RNA side.
- **No 1/k:** significance is a per-pair likelihood ratio against an explicit null, never a uniform 1/copies prior.
- **Honesty as a feature:** the contribution is a **reduction** (two biological constants → one statistical level + a scoring model) and a **schema unification** (one significance logic spanning genome segdup definition, RNA family grouping, per-read assignment), **not** a zero-threshold claim. Naming the scoring scheme and *R* as irreducible is exactly the rigor Canzar rewards.

---

## 3. Graph formalization and the three-way relationship

### 3.1 The segdup graph (three layers)

All on `bench/genome_family_def.py` + `final.bed`. The combinatorial skeleton is **identical** to the RNA read-conflict definition (nodes = loci, edge = a homology tie, family = connected component |C| ≥ 2; `conflict_families()` in `read_conflict.rs:205–239` is the *same* union-find). **The only thing that differs across the two definitions is the EDGE ORACLE** — genomic self-alignment (read-independent) vs expressed-read ambiguity (read-derived). That shared skeleton is what lets the two be overlaid without changing the algebra.

**Layer 0 — raw homology relation (genomic intervals).** Let *S* be the SEDEF self-alignment. Each row is an ordered pair of intervals (I, J) with ℓ(I,J) = col 12, ρ(I,J) = col 21, d = col 23. Define the **Bailey-segdup predicate**

$$\mathrm{SD}(I,J) \iff \ell(I,J)\ge 1000 \;\wedge\; \rho(I,J)\ge 0.90 \;\wedge\; \neg\mathrm{TE}(I,J) \;\wedge\; I\neq J,$$

symmetric and irreflexive ⇒ undirected graph $G_{\mathrm{SD}}^0$. **Empirical caveat (load-bearing):** `final.bed` on disk is **not** pre-filtered to `SD(·)` (only 66,142/253,029 = 26.14% satisfy ℓ ≥ 1 kb ∧ ρ ≥ 0.90), and `load_sedef_pairs` imposes no filter, so the edge oracle **as-coded** is a *superset* of Bailey-segdup. Imposing `SD(·)` is precisely what dissolves the Alu/repeat over-merge — the operational failure of Bailey's TE-exclusion clause.

**Layer 1 — annotated-gene graph (Catalog A).** Nodes $V_A$ = NCBI RefSeq loci. Edge (COVER = 0.50):

$$\{u,v\}\in E_A \iff \exists\,(I,J)\in S:\ \mathrm{SD}(I,J)\wedge \mathrm{cov}(u,I)\ge\tfrac12\wedge\mathrm{cov}(v,J)\ge\tfrac12,\quad \mathrm{cov}(g,X)=\tfrac{|g\cap X|}{|g|}.$$

A multi-copy gene family = connected component with |C| ≥ 2 (union-find). Recovers MAGEA1/MAGEA12/CSAG1, CEACAM5/6/7, PRSS1/2/3, KRAB-ZNF, IFITM1/2, ULBP1/3 **from the genome alone**. Cross-chromosome pairs are first-class (edge emission sharded by side-A contig, merged in one global union-find).

**Layer 1′ — annotation-free duplication-block graph (Catalog B).** Nodes $V_B$ = segdup-footprint **units** (maximal merged SEDEF footprints per contig). Edge $\{U,W\}\in E_B$ iff ∃ (I,J) ∈ S with I ⊆ U ∧ J ⊆ W. A **duplication block** = component ≥ 2 units, *then* labeled by ≥ 50%-covered genes; the graph is built **without** the annotation (the read-independent analog of de-novo family discovery). Granularity: most A-families sit inside a single B-block (B is coarser, merging adjacent paralog families sharing one duplication unit).

### 3.2 Three notions on one ground set

Put all three on the **same** ground set: unordered locus pairs {u,v} over the `genome_family_def` node set *L*. Two of the three carry an explicit operational criterion ($E_a$ via `SD(·)`, $E_c$ via the code predicate); **$E_b$ is the loosely-defined member** and must be pinned down — and, to make the one theorem true, pinned down **asymmetrically**.

$$
\begin{aligned}
E_a\ (\textbf{genomic segdup}) &= \{\{u,v\}: \exists(I,J)\in S,\ \mathrm{SD}(I,J),\ \mathrm{cov}(u,I)\ge\tfrac12,\ \mathrm{cov}(v,J)\ge\tfrac12\}\\
E_b\ (\textbf{exonic homology, asymmetric}) &= \{\{u,v\}: \exists\ \text{a segment that is exonic/expressed at }\ge 1\ \text{endpoint and}\\
&\qquad\quad \text{genomically homologous (alignable, divergence}\le \text{de\_max) at the other}\}\\
E_c\ (\textbf{read-ambiguity}) &= \{\{u,v\}: \ge\!\min\text{-reads reads place on BOTH }u,v\ \text{with divergence the gate cannot significantly distinguish}\}
\end{aligned}
$$

$E_c$ = `conflict_edges` in `read_conflict.rs` (de-tie or sig-tie). **Two honesty points about $E_b$:** (i) it carries **no genome-wide significance gate of its own** (unlike $E_a$'s E-value-thresholded `SD(·)` or $E_c$'s `min_p`), so it is the **loosest** of the three circles — the single clean containment $E_c \subseteq E_b$ below is robust **partly because** $E_b$'s target set is defined permissively. (ii) It is **asymmetric on purpose**: the shared segment need be exonic/expressed only at *one* endpoint and merely *genomically* homologous at the other. This is forced by the code: the de-tie predicate (`read_conflict.rs:86–106`, `de_tied`/`sig_tied`) uses **only** per-placement divergence `de` and `aln_len` and **never checks that the read's placement at the partner locus overlaps an exon/expressed region** (grep-confirmed: the only `exon` strings are doc comments). A read from expressed locus *u* can multimap onto a homologous but **intronic / intergenic / unexpressed-paralog** stretch of *v* — exonic at *u*, not at *v*. A *symmetric* "exonic at both" $E_b$ would therefore **not** contain $E_c$; the asymmetric definition is what makes the theorem hold.

### 3.3 The lattice: ONE containment, two genuine incomparabilities

> **The naïve claim "$E_c \subseteq E_a$" (read-conflict family ⊆ segdup family) is FALSE** — and is broken by a one-grep witness in our own data (APOBEC3D/F, below).

**(1) $E_c \subseteq E_b$ — THE ONE THEOREM (under the asymmetric $E_b$).** A read can de-tie *u* and *v* only if it **aligns to both** with comparable divergence ≤ de_max (`sig_tied`/`de_tied`, `read_conflict.rs:86–106`). The read's sequence is expressed (it *is* an RNA read at its origin locus) and is homologous to *v*'s genomic sequence at the placement — i.e. exonic/expressed at ≥ 1 endpoint and genomically homologous at the other ⇒ {u,v} ∈ $E_b$ **by the asymmetric definition of §3.2**. So **read-ambiguity ⟹ (asymmetric) exonic homology, unconditionally** — *with the caveat that the shared stretch need not be exonic at the partner locus* (the code enforces no partner-exon overlap). Under a *symmetric* "exonic at both" reading the containment would **fail**; the theorem is exactly as strong as $E_b$'s asymmetry. Refinement inside $E_c$: the significance edge $E_c^{\mathrm{sig}} \subseteq E_c^{\mathrm{de}}$ (machine-checked, `sig_edge_is_a_refinement_of_de_tied`; sig can only split/shrink families — real GGO 81 → 71). **This is the only ⊆ in the lattice.**

**(2) $E_a$ vs $E_c$ — INCOMPARABLE (neither ⊆).**
- **$E_c \not\subseteq E_a$ — VERIFIED WITNESS: APOBEC3D ↔ APOBEC3F.** They form a read-conflict edge (`family_def_readconflict.tsv`: `n_conflict = 6`, `mut_cov = 0.915`, 35 / 50 reads) ⇒ ∈ $E_c$. Yet the **only** SEDEF self-pair covering both genes (`NC_086018.1:37,019,503–37,048,982 ↔ 37,033,411–37,085,045`, aln = 54,107) sits at **fracMatch = 0.884 < 0.90**, so it **fails the Bailey identity floor** ⇒ {APOBEC3D, APOBEC3F} ∉ $E_a$ (it survives only in the raw 50%-floor `final.bed` superset, not in `SD(·)`). This is a clean **sub-threshold-identity paralog**: read-ambiguity does not respect the 90% cliff. (A second, *structural* route — a **processed retrocopy**: single-exon, sub-kb exons ⇒ no contiguous ≥ 1 kb genomic block, yet reads cross-map — holds **in principle** and motivates our intronlessness-keyed retrocopy filter, but we have **not** byte-verified a clean single-exon instance in the gorilla T2T annotation. **RABL2 is *not* such a witness** — see (3) below — so the structural retrocopy break is currently *expected-in-principle, not isolated in our data*: an open verification item.)
- $E_a \not\subseteq E_c$: **(dominant mass)** the overwhelming majority of SEDEF SDs cover no expressed locus, so they are trivially ∉ $E_b$ ⊇ ∉ $E_c$; **(graded)** among expressed SDs, a **divergent/old** one carries enough PSVs that reads place uniquely ⇒ ∉ $E_c$ — the regime the break-test probes (bounded; §4.3).

**(3) $E_a$ vs $E_b$ — INCOMPARABLE.**
- $E_a \not\subseteq E_b$: an **unexpressed / non-genic** SD is genomic homology with no transcript ⇒ ∉ $E_b$. Also, COVER = 0.50 is on **genomic span** (mostly intronic), so an SD can cover ≥ 50% of a gene's span while overlapping few exons ⇒ genomic edge without exonic homology (the **span-vs-exon gap**, a structural honesty point about $E_a$).
- $E_b \not\subseteq E_a$: two paralogs sharing **one conserved exon/domain** have (partial) $E_b$ but no ≥ 1 kb / ≥ 90% SD covering both genes ≥ 50% ⇒ ∉ $E_a$ (verified domain-sharers below, all 0-conflict).

**Lattice summary (precise):**

$$\boxed{\,E_c^{\mathrm{sig}}\subseteq E_c^{\mathrm{de}}\subseteq E_b^{\mathrm{asym}}\,}\qquad\text{(the only chain; $E_b$ asymmetric — exonic at $\ge1$ endpoint)}$$
$$E_a\not\subseteq E_b,\quad E_b\not\subseteq E_a,\quad E_a\not\subseteq E_c,\quad E_c\not\subseteq E_a\qquad\text{(four non-containments; $E_c\not\subseteq E_a$ verified via APOBEC3D/F)}$$
$$\text{Triple core }E_a\cap E_b\cap E_c=\{\text{recent, gene-bearing, expressed, exonically near-identical SDs}\}\supseteq\{\text{RABL2A/B}\}\supseteq\text{the K-frontier.}$$

**Conditional partial implication:** $E_a$ restricted to gene-bearing SDs whose copies are **expressed** and whose footprint **covers exons** ⇒ $E_b$; without those conditions it does not.

**Over-merge is NOT cured by `SD(·)` alone.** §1.4 / §3.4 attribute the repeat-array over-merge to the missing `SD(·)` filter — true for the *repeat-bridged* edges (TRNAV-CAC, rRNA, DNFAM0), which the TE/identity floors remove. But connected-components **single-linkage over-merges transitively even after `SD(·)` is imposed**: a mosaic SD sharing legitimate ≥ 90% block-1 with family X and block-2 with family Y chains X–Y through two *individually valid* SD edges. This is precisely why the pipeline carries a **≥ 2-distinct-loci / homology-component refinement** ("raw over-merges via Alu-bridge → refine"); `SD(·)` is **necessary but not sufficient**, and the formal components should be read as *post-refinement*, not raw union-find.

So the correct picture is **not** a subset chain segdup ⊇ exonic ⊇ read. It is: a read-ambiguity edge always lives inside (asymmetric) exonic homology ($E_c \subseteq E_b$), while genomic segdup ($E_a$) is an **orthogonal third circle** that overlaps both but contains neither and is contained by neither — with APOBEC3D/F (∈ $E_c$, ∉ $E_a$) and the unexpressed SD mass (∈ $E_a$, ∉ $E_c$) as the two witnessed escapes.

### 3.4 Named counterexamples (each witnessed in our data)

- **SUB-THRESHOLD-IDENTITY PARALOG — $E_c \setminus E_a$ (VERIFIED).** **APOBEC3D ↔ APOBEC3F** (both on `NC_086018.1`, ~7 kb apart): a read-conflict edge (`family_def_readconflict.tsv`: `n_conflict = 6`, `mut_cov = 0.915`) ⇒ ∈ $E_c$ (and ∈ $E_b$). The only SEDEF pair covering both genes is **88.4% identity** (`37,019,503–37,048,982 ↔ 37,033,411–37,085,045`, aln 54,107) — **below the Bailey 90% floor** ⇒ ∉ $E_a$ under `SD(·)`. This is the **one-grep break of the naïve nesting**: genomic homology exists but falls just under the arbitrary identity cliff, while reads still cannot fully distinguish the copies. **The retrocopy route** ($E_c \setminus E_a$ via single-exon, sub-kb-exon loci with no contiguous ≥ 1 kb genomic block) is *expected in principle* and motivates our intronlessness-keyed retrocopy filter (Bailey: "excluded genes that showed no evidence of intron-exon splicing to avoid… processed pseudogenes"), but is **not byte-verified** in the gorilla T2T annotation here — flagged open.
- **TRIPLE-CORE, RESOLVABLE — $E_a \cap E_b \cap E_c$ (VERIFIED): RABL2A ↔ RABL2B.** Both are ~14–16 kb **multi-exon** genes (`RABL2A` NC_073235.2:15,131,653–15,147,533; `RABL2B` NC_086018.1:48,818,440–48,832,011), **not** retrocopies. They form a **24,995 bp, 98.6%-identity interchromosomal SD** (`final.bed` NC_073235.2:15,128,169–15,153,011 ↔ NC_086018.1:48,812,020–48,835,563, jck 0.014, cov ≈ 1.0) ⇒ ∈ $E_a$ (and recovered together in GFAM0/GFAM1); share exonic sequence ⇒ ∈ $E_b$; and carry a strong read-conflict edge (`n_conflict = 190`, 199/196 reads) ⇒ ∈ $E_c$. So RABL2 sits in all three circles. **Correction note:** an earlier draft (and `denovo_pipeline.rs:951,1156`, "RABL2B retrocopy ~70% genomic") filed RABL2B as a retrocopy ∉ $E_a$; that is **genome-false for this assembly** (a 25 kb / 98.6% genomic SD, two full multi-exon genes) and is **retracted**. RABL2 is triple-core but **resolvable** (≈ 350 PSVs over 25 kb; copy-assignment "RABL2 2cp 100%"), so it is *not* at the unresolvable K-frontier limit — it is the regime where the conflict edge fires yet χ(H) is still 2-colorable.
- **DIVERGENT-RESOLVABLE SEGDUP — $E_a \setminus E_c$.** An old SD that *is* a ≥ 1 kb / ≥ 90% duplication ($E_a$) but whose copies carry enough PSVs that reads place uniquely ⇒ ∉ $E_c$. This is the gradient the break-test probes: Spearman(SEDEF %identity, read-PSV-density) = **−0.443** (n = 154, p = 9.1e-09; `genome_rna_overlay_readcontent.json`, *not* `…_breaktest.py`), i.e. **r² ≈ 0.20** — a weak-to-moderate monotone trend (~80% of variance unexplained), and on an axis that is **read-gated but reference-seeded** (PSV columns from a minimap2 `asm20` self-alignment + ≥ 2 alleles + ≥ 3 reads). See §4.3 for the bound.
- **UNEXPRESSED / NON-GENIC SEGDUP — $E_a \setminus E_b$ (hence \ $E_c$).** The overwhelming majority of the 253,029 SEDEF pairs cover no expressed locus (intergenic/intronic/silent) — the **dominant** $E_a \setminus E_c$ mass. Concretely, the span-vs-exon gap: an SD covering ≥ 50% of a gene's mostly-intronic span while overlapping few exons is in $E_a$ but generates little/no $E_b$.
- **SINGLE-SHARED-EXON / DOMAIN-SHARER — partial $E_b$, ∉ $E_a$, ∉ $E_c$ (VERIFIED).** Genes sharing one conserved domain (`family_def_readconflict.tsv` domain-sharer rows ASDURF/ASNSD1, CASP8/FLACC1, CCDC188/ZDHHC8, CDPF1/PPARA, CREB1/METTL21A, GCA/KCNH7, GPR39/LYPD1 — **all `n_conflict = 0`**) have partial $E_b$ but no ≥ 1 kb / ≥ 90% SD covering both ≥ 50% ⇒ ∉ $E_a$; a read over the single shared domain extends into unique flanks and places uniquely ⇒ no de-tie ⇒ ∉ $E_c$. **This is the case the read-conflict criterion was *designed* to exclude and that a naïve sequence-similarity family definition wrongly merges.**
- **K-FRONTIER (UNRESOLVABLE triple core) — $E_a \cap E_b \cap E_c$ with χ(H) collapsed.** Exonically-identical high-copy co-located copies (RFPL4A family, DAZ on the Y AZFc palindrome): recent genomic SD ($E_a$), identical exons ($E_b$), reads that genuinely **cannot** resolve ($E_c$; locus refs NM:i:0, ~0/386 reads resolve). Distinct from RABL2: here PSV density → 0, so χ(H) cannot separate the copies — the hard core where the gate must **abstain**.
- **OVER-MERGE WITNESS — $E_a$ inflation when `SD(·)` is dropped.** DNFAM0 = 728 members, TRNAV-CAC ×173, rRNA ×70 are repeat-bridged pseudo-edges that survive **only** because the on-disk `final.bed` (and `load_sedef_pairs`) omits the TE-exclusion and ℓ/ρ floors; re-imposing `SD(·)` removes them. **But** note §3.3: even with `SD(·)`, single-linkage still over-merges *transitively* through mosaic SDs — `SD(·)` is necessary, the ≥ 2-distinct-loci/homology-component refinement is the rest. **Concrete evidence that the Bailey thresholds are load-bearing, not decorative.**

---

## 4. Connection to our framework

### 4.1 MCC = χ(H) lives in $E_c$, not $E_a$

Copy-assignment decomposes **exactly** over connected components of the **read-conflict** graph $G_c = (L, E_c)$ (`read_conflict.rs:1–20`: "reads never cross-map outside their component, so the assignment problem decomposes exactly with no information lost"). Within a component, *H* is the read-conflict **hypergraph** whose hyperedges are reads (a read = the set of loci it is de-tied across), and the multi-copy partition cover number **MCC = χ(H)** (copy-assignment theory, Lemma 1).

The segdup graph $G_a = (L, E_a)$ is a **different** graph on the **same** nodes — the read-*independent* oracle. The honest three-way relationship says: **assignment is governed by χ(H) over $E_c$, not over $E_a$.** Using $E_a$ as the decomposition unit would be wrong on both sides — it would **drop retrocopies** in $E_c \setminus E_a$, and **needlessly merge resolvable divergent SDs** in $E_a \setminus E_c$.

### 4.2 The K-frontier *is* the triple core

The triple core $E_a\cap E_b\cap E_c$ is not uniform: it spans a **resolvability gradient**. At moderate identity (RABL2A/B, ρ = 0.986 over 25 kb ⇒ ~350 PSVs) the conflict edge fires but χ(H) is still 2-colorable — the copies *are* recovered (assignment 100%). As genomic identity ρ → 1 (most recent SD), PSV density → 0, reads carry no distinguishing column, and **MCC becomes unrecoverable from RNA** (χ(H) collapses copies into one color class; locus refs NM:i:0). The latter limit is the K-frontier:

$$\text{K-frontier}=\{\,(I,J)\in E_a\cap E_b\cap E_c:\ \rho\to 1,\ \text{PSV density}\to 0\,\}=\text{exactly where the gate must ABSTAIN (assign-or-abstain, no }1/k).$$

This is why the K-frontier is genuinely RNA-unresolvable: it is the high-ρ **limit** of the triple core where $E_c$ is forced and χ(H) cannot separate — RFPL4A / DAZ, not RABL2.

### 4.3 The break-test *probes* $E_a \setminus E_c$ — a bounded trend; the decisive guard is now RUN (read gate **inert**)

The SEDEF-identity-vs-read-PSV-density correlation is **Spearman −0.443** (n = 154, p = 9.1e-09), and its actual source is **`genome_rna_overlay_readcontent.json`** (the `psv_per_kb` axis), **not** `genome_rna_overlay_breaktest.py` (an earlier citation error). Read precisely, this is **r² ≈ 0.20** — a weak-to-moderate monotone trend (~80% of variance unexplained), directionally consistent with $E_c \approx E_a \cap \{\text{higher identity / lower PSV density}\}$: lower SEDEF identity (older SD) ⇒ more PSVs ⇒ reads resolve ⇒ edge in $E_a$ but out of $E_c$.

**Three honest bounds (do not over-read it as a clean map of $E_a \setminus E_c$):**
1. **Effect size is modest.** −0.443 ⇒ r² ≈ 0.20; the correlation is real and **robust to confounds** (survives partial correlation controlling reads/copies/length, ρ = −0.412; n_copies==2 stratum, −0.456; a join-permutation null, p = 2e-4), but it explains a minority of the variance.
2. **The axis is read-GATED but reference-SEEDED.** The PSV columns are defined by a **minimap2 `asm20` self-alignment of the reference copy fasta**, gated by ≥ 2 alleles and ≥ 3 covering reads — themselves hard thresholds, and substrate that **shares reference bytes with SEDEF**. The provenance note in the artifact is explicit ("psvs/K = read-GATED, reference-SEEDED (asm20)"), and the **purest read-BASE axis (assignment rate, no reference-alignment preset) is weak/non-significant: ρ = −0.138, p = 0.087**. So the strong number is *not* a fully read-discovered corroboration.
3. **The decisive refutation guard has now been run — the read gate is INERT (2026-06-30, `genome_rna_overlay_breaktest.json`, 92 min over the 154 families).** It recomputes the PSV columns **without** the read-coverage gate — a pure `asm20` reference-divergence count — and re-correlates with SEDEF identity. The read gate keeps only **36%** of the asm20-divergent columns (9,089 of 30,893), yet removing it barely moves the correlation: **gated ρ = −0.443 → ungated ρ = −0.420** (p = 6e-08), and the gated vs ungated per-family PSV counts correlate at **ρ = +0.908**. So the `psv_per_kb` axis is a **reference quantity** (minimap2 `asm20` self-alignment), and the read requirement adds essentially nothing: **the −0.44 headline is a two-reference-aligner concordance (SEDEF lastz vs minimap2 asm20), carrying no read information.** This is consistent with the purest read-base axis (assignment rate) being non-significant (ρ = −0.138). §4.3 is therefore **not** evidence that $E_c$ carries read content beyond $E_a$; it is a (real, robust) *within-$E_a$* concordance of two reference aligners on the identity gradient.

The verified non-containment $E_c \not\subseteq E_a$ therefore rests **entirely** on the **APOBEC3D/F witness** (§3.3/§3.4), **not** on the break-test. The break-test, now completed, in fact **confirms** that its $E_a \setminus E_c$ "gradient" is a reference–reference concordance — so it corroborates $E_a$'s internal divergence structure via a second aligner, and says **nothing** read-derived about $E_c$. Net: our two attempts at a read-content cross-modal O1 corroboration (the POA-substring overlay `GENOME_RNA_OVERLAY.md` and this PSV-resolvability overlay) **both reduce to shared reference substrate**; a genuinely read-independent O1 corroboration must come from read-**discovered** divergence (read-consensus POA, not asm20-seeded) or from structural witnesses like APOBEC3D/F, not from either overlay.

### 4.4 Net role of each object

- **$E_a$ (segdup graph)** = the read-independent **genome oracle** that validates our families **non-circularly** (Catalogs A/B recover MAGEA / CEACAM / PRSS from the genome alone).
- **$E_b$ (asymmetric exonic homology)** = the **loosest** circle and the necessary **envelope** that read-ambiguity lives inside ($E_c \subseteq E_b$, the one theorem — exonic at ≥ 1 endpoint, genomically homologous at the other).
- **$E_c$ (read-conflict graph)** — refined by the significance gate ($E_c^{\mathrm{sig}} \subseteq E_c$, `min_p = ε^δ`, Thm 4) — = the actual **carrier of χ(H) = MCC**.
- **Triple core** = their intersection (RABL2A/B verified inside it; the K-frontier its unresolvable limit); the **break-test** gives a bounded probe of $E_a \setminus E_c$.

**The formal contribution** is the lattice $E_c^{\mathrm{sig}} \subseteq E_c \subseteq E_b^{\mathrm{asym}}$ with $E_a$ **incomparable** — replacing the false "read-conflict ⊆ segdup" nesting with an honest Venn whose only containment is *read-ambiguity ⟹ (asymmetric) exonic homology*, with **both** non-containments against $E_a$ witnessed in our own data (APOBEC3D/F ∈ $E_c \setminus E_a$; the unexpressed-SD mass ∈ $E_a \setminus E_c$), and identifying $E_a \cap E_c$ as the only region where the genomic and expressed oracles agree.

---

## Appendix A — verification log (2026-06-30)

| check | command / source | result |
|---|---|---|
| row count | parse `final.bed`, skip `#` | 253,029 data rows |
| identity floor | col 21 stats | min 0.5071, p1 0.692, p5 0.770, median 0.8696, max 1.0; **< 0.90 = 72.41%** |
| length | col 12 stats | min 900, median 4,301, p99 123,849, max 3,261,496; < 1 kb = 3.43% |
| divergence | col 23 stats | min ≈ 0, median 0.1433, p95 0.275, p99 0.397, max 0.803 |
| geography | cols 1,4,10 | inter 84.28%, intra 15.72%, inverted 50.75%, self-pairs 1,939 |
| Bailey filter | ρ ≥ 0.90 ∧ ℓ ≥ 1000 | **66,142 (26.14%)** |
| + TE exclusion | + uppercaseMatches/aln_matches ≥ 0.50 | **27,623 (10.92%)**; 54.92% majority-repeat |
| row-0 encoding | cross-checks | fracMatch = 5256/5543; jck = −¾ln(1−4⁄3·0.0518); transitions+transversions = mismatch; matchB+mismatchB = alnB; max_len = endA−startA — all pass |
| consumer audit | `genome_family_def.py:79–93` | keeps only `f[0..5]`; **no ρ/ℓ/TE filter**; COVER = 0.50 (line 44) ⇒ repeat-bridged over-merge is a removable artifact (transitive single-linkage over-merge is **not** — needs refinement) |
| Bailey anchors | `bailey_famcn.md` | "≥ 90% ≥ 1 kb"; R² = 0.96 **WGS-depth↔CN** (a read-depth-vs-reference statistic, *not* assembly self-multiplicity); 16.5% → 5.2%; "last ~40 million years"; "excluded genes that showed no evidence of intron-exon splicing" |
| RABL2A↔RABL2B | `final.bed`, `genes.bed`, `genome_families_*.tsv` | **triple-core, NOT a retrocopy**: 24,995 bp / 98.6% SD ∈ $E_a$ (GFAM0/GFAM1); both ~14–16 kb multi-exon genes; `n_conflict = 190` ⇒ ∈ $E_c$ |
| APOBEC3D↔APOBEC3F | `final.bed`, `family_def_readconflict.tsv` | **$E_c \setminus E_a$ witness**: `n_conflict = 6` ⇒ ∈ $E_c$; only covering SEDEF pair = **88.4% < 90%** ⇒ ∉ $E_a$(`SD·`) |
| read-conflict predicate | `read_conflict.rs:86–106` | `de_tied`/`sig_tied` use only `de`, `aln_len`; **no partner-exon check** ⇒ $E_b$ must be asymmetric for $E_c \subseteq E_b$ |
| break-test source | `genome_rna_overlay_readcontent.json` | ρ = **−0.443** (psv_per_kb, n = 154, **r² ≈ 0.20**); read-base axis (assign rate) **−0.138 NS** |
| **ungated guard (RUN)** | `genome_rna_overlay_breaktest.json` (2026-06-30, 92 min, n = 154) | **read gate INERT**: gated ρ −0.443 → **ungated ρ −0.420** (p 6e-08); gated/ungated PSV counts ρ **+0.908**; gate keeps 36% of asm20 columns (9,089/30,893) ⇒ the −0.44 axis is reference (asm20) concordance, no read info |

**Definitional one-liner for the thesis body:** A *segmental duplication* is **operationally** a SEDEF self-alignment pair with ℓ ≥ 1 kb, ρ ≥ 90%, ¬TE — **interpretively** (§2, computed nowhere in the pipeline) a pair whose maximal local self-alignment is significant above the haploid genome's unique-alignment background, excluding high-copy interspersed repeats. This is **axis (a)**; our families are the connected components (|C| ≥ 2, post-refinement) of the `SD(·)`-thresholded graph. $E_a$ is a third, **incomparable** circle to (asymmetric) exonic homology ($E_b$) and read-ambiguity ($E_c$); the only clean containment is read-ambiguity ⟹ exonic homology ($E_c \subseteq E_b$), and *both* directions of $E_a$-incomparability are witnessed in our data.
