# A Formal Definition of RNA-Level (Transcript) Multi-Copy Gene Families for the Gorilla Pan-Transcriptomic Thesis

**Scope.** This note fixes, for the rest of the thesis, what a **multi-copy gene family at the RNA level** is, as **one canonical combinatorial object** parallel to the three companion notes, and it **corrects a structural flaw in the prior RNA family definition.** The prior object (the `## family_definition_formal` section of `bench/FAMILY_DEF.md`) defined the RNA family as a connected component of the **read-conflict graph** $E_c$ — an edge being a read that *cannot distinguish two loci* (a de-tie / MAPQ-0 ambiguity). $E_c$ is an **ambiguity oracle, not a homology oracle**: by construction it links only loci `minimap2` *fails* to resolve (collapsed / near-identical copies), so a genuine multi-copy family whose copies are divergent enough that every read maps uniquely produces **no de-ties, fragments into singletons, and is invisible as a family** — the RNA analog of the globin problem (MB and HBB are one family that $E_c$ could never group, because no read is shared). This note **reframes the RNA family as the transcript-homology component $E_r$** (which *includes* the read-resolvable copies), and **demotes $E_c$ + the SUN 3-tier resolvability ladder to the within-family O2 (copy-assignment) structure.** The payoff is that all four levels become **uniform homology oracles** — $E_a$ (genomic) / $E_b$ (exonic) / $E_r$ (transcript) / $E_p$ (protein) — and read-ambiguity is correctly demoted to O2. This composes with, and must not contradict:

- `bench/SEGDUP_DEFINITION_FORMAL.md` — the segdup predicate $\mathrm{SD}(\cdot)$, the genomic oracle $E_a$, the **asymmetric** exonic oracle $E_b$, the read-ambiguity oracle $E_c$, and the lattice $E_c^{\mathrm{sig}}\subseteq E_c\subseteq E_b^{\mathrm{asym}}$ with $E_a$ incomparable. **The one clean containment of that note, $E_c\subseteq E_b^{\mathrm{asym}}$, moves down one level here to $E_c\subseteq E_r^{\mathrm{asym}}$** (same ground set $V_R$).
- `bench/DNA_FAMILY_DEFINITION_FORMAL.md` — the DNA-level family: the cohesion-refinement operator $R$ ($\gamma$-quasi-clique + $\ge 2$-distinct-loci) applied to a raw $E_a$ component. **This note reuses that $R$** (via the pipeline's `refine_families`).
- `bench/PROTEIN_FAMILY_DEFINITION_FORMAL.md` — the protein-level family $E_p$, the fourth oracle, same skeleton, same operator $R$.

The combinatorial **skeleton is shared** across all four (nodes = loci, edge = a homology tie, family = $R$-refined connected component with $|C|\ge 2$, $\Lambda(C)\ge 2$ distinct loci); **the only thing that differs is the edge oracle.** What this note adds is the **transcript oracle $E_r$** and the formal **O1/O2 split**: O1 (family definition) = homology; O2 (copy-assignment) = the read-conflict / SUN / gate resolvability layered *inside* each fixed family. This **partially undoes** the prior "one significance criterion for O1+O2" unification (`project_family_def_readconflict`; §7).

All empirical statements were re-verified on 2026-06-30 against `bench/denovo_families.tsv` (the transcript-homology catalog: **1,130 multi-copy families over 3,636 expressed de-novo loci**), `bench/sun_identifiability.tsv` (154 validated multi-copy families / 412 copies), `bench/psv_graph_genomewide.json` (154 co-located families), `/home/juanfra/winloci_scratch/GGO_mm.bam` (the `-N50`, 63%-secondary map) and `/home/juanfra/winloci_scratch/denovo_transcripts.meta.tsv`. Load-bearing counts are inlined (Appendix A).

> **One-line thesis claim.** An RNA-level multi-copy gene family is a **cohesive community of the transcript-homology graph** $G_R=(V_R,E_r)$ satisfying the **same cohesion certificate** as the DNA and protein levels: a block $C=R(\kappa)$ that the operator $R$ carves from a connected component $\kappa$ of $G_R$, a $\gamma$-quasi-clique ($\rho_{\mathrm{in}}(C)\ge\gamma=0.20$) with **$\ge 2$ distinct expressed loci** ($\Lambda(C)\ge 2$). The edge oracle $E_r$ is **significant spliced-transcript self-alignment** (interpretively a Karlin–Altschul E-value $<\alpha_r$ against the transcriptome chance-match null; operationally the shipped POA contiguous-core reciprocal coverage $\ge\gamma_{\text{core}}=0.13$). **This family includes every copy whose transcript is homologous, whether or not reads can resolve it.** Read-ambiguity ($E_c$) and the SUN 3-tier ladder are **not** the family definition; they are the **within-family O2 resolvability structure**, and $E_c\subseteq E_r^{\mathrm{asym}}$ (the RNA analog of $E_c\subseteq E_b^{\mathrm{asym}}$, machine-checked below). Under the old $E_c$-as-definition, **~30% of multi-copy homology families and ~1/4 of copies formed no within-family de-tie edge and were silently dropped as "easily solvable"** (§2 — this de-tie edge-drop is the literal O1-drop). *One adjudication fixed here up front: the **82–96% SUN-Tier-1 / MAPQ$>0$ figures** quoted later measure **per-copy O2 resolvability** (how many individual copies reads can route), a **distinct** quantity that does **not** bound the O1 family-drop — measured, **132/132 = 100%** of the co-located fully-Tier-1 families still form a de-tie edge and SURVIVE $E_c$, so the honest O1-drop is the de-tie **29.8%**, never the ~86% one gets by conflating the two (§5).* This reframe recovers the dropped ~30% at O1 while keeping the collapsed ~18% as the hard O2 core where the gate abstains.

---

## 1. The flaw in the read-conflict definition (why $E_c$ cannot be the family)

The prior RNA family (`FAMILY_DEF.md` §`family_definition_formal`) is the connected component of the read-conflict graph $G_c=(V_R,E_c)$, $|C|\ge 2$, where

$$
\{u,v\}\in E_c \iff \big|\{\,r: r \text{ de-ties } u,v\,\}\big|\ge \mathrm{MIN\_READS}=3,\qquad
r\text{ de-ties }u,v \iff |de_u-de_v|\le\Delta=0.005\ \wedge\ \max(de_u,de_v)\le\mathrm{DE_{max}}=0.05,
$$

(`src/rustle/vg_family/read_conflict.rs:86–106`, `de_tied`/`sig_tied`). This is a beautiful **O2** object — it is exactly the unit over which copy-assignment decomposes (`MCC = χ(H)`, Lemma 1) — but it is the **wrong object for O1 (defining a family)**, for one structural reason:

> **$E_c$ is an ambiguity oracle, not a homology oracle.** An edge exists **iff `minimap2` cannot resolve the two loci**. A multi-copy family whose copies are **divergent enough to each map uniquely** produces **no de-tie edge**, fragments into singletons, and vanishes as a family — even though its copies are unambiguously homologous.

This is not a hypothetical: it is the RNA mirror of the globin case that motivated the protein oracle (`PROTEIN_FAMILY_DEFINITION_FORMAL.md` §4.1). MB↔HBB are one family down to 24% aa identity, yet no read is ever shared between them, so $E_c$ (and any read-ambiguity definition) can never group them. On the GGO substrate the effect is **quantitatively dominant**: **99.49%** of primary reads over $E_r$ loci are MAPQ$>0$ (uniquely resolvable), and genome-wide **99.96%** of primary reads are MAPQ$>0$ (Appendix A). So the read-ambiguity **mass** — the entire node set $E_c$ can ever touch — is a small co-located core (~7–18% of reads, §2, §5). This per-read/per-copy resolvability is **not** itself the O1-drop: the family boundary is set by whether $E_c$ forms **any** within-family de-tie **edge**, and the literal count of families with no such edge is the de-tie **29.8%** of §2 — not the 82–96% one gets by (wrongly) equating it with per-copy resolvability.

The fix is to define the family by **transcript homology**, and to recognize $E_c$ for what it is: a **sub-structure inside** the homology family that says which copies the reads happen to confuse. That is an O2 statement, and §3–§5 make it precise.

---

## 2. The missed class, quantified (GGO)

Because the genome-wide $E_c$ catalog is **not materialized** (build-vs-run gap: `bench/family_def_readconflict.tsv` is a 13-row demo), $E_r\setminus E_c$ cannot be read off by diffing two catalogs. Instead we ask, inside each transcript-homology ($E_r$) family, whether $E_c$ would emit **any** within-family de-tie **edge**: an $E_r$ family with **no** such edge is one $E_c$ drops entirely (it fragments into singletons). This is **edge-formation**, *not* per-copy MAPQ-resolvability — a resolvable copy can still de-tie, so the two must not be conflated (§5). Substrate: `bench/denovo_families.tsv` (the $E_r$ catalog) + `denovo_transcripts.meta.tsv` (per-locus coordinates) + `GGO_mm.bam`. Scripts (run from repo root): `bench/rna_reframe_measure_detie.py`, `bench/rna_reframe_measure_sigtie.py`, `bench/rna_reframe_measure_ec_er.py` (each writes its rows to `bench/rna_reframe_*_rows.tsv`, override with `RNA_REFRAME_OUT`).

**Evaluable set.** `denovo_families.tsv` has **1,130 multi-copy families / 3,636 expressed loci**. Excluding the single repeat/ZNF over-merge blob **DNFAM0 (728 loci** — the RNA analog of the DNA 1,547-gene 22q11 blob, dissolved by the $\gamma$-quasi-clique refinement $R$, not a real family), the evaluable core is **1,129 families / 2,908 loci**; restricting to loci with fetchable BAM coordinates gives the **1,126 families / 2,902 copies** the resolvability measurement ran on.

**Dropped by the shipped $E_c$ criterion** (reproducing `read_conflict.rs::de_tied` bit-for-bit: per-locus best/min `de`, pair de-tied iff $|de_i-de_j|\le0.005 \wedge \max\le0.05$, edge iff $\ge3$ de-tie reads; family EC-visible iff $\ge1$ edge):

| criterion | families dropped (no within-family edge) | copies dropped |
|---|---|---|
| **de-tie** (shipped default, $\Delta=0.005$, $\mathrm{DE_{max}}=0.05$, min 3) | **335 / 1,126 = 29.8%** | **741 / 2,902 = 25.5%** |
| **sig-tie** ($E_c^{\mathrm{sig}}$, $\varepsilon=0.001$, $\alpha=10^{-3}$) | **370 / 1,126 = 32.9%** | **828 / 2,902 = 28.5%** |

Sig is stricter than de (consistent with $E_c^{\mathrm{sig}}\subseteq E_c^{\mathrm{de}}$; it drops slightly more). So **~30% of multi-copy homology families / ~1/4 of copies are silently reclassified as non-families** by the read-conflict definition.

**The drop is size-graded** — exactly the signature of the flaw (small families are more likely to be all-divergent; large ones almost always contain one near-identical pair that de-ties):

| family size | fraction dropped |
|---|---|
| 2 | 284 / 855 = **33%** |
| 3 | **24%** |
| 4–5 | ~15% |
| $\ge 7$ | **0%** |

**The dropped class is genuine paralogy.** Dominant symbol roots among dropped families: **ZNF, ZBTB, TRIM, KRAB-ZNF, TMEM, PNMA, SLC4A**; **266/335** dropped families contain protein-coding genes, and critically **48 are unnamed / novel de-novo families** — invisible to any annotation *and* to $E_c$. These are the globins-style cases: divergent-but-homologous copies each mapping uniquely.

**A per-locus MAPQ=0 (hard-tie) reading measures per-copy O2 resolvability, NOT the O1-drop — and must not be read as a larger drop.** The shipped $E_c$ edge is the **de-tie** predicate ($|de_i-de_j|\le0.005$), *not* MAPQ=0 (both placements exactly co-best); a de-tie is far more permissive than a hard tie, so the O1-drop is the de-tie **29.8%** above, not larger. For O2 context only: **0.51%** of primary reads over $E_r$ loci are MAPQ=0 (907,096 reads; 99.49% uniquely resolvable), genome-wide **99.96%** are MAPQ$>0$ (652,221), and the ambiguity concentrates **14×** in the co-located subset (7.4% MAPQ=0 over the SUN catalog, 72,667 reads). This says $E_c$'s reachable node set is a small co-located core; it does **not** say $E_c$ drops that much at O1. *(An earlier helper, `rna_reframe_measure_ec_er.py`, classified $E_r$ families by a $\ge5\%$-MAPQ0 per-locus proxy and reported a **96.6% "EC_DROPPED"**; that proxy is per-locus **hard-ambiguity**, **not** de-tie edge-formation, and **overstates** the O1-drop ~3× — the faithful edge measure is `rna_reframe_measure_detie.py` = **29.8%**. Its output rows label the proxy `EC_DROPPED`, so read it as an O2 hard-ambiguity screen, not the O1 family-drop.)*

**Verdict: the flaw is real and moderate-to-large, not a corner case.** Roughly one third of all multi-copy families are lost — including named paralog groups and 48 annotation-invisible novel families.

---

## 3. The canonical object: $E_r$, same certificate and operator $R$ as DNA/protein

### 3.1 The transcript-homology oracle $E_r$

**Node set** $V_R$ = the de-novo **expressed** loci — the intron-junction-collapsed representative transcripts the assembler emits ($\ge\mathrm{MIN\_READS}=3$ reads, all-canonical junctions; `bench/denovo_assemble_gate.py`$\to$`denovo_families.py`). This is the **same** node set the read-conflict graph $E_c$ uses (both are "de-novo expressed loci", `read_conflict.rs`), which is exactly what makes $E_c\subseteq E_r$ a clean **same-ground-set** containment (contrast $E_c\subseteq E_b$, which crossed to the larger 34,114-locus genome node set).

**Edge, interpretive form** (parallel to $E_a$'s Karlin–Altschul E-value and $E_p$'s reciprocal-whole-protein E-value). For expressed loci $u,v$ with exon-union (spliced) transcripts $t_u,t_v$,

$$
\{u,v\}\in E_r \iff \exists\ \text{a maximal local alignment } a: t_u\!\leftrightarrow\!t_v \ \text{with}\ E(a)=K\,|T|^2 e^{-\lambda S(a)}<\alpha_r,
$$

i.e. the spliced-transcript self-similarity is significant above the transcriptome's chance-match background — the **same significance schema** as the other three homology oracles, run against the **transcript** chance-match null.

**Edge, shipped operational form** (what actually runs; hard-floored exactly as SEDEF hard-floors $E_a$ and mmseqs hard-floors $E_p$): the strand-symmetric POA contiguous-core reciprocal coverage

$$
\mathrm{core\_recip}(t_u,t_v)\ge\gamma_{\text{core}}=0.13\qquad(\texttt{denovo\_families.py:T\_CORE},\ \text{canonical/RC POA retry}),
$$

read as "$\ge$13%-length contiguous identical block (a shared conserved exon even if flanks diverge)."

**Honesty (parallel to the segdup / protein notes).** The E-value is the interpretive reframe **computed nowhere**; the shipped oracle is the $\gamma_{\text{core}}=0.13$ hard floor, an arbitrary threshold **named as such** — the RNA analog of $E_a$'s $\rho\ge0.90 / \ell\ge1$ kb and $E_p$'s $c_p=0.50$. The L1 circularity caveat (`denovo_assemble_gate.py:58` builds the transcript sequence from `GGO.fasta` reference substrings) bears on corroboration **independence**, not on the definition: for *defining* families by expressed-locus homology, reference-substring transcripts are exactly the right object.

### 3.2 Family = $\gamma$-quasi-clique-refined $E_r$ component

> **Definition (RNA-level multi-copy gene family).** Let $\gamma\in(0,1]$. An **RNA-level multi-copy gene family** is a vertex set $C\subseteq V_R$ such that
> $$
> C\in R(\kappa)\ \text{ for the connected component }\kappa\text{ of }G_R=(V_R,E_r)\text{ containing it,}\qquad\text{and}\qquad\boxed{\ \Lambda(C)\ge 2\ }
> $$
> where $R$ is the **shared** cohesion-refinement operator (`genome_family_def.refine_families`, the same $R$ imported by the DNA and protein notes), so $C$ is a **$\gamma$-quasi-clique** ($\rho_{\mathrm{in}}(C)\ge\gamma=0.20$ on the $E_r$-subgraph, or $|C|\le2$), and $\Lambda(C)$ is the **$\ge2$-distinct-loci** multi-copy count (reciprocal-$\ge50\%$-overlap locus collapse). The predicate "$\gamma$-quasi-clique with $\Lambda\ge2$ inside a raw $E_r$ component" is the **cohesion certificate**, and **it — not any one partition — is the canonical object**: exact max-$\gamma$-quasi-clique partition is NP-hard, so $R$ emits one seed-fixed byte-identical **witness**.

$R$ is required because single-linkage over the raw homology graph **over-merges via repeat-bridges exactly as on DNA**: `denovo_families.tsv` **DNFAM0 = 728 members chained chr1..chrY** is the RNA analog of the 1,547-gene 22q11 DNA blob and the TRNAV / rRNA segdup arrays. So $E_r$ inherits the identical *same skeleton, same operator, different oracle* treatment. This is **bit-for-bit the DNA/protein skeleton** with the edge oracle swapped:

| | RNA (this note) | DNA | protein |
|---|---|---|---|
| ground set | $V_R$ = expressed loci | 34,114 gene loci | 22,614 protein-coding loci |
| edge oracle | $E_r$ transcript homology (spliced; $\alpha_r$ / core_recip $\ge0.13$) | $E_a$ genomic segdup (SEDEF) | $E_p$ reciprocal whole-protein ($E\le\alpha_p$, cov $\ge c_p$) |
| raw over-merge | repeat-bridge blob (DNFAM0 = 728) | transitive 22q11 blob (1,547) | domain-superfamily hubs (GPCR/ZNF) |
| refinement | $\gamma$-quasi-clique $R$ + $\Lambda\ge2$ | $\gamma$-quasi-clique $R$ + $\Lambda\ge2$ | $\gamma$-quasi-clique $R$ + $\Lambda\ge2$ |
| copy count | $K_R=\Lambda(C)$ | $K_{\mathrm{DNA}}=\Lambda(C)$ | $\Lambda(C)$ |

**Crucially — the whole point of the fix: the $E_r$ family includes read-resolvable copies.** A homology family whose copies `minimap2` places **uniquely** — divergent-but-homologous paralogs that never share a read (the globins are the limiting case) — is a *bona fide* $E_r$ family: it has $E_r$ edges (the transcripts are homologous) even though it has **no** $E_c$ edge, so under $E_c$-as-definition it fragments into singletons and vanishes. That **edgeless-$E_c$** class is the **29.8%** de-tie O1-drop of §2. *(Caveat, measured — do NOT read per-copy SUN-Tier-1 resolvability as "no $E_c$ edge": a single distinguishing SUN moves $de$ by only ~$1/\text{read\_len}\le\Delta=0.005$, so SUN-resolvable copies still de-tie. On the co-located SUN catalog **132/132 = 100%** of fully-Tier-1 families still form a de-tie edge and would SURVIVE $E_c$ — `rna_reframe_validate.py` §D. "Fully-Tier-1" is thus an O2 per-copy fact, not the O1-drop; §5.)* The object is **read-independent of ambiguity** — it depends only on the assembled transcript sequences being homologous. $E_r$ is therefore a genuine homology oracle sitting alongside $E_a/E_b/E_p$, not an ambiguity oracle. (Observed GGO refine goes ~1,130 raw $\to$ homology-component $\wedge\ \ge2$-distinct-loci refined, the same discipline that yielded the 157-family / 54-cross-chrom refined catalog in `COPY_ASSIGN_RECOMPUTE`.)

---

## 4. $E_c\subseteq E_r$ — the one containment (asymmetric, machine-checked)

> **$E_c\subseteq E_r^{\mathrm{asym}}$ holds, logically and empirically, provided $E_r$ is the permissive (local / exonic) homology oracle** — the RNA analog of the segdup note's sole clean containment $E_c\subseteq E_b^{\mathrm{asym}}$.

**Logical.** A de-tie requires the read to align to **both** loci with $de\le\mathrm{DE_{max}}=0.05$ and comparable divergence, i.e. the read (transcribed, exonic at its origin) is $\ge95\%$ identical to both loci over its span $\Rightarrow$ a significant **local** transcript alignment exists between $u$ and $v$ $\Rightarrow \{u,v\}\in E_r$ by any permissive/local definition of transcript homology. This is exactly the $E_c\subseteq E_b$ argument moved from the genome node set to the shared $V_R$.

**The asymmetry (honest, forced by the code).** A read can de-tie two loci over a **single shared exon**: it is exonic/expressed at its origin locus and merely homologous to a partner exon at the other. The read-conflict predicate (`read_conflict.rs:86–106`) checks only per-placement `de` and `aln_len` — it **never** verifies that the partner placement overlaps an exon or that the two *full* transcripts are homologous. So an $E_c$ edge can join two loci whose complete transcripts share **only** that one exon and diverge everywhere else. Under a **symmetric/strict** $E_r$ ("the two transcript sequences align significantly over their length" / high reciprocal coverage), that pair would fail $E_r$ and $E_c\not\subseteq E_r$. Therefore, to make the containment hold, $E_r$ must be defined **permissively and asymmetrically**, exactly as $E_b$ is:

$$
\boxed{\,E_c^{\mathrm{sig}}\subseteq E_c\subseteq E_r^{\mathrm{asym}}\,}\qquad
E_r^{\mathrm{asym}}\ \text{edge} = \exists\ \text{a significant \emph{local} transcript alignment (}\ge1\text{ shared exon / significant segment),}
$$

exonic at $\ge1$ endpoint and homologous at the other, regardless of what fraction of either transcript the shared segment spans — **not** full-length reciprocal homology. The theorem is exactly as strong as $E_r$'s asymmetry, precisely mirroring $E_b$'s deliberate one-sidedness.

**Empirical check** (`rna_reframe_validate.py` §B–§D, mapping the 154 co-located ambiguous SUN families onto the denovo $E_r$ catalog by best coordinate overlap, then reading `denovo_family_edges.tsv` and running the faithful `read_conflict.rs` de_tied probe on each split):

- Among the **132** families where $E_r$ assigns a multi-copy label, **129/132 = 97.7%** have **all** their co-located ambiguous copies inside **one** $E_r$ family (the shipped **symmetric** core_recip$\ge0.13$ operational oracle). 22 more carry no multi-copy $E_r$ label at all — an $E_r$ **under-detection** issue (unevaluable; distinct from an $E_c\setminus E_r$ homology gap), *not* an $E_c$-drop.
- The **3** that split across $\ge2$ $E_r$ families are **operational-oracle shared-exon leaks, NOT "EDGE_LINKED"** (an earlier draft's claim, **retracted**). `denovo_family_edges.tsv` — whose families *are* its connected components — carries **0** cross-family `core_recip` edges touching these loci, yet the loci **genuinely de-tie**: the faithful predicate finds real $E_c$ edges — fam16 DNFAM6~DNFAM612 **20** reads (LOC101143406/144261/129525008, tandem paralogs sharing one conserved exon); fam34 DNFAM1120~DNFAM143 **121** (+DNFAM0 13/10; collapsed Tier-3 LOCs, $n_{\rm psv}=0$); fam65 DNFAM10~DNFAM52 **65** (ZNF425 + paralogs on separate chromosomes sharing the ZNF/KRAB domain). Their only shared homology is a **single conserved exon/domain below the $\gamma_{\text{core}}=0.13$ reciprocal-coverage floor**, so the shipped **symmetric** operational oracle **misses** them (~**3/132 = 2.3%** operational leaks).
- These 3 are **absorbed by the DEFINITIONAL permissive-local $E_r^{\mathrm{asym}}$** — each de-tie read (de$\le0.05$ to both loci) is itself a witness of a significant *local* transcript alignment — so the theorem $E_c\subseteq E_r^{\mathrm{asym}}$ **stands**; only the operational 0.13 oracle leaks. **$R$ cannot repair them** — $R$ only *splits* connected components, never merges — so the earlier "fixed by $R$" claim is retracted (backwards: $R$ would further fragment, not rejoin).

So the containment reads two-tier, honestly: the **definitional/permissive** $E_r^{\mathrm{asym}}$ contains **every** witnessed $E_c$ edge (theorem unrefuted); the **operational** core_recip$\ge0.13$ oracle contains **129/132 = 97.7%**, leaking the **3** single-conserved-exon/domain de-ties above. **Residual the definition must name:** because $\gamma_{\text{core}}=0.13$ is a *reciprocal-coverage floor*, a de-tie over an exon/domain covering $<13\%$ of a large transcript escapes the symmetric operational oracle (the 3 witnessed cases) — which is exactly why the *definitional* $E_r$ edge is specified as **permissive local homology** (segment-significance / E-value over the shared local alignment), not the 0.13 reciprocal floor.

---

## 5. $E_c$ + the SUN 3-tier ladder are the within-family O2 structure (demoted from O1)

Fix an O1 family $C$ (an $E_r$ component, $|C|=K_R=\Lambda(C)$ copies). The O2 resolvability structure lives **inside** $C$ and has three coupled layers:

1. **$E_c$ restricted to $C$ (the read-conflict subgraph $E_c|_C$).** Because $E_c\subseteq E_r$ (§4), $E_c|_C$ is a subgraph of the homology family. It partitions the copies by read-resolvability: copies with **no** $E_c|_C$ edge are read-resolvable (MAPQ$>0$ / unique; "resolvable already MAPQ$>0$"); copies joined by $E_c|_C$ edges are the collapsed / MAPQ-0 confusable core. The **read-conflict components** of $E_c|_C$ are the exact-decomposition units of assignment (property P2: no de-tied read crosses two $E_c$ components $\Rightarrow$ assignment separates with no information lost).
2. **$\mathrm{MCC}=\chi(H)$** over the read-conflict hypergraph $H$ on $C$ (Lemma 1) — a strictly O2 quantity: $\mathrm{MCC}=\chi(H)\le K_R=\Lambda(C)$. Equality holds when reads resolve every copy; the **K-frontier** is the strict case $\chi(H)<\Lambda(C)$.
3. **SUN 3-tier per-copy ladder** (`bench/sun_identifiability.*`, `THEORY.md` §5·SUN) stratifies each copy: **Tier-1** SUN-identifiable (single-read gate-deterministic, $|N(r)|=1$); **Tier-2** hap-vector-unique-only (needs $\ge2$-PSV co-observation); **Tier-3** collapsed/frontier (NM:i:0, $\chi(H)$ collapses the copies, gate certifies $\min_p=1$ = tied, assign-or-abstain abstains).

So the demotion is formal:

$$
\mathrm{O2}(C)=\big(\,E_c|_C,\ \chi(H),\ \mathrm{SUN\text{-}tier}(\cdot),\ \text{assign-or-abstain gate }\min_p=\varepsilon^\delta\,\big)
$$

is a structure **attached to** each fixed O1 family, never the O1 boundary. Using $E_c$ as the family *definition* — the prior flaw — is the operation "take the whole catalog, keep only the non-trivial $E_c$-components, discard every copy whose $E_c|_C$ is edgeless." That discards the resolvable copies (the *good* O2 outcome) **and** the families they constitute.

**O2 copy-resolvability inside these families** (`sun_identifiability.tsv`, 154 validated multi-copy homology families / 412 copies) — a **per-copy** quantity that does **NOT** measure the O1 family-drop (the O1-drop is the de-tie **29.8%** of §2):

| SUN tier | copies | share | O2 read-resolvability (per copy) | O1 de-tie-edge status (measured) |
|---|---:|---:|---|---|
| Tier-1 | 338 | **82.0%** | resolvable (MAPQ$>0$ / SUN-identifiable) | **still forms a de-tie edge — $E_c$-VISIBLE, NOT edgeless** (132/132 fully-Tier-1 families, §D) |
| Tier-2 | 1 | 0.2% | hap-vector-unique-only | marginal |
| Tier-3 | 73 | **17.7%** | collapsed / K-frontier (the hard O2 core) | $E_c$-visible |

**132/154 families (85.7%) are fully-Tier-1** — meaning **every copy is individually read-resolvable** (per-copy O2 resolvability). This is a statement about *copy-assignment*, **not** the O1 family-drop, and the two must not be conflated. "Fully-Tier-1" does **not** imply an edgeless $E_c|_C$: a single distinguishing SUN moves $de$ by only ~$1/\text{read\_len}\le\Delta=0.005$, so SUN-resolvable copies still de-tie. **Measured** (`rna_reframe_validate.py` §D): **132/132 = 100%** of these fully-Tier-1 co-located families **do** form a within-family de-tie edge and would **SURVIVE** an $E_c$-as-O1 definition — **0** vanish. So on this *known-ambiguous, co-located* catalog the $E_c$ O1-drop is ~**0%**, far *below* the full-catalog **29.8%** (as expected — the co-located core de-ties the most). The earlier inference "fully-Tier-1 $\Rightarrow$ edgeless-$E_c$ $\Rightarrow$ vanish $\Rightarrow$ $E_c$ captures $\approx1/8$ / silently drops ~86%" is therefore **retracted**: it is empirically **inverted** (the §4 fam16 split is itself fully-Tier-1 yet carries a dense de-tie subgraph and would survive). The **sole** correct O1-drop headline is the de-tie **29.8%** (§2); 82–96% is per-copy O2 resolvability. Demoting $E_c$+SUN to O2 still recovers all 154 SUN families at O1 while keeping the collapsed **17.7%** Tier-3 as the hard O2 core where the gate abstains. (`psv_graph_genomewide.json`: 154 co-located families, **82.5% per-copy-resolvable / 17.5% frontier(K=0)** — again a per-copy O2 split, not an O1-drop.)

---

## 6. The updated 4-level lattice (all four are HOMOLOGY oracles; $E_c$ is an O2 substructure)

The four **O1** oracles are the abstraction-nested homology levels

$$
E_a\ (\text{genomic DNA})\ \text{—}\ E_b\ (\text{exonic DNA})\ \text{—}\ \boxed{E_r\ (\text{spliced transcript})}\ \text{—}\ E_p\ (\text{protein}),
$$

all **read-independent-of-ambiguity** and pairwise **incomparable** (each has a witnessed set-difference in both directions). $E_c$ (read-ambiguity) **leaves the O1 lattice** and becomes an **O2 sub-structure of $E_r$.**

**$E_r$ vs $E_b$ — a symmetric refinement on a shared axis, incomparable overall (NOT equality).** $E_b$ (`SEGDUP_DEFINITION_FORMAL.md` §3.2) is exonic homology defined **asymmetrically on the genomic substrate**: exonic/expressed at $\ge1$ endpoint AND genomically (DNA-)homologous at the other, over the full 34,114-locus set. $E_r$ is transcript homology defined **symmetrically on the spliced substrate**: both endpoints are expressed de-novo loci whose exon-union transcripts align. They differ on two independent axes:
- **$E_b\setminus E_r\ne\varnothing$ ($E_b$'s asymmetry).** An $E_b$ edge whose partner endpoint is intronic / intergenic / an unexpressed paralog is exonic at only one end; its partner is not an expressed $V_R$ node, so it has no $E_r$ counterpart. On the shared *expressed-both, exonic-at-both* core, $E_r\subseteq(E_b$ projected onto $V_R)$: **$E_r$ is the symmetric ("exonic at BOTH") core of the asymmetric $E_b$.**
- **$E_r\setminus E_b\ne\varnothing$ ($E_r$'s spliced substrate).** $E_a/E_b$ align *genomic* DNA; $E_r$ aligns *spliced* transcripts (introns removed). A **retrocopy↔parent** pair aligns cleanly transcript-to-transcript (both spliced products) yet the intronless retrocopy is one genomic block while the parent's homologous exons are scattered across introns — so the contiguous genomic homology of $E_a/E_b$ fragments or fails. Witnessed in our data: EEF1A1/CNN2 retrocopies, OCLN~SEPTIN7, BCAS4~CCDC30 read-through — $E_c/E_r$ families that are retrocopy / processed-pseudogene / read-through, **not** genomic segdups.

So $E_r\ne E_b$ and $E_r\not\subset E_b$ cleanly: $E_r$ is the symmetric refinement of $E_b$'s exonic homology on the shared expressed core, **plus** a spliced-homology superset direction (retrocopies) that genomic homology misses.

**The clean-containment chain moves down one level from $E_b$ to $E_r$:**

$$
\boxed{\,E_c^{\mathrm{sig}}\subseteq E_c\subseteq E_r^{\mathrm{asym}}\,}\qquad\text{(the new sole unconditional containment; supersedes }E_c\subseteq E_b^{\mathrm{asym}}\text{ as the tightest)}
$$

because $E_c$ and $E_r$ **share the node set $V_R$** and $E_r$ sits *between* $E_c$ and $E_b$: $E_c\subseteq E_r\subseteq_{\mathrm{proj}}E_b$. The $E_c\subseteq E_r$ containment is itself **asymmetric** (§4), exactly as $E_c\subseteq E_b$ was.

**Pairwise incomparabilities (all four homology oracles):**

$$
E_r\ \text{vs}\ E_a:\ \text{retrocopy }E_r\setminus E_a\ +\ \text{unexpressed-SD }E_a\setminus E_r;\qquad
E_r\ \text{vs}\ E_b:\ \text{as above};\qquad
E_r\ \text{vs}\ E_p:\ \text{non-coding lncRNA-dup }E_r\setminus E_p\ +\ \text{twilight-zone globin }E_p\setminus E_r.
$$

**Retained conditional containments:** $E_r\subseteq_{\mathrm{proj}}E_b$ (expressed-both symmetric core); $E_c\cap\{\text{coding-both}\}\subseteq E_p$.

**The whole lattice, with O1/O2 separated:**

$$
\underbrace{E_a\ \text{—}\ E_b\ \text{—}\ E_r\ \text{—}\ E_p}_{\text{O1: homology oracles (pairwise incomparable)}}\qquad\Big|\qquad \underbrace{E_c^{\mathrm{sig}}\subseteq E_c\subseteq E_r^{\mathrm{asym}},\ \chi(H)\le\Lambda,\ \text{SUN 1/2/3}}_{\text{O2: resolvability, inside a fixed }E_r\text{ family}}.
$$

---

## 7. O1/O2 separation, and what this changes vs the prior unification

### 7.1 The clean separation

- **O1 — family definition = HOMOLOGY.** A multi-copy family at level $\ell\in\{E_a,E_b,E_r,E_p\}$ is a $\gamma$-quasi-clique-refined connected component ($|C|\ge2$, $\Lambda(C)\ge2$) of the level-$\ell$ homology graph. The RNA O1 object is the **$E_r$ transcript-homology family**; it **includes every homologous copy, whether or not reads resolve it** — the per-copy-resolvable copies (SUN-Tier-1, ~82% of copies) are first-class O1 members, which is a resolvability fact, not an O1-drop. O1 answers *"are these expressed loci copies of one another?"* via sequence homology against a **chance-match** null.
- **O2 — copy-assignment = RESOLVABILITY.** Inside a *fixed* O1 family $C$, O2 = $(E_c|_C,\ \mathrm{MCC}=\chi(H)\le K_R=\Lambda(C),\ \text{SUN 3-tier},\ \text{gate }\min_p=\varepsilon^\delta\ge\alpha)$. O2 answers *"can reads be routed to individual copies?"* against the sequencing-**error** floor null $\varepsilon$. The **K-frontier** ($\chi(H)<\Lambda(C)$; Tier-3 / NM:i:0) is the O2-unresolvable limit where the gate abstains — genuinely RNA-unresolvable *assignment*, **not** a non-family.
- **Orthogonal, not nested.** A fully-resolvable family is one O1 object with a *trivial* (edgeless-$E_c$, all-Tier-1) O2 — it both **is** a family and is **fully assignable**. A K-frontier family is one O1 object with a *degenerate* O2 (χ(H) collapses, abstain). The old $E_c$-as-O1 definition conflated the two: it demanded a *non-trivial O2* as the price of O1 membership, so families with easy (good) O2 were denied O1 status and dropped. One line: **O1 is what the transcriptome duplicated ($E_r$ homology); O2 is what the reads can tell apart ($E_c$/SUN/$\chi(H)$).**

### 7.2 This partially undoes the prior "one significance criterion for O1+O2"

The prior unification (`project_family_def_readconflict` / `CONFLICT_EDGE_UNIFICATION`) collapsed O1 and O2 onto a **single** significance test: it redefined the family **edge** to be the assignment gate's own significance-tie $E_c^{\mathrm{sig}}$ (link two loci iff $\varepsilon^\delta\ge\alpha$ = "cannot significantly resolve them"), removing the last hand-set constant ($\Delta=0.005$) and declaring "one significance criterion for O1+O2." That made the **family boundary identical to the read-resolvability boundary** — i.e. it **baked the flaw into the formalism**: a divergent-but-resolvable copy pair fails $\varepsilon^\delta\ge\alpha$, is not an edge, and the family fragments. The measured GGO 81$\to$71 families / 205$\to$176 copies "narrowing" *was* this — dropping the pairs "significantly resolvable under the gate."

This reframe re-reads that result and restores the separation:

1. **O1 is no longer the same criterion as O2.** O1 = $E_r$ transcript homology ($\alpha_r$ / core_recip $\ge0.13$, against the transcript **chance-match** null). O2 = the significance-tie / read-conflict / SUN ($\varepsilon^\delta\ge\alpha$, against the sequencing-**error** null). Two different tests, two different nulls — exactly the *"same schema, two nulls (no common $\alpha$)"* structure the segdup note §2.2 already established for $E_a$ (genomic-homology null) vs $E_c$ (read-error null). The unification was the anomaly that had $E_c$ uniquely serving as both homology and definition; the reframe realigns RNA with the DNA/protein/segdup treatment.
2. **What of the unification STANDS (it was right about O2).** Its real content — that the family-edge *refinement* $E_c^{\mathrm{sig}}$ should use the same $\varepsilon^\delta\ge\alpha$ criterion as the assignment gate — is correct **as an O2 statement**: $E_c^{\mathrm{sig}}\subseteq E_c$ is now an *O2-internal* refinement of the read-conflict subgraph, sharing the gate's null. That is **kept**. What was overclaimed was extending that single O2 criterion up to O1 (family *definition*).
3. **Re-read of 81$\to$71.** The 10 "lost" families are still **O1 $E_r$ homology families**; they did not disappear, they merely have a fully-resolvable O2 (edgeless-$E_c^{\mathrm{sig}}|_C$). That is the *correct* reading and exactly the easy-family mass the $E_c$-as-definition silently dropped. The number that looked like a principled tightening of the *family count* is re-interpreted as a *within-family O2 sharpening* (fewer confusable pairs), leaving O1 unchanged.

**Net:** keep "one significance **schema** (over-similarity above a noise floor)," but replace "one criterion for O1+O2" with **"two nulls — transcript-homology-background for O1 ($E_r$), read-error for O2 ($E_c^{\mathrm{sig}}$/SUN/gate)."**

### 7.3 Implementation implication

The **shipped family catalog should be the homology / `denovo_families` grouping ($E_r$), not the read-conflict graph.** `src/rustle/vg_family/read_conflict.rs` is **unchanged** — it correctly stays the **O2** carrier (the de-tie / sig-tie graph, `MCC=χ(H)`, the exact-decomposition unit). The family-*definition* output should be the $\gamma$-quasi-clique-refined `denovo_families` catalog (with $R$ applied to dissolve DNFAM0), and the read-conflict graph should be consumed **within** each such family as the O2 resolvability structure. Concretely: O1 emits the ~1,130 (post-$R$) transcript-homology families; O2 runs `read_conflict` + the SUN ladder + the assign-or-abstain gate *inside* each, reporting $\chi(H)\le\Lambda$ and abstaining on the Tier-3 / K-frontier ~18%.

### 7.4 The copy *number* re-enters O1: $\chi(H)$ counts copies, including the unassignable ones (count vs. assign)

The demotion in §5–§7.3 moved the whole read-conflict *object* to O2. That is right for **assignment** but goes one
step too far for **counting**: it discards the family-level fact that the conflict graph also tells us **how many
copies a family has** — a quantity that is well-defined **even for copies the reads cannot assign**. We therefore
**re-promote the count** (not the assignment) to O1, and enrich the O1 family object accordingly. This is not a
reversal of the reframe: membership (*which* loci) stays the $E_r$ transcript-homology component; we add the
conflict-derived **multiplicity** as a family-level property, sitting beside $\Lambda(C)$.

The distinction is theorem-level (`THEORY.md` §3 Remark, Lemma 1 / Theorem 2 / Theorem 4):

- **Counting** copies $=\mathrm{MCC}=\chi(H)$ (chromatic number of the copy-conflict graph $H$). Needs **only** the
  conflict structure; defined even when **no read can be assigned**. On the **de-tied / significance-gated
  (copy-consensus)** graph and in the **full-span** regime, $\chi(H)$ is a **lower bound** on the true copy number
  (on the *raw* allele-disagreement read graph the error edges inflate the colouring to $\approx3\times K$, so the
  gate / copy-consensus is a *precondition* for the bound — `THEORY.md` §3 Lemma).
- **Assigning** is strictly harder, and splits into two theorems the family-level shorthand must keep **distinct**:
  **(i) single-read** assignability of a read $r$ ($\lvert N(r)\rvert=1$) needs a single-position Strong-Separation
  witness — a **SUN** — which is **Theorem 4(ii) / the Bridge** (the per-read certificate, and exactly the Tier-1
  property, a *per-copy* fact); **(ii) family-wide unique-cover** recovery needs **Strong Separation** (Theorem 2), a
  **coverage-dependent** condition. The **SUN / Tier-1** core is therefore **not** the same object as Strong
  Separation (an all-Tier-1 family with short reads can still fail Strong Separation) — the per-read assign
  certificate is **Theorem 4**, and Theorem 2 is reserved for the family-wide **cover** uniqueness.

So a family is **counted** as $\chi(H)$ copies while only the Tier-1 subset is single-read **assignable** — the
Tier-2 copies (distinct hap-vector, no single-position SUN witness) are **counted but not assignable**, and they
**matter**: they are real copies the family has, that the per-read gate must abstain on. The **enriched O1 family
object** is therefore the triple

$$
\text{family}\ =\ \big(\ \underbrace{C=R(\kappa)\subseteq V_R}_{\text{which loci: }E_r\text{ homology, }|C|=\Lambda(C)}\ ,\ \ \underbrace{K=\chi(H)}_{\text{how many copies (conflict-derived, incl. unassignable)}}\ ,\ \ \underbrace{(n_{\text{res}},\,n_{\text{c.u.}},\,\text{coll})}_{\text{3-tier hierarchy}}\ \big),
$$

with the **SUN 3-tier ladder** stratifying the copies:

| tier | definition | counted by $\chi(H)$? | single-read assignable? | who counts it |
|---|---|:---:|:---:|---|
| **Tier-1** RESOLVABLE | SUN-identifiable / Strong-Sep | yes (own color) | **yes** | $\chi(H)$ — the nice O2 core |
| **Tier-2** DISTINGUISHABLE-BUT-UNASSIGNABLE | unique hap-vector, **no** SUN | yes (own color) | **no** (needs $\ge2$-PSV co-observation) | $\chi(H)$ — the advisor's "important" copies |
| **Tier-3** IDENTICAL / collapsed | share a hap-vector, no distinguishing PSV | **no** ($\chi(H)$ collapses them) | no | read-**depth** / DNA (K=0 frontier, parCN / O4) |

**Exact per-family metrics** (`bench/copy_number_catalog.py` → `bench/copy_number_catalog.tsv/.json`; all identities
machine-checked, all 154 families):

- $n_{\text{loci}}$ — distinct co-located reference copies fed to the conflict graph $=\text{Tier1}+\text{Tier2}+\text{Tier3}$.
- $\chi(H)=K$ — number of distinct copy hap-vectors (colors) $=\text{Tier1}+\text{Tier2}+(\#\text{non-singleton groups})$.
- $n_{\text{resolvable}}=\text{Tier1}$ — counted **and** single-read assignable (Strong-Sep, the O2 core).
- $n_{\text{counted\_unassignable}}=\chi(H)-n_{\text{resolvable}}=\text{Tier2}+(\#\text{non-singleton groups})$ — copies $\chi(H)$ counts but a single read cannot pin (Tier-2 + one representative per collapsed group).
- $\text{collapsed\_excess}=\sum_{\text{non-singleton groups}}(\text{size}-1)=\text{Tier3}-(\#\text{non-singleton groups})$ — **identical copies $\chi(H)$ misses** (need depth/DNA).
- $\text{true\_copy\_lower\_bound}=\chi(H)+\text{collapsed\_excess}=n_{\text{loci}}$ (reference-resolved regime).

**Honest bound chain** (the crux; note the direction — $\chi(H)$ *under*-counts here because minimap2 pre-splits
distinct reference loci, so identical loci collapse):

$$
\underbrace{n_{\text{resolvable}}}_{\text{Tier-1}}\ \le\ \underbrace{\chi(H)}_{\text{conflict count}}\ \le\ \underbrace{n_{\text{loci}}=\text{true\_copy\_lower\_bound}}_{\chi(H)+\text{collapsed\_excess}}\ \le\ \text{true copy number},
$$

last gap $=$ reference-**absent** copies (O4). The advisor's dual statement $n_{\text{loci}}\le\chi(H)$ is the
complementary **reference-collapsed** regime (one reference locus hiding several hap-vectors — the Vollger/SDA
collapsed-segdup case; `reference_sda_vollger`), where the reads reveal *more* copies than the reference shows; the
unified statement is $\max(n_{\text{loci}},\chi(H))\le\text{true copy number}$, and the O1 object carries **both**
numbers plus $\text{true\_copy\_lower\_bound}$.

**Genome-wide (GGO, 154 co-located families / 412 copies):** the **read-observed** conflict count (`psv_graph_genomewide.json`, the directly-plumbed $\chi(H)$ of the complete-multipartite Lemma, `THEORY.md` §3) is $\sum\chi(H)=354=322$ **read-level** singleton parts $+\,32$ collapsed parts, so $354+58=412=\sum n_{\text{loci}}$ and **58** identical (Tier-3) copies go uncounted across **30/154** strictly-under-counting families (the read-level singleton count **322** is a read-sampling quantity and is **not** the tier count: the SUN tiers are defined on the copy-consensus graph, where Tier-1 $\uplus$ Tier-2 $=339$, differing from 322 by 17 over the 11 divergent families); the depth-independent **copy-consensus** count (`copyonly_K`, authoritative for the O1 count — see *Which conflict graph* below) is $\sum\chi(H)=361=338\ (\text{Tier-1})+1\ (\text{Tier-2})+22\ (\text{collapsed-group reps})$, missing **51** Tier-3 copies. Both obey the chain against $\sum n_{\text{loci}}=\sum\text{true\_copy\_lower\_bound}=412$: the identical copies counted only by reference-locus multiplicity, not by the conflict count, are the exact Tier-3 mass that needs read-depth/DNA. Exemplars: **fam 0** (GSTM2) $n_{\text{loci}}=\chi(H)=7$, all Tier-1 (fully counted **and** assignable); **fam 1** $n_{\text{loci}}=7,\ \chi(H)=1$ (all identical, collapsed\_excess $=6$ — $\chi(H)$ under-counts $7\!\to\!1$); **fam 42** $n_{\text{loci}}=\chi(H)=8$ but $n_{\text{resolvable}}=7$, $n_{\text{counted\_unassignable}}=1$ — the single Tier-2 copy the family *has* and *counts* but no single read can assign (the advisor's case in one row).

**Which conflict graph defines the O1 count.** $\chi(H)$ is taken on the **copy-consensus** conflict graph
(`sun_identifiability` `copyonly_K`): it counts copies distinguishable **in principle** from the assembled copy
sequences, independent of read-sampling depth, and is **tier-consistent** (singleton groups $=$ Tier1+Tier2). The
**read-level** graph (`psv_graph_genomewide.json` `K`) is the **assignment-realized** count — what *this* read set
resolves — and is carried as a cross-check. They agree on **143/154** families and diverge on **11**, in both
directions (the count-vs-assign split made concrete at the graph level): **fam 22** the read-level graph *over*-splits
($K=4$ vs copyonly $1$; read-level PSVs from within-copy noise not fixed between copy consensuses), and **10**
families (IKBKG/30, 46, 68, 94, 96, 101, 145, 151, 168, 195) the copy-consensus splits *more* (copyonly $2$ vs
read-level $1$; the copies carry consensus SUN but the reads never co-observed the distinguishing column). For the
O1 **count** the copy-consensus graph is authoritative; the read-level divergence is an O2 read-sampling fact.
**Reconciliation:** the read-observed $\sum\chi(H)=354$ and copy-consensus $\sum\chi(H)=361$ differ by exactly $+7$
over these 11 families (10 consensus-over-splits at $+1$ each, minus fam 22's read-over-split of $-3$); both
instantiate the complete-multipartite Lemma of `THEORY.md` §3 ($\chi(H)=n_{\text{groups}}$, and the
count-needs-only-conflict Proposition beside it), so the choice between them is a graph-substrate choice, not a
change to the $O1$ invariant.

---

## Appendix A — verification log (2026-06-30)

| check | source / command | result |
|---|---|---|
| $E_r$ catalog size | `denovo_families.tsv` (multi-copy) | **1,130 families / 3,636 loci**; excl DNFAM0 blob **1,129 / 2,908**; evaluable (fetchable coords) **1,126 / 2,902** |
| DNFAM0 blob | `denovo_families.tsv` row 2 | **728 members** chained chr1..chrY — RNA analog of the DNA 1,547-gene blob; dissolved by $R$ |
| **O1-DROP: $E_c$ de-tie** (headline) | `bench/rna_reframe_measure_detie.py` (reproduces `read_conflict.rs::de_tied`, $\Delta$0.005 / $\mathrm{DE_{max}}$0.05 / min 3, primary+secondary `de:f`; rows $\to$ `bench/rna_reframe_detie_rows.tsv`) | **335/1,126 = 29.8% families**, **741/2,902 = 25.5% copies** dropped (**no within-family de-tie edge** — the literal O1-drop) |
| O1-DROP: $E_c^{\mathrm{sig}}$ | `bench/rna_reframe_measure_sigtie.py` ($m_x=de_x\!\cdot\!\mathrm{aln\_len}_x$, tied iff $\varepsilon^{|m_i-m_j|}\ge\alpha$, $\varepsilon$1e-3/$\alpha$1e-3) | **370/1,126 = 32.9% families**, **828/2,902 = 28.5% copies** ($E_c^{\mathrm{sig}}\subseteq E_c$, drops more) |
| size-grading | `bench/rna_reframe_measure_detie.py` | size-2 284/855 = **33%**; size-3 **24%**; size-4/5 ~15%; size $\ge7$ **0%** |
| dropped-class content | symbol roots of dropped families | ZNF/ZBTB/TRIM/KRAB-ZNF/TMEM/PNMA/SLC4A; **266/335** contain protein_coding; **48** unnamed/novel |
| O2 (NOT O1-drop): MAPQ hard-tie | `bench/rna_reframe_measure_ec_er.py` on `GGO_mm.bam` (a $\ge5\%$-MAPQ0 **per-locus** proxy; its `EC_DROPPED`=96.6% is O2 hard-ambiguity, **overstates O1 ~3×** — use the de-tie row above) | $E_r$ loci **0.51% MAPQ=0** (907,096 reads); SUN co-located **7.4%** (72,667; **14×**); genome-wide **0.04%** (652,221) $\Rightarrow$ **99.96% MAPQ$>0$** — per-copy resolvability |
| $E_c\subseteq E_r$ containment | `bench/rna_reframe_validate.py` §B–§C (154 SUN families $\to$ denovo $E_r$ by best overlap; reads `denovo_family_edges.tsv` + faithful de_tied probe) | **129/132 = 97.7%** all co-located copies in ONE $E_r$ family (operational core_recip$\ge0.13$); **3** split = **operational shared-exon leaks** (0 cross-family core_recip edges; de-tie reads **20 / 121 / 65**; below the 0.13 floor; **not** $R$-fixable; absorbed by permissive $E_r^{\mathrm{asym}}$ — theorem stands); **0** EDGE_LINKED |
| SUN tiers (O2 per-copy) | `sun_identifiability.tsv` | 412 copies: Tier-1 **338 (82.0%)** / Tier-2 **1** / Tier-3 **73 (17.7%)** — per-copy resolvability |
| fully-Tier-1 families (O2, NOT O1-drop) | `sun_identifiability.tsv` + `bench/rna_reframe_validate.py` §D | **132/154 = 85.7%** fully-Tier-1 = **per-copy resolvable**; measured **132/132 = 100%** still form a de-tie edge ($E_c$-VISIBLE, **SURVIVE** $E_c$; 0 vanish) $\Rightarrow$ this is **not** an O1-drop |
| co-located resolvability (O2) | `psv_graph_genomewide.json` (154 families) | **82.5% per-copy-resolvable / 17.5% frontier(K=0)** |
| shipped oracle | `denovo_families.py:T_CORE` | core_recip $\ge$ **0.13** (canonical + RC POA retry) — the hard-floored $E_r$ |
| read-conflict predicate | `read_conflict.rs:86–106` | `de_tied`/`sig_tied` use only `de`, `aln_len`; **no partner-exon check** $\Rightarrow$ $E_r$ must be asymmetric for $E_c\subseteq E_r$ |
| operator $R$ | `genome_family_def.refine_families` | shared $\gamma$-quasi-clique ($\gamma=0.20$) + $\Lambda\ge2$; imported by DNA and protein notes |
| **O1 copy number** $K=\chi(H)$ (§7.4) | `bench/copy_number_catalog.py` (all-tier columns) + `bench/family_copy_number.py` (hap-vector verification companion) → `copy_number_catalog.tsv/.json`, `family_copy_number.tsv/.json` (all 154 fams, identities machine-checked, byte-deterministic, numerically identical on shared columns) | read-observed (`psv_graph`) $\sum\chi(H)=\mathbf{354}=322$ **read-level** singleton parts $+32$ collapsed parts, **58** Tier-3 uncounted, **30/154** strict; copy-consensus (`copyonly_K`, O1-authoritative) $\sum\chi(H)=\mathbf{361}=338$ Tier-1 $+1$ Tier-2 $+22$ collapsed-group reps, **51** uncounted (read-level singleton **322** $\ne$ copy-consensus Tier-1$\uplus$Tier-2 $=339$; tiers are defined on copy-consensus); $\sum n_{\text{loci}}=\sum\text{true\_copy\_lower\_bound}=\mathbf{412}$; reference-resolved chains $338\le354\le412$ and $338\le361\le412$; **22** families carry counted-but-unassignable copies, **21** need depth/DNA; the two substrates differ by $+7$ over 11 divergent families |
| copy count: consensus vs read-level | `sun_identifiability` `copyonly_K` vs `psv_graph_genomewide.json` `K` | agree **143/154**; diverge **11** (fam 22 read over-splits $4$ vs $1$; 10 fams consensus over-splits $2$ vs $1$) — count-vs-assign at the graph level; O1 count uses copy-consensus |

**Definitional one-liner for the thesis body.** An *RNA-level multi-copy gene family* is a **cohesive community of the transcript-homology graph** $G_R=(V_R,E_r)$ meeting the **same cohesion certificate** as the DNA/protein levels: a block $C=R(\kappa)$ that is a $\gamma$-quasi-clique ($\rho_{\mathrm{in}}\ge\gamma=0.20$) with $\ge2$ distinct expressed loci, where $E_r$ = significant spliced-transcript self-alignment ($\alpha_r$ / core_recip $\ge0.13$). It is the **fourth homology oracle** ($E_a$—$E_b$—$E_r$—$E_p$), **includes read-resolvable copies**, **carries a family-level copy number $K=\chi(H)$** (the conflict-derived count on the de-tied / copy-consensus graph — a lower bound on the true copy number that counts even the Tier-2 *unassignable* copies; reference-resolved chain $n_{\text{res}}\le\chi(H)\le n_{\text{loci}}=\text{true\_copy\_lower\_bound}\le\text{true}$, general form $\max(n_{\text{loci}},\chi(H))\le\text{true}$; §7.4, `bench/copy_number_catalog.py` / `bench/family_copy_number.py`), and **contains the read-conflict graph as its within-family O2 *assignment* substructure** ($E_c^{\mathrm{sig}}\subseteq E_c\subseteq E_r^{\mathrm{asym}}$, the sole clean containment, now the tightest since it shares $V_R$) — the **count** is O1, the **assignment** is O2. Read-ambiguity ($E_c$) and the SUN 3-tier ladder are **O2 resolvability for *assignment*, not the O1 family *boundary*** — but the ladder's per-tier **counts** ($n_{\text{res}}$, $n_{\text{c.u.}}$, collapsed\_excess) **are** O1 attributes of the enriched family object (§7.4). This is the correction that recovers the **~30% of multi-copy families / ~1/4 of copies** the prior $E_c$-as-definition silently dropped as easily-solvable.
