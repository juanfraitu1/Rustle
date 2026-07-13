# Five Concepts, Four Baselines: Paralog, Segmental Duplication, Multi-Copy Gene Family, Expansion, and Reference-Absent Copy

*A defense-grade formalization for the Rustle gorilla pan-transcriptomic pipeline. Each concept is defined against the reference frame it is actually measured against, mapped to our formal objects, and fenced with explicit non-implications. Every threshold and citation below survived an adversarial (Canzar-skeptic) pass; loose usages are disclosed, not hidden.*

---

## 0. Intro — these are NOT synonyms; they are four measurement axes

"Paralog," "segmental duplication," "multi-copy gene family," "expansion," and "reference-absent copy" are routinely used interchangeably in duplication genomics. They are not interchangeable. Each answers a different question and is measured against a **different baseline**. The five concepts collapse onto **four distinct measurement axes**, distinguished by their reference frame:

1. **Gene-phylogeny axis — PARALOG (vs ortholog).** Baseline = the reconciled gene tree. An *evolutionary relation on a pair of genes* (their last common ancestor is a duplication node). Not a set, not a count, not a structure.
2. **Genome-self-alignment axis — SEGMENTAL DUPLICATION (SD).** Baseline = one genome aligned to itself. A *DNA/structural block state* (a ≥1 kb, ≥90%-identical duplicated block within one assembly).
3. **Single-genome count axis — MULTI-COPY GENE FAMILY.** Baseline = one genome. A *within-genome membership + count state* (≥2 transcribed homologous copies at ≥2 loci) — **our primary object O1**.
4. **Comparative axis — the fourth axis carries the two concepts that require a baseline OUTSIDE the single genome, and they are orthogonal to each other:**
   - **EXPANSION.** Baseline = a **reconstructed ancestral / outgroup copy count**. A *cross-lineage, directional change* (copy gain along a branch). Needs ≥2 species + an outgroup to polarize.
   - **REFERENCE-ABSENT COPY (O4).** Baseline = a **specific reference assembly**. A *sample-vs-assembly gap* within one individual (a copy the reference does not host).

The organizing claim of this document is the **baseline lattice** (§6): the concepts differ primarily in what they are measured against, and most confusion in the field (and every reviewer trap for this thesis) comes from silently swapping one baseline for another. The student's specific question — *how is an expansion different from a reference-absent copy?* — is the case where the two comparative baselines (ancestor vs assembly) are most often conflated; §7 makes that distinction airtight.

**Oracle vocabulary used throughout.** Our family skeleton (nodes = loci, edge = a homology tie, family = a cohesion-refined connected component with ≥2 loci) is shared across **four homology edge oracles** that differ only in the edge test:

| oracle | edge test |
|---|---|
| **E_a** | genomic DNA / segmental-duplication homology (SEDEF) |
| **E_b** | exonic homology |
| **E_r** | spliced-transcript homology (exon-sum nt identity + reciprocal coverage) — **the O1 oracle** |
| **E_p** | reciprocal whole-protein homology |

**E_c** (read-conflict / MAPQ-0 **ambiguity**) is a **fifth, distinct oracle** and is NOT a homology oracle: it measures identifiability-relevance, not descent. It was deliberately **demoted** from the family definition (O1) to the within-family copy-assignment structure (O2). Standardize on `{E_a, E_b, E_r, E_p}` homology **plus** E_c ambiguity; the oracle set does not move.

---

## 1. Paralog (with the ortholog contrast)

### One line
Two homologous genes are **paralogs** iff their divergence traces to a gene-**duplication** event (their last common ancestor in the gene tree is a duplication node), as opposed to **orthologs** whose divergence traces to a **speciation** event. It is a **pairwise, evolutionary relation** of any age, location, or mechanism — read off the gene phylogeny, not off sequence identity.

### Formal definition

**Scope assumption (name it).** Fix the standard duplication/loss (DL) reconciliation with **no horizontal gene transfer / no xenology**. Under this assumption every internal node of the gene tree is speciation (σ) or duplication (δ) — a binary partition — and the within-genome-collapse theorem below holds. Under HGT a same-genome pair can be *xenologous* (σ-lca via a transfer node; Gray & Fitch 1983), which is out of scope for gorilla nuclear genes but is **stated, not silently assumed** (this is the load-bearing step, so the condition is named).

**Estimand vs estimator.** Paralogy is defined by the **true** event type at the **true** last common ancestor in the true gene genealogy (the estimand):
- **ORTHOLOG(x, y)** ⇔ the lca of x, y is a **speciation** event.
- **PARALOG(x, y)** ⇔ the lca of x, y is a **duplication** event.

The rooted **reconciled** gene tree *T* (leaves = extant genes; internal nodes labeled σ/δ by reconciliation against the species tree *S*) is our **estimator** of that genealogy. Reconciliation is error-prone — it fails precisely on **recent** paralogs, which is our target regime — so the σ/δ labels are inferred, not given. Homology is the ambient equivalence closure; ortholog/paralog **partition** the homolog pairs (complementary and mutually exclusive on homologs, under no-HGT).

**Structural properties.**
- It is a **relation on pairs**, not a property of a gene or a set. It is symmetric but **not transitive** and **not an equivalence relation** (one-to-many and many-to-many orthology are the norm after duplication).
- **Within-genome collapse.** Under no-HGT, any two distinct homologous genes co-resident in **one** genome necessarily have a δ lca (a σ node's children segregate into different species). Hence single-species homolog pairs are **all paralogs by construction**, and Fitch's actual difficulty — discriminating σ vs δ — **does not arise with one genome**.
- **Sub-typing** (Sonnhammer & Koonin 2002): relative to a chosen reference speciation node, *in-paralogs* (duplication after that speciation) vs *out-paralogs* (before it); *ohnologs* = whole-genome-duplication paralogs (Ohno 1970; term Wolfe 2000). Each refinement requires designating a reference σ node — i.e. ≥2 species — so all are undefined from a single genome.

**What a FAMILY adds.** A gene family is a **set**, not a pairwise relation: a chosen grouping of homologs (a clade of *T*, or a similarity cluster). Converting Fitch's node-exact edge relation into a family needs an **extra decision** that paralogy does not: *where* to cut *T* (superfamily vs subfamily) or which clustering threshold to apply. A family therefore (i) collapses the directional, non-transitive, edge-labeled structure into flat membership, and (ii) is **not uniquely determined by the Fitch relations alone** (nested families at different levels). In one genome a family ≈ a monophyletic set of mutual paralogs; across species a family mixes orthologs and paralogs.

### Baseline
The **gene phylogeny** — specifically the true gene genealogy (estimand), approximated by a rooted gene tree reconciled against the species tree (estimator). Paralogy is a statement about the **type of one internal node** (the lca). This needs **neither a second species nor an assembly** — only tree topology and node labels.

### Map to our objects
- **E_r** is our operational proxy for the **within-genome paralog SET**. Because our substrate is a single genome (mGorGor1 gorilla), within-genome collapse applies: every homologous locus pair inside one E_r component is in principle a paralog pair, so an **E_r family (O1) targets a paralog group** — and we can define families **without a gene tree** precisely because the σ-vs-δ problem is moot in one genome. Only the grouping remains.
- **O1** recovers paralog-group *membership* but **not** the pairwise Fitch labels, the tree, ancestry, or duplication order/dates.
- **E_c** is **incomparable** with paralogy: paralogy ⊄ E_c (resolvable true paralogs like RFPL1/2/3 carry zero read-conflict) **and** E_c ⊄ paralogy (E_c fires via shared repeats / exonized TEs that are not duplication edges). E_c measures identifiability, not descent. (Our documented within-stack containment is E_c ⊆ E_b — a de-tie read witnesses local/exonic homology — which is about E_b, not about the true paralogy relation.)
- The four homology oracles E_a/E_b/E_r/E_p are **pairwise incomparable** (each has a witnessed set-difference both ways; `bench/family_definition_formal.md §6`). Informally their **recall** of true paralogy differs — **E_p uniquely reaches twilight-zone paralogy** (the E_p \ E_r globin case; documented witness **MB ↔ HBB at ~24% aa identity**), E_r misses ancient divergent paralogs (blind spot), E_c only touches the small co-located near-identical core — but this is a **recall heuristic, not a containment order**.
- **χ(H) = MCC** (per-genome copy number) counts copies that are mutually paralogous, but it is a **count, not the relation**; the copies it enumerates are paralogs even at the K=0 frontier where they are sequence-identical.
- **O2** assigns reads to paralogous copies; it presupposes the paralog set but never labels edges σ/δ. **O4** copies are additional paralogs of the in-reference copies. **SUN** (Sudmant 2010, SUN ⊊ PSV) is a single column where one paralog's base is private — the per-copy marker that distinguishes paralogs.
- The **ortholog** contrast is what the expansion leg needs: expansion compares copy counts of **orthologous** families across lineages. Paralog = within-genome multiplicity; ortholog = the cross-genome correspondence along which that multiplicity is compared.

### Non-implications
Paralogy does **not** imply: (1) **co-location** (paralogs may be tandem, dispersed, or trans-chromosomal). (2) **high sequence identity** (globins 24–40% aa id are paralogs; conversely high identity is not sufficient — near-identical sequences may be **alleles**, the Bailey PSV-vs-allele problem). Identity is not the definition; tree topology is. (3) **same function** (neo-/sub-functionalization is routine; the "ortholog conjecture," Gabaldón & Koonin 2013, is exactly the open debate that function-conservation is not entailed). (4) **within-genome-only** (out-paralogs live in different species; the "all same-genome homologs are paralogs" implication is one-directional). (5) **resolvability/observability** (K=0 identical unresolved copies are still paralogs). (6) **"paralog" = "copy"** (a copy is a within-individual count; a paralog is an evolutionary edge; two copies are one paralog *pair*).

### Our-usage caveats (load-bearing)
E_r is **not** a gene tree and we do **not** infer duplication nodes: we never reconcile against a species tree, never label σ/δ, never root/order/date duplications, never assign ancestry. E_r is a similarity/homology graph thresholded by exon-sum identity + reciprocal coverage and refined by a γ-quasi-clique operator — a **noisy proxy for the paralog SET**, not the Fitch relation. Consequences we disclose: (i) **E_r ⊋ paralogy** in places (shared repeats / exonized TEs / domain bridges link non-paralogous loci — the AMY2A–ZNF91 over-merge class), so an E_r edge is not a certified duplication edge; (ii) **E_r ⊊ paralogy** in places (divergent paralogs whose reads don't cross-map are missed — the recall gap; ancient divergent paralogs are our blind spot and belong to E_p); (iii) we legitimately **avoid** the σ-vs-δ problem **only** because we are single-species — a scope statement, not a solution: **the moment the thesis makes any cross-species (expansion/orthology) claim, a gene tree or an explicit orthology map becomes required and E_r alone is insufficient**; (iv) "family" (O1) is a threshold/clustering decision layered on E_r, inheriting both arbitrary-level ambiguity and Canzar's standing suspicion of similarity thresholds — defended as identifiability-scoped grouping, **not** as reconstructed paralogy.

### Literature
Fitch 1970, *Syst. Zool.* 19(2):99–113 (canonical ortholog/paralog); Fitch 2000, *Trends Genet.* 16(5):227–231; Sonnhammer & Koonin 2002, *Trends Genet.* 18(12):619–620 (in-/out-paralog); Koonin 2005, *Annu. Rev. Genet.* 39:309–338; Gabaldón & Koonin 2013, *Nat. Rev. Genet.* 14:360–366 (ortholog conjecture); Gray & Fitch 1983 (xenolog); Ohno 1970 / Wolfe 2000 (ohnolog). For the family/clustering contrast: Vilella et al. 2009 (Ensembl Compara/TreeBeST reconciliation, emits Fitch labels) and Emms & Kelly 2019 (OrthoFinder orthogroups) — both noted to **fail on recent near-identical paralogs**, the exact regime our long-read method targets.

---

## 2. Segmental Duplication (SD)

### One line
A duplicated genomic **block** — a contiguous stretch **≥1 kb** present in **≥2 near-identical (≥90%) copies within a single genome** — a DNA/structural object carrying flanking context, detected by genome self-alignment (WGAC) confirmed by read depth (WSSD), and the substrate on which paralogous sequence variants (PSVs) live.

### Formal definition
Fix one (ideally haploid) assembled genome *G*. **Bailey-2002 predicate:** a pair of intervals (A, B), A ≠ B, forms a segmental-duplication edge **SD(A, B)** iff there exists a maximal local self-alignment *a*: A ↔ B with (i) length ≥ 1 kb, (ii) identity ≥ 90% (Bailey's raw sets ~91–100%), (iii) copy multiplicity ≥ 2 in *G* (implied by the pair predicate), and (iv) neither interval is explained by a high-copy interspersed-repeat/TE family (genomic copy number below a **TE cutoff R** — the sole irreducible detection knob; TEs are RepeatMasked before and re-filtered after). Any column difference between the two copies is a **PSV**, categorically distinct from an allelic SNP.

**Detection = two-method intersection SD = WGAC ∩ WSSD.** WGAC = **BLAST-based** assembly self-alignment (keep pairs ≥1 kb, ≥90%). WSSD = WGS read depth in ~5 kb windows, >3 SD above the unique-region mean (depth tracks copy number at R² = 0.96). WGAC-alone flagged 16.5% of the human genome, but ~4/5 of the >98% alignments were allelic/mis-assembly artifacts; WSSD intersection dropped it to **~5% of the euchromatic genome (130.5 Mb)**. **SEDEF** (Numanagić et al. 2018) / **BISER** are the modern **minimizer-era** WGAC (seed → chain → align on one assembly, read-independent) — WGAC only.

**In our pipeline, E_a = the SD edge oracle over loci:** an edge (gene *u*, gene *v*) iff one SEDEF pair (`final.bed`) covers side-A of *u* and side-B of *v* each ≥50% (COVER = 0.5). `final.bed` cols: 12 = aln_len, 21 = fracMatch (gap-excluded identity), 23 = Jukes–Cantor divergence. Our `final.bed` is a **SEDEF candidate SUPERSET** at SEDEF's default `--max-error 0.5` (~50% identity / 1 kb) floor, **NOT** the textbook ≥90% set: raw n = 253,029 pairs (median identity 0.870, 72.4% below 90%); applying the Bailey predicate (fracMatch ≥ 0.90 & aln_len ≥ 1000) survives 66,142 (26.1%); + TE-exclusion survives 27,623 (10.9%). **Soto's SD98 = the identity ≥ 0.98 sub-band of this object.**

### Baseline
The **genome self-alignment** — a locus is an SD iff it is duplicated **elsewhere in the same single genome/haplotype**. Reference frame is **within one assembly** (intra-genomic, one individual). Two implicit sub-baselines: (1) the ≥90% identity criterion is a **recency** baseline against the neutral primate substitution rate (~90% id ≈ last ~40 My), which excludes ancient whole-genome-duplication ohnologs (2R, ~450–500 Mya); (2) the >3-SD depth cut is against the unique-region read-depth mean. **Not** a cross-species/outgroup baseline (that is expansion) and **not** a sample-vs-reference baseline (that is O4). On a haploid T2T the allele-vs-paralog ambiguity (2002's deepest problem, which needed WSSD depth: ~1× = allele, ~2× = duplication) is resolved **by construction** — no second haplotype exists for an allelic overlap, so every self-alignment column difference is a PSV, **modulo assembly correctness**.

### Map to our objects
E_a is exactly our DNA/segdup edge oracle, used only as an **external, corroborating** oracle — **our family (O1) is an E_r component, not an E_a component.** Lattice relations (on the Bailey-thresholded predicate): the only clean containment is **E_c ⊆ E_b** (read-ambiguity ⊆ asymmetric exonic homology; E_b is the exonic-homology oracle, the same transcribed/exonic axis as E_r — stated here to avoid an undefined symbol). **E_a is an incomparable third circle:** E_c \ E_a witnessed by **APOBEC3D/F** (in E_c but SEDEF identity 88.4% < 90% ⇒ not an SD); E_a \ E_c = the large unexpressed/intergenic SD mass; **RABL2A/B** sits in all three (verified triple-core 25 kb / 98.6% SD). So the naïve "read-conflict ⊆ segdup" is **FALSE**. **χ(H) = MCC lives over E_c, not over E_a**; the K-frontier is the high-correlation limit of the SD triple-core. SD is the substrate of the per-copy markers: **Bailey coined PSV**; **SUN** (Sudmant 2010, SUN ⊊ PSV) = a single-position private allele of one SD copy = the Strong-Separation single-read witness and the parCN marker.

**Cross-validation against SD98 is deliberately ONE-SIDED (and passed):** **DAZ (99.63%), RBMY (99.74%), TSPY (99.45%)** are literally SD98 segdups Soto's pipeline clusters and we recover them. **GSTM (95.29%)** and **MAGEA (94.85%)** are SD but **below SD98**. **PCDHB (88.86%)** is **below the 90% Bailey SD floor entirely — NOT a segmental duplication at all (∉ E_a)** — yet a real family we resolve by the PSV gate: **PCDHB is the single strongest witness that RNA + PSVs reach past the DNA identity cutoff** (do not lump it with the below-SD98 cases). This asymmetry is one-sided: SD / SD98 membership **confirms** a family; non-membership does **not** refute one — never invert it into a specificity claim.

SD is orthogonal to O4: an SD is **in** the assembly by definition, whereas O4 is what is **missing/collapsed from** it — but a **collapsed SD** (mis-assembled merge, invisible to depth-free SEDEF-on-T2T) is exactly where O4 lives.

### Non-implications
SD does **not** mean gene family (a block may contain a whole gene, a partial gene, or no gene). SD does **not** mean expansion (within-genome state, no polarization). SD does **not** mean reference-absent (SD is present in the assembly; O4 is absent from it). High SD identity does **not** mean allele (on a haploid T2T, self-alignment differences are PSVs by construction). Being an SD is **not necessary** to be a real family (ancient paralogs below ~90%/98% leave no SD signature — globins, and GSTM/MAGEA below SD98, PCDHB below the SD floor). SD is **not** read-conflict (E_a incomparable to E_c). `final.bed` is **not** the Bailey SD set as-shipped (it is a ~50%-identity candidate superset until the ≥1 kb / ≥90% / not-TE predicate is applied).

### Our-usage caveats
(1) **Loose load:** `genome_family_def.py::load_sedef_pairs` keeps only the 6 coordinate columns and applies **no identity/length/TE filter**, so E_a is currently built on all 253,029 pairs at the ~50% floor — this manufactures repeat-array over-merges (TRNAV-CAC ×173, an rRNA cluster ×70, DNFAM0 with 728 members) that are **removable artifacts of not applying SD(·)**, not properties of segdups. The thesis must state the predicate (≥1 kb, ≥90%, not-TE) as **part of the definition**, or use the opt-in `--bailey-sedef` filter (which, honestly, also drops real divergent families CEACAM/KRAB-ZNF/PRSS/IFITM/ULBP off the 90% cliff). (2) SD is only a corroborating oracle; "read-conflict families ⊆ segdups" is FALSE (E_a ⟂ E_c). (3) The SD98 check is one-sided (confirm-only). (4) We **cannot run WSSD/famCN/parCN** — no gorilla WGS depth — so our SEDEF object is **WGAC-only** (self-alignment structure, read-independent), carrying no per-individual copy-number term. (5) Two depth-only artifacts survive on SEDEF-on-T2T: a collapsed segdup is invisible (where O4 and copy-vs-allele stay open, needing DNA parCN or the mGorGor1 second haplotype), and a false-duplication/haplotype-switch mis-assembly can manufacture a spurious self-pair. (6) Overlaying SEDEF against RNA families is **not independent**: our `denovo_transcripts.fa` are reference substrings spliced at read-defined coords, and SEDEF self-aligns the same reference bytes — shared-reference-substrate corroboration, not read-independent (genuine independence needs read-consensus PSVs (A1/SDA) or the DNA second haplotype).

### Literature
Bailey et al. 2002, *Science* 297:1003–1007 (coined SD ≥1 kb/≥90%, WGAC ∩ WSSD, and PSV). Numanagić I, Gokkaya AS, Zhang L, Berger B, Alkan C, **Hach F** 2018, *Bioinformatics* 34:i706–i714 (SEDEF, modern minimizer WGAC). Soto/Dennis et al. 2025, *Cell* 188:5363–5383, DOI 10.1016/j.cell.2025.06.037 (SD98 ≥98%; famCN via WSSD, parCN via QuicK-mer2). Sudmant et al. 2010, *Science* 330:641 (SUN). Fitch 1970, *Syst. Zool.* 19:99 (the evolutionary vocabulary SD instances).

---

## 3. Multi-copy gene family (RNA-level; our primary object O1)

### One line
A descriptive **state of a single genome** — two or more transcribed, homologous gene copies (paralogs) at **≥2 distinct loci** — formalized as a γ-quasi-clique-refined connected component of the transcript-homology graph, with per-genome copy number reported as **two independent censored lower bounds, Λ(C) and χ(H) = MCC**.

### Formal definition
Let *V_R* = the set of de-novo **expressed** loci (representative spliced transcripts; ≥ MIN_READS = 3 reads, all-canonical junctions). Let *G_R* = (*V_R*, E_r) be the transcript-homology graph, where {u, v} ∈ E_r iff a significant local spliced-transcript self-alignment exists between *t_u* and *t_v*.

**Definitional predicate (own the computed threshold).** Operationally, **{u, v} ∈ E_r iff POA contiguous-core reciprocal coverage ≥ γ_core = 0.13** — a named, disclosed threshold. The Karlin–Altschul form *E(a) = K·|T|²·e^{−λS(a)} < α_r* is an **interpretive analogy that is evaluated nowhere**; its i.i.d./ungapped assumptions do **not** hold for spliced transcripts, and |T|² is the all-vs-all search space, not the per-pair one. We present the 0.13 coverage floor as the definitional predicate and do not dress it as significance.

A **multi-copy gene family** is a vertex set *C* ⊆ *V_R* with *C* ∈ *R*(κ) for the connected component κ containing *C*, where *R* is the cohesion-refinement operator: *C* is a **γ-quasi-clique** (internal edge density ρ_in(C) ≥ γ = 0.20, or |C| ≤ 2) **and** Λ(C) ≥ 2 distinct loci (reciprocal-≥50%-overlap locus collapse). The **cohesion certificate** — "γ-quasi-clique ∧ Λ ≥ 2 inside a raw E_r component" — not any one partition, is the canonical, seed-fixed, permutation-invariant object (exact max-γ-quasi-clique partition is NP-hard, so *R* emits one byte-identical witness). *(Note: the |C| ≤ 2 escape means the density test does no work for size-2 families — the Λ(C) ≥ 2-loci predicate is the sole gate there, so the cohesion certificate is not overclaimed for the most common size-2/3 families.)*

**Per-genome copy number = two independent censored lower bounds; neither dominates:**
- **Λ(C)** = the number of distinct expressed loci.
- **χ(H) = MCC** = the chromatic number of the **read-conflict graph H** (vertices = **reads** within *C*; a copy = an independent set / color; **MCC = minimum clique cover of the compatibility complement H̄**, Lemma 1) = the minimum number of copies that explains the ambiguous reads.

The relation between them is **not** a total order:
- **χ(H) < Λ(C)** whenever copies are **uniquely-mapping-divergent** (H is edgeless, no reads conflict them — **RBMY: 6 loci but χ(H) = 1**) **or** reference-resolved-but-transcript-identical (**K=0 frontier**, collapse to one color).
- **χ(H) > Λ(C)** whenever a **reference-collapsed / O4** copy piles reads at a single assembled locus that PSVs split into ≥2 colors — exactly the "O4 lifts χ(H) toward the true count" story.

The reported per-genome copy number is therefore **max(Λ(C), χ(H))**, and both are **≤ the true count** (transcribed-only, resolvable-only). This matches the memory lower bound **max(n_loci, χ(H)) ≤ true (338 ≤ 361 ≤ 412)**. **MCC = χ(H) holds on the de-tied / error-free read-conflict graph (Thm4 gate), NOT the raw allele-disagreement graph** (error-inflated; coloring median ~3× K).

The combinatorial **skeleton** (nodes = loci, edge = a homology tie, family = *R*-refined component with |C| ≥ 2 and Λ(C) ≥ 2) is shared across the **four homology oracles** E_a genomic / E_b exonic / E_r transcript / E_p protein; **E_c** (read-conflict / MAPQ-0 ambiguity) is a **fifth, distinct oracle demoted to the O2 within-family structure — not one of the four.**

### Baseline
**One genome** — a single-individual snapshot (the sequenced gorilla's testis transcriptome mapped against its own **mGorGor1** assembly; RNA donor = the mGorGor1 individual). No ancestor, no outgroup, no second species, no sample-vs-assembly comparison. A purely **within-genome homology + count state**, licensing **no** evolutionary or cross-species claim. (Contrast: expansion is measured against an ancestral/outgroup count and needs ≥2 species; O4 is measured against the reference assembly. Multi-copy-ness is orthogonal to both.)

### Map to our objects
This **IS O1**, our primary object; its oracle is **E_r**. **E_c is not the family definition:** it was demoted to the within-family O2 structure because as a definition it silently drops divergent uniquely-mapping families (the globin problem — MB ↔ HBB share no read). **E_c ⊆ E_b** is the sole clean containment (a read conflict implies **exon-level** shared sequence; it does **not** imply full-transcript E_r homology — the 0.13 floor's ~2.3% single-conserved-exon leak is exactly the E_c ⊄ E_r gap). **χ(H) = MCC** (Lemma 1) is an O2 quantity computed inside each fixed family; copy-assignment (O2) = max-weight facility location, assign-or-abstain, no 1/k. **SUN** (Sudmant 2010, SUN ⊊ PSV) = the per-copy single-read identifiability witness inside the family (3-tier ladder: SUN-taggable / hap-unique-only / K-frontier). **O4** is an orthogonal axis (sample-vs-assembly) that **interacts**: multi-copy families are the most collapse-prone, so O4 recovery (hidden_copy detector: 18 collapsed flags, 4 endogenous MHC copies) lifts χ(H) toward the true genomic count.

### Non-implications
(1) **Not** an expansion (multi-copy ≠ expanded; **RBMY = 6 copies in both gorilla and human** — multi-copy, not expanded). (2) **Not** a segmental duplication (SD is a proper subset; the family is broader — globins MB↔HBB ~24% aa id; CEACAM5/6/7, KRAB-ZNF, PRSS, IFITM, ULBP all <90% and dropped by an SD filter but captured by E_r/E_p). (3) **Not** reference-absent (O4 is a sample-vs-assembly gap; the family is a within-genome state — independent axes). (4) **Not** orthology (members are paralogs). (5) **Read-ambiguity (E_c) is not required for membership** (a divergent uniquely-mapping copy is a bona fide member; E_c-as-definition dropped ~30% of families / ~25% of copies, size-graded size-2 33% → size ≥7 0% — the exact flaw the E_r reframe corrects). (6) **χ(H) is not the true genomic copy number** — a functional lower bound.

### Our-usage caveats
**χ(H) = MCC is a functional lower bound**, loose three ways: (a) **transcribed-only** (silent/pseudogene/tissue-unexpressed copies invisible; single testis library); (b) **resolvable-only** (K=0 collapsed identical-transcript copies collapse to one color, χ(H) < Λ(C) at the K-frontier); (c) node set = **expressed** de-novo loci only. The DNA oracle confirms the direction: mGorGor1 diploid asm_hapCN ≥ K_read for 62/62 A1 families. The E_r operational oracle uses named arbitrary thresholds **γ_core = 0.13** and **γ = 0.20** — disclosed, not hidden (Canzar distrusts thresholdless graph-merging). Transcript sequences are reference substrings spliced at read-defined coords (L1 circularity: correct for **defining** families, but "read-independent" corroboration is not achievable from this substrate alone — RNA and SEDEF share the reference). Λ(C) locus collapse is ambiguous in tandem arrays (locus↔gene mapping brackets the count). **Precision/recall depends on the oracle** (state it beside each number): vs the **diploid-gold** oracle P ≈ 0.85 / R ≈ 0.84; vs the **Compara** oracle refined precision ≈ 0.94 / recall ≈ 0.52 (identity-graded) — the 0.84-vs-0.52 gap is different oracles, not a contradiction. Residual FPs (repeat-bridge hubs, GSTM2-domain, MAGE-cardinality) are DNA-boundary and RNA-invisible; recall is expression- and resolvability-bounded, not method failure.

### Literature
Fitch 1970, *Syst. Zool.* 19:99–113; Ohno 1970 (*Evolution by Gene Duplication*); Nei & Rooney 2005, *Annu. Rev. Genet.* 39:121–152 (canonical multigene-family review, birth-and-death); Bailey et al. 2002, *Science* 297:1003 (SD, PSV); Soto/Dennis 2025, *Cell* 188:5363 (SD98 + famCN, the DNA/CN analog); Sudmant et al. 2010, *Science* 330:641 (SUN). Our objects: `bench/family_definition_formal.md`, `bench/DNA_FAMILY_DEFINITION_FORMAL.md`, `bench/THEORY.md` (Lemma 1: MCC = χ(H)).

---

## 4. Gene family expansion

### One line
A comparative, **directional** claim that a family's copy number **increased along a specific lineage relative to its ancestor** — a cross-lineage change, not a within-genome state and not a sample-vs-assembly gap.

### Formal definition
**Objects.** (1) A gene family *g* — a set of homologous loci partitioned into orthologs/paralogs by Fitch (1970). (2) A rooted species phylogeny *T* = (*V*, *E*), leaves = extant species, **root fixed by an outgroup**. (3) Observed per-genome counts **n_g: leaves → ℕ** (integers) and **reconstructed ancestral counts n̂_g: V → ℝ_{≥0}** at internal nodes under a gene birth–death process — canonical estimator CAFE ML (Hahn et al. 2005; De Bie et al. 2006), or Dollo/parsimony as a coarse fallback. *(Typing matters: leaf counts are integers; ancestral reconstructions are real-valued expectations.)* (4) A **cross-species orthology map φ** that identifies the family unit across genomes — **required for any Δ to be well-typed** (see the correspondence gap below).

**Predicate (test form, not a bare inequality).** For a branch *b* = (u → v):

> **EXPANSION(g, b)** ⇔ a **significant** net copy gain along *b* under the birth–death null (CAFE branch *p* < α), i.e. n̂_g(v) > n̂_g(u) with the null rejected; **CONTRACTION** symmetric.

Magnitude = n̂_g(v) − n̂_g(u) (or log-ratio). A "lineage-specific expansion of *g* in *L*" = EXPANSION(g, b) for the branch subtending *L* (terminal = species-specific, internal = clade-specific). *(A bare n̂_g(v) > n̂_g(u) on real-valued ML expectations is satisfied by any ε > 0, so the significance test — an unavoidable, canonical threshold here — is part of the predicate, not bolted on.)* The **outgroup + rooted tree + ancestral reconstruction** is precisely what **polarizes** a raw between-taxon inequality n_g(A) ≠ n_g(B) into directed *gain-in-A vs loss-in-B*.

**Minimal operational form (2 taxa + outgroup, parsimony).** Taxa {A, B}, outgroup O, rooted 3-taxon tree; n̂_g(MRCA(A, B)) = the Fitch/Dollo parsimony reconstruction from {A, B, O}; **expansion on the A-branch iff that reconstruction < n_g(A)** (no approximate-equality — parsimony gives the ancestral value directly). Two taxa **without** an outgroup give only the unpolarized difference (necessary but not sufficient for an expansion call).

### Baseline
The **reconstructed ancestral copy number** n̂_g(u) at the parent node (equivalently, in a 2-taxon+outgroup contrast, the outgroup/ancestral count). Measured against a **point on the phylogeny**, one species-tree branch back in time. Explicitly **not** the reference assembly (that is O4) and **not** zero/absence. Requires **≥2 species + an outgroup** to instantiate and polarize.

### Map to our objects
The family *g* is our **O1** object (a within-genome E_r component). Our per-genome count is **χ(H) = MCC**. Our expansion statistic is **Δχ(H) = χ_A(H) − χ_B(H)** for the **φ-corresponded** family.

**The load-bearing dependency (disclose it):** χ(H) and E_r are **strictly within-genome** objects. Nothing in our oracles identifies the gorilla family unit with the human one — the correspondence **MAGEA_gorilla ↔ MAGEA_human** is **imported from gene annotation / Fitch orthology applied externally**, i.e. φ is currently **supplied by annotation, not inferred by us**. This is a **stronger** borrowing than just polarization: not only *which lineage changed* but the *orthogroup correspondence itself* is borrowed. Named, not hidden.

**Δχ(H) is a stack of two lower bounds:**

> Δχ(H) ≤ Δ(transcribed-resolvable copies) [χ(H) = MCC is a lower bound; K=0 exonically-identical copies collapse to Tied] ≤ Δ(famCN) [famCN is DNA read-depth and also counts silent pseudogenes our transcribed oracle never sees].

So the defensible statistic is precisely: **Δχ(H) is a conservative lower bound on Soto's famCN-based expansion magnitude.** **SUN** (Sudmant 2010) is the per-copy witness making expanded copies single-read-taggable (χ(H) T1 tier). **O4** is orthogonal to expansion but coupled operationally (§7).

**Demonstrated contrasts (borrowed polarization, disclosed):** MAGEA **2 (gorilla) → 11 (human)** and TSPY **5 → 33** (identical binary/recipe; `HUMAN_CROSSSPECIES.md`, `SOTO_CONCORDANCE.md`), against conserved control **RBMY 6 = 6**.

### Non-implications
(1) **Not** high copy number (a large family can be entirely ancestral — expansion is a *change on a branch*). (2) **Not** a reference-absent copy (cross-lineage-vs-ancestor, not within-species-vs-assembly). (3) **Not** inferable from a single genome, nor from a bare two-taxon inequality (polarization needs the outgroup + rooted tree + reconstruction). (4) **Not** an adaptive/functional claim (copy gain ≠ positive selection / expression gain / neofunctionalization — needs dN/dS, expression, dosage). (5) **Not** contraction (directional; mis-polarization flips them). (6) **Not** a statement about a specific paralog's identity (a count over the ortholog group; per-copy identity needs parCN/SUN). (7) **Not** identical to Soto's "duplication": the great-ape-max **2.5** cut classifies whether the human gain came from a **near-single-copy** ancestor (great-ape max < 2.5, "expansion") **vs an already-duplicated** ancestor (> 2.5, "duplication") — a split **our two-species contrast cannot make** because we have no ancestral estimate.

### Our-usage caveats
Our usage is a **two-species contrast, not a polarized phylogenetic expansion**, loose four ways. (i) χ(H) is a **lower bound**, expression- and resolvability-gated (only transcribed copies; K=0 identical copies collapse to Tied), so Δχ(H) is a lower bound on true magnitude and **not** the same statistic as Soto's DNA famCN (which counts silent pseudogenes). (ii) We have **only gorilla + human — no outgroup, no in-pipeline ancestral reconstruction, and φ from annotation** — so MAGEA 2 vs 11 and TSPY 5 vs 33 are **differences**; calling them "human-lineage expansions" **borrows polarization** from external phylogeny (Soto's ape-outgroup design), not from our data. To make the claim ours end-to-end we would need ≥1 outgroup (orangutan/macaque), a CAFE-style branch test, and φ computed by our own cross-genome E_r oracle. (iii) Cross-species comparison mixes tissues/libraries/depths (gorilla testis vs a single human testis library); low expression can **hide** an expansion as a false contraction (PCDHB 16 copies, DAZ 4 copies under-expressed and unresolved in the human library) — an ascertainment confound distinct from biology. (iv) Assembly collapse can **fake** (over-collapse the ancestor) or **hide** (under-count the derived) an expansion unless the true count is recovered — exactly why O4 is coupled to expansion (Soto's rationale). **Net defensible claim: "χ(H) differs across species, a conservative lower bound consistent with the annotated human-specific expansion (polarized by Soto's ape-outgroup design)" — NOT "we independently inferred and polarized a lineage-specific expansion."**

### Literature
Fitch 1970, *Syst. Zool.* 19:99–113; Bailey et al. 2002, *Science* 297:1003 (WSSD CN, famCN/parCN ancestor); Hahn et al. 2005, *Genome Res.* 15:1153 and De Bie et al. 2006, *Bioinformatics* 22:1269 (CAFE ML birth–death); Demuth & Hahn 2009, *BioEssays* 31:29 (expansion/contraction concept); Soto/Dennis 2025, *Cell* 188:5363–5383, DOI 10.1016/j.cell.2025.06.037 (human-specific expansions; call = median human famCN > max great-ape famCN across chimp Clint, bonobo Mhudiblu, **gorilla Kamilah**, orangutan Susie; duplication-vs-expansion split at great-ape max 2.5); Sudmant et al. 2010, *Science* 330:641 (SUN).

---

## 5. Reference-absent copy (our objective O4)

### One line
A gene-family copy **transcribed in the sequenced individual's genome but not hosted by the reference assembly** (collapsed, sequence-divergent-in-place, or wholly absent) — inferred from reads that **over-determine** the modeled reference copies, and **detected-and-flagged, never placed**.

### Formal definition
Fix a transcribed multi-copy family *F* (our O1). **Baseline = a specific assembly A** (here the haploid T2T gorilla, **mGorGor1**; the in-pipeline accession is GCF_029281585.2 — the exact accession should be confirmed against the assembly actually mapped against). To keep the comparison **ploidy-consistent**, count distinct paralog **loci at haploid resolution**, quotienting out allelic (het) variation:
- n_loci(*F*, A) = the number of **modeled** *F*-loci A hosts (the slots the copy-assignment model has). *(Defining n_ref via "resolvable loci" would beg the identifiability question O4 adjudicates — so it is defined as annotated/modeled slots, keeping "resolvable" out of the baseline.)*
- n_loci(*F*, individual) = the true number of distinct *F* paralog loci in the individual.

O4 partitions into **two typed predicates** (the original single "Δ > 0" conflated a cardinality gap with a content gap):

> **(COUNT-ABSENT)** Δ_cnt(F) = n_loci(F, individual) − n_loci(F, A) > 0 — **collapsed or wholly-absent copies: A has too few slots.**
>
> **(SEQUENCE-ABSENT)** ∃ a copy *c* of *F* whose sequence is not represented in A within het divergence (~0.5%), even though n_loci is matched — **divergent-in-place / mis-assembled paralog** (the MHC candidates are largely this case, not count-absent).

Neither n_loci(F, individual) nor *c* is directly observable. **The CONCEPT operationalizing COUNT-ABSENT** is the local copy read-conflict graph *H_l* at a reference locus *l*, built over **primary reads only** (paralog-bleed firewall): since **χ(H) is a lower bound** on true copy number (χ = MCC, our O2 object), the **sufficient condition χ(H_l) > n_loci(l) implies Δ_cnt > 0**. This direction is valid and self-consistent with the lower-bound framing.

**OPERATIONAL SURROGATE (explicitly a heuristic for χ, NOT a computed chromatic number):** `hidden_copy.rs::detect_hidden_copy` flags *l* when **≥12 alt columns** with alt-fraction in **[0.20, 0.60]** co-segregate across **≥5 reads**. These cutoffs are **calibrated, not derived**; they approximate "a coherent second haplotype the single reference slot cannot host." The approximation is **imperfect** — ordinary heterozygosity, gene conversion, and linked hets can also reach the balanced-column count (this is why het dominates the false-positive mass; see caveat b). **The detector reports column counts, not a chromatic number** — "χ(H_l) > n_loci" is the concept it approximates, and the two must be labeled concept-vs-implementation, not equated. **SEQUENCE-ABSENT** is instead evidenced by assembling the alt-haplotype reads and finding the consensus maps back at **3–20% divergence** (≫ het ~0.5%), too divergent to be that copy or its allele (the 4 MHC candidates).

**Ploidy note baked into the arithmetic (the root of the het wall):** n_loci(F, individual) is a **diploid** reality mapped against a **haploid** A. Counting at mismatched ploidy would score ordinary heterozygosity (one diploid locus, two alleles, one haploid slot) as reference-absent by construction — which is exactly why both counts are defined at **haploid locus resolution** and why the downstream het funnel (divergence > 3% ≫ het ~0.5%) is **necessary, not optional**. RNA cannot observe ploidy; this is the honest limit.

Discipline: **detect-and-flag only** — each O4 call is "the reads imply an unmodeled copy here," **never a placed sequence**.

### Baseline
The **reference assembly A and its completeness** — a **sample-vs-assembly, within-species (indeed within-individual) comparison**. A copy reference-absent vs an **incomplete** reference (GRCh38-style) can be reference-**present** vs a **complete** one (T2T): accordingly the **divergent-unmapped route (population B) is a real lever only against incomplete references and is DRY on complete T2T** (0.13% unmapped, 79% already 99.7%-present, 1 residual hit), so on T2T the yield is **entirely the mapped pile** (collapsed + divergent-in-place). This baseline is **distinct from and orthogonal to** an expansion's baseline (ancestral/outgroup-species count; cross-lineage; needs ≥2 species + outgroup).

### Map to our objects
**O4 is O2 turned outward:** O2 (max-weight facility location, assign-or-abstain, no 1/k) assigns reads to the copies **A has**; O4 flags copies **A lacks** via the local conflict graph demanding χ(H_l) > n_loci(l). Evidence type = **PSV/SUN over-dispersion**: excess balanced columns are **SUN-consistent** witnesses of an **unenumerated** copy (SUN ⊊ PSV, Sudmant 2010 — a *genuine* SUN only once the copy is enumerated; here it is SUN-**like** evidence inferring an unenumerated hidden copy). Discipline = detect-and-flag.

### Non-implications
(1) **Not an expansion** (different baseline/axis; unpolarizable without ≥2 species + outgroup). (2) **Not a placed/assembled copy** (detect-and-flag). (3) **Not a resolved copy-vs-allele call** — the diploid-individual vs haploid-A ploidy mismatch makes a hyperdivergent **allele** indistinguishable from a distinct **copy** from RNA alone (the het-vs-copy wall). (4) **Not a novel gene/contaminant** (endogeneity via proteome BLAST; population B contaminant-screened). (5) **Not a functional claim** beyond ORF/expression. (6) **Not assembly-independent** ("absent" is relative to a specific A and its completeness).

### Our-usage caveats
(a) **Strict lower bound:** near-identical collapsed copies with no PSV (K=0 / NM:i:0 floor, the same identifiability floor as χ(H)) and low-expression/low-depth copies are invisible; the ≥12-column requirement also drops them. (b) **The raw flag is a non-specific SCREEN dominated by heterozygosity** (measured FP bound **828/11,206 = 7.39%** on single-copy genes vs ~7.93% genome background) — a direct consequence of the diploid-vs-haploid ploidy mismatch, not merely noise. Specificity comes **only** from the downstream funnel (divergence > 3% ≫ het ~0.5%; A→G + T→C editing filter; protein BLAST), collapsing **1,015 flags → 73 candidates → 15 dispersed + 4 protein-confirmed MHC**. (c) **Circularity (adversarial review #4):** the "divergence ⇒ copy" promotion is circular because the alt-haplotype reads are **selected** for reference-divergent columns, so the assembled consensus is divergent **by construction** — flagged, not resolved. (d) The **4 MHC + 15 dispersed are candidates**; copy-vs-allele needs **DNA parCN (QuicK-mer2)**, installed and verified-to-run but **never completed** (~40 GB index RAM, no gorilla DNA WGS — only RNA IsoSeq on disk). The read-free structural surrogate (haploid T2T n_loci + SEDEF partner) resolves only ~44% (26 copies/18% + 38 alleles/26%) and confirms ~45% genuinely DNA-ambiguous — so "copy-vs-allele needs DNA" stands for the MHC bulk. (e) Population-B divergent-unmapped route **DRY on complete T2T**. (f) **No external validation yet** (no DNA, no second individual, no assembled full-copy locus; BLAST confirms endogeneity, not copy status).

### Literature
Bailey et al. 2002, *Science* 297:1003 (SD; WSSD depth CN — ancestor of collapse detection; coined PSV = the copy-vs-allele problem); Sudmant et al. 2010, *Science* 330:641 (SUN); Soto/Dennis 2025, *Cell* 188:5363–5383, DOI 10.1016/j.cell.2025.06.037 (213 human-specific families / 1,002 paralogs, **37% (668/1,002) in GRCh38-missing/erroneous regions**; famCN/parCN; SD98 short-read SNV sensitivity **0.85%** = the long-read + collapse-aware rationale and the expansion↔collapse interaction source); Vollger et al. 2019 (SDA — PSV correlation-clustering of collapsed segdups from long reads; closest prior art; shares the K=0 floor); QuicK-mer2/fastCN (Eichler lineage — the DNA parCN resolver for copy-vs-allele); Fitch 1970, *Syst. Zool.* 19:99 (the paralog is the counted unit).

---

## 6. THE BASELINE LATTICE (the organizing idea)

Every concept above is a comparison, and each is measured against a **different reference frame**. Confusing two concepts is almost always the act of silently swapping one baseline for another. Stated exactly:

| Concept | Kind of object | **Baseline (what it is measured AGAINST)** | Species / genomes needed | Polarizable? |
|---|---|---|---|---|
| **Paralog** (vs ortholog) | pairwise evolutionary relation | **The gene phylogeny** — true gene genealogy (estimand), approximated by the reconciled gene tree; label = type of the lca node (δ vs σ) | 1 genome suffices to *know* two same-genome homologs are paralogs (within-genome collapse); ≥2 species to *label* ortholog vs the σ/δ discrimination | n/a (a relation, not a change) |
| **Segmental duplication** | DNA/structural block state | **The genome self-alignment** — duplicated ≥1 kb, ≥90% id **elsewhere in the same assembly** (WGAC ∩ WSSD) | 1 genome (1 individual/haplotype) | No |
| **Multi-copy gene family (O1)** | within-genome membership + count state | **One genome** — ≥2 transcribed homologous copies at ≥2 loci; count = max(Λ(C), χ(H)) | 1 genome | No |
| **Expansion** | cross-lineage directional change | **A reconstructed ancestral / outgroup copy count** — one species-tree branch back in time | **≥2 species + an outgroup** (to reconstruct + polarize) | **Yes — requires it** |
| **Reference-absent copy (O4)** | sample-vs-assembly gap | **A specific reference assembly** (and its completeness) — true individual count/content minus what A hosts | 1 species / 1 individual | No |

**Two structural facts about the lattice.** (i) The **four homology oracles E_a/E_b/E_r/E_p are pairwise incomparable** (each has witnessed set-differences both ways); the only clean containment in the whole stack is **E_c ⊆ E_b** (read-ambiguity ⊆ exonic homology). There is **no total order** of "closeness to paralogy" — only a per-family recall heuristic (E_p uniquely reaches twilight-zone paralogy). (ii) The five baselines are genuinely different reference frames: a paralog needs no assembly and no second species; an SD needs no second species; a multi-copy family needs no ancestor and no assembly comparison; an expansion needs an ancestor/outgroup; an O4 copy needs an assembly. **Naming the baseline first disambiguates every one of these concepts.**

---

## 7. EXPANSION vs REFERENCE-ABSENT (O4) — the crisp distinction

This is the student's specific question. The two are the two concepts on the **comparative axis**, and they are the pair most often conflated — because both are "there is more here than the naïve baseline shows." They are **different baselines, orthogonal, but operationally coupled.** Get all three parts right.

### 7.1 Different baselines
- **Expansion** is measured against a **reconstructed ancestral / outgroup copy count** — a point on the species phylogeny, **one branch back in time**. It is **cross-lineage** and requires **≥2 species + an outgroup** to polarize *which lineage changed*. Δ over species. Answers: *did this lineage gain copies relative to its ancestor?*
- **Reference-absent (O4)** is measured against a **specific reference assembly** — the same individual's own genome vs the bytes in A. It is **within one species (indeed one individual)** and requires **no second species**. Δ between sample and assembly. Answers: *does the reference host every copy the reads imply?*

One is a fact about **evolution across lineages**; the other is a fact about **assembly completeness for one sample**. Different reference frames, different data requirements, different failure modes.

### 7.2 Orthogonal — each can occur without the other
- **Reference-absent yet NOT expanded:** an assembly collapses two **near-identical copies** of a family whose copy count is **identical across all apes** (no lineage gain). The gap is an assembly artifact, not an evolutionary event. → O4 fires, expansion does not.
- **Expanded yet NOT reference-absent:** a lineage-specific expansion whose every copy is **correctly resolved** by a complete T2T assembly. The evolutionary event is real and the assembly is complete. → expansion holds, O4 fires nowhere.

These are the two independent axes; neither implies the other.

### 7.3 But they interact — and O4 is the TOOL that makes expansion robust to collapse
This is Soto's rationale and the reason the two concepts are studied together. **Expanded families are precisely the most collapse-prone:** recent duplication → high sequence identity → assembly and short-read collapse. So the assembly is **most likely to under-count exactly the families where the expansion is**. Concretely (Soto 2025): **37% (668/1,002) of human-specific paralogs fall in GRCh38-missing/erroneous regions**, and **short-read SNV sensitivity is 0.85% in SD98**. Without collapse-aware, long-read true-count recovery, the expansion counts are **biased downward exactly where the expansion is** — assembly collapse can **fake** an ancestor (over-collapse the baseline) or **hide** the derived count (under-count the tip), either of which corrupts the Δ.

**O4 supplies the true-count recovery** (on the RNA/individual side; DNA-side it is famCN/parCN via WSSD/QuicK-mer2). By flagging copies the assembly does not host — χ(H_l) > n_loci(l), the reads over-determining the reference slot — O4 **de-biases χ(H) upward toward the true count**, making the expansion contrast **robust to assembly collapse rather than an artifact of it**. In our object terms: expansion = Δχ(H) across species; O4 = the mechanism that keeps each per-species χ(H) from being silently depressed by collapse. Orthogonal in definition, coupled in practice.

### 7.4 The one-sentence version
**Expansion asks "did this lineage gain copies vs its ancestor?" (baseline = ancestor/outgroup, needs ≥2 species). Reference-absent asks "does the reference host every copy this individual has?" (baseline = the assembly, needs one individual). They can each occur alone, but because expanded families are the most collapse-prone, O4 (true-count recovery) is the tool that keeps the expansion measurement honest.**

---

## 8. WHERE OUR PIPELINE SITS — honest mapping

| Concept | Our object | What we actually measure | Honest status |
|---|---|---|---|
| **Paralog** | E_r (proxy for within-genome paralog SET) | paralog-group **membership**, not the Fitch relation | We recover the *set* approximately; we do **not** recover the ortholog/paralog *relation*, the gene tree, ancestry, or duplication history. |
| **Segmental duplication** | E_a (SEDEF, corroborating oracle only) | a within-genome DNA block state; **WGAC-only** (no WSSD/famCN — no gorilla WGS depth) | Structural corroboration, **one-sided** (SD98 membership confirms, non-membership does not refute). Not our family definition. |
| **Multi-copy gene family** | **O1 (primary object)** = γ-quasi-clique-refined E_r component, ≥2 loci | per-genome copy number as **max(Λ(C), χ(H) = MCC)**, two independent censored lower bounds | **χ(H) = transcribed, resolvable copies = a functional LOWER BOUND** on true genomic copy number. This is our headline object. |
| **Expansion** | Δχ(H) across species, for an **annotation-supplied** orthology map φ | a two-species **difference** in χ(H) (MAGEA 2→11, TSPY 5→33 vs RBMY 6=6) | We defend "χ(H) **differs** across species, a conservative lower bound **consistent with** the annotated human-specific expansion" — **not** an independently inferred and polarized lineage-specific expansion (no outgroup, no CAFE, φ borrowed). |
| **Reference-absent copy** | **O4** = χ(H_l) > n_loci(l) detector (hidden_copy) | a **sample-vs-assembly** flag: reads over-determine the reference slot | **Detect-and-flag, never placed.** A non-specific screen (het-dominated 7.39% FP) funneled to 15 dispersed + 4 MHC **candidates**; copy-vs-allele **needs DNA** (parCN never completed). |

**The through-line.** All five sit on one **identifiability** spine: **χ(H) = MCC** is the load-bearing quantity. Inward it counts resolvable copies (O1's copy number); it is a **provable lower bound** (Lemma 1: MCC = χ(H) on the de-tied graph; both Λ(C) and χ(H) ≤ true count). Across species its Δ is our expansion signal (a conservative lower bound on famCN). Turned outward, χ(H_l) > n_loci(l) is the O4 reference-absent flag. What we **honestly claim**: a single-genome, transcribed, resolvable **copy catalog with a per-family lower-bound copy number**, plus **corroboration** (SD/annotation) and **screens** (expansion difference, O4 flag) that are explicitly disclosed as borrowing polarization/correspondence/DNA where we lack it.

---

## 9. What to tell the advisor

*"These are four measurement axes, not synonyms, and the whole design turns on naming the baseline. A **paralog** is a pairwise relation read off the gene tree — and because we work in one gorilla genome, every within-family locus pair is a paralog by construction, so we define families (**E_r**, our O1) without ever inferring a duplication node; we recover the paralog **set**, not the Fitch relation, and we say so. A **segmental duplication** is a within-genome DNA block (**E_a**), used only as one-sided corroboration — our oracles E_a/E_b/E_r/E_p are pairwise incomparable, so I never claim read-conflict families are a subset of segdups, and PCDHB (88.9%, below the Bailey floor) is my proof that RNA+PSVs reach past the DNA cutoff. Our primary object is the **multi-copy gene family**, whose per-genome copy number I report as two independent censored lower bounds, Λ(C) and **χ(H) = MCC** — provably a lower bound, never overstated as the true count. The two concepts I keep strictly apart are **expansion** and **reference-absent (O4)**: expansion is cross-lineage vs a reconstructed ancestor and needs ≥2 species plus an outgroup to polarize — I have gorilla+human only, so my MAGEA 2→11 and TSPY 5→33 are **differences** whose polarization **and orthogroup correspondence are borrowed from annotation**, and I state exactly that; O4 is a within-individual sample-vs-assembly gap, a detect-and-flag screen dominated by heterozygosity until the divergence/editing/BLAST funnel, with copy-vs-allele still owed to DNA I don't yet have. They are orthogonal but coupled: expanded families are the most collapse-prone, so O4 is the tool that keeps χ(H) — and therefore the expansion Δ — honest against assembly collapse. Every threshold (γ_core = 0.13, γ = 0.20, the ≥12-column O4 heuristic) is named and owned, never dressed as significance; no circularity in the family definition; the borrowings are disclosed, not hidden."*

---

*Citations verified against project memory and the entries' skeptic-approved forms. Soto/Dennis 2025 confirmed as Cell 188:5363–5383, DOI 10.1016/j.cell.2025.06.037. Gorilla substrate = mGorGor1 (RNA donor == assembly individual, confirmed); the T2T accession (GCF_029281585.2) should be pinned against the exact assembly mapped against before this goes in a written defense. No invented citations.*