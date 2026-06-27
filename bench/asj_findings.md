# Modeling & quantifying allele-specific junctions (genome-wide)

**The idea (advisor interest):** a splice junction is *allele-specific* when its usage depends on the
allele a molecule carries. This is a quantifiable allele/copy phenomenon — not assembly. **Long HiFi
reads make it direct:** a single read observes BOTH its heterozygous-SNP allele AND which junctions it
splices, so allele→junction linkage is per-molecule — no statistical phasing (the short-read sQTL
problem). **Substrate = the het loci that *confounded* copy-detection in task (c):** heterozygosity is
exactly the phasing signal here. The confound becomes the feature.

## Method (bench/allele_specific_junctions.py → asj_aggregate.py → asj_verify.py)
Per gene: call a balanced het **anchor SNP** from the pileup (biallelic, freq∈[0.35,0.65], ≥12×, ≥5
reads/allele) → partition reads by the allele they carry. Per junction (seen ≥3×), among reads carrying
each allele AND spanning the junction, compute **PSI = used/spanning**; **Fisher exact** on the 2×2
(allele × uses-junction); effect = |ΔPSI|. **BH-FDR** genome-wide; **ASJ = q<0.05 AND |ΔPSI|≥0.3**.

### Definition (PSI / ΔPSI, as implemented in `vg_family/allele_specific_junctions.rs`)

Fix a balanced heterozygous **anchor SNP** with alleles $x,y$, and partition the reads covering it by
the allele each molecule carries (per-molecule, from the same HiFi read — no statistical phasing). For
a candidate junction $j=(d,a)$ (donor $d$, acceptor $a$) and allele $x$, define over the reads carrying
allele $x$:

$$S_x(j)=\#\{\,r : \mathrm{ref\_start}(r)\le d \ \wedge\ \mathrm{ref\_end}(r)\ge a\,\}\qquad\text{(spanning — denominator)}$$
$$U_x(j)=\#\{\,r \text{ spanning} : (d,a)\in \mathrm{intron\_chain}(r)\,\}\qquad\text{(using — numerator)}$$
$$\boxed{\ \mathrm{PSI}_x(j)=\frac{U_x(j)}{S_x(j)}\in[0,1]\ }\qquad
\Delta\mathrm{PSI}(j)=\big|\,\mathrm{PSI}_x(j)-\mathrm{PSI}_y(j)\,\big|.$$

$\mathrm{PSI}_x(j)$ is the fraction of allele-$x$ molecules that *span* the junction locus which actually
*splice out* $j$. A junction is called **allele-specific** when the Fisher exact test on

$$\begin{pmatrix} U_x & S_x-U_x\\[2pt] U_y & S_y-U_y\end{pmatrix}$$

passes genome-wide Benjamini–Hochberg FDR ($q<q^\*$) **and** $\Delta\mathrm{PSI}\ge\Delta\mathrm{PSI}_{\min}$
(defaults $q^\*=0.05$, $\Delta\mathrm{PSI}_{\min}=0.30$). Guards: $\ge$ `min_span` ($5$) spanning reads
**per allele**, junction seen $\ge$ `min_j` ($3$) times; junctions constitutive in *both* alleles
($\mathrm{PSI}\ge0.98$ on both, or $\le0.02$ on both) are skipped as un-testable.

This is a **junction-level, per-allele** PSI (inclusion ratio of one splice junction among the
molecules of that allele which span it) — same used/(used+skipped) family as the classic
percent-spliced-in, but the unit is a junction and the denominator is allele-specific spanning reads,
which is what makes it a true per-molecule allele→junction linkage rather than a phased/sQTL estimate.

## Result
- **7,898** phaseable het genes; **74,674** alternatively-spliced junctions tested.
- **Headline (the genetic core): 120 allele-specific junctions with a TRANSVERSION anchor** (FDR q<0.05,
  |ΔPSI|≥0.3, 59 genes) — unambiguously **genetic** (a transversion anchor cannot be an A→I RNA-edit
  site), so these are allele-specific by construction, not editing-coupled.
- **475 total candidates** across 235 genes (median |ΔPSI| **0.64**; **213** ≥0.7; **56** full switches
  ΔPSI=1.0) — but the other **355 have transition (A/G, C/T) anchors** that *could* be RNA-edit-coupled
  splicing (real, but not *genetic* allele-specificity). Lead with the 120; the 475 is the full
  candidate set including the edit-confoundable transitions.

## Verification (deterministic + mechanistic — stronger than an LLM pass here)
- **Collapsed-paralog masquerade — partially controlled, honestly bounded.** The `frac_mq0<0.1` filter
  removed **0/475** (the called anchors all sit in uniquely-mappable flank), so it is **not a binding
  control** — it cannot, by construction, separate a 50/50 *het* from a 50/50 *two-copy PSV* whose flank
  is uniquely mappable. ~36% of the called genes are `LOC*` paralog loci where this ambiguity is live.
  The clean separator is the existing `scan_gene_copy_specific_junctions` (copy-specific vs het),
  which should be run on the LOC* loci and reported as a **separate within-gene-het vs paralog-locus
  count** (TODO). The **120 transversion / ~20 splice-proximal** core is unaffected (mechanism is
  base-level, see below).
- **Editing controlled (this is the load-bearing control):** transversion anchors (120) cannot be
  RNA-edit sites; the 355 transition anchors are flagged as edit-confoundable.
- **Mechanism — textbook splice-site variants.** 20 high-confidence ASJ are splice-proximal (anchor
  ≤100bp from the junction); the cleanest sit **on the canonical splice dinucleotide**, with the
  per-molecule split confirming disruption:
  - **PSMD2**: SNP at the splice **donor** (GT-AG). allele G → 14/14 splice (PSI 1.0); allele T → 0/18 (PSI 0).
  - **DAXX**: SNP at the splice **acceptor** (CT-AC, − strand). allele C → 9/9; allele A → 0/10.
  - Balanced het depth on both alleles rules out reference error; the SNP creates/destroys the splice motif.
- **Long-range allele-specific transcript structure** (reported separately): LOC101141065 and MYO18B —
  one anchor (itself a splice-site SNP) flips a whole *set* of junctions up to 15–100 kb away (ΔPSI=1.0),
  i.e. the allele selects an entire transcript structure. Real (uniquely mapped, anchor at a splice site)
  but a distinct phenomenon from local splice-disruption; long-span ones are also lower-power.
  > ⚠ **Chimera-not-excluded class.** The ASJ path has **no RT/template-switch guard**. The ~23 distal,
  > full-switch (|ΔPSI|=1.0, anchor far from the junction) calls — 17 of them `LOC*` loci, 7 reaching the
  > high-confidence column — are exactly the class an RT/template-switch chimera (or trans-association)
  > could mimic, since a single molecule carrying both the distal allele and the switched junction is the
  > artifact's own signature. These should be tagged **chimera-not-excluded** in the TSV (and could carry
  > the microhomology-signature flag that already exists in `copy_assign`); the local splice-proximal core
  > (120 transversion / ~20 dinucleotide) is unaffected because its mechanism is base-level at the junction.

## Honest caveats
- **Single-anchor phasing:** one balanced het SNP per gene; a junction is tested only among reads
  spanning both the anchor and the junction (long reads help). Multi-SNP haplotype phasing adds reach/power.
- **Transition anchors** (A/G, C/T; the non-transversion 355) could include RNA-edit sites — edit-coupled
  splicing is real but not *genetic* allele-specific; treat the 120 transversion calls as the genetic core.
- **Coverage-limited** (≥12× anchor, ≥5 spanning/allele); deeper data finds more. Long-range full-switches
  need rare long-spanning reads → lower power, wider CIs.
- **Scope:** RefSeq gene loci on the 25 large contigs.

## Why this fits the project
Quantifiable, principled (FDR + |ΔPSI| + per-molecule mechanism), **allele/copy-resolution not assembly**,
and a genuine **long-read advantage** (per-molecule allele→junction linkage). It reuses the task-(c)
machinery wholesale (de-novo allele calling + the het loci) — turning the copy-detection confound into
a first-class measurement. Natural extensions: multi-SNP haplotype phasing; strand-aware editing
demoter to reclassify the transition anchors; tie ASJ to the PSV/copy axis (copy-specific junctions at
real paralog loci); cross-reference splice-proximal hits against known splice-variant catalogs.

## Reproduce
- `MINIFORGE python bench/allele_specific_junctions.py --chrom <C>` (per-chrom; `--region` for one locus)
- `python3 bench/asj_aggregate.py` (BH-FDR + calls → asj_results.md, asj_calls.tsv)
- `python3 bench/asj_verify.py` (confound control → asj_calls_verified.tsv)
- `python3 bench/asj_evidence.py --row N` (per-call per-molecule evidence + splice motifs)
- `python3 bench/asj_fig.py` (figure)
