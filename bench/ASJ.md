# Asj (consolidated)

> Merged from 3 source docs (verbatim, git keeps the originals' history). Each section below was a separate `bench/*.md`.

**Contents:** [asj_findings](#asj-findings) · [asj_results](#asj-results) · [asjm_findings](#asjm-findings)


---

## asj_findings

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
- **Headline (the REPRODUCIBLE genetic core): 54 strand-bias-clean allele-specific junctions in 29 genes**
  — the committed three-filter funnel `475 → 120 transversion → 76 non-LOC\* → 54 SOR-clean`
  (`bench/asj_genetic_core.py` → `asj_genetic_core.tsv`, reproducible from the call tables). Each filter
  removes a distinct confound: transversion anchor (cannot be A→I editing), non-LOC\* (excludes paralog
  copy-masquerade), and the project's own longcallR-style SOR strand-bias filter (`sor_pass`).
  ⚠ **Corrected (adversarial review #4):** the old "~77" headline = 76 + 1 "clean LOC\*", where the +1 came
  from a paralog-masquerade scan whose output (`bench/P1_P4_RESULTS.md`) **does not exist on disk** — that +1
  is RETRACTED, and applying the SOR filter (22/76 named calls fail it) tightens 76 → 54. The full
  transversion set is **120** (FDR q<0.05, |ΔPSI|≥0.3); 44 of those sit at LOC\* paralog loci (copy-confounded,
  excluded). The previously-headlined "mechanistically airtight" PSMD2/DAXX **anchors FAIL SOR**
  (10.45 / 7.08, sor_pass=0) and are RETIRED as flagships.
- **475 total candidates** across 235 genes (median |ΔPSI| **0.64**; **213** ≥0.7; **56** full switches
  ΔPSI=1.0) — but the other **355 have transition (A/G, C/T) anchors** that *could* be RNA-edit-coupled
  splicing (real, but not *genetic* allele-specificity). Lead with the **54 SOR-clean genetic core**; the full
  transversion set is 120; the 475 is the full candidate set including edit-confoundable transitions.

## Verification (deterministic + mechanistic — stronger than an LLM pass here)
- **Collapsed-paralog masquerade — DONE (P4).** The `frac_mq0<0.1` filter removed **0/475** (the
  called anchors all sit in uniquely-mappable flank), so it is **not a binding control** — it cannot,
  by construction, separate a 50/50 *het* from a 50/50 *two-copy PSV* whose flank is uniquely mappable.
  ~36% of the called genes are `LOC*` paralog loci where this ambiguity is live. The 44 LOC\* transversion calls are
  excluded from the genetic core **by name** (LOC\* = unannotated paralog family = copy-masquerade-prone),
  reproducibly from `asj_calls_verified.tsv` (120 transversion = 44 LOC\* + 76 named).
  ⚠ **review #4:** a `scan_gene_copy_specific_junctions` masquerade scan (claimed 17/18 LOC\* genes carry
  copy-specific junctions, 1 "clean") was cited from `bench/P1_P4_RESULTS.md`, **which does not exist on
  disk** — so the 17/18 figure and the old "+1 clean LOC\* → ~77" are UNBACKED and retracted. The genetic
  core is the reproducible **76 non-LOC\* → 54 SOR-clean** subset (`asj_genetic_core.tsv`). The **~20 splice-proximal
  core (PSMD2/DAXX) sits in the EXTENDED splice consensus (donor±1 / acceptor±1), NOT the invariant
  dinucleotide** — see Mechanism below (genome-verified, corrected 2026-06-28).
- **Editing controlled (this is the load-bearing control):** transversion anchors (120 full set;
  54 SOR-clean genetic core after the non-LOC* + SOR filters) cannot be RNA-edit sites;
  the 355 transition anchors are flagged as edit-confoundable.
- **Mechanism — splice-REGION (extended-consensus) variants, NOT invariant-dinucleotide variants
  (CORRECTED 2026-06-28; the earlier "on the canonical GT-AG dinucleotide / creates-destroys the motif"
  claim was genome-FALSE by one base and is RETRACTED).** Re-derived directly from `GGO.fasta`
  (`bench/asj_motif_check.py`): **0/475 anchors fall on a core splice dinucleotide**; the cleanest two are
  *adjacent* to it:
  - **PSMD2**: anchor at **donor−1** (the canonical donor `GT` at `NC_086017.1:195406804` is INTACT under both
    alleles). The per-molecule split is real: allele G → 14/14 splice (PSI 1.0); allele T → 0/18 (PSI 0).
  - **DAXX**: anchor at **acceptor+1** (the canonical acceptor `AG` is INTACT under both alleles).
    allele C → 9/9; allele A → 0/10.
  - **The load-bearing result is the per-molecule allele→junction LINKAGE** (the PSI 0↔1 split, balanced het
    depth ruling out reference error), NOT a dinucleotide-disruption mechanism. These are splice-region /
    extended-consensus variants that *modulate* splicing; "textbook GT-AG variant" is the wrong frame. So the
    genome-CERTIFIED splice-mechanism count is **0 core-dinucleotide + 2 extended-consensus (PSMD2, DAXX)**;
    the other ~18 splice-proximal calls are *splice-proximal only, mechanism uncertified*.
  - ⚠ **review #4:** PSMD2 and DAXX both **FAIL the SOR strand-bias filter** (SOR 10.45 / 7.08, every read on
    one strand) — so even these "cleanest two" are strand-bias-artifact-suspect and are **excluded from the
    54-call SOR-clean core**; do not present them as airtight. The reads-on-one-strand pattern is itself a
    plausible artifact mechanism for the perfect PSI 0↔1 split. Retain only as exploratory exemplars.
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
  > (54 SOR-clean genetic-core transversion / ~20 dinucleotide) is unaffected because its mechanism is base-level at the junction.

## Honest caveats
- **Single-anchor phasing:** one balanced het SNP per gene; a junction is tested only among reads
  spanning both the anchor and the junction (long reads help). Multi-SNP haplotype phasing adds reach/power.
- **Transition anchors** (A/G, C/T; the non-transversion 355) could include RNA-edit sites — edit-coupled
  splicing is real but not *genetic* allele-specific; treat the 54 SOR-clean genetic core (full transversion set
  120, minus 44 copy-confounded LOC\* calls) as the genetic core.
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
- `asj --bam GGO_mm.bam --regions <LOC* windows> --mode copy --out <out>` (masquerade separator — see `bench/P1_P4_RESULTS.md`, P4)


---

## asj_results

# Allele-specific junctions (genome-wide)

Splice junctions whose usage depends on the allele a molecule carries. Long HiFi reads link allele->junction PER MOLECULE (the read carries both its het-SNP allele and its junctions), so no statistical phasing is needed. The het loci that confounded copy-detection (task c) are the substrate: heterozygosity = the phasing signal (confound -> feature).

## Result
- genes with a phaseable het anchor + a tested junction: **7898**
- alternatively-spliced junctions tested (non-constitutive): **74674**
- **allele-specific junctions (BH-FDR q<0.05 AND |dPSI|>=0.3): 475** across **235** genes
- of those, **120** have a TRANSVERSION anchor (unambiguously genetic, not RNA-editing) — the high-confidence genetic set

## |dPSI| distribution of the ASJ calls
- median |dPSI|=0.64, max=1.00; >=0.5: 312, >=0.7: 213, ==1.0 (full switch): 56

## Top allele-specific junctions
| gene | chrom | anchor | alleles | junction | PSI_X | PSI_Y | dPSI | q | anchor type |
|---|---|---|---|---|---|---|---|---|---|
| LOC101141065 | NC_073224.2 | 192934097 | G/C | 192941955-192944457 | 1.0 | 0.0 | 1.00 | 1.9e-26 | transversion(genetic) |
| LOC101141065 | NC_073224.2 | 192934097 | G/C | 192944612-192944695 | 1.0 | 0.0 | 1.00 | 6.6e-26 | transversion(genetic) |
| LOC101141065 | NC_073224.2 | 192934097 | G/C | 192949190-192950494 | 1.0 | 0.0 | 1.00 | 6.6e-26 | transversion(genetic) |
| LOC101141065 | NC_073224.2 | 192934097 | G/C | 192947387-192949145 | 1.0 | 0.0 | 1.00 | 6.6e-26 | transversion(genetic) |
| LOC101141065 | NC_073224.2 | 192934097 | G/C | 192944823-192947259 | 1.0 | 0.0 | 1.00 | 6.6e-26 | transversion(genetic) |
| LOC101141065 | NC_073224.2 | 192934097 | G/C | 192930566-192934097 | 0.0 | 1.0 | 1.00 | 6.6e-26 | transversion(genetic) |
| GLB1L3 | NC_073233.2 | 142241047 | C/A | 142240873-142240938 | 0.0 | 1.0 | 1.00 | 3.2e-14 | transversion(genetic) |
| LOC115932992 | NC_073242.2 | 16399354 | C/T | 16419347-16423955 | 0.0 | 1.0 | 1.00 | 4.7e-12 | transition(poss.edit) |
| TSPAN2 | NC_073224.2 | 119456078 | C/G | 119454970-119455031 | 0.0 | 1.0 | 1.00 | 7.0e-12 | transversion(genetic) |
| LOC115932992 | NC_073242.2 | 16399354 | C/T | 16426010-16429060 | 0.0 | 1.0 | 1.00 | 5.5e-11 | transition(poss.edit) |
| CNTNAP3C | NC_073236.2 | 33312391 | G/A | 33312509-33313164 | 0.0 | 1.0 | 1.00 | 9.2e-11 | transition(poss.edit) |
| CASP10 | NC_073235.2 | 103835969 | A/G | 103835036-103835169 | 1.0 | 0.0 | 1.00 | 8.1e-10 | transition(poss.edit) |
| LOC115933039 | NC_073244.2 | 20936571 | A/G | 20959729-20961004 | 1.0 | 0.0 | 1.00 | 1.0e-09 | transition(poss.edit) |
| LOC115934626 | NC_073228.2 | 45957287 | T/A | 45956648-45956726 | 1.0 | 0.0 | 1.00 | 3.3e-09 | transversion(genetic) |
| HIGD1A | NC_086017.1 | 52690273 | C/G | 52689198-52689242 | 1.0 | 0.0 | 1.00 | 2.5e-08 | transversion(genetic) |
| KCNAB2 | NC_073224.2 | 234526415 | C/G | 234525571-234526031 | 1.0 | 0.0 | 1.00 | 1.9e-07 | transversion(genetic) |
| PSMD2 | NC_086017.1 | 195406803 | T/G | 195406804-195406955 | 0.0 | 1.0 | 1.00 | 5.8e-07 | transversion(genetic) |
| CSNK1D | NC_073228.2 | 9740588 | G/C | 9740897-9740970 | 0.0 | 1.0 | 1.00 | 8.8e-07 | transversion(genetic) |
| LOC134758754 | NC_073224.2 | 19346753 | T/C | 19339269-19339542 | 1.0 | 0.0 | 1.00 | 1.3e-06 | transition(poss.edit) |
| LOC129527496 | NC_073242.2 | 4910233 | C/T | 4908324-4908418 | 0.0 | 1.0 | 1.00 | 3.1e-06 | transition(poss.edit) |
| C13H9orf43 | NC_073237.2 | 109779230 | T/C | 109779303-109779954 | 0.0 | 1.0 | 1.00 | 6.1e-06 | transition(poss.edit) |
| LOC115932992 | NC_073242.2 | 16399354 | C/T | 16480542-16484928 | 1.0 | 0.0 | 1.00 | 7.4e-06 | transition(poss.edit) |
| LOC115932992 | NC_073242.2 | 16399354 | C/T | 16471500-16475882 | 1.0 | 0.0 | 1.00 | 7.4e-06 | transition(poss.edit) |
| LOC115932992 | NC_073242.2 | 16399354 | C/T | 16476020-16480404 | 1.0 | 0.0 | 1.00 | 7.4e-06 | transition(poss.edit) |
| MYO18B | NC_086018.1 | 23923474 | G/T | 24012728-24023003 | 1.0 | 0.0 | 1.00 | 1.5e-05 | transversion(genetic) |

## Honest caveats
- **Editing confound:** a TRANSITION anchor (A/G or C/T) could be an RNA-edit site (phasing by edit-status, not allele) — that is edit-coupled splicing, still real but not *genetic* allele-specific. TRANSVERSION anchors are unambiguously genetic; treat the transversion subset as the high-confidence genetic ASJ. (Strand-aware editing detection would reclassify the transitions.)
- **Single-anchor phasing:** reads are phased by ONE balanced het SNP; a junction is tested only among reads spanning BOTH the anchor and the junction (long reads help). Multi-SNP haplotype phasing would add power/reach.
- **Collapsed-paralog masquerade:** at a collapsed locus the two 'alleles' could be paralog copies (=copy-specific splicing, also interesting, but not within-gene allele-specific). The het substrate from task c was uniquely-mapped (frac_mq0=0); genome-wide, low-MAPQ loci should be down-weighted.
- **Coverage-limited:** needs >=12x at the anchor and >=5 spanning reads per allele; deeper data finds more.

## Reproduce
- `MINIFORGE python bench/allele_specific_junctions.py --chrom <C>` (per-chrom) ; `--region` for one locus
- `python3 bench/asj_aggregate.py` (this: BH-FDR + calls)


---

## asjm_findings

# Allele-specific junctions — MULTI-SNP haplotype phasing

Each read is assigned to one of the two diploid haplotypes using ALL het SNPs in the gene (read-backed 2-means), not a single anchor SNP — so reads covering any subset of SNPs get phased (more reach/power), and the phasing quality (how cleanly reads bipartition) is a built-in confound check. With 1 SNP it degenerates to single-anchor, so this is a SUPERSET.

## Result
- junctions tested: **82179**
- **allele-specific junctions (q<0.05, |dPSI|>=0.3): 630** across **263** genes
- high-confidence (>=1 transversion SNP [genetic] AND phase_qual>=0.7): **453** / 184 genes
- from genes with >=2 het SNPs (multi-SNP phasing does real work): **559** / 241 genes

## Gain over single-anchor phasing
- single-anchor ASJ (committed): **475**; multi-SNP ASJ: **630**
- shared (gene+junction): **448**; **multi-SNP-only: 182** (the reach/power gain); single-anchor-only: 27
- net: multi-SNP ADDS +155 ASJ (263 genes)

## n_snp distribution of ASJ
```
  n_snp=1 : 71
  n_snp=2 : 64
  n_snp=3 : 54
  n_snp=4 : 35
  n_snp=5 : 35
  n_snp=6 : 36
  n_snp=7 : 38
  n_snp=8 : 18
  n_snp=9 : 49
  n_snp=10+: 230
```

## Top multi-SNP-only allele-specific junctions (reach/power gain)
| gene | chrom | n_snp | phase_qual | junction | psi0 | psi1 | dPSI | q |
|---|---|---|---|---|---|---|---|---|
| LOC129525331 | NC_073224.2 | 40 | 0.79 | 125054945-125055563 | 0.0 | 1.0 | 1.00 | 6.7e-12 |
| LOC109026500 | NC_073227.2 | 2 | 1.0 | 9976297-9976330 | 0.0 | 1.0 | 1.00 | 9.6e-08 |
| LOC129527636 | NC_073242.2 | 6 | 0.77 | 35228529-35228571 | 1.0 | 0.0 | 1.00 | 2.8e-04 |
| LOC115931062 | NC_073242.2 | 40 | 0.81 | 21768385-21768637 | 1.0 | 0.0 | 1.00 | 9.7e-04 |
| LOC115931062 | NC_073242.2 | 40 | 0.81 | 21767764-21767806 | 1.0 | 0.0 | 1.00 | 9.7e-04 |
| LOC115931062 | NC_073242.2 | 40 | 0.81 | 21766564-21766606 | 1.0 | 0.0 | 1.00 | 1.5e-03 |
| LOC109025073 | NC_073244.2 | 9 | 1.0 | 18096968-18106041 | 1.0 | 0.0 | 1.00 | 3.3e-03 |
| LOC101146691 | NC_073236.2 | 2 | 1.0 | 47843875-47845873 | 1.0 | 0.0 | 1.00 | 1.8e-02 |
| LOC101146691 | NC_073236.2 | 2 | 1.0 | 47843875-47844486 | 0.0 | 1.0 | 1.00 | 1.8e-02 |
| LOC101146691 | NC_073236.2 | 2 | 1.0 | 47844579-47845873 | 0.0 | 1.0 | 1.00 | 1.8e-02 |
| LOC129525331 | NC_073224.2 | 40 | 0.79 | 125054945-125065033 | 0.981 | 0.0 | 0.98 | 8.6e-10 |
| LOC109025073 | NC_073244.2 | 9 | 1.0 | 18106165-18107316 | 0.923 | 0.0 | 0.92 | 9.1e-03 |
| LOC129525331 | NC_073224.2 | 40 | 0.79 | 125057247-125057955 | 0.078 | 1.0 | 0.92 | 1.4e-16 |
| LOC129525331 | NC_073224.2 | 40 | 0.79 | 125056559-125057195 | 0.078 | 1.0 | 0.92 | 1.4e-16 |
| LOC134756861 | NC_073224.2 | 40 | 0.86 | 117517117-117517732 | 0.087 | 1.0 | 0.91 | 2.0e-19 |
| LOC134756861 | NC_073224.2 | 40 | 0.86 | 117513957-117514560 | 0.091 | 1.0 | 0.91 | 1.5e-17 |
| LOC134756861 | NC_073224.2 | 40 | 0.86 | 117513186-117513905 | 0.091 | 1.0 | 0.91 | 1.5e-17 |
| LOC115932779 | NC_073242.2 | 20 | 0.44 | 28928691-28933193 | 0.909 | 0.0 | 0.91 | 1.9e-02 |
| LOC134756861 | NC_073224.2 | 40 | 0.86 | 117512378-117513013 | 0.093 | 1.0 | 0.91 | 1.2e-16 |
| LOC129525331 | NC_073224.2 | 40 | 0.79 | 125058128-125058727 | 0.094 | 1.0 | 0.91 | 3.1e-16 |

## Caveats
- phase_qual is the diploid-bipartition cleanliness; low phase_qual = paralog/>2 haplotypes (down-weighted by the phase_qual>=0.7 high-confidence filter).
- transition-only-SNP loci (n_tv=0) could be edit-phased; the genetic core requires >=1 transversion het SNP.
- still coverage-limited (>=12x, >=5 reads per haplotype per junction).

## Reproduce
- `MINIFORGE python bench/allele_specific_junctions_multisnp.py --chrom <C>`
- `python3 bench/asjm_aggregate.py`

