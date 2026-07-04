# Applied at the DNA level, is the multi-copy family definition MORE SPECIFIC? — No. It's a recall/specificity trade.

**Question.** The four formal family objects share one skeleton (the same γ-quasi-clique refine + ≥2-loci
predicate) and differ only in the **edge oracle**. If we swap the RNA transcript-homology oracle (E_r) for
the DNA genome-segdup oracle (SEDEF `final.bed`, the catalog `genome_families_refined.tsv`), would the
resulting families be more specific?

**Answer: no — RNA is ~2.4× *more* specific; DNA is *more complete*. The two are a recall-vs-specificity
trade on two INCOMPARABLE homology oracles**, not a containment.

## Premise correction (load-bearing)

An earlier framing claimed `E_c ⊆ E_b(segdup)`, i.e. DNA-segdup ⊇ RNA-read-conflict → "DNA more complete."
**That is false.** The segdup oracle is **incomparable** to the read-conflict oracle — verified by
**APOBEC3D/F**: a read-conflict pair (∈ E_c) whose only-covering SEDEF pair is 88.4% < 90%, so ∉ strict
segdup. The segdup oracle *sees* copies RNA can't (silent/unexpressed SDs) *and misses* copies RNA has
(sub-90% expressed paralogs, retrocopies whose exons scatter across introns). So it was always a trade.

## Specificity — measured on a NON-CIRCULAR truth (protein E_p)

The DNA copy-number / segdup oracle is **circular** for judging a DNA definition (same genomic substrate),
so specificity is scored against **protein E_p** (reciprocal whole-protein mmseqs — a different molecule
from both the genomic DNA behind segdups and the reads behind RNA), applied identically to both catalogs
(mega-families PRFAM0 GPCR / PRFAM1 ZNF excluded).

| specificity (non-circular protein E_p) | DNA (segdup) | RNA (transcript) |
|---|---|---|
| edge precision, evaluable pairs | **0.358** | **0.849** ← 2.4× |
| gene-weighted dominant protein-family purity | 0.714 | 0.876 |
| max distinct protein families merged into ONE family | **22** | **3** |
| families chaining ≥5 distinct protein families | 26 | **0** |
| total false-merge pair mass | **75,193** | **902** ← ~83× smaller |

The one metric DNA "wins" (cDNA-evaluable precision 0.956 vs 0.865) is a **truth-file artifact** — that
cDNA edge file is a 94.5%-positive candidate list, not a negative-labelled test.

## Completeness — DNA wins on the silent axis, but not as a criterion

| completeness (vs cDNA-loose paralog set) | DNA | RNA |
|---|---|---|
| genomic paralog pairs captured | 11,726 (67%) | 418 (2.4%) ← DNA 28× |
| cDNA-truth genes covered | 3,510 (64%) | 592 (11%) |
| **recall on the shared *reachable* substrate** | 0.854 | **0.931** ← RNA |

The last row closes the circularity: **on the same expressed node set, RNA groups paralogs *better***. So
RNA's 28× lower absolute recall is a **transcriptome ceiling, not a criterion failure** — 96.1% of the
2,895 DNA-only-family genes never appear in RNA at all (GH1 pituitary, DRD5 brain, KIR2DL4 NK-cell,
ribosomal retrocopies).

## Why DNA is less specific — a segdup is a genomic feature, not a gene family

DNA's 75,193 false merges decompose into the classes RNA's whole-transcript `core_recip` filters out:

| DNA FP class | % of FP mass | mechanism | RNA has it? |
|---|---|---|---|
| **repeat-bridge / mega-blob** | **79.8%** (33 fam) | Alu/LCR transitive chaining (GRFAM0 = 189-gene 22q11 DiGeorge blob spanning 22 protein families) | **0 — E_r can't chain via Alu/intron** |
| non-coding tandem array | 15.0% (244 fam) | tRNA/snRNA/lncRNA arrays | ≈0 — spliced-protein oracle never assembles them |
| isolated single-domain share | 5.3% (152 fam) | distinct genes sharing one exon/domain (β-defensin+ZNF, LHB+NTF4) | yes, but ~13× fewer (307 pairs) |

## Bottom line for the thesis

DNA and RNA are **complementary, incomparable oracles.** DNA is the natural genomic copy-number **truth**
and more complete (it sees silent copies), **but a different *object*** — genomic duplications, not
transcribed gene families. RNA (E_r) is the more specific definition of the thing you actually study
transcriptomically, and it is **better on the shared substrate**. This is exactly why DNA belongs in the
pipeline as the **validation layer**, not the definition — moving the definition to DNA would ~83× the
false-merge mass to gain copies unobservable in RNA.

Scripts: `family_dna_vs_rna_headtohead.py`, `family_fp_mechanism.py` (+ `.json`). Deterministic
(`PYTHONHASHSEED=0`). Adversarially verified (truth non-circular, comparison fair, FP classes spot-checked).
