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

## Result
- **7,898** phaseable het genes; **74,674** alternatively-spliced junctions tested.
- **475 allele-specific junctions (FDR q<0.05, |ΔPSI|≥0.3) across 235 genes.**
- Strong effects: median |ΔPSI| **0.64**; **213** ≥0.7; **56** full switches (ΔPSI=1.0).
- **120** have a **transversion** anchor → unambiguously **genetic** (not RNA-editing); 59 genes.

## Verification (deterministic + mechanistic — stronger than an LLM pass here)
- **Collapsed-paralog masquerade ruled out:** all 475 are **uniquely mapped** (frac_mq0<0.1) — the
  anchors are genuine heterozygous sites, not paralog copies. (0 collapsed.)
- **Editing controlled:** transversion anchors (120) cannot be RNA-edit sites; transitions flagged.
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
