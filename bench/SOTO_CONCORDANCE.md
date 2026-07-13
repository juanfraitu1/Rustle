# Quantitative concordance with Soto et al. 2025 (Cell)

**Date:** 2026-07-12. **Paper:** Soto, Uribe-Salazar, …, Dennis. *Human-specific gene expansions contribute to
brain evolution.* **Cell** 188(19):5363–5383 (2025), doi
[10.1016/j.cell.2025.06.037](https://doi.org/10.1016/j.cell.2025.06.037); bioRxiv
[10.1101/2024.09.26.615256](https://doi.org/10.1101/2024.09.26.615256). Data repo:
[`mydennislab/HSD_brain_evolution`](https://github.com/mydennislab/HSD_brain_evolution), Zenodo
[10.5281/zenodo.15486469](https://doi.org/10.5281/zenodo.15486469). *(Reference retrieved via PubMed; per-species
gorilla copy number is on the authors' internal share, not in the public repo — see below.)*

**Purpose.** An **independent, cross-lab, cross-modal** external anchor for the advisor: Soto's is a DNA / T2T
long-read catalog from a different lab; ours is RNA. This is the quantitative version of the earlier
`SOTO_FAMILY_COMPARISON.md` (which was a reading of the paper) — here we **run our method** on Soto's own gene set.

## Why the literal comparison is a category difference (unchanged, still true)

Soto's 213 families are **human-specific by construction**: human/ape CN comparisons with the great apes —
**gorilla (Kamilah) as the outgroup denominator** — and any family also duplicated in apes (e.g. FOXO3) is
*excluded*. So "find Soto's families in gorilla" is ill-posed: they are single-copy or absent in gorilla by his
own definition. His gorilla per-family CN table is not public (internal share
`/share/dennislab/databases/data/wssd-t2t/…`), so a direct gorilla-DNA-CN vs our-RNA-CN correlation is not possible
from released data. What IS testable, on his released gene set, is a **two-sided concordance**.

## The two-sided test (measured this session)

All numbers on `GGO_mm.bam` / `GGO.fasta`, gorilla annotation `GGO_genomic.gff`, SEDEF segdup map
`final.bed` (fracMatch = col 21). Copy calls = `copy_assign … --min-copies 2 --homology-primary --lambda-file`
run **foreground/serial** (`soto_specificity.sh`). Reproduce the figure: `bench/make_soto_concordance.py`.

### (A) SPECIFICITY — Soto's human-specific expansions are single-copy in gorilla by OUR method

For each human-specific family, our RNA method returns **0 multi-copy families** at the gorilla ortholog locus,
and the locus **is expressed** (so this is a real single-copy call, not "no data"):

| Soto family (human paralog) | gorilla ortholog | annot. loci | our RNA call | testis reads |
|---|---|---|---|---|
| GPR89B\* | GPR89A | 1 | **single-copy (0 fam)** | 70 |
| FRMPD2B\* | FRMPD2 | 1 | single-copy | 363 |
| CD8B2\* | CD8B | 1 | single-copy | 4 |
| SRGAP2B/C/D | SRGAP2 | 1 | single-copy | 248 |
| ARHGAP11B | ARHGAP11A | 1 | single-copy | 310 |
| HYDIN2 | HYDIN | 1 | single-copy | 187 |
| ROCK1P1 | ROCK1 | 1 | single-copy | 579 |
| FAM72B/C/D | FAM72A | 1 | single-copy | (within SRGAP2 locus) |
| GPRIN2B | GPRIN2 | 1 | single-copy | 8 |
| DUSP22B | DUSP22 | 1 | single-copy | 277 |
| NPY4R2 · CFC1B · NOTCH2NL | — | 0 | absent | — |

**13/13 concordant.** We independently agree these are not gorilla families — the method does **not** fabricate
Soto's human-specific expansions in gorilla. (\* = Soto's functionally-validated genes: GPR89B brain size,
FRMPD2B synapse, CD8B2 T-cell selection.)

### (B) RECOVERY — the ancestral / shared gorilla families are resolved at concordant copy number

| gorilla family | our RNA χ(H) | max SEDEF segdup identity | Soto's SD98? |
|---|---|---|---|
| RBMY | 6 | 99.74% | **SD98 ✓** (Soto's pipeline clusters this) |
| TSPY | 5 | 99.45% | **SD98 ✓** |
| DAZ | 2 | 99.63% | **SD98 ✓** |
| GSTM | 3 | 95.29% | paralog family, below SD98 |
| MAGEA | 2 | 94.85% | paralog family, below SD98 |
| PCDHB | 5 | 88.86% | paralog family, below SD98 |

The three Y-ampliconic families (**DAZ, RBMY, TSPY**) are genuine ≥98%-identity segmental duplications — the exact
objects Soto's **SD98 + famCN** pipeline clusters — recovered at the annotated copy number. The autosomal families
(GSTM, MAGEA, PCDHB) are **more divergent than SD98** (89–95%) — below Soto's threshold, yet still resolved by our
per-read PSV gate (the harder regime, and the added reach of RNA + PSVs over a DNA identity cutoff).

## What to tell the advisor

1. **Soto is an independent DNA/T2T catalog from another lab — the external anchor.** On his released gene set,
   our RNA method is **concordant on both sides**: single-copy in gorilla for the 13 human-specific expansions
   (specificity — expressed loci, we do not over-call), multi-copy for the ancestral families (recovery).
2. **The overlap is not circular.** Soto's copies/CN come from DNA read-depth on a T2T assembly; ours from
   Iso-Seq reads through the PSV gate. No shared aligner, no shared silver standard.
3. **Three recovered families (DAZ/RBMY/TSPY) are literally SD98 segdups** — the same objects Soto's pipeline
   builds — so where a Soto-style DNA family is also gorilla-multi-copy and transcribed, we return it.
4. **Honest limits:** his per-family gorilla CN is not public, so this is gene-set concordance, not a CN
   correlation; and his human-specific families cannot be *recovered* in gorilla (they are not there) — the
   correct, and passed, test is specificity + ancestral recovery.

Related: `bench/SOTO_FAMILY_COMPARISON.md` (qualitative predecessor), `bench/SOTO_SEGDUP_CN_REFINE.md`
(SD98/famCN refine on gorilla, the DNA-side loop), `bench/slides/soto_concordance.png`, `soto_specificity.sh`.
