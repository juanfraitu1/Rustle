# Can we recover Soto 2025's families? Why the comparison is a category difference, and what the right test is

**Date:** 2026-07-10. For the advisor's question: "Soto used DNA + segmental duplications in human data, but it
should work fine here — can we find the families Soto mentions?" **Answer: Soto's families cannot exist in gorilla
by Soto's own definition, so the literal comparison is a category error. But the underlying question — "does the
method recover real gorilla multi-copy families?" — has a clean positive answer, on a set Soto's method would
itself select.** Everything below is quoted from `soto_2025.md` or measured on `GGO_mm.bam` / the gorilla annotation.

## 1. Soto's "family" is a different object from ours — DNA, not RNA

Soto's operational definition (STAR Methods, "Gene family clustering"):

> "SD98 genes were grouped into gene families based on shared exons … refined using … whole-genome shotgun
> sequence detection (WSSD) (famCN) … groupings where the mean absolute deviation of the CN was less than one
> were selected."

and SD98 itself:

> "SD98 regions were defined as an SD with ≥ 98% sequence identity to another locus in the T2T-CHM13 genome."

So a Soto family is a **DNA segmental-duplication object**: ≥98% genomic-identity duplicated segments, clustered by
shared exons and confirmed by read-depth copy number (famCN/WSSD). It is our **E_a (DNA) axis**, not E_r (transcribed
homology), E_c (read-conflict), or E_p (protein). There is no expression requirement, no clique/χ(H) formalism, no
copy-assignment. **Expression is not required for membership** — the paper states only **455 / 1,002 (45%)** of the
paralogs are even brain-expressed (TPM ≥ 1), and unprocessed pseudogenes count toward forming a family. Our object
is the opposite: a family exists for us **only if it is transcribed**.

## 2. Soto's families are HUMAN-SPECIFIC by construction — so they cannot be gorilla families

The paper's thesis is **"213 human-specific gene families."** Human-specificity is established by using the great
apes — gorilla among them — as the **outgroup denominator**:

> "Human-specific and -expanded gene families were identified using CN comparisons between humans and nonhuman
> great apes with previously published WSSD (famCN) CNs from humans (SGDP n = 269) and four nonhuman great apes,
> including … gorilla (Kamilah) …"

Gorilla appears in the entire paper **exactly twice**, both as a comparator: Kamilah's read depth is one of four
ape values used to threshold "duplicated in humans vs not," and gorGor6 is one reference in a dN/dS alignment.
**There is no gorilla-specific family, count, or result anywhere in Soto 2025.**

The decisive detail is the paper's own **exclusion rule**. FOXO3 — CN-constrained and brain-expressed — was
**removed** from the human-expansion list precisely because it is "also duplicated in other great apes at similar
CN." So **any family that is a family in gorilla was deliberately excluded from Soto's list.** Soto's 213 families
and "gorilla multi-copy families" are therefore **disjoint by construction.**

This is confirmed in the gorilla annotation. Every named Soto family is single-copy or absent:

| Soto family | human status | gorilla annotation | testis reads |
|---|---|---|---|
| SRGAP2 (SRGAP2C) | human-specific | **1** copy (ancestral SRGAP2 only) | 195 |
| ARHGAP11 (ARHGAP11B) | human-specific | **1** (ARHGAP11A only) | — |
| NOTCH2NL | human-specific | **0** | 0 |
| NPY4R (NPY4R2) | human-specific | **0** | 0 |
| FAM72 | human-specific | **1** | — |
| CD8B (CD8B2) | human-specific | **1** | 4 |
| FRMPD2 (FRMPD2B, ~2.3 mya) | human-specific | **1** | 12 |
| GPRIN2 | human-specific | **1** | 5 |
| HYDIN | expanded | **1** | 108 |

SRGAP2C, ARHGAP11B, NOTCH2NL do not exist in gorilla — the gorilla has only the ancestral single copy. **A
single-copy gene is not a family**: there is nothing to duplicate-assign. Asking our method to "find SRGAP2 as a
family in gorilla" is asking it to resolve copies of a gene that has one copy. The negative result is a fact about
gorilla biology (and Soto's definition), not about the method.

## 3. The short-read limitation Soto describes IS the rationale for our method

Soto's argument for why this is hard with short reads is, verbatim, the collapsed-copy problem our
copy-assignment gate is built for:

> "paralog-specific variants are mistaken for SNPs when reads from both paralogs map to a single collapsed locus,
> resulting in false mid-frequency alleles."

with the sensitivity numbers: short-read SNV sensitivity **0.85% in SD98** vs 87.6% in non-SD, because only 10.86%
of SD98 is short-read-accessible. **This is precisely the K=0 / collapsed-copy wall we characterised**, and Soto's
conclusion — that long reads are needed to separate paralog-specific variants from errors — is exactly our
per-read significance gate (PSV vs sequencing error), demonstrated on the DAZ2 recovery. Soto argues for long-read
paralog resolution on DNA; we do it on RNA. **The paper is evidence for our approach, not against it.**

## 4. The right test — and it is positive

The correct gorilla analogue of "a Soto family" is: a locus that **is multi-copy in gorilla** and **is
transcribed in our tissue**. That is the set our audits actually cover, and we recover it:

| gorilla testis multi-copy family | annotated | recovered |
|---|---|---|
| GSTM (globin-problem, uniquely-mapping copies) | 3 | **3/3** |
| MAGEA (inverted duplicate) | 2 | **2/2** |
| DAZ (collapsed, recovered this session) | 2 | **2/2** |
| RBMY (Y ampliconic, distal cluster) | 6 | **6/6** |
| TSPY (Y ampliconic array) | 6 | **5/6** (34 reads) |
| PCDHB, APOBEC, APOL, GBP | 2–6 | recovered |

Notably, **Soto's own method would place the Y-ampliconic families (DAZ, RBMY, TSPY) and MAGE among SD98 segmental
duplications** — they are ≥98%-identity duplications — so where a Soto-style DNA family is also gorilla-multi-copy
and testis-expressed, we recover it. The families we cannot recover are the ones that (a) are human-specific
(don't exist in gorilla), (b) are brain-only (not in testis), or (c) are silent pseudogene copies (not transcribed
by definition) — three properties that follow from Soto's DNA-and-brain design, not from a limit of the method.

## What to tell the advisor

1. **Soto's family is a DNA segmental-duplication object; ours is a transcribed-homology object.** Different axes.
   His includes silent pseudogenes and requires no expression; ours requires transcription. This is stated in his
   own Methods (SD98 + famCN, 45% expressed).
2. **His 213 families are human-specific — defined by using gorilla as the outgroup.** By his own exclusion of
   ape-shared expansions (FOXO3), his list is disjoint from gorilla families by construction. Every named family
   is single-copy or absent in the gorilla annotation. So the families cannot be found in gorilla — not because
   the method fails, but because they are not there.
3. **The short-read failure Soto quantifies (0.85% SD98 sensitivity, paralog collapse → false mid-frequency
   alleles) is exactly the problem our long-read copy-assignment gate solves.** His paper motivates our method.
4. **On the correct test — gorilla multi-copy, testis-expressed families — we recover them** (GSTM, MAGEA, DAZ,
   RBMY, TSPY, PCDHB, …), and several of those are SD98 segmental duplications Soto's own pipeline would cluster.

If the advisor wants the direct methodological equivalence: run **Soto's pipeline (SD98 + famCN) on the gorilla
genome** to get the gorilla DNA family catalog, then show our RNA families are the transcribed subset of it. That
is a DNA-side experiment (SEDEF/BISER segdup calling on the gorilla assembly, which the project has scaffolding
for), not an RNA one — and it is the honest way to close the DNA↔RNA loop he is pointing at.

Related: `reference_soto_2025_hsd_brain`, `reference_bailey_2002_segdups` (WSSD/famCN ancestor), `bench/YAG_VS_ISOCON.md`,
`bench/DAZ2_RECOVERED.md`, `bench/FAMILY_ARTIFACT_AUDIT.md`.
