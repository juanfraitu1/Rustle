# Protein-level validation — the orthogonal ground truth for family definition

The cDNA-homology over-merge could never be adjudicated from the inside (0 trustworthy positives in
the big blobs). Protein homology is the orthogonal validator the project flagged as needed. Built from
`proteins.fa` (22,614 representative proteins, gene-name keyed — directly comparable to the cDNA
graph) + mmseqs all-vs-all (`prot_ava.m8`, 982,949 hits). Protein paralog edge = fident≥0.30 ∧
reciprocal-cov≥0.50 ∧ e≤1e-5 → **199,828 protein-paralog gene-pairs** = the positives we lacked.
Scripts: `family_def_protein_validate.py`, manifest at `protein_family_manifest.tsv`.

## 1. The over-merge is CONFIRMED and EXPLAINED

| cDNA component | n | coding % | non-coding | distinct protein orthogroups |
|---|--:|--:|---|--:|
| comp-238 (the resistant blob) | 238 | 38% | **62% lncRNA** | **10** (+148 lncRNA) |
| comp-201 | 201 | 71% | 29% lncRNA | 17 |
| comp-140 | 140 | 51% | 49% lncRNA | 10 |
| comp-47 | 47 | **100%** | — | **1** (a REAL family) |
| comp-41 | 41 | **100%** | — | **1** (a REAL family) |

The dense comp-238 blob that resisted every junction/topology lever is a **non-coding repeat
over-merge** — 62% lncRNA sharing a repeat element, plus 90 coding genes that form 10 *distinct*
protein families. comp-47 and comp-41 are each a single protein orthogroup = genuine large real
families. Big over-merge components are 30% lncRNA vs 16% in small families.

## 2. CORRECTION — co-threading beats cov_min (only visible with ground truth)

Recovering the protein-orthogroup partition inside the over-merged components (pairwise prec/recall):

| de-bridge weight | precision | recall | F1 (bar .30/.50) | F1 (bar .40/.60) |
|---|--:|--:|--:|--:|
| none (raw over-merge) | 0.28 | 1.00 | 0.44 | 0.42 |
| cov_min | 0.70 | 0.60 | 0.64 | 0.62 |
| **frac_ref (co-threading)** | 0.71 | 0.74 | **0.72** | **0.70** |
| product (cov_min×frac_ref) | 0.64 | 0.56 | 0.60 | 0.58 |

**co-threading recovers real families better than cov_min** — cov_min's smaller max-component (43 vs
57) was OVER-FRAGMENTATION (recall 0.60: it splits real protein families), a failure mode invisible
without ground truth. Default reverted to `frac_ref`. (This corrects the prior call to switch to
cov_min, which was based on max-comp reduction alone.)

## 3. DELIVERABLE — protein-confirmed family definition

Keep cDNA-homology edges that are ALSO protein-homologous (coding ∧ reciprocal-cov≥0.50):

- **837 families** (vs cDNA over-merged 1216), **max comp 238 → 47** (a real orthogroup), 2,944 genes.
- 58% of cDNA edges are protein-confirmed; **4/5 named real families** preserved.
- Dropped 1,775 genes with no protein-confirmed paralog: 772 lncRNA (non-coding repeat bridges,
  correctly excluded), 845 protein_coding (cDNA bridges in non-coding regions / sub-bar), and — as a
  KNOWN limitation — 103 V_segment + 13 C_region (immunoglobulin/TCR segments are short and fail the
  reciprocal-cov bar) + 33 rRNA.

This is the defensible "real protein-coding family" claim: members are coding, expressed (cDNA
homology), AND whole-protein homologous — and it breaks the over-merge to a real-family max-comp on
PRINCIPLE, not heuristic, with trustworthy positives.

## 4. Known limitations (separately handleable)
- **Immune-receptor families** (IG/TCR V/C/D/J segments) are missed — short peptides below the
  reciprocal-cov bar. Need a segment-aware bar.
- **Partial paralogs** (GGT1~GGTLC2) are missed — protein homology sub-bar (GGTLC2 covers only part
  of GGT1). The recurring partial-duplication boundary.
- **Non-coding (lncRNA) families** are a separate phenomenon (repeat-bridged) — not protein families;
  handle on the cDNA axis if of interest.

## Status / recommendation
- Protein orthogroups = the family ground truth for coding genes; use them directly.
- The cDNA de-bridge (co-threading `frac_ref`, NOT cov_min) is the proxy where protein is unavailable
  (F1 0.72 vs protein truth).
- lncRNA repeat over-merges are the bulk of the residual blob; they are not protein families.
- Next: segment-aware bar for IG/TCR; decide the partial-paralog policy; optionally a protein-level
  community detection for any protein orthogroup that itself over-merges via shared domains.

---

## ADDENDUM — segment-aware recovery (Task 1) + protein community detection (Task 2)

Script: `family_def_immune_and_protein.py`.

### Task 1 — immune-receptor families recovered (locus constraint)
The 116 IG/TCR segments (V/C) have **0 representative proteins** (not translated until VDJ
recombination), so no protein bar can reach them; and their shared V-DOMAIN over-merges segments
ACROSS loci (one naive cDNA component spanned 6 chromosomes). Fix: segment-biotype cDNA homology +
**SAME-LOCUS constraint** (same chrom, ≤3 Mb) — the inverse of the disjoint-loci rule for dispersed
paralogs. Result: **30 immune families, all 30 single-chromosome compact tandem arrays** (0.06–0.97 Mb
at the IGH/IGK/IGL/TRA/TRB loci). Recovers a whole real-family class the protein definition dropped.

### Task 2 — protein orthogroups DO domain-over-merge; the cDNA bar already handles it
The *pure* protein-orthogroup graph (fident≥0.30) over-merges ancient domain families: ZNF (569,
density 0.34), OR7D (501), GPCR, HOX, kinases — **10/78 large orthogroups** split under protein-recip-
cov-weighted Louvain. BUT the cDNA∩protein deliverable is already clean, because the strict cDNA
id≥0.90 bar splits each ancient domain blob into its RECENT sub-duplications:

| family | protein-alone | cDNA∩protein (recent paralogs) |
|---|---|---|
| ZNF | one 569-gene blob | 43 genes → families of 4,3,2,2,… |
| OR / GPCR | 501 / 231 blobs | 2 genes (the rest are cDNA-divergent ancient) |
| KRT | 78 blob | 9 genes → families of 2,2,2 |

So the two definitions capture different things and that is the right design:
- **cDNA∩protein = RECENT coding paralog families** (cDNA still ≥0.90) — clean, no domain over-merge,
  and exactly the recent-duplication / copy-assignment regime the thesis targets. 837 families, max 47.
- **protein-alone = all homology incl. ancient domain families** — over-merges via domains; needs
  community detection (the Task 2 tool) if ancient families are wanted.

## CONSOLIDATED family definition (layered)
1. **Recent coding paralog families**: cDNA-homology (id≥0.90) ∩ protein-homology (recip-cov≥0.50) →
   **837 families**, max comp 47, retrocopy- and over-merge-free.
2. **Immune-receptor families**: segment biotype + cDNA homology + same-locus → **30 families**.
3. **Ancient domain families** (ZNF/OR/GPCR/HOX): protein-level only; appear as many recent
   sub-families in (1); recoverable as super-families via protein-graph community detection if wanted.

Total **867** trustworthy families, each with an explicit, biotype-appropriate criterion — replacing
the 1,216-family / 238-blob cDNA over-merge.
