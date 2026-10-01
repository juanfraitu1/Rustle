# Pre-registration: does Soto's copy-number cut follow duplicon boundaries? (2026-09-30, KEY=cnduplicon)

Written before any statistic below was computed. Species: human CHM13 v1.0 only. Nothing in `src/` or `tools/` changes.
Script: `bench/soto_m2/soto_m2_duplicons.py` (written after this file is committed; the commit hash of this file is
recorded in the result section).

## 1. Question

Soto's families are sequence families (genes sharing ≥ 98%-identical exons) cut where copy numbers differ by 2 or more
(`docs/SOTO_M2_MEETING_EVIDENCE_2026-09-30.md`; 33 sequence families hold 87 Soto families). Segmental duplications are
mosaics of duplicons (ancestral duplication units, DupMasker). **Does Soto's cut group genes whose exons sit on the same
duplicons?** If yes, the copy-number gate behaves like a duplicon split, and Soto's families reconcile with duplicons as
units. If no, the cut is finer than duplicon structure.

## 2. Inputs (frozen)

- Duplicons: `chm13.draft_v1.0_plus38Y_dupmasker_colors.bed` from Vollger et al. 2022, Zenodo 4726156 (`data.tar`, member
  fetched by byte range, 13,114,022 bytes; the source Soto cite for DupMasker). Coordinates CHM13 v1.0. Segments may overlap
  (48,745 overlapping neighbours): every overlapping duplicon counts.
- Genes and exons: `sd98_gene_exons.tsv` (CAT v4), Soto's table S1C (`bench/soto/soto_famCN_S1C.tsv`), and the exon map-back
  edges `bench/soto/shared_exons_5154_exon_mapback.tsv`.
- Sequence families: `soto_replication.py`'s `pair_families(..., gate=False)` + `collapse_cover`, exactly as in `nesting`.
- Split: `bench/soto/soto_split_2026-09-29.tsv` (frozen gene → dev / heldout).

## 3. Units and statistic (fixed now, no free constant)

- **Clusters:** sequence families that hold ≥ 2 Soto families with ≥ 2 clean members each (expected 33 clusters / 87 Soto
  families; the script asserts it). **Genes:** the clean members (genes Soto places in one family) of those Soto families.
- **Duplicon composition of a gene:** for every duplicon ID, the number of the gene's exonic bases it covers.
- **Similarity of two genes:** weighted Jaccard of their compositions, Σ min / Σ max over duplicon IDs (0 if either is empty).
- **Per-cluster statistic Δ:** mean similarity of gene pairs inside one Soto family minus mean similarity of gene pairs in two
  different Soto families, over the cluster's genes.
- **Pooled statistic:** the mean of Δ over clusters (each cluster weight 1).
- **Null:** shuffle Soto family labels among each cluster's genes (family sizes kept), all clusters at once; 10,000 shuffles,
  seed 20260930. One-sided p = (1 + #{null ≥ observed}) / 10,001. Per-cluster p the same way.

## 4. Decision rule (fixed now)

- **FOLLOWS:** pooled p < 0.01 and Δ > 0 in at least 2/3 of clusters.
- **PARTIAL:** pooled p < 0.05 otherwise.
- **DOES NOT FOLLOW:** pooled p ≥ 0.05.
- **Held out:** each cluster is assigned to the split half of its largest Soto family (ties: dev). The rule is applied to
  each half separately; nothing is tuned on dev (there is nothing to tune). **The verdict is the held-out half's**; if the two
  halves disagree, the verdict is reported as SPLIT with both.

## 5. Secondary (descriptive, no verdict)

- **S1C SD Unit arm:** the same test with Soto's own per-gene `SD Unit` labels as a set (genes labelled `.` dropped).
- **Copy number vs duplicons:** Spearman between |ΔfamCN| (S1C) and 1 − similarity over (a) all gene pairs inside the 33
  clusters and (b) pairs across a Soto boundary only.
- **Genes with no exonic duplicon:** counted; the primary test is rerun without them.

## 6. Seen before writing this file (disclosed)

- 2,259 of 2,334 S1C genes have exonic duplicon overlap; the dominant exonic duplicon is among S1C's `SD Unit` labels for
  1,055 of 1,402 labelled genes.
- Six genes looked at by hand: NPIPB3 and NPIPB4 share their dominant duplicon (SD9443) while NPIPB5's is SD9622; all three are
  in different Soto families. NPIPA1's dominant is SD9613. This hints against FOLLOWS inside NPIP; no other cluster was looked at.
- In S1C, the `SD Unit` values of all six NPIP-side Soto families come from their shared lncRNA members (SD9449, SD9479,
  SD9481); the NPIP genes themselves are labelled `.`.

## 7. Result

(Filled in after the run, below this line, without editing anything above.)
