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

Run on 2026-09-30 after this file was committed (`bbbdda25`, sha1 of the file at that commit
`a9b5dca5f006d2157648307bc2668c92918fd4ed`); `bench/soto_m2/soto_m2_duplicons.py`, 56 s, light lock.

**VERDICT (§4): FOLLOWS** — on the held-out half and the dev half alike.

| set | clusters | pooled Δ | pooled p | clusters with Δ > 0 | verdict |
|---|---|---|---|---|---|
| held-out | 16 | +0.2853 | 0.0001 | 0.88 | FOLLOWS |
| dev | 17 | +0.3276 | 0.0001 | 0.88 | FOLLOWS |
| all (descriptive) | 33 | +0.3071 | 0.0001 | 0.88 | FOLLOWS |

Secondary: Soto's own `SD Unit` labels (28 clusters with ≥ 3 labelled genes) Δ +0.2468, p 0.0001, 0.68 → FOLLOWS; without the
5 genes lacking an exonic duplicon, Δ +0.3121, p 0.0001. Copy-number gap vs duplicon dissimilarity, Spearman +0.577 over all
8,565 gene pairs inside the clusters, +0.313 over the 2,834 pairs across a Soto boundary.

Clusters that do not follow (Δ ≤ 0.05 or p > 0.25): ID_69/ID_76 (−0.029), ID_163/ID_191 (−0.199), ID_96/ID_97 (−0.001),
ID_271/ID_272 (0.000), ID_172/ID_184 (+0.016), TBC1D3 ID_468/ID_469 (+0.049). Strongest: FAM90A/FAM86 ID_356/ID_355 (+0.818),
DUX4 cluster (+0.737), FGF7P (+0.621), FRG1 (+0.595). The NPIP cluster (+0.375, p 0.0001) includes SMG1P (ID_41, the adjacent
LCR16u module); the one-gene NPIPB3/B4/B5 families are outside this test (it needs ≥ 2 clean members per family).

**Reading.** Where Soto cut a sequence family by copy number, the pieces sit on different duplicons far more than chance, and
copy-number gaps grow with duplicon differences: different duplicons carry different copy numbers, so the copy-number gate
largely acts as a duplicon split. The exceptions (TBC1D3 among them) are cuts inside one duplicon composition. This supports
"duplicons as units" as the bridge between Soto's families and sequence families; it does not test the reverse (whether every
duplicon boundary is a family boundary).
