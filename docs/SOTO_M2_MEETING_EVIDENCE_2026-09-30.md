# Soto vs our families: reproduced and extended (2026-09-30)

These are the second-machine tasks of `docs/HANDOFF_SECOND_MACHINE_2026-09-30.md`, run in an isolated clone on this
machine (`/mnt/linuxdisk/home/juanfraitu/rustle_m2`, branch `machine2/soto-evidence`). They extend
`docs/SOTO_VS_OURS_MEETING_2026-09-30.md`. Soto's Table S1C is the truth throughout. Species: human CHM13 v1.0.

## 0. Reproduction: both commands match exactly

- `nesting`: 444 Soto families; sequence only 440/444 (99.1%) and 394/398 (99.0%); 33 clusters hold 87 Soto families;
  exceptions ID_192, ID_347, ID_482, ID_62. **Match.**
- `ladder` (ALL ARI / exact): 0.7307/345, 0.9198/373, 0.8855/375, 0.9277/411, 0.9698/479 (held-out 0.9681/263). **Match.**

## 1. Objective 1: Soto's families sit inside ours (task A2 adds a condition)

| edge set | Soto families inside one of our sequence clusters | our clusters that are unions of whole Soto families |
|---|---|---|
| exon map-back (the edge set of Soto's released code; gives ARI 0.97) | **440/444 (99.1%)** | 394/398 (99.0%) |
| SEDEF-projected exons (their Methods prose, literal; gives ARI 0.71) | 346/444 (77.9%) | 287/446 (64.3%) |

**Condition:** the 99.1% holds on the exon map-back graph. On the sparser SEDEF-projected edges, 98 Soto families span two or
more components (for example AMY1 vs AMY2A in ID_131, and USP17L4 in ID_107). There, Soto joins genes that those edges do not
join, so the difference comes from missing edges, not from a conflict. ⚠ "Our clusters" here are components of the sequence
graph built with Soto's own ≥ 98% shared-exon rule, not the Rustle O1 catalog. For the pipeline catalog, the numbers are
in `bench/SOTO_AS_A_REFINEMENT.md`: L1/L2 nest 90%/80%, all 6 NPIP-side Soto families at L1/L2, and 80.3%/76.4% held out on chr5/7/21.

## 2. Objective 2: under the right conditions we find what they find

Unchanged from the meeting sheet (reproduced above): with Soto's released-code choices and **our own** copy numbers
(268 SGDP samples, their gene-body ∩ SD98 interval), ARI 0.9277 with 411/491 families exact (held-out 0.9343). With their
published copy numbers, ARI 0.9698 with 479/491 exact. Without FAM90A (ID_356) these are 0.925 and 0.965.

## 3. Objective 3: Soto's families are smaller than homology and incomplete

`python3 bench/soto_m2/soto_m2_audit.py narrower ...` (B3):

1. **Single-gene "families" of near-identical paralogs.** In 47 Soto families, only one member belongs to that family alone;
   the other members are genes Soto also places in other families. For 40 of the 47 (85.1%), that lone gene shares
   ≥ 98%-identical exons with an exclusive member of another Soto family. The clearest case is NPIPB3, NPIPB4 and NPIPB5.
   Each is its own family (ID_151, ID_152, ID_153), and each family is padded with the same five shared genes
   (AC126755.6, AP001120.2, MSTRG.2119, PDXDC2P-NPIPB14P, PKD1P6-NPIPP1). NPIPB3 and NPIPB4 are 98.6% identical. Other lone
   genes: GOLGA8B (ID_79), GOLGA6L3 (ID_89), USP17L8 (ID_108), PCMTD2 (ID_183), NIPA1 (ID_189), ANKRD20A4P (ID_281).
2. **The copy-number gaps behind those splits do not replicate.** Soto separates NPIPB3/B4/B5 at famCN 26.5/29.1/22.2 (gaps of 2.6
   to 6.9 copies). Our re-measurement over the same interval reads 31.7/36.9/33.4, which groups them differently. Across
   all 1,793 genes the two measurements agree closely (median relative difference 0.8-2.3%). However, at famCN 20-35, **21%** of genes
   differ by ≥ 2.6 copies, which is the size of the gap that splits NPIPB3 from NPIPB4. NPIPA1 reads 9.0 in S1C against 26.4 in ours.
3. **23% of Soto's multi-gene families are strictly smaller than their sequence component** (102/444). Together they leave out
   1,336 gene memberships (median 6 per family). Inside our clusters, 971 distinct ≥ 98% exon-sharing pairs join two different
   Soto families (median copy-number gap 14.0 S1C, 13.6 ours). The split is by copy number, as designed.
4. **83 real paralog pairs missed** (09-28 attribution, reproduced from `attributed_pairs.tsv`): 65 have no edge in their
   graph (median identity 0.839; GOLGA6/GOLGA8 on chr15 and NPIP A/B clades on chr16), and 18 are split by copy number (median
   identity 0.980, all NPIP). The most identical pair kept apart is NPIPB12/NPIPB13 at 0.991 (ID_154 vs ID_155).
5. **Pseudogene composition (B1)** is a pure tabulation of S1C's `Biotype` column. It matches S1C's own per-family `No. Protein
   Coding` column (0 mismatches):

   | counting | families | entirely pseudogene | no protein-coding member | pseudogene members |
   |---|---|---|---|---|
   | 491 real families (the 114 `Unassigned_*` singletons excluded) | 491 | **183 (37.3%)** | 253 (51.5%) | 1,420/2,458 (57.8%) |
   | as in the 09-10 audit (singletons counted) | 605 | 217 (35.9%) | 287 (47.4%) | 1,454/2,572 (56.5%) |

   The 09-10 numbers reproduce exactly **only** when the `Unassigned_*` singletons are counted as families. Quote the 491 row.
6. **Fragments next to full-length members (B2).** Using merged exonic bp from `sd98_gene_exons.tsv`, all 491 families have ≥ 2
   members with coordinates. 154/491 (31.4%) hold a member < 20% the size of their largest member. 45/491 (9.2%) also hold a
   second member ≥ 80%, meaning two full copies plus a fragment. Worst ratios: ID_280 FAM153CP / CR392039.2 126×, ID_215 PDE4DIP /
   AC239860.2 119×, and ID_150 NOMO3 / MIR3179-1 110×. ⚠ The 09-10 figure "147 of 420 (35.0%)" does **not** reproduce: its script
   is gone, and its 420 denominator cannot be rebuilt. Quote 154/491 instead.

## 4. What can be said

- **Supported:** Soto's families are copy-number refinements of sequence families (99.1% nesting on their released-code graph).
  With our own copy numbers we rebuild them to 0.93 ARI and 411/491 exact. Where they are finer, the cuts rest on copy-number gaps
  that a second measurement does not always reproduce (NPIPB3/4/5), and they leave near-identical paralogs in single-gene families.
- **Not supported: "ours is better" as an accuracy claim.** Soto's table is the truth in every number here, and the two objects
  answer different questions (homology groups vs copy-number-coherent expansion groups). The defensible form is an
  **information** claim. Soto = ours + a copy-number cut: we can derive theirs from ours (0.93 using our own copy numbers). Ours cannot be derived
  from theirs, because their table records no links between families. The information about which Soto families belong together
  (971 cross-family pairs sharing ≥ 98%-identical exons, and 83 paralog pairs at median identity 0.858) is lost.
- The 99.1% needs its condition stated: on the literal-Methods edges it is 77.9%, and on our Rustle O1 catalogs it is 80-90%.

## 5. Regenerate

```
python3 bench/soto/soto_replication.py genesets --out-eligible elig.tsv --out-full full.tsv
python3 bench/soto/soto_replication.py nesting --shared bench/soto/shared_exons_5154_exon_mapback.tsv \
    --geneset elig.tsv --full-geneset full.tsv --famcn-ours famcn_ours_allwssd.tsv
python3 bench/soto/soto_replication.py nesting --shared bench/soto/shared_exons_2334_finalv1_native.tsv \
    --geneset elig.tsv --full-geneset full.tsv --famcn-ours famcn_ours_allwssd.tsv            # A2
python3 bench/soto_m2/soto_m2_audit.py biotype                                                # B1
python3 bench/soto_m2/soto_m2_audit.py fragments --exons sd98_gene_exons.tsv                  # B2
python3 bench/soto_m2/soto_m2_audit.py narrower --geneset elig.tsv --full-geneset full.tsv \
    --famcn-ours famcn_ours_allwssd.tsv --pairs attributed_pairs.tsv                          # B3
```

Inputs outside the repository: `famcn_ours_allwssd.tsv` and `sd98_gene_exons.tsv`
(`/mnt/linuxdisk/home/juanfraitu/winloci_data/soto_replication/`), and `attributed_pairs.tsv`
(09-28 scratchpad `soto_attr/`, copied to `m2data/` in the clone). Each runs in under 2 s.
