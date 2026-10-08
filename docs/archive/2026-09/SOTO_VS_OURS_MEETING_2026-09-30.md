# Soto 2025 and our families: what is measured (meeting sheet, 2026-09-30)

Every number below comes from a file in this repository and can be regenerated from a clone plus two small tables (§5).
Soto's family table (S1C) is the truth throughout, so every "agreement" is agreement with Soto, not an independent
measurement. Species: human CHM13 v1.0, Soto's own 2,334-gene universe.

## 1. Soto's families sit inside ours (Soto is a refinement of our sequence families)

Our sequence-only families are the connected components of the exon map-back graph (≥ 98% identical exons shared
between genes, Soto's own edge rule), with no copy-number test. Soto's families are those components split by copy number.

| our clusters built with | Soto families (444 with ≥ 2 clean genes) wholly inside one of ours | our clusters that are exact unions of whole Soto families |
|---|---|---|
| **sequence only** | **440 / 444 (99.1%)** | **394 / 398 (99.0%)** |
| + our own copy numbers (268 samples, Soto's interval), MAD < 1 | 401 / 444 (90.3%) | 393 / 456 (86.2%) |
| + Soto's published copy numbers, MAD < 1 | 433 / 444 (97.5%) | 433 / 450 (96.2%) |

Soto is finer than ours where they differ: of our 398 sequence clusters, 33 hold 87 Soto families (the largest holds 6; the NPIP cluster
holds 4, ID_154, ID_149, ID_155 and ID_41); the other clusters hold one Soto family or none.

The four exceptions are ID_347 (DUX4; marked "Manual merge" in Soto's table), ID_482 (UBTFL; also "Manual merge"), ID_192
(3 genes) and ID_62 (RPL23AP87 split off a 14-gene family).

Earlier, on our real pipeline catalogs (`bench/SOTO_AS_A_REFINEMENT.md`; small numbers, read with n): all six NPIP-side
Soto families lie inside our single NPIP family at levels L1 and L2 (6 / 6); of the 10 Soto families with ≥ 2 genes in that node
set, 9 (L1) and 8 (L2) nest; held out on chr5/7/21, 80.3% and 76.4% of 76 and 72 Soto families nest; Soto families of
> 3 genes nest less often than smaller ones (47.1% vs 89.8%).

## 2. Given the right conditions we find what they find

The same edges and rule, with more of Soto's own choices switched on (`soto_replication.py ladder`):

| recipe | copy numbers | ARI | exact families (of 491) |
|---|---|---|---|
| literal recipe (SD98 regions mapped back, exons projected; component split) | Soto's published | 0.7096 | 235 |
| released-code choices (exon map-back, per-pair test) | none (sequence only) | 0.7307 | 345 |
| released-code choices | our own (268 SGDP samples, their gene-body ∩ SD98 interval) | **0.9277** | **411** |
| released-code choices | Soto's published (Table S1C) | 0.9698 | 479 |

Conditions stated plainly: the step from 0.93 to 0.97 is agreement with their published copy numbers, not an independent
measurement; one family (FAM90A, ID_356) swings the ARI by about 0.035, so quote exact families beside ARI (without it: 0.925 and
0.965). Their released family-building loop, run exactly as released, scores ARI 0.82-0.87 (it depends on hash order); the 0.97 uses the
loop completed as it evidently intends. Independent recheck of the 0.97: CONFIRMED WITH CORRECTIONS
(`docs/archive/2026-09/PREREG_soto_reconciliation_2026-09-29.md`).

## 3. Where Soto's families are smaller or incomplete (by sequence homology) — the measured evidence

1. **Real paralog pairs Soto's families miss.** Of 247 pairs Soto keeps apart that our graph joins, 83 carry a direct sequence
   edge at median identity 0.858 (37 at ≥ 0.90, 30 at ≥ 0.95). Causes, per pair: 65 have no shared-exon
   edge in Soto's edge set (44 on chr15: GOLGA6 vs GOLGA8 and a GOLGA8A isolation; 21 on chr16: NPIP A-clade vs B-clade);
   18 (all chr16 NPIP) have an edge but are split by their copy-number grouping (NPIPA1 vs NPIPA7: copy number 9.05 vs 47.58).
   (`docs/archive/2026-09/SOTO_REPLICATION_STATUS_2026-09-28.md` §3; independently re-derived, 83 / 83.)
2. **Sequence-identical links that Soto's copy-number split separates.** About 700 cross-family exon links at median identity
   0.9948 join genes Soto places in different families (median copy-number gap 15.1) (register 1160).
3. **Their truth is a cover, not a partition.** 148 of 2,333 genes (6.4%) belong to two or more of Soto's own families.
4. **Their table bundles pseudogenes and fragments as gene families** (audit of S1C itself, `project_soto_family_pseudogene_fragment_audit`,
   2,572 members / 605 families): 56.5% of members are pseudogene-biotype; 217 of 605 families (35.9%) are entirely pseudogene;
   147 of 420 size-comparable families (35.0%) put a fragment (< 20% of the largest member's exonic bases) next to a full-length member
   (worst ratio 100-126x, ANKRD20A1 / ANKRD20A3P).

## 4. What this does and does not show

- **Shown:** Soto's families are a copy-number refinement of sequence families that we reproduce (nesting 99.1% sequence-only;
  0.93 with our own copy numbers); their families miss real paralog pairs and include fragment and pseudogene members.
- **Not shown, and not claimable from this data:** that ours is "better". The objectives differ: Soto splits by copy number to
  study human-specific expansions; we define homology groups. The accuracy comparison uses Soto's table as truth, so it cannot
  rank the two. Our own RNA-level families also cover only expressed loci, and our 0.93 still uses a copy-number gate.
- **Safe phrasing:** "Soto's 444 multi-gene families nest inside our sequence families (99.1%); with their copy-number step
  (our own copy numbers: 0.93 ARI, 411 / 491 exact; theirs: 0.97, 479) we reproduce them. Where they differ, theirs are narrower
  than sequence homology: 83 real paralog pairs are missing, about 700 near-identical exon links are split by copy number, and a
  third of their families are entirely pseudogene."

## 5. Regenerate (two minutes, no heavy data)

Needs a clone and `famcn_ours_allwssd.tsv` (0.2 MB; `famcn_ours_all.tsv`, 0.2 MB, for the 10-sample row).

```
python3 bench/soto/soto_replication.py genesets --out-eligible elig.tsv --out-full full.tsv
python3 bench/soto/soto_replication.py nesting --shared bench/soto/shared_exons_5154_exon_mapback.tsv \
    --geneset elig.tsv --full-geneset full.tsv --famcn-ours famcn_ours_allwssd.tsv
python3 bench/soto/soto_replication.py ladder  --shared bench/soto/shared_exons_5154_exon_mapback.tsv \
    --geneset elig.tsv --full-geneset full.tsv --famcn-ours famcn_ours_allwssd.tsv \
    --famcn-ours10 famcn_ours_all.tsv --split bench/soto/soto_split_2026-09-29.tsv --drop-family ID_356
```
