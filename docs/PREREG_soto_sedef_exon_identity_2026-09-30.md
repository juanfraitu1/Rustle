# Pre-registration: SEDEF links restricted to exon pairs that are themselves >= 98% identical (2026-09-30, KEY=sedefexon)

Written before any exon-level identity was computed. Human CHM13 v1.0, Soto's 2,334-gene table. Nothing in `src/` or `tools/`
changes. Script: `bench/soto_m2/soto_m2_sedef_exonid.py`, written after this file is committed.

## 1. Question

KEY=unionedges found that adding SEDEF-projected exon links to the exon map-back graph breaks 82 of Soto's exact families (S1C copy
numbers, all genes: 479 -> 397). SEDEF rows are kept at >= 98% identity, but that identity is the whole duplication alignment's; the
exons carried through it are never compared. **If a SEDEF link is kept only when the exon pair itself is >= 98% identical, does the
union come back toward the map-back result?** If yes, the damage came from exon pairs below Soto's own 98% scope, and SEDEF and the
map-back agree at the exon level. If no, SEDEF links at exon-level 98% still join genes Soto keeps apart.

## 2. Inputs (frozen)

- SEDEF: the native CHM13 v1.0 table `winloci_data/soto_replication/final_v1_clean.bed` (rows with identity >= 0.98, as in
  `soto_replication.py edges --native-v1`). Only this table is used (the v2.0 table is not: the native arm is the one without a
  liftover, and `soto_replication.py` records that its union with the v2.0 edges scores identically in the replication chain).
- Genome: `winloci_data/soto_replication/t2t-chm13-v1.0.fa.gz`. Genes and exons: `cat_v4.bed`, gene set `full.tsv` (2,334 genes).
- Everything else as KEY=unionedges: map-back edges, S1C, our copy numbers, the frozen dev / held-out split, `classify` unchanged.

## 3. Arms and measures

- **A:** map-back edges alone.
- **B:** A plus every native-v1 SEDEF link, re-derived by this script with `soto_replication.py`'s own CIGAR walk. **Correctness gate:**
  B's SEDEF link set must equal the frozen `bench/soto/shared_exons_2334_finalv1_native.tsv` exactly, else stop.
- **C:** A plus the native-v1 SEDEF links that have at least one projection with **exon-pair identity >= 0.98**: for an exon E (the part
  inside the SEDEF row's span) projected through the row's CIGAR onto interval P on the other side, identity = 1 - (global edit
  distance between E and P) / max(|E|, |P|), with P reverse-complemented when the row's second side is on '-'. Edit distance: edlib,
  global mode.
- Per arm, with S1C and our copy numbers and the four biotype filters: exact families (all, dev, held-out), ARI, recovered / broken
  vs A. Also: number of SEDEF links and of SEDEF-only links (not in A) in B and C, and the identity distribution of B's links.

## 4. Decision rule (primary: S1C copy numbers, all genes)

- **CONVERGES:** C's exact count is above the midpoint of A and B, (exact_A + exact_B) / 2.
- **DOES NOT CONVERGE:** otherwise.
Recovered and broken families are reported beside it, with no verdict attached.

## 5. Seen before (disclosed)

- KEY=unionedges (both SEDEF tables, duplicates collapsed): 397 exact, 0 recovered, 82 broken; 338 extra genes, 283 with famCN at
  median |difference| 0.49 from the nearest member, 325 in another Soto family. The pair rule is sensitive to repeated edges.
- No exon-level identity of any SEDEF link has been computed or looked at.

## 6. Result

(Filled in after the run, below this line, without editing anything above.)
