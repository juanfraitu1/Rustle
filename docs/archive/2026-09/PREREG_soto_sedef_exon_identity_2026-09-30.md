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

Run on 2026-09-30 after this file was committed (`cdaa8338`, sha1 `80ca836f39ea7fa5f307a76f545560956ee7a7b9`);
`bench/soto_m2/soto_m2_sedef_exonid.py` (miniforge python, edlib), 14 s, light lock. Output: `docs/archive/2026-09/SOTO_SEDEF_EXON_IDENTITY_2026-09-30.md`;
per-link identities: `docs/SOTO_SEDEF_EXON_IDENTITY_LINKS_2026-09-30.tsv`. Implementation choice not fixed above: B and C collapse duplicate
edges, as KEY=unionedges did.

Correctness gate passed: the 4,384 re-derived SEDEF links equal the frozen native-v1 edge file (3,628 rows used, 57,456 exon projections).

**VERDICT (section 4): DOES NOT CONVERGE.** S1C copy numbers, all genes: A 479, B 397, C 406 exact (midpoint 438); C recovers 0 and
breaks 73 of A's exact families (held-out 263 -> 214). Every other setting goes the same way (ours, all genes: 411 / 345 / 353, 58 broken).

Why: the exon pairs are already >= 98% identical. Exon-pair identity of all 4,384 SEDEF links, quartiles 0.994 / 1.000 / 1.000; 4,066
(92.7%) are >= 0.98. Of the 1,012 SEDEF-only links (not in the map-back), 856 are >= 0.98 (quartiles 0.986 / 0.996 / 1.000); the filter
removes only 156.

Descriptive, after the verdict: SEDEF-only links mostly involve a non-coding gene (813 of the 856 high-identity ones; 43 join two coding
genes) against 2,492 of 3,372 map-back links; gene strand does not separate them (opposite-strand pairs 46% vs 47%). The damage passes
mostly through those genes: with pseudogenes and lncRNAs removed, B and C each break 12 families.

**Reading.** SEDEF links are not low-identity noise: they are >= 98%-identical exon pairs that Soto's exon map-back does not report
(the map-back queries only exons fully inside a merged SD98 region, keeps up to 50 secondary hits, and needs >= 99% cover of the target
exon). Soto's families are therefore defined by their map-back procedure, not by exon identity alone: at the same 98% exon identity,
SEDEF's links join genes Soto keeps apart. The fair comparison stays the map-back graph.
