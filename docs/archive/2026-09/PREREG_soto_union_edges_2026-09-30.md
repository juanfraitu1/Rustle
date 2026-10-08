# Pre-registration: does adding SEDEF links to the exon map-back graph reproduce more of Soto's families? (2026-09-30, KEY=unionedges)

Written before the union graph was built or scored. Human CHM13, Soto's 2,334-gene table (S1C). Nothing in `src/` or `tools/`
changes. Script: `bench/soto_m2/soto_m2_union.py`, written after this file is committed.

## 1. Question

With Soto's own copy numbers, our reconstruction (exon map-back edges, Soto's per-pair MAD < 1 rule) rebuilds 479 of Soto's 491
families exactly. Three of the 12 misses (ID_62, ID_192, ID_347) are attributed to `our_edges`: a SEDEF-projected exon link joins
the pieces where no map-back link does. **Does the union graph (map-back plus SEDEF-projected links) recover those families without
merging families Soto keeps apart?**

## 2. Inputs (frozen)

Exactly the inputs of `soto_m2_families.py` (meeting page run): map-back edges `bench/soto/shared_exons_5154_exon_mapback.tsv`;
SEDEF-projected edges `bench/soto/shared_exons_2334_finalv1_native.tsv` + `shared_exons_2334_finalhuman.tsv`; S1C; our copy numbers
`famcn_ours_allwssd.tsv`; gene sets `elig.tsv` / `full.tsv`; the frozen split `bench/soto/soto_split_2026-09-29.tsv` (a family's half =
the half of most of its clean members, as in `soto_m2_loosen.py`).

## 3. Arms and measures

- **Baseline:** map-back edges. **Union:** map-back edges plus every SEDEF-projected edge (duplicates collapsed).
- Each arm classified by `soto_m2_families.classify` (unchanged), with S1C copy numbers and with ours, for the four biotype filters.
- Per arm: exact families (all, dev, held-out), ARI, and per family the change baseline -> union: **recovered** (not exact -> exact)
  and **broken** (exact -> not exact).

## 4. Decision rule (primary setting: S1C copy numbers, all genes)

- **UNION HELPS:** exact rises, no baseline-exact family is broken, and held-out exact does not drop.
- **UNION HURTS:** at least one family is broken and exact does not rise.
- **MIXED:** exact rises but at least one family is broken.
- **NO CHANGE:** nothing recovered, nothing broken.

The other seven settings (our copy numbers; the three filters) are reported beside it, without a verdict.

## 5. Seen before (disclosed)

- With S1C copy numbers, all genes: 479 exact; the 12 misses are 7 `soto_rule`, 1 manual merge (UBTFL), 3 `our_edges` (ID_62, ID_192,
  ID_347), 1 `unassigned` (ID_401). With our copy numbers, 19 misses are `our_edges`.
- The `our_edges` attribution itself is computed on the union graph (`components_over(edges + sedef)`), so for those families the union
  joins the pieces by construction; what is unknown is whether the copy-number rule then reproduces them, and what else the union merges.
- Nesting on SEDEF-projected edges alone is 77.9% (346/444) and their ARI 0.71; nothing about the union has been computed.

## 6. Result

(Filled in after the run, below this line, without editing anything above.)

Run on 2026-09-30 after this file was committed (`427da056`, sha1 `de32c7946a723c9de791ed4eb8af63d33a379812`);
`bench/soto_m2/soto_m2_union.py`, 3 s, light lock. Full output: `docs/archive/2026-09/SOTO_UNION_EDGES_2026-09-30.md`.

**VERDICT (section 4): UNION HURTS.** With S1C copy numbers and all genes, the union graph (12,231 map-back + 1,019 SEDEF-only
edges) recovers 0 families and breaks 82: exact 479 -> 397, held-out 263 -> 211, ARI 0.9698 -> 0.9304. The broken families mostly
become `soto_smaller` (1 to 20 extra genes each): SEDEF links join genes that Soto's table keeps apart, i.e. links Soto's code did not
use. Every other setting goes the same way (ours, all genes: 411 -> 345, 66 broken, 0 recovered; filtered settings 12-70 broken, 0
recovered).

Found while checking the three target families (descriptive, after the verdict): the pre-registered union collapses duplicate
edges; the two SEDEF files share 4,174 of their 8,576 rows. Concatenating the files without collapsing gives the same exact count
(397, ARI 0.9350) but recovers ID_62 and ID_192, while the collapsed union recovers neither; ID_347 stays `mixed` (+58 genes) in
both. The map-back graph alone has no duplicates and gives 479 under three random edge orders, so the reconstruction is not
affected; the pair rule's sensitivity to repeated edges matters only once SEDEF links are added. Same verdict either way.

**Reading.** The fair comparison is the map-back graph alone (Soto's released code). SEDEF-projected links describe the duplication
blocks, not Soto's families: adding them merges 82 of Soto's exact families with neighbours.
