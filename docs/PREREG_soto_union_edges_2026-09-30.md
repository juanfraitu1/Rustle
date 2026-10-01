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
