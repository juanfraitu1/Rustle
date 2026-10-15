# Soto's families as our families cut by copy number: T1 and T2 (2026-10-14)

**Pre-registered before any T1 or T2 statistic is computed.**
- Spec: `docs/superpowers/specs/2026-10-14-soto-subset-cn-design.md`.
- Plan: `docs/superpowers/plans/2026-10-14-soto-subset-cn.md`.
- Code: `bench/soto_subset/` at commit `25006d97`. Every number below comes from it, run with the defaults (1,000 null
  draws, 1,000 bootstrap resamples, seed 7):

      bash tools/rlock.sh light python3 bench/soto_subset/run_subset.py guided --out /mnt/linuxdisk/tmp/o1_soto_subset/guided
      bash tools/rlock.sh light python3 bench/soto_subset/run_subset.py denovo --sample human_A119b --out /mnt/linuxdisk/tmp/o1_soto_subset/human_A119b
      bash tools/rlock.sh light python3 bench/soto_subset/run_subset.py denovo --sample human_testis --out /mnt/linuxdisk/tmp/o1_soto_subset/human_testis

## Question

- **T1, Soto ⊂ ours:** are Soto's families (Soto 2025, Table S1C) a refinement of our O1 multi-copy gene families?
- **T2, Soto ≈ ours ∩ CN:** does cutting ours by copy number, with Soto's own pair rule, rebuild them?

Human only. Nothing is tuned on Soto's table. The pair rule is Soto's (|Δ| < 2), and the nulls, draw counts and seed
are fixed above.

## Inputs (frozen)

| What | Path |
|---|---|
| Soto's families and the S1C copy number | `bench/soto/soto_famCN_S1C.tsv` |
| Assembly copy number (HG002 diploid, `fam_HG002`) | `/mnt/linuxdisk/tmp/rustle_figures_dev/soto_parcn_asm/results_s1e.tsv` |
| Read-depth copy number (`famCN_sotoiv`, 268 SGDP genomes) | `/mnt/linuxdisk/home/juanfraitu/winloci_data/soto_replication/famcn_ours_allwssd.tsv` |
| Soto's genes on CHM13 v2.0 (de novo mapping, coordinate check) | `/mnt/linuxdisk/tmp/o1_close/lift/soto_truth.v2.tsv` |
| Guided families (F_O, the 10-09 run) | `/mnt/linuxdisk/tmp/o1_layers/v2/dupmasker/ours.hierarchy.members.tsv`, `ours.loci.tsv`; loci `/mnt/linuxdisk/tmp/ourhom/loci.gff3` |
| De novo catalogs (exact-mcl; the same partitions as the shipped `<sample>.fam.clusters.tsv`) | `/mnt/linuxdisk/tmp/o1_close/hier_shipped/<sample>/ours.hierarchy.members.tsv`, `ours.loci.tsv`; loci `/mnt/linuxdisk/tmp/rustle_figures/runs/<sample>/<sample>.fam.loci.gff3` |
| Soto's recipe | `bench/soto/soto_replication.py cluster --pair-mad` on `bench/soto/shared_exons_5154_exon_mapback.tsv` |

## Known before running (recorded now)

- **T1 on the guided universe is not blind, and strict containment already failed.** On F_O (the 10-09 run, at
  member level, over the 400 evaluable families), 313 (78.2%) of Soto's families are fully contained in one of our
  families and 346 (86.5%) in one superfamily; 215 equal one of our expansions. That prereg's H1 (≥ 90% contained)
  was NOT met, and it concluded that Soto's families are not a strict subset of ours
  (`docs/PREREG_family_hierarchy_flag_2026-10-09.md`). T1 asks a different question: containment beyond random
  partitions with our block sizes, and which partition is the coarser. A T1 pass therefore does not make Soto a strict
  subset (ruling R7). The de novo catalogs are the replication.
- **T1a is expected to pass easily.** Plan 1's global null was passed by design (register 1277), and random
  partitions with our block sizes rarely contain a Soto family. T1a only rules out containment by chance; the
  informative T1 quantities are the containment rate, the direction D and the violation table.
- **Older catalogs nest too.** On the older RNA lattice catalogs, 80-90% of Soto's families nest in ours at L1/L2.
  On held-out chr5/7/21 it was 61/76: 53/59 families of ≤ 3 genes and 8/17 larger ones
  (`bench/SOTO_AS_A_REFINEMENT.md`).
- **Copy number already raised agreement on our edges.** In an exploratory, unregistered run (2026-10-07), our F_O
  graph's edges gated by Soto's pair rule gave ARI 0.7142 / 0.7849 / 0.8263 (no copy number / read depth / S1C) at
  identity ≥ 0.70. This study cuts our families rather than our edges, and adds assembly copy number.
- **Soto's recipe rungs** (registers 1162 and 1169; reproduced 2026-10-14), over the 491 families: no copy number
  0.7307 / 345 exact, read depth 0.9277 / 411, S1C 0.9698 / 479.
- **Assembly vs S1C:** assembly copy number agrees with S1C at Spearman 0.617 (human; register 1172).
- **No held-out half is needed.** S1C's family-hash halves are both spent for copy-number split rules (2026-10-07).
  Nothing here selects a rule, so no half is needed; the de novo catalogs replicate.
- **FAM90A.** ARI over Soto's families swings ±0.035 with FAM90A (register 1170). Hence exact families are a
  co-requirement of T2's bar, and ARI is also reported without FAM90A.
- **Input counts** (`inputs.json`; inputs, not outcomes):

| | guided | A119b | testis |
|---|---|---|---|
| members | 2,334 genes | 1,396 loci | 332 loci |
| Soto families with ≥ 2 members | 491 | 293 | 105 |
| single-label / multi-label members | 2,185 / 149 | 1,190 / 206 | 300 / 32 |
| members on our families | 1,796 | 1,115 | 238 |
| members with assembly / read-depth / S1C copy number | 1,534 / 2,334 / 1,793 | 1,175 / 1,396 / 1,196 | 288 / 332 / 292 |
| common set (a value in all three) | 1,534 | 1,175 | 288 |
| Soto families with ≥ 2 members in the common set | 357 | 257 | 85 |
| FAM90A members in the common set | 2 | 0 | 0 |

## Definitions

- **Guided universe.** Members are Soto's 2,334 genes. A gene's locus is its gene line in the F_O input GFF (by CAT
  v4 gene id), folded through `ours.loci.tsv`.
  - Genes on loci of one F_O family share a block.
  - Genes on one unclustered locus share a block.
  - Every other gene is our singleton.
- **De novo universe.** Members are the catalog loci reached by a gene of Soto's 491 families: an exon shares ≥ 1 bp,
  on the same strand, with a locus exon, and a folded locus counts as its representative.
  - Blocks are the catalog's families. A reached locus in no family is our singleton.
  - A locus carries the families of every gene that reaches it.
- **Soto's labels.** `Unassigned_*` genes are Soto singletons. A member with two or more labels is *multi-label*.
- **Copy number of a member** = the median over its genes that have a value.
- **Common set** = the members with a value in all three sources.

## T1 — Soto ⊂ ours

- **Containment.** A Soto family with ≥ 2 members is contained when all its members lie in one block of ours. The
  statistic is the number contained.
- **Size null.** Random partitions of the members with our exact block-size multiset, singletons included: a uniform
  permutation of the block labels, 1,000 draws, seed 7. p = (1 + draws ≥ observed) / 1,001.
- **T1a bar:** p < 0.001.
- **Direction.** Over the single-label members, D = H(Soto | ours) − H(ours | Soto), in bits.
  - The 95% percentile interval comes from 1,000 bootstrap resamples of *units* (ruling R1), seed 7.
  - A unit is a connected component of the members linked by a shared Soto label or a shared block of ours.
  - A unit drawn twice enters as two disjoint copies.
- **T1b bar:** the interval lies above 0.
- **Reported, no bar:**
  - the rate over all 491 families, and over the families with ≥ 2 members (R6);
  - plan 1's rate: families with ≥ 2 members on clustered members, contained within that restriction;
  - containment for families of ≤ 3 and > 3 members;
  - every violation, classified by the first rule that makes it contained: **coverage** (our singletons removed),
    then **cover** (multi-label members also removed), else **conflict**.

## T2 — Soto ≈ ours ∩ CN

- **The intersection.** Inside each block of ours, members are joined when |CN_a − CN_b| < 2 (Soto's per-pair
  MAD < 1). The components are the intersection: runs of the block's sorted values, split at gaps ≥ 2.
- **Scores,** on the same members:
  - **exact:** distinct Soto member sets, restricted to the single-label members with ≥ 2 of them, that equal a
    component (R3);
  - **ARI** over the single-label members, also reported without FAM90A (`ID_356`);
  - **Hungarian Jaccard:** the mean one-to-one best-match Jaccard over Soto's families with ≥ 2 members, all members
    included.
- **Null.** Copy-number values are permuted among the members of each block of ours, 1,000 draws, seed 7, with the
  same p formula. It is computed for each source.
- **T2 bar:** with assembly copy number on the common set, the exact count and ARI both beat their null at p < 0.001.
- **Reported, no bar:**
  - read depth and S1C on the common set, with their nulls;
  - each source on its own coverage;
  - our blocks alone;
  - the ceiling: Soto's partition intersected with ours, on single-label members (R5);
  - Soto's recipe with each source, guided only (R4), on its full gene set and on the common set.

## Rulings (made before any statistic)

- **R1 — bootstrap unit.** T1b resamples units, not Soto's blocks. Two copies of a Soto block inside one block of
  ours would fake a split and inflate H(Soto | ours), biasing toward the bar. D decomposes exactly over units.
- **R2 — a fourth file.** `universes.py` holds the universes beside the spec's three files.
- **R3 — exact over single-label members.** Exact is counted over single-label members, as Soto's own scorer does:
  479/491 excludes the 149 multi-family genes.
- **R4 — recipe rows on guided only.** Soto's recipe is gene-level.
- **R5 — the ceiling on single-label members.** The meet of a cover is not a partition.
- **R6 — two rates.** The rate over all 491 families (families with < 2 members count as not contained) and plan 1's
  rate are both reported.
- **R7 — no strict-subset wording.** The spec's statement "Soto's families are a refinement of ours" would claim strict
  containment, which the 10-09 prereg already refuted (78.2% < 90%) and which T1's null-based bars cannot establish.
  Readings 1 and 2 therefore state containment beyond chance with the contained count and the listed exceptions.
  Added while executing the plan, before any statistic.

## Readings (fixed now)

"T1 holds" means T1a and T1b both pass. A de novo "fail" means a bar fails on either catalog.

| Guided | De novo | Reading |
|---|---|---|
| T1 and T2 hold | hold on both catalogs | O1's statement: "Soto's families lie inside our multi-copy gene families far beyond chance, with ours the coarser partition (c of f contained; the exceptions are listed), and cutting ours by copy number rebuilds them (exact x, ARI y; ceiling z)" (R7). Replaces the failed hierarchy claim (register 1277) as the Soto sentence. |
| T1 and T2 hold | fail | The same statement, labelled as shown on the guided universe only. |
| T1 holds, T2 fails | — | Soto ⊂ ours is claimed. The copy-number part is claimed only if read depth passes T2's bar on the same set, labelled human-only and population-based. |
| T1 fails | — | No subset claim; T2 is reported, not claimed. The size and violation tables say where Soto and we disagree. |

## Wiring checks (the run stops, and nothing is read, if one fails)

- The gene-id join agrees with coordinates: 0 disagreements.
- Soto's recipe rungs reproduce: 0.7307 / 345, 0.9698 / 479 and 0.9277 / 411.
- Plan 1's coverage reproduces: guided, 1,796 genes on our families; A119b, 1,211 / 155 / 1,092 Soto records on
  clustered loci / only unclustered / none; testis, 286 / 100 / 2,072.

## Results

Added after the runs. Nothing above this heading changes once a statistic has been seen.

### Guided (2026-10-14)

Run: commit `6d5be546` (code `25006d97`); wall 19.72 s, peak 164,044 KB. Wiring checks: all passed (gene-id join 1,710 / 0 / 121;
recipe rungs 0.7307 / 345, 0.9698 / 479, 0.9277 / 411 reproduced; plan 1's coverage reproduced in `inputs.json`).

**T1**

| Quantity | Value |
|---|---|
| Soto families with ≥ 2 members (`families_ge2`) | 491 |
| contained (`contained`); rate over 491 (`rate_all`); over ≥ 2 (`rate_ge2`) | 258; 0.5255; 0.5255 |
| size null mean / max (`null.mean`, `null.max`) | 0.709 / 4 |
| p (`p`) → **T1a** | 0.000999 (1/1001) → **PASS** |
| H(Soto \| ours) / H(ours \| Soto) | 0.2860 / 0.6855 |
| D [95% interval] (`D`, `D_ci95`) → **T1b** | -0.3995 [-0.6154, -0.1945] → **FAIL** |
| ≤ 3 members: contained / families (`by_size.le3`) | 176 / 271 |
| > 3 members: contained / families (`by_size.gt3`) | 82 / 220 |
| plan 1's rate: contained / evaluable (`plan1`) | 315 / 402 |
| violations: coverage / cover / conflict (`violation_classes`) | 139 / 47 / 47 |

The conflict families (from `violations.tsv`, `class` = conflict), as `label (size, blocks)`: ID_100 (19, 13), ID_104 (9, 4), ID_112 (8, 7), ID_12 (5, 2), ID_136 (6, 4), ID_14 (10, 6), ID_142 (15, 2), ID_147 (4, 2), ID_155 (6, 4), ID_156 (6, 3), ID_160 (6, 3), ID_175 (13, 2), ID_182 (8, 3), ID_187 (5, 4), ID_188 (5, 3), ID_208 (7, 5), ID_211 (4, 2), ID_212 (11, 8), ID_215 (11, 2), ID_22 (11, 3), ID_226 (6, 3), ID_24 (8, 4), ID_241 (3, 2), ID_244 (9, 4), ID_257 (5, 2), ID_261 (8, 3), ID_263 (4, 3), ID_266 (5, 3), ID_267 (6, 2), ID_269 (3, 2), ID_270 (6, 2), ID_28 (29, 21), ID_280 (11, 7), ID_283 (8, 2), ID_284 (13, 3), ID_306 (4, 2), ID_324 (6, 2), ID_35 (41, 7), ID_39 (7, 3), ID_393 (5, 2), ID_400 (6, 3), ID_443 (5, 2), ID_47 (4, 2), ID_63 (11, 4), ID_65 (15, 5), ID_78 (9, 7), ID_8 (14, 4).

**T2 (common set, 1,534 members)**

| Copy number | exact | ARI | ARI w/o FAM90A | Hungarian J | null exact mean (max) | null ARI mean (max) | p exact | p ARI | bar |
|---|---|---|---|---|---|---|---|---|---|
| assembly (T2's bar) | 186 | 0.8274 | 0.8275 | 0.7694 | 163.195 (171) | 0.6836 (0.7303) | 0.000999 (1/1001) | 0.000999 (1/1001) | **PASS** |
| read depth | 206 | 0.8685 | 0.8686 | 0.7890 | 173.116 (181) | 0.7026 (0.7464) | 0.000999 (1/1001) | 0.000999 (1/1001) | (reported; passes) |
| S1C (circular upper bound) | 232 | 0.9166 | 0.9167 | 0.8305 | 186.865 (195) | 0.7336 (0.7788) | 0.000999 (1/1001) | 0.000999 (1/1001) | (reported; passes) |
| our blocks alone (`ours_alone`) | 183 | 0.7890 | 0.7889 | 0.7411 | — | — | — | — | — |
| ceiling: Soto ∩ ours (`ceiling`) | 258 | 0.9425 | 0.9425 | 0.8727 | — | — | — | — | — |

**T2, each source on its own coverage (`own`)**

| Copy number | members | exact | ARI | ARI w/o FAM90A | Hungarian J |
|---|---|---|---|---|---|
| assembly | 1534 | 186 | 0.8274 | 0.8275 | 0.7694 |
| read depth | 2334 | 203 | 0.7835 | 0.7556 | 0.6706 |
| S1C | 1793 | 236 | 0.9340 | 0.9187 | 0.8113 |

**Soto's recipe with each source (`recipe.json`)**

| Copy number | full gene set ARI / exact of 491 | common set exact | common set ARI | common set ARI w/o FAM90A | common set Hungarian J |
|---|---|---|---|---|---|
| none | 0.7307 / 345 | 266 | 0.7465 | 0.7469 | 0.8312 |
| S1C | 0.9698 / 479 | 348 | 0.9894 | 0.9894 | 0.9902 |
| read depth | 0.9277 / 411 | 299 | 0.9248 | 0.9250 | 0.9293 |
| assembly | 0.7410 / 337 | 258 | 0.8656 | 0.8657 | 0.8694 |

Guided verdicts: T1a PASS, T1b FAIL, T2 PASS; read depth on the common set passes T2's bar.

**Post hoc (descriptive, not pre-registered; added after the guided run).**
- Why T1b fails: 495 of the 2,185 single-label genes are our singletons (genes F_O places in no family). They split Soto's
  families and raise H(ours | Soto). On the 1,682 single-label genes that sit in our families, H(Soto | ours) = 0.3715,
  H(ours | Soto) = 0.1809 and D = +0.1906 (ours the coarser). That cut conditions on our own catalog, so it is not a test,
  and no interval was computed. Of the 233 families not contained, 139 are coverage violations.
- T2's null sits below our uncut blocks: the shuffle null's mean (163.2 exact, ARI 0.6836) is lower than
  our blocks alone (183 exact, ARI 0.7890). The bar therefore compares the copy-number cut with random cuts,
  not with leaving our families uncut. Against uncut, assembly copy number gains +3 exact and +0.0384 ARI.

### De novo (2026-10-14)

Run: commit `6d329ed9` (code `25006d97`); A119b wall 15.99 s, peak 319,368 KB; testis wall 6.42 s, peak 183,876 KB.
Wiring checks: all passed (gene-id join 1,710 / 0 / 121; plan 1's coverage reproduced in `inputs.json`).

**T1**

| Quantity | A119b | testis |
|---|---|---|
| Soto families with ≥ 2 members | 293 | 105 |
| contained; rate over 491; over ≥ 2 | 80; 0.1629; 0.2730 | 54; 0.1100; 0.5143 |
| size null mean / max | 0.192 / 4 | 0.291 / 3 |
| p → **T1a** | 0.000999 (1/1001) → **PASS** | 0.000999 (1/1001) → **PASS** |
| H(Soto \| ours) / H(ours \| Soto) | 0.2237 / 1.2667 | 0.1034 / 0.3555 |
| D [95% interval] → **T1b** | -1.0430 [-1.2541, -0.8029] → **FAIL** | -0.2520 [-0.3903, -0.1052] → **FAIL** |
| ≤ 3 members: contained / families | 59 / 126 | 45 / 83 |
| > 3 members: contained / families | 21 / 167 | 9 / 22 |
| plan 1's rate: contained / evaluable | 105 / 247 | 69 / 84 |
| violations: coverage / cover / conflict | 108 / 34 / 71 | 44 / 4 / 3 |

Conflict families: A119b ID_100 (31, 13), ID_113 (19, 4), ID_116 (8, 2), ID_127 (4, 2), ID_132 (6, 2), ID_136 (14, 10), ID_14 (12, 8), ID_142 (23, 9), ID_148 (5, 2), ID_149 (12, 6), ID_151 (11, 8), ID_152 (11, 7), ID_153 (10, 7), ID_155 (8, 4), ID_156 (7, 3), ID_160 (5, 4), ID_163 (14, 13), ID_175 (11, 4), ID_176 (6, 3), ID_178 (7, 5), ID_179 (9, 4), ID_181 (10, 7), ID_182 (14, 9), ID_187 (9, 7), ID_188 (8, 6), ID_203 (4, 2), ID_208 (7, 5), ID_211 (9, 4), ID_212 (14, 11), ID_214 (14, 12), ID_215 (13, 7), ID_218 (13, 10), ID_231 (5, 3), ID_233 (15, 4), ID_24 (6, 2), ID_243 (8, 4), ID_244 (7, 3), ID_245 (15, 7), ID_257 (5, 3), ID_261 (7, 4), ID_262 (12, 3), ID_270 (6, 5), ID_279 (6, 2), ID_280 (35, 18), ID_284 (26, 8), ID_31 (14, 3), ID_35 (45, 13), ID_363 (9, 2), ID_368 (8, 3), ID_369 (4, 2), ID_374 (5, 3), ID_391 (10, 2), ID_393 (6, 3), ID_395 (18, 7), ID_396 (4, 3), ID_397 (9, 5), ID_399 (8, 3), ID_400 (8, 7), ID_403 (2, 2), ID_41 (25, 9), ID_455 (4, 2), ID_468 (21, 3), ID_471 (4, 3), ID_49 (5, 3), ID_56 (6, 4), ID_60 (7, 3), ID_62 (15, 9), ID_65 (13, 6), ID_78 (11, 8), ID_8 (8, 4), ID_99 (27, 14); testis ID_154 (7, 2), ID_395 (5, 2), ID_407 (6, 3).

**T2 (common set: A119b 1,175 loci, testis 288 loci)**

| Copy number | catalog | exact | ARI | ARI w/o FAM90A | Hungarian J | null exact mean (max) | null ARI mean (max) | p exact | p ARI | bar |
|---|---|---|---|---|---|---|---|---|---|---|
| assembly (T2's bar) | A119b | 51 | 0.4022 | 0.4022 | 0.5361 | 47.982 (52) | 0.3674 (0.3873) | 0.006993 (7/1001) | 0.000999 (1/1001) | **FAIL** |
| assembly (T2's bar) | testis | 40 | 0.7574 | 0.7574 | 0.7530 | 39.243 (40) | 0.7312 (0.7708) | 0.243756 (244/1001) | 0.044955 (45/1001) | **FAIL** |
| read depth | A119b | 55 | 0.4414 | 0.4414 | 0.5425 | 50.035 (53) | 0.3730 (0.3941) | 0.000999 (1/1001) | 0.000999 (1/1001) | (reported; passes) |
| read depth | testis | 37 | 0.7723 | 0.7723 | 0.7438 | 37.003 (38) | 0.7247 (0.7858) | 1.000000 (1001/1001) | 0.012987 (13/1001) | (reported; does not pass) |
| S1C | A119b | 57 | 0.4992 | 0.4992 | 0.5600 | 50.460 (56) | 0.4111 (0.4298) | 0.000999 (1/1001) | 0.000999 (1/1001) | (reported; passes) |
| S1C | testis | 40 | 0.8305 | 0.8305 | 0.7634 | 38.262 (40) | 0.7633 (0.8286) | 0.005994 (6/1001) | 0.000999 (1/1001) | (reported; does not pass) |
| our blocks alone | A119b | 51 | 0.4621 | 0.4621 | 0.5428 | — | — | — | — | — |
| our blocks alone | testis | 40 | 0.7861 | 0.7861 | 0.7611 | — | — | — | — | — |
| ceiling | A119b | 63 | 0.5132 | 0.5132 | 0.6314 | — | — | — | — | — |
| ceiling | testis | 41 | 0.8409 | 0.8409 | 0.8081 | — | — | — | — | — |

**T2, each source on its own coverage**

| Copy number | catalog | members | exact | ARI | ARI w/o FAM90A | Hungarian J |
|---|---|---|---|---|---|---|
| assembly | A119b | 1175 | 51 | 0.4022 | 0.4022 | 0.5361 |
| assembly | testis | 288 | 40 | 0.7574 | 0.7574 | 0.7530 |
| read depth | A119b | 1396 | 62 | 0.4255 | 0.4253 | 0.5236 |
| read depth | testis | 332 | 39 | 0.7549 | 0.7549 | 0.6827 |
| S1C | A119b | 1196 | 58 | 0.4986 | 0.4984 | 0.5603 |
| S1C | testis | 292 | 41 | 0.8312 | 0.8312 | 0.7631 |

De novo verdicts: A119b T1a PASS, T1b FAIL, T2 FAIL; testis T1a PASS, T1b FAIL, T2 FAIL.

**Post hoc (descriptive, not pre-registered).** As on the guided universe, most violations are coverage violations
(A119b 108 of 213, testis 44 of 51). Cutting our blocks by assembly copy number does not beat leaving them uncut:
A119b 51 exact / ARI 0.4022 against 51 / 0.4621 uncut, and testis 40 / 0.7574 against 40 / 0.7861.

### Verdict (2026-10-14)

Reading 4 of the readings table: T1 fails on the guided universe, so no subset claim is made, and T2 is reported, not
claimed.
- **T1a (guided): PASS.** 258 of 491 Soto families are contained, against a size-null mean
  of 0.709 (max 4); p = 0.000999.
- **T1b (guided): FAIL.** D = -0.3995 [-0.6154, -0.1945], below 0: H(Soto | ours) = 0.2860 is smaller than
  H(ours | Soto) = 0.6855, so ours is not the coarser partition on this universe.
- **T2 (guided, assembly copy number, common set of 1,534 genes): PASS, reported, not claimed.**
  - Exact 186, against a null mean of 163.2 (max 171); p = 0.000999.
  - ARI 0.8274, against a null mean of 0.6836 (max 0.7303); p = 0.000999.
- **Replication: both de novo catalogs pass T1a and fail T1b and T2.**
  - T1b: A119b D = -1.0430 [-1.2541, -0.8029]; testis D = -0.2520 [-0.3903, -0.1052].
  - T2, assembly: A119b p exact 0.006993, p ARI 0.000999; testis p exact 0.243756, p ARI 0.044955.
- **Read depth passes T2's bar on the guided common set** (exact 206, ARI 0.8685). Reading 3 would use
  it, but only when T1 holds.

**What this says.**
- Soto's families sit inside ours far beyond chance (T1a, every universe), but not as a refinement.
- On the guided universe our catalog leaves 495 of Soto's single-label genes in no family. Of the 233 uncontained
  families, 139 are contained once those singletons are removed, 47 once Soto's multi-family genes are removed too,
  and 47 are genuine conflicts.
- On the genes we do group, ours is the coarser partition (post hoc, D = +0.1906).
- The copy-number cut beats random cuts on the guided universe, but gains only
  +3 exact and +0.0384 ARI over leaving our families uncut, and fails on both de novo catalogs.
