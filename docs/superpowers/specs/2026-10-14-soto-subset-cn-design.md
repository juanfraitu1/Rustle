# Soto's families as our families cut by copy number — design (2026-10-14)

## Intent (agreed with the user)

**The user's words (2026-10-14):** "going back to O1, if I want to show that I can replicate Soto but their families
are different from multi-copy gene families that I have, what can we do to show if Soto could be a subset or if we
could also have a CN-like measurement"; then "I think we can try that"; design sections 1-2 approved.

**The claim to test:** Soto's families are a refinement of our O1 families, obtained by cutting them with copy number.
- **Soto ⊂ ours:** each Soto family lies inside one of our families.
- **Soto ≈ ours ∩ CN:** our families intersected with copy-number classes rebuild Soto's families.

If both hold, the thesis can say: we replicate Soto, and Soto's families are the copy-number-resolved, young subset of
our multi-copy gene families. They are not a competing definition.

**Decisions (user):**
1. **Approach:** containment plus intersection. Not a unit cover (human-only, more complex) and not a lattice-level
   pick (threshold-driven).
2. **Universe:** both. Primary is the guided families on Soto's own genes. The replication, on a held-back substrate,
   is the shipped de novo catalogs.
3. **Copy number:** assembly copy number is primary (the user's standing "assembly everywhere" decision from the
   layered spec). Population read depth is reported beside it; Soto's own S1C is reported as a circular upper bound.

**Corrections from the fact check (2026-10-14, before anything ran).** Each changes the approved section 2.
1. **CN tables are joined by CAT v4 gene id, not by coordinates.** All three tables carry that id (an identifier, not a
   gene name). The assembly table is in CHM13 v2.0 coordinates and Soto's genes in v1.0, so a coordinate join would
   need a lift. The id join was checked against coordinates: every S1E region overlaps its gene's lifted body
   (1,710 rows, 0 disagree; 121 rows belong to genes absent from the lifted truth table).
2. **T2 is scored on the members that have a copy-number value.** The approved rule, "a gene without CN never cuts a
   pair", fails inside a family: such a gene pairs with every other member, so it joins every copy-number class into
   one. Soto's own code also drops a pair where neither gene has a value.
3. **T2's bar adds the exact-family count next to ARI.** ARI swings ±0.035 on whether FAM90A (56 genes) stays whole
   (register 1170, a recorded metric trap).
4. **The guided universe is all 2,334 Soto genes, our singletons included.** Plan 1 counted only genes on clustered
   loci, which left 400 evaluable families. Recall is not conditioned on our catalog here; the plan-1 style rate is
   reported beside it.

## What exists (verified 2026-10-14)

- **Soto's families (S1C, `bench/soto/soto_famCN_S1C.tsv`):**
  - 2,334 genes, 491 families (every one with ≥ 2 genes), and 114 `Unassigned_*` labels of one gene each. Unassigned
    genes are Soto singletons here.
  - It is a cover: 149 of 2,334 genes (6.4%) are in two or more families.
- **Earlier containment study** (`bench/SOTO_AS_A_REFINEMENT.md`; older RNA lattice catalogs, descriptive):
  - 80-90% of Soto families lie inside one of ours at L1/L2, including all 6 NPIP-side families.
  - Held-out chr5/7/21: 61/76 (80.3%). On those 76 there is a size effect: families of ≤ 3 genes nest 53/59 (89.8%),
    larger ones 8/17 (47.1%).
- **Soto's recipe with each copy number** (all 491 families; ladder of `docs/archive/2026-09/PREREG_soto_famcn_allwssd_2026-09-29.md`,
  registers 1162, 1168-1170):

  | Copy number used | ARI | Exact |
  |---|---|---|
  | none (sequence only) | 0.7307 | 345 |
  | 268-genome read depth, all merged exons | 0.8855 | 375 |
  | 268-genome read depth, Soto's interval (`famCN_sotoiv`) | 0.9277 | 411 |
  | S1C | 0.9698 | 479 |

  Copy number takes agreement from 0.73 to 0.97.
- **Copy-number tables on disk.** All three are keyed by CAT v4 gene id (`CHM13_G…` / `LOFF_G…`).

  | Source | File and column | Soto genes with a value |
  |---|---|---|
  | assembly | `/mnt/linuxdisk/tmp/rustle_figures_dev/soto_parcn_asm/results_s1e.tsv`, `fam_HG002` | 1,534 |
  | read depth | `winloci_data/soto_replication/famcn_ours_allwssd.tsv`, `famCN_sotoiv` | 2,334 |
  | Soto's own | `bench/soto/soto_famCN_S1C.tsv`, `Median famCN` | 1,793 |

  - Assembly: one row per S1E region (gene ∩ SD98, CHM13 v2.0). `fam_HG002` is the HG002 diploid sum under
    QuicK-mer2's 30-mer rule, and a region with < 100 usable 30-mers is unresolved. A gene's value is the median over
    its resolved rows.
  - The common set (a value in all three) is the 1,534 assembly genes. 357 of the 491 Soto families keep ≥ 2 genes in
    it, and FAM90A keeps 2 of its 56.
  - Assembly copy number agreed with S1C only moderately before (Spearman 0.617, human; register 1172).
- **Soto's copy-number rule** (`bench/soto/soto_replication.py`, `pair_kept`):
  - A pair is kept when at least one member has a value and the MAD of the values present is < 1.
  - The MAD is the unscaled median absolute deviation. For two genes it is |Δ| / 2, so the pair is kept iff |Δ| < 2.
  - A pair with no value is dropped.
- **Our families:**
  - Guided: F_O, the families over Soto's 2,334 genes (input `/mnt/linuxdisk/tmp/ourhom/loci.gff3`).
    - Files: `/mnt/linuxdisk/tmp/o1_layers/v2/dupmasker/ours.clusters.tsv` (families of ≥ 2 loci) and
      `ours.loci.tsv` (93 annotations folded into overlapping loci). Byte-identical under today's default to the
      registered 10-09 run.
    - Plan 1's key mapping places 1,937 of the 2,458 Soto family records on a family member; the rest are our
      singletons.
  - De novo: the shipped catalogs, `/mnt/linuxdisk/tmp/rustle_figures/runs/<sample>/<sample>.fam.clusters.tsv`.
    - A119b has 2,150 families and testis 342. Both are exact-mcl and the same partitions as
      `/mnt/linuxdisk/tmp/o1_close/hier_shipped/<sample>/`.
    - Soto genes reach their loci by plan 1's exon + strand rule (`/mnt/linuxdisk/tmp/o1_close/lift/soto_truth.v2.tsv`,
      CHM13 v2.0).

## 1. The tests

**Universes.**
- **Guided (primary).**
  - Members are Soto's 2,334 genes.
  - Our blocks: genes on loci of the same F_O family, or on the same locus, share a block. Every other gene is our
    singleton.
- **De novo (replication).**
  - Members are the catalog's loci reached by a gene of Soto's 491 families.
  - Our blocks are the catalog's families restricted to those loci; a reached locus in no family is our singleton.
  - A locus carries the Soto families of every gene that reaches it.
  - Soto genes that reach no locus are reported as catalog coverage. They are not members.

**T1 — Soto ⊂ ours.**
- **Containment.** A Soto family with ≥ 2 members is contained when all its members are in one block of ours.
  - The statistic is the number contained.
  - Null: random partitions of the members with our exact block-size multiset, singletons included. 1,000 draws,
    seed 7, and p = (1 + draws ≥ observed) / 1,001.
  - **T1a bar:** p < 0.001.
- **Direction.** Over the members with exactly one Soto label:
  - compute H(ours | Soto) and H(Soto | ours), conditional entropies in bits;
  - D = H(Soto | ours) − H(ours | Soto);
  - 95% percentile interval from 1,000 bootstrap resamples of Soto's blocks, seed 7.
  - **T1b bar:** the interval lies above 0 (ours is the coarser of the two).
- **Reported, no bar:**
  - the containment rate over all 491 families, and over plan 1's evaluable families (≥ 2 members on clustered loci);
  - containment by Soto family size (≤ 3 vs > 3 members);
  - every violation, listed and classified by the first rule that makes it contained:
    1. **coverage** — contained once our singletons are removed;
    2. **cover** — contained once Soto's multi-family members are also removed;
    3. **conflict** — otherwise.

**T2 — Soto ≈ ours ∩ CN.**
- **Copy-number universe.** These are the members with a value in the source.
  - A de novo locus takes the median of the values of the Soto genes reaching it.
  - The primary comparison uses the common set (a value in all three sources). Each source on its own coverage is
    also reported.
- **The intersection.** Inside each block of ours, two members are joined when |CN_a − CN_b| < 2 (Soto's per-pair
  rule). The connected components are the intersection.
- **Scored against Soto,** on the same members:
  - exact: Soto families with ≥ 2 members that equal a component, counted over distinct member sets;
  - ARI over the members with exactly one Soto label, also reported without FAM90A (ID_356);
  - a cover-aware variant: the mean one-to-one (Hungarian) best-match Jaccard over the Soto families with ≥ 2 members.
- **Null:** copy-number values permuted among the members of each block of ours. This keeps each block's values and
  breaks the member-to-value link. 1,000 draws, seed 7, the same p formula, computed for each source.
- **T2 bar:** with assembly copy number on the common set, both the exact count and ARI beat their null at
  p < 0.001.
- **Reported, no bar:**
  - read depth and S1C, each with its null;
  - our blocks alone;
  - Soto's own recipe with each source (`soto_replication.py cluster --pair-mad`; its full gene set, as registered);
  - the ceiling: the intersection of Soto's partition with ours, scored the same way. It is the best any
    copy-number cut of our blocks can do.

**Replication.** T1 and T2 are judged on the guided universe. The same bars are then run on the de novo A119b and testis
catalogs. The claim **replicates** when every bar holds on both.

**Readings, fixed now.** "T1 holds" means T1a and T1b both pass. A de novo "fail" means a bar fails on either catalog.

| Guided | De novo | Reading |
|---|---|---|
| T1 and T2 hold | hold on both catalogs | O1's statement: "Soto's families are a refinement of our multi-copy gene families; cutting ours by copy number rebuilds them (exact x, ARI y; ceiling z)". Replaces the failed hierarchy claim (register 1277) as the Soto sentence. |
| T1 and T2 hold | fail | The same statement, labelled as shown on the guided universe only. |
| T1 holds, T2 fails | — | Soto ⊂ ours is claimed. The copy-number part is claimed only if read depth passes T2's bar on the same set, labelled human-only and population-based. |
| T1 fails | — | No subset claim; T2 is reported, not claimed. The size and violation tables say where Soto and we disagree. |

## 2. Components (Python, `bench/soto_subset/`, test-first)

- **`subset_lib.py`** (pure functions, no I/O):
  - containment count and the size-preserving null;
  - conditional entropies and the bootstrap interval;
  - the intersection under the pair rule and the copy-number shuffle null;
  - ARI, the exact distinct-set count and the cover-aware Hungarian Jaccard;
  - the violation classifier.
- **`cn_tables.py`:** one table per copy-number source, gene id → value.
  - Assembly: the median `fam_HG002` over the gene's resolved S1E rows. Read depth: `famCN_sotoiv`. S1C: as published.
  - Each Soto gene is accounted for as valued, unresolved or absent, and the counts are reported.
  - The coordinate check of correction 1 runs every time and fails loudly on a disagreement.
- **`run_subset.py`:** runs T1 and T2 for one universe. It writes `t1.json` and `t2.json` under
  `/mnt/linuxdisk/tmp/o1_soto_subset/<universe>/`.
- **Reuse:** plan 1's mappings (`bench/o1_hier/soto_truth.py`, `hier_match.py`: key mapping with the fold table for
  F_O, exon + strand for de novo), and `bench/soto/soto_replication.py` for the recipe rows.
- **Tests:**
  - toy partitions with hand-computed containment, entropies, intersection components and violation classes;
  - seeded nulls that are deterministic;
  - a gene-id join fixture, including an unresolved gene and a coordinate disagreement.

## 3. Order

1. `subset_lib` and `cn_tables` with their tests.
2. The pre-registration `docs/PREREG_soto_subset_cn_2026-10-14.md`, copying section 1's bars and readings, committed
   before any real run.
3. The guided runs, then the de novo replication.
4. Results into the prereg, then register rows. `docs/THESIS_OBJECTIVES.md`'s O1 line changes only if the readings
   say so.

## Out of scope

- Gorilla: it needs Soto projected by Liftoff (a later, descriptive extension).
- RNA expressed copy number (later).
- Any change to O1's family definition or defaults.
- The unit cover, and phases 3-5 of the layered spec. This study can inform phase 3's S2 but does not run it.

## Results (2026-10-14)

Reading 4 of section 1's readings table: Soto as a refinement is not claimed. T1b failed on the guided universe and both
de novo catalogs, while T1a passed, with 258 of 491 Soto families inside one of our families against a null mean of 0.7. T2 passed on the guided universe but is not claimed, and it failed on both de novo catalogs. The wording of readings 1 and 2 was
qualified before any statistic (prereg ruling R7). Full results: `docs/PREREG_soto_subset_cn_2026-10-14.md` (Results,
Verdict); register rows 1281 and 1282.
