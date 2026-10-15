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
