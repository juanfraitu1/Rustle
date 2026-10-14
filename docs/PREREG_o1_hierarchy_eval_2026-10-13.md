# PREREG — O1 hierarchy-aware evaluation, human (2026-10-13)

Spec: `docs/superpowers/specs/2026-10-13-close-o1-design.md` section A, with its amendments. Plan:
`docs/superpowers/plans/2026-10-13-close-o1-phase-0A.md` (Tasks 5-6). Where this prereg and the spec or plan differ,
this prereg governs.

**Question.** Are Soto 2025's families nodes of our family hierarchy? The flat one-to-one match against our MCL
families counts every Soto family that sits one level down as a miss.

**Catalogs.**
- `<S>` is `human_A119b` (dev) or `human_testis` (held-out). `R = /mnt/linuxdisk/tmp/rustle_figures/runs/<S>`.
- Both are default catalogs, re-clustered from their cached all-vs-all with the MCL two-node fix ("fix 0",
  `docs/PREREG_mcl_two_node_2026-10-13.md`).
- Binary: `/mnt/linuxdisk/tmp/o1_close/mcl_families.after` (md5 db55a19f8f165587563212c276bcfbca). The pre-fix
  `mcl_families.before` (md5 9914dee365a041347275f2fafba25e07) is used only by the Secondary.
- Lifted inputs (md5): `dupmasker_colors.v2.bed` f5e1146faf778a6d4f52e5bf0c172284; `sedef_pairs.v2.bed`
  7bd521edff9bd1e27b8141027da37a01; `soto_truth.v2.tsv` 60b77a995216237bbc33738dd78e212f.
- Before the run, `md5sum` must reproduce all five values above, and `git diff 8f0515ef -- bench/o1_hier/hier_match.py`
  must be empty. A mismatch stops the run.

**Hierarchy.** `L = /mnt/linuxdisk/tmp/o1_close/lift`; output `O = /mnt/linuxdisk/tmp/o1_close/hier/<S>`.
- Command:
  `mcl_families.after --paf $R/<S>.fam.loci.paf --gff $R/<S>.fam.loci.gff3 --min-exonic-bp 1 --min-shared-exon-frac 0.60 --emit-hierarchy --hierarchy-duplicons $L/dupmasker_colors.v2.bed --hierarchy-sd $L/sedef_pairs.v2.bed --out $O/ours`.
  The nodes are read from `$O/ours.hierarchy.members.tsv`; the fold map is `$O/ours.loci.tsv`.
- **Inputs.**
  - DupMasker colors: `/mnt/linuxdisk/home/juanfraitu/winloci_data/duplicons/chm13.draft_v1.0_plus38Y_dupmasker_colors.bed`.
  - SEDEF: `/mnt/linuxdisk/home/juanfraitu/winloci_data/soto_replication/final_v1_clean.bed`.
  - Both are lifted from CHM13 v1.0 by the T2T chain `/mnt/linuxdisk/home/juanfraitu/winloci_data/soto_replication/v1.0_to_v1.1_rdna_merged.chain`
    (`bench/o1_hier/v1_to_v2.py`).
- **Lift losses and checks.**
  - DupMasker: 6,661 rows dropped (6,306 GRCh38 chrY, 355 acrocentric rDNA).
  - SEDEF: 925 rows dropped. 156 rows (247 sides) span a chain gap and are kept with shifted ends; these are not
    sequence-checked, and most lie within 1 kb of a chromosome end.
  - 1,000 / 1,000 sampled same-length intervals per file are sequence-identical.
- No `--hierarchy-external`: the truth never enters the build.
- Its families must equal the re-clustered catalog:
  `cmp $O/ours.clusters.tsv /mnt/linuxdisk/tmp/o1_close/recluster/after/<S>.fam.clusters.tsv`. A difference stops the
  run.

**Truth.** `$L/soto_truth.v2.tsv` (`bench/o1_hier/soto_truth.py`):
- S1C families proper: 491 families, 2,458 gene rows; `Unassigned_*` excluded.
- Genes are joined to CAT v4 by Gene ID.
- Exons = the union over the gene's transcripts, chain-lifted to v2.0 (0 exons unliftable).

**Mapping.**
- A Soto gene joins a member when one of its exons shares >= 1 bp, on the same strand, with an exon of the member's
  locus or of a locus folded into it.
- Exons come from the catalog's `$R/<S>.fam.loci.gff3`; the fold map is `$O/ours.loci.tsv`.
- The scorer's fallback for a member locus missing from the GFF (its span, either strand) is counted in
  `coverage.member_loci_missing_from_gff`, which is reported. The count covers member and folded-locus keys; absent
  from the JSON means 0. A GFF strand `.` matches either strand.
- A family is EVALUABLE with >= 2 members. The others are reported as coverage.

**Scorer.** `bench/o1_hier/hier_match.py` as of commit 8f0515ef (scipy 1.17.1, numpy 2.4.2):
`hier_match.py $O/ours.hierarchy.members.tsv $L/soto_truth.v2.tsv $O/hier_match.json --map exon --loci-gff $R/<S>.fam.loci.gff3 --fold $O/ours.loci.tsv --null 1000 --seed 7`.
- **Matching.** One-to-one Hungarian (scipy) on -Jaccard between the evaluable truth families and every node: family,
  expansion (including the `<family>.unique` buckets), superfamily and SD complex.
- **Tie policy.**
  - Rows in label order; columns in (level rank family < expansion < superfamily < sd_complex, node id) order.
  - Then canonical placement: the truth families matched to nodes with one identical member set take that set's
    lowest-rank nodes, best Jaccard first, then label (`Scorer._prefer_lower_rank`).
  - Identical sets give identical Jaccards, so the optimum and every number except the reported node and level are
    unchanged.
- **Nulls**, both keeping node sizes and nesting:
  - **global**: one random permutation of all members per round;
  - **within_family**: members permuted inside their own family. Families, superfamilies and SD complexes are
    unchanged; only expansions move.

**Secondary (registered; same statistics, no bars).**
- **Why.** `docs/PREREG_mcl_two_node_2026-10-13.md` (F2) showed that the residue bug split no pair in the shipped
  catalogs: all 6,969 tie-regime pairs in the four default catalogs were already joined (5,127 of them in the two
  scored here).
  - Fix 0 therefore only adds pairs below w* = 0.456885 (A119b +68 families, testis +9).
  - That is a definition change whose adoption is the user's choice (Amendment 1 there).
- **What.** The shipped pre-fix catalogs are also scored. For every two-node component they equal exact-MCL semantics
  (F2: no tie-regime pair was split).
  - They use the same hierarchy inputs, truth, scorer and seed, with `mcl_families.before`.
  - Output: `/mnt/linuxdisk/tmp/o1_close/hier_shipped/<S>`, built with `--out /mnt/linuxdisk/tmp/o1_close/hier_shipped/<S>/ours`.
  - Their families must equal the shipped catalog:
    `cmp /mnt/linuxdisk/tmp/o1_close/hier_shipped/<S>/ours.clusters.tsv $R/<S>.fam.clusters.tsv`. A difference stops
    the Secondary.
- **How reported.** Next to the primary, with no bar. The pair-rule choice must not be made from these scores.

**Bars.** All are computed from the unrounded values in the PRIMARY run's `$O/hier_match.json`.
- `exact` = `hierarchy.exact`, the Hungarian matches with Jaccard 1.
- `rate` = `hierarchy.exact_rate` = `hierarchy.exact` / `hierarchy.n_truth` (the evaluable families).

- **E1:** on each catalog, `null.global.p_exact` < 0.001 AND `null.global.p_mean_jaccard` < 0.001.
- **E1b:** on each catalog, `null.within_family.p_exact` < 0.001 AND `null.within_family.p_mean_jaccard` < 0.001.
- **E2:** |`hierarchy.exact_rate`(testis) - `hierarchy.exact_rate`(A119b)| <= 0.10.
- **p-values:** p = (1 + #{null rounds >= observed}) / (1 + 1000).
  - The smallest is 1/1001 = 0.000999, which passes p < 0.001 only when no null round reaches the observed value.
  - The advisor table prints `p_exact` to 6 dp; `p_mean_jaccard` is quoted from the JSON.

**Reading, fixed now.** The first row that matches applies.
1. **E1 fails on either catalog:** no claim.
2. **E2 fails:** sample-specific; only the per-catalog numbers are reported.
3. **E1b fails (E1 and E2 pass):** the claim stops at the family level and above.
   - Soto's families are matched by our families, superfamilies and SD complexes.
   - The headline rate is (`hierarchy.exact` - `hierarchy.exact_by_level.expansion`) / `hierarchy.n_truth`. Quote it
     with `hierarchy.exact_by_level` and `hierarchy.by_level` excluding `expansion`, next to the flat score.
   - On at least one catalog, the within-family null is not rejected at p < 0.001 on at least one of its two
     statistics. So no claim is made that our expansions recover Soto beyond random within-family splits of the same
     sizes.
4. **E1, E1b and E2 pass:** "Soto's families are nodes of our hierarchy". The O1 headline for Soto becomes, per
   catalog, the exact rate and mean Jaccard, next to the flat score and the coverage.

In every case the headline rate is quoted with `hierarchy.n_truth` (evaluable) out of 491, and with
`hierarchy.exact_of_all`.

**Reported regardless:**
- **Coverage:**
  - records on members / on unclustered loci / without a locus;
  - labels evaluable / one member / none;
  - members in several truth families;
  - `member_loci_missing_from_gff`.
- **Summaries:** the hierarchy and flat summaries: exact, rate, mean Jaccard, macro and micro recall and precision,
  gene recall.
- **Levels:** `exact_any` per level, and where matches land (`hierarchy.by_level`, `hierarchy.exact_by_level`).
- **Unconditioned numbers** over ALL 491 truth families: `hierarchy.exact_of_all`, `macro_gene_recall_all`,
  `micro_gene_recall_all`. A family without a match counts 0, and every truth gene counts, clustered or not.
- **`exact_any.expansion_not_family`:** Soto families equal to an expansion that is not a whole family. These are the
  only matches strictly below the family level. `exact_any.expansion_unique` (equal to a `<family>.unique` bucket) is
  reported separately.
- **Tables:** the per-truth table, and the advisor table (`bench/o1_hier/advisor_table.py`).

**Disclosures.**
- **What was fixed after the F_O probe.** E1b, the exon + strand mapping and the >= 2-member evaluable rule were fixed
  after the plan-time scorer probe on the 10-09 F_O universe (v1.0, guided; 200 null rounds). There, with exon
  mapping:
  - the global null averaged 0.045 exact matches;
  - the within-family null averaged 224 of the observed 231.

  The tie-free exact-expansion count (`exact_any.expansion`) was:
  - 215 with key mapping (the 10-09 Rust placement);
  - 206 with exon + strand mapping;
  - 169 with span-overlap mapping. That mode existed only in the probe; the committed scorer has no span mode.
- **No prior scoring.** No de novo catalog has been scored by this scorer before this prereg.
- **Format-only dry run.** While the plan was written (10-13, before this prereg), the testis hierarchy was built once,
  with the pre-fix binary and scratch copies of the lifted inputs, to check file formats (`/mnt/linuxdisk/tmp/o1_close/probe/`):
  member keys against the loci GFF, and node counts per level. No truth file entered the build, and nothing was
  scored.
- **Precision is a lower bound.** A non-Soto member inside a matched node counts against precision, and Soto's truth
  covers only its SD98 gene set.
- **What the review rounds added.** The unconditioned metrics, the expansion-level split, canonical placement and the
  6-dp printing were added in the scorer's review rounds. They were checked on the 10-09 F_O universe only:
  - Every earlier F_O number was unchanged, except where exact matches land. With key mapping, `exact_by_level`
    (expansion / family / superfamily) is 13 / 228 / 4 instead of 14 / 227 / 4, because an exact match now takes its
    identical family node.
  - `expansion_not_family` = 9 (exon mapping) of the 206 expansion matches.

## Results

Run 2026-10-13 (Task 6 of `docs/superpowers/plans/2026-10-13-close-o1-phase-0A.md`) on the committed prereg (HEAD `fed73e76`); nothing above this heading was edited. Scorer `bench/o1_hier/hier_match.py` at 8f0515ef, python 3.14.4, scipy 1.17.1, numpy 2.4.2, called as `--map exon --null 1000 --seed 7` (the default `--map` is `exon`). The bars below are computed from the UNROUNDED values in the two PRIMARY JSONs (`/mnt/linuxdisk/tmp/o1_close/hier/human_{A119b,testis}/hier_match.json`); the computation is saved in `/mnt/linuxdisk/tmp/o1_close/hier/verdicts.txt`.

### Gates (all passed)

- **Pinned files.** `md5sum` reproduced all five values: `mcl_families.after` db55a19f8f165587563212c276bcfbca, `mcl_families.before` 9914dee365a041347275f2fafba25e07, `dupmasker_colors.v2.bed` f5e1146faf778a6d4f52e5bf0c172284, `sedef_pairs.v2.bed` 7bd521edff9bd1e27b8141027da37a01, `soto_truth.v2.tsv` 60b77a995216237bbc33738dd78e212f.
- **Scorer pin.** `git diff 8f0515ef -- bench/o1_hier/hier_match.py` is empty (the working-tree file equals the pinned commit).
- **Families = the catalog.** Primary: `cmp $O/ours.clusters.tsv /mnt/linuxdisk/tmp/o1_close/recluster/after/<S>.fam.clusters.tsv` is identical for human_testis and human_A119b. Secondary: `cmp /mnt/linuxdisk/tmp/o1_close/hier_shipped/<S>/ours.clusters.tsv $R/<S>.fam.clusters.tsv` is identical for both.
- **`coverage.member_loci_missing_from_gff`** is absent from all four JSONs, i.e. 0: every member and folded-locus key was found in the catalog GFF.
- **Builds** (`--emit-hierarchy` log lines; the log counts the expansions without the `<family>.unique` buckets that `n_nodes.expansion` includes):
  - primary testis: 1133 members, 344 superfamilies, 351 families, 309 expansions, SD complexes 157; primary A119b: 10043 members, 1525 superfamilies, 2218 families, 1584 expansions, SD complexes 598;
  - secondary testis: 1115 members, 335 superfamilies, 342 families, 302 expansions, SD complexes 149; secondary A119b: 9907 members, 1457 superfamilies, 2150 families, 1541 expansions, SD complexes 571.
- **Pre-run timing probe (disclosure).** Before scoring, one probe imported the scorer, loaded A119b, and ran 3 null rounds per kind; no scores or null values were printed or recorded (only elapsed times and peak memory).

Post hoc (added after the run; changes no verdict): **lift caveat.** The lift's sequence check (1,000 random same-length intervals per file, all identical) could not see a cluster at the chromosome ends. Of the 74 lifted SEDEF sides that keep their length and whose v1.0 interval starts within 1 kb of its chromosome's start or ends within 1 kb of its end, 27 differ in sequence between v1.0 and v2.0 (72 and 26 as distinct intervals). They are 74 of 301,428 same-length sides, so a random sample of 1,000 holds 0.25 of them on average. Their effect on the SD-complex layer was not measured. Recomputed by `bench/o1_hier/posthoc_checks.py --lift-ends` (output `lift_ends` in `/mnt/linuxdisk/tmp/o1_close/posthoc_checks.json`, which lists the 27 sides).

### Verdicts (PRIMARY)

| Bar | Rule | A119b (dev) | testis (held-out) | Verdict |
|---|---|---|---|---|
| E1 | `null.global`: `p_exact` < 0.001 and `p_mean_jaccard` < 0.001, on each catalog | p_exact 0.000999000999000999; p_mean_jaccard 0.000999000999000999 | p_exact 0.000999000999000999; p_mean_jaccard 0.000999000999000999 | **PASS** |
| E1b | `null.within_family`: `p_exact` < 0.001 and `p_mean_jaccard` < 0.001, on each catalog | p_exact 0.08291708291708291 (not < 0.001); p_mean_jaccard 0.000999000999000999 | p_exact 0.48951048951048953 (not < 0.001); p_mean_jaccard 0.21878121878121878 (not < 0.001) | **FAIL** |
| E2 | \|`hierarchy.exact_rate`(testis) - `hierarchy.exact_rate`(A119b)\| <= 0.10 | exact_rate = 63/247 = 0.2550607287449393 | exact_rate = 41/84 = 0.4880952380952381 | difference 0.2330345093502988 > 0.10: **FAIL** |

- **E1.** 0.000999000999000999 = 1/1001: no null round reached the observed value in any of the four tests. A119b: observed exact 63 against a global-null mean of 0.001 (max 1); observed mean Jaccard 0.5844 against 0.1934 (max 0.2109). testis: exact 41 against 0.026 (max 2); mean Jaccard 0.7836 against 0.2715 (max 0.2967).
- **E1b.** The within-family null keeps every family and permutes members inside it, so only the expansions move. A119b: observed exact 63 against a null mean of 61.048 (max 64), mean Jaccard 0.5844 against 0.5752 (max 0.5819): the within-family null is not rejected on the exact count (p 0.08291708291708291) but is rejected on the mean Jaccard (p 1/1001), so one of the two statistics passes. testis: exact 41 against 40.569 (max 44), mean Jaccard 0.7836 against 0.7793 (max 0.8007): neither statistic is rejected.
- **E2.** The two exact rates differ by 0.2330345093502988 (exact fraction 4835/20748), more than twice the 0.10 bar; the unconditioned `exact_of_all` goes the other way (A119b 0.1283, testis 0.0835), because the evaluable set is 247 families in A119b and 84 in testis.

### Reading applied: row 2

E1 passes on both catalogs, so row 1 does not match. **E2 fails, so row 2 matches and applies** (the first row that matches applies):

> 2. **E2 fails:** sample-specific; only the per-catalog numbers are reported.

E1b also fails, but row 3 is the case "E1b fails (E1 and E2 pass)" and is not reached. No hierarchy claim for Soto is made across the two catalogs; the per-catalog numbers below stand as reported. The headline rate is quoted with `hierarchy.n_truth` (evaluable) out of 491 and with `hierarchy.exact_of_all`:

| | evaluable of 491 | hierarchy exact | exact rate | exact_of_all | mean Jaccard | flat exact | flat exact rate | flat exact_of_all | flat mean Jaccard |
|---|---|---|---|---|---|---|---|---|---|
| A119b (primary) | 247 | 63 | 0.2551 | 0.1283 (63/491) | 0.5844 | 58 | 0.2348 | 0.1181 | 0.5288 |
| testis (primary) | 84 | 41 | 0.4881 | 0.0835 (41/491) | 0.7836 | 37 | 0.4405 | 0.0754 | 0.7130 |

Post hoc (added after the run; changes no verdict): Descriptive, not a bar: in both catalogs `hierarchy.exact_by_level.family` equals the flat exact count (58 and 37), so the hierarchy's extra exact matches over the flat score are exactly its `exact_by_level.expansion` matches (5 and 4; 1 and 2 of these are twin placements of duplicate Soto families, separated in the post hoc paragraph after the Levels table).

Post hoc (added after the run; changes no verdict): **E2 and set size.** E2's failure stands (its bar compares the unadjusted rates), but half to two thirds of the gap tracks mapped set size (reweighting A119b to testis's size mix leaves 0.079, testis to A119b's mix leaves 0.121; per-size rates in `posthoc_checks.json`, strata 2 / 3 / 4 / 5+), the evaluable truth family's member count in each catalog. The exact rate falls with size in both catalogs: 2 members 43/76 = 0.566 (A119b) and 30/52 = 0.577 (testis); 3 members 10/41 = 0.244 and 6/14 = 0.429; 4 members 4/32 = 0.125 and 3/7 = 0.429; 5 or more 6/98 = 0.061 and 2/11 = 0.182. Testis's evaluable families are smaller: 52 of 84 (0.62) have 2 members against 76 of 247 (0.31) in A119b, and the 77 truth families evaluable in both catalogs have median mapped sizes 5 (A119b) and 2 (testis). Reweighted to testis's size mix (strata 2, 3, 4, 5+), A119b's exact rate is 0.409 against testis's 0.488: a gap of 0.079 instead of the observed 0.233 (by exact size: 0.414, gap 0.074). Restricted to the 77 shared families, the 2-member strata are 12/14 (A119b) and 26/46 (testis). Script `bench/o1_hier/posthoc_checks.py`, output `/mnt/linuxdisk/tmp/o1_close/posthoc_checks.json` (`size_strata`).

### Reported regardless

**Coverage** (records = Soto gene rows). In this and the following tables the two primary catalogs come first and the two Secondary (shipped) catalogs after them, the latter without verdicts:

| | A119b (primary) | testis (primary) | A119b (shipped, secondary) | testis (shipped, secondary) |
|---|---|---|---|---|
| members (`n_members`) | 10043 | 1133 | 9907 | 1115 |
| nodes: family | 2218 | 351 | 2150 | 342 |
| nodes: expansion (incl. `<family>.unique`) | 2863 | 409 | 2785 | 398 |
| nodes: superfamily | 1525 | 344 | 1457 | 335 |
| nodes: sd_complex | 598 | 157 | 571 | 149 |
| truth families (`labels`) | 491 | 491 | 491 | 491 |
| evaluable (>= 2 members) | 247 | 84 | 247 | 84 |
| one member | 55 | 34 | 54 | 34 |
| no member | 189 | 373 | 190 | 373 |
| records | 2458 | 2458 | 2458 | 2458 |
| records on members | 1212 | 287 | 1211 | 286 |
| records on unclustered loci | 154 | 99 | 155 | 100 |
| records without a locus | 1092 | 2072 | 1092 | 2072 |
| records on several members | 298 | 10 | 298 | 10 |
| members in several truth families | 180 | 24 | 180 | 24 |
| `member_loci_missing_from_gff` (key absent = 0) | 0 | 0 | 0 | 0 |

**Summaries** (over the evaluable families; precisions are averages over the matched ones):

| | A119b (primary) | testis (primary) | A119b (shipped, secondary) | testis (shipped, secondary) |
|---|---|---|---|---|
| hierarchy.n_truth | 247 | 84 | 247 | 84 |
| hierarchy.n_matched | 247 | 84 | 247 | 84 |
| hierarchy.exact | 63 | 41 | 63 | 41 |
| hierarchy.exact_rate | 0.2551 | 0.4881 | 0.2551 | 0.4881 |
| hierarchy.mean_jaccard | 0.5844 | 0.7836 | 0.5846 | 0.7852 |
| hierarchy.macro_recall | 0.7073 | 0.9281 | 0.7075 | 0.9305 |
| hierarchy.macro_precision | 0.7810 | 0.8350 | 0.7810 | 0.8350 |
| hierarchy.micro_recall | 0.5999 | 0.9132 | 0.6003 | 0.9170 |
| hierarchy.micro_precision | 0.1920 | 0.7543 | 0.1920 | 0.7543 |
| hierarchy.macro_gene_recall | 0.6720 | 0.6044 | 0.6720 | 0.6044 |
| flat.n_truth | 247 | 84 | 247 | 84 |
| flat.n_matched | 233 | 77 | 233 | 77 |
| flat.exact | 58 | 37 | 58 | 37 |
| flat.exact_rate | 0.2348 | 0.4405 | 0.2348 | 0.4405 |
| flat.mean_jaccard | 0.5288 | 0.7130 | 0.5289 | 0.7146 |
| flat.macro_recall | 0.6439 | 0.8389 | 0.6441 | 0.8412 |
| flat.macro_precision | 0.7506 | 0.8271 | 0.7506 | 0.8271 |
| flat.micro_recall | 0.5123 | 0.8430 | 0.5127 | 0.8465 |
| flat.micro_precision | 0.6454 | 0.7698 | 0.6454 | 0.7698 |
| flat.macro_gene_recall | 0.6236 | 0.5495 | 0.6236 | 0.5495 |

Post hoc (added after the run; changes no verdict): Descriptive, not a bar: A119b's `hierarchy.micro_precision` (0.1920 against 0.6454 for flat) is driven by ONE match. Truth family ID_256 (2 members) has no family-level partner in the flat score; in the primary A119b hierarchy it is matched to superfamily SF166 (2,991 members) at Jaccard 0.0006687 (recall 1.0, precision 0.0006687). Without that match the micro precision is 827/1326 = 0.6237. In the hierarchy every evaluable family gets a partner (`hierarchy.n_matched` 247 of 247 and 84 of 84, against 233 and 77 for flat), so the matched-based precisions and counts are not comparable between the two rows. Testis has no such match (its largest matched node has 13 members; micro precision 0.7543 against 0.7698 flat). `exact_any.any` (63 and 42) against `hierarchy.exact` (63 and 41): in testis, ID_453 and ID_454 are two Soto families with identical member sets that share one expansion node, and the one-to-one matching can give it to only one of them.

**Levels.** `exact_any` (tie-free: truth sets equal to some node of the level, evaluable families only), `hierarchy.by_level` (where the one-to-one matches land) and `hierarchy.exact_by_level`:

| | A119b (primary) | testis (primary) | A119b (shipped, secondary) | testis (shipped, secondary) |
|---|---|---|---|---|
| exact_any.any | 63 | 42 | 63 | 42 |
| exact_any.family | 59 | 39 | 59 | 39 |
| exact_any.expansion | 56 | 41 | 56 | 41 |
| exact_any.superfamily | 57 | 39 | 57 | 39 |
| exact_any.sd_complex | 9 | 5 | 9 | 5 |
| exact_any.expansion_not_family | 4 | 3 | 4 | 3 |
| exact_any.expansion_unique | 14 | 3 | 14 | 3 |
| hierarchy.by_level.family | 170 | 69 | 170 | 69 |
| hierarchy.by_level.expansion | 57 | 12 | 57 | 12 |
| hierarchy.by_level.superfamily | 20 | 3 | 20 | 3 |
| hierarchy.by_level.sd_complex | 0 | 0 | 0 | 0 |
| hierarchy.exact_by_level.family | 58 | 37 | 58 | 37 |
| hierarchy.exact_by_level.expansion | 5 | 4 | 5 | 4 |
| hierarchy.exact_by_level.superfamily | 0 | 0 | 0 | 0 |
| hierarchy.exact_by_level.sd_complex | 0 | 0 | 0 | 0 |

`exact_any.expansion_not_family` (Soto families equal to an expansion that is not a whole family; the only matches strictly below the family level) is 4 in A119b and 3 in testis; `exact_any.expansion_unique` (equal to a `<family>.unique` bucket) is reported separately, 14 and 3. No truth family lands on an SD complex in the one-to-one matching (`by_level.sd_complex` 0), although `exact_any.sd_complex` is 9 and 5.

Post hoc (added after the run; changes no verdict): each of those truth sets is also equal to a family node (the complex is one whole family; checked on the member sets), and the canonical placement gives it the lowest-rank node.

Post hoc (added after the run; changes no verdict): of the 5 (A119b) and 4 (testis) `exact_by_level.expansion` matches, 1 and 2 are canonical-placement twins of duplicate Soto families (the Soto family's member set equals another Soto family's set, and that other family holds the identical whole-family node): primary A119b ID_454 on `MCL1025.x1`, whose set equals family `MCL1025`, held by ID_453; primary testis ID_264 on `MCL202.x1`, whose set equals family `MCL202`, held by ID_262, and ID_481 on `MCL220.x1`, whose set equals family `MCL220`, held by ID_324; so the exact matches strictly below the family level are 4 in A119b (ID_112, ID_200, ID_477, ID_68) and 2 in testis (ID_147, ID_454); primary testis `exact_any.expansion_not_family` = 3 includes ID_453, which is not an exact match (its member set equals ID_454's, and the only node with that set, `MCL46.x2`, went to ID_454); and with the twin placements removed from the exact counts (same denominators) the exact rates are 62/247 = 0.2510 (A119b) and 39/84 = 0.4643 (testis), so the E2 difference is 0.2133, still > 0.10 (check script and output: `/mnt/linuxdisk/tmp/o1_close/hier/twin_check.py`, `twin_check.txt`).

Post hoc (added after the run; changes no verdict): **duplicate nodes.** The one-to-one matching runs over nodes, not over distinct member sets, and one member set often repeats across levels (for example a family with a single expansion, or a superfamily that is one family): the primary A119b hierarchy has 7,204 nodes but 3,452 distinct member sets, primary testis 1,261 nodes and 481 sets. Two truth families can therefore both match one member set through two of its nodes. In the scored runs 13 (A119b) and 6 (testis) member sets carry more than one match, 15 and 6 matches beyond the first on a set; 13 and 4 of these come from truth families with different member sets, the other 2 and 2 from duplicate Soto families. Matched over distinct member sets (one column per set, labelled by its lowest-rank node; the same Hungarian and canonical placement), A119b gives exact 62, mean Jaccard 0.5700, 245 matched and 50 landing on an expansion (as scored 63 / 0.5844 / 247 / 57), and testis 39 / 0.7488 / 81 / 8 (as scored 41 / 0.7836 / 84 / 12). The gain over flat shrinks from +5 to +4 exact and from +0.0557 to +0.0412 mean Jaccard in A119b, and from +4 to +2 exact and from +0.0707 to +0.0358 in testis: matching over distinct sets removes about a quarter (0.26) of A119b's mean-Jaccard gain and about half (0.49) of testis's. The duplicate Soto member-set groups (evaluable families with identical mapped member sets; A119b 4: ID_96 / ID_97, ID_162 / ID_347, ID_255 / ID_256, ID_453 / ID_454; testis 3: ID_262 / ID_264, ID_324 / ID_481, ID_453 / ID_454) explain the twin placements above and the ID_256 -> SF166 match of the Summaries note: in primary A119b, ID_256's set equals ID_255's, ID_255 holds family `MCL86` and ID_254 holds `MCL86.unique`, so the one-to-one matching leaves ID_256 the 2,991-member superfamily. Script `bench/o1_hier/posthoc_checks.py` (it reloads each scored JSON's own inputs and stops unless its re-match equals the scored placements); output `/mnt/linuxdisk/tmp/o1_close/posthoc_checks.json` (`distinct_sets`).

**Unconditioned numbers over ALL 491 truth families** (a family without a match counts 0; every truth gene counts, clustered or not):

| | A119b (primary) | testis (primary) | A119b (shipped, secondary) | testis (shipped, secondary) |
|---|---|---|---|---|
| hierarchy.n_truth_all | 491 | 491 | 491 | 491 |
| hierarchy.exact_of_all | 0.1283 | 0.0835 | 0.1283 | 0.0835 |
| hierarchy.macro_gene_recall_all | 0.3380 | 0.1034 | 0.3380 | 0.1034 |
| hierarchy.micro_gene_recall_all | 0.3470 | 0.0956 | 0.3470 | 0.0956 |
| flat.n_truth_all | 491 | 491 | 491 | 491 |
| flat.exact_of_all | 0.1181 | 0.0754 | 0.1181 | 0.0754 |
| flat.macro_gene_recall_all | 0.3137 | 0.0940 | 0.3137 | 0.0940 |
| flat.micro_gene_recall_all | 0.3104 | 0.0883 | 0.3104 | 0.0883 |

**Nulls** (1,000 rounds, seed 7; observed values are the hierarchy's; p = (1 + #{null rounds >= observed}) / 1001, quoted unrounded):

| | A119b (primary) | testis (primary) | A119b (shipped, secondary) | testis (shipped, secondary) |
|---|---|---|---|---|
| observed hierarchy.exact | 63 | 41 | 63 | 41 |
| observed hierarchy.mean_jaccard | 0.5844 | 0.7836 | 0.5846 | 0.7852 |
| null.global.exact_mean | 0.001 | 0.026 | 0.004 | 0.031 |
| null.global.exact_max | 1 | 2 | 1 | 2 |
| null.global.mean_jaccard_mean | 0.1934 | 0.2715 | 0.1923 | 0.2703 |
| null.global.mean_jaccard_max | 0.2109 | 0.2967 | 0.2105 | 0.2965 |
| null.global.p_exact | 0.000999000999000999 | 0.000999000999000999 | 0.000999000999000999 | 0.000999000999000999 |
| null.global.p_mean_jaccard | 0.000999000999000999 | 0.000999000999000999 | 0.000999000999000999 | 0.000999000999000999 |
| null.within_family.exact_mean | 61.048 | 40.569 | 61.071 | 40.522 |
| null.within_family.exact_max | 64 | 44 | 64 | 44 |
| null.within_family.mean_jaccard_mean | 0.5752 | 0.7793 | 0.5754 | 0.7808 |
| null.within_family.mean_jaccard_max | 0.5819 | 0.8007 | 0.5823 | 0.7988 |
| null.within_family.p_exact | 0.08291708291708291 | 0.48951048951048953 | 0.0939060939060939 | 0.4825174825174825 |
| null.within_family.p_mean_jaccard | 0.000999000999000999 | 0.21878121878121878 | 0.000999000999000999 | 0.2037962037962038 |

**Advisor table** (`bench/o1_hier/advisor_table.py`, as printed, 4 rows: the two primary catalogs, then the two shipped ones; TSV `docs/hierarchy_results/o1_close_human_2026-10-13.tsv`):

| sample | truth_families | evaluable | one_member | no_member | records_on_members | records_on_unclustered_locus | records_without_locus | records_without_member | flat_n_matched | hier_n_matched | flat_exact | flat_exact_rate | hier_exact | hier_exact_rate | hier_exact_of_all | flat_mean_jaccard | hier_mean_jaccard | flat_macro_recall | hier_macro_recall | flat_macro_precision | hier_macro_precision | hier_gene_recall | hier_gene_recall_all_macro | hier_gene_recall_all_micro | lands_family | lands_expansion | lands_superfamily | lands_sd_complex | exact_any_expansion | exact_expansion_not_family | exact_expansion_unique | p_global | p_within_family |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A119b | 491 | 247 | 55 | 189 | 1212 | 154 | 1092 | NA | 233 | 247 | 58 | 0.235 | 63 | 0.255 | 0.128 | 0.529 | 0.584 | 0.644 | 0.707 | 0.751 | 0.781 | 0.672 | 0.338 | 0.347 | 170 | 57 | 20 | 0 | 56 | 4 | 14 | 0.000999 | 0.082917 |
| testis | 491 | 84 | 34 | 373 | 287 | 99 | 2072 | NA | 77 | 84 | 37 | 0.44 | 41 | 0.488 | 0.084 | 0.713 | 0.784 | 0.839 | 0.928 | 0.827 | 0.835 | 0.604 | 0.103 | 0.096 | 69 | 12 | 3 | 0 | 41 | 3 | 3 | 0.000999 | 0.489510 |
| A119b_shipped | 491 | 247 | 54 | 190 | 1211 | 155 | 1092 | NA | 233 | 247 | 58 | 0.235 | 63 | 0.255 | 0.128 | 0.529 | 0.585 | 0.644 | 0.708 | 0.751 | 0.781 | 0.672 | 0.338 | 0.347 | 170 | 57 | 20 | 0 | 56 | 4 | 14 | 0.000999 | 0.093906 |
| testis_shipped | 491 | 84 | 34 | 373 | 286 | 100 | 2072 | NA | 77 | 84 | 37 | 0.44 | 41 | 0.488 | 0.084 | 0.715 | 0.785 | 0.841 | 0.931 | 0.827 | 0.835 | 0.604 | 0.103 | 0.096 | 69 | 12 | 3 | 0 | 41 | 3 | 3 | 0.000999 | 0.482517 |

records_without_member is a --map key column; it is NA in these exon-mapping runs.

Post hoc (added after the run; changes no verdict): **reading the advisor table.** `p_global` and `p_within_family` are the `p_exact` values only; the `p_mean_jaccard` values are in the Nulls table (A119b's within-family `p_mean_jaccard` is 0.000999, from the JSON). `exact_any_expansion`, `exact_expansion_not_family` and `exact_expansion_unique` are tie-free set equality (`exact_any`), not one-to-one matches: in primary testis, ID_453 equals `MCL46.x2` as a set, but its one-to-one match is family `MCL46` (Jaccard 0.6), because `MCL46.x2` went to ID_454, whose set is identical.

**Per-truth tables.** Each JSON carries `per_truth` (one record per evaluable truth family: `label`, `n_records`, `n_members`, `level`, `node`, `exact`, `jaccard`, `recall`, `precision`, `gene_recall`, `flat_node`, `flat_jaccard`) and the `inputs` it was called with:
- primary A119b `/mnt/linuxdisk/tmp/o1_close/hier/human_A119b/hier_match.json`; primary testis `/mnt/linuxdisk/tmp/o1_close/hier/human_testis/hier_match.json`;
- shipped A119b `/mnt/linuxdisk/tmp/o1_close/hier_shipped/human_A119b/hier_match.json`; shipped testis `/mnt/linuxdisk/tmp/o1_close/hier_shipped/human_testis/hier_match.json`.

### Secondary (the shipped pre-fix catalogs; same statistics, no bars, no verdicts)

The Secondary columns of the tables above are the shipped catalogs (A119b 9,907 members / 2,150 families; testis 1,115 / 342), built with `mcl_families.before` from the same hierarchy inputs and scored with the same truth, scorer and seed. They differ from the primary in the catalog sizes, in a few coverage counts, and in ONE evaluable truth family per catalog; no exact count, level count or global-null p-value moves. The same families are evaluable (247 and 84); exact (63 and 41, flat 58 and 37), `n_matched`, `by_level`, `exact_by_level`, `exact_any` and both global-null p-values are identical, and no exact flag flips. The family that changes is ID_214 in A119b (7 mapped members in the primary, 6 in the shipped catalog; Jaccard 0.2857 in the primary against 0.3333 shipped) and ID_391 in testis (5 against 4; 0.6667 against 0.8). That moves mean Jaccard by +0.0002 (A119b) and +0.0016 (testis) in the shipped catalogs, and recall slightly (A119b macro recall 0.7073 against 0.7075, testis 0.9281 against 0.9305). The coverage counts that move do so by 1 each: records on members (A119b 1212 against 1211, testis 287 against 286), records on unclustered loci (154 against 155, 99 against 100), and in A119b one label (ID_373) with one member in the primary and none in the shipped catalog (`one_member` 55 against 54, `no_member` 189 against 190). The nulls' own values move a little (within-family `p_exact` 0.0829 in the primary against 0.0939 shipped in A119b, 0.4895 against 0.4825 in testis; within-family `p_mean_jaccard` in testis 0.2188 against 0.2038; in A119b it is 0.000999 in both). The exact-rate difference between the two shipped catalogs is the same 0.2330 as in the primary. As the prereg says, the pair-rule choice must not be made from these scores.

Post hoc (added after the run; changes no verdict): on 2026-10-13 the user made `exact-mcl` the default pair rule (`docs/PREREG_mcl_two_node_2026-10-13.md`, Results). The choice was made on catalog identity and on MCL's consistency across component sizes, not on these scores; the two arms give the same exact counts (63 / 41). The Secondary arm is therefore the default-rule result.
