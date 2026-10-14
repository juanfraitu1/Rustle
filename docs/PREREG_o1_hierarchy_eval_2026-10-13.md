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
  - Output: `/mnt/linuxdisk/tmp/o1_close/hier_shipped/<S>`.
  - Their families must equal the shipped catalog:
    `cmp /mnt/linuxdisk/tmp/o1_close/hier_shipped/<S>/ours.clusters.tsv $R/<S>.fam.clusters.tsv`.
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
   - The headline rate is (`hierarchy.exact` - `exact_by_level.expansion`) / `hierarchy.n_truth`. Quote it with
     `exact_by_level` and `by_level` excluding `expansion`, next to the flat score.
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
- **Levels:** `exact_any` per level, and where matches land (`by_level`, `exact_by_level`).
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

(pending)
