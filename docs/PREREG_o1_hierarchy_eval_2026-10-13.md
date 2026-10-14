# PREREG — O1 hierarchy-aware evaluation, human (2026-10-13)

Spec: `docs/superpowers/specs/2026-10-13-close-o1-design.md` section A; plan
`docs/superpowers/plans/2026-10-13-close-o1-phase-0A.md`.

**Question.** Are Soto 2025's families nodes of our family hierarchy? The flat one-to-one match against our MCL
families counts every Soto family that sits one level down as a miss.

**Catalogs.**
- Dev: human A119b. Held-out: human testis.
- Both are default catalogs, re-clustered with the MCL two-node fix (`docs/PREREG_mcl_two_node_2026-10-13.md`).

**Hierarchy.**
- Command: `mcl_families.after --paf <S>.fam.loci.paf --gff <S>.fam.loci.gff3 --min-exonic-bp 1
  --min-shared-exon-frac 0.60 --emit-hierarchy --hierarchy-duplicons lift/dupmasker_colors.v2.bed --hierarchy-sd
  lift/sedef_pairs.v2.bed`.
- Inputs: DupMasker colors and SEDEF `final_v1_clean.bed`, lifted from CHM13 v1.0 by the T2T chain (sequence-verified).
- No `--hierarchy-external`: the truth never enters the build.
- Its families must equal the re-clustered catalog (`cmp`).

**Truth.** `lift/soto_truth.v2.tsv`:
- S1C families proper (491; `Unassigned_*` excluded), genes joined to CAT v4 by Gene ID;
- exons = the union over the gene's transcripts, chain-lifted to v2.0.

**Mapping.**
- A Soto gene joins a member when one of its exons shares >= 1 bp, on the same strand, with an exon of the member's
  locus or of a locus folded into it.
- Exons come from the catalog's `.fam.loci.gff3`; the fold map is the run's `ours.loci.tsv`.
- A family is EVALUABLE with >= 2 members. The others are reported as coverage.

**Scorer.** `bench/o1_hier/hier_match.py --map exon --null 1000 --seed 7`.
- One-to-one Hungarian (scipy) on -Jaccard between evaluable truth families and every family, expansion,
  superfamily and SD-complex node.
- Tie policy: rows in label order; columns in (level rank family < expansion < superfamily < sd_complex, node id) order;
  then canonical placement: the truth families matched to nodes with one identical member set take that set's lowest-rank
  nodes, best Jaccard first, then label (`Scorer._prefer_lower_rank`). Identical sets give identical Jaccards, so the
  optimum and every number except the reported node and level are unchanged.
- Nulls, both keeping node sizes and nesting:
  - **global**: one random permutation of all members per round;
  - **within_family**: members permuted inside their own family. Families, superfamilies and SD complexes are
    unchanged; only expansions move.

**Secondary (registered; same statistics, no bars).**
- **Why.** `docs/PREREG_mcl_two_node_2026-10-13.md` (F2) showed that the residue bug split no pair in the shipped
  catalogs: all 6,969 tie-regime pairs were already joined. Fix 0 therefore only adds pairs below w* = 0.456885
  (A119b +68 families, testis +9). That is a definition change whose adoption is the user's choice (Amendment 1 there).
- **What.** The shipped pre-fix catalogs, which equal exact-MCL semantics on these graphs, are also scored. They use the
  same hierarchy inputs, truth, scorer and seed, and `mcl_families.before`. Their families must equal the shipped
  `.fam.clusters.tsv`.
- **How reported.** Next to the primary, with no bar. The pair-rule choice must not be made from these scores.

**Bars** (`exact` = Hungarian matches with Jaccard 1; `rate` = exact / evaluable):
- **E1:** on each catalog, `null.global.p_exact` < 0.001 AND `null.global.p_mean_jaccard` < 0.001.
- **E1b:** on each catalog, `null.within_family.p_exact` < 0.001 AND `null.within_family.p_mean_jaccard` < 0.001.
- **E2:** |rate(testis) - rate(A119b)| <= 0.10.
- p-values are reported to 6 dp: with 1,000 rounds the smallest is 1/1001 = 0.000999, which passes p < 0.001 only when no
  null round reaches the observed value.

**Reading, fixed now:**
- **E1 + E1b + E2 pass:** "Soto's families are nodes of our hierarchy". The O1 headline for Soto becomes, per catalog,
  the exact rate and mean Jaccard, next to the flat score and the coverage.
- **E1 + E2 pass, E1b fails:** the claim stops at the family level and above. Our expansions recover Soto no better
  than random within-family splits of the same sizes.
- **E2 fails:** sample-specific; only the per-catalog numbers are reported.
- **E1 fails:** no claim.

**Reported regardless:**
- coverage (records on members / on unclustered loci / without a locus; labels evaluable / one member / none; members
  in several truth families);
- the hierarchy and flat summaries (exact, rate, mean Jaccard, macro and micro recall / precision, gene recall);
- `exact_any` per level, and where matches land (`by_level`);
- the unconditioned numbers over ALL 491 truth families (`hierarchy.exact_of_all`, `macro_gene_recall_all`,
  `micro_gene_recall_all`; a family without a match counts 0, and every truth gene counts, clustered or not);
- `exact_any.expansion_not_family`: Soto families equal to an expansion that is not a whole family, the only matches
  strictly below the family level. `exact_any.expansion_unique` (equal to a `<family>.unique` bucket) is reported
  separately.
- the per-truth table;
- the advisor table (`bench/o1_hier/advisor_table.py`).

**Disclosures.**
- E1b, the exon + strand mapping and the >= 2-member evaluable rule were fixed after the scorer's regression run on the
  10-09 F_O universe (v1.0, guided). There the global null averaged 0.03 exact matches; the within-family null averaged
  224 of the observed 231; span-overlap mapping cut exact expansions from 215 to 169.
- No de novo catalog has been scored by this scorer before this prereg.
- A non-Soto member inside a matched node counts against precision. Soto's truth covers only its SD98 gene set, so
  precision is a lower bound.
- the unconditioned metrics, the expansion-level split, canonical placement and the 6-dp p-values were added in the
  scorer's review rounds. They were checked on the 10-09 F_O universe only:
  - every earlier F_O number was unchanged, except where exact matches land (`exact_by_level` 13 / 228 / 4 instead of
    14 / 227 / 4 with key mapping), because an exact match now takes its identical family node;
  - `expansion_not_family` = 9 (exon mapping) of the 206 expansion matches.

## Results

(pending)
