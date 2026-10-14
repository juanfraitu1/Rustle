# Closing O1 — design (2026-10-13)

## Intent (agreed with the user)

**The user's words:** "focus on closing O1"; "the bipartite matching is not really showing what we want"; a Hungarian
matching "could be useful"; and the other session's SEDEF port (`/mnt/c/Users/jfris/Desktop/sd/sedef-rs`) should
"derive an unbiased rule for how sedef works, or [be] a dependency of our pipeline and thus not really an external
entry".

**Diagnosis.**
- Today's "bipartite" score is already a Hungarian match (scipy `linear_sum_assignment`), but between our FLAT MCL
  families and Soto's families.
- Soto's families sit one level DOWN in our hierarchy: 215 / 400 are exactly one of our expansions, and 78% lie inside
  one of our families (`project_hierarchy_flag_2026-10-09`).
- So the flat one-to-one match counts every split as misses: recall about 0.49, F about 0.63.

**Decisions** (each chosen by the user):
1. **Evaluation:** hierarchy-aware matching. Truth families are matched one-to-one (Hungarian) to nodes at ANY level
   of our hierarchy, against a label-shuffle null.
2. **`sedef-rs`:** becomes the pipeline's internal SD + duplicon layer, vendored into Rustle. It replaces DupMasker
   (human only) only after validation.
3. **Ape truth:** Soto's families projected by Liftoff, plus the Yoo 2025 lineage-specific expansions.
4. **Sequencing A:** the human evaluation first (on the current DupMasker hierarchy), then `sedef-rs`, then the apes.
5. **Included:** the MCL two-node fix. 57% of genome-wide families have two members, and today their existence is
   decided by floating-point residue (`project_mcl_two_node_fragility`).

Each component below is pre-registered (`docs/PREREG_*.md`) before it runs. Bars are threshold-free where possible.
**Order:** 0 (MCL fix) -> A (scorer + human) -> V (`sedef-rs`) -> B (apes) -> deliverables. The scorer can be built
while 0 is reviewed, but every scored run uses the re-clustered catalogs.

## 0. MCL two-node fix

**Bug.** On a connected 2-node graph, MCL's column converges toward a tie between "self" and "mate". The convergence
test (`|u - v| < 1e-7`) stops at a residual whose SIGN picks the attractor: join or split, depending on the weight
(0.99 splits, 0.98 joins, 0.90 splits, ...).

**Fix.** A connected component of exactly two nodes is one cluster, decided before MCL runs. Three or more nodes are
unchanged (they always join today).
- Location: `annotation_families.rs::mcl`, or its caller in `mcl_families`.
- Red-first tests: the bit-faithful probe weights, and the famsim ladder rungs snp_0.01 / snp_0.1, which failed in both
  modes and both species.

**Re-clustering.** The four default catalogs (human A119b, human testis, chimp, gorilla KB3781) are re-clustered from
their cached `.fam.loci.paf` / `.gff3`. Edges are not recomputed; the families stage is not rerun.

**Bars.**
- Only size-2 components change.
- Every other family is byte-identical.
- The famsim rungs pass.
- The number of families gained or lost is reported per catalog.

## A. Hierarchy-aware scorer and the human evaluation

**Scorer** (`bench/o1_hier/`, Python, unit-tested).
- **Inputs:**
  - the hierarchy flag's tables (`<out>.hierarchy.{members,nodes,edges}.tsv`), from which each node at every level
    (SD complex, superfamily, family, expansion) becomes its set of member loci;
  - the truth: families as sets of genes, mapped to our loci by coordinate overlap, never by name.
- **Coverage:** truth genes that overlap no locus of ours are reported separately as catalog coverage. Recall is never
  conditioned on the prediction.
- **Matching:**
  - one-to-one maximum-weight Hungarian between truth families and all nodes;
  - weight = member-level Jaccard;
  - tie policy named and deterministic (scipy order on a stable node sort);
  - no non-overlap constraint, because the hierarchy is a partial order (expansions may share fragments). The null
    controls depth instead.
- **Per truth family:** the matched node and its level, Jaccard, recall, precision, and exact (Jaccard = 1).
- **Headline:**
  - macro / micro recall and precision, the exact count, and mean Jaccard;
  - where truth lands (truth families matched per level);
  - next to it, today's flat family-level bipartite score on the same data.
- **Null:** 1,000 shuffles of member labels across the hierarchy, keeping node sizes and nesting. Gives the expected
  score and a p-value.

**Human run.**
- **Catalogs:** dev = human A119b; held-out = human testis. Both are default catalogs, after fix 0.
- **Hierarchies:** `mcl_families --emit-hierarchy` with the 10-09 inputs: DupMasker colors
  (`winloci_data/duplicons/chm13.draft_v1.0_plus38Y_dupmasker_colors.bed`) and the CHM13 SEDEF calls the 10-09 run used
  (`winloci_data/soto_replication/final_v1_clean.bed`; corrected from `HSA_sedef_pairs.bed`, whose column 21 is not an
  identity).
  - Both are v1.0 coordinates. The 10-09 run did not lift them: its universe was v1.0.
  - Here they are lifted by the T2T v1.0 -> v1.1 chain, sequence-verified (plan 1, Task 3).
- **Truth:** Soto's 400 families (`/mnt/linuxdisk/tmp/hier_flag/soto_external.tsv`).
- **Bars:**
  - **E1:** the hierarchy-aware score beats the null at p < 0.001 on both catalogs.
  - **E2:** the held-out exact-recovery rate is within 0.10 of the dev rate.
- **Reading, fixed in advance:**
  - E1 and E2 pass: "our hierarchy contains Soto's families" replaces the flat F of about 0.63 as O1's headline.
  - E2 fails: the result is sample-specific and is not claimed.

## V. `sedef-rs` in the pipeline, and its validation

**Packaging.** The `sedef-rs` source (seed, chain, refine, decompose, anchors, kmer, fasta) is copied into Rustle as
`src/sedef/`.
- The Desktop repo stays the user's sandbox; changes are ported deliberately.
- New command: `utilities sd-scan --genome FA --out PREFIX [-t N]`, producing
  `PREFIX.sd.{pairs,segments,families}.tsv`.
- Pipeline stage `sd`: one run per genome, cached by genome hash.
- Converters to the hierarchy flag's existing inputs: an SD-pairs BED, and a duplicon BED (segments labelled by
  duplication family). `mcl_families --emit-hierarchy` is unchanged.

**V1 (faithful SD calls).** Whole-genome runs on CHM13 v2.0 and on mGorGor1 (`GGO.fasta`), compared base-pair-wise
with the SEDEF calls on disk:
- CHM13: `winloci_data/soto_replication/final_v1_clean.bed` (v1.0 SEDEF), lifted by the T2T chain (plan 1, Task 3);
- `winloci_data/GGO_sedef_final.bed`.

Bar: SD recall and precision >= 0.80 on both genomes. Fail: `sedef-rs` is not adopted.

**V2 (can it replace DupMasker for expansions).** The human hierarchy is rebuilt with `sedef-rs` segments as the
duplicon layer and compared with the DupMasker-based hierarchy, using the measures our own duplicon units failed on
10-10 (register 1263):
- expansion-level bipartite agreement >= 0.80 on a held-out family-hash half;
- Soto families recovered as exact expansions no more than 10% below DupMasker's count.

Fail: `sedef-rs` gives the SD level only, and the expansion level stays human-only.

## B. Gorilla and chimp

- **Catalogs:** gorilla KB3781 and chimp, after fix 0.
- **Hierarchies:** built with the `sd` stage on each catalog's assembly; with the SD level only if V2 fails.
- **Truth 1 (Liftoff-Soto):** Soto's genes projected by Liftoff onto each assembly, extra copies included, reusing the
  other session's orangutan setup (`bench/hierarchy/liftoff_cross.py`). Overlapping projections are scored as a cover.
- **Truth 2 (Yoo 2025):** the lineage-specific expansions from `winloci_data/yoo_2025`, scored at the expansion level.
- **Bars:**
  - **E3:** beats the null at p < 0.001 on each ape against Liftoff-Soto.
  - **Reported:** exact and Jaccard next to human's, Yoo expansions recovered exactly, and the coverage of each truth.
- **Caveats, stated:**
  - tissue and read set differ by sample, so differences from human mix sample and species;
  - a missing truth family can be a Liftoff failure, so coverage is reported separately.

## Deliverables

- the MCL fix with its tests, and the re-clustered catalogs;
- the scorer with its tests;
- `src/sedef/`, `utilities sd-scan`, the `sd` stage and the converters;
- pre-registrations and results for 0, A, V and B;
- register rows, and `docs/THESIS_OBJECTIVES.md`'s O1 line updated to the claim the bars support;
- one advisor-facing table: flat vs hierarchy-aware, the level where each truth family lands, and coverage per sample.

## Not in scope

- O2 / O3;
- the paused genome-wide figure 7 scoring;
- cross-species consistency of duplicon IDs (each genome gets its own `sedef-rs` families);
- a non-overlap (antichain) matching variant: NP-hard on the partial order, so the null controls depth instead.

## Amendments from the plan's fact check (2026-10-13, before any scored run)

Recorded in `docs/superpowers/plans/2026-10-13-close-o1-phase-0A.md` ("Rulings").
1. **Truth mapping is exon + strand.** A Soto gene joins a locus when one of its exons shares >= 1 bp, on the same
   strand, with an exon of the locus. Raw span overlap joined antisense and intronic neighbours: on the 10-09 F_O data
   it cut exact expansions from 215 to 169; exon + strand gives 206.
2. **Evaluable truth family = >= 2 members on our loci** (the 10-09 convention). E2's rate is exact / evaluable;
   unconditioned numbers are reported next to it.
3. **E1b is added**, a within-family null: members permuted inside their own family, so only the expansions move.
   - The global null of E1 is passed trivially (0.03 exact matches on F_O).
   - E1b is the test that controls depth. If it fails, the claim stops at the family level.
4. **The pair rule of fix 0 runs after MCL.** This is equivalent, since flow never crosses components. A pair's fate
   also depended on unrelated components, because convergence is tested on all columns at once.

Amendments 5-7 were added on 2026-10-13 after the human evaluation's run and its final review
(`docs/PREREG_o1_hierarchy_eval_2026-10-13.md`, Results). They change no verdict.

5. **Canonical placement replaces "scipy order on a stable node sort"** (section A, Matching, tie policy). The truth
   families matched to nodes with one identical member set take that set's lowest-rank nodes, best Jaccard first, then
   label (`Scorer._prefer_lower_rank`, `bench/o1_hier/hier_match.py` at 8f0515ef). scipy's choice among equal-cost
   columns is not guaranteed.
6. **The diagnosis "Soto's families sit one level DOWN" did not hold.** Of the 215 Soto families that equal one of our
   expansions on F_O (10-09 run, key mapping), 10 lie below a whole family (`exact_any.expansion_not_family`); the other
   205 also equal a whole family. On the de novo catalogs, 4 (A119b) and 2 (testis) exact one-to-one matches lie
   strictly below the family level.
7. **Phase B prerequisites.** Phase B must:
   - match over distinct member sets (one column per set), since identical sets repeat across levels;
   - use the within-family null for E3, since the global null is passed by design;
   - base V2's bar (Soto families recovered as exact expansions) on `expansion_not_family`;
   - stratify cross-catalog bars by mapped set size.
