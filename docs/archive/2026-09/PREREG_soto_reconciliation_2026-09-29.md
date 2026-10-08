# Pre-registration: reconciling our pipeline with Soto 2025 — which named choices close the gap (attribution)

**Written 2026-09-29 (KEY=soto_reconcile) BEFORE any new configuration was computed or scored.** Human CHM13 v1.0 only;
Soto's genome, CAT v4 annotation, gene universe (Table S1C, 2,334 genes / 491 multigene families) and S1C famCN.
Everything here is **concordance with Soto's own tables, not an independent replication** (register 858 / 1085).
Nothing in `src/` or `bench/` is edited; nothing is committed. Analysis code lives under
`/mnt/linuxdisk/tmp/rustle_figures_dev/soto_reconcile/`. This file binds once Amendment 1 records its sha1 and the
instruments' sha1s.

## 0. The question

The replication headline is Soto's literal recipe on native CHM13 v1.0 (**BASE**): ARI 0.7096 / 235 of 491 exact
(median MAD), 0.7039 / 261 (mean MAD) (`docs/archive/2026-09/SOTO_REPLICATION_STATUS_2026-09-28.md` §1.1). This is NOT a search for a
better method. It asks: **which documented choices, set the way Soto's own paper, released code or inputs set them,
move our output to Soto's own family table, and how far?** The answer must read "BASE x → chosen y (held-out z);
remaining gap = a (cover) + b (their rule / manual steps) + c (unreachable by any ≥ 98% exon link) + d (other)".
No continuous parameter is tuned against the truth; there are no per-family choices.

## 1. What was learned today BEFORE this file, from Soto's released code (no configuration scored)

Source: `github.com/mydennislab/HSD_brain_evolution`, fetched 2026-09-29 to `soto_reconcile/soto_code/`
(`A_SD98_regions.md` sha1 736cb4c2, `B_SD98_families.ipynb` sha1 9365476b). The 09-11 trace (§6if) read only the
notebook; `A_SD98_regions.md` had never been read ("not yet read", memory `project_soto_full_replication`). It changes
three readings:

1. **The map-back queries are SD98 EXONS, not SD98 regions.** A §3.1 lines 226-238: exons of every CAT v4 gene
   (all biotypes: the `grep "protein_coding\|unprocessed_pseudogene"` on line 231 is commented out) fully inside a
   merged autosomal SD98 region (`bedtools intersect -wa -f 1`, chrX/chrY/chrM excluded), minus `genes_to_remove.txt`
   (their 71 curated genes), are extracted with `bedtools getfasta -s -name` and mapped with
   `minimap2 -c --end-bonus 5 --eqx -N 50 -p 0.5 -t 64` to the v1.0 genome. The paper's STAR Methods instead say
   "DNA sequences of all SD98 regions were extracted using BEDTools getfasta and mapped back" (fulltext line 1041).
   Our BASE follows the prose (unmerged SD98 units) plus our own CIGAR "projection reading" (the literal region reading
   scores ARI 0.19).
2. **Shared exon = same-strand ≥ 99% cover of an SD98 exon; self = same gene.** A lines 248-263: `bedtools intersect
   -f 0.99 -s -wao -a <SD98 exons> -b <PAF mappings>`; each hit gives the pair (target exon's gene, query exon's gene),
   dropped only if the two gene names are equal (`if($1!=$2)`), and kept if either member is protein-coding /
   unprocessed pseudogene (`grep "protein_coding\|unprocessed_pseudogene"`). Output: one `a,b` pair per line,
   `CHM13.combined.v4.exons.SD98.clusters.tsv`. (BASE removes self-mappings by coordinates, ± 50 bp of the unit.)
3. **Soto's MAD < 1 is applied to each shared-exon PAIR, then families grow through coding genes.** Notebook cell 4
   reads `data/SD98_exon_clusters.txt` as comma-separated lists; cell 7: "We selected **pairwise clusters** based on
   shared exons only if genes had similar overall WSSD copy number"; cell 8 keeps a cluster iff
   `stats.median_abs_deviation(per-gene median famCN) < 1`; cells 10-11 grow a family from every coding / unprocessed
   gene, merging every low-dispersion cluster that contains a coding / unprocessed member already in it; cell 14 only
   REPORTS each family's MAD. For a pair, median-MAD = mean-MAD = |ΔfamCN| / 2, so the gate is |ΔfamCN| < 2 under
   either statistic; a pair with one gene lacking WSSD (non-coding) has MAD 0 and passes. Consequences, all visible in
   S1C: a non-coding gene joins EVERY family it pairs with (the cover: 149 genes, 0 of them coding/unprocessed —
   checked today); a family's own MAD is never gated (the 12 S1C families with Family MAD ≥ 1: ID_8, 22, 28, 35, 41,
   63, 76, 78, 84, 113, 283, 330). The §6if port (`soto_replication.py dennislab`) read "clusters" as connected
   components (ARI 0.4966 on BASE edges, `soto_v1_diff.md` §4): that reading is superseded and is reported only as a
   diagnostic row, not searched.
4. **A second manual step is in S1C itself:** `Family MAD` = "Manual merge" for ID_347 (DUX4 / DUX4L / MIR8078, 46
   genes) and ID_482 (UBTFL1/2/3/5, 4 genes). No rule can be expected to emit these; they are a named wall (§8).
5. **Provenance of the tables asked about:** `winloci_data/soto_replication/shared_exons_1793_final.tsv` and every
   `replicated_families_*` (ids `SEDEFFAM*` / `DENNISFAM*`), and `bench/soto/shared_exons_2334_{finalhuman,
   finalv1_native}.tsv`, are OUR pipeline outputs (09-11 / 09-28). Soto's repository releases no edge, cluster or family
   table (its `data/` holds Tajima's D, pHSD VCFs and zebrafish data; the notebook's inputs `data/SD98_exon_clusters.txt`
   and `data/SD98_WSSD.tsv` are not in it). Soto's own tables are only S1C (`bench/soto/soto_famCN_S1C.tsv`, sha1
   d008a179) and S1E (`soto_parCN_S1E.tsv`).

## 2. Not blind (disclosed)

- The frozen DEV / HELD-OUT split of `PREREG_soto_losses_2026-09-29.md` (`soto_losses/frozen/split.tsv` 49bcbcfe;
  DEV 225 / HELD-OUT 266 families) is reused. **The HELD-OUT half was already scored on 09-29** for BASE, B1, C, ALL
  and readthrough control A (both statistics, 20 null seeds): e.g. BASE 0.8183 / 131 exact (median), C 0.8520, B1
  0.8194. So cells (REGION or REGION+B1) × (GREEDY or DP) × (median or mean) have KNOWN held-out values; B1 and C were
  developed on DEV. The 09-28 per-family diff (union arm, CIGAR arm, BASE) was read on all 491 families.
- The CIGAR edge source (0.6985 / 0.6894 on all 491) and the union (0.7098 / 0.7053) were scored on all 491 on 09-28.
- The curation rule was fitted on all 1,864 genes (both halves) on 09-28.
- NEW, never scored anywhere: the EXON edge source, the PAIR family rule, the "both" statistic, and every cell they
  enter. Only edge / gene COUNTS of the EXON edge set will be looked at before the freeze (Amendment 1), never its
  agreement with S1C.

## 3. The levers (every one justified before scoring)

Labels: **SOTO-CODE** (their released code), **SOTO-PROSE** (their STAR Methods), **INPUT** (a documented input
difference), **OURS** (not supported by their text or code; included because it is an established arm or asked for).

| lever | level | label | justification (quote / line) |
|---|---|---|---|
| **E edge source** | REGION (BASE) | SOTO-PROSE + OURS reading | STAR Methods l.1041 "DNA sequences of all SD98 regions were extracted ... mapped back ... For each SD98 exon, the BEDTools intersect with -f 0.99 ... removing self-mappings"; unmerged SEDEF ≥ 0.98 units of `final_v1.bed`; exon projected through the `cg:Z` CIGAR (our reading, register 946); frozen `soto_losses/frozen/base_edges.tsv` (11,935) |
| | EXON | SOTO-CODE | A §3.1 l.226-263 (quoted in §1.1-1.2); construction in §4 |
| | CIGAR | INPUT | SEDEF's own pairwise alignment instead of a map-back (§6ie / `soto_v1_cigar.md`); `bench/soto/shared_exons_2334_finalv1_native.tsv` (4,384) |
| | UNION | INPUT | REGION ∪ CIGAR (the 09-28 diagnostic arm: map-back plus SEDEF's own ≥ 98% alignments; 12,001 expected) |
| | B1 | OURS | REGION + the CIGAR-walk edges of SD98 rows whose two sides overlap (`soto_losses` fix B1, 38 edges; NOT DISTINGUISHABLE FROM ITS NULL on held-out, register draft 1152). Its stated mechanism (`-p 0.5` against a unit's own self-hit) exists only because BASE maps regions; exon queries do not create it |
| **F family rule** | GREEDY (BASE) | SOTO-PROSE reading + OURS split | l.1041 "groupings where the mean absolute deviation of the CN was less than one were selected": connected components, a component failing MAD < 1 split by our ascending-famCN greedy walk (Soto never states how a failing grouping is split); islands (§6im) and majority attach (§6ii) as in BASE `cluster --full-geneset` |
| | DP | OURS | `soto_losses` fix C (fewest famCN-contiguous MAD < 1 groups; ties least L1 dispersion about the group centre, then earliest cut), frozen code, unchanged tie-break. **Disclosed: C failed its pre-registered mean-arm no-harm clause on held-out 09-29 (−0.0175 ARI, one family ID_356, via this tie-break).** No new tie-break is introduced (it would be chosen after seeing ID_356) |
| | PAIR | SOTO-CODE | notebook cells 4-11 (§1.3): a shared-exon pair is kept iff MAD(per-gene famCN over its members that have one) < 1 (for two coding genes: \|ΔfamCN\| < 2; a pair with one famCN is kept; with none, dropped); families = connected components of the kept coding/unprocessed–coding/unprocessed pairs; every non-coding gene joins every family it has a kept pair with (cover); no islands, no split, no family-level gate |
| **S MAD statistic** (GREEDY, DP only; PAIR is invariant) | median (BASE) | SOTO-CODE | notebook cell 6 `stats.median_abs_deviation` (unscaled) |
| | mean | SOTO-PROSE | l.1041 "mean absolute deviation" |
| | both | OURS | a component / group passes only if median-MAD < 1 AND mean-MAD < 1. The requested "recursive application of MAD < 1 to sub-half minorities": neither their text nor their code recurses, so any recursive rule is OURS; this conjunction is the constant-free stand-in (a mean deviation sees a far sub-half minority that a median cannot). The Soto-supported way below the component is PAIR's per-pair gate |
| **K curation** (NOT searched; BASE already = Soto) | Soto's hand-curated list (S1C `In Table S1` 1,793 eligible) | SOTO (their table) | l.1037 "performed manual curation ... removing redundant and read-through fusion transcripts"; A l.99 "we landed on 71 genes" |
| | rule@0.3 (ablation only) | OURS reconstruction | `soto_replication.py curate --rule contained_outranked --contain 0.3` (`soto_cur_rule.md`; indistinguishable from random pruning, `soto_cur_critique.md`) |
| **U node universe** (NOT searched; BASE = Soto) | 2,334 S1C genes, non-coding members attachable | SOTO-PROSE | l.1041 "SD98 genes associated with other gene features, including lncRNAs and processed pseudogenes, were also assigned a gene family ID" |
| | 1,793 eligible only (ablation only) | OURS (pre-§6ih history) | contradicts the sentence above |
| -p / -N | 0.5 / 50 | SOTO (both) | already their stated values in every arm; no lever (a `-p 0` level would move AWAY from Soto; DEV diagnostic 09-29: +124 edges at 0.40 same-family precision in 25 units) |
| minimap2 version | 2.30 (ours) vs 2.17 (theirs) | INPUT | sensitivity on the EXON edge set only, if the v2.17 release binary runs here; reported, never searched |
| genome | v1.0 without chrY (ours) vs `plus38Y` (theirs) | INPUT | disclosed, not testable (no v1.0 chrY FASTA used here); affects only `-N`/`-p` competition from chrY copies |

**Search space (full factorial, 35 cells ≤ 64):** E ∈ {REGION, EXON, CIGAR, UNION, B1} × [ F ∈ {GREEDY, DP} ×
S ∈ {median, mean, both} (30 cells) + F = PAIR (5 cells; statistic-invariant) ]. **Soto-only subspace (6 cells):**
E ∈ {REGION, EXON} × {GREEDY-median, GREEDY-mean, PAIR}. Every cell is reported (table file), both halves and all.

## 4. EXON edge construction (binding specification)

1. SD98 exons: every block of every `cat_v4.bed` transcript (BED12, gene id = field 19, biotype = field 20) that lies
   fully inside one interval of `winloci_data/soto_replication/sd98_v1.bed` (the merged UCSC v1.0 `$24 >= 0.98`
   regions, 97,797,568 bp) on an autosome; deduplicated on (chrom, start, end, strand, gene id). Probe count (today,
   counts only): 25,971 exons of 5,154 genes (= their line-172 count of autosomal SD98 genes).
2. Query sequences: `samtools faidx` of `t2t-chm13-v1.0.fa.gz`, reverse-complemented for `-` exons (`getfasta -s`).
3. `minimap2 (2.30) -c --end-bonus 5 --eqx -N 50 -p 0.5 -t 5` against an index of the same FASTA built with
   `minimap2 -d` (default parameters = the no-`-x` default, as theirs). Every PAF line (primary or secondary) is a
   mapping (target chrom, tstart, tend, PAF strand).
4. `bedtools intersect -f 0.99 -s -wa -wb -a <SD98 exons, strand = gene strand> -b <mappings, strand = PAF strand>`;
   each hit → pair (gene of the target exon, gene of the query exon); drop if the two gene ids are equal; keep if at
   least one gene's biotype is in {protein_coding, unprocessed_pseudogene, transcribed_unprocessed_pseudogene,
   translated_unprocessed_pseudogene} (their substring grep).
5. Removal of their 71 curated genes (`genes_to_remove.txt`) is applied as an edge filter afterwards (equivalent: each
   query maps independently, each exon × mapping overlap is reported independently); the node universe filter (U)
   then keeps edges with both ends in the universe. Per-gene biotype, not per-transcript (their exon names carry a
   `biotype=` field whose level we cannot see) — disclosed.
6. Heavy steps (index build, mapping) under `tools/rlock.sh heavy`, foreground, shards ≤ 9 min.

## 5. Scoring, selection, null

- Scorer: `soto_losses/frozen/lib.py` `score_half` on the clean (single-family) S1C genes of the half (DEV 1,089 /
  HELD-OUT 1,096 / all 2,185), i.e. `soto_replication.py` `score()` + `bipartite_score()` (scipy tie policy, register
  1045). Per cell × {DEV, HELD-OUT, ALL}: ARI, exact families, pair P / R / F1, `P_guard` and `cross_pairs` (halves),
  bipartite MICRO P / R / F (F = 2PR/(P+R)), MACRO P / R, undetected, predicted families. ALL additionally reports
  **cover-aware exact**: S1C families (all members, cover genes included, 491) equal as sets to a predicted family
  (PAIR's families are covers; partition rules give each gene one family) — the only metric in which PAIR's cover can
  count; not used for selection.
- PAIR is scored as a partition: a non-coding gene in ≥ 2 predicted families goes to the family it has the most kept
  pairs with; tie → the family whose smallest member gene id is smallest.
- A cell reproduces BASE exactly (REGION, GREEDY, median / mean) — asserted against `soto_v1_literal` partitions —
  before anything else runs; DP cells reproduce `soto_losses` C.
- **Selection is on DEV only**: the cell with the highest DEV ARI; tie → higher DEV exact → fewer levers off BASE.
  A cell whose DEV `cross_pairs` exceed 1% of its DEV predicted pairs is ineligible. **Primary = best of the 6
  Soto-only cells ("SOTO-CHOSEN")**; secondary = best of all 35 ("ALL-CHOSEN"). Both are then reported on HELD-OUT
  and on ALL (the latter labelled DEV + HELD-OUT, not a verdict).
- **Random-lever null:** 1,000 cells drawn uniformly with replacement (seed 0) from the 35 (and, separately, from the
  6); distribution of HELD-OUT ARI; report its median, 5-95% range and P(null ≥ chosen) (ties count against the chosen
  cell). Also the Spearman correlation of DEV vs HELD-OUT ARI over the 35 cells (does DEV rank transfer?).

## 6. Clauses (HELD-OUT; BASE held-out median-MAD ARI 0.8183, exact 131)

- **R1 (the Soto-aligned levers explain part of the gap):** SOTO-CHOSEN held-out ARI − 0.8183 ≥ +0.010 AND
  P(35-cell null ≥ SOTO-CHOSEN) ≤ 0.25. Pass → "RECONCILED IN PART by <levers>"; Δ passes but null fails → "NOT
  DISTINGUISHABLE FROM A RANDOM LEVER CHOICE"; Δ fails → "the Soto-aligned levers do not close the gap".
- **R2 (no-harm, reported with R1):** SOTO-CHOSEN held-out exact ≥ 131 − 3 and `P_guard` ≥ BASE − 0.020; failing R2
  qualifies R1 ("trades whole families / precision for ARI").
- **R3 (same conclusion):** number of flagship families (§7) exact on ALL under SOTO-CHOSEN ≥ under BASE.
- The same three are reported for ALL-CHOSEN, labelled OURS-INCLUDED when it uses an OURS level.
- **Guard F-guard:** if the chosen cell's HELD-OUT `cross_pairs` exceed 1% of its predicted pairs, its held-out
  numbers are reported unjudged.

## 7. Ablation and "same conclusion"

- **Ablation** (on HELD-OUT and ALL): from BASE, switch ONE searched lever to the chosen cell's value; from the chosen
  cell, revert ONE lever to BASE. Plus the two unsearched levers as single steps from BASE and from SOTO-CHOSEN:
  K = rule@0.3 (universe S1C ∪ curated eligible, as `soto_cur_rule`; REGION edges recomputed from the frozen
  `soto_v1_literal` PAFs over those genes' exons, EXON / CIGAR filtered from their supersets), and U = 1,793 (no
  islands, no attach; PAIR: coding nodes only). Plus the minimap2 2.17 sensitivity if run.
- **Flagship families** (fixed now from the paper, from S1C gene-name prefixes GPR89, NPY4R, PTPN20, PDZK1, HYDIN, FRMPD2, FAM72, SRGAP2,
  ARHGAP11, CD8B, DUSP22, ROCK1, CFC1, PRSS40, GPRIN2, LRRC37A, LRRC37B, DUX4, USP17L, FAM90A, TBC1D3, NBPF, NPIP,
  NOTCH2NL: the 38 S1C families holding such a gene): the 9 zebrafish-modelled pHSD families (l.191: GPR89 ID_369, NPY4R ID_405, PTPN20 ID_439, PDZK1 ID_425,
  HYDIN ID_214, FRMPD2 ID_366, FAM72 ID_354, SRGAP2 ID_462, ARHGAP11 ID_145), the other pHSD / selection families
  (CD8B ID_312, DUSP22 ID_346, ROCK1 ID_444, CFC1 ID_317, PRSS40 ID_438, GPRIN2 ID_370, LRRC37A/B ID_14, ID_15), the
  three > 50-paralog families named in the Results (l.83: DUX4, USP17L, FAM90A → ID_330, 331, 347, 348, 106, 107, 108,
  356) and the named core-duplicon / conversion families (TBC1D3 ID_341, 468, 469; NBPF ID_396, 397, 400, 401; NPIP ID_149, 151-155; NOTCH2NL ID_400). Status per family on ALL: exact / split (clean members in
  ≥ 2 predicted families or partly unplaced) / merged (matched predicted family holds foreign clean genes) /
  split+merged / undetected, under BASE, SOTO-CHOSEN, ALL-CHOSEN.

## 8. Attribution of the residual (computed on ALL for SOTO-CHOSEN; also for BASE)

G_union = union of all five edge sources (node universe 2,334). For each truth family T (all S1C members, cover genes
included), T's **reachable pieces** = connected components of T's induced subgraph in G_union, restricted to its clean
members.
- **ARI ladder.** P_c = truth with every family cut into its reachable pieces → c = 1 − ARI(P_c) (**unreachable**: no
  ≥ 98% exon link among the family's own members in any source). P_b = P_c with every piece further cut by the chosen
  rule applied to that piece alone (GREEDY / DP: the chosen split with the chosen statistic on the piece's famCN;
  PAIR: components of the piece's kept pairs, coding-only propagation) → b = ARI(P_c) − ARI(P_b) (**their rule wall**).
  d = ARI(P_b) − y (**other**: over-merges, edges one source lacks but another has, split positions, attach choices).
  **a (cover)** = 0 in ARI by the scorer's construction (the 149 cover genes are excluded); its size is reported in the
  cover-inclusive currency (any partition ≤ 409/491 exact; PAIR's cover-aware exact) and causally (re-run the chosen
  cell with the 149 cover genes deleted from its edges: Δ ARI, Δ exact).
- **Exact-family ledger** (491 − exact = c + b + a + d), first applicable cause per non-exact family: c = its clean
  members span ≥ 2 reachable pieces; b = its reachable piece is cut by the chosen rule applied to the piece alone, or
  S1C marks it "Manual merge"; a = it becomes exact when the 149 cover genes are removed from the chosen cell's edges;
  d = the rest. c is also split into acrocentric (chr13/14/15/21/22) vs other; the v2.0 acrocentric assembly wall is
  expected to be 0 here (native v1.0, no liftover) and is reported as such. Known-wall sizes reported beside it:
  12 S1C families with Family MAD ≥ 1, 2 "Manual merge", ~50 families with a member no ≥ 98% alignment reaches
  (`soto_v1_diff.md` §2(c)), cover cap 409.

## 9. Predictions (this author, before any new number)

- EXON builds more edges than REGION within the 2,334 genes: 0.60.
- SOTO-CHOSEN uses F = PAIR: 0.65; uses E = EXON: 0.50.
- R1 passes: 0.50; R2 passes given R1: 0.60; R3 passes: 0.65.
- SOTO-CHOSEN ALL ARI ≥ 0.75: 0.50; ≥ 0.80: 0.30.
- PAIR cover-aware exact > 409 in its best cell: 0.15.
- ALL-CHOSEN uses an OURS level: 0.55.
- c (unreachable) ≥ 0.05 ARI: 0.45.

## 10. Hostile self-review (applied above)

1. **Selection bias.** 35 cells, one picked on DEV; the null (random cell) and the DEV→HELD-OUT rank correlation are
   the controls. With so few cells the null is coarse; P ≤ 0.25 means "top quartile", not significance.
2. **The held-out half is partly spent** (§2). Known held-out values exist for 12 of the 35 cells; the selection rule
   cannot see them (DEV only) but I have seen them. The decisive new cells (EXON, PAIR) are blind.
3. **"Aligned with Soto" is not "better".** A Soto-code lever can lower ARI (e.g. PAIR's cover collapsed to a
   partition). The report says so per lever; the attribution is about direction and size, not wins.
4. **Circularity.** famCN, the universe and the 71 curated genes are Soto's; PAIR uses S1C famCN at the pair level —
   the same famCN Soto used to make S1C. Agreement it produces is concordance by construction, stated as such.
5. **PAIR's cover collapse** (majority tie rule) is our scoring necessity, not Soto's; the cover-aware exact metric
   is reported so the collapse does not hide what the cover achieves.
6. **Per-gene famCN.** Soto's `get_mad` may see several WSSD rows per gene (per SD98 fragment); S1C gives one median
   per gene. The pair gate here uses S1C's value — an approximation, disclosed.
7. **Metric traps checked** (`feedback_metric_traps`): fixed scored universe (no denominator conditioned on the
   prediction); bipartite tie policy named; no hand-picked families (flagship set fixed here from the paper, and the
   headline is all 491); an edge-count-matched null is not used as proof.
8. **One statistic per cell.** "ARI median-MAD and mean-MAD" are the S columns of each (E, F) row; PAIR rows repeat.

## 11. Order, stop rules, machine rules

(1) freeze this file (sha1) → (2) build EXON edges (heavy lock), write the harness, assert BASE / C reproduction and
PAIR on hand-made toy inputs, record counts only → (3) Amendment 1 (instrument sha1s) → (4) one run: DEV scores and
the two selections are written first, then HELD-OUT and ALL, nulls, ablations, flagship, attribution → (5) Outcome.
No lever, level, clause or bar changes after (3); an instrument fix that changes no rule is an Amendment.
Light lock for scoring, heavy lock for minimap2; TMPDIR under `/mnt/linuxdisk`; never `pkill -f`.

## Amendments

### Amendment 1 — the freeze (2026-09-29 10:39, written BEFORE any cell was scored)

**This file's sha1 before this amendment:** `c0f7a4bde904df3f7c0e1dc30a305eb67316868f` (frozen 10:26; byte copy
`soto_reconcile/frozen/PREREG_soto_reconciliation_2026-09-29.pre_amendment1.md`). **Instruments:**
`soto_reconcile/frozen/SHA1SUMS` sha1 `7d719558`: `recon_lib.py` 35c5050e, `run_grid.py` 4ed616a5, `checks.py`
44ec43d2, `build_exons.py` 23aeed57, `build_exon_pairs.py` fda610c8, `gen_exons.py` a674dc05, `nb_exact.py` bf94fc8b;
EXON data `sd98_exons.bed` 4bccb1b1, `sd98_exons.fa` c3b34253, PAF v2.30 a83cb26c / v2.17 17b11561, pair files
d2d36db0 / 039cb8a0; REGION superset for the K ablation `region_k/r2405.project.tsv` 3714d63c; Soto code A 736cb4c2,
B 9365476b; plus the unchanged soto_losses `lib.py` ff0d105a, `split.tsv` 49bcbcfe, `base_edges.tsv` 6cc14eeb,
`cigar_rows.tsv` 2d7175a9, `soto_replication.py` 13ffb698, S1C d008a179.

**What was built and checked (no agreement with S1C was computed for any new cell):**
- EXON (§4): 25,971 SD98 exons / 5,154 genes / 8.99 Mb (their line-172 gene count). Index `minimap2 -d` 86 s /
  12.2 GB; map-back 16 s (v2.30: 215,716 alignments, 24,448 exons mapped); v2.17 release binary (`minimap2-2.17_x64-
  linux`, on-the-fly index) 71 s, 224,016 alignments. `-f 0.99 -s` hits 192,369, same-gene 43,964, both non-coding
  42,705 → 12,231 pairs over 2,351 of the 5,154 genes (v2.17: 12,462 / 2,370). Within the 2,334-gene universe:
  EXON 11,522 edges (v2.17 11,684; 11,483 shared), REGION 11,935, CIGAR 4,384, UNION 12,001, B1 11,973; K universe
  (2,347): 11,628 / 12,059 / 4,443 / 12,125 / 12,097. **Prediction "EXON > REGION in edge count" is therefore already
  false.** B1 for the K ablation uses the 2,334-universe self-overlap CIGAR edges (approximation, as stated in §7's
  spirit; affects only K ablations of B1 cells).
- REGION rebuilt from the frozen `soto_v1_literal` PAFs with `literal_edges.py` over the 2,405-gene superset (S1C ∪
  the 1,864 recomputed eligible genes): its restriction to the 2,334 genes is byte-identical to `base_edges.tsv`
  (11,935, same order); the 2,334-gene rebuild alone also reproduces it byte-for-byte.
- `checks.py`: BASE and C (REGION × GREEDY / DP × median / mean) partitions reproduced by `recon_lib`; PAIR on a toy
  (chaining past a family MAD of 1, a cover gene in two families, a non-coding–non-coding pair dropped, tie rule) as
  specified; PAIR equals an order-preserving port of notebook cells 8-11 on EXON (504 families, 2,210 genes, 158
  genes in ≥ 2 families) and REGION (506 / 2,104 / 151). These counts were seen before scoring (S1C: 491 families,
  149 cover genes — a structural count, not an agreement score).
- **Instrument finding (no rule changed):** the released notebook's closure loop reassigns
  `gene_cluster = list(set(gene_cluster + cluster))` while scanning, so which members get expanded depends on the
  string-hash order. An exact port (`nb_exact.py`, the element at `gene_cluster[i]` re-read after every merge) on EXON
  pairs gives 651-654 families with 625-707 coding / unprocessed genes in more than one family under PYTHONHASHSEED
  0/1/2 (REGION, seed 0: 686 / 642). S1C has **0** coding / unprocessed genes in more than one family. So S1C is not
  the unmodified output of that loop on raw pairs (either their pairs closed fully by chance, or the table was
  reconciled afterwards — "Manual merge" appears twice). PAIR keeps the pre-registered meaning (the full closure the
  loop intends); the exact loop is not added as a level (non-deterministic, contradicts S1C).
- `run_grid.py` implements §5-§8 as written: cells, DEV-first selection (written to `dev_selection.json` before any
  held-out score is computed), null (1,000 draws, seed 0), R1-R3, F-guard, ablations (searched levers; K; U;
  minimap2 2.17 on the EXON cells of the Soto-only set and on the chosen cells), flagship statuses, attribution. The
  rule wall b is additionally split into "admission" (the piece is already cut with the famCN gate switched off:
  eligible-only backbone / islands / leaves) vs "famCN" (cut only by the MAD rule), and "Manual merge" families are
  their own b sub-class — a finer report of the pre-registered b, not a change to it.

**Command:** `bash tools/rlock.sh light python3 frozen/run_grid.py runs/main` (one run).

### Amendment 2 — the run (2026-09-29 10:40; no rule, level, clause or bar changed)

One run, 12.8 s; `sha1sum -c frozen/SHA1SUMS` 28/28 OK (first attempted from the wrong directory, re-run from
`frozen/` immediately after: unchanged). BASE reproduced (asserted). DEV selection written before held-out scoring
(`runs/main/dev_selection.json`). No cell was ineligible (DEV `cross_pairs` ≤ 1%). Cross-check: the SOTO-CHOSEN and
BASE partitions were written out and re-scored with `soto_replication.py score` (identical: 0.9698 / 479, 0.7096 / 235),
and an independent pair count gives 10,937 true of 10,942 predicted clean pairs for SOTO-CHOSEN.

## Outcome (2026-09-29)

**Verdict: RECONCILED — by two of Soto's own code-level choices.** SOTO-CHOSEN = ALL-CHOSEN = **EXON × PAIR** (both
SOTO-CODE: exon-level map-back, per-pair famCN gate with coding-gene propagation). R1 PASS (held-out ARI 0.8183 →
**0.9681**, Δ +0.1498; rank 1 of 35; P(35-cell null ≥) 0.034; 6-cell null 0.169 = 1 of 6), R2 PASS (held-out exact
131 → 263 of 266, `P_guard` 0.845 → 1.000), R3 PASS (flagship exact 13 → 36 of 38). F-guard: 0 cross-half pairs.
DEV → HELD-OUT rank transfer: Spearman 0.71 (35 cells), 0.77 (6 cells). No OURS level is needed or selected.

| | ARI (median / mean MAD) | exact | pair P / R / F1 | MICRO P / R / F | MACRO P / R | undetected | cover-aware exact |
|---|---|---|---|---|---|---|---|
| BASE (REGION × GREEDY), DEV / HO / ALL | 0.5957 / 0.8183 / **0.7096** (mean ALL 0.7039) | 104 / 131 / 235 | .831 / .621 / .711 | .784 / .718 / .749 | .699 / .704 | 113 | 216 |
| **EXON × PAIR**, DEV / HO / ALL | 0.9708 / **0.9681** / **0.9698** (statistic-invariant) | 216 / 263 / **479** | **1.000** / .942 / .970 | **1.000** / .980 / .990 | .999 / .992 | **0** | **481** |

Every cell: `soto_reconcile/runs/main/cells.tsv` (35 cells × DEV / HELD-OUT / ALL). Best cell with any OURS or INPUT
level: held-out 0.8966 is EXON × GREEDY-median (Soto-only); the best non-Soto cell is EXON × DP-median 0.8890 —
B1, UNION, CIGAR, DP and "both" add nothing once E and F are Soto's.

**Ablation (ALL; held-out in brackets).** From BASE: E = EXON alone +0.1263 ARI / +103 exact (+0.0783 / +62); F = PAIR
alone +0.1680 / +78 (+0.0477 / +36); both +0.2601 / +244 (+0.1498 / +132) — sub-additive in ARI, super-additive in
exact families. From EXON × PAIR: E back to REGION −0.0921 / −166 (−0.1021 / −96); F back to GREEDY-median −0.1338 /
−141 (−0.0715 / −70). The MAD statistic is not a lever under PAIR (pair MAD median = mean). Unsearched levers: their
hand-curated list → our rule@0.3 costs −0.0207 / −23 exact at EXON × PAIR (vs +0.0043 / −2 at BASE); 1,793-only
universe −0.1394 / −193 (BASE −0.0608 / −64); minimap2 2.17 (theirs) instead of 2.30 −0.0013 / −9.

**Same conclusion (flagship, 38 families).** BASE exact 13, split 5, merged 12, split+merged 2, undetected 6 →
EXON × PAIR exact 36, split 1 (ID_347, DUX4 "Manual merge"), merged 1 (ID_401 NBPF4/6 + NBPF5P, which S1C leaves
unassigned). All 9 zebrafish pHSD families, CD8B, ROCK1, CFC1, PRSS40, GPRIN2, LRRC37A/B, TBC1D3, NOTCH2NL, all six
NPIP families (ID_149, 151-155), FAM90A (ID_356, 56), USP17L (ID_106, 56) and DUX4L (ID_330, 60) are exact.

**Attribution (ALL).** "BASE 0.7096 → EXON × PAIR 0.9698 (held-out 0.9681); remaining gap 0.0302 = a 0 (cover; the
scorer excludes the 149 cover genes, and EXON × PAIR reproduces 481/491 families exactly INCLUDING cover genes, above
the 409 partition cap) + b 0.0231 (their rule as reconstructed: 7 families the per-pair gate on S1C's one-value-per-
gene famCN cannot join — ID_28, 63, 99, 113, 144, 163, 270 — plus the manual merge ID_347) + c 0.0002 (ID_482 UBTFL,
also a "Manual merge"; no ≥ 98% exon link in any source) + d 0.0070 (ID_62 / ID_192: a non-coding gene in several
predicted families collapsed to the wrong one for scoring; ID_401: joins a gene S1C leaves unassigned)." Exact-family
ledger: 491 − 479 = 12 = a 1 + b 8 + c 1 + d 2. BASE for comparison: gap 0.2904 = a 0 + b 0.1287 + c 0.0002 +
d 0.1615; ledger 256 = 1 + 12 + 1 + 242. Causal cover check at EXON × PAIR: deleting the 149 cover genes' edges
changes ARI by +0.0007 but loses 46 exact families (single-clean-gene families held together by their non-coding
members). Acrocentric wall: 0 (native v1.0).

**The four "known walls", re-measured.** (1) Cover cap 409: a property of partition methods only; Soto's own rule
makes the cover — EXON × PAIR puts 158 genes in ≥ 2 families, including all 149 of S1C's. (2) "12 families violate
their own MAD < 1": not violations — their code never gates a family's MAD; 9 of 12 are exact under EXON × PAIR (ID_28,
63, 113 remain, in b). (3) "~50 families with members no ≥ 98% alignment reaches" (`soto_v1_diff.md`): an artifact of
the region reading — with REGION edges alone 73 families are cut at the alignment stage, with EXON edges alone 2 (ID_347,
ID_482, exactly S1C's two "Manual merge" families), with the union 1. (4) Acrocentric: 0.

**Predictions vs outcome.** EXON > REGION edges (0.60): no (11,522 < 11,935, known at the freeze). SOTO-CHOSEN uses
PAIR (0.65): yes; EXON (0.50): yes. R1 (0.50): yes; R2 | R1 (0.60): yes; R3 (0.65): yes. ALL ARI ≥ 0.75 (0.50) and
≥ 0.80 (0.30): yes. PAIR cover-aware exact > 409 (0.15): yes (481). ALL-CHOSEN uses an OURS level (0.55): no.
c ≥ 0.05 (0.45): no (0.0002). I under-predicted the size of the effect badly.

**Post hoc (after the numbers; not clauses; `soto_reconcile/posthoc/`).** (i) Using the full SD98 node universe (all
5,083 genes minus their 71) instead of S1C's 2,334: identical ARI / exact (cover-aware exact 480); only 19 extra
non-coding genes (15 lncRNA, 4 miRNA) would receive a family — restricting candidates to S1C's gene list does no work.
(ii) Headline counts under EXON × PAIR vs the paper: 504 vs 491 families; 1,672 vs 1,679 paralogs in families (their
notebook prints 1,673); 121 vs 114 singletons; 278 vs 271 families of 2-3 members; the same three families > 50
members (USP17L 56, FAM90A 56, DUX4L 60).

**What this changes.** The Soto replication headline can move from "ARI 0.71, 48% exact" to **"ARI 0.97, 479 of 491
exact (held-out 0.968), with Soto's own released code choices; the difference from the old headline is two named
choices: map SD98 exons (not regions) back to the genome, and apply MAD < 1 to each shared-exon pair then grow
families through coding genes (not to connected components)."** It remains concordance: famCN, the gene universe and
the 71 curated genes are Soto's (register 858 / 1085), and PAIR uses S1C famCN pair by pair. The statements in
`SOTO_REPLICATION_STATUS_2026-09-28.md` §1.1 ("cannot match by construction": cover 409 cap, 12 MAD ≥ 1 families, 50
unreachable families) and the §6if "their real algorithm scores worse" (0.54-0.56; it read clusters as components) are
superseded. Nothing adopted into `bench/` or `src/`; nothing committed; the user decides whether to promote
EXON × PAIR into `soto_replication.py` as the headline arm.

**Register rows (appended 2026-09-29 as 1162-1166; the `register draft 1152` cited in §3 is row 1158).**

| # | date | area | claim | verdict |
|---|---|---|---|---|
| 1162 | 2026-09-29 | Soto replication (reconciliation, held-out by family hash) | Setting the edge source and family rule the way Soto's released code sets them (map SD98 exons back with `-c --end-bonus 5 --eqx -N 50 -p 0.5`; `-f 0.99 -s` exon cover, self = same gene; MAD < 1 per shared-exon pair; families grown through coding / unprocessed genes, non-coding genes joining every family they pair with) reproduces Table S1C | ✅ **RECONCILED** (prereg `PREREG_soto_reconciliation_2026-09-29.md`, c0f7a4bd). DEV-selected among 35 cells; held-out ARI 0.8183 → 0.9681 (rank 1/35, null P 0.034), exact 131 → 263/266; all 491: ARI 0.9698, 479 exact, pair P/R 1.000/.942, MICRO P/R 1.000/.980, 0 undetected, cover-aware exact 481 (> 409 partition cap); flagship 13 → 36/38. Concordance: S1C famCN / universe / curation list are Soto's |
| 1163 | 2026-09-29 | Soto replication (attribution) | The gap between our literal recipe (0.7096) and S1C is the two prose-vs-code choices, not the MAD statistic, SEDEF's own CIGARs, B1, a threshold-free split, curation or version | ✅ from BASE: exon queries +0.126 ARI / +103 exact, pair gate +0.168 / +78, both +0.260 / +244; MAD statistic invariant under the pair gate; B1 / UNION / CIGAR / DP / median∧mean add nothing on top (best such cell held-out 0.889); minimap2 2.17 (theirs) −0.001 / −9 exact vs 2.30; our curation rule instead of their list −0.021 / −23 exact at the reconciled cell |
| 1164 | 2026-09-29 | Soto replication (walls) | The "cannot match by construction" walls of the 09-28 status (cover cap 409, 12 families with Family MAD ≥ 1, ~50 families no ≥ 98% alignment reaches, acrocentric assembly) bound any replication | ⛔ **No — three of four were artefacts of our readings.** Their pair-level rule MAKES the cover (158 multi-family genes incl. all 149 of S1C's) and never gates family MAD (9/12 exact); "unreachable" = region reading (73 families cut by REGION edges, 2 by EXON edges = the two S1C "Manual merge" families ID_347/ID_482); acrocentric 0 on native v1.0. Residual 12 families: 7 per-pair famCN (one S1C value per gene), 2 manual merges, 3 cover-collapse / singleton |
| 1165 | 2026-09-29 | Soto code audit | Soto's released notebook `B_SD98_families.ipynb` computes families as the closure of low-dispersion clusters through coding genes | ⚠ **Its loop is order-dependent**: `gene_cluster = list(set(...))` inside the scan re-orders the list, so an exact port leaves 625-707 coding genes in > 1 family under PYTHONHASHSEED 0/1/2 on our exon pairs; S1C has 0. S1C is not the unmodified output of that loop (or their run happened to close); the intended full closure matches S1C (above). The 09-11 port (clusters read as components) was a misreading: 0.4966 |
| 1166 | 2026-09-29 | Soto code audit | The paper's STAR Methods describe Soto's family construction as run | ⚠ **Two prose-vs-code differences carry the whole replication gap**: prose "DNA sequences of all SD98 regions ... mapped back" vs code maps SD98 exons (`A_SD98_regions.md` l.226-238); prose "groupings where the mean absolute deviation ... less than one" vs code gates each shared-exon pair with a median absolute deviation (`B_SD98_families.ipynb` cells 6-11). S1C also carries 2 manual merges ("Family MAD" = "Manual merge": ID_347, ID_482) the prose does not mention |

## Independent verification (2026-09-29, own code and own exon re-mapping; `figs/soto_reconcile_verify.md`,
## scripts `/mnt/linuxdisk/tmp/rustle_figures_dev/soto_reconcile_verify/`)

**CONFIRMED WITH CORRECTIONS.** Every headline number reproduces exactly (ARI 0.9698, 479/491 exact; held-out 0.9681,
263/266; 5 false pairs / 10,942 = precision 0.9995; cover-aware 481; flagships 36/38; 12 residuals). Split, prereg
freeze (sha1 c0f7a4bd, 10:26, body unchanged) and the 35 configurations check out (30 non-PAIR cells not recomputed;
none within 0.18 on DEV). Their code does both steps as described (A_SD98_regions.md l.226-238, 248-263; notebook cells
7-8); the notebook's printed 1,673 coding genes matches the pair reading (1,672) not components (1,206).
Corrections:
- **Released loop as-is scores ARI 0.815-0.872 (377-398 exact) over 8 hash seeds**: it stops early, wrongly expands
  non-coding genes and merges 87-116 families. "Complete grouping" = the loop's evident intent and S1C's structure,
  but it is a REPAIR worth +0.10-0.15 ARI; say so.
- **Attribution of the gain:** sequence only (no copy-number test) 0.73; S1C copy numbers shuffled within connected
  groups ~0.76; **our own copy numbers (not S1C) 0.92 (held-out 0.93, 373 exact)**; S1C copy numbers 0.97. The last
  ~0.05 ARI / ~106 exact families are agreement with S1C's published values.
- **Universe:** scoring only S1C genes hides errors (all 541 non-coding S1C genes are members by construction).
  Admitting every duplicated-region gene and scoring the 19 extra it places: ARI 0.959, false pairs 5 -> 260.
- **Curation:** the 71 removals come from S1C; adding them back as ordinary genes gives 0.94 / 425 exact.
- **Prose vs code:** only step 1 (exon vs region map-back) is a real text-vs-code difference; the Methods text is
  ambiguous on step 2 and fits the pair reading.
Advisor-safe sentence: "Using Soto et al.'s released code choices (exon map-back, a per-gene-pair copy-number test, and
their family loop completed as it evidently intends) on their genome, annotation, gene list and published copy numbers
reproduces 479 of 491 families (ARI 0.97, 0.97 held-out); sequence alone reaches 0.73 and our own copy numbers 0.92,
so the step to 0.97 is agreement with their published copy numbers, not an independent reconstruction."

**Correction (2026-09-29, later the same day; docs/archive/2026-09/PREREG_soto_famcn_allwssd_2026-09-29.md):** the "our own copy
numbers 0.92" in the advisor sentence above came from a 10-sample WSSD table that turned out to be a favourable draw
(20 random 10-sample draws: ARI 0.870-0.923). With all 268 SGDP samples and Soto's per-gene interval (gene body ∩
SD98) our own copy numbers give **ARI 0.9277 (held-out 0.9343), 411/491 exact**. Use that figure instead of 0.92/373.
