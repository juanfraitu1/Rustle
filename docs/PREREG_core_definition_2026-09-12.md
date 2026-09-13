# PREREG — which contiguous-core definition should `confirm_edge` use? (2026-09-12, written before any score below is computed)

## Question
On real gorilla candidate pairs, does the global poasta contiguous core (production) separate same-family
from different-family rep pairs better than an exact local core (LCS), with or without repeat masking?

## Pairs
Every candidate pair (`candidate_pairs`, production params) in both cells of the thread-free dump
`o1_falsemerge/lcs/{baseline,proposal1}/candidates.tsv` (192 + 207; proposal-#1 pairs whose two sequences are
identical to a baseline pair reuse that pair's values). Pairs with max length > LEN_CAP (20 kb) have no poasta
value in production (it uses LCS there); they are scored for the LCS arms only and reported separately.

## Scores (T_CORE = 0.13 for the edge call; each also scored threshold-free)
- **POA**: production `confirm_edge` value (exact poasta, forward then reverse complement if forward < T_CORE).
- **LCS**: longest common substring / min length, same orientation rule.
- **LCS-masked**: same, after replacing every base covered by a gorilla RepeatMasker (Dfam curated) record with N.

## Labels (independent of all three scores)
- **SAME**: some GGO SEDEF pair has one side overlapping rep A's genomic footprint and, projecting A's
  footprint through that SD (offset within side 1 -> side 2, strand-aware, linear), lands within 2 kb of rep
  B's footprint (or vice versa).
- **DIFF**: no SEDEF pair links the two footprints at all.
- **AMBIG**: an SD links the two loci but the projection misses by > 2 kb (same duplicon, different position);
  excluded from the main score, reported separately.
Footprint = exon blocks for spliced reps, span for single-exon and genomic-span reps.

## Metrics
Edge-level precision, recall, F1 at T_CORE on SAME vs DIFF pairs; ROC AUC of the raw score. Per cell, never
pooled with any human number.

## Decision rule
- If **LCS-masked** has F1 >= POA F1 AND precision >= POA precision - 0.05: recommend replacing poasta in
  `confirm_edge` with LCS-masked (no budget, no leak). Fix the budget only for other poasta callers if any.
- Else if **LCS** (unmasked) meets the same bar: same recommendation, repeat-masking not needed.
- Else: keep poasta; implement the vendored step-limit fix for the leak.
- If SAME or DIFF has fewer than 10 pairs in a cell: that cell is underpowered; report, do not decide from it.

## Known limits (declared now)
- One bounded, NPIP-dense neighbourhood on one substrate; a recommendation here needs a held-out check
  (human is allowed per user 09-12) before any default change.
- SEDEF labels miss very old/divergent duplications (DIFF may contain some true paralogs), and the 2 kb
  projection tolerance is a declared choice.
- Candidate pairs are k-mer pre-filtered, so DIFF pairs are hard negatives by construction.

---
## ADDENDUM A (2026-09-12, written after the gorilla LABEL counts and 23/192 poasta values, before any POA-vs-LCS comparison was looked at)

**Why:** the gorilla NPIP neighbourhood is ~all positional paralogs — baseline labels SAME 189 / DIFF 2 /
SAME_LOCUS 1; proposal #1 DIFF 13 but all > LEN_CAP (no POA arm). Under the rule above both cells are
UNDERPOWERED for precision; they will be reported, not decided from. Implementation detail declared: "linked"
= an SD side overlaps one footprint and the other side lies within 2 kb of the other footprint; footprints that
overlap each other on the genome are labelled SAME_LOCUS and excluded.

**A1 — gorilla, recall only (exploratory):** on SAME pairs <= LEN_CAP, report the 2x2 of POA-pass vs LCS-pass,
stratified by the SD-projected overlap fraction (bp of the shorter footprint projecting onto the other
footprint / shorter footprint bp): [0, 0.13), [0.13, 0.5), [0.5, 1]. Question answered: how often does each
definition miss a core that DNA says is there.

**A2 — human, the deciding arm.**
- Reps: one spliced transcript per gene for the 2,334 Soto SD98 genes (LiftOff v2.0 models,
  `liftoff_v1_to_v2/lifted_v2.gff3`, original CHM13_G genes only, no LOFF extra copies), sequence from
  `chm13v2.0.fa`, reverse-complemented for minus-strand. Genes with an ambiguous Soto family assignment are
  dropped.
- Pairs: production `candidate_pairs` over these reps. Pairs with either gene lacking a Soto family are dropped.
- Labels: SAME = identical Soto family ID. DIFF = different Soto family IDs AND no human SEDEF
  (`final_human.bed`) pair links the two genes' footprints (same "linked" definition as above). Different
  families but SD-linked = AMBIG (reported, excluded). Soto families are an independent DNA-based published
  truth, not derived from any of the three scores.
- Sample: all pairs <= LEN_CAP if <= 600; else a seeded (seed 20260912) random sample of up to 300 SAME + 300
  DIFF. Pairs > LEN_CAP are scored for the LCS arms only.
- LCS-masked on human uses the genome's soft-mask (lowercase -> N), not a RepeatMasker .out (none on disk
  for v2.0); declared as a difference from the gorilla arm.
- Execution: poasta per pair in its own process (8 GB virtual cap, 180 s timeout), at most 4 processes at once.
  Pairs that time out are reported as a count and excluded from POA metrics.
- Decision rule: as above, applied to the human arm. The gorilla arm cannot overturn it (underpowered), but a
  gorilla A1 result where POA misses many DNA-supported cores that LCS finds is reported alongside.

---
## ADDENDUM B (2026-09-13, after §6ja's decision; before any step-2/3/4 number exists) — family-level, threshold, held-out

LCS mode is implemented (`DetectParams::edge_core`, `RUSTLE_EDGE_CORE=lcs`, commit a577d37e). Steps run in order.

**Step 2 — family level (human, development substrate).**
- Subset: connected components of the human candidate-pair graph (12,433 pairs, 241 components), excluding the
  four largest (4,813 / 1,467 / 1,076 / 817 pairs). Shuffle the rest with seed 20260913 and add whole components
  until >= 1,500 candidate pairs. Reps = every rep in the chosen components. Limitation declared: the largest SD
  families are not in the family-level test.
- POA edges: production `confirm_edge` (EdgeCore::Poa, T_CORE 0.13), one process per pair, <= 4 at a time, 8 GB
  cap, 120 s timeout; a timeout counts as NO edge (production would stall), count reported.
- LCS edges: `confirm_edge` EdgeCore::Lcs, T_CORE 0.13.
- Families: shipped `decompose_families(edges, SplitParams::default())` for each mode.
- Truth: Soto family IDs of the subset reps.
- Metrics: ARI over reps; co-membership pair precision / recall / F1; exact matches (a predicted family equal to a
  truth family restricted to the subset, truth families with >= 2 subset reps); number of families and
  singletons.
- Verdict: LCS "holds at family level" iff LCS ARI >= POA ARI AND LCS pair-F1 >= POA pair-F1. Otherwise report
  where edge gains fail to translate.

**Step 3 — LCS threshold (human development pairs: all labelled SAME/DIFF pairs <= LEN_CAP).**
- (a) fraction threshold t on LCS/min_len, F1-maximising over t = 0.01..0.60 step 0.01;
- (b) absolute LCS bp threshold, F1-maximising over 20..2,000 bp step 10;
- (c) theory, Arratia-Waterman null for the longest exact match between unrelated sequences of lengths m, n:
  admit iff LCS >= k*(m,n) = ln((1-p) m n / alpha) / ln(1/p), with p = sum over bases of q_b^2 from the base
  composition of all human reps and alpha = 1 / 12,433 (one expected chance match across the candidate set);
- plus the incumbent t = 0.13.
- Selection: choose (c) if its F1 >= max(F1(a), F1(b)) - 0.02 (a derivable rule is preferred when it costs little);
  else the better of (a), (b). Report all four. Re-run step-2 LCS families at the chosen rule (report only).

**Step 4 — held-out (gorilla, genome-wide, never used for development).**
- Reps: all nodes in `o1_reps/dump/ggo.er._k11_w5.0.reps.fa` (Aug-21 run), EXCLUDING any rep overlapping the 62
  NPIP-seeded SD windows (`o1_fromgenome_sd/npip_seeded_windows2.bed`, the development region). Footprint = span.
- Pairs: production `candidate_pairs` (`fm_pairs`).
- Labels: GGO SEDEF projection exactly as the gorilla arm (SAME / DIFF / AMBIG / SAME_LOCUS, 2 kb tolerance).
- Sample: seeded (20260913) 300 SAME + 300 DIFF <= LEN_CAP; if a class has < 300, take all of it.
- Scores: POA @ 0.13 (one process per pair, <= 4 at a time, 8 GB cap, 180 s timeout; timeouts reported and
  excluded from POA; LCS arms reported with and without them), LCS @ 0.13, LCS @ the step-3 rule.
- Recommendation for a DEFAULT flip (still the user's decision): LCS @ step-3 rule F1 >= POA F1 AND precision
  >= POA precision - 0.05 on the held-out sample. If not met: no flip, report why.
- Human and gorilla numbers are never pooled.

---
## ADDENDUM C (2026-09-13, after §6jb; before any number below exists) — confirmation of LCS @ 0.13 on unseen held-out pairs

**Why:** §6jb's flip criterion targeted the step-3 rule (LCS >= 50 bp), which failed; LCS @ 0.13 met the same bar only
as a report-only row. This is the pre-registered confirmation for LCS @ 0.13 itself.
- Pool: `step4/labeled_pairs.tsv` (gorilla genome-wide held-out, NPIP development windows excluded), SAME and DIFF
  pairs <= LEN_CAP, EXCLUDING every pair in the step-4 sample (`step4/bridge/sample_pairs.list`).
- Sample: seed 20260914, 300 SAME + 300 DIFF (all of a class if fewer).
- Arms: POA @ 0.13 (production `confirm_edge`, one process per pair, <= 4 at a time, 8 GB cap, 180 s timeout;
  timeouts reported and excluded from POA; LCS reported with and without them) vs LCS @ 0.13.
- **Criterion (recommend the default flip to LCS @ 0.13):** on poasta-finished pairs, LCS F1 >= POA F1 AND LCS
  precision >= POA precision - 0.05. Otherwise no flip.
- Limitations declared: edge level only (the family-level check was on human, §6jb step 2); same substrate and
  label source as step 4, different pairs.

---
## ADDENDUM D (2026-09-13, after §6jd; before any number below exists) — does the LCS core narrow the DE NOVO <-> GUIDED gap?

User goal (09-13): two O1 modes — DE NOVO (RNA IsoSeq reps) and GUIDED (annotated GFF) — and the aim is to reduce
the difference between them. §6jd found the default de novo catalog builds E_r edges from minimap2
(`homology_blocks`), not `confirm_edge`, so the LCS default has not touched it.

**Catalogs (existing, provenance disclosed; not re-run with today's code — a de novo rebuild is 2 h 19 min and
background jobs are being killed by the harness):**
- GUIDED: `mcl_ann/gw_units_v3.clusters.tsv` (2026-09-06; mcl_families, GGO_genomic.gff, core_refine + SEDEF, 13,263
  nodes, 2,296 clusters). One row per member locus (cluster_id, chrom, start, end).
- DE NOVO: gorilla homology catalog run of 2026-08-21 (`o1_reps`): node set = the 17,924 reps in
  `dump/ggo.nodes.tsv` + `ggo.er._k11_w5.0.reps.fa`; shipped E_r edges = `dump/ggo.edges.tsv` (4,778).
  The de novo NODE SET is held fixed across all variants; only the edge definition changes.

**De novo variants** (families = shipped `decompose_families(edges, SplitParams::default())`, `fm_decompose`):
- V0 shipped E_r edges (weight = min(1, identity x coverage)). Sanity gate: report V0 vs shipped
  `ggo_reps.families.tsv` co-membership agreement (the shipped catalog adds gates, so exact equality is not required).
- V1 filter: V0 edges that also pass LCS >= 0.13 x min length (forward, else reverse complement).
- V2 replace: production `candidate_pairs` + LCS >= 0.13 (= today's default `detect_edges`) over all 17,924 nodes.
- V3 union: V0 edges plus V2 edges.

**Mapping.** A guided locus is covered by a de novo rep if their spans overlap on the same contig (>= 1 bp). A
locus's de novo family set = families of all overlapping reps that sit in a multi-member family.

**Gap metrics** (guided universe = loci listed in the guided clusters file):
- R_G (fixed denominator): of all locus pairs co-clustered in GUIDED, the fraction whose de novo family sets
  intersect. The denominator depends only on the guided catalog, so it is identical across variants.
- P_G: of all locus pairs (both in the guided universe) whose de novo family sets intersect, the fraction
  co-clustered in GUIDED.
- ARI over guided-universe loci that have >= 1 familied rep (de novo label = family of the max-overlap rep).
- Coverage (reported): loci with any overlapping rep (fixed); loci with a familied rep (per variant).

**Decision.** A variant NARROWS the gap iff R_G > R_G(V0) AND P_G >= P_G(V0) - 0.05. If several do, the
recommendation is the one with the highest R_G. If none do, report that the LCS core does not narrow the gap on
this node set and that node construction (§6j1) remains the dominant difference.

**Limits declared:** catalogs predate today's code (dates above); span-overlap mapping ignores exon structure;
de novo co-familying of a clustered locus with an unclustered one is not counted in P_G; gorilla only.

---
## ADDENDUM E (2026-09-13, after §6je and commit 50736042; before any number below exists) — rebuilt de novo catalog with today's code

**Why:** §6je used an Aug-21 de novo catalog and `decompose_families`; production partitions with
`gamma_quasi_clique_partition` + coverage split + distinct-locus gate. `RUSTLE_ER_UNION_LCS=1` now exists in code.
A genome-wide rebuild (2 h 19 min) does not fit the tool limits, so this check uses the 3-contig substrate
(NC_073241.2, NC_073242.2, NC_073244.2) that fits a foreground call (the same substrate's catalog built in 6 min 22 s).
- De novo catalogs, today's code (commit 50736042): `gw_family_catalog --bam <GGO_ds.bam restricted to the 3
  contigs> --fasta GGO.fasta --homology-primary --threads 4`, run twice: DEFAULT and `RUSTLE_ER_UNION_LCS=1`.
  Families = the emitted `copies.tsv` (family_id, chrom, start, end) — the real production output.
- Guided: `mcl_ann/gw_units_v3.clusters.tsv` restricted to loci on the 3 contigs; clusters with >= 2 such loci.
- Mapping and metrics exactly as Addendum D (span overlap; R_G fixed-denominator guided-pair recall; P_G; ARI;
  node presence), plus the report-only reachable-pair recall and per-cluster mean.
- **Decision:** the union narrows the gap on rebuilt catalogs iff R_G(union) > R_G(default) AND P_G(union) >=
  P_G(default) - 0.05. If met, it is still NOT a default change: the genome-wide rebuild and a second substrate come
  next (hold-a-substrate-back). If not met, report and stop pursuing the union.
- **Limits declared:** 3 contigs, NPIP-dense (part of the core-definition development region); guided catalog is
  2026-09-06 and genome-wide-clustered (cross-contig partners dropped).

---
## ADDENDUM F (2026-09-13, after §6jf; before any number below exists) — second substrate: human

**Why:** hold-a-substrate-back. §6jf's rebuilt catalogs are gorilla, NPIP-dense.
- Reads: `winloci_data/A119b_ds.bam` (human IsoSeq, CHM13 v2.0 coordinates) restricted to **chr15, chr17, chr22**
  (150,146 + 131,286 + 79,881 = 361,313 primary+secondary records by idxstats — matched in size to the gorilla
  3-contig subset; SD-rich; chr16/NPIP deliberately excluded). Genome `winloci_data/chm13v2.0.fa`.
- DE NOVO: `gw_family_catalog --homology-primary --threads 4`, today's code (commit e7bc6d7c or later with no change
  to the pipeline), default and `RUSTLE_ER_UNION_LCS=1`. Families = emitted `copies.tsv`.
- GUIDED: `mcl_families` at its current defaults with `--paf <genes all-vs-all> --gff
  winloci_data/Reference/chm13v2.0_RefSeq_full.gff.gz`; nodes = RefSeq `gene` + `pseudogene` records on the three
  chromosomes, sequences named `CHROM:START-END` (GFF 1-based); PAF = `minimap2 -x asm20 -c -X -N 50 -p 0.1` (the
  gorilla recipe, ledger §6?/9836), run with the query genes split into chunks against one full target index. Chunking
  is verified before use: on one chunk, the chunked records must equal the matching records of an unchunked run over
  the same queries. Guided reference = `<out>.clusters.tsv` (the pre-refinement clusters, same file role as
  `gw_units_v3.clusters.tsv`).
- Metrics and decision exactly as Addendum E (R_G, P_G, ARI; union narrows iff R_G up AND P_G >= default - 0.05).
- If human agrees with gorilla: recommend making the union the default (user's decision) after a genome-wide rebuild
  is run where the tool limits allow it. If human disagrees: no default change; report both.
- Limits declared: three chromosomes; downsampled library; RefSeq-based guided mode (not the gorilla core-refined units).

---
## ADDENDUM G (2026-09-13, user/advisor pivot; before any number below exists) — do the modes recover the literature subclusters of NPIP and TBC1D3?

**Truth** (`o1_falsemerge/lit/lit_truth.tsv`, CHM13 v2.0 RefSeq records):
- NPIP (22 records): level 1 = subfamily NPIPA (A1,A2,A5-A9) vs NPIPB (B1P-B15); level 2 = Dishuck 2025 paralog groups
  A2/3, A6-9, B3-5, B6-9, B12/13, all other copies singletons.
- TBC1D3 (9 protein-coding copies): level 1 = genomic cluster 1 (proximal 37.1-37.44 Mb: B,I,G,H,F) vs cluster 2
  (distal 38.91-39.06 Mb: E,K,D,TBC1D3); level 2 = Guitart/Eichler 2024 phylogenetic groups, mapped by name (declared
  assumption, the paper's per-copy CHM13 table is not on disk): AE = {TBC1D3, TBC1D3E}, CDKL = {TBC1D3D, TBC1D3K},
  B, F, G, H, I singletons.

**Input.** Full `winloci_data/A119b.t2t.bam`, all alignments inside windows = every truth record ±50 kb, plus the two
TBC1D3 cluster spans chr17:37,050,000-37,500,000 and 38,850,000-39,110,000 (never member-exon-restricted — the
subset-BAM trap). One BAM for both families, so cross-family merges are visible.

**Modes.** GUIDED = `mcl_families --paf <all RefSeq gene+pseudogene records inside the windows, all-vs-all
minimap2 -x asm20 -c -X -N 50 -p 0.1> --gff <uncompressed RefSeq restricted to the window chromosomes>` defaults.
DE NOVO default and DE NOVO union = `gw_family_catalog --homology-primary` on the window BAM, today's code, without / with
`RUSTLE_ER_UNION_LCS=1`.

**Level A — emitted families.** Each truth record gets the family of its locus (guided: its cluster via `loci.tsv`
representative; de novo: the family of the max-overlap emitted copy; none = unassigned). Per family and mode: number of
families touching the truth records, unassigned count, ARI and purity of the family labels against level 1 and level 2
(unassigned records count as singletons), and a boundary table (is NPIPA separated from NPIPB? cluster 1 from 2?).

**Level B — within-family identity structure.** Per mode, one sequence per truth record (guided: gene span; de novo:
the max-overlap copy's sequence). All-vs-all `minimap2 -x asm20 -c -X -N 50 -p 0.1`; pair identity = sum(nmatch) /
sum(blocklen) over records (0 if none). UPGMA on 1 - identity. Scores: (i) best 2-way cut vs level 1 (ARI); (ii) cut at
the Guitart criterion, divergence <= 1.5 x allelic variation = 1.5 x 15.3/10 kb = 0.0023, vs level 2 (ARI) for BOTH
families (the same rule applied to NPIP, declared; Dishuck used phylogenetic support, not a cutoff). Records with no
de novo copy are reported and excluded from level B for that mode.

**Reporting.** No pass/fail threshold is set in advance: this is a descriptive correspondence check for the advisor.
What is fixed in advance is the labels, the inputs, the mapping rules, and the two cuts above.
