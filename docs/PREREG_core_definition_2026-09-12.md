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

---
## ADDENDUM H (2026-09-13, after §6jg; before any number below exists) — can alignment COVERAGE (gene structure) separate the literature subclusters?

**Why:** §6jg found identity recovers only the tightest groups; NPIPA/B differ in gene model (Dishuck 2025) and an
earlier review saw coverage, not identity, separate A from B. User: report sensitivity, precision, bipartite matching.

**Units (one sequence per truth record, same truth table as Addendum G):**
- GUIDED: the RefSeq transcript of that gene with the longest spliced length among NM_/NR_ models, else among XM_/XR_
  models; for pseudogenes without transcripts, the exons parented by the gene. Spliced sequence, reverse-complemented
  on the minus strand.
- DE NOVO (default and union separately): POOLED — the union of the exon blocks of ALL emitted copies (any family)
  overlapping the truth record's span, spliced in genomic order, reverse-complemented if the record is on the minus
  strand. A record with no overlapping copy is missing for that mode.

**Pair statistics** (all-vs-all `minimap2 -x asm20 -c -X -N 50 -p 0.1`, all records per pair): identity =
Σnmatch/Σblocklen; coverage = min over the two sequences of (merged aligned bases / length); no alignment -> 0.

**Arms:** distance = (i) 1 - identity, (ii) 1 - coverage, (iii) 1 - coverage x identity. UPGMA. Cuts: k = 2 (root split)
scored against level 1; k = number of distinct level-2 groups among the present records, scored against level 2.
Level A (the emitted families, §6jg mapping) is re-scored in the same metrics for the record.

**Metrics (per family, per mode, per arm, per level):**
- pairwise sensitivity = truth co-member pairs placed together / truth co-member pairs;
- pairwise precision = pairs placed together that are truth co-members / pairs placed together;
- bipartite matching = one-to-one maximum-overlap assignment of predicted clusters to truth groups (Hungarian):
  micro recall = Σ matched overlap / records, micro precision = Σ matched overlap / Σ size of matched clusters,
  macro recall/precision = mean over truth groups (unmatched group = 0), exact matches = truth groups whose matched
  cluster is identical.
Missing records are excluded for that mode and reported.

**Reading, fixed in advance:** a level-1 boundary counts as RECOVERED by an arm iff the k=2 cut gives bipartite exact
matches 2/2 (NPIPA vs NPIPB; TBC1D3 cluster 1 vs 2). Anything less is reported with its numbers, not called recovered.
Descriptive otherwise; no arm is selected post hoc — all arms are reported.

---
## ADDENDUM I (2026-09-13, after §6jh; before any arm below is run) — recovering missing / fragmented members: RNA locus construction vs DNA mode vs guided

**Why:** §6jh follow-up diagnosis. All 31 NPIP/TBC1D3 truth records have >= 15 MAPQ>=1 primary IsoSeq reads, so
de novo losses are pipeline losses. TBC1D3 (chr17:39044723-39055625): 54 skeletons, 51 gate-passing transcripts, 0
reps — TBC1D3-NPEPPSP1 readthrough transcripts share junctions with both genes, the locus collapse joins them, and the
representative is a 2-exon 7-read NPEPPSP1 transcript. NPIP: records are covered by 1-8 de novo copies in 1-6 families
(fragmentation). User: split O1 into RNA and DNA modes and compare both with guided.

**Development substrate (human CHM13, descriptive):** `lit/lit.bam` (19 windows), truth `lit_truth.tsv` (31 records).
Arms, each a full `gw_family_catalog --homology-primary --threads 5` run on the same BAM, one env switch at a time:
- R0 default (existing `dn_default`);
- R1 `RUSTLE_SPLICED_REP=1`; R2 `RUSTLE_LOCUS_JUNCTION_ONLY=1`; R3 `RUSTLE_SHARED_EXON_ISOFORMS=1`;
  R4 `RUSTLE_LOCUS_EXON_UNION=1`; R5 `RUSTLE_LOCUS_GROWTH_EXTENT=1`; R6 `RUSTLE_TIER2_ADMIT=1`;
- R7 readthrough-bridge rule (new, opt-in `RUSTLE_LOCUS_BRIDGE_CUT=1`), defined now:
  on the transcripts that survive the readthrough and mis-chain filters, process spliced transcripts in order
  (n_reads desc, span desc, input index asc), keeping a union-find over junctions (chrom, donor, acceptor) whose
  component support = sum of n_reads of transcripts admitted to it. For transcript T, let D = the distinct existing
  components containing any junction of T with support > T.n_reads. If |D| >= 2, T is a BRIDGE: removed from the
  transcript set (logged); otherwise T's junctions are unioned and T.n_reads is added to the merged support.
  Unspliced transcripts are untouched. Nothing else changes.
- D1 DNA mode: `gw_family_catalog --from-genome-sd` on `HSA_sedef_pairs.bed` restricted to SD pairs with at least one
  side overlapping a lit window, against a chr16/chr17/chr18 FASTA (equal footing with the windowed BAM; the truth
  windows themselves are NOT used as DNA windows, which would be circular);
- G guided (existing `guided.*`).
(A switch that fails to run or exceeds the 390 s foreground window is reported as not run.)

**Scores per arm:** records present (>= 1 emitted copy overlapping the record); fragmentation = mean emitted copies and
mean families per present record; collapse = number of records sharing their best-overlap copy with another record;
Level A (best-overlap family per record) pairwise sensitivity/precision and bipartite micro/macro recall/precision
and exact matches vs level 1 and level 2 — (a) PRIMARY on the records present in every run arm, (b) on all 31 with
missing records as singletons. Level B (identity UPGMA, pooled units as §6jh, k=2 vs L1 and k=#groups vs L2) on the
shared present set for RNA arms; guided with BOTH spliced-transcript and gene-span units, declared now; DNA arm has no
record-specific unit and is scored at presence/fragmentation/collapse/Level A only.

**Selection on development (fixed now):** the candidate RNA rule is the arm among R1-R7 with the most records present
that loses no record present in R0; ties -> lower mean families per present record -> higher Level A L2 bipartite micro
recall on (a). If no arm beats R0 on presence, the candidate is R7 if it recovers TBC1D3, else no candidate.

**Hold-out (gorilla, not looked at for this question):** `rebuild3/sub3.bam` (3 contigs) with the candidate switch,
scored by `rebuild3/score.py` against guided gw_units_v3 exactly as Addendum E. The candidate NARROWS the de novo <->
guided gap iff R_G > cat_default's 0.0736 AND P_G >= 0.3210 - 0.05 = 0.2710; loci_with_family_copy reported. A pass
makes the switch eligible for a default flip (the user's decision); a fail keeps it opt-in and is registered.

---
## ADDENDUM J (2026-09-13, after §6ji; before any arm below is run) — (1) RNA fragmentation levers judged on the gap; (2) copy-level DNA nodes

### J1 — RNA fragmentation
**Why:** §6ji: no arm reduced NPIP fragmentation (2.6 families per record). Diagnosis (R0 `dn_default`): each gene is
one dominant spliced model plus single-exon debris, mostly '+' pieces inside '-' genes; §6as showed the '+' is the
unspliced strand placeholder, which blocks the containment collapse; §6cj/§6ck showed reads link the pieces but the
linked merge runs inside E_r blocks. Both levers were killed or parked under catalog-stability criteria (§6as, §6ck),
never under the de novo <-> guided gap. This re-evaluates them on the gap; §6as's endpoints are also reported.

**Arms** (`lit.bam`, `gw_family_catalog --homology-primary --threads 5`, one change each vs R0):
F1 `RUSTLE_COLLAPSE_UNSTRANDED=1`; F2 `RUSTLE_READ_STRAND=1`; F4 `RUSTLE_LOCUS_LINK_MIN_READS=10` (§6ck's K);
F5 = F1 + F4. Scores exactly as Addendum I (presence, mean copies and families per present record, collapse, Level A
on the shared present set and on all 31, Level B identity on shared units).

**Selection (fixed now):** among F1, F2, F4, F5, the arms that lose no record present in R0 and do not raise collapse
above R0's; the candidate is the one with the lowest mean families per present record, provided it is below R0's
2.57; ties -> higher Level A L2 bipartite micro recall on the shared set. Otherwise no candidate.

**Hold-out (gorilla `rebuild3/sub3.bam`, Addendum E scorer):** NARROWS iff R_G > 0.0736 AND P_G >= 0.2710. Also
reported, not deciding: families and copies vs cat_default; single-copy housekeeping control = number of the §6as
panel genes (bench/negative_control/check_housekeeping.py list) located on the 3 contigs in GGO_genomic.gff that
overlap any emitted family copy, candidate vs cat_default.

### J2 — copy-level DNA nodes (atoms from SD alignment boundaries, edges from the SD alignments)
**Why:** §6ji D1 merged overlapping SD intervals into 170-380 kb nodes (21/31 records collapsed) and spent 96 of 102
minutes re-aligning them with minimap2.

**Nodes (atoms), fixed now:** input = the SD pairs (SEDEF) with both sides on the substrate's contigs [human dev:
`lit/dna/sd_lit.bed` exactly as D1; gorilla hold-out: `GGO_sedef_final.bed` pairs with both sides on the 3 rebuild3
contigs]. Per chromosome, breakpoints = every SD side start and end; elementary segments = intervals between
consecutive breakpoints that lie inside at least one side; signature = set of (pair, side) covering the segment.
Adjacent (touching) segments with identical signatures are merged. Then every segment shorter than 1,000 bp
(`GenomeRepParams::min_block`) is merged into the touching neighbour whose signature has the larger Jaccard similarity
to its own (tie -> left; no touching neighbour -> dropped), repeated shortest-first until none remains below 1,000 bp.

**Edges, fixed now:** identity of a pair = matches / (matches + mismatches) from the SEDEF columns (checked equal to
the file's fracMatch column on every row; any mismatch > 1e-4 aborts). For each SD pair, each atom x overlapping side A
is projected onto side B by linear interpolation of relative position (reversed when the strands differ); every atom
y on side B overlapping the projection receives covered bases on both x and y (interval unions per atom pair, over
all SD pairs, in both directions A->B and B->A). Atom-pair identity = covered-length-weighted mean of the contributing
pairs' identities. Coverage = min(covered_x / len_x, covered_y / len_y). x == y is skipped.

**Families:** the production DNA grouping with those edges instead of minimap2: edges with identity >= 0.80 and
coverage >= 0.50 (the `--from-genome` floors), weight 1.0, `gamma_quasi_clique_partition` gamma 0.20, then
`distinct_locus_reps_grouped` (min_reads 0) and >= 2 copies (new `families_from_edges`, same code path after the
edge list).

**Development (human, lit):** arm D2 scored exactly as Addendum I (presence, fragmentation, collapse, Level A).
Reading fixed now: "the NPIPA/NPIPB separation holds at copy level" iff D2 collapse <= 5 of 31 AND D2 NPIP Level A L1
pairwise precision >= 0.775 (D1's 0.875 - 0.10); otherwise it does not.

**Hold-out (gorilla rebuild3 contigs, Addendum E scorer vs guided gw_units_v3):** DNA atoms are CLOSER to guided than
the RNA default iff R_G > 0.0736 AND P_G >= 0.2710. Descriptive mode comparison; no default changes.

---
## ADDENDUM K (2026-09-13, after §6jj; before any number below exists) — atom edge fix, DNA+RNA hybrid, genome-wide hold-out

### K0 — atom edges from exact aligned blocks (post-hoc fix, disclosed)
**Why:** §6jj's housekeeping control fired once (ATP5F1A). Inspection: its 7.7 kb atom pairs with a 1.6 kb processed
pseudogene; the SEDEF alignment is exons vs retrocopy with intron-sized deletions (1,620 aligned bases), but J2's linear
interpolation between side endpoints scored coverage 1.00 (true ~0.21). This fix was chosen AFTER seeing the rebuild3
hold-out, so rebuild3 is development from here on.
**Rule:** gorilla SEDEF CIGAR (column 33; validated on 8 rows: M aligned, D consumes side A, I consumes side B, side B
reverse-complemented when strand2 = '-', identity over M reproduces fracMatch to 1e-6). Each M run is an exact block
A[pa, pa+n) <-> B genomic [b1+pb, b1+pb+n) ('+') or [b2-pb-n, b2-pb) ('-'); atoms are intersected with blocks and
mapped base-to-base. Covered bases per atom pair = union of mapped aligned bases on each atom; identity = pair fracMatch
weighted by aligned bases; coverage = min over both atoms (J2 form; coverage-of-shorter reported as secondary). Atoms,
floors, gamma and grouping unchanged. A row without a CIGAR aborts the run.
**Development (rebuild3):** report R_G/P_G (any-family and best-overlap), loci with a family copy, housekeeping.

### K1 — DNA+RNA hybrid (atoms as nodes, reads mark expression)
**Expression:** u(interval) = primary (-F 2308) MAPQ >= 1 reads with >= 1 aligned base (M/=/X) inside the interval;
EXPRESSED iff u >= 3 (`GATE_MIN_READS`). MAPQ-0 counts reported alongside, not used.
**Hybrid catalog H:** K0 DNA families restricted to expressed atoms; families with >= 2 expressed atoms are kept.
**Truth:** guided `gw_units_v3` restricted to expressed loci (u over the guided locus span >= 3), clusters with >= 2
expressed loci.
**Compared against that truth:** RNA de novo (dev: rebuild3 `cat_default`; hold-out: `o1_reps/ggo_reps.copies.tsv`,
production genome-wide catalog of 2026-08-21, older code — disclosed), K0 DNA (all atoms) and H. R_G, P_G any-family
(Addendum E scorer) and best-overlap.
**Decision (hold-out only):** H NARROWS the RNA <-> expressed-guided gap iff R_G(H) > R_G(RNA) AND P_G(H) >= P_G(RNA)
- 0.05. Also reported: RNA copies overlapping an expressed atom; expressed guided loci with no expressed atom.

### K2 — genome-wide hold-out
**Substrate:** every gorilla contig EXCEPT NC_073241.2 / NC_073242.2 / NC_073244.2. SEDEF pairs with either side on
those 3 contigs are excluded; guided loci on them are excluded before the >= 2-loci cluster rule; RNA copies on them
are excluded. Expression from `GGO_ds.bam` (the BAM of the RNA catalog).
**Decisions:** (i) K0 DNA atoms are CLOSER to guided than RNA iff R_G(DNA) > R_G(RNA) AND P_G(DNA) >= P_G(RNA) - 0.05
(any-family scorer); (ii) K1 as above. Housekeeping: the 30-gene panel genome-wide, genes in any family, per catalog.
Runs that exceed resources are reported as not run.

---
## ADDENDUM L (2026-09-13, after §6jk; before any number below exists) — (1) atoms + RNA-only families; (2) genome-wide RNA rebuild with today's code

### L2 — genome-wide RNA rebuild (runs first; its catalogs are the RNA side of every hold-out below)
`gw_family_catalog --bam GGO_ds.bam --fasta GGO.fasta --out <p> --homology-primary --threads 4` built at HEAD, with
`RUSTLE_ER_EDGE_DUMP` (same invocation as `o1_reps/run.sh`, 2026-08-21): arm RD default, then arm RU with
`RUSTLE_ER_UNION_LCS=1`. A memory guard kills a run (by PID) if MemAvailable + SwapFree < 2 GB; a killed or failed arm
is reported as not run. Scored on the hold-out contigs (all except the 3 rebuild3 contigs) with
`bench/score_vs_guided.py`.
Decisions: (a) K2 and K1 re-evaluated with RD in place of the 08-21 catalog, same bars (DNA closer iff R_G(DNA) >
R_G(RD) and P_G(DNA) >= P_G(RD) - 0.05; hybrid narrows iff the same holds for H against expressed guided); (b) the
Addendum E union question genome-wide: RU NARROWS iff R_G(RU) > R_G(RD) AND P_G(RU) >= P_G(RD) - 0.05 (full guided).

### L1 — hybrid atoms plus RNA-only families
**Why:** §6jk: 22.6% of expressed guided loci have no atom and 36.5% of RNA copies lie outside atoms.
**Arms** (H = §6jk hybrid, expressed K0 atoms; R = the RNA catalog: dev rebuild3 `cat_default`, hold-out RD):
- U1 = H + R's families restricted to copies overlapping NO atom (any expression), kept with >= 2 such copies;
- U2 = H + R's families restricted to copies overlapping no EXPRESSED atom, kept with >= 2 such copies;
- U3 = H + every R family, whole, that has >= 1 copy overlapping no expressed atom.
Family ids are kept disjoint (R families prefixed). Truth = expressed guided loci (K1 definition).
**Selection on development (rebuild3, fixed now):** among U1-U3, those with P_G(any) >= P_G(H) - 0.05; the candidate
is the one with the highest R_G(any); ties -> higher P_G. If none qualifies, no candidate.
**Hold-out decision (23 contigs, R = RD):** the candidate NARROWS the hybrid's gap iff R_G(U) > R_G(H) AND P_G(U) >=
P_G(H) - 0.05 (any-family); best-overlap reported. If RD is not run, the hold-out uses the 08-21 catalog, disclosed.

**L1 ruling (2026-09-13, after development scoring, before the hold-out):** U1 and U2 tie on development in R_G(any)
0.6383 and P_G(any) 0.6375 (H 0.6356 / 0.6365; U3 P_G 0.4066 fails the precision guard). The fixed tie-break (higher
P_G) does not separate them. Ruling: the candidate is U1, the arm with fewer copies (854 vs 864; U1's exclusion rule is
the stricter one). Cost if wrong: U2 is not evaluated on the hold-out; its hold-out numbers are reported as secondary.

**ADDENDUM L — NOT RUN (2026-09-13, user decision).** The genome-wide RD rebuild was stopped ~20 min in (killed by
PID) and L1/L2 were cancelled after the user redefined scope: no genome-only discovery mode; de novo = cluster loci that
already have aligned reads; guided = start from annotation (even minimal), find new candidate loci, judge them vs ground
truth in width, derive under/over-merge rules. Only L1's development numbers exist (rebuild3; reported, no decision).

---
## ADDENDUM M (2026-09-13, after §6jl; before any number below exists) — guided mode leave-out on NPIP and TBC1D3 (descriptive first look)

**Scope (user):** guided mode starts from a minimal annotation, finds new candidate loci, and is judged against ground
truth in width (locus boundaries AND family breadth). NPIP and TBC1D3, human CHM13.

**Truth:** `lit/lit_truth.tsv` (NPIP 22 records, TBC1D3 9); family = NPIP / TBC1D3; true locus width = RefSeq gene span.

**Seed sequence per record:** the §6jh guided unit (longest NM_/NR_ transcript, else XM_/XR_, else gene-parented exons,
spliced, transcript orientation); a record with no model (NPIPB14P) uses its gene span. All 31 units are aligned once:
`minimap2 -c -x splice -N 100 -p 0.1 -t 4` against the prebuilt `chm13v2.0.fa.mmi` (k15 w10, not the splice preset's
w5 — disclosed). Hit identity = nmatch / block length; query coverage = aligned query fraction.

**Leave-out:** per family, records sorted by name and shuffled with `random.Random(1000 * replicate + family_index)`,
5 replicates. PRIMARY level: keep ceil(n/2) as seeds (NPIP 11, TBC1D3 5), hide the rest. SECONDARY level: keep 1.

**Candidate loci (baseline rule = the DNA E_r floors):** hits from seed queries with identity >= 0.80 and query
coverage >= 0.50, whose target span overlaps no seed record's gene span (either family). Overlapping hit spans are
single-linkage clustered; each cluster's candidate = its highest-nmatch hit (span = that hit's target span, family =
that hit's seed family).

**Classification:** a candidate overlapping a hidden record is matched to the hidden record it overlaps most; otherwise
it is `other_gene` if it overlaps any gene/pseudogene in `chm13v2.0_RefSeq_full.gff.gz`, else `unannotated`.
A hidden record is RECOVERED if a candidate overlaps it; its predicted family is that candidate's family.

**Scores (per family, per level, mean and SD over replicates; also pooled):**
- Breadth: sensitivity = hidden recovered with the correct family / hidden; under-merge = hidden missed + hidden
  recovered into the wrong family; precision = candidates whose matched hidden record is in their family / candidates
  assigned to the family; over-merge = candidates of the family that match another family's record, `other_gene`
  (named) or `unannotated`. Pairwise sensitivity/precision and bipartite matching over items = hidden records (missed
  ones as singletons) plus unmatched candidates (each its own truth singleton).
- Width (hidden records recovered into the correct family): span Jaccard; 5' and 3' boundary offsets in bp relative to
  the truth strand (positive = extends beyond the truth boundary); truncated = candidate covers < 0.90 of the truth
  span; overextended = candidate extends beyond the truth span by > 0.10 of its length.
Descriptive; no rule is tuned or selected here. Rules against under/over-merges are pre-registered separately.

---
## ADDENDUM N (2026-09-13, after §6jm; before any number below exists) — guided-mode rules: iterative expansion and 5' width extension (NPIP, TBC1D3)

**Why:** §6jm: with one seed, NPIP expansion stops at subfamily lines (under-merge); width error is 5' truncation.
Diagnosis before this addendum: 9/31 records' own seed transcript leaves the 5' end of their own gene span uncovered
(gene span = union of isoforms; NPIPB9 0.564, TBC1D3G 0.738, NPIPB5 0.777), so part of the truncation is definitional.

**Common to every arm:** Addendum M's truth, leave-out sets (same RNG, levels keep-50% and keep-1, 5 replicates),
floors (identity >= 0.80, query coverage >= 0.50), blocking (hits overlapping a seed gene span are ignored), single-
linkage clustering with the highest-nmatch hit deciding the candidate's family, classification, and all Addendum M
scores. Added score (secondary): family-named precision, counting `other_gene` candidates whose overlapping RefSeq gene
name or description names the family ("nuclear pore complex-interacting protein"/NPIP; "TBC1 domain family member
3"/TBC1D3) as members.

**Arms:**
- M0 — Addendum M (reference; recomputed by the new script, must reproduce §6jm).
- W1 gene-span projection — extra queries: each seed's gene span (genomic, transcript strand), `minimap2 -c -x asm20
  -N 100 -p 0.1`; a candidate's span becomes the union of its best transcript hit span and the passing gene-span hits
  of the same seed gene that overlap it. Candidate set = M0's (W1 changes width only).
- W2 isoform union — queries: every exon-bearing annotated transcript of each seed gene (139 over 28 genes; NPIPB10P/
  NPIPB1P gene-exons, NPIPB14P gene span), splice mode as M. Candidates are clustered from all passing isoform hits; a
  candidate's span = union of the passing hits in its cluster from the best hit's seed gene.
- I iterative expansion — round 0 = M0. Each new candidate yields a query: the target sequence of its best hit's
  aligned exon blocks (cg CIGAR runs between N, in the hit's transcript orientation), aligned as M. Passing hits that
  overlap no seed gene and no existing candidate form new candidates (single-linkage, best hit) and inherit the family
  of the candidate whose query produced the best hit. Repeat to a fixed point or 10 rounds. Width per candidate = W0.
- I+W2 — round 0 = W2, then iteration as I.
Width ceilings reported per arm's width rule: each record's self-projection coverage of its own gene span.

**Readings (fixed, descriptive, no selection or tuning):** (1) iteration CROSSES the NPIP subfamily barrier iff keep-1
NPIP sensitivity > 0.50 with family-named precision 1.000; (2) a width rule REDUCES truncation iff its truncated count
is lower than M0's at both levels and its overextended count rises by less than its truncated count falls; (3) any
cross-family assignment or non-family-named `other_gene` candidate is an over-merge and is listed.
Rules are developed and read on these two families only; confirmation on another substrate is a later addendum.

---
## ADDENDUM O (2026-09-13, after §6jn; before any number below exists) — the literature's two levels: G1 family by gene body, G2 subfamilies from exon-masked gene-body identity (NPIP, TBC1D3)

**Why:** Dishuck 2025 defines NPIP copies by aligning the ~19 kb gene body (>= 80% identity, >= 15 kb aligned) and
paralogs/subfamilies as ML clades on the MSA with exons, VNTRs and poorly aligned regions removed ("15 kbp of intronic
sequence"). §6jn: transcript-level under-merge; gene-body alignment reaches 104/105 NPIPA->NPIPB pairs.

**Common:** Addendum M truth, leave-out sets (same RNG; keep 50% and keep 1; 5 replicates), blocking by seed gene
spans, clustering (single-linkage, highest total nmatch decides family), classification and all M/N scores incl.
family-named precision. M0 is recomputed as the reference.

**G1 — family by gene body.** Query = each seed's gene span (genomic, transcript orientation), the existing
`genespan.paf` (`minimap2 -c -x asm20 -N 100 -p 0.1`, k15/w10 index). Records of one query on one chromosome and strand
are CHAINED greedily in target order when the target gap to the chain end is <= L_q (query length), query order is
consistent with the strand, and the chain's target span stays <= 2 L_q. Chain identity = sum nmatch / sum block length;
aligned = union of query intervals; L_t = extrapolated target span (below). ACCEPT iff identity >= 0.80 AND aligned >=
0.50 x min(L_q, L_t). Width, two variants: G1-clip = chain target span; G1-extrap = chain target span extended by the
unaligned query ends (strand-aware) — i.e. the projection of the seed's first and last exon boundaries (= gene span
ends). Breadth is identical for both; width is scored for both. Subfamily reach (NPIPA/NPIPB of recovered records) is
reported.

**G2 — subfamilies inside a G1 family.** Members = seeds + G1 candidates of that family (per replicate and level) and,
as the no-leave-out REFERENCE, all truth records of the family. Member sequence = its gene body with exons removed:
seeds (and reference records) lose their annotated union exons (all isoforms; gene-parented exons; NPIPB14P has none
and is not masked); a candidate's region is its G1-extrap span, and the seed's union exons are projected base-to-base
through the chain CIGAR and removed (exons outside aligned blocks cannot be projected and stay). Pieces are
concatenated in transcript orientation. Pairwise identity = sum nmatch / sum block length over all records of
`minimap2 -c -x asm20 -X -N 50 -p 0.1` all-vs-all (both directions pooled). Partition = `bench/identity_gap.py`'s rule
(largest interior gap, outer 10% ignored; gauss/beta/smooth nulls, 10,000 draws, worst p governs): p < 0.05 -> connected
components of pairs with identity >= the gap midpoint; otherwise one group.
Scored on truth-record members vs level 1 (NPIPA/NPIPB; TBC1D3 cluster1/2) and level 2: pairwise sensitivity/precision
and bipartite micro R/P and exact matches.

**Readings (fixed):** (1) G1 ADDRESSES under-merge iff keep-1 NPIP sensitivity > 0.50 with family-named precision >=
0.95 and 0 cross-family assignments. (2) A G1 width variant REDUCES truncation iff (Addendum N reading 2) vs M0.
(3) G2 RECOVERS NPIPA/NPIPB in a run iff p < 0.05 and bipartite exact 2/2 vs level 1 on its truth-record members;
reported for the reference and as a count over leave-out runs. (4) G2 is CORRECT on TBC1D3 iff it does NOT split
(p >= 0.05), since Guitart's clusters are positional (§6gw). Descriptive; no tuning.

---
## ADDENDUM P (2026-09-13, after §6jo; before any number below exists) — P1 subfamilies as tree clades; P2 hybrid width (NPIP, TBC1D3)

### P1 — G2 as a tree (closer to Dishuck 2025)
**Inputs:** exactly Addendum O's G2 member sequences (exon-masked gene bodies; reference = all truth records per family,
leave-out = seeds + G1 candidates, 10 runs per family), regenerated by the same code; candidate members are labelled by
their matched hidden record as in O.
**Tree:** `mafft --retree 2` (FFT-NS-2, the paper's mode), then columns with gaps in > 50% of sequences are removed (an
automated stand-in for the paper's visual trimming; VNTRs are not masked — disclosed), then `iqtree3 -m MFP -B 1000
-alrt 1000 -T 4 --seed 1`. If one tree exceeds 10 minutes, the model is fixed to GTR+F+R4 for all runs (disclosed).
**Supported split:** SH-aLRT > 75 (the paper's criterion); UFBoot >= 95 reported alongside.
**Clade recovery:** a literature group G (>= 2 truth members present) is RECOVERED iff the unrooted tree has a supported
split whose truth-member side equals G exactly (candidate leaves without a truth label are ignored on either side).
Groups: NPIP level 1 NPIPA|NPIPB (one split); level 2 A6-9, B3-5, B6-9, B12/13; the paper's named "human-specific
NPIPB subfamily" {B3, B4, B5, B11, B12, B13}. TBC1D3: AE {TBC1D3, E} and CDKL {D, K}; cluster1|cluster2 is POSITIONAL.
**Partition for sensitivity/precision/bipartite (level 1):** the two child clades of the midpoint-rooted tree.
**Readings (fixed):** (1) NPIPA/NPIPB recovered in a run iff its split is supported — counted over the reference and
the 10 leave-out runs; (2) the recovered fraction of the level-2 and named groups is reported per run; (3) TBC1D3 is
CORRECT in a run iff no supported split has cluster 1 (or cluster 2) as its truth side.

### P2 — hybrid width (transcript boundary where a paralog has one, gene body otherwise)
**Candidates:** Addendum O's G1 candidates (same leave-out sets). **Width:** transcript hits = Addendum M passing hits
(identity >= 0.80, query coverage >= 0.50) from any seed of the candidate's family, not overlapping a seed gene, whose
target span overlaps the candidate's G1-clip span. H = the highest-nmatch such hit's span if any, else the G1-clip span.
H-union (secondary) = the union of all such hit spans, else G1-clip. Reported: fraction of candidates with a transcript
boundary; width on all recovered records and on the records M0 also recovers (like-for-like).
**Reading (fixed):** H IMPROVES on G1 iff its truncated + overextended count on all recovered records is lower than
G1-clip's at both levels; its like-for-like Jaccard vs M0 is reported.

---
## ADDENDUM Q (2026-09-13, while Addendum P's MAFFT run is still going and before any P1 or P2 number is seen) — subfamily trees from a reference-projected alignment

**Why (user):** MAFFT is not recommended for long gene bodies — a forced global collinear alignment of 15-50 kb
repeat-containing, differently sized loci is fragile and scales poorly. P1 is finished as registered as a
literature-matching reference point; Q is the scalable alternative, run on the same members.

**Members and sequences:** exactly P1's per-run members (exon-masked gene bodies, same FASTA files; reference runs and
10 leave-out runs per family).
**Reference member (fixed rule, no truth):** all-vs-all `minimap2 -c -x asm20 -X -N 50 -p 0.1`; the reference is the
member with the largest total aligned query+target bases summed over its pairs (ties -> name order).
**Projection:** every other member is aligned to the reference with `minimap2 -c -x asm20 -N 50 -p 0.1` (member as
query). Records are applied in decreasing alignment score (AS); each reference column takes the member base from the
first record that aligns it (CIGAR walk: M copies the member base, D leaves a gap, I is dropped; reverse-strand records
use the reverse-complemented member). Unfilled columns are gaps. Columns covered in <= 50% of members are removed (P1's
rule). The reference row is its own sequence.
**Tree and scores:** exactly P1 (`iqtree3 -m MFP -B 1000 -alrt 1000 -T 4 --seed 1`, SH-aLRT > 75, clade recovery,
positional check, midpoint-root level-1 partition, readings 1-3 of P1).
**Comparison with P1:** per (run, literature group) cell, agreement of the RECOVERED / not-recovered call; per run,
alignment columns kept and wall time. Descriptive.

---
## ADDENDUM R (2026-09-13, after §6jp; before any AMY number exists) — the guided pipeline on the amylase family (human CHM13), rules unchanged

**Why (user):** test the NPIP/TBC1D3-derived guided pipeline on all AMY genes.
**Literature level (read before scoring):** Bolognini et al. 2024 Nature (PMC11485256) classify copies as AMY1 (salivary),
AMY2A and AMY2B (pancreatic), AMY2Ap (partial AMY2A lacking ~4.5 kb of the 5' end) and AMYP1; their haplotype trees use
unique FLANKING sequence, not gene copies. Yilmaz et al. 2024 Science: AMY2B, AMY2A and AMY1 carry 23, 23 and 36 fixed
coding variants unique to each type. ⟹ the amylase subclusters are GENE TYPES, not intronic clades as for NPIP.

**Truth (fixed now):** every gene/pseudogene in `chm13v2.0_RefSeq_full.gff.gz` named AMY* or described as amylase on
chr1p21 (103.3-103.9 Mb) — 12 loci (AMY2B, AMY2A, AMY1A, AMY1B, AMY1C, 4 "alpha-amylase 1B" LOCs, 2 "pancreatic
alpha-amylase-like" LOCs, AMYP1); width = RefSeq gene span. Level 1 = salivary AMY1 (AMY1A/B/C + the four
"alpha-amylase 1B" LOCs) | pancreatic AMY2 (AMY2A, AMY2B, the two "pancreatic alpha-amylase-like" LOCs) | AMYP1.
Level 2 = AMY1 | AMY2A | AMY2B | pancreatic-like | AMYP1. Types of the LOCs come from RefSeq product names (name
trap disclosed, §6c); LOC124905662's span (36 kb) contains LOC128966568 — a known overlap, not corrected.
MGAM/MGAM2 (alpha-amylase domain, different family) are a named over-merge check.

**Pipeline and scores (identical rules, no tuning):** seeds/units as Addendum M (longest curated transcript, spliced),
isoforms and gene spans as Addendum N, leave-out keep 50% / keep 1 with 5 replicates (same RNG rule); arms M0 (Addendum
M), G1 family by gene body (Addendum O), hybrid width H (Addendum P2), subfamily clades from the reference-projected
exon-masked gene-body tree (Addendum Q; supported = SH-aLRT > 75); family-named precision uses "amylase|AMY".
Literature groups scored as clades: level 1 AMY1 (7) and pancreatic AMY2 (4); level 2 pancreatic-like pair.
**Readings (fixed):** (1) G1 breadth at keep 1: sensitivity and over-merges (MGAM/MGAM2 or other non-amylase genes);
(2) H vs G1 width by Addendum P2's rule; (3) AMY1 and pancreatic clades recovered counts over the reference and 10
leave-out runs. Descriptive.

---
## ADDENDUM S (2026-09-13, after §6jq; before any number below exists) — S1 candidates from either finder; S2 subfamily trees on exons and on introns (NPIP, TBC1D3, AMY)

**AMY truth v2 (evidence gathered before this addendum, from sequence and reads, not from pipeline scores):**
LOC124905662 (XM_047443612.1, Gnomon model, 36.2 kb, '-'): exons 2-8 (chr1:103,575,076-103,581,258) = AMY2A exons 4-10 at
99.6-100% identity (blastn); exon 1 (103,611,079-103,611,298) matches AMY2A exon 3 (100% over 198 bp) and lies 29.8 kb
upstream across the entire opposite-strand AMY1B-like LOC128966568; no A119b.t2t.bam (testis IsoSeq) read has the
103,581,258 -> 103,611,079 junction (3 reads touch exon 1's block, 3 exon 2's; testis is not the expressing tissue).
LOC124905664 (6.2 kb) = AMY2A exons 4-10 as well. Both match Bolognini's AMY2Ap (partial AMY2A lacking the 5' end).
v2 = v1 with LOC124905662's span set to its exons 2-8 block and level 2 of both LOCs = AMY2Ap; units/isoforms stay the
annotation's transcripts; LOC124905662's gene-span query is re-aligned for the new span. **v1 (Addendum R) and v2 are
both run and both reported.**

### S1 — candidates from either finder
Hits = Addendum M passing transcript hits (units) ∪ Addendum O passing gene-body chains, from seeds, ignoring any that
overlap a seed gene span. Single-linkage clusters over the hits' target spans (transcript hit span; chain clip span).
Family = the seed family of the cluster's highest-nmatch TRANSCRIPT hit if the cluster has one, else of its
highest-nmatch chain. Width = H (Addendum P2): that transcript hit's span if present, else the chain's clip span.
Scores: Addendum M/O breadth and width, next to M0 and G1 recomputed.
**Reading (fixed):** the union ADDRESSES under-merge across families iff, at keep 1, its sensitivity is >= max(M0, G1) -
0.02 for NPIP, TBC1D3 and AMY (v1 and v2), with family-named precision >= 0.95 and 0 cross-family assignments.

### S2 — subfamily trees on two sequence classes
Members = seeds + S1 candidates of the family (leave-out) or all truth records (reference). Per member, two sequences:
- EXON class: truth/seed records — annotated union exons spliced in transcript orientation; candidates — the target exon
  blocks of the cluster's best transcript hit (cg runs between N, hit orientation) if present, else the seed's union
  exons projected through the best chain's CIGAR, spliced;
- INTRON class: truth/seed records — gene span minus union exons (Addendum O); candidates — the best chain's extrap span
  minus projected seed exons if a chain exists, else the best transcript hit's span minus its exon blocks.
Each class: Addendum Q's reference-projected alignment (reference = most total aligned bases) + `iqtree3 -m MFP -B 1000
-alrt 1000 -T 4 --seed 1`; supported = SH-aLRT > 75. Groups: NPIP (P1's six), TBC1D3 (AE, CDKL, positional check),
AMY v1 (AMY1, AMY2, pancreatic-like) and v2 (AMY1, AMY2, AMY2Ap).
**Readings (descriptive):** per family and class, recovered counts over the reference and leave-out runs, plus
"either class"; expectation stated in advance from the literature — NPIP subfamilies from introns, AMY types from exons.

---
## ADDENDUM T (2026-09-13, after §6jr; before any number below exists) — guided fixes F1-F4, one tool, re-run on NPIP/TBC1D3/AMY v2

**Why (user):** apply all necessary fixes so the approach works in guided (this addendum) and de novo (next).
**Tool:** `bench/guided_pipeline.py` — builds seed units and gene-body queries from the GFF, aligns them, runs the
leave-out, scores, and builds trees; no exec-reuse of earlier scripts. Its seed units must be identical to the units
used in §6jm-§6jr (checked; a mismatch aborts).

**F2 — gene-body query = CDS envelope.** Per gene: first to last coding base over all its transcripts' CDS features; if
no CDS, the exon envelope of the seed unit's transcript; if none, the gene span. Genomic sequence in transcript
orientation; `minimap2 -c -x asm20 -N 100 -p 0.1`. Chaining, identity >= 0.80 and aligned >= 0.50 x min(L_q, L_t),
extrapolation and exon projection exactly as Addendum O (exons = union exons clipped to the envelope). Truth widths
stay the RefSeq gene span; blocking stays seed gene spans.

**F1 — candidate construction.** (1) Transcript hits (Addendum M floors) are leader-clustered: in decreasing nmatch,
a hit joins the first leader whose span it overlaps by >= 0.50 of BOTH spans (reciprocal), else becomes a leader.
(2) Gene-body chains likewise, on clip spans. (3) A transcript leader and a chain leader are the same locus iff their
spans overlap reciprocally >= 0.50, or the transcript leader lies >= 0.90 inside the chain's extrapolated span and that
chain contains no other transcript leader; components by union-find. Each component is one candidate: family and width
from its highest-nmatch transcript leader if any (hybrid width), else from its highest-nmatch chain (clip span).
Reported next to M0 and G1 as in §6jr (their single-linkage construction unchanged) and to §6jr's U.

**F3 — tree reference member.** The member aligned (>= 1 record in the all-vs-all) to the most other members; ties ->
most total aligned bases -> name. The all-gap guard stays.

**F4 — subfamily step.** Exon and intron classes as Addendum S (intron class = gene-body unit minus exons: seeds and
reference records use the F2 envelope minus union exons; candidates their chain's extrapolated span minus projected
exons, else transcript span minus exon blocks). Reference-projected alignment + IQ-TREE, SH-aLRT > 75; clades reported per
class and "either".

**Bars (fixed):** (B1) for NPIP, TBC1D3 and AMY v2, U sensitivity >= max(M0, G1) - 0.02 at BOTH keep 50% and keep 1,
family-named precision >= 0.95, 0 cross-family; (B2) either-class clade recovery counts >= §6jr's for every literature
group (NPIP six groups; AMY v2 AMY1, AMY2, AMY2Ap), counting reference and leave-out runs; (B3) TBC1D3 positional split
not supported in the exon tree in all runs. A failed bar is reported with its mechanism; nothing is tuned here.

---
## ADDENDUM U (2026-09-13, after Addendum T's results; post-hoc fixes, disclosed) — chain-first candidates; gene-span intron sequence

**T outcome that motivates U:** B1 PASS on NPIP, TBC1D3, AMY v2 (U sensitivity 1.000 at both levels, no over-merge);
B3 PASS; **B2 FAIL** — either-class clade counts below §6jr for NPIP A6-9 (8 vs 10), B3-5 (9 vs 11), B6-9 (10 vs 11),
B12/13 (10 vs 11) and AMY AMY1 (8 vs 10), AMY2 (3 vs 5). Diagnosed: (i) DUPLICATE candidates — transcript hits spanning
22-73 kb whose first exon aligns to another copy (AMY) and partial transcript hits (NPIP) are not reciprocal-overlap
merged with the locus's gene-body chain, so one locus enters trees twice (NPIP keep-50% candidates 25.6 vs 15.6; AMY
keep-1 12.8 for 11 hidden); (ii) the INTRON class was built on the CDS envelope (F2), shorter than the gene span §6jr used.

**U1 — chain-first construction (replaces F1).** Gene-body chains (F2 CDS-envelope queries) are leader-clustered by
reciprocal overlap >= 0.50 on clip spans; each chain leader is a locus. Each transcript hit (decreasing nmatch) is
ATTACHED to the chain leader whose extrapolated span contains >= 0.90 of it (largest overlap, then highest nmatch);
a transcript hit that overlaps any chain leader's extrapolated span but is not contained is DISCARDED (mis-chained or
partial); transcript hits overlapping no chain leader are leader-clustered by reciprocal overlap and form
transcript-only loci. Candidate width and family: the best attached transcript hit if any, else the chain (clip span);
transcript-only loci use their leader.
**U2 — intron sequence from the gene span.** Seeds' annotated gene spans are also aligned (`asm20`, as `genespan.paf`
of §6jo) and chained (Addendum O). Intron class: seeds/reference records = annotated gene span minus union exons (§6jr);
candidates = the extrapolated span of the gene-span chain overlapping the candidate most, minus that seed's union exons
projected through it; if none, the CDS-envelope chain's; if none, transcript span minus its exon blocks. Exon class as T.
F2 (finding), F3 (tree reference) unchanged.
**Bars:** B1, B2, B3 exactly as Addendum T, plus **B4**: hidden records overlapped by more than one candidate <= 5% of
recovered hidden records at both levels. Families NPIP, TBC1D3, AMY v2. Development only — the same families shaped
these fixes, so a pass here is not validation; the de novo hold-out follows.

**Addendum U amendment U1' (2026-09-13, after the first U run on AMY v2 only; disclosed):** that run gave B1 PASS but
B4 FAIL at keep 50% (0.233): each duplicate was a complete chain and a PARTIAL chain from a more divergent seed at one
locus, reciprocal overlap < 0.50. Chain leader clustering now also joins a chain whose clip span lies >= 0.90 inside a
leader's extrapolated span (the same containment rule transcript hits use). Nothing else changes; NPIP/TBC1D3 have not
been run under U yet; all three families are re-run under U1'.

**Addendum U amendment U1'' (2026-09-13, after U1' results on all three families; disclosed; LAST change before the
de novo step):** U1' gave NPIP/TBC1D3 B1 PASS, B4 PASS, B3 PASS, B2 NPIP 5/6 groups >= §6jr (A6-9 9 vs 10); AMY v2 B1
FAIL at keep 50% (0.933 vs 0.967: LOC124905662's mis-joined 36 kb CDS envelope contains LOC128966568 and absorbs it),
B2 PASS, B4 PASS. Width regression: TBC1D3 transcript hits (~10.9 kb incl. UTRs) exceed the 9 kb CDS-envelope chain, are
never "inside" it and were discarded (keep-50% Jaccard 0.969 -> 0.806). Change: a transcript hit is also ATTACHED to a
chain leader when that chain's clip span lies >= 0.90 inside the transcript hit AND the transcript hit overlaps no other
chain leader; otherwise the U1' rules apply. All three families are re-run once under U1''; results are reported
whatever they are, and no further guided change is made in this round.

---
## ADDENDUM V (2026-09-13, after the guided U1'' runs; before any de novo number below exists) — de novo: gene-body edge union (D1) and exon/intron subfamily clades (D2)

**Scope reminder (user):** de novo = cluster loci that already have aligned RNA reads; its fair truth is read-supported
loci. Guided lessons carried over: NPIP is held together by gene body, AMY by coding sequence; subfamilies are clades of
intron (NPIP) or exon (AMY) trees.

**D1 (Rust, opt-in, default off, byte-identical when unset):** `RUSTLE_ER_UNION_GENOMIC_SPAN=1` adds to the exon-sum
E_r edges the E_r edges computed on each rep's genomic span (the existing `homology_genomic_span` substrate, same floors,
coverage of the shorter), before the gamma partition; existing pairs keep their metrics. 2 new unit tests.
**D2:** `bench/denovo_subfamilies.py` — for each emitted family containing a truth record: members = the family's
copies; each truth record is represented by its best-overlapping copy of that family (other copies unlabelled); exon
sequence = the copy's spliced sequence from `copies.fa`; intron sequence = the copy's genomic span minus its exon blocks
(transcript orientation); trees and clade calls exactly as `bench/guided_pipeline.py` (reference-projected alignment,
F3 reference, IQ-TREE MFP + 1000 UFBoot/SH-aLRT, SH-aLRT > 75, clades per class and either).

**Development (human CHM13 IsoSeq A119b.t2t.bam):** (a) the 19 NPIP/TBC1D3 windows (`lit/lit.bam`), (b) an amylase window
(chr1:103,300,000-103,900,000 from `A119b.t2t.bam`, `amy.bam`). Arms: DN0 `gw_family_catalog --homology-primary
--threads 5` (current code); DN1 = DN0 + `RUSTLE_ER_UNION_GENOMIC_SPAN=1`. Scores: records present, mean families and
copies per present record, collapse; family level (best-overlap family per record vs family truth; AMY v2 truth):
pairwise sensitivity/precision and bipartite micro R/P; D2 clade calls per literature group (NPIP six, TBC1D3 AE/CDKL +
positional, AMY AMY1/AMY2/AMY2Ap). Descriptive; expectation stated now: DN1 lowers NPIP families per record without
raising collapse.
**Hold-out (gorilla rebuild3, `sub3.bam`):** DN0 vs DN1 scored by `bench/score_vs_guided.py` against EXPRESSED guided
loci (`hybrid/expr_c3.tsv`, u >= 3): DN1 NARROWS iff R_G(DN1) > R_G(DN0) AND P_G(DN1) >= P_G(DN0) - 0.05 (any-family);
full guided and best-overlap reported. A failed or resource-limited arm is reported as such. No default changes.

**Addendum V run note (disclosed):** the DN1 gorilla hold-out was launched twice by this session around a terminal
restart (outputs `cat_dn1_span`, started ~23:02, and `cat_span_union`, ~23:04, identical settings); the later one was
killed by PID at ~23:20 before finishing. `cat_dn1_span` is the hold-out catalog.

---
## ADDENDUM W (2026-09-13, before any number below exists) — did some multi-copy gene families arise from duplicons? (human CHM13)

**Question (user):** the conclusion sought is that some multi-copy gene families arose from duplicons. The test must be
able to say no: a classifier that also calls retrotransposition families "duplicon-borne" does not support it.

**Families (truth = CHM13 RefSeq gene/pseudogene names; regexes fixed now):**
- Literature core-duplicon families (review PMC6920530: "core or seed duplicons shared between all copies"): NPIP
  `^NPIP[AB]\d+P?$`, TBC1D3 `^TBC1D3[A-Z]?$`, GOLGA8 `^GOLGA8[A-Z]*P?\d*$`, LRRC37A `^LRRC37A\d*P?$`, RGPD `^RGPD\d+$`,
  NBPF `^NBPF\d+P?$`, SPATA31 `^SPATA31[A-Z]\d*P?\d*$`, PMS2 `^PMS2(P\d+)?$`, TRIM51 `^TRIM51[A-Z]*P?\d*$`, GUSB
  `^GUSB(P\d+)?$`.
- Negative controls, retrotransposition families (a parent gene with processed pseudogenes): GAPDH `^GAPDH(P\d+)?$`,
  PPIA `^PPIA(P\d+)?$`, EEF1A1 `^EEF1A1(P\d+)?$`.
- Reported, not in the reading: AMY (`amy_lo_v2/truth.tsv`, 12 loci).
Members on unplaced/alt contigs are dropped; families capped at 40 members + reference (name-sorted, disclosed).

**Unit and pairs:** member region = gene span extended by one gene-span length on each side (clipped to the contig).
Reference member = the protein-coding member whose gene span is closest to the family's median protein-coding span
(ties -> name); for the retro families this is the parent gene (the only protein-coding member). Pairs = reference vs
each other member. Alignment: `minimap2 -c -x asm20 -N 50 -p 0.1` (member region as query, reference region as
target); records chained as Addendum O; the pair's chain = the chain overlapping both gene spans with most aligned
bases; ALIGNED iff chain identity >= 0.80.

**Co-duplicated flank:** for an aligned pair, the chain's extension beyond each gene boundary is measured on both
members in their own coordinates, sides matched through the chain orientation. The pair is DUPLICON (co-duplicated
flank) iff on at least one matched side both members' chains extend >= 1,000 bp (`GenomeRepParams::min_block`) beyond
their gene boundary; otherwise GENE-ONLY. Reported alongside: member intron count (retro signature = member 0 introns,
reference >= 1), chain identity, extension lengths.

**Family call:** DUPLICON-BORNE iff >= 50% of its aligned pairs are DUPLICON pairs (and >= 2 aligned pairs); else NOT.
**Core (descriptive):** reference-region bases covered by the chains of >= 90% of aligned members; core length and
core length / reference gene span.

**Reading (fixed):** the test SUPPORTS "some multi-copy gene families arose from duplicons" iff >= 7 of the 10
literature core-duplicon families are DUPLICON-BORNE AND <= 1 of the 3 retrotransposition families is DUPLICON-BORNE.
If >= 2 retro families are called DUPLICON-BORNE the classifier does not discriminate and the test is uninformative.
