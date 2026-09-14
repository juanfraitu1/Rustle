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

---
## ADDENDUM X (2026-09-13, after Addendum W's NOT SUPPORTED result; before any number below exists) — genomic copy vs RNA-mediated copy, judged on HELD-OUT families

**Why:** W's flank rule scored partial gene copies inside duplicated segments as gene-only (PMS2P, GUSBP, LRRC37A*P, NBPF*P)
and 7 GOLGA8 pairs sat just under 1 kb. The biological distinction the data separate is a GENOMIC copy (shares
non-exonic sequence — introns or flanks — with its source) vs an RNA-MEDIATED retrocopy (shares only spliced exons).
The rule below was chosen after seeing W, so W's 13 families are development only; the reading uses families not
seen before.

**Rule (fixed now):** units, pairs, chaining and the identity floor exactly as W (region = gene span +/- one gene
length; reference = median-length protein-coding member; asm20 chains; aligned iff chain identity >= 0.80). For an
aligned pair, SHARED NON-EXONIC = aligned bases (M runs of the chain's CIGARs) whose member position lies outside the
member's union exons AND whose reference position lies outside the reference's union exons. Union exons = all exon
features of the gene's transcripts, else the gene's own exon features, else (no exon annotation at all) the whole gene
span counts as exonic. Pair is GENOMIC iff shared non-exonic >= 1,000 bp; else RNA-LIKE. Family DUPLICON-DERIVED iff it
has >= 1 aligned pair and >= 50% of its aligned pairs are GENOMIC. W's flank classes are also reported.

**Held-out families (CHM13 RefSeq names; regexes fixed now):**
- Segmental-duplication families (not in the core-duplicon list): NOTCH2NL `^NOTCH2NL[A-Z]?$|^NOTCH2$`, SRGAP2
  `^SRGAP2[A-D]?$`, ARHGAP11 `^ARHGAP11[AB]$`, HYDIN `^HYDIN\d?$`, FAM72 `^FAM72[A-D]$`, SMN `^SMN[12]$`, SERF1
  `^SERF1[AB]$`, ZNG1 `^ZNG1[A-F]$`, NCF1 `^NCF1[BC]?$`.
- Retrotransposition families: RPL21 `^RPL21(P\d+)?$`, HMGB1 `^HMGB1(P\d+)?$`, NPM1 `^NPM1(P\d+)?$`, RPS2
  `^RPS2(P\d+)?$`, TUBB `^TUBB(P\d+)?$`. Cap 40 members + reference as W.
**Development re-score (reported, not deciding):** W's 10 core-duplicon families, 3 retro families and AMY.

**Reading (fixed, held-out only):** SUPPORTED — "multi-copy gene families include duplicon-derived families (shared
genomic context) distinguishable from retrotransposition-derived families (exon-only homology)" — iff >= 7 of the 9
held-out segmental-duplication families are DUPLICON-DERIVED AND <= 1 of the 5 held-out retrotransposition families is.
If >= 2 retro families are DUPLICON-DERIVED the rule does not discriminate; otherwise NOT SUPPORTED.

---
## ADDENDUM Y (2026-09-14, after Addendum X was SUPPORTED; before any number below exists) — a per-copy MECHANISM layer: genomic copy, retrocopy, TE-derived gene (human CHM13)

**Why (user):** define the retrotransposition-derived copies positively instead of as "not genomic", using the repeat
annotation and the soft-masked genome as another layer. Two distinct things are named: RETROCOPIES (a gene copied through
its mRNA by L1 machinery) and TE-DERIVED genes (the gene's own coding sequence is a domesticated retrotransposon).

**Inputs (fixed):** CHM13 v2.0 soft-masked FASTA (lowercase = repeat-masked); UCSC hs1 RepeatMasker `.out`
(`winloci_data/rmsk/hs1.repeatMasker.out.gz`, md5 5cc807e8, matches UCSC md5sum.txt); `chm13v2.0_RefSeq_full.gff.gz`.

**Parent (per family):** the protein-coding member with the most annotated introns (ties: name); its mRNA = the exon
sequence of its transcript with the most exons (ties: longest). Every other member (cap 40, as W) is classified against
the parent, in this order:
1. **GENOMIC** iff the region chain (W/X regions and chaining; chain identity >= 0.80) shares >= 1,000 aligned bases that
   are non-exonic in both copies (X) AND uppercase (not soft-masked) in both copies.
2. **RETROCOPY** evidence, from the parent mRNA aligned to the member region with `minimap2 -x splice -uf -c` (record =
   the one overlapping the member gene with the most aligned query bases; identity matches/block >= 0.80; >= 100 aligned
   query bases):
   - parent junctions scored = those whose parent intron >= 70 bp and that lie >= 20 query bases inside the record;
     genomic distance between the member bases aligned to query J-10 and J+10 minus 20: <= 30 bp = LOST, >= 70 bp = RETAINED,
     otherwise ambiguous (not counted);
   - POLYA: the record reaches within 100 bp of the mRNA 3' end, and a 15-bp window with >= 12 A (T on the minus strand)
     starts inside the 60 bp from 10 bp inside to 50 bp beyond the 3' alignment end (transcript orientation);
   - TSD: an exact direct repeat >= 10 bp (no base > 60% of it) shared between genomic windows [L-40, L+10) and
     [R-10, R+40), where L / R are the left / right genomic boundaries of the insertion (alignment ends, extended to the far
     edge of the poly(A) window when POLYA holds).
   RETROCOPY iff not GENOMIC and [(LOST >= 1 and RETAINED = 0 and (POLYA or TSD)) or (no junction scorable and POLYA and TSD)].
3. **UNRESOLVED** otherwise (region chain or mRNA record passed but neither rule held). Members with neither passing are
   UNASSESSED and excluded from family fractions.

**Chance-rate gate (computed first, fixed seed 13):** 1,000 random uppercase-anchored positions on chr1-22,X, a pseudo-
insertion of 1,000 bp on a random strand; POLYA and TSD evaluated exactly as above. A hallmark whose chance rate > 10% is
dropped from the RETROCOPY rule (the other hallmark alone then applies; if both drop, only junction evidence applies and
no-junction members cannot be RETROCOPY).

**TE-DERIVED (per protein-coding member, independent of 1-3):** >= 50% of its union CDS bases (all transcripts) overlap
RepeatMasker records whose class (text before "/", "?" removed) is LINE, SINE, LTR or Retroposon.

**Family calls:** RETRO-DERIVED iff >= 50% of assessed members are RETROCOPY; GENOMIC-DERIVED iff >= 50% GENOMIC.
**Report-only:** for RETROCOPY members, whether one shares >= 1 kb uppercase flank with a sibling RETROCOPY (asm20,
identity >= 0.80) = "retrocopy, then genomic duplication".

**Families (regexes fixed now; names checked to exist, no alignments run):**
- RETRO positives (held out, 16): processed-pseudogene families RPL7 `^RPL7(P\d+)?$`, RPL23A `^RPL23A(P\d+)?$`, HNRNPA1
  `^HNRNPA1(P\d+)?$`, PTMA `^PTMA(P\d+)?$`, TPT1 `^TPT1(P\d+)?$`, NACA `^NACA(P\d+)?$`, KRT8 `^KRT8(P\d+)?$`, FTH1
  `^FTH1(P\d+)?$`; retrogene pairs PGK `^PGK[12]$`, GLUD `^GLUD[12]$`, UTP14 `^UTP14[AC]$`, TAF1L `^TAF1L?$`, RPL10L
  `^RPL10L?$`, PABPC `^PABPC[13]$`, POU5F1 `^POU5F1B?$`, NANOG `^NANOG(P8)?$`.
- SD negatives (13): Addendum X's 9 SD families (seen under X, new to this layer) + new DEFB4 `^DEFB4[AB]$`, GTF2H2
  `^GTF2H2C?$`, SPDYE `^SPDYE\d+[A-Z]?$`, USP17L `^USP17L\d+$`.
- TE-derived groups: ERV-env `^ERVW-1$|^ERVFRD-1$|^ERVV-[12]$|^ERVH48-1$|^ERVMER34-1$|^ERVK3-1$|^ERV3-1$` (deciding);
  Ty3/gypsy-derived `^PEG10$|^RTL\d+[A-Z]?$|^ARC$|^ASPRV1$|^NYNRIN$`, PNMA `^PNMA\d+[A-Z]?$`, L1TD1 (reported; expected
  LOW before looking — ancient domestications are too diverged for RepeatMasker).
- Development (reported, not deciding): W's 10 core + 3 retro families, X's 5 retro families, AMY.

**Readings (held out, each fixed):**
- R1 RETROCOPY: SUPPORTED iff >= 12/16 RETRO positives RETRO-DERIVED AND <= 1/13 SD negatives RETRO-DERIVED.
- R2 GENOMIC with repeats excluded: SUPPORTED iff >= 10/13 SD negatives GENOMIC-DERIVED AND <= 1/16 RETRO positives.
- R3 TE-DERIVED: SUPPORTED iff >= 6/8 ERV-env genes TE-DERIVED AND <= 5% of protein-coding members of all SD, core and
  retro families (dev + held out) TE-DERIVED.
Each reading is reported separately; a failed reading is reported as failed, not re-tuned on these families.

---
## ADDENDUM Z (2026-09-14, after Addendum Y's R1 failed; before any number below exists) — retrocopy = parent-intron loss, on FRESH held-out families

**Why:** Y's RETROCOPY rule failed on its hallmark requirement (poly(A) gated on reaching the annotated 3' end; hallmarks
decay with age), while intron loss alone looked near-perfect POST HOC (§6jv). That observation was made on Y's families,
so every family below is new to W, X and Y; Y's 16 retro + 13 SD families become development.

**Rule (fixed now):** everything as Addendum Y (parent, regions, GENOMIC with soft-masked bases excluded, parent mRNA,
`minimap2 -x splice -uf`, record selection, mRNA identity >= 0.80 and >= 100 aligned query bases, junction scoring:
parent intron >= 70 bp, junction >= 20 query bases inside the record, <= 30 bp = LOST, >= 70 bp = RETAINED), except:
**RETROCOPY iff not GENOMIC AND the mRNA record passes AND LOST >= 1 AND RETAINED = 0.** No scorable junction
(including every intronless parent) -> UNRESOLVED; intronless parents are outside this rule's scope. Poly(A) (Y's
definition, and also evaluated without the 3'-reach gate) and TSD are REPORTED as age/confidence annotations only.
Family RETRO-DERIVED iff >= 50% of assessed members are RETROCOPY; GENOMIC-DERIVED iff >= 50% GENOMIC.

**Fresh families (regexes fixed now; names checked to exist, no alignments run):**
- Retro positives (16): processed-pseudogene families RPL5 `^RPL5(P\d+)?$`, RPS3A `^RPS3A(P\d+)?$`, RPL31
  `^RPL31(P\d+)?$`, EEF1B2 `^EEF1B2(P\d+)?$`, HMGN2 `^HMGN2(P\d+)?$`, YBX1 `^YBX1(P\d+)?$`, CYCS `^CYCS(P\d+)?$`,
  KRT18 `^KRT18(P\d+)?$`; literature retrogenes GK `^GK2?$`, CETN `^CETN[12]$`, NAP1L `^NAP1L[123]$`, MKRN
  `^MKRN[13]$`, PDHA `^PDHA[12]$`, CSTF2 `^CSTF2T?$`, UBL4 `^UBL4[AB]$`, FAM50 `^FAM50[AB]$`.
- SD negatives (10): GTF2IRD2 `^GTF2IRD2B?$`, PRAMEF `^PRAMEF\d+$`, SPANX `^SPANX[A-D]\d?$`, FCGR3 `^FCGR3[AB]$`,
  CFHR `^CFHR[1-5]$`, RH `^RHD$|^RHCE$`, HP `^HPR?$`, CYP2D `^CYP2D[67]$`, OPN1 `^OPN1[LM]W\d?$`, CGB `^CGB\d+$`.
- Development (reported): Y's 16 retro positives and 13 SD negatives under this rule.

**Reading (fixed):** SUPPORTED iff >= 12/16 fresh retro positives RETRO-DERIVED AND <= 1/10 fresh SD negatives
RETRO-DERIVED. Also reported, not deciding: GENOMIC-DERIVED counts (expected SD high, retro 0) and member-level
junction-loss among SD members. TE-DERIVED is not re-tested (reported as CDS repeat coverage only).

---
## ADDENDUM AA (2026-09-14; before any split tree below is computed) — subfamilies from variation-graph bubbles: a threshold-free split-dominance tree vs IQ-TREE

**Why (user):** merge the graph view with the tree. In a family variation graph every bubble splits the copies (paths)
in two; a pairwise-compatible split set is exactly one tree (Buneman 1971; splits-equivalence theorem), and incompatible
splits are conversion/mosaic/homoplasy signal.

**Inputs (fixed, already on disk, produced by the §6js guided run):** the 66 reference-projected alignments
`lit/guided_t/tree_t/*.proj.fa` (NPIP, TBC1D3: ref + half x5 + keep1 x5, exon and intron) and `lit/amy_lo_t/tree_t/*.proj.fa`
(AMY, same design), their IQ-TREE `.treefile`s, `t.candidates.tsv` and `truth.tsv`. Seen before registering: only
column-completeness counts (e.g. ref NPIP exon has 0 gap-free columns, ref NPIP intron 162 gap-free informative columns);
no split or tree was computed.

**Rule DT (fixed):**
1. Bubble = alignment column with no gap/N in any leaf and exactly two nucleotide states, each carried by >= 2 leaves
   (parsimony-informative biallelic). Columns with gaps, N, >2 states or a singleton state are not used (gaps are not
   indels here: in a projected alignment they mix deletion and missing sequence).
2. Split of a bubble = the two leaf sets; support(S) = number of bubbles inducing S.
3. S and T are incompatible iff all four intersections of their sides are non-empty.
4. **KEEP S iff support(S) > support(T) for every split T incompatible with S.** (Kept splits are pairwise compatible —
   two kept incompatible splits would each need more support than the other — hence display as one unrooted tree.)
   No identity cut, no support threshold, no substitution model.
5. A literature group is RECOVERED by DT iff a kept split restricted to truth-labelled leaves equals the group or its
   complement (same test as `clade_calls`); RECOVERED by IQ-TREE iff the same holds for a split with SH-aLRT > 75.
6. Fewer than 4 leaves or zero bubbles -> not treed (group counts as not recovered for that class).

**Labels:** reference runs use leaf names; leave-out runs label seeds by name and `cand{i}` by its hidden-truth match
from `t.candidates.tsv` (i = row order within that level/rep, class hidden, match in the same family).
**Gate (must pass before any DT number is read):** re-deriving the IQ-TREE calls from the treefiles with these labels
reproduces every per-run call printed in `guided_t/t.out` and `amy_lo_t/t.out`; otherwise stop and report.

**Reported per run and summed:** literature groups recovered (exon, intron, either), TBC1D3 positional split (kept or
not), number of kept non-trivial splits vs IQ-TREE supported splits, LITERATURE-CONFLICTING splits (kept DT / supported
IQ splits incompatible, on truth-labelled leaves, with any present literature group), and the CONFLICT INDEX (fraction of
bubbles whose split is incompatible with >= 1 kept split).

**Reading (fixed):** DT is SUPPORTED as the subfamily definition iff, summed over all runs and literature groups, DT
"either"-class recoveries >= IQ-TREE "either"-class recoveries AND DT literature-conflicting splits <= IQ-TREE
literature-conflicting splits. Otherwise NOT SUPPORTED. Per-family and per-class differences are reported, not deciding.

---
## ADDENDUM AB (2026-09-14; before any number below exists) — the shared definition instantiated DE NOVO, gorilla hold-out; and the (a)/(b)/(c) gap decomposition

**Why (user):** one definition for both modes and both levels (`docs/seeded_family_definition.md` §0★★). De novo = read
loci consolidated to gene level, then the SAME edge rule as guided; families = connected components.

**Substrate (held out; never used to build the guided rule):** gorilla `rebuild3/sub3.bam` contigs NC_073241.2,
NC_073242.2, NC_073244.2; genome `GGO.fasta` restricted to those contigs. Nodes come from the DN0 run's own node dump
(`rebuild3/dump_default/ggo.nodes.tsv`: 2,469 read-supported locus representatives with exons; representative sequence =
genome exon sequence in transcript orientation — checked, no alignment run).

**Arms (fixed):**
- **AB1 (edge rule only):** nodes = the 2,469 dump nodes unchanged.
- **AB2 (full definition, PRIMARY):** nodes consolidated to gene level: (i) cut each node's exon chain at every intron
  longer than 271,359 bp (P99.9 of 1,092,233 annotated GGO_genomic.gff intron lengths, computed before any alignment);
  pieces with exon sum < 100 bp are dropped; (ii) same-strand pieces whose exon blocks overlap by >= 1 bp are merged
  (connected components) into one locus: exons = union, n_reads = sum, representative = the member piece with most reads
  (ties: longest exon sum).
**Edges (both arms, identical to guided `bench/guided_pipeline.py` finders):**
- EXON edge u->v: u's representative transcript, `minimap2 -x splice -uf -N 50 -p 0.1` onto the three contigs; a hit with
  identity (matches/block) >= 0.80 and query coverage >= 0.50 (`transcript_hits`) whose aligned exon blocks overlap v's
  exon blocks by >= 1 bp.
- GENE-BODY edge u->v: u's gene body (first exon start .. last exon end, genome forward), `minimap2 -c -x asm20 -N 50
  -p 0.1` onto the three contigs, chained exactly as `gene_body_chains` (identity >= 0.80, aligned >= 0.50 of
  min(query, extrapolated target span)); a chain whose target interval overlaps >= 1 exon block of v.
- Both: hits/chains overlapping u's own span are ignored; when u and v both have >= 2 exons, the implied transcript
  orientation on v must equal v's strand. Edges are symmetrised.
**Families:** connected components with >= 2 distinct loci; emitted as `copies.tsv` (family_id, chrom, start, end, ...).

**Primary reading (fixed):** `bench/score_vs_guided.py --contigs include:NC_073241.2,NC_073242.2,NC_073244.2 --expr
hybrid/expr_c3.tsv`, any-family. AB2 SUPPORTED iff R_G > 0.4710 AND P_G >= 0.3386 (i.e. beats the shipped opt-in LCS
union's recall, R_G 0.4710, with precision no more than 0.05 below its 0.3886). AB1 is read with the same bar,
reported, not deciding. Also reported: best-overlap metrics, full-guided metrics, DN0 / union / DN1 rows.

**Gap decomposition (reported, not deciding), for DN0, LCS union, DN1, AB1, AB2:** universe = the 741 co-clustered
pairs of expressed guided loci. Guided edges between loci on the three contigs are re-derived from
`mcl_ann/allgenes_gw.asm20.paf` at the guided catalog's floors (identity >= 0.70, >= 300 aligned bp, aligned >= 0.30 of the
longer interval; gene intervals mapped to the locus they overlap most). Each pair is RECOVERED if the catalog puts the
two loci in a shared family (any-family); otherwise:
(c) UNEXPRESSED BRIDGE if the two loci are in different components of the guided edge graph induced on the EXPRESSED loci
of their cluster; else (a) MISSING NODE if either locus overlaps no node of that catalog's node set (DN0/DN1/AB1: the DN0
dump; union: its own dump; AB2: consolidated loci); else (b) MISSING EDGE. Counts and fractions of 741 are reported.

---
## ADDENDUM AC (2026-09-14; before any number below exists) — the two levers from §6jz: seed-style grouping and missing nodes

**Seen before registering (diagnosis only, no arm computed):** 58 of the 304 expressed guided loci have no DN0 node; all 58
carry >= 3 primary MAPQ >= 1 reads and 54 carry >= 3 spliced reads (median locus span 27 kb) — the reads exist, the
pipeline dropped the loci.

**Lever G — leader neighbourhoods instead of components (the guided rule without annotation):** nodes are ordered by
n_reads (desc), then edge degree (desc), then chromosome and start. Walking that order, an unassigned node with >= 1
unassigned neighbour becomes a LEADER and its family = the leader + its unassigned direct neighbours; those are marked
assigned. A node joins a family only as a direct neighbour of that family's leader (no chaining through members). Families
with >= 2 loci are emitted. Edges are exactly Addendum AB's.

**Lever N — read-locus nodes for loci the pipeline dropped:** from `sub3.bam` on the three contigs, primary mapped
reads with MAPQ >= 1 (the expression gate's reads); exon blocks = CIGAR segments between N operations; transcript strand
= read orientation, flipped when `ts:A:-`. Same-strand reads whose exon blocks overlap by >= 1 bp are grouped (connected
components); a group with >= 3 reads becomes a candidate locus whose exons are the bases covered by >= 2 reads' blocks
(merged; exon sum >= 100 bp required). A candidate is ADDED as a node iff its exons overlap no AB2 node's exons (either
strand). Its representative transcript is its exon sequence and its body runs from first to last exon; edges are
computed for the added nodes with Addendum AB's alignments and rules (queries not yet aligned are aligned the same way).
No annotation or guided interval is used to build nodes.

**Arms:** AC1 = AB2 nodes + lever G; AC2 = AB2 nodes + lever N + components; **AC3 (PRIMARY) = AB2 nodes + lever N +
lever G.**
**Reading (fixed):** `bench/score_vs_guided.py` on expressed guided loci (as AB), any-family. AC3 SUPPORTED iff R_G >
0.4710 AND P_G >= 0.3386. AC1 and AC2 are read with the same bar and reported, not deciding. Also reported: best-overlap
and full-guided metrics, number of added nodes, how many of the 58 node-less expressed guided loci gain a node, and the
gap decomposition (missing-node-first order) for AC1-AC3.

---
## ADDENDUM AD (2026-09-14; before any number below exists) — (1) is nucleotide unreachability of old copies synonymous-codon saturation? (2) intermediate grouping; (3) truth granularity

### AD-1 Codon divergence of unreachable retrogenes (human CHM13)
**Question (user):** are the copies that nucleotide alignment cannot reach (Addendum Z UNASSESSED) diverged at synonymous
codon positions while the protein is conserved — and when two of three codon bases differ, can degeneracy still keep the
amino acid?
**Pairs (parent = protein-coding member with most introns):** UNREACHED (Z unassessed): GK-GK2, CETN2-CETN1,
NAP1L1-NAP1L2, NAP1L1-NAP1L3, CSTF2-CSTF2T, UBL4A-UBL4B, FAM50A-FAM50B (7). REACHED (nucleotide mRNA record passed in Y/Z):
MKRN1-MKRN3, PDHA1-PDHA2, UTP14A-UTP14C, PABPC1-PABPC3, POU5F1-POU5F1B, PGK1-PGK2, GLUD1-GLUD2, TAF1-TAF1L, RPL10-RPL10L (9).
**Method (fixed):** CDS = the transcript with the longest CDS (RefSeq CDS features, CHM13 soft-masked genome, upper-cased);
translate (standard code); protein pair aligned with MAFFT (`--auto`); codon alignment by back-translation; over codons
aligned in both: protein identity, identity at codon positions 1, 2, 3, CDS nucleotide identity, fraction of codons
differing at 0/1/2/3 positions and the synonymous share of 2- and 3-difference codons; dN and dS by PAML `yn00`
(Yang-Nielsen). Protein-level reach: `miniprot` of the parent protein onto the member region (gene +/- one gene length);
REACHED-BY-PROTEIN iff the best alignment covers >= 50% of the parent protein; PROTEIN INTRON LOSS iff that alignment has
no intron and its protein span covers >= 1 parent CDS junction.
**Reading (fixed):** SYNONYMOUS SATURATION explains nucleotide unreachability iff in >= 5/7 UNREACHED pairs: identity at
position 3 < both positions 1 and 2, dS >= 1.0, and protein identity >= 0.70. PROTEIN LEVEL RESCUES iff >= 6/7 UNREACHED
pairs are REACHED-BY-PROTEIN with PROTEIN INTRON LOSS. The REACHED pairs are reported as the contrast; the synonymous share
of multi-difference codons is reported for all pairs.

### AD-2 Intermediate grouping on AC nodes and edges (gorilla hold-out)
- **AD2a (PRIMARY): 2-edge-connected components** — remove every bridge (an edge whose removal disconnects its component),
  then take components with >= 2 loci. No parameter; every emitted family has edge connectivity lambda >= 2.
- **AD2b: triangle-supported leaders** — leaders as Addendum AC; the family = leader + free direct neighbours + free nodes
  adjacent to >= 2 nodes of that star; one pass per leader.
**Reading (fixed):** as AB/AC: SUPPORTED iff any-family R_G > 0.4710 AND P_G >= 0.3386 on expressed guided loci. AD2a
decides; AD2b is read with the same bar and reported. Gap decomposition (missing-node first) reported.

### AD-3 Truth granularity (reported, not deciding)
A COARSE truth merges guided MCL clusters whose loci (expressed or not, on the three contigs) are joined by >= 1 guided
edge re-derived as in Addendum AB. DN0, LCS union, AB2, AC1-AC3, AD2a, AD2b are scored against it (expressed loci, same
scorer). Reported: coarse cluster count and pairs, and each catalog's R_G/P_G under both truths.

---
## ADDENDUM AE (2026-09-14; before any number below exists) — (A) confirm triangle-supported leaders on a fresh gorilla substrate; (B) retrocopy rule with a coding-sequence, gap-excluded identity and a protein-level fallback, on fresh human families

### AE-A Fresh substrate for the grouping clause
**Substrate (fixed now, chosen by guided-locus and read counts only, before any de novo run):** gorilla contigs
NC_073233.2, NC_073240.2, NC_073238.2 (1,345 `gw_units_v3` loci; 429,036 reads in `GGO_ds.bam`), none used in
Addenda AB-AD. `subF.bam` = `samtools view -b GGO_ds.bam <3 contigs>` (the sub3.bam recipe).
**Catalogs:** DN0 and LCS union built with the same binary as the V/AB work (`lit/dn_v/gw_family_catalog.54154909
--homology-primary --threads 4`; union = `RUSTLE_ER_UNION_LCS=1`), each with `RUSTLE_ER_EDGE_DUMP`. Expression of the
guided loci on `subF.bam` with `bench/interval_expression.py`; expressed = u >= 3.
**Definition arms (code unchanged from AB/AC/AD except the contig list):** AB2 consolidation of the DN0 dump nodes, AC
read-locus nodes (same rule), AB edges; groupings: components (AC2), leaders (AC3), bridge-split (AD2a), **triangle-
supported leaders (AD2b, DECIDING)**.
**Reading (fixed, relative to the union on the same contigs):** triangle-supported leaders CONFIRMED iff any-family
R_G(triangle) > R_G(union) AND P_G(triangle) >= P_G(union) - 0.05, on expressed guided loci (u >= 3, clusters >= 2
loci). Other arms, best-overlap and full-guided metrics, and the gap decomposition are reported.

### AE-B Retrocopy rule, fixed metric (human CHM13)
**Rule (fixed):** as Addendum Z (parent, regions, GENOMIC first with soft-masked bases excluded) except:
1. NUCLEOTIDE LEVEL: query = the parent's coding sequence (longest-CDS transcript, CDS segments only), `minimap2 -x splice
   -uf -N 50 -p 0.1` onto the member region; record = the one overlapping the member gene with most aligned query bases;
   identity = gap-excluded (aligned M bases minus mismatches, over aligned M bases; mismatches = NM minus inserted and
   deleted bases); passes iff identity >= 0.80 and >= 100 aligned query bases. Junctions = CDS segment boundaries; LOST /
   RETAINED exactly as Z.
2. PROTEIN LEVEL, only when the nucleotide record does not pass: `miniprot` of the parent protein onto the member region;
   passes iff the best alignment covers >= 50% of the parent protein. A parent CDS junction j (amino-acid coordinate)
   lying >= 7 aa inside the alignment is RETAINED iff an intron operation (N/U/V) starts within 5 aa of j, else LOST.
3. RETROCOPY iff not GENOMIC and the first passing level has LOST >= 1 and RETAINED = 0; UNRESOLVED if a level passes
   otherwise; UNASSESSED if neither passes. Family RETRO-DERIVED iff >= 50% of assessed members are RETROCOPY.
**Fresh families (regexes fixed now; names checked, no alignments run):**
- Retro positives (13): processed-pseudogene families RPL13 `^RPL13(P\d+)?$`, RPS6 `^RPS6(P\d+)?$`, RPL35A
  `^RPL35A(P\d+)?$`, EEF1G `^EEF1G(P\d+)?$`, ACTG1 `^ACTG1(P\d+)?$`, HSPA8 `^HSPA8(P\d+)?$`, RPL12 `^RPL12(P\d+)?$`;
  literature retrogenes SET `^SET$|^SETSIP$`, ELOA `^ELOA2?$`, RBMX `^RBMX(L[123])?$`, PPP1R2 `^PPP1R2B?$`, CDY
  `^CDYL$|^CDY\d[AB]?$`, TAF7 `^TAF7L?$`.
- SD negatives (8): PSG `^PSG\d+$`, CT45A `^CT45A\d+$`, GAGE `^GAGE\d+[A-Z]?$`, XAGE1 `^XAGE1[AB]$`, TSPY `^TSPY\d+$`,
  CSAG `^CSAG[123][A-C]?$`, CCL4L `^CCL4L\d$|^CCL4$`, POTE `^POTE[A-M]$`.
- Development (reported): Addendum Z's fresh families (16 retro, 10 SD) under this rule.
**Reading (fixed):** SUPPORTED iff >= 10/13 fresh retro positives RETRO-DERIVED AND <= 1/8 fresh SD negatives
RETRO-DERIVED. Reported: which level assessed each member, and the Z retrogenes that were UNASSESSED (GK2, CETN1, NAP1L2/3,
CSTF2T, UBL4B, FAM50B) under this rule.

---
## ADDENDUM AF (2026-09-14; before any number below exists) — (1) Rust port gate; (2) human substrate; (3) missing nodes

### AF-1 Rust port (engineering gate, not a scientific test)
`RUSTLE_SHARED_DEFINITION=1` in `detect_homology_catalog_genome_wide` replaces the E_r / γ-QC / coverage-split / distinct-
locus stage with the Python prototype's construction (`bench/denovo_shared_def.py`): reps -> AB2 consolidation -> AC read-
locus nodes (from the same MAPQ >= 1 primary reads) -> exon + gene-body edges (minimap2 `-x splice -uf` / `-x asm20`,
`-N 50 -p 0.1`, onto the BAM's contigs) -> triangle-supported leaders -> families (>= 2 loci). Unset = byte-identical.
**Gate:** (i) all library tests pass; (ii) on both gorilla substrates the opt-in Rust catalog's family memberships (as
sets of copy coordinates) equal the Python `ad2b.copies.tsv` exactly; (iii) a default Rust run on the fresh substrate
reproduces `catF_default.copies.tsv` byte-for-byte. A failed gate is reported with the differing families, not tuned.

### AF-2 Human substrate (third substrate; never pooled with gorilla)
`human2/sub_h3.bam` (A119b_ds, chr15/chr17/chr22), CHM13 v2.0, guided truth `human2/guided.clusters.tsv` (RefSeq
gene+pseudogene MCL, built 09-13 before any de novo human catalog). DN0 and LCS union built with binary 54154909 as for
gorilla; expression of guided loci on `sub_h3.bam` with `interval_expression.py`, expressed = u >= 3. Arms: shared
definition with AC read-locus nodes and triangle-supported leaders (PRIMARY), components and leaders reported.
**Reading (fixed):** CONFIRMED on human iff R_G(triangle) > R_G(union) AND P_G(triangle) >= P_G(union) - 0.05 (expressed
guided, any-family).

### AF-3 Missing nodes — split chained read groups
**Diagnosis seen before registering (development substrates):** of the expressed guided loci still without a node after
AC, 14/15 (substrate 1) and 16/31 (substrate 2) belong to a read group whose depth >= 2 exons overlap an existing node
elsewhere (group span median 133 / 76 kb), so AC rejected the whole group; 1/15 and 13/31 sit in groups of < 3 reads.
**Rule (fixed):** within each AC read group, depth >= 2 exon segments are computed as in AC; two segments are LINKED iff
>= 2 reads of the group have exon blocks overlapping both; each connected component of linked segments is a sub-locus
whose supporting reads are the group's reads overlapping any of its segments; a sub-locus with >= 3 supporting reads and
>= 100 bp is added as a node iff none of its segments overlaps an AB2 node's exons (n_reads = supporting reads). Edges
and triangle grouping unchanged.
**Reading (fixed, on human — AF-3 arm vs the AF-2 primary arm):** SUPPORTED iff missing-node pairs decrease AND R_G(split)
>= R_G(AC) AND P_G(split) >= P_G(AC) - 0.05. The two gorilla substrates are reported as development.

---
## ADDENDUM AG (2026-09-14; before any catalog is scored against the truth below) — the RNA-level ground truth of clause 6, and the goal bars

**Why (user goal):** a definition with very high sensitivity, precision and bipartite matching against ground truth.
Clause 6 of the definition (written in §0★★ before any de novo scoring) says the RNA-level family is the DNA-level
family graph INDUCED on expressed copies, and that a pair joined only through unexpressed copies is a level
difference, not a de novo error. Until now de novo catalogs were scored against DNA-level clusters restricted to
expressed loci, which counts those pairs as misses. The earlier attempt to derive the induced graph (§6jz) was
unreliable because its edges used span denominators.

**Guided graph (fixed):** `mcl_families --paf allgenes_gw.asm20.paf --gff GGO_genomic.gff --min-exonic-bp 1 --dump-graph`
reproduces the `gw_units_v3` graph exactly (13,263 nodes, 64,336 edges = its params certificate). Human: the same
command on `human2/genes.asm20.paf` with `human2/refseq_c15_17_22.gff`, checked against `human2/guided.clusters.tsv`'s
params certificate (nodes/edges) before use.
**RNA-level truth (fixed):** within each guided cluster (the scorer's contig set), take its expressed loci (u >= 3);
map each locus to the graph node(s) overlapping it on the same contig (a cluster locus may fold several annotation
records; any overlapping node counts); two expressed loci are joined iff connected by a path of graph edges whose
every node maps to an expressed locus of that cluster; RNA truth clusters = connected components (>= 2 loci).
**Metrics (fixed, both truths reported):** best-overlap assignment of each truth locus to a catalog family;
pairwise sensitivity and precision; bipartite micro recall and precision (`guided_pipeline.bipartite`), and their F.
**Goal bars (declared now, to be met on a HELD-OUT substrate before claiming them):** pairwise sensitivity >= 0.90,
pairwise precision >= 0.90, bipartite F >= 0.90, on the RNA-level truth. A rule change is adopted only after passing its
own pre-registered bar on a substrate not used to design it.

---
## ADDENDUM W (2026-09-14, after §6jt; before any simulated read exists) — can de novo find all amylase loci when every locus is expressed? (simulated IsoSeq)

**Why (user):** in real testis IsoSeq the amylase window was sparse (3,205 primary reads) and de novo missed loci (DN0 9/12,
DN1 11/12). Simulation removes expression as the limit and asks whether node construction + E_r finds every copy.
**Truth:** AMY v2 (12 loci, `lit/amy_lo_v2/truth.tsv`).
**Source transcripts (fixed now):** every exon-bearing RefSeq transcript of the 12 loci (spliced, transcript orientation),
except LOC124905662, whose model is its exons 2-8 (v2 correction, §6jr); AMYP1 (no model) = the exon blocks of the best
`minimap2 -x splice` alignment of AMY2A NM_000699.4 inside the AMYP1 span with identity >= 0.80, else its unspliced span.
**Reads:** `bench/sim_reads.simulate_reads` (substitution 0.003, indel 0.0008, end-truncation up to 30%, deterministic
seeds), reads >= 300 bp; depth PRIMARY 40 reads per locus split evenly over its transcripts, SECONDARY 10 per locus.
Named `SIMAMY|<locus>|<transcript>|<i>`. Aligned `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes` to
CHM13 (prebuilt k15/w10 index, disclosed), sorted and indexed.
**Arms:** DN0 `gw_family_catalog --homology-primary` and DN1 = DN0 + `RUSTLE_ER_UNION_GENOMIC_SPAN=1`, at both depths.
**Scores:** `bench/denovo_subfamilies.py` (presence, families and copies per record, collapse, family-level pairwise
sensitivity/precision and bipartite, exon/intron clades AMY1, AMY2, AMY2Ap); read level: fraction of primary alignments
on the source locus and MAPQ distribution per locus type.
**Readings (fixed):** an arm FINDS ALL AMYLASE LOCI iff 12/12 present, 0 collapsed, and family-level pairwise sensitivity
and precision both 1.000 (all 12 in one family, nothing else in it); otherwise the missing/split loci are listed with
their read-level placement. Descriptive otherwise.

---
## ADDENDUM AH (2026-09-14; before any held-out number below exists) — the node rule and guided minimal annotation on held-out substrates

**Code frozen at commit aa869a08.** Development (substrate 1, ledger §6kg-§6ki) chose:
**(D) de novo definition** = `bench/read_gene_nodes.py` nodes (strong 5'/3' end-site clipping, W = 50, FRAC = 0.5, >= 2-read
linkage split) + the guided catalog's construction (`bench/node_graph_mcl.py`: all-vs-all `minimap2 -x asm20 -c -X -N 50
-p 0.1` of node spans, chunks <= 3 Mb; `mcl_families --min-exonic-bp 1`; for human the human guided catalog's params,
i.e. `--min-exonic-bp 0`).
**(G) guided minimal annotation** = `bench/guided_min.py`: a random 50% of annotated genes on the substrate (seed 1) + candidates
from seed gene-body chains (identity >= 0.70, aligned >= 0.30 of min(query, extrapolated span), 1 round) + the same
construction (seed-seed pairs from the annotation PAF; human `--min-exonic-bp 0`).

**Held-out substrates:** gorilla s3 = NC_073230.2 + NC_073228.2 (never used in AB-AG or development); human chr15/17/22 (used
for AF-2 triangle confirmation only; never for the node rule or the guided minimal design).
**Truths:** (D) the RNA-level truth of Addendum AG (`lit/truth_rna/rna_s3.clusters.tsv`, `rna_h.clusters.tsv`); (G) the DNA-level
guided clusters on the substrate (all loci in clusters >= 2).
**Comparators:** (D) triangle-supported leaders (the confirmed shared definition; gorilla s3 from the Rust opt-in binary
`gw_family_catalog.sd_wip`, human from `lit/sharedef_H/tri.copies.tsv`); (G) the 50% seeds without candidates, and all
annotated genes (ceiling).
**Metrics:** `bench/rna_truth.py score` — pairwise sensitivity/precision and bipartite R/P/F.

**Readings (fixed, each substrate separately, never pooled):**
- D SUPPORTED iff bipartite F(D) > bipartite F(triangle) AND pairwise precision(D) >= pairwise precision(triangle) - 0.05.
- G SUPPORTED iff bipartite F(G) >= bipartite F(seeds only) + 0.10 AND pairwise precision(G) >= pairwise precision(seeds only) - 0.05.
- GOAL bars (Addendum AG), reported for every arm: pairwise sensitivity >= 0.90, pairwise precision >= 0.90, bipartite F >= 0.90.

---
## ADDENDUM X (2026-09-14, after §6ju; before any number below exists) — provenance-typed de novo catalog: typed edges, RNA-typed boundaries, presence categories

**Why (user):** separate DNA and RNA evidence. Feasibility check done before this addendum (no scoring): the three
readthroughs that cost loci (RNPC3->AMY2B, LOC124905668->AMYP1, TBC1D3->NPEPPSP1) contain NO pair of SD-paralogous atoms
— duplication blocks are larger than genes, so annotation-free DNA structure cannot place gene boundaries; boundaries
stay RNA evidence (de novo) or annotation (guided). User chose: typed edges, RNA-typed boundaries, presence categories.

**Edge types (provenance).** TX = E_r edges on the rep's assembled exon sequence (production substrate; `<dump>.edges.tsv`);
SPAN = E_r edges on the rep's genomic span (`<dump>.call2.edges.tsv`, same floors); pairs are typed TX-only, SPAN-only,
or BOTH. **Restricted DNA admission (new, Rust, opt-in `RUSTLE_ER_UNION_GENOMIC_SPAN=restricted`):** TX edges are
partitioned first (same gamma_quasi_clique_partition, same gamma); a SPAN-only pair is ADMITTED unless its two reps lie
in two different TX blocks that both have >= 2 reps (DNA may attach an RNA-isolated rep to a family, never merge two
RNA-established families); then the usual partition runs on TX ∪ admitted SPAN. `=1` keeps the full union (§6jt).
**RNA-typed boundary:** `RUSTLE_LOCUS_BRIDGE_CUT=1` (§6ji), reported as an RNA-evidence arm.

**Arms:** DN0 (default), DN1 (full union), DN1r (restricted), DN1r+B (restricted + bridge cut), DN0+B (bridge only), each
with `RUSTLE_ER_EDGE_DUMP`. Development: human lit windows (`lit.bam`) and the amylase window (`amy.bam`). Hold-out:
gorilla rebuild3 `sub3.bam` — DN1r first, then DN1r+B (each ~3.5 h; a run that fails or is killed is reported as not
run); DN0 and DN1 hold-out catalogs exist (§6jt); DN0+B = §6ji's `cat_bridge`.

**Scores (`bench/provenance_eval.py`):**
1. Per-type edge precision: an edge between reps whose best-overlapping truth loci are both known is TRUE if those loci
   share a truth family (human windows: truth family; gorilla: guided gw_units_v3 cluster), FALSE otherwise; edges with
   an endpoint on no truth locus are unscored. Reported for TX-only, SPAN-only admitted, SPAN-only rejected, BOTH.
2. Presence categories per truth locus (human truth tables; gorilla: guided loci on the 3 contigs): R = an RNA locus of its
   own (a rep whose span lies >= 50% inside the locus); T = readthrough-only (no R rep, but >= 3 primary MAPQ>=1 reads
   with aligned bases in the locus — `interval_expression.py` counts; for gorilla `hybrid/expr_c3.tsv`); N = < 3 such
   reads. Reported: counts per category and the fraction of each category in an emitted family, per arm.
3. Family level and clades as §6jt (`denovo_subfamilies.py`) on the human windows.
**Hold-out decision (fixed):** DN1r (and separately DN1r+B) NARROWS the de novo <-> expressed-guided gap iff R_G > DN0's
0.3887 AND P_G >= 0.3345 - 0.05 = 0.2845 (`score_vs_guided.py --expr hybrid/expr_c3.tsv`, any-family). Expected in
advance: the §6jt precision loss sits in SPAN-only edges that merge two TX blocks (the ones restricted mode rejects).

---
## ADDENDUM AI (2026-09-14; before any number below exists) — is the guided definition robust to the annotation it starts from?

**Why:** the guided definition reproduces its truth on the same annotation by construction (tautological). A
non-tautological test gives it an INDEPENDENT annotation of the same genome and scores against the RefSeq-derived truth.
**Annotation (fixed):** T2T consortium CAT/Liftoff GENCODE-based CHM13 annotation
`chm13.draft_v2.0.gene_annotation.gff3` (downloaded 09-14, md5 prefix 7e946417), every `gene` on
chr15/chr17/chr22 (7,388); exon union per gene from its transcripts' exons.
**Construction (fixed, identical to the human RefSeq guided catalog):** gene spans all-vs-all `minimap2 -x asm20 -c -X -N 50
-p 0.1` (3 Mb chunks, `bench/node_graph_mcl.py`), `mcl_families` defaults (identity >= 0.70, cov_longer >= 0.30 on
exonic length, >= 300 bp, MCL I = 2.8, `--min-exonic-bp 0`).
**Truth:** `human2/guided.clusters.tsv` (RefSeq gene+pseudogene MCL), DNA level (all loci in clusters >= 2); RNA level
(`lit/truth_rna/rna_h.clusters.tsv`) reported.
**Reading (fixed):** the guided definition is ANNOTATION-ROBUST iff, on the DNA-level truth, pairwise sensitivity >= 0.90,
pairwise precision >= 0.90 and bipartite F >= 0.90 (the goal bars). Otherwise NOT, with the failing metric reported.

---
## ADDENDUM AJ (2026-09-14; after AI's NOT, before any number below exists) — the exon-to-exon edge clause on both annotations, held out

**AI result that motivates this (fixed, already reported):** GENCODE/CAT vs RefSeq truth on chr15/17/22 = pairwise
sensitivity 0.461, pairwise precision 0.411, bipartite F 0.666 — NOT annotation-robust. Diagnosis (post hoc, looked at):
71% of false pairs sit in ONE GENCODE family (266 genes: 114 protein-coding, 102 lncRNA), and the RefSeq TRUTH's own
largest cluster has the same shape (202 loci: 96 protein-coding, 89 lncRNA; MCL1 joins TBC1D3 with KRT). Both were built
at `--min-exonic-bp 0`, i.e. WITHOUT §6dt's exon-to-exon clause, which the gorilla truth (`gw_units_v3`) uses and which
removed the gorilla repeat clique (0/33). No E1 human number has been computed.

**Rule E1 (no new constant; the gorilla truth's setting):** `mcl_families --min-exonic-bp 1` (with the default
`--exonic-both-sides`): an edge needs one alignment record mapping exon bases of one gene onto exon bases of the other, and
the pair's merged alignments must cover >= 1 exonic base of the longer gene. Everything else unchanged (identity >= 0.70,
cov_longer >= 0.30, >= 300 bp, I = 2.8, prune 1e-9). E0 = the AI construction (`--min-exonic-bp 0`).

**Substrates.** Development (looked at): chr15/chr17/chr22 — RefSeq truth = `mcl_families` on `human2/genes.asm20.paf` +
`human2/refseq_c15_17_22.gff`; GENCODE = `lit/annot_gencode/all.paf` + `nodes.gff`. **Hold-out (never looked at): chr16,
chr19, chr20** — RefSeq = gene+pseudogene records of `winloci_data/Reference/chm13v2.0_RefSeq_full.gff.gz` (exon union
from exon lines' `gene=`, as `bench/guided_min.py load_genes`); GENCODE = gene records of the CAT GFF (exon union from
transcripts, as AI); BOTH annotations through the identical `bench/node_graph_mcl.py` prep (3 Mb chunks) / align / mcl.

**AJ-1 (decision, held-out substrate; development reported):** under E1 the definition is ANNOTATION-ROBUST iff GENCODE_E1
scored against RefSeq_E1 (`bench/rna_truth.py score`, all loci in truth clusters >= 2) has pairwise sensitivity >= 0.90,
pairwise precision >= 0.90 and bipartite F >= 0.90. E0 vs E0 reported alongside.

**AJ-2 (guard: E1 repairs the truth, it does not degenerate it; each substrate):** HGNC (`hgnc_complete_set.txt`) by
exact symbol. Over the RefSeq loci, H = pairs whose two gene symbols share >= 1 `gene_group_id`; X = pairs where both
symbols have a non-empty `gene_group_id` and share none. A locus's symbols = Names of the RefSeq records it contains
(development: records whose coordinates equal a clusters.tsv row; hold-out: the node's record). REPAIR iff
|pairs(RefSeq_E1) ∩ H| >= 0.90 × |pairs(RefSeq_E0) ∩ H| AND |pairs(RefSeq_E1) ∩ X| <= 0.50 × |pairs(RefSeq_E0) ∩ X|.
If AJ-2 fails on the hold-out, AJ-1 is not read as support whatever its value.

**AJ-3 (report only, development):** the human arms of AH (de novo D = read nodes + MCL; guided minimal G) rebuilt under E1
from their existing PAFs and scored against RefSeq_E1 at DNA and RNA level, with the goal bars. Because the truth changed,
no AJ-3 number is a confirmation; a gain needs its own held-out test.

**Expected in advance:** E1 removes most cross-group pairs in the largest clusters; the hold-out chr19 KRAB-ZNF
superfamily is the likeliest place for AJ-1 to fail (dense, graded homology where MCL cuts depend on the node set).

---
## ADDENDUM AK (2026-09-14; user chose the adjudicated truth; before any adjudicated number exists) — an adjudicated two-annotation ground truth, evidence calibrated on agreement pairs, an independent third annotation held out

**Why:** AI/AJ/§6km — two expert annotations through the identical construction agree only at bipartite F ~0.8, and the
residual is how each annotation models the same genes. User decision (09-14): measure the goal against an ADJUDICATED truth.
Pairs both annotations group are true. Pairs only one groups are settled by independent evidence or left unscored.

**Construction rule for every annotation: E1** (`mcl_families --min-exonic-bp 1`, defaults otherwise), via
`bench/node_graph_mcl.py`. The one exception is development RefSeq (`lit/aj_dev/refseq_e1`, same construction).

**Joint loci.** Records of RefSeq and GENCODE/CAT (`bench/annotation_nodes.py`) whose exon unions share >= 1 bp on a contig
are one locus (union-find; the `--merge-overlapping-loci` rule). Locus exons = union of all member records' exons.

**Opinions.** A record's family comes from `<tag>.clusters.tsv`; a record folded away by the construction takes its
representative's family via `<tag>.loci.tsv`. For locus u, F_A(u) = the families of u's records in annotation A.
- A = SAME for (u, v) iff F_A(u) ∩ F_A(v) ≠ ∅.
- A = DIFF iff both loci have A records and do not share a family.
- A = NONE otherwise.

**Pair status** (only pairs with at least one SAME are enumerated; every other pair is FALSE):
- AGREED TRUE: SAME in both annotations.
- DISPUTED: SAME in exactly one. Adjudicated as:
  - TRUE if any evidence (below) holds;
  - FALSE if both loci are coding (a CDS in either annotation) and no evidence holds;
  - UNSCORED otherwise.

**Evidence (no new constants; the definition's identity 0.70 and coverage 0.30; SEDEF's 1 kb resolution):**
- (P) Protein. The longest CDS of u, over both annotations, is translated and aligned with miniprot 0.18 to the genomic span
  of v. P holds iff an alignment has Identity >= 0.70 and covers >= 0.30 of the protein, in either direction.
- (S) Segmental duplication. There is a UCSC hs1 SEDEF pair (`soto_replication/sd_v1.bed`) with exon bases of u in one side
  and exon bases of v in the other. The linear projection of u's exonic interval within that side must land within
  |len_A − len_B| + 1 kb of v's exonic interval in the other side.

**AK-0 (evidence validity; each substrate; gates AK-1).**
- Evidence sensitivity: among AGREED TRUE pairs with both loci coding, (P or S) holds for >= 0.80.
- Evidence false-positive rate: among coding pairs whose genes share an HGNC `gene_group_id` but are DIFF in both
  annotations, (P or S) holds for <= 0.20.
- Coverage: <= 0.50 of DISPUTED pairs end UNSCORED.

If any of the three fails, the adjudicated truth is NOT VALID on that substrate and AK-1 is not read.

**Truth objects.** TRUE pairs; clusters = connected components of TRUE pairs. Levels:
- DNA: all loci.
- RNA: loci with u >= 3 (`interval_expression.py` counts on locus spans), clusters recomputed on the TRUE pairs among expressed
  loci. This is reported only.

**Scoring (`bench/adjudicated_truth.py score`).** A method's copies are assigned to loci by best overlap.
- Predicted co-family pairs count TP if TRUE and FP if FALSE; UNSCORED pairs are ignored.
- Pairwise sensitivity = TP / |TRUE|; precision = TP / (TP + FP).
- Bipartite R/P/F over the loci in truth clusters (`guided_pipeline.bipartite`).

**Substrates.** Development: chr15/chr17/chr22. **Fresh hold-out, never used: chr5, chr7, chr21** (6,678 RefSeq genes).

**Arms.**
- **ENSEMBL**: Ensembl rapid-release CHM13 gene set `Homo_sapiens-GCA_009914755.4-2022_07-genes.gff3.gz` (md5 5508409d; gene,
  ncRNA_gene and pseudogene records; exon union from their transcripts; contigs renamed N -> chrN) under E1. NOT used to build
  the truth. Caveat fixed in advance: Ensembl derives part of its models from GENCODE and paralogue mapping.
- REFSEQ_E1 and GENCODE_E1: inputs to the truth — reported, never confirmations.
- G: guided minimal, 50% RefSeq seeds (seed 1) + candidates, under E1. Its seeds come from a truth input — reported.
- D: de novo read gene nodes (`bench/read_gene_nodes.py`) + construction under E1, RNA level. Hold-out BAM =
  `A119b_ds.bam` chr5/7/21.

**AK-1 (decision, hold-out; development reported).** THE GOAL IS MET for the definition on an independent annotation iff AK-0
passes AND ENSEMBL at DNA level has pairwise sensitivity >= 0.90, pairwise precision >= 0.90 and bipartite F >= 0.90.
The goal bars are reported for G (DNA) and D (RNA), with "not independent of the truth" stated for G.

---
## ADDENDUM AL (2026-09-14; after AK-0 NOT VALID on development; before any method is scored against the revised truth) — edge-level adjudication with nucleotide evidence; the truth construction frozen

**AK development result (fixed, reported):** AK-0 NOT VALID on chr15/17/22.
- Protein evidence sensitivity on agreed coding pairs was 0.362 (bar 0.80); disputed pairs unscored 0.746 (bar 0.50);
  hard-negative rate 0.005.
- Two causes:
  - Amino-acid Identity is not on the definition's nucleotide scale: agreed direct edges align at 85-95% protein coverage but
    0.25-0.58 identity (GOLGA6L6 vs GOLGA6L2 / GOLGA8).
  - Pair-level adjudication tests direct homology for co-membership that is transitive.
- The chr5/7/21 hold-out was NOT touched.

**Revised truth construction (implemented in `bench/adjudicated_truth.py` build, commit to follow; frozen here):**
1. **Joint loci** as AK.
2. **Edges:** each annotation's E1 graph (`<tag>.graph.tsv`), mapped to joint loci.
   - AGREED edges are present in both graphs.
   - DISPUTED edges are present in one.
3. **Evidence X (u -> v):** dc-megablast (BLAST+ `blastn -task dc-megablast -lcase_masking -evalue 1e-5`, soft-masked CHM13
   v2.0) of u's exon-union sequence (both annotations) against v's genomic span.
   - Non-overlapping HSPs on the query, greedy by nident, must sum to >= 300 aligned bp at identity >= 0.70 and cover >= 0.30
     of u's exonic length.
   - Either direction counts.
   - These are the definition's own constants; 1e-5 is the BLAST homology convention.
4. **Evidence S:** SEDEF as AK.
5. **Graphs.** STRICT = agreed + disputed with (X or S). PERMISSIVE = agreed + all disputed.
   - Edge weight = mean of the annotations' weights.
   - MCL (`bench/mcl_port.py`, I = 2.8, prune 1e-9) on each.
6. **Pair status:**
   - TRUE if co-clustered in both;
   - UNSCORED if co-clustered in exactly one;
   - FALSE otherwise.
7. **Truth clusters** = connected components of TRUE pairs.

**Validity gate AL-0 (each substrate; development result: 0.861 / 0.006 / 0.198 = VALID):**
- evidence (X or S) holds for >= 0.80 of agreed edges;
- evidence holds for <= 0.20 of HGNC hard negatives (coding loci sharing a `gene_group_id`, annotated in both, no edge in
  either graph);
- UNSCORED <= 0.50 of TRUE + UNSCORED.

**Scoring:** as AK. Loci are assigned to method families by best overlap; TP on TRUE, FP on FALSE, UNSCORED ignored;
bipartite on truth clusters.

**Arms and substrates** as AK (ENSEMBL, REFSEQ_E1, GENCODE_E1, G, D). Development chr15/17/22 is reported; the decision is on
the **fresh hold-out chr5/chr7/chr21**.

**AL-1 (decision, hold-out).** THE GOAL IS MET for the definition on an independent annotation iff AL-0 is VALID on the
hold-out AND ENSEMBL (DNA level) has pairwise sensitivity >= 0.90, pairwise precision >= 0.90 and bipartite F >= 0.90.
- G (DNA) and D (RNA, expressed loci u >= 3 on joint-locus spans) are reported against the bars.
- REFSEQ_E1 / GENCODE_E1 are reported as truth inputs.
