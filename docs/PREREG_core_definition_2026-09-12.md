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
