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
