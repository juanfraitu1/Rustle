# PREREG — what junction-jitter tolerance should RUSTLE_JUNCTION_FUZZ_BP use? (2026-09-15, written before any gffcompare/SQANTI3 score from this feature is computed)

## Question
On real chr20 IsoSeq reads, how far do individual reads' own CIGAR-derived splice junction coordinates scatter around the nearest real annotated RefSeq intron boundary -- and what tolerance, chosen by ONE fixed rule stated here, should `merge_fuzzy_skeletons` (docs/superpowers/specs/2026-09-15-fuzzy-junction-merge-design.md) use?

## Method
`bench/measure_junction_jitter.py` against `/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20/chr20.bam` and the same substrate's `chr20_ref.gtf`. Every spliced primary read's own junction is matched to the nearest annotated intron by donor-site distance (capture window 500bp); the signed offset at both donor and acceptor is pooled into one distribution.

## Result (initial run)
reads scanned: 25341
junctions seen: 101697
junctions matched to an annotated intron within 500bp: 100479
pooled donor+acceptor |offset| samples: 200958
median |offset|: 0
90th percentile |offset|: 0
95th percentile |offset|: 672

## Methodology correction (2026-09-15, before Task 5 runs)

**Issue identified:** The initial pooled statistic (p90 = 0bp) was methodologically flawed for the actual merge test's AND-condition semantics.

**Why:** `merge_fuzzy_skeletons` requires BOTH the donor AND acceptor site to be within `tolerance_bp`. However, pooling donor and acceptor offsets into one flat distribution was incorrect because:
- Donor-only p90 = 0bp (artificially tight, since `nearest()` selects introns BY donor-site distance — this measures search accuracy, not independent jitter)
- Acceptor-only p90 = 672bp (real, independent scatter)
- Pooling both drags the combined p90 down to 0bp, hiding the acceptor's real tolerance requirement

**Correct statistic:** Per-junction max(donor_offset, acceptor_offset), then the 90th percentile across junctions. This directly matches the AND-condition logic.

## Result (corrected)
Script re-run with corrected per-junction max aggregation:

reads scanned: 25341
junctions seen: 101697
junctions matched to an annotated intron within 500bp: 100479

Per-axis breakdown (for context):
- donor-only |offset| samples: 100479, median: 0, p90: 0, p95: 0
- acceptor-only |offset| samples: 100479, median: 0, p90: 672, p95: 4374

Per-junction max(|donor_offset|, |acceptor_offset|) [correct for AND-semantics]:
- samples: 100479
- median: 0
- 90th percentile: 672
- 95th percentile: 4374

## Decision (fixed now, before Task 5 runs)
`RUSTLE_JUNCTION_FUZZ_BP` = 672 (the per-junction max 90th percentile from the corrected statistic).

**Provenance, stated accurately (corrected 2026-09-16, final whole-branch review).** The design spec's
ORIGINAL pre-registered rule (`docs/superpowers/specs/2026-09-15-fuzzy-junction-merge-design.md`, "Tolerance
selection" step 2) was the 90th percentile of the POOLED donor+acceptor offset distribution from step 1
("Pool these into a real empirical distribution") — which gave **0bp** (see "Result (initial run)" above).
The per-junction-max statistic was NOT pre-registered from the start; it is a disclosed, post-hoc
methodology correction (see "Methodology correction" above) — the same "post-hoc fixes, disclosed" pattern
this project already uses elsewhere (e.g. `docs/PREREG_core_definition_2026-09-12.md`'s Addendum U and its
amendments U1'/U1''). It was substituted AFTER the pooled result (0bp) was seen but, crucially, BEFORE any
gffcompare/SQANTI3 score from this feature was ever computed. That timing — not a claim that the statistic
itself was chosen unchanged from the outset — is what makes the correction legitimate: no score from this
feature existed yet to have motivated picking a number that would merge more aggressively. Not changed based
on Task 5's result.
