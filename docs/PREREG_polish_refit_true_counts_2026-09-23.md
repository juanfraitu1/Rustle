# Pre-registration — re-fit the assembly polish on TRUE read counts (`--keep-coordinate-duplicates`)

**Written 2026-09-23 (§6z8), before any grid point is scored.** Follows the miss taxonomy of
`docs/PREREG_external_tool_bakeoff_2026-09-22.md` §"What is missing" and r1058.

## Why this and nothing else

The taxonomy of every multi-exon reference transcript we miss (human chr20/21/22, 7,446 misses; gorilla
genome-wide, 3,000 sampled of 69,987) says:

- **≥ 95% have at most ONE read carrying the exact chain** (human: 86-96% zero reads in every structural class;
  gorilla 2.8% of misses have ≥ 2). Pass-1 needs two. The tools match only 4-18% of those classes, and only
  isoseq collapse, at singleton level. Unreachable by any rule short of emitting single-read chains.
- The **≥ 2-read residue is 384 human / ~1,900 gorilla (5.1% / 2.8% of misses)**, and it is not structural:
  303 / 74 of them are in our RAW output and were removed by the polish (fraction 170, ISM 110, both 23 on
  human); the remaining 81 human are the coordinate de-duplication (r1058). "Super" (we emit the longer chain,
  the reference is a contained shorter isoform) is the only structural class with ≥ 2-read mass, and it IS the
  ISM collapse.
- Two candidate rules on that residue were priced and fail: an absolute-support exemption from the fraction
  rule (5 reads: +33 chains for +330 transcripts, precision 14.5 → 14.3) and an ISM exemption for 3′-shorter
  sub-chains that pass an internal-priming test (reference-true ones show the priming signature 29% vs 54%
  for the rest: ~+40 chains for ~+1,200 transcripts). r861 already refuted the 3′-side ISM exemption without
  the priming test.

What remains is the **interaction found in r1058**: with duplicates kept, the RAW output gains 81 chains but the
polish — every dial fitted on de-duplicated counts (§6p8–§6q4) — removes 228 more than that. The dials may
simply be wrong for true counts. This is a re-fit, not a new rule.

## Arms

Binary = current (`--keep-coordinate-duplicates` on, both AS-table envs unset). Grid on the polish dials:
`--polish-isoform-fraction ∈ {0.02, 0.03, 0.05}` × `--polish-ism-ratio ∈ {0.7, 0.5, 0.35}` ×
`--polish-mono-quantile ∈ {0.82, 0.90}` — 18 points, `--polish-mono-shadow` and `full` as shipped.
Baselines already on disk: shipped dials on de-duplicated counts (2,295 chains, 23.3/16.0) and shipped dials on
true counts (2,148, 21.8/13.7).

## Substrates and selection rule

- **Development: human A119b chr20/21/22** (RefSeq subset, 10,851 transcripts). **Selection rule, fixed now:**
  the grid point with the most matching intron chains among those whose intron-chain precision is ≥ 16.0 (the
  de-duplicated shipped value). If none reaches 16.0, the point maximising `chains × precision`.
- **Held out: gorilla `GGO_mm.bam`, all 26 contigs**, at the ONE selected point, no re-tuning.

## Bar — judged on gorilla, against the SHIPPED setting on de-duplicated counts (25,829 chains, 27.0/33.1)

| outcome | verdict |
|---|---|
| chains ≥ 25,829 **and** intron-chain precision ≥ 33.1 | ⭐ **DOMINATES** — true counts + re-fit beat the shipped polish; flip both defaults |
| one of the two, the other within 1.0 pt / 1% | ⚠ **TRADE** — report, user's call |
| neither | ⛔ **NO** — the de-duplication key stays the default duplicate filter |

**Predicted, before looking:** ⚠. The raw gain is +3% chains (human) and the polish under true counts sits at
36.3 precision on gorilla, so a looser fraction or ISM ratio should buy chains back; but r863 found the
testis-panel dials have no single feasible point across chromosomes, and this library is 6× deeper, so I
expect the selected point to recover the chains at a precision between 33.1 and 36.3 — a trade, not a
dominance.

I will not change the grid, the selection rule, the substrates or the bar after seeing any number.

---

# OUTCOME (2026-09-23) — ⛔ **NO: every true-count grid point is dominated by the de-duplicated shipped setting.**

Human chr20/21/22, `--keep-coordinate-duplicates`, 18 points (mono quantile 0.82 vs 0.90 never changed a
chain count, so the table collapses to 9):

| fraction | ISM ratio | transcripts | intron chain SN / PR | matching chains |
|---|---|---|---|---|
| 0.02 | 0.70 (shipped) | 17,296 | 21.8 / 13.7 | 2,148 |
| 0.02 | 0.50 | 18,351 | 22.0 / 13.0 | 2,165 |
| 0.02 | 0.35 | 19,032 | 22.2 / 12.5 | 2,182 |
| 0.03 | 0.70 | 15,714 | 20.5 / 14.3 | 2,017 |
| 0.03 | 0.50 | 16,660 | 20.7 / 13.5 | 2,032 |
| 0.03 | 0.35 | 17,280 | 20.8 / 13.1 | 2,048 |
| 0.05 | 0.70 | 13,578 | 18.6 / 15.3 | 1,828 |
| 0.05 | 0.50 | 14,358 | 18.7 / 14.4 | 1,837 |
| 0.05 | 0.35 | 14,902 | 18.8 / 14.0 | 1,853 |
| **de-duplicated counts, shipped dials** | | **15,876** | **23.3 / 16.0** | **2,295** |

No point reaches precision 16.0, so the selection rule falls to `chains × precision`, which picks the shipped
dials (2,148 × 13.7). That arm is the one already run on gorilla in r1058: **24,733 chains at 36.3 vs the shipped
25,829 at 33.1 — chains 4.2% below, outside the 1% trade band ⇒ ⛔ NO.** The whole true-count frontier lies
below the de-duplicated point: loosening the ISM ratio buys +34 chains for −1.2 pts, tightening the fraction
buys +1.6 pts for −320 chains, and nothing recovers the 147 chains the true counts cost.

## What this settles

Coordinate-identical reads do not add independent evidence. A 2-read chain after de-duplication (two molecules
with different ends) is a better support criterion than two reads with identical ends, at every dial — which
is what one expects if identical-coordinate reads are largely PCR duplicates of one molecule (Iso-Seq
libraries are amplified; the identical rate here is 25-39% of primary records). **The `(chrom, start, end,
chain)` key is therefore a legitimate duplicate filter, not only a double-fetch guard, and stays the default.**
`--keep-coordinate-duplicates` stays as the opt-in. Consequence for the tool comparison (r1065): the 389 / 237
consensus chains all three tools emit and we do not are supported only by coordinate-identical reads — the
tools count PCR copies as support; we do not. Predicted ⚠ (a trade); measured ⛔ (dominated) — the prediction
overestimated what the polish could buy back.
