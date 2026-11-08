# PREREG row 1581: add weak SD-only edges on top of SD damping (`--sd-add-edges W`), held out on HG002

Written 2026-11-08, before any binary run of the flag. The setting is chosen on A119b (dev) by the rule below and appended
here before the held-out jobs start. Companion of row 1580 (`docs/PREREG_sd_damp_span_2026-11-08.md`); 1581 is judged
against 1580's arm, not against plain MCL.

## Question

On top of `--sd-damp 0.25` (row 1580), does an edge of small weight W between loci that one SD pair links but that have
no homology edge buy recall more cheaply than lowering MCL inflation does?

## Rule

Same SD link rule as row 1580 (`SdPairs::holds`: the copy holds >= half of the node's exon bases and >= half of its span).
Every SD-linked node pair without a homology edge gets an edge of weight W (`add_sd_edges`), applied after damping.

## Provenance (dev, A119b, numpy, post hoc)

`r1571/clu/recall_sd2.py`: of the 633 truth pairs the damped arm misses, 351 have no graph path (16 SD-linked) and 282
are MCL splits (110 SD-linked). W = 0.15 at I = 2.8 / 3.2 / 3.6: .890/.735, .904/.730, .910/.726, about +2-2.5
precision over the damped-only curve at equal recall. Merging clusters by SD links instead was below that curve (closed).

## Setting (fixed on A119b binary runs r1581dev before any held-out run)

Candidates: W in {0.1, 0.15} x I in {2.8, 3.2, 3.6}, all with `--sd-damp 0.25`. Baseline: the binary damped-only curve
(`--sd-damp 0.25` at I = 1.6, 2.0, 2.4, 2.8, 3.2, 3.6), `pairscore.py` on PREFIX.fam.copies.tsv. Choose the candidate
with the largest precision gain over that curve at its recall (linear interpolation between the two bracketing points)
among candidates whose recall is >= the damped default's (I = 2.8). If no candidate's recall lies inside the curve's
recall range, or none gains, stop: no held-out run.

## Held-out

HG002 (decisive) and gorilla OR6737 (descriptive, 88 truth pairs), exactly as row 1580: same SEDEF files, same runs,
`r1580/heldout.sh`, extended with the chosen W at every inflation of the curve.

## Bar (HG002)

PASS when both hold:

1. Precision of the chosen arm >= the damped-only curve's precision at its recall (bracketing interpolation) + 0.03.
2. Recall of the chosen arm >= recall of the damped default (`--sd-damp 0.25`, I = 2.8).

If its recall lies outside the damped curve's recall range: FAIL. Row 1581 is judged even if row 1580 fails; if 1580
fails, a 1581 PASS still only says the extra edges beat lower inflation on the damped graph.

## What follows

- PASS: `--sd-add-edges` documented as opt-in next to `--sd-damp` (defaults are the user's call).
- FAIL: register; the flag is removed or stays experimental at the user's choice.

## Setting chosen (2026-11-08, A119b binary r1581dev, before any held-out run)

Refactor check: `--sd-damp 0.25` I = 2.8 copies and clusters byte-identical to r1578span (shared `sd_linked_pairs`).
Damped-only curve (pairscore, copies): I = 3.6 .943/.688, 2.0 .909/.716, 3.2 .943/.716, 2.8 .925/.718 (default),
2.4 .914/.722, 1.6 .883/.736. The curve is not monotone in inflation (I = 3.2 beats I = 2.0 at equal recall); the
interpolation uses the points sorted by recall, as registered.

| candidate | precision / recall | damped curve at that recall | gain |
|---|---|---|---|
| **W 0.15, I 3.2** | **.920 / .733** | .890 (.914 at .722 to .883 at .736) | **+0.030** |
| W 0.15, I 2.8 | .906 / .733 | .890 | +0.016 |
| W 0.15, I 3.6 | .922 / .725 | .907 | +0.015 |
| W 0.10, I 2.8 | .916 / .726 | .905 | +0.011 |
| W 0.10, I 3.2 | .918 / .722 | .914 | +0.004 |
| W 0.10, I 3.6 | .922 / .711 | recall < .718 | excluded |

Chosen: **`--sd-damp 0.25 --sd-add-edges 0.15 --inflation 3.2`**. Its dev gain equals the bar (+0.030), so the held-out
has no margin. Held-out arms: `add_I*` in `r1580/heldout.sh SAMPLE 0.25 0.15` (verdict arm `add_I3.2`).
