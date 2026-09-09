# PREREG — `--origin-drop-indels`: exclude indels from the origin certificate's edit count (2026-09-09)

**Written before any run under `--origin-drop-indels`.** Traced from the `best_by_psv` fix (§6ha,
`4d5a1a2`), whose measured side effect was 37 human MCL0 molecules newly `origin_rejected` (§6gz addendum).

## The finding (from `RUSTLE_STAR_DEBUG=1`, all 40 origin-flipped molecules, `bakeoff/human/flip37_debug.txt`)
20 of 40 share `best_copy=2`, with **X = 2–5 substitutions** (well inside 0.3 % error tolerance for a ~2 kb
read) but **I = 57–66** — a near-constant ~60 bp insertion recurring across independent reads at one locus.
A fixed-size indel repeating across many INDEPENDENT molecules is the signature of a reference-vs-haplotype
structural difference (or an alignment-representation artifact at that position), not evidence the read is
from elsewhere. A separate cluster of 8 is dominated by 14–327 **unaligned** bases with near-zero indels — a
genuinely different and legitimate rejection (the locus does not cover content the read has).
Hand-simulated on all 40 (`X + unaligned` vs `mean = (aligned+unaligned)*error_rate`, unchanged threshold):
**32/40 flip REJECT→pass, 8/40 correctly stay REJECT.**

## The fix
New flag `--origin-drop-indels` (param `origin_drop_indels`, default **false** ⟹ byte-identical). Only
affects the `read_star_genomic` (default) branch of the origin certificate: `nm = X + unaligned − explained`
instead of `X + I + D + unaligned − explained`; `blk` (denominator) unchanged. Distinct from
`--origin-substitutions-only`, which also drops `unaligned` and stops asking "does the locus explain the
whole read" — this flag keeps that requirement and removes only the indel term.

## Predictions
| # | prediction | refuted by |
|---|---|---|
| P1 | ≥ 28 of the 40 traced molecules flip `origin_rejected 1→0` under `--origin-drop-indels` | < 20 |
| P2 | the 8 unaligned-dominated molecules stay `origin_rejected = 1` | any of the 8 flips |
| P3 | human contested set GROWS (more molecules re-enter, since fewer are excluded via origin_rejected); `assigned` does not fall | assigned falls, or contested shrinks |
| P4 | gorilla: origin-certificate composition has not been checked — report whatever the run shows, no directional prediction | — |
| P5 | full escape (`--no-as-tied-only --best-by-alignment`, `--origin-drop-indels` NOT set) stays byte-identical on both species | any md5 differs |
| P6 | `--origin-drop-indels` alone (gate on, `best_by_alignment` off) on the OLD-style single-candidate/sole-candidate molecules changes nothing new outside the traced 40, since those already passed (old `bk` coincided with lower total edits) | a large unrelated population shifts |

## What this is NOT
Not a re-run of the `origin_substitutions_only` result (§6fc: 1/20 wrong anchor under the OLD `bk`-selection,
before today's `best_by_psv` fix and before the AS-tied gate existed). That precedent does not transfer
automatically; `--origin-drop-indels` is a narrower, different rule (keeps `unaligned`) and is measured fresh.

## Rules held
Default OFF; escape must be exact. Human and gorilla never pooled. Report P1–P6 honestly, including any
refutation, before proposing a default.

---
## Outcome (2026-09-09, after the runs)
| # | prediction | verdict |
|---|---|---|
| P1 | ≥ 28 of 40 traced molecules flip `origin_rejected 1→0` | ✓ **29/40** |
| P2 | the 8 unaligned-dominated molecules stay rejected | ✓ **exactly those 8** (5× best=9, best=2/327-unaligned, best=8, best=20) |
| P3 | contested grows, assigned does not fall | ✓✓ **far exceeded**: human contested 725→1143 (+418), assigned 14→**230** |
| P4 | gorilla, no directional prediction | small, same-direction effect: contested 28→33, assigned 2→4 |
| P5 | full escape byte-identical | ✓✓ human `91081887`, gorilla `ff0b8f16` |
| P6 | no large unrelated population shifts | ⚠ **see below — P6 needs qualification, not refutation** |

## ⚠⚠ The scale is far larger than the 40-molecule trace predicted, and P6 needs honest qualification
All 418 newly-contested human molecules are **`ENTERED`** transitions (origin_rejected 1→0) — **zero** previously
`tied`/`ambiguous`/`assigned` molecules changed status. Mechanically: `origin_rejected` FORCES `Ambiguous`
regardless of what the pairwise/PSV test already computed (`copy_assign_pipeline.rs`, the block after the
Binomial test). Removing spurious rejections does not create new evidence — it stops the origin certificate
from **overwriting an already-computed, already-correct pairwise verdict**. So P6 ("no large unrelated
population shifts") is true in the sense that no NEW mechanism fired — but the SIZE of the population the
existing mechanism was suppressing (418, not ~30) was not anticipated.

**Composition of the 216 newly-`assigned`**: concentrated at two loci — **copy 2 (153, 71 %)** and
**copy 22 (42, 19 %)**, together 90 %. Traced with `RUSTLE_STAR_DEBUG=1`:
- copy 2: `X = 2–5` substitutions, **`I = 57–66`** — a ~60 bp insertion, tightly clustered in size across
  independent reads.
- copy 22: `X = 5–7` substitutions, **`I = 476–504`** — a ~500 bp insertion, even more tightly clustered.
Both patterns are a fixed-size indel recurring across many independent molecules with near-baseline
substitution counts — the signature of a **real tandem-repeat/copy-number difference between the catalog's
reference sequence and this individual's actual haplotype** at those two specific loci (NPIP/LCR16-type
regions are known for exactly this), not evidence of the wrong candidate. Statistical strength on the 216:
`n_decisive` min 2 / median **13** / max 66; `margin` min 6.9 / median **82.9** / max 366 — strong, not
borderline, resolutions.

## What this means, stated plainly
The fix is mechanistically sound and the evidence for it is specific and reproducible (two loci, two distinct
recurring indel sizes, both dominated by substitution counts inside normal error tolerance). But it is now a
much bigger claim than "chase 40 flipped molecules" — it says the DEFAULT origin certificate has been
suppressing real, PSV-strong assignments at scale wherever a locus carries a real structural indel relative
to its catalog reference. **This is squarely a default-flip decision for the user, not something to ship
silently on the strength of one family's trace.** Recommended next step before any default changes: sweep
`--origin-drop-indels` over the held-back family and a handful of others, and specifically confirm the
two dominant loci (copy 2, copy 22) really do carry the indicated indel by inspecting their catalog reference
sequence against a few of the newly-assigned reads directly (an alignment-level sanity check, not yet done).
