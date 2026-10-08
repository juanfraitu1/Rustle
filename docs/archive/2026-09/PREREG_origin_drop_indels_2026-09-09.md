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

---
## The requested next step (2026-09-09): direct sequence confirmation + sweep to other families

### ⭐⭐⭐ Direct sequence confirmation — the indels are real, not alignment noise
Fetched the exact genomic locus each origin certificate tested (`locus_start`–`locus_end`, human CHM13),
aligned 2 example reads per dominant locus against it directly (`minimap2 -ax splice -uf -k14 --eqx`), read
the raw CIGAR by eye:
- **Copy 2** (chr16:14,938,074–14,953,180): two independent reads both show `...358=57I...` — a 57 bp
  insertion at the **identical alignment position**, everything else clean matches/mismatches.
- **Copy 22** (chr16:75,773,416–75,821,281): two independent reads show `...8=492I196N5=1X12=1X11=1X3=5D2=
  2D4=8D141=1X25=...` / `...14=1X8=476I176N...` — a ~480–500 bp insertion at the identical position,
  followed by an **identical downstream CIGAR** (same intron, same run of matches/mismatches) between the two
  independent reads.
- **Specificity check**: the copy-2 reads, aligned instead against copy 22's locus, show scattered mismatches
  and large soft-clips — they do not fit there. Confirmed both ways.
⟹ Two independent, recurring, position-identical indels of different characteristic sizes (57 bp, ~490 bp),
present across multiple independent molecules with clean matching sequence immediately before and after —
this is unambiguous evidence of real individual-specific structural variants relative to the CHM13 reference
used to build the catalog (very plausibly VNTR/tandem-repeat copy-number differences, common at these
NPIP/LCR16 loci), not evidence the reads come from elsewhere. **Not alignment noise, not a coincidence.**

### The sweep to other families — one uninformative, one small-and-consistent
| family | copies | AS-tied | contested (base) | contested (+odi) | reading |
|---|---|---|---|---|---|
| MCL1_073242 (produced the fix) | 80 | 716 | 725: 14/409/302 | 1,143: 230/531/382 | the large effect |
| **MCL2_073244 (held-back, trap 15)** | 64 | **1** | 0/0/0 | 0/0/0 | ⚠ **genuinely uninformative** — too sparse an AS-tied population to test either way; NOT a validation, NOT a refutation |
| MCL7_073242 (fresh, unexamined before today) | 32 | 99 | 9: 0/9/0 | 11: 0/10/1 | small, same-direction, non-alarming: +2 contested, `assigned` stays 0, no new false calls |

⚠⚠ **The held-back family did not validate or refute this**, and that must be said plainly rather than
glossed over — trap 15 exists to catch overfitting, and it cannot do that job on a substrate with no
population to overfit. MCL7 is a genuine (if small) second data point and behaved consistently: it gained
contested molecules, gained none falsely to `assigned`, and produced no alarming shift.

## Revised recommendation
The mechanism is now confirmed at the sequence level, not just the statistical level — the two dominant loci
in the human family really do carry the indicated indels. The gorilla sweep adds one small consistent
corroboration (MCL7) and one uninformative null (MCL2). **Still not proposing a default flip**: no OUT-of-family
population beyond MCL7's small one has been checked, and MCL7's own sequence-level confirmation (does its own
newly-resolved molecule really carry a clean recurring indel, the same direct check done above for human) has
not been done. That remains open before any default decision.

---
## 2026-09-09, later: MCL7 sequence check, the "are the tied copies actually similar?" audit, and a blind spot

### MCL7's two entering molecules — the rule confirmed on a DIFFERENT mechanism (`scratchpad/sanity/mcl7_*`)
Both are clean by substitution and were rejected by the tightness of the Binomial null on short reads plus
scattered 1-bp indels (typical HiFi homopolymer errors), not by a large SV:
- `SRR27438212.8747054` → copy 1: `667=1D517=1D987=1D271=1D98=1X90=1I14=1I63=…23S` — **X = 1 over 4.6 kb**,
  eleven scattered 1-bp I/D, a 23 bp soft-clip; pre-fix edits 35 vs mean 13.8 (z 5.7, REJECT).
- `SRR27438212.5998780` → copy 6: `60=1X21=1X53=2X8=1X2=1X289=1I184=` — **X = 7 over 623 bp**, one 1-bp I;
  pre-fix edits 7 vs mean **1.9** (z 3.8, REJECT). Against copy 1 it is a mess of soft-clips and mismatches —
  specificity holds.
⟹ Human MCL0 exercised the large-recurring-SV case (57 bp, ~490 bp); MCL7 exercises the short-read /
scattered-1-bp-indel case. Same rule, both consistent with "indels are not origin evidence".

### ⭐⭐ The audit the user asked for: are AS-tied molecules tied between copies that are actually similar?
Two granularities, because whole-locus identity is the wrong one (a read covers ~1–3 kb, usually exonic;
copies can be 95 % identical over a 48 kb locus with diverged introns yet near-identical over a read):
- whole-locus, 26 human copies, minimap2 asm20 all-vs-all: 325 pairs, identity 0.936–0.998, median 0.960.
- **read-footprint identity at the read's WORST tied placement** (`=`/`X` from the `--eqx` CIGAR), 216
  newly-assigned: **median 0.9926, min 0.9601**; substitutions at the worst tied placement median 20, max 69.
  **190 of 216 ≥ 0.97**: the tied copies are genuinely near-identical over the read and PSV columns resolve a
  fine distinction. **26 < 0.97** (copy 22 ↔ 23/24/25, whole-locus 0.9486): NOT near-identical over the read.

### ⭐⭐⭐ And what the tie actually IS for most of them — a finding that reframes the gate
Classifying the 216 by whether the ASSIGNED copy is even inside the AS-tie set:
| | n |
|---|---|
| assigned copy **BELOW the AS tie** (aligner under-scored it), footprint ≥ 0.97 | **180** |
| assigned copy BELOW the AS tie, footprint < 0.97 | 26 |
| assigned copy IN the AS-tie set | 10 |
Anatomy of a typical one (`…4369/ccs/47_2588`): assigned copy 2 placement **AS 2354, X = 5, I = 58,
primary = true**; the AS tie is at copy 6, **AS 2424, X = 15, I = 1**. The real 58 bp insertion costs ~50–70 AS
in gap penalties, so the aligner's BEST-scoring placements are the WRONG copies (6/7/8 — which ARE near-
identical to each other, so the tie among them is real), and the right copy — which minimap2's chaining stage
itself picked as PRIMARY — sits below the tie. Read-star ignores AS, realigns against every family copy, and
finds copy 2 on substitutions (X 5 vs 15–22). **The assignment is justified; the tie is the aligner's, among the
wrong copies.** For the 26 low-footprint cases (copy 22, ~490 bp insertion) the same story, more extreme.

### ⚠⚠ A blind spot in the AS-tied gate, quantified (`bakeoff/human/blindspot.txt`)
If the true copy is under-scored by an SV gap penalty and the AS-best lands **uniquely** at ONE wrong copy,
there is no tie, the gate calls it a "clear best" unique mapper, and O2 never sees it. Counted over all
molecules with a placement in a copy, human MCL0:
| class | n |
|---|---|
| AS-unique, agrees with primary (true unique mapper) | 13,199 |
| AS-tied, primary IN the tie | 6,958 |
| **AS-unique at a copy ≠ primary copy — gate skips, aligner disagrees with itself** | **2,278** |
| AS-tied, primary BELOW the tie (the 206 above live here) | 666 |
| AS-unique outside every copy, primary inside one | 111 |
Of the 2,278: **1,456 (64 %) carry ≥ 50 bp of insertion in their PRIMARY placement** (median 63, q90 401) — the
same SV-under-scoring mechanism. `--as-tie-ratio` does NOT cover it: AS(primary)/AS(best) median 0.978, so
0.98 admits only 42 %, 0.95 admits 79 %, 0.90 admits 94 %. ⚠ These reads are not mis-assigned by O2 — they are
invisible to it — but they are exactly the ambiguous reads O2 exists for. **Design question for the user, not
implemented**: admit to the gate a molecule whose PRIMARY copy ≠ its AS-best copy (the aligner's two stages
disagree), in addition to AS ties. Human MCL0 would gain up to 2,278 candidates; their read-star outcome is
unmeasured.

### Standing recommendation
`--origin-drop-indels` is confirmed at the sequence level on both mechanisms (large recurring SV; short-read
scattered indels) and on three families (large effect / uninformative / small consistent). Still the user's
default decision. The blind spot above is the more consequential open item.
