# Pre-registration — does the containment escape improve FAMILIES, not just admission?

**Written 2026-09-22, §6x5, before any family is scored.** §6x4 measured a ceiling; this is the adoption
test it explicitly deferred: *"Admission is a gate, not precision."*

## What was built

`mcl_families --min-cov-shorter <C>` (`annotation_families.rs`, deferred path — the one
`--min-exonic-bp 1` takes). A pair whose `cov_longer` fails ALSO passes if the alignment covers ≥ C of the
**shorter** gene's exonic length; the edge weight then uses `cov_shorter`. Default `0.0` = OFF.
**Verified byte-identical when unset** (`clusters.tsv` and `loci.tsv` both `cmp`-clean against the shipped
chr16 de novo catalog with the new binary). New params keys appended LAST (r936).

⚠ The shipped exon conjunct (`--min-exonic-bp 1 --min-shared-exon-frac 0.60`) applies to escaped pairs
unchanged — a **stronger** guard than the `exonic > 0 both sides` §6x4 measured the ceiling with.

## Arms · substrates · truth

**C = 0.90 primary** (§6x4: 216/679 victims at 2.1× hub), **C = 0.70 reported** (234 victims, 2.2×).
Baseline is the shipped catalog.

- Development: chr16 de novo + guided.
- **Held out: chr2 / chr8 / chr10** — the verdict is taken here, with no re-tuning.
- Truth: **protein-family referee** (neutral, §6u5) primary; Soto reported as context and labelled
  SD-circular (r741/r1085).

## The bar — committed now

| outcome | verdict |
|---|---|
| held-out F **up on ≥2 of 3** chromosomes, none down by >0.02, largest component ≤ **2.5×** baseline | ⭐ **ADOPT** |
| F within ±0.02 everywhere but sensitivity up and precision down by ≤0.02 | ⚠ **PARTIAL** — a recall knob, document and leave OFF |
| F down on the held-out set, or the hub guard fails | ⛔ **NO** |

**Predicted, before looking — ⚠ PARTIAL: sensitivity up, precision down, F roughly flat.** The escape
fires only where `cov_longer` already failed, and §6x3/r1002 showed those victims align along ~100% of
their own length, so they should be real paralogs and recall should rise. The risk is entirely on
precision: chr16 gains 365 nodes and 805 edges, and MCL routes around nothing — more edges means more
merging. ⚠**r911 is the standing warning**: gains measured at one layer did not survive being combined at
another ("MCL is what makes the conjunct affordable; gains are NOT additive").

I will not change the arms, the substrates, the truth or the bar after seeing any number.

---

# OUTCOME (2026-09-22) — ⭐ **CLEARS THE BAR. Small, consistent, and precision never falls — with one documented exception.**

`--min-cov-shorter` verified **byte-identical when unset** (`clusters.tsv` + `loci.tsv` `cmp`-clean against
the shipped chr16 de novo catalog, new binary), no new compiler warnings, new params keys appended last.

## Held out — chr2 / chr8 / chr10, protein-family referee (PRIMARY)

| chrom | arm | sens | prec | **F** | largest cluster |
|---|---|---|---|---|---|
| chr2 | A0 | 0.135 | 0.937 | 0.236 | 23 |
| | C=0.90 | 0.137 | **0.938** | **0.239** | 22 (0.96×) |
| | C=0.70 | 0.139 | **0.938** | **0.243** | 26 (1.13×) |
| chr8 | A0 | 0.279 | 0.989 | 0.436 | 46 |
| | C=0.90 | 0.283 | 0.989 | **0.440** | 46 (1.00×) |
| chr10 | A0 | 0.111 | 0.875 | 0.197 | 28 |
| | C=0.90 | 0.111 | 0.875 | 0.197 | 28 (1.00×) |

⭐**F up on 2 of 3, flat on the third, none down; precision up or equal in every cell; hub guard passes at
0.96–1.13× against a 2.5× bar.** That is the ⭐ ADOPT row as written.

## Development — improves in BOTH modes, which is what "all modes" asked for

| arm | sens | prec | **F** |
|---|---|---|---|
| chr16 de novo A0 | 0.121 | 0.949 | 0.214 |
| chr16 de novo C=0.90 | 0.131 | **0.952** | **0.230** |
| chr16 guided A0 | 0.186 | 0.966 | 0.312 |
| chr16 guided C=0.90 | 0.190 | **0.967** | **0.317** |
| chr16 guided C=0.70 | 0.193 | **0.967** | **0.322** |

## ⚠ The exception, stated plainly

**NPIP, chr16, guided mode gets WORSE**: Soto F 0.833 → 0.800 (sensitivity 0.750 → 0.700, collapsed 0 → 1);
on the §6x1 union truth U2 it is flat at 0.721 but collapse rises 1 → 3. **NPIP de novo improves** (Soto
0.727 → 0.750 at **precision 0.923 → 1.000**, collapsed 2 → 1; U2 0.610 → 0.621). So the escape helps the
de novo mode on the thesis's own family and hurts the guided mode there, while helping guided on the
chromosome-wide referee. ⚠**Do not quote the chromosome-wide gain as an NPIP result.**

## Effect size — reported honestly

The held-out gains are **+0.003 to +0.007 F**, i.e. one or two gene-pairs. They are consistent in sign
across arms, truths and both modes, and precision never falls, which is why they clear the bar — but this
is a small effect and must not be presented as a headline improvement. Soto (context, SD-circular) agrees
on direction for chr2 (0.778 → 0.804) and chr8 (flat) and **disagrees on chr10** (0.750 → 0.739).

## Recommendation

**Clears ⭐ ADOPT on the pre-registered bar, but the default stays OFF pending the user's call** — the same
handling as `RUSTLE_JUNCTION_MAJORITY` (§6m8), and for the same reason: a real, guard-clean, held-out-
positive effect that is small and has one known regression (NPIP guided) deserves an explicit decision,
not a silent default flip.

## Prediction scorecard — wrong, and in the informative direction

I predicted ⚠ PARTIAL: *"sensitivity up, precision down, F roughly flat … the risk is entirely on
precision: chr16 gains 365 nodes and 805 edges, and more edges means more merging."* **Precision rose in
every cell and the hub guard never came close to firing.** The reason is the one §6x4's outcome had
already flagged and I failed to carry into the prediction: the ceiling was measured with
`exonic > 0 on both sides`, but the shipped pipeline applies `--min-shared-exon-frac 0.60`, which is far
stronger. The real implementation admitted **+37 nodes on chr16, not +365** — an order of magnitude less
than the ceiling. ⭐**A Python ceiling measured with a weaker guard than the shipped one overstates both
the gain and the risk by ~10×; it is a feasibility screen, not a forecast.**

## ADDENDUM — the third mode, run to complete the goal's "in all modes"

Semi-guided (SD-region nodes, §6x0) uses the same `mcl_families` gate, so the flag applies there too.
It was not in the pre-registered arms; run afterwards and reported as an addendum, not as a bar test.

| arm | referee F | referee prec | NPIP Soto F |
|---|---|---|---|
| semi `exonic` A0 | 0.209 | 0.771 | 0.750 |
| semi `exonic` C=0.90 | 0.207 | **0.725** | 0.750 |
| semi `whole` A0 | 0.210 | 0.973 | 0.750 |
| semi `whole` C=0.90 | **0.175** | **0.833** | **0.645** |

⛔**In semi-guided mode the escape is HARMFUL** — the only mode where precision falls. ⭐**The mechanism
is clear and it bounds where the flag is valid**: an SD region has no intrinsic gene boundary, so a small
region sitting fully inside a large one is the *normal* state of that node set, not evidence of paralogy.
Containment is only informative when a node is a **gene-like unit** — an assembled locus or an annotated
gene body. ⟹ the flag must never be enabled with `--from-genome-sd` nodes.

**Final picture across the three modes:** de novo ⭐ (F .214→.230, precision up) · guided ⭐ chromosome-wide
(.312→.322, precision up) but ⛔ on NPIP specifically (.833→.800) · semi-guided ⛔ (precision .973→.833).

## ADDENDUM 2 — two further arms run against the session goal, both closing negatively

**1. Positional-overlap guard (r1010) — REFUTED and reverted.** Hypothesis: the semi-guided regression is
caused by SD regions overlapping each other (measured 87.2% vs de novo 34.5% / guided 22.4%), so the escape
should be refused for overlapping pairs. Built and scored on every mode: **no effect on semi-guided**
(0.175 either way) and it **HURT de novo, F .230 → .199, below the .214 baseline**. In a duplicated region
two overlapping de novo loci are frequently genuine tandem copies. Reverted; unset byte-identity and
pre-guard output both re-verified.

**2. The false-positive side is NOT solved (r1011).** Rust-to-Rust on the dumped pre-MCL graph, the escape
gains **37 nodes and loses 0**, of which **12 are §6w2 collateral evictions = 3.0% of the 406**. Splitting
recovered 0% and boundary pull-in is negative, so the cost over-merge imposes on innocent partners remains
**~97% unaddressed**. ⚠I first reported 13.3% from a Python gate that omits the shipped exon conjunct;
the same reimplementation implied 137 nodes were "lost", which a monotone operator cannot do — **checking
the monotone invariant ("lost must be 0") caught it immediately.**

## Honest status against the session goal

| | state |
|---|---|
| node FALSE NEGATIVES | ⚠ **partly improved.** Cause identified (r1002), both shrink remedies closed (r1001, §6w6), `--min-cov-shorter` ships behind a flag: held-out F up 2 of 3, precision up-or-equal everywhere, +0.003–0.007 |
| node FALSE POSITIVES (over-merge) | ⛔ **not solved.** Diagnosed and bounded; every remedy tried is closed; the escape treats 3.0% of the collateral cost |
| "in all modes" | ⛔ **not achieved.** de novo ⭐ · guided ⭐ chromosome-wide but ⛔ on NPIP · semi-guided ⛔ harmful, and r1010 shows the regression is not patchable |
| default | OFF, uncommitted — the user's call, per the standing rule that commits happen only when asked |

## ADDENDUM 3 (§6x6) — C is not a fitted threshold, and 0.70 beats the pre-registered 0.90

Swept on **development only** (chr16; the held-out set was not used to choose C):

| C | de novo F | guided F | largest de novo cluster |
|---|---|---|---|
| OFF | 0.214 | 0.312 | 26 |
| 0.40 – 0.80 | **0.230** | **0.322** | 33 |
| 0.90 – 0.95 | 0.230 | 0.317 | 33 |

⭐⭐**Identical to three decimals across a 2× range of C.** A parameter whose output is constant over its
usable range is not a tuned constant — the specific property the advisor asks for. **C = 0.70 (midpoint of
the flat region) dominates the 0.90 chosen from §6x4's ceiling**: held-out chr2 0.243 vs 0.239, chr8 0.440,
chr10 0.197, guided dev 0.322 vs 0.317. ⚠It also means the flag cannot be tuned for more — this is its
ceiling on this substrate.
