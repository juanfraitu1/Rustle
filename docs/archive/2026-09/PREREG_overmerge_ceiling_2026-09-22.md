# Pre-registration — what is the CEILING of repairing node-level over-merge?

**Written 2026-09-22, §6x7, before the arm is built.** The session goal's false-positive half has been
bounded (r1015: 228 fused loci, 28/306 referee genes = 9.2%) and characterised (r1016: 49.6% trimmable
passenger, 34.2% two comparable genes) but every operator is refuted or carries an unsized cost. I
previously called sizing that trade "a design decision, not a measurement." **That was wrong — it is
exactly a measurement**, and this is it.

## The arm

**PERFECT ORACLE REPAIR.** Every de novo locus containing ≥2 fully-contained annotated genes is replaced
by one node per constituent gene (clipped to the locus); every other locus is passed through unchanged.
Then the shipped pipeline runs end to end: `minimap2 -c -N 50 -p 0.1 -x asm20 -X` all-vs-all →
`mcl_families --min-exonic-bp 1 --min-shared-exon-frac 0.60`.

⚠⚠**This is an ORACLE and can never ship** — it consumes the annotation to define de novo nodes, which
the de novo mode must not do. It exists only to answer: *if node over-merge were solved PERFECTLY, by an
oracle, how much family score would that buy?* That number bounds every possible trim, split or cut rule,
including the 49.6% trim whose cost I deferred.

## Metrics and the bar — committed now

Protein-family referee on chr16, plus NPIP (Soto and U2), against the shipped de novo baseline
(F 0.214 referee, NPIP Soto 0.727 at sensitivity 0.600).

| outcome | verdict |
|---|---|
| referee F **+0.05 or more**, or NPIP sensitivity 0.600 → **≥0.750** (guided's level) | ⭐ **THE FP HALF IS WORTH PURSUING** — a real rule is worth its FN cost |
| +0.02 to +0.05 | ⚠ **MARGINAL** |
| **< +0.02**, or NPIP sensitivity does not move | ⛔ **THE FP HALF IS CLOSED** — perfect repair buys too little to justify any lossy approximation of it |

**Predicted, before looking — ⭐ on NPIP, ⚠ on the referee.** r1015 showed the FP cost is concentrated on
exactly the thesis family (`NPIPB4`→`LOC112268174`, `NPIPB5`→`SMG1P1`, `NPIPB12`→`SMG1P2`), so a perfect
repair should lift NPIP sensitivity toward guided's 0.750 while moving the chromosome-wide referee much
less, since only 28 of 306 referee genes are affected.

⚠ Splitting a locus also costs each piece its own alignments (§6w6's mechanism), so the oracle is not a
pure gain and the measured ceiling may be **below** the naive 9.2%. That is itself worth knowing.

I will not change the arm, the metrics or the bar after seeing any number.

---

# OUTCOME (2026-09-22) — ⛔ **THE FP HALF IS CLOSED. Perfect oracle repair buys +0.021 and makes NPIP WORSE.**

2,550 de novo loci → **2,885 oracle-split nodes** (233 fused loci expanded), same
`minimap2 -x asm20 -c --eqx -P -t 4` as the baseline, same `mcl_families` config.

| arm | referee F | referee prec | collapsed | NPIP Soto F | NPIP sens | NPIP U2 F |
|---|---|---|---|---|---|---|
| A0 shipped | 0.214 | 0.949 | 12 | 0.727 | 0.600 | 0.610 |
| **`--min-cov-shorter 0.70`** | **0.230** | 0.952 | 12 | **0.750** | 0.600 | 0.621 |
| **ORACLE perfect split** | **0.235** | 0.953 | **5** | **0.687** | **0.550** | 0.632 |

⛔**Verdict per the bar: NO.** Referee F +0.021 (the ⚠ band), and **NPIP sensitivity does not rise — it
FALLS, 0.600 → 0.550**, which the bar names as a ⛔ condition outright.

⭐⭐**The decisive comparison: the shipped flag delivers +0.016 of the oracle's +0.021 = 76% of the ENTIRE
achievable gain from perfect node-FP repair** — with no oracle, no annotation, and while *improving* NPIP
(0.750 vs the oracle's 0.687). **There is at most +0.005 F left in the false-positive half, and reaching
it costs the thesis family.**

## Why perfect repair loses on NPIP — §6w6's mechanism, now confirmed on the FP side

The pre-registration flagged it and the measurement confirms it: *"splitting a locus also costs each piece
its own alignments."* The oracle correctly recovers `NPIPB4`, `NPIPB5`, `NPIPB12` as their own nodes
(collapsed 12 → **5**, the cleanest collapse count of any arm), but **the pieces then fail the coverage
gate that the fused locus passed**, so NPIP sensitivity drops. ⭐**Over-merge is simultaneously the cause
of NPIP's collapse AND the reason its members are admitted at all** — repairing it trades one for the
other, which is why no operator in this line ever came out ahead.

## Prediction scorecard — wrong on the half I was most confident about

I predicted **⭐ on NPIP**: r1015 had shown the FP cost concentrated on exactly that family
(`NPIPB4`→`LOC112268174`, `NPIPB5`→`SMG1P1`, `NPIPB12`→`SMG1P2`), so perfect repair "should lift NPIP
sensitivity toward guided's 0.750." It fell to 0.550. ⭐**Knowing which genes an error destroys does not
tell you that removing the error recovers them** — the same loci carried the edges. That inference gap is
the single most expensive one in this whole line, and it is why four operators were built before anyone
measured the ceiling.

⚠ Also caught mid-run: I first built the oracle PAF with `-N 50 -p 0.1 -X --secondary=yes` (48,176 records)
when the baseline used `-x asm20 -c --eqx -P` (127,773). **An arm is not comparable until its aligner
invocation is copied from the baseline's own log**, not from the config a different arm happened to use.

> **Scorer (2026-09-22 port):** every `bench/mode_family_score.py` number above is reproduced byte-for-byte by `target/release/family_score` (same flags; 732/732 parity, scipy's assignment tie-breaking included — r1045/r1046). The Python was retired in §6z2.
