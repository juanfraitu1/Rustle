# Pre-registration — a CORROBORATED union rep: widening only where the evidence is real

**Written 2026-09-22, §6z0, before any arm is run.** User: *"can we try to implement that or implement it
for the cases that actually need it?"*

## The prior this must beat, and why it is not a re-proposal

**Register 303 refuted the plain union rep** (`RUSTLE_LOCUS_EXON_UNION=1`): Soto recall **65.5% → 44.8%**,
pure families **135 → 72**, with the cause named — *"widened reps inflate their own coverage denominator
and fall below the 0.50 floor."*

⭐**§6y9/r1040 says WHY the union was the wrong target.** On chr20 (505 multi-transcript loci), of the
3,483 exons a non-rep transcript contributes, **63.5% appear in exactly ONE transcript**. The existing
`RUSTLE_LOCUS_UNION_MIN_READS` floor (default 3) does not separate them — it admits **74.8%**, because a
single chain can carry many reads. So r303 widened reps with a set that is **mostly singleton noise**.

**The change under test** (`RUSTLE_LOCUS_UNION_MIN_TX`, default 0 = off, byte-identical): filter per EXON
rather than per transcript — a merged interval survives only if ≥ k member transcripts place an exon on
it. The rep's own exons are always kept, so a widened rep is never smaller than today's single-chain rep.

## Arms

| arm | flags |
|---|---|
| **B** baseline | single-chain rep (shipped) |
| **U** plain union (r303's refuted arm) | `RUSTLE_LOCUS_EXON_UNION=1` |
| **C2** corroborated | `RUSTLE_LOCUS_EXON_UNION=1 RUSTLE_LOCUS_UNION_MIN_TX=2` |
| **C3** stricter | `RUSTLE_LOCUS_EXON_UNION=1 RUSTLE_LOCUS_UNION_MIN_TX=3` |

⚠**U is mandatory.** Without it I cannot tell whether C2 beat the baseline or merely beat r303's noise;
r303's numbers are from a different pipeline generation and are NOT a valid comparator.

## Metrics · substrate · bar

Per arm: rep exonic bp (did widening happen at all) · **graph nodes / edges** (r303's failure mode is a
coverage-denominator collapse, which shows up as edges lost) · family pairwise precision / recall / F
against the protein referee and Soto.

Development **chr20** (the substrate r1040 measured, and the only one with a de novo GTF on disk).
Held out **chr16** — a different chromosome and the one every other arm this session used.

| outcome | verdict |
|---|---|
| C2 or C3 beats **B** on held-out F, and beats **U** | ⭐ **ADOPT** — corroboration is what r303's union lacked |
| C2/C3 beats U but not B | ⚠ **PARTIAL** — the diagnosis was right, the remedy still loses |
| C2/C3 ≤ U, or edges collapse as in r303 | ⛔ **NO** — widening fails for a reason corroboration does not fix |

**Predicted, before looking — ⚠ PARTIAL.** r1040's diagnosis is well measured, so C2 should clearly beat U.
But r303's mechanism is a **denominator** effect, and a corroborated union still enlarges the rep — just
less. §6x3/r1002 showed the same denominator is what evicts 679 chr16 loci. I expect C2 to recover most of
U's loss without clearing the baseline. ⚠If C2 does clear B, the natural follow-up is C2 **plus**
`--min-cov-shorter`, which exists precisely to absorb denominator inflation — but r911 warns gains are not
additive, so that is a separate arm, not a claim.

I will not change the arms, the substrates, the metrics or the bar after seeing any number.
