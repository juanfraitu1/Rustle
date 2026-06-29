# Would resolving multimapper ambiguity improve families and their isoforms?

*Bridging advisor interests #1 (family definition) and #2 (copy assignment): the achievable improvement from
resolution is bounded by RNA-level identifiability, and on the GGO panel the bound bites to ~0.*

Measured on `GGO.bam` over a 7-family panel + the chrY mega-array (workflow `wf_ad5ebaa7-0fb`, 4 measurements
+ synthesis). The answer separates into two questions that resolve differently.

## Families: detection unchanged, completeness +3 copies

- **Detection (membership) is unimproved — by construction.** The conflict graph is built from the ambiguity
  itself; resolution decides *where* an ambiguous read goes, not *which* reads are ambiguous. Family membership
  is defined prior to resolution.
- **Completeness improves via starved-copy recovery, but modestly.** A copy that loses the arbitrary MAPQ-0
  primary-flag battle to a tied sibling has ~no primary reads and is invisible to primary-only assembly.
  Reclaiming its secondaries (= resolution) recovers it. Found **3 starved copies** (gate: primary < 3, secondary
  ≥ 3): `LOC115934278` (1 prim / 23 sec; sibling AK6 hoards 848 prim) + 2 chrY-array loci
  (`22450601` 2/57, `22474453` 2/91). The 2 chrY copies are the strongest case — their reclaimed secondaries are
  100% MAPQ-0, ~90% spliced, with the arbitrary primary parked at a tied in-array sibling, so resolution recovers
  **both the copy and its multi-intron isoform**. (Lower bound: gap-clustering understates the ≥12-copy array as
  8 loci, and minimap2's N=5 cap drops secondaries past the 6th placement.)

## Isoforms: two clean regimes

| regime | families | ambiguous (primary-MAPQ-0) mass | resolution corrects |
|---|---|---|---|
| **separate-chrom / well-separated** | RABL2, AK6, CCDC196, MAGEA_p2 | 0–12% | **~0** — minimap2's `de` already de-ties them (RABL2: 0/195 MAPQ-0, primary on the better-fit copy 195/195) |
| **co-located tandem** | MAGEA_p1/p3/p4 | **95–99%** | the coin-flip label on **474/495 spliced reads** |

In the co-located MAGEA pairs, 95–98% of family reads are spliced **and** primary-MAPQ-0 **and** exactly
de-tied (de-gap = 0.0000) — nearly every isoform is pinned to a coin-flip copy. That is the headroom resolution
would, in principle, correct.

## The bound — and why it is the whole story

**0 of the panel's 494 ambiguous reads lie in a K≥2 (PSV-bearing) family.** The two levers are **disjoint**:

- The ambiguous mass that **needs** resolving (488/494 = 98.8% MAGEA co-located) is **K=0**: per-read edit
  distance is *identical* against both copies (NM_A == NM_B, 1–3 mismatches = HiFi error floor; 0/311 reads
  differ between copies on MAGEA_p3). The copies are **sequence-identical over the transcribed exonic region** —
  the 67–72% genomic divergence is entirely intronic/intergenic. No PSV column can adjudicate; attribution is
  intrinsically arbitrary.
- The copies that **can** be resolved (chrY AS-gap median 259; RABL2; AK6; CCDC196) are already MAPQ>0 → already
  correctly placed by minimap2. No ambiguity remains to resolve.

> **Where resolution is possible it is unnecessary; where it is needed it is impossible.**
> Achievable-improvement ceiling on this panel's ambiguous mass = **0%** (resolvable 0/494; stuck 494/494).

## What this means (not a defeat — a characterization)

This is the identifiability theorem made concrete, and a precise statement of **when copy-resolution helps**:
**co-located AND exonically-divergent (K≥2)** families. The GGO panel contains no such family — co-located ⇒
exonically identical (K=0), divergent ⇒ already separated (MAPQ>0) — so the two conditions never co-occur here.

- The improvement is real and was demonstrated in the regime it targets: the **sim5x K-ladder** (synthetic
  co-located copies, PSV ladder K=0..8) shows K≥2 → 100% correct assignment. The win exists; it is gated by K.
- The value for the thesis is the **bound itself**: the achievable gain from resolution is a measurable function
  of RNA-level identifiability (K = exonic PSV count), predictable a priori per family. The same K governs
  detection, resolution, and isoform separation — the through-line to the identifiability theorem.
- The **joint-inference loop** (detect → resolve → recover copies/isoforms → refine loci → re-detect) is virtuous
  only in the co-located-AND-divergent regime; on GGO it stalls at the resolve step because the reclaimable mass
  is K=0. The one place it closes (the 2 chrY starved copies) is granularity- and cap-limited — a trickle, not a
  cascade.

**Bottom line.** Yes, in principle: family completeness gains 3 starved copies, and 474/495 co-located spliced
reads carry a coin-flip label resolution would fix. But the achievable gain on this real panel is ~0%, because
the ambiguity that needs resolving is RNA-level-unresolvable (K=0) and the resolvable copies are already
resolved. A real win requires a co-located-AND-divergent family — present on the sim ladder, absent from this
GGO panel.
