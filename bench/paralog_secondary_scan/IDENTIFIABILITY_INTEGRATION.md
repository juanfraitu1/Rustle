# Including minimizer-identifiability ("malleability") in rustle

Reference-only score predicting whether a paralog copy's reads are *assignable*
(separable copies) or at the *sequence floor* (tied / MAPQ-0 / non-identifiable).
This is the principled gate that distinguishes real RABL2-type copy wins from
DAZ3-type phantoms — computable BEFORE looking at a single read.

## The score (validated heuristic)
`disc_frac(copy) = |minimizers(copy) not in ANY sibling copy| / |minimizers(copy)|`
using rustle's exact scheme (FNV-1a, k=15, w=10, fwd). Anchors: RABL2 ~0.58 → 0% tied
(IDENTIFIABLE); DAZ ~0.07 over a 2.9kb identical core → ~19% tied (SEQUENCE_FLOOR).
Standalone scorer: `minimizer_identifiability.py` (family.fasta | --genome+--bed).

## Why it's nearly free in rustle
rustle ALREADY computes `vg_family::family_graph::minimizers(seq,15,10)` and
`minimizer_jaccard` when pairing homologous exons across copies (vg.rs:3709, the
0.30 merge-Jaccard). The same minimizer sets give disc_frac at ~zero extra cost.

## Three integration points
1. **Compute at family build** (vg.rs, where copies/exons are already paired):
   for each copy, disc_frac vs the union of sibling-copy minimizers. Store on the
   copy/CopyCertificate struct alongside the existing own_ev/strict/tied fields.

2. **Reference-prior to the decisive gate** (`RUSTLE_VG_DECISIVE_GATE`, pipeline.rs):
   today the gate is purely empirical (own_ev>=2, strict-frac>=0.5, primary-cov escape).
   Add the prior:
     - `SEQUENCE_FLOOR` copy  -> require MUCH stronger read evidence to emit, or suppress
       (the reads CANNOT carry assignment signal — emitting is fabrication risk; this is
        exactly the DAZ3 phantom failure mode the gate keeps tripping on).
     - `IDENTIFIABLE` copy    -> safe to emit on weaker read evidence (assignment is trustworthy).
   Turns the gate from "tuned thresholds" into "evidence required scaled by what the
   sequence can possibly support."

3. **Output annotation** (the "include the information" ask):
   - GTF: per-emitted-copy attribute `identifiability "0.61"; copy_tier "IDENTIFIABLE";`
   - rescue.tsv (`RUSTLE_VG_RESCUE_REPORT`): new columns `disc_frac` + `tier` next to
     the existing decisive-evidence certificate. Reviewer instantly sees which rescued
     copies are real vs floor.

## How it answers "more copies"
It does NOT fabricate copies. It tells you WHICH missed copies are worth chasing
(IDENTIFIABLE/BORDERLINE — genuine recall headroom) vs which are at the floor
(SEQUENCE_FLOOR — recovering them = inventing transcripts on tied reads). Pairs with
the decisive gate to raise copy *recall* without cratering *precision*.

## Self-test (synthetic, reproduces anchors)
| copy | divergence | disc_frac | tier |
|---|---|---|---|
| near-identical A | 0.2% | 0.013 | SEQUENCE_FLOOR |
| near-identical B | 0.2% | 0.025 | SEQUENCE_FLOOR |
| divergent C | 8% | 0.609 | IDENTIFIABLE |
