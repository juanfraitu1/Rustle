# PREREG — the conjoined read-through as an O1 object (2026-09-07)

**Written before the code.** md5 in `mcl_ann/adj/readthrough/PREREG.md5`.

## Why
Of the 71 NPIP reads rejected for running past their locus (§6fp, register row 723): **48 are real
read-throughs** — canonical splice sites, the intron shared by ≥ 3 molecules, the far end landing in another
catalog unit; 40 of those are one event, NPIP unit 2 → MCL27:0 through a single 15-kb intron. 15 are chaining
artefacts and 8 are terminal exons the cross-family clip severs. **The majority is biology the catalog has no
object for**, so it arrives as rejection. The alternative considered and refused today was feeding annotated
junctions to the aligner (`--junc-bed`): it injects the annotation into the read layer, biases near-ties toward
better-annotated copies, and suppresses exactly what O3 looks for.

## The object
`mcl_families --emit-readthrough-units` (**default OFF**; with it off every output byte is unchanged).
A read-through unit is emitted for an ordered pair of emitted units (A, B) on one contig when

1. ≥ `--min-reads` distinct primary molecules each have an aligned block in A's emitted chain **and** an
   intron whose far end falls inside B's emitted chain, and
2. those molecules share **one and the same** intron (identical coordinates), and
3. that intron is **canonical** — `GT..AG`, `GC..AG` or `AT..AC` read on the unit's strand.

Row: `member_status = readthrough`, `source = readthrough`, exon chain = A's chain ∪ B's chain, span from A's
first exon to B's last, `n_reads` = the joining molecule count. It is a unit, not a member: `bench/o1_eval.py`
counts it with the candidates, never in the specificity denominator. `copy_assign` sees it as an ordinary
candidate, so a molecule spanning both loci can be assigned to it instead of rejected.

## Predictions (gorilla 3-contig, `rna_units_v12` inputs)
| # | prediction |
|---|---|
| **P0** | with the flag off, `units.tsv` is **byte-identical** to the current build |
| **P1** | the NPIP unit 2 → MCL27:0 event is emitted as exactly **one** read-through unit |
| **P2** | on the three contigs, **fewer than 40** read-through units in total — this is a rare object, not a second catalog |
| **P3** | O1 metrics on the 26 LCR16a truth loci are **unchanged** (sensitivity 26/26, specificity 25/25): read-through units enter as candidates, never as members |
| **P4** | rerunning O2 on NPIP with the read-through unit present, **≥ 30 of the 48 real read-through reads** stop being rejected (they are assigned to the read-through unit or abstain with it as a candidate) |
| **P5** | no molecule already assigned to a copy changes its assignment because of the new unit — the object absorbs rejects, it does not move accepted reads |

## Interpretation fixed in advance
- P4 holding with P5 holding ⟹ the object is the right shape and ships behind its flag for the user's decision.
- P5 failing ⟹ the read-through unit competes with real copies and must NOT ship as a candidate; it would then
  be a reporting-only row.
- P2 failing ⟹ the ≥ 3-molecule canonical-intron rule is too loose; report the count and do not tune it silently.
