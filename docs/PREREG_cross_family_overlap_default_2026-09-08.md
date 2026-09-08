# PREREG — flipping the cross-family exon-overlap rule ON by default (2026-09-08)

**Written before the flip.** md5 in `mcl_ann/adj/xfdefault/PREREG.md5`.

## The change
`mcl_families --no-cross-family-exon-overlap` becomes the **default**, with `--allow-cross-family-exon-overlap`
as the escape that reproduces every catalog built before today. The rule: two units of DIFFERENT families may
not claim the same exon bases; contested bases go to the unit whose own annotated member span contains them,
else to the unit with more reads. (The own-donor re-check on read-through junctions is already unconditional
and is not part of this flip.)

**Why.** Cross-family overlap manufactures read-throughs out of ordinary introns: `MCL1:3`'s chain overran
2 kb into `MCL27:0`, whose FIRST EXON is that overlap, so every read splicing across `LOC129527585`'s own
first intron looked like a molecule joining two loci. **32 of 42 guarded read-throughs were annotated introns
of one gene** (§6ge, register 754).

## Measured before the flip (§6gm), on the current OFF/ON pair
| | gorilla 3-contig | human Soto slice |
|---|---|---|
| units trimmed | 31 | 41 |
| unit rows changed (cols 1–9) | 60 of ~1,000, 8 families | — |
| NPIP sensitivity / specificity | 26/26, 25/25 unchanged | — |
| Soto sensitivity / specificity | — | 0.909 unchanged; 0.649 → 0.663 |
| Soto size in-band | — | 0.82 → 0.81 (−3 of 319) |
| gorilla read-throughs | 42 → 14; annotated introns among them 32 → 1 | — |

## Predictions for the flipped build
| # | prediction |
|---|---|
| **P1** | with `--allow-cross-family-exon-overlap`, output is **byte-identical** to the current default. If not, the flip is not a flip and must be reverted |
| **P2** | gorilla NPIP: sensitivity **26/26**, specificity **25/25** — unchanged |
| **P3** | human Soto slice: sensitivity **0.909** unchanged, specificity **≥ 0.649** (not worse than OFF) |
| **P4** | human chr16+18 NPIP: sensitivity **26/26**, specificity **≥ 0.963** — this substrate has NOT been tested with the rule and is the one that could surprise |
| **P5** | gorilla read-throughs **≤ 20**, of which **≤ 3** are annotated introns |

## Interpretation fixed in advance
- Any prediction failing ⟹ the default stays OFF and the failure is reported, not tuned around.
- P4 is the real test: the other three substrates were used to design the rule, chr16+18 was not.
