# Unified gene-conversion vs RT/template-switch discriminator

**Why.** Both gene conversion and RT/template switching produce the *same observable* in copy-assignment:
one molecule (or copy) whose PSV-allele pattern is a **mosaic of two copies** — a cross-copy recombinant.
This is the object that breaks Strong Separation / sits on the K≥3 recombination frontier
(`copy_assignment_theory.md` §5). But **one is biology to report, the other is a library artifact to
filter**, and they look identical per-read. The mosaic detector previously confirmed a conversion on
**recurrence alone** (`ConversionEvent.confirmed = recurs ≥3 molecules, tight dispersion`) — and that is
*insufficient*: a sequence-driven **microhomology hotspot** makes RT-switches *recur*, so they pass the
recurrence gate. The fix adds two orthogonal legs.

## The three-leg discriminator (`mosaic.rs::classify_event`)

| leg | signal | meaning |
|---|---|---|
| **recurrence** | `ev.confirmed` (already there) | breakpoint recurs across independent molecules |
| **microhomology** | `genome::breakpoint_microhomology(chrom, br_lo, br_hi, 6..12)` | a direct repeat flanking the breakpoint = the RT/template-switch signature |
| **DNA support** | catalog lookup (`Option<bool>`) | gene conversion is heritable → in matched DNA; RT switch is RNA-only |

Truth table (evaluation order matters — the artifact rule fires *before* the conversion rule):

```
microhomology && !dna_present                    -> RtSwitchArtifact   // template signature, not DNA-rescued
ev.confirmed  && !microhomology && !dna_absent   -> GeneConversion     // recurrent, no signature, DNA not contradicting
!ev.confirmed && !microhomology                  -> ChimeraSuspect     // sporadic one-off
otherwise                                        -> Ambiguous          // mh∧dna conflict, or recurrent-but-DNA-absent
```

Key semantics: **DNA is a veto, not a requirement.** `dna_supported = None` (unchecked) does **not**
block a `GeneConversion` call (so the two cheap legs ship without the DNA catalog wired); `Some(false)`
(checked and ABSENT from the genome — contradicts heritability) downgrades to `Ambiguous`. We never
upgrade to `GeneConversion` on a microhomology breakpoint.

## What's wired

- **`mosaic.rs`** — `Classification` enum + pure `classify_event` (6 unit tests on the truth table).
- **`genome.rs`** — `breakpoint_microhomology` (scans `is_rt_switch`'s direct-repeat test over a k-range;
  widens past the old fixed 8 bp) with a **low-complexity guard**: a homopolymer / dinucleotide-repeat
  window (`AAAAAA`, `ATATAT`) matches a direct repeat almost everywhere, so without a filter a real
  conversion near a simple repeat would be wrongly demoted to an artifact (the error direction
  *suppresses real conversions*). The matched window must carry ≥3 distinct bases. 2 unit tests.
- **`copy_assign_pipeline.rs`** — `FamilyDetail.conversion_class` + `classify_conversions(detail, genome,
  chrom, dna_closure)`.
- **`denovo_pipeline.rs`** (the de-novo family / `copy_assign` path) — classifies each confirmed event
  with the **live microhomology leg** and reports `… N confirmed -> G gene_conversion, R rt_switch_artifact`.
  DNA leg passed as `None` here (catalog leg = follow-up).
- **`ConversionEvent.chrom`** — each event now carries its chromosome (set in `aggregate_family`, which
  also makes chrom part of the **cluster key**: breakpoints at coincidentally-similar positions on
  DIFFERENT chromosomes — multi-chrom paralog families like RABL2A/RABL2B — no longer wrongly merge).
- **production `--vg` emit path** (`pipeline.rs`, `RUSTLE_VG_MOSAIC_EMIT`) — now **live**: when a genome
  is loaded, an event whose breakpoint carries the microhomology signature is classified
  `RtSwitchArtifact` and **suppressed from emission** (not promoted to a gene-conversion isoform);
  `GeneConversion`/`Ambiguous` still emit. Recurrence-alone could not tell these apart; the genome
  microhomology leg can. Reports `[VG-MOSAIC-EMIT] suppressed N RT/template-switch artifact event(s)`.

Default-off preserved: the whole mosaic pass is opt-in (`RUSTLE_VG_MOSAIC_ON`); `classify_conversions`
only runs when events exist, and `ConversionEvent.confirmed` is untouched — so the discriminator is
purely additive metadata.

## Measurement — ground-truth confusion matrix (both directions, real code path)

`mosaic::tests::ground_truth_conversion_vs_rt_switch_confusion_matrix` drives the **full real path**
(`detect_mosaic → aggregate_family → genome.breakpoint_microhomology → classify_event`) over a
constructed genome with two recurrent recombinant loci:

| planted locus | breakpoint sequence | classified as |
|---|---|---|
| gene conversion | flanks **differ** (no direct repeat) | **GeneConversion** ✓ |
| RT/template switch | breakpoint at an **exact 8 bp direct repeat** | **RtSwitchArtifact** ✓ |
| (same conversion, DNA absent) | — | **Ambiguous** (DNA veto) ✓ |

The RT-switch case is the load-bearing one: it is **recurrent** (5 molecules → `confirmed = true`), so
recurrence alone would have mis-called it a gene conversion; the microhomology leg correctly reclassifies
it as an artifact. 627-test suite green.

## Adversarial review

A 3-lens review (truth-table soundness / additive-off / microhomology correctness) confirmed
**additive-off holds** and the **truth table is sound**, with **zero blockers**. It flagged (and the
synthesis rejected) a proposed "abort if the breakpoint bracket spans >1 kb" guard as a misread —
`is_rt_switch` only reads two short windows *ending at* each endpoint, never the span between, so an
intron between the two flanking PSVs does not corrupt the test; the guard would have suppressed real
signal. The one real, non-blocking issue it surfaced — low-complexity windows over-calling microhomology
(direction: suppresses real conversions) — is fixed by the ≥3-distinct-base entropy guard above.

## Status / follow-ups

- **`copy_assign` / `detect_families` path: live** (microhomology leg wired, classification reported).
- **Production `--vg` GTF-emit path: live** — `ConversionEvent` now carries `chrom`, so the emit stage
  runs the microhomology discriminator and suppresses RT/template-switch artifacts from gene-conversion
  isoform emission. Default-off preserved (gated on `RUSTLE_VG_MOSAIC_ON` + `RUSTLE_VG_MOSAIC_EMIT`;
  byte-identical no-op verified on chr19: MOSAIC on vs off = 0 structural diff).
- ⚠ **Two surfaces still label on recurrence alone** (the discriminator gates *emission*, not these
  labels): the `[VG-MOSAIC]` stderr *report* line in `detect_and_report_mosaics`, and the GTF **attribute**
  `gene_conversion "confirmed"` written in `gtf.rs` — both read `ConversionEvent.confirmed` (recurrence),
  not the microhomology classification. So a recurrent RT-switch is *suppressed from isoform emission* but
  would still be *tagged* `confirmed` on its host transcript if `MOSAIC_EMIT` is off. To make the labels
  consistent, thread `classify_event` into the `family_verdict.conversion` so the attribute can read
  `rt_switch_artifact` / `ambiguous` (loose end L21). The actual emitted gene-conversion *isoforms* are
  already correctly filtered; only the descriptive attribute/report lag.

## DNA heritability leg — measured (prototype), NOT wired as a veto

`dna_support.py` prototypes the DNA leg against the T2T PSV catalog
(`/home/juanfra/winloci_scratch/dna_catalog/DNFAM*.json`) — measure-first before any Rust port. The
catalog is **ref0-centric**: `matrix[label][pos] = base`, `pos = genomic − ref0_start`, sparse
(absent = matches ref0). A heritable gene conversion ⇒ a DNA copy whose allele vector is itself a
**mosaic** of two others (verified convention against the FASTA).

Findings:
- **(A) The signal is real:** **42% of multi-copy families (17/40)** carry a DNA mosaic
  (historical-conversion signature) — e.g. `DNFAM39 L0 = L2 | L3`. So DNA *can* corroborate a conversion.
- **(B) Coverage is low:** ref0 intervals cover only **2.88% of the genome**; and within a family only
  breakpoints near the *ref0* copy are projectable (non-ref0 copies would need re-alignment).
- **(C)/(D) "absent" is unreliable:** DNA mosaics are *localized*, so families that HAVE one return
  "absent" at almost every position. Sparse + ref0-centric + position-specific.

**Conclusion (this is why measure-first mattered):** catalog **"absent" is weak negative evidence** —
wiring it as a **veto** (`Some(false)` → downgrade GeneConversion → Ambiguous) would **wrongly suppress
real conversions**, the exact failure direction we guard against. So the catalog DNA leg must be
**positive-only**: return `Some(true)` (a DNA mosaic coincides with the breakpoint) or `None`, **never
`Some(false)`**. `classify_event`'s truth table stays general (it can still veto a *reliable* absence
source); the catalog-backed closure is restricted to positive corroboration. The prototype
(`dna_support()`) now returns only `supported` / `unchecked`.

Given the low coverage and coordinate fragility, the DNA catalog is best used as an **offline
corroboration** (annotate confirmed conversions that coincide with a known DNA mosaic as
"DNA-corroborated"), not an inline Rust catalog loader. The **microhomology leg remains the primary,
shipped discriminator**; the DNA leg adds rare but strong positive confirmation where it applies.
