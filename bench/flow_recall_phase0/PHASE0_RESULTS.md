# Phase 0 flow-recall diagnostic — RESULTS

**Date:** 2026-06-09
**Spec:** `docs/superpowers/specs/2026-06-09-parse-trflong-flow-recall-scope.md`
**Plan:** `docs/superpowers/plans/2026-06-09-flow-recall-phase0.md`
**Input:** 668 st_only misses (annotated isoforms ST recovers, rustle misses), all 668 captured
cross-tool (baseline `-L`, per-locus slices). No OOM (18 GB free throughout).

## Headline: literal gate = PROCEED, but the honest verdict for the parse_trflong work is STOP/REDIRECT

The `gate_report` script printed **PROCEED** (ceiling 668 ≥ 150, discriminator AUC 0.78).
**That verdict is overstated** — it used the full 668 st_only as the recall ceiling, which
conflates three different mechanisms. The attribution shows the `parse_trflong` flow-generation
work targets only **59** of them.

## (a) Generation census — 416 generated / 252 non-generated

Of the 668 misses, **62% (416) ARE generated** by rustle's flow (the exact ref chain appears
in `path_extracted`) and then die in a downstream filter. Only **38% (252) are true
non-generation.** (The container-only probe had estimated ~50% non-gen; genome-wide it's 38%.)

| category | generated | non_generated |
|---|---:|---:|
| partial_multi | 188 | 80 |
| altsplice_1junc | 100 | 41 |
| no_overlap | 63 | 11 |
| container | 43 | 43 |
| contained | 16 | 12 |
| no_shared_introns | 6 | 5 |
| single_exon | 0 | 59 |
| partial_disjoint_1junc | 0 | 1 |

(single_exon is entirely non-generated, as expected — rustle emits no single-exon tx.)

## (b) Recall-oracle ceiling + precision cost

- Current rustle FSM (genome-wide): **24,373**. Ceiling if rustle adopted ST's extraction
  exactly: **25,041 (+668)**.
- **Precision cost: 1,713** ST-extracted non-annotated chains at the miss loci would be added
  as FPs. **Crude recall:FP ratio = 668:1713 = 0.39** — naively matching ST's extraction adds
  ~2.6 spurious chains per real chain recovered. Net-F1-negative without a discriminator.

## (c) Separability — moderate discriminator exists

Among ST-extracted chains rustle misses (real=316 annotated, spurious=2152 unannotated):

| feature | AUC(real > spurious) |
|---|---:|
| cov | **0.782** |
| longcov | 0.765 |
| entry_abund | 0.765 |
| nexons | 0.575 |

Coverage separates real from spurious at **AUC ~0.78** — a *moderate* discriminator (clears
the 0.70 bar but is far from clean; any threshold still mixes substantial FPs). This is more
hopeful than the prior precision-direction finding (which found no discriminator via
min_jct_mm), but 0.78 is not strong enough to make naive adoption net-positive on its own.

## (d) Attribution of the 252 non-generated — THE decisive number

| cause | n | % of non-gen |
|---|---:|---:|
| **graph_missing** (a ref intron absent from rustle's junction_accept) | 193 | 77% |
| **flow_enumeration** (junctions present, flow didn't walk the chain) | 59 | 23% |

**Only the 59 flow_enumeration cases are the `parse_trflong` seed/extension target.** The 193
graph_missing cases are a *different, deeper* lever (junction acceptance / graph construction).

## Synthesis — the 668 gap decomposes into THREE mechanisms, not one

| mechanism | n | % | lever | parse_trflong? |
|---|---:|---:|---|---|
| generated-then-filtered | 416 | 62% | downstream filters (retained_intron / isofrac) | NO |
| graph_missing | 193 | 29% | junction acceptance / graph construction | NO |
| flow_enumeration | 59 | 9% | parse_trflong seed_order_st / update_abundance_st | **YES** |

## Verdict for the parse_trflong flow work: REDIRECT (do not commit the rewrite)

- The `parse_trflong`-specific ceiling is **~59 chains (9% of the gap)** — **below** the 150
  materiality bar the spec pre-registered. The literal PROCEED was an artifact of the gate
  using 668 (total st_only) instead of 59 (flow-attributable) as the ceiling.
- Completing `seed_order_st` + `update_abundance_st` is a high-risk core-flow rewrite for a
  ~59-chain ceiling. Not worth it on these numbers.

**What the diagnostic redirected us toward (the real recall mass):**
1. **Generated-then-filtered (416, 62%)** — the largest lever, and unlike chr19 it now has a
   *moderate* discriminator (cov AUC 0.78 on this population). This is the filter side
   (retained_intron/isofrac) we found hard before — but the 0.78 signal is worth a fresh,
   recall-framed look: can a cov-gated rescue of generated-but-killed annotated chains net out
   F1-positive? Most tractable next probe.
2. **graph_missing (193, 29%)** — a newly-surfaced lever: rustle's junction acceptance / graph
   construction drops junctions ST keeps. Distinct from both filters and parse_trflong; its own
   precision risk (accepting more junctions → more FPs). Worth a separate characterization.
3. **flow_enumeration (59, 9%)** — the parse_trflong target; smallest. Park unless 1 & 2 are
   exhausted.

The diagnostic succeeded: it prevented committing a high-risk flow rewrite that addresses only
9% of the recall gap, and re-pointed at two larger, better-characterized levers.

## ⚠ CORRECTION (2026-06-09): attribution had a junction-coordinate bug — parse_trflong verdict REVERSED

The Phase 0 attribution (d) used the raw `junction_accept` (start,end) to test whether a ref
intron is in rustle's graph. But `junction_accept` logs the **exon-flanking** convention
`(intron_first-1, intron_last+1)`, NOT the intron-interior coords `ref_chain` emits. The raw
comparison NEVER matched a multi-exon intron → every multi-exon non-generated miss fell to
`graph_missing`, and the only "flow_enumeration" hits were the **vacuous 0-intron single-exon**
cases. So the reported "graph_missing 193 / flow_enumeration 59" was an artifact.

**Corrected attribution of the 252 non-generated** (`attribute.py` fixed to match both coord
conventions + separate single-exon; verified end-to-end on NC_073224.2 XM_055373957.2: all 6
junctions accepted, chain absent from path_extracted):

| corrected cause | n | was (buggy) |
|---|---:|---|
| **flow_enumeration** (junctions in graph, flow didn't walk) | **190 (28% of 668)** | 59 |
| single_exon (0 introns) | 59 | (mislabeled flow_enum) |
| graph_missing (junction truly absent) | 3 | 193 |

**This REVERSES the parse_trflong REDIRECT.** The flow-path-enumeration target is **190 chains
(28% of the gap)** — comfortably ABOVE the 150 materiality bar, not 59 (9%). The gate's
recall-ceiling condition now PASSES. `graph_missing` is negligible (3), NOT a lever.

**Corrected 668 landscape:** 416 generated-then-filtered (62%) · 190 flow_enumeration (28%,
parse_trflong) · 59 single_exon (9%) · 3 graph_missing.

**Status of parse_trflong:** re-opened as a material lever (190). The recall ceiling passes; the
OPEN question is the net-F1 feasibility — the precision cost of *generating* those 190 (how many
spurious chains the completed seed/extension would also produce). That generation-precision-cost
was never measured (the recall-oracle +1713-FP figure was for adopting ST's FULL extraction, a
superset). That measurement is the real Phase-0.5 gate before committing the rewrite.

`graph_missing_anatomy.py` is superseded (built on the buggy 193). Corrected tooling:
`attribution_fixed.py`, fixed `attribute.py`.

## FOLLOW-UP (2026-06-09): filter-rescue lever (416) — FALSIFIED

The Phase 0 separability AUC 0.78 was measured on the WRONG population (ST-extracted chains
rustle misses). The filter-rescue operation actually acts on **rustle's OWN generated-then-killed
chains**. `rescue_feasibility.py` measured the correct population (per locus: rustle
path_extracted = generated; rustle r.gtf = kept; killed = generated − kept; label real if
annotated):

- real-killed (annotated, want to rescue): **680** · spurious-killed (correctly killed): **5855**
- **cov AUC(real > spurious) = 0.568** — essentially random. cov does NOT separate rustle's real
  killed chains from spurious ones. (The 0.78 was an artifact of ST's extraction quality.)
- cov-threshold sweep is net-F1-NEGATIVE everywhere: t≥3 rescues 147 real but adds 495 FP
  (1:3.3); best ratio (t=50) rescues only 3.

**Filter-rescue is DEAD** — same indistinguishability rock as chr19 / pathpat_flow_parity_scope,
now confirmed genome-wide on the correct population. Two of the three redirect levers
(parse_trflong, filter-rescue) are now no-gos. Only `graph_missing` (193) is untested — but the
recurring pattern (real/spurious indistinguishable by available features) predicts it is also a
Sn/Pr reshape. The one orthogonal, non-reshape recall win remains strand-aware bundling
(`RUSTLE_STRAND_PURE_MINORITY`).

Artifacts: `gen_census.jsonl`, `attribution.jsonl`, `separability_rows.jsonl`,
`rescue_feasibility.py`; harness in this directory; cache under `cache/` (gitignored).
