# Ideal-scenario chromosome simulation — in theory it does work

**§6v1, 2026-09-21.** Pre-registration `docs/PREREG_ideal_chromosome_sim_2026-09-21.md` (md5
`ff226f41`), committed `2a809b76` before any arm was scored. Tool `bench/ideal_chromosome_sim.py`.
⚠HUMAN substrate (A119b/CHM13 chr16) — do not pool with gorilla.

## Question

*"Simulate an ideal scenario using 1 chromosome with some multi-copy gene families but also some others
that are single copy, to prove that in theory it should work. We just need a way to remove unspecific
phenomena like readthroughs and ensure the reads are complete."* — establish the CEILING, and attribute
the real-data gap.

## Setup

chr16, **4,411 annotated transcripts over 1,443 genes**, 10 reads per transcript, ends jittered +-0-30 bp
(mandatory — identical reads collapse under dedup, §6n0), err 0.001. Aligned with the shipped
`-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes` and assembled through the shipped
`--assemble-only --assembly-polish full` path, unchanged.

Strata fixed from the SHIPPED chr16 homology graph before any arm ran: **MULTI-COPY** = the gene has >= 1
homology edge (297 genes); **SINGLE-COPY** = none (1,110). Universe = the 1,407 simulated genes that carry
an annotation record — an INPUT, so nothing is conditioned on any arm's output.

The readthrough rate is measured, not invented: **33,058 of 439,985 (7.51%)** real primary MAPQ-60 chr16
reads hit the exons of >= 2 distinct annotated genes at >= 25 bp. A_rt injects molecules splicing a gene's
transcript to a neighbouring same-strand gene's at that rate (3,014 of 43,700 reads = 6.9% realised).

## Result

| arm | stratum | raw endpoint | unscoreable | **genuine over-merges** | corrected |
|---|---|---|---|---|---|
| **A_ideal** | ALL | 0.8067 | 15.8% | **12** | **0.9578** |
| | MULTI-COPY | 0.6364 | 28.6% | 12 | 0.8915 |
| | SINGLE-COPY | 0.8523 | 12.3% | **0** | **0.9723** |
| A_trunc | ALL | 0.8188 | — | — | — |
| A_rt | ALL | 0.6503 | 12.6% | **188** | 0.7439 |
| | MULTI-COPY | 0.5556 | 22.9% | 46 | 0.7205 |
| | SINGLE-COPY | 0.6757 | 9.8% | 142 | 0.7493 |
| real data | ALL | 0.4961 | 13.9% | 139 | 0.5759 |

"Unscoreable" = the gene's node is shared ONLY with an annotated gene that overlaps it. Two annotated
genes sharing exonic sequence cannot be separated into one-to-one nodes by ANY method, so these are a
truth ceiling, not an error. "Corrected" excludes them; it is **post-hoc** and is reported beside the
pre-registered raw number, never instead of it.

## Verdict against the pre-registered bars

- **Bar 1 (A_ideal >= 0.95 overall and >= 0.90 multi-copy) FAILS as written**: 0.8067 / 0.6364.
- **Bar 2 (A_ideal < 0.80 ⟹ algorithmic defect) NOT triggered**: 0.8067.
- **Bar 3 (single-copy >= 0.95) FAILS as written**: 0.8523.
- **Bar 4 (attribution) PASSES decisively** — see below.

⚠**The bars were mis-set, and that is a prereg design flaw worth recording.** The endpoint has a
built-in ceiling of roughly **0.84** on this substrate, because 15.8% of the universe consists of
annotated genes that overlap another annotated gene. A bar of 0.95 on a metric that cannot exceed ~0.84
was unreachable by construction. This is the same ~0.8 ground-truth ceiling already documented in
§6kl/§6km, arrived at independently.

## What the simulation actually shows

⭐**With full-length reads and no readthrough the algorithm makes almost no errors**: **12 genuine
over-merges across 1,407 genes**, and **ZERO in the single-copy stratum**. Corrected for the truth
ceiling that is 0.9578 overall, 0.9723 single-copy, 0.8915 multi-copy. **So yes — in theory it works,
and multi-copy families are recovered as separate nodes.**

⭐**Readthrough alone reproduces the real pathology, and then some.** Injecting it at the measured 7.51%
takes genuine over-merges from **12 to 188** — a 15.7x increase — against **139** in the real data. The
injected rate slightly OVERSHOOTS the real defect, so readthrough is not merely sufficient to explain the
over-merge found in §6u7, it is more than sufficient. Nothing else needs to be invoked.

⭐**5' truncation costs nothing on this endpoint** (A_trunc 0.8188 vs A_ideal 0.8067, and the identical
multi-copy rate 0.6364). ⚠This does NOT contradict §6n0, where truncation cost 6 of 26 NPIP copies:
that endpoint was transcript COMPLETENESS, this one is where the node SITS. Truncation shortens models
without moving them.

## Consequence

The gap between the ideal ceiling (0.9578 corrected) and real data (0.5759 corrected) is attributed:
**readthrough carries it**. This is the third independent line pointing at the same conclusion as §6u7
(the defect is over-merge, not fragmentation) and §6v0 (full-length-ness corroborates bridges rather
than discriminating them) — and it is the first to show the algorithm is clean when the phenomenon is
removed. ⟹ **the open lever remains a readthrough-aware node SPLIT**, and it now has a measured ceiling
to aim at.
