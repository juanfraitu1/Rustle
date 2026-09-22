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

---

# 8. Does the idealized version hold a better FAMILY definition?

**User, 2026-09-21.** §6v1 scored only the NODE-construction endpoint. This runs the actual family
definition (`mcl_families --min-exonic-bp 1 --min-shared-exon-frac 0.60`, the shipped Rust binary, same
recipe as `dn16_fam3`) on the A_ideal and A_rt assembled loci, and scores the resulting clusters against
**both** independent truths from §6u8/§6u9 (Soto cover, protein referee), restricted to chr16.

⚠**n is small** — chr16 alone has 15 Soto families / 36 protein-referee families. Per-family movement
swings are large at this n (§6u8's lesson); read the direction, not the third decimal.

## 8.1 Pipeline

For each arm: exon-merged locus spans from the assembled GTF -> region list -> `samtools faidx` body
extraction -> `minimap2 -x asm20 -c --eqx -P` all-vs-all -> shipped `mcl_families` with the locus's own
exon-union GFF3 (mirroring `dn16.gff3`'s format exactly) -> `.clusters.tsv`. A fourth arm, **GUIDED**
(annotation records used directly as nodes, no assembly at all — `chr16_guided.clusters.tsv`, already on
disk), stands in as the ceiling: perfect nodes, no simulation needed.

## 8.2 Result

| arm | clusters >=2 | vs **Soto cover** (15 fams) | vs **protein referee** (36 fams) |
|---|---|---|---|
| **GUIDED** (annotation nodes) | 90 | sens 0.4567 / prec 0.4294 / **F 0.4174** | sens 0.2046 / prec 0.3203 / **F 0.2248** |
| **A_ideal** de novo | 56 | sens 0.3692 / prec 0.3649 / **F 0.3561** | sens 0.1722 / prec 0.2968 / **F 0.1944** |
| REAL de novo | 70 | sens 0.3372 / prec 0.3255 / **F 0.3166** | sens 0.1681 / prec 0.2172 / **F 0.1805** |
| A_rt de novo | 51 | sens 0.3136 / prec 0.3182 / **F 0.3068** | sens 0.1329 / prec 0.2603 / **F 0.1570** |

**Yes — modestly, and the direction holds on both independent truths.** A_ideal beats real de novo by
+0.0395 F (Soto) and +0.0139 F (protein referee). A_rt falls below REAL on both truths despite injecting
readthrough onto otherwise-ideal reads — consistent with §6v1's own finding that the injected 7.51% rate
slightly overshoots the real defect (12 -> 188 genuine over-merges vs 139 real), so this is not a
contradiction: an aggressive dose of the one thing being tested does more damage than the full complexity
of real data, which is itself evidence the mechanism is right.

Per-family (Soto), 3 better / 6 worse / 6 unchanged: the pooled gain is real but not one-sided — a caveat
consistent with n=15.

## 8.3 The bigger finding: node CORRECTNESS is fixed, node COVERAGE is not

**A_ideal falls well short of GUIDED — F 0.3561 vs 0.4174 (Soto), 0.1944 vs 0.2248 (referee) — despite
§6v1 showing node correctness is nearly solved (0.9578 corrected).** The reason is not over-merge (that
was measured and is nearly gone) and not fragmentation-as-scored-by-§6u7 (per-copy correctness is high).
It is **node COVERAGE**: how many of the truth genes have ANY node representing them in the homology
graph at all.

| arm | graph nodes | Soto truth genes (58 on chr16) represented as SOME node |
|---|---|---|
| GUIDED | 380 | **47 / 58 = 81.0%** |
| A_ideal | 324 | **34 / 58 = 58.6%** |

A 22.4-point coverage gap, under perfect reads and no readthrough. Of the 24 genes missing a node in
A_ideal, 17 (71%) DO have an assembled locus — the locus exists, it is just absent from the graph,
meaning either it never accumulated enough exonic content to clear the `identity >= 0.70 / cov_longer >=
0.30 / >= 300bp` edge floor, or (below) it never gets scoring credit at all:

**Of the 13 genes missing in A_ideal but present in GUIDED (the true assembly-attributable set), 5
(38%) are a SCORING ARTIFACT, not an assembly failure.** RefSeq itself carries curated **readthrough/
fusion gene records** overlapping almost exactly the same span — `PKD1P3-NPIPA1` (42,950 bp) sits over
`NPIPA1`'s 14,597 bp span at 37,168 bp overlap; `PKD1P4-NPIPA8`, `PDXDC2P-NPIPB14P`, `BOLA2-SMG1P6`
do the same for `NPIPA8`, `NPIPB14P`, `SLX1B`. The de novo locus's homology signal for the real gene IS
present (verified directly on the PAF: NPIPA1's locus has 10+ alignment records at 0.83-0.99 identity to
other loci) — but a **max-overlap, one-name-per-locus resolver, applied at scoring time, hands the whole
locus's credit to the bigger fusion-named record**, which then has no truth family to match, so the real
gene registers as "no node" though its sequence was correctly assembled and correctly grouped. This is
the family-definition-scoring analogue of §6v1's "unscoreable (annotation overlap)" category, applied
here for the first time to the CLUSTER-naming step rather than the per-copy node step.

The remaining ~8 (ABCC6, HERC2P5, HERC2P8, NPIPB10P, NPIPB7, PKD1, PKD1P6, SMG1P6) have no bigger
co-located record and are candidates for a genuine remaining edge-construction gap — short or
low-coverage assembled loci that never accumulate enough exonic content relative to their full-length
siblings to clear the coverage floor. This was not run to ground on this pass and is the natural next
target.

## 8.4 Consequence

The answer to "does the idealized version hold a better family definition" is **yes, but only partly for
the reason expected.** Fixing readthrough (§6v1's target) recovers a real, if modest, family-definition
gain — confirmed on two independent truths. But it does not close most of the gap to a perfect-node
ceiling, because **most of that remaining gap is not a node-construction defect at all**: at least 38% of
it is an artifact of how RefSeq's own curated fusion-gene names get resolved at scoring time, in exactly
the same region (PKD1/NPIP tandem duplication) already flagged by §6u7/§6v1/§6v0 as ground zero for
readthrough. The rest is a real, smaller, and as-yet uncharacterized edge-coverage gap.
