# Pre-registration — adapting `RUSTLE_ER_COVERAGE_LONGER_FLOOR` to the CURRENT definition

**Written 2026-09-22, §6y7, before any floor is swept.** User: *"can we adapt the flag to the current
definition?"*

## What the adaptation actually is

§6y6/r1032: the E_r flag is unreachable from the shipped catalog's driver. The current definition lives in
`mcl_families`/`annotation_families.rs`. Its faithful analogue is a **CONJUNCT**, not an escape:

| layer | longer side | shorter side |
|---|---|---|
| E_r (RNA, the flag's home) | `RUSTLE_ER_COVERAGE_LONGER_FLOOR` (**added** clause) | `min_coverage` 0.50 (shipped) |
| **DNA, current definition** | `cov_longer` ≥ **0.30** (shipped) | **`min_shared_exon_frac` ≥ 0.60** (shipped) |

⚠⚠**The current definition is ALREADY two-sided** — it charges the longer side on alignment span and the
shorter side on **exon-to-exon overlap as a fraction of the smaller gene's exonic length**. So BLAST's
`qcovs AND scovs` shape is structurally present, with a different quantity on each side.

⚠**Semantics differ from `--min-cov-shorter` (§6x4, shipped this session), which is an OR-escape and
LOOSENS.** This is an AND-conjunct and TIGHTENS. They are opposite directions and must not be conflated.

## The question this actually tests

**Does an ALIGNMENT-COVERAGE floor on the shorter side remove anything the EXON-SHARING floor does not?**
If the two clauses are near-redundant the adaptation is a no-op and the definition already has it.

## The arm

Edge set = the shipped `--dump-graph` output (post every shipped conjunct). Each surviving edge gets
`cov_shorter = merged alignment on the shorter locus / that locus's exonic length`, computed the deferred
way (validated in §6x4 against `node_graph_admission.py` at 866/864/679 exactly). A conjunctive floor
`S ∈ {0.30, 0.50, 0.70}` then REMOVES edges below it; families are rebuilt with `mcl_port` on the
identical graph minus those edges.

⚠**Comparator is `mcl_port` on the unfiltered graph, never the shipped Rust F** (r917).

## Metrics · substrates · bar

Pairwise precision / recall / F against the **protein-family referee** (neutral) and **Soto** (reported,
SD-scoped). Plus **edges removed** and **redundancy**: what fraction of removed edges the exon conjunct
would have removed at a stricter `min_shared_exon_frac`.

Development **chr16**; held out **chr2 / chr8 / chr10**, no re-tuning.

| outcome | verdict |
|---|---|
| a floor raises held-out F by **≥0.02** on the referee without losing >0.02 on Soto | ⭐ **ADAPT IT** — add the conjunct to `mcl_families` |
| F within ±0.02 but precision up ≥0.02 at recall cost ≤0.02 | ⚠ **PARTIAL** — a precision knob, document and leave off |
| **F down, or <1% of edges removed at every floor** | ⛔ **NO** — redundant with the exon conjunct; the definition already has its two-sided rule |

**Predicted, before looking — ⛔ REDUNDANT.** `min_shared_exon_frac 0.60` already demands that the best
record's exon-to-exon overlap cover 60% of the smaller gene's exons; a pair passing that can hardly have a
low alignment coverage on the shorter side. I expect very few edges removed and no F movement. If instead a
floor removes a meaningful slice AND improves held-out F, the two clauses are measuring different things
and the adaptation is real.

I will not change the arm, the floors, the truth or the bar after seeing any number.

---

# OUTCOME (2026-09-22) — ⛔ **NO. It removes real edges and costs recall without buying precision.**

Held-out pairwise F against the protein referee, by shorter-side floor `S` (span/span):

| S | chr2 | chr8 | chr10 | chr16 (dev) | edges removed |
|---|---|---|---|---|---|
| none | **0.417** | **0.934** | **0.206** | **0.881** | — |
| 0.30 | 0.414 | 0.934 | 0.206 | 0.879 | 0.3–1.6% |
| 0.50 | 0.414 | 0.934 | 0.206 | 0.879 | 1.0–3.7% |
| 0.70 | 0.409 | 0.934 | 0.206 | 0.876 | 1.8–8.7% |
| 0.90 | 0.382 | **0.760** | 0.185 | 0.882 | 9.9–22.0% |

⛔**Every held-out cell is flat or down**, and precision does not rise to compensate — chr2 .975→.972,
chr8 .968→.954, chr10 .857→.842 at S=0.90, all *falling*. chr8 at 0.90 is the clearest: removing **11.8%**
of edges costs recall **0.903 → 0.631**. Only chr16 (development) shows the hoped-for shape
(precision .968→.981, F .881→.882) and it does not transfer.

## ⭐ Why — the same mechanism that refuted the SYMMETRIC variant at the RNA layer

The removed edges carry real pairs. `denovo_pipeline.rs:4900-4926` already recorded the cause for E_r:
**only 134/171 NPIP true pairs can reach 0.50 on the longer axis at all** (`NPIPB8-NPIPB2` caps at 0.215)
*"because the duplicated unit is size-invariant while annotated spans are not."* The identical argument
applies to the shorter side at the DNA layer: a genuine paralogue pair routinely has low span coverage on
one side, so a span-coverage floor deletes true edges wherever the two spans disagree in length.

⭐**And the definition already has its shorter-side conjunct — on the right quantity.**
`min_shared_exon_frac 0.60` charges exon-to-exon overlap as a fraction of the **smaller gene's exonic
length**, and it is strongly binding: it rejects **965** pairs on guided chr16, **3,874** on de novo chr16,
**834** on chr2. Exon sharing is comparatively size-invariant; span coverage is not. **The adaptation the
user asked for is present, and it is on the better quantity.**

## ⚠⚠ A units bug of my own, and it produced a FALSE no-op I nearly reported

My first pass computed `cov_shorter` as *genomic alignment span ÷ exonic length* — mismatched units. The
quantity ran to a **median of 6.87 and a max of 2,416** on chr16, so every edge trivially cleared any floor
below 1.0 and the sweep reported **0 edges removed at S=0.30/0.50 and 1 at 0.70**, i.e. "redundant, a
no-op". With matched span/span units the same floors remove **1.6% / 3.7% / 7.0%**. ⭐**The tell was in the
printed distribution — a "coverage fraction" with a median of 6.87 is not a fraction** — and the
`annotation_families.rs` comments warn about exactly this ("the span numerator and the exonic denominator
are in DIFFERENT UNITS"). **Print the distribution of any ratio before thresholding it.**

## Prediction scorecard

Predicted ⛔ REDUNDANT — *"a pair passing `min_shared_exon_frac 0.60` can hardly have low alignment coverage
on the shorter side; I expect very few edges removed and no F movement."* **The verdict is right and the
reasoning is wrong.** The floors remove 2–22% of edges, not "very few"; they are not redundant, they are
*harmful*. My prediction happened to match the first (buggy) run, which is exactly how a units error
survives review — it confirmed what I already expected.

> **Generator (2026-09-22 consolidation):** `python3 bench/edge_probes.py cov-shorter ...` — the original `bench/cov_shorter_conjunct.py` was folded in verbatim (§6z2); the register rows above cite this file.
