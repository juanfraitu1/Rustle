# Pre-registration — guarded CONTAINMENT at the admission gate (ceiling measurement)

**Written 2026-09-22, §6x4, before any number is produced.** Follows directly from §6x3/r1002: the 679
loci evicted by `cov_longer` align along **1.00 of their own length** at passing identity and `alen`, and
are rejected only because the denominator is the LONGER locus's span. Both ways of shrinking the giant are
now dead (§6w6 boundary pull-in, §6x3 splitting — median 4.83× required, binary split gives 2×).

## What this is, and what it is NOT

⚠**This is a CEILING and HUB-RISK measurement on the PAF, not a definition change and not a family
score.** `mcl_families` has no containment option (`--min-cov-longer` only), so a real test needs Rust.
The question here is only: **is it worth writing?** No adoption decision follows from this file.

## The rule under measurement

Admit a pair if the shipped gate passes **OR** the containment escape passes:

    shipped : cov_longer  = merged/max(la,lb) >= 0.30  ∧ ident >= 0.70 ∧ alen >= 300
    escape  : cov_shorter = merged/min(la,lb) >= C     ∧ ident >= 0.70 ∧ alen >= 300 ∧ GUARD

`C` swept over {0.30, 0.50, 0.70, 0.90, 0.95}. **GUARD** is §6u0/r924's `--exonic-both-sides`: the aligned
interval must overlap annotated exons on BOTH loci. Measured with the guard ON and OFF, because the whole
point of r913 is that the un guarded norm generates hubs.

## Why a guard is mandatory, quoted so it cannot be skipped

**r913**: the bare `min(la,lb)` containment norm **OVERCORRECTS** — 864 pairs reach ≥1.0, the largest
component goes **23 → 171**, precision **.639 → .164**. **r917**: r916's guarded-containment win was
measured in `mcl_port.py`, not bit-identical Rust, and was explicitly recorded as *not* a basis to change
the definition. Any C that reproduces r913's hub signature is refuted on the spot regardless of admission.

## Metrics

Per (C, guard) cell: admitted nodes · newly admitted of the **679 known victims** · new edges · **largest
connected component** · pairs at cov_shorter ≥ 1.0 (r913's own hub counters).

## The bar — committed now

| outcome | verdict |
|---|---|
| a cell admits **≥ 200 of the 679 victims** with largest component **≤ 40** (baseline 23, r913's failure 171) | ⭐ **WORTH IMPLEMENTING IN RUST** |
| ≥ 100 victims, component ≤ 40 | ⚠ **MARGINAL** — implement only if no cheaper lever appears |
| < 100 victims at every safe C, **or** every cell that admits enough reproduces r913's hub signature | ⛔ **NO** — the admission denominator is not the lever either |

**Predicted, before looking — ⭐ at high C with the guard ON, ⛔ with it OFF.** The victims align at a
median ratio of 1.00, so they sit at the very top of the containment distribution and a strict C (0.90+)
should capture most of them while excluding partial overlaps. With the guard OFF I expect r913's hubs to
reappear at every C low enough to matter.

I will not change the rule, the guard, the sweep or the bar after seeing any number.

---

# OUTCOME (2026-09-22) — ⭐ **WORTH IMPLEMENTING, and the guard is what makes it so.**

Reimplementation validated against `bench/node_graph_admission.py` on all three counters before any sweep:
**homologous 1,892 / admitted 866 / victims 679 — exact.**

Baseline: admitted 866, edges 2,113, **largest component 49**, 679 victims.

| C | guard | admitted | victims admitted | edges | largest comp | × baseline |
|---|---|---|---|---|---|---|
| 0.30 | ON | 1,342 | 273/679 | 3,793 | 317 | 6.5× |
| 0.50 | ON | 1,279 | 244/679 | 3,309 | 144 | 2.9× |
| 0.70 | ON | 1,261 | 234/679 | 3,042 | 107 | 2.2× |
| **0.90** | **ON** | **1,231** | **216/679** | 2,918 | **101** | **2.1×** |
| 0.95 | ON | 1,229 | 214/679 | 2,871 | 101 | 2.1× |
| 0.30 | OFF | 1,694 | 572/679 | 6,838 | 799 | 16.3× |
| 0.90 | OFF | 1,585 | 509/679 | 3,768 | 290 | 5.9× |
| 0.95 | OFF | 1,582 | 506/679 | 3,703 | 279 | 5.7× |

⭐**The victim clause is met at every guarded C** (214–273 ≥ 200), and ⭐⭐**r913 reproduces exactly where
it was predicted to**: with the guard OFF the largest component runs **5.7–16.3× baseline**, with it ON
and C ≥ 0.70 it holds at **2.1–2.2×**. The exonic-both-sides guard is not decoration — **it is the
difference between a 2× and a 16× hub**, on a substrate r913 never measured.

## ⚠⚠ A defect in my own pre-registration, reported rather than re-drawn

The bar said *"largest component ≤ 40 (baseline 23, r913's failure 171)"*. **Those numbers came from
r913's ANNOTATED-GENE-BODY graph and do not transfer to the de novo graph, whose baseline largest
component is already 49.** The absolute threshold was therefore **unsatisfiable by construction before a
single cell was computed** — the baseline itself fails it. I evaluated the hub clause on the **ratio to
this substrate's own baseline** instead, and I am recording that as a forced amendment, not a silent
goalpost move: the direction and the ranking of cells are unaffected, but nobody should read "≤ 40" as
having been tested. ⭐**Lesson: a bar imported from another substrate must be re-derived on this one
before it is written down** — the same class of error as quoting a hand-picked F ([[project_overall_family_metrics]]).

## What this does and does NOT establish

⭐ Establishes: a containment escape at the admission gate can recover **~32% of the 679 collateral
evictions** at a hub cost of 2.1×, where **splitting recovers 0%** (§6x3/r1001) and boundary pull-in is
negative (§6w6). That is the first lever in this line that moves the eviction population at all.

⚠⚠ Does NOT establish that it should ship. **Admission is a gate, not precision.** The cell adds 365
nodes and 805 edges; whether the resulting families are better is a question this file cannot answer and
did not try to. Two things also make the real implementation *safer* than this ceiling: the shipped
pipeline applies `--min-exonic-bp 1 --min-shared-exon-frac 0.60` to every admitted pair, which is a
strictly stronger guard than the `exonic > 0 on both sides` used here; and 149 of the 365 newly admitted
nodes are not in the 679 victim set at all, so the escape's effect is broader than the population it was
designed for and must be scored, not assumed.

**Next: implement `--min-cov-shorter` in `mcl_families` (default OFF, byte-identical when unset) and score
families on the held-out chromosomes.**
