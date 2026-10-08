# Pre-registration — a depth-adaptive locus boundary (fixed QUANTILE instead of fixed RANK)

**Written 2026-09-22, §6w6, before any arm is run.** User: *"should there be a depth based calculation we
use to improve metrics?"* → *"ok lets do that."*

## The defect, established twice independently

The transcript boundary is the **k-th most extreme read start/end**, `k = min_terminal_support`, default
**2, the same at both ends** (`denovo_assemble.rs:184`, `:325-345`). With *n* reads the 2nd-most-extreme
sits ~`1.5/n` into the tail — **75% of the way in at n=2, but 5% at n=30**. A fixed rank is a moving
quantile, and the drift is measurable in two places that share no metric and no population:

| | low depth | high depth |
|---|---|---|
| **r975** (§6w3, all single-gene loci, boundary error) human median `d5` | +19 bp (n=2) | **−40 bp** (n≥30) |
| same, gorilla NC_073244.2 | +86 bp | +12 bp |
| **r978** (§6w5, Soto families genome-wide, WIDTH) under-capture `ratio<0.5` | **28.2%** (n=2) | 11.2% (n=10-29) |
| same, over-extension `ratio>2` | 12.7% | **25.5%** (n=30-99) |

Deep loci overshoot, shallow loci fall short — on both metrics, both directions, two species.

## The rule

    k = clamp(round(q * n), 1, n)          # q fixed; today's behaviour is k = 2, i.e. q = 2/n

At n=2 this gives k=1 (the outermost read → the locus GROWS, which is the direction r978 says low-depth
loci need); at n=100 and q=0.05 it gives k=5 (the boundary pulls IN, which is the direction high-depth
loci need). One parameter, and it moves both failure modes the right way.

⚠ Implementation note: `groups` currently truncates to the k most extreme values (`:320-331`), so the
full order statistic is not retained. `allpos` already collects every start/end per group but only when
`snap.is_some()` (`:302`, `:313-320`). The arm enables that collection unconditionally under the new mode
and takes the exact k-th order statistic — no approximation.

## Why the score is GRAPH ADMISSION, not boundary error

Per **r971/§6w2**: `cov_longer`'s denominator is the LONGER locus's whole span, so an over-long locus
evicts its correctly-assembled partners from the family graph. On chr16, **1,892 loci have real homology
and only 864 become graph nodes; 99.8% of the 1,026 dropped fail `cov_longer` alone, and 39.6% are
evicted by a partner holding ≥2 whole genes.** r978's over-extension is therefore feeding r971's
eviction, and boundary error in bp is the wrong place to look for the payoff.

⚠ Expected size, stated up front so a small number is not spun as a win: the depth drift is
**~60 bp (human) / ~74 bp (gorilla)** against a **~150-250 bp irreducible scatter** (r973/r974). This
tightens depth-dependence; it cannot touch the per-library bias or the scatter.

## Substrates — parameter on human, verdict on gorilla

- **Development (q chosen here, and ONLY here): human chr16**, A119b (`/mnt/linuxdisk/tmp/regress/dn16*`
  recipe, rebuilt under each arm).
- **Held out: gorilla NC_073244.2** (`GGO_mm.bam`), different species AND library. Reported with **no
  re-tuning** ([[feedback_hold_a_substrate_back]]). ⚠Never pooled with human.

## Arms — the control is what makes the result attributable

| arm | k |
|---|---|
| **A0** shipped | k = 2 fixed |
| **A1** adaptive | k = clamp(round(q·n), 1, n), q swept on human over {0.02, 0.05, 0.10, 0.15, 0.20} |
| **C1** control | k = 3 fixed |
| **C2** control | k = 1 fixed |

⚠⚠ **C1/C2 are not optional.** If a fixed k=3 or k=1 reproduces A1's gain, then the gain is a constant
offset and **adaptivity bought nothing** — which is exactly the claim under test. A1 must beat the better
of C1/C2 to be called a depth effect at all.

## Metrics

**Primary:** graph admission — loci with ≥1 non-self PAF record that clear
`identity ≥ 0.70 ∧ cov_longer ≥ 0.30 ∧ alen ≥ 300bp` and reach the family graph
(`bench/node_graph_admission.py`, whose reimplementation matched the shipped graph 866 vs 864).

**Guards — all must hold, because a rule that merely SHRINKS every locus inflates admission for free:**
1. median exonic Jaccard vs Soto truth must not fall by >0.01 (`bench/soto_family_locus_fidelity.py`);
2. width-ratio in-band fraction `[0.8, 1.25]` must not fall;
3. total emitted transcripts must not fall by >5%;
4. matching intron chains vs RefSeq (gffcompare) must not fall.

**Secondary (reported, not decisive):** median |d5|, |d3|, and the depth-stratified under-capture /
over-extension table from r978 — which must show the two rates CONVERGING if the mechanism is real.

## The bar — committed now

Judged on **held-out gorilla**, at the single q chosen on human:

| outcome | verdict |
|---|---|
| admitted nodes **+3% or more** over A0, **every guard holds**, and A1 beats the better of C1/C2 | ⭐ **ADOPT-WORTHY** |
| +1-3%, guards hold, beats controls | ⚠ **PARTIAL** — real but small |
| <+1%, or any guard fails, or a fixed-k control matches it | ⛔ **NO** |

⚠ **A depth-adjusted or depth-weighted headline METRIC is explicitly NOT part of this** and is rejected in
advance: read count is a property the construction assigns, so conditioning a score on it is the
"denominator conditioned on the prediction" trap ([[feedback_metric_traps]]), and changing k moves genes
BETWEEN depth strata, so within-band before/after comparisons would not compare the same populations.
Depth stratification is used for DIAGNOSIS only.

I will not change the arms, the substrates, the guards or the bar after seeing any number. q is fitted on
human chr16 only.

---

# OUTCOME (2026-09-22) — ⛔ **NO. Refuted on development; the held-out substrate was NOT spent.**

Development, human chr16, A119b, all arms end-to-end (assemble → locus bodies → all-vs-all asm20 →
graph admission via `bench/node_graph_admission.py`):

| arm | transcripts | loci | with homology | **ADMITTED** |
|---|---|---|---|---|
| **A0 shipped (k=2)** | 10,093 | 2,544 | 1,904 | **867** |
| A1 q=0.02 | 10,093 | 2,544 | 1,894 | 865 |
| A1 q=0.05 | 10,091 | 2,539 | 1,887 | 861 |
| A1 q=0.10 | 10,090 | 2,533 | 1,872 | 854 |
| A1 q=0.15 | 10,088 | 2,532 | 1,866 | 852 |
| A1 q=0.20 | 10,088 | 2,527 | 1,859 | **840** |
| C k=1 (control) | 10,093 | 2,544 | 1,904 | 867 |
| C k=3 (control) | 10,131 | 2,580 | 1,874 | 821 |

**Every quantile arm is WORSE than shipped, monotonically in q.** Not one reaches baseline, let alone the
pre-registered +3%. Per the bar this is ⛔ **NO**, and **gorilla was deliberately not run** — spending a
held-out substrate on a rule that fails development is exactly what holding it back is for.

## Why it fails — the prediction was directionally wrong about the cost

The rule does what it was designed to do: raising q raises k for deep loci, which pulls their boundaries
IN. §6w2/r971 predicted that should ADMIT more nodes, because a shorter giant stops inflating
`cov_longer`'s denominator for its partners. What actually happens is visible in the **`with homology`**
column, which falls in lockstep (1,904 → 1,859): **shrinking a locus also costs it its own alignments.**
A shorter body has less sequence to align, so it drops below `alen ≥ 300bp` / the coverage floor and
leaves the graph itself. The gain to partners is real but smaller than the loss to the shrunken locus,
so the net is negative at every q tested.

⭐**The shipped k=2 sits at a local optimum on this metric**, approached from three directions:
k=1 is a **no-op**, k=3 costs 46 admitted nodes, and the quantile costs 2-27.

## Incidental, and worth knowing

⚠**`RUSTLE_TERMINAL_K=1` is byte-identical to the default on FULL chr16 with the polish** (verified by an
independent re-run: `verify_k1.gtf == A0.gtf`), **but NOT on `chr16:1-30,000,000`**, where k=1 vs k=2
differ both with and without polish (11,381 vs 11,350 transcripts unpolished). So *k=1's no-op status is
a property of the region+polish combination, not of k=1* — a reminder that an assembly arm validated on a
sub-region can behave differently on the whole contig, and that the polish absorbs some boundary changes
entirely.

## What survives

r975/r978's **diagnosis** is untouched — a fixed rank IS a moving quantile, and the depth drift is real
and measured on two metrics and two species. What is refuted is the **remedy**: correcting the drift at
the boundary costs more alignment than it recovers admission. Any future attempt must shrink the giant
**without shortening the sequence that earns its own edges** — i.e. a node SPLIT, not a boundary pull-in.
