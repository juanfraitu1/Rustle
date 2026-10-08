# Pre-registration — an empirical rule for the missing 5′ of a constructed locus

**Written 2026-09-21, §6w3, before any deficit is measured.** User goal: *"find an empirical rule to
account for missing 5′ in representative constructed locus from reads."*

## What is already established, and what it licenses

- **r806** — the de novo exon-sum does **not** systematically under-represent a copy's width (median
  relative deficit −0.1% on both species). ⭐But: *"~100% of whatever deficit remains is at the 5′ end;
  the 3′ end is exact (< 30 bp) on every tool tested, including flair/StringTie/isoseq."* **That
  asymmetry is the thing this pre-registration is about**, and the median-zero result means the rule must
  be judged on the *tail*, not the median, or it will look like a no-op.
- **r912** — `RUSTLE_TSS_SNAP` (5′ snap to a read-start peak) is *"opt-in for absence of benefit, not
  harm"*: the claimed "makes size agreement worse" was itself refuted (paired sign test p = 0.69). So a
  5′ extension is neither proven to help nor proven to hurt — it is unresolved, not refuted.
- **r396** — trimming a locus to its read-supported core made truncation **worse** (64% → 75%). The
  correction must therefore be an *extension* question, not a trimming one.
- **§6p4** — CDS recovery is 99.8% but TBC1D3's 5′UTR is 75.3% vs 3′UTR 100% ⟹ *5′/TSS-based rules are
  the unsafe ones*. This is a warning about the rule's downstream use, not a reason not to measure it.
- **The boundary rule itself** (`denovo_assemble.rs:184, 342-345`): the transcript extent is the **k-th
  most extreme read start/end**, `k = min_terminal_support`, default **2** — *the same k at both ends*.
  Zero register hits for `min_terminal_support` or for an asymmetric k. Given r806 says the 3′ end is
  already exact at k=2 and the 5′ is not, **an asymmetric k is the untried rule this design is built to
  test.**
- `boundary_gap_left/right` in the GTF is an **outlier flag**, not a deficit measure
  (`copy_assign.rs:4105-4147`: farthest start/end bucket vs the next, gated on read fraction).

## Substrate — development and held-out, never pooled

- **Development: human** chr16 de novo, `/mnt/linuxdisk/tmp/regress/dn16.gtf` vs `chr16.genes.gff`.
- **Held out: gorilla** `NC_073244.2` de novo, `ggo_polished.gtf` vs `winloci_data/GGO_genomic.gff`.
  Different species AND different library ([[feedback_hold_a_substrate_back]]). The rule's parameter is
  chosen on human only and reported on gorilla with **no re-tuning**.

## Population — frozen

De novo loci whose span contains **exactly one** annotated gene at ≥50% of that gene's own span, and
**no second gene** meeting the same bar, **and** whose representative has ≥2 exons.

⚠ The single-gene restriction is not cosmetic: a readthrough-fused locus's "5′ deficit" is **undefined**
(§6v8 showed a fused locus spans unrelated genes), so including them would measure over-merge, not 5′
truncation. Spliced-only because the terminal-exon geometry the rule uses is meaningless on a stub, and
r1063 already showed stub reps are a separate population (median ratio 0.25 vs 0.77).

## Measure — strand-aware, signed, both ends

Using the **annotated gene's** strand:

    + strand:  d5 = locus.start − gene.start      d3 = gene.end − locus.end
    − strand:  d5 = gene.end   − locus.end        d3 = locus.start − gene.start

`d5 > 0` = the locus falls SHORT at the 5′ end (the phenomenon). `d5 < 0` = it overshoots.

Truth is ambiguous at the 5′ end because RefSeq carries several transcripts per gene, so **both are
reported and neither is chosen after the fact**: `d5_repr` against the gene record's own terminus, and
`d5_union` against the most extreme 5′ end over all that gene's transcripts.

## Candidate rules — all functions of de-novo-observable quantities only

| rule | form |
|---|---|
| **R0** null | no correction (the comparator) |
| **R1** constant | extend the 5′ by the development-set median `d5` |
| **R2** depth-scaled | extend by a function of the locus's own `reads` count |
| **R3** asymmetric k | take the 5′ boundary at the **1st** most extreme read start (k=1) while the 3′ stays k=2 |
| **R4** quantile | 5′ boundary at the q-th percentile of read starts, q chosen on development only |

## The bar — committed now

**Gate 0 (does the phenomenon even reproduce here?)** median `d5` must exceed median `d3` by **≥ 100 bp**
on the development set. If it does not, r806's asymmetry does not hold on chr16 and I report that and
stop — no rule is fitted to a phenomenon that is not there.

Given gate 0 passes, a rule is **adopted-worthy** only if, on the **held-out gorilla** substrate:

| outcome | verdict |
|---|---|
| median \|d5\| falls by **≥ 30%** vs R0 **and** median \|d3\| does not worsen **and** the 5′ **overshoot rate stays ≤ 25%** | ⭐ **RULE FOUND** |
| median \|d5\| falls but overshoot > 25%, or \|d3\| worsens | ⚠ **PARTIAL** — trades one error for another |
| median \|d5\| falls < 30%, or it does not transfer from human to gorilla | ⛔ **NO** |

⚠ **The overshoot rate is not optional.** An extension rule that reaches past the true 5′ is not a fix,
it is over-merge fuel — it is the same failure r396 produced in the trimming direction, and §6w2 just
showed that lengthening a locus inflates `cov_longer`'s denominator and evicts its partners from the
graph. **Any 5′ extension must be reported with its effect on locus length**, because length is now known
to be load-bearing.

I will not change the population, the strand convention, the two truth variants, or the bars after
seeing the numbers. Parameters for R1/R2/R4 are fitted on human chr16 ONLY.

---

# OUTCOME (2026-09-21, `bench/five_prime_deficit.py`)

⚠**Numbers are reported per substrate and NEVER pooled** (human chr16 n=4,606 · gorilla NC_073244.2
n=2,929).

## Gate 0 fails as written — and the reason is the finding

| | human chr16 (dev) | gorilla NC_073244.2 (held out) |
|---|---|---|
| median `d5` (signed) | **1 bp** | **67 bp** |
| median `d3` (signed) | 1 bp | 1 bp |
| **gate 0** (median d5 − d3, bar ≥100) | **0 bp ⛔** | 66 bp ⛔ |
| loci short at 5′ (`d5>0`) | **50.0%** | 76.7% |

On human the 5′ end is **not systematically short at all** — it overshoots exactly as often as it falls
short. Gate 0 was operationalised as a *bias* test, and by that test the phenomenon is absent on the
development substrate. **But the phenomenon r806 described is real and large — it is dispersion, not
bias:**

| | human chr16 | gorilla |
|---|---|---|
| **median \|d5\|** | **229 bp** | **136 bp** |
| **median \|d3\|** | **22 bp** | **5 bp** |
| ratio | **10×** | **27×** |
| 5′ within ±30 bp | **26.1%** | **25.1%** |
| 3′ within ±30 bp | 51.4% | 65.4% |

⟹ **r806 reproduces on both substrates, as a variance asymmetry.** The 5′-within-±30bp rate is
**26.1% vs 25.1%** across two species and two libraries — the most stable number in this measurement.

## Consequence: R1, R2, R3 and R4 are all dead on arrival

Every pre-registered candidate is an **extension** rule — it shifts the 5′ boundary outward. Against a
distribution already centred on zero (human: 50.0% short), shifting it outward converts ~50% short into
~75% overshoot. **You cannot correct a zero-bias, high-variance error with an offset**, and §6w2 has just
shown that lengthening a locus is not free: it inflates `cov_longer`'s denominator and evicts partners
from the graph. R3 (k=1) would make overshoot strictly worse.

## The bias that does exist does NOT transfer

Human median `d5` = **+1 bp** (50.0% short) · gorilla = **+67 bp** (76.7% short). A constant fitted on
gorilla is a no-op-to-harmful on human and vice versa. **Per-library, not per-method** — which is exactly
what holding gorilla back was for.

## Over half the 5′ "deficit" is not an error

On human, **55.3%** of loci have `|d5|` within the gene's **own annotated TSS spread** (median spread
**589 bp** across a gene's transcripts). Genes with a single annotated transcript score *worse*
(median |d5| **1,357 bp**, n=1,254) than multi-transcript genes (**154 bp**, n=3,352) — single-transcript
RefSeq records are largely Gnomon/LOC predictions whose own boundaries are unreliable, so part of the
residual is truth-side, not assembly-side.

## The one component that IS a code fix: k is a fixed RANK, so it is a MOVING QUANTILE

The boundary is the **k-th most extreme** read start with **k = 2, fixed, at both ends**
(`denovo_assemble.rs:184, 342-345`). With *n* reads, the 2nd-most-extreme sits at roughly `1.5/n` into
the tail — **75% of the way in at n=2, but only 5% in at n=30.** The measured drift matches exactly, on
both substrates, monotonically:

| reads | human median d5 | human %short | gorilla median d5 | gorilla %short |
|---|---|---|---|---|
| 2 | +19 | 58.7% | +86 | 81.0% |
| 3–4 | +8 | 53.2% | +108 | 82.8% |
| 5–9 | 0 | 48.1% | +65 | 80.3% |
| 10–29 | −14 | 41.1% | +52 | 72.6% |
| ≥30 | **−40** | 34.3% | **+12** | 56.3% |

Deep loci **overshoot** (human −40 bp) and shallow loci fall short — the signature of a fixed rank
sampling further into the tail as depth grows. Swing: **~60 bp human, ~74 bp gorilla.**

## THE RULE

The 5′ boundary error decomposes into three parts, and only one of them is fixable in code:

1. **Irreducible scatter — carry it, do not correct it.** Treat the 5′ terminus as **±150–250 bp** and
   the 3′ as **±5–25 bp**. This transfers across species (26.1% / 25.1% within ±30 bp) and is the first
   *number* behind §6p4's standing "CDS-based rules are safe, 5′/TSS-based rules are not". Any downstream
   rule that reads the 5′ boundary at finer resolution than ~200 bp is reading noise.
2. **Per-library bias — estimate per library, never inherit.** +1 bp human vs +67 bp gorilla.
3. **Depth drift (~60–75 bp) — the only code fix: make k a fixed QUANTILE of read starts, not a fixed
   rank.** Untried (zero register hits for `min_terminal_support`). ⚠Expected gain is ~60 bp against a
   ~200 bp scatter, so it tightens the depth-dependence without moving the bulk — worth doing for
   principle, not for a big metric jump, and it must be re-scored on **graph admission** per §6w2/r971,
   not on boundary error alone.

⚠**Not adopted, not implemented** — this is the measurement and the rule it licenses. Item 3 needs its
own pre-registered arm with an end-to-end score before any default changes.

> **Generator (2026-09-22 consolidation):** `python3 bench/locus_probes.py five-prime ...` — the original `bench/five_prime_deficit.py` was folded in verbatim and verified identical on its documented inputs (§6z3); the register rows above cite this file.
