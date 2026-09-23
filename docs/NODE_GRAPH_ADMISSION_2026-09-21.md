# What is missing for the node definition: it is not the boundary, it is graph admission

**§6w2, 2026-09-21.** User goal: *"analyse what is missing for improving the node definition."*
Substrate: **human** A119b IsoSeq vs CHM13, chr16 de novo (`/mnt/linuxdisk/tmp/regress/dn16.*`).
⚠Never pooled with gorilla.

## Summary

Every node-construction lever tried since §6u7 has attacked **where a locus's boundaries fall**. Measured
end to end on chr16, that is not where the loss is. Of the 2,550 assembled loci, **1,892 have real
homology to another locus, and only 864 become nodes in the shipped family graph.** The other **1,026
(54.2% of everything with homology) are evicted by ONE gate**, and it is not the one the recent work has
been aimed at.

| stage | loci |
|---|---|
| assembled loci (chr16) | 2,550 |
| ...with any non-self PAF alignment | **1,892** |
| ...clearing `identity ≥ 0.70 ∧ cov_longer ≥ 0.30 ∧ alen ≥ 300bp` | 866 *(reimplementation)* |
| ...actually in the shipped graph `dn16.graph.tsv` | **864** |

The reimplementation lands within **2 loci of the shipped graph (866 vs 864)**, so the decomposition
below is faithful to the shipped rule, not to an approximation of it.

## Finding 1 — the binding gate is `cov_longer`, alone

Reason each of the 1,026 dropped loci fails, taken from that locus's own best-coverage pair:

| reason | loci | share |
|---|---|---|
| **`cov_longer < 0.30`** | **1,024** | **99.8%** |
| `cov_longer < 0.30` + `identity < 0.70` | 2 | 0.2% |
| identity alone / length alone | **0** | 0% |

**Identity never rejects anything and neither does the 300 bp floor.** The entire node→graph loss is one
coverage threshold. This sharpens §6u4's "identity/cov gate 22.7% of FNs" into a specific, single knob,
and it confirms §6v2's observation ("17 of 24 missing genes DO have an assembled locus") at full scale
rather than on a 24-gene sample.

## Finding 2 — over-merge's dominant cost is collateral, not self-inflicted

`cov_longer = aligned_union_on_the_longer_locus / length_of_the_longer_locus`
(`annotation_families.rs:407` for the immediate path; the deferred path used by the shipped
`--min-exonic-bp 1` config unions merged intervals and divides by the same denominator).

The denominator is **the longer locus's whole span**. So when a readthrough fuses two genes into one
giant locus, that locus does not only mis-name itself — **it inflates the denominator for every
correctly-assembled partner it aligns to, and evicts those partners from the graph.**

Measured on the 1,026 dropped loci:

| | |
|---|---|
| best partner is LONGER than the locus itself | **681 / 1,026 = 66.4%** |
| median partner/self length ratio | **3.92×** |
| median length of the evicted locus itself | 7,824 bp (ordinary) |
| of those 681, culprit partner contains **≥ 2 whole annotated genes** | **406 / 681 = 59.6%** |

⟹ **406 of 1,026 graph-admission losses (39.6%) are caused by an over-merged partner**, not by anything
wrong with the evicted locus.

Put beside §6u7's first-order count — **172 genuinely-fused gene pairs** — the collateral damage is
**~2.4× larger than the over-merge itself**. The load is diffuse, not a few monsters: **380 distinct
culprit loci, top-10 responsible for only 15.6%** of evictions.

## Consequence — the node-split work was aimed at the smaller half of its own payoff

Every node-split attempt (registers 845/846/937-940/948-949/967/968) was scored on whether the FUSED
locus gets split and correctly named. That is the first-order effect. This measurement says the
second-order effect is bigger: **a successful split shrinks the denominator and re-admits the innocent
partners.** A split trigger that looked worthless on per-copy correctness could still be worth several
hundred re-admitted nodes — and no arm so far has been scored that way.

⚠ This does **not** license changing the norm. `min(la,lb)` is already refuted (**r913**: long genes
become hubs, 864 pairs ≥ 1.0, largest component 171 vs 23, precision .164 vs .639), and under the faithful
accumulation 74.9% of these loci would clear a `min()` floor — which is exactly the hub blow-up r913
measured, not a fix. `RUSTLE_COLLAPSE_EXONIC` is likewise refuted as a default (**§6ab**, hub fusion).
**The denominator is not the bug; the giant locus in the denominator is.**

## Finding 3 — a positional partner-discontinuity probe: real signal, narrow reach

Every refuted split trigger used **read-level** evidence inside the locus, and all failed the same way: a
readthrough molecule is a genuine, abundant, full-length FLNC transcript, so the bridge *is* the
population. Graph-structural splits were tried on graph **topology** (r300 bridges 0.3%, r301 connectivity
*inverts*, r827 λ≥2) but never on **where along a locus's own axis its partners align**.

    discontinuity(L) = min over cut c of  Jaccard( partners left of c , partners right of c )

with both sides required to carry ≥5 *exclusive* partners and the cut confined to [0.15, 0.85].
⚠The naive form is degenerate — without those floors the minimum always lands at an extreme cut with two
partners on one side. Length-matched, because a longer locus has more partners and more chances to find a
low-J cut:

| length band | 1 gene inside | ≥2 genes inside |
|---|---|---|
| <20 kb | 0.366 (n=13) | 0.516 (n=3) ⚠n |
| 20–50 kb | 0.356 (n=15) | **0.127** (n=14) |
| 50–100 kb | 0.400 (n=8) | **0.105** (n=11) |
| ≥100 kb | n<3 | 0.038 (n=14) |

Within a band the separation is ~3×, and it is not a length artifact. **But the reach is narrow and it
fails on the canonical cases**: only **145 of 2,550 loci (5.7%)** carry enough partners to be scoreable;
the `NPIPB4`+`RRN3P1` 65 kb super-locus scores **0.556** (no signal), `CDR2`'s 91.6 kb locus is **not
scoreable at all**, and only the `PKD1P6` 190 kb locus separates (0.007). **Population signal, not a
per-locus decision rule** — the same shape as §6u3's neighbourhood Jaccard.

## What this rules out, and what it leaves

Dead, do not re-propose: the `min()` norm (r913) · `COLLAPSE_EXONIC` as default (§6ab) · every read-level
split trigger (967/968/948/949) · depth valley (the rule may not be quoted — `c=0.20` is the argmax
against the truth denominator; its null was matched on edge count only) · graph-topology splits
(r300/r301/r827) · scoring-layer multi-labelling (964, §6v9).

Left standing, in priority order implied by the numbers:
1. **Re-score the existing split triggers on graph admission, not per-copy correctness.** r968's
   turnover trigger moved per-copy +0.78pp and was judged 4× short of its bar; its effect on the 406
   collateral evictions was never measured. This is a re-scoring of work already done, not new machinery.
2. **ORF architecture of the assembled transcript** — still zero register hits as a *split* criterion
   (r351/r746/r762 use ORFs only at the family/edge layer). ⚠Weak exactly where it matters: NPIP-type
   families are largely pseudogenes with no ORF to find, and 35.9% of Soto families are all-pseudogene.
3. **Internal 3′-end (TES) pile-up** — r543 measures TES at ~5× the signal of TSS and r750 guarantees
   every FLNC read end is a genuine transcript 3′ end; `DenovoTranscript` already carries `tes` from
   `sharp_tes` (`denovo_assemble.rs:161`). ⚠Inherits the route-3/T20 tautology risk: those same read ends
   set the locus span in the first place.

## Reproduce

```sh
cd /mnt/linuxdisk/tmp/regress          # dn16.paf, dn16.graph.tsv, chr16.genes.gff

# Findings 1 and 2 — gate decomposition and the collateral-eviction measurement
python3 bench/node_graph_admission.py --paf dn16.paf --graph dn16.graph.tsv --gff chr16.genes.gff

# Finding 3 — positional partner discontinuity, length-matched
python3 bench/locus_probes.py partner-discontinuity --paf dn16.paf --gff chr16.genes.gff \
  --locus chr16:22357769-22422849 --locus chr16:21779645-21871278 --locus chr16:14975608-15166427
```

`cov_longer` must be computed the **deferred** way (union of merged intervals on the longer locus divided
by that locus's length), because the shipped `--min-exonic-bp 1` config takes the deferred path and
accumulates every record for a pair. A best-single-record approximation **overstates** the effect
(78.0% vs 66.4% longer-partner; ratio 9.79× vs 3.92×) — the agreement check (866 vs the shipped 864) is
what licenses the decomposition and should be re-checked on any other chromosome before quoting it.

> **Generator (2026-09-22 consolidation):** `python3 bench/locus_probes.py partner-discontinuity ...` — the original `bench/partner_discontinuity.py` was folded in verbatim and verified identical on its documented inputs (§6z3); the register rows above cite this file.
