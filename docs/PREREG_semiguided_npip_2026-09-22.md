# Pre-registration — does a SEMI-GUIDED (DNA-region) node set beat de novo and guided on NPIP?

**Written 2026-09-22, §6x0, before any arm is scored.** User goal: *"implement the semi-guided mode for
NPIP, check if it improves precision, sensitivity and bipartite matching."*

## What already exists, and what this is not

⚠ **The mode is already implemented.** `gw_family_catalog --from-genome-sd <SD-pairs BED>` is
GENOME-ONLY mode: "every interval on either side of every SD pair becomes a candidate window
(annotation-free)". It is proposal #2 from §6j5. Nothing is being built here; this is a measurement.

⚠ **Region nodes have been scored once before, and the result is a warning:**
- **r817 / §6ji arm D1** (`--from-genome-sd`, SEDEF): NPIP L1 precision **0.875**, reproducing Dishuck's
  B6–9 exactly — **but 21 of 31 truth records COLLAPSED** (shared their best copy with another record),
  against 0 collapsed for guided and for the de novo default.
- **r817 verdict**: the separation "came from flanking duplicon context in the region nodes", and with
  copy-level SD atoms it **did not hold** — all 22 NPIP genes fell in one family, L1 precision **0.545,
  identical to guided**.

So the live question is not "does it separate NPIPA from NPIPB" (that was context, not method) but
**whether a region node set survives a one-to-one bipartite score at all**, given it is known to collapse
records.

## Substrate — human, and only human

**Human chr16**, CHM13, because the de novo and guided arms already exist there and are the comparators:
- de novo `/mnt/linuxdisk/tmp/regress/dn16.*` (2,550 loci)
- guided `/mnt/linuxdisk/tmp/regress/chr16_guided.*` (annotated gene bodies)
- semi-guided: **2,465 unique chr16 SD intervals** from `HSA_sedef_pairs.bed` (42.8 Mb, median 7,117 bp,
  max 564,039 bp).

⚠**Never pooled with gorilla.** The pre-existing `o1_fromgenome_sd/npip_chroms_sd.bed` is GORILLA
(`NC_073224.2`) and is deliberately NOT used.

## What changes, and what must not

**Only the node interval set changes.** All three arms then run the identical downstream:
same `minimap2 -x asm20 -c -X -N 50 -p 0.1 --secondary=yes` all-vs-all, same
`mcl_families --min-exonic-bp 1 --min-shared-exon-frac 0.60`, same clustering.

⚠ **The exon conjunct still consults the annotation in every arm**, because a DNA region has no exons of
its own — `mcl_families --gff` computes exonic content for whatever interval it is given. That is why
this mode is *semi*-guided and not genome-only: **the nodes are annotation-free, the edge filter is not.**
Stating it up front so the arm is not later mistaken for a pure genome-only result.

## Truth and metrics

Truth = **Soto family IDs** (`bench/soto/soto_famCN_S1C.tsv`), restricted to the NPIP family, plus the
all-family chr16 number as context. ⚠**Symbol-root truth is VOID** (r902 — it fails on the development
set) and is not used.

Reported for every arm, per the standing reporting rule:
**sensitivity · precision · one-to-one bipartite matching F**, with collisions counted as misses.
Also reported because r817 says it is the failure mode: **how many truth genes COLLAPSE onto a shared
node**.

## The bar — committed now

| outcome | verdict |
|---|---|
| semi-guided bipartite **F beats BOTH** de novo and guided on NPIP, and collapse ≤ guided's | ⭐ **BETTER** — a real third mode |
| F within ±0.02 of the better of the two, collapse ≤ guided's | ⚠ **TIED** — no reason to prefer it |
| F below both, **or** collapse materially worse (the r817/§6ji failure reproducing) | ⛔ **NO** |

**Predicted outcome, stated before looking**: ⛔. SD intervals have a median of 7.1 kb and a max of
564 kb against a typical NPIP gene of ~15–30 kb, so one region will routinely swallow several genes.
Bipartite matching is one-to-one, so every swallowed gene beyond the first is a miss — which is exactly
the 21-of-31 collapse §6ji measured. If this prediction is wrong and F improves, that is the interesting
result and it must be explained by something other than "more sequence per node".

I will not change the node set, the truth, the metrics or the bar after seeing any number.

---

# OUTCOME (2026-09-22) — ⛔ **NO. Semi-guided beats de novo but loses to guided; guided stays best.**

Human chr16, Soto truth, one-to-one bipartite, identical downstream in every arm.
Two semi-guided variants because the rule forced a choice (see "correction" below): `exonic` = the
region's exons are the annotated exons falling inside it; `whole` = the region is one exon.

**NPIP (2 Soto families, 20 genes)**

| arm | nodes | clusters | sens | prec | **F** | collapsed |
|---|---|---|---|---|---|---|
| de novo | 2,550 loci | 7 | 0.600 | 0.923 | 0.727 | 2 |
| **guided** | annotated gene bodies | 2 | **0.750** | 0.938 | **0.833** | **0** |
| semi-guided `exonic` | 2,465 SD regions | 6 | 0.600 | **1.000** | 0.750 | 2 |
| semi-guided `whole` | 2,465 SD regions | 10 | 0.600 | **1.000** | 0.750 | 1 |

**All chr16 families (20 families, 71 genes)**

| arm | sens | prec | **F** |
|---|---|---|---|
| de novo | 0.465 | 0.971 | 0.629 |
| **guided** | 0.718 | 0.981 | **0.829** |
| semi-guided `exonic` | 0.535 | 0.745 | 0.623 |
| semi-guided `whole` | 0.535 | 1.000 | 0.697 |

⛔**Verdict: NO.** Semi-guided lands **between** the two existing modes on NPIP (F 0.750 vs de novo 0.727
and guided 0.833) and does not reach guided on any cohort. It is not a reason to add a third mode.

⭐**But the one thing it is best at is real and reproduces r817**: **precision 1.000 on NPIP in both
variants**, the highest of any arm — r817 saw the same shape (region nodes, NPIP L1 precision 0.875 vs
guided's 0.545). **Its ceiling is sensitivity, fixed at 0.600 — identical to de novo's** — because an SD
interval only exists where the genome is self-similar, so a family member in unique sequence has no node.
Guided's 0.750 comes from annotation covering members SD calls never propose.

## Two corrections to my own pre-registration

⚠**1. The predicted failure mode did NOT dominate.** I predicted coarse regions (median 7.1 kb, max
564 kb) would swallow several genes each and be punished by one-to-one matching. Measured collapse is
**2 (exonic) / 1 (whole) against de novo's 2 and guided's 0** — semi-guided is no worse than de novo
here. §6ji's 21-of-31 collapse did not reproduce on this substrate.

⚠⚠**2. My collapse metric was 0 BY CONSTRUCTION in the first run and I nearly reported it.** It counted
truth genes sharing a locus by iterating **locus → gene**, but the max-overlap resolver gives every locus
exactly one gene, so the count can never exceed 0. It has to be computed **gene → its best-covering
locus**. Fixed before any number above was recorded; the first (all-zero) column was discarded.

⚠**3. The prereg's mechanical assumption was wrong**: I wrote that `mcl_families --gff` "computes exonic
content for whatever interval it is given". It does not — it **joins node headers against annotation
records**, so SD regions matching no record made every node fall back to span and produced **0 clusters**
(register 899's "silent 0 nodes"). Both variants above exist because the fix required inventing an exon
model for a DNA region, which is itself a choice the mode has to make and the prereg had not anticipated.
