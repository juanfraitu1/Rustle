# Pre-registration — do the §6w0 split triggers recover the 406 COLLATERAL evictions?

**Written 2026-09-22, §6x3, before any arm is scored.** User goal: *"improve node definition to avoid
false positives and false negatives in all modes."* This executes **priority 1 of
`docs/NODE_GRAPH_ADMISSION_2026-09-21.md`** verbatim: *"Re-score the existing split triggers on graph
admission, not per-copy correctness. r968's turnover trigger moved per-copy +0.78pp and was judged 4×
short of its bar; its effect on the 406 collateral evictions was never measured. This is a re-scoring of
work already done, not new machinery."*

## Why this is the highest-value open item, in one paragraph

§6w2/r970-972: on chr16, **1,892 de novo loci have real homology and only 864 become graph nodes**;
**99.8% of the 1,026 dropped fail `cov_longer < 0.30` ALONE**, and **406 of them (39.6%) are evicted by a
partner holding ≥2 whole genes**. So a node-level FALSE POSITIVE (over-merge) manufactures node-level
FALSE NEGATIVES in innocent third parties, at ~2.4× its own count — `cov_longer`'s denominator is the
LONGER locus's whole span. §6w6 then refuted boundary pull-in as the remedy (shrinking a locus costs it
its own alignments) and concluded: *"Any future attempt must shrink the giant WITHOUT shortening the
sequence that earns its own edges — i.e. a node SPLIT, not a boundary pull-in."* The split arms exist and
have never been scored on this metric.

## Arms — all three already on disk, nothing is rebuilt

| arm | trigger | reads flagged | artefacts |
|---|---|---|---|
| **A0** baseline | none | — | `/mnt/linuxdisk/tmp/regress/dn16.{paf,graph.tsv}` |
| **T1** chimeric bridge (r967) | disjoint-span per primary read | 63 | `/mnt/linuxdisk/tmp/idealsim/chr16_t1split_loci.{paf,graph.tsv}` |
| **T2** read-identity turnover (r968) | read set turns over across a junction | 1,926 | `…chr16_t2split_loci.{paf,graph.tsv}` |

Scored with `bench/node_graph_admission.py` unchanged (its reimplementation matched the shipped graph
866 vs 864), `cov_longer` on the **deferred** path.

## ⚠⚠ The decomposition that makes this falsifiable

A split produces two SHORTER loci. A shorter locus is a smaller `cov_longer` denominator whenever it is
the longer partner, **so admission can rise purely from fragmentation with no innocent partner rescued.**
Admitted-node count alone therefore cannot answer the question. Loci are identified by
`chrom:start-end`, so the arms decompose exactly:

- **PRESERVED locus** = identical `chrom:start-end` in A0 and the split arm — untouched by the split.
- ⭐**RECOVERED COLLATERAL** = a PRESERVED locus **dropped in A0 and admitted in the split arm**. This is
  the mechanism under test and the only number that counts.
- ⚠**SPLIT-PIECE GAIN** = an admitted locus whose coordinates do not exist in A0. This is the free-lunch
  component and is reported separately, never pooled into the headline.
- ⛔**COLLATERAL LOSS** = a PRESERVED locus admitted in A0 and dropped in the split arm. Netted against
  recovery; a split that evicts as many innocents as it rescues has done nothing.

## The bar — committed now

Judged on **net preserved-locus admission** (recovered − lost), against §6w2's 406 collateral evictions:

| outcome | verdict |
|---|---|
| net ≥ **+40** preserved loci (≈10% of the 406), split-piece gain reported but not required | ⭐ **DIRECTION VALIDATED** — node splitting is the right remedy, pursue it |
| net **+10 to +40** | ⚠ **PARTIAL** — real mechanism, too small as specified |
| net **< +10**, or the total admission gain is ≥90% split pieces | ⛔ **NO** — admission rose by fragmentation, not by relieving over-merge |

**Predicted, before looking — ⛔, and the reason is already on record.** r967 and r968 BOTH report that
**zero of the flagged reads fell within any of the three known over-merge sites** (NPIPB4/RRN3P1's shared
locus, CDR2's 91 kb engulfing locus, PKD1P6-NPIPP1). A trigger that does not fire at the over-merge sites
cannot relieve those sites' denominators, so the collateral evictions they cause must survive. If the
prediction is wrong and net recovery is large, then the collateral evictions are driven by *many small*
over-merges rather than the three big ones — which would itself be the finding, and would redirect the
whole line away from the giants.

⚠**This measures the de novo mode only** — guided and semi-guided nodes are annotation/SD intervals and
are not produced by the assembler, so no read-level split applies to them. The "all modes" half of the
goal is addressed separately by whatever this establishes about the shared `cov_longer` denominator.

I will not change the arms, the decomposition, or the bar after seeing any number.

---

# OUTCOME (2026-09-22) — ⛔ **NO, on both clauses of the bar. And splitting is STRUCTURALLY incapable, not merely mis-tuned.**

A0 reproduced §6w2 exactly (2,550 loci → 1,892 homologous → 866 admitted / 864 shipped; 1,026 dropped,
99.8% on `cov_longer` alone; 406/1,026 = 39.6% collateral).

| arm | admitted | preserved loci | RECOVERED | LOST | **NET** | new coords (split pieces) |
|---|---|---|---|---|---|---|
| A0 | 866 | — | — | — | — | — |
| T1 chimeric bridge | 869 | 1,841 | 2 | 4 | **−2** | 32 |
| T2 turnover | 885 | 1,782 | 7 | 7 | **0** | 71 |

⛔**Every admitted node either arm gains is a NEW coordinate.** Net preserved-locus admission is 0 (T2)
and −2 (T1) against a +40 bar — and the second ⛔ clause fires too, since 100% of the gain is split
pieces. The eviction structure is untouched: dropped share 54.2% → 55.1%, collateral share 39.6% → 40.1%.
**The +19 admitted nodes are fragmentation, not rescue.**

## The failure is double, and only the second one is interesting

1. **Targeting.** Of A0's **380 distinct culprit loci**, **91.3% (T2) / 95.0% (T1) survive with their
   coordinates completely unchanged** — the triggers fire almost entirely somewhere else.
2. ⭐**Mechanism.** Where a culprit IS split (T2: 33 culprits, 56 victims still present), **1 of 56
   victims is rescued — 1.8%.** Splitting the giant does not relieve its victims.

## ⭐⭐ Why — the required shrink is quantified, and a split cannot deliver it

For each of the **679 loci evicted by a longer partner whose identity and `alen` already pass**, the
culprit must fall below `overlap / 0.30` for the victim to be admitted:

| | |
|---|---|
| median culprit length | **64,113 bp** |
| median length it must fall below | **13,870 bp** |
| ⭐**median shrink factor required** | **4.83×** |
| 25th / 75th / 90th pct | 2.23× / 12.30× / 22.90× |
| victims reachable by a <2× shrink | **147 = 21.6%** |
| victims needing **>5×** | **335 = 49.3%** |

**A binary split delivers at most 2×.** So even a perfectly targeted, perfectly placed split reaches at
most the 21.6% of victims that need under 2× — and half of them need more than 5×, which no split into
two pieces can ever reach. ⛔**This closes read-level splitting as a remedy for collateral eviction**, and
it closes it by arithmetic rather than by one more refuted trigger. Combined with §6w6 (boundary pull-in
costs the locus its own alignments) **both ways of shrinking the giant are now dead.**

## ⭐⭐⭐ What the same numbers say the answer IS

The victims are **not partial alignments**: median overlap 4,161 bp against a median victim length of
5,208 bp, and the **median per-pair ratio is 1.00** — the evicted locus aligns along essentially its
*entire* length. Nothing is wrong with the victim, the alignment, or the boundary. It is rejected only
because `cov_longer` divides by the **longer** locus's span.

⭐**That is the containment/asymmetry problem of §6t7–§6u0, appearing at the ADMISSION layer, where its
solution has never been applied.** §6u3/r915 measured the blind spot (25% of true pairs are asymmetric,
jaccard/ochiai/dice all score 0.0% on them, containment 91.5%); §6t7/r916 found **guarded containment**
beats the shipped metric (F .6648 vs .6394, precision .586 → .650); §6u0/r924 established the guard that
makes it safe (`--exonic-both-sides`, exonic fraction .367 TRUE vs .107 FALSE). All of that was developed
and scored on **edge weights inside an existing graph**. `cov_longer ≥ 0.30` at admission still uses the
bare longer-span denominator — the exact norm r913 showed becomes a hub generator when it is
**un**guarded.

⚠**Not a licence to swap in `min(la,lb)`**: r913 refuted the unguarded containment norm (864 pairs ≥1.0,
largest component 171 vs 23, precision .164 vs .639), and r917 recorded that r916's guarded win was
measured in `mcl_port.py`, not bit-identical Rust, so it is not yet a basis to change the definition.
The next arm is **guarded containment at admission**, pre-registered separately and scored on graph
admission plus a hub guard.

## Prediction scorecard

⛔ predicted, ⛔ measured — but **my stated reason was wrong.** I predicted failure because r967/r968 report
that no flagged read fell in the three known over-merge sites. A0's own output contradicts the premise:
**380 distinct culprits, top-10 share 15.6% — "diffuse, not a few monsters."** The three giants were never
the target. The real reasons are the two above, neither of which I anticipated, and the arithmetic one
(4.83× median) is the durable finding.

---

# ADDENDUM (§6x6) — the premise of this whole line was wrong

After the split arms, the containment escape, and a refuted overlap guard, one framing was still untested:
attack the **denominator locally** rather than the giant's **sequence**. Measured with the annotation as an
ORACLE (never a pipeline input): give each collateral eviction a denominator equal to the ONE constituent
gene its alignment lands in.

**Ceiling 13.0% (47 of 362) — and the reason kills the premise.** Median fused locus **125,450 bp**, median
constituent gene **97,026 bp**, **median denominator shrink 1.02×**. ⭐⭐**There is no inflated denominator
to recover: the "fused giant" is ~98% one genuinely large gene.** The victim (median 7,824 bp) fails
because 7,824/97,026 = 0.08.

⛔⛔**So §6w2/r971's causal claim is wrong.** `gene_counter >= 2` finds loci that *contain* two or more
annotated genes, but the alignment lands in one that dominates the span. **Node false positives
(over-merge) and node false negatives (eviction) are largely INDEPENDENT.** That is why every remedy
failed: §6w6 pull-in, r1001 splitting (4.83× needed), r846 node cut, r1012 local denominator (1.02×
available) are four independent refutations of a single mis-attributed cause.

⭐**The correctly-shaped fix is to normalise by the SHORTER locus** — which is exactly `--min-cov-shorter`,
and why it is the one lever in this line that moved anything.
