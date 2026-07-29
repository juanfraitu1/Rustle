# Over- and under-merging in the DNA family partition — diagnosis and levers

**Date:** 2026-07-27. Subject: `gw_family_catalog --from-genome` vs Soto 2025 (83 families / 362 members, CHM13v2.0).

## 1. The problem the member metric hides

Member sensitivity/precision are **partition-blind**. If two Soto families fuse into one DNA family, every
member is still inside a ≥2-copy family → still counted as recovered; and no off-catalog locus is created, so
precision barely moves either. Over-merging is therefore **invisible** to the headline numbers.

Measuring the partition instead (`bench/soto/partition_score.py`):

| metric | value |
|---|---|
| member recall | **361/362 = 99.7%** |
| Soto families with ≥1 member recovered | 82/83 |
| DNA families emitted (≥2 copies) | 69 |
| **clean 1:1** (DNA family = exactly one Soto family) | **38** |
| **over-merge** (one DNA family fuses ≥2 Soto families) | **18 families, fusing 48 Soto families** |
| **under-merge / split** (one Soto family across ≥2 DNA families) | **13 families** |
| homogeneity (DNA families that are "pure") | 74% |
| **completeness** (Soto families reproduced intact) | **46%** |

So: memberships reproduce at 99.7%, but **only 46% of Soto's families are reproduced as families.** Both numbers
must be quoted together — the first alone overstates the replication.

## 2. Mechanism — why families fuse (over-merge)

A homology edge forms when two loci share an aligned block clearing `min_identity` **and** `min_coverage`
(fraction of the **shorter** locus covered; default **0.50**). In segmental-duplication-rich regions
(pericentromeric / subtelomeric), **unrelated families share a duplicated block** — a mobile duplicon, an Alu,
a shared protein domain. Pairwise sequence cannot distinguish "same family" from "shares a duplicon", so the
quasi-clique fuses them. The code comments already name `min_coverage` as the intended defense
("Repeat bridges are held off by `min_coverage`").

Observed composition of over-merged families (baseline run):

| DNA family | Soto families fused | copies | on a Soto member | EXTRA (discovered) loci |
|---|---:|---:|---:|---:|
| GWFAM13 | 6 (FAR2P, ANKRD20, TEKT4P2, AC119751.3, AL669942.1, CR381670.1) | 53 | 13 | 40 (**75%**) |
| GWFAM5 | 5 (FAM72B, SRGAP2C, AC244669.2, FAM153B, POM121) | 34 | 9 | 25 (74%) |
| GWFAM57 | 4 (PMS2P4, POM121, NSUN5P2, SPDYE1) | 39 | 26 | 13 (33%) |
| GWFAM10 | 2 | 31 | 6 | 25 (**81%**) |

Over-merged families are **61%** EXTRA (non-Soto discovered) loci vs **53%** for pure families — and the
largest fusions are 74–81% extra. Suggestive that discovered paralog loci act as **bridges** between families;
decisive test = restrict the node set to the member windows only (below).

## 3. Mechanism — why families split (under-merge)

The 13 split families track two properties, both consequences of the same coverage rule:

| | split families | clean families |
|---|---|---|
| median(min member paralog identity) | 0.969 | 0.980 |
| **median size spread (max/min member bp)** | **5.8×** | 2.7× |

Examples: ID_215 spans **1,165 bp → 239,362 bp** members (205×), ID_481 **103 bp → 15,899 bp**, ID_182
**3,453 bp → 76,044 bp** at 0.889 identity.

Because `min_coverage` is a fraction of the shorter locus, **two large members that share only part of their
span fail the floor** and no edge forms → the family splits. Divergent members (identity below the floor) split
for the complementary reason.

## 4. The core tension

**The same knob moves the two error modes in opposite directions:**

- Raise `min_coverage` → kills duplicon bridges (**less over-merge**) → but more large/partial-overlap members
  fail the edge (**more splits**).
- Lower it → reconnects partial members (**fewer splits**) → but admits shared-block bridges (**more over-merge**).

A single global threshold therefore cannot fix both. This is the classic segmental-duplication
family-definition problem (Bailey 2002 → Eichler-lab work): **pairwise sequence similarity is not transitive
in the way family membership is.**

## 5. Levers

**A. Global thresholds (measured — `bench/soto/merge_lever_sweep.sh`).** `min_coverage`
(`RUSTLE_GENOME_MIN_COVERAGE`), `min_identity` (`--min-identity`), γ quasi-clique density
(`RUSTLE_GENOME_GAMMA`). These trade over-merge against split per §4; the sweep quantifies the curve.

**B. Restrict the node set (decisive test for §2).** Group only the given member loci — no discovered extras.
This is also the closest analog to Soto's own setup (they group *given* SD98 loci). If the fusions largely
disappear, the bridges are the discovered loci and the fix is to admit extras only when they do not join two
otherwise-disconnected blocks.

**C. Structure-aware de-bridging — ORTHOGONAL to the threshold trade, and already implemented in this
codebase but NOT wired into the `--from-genome` path:**
- `multi_repeat_bridge` — cut families joined only through a shared repeat unit.
- `recombinant_split` — split at articulation points (a single locus holding two blocks together).
- `bridge_detector` / `catalog_overlaps` — exon-graph bridge and shared-sequence diagnostics.

These act on **graph structure and repeat content**, not on the pairwise coverage of every edge, so they can
cut a bridge *without* penalizing legitimate partial-overlap edges — i.e. they attack over-merge without
causing the splits that raising `min_coverage` causes. **Highest-value next step.**

**D. Copy-number concordance (Soto's actual splitter) — irreducible here.** Soto refine shared-exon groups by
**famCN/parCN** (median read-depth copy number from WGS; group only when copy-number mean-absolute-deviation
< 1). That is an orthogonal modality: two families that share a duplicon have *different copy-number profiles*
even when their sequence bridges. We have no WGS depth for this benchmark, so this lever is unavailable — and
no amount of sequence-threshold tuning substitutes for it. This is the honest ceiling of a sequence-only
partition.

**E. External duplication-unit annotation.** Soto also anchor groups to DupMasker ancestral duplication units.
Equivalent to D in spirit: an orthogonal signal about *which* duplication a block belongs to.

## 6. Honest statement for the write-up

> Using the same pairwise-alignment + γ-quasi-clique engine on the genome, we reproduce Soto's family
> **memberships** (361/362 = 99.7% of paralogs) but only **46%** of their family **partition** — the shortfall
> is over-merging in SD-rich regions where unrelated families share a duplicated block, plus splits of families
> whose members differ greatly in length. Sequence-only thresholds trade one error for the other; separating
> them requires either structure-aware de-bridging (implemented, not yet wired) or the copy-number modality
> Soto used and this benchmark lacks.

---

## 7. RNA side: fragmentation diagnosis and the genomic-span fix (2026-07-27)

### RNA vs DNA have OPPOSITE failure modes (same metric)

| catalog | recall | homogeneity | over-merge | splits | completeness |
|---|---:|---:|---:|---:|---:|
| DNA `--from-genome` | 99.7% | 74% | 18 (fusing 48) | 13 | 46% |
| RNA `--cross-chrom --refine` | 65.5% | **91%** | 14 (fusing 28) | **45** | 25% |
| RNA `--homology-primary` | 77.3% | 86% | 32 (fusing 73) | 52 | 15% |
| RNA per-chrom conflict | 22.1% | **92%** | 6 (fusing 12) | 15 | 50% |

**RNA is markedly PURER (86–92% vs 74%) but far more FRAGMENTED (45–52 splits vs 13).** The intron-chain /
splicing constraint is real and free: the sequence that bridges unrelated families at the DNA level
(duplicons, intronic repeats, flanks) is exactly what splicing discards, so RNA reps cannot form those
spurious edges. **One mechanism, opposite effects: splicing removes both the sequence that DISTINGUISHES
near-identical copies (→ K=0) and the sequence that spuriously BRIDGES families (→ no over-merge).**

### Why RNA fragments — measured

Across the 45 split families, 1831 pairs should be one family but landed in different RNA families:
- **75% have NO exon-sum alignment at all** (even with the sensitive `-k11 -w5` tier);
- of the 449 that align: median identity **0.900** but median coverage-of-shorter **0.119**;
- relaxing the coverage floor to 0.05 recovers only **20%**.

It is **not** a threshold problem: the copies assemble **disjoint exon subsets** (one gets exons 1–2, its
sibling 7–8), so comparing *assembled* sequence has nothing to compare. Corroborating: split families have
**7.0×** median exon-count spread vs **1.5×** for clean ones, and **62%** of copies in split families are
single-exon fragments.

### The fix: compute the E_r edge on the GENOMIC SPAN of each RNA-detected locus

Same loci, genomic spans instead of exon-sums: median identity **0.900 → 0.973**, median coverage
**0.119 → 0.982**. Per split family, fully re-linked: exon-sum 5/45 (11%) → genomic 13/45 (29%) at cov≥0.30.

Simulated re-partition of all 863 RNA loci (connected components — merge-prone vs the engine's γ-quasi-clique):

| substrate | recall | splits | over-merge | homogeneity | completeness |
|---|---:|---:|---:|---:|---:|
| exon-sum (current) | 65.5% | 45 | 14 | 91% | 25% |
| genomic id≥0.90 cov≥0.50 | 64.9% | 32 | 11 | 90% | **43%** |

**Real-data validation — 4 chromosomes** (Soto regions only, `--homology-primary`, flag OFF vs ON):

| chrom | recall | families | clean 1:1 | over-merge | splits | **purity** | **completeness** |
|---|---|---:|---:|---|---:|---:|---:|
| chr7  | 77→79% ↑ | 16→13 | 1→1 | 6(15)→6(13) | 6→6 = | 62→**54** ↓ | 14→14 = |
| chr9  | 73→81% ↑ | 32→23 | 1→2 | 3(6)→3(6)   | 3→2 ↑ | 91→**87** ↓ | 14→25 ↑ |
| chr16 | 65→77% ↑ |  5→3  | 1→2 | 1(2)→1(3)   | 1→0 ↑ | 80→**67** ↓ | 33→40 ↑ |
| chr15 | 86→89% ↑ | 24→26 | 2→3 | 9(21)→8(19) | 7→7 = | 62→69 ↑ | 17→23 ↑ |

**Pooled (166 members):** member recall **78.3% → 83.1% (+4.8 pts)**; clean 1:1 5→8; splits 17→15;
over-merge 19 fams fusing 44 → 18 fusing 41; families **77 → 65** (fewer = more merging).

⚠**CORRECTION.** An earlier note claimed, from chr15 alone, that "every metric improves or holds". **That does
NOT generalize.** Across 4 chromosomes:

- **member recall improves 4/4** (robust, +4.8 pts pooled);
- **completeness improves 3/4, never worse**;
- **splits never worse** (2 improve, 2 flat);
- **but PURITY WORSENS in 3/4** (chr7 62→54, chr9 91→87, chr16 80→67). chr15 (62→69) was the exception.

So there is a **real trade, not a free lunch**: genomic spans buy recall and completeness at a cost in purity —
they reintroduce exactly the intronic/duplicon bridging sequence that splicing removes for free (§7 opening).
The worst purity loss is chr7 (SPDYE / GTF2I / POM121 / PMS2P — the most duplicon-dense region tested), which
also gained no completeness; testable prediction: **purity loss scales with regional duplicon density.**

**Verdict: keep `--homology-genomic-span` OPT-IN, not default.** It is the right lever when completeness/recall
matter more than partition purity (e.g. recovering fragmented families), and the wrong one in duplicon-dense
regions. A targeted variant — apply genomic spans only to loci whose assembly is incomplete (single-exon /
low exon-count fragments), keeping exon-sum edges elsewhere — would plausibly capture the completeness gain
without the global purity cost, and is the natural next experiment.

Shipped as `--homology-genomic-span` (opt-in; default OFF keeps every existing catalog byte-identical).
Needs **no read depth / copy number** — RNA decides WHICH loci are expressed, the reference supplies
complete sequence for WHAT-GROUPS-WITH-WHAT.

### §7b. Four remediation levers measured — none is a silver bullet (2026-07-28)

The split-away pieces of the 45 fragmented RNA families were used as the test set (233 loci). "Re-linked"
= the piece gains an edge/assignment back to its family.

| lever | families with ≥1 piece re-linked | measured cost | status |
|---|---:|---|---|
| exon-sum pairwise (**current default**) | 38% | — | shipped |
| lower the coverage floor (0.50→0.05) | ~+9% of *pairs* only | none | rejected: 75% of pairs have NO alignment at all |
| **projection** (family consensus → genome, spliced) | **49%** | not measured; spliced so no intron bridging expected | machinery exists (`--project-all-families`) but writes a separate file — hits are treated as copy-NUMBER, not membership |
| **genomic span** (`--homology-genomic-span`) | **58%** | **purity ↓ in 3/4 chromosomes** | shipped, opt-in |

Two approaches were **ruled out by code inspection** (minutes, not hours):

- **Genomic spans targeted at fragments.** Single-exon copies have genomic span ≈ their own length
  (median span/seq = **1.00×** vs 3.82× for multi-exon), so the substrate switch adds them no new sequence.
  The measured gain must come from MULTI-exon copies — the opposite of the intended target.
- **`family_rescue` (borrow-strength).** `thin_loci` starts with `if r.introns.is_empty() { continue; }`
  — unspliced reads are skipped — but **66% of RNA copies are single-exon**. It addresses multi-exon loci
  below the ≥3-read gate (a RECALL lever), not the fragment/split population.

**Why every failed lever failed the same way:** they all compare *the fragment's own assembled sequence* to
something, and a fragment covering exons 7–8 has nothing to offer that comparison. Projection is the only
one that inverts the direction (the family searches for the fragment), which is why it beats the current
default without touching purity.

**Conclusion.** RNA family fragmentation is **partially remediable, not solvable** by these levers:
38% → 49–58% of families gain a re-link, each with a trade. This is a real identifiability-adjacent limit of
the substrate, not a tuning failure — and it is the honest counterpart to the DNA side's over-merging.

### §7c. Projection edge tier — prediction refuted, and why that matters (2026-07-28)

Simulated a third edge substrate: **edge(A,B) when rep A's SPLICED consensus projects onto rep B's genomic
locus.** Asymmetric by construction — query contributes exons only, target contributes the complete locus.
**Predicted:** fragments get found (target is complete) without intron-to-intron bridging (query is spliced),
so completeness should rise *without* the purity cost genomic spans pay.

| substrate | recall | splits | purity | completeness |
|---|---:|---:|---:|---:|
| exon-sum (current default) | 65.5% | 45 | **91%** | 25% |
| genomic span id≥0.90 | 64.9% | 32 | 90% | **43%** |
| projection id≥0.90 | 65.5% | **22** | 84% | 33% |
| projection id≥0.95 | 65.2% | **22** | 85% | 38% |
| projection id≥0.98 | 62.2% | 29 | 89% | 40% |

⚠**The prediction is REFUTED.** At matched completeness (~40%), projection (purity 89%) and genomic spans
(purity 90%) sit on the *same* trade curve. The spliced query does not avoid the bridging: a spliced
consensus still aligns to a paralogous exon in an unrelated family, and the complete target locus supplies
the rest.

**Why this is the more valuable result.** Two structurally different mechanisms — symmetric
genomic-vs-genomic alignment, and asymmetric spliced-vs-genomic projection — land on the **same
completeness/purity frontier**. That is evidence the trade is **intrinsic to the problem**, not an artifact
of a particular edge substrate: any rule powerful enough to reach a fragment covering exons 7–8 is also
powerful enough to connect two families that share those exons. Escaping it needs an **orthogonal signal**
(copy number / parCN), not a better sequence rule — the same conclusion the DNA side reached independently.

**Secondary finding:** projection is the strongest *split* reducer measured (45 → **22**, vs 32 for genomic
spans), so it is the better choice if de-fragmentation specifically is the goal and ~6 points of purity is
acceptable.

**Decision: do not implement.** `--homology-genomic-span` is already shipped and validated on 4 chromosomes;
projection is an equivalent point on the same curve, so a second lever adds engineering without new
capability. Recorded here so the equivalence is documented rather than re-derived.

## 8. Bipartite size matching — do predicted copies have the RIGHT SIZE? (2026-07-28)

Requested by the advisor. Every number up to here scores *placement* — is a truth member covered by some
predicted copy? That is a many-to-many test and it is blind to size: one over-extended prediction spanning
two members scores the same as two correctly-sized ones.

**Method** (`bench/soto/bipartite_size_match.py`). Force an optimal **1:1** assignment between predicted
copies and Soto truth members — Hungarian algorithm (`scipy.optimize.linear_sum_assignment`) maximising
total **reciprocal overlap**, `overlap / max(len_pred, len_true)`. Reciprocal (not raw) overlap is the point:
it reaches 1.0 only when the two intervals *coincide*, so the objective rewards getting the size right rather
than merely touching the right place. The forced 1:1 then makes the leftovers explicit:

| pattern | reading |
|---|---|
| matched, ratio ≥ 2× | prediction OVER-EXTENDED — over-merge signature |
| matched, ratio ≤ 0.5× | prediction TRUNCATED — fragmentation signature |
| unmatched prediction | spurious / extra copy |
| unmatched truth | missed member |

### Result — RNA (`definitive.copies.tsv`, 863 copies vs 362 members)

```
MATCHED 1:1        232  (64.1% of truth)
unmatched truth    130   (missed)
unmatched pred     631   (extra)

size ratio (pred/true) over matched pairs
  median             0.54       median log2 bias  -0.88   (UNDER-prediction)
  IQR                0.19 - 1.00
  within 2x          102  (44% of matched)
  OVER-EXTENDED >=2x  17
  TRUNCATED    <=0.5x 113
```

**RNA's dominant size error is TRUNCATION, not over-merging** — 113 truncated vs 17 over-extended, median
0.54×. The extremes are large: SRGAP2C 5.6 kb predicted vs 208 kb true (0.03×), GUSBP1 0.03×, NOTCH2 0.03×.

This is the same fragmentation already diagnosed in §7, now measured on the size axis instead of the
partition axis, and it is *mechanistically expected*: a truth member's span is the whole genomic locus,
while an RNA copy spans only the transcribed, covered, assembled portion. It is worth stating plainly for
the write-up — the advisor's concern is over-merging, and on RNA the measured bias runs the other way.

### ⚠ The DNA row is CIRCULAR — do not quote it

`--from-genome` scores 348/362 matched, median ratio **1.00**, IQR 1.00–1.00, 0 over-extended, 0 truncated.
That is not a result. DNA mode seeds its rep nodes **from the Soto member windows**, so predicted intervals
are the truth intervals by construction and the size ratio is 1.00 by definition. The metric is meaningful
only for RNA (and, separately, for DNA mode's *discovered extra* loci, which have no truth interval to
compare against). Reporting the DNA row as a 100% size agreement would be exactly the circularity already
corrected once in `project_soto_honest_pr_artifacts`.

### What it does NOT say

The 631 unmatched RNA predictions are not 631 false positives — the catalog covers the whole Soto regions,
including real paralogs outside Soto's 80-family subset (§6, and the over-detection audit in
`project_soto_recall_recompute`). The 1:1 constraint also *forces* leftovers whenever a truth member is
genuinely covered by two fragments; that is counted here as one match plus one "extra", which is the
intended reading (fragmentation), not an independent error.

### 8b. Control — is the "truth" a duplicon rather than a gene? (tested, and it is not the explanation)

Fair objection: Soto's member intervals come from the SD98 self-alignment, so they are **duplication blocks**,
not gene models. A block can be much larger than the transcribed gene inside it, and an RNA copy can only
ever span the gene — so a "0.03x truncated" call might be comparing a transcript to a duplicon.

`bench/soto/gene_preferred_truth.py` rebuilds the truth: where a NAMED RefSeq gene matching the member's gene
label overlaps the member interval, the gene's span replaces the block; otherwise Soto's block is kept.

```
members 362 | matched a named annotation 222 | no named match (Soto kept) 140
interval CHANGED 222  (shrunk 99, grown 123)
block->gene size ratio: median 1.00   p10 0.77   p90 1.82
largest shrinks: GTF2IP12 49129->5058 (0.10x), CTSLP3 6137->977 (0.16x), HERC2P9 95778->30825 (0.32x)
```

Re-scoring against the gene-preferred truth barely moves anything:

| | vs Soto blocks | vs gene-preferred |
|---|---|---|
| median size ratio | 0.54 | 0.54 |
| TRUNCATED <=0.5x | 113 | 112 |
| OVER-EXTENDED >=2x | 17 | 16 |

Individual intervals *do* shift (the effect is real for a minority — the p10/p90 spread is 0.77-1.82), but the
median block/gene ratio is 1.00, so the correction cancels in aggregate. And the extreme truncations survive
with annotations attached: SRGAP2C gene = 207.9 kb (block 208.1 kb), NOTCH2 = 158.1 kb (block 189.2 kb),
GUSBP1 = 231.3 kb (block 346.6 kb). These are genuinely large, intron-rich genes; a 5.6 kb prediction for
SRGAP2C is a fragmentary assembly, not a units mismatch.

**Conclusion: the truncation is real.** The size axis independently confirms the fragmentation diagnosed in
§7 from the partition axis. The gene-preferred BED is kept as the fairer default truth for future size
scoring even though it does not change the headline.

## 9. Conditional TSS extension — implemented, measured, and NOT recommended (2026-07-28)

Requested as a way to give copies "the adequate size" and so reduce over/under-merging. Implemented as
`snap_boundary` + the `RUSTLE_TSS_SNAP` opt-in in `denovo_assemble.rs`: where a skeleton's read starts (or
ends) are SHARPLY PEAKED — >= 30% inside one 400 bp window, the criterion validated in §8 of
`bam_tie_signals.md` — the boundary snaps to that peak's outer edge instead of the k-th-read quantile.
Default OFF, so every existing catalog is byte-identical.

**Two implementation bugs found and fixed during validation** (both worth recording, since the first would
have shipped silently):

1. Snapping to the peak's **median** orphaned every read in the peak's first half — a guaranteed ~200 bp
   truncation on every snap. On chr1 it moved 8/43 boundaries (max 756 bp); the corrected **edge** rule moves
   2/43, all <= 20 bp. Fixed.
2. A doc comment claimed the snap "never shortens". It does — pulling an outlier-dragged boundary inward IS
   shortening, and that is the whole mechanism. Comment corrected to state the cost.

**Measured effect** — same binary, same BAMs, same `--refine`, boundary rule the only difference, scored by
the §8 bipartite matcher against the §8b gene-preferred truth (chr1+chr7+chr15+chr16, 94 vs 95 copies):

| | snap OFF | snap ON |
|---|---|---|
| matched 1:1 | 23 | 23 |
| median size ratio | 0.35 | 0.35 |
| **within 2x** | **10 (43%)** | **8 (35%)** |
| **TRUNCATED <=0.5x** | **13** | **15** |
| size-ratio IQR | 0.08 – 0.71 | 0.07 – 0.60 |

**The extension makes size agreement WORSE, and the reason is structural rather than statistical.** A snap
can only ever pull a boundary *inward* (it fires when the quantile sits outside the peak), so it can only
shorten. RNA copies are already systematically too short — median 0.54x in §8, 0.35x on this subset. Making
them shorter moves them away from the truth. The lever pushes on the wrong side of the error.

It also has a detection cost: on chr7 it drops `DN_chr7_74502675_10` (10 exons, 105 reads) entirely, because
the shortened exon-sum falls below `--refine`'s cov >= 0.50 gate and the copy loses its family edge. chr16
gains 2 copies (15 -> 17), so the copy-count effect is close to a wash; the size effect is not.

**Status: shipped as opt-in, default off, not recommended.** The code and its tests stay because the
measurement is the deliverable — the honest answer to "should the exon-sum encode TSS?" is now backed by
three independent results: only 3/40 copy pairs have a genuinely distinct TSS (`bam_tie_signals.md` §8),
snapping the boundary degrades size agreement (this section), and the size error that actually matters is
truncation from fragmentary assembly (§8), which no boundary rule can fix.
