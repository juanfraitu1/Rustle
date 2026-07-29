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
MATCHED 1:1        215  (59.4% of truth, containment >= 0.50)
grazing pairs       17   (overlap > 0 but neither interval 50% covered -- NOT counted)
unmatched truth    147   (missed)
unmatched pred     648   (extra)

size ratio (pred/true) over matched pairs
  median             0.54       median log2 bias  -0.89   (UNDER-prediction)
  IQR                0.17 - 0.97
  within 2x           99  (46% of matched)
  OVER-EXTENDED >=2x  12
  TRUNCATED    <=0.5x 104
```

CORRECTION (2026-07-28, found by adversarial review). The first version of this table admitted ANY overlap
> 0 as a 1:1 match. That let a GRAZE count as a match: prediction chr1:144357575-144366404 clips truth member
PPIAL4D (chr1:144366357-144366853) by 47 bp -- 0.5% of the prediction, 9.5% of the member -- and was reported
BOTH as a matched member and as a 17.8x "over-merge signature", because the ratio uses full interval lengths
regardless of how little lies in the intersection. 18 of the original 232 pairs had NEITHER interval 50%
covered; the smallest admitted overlap was 33 bp. The matcher now requires the intersection to cover >= 50%
of one interval (`--min-contain`, default 0.50). Superseded figures were: matched 232 (64.1%), within-2x 102,
over-extended 17, truncated 113.

The conclusion is unaffected and the median is threshold-stable:

| --min-contain | matched | median ratio | within 2x | over | trunc |
|---|---|---|---|---|---|
| 0.00 (old) | 232 | 0.54 | 104 | 16 | 112 |
| 0.25 | 229 | 0.54 | 102 | 16 | 111 |
| **0.50 (default)** | **215** | **0.54** | **99** | **12** | **104** |
| 0.75 | 201 | 0.52 | 90 | 12 | 99 |

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

**Tail risk of the substitution, and why it does not matter here.** The corrections are NOT one-directional
as the premise assumed: 99 shrink, 122 grow, 1 moves without resizing. Some growths are extreme and are
genuine RefSeq model differences rather than matching bugs -- H3P4's Soto block is 417 bp (chr1:121127576-
121127993) while RefSeq annotates a 58.9 kb *pseudogene* whose SD-homologous core that block is (141x);
UBE2Q2P8 is 103 bp vs a 7.6 kb pseudogene (74x). Which interval is "the truth" is then a modelling choice,
not a fact. The report now prints BOTH tails so the premise cannot confirm itself.

Because it is a choice, every conclusion drawn from it was re-run under both truths. They agree: the §8
headline is unchanged (median 0.54 either way, table above), and the §9 paired test is *identical* under both
(matched 21, 6 changed, 4 worse / 2 better, sign-test p = 0.688 for Soto blocks AND for gene-preferred). No
conclusion in this document depends on the truth definition.

**Conclusion: the truncation is real.** The size axis independently confirms the fragmentation diagnosed in
§7 from the partition axis. The gene-preferred BED is kept as the fairer default truth for future size
scoring even though it does not change the headline.

## 9. Conditional TSS extension — implemented, measured, and left OFF (2026-07-28)

Requested as a way to give copies "the adequate size" and so reduce over/under-merging. Implemented as
`snap_boundary` + the `RUSTLE_TSS_SNAP` opt-in in `denovo_assemble.rs`: where a skeleton's read starts (or
ends) are SHARPLY PEAKED — >= 30% inside one 400 bp window, the criterion validated in §8 of
`bam_tie_signals.md` — the boundary snaps to that peak's outer edge instead of the k-th-read quantile.
Default OFF, so every existing catalog is byte-identical.

**Measured effect** — same binary, same BAMs, same `--refine`, boundary rule the only difference, scored by
the §8 bipartite matcher against the §8b gene-preferred truth (chr1+chr7+chr15+chr16, 94 vs 95 copies). The
matcher selects the SAME 23 truth members in both catalogs, so the comparison can and must be PAIRED:

| | snap OFF | snap ON |
|---|---|---|
| matched 1:1 | 23 | 23 |
| median size ratio | 0.35 | 0.35 |
| within 2x | 10 (43%) | 8 (35%) |
| TRUNCATED <=0.5x | 13 | 15 |
| size-ratio IQR | 0.08 – 0.71 | 0.07 – 0.60 |

**Paired analysis (the honest test):** only **6 of 23** matched members change at all — 4 worse, 2 better.
Sign test p = **0.69**; Wilcoxon signed-rank on |log2 size error| p = **0.16**. Mean |log2 error| 2.065 (off)
vs 2.182 (on). The six are three sibling pairs from three families: AMY1B/AMY1C, POM121/POM121C,
NSUN5P1/NSUN5P2. **There is no detectable effect on size agreement in either direction.**

### CORRECTION (2026-07-28) — two claims in the first version of this section were wrong

Both were caught by adversarial review, and both had been asserted confidently:

1. **"A snap can only ever pull a boundary inward, so it can only shorten."** FALSE. `snap_boundary`'s near
   branch returns `fallback.min(peak)` for a start and `fallback.max(peak)` for an end, so when the peak lies
   OUTSIDE the quantile the boundary moves outward and the transcript gets LONGER. Measured over the four
   chromosomes: **12 copies lengthen, 15 shorten.** This was the load-bearing mechanism for the original
   "pushes on the wrong side of the error" conclusion, and it does not exist.
2. **"The extension makes size agreement WORSE."** NOT SUPPORTED. The unpaired 10-vs-8 and 13-vs-15 counts
   looked like a trend, but the paired test above (p = 0.69) shows the entire difference is 4 vs 2 discordant
   pairs among 6 changes in 3 families. Quoting a within-2x drop of "43% → 35%" from n=23 without pairing was
   the error; the same data, paired, says nothing.

The corrected reading is weaker and more boring: **on this benchmark the snap is close to a no-op.** It
touches 6 of 23 scored members, moves size agreement neither way at p = 0.69, and the copy-count effect is
roughly a wash (chr7 loses `DN_chr7_74502675_10`, 10 exons / 105 reads, when the changed exon-sum falls below
`--refine`'s cov >= 0.50; chr16 gains 2, 15 → 17).

**Status: shipped as opt-in, default off, not enabled.** The reason is now *absence of demonstrated benefit*,
not demonstrated harm — a weaker but defensible basis for leaving a lever off. What genuinely holds against
encoding TSS in the exon-sum is the other two results, neither of which this correction touches: only 3/40
copy pairs have a distinct TSS at all (`bam_tie_signals.md` §8), and the size error that actually matters is
truncation from fragmentary assembly (§8), which no boundary rule can fix.

**Implementation bugs found while validating** (recorded because the first would have shipped silently):

1. Snapping to the peak's **median** orphaned every read in the peak's first half — a guaranteed ~200 bp
   truncation on every snap. On chr1 the median rule changed 8 of 43 boundaries (max 756 bp); the corrected
   **edge** rule changes 5 of 43, max 545 bp. (An earlier draft said "2/43, all ≤ 20 bp" — that came from a
   partial comparison keyed on an unstable column, and is corrected here.)
2. A doc comment claimed the snap "never shortens"; the replacement, added when that was corrected, claimed
   it "necessarily SHORTENS". Both are wrong — it does both. Source and unit-test comments reconciled.

## 10. `--cross-chrom` chi(H) guard was inert — fix and impact measurement (2026-07-28)

Found while chasing a failing integration test that `cargo test --lib` never runs.

**The bug.** `detect_conflict_catalog_genome_wide_xchrom` never populated `DenovoTranscript::
distinguishing_uniq` (its `reps` binding was not even `mut`), while the two other read-based catalog paths do
(`denovo_pipeline.rs:2122`, `:2457`). With the field 0 on every rep, `reads_distinguish(0, 0, ..)` is always
false, so `distinct_locus_reps` collapsed ANY overlapping same-strand pair **unconditionally**. The
read-evidence guard added in 9e887b4 — the thing that makes co-located copies merge only when reads genuinely
cannot tell them apart, i.e. the chi(H) identifiability argument in code — was **permanently inert on the
`--cross-chrom` path**. It became visible only when 121b7ea made `--refine` default-on.

**Why it went unnoticed for weeks.** Two test defects, one masking the other:
- `default_cross_chrom_output_is_unchanged` pointed `--out` at its own **tracked golden files**, so every run
  overwrote the baseline it was meant to compare against and then asserted only `lines().count() > 1`. It
  could not fail. The committed golden still contained the same-chrom family `GWFAM1` and would have caught
  the regression on the first real comparison.
- `cross_chrom_catalog_emits_same_chrom_family` did fail — but it is an integration test, and the routine
  suite command is `cargo test --lib`, which does not run `tests/`.

**Evidence the fix is right.** With the wiring restored, the committed golden `out_default.*` — which
predates the regression — is reproduced **byte-for-byte**, and both integration tests pass.

**Impact on the Soto benchmark: NEUTRAL.** Measured on the real cross-chrom pass (`--cross-chrom --refine`,
the 18 cross-chrom families, 1.9 GB BAM), against the cached pre-fix run:

| | pre-fix | post-fix |
|---|---|---|
| xchrom copies | 284 | 285 |
| xchrom families | 85 | 86 |
| loci only in one side | 14 old | 15 new |
| **combined member recall** | **95/362** | **95/362** |

The member score is **identical** — verified by scoring the pre-fix combined catalog with the same script as
a control, which returns the same 95/362. So the fix is a correctness repair with no measured movement in
either direction: it does not rescue members, and — contra the concern that restoring the guard would
inflate copy counts by splitting duplicate transcript models — it adds exactly one copy and one family.

⚠ **Separate pre-existing discrepancy, NOT caused by this fix.** The per-chrom combined catalog
(`perchrom_catalog.copies.tsv`) scores 95/362 = 26.2% under `soto_cache_score.py`, while
`definitive.copies.tsv` scores 237/362 = 65.5% and the committed headline is 276/362 = 76.2%. By the score
script's own validation line ("recall near 76.2% => cache+recipe reproduce the genome-wide result; far off =>
not comparable") the per-chrom recipe is **not comparable** to the headline. Both numbers predate this
session's work. The headline figures quoted elsewhere derive from the `definitive` recipe; the per-chrom
cache should not be used for recall claims until that gap is explained.

## 11. The per-chrom recipe scored 26.2% against a 76.2% headline — root cause and fix (2026-07-28)

Flagged in §10 as a pre-existing discrepancy; chased down here. Two independent causes, both real.

### Cause 1 (primary): the recipe ran a DIFFERENT ALGORITHM from the headline

The headline `definitive` catalog — 863 copies, 280 families, the 76.2% member figure — is **one
`--cross-chrom` pass over the whole `soto_regions.bam`**: Louvain community split on a *global* conflict
graph (`detect_conflict_catalog_genome_wide_xchrom`). `recompute_perchrom.sh` instead ran
`gw_family_catalog`'s **default** mode on each chromosome — per-chrom connected components emitted as "clean
families (same-strand, disjoint-loci, >= 2 copies)" (`gw_reps_and_catalog`) — and its header asserted this
"reproduces the genome-wide within-chrom detection EXACTLY".

That assertion was false, and it is a much more conservative algorithm. Measured per chromosome
(default -> `--cross-chrom`, with `definitive`'s own per-chromosome breakdown for reference):

| chrom | default | --cross-chrom | definitive |
|---|---:|---:|---:|
| chr1 | 45 | 180 | 163 |
| chr15 | 18 | 132 | 123 |
| chr9 | 6 | 98 | 97 |
| chr7 | 19 | 65 | 65 |
| chr10 | 8 | 29 | 30 |

Per-chrom units now run `--cross-chrom`. Within a single-chromosome BAM that path can only form same-chrom
families, so the parallel split stays valid; genuinely cross-chrom families keep their separate pass.

**Result: 95/362 = 26.2% -> 257/362 = 71.0%** (1107 copies vs the old 506). For reference `definitive`
itself scores 237/362 = 65.5% on the same scorer, so the corrected split is not merely comparable — it is
slightly *ahead* of the single-pass catalog, which is expected: isolating a chromosome removes the global
graph's opportunity to absorb its loci into cross-chrom communities.

### Cause 2: the 76.2% headline is a 4-LEG UNION; only one leg is being rebuilt

`soto_cache_score.py` unions four detection legs. Only the first is recomputed here; the other three are
cached artifacts that are **absent from the cache directory entirely**:

```
leg RNA-split    : 1107 copies      <- the rebuilt catalog
leg protein-tail : MISSING          soto_gw_prot.copies.tsv
leg projection   : MISSING          soto_pall.allproj.tsv
leg expr-collapse: MISSING          soto_ce.expressed_collapsed.tsv
```

So the residual 71.0% vs 76.2% gap is **not** a recipe defect: it is one leg measured against a four-leg
union. Any recall number from this cache is an RNA-split-only number and must be labelled as such.

### Cause 3 (found while fixing): the runner reported success for a total failure

`_detect_unit.sh` captured `rc=$?` *after* an `echo`/`wc` pipeline, so a binary that died instantly reported
`rc=0`; and stale `*.copies.tsv` from the previous run survived the failure, so the combine step re-globbed
them. A recompute in which **all 21 units failed on a missing FASTA** therefore printed "all units done in
17s" and reproduced the previous run's numbers exactly. Now: the outputs are removed before each run, the
binary's rc is captured immediately and propagated, and the hardcoded FASTA path (a manual WSL mount that is
usually absent) falls back to the Desktop copy and exits 2 if neither exists.

**Standing caution:** a benchmark harness that cannot fail is worse than no harness. All three causes here
were silent — a false comment, a missing-leg union printed as "MISSING" but not as an error, and a runner
that swallowed 21 consecutive crashes.

## 12. What the truncation actually is — spliced vs unspliced (2026-07-28)

§8 reported "RNA copies are systematically truncated, median 0.54x". That number is a **mixture of two
populations with opposite behaviour**, and quoting it unstratified is misleading:

| matched pairs | n | median ratio | within 2x | truncated <=0.5x | over >=2x |
|---|---:|---:|---:|---:|---:|
| **SPLICED (>= 2 exons)** | 127 | **0.77** | 84 (66%) | 36 | 7 |
| **UNSPLICED (1 exon)** | 91 | **0.17** | 15 (16%) | 72 | 4 |

67% of all truncated pairs are single-exon copies (vs 17% of well-sized ones), and the ratio degrades
monotonically with exon count: 1 exon -> 0.133, 2 -> 0.263, 3 -> 0.096, 4 -> 0.268, 5+ -> 0.361.

### Why unspliced copies cannot span a locus — and why that is structural, not a bug

Unspliced copies are produced by `cluster_unspliced` in `pass1_skeletons_robust`: single-linkage
**span-overlap** clustering of unspliced reads per chromosome. A cluster's genomic extent is therefore
bounded by **read length** (~1-10 kb), not by gene length, no matter how long the locus is. Nothing can join
two unspliced read-clusters 50 kb apart, because an unspliced read carries no intron chain — the very
evidence that lets spliced reads reach across a locus. GUSBP1's 231 kb annotated span is covered by **11
separate ~9 kb unspliced clusters**; CNTNAP3P2's 238 kb by 19.

Merging them by proximity is possible but is exactly the readthrough/over-merge failure the pipeline
deliberately guards against (the readthrough filter already drops single-exon transcripts engulfing >= 5
junctions, `project_readthrough_filter`).

### The hypothesis this refutes, and the one it supports

- **"Spliced RNA inherently cannot span the whole locus"** — REFUTED, and inverted. A spliced transcript's
  genomic span *includes* its introns, so splicing does not shrink it. The spliced copies are the
  well-sized population (0.77). Full-length reads demonstrably reach the whole locus: at SRGAP2C the longest
  primary read spans 206,964 bp against a 207,899 bp gene; SRGAP2 260,860 vs 260,886.
- **"The truncation is inherent"** — SUPPORTED, but for the unspliced population only, and for a different
  reason: no intron chain means no evidence linking distant read clusters.

⚠ Two measurement cautions found while establishing this:
- Max-read-span is contaminated by MIS-CHAINS: NPIPB2 shows a read spanning 1270% of its gene, GTF2IRD2B
  680%, ULK4P2 750%. Only spans near 100% (SRGAP2/2C/2D) are evidence of genuine full-length alignment.
- Read presence over an unincluded exon does not prove the reads BELONG to that copy — in SD regions they may
  be homology-shadow multimappers. `bench/soto/truncation_cause.py` establishes the locus is covered
  (65/65 truncated pairs have read-supported exons we omitted), not that the reads are assignable.

**Recommendation: report size agreement stratified.** The spliced number (0.77 median, 66% within 2x) is the
meaningful one for "do our transcript models have the right size"; the unspliced number measures a read-cluster
locus against a gene model, which is a category mismatch rather than a pipeline error.

## 13. Empirical rules for merging unspliced clusters — tested against ground truth (2026-07-28)

§12 left the question: can unspliced read-clusters be merged into one locus without reintroducing
readthrough over-merge? `bench/soto/unspliced_merge_rules.py` tests candidate rules against the Soto members
(POSITIVE = two unspliced copies inside the SAME member, NEGATIVE = inside DIFFERENT members, same
chromosome, gap < 300 kb so the decision is non-trivial). 400 positives / 141 negatives.

**A structural constraint first.** `cluster_unspliced` is single-linkage span-overlap clustering, so if any
read overlapped both clusters they would already BE one cluster. No unspliced read can bridge two clusters —
that witness is unavailable *by construction*. The evidence must come from reads the unspliced path never
used: **spliced** reads with aligned blocks in both clusters.

### Headline (always-merge baseline precision = 0.739)

| rule | precision | recall | F1 |
|---|---:|---:|---:|
| SPLICED-BRIDGE >= 1 | 0.924 | 0.335 | 0.492 |
| GAP-COVERAGE >= 3 | 0.734 | 0.958 | 0.831 |
| GAP <= 50 kb (naive proximity) | 0.938 | 0.605 | 0.736 |
| **bridge >= 1 OR gap <= 50 kb** | **0.925** | **0.682** | **0.786** |

### The gap-matched control, which is what actually decides this

Positives are pairs INSIDE one member (gap <= member length) while negatives are BETWEEN members, so a raw
gap rule partly encodes "member length vs inter-member distance" rather than same-locus evidence. Holding
the gap constant:

| gap bin | n+ / n- | base rate | bridge >= 1 precision | gapcov >= 3 precision |
|---|---|---:|---:|---:|
| 0-10 kb | 88 / 3 | 0.967 | 0.938 (**worse than base**) | 0.986 |
| 10-50 kb | 154 / 13 | 0.922 | 0.967 | 0.922 |
| **50-150 kb** | 143 / 67 | **0.681** | **0.882** | 0.681 |
| 150-300 kb | 15 / 58 | 0.205 | 0.333 (n+ = 15, unreliable) | 0.205 |

Three findings, two of them negative:

1. **GAP-COVERAGE carries literally zero information.** Its precision equals the base rate to three decimals
   in EVERY bin (0.986 / 0.922 / 0.681 / 0.205). In SD regions everything is covered, so "is the intervening
   region covered?" — the intuitive signal — discriminates nothing. Worth recording precisely because it is
   the first thing one would reach for.
2. **The spliced-bridge certificate is real but NARROW.** It only beats the base rate where proximity is
   weak, the 50-150 kb regime (0.681 -> 0.882), and even there recall is 0.21. Below 10 kb it is *worse*
   than the base rate (0.938 vs 0.967) — at short range a spliced bridge is as likely to be a readthrough
   transcript as evidence of one locus.
3. **The naive distance rule is hard to beat**, but partly by construction: Soto members are compact, so
   "same member" and "small gap" are nearly the same statement on this benchmark. That is a property of the
   benchmark, not a general law, and it would fail wherever members are adjacent.

### Recommendation

`bridge >= 1 OR gap <= 50 kb` is the best precision-preserving rule found (0.925 / 0.682), but it still
over-merges 7.5% of pairs, which is worse than the current behaviour of not merging at all when the goal is
partition purity. The defensible use of the certificate is **narrow and additive**: extend merging into the
50-150 kb regime where proximity alone is only 68% precise, and abstain elsewhere — an assign-or-abstain
shape consistent with the rest of the method, rather than a global threshold.

**Not implemented.** The measurement is the deliverable: it establishes that the intuitive signal (gap
coverage) is worthless, that the principled certificate buys accuracy only in one band, and that a naive
distance rule's apparent strength here is partly a benchmark artifact.

## 14. Does the method find the FULL locus of a complete gene? No — and this qualifies the recall headline (2026-07-28)

Member recall counts a member as recovered when **any** predicted copy overlaps it. That is a locus-*overlap*
test, not a locus-*recovery* test. `bench/soto/locus_completeness.py` asks the stronger question: taking the
UNION of every predicted copy overlapping a member (so fragmentation is forgiven), what fraction of that
gene's ANNOTATED EXONIC bases do we cover? Stratified by `gene_biotype`, since "complete gene" is the
question — of the 252 Soto members carrying a named annotation, 126 are protein_coding.

### Complete genes (protein_coding, n = 105 with annotated exons in the member interval)

| | |
|---|---|
| detected at all (>= 1 overlapping copy) | **76/105 = 72%** |
| median exon coverage, over detected | **0.69** |
| **>= 90% of annotated exons covered** | **18/76 = 24% of detected, 17% of all** |
| >= 75% | 31/76 = 41% of detected |
| >= 50% | 50/76 = 66% of detected |
| spliced copies only | median 0.53; >= 90% in 13/76 |

By biotype: transcribed_pseudogene 81% detected / median 0.61; pseudogene 42% detected / median 0.77.

**Isoform-fairness control.** Scoring against the union of all isoforms' exons is a harsh bar — no single
transcript expresses it. Re-scored against the **best-matching single transcript** (median 4 isoforms/gene):
median coverage **0.67**, >= 90% in 36%. So the union bar accounts for only ~12 points at the 90% threshold
and does not explain the result. (This control joins annotation slightly differently and so runs over 90
detected members rather than 76; the distributions are the comparable quantity, not the denominators.)

### The qualitative finding underneath: some "detections" are purely INTRONIC

Several protein_coding members are scored as recovered while covering **zero** annotated exons:

```
SRGAP2     2 copies, exon_cov 0.00   both clusters lie inside ONE ~93kb intron
                                     (13,068bp/95 reads and 6,696bp/3 reads; nearest exon 1.7-9.7kb away)
NPIPB2     1 copy,   exon_cov 0.00   1,356bp cluster, nearest exon 2.7kb away, intronic
NOTCH2     1 copy,   exon_cov 0.03      NOTCH2NLC 0.03    NOTCH2NLA 0.04    LIMS4 0.05
```

These are unspliced clusters sitting in intronic sequence — pre-mRNA or unannotated intronic transcription —
which overlap the gene's *span* and therefore satisfy the overlap-based recall test while contributing
nothing to the gene's exon structure. This is the same population identified in §12 as driving the size
truncation, now shown to be positionally intronic rather than merely short.

### Honest statement

> Member recall (76.2% / 65.5% / 71.0% depending on recipe, §11) measures whether a Soto member's locus is
> **touched** by a predicted copy. It is not a statement that the gene was reconstructed. For complete
> (protein-coding) genes, 72% are touched, but only **17%** have >= 90% of their annotated exons covered even
> when all overlapping copies are unioned, and the median covered fraction is **0.69**. Some members counted
> as recovered are represented solely by intronic unspliced clusters covering none of the gene's exons.

Both numbers should be quoted together, exactly as §1 argued for member recall vs partition completeness.
The gap between them is not a defect to be hidden — it is the honest scope of what an RNA-only method
recovers — but presenting overlap-recall alone overstates it.

## 15. Can intronic clusters be excluded WITHOUT annotation? Tested — no (2026-07-28)

§14 ended by suggesting a copy should overlap an annotated exon before counting toward member recall. That is
legitimate as a **scoring** change — an evaluation may use information the method cannot — but it does not
fix the pipeline, and gating *detection* on a GFF would destroy the annotation-independence that is the whole
point of the approach. So: can the reads alone tell that a cluster is intronic?

The reads do carry a direct observation of splicing: a spliced read's `N` CIGAR gap **is** an intron. Two
annotation-free rules were tested on the 254 unspliced clusters that lie inside a Soto member (label, from
annotation and used ONLY for grading: intronic = 0 annotated exonic bases; base rate 0.374):

| rule (reads only) | best precision | recall |
|---|---:|---:|
| >= k spliced reads whose intron CONTAINS the cluster | 0.482 (k=10) | 0.705 |
| intron-witnesses OUTVOTE exonic-block witnesses among spliced reads (>=95%, >=5 votes) | **0.528** | 0.695 |
| *(base rate — always predict intronic)* | *0.374* | *1.000* |

Both are informative but neither is usable: the best rule reaches 0.53 precision, i.e. **roughly half of the
clusters it would discard are not intronic**. The second rule is the more principled one — it excludes
unspliced reads (which created the cluster, so counting them is circular) and makes it a contest between
mature-transcript evidence for and against the interval being retained — and it gains only ~0.05 over the
cruder version.

**Why it fails is itself the finding.** In segmental duplications a spliced read that skips this interval may
belong to a *paralogous copy* in which the interval genuinely is intronic, while for this copy it is exonic.
Read-level splice evidence is therefore not copy-specific — the same multimapping ambiguity the whole method
is about, reappearing at the level of "is this base transcribed here?".

⚠ One caveat that cuts in the rule's favour and cannot be resolved here: a cluster with zero *annotated*
exonic bases may be a genuine **unannotated** exon rather than intronic pre-mRNA, in which case some of the
counted false positives are correct calls and the annotation label is the thing at fault. Settling that needs
independent evidence (CAGE, poly(A), or cross-tissue reproducibility), not more of this BAM.

### Consequence

There is currently **no reliable annotation-free way to exclude intronic clusters from the catalog.** The
overstatement §14 identified therefore cannot be fixed inside the method's own constraints — only in the
evaluation, where annotation is fair game. This is an argument for reporting overlap-recall and exon-coverage
side by side (§14) rather than for adding a filter: the filter would be wrong half the time, and it would
make the method depend on exactly the annotation the thesis claims not to need.

For the write-up this is a limitation worth stating plainly rather than engineering around: **an RNA-only,
annotation-free method cannot in general distinguish intronic pre-mRNA signal from exonic signal at a
duplicated locus**, and that is part of why exon-level recovery (17% of complete genes at >= 90%) sits so far
below locus-overlap recall (72%).

## 16. In-scope resolution: a certificate that decides a minority, and abstains on the rest (2026-07-28)

Scope for this project is **T2T CHM13v2.0 + long RNA-seq reads only** — no CAGE, poly(A), or cross-tissue
data. §15's suggestion to settle the unannotated-exon question with orthogonal assays is therefore withdrawn.
Within scope, one signal remained untested.

§15 tested only NEGATIVE evidence (reads splicing OVER a cluster) and got 0.53 precision, because in an SD a
read skipping the interval may belong to a paralog. POSITIVE evidence is not symmetric: a junction whose
**acceptor lands on the cluster's start**, or whose **donor lands on its end**, is direct observed proof that
this interval is retained in a mature transcript. Intronic pre-mRNA has no reason to generate junctions at
its own boundaries.

Measured on the same 254 clusters (TOL = 25 bp; exonic base rate 0.626):

| rule (reads only) | precision | recall | n fired |
|---|---:|---:|---:|
| acceptor >= 1 **AND** donor >= 1 (both boundaries) | **1.000** | 0.025 | 4 |
| acceptor >= 3 OR donor >= 3 | 0.882 | 0.094 | 17 |
| acceptor + donor >= 10 | 0.889 | 0.050 | 9 |
| acceptor >= 1 OR donor >= 1 | 0.750 | 0.170 | 36 |

**This is a certificate, not a classifier.** Bilateral junction support is perfect on this benchmark (4/4)
and unilateral support at >= 3 reads is 0.88 — but together they decide only ~10-17% of clusters. The
structural reason mirrors §13: these clusters were built *from unspliced reads*. Had junctions abutted them,
spliced reads would have formed a spliced skeleton there instead, and they would not be in this population at
all. The evidence is scarce for the same reason the cluster exists.

### Where this leaves the question, honestly

Within the project's data, the three available signals are:

| evidence | direction | precision | coverage |
|---|---|---:|---:|
| junction abuts a boundary (§16) | exonic | 0.88 - 1.00 | ~10% |
| spliced reads outvote exonic blocks (§15) | intronic | 0.53 | ~70% |
| annotation overlap (§14) | either | — | out of method scope |

So a **minority can be decided with a read-level certificate, and the majority cannot be decided at all from
T2T + long RNA reads.** That is not a gap to be engineered away — it is the same identifiability argument the
method already makes for copy assignment (K=0, `project_k0_frontier_unresolvable`), reappearing on a new
axis: *is this base transcribed in this copy?* The honest move is the one the method already takes elsewhere
— **assign where a witness exists, abstain otherwise** — and to report exon-coverage alongside
overlap-recall (§14) rather than filtering the catalog on evidence that is wrong half the time.

**Suggested framing for the write-up:** the intronic-cluster population is not a bug but a second
identifiability frontier. Copy assignment has K=0 (copies that reads cannot separate); locus reconstruction
has this (intervals whose transcriptional status reads cannot determine). Both are properties of the data,
both are stated as theorems about what is knowable rather than as pipeline failures, and both are why
exon-level recovery (17% of complete genes at >= 90%) sits below locus-overlap recall (72%).

## 17. How many of Soto's 83 families does RNA mode find? (2026-07-28, verified)

Independently recomputed by a separate agent from scratch before reading the script; all counts reproduced
exactly. Catalog = the corrected per-chrom recipe (§11). Two corrections were applied after verification:

- **`family_id` is only unique WITHIN a per-chrom unit** (each unit numbers from `GWFAM0`), so the
  concatenated catalog merges e.g. chr1:GWFAM0 with chr15:GWFAM0. Counts are computed per unit and unioned.
  This does not change FOUND (59) but raises the isoform-qualified count 41 -> 48.
- **7 of the 83 Soto "families" have exactly ONE member** and can never be found *as a family*. The scorable
  denominator is **76**. (Of those 7: 1 is detected, 6 are not.)

### The answer

| standard | count | of 76 scorable |
|---|---:|---:|
| **FAMILY FOUND** — >= 2 members grouped into ONE predicted family | **59** | **78%** |
| **+ ISOFORM REQUIRED** — >= 1 grouped member carries a spliced (>= 2 exon) copy | **48** | **63%** |
| **COMPLETE** — one predicted family covers ALL the family's members | **16** | **21%** |
| singleton only (1 member covered, no family formable) | 6 | 8% |
| missed entirely (0 members covered) | 11 | 14% |

### The isoform requirement matters, and removes the flagship families

Requiring a real isoform — not merely an overlapping copy — drops 11 families, and they are exactly the
ones a reviewer will ask about, because they are represented **only by unspliced clusters** (§14: those can
be entirely intronic):

```
ID_400  NOTCH2, NOTCH2NLA, NOTCH2NLB, NOTCH2NLC     ID_462  SRGAP2B, SRGAP2C, SRGAP2D
ID_395  RGPD1-4                                     ID_226  H2BC18, H2BP1, H3-2, H3C13, H3P4
ID_212  AC244669.1, AC245100.8, LINC00869           ID_448  SEC22B, SEC22B3P
ID_22   AC006453.1, AC027612.2, AL356585.3, BAGE2   ID_24   LSP1P4, LSP1P5
```

NOTCH2NL and SRGAP2 are the canonical human-specific expansions. Reporting them as "found" on the strength of
an intronic unspliced cluster would not survive scrutiny; **63% is the number to quote**, not 78%.

### Verified caveats on the surrounding numbers

An independent audit upheld these, and they qualify figures used elsewhere in this document:

1. **The exon-recovery denominator (§14, §17) is gene-SYMBOL string equality.** The Soto BED uses
   Ensembl-style symbols (AL669831.1, MST1L) while CHM13 RefSeq annotates the same loci as LOC100288069,
   LOC124905698. 124 of 362 members overlap genuine annotated exons but are silently excluded, leaving 200
   judgeable rather than 324. Substituting interval overlap moves exon-formable from 23/62 = 37% to
   36/77 = 47% (definitive) and 26/62 = 42% to 39/77 = 51% (per-chrom). Interval overlap is itself an upper
   bound (it would credit an unrelated overlapping gene's exons); the correct repair is alias reconciliation.
2. **The exon COMPLETE column was inflated** by using judgeable members as the denominator: 8 of 11
   per-chrom exon-COMPLETE families are complete only over the symbol-matched subset. With a fair
   denominator COMPLETE falls 11 -> 3. FOUND is unaffected.
3. **"1107 copies" double-counts 140 loci** (967 distinct) — the per-chrom and xchrom passes both report
   copies in the overlapping regions and the combine does not dedup. `definitive`'s 863 are all distinct, so
   "863 vs 1107" is not apples-to-apples.
4. **definitive vs per-chrom is confounded by binary date**: definitive was built 07-23, the per-chrom units
   07-28, with source changes in between. The 69% -> 71% delta mixes the recipe fix with 5 days of code
   changes and should not be attributed solely to the recipe.

## 18. Applying the RNA lessons to the DNA side (2026-07-28)

### What transfers, and what does not

| lesson from the RNA work | applies to DNA? |
|---|---|
| scorable denominator is 76, not 83 (7 single-member families) | **yes** — same benchmark |
| strict grouping: >= 2 members in ONE predicted family, not merely "touched" | **yes** |
| **purity guard** — do not credit a prediction that also absorbs OTHER Soto families | **yes, and it is decisive** |
| isoform requirement (>= 1 spliced copy) | **no** — no splicing in DNA |
| exon-coverage / truncation analysis (§12, §14) | **no** — `--from-genome` seeds nodes FROM the Soto windows, so span and coverage are circular by construction (§8) |
| member-detection recall | **no** — circular for the same reason (~100% by construction) |

The grouping IS fair: which seeded windows the method joins into a family is its own output and was never
given to it. So family-level partition verdicts are directly comparable between the two modes.

### Like-for-like, both modes, same definitions (denominator 76)

| | RNA (per-chrom corrected) | DNA (`--from-genome`) |
|---|---:|---:|
| predicted families absorbing >= 2 Soto families (over-merged) | 21 | 18 |
| COMPLETE — one predicted family covers all members | 16 (21%) | **66 (87%)** |
| **COMPLETE — and that family is PURE** | **15 (20%)** | **35 (46%)** |
| FOUND — >= 2 members in one predicted family | 59 (78%) | **76 (100%)** |
| **FOUND — and that family is PURE** | **51 (67%)** | **43 (57%)** |
| (RNA only) FOUND + >= 1 member carries an isoform | 48 (63%) | n/a |

### What this says

**The purity guard is to DNA what the isoform guard is to RNA** — the correction that stops a headline being
carried by artifacts. It barely touches RNA's COMPLETE (21% -> 20%) but roughly halves DNA's (87% -> 46%),
and DNA's perfect 100% FOUND becomes 57%. That is the §1/§7 finding — DNA over-merges (homogeneity 74%),
RNA is purer (91%) — arriving again through an independent route.

The honest two-line summary of the two modes:

> **DNA reconstructs whole families better** (46% complete-and-pure vs RNA's 20%), because it has the
> intronic and flanking sequence that makes homology edges easy to form. **RNA forms pure families more
> often** (67% found-and-pure vs DNA's 57%), because splicing discards exactly the duplicon sequence that
> spuriously bridges unrelated families at the DNA level. Neither is uniformly better, and both are far below
> their unguarded numbers (DNA 100%/87%, RNA 78%/21%).

⚠ Do not present DNA's 361/362 = 99.7% member recovery beside RNA's member recall as if they measure the same
thing. DNA's is ~100% by construction (windows seeded from the truth); RNA's is a detection result. The
family-partition rows above are the only rows where the two modes are genuinely comparable.
