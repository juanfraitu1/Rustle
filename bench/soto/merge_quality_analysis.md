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
