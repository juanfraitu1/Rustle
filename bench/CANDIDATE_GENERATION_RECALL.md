# Candidate-generation recall: are the 7 missed multi-copy genes a pre-filter bottleneck?

**Question.** `bench/GRAPH_TO_GRAPH_FAMILY_ALN.md` found the family-definition recall gap (48/57
diploid-oracle multi-copy genes recovered) is an **edge-creation** problem: 7 of the 9 missed genes
(**LOC101141440, LOC115930538, LOC115934629, LOC129534585, TNK2, UBE2Q2P16, ZNF425**) were reported to
have *zero* candidate edges. The shipped candidate-generation (`bench/denovo_families.py`) proposes a
pair only if the two de-novo transcripts share ≥ `K_SHARE=6` exact canonical **18-mers**; POA
`core_recip ≥ 0.13` then decides membership. **Is the miss a pre-filter that is too coarse (recoverable
by shorter k / minimizers), or is it expression / divergence-limited (needs more RNA / DNA)?**

Each missed gene is classified into one of four causes — the task's three plus **REP-COLLAPSE/
PROJECTION**, the honest fourth surfaced by the data:

1. **PRE-FILTER-BOTTLENECK** *(the only loosening-recoverable class)* — an expressed paralog pair passes
   POA scoring but the exact-18-mer screen did **not** propose it; a shorter-k / minimizer screen
   proposes it and it passes → recoverable from RNA.
2. **EXPRESSION-LIMIT** — the family has no distinct expressed paralog copy → no edge possible from this
   IsoSeq regardless of the pre-filter.
3. **DIVERGENCE-LIMIT** — a candidate pair exists (proposed at k18 or by loosening) but the family's
   *cleanest* cross-copy pair **fails** POA scoring (too divergent) → DNA-bound.
4. **REP-COLLAPSE/PROJECTION** — the family **does** contain a POA-passing, 18-mer-proposed cross-copy
   member pair, so candidate-generation **and** scoring both already succeed for the family; the gene's
   annotation region is lost **downstream** (intron-junction rep-collapse merges the gene's copy into a
   readthrough, `gene_of` max-overlap labels the linkable locus as a neighbour gene). Not a
   candidate-generation-k problem; loosening cannot help (the edge is already proposed at k18).

Artifacts (deterministic, `PYTHONHASHSEED=0`, `minimap2 -t1`, Bio.Align global NW):
`bench/candidate_generation_recall.py`, `.tsv` (one row per cross-copy candidate edge, md5
`beee5f94…`), `.json` (md5 `c3eb0778…`, byte-identical on `--json` re-run).

---

## Method fix: score the ACTUAL family members (not a span-derived rep)

An earlier revision re-derived each copy's representative as the **max-n_reads de-novo locus overlapping
the copy's genomic span**. On rep-collapsed loci that picks the *wrong* transcript — a huge readthrough
that is **not a family member**, or a member of a **different** family — and mis-scores the pair. Concrete
failure it produced:

- **LOC101141440**: it scored the 129 kb readthrough `DN_NC_073232.2_99594354_7` (in **no** validated
  family) against `DN_NC_073232.2_92169976_5` (a member of family **154**, not the gene's family 21) →
  core_recip **0.107** (243 fragmented POA blocks — the readthrough signature) → mislabelled
  **DIVERGENCE-LIMIT**. The **real** family-21 members score core_recip **0.926** (2 blocks, 99.75% core
  identity).

The fix (`cluster_members` / `best_cluster_pair`): cluster the **actual `validated_families` member loci**
into distinct genomic copies (gap > 100 kb), keep the top-8 members per copy by n_reads, and score **every
cross-copy member pair** with the shipped machinery. Classification is driven by the family's strongest
edge (`fam_best`); the gene's-own-copy edge is reported separately (`gene_to_paralog_edge` column,
`gene_own_copy_links` in JSON) so a family recoverable from its *other* copies while the gene's *own*
transcript is a divergent/short-span paralog is flagged honestly. Every `validated_families` member is an
expressed de-novo transcript by construction, so **every paralog copy IS expressed** (EXPRESSION-LIMIT is
structurally 0 here).

---

## Per-gene diagnosis

| gene | family | copies | family-best cross-copy pair (core_recip / aln_frac / shared-18mer) | gene's own copy links? | cause |
|---|---|---|---|---|---|
| LOC101141440 | 21 | 2 | 0.926 / 0.991 / **2019** — PASS, 18-mer-proposed | yes (0.926) | **REP-COLLAPSE/PROJECTION** |
| LOC115930538 | 34 | 7 | 0.997 / 0.997 / **5311** — PASS, 18-mer-proposed | **no** (gene copy 0.12) | **REP-COLLAPSE/PROJECTION** |
| LOC115934629 | 32 | 2 | 0.850 / 0.850 / **3725** — PASS, 18-mer-proposed | yes (0.850) | **REP-COLLAPSE/PROJECTION** |
| LOC129534585 | 7 | 4 | 1.000 / 1.000 / **1401** — PASS, 18-mer-proposed | yes (1.000) | **REP-COLLAPSE/PROJECTION** |
| ZNF425 | 65 | 4 | 0.999 / 1.000 / **1237** — PASS, 18-mer-proposed | **no** (gene copy 0.029) | **REP-COLLAPSE/PROJECTION** |
| TNK2 | 136 | 2 | 0.053 / 0.0 / 90 — FAIL (2-member family) | no (0.053) | **DIVERGENCE-LIMIT** |
| UBE2Q2P16 | 175 | 2 | 0.023 / 0.0 / 0 — FAIL (2-member family) | no (0.023) | **DIVERGENCE-LIMIT** |

**Tally: 0 PRE-FILTER-BOTTLENECK / 0 EXPRESSION-LIMIT / 2 DIVERGENCE-LIMIT / 5 REP-COLLAPSE/PROJECTION.**

- **5/7 are REP-COLLAPSE/PROJECTION.** Each family contains a cross-copy member pair that passes POA
  (core_recip 0.85–1.0) **and** shares thousands of exact 18-mers (2019/5311/3725/1401/1237) — so the
  shipped 18-mer screen **already proposes** it and POA **already accepts** it. Candidate-generation is
  *not* the bottleneck. The gene's annotation region is missed downstream. Two sub-cases:
  - *gene copy links directly* (LOC101141440, LOC115934629, LOC129534585): the gene's own transcript is
    itself a passing, 18-mer-proposed edge (0.93 / 0.85 / 1.0) — pure rep-collapse / `gene_of` labelling.
  - *family recovered from its OTHER copies* (LOC115930538, ZNF425): the family's near-identical copies
    form the edges, while the gene's **own** transcript is a divergent/short-span paralog (best gene edge
    0.12 / 0.029, below the 0.13 gate). Even here loosening is irrelevant: LOC115930538's gene copy is
    already 18-mer-proposed (279 shared) and POA rejects it; ZNF425's gene copy fails POA (0.029) and is
    not proposed even at k11.
- **2/7 are genuine DIVERGENCE-LIMIT** (TNK2, UBE2Q2P16): 2-member families whose single possible
  cross-copy RNA pair fails the POA gate (core_recip 0.05 / 0.02, both < the 0.13 shipped line and the
  0.19/0.24 task line — threshold-robust). TNK2 is already 18-mer-proposed (90 shared) and POA correctly
  rejects; UBE2Q2P16 has 0 shared 18-mers, loosening to k=11 proposes it (41 shared) but POA still fails
  (0.023). DNA-bound.

**Correction vs the first-pass narrative.** The earlier "5/7 DNA-bound divergence" was wrong: 3 of those
5 (LOC101141440, LOC115930538, ZNF425) were mis-scored by the span-derived rep against non-member /
wrong-family readthroughs. Scoring the real members flips them out of divergence — their families contain
99.7–99.9%-identical expressed members. Honest split: **~2/7 divergence, ~5/7 downstream
projection/rep-collapse.** The bottom line is unchanged: none is a candidate-generation-k problem.

---

## Would loosening (shorter k / minimizers) recover the pre-filter-limited genes?

**Recoverable by loosening = 0/7.** PRE-FILTER-BOTTLENECK is the only loosening-recoverable class and it
is empty. For **all** 5 REP-COLLAPSE genes the family's linking edge is already proposed at k18 (thousands
of shared 18-mers); for the 2 DIVERGENCE genes the pair is already proposed (TNK2: 90) or fails POA even
when loosening does propose it (UBE2Q2P16 at k11). No `(shipped_proposes=0, loosened_proposes=1,
passes_scoring=1)` row exists in the TSV.

### Genome-wide candidate blow-up + POA compute (`--cost`, 150-locus sample, seed 0, full 94,192-transcript scan)

| k | fan-out multiplier vs k18 | median partners/locus | genome-wide candidate pairs (raw) | extra POA vs k18 |
|---|---|---|---|---|
| 18 | 1.0× | 30 | ~91.7 M | — |
| 15 | 1.6× | 120 | ~143.0 M | +51.3 M (~133 CPU-days) |
| 13 | 3.3× | 2,694 | ~299.0 M | **+207.2 M (~538 CPU-days)** |
| 11 | 36.9× | 75,545 | ~3.38 B | +3.29 B (~8,546 CPU-days) |

At `0.2243 s`/POA-pair single-thread. At k=11 the median locus shares ≥6 k-mers with ~75k of 94k
transcripts (near-total saturation — matches `denovo_families.py`'s own note that k=13/16 are "too
saturated"). All of this extra compute buys **0** recovered genes.

### FP / precision impact (`--fpprobe`, 80 random targets scanned against all 94,192 transcripts, seed 0)

Loosening changes **only k**; the downstream POA gate (`core_recip ≥ 0.13`) still decides every edge. The
question is whether the extra k13-only candidate pairs survive the gate and become over-merges.

- pools: base k18 = 157,618; loosened k13 = 516,059; **extra (k13-only) = 358,441**.
- **extra (k13-only) POA survival at the shipped gate = 0/300 = 0.0 %** (max core_recip observed
  **0.0706**, far below 0.13) — *below* the k18 baseline survival of **1.33 %** (4/300). The loosened-only
  pairs are the more-divergent tail; POA rejects all sampled.
- Expected extra over-merge edges over the ~207 M raw k13-extra pairs: point estimate **0**; rule-of-3
  95%-CI upper bound ≤ **~2.07 M**. The FP over-merge set does **not** explode — POA holds precision.

So loosening is **recall-SAFE** (precision holds) but **recall-USELESS** (0 genes recovered).

---

## Verdict

**The family recall gap is NOT a fixable candidate-generation bottleneck. Loosening the pre-filter is not
worth adopting.**

- **Recall gain from loosening = 0/7.** The recall gap is not candidate-generation-k-bound — every family
  edge that matters is already proposed at k18.
- **~2/7 are genuinely DNA-bound divergence** (TNK2, UBE2Q2P16): the cleanest cross-copy RNA pair fails
  POA; recovery needs the paralog copies to diverge less (DNA), not a shorter k.
- **~5/7 need a different downstream fix, not a shorter k** — the family's copies are near-identical,
  expressed, and already 18-mer-proposed + POA-passing; the gene's annotation region is lost to
  intron-junction rep-collapse and `gene_of` max-overlap projection (and, for 2 of the 5, the gene's own
  transcript is itself a divergent/short-span paralog while the family is recovered from its other copies).
- **Cost of loosening:** 3.3× candidate/POA load at k13 (~+207 M pairs, ~538 CPU-days), 36.9× at k11
  (~+3.3 B pairs, ~8,546 CPU-days, near-saturation). Precision is safe (0/300 k13-only pairs survive the
  POA gate, max core_recip 0.0706 ≪ 0.13; expected extra over-merges 0, 95%-CI ≤ ~2.07 M) — which only
  confirms the extra compute is pure waste.

**Return values.** Diagnosis counts = `{PRE-FILTER-BOTTLENECK: 0, EXPRESSION-LIMIT: 0, DIVERGENCE-LIMIT:
2, REP-COLLAPSE/PROJECTION: 5}`. Recall-recoverable-by-loosening = **0/7**. Cost = 3.3× (k13) / 36.9×
(k11) candidate fan-out, ~538 / ~8,546 extra CPU-days, FP-safe (0/300 loosened-only survive the gate).
Verdict = **loosening not worth it; the recall gap is ~2/7 divergence (DNA-bound) + ~5/7 downstream
rep-collapse/projection, 0/7 candidate-generation-k-bound.**
