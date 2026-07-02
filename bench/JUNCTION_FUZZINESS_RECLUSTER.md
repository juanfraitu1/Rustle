# Junction fuzziness (FLAIR-style wobble tolerance) vs the shipped exact locus definition

**Question.** The shipped locus definition (`bench/denovo_families.py:149-183`) clusters transcripts
into gene-loci by **union-find on identical `(chrom, donor, acceptor)` intron triples** — no junction
fuzziness, no FLAIR-style snapping. minimap2 spliced junctions wobble ±a few bp (worse at low
coverage / divergent paralogs), so a low-coverage or divergent paralog copy's reads could in
principle **fragment** under exact matching into partial loci too short to reach the family stage —
masquerading as an expression limit. Does adding a **junction wobble tolerance** (snap donor/acceptor
sites to a consensus within ±N bp, then union-find; sweep N=0,2,5,10,20) recover any real missed
paralog locus, and at what copy-collapse cost?

**Answer (one line).** Junction fuzziness is a **no-op**: it recovers **0** missed loci and costs
**0** distinct-copy collapses at every N. **best_N = 0 (keep the shipped exact union-find).** The
family recall gap is **NOT** a fixable clustering-fuzziness artifact — it lives downstream.

Reproducer (deterministic, `PYTHONHASHSEED=0`):
- `bench/junction_fuzziness_recluster.py`
- `bench/junction_fuzziness_recluster.json` (sha256 `ea6fb1b8…b95b37`)
- `bench/junction_fuzziness_recluster.tsv` (sha256 `7d606e37…e33e8`)

Two `PYTHONHASHSEED=0` runs produce byte-identical JSON and TSV.

---

## Method

Operates at the **read level**: each read's CIGAR `N` ops → `(chrom, donor, acceptor)` junctions.
`N=0` reproduces the shipped exact identical-junction union-find **rule** (`denovo_families.py:164`,
key `(chrom, donor, acceptor)`; single-exon reads stay singletons). `N>0` = FLAIR-style
support-anchored greedy consensus snapping: the highest-support unassigned donor (resp. acceptor)
site absorbs every unassigned same-chromosome site within ±N, then union-find on the snapped keys.
RNA-only; the DNA oracle supplies target/collapse **truth only**.

**Caveat (disclosed).** The baseline is a **read-level** reproduction of the union-find *rule*, one
level below the shipped assembled-transcript collapse. It therefore probes only the fragmentation
introduced by junction **wobble** (the read→reference spliced alignment — exactly where wobble
lives), not any fragmentation introduced by the assembly/skeleton step itself. This is the right
level to isolate the fuzziness lever; it does not exercise assembly-stage fragmentation.

**Two fixes applied to the script during this review:**
1. **Chromosome is now part of the junction key** (`snap_map`/`cluster` are chrom-aware; snapping is
   strictly within-chromosome). Previously the key dropped chrom — inert on this data (all recall
   regions and all selected collapse families are single-chromosome) but a latent bug that could
   spuriously union two loci by a cross-chromosome coordinate coincidence when reads are pooled.
2. **The collapse test now runs TWO passes** (see below). The original single pass pre-merged any
   copies <5 kb apart into one block (`BLOCK_GAP=5000`) and so **structurally hid the
   adjacent-tandem configuration — the only real collapse risk.** Its headline "tightest inter-copy
   gap = 9021 bp" is a **floor imposed by that pre-merge**, not an empirical minimum. An
   adjacent-tandem stress pass (no 5 kb pre-merge) was added.

---

## (1) RECALL — do wobble-fragmented paralog loci reappear under fuzzy?

**23 recall targets** = the 7 missed diploid-oracle genes' own annotation loci + each of their
expressed paralog copy loci (`bench/candidate_generation_recall.tsv`). **20 / 23 already form an
admitted multi-read locus under EXACT (N=0).** The 3 non-admitted are:

| target | reads | spliced | why unadmitted at exact |
|---|---|---|---|
| `LOC115930538` own (NC_073240.2:21,123,901-21,134,037) | 1 | 1 | single read — nothing to cluster |
| `LOC129534585` own (NC_073224.2:15,805,110-15,833,058) | 1 | 0 | single **single-exon** read — no junction to snap |
| `LOC115934629` own (NC_073228.2:45,754,758-45,835,592) | 18 | 10 | 10 spliced reads, 10 **distinct** junction chains |

**Named missed loci recovered under fuzzy: NONE** (`recovered_loci = []`, `genes_recovered = []`).
**copies_rescued (unadmitted → admitted) = 0 at every N.** No clean locus emerges at N>0 that did not
exist at N=0.

The one genuine spliced candidate, **`LOC115934629` own**, stays `multi=0 / 10 loci` even at N=20.
Its 10 spliced reads carry 10 distinct junction chains — e.g. three reads share donor `45782046`
but their acceptors are `45802605 / 45803510 / 45833543` (>900 bp apart). These are **real isoform
diversity, not ±20 bp alignment wobble**, so no amount of snapping unites them. Not a fuzzy-fixable
fragment.

**Reads / fragments recovered per N (summed over all 23 targets):**

| N | spliced orphans recovered | fragment-merges | newborn 2-read loci | copies_rescued | total loci |
|---|---|---|---|---|---|
| 0 (exact) | — | — | — | 0 | 350 |
| 2 | 2 | 1 | 1 | 0 | 349 |
| 5 | 3 | 2 | 1 | 0 | 348 |
| 10 | 2 | 2 | 0 | 0 | 348 |
| 20 | **3** | **3** | 0 | 0 | 347 |

Negligible and cosmetic: at N=20, 3 reads and 3 merges across **all** targets; the lone "newborn"
2-read locus (N=2/5) is two orphaned spliced reads of the `LOC115930538` paralog coalescing — a copy
**already admitted** at exact. Total loci barely move (350 → 347).

**Mechanism (the decisive finding).** Under exact matching the spliced reads of each paralog copy
*already coalesce* into 1–few loci; there is essentially **no wobble-fragmentation to recover**. The
large per-region locus counts are dominated by **single-exon reads**, which have no junction and are
structurally immune to snapping — e.g. the `LOC101141440` paralog: 180 reads → 64 loci, but 61 are
single-exon singletons; the 119 spliced reads already form essentially one locus (unchanged 64→64 at
every N). `LOC129534585`, `LOC115934629` paralog, `TNK2`, `ZNF425`, `UBE2Q2P16` copies each hold
their locus count across all N.

---

## (2) COPY-COLLAPSE — do distinct diploid-oracle copies merge under fuzzy?

Collapse is tested **twice** over `validated_families.tsv` ∩ `diploid_cn_oracle.tsv` `multi_copy`
families (DNA-confirmed distinct copies):

**Pass A — distant-copy** (top-8 families by tightest gap, `BLOCK_GAP=5000` pre-merge). Tightest
surviving gap **9021 bp** (a *floor* from the 5 kb pre-merge). Fuzzy-induced collapse **= 0 at every
N**. The only collapse component (**FAM26**, NC_073229.2) is present already at **N=0/exact**
(`exact_present_collapse_pairs=1`) — an identical-junction share the shipped pipeline also merges,
constant across all N, **not** snapping-caused.

**Pass B — adjacent-tandem stress** (no 5 kb pre-merge; every non-overlapping DN-locus group is a
candidate copy-block; the 6 `multi_copy` families with ≥2 same-chrom blocks <5 kb apart). This is the
configuration the user flagged as the **only** real collapse risk and that Pass A hides. Each
fine-block is a *superset-of-a-true-copy* candidate, so any fuzzy cross-block merge here is an **upper
bound** on true distinct-copy collapse.

| family | tightest gap | reads (spliced) | fuzzy-induced collapse (N=2..20) | note |
|---|---|---|---|---|
| FAM102 | **2 bp** | 3217 (140) | **0** | 1 within-block frag-merge at N≥5 (not cross-copy) |
| FAM177 | 221 bp | 144 (99) | **0** | 1 collapse **present at N=0/exact**, constant across N |
| FAM42 | 644 bp | 243 (242) | **0** | 1 collapse **present at N=0/exact**, constant across N |
| FAM1 | 959 bp | 1386 (1222) | **0** | 1 within-block frag-merge at N≥2 (not cross-copy) |
| FAM161 | 1436 bp | 52 (24) | **0** | — |
| FAM176 | 2917 bp | 8 (8) | **0** | — |

**Fuzzy-induced distinct-copy collapse = 0 at every N in BOTH passes**, including the 2 bp-adjacent
FAM102 (3217 reads). The two tandem collapse components (FAM177, FAM42) are **exact-present (N=0)**
identical-junction shares — shipped-exact behaviour, not a cost of fuzziness — and are constant across
N. **Tandem vs genuine split: 0 genuine fuzzy-caused merges; 0 adjacent-tandem fuzzy-caused merges;**
the only collapses (FAM26 distant, FAM177/FAM42 tandem) are all N=0/exact.

**Why it is structurally 0.** Snapping merges junctions only when their **absolute** coordinates lie
within ±N. Tandem copies are 7–70 kb long, so a copy's corresponding junctions sit ~one copy-length
(thousands of bp) apart — far outside ±20 bp — even when the two copies' *loci* abut at 2 bp. Distant
copies (median spacing ~554 kb) are orders of magnitude further. Snapping **cannot** collapse distinct
copies on this substrate.

---

## (3) NET trade-off + best wobble tolerance

`net_score = copies_rescued − fuzzy_induced_distinct_copy_loss = 0 − 0 = 0` at **every** N
(2, 5, 10, 20). No recall gain, no collapse cost — a **wash**.

**best_N = 0** (the shipped exact union-find). Adding a wobble tolerance buys nothing and costs
nothing.

---

## Honest verdict

**The family recall gap is NOT a fixable junction-fuzziness artifact — keep the exact
`(chrom, donor, acceptor)` union-find.** Every expressed paralog copy of the 7 missed oracle genes is
**already admitted** as a locus under exact matching. The sole unadmitted spliced target
(`LOC115934629` own) has genuinely distinct junction chains, not ±20 bp wobble, so fuzziness does not
rescue it. And distinct DNA-confirmed copies — even 2 bp-adjacent tandems — sit far too many bp apart
in junction coordinates for ±N snapping to collapse them. Fuzzy clustering here **trades nothing for
nothing**: it neither fixes recall nor costs copy resolution.

**Where the recall gap actually lives.** Two places, both outside the locus-clustering step:
1. **Downstream family edge-creation / exact-k-mer prefilter** — per
   `bench/candidate_generation_recall.tsv`, the missed genes' loci exist but fail the family-linkage
   stage (low shared-18-mer / minimizer share against their paralog).
2. **Single-exon read degradation** — a large fraction of paralog-copy reads are single-exon (no
   junction to snap or share), so they are structurally invisible to *any* junction-based clustering,
   fuzzy or exact.

**Connection to the other stages.** The shipped exact-junction union-find already absorbs the
degradation / alt-TSS / alt-TES variation the fuzziness lever was hypothesized to fix: reads that
differ only in their 5′/3′ ends but **share any internal intron** are unioned into one locus
regardless of terminal-exon truncation (that is precisely why the paralog copies coalesce under
exact). Junction fuzziness is the one remaining untested lever on top of that — and this experiment
shows it is inert. The exact union-find is already near-optimal at *admitting the spliced loci*; the
recall problem is a family-edge / prefilter problem and a single-exon-read problem, not a
junction-wobble problem.

---

### Result summary
- **Missed loci recovered under fuzzy:** NONE (0 genes, 0 copies_rescued at every N).
- **Reads/fragments recovered:** ≤3 reads / ≤3 merges genome-wide at N=20 (cosmetic, all in
  already-admitted copies).
- **Copy-collapse cost:** 0 fuzzy-induced at every N, in both the distant-copy pass and the
  adjacent-tandem stress pass (tightest tandem 2 bp, FAM102, 3217 reads); the only collapses (FAM26,
  FAM177, FAM42) are N=0/exact, not fuzzy-caused.
- **Best wobble tolerance:** **N = 0** (keep exact).
- **Verdict:** WASH / not a clustering-fuzziness artifact — keep the shipped exact union-find.
