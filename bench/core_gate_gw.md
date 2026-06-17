# Genome-wide contiguous-core family-gate measurement
**Site measured:** rustle's de-novo cross-copy singleton-exon merges (`merge_singletons_by_sequence`) — the gate's actual call site.
**Threshold:** T = 0.13 on contiguous-core coverage (`core_cov`).
**Gate status:** DEFAULT-OFF. This is a *measurement* of what the gate would do, not a behavior change.

## Headline
- **25 chromosomes** parsed (22 carry cross-copy merge pairs; 3 are empty).
- **1876** cross-copy merge pairs total.
- **155 (8.3%) DROP** at core_cov < 0.13.
  - 145 have core_cov < 0.05 (deep lobe).
  - 10 sit in the valley [0.05, 0.13).
- **1721 (91.7%) KEEP** at core_cov ≥ 0.13.
- **616** pairs have core_cov ≥ 0.95 (near-identical true copies).
- `NC_073235.2` = **5 / 67** drops (matches the prior pipeline test).

## Bimodality — why 0.13 is a principled threshold
The core_cov distribution is **bimodal**: a heavy mass near 1.0 (near-identical true copies; 616 pairs ≥ 0.95) and a second lobe down near 0.0 (145 pairs < 0.05). Between them is a sparse **valley** spanning roughly [0.05, 0.13) holding only 10 pairs. T = 0.13 is placed inside that valley, so the decision boundary does not slice through a dense region — small threshold perturbations move very few pairs. This is exactly the kind of separated, non-arbitrary cut a gate wants: the data, not the knob, picks the split.

## Why jaccard alone fails and core fixes it (the danger zone)
At the same numeric cut, **jaccard alone keeps all 1876 pairs** (every pair has jaccard ≥ 0.13) — i.e. a jaccard-only gate at 0.13 would never drop anything here. Yet **106 pairs have high jaccard (≥ 0.50) but low contiguous-core (< 0.13)**, and 64 of those have jaccard ≥ 0.80. These are **domain-sharers**: two copies that share a large fraction of *k*-mers/segments (high jaccard) but whose shared material is *not contiguous core* — it is a repeated domain, a shared exon cassette, or scattered homology rather than a single colinear copy body. Jaccard cannot tell a tandem-repeated domain from a real copy; contiguous-core coverage can, because it demands the overlap form one run.

Worst offenders (jaccard, core_cov) — jaccard says "merge", core says "don't":

| jaccard | core_cov |
|---|---|
| 0.997 | 0.004 |
| 0.995 | 0.004 |
| 0.995 | 0.004 |
| 0.984 | 0.004 |
| 0.971 | 0.008 |
| 0.971 | 0.008 |

## Per-chromosome
| chrom | pairs | drop (core<0.13) | drop %% |
|---|---|---|---|
| NC_073228.2 | 135 | 54 | 40.0% |
| NC_073224.2 | 334 | 24 | 7.2% |
| NC_073246.2 | 28 | 17 | 60.7% |
| NC_073244.2 | 69 | 15 | 21.7% |
| NC_073238.2 | 135 | 14 | 10.4% |
| NC_073240.2 | 219 | 6 | 2.7% |
| NC_073235.2 | 67 | 5 | 7.5% |
| NC_073236.2 | 86 | 3 | 3.5% |
| NC_073229.2 | 9 | 3 | 33.3% |
| NC_073248.2 | 395 | 2 | 0.5% |
| NC_073230.2 | 85 | 2 | 2.4% |
| NC_073247.2 | 42 | 2 | 4.8% |
| NC_073243.2 | 9 | 2 | 22.2% |
| NC_073242.2 | 7 | 2 | 28.6% |
| NC_086017.1 | 6 | 2 | 33.3% |
| NC_073237.2 | 48 | 1 | 2.1% |
| NC_086018.1 | 11 | 1 | 9.1% |
| NC_073241.2 | 114 | 0 | 0.0% |
| NC_073233.2 | 67 | 0 | 0.0% |
| NC_073227.2 | 6 | 0 | 0.0% |
| NC_073231.2 | 2 | 0 | 0.0% |
| NC_073234.2 | 2 | 0 | 0.0% |
| NC_073232.2 | 0 | 0 | 0.0% |
| NC_073239.2 | 0 | 0 | 0.0% |
| NC_073245.2 | 0 | 0 | 0.0% |

## Honest scope & caveats
- **What this measures.** rustle's *de-novo cross-copy singleton-exon merges* (`merge_singletons_by_sequence`). That is the gate's literal call site, so these counts are the population the gate would actually act on.
- **Per-chrom / within-chrom families only.** The OOM-safe genome-wide protocol runs one chromosome at a time, so only *within-chromosome* families are grouped. Cross-chromosome paralogs (e.g. RABL2A/RABL2B) are not co-grouped here — but those are the coordinate-separated case this gate does **not** target anyway; the gate is about distinguishing true copies from domain-sharers inside a candidate family, not about linking copies across chromosomes.
- **`would_gate` trace field is vacuous — ignored.** Every trace line reads `would_gate=false` because it was computed against the active threshold 0.0 (gate-OFF during the trace run). We do **not** trust that field; every drop/keep number above is recomputed directly from `core_cov` against T = 0.13.
- **Gate is default-off.** Nothing here changes rustle's emitted output; this is a measurement of the gate's *would-be* effect.

## Reproduce
```
/home/juanfra/miniforge3/bin/python bench/core_gate_gw.py
```
Inputs: `/tmp/gw/cg_NC_*.trace` (deterministic parse).

## Genome-wide verdict

**What is solidly established (independently re-verified).** Reparsing all 25
`/tmp/gw/cg_NC_*.trace` files (1876 `[CORE_TRACE]` lines) reproduces every
headline exactly: **1876 pairs, 155 (8.3%) drop at core<0.13** (145 deep-lobe
<0.05 + 10 in the [0.05,0.13) valley), 616 near-identical (≥0.95),
**NC_073235.2 = 5/67**. The `core_cov` distribution is genuinely bimodal and
T=0.13 sits in a sparse valley (only 10 of 1876 pairs in [0.05,0.13)), so the
cut is data-placed, not knob-tuned. Core does separate domain-sharers from true
copies that jaccard cannot: **106 pairs (64 at jaccard≥0.8) have high jaccard
but core<0.13** — these would be merged by jaccard, dropped by core. The
`would_gate` trace field is confirmed vacuous (all 1876 = `false`, computed
against gate-OFF 0.0) and is correctly ignored; all decisions are recomputed
from `core_cov`.

**Downstream effect is benign-cleaning, NOT recall loss — confirmed by my own
coordinate-keyed diff.** On the worst-case contig NC_073228.2 (gate fires on
exactly 54/135 pairs), an independent re-keying of every transcript by
(strand, exon-chain) finds **0 chains lost, 0 chains gained, 4989=4989 in both**
(genuine identity, not cancellation). The only output change is **2 attribute
edits at one locus (RSTL.262)**: two guide-backed transcripts (guide:STRG.440.1
and .5) shed a spurious `copy_status "novel"` + `rescue_class
"strand_pure_minority"` tag while keeping identical exon coordinates and
identical cov/FPKM/TPM (verified byte-equal). The transcripts stay emitted.
gffcompare is unchanged both ways: **FSM 1733=1733, intron-chains 1716=1716,
loci 945=945, Sn/Pr 23.3/34.7=23.3/34.7**. So the answer to the downstream
question is the same as NC_073235.2: **benign-cleaning, no real isoform lost** —
the gate removes a false novel-copy attribution and nothing else.

**Honest deflation — do not overclaim.**
- **The downstream payoff is essentially negligible in magnitude.** Even on the
  worst-case contig with 54 gated merges (11x NC_073235.2's 5), the entire
  output effect is **2 attribute tags on 1 of 1497 genes**. Drop count does not
  scale to output churn — the gated cross-copy fusions simply did not survive to
  the final flow-extracted transcripts. The gate's value is precision *hygiene*
  on copy attribution, **not a headline recall/precision metric move** (Sn/Pr
  literally do not change).
- **The "jaccard alone drops nothing" line is true but tautological.** Every one
  of the 1876 pairs already passed an upstream jaccard prefilter (observed
  **min jaccard = 0.300** across all pairs; 0 pairs below 0.13). So this is not a
  fair head-to-head of two free-standing gates at 0.13 — it is "among
  jaccard-selected candidates, core rejects 155 that jaccard's own scale would
  never reach." The substantive finding (core catches 106 high-jaccard
  domain-sharers) stands; the rhetorical framing slightly oversells it.
- **Scope is narrow and must stay stated.** This measures only de-novo
  *within-chromosome* cross-copy singleton-exon merges (`merge_singletons_by_
  sequence`), per the OOM-safe one-chrom-at-a-time protocol. Cross-chromosome
  paralogs are not co-grouped (and are not this gate's target). The genome-wide
  drop measurement is from traces; the downstream-effect proof is **one contig
  (NC_073228.2) + the prior NC_073235.2** — not all 22 contigs run gate-on/off.
- **Net.** The gate is **strictly safe** (provably loses nothing real on the two
  contigs tested) and **principled** (valley-placed, separates a class jaccard
  cannot), but its demonstrated benefit at scale is **small precision hygiene**,
  not a metric win. Keep it default-off as a measured, optional precision lever;
  do not pitch it as a recall or F1 improvement.
