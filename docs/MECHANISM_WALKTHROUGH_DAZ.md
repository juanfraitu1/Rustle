# Mechanism walkthrough — real DAZ reads, end to end

**Question this answers (for the advisor):** for a real multi-copy locus, (1) how do we use the `de` tag to
decide which copy a read belongs to, with every number shown; and (2) once a read is assigned, is it accepted
into a bundle, and is max-flow applied to assemble the copy's transcripts?

**Everything below is reproducible** from the gorilla IsoSeq BAM:
`samtools view GGO.bam NC_073248.2:42700000-43000000` (the DAZ cluster on chrY).

## The locus (real, annotated)
Two paralogous copies, ground truth known:
- **DAZ1** — NC_073248.2:42,783,133–42,859,657, **− strand** (reads align reverse).
- **DAZ3 / LOC129530216** — 42,879,918–42,945,552, **+ strand** (reads align forward).

They are near-identical inverted amplicons, so most reads **multi-map** to both copies (MAPQ 0). The job: decide,
per read, which copy it came from — then assemble each copy.

---

## Step 1 — Read-level copy assignment with `de`

minimap2 puts three relevant numbers on each alignment of a read:
- `NM:i:` — raw edit distance (substitutions **+ every indel base**).
- `de:f:` — **gap-compressed** per-base divergence (each indel counts **once**, not per base).
- `AS:i:` — alignment score.

We work in **events** = `de × aligned_length` (the gap-compressed count of distinguishing differences). The read
belongs to the copy with **fewer events**; the discriminant is `ΔEvents = events(sibling) − events(thiscopy)`.

### Three real reads (every value straight from the BAM)

**(a) Decisive — read …46596327** (all metrics agree → DAZ1)

| | `de` | `NM` | `AS` | aligned len | **events = de·len** |
|---|---|---|---|---|---|
| DAZ1 | 0.0004 | 1 | 1686 | 2359 | **1** |
| DAZ3 | 0.0063 | 14 | 1519 | 2205 | **14** |

`ΔEvents = 14 − 1 = +13` → **DAZ1** (13 distinguishing events fewer at DAZ1).

**(b) The one that proves why `de`, not `NM` — read …53870980**

| | `de` | `NM` | `AS` | aligned len | **events = de·len** |
|---|---|---|---|---|---|
| DAZ1 | 0.0278 | **88** | 1588 | 2865 | **80** |
| DAZ3 | 0.0229 | **565** | 1381 | 2370 | **54** |

Raw `NM` says **DAZ1** (88 ≪ 565) — but DAZ3's 565 is **mostly indel slippage** (homopolymer 1-base insertions).
Gap-compress it and DAZ3 is the **better** fit: `ΔEvents = 54 − 80 = −26` → **DAZ3**. Raw edit distance would
**mis-assign this read**; `de` corrects it. *(Across this locus, `NM` is ~87% indel bases, and **8 of 177**
two-copy reads are assigned to a different copy under `de` than under `NM`.)*

**(c) Non-identifiable — read …24709554** (the honest floor)

| | `de` | `NM` | `AS` | aligned len | **events** |
|---|---|---|---|---|---|
| DAZ1 | 0.0000 | 0 | 2518 | 2518 | **0** |
| DAZ3 | 0.0000 | 0 | 2518 | 2518 | **0** |

`ΔEvents = 0`. The read covers a region **identical** between the copies — there is no evidence. We do **not**
guess; it stays split by the prior (Step 2).

---

## Step 2 — From evidence to a soft assignment (the EM)

We never hard-assign. Each read's per-copy events become a **likelihood**, and a posterior **responsibility**
`γ` (the EM E-step is `posterior = softmax(log_score + log_prior)`, `src/rustle/vg.rs:2436`). Using the event
model (each distinguishing event carries log-odds `ln((1−ε)/ε) ≈ 2.94` at ε=0.05):

| read | ΔEvents | **γ(DAZ1)** | **γ(DAZ3)** |
|---|---|---|---|
| (a) …46596327 | +13 | **1.000** | 0.000 |
| (b) …53870980 | −26 | 0.000 | **1.000** |
| (c) …24709554 | 0 | **0.500** | 0.500 |

The tie (c) splits 50/50 by the prior — apportioned, never fabricated. The EM iterates (the per-copy prior
updates from the total responsibility mass) until `γ` converges. *(On real DAZ: `[VG-FP-EM] converged in 7
iter, 179 reads adjusted`.)*

`de` is the **interpretable, by-hand** form of the same evidence the production fingerprint-EM scores from the
copies' distinguishing positions — both point the same way.

---

## Step 3 — Reads accepted into bundles, then **max-flow** assembly

This is the part after assignment.

1. **Both bundles keep the read.** A multimapper sits in **DAZ1's bundle and DAZ3's bundle** simultaneously
   (initial weight `1/NH`). Nothing is discarded.
2. **The EM rewrites the read's weight in each bundle to its responsibility γ** — in place:
   `bundles[bi].reads[ri].weight = γ` (`src/rustle/vg.rs:2480`, and `:1351` for the graph-compat EM).
   - read (a): weight **1.0** in DAZ1's bundle, **0.0** in DAZ3's.
   - read (b): weight **0.0** in DAZ1's bundle, **1.0** in DAZ3's.
   - read (c): weight **0.5** in each.
3. **Each copy's bundle becomes a splice graph**, where every node/edge **coverage = Σ of the contributing
   reads' (now-reweighted) weights**. So a copy's graph carries only the evidence the EM assigned to it; the
   spillover reads (weight ≈ 0 in the sibling) add ≈ no coverage there.
4. **Max-flow extracts the transcripts.** The assembler runs network max-flow on each copy's splice graph —
   `push_max_flow_seeded_full` / `long_max_flow_seeded_with_used_pathpat` (`src/rustle/path_extract.rs:11–13`):
   it repeatedly pushes the **heaviest source→sink path** (an isoform), subtracts that flow, and repeats, so
   each copy's isoforms are the high-flow paths through its **reweighted** graph.

**Net effect:** DAZ1's graph carries the full weight of the DAZ1-decisive reads → DAZ1's isoforms are recovered
at full abundance. DAZ3's graph carries only its genuinely-assigned reads → DAZ3 is assembled at its **honest**
(low) abundance, and the ~150 spillover reads that *align* to DAZ3 but belong to DAZ1 contribute ≈ 0 flow, so
they do **not** fabricate phantom DAZ3 isoforms.

---

## The whole chain, on read …53870980 (the hard one)

```
BAM tags         de(DAZ1)=0.0278  de(DAZ3)=0.0229   (NM 88 vs 565 — misleading: indel-inflated)
  ↓ events = de·len
events           DAZ1 = 80        DAZ3 = 54
  ↓ ΔEvents = 54 − 80 = −26  (DAZ3 fits better)
EM responsibility γ(DAZ3) = softmax(...) = 1.000
  ↓ reads[ri].weight = γ
bundle weight    DAZ3 bundle: weight 1.0   |   DAZ1 bundle: weight 0.0
  ↓ coverage = Σ weights  →  splice graph
max-flow         heaviest source→sink path through DAZ3's graph
  ↓
result           this read's evidence builds a DAZ3 isoform, not a DAZ1 one
```

Every value is checkable: the tags from `samtools view`, the EM/weight step at `vg.rs:2436/2480`, the max-flow
at `path_extract.rs:11`. The point for the advisor: **the copy decision is a concrete per-read calculation on
the `de` tag (not a black box), it is soft (ties are apportioned, not invented), and the assigned weight is
exactly what the max-flow assembler integrates into each copy's transcripts.**
