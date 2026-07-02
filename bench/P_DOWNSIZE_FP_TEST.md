# Is `-p 0.1` a false-positive source? Would `-p 0.8` remove noise?

**2026-07-01.** GGO_mm.bam was aligned `minimap2 -ax splice:hq -uf --eqx -Y -N50 -p 0.1 --secondary=yes`. Question: does the permissive `-p 0.1` (min secondary-to-primary score ratio) inflate false positives (spurious multi-mappings → spurious read-conflict edges), and would "downsizing" to `-p 0.8` (minimap2 default) clean it? Scripts: `bench/p_downsize_fp_test.py` (+ `.tsv`/`.json`). Deterministic (PYTHONHASHSEED=0, crc32 sampling, sorted writes; two byte-identical runs).

## Key fact that makes this testable without re-aligning

`-p` **only filters secondary alignments — it never changes a primary or re-aligns anything.** Verified empirically: a 200-read sample re-mapped at `-p 0.1` vs `-p 0.8` has **200/200 identical primaries** (pos/AS/CIGAR). So `-p 0.8` is faithfully simulated by **post-filtering** the existing BAM — no cluster re-align.

**Proxy validation** (200-read sample, re-mapped with the exact original params at `-p 0.8`):
- AS-ratio ≥ 0.8 post-filter reproduces real `-p 0.8`: **recall 0.936, precision 0.958, 92.5% exact per-read set match**.
- `s1`-ratio ≥ 0.8 (chaining score — what `-p` literally thresholds): **recall 0.984, precision 1.00, 99% exact** — the faithful predictor. (The family correspondence below uses AS-ratio; the ~5% proxy gap does not change any conclusion.)

## 1. Noise volume — at the RAW level, `-p 0.1` *is* a large wide-net

Genome-wide (1/32 crc32 subsample, 137,742 reads / 195,981 secondaries):
- **`-p 0.8` would drop 90.8% of all secondaries** (s1-ratio: 92.7%); `-p 0.5` drops 71.4%.
- The AS-ratio histogram piles at 0.1–0.4 (49,916 secondaries in [0.1,0.2) alone); only ~9% survive ≥0.8.

So as a **raw** conflict-graph input, `-p 0.1` admits an enormous low-ratio secondary mass and inflates raw multi-mapping degree/edges. Taken alone, yes — it is a noisy wide-net.

## 2. …but the de-tie gate already filters it — and on the right axis

Co-located families (165 families with multimappers; 174,749 in-family secondary placements). Raw multi-mapping degree **3.13 → 2.10** per read under `-p 0.8`; 516 of 2,491 raw de-tie edges would drop. But the gate is **stricter inside families than `-p 0.8`**: `-p 0.8` drops only **48.0%** of in-family secondaries vs the **gate's 65.1%** — because sibling copies are near-identical (high score-ratio, survive `-p 0.8`) yet the gate still cuts them on divergence (`|Δde|`).

### Gate-vs-`-p` correspondence (2×2 on in-family secondaries)

| | gate-reject | gate-keep |
|---|---|---|
| **`-p0.8` drop** | 54,524 (junk both kill) | **29,443 GENUINE-TIE LOSS** |
| **`-p0.8` keep** | **59,224 (junk gate catches, `-p0.8` misses)** | 31,558 (both keep) |

- Of what `-p 0.8` drops: **64.9% is junk the gate also rejects**, but **35.1% (29,443 placements = 291 gate-kept de-tie edges = 13.8% of all 2,112 gate edges) are GENUINE ties the gate keeps and `-p 0.8` would destroy.**
- **Neither filter subsets the other.** The gate additionally rejects **59,224** divergent secondaries that `-p 0.8` *keeps* (high score-ratio but `|Δde|>0.005` or `de>0.05`). The gate filters on `de` (divergence) — the correct axis — catching junk `-p 0.8` lets through.
- The final **gate** edge set (**2,112**) is *larger* than the `-p 0.8` edge set (**1,975**) — and that largeness is itself the evidence that `-p 0.8` destroys genuine edges the gate keeps.

### What the genuine loss actually is

The 13.8% of de-tie edges `-p 0.8` would destroy are **near-identical (≤5% divergent) co-located copies** — high score-ratio siblings — concentrated in a few large paralog arrays (**~65% is one 47-copy GSTM2 family**; also LOC109025447-T3, and other SUN Tier-1 *and* Tier-3 families). These are exactly the collapsed / near-identical copy ties the copy-assignment targets.

Sharper point: **genuinely divergent (≈70–80% identity) multimappers are rejected by the *gate* itself** (`DE_MAX=0.05`), independent of `-p`. So the **gate — not `-p` — is the binding constraint on divergent-copy ambiguity.** `-p 0.8`'s damage is specifically to the *near-identical* ambiguous ties, the hardest and most valuable cases.

## Verdict

- **Is `-p 0.1` an FP source?** At the raw-graph level, yes — it admits ~91% low-ratio secondaries and inflates raw degree/edges.
- **Does the gate already make it safe?** Effectively **yes.** The shipped de-tie gate (`DELTA=0.005`, `DE_MAX=0.05`) rejects 64.9% of what `-p 0.8` would drop *plus* 59,224 more divergence-junk, on `de` rather than raw score. The final gated edge set is not inflated by the `-p 0.1` junk — the gate does `-p`'s job and does it better.
- **Should we downsize to `-p 0.8`?** **No — net harmful.** It is mostly redundant with the gate (54,524 overlap), catches *less* junk than the gate (misses 59,224 divergent placements the gate rejects), and would additionally **destroy 29,443 genuine gate-kept ties / 291 de-tie edges (13.8%)** — the near-identical collapsed ambiguity copy-assignment exists to resolve.

**Bottom line: `-p 0.1` was a safe wide net; the de-tie/significance gate is the real filter, and it filters on the correct axis (`de`). Keep `-p 0.1 + gate`; do not adopt `-p 0.8`.** This confirms the "wide net + gate" design rather than refuting it.
