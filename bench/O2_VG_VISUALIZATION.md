# O2 copy-assignment VARIATION GRAPH — materialized object + figures

`bench/o2_vg_visualization.py` turns the thesis's distinctive **Canzar-flip** framing of copy
assignment (objective **O2**) into a single, concrete, drawable object: the **per-family variation
graph**. It reuses the real catalog machinery and the real per-read gate — nothing here is
synthetic — and renders publication-legible figures for three flagship families that span the range
from *fully resolvable* to the *K-frontier identifiability floor*.

---

## The Canzar-flip framing (what the graph IS)

In classic multi-mapping, a read that maps equally well to *k* paralog copies is a problem: minimap2
leaves it at MAPQ 0 and the naive fix is to hand it out `1/k`. The thesis flips this: a multi-copy
family is **one variation graph**, and the multi-mapping reads are **shared evidence threaded through
it**.

| VG object | Biological meaning |
|---|---|
| **Backbone** (shared consensus path) | the sequence all copies agree on |
| **Bubble** | a PSV — a position where copies differ (≥2 alleles) |
| **SUN bubble** (⭐, yellow/bold allele) | a bubble where **one copy's allele is private** (Sudmant 2010 Singly Unique Nucleotide) → a single read over it **pins that copy deterministically** (a Strong-Separation witness, THEORY.md §5) |
| **Copy-path** | each copy = a colored path threading its allele at every bubble |
| **Read** | a partial path; **assigned** to the copy-path it is consistent with, via the significance gate |
| **TIED read** | a read that spans **no distinguishing bubble** → genuinely unresolvable (the K=0 / K-frontier floor) |

The per-read gate is the theorem made operational: `copy_assign.assign_read` computes, over the
bubbles a read spans, `n_decisive` (# bubbles where the candidate copies actually differ — the
identifiability gate) and a **log-likelihood-ratio margin** over the runner-up copy. A read is
**ASSIGNED** iff `n_decisive ≥ 1` **and** `margin ≥ MARGIN (2.0)` — this is the per-read `min_p`
certificate of THEORY.md Thm 4 / Thm 7. `n_decisive == 0` ⟹ **TIED**. Family-level, the number of
distinct copy-paths through the bubbles is `K = χ(H)` (the graph's minimum copy cover): `K = #copies`
⟹ fully resolvable; `K < #copies` ⟹ some copies collapse; `K = 1` (or 0 bubbles) ⟹ the whole family
is the floor.

Until now the theory *described* this VG; there was no single materialized VG object drawn. This is
that object.

---

## The object (reusable API)

```python
from o2_vg_visualization import materialize_family, render_family, render_graph_view
vg = materialize_family(39)            # dict; cached to /home/juanfra/winloci_scratch/o2vg/fam39.json
render_family(vg, "fig.png")           # matplotlib backbone+bubbles+paths+reads (the deliverable)
render_graph_view(vg, "fig_graph.png") # optional networkx bubble/source-sink view
```

`vg` is a JSON-able dict:

- `copies[]` — `{name, chrom, start, end, gene, hap_group, tier(1/2/3), n_sun, sun_cols[], n_assigned_reads}`
- `bubbles[]` — `{col, pos, hap{copy→allele}, alleles{allele→[copies]}, is_sun, sun_copy, n_read_support}`
- `reads[]` — `{name, best_copy|None, status(assigned|ambiguous|tied|no_cover), margin, n_span, n_decisive, span_cols}`
- scalars — `n_copies, backbone_len, K, cls, group_sizes, n_bubbles, n_sun_bubbles, tiers{}, assignment{}, per_copy{}`

**Per-copy 3-tier ladder** (recomputed here, internally consistent — *not* joined from the possibly
stale `sun_identifiability.tsv`): **T1** = has ≥1 private (SUN) bubble → single-read resolvable;
**T2** = unique hap-row but no SUN → needs a read co-observing ≥2 bubbles (linkage); **T3** = hap-row
identical to another copy → collapses (the frontier).

---

## Machinery reused (all real data, no re-derivation)

- **Copies, backbone, bubbles**: `bench/psv_graph_genomewide.py` helpers (`dedup_copies`,
  `column_alleles`, `sh`) over `/home/juanfra/winloci_scratch/validated_families.tsv` — the same 154
  validated co-located families as `psv_graph_genomewide.json`. Longest copy = backbone; other copies
  aligned `minimap2 -ax asm20`; reads aligned `minimap2 -ax splice:hq`; a bubble = a backbone column
  where all copies are defined, ≥2 distinct alleles, and ≥3 reads support it.
- **Reads + assignment = the REAL gate**: every read is threaded and classified by
  `copy_assign.assign_read` (the exact `n_decisive` identifiability gate + LLR `MARGIN=2.0`), reading
  its base at each bubble from the real `GGO_mm.bam` pileup. The figures' read colors are that gate's
  output, not a re-implementation.
- **SUN / tier** recomputed from the materialized bubbles (a copy's allele private at a column).

Deterministic: `PYTHONHASHSEED=0`, no RNG in layout; re-rendering a family from cache is byte-identical
(verified: fam39 PNG md5 stable across runs).

---

## The three flagship families (verified against the data)

> Note on family identity: the validated catalog groups paralogs by homology, so "RABL2" is a
> **5-copy** family (RABL2A + RABL2B + 3 LOC paralogs), and "ANKRD18" a 6-copy family, etc. Genes and
> tiers below are the freshly materialized values.

### 1. `fam39` RABL2 — **RESOLVABLE end** → `fig_o2_vg_fam39_RABL2.png`
5 copies on **5 different chromosomes** (RABL2A, RABL2B, LOC134756389, LOC101142457, LOC129523511).
**235 PSV bubbles, 195 SUN**; **K = 5** (fully resolvable), **all 5 copies Tier-1**. Every copy has a
private-allele bubble, so the copy-paths never converge. Of **222 reads**, the real gate **ASSIGNED
191 (86%)** spread across all 5 copies (41/31/52/19/48), 4 ambiguous, **0 tied**, 27 no-PSV-cover.
→ *This is the clean case: many SUN bubbles, five distinct colored paths, reads cleanly colored by
their assigned copy.*

### 2. `fam1` RGPD8 — **K-FRONTIER FLOOR** → `fig_o2_vg_fam1_RGPD8.png`
7 copies (RANBP2-paralog / RGPD cluster, one chromosome). **0 read-supported PSV bubbles → K = 1**,
all 7 copies Tier-3 (identical hap-rows). All 7 copy-paths **collapse onto the backbone** — there is
no distinguishing position, so **all 2075 reads are TIED** (0 assignable, even in principle). `MCC =
χ(H) = 1`.
→ *This is the identifiability floor: no bubble to thread, every read gray.*

### 3. `fam22` ANKRD18 — **MID / within-graph frontier** → `fig_o2_vg_fam22_ANKRD18.png`
6 copies (ANKRD18A, ANKRD18B + 4 LOC paralogs). **61 bubbles, 50 SUN**; **K = 4 of 6**, group sizes
`[2,2,1,1]`: **2 copies resolve** (ANKRD18B T1/47 SUN, LOC101142783 T1/3 SUN) and **4 copies collapse
into 2 Tier-3 pairs** (paths converge). The gate **ASSIGNS 77 reads (9%) only to the two resolvable
copies** (51 + 26); the 127 ambiguous reads are the collapsed-pair reads that span a bubble but can't
win the margin; 613 no-PSV-cover.
→ *One graph shows both regimes at once: resolvable paths get colored reads, collapsed paths do not.*

**Honesty / caveats.**
- `no_cover` reads (partial reads that reach **no** PSV column) are a *coverage* artifact, kept
  distinct from a true identifiability tie. They are large for `fam22` (75%) because its backbone
  (copy0 LOC115933254, 122 kb) is far longer than the other copies, so many reads fall outside the
  shared, PSV-bearing region.
- For legibility the figures show a deterministic subset of bubbles (≤30, SUN-prioritized) and reads
  (≤60, proportional by status); **all counts in the titles/annotation are the true totals**.
- Assignment "accuracy" is not claimed here — the figure shows the gate's *decisions* (assigned /
  ambiguous / tied) on real reads; the non-circular accuracy question lives in `copy_assign.py`
  crossval / `SIM_GROUND_TRUTH.md`.
- Family 39's per-copy SUN counts here (RABL2A 19, RABL2B 92) differ slightly from earlier notes
  (~25 / ~153) because the backbone/alignment frame differs; the classification (all Tier-1) is
  unchanged.

---

## Files

| File | What |
|---|---|
| `bench/o2_vg_visualization.py` | materialize + render module (reusable API + CLI) |
| `bench/fig_o2_vg_fam39_RABL2.png` | RESOLVABLE flagship (matplotlib, required deliverable) |
| `bench/fig_o2_vg_fam1_RGPD8.png` | K-FRONTIER FLOOR flagship |
| `bench/fig_o2_vg_fam22_ANKRD18.png` | MID / within-graph-frontier flagship |
| `bench/fig_o2_vg_fam{39,1,22}_*_graph.png` | optional networkx source→sink bubble views |
| `/home/juanfra/winloci_scratch/o2vg/fam{39,1,22}.json` | cached materialized VG objects |

## Run

```bash
PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/o2_vg_visualization.py            # all 3 flagships
PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/o2_vg_visualization.py 39 --fresh # one family, ignore cache
```
Materializing a family from scratch (`--fresh`) re-runs the copy/read alignments (minutes for the
large families); once cached, re-rendering is instant.
