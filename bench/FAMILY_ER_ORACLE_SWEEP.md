# E_r edge-oracle sweep: an operating point that keeps the divergent-paralog recall and cuts the over-merges

**Scripts / artifacts** — `bench/family_er_oracle_sweep.py` (deterministic, PYTHONHASHSEED=0 pinned in-script),
`bench/family_er_oracle_sweep.tsv` (one row per operating point), `bench/family_er_oracle_sweep.json`
(full curve + recommended point + tautology note + disproportionate-drop Fisher test + named
survivors/casualties). Reuses `bench/family_er_pr.py` (E_r edge set, FP split, E_p truth, γ-quasi-clique
refine) and `bench/genome_family_def.py` (the refine operator). Reproduce:
`PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/family_er_oracle_sweep.py`. TSV and JSON are
byte-identical across independent launches.

## The problem

The reframed RNA family `E_r` (transcript-homology, `core_recip>=0.13`, γ-quasi-clique refined) is a
**precision/recall trade-off**, not a strict win. Its 10,755 refined gene-pair edges split into:

| class | n | what it is |
|---|---|---|
| **TP** | 450 | in the cDNA-90% truth (`in_dna_loose`) |
| **truthbar** | 6,142 | protein-homologous divergent paralogs *below* the cDNA-90% bar — **the E_r win** (the truth is incomplete, these are real) |
| **genuine** | 4,163 | real over-merges — unrelated genes bridged by the loose 0.13-core POA homology (AMY2A–ZNF91, RPL14–ZNF669) |

Baseline genuine edge-precision `1 − 4163/10755 = 0.613`; family-level over-merge (PRFAM-mixing blocks)
15.9% (= the gene-graph analog of the 13.4% headline in `FAMILY_ER_PRECISION_RECALL.md`). The loose core
edge bridges unrelated genes; the task is to gate those edges to cut the over-merges **without** dropping
the divergent-paralog recall the reframe exists for.

## Two honesty caveats that gate every number below

**(A) The protein route's "100% recall preservation" is a TAUTOLOGY, not a discovery.** `family_er_pr.classify_fp`
*defines* `truthbar-FP := protein_homologous AND id_ok AND cov_ok`, so **truthbar ⊆ in_ep by construction**
(verified: 0 violations — a violation would be a bug, not evidence). Any operating point that keeps every
`in_ep` edge keeps every truthbar edge *definitionally*. Read "the protein route preserves 100% of the E_r
win at zero cost" as **"zero cost to the protein-homologous PROXY of the win."** The genuinely non-circular
axis measured here is the recovery of the **97 protein-missed TP** (real per the cDNA truth, yet `in_ep=0`)
via the protein-independent POA core signal. Because any non-protein-homologous divergent paralog is
definitionally binned `genuine` and cut, the **`genuine` count is an upper bound** on true over-merges and
the recall of non-protein paralogs is **unmeasured**; the 97 protein-missed TP prove that hidden category
is non-empty and bound how large it can be.

**(B) The block over-merge rate has two faces — both are reported, never just the flattering one.**

| metric | what it counts | recommended point | why it differs |
|---|---|---|---|
| `genuine_rate_block` (PRFAM-mixing) | a refined block spans ≥2 distinct non-mega PRFAM | **2.75%** | matches the MD headline but is **blind to LOC\*/no-PRFAM bridges** |
| `genuine_rate_block_edgehost` (direct) | a block still co-hosts both endpoints of a surviving genuine edge | **11.2%** | the honest family-contamination measure |

At the recommended point **76/86 (88%)** of residual genuine over-merges involve a LOC\*/no-PRFAM gene, so
the PRFAM-mixing metric can't see them — that is exactly why edgehost (11.2%) is ~4× higher than
PRFAM-mixing (2.75%). Family contamination at the recommended point is **~11%, not 2.75%**.

## The knobs

1. `core_recip` threshold `t` on the per-edge POA transcript-homology core (from `denovo_family_edges.tsv`,
   DN-locus pairs mapped to genes by the same floored `OV_FLOOR=0.20` projection as `family_er_pr`; a gene
   pair's core = **max** over its direct co-blocked DN edges; transitive-only pairs have no direct edge →
   `core=None` → fail any `t`).
2. `protein` gate (`in_ep`): reciprocal whole-protein mmseqs edge OR same non-mega PRFAM.
3. `id` floor (proxy) — shown to be the **wrong** knob (destroys the divergent-paralog recall).
4. **combined** `in_ep OR core_recip>=t_high` — protein confirms coding paralogs; a high transcript-core
   confirms the non-coding/protein-missed real families the pure protein gate drops; unrelated bridges
   (no protein AND low/no core) are cut. **This is the recommended family.**

## The trade-off curve (key operating points)

| gate | t_high | n_edges | TP | TPnoEp /97 | **ncTP /31** | truthbar | genuine | g_edge | **g_blk prfam / edgehost** | recall .90–.95/.95–.99/≥.99 |
|---|---|---|---|---|---|---|---|---|---|---|
| none | – | 10755 | 450 | 97 | 31 | 6142 | 4163 | 0.387 | **0.159 / 0.487** | 0.951/0.923/0.972 |
| protein | – | 6506 | 353 | 0 | 0 | 6142 | 11 | 0.002 | **0.003 / 0.020** | 0.738/0.731/0.832 |
| id≥0.80 | – | 478 | 450 | 97 | 31 | **9** | 19 | 0.040 | 0.098/0.078 | 0.951/0.923/0.972 |
| combined | 0.30 | 7091 | 417 | 64 | 20 | 6142 | 532 | 0.075 | **0.111 / 0.333** | 0.869/0.815/0.958 |
| combined | 0.50 | 6721 | 386 | 33 | 11 | 6142 | 193 | 0.029 | **0.062 / 0.213** | 0.836/0.754/0.909 |
| **combined** | **0.70** | **6600** | **372** | **19** | **6** | **6142** | **86** | **0.013** | **0.028 / 0.112** | **0.803/0.739/0.881** |
| combined | 0.80 | 6558 | 362 | 9 | 2 | 6142 | 54 | 0.008 | 0.019 / 0.080 | 0.738/0.731/0.874 |

Two structural facts from the sweep: (i) the **id gate is wrong** — at id≥0.80 truthbar collapses 6142→9
(it destroys the divergent-paralog recall the reframe is for); (ii) **truthbar is invariant (6142) for every
protein-containing gate** — the tautology of caveat (A), not a discovery.

## Recommended operating point

> **`E_r edge = core_recip>=0.13 AND (in_ep OR core_recip>=0.70)`** — gate = **combined, t_high = 0.70**.

This is the **task-literal knee**: `max(TP + truthbar)` subject to PRFAM-mixing block over-merge ≤ 3%. Among
all gated points with block ≤ 3% it maximizes TP+truthbar (372+6142 = 6514, vs combined-t0.80 6504, vs pure
protein 6495).

**Headline numbers**

- **Precision.** Genuine edge-precision **0.613 → 0.987** (genuine over-merge 4163 → 86, a **97.9% cut**).
  This is **above the old read-conflict-refined `E_c` 0.94** — the recommended `E_r` recovers the precision
  the reframe had lost (0.61) *and then some*, while keeping the divergent-paralog reach `E_c` never had.
- **Over-merge, both metrics.** PRFAM-mixing block over-merge **15.9% → 2.75%**; direct edge-host block
  over-merge **48.7% → 11.2%**. Cite both; real family contamination is ~11%.
- **Divergent-paralog recall (the win).** truthbar **6142/6142 = 100% retained** — but this is
  **TAUTOLOGICAL** (truthbar ⊆ in_ep), so the honest phrasing is "100% of the *protein-homologous proxy* of
  the win, at zero cost."
- **Protein-missed TP recovered by the core route (the non-circular win).** **19/97**, of which **6/31**
  non-coding TP.
- **Reachable recall** (on the reachable cDNA-truth subset): overall **0.811** (from 0.949 baseline);
  id-stratified 0.803 / 0.739 / 0.881 for .90–.95 / .95–.99 / ≥.99. High-identity paralogs are mostly
  retained; the loss is concentrated in the low-id / protein-missed tail.

**Named over-merges it CUTS** (genuine, `in_ep=False`, low/no core): **AMY2A–ZNF91, RPL14–ZNF669**,
AMY2A–ZNF141, DOCK1–LOC129529082, THNSL1–LOC109023154, ZNF232–LOC101153526.

**Named paralogs it KEEPS**: ENO2–ENO3 (core 0.628), OGDH–OGDHL (core 0.551), **ZNF234–ZNF45**,
ZNF300–ZNF585A, ZNF529–ZNF665, ZNF616–ZNF675 (all protein-homologous → truthbar → kept).

## Honest residual cost (do not omit)

1. **Non-coding families are dropped disproportionately.** All 31 real non-coding TP are protein-missed, so
   the core route is their *only* path — and at t_high=0.70 it recovers just **6/31**. Non-coding TP drop
   rate **0.806 vs coding 0.127** — a **6.4× higher** drop rate (Fisher **OR = 28.8, p = 8.3e-16**). Named
   non-coding TP lost: **RFPL3–LOC134758217, FRG1–LOC129531431**, and 23 more (mostly lncRNA–protein_coding
   and lncRNA–lncRNA pairs). If non-coding families matter, this point is the wrong pick.
2. **The old `noncoding_families_retained=20` headline was inverted and is retracted.** Of the non-coding-
   *dominant* blocks at the recommended point (15), only **4 host a real TP edge**; the rest are genuine-only
   over-merges (LOC/pseudogene bridges). The report now tracks **real non-coding TP retention (pair-level,
   6/31)**, not a block-dominant count that was 75%+ over-merges.
3. **The `genuine` count is an upper bound** (caveat A): non-protein-homologous divergent paralogs are
   swept into it and cut, and their recall is unmeasured.
4. **TP cost.** 78 of 450 TP dropped (64 of the 97 protein-missed still cut: 7 transitive-unrecoverable +
   the rest below core t_high=0.70).

## Alternatives (explicit value judgments, not the recommendation)

| choice | gate | recovers | block prfam / edgehost | when to pick |
|---|---|---|---|---|
| **precision-first** | pure protein | 0/97, 0/31 nc | 0.003 / 0.020 | if only protein-homologous families matter; max precision |
| **recommended** | combined t=0.70 | 19/97, **6/31 nc** | **0.028 / 0.112** | task knee: precision ≤3% block-rate, keep all divergent paralogs |
| **edge-budget** | combined t=0.50 | 33/97, 11/31 nc | 0.062 / **0.213** | more recovery, but PRFAM-mixing block 6.15% **> the 3% target** — an *edge*-genuine-rate≤3% budget, **not** the block knee |
| **non-coding-preserving** | combined t=0.30 | 64/97, **20/31 nc** | 0.111 / 0.333 | if lncRNA/non-coding paralogs are the goal; **fails** the precision target |

The combined-t=0.50 point (previously mis-headlined as "the knee") is a looser **edge**-genuine-rate≤3%
budget whose **block** over-merge (6.15% PRFAM / 21.3% edge-host) exceeds the task's ≤3% block target by
>2×; it is reported here only as a recovery-favoring alternative.

## Bottom line vs the old measurements

| definition | genuine edge-precision | family-level over-merge | divergent-paralog reach |
|---|---|---|---|
| old read-conflict-refined `E_c` | 0.94 | — | **no** (misses divergent paralogs) |
| reframed `E_r` (ungated, `core≥0.13`) | 0.613 | 15.9% prfam / 48.7% edge-host | yes (6142) |
| **`E_r` + recommended gate (`in_ep OR core≥0.70`)** | **0.987** | **2.75% prfam / 11.2% edge-host** | **yes (6142, protein-homologous proxy)** |

The one concrete oracle change — replace the flat `core_recip>=0.13` edge with
**`core_recip>=0.13 AND (in_ep OR core_recip>=0.70)`** — restores precision *past* the old read-conflict
0.94 while keeping the divergent-paralog reach the reframe was built for. The residual, stated plainly: it
keeps that reach only for **protein-homologous** paralogs (tautological), buys back just **19/97** of the
protein-missed TP and **6/31** of the real non-coding families, and its true family contamination is the
**11.2% edge-host** figure, not the 2.75% PRFAM-mixing one.
