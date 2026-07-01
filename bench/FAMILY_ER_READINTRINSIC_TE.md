# Can we detect TE-ness DIRECTLY FROM THE READS? (self-contained, no RepeatMasker)

**Scripts / artifacts** — `bench/family_er_readintrinsic_te.py` (deterministic, `PYTHONHASHSEED=0` pinned
in-script; re-execs itself if unset), `bench/family_er_readintrinsic_te.tsv` (one row per gene-pair edge:
degree, read-multimapping, k-mer, shared-length features + `arbitration` flag),
`bench/family_er_readintrinsic_te.json` (per-feature AUC by core band + complementarity + gate sweep +
Pareto + honest caveat). Feature caches: `bench/denovo_family_edges.tsv` (de-novo homology graph → degree),
`bench/ri_readmm_cache.tsv` (per-locus minimap2 `-N50` alignment multiplicity), `bench/ri_kmer_cache.tsv`
(canonical 21-mer multiplicity over `GGO.fasta`), `bench/ri_sharedlen_cache.tsv` (aligned shared length).
Reuses the committed gate machinery `bench/family_er_te_gate.py` / `bench/family_er_oracle_sweep.py`
(`compute_point`, `prepare`) and the label cache `bench/edge_bridge_mask.tsv` (`cls` TP/truthbar/genuine,
`in_ep`, `core`, `bridge_mask` = the soft-mask feature). Reproduce:
`PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/family_er_readintrinsic_te.py`.
**Byte-identical across independent launches** — md5 `.tsv` `6508a83c…`, `.json` `4174010e…`; the heavy
features are precomputed to integer caches with deterministic builders and the gate re-refine is seed-0.

## The question

The committed soft-mask gate (`bench/FAMILY_ER_TE_GATE.md`, commit `e472671`) established that the
RepeatMasker soft-mask is a **real but SOFT** signal: on the arbitration population (`in_ep=0`, direct
`core` edge, mostly gorilla `LOC*–LOC*`) it separates genuine-junk over-merges from protein-missed real TP
at only **AUC 0.694** at the recommended operating band (`core≥0.50`), and **0.777** over the whole re-admit
region (`core≥0.13`), because *both* classes draw on repeat-rich `LOC*` loci. The soft-mask needs an external
RepeatMasker library. **Can a signal read straight off the reads (no library) beat 0.694 — and does it help
the E_r over-merge oracle?**

## The arbitration population (identical to the soft-mask gate)

`in_ep=0` AND a direct `core` edge, genuine-FP vs protein-missed real TP:
**1853 edges = 1763 genuine + 90 protein-missed TP**. This is the exact population where the soft-mask scored
0.694/0.777, so every AUC below is head-to-head comparable.

## The three read-intrinsic signals (all self-contained, no RepeatMasker)

1. **Homology-graph DEGREE** *(topological, read-derived)* — how many other de-novo loci each bridged locus
   links to at `core_recip ≥ 0.13` in the de-novo homology graph (`denovo_family_edges.tsv`). A TE bridge's
   shared sequence connects to MANY loci (repeat hub); a real paralog connects only to its family. Purely
   from the reads — the graph is de-novo, no annotation. Per edge = max/mean of the two endpoints' degree.
   **Not circular with the gate**: the edge `core` is the *weight*; degree is a node-level neighbor *count*.
   Within the arbitration population Spearman(core, degree) ≈ −0.23 (near-orthogonal), and AUC is measured
   *within* core bands (core held roughly constant) yet degree still discriminates. Labels come from
   protein-homology/truth, not from degree.
2. **Read multi-mapping DEGREE** *(the most literal "from the reads")* — count each read's genome-wide
   alignment multiplicity in `GGO_mm.bam` (`minimap2 -ax splice:hq -N 50 -p 0.1 --secondary=yes`, so it CAPS
   at ~50). Per bridge locus = median/mean/max over its reads.
3. **Genome k-mer MULTIPLICITY** — canonical 21-mer multiplicity of the bridge exonic sequence against the
   whole `GGO.fasta`. High = repeat. Mask-agnostic and catches novel/gorilla-specific TEs a human-biased
   RepeatMasker library would miss. (Honest scope note: this one is *assembly*-intrinsic, not strictly
   read-set-intrinsic, but needs no external annotation library.)

## AUC(genuine > TP) by core band vs the soft-mask baseline

| signal | c≥0.13 | c≥0.30 | c≥0.50 **(operating band)** | c≥0.70 |
|---|---|---|---|---|
| **soft-mask** (baseline, needs RepeatMasker) | 0.777 | 0.705 | **0.694** | 0.675 |
| **homology-graph degree** (max) | **0.830** | 0.734 | 0.650 | 0.762 |
| read multi-mapping (best variant) | 0.686 | 0.668 | 0.592 | 0.650 |
| genome k-mer (best variant) | 0.732 | 0.680 | 0.619 | 0.557 |
| *shared-length* (caveat fix, AUC TP>genuine) | 0.821 | 0.640 | 0.521 | 0.482 |

`n_gen / n_tp` per band: 1763/90, 521/64, 182/33, 75/19 (population shrinks with core → high-band AUCs are
noisy; degree's non-monotone 0.83/0.73/0.65/0.76 is this shrink, not a real dip-and-recover).

## Does any read-intrinsic signal BEAT the soft-mask?

- **At the actual operating band (`core≥0.50`, where soft-mask = 0.694): NO single read-intrinsic signal
  beats it.** Degree 0.650, read-mm 0.592, k-mer 0.619 all LOSE. State this against the *correct* baseline:
  the 0.694 anchor, not the 0.777 broad-region number.
- **Homology-graph degree beats the soft-mask only on the BROAD re-admit region** `core≥0.13`
  (0.830 > 0.777) and at `core≥0.70` (0.762 > 0.675) — real, but not at the deployed operating point.
- **Read multi-mapping is near-chance (0.686), and the mechanism is not what one would guess.** It fails NOT
  because multiplicity is low — it fails because **BOTH classes are high-multiplicity**: per-edge median
  `mm_median` ≈ **11** (mean 21), and **~93% of arbitration edges saturate the `-N50` secondary cap**
  (~70% of individual loci hit the cap). The capped count cannot separate genuine from TP because both pile
  at the ceiling. (The earlier "median multiplicity ~1" explanation was wrong and has been corrected.)
- **K-mer multiplicity is intermediate (0.732 @ `core≥0.13`)** — a genuine repeat signal, weaker than degree.

## The load-bearing result: soft-mask + degree COMPLEMENT (they do not replace each other)

The equal-weight **rank-average of soft-mask and degree** (no weights fitted) **exceeds the soft-mask at
EVERY band, including the operating band**:

| band | soft-mask | degree | **combined rank-avg** | Spearman(soft-mask, degree) |
|---|---|---|---|---|
| c≥0.13 | 0.777 | 0.830 | **0.841** | 0.500 |
| c≥0.30 | 0.705 | 0.734 | **0.754** | 0.622 |
| c≥0.50 **(operating)** | 0.694 | 0.650 | **0.740** | 0.227 |
| c≥0.70 | 0.675 | 0.762 | **0.788** | 0.034 |

Spearman ≈ 0.50 → the two signals are only moderately correlated (partly orthogonal), which is exactly why
their combination helps even at the band where degree alone loses. **Degree complements soft-mask; it does
not replace it** — a degree-ONLY gate (no soft-mask) cannot even match the soft-mask recommended operating
point (`degree_only_matches_softmask_recommended = None`).

## The combined self-contained gate through the committed `compute_point` machinery

> **`E_r edge valid  ⟺  in_ep  OR  (core ≥ t  AND  bridge_mask < m  AND  degree < dmax)`**

Grid: `t ∈ {0.13…0.70}`, `m = 0.10` (fixed from the committed soft-mask gate), `dmax ∈ {3,4,5,6,8,10,15,
20,30,∞}` (∞ = degree OFF, reduces to the soft-mask gate). Pareto-beat of the soft-mask recommended point =
`TPpm ≥ ref` AND `edge-host block ≤ ref`, strictly better on one.

| gate | TP (all) | **TPpm** (protein-missed) | ncTP | genuine over-merge | edge-host block | prfam genuine-rate |
|---|---|---|---|---|---|---|
| no gate | 450 | 97 | 31 | 4163 | 0.487 | 0.159 |
| protein-only (`in_ep`) | 353 | **0** | 0 | 11 | 0.020 | 0.0025 |
| **soft-mask recommended** (`t=0.50, m=0.10`) | 374 | **21** | 7 | 55 | 0.0892 | 0.0282 |
| core-knee (`t=0.70`, no mask/deg) | 372 | 19 | 6 | 86 | 0.1121 | 0.0275 |
| **combined `t=0.50, m=0.10, dmax=5`** | 374 | **21** | 7 | **40** | **0.0731** | **0.0259** |

**The combined gate Pareto-beats the soft-mask-only gate:** at **zero real-TP cost** (TPpm 21→21, ncTP 7→7,
truthbar 6142 unchanged) it removes **15 more genuine over-merges (55 → 40)** and lowers the edge-host block
rate **0.0892 → 0.0731** and prfam genuine-rate **0.0282 → 0.0259**. `dmax=5` is **not** a knife-edge: `dmax
∈ {5,6}` give the identical optimum (width-2 plateau), `dmax=4` drops TPpm to 15 (loses 6 real TP), so
`dmax=5` is the *lowest* cap that preserves all real TP. **Honest caveat:** `dmax` is chosen **in-sample**
(single free parameter, coarse grid, no held-out split), so the modest 55→40 gain carries a small overfit
risk. The gain is genuine and `compute_point`-verified, but modest.

## The HONEST high-copy-real-gene-family caveat

**Degree/multiplicity flag ANY high-copy locus as TE, including REAL high-copy gene families** — KRAB-ZNF
(~400 paralogs), olfactory receptors (~1000). So **degree alone cannot separate a 400-copy zinc-finger family
from a TE**. In the arbitration population the top-decile-degree protein-missed TP are exactly these:
`MMP24OS–SPC25` (deg 51), `KRBOX4–ZNF674` (deg 50), `LOC109024259–ZNF224` (deg 22). Two mechanisms keep the
combined gate honest:

1. **Protein-first ordering.** `in_ep` is evaluated *first*, so annotated high-copy families (ZNF/OR) that
   ARE protein-homologous are kept by the `in_ep` branch and never reach the degree branch. Empirically
   **37/41 (90%) of top-decile-degree TP are rescued by `in_ep`**; degree only ever judges the
   protein-*missed* residual.
2. **At the operating band the `softmask<0.10` branch is the real protection, NOT shared-length.** All 9
   high-degree protein-missed TP have `softmask ≥ 0.235`, so the `softmask<0.10` branch excludes every one of
   them regardless of degree. Shared-segment length — a real paralog shares a LONG contiguous gene body, a TE
   bridge a SHORT exonized-TE fragment — separates the two **only on the broad re-admit region** `core≥0.13`
   (AUC TP>genuine 0.821; median shared-len genuine 439 bp vs TP 1493 bp). It is a **sign-inversion at high
   core**: at `core≥0.50` it is at chance (0.521; medians 1601 vs 1593 bp) and at `core≥0.70` it **inverts
   below chance** (0.482; genuine bridges share *more*, 1944 vs 1522 bp). **Shared-length does not rescue the
   high-copy-family confound at the operating band** and is not a term in the deployed gate; it is disclosed
   as a broad-region-only separator.

## Bottom line

**Can we detect TE-ness directly from the reads? Partly — one of the three signals is genuinely useful, and
the useful one is the topological, not the literal, multiplicity signal.** The two *literal* multiplicity
signals are the weak ones: read multi-mapping is near-chance (0.686 — both classes saturate the `-N50` cap,
median multiplicity ~11, not ~1) and genome k-mer multiplicity is intermediate (0.732). The winner is
**homology-graph DEGREE**, a read-derived *graph* signal: it is the strongest read-intrinsic detector
(0.830 @ `core≥0.13`) and it is the one that helps the E_r oracle.

**Is a read-only TE detector viable, and does it help the oracle?** As a *replacement* for the soft-mask, no:
no single read-intrinsic signal beats the soft-mask at the 0.694 operating band, and a degree-only gate
cannot match the soft-mask operating point. As a *complement*, yes and modestly: soft-mask + degree
rank-average exceeds the soft-mask at **every** band (0.740 > 0.694 at the operating band), and the combined
gate `in_ep OR (core≥0.50 AND softmask<0.10 AND degree<5)` **Pareto-beats** the soft-mask-only gate — 15
fewer genuine over-merges (55→40) and lower edge-host block (0.089→0.073) at **zero** real-TP cost. The
honest limit is that **degree/multiplicity confounds real high-copy gene families with TEs**, resolvable only
by protein-homology (`in_ep`, evaluated first) plus — on the broad re-admit region only — shared-segment
length; at the operating point the residual is held by `in_ep` + `softmask<0.10`, so the read-only signal
adds a real but bounded second line of defense rather than eliminating the external mask.
