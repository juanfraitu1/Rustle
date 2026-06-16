# Copy-split in-silico validation (collapsed-tandem regime + identifiability boundary)

`bench/copy_split_sim.py` -> `bench/copy_split_sim.png`

Definitive controlled simulation that the joint **read-coherence + PSV copy-split**
separates co-located / collapsed-tandem paralog copies that read-coherence alone
cannot, and an empirical map of the **identifiability boundary** (the thesis
identifiability theorem). The split logic is a faithful line-by-line port of
`split_readchain_by_psv` (+ `select_identifiable_copies`, `read_consistent_with`,
`distinct_columns`, `consensus_haplotype`) from
`src/rustle/vg_family/copy_split.rs`.

## Model

- **N copies share ONE intron chain** (the collapsed-tandem regime). Read-coherence
  (the structural axis) keys transcripts by exact ordered intron chain, so the whole
  family collapses to a single transcript — this is the adversarial case where
  structure alone is blind.
- **P PSV columns** carry the only separating signal. Each copy gets a fixed allele
  vector, rejection-sampled so any two copies differ at >= K columns (the
  identifiability premise).
- **Reads** are simulated per copy (`C` reads/copy). Each read covers a random
  **contiguous sub-interval** of the transcript (models finite read length L /
  partial spanning): it observes only the PSVs inside its span (others = `None`).
  Each observed allele is corrupted at per-base error `e` (HiFi ~ 0.003).
- A copy is **recovered** iff enough error-free reads spanning >= K of its columns
  land on its exact allele vector to reach `min_reads`, and it stays >= K-distinct
  from every already-admitted copy (greedy `select_identifiable_copies`).

## Determinism

Fully deterministic. `random.seed` / `numpy.random.seed` are set from **fixed
integer constants only** (`MASTER_SEED=20260616`, `DEMO_SEED=11`,
`SWEEP_SEED_BASE=700`); every sweep cell re-seeds from a constant + cell offset.
Never seeded from the wall clock. Verified byte-identical across reruns:
- stdout: identical
- `copy_split_sim.png` md5 stable = `1e1e2a32b9032f6f0ea46e757ef5603a`

## (1) Single-config demo  (N=3, P=8, C=10, e=0.003, K=2, min_reads=2)

```
read-coherence ALONE -> 1 transcript  (the whole family collapses)
read-coherence + PSV -> 3 identifiable copies
truth = 3 copies; fraction recovered = 1.000
  copy reads: 10, 10, 9   (all identifiable=True)
```

The PSV axis recovers all 3 copies where structure alone gives 1.

## (2) Sweep — identifiability boundary

N=3, P=8, K=2, e=0.003, min_reads=2, 40 seeded replicates/cell.
Sweep over read span (12.5% .. 100% of PSVs) x coverage (1 .. 32 reads/copy),
plus an error-`e` sweep (0 .. 0.20) at fixed generous geometry.

- **Panel A** heatmap: fraction of copies recovered. Sharp navy(merged) ->
  green(recovered) transition; white contour = 0.99 boundary.
- **Panel B**: recovery vs the theorem's order parameter = per-copy clean
  (>= K-spanning, error-free) support for the candidate vector. A clean sigmoid:
  ~0 below `min_reads`, snapping to 1.0 once support clears the candidate
  threshold (green dashed). This is the empirical identifiability boundary.
- **Panel C**: error floor. Flat at 1.0 through HiFi e=0.003; degrades as `e`
  rises (errors fragment support off the exact true vector below `min_reads`).

### Where the boundary sits (key numbers)

- Full recovery at **full span** needs **coverage >= 3 reads/copy**.
- Full recovery at **max coverage** needs **span >= 62% of the PSVs**
  (a read must reliably span >= K columns above the error floor).
- **Error floor**: recovery = **1.000 at e=0.003 (HiFi)**; holds to ~e=0.01,
  then degrades, falling to **0.550 at e=0.20**.

## Interpretation (the theorem, empirically)

Recovery -> 1.0 exactly when the expected count of error-free reads per copy that
span >= K distinguishing columns crosses `min_reads`; below that the chain stays
**merged** (non-identifiable). Both the geometric axis (coverage x span must put
>= min_reads clean reads on each copy's >= K columns) and the error axis (per-base
`e` erodes that clean support) reproduce the predicted boundary. This is the
in-silico validation of the read-coherence + PSV copy-split and its
identifiability bound.

---

## Verdict

**The method works in-silico in its target regime, the recovery boundary matches
the identifiability theorem, and the real-data win is blocked only by coverage —
not by the algorithm.**

**(1) Real GGO co-located families — coverage-limited, valid NEGATIVE.**
Both real autosomal co-located/tandem clusters were tested (`bench/copy_split_colocated.md`)
and recovered **0 splits (K=2 and K=3, min_reads=2)** — but the cause is sparsity,
not a method failure. The two clusters bracket the failure modes:
- **LOC101144552 cluster**: putative copies are **divergent paralogs (59-73 % pairwise
  id)**, not a tight family, with only **29 real primaries** (212 alignments, but
  **182 = 86 % are uninformative MAPQ-0 cross-mapped secondaries**), and **0** shared
  intron chains. Structure- AND sparsity-limited.
- **LOC101132628 cluster (736/737)**: the genuine clean tandem pair — **99.06 %
  genomic id, 35 direct PSVs (18 exonic)** — but only **10 alignments total / 7 primary /
  6 PSV-spanning** over the whole cluster. The one real multimapper event (3 molecules,
  PRIMARY-MAPQ60 on 737 / SECONDARY-MAPQ0 on 736) **is** PSV-resolvable (**8 decisive
  exonic columns, 34/34 vs 26/34**) — but 3 reads is far below phasing depth, so this is
  per-read **assignment**, not a copy **split**.
This matches the task's prior: the deep-tandem families that need this method (DAZ
collapsed at 3.4x; RBMY unannotated) are exactly the ones GGO does not cover at depth.
GGO marks the boundaries of the multimapping well, it does not contain the target case.

**(2) Simulation — the split DOES recover what read-coherence merges, and the
boundary matches the theorem.** In the collapsed-tandem regime (N copies sharing ONE
intron chain), read-coherence alone collapses to **1 transcript**; read-coherence + PSV
recovers all copies (demo N=3 -> **3 copies, fraction 1.000**). The faithful port of
`split_readchain_by_psv` recovers a copy **iff >= min_reads error-free reads span >= K of
its distinguishing PSV columns** — exactly the identifiability condition. The sweep
confirms the boundary empirically:
- full span needs **coverage >= 3 reads/copy**;
- max coverage needs **span >= 62 % of PSVs** (a read must reliably span >= K columns);
- **error floor: recovery = 1.000 at HiFi e=0.003**, holding to ~e=0.01, degrading to
  **0.550 at e=0.20**.
Panel B is the clean order-parameter sigmoid: recovery snaps to 1.0 once clean
>= K-spanning support per copy clears the `min_reads=2` candidate threshold, and is ~0
(merged) below it. **Recoverable iff enough reads span >= K PSVs above the error floor** —
the theorem, reproduced.

**(3) Thesis implication.** The read-coherence + PSV copy-split is **validated in its
target regime in-silico**, with a recovery boundary that **matches the identifiability
theorem on both the geometric (span x coverage) and error axes**. The real-data
demonstration is **gated on coverage, not on the method**: it needs a **deep co-located
long-read dataset** in the tight-tandem regime (e.g. **testis HiFi for DAZ/RBMY**, where
the high-identity arrays are expressed at depth). GGO usefully **marks the boundaries** of
where the method can and cannot fire — RABL2 copies are too separated (separate coordinate
frames, no within-chain collision), and its co-located tandems (LOC101144552/LOC101132628)
are either too divergent to be a real family or too shallow (6-7 informative reads) to
phase. The honest position: **method correct and boundary-matched in silico; real win
awaits a depth-adequate co-located dataset.**
