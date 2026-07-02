# A1 SEDEF re-corroboration — full `final.bed`, depth+size-controlled

**Question.** A1 (`bench/A1_READ_CONSENSUS_O1.md`, commit `396422d`) established a genuinely
read-derived copy signal (SDA / Vollger φ-correlation, no asm20) on the collapsed frontier, but its
**genome-corroboration leg** was downgraded to *"consistent-not-established"*: the read copy count
`K_read` vs a **basic** SEDEF partner count gave Spearman ρ = 0.443, p = 0.044, perm-p = 0.052
(straddles the null), and `K_read` is ~80 % a read-**depth** proxy
(ρ(K_read, SEDEF | depth) = 0.40, p = 0.08). The user's point: the **full** SEDEF map is already
computed (`/mnt/c/Users/jfris/Desktop/final.bed`, 253 029 pairs); A1 used only a basic partner count.
Can the corroboration be pushed **above the null** using the full `final.bed`, with confounds properly
controlled?

**Artifacts.** `bench/a1_sedef_recorroborate.py` → `bench/a1_sedef_recorroborate.tsv` / `.json`.
Run: `PYTHONHASHSEED=0 python bench/a1_sedef_recorroborate.py`. The read side
(`K_read`, `med_depth`, `n_read_psv`, `mean_phi`) is read from the cached A1 `rows.json` — **no SDA
re-run**; only the SEDEF (genome) side is recomputed from `final.bed`.

---

## 1. Method

- **GGO-only parse.** `final.bed` has 34 tab columns (BEDPE: `c1,s1,e1,c2,s2,e2,…,identity(col20),
  …,divergence(col22),…`, verified against the example row). Each GGO–GGO pair contributes two
  features (one indexed on each side). **Both** `c1` and `c2` must be a GGO nuclear chrom (25 scaffolds
  > 1 Mb: 23 autosomes `NC_0732xx.2` + X `NC_086017.1` + Y `NC_086018.1`); the 36 MT rows
  (`NC_011120.1`) are dropped. **252 993 / 253 029** pairs kept. The 7 perfect-identity (`identity=1.0`,
  `div=0`) pairs are **kept** (they are the most-collapsed segdups — exactly the frontier signal).
- **Axis 1 — expected-CN proxy.** For each collapsed backbone locus, `expected_CN` = (# **distinct**
  merged GGO SEDEF partner regions overlapping it) + 1 — a read-independent copy-number estimate.
- **Axis 2 — divergence (cleaner).** `sedef_div` = median SEDEF divergence of the locus's overlapping
  segdup partners — a read-independent measure of how far the copies diverged. Sidesteps the
  copy-vs-allele confound because it corroborates **divergence**, not copy number.
- **Population.** Collapsed-core = A1 families with `collapsed_excess > 0`, **n = 21** (18 `no_psv`,
  3 `partial`). All 21 overlap ≥ 1 GGO segdup, so the divergence axis is fully populated.
- **Confound control.** Every axis is reported raw and as a **rank-based partial Spearman** controlling
  (i) `med_depth`, (ii) locus `size`, (iii) **both**, each with a Monte-Carlo label-permutation
  null (B = 20 000) plus a **size-matched shifted-genomic null** (B = 5 000: recompute the genome
  variable at random same-size windows on the same chrom).

### Why `size` had to be added (the fix over A1)

A1 controlled only `med_depth`. But **locus size drives both sides** while being uncorrelated with
depth — so a count-vs-count agreement can be a pure size artifact:

| pair | Spearman ρ | p |
|---|---|---|
| size ↔ expected_CN | **+0.492** | 0.024 |
| size ↔ n_read_psv | **+0.477** | 0.029 |
| size ↔ med_depth | −0.012 | 0.96 |
| **size ↔ sedef_div** | **+0.001** | 1.00 |
| K_read ↔ med_depth | +0.426 | 0.054 |
| n_read_psv ↔ med_depth | +0.385 | 0.085 |

The last two rows confirm A1's depth-proxy finding. The decisive rows are the first and last-of-genome:
`expected_CN` (a **count**) scales with size, but `sedef_div` (a **rate**) is size-**immune**
(ρ = 0.001). This is exactly why the divergence axis is the more defensible corroboration.

---

## 2. Does the full `final.bed` help Axis 1? — No (it is weaker)

| genome copy-number variable vs `K_read` | Spearman ρ | p |
|---|---|---|
| **cached** basic partner count (= A1 headline) | 0.443 | 0.044 |
| **full** `final.bed` distinct-region `expected_CN` | **0.369** | 0.099 |

The A1 headline reproduces exactly (0.443) from the cached count; the full-map distinct-region count is
**weaker, not stronger**. The basic count correlated better only because raw fragmented partner counts
scale even harder with size/repeat content — i.e. it was **more** confounded, not more signal.

**Axis 1 under confound control (n = 21):**

| control | partial ρ | analytic p | perm-p |
|---|---|---|---|
| raw | 0.369 | 0.099 | 0.099 |
| \| depth | 0.298 | 0.201 | 0.189 |
| \| size | 0.148 | 0.535 | 0.519 |
| **\| depth + size** | **−0.005** | 0.983 | 0.982 |
| size-matched shifted null (raw / \|depth) | — | — | 0.096 / 0.207 |

Controlling depth **and** size, the `K_read`↔`expected_CN` correlation **vanishes to zero**
(ρ = −0.005). The genome copy-number proxy does **not** corroborate `K_read` once size and depth are
removed. Using the full `final.bed` cannot rescue this — it is a confound problem, not a
data-resolution problem. (Dedup n = 20 raw rises to 0.449 by dropping one discordant point, but its
depth+size partial is 0.186, p = 0.46 — same conclusion.)

---

## 3. Divergence axis — the cleaner corroboration, but still marginal

**Axis 2: `n_read_psv` vs `sedef_div` (n = 21):**

| control | partial ρ | analytic p | perm-p |
|---|---|---|---|
| raw | 0.398 | 0.074 | 0.074 |
| \| depth | 0.366 | 0.113 | 0.104 |
| \| size | 0.453 | 0.045 | **0.040** |
| **\| depth + size** | **0.427** | 0.068 | **0.051** |
| dedup n = 20, \| depth+size | 0.422 | 0.081 | — |
| size-matched shifted null (raw / \|depth) | — | — | 0.458 / 0.529 |

The divergence axis **survives** joint depth+size control at ρ ≈ 0.43 — whereas Axis 1 collapsed to 0.
This is the qualitative win: because `sedef_div` is a size-immune **rate**, controlling size *sharpens*
the divergence signal (partial | size = 0.45, perm-p = 0.040) instead of destroying it. Direction is as
hypothesised: **more-divergent segdup partners → more read PSVs**.

**But it is still marginal, not above the null.** The strongest depth+size-controlled result sits right
at the boundary (perm-p = 0.051), and — critically — the **size-matched shifted-genomic null does not
support spatial specificity** (0.46): divergence measured at an arbitrary same-size window predicts a
locus's `n_read_psv` about as well as divergence at its own locus. (This shifted null is itself weak for
a rate variable — random windows are unconditioned on carrying a partner — so the label-permutation is
the primary null; but its failure is an honest negative, not a strength.)

**Axis 2b (`mean_phi` vs `sedef_div`)** is non-significant throughout (raw ρ = 0.213, perm-p = 0.35;
depth+size ρ = 0.203). Read linkage *coherence* does not track partner divergence — consistent with φ
being a within-copy consistency measure rather than a between-copy divergence measure.

---

## 4. Multiple testing

Three primary axes tested (`raw` perm-p): Axis 1 = 0.099, Axis 2 = **0.074**, Axis 2b = 0.353. The best
is the divergence axis at perm-p = 0.074; **Bonferroni × 3 = 0.223**. Nothing clears α = 0.05 after
correction. The single sub-0.05 value in the whole analysis (Axis 2 partial | size, perm-p = 0.040) is
one cell of a 3-axis × 4-control grid and does not survive multiplicity.

## 5. Honest limits (n and reference-derivation)

- **n ≈ 21 caps everything.** 7/21 rows (6 distinct loci; 6/20 after dedup) have `n_read_psv = 0` (no read heterogeneity) and one pair
  (families 169/170, the same `LOC101124683` locus, nested backbones — a **pseudo-replicate** now
  flagged in the TSV and dropped in the dedup variant) is a tied zero point, leaving ~14–15 informative
  points. No SEDEF re-parse can add statistical power here.
- **SEDEF is reference-derived — shared substrate.** `final.bed` is computed from the same T2T assembly
  to which the reads were mapped for recruitment. Both the genome side (`expected_CN`, `sedef_div`) and
  the read *recruitment* (which MAPQ-0 reads pile up) trace to that reference. Only the read
  *partition / het* is read-derived. So even a clean above-null result would be a partly
  shared-substrate corroboration — a structural ceiling on this cross-modal test.

---

## 6. Verdict

Using the full `final.bed` does **not** push A1's genome corroboration above the null. The
expected-CN (copy-number) axis is actually **weaker** with the full map (ρ 0.369 vs cached 0.443) and
**collapses to ρ = −0.005** once locus size and depth are jointly controlled — the A1 headline 0.443 was
itself a size+depth artifact, not copy-count agreement. The **divergence axis is genuinely cleaner and
more defensible**: its genome side is a size-immune rate, it sidesteps the copy-vs-allele confound, and
it *survives* joint depth+size control at ρ ≈ 0.43 in the hypothesised direction (more-divergent
partners → more read PSVs). But it remains **marginal** — best depth+size perm-p = 0.051 (boundary), it
fails Bonferroni (0.223), and it fails the size-matched shifted-genomic spatial-specificity null (0.46)
— capped by the collapsed-n ≈ 21 (~14 informative points, one pseudo-replicate) and by SEDEF sharing the
reference substrate. **A1's bottom line is unchanged: "consistent-not-established" stands.** The one
substantive update is qualitative: the **divergence axis, not copy number, is the corroboration to lead
with** — it is confound-robust where copy number is not, even though neither clears the null.

---

### Determinism & reproducibility
`PYTHONHASHSEED=0` asserted; every permutation uses an explicit `np.random.default_rng(SEED+offset)`,
`SEED = 1234`, `B_perm = 20000`, `B_shift = 5000`; rows sorted before write. Two independent full runs
produced **byte-identical** `.json`, `.tsv`, and stdout. `expected_CN` = distinct merged GGO partner
regions + 1; `sedef_div` = median partner divergence (`final.bed` col 22); GGO-only, MT dropped;
`identity ∈ (0,1]`, `divergence ∈ [0,1]` (perfect-identity segdups retained).

### Independently verified (commit `b9915fc`)
Re-run + independent re-derivation (this pass) confirm the committed result. (1) The script regenerates
**byte-identical** `.tsv`/`.json` under `PYTHONHASHSEED=0` (md5 `51fe243c…` / `d44811de…`). (2) All
headlines reproduce from an **independent** partial-Spearman + permutation implementation on the TSV:
Axis 1 raw ρ = 0.369 (p = 0.099), partial | depth+size = **−0.005** (vanishes); Axis 2 raw ρ = 0.399,
partial | size = 0.453 (perm-p ≈ 0.04), partial | depth+size = **0.427** (perm-p ≈ 0.051, boundary);
Axis 2b NS; confounds size↔expected_CN = 0.492, size↔n_read_psv = 0.477, size↔sedef_div = 0.001;
cached 0.443 vs full 0.369; Bonferroni × 3 = 0.223. (3) The genome side (`expected_CN`, `sedef_div`) was
**recomputed from raw `final.bed`** independently and matches exactly (incl. the large-partner case
fam14 CN = 55, div = 0.0827, 75 pairs); GGO–GGO filtering keeps 252 993 / 253 029, drops all 36 MT rows —
the very first data row is an MT row and is correctly dropped (no parse leak); identity col 20 is written
but never used in any statistic. Only fix this pass: §5 prose `6/21 loci` → `7/21 rows (6 distinct loci)`
(immaterial; ~14 informative points unchanged). **Verdict unchanged: "consistent-not-established"; lead
the corroboration with the divergence axis, not copy number.**
