# Unbiased differentiators of REAL family vs OVER-MERGE/bridge edges

**Question.** Beyond the two known levers (junction-position co-threading; same-strand filter),
are there OTHER *natural* differentiators — especially graph-topology ones (edge betweenness) —
discoverable label-free on the full 12,212-edge both-spliced cDNA-homology population? Found via
label-free intrinsic bimodality, oriented only with sparse independent labels (12 compara_pos
positives; antisense/overlap negatives), and verified against the two killer confounds
(component-size and node-degree) plus a whole-component grouped null.

All numbers below were re-measured with sklearn 1.5.2 in the sqanti3 env on
`bench/family_def_features.tsv` (scripts: `family_def_bimodality.py`, `family_def_edgebetw_probe.py`,
`family_def_confound_resolve.py`, `family_def_residual_bimodality.py`, `a2_unsupervised_structure.py`,
`a2_edgebetw_orient.py`, `a2_confound_confirm.py`, `a2_refute.py`, `family_def_refute_edgebetw*.py`,
`family_def_refute_lenmax.py`).

---

## BOTTOM LINE

**No new differentiator beats co-threading once the component-size / node-degree confound is
removed.** The two candidates that looked new both fail the orientation+confound test:

- **edge betweenness (`edge_betw`) is node DEGREE in disguise**, not an independent bridge signal.
  Its headline AUC vs compara (0.683 raw, 0.939 in an earlier unstratified run) is a near-pure
  component-size artifact: corr(edge_betw, log component-size) = **−0.71**, spearman(edge_betw,
  deg_min) = **−0.78**. Residualize on size → AUC **0.531** (chance). Inside the big over-merge blob
  the raw bridge tail *does* track artifact edges (AUC overlap **0.78**, antisense **0.664**), but
  **degree alone does it better** (deg_min AUC overlap **0.871**, antisense **0.862**) and
  degree-residualized edge_betw collapses to **0.574 / 0.282** — i.e. it adds ≈0 beyond degree.
  And it **cannot be oriented by any trustworthy positive: 0/12 compara_pos and 0/5 panel_real lie
  in any component ≥50 edges.** Verdict: confounded, strictly dominated by degree. Not worth adding.

- **`len_max` is the most intrinsically bimodal feature** (BC=0.866, silhouette=0.799,
  dBIC 1→2 = +28,529, survives size-residualization BC_resid≈0.86) and is genuinely new and
  non-circular — **but it does not orient.** It is a short-cDNA vs long-cDNA axis, not real-vs-bridge:
  AUC vs compara only **0.317** (anti-oriented vs the artifact-negative scale), and **10/12 compara
  positives sit in the LOW-len mode.** Bimodal but label-irrelevant; would be a confound if used as a
  family filter.

**The strongest CLEAN lever remains co-threading** (`co_junc`/`contig`/`frac_mem`/`frac_ref`,
AUC ≈ 0.69–0.70 vs compara), which is **size-independent** (corr(co_junc, log-size) = −0.10) and
**orientable without circularity**. Same-strand (antisense=0) stays the orthogonal second lever.

---

## Ranked verdict table

| Feature | Intrinsic bimodality | AUC vs compara (raw) | After size/degree control | Orients? | Verdict | Confound / note |
|---|---|---|---|---|---|---|
| **co_junc** (co-thread) | BC 0.59, sil 0.71 | 0.695 | 0.70 (size corr −0.10, edge_betw corr ≈0) | yes | **KNOWN clean baseline** — best clean lever | none; the bar to beat |
| **contig / frac_mem / frac_ref** (co-thread) | BC 0.53–0.88 | 0.692–0.695 | size-independent | yes | **KNOWN** baseline | circular only for community-derived labels (not used) |
| **antisense / overlap** (strand/co-loc) | n/a (binary) | independent NEG | — | NEG by def | **KNOWN** 2nd lever | orthogonal to co-threading + topology |
| **edge_betw** (topology bridge) | BC 0.94 (size-driven) | 0.683 (0.939 unstratified) | **0.531** size-resid; **0.574/0.282** degree-resid in-blob | no (0/12 pos in blobs) | **CONFOUNDED — degree in disguise** | comp-size (−0.71), node degree (spearman −0.78); deg_min alone wins (0.87/0.86) |
| **common_nbr** (topology) | BC 0.62 | 0.914 → | **0.565** size-resid | weak | **CONFOUNDED** | component size (corr +0.64) |
| **deg_max / deg_min** (topology) | BC 0.45 | 0.932 / 0.926 → | 0.728 / 0.626 size-resid | weak | **CONFOUNDED** | component size (corr +0.81 / +0.75); also anti-orients in-blob |
| **jaccard_nbr** (topology) | BC 0.65 | 0.512 (perm-p 0.87) | 0.488 | no | **NOISE** | size-driven bimodality, no label signal |
| **len_max** | **BC 0.87, sil 0.80** (most bimodal) | 0.630 / **0.317** | survives resid but still doesn't orient | **no** | **NEW bimodal but NON-discriminative** | cDNA length (spearman +0.64 len_min / −0.67 len_ratio; −0.34 seq_id) |
| **len_min / njunc_*** | BC 0.49–0.60 | 0.50–0.71 | — | mixed | bimodal, weak/none | label-irrelevant bimodality |
| **log_dist** | BC 0.28, sil 0.76 | 0.787 (same_chrom only) | — | yes | KNOWN strand-genomic | only defined on same-chrom edges |

---

## (1) NEW genuine differentiators

**None.** No feature is simultaneously (label-free bimodal) AND (orients to compara) AND
(non-circular) AND (not a pure size/degree confound). The closest structural candidate
(edge_betw within blobs) fails the degree test; the closest bimodal candidate (len_max) fails
orientation.

## (2) Confirmed-but-known (co-threading / strand)

- **Co-threading** (`co_junc`, `contig`, `frac_mem`, `frac_ref`): AUC ≈ 0.69–0.70 vs compara,
  size-independent (|corr| ≤ 0.10), essentially uncorrelated with edge_betw — the clean lever.
- **Strand / co-location** (`antisense`, `overlap`): the independent negative axis; orthogonal to
  both co-threading and topology by construction.
- `log_dist` (closer = more real) confirms but is restricted to same-chromosome edges.

## (3) Confounded / circular (downgraded by verifiers)

- **edge_betw, common_nbr, deg_min, deg_max, jaccard_nbr** — all collapse after residualizing on
  **log component-size** (the dominant confound: small components mechanically force betweenness→1.0
  and degree→low). The component-size axis is itself partly the co-threading community-detection
  output, so "small comp = real" risks the documented circularity — another reason to distrust the
  raw topology AUCs.
- **len_max** — bimodal and non-circular but confounded with **cDNA length** (and indirectly with the
  known anti-discriminative `seq_id`, whose own AUC 0.659 ≥ len_max 0.630). Does not separate
  real-vs-bridge.

---

## Edge betweenness specifically — VERDICT WITH NUMBERS

**Independent alignment-free bridge signal, or just degree in disguise? → DEGREE in disguise.
Do NOT add it as a standalone differentiator.**

| Test | Number |
|---|---|
| corr(edge_betw, log component-size) | **−0.71** |
| spearman(edge_betw, deg_min / deg_max) | **−0.78 / −0.77** |
| AUC vs compara, raw | 0.683 (0.939 in unstratified run) |
| AUC vs compara, **residualized on log size** | **0.531** (chance) |
| Whole-group permutation null (positives span ~10 groups) | **p = 0.176** (the global edge-perm p=5e-5 ignores grouping → over-optimistic) |
| compara_pos / panel_real inside components ≥50 edges | **0 / 12 and 0 / 5** (cannot be oriented by any real positive in the blob) |
| Within-blob raw AUC: overlap / antisense | 0.78 / 0.664 |
| Within-blob **degree-residualized** AUC: overlap / antisense | **0.574 / 0.282** (adds ≈0; antisense reverses) |
| Within-blob **deg_min alone** AUC: overlap / antisense | **0.871 / 0.862** (degree dominates) |
| Within-blob grouped-block permutation p: overlap / antisense | **0.001 / 0.418** (only overlap even significant) |
| Within-blob single-feature GroupKFold (raw ranking transfer): overlap / antisense | 0.86 / 0.77 (transfers as a *ranking* — but it is the degree ranking; marginal contribution beyond degree ≈ +0.00) |

Reading: the *raw* within-blob bridge tail does correlate with co-location/antisense artifact edges
and transfers across held-out components — but only because those artifact edges are the low-degree
periphery of the giant component. Once degree is controlled, edge_betw's contribution vanishes
(overlap +0.00, antisense reverses below chance). It is non-circular (independent of co-threading:
corr(edge_betw, co_junc) ≈ −0.06, corr(edge_betw, frac_ref) ≈ −0.06) but **strictly dominated by
node degree**, which is cheaper to compute. **If anything is worth a probe inside blobs, it is
`deg_min` (peripheral low-degree edges = artifact candidates), not edge betweenness** — and even that
is unoriented by any real positive and should be treated as a heuristic, not a differentiator.

---

## How a (hypothetical) new differentiator would complement the known levers

If a clean topology signal had survived, its value would be **orthogonality + cost**: topology is
computed from the homology graph alone (no per-read realignment), whereas co-threading needs the
read-to-locus threading and strand needs orientation. The measured orthogonality holds
(edge_betw ⟂ co_junc, |corr| ≈ 0.06) — the problem is not independence, it is that the only
size-independent residual signal is **degree**, which is a confound rather than a real-vs-bridge
axis, and which has **no trustworthy positive inside the blob to orient it**. So the cheap
alignment-free graph features buy nothing over co-threading here.

---

## Recommended next step

Treat topology as a confound, not a lever: keep co-threading (`co_junc`/`contig`) + same-strand as
the family-edge gate, and — if a blob-internal heuristic is wanted — try **dropping low-`deg_min`
peripheral edges inside components ≥50** as a *candidate-flagging* step only, validated against
overlap/antisense (not compara, since 0 positives live there). To find a genuinely new differentiator,
the bottleneck is **labels inside the blobs**: get a handful of trustworthy positives (e.g. Compara
paralog pairs or curated tandem arrays) that actually fall in components ≥50 edges, otherwise every
in-blob signal is orientable only by artifact-negatives and indistinguishable from degree.

---

## ADDENDUM — the de-bridge test the AUC framing missed (cov_min, jaccard_nbr)

The verification above ranked features by "AUC vs the 12 Compara positives". That test is
*uninformative inside the over-merge blobs* (0/12 Compara positives and 0/5 panel reals land in any
component ≥50 edges), so it cannot adjudicate the actual goal — de-bridging the blobs. Re-ran the
PRACTICAL test instead: plug each cheap feature in as the community-detection edge weight (same
weighted-Louvain as co-threading) and measure max-component reduction + real-family preservation.

| weight (cost) | families | max comp | reals intact | anchor sep (238/140) |
|---|---:|---:|:--:|:--:|
| base (over-merged) | 1216 | 238 | — | — |
| frac_ref — co-threading (needs minimap2) | 1286 | 57 | 5/5 | 3/3, 2/2 |
| **jaccard_nbr — graph topology (FREE)** | 1263 | 58 | 5/5 | 1/1, 1/1 |
| **cov_min — reciprocal coverage (FREE, already computed)** | 1280 | **44** | 5/5 | 3/3, 2/2 |
| common_nbr — graph topology (free) | 1252 | 68 | 5/5 | — |

**`cov_min` de-bridges BETTER than co-threading (max comp 238→44 vs 57), keeps all 5 real families,
separates the unrelated anchors (3/3, 2/2) — at ZERO extra compute** (it is just `min(cov_a,cov_b)`
from the existing all-vs-all; `dna_homology()` already loads both, we only used the `max` for the
gate). `jaccard_nbr` (pure graph topology, no alignment) matches co-threading (58 vs 57).

Mechanism (Canzar-clean): real paralogs cover each other RECIPROCALLY (both cDNAs mostly align);
domain-sharers/bridges cover ASYMMETRICALLY (the shared domain is most of the small gene but little
of the large one → low cov_min). This is the continuous form of the airtight panel's 3rd axis
(recip-cov ≥ 0.30), now used as a de-bridge weight rather than a gate.

**Honest scope:** these three agree (Spearman +0.60–0.65) — they converge on the same dense-subgraph
structure, so they are NOT an independent validation axis (the synthesizer's point stands: nothing
adjudicates whether the LOC blob sub-communities are *real* families; that needs orthogonal
protein/Compara evidence). The win is practical, not epistemic: **cov_min is a cheaper, alignment-free
de-bridge weight that is at least as good as co-threading** — and ~35–40% of its rank variance is
independent of co-threading, so multiplying the two weights is worth trying.

**Updated recommendation:** switch the community-detection weight to `cov_min` (or `cov_min × frac_ref`)
— cheaper and a tighter de-bridge — keeping same-strand as the orthogonal filter. The remaining
bottleneck is unchanged: trustworthy positives inside the big components (protein-level) to confirm
the sub-communities are real families.
