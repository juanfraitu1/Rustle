# Graph-to-graph family alignment: does comparing whole isoform/exon-junction GRAPHS recover paralog families the single-transcript `aln_frac` misses?

**TL;DR — KEEP the shipped `aln_frac`; adopt graph-to-graph (GGA) only as the formal, VG-native
restatement of O1, NOT as a scoring replacement.** GGA is genuine graph structure (not `aln_frac`
relabeled — proven below), and it uniquely recovers 113 divergent-isoform paralog edges (CSNK1D/CSNK1E,
GNB1/GNB2, OGDH/OGDHL, …) that both `aln_frac` and an all-isoform proxy miss. But on the shipped
discrimination workhorse it is **worse** than `aln_frac` (held-out CV AUC 0.793 vs 0.835), its recall
lift is **indiscriminate** (real-edge rate 46.2% ≈ over-merge rate 42.7%), it recovers **0 of the 9**
missed diploid-oracle genes (re-scoring existing edges cannot CREATE the missing ones), and there is **no
matched-precision Pareto gain** (recall = 0.0 at precision ≥0.70 for every score). The value is in
**formalization elegance**, not in recall or FP. Fully consistent with the honest prior that everything
structural is redundant with `aln_frac`.

Artifacts (all absolute):
- `/mnt/c/Users/jfris/Desktop/Rustle/bench/graph_to_graph_family_aln.py` — build + analyze
- `/mnt/c/Users/jfris/Desktop/Rustle/bench/graph_to_graph_family_aln.tsv` — 5571 E_r candidate gene-pairs
- `/mnt/c/Users/jfris/Desktop/Rustle/bench/graph_to_graph_family_aln.json` — summary

Determinism: `PYTHONHASHSEED=0` re-exec, `minimap2 -t1`, deterministic `edlib`, sorted output, 5
component-level CV folds seed=0. `--analyze-only` re-runs byte-identical (md5 `7cf9a09…`).

---

## What was built (RNA-only; provenance identical to the shipped `aln_frac`)

Per gene, the **isoform repertoire** = all de-novo loci (`DN_*`) projecting to it (`family_er_pr`
max-overlap, `OV_FLOOR=0.20`); the shipped bridge representative(s) are always retained (so the graph
score is a strict superset of the single-rep `aln_frac`), then filled to `K_ISO=8` by `n_reads`
(dominant expression). Per gene the **exon-junction graph** = merged-union canonical exons (nodes) +
junctions used by any isoform (edges); node sequences are spliced genome bases from `GGO.fasta` — the
**same provenance and same `minimap2 -cx asm20` metric** as `ri_sharedlen_universal.tsv` (the shipped
`aln_frac`). Three comparable scores per pair:

| score | what it is |
|---|---|
| `aln_frac` | SHIPPED single-representative longest-shared-block / shorter (the workhorse). |
| `all_isoform_alnfrac` (proxy) | SAME metric, MAX over all isoform-A × isoform-B pairs. Isolates "does trying every isoform help." |
| `gga_score` | GRAPH-native: weighted mean of exon-homology (edlib infix id≥0.7, fwd/rc) + junction concordance + maximal shared colinear path. |

Substrate: **5571 candidate gene-pairs** — 416 TP (DNA-confirmed), 3382 truthbar-FP (truth-bar
divergent paralogs), 1773 genuine-FP (over-merges, no truth).

---

## 1. GGA-vs-`aln_frac` correlation — related, not a duplicate, not independent

| pair | Pearson | Spearman |
|---|---|---|
| `gga` vs `aln_frac` | **0.433** | **0.367** |
| `all_isoform` vs `aln_frac` | 0.150 | 0.029 |
| `gga` vs `all_isoform` | 0.171 | — |

`gga` is **moderately** correlated with `aln_frac` — not a relabel, not independent. The all-isoform
proxy is essentially **uncorrelated** with `aln_frac` because it **saturates**: mean `all_isoform`=0.883
and **83.8% of pairs hit ≥0.99** (some isoform pair always aligns fully).

---

## 2. Does full-GGA beat the all-isoform proxy? — TWO honest framings, and the framing is decisive

**On raw lift-COUNT: NO — structure "loses" 10:1.** Of the 2873 real paralog edges (TP + truthbar-FP)
with shipped `aln_frac`<0.24, the all-isoform proxy lifts **2523 (87.8%)** above 0.24 but `gga` lifts only
**1326 (46.2%)**; all-isoform-only recovers **1310**, `gga`-only just **113**. And adding junction+path
structure *beyond bare exon homology* net-**suppresses** recovery at the shipped line: it **helps 19,
hurts 523**. So on count, "try every isoform pair" wins and graph structure dilutes.

**But that win is indiscriminate saturation, and here is the decisive correction.** The proxy "wins" on
count only by aligning *some* short isoform fully for almost everything — inside its own saturated set it
is **random** (AUC 0.480). Graph STRUCTURE is exactly what **restores the discrimination the proxy
destroys**: within the 4668 pairs where `all_isoform`≥0.99,

| within saturated set (n=4668, 83.8%) | mean `gga` |
|---|---|
| TP (n=342) | **0.663** |
| truthbar-FP (n=2913) | 0.260 |
| genuine-FP (n=1413) | 0.322 |

AUC(`gga`, TP vs genuine-FP) within the saturated set = **0.812**, versus AUC(`all_isoform`) = **0.480**
(random) there. **A single-transcript / longest-block score cannot separate TP from FP where every
isoform pair already aligns fully — GGA does.** This is the proof that `gga` is *genuine graph structure,
NOT `aln_frac` relabeled.* So: structure DOES add discriminative signal over the proxy
(`gga_structure_beats_all_isoform_proxy_on_discrimination=True`); it merely loses the *lift-count* race
to indiscriminate saturation (`..._on_liftcount=False`). The honest bottom line is **"structure does not
BEAT `aln_frac`," not "structure is inert."**

---

## 3. How many real missed families does GGA recover, and at what FP cost?

### (a) The 9 diploid-oracle genes: **0/9 recovered** — recall is edge-GENERATION-bound, not score-bound

GGA re-scores *existing* E_r edges; it cannot CREATE them. **7/9 missed oracle genes have ZERO candidate
edges** — LOC101141440, LOC115930538, LOC115934629, LOC129534585, TNK2, UBE2Q2P16, ZNF425 — so there is
nothing to lift. The other 2 have edges only into over-merge blobs: **LOC129523567** (KRAB-ZNF, 20 edges,
8 real, best `gga`=0.483 but every real edge <0.24; its high-`gga` edges are all genuine-FP ZNF-blob
edges — lifting them re-merges the blob) and **LOC115932259** (1 genuine-FP edge only). Closing the
9-gene gap needs a different edge-generation step, not a better score.

### (b) Low-`aln_frac` real edges: structure-unique recoveries are genuine paralog families (named)

`gga` uniquely lifts **113 real paralog edges** where the all-isoform proxy also fails
(`gga_beyond_all_isoform`=113 vs `all_isoform_beyond_gga`=1310 — structure loses the count 10:1, but
these 113 are the *only* place graph structure earns recall). All have shipped `aln_frac`=0.0 (different
dominant isoforms make the single representative under-score) → `gga` 0.50–0.71. Named:

> **CSNK1D|CSNK1E** (0.71), **OGDH|OGDHL** (0.69), **CAMK2B|CAMK2G** (0.69), **GNB1|GNB2** (0.65),
> **ENO2|ENO3** (0.63), **PPP1CA|PPP1CC** (0.61), **RPS6KA1|RPS6KA2** (0.60), CBR1|CBR3, SLC9A6|SLC9A7,
> CPNE6|CPNE7, RARA|RARB, CKB|CKM, MAPKAPK2|MAPKAPK3, CLTC|CLTCL1, MAGOH|MAGOHB, ATP1A1|ATP1A3,
> CYTH2|CYTH4, MAP2K1|MAP2K2, PRPS1|PRPS2, DVL1|DVL3, ARRB1|ARRB2, MEIS1|MEIS2, PUM1|PUM2, LIMK1|LIMK2,
> EHD3|EHD4, GRK2|GRK3, RYR2|RYR3, ELMO1|ELMO2, RPL39|RPL39L, MYL6|MYL6B, …

(Only 2 of the 113 are DNA-confirmed TP — LOC129529768|LOC129529807, LOC129529767|LOC129529807; the rest
are truthbar-divergent-paralogs, the legitimate recall target defined by the task.) Where the all-isoform
proxy ALSO lifts, `gga` agrees on more genuine families — FMO1|FMO4, DPEP2|DPEP3, ACTN1|ACTN3,
SAMD9|SAMD9L, RAC2|RAC3, FRG1|LOC101142904 (`gga`=1.0).

### (c) FP cost: the lift is INDISCRIMINATE

| lever | real-edge lift-rate | genuine-FP lift-rate |
|---|---|---|
| `gga` ≥ 0.24 | **0.462** (1326/2873) | **0.427** (436/1022) |
| all-isoform ≥ 0.24 | 0.878 | 0.850 |
| structure-unique (`gga`≥0.24 ∧ proxy<0.24) | 0.0393 (113 real) | 0.0274 (28 FP) → 1.4× |

`gga` lifts real edges and genuine over-merges at **near-identical rates** (46% vs 43%). Named over-merges
it inflates: **MAGEB1|MAGEB3**, **ZNF70|ZNF79**. The residual over-merges (GSTM2/MAGE/KRAB-ZNF blobs) ARE
homologous — exactly the DNA-bound prior. Even the structure-unique lift is only 1.4× enriched for real.

---

## 4. Held-out component-CV AUC (5 folds, seed=0) — GGA is WORSE than `aln_frac` on the workhorse

**Workhorse task (TP-strict vs genuine-FP; comparable to the shipped ~0.83):**

| score | held-out CV AUC |
|---|---|
| **`aln_frac`** | **0.835 ± 0.046** |
| `gga` | 0.793 ± 0.014 |
| `exon_hom` | 0.777 ± 0.014 |
| `all_isoform` (proxy) | 0.516 ± 0.068 |

Graph structure is **slightly worse** than the single representative; no FP-side gain. (The proxy alone
collapses to random — again the saturation artifact.)

**Real (TP+truthbar) vs genuine-FP — everything near-random:** `gga` 0.535±0.065, `exon_hom` 0.527,
`all_isoform` 0.506, `aln_frac` **0.406** (below random — `aln_frac` actively inverts here because it
scores truthbar paralogs median 0.0 while over-merges mean 0.275). `gga` is the *best* of a near-random
field, but 0.535 is operationally useless — `exon_hom` (0.527) and `all_isoform` (0.506) are also ≈0.5,
so `gga` has no distinguishing FP-discrimination advantage.

---

## 5. NET Pareto — no operating point where GGA gives higher recall than `aln_frac`

Matched-precision recall (max recall at precision ≥ target, by thresholding each score):

| target precision | `aln_frac` | `gga` | `all_isoform` |
|---|---|---|---|
| ≥0.70 | 0.0008 | 0.0 | 0.0 |
| ≥0.80 | 0.0 | 0.0 | 0.0 |
| ≥0.90 | 0.0 | 0.0 | 0.0 |

**Recall = 0.0 at prec ≥0.70/0.80/0.90 for every score, in both framings.** The ceiling is not
`gga`-specific: the **maximum precision anywhere on the `aln_frac` PR curve is only 0.667** (base
precision 0.682; ~101 genuine over-merges have `aln_frac`≥0.99, as homologous as the 152 TP at ≥0.99). No
RNA-homology score reaches those precisions on this candidate set. Redundancy: ρ(`gga`,`aln_frac`)=0.367
(moderate). **No Pareto recall gain.**

---

## Honest verdict — is the value in RECALL, FP, or FORMALIZATION?

**FORMALIZATION only.** Graph-to-graph alignment does **not** make a measurable difference for recall or
for FP:

- **RECALL — no.** It recovers **0/9** oracle genes (7/9 have no edge to re-score — recall is
  edge-generation-bound, not score-bound); its recall lift is **indiscriminate** (real-rate 46% ≈
  over-merge-rate 43%); and there is **no matched-precision Pareto gain** over `aln_frac`. The 113
  structure-unique divergent-paralog recoveries (CSNK1D/CSNK1E, GNB1/GNB2, …) are real and confirm the
  hypothesis's kernel — different dominant isoforms make the single representative under-score — but they
  arrive with a matching flood of over-merges, so they are not usable recall.
- **FP — no.** On the workhorse, held-out CV AUC(`gga`)=0.793 < AUC(`aln_frac`)=0.835 (GGA is *worse*);
  on real-vs-over-merge every score is near-random (best 0.535). The residual over-merges are
  DNA-bound homologs no RNA score separates — exactly the stated prior.
- **FORMALIZATION — yes.** GGA is genuine graph structure (proven: within the saturated proxy set it
  ranks TP above genuine-FP at AUC 0.812 where the single-transcript proxy is random 0.480 — NOT
  `aln_frac` relabeled). Its defensible value is a clean **VG-native restatement of O1** (a gene is a
  variation graph; family = graph-homology component) that unifies with the O2 copy-assignment layer
  (copies = PSV-paths through one shared VG), not a scoring improvement.

**RECOMMENDATION: keep the shipped `aln_frac` as the operational O1 score. Adopt graph-to-graph only as
the formal definition** — the elegant VG-native statement of what `aln_frac` operationally approximates.
Fully consistent with the honest prior: everything structural we have tried (exon-containment, the L1
splice lever, the ~108-feature panel, and now GGA) is redundant with `aln_frac`; the only measurable
levers remain `aln_frac` + VG-multiplicity, and the FP residual is DNA-bound.
