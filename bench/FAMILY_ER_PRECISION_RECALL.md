# Precision / Recall of the reframed RNA family definition $E_r$

**Scope.** Re-measures precision/recall for the *reframed* RNA multi-copy family definition
$E_r$ (transcript-homology component, $\gamma$-quasi-clique refined; `bench/family_definition_formal.md`),
replacing the old read-conflict $E_c$ measurement. Two independent truths are used:

- **cDNA-90% truth** (`bench/family_def_dna_pr_edges.tsv`, `in_dna_loose` = all-vs-all cDNA id $\ge$90% cov $\ge$30%).
  STRICT — structurally **under-counts divergent paralogs** below the 90% nt bar.
- **protein truth $E_p$** (`bench/protein_families_refined.tsv` + raw `aln.m8`): reciprocal whole-protein
  mmseqs homology, **no identity floor**, reaches the 40–70% protein-id twilight zone
  (globins / CEACAM / KRAB-ZNF). Used ONLY for the precision-side FP split — never mixed into recall.

Reproduce: `PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/family_er_pr.py`
→ `bench/family_er_pr.{tsv,json}`. **Byte-deterministic** across launches (hash seed pinned; the refined
block *set* is seed-invariant, the witness partition is pinned to `seed=0`).

---

## 1. Headline numbers

| quantity | value |
|---|---|
| refined $E_r$ gene-pairs | **10,755** |
| TP vs cDNA-90% truth | 450 (raw cDNA precision only 0.042 — see §3, the truth is incomplete, this is **not** the honest metric) |
| **genuine precision** (truth-bar-misses counted correct) | **0.613** |
| **GENUINE over-merge (headline FP)** | **4,163 edges = 38.7% of refined edges** |
| — **family-level over-merge (lead with THIS)** | **115 / 858 refined blocks $E_p$-impure = 13.4%** |
| real divergent paralog the truth misses (**not** an error) | 6,142 edges (59.6% of all FP; 99.3% a direct reciprocal-protein edge) |
| **absolute recall** (expressed & cross-map) | $E_r$ **0.150** vs $E_c$ **0.552** (E_r LOWER — footprint-bound) |
| **reachable recall** (both loci in footprint) | $E_r$ **0.949** vs $E_c$ **0.686** (**+0.26**, criterion-isolating, conditional) |

**One-line takeaway.** $E_r$ is a **recall/precision trade-off vs the old $E_c$-refined**, not a domination:
it recovers the divergent low-identity paralogy $E_c$ structurally drops (id 0.90–0.95 reachable recall
**0.951 vs 0.328**) at the cost of genuine precision **0.613** (vs the old $E_c$-refined 0.94). The looser
0.13-core homology grouping buys the divergent paralogs and pays in over-merges.

---

## 2. Precision: raw vs refined $E_r$, vs both truths

| | gene-pairs | cDNA-90% precision | $E_p$ precision (incl. divergent paralogs) | genuine precision |
|---|---|---|---|---|
| **RAW** $E_r$ (single-linkage, core_recip$\ge$0.13) | 170,411 | 0.003 | 0.243 | 0.243 |
| **REFINED** $E_r$ ($\gamma$-quasi-clique, $\ge$2 loci) | 10,755 | 0.042 | 0.614 | **0.613** |

Refinement $R$ is load-bearing: it collapses the 170k raw pairs (dominated by the DNFAM0 728-locus
repeat/ZNF blob, the RNA analog of the 22q11 DNA blob) to 10.8k, lifting genuine precision 0.243 → 0.613.

The **cDNA-90% precision (0.042) is deliberately NOT the headline** and no "effective $\approx$1.00" claim
is made: the old "effective 1.00" was **circular** (an id$\ge$0.80 pre-filter forced residual FP to have
id$\ge$0.80). Here precision is stated against an *independent* protein arbiter, and the FP are split, not
laundered.

---

## 3. The honest FP split (the deliverable)

For every refined-$E_r$ FP edge (`in_er & !in_dna_loose`), classify with the **independent** protein truth:

```
protein_homologous := reciprocal-whole-protein mmseqs edge (evalue<=1e-5, both qcov,tcov>=0.5)
                      OR same non-mega PRFAM
real_divergent_paralog := protein_homologous AND cdna_id>=0.50 AND min_cov>=0.60
genuine_over_merge     := NOT real_divergent_paralog     <-- the real false positive
```

| bucket | count | what it is |
|---|---|---|
| **REAL divergent paralog** (truth-bar miss) | **6,142** (59.6% of FP) | protein-homologous below the cDNA-90% bar — the cDNA truth is *incomplete*, not $E_r$ wrong |
| — tier-1 direct reciprocal protein edge | 6,099 (99.3%) | prot pident median 47%, whole-protein min-cov median 0.84 (all $\ge$0.50) |
| — tier-2 same non-mega PRFAM | 43 | |
| **GENUINE over-merge** (HEADLINE FP) | **4,163** (38.7% of edges) | real errors — 4,198-of-4,209 (pre-floor) have **no protein edge at all** |
| — both protein-coding | 2,363 | |
| — of which mega-blob co-family (KRAB-ZNF/GPCR, weak) | 759 (ZNF–ZNF 420) | **kept in the over-merge count** (conservative; arguably weak real paralogs) |

**Why the split is honest / non-circular.** The $E_p$ arbiter is independent of both the cDNA-90% truth and
the $E_r$ RNA construction. 99.3% of the 6,142 truth-bar FP carry a *direct* reciprocal whole-protein edge
with substantial coverage (recomputed straight from `aln.m8`: 0/6,099 fail the reciprocal $\ge$0.5-coverage
test). Named exemplars, all reciprocal substantial-coverage homology the cDNA-90% bar drops:
ENO2–ENO3 (83.5% / cov 0.98), OGDH–OGDHL (77% / 0.94), GALNT14–GALNT16 (57% / 0.90), ZNF234–ZNF45 (55.5% / 0.96).
Conversely the genuine-over-merge bucket is genuinely errors: AMY2A–ZNF91, AMY2A–ZNF141, RPL14–ZNF669 have
**empty `aln.m8` rows** (unrelated-gene bridges surviving the loose 0.13 POA homology). The guard bites
CONSERVATIVELY toward over-merge (the 11 protein-homologous pairs with whole-protein min-cov$<$0.6 stay in the
error bucket; the 759 mega-blob ZNF pairs stay in the error bucket) — the opposite of laundering.

### Headline the family-level rate, not the clique-inflated edge rate
The 38.7% edge rate is **clique-inflated**: a single over-merged block contributes $\binom{k}{2}$ FP edges
(the top-6 residual blobs alone are 57% of all genuine-FP edges). The honest family-over-merge rate is
**block-level: 115 / 858 refined blocks are $E_p$-impure = 13.4%** (mix $\ge$2 distinct non-mega protein
families). *Lead with 13.4%.*

**Disclosed floor, never the headline:** the "clean both-protein-coding, non-blob genuine over-merge = 1,604"
is a *lower bound* only. The over-merge count is **4,163** (edge) / **13.4%** (block); 1,604 must never be
quoted as the over-merge count.

### Min-overlap floor (precision-only hardening, this revision)
The DN-locus→gene projection now applies a **min-overlap floor** (`OV_FLOOR=0.20`): a DN→gene map is kept
iff the overlap covers $\ge$20% of the DN span **OR** $\ge$20% of the gene span. This removes incidental
$<$20%-clip attributions (e.g. a DN locus grazing ALB by 0.18 of the gene) that otherwise attribute an
unrelated gene into a block and inflate the over-merge count. Effect: dropped 26 spurious locus maps,
removing **46 genuine-over-merge edges (4,209 → 4,163)**; it touched **0 TP and 0 truth-bar edges** — a
surgical precision-only fix. The pre-floor **4,209** is retained in the JSON
(`overlap_floor.genuine_over_merge_prefloor`) as a **conservative upper bound**.

---

## 4. Recall (identity-stratified) and the gain over $E_c$

Recall is measured for **both** $E_r$ and $E_c$ against the **same** cDNA-90% truth (`in_dna_loose`).
$E_p$ is *not* used here (no cross-contamination).

| subset | n | $E_r$ recall | $E_c$ recall | note |
|---|---|---|---|---|
| **(a) expressed & cross-map** (task-specified) | 2,116 | **0.150** | **0.552** | ABSOLUTE — $E_r$ **lower** |
| (b) expressed only | 4,257 | 0.106 | 0.275 | ABSOLUTE — $E_r$ lower |
| **(c) reachable** (both loci in $E_r$ footprint) | 334 | **0.949** | **0.686** | CONDITIONAL — isolates the criterion, +0.26 |

**Identity strata (the smoking gun for the divergent-paralog recovery):**

| id bin | (a) $E_r$ / $E_c$ | (c) reachable $E_r$ / $E_c$ |
|---|---|---|
| 0.90–0.95 | 0.190 / 0.367 | **0.951 / 0.328** |
| 0.95–0.99 | 0.128 / 0.405 | 0.923 / 0.608 |
| $\ge$0.99 | 0.160 / 0.776 | 0.972 / 0.909 |

**How to read the gain honestly (guardrails).**
- On **both absolute** subsets $E_r$ recall is **BELOW** $E_c$ — because $E_r$'s de-novo footprint reaches only
  **13.7%** of the 5,517 cDNA-truth genes (refined-block representation; the raw-locus footprint is 14.2%).
  Cite the **13.7%** refined figure when discussing refined recall.
- The **+0.26** appears **only** under the reachable conditioning, which drops **1,782** truth pairs. Of those,
  $E_c$ loses **940 correctly-grouped TP** (pairs it groups that lie *outside* $E_r$'s footprint) while $E_r$
  loses **0** by construction (`in_er` requires both genes co-blocked ⟹ both represented ⟹ reachable). So (c)
  measures the **family-definition CRITERION** (homology $E_r$ vs read-conflict $E_c$) on a convenient
  denominator — it is **never** absolute recall superiority.
- **Rule for the thesis body:** the +0.26 gain / the "0.95 vs 0.33" low-id win must ALWAYS travel with the
  attached caveat *"conditional on both loci in $E_r$'s 13.7% de-novo footprint; on the full expressed &
  cross-mapping universe $E_r$ absolute recall is lower (0.150 vs 0.552)."* The reachable n is small
  (n=334; id 0.90–0.95 bin n=61) — note the modest n.

**"Recover ~30% of what $E_c$ dropped" — corroborated within the footprint.** On the reachable subset $E_c$
misses 31.4%; $E_r$ closes 26.3 of those points (84% of $E_c$'s reachable misses), and the recovered mass sits
at low identity (id 0.90–0.95: $E_c$ 0.328 → $E_r$ 0.951) — exactly the divergent paralogy $E_c$ structurally
drops, consistent with the $E_r$ reframe's stated purpose.

---

## 5. Consistency with the old $E_c$ note (no contradiction)

The old "$E_c$-refined: precision 0.94 / 79 residual FP" was a **different, tighter edge set** (the 1,367-edge
homology-filtered read-conflict catalog), and its "effective $\approx$1.00" was circular (id$\ge$0.80 prefilter).
This note does **not** re-make that argument; it reports $E_r$'s own genuine precision **0.613** against an
independent protein arbiter. The same truth file underlies both eras (verified: $E_c \cap$ dna_loose = 1,822
= the old "RAW 2,829 TP1822"; $E_c$ recall over all dna_loose = 1,822/17,410 = 0.105 = old "R.10"). The new
$E_r$ artifact is *consistent with*, not contradictory to, the old note.

---

## 6. Determinism

`family_er_pr.{tsv,json}` are **byte-identical across launches** (`PYTHONHASHSEED=0` pinned in-script;
`seed=0` witness partition). Arithmetic closes exactly: TP 450 + truth-bar 6,142 + genuine 4,163 = 10,755
refined edges. Multi-seed wobble (documented, expected — the $\gamma$-quasi-clique certificate is
seed-invariant, the witness partition is not): block counts 858/857/859, pair counts 10,755/10,661/11,240
across seeds 0/1/2.
