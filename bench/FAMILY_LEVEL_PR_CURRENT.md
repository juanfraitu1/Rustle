# Gene-family-level (block-level) precision / recall of the CURRENT DEFAULT family definition

**Current default** = `bench/family_rna_refine.tsv` — the RNA-only recall-preserving refined oracle:
`core_recip >= 0.19 AND aln_frac >= 0.24` -> gamma-quasi-clique(gamma=0.20) -> allele-demote.
**607 multi-copy families** (committed `cc472a4`). The catalog is **loaded, not re-derived**; all
numbers come from the shipped eval machinery (`graph_def_refine_sweep.eval_partition`,
`rna_only_edge_oracle.oracle_residuals`) so they match the committed catalogs.

**Legacy** = shipped gamma-quasi-clique on the wider `core_recip >= 0.13` seed (gamma=0.20), **858 families**.
The legacy partition reproduces the committed `graph_def_refine_sweep.json` baseline **exactly**
(858 families / 115 impure / 50 oracle-mapped / 12 residual-FP / 52 recovered = OK cross-check).

Reproduce: `PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/family_level_pr_current.py`
Outputs (byte-identical across runs): `bench/family_level_pr_current.{json,tsv}`.
All three primary metrics are computed **at the family/block level** (they loop over blocks), NOT
at the edge level.

---

## The three truths, clearly labelled

### Truth 1 — PROTEIN E_p purity (block precision; conservative lower bound)
A family block is **PURE** iff its member genes span at most one non-mega protein family E_p;
**IMPURE = over-merge**. `Precision_Ep = pure_blocks / total_blocks`.

| | pure / total | impure | **P_Ep** |
|---|---|---|---|
| **CURRENT** | 528 / 607 | 79 | **0.870** |
| LEGACY | 743 / 858 | 115 | 0.866 |

Essentially flat (~0.87). **Caveat:** this is a *conservative lower bound*. E_p splits divergent
paralogs into separate PRFAMs, so genuinely-real divergent-paralog families (e.g. the APOBEC3
locus) get counted as "impurity." Real precision is higher than 0.87.

### Truth 2 — DNA-loose cDNA (`in_dna_loose`)
Gene-pair truth from `family_er_pr.tsv` (same values as `family_def_dna_pr_edges.tsv`), cDNA
id >= 0.90 & cov >= 0.30.

| metric (block-level) | CURRENT | LEGACY |
|---|---|---|
| **Component-recovery recall of real cDNA families** | **182 / 187 = 0.973** | 195 / 200 = 0.975 |
| Block over-merge precision (blocks with NO genuine over-merge pair) | 0.424 (249/432 carry a genuine pair) | 0.502 (350/703) |
| Pair-projection recall (shipped, *pair-level not block-level*) | 0.922 | 0.978 |
| Pair-projection precision (shipped, *pair-level not block-level*) | 0.472 | 0.502 |

The headline is the **block-level component-recovery recall of real cDNA families (~0.97)**: the
current default still fully co-members essentially every real cDNA family. **Caveat:** the block
over-merge precision (0.42) is a *pessimistic lower bound* — `in_dna_loose` requires cDNA
id >= 0.90, so divergent paralogs (real families) are labeled "genuine over-merge." The
pair-projection numbers are reported for continuity only and are explicitly *not* block-level
(they are dominated by megablob / clique-expansion pairs).

### Truth 3 — INDEPENDENT DIPLOID DNA ORACLE (the gold, assembly-independent)
The only assembly-independent scalar (`asm_hapCN` = maternal + paternal diploid copy number).
`P_oracle = 1 - (multifam + oversize + allele) / oracle_mapped`;
`R_oracle = oracle_multicopy_genes_recovered / oracle_multicopy_genes`.

| metric | CURRENT | LEGACY |
|---|---|---|
| oracle-mapped families | 48 | 50 |
| FP (allele / oversize / multifam) | **0 / 3 / 4** | 2 / 4 / 6 |
| flag-instances / distinct FP blocks | 7 / **6** | 12 / 11 |
| **P_oracle (task formula 1 - flags/mapped)** | **1 - 7/48 = 0.854** | 1 - 12/50 = 0.76 |
| **P_oracle (dedup 1 - distinct/mapped)** | **1 - 6/48 = 0.875** | 1 - 11/50 = 0.78 |
| **R_oracle (multi-copy genes)** | **48 / 57 = 0.842** | 50 / 57 = 0.877 |

**Current is materially cleaner on the gold oracle: residual FP 12 -> 7 (precision 0.76 -> 0.85),
allele FPs 2 -> 0.** The cost is a small recall drop (0.88 -> 0.84): current over-splits exactly
**two** multi-copy oracle genes vs legacy — `LOC115932259` and `LOC129523567`. The latter was the
*worst* legacy oversize FP (dl=54, 13.5x diploid), which the current default over-splits away —
i.e. it trades one recall point to kill the single largest over-merge.

---

## Enumerated residual false-positive set we STILL get (current default, named)

### (a) DNA-CONFIRMED over-merges — 6 distinct blocks (7 flag-instances)

**Multifam (block merges >= 2 distinct diploid-oracle genes):**

| family_id | dominant_gene | distinct_loci | offending oracle genes (+ context) |
|---|---|---|---|
| **9** | LOC115930164 | 15 | **GSTM2 + LOC115930164 + LOC115930576** (+9 GSTM-cluster LOCs) — GSTM2-hub #1 |
| **13** | SEC22B | 3 | **GSTM2 + LOC101129940** (+ SEC22B) — GSTM2-hub #2 (the other GSTM2 splinter) |
| **333** | LOC115933254 | 10 | **FOXO1 + LOC115933254** (+ ANKRD18A/B, LOC101142783/147259, LOC129524732) |
| **550** | MAGEA12 | 7 | **LOC129529978 + LOC129529986** inside the MAGEA1/MAGEA4/MAGEA12 block — MAGE-A sub-cluster (also oversize) |

**Oversize (distinct RNA loci > 1.5x diploid DNA CN = collapsed / inflated copies):**

| family_id | dominant_gene | dl | dipCN | ratio | offending genes |
|---|---|---|---|---|---|
| **550** | MAGEA12 | 7 | 2 | 3.5x | LOC129529978 + LOC129529986 (same block as multifam #4) |
| **37** | SVBP | 11 | 4 | 2.75x | **MPHOSPH8** (+ FBXO17, SVBP, 7 LOCs) |
| **109** | LOC134758618 | 10 | 6 | 1.67x | **LOC134758618 + NLGN1** |

### (b) Residual allele-as-copy: **0**
The demote step (`balanced_frac >= 0.90 & copy_like <= 0.10`) removed both shipped allele FPs
(DHRSX, LOC129530050); neither recurs.

### (c) E_p-impure blocks NOT reachable by the DNA oracle: **73 of 79**
The diploid oracle covers only 48/607 families, so it cannot adjudicate these. Salient members:

- **fam 17** (dl=28, **16 protein families**): DCP1B / GIMAP4 / KLHL33 / LIPA / ... / UCHL1 — a
  repeat-bridge mega-hub, the clearest genuine over-merge.
- fam 78 (dl=16, 3 fams: LOC129532709 / LOC134759133 / NDUFA8), fam 68 (KLHL3 + SLC38A6),
  fam 91 (SERPINA3 + LOCs), fam 120 (ARL17A + LRRC37A3 cluster), fam 340 (RGPD8 cluster + LOC),
  fam 421 (APOL2/3/4 + PP2D1 + RAB5A).
- **Probably REAL divergent-paralog families mis-flagged** because E_p splits them into separate
  PRFAMs (do NOT count against biological precision): **fam 603 APOBEC3C + APOBEC3D/F** (the
  APOBEC3 locus), **fam 164 AKR1C3 + IL20RB**, **fam 42 RHD + SDHD**.

---

## One-line takeaway
On the gold diploid oracle the current default improves precision **0.76 -> 0.85** (residual FP
**12 -> 7**, allele FP **2 -> 0**) at a small recall cost **0.88 -> 0.84**; protein E_p purity is
flat at **~0.87** (conservative lower bound); cDNA family recovery is **~0.97**. The only
DNA-confirmed FPs left are **6 blocks** — the two GSTM2-hub splinters (fam 9, fam 13), the
FOXO1 + LOC115933254 merge (fam 333), the MAGE-A sub-cluster LOC129529978/986 (fam 550), and two
collapsed-copy oversizes MPHOSPH8 (fam 37) and LOC134758618 + NLGN1 (fam 109) — plus 73 E_p-impure
blocks the DNA oracle cannot reach (led by the fam-17 16-family repeat-bridge hub).

---

**Determinism:** `PYTHONHASHSEED=0`; `family_level_pr_current.{json,tsv}` byte-identical across
repeated runs; legacy partition reproduces the committed `graph_def_refine_sweep.json` exactly.

**Files (absolute):**
- `/mnt/c/Users/jfris/Desktop/Rustle/bench/family_level_pr_current.py`
- `/mnt/c/Users/jfris/Desktop/Rustle/bench/family_level_pr_current.tsv` (per-family flags:
  `ep_impure`, `fp_multifam`, `fp_oversize`, `fp_allele`, `oracle_genes`, `diploid_CN`)
- `/mnt/c/Users/jfris/Desktop/Rustle/bench/family_level_pr_current.json` (full P/R table + named FP roster)
