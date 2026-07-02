# Per-block graph-structural FP detector — can a family's whole-block graph signature triage over-merges?

**Code:** `bench/block_graph_fp_detector.py` · **Outputs:** `bench/block_graph_fp_detector.tsv` (606 families × features + `fp_score` + `is_impure` + `oracle_fp`), `bench/block_graph_fp_detector.json`.
**Determinism:** `PYTHONHASHSEED=0`, seed=0 folds/Louvain/Fiedler; byte-identical across two runs (json md5 `b09545c9…`, tsv md5 `d9d8022e…`).
**Truth reproductions (exact vs shipped `bench/family_level_pr_current.json`):** E_p-impure = **80/606**; diploid-oracle FP = **6 distinct blocks**; oracle-mapped = **48**.

## The question

Prior work established two negatives at the **per-EDGE** level: (a) node-multiplicity catches the repeat-hub class (shipped as the fam17 repeat gate, `family_rna_refine.py`); (b) connectivity/topology (k-truss / k-core / betweenness / community / spectral, `graph_def_refine_sweep.py`) is **falsified** for the general over-merge because over-merges are DENSE not sparse (triangle-support AUC inverted 0.28, betweenness size-confounded). The new, untested angle here is not an edge cutter but a per-**BLOCK** confidence / triage score: does a family's *whole-graph* signature (multiplicity + connectivity + density + **cardinality**) predict that the block is an over-merge, so we can rank families for DNA review?

**Design (RNA-only detector; DNA = validation only).** Train on **E_p-impurity** (block spans ≥2 non-mega protein families; 80/606 positives; partly reference-derived but large enough to fit). Validate flags on the **diploid CN oracle** (6 DNA-confirmed FPs over 48 oracle-mapped families — too small to train, used only to validate). Features are derived from the induced denovo-family subgraph, `core_recip` edge weights, distinct-loci, and VG minimizer node-multiplicity (`vg_repeat_catalog.tsv`); the protein-family COUNT (the E_p label) is explicitly excluded from the feature set, and the DNA oracle enters only in precision@K.

## (1) The size confound is real and dominant

E_p-impurity rate rises monotonically with block size:

| block size (distinct loci) | 2 | 3–4 | 5–9 | ≥10 |
|---|---|---|---|---|
| impure rate | 8.5% (40/470) | 23.8% (20/84) | 35.1% (13/37) | 46.7% (7/15) |

So any raw-AUC signal is mostly "bigger blocks are more often impure." Two independent size controls follow: OLS residualization, and nested CV against a size-only baseline.

### Per-feature AUC(is_Ep_impure): raw vs LINEAR-size-resid vs QUADRATIC-size-resid

Residualized on a linear basis `[log n_members, log distinct_loci]` and, to expose linear-residualizer artifacts, on a **quadratic** basis `[logn, logdl, logn², logdl², logn·logdl]`. A ratio feature (e.g. `mem_per_locus = n_members/distinct_loci`, or `frac_in_giant`) is a *nonlinear* function of the size variables, so a linear residualizer cannot remove it and produces a **spurious** residual-AUC; the quadratic basis removes it and the apparent signal collapses to ~0.5.

| feature | raw | lin-resid | **quad-resid** | reading |
|---|---|---|---|---|
| n_genes / distinct_loci / n_members | 0.75 / 0.67 / 0.66 | 0.52 / 0.39 / 0.38 | 0.64 / 0.50 / 0.51 | = size, no residual graph signal |
| **med_mult** (VG node mult) | 0.67 | 0.641 | **0.655** | MULTIPLICITY survives BOTH bases |
| **max_mult / max_min_mult** | 0.69 / 0.68 | 0.571 / 0.558 | **0.566 / 0.550** | MULTIPLICITY survives |
| mem_per_locus | 0.48 | 0.665 | **0.329** | **ARTIFACT** — collapses/inverts under quad basis |
| frac_in_giant | 0.50 | 0.655 | **0.463** | **ARTIFACT** — collapses to ~0.5 under quad basis |
| edge_density | 0.33 | 0.361 | 0.367 | **INVERTED** (impure = denser) |
| fiedler / mincut / clustering | 0.35 / 0.43 / 0.57 | 0.350 / 0.352 / 0.424 | 0.376 / 0.391 / 0.395 | INVERTED / null |
| loci_per_seqclust (CARDINALITY) | 0.49 | 0.391 | 0.394 | does NOT predict E_p beyond size |

Two readings. **First**, the connectivity/density inversion confirms and extends the prior per-edge "over-merges are DENSE not sparse" finding to the whole-block level: impure blocks are denser, better-connected, higher-Fiedler, so connectivity is anti-predictive. **Second**, the two features that residualize *highest* under the linear basis (`mem_per_locus` 0.665, `frac_in_giant` 0.655) are raw-null features (raw AUC ~0.48–0.50) whose apparent residual signal is a **linear-residualizer artifact of a nonlinear size ratio** — the quadratic basis makes both vanish (0.329, 0.463). The surviving size-residual signals are all weak (quad-resid AUC 0.55–0.66) and the **only mechanistically interpretable size-independent lever is VG node MULTIPLICITY** (med/max/max_min ≈ 0.55–0.66), which is exactly the already-shipped fam17 repeat-hub feature — no new lever.

### CV logistic (5-fold random split, seed 0, pooled out-of-fold AUC)

| model | OOF AUC |
|---|---|
| size-only | 0.691 |
| graph-only | 0.762 |
| size + graph | 0.755 (Δ = **+0.064** over size) |

graph-only ≈ size+graph and the lift is small (+0.064) and small-n-noisy; it is carried by multiplicity.

## (2) Which FP CLASSES the signature flags

- **fam17 repeat-hub → YES, multiplicity flags it (held-out rank 7/606).** Block 17 (20 genes / 10 protein families) has `max_min_mult=26` (≥ the shipped gate's 20), `max_mult=497`, `frac_repeat_bridged=0.93`, fp_score 0.90. The per-EDGE gate did not dissolve it (the block stays connected through non-repeat edges) — the per-BLOCK multiplicity signature still catches it. But block 17 is **E_p-impure, NOT a DNA-confirmed FP** (`oracle_fp=0`); these repeat-bridged ZNF blobs may be real KRAB-ZNF families. Broadly, 74/80 E_p-impure blocks are not DNA-named FPs — the E_p proxy massively over-flags real divergent-paralog families.
- **MAGE cardinality → NO.** Block 549 (MAGEA1/4/12 + LOC129529978, the oracle multifam+oversize FP) is **E_p-PURE** (`n_protein_families=1`), so the E_p-trained score buries it at held-out **rank 563/606** (fp_score 0.014). And `loci_per_seqclust` is only **2.33** (not in the top-10, whose minimum is 4.38): **MAGE-A are DIVERGENT paralogs, not near-identical tandems**, so the "many near-identical loci" cardinality premise does not hold. The top-cardinality blocks are ZNF / tubulin / UBC families that are NOT confirmed FPs (1/10 oracle-FP). Cardinality does not flag MAGE.
- **GSTM2 domain-hub → no clean separation.** 5 GSTM2 blocks: FP block 9 (rank 75) and FP block 13 (rank 52) sit above the pure non-FP block 23 (rank 417) and above non-FP block 14 (rank 222) — but **FP and non-FP GSTM2 blocks interleave**: FP block 9 (rank 75) ≈ non-FP block 10 (rank 76). Its dense / high-clustering signature looks like a real dense family; the score gives no clean separation between FP and non-FP GSTM2.

## (3) Triage precision@K vs the DIPLOID ORACLE (held-out OOF score; 6 confirmed FPs)

| top-K | oracle-FP captured | recall | E_p-impure in top-K |
|---|---|---|---|
| 10 | **0/6** | 0.00 | 9/10 |
| 20 | 1/6 | 0.17 | 13/20 |
| 30 | 2/6 | 0.33 | 18/30 |
| 50 | **2/6** | **0.33** | 27/50 |

The score concentrates E_p-impure blocks (its training target) but does **NOT** concentrate the DNA-confirmed FPs. Held-out oracle-FP ranks are scattered — FOXO1-partner (fam 332) #19, SVBP (fam 37) #22, SEC22B (fam 13) #52, GSTM2 (fam 9) #75 — and the two single-protein-family over-merges are buried at the BOTTOM: **MAGE (fam 549) #563, LOC134758618 (fam 110) #597**. This is a structural blind spot: E_p supervision cannot see within-family cardinality/oversize FPs, not merely small-n.

## (4) Honest verdict

**A per-block graph signature is NOT a usable general FP triage / confidence score.**

1. **Size dominates.** Once you control for block size (log n_members + log distinct_loci), connectivity, density, and cardinality fall to null or invert; the two features that appear to survive the linear control (`mem_per_locus`, `frac_in_giant`) are nonlinear-size-ratio artifacts that vanish under a quadratic size basis. "Big/dense blocks are impure" is essentially the whole raw signal.
2. **The only size-independent lever is MULTIPLICITY — and it is already shipped.** Med/max VG-node multiplicity (quad-resid AUC 0.55–0.66) survives both residualizers and flags the repeat-hub class (fam17, ZNF blobs), reproducing the fam17 gate at block level. It adds nothing new and flags E_p-impure blocks that are largely NOT DNA-confirmed over-merges.
3. **It misses the DNA truth.** Only 2/6 oracle FPs in the top-50; the two single-family cardinality/oversize FPs (MAGE, LOC134758618) are actively ranked last because E_p-supervision is blind to within-family cardinality. GSTM2's domain-hub interleaves with real dense families. Cardinality (loci/seq-cluster) does not flag MAGE because MAGE-A paralogs are divergent, not near-identical.

**Caveats built in.** (i) Small oracle-n — 6 DNA FPs over 48 mapped families is a validation-only sample, never trained on. (ii) Size confound — controlled two ways (OLS residualization + size-only CV baseline). (iii) E_p proxy is partly reference-derived and over-flags (74/80 impure blocks are not DNA-named FPs). Net: the block graph confirms and extends the earlier "dense-not-sparse" falsification to the whole-block level; the single defensible signal (multiplicity) is the already-deployed repeat gate; and no block-graph score triages the DNA-confirmed over-merges — the confirmed FP classes (single-family cardinality, oversize, domain-hub) are precisely the ones the graph signature cannot see. **DNA/CN review remains necessary; a graph-only confidence score would be a mislabeled size-and-repeat heuristic.**
