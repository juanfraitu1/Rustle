# Scope: a divergent-paralog detection track (beyond the contiguous-core horizon)

_Status: SCOPE / not started. Author handoff doc. Builds on `compara_validation.py`, `family_detection_validation.py`, `denovo_family_pipeline.md`._

## 1. Motivation — the gap, with evidence

The shipped family detector defines a family by POA **contiguous-core coverage ≥ `t_core` (0.13)** — the longest *contiguous exact* shared run / min(len). This is tuned for **recent tandem duplicates** and is correctly conservative.

It is blind to **divergent paralogs**. Measured on DSFAM45 (`diag_prefilter_homology`, NC_073234.2:48446103-49179309):

| pair | shared 18-mers | longest contiguous run | core cov (poasta == LCS) |
|---|---|---|---|
| [3]×[4] (≈2 kb copies) | **374** | ~76 bp | **0.036** |
| unrelated pair | 2–9 | — | 0.003–0.01 |

374 shared k-mers vs 2–9 for unrelated ⇒ unmistakable common descent; but the longest contiguous exact run is ~76 bp, so contiguous-core (poasta **and** the alignment-free LCS agree — not an artifact) puts it at 0.036, far under 0.13. **These are real paralogs past the recent-duplicate horizon.**

**Why it's worth detecting (the through-line):** divergent paralogs carry MANY fixed differences → by the identifiability theorem they are the *easy, high-PSV* regime for copy-assignment (the advisor's #2 interest). The families we currently can't *detect* are the ones we could most cleanly *assign*. Detection is the gate; assignment is the payoff.

## 2. The central question (sharp, falsifiable)

> Among de-novo locus pairs the contiguous-core criterion **rejects**, can a sequence feature — validated against an **independent, LOC-covering** paralog truth set — separate TRUE divergent paralogs from FALSE domain/repeat-sharers, at a **calibrated FDR**?

Reject criterion to beat: a clean, ROC-calibrated edge predictor with reported FDR (no arbitrary threshold) that recovers DSFAM45-type families without admitting domain-sharers.

## 3. The make-or-break — the truth set — is largely already in place

The original sin (`family_detection_validation.py`) was **circularity**: a minimizer-Jaccard clustering scored against a minimap2-built universe — both reward shared subsequence. `compara_validation.py` de-circularized it with Ensembl Compara (gene-tree + species-tree, never sees our alignments) but hit a **coverage wall: 154/195 genes are LOC-only → unmappable to Ensembl → only 12 checkable pairs.** Too small to calibrate a threshold.

**The fix we can now afford — OrthoFinder on `proteins.fa`:**
- `/home/juanfra/winloci_scratch/proteins.fa` = the RefSeq gorilla proteome. RefSeq annotates **LOC genes with protein products (XP_)**, so OrthoFinder clusters them too → the truth set **covers the LOC paralogs Compara could not** (the majority, and exactly the interesting copies).
- Self-contained: no Ensembl-ID mapping, no gorGor↔T2T liftover.
- Independent of a DNA k-mer feature (DIAMOND protein all-vs-all → MCL → gene trees).
- The advisor-endorsed modality ("need protein-level (OrthoFinder) for a real claim").

**Truth plan:** OrthoFinder orthogroups = PRIMARY truth (LOC-covering, hundreds–thousands of paralog pairs). Compara (already cached, `compara_paralog_relation.json`) = independent CROSS-CHECK on the mappable named subset (guards against an OrthoFinder-specific bias). Granularity caveat carries over: orthogroups are COARSER than copy-families ⇒ score **precision** (our pairs that the truth also relates), not recall.

## 4. Feature candidates (test both; lead with the independent one)

| feature | independence vs protein truth | sensitivity to divergence | cost |
|---|---|---|---|
| **DNA k-mer containment / Jaccard** on the de-novo nt seq | FULL (cross-modality) | moderate | cheap (have the k-mers already) |
| **Protein-ORF identity** (translate longest ORF → AA → identity/k-mer) | PARTIAL (both protein; mitigate: pairwise vs OrthoFinder's tree-clustering) | high | needs ORF-finding on de-novo tx |
| shared-k-mer **count**, length-normalized | FULL | moderate (the DSFAM45 signal: 374 vs ~5) | cheap |

Lead claim = **DNA/k-mer feature vs OrthoFinder truth** (clean cross-modality). Report protein-ORF as the sensitive-but-caveated upper bound. If the cheap DNA feature already separates on independent truth, that is the strongest, simplest result.

## 5. Phased plan — gates, kill-criteria, effort

**Phase 0 — Truth set + labeled pairs (≈1–2 d). THE gate.**
- Run OrthoFinder (or a direct DIAMOND-all-vs-all + MCL if OrthoFinder is too heavy) on `proteins.fa` → orthogroups.
- Map each de-novo locus → overlapping RefSeq gene (`GGO_genomic.gff`) → protein → orthogroup.
- Build labeled pair set: POSITIVE = same orthogroup; NEGATIVE = length/GC-matched different-orthogroup pairs.
- **GATE:** ≥ a few hundred labeled paralog pairs INCLUDING LOC genes. If not (truth too sparse even with OrthoFinder) → fall back to Compara-only, smaller claim, and flag the reduced scope. **KILL if no independent truth materialises.**

**Phase 1 — Is there signal? (≈1 d).**
- Feature AUC vs the independent label, restricted to the contiguous-core-REJECT population (cov < `t_core`) — i.e. exactly the pairs the shipped detector misses.
- **KILL if best AUC < ~0.7 on independent truth** → divergent paralogs are sequence-indistinguishable from domain-sharers (an honest, publishable negative; stop).
- Decide the feature.

**Phase 2 — Calibration / principled edge (≈1 d).**
- ROC-calibrate the operating threshold on the truth; report FDR at the operating point. NO arbitrary cutoff.
- Advisor-principled option: emit the feature as a weighted edge and let the EXISTING weighted-modularity Louvain (`family_split.rs`) decompose — the threshold folds into the objective.

**Phase 3 — Integration as a SEPARATE, gated track (≈2–3 d).**
- New `divergent_family_detect` path running ALONGSIDE (not replacing) contiguous-core: same candidate_pairs prefilter, but the REJECTED-by-core pairs get the validated divergent edge criterion; reuse `decompose_families`.
- Opt-in (`RUSTLE_DIVERGENT_PARALOGS`), default-OFF → byte-identical baseline (the project's standard for every lever). TDD, adversarial review.
- Apply to DSFAM45 → does it recover the family at the reported FDR?

**Phase 4 — The payoff: does it help copy-assignment? (≈1–2 d).**
- On the recovered divergent families, run the existing copy-assignment + measure resolvable-PSV / identifiability. Hypothesis: divergent paralogs are HIGH-identifiability (many PSVs) → clean assignments. This closes the loop to the advisor's interest and is the real success bar.

## 6. Principled framing (for Canzar)

- **No arbitrary threshold:** the edge criterion is ROC-calibrated on an independent biological truth with a reported FDR; ideally folded into the existing weighted-modularity objective so the cutoff is emergent, not hand-set.
- **Precision, not recall:** truth is coarser (superfamilies) ⇒ we score "of the pairs we join, what fraction the independent truth also relates," and intentionally do NOT penalise correctly splitting ancient superfamilies into copy-clusters (the established framing).
- **Through-line:** detection feeds the identifiability story — divergent paralogs are the high-PSV, cleanly-assignable regime. The track expands the set of families we can *resolve*, which is the thesis core.

## 7. Risks / open decisions for the user

1. **Truth source priority.** OrthoFinder-primary (more setup, breaks the 12-pair ceiling, the real claim) vs Compara-only (already cached, but only 12 mappable pairs → underpowered). *Recommend OrthoFinder-primary.*
2. **Feature priority.** DNA-independent (cleanest claim) vs protein-ORF (most sensitive, needs de-novo ORF-finding, semi-circular caveat). *Recommend lead DNA, report protein.*
3. **Success bar.** Phase 3 (recover DSFAM45-type families at a calibrated FDR) — or push to Phase 4 (show the recovered families assign cleanly), which is the advisor-facing payoff. *Recommend committing to Phase 4 as the bar.*
4. **Honest-negative acceptance.** The kill-criteria can end in "no separable feature" — a real result. Agree up front that a clean negative is an acceptable outcome (it bounds the contiguous-core decision as provably optimal on sequence alone).

## 8. First concrete step

Phase 0, sub-step 1: get orthogroups from `proteins.fa` (OrthoFinder, or DIAMOND+MCL) and map de-novo loci → genes → orthogroups via `GGO_genomic.gff`. That single artifact (a labeled paralog-pair table covering LOC genes) is the gate the entire track hinges on; everything after is cheap by comparison.
