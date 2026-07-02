# The RNA-only refined edge oracle for the multi-copy family definition

**Thesis-critical division of labor.** The RNA multi-copy family definition is
DERIVED FROM RNA. DNA/protein truth is used only to (i) LEARN the operating point
(threshold calibration), (ii) VALIDATE the result (scoring), and (iii) name the
IRREDUCIBLE residual. **No DNA, protein, or genome-annotation value ever gates an
edge at inference.** This document reports the RNA-only feature set, the learned
interpretable rule, the held-out performance vs the shipped `core_recip>=0.13`
threshold, the RNA-only refined family catalog, the over-merge FP it removes vs the
honest irreducible residual (DNA-validated on BOTH tails), and the allele demotions.

Everything is deterministic: `PYTHONHASHSEED=0`, fixed CV folds by seed, sorted
writes. Re-runs are byte-identical (`bench/rna_only_edge_oracle.tsv`/`.json`); the
feature table md5 is preserved (`3806b13c…`).

---

## 1. RNA-only feature set (no DNA at inference)

The edge unit is a gene-pair asserted to sit in the same shipped refined-E_r family
(`edge_bridge_mask.tsv`, N=10755). The **inference** decision reads exactly two
RNA-derived features; a third read-derived feature was tested and dropped:

| feature | provenance (RNA) | role |
|---|---|---|
| **core_recip** | whole-transcript reciprocal homology weight = POA reciprocal contiguous-core over `denovo_transcripts.fa` (de-novo **read-assembled** transcript consensus, annotation-free) | universal-coverage **precision gate** |
| **aln_frac** | longest shared spliced-exon block / shorter transcript; RNA splice STRUCTURE from the de-novo assembly (`ri_sharedlen_universal.tsv`) | **workhorse** separator |
| ~~mm read-degree~~ | genome-wide read multi-mapping degree of the bridge locus | tested, **dropped** (held-out ΔF1 +0.0018) |

**Allele DEMOTE** (family-level, orthogonal to the edge rule) reads only two
read-derived quantities: `balanced_frac`, `copy_like` (minor-allele-fraction from
the read-consensus O1 signal, `a1_read_consensus_o1.tsv`).

**Hard guards (asserted at runtime, all pass):**
- Inference feature set `{core_recip, aln_frac, mm_log}` is `⊆` the RNA whitelist and
  **disjoint** from `{in_dna_loose, in_ep, ep_tier, cls, abl_bridge_mask, allele_*}`.
  `no_dna_in_decision=True`, `edge_decision_features=[core_recip, aln_frac]`,
  `demote_features=[balanced_frac, copy_like]`.
- The genome/library soft-mask (`abl_bridge_mask`) is emitted ONLY as a flagged
  ablation column; it is never in a decision.
- `in_dna_loose` (reference-cDNA), `in_ep` (whole-protein), and the diploid CN oracle
  enter only as the training TARGET `y` (calibration) and via `eval_partition` to
  SCORE (validation).

### aln_frac provenance — the one honest exposure, now bounded

`aln_frac` uses RNA splice structure but fetches the exon **bases** from the genome
assembly. A strict RNA-primary examiner is entitled to ask whether the dominant
separator hinges on genome bases. It does not. We recomputed the **identical**
longest-shared-block fraction from the de-novo **read-consensus** POA transcripts
(`denovo_transcripts.fa`, the same source as `core_recip`; NO genome bases,
`ri_sharedlen_readcons.tsv`) and compared:

| provenance check (n=5571 core-present edges) | value |
|---|---|
| Spearman ρ(genome-base aln_frac, read-consensus aln_frac) | **0.9941** |
| KEEP/CUT gate agreement at `aln_frac>=0.72` | **0.9975** (5557/5571) |
| held-out AUC vs `in_dna_loose` — genome bases | 0.9151 |
| held-out AUC vs `in_dna_loose` — read consensus | **0.9149** |

The purity-max catalog `rna_only_ruleF1_readcons` (aln from read consensus) yields
**420 families vs 421** for the genome-base rule, with identical residual counts.
Independently, the pure read-POA feature `core_recip` alone confirms the separation
at held-out AUC **0.815**. The separation is a property of the transcript, not the
base source.

---

## 2. The learned interpretable RNA-only rule

Component-level 5-fold CV: whole connected components (669 families, union-find over
all 10755 edges) are held out — no family spans folds, no gene appears in >1 fold
(asserted). Population = the **5580 core-present direct edges** the gamma operator
thresholds (416 TP / 1774 genuine over-merge / 3390 truthbar divergent-paralog).

**Deployable rule (F1-optimal; fold-stable, median == final):**
> **KEEP an E_r edge iff `core_recip >= 0.26` AND `aln_frac >= 0.72`, else CUT.**

**Recall-preserving rule (max train precision s.t. recall ≥ 0.90):**
> **KEEP iff `core_recip >= 0.19` AND `aln_frac >= 0.24`.**

### Held-out performance (out-of-fold; each edge scored by a rule trained on folds NOT containing its family)

| | AUC | precision | recall | F1 | TP kept | genuine over-merge cut | truthbar cut |
|---|---|---|---|---|---|---|---|
| aln_frac (alone) | **0.915** | — | — | 0.439 | — | — | — |
| core_recip (alone) | 0.815 | — | — | 0.365 | — | — | — |
| logistic core+aln+mm | 0.904 | — | — | — | — | — | — |
| **rule F1** (0.26 / 0.72) | — | **0.416** | 0.577 | **0.483** | 240/416 | **1587/1774 = 89.5%** | 3240/3390 = 95.6% |
| **rule recall-preserving** (0.19 / 0.24) | — | 0.242 | **0.873** | 0.378 | 363/416 | — (catalog ge-cut 0.669) | — |

**Honest framing (not the x5.6 headline).** The shipped `core_recip>=0.13` is the
edge-**creation** threshold: on this population it keeps 100% of edges, so its
"precision" 0.0746 is the base rate, NOT an operating point — the x5.6 factor is
generous and is footnoted, not headlined. The real, defensible contribution is:
- **aln_frac is the workhorse** (held-out AUC 0.815 core-only → **0.915** core+aln);
  a smooth logistic does **not** beat it out-of-fold (0.904 < 0.915), and the
  read-degree veto adds nothing (ΔF1 +0.0018, dropped).
- **core_recip is a precision GATE** on top of aln (held-out F1 0.365 core-only,
  0.439 aln-only → **0.483** core+aln; precision 0.352 → 0.416 at matched recall).
- Standardized logistic coefficients: `aln_frac +1.23` (dominant), `core_recip +0.30`,
  `mm_log −0.10` — consistent with the AND-rule.

---

## 3. RNA-only refined family catalog

Apply the rule at inference (RNA features only) → re-threshold the E_r DN-graph →
rebuild families with the UNCHANGED shipped operator (connected components →
γ-quasi-clique refinement γ=0.20 seed=0 → ≥2-distinct-loci predicate) → DEMOTE
allele-as-copy. Scored with the sweep's own `eval_partition` + the
assembly-independent diploid CN oracle. `tp_ret`=DNA-real edge retention;
`ge_gcut`=genuine over-merge edge-cut (the honest per-edge number);
`tb_ret`=truthbar/divergent-paralog retention; `recov`=validated oracle genes
retained multi-copy; `undersz`=over-SPLIT tail (diploid CN > 1.5× RNA loci).

| catalog | n_fam | tp_ret | ge_gcut | tb_ret | orac | recov | allele | oversz | undersz | multifam |
|---|---|---|---|---|---|---|---|---|---|---|
| gamma_shipped_0.20 | 858 | 0.978 | 0.316 | 0.150 | 50 | 52 | 2 | 4 | 30 | 6 |
| **rna_only_ruleF1** | 421 | 0.674 | **0.899** | 0.007 | 45 | 47 | 2 | **2** | 30 | **4** |
| **rna_only_ruleF1_demoted** | 418 | 0.674 | 0.900 | 0.007 | 42 | 44 | **0** | 2 | 30 | 4 |
| rna_only_heldout_ruleF1 (non-circ) | 418 | 0.654 | 0.893 | 0.007 | 44 | 46 | 2 | 2 | 30 | 4 |
| **rna_only_recall_preserving** | 610 | **0.922** | 0.669 | 0.019 | 51 | 51 | 2 | 3 | 30 | 4 |
| rna_only_ruleF1_readcons (purity-max) | 420 | 0.674 | 0.899 | 0.006 | 45 | 47 | 2 | 2 | 30 | 4 |

The held-out catalog ≈ the deploy catalog (oversize 2, multifam 4 in both) → the
removal is **not** an artifact of fitting thresholds to these families.

---

## 4. Over-merge FP removed (RNA-only) vs the irreducible residual

Of the 12 named DNA-confirmed residual FP from `GRAPH_DEF_REFINE_SWEEP.md`
(2 allele + 4 oversize + 6 multifam), the RNA-only definition **removes 6** — all
out-of-fold:

- **Edge rule shatters** `LOC129523567` (54-locus ZNF/repeat-bridge blob, aln_frac=0),
  the `MAGEA9` 18-locus cluster, `LOC134758618`, `LOC101142904+LOC129526550`, and one
  `GSTM2`-hub instance (3→2 spanned oracle genes).
- **DEMOTE removes both alleles** (allele_as_copy 2 → 0).

**Irreducible residual (survives RNA-only + demote; DNA-validated as the honest RNA limit):**
- `MPHOSPH8` (dl=11 vs diploid CN 4 = 2.75×): same-gene locus multiplicity, no
  cross-gene bridge for the edge rule to cut.
- `LOC129529978+LOC129529986` (dl=7 vs 2 = 3.5×): residual MAGE sub-cluster left after
  the MAGEA9 bridge was severed.
- `GSTM2+LOC115930164+LOC115930576` (dl=14), `GSTM2+LOC101129940` (dl=3): GSTM paralogs
  aligning near-FULL length (core 0.79–0.97, aln_frac ≈ 1.0) — DNA says separate/fewer
  loci, but RNA shares a long transcript body. **Only the DNA-CN/protein arbiter
  separates them; this is the honest RNA limit.**
- `FOXO1+LOC115933254` (dl=10, aln_frac 0.88).

### Symmetric over-split bound (both tails, assembly-independent)

Precision is bought by aggressively cutting divergent-paralog edges, which narrows
the E_r homology definition. We bound that cost on the assembly-independent diploid
CN oracle (`asm_hapCN`), the ONLY assembly-independent quantity:

- **Held-out CUT edges touching a diploid MULTI-COPY locus (dip≥2 = a potential real
  copy split off):** TP **28/176**, truthbar **35/3240**, genuine **61/1587** (oracle
  coverage is small because only 60 genes are in the curated oracle; every covered cut
  edge is dip≥2).
- **Validated oracle-gene recovery:** shipped **52 → RNA-only 44** (Δ **−8**): eight
  DNA-validated families are over-split out of the multi-copy catalog by the F1 rule.
  The recall-preserving rule keeps **51** (Δ −1).
- **Undersize families (dip > 1.5× RNA loci):** **30 → 30** — the rule does NOT worsen
  the pre-existing K=0 / collapsed-copy floor (RNA genuinely under-counts high-CN
  segdups like GSTM2 dipCN=38 → ~14 RNA loci; that floor exists in the shipped catalog
  and is unchanged).

Net: the F1 rule removes 6/12 over-merge FP at an honest cost of ~176 divergent-paralog
edges and 8 validated families; the recall-preserving rule removes fewer over-merges
(oversize 4→3, multifam 6→4) at near-full recall (tp_ret 0.922, recovery 51).

---

## 5. Allele demotions (RNA read signal only)

`balanced_frac >= 0.90 AND copy_like <= 0.10` (minor-allele fraction ~0.5 = diploid
het, not ~1/K = copy) fired on exactly **3 families**, all with identical pure-het
signatures (`balanced_frac=1.000, copy_like=0.000`):

| gene | med_maf | verdict |
|---|---|---|
| **DHRSX** | 0.464 | DNA-confirmed (asm_hapCN=1) — removed |
| **LOC129530050** | 0.410 | DNA-confirmed (asm_hapCN=1) — removed |
| LOC115932259 | ~0.5 | novel, DNA-unvalidated — disclosed |

`oracle_allele_as_copy` 2 → 0. The edge rule cannot touch these (both loci project
to ONE gene → within-gene edge, always kept), so the allele-as-copy FP genuinely
needs this separate RNA read signal.

---

## 6. Recommended RNA-only refined definition

**Recommendation: adopt the two-feature RNA-only gate as the refinement layer, at the
RECALL-PRESERVING operating point by default, and expose the F1 point as a
high-precision variant. Do NOT replace the shipped `core_recip>=0.13` wholesale.**

Rationale, advisor-aligned (interpretable rule + honest residual):
1. `core_recip>=0.13` is the edge-**creation** threshold and stays — it defines which
   edges exist. The RNA-only rule is a **refinement gate** layered on top, not a
   replacement of the creation threshold.
2. The E_r homology family is *defined* to include divergent paralogs. The F1 point
   (0.26/0.72) trades 42% of DNA-real TP edges and 95.6% of truthbar edges for
   over-merge precision — a genuine **re-narrowing** of the definition toward
   high-identity recent paralogs. That is a legitimate high-precision variant but
   should not be the default definition.
3. The **recall-preserving** point (0.19/0.24) preserves divergent-paralog recall
   (tp_ret 0.922, recovery 51 ≈ shipped 52, undersize unchanged) while still removing
   the worst over-merges (genuine edge-cut 0.316→0.669, oversize 4→3, multifam 6→4).
   This is the honest default: cut the dense repeat/TE-bridge blobs, keep the paralogs.
4. Layer the DEMOTE step (balanced_frac/copy_like) unconditionally — it is orthogonal,
   removes the allele-as-copy FP, and costs nothing on the copy families.

DNA/protein remain **calibration + validation only**: they picked the thresholds and
scored the result; they never gate an edge. The purity-max read-consensus catalog and
the pure-read `core_recip` feature confirm the result does not hinge on genome bases.

---

## 7. Reproducibility

Deterministic (`PYTHONHASHSEED=0`, seed=0 folds and community detection, sorted
writes); byte-identical across re-runs.

- `bench/rna_only_edge_oracle.py` — feature + LEARN (component-CV) + APPLY stages,
  provenance concordance, symmetric over-split bound, guards.
- `bench/ri_build_sharedlen_readcons.py` → `bench/ri_sharedlen_readcons.tsv` —
  read-consensus aln_frac (no genome bases) for the provenance check.
- `bench/ri_build_sharedlen_universal.py` → `bench/ri_sharedlen_universal.tsv` —
  leakage-free universal aln_frac (in_ep-independent).
- `bench/rna_only_edge_features.tsv` (md5 `3806b13c…`, unchanged) — per-edge feature
  table with DNA/protein LABELS clearly separated from RNA features.
- `bench/rna_only_edge_oracle.tsv` / `.json` — the 6 catalog rows + full per-example
  tracking, provenance, over-split bound, and guards.
