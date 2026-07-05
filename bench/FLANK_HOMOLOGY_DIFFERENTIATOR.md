# Flanking-homology differentiator — a new orthogonal axis with modest lift

**Question.** The settled FP boundary says the residual single-domain over-merges are
separable only by DNA copy-number (`FP_EXCLUSION_DISCRIMINATORS.md`).  Those five axes
all look *inside* the transcript (nt/protein/TE/topology/architecture).  Does sequence
*outside* the transcript — upstream/downstream flanking regions — carry an orthogonal
signal?  Whole-locus duplications should share regulatory/context sequence; domain-shares
should not.

**Answer.** **Yes, but modestly.** Flank homology is a genuine new differentiator: it
adds **+0.016 to +0.024** OOF AUC over the shipped `core_recip + aln_frac` baseline on
the 5,571 refined-`E_r` direct edges, and it reaches **AUC ~0.76–0.80** on the
irreducible residual roster.  However, it is **not a clean hard gate**: coverage is sparse
(~17–21% of edges have a detectable flank hit), many real paralogs have no detectable
flank homology, and confirmed over-merges such as the **GSTM2 hub also show high flank
homology**.  The raw signal is partly repeat-driven, but masking repeats does not remove
the additive lift.

Code: `bench/flank_homology_differentiator.py`  
Outputs: `bench/flank_homology_edge*.tsv`, `.json`  
Deterministic: `PYTHONHASHSEED=0`; byte-identical re-runs for a fixed flank size.

---

## 1. Method

**Population.** The same 5,571 refined-`E_r` direct edges used in `rna_levers_explore`
(joined via `ri_sharedlen_universal.tsv` → representative `dn_a`/`dn_b` loci).

**Flank extraction.** For each de-novo locus, extract `F` bp upstream and `F` bp
downstream of the transcript, oriented 5'→3' relative to transcription (reverse-complement
for `-` strand genes).  Contig boundaries are clipped.  Tested `F = 1000, 2000, 5000`.

**Homology.** Build separate FASTA files for all upstream and all downstream flanks,
run `minimap2 -x asm20 -c -X -t 4` all-vs-all for each orientation, and parse PAF.  For
each edge we keep the best hit by matching bases and compute:

- `up_aln_frac` / `down_aln_frac` — reciprocal alignment span / flank length
- `up_id` / `down_id` — alignment identity
- `max_flank_aln_frac`, `mean_flank_aln_frac`, `min_flank_aln_frac`, `best_flank_id`

Missing alignments are encoded as `NA` in the TSV and imputed as **0** for the additive
logistic evaluation (no detectable flank homology = no homology).

**Evaluation.** Out-of-fold logistic AUC over the same component folds stored in
`rna_levers_explore.tsv`.  Baseline = `core_recip + aln_frac`.  Labels =
`in_dna_loose` (DNA/cDNA real edge, validation-only).  Residual-crack test on the
irreducible roster (MAGE/GSTM2/etc.).

---

## 2. Results

### 2.1 Lift over the shipped `core_recip + aln_frac` baseline

| flank | edges with hit | `max_flank_aln_frac` AUC | OOF AUC `core+aln+max_flank` | lift |
|---|---:|---:|---:|---:|
| 1 kb | 864 / 5571 (15.5%) | 0.633 | 0.9199 | **+0.0160** |
| 2 kb | 939 / 5571 (16.9%) | 0.663 | 0.9216 | **+0.0178** |
| 5 kb | 1075 / 5571 (19.3%) | 0.699 | 0.9272 | **+0.0234** |
| 2 kb masked | 695 / 5571 (12.5%) | 0.571 | 0.9224 | **+0.0185** |

Baseline `core_recip + aln_frac` OOF AUC = **0.9038**.  Lift is **robust across flank
sizes** and slightly strengthens with longer flanks.  Soft-masking repeats lowers the
univariate AUC (0.663 → 0.571) but leaves the additive lift essentially unchanged,
suggesting the logistic value is not purely repeat-driven.

For comparison, the best RNA lever in `rna_levers_explore.md` (splice/structural) gave
**+0.0089** lift.  Flank homology roughly **doubles** that.

### 2.2 Residual-crack test (irreducible MAGE/GSTM2 roster)

| flank | `max_flank_aln_frac` AUC | n_pos (real) | n_neg (over-merge) |
|---|---:|---:|---:|
| 1 kb | **0.803** | 17 | 22 |
| 2 kb | 0.759 | 18 | 24 |
| 5 kb | 0.765 | 18 | 26 |
| 2 kb masked | 0.705 | 18 | 16 |

On the hard residual, flank homology is a noticeably stronger separator than on the full
population.  This is the regime where `aln_frac` alone cannot split real paralogs from
domain-shares.

### 2.3 Distribution by class (2 kb unmasked)

| `max_flank_aln_frac` | TP (real) | truthbar | genuine (over-merge) |
|---|---:|---:|---:|
| NA | 96 | 3259 | 1277 |
| 0.9–1.0 | **235** | 28 | 280 |
| 0.7–0.9 | 37 | 21 | 53 |
| 0.5–0.7 | 22 | 24 | 52 |
| 0.3–0.5 | 18 | 34 | 56 |
| 0.1–0.3 | 8 | 16 | 55 |

When flank homology is **high (>0.9)**, it is enriched in real edges (235/320 = 73% of
edges with a hit are real).  But **280 genuine over-merge edges also have high flank
homology**, so a simple threshold would misclassify many of them.

### 2.4 Mechanistic spot-checks

- **GSTM2 hub.** Many GSTM2-over-merge edges show `max_flank_aln_frac` ≥ 0.9 (e.g.
  GSTM2–LOC101129940 = 0.94, GSTM2–LOC109023809 = 1.0).  The GSTM2 copies share not
  only the GST domain but also adjacent sequence, consistent with a real recent
  segmental duplication.  This is why flank homology does **not** dissolve the GSTM2
  residual — the copies look like whole-locus paralogs at the sequence level.

- **MAGE sub-cluster.** Real MAGE edges (e.g. LOC129529976–LOC129529986 = 1.0) and some
  truthbar/cardinality edges show high flank homology, while other MAGEA edges have NA.
  The signal is mixed, matching the sub-cluster's complex expansion history.

- **LOC134758618 / FOXO1-partner.** Some confirmed over-merge edges have low/moderate
  flank homology (0.35–0.50), the regime where the feature is genuinely informative.

---

## 3. Interpretation and caveats

**What flank homology does well.** It directly tests the *whole-locus duplication*
mechanism: real paralogs share context, domain-shares do not.  It adds the largest lift
of any RNA lever tried so far and is strongest exactly on the high-`aln_frac` residual
where the shipped oracle struggles.

**Why it is not a clean gate.**

1. **Sparse coverage.** Only ~15–20% of edges have a detectable flank hit.  For the other
   ~80% the feature is effectively 0, so it cannot directly cut most over-merges.
2. **Confirmed FPs can be whole-locus duplications at the DNA level.** The GSTM2 hub and
   some MAGEA edges have homologous flanks, so the feature flags them as *real* even
   though the DNA/protein truth splits them into distinct families.  This is not a bug —
   it means the RNA+flank signal and the DNA-family boundary disagree on those cases.
3. **Repeat contribution.** Univariate AUC drops ~0.09 when repeats are soft-masked,
   indicating a substantial fraction of raw flank homology is repeat/TE content.  The
   existing repeat gates handle repeat bridges *inside* transcripts; flank repeats are a
   separate, weaker axis.
4. **No natural threshold.** Unlike the antisense/reciprocal-overlap rule (a geometric
   impossibility at ≥50%), flank homology is a continuous similarity score with overlap
   between real and false-positive distributions.

**Relation to the settled FP boundary.** This is a **sixth axis**, genuinely orthogonal
to the five in `FP_EXCLUSION_DISCRIMINATORS.md`, but it does **not** break the
DNA-CN-bound residual.  It clips part of it and adds probabilistic signal, much like the
splice/structural lever in `RNA_LEVERS_EXPLORE.md`.

---

## 4. Conclusion and recommendation

Flank homology is a **real, reproducible, orthogonal differentiator** that improves the
edge classifier.  The cleanest operational value is as a **soft feature** in an edge-level
model, not as a hard threshold gate.  Adding it to the oracle would give ~+0.02 AUC and
help triage the high-`aln_frac` residual.

**Recommendation.** If the goal is to improve the edge oracle, flank homology is worth
incorporating as a continuous feature (best form: `max_flank_aln_frac` over ~2–5 kb).  If
the goal is a new *hard gate* that removes FPs without collateral, flank homology is
**not** it — the GSTM2/MAGE residuals remain.

**Next experiments to consider:**
- Cross-species ortholog-pattern consistency (Ensembl Compara) — do all family members map
  to one human gene family?
- Genome-wide intron-phase / splice-site conservation (panel was promising in
  `family_def_splicing_signals_a4.py`).
- Combining flank homology with the spectral sub-structure features from
  `TP_FP_DIFFERENTIATOR_PANEL.md` in a single edge model.
