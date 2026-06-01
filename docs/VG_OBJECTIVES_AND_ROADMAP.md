# Multi-Copy Gene Family (VG/EM) — Objectives & Roadmap

The thesis's core objective is NOT to reproduce StringTie (that's the StringTie-faithful *floor*, proven
via `-G`: 100% Sn / >99% Pr — see `STRINGTIE_PARITY_CONSOLIDATED.md`). It is to do what StringTie
*cannot*: assemble the copies and isoforms of multi-copy gene families. The four research objectives
(documented in `docs/ALGORITHMS.md §5–12`) and the concrete avenues to attain them:

## Obj 1 — Cross-family information borrowing during assembly
Let the assembler borrow exon structure / junction evidence / coverage from paralog loci so each copy
is reinforced by family-wide evidence (a read supporting exon E in copy A is evidence for E in copy B).
- **State:** assembly is still PER-BUNDLE; the EM only reweights reads. The architectural enabler
  (a family-bundle merged graph) is not built.
- **Hard blocker (demonstrated, GOLGA6L7):** when paralogs are <~85kb apart, the bundle-builder MERGES
  them into ONE bundle, so VG (which acts at bundle boundaries) can't separate them → L7_2 gets 0
  predictions. This is the most common tandem-array structure, so it gates the whole objective.
- **Avenues:**
  1. **Intra-bundle paralog splitter (highest value):** detect paralog sub-structure within a merged
     bundle (repeat-unit periodicity / k-mer self-similarity / coverage-shoulder) and split or
     sub-cluster, then assemble per-paralog. Target example: GOLGA6L7 (recover L7_2).
  2. **Family-merged graph + shared-evidence flow:** build the conceptual family VG (collapse
     near-identical exon classes), run flow with shared-node capacity from all copies. Target: any
     family where a copy is starved of unique evidence.

## Obj 2 — Discovery of copies absent from the reference
If a copy is missing from the reference, its reads have no locus to map to and are lost.
- **State:** scaffold exists (`--vg-scan-novel-loci`, `discover_novel_copies`, HMM rescue in
  `vg_hmm/rescue.rs`); UNVALIDATED on a known absent-copy case.
- **Avenues:**
  1. **Validate on a constructed absent-copy case:** mask one copy out of the reference, confirm the
     k-mer/HMM profile scan recovers its reads + assembles it. Measurable, decisive.
  2. **Unmapped-read rescue:** pull reads that map poorly everywhere but match a family's k-mer profile;
     assemble a de-novo copy. Target: a family with a known reference-absent paralog.
  3. **Distinguish novel copy vs artifact:** the k-mer-Jaccard / POA-identity filters (just fixed for
     inversions) are the discriminator — extend to score candidate novel loci.

## Obj 3 — Novel isoforms and structural variants per copy
Per-copy splice graph built from post-EM reads, so unusual combinations unique to one paralog survive.
- **State:** the EM-reweighting prerequisite now works (incl. INVERTED families after the 2026-05-31
  canonical-k-mer fix — DAZ3 recovered as 5 isoforms). Standard `path_extract` assembles per-copy from
  the reweighted coverage. **Structural variants (inversions, large indels) NOT handled at graph level.**
- **Avenues:**
  1. **Validate per-copy isoform completeness:** for DAZ3 (5 isoforms recovered), confirm against the
     DAZ repeat-array structure how many isoforms are real vs over-enumerated (needs a curated truth).
  2. **Structural-variant handling:** parse soft-clips / supplementary alignments per copy to capture
     inversions / tandem expansions / TE-insertion exons. Research-grade; the DAZ inversion is the
     entry example.

## Obj 4 — Accurate read-to-copy assignment in ambiguous cases
A hierarchy of signals: (1) junction pattern, (2) diagnostic SNPs/INDELs (`--vg-snp`), (3) exon-length,
(4) copy-level EM prior; propagate fractional weight for genuinely ambiguous reads.
- **State:** signals 1–3 implemented and now **validated end-to-end at 100% decisive accuracy** on the
  synthetic fixture (`test_data/synthetic_family/`, 28 multi-mappers, oracle Obj 4). Reaching that took
  fixing three separable blockers found via the oracle (2026-05-31):
  1. **Site pollution (rustle bug, fixed):** `build_exon_fingerprints` compared *same-copy* read
     fragments (ragged terminal exons give one copy several variable-length versions) as if they were
     distinct copies → 132 spurious "diagnostic" sites + 29× multiplicity → assignment collapsed to
     noise. Fixed by collapsing `per_copy_sequences` to one representative (longest) per *distinct*
     CopyId and requiring ≥2 distinct copies (132→12 real sites).
  2. **Genome not threaded without `--vg-snp`:** the fingerprint-EM (`--vg-no-hmm`) needs
     `read.mismatches`, which are only extracted when the genome reaches BAM ingestion — gated on
     `--vg-snp`. Without it every placement scored as a perfect self-match → 0.5/0.5. The oracle scorer
     now passes `--vg-no-hmm --vg-snp` (the design-spec invocation).
  3. **Over-strict gate:** the gap-gate (`10.0` log-units) is calibrated for reads covering *thousands*
     of sites; on sparse coverage it silenced correct decisions (a 2-net-site / ~80:1 signal is only
     ~4.4 log-units). Now **evidence-adaptive**: `eff_gap = min(gap_threshold, per_site_gap·n_sites)`
     (`RUSTLE_VG_EM_SCORE_GAP_PER_SITE`, default 0.25) — unchanged in the abundant-site regime,
     proportional when sparse. Isolated to `--vg-no-hmm`; the default HMM-EM (DAZ) is untouched.
  **Missing: a per-transcript `copy_assignment_confidence` GTF attribute** (copy_id is set; the EM weight
  gap is computed but not yet emitted as an attribute).
- **Honest finding (this session):** on DAZ the *iterative* EM converges trivially (delta=0) — the
  discriminating signal is the pileup-depth PRIOR (53/47), not iteration. So in the
  identical-compatibility regime it's "pileup-weighted assignment," and the value is honest uncertainty
  quantification, not iterative refinement.
- **Avenues:**
  1. **Emit `copy_assignment_confidence`** per transcript (from the EM weight gap / n_sites_covered) —
     turns "copy recovered" into "copy recovered with quantified confidence." Cheap, high thesis value.
  2. **Calibrate the per-copy expression split** end-to-end: confirm DAZ1≈53% / DAZ3≈47% materializes in
     the per-copy coverage, not just the prior. Validates the EM is doing the right thing.

## Cross-cutting enabler — a MEASURABLE multi-copy benchmark (the missing oracle)
The parity work is measurable (gffcompare vs StringTie). The multi-copy work's whole point is going
BEYOND StringTie, so there is no off-the-shelf gold standard — which means "are we attaining the
objectives" is currently anecdotal (DAZ, GOLGA6L7). **The foundational avenue is a multi-copy benchmark
+ oracle:**
- **Synthetic ground-truth** (extend `test_data/synthetic_family/`): planted copies with known
  per-copy isoforms + known read-to-copy truth → measure recovery Sn/Pr AND assignment accuracy per
  objective.
- **Curated real families** (DAZ1/DAZ3, GOLGA6L7, AMY, a reference-absent case): a small panel with
  expected outcomes, run on every change. This is the multi-copy analogue of the parity oracles — the
  user-endorsed "oracle for the flow-related steps" direction, applied to the thesis objectives.

## Prioritized recommendation
1. **Obj 4 confidence + Obj 3 calibration (builds directly on today's DAZ EM fix; cheap, high thesis
   value):** emit `copy_assignment_confidence`, validate the 53/47 split.
2. **The multi-copy benchmark/oracle (foundational — makes all objectives measurable).**
3. **Obj 1 intra-bundle splitter (unlocks the largest family class — GOLGA6L7).**
4. **Obj 2 absent-copy validation (the reference-bias differentiator).**

Parity track (on the side, user-endorsed): lever #3 (checktrf rescue-vs-drop, the 2B-shape), the
flow-step oracle, and the `RUSTLE_CSR_FOLD` default-flip decision.
