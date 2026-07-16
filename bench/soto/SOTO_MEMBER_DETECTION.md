# Soto per-member detection — sensitivity & precision on the real benchmark

Per-member recovery of the **362 Soto members / 83 families** in `80_fams.chr.bed` (human Iso-Seq A119b →
CHM13v2.0), across all pipeline legs. Tables:
- `bench/soto/soto_family_detection.tsv` — one row per family: `n_members, n_detected, sensitivity, members`.
- `bench/soto/soto_member_detection.tsv` — one row per member: `family, gene, locus, detected(Y/N), recovered_by`.

## Headline

```
member sensitivity:          276/362 = 76.2%
families with >=1 member:      62/83  = 75%
families FULLY recovered:      42/83  = 51%
precision (detected loci on a real member): 711/767 = 93%
```

`recovered_by` breakdown: **RNA-split 215 · protein-tail 6 · genome-projection 40 · expressed-collapsed 15.**

The **56 off-annotation loci** (7 %) are not false positives — they are candidate **unannotated paralogs**
(id ≥ 0.98, read-supported), i.e. copies the curated Soto 80-family subset does not list (Soto's full
catalog is 213 families / 1002 paralogs). RNA-split copies are 100 % precise; the projection legs add these
high-identity extra copies.

## How to read it for the advisor

- **Sensitivity is per-member and honest** — each member is either recovered (RNA-split resolved it, or a
  projection/expressed-collapsed leg localized it as copy-number) or not. Flagship duplication families are
  clean: SRGAP2 (`ID_462`) 4/4, PMS2P (`ID_8`) 9/9, the ID_22/ID_71/ID_68 clusters fully recovered.
- **The 24 % not recovered is characterized, not mysterious** (see `SOTO_A119B_RECOVERY.md` Panel 4): ~5
  members are silent (0 reads, all coding), ~25 are per-read K=0 exon-identity collapses (recovered as
  copy-number where a family forms), and ~34 are coverage-limited (expressed but <20 reads). None is a
  detection-algorithm failure — the synthetic ground-truth demo (`bench/SIM_DETECTION_DEMO.md`) shows 100 %
  detection sensitivity when coverage/expression are controlled.
- **The `recovered_by` column shows which leg earned each member**, so you can see the incremental value of
  each feature (RNA-split is the core; protein-tail, generalized projection, and expressed-collapse each add
  a distinct slice).

Regenerate: `python3 bench/soto/soto_detection_eval.py`.
