# Clair3-RNA editing filter for copy-assignment PSVs — design

**Date:** 2026-06-26 · **Status:** approved, implementing · **Scope:** `copy_assign.rs` + `copy_assign_pipeline.rs` + binary.

## Goal

Stop A-to-I RNA-editing sites from being trusted as copy-distinguishing PSVs in the significance gate. At an
A/G PSV that is actually an editing site, `εⱼ = e/3` is anti-conservative (a true-A read shows G at the
editing rate, ≫ e/3), so a c(A)-read edited to G can fake support for a G-copy. Approach **B** (detected):
flag *true* editing columns by their signature and downweight only those — preserving genuine A/G paralog SNVs.

PSV alleles are stored in **transcription (mRNA) orientation**, so editing is uniformly **A→G** in this frame
(no genomic-strand T↔C handling needed).

## Detector (`copy_assign_pipeline::detect_editing_columns(reads_obs, copies) -> Vec<bool>`)

Per PSV column, returns a flag. A column is **editing-flagged** iff:
1. it is an **A↔G** column — the copy-consensus alleles at it are exactly `{A, G}`; and
2. some copy shows **within-copy A/G heterogeneity** — provisionally assign each read to its argmax copy
   (match count over ALL columns; real PSVs dominate the minority of editing columns), then for each copy
   count `(n_A, n_G)` at the column among its reads; flag if any copy has `minor = min(n_A,n_G)` with
   `minor ≥ EDIT_MIN_READS (=2)` and `minor/(n_A+n_G) ≥ EDIT_MIN_FRAC (=0.05)`.

A real A/G PSV → each copy monomorphic (minor ≈ 0, ≪ EDIT_MIN_FRAC) → **not** flagged; the ~3e-4 sequencing
error is far below EDIT_MIN_FRAC so it never trips. Editing → the edited copy carries both A and G → flagged.

## Handling (`copy_assign::assign_read_editing`)

`assign_read(read, copies, p)` becomes a thin wrapper for `assign_read_editing(read, copies, p, &[])` (so the
~15 existing call sites/tests are untouched). The core gains a `editing_cols: &[bool]` parameter; in the
significance loop, for a distinguishing column `j` with `editing_cols[j] == true` and `p.rna_editing_filter`,
`εⱼ = max(e/3, p.edit_rate)` (default `edit_rate = 0.2`) — downweighting the flagged column to near-
uninformative for the certificate, in both directions (a flagged column is not a reliable genomic variant).
The likelihood *ranking* is unchanged (the column is still mild evidence); only the certificate stops
trusting it. Empty `editing_cols` ⇒ byte-identical to the current gate.

## Wiring

`assign_family_detailed` runs a pre-pass building the per-read `psv_obs` matrix (`fill_psv_obs`), computes
`editing_cols = detect_editing_columns(...)` once per family (when `rna_editing_filter`), then the main loop
calls `assign_read_editing(..., &editing_cols)` for both the combined and psv-only assignments.

## Params / CLI

- `AssignParams`: `edit_rate: f64 = 0.2`, `rna_editing_filter: bool = true` (default-on — it only acts on
  detected-het A↔G columns). Binary: `--edit-rate <f>`, `--no-editing-filter`.
- Detector thresholds `EDIT_MIN_FRAC = 0.05`, `EDIT_MIN_READS = 2` as documented constants.

## Tests

- `detect_editing_columns`: 2-copy family, column 0 = real A/G PSV (each copy monomorphic) → not flagged;
  column 1 = editing (a copy carries ~30% G) → flagged.
- `assign_read_editing`: a read resolvable only via a flagged column drops from Assigned to abstain; an
  unflagged read is unchanged; empty `editing_cols` == current behavior.
- Regression: full suite (currently 655) stays green.

## Caveat

The detector needs enough reads per copy to see the heterogeneity; very low-coverage copies may miss editing
(falls back to the uniform-error model). Documented, acceptable.
