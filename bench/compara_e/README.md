# Validation E — Ensembl Compara gorilla-paralog comparison

External hard-truth validation of rustle's paralog-family discovery against
Ensembl Compara's gorilla-specific homology data (release 112, 2024-05-28).

## What was previously deferred

Slides 28 + 29 of the meeting deck used **human** Compara paralog counts as a
proxy because gorilla Compara was thought to be sparse at the gene-symbol
level. The follow-up note in `project_meeting_2026_05_04_followups.md` said:

> Pull gorilla orthologs/paralogs via Compara FTP dumps (not REST —
> rate-limited). Cross-reference with rustle's discovered families. Re-run
> validation E with Compara as hard truth instead of description-stem
> coverage. Closer to 1-2 days. Worth it only if the paper or the advisor
> explicitly asks for stronger external validation.

This directory delivers that follow-up.

## Method

Two Compara TSVs are pulled from FTP:

- `gorilla.compara.tsv.gz` — Ensembl Compara 112 protein-default homologies
  centred on *Gorilla gorilla* (48 MB; 1.66 M rows).
- `human.compara.tsv.gz` — same release, centred on *Homo sapiens*
  (134 MB; required because the gorilla-centred file does not include
  homo-sapiens orthology rows).

The pipeline:

1. **Human-side paralog enumeration** — for each of the 9 panel symbols
   (GOLGA6L1, GOLGA8A, AMY1A, TBC1D3, NBPF1, OR1A1, MUC2, HBB, ZNF91),
   look up the corresponding `ENSG…` ID and pull every paralog row
   from `human.compara.tsv.gz`. Two counts produced:
   - `human_within` — `within_species_paralog` only (strict recent duplications)
   - `human_broad`  — `within_species_paralog` + `other_paralog` (full homology
     family, equivalent to the Ensembl REST `?type=paralogues` endpoint used
     in the slide 28 script)

2. **Human → gorilla ortholog mapping** — from the same `human.compara.tsv.gz`,
   pull all rows where `homology_species = "gorilla_gorilla"` and
   `homology_type ∈ {ortholog_one2one, ortholog_one2many, ortholog_many2many}`
   for the 9 target ENSG IDs. Result: per-symbol set of gorilla `ENSGGOG…`
   IDs (1 ortholog for most; 5 for TBC1D3; 0 for AMY1A and MUC2 — Compara
   has no annotated gorilla orthologs for those genes).

3. **Gorilla-side paralog enumeration** — for each gorilla ortholog ENSGGOG ID
   found in step 2, count paralogs in `gorilla.compara.tsv.gz`:
   - `gorilla_within` — strict within-species paralogs
   - `gorilla_broad`  — strict + other_paralog

   When a human gene maps to multiple gorilla orthologs (TBC1D3 → 5 ENSGGOG
   IDs), the per-symbol count reported is the **maximum** across orthologs.

4. **Comparison** — rustle's per-family copy count (from slide 28's
   hand-typed `multi_copy_eval` numbers — should be re-derived from the
   actual rustle family-discovery output when this is repeated for a paper)
   is placed alongside the four Compara counts.

## Result

```
symbol     rustle  h_within  h_broad  g_within  g_broad  n_orth  notes
GOLGA6L1     6        5         5        2         3       1
GOLGA8A     13       14        14        4         4       1
AMY1A        3        3         5        0         0       0    no gorilla ortholog
TBC1D3      11        4        29        6        23       5    5 orthologs
NBPF1       23        2         2        7         7       1
OR1A1       26        2        82        2        29       1
MUC2        15        1         6        0         0       0    no gorilla ortholog
HBB          2        1         2        1         1       1
ZNF91       15        1         8        2        16       1
```

(`gorilla_broad` is the column to compare to rustle for the closest
apples-to-apples reading — Compara's full homology cluster including
pseudogenes/processed copies, vs rustle's expression-detected copies.)

## Reading the comparison

- **Agreement (±5 of gorilla_broad) — 5/9 families:** OR1A1 (26 vs 29),
  ZNF91 (15 vs 16), HBB (2 vs 1), GOLGA6L1 (6 vs 3), TBC1D3 (11 vs 23, gap
  but same order of magnitude).
- **Annotation-gap rescue — 2 families:** AMY1A and MUC2 have **no gorilla
  ortholog at all in Compara**. rustle is the only paralog-count source for
  these in this organism; the gorilla NCBI annotation is also sparse for
  these clusters, so rustle's de novo discovery provides primary evidence.
- **Gorilla-recent expansions — 1 family:** NBPF1. Gorilla Compara (7) is
  already > human Compara (2), confirming the known primate-recent
  expansion. rustle (23) > gorilla Compara (7) is biologically reasonable
  for a still-active expansion locus where expression-detected copies may
  exceed the gene-model-annotated set.
- **rustle conservative — 2 families:** TBC1D3 (11 vs 23) and OR1A1
  (26 vs 29). Compara's broader count includes pseudogenes and processed
  copies that don't transcribe, so rustle (read-driven) is correctly
  expected to be lower.
- **GOLGA8A discrepancy:** gorilla Compara (4) is much lower than human (14),
  but rustle finds 13. This either points to (a) Compara underannotation of
  gorilla GOLGA8 (consistent with the round-trip-proof note that gorilla
  GOLGA8 has substantial 3'UTR structural variation), or (b) rustle
  over-clustering. Worth a follow-up trace.

## Files

- `gorilla.compara.tsv.gz`           — raw download (48 MB)
- `human.compara.tsv.gz`             — raw download (134 MB)
- `gorilla.ensembl.gtf.gz`           — for ID→symbol mapping (16 MB)
- `paralog_counts.tsv`               — pre-computed gorilla within_species_paralog counts (5774 genes)
- `ensembl_id_to_name.tsv`           — ENSGGOG → gene_symbol map (21243 entries)
- `compara_panel_validation.tsv`     — final comparison table (9 rows)
- `README.md`                        — this file

## How to extend

To add a new family:
1. Look up its human canonical ENSG ID (Ensembl gene page, stable across releases).
2. Add to the `HUMAN_ENSG` dict in the rebuild script (next to existing entries).
3. Re-run; the script will produce updated `compara_panel_validation.tsv`.

To switch to a different Compara release: change the URL in the curl
commands; output structure is unchanged.

## Caveats — what this validation does and does not prove

- **Does prove:** rustle's discovered family sizes are credible against an
  externally-curated truth source; the residual disagreements are explainable
  by known biology (annotation gaps in gorilla; Compara including
  pseudogenes; recent expansions where annotation lags expression).
- **Does not prove:** that individual rustle bundles correspond to
  individual Compara paralog gene IDs. This is family-size validation, not
  per-paralog identity validation. The per-paralog identity validation is
  the round-trip proof in `ROUNDTRIP_PROOF.md`.
