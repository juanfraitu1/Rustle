# Exon-colinearity family merge (`family_merge_colinear.py`)

## What this is

A **post-demote, default-ON** stage in `bench/family_rna_refine.py` that merges catalog families which
are split by the gamma-quasi-clique refine at weak edges, but which still share a short colinear exon
block between any pair of loci.  It is the 5th default-ON stage after repeat-hub gate, recombinant-split
gate, multi-repeat-bridge gate, and antisense/reciprocal-overlap gate.

Motivating cases:
- **MAGEA** — MAGEA9 sat in a sibling family because the Xq28 array edges are weak after gamma refine.
  The merge recovers it into the main MAGEA family (fam496, 13 copies).
- **GSTM** — the divergent GSTM1/2/4/5 sub-tandem was split across families; the merge rejoins them
  (fam17, 6 copies).

The risk is domain-share over-merges:
- **ANKRD18 + ANKRD36C** share a few exons in colinear order but are different ankyrin-repeat proteins.
  A fixed "≥2 colinear exons" rule would merge them.

## Rule

A pair of catalog blocks is considered if:
1. they share a dominant chromosome, AND
2. either they share an annotated gene symbol, OR their spans are within `WINDOW_BP = 5 Mb`, OR there
   is at least one raw homology edge (from the pre-gate `denovo_family_edges.tsv` graph) connecting them.

They merge iff:
- the best strict-LIS colinear shared-exon count between any two loci is `>= MIN_COLINEAR = 2` at
  `ID_THRESH = 0.70`, AND
- the bridge is not blocked by the antisense/reciprocal-overlap or repeat-hub gates, AND
- **window-only pairs** (no shared gene symbol and no raw edge) additionally pass an **adaptive
  adjacent-junction floor**.

### Adaptive adjacent-junction floor

For window-only pairs we count **adjacent preserved splice junctions**: consecutive matched exon pairs
`(i, i+1)` in family A that map to consecutive exons `(j, j+1)` in family B.  The required floor is
`min(MIN_ADJACENT_JUNCTIONS, colinear_exons - 1)`:

- a 2-exon colinear hit needs **1** adjacent junction,
- a ≥3-exon colinear hit needs **2** adjacent junctions.

This blocks ANKRD18 + ANKRD36C (`col=3, junc=1`) while keeping GSTM1/2/4/5 (`col=2, junc=1`) and
MAGEA9 (`col=3, junc=2`, supported by a raw edge anyway).

## Ablations

- `--no-colinear-merge` / `RUSTLE_NO_COLINEAR_MERGE=1` disables ONLY this stage and recovers the
  pre-merge catalog (`548029ad`, 566 families) byte-identically.
- The stage composes with all other `--no-*` flags and with `--high-precision`.

## Files

- Implementation: `bench/family_merge_colinear.py`
- Integration: `bench/family_rna_refine.py` (post-demote call)
- Diagnostics: `bench/diag_family_merge_candidates.py`

## Metrics impact (current default)

| metric | pre-merge (`--no-colinear-merge`) | post-merge (default) |
|--------|-----------------------------------|----------------------|
| families | 566 | 553 |
| merge edges | 0 | 19 |
| R_oracle | 51/57 = 0.8947 | 51/57 = 0.8947 |
| P_Ep | 0.8940 | 0.8879 |
| P_oracle(dedup) | 0.9170 | 0.9149 |
| distinct FP blocks | 4 | 4 |
| MAGEA recall | 0.85 (11 cp) | **0.87 (13 cp)** |
| GSTM recall | 0.27 (3 cp) | **0.50 (6 cp)** |
| ANKRD18 recall | 1.00 → 0.89 (ANKRD36C pulled in) | **1.00** (ANKRD36C blocked) |

The stage is **recall-neutral at the oracle level** (R_oracle held), recovers two split known families,
and restores ANKRD18 precision by blocking the ANKRD36C domain-share.  P_Ep dips slightly because the
stage adds a few weak merges, but P_oracle(dedup) and distinct-FP blocks remain at the pre-merge level.

## Determinism

`PYTHONHASHSEED=0`, sorted iteration, deterministic family-id assignment.  The default catalog md5 is
`991913da` (553 families); `--no-colinear-merge` recovers `548029ad` byte-identically.
