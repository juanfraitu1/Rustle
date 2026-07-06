# Exon-colinearity family merge and split stages (`family_merge_colinear.py`)

## What this is

A set of **post-demote** stages in `bench/family_rna_refine.py` that adjust the gamma-quasi-clique
refined families using exon colinearity:

1. **Colinear family merge** (default-ON, 5th stage) — merges catalog families split by gamma refine at
   weak edges but still sharing a short colinear exon block.
2. **Divergent same-chromosome merge** (default-ON, 6th stage) — a second pass with a lower exon-identity
   threshold but a longer colinear backbone, aimed at recent duplicons such as **HERC2**.
3. **Cross-chromosome domain-bridge split** (default-OFF) — an experimental opt-in that splits
   multi-chromosome families unless cross-chromosome components share strong homology or a curated
   gene-symbol link.  It is **off by default** because it over-splits real singleton cross-chromosome
   paralogs (R_oracle drops 51 → 46).

Motivating cases for the colinear merge:
- **MAGEA** — MAGEA9 sat in a sibling family because the Xq28 array edges are weak after gamma refine.
  The merge recovers it into the main MAGEA family.
- **GSTM** — the divergent GSTM1/2/4/5 sub-tandem was split across families; the merge rejoins them.

The risk is domain-share over-merges:
- **ANKRD18 + ANKRD36C** share a few exons in colinear order but are different ankyrin-repeat proteins.
  A fixed "≥2 colinear exons" rule would merge them.

## Colinear merge rule

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

## Divergent same-chromosome merge rule

A second pass for recent duplicons that have diverged below the 0.70 exon-identity floor while keeping
a longer colinear backbone (`ID_THRESH = 0.55`, `MIN_COLINEAR = 4`, `MIN_ADJACENT_JUNCTIONS = 3`,
`WINDOW_BP = 10 Mb`).

A pair of same-chromosome blocks merges iff:
- they share at least one annotated gene symbol,
- the smaller block's non-LOC gene symbols are a subset of the larger block's non-LOC symbols (so a
  GSTM2 block cannot pull in an unrelated named gene such as TRIP13),
- the colinear backbone passes the identity/exon/junction thresholds,
- the pair passes the antisense/reciprocal-overlap and repeat-hub gates.

This recovers the **HERC2** duplicon cluster (8 loci in one family) without reopening domain-share
over-merges.  It also merges additional GSTM2-related fragments on the same chromosome; the resulting
GSTM2 family is large but does not increase the diploid-oracle distinct-FP count.

## Cross-chromosome domain-bridge split rule

For each multi-chromosome family, the gate challenges every cross-chromosome link where the smaller
chromosome-side component has exactly one locus.  The link is kept only if there is strong colinear
support (`id >= 0.80`, `col >= 4`, `junc >= 3`), an identical non-LOC gene symbol, or a shared
 gene-symbol root.

Despite these keep-rules, the gate splits 79 families in the opt-in mode and drops diploid-oracle
recall from 51/57 to 46/57 by removing real singleton cross-chromosome paralogs.  It is therefore
**default-OFF** and available as `--crosschrom-split` / `RUSTLE_CROSSCHROM_SPLIT=1` for experiments
aimed at singleton domain-bridge cases such as LOC129527827→ZNF92.

## Ablations

- `--no-colinear-merge` / `RUSTLE_NO_COLINEAR_MERGE=1` disables the colinear merge **and** the
  downstream divergent/crosschrom stages (they depend on the merge pass) and recovers the pre-merge
  catalog (`548029ad`, 566 families) byte-identically.
- `--no-divergent-merge` / `RUSTLE_NO_DIVERGENT_MERGE=1` disables only the divergent merge and recovers
  the pre-divergent catalog (`991913da`, 553 families) byte-identically.
- `--crosschrom-split` / `RUSTLE_CROSSCHROM_SPLIT=1` **enables** the cross-chromosome split (default OFF).
- The stages compose with all other `--no-*` flags and with `--high-precision`.

## Files

- Implementation: `bench/family_merge_colinear.py`
- Integration: `bench/family_rna_refine.py` (post-demote call)
- Diagnostics: `bench/diag_family_merge_candidates.py`

## Metrics impact (current default)

| metric | pre-merge (`--no-colinear-merge`) | post-merge (default) |
|--------|-----------------------------------|----------------------|
| families | 566 | 551 |
| colinear merge edges | 0 | 19 |
| divergent merge edges | 0 | 2 |
| crosschrom families split | 0 | 0 (default OFF) |
| R_oracle | 51/57 = 0.8947 | 51/57 = 0.8947 |
| P_Ep | 0.8940 | 0.8893 |
| P_oracle(dedup) | 0.9167 | 0.9111 |
| distinct FP blocks | 4 | 4 |
| MAGEA recall | 0.85 (11 cp) | **0.87 (13 cp)** |
| GSTM recall | 0.27 (3 cp) | **0.50 (6 cp)** |
| ANKRD18 recall | 1.00 → 0.89 (ANKRD36C pulled in) | **1.00** (ANKRD36C blocked) |

The colinear + divergent stages are **recall-neutral at the oracle level** (R_oracle held), recover
MAGEA and GSTM, and keep ANKRD18 precision by blocking the ANKRD36C domain-share.  P_Ep dips slightly
because the stages add weak merges, but P_oracle(dedup) and distinct-FP blocks remain at the pre-merge
level.

## Determinism

`PYTHONHASHSEED=0`, sorted iteration, deterministic family-id assignment.  The default catalog md5 is
`de430908` (551 families); `--no-colinear-merge` recovers `548029ad` byte-identically and
`--no-divergent-merge` recovers `991913da` byte-identically.
