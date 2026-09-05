# Roadmap: making O1 slide into O2 (2026-09-05)

**Principle (docs/O1_O2_COMPOSITION.md):** O1's vertex set IS O2's path set, and — extended today — O1's
pairwise alignments ARE O2's PSV columns. O2 must consume what O1 emits and re-derive nothing.
Canonical catalogs: `mcl_ann/rna_bp1_p9` (3 contigs) and `mcl_ann/gw_bp1_p9` (genome-wide), prune 1e-9.
Every prior plan file is archived under `docs/archive/` (superseded, not wrong); the ledger `o1_ledger.md`
and the register `NEGATIVE_RESULTS_REGISTER.md` remain the record of why.

## The chain, as measured (ledger §6ee–§6ep)
| stage | what exists | flag / file | status |
|---|---|---|---|
| node = locus | exon-union-overlap components, longest exon-union as representative | `mcl_families --merge-overlapping-loci` | OFF; attribution semantics reconstruct the duplication BLOCK (§6eg) — ship with representative-only semantics + the CHM13 landing as co-signatory |
| core = the family's shared duplicon | SEDEF depth ≥ half the family; trim chimeras to the hull | `--core-refine --sedef` → `.cores.tsv` | OFF; TRIM arm accepted, DROP arm is a flag (§6ei, row 677) |
| unit = core + read-supported exon chain | majority-splice rule, GFF fallback | `--emit-units --bam --fasta` → `.units.tsv/.fa/.regions` (+ `core_hull`) | OFF; faithful (§6en) |
| O2 input | `copy_assign --families --copies-fa` reads `.units.tsv` directly | contract in `catalog_input.rs` | works; `core_hull` optional column (§6ep) |
| O2 columns | spliced star projection (default) / hull-span chained minimap2 (`--psv-genomic`) | | pairwise correct; family-level star projection is the open defect (§6ep) |
| O2 decision | assign-or-abstain, `min_p < α/(n−1)` | | unchanged; NPIP = abstention certificate |
| truth | junction-anchored reads (`bench/o2_junction_truth.py`, `bench/o2_score_fixed_truth.py`) | | valid on genuinely different junctions; UNAUDITED on near-identical paralogues (row 686) |

## Steps, in importance × effort order
| # | step | why it matters | effort | acceptance |
|---|---|---|---|---|
| **S1** | **Audit the junction truth** on near-identical paralogues: for every low-MAPQ anchored read whose O2 call disagrees, realign the read (splice-aware) to BOTH hulls; if it carries the same splice on the chosen paralogue at comparable edit distance, the anchor was an annotation gap, not a truth | every later acceptance is scored on this truth | ½ session | a table of the disagreements with a verdict each; the truth set re-issued with the invalid anchors removed |
| **S2** | **The unit is the catalog row**: flip the three defaults (locus rule with representative-only semantics; core refinement TRIM arm; unit emission) so `mcl_families` writes a unit table with cluster, core hull, chain source, SD depth, curated-repeat fraction, CHM13 landing, nearest-paralogue identity | O2's input file IS O1's output; the bridge scripts retire | 1 session | byte-identity controls for each flag; NPIP/tandem/MCL4/MCL6/ZNF acceptance sets unchanged; `copy_assign` runs on the emitted table with no script in between |
| **S3** | **Abstention forecast column**: per unit, max identity to any paralogue in its family (from O1's edges); per family, the count of forecast-resolvable copies | O2's `tied` becomes a prediction the catalog made before reads (ρ −0.55 on NPIP; ≥0.99 ⟹ share 0) | ½ session | column emitted; on NPIP the forecast matches the observed tie structure |
| **S4** | **O1's alignments are O2's columns (O2-8d)**: replace the star projection with a true multiple alignment of the hulls (`poa_msa_with_costs` takes n sequences) or a transitive-consistency filter over the pairwise alignments | the last place O2 re-derives homology, and the source of the min_p-0 wrong calls | 1–2 sessions | 117 anchored reads (post-S1): 0 confident disagreements; pairwise column counts stay ≈ direct NM |
| **S5** | **Sweep**: `copy_assign` over every SD-evidenced family of the unit catalog (artefacts excluded), reporting assigned / tied / forecast-vs-observed per family | the catalog-wide O2 result | 1 session (memory fixed, §6em) | table + the per-family forecast agreement |
| **S6** | **Write the composition chapter** with NPIP as the worked example (LCR16a core, EIF3C and chimeras excluded, read-chain units, abstention certificate) and the limitations (node boundaries, library-silent elements, slice conditioning) | | 1 session | text |

Dependencies: S1 gates S4; S2 gates S5; S3 is independent. Order of execution: S1 → S2 → S3 → S4 → S5 → S6.

## Standing decisions recorded
- Locus rule: representative-only semantics (attribution merges LCR16a with LCR16u, §6eg).
- Prune 1e-9 (§6ed). Inflation-insensitive 2.0–4.0 under it (§6ec).
- Cell-B artefacts (MCL0/NEW101, MCL15, MCL28) and MCL4 are never O2 inputs.
- O2's quotable NPIP statement until S4: abstention on the LCR16a cores; per-read calls not reliable.
