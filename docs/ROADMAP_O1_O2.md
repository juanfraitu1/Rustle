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
| truth | junction-anchored reads, AUDITED (`bench/o2_truth_audit.py`, §6eq): 62 valid of 117 | `o2scale/truth_valid.tsv` | valid; ambiguous anchors removed |

## Steps, in importance × effort order
| # | step | why it matters | effort | acceptance |
|---|---|---|---|---|
| **S1 ✓ (§6eq)** | ~~Audit the junction truth~~ DONE: 55/117 anchors ambiguous; truth re-issued as 62 valid anchors (`o2scale/truth_valid.tsv`). ⚠ Exact rescoring (§6eq addendum, row 687): GFF-model spliced 5/5; faithful units 0/3 in every column mode — all on the nested pair 9/10 ⟹ units must be per LOCUS (S2 first) on near-identical paralogues: for every low-MAPQ anchored read whose O2 call disagrees, realign the read (splice-aware) to BOTH hulls; if it carries the same splice on the chosen paralogue at comparable edit distance, the anchor was an annotation gap, not a truth | every later acceptance is scored on this truth | ½ session | a table of the disagreements with a verdict each; the truth set re-issued with the invalid anchors removed |
| **S2 ✓ (§6er)** | ~~The unit is the catalog row~~ DONE: defaults flipped with escape hatches (byte-identical), `mcl_ann/rna_units_v1` = 751 units with sd_depth/core_bp/nearest_ident/rep_frac; NPIP 29-locus cluster, LCR16u separate; O2 on it: 0 wrong / 0 right on 62 valid anchors (all ambiguous), MAPQ-60 placement 75/76. Original:: flip the three defaults (locus rule with representative-only semantics; core refinement TRIM arm; unit emission) so `mcl_families` writes a unit table with cluster, core hull, chain source, SD depth, curated-repeat fraction, CHM13 landing, nearest-paralogue identity | O2's input file IS O1's output; the bridge scripts retire | 1 session | byte-identity controls for each flag; NPIP/tandem/MCL4/MCL6/ZNF acceptance sets unchanged; `copy_assign` runs on the emitted table with no script in between |
| **S3 ✓ (column, §6er)** | ~~Abstention forecast column~~ `nearest_ident` emitted (313/751 units ≥0.99); forecast-vs-observed check folded into S5. Original:: per unit, max identity to any paralogue in its family (from O1's edges); per family, the count of forecast-resolvable copies | O2's `tied` becomes a prediction the catalog made before reads (ρ −0.55 on NPIP; ≥0.99 ⟹ share 0) | ½ session | column emitted; on NPIP the forecast matches the observed tie structure |
| **S4 (re-scoped, §6eq addendum)** | ~~O2-8d~~ the spliced star projection is accepted on GFF-model cores (5/5); on read-chain units it must be re-scored AFTER S2 (locus units); genomic modes stay OFF and off the critical path. Remaining S4 = consistency check of the spliced columns on the sweep (S5), not a rewrite. Original: O1's alignments are O2's columns (O2-8d): replace the star projection with a true multiple alignment of the hulls (`poa_msa_with_costs` takes n sequences) or a transitive-consistency filter over the pairwise alignments | the last place O2 re-derives homology, and the source of the min_p-0 wrong calls | 1–2 sessions | 117 anchored reads (post-S1): 0 confident disagreements; pairwise column counts stay ≈ direct NM |
| **S5 ✓ (§6es)** | ~~Sweep~~ DONE: 87 families, 16,723 unit reads: 29% assigned / 7% tied / 64% ambiguous; MAPQ-60 placement 95%; forecast column does NOT predict (row 688). Original:: `copy_assign` over every SD-evidenced family of the unit catalog (artefacts excluded), reporting assigned / tied / forecast-vs-observed per family | the catalog-wide O2 result | 1 session (memory fixed, §6em) | table + the per-family forecast agreement |
| **S6 (draft, `docs/CHAPTER_O1_O2_COMPOSITION_DRAFT.md`)** | ~~Write the composition chapter~~ DRAFTED: sections 1–7 with the measured numbers; needs the advisor framing decision (block vs family stated up front) and figures. Original: with NPIP as the worked example (LCR16a core, EIF3C and chimeras excluded, read-chain units, abstention certificate) and the limitations (node boundaries, library-silent elements, slice conditioning) | | 1 session | text |

All steps done (S1–S5) or drafted (S6). O2 hygiene closed (§6et): `in_copy` + `catalog_copy_idx` columns, in-chain read counts, sweep v2 35/35. Canonical unit catalog: `mcl_ann/rna_units_v2`. Framing settled: the thesis defines the gene family (chapter §1).

## Standing decisions recorded
- Locus rule: representative-only semantics (attribution merges LCR16a with LCR16u, §6eg).
- Prune 1e-9 (§6ed). Inflation-insensitive 2.0–4.0 under it (§6ec).
- Cell-B artefacts (MCL0/NEW101, MCL15, MCL28) and MCL4 are never O2 inputs.
- O2's quotable NPIP statement until S4: abstention on the LCR16a cores; per-read calls not reliable.

## After S6 — the sweep's worst families (§6eu, 2026-09-05)
Three defects from adjudicating MCL32 (ZSCAN5) and MCL38 (ZNF875 / ZNF569-like): **D1** hull-clipped unit chains
(fixed, `mcl_families --units-follow-reads`, OFF; bounded by the annotated locus — the unbounded form engulfs, row 690);
**D2** the read-support PSV filter deletes the columns of unexpressed paralogues (row 689; `RUSTLE_PSV_READFILTER=0`,
default still ON); **D3** a molecule's secondary records scored as independent observations (row 691, OPEN — needs its
own PREREG; the sound form scores each molecule's sequence once). Paired sweep table in §6eu; default flips are
decisions, not patches.
