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

## Definition switch (§6ev, 2026-09-05)
Both O1 definitions scored on Soto (one truth, one scorer): `docs/O1_DEFINITION_SWITCH.md`. Gaps G1–G5 tracked
there with owners; the O1 restatement (§5 of that file) awaits the user's decision; the old definition stays
opt-in. Next in order: G5 (RNA admission of unannotated loci — the unbuilt pivot stage), G3 (threshold
sensitivity on the anchors), Q9 (NPIPA/B), D3 pre-registration, the two §6eu default flips.

## G5 / Q9 closed, min_size flipped (§6ew–§6ex, 2026-09-05)
Q9 answered on the thesis definition (A/B is not a cut in the gorilla catalog; landing per locus is the certificate).
G5 decomposed: the node's shortfall was the size floor (19/54) and MCL fragmentation (9/54), annotation gaps 11/54;
RNA admission prototyped (`bench` only, recovers 2/11). **min_size default 2** (user); canonical 3-contig catalog =
`mcl_ann/rna_units_v4` (268 clusters, 432 units; every ≥3 cluster identical to v3). Open next: MCL fragmentation
(row 694), G3 sensitivity table, D3 pre-registration, the two §6eu default flips, a duplicon-name column per family.

## Row 694 root-caused (§6ey, 2026-09-05)
Not MCL: the fold-first locus rule loses every record that overlaps another family's record on exon bases, and
exon-less pseudogene records can carry no edge. Fix = three OFF flags (`--fold-within-clusters`, `--exonless-span`,
`--exonic-both-sides`); Soto 50 % floor: detection 0.751 → 0.798, band precision 0.949 → 0.943 (CIs overlap),
recall|both 0.942 → 0.949; gorilla `rna_units_v6`: NPIP 44/44 records one cluster (32 loci), LCR16u separate,
anchors intact, the L1 artefact blob dissolves (64/104 unclustered). **Default flips = user decision.** Then:
G3 sensitivity table, D3 pre-registration, the two §6eu flips, duplicon-name column, CHM13 landing for the 3 new NPIP loci.

## Defaults flipped, queue executed (§6ez, 2026-09-05)
Five defaults flipped (user); canonical catalog `mcl_ann/rna_units_v8` (+ `blocks.tsv`); O2 sweep on it (paired
35: 28.2 % → 33.7 % assigned @ 99.5 %); the 3 recovered NPIP loci landed (2 NPIP, 1 NPIPB6/EIF3CL block record);
G3 closed (thresholds = point in a plateau); D3 pre-registered (implementation next); duplication blocks emitted
(direct SD partners per family). Open: D3 implementation (O2-9), the MCL15/MCL28 residual artefacts, genome-wide
rebuild of the catalog under the new defaults, the human-side residual (PMS2P14-type exon-less block spans).

## D3 implemented (§6fa, 2026-09-05)
`copy_assign --molecule-observations` (OFF): one observation per molecule ("the read is the star": minimap2 of the
read to every copy's unit, columns = the read's own differing positions) + an ORIGIN CERTIFICATE (edits to the best
unit ~ Binomial(aligned, error_rate), else no candidate explains the read). Pooling the BAM records was refuted
(row 700). NPIP: 91 → 175 assigned, 0 wrong on the 62 anchors, and **285/433 low-MAPQ reads explained by no unit**
(row 702) — the abstention is two thirds node incompleteness. Paired 35: assigned share unchanged (33 %),
ambiguous 59 % → 17 %, tied 7 % → 49 %, agreement 99.6 %. **Default flip = user decision.** Next: the unexplained
reads are an O1/O3 item (units that lack their content; copies the catalog lacks) — count them per family and
feed the RNA-admission prototype; then the genome-wide rebuild.

## Read-star default, unit merge, unexplained reads located, genome-wide rebuild (§6fb, 2026-09-05)
`copy_assign --molecule-observations` ON (user); `origin_rejected` column. O1: `mcl_families --merge-overlapping-units`
ON (units sharing exon bases are one locus; row 703 — the record-level 69 % was nested-unit inflation); canonical
3-contig catalog **`rna_units_v9`** (349 units, 40 merged). The certificate's rejected reads locate 62 candidate
loci with ≥ 3 reads (`docs/o2_unexplained_v8_loci_2026-09-05.tsv`): genes that are the reads' origin and not
units of the family (dropped members, neighbouring clusters) — O1's roster gaps named by O2. RNA admission finds
an unannotated expressed LCR16u copy (NC_073242.2:35.99 Mb, 123 reads). Genome-wide catalog under the new
defaults: `mcl_ann/gw_units_v1` (`bench/gw_rebuild.sh`). Next: fold the 62 loci back into O1 (roster admission
by read origin), score gw_units_v1's anchors (tandem array, NPIP), the sub-30 %-score competitor question
(`-p 0.3`: a copy below it is not a candidate — state or test).

## §6fc — pairwise read-star is the O2 default; admission, floor and junction forms measured (2026-09-05)
Final default: read-star with pairwise certificates, every edit counted in the origin certificate, no junction
term, `-p 0.3`. Paired 35: 32.1 % assigned @ 99.90 % (record-level 33.7 % @ 99.51 %), MAPQ<60 assigned 132 → 232;
NPIP 91 → 316 assigned @ 201/201, 0 wrong anchors. Opt-ins measured and shipped OFF: `--origin-substitutions-only`
(row 707), `--read-star-junctions` (row 708), `RUSTLE_STAR_P` (row 704). Roster admission: safe, modest (row 706).
Open: the certificate's structure-vs-origin dichotomy (a per-read isoform model would lift the LCR16u rejections
without the substitution-only precision cost); the 60-member L1 artefact class genome-wide (classified, row-free);
score `gw_units_v1` with O2 genome-wide (GGO_ds.bam).

## §6fd — genomic read-star is the O2 default (2026-09-05)
Read vs each candidate's LOCUS (unit extent padded by the longest read), splice-aware; certificate = X + I + D +
unaligned bases over the read length; pairwise columns from genomic bases. PREREG held (0 wrong anchors, 100 %
placement agreement on the paired 35, 55.7 % assigned, MAPQ<60 386/1,158). Rows 709–710 correct the isoform
reading: the rejections were unit EXTENT. Open: the 10× time cost (minimap2 -x splice per family; batch across
families or index once per contig); O1 unit extents that exclude expressed exons (the anchor read's locus 4–16 kb
past its unit); genome-wide O2 on `gw_units_v1`.
