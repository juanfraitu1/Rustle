# PREREG — the SD-core definition on HUMAN NPIP (2026-09-07)

**Written before the substrate was built.** md5 in `soto_mcl/npip_hsa/PREREG.md5`.
⚠ Human and gorilla numbers are never pooled. This is a separate substrate with its own truth.

## Why a new substrate
The human catalog on disk (`soto_mcl/mcl_v12`) was built **without the core rule** (`core_refine false`,
no SEDEF) and over the **Soto slice**, whose PAF holds 696 gene spans — 13 of the 22 NPIP-named records,
the other 9 never proposed. Neither the definition nor the substrate is the one being asked about.

## Substrate
CHM13 v2.0, **chr16 and chr18** (every NPIP-named record lies there). Gene and pseudogene spans from
`hsa.gff` (CHM13 RefSeq) → FASTA → all-vs-all `minimap2 -x asm20 -c -X -N 50 -p 0.1 -t 4` →
`mcl_families --min-exonic-bp 1 --merge-overlapping-loci --core-refine --core-from-paf --emit-units
--bam soto_adj/soto.bam --fasta chm13v2.0.fa`. Every threshold at the gorilla defaults.
⚠ **No human SEDEF exists on this machine**, so the cores come from the run's own alignments
(`--core-from-paf`, §6fo — on gorilla NPIP the two routes agree on 30 of 32 members).
Reads are the human A119b IsoSeq already aligned to CHM13 (5,193 primaries in chr16:14–31 Mb).

## Truth
The **22 records whose gene symbol begins with NPIP** in CHM13 RefSeq: 21 on chr16, `NPIPB1P` on chr18.
This is a curated, annotation-derived truth, independent of anything the method computes — the opposite
failure from `lcr16a.bed`, whose circularity was caught today (register 727). ⚠ Its own bias is the
gene-symbol trap: an unnamed NPIP copy is invisible to it, so **sensitivity here is against named copies
only** and a unit outside the 22 is not necessarily false.

## Metrics (`bench/o1_eval.py`, unchanged)
Sensitivity (truth records rediscovered), specificity over MEMBERS with candidates counted separately,
and Hungarian 1:1 bipartite matching reported on the three interval forms (core hull, unit, locus extent):
median truth coverage and the in-band 0.5–2× fraction. Plus the selection-free count of clusters holding
≥1 NPIP record.

## Predictions
| # | prediction |
|---|---|
| **P1** | the 22 records fall into **≤ 2** clusters |
| **P2** | sensitivity **≥ 18/22** (members plus dropped candidates) |
| **P3** | precision of the dominant cluster **≥ 0.85** |
| **P4** | core-hull bipartite in-band **≥ 0.80** of matched pairs |
| **P5** | NPIPA and NPIPB records are **not** separated into different clusters by the core rule (Q9: they are one family with a shared core) |

## Interpretation fixed in advance
- P1+P5 holding ⟹ the definition transports to human NPIP and answers Q9 the same way it does in gorilla.
- P2 failing ⟹ report which records were never proposed (outside the annotation's reach) separately from
  those proposed and rejected; those are different failures.
- Any unit outside the 22 is reported as an unnamed candidate, never silently as a false positive.
