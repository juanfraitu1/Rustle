# Row 1854: would IsoQuant's rules for 2-exon models separate Rustle's wrong single-intron chains? (pre-registered 2026-11-21)

User (2026-11-21): try it, after row 1853 (single-intron chains are the one group where IsoQuant is ahead: mouse ES 0.436
against 0.360, fly simulated 0.991 / 0.965, mouse simulated 0.972 / 0.949; not on the three ape libraries). Measurement only,
no assembler code. Workspace `r1854_single_intron` (`one.py`). Committed before the script was written or run.

## IsoQuant's rules, read from the source before this row (4.1.0, `model_filter.py`, preset `default_pacbio`)
For a novel model of at most two exons: (a) the mean mapping quality of its reads must reach `simple_models_mapq_cutoff` = 30;
(b) a lone intron (a component with one intron) needs max(2, `min_mono_count_rel` 0.005 x the highest coverage of the
component it overlaps) reads. A polyA tail is NOT required for 2-exon models in this preset (seen in the ES run's log).

## Data
Single-intron chains of today's `final.gtf` (`r1849_dup_count_rerun/NAME/dc`) on the six datasets with an IsoQuant run: mouse
ES, fly and mouse simulated (IsoQuant ahead), chimpanzee, bonobo, KB3781 (Rustle ahead). CORRECT = exact chain in the truth.

## Measured, per chain
From the BAM, the primary records whose only intron is the chain's: their number and mean MAPQ. From `final.gtf`: the highest
`reads` of a multi-intron transcript whose span overlaps the chain's intron on the same strand, and on the other strand.
Signatures (a chain HOLDS one = the rule would drop it):
| | holds when |
|---|---|
| Q30 | mean MAPQ < 30 (IsoQuant's rule) |
| Q50 | mean MAPQ < 50 |
| REL | `reads` < 0.005 x the strongest overlapping multi-intron transcript on the same strand (IsoQuant's ratio) |
| REL02 | the same at 0.02 (Rustle's isoform fraction, across gene_ids) |
| ANTI | overlaps a multi-intron transcript on the other strand only, with fewer reads than it |
Beside them: precision by mean-MAPQ bin (< 20, 20-30, 30-50, 50-59, 60), and the share of wrong single-intron chains that
IsoQuant also reports.

## Reading, fixed now
- A signature **separates** on a dataset when the chains holding it are at most half as precise as all single-intron chains
  there, they are at least 20% of the wrong single-intron chains, and at most 2% of the correct single-intron chains hold it.
- **Candidate** = separates on mouse ES, fly simulated and mouse simulated, and holds at most 2% of the correct single-intron
  chains on each of the three ape libraries. A candidate becomes a build row (opt-in switch, judged on the twelve datasets).
- Otherwise: no port; the single-intron gap is reported with what the wrong chains are.

## Expectation
Q30 will hold few chains: Rustle's inputs are primary alignments and most have MAPQ 60, so I expect it to miss the 20% bar.
REL and ANTI are the likelier ones on mouse ES; on the simulated sets I do not know what the wrong single-intron chains are
(there is no noise in the simulation), which is what this row should show.

## Limits
- Rustle's `reads` and the BAM count can differ (junction snap, collapse); the BAM count is used for MAPQ only.
- IsoQuant's "component" is its intron graph's; the overlap test here is a stand-in.
- Single-intron chains are 7% of the ES output: the ceiling of any rule here is about 0.006 of ES precision.
