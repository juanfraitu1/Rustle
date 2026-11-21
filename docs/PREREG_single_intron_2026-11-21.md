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

## Result (2026-11-21): no candidate; the mapping-quality rule works on mouse ES and fly, and would cost true chains on the apes

Jobs of `r1854_single_intron/jobids.txt`, six runs exit 0; `one.NAME.txt`, per chain `cand.NAME.tsv`.

Single-intron chains (correct / all) and the two IsoQuant rules: chains holding the signature (their precision), share of the
wrong and of the correct single-intron chains that hold it.
| dataset | single-intron chains | Q30: mean MAPQ < 30 | REL: under 0.005 of the transcript over it |
|---|---|---|---|
| mouse ES | 1,044 / 2,902 = 0.360 | 399 (0.045); 20.5% of wrong, 1.7% of correct: SEPARATES | 80 (0.100); 3.9%, 0.8% |
| fly simulated | 1,394 / 1,445 = 0.965 | 40 (0.150); 66.7% of wrong, 0.4% of correct: SEPARATES | 2 (1.000) |
| mouse simulated | 1,982 / 2,089 = 0.949 | 40 (0.575); 15.9% of wrong, 1.2% of correct: no | 47 (0.872); 5.6%, 2.1% |
| chimpanzee | 1,058 / 3,870 = 0.273 | 76 (0.553); 1.2% of wrong, 4.0% of correct | 3 |
| bonobo | 729 / 3,455 = 0.211 | 92 (0.457); 1.8% of wrong, 5.8% of correct | 5 |
| KB3781 | 959 / 3,747 = 0.256 | 33 (0.394); 0.7% of wrong, 1.4% of correct | 240 (0.021); 8.4%, 0.5% |
Q50, REL02 and ANTI separate nowhere except Q50 on fly (in the files).

- **Not a candidate.** Q30 separates on mouse ES and fly simulated but not on mouse simulated (the chains it holds are 0.575
  correct, more than half the group's 0.949, and 16% of the wrong), and on chimpanzee and bonobo it holds 4.0% and 5.8% of
  the CORRECT single-intron chains (bar 2%).
- **On the apes low mapping quality marks true chains, not wrong ones:** single-intron chains with mean MAPQ under 30 are
  0.55 / 0.46 / 0.39 correct on chimpanzee / bonobo / KB3781, against 0.27 / 0.21 / 0.26 for the group. A read that maps
  ambiguously there is typically in a duplicated gene, which is the kind of transcript this project is about; IsoQuant's
  cutoff removes them.
- On mouse ES the rule would remove 381 wrong and 18 correct single-intron chains: precision of the whole output 0.4100 ->
  0.4135. On fly simulated 34 wrong and 6 correct. That is the whole of IsoQuant's lead in this group on those two
  datasets (ES: 399 of the 616 chains IsoQuant does not report).
- What the wrong single-intron chains are, mouse ES: 1,858, of which 967 are also reported by IsoQuant; 0.076 correct when
  antisense to a multi-intron transcript only (303 chains), 0.311 inside one on the same strand, 0.445 with none over it;
  2-read chains 0.174. On the apes 67-87% of the wrong single-intron chains are also in IsoQuant.
- Expectation: wrong about Q30 holding too few chains on mouse ES (it holds 399, 20.5% of the wrong, just over the bar);
  REL and ANTI, which I thought likelier, do nothing.

**Decision:** no port. A mapping-quality cutoff for 2-exon models is library-dependent: useful on mouse ES and fly, harmful
where low mapping quality means a duplicated gene. With rows 1765, 1838 and 1853 this closes IsoQuant as a source of rules.
