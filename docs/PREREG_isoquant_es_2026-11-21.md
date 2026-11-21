# Row 1853: where is IsoQuant ahead on mouse ES, and is any of it a rule to port? (pre-registered 2026-11-21)

User (2026-11-21): understand why IsoQuant is much better on mouse ES cells and whether anything should be ported. Known:
IsoQuant reports 28,508 chains there at 0.472, Rustle 41,641 at 0.410; Rustle cut by reads to IsoQuant's correct count is at
0.611 (row 1850); row 1838 read IsoQuant's rules from its source and found "operating point" on six datasets, on the output
of 2026-11-20 morning. Tried since from that list: the end peak for contained chains (rows 1840, 1841, failed), the junction
merge (ratio 20, now default). Not looked at: the tools side by side by transcript structure, and the 2-exon rule.
Measurement only. Workspace `r1853_isoquant_es` (`iq.py`). Committed before the script was written or run.

## Data
Today's `final.gtf` (`r1849_dup_count_rerun/NAME/dc`) and IsoQuant de novo (`apetools/*/iq_dn`) on the six datasets that have
an IsoQuant run: mouse ES, chimpanzee, bonobo, KB3781 (real), fly and mouse simulated. Exact intron chains, truths of row 1849.

## Measured, per dataset
1. Sets: both tools, Rustle only, IsoQuant only: chains, correct, precision; median Rustle `reads` of the first two.
2. Per intron-count bin (1, 2-3, 4-9, >= 10): reported, correct, precision for each tool.
3. Rustle-only chains: precision by `reads` (2, 3-5, 6-20, > 20), and the share with a single intron.
4. IsoQuant-only chains that are correct: in Rustle's `x.gtf` (removed by the final filter) or absent from it; by intron bin.
5. Rustle cut by reads within each intron bin to IsoQuant's correct count in that bin: precision beside IsoQuant's.

## Reading, fixed now
- **A pointer for a port** = an intron bin in which, with Rustle cut (item 5) to IsoQuant's correct count in that bin,
  IsoQuant's precision is still higher by at least 0.03, on mouse ES and on at least two of the three other real libraries.
  A pointer names where IsoQuant does something a read cut does not; the rule behind it is then read from its source and
  becomes its own registered row.
- No pointer = IsoQuant is more precise on mouse ES because it reports fewer chains, not because of a rule Rustle lacks
  (row 1838's verdict, confirmed on today's output and by structure).
- Reported whatever the verdict: how much of the gap in raw precision is the Rustle-only set, and what it is made of.

## Expectation
No pointer in the bins of 2 or more introns. I am unsure about single-intron chains: IsoQuant has a rule for novel 2-exon
models (strand from canonical sites, a mapping-quality cutoff) and Rustle has none; row 1851 put Rustle's single-intron
chains at 0.360 on ES.

## Limits
- IsoQuant ran with its defaults in an earlier row; its output carries no read count comparable to Rustle's, so the matched
  cut is on Rustle's side only.
- Real-library truths are incomplete (row 1851); the simulated sets are shown as the complete-truth contrast.
