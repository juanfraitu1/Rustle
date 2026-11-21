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

## Result (2026-11-21): no pointer; IsoQuant reports fewer chains, it is not more precise at equal output

`r1853_isoquant_es/iq.txt` (`iq.py`).

**Mouse ES, the three sets:** both tools 12,648 / 23,240 = 0.544 (median 17 Rustle reads); Rustle only 4,425 / 18,401 = 0.240
(median 5 reads); IsoQuant only 795 / 5,268 = 0.151. The raw gap (0.410 against 0.472) is the size of the tool-only sets:
Rustle adds 18,401 chains to the shared ones, IsoQuant 5,268. Rustle-only chains by reads: 2 reads 0.087, 3-5 0.121, 6-20
0.180, over 20 0.573 (4,640 chains, 2,659 correct, which IsoQuant does not report at all).

**By intron count, mouse ES** (reported, precision; Rustle cut by reads to IsoQuant's correct count in the bin):
| introns | Rustle | IsoQuant | Rustle at IsoQuant's correct count | IsoQuant minus that |
|---|---|---|---|---|
| 1 | 1,044 / 2,902 = 0.360 | 996 / 2,286 = 0.436 | no cut possible (0.360) | +0.076 |
| 2-3 | 2,872 / 5,872 = 0.489 | 2,516 / 5,156 = 0.488 | reads >= 3: 0.564 | -0.076 |
| 4-9 | 6,604 / 14,746 = 0.448 | 5,194 / 10,707 = 0.485 | reads >= 10: 0.670 | -0.185 |
| >= 10 | 6,553 / 18,121 = 0.362 | 4,737 / 10,359 = 0.457 | reads >= 19: 0.641 | -0.184 |
The same column on the other datasets (IsoQuant minus Rustle's cut; 1 / 2-3 / 4-9 / >= 10 introns): chimpanzee -0.046 /
-0.021 / -0.076 / -0.026; bonobo -0.032 / -0.019 / -0.083 / -0.065; KB3781 -0.157 / -0.180 / -0.195 / -0.159; fly simulated
(no cut) / +0.004 / -0.009 / -0.039; mouse simulated (no cut) / (no cut) / -0.003 / -0.013.

- **No pointer** (IsoQuant ahead by 0.03 at matched count on ES and on two of the three other real libraries): in every bin
  of two or more introns Rustle at IsoQuant's correct count is MORE precise on all four real libraries, by 0.02-0.20.
- **The one place IsoQuant is ahead is single-intron chains, off the apes:** mouse ES 0.436 against 0.360 with nearly the
  same correct count (996 / 1,044); fly simulated 1,396 correct at 0.991 against 1,394 at 0.965; mouse simulated 2,074 at
  0.972 against 1,982 at 0.949 (more correct AND more precise). On chimpanzee, bonobo and KB3781 Rustle is ahead in that bin.
  It misses the registered bar (needs two other real libraries) and is 7% of the ES chains, so it is an observation: about
  570 fewer wrong single-intron chains on ES would move its precision by about 0.006.
- IsoQuant-only correct chains: 795 on ES, 727 of them absent from Rustle's `x.gtf` as well (not a filter loss); on mouse
  simulated 1,905, spread over all intron bins. These are chains Rustle does not build, 4-8% of IsoQuant's correct ones.
- Expectation: held (no pointer at two or more introns; the doubt about single-intron chains was the right place to look).

**Reading (mine).** IsoQuant looks better on mouse ES because it stops earlier: it reports 28,508 chains where Rustle
reports 41,641, and the 18,401 chains only Rustle has are mostly low-count (median 5 reads, 0.24 annotated). At the same
number of correct chains Rustle is 0.14 more precise overall and ahead in every multi-intron bin. Nothing in IsoQuant's
multi-intron handling is worth porting. The single-intron bin is the only structural difference in IsoQuant's favour, on
mouse and fly only; IsoQuant's rule there (novel 2-exon models need a strand from canonical sites and a mapping-quality
cutoff) is the one item of row 1838's list that has no Rustle counterpart and was never tested.
