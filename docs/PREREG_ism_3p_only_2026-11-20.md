# Row 1843: the sub-chain collapse spares chains cut at the 3' side only (pre-registered 2026-11-20)

User (2026-11-20): implement the exemption for chains cut only at the 3' side. Lead = row 1841's post hoc table: among the
chains the collapse removes, those no longer chain reaches upstream of (they only END early) are 0.22-0.34 correct on the
mouse libraries and 0.96-0.98 on fly and mouse simulated, against 0.04-0.07 and 0.13-0.34 for the 5'-cut ones; on the apes
0.010-0.023. That table was seen before this row; this row is the registered test of the rule built from it. Workspace
`r1843_ism_3p_only`. Committed before the code was written and before any run.

## The rule, fixed
Switch `copy_assign --polish-ism-ratio-3p R` (`RUSTLE_POLISH_ISM_RATIO_3P=R` in the driver; unset = off = today's output,
byte-identical). A multi-exon transcript Y of the polish input is *3'-only* when its chain is a contiguous part of at least
one longer chain of its contig and strand and NO such container has introns upstream of Y's first intron (5' side by the
strand). For a 3'-only Y the collapse uses the read ratio R in place of `--polish-ism-ratio` (0.3): Y is kept against a
container X when reads(Y) >= R x reads(X). R = 0 spares it outright. Every other polish step still judges Y; chains cut at
the 5' side or at both sides are judged as today.
Not the existing `--polish-ism-3p` (refuted, rows 861 / 1064): that one spares every sub-chain not flush with the
container's 3' end, pieces from the middle included (row 1841: 0.05-0.09 correct on mouse, 0.002 on apes), and does not ask
whether another container reaches upstream.

## Arms and data
`base` (today's defaults), `r3p0` (R = 0), `r3p01` (R = 0.1) on the twelve datasets of rows 1837-1842; `r0` of row 1839 (the
collapse off) as the ceiling, not rerun. Stage `assemble`, `RUSTLE_PTC=0`, `RUSTLE_ONEREAD_TIER=0`, scorer of row 1839.
Reused datasets, said. R = 0 and 0.1 chosen a priori (the exemption asked for, and row 1800's non-primate ratio).

## Reading, fixed now (`x.final.gtf`; `x.gtf` beside it)
- **Pays** on a dataset as in rows 1839 / 1840: the marginal precision of the chains the arm adds is at least the dataset's
  `base` precision, or precision falls by at most 0.005 while correct chains rise by at least 1%.
- **Harmless** on a dataset: precision falls by at most 0.005 and correct chains do not fall.
- **PASS** = an arm pays on at least 6 of 7 non-primate and 4 of 5 primate datasets.
- **NON-PRIMATE CANDIDATE** = an arm pays on at least 6 of 7 non-primate datasets and is harmless on all 5 primate ones.
- Anything else = FAIL, with the table. A PASS or a NON-PRIMATE CANDIDATE makes the switch a candidate; a default is the
  user's call and is then checked on one dataset not used here.
- Parity: without the switch the new binary's `x.gtf` and `x.final.gtf` are byte-identical to row 1840's `base` (md5).
- Reported: chains added (correct / wrong) per arm, their share of the ceiling's correct chains, and the precision of the
  added chains against row 1841's 3'-only figures (does the rule reach the class the table described).

## Expectation, written before the build (from row 1841's table, `x.gtf` arithmetic)
`r3p0` pays on fly and mouse simulated (added chains more precise than the assembly) and on mouse ES by the second clause
(+2.5% correct, precision about -0.003). I expect it NOT to pay on zebrafish (+0.6% correct), heart (+0.9%) and adrenal
(precision about -0.015), so 3 or 4 of 7: FAIL as registered, unless the final filter removes the wrong additions. `r3p01`
adds less and should lose less precision; I am unsure it reaches 1% anywhere it matters. On the apes both arms add wrong
chains at 0.01-0.02; harmless on precision only where the 3'-only class is small against the assembly.

## Limits
- "No container reaches upstream" is judged on the polish input, which holds more chains than row 1841's output files: the
  class is narrower than the one measured there.
- A 3'-only chain is still exposed to the isoform fraction, the floors and the final filter, which were set with the
  collapse at 0.3.
- Simulated reads have idealised 3' ends and no internal priming, which is the artifact this class would carry on real
  libraries; real-library truths are incomplete.

## Result (2026-11-20): FAIL; pays on 3 of 7 non-primate and 0 of 5 primate, as expected beforehand

Jobs 4490354-4490366, 36 arms exit 0; `r1843_ism_3p_only/sep.NAME.txt` (`sep.py`). Code `5e5b80c2`, binaries in
`~/ws/rustle_target_b`. Parity: `base` is byte-identical to row 1840's `base` (md5 of `x.gtf` and `x.final.gtf`; fly, mouse
ES, A119b). Full test suite 1,586 passed, 0 failed.

`x.final.gtf`, arm against `base`: correct added / wrong added, precision base -> arm.
| | `r3p0` (spared) | | `r3p01` (ratio 0.1) | |
|---|---|---|---|---|
| fly simulated | +14 / 0, 0.964 -> 0.964 | PAYS | +12 / 0 | PAYS |
| mouse simulated | +545 / +6 (3 correct lost), 0.953 -> 0.954 | PAYS | +320 / +5 | PAYS |
| zebrafish simulated | +27 / +27, 0.891 -> 0.890 | no (+0.2% correct) | +20 / +14 | no |
| mouse ES | +290 / +1,039, 0.347 -> 0.344 | PAYS (+1.6% correct) | +161 / +499, 0.346 | no (+0.9%) |
| mouse heart | +91 / +174, 0.577 -> 0.574 | no (+0.7%) | +63 / +101, 0.575 | no |
| mouse adrenal | +175 / +469, 0.578 -> 0.570 | no | +110 / +264, 0.574 | no |
| mouse cortex | +142 / +304, 0.487 -> 0.485 | no (+0.9%) | +103 / +186, 0.486 | no |
| OR6737 | +47 / +1,741, 0.377 -> 0.369 | no | +39 / +1,000, 0.372 | no |
| A119b | +188 / +3,946, 0.219 -> 0.215 | no | +148 / +2,878, 0.216 | no |
| chimpanzee | +30 / +1,479, 0.416 -> 0.404 | no | +20 / +973, 0.408 | no |
| bonobo | +40 / +1,487, 0.337 -> 0.328 | no | +29 / +894, 0.332 | no |
| KB3781 | +48 / +2,593, 0.373 -> 0.360 | no | +34 / +1,503, 0.366 | no |

- **`r3p0`: 3 of 7 non-primate, 0 of 5 primate; not harmless on four of the five apes** (precision -0.008 to -0.013; A119b
  -0.004). **`r3p01`: 2 of 7 and 0 of 5.** Neither PASS nor NON-PRIMATE CANDIDATE.
- The rule reaches the class row 1841 described: in `x.gtf` the chains it adds are 0.955 / 0.984 / 0.403 correct on the
  simulated sets, 0.24 / 0.34 / 0.22 / 0.34 on the mouse libraries, 0.010-0.023 on the apes, the figures of that table, and
  the rest of what ratio 0 adds is 0.04-0.08 on mouse and 0.004-0.006 on apes. The class separation is real; the class is
  just not precise enough on real libraries: 2-3 wrong chains per correct one on mouse (assembly: about 1 per correct on
  heart and adrenal), 20-55 on apes.
- **After the final filter the exemption equals the collapse switched off.** In `final.gtf` `r3p0` adds what ratio 0 adds
  (mouse ES +290 / +1,039 against +282 / +1,010; share of the ceiling's correct chains 0.94-1.00 on all twelve): the 5'
  truncation rule of the filter already removes the 5'-cut chains, so what a lower ratio lets through to `final.gtf` IS the
  3'-only class. Rows 1839 and 1843 measured the same additions.
- Expectation: held on every dataset named (pays fly, mouse simulated, ES; not zebrafish, heart, adrenal); the final filter
  did not remove the wrong additions.

**Decision:** not a candidate for a default. The switch stays in the code as an opt-in (off by default, parity checked); it
is the precise form of "more recall on mouse-like libraries" (mouse simulated +2.0% correct at 0.99, mouse ES +1.6% at 0.22)
and removing it is the user's call. The shared `~/rustle_target` was not rebuilt with it.
