# Row 1819: replace a weak junction by its strong neighbour instead of dropping the chain (pre-registered 2026-11-19)

User (2026-11-19), after row 1818 (an intron within 10 bp of a true one is 18-27% of the wrong chains left on complete truth, 22 of
28 on SIRV; dropping such chains costs 0.6-2.9% of the correct ones): how about replacing junctions?
**What exists.** The default polish junction snap (row 1800) already replaces junction K by neighbour J, but only when fewer than
two of K's reads align better at K than at J; it uses no support ratio (its design note: ratio merges cost chains on three
development contigs). The chains left are those the read test does not move. Registered before any replacement was computed.

**Rule tested, as a post-step on `PREFIX.gtf` at today's defaults.** Support of an intron = summed `reads` of the transcripts
of the assembly that carry it (contig, strand). Intron K of a chain is replaced by the best-supported intron J of the same contig
and strand with both ends within 10 bp of K's, when support(J) >= R x support(K). Chains that become identical collapse.
R in 2, 5, 10, 20, 50. No annotation, no new BAM pass.
**Development:** zebrafish simulated and OR6737. R* = the R with the fewest wrong chains (summed over the two) among those that
keep at least 0.995 of the correct chains on both; none = NO CANDIDATE.
**Check, R* frozen,** on fly simulated, mouse simulated, mouse ES, mouse heart, mouse adrenal, A119b, chimpanzee testis; SIRV shown.
ACROSS THE BOARD = on all seven the correct chains are at least 0.995 of before and precision is higher, and on at least five
the wrong chains fall by at least 3%. Reported per dataset: chains changed, wrong -> correct (of them not in the output before =
gained), wrong -> wrong, correct -> wrong, correct -> correct.

**Limits.** Support from the assembly's own read counts, not from the BAM; inside the assembler the read floor and later rules
would see the merged counts. Real alternative sites a few bases apart (NAGNAG) are what the rule destroys; the 0.995 bound is the
guard. A replaced chain is not re-checked for exon validity.

## Result (2026-11-19): R* = 20, ACROSS THE BOARD (a small effect)

Tables: workspace `r1819_junction_replace` (`dev.*.txt`, `c.*.txt`; script `replace.py`).
**Development** (correct chains kept / wrong chains fewer): R 2: zebrafish 0.991 / 2.2%, OR6737 0.948 / 5.5%; R 5: 0.995 / 2.5%,
0.974 / 4.2%; R 10: 0.997 / 2.6%, 0.990 / 3.3%; **R 20: 0.998 / 2.2%, 0.997 / 2.2%**; R 50: 1.000 / 0.8%, 1.000 / 0.6%.
Only R 20 and R 50 keep 0.995 on both; R* = 20.

**Check at R = 20:**

| dataset | correct before -> after (ratio) | precision before -> after | wrong chains fewer | wrong -> correct | wrong -> wrong | correct -> wrong | correct -> another correct |
|---|---|---|---|---|---|---|---|
| fly simulated | 7,082 -> 7,057 (0.996) | 0.926 -> 0.928 | 3.3% | 33 | 0 | 14 | 11 |
| mouse simulated | 27,731 -> 27,614 (0.996) | 0.903 -> 0.907 | 4.9% | 176 | 31 | 52 | 67 |
| mouse ES | 18,058 -> 18,007 (0.997) | 0.331 -> 0.349 | 7.8% | 2,425 | 469 | 2 | 50 |
| mouse heart | 13,287 -> 13,264 (0.998) | 0.546 -> 0.558 | 4.8% | 447 | 89 | 0 | 25 |
| mouse adrenal | 14,650 -> 14,617 (0.998) | 0.540 -> 0.552 | 4.7% | 550 | 48 | 0 | 34 |
| A119b | 39,049 -> 38,889 (0.996) | 0.158 -> 0.159 | 1.2% | 1,383 | 1,428 | 6 | 160 |
| chimpanzee testis | 20,940 -> 20,901 (0.998) | 0.351 -> 0.354 | 1.4% | 381 | 172 | 1 | 40 |
| SIRV (shown) | 53 -> 54 | 0.589 -> 0.675 | 29.7% | 11 | 0 | 0 | 0 |

All seven keep at least 0.995 of the correct chains with higher precision; wrong chains fall by at least 3% on five (not on the
two apes): the bar passes. Almost no chain is gained (0 to 4 per dataset): a wrong chain that becomes correct lands on a chain the
output already had, so the effect is the removal of a duplicate with a shifted site. The cost is real neighbouring sites:
a correct chain moved onto another correct chain (11 to 160 per dataset), and on simulated reads onto a wrong one (14, 52).
**Against dropping** (`wobble`, row 1818): dropping costs 0.6-2.9% of the correct chains for 3-14% of the wrong; replacing at
ratio 20 costs 0.2-0.4% for 1-8%.
**Not done:** the rule inside the assembler (the existing snap is in the binary; a ratio condition there needs a rebuild), and
the combination with the filter on `final.gtf`.

## Port into the assembler (row 1835, criteria written 2026-11-20 before the parity run)

User (2026-11-20): implement it. Opt-in `copy_assign --polish-junction-ratio R` (`RUSTLE_POLISH_JUNCTION_RATIO=R` in the driver;
unset or 0 = off): the rule above applied per contig to the final polished transcripts, through the junction snap's own
re-place and collapse step (`snap_apply`). Differences from `replace.py`, known before the run: a move that would empty an exon
is skipped per junction (the script refused the whole chain); it runs in the `final.gtf` pass too, where the chain score filter
then sees the merged read counts; `reads` of a collapsed transcript is the sum.
Parity run, workspace `r1835_junction_ratio_port`, binaries built into `rustle_target_b` (the shared ones untouched), arms `base`
(new binary, no switch) and `ratio20`, on the nine datasets of the development and check above:
- **P1.** `base` `x.gtf` and `x.final.gtf` are byte-identical to the installed binary's (rows 1831 / 1833 `base`).
- **P2.** `ratio20` on `x.gtf` meets row 1819's across-the-board criteria on its seven check datasets: at least 0.995 of the
  correct chains, higher precision, and wrong chains at least 3% fewer on at least five.
- Shown: `x.final.gtf`, zebrafish and OR6737 (development), gained and lost chains.
P1 and P2 = the switch is committed as opt-in. Default: the user's call.
