# PREREG r1618: the existing options that use incomplete reads, on the missed-isoform measure (2026-11-10)

User (2026-11-10): "can we try the option that tries to model incomplete reads?" Three opt-in modes already exist (October
ledgers; tested then on other measures). Written before any run here.

## Arms (gorilla OR6737, assemble only, default seed pool; default run = `r1614/default`)
- P. `RUSTLE_COMPAT_SUPPORT=1`: a chain below the read floor is admitted when every junction is carried by reads whose chain is
  a contiguous, end-compatible sub-chain of it.
- U. `RUSTLE_DIRECT_TRANSFRAGS=union`: default chains plus maximal transfrags credited with the reads they contain.
- X. `RUSTLE_EXTEND_TRANSFRAGS=1`: admitted chains also emitted extended along the best-linked junction pairs (StringTie-like
  stitching; the only arm that can reach the 398 chains no single read carries).
- S. `RUSTLE_READ_JUNCTION_SNAP=10`: read junctions snap to better-supported ones within 10 bp (the 323 junction-wobble
  misses; not an incomplete-read mode, added because it is the other untested class).

## Measure and bars (as r1616 / r1617)
Exact RefSeq intron chains. PASS = correct gained >= 100 and share correct among added chains >= the share among gorilla
default 2-read chains (0.170, seen). A passing arm is then run held out on A119b against 0.066 before any claim.
Descriptive: misses left per read class (0 / 1 / 2+ full-chain records); overall precision.

## Outcome on gorilla (2026-11-10, `r1618/`; default 73,348 chains, 25,794 correct, 0.352)
| arm | chains | correct | precision | gained (net) | added (net) | share correct of net added | bar |
|---|---|---|---|---|---|---|---|
| P compat | 74,907 | 26,111 | 0.349 | +317 | 1,559 | 0.203 |FAIL (new-chain reading 0.150, see correction below) |
| U union | 79,296 | 26,726 | 0.337 | +932 | 5,948 | 0.157 | FAIL (< 0.170) |
| X extend | 76,139 | 26,424 | 0.347 | +630 | 2,791 | 0.226 | PASS |
| S snap 10 | 68,240 | 24,196 | 0.355 | -1,598 | -5,108 | - | FAIL (loses chains) |
Misses left (0 / 1 / 2+ full-chain records; default 398 / 3,360 / 1,108): P 398 / 3,074 / 1,112; U 398 / 2,526 / 1,150;
X 276 / 3,138 / 1,125; S 390 / 3,124 / 2,776. P and X go to the A119b held-out run (bar 0.066, gained >= 100; default =
`r1617/default`).

## Held out on A119b (2026-11-10, `r1618/h_*`; default `r1617/default` 219,517 chains, 38,715 correct, 0.176)
| arm | chains | correct | precision | new chains | new correct | share (as registered) | lost | lost correct | net gained / net added |
|---|---|---|---|---|---|---|---|---|---|
| P compat | 229,437 | 39,379 | 0.172 | 17,116 | 782 | 0.046 | 7,196 | 118 | +664 / 9,920 = 0.067 |
| X extend | 236,301 | 40,086 | 0.170 | 25,894 | 1,666 | 0.064 | 9,110 | 295 | +1,371 / 16,784 = 0.082 |
- **Both FAIL as registered** (share correct among added chains < 0.066; X by 0.002). Gained >= 100 holds for both.
- Correction to the gorilla table: its "added" column was the NET change in chain count. The modes also drop chains (the
  new chains change locus fractions). Gorilla on the registered reading (new chains only): P 331 of 2,202 = 0.150 (FAIL
  vs 0.170), X 656 of 3,788 = 0.173 (PASS), U 1,160 of 12,351 = 0.094 (FAIL). On the net reading X passes on both samples
  and P passes on both narrowly; the net reading was not the registered one.
- Not adopted. X (`RUSTLE_EXTEND_TRANSFRAGS=1`) is the best opt-in trade: A119b +1,371 correct for -0.006 precision, gorilla
  +630 for -0.005.
