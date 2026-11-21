# Row 1882: does transfrag extension close the short-read testis miss of row 1880? (pre-registered 2026-11-21)

Row 1880: Rustle's ranked list is behind StringTie on one dataset of twenty, Sumatran orangutan testis (`PAB_mm`, -686 correct
chains at StringTie's size, behind along StringTie's whole ranking), and behind at StringTie's top quarter on Bornean testis
(`PPY_mm`). Both are the short-read testis libraries (memory: shorter reads, about three secondary records per read). StringTie
builds chains from partial reads; Rustle's answer to that is transfrag extension, off by default since row 1811 because its
gain was confined to two ape testis libraries (`RUSTLE_EXTEND_TRANSFRAGS=1` turns it on).

## Runs
The row 1880 job (`assemble` + `oneread`, same driver snapshot and binaries) with `RUSTLE_EXTEND_TRANSFRAGS=1`, on `PAB_mm`
and `PPY_mm`, and on two controls where the list is already ahead: chimpanzee testis (`PTR_mm`) and mouse cortex. Workspace
`r1882_extension_short_reads`; scorer `r1880_ranked_vs_rivals/rank.py` unchanged.

## Reading, fixed now
- **CLOSES** = with extension the list is AHEAD of StringTie on `PAB_mm` at StringTie's size and at its three inner points, and
  no rival that row 1880 found outperformed on the four datasets turns to BEHIND.
- **PARTLY** = the `PAB_mm` deficit at StringTie's size shrinks by at least half, no rival turns.
- Otherwise the miss stands as a property of the library and nothing is proposed.
- A CLOSES or PARTLY is a lead only: extension stays off by default (row 1811, user), and a switch that turns it on from a
  library statistic would need its own row with more libraries.

## Result (2026-11-21): PARTLY; and with extension `final.gtf` itself gains on both axes on the three testis libraries

Jobs 4494676-4494679, every stage exit 0; `r1882_extension_short_reads/score.NAME.txt`.
| dataset | `final.gtf` off | `final.gtf` with extension | StringTie at its size, off -> on | top 25% | top 50% | top 75% |
|---|---|---|---|---|---|---|
| Sumatran testis | 54,748 / 15,552 / 0.2841 | 56,157 / 16,327 / 0.2907 | -686 -> +450 | -673 -> -610 | -460 -> -162 | -62 -> +474 |
| Bornean testis | 83,473 / 19,525 / 0.2339 | 85,484 / 20,478 / 0.2396 | +666 -> +2,153 | -235 -> -47 | +326 -> +854 | +842 -> +1,901 |
| chimpanzee testis | 46,969 / 20,309 / 0.4324 | 47,756 / 21,106 / 0.4420 | +778 -> +1,865 | -99 -> +26 | +313 -> +663 | +706 -> +1,284 |
| mouse cortex | 29,767 / 15,490 / 0.5204 | 30,023 / 15,490 / 0.5159 | +1,272 -> +1,121 | +349 -> +351 | +664 -> +655 | +882 -> +845 |

- **PARTLY as registered:** on Sumatran testis the deficit at StringTie's size turns into a lead (-686 -> +450) and no rival
  turns to BEHIND on any of the four; the list is still behind StringTie's top quarter and top half there (-610, -162), so not CLOSES.
- **Not asked by this row, and the larger finding:** with extension `final.gtf` has more correct chains AND higher precision on
  all three testis libraries (+775 / +953 / +797 correct chains, +0.007 / +0.006 / +0.010 precision). Row 1811 called
  extension a recall rule that lowers precision; that was `x.gtf` under the defaults of 2026-11-19. On mouse cortex it adds 256
  chains, none annotated (precision -0.0045).
- Follow-up registered as row 1883: extension at today's defaults on all twenty datasets.
