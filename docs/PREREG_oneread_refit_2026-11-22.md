# Row 1899: the one-read ranking refitted at today's defaults, judged with each species held out (pre-registered 2026-11-22)

Row 1898 (seen before this registration): `tools/data/oneread_r1806.json` was fitted on 2026-11-18 under that day's defaults;
a refit on the single library WTC11 r2 at today's defaults beats it on all thirteen primate datasets at StringTie's size
(Sumatran testis from 686 behind StringTie to 1,168 ahead) and loses on mouse ES and mouse cortex. Neither the stale model nor
a one-library model is right. This row refits the same 23 features and the same design on all datasets at today's defaults
and judges the fit only on data of a species it did not see.

## Fit
- Features and design exactly as `oneread_tier.py` (dumps of row 1898: `r1898_replication_score/feat/NAME.npz`, one-read chains
  the default lacks, standardised inside each library). Label = intron chain in the dataset's truth (complete truth on the
  simulated sets; annotation elsewhere). Logistic fit `chain_rules.fit`.
- Balance: at most 20,000 one-read chains per dataset, drawn with a fixed seed (numpy `default_rng(1899)`).
- **Held-out models:** for each of nine species groups (fly; zebrafish; mouse: simulated + ES + heart + adrenal + cortex; human:
  WTC11 r2 + A119b + HG002; gorilla; chimpanzee; bonobo; Bornean orangutan; Sumatran orangutan) one model fitted on the other
  eight groups. A dataset is always scored by the model that never saw its species.
- **Final model** (only if the verdict allows): the fit on all twenty datasets.

## Measures (`r1899_oneread_refit/lodo.py`)
Per dataset, tier 2 ordered by the frozen score and by its held-out model:
- correct chains of the ranked list at StringTie's size (the 17 datasets where StringTie is larger than `final.gtf`);
- annotated chains in the best 10% of the tier (the `list_sensitive` cut), all 20;
- row 1880's AHEAD / BEHIND against every rival under the held-out order.

## Reading, fixed now
- **ADOPTABLE** = the held-out order gives at least the frozen order's correct chains at StringTie's size on at least 15 of 17
  datasets and at least 3 of the 4 non-primate ones among them; no comparison that row 1880 found AHEAD turns to BEHIND; and on
  the three simulated sets the best 10% of the tier holds at least 95% of the frozen order's correct chains on at least 2 of 3.
- ADOPTABLE makes the all-data model a candidate replacement of `oneread_r1806.json` (the user's call). It changes the order of
  `PREFIX.oneread.ranked.gtf` and of tier 2 of `PREFIX.ranked.gtf`, and with partners the promoted set; never `final.gtf`.
- Not ADOPTABLE: the frozen model stays; the table says which group fails.

## Limits, stated in advance
- Every one of the twenty datasets has been used in this project; "held out" means held out of THIS fit, by species.
- The label on real libraries is the annotation, with the non-human ape sets being largely computed (row 1883).
- A species group with one dataset (fly, zebrafish) is judged by a model fitted on vertebrates or on everything else.

## Result (2026-11-22): not ADOPTABLE (12 of 17; non-primate 1 of 4); primates gain, three mouse libraries lose; StringTie turns on Sumatran testis

`r1899_oneread_refit/lodo.txt` (`lodo.py`), feature dumps of row 1898. Each dataset scored by the model fitted without its species.
| dataset | group | species held out | StringTie correct | frozen order at its size | held-out model | difference | best 10% of the tier, frozen / held-out |
|---|---|---|---|---|---|---|---|
| dme | N | fly | StringTie smaller than `final.gtf` | | | | 180 | 179 ok |
| mmu | N | mouse | StringTie smaller than `final.gtf` | | | | 416 | 440 ok |
| dre | N | zebrafish | StringTie smaller than `final.gtf` | | | | 326 | 369 ok |
| es | N | mouse | 16401 | 17405 (+1004) | 17260 (+859) | -145 WORSE | 442 | 292 |
| heart | N | mouse | 12577 | 13121 (+544) | 13098 (+521) | -23 WORSE | 400 | 366 |
| adrenal | N | mouse | 13741 | 14388 (+647) | 14465 (+724) | +77 ok | 404 | 467 |
| cortex | N | mouse | 14935 | 16207 (+1272) | 16052 (+1117) | -155 WORSE | 563 | 403 |
| wtc | P | human | 20064 | 21134 (+1070) | 21193 (+1129) | +59 ok | 1052 | 1103 |
| A119b | P | human | 31781 | 36480 (+4699) | 37320 (+5539) | +840 ok | 1153 | 1754 |
| GGO_KB3781_fib_mm | P | gorilla | 22087 | 23910 (+1823) | 24077 (+1990) | +167 ok | 667 | 853 |
| HG002_kinnex.t2t | P | human | 27456 | 30147 (+2691) | 29759 (+2303) | -388 WORSE | 1026 | 622 |
| OR6737 | P | gorilla | 23248 | 25407 (+2159) | 25346 (+2098) | -61 WORSE | 1334 | 1233 |
| PAB_AG06213_fib_mm | P | Sumatran | 21510 | 22819 (+1309) | 23088 (+1578) | +269 ok | 718 | 1076 |
| PAB_mm | P | Sumatran | 16981 | 16295 (-686) | 17531 (+550) | +1236 ok | 720 | 1933 |
| PPA_mm | P | bonobo | 17288 | 18644 (+1356) | 18834 (+1546) | +190 ok | 1245 | 1434 |
| PPA_PR00251_fib_mm | P | bonobo | 18780 | 20452 (+1672) | 20907 (+2127) | +455 ok | 1108 | 1546 |
| PPY_AG05252_fib_mm | P | Bornean | 20787 | 22164 (+1377) | 22588 (+1801) | +424 ok | 737 | 1129 |
| PPY_mm | P | Bornean | 20561 | 21227 (+666) | 21942 (+1381) | +715 ok | 1302 | 1914 |
| PTR_AG18354_lcl_mm | P | chimpanzee | 25433 | 27077 (+1644) | 27539 (+2106) | +462 ok | 948 | 1496 |
| PTR_mm | P | chimpanzee | 21532 | 22309 (+777) | 22841 (+1309) | +532 ok | 1849 | 2346 |

- **At StringTie's size the held-out order is at least the frozen one on 12 of 17:** primate 11 of 13 (HG002 -388, OR6737 -61),
  non-primate 1 of 4 (mouse ES -145, heart -23, cortex -155; adrenal +77). Bar 15 of 17 and 3 of 4: **not ADOPTABLE.**
- Simulated sets, best 10% of the tier: 3 of 3 at or above 95% of the frozen order (fly 180 -> 179, mouse 416 -> 440,
  zebrafish 326 -> 369). No comparison that row 1880 found AHEAD turns to BEHIND.
- **Sumatran testis against StringTie turns to AHEAD: -686 -> +550**, with a model that saw no orangutan library. Under the
  held-out order the list would be ahead of StringTie on 20 of 20 datasets.
- The sensitive cut (best 10%) moves the same way and harder: +10% to +168% annotated chains on eleven primate datasets,
  -34% on mouse ES, -28% on cortex, -39% on HG002.
- The failure is the mouse group: with every mouse dataset out of the fit, the model is fitted mostly to ape and human libraries
  and ranks mouse libraries worse than the frozen model, which was itself fitted on mouse ES and heart. That is the direction
  of a primate fit, and the registered bar exists to catch it.

**For context, not a test (`insample.py`):** the model fitted on all twenty datasets (`r1899_oneread_refit/oneread_r1899.json`)
is at least the frozen order on 15 of 17 at StringTie's size (mouse ES -59, HG002 -193; heart +61, adrenal +135, cortex +97;
Sumatran testis -686 -> +749) and on 18 of 20 in the best 10%. Every dataset is in its own fit there, as seven of them are in
the frozen model's.

**Reading (mine).** The frozen ranking is stale and a refit at today's defaults is clearly better on primates and on the
simulated sets; whether it is better on a species it has not seen depends on the species. The honest state: no adoptable
replacement by the bar set in advance. Candidate, the user's call and against that bar: replace `oneread_r1806.json` by the
all-data model. It would close the last StringTie comparison and it is not validated on unseen species.

## Descriptive addendum (2026-11-22, not registered): one-read chains by the read support of their weakest junction

The user's third question ("junctions shared by 1 read"), shown directly; the quantity is already the feature `min_isup` of the
ranking. `r1899_oneread_refit/junc1.NAME.txt` (job `junc1`), one-read chains the default lacks, annotated / all (share):
| dataset | this read only | 2 reads | 3-5 reads | 6+ reads | 6+ reads and no other read on the chain | whole tier |
|---|---|---|---|---|---|---|
| mouse simulated | 3,446 / 4,117 (0.837) | 309 / 462 (0.669) | 175 / 338 (0.518) | 88 / 173 (0.509) | 53 / 111 (0.477) | 4,018 / 5,090 (0.789) |
| zebrafish simulated | 2,545 / 3,427 (0.743) | 216 / 390 (0.554) | 107 / 296 (0.361) | 25 / 93 (0.269) | 14 / 69 (0.203) | 2,893 / 4,206 (0.688) |
| mouse heart | 2,087 / 16,404 (0.127) | 353 / 2,783 (0.127) | 324 / 2,907 (0.111) | 133 / 2,295 (0.058) | 54 / 703 (0.077) | 2,897 / 24,389 (0.119) |
| mouse cortex | 1,837 / 22,763 (0.081) | 425 / 4,512 (0.094) | 390 / 4,821 (0.081) | 256 / 4,464 (0.057) | 75 / 1,016 (0.074) | 2,908 / 36,560 (0.080) |
| WTC11 r2 | 3,797 / 32,830 (0.116) | 1,208 / 7,112 (0.170) | 1,228 / 7,430 (0.165) | 545 / 5,384 (0.101) | 62 / 1,329 (0.047) | 6,778 / 52,756 (0.128) |
| chimpanzee testis | 4,376 / 73,043 (0.060) | 1,888 / 17,598 (0.107) | 2,076 / 19,962 (0.104) | 1,676 / 21,944 (0.076) | 84 / 4,115 (0.020) | 10,016 / 132,547 (0.076) |
| Sumatran testis | 3,275 / 98,751 (0.033) | 1,400 / 22,756 (0.062) | 1,442 / 25,030 (0.058) | 895 / 24,976 (0.036) | 49 / 4,966 (0.010) | 7,012 / 171,513 (0.041) |

- **A junction carried by one read is not a mark of error.** On complete truth (mouse, zebrafish simulated) the one-read
  chains whose weakest junction only that read carries are the MOST often true (0.84 on mouse simulated against 0.51-0.67
  for chains whose junctions other reads carry): they are whole transcripts seen once. A rule that drops single-read
  junctions would remove the best class there.
- On real libraries the relation is not monotone: support of 2-5 reads is best (1.3-1.9 times the single-read class on the
  three primate libraries, no difference on mouse heart), and 6 or more is the WORST (a one-read chain made of heavily used
  junctions in a busy locus is a rare recombination or an error); with no other read on the chain that cell is 1-8%
  annotated. No cell exceeds 0.19 on a real library.
- So no cut on junction support is a rule; what there is to gain from it the ranking takes, and the refit of this row is
  the better use of the same information.
