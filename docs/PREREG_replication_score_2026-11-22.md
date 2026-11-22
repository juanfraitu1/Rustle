# Row 1898: a one-read ranking fitted to replication instead of to the annotation (pre-registered 2026-11-22)

User goal (2026-11-22): "... better rules for the ones only covered by 1 long read or the junctions shared by 1 read."

State. The one-read tier is ranked by `oneread_score` (row 1806): 23 read features, among them the weakest junction's read
support (`min_isup`, `min_isup_rel`), the partial reads on the chain, its share of the locus, read-end concentration. So the
support of single-read junctions is already a feature. Row 1806's verdict was "not separable": on real libraries no fit
reaches the precision of two-read chains, and the label was the ANNOTATION, which cannot judge a real transcript seen once.
Rows 1890 / 1895: two hand-made orders are no better.

New evidence, looked at before this registration (`r1876_whats_left/repl.py`): WTC11 has three replicates of one sample.
Of the 52,756 one-read chains of replicate 2, 15,257 (29%) are seen again in replicate 1 or 3; those are 31.3% annotated, the
rest 5.3%. The frozen score hardly predicts it (best decile 41% seen again, worst 20%). Replication is a label that does not
depend on the annotation, and it is the same quantity row 1808's partner promotion uses with a second library in hand.

**Rule tested:** the same 23 features and the same design as `oneread_r1806.json`, refitted to "seen again in another
replicate" on WTC11 r2. The fitted ranking needs ONE library when applied; the replicates are training data only.

## Fit and data
- Features: `chain_rules_feat.py` over the floor-1 assemblies that exist (rows 1880, 1887), design of `oneread_tier.py`,
  standardised inside each library; dumped by `r1898_replication_score/feats.py`. No new assembly.
- Label on a library L of WTC11: a one-read chain is REPLICATED when its intron chain is in `x.gtf`, `final.gtf` or the
  one-read tier of either other replicate.
- **Model R:** logistic fit (`chain_rules.fit`) on r2's one-read chains, label = replicated.
- **Model A (control, same data, other label):** the same fit with label = annotated (RefSeq on CHM13). It tells whether a gain
  comes from the label or from refitting on this library.

## Measures
1. Inside WTC11, on r1 and r3 (not fitted): share replicated and share annotated in the best 10% and 30% of the tier, frozen
   score against R and A.
2. **Transfer, the bar:** on the datasets of row 1880 where StringTie is larger than `final.gtf` (17), correct chains of the
   ranked list at StringTie's size with tier 2 ordered by R, against the frozen order (`rank.py`'s count). Same for A.
3. The `list_sensitive` cut (best 10% of the tier): annotated chains under each order.

## Reading, fixed now (bars of rows 1890 / 1895)
- **ADOPTABLE** = R gives at least as many correct chains as the frozen order at StringTie's size on at least 15 of the 17
  datasets, on at least 3 of the 4 non-primate ones, and loses no AHEAD of row 1880. WTC11 r2 is fit data and is counted but
  marked; the verdict must also hold without it.
- **LABEL EFFECT** is claimed only if R is ADOPTABLE and A is not.
- If R is not adoptable the frozen ranking stays, and the reading is that single-library features cannot predict which
  one-read chains recur, whatever the label.
- An ADOPTABLE R is a candidate replacement of `tools/data/oneread_r1806.json` (the user's call); it changes only the order of
  the one-read tier, hence `PREFIX.ranked.gtf` and `PREFIX.oneread.ranked.gtf`, never `final.gtf`.

## Limits, stated in advance
- Replicated is not "real": an artefact of the aligner or of the library chemistry recurs too, and a real transcript at one
  molecule per library often does not (replicates 1 and 3 are larger than 2, which helps).
- One sample (human iPS), three replicates. Transfer to the other nineteen datasets is the only protection against fitting it.
- The judging truth on the transfer bar is the annotation again; a ranking that finds real unannotated chains is not rewarded.

## Result (2026-11-22): R is not ADOPTABLE (14 of 17, non-primate 2 of 4); no label effect; the frozen ranking is stale on primates

Feature dumps: jobs of `r1898_replication_score/jobids.txt`, 22 of 22 exit 0; `fit.py` -> `fit.txt`. Fit on WTC11 r2: 52,756
one-read chains, 15,257 replicated, 6,778 annotated.

**1. Inside WTC11, replicates not fitted** (best 10% of the tier; replicated / annotated):
| | frozen | R (replication label) | A (annotation label) | whole tier |
|---|---|---|---|---|
| r1 | 0.393 / 0.197 | 0.446 / 0.244 | 0.406 / 0.268 | 0.225 / 0.108 |
| r3 | 0.353 / 0.185 | 0.446 / 0.245 | 0.423 / 0.274 | 0.216 / 0.104 |
Each model is best on its own label; both beat the frozen score on both labels.

**2 and 3. Transfer** (correct chains of the list at StringTie's size, difference to StringTie in brackets; last column: annotated
chains in the best 10% of the tier, frozen | R | A):
| dataset | group | StringTie | frozen | R | A | best 10% |
|---|---|---|---|---|---|---|
| dme | N | (StringTie smaller than `final.gtf`) | | | | 180 | 60 | 174 |
| mmu | N | (StringTie smaller than `final.gtf`) | | | | 416 | 348 | 438 |
| dre | N | (StringTie smaller than `final.gtf`) | | | | 326 | 224 | 364 |
| es | N | 16401 | 17405 (+1004) | 17177 (+776) WORSE | 17284 (+883) WORSE | 442 | 176 | 303 |
| heart | N | 12577 | 13121 (+544) | 13177 (+600) ok | 13197 (+620) ok | 400 | 408 | 489 |
| adrenal | N | 13741 | 14388 (+647) | 14532 (+791) ok | 14562 (+821) ok | 404 | 505 | 571 |
| cortex | N | 14935 | 16207 (+1272) | 16084 (+1149) WORSE | 16117 (+1182) WORSE | 563 | 452 | 493 |
| wtc (fit data) | P | 20064 | 21134 (+1070) | 21616 (+1552) ok | 21803 (+1739) ok | 1052 | 1518 | 1703 |
| A119b | P | 31781 | 36480 (+4699) | 37211 (+5430) ok | 37972 (+6191) ok | 1153 | 1619 | 2417 |
| GGO_KB3781_fib_mm | P | 22087 | 23910 (+1823) | 24134 (+2047) ok | 24295 (+2208) ok | 667 | 914 | 1071 |
| HG002_kinnex.t2t | P | 27456 | 30147 (+2691) | 29824 (+2368) WORSE | 30825 (+3369) ok | 1026 | 724 | 1722 |
| OR6737 | P | 23248 | 25407 (+2159) | 25868 (+2620) ok | 25996 (+2748) ok | 1334 | 1891 | 1989 |
| PAB_AG06213_fib_mm | P | 21510 | 22819 (+1309) | 23122 (+1612) ok | 23298 (+1788) ok | 718 | 1128 | 1337 |
| PAB_mm | P | 16981 | 16295 (-686) | 17655 (+674) ok | 18149 (+1168) ok | 720 | 2039 | 2527 |
| PPA_mm | P | 17288 | 18644 (+1356) | 18913 (+1625) ok | 19137 (+1849) ok | 1245 | 1511 | 1726 |
| PPA_PR00251_fib_mm | P | 18780 | 20452 (+1672) | 20806 (+2026) ok | 21355 (+2575) ok | 1108 | 1293 | 1915 |
| PPY_AG05252_fib_mm | P | 20787 | 22164 (+1377) | 22752 (+1965) ok | 23075 (+2288) ok | 737 | 1306 | 1621 |
| PPY_mm | P | 20561 | 21227 (+666) | 22347 (+1786) ok | 22671 (+2110) ok | 1302 | 2238 | 2574 |
| PTR_AG18354_lcl_mm | P | 25433 | 27077 (+1644) | 27635 (+2202) ok | 27956 (+2523) ok | 948 | 1616 | 1949 |
| PTR_mm | P | 21532 | 22309 (+777) | 23184 (+1652) ok | 23529 (+1997) ok | 1849 | 2683 | 3023 |

- **R against the frozen order: at least as many correct chains on 14 of 17** (non-primate 2 of 4: heart and adrenal better,
  mouse ES -228 and cortex -123; primate 12 of 13, HG002 -323). Bar 15 of 17 and 3 of 4: **not ADOPTABLE.** Without the fit
  library 13 of 16.
- **A: 15 of 17** (non-primate 2 of 4, primate 13 of 13), and ahead of R on every dataset. So there is **no LABEL EFFECT**:
  the annotation label fits better than replication even when judged across libraries. Replication is not the missing piece.
- **What the control shows instead:** a fit on ONE library at today's defaults beats the frozen model on every primate
  dataset by 1-11% of the list's correct chains at StringTie's size (Sumatran testis 16,295 -> 17,655 with R, 18,149 with A:
  from 686 behind StringTie to 674 / 1,168 ahead), and loses on mouse ES, mouse cortex and (best 10%) fly simulated.
  `oneread_r1806.json` was fitted on 2026-11-18 on seven datasets under the defaults of that day (collapsed read counts,
  before the junction ratio, the spliced share and the count of every read); the read-count features it weights have changed
  scale since. A WTC11-only model is the opposite error: fitted to one human library.
- Not proposed: R or A as they are. **Next (row 1899):** the same design refitted at today's defaults on all datasets pooled,
  judged with each dataset held out of its own fit.
