# PREREG: family recall by age of the duplication, per tool (row 1857, 2026-11-21)

User (2026-11-21, after row 1855): the youngest families have near-identical members, their reads are the near-tie
multimappers, StringTie and the other tools discard secondary alignments, so they should lose those members and Rustle
should find more of them. Row 1855 scored all ages together and cannot show this. Written before any stratified score.

## Measurement
Row 1855's family tables, unchanged (`r1847_family_benchmark/a1/ARM/x.fam.copies.tsv`, A119b, fixed truth set
`union.gtf`). `pairscore.py --level L` per stratum of the Compara duplication node: `Homo sapiens`;
`Homininae,Hominidae,Hominoidea`; `Catarrhini`; `Simiiformes`. Per arm and stratum: true pairs found / true pairs (recall).
Precision is not stratified (a predicted pair has no age). OR6737 is reported beside with the same strata if its protein-hit
truth holds >= 30 pairs in the youngest stratum; no bar there.

## Rule (A119b, arm `rustle` against `stringtie`, `isoquant`, `bambu`, `flair`)
Lead(stratum, tool) = Rustle recall - tool recall.
- SUPPORTED: in the `Homo sapiens` stratum Rustle's recall is the highest of the five, AND its lead over StringTie there is
  >= 0.05 larger than its lead over StringTie in the `Simiiformes` stratum.
- NOT SUPPORTED: otherwise. A stratum with < 30 true pairs is reported and not judged.

## Predictions (mine)
Not supported, or weakly: the tools read the primary alignment, and a near-tie read's primary still lands on one of the
copies, so each copy keeps part of its reads and its transcript; the tools lose depth, not the locus. Rustle's lead on A119b
comes from all ages.

## Limits
One sample judged. Age is the annotation's, not the pair's sequence identity; gene conversion keeps some old pairs
near-identical. The judged Rustle arm seeds loci from primaries too; this measures the outcome, not the mechanism.

## Result of row 1857 (2026-11-21): NOT SUPPORTED
A119b, recall per stratum (true pairs found / true pairs), job 4491242, tables `r1857_family_by_age/a1.{1..4}.tsv`:

| stratum (true pairs) | rustle | stringtie | flair | isoquant | bambu |
|---|---|---|---|---|---|
| Homo sapiens (780) | 0.588 (459) | 0.658 (513) | 0.517 (403) | 0.414 (323) | 0.392 (306) |
| Homininae to Hominoidea (465) | 0.686 (319) | 0.456 (212) | 0.641 (298) | 0.391 (182) | 0.185 (86) |
| Catarrhini (239) | 0.569 (136) | 0.243 (58) | 0.548 (131) | 0.385 (92) | 0.167 (40) |
| Simiiformes (458) | 0.231 (106) | 0.103 (47) | 0.188 (86) | 0.114 (52) | 0.105 (48) |

Rustle is not the highest in the youngest stratum (StringTie leads by 0.070, 54 pairs); its lead over StringTie is -0.070
there against +0.128 in the oldest. Rustle's lead of row 1855 comes from the three older strata. OR6737 (protein-hit truth;
93 / 71 / 108 / 217 pairs): youngest stratum flair 0.731, stringtie 0.677, rustle 0.624, isoquant 0.559, bambu 0.376; rustle
below stringtie in the three younger strata, level in the oldest (0.129 / 0.124). The strata sum to 1,942 on A119b, 6 more
than row 1855's 1,936 (pairs annotated at two nodes).

# PREREG: Rustle with and without secondary alignments (row 1858, 2026-11-21)

User (2026-11-21): "one run with rustle's secondaries totally off". Found while setting it up, before any score: row 1855's
judged Rustle arm was ALREADY assembled from primary alignments only (`ours.log`: "seeding pool: primary"; the default of
`assemble` since row 1812), and its families stage reads no alignments (`RUSTLE_FAMREPS=0`). Row 1855 therefore compared
tools that all used primaries only. The contrast that tests the user's idea is the reverse one: secondaries ON.

## Arms (today's binaries and defaults, A119b and OR6737, nothing else changed)
- `primary`: `assemble --seed-pool primary` (secondaries totally off).
- `good`: `assemble --seed-pool good` (primaries + secondaries within 2% of the molecule's best alignment score).
Then `families` with `RUSTLE_FAMREPS=0` on each, scored as rows 1855 / 1857 (fixed truth set `r1847_family_benchmark/S/union.gtf`,
all ages and the four age strata).

## Rule (A119b)
- SECONDARIES HELP: `good` - `primary` recall in the `Homo sapiens` stratum >= +0.03, and overall precision not lower by
  more than 0.02.
- NO HELP: otherwise. OR6737 reported beside, no bar.

## Prediction (mine)
A small gain, below the bar: a near-tie read's primary already lands on one of the copies.

## Result of row 1858 (2026-11-21): the recall half holds strongly, the precision half fails
Jobs 4491257 / 4491258, tables `r1858_seed_pool_families/{a1,or}/score.{0..4}.tsv` (commit in `commit.txt`).

| A119b | primary | good |
|---|---|---|
| all ages: predicted pairs, precision / recall | 1,183, 0.841 / 0.514 | 1,903, 0.650 / 0.639 |
| Homo sapiens (780) recall | 0.588 (459) | 0.853 (665) |
| Homininae to Hominoidea (465) | 0.667 | 0.733 |
| Catarrhini (239) | 0.540 | 0.573 |
| Simiiformes (458) | 0.223 | 0.216 |

Registered verdict: NO HELP, because overall precision falls by 0.191 (bar 0.02). The recall gain in the youngest stratum is
+0.265 (bar +0.03), 206 more true pairs, and passes StringTie's 0.658 of row 1857; the gain shrinks with age and is nil in the
oldest stratum, as the user's mechanism predicts. My prediction (a small gain below the bar) was wrong on recall.
OR6737: all ages 598 pairs 0.271 / 0.333 -> 2,135 pairs 0.090 / 0.394; youngest stratum (93) 0.699 -> 0.720.
Transcripts: A119b 213,397 -> 225,945; OR6737 68,822 -> 70,155.
