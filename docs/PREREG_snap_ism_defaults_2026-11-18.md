# Row 1800: polish junction snap on by default, ISM ratio re-set on non-primate data (pre-registered 2026-11-18)

User (2026-11-18), after row 1793 (`docs/PREREG_rule_ablation_2026-11-18.md`): do proposals 1 (polish junction snap on) and 2
(re-set the ISM ratio using non-primate data). Registered before any arm below ran. Workspace `r1800_snap_ism`.

**Unit.** Exact intron chains of multi-exon transcripts (`bench/famcmp/curve.py`), ours de novo, `RUSTLE_PTC=0` as in row 1793.
"Old defaults" = the `base` arm of row 1793 (same binaries, built 2026-11-18 04:20).

## Step 1, choose the ISM ratio (diagnostic data of row 1793)

Arms: polish junction snap `reads` on, `--polish-ism-ratio` 0.7, 0.5, 0.4, 0.3, 0.2, 0.1, 0, on fly simulated, mouse simulated,
mouse ES, OR6737, A119b. Read on `PREFIX.gtf`.
A ratio *pays* on a dataset as in row 1793: against ratio 0 it has higher precision and costs at most 5% of ratio 0's correct chains.
**Choice, fixed now:** r* = the largest ratio of the grid that pays on all three non-primate datasets. If 0.7 does so with the snap
on, the ratio stays and proposal 2 is void. If no ratio does, r* = 0 is not adopted: reported, the ratio stays, the user decides.
r* must also keep precision above ratio 0 on OR6737 and A119b; if not, reported and the user decides.

## Step 2, check on data never used for a choice

- **U1, zebrafish simulated:** GRCz11 + Ensembl 112, the row 1779 recipe unchanged (IsoSeqSim call of LRGASP, mouse per-transcript
  counts drawn onto the transcripts, seed 1800), truth = transcripts with >= 1 simulated read.
- **U2, mouse tissue, real:** one ENCODE long-read RNA-seq PacBio library of an adult mouse tissue (not ES cells, accession fixed
  in the Result before mapping), truth GENCODE vM27 (incomplete).
Mapping and assembly as in rows 1779 / 1793. Arms: old defaults; old + snap; snap + r*; snap + ratio 0.

Bars, per dataset:
- **C1, snap.** old + snap against old defaults, `PREFIX.gtf`: precision higher and at most 2% of the old correct chains lost.
- **C2, ratio.** snap + r* against snap + 0.7: more correct chains; and r* pays against snap + ratio 0.
- **C3, the two together under the filter.** `PREFIX.final.gtf` of snap + r* is not dominated by `PREFIX.final.gtf` of the old
  defaults (not both fewer correct chains and lower precision).

**Decision.** Snap becomes the driver default when C1 passes on U1 and U2. r* becomes the default when C2 passes on both.
C3 failing on either stops both changes: reported, the user decides. Gained and lost chains are reported separately for every
comparison. Defaults change in the working tree only; commit on the user's word.

**Limits, stated in advance.** U2 is the same species as one diagnostic dataset (another tissue and library). U1 is simulated, so
it lacks truncated reads (row 1793) and flatters precision. The chain score filter was fitted without the snap; row 1793 saw the
snap cost 3.3% of correct chains in `final.gtf` on mouse simulated. Refitting the filter is not part of this row. The 12 primate
samples are represented by two. Known wall-time cost of the snap: A119b assemble 680 s -> 872 s.

**Fixed before any mapping (2026-11-18).** U2 = ENCODE ENCSR938RPH, adult (2 months) B6CASTF1/J mouse heart, PacBio Sequel II,
Mortazavi lab, replicates ENCFF836SSV + ENCFF031EVU pooled (1,764,827 reads; the F1 hybrid carries a CAST haplotype that differs
from GRCm39). First looked at: ENCSR491KJF (C57 heart, Guigo lab, CapTrap), set aside because its library is unstranded and the
shipped mapping (`-uf`) needs oriented reads. Jobs: preparation 4482464, mapping 4482465 (zebrafish), 4482466 (heart); step 1
grid 4482459-4482463. The simulated FASTA is not kept (written to the job's temporary directory).

## Result (2026-11-18): r* = 0.1; C1, C2, C3 pass on both unused datasets; snap made the default, ratio held for the user

All arms exit 0 (grid 35, check 8). Tables: workspace `r1800_snap_ism/score.*.txt`. Zebrafish: 3,476,946 simulated reads, 19,955
truth transcripts.

**Step 1.** Snap on, `PREFIX.gtf`, correct chains / precision; cost = share of ratio 0's correct chains.

| ratio | fly sim | mouse sim | mouse ES | OR6737 | A119b |
|---|---|---|---|---|---|
| 0 | 7,154 / 0.809 | 30,074 / 0.824 | 19,492 / 0.234 | 26,787 / 0.220 | 41,354 / 0.109 |
| 0.1 | 7,136 / 0.819 (0.3%) | 28,732 / 0.844 (4.5%) | 18,700 / 0.296 (4.1%) | 26,718 / 0.272 | 41,183 / 0.124 |
| 0.2 | 7,119 / 0.828 | 28,138 / 0.858 (**6.4%**) | 18,338 / 0.309 (**5.9%**) | 26,656 / 0.297 | 40,998 / 0.137 |
| 0.3 | 7,089 / 0.836 | 27,777 / 0.869 | 18,132 / 0.316 | 26,597 / 0.314 | 40,772 / 0.147 |
| 0.4 | 7,079 / 0.843 | 27,562 / 0.876 | 17,983 / 0.318 | 26,558 / 0.324 | 40,596 / 0.153 |
| 0.5 | 7,072 / 0.847 | 27,406 / 0.883 | 17,869 / 0.320 | 26,511 / 0.332 | 40,450 / 0.159 |
| 0.7 | 7,053 / 0.855 | 27,144 / 0.894 (9.7%) | 17,708 / 0.325 (9.2%) | 26,419 / 0.350 | 40,078 / 0.170 |
| old defaults | 7,058 / 0.833 | 27,175 / 0.821 | 17,714 / 0.320 | 26,424 / 0.347 | 40,086 / 0.170 |

r* = 0.1: the only ratio within 5% on all three non-primate datasets; precision above ratio 0 on both tuning samples.

**Step 2.** Correct chains / precision.

| arm | zebrafish `gtf` | zebrafish `final.gtf` | heart `gtf` | heart `final.gtf` |
|---|---|---|---|---|
| old defaults | 12,475 / 0.738 | 9,154 / 0.939 | 13,101 / 0.553 | 11,799 / 0.664 |
| old + snap | 12,472 / 0.793 | 8,841 / 0.944 | 13,099 / 0.557 | 11,763 / 0.667 |
| snap + 0.1 | 12,753 / 0.725 | 9,758 / 0.930 | 13,611 / 0.499 | 12,435 / 0.626 |
| snap + 0 | 12,845 / 0.705 | 10,109 / 0.921 | 13,876 / 0.442 | 12,965 / 0.585 |

- **C1 PASS, both.** Snap against old: zebrafish correct +4 / -7, wrong +10 / -1,192; heart +0 / -2, +9 / -178.
- **C2 PASS, both.** 0.1 against 0.7 (snap on): zebrafish correct +282 / -1, wrong +1,598 / -2; heart +513 / -1, +3,217 / -3.
  0.1 against 0: zebrafish correct +0 / -92 (0.7%), wrong -538; heart +1 / -266 (1.9%), wrong +5 / -3,882.
- **C3 PASS, both** (more correct chains, lower precision: not dominated). `final.gtf`, snap + 0.1 against old: zebrafish
  correct +737 / -133, wrong +184 / -43; heart +645 / -9, +1,535 / -70.
- Snap alone under the filter (`final.gtf`, not a bar): zebrafish correct +71 / -384 (3.4% lost, as on mouse simulated in row
  1793); heart +6 / -42.

**Done.** `tools/rustle_pipeline.sh`: `RUSTLE_POLISH_JUNCTION_SNAP` unset = `reads` (`off` restores the old command).
**Not done, held for the user:** `--polish-ism-ratio` stays 0.7. The bars pass, but the bars only ask for more correct chains
than 0.7 and a gain over ratio 0; against the old defaults 0.1 lowers `PREFIX.gtf` precision on every real library (heart 0.553
-> 0.499, mouse ES 0.320 -> 0.296, OR6737 0.347 -> 0.272, A119b 0.170 -> 0.124: A119b gains 1,097 correct chains for about
96,000 chains RefSeq lacks). It is a recall-for-precision trade, the same on primates and non-primates.

**User ruling (2026-11-18).** `--polish-ism-ratio` stays 0.7. It is marked a candidate for change, to be taken up when precision
is next tuned: the measured non-primate choice is 0.1, and the grid above is the starting table (0.3 as a middle value).
