# Which assembler levers are overfit to the primate samples, and what was toned down (18 Nov 2026)

One page. Every threshold of the assembler and both chain scores were set on 12 primate samples. Rows 1779-1810 tested them on
data they were never tuned on. Details and full tables: the pre-registration named per row.

## Current state
- **Toned down (user, 2026-11-18, row 1810):** sub-chain collapse ratio 0.7 -> **0.3**; chain score cut on `PREFIX.final.gtf`
  0.055 -> **0.02**. `RUSTLE_POLISH_ISM_RATIO=0.7` and `RUSTLE_CHAINSCORE_TAU=0.055` give the earlier commands.
- **New default (row 1800):** polish junction snap (`RUSTLE_POLISH_JUNCTION_SNAP` unset = `reads`).
- **Switched off (user, 2026-11-19, row 1812):** transfrag extension (opt-in: `RUSTLE_EXTEND_TRANSFRAGS=1`); the stage `assemble`
  on its own seeds with primary alignments (`--seed-pool good` = before; `all` and the family stages keep the good pool).
- **Filter replaced (user, 2026-11-19, row 1816):** `final.gtf` is filtered by chain score r1813, whose score for a chain depends on
  that chain's reads only (r1769 changed 2-20% of its decisions when the chain set shifted, row 1813), and by the 5' truncation
  rule (row 1815: drop a chain when reads splicing on past its 5' end outnumber the reads matching it by more than 3.22; found
  on two apes, passes a third ape and 5 of 6 non-primate datasets). `RUSTLE_CHAINSCORE_MODEL=r1769`, `RUSTLE_TRUNC5_RATIO=0` = before.
- **Not changed:** every lever listed as general or as a recall trade below.
- Base runs and figures made before 2026-11-18 do not carry these changes.

## The levers
Reading of row 1793: a lever *pays* on a dataset when it raises precision and costs at most 5% of the correct chains (10% for the
three levers meant as large cuts). OVERFIT = pays on both tuning samples (OR6737, A119b), fails on at least 2 of 3 new datasets
(fly simulated, mouse simulated, mouse ES real).

| lever | switch | verdict | evidence | action |
|---|---|---|---|---|
| chain score filter, cut | `tools/data/chainscore_r1769.json` `tau`, `RUSTLE_CHAINSCORE_TAU` | **OVERFIT** | cost 15.3% / 12.9% of correct chains on fly / mouse simulated, 7.8% mouse ES; 1.0% / 4.2% on the tuning samples. It rewards chains for carrying truncated reads; simulated reads have fewer | cut 0.055 -> 0.02 |
| sub-chain collapse ratio ("ISM ratio") | `RUSTLE_POLISH_ISM_RATIO` | **OVERFIT** | cost 9.8% mouse simulated, 9.2% mouse ES; 1.4% / 3.1% on the tuning samples | 0.7 -> 0.3 |
| isoform fraction 0.02 | `--polish-isoform-fraction` | general | pays on all five (cost 2.2-7.1%) | none |
| strict junction gate | `--assembly-junctions strict` | general | pays on all five (cost <= 0.6%) | none |
| retained-intron ratio 10 | `RUSTLE_POLISH_RETAINED_RATIO` | general | pays on all five (cost <= 2.4%) | none |
| sub-chain drop on `final.gtf` | `RUSTLE_FINAL_SUBCHAIN` | general | pays on all five (cost <= 2.1%) | none |
| read floor 2 | `RUSTLE_MIN_CHAIN_READS` | mixed | costs 10.3% on OR6737, 5.5-6.4% elsewhere | none; one-read chains go to separate tiers (below) |
| transfrag extension | `RUSTLE_EXTEND_TRANSFRAGS` | recall trade | row 1811: **PRIMATE-ONLY**, +3.6% / +2.5% correct chains on A119b / OR6737, under 0.4% on five non-primate datasets, 7+ wrong chains per correct one | default off (row 1812) |
| seed pool | `RUSTLE_SEED_POOL` | recall trade | row 1811: under 1% correct chains on all seven datasets (A119b +276 correct, +12,446 wrong) | `assemble` alone: primary (row 1812) |
| chain score r1769 over r1677 | `RUSTLE_CHAINSCORE_MODEL` | mixed | fails OR6737 and fly | none |
| TSS rescue | `RUSTLE_POLISH_TSS` | no effect off A119b | changes `PREFIX.gtf` on A119b only | none |
| polish junction snap | `RUSTLE_POLISH_JUNCTION_SNAP` | was opt-in | raises precision on all five for < 0.2% of correct chains; passed on zebrafish simulated and mouse heart | default on (row 1800) |

The sub-chain collapse uses no annotation: a transcript whose intron chain is a contiguous part of a longer assembled chain is
kept when it has at least R times that chain's reads. Only the value of R was fitted to primates.

## What the toned-down values cost and buy
Sub-chain ratio, snap on, `PREFIX.gtf`, correct chains / precision (row 1800 grid):

| ratio | fly sim | mouse sim | mouse ES | OR6737 | A119b |
|---|---|---|---|---|---|
| 0.7 (before) | 7,053 / 0.855 | 27,144 / 0.894 | 17,708 / 0.325 | 26,419 / 0.350 | 40,078 / 0.170 |
| **0.3 (now)** | 7,089 / 0.836 | 27,777 / 0.869 | 18,132 / 0.316 | 26,597 / 0.314 | 40,772 / 0.147 |
| 0.1 (non-primate choice) | 7,136 / 0.819 | 28,732 / 0.844 | 18,700 / 0.296 | 26,718 / 0.272 | 41,183 / 0.124 |

At 0.3 the ratio still costs 7.6% of ratio 0's correct chains on mouse simulated (limit 5%; only 0.1 is inside it).

Chain score cut, share of correct chains kept / precision of the kept chains (snap on, ratio 0.7 runs; `tau.py`, row 1810):

| cut | fly sim | mouse sim | zebrafish sim | mouse ES | mouse heart | OR6737 | A119b |
|---|---|---|---|---|---|---|---|
| 0.055 (before) | 0.820 / 0.978 | 0.835 / 0.950 | 0.714 / 0.943 | 0.917 / 0.415 | 0.899 / 0.662 | 0.989 / 0.377 | 0.955 / 0.243 |
| **0.02 (now)** | 0.937 / 0.962 | 0.939 / 0.923 | 0.892 / 0.909 | 0.976 / 0.355 | 0.968 / 0.599 | 0.998 / 0.358 | 0.988 / 0.200 |
| no filter | 1 / 0.855 | 1 / 0.894 | 1 / 0.794 | 1 / 0.325 | 1 / 0.557 | 1 / 0.350 | 1 / 0.170 |

Both changes move recall up and precision down on every real library. The two together were not measured as one run on all
datasets; the cut table was made at ratio 0.7.

**Problem found with row 1812, repaired by row 1816: the chain score filter r1769 is unstable.** With extension off and primary seeds it drops 16,480
chains on mouse ES instead of 5,582 (1,405 correct chains lost from `final.gtf`) and 6,729 instead of 3,518 on OR6737, though
`PREFIX.gtf` barely changed. Likely cause: its decile thresholds are taken within the sample and sit on tied feature values.

## One-read chains (rows 1805-1808)
The floor stays at two reads. Opt-in stage `oneread` writes: a strict tier (`PREFIX.oneread.gtf`, row 1766); every one-read chain
ranked by a score from the reads (`PREFIX.oneread.ranked.gtf`; a ranking, not a filter); and, with
`RUSTLE_ONEREAD_PARTNERS=A.bam:B.bam`, promoted chains (`PREFIX.oneread.promoted.gtf`: seen in a second library of the species
and in the best 20% of the ranking; as precise as the default's two-read chains or better on three unused libraries, 2.3-3.1%
more correct chains).

## Limits
Simulated reads lack degradation and nascent RNA, so their precision flatters every tool. Annotations are incomplete, so "wrong"
on real libraries includes real chains. Two tuning samples stood for the twelve. Levers were switched one at a time.

Documents: `PREREG_simulation_transfer_2026-11-17.md` (rows 1779 / 1780), `PREREG_rule_ablation_2026-11-18.md` (1793),
`PREREG_snap_ism_defaults_2026-11-18.md` (1800), `PREREG_one_read_chains_2026-11-18.md` (1805),
`PREREG_one_read_score_2026-11-18.md` (1806), `PREREG_one_read_promote_2026-11-18.md` (1807),
`PREREG_one_read_promote_check_2026-11-18.md` (1808). Workspaces `r1779_simulated_reads`, `r1793_rule_ablation`,
`r1800_snap_ism`, `r1805_one_read` ... `r1810_tone_down` under `/anvme/workspace/v138da13-rustle`.

**Verification (job 4482795, fly simulated).** A run at the new defaults reproduces the tested ratio-0.3 arm (`PREFIX.gtf`, same
md5) and filters at 0.02; with `RUSTLE_POLISH_ISM_RATIO=0.7 RUSTLE_CHAINSCORE_TAU=0.055` both outputs reproduce the earlier run
(same md5). Both changes together on fly: `final.gtf` 5,768 correct / 0.980 before, 6,748 / 0.967 now (truth 9,379).
