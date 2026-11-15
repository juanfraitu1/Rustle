# PREREG row 1669: single-exon feed as an opt-in `families` switch, unseen-sample test (2026-11-13)

Families layer (O1). User request ("re-register it") after row 1665 failed only on a recall bar whose denominator moves.

## What is wired (before the test, opt-in, default unchanged)
`RUSTLE_MONO_FEED=span` in `tools/rustle_pipeline.sh` `families`: input = `$FAM_GTF` plus one single-exon model per pile
of `PREFIX.final.gtf` cluster models that no exon of it overlaps on the same strand (`tools/mono_feed.py`, moved from
bench), written to `PREFIX.fam.monofeed.gtf`. Unset = the same command as before. `flag` still reads the unfed loci.
Driver test `families_feeds_single_exon_cluster_models_only_when_asked` (watched fail, then pass).

Arm choice is post hoc: `span` was picked after row 1665 showed `span` 497 vs `most` 495 placed on A119b. Row 1665's
A119b figures on a fixed set (hits 1,768 -> 1,780, 0 lost, false 192 -> 194) are therefore not a test.

## Test: HG002 (not seen by this rule), frozen before the run
`assemble` (default, layers on) once on `~/work/jfim/hg002_run/hg002.pat.bam`, then `families` with the switch unset and
with `span`, same binaries and driver copy. Truth: `--species hg002`. Fixed truth set = pairscore hit dumps and
false-pair counts (no ratio with a moving denominator).
Validity: baseline loci all present in the arm; added loci reach `fam.loci.gff3`.
Bars: (1) hits >= baseline hits and true pairs lost <= 5; (2) new false pairs <= new true pairs;
(3) members placed (member_funnel) >= baseline + 1.
PASS = the switch stays as a documented opt-in. FAIL = the switch is removed from the driver and the row closes.
Gorilla KB3781 was already run in row 1665: no gain, no harm.

## Result (2026-11-13, job 4468894, `/anvme/workspace/v138da13-rustle/r1669_monofeed_hg`)

HG002, `span`: piles 1,067, added 891. Validity holds (baseline loci missing 0; added loci in gff3 1,782; added copies 29).

| | pred pairs | hit | false | precision | recall | placed | wrong_family | no_family | no_model |
|---|---|---|---|---|---|---|---|---|---|
| base | 256 | 184 | 72 | 0.719 | 0.448 | 124 | 38 | 129 | 145 |
| span | 257 | 184 | 73 | 0.716 | 0.446 | 124 | 39 | 131 | 142 |

Bar 1 PASS (hits 184 = 184, true pairs lost 0). Bar 2 FAIL (new false 1 > new true 0). Bar 3 FAIL (placed 124, needs 125).
**FAIL.** Three members gain a model; none is placed (one lands in a wrong family, two in none).
As registered: `RUSTLE_MONO_FEED` and its driver test are removed; `mono_feed.py` is back in `bench/denovo_eval/`. Row closed.
With row 1665 (A119b +15 placed, KB3781 no gain) the feed helps only on the sample it was found on.
