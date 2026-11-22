# Rows 1900 and 1901 (2026-11-22, user decisions): the refitted one-read ranking is the default; two single-exon pairing opt-ins

Both decided by the user on 2026-11-22 after rows 1897-1899 (asked in the session, both answered).

## Row 1900: `RUSTLE_ONEREAD_MODEL` unset = `r1899`
User: "Replace, keep old as switch", decided AGAINST row 1899's registered verdict (not adoptable: with each species held out
of its fit the refit is at least the old order on 12 of 17 datasets at StringTie's size, primates 11 of 13, mouse 1 of 4).

- `tools/data/oneread_r1899.json`: the 23 features and the design of `oneread_r1806.json`, fitted on twenty datasets at the
  defaults of 2026-11-22 (at most 20,000 one-read chains each; `r1899_oneread_refit/lodo.py`). `RUSTLE_ONEREAD_MODEL=r1806` =
  the ranking before; any other value stops the stage.
- It changes the ORDER of `PREFIX.oneread.ranked.gtf` and of tier 2 of `PREFIX.ranked.gtf`, hence the `list_sensitive` cut and,
  with partner libraries, the promoted set (best 20%). `PREFIX.gtf` and `PREFIX.final.gtf` are untouched.
- Check (workspace `r1900_oneread_model`, jobs 4495015-4495016, Sumatran testis and mouse heart): with `r1806` the ranked tier
  is byte-identical to row 1880's; with the default the list at StringTie's size holds 17,730 correct chains on Sumatran testis
  (+749 against StringTie; was -686) and 13,182 on heart (+605; was +544), the figures of row 1899's in-sample table.
- **What is and is not known.** In sample (every dataset in the fit) the new order is at least the old one on 15 of 17 at
  StringTie's size (mouse ES -59, HG002 -193) and on 18 of 20 in the best 10% of the tier. On a species the model has not
  seen, row 1899 says: better on 11 of 13 primate datasets, worse on three of four mouse libraries (-23 to -155 chains) when no
  mouse data is in the fit. The shipped model has mouse data in it; a new species has no such guarantee.
- With it the ranked list is ahead of StringTie on 20 of 20 datasets of row 1880 in sample (and under the held-out models).

## Row 1901: `RUSTLE_POLISH_MONO_SPARE=1`, `RUSTLE_POLISH_MONO_PAIR=1` (opt-in, both off by default)
User: "Both, opt-in" (row 1897: pairing single-exon transcripts with the multi-exon ones around them).

- `copy_assign --polish-mono-shadow-spare`: the shadow rule drops a single-exon transcript only when it overlaps an exon of a
  multi-exon transcript of ITS OWN strand (it can be a piece of that transcript). One wholly in an intron of its strand's
  gene, or covered only by an opposite-strand gene, is kept.
- `copy_assign --polish-mono-pair`: a multi-exon transcript gets `frag3_models "N"; frag3_reads "R";` for the N single-exon
  transcripts the polish dropped (steps `ism_mono`, `floor`, `shadow`; R reads together) that lie at least 90% inside its last
  exon with their 3' end within 50 bp of its own. A fragment under several isoforms that share the last exon is counted for
  each. An attribute; no transcript changes.
- Finding while building: a single-exon transcript inside a spliced transcript's span is first judged by the sub-chain step
  (`ism_mono`: dropped under 0.3 x the container's reads), so the shadow only ever sees the well-supported ones; the pair tag
  therefore counts fragments from all three single-exon drop steps, far more than row 1897's "shadowed" fragments.
- Check (workspace `r1901_mono_pair_build`, jobs 4495017-4495020; binary `~/ws/rustle_target_b`, `x.gtf`):
  | library | single-exon models, off -> spare | matching the truth | truth found | added models (matching) | multi-exon chains | transcripts tagged, fragment reads |
  |---|---|---|---|---|---|---|
  | mouse ES | 1,932 -> 2,408 | 266 -> 294 | 275 -> 304 | 476 (28, 5.9%) | unchanged | 9,852, 167,279 |
  | mouse heart | 87 -> 146 | 50 -> 59 | 50 -> 59 | 59 (9, 15.3%) | unchanged | 6,297, 57,834 |
  | OR6737 | 312 -> 437 | 76 -> 101 | 76 -> 101 | 125 (25, 20.0%) | unchanged | 10,368, 214,172 |
  | WTC11 r2 | 183 -> 227 | 113 -> 132 | 114 -> 133 | 44 (19, 43.2%) | unchanged | 10,087, 90,501 |
  The added models are row 1897's classes 4 + 5 (467 / 59 / 125 / 44 predicted). With both switches off the new binary's
  `x.gtf` is byte-identical to the shared binary's on all four. Unit test
  `mono_shadow_spare_keeps_the_intronic_and_the_opposite_strand_one_and_pair_tags_the_3p_fragment`; `copy_assign` unit tests
  120 of 120, driver tests 66 of 66.
- **Recall option, not a gain on both axes:** 11-33% more annotated single-exon genes found, single-exon precision lower on all
  four (0.138 -> 0.122, 0.575 -> 0.404, 0.244 -> 0.231, 0.617 -> 0.581).
- **The shared binaries are NOT rebuilt.** `~/rustle_target/release` (2026-11-20) does not know the two flags; other sessions'
  jobs were running on it when this was built. Until it is rebuilt the two switches need `--bin ~/ws/rustle_target_b/release`;
  with the shared binary and a switch set, `copy_assign` stops on the unknown option.

## Not done
- The public pages; the twelve default numbers under `r1899` (only the order of side files changes, no default count does).
- `families` / `assign` with the spared single-exon transcripts (they enter `PREFIX.gtf` when the switch is on).

## Tally of the ranked list against every other assembler under the shipped defaults (2026-11-22, `r1899_oneread_refit/current.py`)
The comparison of row 1880 (list cut at each rival's own output size, correct chains at least the rival's), recomputed with
the shipped ranking `r1899`. Every one of the twenty datasets is in that model's fit, so this is a description of the shipped
state, not a test; the test is row 1899.
| rival | AHEAD | BEHIND | note |
|---|---|---|---|
| StringTie | 20 of 20 | none | Sumatran testis +749 (was -686); the three simulated sets are inside tier 1 and unchanged |
| IsoQuant | 20 of 20 | none | all inside tier 1, unchanged from row 1880 |
| bambu (NDR 1) | 12 of 12 | none | all inside tier 1, unchanged |
| FLAIR | 14 of 15 | HG002, -6,171 | FLAIR is longer than the whole two-tier list there; `RUSTLE_RANKED_TAIL=1` (row 1891, opt-in) makes it +6,718 |
66 of 67 with the default two tiers, 67 of 67 with the opt-in third tier.
