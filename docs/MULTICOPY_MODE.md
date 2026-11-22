# Multi-copy mode: assembling the copies of known gene families

## Current state
(2026-11-21) The mode exists as two driver commands, uncommitted:

```
tools/rustle_pipeline.sh assemble  --seed-pool good --bam B --fasta G --out PREFIX
tools/rustle_pipeline.sh multicopy --seed-pool good --bam B --fasta G --out PREFIX
```

Output: `PREFIX.multicopy.gtf` (= `PREFIX.gtf` + rescued transcripts tagged `multicopy_rescue "1"`) and
`PREFIX.multicopy.tsv` (one line per rescued locus: gene, family, copies in the family, copies from the base assembly).
Since row 1871 every transcript also carries tie-adjusted counts and its tie group, and `PREFIX.multicopy.groups.tsv`
lists the groups (see "Counts and tie groups").
No other product changes. Opt-in; not part of `all`; no default was changed.

On the named multi-copy families (mainstays + Y-ampliconic) it models more annotated copies than StringTie, FLAIR and
IsoQuant on six testis samples, at StringTie's junction precision or better. It was confirmed on four samples it was not
designed on (row 1867). Whether a copy modeled from tied reads is expressed: the mode gives the tie group and its range of
expressed copies; the base-level verdicts are still a standalone script.

Goal this serves: `GOAL_multicopy_families_2026-11-19.md`. Pre-registrations with every number:
`PREREG_family_benchmark_2026-11-20.md`, `PREREG_family_by_age_2026-11-21.md`, `PREREG_family_by_family_2026-11-21.md`,
`PREREG_copy_evidence_2026-11-21.md`, `PREREG_multicopy_mode_2026-11-21.md`, `PREREG_copy_psv_2026-11-21.md`.

## What the mode does
1. **Near-tie secondaries seed loci** (`assemble --seed-pool good`): besides primaries, a secondary alignment within 2% of
   the molecule's best alignment score (AS) can start a locus. `assemble` on its own seeds from primaries only (row 1812).
2. **Rescue of copies no two reads agree on** (stage `multicopy`, `tools/multicopy_rescue.py`):
   - assemble again with the read floor at 1 (`PREFIX.multicopy.d/f1.gtf`);
   - take its spliced transcripts at places `PREFIX.gtf` leaves empty on that strand;
   - run the families stage on `PREFIX.gtf` + those candidates;
   - keep a candidate only when its locus is a copy of a family of >= 2 loci.

Switches: `RUSTLE_MULTICOPY_F1=FILE` reuses a floor-1 GTF of the same BAM and seeding pool. `RUSTLE_MULTICOPY_GATE`
(unset = 0) = 1 adds the pooling-guard bypass of row 1866: 2.4 x the candidates on OR6737, not judged for this mode.

## Counts and tie groups (row 1871)
Identical copies each get the same transcripts, and the assembler's `reads` counts the same molecules at each of them
(72 at each of three TSPY10 copies). The stage therefore runs `tools/eq_classes.py` on its GTF (ties on minimap2's `ms`,
exact band; `RUSTLE_EQ_WINDOW` for the wide band; `RUSTLE_MULTICOPY_EQ=0` skips the step) and adds per transcript:

| attribute | meaning |
|---|---|
| `reads_placed` | molecules whose best record is at this locus alone |
| `reads_max` | placed + every molecule tied here |
| `reads_em` | tied molecules shared out by EM; sums to the molecules of the group |
| `reads_tieadj` | `reads` x `reads_em` / `reads_max` of the locus (an approximation per transcript) |
| `tie_group`, `tie_group_loci` | the group of loci the reads do not separate, and its size |

`PREFIX.multicopy.groups.tsv`: per group its loci, placed and tied molecules, `min_expressed` (loci with >= 2 placed
molecules, at least 1) and `max_expressed` (all), and the members with their counts.

- A119b: 679 groups, 2,234 loci, between 1,798 and 2,234 expressed. OR6737: 213 groups, 540 loci, 374 to 540.
- TSPY on A119b: one group of 39 chrY loci (29 annotated copies, all 14 IDENTICAL ones), 165 molecules; `reads_em` sums to
  164.7, `reads_max` to 1,443; expressed copies 8 to 39.
- The EM gives 0.0 to a copy with no placed molecule when a tied partner has some (the three identical TSPY10 copies).
  That is one end of the range 0..`reads_max`: quote the group's range, not `reads_tieadj`, for such copies.
- The registered count check failed as written (0.55 of my reference count: it included 79 molecules that score better
  elsewhere); the conservation inside the group is post hoc.
- Cost: 47 min on A119b (the stage goes from 7 to 54 min), 2 min on OR6737.

**Is the EM share right? Tested where the origin is known** (rows 1873, 1874; `PREREG_tied_truth_2026-11-21.md`).
Measure: over tied reads, the mean share a rule gives the copy the read came from.

| test | EM (current) | equal split | |
|---|---|---|---|
| real A119b reads cut in pieces that become tied (297,936) | 0.507 | 0.396 | EM, in 184 of 244 groups |
| planted, every copy 60 reads | 0.279 | 0.337 | equal split |
| planted, one copy 600, others 20 | 0.313 | 0.320 | none |
| planted, half the copies silent | 0.305 | 0.338 | none |

No rule wins by the registered decision: EM stays. With half the copies silent, EM calls 10 of 63 silent copies
expressed and 6 of 70 expressed copies silent; the equal split 25 and 0. Neither identifies a silent identical copy.
The planted count fell outside [`reads_placed`, `reads_max`] for 4-8 of 133 copies (by a few reads): the range is a
guide, not a guarantee.

**Other near-tie cases** (row 1875, planted reads). Gene + processed pseudogene (150 genes, 417 pseudogenes): under `ms`
no read of a gene ties with a silent pseudogene and none is called expressed; under AS 3.1% tie and 1-5 are called.
Tandem pairs tie three times as often as dispersed pairs (0.29-0.31 of reads against 0.10); no alignment runs across
two copies.

**Known flaw, fix proposed, not applied** (rows 1875, 1877). `eq_classes.py` pools every tied record that lies on no
locus into one EM member `outside`, which takes tied reads from real copies (any copy the assembly lacks plays that
part). With one `outside` member per class: real pieces 0.507 -> 0.542, scenarios 0.279 -> 0.340, 0.313 -> 0.365,
0.305 -> 0.376; expressed copies called silent 15 -> 4, 15 -> 5, 6 -> 1. The tool also feeds stage `eqclasses`: the
change waits for the user.
- These ties are on `ms`; the OWN / TIE ONLY classes below are on AS. They need not agree copy by copy.

**Aligner settings `-p 0.9 -N 1000` change no member** (row 1885, `PREREG_minimap_p09_n1000_2026-11-22.md`): on planted
family reads and on the family reads of A119b and OR6737 the same copies are reached and modeled as with `-p 0.1 -N 50`,
with 16-80% of the secondaries. Whole OR6737 (row 1888): only 1.8% of the usual BAM's secondaries are near-tie; `-N 1000 -p 0.5`
keeps 99.7% of them in a 7.0 GB BAM (usual 11.7), `-p 0.8` 98.3%, `-p 0.9` 92.7%, `-p 0.95` 76.7%; named copies, family pairs and
intron chains are those of the usual BAM up to `-p 0.9`. Run directly (row 1894): `-N 50 -p 0.8` gives the usual named copies and
family pairs with 8% of the secondaries (BAM 5.8 GB). Two-pass alignment (`--write-junc`, `--pass1`) changes nothing measurable
and cannot be combined with `--eqx`. `-Y` already puts the sequence in secondary records. The mapping recipe is unchanged (the user's call; one sample).

**Slim mapping mode** (row 1902, opt-in): `sbatch --export=ALL,MAP_MODE=slim -J TAG genomes/map_reads.sbatch READS` maps with
`-p 0.8` instead of `-p 0.1` and writes `reads/TAG/TAG_slim.bam`: 8% of the secondaries, 98.5% of the near-tie ones, half the
BAM, the same named copies, family pairs and transcripts on OR6737. `--eqx` is kept: assembly, families and this mode do not
read =/X, but stage `flag` (PSV tallies of the missing-copy test) and the base-level copy test do. Do not use a slim BAM for
work that needs low-scoring secondaries: the join census of row 1881, and finding silent family members (row 1909: with
`-p 0.1` an expressed gene's reads reach 292 of 417 silent processed pseudogenes, with `-p 0.8` 20).

**The aligner's cap of 50 secondaries is not a limit here.** A119b: 367,525 molecules (1.7%) reach the cap; 8,896 of them
(2.4%) have all records within 2% of the best AS, none on a named-family copy; 81% have at most 5 near-tie records.

## Joins between tandem copies (row 1881, `PREREG_tandem_joins_2026-11-21.md`)
- Near-tied alignments do not join tandem copies more often than ordinary ones on A119b: 0.078% of near-tie records at
  2,319 tandem pairs cover both genes, 0.74% of untied primaries; ordinary neighbours 0.23%.
- On planted reads with no readthrough, secondaries do join loci (3.5% of molecules have such a secondary; primaries
  1-2 reads), and Rustle writes 1 and 0 joined transcripts with secondaries, 0 with primaries; StringTie 0.
- A119b, distinct transcripts joining a tandem pair: Rustle 13 (primaries), 15 (secondaries), StringTie 49.
- Known error of the mode: `DN_chr5_17456902_1`, one 113 kb single-exon model over the TAF11L array (946 reads), made
  by near-tie secondaries covering an array of single-exon copies end to end. Open.

## Why each half is there
| finding | row | number |
|---|---|---|
| With primaries only Rustle models fewer named copies than StringTie | 1859 | A119b 93 against 107 of 110 |
| The loss is at seeding, not at a filter | 1861 | copies flat across stages (91 -> 93; with secondaries 122 -> 123) |
| Secondaries recover the youngest families, by an amount that falls with age | 1858 | A119b pair recall: human-specific 0.588 -> 0.853, great ape 0.667 -> 0.733, Old World monkey 0.540 -> 0.573, older 0.223 -> 0.216 |
| Genome-wide, secondaries cost family-pair precision | 1858 | A119b 0.841 -> 0.650, OR6737 0.271 -> 0.090 |
| Some copies have own reads and no transcript: every read on a different intron chain | 1863, 1865 | 13 of 14 such copies; PLA2G4B 22 reads, 22 chains |
| The one-read tiers do not fix it; floor 1 everywhere costs junction precision | 1865 | 1-2 of 14 copies; A119b 0.49 -> 0.39 |
| Floor 1 only at empty places fixes it | 1865 (post hoc), 1867 | 11 of 14 copies; confirmed on unused samples |
| Genome-wide the rescued transcripts are rarely annotated, hence the confinement to families | 1867, 1869 | chain precision 0.03-0.10 against 0.18-0.37; confinement keeps 9-11% of them and 83 of 90 copies |

## Result on the named families
Copies with a spliced transcript / junction precision (share of introns on annotated exon boundaries).
Families: human `MAINSTAYS` + `YAGS` on A119b, `GGO_MAINSTAYS` on OR6737 (`bench/denovo_eval/famtruth.py`), the nineteen
of both lists found in each ape RefSeq GFF by symbol or description (`r1867_multicopy_judge/targets.py`).

| sample | Rustle, primaries only | Rustle, secondaries | mode (stage) | StringTie | FLAIR | IsoQuant |
|---|---|---|---|---|---|---|
| chimp testis | 94 / 0.81 | 129 / 0.81 | 155 / 0.76 | 121 / 0.73 | 89 / 0.70 | 81 / 0.80 |
| bonobo testis | 86 / 0.79 | 139 / 0.77 | 177 / 0.72 | 126 / 0.74 | 86 / 0.70 | 95 / 0.79 |
| Bornean orangutan testis | 165 / 0.77 | 176 / 0.74 | 182 / 0.73 | 173 / 0.66 | 165 / 0.50 | 158 / 0.81 |
| Sumatran orangutan testis | 135 / 0.79 | 174 / 0.75 | 187 / 0.73 | 167 / 0.70 | 152 / 0.55 | 129 / 0.80 |
| the four together | 480 / 0.785 | 618 / 0.760 | 701 / 0.732 | 587 / 0.699 | 492 / 0.561 | 463 / 0.801 |
| A119b, human testis (development) | 93 / 0.57 | 123 / 0.49 | 126 / 0.49 | 107 / 0.47 | 98 / 0.35 | 85 / 0.65 |
| OR6737, gorilla testis (development) | 25 / 0.61 | 27 / 0.58 | 34 / 0.66 | 28 / 0.59 | 23 / 0.55 | 23 / 0.56 |

Rows 1859, 1861, 1867, 1869. The largest gains over StringTie are TSPY, RBMY and NPIP. Behind StringTie: TBC1D3 on bonobo
and Sumatran orangutan, GOLGA8 on bonobo, by one or two copies.

## Is a modeled copy expressed
A secondary alignment places one read on every near-identical copy, so "modeled" is not "expressed". Two standalone
scripts on the named families (not in the stage):

**By alignment score** (`r1863_copy_evidence/ev2.py`, `tiegroups.py`; row 1863). A copy is OWN (>= 2 spliced molecules
score best there alone), TIE ONLY, or BORROWED (its reads score better elsewhere). A119b, mode's assembly: 85 OWN,
38 TIE ONLY, 9 BORROWED; expressed copies between 88 and 123. TSPY 6 to 32: the same ~200 molecules lie on every copy.
OR6737: all 36 OWN. Fifteen copies found with primaries only are TIE ONLY too: any tool reporting them rests on the
aligner's arbitrary primary.

**By bases that distinguish the copies** (`r1870_copy_psv/psv.py`; row 1870). Per copy: reads that match it at a base
where they mismatch every other family copy they align to (support), the reverse (contradiction), and copy pairs no
shared read tells apart (identical). Calibration on A119b: 81 of 85 OWN copies are SUPPORTED (95%); OR6737 29 of 36.
The 38 TIE ONLY copies: 7 SUPPORTED (BPY2B, CT47A5, CT47A6, HSFY1, HSFY2, two XAGE1A), 14 IDENTICAL (all TSPY: no
sequence-based method can separate them on these reads), 17 not needed (12 TSPY, CT47A1-4, XAGE1B). All copies:
88 SUPPORTED, 14 IDENTICAL, 19 not needed, 2 unresolved.

These are read-derived paralogue-specific positions, not a PSV catalogue: only positions under reads, substitutions only,
never checked against an alignment of the copies. "Not needed" is not proof of absence (TSPY2 has 12 own reads by score
and lands there).

## StringTie and secondary alignments
`stringtie -L` 3.0.3 reads secondary records: on chimp testis 62,418 transcripts with them, 62,494 without, 61,268
shared; the named families are identical either way (row 1870). Every StringTie assembly compared here was run on BAMs
that hold the secondaries. IsoQuant, FLAIR and bambu were not checked.

## Planted truth across tools (row 1903, `O2_FAMILY_BENCHMARK.md`)
On planted reads of the 133 named copies every tool reaches 0.87-0.99 of the expressed copies; the mode has no lead
(0.970 / 0.970 / 0.971 against StringTie 0.970 / 0.977 / 0.986). With 63 copies silent the mode gives a transcript to 42 of
them, StringTie to 25, every tool to at least 18. Rustle reproduces more planted isoforms (0.74-0.79 against 0.64-0.65)
at lower isoform precision. The table of real samples below counts copies modeled; it does not show them expressed.

## Wording
Supported: "On known multi-copy families, Rustle's multi-copy mode assembles transcripts at more annotated copies than
StringTie, FLAIR or IsoQuant run on the same alignments (701 against 587 for StringTie over four ape testis samples), at
equal or better junction precision. Within Rustle the gain comes from letting near-tie secondary alignments seed loci."

Not supported:
- "more expressed paralogs": for copies resting on tied reads, give the range or the base-level verdicts.
- "Rustle by default": with primaries only it is behind StringTie on these families (480 against 587).
- "StringTie discards secondaries": it reads them and does not turn them into copies.
- "better families than every tool genome-wide": row 1855 failed its bar on OR6737 recall (0.320 against StringTie 0.353),
  and with primaries only StringTie leads the human-specific pairs on A119b (0.658 against 0.588, row 1857).

## Limits
- A copy counts with any spliced transcript mapped to it, complete or not; annotated copies only.
- The rescued transcripts are single-read chains: 3-7% of their intron chains are in RefSeq.
- The four ape samples judged the unconfined candidate (row 1867) and then checked the confined stage (row 1869): the final
  code has no sample it has not seen.
- The expressed-or-not scripts ran on A119b and OR6737 only; reference copies are CHM13's, not the sample's; indels unscored.
- A119b and OR6737 were development samples for the family definition and for this mode.

## Open
- A PSV catalogue from an alignment of the copies (with indels), replacing the read-derived positions; `assign` already
  builds PSVs for read placement.
- The base-level verdicts (supported / identical / not needed) as an output of the stage (tie groups: done, row 1871).
- The 47 minutes of the counting step on a human sample.
- One unused testis sample for the final code.
- Commit: `tools/rustle_pipeline.sh` (stage `multicopy`; the file also carries another session's edits),
  `tools/multicopy_rescue.py`, `tools/test_selftests.py`, the six pre-registrations, register rows 1855-1877 (not 1872, 1876) and 1881, `PREREG_tandem_joins_2026-11-21.md`, `PREREG_tied_truth_2026-11-21.md`, this file.

## Where things are (Fritz workspace `/anvme/workspace/v138da13-rustle`)
| row | directory | what |
|---|---|---|
| 1855 | `r1847_family_benchmark` | one family definition on five tools' GTFs |
| 1857 | `r1857_family_by_age` | recall by duplication age |
| 1858 | `r1858_seed_pool_families` | primaries only against near-tie secondaries |
| 1859 | `r1859_family_by_family` | named families, five tools |
| 1861 | `r1861_multicopy_mode` | stage-by-stage diagnosis, `why.py` |
| 1863 | `r1863_copy_evidence` | own / tie / borrowed, tie groups |
| 1865 | `r1865_multicopy_oneread` | one-read tiers, the post hoc rescue |
| 1867 | `r1867_multicopy_judge` | four unused ape testis samples |
| 1869 | `r1869_multicopy_stage` | the driver stage, confined |
| 1870 | `r1870_copy_psv` | base-level copy evidence, StringTie with and without secondaries |
| 1871 | `r1871_multicopy_counts` | tie-adjusted counts, tie groups, the cap of 50 secondaries |
| 1873 | `r1873_tied_truth_clip` | EM against equal split on real reads cut in pieces (`parents.py`, `score.py`) |
| 1874 | `r1874_tied_truth_sim` | the same on planted reads (`sim.py`, `calls.py`) |
| 1875 | `r1875_tied_truth_cases` | gene + pseudogene, tandem and dispersed pairs (`plan.py`, `cases.py`) |
| 1877 | `r1877_outside_unpooled` | one no-locus member per class (`score_unpool.py`) |
| 1881 | `r1881_tandem_joins` | joins between tandem copies in reads, assemblies, planted reads |
| 1885 | `r1885_minimap_p09_n1000` | minimap2 -p 0.9 -N 1000 against the usual settings |
| 1894 | `r1894_minimap_two_pass` | two-pass minimap2, recipe -N 50 -p 0.8 (`arm.sbatch`, `junc.py`, `compare.sh`) |
| 1888 | `r1888_multimapper_limit` | OR6737 realigned whole at four values of -p (`arm.sbatch`, `stats.py`, `compare.sh`) |
