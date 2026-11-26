# Multi-copy mode: assembling the copies of known gene families

## Current state
(2026-11-23, evening) The mode is the O2 stage of the default pipeline (user, 2026-11-22; committed 074557f7; the
labels' EM carries the alignment scores since e09fbf02; the sequence count of row 1926 is an opt-in):

```
tools/rustle_pipeline.sh all --bam B --fasta G --out PREFIX          # assemble (--seed-pool good), families, o2, flag
tools/rustle_pipeline.sh o2  --bam B --fasta G --out PREFIX          # the same step alone, after assemble + families
```

Stage `o2` runs stage `multicopy` (near-tie secondaries seed the family copies, floor-1 rescue inside families ->
`PREFIX.multicopy.gtf`, rescued transcripts tagged `multicopy_rescue "1"`, `PREFIX.multicopy.tsv`), then the provenance
labels (`o2_status` expressed / shared / homology_only from the EM over the tied reads; nothing dropped) ->
`PREFIX.o2.gtf`, `PREFIX.o2.em.tsv`, then the family report `PREFIX.o2.families.tsv`, `PREFIX.o2.copies.tsv`.
`RUSTLE_O2_SEQUENCES=1` (opt-in, register 1926) adds `PREFIX.o2.sequences.tsv` and the report column `tied_sequences`
(IsoCon's test ported to the genome alignments of the tied reads; sequences, not copies: alleles count).
`RUSTLE_O2_EXTEND=0` turns the extension off; `RUSTLE_O2_DROP=1` drops the homology_only loci; `RUSTLE_O2_EM_WEIGHT=0`
puts the labels back on the counts-only EM (default since 2026-11-23, register 1921: scores inside the EM); `PREFIX.gtf` is never
changed. `multicopy` can still be called alone; with `RUSTLE_MULTICOPY_EQ=1` it also writes tie-adjusted counts and
`PREFIX.multicopy.groups.tsv` (see "Counts and tie groups"; inside `o2` that step is off, the labels replace it).
Validation and every planted number: `O2_FAMILY_BENCHMARK.md`.

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
CORRECTED 2026-11-25 (`docs/MULTIMAPPERS_ACROSS_ASSEMBLERS_2026-11-25.md`): in long-read mode StringTie DROPS secondary and
supplementary alignments before its splice graph (`rlink.cpp:713`); they only move bundle boundaries, which is why on chimp
testis the output differs slightly (62,418 transcripts with them, 62,494 without, 61,268 shared; named families identical;
row 1870), and on planted reads not at all. FLAIR, IsoQuant and bambu were checked the same day, in source and by running:
none uses a secondary alignment; all four give the same assembly with and without them.

## Planted truth across tools (`O2_FAMILY_BENCHMARK.md`; rows 1917, 1918, reads with ragged ends)
- Isoforms: Rustle reproduces the most of the genome-based tools at every coverage (0.82 against 0.66-0.82 at 60 reads
  per copy, 0.67 against 0.31 at 3; gorilla's own families 0.96-0.98 against 0.85-0.94).
- Copies: near-tie secondaries add expressed copies only at about 3 reads per copy (the mode 125 of 133, StringTie 95);
  from 5 reads up every tool reaches them. On copies with ordinary introns Rustle finds 267 of 267; it misses
  single-exon models and lift-over models with non-canonical introns.
- False copies: taking every near tie gives 24 false units with half the copies silent, the other tools 9-16, the mode
  with the EM test 3. The default labels instead of dropping: 86-92 units available against 75-76, 0-2 expressed units
  mislabelled.
- The earlier planted numbers (rows 1903-1911) were made with reads whose ends coincide, which Rustle counts as
  duplicates; they are superseded for Rustle against other tools.
The table of real samples above counts copies modeled; it does not show them expressed.

## Wording
Supported, real samples: "On known multi-copy families Rustle's multi-copy mode has transcripts at more annotated copies
than StringTie, FLAIR or IsoQuant run on the same alignments (701 against 587 for StringTie over four ape testis
samples), at equal or better junction precision. Within Rustle the gain comes from near-tie secondaries seeding loci."
Supported, planted reads: the three statements of `O2_FAMILY_BENCHMARK.md` (most isoforms; fewest false copies with the
provenance test; most members available with labels), and copies gained at about three reads per copy.

Not supported:
- "more expressed paralogs" / "the extra copies are expressed": give the label (`o2_status`), the tie-group range or the
  base-level verdicts.
- "Rustle finds more copies" in general: on planted reads only at very low coverage.
- "Rustle by default with primaries only": on real samples it is behind StringTie on these families (480 against 587).
- "StringTie uses secondaries" (this file said so until 2026-11-25): in long-read mode it drops them before building; they only move bundle boundaries. Say "none of the four comparison assemblers uses a secondary alignment as evidence".
- "better families than every tool genome-wide": row 1855 failed its bar on OR6737 recall (0.320 against StringTie 0.353),
  and with primaries only StringTie leads the human-specific pairs on A119b (0.658 against 0.588, row 1857).

## Limits
- A copy counts with any spliced transcript mapped to it, complete or not; annotated copies only.
- The rescued transcripts are single-read chains: 3-7% of their intron chains are in RefSeq.
- The four ape samples judged the unconfined candidate (row 1867) and then checked the confined stage (row 1869): the final
  code has no sample it has not seen.
- The expressed-or-not scripts ran on A119b and OR6737 only; reference copies are CHM13's, not the sample's; indels unscored.
- A119b and OR6737 were development samples for the family definition and for this mode.

## Tried and not adopted
- `RUSTLE_MULTICOPY_MONO=1` (row 1919, opt-in, off): single-exon models of a floor-free run kept only as copies of a
  family of >= 2 loci (`multicopy_rescue.py candidates --mono`, driver env `RUSTLE_POLISH_MONO_QUANTILE`). Planted
  gorilla: expressed single-exon copies 101 -> 119 of 147 (FLAIR 146). Real OR6737: the 728 rescued models are on an
  annotated exon 13% of the time against 56% for the default's single-exon models. Family membership does not pick out
  real single-exon genes. `PREREG_mono_family_rescue_2026-11-23.md`.
- Pairwise decisive-read verdicts (`tools/copy_verdicts.py`, row 1904): 42 -> 40 false copies; the joint EM test is used.
- Restoring the single-exon models a primaries-only assembly passes (row 1920, `multicopy_rescue.py mono_primary`):
  closes the planted gap (119 of 147 against 98), but real libraries have almost nothing to restore (OR6737 one model;
  A119b 224, 6 annotated). Not built. The single-exon gap to FLAIR on planted reads is a benchmark figure.
- Overruling minimap2's choice from the read itself (rows 1922, 1923): its wrong calls are exact ties identical in
  every tag; error-aware re-scoring moves nothing on KB3781 and harms HG002. Closed.
- A variant counter of my own (`sample_psv.py`, row 1924) and running IsoCon on the tied groups (`isocon_groups.sh`,
  row 1925): superseded by the ported test of row 1926.

## Open
- The sequence count (`RUSTLE_O2_SEQUENCES=1`) on real reads: HiFi errors, no truth; alleles against copies (row 1927);
  a sample that lacks or gains a copy.
- Isoform precision with near-tie secondaries; strongly supported non-canonical junctions; the EM test's one to two
  units and its brittleness at two-read loci on real data.
- One unused real testis sample for the final code; the extension's run time on human (A119b three pools: 2.4 h).
- Known error: `DN_chr5_17456902_1`, a 113 kb single-exon model over the TAF11L array made by near-tie secondaries (row 1881).

## Alignment scores inside the EM (row 1921, 2026-11-23): confirmed, small, the default
The default EM of the labels is counts-only: the `ms` scores decide which loci a read is exactly tied between, and a read
with a unique best is placed outright. `em_filter.py --weighted 1` makes every alignment within 10 `ms` a candidate with
weight exp(-gap) and lets abundance x weight share the read. On planted human families it removes one silent unit wrongly
labelled `expressed` per run on average (five seeds 23 -> 21 mislabelled silent units, three unseen seeds 12 -> 9) and
changes no expressed unit, the 10-fold unequal-expression scenario included; gorilla unchanged. A softer weight (T = 5)
demotes real expressed units and is rejected. What it does not touch: silent units labelled `shared`, which are exact
ties over the whole read. DEFAULT of stage `o2` since 2026-11-23 (user); `RUSTLE_O2_EM_WEIGHT=0` = the counts-only EM
(benchmark runner: `O2_EM_WEIGHT`). Planted results quoted before this date used the counts-only EM.
On real gorilla OR6737 (no truth) it changes 18 of 18,868 `expressed` loci to `homology_only`, all with 1-4 reads; a
locus with exactly two reads and a similar relative one substitution away is where it is brittle.
`PREREG_weighted_em_2026-11-23.md`.

## Sequences among copies the reference cannot separate (row 1926, opt-in)
`RUSTLE_O2_SEQUENCES=1` runs `bench/o2_families/recurrence_test.py` on the stage's assembly, BAM and score table:
loci joined by exactly tied reads form a group; its tied reads are read at loci of the group, split by intron chain;
differences from the reference that recur are candidate sequences; IsoCon's statistical test (reads carrying every
separating variant, per-read error probabilities, Raghavan's bound, iterative merge at 0.01) decides which are real.
-> `PREFIX.o2.sequences.tsv` and the column `tied_sequences` of the report. On planted samples with per-copy variants
the count is exact in 0.80 of (unit, isoform) cells on unseen seeds, against 0.57 for IsoCon on the same groups and
0.50 for the best constant; with no variants it never over-counts. It counts sequences: two alleles of one copy are two,
and on a planted heterozygous sample the count is exact on haplotypes in 0.53 of cells and over-counts copies in 140 of
199 (row 1927); it never says which copy. Real OR6737: 392 sequences at 194 loci, 116 of 379 families.
Real samples with both haplotypes (rows 1928-1931): the non-reference sequences are in the individual's own DNA for
476 of 648 on HG002 (0.735; 0.15 nowhere; hypervariable loci such as the immunoglobulin lambda locus are not tested and
are listed in `PREFIX.o2.sequences.skipped.tsv`) and 52 of 107 on KB3781 (0.49; 0.44 nowhere, unexplained and specific to
that sample). Read the column as recurrent sequences with that confirmation rate, not as copies. On a deep human
library the test adds about an hour.
`PREREG_isocon_port_2026-11-23.md`.

## One floor-1 assembly for two stages (row 1937, default since 2026-11-24)
The multi-copy stage takes its candidates from a one-read-floor assembly made WITH the pooling-guard bypass
(`RUSTLE_MULTICOPY_GATE` unset = 1; 0 = before), which is the run the one-read tiers make. Whichever of the two stages
runs second reuses the other's GTF when the recorded settings are equal (`RUSTLE_F1_SHARE=0` = no sharing). Judged on
planted reads (units and labels unchanged; the extra rescued silent units are all `homology_only`) and on OR6737 /
KB3781 (annotated copies 34 -> 37 and 27 -> 28 at equal junction precision; rescued loci 445 -> 1,688 and 161 -> 514).
Saves one whole-genome assembly per sample (about 1.9 h of 5 on HG002). `PREREG_shared_floor1_2026-11-24.md`.

## One pass for the tie-based steps (`tools/tie_table.py`, 2026-11-24)
The labels and the sequence test both need, for every read with a second alignment within 10 `ms` of its best, the loci
its alignments fall on. That pass over the BAM is now made once, by whichever step runs first, and written as
`PREFIX.o2.ties.tsv.gz` (read, locus, `ms` gap to the read's best, best `ms`) with `PREFIX.o2.ties.tsv.gz.placed.tsv`
(per locus, the reads placed there with no rival); both steps take `--ties FILE`, and a table made from other inputs
is ignored and rewritten. The pass stops when the score table's header records another BAM size than the BAM given
(`bam_bytes`; 2026-11-25, after a shared table was overwritten for a reduced BAM). Outputs unchanged (workspace `tie_table_2026-11-24`): labels identical without the table,
while writing it and when reading it, and identical to the stage's earlier tables on OR6737 and KB3781; sequences
identical. Time: labels 142 s (OR6737) / 269 s (KB3781) when they make the pass, 8 s / 7 s when they read the table; the
sequence test 6 s / 11 s reading it, against 157 s / 282 s when it made its own pass.

## Which transcripts a primary alignment carries (default tag, 2026-11-25, register 1952)
Stage `o2` tags every transcript line of `PREFIX.o2.gtf` with `support "primary"` (at least one primary alignment has its
intron chain; for a single-exon transcript, lies inside it) or `support "tied_only"` (none: the transcript rests on
near-tied secondary alignments alone). Nothing is dropped. `RUSTLE_O2_SUPPORT=0` switches it off. It costs the provenance
pass over the BAM and writes that stage's tables (`PREFIX.provenance.*`); tool `bench/o2_families/support_tier.py`.
Why (row 1952, three fresh planted seeds): at 60 reads per copy the `tied_only` transcripts are 0% exact (0 of 364) and
the `primary` ones have isoform precision 0.80-0.82 at the same recall (whole mode 0.65; FLAIR 0.73-0.77); at 3 reads the
`primary` ones keep 0.59-0.61 isoform recall at 0.81-0.84 precision and the `tied_only` ones add 13-15% of the exact
isoforms at 0.42-0.48. On gorilla OR6737 2,455 of 73,173 transcripts are `tied_only`. Junction precision at annotated
copies cannot tell the tiers apart (a paralog's isoform laid onto a copy has that copy's junctions): do not use it to
argue that `tied_only` transcripts are expressed.

**DEFINITION CHANGED 2026-11-26: contested primaries no longer count; see "The layers without contested primaries" at the end of this file. The numbers in this section are for the definition before.**

**Three layers and their products (user 2026-11-25).** The same step also writes `primary_reads "N"` (primary alignments
that carry the chain) and `layer`: `primary` (N >= 2: stands without tied reads), `tied_support` (N = 1: tied reads lifted
it over the read floor), `tied_only` (N = 0). Products:
- `PREFIX.o2.primary.gtf`: the primaries-only product (layer `primary`). On planted reads it reproduces a primaries-only
  assembly (82 / 69 / 72 exact isoforms at 3 reads per copy against 83 / 70 / 72) without a second run.
- `PREFIX.o2.tied.gtf`: what the near-tied multimappers add (layers `tied_support` and `tied_only`, tagged). At 3 reads per
  copy `tied_support` doubles the exact isoforms (+60 to +74 at precision 0.72) and `tied_only` adds 22-24 at 0.45; at 60
  reads per copy the two hold nothing correct.
- `PREFIX.o2.layers.tsv`: per layer the transcripts, loci and labels, the loci each tied layer adds, and what happens to
  the near-tied reads (in which layer's transcript their primary sits, or why no transcript carries them). Gorilla OR6737:
  of 16,275 near-tied reads 64% sit in a `primary`-layer transcript, 4% in a tied layer, and 25% are carried by no
  transcript (unspliced with no model there, under the read floor, non-canonical junction). Since 2026-11-25 the table
  also looks each such read up at its TIED places (the tie table of the stage): a read whose tied alignment sits at
  another assembled locus is listed apart, and the last two rows give the reads carried at NO tied locus and their share,
  which is the unexplained part (OR6737: 1,523 of 16,275, 9.4%; the primary of a tied read is the aligner's arbitrary
  choice, so counting by the primary alone overstated it at 25%). Described the same day
  (`docs/PREREG_support_tier_2026-11-25.md`, last section): 55% of those have their TIED alignment on an assembled transcript
  elsewhere, leaving 11% of near-tied reads carried at no tied place; the unspliced ones without a model are mostly single
  reads tied between intronic or unannotated places, plus one clear class, intronless copies of ribosomal-protein genes
  tied with their spliced parent (207 reads).

## The EM's share of each tied read (opt-in, 2026-11-25)
`RUSTLE_O2_SHARES=1` (or `em_filter.py --shares FILE.tsv.gz`, with the weighted EM) writes `PREFIX.o2.shares.tsv.gz`: one
row per tied read and candidate locus with the `ms` gap and the share of the read the EM gives that locus (shares of a
read sum to 1). OR6737: 14,364 tied reads, 35,088 rows; labels and GTF unchanged by the option. A share is how the EM
divides a read between loci it cannot tell apart, not where the read came from.

## Where each read went (opt-in stage `provenance`, 2026-11-23)
`tools/rustle_pipeline.sh provenance --bam B --fasta G --out PREFIX` (after `assemble`; uses `PREFIX.o2.gtf` when it is
there) writes, with `tools/read_provenance.py`:
- `PREFIX.provenance.reads.tsv.gz`: one row per primary alignment: the transcript the read supports and how (exact;
  exact after a junction replacement, old > new; a shorter read of one transcript; a piece of several; unspliced; none
  and why), whether the short-exon repair rewrote it, whether it is a near-tied multimapper, its gene's `o2_status`, and
  for a tied read the sequence it supports.
- `PREFIX.provenance.transcripts.tsv`: the reads behind every transcript by kind, beside the assembler's own count.
- `PREFIX.provenance.family_reads.tsv.gz`: the rows of multi-copy families and near-tied reads.
It is RECONSTRUCTED from the BAM and the assembly: the assembler keeps no per-read record. Exact: the alignment, the
match to a transcript, the rewrite, the tie flag, the label, the sequence. Inferred from the outcome: junction
replacements, and the reason a chain gave no transcript (non-canonical junction; under the read floor; "removed by a
later rule", not resolved further). A tied read is never given a copy of origin.
Check against the assembler's read counts (assemblies from primary alignments; workspace `provenance_2026-11-23`):
planted human families 260 of 267 transcripts equal, all within 10%; gorilla OR6737 55,061 of 68,847 equal (0.80),
67,267 within 10% or one read (0.98), 219 with no read found. OR6737, 4.4 M reads: 0.69 exact, 0.015 exact after a
replacement, 0.035 a shorter read of one transcript, 0.062 a piece of several, 0.055 unspliced, 0.145 none. About one
minute on 72 threads.

## The polishing rules in Python (`tools/polish_rules.py`, 2026-11-24)
A line-by-line port of the assembler's polishing routine (`polish_drop_set_with`, `src/copy_assign_bin.rs`), written to
be read: one short block per rule, in the assembler's order (sub-chain collapse `ism`, single exon inside a spliced span
`ism_mono`, single-exon read floor `floor`, `shadow`, isoform fraction `frac`, retained intron `ret`), and a table naming
the rule that removed each transcript. `RUSTLE_ASSEMBLY_POLISH=none` makes `assemble` write the transcripts before the
rules (it also switches off junction snapping, the 5' rescue and junction replacement, which the assembler refuses
without its polishing: the file is the transcripts BEFORE snapping).
Parity, those steps off on both sides (workspace `polish_parity_2026-11-24`): planted sets 267 = 267 and 301 = 301
transcripts kept; gorilla OR6737 every rule's count equal to the assembler's (sub-chain 58,333; floor + shadow 28,772;
fraction 28,907; retained 2,392), the kept sets equal but for 92 single-exon transcripts removed by the internal-priming
test, a separate rule that needs the reads. One fact the port brought out: the assembler applies its per-gene rules
BEFORE it regroups genes (ids ending `.rgN`, `.fusN`); the script judges each transcript in the gene it had then.
Against a normal run (snapping, spliced share, rescue, replacement on): 70,124 of its 70,181 transcripts are kept by the
Python rules too, which keep 1,566 more.
With `read_provenance.py --unpolished U.gtf --reasons R.tsv` the reads whose chain "passed the floor but has no
transcript" get the rule by name: on OR6737 108,044 of 121,013 such reads (frac 96,115; ret 10,928; ism 1,001); 12,969
stay "removed by a later rule". These reasons come from the rules applied before snapping, so they are close, not
exact.
**Exact reasons from the assembler (2026-11-24, user: rebuild the shared binary).** `copy_assign --polish-reasons FILE`
writes one row per transcript the polishing removed (name, contig, strand, exons as they stood when the rules ran,
reads, rule: `spliced_share`, `ism`, `ism_mono`, `floor`, `shadow`, `frac`, `ret`, `mono_primed`, `collapsed_or_other`).
The assemble stage asks for it as `PREFIX.polish_reasons.tsv` (`RUSTLE_POLISH_REASONS=0` = no table; a binary without
the option is detected and left alone) and stage `provenance` reads it. No other output changes: `PREFIX.gtf` and
`PREFIX.final.gtf` byte-identical on a planted set and on OR6737. OR6737: 118,353 removed transcripts, their counts
per rule equal to the assembler's log (sub-chain 58,440; floor + shadow 28,642; fraction 26,782; retained 2,375; plus
spliced share 1,134, internal priming 89, other 891). Provenance with it: of the reads whose chain passed the floor
and has no transcript, fraction 93,534, retained 10,868, spliced share 4,568, sub-chain 943; 11,103 remain unnamed
(their chain was never a transcript: merged by snapping, or removed before the rules).

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
| 1903 | `r1903_o2_family_benchmark` | O2 benchmark, seven tools on planted reads (old reads); driver checks |
| 1904-1905 | `r1904_copy_verdicts`, `r1905_em_filter` | pairwise verdicts; the EM provenance test |
| 1906-1908, 1911 | `r1906_yag_cases`, `r1907_o2_replication`, `r1908_depth_sweep`, `r1911_label_validation` | units, IsoCon, replication, depth, labels (old reads) |
| 1909, 1910, 1913-1915 | `r1909_silent_members`, `r1910_member_discovery`, `r1913_members_real`, `r1914_read_relatedness`, `r1915_members_blocks` | silent relatives, the member layer |
| 1916-1918 | `r1916_isoform_loss`, `r1917_ragged_reads`, `r1918_rerun_ragged` | the read-end artefact and the corrected reruns |
| - | `annot_related` | relatedness annotation of CHM13, NCBI gene_group |

## Which rules use multi-mapping reads: the state on 2026-11-25
Two bands are in use. They are different thresholds on different scores and must not be mixed in wording.

| use | band | status |
|---|---|---|
| Seeding the assembly (stage `assemble` in `all`, pool `good`) | secondary AS >= 0.98 x the read's best AS | default |
| Tie table, labels by the EM, `support` / `layer` tags (stage `o2`) | secondary `ms` within 10 of the primary's | default |
| Distinct sequences among tied copies (`RUSTLE_O2_SEQUENCES=1`) | `ms` within 10 | opt-in |
| Relatives list (stage `members`, with `RUSTLE_MEMBERS_REPEATS`) | >= 0.9 of the best score | opt-in; a list of relatives, never members |

Every other rule tried on the secondaries was not adopted:

| rule | row | why not |
|---|---|---|
| Pool `ms` within 10 for seeding | 1942, 1942b | more precise, fewer copies; did not confirm on unused samples |
| Pool of exact `ms` ties only | 1951, 1958 | helps only at high depth; loses true transcripts at 3 reads |
| Pool with the intron allowance | 1953 | harmless on real data, adds nothing correct |
| Pools at 0.95 and 0.90 of the best AS | 1959 | adds transcripts with no read of their own |
| Anchored pools at 0.9 (another session) | 1961 | planted: safe with the intron condition, no true transcript gained; real arms unfinished when this was written |
| Per-read re-placement, error-aware re-scoring | 1922, 1923 | nothing inside one read decides a tie |
| A MAPQ filter | `coin_toss_2026-11-25` | removes every copy that has an identical twin |
| Members or silent copies from secondaries | 1915, 1955-1957, 1960 | secondaries mark homology, not membership or expression |
| Ties from a graph instead of secondaries | 1943-1945 | the graph gives the places, not the band |
| Label `processed_copy_of` (intronless against spliced, >= 0.9) | 1962 | 13.8% of places on a processed copy on A119b, reach 5.4% |
| Seeding in two steps: uniquely placed primaries first, tied reads conditioned on them | 1963 | the unique core has no made-up copy but a third of the isoforms; the conditioned second step loses 3-5 true isoforms and removes almost no made-up copy |

Census behind row 1962 (OR6737, workspace `multimapper_rules_2026-11-25/census.or.txt`, a description on one sample):
105,396 of 4.40 million reads have a secondary at 0.9 or more of the primary's `ms` (7,076 ties, 9,194 within 10,
89,126 below that). Three more patterns were counted and none is a rule: a join of two RefSeq genes has a tied or
near-tied alternative that stays in one gene for 198 of 14,652 reads (1.4%; so read-throughs are not a mapping
ambiguity), an alternative chain inside the same locus exists for about 440 reads, and of 236 tie-partner places on no
RefSeq exon 10 are in the CAT annotation.

### Two products side by side: unique primaries, and the near-tie pool (2026-11-25, description)
Workspace `two_products_2026-11-25` (job 4509336, `run.sbatch` is the recipe). Product 1 = primary alignments with
mapping quality 1 or more, no secondaries (`samtools view -b -F 0x100 -q 1`, then `assemble --seed-pool primary`; the
assembler has no mapping-quality option). Product 2 = today's assembly with the 0.98 pool. Mapping quality 0 removes
the exact ties only; near ties that are not exact stay in product 1.

| planted sets, seeds 1967-1969 | unique primaries | all primaries | near-tie pool 0.98 | multi-copy mode |
|---|---|---|---|---|
| all expressed, 60 reads: exact isoforms of 239 (precision) | 108 (0.86-0.89) | 195-197 (0.81-0.83) | 197 (0.65-0.66) | 197 (0.65-0.66) |
| half silent, 60 reads: exact isoforms (precision) | 56-57 (0.85-0.89) | 102-106 (0.62-0.65) | 103-105 (0.42-0.49) | 103-105 (0.42-0.49) |
| half silent: silent copies with a transcript, of 63 | 0 / 0 / 0 | 23 / 25 / 24 | 44 / 47 / 39 | 45 / 48 / 40 |
| all expressed, 3 reads: exact isoforms (precision) | 44-47 (0.90-0.96) | 70-83 (0.92-0.95) | 137-148 (0.77-0.82) | 162-168 (0.73-0.76) |
| all expressed, 3 reads: expressed copies found, of 133 | 47-50 | 73-82 | 110-113 | 128 |

Real OR6737: 11,291 of 4,409,946 primaries have mapping quality 0 (0.26%; 42% in the planted sets). Transcripts 68,640
(unique primaries), 68,847 (all primaries), 70,181 (pool 0.98), 71,443 (mode). Annotated copies of the four named
families modelled / junction precision: 25 / 0.61, 25 / 0.61, 27 / 0.58, 34 / 0.66.
Reading: on planted families product 1 has no made-up copy and about half of the isoforms; on a real sample it is the
all-primaries assembly within 0.3%, because exact ties are rare there. The near ties earn their place at low depth
(3 reads per isoform: 2 to 3 times the exact isoforms of either primaries product).

### The depth effect on real reads (row 1964, 2026-11-26)
OR6737 and KB3781 downsampled to 5 / 10 / 25%, truth = the full-depth primaries-only assembly
(`PREREG_real_downsample_2026-11-26.md`). The 0.98 pool adds 0.3-0.5% confirmed chains at every depth, and 0.52 / 0.51
of what it adds is confirmed at 5%, 0.29 / 0.30 at 25% (primaries themselves: 0.90-0.92). The mode adds 2.8% / 1.2% at
5% and 1.0% / 0.5% at 25%, a third to a sixth of it confirmed. The planted "two to three times the isoforms at 3 reads"
is a statement about multi-copy families, not about a whole sample.

### The layers without contested primaries (default since 2026-11-26, user)
The primary alignment of a read with a rival within 10 `ms` was placed by the aligner's tie-break. It no longer counts
as primary support: `primary_reads "N"` is the number of UNCONTESTED primaries that carry the chain, the new attribute
`contested_primaries "M"` holds the ones set aside, and the layers follow N as before (`primary` N >= 2, `tied_support`
N = 1, `tied_only` N = 0). A transcript that only contested primaries carry is now `tied_only` and is written to
`PREFIX.o2.tied.gtf`. The assembly itself is unchanged: same transcripts, other tags. Tool
`bench/o2_families/support_tier.py` (with `--reads`, as the stage calls it); check in workspace
`contested_primaries_2026-11-26` (job 4509358; self-tests and `cargo test` pass).

| | layer `primary` before | layer `primary` now |
|---|---|---|
| OR6737: transcripts / loci | 68,657 / 19,097 | 68,121 / 18,860 |
| OR6737: tied_support / tied_only transcripts | 1,133 / 1,653 | 874 / 2,448 |
| OR6737: annotated copies of the four named families (junction precision) | 25 (0.58) | 21 (0.59) |
| planted, half silent, 60 reads: made-up units | 8 / 11 / 6 | 0 / 0 / 0 |
| planted, half silent, 60 reads: exact isoforms (precision) | 102-105 (0.64-0.65) | 32-34 (0.81-0.92) |
| planted, all expressed, 60 reads: exact isoforms (precision) | 196-197 (0.81-0.83) | 65 (0.83-0.88) |
| planted, all expressed, 3 reads: exact isoforms (precision) | 69-82 (0.93-0.95) | 25-31 (0.86-0.97) |

What this changes in the wording: on a real sample 0.37% of the primaries are contested and the `primary` layer loses
536 transcripts and 237 loci (0.8% and 1.2%). In multi-copy families the layer is now the part no tie-break touched:
no made-up copy on planted reads, and a third of the isoforms. The tied layers are no longer "almost all wrong": at 60
reads `tied_only` held 0 exact isoforms of 58 before and holds 130 of 224 now (seed 1967), because true isoforms of
copies with an identical twin sit there. Row 1952's sentences about the tied layers describe the definition before.
