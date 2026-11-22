# O2: known gene families across tools, on planted and real reads

## Current state
(2026-11-22) O2 was replaced on the user's word (register 1903). It is no longer per-read copy assignment
(`assign`, `copy_assign`; record in `LEDGER_O2.md`, kept as legacy). O2 now asks, for every assembler on the same
alignments: **which copies and isoforms of known multi-copy families does it report, and are they the expressed ones?**
Two arms:
- **Planted reads** (`bench/o2_families/`): reads written from annotated transcripts of each family copy, so the copy and
  isoform of origin of every read is known. This arm has truth for silent copies.
- **Real samples** (`bench/denovo_eval/toolcompare.py`): annotated copies with a transcript, per tool. More realistic,
  no truth for which copies are expressed.
In the driver (`tools/rustle_pipeline.sh`), `all` = assemble, families, **o2**, flag: the stage `o2` took the place of
`assign`. Since 2026-11-22 (user: "lets readd it by default") `o2` first runs the multi-copy extension (near-tied
multimappers seed the family copies) and the provenance labels (`o2_status` expressed / shared / homology_only; nothing
dropped) -> `PREFIX.o2.gtf`, `PREFIX.o2.em.tsv`, then the family report on that. `RUSTLE_O2_EXTEND=0` turns the extension
off; if it fails the report is made on `PREFIX.gtf` and the stage succeeds (checked: default, off, forced failure;
27 dispatch tests pass). `PREFIX.gtf` is never changed. `assign` is in NO plan of `all` any more (not with `--candidates` or `--legacy-catalog`, and there is no switch); the
stage is kept and runs when called by name (`rustle_pipeline.sh assign ...`). Dispatch tests updated (`tests/pipeline_driver_dispatch.rs`). Uncommitted.
What a run of `all` no longer writes by default: `PREFIX.assign.*` (assignments, quant, families, PSV evidence).
`tools/check_run_tables.py` and `tools/rustle_samples.sh` (which calls `assign` itself) need the legacy stage.

## How to run
```
tools/rustle_pipeline.sh all     --bam B --fasta G --out PREFIX     # ... families, o2, flag
tools/rustle_pipeline.sh o2      --bam B --fasta G --out PREFIX     # PREFIX.o2.families.tsv, PREFIX.o2.copies.tsv
tools/rustle_pipeline.sh o2bench --gff ANNOT.gff3.gz --bam B --fasta G --out PREFIX   # planted benchmark, PREFIX.o2bench.tsv
```
`o2` (`tools/o2_family_report.py`, seconds, no BAM pass): per family of `families` its copies, copies with a spliced
transcript, transcripts, isoforms and reads; when `PREFIX.eq.loci.tsv` (stage `eqclasses`) or
`PREFIX.multicopy.eq.loci.tsv` exists, also reads placed on one copy, tie groups and the range of expressed copies.
The `reads` column is the assembler's count per transcript, summed: a read tied between identical copies is counted at
each of them (row 1871; 7,980 planted reads give 48,083 in one family), so read it per copy, never as a family total; the
tie-group columns are the ones that count each read once.
Verified 2026-11-22: 27 dispatch tests pass; `all` on the planted BAM runs assemble, families, o2 (17 families, 148
copies) (`r1903_o2_family_benchmark/e2e_new`; the legacy switch tried there was removed later the same day).
`o2bench` is opt-in and human-only for now. The pieces, by hand:
```
python3 bench/o2_families/plant.py ANNOT.gff3.gz GENOME.fa PLANT_DIR            # loci.gtf, plan.tsv, {E,D,S}.fq + truth
bash bench/o2_families/run_tools.sh PLANT_DIR E GENOME.fa OUTDIR [THREADS]      # align, every tool, score.tsv, reads.txt
python3 bench/o2_families/score.py PLANT_DIR E --families LABEL=GTF ...         # any assemblies against the truth
python3 bench/o2_families/tied_share.py LOCI.gtf BAM AS_TABLE TRUTH.tsv OUT     # reads: placed / tied / sharing rules
```
`plant.py`: members of `MAINSTAYS` + `YAGS` (`famtruth.py`), up to `--isoforms` (2) annotated isoforms per copy, reads
2 : 1; scenarios E (every copy `--reads`), D (one copy 10 x, the others a third), S (half the copies of a family silent);
`--seed`; errors and 5' truncation. `run_tools.sh` runs FLAIR, IsoQuant, bambu when `O2_FLAIR_PATH`,
`O2_ISOQUANT_MANIFEST`, `O2_BAMBU_MANIFEST` are set (Fritz values in `r1903_o2_family_benchmark/run.sbatch`); `MAP_P=0.8`
= the slim mapping mode. `plant.py --selftest`, `score.py --selftest`.

Measures, always in pairs:
- copy recall = expressed copies with a transcript; copy precision = called copies that are expressed;
- isoform recall = planted isoforms reproduced exactly at their copy; isoform precision = transcripts on planted copies
  that are such an isoform;
- joined transcripts (two copies), transcripts elsewhere; at read level: placed on the true copy / another / tied.

## Result on planted reads (row 1903; 133 copies, 239 isoforms, CHM13, the usual alignment)
Copy recall / copy precision / isoform recall / isoform precision.

| tool | E: all expressed (7,980 reads) | D: one dominant (11,360) | S: 63 of 133 silent (4,200) |
|---|---|---|---|
| Rustle, primaries only | 0.970 / - / 0.766 / 0.803 | 0.910 / - / 0.669 / 0.788 | 0.971 / 0.739 / 0.719 / 0.575 |
| Rustle, near-tie secondaries | 0.970 / - / 0.791 / 0.677 | 0.947 / - / 0.724 / 0.603 | 0.971 / 0.624 / 0.750 / 0.412 |
| Rustle, multi-copy mode | 0.970 / - / 0.791 / 0.677 | 0.970 / - / 0.741 / 0.592 | 0.971 / 0.618 / 0.750 / 0.409 |
| StringTie -L | 0.970 / - / 0.649 / 0.791 | 0.977 / - / 0.636 / 0.776 | 0.986 / 0.734 / 0.648 / 0.606 |
| FLAIR | 0.887 / - / 0.661 / 0.790 | 0.887 / - / 0.640 / 0.712 | 0.871 / 0.772 / 0.648 / 0.656 |
| IsoQuant | 0.955 / - / 0.615 / 0.803 | 0.955 / - / 0.611 / 0.785 | 0.943 / 0.750 / 0.594 / 0.628 |
| bambu (NDR 1) | 0.970 / - / 0.774 / 0.661 | 0.970 / - / 0.715 / 0.684 | 0.971 / 0.739 / 0.719 / 0.503 |

Silent copies given a transcript in S (of 63): mode 42, Rustle with secondaries 41, StringTie 25, Rustle primaries 24,
bambu 24, IsoQuant 22, FLAIR 18. By family the mode's extra false copies are NPIP (4 against 0), TBC1D3 (4 against 0),
TSPY (16 against 12), GOLGA8 (3 against 1).

What it says:
- **Every tool reaches nearly every expressed copy** (0.87-0.99): the aligner's toss between identical copies spreads
  primaries over all of them. The mode has no copy-recall advantage on planted reads; StringTie is one copy ahead in D
  and S.
- **Every tool reports silent copies** (18-25 of 63 with primaries alone), because that same toss puts primaries on
  silent identical copies. Near-tie secondaries add 17 more.
- **Rustle reproduces the most planted isoforms** (0.74-0.79 against 0.61-0.72 for the other tools except bambu), and
  with secondaries pays for it in isoform precision (0.41-0.68 against 0.61-0.79 for StringTie).
- Registered claims: C1 (mode's copy recall >= every tool) FAILS in D and S by one copy; C2 (copy precision within 0.05
  of StringTie in S) FAILS, 0.618 against 0.734; C3 holds on isoform recall and FAILS on isoform precision.

## A second genome: gorilla OR6737, its own families (unregistered, the driver's `o2bench` check, 2026-11-22)
Planted at the 1,328 copies of the run's own 379 families (the run's own transcripts as isoforms, 30 reads per copy; no
annotation). Copy recall / copy precision / isoform recall / isoform precision:

| tool | E: all expressed | S: 616 of 1,328 silent |
|---|---|---|
| Rustle, primaries only | 0.853 / - / 0.723 / 0.912 | 0.836 / 0.895 / 0.693 / 0.814 |
| Rustle, near-tie secondaries | 0.895 / - / 0.776 / 0.884 | 0.882 / 0.724 / 0.737 / 0.658 |
| Rustle, multi-copy mode | 0.907 / - / 0.789 / 0.871 | 0.892 / 0.721 / 0.750 / 0.638 |
| StringTie | 0.888 / - / 0.853 / 0.972 | 0.889 / 0.892 / 0.859 / 0.893 |
| FLAIR | 0.950 / - / 0.893 / 0.975 | 0.938 / 0.892 / 0.890 / 0.904 |
| IsoQuant | 0.845 / - / 0.818 / 0.985 | 0.846 / 0.907 / 0.827 / 0.916 |
| bambu | 0.889 / - / 0.904 / 0.730 | 0.892 / 0.892 / 0.906 / 0.687 |

Here Rustle leads on nothing: FLAIR has the highest copy recall, and every other tool reproduces more of the planted
isoforms than Rustle, although the planted isoforms are Rustle's own models. The isoform-recall lead of row 1903 (human
named families, annotated isoforms, 60 reads per copy) does not carry over to this set (single-exon and low-count copies,
30 reads). One run, one seed, not registered: a result to explain, not to quote as a ranking.
`r1903_o2_family_benchmark/or/g.o2bench.tsv`.

## Can evidence decide which copies to keep? (rows 1904, 1905; `PREREG_copy_verdicts_2026-11-22.md`)
The advisor's point: other tools ignore the near-tied multimappers or take them all; ours should use them with evidence.
What the planted reads say:
- **Per read there is nothing to decide in an exact tie.** Exact ties are 41% of the planted reads; the primary sits on the
  true copy for 31-41% of them (a toss among 3-5 copies); fewest edits single out the true copy in 19 of 3,309; the legacy
  `assign` stage assigns 1 of 3,125. Near ties (rival 1-10 `ms` below) are 96-99% right as they are.
- **So the decision is per copy.** Half silent (63 silent, 70 expressed), silent copies given a transcript: ignoring or
  tossing 18-25 (every other tool), taking all near ties 42 (the mode).
- **Pairwise evidence fails** (row 1904, `tools/copy_verdicts.py`): 42 -> 40. A silent copy that matches two expressed
  copies at different places is never refuted by either.
- **The joint test works halfway** (row 1905, `bench/o2_families/em_filter.py`: keep a locus only if the EM leaves it
  >= 2 reads): 42 -> 20, copy precision 0.767 (StringTie 0.734, FLAIR 0.772), isoform recall still 0.74 against 0.59-0.72;
  but it loses 2 expressed copies in S and 10 in D (low-count copies beside a dominant one). Not validated, not a default.

## The cases one by one: identical copies, recently diverged copies, the Y families, IsoCon (row 1906)
`PREREG_yag_cases_2026-11-22.md`. Copies with identical planted isoforms form one UNIT (no read can separate them, for
any method): the 133 copies are 108 units (TSPY 35 copies / 20 units, CT47A 7 / 3, BPY2 3 / 1, PRY 2 / 1). Scoring by
unit removes the calls nobody could avoid; divergence classes show the recently diverged copies directly.
Fresh seed, half silent (62 expressed units, 46 silent):

| tool | units found | false units | unit precision | silent copies 1-5 edits from an expressed one, called (of 14) |
|---|---|---|---|---|
| ignore near ties (Rustle, primaries) | 58 | 12 | 0.829 | 6 |
| take all near ties (Rustle mode) | 59 | 30 | 0.663 | 14 |
| **mode + EM test** | 59 | 4 | 0.937 | 1 |
| StringTie | 59 | 12 | 0.831 | 7 |
| FLAIR | 58 | 12 | 0.829 | 6 |
| IsoQuant | 58 | 10 | 0.853 | 6 |
| bambu | 59 | 12 | 0.831 | 6 |
| IsoCon | 62 | 9 | 0.873 | 4 |

- The near ties used with the statistical test give the fewest false units of all tools at StringTie's unit recall, and
  call 1 of 14 recently diverged silent copies where the others call 4-7 and taking all calls 14 (claims Y1, Y3 hold).
- IsoCon finds the most units (62 of 62; 107 of 108 with all expressed) and every expressed sequence exactly; on the Y
  families it finds 26 of 26 units against 24 for mode + EM and 25 for StringTie (claim Y2 fails on recall).
- The cost of the EM test seen at copy level (row 1905) is one unit at unit level.
- Joined transcripts: 1-6 per tool, no difference.
`bench/o2_families/units.py PLANT_DIR GENOME.fa`, `units_score.py PLANT_DIR SC [--families ...] LABEL=GTF ...`,
`isocon_arm.sh`, `em_filter.py`.

## Replication of the EM test (row 1907): five more seeds and gorilla
Half silent, false units: mode + EM 2, 2, 4, 2, 2; the best other genome-based tool 9, 7, 12, 11, 6; the mode alone
26-36; IsoCon 7-12. Unit precision 0.93-0.97 against 0.80-0.91. Recently diverged silent copies called: 0-1 against 3-6
for StringTie. Gorilla OR6737's own families (1,227 units): 6 false units against StringTie's 13 and the mode's 178, at
607 of 684 units found (StringTie 609; FLAIR 682 with 40 false).
Cost: one to two units fewer than the mode in a scenario (two on three of the five seeds), which the registered rule did
not allow: by that rule the test is NOT validated and the extension stays opt-in. With labelling instead of dropping the
cost is one to two expressed units tagged `homology_only`.
Real gorilla run of the wired extension: 20,120 loci labelled expressed 18,868, shared 46, homology_only 1,206.

## Validation of the default: the labels (row 1911)
The default drops nothing and labels. Five fresh seeds, half silent: units with a transcript 84-89 against 68-74 for the
best other genome-based tool; expressed units wrongly labelled `homology_only` 0-2 of 59-62; silent units with a
transcript labelled `homology_only` 79-93% (four of five seeds at or above the 80% bar). Gorilla OR6737's own families:
786 units available against FLAIR 722 and StringTie 622; 1 of 608 expressed units mislabelled; 172 of 178 silent ones
labelled `homology_only`. Validated by the rule fixed before the data (`label_score.py PLANT_DIR SC LABELLED.gtf ...`).

## Depth: where the near ties add copies (row 1908)
Every run above gave each copy 20-60 reads, enough for the aligner's toss to leave primaries on every identical copy.
With fewer reads the toss leaves copies short. Expressed copies with a transcript (of 133; identical copies count each,
as a read that fits several identical loci belongs to all of them):

| reads per copy | Rustle, primaries | StringTie | FLAIR | IsoQuant | bambu | Rustle mode | mode + EM test |
|---|---|---|---|---|---|---|---|
| 3 | 44 | 91 | 52 | 76 | 77 | 121 | 115 |
| 5 | 63 | 117 | 106 | 112 | 113 | 128 | 126 |
| 10 | 97 | 124 | 114 | 125 | 126 | 126 | 124 |
| 20 | 117 | 126 | 118 | 127 | 129 | 129 | 127 |
Half silent, expressed copies found of 70 / false units: 3 reads: primaries 13 / 0, StringTie 40 / 2, mode 57 / 27,
mode + EM 41 / 0; 5 reads: 29 / 1, 59 / 2, 65 / 30, 60 / 1; 10 reads: 44 / 2, 64 / 5, 66 / 31, 66 / 2; 20 reads: 61 / 5,
66 / 7, 66 / 30, 66 / 2.

At 3 reads per copy the mode finds 30 more expressed copies than StringTie and 77 more than Rustle on primaries alone;
with the EM test 24 more than StringTie, with 0-2 false units at every depth (the mode alone: 27-31). The lead is gone by
10 reads per copy. This is the regime of real samples and the likely reason the mode models more annotated copies there.

## Available and labelled, not dropped (user and advisor, 2026-11-22)
The advisor's position: when one member of a family is found, the others should be available too; the tool should not
declare a member "not expressed" and hide it. Since then the extension LABELS instead of dropping
(`em_filter.py --label`, the driver's default under `RUSTLE_O2_EXTEND=1`; `RUSTLE_O2_DROP=1` = the dropping form scored in
rows 1905-1908): every transcript stays, tagged `o2_status` = `expressed` (reads placed there alone), `shared`
(indistinguishable from its partners on these reads) or `homology_only` (reads reach it, all better explained elsewhere).
Seed 1906, half silent, members with a transcript of 133 (expressed + silent): the mode 114 (67 + 47), StringTie 95,
bambu 95, Rustle primaries 92, IsoQuant 92, IsoCon 91, FLAIR 82. Of the mode's 47 silent members, the label says
`homology_only` for 28 and keeps 19, of which 17 are identical to an expressed copy and 2 are wrong; of its 67 expressed
members 1 is labelled `homology_only`. So availability and the expression label are two outputs: the planted truth
scores the label, not the presence. Not claimed: that every member of a family is expressed when one is (row 1870, real
A119b reads: CT47A1-4 have 125-131 reads each and none supports them at the bases that distinguish them).

## Silent relatives are visible in the low-scoring secondaries (row 1909, `PREREG_silent_members_2026-11-22.md`)
| silent loci, reached from the reads of their expressed relative | n | within 2% of the best | >= 0.9 | >= 0.8 | >= 0.5 | any score | slim mode `-p 0.8`, any score |
|---|---|---|---|---|---|---|---|
| processed pseudogenes, from their gene's reads | 417 | 7 | 14 | 31 | 119 | 292 (70%) | 20 |
| genes, from their pseudogenes' reads | 150 | 1 | 2 | 7 | 64 | 122 (81%) | 13 |
| silent copy of a tandem pair | 150 | 49 | 74 | 82 | 118 | 141 (94%) | 84 |
| silent copy of a dispersed pair | 128 | 10 | 35 | 54 | 82 | 96 (75%) | 53 |

The advisor's statement holds on planted reads: an expressed gene's reads reach 70% of its silent processed pseudogenes
and expressed pseudogenes reach 81% of their silent gene, almost never through near ties. These records are what the
slim mapping mode discards. Specificity is not measured (over a thousand other places are reached too).

## The member-discovery layer (row 1910, `PREREG_member_discovery_2026-11-22.md`)
`tools/rustle_pipeline.sh members --bam B --fasta G --out PREFIX` (opt-in; `tools/member_discovery.py`): for every
assembled locus, the places its reads also align to, at any score down to 0.3 of the read's best, with >= 5 reads ->
`PREFIX.members.tsv` and `PREFIX.members.gtf` (models tagged `member_candidate`). Where relatives are, expressed or not.
On planted reads with every relative silent and a control set of 150 genes without relatives (second draw, thresholds
fixed on the first): 6% of control genes get a place; 90% of the 1,350 places lie on a listed relative (969) or an
annotated pseudogene (246); reached: 72% of recent paralogues, 48% of named processed pseudogenes, 1% of older
paralogues. Near ties alone reach a third as many relatives. Needs the usual `-p 0.1` mapping.
On a real sample it is NOT validated (row 1912, gorilla OR6737 against RefSeq): 53-58% of its places lie on an annotated
paralogue or pseudogene, 42-47% on another gene or unannotated sequence (planted: 10%); 21% of recent paralogues of
expressed genes are reached (planted: 72%); control genes get places as often as genes with relatives, the control being
unable to exclude pseudogene parents. Treat `PREFIX.members.gtf` as candidates to check, not as members.

## Result on real samples (rows 1859-1869; annotated copies with a spliced transcript / junction precision)
| sample | Rustle, primaries | Rustle, secondaries | mode | StringTie | FLAIR | IsoQuant |
|---|---|---|---|---|---|---|
| four ape testis samples together | 480 / 0.785 | 618 / 0.760 | 701 / 0.732 | 587 / 0.699 | 492 / 0.561 | 463 / 0.801 |
| A119b, human testis | 93 / 0.57 | 123 / 0.49 | 126 / 0.49 | 107 / 0.47 | 98 / 0.35 | 85 / 0.65 |
| OR6737, gorilla testis | 25 / 0.61 | 27 / 0.58 | 34 / 0.66 | 28 / 0.59 | 23 / 0.55 | 23 / 0.56 |

## How the two arms fit together
On real samples the mode models more annotated copies than any other tool. The planted arm shows what that number can
and cannot mean: where copies are identical, a tool that follows near-tie secondaries reports the silent ones too
(42 of 63 against 25 for StringTie), and no tool separates them. So "more copies modeled" is supported; "more expressed
copies found" is not, and on planted reads it is false. What the planted arm does support for Rustle is isoform recall.

## Wording
- Supported: "On known multi-copy families Rustle reproduces more of the planted isoforms than StringTie, FLAIR or
  IsoQuant (0.74-0.79 against 0.59-0.66), and on real samples its multi-copy mode has transcripts at more annotated
  copies."
- Not supported: "Rustle finds more expressed copies"; "the extra copies are expressed"; any copy-recall advantage on
  planted reads.
- Always give copy precision in S beside any copy count, and isoform precision beside isoform recall.

## Limits
- Planted reads: the sample equals the reference (no private variants, no copy the reference lacks), simple error and
  truncation models, three scenarios of our choosing, 60 reads per copy; a scenario, not a mimic of real data
  (rows 1797, 1809). One seed.
- Isoform truth = annotated isoforms; "exact" = the intron chain, ends not scored.
- Human named families only (`plant.py` reads the CHM13 loader); other species need a member list.
- Real-sample arm: annotated copies only, any spliced transcript counts, no truth on expression.

## Related rows
1855-1869 real samples per tool and the mode; 1871 tie groups and counts; 1873, 1874, 1877 sharing tied reads with known
origin; 1875 gene + pseudogene, tandem and dispersed pairs; 1881 joins between tandem copies; 1885, 1888, 1894, 1902
aligner settings and the slim mapping mode. `MULTICOPY_MODE.md`, `PREREG_tied_truth_2026-11-21.md`,
`PREREG_tandem_joins_2026-11-21.md`, `PREREG_minimap_p09_n1000_2026-11-22.md`, `PREREG_o2_family_benchmark_2026-11-22.md`.
Workspace `r1903_o2_family_benchmark` (`plant/`, `{E,D,S}/score.tsv`, `run.sbatch`).
