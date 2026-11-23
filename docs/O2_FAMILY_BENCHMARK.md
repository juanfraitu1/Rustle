# O2: known gene families across tools, on planted and real reads

## Current state
(2026-11-23) O2 was replaced on the user's word (register 1903) and confirmed complete as it stands (2026-11-22).
It is no longer per-read copy assignment (`assign`, `copy_assign`'s assignment mode; record in `LEDGER_O2.md`, kept as
dormant code that runs only when called by name). O2 asks, for every assembler on the same alignments: **which copies and
isoforms of known multi-copy families does it report, and are they the expressed ones?**

- **Default pipeline.** `tools/rustle_pipeline.sh all` = assemble, families, **o2**, flag. Stage `o2` extends the assembly
  to the multi-copy families with the near-tied multimappers (stage `multicopy`), labels every locus by a provenance
  test (`o2_status` expressed / shared / homology_only; nothing is dropped) -> `PREFIX.o2.gtf`, `PREFIX.o2.em.tsv`, and
  writes the per-family report `PREFIX.o2.families.tsv`, `PREFIX.o2.copies.tsv`. `RUSTLE_O2_EXTEND=0` = report only;
  `RUSTLE_O2_DROP=1` drops the homology_only loci; a failure of the extension leaves the report on `PREFIX.gtf` and the
  stage still succeeds. `PREFIX.gtf` is never changed. `assign` is in no plan of `all`.
- **Validation.** Planted reads of known origin, opt-in stage `o2bench` and `bench/o2_families/` (planting is the
  validation, not a per-run step: user).
- **Committed** 074557f7 (2026-11-22). Since then, uncommitted: ragged read ends in `plant.py`, `member_discovery.py`
  rebuilt, `related_gtf.py`, the records of rows 1913-1918.

What holds, on planted reads with ragged ends (rows 1917, 1918) unless said:

| claim | status |
|---|---|
| Rustle reproduces the most isoforms of the genome-based tools | holds at every depth: 0.82 against 0.66-0.82 at 60 reads per copy, 0.67 against 0.31 at 3; gorilla 0.96-0.98 against 0.85-0.94 |
| Near ties add expressed copies | only at about 3 reads per copy (125 of 133 against StringTie 95); no lead from 5 reads up |
| Taking every near tie is wrong | 24 false units (half silent); ignoring them leaves 9-16 for every tool |
| The EM provenance test fixes that | 3 false units; recently diverged silent copies called 1 against StringTie 9 |
| Labels make more members available, accurately | 86-92 units against 75-76; 0-2 expressed units mislabelled; ten seeds and gorilla |
| Rustle finds every copy with ordinary introns | 267 of 267 (IsoCon 245) |
| Rustle finds the most copies overall | does NOT hold: it misses single-exon models and lift-over models with non-canonical introns |
| Low-scoring secondaries reveal silent relatives | holds on planted reads (70-81%); as a member list on real data not usable yet |

Everything above is on planted reads where the sample equals the reference.

## How to run
```
tools/rustle_pipeline.sh all     --bam B --fasta G --out PREFIX     # ... families, o2, flag
tools/rustle_pipeline.sh o2      --bam B --fasta G --out PREFIX     # extension + labels + family report
tools/rustle_pipeline.sh o2bench [--gff ANNOT.gff3.gz] --bam B --fasta G --out PREFIX   # planted benchmark (opt-in)
tools/rustle_pipeline.sh members --bam B --fasta G --out PREFIX     # relatives through all secondaries (opt-in)
```
`o2bench` with `--gff`: the named human families (famtruth environment, CHM13); without: the run's own families, any
genome. The pieces (`bench/o2_families/`):
```
plant.py ANNOT.gff3.gz GENOME.fa DIR [--seed N --reads 60 --isoforms 2 --scenarios E,D,S --ragged 20]
plant.py --from-run PREFIX.gtf PREFIX.fam.copies.tsv GENOME.fa DIR
run_tools.sh DIR SC GENOME.fa OUT [THREADS]      # one BAM, every tool, score.tsv, units.tsv, reads.txt
score.py DIR SC [--families] LABEL=GTF ...       # copy and isoform recall / precision
units.py DIR GENOME.fa ; units_score.py DIR SC [--families F1,F2] LABEL=GTF ...
em_filter.py ASSEMBLY.gtf BAM AS_TABLE OUT [--label]   ; label_score.py DIR SC LABELLED.gtf OTHER=GTF ...
isocon_arm.sh DIR SC GENOME.fa OUT ; tied_share.py LOCI.gtf BAM AS_TABLE TRUTH.tsv OUT
```
Scenarios: E every copy `--reads`; D one copy 10 x, the others a third; S half the copies of a family silent.
`run_tools.sh` runs FLAIR, IsoQuant, bambu when `O2_FLAIR_PATH`, `O2_ISOQUANT_MANIFEST`, `O2_BAMBU_MANIFEST` are set
(Fritz values in any `r19xx/*.sbatch`). Self-tests: `plant.py`, `score.py`, `tools/o2_family_report.py`,
`tools/member_discovery.py`, `tools/copy_verdicts.py --selftest`.

## Rules for reading a planted result (each learned the hard way)
- **Ragged ends.** Until row 1917 the full-length reads of a transcript shared both ends; Rustle counts coordinate
  duplicates once toward its 2-read floor (row 1850) and lost a third of the minor isoforms (found 32% with one distinct
  start-end pair, 99.6% with two). `plant.py` trims 0-20 bases off each end by default. Rows 1903-1911 are superseded for
  every Rustle-against-other-tool number; their findings on false copies, the EM test and the labels repeated.
- **Units.** Copies whose planted isoforms are identical form one unit (133 human copies = 108 units; TSPY 35 -> 20):
  no read separates them for any method. Count false calls by unit; count recall per copy when the question is
  availability (a read that fits several identical loci belongs to all of them: user).
- **Transcript classes.** 18 of the 133 human copies come from lift-over models with non-canonical or zero-length
  introns that real reads mostly do not show (3 of 18 on A119b), 9 are single-exon. Give the canonical class beside
  totals.
- **Family-only BAMs distort quantile rules.** Rustle's single-exon floor is the 0.82 quantile of the multi-exon read
  counts of the contig; where tied reads count at every copy it reaches hundreds of reads.
- **Sample = reference.** No private variant, no copy the reference lacks; a scenario, not a mimic (rows 1797, 1809).
- Fresh seed for every confirmation; a seed that has been looked at is development data.

## Results on planted reads, ragged ends

### All tools (row 1917; copy recall / isoform recall / isoform precision)
| tool | human named families, all expressed | human, half silent | gorilla own families, all expressed |
|---|---|---|---|
| Rustle, primaries only | 0.970 / 0.820 / 0.800 | 0.957 / 0.806 / 0.592 | 0.978 / 0.982 / 0.950 |
| Rustle mode (all near ties) | 0.970 / 0.824 / 0.629 | 0.957 / 0.815 / 0.449 | 0.959 / 0.964 / 0.909 |
| mode + EM test | 0.962 / 0.820 / 0.641 | 0.943 / 0.798 / 0.576 | 0.956 / 0.960 / 0.912 |
| StringTie | 1.000 / 0.669 / 0.773 | 1.000 / 0.661 / 0.582 | 0.883 / 0.848 / 0.973 |
| FLAIR | 0.992 / 0.816 / 0.741 | 1.000 / 0.734 / 0.532 | 0.999 / 0.942 / 0.958 |
| IsoQuant | 0.940 / 0.590 / 0.815 | 0.914 / 0.565 / 0.569 | 0.841 / 0.814 / 0.984 |
| bambu | 0.970 / 0.745 / 0.622 | 0.957 / 0.710 / 0.458 | 0.888 / 0.902 / 0.728 |

On gorilla at 30 reads per copy Rustle with primaries only is slightly ahead of the mode (the secondaries add noise
there); with secondaries its isoform precision is the weak point on the human set.

### Depth (row 1918, seed 1918; expressed copies of 133, isoform recall in brackets)
| reads per copy | Rustle primaries | Rustle mode | mode + EM | StringTie | FLAIR | IsoQuant | bambu |
|---|---|---|---|---|---|---|---|
| 3 | 81 (0.34) | 125 (0.67) | 120 (0.65) | 95 (0.31) | 45 (0.18) | 79 (0.31) | 80 (0.32) |
| 5 | 115 (0.60) | 128 (0.74) | 117 (0.67) | 121 (0.52) | 107 (0.43) | 110 (0.48) | 114 (0.51) |
| 10 | 129 (0.73) | 129 (0.80) | 127 (0.79) | 132 (0.64) | 126 (0.64) | 125 (0.57) | 128 (0.69) |
| 20 | 129 (0.80) | 129 (0.82) | 127 (0.81) | 132 (0.66) | 131 (0.78) | 126 (0.60) | 129 (0.70) |

Half silent at 3 reads: expressed copies found the mode 67 of 70 with 22 false units, mode + EM 52 with 0, StringTie 39
with 1.

### Units, recently diverged copies, IsoCon (row 1918, seed 1919; half silent: 62 expressed units, 46 silent)
| tool | units found | false units | unit precision | silent copies 1-5 edits from an expressed one, called |
|---|---|---|---|---|
| Rustle, primaries only | 60 | 11 | 0.845 | 5 |
| Rustle mode | 60 | 24 | 0.714 | 17 |
| mode + EM test | 60 | 3 | 0.952 | 1 |
| StringTie | 62 | 13 | 0.827 | 9 |
| FLAIR | 62 | 16 | 0.795 | 13 |
| IsoQuant | 57 | 9 | 0.864 | 5 |
| bambu | 59 | 11 | 0.843 | 5 |
| IsoCon | 62 | 12 | 0.838 | 8 |

Y-ampliconic families (25 expressed units): mode + EM 23 with 3 false, StringTie 25 with 8, IsoCon 25 with 9.
Expressed copies found by the class of the planted transcript (three scenarios pooled):

| class | copies | Rustle | StringTie | FLAIR | bambu | IsoQuant | IsoCon |
|---|---|---|---|---|---|---|---|
| canonical introns | 267 | 267 | 267 | 264 | 267 | 262 | 245 |
| single exon | 22 | 17 | 20 | 21 | 17 | 12 | 22 |
| non-canonical or 0-length intron | 47 | 40-42 | 47 | 47 | 41 | 41 | 46 |

IsoCon has no model of introns, so splice-site and single-exon rules cannot remove a copy: that is its whole lead. On
ordinary copies it finds fewer, merges copies one edit apart, and is the weakest tool when expression is uneven
(115 of 133 copies). Its reference-free side is not exercised here.

### Labels, the default (rows 1911, 1918)
Five fresh seeds with ragged ends, half silent: units with a transcript 92 / 90 / 88 / 91 / 86 against 75-76 for the best
other genome-based tool; expressed units labelled `homology_only` 0-1 (1-2 with all copies expressed); silent units with
a transcript labelled `homology_only` 0.91, 0.81, 0.79, 0.91, 0.82. Gorilla OR6737's own families: 860 units available
against FLAIR 725; 1 of 664 expressed units mislabelled; 193 of 196 silent ones labelled right. Validated by its
pre-registered rule twice (row 1911 on the old reads, row 1918 on the corrected ones).
With all copies expressed the labelled extension has 104 of 108 units, StringTie and FLAIR 108.

### What a tied read can decide (row 1904 step 0, old reads; not re-run)
Exact ties are 41% of planted reads; the primary sits on the true copy for 31-41% of them; fewest edits single out the
true copy in 19 of 3,309; the legacy `assign` stage assigns 1 of 3,125. Near ties (rival 1-10 `ms` below) are 96-99%
right as placed. So the decision is per copy, not per read. Pairwise decisive-read verdicts (`tools/copy_verdicts.py`)
fail (42 -> 40 false copies): a silent copy that matches two expressed copies at different places is refuted by
neither. The joint EM test (`em_filter.py`) is what works.

## Silent relatives and the member layer
**Low-scoring secondaries carry them** (row 1909, planted): a silent locus reached by >= 2 reads of its expressed relative.

| silent loci | n | within 2% of the best | any score | slim mapping `-p 0.8` |
|---|---|---|---|---|
| processed pseudogenes, from their gene's reads | 417 | 7 | 292 (70%) | 20 |
| genes, from their pseudogenes' reads | 150 | 1 | 122 (81%) | 13 |
| silent copy of a tandem pair | 150 | 49 | 141 (94%) | 84 |
| silent copy of a dispersed pair | 128 | 10 | 96 (75%) | 53 |

**The layer** (`members`, `tools/member_discovery.py`; rows 1910-1915): for every assembled locus the places its reads'
secondaries also align to (>= 5 reads, AS >= 0.3 of the best), a place = hits whose aligned blocks overlap.
- Planted (row 1915): 6-9% of control genes without relatives get a place; 93-95% of places on a known relative or
  annotated pseudogene; 75-77% of recent paralogues, 38-50% of pseudogenes reached.
- Real human A119b (row 1915; 54 min, 7.6 GB with `--as-table`, `--exclude` satellites): recent paralogues 64%,
  pseudogenes 56%, controls 12%, but 48% of places on a relative or pseudogene and 46% on no annotated exon; as gene
  pairs 0.23 precision / 0.74 recall. NOT usable as a member list on real data. Places still chain (median 122 kb).
- Read-derived relatedness that works: the `families` stage, 0.82 / 0.82 against Compara on A119b; union with the
  cross-mapping links 0.27 / 0.92; agreed by both 0.82 / 0.64.
- Old paralogues are out of scope: their reads place uniquely and they are separate loci (user).
- Rows 1910, 1913, 1914 (places merged by span) and 1912 (gorilla, same rule) are superseded.
**Relatedness as annotation**: `bench/denovo_eval/related_gtf.py ANNOT.slim.gff3.gz NCBI_gene_group.gz OUT` ->
workspace `annot_related/chm13.related.gtf` (`related_group`, `paralogs_recent`, `pseudogene_parent`, `gene_group`,
`copy_of`) and `.relations.tsv` (recent paralogues 4,696, older 133,398, pseudogene parents 12,826 from NCBI Gene + 245
by name, copies of one source gene 19,006). No public GTF carries this. NCBI has pseudogene parents for human only
(gorilla 7). The two largest groups are snRNA copies, not families.

## Aligner settings (rows 1885, 1888, 1894, 1902)
Only 1.8% of the secondaries of a `-p 0.1` BAM are near ties. `-N 50 -p 0.8` (the slim mapping mode,
`MAP_MODE=slim` of the workspace's `genomes/map_reads.sbatch`) gives the same named copies, family pairs and
transcripts on gorilla OR6737 with 8% of the secondaries and half the BAM. It is NOT equivalent for silent relatives
(table above). `-N 25` would cut near ties of families with more than 25 copies. Two-pass alignment changes nothing and
cannot be combined with `--eqx`; `--eqx` stays (stage `flag` and the base-level copy test read =/X); `-Y` already
writes the sequence in secondary records.

## Result on real samples (rows 1859-1869; annotated copies with a spliced transcript / junction precision)
| sample | Rustle, primaries | Rustle, secondaries | mode | StringTie | FLAIR | IsoQuant |
|---|---|---|---|---|---|---|
| four ape testis samples together | 480 / 0.785 | 618 / 0.760 | 701 / 0.732 | 587 / 0.699 | 492 / 0.561 | 463 / 0.801 |
| A119b, human testis | 93 / 0.57 | 123 / 0.49 | 126 / 0.49 | 107 / 0.47 | 98 / 0.35 | 85 / 0.65 |
| OR6737, gorilla testis | 25 / 0.61 | 27 / 0.58 | 34 / 0.66 | 28 / 0.59 | 23 / 0.55 | 23 / 0.56 |

Copies modeled, not copies shown expressed. The planted depth sweep is the likely reason for the lead: real samples
hold many copies at a few reads, where near ties seed loci that primaries leave under the floor.

## Wording
- Supported: "On planted reads of known multi-copy families Rustle reproduces more isoforms than StringTie, FLAIR,
  IsoQuant or bambu at every coverage; used with its provenance test the near-tied multimappers give the fewest false
  copies of any tool, and labelling instead of dropping makes the most family members available with their evidence."
- Supported with its condition: "near ties recover expressed copies at about three reads per copy."
- Not supported: "Rustle finds more expressed copies" in general; "the extra copies on real samples are expressed";
  "reads can be assigned to their copy of origin" for exact ties; any member list from `members` on real data.
- The advisor's position, as implemented: a found member stays available (label, not drop); every member being
  expressed when one is, is not claimed (row 1870: CT47A1-4 on A119b have 125-131 reads and none supports them at the
  bases that distinguish them).

## Single-exon copies (row 1919, `PREREG_mono_family_rescue_2026-11-23.md`)
Tried: floor-free single-exon models kept only as copies of a family of >= 2 loci (`RUSTLE_MULTICOPY_MONO=1`, opt-in).
Planted gorilla own families, expressed single-exon copies found: the mode 101 -> 119 of 147 (FLAIR 146, Rustle
primaries 124, StringTie 0); with the EM test 8 silent ones called against FLAIR's 11. Real OR6737: the 728 rescued
models are on an annotated exon 13% of the time, the default's 534 single-exon models 56%, StringTie's 51%, FLAIR's
14,850 22%. Family membership does not select real single-exon genes; not proposed. FLAIR's lead on planted single-exon
copies is bought with tens of thousands of unannotated single-exon models on real data.

## Open
- The four human copies and the single-exon class Rustle misses: single-exon copies inside families; the HSFY2-type
  fragments; whether a strongly supported non-canonical junction should pass (GOLGA6L2, GOLGA8Q on real reads).
- Isoform precision with near-tie secondaries (0.63 against 0.74-0.80).
- The EM test's cost of one to two units; the unit of `members` on real data (links per source and target locus
  instead of connected places) and what its unannotated half is.
- A planted set whose sample differs from the reference (private variants, a copy the reference lacks): the only place
  IsoCon's and any PSV-based method's real advantage could show.
- Gorilla pseudogene parents (no public table).

## History (superseded numbers, kept in the pre-registrations)
| row | what | fate |
|---|---|---|
| 1903 | seven tools on planted named families, old reads | superseded by 1917 for Rustle against other tools |
| 1904, 1905 | pairwise verdicts fail; EM test halves false copies, copy-level bars fail | stands; unit-level in 1906 / 1918 |
| 1906 | units, divergence, IsoCon, old reads | superseded by 1918 part 2 |
| 1907 | EM test replication, dropping form | superseded by the label form (1911, 1918) |
| 1908 | depth sweep, old reads | superseded by 1918 part 1 (lead only at 3 reads) |
| 1911 | labels validated, old reads | repeated in 1918 part 3 |
| 1910, 1912-1914 | member layer with span-merged places | superseded by 1915 |
| 1916 | the isoform loss traced to coinciding read ends | the correction |

Pre-registrations: `PREREG_o2_family_benchmark_2026-11-22.md` (1903), `PREREG_copy_verdicts_2026-11-22.md` (1904, 1905),
`PREREG_yag_cases_2026-11-22.md` (1906-1908, 1911, 1916-1918), `PREREG_silent_members_2026-11-22.md` (1909),
`PREREG_member_discovery_2026-11-22.md` (1910, 1912-1915), `PREREG_tied_truth_2026-11-21.md` (1873-1877),
`PREREG_tandem_joins_2026-11-21.md` (1881), `PREREG_minimap_p09_n1000_2026-11-22.md` (1885, 1888, 1894),
`MULTICOPY_MODE.md` (the mode itself). Workspace `/anvme/workspace/v138da13-rustle/r19NN_*`, one directory per row;
`annot_related/` for the relatedness annotation.

### Restoring the primaries' single-exon models (row 1920, 2026-11-23): not built
Adding to the mode the single-exon models a primaries-only assembly passes closes the planted gap (119 of 147 against the
mode's 98 and the primaries' 120; silent called 5 against FLAIR's 10) but on real libraries there is almost nothing to
restore: one model on OR6737, 224 on A119b of which 6 are annotated (0.027, bar 0.058: the registered claim fails).
The single-exon gap to FLAIR on planted reads is a benchmark figure, not lost genes. `PREREG_mono_family_rescue_2026-11-23.md`.

### Alignment scores inside the EM (row 1921, 2026-11-23): confirmed, small, opt-in
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
