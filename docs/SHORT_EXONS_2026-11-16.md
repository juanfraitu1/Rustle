# Alternative splicing classes, retained introns and short exons: what was measured and the `shortexons` stage, state on 2026-11-16

Summary of the sessions of 2026-11-15 / 16 that started from the user question "are we modelling correctly all
alternative splicing events in comparison with those expected by the annotation? Would using vg help". Register rows
1742-1746, 1753, 1754, 1757. Figures are copied from the pre-registrations named in each section, which hold the bars
and the full tables.

## Outcome in five lines

- Cassette exons, mutually exclusive exons and exons of 16-30 bp are emitted as often as the control class. Two gaps
  were found: retained introns and exons of 15 bp or less.
- Retained introns: the retained-intron filter (`--polish-retained-ratio 10`) is the whole cause; relaxing it fails
  the cost bar. It stays at 10 (user, 2026-11-15); `RUSTLE_POLISH_RETAINED_RATIO` is the opt-in switch.
- Exons of 15 bp or less: minimap2 drops them (the skip reads carry the exon's bases 74-96% of the time on eight
  samples; from 16 bp up, under 2% on the four samples where that was read). This is the one splicing case where a graph-style step has a measured target.
- Two annotation-free rules rewrite those reads through the exon. Together they reach 70-89% of the dropped annotated
  exons; where they reach, false skip chains nearly vanish. Whole-assembly precision does not move consistently.
- Shipped as the opt-in stage `shortexons` (user, 2026-11-16). Not a default: it would be promoted only if another
  rule removes the wrong new chains it lets in. No bar defines that yet.

## 1. Event classes (row 1744, `PREREG_alt_events_missing_classes_2026-11-15.md`)

Share of annotated pairs with both forms emitted, `bench/famcmp/alt_events.py`: cassette exons <= 30 bp and mutually
exclusive exons at the cassette control; retained introns 0.744 (OR6737) / 0.682 (A119b) against 0.975 / 0.912 (FLAIR
0.953 / 0.906). Short-exon alignment loss was not flagged as registered and is concentrated at <= 15 bp post hoc.

## 2. Retained introns (row 1745, `PREREG_retained_intron_rule_2026-11-15.md`)

| sample | filter | retained-intron pairs, both emitted | RefSeq chain precision | chains added (RefSeq-correct) |
|---|---|---|---|---|
| OR6737 | 10 (default) | 0.744 | 0.3517 | |
| OR6737 | 30 | 0.953 | 0.3458 | 1,458 (71) |
| OR6737 | off | 1.000 | 0.3388 | 3,033 (87) |
| A119b | 10 (default) | 0.682 | 0.1764 | |
| A119b | 30 | 0.812 | 0.1722 | 6,641 (236) |
| A119b | off | 0.932 | 0.1660 | 15,916 (355) |

Attribution PASS (filter off puts the class at the control on three samples), cost FAIL for both arms (30 misses the
0.005 bar by 0.0009 on OR6737). RefSeq cannot say which unannotated retained-intron chains are real, so the cost is an
upper bound.

## 3. Exons of 15 bp or less

| row | question | result |
|---|---|---|
| 1746 | do reads that skip an annotated short exon carry its bases? (`skip_realign.py`) | yes: 0.853 / 0.768 / 0.794 at <= 15 bp (OR6737 / A119b / KB3781), 0.018 / 0.017 / 0.000 at 16-30 bp |
| 1753 (pre-registered as 1751) | junction rule: rewrite the read through an exon another read of the locus spells (`skip_fix.py`) | FAIL on reach, 0.26-0.52: most dropped exons are dropped in every read |
| 1754 | intron search: match the read's extra junction bases to one splice-flanked site of the intron, used only when sites x 4^-L <= 0.01 (`skip_search.py`) | FAIL on precision (HG002 0.908 vs 0.95; single-read finds 2-3 bp off); reach with the junction rule 0.74-0.89 |
| 1757 | the same with a floor of 2 reads, on PTR / PPA / PPY / PAB | NOT TESTABLE on precision (one sample with 20 judged units; pooled 68 of 71 post hoc); reach 0.70-0.82 |

The two rules are complementary. The search finds exons of 4-10 bp that every read drops; the junction rule covers
11-15 bp, where the reads that drop the exon are the noisy ones and fail the search's one-edit limit. Exons of 3-5 bp
in long introns stay out of reach (chance matches). Decoys (made-up exons, reversed introns) fire at 0.0001-0.005.

Effect on the assembly, panel of reads within 100 kb of a fired exon (rows 1754 and 1757; decides nothing as
registered):

| | OR6737 | KB3781 | A119b | HG002 | PTR | PPA | PPY | PAB |
|---|---|---|---|---|---|---|---|---|
| skip chains where only the exon form is annotated | 72 -> 0 | 57 -> 0 | 210 -> 20 | 112 -> 13 | 51 -> 1 | 47 -> 2 | 39 -> 2 | 27 -> 3 |
| chains with the exon there | 45 -> 107 | 66 -> 116 | 140 -> 326 | 104 -> 200 | 25 -> 74 | 31 -> 73 | 51 -> 85 | 23 -> 48 |
| new / lost RefSeq-correct chains | 14 / 20 | 14 / 14 | 96 / 18 | 73 / 11 | 19 / 8 | 11 / 7 | 6 / 16 | 8 / 5 |
| panel precision | -0.0016 | -0.0006 | +0.0017 | +0.0017 | +0.0049 | -0.0010 | -0.0016 | +0.0005 |

Rows 1754's four samples ran without the read floor, row 1757's four with it. Most lost correct chains are annotated
skip forms at exons where both forms are annotated: their reads moved to the exon form and the skip chain fell under
the emission threshold. The new chains are 3-28% RefSeq-correct; those are the "wrong calls" a promotion would have to
remove.

## The stage

| product | status | command |
|---|---|---|
| `PREFIX.shortexons.gtf` (+ `.final.gtf`): default assembly of the reads after the rewrite; `PREFIX.shortexons.units.tsv` (exon, rule, read counts); `PREFIX.shortexons.reads.tsv` (read, old and new CIGAR) | opt-in, named only, never part of `all` | `bash tools/rustle_pipeline.sh shortexons --bam B --fasta G --out PREFIX [--threads N]` |

- `bench/famcmp/short_exons.py`; needs python3 with pysam. No annotation. `PREFIX.gtf` is never changed.
- The inner assemble is the default one: `RUSTLE_*` switches apply, command-line assemble options do not.
- While it runs, `PREFIX.shortexons.d` holds a second copy of the BAM (removed at the end and on failure;
  `RUSTLE_SHORTEXONS_KEEP=1` keeps it). Unmapped reads are left out of that copy.
- Checked through the driver: PTR gives the same 294 units and 3,280 rewritten reads as row 1757 (151 s, 4.6 GB BAM);
  A119b 1,477 units, 29,989 reads, 2,324 s on 72 cores, work directory peak about 180 GB (96 GB BAM), 529 new chains
  (96 RefSeq-correct), 449 lost (18 correct) against an earlier default assembly. Peak memory was not captured.

## Not done

- No rule for the wrong new chains, and no bar for "promotable".
- Exons of 3-5 bp; reads whose skip junction the aligner also shifted; supplementary and secondary alignments.
- A fibroblast or cell-line sample with the read floor (all four floor samples are testis).
- vg itself was not run: the junction rule is the graph use (a path other reads show), and it reaches the minority.

## Also in these sessions

- `trunc5` against CAGE (rows 1742 INVALID, 1743 valid with 1 kb flank backgrounds): the stage drops a CAGE-supported
  minority (0.048 A119b / 0.125 HG002 of dropped starts in a peak, against 0.27 / 0.49 kept): tag, not drop.
- `capstarts` refuses both gorilla libraries because of primer trimming: OR6737's 5' primer ends in GGGG, so the cap G
  is trimmed away; KB3781 keeps one constant G and fails the stage's self-check.
