# How many multi-placement reads can O2 assign, and how many does it move off the aligner's primary? 2026-10-12

**Pre-registered before the script runs.** User (2026-10-12): "how common are assignable reads to copies, do we have a lot of cases in
which it is impossible or do we have enough cases of a decent number of reads correctly assigned that StringTie and FLAIR cannot do
since they don't use secondary alignments".

**Population** (per species, benchmark BAMs `GGO_mm.bam` / `A119b.t2t.bam`, all chromosomes): molecules with another alignment within
the good window — `as_table`'s second-best AS >= 0.98 x best AS. **Sample**: 20,000 molecules per species, uniform, seed 1
(`random.Random(1).sample` over the population in table order). **Per molecule**: the primary record vs every non-supplementary
record with AS >= 0.98 x best (O2's certificate exactly as `bench/tool_difference/s1_o2.py`). Classes, first that applies:
**moved** (some rival certified over the primary: the aligner's — and so StringTie's / FLAIR's — placement is wrong and O2 can
correct it), **confirmed** (the primary certified over every rival), **identical** (no distinguishing column against any rival:
no read-level method can assign it), **undecided** (otherwise). Reported with 95% Wilson intervals and extrapolated to the
population; also split by the primary's AS rank (true tie primary = best = second; primary below the best; primary best with a close
rival), the mosaic share, and, for moved molecules, whether the winning placement is spliced where the primary is not (the
retrocopy direction). Descriptive; nothing adopted. Expected (written now): identical and confirmed dominate; moved is a small share.

**Instrument** (committed before running): `bench/tool_difference/o2_population.py` sha1 1a1180dcdf3b5358c2d254fa67a155da8c69dfe1 (2 tests).

## Result (2026-10-12)

Outputs `/mnt/linuxdisk/tmp/asm_gap/op/{gorilla,human}.o2pop.{txt,molecules.tsv}`.

| | gorilla | human |
|---|---|---|
| molecules / multi-placement population | 4,404,440 / 49,444 (1.1%) | 21,391,966 / 532,876 (2.5%) |
| **moved** (O2 certifies a rival over the aligner's primary) | 1.0% [0.9, 1.2] ~ 504 reads | 2.7% [2.5, 3.0] ~ 14,494 reads |
| **confirmed** (primary certified over every rival) | 76.0% ~ 37,565 | 62.9% ~ 335,046 |
| **identical** (no distinguishing column: unassignable by any read-level method) | 8.2% ~ 4,035 | 10.2% ~ 54,380 |
| **undecided** (too few columns) | 14.8% ~ 7,340 | 24.2% ~ 128,956 |
| true ties (primary = best = second): moved / confirmed / identical / undecided | 0.6 / 12.3 / 59.2 / 28.0% | 1.6 / 3.2 / 48.1 / 47.1% |
| primary below the best AS | 0.5 / 96.7 / 0.8 / 2.0% | 4.5 / 70.8 / 0.8 / 23.9% |
| primary best, close rival | 1.4 / 78.5 / 0.9 / 19.2% | 2.4 / 77.2 / 2.7 / 17.7% |
| mosaic (>= 2 columns for each of a pair) | 18.8% | 39.7% |
| moved where the winner has more introns (retrocopy reversal) | 8 / 204 | 24 / 544 |

**Reading.** Most multi-placement reads are assignable and the aligner's primary is already right (confirmed 76 / 63%); a further
8-10% are truly unassignable (identical copies) and 15-24% lack evidence. The reads StringTie and FLAIR misplace and O2 can move are
1.0% / 2.7% of the multi-placement population — about 500 gorilla and 14,500 human reads genome-wide (0.01% / 0.07% of all reads).
"Primary below the best AS" is confirmed 97% / 71% of the time: the aligner's primary, not the best-AS placement, is usually the
right one (the retrocopy pattern of register 1270). Mosaic reads are common (19 / 40%): pairs where the read carries each copy's
allele at >= 2 columns — the O3 regime (variants absent from both reference copies). Expected outcome (written before) held.
Open: where the moved reads land (how many copies get >= 2 of them, i.e. copies the tools cannot assemble) needs the full population.

## Amendment 1 (2026-10-12, user: "yes lets do that"; before any full-population number)

The same classification over the FULL population (`--n 0`; names processed in batches of 100,000 per BAM pass to bound memory;
identical per-molecule rule), recording for every moved molecule the winning placement (chromosome, aligned span). Aggregated by
`bench/tool_difference/o2_moved_loci.py`: each winning placement is mapped to the reference genes whose exons its aligned span
overlaps (none = unannotated); reported: genes receiving >= 1 / 2 / 5 / 10 moved reads, and for genes with >= 2, whether OP_good,
StringTie and FLAIR (benchmark GTFs) emit any transcript touching the gene's exons. A gene with >= 2 moved reads and no StringTie /
FLAIR transcript is a copy the tools miss and O2's reads could support. Descriptive; nothing adopted.

**Instruments for amendment 1** (committed before the full run; the batched script reproduces the 20,000 sample exactly with 3 batches): `o2_population.py` sha1 f8d7435e304215df6c0c510c87179c3140a1ebec, `o2_moved_loci.py` sha1 82bf3325f87839dca3452370d16abe704fa1c958 (4 tests).
