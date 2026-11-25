# How the comparison assemblers treat multi-mapping reads, and what Rustle does differently (2026-11-25)

Question (user, for the advisor): say with foundations what Rustle does differently from, or like, the other assemblers
to use secondary alignments, above all near-tied multimappers. Method: each tool's SOURCE was read at the lines that
decide (installed version and a clone of the upstream repository, workspace `assembler_src_2026-11-25/`), and each
tool was RUN on the same planted reads with and without the secondary alignments (job 4504616, `check/run.sbatch`).
Nothing here is taken from a paper or a manual alone.

## 1. What the aligner hands over
All tools read one BAM: `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`. Per read: one primary
alignment and up to 50 secondary alignments whose chaining score is at least 0.1 x the best. minimap2 ranks on `ms` and
breaks exact ties by a hash, so at identical copies the primary is at an arbitrary copy (row 1922: wrong unit for 12.6% of
planted reads). It writes no `NH` tag. Its mapping quality is low where copies tie: on the planted set of the test (half the
copies silent) 46% of primary alignments have MAPQ 0 and 66% are below 30.

### 1b. Where those properties come from in minimap2's own source (2.30; `hit.c`, `format.c`, `minimap.h`)
- **Which alignments become secondary:** after alignment the hits of a read are sorted by score, and a hit that overlaps a
  better one on the read is its secondary. It is kept when its score is at least `-p` x the primary's and fewer than `-N`
  are kept already (`mm_select_sub`, `hit.c:265`). Our `-p 0.1 -N 50` therefore keeps almost everything, down to a tenth
  of the best score; minimap2's own defaults (`-p 0.8 -N 5`) would keep few. This is why 98% of our secondary
  alignments are far from ties and every tool is right to ignore THOSE.
- **Which of two equal hits is primary:** the sort key is `score << 32 | hash` (`mm_hit_sort`, `hit.c:203`), where the hash
  is derived from the hit's anchors and a per-read seed (`hit.c:64`). Equal scores are ordered by that hash: deterministic
  for a given read and index, arbitrary with respect to biology. This is the "coin toss".
- **The score that is compared** is `dp_max`, the best-scoring segment of the alignment; it is written as the tag `ms`
  (`dp_max0`), while `AS` is the full dynamic-programming score `dp_score`, which in splice mode also pays for introns
  (`format.c:403`, `minimap.h:105-106`). Hence "say which score": two hits can tie on `ms` and differ on `AS`.
- **Mapping quality** is computed for the primary only from the gap to the best alternative: with `x = dp_max2 * subsc /
  (dp_max * score0)` it is proportional to `(1 - x^2)`, capped by a second formula `6.02 x identity^2 x (dp_max - dp_max2)`,
  then reduced by `4.343 x ln(number of suboptimal hits + 1)` (`mm_set_mapq2`, `hit.c:444-482`). An exact tie gives 0; every
  secondary alignment gets 0 (`hit.c:482`). So a tool that filters on mapping quality (IsoQuant) is filtering on exactly
  this gap, and a tool that ignores it (StringTie, bambu, FLAIR at threshold 0) treats a tied primary like a unique one.

## 2. Tool by tool (version; the deciding lines)

**StringTie 3.0.3, long-read mode (`-L`)** (`stringtie.cpp`, `rlink.cpp`; the modified copy `~/stringtie_modded` is 3.0.1 with
logging added and decides the same way)
- Secondary and supplementary alignments are DROPPED before anything is built: `rlink.cpp:713` upstream
  (`if((longreads || ...) && (flags & BAM_FSECONDARY || flags & BAM_FSUPPLEMENTARY)) return;`; `rlink.cpp:579` in the modified
  copy, which counts them as `dropped_secondary_or_supplementary`).
- They are still READ by the main loop, which decides where a bundle (the stretch of reads processed together) starts and
  ends before that line is reached (`stringtie.cpp:945-1131`, then `processRead` at 1166). So a secondary alignment can
  join or lengthen bundles without adding a read to any of them.
- No use of mapping quality. Read weight = 1 / `NH` (`rlink.cpp:719` in the copy; `-u` switches it off), and a bundle whose
  multi-hit share exceeds `-M` (default 1 = never) is skipped; minimap2 writes no `NH`, so every alignment weighs 1 and the
  filter never fires.
- In SHORT-read mode (no `-L`) secondary alignments are used as reads at full weight (same BAM: 493 transcripts with them,
  191 without). That mode is not the one compared.

**FLAIR 3.0.1** (`flair_align.py`, `flair_collapse.py`, `count_sam_transcripts.py`)
- Genome alignments: primary only. FLAIR's own aligner call sets `--secondary=no` (`flair_align.py:150`), and its filter
  keeps `read.is_mapped and not read.is_secondary` (`flair_align.py:187`); the benchmark wrapper calls that same filter.
  Supplementary alignments are removed by default. Mapping quality threshold 0.
- Then every read is RE-ALIGNED to the first-pass isoform sequences with `minimap2 -N 4` (`flair_collapse.py:326`) and given
  to ONE isoform: candidates are sorted by alignment score, matched bases, clipping, transcript length and then NAME, and the
  first is taken (`count_sam_transcripts.py:241-255`). A read that fits two copies' isoforms equally goes to the one whose
  name sorts first.
- An option exists for exactly this case: `--allow_paralogs`, "allow reads to be assigned to multiple paralogs with
  equivalent alignment" (`count_sam_transcripts.py:56`, passed on at `flair_collapse.py:347`): the read then counts for every
  tied isoform. Off by default; not used in our comparisons so far.
- An isoform needs 3 supporting reads (`--support`, default 3).

**IsoQuant 4.1.0** (`isoquant.py`, `isoquant_lib/alignment/alignment_processor.py`, `isoquant_lib/model_construction/`)
- Secondary alignments are not used unless `--use_secondary` (`isoquant.py:289-291`, default off); supplementary ones never
  (`alignment_processor.py:385-387`).
- Mapping quality is used three times: an alignment with 1 or 2 exons and MAPQ < 1 is dropped
  (`alignment_processor.py:399-401`, `--simple_alignments_mapq_cutoff` 1); a novel model with 1 or 2 exons is deleted when
  the MEAN MAPQ of its reads is below 30 (`model_filter.py:72-75` and `112-116`, cutoff set at `isoquant.py:1306`);
  `--min_mapq` (unset) would drop all alignments below it.
- A read flagged `multimapper` builds nothing: it is left out of the intron collection, the intron graph and the path
  counts (`intron_graph.py:56`, `280`, `629`; `intron_path.py:28`). With `--use_secondary` a read with several consistent
  placements is flagged so (`assignment/multimap_resolver.py`, strategy `take_best`).

**bambu 3.12.1** (`R/prepareDataFromBam.R`)
- Reads the BAM with `scanBamFlag(isSecondaryAlignment = FALSE)` (lines 31-32): secondary alignments are never loaded.
  Supplementary alignments are loaded. No mapping-quality filter. One read, one place; nothing models a tie.

**IsoCon 0.3.3**
- Reference-free: it never sees a genome alignment, so "secondary alignment" does not exist for it. Reads are clustered
  by sequence, and a statistical test decides whether two near-identical clusters are different transcripts. Its test is
  what Rustle's sequence count ports (row 1926).

## 3. The direct test (planted human families, half the copies silent; 4,200 primary + 68,929 secondary alignments)

| tool | transcripts on the full BAM | on primaries only | same structures |
|---|---|---|---|
| StringTie `-L` | 156 | 156 | 156 |
| FLAIR | 203 | 203 | 203 (its read table is byte-identical) |
| IsoQuant | 123 | 123 | 123 |
| bambu (NDR 1) | 198 | 198 | 198 |

Every one of the four gives the same assembly whether the secondary alignments are in the file or not.
Their own switches, on the same reads:
- FLAIR `--allow_paralogs`: 216 transcripts, all 203 of the default plus 13.
- IsoQuant `--use_secondary`: 30 transcripts instead of 123 (12 in common). Turning secondaries on flags most family reads
  as multimappers, and those build nothing. (Registered afterwards on nine fresh sets, row 1947: section 8.)

On a real sample StringTie's output does change slightly with the secondaries in the file (row 1870, chimpanzee testis:
62,418 transcripts with, 62,494 without, 61,268 shared, named families identical). The source explains it: they move
bundle boundaries, they are not evidence.

Context, ONE planted set, scored per copy (the registered comparisons are rows 1917-1918, not this table):

| tool | expressed copies found (of 70) | silent copies called (of 63) | isoform recall | isoform precision |
|---|---|---|---|---|
| StringTie | 70 | 30 | 0.669 | 0.568 |
| FLAIR | 70 | 30 | 0.831 | 0.551 |
| FLAIR `--allow_paralogs` | 70 | 31 | 0.839 | 0.523 |
| IsoQuant | 64 | 24 | 0.581 | 0.615 |
| IsoQuant `--use_secondary` | 19 | 4 | 0.121 | 0.500 |
| bambu | 67 | 26 | 0.718 | 0.464 |
| Rustle, tied secondaries seeding | 67 | 48 | 0.815 | 0.432 |

Note what the middle column shows: tools that never read a secondary alignment still model 24-30 silent copies, because
the PRIMARY of a tied read lands on a silent copy by the aligner's coin toss. Ignoring secondaries does not avoid the
problem; it hides which models rest on tied reads.

## 4. What Rustle does (where; all in this repository)
1. **Admits a secondary alignment as building evidence only when the aligner could not separate it from the best**: AS >=
   0.98 x the read's genome-wide best, from one table per BAM (`src/family.rs`, `as_tie_keep_with`; driver `--seed-pool
   good`). Every other secondary is ignored, as in the other tools. On real BAMs that is 1.8% of secondary alignments.
   A stricter band on `ms` was tested: more precise, fewer copies on the deepest sample, not adopted (rows 1942, 1942b).
2. **Rescues family copies no two reads agree on**, inside families only (stage `multicopy`).
3. **Labels instead of dropping**: an EM over the reads tied between loci, with the alignment scores inside it, marks each
   locus `expressed`, `shared` or `homology_only` (`bench/o2_families/em_filter.py`). Nothing is removed; the label says
   what the reads can support. Per-read shares can be written (`PREFIX.o2.shares.tsv.gz`).
4. **Counts the distinct sequences** expressed among copies the reference cannot separate, with IsoCon's test on the tied
   reads (`bench/o2_families/recurrence_test.py`, opt-in).
5. **Records the ties**: `PREFIX.o2.ties.tsv.gz` lists, per tied read, every locus it can sit at with its score gap; the
   provenance table marks near-tied reads.
It does NOT assign a tied read to one copy, and it does not use mapping quality anywhere.

## 5. Same and different

| | StringTie `-L` | FLAIR | IsoQuant | bambu | IsoCon | Rustle |
|---|---|---|---|---|---|---|
| secondary alignments as evidence | no (dropped) | no | no (opt-in breaks it) | no (not loaded) | no genome | near ties only |
| supplementary alignments | dropped | dropped | dropped | kept | - | not used |
| mapping quality | not used | threshold 0 | drops MAPQ-0 short alignments; needs mean MAPQ 30 for 1-2 exon models | not used | - | not used |
| a read tied between copies | counted at the primary's copy | one isoform, by name on a tie (all with `--allow_paralogs`) | builds nothing when flagged | counted at the primary's copy | by sequence | seeds every tied copy; shared out by an EM for the labels |
| says which models rest on tied reads | no | no | no | no | - | yes (labels, tie table) |
| counts sequences among identical copies | no | no | no | no | yes (its purpose) | yes (IsoCon's test, ported) |

**Like the others:** Rustle also ignores the large majority of secondary alignments (98% on real BAMs), builds on primary
alignments first, and cannot say which copy a tied read came from.
**Unlike the others:** it is the only one of the six that (a) uses the near-tied part of the secondary alignments as
evidence to build, (b) reports per locus how much of its support is tied, and (c) keeps a per-read record of the ties.
IsoQuant is the only other tool that reacts to ambiguity at all, and it does so by excluding: low mapping quality removes
short alignments and short models, which is why it under-calls family copies (64 of 70 here; 0.59 isoform recall in row 1917).
FLAIR has the nearest idea (`--allow_paralogs`), but applied after its first-pass isoforms were built from primary
alignments only, so a copy that received no primary gets no isoform to share reads with.

## 6. What may be said, and what may not
- SAY: "With the same alignment file, StringTie (long-read mode), FLAIR, IsoQuant and bambu produce the same assembly whether
  the secondary alignments are present or not; we checked in their source and by running them."
- SAY: "Rustle uses the secondary alignments the aligner could not separate from the best one; the others use none."
- SAY: "Ignoring secondary alignments does not remove the problem: the primary of a tied read is an arbitrary choice, and
  each tool models 24-30 silent copies from it on planted reads."
- DO NOT SAY "StringTie uses secondary alignments" (row 1870's wording, corrected below) nor "StringTie discards them
  without effect": they change bundle boundaries.
- DO NOT SAY the other tools "ignore multimappers" for IsoQuant: it filters on mapping quality, which is a reaction to them.
- DO NOT SAY Rustle finds more expressed copies in general (rows 1917-1918: only at about 3 reads per copy).

## 7. Correction of an earlier record
Row 1870 concluded "StringTie USES secondaries" from a 2% difference in its output with and without them. The measurement
stands; the reading was wrong. The source drops secondary and supplementary alignments in long-read mode before the
splice graph, and the difference comes from bundle boundaries. Statements in `docs/MULTICOPY_MODE.md` are corrected.

## 8. The items first listed as not done, now closed
- **FLAIR `--allow_paralogs` as a registered arm (row 1946, three fresh seeds):** at 3 reads per copy it lifts FLAIR from
  44 / 40 / 40 to 54 / 48 / 47 expressed units of 108 and adds 5 false units at 60 reads; Rustle's mode finds 104 / 102 / 104,
  and its labels leave 4 / 3 / 4 silent units called against FLAIR's 19 / 15 / 17. FLAIR needs 3 supporting reads, so the
  3-read rows sit at its own floor; the comparison that isolates the option is FLAIR against FLAIR.
- **IsoQuant `--use_secondary` as a registered arm (row 1947, the same nine sets):** it finds about a third of the units it
  finds by default (102 / 99 / 100 -> 35 / 33 / 37 of 108; isoform recall 0.60 -> 0.14). Reading the secondary alignments
  makes IsoQuant flag the family reads as multimappers, and flagged reads build nothing.
- **bambu's supplementary alignments (row 1947, gorilla OR6737):** 0.8% of primary alignments genome-wide, 1.3% at the
  family loci; removing them changes 49 of about 1,870 bambu models there (30 lost, 19 gained). A small difference, and a
  supplementary alignment is a second piece of a split read, not a second placement.
- **The aligner side:** section 1b, from minimap2's source.

So the two tools that HAVE a switch for multi-mapping reads use it in opposite ways: FLAIR counts a tied read at every tied
isoform and gains a little; IsoQuant excludes reads with several placements and loses most of the family. Neither labels
what rests on ties.

## 9. Still open
- The other tools were run at their benchmark settings; a tuned FLAIR (`--support 2` with `--allow_paralogs`) was not tried.
- All direct runs here are on planted human families plus one gorilla sample for bambu; no real-sample arm of FLAIR
  `--allow_paralogs` was run.
