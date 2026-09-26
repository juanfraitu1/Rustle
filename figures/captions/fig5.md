# Figure 5: Copy assignment by copy identity (simulation), and the hard-locus benchmark (real reads)

**Status.** Panels a and b are laid out for all six samples, and panel c for the two samples on which the lab ran
StringTie 3.0.1, FLAIR 3.0.1 and IsoSeq collapse (human A119b, gorilla OR6737 testis), each genome-wide
(`docs/PREREG_genome_wide_copy_assignment_2026-09-25.md`, experiments A and B). Until those runs finish, the rendered
figure shows the development tables: in a and b, human A119b with its chr16 catalog and gorilla OR6737 with its chr20
(NC_073244.2) catalog; in c, only the chr16 NPIP benchmark (the inset). Every number below is from the development
tables. The pre-registration fixes, before the genome-wide numbers exist, the bars that decide each claim.

**Claim (development tables).** The simulation is scored against each read's source copy.

1. **Scored within the read's source family, the copy-assignment test made no wrong call where it assigned.** In
   human A119b (chr16 catalog) it assigned 163 of the 1,263 MAPQ-0 reads, all to the source copy. They come from 21
   source copies (95% interval resampling copies, 0.85 to 1): 153 through a decisive site, 10 from one copy whose
   family had a single candidate copy. It assigns most MAPQ-0 reads from copies 99 to 99.5% and 98 to 99% identical
   to their most similar copy (fraction assigned 0.77 and 0.80), where nearly every such read has a decisive site
   (56 of 57, 51 of 51). On those bands the aligner's primary alignment is on the source copy for 0.66 and 0.53 of
   the reads it places. Where no decisive site exists, the test leaves the read unassigned: of the 858 MAPQ-0 reads
   from copies identical over the aligned segment, 4 have a decisive site and 1 is assigned.
2. **That reading needs the simulation to pick the source family, and the outputs a user can apply do not reproduce
   it.** The default output tests each family on its own. Read as any assigned result, it is correct on 60 of its
   648 assigned MAPQ-0 reads (0.09), below the aligner's 0.50. The union test assigns none of the 1,263. For 1,255
   of them no test could: they have an NM-identical twin, an alignment at another locus with the best alignment
   score and the same number of mismatches and indels as the source-copy alignment (Fig. 4). Assignment of these
   reads is left undecided, not solved. In gorilla OR6737 (chr20 catalog) no reading assigns any of the 30 MAPQ-0
   reads, so that sample tests neither correctness nor the fraction assigned.
3. **At the real NPIP locus the fraction carried depends on the stratum.** The hard set is defined by Rustle's own
   copy-assignment step, and Rustle's arm here is that step's transcripts (`copy_assign --families --gtf`), not the
   assembler of Figs 1–3. On all hard molecules Rustle's transcripts carry 0.858, IsoSeq collapse 0.732, FLAIR 0.473
   and StringTie 0.441. On the contested stratum IsoSeq collapse carries 0.813 and Rustle 0.692, with 3,403
   transcripts overlapping the copies against Rustle's 911.

## Caption

**Figure 5 | Copy assignment of the reads the aligner cannot place, by copy identity (a, b), and the hard-locus
benchmark on real reads (c).**

**a, b,** One small multiple per sample, species kept together and never pooled (a: human A119b, human testis,
chimpanzee; b: gorilla OR6737, gorilla KB3781, orangutan); axes are shared. The rows are the simulated reads of
Fig. 4 whose primary alignment has MAPQ 0, first all of them, then one row per identity band of the source copy: its
identity to its most similar directly aligned copy of the family, where 100%\* is identity over the aligned segment
(at least half of the shorter copy), not full-length identity. The numbers at the right give the reads and the
distinct source copies of each row. Left axes: the **fraction assigned**. The pale bar behind the points is the share
of the row's reads with at least one decisive site in their source family (a PSV or splice junction, covered by the
read, at which the candidate copies differ), which is what the test can use. Right axes: the **fraction correct
among assigned**, with a 95% interval **resampling source copies**, because reads from one copy are not independent
(2,000 resamples, seed 20260925). When every assigned copy has the same ratio, the resampling is degenerate and the
interval is the Wilson interval at the pooled ratio with the copies as the sample size. With fewer than 5 assigned
source copies no interval is drawn, and the note under the panel says so. Four readings, each a filled marker:
- **Copy-assignment test, scored in the source family\* (circle):** the test's own correctness; the simulation picks
  the source family's result, which a user cannot do.
- **Default output, any assigned result (plus sign):** what a user of the default output sees. A read claimed for
  two loci counts as assigned and not correct.
- **Union test (triangle):** one test per read over all its candidate alignments, across families and outside the
  catalog, so a read is assigned at most once.
- **Aligner's primary alignment (grey cross):** the catalog copy that the primary alignment overlaps most. A primary
  outside every catalog copy counts as not placed. This is the pre-registered baseline.

The **copy-assignment test** assigns a read to its most likely copy only if, against every other candidate copy,
sequencing error alone would explain its agreement with that copy at the distinguishing sites with probability below
0.001/(n − 1), where n is the number of candidate copies; otherwise the read is left unassigned.

**Development tables shown.** Human A119b, chr16 catalog (the development substrate of the copy-assignment rules):
1,263 MAPQ-0 reads from 99 copies. The test assigns 163 (fraction assigned 0.13), all correct, from 21 copies
(interval 0.85 to 1). By band it assigns 1 of 858 at 100%\* (51 copies; 4 reads with a decisive site; 1 copy
assigned, no interval), 56 of 225 at 99.5 to 100% (0.25; 139 with a decisive site; 7 of 26 copies, interval 0.65 to
1), 44 of 57 at 99 to 99.5% (0.77; 56 with a decisive site; 5 of 7 copies, interval 0.57 to 1), 41 of 51 at 98 to
99% (0.80; 51 with a decisive site; 6 of 7 copies, interval 0.61 to 1), and 21 of 72 below 98% (0.29; 30 with a
decisive site; 2 of 8 copies, no interval). The default output assigns 648 MAPQ-0 reads (0.51), mostly through other
families' results; 60 are correct (0.09; interval 0.03 to 0.18), 432 are wrong and 156 are claimed for two loci. The
union test assigns none: 1,247 are left unassigned and 16 have no result. The aligner places 1,184 on a catalog copy,
596 of them on the source copy (0.50; interval 0.47 to 0.54), and 79 outside every catalog copy. Gorilla OR6737,
chr20 (NC_073244.2) catalog, not used to develop the copy-assignment rules (scored once as held out, 2026-09-23): 30
MAPQ-0 reads from three copies (20 at 100%\* from two copies, 10 at 99.5 to 100% from one). The 10 have a decisive
site, and all 30 have an NM-identical twin. No Rustle reading assigns any of them. The aligner places 13 on a catalog
copy, 9 of them on the source copy (0.69; three copies, no interval), and 17 outside every catalog copy.

**c,** Real reads carry no copy truth. What is scored is whether a method's transcripts carry a molecule's exact
intron chain inside a catalog copy (the derived copy call of `score.py bakeoff-calls`: the copy is the one the
transcript overlaps most; for a single-exon molecule the exact chain means span containment). **Rustle's arm is the
copy-assignment step's own transcripts (`copy_assign --families --gtf`)**, not the assembler of Figs 1–3. The other
methods are the lab's StringTie 3.0.1, FLAIR 3.0.1 and IsoSeq collapse on the same BAMs. **The hard set is defined by
Rustle's copy-assignment step:** the molecules its gate admits (at least two alignments within 50 kb of a catalog copy
on one contig, the runner-up alignment score equal to the best) whose primary alignment lies inside a copy, a
universe that is the same for every method. Three groups: all hard molecules; contested (at least two candidate
copies, and the best candidate passes the origin check); and molecules whose exact chain is carried by at least two
reads. Genome-wide panels (human A119b, gorilla OR6737): every family with at least two copies in the sample's
genome-wide catalog; each point is one family with at least 20 hard molecules (the fraction of its hard molecules a
method carries), and the black bar is the fraction over all hard molecules of the genome. **Inset:** the chr16 NPIP
family (26 copies, human A119b), the development benchmark:

| Group | n | Rustle | StringTie | FLAIR | IsoSeq collapse |
|---|---:|---:|---:|---:|---:|
| All hard molecules | 4,115 | 0.858 | 0.441 | 0.473 | 0.732 |
| Contested | 1,046 | 0.692 | 0.327 | 0.359 | 0.813 |
| Chain carried by ≥ 2 reads | 3,688 | 0.957 | 0.488 | 0.507 | 0.737 |

The chain-carried-by-≥ 2-reads group was added after the NPIP pre-registration's prediction P5 failed (post hoc
there); it is pre-registered for the genome-wide run. Transcripts overlapping the NPIP copies (in the table, not
drawn): Rustle 911, StringTie 486, FLAIR 917 and IsoSeq collapse 3,403; the shares whose exact intron chain (ends not
compared) is carried by at least one read are 0.99, 0.75, 0.99 and 0.84. Rustle's value is close to 1 largely by
construction, because it assembles its chains from reads, so the share flags transcripts no read carries and is not
a precision. Transcripts overlapping two or more copies (conflation): 1, 2, 11 and 57.

## Methods

- **a, b.** The runs of Fig. 4 (`figures/_o2.py`): the simulation (`bench/sim.py copies`, seed 20260925),
  `copy_assign --families` with and without `--union-certificate`, then `python3 bench/score.py reads --per-read`.
  The build requires the table's tallies to equal the scorer's totals. Intervals: `_o2.copy_level_ci`. The
  decisive-site share (`n_family_psv`) and the twin count (`n_nm_twin`) come from the same per-read records.
- **c, NPIP inset.** One current run writes both the transcripts and the rows that define the hard set:
  `copy_assign --bam hsa16.bam --fasta chm13v2.0.fa --region chr16:11963320-80438591 --families copies16.tsv
  --copies-fa copies16.fa --gtf --origin-drop-indels --threads 4`, default settings, no `RUSTLE_*` variables (88 s,
  3.1 GB). `hsa16.bam` is `samtools view -b -M -L copyregions.bed A119b.t2t.bam` (copy spans ± 100 kb). Each method is
  scored with `score.py bakeoff-calls <gtf> hsa16.bam copies16.tsv --fuzz <0|5> --tx-out …`, then `bakeoff-compare
  --tx-support …`, once as is and once with `--min-mult 2 --bam hsa16.bam`. Junction tolerance: 0 bp for Rustle,
  StringTie and FLAIR; 5 bp for IsoSeq collapse, which itself runs `--max-fuzzy-junction 5`.
- **c, genome-wide.** `fig_assign_accuracy.ensure_hard_gw`: the catalog's families with at least two copies, in
  shards of read-connected families (`_o2.plan_shards` over the sample's BAM); per shard, `copy_assign --bam <full
  BAM> --regions <every contig carrying a copy> --families <shard> --copies-fa <shard> --gtf --origin-drop-indels`;
  a subset BAM of the shard's copy spans ± 100 kb (the `hsa16.bam` recipe); the lab GTFs restricted to the
  transcripts that overlap the shard's copies (`bakeoff-calls` ignores the others, so no call changes); copy ids
  renumbered across families (`bakeoff-calls` keys copies by `copy_idx`, which repeats across families); then the
  same `bakeoff-calls` / `bakeoff-compare` commands. Per-family and pooled fractions are tallied per molecule by
  `bakeoff-compare`'s own definitions and checked shard by shard against its printed totals. A molecule's family is
  the family whose copy its primary alignment overlaps (`primary_local = 1`). Validation (2026-09-25): run on the
  NPIP catalog alone, this pipeline reproduced the inset exactly (4,115 hard molecules; 0.858, 0.441, 0.473, 0.732;
  contested 1,046).
- **Commands.** `python3 figures/make.py data fig5` (heavy, bounded: one heavy step per call, run under
  `flock /mnt/linuxdisk/tmp/rustle_heavy.lock` and repeat until it finishes; `--set o2_scope=dev` rebuilds the
  development tables; `--set fig5c_samples=none` skips the genome-wide panel c). `python3 figures/make.py plot fig5`
  writes `figures/out/fig5_assign_accuracy.{pdf,png,svg}`.
- **Tables.** `figures/data/fig5_assign_accuracy_bands.tsv` (sample × reading × band), `fig5_hard_locus.tsv` (every
  group `bakeoff-compare` prints for NPIP; the pooled groups genome-wide), `fig5_hard_locus_transcripts.tsv`, and,
  once the genome-wide run exists, `fig5_hard_locus_families.tsv` (sample × family × method).

## Caveats

- The hard set and its groups are defined by Rustle's copy-assignment step, not by the BAM alone. Panel c therefore
  measures how well each method's transcripts represent the molecules that step flags; it does not measure per-read
  copy correctness (panels a and b do, in simulation).
- Rustle's arm runs with the current `copy_assign` defaults, which postdate the NPIP pre-registration; the table notes
  list the differences. The other methods are fixed recorded outputs.
- The composite quoted in older notes (0.846 / 0.574 / 0.502 / 0.755) came from a retired scorer on 3,633 molecules.
  It is a different metric and is not shown.
- The circle in a and b uses the simulation to choose the source family's result. Only the default output and the
  union test are available to a user.
- The simulation has no readthrough molecules, no intron retention and no reads from copies missing from the catalog,
  so the correctness in a and b is an upper bound for real reads.
- Panel c's Rustle arm is `copy_assign --gtf`, a different product from the assembler of Figs 1–3.

---

# Figure 5s (supplement): Rustle and the alignment-score margin rule

**Status.** Laid out for all six samples (a, b: simulation) and for human A119b and gorilla OR6737 (c: real reads),
genome-wide (`docs/PREREG_genome_wide_copy_assignment_2026-09-25.md`, Amendment 2, experiment C, written before any
number below existed). Until the genome-wide runs finish, a and b show the development simulations of Figs 4 and 5
(human A119b with its chr16 catalog, gorilla OR6737 with its chr20 (NC_073244.2) catalog) and c is empty. The claims
are decided on the genome-wide runs only. The rule is attributed to the Eichler lab as the user states it; the
citation is still to be added.

**Claim (development tables).** The margin rule assigns a read only when its best alignment score beats every other
alignment of the read by at least T. Where the rule assigns, Rustle almost always puts the read in the same place,
and Rustle also places reads the rule discards.
- Human A119b (chr16 catalog, 28,453 simulated reads). At T = 10 the rule assigns 26,669 reads (fraction assigned
  0.937), with fraction correct 0.9993. Rustle puts 26,660 of them in the same place (0.9997).
  - The other 9 are exceptions: the rule is right on 4 and Rustle on 2.
  - Rustle also keeps 522 reads that the rule discards, at the aligner's placement (MAPQ above 0 but a margin
    below 10). 516 of them are correct (0.989; 95% interval resampling source copies, 0.969 to 1). This is just
    under the pre-registered bar of 0.99, and most of these reads come from copies 99.5 to 100% identical to their
    closest copy (258 reads, 253 correct).
  - The union test adds no MAPQ-0 read. Scored within the source family, the test adds 163 reads, all correct.
- Gorilla OR6737 (chr20 catalog, 11,448 reads). The rule assigns 11,386 (0.995) and is correct on 11,218 (0.985).
  - All 168 exceptions are reads the rule sends to an unspliced alignment elsewhere in the genome, where Rustle
    keeps the aligner's spliced primary alignment on the source copy. The truth favours Rustle in all 168.
  - Cause: the spliced alignment's score includes the intron penalties, so an unspliced copy of the transcript
    can outscore it. Example: AS 515 for the spliced primary alignment against AS 545 for the unspliced one. The
    168 reads come from two source copies.
  - Rustle also keeps 32 reads (all correct, from one copy) that the rule discards.
- T is a convention. Human: the rule assigns 0.956 at T = 1 and 0.903 at T = 20 (fraction correct 0.9990 and
  0.9994). Gorilla: 0.997 at T = 1 and 0.992 at T = 20.

## Caption

**Figure 5s | Rustle and the alignment-score margin rule, on simulated reads (a, b) and real reads (c).**
- **The margin rule** assigns a read to its best-scoring alignment only if no other alignment of the read, anywhere
  in the genome, scores within T alignment-score (AS) units. A read with no other alignment is assigned. Otherwise
  the read is discarded. The bars use T = 10. Black ticks show where the rule's part of each bar would end at T = 1
  and at T = 20. The rule counts every mapped primary and secondary alignment (supplementary alignments are other
  pieces of the same read, not alternative placements); a missing AS counts as 0.
- **Rustle.**
  - Reads whose primary alignment has MAPQ above 0 stay at that alignment: the copy-assignment step leaves them to
    the aligner (Fig. 4).
  - Reads with MAPQ 0 get the copy-assignment result. The union test (one test per read over all its alignments,
    across families and outside the catalog) is what a user can apply, and it is what the bars show. The test scored
    within the read's source family is known only in simulation (asterisk) and is given as a number in a.
- **Placement.** A placement's copy is the catalog copy its alignment overlaps most, or "outside the catalog".
  - **Correct** (simulation only): the source copy, or a catalog copy at the same locus (overlap of at least 50% of
    the shorter span). Another copy, or a place outside the catalog, is wrong for both rules.
  - **Same place:** the same alignment, or the same catalog copy or locus.

**a,** Every simulated read of each sample (reads from every catalog copy of at least 300 bp in a family of at least
two copies, mapped to the whole genome), split by what the two rules do with it:
- both assign it to the same place (grey);
- the rule assigns it and Rustle places it elsewhere or leaves it unassigned (orange; the columns say which rule the
  truth favours);
- only Rustle assigns it, from the aligner's placement (dark blue) or from the copy-assignment test (light blue),
  with wrong ones in red;
- neither assigns it (light grey).

The columns give the reads the rule assigns and its fraction correct among assigned. Then the share Rustle places
the same, the exceptions, the reads only Rustle assigns with their fraction correct, and the test scored within the
source family\*. Samples are grouped by species and never pooled.

**b,** The same split per identity band of the source copy: its identity to its most similar directly aligned copy
in the family. 100%\* means identity over the aligned segment, which covers at least half of the shorter copy; it is
not full-length identity.

**c,** Real reads carry no copy truth, so this panel shows agreement and coverage only, at T = 1, 10 and 20.
- **Reads compared:** those whose primary alignment overlaps a copy of a family with at least two copies in the
  sample's genome-wide catalog (the reads of Fig. 5c).
- **The rule** uses the best and second-best score over every alignment of the read in the whole BAM.
- **Rustle** uses the primary alignment's copy (MAPQ above 0) or the union test's result (MAPQ 0).

## Methods

- **Simulation (a, b).** The runs of Figs 4 and 5, with no extra alignment or assignment run.
  - `python3 bench/score.py eichler --sim sim --catalog cat.copies.tsv --default o2 --union u2 --per-read <tsv>`
    reads every non-supplementary alignment of `sim.bam` (`-F 2048`, unmapped records counted as unassigned). It
    computes the best and second-best AS per read and joins Rustle's answer: the primary alignment for MAPQ above 0;
    for MAPQ 0, the assigned, not origin-rejected results of the default run (source family) and of the
    `--union-certificate` run (any family).
  - `figures/_o2.py` checks every read against Figs 4 and 5's own records: the same reads, the same aligner
    placement, and the same MAPQ-0 verdicts as `score.py reads --per-read`. The build stops on any difference.
  - Intervals: `_o2.copy_level_ci`, resampling source copies.
  - Commands: `python3 figures/make.py data fig5` (with the rest of Fig. 5), or, from runs that already exist and
    without running anything heavy, `python3 figures/fig_assign_accuracy.py margin-rule [--set o2_scope=dev]`.
    `python3 figures/make.py plot fig5` writes `figures/out/fig5s_margin_rule.{pdf,png,svg}`.
- **Why not `copy_assign --eichler-margin`.** Those columns see only the alignments inside the swept windows of one
  contig. A read whose only rival lies on another contig therefore looks rival-free and is assigned. Under
  `--no-as-tied-only` they also count supplementary records as rivals. The rule is therefore computed from the BAM
  genome-wide (pre-registered, Amendment 2).
- **Real reads (c).**
  - On each shard of Fig. 5c's genome-wide run, one extra `copy_assign --bam <BAM> --fasta <genome> --regions
    <shard regions> --families <shard> --copies-fa <shard> --union-certificate --threads 4` run is made, with current
    defaults and every `RUSTLE_*` variable unset.
  - Then `python3 bench/score.py eichler --real --bam <BAM> --catalog <the shards' families>
    --as-table runs/<sample>/<sample>.molecules.tsv --union <shard prefixes> --cache <dir>`. The genome-wide best and
    second-best AS come from the sample's `as_table` (one pass over the whole BAM). The alignments over the copies
    show where the best one lies; the pass over the BAM is cached per contig and bounded per call.
- **Tables.**
  - `figures/data/fig5s_margin_rule_sim.tsv`: sample × identity band × T × reading × group of reads (the groups
    drawn in a and b).
  - `fig5s_margin_rule_accuracy.tsv`: fraction assigned and fraction correct among assigned for the rule at each T,
    for Rustle, and for the reads only Rustle assigns, with intervals.
  - `fig5s_margin_rule_real.tsv`: panel c, once run.

## Caveats

- The comparison is with the rule as stated. A pipeline that ranks alignments by another score, or that treats
  spliced and unspliced alignments differently, would behave differently on the reads where an unspliced alignment
  scores higher.
- Most of Rustle's extra reads are the aligner's own placements, not the copy-assignment test's. The copy-assignment
  step leaves every read with MAPQ above 0 to the aligner, including reads whose best and second-best scores differ
  by less than T, or are equal.
- The union test assigns almost no MAPQ-0 read in simulation (Fig. 4), so the test's own contribution shows only
  when it is scored within the source family, which needs the simulation.
- Reads from one copy are not independent; intervals resample source copies. The simulation has no readthrough, no
  intron retention and no reads from copies missing from the catalog.
- In the development tables, alignments to copies of families on other contigs count as outside the catalog, for
  both rules.
