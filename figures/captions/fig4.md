# Figure 4: What happens to the simulated reads the aligner cannot place

**Status.** The figure is laid out for all six samples, each simulated from its own genome-wide copy catalog
(`docs/PREREG_genome_wide_copy_assignment_2026-09-25.md`, experiment A). Until those runs finish, the rendered figure
shows the development tables: human A119b with its chr16 catalog and gorilla OR6737 with its chr20 (NC_073244.2)
catalog. The other four samples are marked "not run yet". Every number below is from the development tables. The
genome-wide numbers replace them, and the pre-registration fixes, before they exist, the bars that decide each claim.

**Claim (development tables).** Scored within the read's source family, which is known only in simulation, the
copy-assignment test made no wrong call where it assigned. In human A119b with the chr16 catalog it assigned 163 of
the 1,263 reads whose primary alignment has MAPQ 0, all 163 to the source copy. They come from 21 source copies: 153
through a decisive site, and 10 from one copy whose family had a single candidate copy among the read's alignments.
Another 127 MAPQ-0 reads have a decisive site in their source family and were left unassigned. Most MAPQ-0 reads (858
of 1,263) come from copies identical to a sibling over the aligned segment, and 1 of those 858 is assigned.

That reading uses the simulation to pick the source family's result, and the outputs a user can apply do not
reproduce it. The default output, read as any assigned result, is correct for 60 of the 648 MAPQ-0 reads it assigns.
The union test assigns none. 1,255 of the 1,263 MAPQ-0 reads have an NM-identical twin: an alignment at another locus
with the best alignment score and the same number of mismatches and indels as the alignment at the source copy, so no
base of the read separates the two places. Of the other 8, the source copy holds the best score alone for 3, the
source-copy alignment scores below the best for 3, and 2 have no alignment on the source copy. Assignment of these
MAPQ-0 reads is left undecided, not solved. In gorilla OR6737 (chr20 catalog), only 30 reads have MAPQ 0, from three
copies, and none is assigned, so that sample tests neither correctness nor the fraction assigned.

## Caption

**Figure 4 | Fate of the simulated reads that the aligner cannot place (MAPQ 0), one bar per sample.** Reads were
simulated from the spliced sequence of every catalog copy of at least 300 bp in a family of at least two copies. The
read name records the source copy, which is the truth for every number here. Reads were mapped to the whole genome
with minimap2 2.31 and assigned by Rustle's copy-assignment step (`copy_assign --families`). Samples are grouped by
species and never pooled: human A119b and human testis (T2T-CHM13 v2.0), chimpanzee (mPanTro3), gorilla OR6737 testis
and gorilla KB3781 fibroblast (mGorGor1, GCF_029281585.2), and orangutan (mPonPyg2).

Each bar splits one sample's MAPQ-0 reads into five fates. They are scored within the read's source family, which is
known only in simulation (asterisk):
- **Assigned to the source copy:** the copy-assignment test assigned the read to its source copy, or to a catalog
  copy at the same locus (overlap of at least 50% of the shorter span).
- **Assigned to another copy:** a wrong call. Its count is printed in bold for every sample, even when it is 0.
- **At least 1 decisive site, left unassigned:** the read covers a decisive site, but the test did not reach
  significance against every other candidate copy.
- **No decisive site, left unassigned:** no position of the read distinguishes the candidate copies.
- **No result in the source family:** the default output has no result for the read in its source family.

A **decisive site** is a paralogous sequence variant (PSV, a position at which the copies of a family carry
different bases) or a splice junction, covered by the read, at which the candidate copies differ. The
**copy-assignment test** assigns a read to its most likely copy only if, against every other candidate copy,
sequencing error alone would explain its agreement with that copy at the distinguishing sites with probability below
0.001/(n − 1), where n is the number of candidate copies. Otherwise the read is left unassigned.

Columns to the right of the bars: simulated reads; MAPQ-0 reads and their share; wrong / assigned in the
source-family reading; the **default output** (one result per read and family, each family tested on its own, so a
read can be assigned in more than one family) read as a user would, correct / assigned; the **union test**
(`--union-certificate`, one test per read over all its candidate alignments, across families and outside the
catalog), reads assigned; and the MAPQ-0 reads with an NM-identical twin.

**Development tables shown.** Human A119b, chr16 catalog: 28,453 simulated reads, 1,263 at MAPQ 0 (4.4%). Of these,
163 were assigned to the source copy, 0 to another copy, 127 have a decisive site and were left unassigned, 957 have
no decisive site, and 16 have no result. Gorilla OR6737, chr20 (NC_073244.2) catalog: 11,448 reads, 30 at MAPQ 0
(0.3%), of which 10 have a decisive site and 20 do not; none is assigned. The chr16 catalog is the development
substrate of the copy-assignment rules. The chr20 catalog was not used to develop them; it was scored once, as the
held-out substrate, on 2026-09-23.

**Figure 4s (supplement) | Every simulated read of each sample, by set membership.** One UpSet plot per sample,
species in rows. Each read belongs to exactly one column (its exact combination of sets). Bars count reads and are
stacked by the identity of the source copy to its most similar directly aligned copy in the family: 100%\*, 99.5 to
100%, 99 to 99.5%, 98 to 99% and below 98%. **100%\* means identity over the aligned segment, which covers at least
half of the shorter copy; it is not full-length identity.** This is why some reads from that band have MAPQ > 0.
The column of MAPQ > 0 reads is drawn on a broken axis: both parts share one scale. Columns beyond the 12 largest are
pooled as "other". The sets are:
- **MAPQ > 0:** the primary alignment has mapping quality above 0. The copy-assignment step leaves these reads to
  the aligner. MAPQ > 0 does not mean the read has no other placement.
- **MAPQ 0:** the aligner could not choose a placement.
- **Not tested (no result for the read):** no result for the read in the default output.
- **No decisive site over all placements:** the union run's result for the source family has no decisive site.
- **≥ 1 decisive site in the source family:** the default run's result for the source family has a decisive site.
- **Assigned in the source family\*:** that result is "assigned" and the read passed the origin check.
- **… to the source copy.**

## Methods

- **Catalogs.** Genome-wide: each sample's catalog from the run cache
  (`python3 figures/make.py runs --sample <id> --stage catalog`, `gw_family_catalog --piecewise`). Development: human
  chr16 (`chr16_arm/on.copies.*`) and gorilla chr20 (`hom_c234.copies.*`). The catalog handed to `copy_assign` drops
  each copy that was neither simulated nor covered by a simulated primary alignment (`_o2.derive_catalog`: 18 human,
  0 gorilla), because `--families` refuses a copy with no reads.
- **Simulation.** `python3 bench/sim.py copies <catalog.tsv> <catalog.fa> <splice.mmi> sim 20260925`: each copy gets
  min(100, max(10, n_reads)) reads with 0.1% substitutions, 0.03% deletions and 0.03% insertions, up to 10% trimmed
  from one end, then 0 to 30 bp from each end. Reads are mapped with the shipped command
  (`minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`) in read-disjoint parts (each part a separate
  minimap2 run over the same index; minimap2 maps each read independently, so the merged records are those of one
  run). Genome-wide, the number of parts is fixed from the catalog before any read is simulated (17,000 reads per
  part for human, 4,400 for the apes), the first call only simulates (`--simulate-only`), later calls reuse the
  FASTQ (`--reuse-fastq`), and the closest-sibling all-vs-all is skipped (`--sibling none`; no figure uses it). The
  build checks that `sim.bam` holds one primary or unmapped record per simulated read.
- **Assignment.** `copy_assign --bam sim.bam --fasta <genome> --regions <every contig> --families cat.copies.tsv
  --copies-fa cat.copies.fa`, then the same with `--union-certificate`, current defaults, every `RUSTLE_*` variable
  unset. Genome-wide, the runs are split into shards of whole read-connected families (`_o2.plan_shards`): two
  families are in one shard when their read windows (copy ± 50 kb) overlap, when one read has alignments in windows
  of both, or when a family with copies on several contigs meets a single-contig family there. Check V3 (2026-09-25,
  human chr16 simulation, one frozen binary): two shards gave the same rows as one run for the default output
  (3,767 rows), the union run (3,767 rows) and the union side file (993 rows), and `score.py reads` printed identical
  tables. On chr16, 287 of the 290 families fall in one component, so a genome-wide component too large for one
  10-minute call is run as a whole on the cluster (pre-registered stop rule).
- **Scoring.** `python3 bench/score.py reads --catalog cat.copies.tsv sim o2|u2 --per-read <tsv>`.
  `figures/_o2.py` joins the primary alignments (truth, aligner placement) with those per-read verdicts. The build
  stops unless the table's MAPQ-0 tallies equal the scorer's own totals (OWN on the default run, ANY on both runs).
- **Commands.** `python3 figures/make.py data fig4` (heavy; genome-wide: bounded, one heavy step per call, run under
  `flock /mnt/linuxdisk/tmp/rustle_heavy.lock` and repeat until it finishes; `--set o2_scope=dev` rebuilds the
  development tables). `python3 figures/make.py plot fig4` writes `figures/out/fig4_assignability.{pdf,png,svg}` and
  the supplement `figures/out/fig4s_assignability_upset.{pdf,png,svg}`.
- **Table.** `figures/data/fig4_assignability_upset.tsv`: one row per sample, exact set combination, default-output
  and union verdicts, twin state, aligner outcome and identity band. Its notes give the copy counts, the source copies
  of the assigned reads and the twin tally. The fates of the main figure are sums of its rows (`_o2.fate_of`).

## Caveats

- Reads from one copy (10 to 100 per copy) are not independent. The 163 assigned human reads come from 21 copies, so
  Fig. 5 draws its intervals resampling source copies.
- The fates use the simulation to choose the source family's result. They measure the test itself, not what a user
  of the default output sees; the right-hand columns and Fig. 5 give the outputs a user can apply.
- The twin measure also uses the truth (which alignment lies on the source copy). It explains the union test's
  abstention; a user cannot compute it.
- The simulation has no readthrough molecules, no intron retention and no reads from copies missing from the catalog,
  so the correctness shown is an upper bound for real reads.
- Human testis was mapped without `-uf`; its catalog carries that caveat (the simulated reads are mapped with the
  shipped command).
- The development tables use one contig's catalog per sample, mapped to the whole genome; alignments to copies on
  other contigs count as outside the catalog. The genome-wide catalogs remove that restriction.
- Supplementary Fig. 5s compares every simulated read's fate, not only the MAPQ-0 reads', with the alignment-score
  margin rule (a read is assigned to its best alignment only if no other alignment scores within T units).
