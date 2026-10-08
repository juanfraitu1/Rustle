# Figure 8: Loci in the Liftoff framework

**Status.** Pre-registered in `docs/archive/2026-09/PREREG_liftoff_loci_2026-09-25.md` before any comparison existed. The Liftoff
baseline is being computed genome-wide, one call per block of records (human 36 calls, gorilla 28, chimpanzee 27,
orangutan 27; about 5 minutes each); the only number filled in below is the check against a single run (amendment
3). `python3 figures/fig_loci.py summary` prints every number this caption quotes, with
the pass/fail of claims L1, G1, D1 and F1. Until a table exists, its panel says "not built"; until a sample's default
families (with their copy table) or missing-copy flags exist, its row says so. Amendment 5 (2026-09-25, before any
comparison): claims C3 and F1 apply to Rustle's default de novo families (the families stage) in place of the legacy
copy catalog, which is scored by the same code as a labelled secondary row of `fig8_samples` and not drawn.

**Claim (to be filled from the tables; the bars were fixed before they existed).**
- L1: Liftoff's self-lift re-places at least 99.0% of each genome's gene and pseudogene records at their own
  position, so it is a usable baseline.
- G1: Rustle's guided candidate search finds at least 90% of Liftoff's extra copies (exon union of at least 200 bp,
  outside every annotated record).
- D1: in every sample, the de novo loci cover Liftoff's read-supported extra copies about as often as its
  read-supported annotated loci (at most 0.10 lower).
- F1: in every sample, at least 90% of Liftoff's (record, extra copy) pairs whose two loci are both covered by members
  of Rustle's default families lie in one family.

## Caption

**Figure 8 | Rustle's loci compared with Liftoff's, using Liftoff's own matching criterion.** Liftoff v1.6.3
(Shumate and Salzberg 2021) lifts an annotation from one genome to another. Here each genome's own RefSeq annotation is
lifted onto the same genome with `-copies`, which re-places every gene and pseudogene record and then searches the
whole genome for extra copies of each record. This is the annotation-guided locus baseline: it uses the annotation
and the genome, and no reads. Four genomes, never pooled: human T2T-CHM13 v2.0 (RefSeq full annotation), gorilla
mGorGor1 (GCF_029281585.2), chimpanzee mPanTro3 (GCF_028858775.2) and orangutan mPonPyg2 (GCF_028885625.2). The six
RNA samples are scored against their own genome's Liftoff loci: human A119b and human testis; gorilla OR6737 (testis)
and gorilla KB3781 (fibroblast cell line); chimpanzee PTR; orangutan PPY.

**Liftoff settings.** `-copies -sc 0.95`, with Liftoff's defaults `-a 0.5` (a record is placed only if at least 50%
of its exon bases align), `-s 0.5` (and its exon columns are at least 50% identical) and `-overlap 0.1`; `-f
pseudogene` adds pseudogene records to Liftoff's default `gene` records, so the records are exactly the loci of
Rustle's guided mode. An extra copy is an alignment of the record's whole gene body (introns included) that covers it
end to end, whose exons are at least 95% identical to the record's (`-sc`; Liftoff's default 1.00 admits only exact
copies), and that overlaps no other placed record on the same strand, so an extra copy is never annotated. The search
covers the whole genome, other chromosomes included, and keeps at most 51 end-to-end alignments per record. Liftoff
aligned with its own minimap2 (2.24). To fit the machine's 10-minute limit per step, Liftoff ran once per contig's
records, each time against the whole genome; overlaps between calls were resolved with Liftoff's own rule (an extra
copy loses to an annotated record; of two overlapping extra copies, the more identical one is kept). This is an
approximation of one Liftoff run. A check on chr20–22 (3,750 placed records), with blocks of at most 600 records,
agreed with one run on 3,730 of 3,750 placements and on 438 of 449 extra copies; every difference lay in the arrays of
identical copies on the short arms of chr21 and chr22, where equivalent positions can be exchanged. The pre-registered
bar (identical placements and at least 98% of extra copies) was missed narrowly, so the figure says "approximation"
and a single run per genome is requested on a larger machine.

**Matching.** Liftoff reports, for each placed locus, `coverage` (the fraction of the record's exon bases aligned)
and `sequence_ID` (the fraction of identical columns). Every comparison here uses that criterion: a reference locus is
**found** by a Rustle locus when the Rustle locus's exons cover at least 50% of the reference locus's exon bases. In
one genome the covered bases are the same bases, so Liftoff's sequence_ID of such a match equals its coverage and the
two criteria are one number. A single Rustle locus must reach 50%; strand is not required.

**a** Liftoff's baseline, per genome. Row labels: annotated records re-placed at their own position (same contig,
spans overlapping by at least half). Bars: extra copies with an exon union of at least 200 bp, split by identity to the
record (sequence_ID at least 0.95, 0.98, 0.99 and 1.00; the darker segments are subsets of the lighter ones); the
number in parentheses counts extra copies shorter than 200 bp (mostly small RNAs), which Iso-Seq reads cannot show and
which are never scored for Rustle.

**b** Like for like: both methods start from the annotation and the genome and use no reads. Rustle's guided mode adds
to the annotated records the candidate loci found by aligning each record's longest transcript (minimap2 `-x splice`)
and its coding-sequence envelope (`-x asm20`) to the genome, keeping hits of at least 80% identity over at least half
the query that overlap no annotated record. Here every record is a seed. Top bar of each genome: Liftoff's extra
copies (at least 200 bp, outside every annotated record), found or not by a Rustle candidate. Bottom bar: Rustle's
candidates of at least 95% identity (the `-sc` of panel a), at a Liftoff extra copy or not. Numbers: k of n (95%
interval).

**c** Rustle de novo loci (the genome-wide assembly of each sample: reads and genome, no annotation), scored against
Liftoff's loci as a reference, not as a competitor. Only reference loci with read support are scored: at least 2 reads
whose primary alignment has an aligned block on the locus's exons. Filled circle: annotated records re-placed in place;
open diamond: extra copies (sequence_ID at least 0.95); grey square: the fraction of all de novo loci that lie at a
Liftoff locus (the others are outside every Liftoff locus, which is not an error). Row labels give n (annotated · extra
copies).

**d** The same for the member loci of Rustle's default de novo families (each sample's genome-wide `families` stage:
reads → loci seeded with secondary alignments within 2% of the best score → one representative per locus, its
most-read transcript → families of loci whose genomic sequences align at ≥ 70% identity and share ≥ 60% of the smaller
locus's exonic sequence, Markov clustering). A member locus is represented by its representative's exons (the copy
table copy assignment uses), so a reference locus counts as found only when a family member's representative covers
it; loci outside every family are not in this set (panel c has them).

**e** For each Liftoff (record, extra copy) pair whose two loci are both covered by family members: the fraction in
one default family. The vertical line is the pre-registered bar (0.90). Row labels give k of n pairs. (Figures 6s and
7 use the same pairs as a family reference, restricted to pairs whose two loci both have read support, with every
such pair in the denominator.)

**f** Missing-copy flags. For every annotated record with at least 10 reads, Rustle tests whether the reads hide a
copy the reference lacks. The markers give the fraction of records that have a Liftoff extra copy among all scanned
records (open circle), records whose test fired (triangle) and reference-absent candidates (open diamond). Liftoff's
extra copies are in the reference, so a candidate whose record has one is listed for review when the flag's own search
for another home of its consensus did not land on that copy (row labels).

Intervals: 95%, from 2,000 bootstrap resamples of source records (a record and its extra copies are one unit; random
seed 20260925). Exposure: nothing here was scored before. The family rules were developed on human chr16 and gorilla
chr20 (NC_073244.2), and the guided search's thresholds on three human families (NPIP, TBC1D3, AMY); panels c-e are
also reported without those contigs (table `fig8_samples`, scope `minus_dev`). The legacy copy catalog (`catalog`
rows of `fig8_samples`) is a comparison with no claim.
