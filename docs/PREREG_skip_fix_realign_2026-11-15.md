# PREREG row 1751 (2026-11-15): de novo local realignment of reads that skip an exon of <= 15 bp

Row 1746: reads that skip an annotated exon of <= 15 bp carry its bases 77-85% of the time (16 bp and up: under 2%), so
minimap2 drops those exons. That row used the annotation to know where the exons are. This row tests a rule that does not:
the exon must be spelled by another read of the locus (a path of the locus's splice graph). Written and committed before
the rule was run on any sample; the per-unit tables of row 1746 (OR6737, A119b, KB3781) were already seen, HG002 was not.

## Rule (`bench/famcmp/skip_fix.py`, nothing tuned)

- Candidate: an intron (d, a) of a primary read, where another primary read has the consecutive introns (d, x), (y, a)
  with 0 < y - x <= 15.
- Test, per read with the intron (d, a): row 1746's test (25 bp window, exon form better by >= 2 edits), the window
  clamped to the read's own flanking exons. Several candidate exons: the best fit.
- Fix: the read's CIGAR is rewritten through the exon (start, end and sequence unchanged). Primary alignments only.

## Step 0: read level (`skip_fix.py scan`), four samples

- **Gate (decoys):** up to 300 plain introns per contig (>= 5 reads, not candidates) get a made-up exon of a real
  candidate's length, cut from the middle of the intron; up to 20 reads each. Reads fired / reads tested <= 0.005 on
  each sample, else INVALID. This is the rate at which stray bases at a junction pass for an exon.
- **Reach:** of row 1746's units where >= half the skip reads carry the exon (30 OR6737, 79 A119b, 36 KB3781; HG002's
  list is made in this job by `skip_realign.py`), the share where the rule fixes >= 1 read. Bar: >= 0.75 on HG002 and
  on >= 2 of the 3 seen samples.
- Reported, no bar: fired units by annotation (only the exon form, both, only the skip, neither). "Neither" is the
  count of unannotated short exons lost the same way, not measured so far.

## Step 1: assembly (read only if step 0 passes the gate and the reach bar)

Panel BAM = all reads overlapping a fired unit +- 100 kb; arm `base` = that BAM, arm `fix` = the same with the CIGARs
rewritten. Both assembled in this job with the default command (same binary), scored on the panel.

- **S1:** at fired units where the annotation has only the exon form, emitted chains with the skip intron fall by >= 50%.
- **S2:** new RefSeq-correct chains >= lost RefSeq-correct chains (`chain_prec.py`).
- **S3 (no harm):** RefSeq chain precision of `fix` >= `base` - 0.0005.

PASS = S1-S3 on HG002 and on >= 2 of the 3 seen samples. A PASS makes the rule a candidate for an opt-in stage, which
is the user's decision; no default changes in this row. Check of the rewrite: records differing between the two BAMs =
reads logged as rewritten.

## Limits

Only an exact (d, a) skip intron is a candidate: a read whose junction the aligner also shifted is not reached. An exon
no read spells is not reached. With a 2-edit margin an exon of 2-5 bp fires on few differing bases; the decoy gate
measures that. Panel scores are not genome scores (loci cut at the panel edge are cut the same way in both arms). The
class is small (row 1746: 30-80 annotated units per sample), so whole-GTF figures will barely move.
