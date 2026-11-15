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

## OUTCOME (2026-11-15, job 4472291, `r1751_skipfix/{or,kb,a,hg}.{scan,prec,score}.txt`): FAIL on reach

Register row 1753: another session took 1751 while this ran. The reach denominators are distinct exons (26 / 27 / 34;
the 30 / 36 / 79 above counted chain x exon units). Rewrite check: differing records = logged on all four samples.

| | OR6737 | KB3781 | A119b | HG002 (unseen) |
|---|---|---|---|---|
| decoy fire rate (gate <= 0.005) | 0.0008 | 0.0003 | 0.0011 | 0.0006 |
| row 1746 exons / candidates / fired | 26 / 9 / 9 | 27 / 14 / 14 | 34 / 9 / 9 | 27 / 12 / 12 |
| reach (bar 0.75) | 0.35 | 0.52 | 0.26 | 0.44 |

The gate is valid; the reach bar fails on every sample. Every exon that is a candidate fires. The others are not
candidates: no primary read has both flanking junctions, so the exon is dropped in every read and there is no path in
the locus's junctions to realign to.

Post hoc (step 1 was to be read only on a pass; recorded, decides nothing):

| | OR6737 | KB3781 | A119b | HG002 |
|---|---|---|---|---|
| fired units (exon form only / both / skip only / neither annotated) | 16 / 15 / 21 / 34 | 21 / 15 / 16 / 72 | 28 / 25 / 20 / 77 | 22 / 19 / 25 / 104 |
| reads rewritten | 994 | 8,568 | 2,499 | 8,452 |
| S1 skip chains at exon-only units, base -> fix | 9 -> 0 | 8 -> 0 | 11 -> 4 | 14 -> 2 |
| S2 new correct / lost correct chains | 0 / 7 | 1 / 3 | 0 / 3 | 0 / 2 |
| S3 panel precision, base -> fix | 0.3323 -> 0.3310 | 0.3449 -> 0.3455 | 0.1510 -> 0.1510 | 0.2343 -> 0.2344 |

Where it reaches, the rule removes the false skip chains, but it loses more RefSeq-correct chains than it adds on all
four samples (why was not examined). No stage built. Reaching the remaining exons needs a search of the intron for the
bases the read carries at the junction; not tried. The panel BAMs were deleted; `*.fix.bam.fixes.tsv` keeps every
rewritten read with its old and new CIGAR.
