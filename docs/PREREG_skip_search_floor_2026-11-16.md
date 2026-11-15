# PREREG row 1757 (2026-11-16): short-exon intron search with a floor of 2 reads, on four unseen samples

Row 1754 failed on precision (HG002 0.908 vs 0.95); 8 of its 12 disagreeing finds rested on one read. This row adds
one thing: a found exon needs >= 2 reads. Everything else is row 1754's rule (`bench/famcmp/skip_search.py`, 8th
argument 2) with row 1751's junction rule (`skip_fix.py`). User request 2026-11-16. Written and committed before the
rule was run on the test samples.

Seen, decides nothing: under the floor the four earlier samples give 17 / 18 (OR6737), 15 / 17 (KB3781), 99 / 99
(A119b), 63 / 65 (HG002), and the floor leaves out about 60% of the found units.

## Test samples (none has been run through rows 1746, 1751 or 1754)

Chimp PTR, bonobo PPA, orangutans PPY and PAB: testis reads `reads/TAG/TAG_mm.bam` on their own T2T genomes
(`genomes/TAG/TAG.fa`), RefSeq as `ptr_ppa_rules/{PTR,PPA}.refseq.gtf` and `pongo_asm/{PPY,PAB}.refseq.gtf`. The job
first makes each sample's list of dropped annotated exons (`skip_realign.py`, row 1746) and row 1751's units.

## Step 0

- **Gate:** reversed-intron decoy reads / found reads <= 0.05 on each sample (reads, before the floor), else INVALID.
- **Precision:** found units (>= 2 reads) whose outer junctions flank an annotated exon of <= 15 bp: the annotated
  exon in >= 0.95, on every sample with >= 20 such units; at least 2 samples must have 20, else NOT TESTABLE.
- **Reach:** dropped annotated exons (>= half the skip reads carry the exon) reached by the search or by row 1751's
  rule: >= 0.75 on >= 3 of the 4 samples (samples with < 10 such exons do not count).

## Step 1 (read only if step 0 passes; panel and arms as rows 1751 / 1754)

- **S1:** emitted chains with the skip intron at fired units where only the exon form is annotated fall by >= 50%.
- **S2:** lost RefSeq-correct chains that are not the skip form of a unit with both forms annotated <= new correct
  chains (`skip_fix.py score ... ref=`).
- **S3:** panel precision of `fix` >= `base` - 0.0005.

PASS = step 0 and S1-S3 on >= 3 of 4 samples. A PASS makes the two rules together a candidate for an opt-in stage;
building it is the user's decision. No default changes in this row.

## Limits

RefSeq is thinner on the orangutans (rows 1748 / 1750), so fewer units can be judged there and S2 / S3 undercount
correct new chains. All four are testis; no fibroblast sample. Exons of 3-5 bp stay out of reach.
