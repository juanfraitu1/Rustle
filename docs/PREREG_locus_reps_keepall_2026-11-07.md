# PREREG row 1577: multiple representatives, every attached group kept (no merge-back)

## Why
Row 1575 traced the 25 absorbed genes: 17 get a node, 10 reach a family, 4 survive the attach. The merge-back (built for
row 1571's positional halves of one gene) keeps only the most-read group when two groups of one locus share a family
without aligning: 6 genes lost (NPIPB7 to its overlapping lncRNA group, TCP11X1/X2, FAM25G, OR4M2, ZNF324B). Exon-disjoint
groups with their own reads are separate transcription units, so F2 keeps them all. The other 8 are not attach losses:
they align (UGT2B17 only to intronic repeats; ZNF585B to ZNF585A at 0.91 over 3.4 kb) but fail the edge rule = the
"aligned but not linked" class (row 1575's self-hit guess was wrong; verified 2026-11-07).
Development evidence on A119b AND HG002 (HG002 was row 1576's held-out; it is now seen). A third sample is the held-out.

## Arm F2
Row 1575's F runs as they are (A119b `/anvme/workspace/v138da13-rustle/r1571/F`, HG002 `.../r1576/F`), no piece_merge,
`piece_attach.py --keep-all` onto each sample's arm A.

## Bars per sample (as row 1571)
1. Absorbed pairs still missed <= half of A's.  2. Precision >= A - 0.01.  3. Recall >= A.  4. Copies <= 1.05 x A.

## Prediction
A119b: bar 1 FAIL (at most 10 of 25 genes), bar 2 FAIL (as F); HG002: bar 1 FAIL or borderline, bars 2-4 PASS.
