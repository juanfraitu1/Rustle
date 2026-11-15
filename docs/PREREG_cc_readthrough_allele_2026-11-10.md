# PREREG r1623 / r1624: step 0 for the two parked correlation-clustering uses (2026-11-10)

User (2026-11-10): "can we try one or more of the other possible uses for correlation clustering?" The parked ones are tandem
readthrough cuts and allele vs paralog. r1620 and r1621 both showed that a gain needs repel evidence that exists and points
the right way, so each use gets a counting step first and a clusterer only if it survives. Written before any count.

## r1623 tandem readthrough: does the repel evidence exist, and only at true boundaries?
In data terms: inside one assembled locus, two reads ATTRACT when their aligned spans overlap and REPEL when they do not. A
clusterer would cut where reads confined to each side outweigh the reads that cross.
- Cases: the 69 spliced fused loci of `~/fused/tandem_fused.tsv` (A119b; the 4 `mono_blob` rows are O2's), first two genes,
  RefSeq gene spans. Reads: primary alignments of `A119b.t2t.bam` on the genes' strand. Left gene = smaller start.
  confined-left = read overlaps the left gene and ends at or before the right gene's start; confined-right the mirror;
  crosser = overlaps both genes and is neither. Gene spans that overlap are counted as "no boundary".
- **Cut statistic (frozen):** both confined counts >= 2 and crossers < the smaller confined count.
- Control: RefSeq multi-exon genes on no fused row with >= 20 same-strand primary alignments, the first 1,000 by md5 of the
  gene name. A control gene is a FALSE CUT if the statistic holds at ANY of its intron positions (confined-left = read ends
  at or before the intron start, confined-right = starts at or after its end), since a clusterer picks its own boundary.
- **Stop rule:** STOP if cut loci <= 10 of 69 (rule T reached 6) or control false cuts > 5%. Otherwise register the clusterer.

## r1624 allele vs paralog: closed at step 0 from the record
In data terms: reads piled on one copy ATTRACT when they agree at the copy's variable columns and REPEL when they differ; two
clusters appear both for two haplotypes of one copy and for a collapsed second copy.
- What would tell them apart must be written first, and the only quantity available from RNA reads alone is polymorphism of
  the copy's own confidently placed reads at those columns. That is the r1541 allele-site guard (FAIL: fixed 4 wrong,
  abstained on 12 right; panel 0 of 9). Linking columns by read does not add a new quantity.
- Size of the class where truth exists (KB3781, haplotype truth): one family, 12 reads (MCL183, the 1413 trap).
- **STOP, no run.** Reopen only with a second signal (parental or DNA haplotypes) or a sample where the class is larger.

## Outcome r1623 (2026-11-10, `bench/famcmp/readthrough_step0.py`, `r1623.step0.tsv` in the workspace): STOP
- Fused spliced loci 69: gene spans overlap at 3; the cut statistic holds at **59 of 69**. The repel evidence exists: most
  fused loci have confined reads on both sides that outnumber the crossers.
- Control, 1,000 single genes with >= 20 reads: the statistic holds at some intron of **266 (0.266; bar 0.05)**. Reads confined
  to the two ends of one gene (partial reads, alternative first and last exons) look the same as two genes. STOP as registered.
- Same wall as r1516/r1517 from the other side: the signal reaches the true boundaries and cannot tell them from the inside
  of a gene. No clusterer written.
