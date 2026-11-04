# PREREG r1533: tied-read loci + graph threading, a prototype on the panel (2026-11-04)

From r1532: a family graph resolves tied reads only when every placement of those reads is a path (closed world), and the
catalog has no copy at 9 of the 10 reachable members. The closed world that needs no catalog and no truth: the TIED READS'
OWN PLACEMENTS. Every record of a tied read (AS >= 0.98 x best, a1.molecules.tsv) is a candidate locus; the paths are those
loci; only tied reads enter (uniquely placed reads never do, which is what kept FKBP1C out in r1532).

## Prototype `bench/famcmp/thread_tied.py` (Python; a Rust port only after a PASS)
1. Reads: tied reads with any record at the panel regions (r1532's 28 truth families' member spans: used only to CHOOSE the
   reads for the fast panel, never as loci). All their records come from one `samtools view -N` pass over the BAM.
2. Loci: records clustered by (chrom, strand) and span overlap; a locus's path sequence = the genomic union of its records'
   aligned blocks, strand-normalised (reads align contiguously; constraint 3 of r1532).
3. Graph: loci connected by shared reads form one seqwish graph (minimap2 asm20 all-vs-all, k 19); GraphAligner -x vg
   threads the tied reads. A read is ASSIGNED to the single locus whose private nodes its best alignment covers; a read
   covering private nodes of two loci, or none, stays unassigned (the method is silent where the graph cannot decide).
4. Model: a locus with >= 2 assigned reads gets one model = merged aligned blocks of its assigned reads' records AT THAT LOCUS
   (no borrow strength: only its own reads). Output: a GTF of new models; the base GTF is untouched (additive layer).
Scoring: family_panel.py base = fusedgw/R1/x.gtf vs arm = base + new models, full A119b BAM + a1.molecules.tsv.

## Bars (all must hold; fixed before running)
- reach: a same-strand model at >= 6 of the 10 strict-reached members of r1532 (PRAMEF7, FRG2, FRG2B, OR4F17, OR4F5, TRIM64B,
  LGALS9B, MBD3L3, NPY4R, OR4F21)
- specificity: 0 new models at the 18 no-private-node members of r1532 (exon-identical siblings)
- panel: fused families <= base (49, definition T); complete >= base (209) [additive, so >= is the floor]; T-only hit rate >= 0.678
- new models: >= 0.95 overlap a CAT gene on the same strand (any biotype; models at unannotated loci count against)
- wall: <= 60 min for the whole prototype on the panel's reads
PASS -> a genome-wide row (all tied reads) with the same bars plus chains kept = 1.0 and the one-chain / r1521 bars; a Rust
port is only then on the table. FAIL on any bar -> the route is recorded at its measured size (10 of 1,136 at best) and closed.
