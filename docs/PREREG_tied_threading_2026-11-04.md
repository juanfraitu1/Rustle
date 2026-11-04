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

## RESULT (job 4441399, ~/r1533/run.4441399.out; products /anvme/workspace/v138da13-rustle/thread/)
Panel: 19 families, 905 tied reads, 3,156 tied placements -> 271 loci in 88 connected groups; 581 reads unassigned (no single
locus with private nodes covered), 29 loci under 2 reads, **31 models**; prototype 526 s.
| bar | value | verdict |
|---|---|---|
| reach: >= 6 of the 10 strict-reached members gain a model | **1** (MBD3L3) | FAIL |
| specificity: 0 models at the 18 exon-identical members | **3** (DEFB103B, DEFB105A, TAF11L5) | FAIL |
| fused <= 49 | 49 | PASS |
| complete >= 209 | 209 (T), 210 -> 212 (P) | PASS |
| T-only hit rate >= 0.678 | 0.690 | PASS |
| new models on a CAT gene same strand >= 0.95 | 25 of 31 = 0.806 | FAIL |
| wall <= 3600 s | 526 s | PASS |
**Verdict: FAIL as registered. Route closed at its measured size.**

### After the numbers (checked; ~/r1533/check_spec.out, check_loci.out)
- Assignment split over the 905 reads' best alignments: 342 cover private nodes of exactly one locus, 305 cover two or more
  loci, 250 cover none. The method was confused (305) more than silent (250).
- Locus construction is the defect: overlap-clustering of placements chains across tandem clusters into mega-loci. Model span
  median 13.9 kb, 17 of 31 over 10 kb, 7 over 100 kb (chr7:11839-252847, chr1:10799-175488); 13 of 31 models span >= 2 CAT
  genes. The TAF11L5 (50 kb, 12 kb) and DEFB105A (25 kb) models are such mega-loci; the DEFB103B model (1.4 kb) is a genuine
  assignment to an exon-identical sibling, i.e. the private node there is not real separating sequence.
- Not checked: whether a locus definition that does not merge across genes would raise the reach above 1 of 10; the
  step-0 graphs on truth exon unions (r1532) already bound it at 10, with 2-25 reads each.
