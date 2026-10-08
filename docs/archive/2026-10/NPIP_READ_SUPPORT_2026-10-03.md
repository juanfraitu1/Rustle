# Do the reads at an NPIP copy support its transcript? Junction-level check of the A119b primaries at the 25 CAT/Liftoff copies (2026-10-03)

Prompted by the NPIP Read Pools page (artifact F3gJty4egn598SCZ9RBiM1): at NPIPB4 the de novo loci are 2-exon pieces whose "exons" span
annotated introns, so the reads behind them are unspliced or barely spliced. Data: `A119b.t2t.bam` (CHM13 v2.0), same-strand primary
records inside each copy span; scripts inline in the session, outputs `/mnt/linuxdisk/tmp/psv_ceiling/npip*_*.out`.

## NPIPB4 (901 primaries, median read 2,270 bp, 21 CAT exons)
- Unspliced 424 (47 %): 212 inside a single annotated exon (3' fragments), 140 mixed exon/intron, 72 intronic. One junction 219 (24 %):
  only 12 of them match a CAT intron. >= 2 junctions 258; **matching >= 2 CAT introns 144 (16 %)**, >= 2 RefSeq introns 161 (18 %).
- The two annotations disagree with each other: CAT and RefSeq NPIPB4 share **2 of 20 introns**; RefSeq's single model has introns of
  1, 1, 2, 2 bp (frameshift gaps), CAT has introns of 13, 19, 38, 61 bp (lift artefacts). Neither is a usable "transcript to support".
- Annotation-free: a junction seen in >= 3 reads of the copy is "supported". NPIPB4 reads with >= 2 supported junctions: **24 %**
  (median 6 reads per supported junction). The recurrent junctions (22366126-22369656: 60 reads; 22337012-22337742: 49; 22336478-22336852: 39)
  are real and partly CAT's, but no read carries the whole chain: reads are 5' or 3' fragments of a 46 kb gene.
- The primaries-only "own node" of NPIPB4 on the Read Pools page (DN_chr16_22349614_1) is a single 6.6 kb exon over three annotated exons
  and their introns: an unspliced stub, the T5 trap ("one giant exon") in the flesh. The copy keeps a node, but not a transcript.

## All 25 copies (11,471 primaries)
| | unspliced | >= 2 supported junctions | MAPQ 0 overall | MAPQ 0 among unspliced | MAPQ 0 among >= 2-junction reads |
|---|---|---|---|---|---|
| all copies | 16 % | **68 %** | 13 % | **27 %** | **8 %** |
Per copy (unspliced / >= 2 supported): NPIPB2 5/86, NPIPB14P 4/92, NPIPA9 6/88, NPIPB6 2/87 ... **NPIPB4 47/24, NPIPB13 15/24,
NPIPB10P 36/40, LOC124907807 30/45, LOC128966608 29/54**. Tied copies: NPIPA7 63 % MAPQ 0 (its >= 2-junction reads 63 % tied), NPIPA8 82/82,
LOC124907808 63/77: there the spliced reads tie too, because the sibling copies share the exon structure (A6-A9 group).

## Reading
1. The worry is right for NPIPB4 and a handful of copies (B13, B10P, the two LOC12490780x, LOC128966608): most of their reads are
   unspliced or single-junction fragments and do not support any transcript model; a locus built from them is a stub.
2. It is not the norm: at 19 of 25 copies, 65-92 % of primaries carry >= 2 read-supported junctions.
3. Junction count is a specificity signal: unspliced reads are tied 3.4x more often than >= 2-junction reads (27 % vs 8 %). But at copies
   whose siblings share the structure (A7/A8, LOC124907808) spliced reads tie just as much: structure separates copies only where
   structure differs; otherwise only PSVs can.
4. "Supported" must be defined from the reads (junction shared by >= k reads), not from the annotation: at NPIPB4 both annotations are
   wrong and disagree with each other.
Suggested use (not implemented): report per copy the fraction of reads with >= 2 read-supported junctions as the "structure-specific" read
pool; reads outside it count toward a copy only through PSV evidence (O2) and never establish a transcript (O1 node = stub). The
assembler's floor 2 + strict junctions already enforce most of this for transcript models; the gap is that unspliced stubs still become
nodes (register T5, row 1225).
