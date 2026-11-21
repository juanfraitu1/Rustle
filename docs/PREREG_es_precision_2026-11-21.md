# Row 1851: why is precision low on the mouse ES library? (pre-registered 2026-11-21)

User (2026-11-21): explore why we miss on precision on mouse ES cells. Today's `final.gtf` there: 17,074 correct of 41,679
reported (0.410); mouse heart, adrenal and cortex on the same annotation 0.52-0.62; IsoQuant on ES 13,443 at 0.472. Measurement
only, GTF level plus the genome for splice motifs; nothing is built. Workspace `r1851_es_precision` (`why.py`). Committed
before the script was written or run. Already known and not repeated: Rustle cut to IsoQuant's count is at 0.611 (row 1850).

## Data
`r1849_dup_count_rerun/{es,heart,adrenal,cortex}/dc/x.final.gtf` (today's default), mouse contigs, GENCODE vM27 as scored
(`genomes/MMU/MMU.ref.gtf`); StringTie and IsoQuant de novo on the ES library (`apetools/MMU_ES_mm`). The three tissues are the
contrast: same annotation, same assembler, higher precision.

## Measured
1. **Classes of the wrong chains** (first that fits; same strand): PIECE = contiguous part of an annotated chain; SHIFT = every
   unannotated intron lies within 10 bp of an annotated one at both ends; COMBINATION = all introns annotated, not a piece;
   NOVEL SITE = some intron annotated, some not (split: the novel intron shares a donor or acceptor with the annotation, or
   neither); NO ANNOTATED INTRON (split: inside an annotated gene's span on its strand, antisense only, intergenic).
   Per class: chains, share of the wrong, median `reads`, share with every intron GT-AG / GC-AG / AT-AC, and on ES the share
   also reported (exact chain) by StringTie, by IsoQuant, by either.
2. **Precision by `reads`** (2, 3-5, 6-20, 21-100, > 100) and by intron count (1, 2-3, 4-9, >= 10), ES and the three tissues.
3. **Agreement:** of Rustle's wrong ES chains the share another tool reports; the same for IsoQuant's and StringTie's wrong
   chains (reported by one of the other two); precision of each tool on the chains at least one other tool also reports.
4. **Depth:** reported, correct and precision when ES is cut to the chains with at least as many reads as give heart's
   reported count (is ES just deeper into the tail?).

## Reading, fixed now
- A class is **likely real and unannotated** when at least half of its ES chains are reported by another tool AND at least
  95% are canonical; **likely an assembler error** when under 20% are reported by another tool; otherwise undecided.
- "ES is deeper" explains the gap to the tissues when ES cut to heart's reported count reaches heart's precision within 0.03.
- "Annotation, not assembler" is supported when the three tools' wrong sets agree at similar rates (within 10 points) and
  Rustle's precision on shared chains is within 0.05 of IsoQuant's.
Everything else is reported as a table for the user; any rule that follows is its own registered row.

## Expectation
ES cells express many unannotated transcripts (repeat-driven and lncRNA starts), and the library is the deepest mouse one, so
I expect NO ANNOTATED INTRON and NOVEL SITE to be the largest classes with high canonical shares and moderate agreement, and
the low-read bins to carry most wrong chains. I expect the depth cut to close part of the gap to heart, not all of it.

## Limits
- Agreement between tools is evidence of a real molecule in the reads, not of a real transcript: all three read the same
  alignments, so a shared artefact (a mis-spliced read population) also agrees.
- Exact-chain sharing understates agreement for chains that differ only at one site.
- One library; the three tissues come from ENCODE and differ in protocol and depth.

## Result (2026-11-21): the gap to the tissues is depth; most "wrong" ES chains are canonical and seen by another tool

Job 4491224 (`r1851_es_precision/why.txt`, 82 s). Mouse contigs, today's `final.gtf`: ES 17,073 correct of 41,641 (0.410).

**1. Depth explains the gap to the tissues (registered reading: within 0.03).**
| tissue | reported, precision | ES cut to that reported count |
|---|---|---|
| heart | 21,160, 0.611 | reads >= 10: 21,153 reported, 13,220 correct, 0.625 |
| adrenal | 23,018, 0.615 | reads >= 9: 22,048 reported, 13,475 correct, 0.611 |
| cortex | 29,767, 0.520 | reads >= 5: 27,952 reported, 14,927 correct, 0.534 |
ES is not a harder library: at the same number of reported chains it is as precise as each tissue. It reports nearly twice as
many chains as heart because it is deeper, and the extra ones are the low-count tail.
Precision by `reads`, ES / heart: 2 reads 0.133 / 0.313; 3-5 0.189 / 0.453; 6-20 0.312 / 0.682; 21-100 0.560 / 0.869; > 100
0.889 / 0.934. A chain with a given read count is worth less in the deeper library (correct chains: median 58 reads on ES,
14 on heart).

**2. Classes of the 24,568 wrong ES chains.** All eight classes are 100% canonical (GT-AG / GC-AG / AT-AC).
| class | chains | share of wrong | median reads | in StringTie | in IsoQuant | in either | reading |
|---|---|---|---|---|---|---|---|
| novel site, shares a splice site with the annotation | 10,153 | 41.3% | 4 | 0.634 | 0.397 | 0.716 | likely real, unannotated |
| combination of annotated introns | 7,219 | 29.4% | 9 | 0.711 | 0.457 | 0.769 | likely real, unannotated |
| novel site, neither end annotated | 1,607 | 6.5% | 4 | 0.619 | 0.445 | 0.714 | likely real, unannotated |
| piece of an annotated chain | 1,565 | 6.4% | 6 | 0.548 | 0.558 | 0.716 | likely real, unannotated |
| shifted site (<= 10 bp) | 1,365 | 5.6% | 7 | 0.246 | 0.145 | 0.284 | undecided |
| no annotated intron, inside a gene | 1,105 | 4.5% | 3 | 0.483 | 0.526 | 0.666 | likely real, unannotated |
| no annotated intron, intergenic | 1,042 | 4.2% | 3 | 0.705 | 0.616 | 0.809 | likely real, unannotated |
| no annotated intron, antisense | 512 | 2.1% | 3 | 0.496 | 0.492 | 0.654 | likely real, unannotated |
Seven of eight classes meet "likely real and unannotated" (at least half reported by another tool, at least 95% canonical).
The shifted-site class is the one the other tools mostly do not report (28%). The class mix is the same on heart, adrenal and
cortex (novel site sharing a splice site 31-36%, combination 35-40%), so ES has no error class of its own.
Expectation: wrong on the composition (I expected chains with no annotated intron to lead; they are 11%); 71% of the wrong
chains are alternative uses of annotated genes: a new site next to a known one, or known introns in a new order.

**3. Agreement between tools (exact chains).**
| tool | reported | precision | wrong chains also in another tool | precision on chains another tool also has | precision on chains only it has |
|---|---|---|---|---|---|
| Rustle | 41,641 | 0.410 | 17,398 of 24,568 (0.708) | 15,952 / 33,350 = 0.478 | 1,121 / 8,291 = 0.135 |
| StringTie | 44,591 | 0.368 | 0.609 | 0.476 | 813 / 11,837 = 0.069 |
| IsoQuant | 28,508 | 0.472 | 0.828 | 0.510 | 468 / 3,052 = 0.153 |
Registered reading "annotation, not assembler": half met. Precision on shared chains is within 0.05 (Rustle 0.478, IsoQuant
0.510, StringTie 0.476): where the tools agree they hit the same annotation ceiling of about a half. The wrong-set agreement
rates are NOT within 10 points (0.61 / 0.71 / 0.83). The difference between the tools is the size of the tool-only tail:
Rustle reports 8,291 chains no other tool has, at 0.135; IsoQuant 3,052 at 0.153.

**Reading (mine).** Three things make up 0.410, in this order: (a) the annotation: even chains two tools agree on match
GENCODE only half the time; (b) depth: ES reports twice heart's chains and the tail is poor; (c) Rustle's own tail: 8,291
Rustle-only chains hold 7,170 of its wrong ones and 1,121 correct ones (6.6% of the correct). Nothing here is an ES-specific
assembler fault. Levers this points to, none built, each its own registered row: a read floor that scales with library depth
(the fixed 2-read floor means something else at 8.7 M reads than at heart's depth), and scoring against an expression-aware
or agreement-based yardstick on real libraries instead of the whole annotation.
Limits as registered: agreement shows a molecule in the reads, not a real transcript.
