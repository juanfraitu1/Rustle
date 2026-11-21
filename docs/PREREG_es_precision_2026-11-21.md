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
