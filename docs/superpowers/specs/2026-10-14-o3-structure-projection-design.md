# O3: exon structure for reference-absent candidates (structure projection) — design

Date 2026-10-14. Status: implemented in `candidate_copies` (default products `structure.tsv`, `contigs.gtf`); held-out test PASSED (S1 0.946, S2, S3 0.995) — see
`docs/PREREG_o3_structure_projection_2026-10-14.md`.

## 1. Goal: what "O3 works" has to show

O3 = detection from RNA of an expressed copy the reference lacks (Amendment 48). It is proven when three claims carry measured
numbers, each against truth at three levels of realism:

| claim | measure |
|---|---|
| **detection** | the copy is flagged (recall by divergence, expression, structure) |
| **reconstruction** | the candidate's sequence is the copy's transcript (identity x coverage; ≥ 0.999 = exact) — IsoCon's job, IsoCon as the comparator |
| **structure** | the candidate carries the copy's exon structure (exon boundaries vs truth, ±5 bp and exact) — new |

plus specificity (nothing flagged when nothing is missing). Levels: (1) famsim worlds, exact truth (laptop); (2) real reads with an
erased, annotated copy (bed H, main machine); (3) a fresh individual with haplotype assemblies (HG002, cluster + main machine).

This document covers **structure**, the missing piece: `o3_export.py` and `candidate_copies` emit every candidate as a single-exon
sequence. Read correction before clustering (IsoCon's step 2, Amendment 44's diagnosis of the 14/40 exact reconstructions) is the
next step and gets its own design.

## 2. Idea

A candidate's reads are partly unmapped and carry no splice alignment, but the candidate always has a family paralog in the
reference. Its union (the exon union of its clusters, `contigs.fa`) is aligned to the genome with the spliced preset the stage
already uses; minimap2 places the paralog's introns into the alignment. The paralog's exon structure is thereby projected onto the
candidate, and wherever the candidate differs (an exon the paralog lacks, a moved exon) the alignment shows a gap that the rule turns
into a novel segment rather than forcing the paralog's model.

## 3. The rule (frozen 2026-10-14; reference implementation `bench/o3_structure/proto.py`, md5 511562b804907e1f3cb54b008000dc68)

Per flagged candidate, from its genome records (`minimap2 -c --eqx -x splice:hq -uf -N 20`, the union against the genome index):

1. **Frame.** The primary record with the best identity x coverage defines contig and strand; a `-` primary is handled on the reverse
   complement of the union (boundaries mapped back). Records kept: primary/supplementary (`tp:A:P`), same contig, same strand.
2. **Home copy.** The family copy (copies table) overlapping the best record most.
3. **Solid blocks.** Per record, matched runs (`=`/`X` merged) whose longest exact-match run is ≥ **ANCHOR = 8** bp. Shorter runs
   between gaps are chance agreement inside an insertion.
4. **Exons and junctions.** Consecutive solid blocks with an intron (`N`) between them are separated by a junction; without one the
   gap is an indel inside an exon (no boundary).
5. **End exons.** Dropped (their stretch becomes terminal overhang) while (a) they span < **SEG = 50** union bases, or (b) they lie
   outside the home copy ± **SLACK = 2,000** bp and the hop's far side is not within **EDGE_TOL = 5** bp of an annotated exon edge of a
   family copy (a hop onto another copy at a real exon junction is kept; a hop inside an exon is not).
6. **Junction normalisation.** Query gap [qL, qR) vs genome gap [tL, tR) between the two solid blocks: ungapped extension
   (match +1, mismatch −3) of the query gap's head against the genome gap's head (k1), then of the remaining tail against the
   genome gap's tail (k2). Residual [qL+k1, qR−k2): **≥ MIN_SEG = 20** bp → a novel segment, two boundaries; otherwise one boundary at
   qL+k1. This undoes the aligner's shift of an intron into a flanking exon and collapses a shredded novel exon.
7. **Segment pass.** Novel segments ≥ SEG, uncovered stretches ≥ SEG between records, and terminal overhangs ≥ SEG are re-aligned on
   their own (same preset, same contig, same strand, terminal ones only on the home copy ± SLACK): their introns and the breakpoints
   between their records are boundaries; a terminal overhang that lands on the home copy makes its inner edge a boundary.
8. **Breakpoints.** Between two kept records: the edges of an uncovered stretch ≥ MIN_SEG, else the next record's start.

Output: internal boundaries in union coordinates → exons of the candidate in its own coordinates.

## 4. Development record (seeds 11, 12, 13; nothing else used)

Worlds: `bench/o3_structure/world.py` (A, B in the reference; E absent): divergence ladder 0.5-8%, exon loss / gain / duplication /
shuffle, 5' truncation, resized intron, E at 8-16 reads; seed 13 with mixed strands (`--strands`). Read variants E0 (jitter 30) and E2
(30% of reads 5'-truncated). Iterations, each prompted by a dev failure: secondary records excluded (tp:A:S is a secondary, not a
supplementary); insertion-beside-intron normalisation (gained exon swallowing 38 bp of the previous exon); terminal stubs (a 39-bp end
spliced onto the next copy 31 kb away); solid-block junctions (a gained exon shredded into 1-6-bp matches around an intron); the home-copy
hop rule with the annotated-edge exception; exact-run solidity; frame from the primary, not the best secondary; terminal overhangs re-aligned.
Final dev scores (all boundaries within ±5 bp / exact): s11 E0 25/25 (24), E2 25/25 (24); s12 E0 24/24 (17), E2 24/24 (17);
s13 E0 24/25 (22), E2 24/25 (22) — 146/148 within ±5. Not flagged at all: E at 0.5% (linked as an allele by design, delta = 0.958%), one
E at 1% (seed 12).

## 5. Implementation (Rust, `candidate_copies`)

- `PafHit` gains the CIGAR (`cg:Z:`) and `tp`; the genome call already has `-c --eqx`.
- A pure function `project_structure(union, records, genome_contig_seq, copies_of_family, realign)` in `rustle::candidates::detect`
  implementing §3, unit-tested on hand-built records (every rule, both strands) and parity-tested against the prototype on the dev
  worlds' flagged candidates.
- One extra batched minimap2 call for all segments of all flagged candidates.
- New products: `<out>.structure.tsv` (candidate, strand frame, home copy, boundaries, exon blocks, per-junction kind: intron /
  novel_segment / breakpoint, genomic intron coordinates where known) and `<out>.contigs.gtf` (each candidate a multi-exon transcript on
  its own sequence, as `o3_export.py`'s unplaced GTF). Existing products unchanged byte for byte (identity check).

## 6. Not in scope here

Read correction / IsoCon-style clustering (next design); bed H structure scoring (main machine, after this passes); HG002.
