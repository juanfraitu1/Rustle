# PREREG r1532: family graph threading of multimapper reads: a ceiling measurement (user, 2026-11-04)

The user asks for family-wide assembly by threading reads through a family variation graph (the thesis model, LEDGER_GENERAL §18:
graph = model, linear alignment = implementation). On the r1521 panel the multimapper misses are 35 of the 212 expressed
no-model members (28 T-only + 7 MAPQ0-share >= 0.5); the other 177 have uniquely placed primaries and are a single-exon /
unspliced-emission matter, not a multimapper one. So the ceiling of ANY multimapper method on this panel is 35 members,
most with 2-5 tied reads in total (FRG2, PRAMEF7, SPANXC): a PASS means 2-read models at those genes.

## Constraints (binding, from the register and the family_gfa design)
1. Paths come from TRUTH loci here (15 of the 28 T-only misses have no catalog copy, r1528 provenance pre-check): this row is a
   ceiling measurement; a tool would first need O1 to find the locus. It must never be read as a method result.
2. Not borrow strength (register L647-648): a copy's model uses only reads threaded to its own path, never pooled family reads.
   Not a family definer: graphs are built from known members only.
3. Paths are each member's EXON-UNION sequence, strand-normalised (as r1529 step 0), so a cDNA read aligns contiguously; a
   skipped exon is a deletion, not an intron gap (GraphAligner has no spliced mode).

## Step 0 (fixed before any graph is built)
- Members: the 35 multimapper misses, in their truth families; paths = exon-union sequence of every expressed (T >= 2) member of
  each family.
- Graph: one seqwish graph per family (minimap2 -cx asm20 -X all-vs-all, seqwish -k 19; ~/runs/family_graphs/build.sh pattern),
  pixi graph env with the default env on PATH.
- Reads: every record at every member span of the family in the A119b BAM, primaries and secondaries, deduped by name, SEQ
  taken from a record that carries it.
- Threading: GraphAligner (-x vg) per family. A read REACHES a member when its best alignment covers >= 1 node private to that
  member's path (on no other path). Members are classed: reached (>= 2 reaching reads) / path but no private node covered /
  no private node in the graph (identical to a sibling at the exon level).
- **Stop rule: fewer than 9 of 35 members reached (25%, as r1529) -> closed, ceiling recorded, no arm.** >= 9 -> the arm is a
  new row (DNA-first loci + threading + per-path assembly), designed then.
- Wall: minutes; interactive.
