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

## Step 0 RESULT (~/r1532/step0_thread.py -> step0.out; graphs in ~/r1532/fam/<family>/; 65 s)
| class (registered rule) | members |
|---|---|
| reached (>= 2 reads cover a private node) | **12** (stop rule < 9 not fired) |
| path has private nodes, none covered | 5 |
| no private node: identical to a sibling at the exon level | 18 (SPANXC, CDY2A, DEFB103A/B, DEFB105A/B, DEFA1, PPIAL4E/F, H2AB2/3, 5 TAF11L5 copies, TRIM49D2, MBD3L4) |
The 18 are a ceiling for any transcript-level graph: nothing in their spliced sequence separates them from a sibling.

### After the numbers (labelled; the registered rule is not rewritten)
The registered reach rule counts reads whose LINEAR primary is unique at a sibling or elsewhere (check_strict2.out): FKBP1C 1193
such reads (1153 with a unique primary at FKBP1A; the graph has no path for the chr15/chr19 pseudogenes their homology points
to), FRG2B 100 (primaries on chr16/20/22), LGALS9B 47, TRIM64B 24, TRIM49D1 21. A family graph built from the truth family's
expressed members is not a closed world: reads from homologs absent from the graph land on the nearest path and cover its
"private" nodes. Strict count (reaching reads with an on-strand primary at the member, or AS-tied in a1.molecules.tsv):
**10 of 35 reached** (FKBP1C and TRIM49D1 drop), still above 9. Tied reads resolved per member: OR4F21 25, OR4F5 16, LGALS9B
14, OR4F17 13, FRG2B 12, FRG2 9, TRIM64B 8, MBD3L3 2, NPY4R 2, PRAMEF7 2 (the last three on 1-bp private nodes: single-SNV
calls, marginal). So threading is sound only for reads whose linear evidence is already tied (O2's scope) or on a graph that
holds every genomic copy of the family, i.e. O1's catalog: the sketch finds the members, the graph resolves the reads.
Size of the route before any arm: 10 members of 1,136, with 2-25 tied reads each; the best cases are single-exon OR models.
The arm is a new row; it adopts the strict rule.

### Gating check for the arm (~/r1532/gate_catalog.out)
Strict-reached members with a catalog copy at 0.60: **1 of 10** (OR4F17); the siblings the tied reads sit on have a copy in 3
of 10 (OR4F17, LGALS9, OR4F16), never in the same catalog family as the member. These families have no model on any member,
so O1 has no locus to sketch from: the arm's first component is LOCUS DISCOVERY, not threading.
