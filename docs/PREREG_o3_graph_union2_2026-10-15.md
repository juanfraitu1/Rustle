# Union-path graph O2 with padded path ends, and family-pooled O3, on a fourth handpicked family set (2026-10-15)

Pre-registration, written before any world of the new held-out set was built. Code: branch `o3-graph` (`bench/o3_graph/`),
the commit that adds this file.

**User goal (2026-10-15):** "improve the simulation and use handpicked families (maybe soto or NPIP or other mainstays) to
improve the graph based O2 and O3".

Round 3. Round 2 (`docs/PREREG_o3_graph_union_2026-10-15.md`) claimed union O2 in part and graph O3 without the removed-copy
clause. Each failure pointed at one fix:
- **U4** (0.981 on paralog-structure reads): 5 AMY1 reads ran 16 bp past the end of union O2's own projection of their
  structure. Those bases counted as distance against the true copy, more than the near-identical copies' PSVs.
- **G3** (2 of 3 removed copies): PSG's removed copy was split over four hosts (27, 16, 10, 8 reads). Each host's cluster was
  pure, but none held 80% of the copy's reads.

## Methods under test

**O2, primary: padded union O2** (`o2graph.py run WORLD unionp 2`). Round 2's union O2 (every family isoform structure
threaded through every copy), with each path's two transcript ends extended by 50 bp of the copy's own genome. A read
running past an annotated end is compared with each copy's real sequence there. It is neither penalised nor excused for
the annotation.

A first version, which simply did not count read bases past a path end ("unione"), fixed AMY. But it also dropped real end
evidence: on dev-3 CROSS + REMOVED it fell from x1.07 to x1.03 and made 2 wrong calls. Padding replaced it.

**O3, primary: H** (`detect.py`). As round 2's G (cluster only the reads union O2 cannot explain, the mosaic rule, branch
calls, the genome check), but on padded union O2's unexplained reads POOLED over all of the family's hosts. The pool is
aligned to each host and clustered there, so an absent copy that copy_assign spreads over several hosts clusters once.

Reported beside them: the current O2 and CUR; round 2's union O2 and G; path O2, locus O2; S23, S123m; padded union O2's
gains split into PSV-separated (the runner-up copy can carry the read's structure) and structure-separated (`o2graph.py
gains`; round 2: 35% / 65%).

## Worlds

As round 2, seed 41: INTACT, OUTSIDE, REMOVED (special copy hard-masked, own index), BRANCH (novel 150-bp exon), and CROSS with
`--projection liftover` (paralog structures lifted through the loci's genomic alignment, independent of union O2's own
spliced projection).

## Families (held out, picked before any run)

X-linked cancer-testis families SSX, MAGEA, SPANX; UGT2B; the opsins OPN1LW/MW; CGB; KIR; CYP4F. TRIM, LILR, DEFB10 and CSAG
were named as spares, in order, should one form no family; none was needed.

| family | copies / isoforms |
|---|---|
| SSX | 5 / 10 |
| MAGEA | 5 / 16 |
| SPANX | 4 / 4 |
| UGT2B | 3 / 8 |
| OPN1 | 3 / 9 (OPN1LW, OPN1MW, OPN1MW2) |
| CGB | 6 / 10 |
| KIR | 8 / 16 |
| CYP4F | 2 / 8 |

`/mnt/linuxdisk/tmp/o3_graph/families3.json` (its earlier sets are identical to `families2.json`).

## Known before running (development: dev, dev-2 = register 1298's held-out set, dev-3 = round 2's held-out set)

O2, padded union O2 (margin 2) against the current O2 and round 2's union O2; correct / assigned, wrong-family calls counted:

| set | worlds | current | **padded union O2** | union O2 (round 2) |
|---|---|---|---|---|
| dev | INTACT, REMOVED, OUTSIDE, BRANCH, CROSS (spliced), CROSS (liftover) | 10,939 / 10,957 | **11,315 / 11,315 (x1.034)** | 12,391 / 12,391 (x1.133) |
| dev-2 | INTACT, REMOVED, OUTSIDE, BRANCH, CROSS (spliced) | 11,530 / 11,657 | **12,701 / 12,701 (x1.102)** | 12,701 / 12,701 (x1.102) |
| dev-3 | INTACT, CROSS (liftover), OUTSIDE, REMOVED, BRANCH | 6,678 / 6,678 | **7,106 / 7,106 (x1.064)** | 7,170 / 7,175 (x1.074) |

- Paralog-structure reads, padded union O2: dev 883 / 883, dev-2 611 / 611, dev-3 261 / 261 (round 2's AMY1 errors gone).
- Padding costs most on dev (x1.133 -> x1.034). There copies' annotated isoforms end at different points, and the reads,
  trimmed inside their own isoform, ran past a paralog's shorter path end. That counted as a difference, but it is the
  annotation again, not sequence. Padding lets each copy's genome decide those bases, so the honest gain is family-set
  dependent (x1.03-x1.10) -- hence V2's 2% below (5% and 10% reported).
- Gain split, padded union O2 on dev-3: separated by PSVs 102, by structure 342 (the runner-up cannot carry the read's
  structure under the projection rules).

O3 (families with a detection or a false call; detectable = median `de` > delta):

| set | | CUR | S23 | G | **H** |
|---|---|---|---|---|---|
| dev | BRANCH detected | 0 / 8 | 8 / 8 | 8 / 8 | **8 / 8** |
| dev | REMOVED detected / detectable | 3 / 4 | 3 / 4 | 4 / 4 | **4 / 4** |
| dev | false-call families, INTACT + OUTSIDE + both CROSS | 4 / 32 | 0 | 0 | **0** |
| dev-2 | BRANCH | 0 / 8 | 8 / 8 | 8 / 8 | **8 / 8** |
| dev-2 | REMOVED | 0 / 2 | 1 / 2 | 2 / 2 | **2 / 2** |
| dev-2 | false-call families of 40 | 13 | 0 | 0 | **0** |
| dev-3 | BRANCH | 0 / 8 | 8 / 8 | 8 / 8 | **8 / 8** |
| dev-3 | REMOVED | 2 / 3 | 2 / 3 | 2 / 3 | **3 / 3** (PSG, pooled) |
| dev-3 | false-call families of 40 | 3 | 0 | 0 | **0** |

## Bars (held out: the 8 families above, seed 41, five worlds; wrong-family calls count as wrong)

| id | bar |
|---|---|
| V1 | padded union O2 precision >= 0.995 in every world (reads of catalogued copies) |
| V2 | padded union O2 correct >= 1.02 x the current O2's, pooled over the five worlds |
| V3 | padded union O2 correct >= the current O2's in every world |
| V4 | padded union O2 precision >= 0.99 on CROSS's paralog-structure reads |
| V5 | padded union O2 assigns no more reads of non-catalogued copies than the current O2, pooled |
| H1 | H detects >= 6 of the 8 BRANCH absent copies |
| H2 | H has a false call in <= 1 of the 24 INTACT, OUTSIDE and CROSS families |
| H3 | H detects every detectable REMOVED copy CUR detects, and >= 0.75 of the detectable ones; needs >= 3 detectable |
| H4 | H detects at least as many absent copies (REMOVED + BRANCH) as S23, with no more false-call families |

Detectable (H3): the removed copy's reads have a median primary `de` on the masked genome > delta. Detection and false
calls use `score.py`'s rule (register 1294).

## Readings (fixed now)

1. **V1-V5 pass:** claimed, "on simulated reads of handpicked real families, O2 on a family graph that threads every isoform
   structure through every copy, each path ending in its copy's own genome, assigns more reads correctly than the current
   O2 (>= 2%, the figure reported) at >= 99.5% precision, including when isoforms are annotated on another paralog".
2. **H1, H2, H4 pass and H3 passes:** claimed, "graph O3 that clusters the family's O2-unexplained reads together detects
   novel-exon absent copies and >= 75% of detectable removed copies, with at most one false call over intact,
   uncatalogued-paralog and paralog-annotation families". H3 not measurable: claimed without the removed-copy clause.
3. Whether the O2 step adds detections: H detects more than S23 (H4 holding) -> "the O2-informed pooling adds detections on
   held-out families"; equal -> "it adds nothing on this set" (reported, not a failure).
4. Each claim stands or falls separately. A failed bar is reported by name.

Scope either way: simulated reads (curated CAT isoforms, HiFi errors) on the real CHM13 genome. Real data's other species,
fragments and allele background are not modelled.

## Results

Added after the runs. Nothing above this heading changes once a statistic has been seen.
