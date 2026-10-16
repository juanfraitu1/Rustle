# Union-path graph O2 and union-O2-informed O3 on a fresh set of handpicked families (2026-10-15)

Pre-registration, written before any world of the new held-out set was built. Code: branch `o3-graph` (`bench/o3_graph/`),
the commit that adds this file.

**User goal (2026-10-15):** "improve the simulation and use handpicked families (maybe soto or NPIP or other mainstays) to
improve the graph based O2 and O3".

This follows `docs/PREREG_o3_graph_real_2026-10-15.md` (register 1298-1299). There, path O2 (each read against every curated
isoform of every copy) gained +14.3% only because each copy's catalog entry listed exactly the isoforms it expressed. When
copies also expressed isoforms annotated on a paralog (world CROSS), it named the annotated paralog: precision 0.42 on those
reads, 0.946 overall, and fewer correct calls than the current O2. On O3, the O2 step added nothing (S4 = S23).

## Methods under test

**O2** (`o2graph.py`; the reads O2 places on a family = copy_assign's rows):

| method | what it does |
|---|---|
| current O2 | copy_assign `--no-as-tied-only`: rows with status `assigned` |
| **union O2** | the family graph threads EVERY family isoform structure through EVERY copy. A copy's paths: its own curated isoforms, plus each paralog isoform spliced-aligned to its locus (kept with >= 95% coverage, the same exon count, canonical junctions, an intron chain new to it). Then as path O2: edits + unaligned bases to the copy's best path; assigned when the next copy is >= MARGIN farther, abstains when tied or thin, `unexplained` beyond delta. **MARGIN = 2.** Only the sequence separates copies, not which copy the annotation lists a structure on. |
| path O2 | each copy's own curated isoforms only (register 1298), reported |
| locus O2 | the annotation-free baseline: the same decision on spliced alignments to each copy's genomic locus (any canonical junction), reported |

**O3** (`detect.py`):

| arm | what it does |
|---|---|
| CUR | `copy_assign --psv-clusters` |
| S23 | branch calls (reads on no path of the family's vg graph) + PSV calls, kept only when the genome does not explain their reads |
| S123m | S123 with the MOSAIC rule: a cluster no single copy explains whose NEAREST copy differs from its consensus at > delta of its bases (GOLGA8's removed copy, register 1298, carries alleles of several copies) |
| **G** | union-O2-informed O3: the clustering run only on the reads union O2 calls `unexplained`, with the mosaic rule, plus the branch calls, kept by the genome check |

## Worlds (`realworld.py`)

As register 1298 (handpicked CAT families on the real CHM13 genome, reads from curated isoforms, anonymised names), seed 31:
INTACT, OUTSIDE, REMOVED (special copy hard-masked, own index), BRANCH (novel 150-bp exon, 0.2% SNPs), and **CROSS with
`--projection liftover`**. In this CROSS, every copy also expresses up to 2 paralog-annotated structures, lifted onto it through
the two loci's genomic alignment (minimap2 asm20). That projection is independent of the spliced transcript alignment union O2
itself uses. On development, 45 of its 66 structures matched the spliced projection exactly, 12 differed in exon ends, and 9
were different choices.

## Families (held out, picked before any run)

Soto's human-specific duplications FAM72, ARHGAP11 and TCAF2; copy-number mainstays AMY, SULT1A, FCGR, APOBEC3 and PSG.
ARHGAP11 and APOBEC3 form no family at the build's 90% identity linkage (built, no world). GTF2H2 (the SMA region) and LRRC37A
(17q21.31) were added in their place, still before any run.

| family | copies / isoforms |
|---|---|
| FAM72 | 4 / 11 |
| TCAF2 | 2 / 6 (TCAF2, TCAF2P1) |
| AMY | 9 / 30 |
| SULT1A | 4 / 16 |
| FCGR | 4 / 16 |
| PSG | 7 / 28 |
| GTF2H2 | 2 / 8 |
| LRRC37A | 3 / 8 |

`/mnt/linuxdisk/tmp/o3_graph/families2.json` (its dev and first held-out sets are identical to `families.json`).

## Known before running (development)

Development = the 8 original development families, plus the 8 families of register 1298's held-out set. Those are seen now
and count as development ("dev-2").

O2, union O2 (margin 2) against the current O2; correct / assigned, wrong-family calls counted:

| set | worlds | current | union O2 | path O2 | locus O2 |
|---|---|---|---|---|---|
| dev | INTACT, REMOVED, OUTSIDE, BRANCH, CROSS (spliced) | 8,703 / 8,712 | **9,944 / 9,944 (x1.143)** | 9,967 / 10,203 | 8,612 / 8,612 |
| dev | CROSS (liftover) | 2,236 / 2,245 | **2,447 / 2,447 (x1.094)** | 2,222 / 2,443 | 2,247 / 2,247 |
| dev-2 | INTACT, REMOVED, OUTSIDE, BRANCH, CROSS (spliced) | 11,530 / 11,657 | **12,701 / 12,701 (x1.102)** | 12,707 / 12,865 | 11,285 / 11,436 |

- Paralog-structure reads: union O2 is right on every one it assigns (dev 490 / 490, liftover 395 / 395, dev-2 611 / 611).
  Path O2 gets 67 / 303, 66 / 287 and 115 / 273.
- Locus O2 gains nothing (x0.99, x0.98). The gain is not genome alignment as such but the family's structures with their PSVs.
- Over dev-2's four original worlds alone union O2 is +9.7%, below register 1298's 10% bar; hence U2's 5% below.

O3, families with a detection or a false call:

| set | scenario | CUR | S23 | S123m | G |
|---|---|---|---|---|---|
| dev | BRANCH detected | 0 | 8 | 8 | 8 |
| dev | REMOVED detected (4 detectable) | 3 | 3 | 3 | **4** (GAGE: the O2 restriction) |
| dev | false-call families, INTACT + OUTSIDE + CROSS | 4 | 0 | 0 | 0 |
| dev-2 | BRANCH detected | 0 | 8 | 8 | 8 |
| dev-2 | REMOVED detected (2 detectable) | 0 | 1 | 2 | **2** (GOLGA8: the mosaic rule) |
| dev-2 | false-call families of 32 | 13 | 0 | 0 | 0 |

## Bars (held out: the 8 families above, seed 31, five worlds; wrong-family calls count as wrong)

| id | bar |
|---|---|
| U1 | union O2 precision >= 0.995 in every world (reads of catalogued copies) |
| U2 | union O2 correct >= 1.05 x the current O2's, pooled over the five worlds |
| U3 | union O2 correct >= the current O2's in every world |
| U4 | union O2 precision >= 0.99 on CROSS's paralog-structure reads (register 1299's X1, which path O2 failed) |
| U5 | union O2 assigns no more reads of non-catalogued copies than the current O2, pooled |
| G1 | G detects >= 6 of the 8 BRANCH absent copies |
| G2 | G has a false call in <= 1 of the 24 INTACT, OUTSIDE and CROSS families |
| G3 | G detects every detectable REMOVED copy CUR detects, and >= 0.75 of the detectable ones; needs >= 3 detectable |
| G4 | G detects at least as many absent copies (REMOVED + BRANCH) as S23, with no more false-call families |

Detectable (G3) is as before: the removed copy's reads have a median primary `de` on the masked genome > delta. Detection and
false calls use `score.py`'s rule (register 1294).

## Readings (fixed now)

1. **U1-U5 pass:** claimed, "on simulated reads of handpicked real families, O2 on a family graph that threads every isoform
   structure through every copy assigns >= 5% more reads correctly than the current O2, at >= 99.5% precision, including
   when isoforms are annotated on another paralog".
2. **G1, G2, G4 pass and G3 passes:** claimed, "graph O3 detects novel-exon absent copies and >= 75% of detectable removed
   copies, with at most one false call over intact, uncatalogued-paralog and paralog-annotation families". G3 not measurable:
   claimed without the removed-copy clause.
3. Whether the O2 restriction and the mosaic rule ADD detections: G detects more than S23 (G4 holding) -> "they add
   detections on held-out families"; equal -> "they add nothing on this set" (reported, not a failure).
4. Each claim stands or falls separately. A failed bar is reported by name.

Reported beside the bars: path O2 and locus O2; G against S123m, S123 and S23 per family; union O2 against register 1298's 10%;
union O2's gains split into reads whose runner-up copy also carries the read's structure (separated by PSVs) and reads whose
runner-up copy cannot carry it (separated by structure).

Scope either way: simulated reads (curated CAT isoforms, HiFi errors) on the real CHM13 genome. Real data's other species,
fragments and allele background are not modelled.

## Results

Added after the runs. Nothing above this heading changes once a statistic has been seen.
