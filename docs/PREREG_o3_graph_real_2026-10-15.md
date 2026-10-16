# Graph-based O2 and O3 on handpicked real families (simulation v2, 2026-10-15)

Pre-registration, written before any held-out family's world was built. Code: branch `o3-graph`, commit `3dcf98a3`
(`bench/o3_graph/`).

**User goal (2026-10-15):** "improve the simulation and use handpicked families (maybe soto or NPIP or other mainstays) to
improve the graph based O2 and O3".

This follows `docs/PREREG_o3_graph_sim_2026-10-15.md` (register 1297), whose worlds were famsim copies of single curated
genes in small genomes.

## Simulation v2 (`families.py`, `realworld.py`)

**Families.** Handpicked mainstay multi-copy families, taken by CAT gene-name pattern from CHM13 v2.0's CAT/Liftoff
annotation (the project's human default):
- copies linked when their canonical transcripts align at >= 90% identity (largest linked component);
- span <= 60 kb, non-overlapping, at most 10 copies (seeded);
- up to 4 curated isoforms per copy with distinct intron chains.

The split was fixed before any run:

| set | families (copies / isoforms) |
|---|---|
| development | SMN 2 / 8, NCF1 3 / 8, PAGE 3 / 8, CT45A 9 / 16, GAGE 10 / 12, POTE 10 / 28, GSTA 5 / 11, TP53TG3 4 / 13 |
| **held out** | NPIPA 6 / 17, NPIPB 10 / 24, TBC1D3 10 / 20, GOLGA8 10 / 24, GOLGA6L 4 / 5, CT47A 7 / 7, XAGE 3 / 12, SPATA31 7 / 14 |

NPIP, TBC1D3 and GOLGA8 are Soto's mainstays.

**Worlds.**
- The real CHM13 genome; reads mapped genome-wide with the baseline's minimap2 command.
- 60 reads per expressed copy over its isoforms: the canonical "main" isoform weighs 2.0, the others 0.5.
- HiFi errors from `bench/sim.py` (0.1% substitutions, 0.03% indels) and 0-30 bp end jitter.
- **Read names are a seeded shuffle and carry no truth** (review M5 of register 1297).
- The catalog's representative is each copy's main isoform.
- One SPECIAL copy per family (seeded) plays the scenario's part.

| scenario | what the special copy is |
|---|---|
| INTACT | an ordinary copy: every copy is in the genome and the catalog |
| REMOVED | hard-masked in the genome and dropped from the catalog; its reads still exist |
| OUTSIDE | in the genome, not in the catalog (an uncatalogued paralog) |
| BRANCH | the source of an extra ABSENT copy: its isoforms with a novel 150-bp exon at the middle of the main isoform's middle intron (spliced into every isoform with an intron across that point) and 0.2% exonic SNPs |

## Methods under test

**O2** (`o2graph.py`; the reads O2 places on the family = copy_assign's rows):

| method | what it does |
|---|---|
| current O2 | copy_assign `--no-as-tied-only`: rows with status `assigned` |
| **path O2** | each read aligned (minimap2 -x map-hifi, all secondaries) to EVERY curated isoform path of every catalogued copy. Distance to a copy = edits + unaligned bases against its best-fitting isoform. Assigned to the closest copy when the next is >= MARGIN farther; abstains when tied or thin, and as `unexplained` when even the closest copy is > delta away. **Primary: MARGIN = 2**; MARGIN = 1 reported. |
| path O2, representatives only | the ablation: each copy's representative path only (MARGIN 2) |

**O3** (`detect.py`; arms as in register 1297, plus S4):

| arm | what it does |
|---|---|
| CUR | `copy_assign --psv-clusters` |
| S123 | all-isoform explanation + branch calls + genome check |
| **S4** | the S1 clustering run only on the reads path O2 calls `unexplained`, plus the branch calls, kept by the genome check |

## Known before running (recorded now)

**Development results**, final code. O2, pooled over the four worlds (10,106 reads of catalogued copies, 1,449 of
non-catalogued ones):

| | assigned | correct | wrong | non-catalogued reads assigned |
|---|---|---|---|---|
| current O2 | 6,484 | 6,484 | 0 | 660 |
| path O2, margin 2 | 7,744 | 7,744 | 0 | 217 |
| path O2, margin 1 | 9,114 | 9,114 | 0 | 226 |
| path O2, representatives only, margin 2 | 6,502 | 6,348 | 154 | 170 |

O3 (8 families per scenario):

| | CUR | S123 | S4 |
|---|---|---|---|
| BRANCH detected | 0 | 8 | 8 |
| REMOVED detected | 3 | 3 | 4 |
| OUTSIDE false-call families | 4 | 0 | 0 |
| INTACT false-call families | 0 | 0 | 0 |

REMOVED's 8 removed copies include 4 that are detectable by O3's definition (the reads' median primary `de` on the masked
genome > delta: CT45A, GAGE, GSTA, PAGE). The 4 within delta are NCF1, POTE, SMN and TP53TG3; no arm calls them.

**Found and fixed on development:**
- copy_assign needs a family's span on a chromosome in one region, so regions are family hulls.
- Node-based graph O2 (allele nodes of a free GraphAligner alignment) makes artifactual calls. In GAGE it made 74 wrong of 289;
  the reads match both copies equally well linearly. Normalising the graph made it worse (173 wrong). Path-constrained
  alignment replaces it.
- A novel exon goes into every isoform whose intron spans it (not only the exact junction).
- The `novel` truth flag must look at both strands.
- S4 was added after GAGE's removed copy formed an impure cluster (59 of its reads with 48 of the host's, purity 0.55).
- BRANCH is scored on the absent copy's reads that carry the novel exon. Its other reads are the source copy's sequence
  within 0.2% SNPs, so no method can tell them apart; the all-reads result is reported beside it.

## Bars (held out: the 8 held-out families, seed 21; pooled over INTACT, REMOVED, OUTSIDE, BRANCH)

| id | bar |
|---|---|
| O2a | path O2 (margin 2) precision over reads of catalogued copies >= 0.995 |
| O2b | path O2 (margin 2) correct assignments >= 1.10 x the current O2's |
| O2c | path O2 (margin 2) assigns no more reads of non-catalogued copies than the current O2 |
| O2d | the representatives-only ablation makes more wrong assignments than all-isoform path O2 (step 1 is what keeps it precise) |
| O3a | S4 detects >= 6 of the 8 BRANCH absent copies |
| O3b | S4 detects every detectable REMOVED copy that CUR detects, and >= 0.75 of the detectable ones; needs >= 3 detectable |
| O3c | S4 has a false call in <= 1 of the 16 INTACT and OUTSIDE families |
| O3d | S4 detects >= 5 more BRANCH copies than CUR, and has no more false-call families than CUR |

Detectable (O3b) means the removed copy's reads have a median primary `de` on the masked genome > delta. That is O3's
definition: nothing in the genome explains the reads within delta.

## Readings (fixed now)

1. **O2a-O2d pass:** claimed, "on simulated reads of handpicked real families, O2 on every curated isoform path assigns
   >= 10% more reads correctly than the current O2 at >= 99.5% precision, and the representatives-only version is less
   precise".
2. **O3a-O3d pass:** claimed, "O2-informed O3 (S4) detects absent copies carrying a novel exon, keeps the current method's
   detections of sequence-divergent absent copies, and makes at most one false call over intact and uncatalogued-paralog
   families".
3. Each claim stands or falls separately. A failed bar is reported by name.

Scope either way: simulated reads (curated isoforms, HiFi errors) on the real genome. Real data's other species, fragments
and allele background are not modelled.

## Results

Added after the runs. Nothing above this heading changes once a statistic has been seen.

## Results (2026-10-15)

Run: code `549b5d1c`, the frozen copy_assign binary. Held-out worlds: the 8 held-out families, seed 21. Outputs in
`/mnt/linuxdisk/tmp/o3_graph/real/heldout/` (one directory per scenario, `o3_score.json`, `o2_*.json`).

**Disclosure: a scorer bug fixed after the first scoring.** Every held-out world carries the truth column `novel`, so
`score.py` took each REMOVED copy's source to be its novel-exon reads, which are none. That turned REMOVED's detections into
"no source" and its calls into false calls. The prereg scores BRANCH on novel-exon reads and REMOVED on all of the copy's
reads (register 1294), so the implementation was corrected to that rule (test added). Development's REMOVED world predates
the column and was unaffected. Only REMOVED's numbers changed.

**O2** (pooled over the four worlds: 12,766 reads of catalogued copies, 1,443 of non-catalogued ones):

| | assigned | correct | wrong | non-catalogued reads assigned |
|---|---|---|---|---|
| current O2 | 8,720 | 8,720 | 0 | 426 |
| **path O2, margin 2** | 9,964 | 9,964 | 0 | 181 |
| path O2, margin 1 | 10,650 | 10,650 | 0 | 213 |
| path O2, representatives only, margin 2 | 8,148 | 8,136 | 12 | 234 |

Per family (correct, current -> path O2 margin 2):

| family | reads | current | path O2 |
|---|---|---|---|
| GOLGA6L | 840 | 408 | 840 |
| NPIPA | 1,324 | 1,040 | 1,300 |
| NPIPB | 2,292 | 1,548 | 1,754 |
| TBC1D3 | 2,288 | 2,072 | 2,164 |
| GOLGA8 | 2,288 | 2,120 | 2,152 |
| SPATA31 | 1,564 | 1,532 | 1,564 |
| XAGE | 610 | 0 | 190 |
| CT47A | 1,560 | 0 | 0 (copies identical over their reads) |

**O3** (8 families per scenario):

| | CUR | S12 | S123 | **S4** |
|---|---|---|---|---|
| BRANCH detected | 0 | 8 | 8 | **8** |
| BRANCH false-call families | 3 | 3 | 0 | **0** |
| INTACT false-call families | 3 | 3 | 0 | **0** |
| OUTSIDE false-call families | 4 | 4 | 0 | **0** |
| REMOVED detected (2 detectable) | 0 | 1 | 1 | **1** |
| REMOVED false-call families | 3 | 3 | 0 | **0** |

- BRANCH holds on all of the absent copy's reads too: S4 8 / 8, CUR 0.
- 6 of the 8 removed copies are within delta of the masked genome by their reads' median `de`: CT47A 0.0016, NPIPB 0.0014,
  XAGE 0.0021, NPIPA 0.0026, SPATA31 0.0057, GOLGA6L 0.0093.
- The 2 detectable ones are TBC1D3 (0.0141; S4 detects it, CUR does not) and GOLGA8 (0.0143; missed by every arm).

| bar | value | |
|---|---|---|
| O2a | precision 1.0 | pass |
| O2b | 9,964 vs 8,720: +14.3% | pass |
| O2c | 181 vs 426 | pass |
| O2d | representatives only: 12 wrong; all isoforms: 0 | pass |
| O3a | BRANCH 8 / 8 | pass |
| O3b | 2 detectable removed copies (needs >= 3) | **not measurable** |
| O3c | 0 of 16 | pass |
| O3d | BRANCH +8; false-call families 0 vs 13 | pass |

### Verdict (2026-10-15)

**O2: reading 1, claimed** (simulated reads of handpicked real families). O2 on every curated isoform path assigns 14.3% more
reads correctly than the current O2, at precision 1.0. It also places fewer reads of non-catalogued copies (181 vs 426).
The representatives-only version makes wrong assignments (12 held out, 154 in development), so the all-isoform content is
what keeps it precise. Margin 1 assigns more still (10,650, no error) but is not the pre-registered primary.

**O3: claimed in part.** O3a, O3c and O3d pass. O2-informed O3 (S4):
- detects every absent copy carrying a novel exon (8 / 8; the current method 0);
- makes no false call over intact, uncatalogued-paralog, removed and branch families (0 / 32; the current method 13 / 32).

The REMOVED clause (O3b) cannot be tested on this set: 6 of its 8 removed copies are within delta of a remaining paralog, as
these young families' copies are. On the two that are not, S4 detects TBC1D3 and the current method neither.

**What the miss shows (post hoc).** GOLGA8's removed copy is a MOSAIC: its 58 reads form a pure cluster with 23 alt sites,
22 of them carried by some other copy, so no copy explains the cluster and too few sites are unexplained (`unresolved`). Its
reads are 1.43% from anything in the genome. A rule on the cluster's distance to its nearest single copy, not on
copy-private alleles, would call it. That is a method change for a new prereg on another family set.

**Scope:**
- simulated reads (curated CAT isoforms, HiFi errors) on the real CHM13 genome;
- real data's other species, fragments and allele background are not modelled;
- the special copies are one per family (seeded).

## Amendment 1 (2026-10-15): a catalog that is not copy-specific (the CROSS world)

Written after the whole-branch review and before the held-out CROSS world was built. Code: `9b6b98df`.

**Why.** The review found that much of path O2's gain is decided by the annotation rather than the sequence (review C1;
the corrections are below, under "Review corrections"). Every read here is simulated from one of its copy's own curated
isoforms, and each copy's catalog entry lists exactly those. Real annotation is not that copy-specific: an isoform
annotated on one paralog can be expressed by another. Then path O2 meets a read whose exon structure only a paralog's
catalog entry carries.

**World CROSS** (`realworld.py --scenario CROSS`):
- as INTACT: every copy is in the genome and the catalog, and its own-isoform reads are INTACT's;
- in addition, every copy expresses up to 2 isoform structures (seeded) annotated only on a paralog;
- such a structure is the paralog isoform's sequence spliced-aligned to the copy's own locus (minimap2 splice:hq). It is
  kept when the primary alignment covers >= 95% of the isoform with the same exon count, canonical junctions and an intron
  chain none of the copy's own isoforms has (junctions within 10 bp count as the same);
- the transcript is the copy's own genomic sequence over those exons, at the weight of a non-main isoform;
- the catalog still lists each copy's own isoforms only. Truth isoform `proj:<paralog>|<isoform>`.

**Known before running (development, seed 11).** 3,539 reads of catalogued copies, 771 of them from 67 paralog-annotated
structures over 46 copies:

| | assigned | correct | precision | paralog-structure reads: assigned / correct |
|---|---|---|---|---|
| current O2 | 2,228 | 2,219 | 0.996 | 400 / 391 |
| path O2, margin 2 | 2,459 | 2,223 | 0.904 | 303 / 67 |
| path O2, margin 1 | 2,864 | 2,620 | 0.915 | 361 / 117 |
| path O2, representatives only | 1,945 | 1,817 | 0.934 | 152 / 67 |

229 of path O2's 236 wrong calls name the paralog that carries the annotation. O3: every arm 0 / 8 false-call families.

**Bars (held out: the 8 held-out families, seed 21, world CROSS; wrong-family assignments count as wrong):**

| id | bar |
|---|---|
| X1 | path O2 (margin 2) precision over the reads of paralog-annotated structures >= 0.99 |
| X2 | path O2 (margin 2) precision over all reads of catalogued copies >= 0.995 (O2a's bar) |
| X3 | path O2 (margin 2) correct assignments >= the current O2's |
| X4 | S4 has a false call in <= 1 of the 8 families |

**Readings (fixed now):**
1. X1-X3 pass: O2's claim extends to a catalog that is not copy-specific.
2. X1 or X2 fails: O2's claim is scoped to a complete, copy-specific isoform catalog, and path O2 is not safe on real
   annotation as it is. A fix (every family structure threaded through every copy) needs its own prereg on a new family set.
3. X4 is read on its own: it pass/fails O3's specificity under such a catalog.
