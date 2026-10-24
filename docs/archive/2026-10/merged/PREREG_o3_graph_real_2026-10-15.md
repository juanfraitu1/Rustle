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

*Superseded by "Verdict, revised after the review" at the end of this document.*

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

## Review corrections (whole-branch review, 2026-10-15; post hoc)

A fresh reviewer re-ran every number above and all of them reproduce. The prereg timing, the absence of leakage (anonymised
names; no method reads the truth) and the simulation's construction check out. The corrections below are post hoc; the bars
and their values above stand as computed.

**C1. About 42% of path O2's gain is decided by the annotation, not the sequence.**
- `o2graph.py seqsupport` aligns every hosted read to each catalogued copy's GENOMIC locus (spliced, +-1 kb). A read is
  genome-separated when its own copy's locus explains it at least MARGIN (2) edits better than any other copy's locus.
- 9,330 of the 12,766 reads of catalogued copies are genome-separated. On them path O2 makes 9,326 correct calls and the
  current O2 8,602: **+8.4%, below O2b's 10% bar.**
- Of path O2's 1,244 additional correct calls, 724 are genome-separated. The other 520 are reads another copy's locus
  explains within 1 edit; only the curated isoforms listed for each copy separate them. Every read is simulated from its
  copy's own curated isoforms, so those calls are correct by construction.
- Margin 1's extra calls are all of that kind (genome-separated correct: 9,326 at both margins).
- A read that aligns to one copy's isoform paths only is assigned without the margin test (`decide`); TBC1D3's gain is of
  this kind (reviewer).
- Representatives only (an incomplete catalog): 7,864 genome-separated reads correct, fewer than the current O2's 8,602.
- Hence Amendment 1 (above) and its results (below).

**I1. S4's O2 step adds nothing on the held-out set.** S4 equals S23 (branch calls + genome check, without the isoform step
or path O2) in every held-out family and scenario. Every BRANCH detection and the TBC1D3 REMOVED detection are branch calls;
every removed false call is removed by the genome check. Development's one difference (REMOVED GAGE, S4 only) did not
replicate. The O3 result belongs to the branch calls and the genome check, not to "O2-informed O3".

**I2. The current method's 13 false-call families are mostly another family's reads.**
- 12 of the 13 are GOLGA6L, NPIPA and NPIPB in each of the four worlds: reads of a sister family placed on them. Splitting
  families by gene-name pattern creates those sisters (NPIPA / NPIPB, GOLGA8 / GOLGA6L).
- Only TBC1D3 in OUTSIDE is the designed failure. There the current method scans the uncatalogued paralog only because the
  region hull includes the special copy (M2, below).
- The current method plus the genome check (S3) is also 0 / 32.
- copy_assign's opt-in `--psv-family-best` (register 1296) targets exactly such reads; it was not run here.
- So O3c / O3d's false-call comparison shows what the genome check does. It is not an advantage of S4 over a
  genome-checked current method.

**I3. Fixed:** S4 keyed path O2's unexplained reads by read only. A read unexplained in a sister family it was hosted in
entered its own family's S4 pool (about 550 reads per world). It is now keyed by (read, family), with a test. Re-run on all
four held-out worlds, every call is identical.

**I4. O2c (181 vs 426) comes mostly from BRANCH.** Reads of non-catalogued copies assigned, current vs path O2: BRANCH 240 vs 9,
REMOVED 119 vs 86, OUTSIDE 67 vs 86 (path O2 places more), INTACT 0 vs 0. In BRANCH the 150-bp novel exon puts the reads
past delta.

**M1 (graded up: a reported number was wrong). Fixed:** wrong-family assignments were not scored. The current O2 makes 56
(NPIPA 40, NPIPB 16), so its "wrong 0" was wrong: its precision is 0.9936. Path O2 makes none (1.0). Counting reads of
non-catalogued copies too: path O2 0.982, current 0.948.

**M3 (graded up: the O3b reading rests on it).** GOLGA6L's removed copy has median `de` 0.00935, against delta 0.00958. Had it
been above delta, O3b would have been measurable (3 detectable) and failed (S4 1 of 3).

**M4 (graded up: the only REMOVED detection).** S4's TBC1D3 detection is one branch call holding 50 of 60 source reads (0.83
against the 0.8 bar). The current method made two pure calls, of 17 and 42 reads, split across host copies; neither
reaches the bar.

**M5. Corrected:** CT47A's copies are not identical over their reads. At margin 1 path O2 assigns 480 CT47A reads correctly; at
margin 2, none.

**M7. Development provenance.**
- Development INTACT, OUTSIDE and REMOVED were built (18:59-19:17) before `realworld.py`'s final version (19:34).
- Development INTACT's path-O2 output predates the final `o2graph.py` by 24 s.
- Rebuilt with the final code, the three worlds' reads, truth, catalog and regions are byte-identical. Path O2 reproduces
  development INTACT's output exactly.

**Deferred (minor, not fixed):**
- M2: regions are built from every copy, the special one included, so copy_assign scans the OUTSIDE paralog (see I2).
- M6: the representatives-only run overwrites path O2's `graph/o2paths/*.tx.fa` intermediates; outputs are unaffected.
- M8: `families.py`'s tie-breaker does nothing.

### Amendment 1 results (held out, world CROSS, 2026-10-15)

Run: code `9b6b98df`, prereg amendment `99932167`. 4,372 reads of catalogued copies; 940 of them come from 79
paralog-annotated structures over 57 copies (475 paralog isoforms projected). CT47A, GOLGA6L and XAGE (9 reads) have
almost none.

| | assigned | correct | wrong | precision | paralog-structure reads: assigned / correct |
|---|---|---|---|---|---|
| current O2 | 2,881 | 2,810 | 71 (all wrong-family: NPIPA 58, NPIPB 13) | 0.975 | 538 / 538 |
| path O2, margin 2 | 2,901 | 2,743 | 158 | 0.946 | 273 / 115 |
| path O2, margin 1 | 3,167 | 2,916 | 251 | 0.921 | 367 / 116 |
| path O2, representatives only | 2,250 | 2,232 | 18 | 0.992 | 95 / 83 |

- Genome-separated reads (3,036): correct, current 2,763 and path O2 2,535.
- Path O2 calls 1,350 reads `unexplained`; S4 forms no call from them once the genome check runs.
- O3: S4, S123, S23 and S3 have 0 false-call families of 8. CUR, R0, S1, S12 and S4n have 3 (GOLGA6L, NPIPA, NPIPB, the
  sister-family reads of I2).

| bar | value | |
|---|---|---|
| X1 | 115 / 273 = 0.42 | **fail** |
| X2 | 2,743 / 2,901 = 0.946 | **fail** |
| X3 | 2,743 vs 2,810 | **fail** |
| X4 | 0 of 8 | pass |

**Reading 2.** When a copy expresses a structure annotated only on a paralog, path O2 names the annotated paralog. Its
claim is scoped to a complete, copy-specific isoform catalog, and it is not safe on real annotation as it is.

### Verdict, revised after the review (2026-10-15)

**O2: claimed only for a complete, copy-specific isoform catalog; not for real annotation.**
- O2a-O2d pass as computed, but only on the simulation's premise that each copy's catalog entry lists exactly the isoforms
  it expresses.
- Counting only reads the genome itself separates, path O2's gain is +8.4%, below O2b's 10%.
- When copies also express paralog-annotated structures (Amendment 1), path O2 is wrong on 58% of the reads it assigns
  from them. Its overall precision drops to 0.946, and it assigns fewer reads correctly than the current O2 (X1-X3 fail).
- With representatives only (an incomplete catalog), it gets fewer genome-separated reads right than the current O2.
- The current O2's precision is 0.9936 with wrong-family calls counted (0.975 in CROSS; all of them sister-family calls).
- Next, under its own prereg on a new family set: thread every family structure through every copy (each copy's paths =
  its own isoforms + every paralog isoform projected onto its locus), so that only the sequence separates copies. Its
  ceiling on this set is the genome-separated gain (+8.4%).

**O3: the branch calls and the genome check, not O2-informed O3.**
- BRANCH 8 / 8 (current 0) holds, from the branch calls.
- No false call over 32 families (and 0 / 8 in CROSS), from the genome check. The current method plus the genome check
  does as well; the current method's own false calls are mostly sister-family reads.
- REMOVED is not measurable: 6 of 8 removed copies are within delta, one of them by a hair (0.00935 vs 0.00958).
  TBC1D3's detection is marginal (0.83 vs 0.8).
- The GOLGA8 mosaic rule remains the next O3 step (new prereg).

Scope unchanged: simulated reads (curated CAT isoforms, HiFi errors) on the real CHM13 genome.
