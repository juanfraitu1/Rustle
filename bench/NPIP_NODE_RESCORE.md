# NPIP node arms re-scored: SCORER v2 (one-to-one copy mapping) + GATE v2 (agreement with guided A1)

DEVELOPMENT ONLY. NPIP is the development family, on copy-centred windows. Nothing here counts as validation.
Nothing committed, `src/` untouched. Pre-registered in `sd_rescore/DECLARATIONS.txt` (08:35:28, md5 checked). The G4
node list was written at 08:36 and results at 08:38. No addenda.

## Read first

The NPIP columns use triangle families at FAMILY level, scored with SCORER v2. The gate columns use GATE v2 with
singletons (167 A1 genes). A gate pass requires R >= 0.740 and F >= 0.818 (G0 minus 0.02).

| arm | v2 R | v2 F | v2 P strict | v2 F strict | MISSING (by collision) | M1 F strict (old) | M1 copies sharing a node | gate R | gate F | gate P strict | gate |
|---|---|---|---|---|---|---|---|---|---|---|---|
| A1 | 1.000 | 1.000 | 0.675 | 0.806 | 0 (0) | 0.806 | 0 | 1.000 | 1.000 | 1.000 | PASS (sanity) |
| G0 | 1.000 | 1.000 | 0.290 | 0.450 | 0 (0) | 0.450 | 0 | 0.760 | 0.838 | 0.483 | PASS (ref) |
| F1 | 1.000 | 1.000 | 0.300 | 0.462 | 0 (0) | 0.462 | 0 | 0.760 | 0.838 | 0.492 | PASS |
| G3 | 1.000 | 1.000 | 0.321 | **0.486** | 0 (0) | 0.486 | 0 | 0.766 | 0.839 | 0.527 | PASS |
| G4 | 1.000 | 1.000 | 0.321 | **0.486** | 0 (0) | 0.486 | 0 | 0.760 | 0.830 | 0.527 | PASS |
| G1a | 0.704 | 0.826 | 0.339 | 0.458 | 7 (7) | 0.581 | 19 | 0.551 | 0.652 | 0.544 | fail |
| G1b | 0.630 | 0.773 | 0.586 | 0.607 | 8 (8) | 0.818 | 19 | 0.503 | 0.618 | 0.618 | fail |
| G2a | 0.778 | 0.875 | 0.280 | 0.412 | 5 (5) | 0.482 | 19 | 0.545 | 0.636 | 0.442 | fail |
| G2b | 0.667 | 0.800 | 0.400 | 0.500 | 6 (6) | 0.659 | 19 | 0.515 | 0.616 | 0.462 | fail |

**Decision: NO CANDIDATE.** Three arms pass the gate: F1, G3 and G4. G3 and G4 tie for the best v2 F strict at
0.486, which is only 0.036 above G0's 0.450 (the bar is 0.05).

Checks that passed:
- A1 scores 1.000 on every gate metric, with and without singletons, so the gate is not broken.
- The frozen M1 numbers match the stored `sd_readgroup/results.json` and `sd_fold` values exactly.
- Every recomputed family list matches the stored one.

## G4 build (spliced-read strand on shipped reps)

| step | count |
|---|---|
| reads (primary, MAPQ >= 1) | 52,917 |
| unspliced reads: strand changed / kept / no spliced read overlapping | 1,486 / 8,442 / 501 (same as G1b, checked in code) |
| single-exon reps: + to - / stay + with reads / stay + with no reads | 225 / 134 / 12 (G3 was 215 / 144 / 12) |
| nodes after consolidate / read-locus nodes added / total | 442 / 24 / 466 |
| edges / triangle families / genes in families / largest family | 1,587 / 59 / 301 / 84 |
| new queries mapped (tx / body) | 26 (65.8 kb) / 1 (343 kb) |

## SCORER v2, all levels (triangle families; components shown as F / F strict)

| arm | level | R | P | F | P strict | F strict | pair sens | pair prec strict | comp F | comp F strict |
|---|---|---|---|---|---|---|---|---|---|---|
| A1 | FAMILY | 1.000 | 1.000 | 1.000 | 0.675 | 0.806 | 1.000 | 0.450 | 1.000 | 0.540 |
| A1 | SUB-1 | 0.704 | 0.704 | 0.704 | 0.475 | 0.567 | 1.000 | 0.255 | 0.704 | 0.380 |
| A1 | SUB-2 | 0.148 | 0.148 | 0.148 | 0.100 | 0.119 | 1.000 | 0.031 | 0.148 | 0.080 |
| G0 | FAMILY | 1.000 | 1.000 | 1.000 | 0.290 | 0.450 | 1.000 | 0.082 | 1.000 | 0.196 |
| G0 | SUB-1 | 0.704 | 0.704 | 0.704 | 0.204 | 0.317 | 1.000 | 0.047 | 0.704 | 0.138 |
| G0 | SUB-2 | 0.148 | 0.148 | 0.148 | 0.043 | 0.067 | 1.000 | 0.006 | 0.148 | 0.029 |
| F1 | FAMILY | 1.000 | 1.000 | 1.000 | 0.300 | 0.462 | 1.000 | 0.088 | 1.000 | 0.205 |
| F1 | SUB-1 | 0.704 | 0.704 | 0.704 | 0.211 | 0.325 | 1.000 | 0.050 | 0.704 | 0.144 |
| F1 | SUB-2 | 0.148 | 0.148 | 0.148 | 0.044 | 0.068 | 1.000 | 0.006 | 0.148 | 0.030 |
| G3 | FAMILY | 1.000 | 1.000 | 1.000 | 0.321 | 0.486 | 1.000 | 0.101 | 1.000 | 0.234 |
| G3 | SUB-1 | 0.704 | 0.704 | 0.704 | 0.226 | 0.342 | 1.000 | 0.057 | 0.704 | 0.165 |
| G3 | SUB-2 | 0.148 | 0.148 | 0.148 | 0.048 | 0.072 | 1.000 | 0.007 | 0.148 | 0.035 |
| G4 | FAMILY | 1.000 | 1.000 | 1.000 | 0.321 | 0.486 | 1.000 | 0.101 | 1.000 | 0.226 |
| G4 | SUB-1 | 0.704 | 0.704 | 0.704 | 0.226 | 0.342 | 1.000 | 0.057 | 0.704 | 0.159 |
| G4 | SUB-2 | 0.148 | 0.148 | 0.148 | 0.048 | 0.072 | 1.000 | 0.007 | 0.148 | 0.033 |
| G1a | FAMILY | 0.704 | 1.000 | 0.826 | 0.339 | 0.458 | 0.487 | 0.111 | 0.851 | 0.237 |
| G1a | SUB-1 | 0.556 | 0.750 | 0.638 | 0.263 | 0.357 | 0.508 | 0.066 | 0.667 | 0.188 |
| G1a | SUB-2 | 0.333 | 0.375 | 0.353 | 0.148 | 0.205 | 0.458 | 0.007 | 0.346 | 0.103 |
| G1b | FAMILY | 0.630 | 1.000 | 0.773 | 0.586 | 0.607 | 0.390 | 0.260 | 0.826 | 0.345 |
| G1b | SUB-1 | 0.481 | 0.722 | 0.578 | 0.433 | 0.456 | 0.382 | 0.144 | 0.596 | 0.252 |
| G1b | SUB-2 | 0.333 | 0.375 | 0.353 | 0.180 | 0.234 | 0.250 | 0.011 | 0.314 | 0.139 |
| G2a | FAMILY | 0.778 | 1.000 | 0.875 | 0.280 | 0.412 | 0.598 | 0.076 | 0.898 | 0.175 |
| G2a | SUB-1 | 0.593 | 0.727 | 0.653 | 0.211 | 0.311 | 0.603 | 0.043 | 0.680 | 0.134 |
| G2a | SUB-2 | 0.296 | 0.320 | 0.308 | 0.099 | 0.148 | 0.500 | 0.004 | 0.302 | 0.062 |
| G2b | FAMILY | 0.667 | 1.000 | 0.800 | 0.400 | 0.500 | 0.436 | 0.111 | 0.875 | 0.227 |
| G2b | SUB-1 | 0.519 | 0.737 | 0.609 | 0.286 | 0.368 | 0.442 | 0.064 | 0.653 | 0.172 |
| G2b | SUB-2 | 0.333 | 0.375 | 0.353 | 0.110 | 0.165 | 0.292 | 0.005 | 0.302 | 0.084 |

Copies left MISSING under v2, by arm:
- G1a: NPIPA1, A7, A8, B8, B10P, B12, LOC124907834
- G1b: NPIPA7, A8, LOC128966608, B7, B8, B10P, B12, LOC124907834
- G2a: NPIPA7, A8, B8, B10P, LOC124907834
- G2b: NPIPA7, A8, LOC128966608, B8, B10P, LOC124907834

Every missing copy is a collision: its best node was already taken by a neighbouring copy.

## GATE v2, whole substrate (167 A1 genes on chr16/17/18, 23 A1 families)

| arm | with singletons: R | P | F | P strict | F strict | without: R | P | F | P strict | F strict | A1 genes unmatched | nodes | families | genes in families | largest family |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A1 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0 | 167 | 23 | 119 | 40 |
| G0 | 0.760 | 0.934 | 0.838 | 0.483 | 0.591 | 0.664 | 0.898 | 0.763 | 0.391 | 0.492 | 27 | 560 | 68 | 364 | 93 |
| F1 | 0.760 | 0.934 | 0.838 | 0.492 | 0.598 | 0.672 | 0.909 | 0.773 | 0.396 | 0.498 | 28 | 539 | 65 | 348 | 90 |
| G3 | 0.766 | 0.928 | 0.839 | 0.527 | 0.624 | 0.689 | 0.911 | 0.785 | 0.439 | 0.536 | 30 | 483 | 60 | 313 | 84 |
| G4 | 0.760 | 0.914 | 0.830 | 0.527 | 0.623 | 0.681 | 0.890 | 0.771 | 0.438 | 0.533 | 30 | 466 | 59 | 301 | 84 |
| G1a | 0.551 | 0.800 | 0.652 | 0.544 | 0.548 | 0.429 | 0.773 | 0.551 | 0.402 | 0.415 | 67 | 259 | 22 | 144 | 56 |
| G1b | 0.503 | 0.800 | 0.618 | 0.618 | 0.554 | 0.361 | 0.768 | 0.491 | 0.462 | 0.406 | 83 | 191 | 16 | 99 | 29 |
| G2a | 0.545 | 0.765 | 0.636 | 0.442 | 0.488 | 0.412 | 0.710 | 0.521 | 0.302 | 0.349 | 57 | 379 | 34 | 225 | 75 |
| G2b | 0.515 | 0.768 | 0.616 | 0.462 | 0.487 | 0.370 | 0.710 | 0.486 | 0.314 | 0.340 | 72 | 297 | 25 | 174 | 45 |

## Per-copy tables (v2 mapping; "shared" = the old M1 node also served another copy)

**G0, G3, G4** (G3 and G4 tie as the best gate-passing arms). All 27 copies are matched, all fall in GWFAM0, and none
share an M1 node. The v2 node equals the M1 node for every copy. The table shows bp of copy exons covered (out of the
copy's exon bp), and the node strand is marked where it differs by arm.

| copy | str | copy bp | G0 node (bp) | G3 node (bp) | G4 node (bp) |
|---|---|---|---|---|---|
| NPIPB2 | - | 1702 | 11963703-11968998 + 2ex read-locus (416) | same as G0 (416) | **11977893-11978376 + 1ex (63)** |
| NPIPA2 | + | 4517 | 14749588-14763868 + 7ex (2827) | = | = |
| NPIPA1 | + | 1084 | 14938814-14953123 + 7ex (1021) | = | = |
| PKD1P6-NPIPP1 | - | 5391 | 15115422-15132313 - 10ex (1898) | = | = |
| NPIPA5 | - | 3741 | 15368407-15382704 - 7ex (2915) | = | = |
| NPIPA6 | + | 1570 | 16344720-16359014 + 8ex (1088) | = | = |
| NPIPA7 | + | 1264 | 16391852-16406187 + 7ex (1093) | = | = |
| NPIPA8 | - | 1588 | 18325177-18339582 - 7ex (1080) | = | = |
| NPIPA9 | - | 2869 | 18372354-18376104 - 2ex (704) | = | = |
| NPIPB3 | - | 3674 | 21337422-21344055 - 3ex (2989) | = | = |
| LOC128966608 | - | 3557 | 21692009-21693757 **+** 1ex (129) | - (129) | - (129) |
| NPIPB4 | + | 3634 | 22382178-22387825 + 2ex (184) | = | = |
| NPIPB5 | + | 9237 | 22788631-22792526 + 1ex (3031) | = | = |
| NPIPB6 | - | 7025 | 28625753-28628690 **+** 1ex (2635) | - (2635) | - (2635) |
| NPIPB7 | - | 2952 | 28736995-28771999 - 8ex (1274) | ...-28772246 (1274) | ...-28772246 (1274) |
| NPIPB8 | + | 1293 | 28935124-28939424 + 3ex (781) | = | = |
| NPIPB9 | + | 4193 | 29048645-29050801 + 1ex (520) | = | = |
| NPIPB10P | + | 880 | 29328121-29333923 + 1ex (441) | = | = |
| NPIPB11 | - | 6117 | 29663364-29679826 - 7ex (3236) | = | = |
| NPIPB12 | - | 3737 | 29765341-29767376 **+** 1ex (1783) | - (1783) | - (1783) |
| LOC124907834 | - | 3244 | 30507437-30523997 - 7ex (3120) | ...-30528300 (3120) | ...-30528300 (3120) |
| NPIPB13 | - | 6144 | 30609406-30625966 - 7ex (3108) | ...-30628027 (3108) | ...-30628027 (3108) |
| NPIPB14P | - | 2418 | 75785725-75790348 - 2ex (1574) | = | = |
| NPIPB15 | + | 4340 | 80195301-80209898 + 7ex (1605) | = | = |
| LOC124907808 | + | 4981 | 80319739-80324246 + 2ex (2777) | = | = |
| LOC124907807 | + | 4571 | 80424045-80438604 + 7ex (1805) | = | = |
| NPIPB1P | - | 981 | chr18:11781592-11796212 - 7ex (810) | = | = |

In G4 the 2-exon read-locus node over NPIPB2 disappears. With the restranded reads, a new 4-read `+` read-locus node
appears at 11951326-11952692 instead, and NPIPB2 falls back to a 63 bp single-exon placeholder. The NPIP family
numbers do not change.

**G1b** is shown for contrast (the old M1 winner, which fails the gate). It lists every copy that is missing or out of
the family, plus the other copies whose node also served a neighbour under M1.

| copy | v2 node (bp) | v2 family | M1 node |
|---|---|---|---|
| NPIPA2 | 14740775-15167445 + 85ex (4472) | GWFAM0 | same node, shared with NPIPA1 |
| NPIPA1 | 14948058-14988682 - 11ex (389) | GWFAM2 (out) | the NPIPA2 node |
| NPIPA6 | 16328099-16438449 + 50ex (1569) | GWFAM0 | shared with NPIPA7 |
| NPIPA7 | none | MISSING (collision) | the NPIPA6 node |
| NPIPA9 | 18147539-18493395 - 68ex (2869) | GWFAM0 | shared with NPIPA8 |
| NPIPA8 | none | MISSING (collision) | the NPIPA9 node |
| NPIPB3 | 21314626-21842398 - 134ex (3668) | GWFAM0 | shared with LOC128966608 |
| LOC128966608 | none | MISSING (collision) | the NPIPB3 node |
| NPIPB5 | 22216076-22958833 + 126ex (9229) | GWFAM0 | shared with NPIPB4 |
| NPIPB4 | 22233283-22365456 - 25ex (207) | GWFAM2 (out) | the NPIPB5 node |
| NPIPB6 | 28623009-28772328 - 43ex (6325) | GWFAM0 | shared with NPIPB7 |
| NPIPB7 | none | MISSING (collision) | the NPIPB6 node |
| NPIPB9 | 28772052-29800624 + 210ex (4105) | GWFAM0 | shared with NPIPB8 and NPIPB10P |
| NPIPB8, NPIPB10P | none | MISSING (collision) | the NPIPB9 node |
| NPIPB11 | 29663336-29949744 - 85ex (4550) | GWFAM0 | shared with NPIPB12 |
| NPIPB12 | none | MISSING (collision) | the NPIPB11 node |
| NPIPB13 | 30437930-30741477 - 111ex (6144) | GWFAM0 | shared with LOC124907834 |
| LOC124907834 | none | MISSING (collision) | the NPIPB13 node |

## Readings

1. **Defect 1 is confirmed and fixed.** Under v2, G1b falls from F strict 0.818 to 0.607. R drops to 0.630 because 8
   copies are missing and 2 more land outside the family. Its 9 shared nodes (the M1 nodes serving 19 copies) span
   0.11-1.03 Mb, cover 2-3 copies each, and have 43-210 exons. The read-group arms (G1a/G1b/G2a/G2b) lose 5-8 copies each, all by collision. The catalog arms (A1,
   G0, F1, G3, G4) do not change at all: none has two copies sharing a node, so v2 and M1 give the same numbers.
2. **Defect 2 is fixed.** GATE v2 gives A1 1.000 by construction, and it separates the arms clearly. The catalog arms
   score 0.760-0.766 R and 0.830-0.839 F. The read-group arms score 0.503-0.551 R and 0.616-0.652 F. They fail because
   whole-substrate genes merge into long read-group nodes and 57-83 of 167 A1 genes go unmatched.
3. **The G1b strand fix does not rescue shipped reps (G4).** Moving the strand evidence from raw reads to
   spliced-restranded reads flips 10 more single-exon reps to `-` (225 vs 215). Even so, the NPIP numbers are
   identical to G3, and substrate-wide agreement is slightly lower (F 0.830 vs 0.839). G4 adds nothing over G3.
4. **The best honest gain is G3's strand fix at +0.036 F strict**, which is below the 0.05 bar. All of that gain is
   precision: 84 vs 93 nodes in the matched family. G3 also raises gate P strict over G0 (0.527 vs 0.483), so the
   gain is not NPIP-only. That P strict is not the best of any arm: G4 is marginally higher (0.52697 vs 0.52675), and
   the gate-failing read-group arms G1a (0.544) and G1b (0.618) are higher still. The gap to A1 (0.806) is still 0.32.

## Provenance

- Repo `c45c3261`, read-only. `shared_definition.rs` md5 5d16fc38. The `lib.py` logic is the verified mirror
  (`npip_ladder/verify`), copied through `sd_fold/charz/lib.py`.
- Data paths are relative to `/mnt/linuxdisk/home/juanfraitu/`:
  - `sd_rescore/DECLARATIONS.txt` (+ `.md5`)
  - `sd_rescore/scripts/build_g4.py` → `g4.pkl`, `logs/build_g4.log`
  - `sd_rescore/scripts/rescore.py` → `results.json`, `logs/rescore.log`
- Mapping used only the prebuilt indexes, for 27 new md5s:
  - `minimap2 2.30 -c -N 50 -p 0.1 -x splice -uf -t 4 npip_ladder/idx/target.splice.mmi` for `sd_rescore/map/tx_chunk0`
  - `-x asm20 ... target.asm20.mmi` for `body_chunk0`
- Inputs reused: `npip_ladder/verify/{nodes,edges}.pkl`, `npip_ladder/union/*.paf`, `sd_fold/map/new_*.paf`,
  `sd_readgroup/{arms,families}.pkl` + `map/new_*.paf`, `sd_fold/charz/{variants,results}.pkl`.
- Asserts that passed:
  - Truth exon_bp equals `npip_per_member/wholechr16/truth_exons.tsv`.
  - Recomputed families equal the stored ones (G0 G3 G1a G1b G2a G2b F1; A1 is stored).
  - Frozen M1 core rows equal the stored ones (triangle and components, 3 levels) for 7 arms plus F1.
  - The G4 read restrand counts equal G1b's.
  - A1 gate = 1.000.

## Caveats

- The threshold of v2 R >= G0's R is not binding here. G0 and every gate-passing arm have R = 1.000.
- The Hungarian mapping picks the copy-node pairing that maximises total bp. When one node covers several copies, it
  goes to whichever copy gives the largest total, which is not always the "right" copy (see NPIPA1 vs NPIPA2 in
  G1b). Missing copies are penalised in R, not in P.
- The gate's reference is only as good as A1. A1 is limited to RefSeq genes with >= 3 reads, and its edges and
  triangle families have their own errors, so "agreement" is not "truth". The gate is relative (G0 minus 0.02), so it
  asks "no worse than shipped", not "good". G0 itself agrees only 0.760 / 0.838.
- "Without singletons" P strict counts a node matched to a singleton A1 gene as a non-reference member. This was a
  declared choice. It lowers P strict for every arm except A1 and G0: G0 is unchanged at 0.391, and the other arms
  drop by 0.002-0.038.
- NPIP is a development family on copy-centred windows. G3's +0.036 and the G1b collapse both need a separate,
  pre-registered held-out test before anyone quotes them.

## Verification (independent recompute)

An independent verifier checked this report. It did not read the builder's scripts. Its code and logs are in
`sd_rescore/verify/` (`v_g4.py`, `v_score_v2.py`, `v_extras.py`, `v_samtools_g4.py`, `v_alt_without.py`, plus `*.log`).
It reuses the already-verified mirrors in `npip_ladder/verify/` and `sd_readgroup/verify/`.

1. **Declarations came first.** `DECLARATIONS.txt` (08:35:28) is older than `g4.pkl` and the map FASTAs (08:36:15),
   the PAFs (08:37) and `results.json` (08:38:50). Its md5 (d605c100…) matches `DECLARATIONS.md5`. There are no
   addenda.
2. **G4 nodes, rebuilt from scratch.** The verifier started from the ladder-verified reads and the dump, applied the
   declared restrand rules, then ran consolidate and read-locus. The result is 466 nodes (442 + 24), and 0 rows differ
   from `g4.pkl` on chrom, strand, n_reads, exons or rep_exons. Read restranding gives 1,486 / 8,442 / 501, identical to
   the verified G1b list. Rep flips give 225 / 134 / 12, against G3's 215 / 144 / 12.
3. **G4 edges and families, recomputed from the PAFs.** The new map FASTAs hold 26 tx queries (65,788 bp) and 1 body
   query (343,201 bp). Every md5 is correct, none was already in an older source, and together they are exactly the
   queries G4 needs. The recompute gives 1,587 pairs, 59 triangle families, 301 nodes in families, and a largest
   family of 84. F1's pairs and its triangle and component families equal `sd_fold/charz/results.pkl`.
4. **SCORER v2 and GATE v2, reimplemented.** The verifier wrote its own Hungarian assignment: items × (nodes + one dummy
   per item), integer cost -w·BIG + node row, zero-weight pairs forbidden. It then recomputed all 9 arms:
   - All 82 numeric leaves per arm match `results.json` to full precision (v2 and M1; triangle and components; 3
     levels; the gate with and without singletons).
   - Pair sensitivity and pair precision strict also match.
   - The missing and collision lists, copies sharing an M1 node, unmatched A1 genes, node, family and largest-family
     counts, and all per-copy nodes, bp and families match. The only differences are cosmetic label formats.
   - Frozen M1 numbers equal the stored `sd_readgroup/results.json`, `sd_fold` F1 and ladder A1 values.
   - A1 scores 1.000 on every gate metric in both versions and passes.
   - The decision is reproduced: the gate-passing arms are F1, G3 and G4. The top v2 F strict is 0.486486, tied by G3
     and G4. That is 0.0365 above G0, below the 0.05 bar, so there is **NO CANDIDATE**.
   - Every missing copy's M1 node is matched to another copy (7/7, 8/8, 5/5, 6/6).
5. **samtools spot-check of G4 restranding.** The verifier read `reads.bam` with `-F 2308 -q 1` and applied the ts:A:-
   flip and the spliced-read restrand. It then took a 2/3 majority on the dump single-exon rep for LOC128966608,
   NPIPB6 and NPIPB12 (→ -) and for NPIPB5 and NPIPB9 (→ +). All 5 agree with the G4 node strand and the copy strand.

Sentences corrected. None of the numbers were wrong.
- Reading 1 said G1b's "nodes span 0.1-1 Mb and cover 2-3 copies each". That is true only of its 9 shared M1 nodes,
  which serve 19 copies. The other 8 copies each have their own node.
- Reading 4 said "G3 also gives the best gate P strict". It does not. G4 is marginally higher (0.52697 vs 0.52675), and
  G1a (0.544) and G1b (0.618) are higher still.
- A caveat said the without-singletons P strict choice "lowers every arm except A1". G0 is unchanged (0.3911 either
  way).
