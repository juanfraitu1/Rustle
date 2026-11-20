# Row 1833: every read counted (no collapse of coordinate duplicates), tested off the spike-ins (pre-registered 2026-11-20)

User (2026-11-20): test the read-counting rule; bar accepted as proposed. This is arm K of the draft
`docs/PREREG_sirv_rules_2026-11-18.md` (row 1801, never run on whole datasets), taken alone at today's defaults. The draft's
other arms (spliced share, shadow support) are not part of this row; its extension arm is moot (extension is off since row 1812).

**The lever.** By default primary reads with the same contig, start, end and intron chain count once. It was kept as a default
for primate recall (`docs/archive/2026-09/PREREG_primary_dedupe_2026-09-22.md`: gorilla -1,096 matching chains at +3.2 points
of precision without it). On the SIRV contigs of the mouse ES library it is the main cause of wrong chains (71 -> 17 with
`--keep-coordinate-duplicates`, diagnosis of 2026-11-18 at the defaults of that day): deep loci have many identical reads, so the
0.02 isoform fraction is computed on collapsed counts and stops removing errors.

**Arms.** `base` (today's defaults) and `kd` (the driver with `--keep-coordinate-duplicates`); stage `assemble`, `RUSTLE_PTC=0`.
Workspace `r1833_count_all_reads`. Registered before either arm ran.

**Data.** Bars on U1 zebrafish simulated and U2 mouse heart (ENCSR938RPH), as the draft fixed: neither was used to diagnose or
choose this lever (both served other rows as check data). Shown without a bar: fly and mouse simulated, mouse ES (mouse contigs,
and its SIRV contigs apart), OR6737, A119b. Zebrafish runs on the re-mapped BAM of row 1830.

**Unit.** Exact intron chains (`bench/famcmp/curve.py`), `PREFIX.gtf`; `PREFIX.final.gtf` shown. Gained and lost chains
separately (correct and wrong), never the net alone.

**Bar D3, per dataset:** `kd` against `base` on `PREFIX.gtf`: precision higher, and at most 5% of base's correct chains lost
(row 1793's "pays").
**Decision.** D3 on both U1 and U2 = candidate default; the user rules, and the known primate cost is put beside it. Either
fails = reported, nothing proposed.

**Limits, stated in advance.** U1 is simulated: IsoSeqSim draws read ends from completeness tables, so identical coordinates
may be rarer or commoner than in a real library, and the result there says little about PCR or deep-locus duplicates. U2 shares
the species of the library the diagnosis came from. The chain score filter and the 5' truncation rule were fitted on collapsed
counts, so `final.gtf` under `kd` is outside their fit. The `reads` attribute changes meaning under `kd` (molecules instead of
distinct coordinates), which later stages read. Two primate samples stand for twelve.

## Result (2026-11-20): D3 PASSES on both check datasets; A119b gets worse on both axes

Jobs 4484224-4484230, both arms exit 0 everywhere; `r1833_count_all_reads/score.NAME.txt`. Repo `d83140b4`, binaries of
2026-11-20. `PREFIX.gtf`, `kd` against `base`:

| dataset | precision | correct chains (net) | correct gained / lost | wrong gained / lost | |
|---|---|---|---|---|---|
| **U1 zebrafish simulated** | 0.813 -> 0.834 | 12,574 -> 12,542 (-0.3%) | +229 / -261 | +93 / -499 | **D3 PASS** |
| **U2 mouse heart** | 0.546 -> 0.568 | 13,234 -> 13,093 (-1.1%) | +43 / -184 | +200 / -1,244 | **D3 PASS** |
| fly simulated | 0.926 -> 0.939 | 7,082 -> 6,913 (-2.4%) | +73 / -242 | +13 / -130 | shown |
| mouse simulated | 0.903 -> 0.925 | 27,731 -> 25,963 (-6.4%) | +320 / -2,088 | +77 / -975 | shown (over the 5% bound) |
| mouse ES, mouse contigs | 0.331 -> 0.372 | 18,058 -> 17,378 (-3.8%) | +24 / -704 | +379 / -7,481 | shown |
| SIRV contigs of mouse ES | 0.589 -> 0.862 | 53 -> 50 | 0 / -3 | +1 / -30 | shown (the diagnosis data) |
| OR6737 | 0.321 -> 0.346 | 25,792 -> 24,695 (-4.3%) | +142 / -1,239 | +1,883 / -9,638 | shown |
| A119b | 0.158 -> 0.139 | 39,049 -> 36,277 (-7.1%) | +874 / -3,646 | +59,290 / -42,185 | shown |

`PREFIX.final.gtf` moves the same way: zebrafish 0.891 -> 0.908 (12,381 -> 12,393), heart 0.577 -> 0.600 (-0.9%), fly 0.964 ->
0.971 (-2.4%), mouse simulated 0.953 -> 0.968 (-5.8%), mouse ES 0.347 -> 0.392 (-3.3%), SIRV 48 / 0.632 -> 49 / 0.875, OR6737
0.377 -> 0.408 (-4.1%), A119b 0.219 -> 0.181 (-6.8%).

**Verdict by the registered rule:** candidate (D3 on U1 and U2). **What stands beside it:**
- It raises precision on seven of eight views and removes 20% of the wrong chains on mouse ES (7,102 net) and 14% on OR6737.
- It is over the 5% loss bound on mouse simulated, and on A119b it is worse on both axes: 7.1% fewer correct chains AND 17,105
  more wrong chains. The earlier "primate recall" reading is not what this shows: OR6737 behaves like the mouse libraries;
  A119b alone, the deepest library, floods (counting every read lifts many one-coordinate chains over the 2-read floor).
- So the collapse is not a primate fit; it protects one deep library. A default for every library would trade A119b for the rest.
Not examined: why A119b floods (PCR duplicates of that library, or depth alone), and whether a read floor that scales with
depth would let `kd` pass there. The driver has no switch for the flag (arms patched a driver copy).
