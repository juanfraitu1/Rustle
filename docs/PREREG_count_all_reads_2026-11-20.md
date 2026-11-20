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
