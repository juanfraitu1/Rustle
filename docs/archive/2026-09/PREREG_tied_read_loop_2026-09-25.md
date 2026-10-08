# Pre-registration: the closed loop (families -> assign tied reads -> give each read to its copy -> re-assemble)

**Written 2026-09-25 16:20, before any pass-2 assembly, any home table or any real-read union-certificate count
exists.** User decision of 2026-09-25 16:00 (recorded in the figure-phase notes, item 3): "Build the CLOSED LOOP,
pre-registered: families -> assign tied multimappers (cross-family test) -> give each read to its assigned copy ->
re-assemble -> new loci that help assembly; judged like StringTie/FLAIR de novo and in the Liftoff framework; verdict
on held-out samples." Figure 9 (`figures/fig_loop.py`, `figures/captions/fig9.md`) shows the result.

## 0. What existed and what was seen before this file

- **No pass-2 assembly has ever been run.** Until this work no assembler pass consumed an assignment:
  `copy_assign`'s own assembly is built from the pre-assignment read pool (`src/bin/copy_assign.rs`, the
  `--gtf` comment "independent of the assignment"), and the driver's `assemble` stage runs before `assign`
  (`tools/rustle_pipeline.sh`). Register r847 (§6m5, "O2 contributed 0 junctions") was therefore a null by
  construction, not a measurement of redistribution. Its NPIP ceiling (perfect placement of the existing records:
  85 -> 91 of 106 junctions, 2/9 -> 4/9 complete chains, `bench/LOCUS_ASSEMBLY_NPIP.md`) is the only prior bound.
- **Seen before this file (disclosed; these runs used the legacy catalog, not the default families):** the union
  certificate on the Fig. 4 development simulations assigns **0** molecules: human chr16 993 molecules scored, 991
  left unassigned (`tied`), 2 `ambiguous`; gorilla NC_073244.2 30 `tied`
  (`${work}/o2sim/{human,gorilla}/u2.union_certificate.tsv`; r1103, r1109: 1,255 of 1,263 MAPQ-0 human reads have an
  NM-identical twin). **No union-certificate run on real reads exists for any sample** (checked with `find` over
  `${work}` and `rustle_figures_dev`: every `union_certificate.tsv` is a simulation).
- Closest refuted neighbours (checked in `docs/NEGATIVE_RESULTS_REGISTER.md`), and how this differs:
  - r1102 (§6zl, SCOPED within-family secondary pool): complete members 58 = GOOD = ALL on NC_073244.2. SCOPED
    **added** a read to every member it touched; the loop gives an assigned read to **one** copy and **removes** its
    other records, the primary included. Nothing in the register tests that subtraction.
  - r1059 (all secondaries seed loci): the loop never adds an unassigned echo.
  - r338 / r435 / r980 (`--tied-seed`): the loop seeds nothing without a certificate verdict.
  - r579 / r1089 (EM, 1/k): **no EM, no fractional reads.** A read is moved only on an `assigned` verdict; every other
    read keeps a fixed rule (the base pool of its arm).
  - r395 (bounding loci by family siblings, "a structural feedback loop cements the error"): **families are frozen
    from pass 1** (§1.1); pass-2 loci never re-derive the families used in this claim.
  - r1092 / r1093 (the per-family table has no cross-family arbitration; 904 foreign `assigned` rows): **only the
    union certificate's verdicts are used** (§1.2); the per-family `assignments.tsv` is never read.
  - r1101 / r1026 (family-level losses are mostly missing edges): the loop can change loci and chains; no family-level
    claim is made here.
- Predicted before looking: ⛔ or a small effect (r1102; the union assigns 0 on both development simulations; tied
  molecules are 2.5% of human and 1.1% of gorilla reads at the 0.98 rule, `docs/PREREG_locus_read_pool_2026-09-22.md`).

## 1. The loop

### 1.1 Pass 1 and the families (frozen)

- **Pass-1 assembly** = the shipped `assemble` stage (`tools/rustle_pipeline.sh assemble`): `copy_assign
  --assemble-only --genome-wide --assembly-junctions strict` + the shipped polish, loci seeded from primaries plus
  secondaries with AS >= 0.98 x the molecule's genome-wide best AS (`as_table`). Called **GOOD** below.
- **Families** = the default de novo family definition (user decision 2026-09-25 16:00, item 1): the `families` stage
  on the pass-1 GTF (`mcl_families --from-gtf <pass1.gtf> --min-exonic-bp 1 --min-shared-exon-frac 0.60`, MCL 2.8),
  with `--bam` and `--fasta` so that it writes the per-locus **units** in the `copy_assign --families/--copies-fa`
  contract (`<prefix>.fam.units.tsv/.fa`). Copy assignment consumes these units (the same families). The legacy
  `gw_family_catalog` copies may be used only for the development mechanics of §5, never for a verdict number.
- The families and units are computed once, from pass 1, and are not recomputed from any pass-2 output.

### 1.2 Assignment (union verdicts only)

`copy_assign --bam B --fasta G --regions <every contig whole> --families UNITS.tsv --copies-fa UNITS.fa
--union-certificate --out <prefix>.loop` (current defaults: the AS-tied gate at an **exact** tie,
`--as-tie-ratio 1.0`; placements on other contigs are invisible to it and such molecules are "left as today").
Only `<prefix>.loop.union_certificate.tsv` is read:
- `assigned_family`: `winner` = `family_id:copy_idx`, resolved to that unit's row (`chrom`, `start`, `end`);
- `assigned_outside`: `winner` = `outside:chrom:start1-end`, the pseudo-copy over the placement's aligned blocks
  (1-based start in the label; converted to 0-based half-open);
- `tied`, `ambiguous`, `no_result`: not moved.

Tie widths differ and are named everywhere: assignment uses the exact tie (1.0); pass-1 seeding uses 0.98.

### 1.3 The home table

`bench/loop_home.py union --union U.tsv --families UNITS.tsv --bam B --out HOME.tsv`: one row per assigned molecule,
`read_name, chrom, start, end` (0-based half-open home span) + provenance columns. The builder checks, in the BAM, that
the molecule has at least one mapped non-supplementary record with an aligned block (M/=/X) overlapping its home span.
A molecule with **no** record at home is **not** put in the table (it keeps its arm's base rule) and is counted as
`H` (its placement would need a lifted record; deferred, never faked).

### 1.4 The pass-2 filter (`RUSTLE_READ_HOME_TABLE=HOME.tsv`, opt-in)

For a molecule in the table, a mapped non-supplementary record enters the assembler's read pool **iff** one of its
aligned blocks overlaps the home span on the same contig by >= 1 bp, whether it is primary or secondary and whatever
its AS. Every other record of that molecule is dropped, **the primary included**. Molecules not in the table follow
the arm's base rule unchanged. The filter touches only the assembler's read pool (pass-1 skeletons); it never changes
the reads a copy assignment sees. Unset, every output is byte-identical (proved with `cmp`, §5).

### 1.5 Arms (same binary, BAM, genome, polish, junction mode and best-AS table in every arm)

| arm | base rule for molecules not in the table | table |
|---|---|---|
| **P** | primaries only (`--no-seed-secondaries`) | none |
| **GOOD** | primaries + secondaries >= 0.98 x genome-wide best AS (the shipped default) | none |
| **LOOP-S** | GOOD | union home table (§1.3) |
| **LOOP-P** | P | union home table |
| **ORACLE** (simulation only) | GOOD | every simulated molecule's home = its true source copy's span |

LOOP-S is the headline: it cannot lose GOOD's transcripts built from abstained reads, because an abstained read keeps
GOOD's placements. LOOP-P is the "redistribute only" reading (descriptive).

## 2. Gate G0 (a stop rule, not a result)

Per verdict sample, from the LOOP-S pass-2 run's own counters (records, not scores; printed by the filter):
- `A` = molecules assigned by the union (family + outside); `H` = assigned molecules with no record at home;
- `R_drop` = records GOOD admits that the filter removed (a record GOOD would not admit anyway is counted
  separately and does not count); `R_add` = records the filter admitted at home that GOOD would not admit.
  Both are counted before the assembler's coordinate de-duplication.

Each record lands in at most one pass-1 read group, so at most `R_drop + R_add` pass-1 read groups differ between GOOD
and LOOP-S. **If `R_drop + R_add < 20` on a sample, that sample is ⛔ "bounded" (no scoring); if every verdict sample is
bounded, the loop is ⛔ as a bounded negative with these counts.**

## 3. Metrics (per sample and arm; species never pooled; every contig the sample's annotation covers)

- **M1** intron-chain sensitivity and precision against the sample's RefSeq annotation, gffcompare 0.12.10, exactly as
  Fig. 1 (`figures/assembly.py`).
- **M2** paired comparison on the same reference chains, LOOP-S vs GOOD (and LOOP-P vs P): chains with >= 2
  exact-chain reads (the Fig. 1 paired rule), Tango 95% interval, exact McNemar p.
- **M3** Fig. 3's tie bins (0, (0, 0.1], (0.1, 0.5], (0.5, 0.9], > 0.9; tie = second AS >= 0.98 x best, genome-wide)
  and the no-primary facet; matched chains per arm; read-sharing groups (Fig. 3) are the unit; every changed chain is
  listed with its read-sharing group, and a CGB-like / RFPL4A-like array is named when it carries a change (r1110,
  r1116).
- **M4** loci in the Liftoff framework (`docs/PREREG_liftoff_loci_2026-09-25.md` §3-§4, C2): the fraction of
  read-supported Liftoff reference loci found by an arm's de novo locus (cov >= 0.5 of the reference's exon bases), for
  annotated loci in place, moved, and extra copies (sequence_ID >= 0.95 / 0.98 / 0.99 / 1.00), and the fraction of the
  arm's de novo loci at a Liftoff locus. **Two fixed universes**, the same for every arm:
  **U1** = C2's (>= 2 reads whose primary alignment, `-F 2308`, has an aligned block on the exon union);
  **U2** = >= 2 reads with a candidate placement (the primary, or a secondary with AS >= 0.98 x the molecule's
  genome-wide best AS from `as_table`) with an aligned block on the exon union. U1 excludes exactly the copies whose
  reads are tied secondaries, which the loop targets; U2 counts them.
- **M5** (simulation only) per simulated source copy (>= 1 simulated read): multi-exon copies, an arm transcript with
  the copy's exact intron chain (gffcompare `=` against a GTF of the simulated copies); single-exon copies, an arm
  locus covering >= 0.5 of the copy's exon bases. Sensitivity over copies; precision = arm transcripts matching some
  simulated copy / arm transcripts overlapping a simulated copy's span.

Intervals: M1-M3 as in Figs 1 and 3 (chains; read-sharing groups or genes resampled); M4 as in
`PREREG_liftoff_loci` (source records resampled, seed 20260925); M5 resampling source copies.

## 4. Claims and bars

| id | claim | bar | if it fails |
|---|---|---|---|
| LOOP1 | giving union-assigned reads to their copy and re-assembling helps de novo assembly where reads tie | per verdict sample PASS iff (a) M3: LOOP-S minus GOOD, net matched reference chains in the tie > 0.5 bins > 0, with the gained chains in >= 2 distinct read-sharing groups; (b) M1 precision of LOOP-S >= GOOD - 0.2 points; (c) M4 fraction of de novo loci at a Liftoff locus of LOOP-S >= GOOD - 0.2 points. ⭐ if >= 3 of the 4 verdict samples PASS; ⛔ if <= 1 PASS or G0 bounds every sample; ⚠ (mixed) at 2 | reported as ⛔ with the counts; no rule is changed after a verdict number |
| LOOP2 | a perfect placement of the reads (ORACLE) would change assembly | on the simulation of each verdict sample (the Fig. 4 genome-wide simulation of that sample, when built), M5 exact-chain sensitivity ORACLE - GOOD >= +1.0 point **or** >= 5 more copies recovered | "the placement ceiling is closed": no assignment rule can help assembly on that simulation |
| LOOP3 | (descriptive) LOOP-P vs P and vs GOOD; Liftoff D1 on U1 and U2 for every arm; A, H, R_drop, R_add per sample | none | - |

Samples: **verdict** = human_testis, gorilla_KB3781, chimp_PTR, orangutan_PPY (never used to develop the assembler,
the seeding rule or copy assignment; KB3781 was used for O3 only, §6ze, which is disclosed). **Development** = human
A119b (copy-assignment rules on chr16) and gorilla OR6737 (seeding verdict contig NC_073244.2): reported, never part
of a verdict. Species and samples are never pooled.

## 5. Development and validation before the verdict (mechanics only)

- **V-byte.** With `RUSTLE_READ_HOME_TABLE` unset, the pass-2 binary's GTF is `cmp`-identical to the pre-change
  binary's on a slice, on the streaming path (GOOD with a best-AS table) and on the buffered path (GOOD without a
  table, and P). With an EMPTY table file the filter is inert (same GTF).
- **V-paths.** On the development simulation (human chr16, `${work}/o2sim/human`), the ORACLE table gives the same
  GTF through the streaming path (best-AS table loaded) and the buffered path (`--materialize-reads`).
- **Development run D1:** every arm on the human chr16 development simulation (union verdicts from the existing
  `u2.union_certificate.tsv`, legacy catalog; ORACLE from the simulation's truth). Its numbers go into the amendments
  below with the date; they are development numbers and cannot change §1-§4.
- Unit tests of the home-table parser and the keep/drop rule (`cargo test`).

## 6. Stop rules and cost

- Every call <= 10 min and <= 20 GB under `/mnt/linuxdisk/tmp/rustle_heavy.lock`. The genome-wide union assignment is
  the expensive step (the `assign` stage: est. 6 min for chimp_PTR to 5 h for human_A119b in one process); if it does
  not fit, it is sharded by families as in `docs/PREREG_genome_wide_copy_assignment_2026-09-25.md` (validated shards),
  and if one shard does not fit it stops for the cluster. Pass-2 assembly costs about one pass-1 assembly.
- If the default families' units (§1.1) do not exist for a verdict sample, its verdict waits; it is never computed on
  the legacy catalog.
- Numbers are new register rows; r847 and r1102 are not overwritten.

## Amendments

**Amendment 1 (2026-09-25 17:15): V-byte, V-paths and development run D1 (human chr16 simulation; legacy catalog).
Development numbers; §1-§4 unchanged.**
- *V-byte (PASS).* On the development simulation BAM (`${work}/o2sim/human/sim.bam`, 164,300 records, 135,847
  secondary), the pre-change and post-change `copy_assign` (built from the same tree, differing only by the filter)
  give `cmp`-identical products (GTF, assignments, families, famcn, params, quant) on four paths: GOOD streaming (best-AS
  table), GOOD buffered (`--materialize-reads`), P streaming, P buffered; with an EMPTY table (header only) GOOD
  streaming and GOOD buffered are also identical to the unset runs. The same holds on the real human A119b chr16
  slice (`chr16_arm/chr16.bam`, 1,769,085 records; best-AS table of the slice): GOOD streaming and P, pre vs post, and
  GOOD with an empty table, every product `cmp`-identical (11,435 transcripts; 15-19 s, 0.3 GB per run).
- *V-paths (PASS).* ORACLE through the streaming path and through the buffered path: identical GTF apart from the
  `matched_reads` attribute (the known streaming difference, r1086), and identical filter counters (kept at home
  28,446, R_add 1,518, R_drop 8,554, away and outside the pool 125,680).
- *D1 (M5, `rustle_figures_dev/loop/d1_sim_human/m5.tsv`).* Union verdicts on this simulation: 991 `tied`, 2
  `ambiguous`, 0 assigned (as disclosed in §0), so LOOP-S and LOOP-P are byte-identical to GOOD and P. ORACLE home table:
  28,453 simulated molecules, 2 with no record at home (left out). Multi-exon copies with the exact chain (of 217):
  P 191, GOOD 188, ORACLE 191. Single-exon copies covered (of 1,181): P 226, GOOD 247, ORACLE 209. Transcripts
  matching a copy / at a copy: P 401/435, GOOD 393/493, ORACLE 386/409. Read from the logs: ORACLE's single-exon loss
  is the polish's data-adaptive single-exon read floor (28 reads on chr16) cutting copies whose own 10-23 simulated reads
  no longer carry the echo reads of their siblings. Reported, not tuned: no rule changes (the polish is the shipped
  default and part of every arm).
