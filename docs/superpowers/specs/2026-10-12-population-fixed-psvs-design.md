# Population-fixed PSVs for O2 — design (2026-10-12)

## Intent (agreed with the user)

O2 should assign at least some reads from a man whose genome is NOT the reference (A119b testis Iso-Seq vs CHM13's
HG002 Y) and say how far each call can be trusted. A matched genome + reads from the same man stays the ideal; this is
the fallback. Motivation: `docs/archive/2026-10/merged/PREREG_yag_o2_2026-10-12.md` — on the YAGs O2 is calibrated but split-half concordance on
real testis reads was 0.544, plausibly because the donor's copies differ from the reference man's (gene conversion
between palindrome arms varies between men). Scope chosen by the user: sub-project **B (population-fixed PSVs) first**,
**human only**, population = **Hallast et al. 2023's 43 Y assemblies**, usage = **both** a portability label and an
opt-in strict mode, fixed-ness by a **hybrid** k-mer + alignment check. Sub-project A (assign to copy groups under a
flag) is separate and later.

## 1. Data

- The 48 per-man chrY assemblies `verkko_1.0_release/chrY/*.chrY.fasta` of
  `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/working/20230412_sigY_assembly_data/` (2.78 GB),
  each checked against `MANIFEST_20230412_sigY_assembly_data.tsv` (md5).
- Stored at `/mnt/c/Users/jfris/Desktop/data/hgsvc_sigY/` (C: has room; `/mnt/linuxdisk` has ~5.6 GB free), with
  `SOURCE.md` (URL, date, manifest, sample ids). Not downloaded: whole-genome, chrX and region-subset files.

## 2. Population status of each PSV site (`bench/popfix/`)

- **PSV sites:** the decisive columns of every catalog family with copies on chrY, from `utilities psv-discover`
  (default engine), which gives each copy's genome position and base per column. A column involving a copy not on chrY
  is `NOT_ASSESSED` (Y assemblies only).
- **k-mer layer (every column):** man m *carries* copy c's allele iff his assembly contains at least one of the 31
  31-mers spanning c's position in c's reference sequence with that allele; m *covers* the column iff his assembly
  contains c's flanking 31-mers starting 40 bp left and right of the position (both sides). Column `FIXED_K` iff at least
  half the men cover it and every covering man carries every distinct allele of the column; `POLYMORPHIC_K` if some
  covering man lacks an allele. Counting: the existing `kscan` (candidate k-mers only) per assembly.
- **Alignment layer (unambiguous copies only):** a copy qualifies if its locus +- 10 kb flanking 31-mers are unique in
  CHM13. The locus + flanks is aligned to each man's Y (`minimap2 -cx asm5 --cs`); his ortholog = the single alignment
  covering both flanks; his base at each PSV position is read through the CIGAR. Column `FIXED_A` iff every man with an
  ortholog for each of its qualifying copies has the reference base there and at least half the men have one;
  `POLYMORPHIC_A` otherwise. Where a column has an alignment verdict it overrides the k-mer verdict.
- **Statuses:** `FIXED` (A or K), `POLYMORPHIC`, `UNCOVERED` (< half the men cover / have an ortholog), `NOT_ASSESSED`.
- **Output** `psv_population.tsv`: one row per (chrom, position) site: family, status, method, n_men, n_cover,
  n_carry, n_mismatch.

## 3. O2 integration (`src/family.rs`, `src/copy_assign_bin.rs`)

- The read-star stream already computes, per candidate, the target offset of every read position (`per_pos`); it is
  turned into a genome position per (candidate, column) — genomic form: locus start + offset; unit form: `exon_map`.
  Kept in `ReadStarObs` and in the star dump.
- A column counts as fixed for the pair (best, rival) iff both copies' positions are `FIXED` in the table.
- `--psv-population FILE`: load the table; write `<out>.population_tier.tsv` per contested molecule: normal status /
  copy, the pairwise certificate re-run on the fixed columns only (status / copy), tier `PORTABLE` (fixed-only assigns the
  same copy) | `REFERENCE_ONLY` (assigned only with non-fixed columns) | `UNASSIGNED`. A separate file: no change to
  `assignments.tsv` or its readers.
- `--psv-fixed-only` (requires the table): the fixed-only certificate IS the decision (`assignments.tsv` = portable
  calls only).
- Without the flags every output is byte-identical to today.
- Tests: position mapping (both forms), masking, tier logic, flag validation.

## 4. Validation (pre-registered before any run: `docs/archive/2026-10/merged/PREREG_popfix_2026-10-12.md`)

- (a) Table sanity (descriptive): status counts by YAG family and region class; expectation: palindrome-arm PSVs often
  `POLYMORPHIC`, unique-flank PSVs mostly `FIXED`; k-mer vs alignment agreement where both apply.
- (b) Simulation (the 10-06 truth-labelled chrY reads) with `--psv-fixed-only`: 0 wrong among assigned and 0 of the
  undecidable assigned; portable count reported.
- (c) Real A119b testis, YAG families: split-half concordance on fixed PSVs only >= 0.98 with >= 30 molecules whose two
  halves both assign — the test of the donor-mismatch explanation (if it stays low, the disagreement is O2's evidence,
  not the donor); molecules with a `PORTABLE` single-copy assignment reported.
- (d) No regression: outputs byte-identical without the flags on the dev subset; full suite passes (3 known
  `pipeline_driver_dispatch` failures excepted).

## Order

Download + table (bench, no Rust) -> prereg -> Rust (TDD) -> runs (a)-(d).

## Not in scope

Gorilla (no population Y assemblies), autosomes (would need HPRC/HGSVC whole genomes), copy-group assignment
(sub-project A).
