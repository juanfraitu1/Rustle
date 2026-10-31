# Multi-copy assembly mode (FAMASM) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A Python prototype that pools a family's reads across copies, defines chains family-wide, and emits a chain at a copy only where O2-certified reads support it, layered on the unchanged base GTF; scored on dev A119b and held-out HG002 under the spec's bars M1-M4.

**Architecture:** One script `bench/stgaps/famasm.py` with five subcommands (`pool`, `lift`, `admit`, `gtf`, `score`), each writing an inspectable TSV; the copy map and chain parsing come from `bench/stgaps/fw_layer.py`; certificates come from the shipped assign stage (`tools/rustle_pipeline.sh assign`, product `PREFIX.assignments.tsv`, status `assigned|ambiguous|tied`); scoring reuses `stratify.py`, `observed.py`, `mainstay_score.py`, `pairscore.py` and a `mcl_families --from-gtf` rerun on the FAMASM GTF.

**Tech Stack:** Python 3 + pysam (pixi env `~/.pixi/detached-envs/rustle-7434162673922408322/envs/default/bin`), minimap2, samtools, gffcompare (to install), the Rust binaries in `~/rustle_target/release/`, Slurm `singlenode`.

**Spec:** `docs/superpowers/specs/2026-10-31-multicopy-assembly-design.md`

## Global Constraints
- BASE GTF lines are never altered or removed; FAMASM = BASE + admitted chains absent from BASE at that copy (same contig, exact intron chain).
- Admission needs reads CERTIFIED (`status == assigned`, `assigned_copy == c`) to copy c: >= 2 whole-chain reads (`via whole`) or compatible sub-chains covering every junction with >= 2 reads each (`via compat`). `ambiguous` / `tied` / uncertified reads count for no copy.
- Pool = primary records plus secondary records with `AS >= 0.98 x` the molecule's genome-wide best AS (`PREFIX.molecules.tsv`, columns `name best_as second_as ...`); one molecule once per family.
- Family chain floor = 2 molecules (pass-1 rule); compatible support is used in admission only, never for chain discovery.
- Single-copy loci and one-copy families: byte-identical to BASE (checked by `diff` on the GTF restricted to those loci).
- Bars (spec): M1 multi-copy recall +1.0 point; M2 multi-copy precision >= BASE - 0.5; M3 single-copy |delta| <= 0.1 and every BASE line kept; M4 mainstays recalled >= BASE, unrelated foreign <= BASE, pair precision and recall >= BASE - 0.01. Dev A119b reported; held-out HG002 Kinnex decides. Prereg row in `docs/NEGATIVE_RESULTS_REGISTER.md` and section in `docs/LEDGER_O1.md` BEFORE any dev number (Task 7).
- Every file under `~/work/jfim/runs/famasm_{A119b,HG002}/`; commits end with the session attribution lines; `git fetch -q && git rebase -q origin/main` before push.

## Review Focus
1. A read certified to copy c whose certified record's chain has a junction that does not lift: the chain stays copy-private and may only be admitted at c (Task 3 test `private_junction_admits_only_at_own_copy`).
2. A family whose representative copy has no PAF hit for one member copy (lift table missing): that copy's reads are private, nothing crashes, the copy is reported in `lift.unmapped.tsv` (Task 3 test `missing_lift_is_private`).
3. Two admitted chains that invert to the same genomic chain at c (two family chains differing only in private junctions): emitted once (Task 4 test `inverse_lift_dedups`).
4. A read with a primary on copy A but certified to copy B: its chain is taken from the record on B when one exists in the pool, else the read is skipped and counted in `admit.stats.tsv` as `certified_without_record` (Task 4 test `chain_comes_from_certified_record`).
5. BASE already has the chain at c (exact intron chain): not appended, counted `already_in_base` (Task 5 test `no_duplicate_of_base`).

---

### Task 1: Scoring assets and tools (no product of the mode yet)

**Files:**
- Create: `~/work/jfim/runs/famasm_A119b/ref/` and `~/work/jfim/runs/famasm_HG002/ref/` (reference transcripts, PAF, strata)
- Modify: `~/work/jfim/tools/isocon/pixi.toml` (add `gffcompare`, channel `bioconda`) or a new `~/work/jfim/tools/gff/` pixi env if the isocon pins conflict

**Interfaces:**
- Produces: `REF.gtf` (human `A119b.chm13v2.0_RefSeq.gtf.gz` as is; HG002 `gffread ~/data/hg002/genes.gff3 -T -o ref/hg002.gtf`), `reftx.tsv` + `reftx.fa` (`bench/stgaps/ref_tx_fasta.py REF.gtf GENOME reftx.fa reftx.tsv`), `reftx.paf` (`minimap2 -x splice -N 50 --secondary=yes -c GENOME reftx.fa`), `strata.tsv` from `bench/stgaps/stratify.py` (its OUTDIR product naming the stratum of each reference transcript; read the script for the exact file name and record it in `ref/README`).

- [ ] **Step 1:** `cd ~/work/jfim/tools/isocon && pixi add -c bioconda gffcompare`; run `pixi run gffcompare --version`. Expected: a version line. If the solve fails, create `~/work/jfim/tools/gff` with `pixi init` + `pixi add -c bioconda -c conda-forge gffcompare gffread` instead.
- [ ] **Step 2:** Build the four reference products per sample with the commands in Interfaces (Slurm job `ref.sbatch` in each run dir, 4 h, 16 threads); verify `wc -l reftx.tsv` > 100,000 (human) and > 50,000 (HG002) and the PAF is non-empty.
- [ ] **Step 3:** Run `stratify.py reftx.tsv reftx.paf REF.gtf ref/strata base=PREFIX.families.gtf` for both samples; record multi-copy and single-copy recall / precision of BASE in `ref/base_strata.txt`. These BASE numbers are the comparators of M1-M3.
- [ ] **Step 4:** Commit nothing in the repo (assets live under `~/work`); write `ref/README` listing every command used.

### Task 2: `famasm.py pool`

**Files:**
- Create: `bench/stgaps/famasm.py`
- Test: `bench/stgaps/test_famasm.py`

**Interfaces:**
- Consumes: `fw_layer.read_copies(path) -> {(fam, idx): (chrom, start, end, n_reads)}`, `fw_layer.chain_of(rec) -> tuple[(intron_start, intron_end)]` (0-based genomic), `PREFIX.molecules.tsv`.
- Produces: `pool(copies_path, bam_path, molecules_path, out_dir) -> None` writing `out_dir/pool/<fam>.tsv` with header `read\tcopy_idx\tis_primary\tAS\tchain` (chain via `fw_layer.fmt`, empty for an unspliced record) and `out_dir/pool.stats.tsv` (`family\tmolecules\trecords\tsecondary_kept`). Helper `good_secondary(as_: int, best_as: int) -> bool` = `as_ >= 0.98 * best_as`.

- [ ] **Step 1: Failing tests** in `test_famasm.py`: `test_good_secondary_threshold` (`good_secondary(98, 100)` True, `good_secondary(97, 100)` False); `test_pool_keeps_one_record_per_copy_and_drops_bad_secondaries` on a tiny BAM built in the test with pysam (two copies on one contig, read r1 primary on copy 1 AS 100 and secondary on copy 2 AS 99 -> both rows; read r2 secondary on copy 2 AS 80 -> dropped; a supplementary record -> dropped).
- [ ] **Step 2:** Run `python3 -m pytest bench/stgaps/test_famasm.py -k "good_secondary or pool" -v`. Expected: FAIL (module missing).
- [ ] **Step 3:** Implement `pool` and `good_secondary` in `famasm.py`; regions per copy come from `read_copies`; use `samtools`-free pysam `fetch(chrom, start, end)` with `-F 0x904` semantics (`is_unmapped`, `is_supplementary`, `is_qcfail` skipped; secondaries kept when `good_secondary`). A record is assigned to the copy whose locus it overlaps most.
- [ ] **Step 4:** Run the tests. Expected: PASS.
- [ ] **Step 5:** Commit: `git add bench/stgaps/famasm.py bench/stgaps/test_famasm.py && git commit -m "famasm pool: family read pools with good secondaries"`.

### Task 3: `famasm.py lift`

**Files:**
- Modify: `bench/stgaps/famasm.py`
- Test: `bench/stgaps/test_famasm.py`

**Interfaces:**
- Consumes: `fw_layer.map_(copies_path, fa_path, out)` (writes `reps.fa`, `copies.fa`), the PAF `minimap2 -x asm20 -c --cs copies.fa reps.fa`, `fw_layer.load_lift(paf)`, `fw_layer.fam_chain(copies, reps, lifts, fam, idx, chain) -> tuple | None`.
- Produces: `lift(copies_path, fa_path, out_dir) -> None` writing `out_dir/lift/<fam>.tsv` with header `read\tcopy_idx\tis_primary\tfam_chain\tprivate` where `fam_chain` is `fw_layer.fmt` of the lifted chain, and `private` is `1` when any junction failed to lift (then `fam_chain` carries the genomic chain prefixed `c<idx>:`); `out_dir/lift.unmapped.tsv` (`family\tcopy_idx`) for copies with no lift table. Also `inverse_lift(lifts, copies, reps, fam, idx, fam_chain) -> tuple | None` (a family junction back to copy idx's genomic coordinates; `None` when any junction has no aligned column).

- [ ] **Step 1: Failing tests**: `test_private_junction_admits_only_at_own_copy` (a chain with one unliftable junction -> `private == 1`, `fam_chain` starts with `c`), `test_missing_lift_is_private` (copy absent from the lift dict -> every read private, row in `lift.unmapped.tsv`), `test_inverse_lift_roundtrip` (a representative chain lifted then inverted at a copy with a synthetic `+` and a `-` strand lift block returns the original genomic chain). Build lift dicts directly in the `load_lift` return shape `{(fam, idx): (q_starts, blocks)}`.
- [ ] **Step 2:** Run `pytest -k "private or missing_lift or inverse"`. Expected: FAIL.
- [ ] **Step 3:** Implement `lift` (calls `map_`, runs minimap2 via `subprocess`, `load_lift`, applies `fam_chain` to every `pool` row) and `inverse_lift` (search `blocks` by target coordinate: for a block `(q_lo, q_hi, t0, s)` the target range is `[t0, t0 + (q_hi - q_lo))`; invert the `+`/`-` formulas of `fw_layer.lift`).
- [ ] **Step 4:** Run the tests. Expected: PASS.
- [ ] **Step 5:** Commit: `"famasm lift: family junctions and inverse lift"`.

### Task 4: `famasm.py admit`

**Files:**
- Modify: `bench/stgaps/famasm.py`
- Test: `bench/stgaps/test_famasm.py`

**Interfaces:**
- Consumes: `lift/<fam>.tsv`, `PREFIX.assignments.tsv` (columns `read_name family_id assigned_copy status`; certified = `status == "assigned"`), `inverse_lift`.
- Produces: `admit(out_dir, assignments_path, min_reads=2) -> None` writing `out_dir/admit/<fam>.tsv` with header `copy_idx\tfam_chain\tgenomic_chain\tvia\tcertified_reads` (`via` in `whole|compat`) and `out_dir/admit.stats.tsv` (`family\tfamily_chains\tadmitted_whole\tadmitted_compat\tcertified_without_record\tuninvertible`). Pure functions: `family_chains(rows, min_reads=2) -> set[chain]` (exact family-chain groups with >= min_reads distinct molecules over ALL copies, private chains excluded); `admit_at_copy(chain, certified_rows_at_c, min_reads=2) -> str | None` returning `"whole"`, `"compat"` or `None`, where compat = every junction of `chain` is covered by >= min_reads certified reads whose own chain is a contiguous sub-chain of `chain` (its junctions are a consecutive run of `chain`'s junctions) and whose record's first/last exon does not extend past the sub-chain's flanking junctions into an intron of `chain` (end-compatible).

- [ ] **Step 1: Failing tests**: `test_family_chain_floor_counts_molecules_over_copies` (one read at copy 1 and one at copy 2 with the same family chain -> in `family_chains`; a chain with one read -> not), `test_admit_whole_needs_two_certified` (two `assigned` reads at c with the whole chain -> `whole`; one assigned + one `tied` -> `None`), `test_admit_compat_covers_every_junction_twice` (chain J1-J2-J3; reads with J1-J2 x2 and J2-J3 x2 -> `compat`; J1-J2 x2 and J3 x1 -> `None`), `test_chain_comes_from_certified_record` (a read certified to copy 2 with records on copies 1 and 2 uses the copy-2 chain; with a record on copy 1 only -> counted `certified_without_record`), `test_inverse_lift_dedups` (two family chains inverting to the same genomic chain at c -> one output row, `certified_reads` summed).
- [ ] **Step 2:** Run `pytest -k "family_chain or admit or certified_record or dedups"`. Expected: FAIL.
- [ ] **Step 3:** Implement `family_chains`, `admit_at_copy`, `admit`. A private chain (`c<idx>:` prefix) is a family chain of its own copy only and never lifts elsewhere.
- [ ] **Step 4:** Run the tests. Expected: PASS.
- [ ] **Step 5:** Commit: `"famasm admit: per-copy certified admission, whole and compat"`.

### Task 5: `famasm.py gtf` and the single-copy identity check

**Files:**
- Modify: `bench/stgaps/famasm.py`
- Test: `bench/stgaps/test_famasm.py`

**Interfaces:**
- Consumes: `admit/<fam>.tsv`, BASE GTF (`fw_layer.gtf_tx(path)` gives per-transcript exons and chains), `read_copies`.
- Produces: `gtf(out_dir, base_gtf, copies_path, out_gtf) -> None`: FAMASM = every BASE line verbatim, then one transcript per admitted `(copy, genomic_chain)` absent from BASE on that contig, exons = the chain's complement inside the copy locus `[locus_start, locus_end]` trimmed to the certified reads' outermost aligned bases at that copy (recorded in `admit/<fam>.tsv` as `span_start\tspan_end`, add those two columns in Task 4's writer), attributes `gene_id "FAMASM_<fam>_<idx>"; transcript_id "FAMASM_<fam>_<idx>_<n>"; reads "<certified_reads>"; famasm "1"; family "<fam>"; via "<via>"; certified_reads "<n>"`. Also `gtf.stats.tsv` (`appended\talready_in_base`).

- [ ] **Step 1: Failing tests**: `test_no_duplicate_of_base` (an admitted chain equal to a BASE chain on the contig -> not appended, `already_in_base` 1), `test_base_lines_verbatim` (every BASE line appears unchanged, in order, before any FAMASM line), `test_single_copy_untouched` (a BASE GTF with a one-copy family and no admissions -> output == input byte for byte).
- [ ] **Step 2:** Run `pytest -k "duplicate or verbatim or untouched"`. Expected: FAIL.
- [ ] **Step 3:** Implement `gtf`.
- [ ] **Step 4:** Run the tests. Expected: PASS. Then run the whole file: `python3 -m pytest bench/stgaps/test_famasm.py -v`. Expected: all PASS.
- [ ] **Step 5:** Commit: `"famasm gtf: BASE + admitted chains, byte-identical elsewhere"`.

### Task 6: `famasm.py score` and the driver for one sample

**Files:**
- Modify: `bench/stgaps/famasm.py`
- Create: `bench/stgaps/famasm_run.sh`

**Interfaces:**
- Produces: `famasm_run.sh SAMPLE PREFIX BAM GENOME RUN_DIR SPECIES` running, in order: `tools/rustle_pipeline.sh assign --bam BAM --fasta GENOME --out PREFIX --threads 72 --bin ~/rustle_target/release` (skipped when `PREFIX.assignments.tsv` has > 1 line), `famasm.py pool|lift|admit|gtf`, `stratify.py` on `base=PREFIX.families.gtf famasm=RUN_DIR/famasm.gtf`, `observed.py count/score`, `mcl_families --from-gtf RUN_DIR/famasm.gtf --fasta GENOME --threads 72 --min-exonic-bp 1 --min-shared-exon-frac 0.60 --emit-units --out RUN_DIR/famasm.fam` (realignment; `RUSTLE_CACHE_DIR=PREFIX.cache`), `mainstay_score.py --species SPECIES` and `pairscore.py --species SPECIES` against `PREFIX.fam.copies.tsv`; then `famasm.py score RUN_DIR` which prints the M1-M4 table (`verdict.txt`) from `ref/base_strata.txt`, the famasm strata file, `mainstays.txt`, `pairs.txt`, plus the descriptive counts from the stats files.
- Single-copy byte check: `famasm.py score` also diffs BASE and FAMASM restricted to contigs/loci outside multi-copy families and reports `single_copy_identical yes|no`.

- [ ] **Step 1:** Write `famasm_run.sh` (Slurm-ready: `#SBATCH -p singlenode -t 12:00:00`, the PATH/env exports of `~/work/jfim/runs/absorb_A119b/job.sbatch`).
- [ ] **Step 2:** Implement `score` (pure table assembly; the bar arithmetic as in Global Constraints; prints PASS/FAIL per bar and the overall verdict).
- [ ] **Step 3:** Smoke-run on A119b restricted to one contig: `famasm.py pool ... --contig chr16` (add a `--contig` option to `pool` that limits the copies table to that contig) and check `admit.stats.tsv` has NPIP rows; do not look at any score.
- [ ] **Step 4:** Commit: `"famasm score + run driver"`.

### Task 7: Pre-registration, dev run, held-out run

**Files:**
- Modify: `docs/LEDGER_O1.md` (section `## PREREG FAMASM (2026-10-31)`), `docs/NEGATIVE_RESULTS_REGISTER.md` (next row after 1478)

- [ ] **Step 1:** Write the prereg from the spec's Measurement section verbatim (bars M1-M4, dev A119b, held-out HG002 with the disclosure of its earlier use today), BEFORE any dev number is seen; commit and push.
- [ ] **Step 2:** `sbatch` the dev run: `famasm_run.sh A119b ~/runs/human_A119b/a1 ~/work/jfim/rustle_data/reads/A119b.t2t.bam ~/data/chm13/chm13v2.0.fa ~/work/jfim/runs/famasm_A119b human`. Record `verdict.txt` in the ledger as `### FAMASM dev result`.
- [ ] **Step 3:** Regardless of dev (reported, not judged), `sbatch` the held-out: `famasm_run.sh HG002 ~/runs/hg002_kinnex/hg ~/work/jfim/rustle_data/reads/HG002_kinnex.t2t.bam ~/data/hg002/hg002.pat.fa ~/work/jfim/runs/famasm_HG002 hg002`. Record `### FAMASM RESULT` with the M1-M4 table and the verdict; update the register row; update memory (`locus-absorb.md` or a new `famasm.md`).
- [ ] **Step 4:** If PASS: open the Rust port as the next plan (`--multicopy` on the assemble stage). If FAIL: the prototype stays in `bench/stgaps/`, nothing else changes.
