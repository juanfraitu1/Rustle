# Family HMM layer Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** An opt-in bench script that builds a DNA and a protein profile HMM per multi-copy family, labels families with Pfam domains, and reports a pre-registered concordance table against the MCL families.

**Architecture:** One Python module `bench/famhmm/famhmm.py` (pure helpers + subprocess calls to HMMER/MAFFT) with subcommands build / scan / concord / all; unittest file beside it in the `bench/unmapped_rescue` style; nhmmer tblout parsing reused from `bench/unmapped_rescue/profiles.py`. Nothing in the driver, the Rust binaries or shipped outputs changes.

**Tech Stack:** Python 3 stdlib (subprocess, csv, argparse, tempfile), HMMER 3.4 (hmmbuild, hmmpress, nhmmer, hmmscan, hmmemit, esl-translate) and MAFFT from the pixi `bench` env.

**Spec:** docs/superpowers/specs/2026-11-05-family-hmm-design.md

## Global Constraints

- Run with `pixi run -e bench python3 -B ...`; tests skip (not fail) when `hmmbuild` is not on PATH.
- Outputs only under `--out` (workspace `/anvme/workspace/v138da13-rustle/famhmm/<sample>`), never home; no file under PREFIX is written.
- Defaults from the spec: `--min-copies 2`, `--min-orf-aa 50`, Pfam domain E <= 1e-5, concordance hit E <= 1e-10 and >= 50% profile cover, rule thresholds 0.95 and 0.95.
- Pfam absent and not `--no-pfam`: print the download command, exit 2.
- First real run: the 57-family panel `bench/o2levers/panel.families.txt` with a 1 h wall; genome-wide only after.
- Copy header format `FAM|idx|chrom:start-end|strand|nexon=N`; loci header `chrom:start-end`; a copy's locus key is field 3 of its header.

## Review Focus

1. A family whose copies are all on the minus strand or mixed strands: MAFFT alignment must use the FASTA as written (copies are already transcript-oriented); test pins that a reverse-complemented copy is still recovered by nhmmer (hits on `-` strand count).
2. A family with one copy with an ORF and one without: protein profile must be skipped with `prot_status=none`, DNA profile still built (Task 2 test).
3. Two profiles hitting the same locus: `n_top_own` must use the best bit score across ALL families, not the own family's best (Task 4 test).
4. A family in `--only-families` but absent from copies.tsv: skipped with a stderr line, not a crash (Task 2 test).
5. Zero families with >= 2 copies: concord writes a header-only table and the rule prints `FAIL (no families)` (Task 4 test).

---

### Task 1: ORF selection and FASTA helpers

**Files:**
- Create: `bench/famhmm/__init__.py` (empty), `bench/famhmm/famhmm.py`
- Test: `bench/famhmm/test_famhmm.py`

**Interfaces:**
- Produces: `read_fasta(path) -> list[(header, seq)]`; `copies_by_family(copies_fa) -> dict[str, list[(header, seq)]]` keyed by the header's first `|` field; `longest_orf(records) -> (copy_header, frame:int, nt_offset:int, aa_seq) | None` where records are the `esl-translate` FASTA records for ONE copy (headers like `orf1 source=X coords=10..129 length=40 frame=2 desc=...`), picks the longest, returns frame and 0-based nt offset parsed from `coords=`; `esl_translate(copy_fa, min_aa) -> str` returns the raw FASTA text from `esl-translate -l min_aa --watson`.

- [ ] **Step 1: Write failing tests** `test_copies_by_family_groups_on_first_field`, `test_longest_orf_picks_longest_and_parses_coords` (records with lengths 40 and 60, frame=2 coords=10..129 -> (header, 2, 9, aa)) and `test_longest_orf_none_when_empty`.
- [ ] **Step 2: Run** `pixi run -e bench python3 -B -m unittest bench/famhmm/test_famhmm.py -v` Expected: FAIL, import error.
- [ ] **Step 3: Implement the four functions** in `famhmm.py`; `esl_translate` uses `subprocess.run([... "esl-translate", "-l", str(min_aa), "--watson", copy_fa], capture_output=True, text=True, check=True).stdout`.
- [ ] **Step 4: Run tests** Expected: PASS.
- [ ] **Step 5: Commit** `git add bench/famhmm && git commit -m "famhmm: FASTA and ORF helpers"`

### Task 2: build

**Files:**
- Modify: `bench/famhmm/famhmm.py`
- Test: `bench/famhmm/test_famhmm.py`

**Interfaces:**
- Consumes: Task 1 helpers.
- Produces: `build(prefix, out, only=None, min_copies=2, min_orf_aa=50, threads=1) -> list[dict]` one dict per family `{family, n_copies, prot_status: "ok"|"none"}`; writes `out/dna/FAM.hmm`, `out/prot/FAM.hmm` (when >= 2 ORF copies), `out/orf.tsv` (columns family, tid, frame, nt_offset, aa_len), `out/build.tsv` (the returned dicts), concatenated + pressed `out/families.dna.hmm`, `out/families.prot.hmm`. Helper `build_profile(fasta, hmm_out, name, alphabet)` runs `mafft --auto --quiet` then `hmmbuild --<alphabet> -n name`.
- Family selection: families in `copies_by_family` with `len >= min_copies`, intersected with `only` (a file of family ids, one per line) when given; ids in `only` absent from the table go to stderr as `skip FAM: not in copies`.

- [ ] **Step 1: Write failing tests** (skipped when `shutil.which("hmmbuild")` is None): `test_build_toy_family` writes a 3-copy toy family (one 180-nt ORF-bearing sequence with planted frame 1 at offset 0 plus two 2-base-mutated copies, 1 copy reverse-complemented in the FASTA is NOT needed here) and asserts `out/dna/FAM1.hmm`, `out/prot/FAM1.hmm`, `families.dna.hmm.h3m` exist and `orf.tsv` has 3 rows with frame 1 offset 0; `test_build_prot_none_when_single_orf` (3 copies, only one with an ORF -> prot_status none, no `prot/FAM1.hmm`); `test_build_only_missing_family_warns` (only=["NOPE"] -> empty result, stderr contains `skip NOPE`).
- [ ] **Step 2: Run** Expected: FAIL, `build` not defined.
- [ ] **Step 3: Implement `build`, `build_profile`** ; per-family temp files under `tempfile.TemporaryDirectory`; hmmpress with `-f`.
- [ ] **Step 4: Run tests** Expected: PASS.
- [ ] **Step 5: Commit** `git commit -am "famhmm: build DNA and protein profiles per family"`

### Task 3: scan (Pfam labels)

**Files:**
- Modify: `bench/famhmm/famhmm.py`
- Test: `bench/famhmm/test_famhmm.py`

**Interfaces:**
- Consumes: `out/families.prot.hmm` from Task 2.
- Produces: `parse_domtblout(lines, max_evalue=1e-5) -> dict[query, list[(ali_from, pfam_name, ievalue)]]` (columns: target name f[0], query name f[3], i-Evalue f[12], ali from f[17]); `scan(out, pfam, max_evalue=1e-5) -> None` runs `hmmemit -c out/families.prot.hmm` to `out/prot.consensus.fa`, `hmmscan --domtblout out/pfam.domtblout --cpu T pfam out/prot.consensus.fa`, writes `out/pfam_labels.tsv` (family, domains joined by `;` in ali_from order, best_evalue, prot_status from build.tsv). `pfam_path(arg) -> str|None`: arg, else `$WS/pfam/Pfam-A.hmm` with `WS=/anvme/workspace/v138da13-rustle`; missing -> print `wget -P $WS/pfam https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz && gunzip ... && hmmpress ...` and `sys.exit(2)`.

- [ ] **Step 1: Write failing test** `test_parse_domtblout_orders_domains_by_position` with two synthetic domtblout lines for one query (ali_from 120 then 5, one above the E cut) -> ordered list of the names under the cut.
- [ ] **Step 2: Run** Expected: FAIL.
- [ ] **Step 3: Implement `parse_domtblout`, `scan`, `pfam_path`**.
- [ ] **Step 4: Run tests** Expected: PASS.
- [ ] **Step 5: Commit** `git commit -am "famhmm: Pfam domain labels per family"`

### Task 4: concord and the decision rule

**Files:**
- Modify: `bench/famhmm/famhmm.py`
- Test: `bench/famhmm/test_famhmm.py`

**Interfaces:**
- Consumes: `bench/unmapped_rescue/profiles.parse_tblout` (sys.path insert of that dir); `out/families.dna.hmm`; copy header field 3 as locus key.
- Produces: `profile_lengths(hmm_path) -> dict[name, int]` (LENG lines); `concord_rows(hits, members, lengths, min_cover=0.5) -> list[dict]` where `hits` is parse_tblout output (already E-filtered at 1e-10), `members = {family: set(locus keys)}`; row keys family, n_copies, n_recovered, n_extra, n_top_own; a hit counts when `(end-start+1) >= min_cover * lengths[family]`; `n_top_own` counts member loci whose max-score counted hit over ALL families belongs to the own family. `decide(rows) -> (str, float, float)` returns `("CONCORDANT"|"FAIL", rec_frac, top_frac)` with 0.95/0.95, `("FAIL (no families)", 0, 0)` on empty rows. `concord(prefix, out, threads) -> None` runs `nhmmer --tblout out/loci.tblout --cpu T out/families.dna.hmm PREFIX.fam.loci.fa`, writes `out/concordance.tsv` plus footer `# rule: VERDICT recovered=R top_own=T`.

- [ ] **Step 1: Write failing tests** (pure, no HMMER): `test_concord_rows_counts_recovered_extra_top_own` with two families, one shared locus where FAM2 outscores FAM1 (so FAM1 n_top_own excludes it); `test_concord_rows_cover_filter` (a 100-bp hit on a 400-length profile not counted); `test_decide_thresholds` (rows giving 0.96/0.94 -> FAIL, 0.95/0.95 -> CONCORDANT, [] -> FAIL (no families)). HMMER-dependent: `test_concord_toy_end_to_end` builds the Task 2 toy family where copy 3 is written reverse-complemented into a loci FASTA, runs `concord`, asserts n_recovered 3, n_top_own 3, footer CONCORDANT.
- [ ] **Step 2: Run** Expected: FAIL.
- [ ] **Step 3: Implement the four functions**.
- [ ] **Step 4: Run tests** Expected: PASS.
- [ ] **Step 5: Commit** `git commit -am "famhmm: concordance table and pre-registered rule"`

### Task 5: CLI, `all`, panel run

**Files:**
- Modify: `bench/famhmm/famhmm.py` (argparse `main`), `docs/superpowers/specs/2026-11-05-family-hmm-design.md` (append a `## Result` line after the run)

**Interfaces:**
- CLI exactly as the spec: `famhmm.py {build,scan,concord,all} --prefix P --out D [--only-families F] [--threads N] [--no-pfam] [--pfam PATH] [--min-copies 2] [--min-orf-aa 50]`; `all` = build, scan unless `--no-pfam`, concord.

- [ ] **Step 1: Implement `main`** and `if __name__ == "__main__": main()`.
- [ ] **Step 2: Smoke** `pixi run -e bench python3 -B bench/famhmm/famhmm.py all --prefix ~/runs/human_A119b/a1 --out /anvme/workspace/v138da13-rustle/famhmm/A119b_panel --only-families bench/o2levers/panel.families.txt --threads 8 --no-pfam` on a compute node under `timeout 3600` (sbatch singlenode; the storage hook gates it). Expected: `concordance.tsv` with rows for the panel families present in copies.tsv and a footer verdict; record wall time.
- [ ] **Step 3: Append `## Result`** to the spec: panel verdict, recovered/top_own fractions, wall, Pfam skipped or not. A FAIL is written as FAIL.
- [ ] **Step 4: Commit** `git commit -am "famhmm: CLI; panel result"`

Genome-wide and the Pfam download are separate decisions after the panel wall time is known.
