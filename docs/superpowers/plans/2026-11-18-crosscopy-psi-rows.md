# Cross-copy splicing, part 1: rows for the copies the `psi` report leaves out. Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `PREFIX.fam.psi.tsv` gains a row for every copy that has the junction's exons and reads at its anchor but
shows the junction with fewer than 3 reads (including none), and five columns that state each copy's range and gap.

**Architecture:** `psi_bunches.py` gets an optional dump of every intron's raw counts. `psi_family.py` projects a
junction from the copy that uses it to the other members through the family's pair alignments
(`bench/denovo_eval/exonshare.py`), looks the projected junction up in the dump, and appends rows and columns. The
driver's `psi` stage passes the two new inputs when they exist.

**Tech Stack:** Python 3 (stdlib + pysam already used by `psi_bunches.py`), pytest, bash driver. Run tests with the
pixi env first on PATH: `~/.pixi/detached-envs/rustle-7434162673922408322/envs/default/bin`.

**Spec:** `docs/superpowers/specs/2026-11-18-cross-copy-splicing-design.md` (sections 1 and 2).

## Global Constraints
- `PREFIX.psi.tsv` stays byte-identical. In `PREFIX.fam.psi.tsv` the rows written today keep their 26 columns
  byte-identical and their order; added rows come after them; new columns are appended to every row.
- Opt-in report only: nothing becomes part of `all`, nothing else reads these tables.
- No sample name, path or species constant in `tools/`.
- Shared repository: never `git add -A`; `git diff --cached --stat` before every commit; commit only this plan's paths,
  and only when the user asks for the commit (project rule). The "Commit" steps below mean: stage by explicit path and
  report; the commit itself waits for the user's word.
- `tools/rustle_pipeline.sh` is read by running jobs: write the new driver to a temp file, `bash -n` it, then
  `os.replace` / `mv -f` it into place.
- Eligible copy = at least `--min-bunch` (default 20) certain reads; silent = no read at all in the wide band.

## Review Focus
1. Members on opposite strands (PAF strand `-`): the projected junction must come back with start <= end on the target
   and be looked up under the target member's strand.
2. A junction end that falls in an alignment gap, or whose two ends project through different records: no row, counted
   as unlocated in the log, no exception.
3. `psi_family.py` run without `--introns` / `--paf` (older run directories, the driver before `exonmap`): no added
   rows, the new columns still present, exit 0.
4. A key with one eligible copy: `gap` is `NA`, `n_eligible` is 1.
5. A member that already has a row under the key in another bunch: no second, added row for it.

---

### Task 1: raw intron counts from `psi_bunches.py`

**Files:**
- Modify: `tools/psi_bunches.py` (new function beside `bunch_rows`; one argument in `main`)
- Test: `tools/test_psi_bunches.py`

**Interfaces:**
- Produces: `INTRON_HEADER = "chrom strand start end placed max placed_w max_w".split()`;
  `intron_rows(cnt) -> list[list]`, one row per intron of `tally`'s result with `max_w >= 1`, sorted by
  (chrom, start, end, strand); CLI flag `--introns OUT.tsv` (written through `OUT.tsv.tmp` + `os.replace`).

- [ ] **Step 1: Write the failing test** in `tools/test_psi_bunches.py`

```python
def test_intron_rows_keep_what_bunches_drop():
    cnt = {A: [12, 14, 11, 14], B: [2, 3, 2, 3], C: [0, 1, 0, 1], X: [5, 5, 5, 5], ("c1", "-", 5, 9): [0, 0, 0, 0]}
    rows = intron_rows(cnt)
    assert [tuple(r[:4]) for r in rows] == [("c1", "+", 101, 200), ("c1", "+", 101, 300), ("c1", "+", 251, 300), ("c2", "+", 11, 90)]
    assert rows[2][4:] == [0, 1, 0, 1] and len(INTRON_HEADER) == len(rows[0]) == 8   # C: under the 3-read floor, kept here
```

- [ ] **Step 2: Run it, expect an ImportError** — `python3 -m pytest tools/test_psi_bunches.py -q`
- [ ] **Step 3: Implement `intron_rows` and `--introns`** (the dump is written after `PREFIX.psi.tsv`, from the same `cnt`)
- [ ] **Step 4: Tests pass** — `python3 -m pytest tools/test_psi_bunches.py -q` → all passed
- [ ] **Step 5: Byte check on real data.** In a scratch directory:
  `python3 tools/psi_bunches.py $W/psi_2026-11-16/or.gtf $W/reads/GGO/GGO_mm.bam $W/reads/GGO/GGO_mm.bam.as_table.tsv a.tsv --contigs NC_073224.2 --threads 8`
  and the same with `--introns a.introns.tsv` into `b.tsv` (`W=/anvme/workspace/v138da13-rustle`); `cmp a.tsv b.tsv`
  → identical; `a.introns.tsv` has more rows than `a.tsv`.
- [ ] **Step 6: Stage** `tools/psi_bunches.py tools/test_psi_bunches.py`

### Task 2: projection of a junction between members

**Files:**
- Modify: `tools/psi_family.py`
- Test: `tools/test_psi_family.py`

**Interfaces:**
- Consumes: `blocks_from_cigar(cg, qs, qe, ts, strand)`, `seg_map(blocks, strand, from_q)`, `project(segs, a, b)` from
  `bench/denovo_eval/exonshare.py` (import through `sys.path` as `tools/` scripts already do for `bench/`). PAF names
  are member keys `CONTIG:START-END` (START 1-based); PAF coordinates are 0-based within the member's span.
- Produces:
  - `load_pair_records(paf, members) -> {(src, dst): [segs, ...]}`: for every PAF record whose two names are both in
    `members` and differ, the segment maps in both directions (genomic 1-based coordinates on both sides), each list
    ordered by the record's matching bases (PAF column 10), most first.
  - `project_intron(records, src, dst, start, end) -> (start2, end2) | None`: `start`, `end` = first and last intron
    base (1-based closed) on `src`. Projects the two exonic bases flanking the intron (`start - 1`, `end + 1`) through
    the first record of `records[(src, dst)]` that maps both; the junction on `dst` is the open interval between the
    two images. `None` when no record maps both bases or the interval is shorter than 20 bases.

- [ ] **Step 1: Write the failing tests** (a PAF written to `tmp_path`; members `c1:1-500` and `c2:1-500`)

```python
def test_project_intron_plus_minus_and_gap(tmp_path):
    paf = tmp_path / "p.paf"   # c1:1-500 vs c2:1-500, 500 bases, minus strand, one 500= block; then a second pair with a gap
    paf.write_text("c1:1-500\t500\t0\t500\t-\tc2:1-500\t500\t0\t500\t500\t500\t0\tcg:Z:500=\n"
                   "c1:1-500\t500\t0\t500\t+\tc3:1-490\t490\t0\t490\t490\t500\t0\tcg:Z:95=10I395=\n")
    R = load_pair_records(str(paf), {"c1:1-500", "c2:1-500", "c3:1-490"})
    assert project_intron(R, "c1:1-500", "c2:1-500", 101, 200) == (301, 400)      # minus: base x maps to 501 - x
    assert project_intron(R, "c2:1-500", "c1:1-500", 301, 400) == (101, 200)      # and back
    assert project_intron(R, "c1:1-500", "c3:1-490", 201, 300) == (191, 290)      # past the 10-base insertion in c1
    assert project_intron(R, "c1:1-500", "c3:1-490", 101, 200) is None            # exonic base 100 sits in the insertion
    assert project_intron(R, "c1:1-500", "c9:1-9", 101, 200) is None              # no record for the pair
```

- [ ] **Step 2: Run, expect ImportError** — `python3 -m pytest tools/test_psi_family.py -q`
- [ ] **Step 3: Implement both functions**
- [ ] **Step 4: Tests pass**
- [ ] **Step 5: Stage** `tools/psi_family.py tools/test_psi_family.py`

### Task 3: added rows and the five new columns

**Files:**
- Modify: `tools/psi_family.py` (`key_rows`, `HEADER`, `main`, module docstring)
- Test: `tools/test_psi_family.py`

**Interfaces:**
- Consumes: Task 1's dump (`INTRON_HEADER` columns), Task 2's two functions, `psi(p, m, sp, sm)` from
  `tools/psi_bunches.py`.
- Produces:
  - `HEADER` = today's 26 names + `silent n_eligible range_lo range_hi gap`.
  - `key_rows(rows, cols, min_bunch=20, introns=None, records=None) -> list[list]`; with `introns` (dict
    `{(chrom, strand, start, end): [placed, max, placed_w, max_w]}`) and `records` it appends added rows.
  - CLI: `--introns PSI.introns.tsv --paf FAM.loci.gated.paf` (both optional; added rows need both).
  - Log: today's line unchanged, then one line
    `added rows N (silent S) for K keys; M candidate copies not located; eligible copies counted with added rows`.
- Rules the code must follow (spec section 2):
  - `range_lo`, `range_hi`, and the reads `n` of a copy under a key are exactly the values `key_rows` already computes
    in `members[key][member]`; they are printed on every row of that member under the key.
  - Candidate for an added row: for each key and each member U that has a row under it, every other member S of the
    family that owns both of the key's columns and has no row under the key. U's junction = U's row with the most
    `placed_w` in its deepest bunch; project it to S with `project_intron`.
  - S's anchor site = the projected junction's donor or acceptor end, by the key's `anchor` and S's strand. Its reads
    `n` = sum of `placed_w` over the dump's introns sharing that site on S's strand. `n < min_bunch`: no row.
  - The added row carries S's projected coordinates; `placed`, `max`, `placed_w`, `max_w` from the dump (0 when
    absent); the PSI columns from `psi()` over all the dump's introns at that site; `locus` `.`, `in_gtf` 0;
    `bunch_placed_w` = `n`; `n_copies` and `differs` 0; `silent` 1 when `max_w` is 0, else 0.
  - If the dump shows `max_w >= 3` for S's projected junction and S still has no row (its site holds one intron),
    the row is added the same way (`silent` 0).
  - Several U for one (key, S): use the U with the most reads; one added row per (key, S).
  - `n_eligible` = copies under the key with `n >= min_bunch`, added rows included. `gap` = this copy's `range_lo`
    minus the lowest `range_hi` among the other eligible copies, 4 decimals; `NA` when this copy is not eligible or
    has no eligible partner.

- [ ] **Step 1: Write the failing tests.** Reuse `setup` and `row` of the test file.

```python
def test_old_columns_unchanged_and_new_ones(tmp_path):
    cols = setup(tmp_path)
    rows = [row("a", "donor", "c1", "+", 101, 200, 18, 0.9, 0.95), row("a", "donor", "c1", "+", 101, 400, 2, 0.05, 0.1),
            row("b", "donor", "c2", "-", 301, 400, 4, 0.1, 0.3), row("b", "donor", "c2", "-", 101, 400, 26, 0.7, 0.9)]
    out = key_rows(rows, cols, min_bunch=20)
    assert len(out) == 4 and all(len(r) == len(HEADER) == 31 for r in out)
    by = {(r[1], r[5]): r for r in out}
    a = by[("c1:1-500", "C1")]
    assert a[24:26] == [2, 1] and a[26:31] == [0, 2, 0.9, 0.95, 0.6]            # 0.9 - 0.3
    assert by[("c2:1-500", "C1")][30] == -0.85                                   # 0.1 - 0.95: below its partner

def test_silent_copy_gets_a_row(tmp_path):
    cols = setup(tmp_path)                                                       # B is A reversed: x on c1 = 501 - x on c2
    (tmp_path / "p.paf").write_text("c1:1-500\t500\t0\t500\t-\tc2:1-500\t500\t0\t500\t500\t500\t0\tcg:Z:500=\n")
    rec = load_pair_records(str(tmp_path / "p.paf"), {"c1:1-500", "c2:1-500"})
    rows = [row("a", "donor", "c1", "+", 101, 200, 18, 0.6, 0.6), row("a", "donor", "c1", "+", 101, 400, 12, 0.4, 0.4)]   # A skips C1 at 0.4
    intr = {("c2", "-", 301, 400): [40, 40, 40, 40]}                             # B: only C0->C1 from that donor, never the skip
    out = key_rows(rows, cols, min_bunch=20, introns=intr, records=rec)
    add = [r for r in out[2:]]
    skip = next(r for r in add if r[5] == "C2")
    assert skip[1] == "c2:1-500" and skip[10:12] == ["101", "400"] and skip[26] == 1 and skip[28:30] == [0.0, 0.0]
    incl = next(r for r in add if r[5] == "C1")                                  # B's constitutive junction: PSI 1, not silent
    assert incl[26] == 0 and incl[28:30] == [1.0, 1.0] and incl[18] == 40
    a_skip = next(r for r in out[:2] if r[5] == "C2")
    assert a_skip[27] == 2 and a_skip[30] == 0.4                                 # eligible with the added row; gap 0.4 - 0
    assert key_rows(rows, cols, min_bunch=20, introns={("c2", "-", 301, 400): [5, 5, 5, 5]}, records=rec)[2:] == []   # 5 reads: not eligible
```

- [ ] **Step 2: Run, expect failures** (31 columns, missing parameters)
- [ ] **Step 3: Implement** (the existing loop that fills `members` stays the single source of a copy's range)
- [ ] **Step 4: Tests pass**, including the untouched `test_keys_and_differs` after widening its length assertions to
  the new `HEADER`
- [ ] **Step 5: Rewrite the module docstring**: the sentence "a copy that never uses the junction has no row at all" is
  replaced by the added-row rule and the five columns, each in one line.
- [ ] **Step 6: Stage** `tools/psi_family.py tools/test_psi_family.py`

### Task 4: the `psi` stage passes the new inputs; real-data check

**Files:**
- Modify: `tools/rustle_pipeline.sh` (`stage_psi`, the header lines of `psi`), `README.md` (the `psi` row)

**Interfaces:**
- Consumes: Tasks 1-3's CLI flags.
- Produces: `PREFIX.psi.introns.tsv` written by `stage_psi` only when `PREFIX.fam.exonmap.columns.tsv` exists and is
  not older than `PREFIX.fam.loci.gff3` (the condition the stage already tests); `psi_family.py` called with
  `--introns` and, when `PREFIX.fam.loci.gated.paf` exists, `--paf`. Without the PAF the stage says so in its `say`
  line and writes no added rows.

- [ ] **Step 1: Edit a copy of the driver, `bash -n`, replace atomically.** The stage's failure handling and messages
  keep their form.
- [ ] **Step 2: Real-data check on the development sample.** A fresh directory in the workspace with the same symlinks
  as `psi_2026-11-16` (`or.gtf`, `or.fam.clusters.tsv`, `or.fam.loci.gff3`, `or.fam.loci.gated.paf`) and a copy of
  `or.fam.exonmap.columns.tsv`; run the two driver stages as `psi_2026-11-16/or2.sbatch` does (Slurm job, 2 h limit,
  genome from the workspace copy `genomes/GGO/GGO.fa`). Expected:
  `cmp <(cut -f1-26 NEW/or.fam.psi.tsv | head -n 3653) <(cut -f1-26 psi_2026-11-16/or.fam.psi.tsv)` identical;
  `cmp NEW/or.psi.tsv psi_2026-11-16/or.psi.tsv` identical; the log's first line equals the stored
  `or.fam.psi.log`; the second reports added rows > 0. Record the second line in the plan's progress note.
- [ ] **Step 3: Driver tests** — run the repository's driver dispatch test the way the famwise stage's commit did
  (`git show --stat 59639df4` names it) → passes.
- [ ] **Step 4: README row**: one clause added to the `psi` row (added rows for copies under the read floor; columns
  `silent`, `n_eligible`, `range_lo`, `range_hi`, `gap`).
- [ ] **Step 5: Stage** `tools/rustle_pipeline.sh README.md`; report the staged paths with `git diff --cached --stat`.
