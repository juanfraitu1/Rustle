# Cross-copy splicing, part 2: the tests against sequence cause (A) and read origin (C). Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Score the differences `PREFIX.fam.psi.tsv` reports between the copies of a family against the splice-site
sequence (A) and against the origin of the reads (C), and count which assembler emits the copy-specific junctions.

**Architecture:** Scorers in a new `bench/crosscopy/`, all reading calls from one function over `PREFIX.fam.psi.tsv`.
Site strength is the PWM code of row 1789 moved into the repository. Origin reuses `bench/famwise_eval/records.py` and
`bench/o2_truth/neutral_truth.py` unchanged. A pre-registration is written and fixed before any judging sample is
opened; the development sample is gorilla OR6737.

**Tech Stack:** Python 3 (stdlib, pysam), pytest, Slurm job scripts in the workspace, samtools, minimap2.

**Spec:** `docs/superpowers/specs/2026-11-18-cross-copy-splicing-design.md` (sections 1, 3, 4, 5, 6 and "Rulings").
Part 1 (`docs/superpowers/plans/2026-11-18-crosscopy-psi-rows.md`) must be finished first.

## Global Constraints
- Classes: LARGE gap > 0.2; SMALL 0 < gap <= 0.2; SAME = ranges overlap and the eligible copies' `psi_placed_w` within
  0.05; COPY-SPECIFIC = LARGE with `range_hi` of S <= 0.05. Eligible = >= 20 certain reads (`--min-bunch`).
- Counted unit = (reference genome, U's junction, S). A unit seen in two samples of one genome counts once, with the
  numbers of the sample holding more certain reads at U. A unit that is a call on OR6737 counts in no bar.
- A: tau = 1 bit (0.5 and 2 printed beside). Bar A1: AGREES share among AGREES + OPPOSES, Wilson 95% lower bound > 0.5,
  at least 30 units, else undecidable.
- C: judged kinds `one` and `several`; `outside` printed beside (ruling 2). Bar C1: SUPPORTED / judged >= 0.95 with
  Wilson lower bound >= 0.90. Bar C2: >= 0.90 of judged LARGE calls hold; a call is judged with >= 20 SUPPORTED reads in
  the bunch at U and at S; fewer than 20 judged calls = undecidable. HG002 reads of every hash residue (ruling 1).
- No count on a judging sample before the pre-registration is committed. OR6737 is seen and debugs everything.
- No annotation is read by any call or verdict. No path, sample or species name in `bench/crosscopy/*.py`; they live in
  the workspace job scripts.
- Inputs from the workspace (`/anvme/workspace/v138da13-rustle`, below `$W`), not atuin; sort and alignment temps in
  `$TMPDIR`; nothing new under `$HOME`. `PYTHONHASHSEED=0` in every job.
- Shared repository: commit by explicit path, `git diff --cached --stat` first, only when the user asks. Register rows:
  claim each number by creating its workspace directory, re-check the register tail right before writing the row.
- Do not write into `artifact_refresh_2026-11-16/` (published pages are rebuilt from it): link its files into new run
  directories.

## Review Focus
1. A projected site within 30 bases of a contig end: `windows` returns a short string and the score is NaN. The call is
   UNLOCATED, not AGREES / OPPOSES, and nothing raises.
2. A truth read whose origin index lies beyond its 12 written loci (`kind` unjudged in `records.load_truth`): in neither
   C1's numerator nor its denominator.
3. The same unit called by a testis and a cell-line sample of one genome, with U and S swapped between them: two units
   (different U's junction), each counted once; never merged into one with a sign flip.
4. An arm GTF with `transcript_id` equal on two contigs, or with gzip input: the junction lookup keys on
   (contig, strand, start, end) and opens `.gz`.
5. A LARGE call whose S is an added silent row and whose U's bunch holds one intron only: A has no competitor at
   either copy, so dA is the difference of the key junction's own site strength; it must not divide or take a max of
   an empty list.

---

### Task 1: calls

**Files:**
- Create: `bench/crosscopy/calls.py`
- Test: `bench/crosscopy/test_calls.py`

**Interfaces:**
- Consumes: `PREFIX.fam.psi.tsv` with part 1's columns (`silent n_eligible range_lo range_hi gap`).
- Produces:
  - `Call = namedtuple("Call", "family anchor col_anchor col_other cls gap u s u_intron s_intron reads_u reads_s lo_u hi_s silent copy_specific")`
    with `u`, `s` member keys, `u_intron` / `s_intron` = `(chrom, strand, start, end)` (1-based closed; `s_intron` is
    S's row with the most `placed_w`, or the added row's), `cls` in `LARGE SMALL SAME`.
  - `load_calls(fam_psi_tsv, min_bunch=20) -> list[Call]`, one per key at most (keys with one eligible copy or with
    neither SAME nor a positive gap give none), after the dedup "per (family, `u_intron`) keep the larger gap",
    ordered by (family, anchor, col_anchor, col_other). A copy's range and reads are read from the rows themselves
    (`range_lo`, `range_hi`, the largest `bunch_placed_w` of its rows under the key; its `psi_placed_w` = the sum over
    the rows of that bunch), never from the `gap` / `n_eligible` columns, so a table cut to a subset of rows still
    gives that subset's calls.
  - `unit(call, genome) -> tuple` = `(genome, *call.u_intron, call.s)`.
  - CLI `calls.py FAM.psi.tsv [OUT.tsv]`: class counts on stdout, the calls as a table when OUT is given.

- [ ] **Step 1: Write the failing tests** (rows built as dicts with the 31 column names; helper in the test file)

```python
def test_classes_dedup_and_units(tmp_path):
    p = write(tmp_path, [   # family, member, anchor, col_a, col_o, chrom, strand, start, end, placed_w, bunch, lo, hi, psi_placed_w, silent
        ("F", "c1:1-500", "donor", "C0", "C2", "c1", "+", 101, 400, 12, 30, 0.4, 0.4, 0.4, 0),
        ("F", "c2:1-500", "donor", "C0", "C2", "c2", "-", 101, 400, 0, 40, 0.0, 0.0, 0.0, 1),
        ("F", "c1:1-500", "acceptor", "C2", "C0", "c1", "+", 101, 400, 12, 25, 0.3, 0.3, 0.3, 0),   # the same intron in its acceptor bunch
        ("F", "c2:1-500", "acceptor", "C2", "C0", "c2", "-", 101, 400, 0, 40, 0.0, 0.0, 0.0, 1),
        ("G", "c3:1-9", "donor", "C5", "C6", "c3", "+", 11, 90, 50, 50, 0.50, 0.52, 0.51, 0),
        ("G", "c4:1-9", "donor", "C5", "C6", "c4", "+", 11, 90, 60, 60, 0.51, 0.55, 0.53, 0),
        ("H", "c5:1-9", "donor", "C1", "C2", "c5", "+", 11, 90, 9, 9, 0.9, 0.9, 0.9, 0),            # 9 reads: not eligible
        ("H", "c6:1-9", "donor", "C1", "C2", "c6", "+", 11, 90, 50, 50, 0.1, 0.1, 0.1, 0)])
    C = load_calls(p)
    assert [(c.family, c.cls) for c in C] == [("F", "LARGE"), ("G", "SAME")]
    f = C[0]
    assert (f.anchor, f.gap, f.u, f.s, f.silent, f.copy_specific) == ("donor", 0.4, "c1:1-500", "c2:1-500", True, True)
    assert f.u_intron == ("c1", "+", 101, 400) and f.s_intron == ("c2", "-", 101, 400)
    assert unit(f, "GGO") == ("GGO", "c1", "+", 101, 400, "c2:1-500")
```

- [ ] **Step 2: Run, expect ImportError** — `python3 -m pytest bench/crosscopy/test_calls.py -q`
- [ ] **Step 3: Implement `calls.py`**
- [ ] **Step 4: Tests pass**
- [ ] **Step 5: Development-sample check.** On the rows the stage wrote before part 1
  (`head -n 3653 NEW/or.fam.psi.tsv > old_rows.tsv`, part 1 Task 4's directory; scratch file),
  `python3 bench/crosscopy/calls.py old_rows.tsv` must give 32 LARGE calls, 7 of them COPY-SPECIFIC (the count of the
  design's review on the stored table). Another number: stop and find the rule that differs. Then run it on the whole
  table and note the totals with added rows in the progress note.
- [ ] **Step 6: Stage** the two files

### Task 2: site strength and the sequence-cause scorer

**Files:**
- Create: `bench/crosscopy/sites.py` (from `$W/r1789_splicesites/ssfeat.py`: `rc`, `windows`, `don`, `acc`, `pwm`,
  `score` and its `selftest`, unchanged in behaviour; plus `sample_pwms`)
- Create: `bench/crosscopy/site_cause.py`
- Test: `bench/crosscopy/test_site_cause.py`

**Interfaces:**
- Consumes: `load_calls` (Task 1); `load_pair_records`, `project_intron` (`tools/psi_family.py`, part 1); the intron
  dump `PREFIX.psi.introns.tsv`.
- Produces:
  - `sites.sample_pwms(gtf, fasta) -> (donor_pwm, acceptor_pwm)` from the introns of the GTF's multi-exon transcripts.
  - `sites.other_end_strength(fasta, pwms, intron, anchor) -> float`: strength of the intron's non-anchor site
    (`anchor` = `donor` or `acceptor` names the SHARED end, so the scored site is the other one); NaN when the window
    is short or holds N.
  - `site_cause.advantage(key_strength, competitor_strengths) -> float` = key minus the highest competitor, or the key
    strength alone when there is no competitor.
  - `site_cause.verdict(dA, tau=1.0) -> "AGREES" | "OPPOSES" | "NO SITE DIFFERENCE"`; NaN anywhere -> `"UNLOCATED"`.
  - CLI `site_cause.py score FAM.psi.tsv PSI.introns.tsv FAM.loci.gated.paf GTF FASTA OUT.tsv`: one row per call with
    `unit` fields, class, gap, `a_u`, `a_s`, `dA`, verdicts at tau 0.5 / 1 / 2, `canonical_lost` (GT or GC donor / AG
    acceptor at U's key site and not at S's), `n_competitors`, `dropped_competitors`, `exon_identity` (matching bases
    over aligned bases of the record used).
  - CLI `site_cause.py pool --dev DEV.tsv GENOME=OUT.tsv [...]`: drops units present in DEV, merges units per genome
    (more `reads_u` wins), prints A1 (count, share, Wilson interval from `bench/famwise_eval/records.wilson`), A2 (site
    difference share LARGE vs SAME, overall and split at `exon_identity` 0.98), per-genome and per-family tables, and
    the share of LARGE units in the three largest families.
- Rules (spec section 3):
  - Competitors at copy X = the dump's introns on X's strand sharing the anchor site of X's key junction with
    `max_w >= 3`, the key junction excluded. The compared set is the union over U and S, each intron carried to the
    other copy with `project_intron`; two introns are the same when their coordinates coincide at either copy.
  - A competitor that does not project is dropped when it holds under 10% of its site's `placed_w`
    (`dropped_competitors` counts them); otherwise the call is UNLOCATED. A key junction that does not project (S
    without `s_intron`): UNLOCATED.
  - `a_x` uses, for every intron of the compared set, its coordinates AT copy X.

- [ ] **Step 1: Write the failing tests**

```python
def test_advantage_and_verdict():
    assert advantage(8.0, [3.0, 6.5]) == 1.5 and advantage(8.0, []) == 8.0
    assert verdict(1.0) == "AGREES" and verdict(-1.0) == "OPPOSES" and verdict(0.99) == "NO SITE DIFFERENCE"
    assert verdict(float("nan")) == "UNLOCATED" and verdict(0.6, tau=0.5) == "AGREES"

def test_broken_acceptor_at_the_silent_copy(tmp_path):
    # two copies of one 3-exon gene on c1 and c2 (plus strand, identical but for c2's exon-2 acceptor AG -> AA);
    # c1 skips exon 2 in 40% of reads, c2 never includes it: the INCLUSION junction is the call, U = c1, S = c2
    fa, psi, dump, paf, gtf = build_two_copies(tmp_path, break_acceptor_on="c2")
    main(["score", psi, dump, paf, gtf, fa, str(tmp_path / "o.tsv")])
    r = read_rows(tmp_path / "o.tsv")
    assert len(r) == 1 and r[0]["verdict_1"] == "AGREES" and float(r[0]["dA"]) >= 1 and r[0]["canonical_lost"] == "1"

def test_site_at_contig_end_is_unlocated(tmp_path):
    fa, psi, dump, paf, gtf = build_two_copies(tmp_path, s_copy_offset=5)      # S's acceptor 5 bases from the contig start
    main(["score", psi, dump, paf, gtf, fa, str(tmp_path / "o.tsv")])
    assert read_rows(tmp_path / "o.tsv")[0]["verdict_1"] == "UNLOCATED"
```

`build_two_copies` (test helper, written in this step): writes a FASTA with two contigs built from one random 3-exon
template with canonical sites, the four tables with hand-set counts, and a PAF of one `=` block.

- [ ] **Step 2: Run, expect ImportError**; also `python3 bench/crosscopy/sites.py --selftest` → fails (no file)
- [ ] **Step 3: Implement `sites.py`, then `site_cause.py`**
- [ ] **Step 4: Tests and self-test pass**
- [ ] **Step 5: Development-sample run** (login node is enough: one genome fetch per contig)
  `python3 bench/crosscopy/site_cause.py score NEW/or.fam.psi.tsv NEW/or.psi.introns.tsv NEW/or.fam.loci.gated.paf NEW/or.gtf $W/genomes/GGO/GGO.fa DEV/or.sites.tsv`.
  Read the UNLOCATED share and the verdict counts; fix defects found here; this is the only sample where the output may
  be inspected before the freeze. Note the counts in the progress note.
- [ ] **Step 6: Stage** the four files

### Task 3: read-origin scorer

**Files:**
- Create: `bench/crosscopy/origin.py`
- Test: `bench/crosscopy/test_origin.py`

**Interfaces:**
- Consumes: `load_calls`; `load_truth(reads_tsv, tie_split)`, `load_records(bam)`, `statement(recs, "ms", 10)`,
  `wilson` from `bench/famwise_eval/records.py` (a `Rec.chain` is the record's intron chain; read `famwise.record_chain`
  for its coordinate convention before comparing it with a call's 1-based intron); `PREFIX.fam.exonmap.columns.tsv` +
  `PREFIX.fam.loci.gff3` through `tools/psi_family.load_columns`.
- Produces:
  - `placed_introns(recs) -> set`: the introns carried by EVERY record of `statement(recs, "ms", 10)` (the wide band's
    `placed_w` rule of `psi_bunches.tally`).
  - `read_verdict(truth, member_exons) -> "SUPPORTED" | "AGAINST" | "OUTSIDE" | None` for a read placed at copy X:
    kind `one` with its origin locus on X's contig and strand sharing an exon base with X's column exons = SUPPORTED;
    kind `one` elsewhere or kind `several` = AGAINST; kind `outside` = OUTSIDE; any other kind = None.
  - CLI `origin.py loci FAM.psi.tsv FAM.loci.gff3 OUT.bed`: the member loci of LARGE and SAME calls (BED, family in
    column 4).
  - CLI `origin.py score FAM.psi.tsv COLUMNS.tsv FAM.loci.gff3 PSI.introns.tsv SEL.bam TRUTH.tsv TIE_SPLIT.tsv OUT_PREFIX`:
    `OUT_PREFIX.reads.tsv` (read, call, copy role U / S, intron, kind, verdict), `OUT_PREFIX.calls.tsv` (per LARGE call:
    SUPPORTED reads in the bunch at U and S, origin-only PSI at U and S, origin gap, HOLDS / FAILS / UNJUDGED), and on
    stdout step-0 counts first (judged reads, judged calls), then C1 (pooled, per kind, per role, per family; OUTSIDE
    share beside) and C2.
- Rules (spec section 4): a call's bunch at X = the dump's introns sharing the anchor site of X's key junction with
  `max_w >= 3`, plus the key junction. Origin-only PSI at X = SUPPORTED reads placed on the key junction over SUPPORTED
  reads placed on any intron of the bunch. HOLDS = origin gap (U minus S) >= half the call's gap.

- [ ] **Step 1: Write the failing tests**

```python
def test_read_verdicts():
    ex = [("c1", "+", 101, 200), ("c1", "+", 301, 400)]                         # copy X's column exons
    one_here = T(kind="one", origin=frozenset({0}), loci=[("c1", 150, 320, "+", 900)])
    one_else = T(kind="one", origin=frozenset({1}), loci=[("c1", 150, 320, "+", 900), ("c2", 10, 500, "+", 890)])
    antisense = T(kind="one", origin=frozenset({0}), loci=[("c1", 150, 320, "-", 900)])
    assert read_verdict(one_here, ex) == "SUPPORTED" and read_verdict(one_else, ex) == "AGAINST"
    assert read_verdict(antisense, ex) == "AGAINST"
    assert read_verdict(T(kind="several", origin=frozenset({0, 1}), loci=one_else.loci), ex) == "AGAINST"
    assert read_verdict(T(kind="outside", origin=frozenset(), loci=one_here.loci), ex) == "OUTSIDE"
    assert read_verdict(T(kind="unjudged", origin=frozenset(), loci=one_here.loci), ex) is None

def test_placed_needs_every_record_in_the_band():
    a = R(chain=((200, 300),), ms=900, primary=True); b = R(chain=((200, 300), (400, 500)), ms=895, primary=False)
    far = R(chain=(), ms=700, primary=False)
    assert placed_introns([a, b, far]) == {(200, 300)}

def test_call_holds_on_origin_reads_only(tmp_path):
    # U: 30 SUPPORTED reads, 15 on the key junction; S: 30 SUPPORTED, 0 on it; plus 10 OUTSIDE reads on the junction at S
    out = score_fixture(tmp_path, u_key=15, u_other=15, s_key=0, s_other=30, s_outside_key=10)
    assert out.calls[0]["verdict"] == "HOLDS" and out.c1 == (60, 60) and out.outside == 10
    out = score_fixture(tmp_path, u_key=15, u_other=15, s_key=14, s_other=16)   # the stated gap 0.5 shrinks to 0.03
    assert out.calls[0]["verdict"] == "FAILS"
    assert score_fixture(tmp_path, u_key=9, u_other=9, s_key=0, s_other=30).calls[0]["verdict"] == "UNJUDGED"
```

`T`, `R`, `score_fixture` are test helpers written in this step (`T` / `R` fill the other fields of `records.Truth` /
`records.Rec` with defaults; `score_fixture` builds the tables and an unsorted BAM with pysam and returns the parsed
result of `origin.score`).

- [ ] **Step 2: Run, expect ImportError**
- [ ] **Step 3: Implement `origin.py`**
- [ ] **Step 4: Tests pass**
- [ ] **Step 5: Stage** the two files

### Task 4: tool table

**Files:**
- Create: `bench/crosscopy/tool_table.py`
- Test: `bench/crosscopy/test_tool_table.py`

**Interfaces:**
- Consumes: `load_calls`, `unit`; `site_cause.py score`'s OUT.tsv for the verdict split (optional argument).
- Produces: CLI `tool_table.py --dev DEV.calls.tsv [--sites GENOME=SITES.tsv ...] GENOME:SAMPLE:FAM.psi.tsv:NAME=GTF[,NAME=GTF...] [...]`;
  `emits(gtf_introns, intron) -> bool` with `gtf_introns(path) -> set[(chrom, strand, start, end)]` (opens `.gz`).
  Output: per arm, over COPY-SPECIFIC units of the samples that have the arm: units, share emitted at U, share emitted
  at S; per pair of arms on their shared samples the reading "behind" (share at U lower by more than 0.05 on >= 30
  units) or "level"; the same split by A's verdict.

- [ ] **Step 1: Write the failing test**

```python
def test_shares_and_reading(tmp_path):
    # 40 copy-specific units; arm x emits 36 at U and 0 at S, arm y emits 30 at U and 4 at S
    t = table_fixture(tmp_path, units=40, x=(36, 0), y=(30, 4))
    assert t["x"] == (40, 0.9, 0.0) and t["y"] == (40, 0.75, 0.1) and t.reading("y", "x") == "behind" and t.reading("x", "y") == "level"
    assert table_fixture(tmp_path, units=29, x=(29, 0), y=(0, 0)).reading("y", "x") == "level"   # under 30 units
```

- [ ] **Step 2: Run, expect ImportError**
- [ ] **Step 3: Implement**
- [ ] **Step 4: Tests pass; development-sample run** with the OR6737 arms (`$W/artifact_refresh_2026-11-16/g.def/or.gtf`,
  `or.final.gtf`, `gc/GGO.stringtie.gtf`, `gc/GGO.flair.isoforms.gtf`, `gc/GGO.isoquant.gtf`, `$W/r1641_bambu/g.dn.gtf`)
  to see the table print; numbers go in the progress note as development counts only.
- [ ] **Step 5: Stage** the two files

### Task 5: pre-registration and freeze

**Files:**
- Create: `docs/PREREG_crosscopy_splicing_2026-11-18.md`
- Modify: `docs/NEGATIVE_RESULTS_REGISTER.md` (three rows, added at scoring time in Task 7), `docs/INDEX.md`

- [ ] **Step 1: Claim three row numbers** (A, C, tool table): `ls $W | tail`, `squeue -u $USER`,
  `tail -5 docs/NEGATIVE_RESULTS_REGISTER.md`, then `mkdir $W/rNNNN_crosscopy_sites $W/rNNNN_crosscopy_origin $W/rNNNN_crosscopy_tools`.
- [ ] **Step 2: Write the pre-registration**, in the form of `docs/PREREG_famwise_origin_2026-11-18.md`: Why; Rule,
  fixed now (the Global Constraints above verbatim, the call rules of spec section 1, A's and C's rules, the sample
  table of Task 6, the development counts from Tasks 1, 2, 4 labelled as seen); Readings, fixed now (A1, A2, C1, C2,
  the tool table's reading, every "undecidable" floor); what each outcome will be said to mean (spec "Risks"). State
  that scripts are the snapshot of Step 4 and that OR6737 counts in no bar.
- [ ] **Step 3: Show the pre-registration to the user and ask for the commit** (the bars are theirs; working rule).
  No job of Task 6 or 7 is submitted before this commit exists.
- [ ] **Step 4: Freeze**: `cp -r ~/Rustle/tools $W/rNNNN_crosscopy_sites/snap/tools` and `bench` beside it; every job
  below runs the snapshot. Regenerate `docs/INDEX.md` through a temp file with
  `RUSTLE_MEMORY_DIR=~/.claude/projects/-home-hpc-v138da-v138da13/memory python3 tools/docs_index.py`.

### Task 6: the report on the eleven judging samples

**Files:** job scripts in `$W/rNNNN_crosscopy_sites/` only.

Samples (run directory = links to `$W/artifact_refresh_2026-11-16/apes/SAMPLE/def/x.gtf` and `x.molecules.tsv*`):

| genome | samples | BAM (`$W/reads/`) | FASTA |
|---|---|---|---|
| GGO | GGO_KB3781_fib_mm | `GGO/GGO_KB3781_fib_mm.bam` | `$W/genomes/GGO/GGO.fa` |
| CHM13 | A119b, HG002_kinnex.t2t | `HSA/A119b.t2t.bam`, `HSA/HG002_kinnex.t2t.bam` | `~/data/chm13/chm13v2.0.fa` |
| PTR, PPA, PPY, PAB | `TAG_mm` and the cell line of each (`PTR_AG18354_lcl_mm`, `PPA_PR00251_fib_mm`, `PPY_AG05252_fib_mm`, `PAB_AG06213_fib_mm`) | `TAG/STEM_mm.bam` | `$W/genomes/TAG/TAG.fa` |

- [ ] **Step 1: One job script** `run.sbatch STEM GENOME BAM FASTA` (singlenode, 72 cores, 4 h): driver stages
  `families`, `exonmap`, `psi` from the snapshot on the linked run directory; human samples source
  `$W/r1662_edgemiss/env.sh` first (the families defaults with the SEDEF BED and exclusion sidecar). Check before
  submitting that every BAM path exists and that `x.molecules.tsv` is the score table of that BAM.
- [ ] **Step 2: Submit the eleven jobs**; wait with one background `until squeue` loop, not by polling.
- [ ] **Step 3: Per sample, verify** exit 0, `x.fam.psi.tsv` present with 31 columns, the two log lines. A failed
  `families` on an ape (missing sidecar or environment) is diagnosed and rerun; a sample that cannot run is reported
  as missing, not dropped silently.

### Task 7: scoring

- [ ] **Step 1: A.** `site_cause.py score` per sample (login node or one small job), then
  `site_cause.py pool --dev DEV/or.sites.tsv GGO=... CHM13=... CHM13=... PTR=... ...` → `A.txt`. Opened once.
- [ ] **Step 2: C truth**, one job per sample (KB3781; HG002) in `$W/rNNNN_crosscopy_origin/`:
  `origin.py loci` → BED; `samtools view -L BED` on the full BAM → read names → `samtools view -N names` over the full
  BAM → `sel.bam` (every record of those reads; sort temps in `$TMPDIR`); an empty read-to-family file (`tie_split.py`
  uses it only to label its printed summary); then exactly the steps of `$W/r1796_famwise_origin/hg_held.sbatch`
  from `neutral_truth.py score` on, with `NT_REPAIR=1 NT_CAP=100000000`, no `NT_SELECT_*`, the haplotype indexes
  `$W/r1776_hg002_room/{mat,pat}.mmi` (HG002) and `$W/o2contested/haps/KB3781.{mat,pat}.mmi` (KB3781), and
  `tie_split.py`. Size the wall time from that job's log (42,572 reads scored in 1,169 s).
- [ ] **Step 3: C score.** `origin.py score ...` per sample → `C.KB3781.txt`, `C.HG002.txt`. The step-0 counts are read
  first; under a floor the bar is undecidable and is written so.
- [ ] **Step 4: Tool table** over all judging samples with the arms on disk: `$W/apetools/SAMPLE/st_dn.gtf`,
  `fl_dn.isoforms.gtf`, `iq_dn/x/x.transcript_models.gtf`; bambu where a de novo GTF exists (`$W/r1641_bambu/a.dn.gtf`,
  `$W/r1766_allmiss/{ptr,ppa}.bambu.gtf`) → `T.txt`.
- [ ] **Step 5: Record.** Result section of the pre-registration (every number with its count and interval, every
  undecidable said as such, the development counts kept apart); three register rows written read-compute-write with
  the row numbers re-checked; an entry in `docs/LEDGER_O1.md`; one line in `docs/ADVISOR_QUESTIONS.md` under the
  advisor's "nothing new" objection. Stage by explicit path, show `git diff --cached --stat`, ask the user for the
  commit.
- [ ] **Step 6: Report to the user** in plain words: A1, C1, C2 with PASS / FAIL / undecidable, the tool table's
  reading, and what may and may not be said to the advisor from them. Offer the advisor page (visual first).
