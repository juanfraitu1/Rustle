# Family-wise layer, part 2: measurement against read origin Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Score the layer and three baseline assemblies against the origin of the reads (E1), bound what any provenance
method could change in an assembly on real and planted reads (E2), and measure what a re-drawn tie-break changes (E3).

**Architecture:** One bench package, `bench/famwise_eval/`. A shared reader turns a truth table and a small BAM into
per-molecule records with both scores. The scorer, the input builder for the three assembly inputs, the region chain
counter, the simulator map and the re-toss comparison are separate scripts on that reader. Jobs live in the workspace.

**Tech Stack:** Python 3 with pysam, pytest for the package, `bench/famsim`, minimap2 2.30, samtools, StringTie 3.0.3,
the driver's `assemble` and `famwise` stages, Slurm.

**Spec:** `docs/superpowers/specs/2026-11-18-family-wise-assembly-design.md` (sections "Scores" and "Measurement").
Part 1 (`docs/superpowers/plans/2026-11-18-famwise-stage.md`) must be finished first.

## Global Constraints
- **Nothing that decides a bar runs before its pre-registration is written and the user has set the bar** (Task 1 is a
  gate). KB3781 develops; the held-out quarter of HG002 reads is opened once, by Task 5.
- `ms` first, AS beside it: every band, stratum and count is computed on both scores and printed with the score named;
  bars are on `ms`. The truth tables are built on AS, so the scorer re-reads both scores from the BAM.
- Truth: `/anvme/workspace/v138da13-rustle/r1785_truth_repair/{kb,hg}/reads.tsv`,
  `r1787_winnowmap_orangutan/pab/reads.tsv`; header `read	stratum	lead	status	picks	origin_loci	loci`; `loci` items
  `i:chrom:start-end:strand:ASn` (0-based half-open spans), at most 12; an index >= 12 in `picks` or `origin_loci` is
  unjudgeable and counted.
- Limits printed in every result file: annotation families only; clear-lead reads sampled (30,000); 12 loci per read;
  kept placements at AS >= 0.8 x best; HG002 is a 25% read sample on another individual's reference; `outside` reads
  give the upper bound on failures of assumption A.
- Baselines are trusted baselines, never "competitors". Read attachment uses one rule for every tool: exact equality of
  (contig, strand, intron chain) between a record and a transcript. No tool's own read map.
- Wording: "reads with primary or secondary alignments", never "primary / secondary reads".
- Paths: `W=/anvme/workspace/v138da13-rustle`; tools on `PATH` from
  `$HOME/.pixi/detached-envs/rustle-7434162673922408322/envs/default/bin` and `$HOME/.pixi/bin` (StringTie); assembler
  binaries `--bin $HOME/rustle_target/release`. No path is hard-coded in a script: all are arguments.
- Jobs: `-p singlenode -t 01:30:00`; sort temps and indexes in `$TMPDIR`, products to the row's workspace directory.
  If `sbatch` is denied by the storage guard, stop and report.
- No commit unless the user asks. Register rows are claimed by `mkdir $W/rNNNN_name` (next free was 1793 on
  2026-11-18; re-check) and written only with the user's go.

## Review Focus
1. A record inside the `ms` band that the truth never judged (AS < 0.8 x best, or beyond its 12 loci): the molecule has
   no verdict for that statement; it is counted and listed, never scored as wrong or right. (Task 3.)
2. A read near-tied on one score only (KB3781: 116 of 1,321 AS near ties are not within 10 on `ms`): it appears in the
   per-score tables under its own stratum and in `one_score_only.tsv`. (Task 3.)
3. A truth read missing from the BAM, or a BAM read without a truth row: counted per file, no crash. (Task 2.)
4. An origin read with no record at its origin copy when the oracle input is built: counted `no_home`, kept at its
   primary alignment. (Task 6.)
5. A tie-split table made under another truth rule (4,249 against 4,231 ties): join on read name, unmatched ties stay
   `tie_unsplit` and are reported, not dropped. (Task 2.)

---

### Task 1: pre-registration (GATE: the user sets the bars)

**Files:**
- Create: `docs/PREREG_famwise_origin_<date>.md`

- [ ] **Step 1: Write the document** in the house form of `docs/PREREG_tie_split_2026-11-18.md`: first paragraph
  "Written before any count of this row" with workspace and inputs; `## Why`; `## Rule, fixed now`; `## Readings, fixed
  now`; `## Prediction`; `## Limits stated beforehand`. It registers three rows:
  - **Row a (E1)**, bars B1 and B2. B1, on the held-out HG002 quarter, molecules near-tied on `ms` (two best records
    within 10), truth kinds `one` and `several` pooled, bands on `ms`: the exact-band statement equals the truth set more often than the primary alignment alone, Wilson 95%
    intervals disjoint, and in more than half of the families with >= 20 judged molecules. B2, on the same reads:
    (i) the share of PHANTOM + SHARED claims is higher among `tied_only` nodes than among `certified` ones, intervals
    disjoint; (ii) `certified` + `best_only` nodes hold >= 90% of the default GTF's TRUE claims; fewer than 30 judged
    `tied_only` nodes = undecidable, said so.
  - **Row b (E2, real reads)**, bar B3, on the regions a tied or moved read touches: advisor's expectation holds if
    R >= 2 for either assembler; ours if R <= 1.1 for both; in between neither is claimed.
  - **Row c (E2, planted)**, bar B4: several-fold cells (median over 5 seeds of recovered pairs, O over P, >= 2, or P
    none and O some) lie only at <= 2 differences per transcript or in the array row, and at least one exists.
  - E3 and the range table are descriptive, stated as such.
  - "Seen already": the AS-proxy table of the design (three samples), the direction of B3 on KB3781 (937 of 100,810
    reads differ between O and P), the tied-loci ratios x1.10-1.22.
- [ ] **Step 2: STOP. Show the user the four bars as drafted and ask for the values.** Do not run Tasks 4-9 before the
  answer is in the document. (Answered 2026-11-18: the user accepts B1-B4 as drafted above, "yes these are good"; the
  document still has to be written with them before any scoring run.)

---

### Task 2: the shared reader (`bench/famwise_eval/records.py`)

**Files:**
- Create: `bench/famwise_eval/records.py`, `bench/famwise_eval/test_records.py`

**Interfaces:**
- Consumes: `tools/famwise.py`: `record_chain`.
- Produces:
  - `Truth = namedtuple("Truth", "stratum status kind origin loci")`; `loci` = `[(chrom, start, end, strand, AS)]`
    (0-based half-open); `origin` = frozenset of indexes into `loci`; `kind` one of `one` (status `decidable`),
    `several` (status `tie`, split class `SEVERAL COPIES`), `one_place`, `tie_unsplit` (status `tie`, no split row),
    `outside`, `unjudged` (any other status, or an origin index >= `len(loci)`).
  - `load_truth(reads_tsv: str, tie_split: str | None = None) -> dict[str, Truth]`
  - `Rec = namedtuple("Rec", "chrom start end strand chain ms AS primary")`
  - `load_records(bam: str) -> dict[str, list[Rec]]`: every mapped, non-supplementary record, read with
    `until_eof=True` (the selected-reads BAMs are unsorted and unindexed).
  - `locus_index(rec: Rec, t: Truth) -> int | None`: the truth locus on the record's contig and strand that overlaps
    it; `None` = the truth never judged this record.
  - `statement(recs: list[Rec], score: str, window: float | None) -> list[Rec]`: `window=None` = the primary record
    alone; else the primary plus every record whose `score` (`"ms"` or `"AS"`) is within `window` of the molecule's best.

- [ ] **Step 1: Failing tests** (`test_records.py`, pytest, truth rows and a tiny unsorted BAM written in `tmp_path`):
  `test_kinds` (one row per kind above, including an origin index 12 with 12 loci -> `unjudged`, and a tie absent from
  the split table -> `tie_unsplit`); `test_locus_index` (overlap on the right strand -> index, other strand -> `None`);
  `test_statement_both_scores` (records ms 200 / 195 / 185, AS 190 / 196 / 150: `statement(r, "ms", 10)` has two
  records, `statement(r, "AS", 10)` has two different ones, `statement(r, "ms", None)` is the primary).
- [ ] **Step 2: Run** `PYTHONDONTWRITEBYTECODE=1 $PY -m pytest bench/famwise_eval -q -p no:cacheprovider` -> FAIL
  (module missing).
- [ ] **Step 3: Implement.** - [ ] **Step 4: Run** -> 3 passed.

---

### Task 3: the origin scorer (`bench/famwise_eval/origin_score.py`)

**Files:**
- Create: `bench/famwise_eval/origin_score.py`, tests appended to `bench/famwise_eval/test_records.py`

**Interfaces:**
- Consumes: Task 2; `famwise.load_nodes` is NOT used (a baseline's `gene_id` is not a locus): claims are keyed by
  `(contig, strand, chain)` from each arm's GTF (`gtf_names.chain` on its exons, multi-exon transcripts).
- Produces: CLI
  `origin_score.py --truth READS.tsv [--tie-split T.tsv] [--families READ_FAMILY.tsv] --bam SEL.bam --arm NAME=GTF [...] [--tx FW.tx.tsv --tx-arm NAME] --out DIR`
  writing, each with a `# limits:` first line:
  - `molecule.tsv`: `score	window	near_tied	kind	n	equals	sound	no_verdict` for score in ms, AS; window in primary,
    0, 5, 10, 20; `near_tied` yes / no = the molecule's two best records are within 10 on THAT score (the truth's own
    stratum column is on AS and is not used for this). A molecule is `no_verdict` for a statement holding a record
    with `locus_index` `None`. First comment lines: reads in the truth, in the BAM, in both.
  - `molecule_family.tsv`: the same for `ms`, windows primary and 0, per family (needs `--families`).
  - `one_score_only.tsv`: `read	near_tied_on	lead_ms	lead_AS` for molecules within 10 on one score and not the other.
  - `claims.tsv`: `arm	status	claims	TRUE	SHARED	PHANTOM	UNJUDGED`; `status` = the layer's status for the arm named
    by `--tx-arm` (joined through the `transcripts` column of `FW.tx.tsv`, band `wide`), `.` for the other arms.
    Support of a transcript = judged molecules (kinds `one`, `several`) with a record within 10 `ms` of their best (or
    primary) whose `(contig, strand, chain)` equals the transcript's. TRUE: a supporter's origin set is exactly that
    record's locus index; else SHARED: a supporter's origin set contains it; else PHANTOM; no supporter: UNJUDGED.
  - `range.tsv`: `band	nodes	held	a_holds	a_n`: on kind `one` molecules only, nodes where (molecules with every band
    record at the node) <= (molecules whose origin record is at the node) <= (molecules with any band record there);
    and the molecules whose origin record is in their band.
  - `bars.tsv`: the quantities B1 and B2 name, with Wilson 95% intervals, and `undecidable` where the registered
    minimum is not met.

- [ ] **Step 1: Failing tests** on hand-made inputs: `test_molecule_table` (one `one` read whose primary is off its
  origin and whose exact band holds both -> primary: equals 0, sound 0; window 0: equals 0, sound 1);
  `test_claims` (three transcripts: supported by a one-origin read from there -> TRUE; supported only by a `several`
  read -> SHARED; supported only by a read from elsewhere -> PHANTOM; a fourth without support -> UNJUDGED; the same
  transcript in two arms gets the same verdict); `test_no_verdict_is_not_wrong`; `test_one_score_only` (a molecule
  within 10 on AS and 15 apart on `ms` is listed, and sits under `near_tied` yes for AS, no for `ms`);
  `test_unmatched_reads_counted` (a truth read absent from the BAM and a BAM read absent from the truth: both counted
  in the comment lines, exit 0).
- [ ] **Step 2: Run** -> FAIL. - [ ] **Step 3: Implement.** - [ ] **Step 4: Run** -> all pass.

---

### Task 4: development run on KB3781 (after the gate)

**Files:**
- Create (workspace): `$W/rNNNN_famwise_origin/kb.sbatch` and outputs under `kb/`

- [ ] **Step 1: Job.** From `$W/r1775_eichler_baseline/mm.sel.bam`: `samtools sort` + `index` into the row directory
  (temps in `$TMPDIR`); `$HOME/rustle_target/release/as_table --bam kb/sel.bam --out kb/sel.as_table.tsv --threads 32`;
  the layer on it with the default run's GTF
  (`eq_classes.py $W/artifact_refresh_2026-11-16/apes/GGO_KB3781_fib_mm/def/x.families.gtf kb/sel.bam kb/sel.as_table.tsv kb/eq --threads 32 --famwise kb/fw`,
  then `famwise.py GTF kb/fw --eq-loci kb/eq.loci.tsv`); then `origin_score.py` with the truth
  `$W/r1785_truth_repair/kb/reads.tsv`, arms `rustle=`the default `x.families.gtf`, `stringtie=$W/apetools/GGO_KB3781_fib_mm/st_dn.gtf`,
  `flair=.../fl_dn.isoforms.gtf`, `isoquant=.../iq_dn/x/x.transcript_models.gtf`, `--tx kb/fw.tx.tsv --tx-arm rustle`.
- [ ] **Step 2: Check the reader against known counts** (not a bar): among molecules near-tied on AS, kind `one` has
  45 and `tie_unsplit` 786; `molecule.tsv` on AS reproduces the design's preview (kind `one`: primary equals 40 of 45;
  kind `tie_unsplit`, window 10: 778 of 786). Then write the `ms` version of that preview table into the design's
  "Preview" section, under the AS one, as the spec's "Scores" section asks.
- [ ] **Step 3: A second scoring** with `--tx` from the whole-sample stage run of part 1 (`run/kb.famwise.tx.tsv`):
  how many nodes change status between the selected reads and all reads. Reported; the registered scoring uses the
  selected reads.
- [ ] **Step 4: Freeze.** Any change to the scorer after this step is listed in the pre-registration as an amendment
  before Task 5.

---

### Task 5: the held-out HG002 quarter (opened once)

**Files:**
- Modify: `bench/o2_truth/neutral_truth.py:24` and `:40` (read selection)
- Create: `bench/o2_truth/tie_split.py` (port of `$W/r1791_tie_split/tie_split.py`)
- Create (workspace): `$W/rNNNN_famwise_origin/hg.sbatch`

**Interfaces:**
- Produces: env `NT_SELECT_OFFSET` (default 0): a read is selected when
  `(int(sha256(name)[:8], 16) - OFFSET) % 100 < PCT`. With `NT_SELECT_PCT=25 NT_SELECT_OFFSET=25` the selection is the
  residues 25-49, disjoint from every read looked at so far.
  `tie_split.py READS.tsv READS.fa PLACEMENTS.fa MAT.mmi PAT.mmi READ_FAMILY.tsv OUT.tsv`: the workspace script
  unchanged except that it imports `neutral_truth` from its own directory and drops `assert not NT.REPAIR` (the truth
  it now reads is made with `NT_REPAIR=1`); a tie for which its own locus test finds fewer than two origin loci is
  written with class `UNSPLIT` and counted, instead of failing the second assert (the reader then keeps it as
  `tie_unsplit`).

- [ ] **Step 1: Test the offset** (`bench/o2_truth/test_select_offset.py`): over 2,000 synthetic names, offsets 0 and
  25 at PCT 25 select disjoint sets of about a quarter each, and offset 0 selects exactly what the old expression did.
- [ ] **Step 2: Job**, copying `$W/r1776_hg002_room/room.sbatch` (select 20 min from the full BAM, minimap2 6 min, score
  6 min) with `NT_SELECT_PCT=25 NT_SELECT_OFFSET=25 NT_REPAIR=1`, then the read-family table
  (`$W/r1790_family_ties/family_ties.py`), then `tie_split.py`, then sort / index / `as_table` / layer / `origin_score.py`
  exactly as Task 4, with `--tie-split`, `--families` and the HG002 arms
  (`$W/artifact_refresh_2026-11-16/apes/HG002_kinnex.t2t/def/x.families.gtf`, `$W/apetools/HG002_kinnex.t2t/...`).
- [ ] **Step 3: Read `bars.tsv` against the registered B1 and B2.** Write the result, pass or fail, into the
  pre-registration under `## Result` with the limits line. Orangutan (`r1787`) is the same job without the split and
  is cut first.

---

### Task 6: three assembly inputs from one alignment (`bench/famwise_eval/origin_bam.py`)

**Files:**
- Create: `bench/famwise_eval/origin_bam.py`, tests in `bench/famwise_eval/test_origin_bam.py`

**Interfaces:**
- Consumes: Task 2.
- Produces: CLI `origin_bam.py --truth READS.tsv [--tie-split T.tsv] --bam SEL.bam --out DIR` writing sorted, indexed
  `DIR/D.bam` and `DIR/O.bam`, `DIR/regions.bed` and `DIR/inputs.tsv`.
  - A read with truth kind `one`, `several` or `tie_unsplit`: `D` keeps all its records; `O` keeps the records whose
    `locus_index` is in its origin set. If none: `no_home`, both keep the primary only. If the kept set has no primary
    record, the kept record with the highest `ms` gets flag 0x100 cleared.
  - Any other read: `D` and `O` keep its primary record only.
  - `regions.bed`: merged truth spans of every read that is tied (`several`, `tie_unsplit`) or moved (`one` with the
    primary off its origin).
  - `inputs.tsv`: reads per case (kept as is, restricted, moved, `no_home`, primary-only), and for `D` the secondary
    records that tie the best `ms` exactly but lie below 0.98 x the best AS (what the default seeding band drops).
  - The three inputs: P = `D.bam` with `--seed-pool primary` (StringTie: `samtools view -b -F 256 D.bam`); D = `D.bam`
    with defaults; O = `O.bam` with `--seed-pool all` (StringTie: `O.bam`).

- [ ] **Step 1: Failing tests** on a tiny BAM: a moved read (secondary promoted, old primary gone from `O`), a
  `several` read (two records kept, one primary), a read without truth (one record in both), a `no_home` read, and
  `regions.bed` holding the first two only.
- [ ] **Step 2: Run** -> FAIL. - [ ] **Step 3: Implement.** - [ ] **Step 4: Run** -> pass.
- [ ] **Step 5: Settle one question before any run:** on one famsim family (Task 8's cell "0 differences, 1:1"), does
  `assemble --seed-pool all` seed a promoted record with MAPQ 0, and does StringTie use it? If either ignores it, set
  MAPQ 60 on promoted records, say so in `inputs.tsv`, and re-run the test.

---

### Task 7: E2 on real reads (`bench/famwise_eval/region_chains.py` and the job)

**Files:**
- Create: `bench/famwise_eval/region_chains.py`, test in `bench/famwise_eval/test_origin_bam.py`
- Create (workspace): `$W/rNNNN_famwise_oracle/kb.sbatch`

**Interfaces:**
- Consumes: `bench/famcmp/chain_prec.py`: `chains(path) -> set[(contig, strand, chain)]`; Task 6.
- Produces: CLI `region_chains.py REF.gtf REGIONS.bed P=GTF D=GTF O=GTF` printing
  `scope	input	chains	refseq_exact` for scope `regions` and `genome`, then `R	scope	value` (RefSeq-exact, O over P) and
  `D_not_O	scope	chains	refseq_exact`. A chain is in `regions` when its span overlaps a BED interval on its contig.

- [ ] **Step 1: Failing test** (three tiny GTFs and a BED: counts, R and `D_not_O` by hand). - [ ] **Step 2-4:**
  implement, pass.
- [ ] **Step 5: Job** on KB3781: `origin_bam.py`; `rustle_pipeline.sh assemble --bam X --fasta $W/genomes/GGO/GGO.fa
  --out DIR/x --threads 36 --bin ...` three times (flags of Task 6); `stringtie -L -p 8 -o DIR/st.gtf X` three times;
  `region_chains.py $W/apetools/GGO.refseq.gtf regions.bed ...` per assembler.
- [ ] **Step 6: Read R on the regions against the registered B3**; write the result into the pre-registration. HG002
  is the same job on the held-out truth of Task 5 and is cut before KB3781.

---

### Task 8: the planted map (`bench/famwise_eval/regime_map.py`)

**Files:**
- Create: `bench/famwise_eval/regime_map.py`, test `bench/famwise_eval/test_regime_map.py`

**Interfaces:**
- Consumes: `bench/famsim` (`make`, `verify`, `align` through its CLI, pattern of `bench/o2_contested/verdict_sim.py:36-43`),
  `chromosome.load_planted`, `Planted.junctions(skip=())`, `Planted.isoform_list()`, `evaluate.gtf_loci`.
- Produces: CLI `regime_map.py --out DIR --bin RELEASE_DIR [--seeds 1,2,3,4,5] [--only CELL,...]` writing
  `DIR/map.tsv`: `cell	differences	expression	reads	seed	input	planted_pairs	recovered	silent_chains	tied_ms	tied_AS`
  and `DIR/cells.tsv` (median of recovered O over P per cell, `several_fold` yes / no).
  - Template: synthetic, exons `[300, 200, 400, 250, 150, 700]` (2,000 bp), introns `[300, 300, 300, 300, 300]`, saved
    once and reused as `{"source": "file"}` so that seeds change reads and SNP positions only. Background random,
    200,000 bp.
  - Two copies A and B, each with `isoforms: [{"skip": [3], "weight": 0.3}]`. B carries
    `{"op": "snp", "rate": k / 2000, "region": "exons"}` for k in 0, 1, 2, 4, 8, 20 (no op at k = 0).
    Expression A / B: 200 / 200, 200 / 20, 200 / 0. Reads: full length, or `trunc5_frac 1.0, trunc5_max 0.7`
    (3' fragments). 36 cells. Array row: 8 copies without ops, 25 reads each, both read kinds: 2 cells.
  - Inputs per cell from the one `reads.bam`: P (`--seed-pool primary`), D (default), O (records overlapping the read's
    planted copy, name field `copy|isoform|i`, spans from `copies.tsv`, promoted as in Task 6; `--seed-pool all`),
    each `assemble` with its own `--out`; then the stage `famwise` on D.
  - `recovered` = planted (isoform, copy) pairs with expression > 0 whose exact junction chain is emitted on the copy's
    span; `silent_chains` = multi-exon chains emitted on a copy with expression 0; `tied_ms` / `tied_AS` = share of
    reads whose two best records tie exactly on each score (score-table columns 7-8 and 2-3).
  - `verify` FAILs of `pair_identity_genomic` at exon-only SNPs are expected and listed; any other FAIL stops the cell.

- [ ] **Step 1: Failing test:** `spec_for(cell)` gives the JSON above for three cells (k = 0 has no op; k = 20 rate
  0.01; array has 8 copies); `several_fold([..])` on hand-made rows (P 0 and O 2 -> yes; medians 3 / 2 -> no).
- [ ] **Step 2-4:** implement, pass.
- [ ] **Step 5: One cell end to end** on the login node (`--only k0_1to1_full --seeds 1`): it runs, the three GTFs
  differ or not as printed, wall time noted. Then the whole map as one job.
- [ ] **Step 6: Read `cells.tsv` against the registered B4.** Draw the map (load the `dataviz` skill first): one panel
  per read kind, differences on x, expression ratio on y, the real two-member pairs placed by
  k = (1 - `identity`) x 2000 from `$W/r1773_sibling_reverse/{kb,or}.siblings.tsv` (columns `family_size`, `identity`,
  `own_primaries`).

---

### Task 9: the re-toss (`bench/famwise_eval/retoss.py`)

**Files:**
- Create: `bench/famwise_eval/retoss.py`, test in `bench/famwise_eval/test_records.py`
- Create (workspace): `$W/rNNNN_famwise_oracle/retoss.sbatch`

**Interfaces:**
- Produces: CLI `retoss.py records A.bam B.bam` (records keyed `(read, contig, start, strand, cigar, ms)`: counts in
  both, only A, only B, and reads whose primary record differs) and
  `retoss.py claims NAME=GTF[,GTF...] [...]` (per arm, chains present in every GTF of the arm against chains in any).

- [ ] **Step 1: Failing test** (two tiny BAMs differing in one primary flag and one record). - [ ] **Step 2-4.**
- [ ] **Step 5: Job** on `$W/r1775_eichler_baseline/kb.sel.fa`: the BAM's own command
  (`minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes -t 72 REF kb.sel.fa`) with `--seed 11`, 12, 13, 14;
  sort, index. Control: `retoss.py records` seed 11 against `mm.sel.bam` must report no difference; if it does, stop
  and report (the control failed, nothing else is read). Then per seed: `assemble` (defaults), `stringtie -L`, and the
  layer with the seed-11 GTF on each seed's BAM.
- [ ] **Step 6: Report** per arm the share of claims not in every seed, the reads whose primary moved, and the rows of
  `fw.tx.tsv` (`lo`, `hi`, `status`) that differ between seeds, each on `ms` with the AS count beside. Descriptive.

---

### Task 10: results into the project record

- [ ] **Step 1:** Each pre-registered row gets its `## Result` text and limits. Ask the user before writing register
  rows (`docs/NEGATIVE_RESULTS_REGISTER.md`, next free number re-checked) and before any commit.
- [ ] **Step 2:** Offer the advisor page (visual first: the map of Task 8 with the real pairs, the claims table of
  Task 5, R of Task 7). It is published only on the user's request.

## Cut order if time runs short
Orangutan; IsoQuant and FLAIR in Task 9; the expression axis of Task 8; HG002 in Task 7. Never cut: Task 1, Task 4,
Task 5, Task 7 on KB3781 with both assemblers, the difference ladder with the array row.
