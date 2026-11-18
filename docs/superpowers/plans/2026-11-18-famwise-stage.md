# Family-wise layer, part 1: the opt-in stage `famwise` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** An opt-in driver stage `famwise` that joins the transcripts of an assembly into family transcripts through tied
molecules and writes per-node count ranges and statuses, without changing any existing product.

**Architecture:** `tools/eq_classes.py` gets one switch that, in its single BAM pass, also writes every tied molecule's
records (locus, node, both scores) and the node of every untied molecule. A new BAM-free `tools/famwise.py` turns those
two tables into family transcripts, ranges, statuses and GTF attributes. The driver gets `stage_famwise` = `eqclasses`
with the switch, then `famwise.py`.

**Tech Stack:** Python 3 with pysam (the pixi env), pytest, bash driver, one Rust dispatch test. No new dependency, no Rust.

**Spec:** `docs/superpowers/specs/2026-11-18-family-wise-assembly-design.md` (sections "The layer", "Scores", "Shape of the
implementation"). Part 2 (measurement E1-E3) is `docs/superpowers/plans/2026-11-18-famwise-evaluation.md`.

## Global Constraints
- Opt-in only: `famwise` is named explicitly and never enters the plan of `all`.
- Every existing product stays byte-identical. `PREFIX.eq.loci.tsv` and `PREFIX.eq.classes.tsv` are identical with and
  without the switch (`cmp`). `eq_classes.py` prints exactly its two summary lines to stdout, nothing after them.
- Bands are on `ms` (minimap2 ranks on it). Every record row also carries AS (user, 2026-11-18).
- Windows: exact = gap 0; wide = gap <= 10 (`RUSTLE_EQ_WINDOW`); the reads table keeps gap <= 20 (`FW_WINDOW = 20`).
  The primary record is in every band whatever its score. `--ratio` together with the switch is refused (exit 2).
- Link rule: one molecule (`--min-link 1`).
- Intron chain convention: `gtf_names.chain`: tuple of `(exon_end, next_exon_start)`, 1-based closed exons. Introns of a
  record come from `cigartuples` op 3 only (`get_blocks()` splits at every `=` / `X` / `D` on these `--eqx` BAMs).
- Python: `PY=$HOME/.pixi/detached-envs/rustle-7434162673922408322/envs/default/bin/python3`. Tests:
  `PYTHONDONTWRITEBYTECODE=1 $PY -m pytest tools/test_famwise.py tools/test_eq_classes.py -q -p no:cacheprovider`.
- No hard-coded paths in code. Outputs of runs go to `/anvme/workspace/v138da13-rustle/famwise_design_2026-11-18/`.
- Shared repository: `tools/eq_classes.py`, `tools/rustle_pipeline.sh` and `README.md` hold uncommitted edits of other
  sessions. Edit the working tree by exact anchors, never revert or reformat. Do NOT commit: the user commits on request
  (then `git add -p` own hunks only, `git diff --cached --stat` first).
- `sbatch` may be denied by the storage guard: do not bypass, report.

## Review Focus
1. A molecule whose primary record is not its best-scoring one (primary `ms` < best): the primary is in the exact band
   by rule, the best record too; the molecule is tied between them. (Task 2, `test_primary_below_best`.)
2. A tied molecule whose tied record was not scanned (`--contigs` panel, e.g. CD99 on X with Y unscanned): it must not
   certify a node. It gets the member `unseen` and its class is open. (Task 3, `test_unseen_rival`.)
3. Two transcripts with the same locus and chain (FLAIR writes 635 such): one node, both lines get the same values.
   (Task 1, `test_load_nodes`; Task 4, `test_gtf_attributes`.)
4. Nothing tied at all, or a GTF with single-exon transcripts only: header-only tables, exit 0. (Task 4, `test_empty`.)
5. `--ratio` with the switch, or `RUSTLE_EQ_RATIO` set for the stage: refused with a message, exit 2. (Tasks 3 and 5.)

---

### Task 1: nodes and record chains (`tools/famwise.py`, first part)

**Files:**
- Create: `tools/famwise.py`
- Test: `tools/test_famwise.py`

**Interfaces:**
- Produces:
  - constants `UNSPLICED = "U"`, `OTHER = "other"`, `OUTSIDE = "outside"`, `UNSEEN = "unseen"`, `FW_WINDOW = 20`
  - `record_chain(pos: int, cigartuples: list[tuple[int, int]]) -> tuple` (`pos` = pysam `reference_start`, 0-based)
  - `load_nodes(gtf: str) -> tuple[dict, dict]`: `nodes {locus: {chain: node_id}}`, `tx {node_id: [transcript_id, ...]}`.
    Locus = `gene_id`. Multi-exon transcripts only. `node_id` = the smallest `transcript_id` sharing (locus, chain).
  - `member_of(locus: str, chain: tuple, nodes: dict) -> str`: `"."` when `locus == OUTSIDE`, `UNSPLICED` when
    `chain == ()`, else `nodes.get(locus, {}).get(chain, OTHER)`.
  - A member is addressed everywhere as the pair `(locus, member)`.

- [ ] **Step 1: Write the failing tests**

```python
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import famwise as fw
from gtf_names import chain

GTF = "\n".join("\t".join(r) for r in [
    ("c1", "x", "exon", "1", "100", ".", "+", ".", 'gene_id "A"; transcript_id "A.2";'),
    ("c1", "x", "exon", "151", "250", ".", "+", ".", 'gene_id "A"; transcript_id "A.2";'),
    ("c1", "x", "exon", "1", "100", ".", "+", ".", 'gene_id "A"; transcript_id "A.1";'),
    ("c1", "x", "exon", "151", "300", ".", "+", ".", 'gene_id "A"; transcript_id "A.1";'),
    ("c1", "x", "exon", "500", "700", ".", ".", ".", 'gene_id "S"; transcript_id "S.1";'),
]) + "\n"

def test_record_chain():
    assert fw.record_chain(0, [(7, 100), (3, 50), (8, 1), (7, 99)]) == ((100, 151),) == chain([(1, 100), (151, 250)])
    assert fw.record_chain(0, [(7, 40), (2, 3), (7, 57)]) == ()          # a deletion is not an intron
    assert fw.record_chain(9, [(4, 5), (7, 10), (3, 20), (7, 10)]) == ((19, 40),)   # soft clip consumes no reference

def test_load_nodes(tmp_path):
    p = tmp_path / "a.gtf"; p.write_text(GTF)
    nodes, tx = fw.load_nodes(str(p))
    assert nodes == {"A": {((100, 151),): "A.1"}} and tx == {"A.1": ["A.1", "A.2"]}   # same chain: one node; S is single-exon
    assert fw.member_of("A", ((100, 151),), nodes) == "A.1"
    assert fw.member_of("A", ((100, 160),), nodes) == fw.OTHER
    assert fw.member_of("A", (), nodes) == fw.UNSPLICED and fw.member_of(fw.OUTSIDE, (), nodes) == "."
```

- [ ] **Step 2: Run to see them fail**

Run: `PYTHONDONTWRITEBYTECODE=1 $PY -m pytest tools/test_famwise.py -q -p no:cacheprovider`
Expected: FAIL, `ModuleNotFoundError: No module named 'famwise'`.

- [ ] **Step 3: Implement the three functions and the constants in `tools/famwise.py`**

Module docstring states the CLI of Task 4. Read the GTF with `gtf_names.opn`; exons per `(gene_id, transcript_id)`;
`chain` from `gtf_names`. Reference-consuming CIGAR ops are 0, 2, 3, 7, 8; only op 3 closes an exon.

- [ ] **Step 4: Run the tests**

Run: same command. Expected: 2 passed.

---

### Task 2: the layer on tables (`tools/famwise.py`, core; no BAM, no files)

**Files:**
- Modify: `tools/famwise.py`
- Test: `tools/test_famwise.py`

**Interfaces:**
- Consumes: constants of Task 1.
- Produces:
  - `bands(records: dict, window: float = 10) -> dict`: `records {read: [(locus, member, gap, primary), ...]}` ->
    `{"exact": {read: frozenset[(locus, member)]}, "wide": {...}}`. A record is in `exact` when `primary or gap == 0`,
    in `wide` when `primary or gap <= window`. A read with no record in a band is absent from that band.
  - `layer(sets: dict, base: Counter, repeat: frozenset = frozenset(), min_link: int = 1) -> dict` with keys
    `placed` (Counter over `(locus, member)`: `base` plus every read whose set has one member),
    `classes` (Counter over frozensets of two or more members),
    `ft` (`{(locus, node): key}` for every node that occurs; `key` = the smallest `(locus, node)` of its component),
    `inside`, `open` (Counters over `key`), `lo`, `hi` (Counters over `(locus, member)`),
    `collision` (set of keys whose component holds two nodes of one locus).
    A member is linkable when its member string is none of `UNSPLICED`, `OTHER`, `"."` and its locus is neither
    `OUTSIDE`, `UNSEEN` nor in `repeat`. Two linkable members are linked when classes holding both sum to `>= min_link`
    molecules. A class is `inside` when every member is linkable and all lie in one component; otherwise it is `open`
    for every component it touches. `lo` = `placed`; `hi` = `lo` + molecules of every class holding the member.
  - `locus_counts(sets: dict, base: Counter) -> tuple[Counter, Counter]`: per locus, molecules with every band record at
    that locus, and that plus molecules with any record there. Members of locus `UNSEEN` are dropped before counting
    (the `eq` tables never saw them); `OUTSIDE` counts as a locus but gets no entry. Must equal `placed` / `max` of
    `eq_classes`.
  - `statuses(exact: dict, wide: dict, all_nodes) -> dict[(locus, node), str]`, first match wins:
    `certified` (wide `lo` >= 1), `best_only` (exact `lo` >= 1), `tied_only` (exact `hi` >= 1), `near_only`
    (wide `hi` >= 1: only records below the molecule's best carry it), `no_band_record`.
  - `min_copies(required: set, class_loci: list[frozenset], cap: int = 200000) -> tuple[int, bool]`: fewest loci that
    contain `required` and meet every set of `class_loci`; `(n, True)` when solved exactly by combinations within `cap`
    tries, else `(len(required) + greedy count of pairwise disjoint unmet sets, False)`, the sets taken by size, then by
    sorted members (still a lower bound).

- [ ] **Step 1: Write the failing tests** (helpers at the top of the test file)

```python
import collections, random
A1, A2, B1, AO = ("A", "a1"), ("A", "a2"), ("B", "b1"), ("A", fw.OTHER)
def L(records, base=(), **kw):
    s = fw.bands(records); b = collections.Counter(dict(base))
    return fw.layer(s["exact"], b, **kw), fw.layer(s["wide"], b, **kw), s, b

def test_two_tied_emitted_chains():
    e, w, s, b = L({"r1": [("A", "a1", 0, 1), ("B", "b1", 0, 0)], "r2": [("A", "a1", 0, 0), ("B", "b1", 0, 1)]})
    assert e["ft"][A1] == e["ft"][B1] and e["inside"][e["ft"][A1]] == 2 and not e["open"]
    assert (e["lo"][A1], e["hi"][A1], e["hi"][B1]) == (0, 2, 2)
    assert fw.statuses(e, w, [A1, B1]) == {A1: "tied_only", B1: "tied_only"}
    assert fw.min_copies(set(), [frozenset("AB")]) == (1, True)

def test_one_emitted_one_not():
    e, w, s, b = L({"r1": [("A", "a1", 0, 1), ("B", fw.OTHER, 0, 0)]})
    k = e["ft"][A1]
    assert e["open"][k] == 1 and e["inside"][k] == 0 and e["hi"][A1] == 1 and e["lo"][A1] == 0
    assert sum(e["placed"].values()) + sum(e["classes"].values()) == 1            # P1: nobody lost or doubled

def test_two_records_at_one_locus():
    e, w, s, b = L({"r1": [("A", "a1", 0, 1), ("A", "a2", 0, 0), ("B", "b1", 0, 0)], "r2": [("A", "a1", 0, 1), ("A", "a2", 0, 0)]})
    assert e["ft"][A1] == e["ft"][A2] == e["ft"][B1] and e["ft"][A1] in e["collision"]
    pl, mx = fw.locus_counts(s["exact"], b)
    assert (pl["A"], mx["A"], pl["B"], mx["B"]) == (1, 2, 0, 1) and e["lo"][A1] == 0   # r2: placed on the locus, on no node

def test_statuses_and_primary_below_best():
    recs = {"r1": [("A", "a1", 0, 1), ("B", "b1", 5, 0)],          # best at A, B within 10
            "r2": [("C", "c1", 3, 1), ("D", "d1", 0, 0)]}          # primary 3 below the best record: tied between both
    e, w, s, b = L(recs, base={("E", "e1"): 1})
    st = fw.statuses(e, w, [A1, B1, ("C", "c1"), ("D", "d1"), ("E", "e1"), ("F", "f1")])
    assert st == {A1: "best_only", B1: "near_only", ("C", "c1"): "tied_only", ("D", "d1"): "tied_only",
                  ("E", "e1"): "certified", ("F", "f1"): "no_band_record"}
    assert s["exact"]["r2"] == frozenset({("C", "c1"), ("D", "d1")})

def test_min_link_repeat_and_unseen():
    recs = {"r1": [("A", "a1", 0, 1), ("B", "b1", 0, 0)]}
    assert L(recs)[0]["inside"] and not L(recs, min_link=2)[0]["inside"] and L(recs, min_link=2)[0]["open"]
    assert not L(recs, repeat=frozenset({"B"}))[0]["inside"]
    e, w, s, b = L({"r1": [("A", "a1", 0, 1), (fw.UNSEEN, ".", 0, 0)]})
    assert e["lo"][A1] == 0 and e["open"][e["ft"][A1]] == 1 and fw.locus_counts(s["exact"], b)[0]["A"] == 1

def test_min_copies():
    assert fw.min_copies(set(), [frozenset("AB"), frozenset("BC")]) == (1, True)
    assert fw.min_copies({"A"}, [frozenset("BC")]) == (2, True)
    assert fw.min_copies(set(), [frozenset("AB"), frozenset("CD"), frozenset("BC")], cap=0) == (2, False)

def test_invariance_and_bounds_random():
    rng = random.Random(7); nodes = [(l, l.lower() + str(i)) for l in "ABCD" for i in (1, 2)]
    recs, origin = {}, {}
    for i in range(300):
        k = rng.sample(nodes, rng.choice((1, 1, 2, 3))); origin[i] = k[0]
        recs[i] = [(l, m, 0, int(j == 0)) for j, (l, m) in enumerate(k)]
    e, w, s, b = L(recs)
    flipped = {i: [(l, m, g, int(j == len(r) - 1)) for j, (l, m, g, p) in enumerate(r)] for i, r in recs.items()}
    assert fw.bands(flipped) == s                                                  # P2: the coin does not matter
    x = collections.Counter(origin.values())
    assert all(e["lo"][n] <= x[n] <= e["hi"][n] for n in nodes)                    # P3 per node
    for key in set(e["ft"].values()):                                              # P5: exact total without open molecules
        comp = [n for n in nodes if e["ft"].get(n) == key]
        if not e["open"][key]:
            assert sum(x[n] for n in comp) == sum(e["placed"][n] for n in comp) + e["inside"][key]
```

- [ ] **Step 2: Run to see them fail**

Run: the pytest command. Expected: 7 new tests FAIL with `AttributeError: module 'famwise' has no attribute 'bands'`.

- [ ] **Step 3: Implement `bands`, `layer`, `locus_counts`, `statuses`, `min_copies`**

Union-find over linkable members, fed by a Counter of member pairs summed over classes. `min_copies`: drop sets already
met by `required`, drop supersets, then `itertools.combinations` over the loci of the remaining sets with growing size,
counting tries against `cap`.

- [ ] **Step 4: Run the tests**

Run: the pytest command. Expected: 9 passed.

---

### Task 3: the switch in `tools/eq_classes.py`

**Files:**
- Modify: `tools/eq_classes.py` (`scan` :71-91, `main` :156-204, the docstring's usage line)
- Test: `tools/test_famwise.py`

**Interfaces:**
- Consumes: `famwise.load_nodes`, `record_chain`, `member_of`, `FW_WINDOW`, `UNSEEN` (imported inside `main`, only with
  the switch, so `psi_bunches.py` and `test_eq_classes.py` keep importing `load_ties`, `floor_fn`, `summarise`, `em`).
- Produces: flag `--famwise FW_PREFIX` writing
  - `FW_PREFIX.reads.tsv`, header `read	locus	member	gap	ms	AS	primary`: one row per kept record of every molecule
    the score table calls tied within `max(window, FW_WINDOW)` (or whose primary is not its best). `gap` = best `ms` of
    the score table minus the record's `ms` (integer); `primary` 1 / 0. A tied molecule with exactly one kept record
    whose `second_ms` (score-table field 7) is within the same window gets one extra row
    `read	unseen	.	<best - second>	<second>	NA	0`. Rows sorted by read, gap, locus, member.
  - `FW_PREFIX.placed.tsv`, header `locus	member	n`: untied molecules by the member of their primary record
    (`outside	.` included), sorted.

- [ ] **Step 1: Reference outputs BEFORE editing (40 s at 2 processes on the login node)**

```bash
W=/anvme/workspace/v138da13-rustle; P=$W/famwise_design_2026-11-18/panel; mkdir -p $P
B=$W/reads/GGO/GGO_KB3781_fib_mm.bam; G=$W/famwise_look_2026-11-16/eq/k.families.named.gtf
$PY tools/eq_classes.py $G $B $B.as_table.tsv $P/base --threads 2 --contigs NC_073247.2,NC_073248.2 > $P/base.log
```
Expected in `base.log`: two lines, `ms exact: ... 35 classes ...` and `ms wide ...`.

- [ ] **Step 2: Write the failing test** (tiny sorted, indexed BAM with pysam, pattern of `bench/stgaps/test_famasm.py:17-48`
  plus `a.set_tag("ms", n)`; the 17-column score table written as text, first line `#as_table`)

```python
import subprocess
EQ = os.path.join(os.path.dirname(os.path.abspath(__file__)), "eq_classes.py")
FWPY = os.path.join(os.path.dirname(os.path.abspath(__file__)), "famwise.py")

def tiny(tmp_path):
    """c1: locus A (exons 1-100, 151-250, +) and locus B (1001-1100, 1151-1250, +), one 2-exon transcript each
    (A.1, B.1). Every record carries AS equal to its ms.
    t1, t2: primary at A, secondary at B, ms 200 / 200. u1: one record at A. w1: primary at A ms 200, secondary at B
    ms 185. s1: unspliced primary at A. Returns (gtf, bam, table) paths."""
    ...  # build with pysam.AlignedSegment; cigars 100=50N100= ; sort + index; table rows best_ms / second_ms / primary_ms

def test_switch_keeps_eq_tables_and_writes_reads(tmp_path):
    gtf, bam, tab = tiny(tmp_path)
    run = lambda out, *x: subprocess.run([sys.executable, EQ, gtf, bam, tab, str(tmp_path / out), "--threads", "1", *x], capture_output=True, text=True)
    a, b = run("off"), run("on", "--famwise", str(tmp_path / "fw"))
    assert a.returncode == b.returncode == 0 and a.stdout == b.stdout and len(a.stdout.splitlines()) == 2
    for t in ("loci", "classes"):
        assert (tmp_path / f"off.{t}.tsv").read_bytes() == (tmp_path / f"on.{t}.tsv").read_bytes()
    rows = [l.split("\t") for l in (tmp_path / "fw.reads.tsv").read_text().splitlines()]
    assert rows[0] == "read locus member gap ms AS primary".split()
    assert ["w1", "B", "B.1", "15", "185", "185", "0"] in rows and sum(r[0] == "t1" for r in rows) == 2
    assert (tmp_path / "fw.placed.tsv").read_text().splitlines()[1:] == ["A\tA.1\t1", "A\tU\t1"]
    assert run("bad", "--famwise", str(tmp_path / "x"), "--ratio", "0.98").returncode == 2

def test_unseen_rival(tmp_path):
    gtf, bam, tab = tiny(tmp_path)   # scan only a contig that holds t1's primary: build the BAM with locus B on contig c2 for this test
    ...  # run with --contigs c1 --famwise; assert a row ["t1", "unseen", ".", "0", "200", "NA", "0"]
```

- [ ] **Step 3: Run to see it fail**

Expected: FAIL, `error: unrecognized arguments: --famwise`.

- [ ] **Step 4: Implement the switch**

In `main`: `--famwise`; with `a.ratio is not None` exit 2 with
`eq_classes: --famwise needs the window band (no --ratio)`. With the switch the load floor and the record floor are
`floor_fn(max(a.window, FW_WINDOW))`; a new global holds the nodes (`None` without the switch). `scan` keeps returning
`(placed, pairs)` records of today's shape when the global is `None`. With it, each pair also carries
`wide` (`not secondary or ms >= best - a.window`), member, gap, `ms`, AS and primary, and `scan` also returns a Counter
of `(locus, member)` for untied molecules. In the merge, a record enters `sets["wide"]` only when `wide`, `sets["exact"]`
only when exact, exactly as today. Write the two files after the existing tables; print nothing.

- [ ] **Step 5: Run the tests, then the panel**

Run: the pytest command. Expected: 11 passed.
```bash
$PY tools/eq_classes.py $G $B $B.as_table.tsv $P/off --threads 2 --contigs NC_073247.2,NC_073248.2 > $P/off.log
$PY tools/eq_classes.py $G $B $B.as_table.tsv $P/on --threads 2 --contigs NC_073247.2,NC_073248.2 --famwise $P/fw > $P/on.log
for t in loci classes; do cmp $P/base.$t.tsv $P/off.$t.tsv && cmp $P/base.$t.tsv $P/on.$t.tsv; done; cmp $P/base.log $P/on.log
```
Expected: no output from `cmp` (all identical).

---

### Task 4: tables, GTF attributes and CLI of `tools/famwise.py`

**Files:**
- Modify: `tools/famwise.py`
- Test: `tools/test_famwise.py`

**Interfaces:**
- Consumes: Tasks 1-3.
- Produces: CLI `famwise.py GTF FW_PREFIX [--eq-loci EQ.loci.tsv] [--window 10] [--min-link 1]`, reading
  `FW_PREFIX.reads.tsv` and `FW_PREFIX.placed.tsv`, writing:
  - `FW_PREFIX.tx.tsv`: `band	node	locus	chrom	start	end	strand	transcripts	ft	lo	hi	status` (one row per node of the
    GTF and band; `transcripts` comma-joined; coordinates of the node's first transcript).
  - `FW_PREFIX.classes.tsv`: `band	reads	size	members	ft	kind` (`members` = sorted `locus|member`, comma-joined;
    `kind` inside / open; `ft` = `.` for a class touching no family transcript).
  - `FW_PREFIX.families.tsv`: `band	ft	nodes	loci	copies	molecules	placed	inside	open	certified	best_only	tied_only	near_only	no_band_record	min_copies	min_copies_exact	collision	rep`
    (`ft` = `FT1`, `FT2`, ... by falling `molecules` = placed + inside, then smallest node; `copies` = sorted loci,
    comma-joined; `rep` = the node with most placed molecules, then smallest locus start, then node id; `min_copies` from
    `min_copies(loci holding a certified or best_only node, loci sets of the band's inside classes)`).
  - `FW_PREFIX.gtf`: the input GTF, every line kept in order; multi-exon transcript lines gain
    ` fw_tx "FT7"; fw_status "certified"; fw_lo "3"; fw_hi "10"; fw_copies "2"; fw_rep "yes";` from the wide band.
  - With `--eq-loci`: `repeat` loci are the rows whose `label` or `label_w` is `repeat`; and a self-check: for every
    row, `placed` / `max` (exact) and `placed_w` / `max_w` (wide) equal `locus_counts`. First mismatch ->
    `famwise: P1 fails at locus L (band): eq placed/max a/b, famwise c/d`, exit 1.
  - Stdout: one line per band: family transcripts of two or more nodes, transcripts in them, molecules inside, open,
    collisions; then the count of nodes per status.

- [ ] **Step 1: Write the failing tests**

```python
def test_end_to_end_tiny(tmp_path):
    gtf, bam, tab = tiny(tmp_path)
    subprocess.run([sys.executable, EQ, gtf, bam, tab, str(tmp_path / "eq"), "--threads", "1", "--famwise", str(tmp_path / "fw")], check=True)
    r = subprocess.run([sys.executable, FWPY, gtf, str(tmp_path / "fw"), "--eq-loci", str(tmp_path / "eq.loci.tsv")], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    fam = list(csv.DictReader(open(tmp_path / "fw.families.tsv"), delimiter="\t"))      # import csv at the top
    top = next(f for f in fam if f["band"] == "exact" and f["ft"] == "FT1")
    assert (top["nodes"], top["copies"], top["inside"], top["open"], top["min_copies"]) == ("2", "A,B", "2", "0", "1")
    tx = (tmp_path / "fw.tx.tsv").read_text()
    assert "\tA.1\t" in tx and "certified" in tx          # u1 is placed on A.1 alone
def test_gtf_attributes(tmp_path):
    ...  # GTF of Task 1 (A.1 and A.2 share a chain): both transcript lines carry the same fw_tx and fw_status; line count unchanged; S.1 untouched
def test_p1_mismatch_exits_1(tmp_path):
    ...  # edit one `placed` value of eq.loci.tsv -> returncode 1, "P1 fails at locus" in stderr
def test_empty(tmp_path):
    ...  # reads.tsv and placed.tsv with headers only, GTF with one single-exon transcript -> rc 0, four tables with headers only
```

- [ ] **Step 2: Run to see them fail.** Expected: `famwise.py` exits with usage error / no `main`.

- [ ] **Step 3: Implement `main` and the writers** (argparse, `if __name__ == "__main__": main()`; coordinates from the
  GTF's transcript exons; rows sorted for determinism).

- [ ] **Step 4: Run the tests, then the panel**

Run: the pytest command. Expected: 15 passed.
```bash
$PY tools/famwise.py $G $P/fw --eq-loci $P/on.loci.tsv | tee $P/fw.log
```
Expected: exit 0 (P1 holds on real reads). The exact line should show the CD99 X / Y family transcripts (preview of the
design: 9 family transcripts, about 13,900 molecules inside) and no row with locus `unseen` in `$P/fw.reads.tsv` beyond
molecules tied to a third contig. The printed numbers stay in `$P/fw.log`; they are not asserted.

---

### Task 5: driver stage, usage text, README row, dispatch test

**Files:**
- Modify: `tools/rustle_pipeline.sh` (`stage_eqclasses` :1118-1136; new `stage_famwise` after it; usage header after
  the `eqclasses` entry at :80-86; dispatch `case` line :1860)
- Modify: `README.md` (table under `### Stages added in November 2026`, after the `eqclasses` row at :89)
- Test: `tests/pipeline_driver_dispatch.rs`

**Interfaces:**
- Consumes: `eq_classes.py --famwise`, `famwise.py` CLI.
- Produces: `tools/rustle_pipeline.sh famwise --bam B --fasta G --out PREFIX [--names ANNOT] [--threads N]`, writing
  `PREFIX.eq.loci.tsv`, `PREFIX.eq.classes.tsv` (as `eqclasses`) and `PREFIX.famwise.{reads,placed,tx,classes,families}.tsv`,
  `PREFIX.famwise.gtf`, log `PREFIX.famwise.log`. Env: `RUSTLE_EQ_WINDOW` (10), `RUSTLE_FW_MIN_LINK` (1).

- [ ] **Step 1: Write the failing dispatch test** in `tests/pipeline_driver_dispatch.rs`, copying
  `altends_needs_the_assembly_and_then_writes_the_altends_gtf` (:470): name
  `famwise_runs_eqclasses_with_the_switch_then_the_layer`; script = variables (`OUT`, `BAM`, `FASTA`, `FAM_GTF`,
  `JOB_THREADS=2`, `NAMES=`, `RP_DRIVER_DIR=tools`, `RUSTLE_EXCLUDE_BED=`) + `function_block("fam_gtf_guard")` +
  `function_block("stage_eqclasses")` + `function_block("stage_famwise")` + the call. Stub `RUSTLE_PYTHON='bash py.sh'`
  where `py.sh` appends its arguments to `calls.txt` and creates the files each tool would write (for `eq_classes.py`:
  `$4.loci.tsv`, `$4.classes.tsv` and, when `--famwise X` is among the arguments, `X.reads.tsv`, `X.placed.tsv`; for
  `famwise.py`: `$2.{tx,classes,families}.tsv`, `$2.gtf`). Assert: exit 0; `calls.txt` has two lines, the first
  containing `eq_classes.py` and `--famwise`, the second `famwise.py` and `--eq-loci`; the six `run.famwise.*` files and
  the two `run.eq.*` tables exist; with `RUSTLE_EQ_RATIO=0.98` exported the exit code is 2 and stderr contains
  `famwise needs the window band`. Also assert `all` is unchanged: the existing test
  `all_runs_assemble_families_assign_flag_and_neither_catalog_nor_candidates` (:120) must still pass untouched.

- [ ] **Step 2: Run to see it fail**

Run: `CARGO_TARGET_DIR=$HOME/rustle_target cargo test --release --test pipeline_driver_dispatch > test_disp.log 2>&1; tail -5 test_disp.log`
Expected: 1 failed (`stage_famwise` not found), the others pass.

- [ ] **Step 3: Implement**

`stage_eqclasses`: add `local fw=(); [ -z "${RP_FAMWISE:-}" ] || fw=(--famwise "$OUT.eq.tmp.fw")` and pass
`"${fw[@]}"` to `eq_classes.py` (its failure path `rm -f "$OUT.eq.tmp".*` then covers the new files).
`stage_famwise`: refuse `RUSTLE_EQ_RATIO` (exit 2, message above); `RP_FAMWISE=1 stage_eqclasses`; the GTF is the one
`eqclasses` used (`$OUT.eq.named.gtf` when `NAMES` is set, else `$FAM_GTF`); run `famwise.py "$gtf" "$OUT.eq.tmp.fw"
--eq-loci "$OUT.eq.loci.tsv" --window "${RUSTLE_EQ_WINDOW:-10}" --min-link "${RUSTLE_FW_MIN_LINK:-1}" > "$OUT.famwise.log" 2>&1`;
on failure say `famwise failed (see $OUT.famwise.log); no other product is affected`, remove `"$OUT.eq.tmp.fw".*`, exit 1;
then `mv -f` each `"$OUT.eq.tmp.fw".X` to `"$OUT.famwise.X"` and `say` the log's last lines.
Dispatch: add `famwise) stage_famwise;;` on the `eqclasses` line. `all` is not touched.
Usage header (after the `eqclasses` entry):
```
#   tools/rustle_pipeline.sh famwise --bam B --fasta G --out PREFIX [--names ANNOT] [--threads N]   OPT-IN, named only:
#     `eqclasses` plus PREFIX.famwise.*: transcripts of different loci that the same tied molecule carries are joined into
#     family transcripts (families.tsv: copy set, each molecule once); per transcript the molecules placed only there and
#     the maximum with every tied one (tx.tsv), and a status: certified, best_only, tied_only, near_only, no_band_record.
#     Bands on minimap2's ms (exact, and within RUSTLE_EQ_WINDOW); every record row also carries AS (reads.tsv, within 20).
#     Reporting only: no read is assigned, no transcript added or dropped, no other product changes (tools/famwise.py).
```
README row:
`| `famwise` | `PREFIX.famwise.families.tsv`, `.tx.tsv`, `.classes.tsv`, `.reads.tsv`, `.gtf` (transcripts joined across loci by tied molecules, each molecule counted once; per transcript a count range and a status; reporting only) | none |`

- [ ] **Step 4: Run the dispatch tests and the Python tests**

Run: the cargo command, then the pytest command. Expected: all dispatch tests pass (one more than before); 15 passed.

---

### Task 6: whole-sample runs and the record of what came out

**Files:**
- Create (workspace, not the repo): `/anvme/workspace/v138da13-rustle/famwise_design_2026-11-18/run/famwise.sbatch`,
  outputs `run/kb.*`, `run/a1.*`

**Interfaces:**
- Consumes: the stage of Task 5.
- Produces: `run/kb.famwise.*` (KB3781) and `run/a1.famwise.*` (A119b), which part 2 reads.

- [ ] **Step 1: Job script** (`-p singlenode -N1 -n1 -c 72 -t 02:00:00`), one sample per job, arguments
  `NAME BAM FASTA DEFAULT_RUN_PREFIX`: link `DEFAULT_RUN_PREFIX.gtf`, `.families.gtf` and `.molecules.tsv` under
  `run/NAME.*` (so nothing is written beside the default run), then
  `bash tools/rustle_pipeline.sh famwise --bam BAM --fasta FASTA --out run/NAME --threads 64 --bin $HOME/rustle_target/release`
  under `/usr/bin/time -v`, with `TMPDIR` left at the node-local default for any temporary file. Default runs: `artifact_refresh_2026-11-16/apes/GGO_KB3781_fib_mm/def/x` and
  `artifact_refresh_2026-11-16/h.def/a1`; BAMs `reads/GGO/GGO_KB3781_fib_mm.bam`, `reads/HSA/A119b.t2t.bam`; FASTA
  paths from the `@PG` line of each BAM.
- [ ] **Step 2: Submit KB3781** (`sbatch`). Expected: exit 0 in about 10 min; `kb.famwise.log` ends with the two band
  lines; P1 self-check passed (no `P1 fails`). If `sbatch` is denied by the storage guard, stop and report.
- [ ] **Step 3: Check against the design's preview** (descriptive, not a bar): exact band, family transcripts of two or
  more nodes (preview at the old tie-group loci: 134 joining 303 transcripts, CD99 13,868 molecules); `cmp` of
  `kb.eq.loci.tsv` against a run of the same command with stage `eqclasses`.
- [ ] **Step 4: Submit A119b** (about 1 h). Check the TSPY region `chrY:8940000-9822000`: family transcripts there
  (preview within 10: 13 plus 2 single, no collision).
- [ ] **Step 5: Append the numbers** of steps 3 and 4, with the commands, to the design's "Preview" section under a line
  `Stage output (date)`, replacing nothing. No register row: this is a product, not a result.

## Not in this plan
The measurement (E1-E3), its pre-registration and the advisor page: part 2. Rust changes, a new default, estimates
inside a range: out of scope by the spec.
