# O2 from the sample's own variant calls — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Place AS-tied multimapping reads using per-copy variant calls made by longcallR on the sample's confidently placed reads. Abstain when nothing separates the copies, and flag reads no reference copy explains. Show it end to end on GTF2IRD2 (HG002) under a pre-registration.

**Architecture:** Three Python files under `bench/o3_graph/`, written in the style of `o2correct.py`:
- `o2vc_rule.py`: the pure decision rule (genotype sets, assign or abstain, foreign alleles, flags).
- `o2vc_data.py`: the BAM, VCF and longcallR side (tied domain, teachers, read-frame columns, depth, caller, VCF reader).
- `o2vc.py`: the CLI (`run`, `score`, `score-worlds`, `copies-from-sets`, `copies-from-genes`).

`run` never reads truth; `score` adds the truth, the arms and the bars.

**Tech Stack:** Python 3 (`/home/linuxbrew/.linuxbrew/bin/python3`, pysam 0.23.3, pytest 9.0.3), longcallR 2.0.0 (bioconda), samtools 1.22.

**Spec:** `docs/superpowers/specs/2026-10-22-o2-variant-calls-design.md` (approved 2026-10-22).

## Global Constraints

**The rule (spec values):**
- **Tied domain:** a read whose alignments within 10 AS of its best (`SLACK = 9`) lie on ≥ 2 copies of one family and on no locus off it. minimap2's primary MAY lie outside that set (reported as `outside`).
- **Teacher:** a read whose primary is on the copy, with every alignment on another copy of its family ≥ 10 AS lower. A tied read never teaches.
- **Caller:** longcallR 2.0.0 at `/home/juanfra/miniforge3/envs/longcallr/bin/longcallR`, flags `-p hifi-masseq --min-mapq 0 --no-bam-output`, every other option default.
- **Known / unknown:** a copy is known at a position when teacher depth is ≥ 6 (`MIN_DEPTH`, the HiFi presets' min-depth). Depth is longcallR 2.0.0's own DP for the hifi-masseq preset (Task 3 ruling, 2026-10-22: primary teachers, read length ≥ 500, divergence < 0.5, any base quality, no HiFi read-end trim, the caller's poly-run rule; `teacher_depth(..., fetch)`). `RnaEdit` and `Dense` records → unknown.
- **Genotype:** the record's GT alleles, whatever its FILTER. No record → {reference base}.
- **Votes:** a read base votes only at base quality ≥ 10 (`MIN_BASEQ`).
- **Read frame:** every comparison is made in the read's own frame, as sequenced (FLNC = transcript sense). A copy aligned opposite to the read has its alleles complemented.

**Discipline:**
- Truth is read only by `score`, never by `run`.
- **Exclusions:**
  - PAR families (`/mnt/linuxdisk/tmp/o3_graph/famident/par_families.json`) are out of every pooled bar and reported apart.
  - Families with a chrM copy (S19MTCO2, S19MTND6, S19MTRNR2L8) are out of the HG002 runs (`--exclude-chrM`).
- **Development** uses only S19SMN2 and S19MSTRG3150 on HG002.
- **Held back until the prereg (Task 8) is committed:** GTF2IRD2, NBPF1, PDPR, every other HG002 family, DAZ, gorilla NPIP, the round-23 worlds and HG00733.
- Never pool human and gorilla numbers.
- **Determinism:** iterate sorted, never over a bare `set` (metric trap: PYTHONHASHSEED).

**Practical:**
- **Location:** worktree `/mnt/linuxdisk/home/juanfraitu/rustle_soto_subset`, branch `o3-fast`; outputs under `/mnt/linuxdisk/tmp/o3_graph/o2vc/`.
- **WSL2:** foreground, serial, no `pkill -f`. The HG002 run and the gorilla extraction run under `tools/rlock.sh heavy` (`RLOCK_TIMEOUT=5400`); small runs under `tools/rlock.sh light`.
- **Tests:** `python3 -m pytest bench/o3_graph/test_o2vc_rule.py bench/o3_graph/test_o2vc_data.py bench/o3_graph/test_o2vc.py -q`, run from the worktree root.
- **Commits:** every commit message ends with:
  ```
  Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01JeH8K18apUtXvWEi2vp84Y
  ```

## Review Focus

1. **A copy aligned opposite to the read** (inverted paralogs, e.g. NPIP, TBC1D3). Alleles must be complemented into the read frame before any comparison, or reads get placed on strand-mismatched evidence. Pinned in Task 1 (`to_read_frame`) and Task 3 (`columns` with a flipped copy).
2. **Secondary records without SEQ, carrying hard clips.** Read-frame positions must count both clip kinds, or every column shifts by the clip length. Pinned in Task 3 (`read_frame_positions`, the `scan` secondary).
3. **VCF coordinates.** POS is 1-based; positions are used 0-based. A REF that disagrees with the genome must stop the run. Pinned in Task 2 (`parse_vcf`) and Task 4 (REF mismatch raises).
4. **A copy with no teachers, or too few.** It must stay unknown: never beaten, never beating, no longcallR call. Pinned in Task 1 (unknown copy) and Task 4 (the caller is not run for a teacherless copy).
5. **PAR and chrM families leaking into pooled numbers.** Pinned in Task 4 (`--exclude-chrM`) and Task 5 (`--exclude` reports PAR apart).

---

### Task 1: The decision rule

**Files:**
- Create: `bench/o3_graph/o2vc_rule.py`
- Test: `bench/o3_graph/test_o2vc_rule.py`

**Interfaces:**
- Produces:
  - `MIN_DEPTH = 6`, `MIN_BASEQ = 10`, `UNKNOWN_FILTERS = frozenset({"RnaEdit", "Dense"})`
  - `Call = namedtuple("Call", "ref alts gt filter dp")`: `alts` a tuple of str; `gt` a tuple of allele indices (0 = ref).
  - `genotype(call: Call | None, depth: int, ref_base: str) -> frozenset | None`
  - `to_read_frame(g: frozenset | None, flip: bool) -> frozenset | None`
  - A column is a `dict(base: str, qual: int, g: {copy: frozenset | None})`, with genotypes already in the read frame.
  - `assign(copies: Iterable[str], cols: list[dict]) -> tuple[str | None, frozenset, int]`: the placed copy or None, the kept set, and the number of separating columns voting for the winner (0 when none wins).
  - `foreign(col: dict) -> str`: `"foreign"`, `"editing"` or `""`.
  - `flags(by_read: dict[str, list[tuple[int, frozenset]]]) -> list[dict]`. The input maps read → [(column id within the read, keys)], each key a tuple `(copy, chrom, pos, genome-strand base)`. Each output is `dict(id=int, alleles=sorted list of keys, reads=sorted list)`.

- [ ] **Step 1: Write the failing tests**

```python
import importlib.util
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("o3_graph_o2vc_rule", os.path.join(HERE, "o2vc_rule.py"))
R = importlib.util.module_from_spec(spec)
spec.loader.exec_module(R)


def col(base, qual=30, **g):
    """The read's base and quality, each copy's genotype in the read's frame (F_0 -> F.0; None = unknown)."""
    return dict(base=base, qual=qual, g={k.replace("_", "."): None if v is None else frozenset(v) for k, v in g.items()})


def test_a_copy_s_genotype_is_its_record_s_gt_alleles_whatever_the_filter():
    C = R.Call
    assert R.genotype(None, 6, "A") == frozenset("A")
    assert R.genotype(C("A", ("G",), (0, 1), "LowQual", 9), 9, "A") == frozenset("AG")
    assert R.genotype(C("A", ("G",), (1, 0), "PASS", 9), 9, "A") == frozenset("AG")
    assert R.genotype(C("A", ("G",), (1, 1), "PASS", 9), 9, "A") == frozenset("G")
    assert R.genotype(C("A", ("G", "T"), (1, 2), "Multiallelic", 9), 9, "A") == frozenset("GT")
    assert R.genotype(C("A", ("G",), (0, 0), "HomRef", 9), 9, "A") == frozenset("A")


def test_a_copy_is_unknown_below_depth_six_or_at_an_editing_or_dense_record():
    C = R.Call
    assert R.genotype(None, 5, "A") is None
    assert R.genotype(C("A", ("G",), (0, 1), "RnaEdit", 30), 30, "A") is None
    assert R.genotype(C("A", ("G",), (0, 1), "Dense", 30), 30, "A") is None


def test_alleles_are_complemented_into_the_read_s_frame_for_a_copy_aligned_opposite_to_the_read():
    assert R.to_read_frame(frozenset("AG"), True) == frozenset("TC")
    assert R.to_read_frame(frozenset("AG"), False) == frozenset("AG")
    assert R.to_read_frame(None, True) is None


def test_a_copy_is_placed_when_a_separating_column_votes_for_it_and_none_for_any_rival():
    assert R.assign(["F.0", "F.1"], [col("C", F_0="C", F_1="T"), col("A", F_0="A", F_1="A")]) == ("F.0", frozenset({"F.0"}), 1)
    three = [col("C", F_0="C", F_1="T", F_2="C"), col("G", F_0="G", F_1="G", F_2="A")]
    assert R.assign(["F.0", "F.1", "F.2"], three) == ("F.0", frozenset({"F.0"}), 2)


def test_contradicting_columns_and_shared_alleles_leave_the_copy_set():
    assert R.assign(["F.0", "F.1"], [col("C", F_0="C", F_1="T"), col("G", F_0="A", F_1="G")]) == (None, frozenset({"F.0", "F.1"}), 0)
    assert R.assign(["F.0", "F.1"], [col("C", F_0="CT", F_1="T")]) == (None, frozenset({"F.0", "F.1"}), 0)


def test_an_unknown_copy_is_never_beaten_and_a_low_quality_base_never_votes():
    assert R.assign(["F.0", "F.1", "F.2"], [col("C", F_0="C", F_1="T", F_2=None)]) == (None, frozenset({"F.0", "F.2"}), 0)
    assert R.assign(["F.0", "F.1"], [col("C", qual=9, F_0="C", F_1="T")]) == (None, frozenset({"F.0", "F.1"}), 0)


def test_a_foreign_allele_needs_every_copy_known_and_none_holding_the_read_s_base():
    assert R.foreign(col("G", F_0="C", F_1="T")) == "foreign"
    assert R.foreign(col("G", F_0="C", F_1=None)) == ""
    assert R.foreign(col("C", F_0="C", F_1="T")) == ""
    assert R.foreign(col("G", qual=5, F_0="C", F_1="T")) == ""
    assert R.foreign(col("N", F_0="C", F_1="T")) == ""
    assert R.foreign(col("G", F_0="A", F_1="AT")) == "editing"          # A>G in the read's (transcript) frame


def _k(pos, b, copy="F.0", chrom="chr1"):
    return (copy, chrom, pos, b)


def test_a_flag_is_two_reads_sharing_two_linked_foreign_alleles_at_different_columns():
    by_read = {"r1": [(10, frozenset({_k(100, "G"), _k(500, "G", "F.1", "chr2")})), (20, frozenset({_k(110, "T")}))],
               "r2": [(11, frozenset({_k(100, "G")})), (21, frozenset({_k(110, "T")}))],
               "r3": [(12, frozenset({_k(100, "G")}))]}
    f = R.flags(by_read)
    assert [x["reads"] for x in f] == [["r1", "r2"]]
    assert f[0]["alleles"] == [_k(100, "G"), _k(110, "T")] and f[0]["id"] == 0


def test_one_shared_foreign_column_seen_on_two_copies_makes_no_flag():
    same = frozenset({_k(100, "G"), _k(500, "G", "F.1", "chr2")})
    assert R.flags({"r1": [(10, same)], "r2": [(11, same)]}) == []


def test_flags_do_not_depend_on_the_order_reads_arrive():
    by_read = {f"r{i}": [(1, frozenset({_k(100 + i % 2, "G")})), (2, frozenset({_k(200, "T")}))] for i in range(6)}
    keys = list(by_read)
    random.Random(7).shuffle(keys)
    assert R.flags({k: by_read[k] for k in keys}) == R.flags(by_read)
```

- [ ] **Step 2: Run the tests and see them fail**

Run: `python3 -m pytest bench/o3_graph/test_o2vc_rule.py -q`
Expected: FAIL (`o2vc_rule.py` does not exist).

- [ ] **Step 3: Implement `o2vc_rule.py`**

**Definitions:**
- A column **separates** c and d when both genotypes are known and disjoint.
- **c beats d** when ≥ 1 separating column with `qual >= MIN_BASEQ` has the read's base in c's genotype, and none has it in d's.
- **`assign`** places c when c beats every other copy. Otherwise it returns the copies no copy beats.
- **`foreign`:** the read's base must be one of `ACGT` with `qual >= MIN_BASEQ`. Every copy must be known, and no copy may hold the base. The column is `"editing"` instead when the base is `G` and every copy's genotype holds `A`.

**`flags` algorithm** (the signature does not determine it):
1. A key is **recurrent** when ≥ 2 distinct reads carry it.
2. For each read, take every pair of recurrent keys from two **different** column ids. A pair is **linked** when ≥ 2 distinct reads give it.
3. Union-find over linked pairs gives the components, each with ≥ 2 keys.
4. A component's reads are the reads carrying its keys at ≥ 2 distinct column ids.
5. Drop components with < 2 reads.
6. Sort components by their smallest key and number them from 0. Iterate sorted throughout.

- [ ] **Step 4: Run the tests and see them pass**

Run: `python3 -m pytest bench/o3_graph/test_o2vc_rule.py -q`
Expected: 10 passed.

- [ ] **Step 5: Commit**

```bash
git add bench/o3_graph/o2vc_rule.py bench/o3_graph/test_o2vc_rule.py
git commit -m "o2vc rule: per-copy genotype sets from variant calls, assign-or-abstain, foreign alleles and linked-allele flags"
```

### Task 2: longcallR install, runner and VCF reader

**Files:**
- Create: `bench/o3_graph/o2vc_data.py` (first part)
- Test: `bench/o3_graph/test_o2vc_data.py`

**Interfaces:**
- Consumes: `o2vc_rule.Call` (load the sibling with an importlib helper, as `o2correct._sibling` does).
- Produces:
  - `LONGCALLR = "/home/juanfra/miniforge3/envs/longcallr/bin/longcallR"`
  - `run_longcallr(teacher_bam: str, genome: str, out_prefix: str, binary: str = LONGCALLR, threads: int = 2, region: str | None = None) -> str` (the VCF path)
  - `parse_vcf(path: str) -> dict[tuple[str, int], Call]`: keys `(chrom, POS - 1)`; handles `.vcf` and `.vcf.gz`; FORMAT keys read by name.
  - `longcallr_version(binary: str = LONGCALLR) -> str`

- [ ] **Step 1: Install longcallR and see what it writes**

1. Install:
   ```bash
   /home/juanfra/miniforge3/bin/mamba create -y -n longcallr -c conda-forge -c bioconda longcallr=2.0.0
   ```
   Then run `/home/juanfra/miniforge3/envs/longcallr/bin/longcallR --version`. Expected: a 2.0.0 version string.
2. Smoke run on the development family SMN2. This is primaries only and not teacher-filtered, which is fine for a smoke test:
   ```bash
   mkdir -p /mnt/linuxdisk/tmp/o3_graph/o2vc/smoke && cd /mnt/linuxdisk/tmp/o3_graph/o2vc/smoke
   samtools view -b -F 2308 ../../hg002proof/chm13.bam chr5:70808816-70837687 > smn2.bam && samtools index smn2.bam
   /home/juanfra/miniforge3/envs/longcallr/bin/longcallR -b smn2.bam -f /mnt/linuxdisk/home/juanfraitu/npip_ladder/idx/target.fa -p hifi-masseq --min-mapq 0 --no-bam-output -t 2 -o smn2 -r chr5:70808816-70837687
   ls
   ```
   Expected: one VCF. Note its exact file name, whether it is gzipped, and the FILTER values present. These fix `run_longcallr`'s return value and the next test.

- [ ] **Step 2: Write the failing tests**

```python
import gzip
import importlib.util
import os
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    s = importlib.util.spec_from_file_location(f"o3_graph_{name}_t", os.path.join(HERE, f"{name}.py"))
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


D, R = _load("o2vc_data"), _load("o2vc_rule")
VCF = ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
       "chr1\t101\t.\tA\tG\t30\tPASS\t.\tGT:GQ:PS:DP:AF:PQ\t0|1:30:101:12:0.50:20\n"
       "chr1\t201\t.\tC\tT,G\t20\tMultiallelic\t.\tGT:GQ:DP:AF\t1/2:20:10:0.40,0.60\n"
       "chr1\t301\t.\tT\tC\t5\tLowQual\t.\tGT:GQ:DP:AF\t0/1:5:8:0.30\n")


def test_parse_vcf_reads_every_record_0_based_with_gt_filter_and_depth(tmp_path):
    for p in (tmp_path / "x.vcf", tmp_path / "x.vcf.gz"):
        (gzip.open if p.suffix == ".gz" else open)(p, "wt").write(VCF)
        v = D.parse_vcf(str(p))
        assert v[("chr1", 100)] == R.Call("A", ("G",), (0, 1), "PASS", 12)
        assert v[("chr1", 200)] == R.Call("C", ("T", "G"), (1, 2), "Multiallelic", 10)
        assert v[("chr1", 300)] == R.Call("T", ("C",), (0, 1), "LowQual", 8)


def test_longcallr_runs_with_the_masseq_preset_min_mapq_0_and_no_bam(monkeypatch, tmp_path):
    seen = {}

    def fake(cmd, **kw):
        seen["cmd"] = cmd
        open(str(tmp_path / "o") + ".vcf", "w").write("#CHROM\n")       # adjust to the name Step 1 observed
        return subprocess.CompletedProcess(cmd, 0, "", "")
    monkeypatch.setattr(D.subprocess, "run", fake)
    out = D.run_longcallr("t.bam", "g.fa", str(tmp_path / "o"), binary="LC", threads=2, region="chr1:1-100")
    assert seen["cmd"] == ["LC", "-b", "t.bam", "-f", "g.fa", "-p", "hifi-masseq", "--min-mapq", "0", "--no-bam-output",
                           "-t", "2", "-o", str(tmp_path / "o"), "-r", "chr1:1-100"]
    assert out == str(tmp_path / "o") + ".vcf"
```

- [ ] **Step 3: Run the tests and see them fail**

Run: `python3 -m pytest bench/o3_graph/test_o2vc_data.py -q`
Expected: FAIL (`o2vc_data.py` does not exist).

- [ ] **Step 4: Implement `parse_vcf`, `run_longcallr` and `longcallr_version` in `o2vc_data.py`**

- `run_longcallr` calls `subprocess.run(cmd, capture_output=True, text=True)` with exactly the list above, adding `-r` only when `region` is given.
- On a non-zero exit it raises `RuntimeError` with the last 20 lines of stderr.

- [ ] **Step 5: Run the tests and see them pass**

Run: `python3 -m pytest bench/o3_graph/test_o2vc_data.py -q`
Expected: 2 passed.

- [ ] **Step 6: Commit**

```bash
git add bench/o3_graph/o2vc_data.py bench/o3_graph/test_o2vc_data.py
git commit -m "o2vc data: longcallR 2.0.0 runner (hifi-masseq, --min-mapq 0) and VCF reader"
```

### Task 3: Tied domain, teachers, read-frame columns, teacher BAM, depth

**Files:**
- Modify: `bench/o3_graph/o2vc_data.py`
- Test: `bench/o3_graph/test_o2vc_data.py`

**Interfaces:**
- Consumes: `o2psv.load_copies`, `o2psv._strand`, `summary.attribute`, `o2correct.span_index` / `overlapping` (load them with the sibling helper).
- Produces:
  - `SLACK = 9`
  - `tied_wide(alns: list[dict]) -> tuple[dict | None, str | None]`: `alns` items are `dict(AS, prim, copy)`. Returns `{copy: its best alignment within SLACK of the best}` and minimap2's primary copy, or `(None, prim)`.
  - `is_teacher(primary_as: int, other_copy_as: list[int]) -> bool`
  - `@dataclass Aln(copy, chrom, start, cigar: tuple, flip: bool, AS, NM, de, ms)`
  - `@dataclass TiedRead(family, copies: dict[str, Aln], prim: str | None, mapq: int, margin: int, seq: str, qual: list[int])`: `seq` and `qual` in the read frame, full length including soft clips.
  - `@dataclass Scan(tied: dict[str, TiedRead], teachers: dict[str, set[str]])`
  - `scan(bam_p: str, copies_p: str, families: set[str] | None = None) -> Scan`
  - `read_frame_positions(start: int, cigar, length: int, flip: bool) -> dict[int, int]`
  - `columns(read: TiedRead, fetch: Callable[[str, int, int], str]) -> list[dict]`: each `dict(p, base, qual, at={copy: (chrom, pos)}, ref={copy: base in read frame})`, only at read positions every tied alignment aligns.
  - `write_teacher_bam(bam_p: str, spans: dict[str, tuple[str, int, int]], teachers: dict[str, set[str]], out_bam: str) -> int`
  - `teacher_depth(teacher_bam: str, chrom: str, positions: Iterable[int]) -> dict[int, int]`

- [ ] **Step 1: Write the failing tests**

**Fixture.** Add `write_bam(path, refs, recs, copies_tsv)` to the test file:
- It writes `refs` ({chrom: seq}) as FASTA with a `.fai` (`pysam.faidx`), a pool-format copies TSV (`copy family chrom start end strand exons`, exons = the span), and a sorted, indexed BAM.
- `recs` items: `dict(name, chrom, start, cigar, seq | None, qual | None, secondary, reverse, AS, NM=0, de=0.0, ms=AS)`.
- A secondary has no SEQ; leave `query_sequence` unset.

**World** (`chrA`, `chrB` = 2,000 random bases, seed 1). Copies:
- F.0 = chrA:100-700 +
- F.1 = chrB:100-700 +
- F.2 = chrB:1000-1600 +
- G.0 = chrA:1000-1600 +

All reads are 600 bp, `600M`. Teachers carry their copy's reference bases.

| read | records | role |
|---|---|---|
| T1–T6 | primary on F.0, AS 600 | teachers of F.0 |
| Y | primary F.0 AS 600, secondary F.1 AS 590 | teacher (margin 10) |
| W | primary F.0 AS 600, secondary G.0 AS 598 | teacher: off-family alignments do not count, as in the corrector |
| X | primary F.0 AS 600, secondary F.1 AS 595 | tied |
| Z | primary F.0 AS 500, secondaries F.1 AS 595 and F.2 AS 590 (reverse) | tied, primary outside the set |
| S | primary F.0, AS 400, only 400 bp long | teacher, but under the caller's 500-bp read-length floor |
| H | primary F.0, AS 600, de 0.08 | teacher, but over the caller's 0.05 divergence limit |

T6 has base quality 5 at chrA:400.

```python
def test_tied_wide_keeps_an_as_tie_whose_primary_lies_below_it():
    alns = [dict(AS=1734, prim=True, copy="F.2"), dict(AS=1819, prim=False, copy="F.5"), dict(AS=1819, prim=False, copy="F.7")]
    per, prim = D.tied_wide(alns)
    assert prim == "F.2" and sorted(per) == ["F.5", "F.7"]
    assert D.tied_wide(alns + [dict(AS=1815, prim=False, copy=None)]) == (None, "F.2")
    assert D.tied_wide([dict(AS=900, prim=True, copy="F.0"), dict(AS=890, prim=False, copy="F.1")]) == (None, "F.0")
    assert D.tied_wide([dict(AS=900, prim=True, copy="F.0"), dict(AS=899, prim=False, copy="G.1")]) == (None, "F.0")


def test_a_teacher_is_a_primary_at_least_ten_as_above_every_other_copy_of_its_family():
    assert D.is_teacher(1000, [990]) and not D.is_teacher(1000, [991]) and D.is_teacher(1000, [])


def test_read_frame_positions_count_both_clip_kinds_skip_insertions_and_introns_and_flip_reverse_alignments():
    assert D.read_frame_positions(100, [(4, 2), (0, 3), (1, 1), (0, 2), (5, 2)], 10, False) == {2: 100, 3: 101, 4: 102, 6: 103, 7: 104}
    assert D.read_frame_positions(100, [(5, 2), (0, 3), (1, 1), (0, 2), (4, 2)], 10, True) == {7: 100, 6: 101, 5: 102, 3: 103, 2: 104}
    assert D.read_frame_positions(100, [(0, 2), (3, 50), (2, 1), (0, 2)], 4, False) == {0: 100, 1: 101, 2: 153, 3: 154}


def test_columns_are_the_read_bases_every_tied_alignment_aligns_with_each_copy_s_base_in_the_read_frame():
    ref = {"chrA": "AAAACGTACG", "chrB": "CGTGCGTTTT"}
    fetch = lambda ch, lo, hi: ref[ch][lo:hi]
    a0 = D.Aln("F.0", "chrA", 4, ((0, 5),), False, 5, 0, 0.0, 5)
    a1 = D.Aln("F.1", "chrB", 1, ((0, 5),), True, 3, 1, 0.2, 3)
    read = D.TiedRead("F", {"F.0": a0, "F.1": a1}, "F.0", 0, 2, "CGTAC", [30] * 5)
    cols = D.columns(read, fetch)
    assert [c["p"] for c in cols] == [0, 1, 2, 3, 4]
    assert cols[2] == dict(p=2, base="T", qual=30, at={"F.0": ("chrA", 6), "F.1": ("chrB", 3)}, ref={"F.0": "T", "F.1": "C"})
    ins = D.TiedRead("F", {"F.0": D.Aln("F.0", "chrA", 4, ((0, 2), (1, 1), (0, 2)), False, 3, 1, 0.2, 3), "F.1": a1},
                     "F.0", 0, 0, "CGTAC", [30] * 5)
    assert [c["p"] for c in D.columns(ins, fetch)] == [0, 1, 3, 4]


def test_scan_finds_the_tied_reads_and_the_teachers(tmp_path):
    bam, copies = make_world(tmp_path)                               # the world above
    s = D.scan(bam, copies)
    assert sorted(s.tied) == ["X", "Z"]
    assert s.tied["Z"].prim == "F.0" and sorted(s.tied["Z"].copies) == ["F.1", "F.2"] and s.tied["Z"].copies["F.2"].flip
    assert s.tied["X"].copies["F.1"].AS == 595 and s.tied["X"].margin == 5 and len(s.tied["X"].seq) == 600
    assert s.teachers["F.0"] == {"T1", "T2", "T3", "T4", "T5", "T6", "Y", "W", "S", "H"}


def test_the_teacher_bam_holds_only_the_teachers_primaries_and_depth_counts_reads_as_the_caller_does(tmp_path):
    bam, copies = make_world(tmp_path)
    s = D.scan(bam, copies)
    out = str(tmp_path / "teach.bam")
    assert D.write_teacher_bam(bam, {"F.0": ("chrA", 100, 700)}, s.teachers, out) == 10
    assert D.teacher_depth(out, "chrA", [110, 400]) == {110: 0, 400: 7}   # T1-T5, Y, W; not S (short), H (de), T6 (qual 5)
```

- [ ] **Step 2: Run the tests and see them fail**

Run: `python3 -m pytest bench/o3_graph/test_o2vc_data.py -q`
Expected: the 6 new tests FAIL; Task 2's 2 still pass.

- [ ] **Step 3: Implement `o2vc_data.py`'s BAM side**

- **`scan`:** two passes over the BAM, like `o2correct.collect()`.
  - Pass 1 builds `dict(AS, prim, copy)` per read, where `copy` comes from `summary.attribute` over the overlapping copies.
  - Pass 2 keeps each tied copy's best record (start, cigar, `flip = rec.is_reverse`, AS, NM, de, ms) and the primary's SEQ and QUAL, reverse-complemented into the read frame when the primary is reverse.
  - Skip a read whose primary SEQ length differs from its CIGAR's query length (both clip kinds counted). Report the count on stderr.
  - Teachers come from pass 1: primary on copy c, not tied, `is_teacher(primary AS, AS on the family's other copies)`.
  - `families` restricts both outputs.
- **`read_frame_positions`:** walks the CIGAR. Ops 0/7/8 map, 1/4/5 advance the read only, 2/3 the reference only. When `flip`, `p = length - 1 - q`.
- **`teacher_depth`:**
  - Amended 2026-10-22 (pre-review correction, controller ruling): reproduces longcallR 2.0.0's hifi-masseq DP from its source (primary, length ≥ 500, float `de` < 0.5, any base quality, no HiFi end trim, poly-run rule needing `fetch`), pinned against the real binary. The filters first written here (de ≤ 0.05, quality ≥ 10, 20-bp trim) were a mistranslation of the usage page.

- [ ] **Step 4: Run the tests and see them pass**

Run: `python3 -m pytest bench/o3_graph/test_o2vc_data.py -q`
Expected: 8 passed.

- [ ] **Step 5: Commit**

```bash
git add bench/o3_graph/o2vc_data.py bench/o3_graph/test_o2vc_data.py
git commit -m "o2vc data: AS-tied domain incl. primaries below the tie, teachers (margin >= 10), read-frame columns, teacher BAM and caller-like depth"
```

### Task 4: `o2vc.py run` and the truth-free case files

**Files:**
- Create: `bench/o3_graph/o2vc.py`
- Test: `bench/o3_graph/test_o2vc.py`

**Interfaces:**
- Consumes: everything Tasks 1–3 produce; `o2psv.substitutions`.
- Produces:
  - `cmd_run(bam_p, copies_p, genome_p, out_dir, families: set[str] | None = None, exclude_chrM: bool = False, longcallr: str = D.LONGCALLR, threads: int = 2) -> int`
  - CLI: `o2vc.py run BAM COPIES.tsv GENOME.fa OUTDIR [--families A,B] [--exclude-chrM] [--threads N]`
  - Files in `OUTDIR`:
    - `calls.tsv`: `read family candidates minimap2 mapq margin outside vc vc_set vc_cols n3 n3_set flag`. `flag` is `"<family>#<id>"` or empty; a flagged read has empty `vc`.
    - `<family>.reads.tsv`: `read copy primary AS NM de ms subs indels ref_for ref_against call_for call_against`.
    - `<family>.columns.tsv`: `read p base qual kind copies`.
      - `kind` is one of `separating_call` (separates under VC but not under N3), `separating_ref` (separates under N3 and VC), `foreign` or `editing`.
      - `copies` is `;`-joined `copy|chrom|pos|ref|genotype|depth` entries.
    - `<family>.flags.tsv`: `flag n_reads reads alleles`.
    - `<family>.teachers.tsv`: `copy teachers records pass lowqual homref multiallelic rnaedit dense teacher_flag_reads`.
    - `vcf/<copy>.vcf`, `teach/<family>.bam`.
    - `run.json`: `longcallr`, `preset`, `min_mapq`, `families`, `depth_agreement: [n with |DP - teacher_depth| <= 1, n records]`.

- [ ] **Step 1: Write the failing test**

**Fixture `make_family_world(tmp_path)`.** It returns an object with `.bam`, `.copies` and `.genome`. Reference bases are
forced at chrA:300 = `A`, :350 = `A` and :420 = `G`. Copies:
- F.0 = chrA:100-700
- F.1 = chrB:100-700, where chrB is identical to chrA, so the reference has no PSV
- F.2 = chrB:1000-1600, with no reads
- M.0 = chrM:100-700 and M.1 = chrA:1000-1600

All reads are 600 bp, `600M`, at 100-700.

| read | records | bases |
|---|---|---|
| 8 teachers on F.0 | primary on F.0 | the reference |
| 8 teachers on F.1 | primary on F.1 | the reference, but `G` at chrB:300 (reference `A` there) |
| X | primary F.0 AS 594, secondary F.1 AS 594 | `G` at 300 |
| Y1, Y2 | primary F.0 AS 588, secondary F.1 AS 588 | `C` at 350 and `T` at 420, both differing from the reference there |
| 2 tied reads | on M.0 / M.1 | — |

**Fake caller.** `fake_caller(prefixes: list, ref: str = "A")` returns a function with `run_longcallr`'s signature. It:
- appends each `out_prefix` to `prefixes`;
- writes the VCF where Task 2 observed longcallR writes it;
- holds one record, only for a `chrB` region: `chrB 301 . {ref} G . PASS . GT:GQ:DP:AF 1/1:30:8:1.00`.

```python
def test_run_places_a_tied_read_by_the_sample_s_call_where_the_reference_alone_cannot_and_flags_a_shared_foreign_pair(tmp_path, monkeypatch):
    w, prefixes = make_family_world(tmp_path), []
    monkeypatch.setattr(O.D, "run_longcallr", fake_caller(prefixes))
    assert O.cmd_run(w.bam, w.copies, w.genome, str(tmp_path / "out"), exclude_chrM=True) == 0
    calls = {r["read"]: r for r in csv.DictReader(open(tmp_path / "out" / "calls.tsv"), delimiter="\t")}
    assert calls["X"]["vc"] == "F.1" and calls["X"]["vc_cols"] == "1" and calls["X"]["n3"] == ""
    assert calls["Y1"]["flag"] == calls["Y2"]["flag"] == "F#0" and calls["Y1"]["vc"] == ""
    assert not any(r["family"] == "M" for r in calls.values())                       # a chrM copy: family left out
    assert not any(p.endswith("F.2") or "/M." in p for p in prefixes)                # F.2 has no teachers: never called
    kinds = {r["kind"] for r in csv.DictReader(open(tmp_path / "out" / "F.columns.tsv"), delimiter="\t") if r["read"] == "X"}
    assert kinds == {"separating_call"}


def test_run_stops_when_a_record_s_ref_disagrees_with_the_genome(tmp_path, monkeypatch):
    w = make_family_world(tmp_path)
    monkeypatch.setattr(O.D, "run_longcallr", fake_caller([], ref="C"))
    with pytest.raises(ValueError, match="chrB:301"):
        O.cmd_run(w.bam, w.copies, w.genome, str(tmp_path / "out"))
```

- [ ] **Step 2: Run the tests and see them fail**

Run: `python3 -m pytest bench/o3_graph/test_o2vc.py -q`
Expected: FAIL (`o2vc.py` does not exist).

- [ ] **Step 3: Implement `cmd_run` and the `run` CLI**

Per family, in sorted order. Skip a family with no tied reads, or one with a chrM copy under `--exclude-chrM`.
1. Write the teacher BAM.
2. Run longcallR once per copy that has ≥ 1 teacher, with `region = f"{chrom}:{max(0, start - 1000)}-{end + 1000}"`.
3. Check every record's REF against the genome (`ValueError(f"{chrom}:{pos0 + 1} REF ...")`).
4. Build each tied read's columns.
5. **VC genotype:** `to_read_frame(genotype(call, call.dp if call else teacher_depth, ref), flip)`.
6. **N3 genotype:** `{ref}` in the read frame wherever VC's is known, else None.
7. Assign under VC and N3, and collect foreign columns (keys in genome strand). Run `flags` per family and clear `vc` for flagged reads.
8. **Teacher side**, reported only:
   - FILTER counts per copy;
   - `flags` over the teachers' own-copy foreign alleles (known positions, base not in the genome-strand genotype).
9. **Tags in `reads.tsv`:** `subs` from `o2psv.substitutions(cigar, NM)`; `indels` from the CIGAR.

- [ ] **Step 4: Run the tests and see them pass**

Run: `python3 -m pytest bench/o3_graph/test_o2vc.py -q`
Expected: 2 passed.

- [ ] **Step 5: Commit**

```bash
git add bench/o3_graph/o2vc.py bench/o3_graph/test_o2vc.py
git commit -m "o2vc run: per-family longcallR on teachers, VC and reference-only (N3) placements, missing-copy flags, truth-free case files"
```

### Task 5: `score` and `score-worlds`

**Files:**
- Modify: `bench/o3_graph/o2vc.py`
- Test: `bench/o3_graph/test_o2vc.py`

**Interfaces:**
- Consumes: `calls.tsv`, `<family>.reads.tsv` and `<family>.columns.tsv` from Task 4; `o2psv._truth`.
- Produces:
  - `cmd_score(out_dir, truth_ps: list[str], mode: str, n4_p=None, panel=None, focus=None, exclude_p=None) -> int`. `mode` is one of:
    - `donor`: `class`/`truth` tables;
    - `world`: `read`/`copy`;
    - `units`: `read`/`unit`, the gorilla truth;
    - `none`: no truth, DAZ.
  - Writes `OUTDIR/score[.<panel>[.<focus>]].json`, plus `<family>.case.tsv` and `<family>.summary.md` for every family scored.
  - `cmd_score_worlds(out_p, dirs: list[str]) -> int`: each dir holds `o2vc/calls.tsv`, `reads.truth.E0.tsv` and `o2correct.tsv`; read names are prefixed by the dir.
  - `evidence(kinds: set[str]) -> str`: `"sample call"` when any deciding column is `separating_call`, else `"reference PSV"`.
  - CLI:
    - `o2vc.py score OUTDIR TRUTH... --mode M [--n4 O2CORRECT.tsv] [--panel A|C|Cg|D|E] [--focus FAMILY] [--exclude PAR.json]`
    - `o2vc.py score-worlds OUT.json DIR...`

**Scoring rules (from the spec):**
- **Decidable:**
  - `donor`: class ∈ {unique, unique_elsewhere} and truth ∈ candidates;
  - `world`: truth ∈ candidates;
  - `units`: no placement truth.
- **Arms per read:**
  - N1 = `minimap2`;
  - N2 = N1 when `mapq >= 1`, else abstain;
  - N3 = `n3`;
  - N4 = the corrector's pick from `--n4` (N1 when the read is absent there);
  - VC = `vc`.
- **Per arm:**
  - `placed`, `right`, `wrong` and `abstained` over decidable reads;
  - `fixes` (placed right, N1 wrong), `breaks` (placed wrong, N1 right) and `both_wrong`.
- **Calibration and flags:**
  - **K1** (`donor` only): VC placements among tied `own_tie` reads.
  - **K2:** flagged tied reads whose truth copy is present (classes unique, unique_elsewhere, own_tie; every read in `world`), over all such tied reads.
  - **C1:** flags with ≥ 2 `no_counterpart` reads, with each flag's class counts listed.
  - **C1g:** flags with ≥ 2 reads of unit `C`.
  - **D2:** decidable reads N1 places right that VC moves or flags.
- **Bars:**
  - **A** (on `--focus`): A1 VC wrong/placed ≤ 0.02, A2 fixes ≥ 10, A3 VC fixes > N3 fixes, K1 ≤ 0.02 and K2 ≤ 0.02 (pooled, excluded families out).
  - **C:** C1 ≥ 1.
  - **Cg:** C1g ≥ 1.
  - **D:** wrong/placed ≤ 0.02 and D2 ≤ 0.02.
  - **E:** E1 ≤ 0.02, E2 fixes > breaks, K2 ≤ 0.02.
  - Families in `--exclude` are scored under `"excluded"` only.

- [ ] **Step 1: Write the failing tests**

**Toy run directory.** Family F, candidates `F.0,F.1`. Truth is given as `class/truth`.

| read | minimap2 | mapq | vc | n3 | flag | truth |
|---|---|---|---|---|---|---|
| r1 | F.0 | 0 | F.1 | | | unique / F.1 |
| r2 | F.0 | 3 | F.1 | F.1 | | unique / F.0 |
| r3 | F.1 | 0 | | | | unique / F.1 |
| r4 | F.0 | 0 | F.0 | | | own_tie / — |
| r5 | F.0 | 0 | | | F#0 | no_counterpart |
| r6 | F.0 | 0 | | | F#0 | no_counterpart |
| r7 | F.0 | 0 | | | F#1 | unique / F.0 |

Write a family G (one decidable read, VC right) into `par.json`.

```python
def test_score_counts_each_arm_on_decidable_reads_and_the_calibration_bars(tmp_path):
    d = make_score_toy(tmp_path)
    assert O.cmd_score(d.out, [d.truth], "donor", panel="A", focus="F", exclude_p=d.par) == 0
    s = json.load(open(os.path.join(d.out, "score.A.F.json")))
    assert s["arms"]["VC"] == dict(placed=2, right=1, wrong=1, abstained=2, fixes=1, breaks=1, both_wrong=0)
    assert s["arms"]["N2"]["placed"] == 1 and s["arms"]["N2"]["right"] == 1
    assert s["K1"] == dict(reads=1, placed=1) and s["K2"] == dict(reads=5, flagged=1)
    assert [f["classes"] for f in s["flags"] if f["id"] == "F#0"] == [{"no_counterpart": 2}]
    assert s["bars"]["A1"]["pass"] is False and s["bars"]["A3"]["pass"] is True and s["bars"]["K1"]["pass"] is False
    assert "G" in s["excluded"] and "G" not in s["families"]


def test_evidence_names_a_sample_call_when_any_deciding_column_is_one():
    assert O.evidence({"separating_call", "separating_ref"}) == "sample call"
    assert O.evidence({"separating_ref"}) == "reference PSV"


def test_score_worlds_pools_worlds_keeping_their_read_names_apart(tmp_path):
    dirs = make_world_toys(tmp_path)          # two worlds, each with read "s1" tied F.0/F.1: VC right in one, wrong in the other
    assert O.cmd_score_worlds(str(tmp_path / "d.json"), dirs) == 0
    s = json.load(open(tmp_path / "d.json"))
    assert s["arms"]["VC"]["placed"] == 2 and s["arms"]["VC"]["wrong"] == 1
```

- [ ] **Step 2: Run the tests and see them fail**

Run: `python3 -m pytest bench/o3_graph/test_o2vc.py -q`
Expected: the 3 new tests FAIL.

- [ ] **Step 3: Implement `cmd_score`, `cmd_score_worlds`, `evidence` and the case files**

- **`<family>.case.tsv`:** `reads.tsv` plus `truth_class`, `truth` and `vc_outcome` (fix, break, both right, both wrong, abstain, flagged).
- **`<family>.summary.md`:**
  - the arm table;
  - how many VC placements rest on a single column (`vc_cols == 1`);
  - for each fix, its `evidence(...)` and whether AS (max), NM (min) or subs (min) alone picks the truth, with ties counting as no pick, as in `o2tags.pick`;
  - three worked reads with their `columns.tsv` rows: the first fix, the first abstention and the first flagged read, by read name.

- [ ] **Step 4: Run the tests and see them pass**

Run: `python3 -m pytest bench/o3_graph/test_o2vc_rule.py bench/o3_graph/test_o2vc_data.py bench/o3_graph/test_o2vc.py -q`
Expected: 23 passed.

- [ ] **Step 5: Commit**

```bash
git add bench/o3_graph/o2vc.py bench/o3_graph/test_o2vc.py
git commit -m "o2vc score: arms N1-N4 and VC on decidable reads, calibration K1/K2, flag truth, panel bars, case files and summaries"
```

### Task 6: Development run on SMN2 and MSTRG3150

**Files:**
- Modify only through TDD fixes (each bug found gets a failing test first, in the file of the code it fixes).

- [ ] **Step 1: Run the development families**

```bash
cd /mnt/linuxdisk/home/juanfraitu/rustle_soto_subset
T=/mnt/linuxdisk/tmp/o3_graph; G=/mnt/linuxdisk/home/juanfraitu/npip_ladder/idx/target.fa
RLOCK_TIMEOUT=5400 bash tools/rlock.sh heavy python3 bench/o3_graph/o2vc.py run $T/hg002proof/chm13.bam $T/o2psv_dev/pool_copies.tsv $G $T/o2vc/dev --families S19SMN2,S19MSTRG3150 --exclude-chrM
```

Expected: exit 0; `calls.tsv` holds both families' tied reads; `run.json` reports longcallR 2.0.0.

- [ ] **Step 2: Check depth agreement**

Read `run.json` `depth_agreement`. Expected: ≥ 95% of records have |DP − teacher_depth| ≤ 1. If not, align `teacher_depth`'s filters with longcallR's source (`src/` of v2.0.0) under TDD, then rerun Step 1.

- [ ] **Step 3: Score the development run**

```bash
python3 bench/o3_graph/o2vc.py score $T/o2vc/dev $T/hg002proof/truth.tsv $T/hg002proof/truth_untied.tsv --mode donor --n4 $T/o2psv_dev/hg002.v2.tsv
```

Expected: completes. Then check by hand:
- every VC placement has `vc_cols >= 1`;
- no flagged read is placed;
- `S19SMN2.columns.tsv` shows the exon-7 c.840 site.

Write the dev numbers (arms, K1, K2, flags) into `docs/O2VC_DEV_2026-10-22.md`. They are quoted in the prereg as already seen.

- [ ] **Step 4: Commit**

```bash
git add docs/O2VC_DEV_2026-10-22.md bench/o3_graph/
git commit -m "o2vc: development run on SMN2 and MSTRG3150 (HG002), the numbers seen before the prereg"
```

### Task 7: Panel inputs

**Files:**
- Modify: `bench/o3_graph/o2vc.py` (two converters)
- Test: `bench/o3_graph/test_o2vc.py`

**Interfaces:**
- Produces:
  - `copies_from_sets(fam_json: str, set_name: str, family: str) -> list[dict]`: copy `f"{family}.{k}"` in the JSON's order, with chrom, start, end and strand. Exons are the sorted, de-duplicated union of its transcripts' exons.
  - `copies_from_genes(genes_tsv: str, family: str) -> list[dict]`: rows `NPIP:i chrom strand s-e` become copy `f"{family}.{i}"`, with the span taken from the exons.
  - CLI:
    - `o2vc.py copies-from-sets FAM.json SET FAMILY OUT.tsv`
    - `o2vc.py copies-from-genes GENES.tsv FAMILY OUT.tsv`
  - Both write the pool format: `copy family chrom start end strand exons`.

- [ ] **Step 1: Write the failing tests**

```python
def test_copies_from_sets_writes_pool_rows_with_the_union_of_transcript_exons(tmp_path):
    p = tmp_path / "f.json"
    p.write_text(json.dumps({"sets": {"seq19": [{"family": "S19X", "copies": [
        {"chrom": "chrY", "start": 10, "end": 90, "strand": "-", "transcripts": {"t1": [[10, 20], [50, 60]], "t2": [[10, 20], [70, 90]]}}]}]}}))
    assert O.copies_from_sets(str(p), "seq19", "S19X") == [dict(copy="S19X.0", family="S19X", chrom="chrY", start=10, end=90,
                                                                strand="-", exons="10-20,50-60,70-90")]


def test_copies_from_genes_names_copies_by_their_index(tmp_path):
    p = tmp_path / "g.tsv"
    p.write_text("gene_id\tchrom\tstrand\texons\nNPIP:0\tNC_1\t+\t100-900\nNPIP:1\tNC_2\t-\t50-60,70-80\n")
    assert [r["copy"] for r in O.copies_from_genes(str(p), "GGONPIP")] == ["GGONPIP.0", "GGONPIP.1"]
    assert O.copies_from_genes(str(p), "GGONPIP")[1]["start"] == 50 and O.copies_from_genes(str(p), "GGONPIP")[1]["end"] == 80
```

- [ ] **Step 2: Run the tests, implement the converters, and see the tests pass**

Run: `python3 -m pytest bench/o3_graph/test_o2vc.py -q`. Expected after the implementation: all pass.

- [ ] **Step 3: Build the panel inputs**

These produce outputs only; no panel is run.

```bash
T=/mnt/linuxdisk/tmp/o3_graph; I=$T/o2vc/inputs; mkdir -p $I
python3 bench/o3_graph/o2vc.py copies-from-sets $T/r19seq/families19cal_full2_all.json seq19 S19DAZ1 $I/daz_copies.tsv
python3 bench/o3_graph/o2vc.py copies-from-genes /mnt/linuxdisk/tmp/o3_npip/copies_genes.tsv GGONPIP $I/npip_ggo_copies.tsv
awk -F'\t' 'NR>1{print $3"\t"$4"\t"$5}' $I/npip_ggo_copies.tsv > $I/npip_spans.bed
L=/mnt/linuxdisk/home/juanfraitu/fibroblasts/GCA_029281585.2_flnc_mm.bam
RLOCK_TIMEOUT=5400 bash tools/rlock.sh heavy bash -c "samtools view -@4 -L $I/npip_spans.bed $L | cut -f1 | sort -u > $I/npip_names.txt && samtools view -@4 -b -N $I/npip_names.txt $L | samtools sort -@4 -o $I/npip.bam - && samtools index $I/npip.bam"
```

Expected:
- `daz_copies.tsv` has 4 copies (DAZ1–4).
- `npip_ggo_copies.tsv` has 25 rows.
- `npip.bam` holds every record of the reads touching NPIP, including their alignments elsewhere.
- `npip_names.txt` contains the unit-C reads of `/mnt/linuxdisk/tmp/o3_npip/realign/units.tsv` (39 reads).

- [ ] **Step 4: Commit**

```bash
git add bench/o3_graph/o2vc.py bench/o3_graph/test_o2vc.py
git commit -m "o2vc: copy-table converters for the DAZ (family sets) and gorilla NPIP (gene table) panels"
```

### Task 8: Pre-registration

**Files:**
- Create: `docs/PREREG_o2_variant_calls_2026-10-22.md`

- [ ] **Step 1: Write the prereg**

It contains:
- **The user's words** (2026-10-22): "it just needs to be something that passes from beginning to end … show that we got the right placement instead of minimap2 or that we found a possibly missing copy based on the variants".
- **The frozen code commit** (`git rev-parse HEAD`).
- **The rule, verbatim** from the spec's method.
- **Arms** N1–N4 and VC.
- **Panels with their exact input paths and Task 9's commands.**
- **The bars, verbatim** from the spec's pass criteria, including the 2026-10-22 amendments: K1 on HG002 only, PAR apart, chrM out.
- **What was seen before:** `docs/O2VC_DEV_2026-10-22.md` (SMN2, MSTRG3150) and the corrector's HG002 per-family numbers (GTF2IRD2 19/0, NBPF1 11/0, PDPR 14/0), which chose the showcase.
- **What was never run:** every other panel.
- **The decision rule:**
  - the showcase passes iff A1–A3, K1 and K2 hold;
  - C-human, C-gorilla, D and E each get their own verdict, never pooled;
  - a failed GTF2IRD2 is a FAIL; the backups are reported, never substituted.

- [ ] **Step 2: Commit**

```bash
git add docs/PREREG_o2_variant_calls_2026-10-22.md
git commit -m "Pre-register O2 from the sample's own variant calls: showcase GTF2IRD2 (HG002), DAZ, missing copy (HG002, gorilla), one genome, HG00733 held out"
```

### Task 9: Run every panel once and score it

Run each command once, in this order, serially. Never rerun a panel after seeing its score. A crash gets its fix under TDD and a prereg amendment that says so, before the rerun.

- [ ] **Step 1: HG002, every family (panels A and C-human, calibration)**

```bash
T=/mnt/linuxdisk/tmp/o3_graph; G=/mnt/linuxdisk/home/juanfraitu/npip_ladder/idx/target.fa; P=$T/famident/par_families.json
RLOCK_TIMEOUT=5400 bash tools/rlock.sh heavy python3 bench/o3_graph/o2vc.py run $T/hg002proof/chm13.bam $T/o2psv_dev/pool_copies.tsv $G $T/o2vc/hg002 --exclude-chrM
for F in S19GTF2IRD2 S19NBPF1 S19PDPR; do
  python3 bench/o3_graph/o2vc.py score $T/o2vc/hg002 $T/hg002proof/truth.tsv $T/hg002proof/truth_untied.tsv --mode donor --n4 $T/o2psv_dev/hg002.v2.tsv --panel A --focus $F --exclude $P
done
python3 bench/o3_graph/o2vc.py score $T/o2vc/hg002 $T/hg002proof/truth.tsv $T/hg002proof/truth_untied.tsv --mode donor --n4 $T/o2psv_dev/hg002.v2.tsv --panel C --exclude $P
```

Expected: `score.A.S19GTF2IRD2.json` (the verdict), the backups' JSONs and `score.C.json`.

- [ ] **Step 2: DAZ (panel B, descriptive)**

```bash
python3 bench/o3_graph/o2vc.py run $T/yag_testis/yag.bam $T/o2vc/inputs/daz_copies.tsv $G $T/o2vc/daz --families S19DAZ1
python3 bench/o3_graph/o2vc.py score $T/o2vc/daz --mode none
```

- [ ] **Step 3: Gorilla NPIP (panel C-gorilla)**

```bash
RLOCK_TIMEOUT=5400 bash tools/rlock.sh heavy python3 bench/o3_graph/o2vc.py run $T/o2vc/inputs/npip.bam $T/o2vc/inputs/npip_ggo_copies.tsv /mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta $T/o2vc/ggo_npip
python3 bench/o3_graph/o2vc.py score $T/o2vc/ggo_npip /mnt/linuxdisk/tmp/o3_npip/realign/units.tsv --mode units --panel Cg
```

- [ ] **Step 4: One genome (panel D, the 18 round-23 worlds)**

```bash
for W in $T/fast/heldout23/*_s13*/{INTACT,CROSS,SILENT}; do
  RLOCK_TIMEOUT=1800 bash tools/rlock.sh light python3 bench/o3_graph/o2vc.py run $W/E0.bam $W/o2psv.copies.tsv $W/genome.ref.fa $W/o2vc || echo "CRASH $W"
done
python3 bench/o3_graph/o2vc.py score-worlds $T/o2vc/d_pooled.json $T/fast/heldout23/*_s13*/{INTACT,CROSS,SILENT}
```

Report the pooled numbers and each world type separately: INTACT and SILENT are subsets of CROSS (register 26 audit).

- [ ] **Step 5: HG00733, held out (panel E)**

```bash
S=$T/donorsim/HG00733/sim
RLOCK_TIMEOUT=3600 bash tools/rlock.sh heavy python3 bench/o3_graph/o2vc.py run $S/chm13.bam $T/o2psv_dev/pool_copies.tsv $G $T/o2vc/hg00733
python3 bench/o3_graph/o2vc.py score $T/o2vc/hg00733 $S/truth.tsv --mode donor --n4 $S/second/o2correct.tsv --panel E --exclude $P
```

- [ ] **Step 6: Write the results into the prereg and commit**

Append `## Results (code <hash>; run once)` to the prereg:
- one table per panel, with each bar's value and verdict;
- the GTF2IRD2 summary's worked reads.

```bash
git add docs/PREREG_o2_variant_calls_2026-10-22.md
git commit -m "O2 variant calls: results, run once (verdicts per panel)"
```

### Task 10: Audit, register, memory, report

- [ ] **Step 1: Freeze and audit**

Copy `$T/o2vc/` (without BAMs) to `$T/o2vc/audit_frozen_<hash>/` and make it read-only. Dispatch the `metric-auditor` agent with:
- the claim (each panel's verdict);
- `docs/PREREG_o2_variant_calls_2026-10-22.md`;
- `bench/o3_graph/o2vc.py` as the scorer;
- the frozen directory.

Fold its findings into the Results as "Audit" notes, without changing a verdict unless the auditor shows a scoring bug. A scoring bug is fixed under TDD and noted as an amendment.

- [ ] **Step 2: Register row**

Append row 1374 to `docs/NEGATIVE_RESULTS_REGISTER.md`. Use the 5-column format `| # | date | what was tried | substrate | verdict |`, with the verdicts and the audit's caveats.

- [ ] **Step 3: Commit**

```bash
git add docs/NEGATIVE_RESULTS_REGISTER.md docs/PREREG_o2_variant_calls_2026-10-22.md
git commit -m "Register 1374: O2 from the sample's own variant calls (audited)"
```

- [ ] **Step 4: Memory and report**

- Update `project_o2_variant_calls_2026-10-22.md` and its `MEMORY.md` line with the verdicts and commits.
- Report to the user:
  - each panel's verdict with its numbers;
  - for GTF2IRD2, the fixes, breaks and abstentions;
  - which evidence decided the fixes;
  - whether the tags alone would have;
  - one worked read.
