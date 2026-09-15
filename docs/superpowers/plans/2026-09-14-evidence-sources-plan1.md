# Evidence Sources (Plan 1 of 2) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Turn every optional input into one internal format, and measure on gorilla what each missing source costs.
- Duplication evidence (D1 supplied SD calls, D2 BISER, D3 minimap2 self-alignment, D4 meryl k-mer copy number) becomes a
  pair table / copy-number track.
- Repeat evidence (R1 RepeatMasker, R2 soft-mask, R3 WindowMasker+DustMasker, R4 meryl high-copy) becomes repeat intervals.
- Build the ape substrates (human-lifted orthologous chromosomes).

**Architecture:** Three new Python modules in `bench/`:
- `dup_evidence.py` for D1–D4;
- `repeat_evidence.py` for R1–R4 and C_max;
- `ape_substrate.py` for Liftoff batches, the orthologous substrate, subsets and manifests.

Plus `evidence_compare.py`, the source-cost report. Each converts a tool's output into a fixed TSV/BED format with a `source`
column. `bench/dna_sd_atoms.py` stays unchanged: pair tables are written back as gorilla-SEDEF-format rows so its CIGAR mode
consumes any pair source. Plan 2 (layers G/P/E, arms, pre-registration AO, orangutan) is written after Task 7's report,
because that report fixes the hold-out source defaults.

**Tech Stack:** Python 3 (pysam, pytest), minimap2, BISER (`envs/biser`), meryl + meryl-lookup (`envs/phasing_eval`),
WindowMasker/DustMasker (`envs/blast`), Liftoff (`envs/liftoff`), samtools.

**Spec:** `docs/superpowers/specs/2026-09-14-evidence-agnostic-family-definition-design.md` (sections 1a, 2 "Shared
preparation", "Evidence-source comparison", 4).

## Global Constraints

- WSL2: long work runs as FOREGROUND calls of <= 10 min (`timeout 590`), resumable, one heavy run at a time. No `nohup`/
  background drivers, no `pkill -f`; kill by PID after `readlink /proc/<pid>/cwd`.
- Large outputs go under `/mnt/linuxdisk/home/juanfraitu/o1_falsemerge/lit/evid/<species>/`, never under `/` or `winloci_scratch`.
- Tool paths:
  - `BISER=/home/juanfra/miniforge3/envs/biser/bin/biser`
  - `MERYL=/home/juanfra/miniforge3/envs/phasing_eval/bin/meryl`
  - `MERYL_LOOKUP=/home/juanfra/miniforge3/envs/phasing_eval/bin/meryl-lookup`
  - `WINDOWMASKER=/home/juanfra/miniforge3/envs/blast/bin/windowmasker`
  - `DUSTMASKER=/home/juanfra/miniforge3/envs/blast/bin/dustmasker`
  - `LIFTOFF=/home/juanfra/miniforge3/envs/liftoff/bin/liftoff`
  - `minimap2`, `samtools` on PATH.
- Genomes:
  - human `/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa`, with RefSeq
    `/mnt/linuxdisk/home/juanfraitu/winloci_data/Reference/chm13v2.0_RefSeq_full.gff.gz` (tabix-indexed);
  - gorilla `/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta` + `winloci_data/GGO_genomic.gff`;
  - orangutan `winloci_data/GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.fna` + `winloci_data/PPY_genomic.gff`;
  - chimpanzee `winloci_data/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna` + `winloci_data/PTR_genomic.gff`.
- Gorilla-only sources:
  - D1 `winloci_data/GGO_sedef_final.bed` (gorilla SEDEF format, CIGAR col 32);
  - R1 `winloci_data/rmsk/GCF_029281585.2.repeatMasker.out.gz`.
- **Orangutan is the hold-out: no task in this plan may read `PPY_genomic.gff` or run anything on the orangutan genome
  except Task 5's Liftoff placement of HUMAN genes and substrate selection** (those read no native annotation). Chimpanzee
  is report-only.
- SD definition constants (D3): aligned >= 1000 bp, identity >= 0.90, and >= 1000 aligned bases outside the repeat intervals.
- Never pool species numbers. Commit after each task, ending messages with
  `Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>`.
- Tests: `python3 -m pytest bench/test_<module>.py -q` from `/mnt/c/Users/jfris/Desktop/Rustle`.

---

### Task 1: Pair table and D1/D2/D3 format adapters

**Files:**
- Create: `bench/dup_evidence.py`
- Test: `bench/test_dup_evidence.py`

**Interfaces:**
- Produces:
  - `PAIR_COLS` (tuple of 11 names);
  - `from_sedef(line: str, fmt: str) -> tuple | None`;
  - `from_biser(line: str) -> tuple | None`;
  - `from_selfpaf(line: str) -> tuple | None` (chunk query names `chrom@offset`);
  - `normalize_cigar(cg: str) -> str`;
  - `to_sedef_gorilla(pair: tuple) -> str`;
  - `write_pairs(pairs, path)`, `read_pairs(path) -> list[tuple]`.
- Pair tuple order: `(chrom_a, start_a, end_a, chrom_b, start_b, end_b, strand_a, strand_b, identity, cigar, source)`.
  - Coordinates are 0-based half-open.
  - CIGAR uses `M/I/D` only, with SEDEF semantics: M consumes both, D consumes side A, I consumes side B, and side B is
    reverse-complemented when `strand_b == "-"`.

- [ ] **Step 1: Write the failing tests**

```python
# bench/test_dup_evidence.py
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dup_evidence as de  # noqa: E402

SEDEF_GGO = "\t".join(["NC_011120.1", "1911", "7456", "NC_073224.2", "125568083", "125573628", "S", "5.2", "+", "-",
                       "5545", "5547", "m=5.2;g=0.1", "2", "2", "5543", "5256", "287", "259", "28", "0.948223",
                       "0.947539", "0.053651", "0.0543876", "4", "5545", "5301", "5023", "5256", "287", "4", "4",
                       "673M1I2734M1I328M1D1289M1D519M", "0.947539"])
BISER = "\t".join(["chrA", "99997", "120000", "chrB", "59997", "80000", "toy:toy", "3.7", "+", "+", "20003", "20003",
                   "20003M", "X=3.7;ID=0"])


def test_from_sedef_gorilla():
    p = de.from_sedef(SEDEF_GGO, "gorilla")
    assert p[:8] == ("NC_011120.1", 1911, 7456, "NC_073224.2", 125568083, 125573628, "+", "-")
    assert abs(p[8] - 0.948223) < 1e-9 and p[9] == "673M1I2734M1I328M1D1289M1D519M" and p[10] == "D1:sedef"


def test_from_biser_identity_from_error():
    p = de.from_biser(BISER)
    assert p[:8] == ("chrA", 99997, 120000, "chrB", 59997, 80000, "+", "+")
    assert abs(p[8] - 0.963) < 1e-9 and p[9] == "20003M" and p[10] == "D2:biser"


def test_normalize_cigar_merges_eq_x():
    assert de.normalize_cigar("10=2X5=3I4=1D") == "17M3I4M1D"


def paf(q, qlen, qs, qe, strand, t, tlen, ts, te, nm, bl, cg):
    return "\t".join(map(str, [q, qlen, qs, qe, strand, t, tlen, ts, te, nm, bl, 60, f"cg:Z:{cg}"]))


def test_selfpaf_side_a_is_target_and_filters():
    # chunk chr1@1000 query 0-2000 aligned to chr2 5000-7000, identity 0.95, canonical (chr1 > chr2 so target is side A)
    p = de.from_selfpaf(paf("chr2@1000", 5000, 0, 2000, "+", "chr1", 9000, 5000, 7000, 1900, 2000, "2000M"))
    assert p == ("chr1", 5000, 7000, "chr2", 1000, 3000, "+", "+", 0.95, "2000M", "D3:selfaln")
    # self hit (same chrom, overlapping) dropped; short dropped; low identity dropped
    assert de.from_selfpaf(paf("chr1@0", 9000, 5000, 7000, "+", "chr1", 9000, 5000, 7000, 2000, 2000, "2000M")) is None
    assert de.from_selfpaf(paf("chr2@0", 5000, 0, 900, "+", "chr1", 9000, 0, 900, 900, 900, "900M")) is None
    assert de.from_selfpaf(paf("chr2@0", 5000, 0, 2000, "+", "chr1", 9000, 0, 2000, 1700, 2000, "2000M")) is None
    # non-canonical orientation (query side sorts first) dropped: the mirror record carries the pair
    assert de.from_selfpaf(paf("chr1@0", 9000, 5000, 7000, "+", "chr2", 5000, 1000, 3000, 1900, 2000, "2000M")) is None


def test_to_sedef_gorilla_roundtrip_identity_check():
    p = ("chrA", 0, 1000, "chrB", 10, 1010, "+", "-", 0.9123, "1000M", "D2:biser")
    f = de.to_sedef_gorilla(p).split("\t")
    m, mm, frac = float(f[16]), float(f[17]), float(f[20])
    assert abs(m / (m + mm) - frac) <= 1e-4 and f[32] == "1000M" and f[8] == "+" and f[9] == "-"
    assert de.from_sedef("\t".join(f), "gorilla")[:8] == p[:8]


def test_pairs_file_roundtrip(tmp_path):
    p = ("chrA", 0, 1000, "chrB", 10, 1010, "+", "-", 0.9123, "1000M", "D2:biser")
    de.write_pairs([p], tmp_path / "x.pairs.tsv")
    assert de.read_pairs(tmp_path / "x.pairs.tsv") == [p]
```

- [ ] **Step 2: Run the tests and verify they fail**

Run: `python3 -m pytest bench/test_dup_evidence.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'dup_evidence'`

- [ ] **Step 3: Write the implementation**

```python
#!/usr/bin/env python3
"""Evidence-agnostic families (spec 2026-09-14 §1a): duplication evidence D1-D4 in one internal format.

Pair table (TSV, header = PAIR_COLS): 0-based half-open intervals of side A and side B, strands, identity, CIGAR (M/I/D;
M both sides, D side A, I side B; side B reverse-complemented when strand_b is '-'), source (D1:sedef, D2:biser, D3:selfaln).
"""
import csv
import re

PAIR_COLS = ("chrom_a", "start_a", "end_a", "chrom_b", "start_b", "end_b", "strand_a", "strand_b", "identity", "cigar", "source")
MIN_SD_BP, MIN_SD_ID = 1000, 0.90


def normalize_cigar(cg):
    ops = [(int(n), "M" if o in "=X" else o) for n, o in re.findall(r"(\d+)([MIDNSHP=X])", cg)]
    out = []
    for n, o in ops:
        if out and out[-1][1] == o:
            out[-1] = (out[-1][0] + n, o)
        else:
            out.append((n, o))
    return "".join(f"{n}{o}" for n, o in out)


def from_sedef(line, fmt):
    f = line.rstrip("\n").split("\t")
    if len(f) < 23 or not f[1].isdigit():
        return None
    if fmt == "gorilla":
        ident, cig = float(f[20]), (f[32] if len(f) > 32 else "")
    elif fmt == "human":
        ident, cig = float(f[22]), ""
    else:
        raise ValueError(fmt)
    return (f[0], int(f[1]), int(f[2]), f[3], int(f[4]), int(f[5]), f[8], f[9], ident, cig, "D1:sedef")


def from_biser(line):
    f = line.rstrip("\n").split("\t")
    if len(f) < 13 or not f[1].isdigit():
        return None
    m = re.search(r"X=([\d.]+)", f[13]) if len(f) > 13 else None
    err = float(m.group(1)) if m else float(f[7])
    return (f[0], int(f[1]), int(f[2]), f[3], int(f[4]), int(f[5]), f[8], f[9], round(1 - err / 100, 9),
            normalize_cigar(f[12]), "D2:biser")


def from_selfpaf(line):
    """minimap2 self-alignment record of a chunk (query `chrom@offset`) against the genome. Side A = target (forward),
    side B = query with the record's strand, so minimap2's CIGAR (D = target, I = query) already has SEDEF semantics.
    Keeps canonical non-self records >= MIN_SD_BP at identity >= MIN_SD_ID."""
    f = line.rstrip("\n").split("\t")
    qc, off = f[0].rsplit("@", 1)
    qs, qe = int(off) + int(f[2]), int(off) + int(f[3])
    tc, ts, te = f[5], int(f[7]), int(f[8])
    if qc == tc and qs < te and ts < qe:
        return None
    nm, bl = int(f[9]), int(f[10])
    if bl < MIN_SD_BP or nm / bl < MIN_SD_ID:
        return None
    if (tc, ts) > (qc, qs):
        return None
    cg = next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), "")
    return (tc, ts, te, qc, qs, qe, "+", f[4], round(nm / bl, 9), normalize_cigar(cg), "D3:selfaln")


def to_sedef_gorilla(p):
    ca, a1, a2, cb, b1, b2, sa, sb, ident, cig, _src = p
    L = sum(int(n) for n, o in re.findall(r"(\d+)([MID])", cig) if o == "M") or max(a2 - a1, b2 - b1)
    m = round(ident * L)
    mm = L - m
    frac = m / (m + mm) if L else 0.0
    row = [""] * 34
    row[0:6] = [ca, str(a1), str(a2), cb, str(b1), str(b2)]
    row[6], row[7], row[8], row[9] = "S", f"{(1 - ident) * 100:.1f}", sa, sb
    row[16], row[17], row[20], row[32], row[33] = str(m), str(mm), f"{frac:.6f}", cig, f"{frac:.6f}"
    return "\t".join(row)


def write_pairs(pairs, path):
    with open(path, "w") as fh:
        fh.write("\t".join(PAIR_COLS) + "\n")
        for p in pairs:
            fh.write("\t".join(map(str, p)) + "\n")


def read_pairs(path):
    out = []
    for r in csv.DictReader(open(path), delimiter="\t"):
        out.append((r["chrom_a"], int(r["start_a"]), int(r["end_a"]), r["chrom_b"], int(r["start_b"]), int(r["end_b"]),
                    r["strand_a"], r["strand_b"], float(r["identity"]), r["cigar"], r["source"]))
    return out
```

- [ ] **Step 4: Run the tests and verify they pass**

Run: `python3 -m pytest bench/test_dup_evidence.py -q`
Expected: `6 passed`

- [ ] **Step 5: Commit**

```bash
git add bench/dup_evidence.py bench/test_dup_evidence.py
git commit -m "Evidence sources: pair table and SEDEF/BISER/self-alignment adapters (spec 2026-09-14 §1a)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 2: Repeat evidence R1–R3 in one interval format

**Files:**
- Create: `bench/repeat_evidence.py`
- Test: `bench/test_repeat_evidence.py`

**Interfaces:**
- Produces:
  - `parse_rmsk(path, contigs=None) -> list[(chrom, start, end, cls)]` (R1; `.out` or `.out.gz`);
  - `lowercase_runs(fasta, contigs=None) -> list[(chrom, start, end, ".")]` (R2);
  - `parse_masker_intervals(path) -> list[(chrom, start, end, ".")]` (WindowMasker/DustMasker `-outfmt interval`, 0-based
    inclusive → half-open);
  - `merge(ivs) -> dict[chrom, list[(s, e)]]`;
  - `write_bed(ivs, path, source)`, `read_bed(path) -> list[(chrom, start, end, cls)]`;
  - `masked_bases(merged, chrom, s, e) -> int`.
- BED columns: `chrom start end class source`. Source is `R1:rmsk` / `R2:softmask` / `R3:windowmasker+dust` / `R4:meryl`.

- [ ] **Step 1: Write the failing tests**

```python
# bench/test_repeat_evidence.py
import gzip
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import repeat_evidence as re_  # noqa: E402

RMSK = """   SW  perc perc perc  query      position in query           matching       repeat              position in  repeat
score  div. del. ins.  sequence    begin     end    (left)    repeat         class/family         begin  end (left)   ID

  311  32.7  3.2  3.9  chrA         11      20   (14400) +  AluY           SINE/Alu              1  155 (1128) 1
  650   5.1  0.0  0.0  chrB          1      78    (13839) +  (CA)n          Simple_repeat          1   78    (0) 2
"""


def test_parse_rmsk_plain_and_gz(tmp_path):
    p = tmp_path / "x.out"
    p.write_text(RMSK)
    g = tmp_path / "x.out.gz"
    with gzip.open(g, "wt") as fh:
        fh.write(RMSK)
    want = [("chrA", 10, 20, "SINE/Alu"), ("chrB", 0, 78, "Simple_repeat")]
    assert re_.parse_rmsk(p) == want and re_.parse_rmsk(g) == want and re_.parse_rmsk(p, {"chrA"}) == want[:1]


def test_lowercase_runs(tmp_path):
    fa = tmp_path / "t.fa"
    fa.write_text(">c1\nACgtaCG\nttA\n>c2\nACGT\n")
    import pysam
    pysam.faidx(str(fa))
    assert re_.lowercase_runs(fa) == [("c1", 2, 4, "."), ("c1", 7, 9, ".")]


def test_masker_interval_parse(tmp_path):
    p = tmp_path / "wm.txt"
    p.write_text(">c1 some description\n0 - 9\n20 - 20\n>c2\n5 - 6\n")
    assert re_.parse_masker_intervals(p) == [("c1", 0, 10, "."), ("c1", 20, 21, "."), ("c2", 5, 7, ".")]


def test_merge_and_masked_bases():
    m = re_.merge([("c", 0, 10, "."), ("c", 5, 20, "."), ("c", 30, 40, ".")])
    assert m == {"c": [(0, 20), (30, 40)]}
    assert re_.masked_bases(m, "c", 15, 35) == 10 and re_.masked_bases(m, "x", 0, 5) == 0


def test_bed_roundtrip(tmp_path):
    ivs = [("c", 0, 10, "SINE/Alu")]
    re_.write_bed(ivs, tmp_path / "r.bed", "R1:rmsk")
    assert re_.read_bed(tmp_path / "r.bed") == ivs
```

- [ ] **Step 2: Run the tests and verify they fail**

Run: `python3 -m pytest bench/test_repeat_evidence.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'repeat_evidence'`

- [ ] **Step 3: Write the implementation**

```python
#!/usr/bin/env python3
"""Evidence-agnostic families (spec 2026-09-14 §1a): repeat evidence R1-R4 as one BED format
(chrom, start, end, class, source); 0-based half-open. R1 RepeatMasker .out (class), R2 lowercase runs of a soft-masked
assembly, R3 WindowMasker + DustMasker intervals, R4 meryl high-copy runs (dup_evidence.py)."""
import bisect
import collections
import gzip
import re

import pysam


def _open(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_rmsk(path, contigs=None):
    out = []
    for line in _open(path):
        f = line.split()
        if len(f) < 11 or not f[0].isdigit():
            continue
        if contigs is not None and f[4] not in contigs:
            continue
        out.append((f[4], int(f[5]) - 1, int(f[6]), f[10]))
    return out


def lowercase_runs(fasta, contigs=None):
    g = pysam.FastaFile(str(fasta))
    out = []
    for c in g.references:
        if contigs is not None and c not in contigs:
            continue
        for m in re.finditer(r"[a-z]+", g.fetch(c)):
            out.append((c, m.start(), m.end(), "."))
    return out


def parse_masker_intervals(path):
    out, chrom = [], None
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            chrom = line[1:].split()[0]
            continue
        a, b = line.split(" - ")
        out.append((chrom, int(a), int(b) + 1, "."))
    return out


def merge(ivs):
    by = collections.defaultdict(list)
    for c, s, e, _ in ivs:
        by[c].append((s, e))
    out = {}
    for c, v in by.items():
        m = []
        for s, e in sorted(v):
            if m and s <= m[-1][1]:
                m[-1] = (m[-1][0], max(m[-1][1], e))
            else:
                m.append((s, e))
        out[c] = m
    return out


def masked_bases(merged, chrom, s, e):
    v = merged.get(chrom)
    if not v or e <= s:
        return 0
    lo, hi = 0, len(v)
    while lo < hi:  # first interval ending after s
        mid = (lo + hi) // 2
        if v[mid][1] <= s:
            lo = mid + 1
        else:
            hi = mid
    tot = 0
    for a, b in v[lo:]:
        if a >= e:
            break
        tot += max(0, min(b, e) - max(a, s))
    return tot


def write_bed(ivs, path, source):
    with open(path, "w") as fh:
        for c, s, e, cls in ivs:
            fh.write(f"{c}\t{s}\t{e}\t{cls}\t{source}\n")


def read_bed(path):
    return [(f[0], int(f[1]), int(f[2]), f[3]) for f in (l.rstrip("\n").split("\t") for l in open(path)) if len(f) >= 4]
```

- [ ] **Step 4: Run the tests and verify they pass**

Run: `python3 -m pytest bench/test_repeat_evidence.py -q`
Expected: `5 passed`

- [ ] **Step 5: Commit**

```bash
git add bench/repeat_evidence.py bench/test_repeat_evidence.py
git commit -m "Evidence sources: repeat intervals R1-R3 (RepeatMasker, soft-mask, WindowMasker/DustMasker)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 3: meryl k-mer database, C_max valley, D4 low-copy and R4 high-copy runs

**Files:**
- Modify: `bench/dup_evidence.py` (append functions + CLI subcommand `meryl`)
- Test: `bench/test_dup_evidence.py` (append)

**Interfaces:**
- Consumes: `repeat_evidence.write_bed`.
- Produces:
  - `valley(hist: dict[int,int], max_count=100000) -> int | None`;
  - `cmd_meryl(genome, outdir, threads)`, which writes `<outdir>/kmers.meryl`, `hist.tsv`, `cmax.txt` (the integer or `NA`),
    `low_copy.bed` (D4 runs) and `high_copy.bed` (R4 runs) when C_max exists.
- `valley` rule: over counts c >= 2 with frequency f(c), find the smallest c whose f(c) is a strict local minimum AND some
  c' > c has f(c') > f(c) (a second mode exists). Otherwise `None`.

- [ ] **Step 1: Write the failing tests**

```python
# append to bench/test_dup_evidence.py
def test_valley_bimodal_and_monotone():
    hist = {1: 1000, 2: 400, 3: 120, 4: 30, 5: 12, 6: 20, 7: 45, 8: 60, 9: 40, 10: 10}
    assert de.valley(hist) == 5
    assert de.valley({1: 1000, 2: 500, 3: 250, 4: 100, 5: 50}) is None
    assert de.valley({1: 5}) is None
```

- [ ] **Step 2: Run and verify it fails**

Run: `python3 -m pytest bench/test_dup_evidence.py -q -k valley`
Expected: FAIL with `AttributeError: module 'dup_evidence' has no attribute 'valley'`

- [ ] **Step 3: Implement**

```python
# append to bench/dup_evidence.py
import os
import subprocess

MERYL = "/home/juanfra/miniforge3/envs/phasing_eval/bin/meryl"
MERYL_LOOKUP = "/home/juanfra/miniforge3/envs/phasing_eval/bin/meryl-lookup"


def valley(hist, max_count=100000):
    cs = sorted(c for c in hist if 2 <= c <= max_count)
    for i in range(1, len(cs) - 1):
        c = cs[i]
        if hist[c] < hist[cs[i - 1]] and hist[c] < hist[cs[i + 1]] and any(hist[d] > hist[c] for d in cs[i + 1:]):
            return c
    return None


def _run(cmd, out=None):
    with (open(out, "w") if out else open(os.devnull, "w")) as fh:
        subprocess.run(cmd, stdout=fh, stderr=subprocess.DEVNULL, check=True)


def cmd_meryl(genome, outdir, threads=4):
    import repeat_evidence as rep
    os.makedirs(outdir, exist_ok=True)
    db = f"{outdir}/kmers.meryl"
    if not os.path.exists(db):
        _run([MERYL, "count", "k=31", f"threads={threads}", "memory=12", str(genome), "output", db])
    hist_path = f"{outdir}/hist.tsv"
    if not os.path.exists(hist_path):
        _run([MERYL, "histogram", db], hist_path)
    hist = {int(a): int(b) for a, b in (l.split()[:2] for l in open(hist_path) if l.strip() and l.split()[0].isdigit())}
    c = valley(hist)
    open(f"{outdir}/cmax.txt", "w").write(f"{c if c is not None else 'NA'}\n")
    if c is None:
        print(f"[meryl] no histogram valley: D4 and R4 not available for {genome}")
        return None
    for name, ops in (("low", ["less-than", str(c + 1), "[", "greater-than", "1", db, "]"]), ("high", ["greater-than", str(c), db])):
        sub = f"{outdir}/{name}.meryl"
        if not os.path.exists(sub):
            _run([MERYL] + ops + ["output", sub])
        bed = f"{outdir}/{name}_copy.runs.bed"
        if not os.path.exists(bed):
            _run([MERYL_LOOKUP, "-bed-runs", "-sequence", str(genome), "-mers", sub, "-output", bed])
    for name, src in (("low", "D4:meryl"), ("high", "R4:meryl")):
        ivs = [(f[0], int(f[1]), int(f[2]), ".") for f in (l.split("\t") for l in open(f"{outdir}/{name}_copy.runs.bed")) if len(f) >= 3]
        rep.write_bed(ivs, f"{outdir}/{name}_copy.bed", src)
    print(f"[meryl] C_max = {c}")
    return c
```

- [ ] **Step 4: Run the tests and verify they pass**

Run: `python3 -m pytest bench/test_dup_evidence.py -q`
Expected: `7 passed`

- [ ] **Step 5: Smoke-test meryl on the BISER toy genome**

Run:
```bash
S=/tmp/claude-1000/-mnt-c-Users-jfris-Desktop/931c208e-8acb-4dd2-aacb-cf92d5ad051f/scratchpad/biser_probe
python3 -c "import sys; sys.path.insert(0,'bench'); import dup_evidence as de; print(de.cmd_meryl('$S/toy.fa', '$S/meryl_toy', 2))"
head -3 $S/meryl_toy/hist.tsv; cat $S/meryl_toy/cmax.txt
```
Then exercise the filter and lookup syntax, which the toy's `NA` skips:
```bash
M=/home/juanfra/miniforge3/envs/phasing_eval/bin/meryl
$M less-than 3 [ greater-than 1 $S/meryl_toy/kmers.meryl ] output $S/meryl_toy/low_test.meryl
/home/juanfra/miniforge3/envs/phasing_eval/bin/meryl-lookup -bed-runs -sequence $S/toy.fa -mers $S/meryl_toy/low_test.meryl -output $S/meryl_toy/low_test.bed
head -3 $S/meryl_toy/low_test.bed
```
Expected:
- `hist.tsv` has count-1 and count-2 rows (the 20 kb duplicated segment), and `cmax.txt` is `NA` (the toy has no second
  mode).
- `low_test.bed` has runs on chrA near 100000-120000 and on chrB near 60000-80000.
- No exception.

If `meryl` rejects the bracketed `less-than ... [ greater-than ... ]` expression, replace it with two steps: `greater-than 1
db output ge2.meryl`, then `less-than C+1 ge2.meryl output low.meryl`. Rerun Step 5.

- [ ] **Step 6: Commit**

```bash
git add bench/dup_evidence.py bench/test_dup_evidence.py
git commit -m "Evidence sources: meryl k-mer database, histogram valley C_max, D4 low-copy / R4 high-copy runs

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 4: D2 BISER runner, D3 self-alignment runner, and the CLI

**Files:**
- Modify: `bench/dup_evidence.py` (append runners + `main()`)
- Test: `bench/test_dup_evidence.py` (append an integration test on the toy genome)

**Interfaces:**
- Consumes: `repeat_evidence.read_bed`, `merge`, `masked_bases`.
- Produces the CLI:
  - `dup_evidence.py biser --genome FA --outdir DIR [--threads 4]`: resumable; writes `DIR/biser.bed` and `DIR/pairs.D2.tsv`.
  - `dup_evidence.py selfaln --genome FA --outdir DIR --repeats BED [--chunk-bp 2000000] [--budget 540]`: aligns missing
    chunks until the budget is spent. When all are done it writes `DIR/pairs.D3.tsv` (records with >= 1000 non-repeat aligned
    bases).
  - `dup_evidence.py d1 --sedef BED --fmt gorilla --contigs c1,c2 --out DIR/pairs.D1.tsv`.
  - `dup_evidence.py meryl --genome FA --outdir DIR`.
  - `dup_evidence.py to-sedef --pairs TSV --out BED`: gorilla-format rows for `dna_sd_atoms.py ... gorilla ... cigar`.
- Produces the helper `nonrepeat_aligned(pair, merged) -> int` (CIGAR walk; min over sides per M block).

- [ ] **Step 1: Write the failing tests**

```python
# append to bench/test_dup_evidence.py
import subprocess

TOY = "/tmp/claude-1000/-mnt-c-Users-jfris-Desktop/931c208e-8acb-4dd2-aacb-cf92d5ad051f/scratchpad/biser_probe/toy.fa"


def test_nonrepeat_aligned_minus_strand():
    import repeat_evidence as rep
    merged = rep.merge([("B", 0, 100, ".")])  # side B's first 100 bp are repeat
    # side B '-' : the first CIGAR block on B maps to B's END; B = 0..1000, block 1000M covers all of B
    pair = ("A", 0, 1000, "B", 0, 1000, "+", "-", 0.95, "1000M", "D3:selfaln")
    assert de.nonrepeat_aligned(pair, merged) == 900
    pair2 = ("A", 0, 1000, "B", 0, 950, "+", "-", 0.95, "100M50D850M", "D3:selfaln")
    # block1 A 0-100 / B 850-950 (no repeat) = 100; block2 A 150-1000 / B 0-850: B repeat 0-100 -> 750
    assert de.nonrepeat_aligned(pair2, merged) == 850


def test_toy_biser_and_selfaln_find_the_duplication(tmp_path):
    out = tmp_path / "ev"
    r = subprocess.run([sys.executable, "bench/dup_evidence.py", "biser", "--genome", TOY, "--outdir", str(out), "--threads", "2"],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    p2 = de.read_pairs(out / "pairs.D2.tsv")
    assert any({x[0], x[3]} == {"chrA", "chrB"} and x[2] - x[1] > 15000 for x in p2)
    (tmp_path / "empty.bed").write_text("")
    for _ in range(3):
        r = subprocess.run([sys.executable, "bench/dup_evidence.py", "selfaln", "--genome", TOY, "--outdir", str(out),
                            "--repeats", str(tmp_path / "empty.bed"), "--chunk-bp", "100000", "--budget", "120"],
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
    p3 = de.read_pairs(out / "pairs.D3.tsv")
    assert any({x[0], x[3]} == {"chrA", "chrB"} and x[2] - x[1] > 15000 and x[8] > 0.9 for x in p3)
```

- [ ] **Step 2: Run and verify they fail**

Run: `python3 -m pytest bench/test_dup_evidence.py -q -k "nonrepeat or toy"`
Expected: FAIL (`nonrepeat_aligned` missing; CLI has no subcommands)

- [ ] **Step 3: Implement**

```python
# append to bench/dup_evidence.py
import argparse
import glob
import sys
import time

BISER = "/home/juanfra/miniforge3/envs/biser/bin/biser"


def run_budget(cmd, budget, stdout=None):
    """Run cmd in its own process group; on budget expiry kill the WHOLE group (BISER workers, Liftoff's minimap2) so no
    orphan survives (WSL crash rule). Returns True if it finished."""
    import signal
    proc = subprocess.Popen(cmd, stdout=stdout or subprocess.DEVNULL, stderr=subprocess.DEVNULL, start_new_session=True)
    try:
        rc = proc.wait(timeout=budget)
    except subprocess.TimeoutExpired:
        os.killpg(proc.pid, signal.SIGKILL)
        proc.wait()
        return False
    if rc != 0:
        raise subprocess.CalledProcessError(rc, cmd)
    return True


def nonrepeat_aligned(p, merged):
    import repeat_evidence as rep
    ca, a1, a2, cb, b1, b2, sa, sb, ident, cig, _ = p
    oa, ob, tot = 0, 0, 0
    for n, o in re.findall(r"(\d+)([MID])", cig):
        n = int(n)
        if o == "M":
            ga = (a1 + oa, a1 + oa + n)
            gb = (b2 - ob - n, b2 - ob) if sb == "-" else (b1 + ob, b1 + ob + n)
            tot += min(n - rep.masked_bases(merged, ca, *ga), n - rep.masked_bases(merged, cb, *gb))
            oa += n
            ob += n
        elif o == "D":
            oa += n
        else:
            ob += n
    return tot


def cmd_biser(a):
    import pysam
    os.makedirs(a.outdir, exist_ok=True)
    bed, tmp = f"{a.outdir}/biser.bed", f"{a.outdir}/biser_tmp"
    if not os.path.exists(bed):
        cmd = [BISER, "-t", str(a.threads), "-o", bed, "--keep-temp", "-T", tmp]
        if os.path.isdir(tmp):
            cmd += ["--resume", tmp]
        if not os.path.exists(str(a.genome) + ".fai"):
            pysam.faidx(str(a.genome))
        if not run_budget(cmd + [str(a.genome)], a.budget):
            print("[biser] budget spent; rerun the same command to resume")
            return
    pairs = [p for p in (from_biser(l) for l in open(bed)) if p]
    write_pairs(pairs, f"{a.outdir}/pairs.D2.tsv")
    print(f"[biser] {len(pairs)} pairs -> {a.outdir}/pairs.D2.tsv")


def cmd_selfaln(a):
    import pysam
    import repeat_evidence as rep
    d = f"{a.outdir}/selfaln"
    os.makedirs(d, exist_ok=True)
    g = pysam.FastaFile(str(a.genome))
    mmi = f"{d}/genome.mmi"
    if not os.path.exists(mmi):
        subprocess.run(["minimap2", "-x", "asm20", "-t", str(a.threads), "-d", mmi, str(a.genome)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    chunks = []
    for c, L in zip(g.references, g.lengths):
        for off in range(0, L, a.chunk_bp):
            chunks.append((c, off, min(L, off + a.chunk_bp)))
    t0 = time.time()
    done = 0
    for i, (c, s, e) in enumerate(chunks):
        paf = f"{d}/c{i:05d}.paf"
        if os.path.exists(paf):
            done += 1
            continue
        if time.time() - t0 > a.budget:
            break
        q = f"{d}/c{i:05d}.fa"
        open(q, "w").write(f">{c}@{s}\n{g.fetch(c, s, e).upper()}\n")
        with open(paf + ".tmp", "w") as fh:
            ok = run_budget(["minimap2", "-x", "asm20", "-c", "-N", "50", "-p", "0.1", "-t", str(a.threads), mmi, q],
                            max(30, a.budget - (time.time() - t0)), stdout=fh)
        if not ok:
            os.remove(paf + ".tmp")
            print(f"[selfaln] chunk {i} did not finish inside the budget; rerun (a chunk that never fits needs -f 0.001)")
            break
        os.replace(paf + ".tmp", paf)
        os.remove(q)
        done += 1
    print(f"[selfaln] {done}/{len(chunks)} chunks aligned")
    if done < len(chunks):
        return
    merged = rep.merge(rep.read_bed(a.repeats))
    pairs = []
    for paf in sorted(glob.glob(f"{d}/c*.paf")):
        for line in open(paf):
            p = from_selfpaf(line)
            if p and nonrepeat_aligned(p, merged) >= MIN_SD_BP:
                pairs.append(p)
    write_pairs(pairs, f"{a.outdir}/pairs.D3.tsv")
    print(f"[selfaln] {len(pairs)} pairs -> {a.outdir}/pairs.D3.tsv")


def cmd_d1(a):
    contigs = set(a.contigs.split(","))
    pairs = [p for p in (from_sedef(l, a.fmt) for l in open(a.sedef) if not l.startswith("#")) if p and p[0] in contigs and p[3] in contigs]
    write_pairs(pairs, a.out)
    print(f"[d1] {len(pairs)} pairs -> {a.out}")


def cmd_to_sedef(a):
    with open(a.out, "w") as fh:
        for p in read_pairs(a.pairs):
            if p[9]:
                fh.write(to_sedef_gorilla(p) + "\n")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("biser")
    p.add_argument("--genome", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--budget", type=int, default=560)
    p = sub.add_parser("selfaln")
    p.add_argument("--genome", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--repeats", required=True)
    p.add_argument("--chunk-bp", type=int, default=2_000_000)
    p.add_argument("--budget", type=int, default=540)
    p.add_argument("--threads", type=int, default=4)
    p = sub.add_parser("d1")
    p.add_argument("--sedef", required=True)
    p.add_argument("--fmt", required=True)
    p.add_argument("--contigs", required=True)
    p.add_argument("--out", required=True)
    p = sub.add_parser("meryl")
    p.add_argument("--genome", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--threads", type=int, default=4)
    p = sub.add_parser("to-sedef")
    p.add_argument("--pairs", required=True)
    p.add_argument("--out", required=True)
    a = ap.parse_args()
    if a.cmd == "meryl":
        cmd_meryl(a.genome, a.outdir, a.threads)
    else:
        {"biser": cmd_biser, "selfaln": cmd_selfaln, "d1": cmd_d1, "to-sedef": cmd_to_sedef}[a.cmd](a)


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run the tests and verify they pass**

Run: `python3 -m pytest bench/test_dup_evidence.py -q`
Expected: `9 passed`.
- If BISER rejects `--keep-temp -T` on a fresh run, drop `--keep-temp` for runs without an existing temp directory and
  rerun.
- If the toy minus-strand test fails, recheck the side-B coordinate mapping against `dna_sd_atoms.py` K0 semantics before
  changing the test.

- [ ] **Step 5: Validate the adapters against `dna_sd_atoms.py` on real gorilla SEDEF rows**

Run:
```bash
L=/mnt/linuxdisk/home/juanfraitu/o1_falsemerge/lit/evid/GGO; mkdir -p $L
python3 bench/dup_evidence.py d1 --sedef /mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_sedef_final.bed --fmt gorilla --contigs NC_073241.2,NC_073242.2,NC_073244.2 --out $L/check.D1.tsv
python3 bench/dup_evidence.py to-sedef --pairs $L/check.D1.tsv --out $L/check.D1.sedef.bed
timeout 590 python3 bench/dna_sd_atoms.py $L/check.D1.sedef.bed gorilla NC_073241.2,NC_073242.2,NC_073244.2 $L/check_rt cigar
awk -F'\t' '$1=="NC_073241.2"||$1=="NC_073242.2"||$1=="NC_073244.2"' /mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_sedef_final.bed | awk -F'\t' '$4=="NC_073241.2"||$4=="NC_073242.2"||$4=="NC_073244.2"' > $L/check.orig.bed
timeout 590 python3 bench/dna_sd_atoms.py $L/check.orig.bed gorilla NC_073241.2,NC_073242.2,NC_073244.2 $L/check_orig cigar
cmp $L/check_rt.nodes.tsv $L/check_orig.nodes.tsv && echo ATOMS_IDENTICAL
python3 -c "
import csv
def E(p): return {(r[0],r[1]):(round(float(r[2]),3),round(float(r[3]),3)) for r in (l.split('\t') for l in open(p)) if r[0].isdigit()}
a,b=E('$L/check_rt.edges.tsv'),E('$L/check_orig.edges.tsv'); print('edges', len(a), len(b), 'identical' if a==b else 'DIFFER')"
```
Expected: `ATOMS_IDENTICAL` and `edges N N identical`.
- If the edges differ only in identity at the 3rd decimal, the cause is the synthetic matches/mismatches rounding. Record it
  in the report, and do not change `dna_sd_atoms.py`.

- [ ] **Step 6: Commit**

```bash
git add bench/dup_evidence.py bench/test_dup_evidence.py
git commit -m "Evidence sources: BISER (D2) and minimap2 self-alignment (D3) runners, D1 import, SEDEF-format export validated on dna_sd_atoms

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 5: Ape substrates from human-lifted orthologous chromosomes

**Files:**
- Create: `bench/ape_substrate.py`
- Test: `bench/test_ape_substrate.py`

**Interfaces:**
- Consumes: `annotation_nodes.py` node tables (Plan 2 uses them).
- Produces:
  - `gff_batches(lines, max_genes) -> list[list[str]]`: every line stays with its top-level gene; a batch closes at
    `max_genes` genes.
  - `select_substrate(counts: dict[(human_chrom, target_chrom), int], min_frac=0.10) -> list[target_chrom]`.
  - CLI `ape_substrate.py lift --species NAME --target FA --out DIR [--max-genes 1500]`: one batch per call until all are
    lifted; human chromosomes chr7, chr15, chr16, chr17.
  - CLI `ape_substrate.py select --out DIR`: writes `DIR/substrate.txt` and the counts matrix `DIR/lift_matrix.tsv`.
  - CLI `ape_substrate.py subset --target FA --out DIR [--native-gff GFF]`: writes `DIR/substrate.fa`; with `--native-gff`
    also `DIR/native.gff` (only for gorilla/chimpanzee in Plan 1) and `DIR/native_subsample50.txt` (gene names, seed 1, via
    `guided_min.load_genes`).
  - `DIR/lifted.gff`: concatenated Liftoff output with an added `human_source=<chrom>` attribute.
- `min_frac = 0.10` is a scoping choice, not a definition constant. A target chromosome is in the substrate iff it receives
  >= 10% of the lifted genes of at least one of the four human chromosomes. The whole matrix is reported.

- [ ] **Step 1: Write the failing tests**

```python
# bench/test_ape_substrate.py
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ape_substrate as asub  # noqa: E402

GFF = [
    "chr7\tRefSeq\tgene\t100\t900\t.\t+\t.\tID=gene-A;Name=A\n",
    "chr7\tRefSeq\tmRNA\t100\t900\t.\t+\t.\tID=rna-A1;Parent=gene-A\n",
    "chr7\tRefSeq\texon\t100\t200\t.\t+\t.\tID=exon-A1;Parent=rna-A1\n",
    "chr7\tRefSeq\tgene\t1000\t1900\t.\t+\t.\tID=gene-B;Name=B\n",
    "chr7\tRefSeq\texon\t1000\t1100\t.\t+\t.\tID=exon-B1;Parent=gene-B\n",
    "chr7\tRefSeq\tpseudogene\t2000\t2900\t.\t+\t.\tID=gene-C;Name=C\n",
]


def test_gff_batches_keep_children_with_gene():
    b = asub.gff_batches(GFF, max_genes=2)
    assert len(b) == 2
    assert [l.split("\t")[8].split(";")[0] for l in b[0]] == ["ID=gene-A", "ID=rna-A1", "ID=exon-A1", "ID=gene-B", "ID=exon-B1"]
    assert b[1] == [GFF[5]]


def test_select_substrate_handles_translocation():
    counts = {("chr7", "t7"): 900, ("chr7", "t3"): 20, ("chr17", "t17"): 600, ("chr17", "t5"): 300,
              ("chr15", "t15"): 700, ("chr16", "t16"): 800, ("chr16", "t1"): 50}
    assert asub.select_substrate(counts) == ["t15", "t16", "t17", "t5", "t7"]
```

- [ ] **Step 2: Run and verify they fail**

Run: `python3 -m pytest bench/test_ape_substrate.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'ape_substrate'`

- [ ] **Step 3: Implement**

```python
#!/usr/bin/env python3
"""Evidence-agnostic families (spec 2026-09-14 §2): ape substrates = target chromosomes orthologous to human chr7/15/16/17,
found by Liftoff placement of human RefSeq genes (batched foreground calls); substrate FASTA; native GFF subset and the
50% native subsample (development/report species only)."""
import argparse
import collections
import glob
import os
import random
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

HUMAN_FA = "/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa"
HUMAN_GFF = "/mnt/linuxdisk/home/juanfraitu/winloci_data/Reference/chm13v2.0_RefSeq_full.gff.gz"
LIFTOFF = "/home/juanfra/miniforge3/envs/liftoff/bin/liftoff"
HUMAN_CHROMS = ("chr7", "chr15", "chr16", "chr17")
TOP = ("gene", "pseudogene")


def attrs(col):
    return dict(kv.split("=", 1) for kv in col.strip().split(";") if "=" in kv)


def gff_batches(lines, max_genes):
    top_of, order, members = {}, [], collections.defaultdict(list)
    for line in lines:
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9:
            continue
        a = attrs(f[8])
        if f[2] in TOP and "ID" in a:
            top_of[a["ID"]] = a["ID"]
            order.append(a["ID"])
            members[a["ID"]].append(line)
            continue
        p = a.get("Parent", "").split(",")[0]
        top = top_of.get(p)
        if top is None:
            continue
        if "ID" in a:
            top_of[a["ID"]] = top
        members[top].append(line)
    batches, cur, n = [], [], 0
    for gid in order:
        cur.extend(members[gid])
        n += 1
        if n == max_genes:
            batches.append(cur)
            cur, n = [], 0
    if cur:
        batches.append(cur)
    return batches


def select_substrate(counts, min_frac=0.10):
    tot = collections.Counter()
    for (h, t), n in counts.items():
        tot[h] += n
    keep = {t for (h, t), n in counts.items() if tot[h] and n / tot[h] >= min_frac}
    return sorted(keep)


def cmd_lift(a):
    os.makedirs(a.out, exist_ok=True)
    bdir = f"{a.out}/batches"
    os.makedirs(bdir, exist_ok=True)
    if not glob.glob(f"{bdir}/*.gff"):
        for h in HUMAN_CHROMS:
            lines = subprocess.run(["tabix", HUMAN_GFF, h], capture_output=True, text=True, check=True).stdout.splitlines(True)
            for i, b in enumerate(gff_batches(lines, a.max_genes)):
                open(f"{bdir}/{h}.{i:03d}.gff", "w").writelines(b)
    todo = [p for p in sorted(glob.glob(f"{bdir}/*.gff")) if not os.path.exists(p.replace(".gff", ".lifted.gff3"))]
    if not todo:
        print("[lift] all batches lifted")
        return
    import shutil
    from dup_evidence import run_budget
    b = todo[0]
    out = b.replace(".gff", ".lifted.gff3")
    ok = run_budget([LIFTOFF, "-g", b, "-o", out + ".tmp", "-u", b.replace(".gff", ".unmapped.txt"), "-dir",
                     f"{a.out}/liftoff_tmp", "-p", str(a.threads), "-m", shutil.which("minimap2"), a.target, HUMAN_FA], a.budget)
    if not ok:
        if os.path.exists(out + ".tmp"):
            os.remove(out + ".tmp")
        print(f"[lift] {os.path.basename(b)} did not finish in {a.budget} s; regenerate batches with a smaller --max-genes")
        return
    os.replace(out + ".tmp", out)
    print(f"[lift] {os.path.basename(b)} done; {len(todo) - 1} batches left")


def cmd_select(a):
    counts = collections.Counter()
    with open(f"{a.out}/lifted.gff", "w") as fo:
        for p in sorted(glob.glob(f"{a.out}/batches/*.lifted.gff3")):
            h = os.path.basename(p).split(".")[0]
            for line in open(p):
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                if f[2] in TOP:
                    counts[(h, f[0])] += 1
                fo.write("\t".join(f[:8] + [f[8] + f";human_source={h}"]) + "\n")
    sub = select_substrate(counts)
    with open(f"{a.out}/lift_matrix.tsv", "w") as fh:
        fh.write("human_chrom\ttarget_chrom\tgenes\n")
        for (h, t), n in sorted(counts.items()):
            fh.write(f"{h}\t{t}\t{n}\n")
    open(f"{a.out}/substrate.txt", "w").write("\n".join(sub) + "\n")
    print(f"[select] substrate: {','.join(sub)}")


def cmd_subset(a):
    sub = open(f"{a.out}/substrate.txt").read().split()
    fa = f"{a.out}/substrate.fa"
    if not os.path.exists(fa):
        with open(fa, "w") as fh:
            subprocess.run(["samtools", "faidx", a.target] + sub, stdout=fh, check=True)
        subprocess.run(["samtools", "faidx", fa], check=True)
    if a.native_gff:
        import guided_min
        s = set(sub)
        with open(f"{a.out}/native.gff", "w") as fh:
            for line in open(a.native_gff):
                if not line.startswith("#") and line.split("\t", 1)[0] in s:
                    fh.write(line)
        genes, _ = guided_min.load_genes(f"{a.out}/native.gff", s)
        names = sorted(genes)
        keep = sorted(random.Random(1).sample(names, round(0.5 * len(names))))
        open(f"{a.out}/native_subsample50.txt", "w").write("\n".join(keep) + "\n")
        print(f"[subset] native genes {len(names)}; subsample {len(keep)}")
    print(f"[subset] {fa}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("lift")
    p.add_argument("--species", required=True)
    p.add_argument("--target", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--max-genes", type=int, default=1500)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--budget", type=int, default=560)
    p = sub.add_parser("select")
    p.add_argument("--out", required=True)
    p = sub.add_parser("subset")
    p.add_argument("--target", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--native-gff")
    a = ap.parse_args()
    {"lift": cmd_lift, "select": cmd_select, "subset": cmd_subset}[a.cmd](a)


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run the tests and verify they pass**

Run: `python3 -m pytest bench/test_ape_substrate.py -q`
Expected: `2 passed`

- [ ] **Step 5: Pilot one Liftoff batch on gorilla (timing)**

Run:
```bash
L=/mnt/linuxdisk/home/juanfraitu/o1_falsemerge/lit/evid/GGO
time timeout 590 python3 bench/ape_substrate.py lift --budget 560 --species GGO --target /mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta --out $L/sub
ls $L/sub/batches | head
```
Expected: `[lift] chr15.000.gff done; N batches left` within 590 s.
- If it times out, delete the partial `.tmp` file, then rerun with `--max-genes 500`: remove `$L/sub/batches/*.gff` first so
  the batches are regenerated.
- Record the per-batch time in the Task 7 report.

- [ ] **Step 6: Lift all batches, select, subset (gorilla, chimpanzee, orangutan)**

Repeat the `lift` command (one foreground call per batch) until it prints `all batches lifted`. Then run:
```bash
python3 bench/ape_substrate.py select --out $L/sub
python3 bench/ape_substrate.py subset --target /mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta --out $L/sub --native-gff /mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_genomic.gff
```

Do the same for chimpanzee:
- `L=.../evid/PTR`
- `--target .../winloci_data/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna`
- `--native-gff .../PTR_genomic.gff`

Do the same for orangutan, **without `--native-gff`**:
- `L=.../evid/PPY`
- `--target .../winloci_data/GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.fna`

Expected: each species prints a substrate of 4–6 chromosomes. Gorilla includes the chromosome(s) carrying human chr17's
t(5;17) segment.

- [ ] **Step 7: Commit (code only; outputs stay on /mnt/linuxdisk)**

```bash
git add bench/ape_substrate.py bench/test_ape_substrate.py
git commit -m "Evidence sources: ape substrates from Liftoff placement of human chr7/15/16/17 genes (batched), native subsets

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 6: Run every available source on the gorilla and chimpanzee substrates

**Files:**
- Modify: none (runs only). Outputs go under `/mnt/linuxdisk/home/juanfraitu/o1_falsemerge/lit/evid/{GGO,PTR}/`.

**Interfaces:**
- Consumes: Tasks 1–5 CLIs.
- Produces, per species directory `L`:
  - `R1.bed` (gorilla only), `R2.bed`, `R3.bed`;
  - `meryl/` (`cmax.txt`, `low_copy.bed`, `high_copy.bed` when C_max exists);
  - `pairs.D1.tsv` (gorilla only), `pairs.D2.tsv`, `pairs.D3.tsv`.

- [ ] **Step 1: Repeat sources (gorilla)**

```bash
L=/mnt/linuxdisk/home/juanfraitu/o1_falsemerge/lit/evid/GGO; G=$L/sub/substrate.fa; C=$(paste -sd, $L/sub/substrate.txt)
python3 -c "
import sys; sys.path.insert(0,'bench'); import repeat_evidence as r
s=set('$C'.split(','))
r.write_bed(r.parse_rmsk('/mnt/linuxdisk/home/juanfraitu/winloci_data/rmsk/GCF_029281585.2.repeatMasker.out.gz', s), '$L/R1.bed', 'R1:rmsk')
r.write_bed(r.lowercase_runs('$G', s), '$L/R2.bed', 'R2:softmask')"
timeout 590 /home/juanfra/miniforge3/envs/blast/bin/windowmasker -mk_counts -in $G -out $L/wm.counts
timeout 590 /home/juanfra/miniforge3/envs/blast/bin/windowmasker -ustat $L/wm.counts -in $G -outfmt interval -out $L/wm.txt
timeout 590 /home/juanfra/miniforge3/envs/blast/bin/dustmasker -in $G -outfmt interval -out $L/dust.txt
python3 -c "
import sys; sys.path.insert(0,'bench'); import repeat_evidence as r
r.write_bed(r.parse_masker_intervals('$L/wm.txt') + r.parse_masker_intervals('$L/dust.txt'), '$L/R3.bed', 'R3:windowmasker+dust')"
wc -l $L/R1.bed $L/R2.bed $L/R3.bed
```
Expected: three non-empty BED files. If a WindowMasker step exceeds 590 s, run it per chromosome
(`samtools faidx $G <chrom>`), then concatenate the interval files.

- [ ] **Step 2: meryl (D4 + R4)**

```bash
timeout 590 python3 bench/dup_evidence.py meryl --genome $G --outdir $L/meryl
cat $L/meryl/cmax.txt
```
Expected: an integer C_max, or `NA`. If `meryl count` exceeds 590 s, rerun the same command: `meryl count` has no resume
(the call restarts it), so lower to `threads=4 memory=8` and count per chromosome, then run `meryl union-sum`.

- [ ] **Step 3: Duplication pairs**

```bash
python3 bench/dup_evidence.py d1 --sedef /mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_sedef_final.bed --fmt gorilla --contigs $C --out $L/pairs.D1.tsv
timeout 590 python3 bench/dup_evidence.py biser --genome $G --outdir $L --threads 4 --budget 560   # repeat until pairs.D2.tsv exists; after any budget stop run `ps -eo pid,args | grep -i biser` and kill leftovers by PID
timeout 590 python3 bench/dup_evidence.py selfaln --genome $G --outdir $L --repeats $L/R2.bed --budget 540   # repeat until pairs.D3.tsv exists
wc -l $L/pairs.D*.tsv
```
Expected: three pair tables.
- If one self-alignment chunk alone exceeds the budget (satellite arrays), add `-f 0.001` to the minimap2 command in
  `cmd_selfaln` for that run.
- Record that change in the Task 7 report.

- [ ] **Step 4: Chimpanzee (report species)**

Repeat Steps 1–3 with `L=.../evid/PTR`, skipping R1 and D1 (they do not exist for chimpanzee).

- [ ] **Step 5: Orangutan duplication and repeat sources (hold-out; no native annotation read)**

Repeat Steps 1–3 with `L=.../evid/PPY`, skipping R1 and D1. These are genome-only computations and are allowed before AO;
Plan 2 consumes them.

---

### Task 7: Source-cost comparison on gorilla (report, then ledger)

**Files:**
- Create: `bench/evidence_compare.py`
- Test: `bench/test_evidence_compare.py`
- Modify: `docs/o1_ledger.md` (append section §6kp)

**Interfaces:**
- Consumes: pair tables, repeat BEDs, meryl outputs.
- Produces:
  - `pair_match(ref_pairs, test_pairs, min_ro=0.5) -> (recall, precision)`: a pair matches if both sides reciprocally
    overlap >= 0.5 in either orientation.
  - `bed_overlap(ref_merged, test_merged) -> (ref_covered_frac, test_in_ref_frac, jaccard)`.
  - CLI `evidence_compare.py --species-dir DIR` prints the report table.

- [ ] **Step 1: Write the failing tests**

```python
# bench/test_evidence_compare.py
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evidence_compare as ec  # noqa: E402


def P(ca, a1, a2, cb, b1, b2):
    return (ca, a1, a2, cb, b1, b2, "+", "+", 0.95, "", "x")


def test_pair_match_reciprocal_and_orientation():
    ref = [P("c1", 0, 1000, "c2", 0, 1000), P("c1", 5000, 7000, "c3", 0, 2000)]
    test = [P("c2", 50, 1050, "c1", 0, 900), P("c1", 5000, 5500, "c3", 0, 500), P("c4", 0, 10, "c5", 0, 10)]
    r, p = ec.pair_match(ref, test)
    assert (r, p) == (0.5, 1 / 3)


def test_bed_overlap():
    import repeat_evidence as rep
    a = rep.merge([("c", 0, 100, ".")])
    b = rep.merge([("c", 50, 150, ".")])
    assert ec.bed_overlap(a, b) == (0.5, 0.5, 50 / 150)
```

- [ ] **Step 2: Run and verify they fail**

Run: `python3 -m pytest bench/test_evidence_compare.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'evidence_compare'`

- [ ] **Step 3: Implement**

```python
#!/usr/bin/env python3
"""Evidence-agnostic families (spec 2026-09-14 §2, evidence-source comparison): D2/D3 pairs vs D1, R2/R3/R4 masks vs R1,
C_max check (R4 bases inside R1), D4 coverage of D1 SD bases.
usage: evidence_compare.py --species-dir DIR"""
import argparse
import collections
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dup_evidence as de  # noqa: E402
import repeat_evidence as rep  # noqa: E402


def _ro(s1, e1, s2, e2):
    ov = min(e1, e2) - max(s1, s2)
    return ov > 0 and ov >= 0.5 * (e1 - s1) and ov >= 0.5 * (e2 - s2)


def _same(p, q):
    a = (p[0], p[1], p[2]), (p[3], p[4], p[5])
    b = (q[0], q[1], q[2]), (q[3], q[4], q[5])
    for x, y in ((b[0], b[1]), (b[1], b[0])):
        if a[0][0] == x[0] and a[1][0] == y[0] and _ro(a[0][1], a[0][2], x[1], x[2]) and _ro(a[1][1], a[1][2], y[1], y[2]):
            return True
    return False


def pair_match(ref, test, min_ro=0.5):
    idx = collections.defaultdict(list)
    for i, q in enumerate(test):
        idx[(q[0], q[1] // 100000)].append(i)
        idx[(q[3], q[4] // 100000)].append(i)

    def cands(p):
        out = set()
        for c, s, e in ((p[0], p[1], p[2]), (p[3], p[4], p[5])):
            for b in range(s // 100000 - 1, e // 100000 + 2):
                out.update(idx.get((c, b), ()))
        return out
    hit_ref, hit_test = 0, set()
    for p in ref:
        m = [i for i in cands(p) if _same(p, test[i])]
        if m:
            hit_ref += 1
            hit_test.update(m)
    return hit_ref / max(1, len(ref)), len(hit_test) / max(1, len(test))


def bed_overlap(ref, test):
    def total(m):
        return sum(e - s for v in m.values() for s, e in v)
    inter = sum(rep.masked_bases(ref, c, s, e) for c, v in test.items() for s, e in v)
    rt, tt = total(ref), total(test)
    return inter / max(1, rt), inter / max(1, tt), inter / max(1, rt + tt - inter)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species-dir", required=True)
    a = ap.parse_args()
    L = a.species_dir
    print("== duplication sources vs D1 (SEDEF): pair recall / precision (reciprocal overlap >= 0.5 both sides)")
    d1 = de.read_pairs(f"{L}/pairs.D1.tsv") if os.path.exists(f"{L}/pairs.D1.tsv") else None
    for k in ("D2", "D3"):
        p = f"{L}/pairs.{k}.tsv"
        if d1 is None or not os.path.exists(p):
            print(f"  {k}: not available")
            continue
        t = de.read_pairs(p)
        r, pr = pair_match(d1, t)
        print(f"  {k}: pairs {len(t)} vs D1 {len(d1)}; recall {r:.3f}; precision {pr:.3f}")
    print("== repeat sources vs R1 (RepeatMasker): R1 bases covered / source bases inside R1 / Jaccard")
    r1 = rep.merge(rep.read_bed(f"{L}/R1.bed")) if os.path.exists(f"{L}/R1.bed") else None
    for k, p in (("R2", f"{L}/R2.bed"), ("R3", f"{L}/R3.bed"), ("R4", f"{L}/meryl/high_copy.bed")):
        if r1 is None or not os.path.exists(p):
            print(f"  {k}: not available")
            continue
        print(f"  {k}: " + " / ".join(f"{x:.3f}" for x in bed_overlap(r1, rep.merge(rep.read_bed(p)))))
    cm = open(f"{L}/meryl/cmax.txt").read().strip() if os.path.exists(f"{L}/meryl/cmax.txt") else "not run"
    print(f"== C_max: {cm}")
    if d1 is not None and os.path.exists(f"{L}/meryl/low_copy.bed"):
        sd = rep.merge([(p[0], p[1], p[2], ".") for p in d1] + [(p[3], p[4], p[5], ".") for p in d1])
        low = rep.merge(rep.read_bed(f"{L}/meryl/low_copy.bed"))
        cov, inside, j = bed_overlap(sd, low)
        print(f"== D4 low-copy runs vs D1 SD bases: SD bases covered {cov:.3f}; low-copy bases inside SDs {inside:.3f}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run the tests and verify they pass**

Run: `python3 -m pytest bench/test_evidence_compare.py -q`
Expected: `2 passed`

- [ ] **Step 5: Run the gorilla report**

Run: `timeout 590 python3 bench/evidence_compare.py --species-dir /mnt/linuxdisk/home/juanfraitu/o1_falsemerge/lit/evid/GGO | tee /mnt/linuxdisk/home/juanfraitu/o1_falsemerge/lit/evid/GGO/compare.txt`
Expected: the four report blocks with numbers ("not available" only where a source legitimately failed, e.g. C_max `NA`).

- [ ] **Step 6: Append ledger §6kp and commit**

Append `docs/o1_ledger.md` section `## §6kp — Evidence sources on gorilla: what each missing source costs (Plan 1)` with:
- the substrate chromosomes and lift matrix per species;
- Liftoff batch timing;
- D2/D3 recall/precision vs SEDEF;
- R2/R3/R4 vs RepeatMasker;
- C_max per species;
- D4 coverage of SD bases;
- any runtime deviations (`-f 0.001`, per-chromosome WindowMasker, meryl per chromosome).

State which sources Plan 2 will use as orangutan defaults: **the best-agreeing of D2/D3 and of R2/R3/R4 on gorilla.** Then:

```bash
git add bench/evidence_compare.py bench/test_evidence_compare.py docs/o1_ledger.md
git commit -m "Evidence sources: gorilla source-cost comparison (D2/D3 vs SEDEF, R2-R4 vs RepeatMasker, C_max, D4) — ledger §6kp

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

## After Plan 1

Write Plan 2 (`docs/superpowers/plans/2026-09-1x-evidence-layers-plan2.md`) from the spec sections 1b, 2–4:
- `protein_projection.py`, `sd_segments.py`, `layers_map.py`, `evidence_arms.py`;
- arms T/1/2a/2b on gorilla, and the gorilla evidence-source sweep over layers;
- pre-registration Addendum AO with the orangutan source defaults from §6kp;
- orangutan hold-out runs;
- chimpanzee report.
