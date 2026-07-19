# DNA Variation-Graph Ceiling vs Soto SDs — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A new bench script `bench/soto/soto_dna_vg.py` that builds, per flagship Soto SD family, a base-level variation-graph GFA from the family's member genomic sequences, colored green (RNA-recovered) / red (DNA-only), so the pictures show DNA = 100% identifiability ceiling and RNA = 76.2% of it.

**Architecture:** Pure core (member extraction → abpoa column-MSA → GFA → node colours) with a thin I/O shell (`main`) that builds the flagships and writes `<fam>.gfa` + `<fam>.colours.csv` + `<fam>.legend.tsv` + `index.md`. Per-family, foreground, serial; abpoa (base-level) is the graph builder; nothing whole-genome, no `copy_assign`.

**Tech Stack:** Python 3 (`pyabpoa` for the MSA — already a dependency; stdlib otherwise). pytest 9.0.3.

**⚠ Split interpreters (verified):** `/home/juanfra/.local/bin/pytest` has **pytest but NOT `pyabpoa`**; `/home/juanfra/miniforge3/bin/python` has **`pyabpoa` but NOT pytest**. So: the module imports `pyabpoa` **inside** `abpoa_msa` (module import stays pyabpoa-free); Tasks 1–4 unit tests are hermetic and run via `.local` pytest (Task 4 monkeypatches `abpoa_msa`); the real-abpoa Task 5 smoke is a **standalone script run via miniforge python**, not a pytest test. Do not install anything.

## Global Constraints

- **Visual artifacts only**: emit `<fam>.gfa` + `<fam>.colours.csv` + `<fam>.legend.tsv` (+ `index.md` captions) per flagship; no catalog-wide number required.
- **Honest ceiling framing**: the VG *represents* the given copies (the ceiling); it does not "detect" families. Every caption carries this.
- **Presence is checked, not assumed**: per family, assert every Soto member appears as a `P`-line; a member missing from the graph is logged honestly, never silently counted as present.
- **`pyabpoa` needs `str` inputs** (bytes silently yield all-`N`). The MSA uses base-level POA so SNP-distinct copies remain distinct paths. `vg msga` is an optional `--builder vg` alternative (deprecated); `minigraph` is rejected (SV-level).
- Colours: green `#1e8e3e` (RNA-recovered), red `#d93025` (DNA-only), grey `#9aa0a6` (shared) — consistent with `copy_graph.rs`.
- Per-family, foreground, serial; member seqs from the 11 MB `soto_members.fa` (never the 3 GB genome, except the optional SLC9B1P1 loci); **no `copy_assign`**; outputs under `/home/juanfra/winloci_scratch/soto_vg/`.

---

### Task 1: Member parsing + genomic-sequence extraction

**Files:**
- Create: `bench/soto/soto_dna_vg.py`
- Test: `bench/soto/test_soto_dna_vg.py`

**Interfaces:**
- Produces:
  - `parse_family_members(bed_lines, family_id) -> list[tuple[str,str,int,int]]` — `(gene, chrom, start, end)` for every member of `family_id`, in file order. BED col4 is `GENE|ID_k`; keep rows whose `ID_k == family_id`.
  - `read_fasta(path) -> dict[str,str]` — header (first token after `>`) → sequence.
  - `member_seq(fa, chrom, start, end) -> str` — look up `f"{chrom}:{start+1}-{end}"` in `fa` (the BED `start` is 0-based; `soto_members.fa` headers are 1-based, so `start+1`). Raises `KeyError` if absent.

- [ ] **Step 1: Write the failing tests**

Create `bench/soto/test_soto_dna_vg.py`:

```python
#!/usr/bin/env python3
"""Tests for the DNA variation-graph ceiling demo (bench/soto/soto_dna_vg.py).
Run: PYTHONHASHSEED=0 pytest bench/soto/test_soto_dna_vg.py -v"""
import os
import sys
sys.path.insert(0, "bench/soto")
import soto_dna_vg as V


def test_parse_family_members():
    bed = [
        "chr1\t100\t200\tSRGAP2C|ID_462\t0\t.",
        "chr7\t68507281\t68522154\tPMS2P4|ID_8\t0\t.",
        "chr1\t300\t400\tSRGAP2|ID_462\t0\t.",
    ]
    m = V.parse_family_members(bed, "ID_462")
    assert m == [("SRGAP2C", "chr1", 100, 200), ("SRGAP2", "chr1", 300, 400)]


def test_member_seq_plus_one_offset():
    fa = {"chr1:101-200": "ACGT", "chr7:1-9": "TTTT"}
    # BED start 100 (0-based) -> header start 101 (1-based)
    assert V.member_seq(fa, "chr1", 100, 200) == "ACGT"


def test_read_fasta(tmp_path):
    p = tmp_path / "m.fa"
    p.write_text(">chr1:101-200\nAC\nGT\n>chr7:1-9\nTTTT\n")
    assert V.read_fasta(str(p)) == {"chr1:101-200": "ACGT", "chr7:1-9": "TTTT"}
```

- [ ] **Step 2: Run to verify failure**

Run: `PYTHONHASHSEED=0 /home/juanfra/.local/bin/pytest bench/soto/test_soto_dna_vg.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'soto_dna_vg'`.

- [ ] **Step 3: Create the module with the parsing + extraction helpers**

Create `bench/soto/soto_dna_vg.py`:

```python
#!/usr/bin/env python3
"""DNA variation-graph ceiling vs Soto SDs (visual demo). Per flagship Soto family, build a base-level
variation graph (abpoa MSA -> GFA) from the member GENOMIC sequences and colour each copy green
(RNA-recovered) / red (DNA-only). The graph REPRESENTS all copies (the DNA ceiling); it does not
"detect" families. See docs/superpowers/specs/2026-07-18-soto-dna-vg-ceiling-design.md."""
import os
import sys
from collections import defaultdict

# ---- paths ----
BED = "bench/soto/80_fams.chr.bed"
MEMFA = "/mnt/linuxdisk/home/juanfraitu/winloci_data/soto_members.fa"
DETECT = "bench/soto/soto_member_detection.tsv"
OUT = "/home/juanfra/winloci_scratch/soto_vg"

# ---- colours (match copy_graph.rs) ----
GREEN = "#1e8e3e"   # RNA-recovered
RED = "#d93025"     # DNA-only (K=0 / silent / coverage)
GREY = "#9aa0a6"    # shared / conserved


def parse_family_members(bed_lines, family_id):
    """(gene, chrom, start, end) for every member of family_id, in file order. BED col4 = GENE|ID_k."""
    out = []
    for ln in bed_lines:
        f = ln.rstrip("\n").split("\t")
        if len(f) < 4:
            continue
        name = f[3]
        if "|" not in name:
            continue
        gene, fam = name.rsplit("|", 1)
        if fam == family_id:
            out.append((gene, f[0], int(f[1]), int(f[2])))
    return out


def read_fasta(path):
    seqs, cur, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if cur is not None:
                seqs[cur] = "".join(buf)
            cur = line[1:].strip().split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if cur is not None:
        seqs[cur] = "".join(buf)
    return seqs


def member_seq(fa, chrom, start, end):
    """Genomic sequence for a BED member. soto_members.fa headers are 1-based (start+1)."""
    return fa[f"{chrom}:{start + 1}-{end}"]
```

- [ ] **Step 4: Run to verify pass**

Run: `PYTHONHASHSEED=0 /home/juanfra/.local/bin/pytest bench/soto/test_soto_dna_vg.py -v`
Expected: PASS (3 tests).

- [ ] **Step 5: Commit**

```bash
git add bench/soto/soto_dna_vg.py bench/soto/test_soto_dna_vg.py
git commit -m "feat(soto-vg): member parsing + genomic-sequence extraction"
```

---

### Task 2: abpoa column-MSA → base-level variation-graph GFA

**Files:**
- Modify: `bench/soto/soto_dna_vg.py`
- Test: `bench/soto/test_soto_dna_vg.py` (add tests)

**Interfaces:**
- Consumes: nothing from Task 1.
- Produces:
  - `abpoa_msa(seqs: list[str]) -> list[str]` — run abpoa, return the aligned gapped rows for the members only (the first `len(seqs)` rows of `msa_seq`; `str` inputs — bytes yield all-`N`).
  - `msa_to_gfa(rows: list[str], names: list[str]) -> tuple[str, dict[str, list[str]]]` — pure. `rows` are equal-length uppercase gapped strings (`-` = gap), one per `names`. Returns `(gfa_text, paths)` where `paths[name]` is the ordered node-id list. Deterministic.

- [ ] **Step 1: Write the failing tests**

Append to `bench/soto/test_soto_dna_vg.py`:

```python
def test_msa_to_gfa_snp_bubble():
    # col0 "A" invariant; col1 variant (C/G/C); col2-3 "GT" invariant
    rows = ["ACGT", "AGGT", "ACGT"]
    gfa, paths = V.msa_to_gfa(rows, ["m1", "m2", "m3"])
    assert paths["m1"] == ["1", "2", "4"]   # A , C , GT
    assert paths["m2"] == ["1", "3", "4"]   # A , G , GT
    assert paths["m3"] == ["1", "2", "4"]
    assert "S\t2\tC" in gfa and "S\t3\tG" in gfa
    assert "S\t1\tA" in gfa and "S\t4\tGT" in gfa
    assert "P\tm2\t1+,3+,4+\t*" in gfa
    # link m2 traverses 1->3 and 3->4
    links = {(l.split("\t")[1], l.split("\t")[3]) for l in gfa.splitlines() if l.startswith("L\t")}
    assert ("1", "3") in links and ("3", "4") in links


def test_msa_to_gfa_indel_skips_gap_member():
    # col1 gap in m2 -> m2 skips that allele node
    rows = ["ACGT", "A-GT", "ACGT"]
    gfa, paths = V.msa_to_gfa(rows, ["m1", "m2", "m3"])
    assert paths["m1"] == ["1", "2", "3"]
    assert paths["m2"] == ["1", "3"]        # skips the gap region
    assert paths["m3"] == ["1", "2", "3"]
    links = {(l.split("\t")[1], l.split("\t")[3]) for l in gfa.splitlines() if l.startswith("L\t")}
    assert ("1", "3") in links              # m2's skip link


def test_msa_to_gfa_all_identical_single_node():
    rows = ["ACGT", "ACGT"]
    gfa, paths = V.msa_to_gfa(rows, ["a", "b"])
    assert paths["a"] == ["1"] and paths["b"] == ["1"]
    assert "S\t1\tACGT" in gfa
    assert not any(l.startswith("L\t") for l in gfa.splitlines())  # no links, one node
```

- [ ] **Step 2: Run to verify failure**

Run: `PYTHONHASHSEED=0 /home/juanfra/.local/bin/pytest bench/soto/test_soto_dna_vg.py -v -k "msa_to_gfa"`
Expected: FAIL — `AttributeError: module 'soto_dna_vg' has no attribute 'msa_to_gfa'`.

- [ ] **Step 3: Implement the abpoa wrapper + the pure converter**

Append to `bench/soto/soto_dna_vg.py`:

```python
def abpoa_msa(seqs):
    """Aligned gapped rows for the members (str inputs; bytes silently yield all-N). Returns the first
    len(seqs) rows of abpoa's column-MSA (excludes any appended consensus row)."""
    import pyabpoa
    aligner = pyabpoa.msa_aligner()
    res = aligner.msa([s.upper() for s in seqs], out_cons=False, out_msa=True)
    rows = list(res.msa_seq)[:len(seqs)]
    return [r.upper() for r in rows]


def msa_to_gfa(rows, names):
    """Column-MSA -> base-level variation-graph GFA. rows: equal-length uppercase gapped strings
    ('-' = gap), one per name. Returns (gfa_text, paths). Deterministic:
      - a maximal run of columns where all rows share the SAME non-gap base -> one shared node (all paths);
      - a maximal run of variant/gap columns -> one node per distinct gap-stripped allele (sorted); a
        member whose slice is all-gap skips the region."""
    assert rows and len({len(r) for r in rows}) == 1, "rows must be non-empty and equal length"
    assert len(rows) == len(names)
    m, L = len(rows), len(rows[0])

    def invariant(j):
        c0 = rows[0][j]
        return c0 != "-" and all(rows[i][j] == c0 for i in range(m))

    # segment columns into maximal same-class runs
    segments, j = [], 0
    while j < L:
        kind = "inv" if invariant(j) else "var"
        k = j + 1
        while k < L and ("inv" if invariant(k) else "var") == kind:
            k += 1
        segments.append((kind, j, k))
        j = k

    nodes, paths, nid = [], {n: [] for n in names}, 0
    for kind, a, b in segments:
        if kind == "inv":
            nid += 1
            node = str(nid)
            nodes.append((node, rows[0][a:b]))
            for n in names:
                paths[n].append(node)
        else:
            allele = {n: rows[i][a:b].replace("-", "") for i, n in enumerate(names)}
            node_of = {}
            for s in sorted(set(v for v in allele.values() if v)):
                nid += 1
                node_of[s] = str(nid)
                nodes.append((str(nid), s))
            for n in names:
                if allele[n]:
                    paths[n].append(node_of[allele[n]])

    links = set()
    for n in names:
        p = paths[n]
        for x, y in zip(p, p[1:]):
            links.add((x, y))

    out = ["H\tVN:Z:1.0"]
    for node, seq in nodes:
        out.append(f"S\t{node}\t{seq}")
    for x, y in sorted(links, key=lambda t: (int(t[0]), int(t[1]))):
        out.append(f"L\t{x}\t+\t{y}\t+\t0M")
    for n in names:
        out.append(f"P\t{n}\t{'+,'.join(paths[n])}+\t*")
    return "\n".join(out) + "\n", paths
```

- [ ] **Step 4: Run to verify pass**

Run: `PYTHONHASHSEED=0 /home/juanfra/.local/bin/pytest bench/soto/test_soto_dna_vg.py -v`
Expected: PASS (Task 1 + Task 2 tests).

- [ ] **Step 5: Commit**

```bash
git add bench/soto/soto_dna_vg.py bench/soto/test_soto_dna_vg.py
git commit -m "feat(soto-vg): abpoa column-MSA -> base-level variation-graph GFA"
```

---

### Task 3: Detection labels → node colours + legend

**Files:**
- Modify: `bench/soto/soto_dna_vg.py`
- Test: `bench/soto/test_soto_dna_vg.py` (add tests)

**Interfaces:**
- Consumes: `GREEN/RED/GREY`, `paths` from `msa_to_gfa`.
- Produces:
  - `load_detection(tsv_lines) -> dict[tuple[str,int,int], tuple[bool,str]]` — `(chrom,start,end) -> (detected, recovered_by)` from `soto_member_detection.tsv` (`family_id,gene,chrom,start,end,detected,recovered_by`; `detected == "Y"`).
  - `node_colour(members_through_node, detected_by_name) -> str` — pure. green if every member through the node is RNA-recovered, red if none is, grey otherwise (mixed/empty).
  - `colours_csv(paths, detected_by_name) -> str` — Bandage CSV `Node,Colour` (one row per node, colored by the members traversing it).
  - `legend_tsv(members, detected_by_name, recovered_by_name) -> str` — `gene\tlocus\tdetected\trecovered_by\tcolour`.

- [ ] **Step 1: Write the failing tests**

Append to `bench/soto/test_soto_dna_vg.py`:

```python
def test_load_detection():
    tsv = [
        "family_id\tgene\tchrom\tstart\tend\tdetected\trecovered_by",
        "ID_462\tSRGAP2C\tchr1\t121194173\t121402237\tY\tRNA-split",
        "ID_26\tSLC9B1P1\tchr16\t32804124\t32821138\tN\t",
    ]
    d = V.load_detection(tsv)
    assert d[("chr1", 121194173, 121402237)] == (True, "RNA-split")
    assert d[("chr16", 32804124, 32821138)] == (False, "")


def test_node_colour():
    det = {"a": True, "b": True, "c": False}
    assert V.node_colour({"a", "b"}, det) == V.GREEN     # all recovered
    assert V.node_colour({"c"}, det) == V.RED            # all DNA-only
    assert V.node_colour({"a", "c"}, det) == V.GREY      # mixed = shared/conserved
    assert V.node_colour(set(), det) == V.GREY           # empty


def test_colours_csv():
    paths = {"a": ["1", "2", "4"], "b": ["1", "3", "4"], "c": ["1", "2", "4"]}
    det = {"a": True, "b": False, "c": True}
    csv = V.colours_csv(paths, det)
    rows = dict(l.split(",") for l in csv.strip().splitlines()[1:])
    assert rows["1"] == V.GREY    # a,b,c -> mixed
    assert rows["2"] == V.GREEN   # a,c -> recovered
    assert rows["3"] == V.RED     # b -> DNA-only
    assert rows["4"] == V.GREY    # a,b,c -> mixed
    assert csv.splitlines()[0] == "Node,Colour"


def test_legend_tsv():
    members = [("SRGAP2C", "chr1", 100, 200)]
    det = {"SRGAP2C": True}
    rec = {"SRGAP2C": "RNA-split"}
    leg = V.legend_tsv(members, det, rec)
    assert leg.splitlines()[0] == "gene\tlocus\tdetected\trecovered_by\tcolour"
    assert "SRGAP2C\tchr1:100-200\tY\tRNA-split\t#1e8e3e" in leg
```

- [ ] **Step 2: Run to verify failure**

Run: `PYTHONHASHSEED=0 /home/juanfra/.local/bin/pytest bench/soto/test_soto_dna_vg.py -v -k "detection or node_colour or colours_csv or legend"`
Expected: FAIL — `AttributeError: ... 'load_detection'`.

- [ ] **Step 3: Implement colours + legend**

Append to `bench/soto/soto_dna_vg.py`:

```python
def load_detection(tsv_lines):
    """(chrom,start,end) -> (detected: bool, recovered_by: str)."""
    out = {}
    for i, ln in enumerate(tsv_lines):
        if i == 0:
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 6:
            continue
        out[(f[2], int(f[3]), int(f[4]))] = (f[5] == "Y", f[6] if len(f) > 6 else "")
    return out


def node_colour(members_through_node, detected_by_name):
    """green if every member through the node is RNA-recovered, red if none is, grey otherwise."""
    flags = {bool(detected_by_name.get(n, False)) for n in members_through_node}
    if flags == {True}:
        return GREEN
    if flags == {False}:
        return RED
    return GREY


def colours_csv(paths, detected_by_name):
    node_members = defaultdict(set)
    for n, p in paths.items():
        for node in p:
            node_members[node].add(n)
    lines = ["Node,Colour"]
    for node in sorted(node_members, key=int):
        lines.append(f"{node},{node_colour(node_members[node], detected_by_name)}")
    return "\n".join(lines) + "\n"


def legend_tsv(members, detected_by_name, recovered_by_name):
    lines = ["gene\tlocus\tdetected\trecovered_by\tcolour"]
    for gene, chrom, start, end in members:
        det = bool(detected_by_name.get(gene, False))
        lines.append("\t".join([
            gene, f"{chrom}:{start}-{end}", "Y" if det else "N",
            recovered_by_name.get(gene, ""), GREEN if det else RED]))
    return "\n".join(lines) + "\n"
```

- [ ] **Step 4: Run to verify pass**

Run: `PYTHONHASHSEED=0 /home/juanfra/.local/bin/pytest bench/soto/test_soto_dna_vg.py -v`
Expected: PASS (all tests through Task 3).

- [ ] **Step 5: Commit**

```bash
git add bench/soto/soto_dna_vg.py bench/soto/test_soto_dna_vg.py
git commit -m "feat(soto-vg): detection labels -> node colours + legend"
```

---

### Task 4: Per-family orchestration + presence check + main + captions

**Files:**
- Modify: `bench/soto/soto_dna_vg.py`
- Test: `bench/soto/test_soto_dna_vg.py` (add tests)

**Interfaces:**
- Consumes: all of Tasks 1–3.
- Produces:
  - `build_family(family_id, members, fa, detection) -> dict` — extract seqs (skip members absent from `fa`, log), abpoa MSA, `msa_to_gfa` (P-line names = gene), colours + legend; returns `{"family_id","n_members","n_present","gfa","colours","legend","missing":[gene,...]}`. The **presence check**: `n_present` = members that appear as a `P`-line; `missing` = members dropped (absent from `fa`) — never silently counted as present.
  - `main()` — build the flagships (`ID_462` SRGAP2, `ID_8` PMS2P, `ID_63`), write `<fam>.gfa`/`.colours.csv`/`.legend.tsv` under `OUT`, and an `index.md` with the honesty caption; print a one-line per-family summary (`family n_present/n_members`).

- [ ] **Step 1: Write the failing tests**

Append to `bench/soto/test_soto_dna_vg.py`:

```python
def test_build_family_presence_and_colours(monkeypatch):
    # isolate build_family's orchestration from abpoa (keeps the test hermetic — no pyabpoa needed).
    # these 3 seqs are equal-length with one SNP, so the identity "MSA" is the correct alignment.
    monkeypatch.setattr(V, "abpoa_msa", lambda seqs: list(seqs))
    fa = {"chr1:101-200": "ACGTACGT", "chr1:301-400": "ACGTACGT", "chr1:501-600": "ACGTTCGT"}
    members = [("g1", "chr1", 100, 200), ("g2", "chr1", 300, 400), ("g3", "chr1", 500, 600)]
    detection = {("chr1", 100, 200): (True, "RNA-split"),
                 ("chr1", 300, 400): (True, "projection"),
                 ("chr1", 500, 600): (False, "")}
    r = V.build_family("ID_test", members, fa, detection)
    assert r["n_members"] == 3 and r["n_present"] == 3 and r["missing"] == []
    # every gene is a P-line
    plines = {l.split("\t")[1] for l in r["gfa"].splitlines() if l.startswith("P\t")}
    assert plines == {"g1", "g2", "g3"}
    # g3 differs (T at col4) -> its unique node is red; shared nodes grey; g1/g2 identical
    assert V.RED in r["colours"] and V.GREY in r["colours"]
    assert "g3\tchr1:500-600\tN" in r["legend"]


def test_build_family_missing_member_logged_not_counted():
    fa = {"chr1:101-200": "ACGT"}   # g2 absent from fa
    members = [("g1", "chr1", 100, 200), ("g2", "chr9", 100, 200)]
    detection = {("chr1", 100, 200): (True, "RNA-split")}
    r = V.build_family("ID_x", members, fa, detection)
    assert r["n_members"] == 2 and r["n_present"] == 1 and r["missing"] == ["g2"]
```

- [ ] **Step 2: Run to verify failure**

Run: `PYTHONHASHSEED=0 /home/juanfra/.local/bin/pytest bench/soto/test_soto_dna_vg.py -v -k "build_family"`
Expected: FAIL — `AttributeError: ... 'build_family'`.

- [ ] **Step 3: Implement orchestration + main**

Append to `bench/soto/soto_dna_vg.py`:

```python
FLAGSHIPS = [("ID_462", "SRGAP2"), ("ID_8", "PMS2P"), ("ID_63", "mixed-recovery")]

CAPTION = (
    "DNA variation graph = the ceiling: all {n} Soto copies present as paths (Soto-corroborated, "
    "independent DNA-read-depth catalog). green = RNA recovered; red = DNA-only (K=0 exon-identity / "
    "silent / coverage). The VG REPRESENTS what is given; it does not 'detect' families. RNA recovers "
    "76.2% of this ceiling genome-wide; the gap is the decomposed identifiability floor, not a method failure."
)


def build_family(family_id, members, fa, detection):
    """Extract member seqs (skip+log those absent from fa), abpoa MSA, GFA, colours, legend, presence check."""
    present, missing, seqs, names = [], [], [], []
    for gene, chrom, start, end in members:
        try:
            seqs.append(member_seq(fa, chrom, start, end))
            names.append(gene)
            present.append((gene, chrom, start, end))
        except KeyError:
            missing.append(gene)
    if not present:
        return dict(family_id=family_id, n_members=len(members), n_present=0,
                    gfa="", colours="", legend="", missing=missing)
    rows = abpoa_msa(seqs) if len(seqs) > 1 else [seqs[0]]
    gfa, paths = msa_to_gfa(rows, names)
    det_by_gene = {g: detection.get((c, s, e), (False, ""))[0] for g, c, s, e in present}
    rec_by_gene = {g: detection.get((c, s, e), (False, ""))[1] for g, c, s, e in present}
    colours = colours_csv(paths, det_by_gene)
    legend = legend_tsv(present, det_by_gene, rec_by_gene)
    n_present = sum(1 for l in gfa.splitlines() if l.startswith("P\t"))   # checked, not assumed
    return dict(family_id=family_id, n_members=len(members), n_present=n_present,
                gfa=gfa, colours=colours, legend=legend, missing=missing)


def main():
    os.makedirs(OUT, exist_ok=True)
    bed = open(BED).read().splitlines()
    fa = read_fasta(MEMFA)
    detection = load_detection(open(DETECT).read().splitlines())
    index = ["# Soto DNA variation-graph ceiling — flagship families\n"]
    for family_id, label in FLAGSHIPS:
        members = parse_family_members(bed, family_id)
        r = build_family(family_id, members, fa, detection)
        base = f"{OUT}/{family_id}"
        open(f"{base}.gfa", "w").write(r["gfa"])
        open(f"{base}.colours.csv", "w").write(r["colours"])
        open(f"{base}.legend.tsv", "w").write(r["legend"])
        miss = f"  (MISSING from graph: {r['missing']})" if r["missing"] else ""
        print(f"{family_id} ({label}): {r['n_present']}/{r['n_members']} copies as paths{miss}")
        index.append(f"## {family_id} — {label}: {r['n_present']}/{r['n_members']} copies present as paths")
        index.append(f"`{family_id}.gfa` + `{family_id}.colours.csv` (Bandage). " + CAPTION.format(n=r["n_present"]))
        if r["missing"]:
            index.append(f"> honesty: {len(r['missing'])} member(s) absent from the graph: {r['missing']}")
    open(f"{OUT}/index.md", "w").write("\n\n".join(index) + "\n")
    print(f"wrote {OUT}/ (gfa + colours.csv + legend.tsv per family + index.md)")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run to verify pass**

Run: `PYTHONHASHSEED=0 /home/juanfra/.local/bin/pytest bench/soto/test_soto_dna_vg.py -v`
Expected: PASS (all unit tests).

- [ ] **Step 5: Commit**

```bash
git add bench/soto/soto_dna_vg.py bench/soto/test_soto_dna_vg.py
git commit -m "feat(soto-vg): per-family orchestration + presence check + main + captions"
```

---

### Task 5: Standalone real-abpoa smoke (miniforge python) + render the flagships

**Files:**
- Create: `bench/soto/smoke_soto_dna_vg.py`

**Interfaces:**
- Consumes: `build_family` + the on-disk `soto_members.fa`, `80_fams.chr.bed`, `soto_member_detection.tsv`. Runs abpoa on one real family via **miniforge python** (the interpreter that has `pyabpoa`); single foreground process, WSL2-safe.

- [ ] **Step 1: Write the standalone smoke script**

Create `bench/soto/smoke_soto_dna_vg.py`:

```python
#!/usr/bin/env python3
"""Real-abpoa smoke for soto_dna_vg. Run via the interpreter that has pyabpoa:
  /home/juanfra/miniforge3/bin/python bench/soto/smoke_soto_dna_vg.py
Builds PMS2P (ID_8) end-to-end and asserts the DNA-ceiling presence claim. Exit 0 = pass.
(A standalone script, not a pytest test, because the pytest interpreter lacks pyabpoa.)"""
import sys
sys.path.insert(0, "bench/soto")
import soto_dna_vg as V

bed = open("bench/soto/80_fams.chr.bed").read().splitlines()
fa = V.read_fasta(V.MEMFA)
detection = V.load_detection(open("bench/soto/soto_member_detection.tsv").read().splitlines())
members = V.parse_family_members(bed, "ID_8")            # PMS2P, moderate ~5-15 kb members
assert len(members) >= 5, members
r = V.build_family("ID_8", members, fa, detection)
assert r["missing"] == [], f"members missing from graph: {r['missing']}"   # the DNA-ceiling claim
assert r["n_present"] == r["n_members"], (r["n_present"], r["n_members"])
kinds = {l[0] for l in r["gfa"].splitlines() if l}
assert {"H", "S", "P"} <= kinds, kinds
assert r["colours"].splitlines()[0] == "Node,Colour" and len(r["colours"].splitlines()) > 1
plines = {l.split("\t")[1] for l in r["gfa"].splitlines() if l.startswith("P\t")}
assert plines == {g for g, _, _, _ in members}, plines
n_nodes = sum(1 for l in r["gfa"].splitlines() if l.startswith("S"))
print(f"[smoke] PMS2P ID_8: {r['n_present']}/{r['n_members']} copies as paths, {n_nodes} nodes  OK")
```

- [ ] **Step 2: Run the smoke via miniforge python (foreground, single abpoa — WSL2 crash rule)**

Run: `cd /mnt/c/Users/jfris/Desktop/Rustle && PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/soto/smoke_soto_dna_vg.py`
Expected: `[smoke] PMS2P ID_8: N/N copies as paths, <k> nodes  OK` (exit 0). If a member is unexpectedly `missing`, do NOT delete it to force success — check its `soto_members.fa` header (`chrom:start+1-end`) and report the mismatch honestly.

- [ ] **Step 3: Render the full flagship set (foreground, miniforge python)**

Run: `cd /mnt/c/Users/jfris/Desktop/Rustle && PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/soto/soto_dna_vg.py`
Expected: three `n_present/n_members` lines (SRGAP2 4/4, PMS2P all present, ID_63 with green+red), files under `/home/juanfra/winloci_scratch/soto_vg/`. SRGAP2 members span ~260 kb; if its graph is large, note it — do not truncate. This generates the actual Bandage artifacts.

- [ ] **Step 4: Run the hermetic unit suite (via .local pytest)**

Run: `cd /mnt/c/Users/jfris/Desktop/Rustle && PYTHONHASHSEED=0 /home/juanfra/.local/bin/pytest bench/soto/test_soto_dna_vg.py -v`
Expected: PASS (all Task 1–4 unit tests; no pyabpoa needed).

- [ ] **Step 5: Commit**

```bash
git add bench/soto/smoke_soto_dna_vg.py
git commit -m "test(soto-vg): standalone real-abpoa smoke + flagship render"
```

---

## Optional follow-ups (not in this plan's core; documented for later)

- **Exon overlay** (`--exon-overlay`): tag graph nodes exonic vs intronic from `A119b.t2t.bam` spliced-read blocks (CIGAR `N`), to make the K=0 exon-identity floor explicit. Adds per-locus BAM fetches.
- **SLC9B1P1 (`ID_26`)**: the "DNA exceeds even Soto's annotation" extreme — its beyond-Soto copies come from `soto_dna_oracle_prototype.tsv` (extract loci from `chm13v2.0.fa`).
- **`--builder vg`**: emit via `vg msga -f <fam>.members.fa` (deprecated) for an authentic-`vg`-branded GFA.

## Notes for the implementer

- **Run pytest from the repo root** so `sys.path.insert(0, "bench/soto")` resolves.
- **pyabpoa needs `str`, not bytes** (bytes yield all-`N`); `abpoa_msa` upper-cases and slices to the member rows.
- **WSL2 crash rule:** abpoa runs foreground, one family at a time; member seqs from the 11 MB `soto_members.fa`; no whole-genome, no `copy_assign`; outputs under `winloci_scratch/soto_vg/`.
- **Honesty:** the presence check is real — `n_present` counts actual `P`-lines; a member missing from `soto_members.fa` is reported in `missing`, never silently counted toward the "100%".
```
