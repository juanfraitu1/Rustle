# VG Paralog-Copy-Recovery Evaluation Protocol — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a deterministic, reproducible genome-wide protocol that measures whether rustle `--vg` recovers RefSeq paralog-family transcript copies StringTie2 misses, using SQANTI3 as the structural matcher plus a tool-independent primary-alignment authenticity guard.

**Architecture:** Six isolated units behind one driver. U1 generates transcriptomes (StringTie2 + rustle arms) once on the full BAM; U2 runs SQANTI3 per arm; U3 builds a tool-independent paralog universe; U4 extracts FSM/ISM recovery sets; U5 applies the primary-support guard; U6 produces the head-to-head report. The expensive runs (U1/U2) feed file artifacts to the cheap deterministic Python scoring units (U3–U6), which are TDD'd on toy fixtures so the logic is verified before any genome-wide run.

**Tech Stack:** Bash orchestration; Python 3.11 (pandas, pysam) for scoring with pytest; SQANTI3 (conda env); StringTie 3.0.1; rustle (`--vg`); samtools; minimap2; gffread.

**Spec:** `docs/superpowers/specs/2026-06-08-vg-copy-recovery-eval-protocol-design.md`

**Data (pinned):** `SCRATCH=/home/juanfra/winloci_scratch` → `GGO.bam` (1.6 G, 380,369 secondary alns), `GGO.fasta` (3.4 G), `GGO_tx.gff` (403 M). rustle binary `target/release/rustle`; StringTie `tools/stringtie/stringtie`.

**SQANTI3 category strings (exact, used in U4):** `full-splice_match` (FSM), `incomplete-splice_match` (ISM), `novel_in_catalog` (NIC), `novel_not_in_catalog` (NNC), `genic`, `antisense`, `fusion`, `intergenic`, `genic_intron`.

---

## File structure (all new, under `bench/copy_recovery_eval/`)

- `config.sh` — all paths, thresholds, tool flags, version pins; emits a config hash.
- `00_install_sqanti3.sh` — clone SQANTI3 + create its conda env + smoke test.
- `10_generate_transcriptomes.sh` — U1: run arms on the BAM → per-arm GTF + `provenance.json`; canonical-sort GTFs.
- `20_run_sqanti3.sh` — U2: `sqanti3_qc.py` per arm → `<arm>_classification.txt`.
- `lib_eval.py` — U3–U6 pure functions (universe, recovery extraction, guard helpers, head-to-head) — importable, no I/O side effects in the core functions.
- `30_build_universe.py` — U3 CLI wrapping `lib_eval.build_universe` → `universe.tsv`.
- `40_recovery_sets.py` — U4 CLI wrapping `lib_eval.recovery_set` → `<arm>_recovery.tsv`.
- `50_authenticity_guard.py` — U5 CLI: primary-support via pysam → `authenticity.tsv`.
- `60_headtohead.py` — U6 CLI wrapping `lib_eval.head_to_head` → `headtohead.tsv` + `REPORT.md`.
- `run_protocol.sh` — driver chaining 00→60.
- `tests/test_lib_eval.py` — pytest unit tests for U3/U4/U6 logic on toy fixtures.
- `tests/test_authenticity.py` — pytest unit test for U5 primary-counting on a toy BAM.
- `tests/fixtures/` — toy `classification.txt`, toy annotation, toy SAM→BAM.

**Why a `lib_eval.py` core:** the set logic (universe, recovery, head-to-head) is the protocol's intellectual content and must be unit-tested in isolation from SQANTI3/BAM I/O. CLIs are thin wrappers.

---

## Task 1: Scaffold + pinned config

**Files:**
- Create: `bench/copy_recovery_eval/config.sh`
- Create: `bench/copy_recovery_eval/README.md`

- [ ] **Step 1: Write `config.sh`**

```bash
#!/usr/bin/env bash
# Pinned configuration for the VG copy-recovery evaluation protocol.
# Source this from every stage script: `source "$(dirname "$0")/config.sh"`.
set -euo pipefail

# ── Data ────────────────────────────────────────────────────────────────────
SCRATCH="${SCRATCH:-/home/juanfra/winloci_scratch}"
BAM="${BAM:-$SCRATCH/GGO.bam}"
GENOME_FASTA="${GENOME_FASTA:-$SCRATCH/GGO.fasta}"
REF_GFF="${REF_GFF:-$SCRATCH/GGO_tx.gff}"

# ── Tools ───────────────────────────────────────────────────────────────────
ROOT="${ROOT:-/mnt/c/Users/jfris/Desktop/Rustle}"
RUSTLE="${RUSTLE:-$ROOT/target/release/rustle}"
STRINGTIE="${STRINGTIE:-$ROOT/tools/stringtie/stringtie}"
SAMTOOLS="${SAMTOOLS:-samtools}"
MINIMAP2="${MINIMAP2:-minimap2}"
GFFREAD="${GFFREAD:-gffread}"
SQANTI3_DIR="${SQANTI3_DIR:-$HOME/tools/SQANTI3}"
SQANTI3_ENV="${SQANTI3_ENV:-sqanti3}"

# ── Output ──────────────────────────────────────────────────────────────────
OUTDIR="${OUTDIR:-$ROOT/bench/copy_recovery_eval/results}"

# ── Arms (1=on) ─────────────────────────────────────────────────────────────
ARM_STRINGTIE=1                 # baseline (headline competitor)
ARM_RUSTLE_VG=1                 # headline method (decisive gate ON)
ARM_RUSTLE_VG_RAW="${ARM_RUSTLE_VG_RAW:-0}"      # optional ablation (gate OFF)
ARM_RUSTLE_PRIMARY="${ARM_RUSTLE_PRIMARY:-0}"    # optional attribution arm (-L primary-only)

# ── rustle --vg headline flags (win stack + decisive gate ON) ───────────────
RUSTLE_VG_ENV="RUSTLE_VG_TANDEM=1 RUSTLE_VG_TANDEM_PRIMARY_JUNCTIONS=1 RUSTLE_VG_DECISIVE_GATE=1 RUSTLE_VG_DECISIVE_GATE_MIN_PRIM=4"
RUSTLE_VG_RAW_ENV="RUSTLE_VG_TANDEM=1 RUSTLE_VG_TANDEM_PRIMARY_JUNCTIONS=1"   # no gate

# ── Thresholds (pinned) ─────────────────────────────────────────────────────
MIN_IDENTITY="${MIN_IDENTITY:-0.90}"   # U3 paralog identity floor
MIN_COV_FRAC="${MIN_COV_FRAC:-0.50}"   # U3 alignment length-coverage floor
GUARD_K="${GUARD_K:-1}"                # U5 min primary reads for authentic
# Restrict expensive runs to one chrom for smoke tests (empty = whole genome):
CHROM_SUBSET="${CHROM_SUBSET:-}"

# ── Provenance ──────────────────────────────────────────────────────────────
config_hash() {
  { echo "$MIN_IDENTITY $MIN_COV_FRAC $GUARD_K $RUSTLE_VG_ENV $RUSTLE_VG_RAW_ENV"; } | sha1sum | cut -d' ' -f1
}
mkdir -p "$OUTDIR"
```

- [ ] **Step 2: Write a one-paragraph `README.md`** describing the protocol, the run order (`run_protocol.sh`), and that `config.sh` holds all pins. Keep it short.

- [ ] **Step 3: Commit**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
git add bench/copy_recovery_eval/config.sh bench/copy_recovery_eval/README.md
git commit -m "feat(copy-eval): pinned config + scaffold for VG copy-recovery protocol

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: U3 — paralog universe (pure logic, TDD first)

The universe is tool-independent: RefSeq genes in families with ≥2 copies, confirmed by sequence identity. Logic split into a pure function (parse pre-computed self-alignment + annotation → universe rows) that is unit-tested, and a CLI that produces the inputs.

**Files:**
- Create: `bench/copy_recovery_eval/lib_eval.py`
- Create: `bench/copy_recovery_eval/tests/test_lib_eval.py`
- Create: `bench/copy_recovery_eval/tests/fixtures/` (fixtures created inline in tests)

- [ ] **Step 1: Write the failing test**

Create `tests/test_lib_eval.py`:

```python
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import lib_eval


def test_build_universe_groups_paralog_families():
    # gene->transcripts map (from annotation)
    gene_tx = {
        "RABL2A": ["rna-A1", "rna-A2"],
        "RABL2B": ["rna-B1"],
        "SOLO":   ["rna-S1"],
    }
    # paralog pairs from self-alignment that PASS identity/coverage thresholds:
    paralog_pairs = [("RABL2A", "RABL2B")]   # SOLO has no paralog
    uni = lib_eval.build_universe(gene_tx, paralog_pairs)
    # Every transcript of RABL2A and RABL2B is in U, with a shared family id.
    tx_in = {r["transcript_id"] for r in uni}
    assert tx_in == {"rna-A1", "rna-A2", "rna-B1"}
    fam = {r["transcript_id"]: r["family_id"] for r in uni}
    assert fam["rna-A1"] == fam["rna-A2"] == fam["rna-B1"]
    # SOLO excluded.
    assert "rna-S1" not in tx_in
    # n_family_copies counts GENES in the family (2 here).
    assert all(r["n_family_copies"] == 2 for r in uni)


def test_build_universe_transitive_closure():
    # A~B and B~C  => one family {A,B,C}
    gene_tx = {"A": ["a"], "B": ["b"], "C": ["c"], "D": ["d"]}
    paralog_pairs = [("A", "B"), ("B", "C")]
    uni = lib_eval.build_universe(gene_tx, paralog_pairs)
    fam = {r["transcript_id"]: r["family_id"] for r in uni}
    assert fam["a"] == fam["b"] == fam["c"]
    assert "d" not in fam
    assert all(r["n_family_copies"] == 3 for r in uni)
```

- [ ] **Step 2: Run to verify failure**

Run: `cd bench/copy_recovery_eval && python -m pytest tests/test_lib_eval.py -q`
Expected: FAIL (`lib_eval` has no `build_universe`).

- [ ] **Step 3: Implement `build_universe` in `lib_eval.py`**

```python
"""Pure scoring logic for the VG copy-recovery protocol (U3, U4, U6).
No file or subprocess I/O in these functions — they take parsed data and
return plain dict rows, so they are unit-testable in isolation."""

def _union_find_families(genes, pairs):
    parent = {g: g for g in genes}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in pairs:
        if a in parent and b in parent:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb
    # group genes by root
    fam = {}
    for g in genes:
        fam.setdefault(find(g), []).append(g)
    return fam  # root_gene -> [genes]


def build_universe(gene_tx, paralog_pairs):
    """gene_tx: {gene_id: [transcript_id,...]}. paralog_pairs: [(geneA,geneB),...]
    that already passed identity/coverage thresholds.
    Returns rows: {transcript_id, gene_id, family_id, n_family_copies} for genes
    in families of >=2 genes. family_id is the deterministic min gene_id in the family."""
    families = _union_find_families(list(gene_tx.keys()), paralog_pairs)
    rows = []
    for _root, genes in families.items():
        if len(genes) < 2:
            continue
        family_id = min(genes)
        n_copies = len(genes)
        for g in sorted(genes):
            for tx in gene_tx[g]:
                rows.append({
                    "transcript_id": tx,
                    "gene_id": g,
                    "family_id": family_id,
                    "n_family_copies": n_copies,
                })
    rows.sort(key=lambda r: (r["family_id"], r["gene_id"], r["transcript_id"]))
    return rows
```

- [ ] **Step 4: Run to verify pass**

Run: `python -m pytest tests/test_lib_eval.py -q`
Expected: 2 passed.

- [ ] **Step 5: Commit**

```bash
git add bench/copy_recovery_eval/lib_eval.py bench/copy_recovery_eval/tests/test_lib_eval.py
git commit -m "feat(copy-eval): U3 paralog-universe set logic (TDD)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: U4 — recovery-set extraction (pure logic, TDD)

Parse a SQANTI3 classification table (already loaded as a list of dict rows) into per-arm recovery: which reference transcripts have an FSM (and separately ISM) query isoform, intersected with U.

**Files:**
- Modify: `bench/copy_recovery_eval/lib_eval.py`
- Modify: `bench/copy_recovery_eval/tests/test_lib_eval.py`

- [ ] **Step 1: Write the failing test** (append to `tests/test_lib_eval.py`)

```python
def test_recovery_set_fsm_ism_within_universe():
    # SQANTI3 classification rows (only the columns we use).
    classif = [
        {"isoform": "q1", "structural_category": "full-splice_match",       "associated_transcript": "rna-A1"},
        {"isoform": "q2", "structural_category": "incomplete-splice_match", "associated_transcript": "rna-A2"},
        {"isoform": "q3", "structural_category": "novel_not_in_catalog",    "associated_transcript": "novel"},
        {"isoform": "q4", "structural_category": "full-splice_match",       "associated_transcript": "rna-S1"},  # not in U
    ]
    universe_tx = {"rna-A1", "rna-A2", "rna-B1"}
    rec = lib_eval.recovery_set(classif, universe_tx)
    # rna-A1 FSM; rna-A2 ISM only; rna-S1 excluded (not in U).
    assert rec["rna-A1"] == {"fsm": True, "ism": False}
    assert rec["rna-A2"] == {"fsm": False, "ism": True}
    assert "rna-S1" not in rec
    assert "novel" not in rec
```

- [ ] **Step 2: Run to verify failure**

Run: `python -m pytest tests/test_lib_eval.py::test_recovery_set_fsm_ism_within_universe -q`
Expected: FAIL (no `recovery_set`).

- [ ] **Step 3: Implement `recovery_set`** (append to `lib_eval.py`)

```python
FSM = "full-splice_match"
ISM = "incomplete-splice_match"


def recovery_set(classification_rows, universe_tx):
    """classification_rows: list of dicts with 'structural_category' and
    'associated_transcript'. universe_tx: set of reference transcript ids (U).
    Returns {ref_transcript: {'fsm': bool, 'ism': bool}} restricted to U."""
    rec = {}
    for row in classification_rows:
        cat = row.get("structural_category", "")
        ref = row.get("associated_transcript", "")
        if ref not in universe_tx:
            continue
        if cat not in (FSM, ISM):
            continue
        cur = rec.setdefault(ref, {"fsm": False, "ism": False})
        if cat == FSM:
            cur["fsm"] = True
        elif cat == ISM:
            cur["ism"] = True
    return rec
```

- [ ] **Step 4: Run to verify pass**

Run: `python -m pytest tests/test_lib_eval.py -q`
Expected: 3 passed (2 prior + this).

- [ ] **Step 5: Commit**

```bash
git add bench/copy_recovery_eval/lib_eval.py bench/copy_recovery_eval/tests/test_lib_eval.py
git commit -m "feat(copy-eval): U4 FSM/ISM recovery-set extraction (TDD)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: U6 — head-to-head logic (pure logic, TDD)

Combine per-arm recovery sets + authenticity flags into the headline: copies recovered by rustle-VG but not StringTie2, split FSM/ISM × authentic/phantom.

**Files:**
- Modify: `bench/copy_recovery_eval/lib_eval.py`
- Modify: `bench/copy_recovery_eval/tests/test_lib_eval.py`

- [ ] **Step 1: Write the failing test** (append)

```python
def test_head_to_head_splits_authentic_and_phantom():
    rustle = {
        "rna-A1": {"fsm": True, "ism": False},   # ST misses, authentic -> WIN
        "rna-B1": {"fsm": True, "ism": False},   # ST misses, phantom   -> phantom
        "rna-C1": {"fsm": True, "ism": False},   # ST also has it       -> not rustle-only
    }
    stringtie = {
        "rna-C1": {"fsm": True, "ism": False},
    }
    authentic = {"rna-A1": True, "rna-B1": False, "rna-C1": True}
    fam = {"rna-A1": "FAM1", "rna-B1": "FAM2", "rna-C1": "FAM3"}
    res = lib_eval.head_to_head(rustle, stringtie, authentic, fam)
    assert res["rustle_only_fsm_authentic"] == ["rna-A1"]
    assert res["rustle_only_fsm_phantom"] == ["rna-B1"]
    assert res["n_win"] == 1
    assert res["n_phantom"] == 1
```

- [ ] **Step 2: Run to verify failure**

Run: `python -m pytest tests/test_lib_eval.py::test_head_to_head_splits_authentic_and_phantom -q`
Expected: FAIL (no `head_to_head`).

- [ ] **Step 3: Implement `head_to_head`** (append)

```python
def head_to_head(rustle_rec, stringtie_rec, authentic, family_of):
    """rustle_rec / stringtie_rec: {ref_tx: {'fsm','ism'}}. authentic: {ref_tx: bool}
    (only meaningful for rustle recoveries). family_of: {ref_tx: family_id}.
    Returns the headline split as sorted lists + counts."""
    def fsm_set(rec):
        return {tx for tx, v in rec.items() if v["fsm"]}
    r_fsm = fsm_set(rustle_rec)
    s_fsm = fsm_set(stringtie_rec)
    rustle_only = r_fsm - s_fsm
    auth = sorted(tx for tx in rustle_only if authentic.get(tx, False))
    phantom = sorted(tx for tx in rustle_only if not authentic.get(tx, False))
    return {
        "rustle_only_fsm_authentic": auth,
        "rustle_only_fsm_phantom": phantom,
        "n_win": len(auth),
        "n_phantom": len(phantom),
        "families_won": sorted({family_of.get(tx, "NA") for tx in auth}),
    }
```

- [ ] **Step 4: Run to verify pass**

Run: `python -m pytest tests/test_lib_eval.py -q`
Expected: 4 passed.

- [ ] **Step 5: Commit**

```bash
git add bench/copy_recovery_eval/lib_eval.py bench/copy_recovery_eval/tests/test_lib_eval.py
git commit -m "feat(copy-eval): U6 head-to-head authentic/phantom split (TDD)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: U5 — primary-support authenticity guard (pysam, TDD on toy BAM)

For a recovered copy, count reads whose **primary** alignment overlaps the copy's exons. Tool-independent (uses SAM flags only).

**Files:**
- Create: `bench/copy_recovery_eval/tests/test_authenticity.py`
- Create: `bench/copy_recovery_eval/50_authenticity_guard.py`

- [ ] **Step 1: Write the failing test**

Create `tests/test_authenticity.py` (builds a tiny BAM with pysam: one primary + one secondary read over a region):

```python
import sys, os, subprocess
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import pysam
import importlib.util

def _load_guard():
    path = os.path.join(os.path.dirname(__file__), "..", "50_authenticity_guard.py")
    spec = importlib.util.spec_from_file_location("guard", path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m

def _make_bam(tmp):
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 10000}]}
    bam = os.path.join(tmp, "toy.bam")
    with pysam.AlignmentFile(bam, "wb", header=header) as out:
        # primary read overlapping 1000-1100
        a = pysam.AlignedSegment(out.header)
        a.query_name = "r1"; a.flag = 0; a.reference_id = 0
        a.reference_start = 1000; a.mapping_quality = 60
        a.cigartuples = [(0, 100)]; a.query_sequence = "A"*100
        a.query_qualities = [30]*100
        out.write(a)
        # secondary read overlapping the SAME region (flag 0x100)
        b = pysam.AlignedSegment(out.header)
        b.query_name = "r2"; b.flag = 0x100; b.reference_id = 0
        b.reference_start = 1010; b.mapping_quality = 0
        b.cigartuples = [(0, 100)]; b.query_sequence = "A"*100
        b.query_qualities = [30]*100
        out.write(b)
    pysam.index(bam)
    return bam

def test_primary_support_counts_only_primary(tmp_path):
    guard = _load_guard()
    bam = _make_bam(str(tmp_path))
    # exons overlapping the reads' region -> 1 primary (r1), secondary r2 ignored
    n = guard.primary_support(bam, "chr1", [(1000, 1100)])
    assert n == 1
    # a region with no reads -> 0
    assert guard.primary_support(bam, "chr1", [(5000, 5100)]) == 0
```

- [ ] **Step 2: Run to verify failure**

Run: `python -m pytest tests/test_authenticity.py -q`
Expected: FAIL (no `50_authenticity_guard.py`).

- [ ] **Step 3: Implement `50_authenticity_guard.py`**

```python
#!/usr/bin/env python3
"""U5: tool-independent authenticity guard. For each recovered copy, count reads
whose PRIMARY alignment overlaps the copy's exons (SAM flags only). A recovery is
authentic iff that count >= K."""
import argparse, csv, sys
import pysam


def primary_support(bam_path, chrom, exons):
    """Count distinct reads with a PRIMARY alignment overlapping any exon."""
    n = 0
    seen = set()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for (start, end) in exons:
            for read in bam.fetch(chrom, start, end):
                if read.is_secondary or read.is_supplementary or read.is_unmapped:
                    continue
                if read.query_name in seen:
                    continue
                seen.add(read.query_name)
                n += 1
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--recovery", required=True, help="rustle-VG recovery.tsv (ref_transcript,family_id,fsm,ism)")
    ap.add_argument("--exons", required=True, help="TSV: transcript_id,chrom,exon_starts(csv),exon_ends(csv)")
    ap.add_argument("--k", type=int, default=1)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    # exon map
    exon_map = {}
    with open(args.exons) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            starts = [int(x) for x in row["exon_starts"].split(",") if x]
            ends = [int(x) for x in row["exon_ends"].split(",") if x]
            exon_map[row["transcript_id"]] = (row["chrom"], list(zip(starts, ends)))

    with open(args.recovery) as f, open(args.out, "w", newline="") as o:
        rows = list(csv.DictReader(f, delimiter="\t"))
        w = csv.writer(o, delimiter="\t")
        w.writerow(["ref_transcript", "family_id", "primary_support", "authentic"])
        out_rows = []
        for r in rows:
            tx = r["ref_transcript"]
            if tx not in exon_map:
                sys.stderr.write(f"[guard] no exons for {tx}, skipping\n"); continue
            chrom, exons = exon_map[tx]
            n = primary_support(args.bam, chrom, exons)
            out_rows.append((tx, r["family_id"], n, "true" if n >= args.k else "false"))
        out_rows.sort(key=lambda x: (x[1], x[0]))   # deterministic
        for row in out_rows:
            w.writerow(row)


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run to verify pass**

Run: `python -m pytest tests/test_authenticity.py -q`
Expected: 1 passed. (Requires `pysam` in the python env; if missing: `mamba install -n base -c bioconda pysam` or use the `phasing_eval` env which has bioconda packages.)

- [ ] **Step 5: Commit**

```bash
git add bench/copy_recovery_eval/50_authenticity_guard.py bench/copy_recovery_eval/tests/test_authenticity.py
git commit -m "feat(copy-eval): U5 primary-support authenticity guard (TDD on toy BAM)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 6: U3 CLI — build universe from real annotation + self-alignment

Wraps `lib_eval.build_universe` with the real input producers: parse the RefSeq annotation for gene→transcript and per-gene sequences, self-align with minimap2, threshold, and emit `universe.tsv` + a per-transcript exon table (reused by U5).

**Files:**
- Create: `bench/copy_recovery_eval/30_build_universe.py`

- [ ] **Step 1: Implement `30_build_universe.py`**

```python
#!/usr/bin/env python3
"""U3 CLI: build the paralog universe from the RefSeq annotation + genome.
Steps: (1) parse annotation -> gene_tx + per-transcript exons; (2) write a
representative-sequence FASTA per gene (longest transcript) via gffread; (3)
self-align with minimap2; (4) keep pairs >= identity/coverage; (5) build_universe.
Outputs universe.tsv and exons.tsv (transcript_id,chrom,exon_starts,exon_ends)."""
import argparse, csv, subprocess, sys, os
sys.path.insert(0, os.path.dirname(__file__))
import lib_eval


def parse_annotation(gff):
    """Return (gene_tx: {gene:[tx]}, exons: {tx:(chrom,[(s,e)])}, tx_gene:{tx:gene})."""
    gene_tx, exons, tx_gene = {}, {}, {}
    with open(gff) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            chrom, feat, start, end, attrs = c[0], c[2], int(c[3]) - 1, int(c[4]), c[8]
            kv = {}
            for field in attrs.replace('"', "").split(";"):
                field = field.strip()
                if not field:
                    continue
                if "=" in field:
                    k, v = field.split("=", 1)
                elif " " in field:
                    k, v = field.split(" ", 1)
                else:
                    continue
                kv[k.strip()] = v.strip()
            if feat in ("transcript", "mRNA"):
                tx = kv.get("transcript_id") or kv.get("ID")
                gene = kv.get("gene_id") or kv.get("gene") or kv.get("Parent")
                if tx and gene:
                    gene_tx.setdefault(gene, []).append(tx)
                    tx_gene[tx] = gene
                    exons.setdefault(tx, (chrom, []))
            elif feat == "exon":
                tx = kv.get("transcript_id") or kv.get("Parent")
                if tx:
                    ch, lst = exons.setdefault(tx, (chrom, []))
                    lst.append((start, end))
    return gene_tx, exons, tx_gene


def write_exons_tsv(exons, path):
    with open(path, "w", newline="") as o:
        w = csv.writer(o, delimiter="\t")
        w.writerow(["transcript_id", "chrom", "exon_starts", "exon_ends"])
        for tx in sorted(exons):
            chrom, lst = exons[tx]
            lst = sorted(lst)
            w.writerow([tx, chrom,
                        ",".join(str(s) for s, _ in lst),
                        ",".join(str(e) for _, e in lst)])


def gene_paralog_pairs(gene_fasta, min_identity, min_cov_frac, minimap2):
    """All-vs-all minimap2 of per-gene representative seqs; return passing gene pairs."""
    paf = subprocess.run(
        [minimap2, "-x", "asm20", "-X", "--no-long-join", "-c", gene_fasta, gene_fasta],
        capture_output=True, text=True, check=True).stdout
    pairs = set()
    for line in paf.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        q, qlen, tname = f[0], int(f[1]), f[5]
        nmatch, alen = int(f[9]), int(f[10])
        if q == tname:
            continue
        ident = nmatch / alen if alen else 0.0
        cov = alen / qlen if qlen else 0.0
        if ident >= min_identity and cov >= min_cov_frac:
            pairs.add(tuple(sorted((q, tname))))
    return sorted(pairs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref-gff", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--min-identity", type=float, default=0.90)
    ap.add_argument("--min-cov-frac", type=float, default=0.50)
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--gffread", default="gffread")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    gene_tx, exons, tx_gene = parse_annotation(args.ref_gff)
    write_exons_tsv(exons, os.path.join(args.outdir, "exons.tsv"))

    # Representative gene FASTA: longest transcript per gene, named by GENE id.
    tx_len = {tx: sum(e - s for s, e in lst) for tx, (_, lst) in exons.items()}
    rep_tx = {}
    for gene, txs in gene_tx.items():
        best = max((t for t in txs if t in tx_len), key=lambda t: tx_len[t], default=None)
        if best:
            rep_tx[best] = gene
    # transcript FASTA from gffread, then rename headers tx->gene
    tx_fa = os.path.join(args.outdir, "tx.fa")
    subprocess.run([args.gffread, "-w", tx_fa, "-g", args.genome, args.ref_gff], check=True)
    gene_fa = os.path.join(args.outdir, "gene_rep.fa")
    with open(tx_fa) as f, open(gene_fa, "w") as o:
        keep = False
        for line in f:
            if line.startswith(">"):
                name = line[1:].split()[0]
                keep = name in rep_tx
                if keep:
                    o.write(f">{rep_tx[name]}\n")
            elif keep:
                o.write(line)

    pairs = gene_paralog_pairs(gene_fa, args.min_identity, args.min_cov_frac, args.minimap2)
    uni = lib_eval.build_universe(gene_tx, pairs)
    with open(os.path.join(args.outdir, "universe.tsv"), "w", newline="") as o:
        w = csv.writer(o, delimiter="\t")
        w.writerow(["transcript_id", "gene_id", "family_id", "n_family_copies"])
        for r in uni:
            w.writerow([r["transcript_id"], r["gene_id"], r["family_id"], r["n_family_copies"]])
    sys.stderr.write(f"[U3] {len(uni)} transcripts in {len({r['family_id'] for r in uni})} paralog families\n")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Smoke test on a small annotation slice**

Build a tiny GFF slice (first ~2000 lines of the real annotation) and run the CLI on it against the genome, confirming `universe.tsv` and `exons.tsv` are produced with the expected header. (Full-genome run happens in the driver, Task 9.)

Run:
```bash
cd bench/copy_recovery_eval && source config.sh
head -2000 "$REF_GFF" > /tmp/ref_slice.gff
python 30_build_universe.py --ref-gff /tmp/ref_slice.gff --genome "$GENOME_FASTA" \
  --min-identity "$MIN_IDENTITY" --min-cov-frac "$MIN_COV_FRAC" \
  --minimap2 "$MINIMAP2" --gffread "$GFFREAD" --outdir /tmp/u3_smoke
head -3 /tmp/u3_smoke/universe.tsv; head -3 /tmp/u3_smoke/exons.tsv
```
Expected: both files exist with correct headers (universe may be empty on a tiny slice — that's fine; the test is that it runs and the exon table is populated).

- [ ] **Step 3: Commit**

```bash
git add bench/copy_recovery_eval/30_build_universe.py
git commit -m "feat(copy-eval): U3 CLI - universe from annotation + minimap2 self-align

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 7: U4 CLI + U6 CLI — recovery sets and report from real SQANTI3 output

**Files:**
- Create: `bench/copy_recovery_eval/40_recovery_sets.py`
- Create: `bench/copy_recovery_eval/60_headtohead.py`

- [ ] **Step 1: Implement `40_recovery_sets.py`**

```python
#!/usr/bin/env python3
"""U4 CLI: read a SQANTI3 *_classification.txt + universe.tsv -> <arm>_recovery.tsv."""
import argparse, csv, sys, os
sys.path.insert(0, os.path.dirname(__file__))
import lib_eval


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--classification", required=True)
    ap.add_argument("--universe", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    with open(args.universe) as f:
        uni_rows = list(csv.DictReader(f, delimiter="\t"))
    universe_tx = {r["transcript_id"] for r in uni_rows}
    fam_of = {r["transcript_id"]: r["family_id"] for r in uni_rows}

    with open(args.classification) as f:
        classif = list(csv.DictReader(f, delimiter="\t"))
    rec = lib_eval.recovery_set(classif, universe_tx)

    with open(args.out, "w", newline="") as o:
        w = csv.writer(o, delimiter="\t")
        w.writerow(["ref_transcript", "family_id", "fsm", "ism"])
        for tx in sorted(rec):
            w.writerow([tx, fam_of.get(tx, "NA"),
                        "true" if rec[tx]["fsm"] else "false",
                        "true" if rec[tx]["ism"] else "false"])


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Implement `60_headtohead.py`**

```python
#!/usr/bin/env python3
"""U6 CLI: combine per-arm recovery.tsv + authenticity.tsv -> headtohead.tsv + REPORT.md."""
import argparse, csv, sys, os
sys.path.insert(0, os.path.dirname(__file__))
import lib_eval


def _load_recovery(path):
    rec, fam = {}, {}
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rec[r["ref_transcript"]] = {"fsm": r["fsm"] == "true", "ism": r["ism"] == "true"}
            fam[r["ref_transcript"]] = r["family_id"]
    return rec, fam


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rustle-recovery", required=True)
    ap.add_argument("--stringtie-recovery", required=True)
    ap.add_argument("--authenticity", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    rustle_rec, fam = _load_recovery(args.rustle_recovery)
    stringtie_rec, fam2 = _load_recovery(args.stringtie_recovery)
    fam.update(fam2)
    authentic = {}
    with open(args.authenticity) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            authentic[r["ref_transcript"]] = r["authentic"] == "true"

    res = lib_eval.head_to_head(rustle_rec, stringtie_rec, authentic, fam)

    h2h = os.path.join(args.outdir, "headtohead.tsv")
    with open(h2h, "w", newline="") as o:
        w = csv.writer(o, delimiter="\t")
        w.writerow(["category", "ref_transcript", "family_id"])
        for tx in res["rustle_only_fsm_authentic"]:
            w.writerow(["rustle_only_fsm_authentic", tx, fam.get(tx, "NA")])
        for tx in res["rustle_only_fsm_phantom"]:
            w.writerow(["rustle_only_fsm_phantom", tx, fam.get(tx, "NA")])

    with open(os.path.join(args.outdir, "REPORT.md"), "w") as o:
        o.write("# VG copy-recovery head-to-head\n\n")
        o.write(f"- Authentic rustle-VG-only FSM recoveries (StringTie2 misses): **{res['n_win']}**\n")
        o.write(f"- Phantom (secondary-only) rustle-VG-only FSM recoveries: **{res['n_phantom']}**\n")
        o.write(f"- Families won: {', '.join(res['families_won']) or 'none'}\n\n")
        o.write("See `headtohead.tsv` for per-transcript rows.\n")
    sys.stderr.write(f"[U6] win={res['n_win']} phantom={res['n_phantom']} -> {h2h}\n")


if __name__ == "__main__":
    main()
```

- [ ] **Step 3: Smoke test the two CLIs on toy inputs**

```bash
cd bench/copy_recovery_eval
printf "transcript_id\tgene_id\tfamily_id\tn_family_copies\nrna-A1\tGA\tFAM1\t2\nrna-B1\tGB\tFAM1\t2\n" > /tmp/uni.tsv
printf "isoform\tstructural_category\tassociated_transcript\nq1\tfull-splice_match\trna-A1\nq2\tfull-splice_match\trna-B1\n" > /tmp/cl_rustle.txt
printf "isoform\tstructural_category\tassociated_transcript\nq9\tfull-splice_match\trna-B1\n" > /tmp/cl_st.txt
python 40_recovery_sets.py --classification /tmp/cl_rustle.txt --universe /tmp/uni.tsv --out /tmp/rec_rustle.tsv
python 40_recovery_sets.py --classification /tmp/cl_st.txt --universe /tmp/uni.tsv --out /tmp/rec_st.tsv
printf "ref_transcript\tfamily_id\tprimary_support\tauthentic\nrna-A1\tFAM1\t7\ttrue\nrna-B1\tFAM1\t9\ttrue\n" > /tmp/auth.tsv
python 60_headtohead.py --rustle-recovery /tmp/rec_rustle.tsv --stringtie-recovery /tmp/rec_st.tsv --authenticity /tmp/auth.tsv --outdir /tmp
cat /tmp/headtohead.tsv
```
Expected: `headtohead.tsv` lists `rna-A1` under `rustle_only_fsm_authentic` (StringTie has rna-B1, so only rna-A1 is rustle-only), and REPORT.md says win=1.

- [ ] **Step 4: Commit**

```bash
git add bench/copy_recovery_eval/40_recovery_sets.py bench/copy_recovery_eval/60_headtohead.py
git commit -m "feat(copy-eval): U4/U6 CLIs - recovery sets + head-to-head report

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 8: SQANTI3 install + U2 runner

**Files:**
- Create: `bench/copy_recovery_eval/00_install_sqanti3.sh`
- Create: `bench/copy_recovery_eval/20_run_sqanti3.sh`

- [ ] **Step 1: Implement `00_install_sqanti3.sh`**

```bash
#!/usr/bin/env bash
# Install SQANTI3 into a dedicated conda env. SQANTI3 ships its own env yml.
source "$(dirname "$0")/config.sh"
set -euo pipefail

if [ ! -d "$SQANTI3_DIR" ]; then
  mkdir -p "$(dirname "$SQANTI3_DIR")"
  git clone https://github.com/ConesaLab/SQANTI3.git "$SQANTI3_DIR"
fi
cd "$SQANTI3_DIR"
# Pin to a known-good release tag for reproducibility.
git fetch --tags
git checkout v5.2.1 2>/dev/null || git checkout v5.2 2>/dev/null || true

# Create env from the repo's yml if not present.
if ! mamba env list | grep -qE "^\s*$SQANTI3_ENV\s"; then
  ENV_YML=$(ls SQANTI3.conda_env.yml environment.yml 2>/dev/null | head -1)
  mamba env create -n "$SQANTI3_ENV" -f "$ENV_YML"
fi

# Smoke test: sqanti3_qc.py --version (or --help) must run in the env.
mamba run -n "$SQANTI3_ENV" python "$SQANTI3_DIR/sqanti3_qc.py" --version 2>&1 | head -3 \
  || mamba run -n "$SQANTI3_ENV" python "$SQANTI3_DIR/sqanti3_qc.py" --help 2>&1 | head -5
echo "[SQANTI3] installed at $SQANTI3_DIR, env $SQANTI3_ENV"
```

- [ ] **Step 2: Run the install + record the version**

Run: `bash bench/copy_recovery_eval/00_install_sqanti3.sh 2>&1 | tail -20`
Expected: clone + env create succeed; the smoke line prints SQANTI3 help/version. If the pinned tag's env fails to solve, fall back to the latest release tag and record whichever installed. **If the install genuinely cannot be made to work in this environment, STOP and report BLOCKED** — do not fake it; the protocol needs a real SQANTI3.

- [ ] **Step 3: Implement `20_run_sqanti3.sh`**

```bash
#!/usr/bin/env bash
# U2: run SQANTI3 QC for one arm's GTF against the RefSeq reference + genome.
# Usage: 20_run_sqanti3.sh <arm_name> <arm_gtf>
source "$(dirname "$0")/config.sh"
set -euo pipefail
ARM="${1:?arm name}"; GTF="${2:?arm gtf}"
OUT="$OUTDIR/sqanti3/$ARM"; mkdir -p "$OUT"

# SQANTI3 prefers a GTF reference; convert the RefSeq GFF once (cached).
REF_GTF="$OUTDIR/ref.gtf"
if [ ! -s "$REF_GTF" ]; then
  "$GFFREAD" "$REF_GFF" -T -o "$REF_GTF"
fi

mamba run -n "$SQANTI3_ENV" python "$SQANTI3_DIR/sqanti3_qc.py" \
  "$GTF" "$REF_GTF" "$GENOME_FASTA" \
  -d "$OUT" -o "$ARM" --report skip --skipORF 2>&1 | tail -5

# The classification file SQANTI3 writes:
ls "$OUT/${ARM}_classification.txt"
echo "[U2] $ARM classification -> $OUT/${ARM}_classification.txt"
```

- [ ] **Step 4: Commit** (the U2 runner is verified end-to-end in Task 9's driver, since it needs a real GTF)

```bash
git add bench/copy_recovery_eval/00_install_sqanti3.sh bench/copy_recovery_eval/20_run_sqanti3.sh
git commit -m "feat(copy-eval): SQANTI3 install + U2 per-arm QC runner

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 9: U1 transcriptome generation + full driver + oracle validation

**Files:**
- Create: `bench/copy_recovery_eval/10_generate_transcriptomes.sh`
- Create: `bench/copy_recovery_eval/run_protocol.sh`

- [ ] **Step 1: Implement `10_generate_transcriptomes.sh`**

```bash
#!/usr/bin/env bash
# U1: generate per-arm transcriptomes from the full BAM (no windows).
# Canonical-sorts each GTF (rustle GTF line order has rayon nondeterminism).
source "$(dirname "$0")/config.sh"
set -euo pipefail
GENDIR="$OUTDIR/gtf"; mkdir -p "$GENDIR"

# Optional chrom subset for smoke testing (CHROM_SUBSET in config).
INBAM="$BAM"
if [ -n "$CHROM_SUBSET" ]; then
  INBAM="$OUTDIR/subset.bam"
  "$SAMTOOLS" view -b "$BAM" "$CHROM_SUBSET" > "$INBAM"
  "$SAMTOOLS" index "$INBAM"
fi

canon_sort() { # deterministic GTF order: by chrom,start,feature
  sort -t$'\t' -k1,1 -k4,4n -k3,3 "$1" > "$1.sorted" && mv "$1.sorted" "$1"; }

PROV="$OUTDIR/provenance.json"
echo "{" > "$PROV"
echo "  \"config_hash\": \"$(config_hash)\"," >> "$PROV"
echo "  \"rustle_commit\": \"$(cd "$ROOT" && git rev-parse HEAD)\"," >> "$PROV"
echo "  \"stringtie_version\": \"$($STRINGTIE --version 2>&1 | head -1)\"," >> "$PROV"
echo "  \"bam\": \"$INBAM\"" >> "$PROV"
echo "}" >> "$PROV"

if [ "$ARM_STRINGTIE" = 1 ]; then
  # De novo (no -G): fair "find it from reads alone" comparison.
  "$STRINGTIE" -L -o "$GENDIR/stringtie.gtf" "$INBAM"
  canon_sort "$GENDIR/stringtie.gtf"
fi
if [ "$ARM_RUSTLE_VG" = 1 ]; then
  env $RUSTLE_VG_ENV "$RUSTLE" --vg --vg-snp --genome-fasta "$GENOME_FASTA" \
    -L "$INBAM" -o "$GENDIR/rustle_vg.gtf"
  canon_sort "$GENDIR/rustle_vg.gtf"
fi
if [ "$ARM_RUSTLE_VG_RAW" = 1 ]; then
  env $RUSTLE_VG_RAW_ENV "$RUSTLE" --vg --vg-snp --genome-fasta "$GENOME_FASTA" \
    -L "$INBAM" -o "$GENDIR/rustle_vg_raw.gtf"
  canon_sort "$GENDIR/rustle_vg_raw.gtf"
fi
if [ "$ARM_RUSTLE_PRIMARY" = 1 ]; then
  "$RUSTLE" -L "$INBAM" -o "$GENDIR/rustle_primary.gtf"
  canon_sort "$GENDIR/rustle_primary.gtf"
fi
echo "[U1] GTFs in $GENDIR"; ls -la "$GENDIR"
```

- [ ] **Step 2: Implement `run_protocol.sh`** (the deterministic driver, no LLM in the loop)

```bash
#!/usr/bin/env bash
# Full protocol driver: U1 -> U2 -> U3 -> U4 -> U5 -> U6.
source "$(dirname "$0")/config.sh"
set -euo pipefail
D="$(dirname "$0")"

bash "$D/10_generate_transcriptomes.sh"

# U3 universe (tool-independent; once).
"$SQANTI3_PY_ENV_PYTHON_OK" 2>/dev/null || true
python "$D/30_build_universe.py" --ref-gff "$REF_GFF" --genome "$GENOME_FASTA" \
  --min-identity "$MIN_IDENTITY" --min-cov-frac "$MIN_COV_FRAC" \
  --minimap2 "$MINIMAP2" --gffread "$GFFREAD" --outdir "$OUTDIR"

# U2 + U4 per arm.
declare -A GTFS=( [stringtie]="$OUTDIR/gtf/stringtie.gtf" [rustle_vg]="$OUTDIR/gtf/rustle_vg.gtf" )
for arm in "${!GTFS[@]}"; do
  [ -s "${GTFS[$arm]}" ] || continue
  bash "$D/20_run_sqanti3.sh" "$arm" "${GTFS[$arm]}"
  python "$D/40_recovery_sets.py" \
    --classification "$OUTDIR/sqanti3/$arm/${arm}_classification.txt" \
    --universe "$OUTDIR/universe.tsv" --out "$OUTDIR/${arm}_recovery.tsv"
done

# U5 authenticity (rustle-VG only).
python "$D/50_authenticity_guard.py" --bam "$BAM" \
  --recovery "$OUTDIR/rustle_vg_recovery.tsv" --exons "$OUTDIR/exons.tsv" \
  --k "$GUARD_K" --out "$OUTDIR/authenticity.tsv"

# U6 head-to-head.
python "$D/60_headtohead.py" \
  --rustle-recovery "$OUTDIR/rustle_vg_recovery.tsv" \
  --stringtie-recovery "$OUTDIR/stringtie_recovery.tsv" \
  --authenticity "$OUTDIR/authenticity.tsv" --outdir "$OUTDIR"

echo "[DONE] see $OUTDIR/REPORT.md and $OUTDIR/headtohead.tsv"
```

Note: remove the stray `$SQANTI3_PY_ENV_PYTHON_OK` line — it was a copy artifact; the Python scripts run in the base/phasing env that has pysam/pandas. Ensure `python` resolves to an interpreter with `pysam` (use `mamba run -n phasing_eval python ...` if base lacks it; set a `PY` variable in config and use it throughout).

- [ ] **Step 3: Smoke-run the whole protocol on a single chromosome**

Pick the chromosome containing a known family (RABL2). Set `CHROM_SUBSET` to it in the environment and run:
```bash
cd bench/copy_recovery_eval
CHROM_SUBSET="<chrom_with_RABL2>" bash run_protocol.sh 2>&1 | tail -30
cat results/REPORT.md
```
Expected: completes end-to-end; `REPORT.md` and `headtohead.tsv` produced. (Determine the RABL2 chromosome from the annotation: `grep -m1 RABL2 "$REF_GFF"`.)

- [ ] **Step 4: Oracle validation**

Confirm the protocol's self-validation oracles on the smoke run (or a targeted region run):
- **RABL2** appears in `headtohead.tsv` under `rustle_only_fsm_authentic` (rustle-VG recovers a copy StringTie2 misses, and it passes the primary-support guard).
- **DAZ3** (if in the subset/region) appears as `rustle_only_fsm_phantom` OR is absent from authentic recoveries (guard flags prim=0).

If either oracle fails, STOP and report — it means a unit is mis-wired (likely U4 associated_transcript joining or U5 exon coordinates). Document the finding; do not paper over it.

- [ ] **Step 5: Commit**

```bash
git add bench/copy_recovery_eval/10_generate_transcriptomes.sh bench/copy_recovery_eval/run_protocol.sh
git commit -m "feat(copy-eval): U1 transcriptome generation + full deterministic driver + oracles

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 10: Genome-wide run + provenance + determinism check

**Files:**
- Modify: `bench/copy_recovery_eval/README.md` (record the headline result + how to reproduce)

- [ ] **Step 1: Full genome-wide run**

With `CHROM_SUBSET` empty (whole genome), run `bash run_protocol.sh`. This is expensive (1.6 G BAM; rustle `--vg` + StringTie2 + 2× SQANTI3). Expect a long runtime; run it in the background and capture logs. Expected artifacts: `results/REPORT.md`, `results/headtohead.tsv`, per-arm classifications, `universe.tsv`, `authenticity.tsv`, `provenance.json`.

- [ ] **Step 2: Determinism check**

Re-run the deterministic scoring stages (U3–U6) on the same U1/U2 artifacts and confirm `headtohead.tsv` is byte-identical:
```bash
cd bench/copy_recovery_eval
cp results/headtohead.tsv /tmp/h2h_a.tsv
python 30_build_universe.py --ref-gff "$REF_GFF" --genome "$GENOME_FASTA" --min-identity "$MIN_IDENTITY" --min-cov-frac "$MIN_COV_FRAC" --minimap2 "$MINIMAP2" --gffread "$GFFREAD" --outdir results
# re-run U4/U5/U6 ... then:
diff /tmp/h2h_a.tsv results/headtohead.tsv && echo "DETERMINISTIC OK"
```
Expected: `DETERMINISTIC OK`. (Note any SQANTI3-level nondeterminism separately if the upstream classification differs across runs.)

- [ ] **Step 3: Record the headline in README + commit**

Write the genome-wide result (n_win, n_phantom, families won, global SQANTI3 category table) into `README.md` with the exact reproduce command and `provenance.json` hash.

```bash
git add bench/copy_recovery_eval/README.md
git commit -m "docs(copy-eval): genome-wide headline result + reproduction record

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Self-review (completed during planning)

**Spec coverage:** Claim → Task 9/10 head-to-head. SQANTI3 matcher → Task 8 (U2). Primary-support guard → Task 5 (U5). Genome-wide no-windows → Task 9 (full BAM, no extraction). Paralog universe (annotation+identity, tool-independent) → Tasks 2+6 (U3). FSM/ISM recovery → Tasks 3+7 (U4). StringTie2-only baseline + optional arms → config.sh + Task 9. rustle gate-on headline + raw ablation → config (`RUSTLE_VG_ENV` vs `RUSTLE_VG_RAW_ENV`). Determinism/provenance → Task 1 config_hash, Task 9 canon_sort + provenance.json, Task 10 determinism check. Oracles (RABL2 authentic, DAZ3 phantom) → Task 9 Step 4. Outputs dir → all units write to `$OUTDIR`. Global SQANTI3 distributions → recorded in Task 10 README (note: add a small column-count helper if needed; the classification files already contain per-isoform categories to tabulate).

**Placeholder scan:** No TBDs. The one stray artifact line (`$SQANTI3_PY_ENV_PYTHON_OK`) in `run_protocol.sh` is explicitly called out in Task 9 Step 2 to be removed and replaced with a `PY` variable — the implementer must wire `PY` (interpreter with pysam) in `config.sh` and use it for all `python` calls. **Implementer: add `PY="${PY:-mamba run -n phasing_eval python}"` to config.sh and replace bare `python` in the driver/CLIs invocations accordingly.**

**Type/interface consistency:** `universe.tsv` columns (transcript_id, gene_id, family_id, n_family_copies) produced by Task 6, consumed by Tasks 7 (U4) + the family map. `<arm>_recovery.tsv` columns (ref_transcript, family_id, fsm, ism) produced by Task 7, consumed by U5 (Task 5) + U6 (Task 7). `exons.tsv` (transcript_id, chrom, exon_starts, exon_ends) produced by Task 6, consumed by U5 (Task 5). `authenticity.tsv` (ref_transcript, family_id, primary_support, authentic) produced by U5, consumed by U6. `lib_eval.build_universe/recovery_set/head_to_head` signatures consistent between tests (Tasks 2–4) and CLIs (Tasks 6–7). SQANTI3 category strings (`full-splice_match`, `incomplete-splice_match`) consistent between U4 logic and tests.

**Known risks flagged for execution:** (1) SQANTI3 install is the top risk — Task 8 Step 2 says report BLOCKED rather than fake it. (2) Annotation attribute parsing (Task 6 `parse_annotation`) depends on the real GFF's attribute keys — the implementer must inspect a few lines of `GGO_tx.gff` and adjust the `transcript_id`/`gene_id`/`Parent` key handling to match (the parser already tries common GFF3/GTF spellings; verify on real data in Task 6 Step 2). (3) SQANTI3 `associated_transcript` ID format must match the annotation transcript IDs used in `universe.tsv`/`exons.tsv` — if SQANTI3 strips/renames IDs, U4's join will silently miss; verify on the smoke run (Task 9) via the RABL2 oracle. (4) Full genome-wide rustle `--vg` runtime/memory on the 1.6 G BAM is unproven at whole-genome scale — Task 10 may need per-chromosome generation concatenated (documented as a deviation if so; chromosome boundaries are not the fragile per-locus windows the spec forbids, but cross-chromosome paralog families would then need the whole-genome run — note this trade-off).
```
