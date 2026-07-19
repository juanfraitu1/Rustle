# Non-Coding-Aware Promotion Track — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Recover credible non-coding / lncRNA reference-absent collapses that the protein-BLASTx funnel drops, by promoting them on a collapse-quality bar into a new additive `gw_noncoding_copies.json` track.

**Architecture:** One new standalone script `bench/promote_noncoding.py` = a pure gate (`classify_noncoding`) + pure label/record helpers + a pure orchestration function (`promote_all`) + a thin I/O shell (`main`) that runs a single `minimap2` re-score of the 734 pre-ORF-gate consensuses in `cons.fa`. Nothing existing is modified or regenerated.

**Tech Stack:** Python 3 (stdlib only in the module; `minimap2` shelled from `main`); pytest 9.0.3. Tests run from the repo root: `PYTHONHASHSEED=0 pytest bench/test_promote_noncoding.py -v`.

## Global Constraints

- **Purely additive:** no existing file is modified; no existing output is regenerated; the coding catalog and all prior outputs stay **byte-identical**.
- **Protein/ORF are labels, never gates.**
- **Every promoted record carries `copy_vs_allele = "candidate-DNA-needed"` and `status = "flagged-reference-divergent-candidate"`, set unconditionally.**
- The gate is a **pure function** with an I/O shell (the module must import without `pysam`/`pyabpoa`/`minimap2`, so unit tests are hermetic). Therefore `MIN_DEPTH` and `orf_aa` are defined *locally* in `promote_noncoding.py` (copying the values/logic from `hidden_copy_scan`/`promote_genomewide`, which pull heavy imports), while `load_sedef`/`sedef_partners` are imported from the stdlib-only `copy_vs_allele_structural`.
- Thresholds are **named module-level constants** (single-edit recall/precision knob): `GID_LO=0.60`, `GID_HI=0.97`, `COV_MIN=0.90`, `MIN_COLS=8`, `AF_LO=0.15`, `AF_HI=0.60`, `MIN_ALT=5`, `MIN_DEPTH=8`.
- Deterministic; single `minimap2` process; **no `copy_assign`**; outputs under `winloci_scratch/refabsent/gw_promoted/`, not `/tmp`.
- `~REF`, `artifact`, `thin`, `not-own-locus` exclusions are **logged** (tally), not silently dropped.

---

### Task 1: Module scaffold + constants + the pure `classify_noncoding` gate

**Files:**
- Create: `bench/promote_noncoding.py`
- Test: `bench/test_promote_noncoding.py`

**Interfaces:**
- Produces: `classify_noncoding(rec: dict) -> tuple[bool, str, str]` returning `(promote, call, reason)`. `rec` keys: `genome_id` (float), `genome_cov` (float), `own_locus` (bool), `alt_cols` (int), `alt_read_fraction` (float), `alt_reads` (int), `n_primary` (int). `call` ∈ {`noncoding-candidate`, `~REF`, `artifact`, `chimera`, `not-own-locus`, `thin`}; only `noncoding-candidate` has `promote=True`.
- Produces module constants `GID_LO, GID_HI, COV_MIN, MIN_COLS, AF_LO, AF_HI, MIN_ALT, MIN_DEPTH`.

- [ ] **Step 1: Write the failing tests**

Create `bench/test_promote_noncoding.py`:

```python
#!/usr/bin/env python3
"""Tests for the non-coding-aware promotion track (bench/promote_noncoding.py).
Hermetic: the module must import without pysam/pyabpoa/minimap2.
Run: PYTHONHASHSEED=0 pytest bench/test_promote_noncoding.py -v"""
import os
import sys
sys.path.insert(0, "bench")
import promote_noncoding as PN


def _rec(**kw):
    base = dict(genome_id=0.93, genome_cov=1.02, own_locus=True,
                alt_cols=82, alt_read_fraction=0.318, alt_reads=61, n_primary=192)
    base.update(kw)
    return base


def test_flagship_promotes():
    promote, call, reason = PN.classify_noncoding(_rec())
    assert promote is True
    assert call == "noncoding-candidate"
    assert reason  # non-empty


def test_artifact_low_identity_rejected():
    promote, call, _ = PN.classify_noncoding(_rec(genome_id=0.28))
    assert promote is False and call == "artifact"


def test_near_reference_rejected():
    promote, call, _ = PN.classify_noncoding(_rec(genome_id=0.99))
    assert promote is False and call == "~REF"


def test_chimeric_low_coverage_rejected():
    promote, call, _ = PN.classify_noncoding(_rec(genome_cov=0.5))
    assert promote is False and call == "chimera"


def test_not_own_locus_rejected():
    promote, call, _ = PN.classify_noncoding(_rec(own_locus=False))
    assert promote is False and call == "not-own-locus"


def test_thin_columns_rejected():
    promote, call, _ = PN.classify_noncoding(_rec(alt_cols=2))
    assert promote is False and call == "thin"


def test_unbalanced_fraction_rejected():
    promote, call, _ = PN.classify_noncoding(_rec(alt_read_fraction=0.05))
    assert promote is False and call == "thin"


def test_low_depth_rejected():
    promote, call, _ = PN.classify_noncoding(_rec(alt_reads=3))
    assert promote is False and call == "thin"


def test_boundaries_inclusive():
    # genome_id == GID_LO promotes; == GID_HI is ~REF (excluded)
    assert PN.classify_noncoding(_rec(genome_id=PN.GID_LO))[0] is True
    assert PN.classify_noncoding(_rec(genome_id=PN.GID_HI))[1] == "~REF"
    # genome_cov == COV_MIN promotes; alt_cols == MIN_COLS promotes
    assert PN.classify_noncoding(_rec(genome_cov=PN.COV_MIN))[0] is True
    assert PN.classify_noncoding(_rec(alt_cols=PN.MIN_COLS))[0] is True
    # alt_read_fraction at both band edges promotes
    assert PN.classify_noncoding(_rec(alt_read_fraction=PN.AF_LO))[0] is True
    assert PN.classify_noncoding(_rec(alt_read_fraction=PN.AF_HI))[0] is True
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `PYTHONHASHSEED=0 pytest bench/test_promote_noncoding.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'promote_noncoding'` (module not created yet).

- [ ] **Step 3: Create the module with constants + the gate**

Create `bench/promote_noncoding.py`:

```python
#!/usr/bin/env python3
"""Non-coding-aware promotion of reference-absent copies: promote credible collapses on a
COLLAPSE-QUALITY bar (divergence band + genome coverage + own-locus + balanced co-segregation)
instead of a protein-BLASTx hit. Protein/ORF are LABELS, not gates. Purely additive: re-scores the
734 pre-ORF-gate consensuses already in gw_promoted/cons.fa and writes gw_noncoding_copies.json.
Every record is a FLAGGED reference-divergent candidate (copy_vs_allele='candidate-DNA-needed');
copy-vs-allele is not resolvable from RNA and needs DNA parCN. See
docs/superpowers/specs/2026-07-18-noncoding-promotion-design.md."""
import os
import sys
import json
import subprocess
from collections import Counter, defaultdict

sys.path.insert(0, "bench")
from copy_vs_allele_structural import load_sedef, sedef_partners  # stdlib-only module

# ---- paths (overridable in main via argparse) ----
MM2 = "/home/juanfra/miniforge3/bin/minimap2"
FASTA = "/home/juanfra/winloci_scratch/GGO.fasta"
CAT = "/home/juanfra/winloci_scratch/refabsent"
OUT = f"{CAT}/gw_promoted"
GFF = "/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_genomic.gff"
SEDEF = next((p for p in ["/mnt/c/Users/jfris/Desktop/final.bed",
                          "/mnt/c/Users/jfris/Desktop/Rustle/final.bed"] if os.path.exists(p)), None)

# ---- the collapse-quality bar (single-edit recall/precision knob) ----
GID_LO, GID_HI = 0.60, 0.97   # genome_id band: <LO = repeat/chimera artifact, >=HI = ~reference allele
COV_MIN = 0.90                # consensus must map near full-length (not a chimeric fragment)
MIN_COLS = 8                  # co-segregation breadth (alt columns)
AF_LO, AF_HI = 0.15, 0.60     # overall alt-read fraction band (balanced collapse)
MIN_ALT = 5                   # minimum alt reads
MIN_DEPTH = 8                 # minimum primary reads (matches hidden_copy_scan.MIN_DEPTH)

# ---- honesty rail (set unconditionally on every promoted record) ----
CANDIDATE = "candidate-DNA-needed"
STATUS = "flagged-reference-divergent-candidate"


def classify_noncoding(rec):
    """Pure collapse-quality gate. rec keys: genome_id, genome_cov, own_locus (bool),
    alt_cols, alt_read_fraction, alt_reads, n_primary.
    Returns (promote: bool, call: str, reason: str). Only call=='noncoding-candidate' promotes."""
    gid = rec["genome_id"]
    cov = rec["genome_cov"]
    if gid >= GID_HI:
        return (False, "~REF", f"genome_id {gid:.3f} >= {GID_HI} (essentially the reference allele)")
    if gid < GID_LO:
        return (False, "artifact", f"genome_id {gid:.3f} < {GID_LO} (repeat/chimera, not a divergent copy)")
    if cov < COV_MIN:
        return (False, "chimera", f"genome_cov {cov:.2f} < {COV_MIN} (partial/chimeric consensus)")
    if not rec["own_locus"]:
        return (False, "not-own-locus", "best genome hit is not the candidate's own locus")
    if rec["alt_cols"] < MIN_COLS:
        return (False, "thin", f"alt_cols {rec['alt_cols']} < {MIN_COLS} (too few co-segregating columns)")
    if not (AF_LO <= rec["alt_read_fraction"] <= AF_HI):
        return (False, "thin", f"alt_read_fraction {rec['alt_read_fraction']:.3f} outside "
                               f"[{AF_LO},{AF_HI}] (unbalanced)")
    if rec["alt_reads"] < MIN_ALT or rec["n_primary"] < MIN_DEPTH:
        return (False, "thin", f"depth alt_reads={rec['alt_reads']} (min {MIN_ALT}) / "
                               f"n_primary={rec['n_primary']} (min {MIN_DEPTH})")
    return (True, "noncoding-candidate",
            f"divergent ({100*(1-gid):.1f}%) full-length (cov {cov:.2f}) balanced collapse "
            f"({rec['alt_cols']} cols, af {rec['alt_read_fraction']:.2f}) at own locus")
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `PYTHONHASHSEED=0 pytest bench/test_promote_noncoding.py -v`
Expected: PASS (10 tests).

- [ ] **Step 5: Commit**

```bash
git add bench/promote_noncoding.py bench/test_promote_noncoding.py
git commit -m "feat(noncoding): pure classify_noncoding collapse-quality gate + constants"
```

---

### Task 2: Label helpers + record assembly with the honesty rail

**Files:**
- Modify: `bench/promote_noncoding.py`
- Test: `bench/test_promote_noncoding.py` (add tests)

**Interfaces:**
- Consumes: `CANDIDATE`, `STATUS` (Task 1).
- Produces:
  - `orf_aa(seq: str) -> int` (longest ATG..stop ORF length in aa, both strands; copied from `promote_genomewide` to keep the module import-light).
  - `coding_potential(orf_aa: int, cons_len: int) -> str` — `"noncoding"` if `orf_aa*3/cons_len < 0.30` else `"coding-capable/novel-protein"`.
  - `load_gff_genes(path) -> dict[str, list[tuple[int,int,str]]]` — chrom → sorted `(start, end, biotype)` for `gene` features; `{}` if path missing.
  - `biotype_of(idx, chrom, start, end) -> str` — `"lncRNA"` / `"protein_coding-gene-body"` / `"pseudogene"` / `"intergenic"` / `"unknown"`.
  - `read_fasta(path) -> dict[str, str]` — header (first token after `>`) → sequence.
  - `write_tsv(records: list[dict], path) -> None` — fixed 21-column TSV mirror.
  - `build_record(cid, flag, hit, cons_len, orf, biotype, sedef_partner, protein, reason) -> dict` — assembles a promoted record; `hit` is `(tgt, pos, gid, gcov)`; sets `copy_vs_allele=CANDIDATE` and `status=STATUS` **unconditionally**.

- [ ] **Step 1: Write the failing tests**

Append to `bench/test_promote_noncoding.py`:

```python
def test_coding_potential_label():
    assert PN.coding_potential(102, 3342) == "noncoding"            # flagship: 9.2%
    assert PN.coding_potential(1028, 4915) == "coding-capable/novel-protein"  # 62.7%
    assert PN.coding_potential(0, 0) == "noncoding"                 # guard div-by-zero


def test_orf_aa_finds_frame():
    # ATG GCA GCA TAA  -> M A A * -> ORF "MAA" length 3
    assert PN.orf_aa("ATGGCAGCATAA") == 3


def test_biotype_of_overlap():
    idx = PN.load_gff_from_lines([
        "NC_1\tX\tgene\t100\t500\t.\t+\t.\tgene_biotype=lncRNA",
        "NC_1\tX\tgene\t2000\t2500\t.\t+\t.\tgene_biotype=protein_coding",
    ])
    assert PN.biotype_of(idx, "NC_1", 150, 200) == "lncRNA"
    assert PN.biotype_of(idx, "NC_1", 2100, 2200) == "protein_coding-gene-body"
    assert PN.biotype_of(idx, "NC_1", 900, 950) == "intergenic"  # annotation loaded, no overlap
    assert PN.biotype_of({}, "NC_1", 150, 200) == "unknown"      # no annotation


def test_read_fasta_roundtrip(tmp_path):
    p = tmp_path / "c.fa"
    p.write_text(">a\nACGT\nACGT\n>b\nTTTT\n")
    d = PN.read_fasta(str(p))
    assert d == {"a": "ACGTACGT", "b": "TTTT"}


def test_build_record_honesty_rail():
    flag = dict(chrom="NC_073236.2", start=139051025, end=139225258,
                n_alt_positions=82, alt_read_fraction=0.3177, n_alt_reads=61, n_primary_reads=192)
    r = PN.build_record("NC_073236.2_139051025", flag, ("NC_073236.2", 139171134, 0.931, 1.02),
                        3342, 102, "lncRNA", None, "not-tested", "reason text")
    assert r["copy_vs_allele"] == "candidate-DNA-needed"
    assert r["status"] == "flagged-reference-divergent-candidate"
    assert r["track"] == "noncoding"
    assert r["divergence"] == 6.9 and r["genome_id"] == 0.931
    assert r["coding_potential"] == "noncoding" and r["biotype"] == "lncRNA"


def test_write_tsv_has_all_columns(tmp_path):
    flag = dict(chrom="NC_1", start=1, end=9, n_alt_positions=10, alt_read_fraction=0.3,
                n_alt_reads=6, n_primary_reads=20)
    r = PN.build_record("NC_1_1", flag, ("NC_1", 1, 0.9, 1.0), 300, 50, "intergenic", None, "not-tested", "x")
    p = tmp_path / "o.tsv"
    PN.write_tsv([r], str(p))
    header = p.read_text().splitlines()[0].split("\t")
    assert header[0] == "cid" and "copy_vs_allele" in header and "status" in header and len(header) == 21
```

Note: `load_gff_from_lines` is a tiny test seam so `biotype_of` can be tested without a GFF file; add it in Step 3.

- [ ] **Step 2: Run to verify failure**

Run: `PYTHONHASHSEED=0 pytest bench/test_promote_noncoding.py -v -k "coding_potential or orf_aa or biotype or read_fasta or honesty or write_tsv"`
Expected: FAIL — `AttributeError: module 'promote_noncoding' has no attribute 'coding_potential'` (etc.).

- [ ] **Step 3: Implement the helpers**

Append to `bench/promote_noncoding.py` (after `classify_noncoding`):

```python
# ---- ORF label (copied from promote_genomewide to keep this module import-light) ----
_COMP = str.maketrans("ACGTN", "TGCAN")
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
from itertools import product as _product
_CODON = {"".join(c): _AAS[i] for i, c in enumerate(_product("TCAG", repeat=3))}


def orf_aa(seq):
    """Longest ATG..stop ORF length in aa, scanning all 3 frames on both strands."""
    best = 0
    for strand in (seq, seq.translate(_COMP)[::-1]):
        for f in range(3):
            prot = "".join(_CODON.get(strand[i:i+3], "X") for i in range(f, len(strand)-2, 3))
            for p in prot.split("*"):
                k = p.find("M")
                if k >= 0:
                    best = max(best, len(p) - k)
    return best


def coding_potential(orf, cons_len):
    """Label (never a gate): 'noncoding' if the ORF spans <30% of the transcript."""
    frac = (orf * 3) / cons_len if cons_len else 0.0
    return "noncoding" if frac < 0.30 else "coding-capable/novel-protein"


def load_gff_from_lines(lines):
    """Build the biotype index from an iterable of GFF lines (test seam + file loader core)."""
    idx = defaultdict(list)
    for ln in lines:
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        bt = "unknown"
        for kv in f[8].split(";"):
            if kv.startswith("gene_biotype="):
                bt = kv.split("=", 1)[1]
        try:
            s, e = int(f[3]), int(f[4])
        except ValueError:
            continue
        idx[f[0]].append((s, e, bt))
    for c in idx:
        idx[c].sort()
    return idx


def load_gff_genes(path):
    """chrom -> sorted [(start, end, biotype)] for gene features; {} if path missing/unreadable."""
    if not path or not os.path.exists(path):
        return {}
    with open(path) as fh:
        return load_gff_from_lines(fh)


def biotype_of(idx, chrom, start, end):
    """Annotation LABEL of a gene overlapping [start,end]; never a gate.
    'unknown' if no annotation loaded, 'intergenic' if loaded but no overlap."""
    if not idx:
        return "unknown"
    ivs = idx.get(chrom, [])
    hits = [bt for (s, e, bt) in ivs if s <= end and e >= start]
    if not hits:
        return "intergenic"
    if any(b == "lncRNA" for b in hits):
        return "lncRNA"
    if any(b == "protein_coding" for b in hits):
        return "protein_coding-gene-body"
    if any(b == "pseudogene" for b in hits):
        return "pseudogene"
    return hits[0]


def read_fasta(path):
    """header (first whitespace-delimited token after '>') -> concatenated sequence."""
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


_TSV_COLS = ["cid", "chrom", "start", "end", "track", "genome_hit", "genome_id", "genome_cov",
             "divergence", "alt_cols", "alt_read_fraction", "alt_reads", "n_primary", "orf_aa",
             "coding_potential", "biotype", "sedef_partner", "protein", "copy_vs_allele",
             "status", "reason"]


def write_tsv(records, path):
    with open(path, "w") as fh:
        fh.write("\t".join(_TSV_COLS) + "\n")
        for r in records:
            fh.write("\t".join(str(r.get(c, "")) for c in _TSV_COLS) + "\n")


def build_record(cid, flag, hit, cons_len, orf, biotype, sedef_partner, protein, reason):
    """Assemble a promoted non-coding record. copy_vs_allele and status are set UNCONDITIONALLY."""
    tgt, pos, gid, gcov = hit
    return dict(
        cid=cid, chrom=flag["chrom"], start=flag["start"], end=flag["end"], track="noncoding",
        genome_hit=f"{tgt}:{pos}", genome_id=round(gid, 3), genome_cov=round(gcov, 2),
        divergence=round(100 * (1 - gid), 1),
        alt_cols=flag["n_alt_positions"], alt_read_fraction=round(flag["alt_read_fraction"], 3),
        alt_reads=flag["n_alt_reads"], n_primary=flag["n_primary_reads"],
        orf_aa=orf, coding_potential=coding_potential(orf, cons_len),
        biotype=biotype, sedef_partner=sedef_partner, protein=protein,
        copy_vs_allele=CANDIDATE, status=STATUS, reason=reason)
```

- [ ] **Step 4: Run to verify pass**

Run: `PYTHONHASHSEED=0 pytest bench/test_promote_noncoding.py -v`
Expected: PASS (all Task 1 + Task 2 tests).

- [ ] **Step 5: Commit**

```bash
git add bench/promote_noncoding.py bench/test_promote_noncoding.py
git commit -m "feat(noncoding): label helpers + build_record honesty rail + fasta/tsv io"
```

---

### Task 3: PAF re-score parser + `promote_all` orchestration + `main` I/O shell

**Files:**
- Modify: `bench/promote_noncoding.py`
- Test: `bench/test_promote_noncoding.py` (add tests)

**Interfaces:**
- Consumes: `classify_noncoding`, `build_record`, `orf_aa`, `biotype_of`, `sedef_partners` (imported).
- Produces:
  - `best_hits_from_paf(paf_text: str) -> dict[str, tuple[str,int,float,float,int]]` — cid → `(target, pos, genome_id, genome_cov, mapq)` for the best hit by `identity*coverage`. PAF columns: 0 qname, 1 qlen, 5 tname, 7 tstart, 9 nmatch, 10 blocklen, 11 mapq.
  - `promote_all(cons, flags, hits, gff, sed, carried, exclude=frozenset()) -> tuple[list[dict], Counter]` — pure orchestration; `own_locus = (target == flag['chrom']) and abs(pos - flag['start']) < 200000`. `exclude` = cids already promoted to the coding catalog (skipped, tallied `already-coding`), so the track surfaces only what the protein gate missed.
  - `main()` — thin I/O shell: union-load the flags JSONs, `read_fasta(cons.fa)`, load GFF + SEDEF + carried protein labels, run one `minimap2 -cx splice:hq --eqx -N1`, then `promote_all`, then write `gw_noncoding_copies.json` + `.tsv` + a tally log.

- [ ] **Step 1: Write the failing tests**

Append to `bench/test_promote_noncoding.py`:

```python
def test_best_hits_from_paf_picks_best():
    # cid q1: two hits; the higher identity*coverage wins (line B)
    paf = "\n".join([
        # qname qlen qs qe strand tname tlen ts te nmatch blocklen mapq ...tags
        "q1\t1000\t0\t900\t+\tNC_1\t9\t500\t1400\t800\t900\t60\ttp:A:P",
        "q1\t1000\t0\t950\t+\tNC_2\t9\t10\t960\t930\t950\t40\ttp:A:S",
    ])
    hits = PN.best_hits_from_paf(paf)
    assert set(hits) == {"q1"}
    tname, pos, gid, gcov, mapq = hits["q1"]
    assert tname == "NC_2"          # 0.979*0.95 > 0.889*0.90
    assert pos == 10 and mapq == 40
    assert round(gid, 3) == 0.979 and round(gcov, 3) == 0.95


def test_best_hits_ignores_short_lines():
    assert PN.best_hits_from_paf("garbage\tline\n\n") == {}


def test_promote_all_flagship_and_exclusions():
    cons = {"C_flag": "ATG" + "GCA" * 200, "C_art": "AAAA", "C_ref": "CCCC"}
    flags = {
        "C_flag": dict(chrom="NC_1", start=1000, end=2000, n_alt_positions=40,
                       alt_read_fraction=0.32, n_alt_reads=30, n_primary_reads=90),
        "C_art": dict(chrom="NC_1", start=5000, end=6000, n_alt_positions=40,
                      alt_read_fraction=0.32, n_alt_reads=30, n_primary_reads=90),
        "C_ref": dict(chrom="NC_1", start=8000, end=9000, n_alt_positions=40,
                      alt_read_fraction=0.32, n_alt_reads=30, n_primary_reads=90),
    }
    hits = {
        "C_flag": ("NC_1", 1050, 0.90, 1.00, 60),   # divergent, own-locus -> promote
        "C_art": ("NC_1", 5050, 0.30, 1.00, 60),    # genome_id < GID_LO -> artifact
        "C_ref": ("NC_1", 8050, 0.99, 1.00, 60),    # genome_id >= GID_HI -> ~REF
    }
    gff = PN.load_gff_from_lines(["NC_1\tX\tgene\t900\t2100\t.\t+\t.\tgene_biotype=lncRNA"])
    promoted, tally = PN.promote_all(cons, flags, hits, gff, {}, {})
    assert [r["cid"] for r in promoted] == ["C_flag"]
    assert promoted[0]["biotype"] == "lncRNA"
    assert promoted[0]["copy_vs_allele"] == "candidate-DNA-needed"
    assert promoted[0]["protein"] == "not-tested"     # carried default
    assert tally["artifact"] == 1 and tally["~REF"] == 1 and tally["noncoding-candidate"] == 1


def test_promote_all_missing_flag_or_hit_tallied():
    cons = {"C_x": "ACGT"}
    promoted, tally = PN.promote_all(cons, {}, {}, {}, {}, {})
    assert promoted == [] and tally["no-flag"] == 1
    promoted, tally = PN.promote_all(cons, {"C_x": dict(chrom="NC_1", start=1, end=2,
        n_alt_positions=40, alt_read_fraction=0.3, n_alt_reads=30, n_primary_reads=90)}, {}, {}, {}, {})
    assert promoted == [] and tally["no-hit"] == 1


def test_promote_all_excludes_already_coding():
    cons = {"C_flag": "ATG" + "GCA" * 200}
    flags = {"C_flag": dict(chrom="NC_1", start=1000, end=2000, n_alt_positions=40,
                            alt_read_fraction=0.32, n_alt_reads=30, n_primary_reads=90)}
    hits = {"C_flag": ("NC_1", 1050, 0.90, 1.00, 60)}   # would promote if not excluded
    promoted, tally = PN.promote_all(cons, flags, hits, {}, {}, {}, exclude={"C_flag"})
    assert promoted == [] and tally["already-coding"] == 1
```

- [ ] **Step 2: Run to verify failure**

Run: `PYTHONHASHSEED=0 pytest bench/test_promote_noncoding.py -v -k "best_hits or promote_all"`
Expected: FAIL — `AttributeError: module 'promote_noncoding' has no attribute 'best_hits_from_paf'`.

- [ ] **Step 3: Implement the parser, orchestration, and I/O shell**

Append to `bench/promote_noncoding.py`:

```python
def best_hits_from_paf(paf_text):
    """cid -> (target, pos, genome_id, genome_cov, mapq) for the best hit by identity*coverage."""
    best = {}
    for line in paf_text.splitlines():
        c = line.split("\t")
        if len(c) < 12:
            continue
        try:
            cid, qlen, tname, tstart = c[0], int(c[1]), c[5], int(c[7])
            nm, al, mapq = int(c[9]), int(c[10]), int(c[11])
        except (ValueError, IndexError):
            continue
        gid = nm / al if al else 0.0
        gcov = al / qlen if qlen else 0.0
        sc = gid * gcov
        if cid not in best or sc > best[cid][5]:
            best[cid] = (tname, tstart, gid, gcov, mapq, sc)
    return {cid: (t, p, gid, gcov, mapq) for cid, (t, p, gid, gcov, mapq, _sc) in best.items()}


def promote_all(cons, flags, hits, gff, sed, carried, exclude=frozenset()):
    """Pure orchestration: classify every consensus, assemble promoted records, tally calls.
    `exclude` = cids already promoted to the coding catalog; they are skipped (tallied
    'already-coding'), so the non-coding track surfaces only what the protein gate MISSED
    rather than duplicating the coding catalog."""
    promoted, tally = [], Counter()
    for cid, seq in cons.items():
        if cid in exclude:
            tally["already-coding"] += 1
            continue
        flag = flags.get(cid)
        if flag is None:
            tally["no-flag"] += 1
            continue
        h = hits.get(cid)
        if h is None:
            tally["no-hit"] += 1
            continue
        tgt, pos, gid, gcov, _mapq = h
        own = (tgt == flag["chrom"]) and abs(pos - flag["start"]) < 200000
        rec = dict(genome_id=gid, genome_cov=gcov, own_locus=own,
                   alt_cols=flag["n_alt_positions"], alt_read_fraction=flag["alt_read_fraction"],
                   alt_reads=flag["n_alt_reads"], n_primary=flag["n_primary_reads"])
        promote, call, reason = classify_noncoding(rec)
        tally[call] += 1
        if not promote:
            continue
        parts = sedef_partners(sed, flag["chrom"], flag["start"], flag["end"]) if sed else []
        sp = f"{parts[0][0]}:{parts[0][1]}" if parts else None
        promoted.append(build_record(
            cid, flag, (tgt, pos, gid, gcov), len(seq), orf_aa(seq),
            biotype_of(gff, flag["chrom"], flag["start"], flag["end"]),
            sp, carried.get(cid, "not-tested"), reason))
    promoted.sort(key=lambda r: r["genome_id"])
    return promoted, tally


def main():
    import argparse
    ap = argparse.ArgumentParser(description="Non-coding-aware reference-absent promotion.")
    ap.add_argument("--cons", default=f"{OUT}/cons.fa")
    ap.add_argument("--out", default=f"{OUT}/gw_noncoding_copies.json")
    ap.add_argument("--threads", default="6")
    a = ap.parse_args()

    # union-load both flag files so every cons.fa cid finds its flag record
    flags = {}
    for fn in [f"{CAT}/genomewide_flags.json", f"{CAT}/genomewide_flags_new.json"]:
        if os.path.exists(fn):
            for f in json.load(open(fn)):
                flags[f"{f['chrom']}_{f['start']}"] = f
    cons = read_fasta(a.cons)
    gff = load_gff_genes(GFF)
    sed = load_sedef(SEDEF) if SEDEF and os.path.exists(SEDEF) else {}
    carried = {}
    rap = f"{OUT}/gw_reference_absent_copies.json"
    if os.path.exists(rap):
        carried = {r["cid"]: r.get("protein") for r in json.load(open(rap))}
    exclude = set()
    disc = f"{OUT}/gw_discriminated.json"
    if os.path.exists(disc):
        exclude = {r["cid"] for r in json.load(open(disc))}   # already-promoted coding copies

    paf = subprocess.run([MM2, "-cx", "splice:hq", "--eqx", "-N1", "-t", a.threads, FASTA, a.cons],
                         capture_output=True, text=True, timeout=3600).stdout
    hits = best_hits_from_paf(paf)

    promoted, tally = promote_all(cons, flags, hits, gff, sed, carried, exclude)
    json.dump(promoted, open(a.out, "w"), indent=1)
    write_tsv(promoted, a.out[:-5] + ".tsv" if a.out.endswith(".json") else a.out + ".tsv")

    bt = Counter(r["biotype"] for r in promoted)
    print(f"non-coding promotion: {len(promoted)} flagged reference-divergent candidates "
          f"(track=noncoding; all copy_vs_allele={CANDIDATE})")
    print(f"  by biotype: {dict(bt)}")
    print("  call tally (incl. exclusions, logged not dropped):")
    for k, v in tally.most_common():
        print(f"    {k:18s} {v}")
    print(f"  wrote {a.out} (+ .tsv)")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run to verify pass**

Run: `PYTHONHASHSEED=0 pytest bench/test_promote_noncoding.py -v`
Expected: PASS (all unit tests through Task 3).

- [ ] **Step 5: Commit**

```bash
git add bench/promote_noncoding.py bench/test_promote_noncoding.py
git commit -m "feat(noncoding): PAF re-score parser + promote_all + main I/O shell"
```

---

### Task 4: Data-gated integration test on the real 734 consensuses

**Files:**
- Test: `bench/test_promote_noncoding.py` (add one data-gated test)

**Interfaces:**
- Consumes: `main` (Task 3) and the on-disk `cons.fa`, `GGO.fasta`, flag JSONs. Runs one real `minimap2` (single foreground process; not `copy_assign`; WSL2-safe).

- [ ] **Step 1: Write the data-gated integration test**

Append to `bench/test_promote_noncoding.py`:

```python
import json as _json

_CONS = "/home/juanfra/winloci_scratch/refabsent/gw_promoted/cons.fa"
_GENOME = "/home/juanfra/winloci_scratch/GGO.fasta"
_DATA_PRESENT = os.path.exists(_CONS) and os.path.exists(_GENOME) and os.path.exists(PN.MM2)

_STRONG_CREDIBLE = [
    "NC_073236.2_139051025",  # flagship lncRNA
    "NC_073236.2_44341131",   # lncRNA
    "NC_073234.2_27052843",   # ETV6 gene-body
    "NC_073230.2_167930813",  # VIPR2 gene-body
    "NC_073231.2_4861523",    # TDRP gene-body
    "NC_073242.2_18569043",   # intergenic, 1028-aa novel ORF
    "NC_073224.2_45327469",   # intergenic -> TRAF5
]
_ARTIFACTS = ["NC_073243.2_30120917", "NC_073231.2_39477379"]  # div ~70/74, genome_id ~0.29/0.26


import pytest


@pytest.mark.skipif(not _DATA_PRESENT, reason="cons.fa / GGO.fasta / minimap2 not present")
def test_integration_real_734(tmp_path):
    out = tmp_path / "gw_noncoding_copies.json"
    sys.argv = ["promote_noncoding", "--cons", _CONS, "--out", str(out), "--threads", "6"]
    PN.main()
    recs = _json.load(open(out))
    cids = {r["cid"] for r in recs}
    # flagship + all 7 strong-credible are promoted
    missing = [c for c in _STRONG_CREDIBLE if c not in cids]
    assert not missing, f"strong-credible missing from promotion: {missing}"
    # the 2 high-divergence repeat/chimera artifacts are excluded
    assert not (set(_ARTIFACTS) & cids), "artifact leaked into promotion"
    # sane yield band (looser bar than the workflow's ~19; report the count)
    print(f"\n[integration] promoted {len(recs)} non-coding candidates; "
          f"biotypes={sorted({r['biotype'] for r in recs})}")
    assert len(recs) >= 15
    # honesty rail holds on every record
    assert all(r["copy_vs_allele"] == "candidate-DNA-needed" for r in recs)
    assert all(r["status"] == "flagged-reference-divergent-candidate" for r in recs)
    assert all(r["track"] == "noncoding" for r in recs)
    # flagship specifics
    fs = next(r for r in recs if r["cid"] == "NC_073236.2_139051025")
    assert fs["coding_potential"] == "noncoding"
    if os.path.exists(PN.GFF):            # biotype needs the (linuxdisk-mounted) annotation
        assert fs["biotype"] == "lncRNA"
```

- [ ] **Step 2: Run the integration test (foreground, single minimap2 — obeys WSL2 crash rule)**

Run: `PYTHONHASHSEED=0 pytest bench/test_promote_noncoding.py::test_integration_real_734 -v -s`
Expected: PASS. The `-s` prints the actual promoted count + biotype set. If `cons.fa`/genome absent, the test SKIPS (still green).

If the test fails because a strong-credible cid is missing:
- Inspect its re-scored hit: `grep '<cid>' <the PAF>` — check whether `genome_id`/`genome_cov` landed just outside a band, and confirm the flag record exists in the union (`no-flag` in the tally means a cid-key mismatch → verify the flag files cover it). Do not loosen a constant to force a single cid without confirming the signal is genuinely credible.

- [ ] **Step 3: Run the full suite once more**

Run: `PYTHONHASHSEED=0 pytest bench/test_promote_noncoding.py -v`
Expected: PASS (all unit tests + the integration test, or the integration SKIPS off-box).

- [ ] **Step 4: Commit**

```bash
git add bench/test_promote_noncoding.py
git commit -m "test(noncoding): data-gated integration on the real 734 consensuses"
```

---

## Notes for the implementer

- **Run every `pytest` from the repo root** (`/mnt/c/Users/jfris/Desktop/Rustle`) so the module's `sys.path.insert(0, "bench")` + `from copy_vs_allele_structural import ...` resolves.
- **Do not import `hidden_copy_scan` or `promote_genomewide`** into `promote_noncoding.py` — they pull `pysam`/`pyabpoa` and would break hermetic unit tests. `MIN_DEPTH` and `orf_aa` are intentionally local (see Global Constraints).
- **WSL2 crash rule:** the only heavy step is the single `minimap2` in `main` (734 small queries vs the genome). It is one foreground subprocess — never `copy_assign`, never backgrounded/nohup/pkill. Outputs live under `winloci_scratch/refabsent/gw_promoted/`.
- **Byte-identical guarantee:** the script only ever writes `gw_noncoding_copies.json`/`.tsv`; it must never write or touch `gw_reference_absent_copies.json`, `gw_discriminated.json`, `cons.fa`, or any coding-catalog file.
```
