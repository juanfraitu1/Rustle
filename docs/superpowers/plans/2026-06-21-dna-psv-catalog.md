# DNA-derived PSV Identifiability Catalog (Phase 1) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A reference-only (T2T DNA) catalog of exonic PSVs per multi-copy-family copy pair → a genome-wide identifiability table (resolvable K≥1 vs genuine-K=0), cross-checked against the RNA census.

**Architecture:** One self-contained Python bench script `bench/dna_psv_catalog.py`. Per family: extract copy intervals from the T2T reference, align each copy to a family reference copy with `minimap2 asm20 --cs`, project to an allele matrix, restrict PSV columns to exons, compute per-pair K. Then a cross-check against `bench/copy_resolution_census.py`'s `classify_pair`, and a deterministic `check()`.

**Tech Stack:** Python 3, pysam (FASTA fetch + tags), minimap2 2.30 (subprocess), GFF parsing (plain text/awk), markdown summary.

## Global Constraints

- Python bench tooling only; no Rust / pipeline integration; reference-only (no read re-alignment in Phase 1).
- Substrate (absolute paths; scripts must not rely on cwd):
  - Families `/mnt/c/Users/jfris/Desktop/Rustle/bench/denovo_families.tsv` (header `family_id\tn_copies\tmembers`; members comma-separated `DN_<chrom>_<start>_<n>`; chrom has dots, e.g. `NC_073247.2`).
  - Reference `/home/juanfra/winloci_scratch/GGO.fasta` (+ `.fai`).
  - Annotation `/home/juanfra/winloci_scratch/GGO_genomic.gff`.
  - Scratch (per-family matrices) `/home/juanfra/winloci_scratch/dna_catalog/`.
- Copy interval = annotated gene span (GFF `gene` overlapping the locus start) if present, else `[max(0,start-2000), start+40000]`.
- Alignment: `minimap2 -cx asm20 --cs=long ref0.fa copy_i.fa` (paf + cs). PSV = substitution column where ≥2 copies differ; K = exonic PSV Hamming per pair. Exonic = ref0 position inside a GFF exon overlapping the ref0 window; no overlapping exon ⇒ `psv_exonic = NA`, verdict `unannotated`.
- Output columns (exact order): `family, copyA, copyB, chromA, chromB, co_located, aln_identity, psv_total, psv_exonic, private_exon_bp, verdict`.
- Verdict ∈ {`resolvable`, `genuine_k0`, `unannotated`}; co_located = same chrom within 2,000,000 bp.
- Mega-family cap: process at most 200 members per family (log dropped count; do not silently truncate).
- Commit trailers (every commit ends with):
  ```
  Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_0193po8vtwu5dbDNComynKAa
  ```
- Work on the current branch; add only the new files this plan names; do not touch unrelated modified files.

---

## Task 1: Copy intervals, per-family alignment, and the allele matrix

**Files:**
- Create: `bench/dna_psv_catalog.py`

**Interfaces:**
- Produces: `parse_member(m) -> (chrom,start)|None`; `gene_span(chrom,start) -> (s,e)|None`; `copy_interval(chrom,start) -> (s,e)`; `extract_seq(chrom,s,e) -> str`; `align_to_ref(ref_seq, qry_seq) -> list[(ref_pos, qry_base)]` (substitution columns, 0-based ref-local); `family_allele_matrix(members) -> (ref0, cols, matrix, meta)` where `matrix[copy_label][ref0_localpos] = base`.

- [ ] **Step 1: Write the alignment + matrix smoke test (real data, MAGEA copies)**

Create `bench/dna_psv_catalog.py` with the imports + helpers below, then a `__main__` smoke that aligns two known MAGEA copies and prints substitution columns. Use two copies that are NOT identical so the test is meaningful — `MAGd1` pair (161266095 vs 161413981) is expected to carry substitutions:

```python
"""DNA-derived PSV identifiability catalog (Phase 1, reference-only). See
docs/superpowers/specs/2026-06-21-dna-psv-catalog-design.md.

Per family: extract each copy's T2T interval, align every copy to the longest (ref0) with
minimap2 asm20 --cs, project substitutions onto ref0 coordinates -> allele matrix -> exonic PSV columns ->
per-pair K (Hamming over exonic PSVs). Cross-check vs the RNA census. Deterministic check() at the bottom."""
import collections, os, re, subprocess, tempfile
import pysam

FAM_TSV  = "/mnt/c/Users/jfris/Desktop/Rustle/bench/denovo_families.tsv"
FASTA    = "/home/juanfra/winloci_scratch/GGO.fasta"
GFF      = "/home/juanfra/winloci_scratch/GGO_genomic.gff"
SCRATCH  = "/home/juanfra/winloci_scratch/dna_catalog"
OUT_TSV  = "/mnt/c/Users/jfris/Desktop/Rustle/bench/dna_psv_catalog.tsv"
OUT_MD   = "/mnt/c/Users/jfris/Desktop/Rustle/bench/dna_psv_catalog_summary.md"
WIN, COLOC_WIN, MAX_MEMBERS = 40_000, 2_000_000, 200

_fa = pysam.FastaFile(FASTA)

def parse_member(m):
    mm = re.match(r"DN_(N[CW]_\d+\.\d+)_(\d+)_(\d+)$", m)
    return (mm.group(1), int(mm.group(2))) if mm else None

def gene_span(chrom, start):
    """Annotated gene span overlapping `start` (first hit), via awk over the GFF."""
    r = subprocess.run(["awk", "-F", "\t",
        f'$1=="{chrom}" && $3=="gene" && $4<={start} && $5>={start} {{print $4"\\t"$5; exit}}', GFF],
        capture_output=True, text=True)
    if r.stdout.strip():
        s, e = r.stdout.split()[:2]; return (int(s) - 1, int(e))   # GFF 1-based inclusive -> 0-based half-open
    return None

def copy_interval(chrom, start):
    g = gene_span(chrom, start)
    return g if g else (max(0, start - 2000), start + WIN)

def extract_seq(chrom, s, e):
    return _fa.fetch(chrom, s, e).upper()

def align_to_ref(ref_seq, qry_seq):
    """minimap2 asm20 --cs=long: return substitution columns as (ref_localpos, qry_base), 0-based on ref_seq."""
    with tempfile.TemporaryDirectory() as td:
        rp, qp = os.path.join(td, "r.fa"), os.path.join(td, "q.fa")
        open(rp, "w").write(">r\n" + ref_seq + "\n"); open(qp, "w").write(">q\n" + qry_seq + "\n")
        out = subprocess.run(["minimap2", "-cx", "asm20", "--cs=long", "-t", "1", rp, qp],
                             capture_output=True, text=True).stdout
    subs = []
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 9: continue
        rstart = int(f[7])                          # target (ref) start, 0-based
        cs = next((x[5:] for x in f if x.startswith("cs:Z:")), None)
        if cs is None: continue
        rpos = rstart
        for op, val in re.findall(r"([:=*+\-])([A-Za-z0-9]+)", cs):
            if op == ":":            rpos += int(val)                 # identical run
            elif op == "=":          rpos += len(val)                 # identical (long form)
            elif op == "*":          subs.append((rpos, val[1].upper())); rpos += 1   # *<ref><qry> substitution
            elif op == "-":          rpos += len(val)                 # deletion from query (consume ref)
            elif op == "+":          pass                             # insertion in query (no ref consume)
    return subs

def family_allele_matrix(members):
    """members = [(label, chrom, start)]. Align all copies to the longest (ref0); return
    (ref0_label, ref0_chrom, ref0_iv, matrix, intervals) with matrix[label][ref0_localpos]=base."""
    ivs = {}
    for label, chrom, start in members:
        s, e = copy_interval(chrom, start); ivs[label] = (chrom, s, e)
    ref0 = max(ivs, key=lambda L: ivs[L][2] - ivs[L][1])
    rc, rs, re_ = ivs[ref0]; ref_seq = extract_seq(rc, rs, re_)
    matrix = collections.defaultdict(dict)        # label -> {ref0_localpos: base}
    for label, (c, s, e) in ivs.items():
        if label == ref0: continue
        subs = align_to_ref(ref_seq, extract_seq(c, s, e))
        for pos, base in subs:
            matrix[label][pos] = base
    return ref0, rc, (rs, re_), matrix, ivs

if __name__ == "__main__":
    # smoke: two MAGEA copies expected to differ (MAGd1 pair)
    members = [("A", "NC_073247.2", 161266095), ("B", "NC_073247.2", 161413981)]
    ref0, rc, riv, matrix, ivs = family_allele_matrix(members)
    nonref = [L for L in ivs if L != ref0][0]
    print(f"ref0={ref0} {rc}:{riv[0]}-{riv[1]} ({riv[1]-riv[0]}bp); "
          f"substitution columns {nonref} vs ref0: {len(matrix[nonref])}")
```

- [ ] **Step 2: Run the smoke**

Run: `python3 bench/dna_psv_catalog.py`
Expected: prints a `ref0=...` line with a **nonzero** substitution-column count (these MAGEA copies are homologous-but-not-identical; if it prints 0, the `cs` parse or interval extraction is wrong — debug before proceeding). Confirm `minimap2`/`pysam` import cleanly.

- [ ] **Step 3: Commit**

```bash
git add bench/dna_psv_catalog.py
git commit -m "$(printf '%s\n' "bench(dna-catalog): per-family copy alignment + allele matrix (minimap2 asm20 --cs)" "" "Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>" "Claude-Session: https://claude.ai/code/session_0193po8vtwu5dbDNComynKAa")"
```

---

## Task 2: Exonic restriction, per-pair K, genome-wide TSV, matrices

**Files:**
- Modify: `bench/dna_psv_catalog.py`
- Create (outputs): `bench/dna_psv_catalog.tsv`, `/home/juanfra/winloci_scratch/dna_catalog/<family>.json`

**Interfaces:**
- Consumes: `family_allele_matrix`, `copy_interval` (Task 1).
- Produces: `exon_intervals(chrom, s, e) -> list[(s,e)]`; `pair_K(matrix, exonic_set, a, b) -> (psv_total, psv_exonic|None)`; `build_catalog() -> rows` (writes TSV + matrices).

- [ ] **Step 1: Write the exonic + per-pair logic and a panel-anchored test**

Add to `bench/dna_psv_catalog.py` (before `__main__`). `exon_intervals` pulls GFF exons overlapping the ref0 window; a PSV column is exonic iff its ref0 *genomic* position is inside an exon:

```python
import json

def exon_intervals(chrom, s, e):
    """GFF exon intervals (0-based half-open) overlapping [s,e]."""
    r = subprocess.run(["awk", "-F", "\t",
        f'$1=="{chrom}" && $3=="exon" && $4<{e} && $5>{s} {{print $4"\\t"$5}}', GFF],
        capture_output=True, text=True)
    out = []
    for line in r.stdout.splitlines():
        a, b = line.split()[:2]; out.append((int(a) - 1, int(b)))
    return out

def _exonic_localpos(rc, riv, exons):
    """set of ref0-LOCAL positions (offset from riv[0]) that fall inside any exon; None if no exon annotation."""
    if not exons: return None
    s0 = riv[0]; ex = set()
    for a, b in exons:
        for p in range(max(a, riv[0]), min(b, riv[1])):
            ex.add(p - s0)
    return ex

def pair_K(matrix, ref_label, exonic, a, b):
    """psv_total, psv_exonic(None if exonic is None). Allele at ref0 positions: ref label = ref base (implicit:
    a position is a substitution col iff some copy differs from ref0). Build the column set, then Hamming."""
    cols = set(matrix.get(a, {})) | set(matrix.get(b, {}))
    if a == ref_label: cols = set(matrix.get(b, {}))
    if b == ref_label: cols = set(matrix.get(a, {}))
    def allele(label, p):
        if label == ref_label: return "REF"          # ref0 carries the reference base everywhere
        return matrix.get(label, {}).get(p, "REF")    # absent substitution => matches ref0
    total = sum(1 for p in cols if allele(a, p) != allele(b, p))
    if exonic is None: return total, None
    exonic_cnt = sum(1 for p in cols if p in exonic and allele(a, p) != allele(b, p))
    return total, exonic_cnt
```

Then a test in `__main__` (replace the Task-1 smoke) over a small panel — the MAGEA de-novo pairs (exon-identical ⇒ exonic K=0) plus a resolvable cross-chrom copy (AK6, K≥1):

```python
PANEL_FAMILIES = {
    # label -> (chrom, start); grouped families for the check
    "MAGEA_dn": [("d0a","NC_073247.2",161251228),("d0b","NC_073247.2",161458538),
                 ("d2a","NC_073247.2",164381222),("d2b","NC_073247.2",164442447),
                 ("d3a","NC_073247.2",164397061),("d3b","NC_073247.2",164426194)],
    "AK6":      [("AK6","NC_073243.2",40875017),("LOC","NC_073227.2",47415194)],
}

def _panel_pair_K(members, a_lbl, b_lbl):
    ref0, rc, riv, matrix, ivs = family_allele_matrix([(l, c, s) for l, c, s in members])
    exonic = _exonic_localpos(rc, riv, exon_intervals(rc, riv[0], riv[1]))
    return pair_K(matrix, ref0, exonic, a_lbl, b_lbl), ref0
```

- [ ] **Step 2: Write `build_catalog()` (genome-wide) + run it**

Add `build_catalog()` that iterates all families, builds the matrix, computes every co-located (and cross-chrom) pair's K, writes the TSV + per-family matrix json:

```python
def build_catalog():
    os.makedirs(SCRATCH, exist_ok=True)
    cols = ["family","copyA","copyB","chromA","chromB","co_located","aln_identity",
            "psv_total","psv_exonic","private_exon_bp","verdict"]
    rows = []
    for line in open(FAM_TSV):
        if line.startswith("family_id"): continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 3: continue
        fid, members = p[0], [x for x in (parse_member(m) for m in p[2].split(",")) if x]
        if len(members) < 2: continue
        labeled = [(f"L{i}", c, s) for i, (c, s) in enumerate(members[:MAX_MEMBERS])]
        if len(members) > MAX_MEMBERS:
            print(f"  [cap] {fid}: {len(members)} members -> processed {MAX_MEMBERS}", flush=True)
        try:
            ref0, rc, riv, matrix, ivs = family_allele_matrix(labeled)
        except Exception as ex:
            print(f"  [skip] {fid}: {ex}", flush=True); continue
        exonic = _exonic_localpos(rc, riv, exon_intervals(rc, riv[0], riv[1]))
        json.dump({"ref0": ref0, "chrom": rc, "iv": riv,
                   "matrix": {k: v for k, v in matrix.items()},
                   "exonic_n": (len(exonic) if exonic else 0)},
                  open(os.path.join(SCRATCH, f"{fid}.json"), "w"))
        labels = [l for l, _, _ in labeled]
        for i in range(len(labels)):
            for j in range(i + 1, len(labels)):
                a, b = labels[i], labels[j]
                ca, cb = ivs[a][0], ivs[b][0]
                coloc = (ca == cb and abs(ivs[a][1] - ivs[b][1]) <= COLOC_WIN)
                total, exo = pair_K(matrix, ref0, exonic, a, b)
                verdict = "unannotated" if exo is None else ("genuine_k0" if exo == 0 else "resolvable")
                rows.append(dict(family=fid, copyA=f"{ca}:{ivs[a][1]}", copyB=f"{cb}:{ivs[b][1]}",
                                 chromA=ca, chromB=cb, co_located=int(coloc), aln_identity="",
                                 psv_total=total, psv_exonic=("NA" if exo is None else exo),
                                 private_exon_bp=0, verdict=verdict))
    with open(OUT_TSV, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    return rows
```
(Restrict pairs to co-located OR cross-chrom-within-family; `aln_identity`/`private_exon_bp` may be left as placeholders `""`/`0` in Phase 1 — they are reported columns, not used by the check. Keep them in the header for Phase 2.)

Run: `python3 -c "import bench.dna_psv_catalog as d; d.build_catalog()" 2>&1 | tail -5`
Expected: completes (minutes); writes `bench/dna_psv_catalog.tsv` with co-located pairs and `verdict` populated.

- [ ] **Step 3: Confirm the panel anchors by hand**

Run: `python3 -c "import bench.dna_psv_catalog as d; print('d0',_:=d._panel_pair_K(d.PANEL_FAMILIES['MAGEA_dn'],'d0a','d0b')); print('d2',d._panel_pair_K(d.PANEL_FAMILIES['MAGEA_dn'],'d2a','d2b')); print('AK6',d._panel_pair_K(d.PANEL_FAMILIES['AK6'],'AK6','LOC'))"`
Expected: MAGEA `d0`/`d2` exonic K = 0 (exon-identical — matches the per-read finding that pair0's divergence is purely intronic and pair2/3 are identical); AK6 exonic K ≥ 1. If `d0`/`d2` exonic K > 0, investigate (the exonic restriction or the cs-substitution parse) before Task 3 — do not proceed on a contradiction.

- [ ] **Step 4: Commit**

```bash
git add bench/dna_psv_catalog.py bench/dna_psv_catalog.tsv
git commit -m "$(printf '%s\n' "bench(dna-catalog): exonic PSV restriction + per-pair K + genome-wide TSV" "" "Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>" "Claude-Session: https://claude.ai/code/session_0193po8vtwu5dbDNComynKAa")"
```

---

## Task 3: RNA cross-check, summary, and the deterministic `check()`

**Files:**
- Modify: `bench/dna_psv_catalog.py`
- Create (output): `bench/dna_psv_catalog_summary.md`

**Interfaces:**
- Consumes: `build_catalog` (Task 2); `copy_resolution_census.classify_pair` (existing).
- Produces: `cross_check(rows) -> dict`; `check()` (run as `python3 bench/dna_psv_catalog.py`).

- [ ] **Step 1: Write the cross-check + summary + check, with panel + concordance asserts**

Add to `bench/dna_psv_catalog.py`. The cross-check joins co-located catalog pairs to the RNA census via `classify_pair` (which returns `frac_same` for the same loci):

```python
import sys
sys.path.insert(0, os.path.dirname(__file__))

def cross_check(rows):
    """DNA-K=0 vs RNA-K0 (frac_same>=0.95) on co-located pairs the census can classify. Returns confusion + lists."""
    import copy_resolution_census as census
    conf = collections.Counter()         # (dna_k0, rna_k0) -> count
    discord = {"dna_k_rna_tied": [], "dna_k0_rna_resolv": []}
    n = 0
    for r in rows:
        if not r["co_located"] or r["psv_exonic"] == "NA": continue
        ca = r["chromA"]; sa = int(r["copyA"].split(":")[1]); sb = int(r["copyB"].split(":")[1])
        try:
            cp = census.classify_pair(ca, min(sa, sb), max(sa, sb))
        except Exception:
            continue
        if cp is None or cp.get("n_xmap", 0) < 3 or "frac_same" not in cp: continue
        dna_k0 = (int(r["psv_exonic"]) == 0)
        rna_k0 = (cp["frac_same"] >= 0.95)
        conf[(dna_k0, rna_k0)] += 1; n += 1
        if dna_k0 != rna_k0:
            (discord["dna_k_rna_tied"] if (not dna_k0 and rna_k0) else discord["dna_k0_rna_resolv"]).append(
                dict(family=r["family"], pair=f"{r['copyA']}~{r['copyB']}",
                     dna_K=r["psv_exonic"], rna_frac_same=round(cp["frac_same"], 3)))
    concordant = conf[(True, True)] + conf[(False, False)]
    return dict(n=n, concordant=concordant, conf=dict(conf), discord=discord,
                concordance=(concordant / n if n else 0.0))

def write_summary(rows, xc):
    coloc = [r for r in rows if r["co_located"] and r["psv_exonic"] != "NA"]
    res = sum(1 for r in coloc if r["verdict"] == "resolvable")
    k0  = sum(1 for r in coloc if r["verdict"] == "genuine_k0")
    na  = sum(1 for r in rows if r["psv_exonic"] == "NA")
    L = ["# DNA-derived PSV identifiability catalog (Phase 1)\n",
         f"- co-located classified pairs: **{len(coloc)}** -> resolvable **{res}** "
         f"({100*res/max(1,len(coloc)):.0f}%), genuine-K=0 **{k0}** ({100*k0/max(1,len(coloc)):.0f}%)",
         f"- pairs with no exon annotation (NA): {na}",
         f"- **cross-check DNA-K=0 vs RNA-K0** on {xc['n']} census-classified pairs: "
         f"concordance **{100*xc['concordance']:.0f}%** (confusion {xc['conf']})",
         f"- discordant DNA-K≥1 ∧ RNA-tied (expressible-but-not-expressed): {len(xc['discord']['dna_k_rna_tied'])}",
         f"- discordant DNA-K=0 ∧ RNA-resolvable (pseudo-K=0 / indel): {len(xc['discord']['dna_k0_rna_resolv'])}\n"]
    open(OUT_MD, "w").write("\n".join(L) + "\n")

def check():
    # panel anchors
    d0 = _panel_pair_K(PANEL_FAMILIES["MAGEA_dn"], "d0a", "d0b")[0]
    d2 = _panel_pair_K(PANEL_FAMILIES["MAGEA_dn"], "d2a", "d2b")[0]
    ak = _panel_pair_K(PANEL_FAMILIES["AK6"], "AK6", "LOC")[0]
    print(f"    panel: MAGEA d0 exonicK={d0[1]} d2 exonicK={d2[1]}  AK6 exonicK={ak[1]}")
    assert d0[1] == 0, f"MAGEA d0 should be exon-identical (K=0), got {d0[1]}"
    assert d2[1] == 0, f"MAGEA d2 should be exon-identical (K=0), got {d2[1]}"
    assert ak[1] is not None and ak[1] >= 1, f"AK6 should be resolvable (K>=1), got {ak[1]}"
    rows = build_catalog()
    for r in rows:
        if r["psv_exonic"] != "NA":
            assert int(r["psv_exonic"]) <= r["psv_total"], f"exonic>total in {r['family']}"
    xc = cross_check(rows); write_summary(rows, xc)
    print(f"    cross-check: n={xc['n']} concordance={xc['concordance']:.2f} conf={xc['conf']}")
    assert xc["n"] >= 50, f"too few census-classified co-located pairs to cross-check: {xc['n']}"
    assert xc["concordance"] >= 0.80, f"DNA/RNA identifiability concordance too low: {xc['concordance']:.2f}"
    print(f"OK  - dna-catalog: {len(rows)} pairs; "
          f"co-located resolvable/K0 split + {100*xc['concordance']:.0f}% DNA/RNA concordance")
    return rows

if __name__ == "__main__":
    check()
```

- [ ] **Step 2: Run the full check**

Run: `python3 bench/dna_psv_catalog.py 2>&1 | tail -6`
Expected: the `panel:` line (MAGEA d0/d2 exonicK=0, AK6 exonicK≥1), the `cross-check:` line, then `OK  - dna-catalog: <N> pairs; ... ~XX% DNA/RNA concordance`, exit 0. Load-bearing asserts: panel anchors, concordance ≥ 0.80, ≥50 cross-checked pairs.

- [ ] **Step 3: If concordance < 0.80, diagnose (do NOT weaken the assert)**

The cross-check is the scientific claim. If concordance is low: (a) confirm `classify_pair` is being called with the same loci the catalog used (the census infers spans from reads — coordinates should still match on start); (b) inspect `xc['discord']` — a *systematic* `DNA-K≥1 ∧ RNA-tied` skew is the real "expressible-but-not-expressed" signal (PSV in an unexpressed exon), not a bug, and is expected to be a minority; a *systematic* `DNA-K=0 ∧ RNA-resolvable` skew means substitution-only PSVs miss indel-driven divergence (note it; `private_exon_bp` is the Phase-2 fix). Only if the discordance is an artifact (coordinate mismatch, cs-parse bug) fix the code; if it is the genuine biology, report it in the summary and bring the number to the reviewer rather than forcing 0.80.

- [ ] **Step 4: Confirm outputs**

Run: `head -1 bench/dna_psv_catalog.tsv; grep -c "genuine_k0" bench/dna_psv_catalog.tsv; sed -n '1,8p' bench/dna_psv_catalog_summary.md`
Expected: the 11-column header; a nonzero genuine-K=0 count; the summary with the resolvable/K=0 split and the concordance line.

- [ ] **Step 5: Commit**

```bash
git add bench/dna_psv_catalog.py bench/dna_psv_catalog.tsv bench/dna_psv_catalog_summary.md
git commit -m "$(printf '%s\n' "bench(dna-catalog): RNA cross-check + summary + deterministic check (Phase 1 complete)" "" "Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>" "Claude-Session: https://claude.ai/code/session_0193po8vtwu5dbDNComynKAa")"
```

---

## Self-Review

**Spec coverage:** copy intervals + per-family `minimap2 asm20 --cs` alignment + allele matrix → Task 1; exonic PSV restriction + per-pair K + genome-wide TSV + persisted matrices → Task 2; RNA cross-check + summary + deterministic check → Task 3. Non-goal (Phase 2 decoder) excluded; matrices persisted for it. All spec sections covered.

**Placeholder scan:** the only intentionally-deferred fields are `aln_identity` (`""`) and `private_exon_bp` (`0`) — explicitly Phase-2 columns kept in the header, not used by any assert; called out in Task 2 Step 2. No TODO/"handle appropriately"/missing code. Every step has runnable code + expected output.

**Type/name consistency:** `family_allele_matrix` returns `(ref0, rc, riv, matrix, ivs)` — consumed with that exact unpacking in `_panel_pair_K` and `build_catalog`. `pair_K(matrix, ref_label, exonic, a, b) -> (total, exo|None)` — `exo is None` ⇔ `unannotated`; used identically in `build_catalog` (verdict) and `check` (panel via `_panel_pair_K`, which returns `(total,exo)` as `d0[1]`=exo). `classify_pair` reused from `copy_resolution_census` returns `frac_same` (matches the census's contract). TSV columns identical between the `cols` list and the `write` loop. `co_located` int(bool) consistent between `build_catalog` and `cross_check`.

**Risk note:** runtime is one minimap2 call per copy (~ total copies, bounded); the MAX_MEMBERS=200 cap is logged. The concordance assert (≥0.80) encodes the headline claim and Task 3 Step 3 explicitly forbids weakening it to mask genuine biology — the discordance is brought to review instead.
