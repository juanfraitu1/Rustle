# parse_trflong Flow Recall — Phase 0 Diagnostic Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the read-only diagnostic harness that measures whether completing rustle's `parse_trflong` flow stubs can recover the 668 annotated isoforms StringTie finds but rustle misses, **without regressing precision** — producing a pre-registered go/no-go decision before any flow code is written.

**Architecture:** A locus-sliced cross-tool harness. For each st_only miss we slice both BAMs around the locus and run rustle + instrumented StringTie once, caching their `path_extracted` and `junction_accept` parity logs. Four analyses read the cache: (a) generation census, (b) recall-oracle ceiling, (c) separability test, (d) mechanism attribution. A final report applies the gate. No genome-wide whole-chromosome runs (OOM-safe); the aggregate Sn/Pr oracle reuses GTFs already on disk.

**Tech Stack:** Python 3 (stdlib only), `samtools`, the release `rustle` binary, the instrumented `tools/stringtie/stringtie` binary, the existing `parity_decisions` logging infra (`RUSTLE_PARITY_LOG` / `STRINGTIE_PARITY_LOG`).

---

## Context the engineer needs

- **Working dir for all paths below:** `/mnt/c/Users/jfris/Desktop/Rustle`.
- **Spec:** `docs/superpowers/specs/2026-06-09-parse-trflong-flow-recall-scope.md` (read it).
- **Inputs already on disk:**
  - `bench/copy_recovery_eval/results_genomewide/st_only_classified.jsonl` — the 668 misses, one JSON/line: `{"chrom","tid","cat",...}`. `tid` is a RefSeq transcript id like `rna-XM_004027824.5`.
  - `bench/copy_recovery_eval/results_genomewide/perchrom/<CHROM>/` — per chrom: `c.bam` (+`.bai`), `ref.gtf`, `rustle.gtf`, `st.gtf`, `r_fsm.txt`, `s_fsm.txt`.
  - `bench/copy_recovery_eval/results_genomewide/results.jsonl` — 26 per-chrom aggregate records.
- **Binaries:** `target/release/rustle` (release build) and `tools/stringtie/stringtie` (instrumented, built `make clean release`).
- **Parity log convention (both tools, same schema):** set `<TOOL>_PARITY_LOG=<file>`, `<TOOL>_PARITY_FILTER_CHROM=<chrom>`, `<TOOL>_PARITY_FILTER_STEPS=path_extracted,junction_accept`. `<TOOL>` = `RUSTLE` or `STRINGTIE`. Each line is `{"step","tool","chrom","start","end","strand","payload":{...}}`. A `{"step":"_meta",...}` line appears first.
- **`path_extracted` payload** (both tools): `{"source","cov","longcov","entry_abund","nexons","introns",...}`. `introns` is a comma-joined string of `firstbase-lastbase` 1-based inclusive intron coords. rustle adds `seed_tf`,`flux`,`raw_flow`; ST adds `tlen`.
- **`junction_accept` payload:** `{"accepted":bool,"nreads",...}` emitted per graph junction; `start`/`end` are the junction's 1-based inclusive intron coords.
- **Tool invocation:** `<bin> -L <slice.bam> -o <out.gtf>`. **Gotcha:** never pass `-o /dev/null` to StringTie — it errors on the temp dir; use a real path. ST writes coverage as a huge `tlen`-inflated `cov`; rely on `introns` for chain identity, not coords.
- **Determinism:** sort everything; do not use timestamps/randomness in outputs.
- **OOM safety:** only ever slice per-locus (a few-kb window, ~150 reads). Never run a tool on a whole-chromosome BAM here.

All new code lives under `bench/flow_recall_phase0/`.

## File structure

- `bench/flow_recall_phase0/lib.py` — shared helpers: GTF parsing, intron-chain extraction, per-locus slice+run+cache, parity-log parsing. One responsibility: turn a `(chrom, tid)` into parsed cross-tool logs.
- `bench/flow_recall_phase0/capture.py` — driver that calls `lib.ensure_locus_logs` over all 668 misses (resumable). Produces the log cache. The only expensive step.
- `bench/flow_recall_phase0/gen_census.py` — analysis (a): generated vs non-generated split.
- `bench/flow_recall_phase0/recall_oracle.py` — analysis (b): Sn/Pr ceiling + precision cost, from GTFs.
- `bench/flow_recall_phase0/separability.py` — analysis (c): feature vectors + separation metric (the gate input).
- `bench/flow_recall_phase0/attribute.py` — analysis (d): graph-missing vs flow-enumeration split via junction presence.
- `bench/flow_recall_phase0/gate_report.py` — aggregates (a)-(d), prints the go/no-go verdict.
- `bench/flow_recall_phase0/tests/` — fixture-based unit tests (no live tool runs; hand-written tiny log/GTF fixtures).

Cache layout: `bench/flow_recall_phase0/cache/<chrom>/<tid>/{r.jsonl,s.jsonl,r.gtf,s.gtf,slice.bam,.done}`.

---

## Task 1: Shared library — GTF + parity-log parsing

**Files:**
- Create: `bench/flow_recall_phase0/lib.py`
- Create: `bench/flow_recall_phase0/tests/test_lib_parse.py`
- Create: `bench/flow_recall_phase0/tests/fixtures/mini.gtf`
- Create: `bench/flow_recall_phase0/tests/fixtures/pe.jsonl`

- [ ] **Step 1: Write the failing test**

Create `bench/flow_recall_phase0/tests/fixtures/mini.gtf`:
```
chrA	x	transcript	100	900	.	+	.	transcript_id "t1";
chrA	x	exon	100	200	.	+	.	transcript_id "t1";
chrA	x	exon	301	500	.	+	.	transcript_id "t1";
chrA	x	exon	801	900	.	+	.	transcript_id "t1";
```
Create `bench/flow_recall_phase0/tests/fixtures/pe.jsonl`:
```
{"step":"_meta","tool":"rustle","payload":{"version":1}}
{"step":"path_extracted","tool":"rustle","chrom":"chrA","start":100,"end":900,"strand":"+","payload":{"cov":3.6,"longcov":3.0,"entry_abund":2.0,"nexons":3,"seed_tf":7,"introns":"201-300,501-800"}}
{"step":"junction_accept","tool":"rustle","chrom":"chrA","start":201,"end":300,"strand":"+","payload":{"accepted":true,"nreads":40}}
```
Create `bench/flow_recall_phase0/tests/test_lib_parse.py`:
```python
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib
FX = os.path.join(os.path.dirname(__file__), "fixtures")

def test_gtf_intron_chain():
    tx = lib.parse_gtf(os.path.join(FX, "mini.gtf"))
    assert tx["t1"]["strand"] == "+"
    assert tx["t1"]["introns"] == ((201, 300), (501, 800))
    assert tx["t1"]["span"] == (100, 900)

def test_chain_str_roundtrip():
    assert lib.chain_str(((201, 300), (501, 800))) == "201-300,501-800"
    assert lib.parse_chain_str("201-300,501-800") == ((201, 300), (501, 800))

def test_parse_path_extracted():
    evs = lib.read_parity(os.path.join(FX, "pe.jsonl"), "path_extracted")
    assert len(evs) == 1
    assert evs[0]["payload"]["introns"] == "201-300,501-800"

def test_parse_junctions():
    js = lib.read_parity(os.path.join(FX, "pe.jsonl"), "junction_accept")
    assert (201, 300) in {(j["start"], j["end"]) for j in js}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd bench/flow_recall_phase0 && python3 -m pytest tests/test_lib_parse.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'lib'`.

- [ ] **Step 3: Write minimal implementation**

Create `bench/flow_recall_phase0/lib.py`:
```python
"""Shared helpers for the Phase 0 flow-recall diagnostic."""
import json, os, re, subprocess

ROOT = "/mnt/c/Users/jfris/Desktop/Rustle"
RG = f"{ROOT}/bench/copy_recovery_eval/results_genomewide"
PERCHROM = f"{RG}/perchrom"
CACHE = f"{ROOT}/bench/flow_recall_phase0/cache"
RUSTLE = f"{ROOT}/target/release/rustle"
ST = f"{ROOT}/tools/stringtie/stringtie"

def parse_gtf(path):
    """transcript_id -> {strand, chrom, exons, introns(tuple), span, splice_sites(set)}"""
    tx = {}
    if not os.path.exists(path):
        return tx
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m:
            continue
        d = tx.setdefault(m.group(1), {"strand": f[6], "chrom": f[0], "exons": []})
        d["exons"].append((int(f[3]), int(f[4])))
    for d in tx.values():
        d["exons"].sort()
        ex = d["exons"]
        d["introns"] = tuple((ex[i][1] + 1, ex[i + 1][0] - 1) for i in range(len(ex) - 1))
        d["span"] = (ex[0][0], ex[-1][1]) if ex else (0, 0)
        d["splice_sites"] = {c for iv in d["introns"] for c in iv}
    return tx

def chain_str(introns):
    return ",".join(f"{a}-{b}" for a, b in introns)

def parse_chain_str(s):
    return tuple(tuple(int(x) for x in p.split("-")) for p in s.split(",") if "-" in p)

def read_parity(path, step):
    out = []
    if not os.path.exists(path):
        return out
    for line in open(path):
        try:
            o = json.loads(line)
        except ValueError:
            continue
        if o.get("step") == step:
            out.append(o)
    return out
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd bench/flow_recall_phase0 && python3 -m pytest tests/test_lib_parse.py -q`
Expected: PASS (4 passed).

- [ ] **Step 5: Commit**

```bash
git add bench/flow_recall_phase0/lib.py bench/flow_recall_phase0/tests/
git commit -m "feat(phase0): lib GTF + parity-log parsing helpers"
```

---

## Task 2: Per-locus capture (slice + run both tools + cache)

**Files:**
- Modify: `bench/flow_recall_phase0/lib.py` (add `ensure_locus_logs`)
- Create: `bench/flow_recall_phase0/capture.py`
- Create: `bench/flow_recall_phase0/tests/test_capture_smoke.py`

- [ ] **Step 1: Write the failing test**

Create `bench/flow_recall_phase0/tests/test_capture_smoke.py`:
```python
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

def test_ref_chain_lookup():
    # XM_004027824.5 on NC_073224.2 is a known container miss; rustle GENERATES it.
    ch = lib.ref_chain("NC_073224.2", "rna-XM_004027824.5")
    assert ch is not None
    span_lo, span_hi, introns = ch
    assert len(introns) == 10            # 11 exons -> 10 introns
    assert span_lo < span_hi
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd bench/flow_recall_phase0 && python3 -m pytest tests/test_capture_smoke.py -q`
Expected: FAIL — `AttributeError: module 'lib' has no attribute 'ref_chain'`.

- [ ] **Step 3: Write minimal implementation**

Append to `bench/flow_recall_phase0/lib.py`:
```python
def ref_chain(chrom, tid):
    """Return (span_lo, span_hi, introns_tuple) for an annotated transcript, or None."""
    ex = []
    path = f"{PERCHROM}/{chrom}/ref.gtf"
    if not os.path.exists(path):
        return None
    for line in open(path):
        f = line.split("\t")
        if len(f) < 9 or f[2] != "exon" or tid not in f[8]:
            continue
        ex.append((int(f[3]), int(f[4])))
    if not ex:
        return None
    ex.sort()
    introns = tuple((ex[i][1] + 1, ex[i + 1][0] - 1) for i in range(len(ex) - 1))
    return ex[0][0], ex[-1][1], introns

def _run_tool(binpath, log_var, slice_bam, chrom, out_gtf, log_file):
    env = dict(os.environ)
    env[log_var] = log_file
    env[log_var.replace("LOG", "FILTER_CHROM")] = chrom
    env[log_var.replace("LOG", "FILTER_STEPS")] = "path_extracted,junction_accept"
    subprocess.run([binpath, "-L", slice_bam, "-o", out_gtf], env=env,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def ensure_locus_logs(chrom, tid, pad=8000):
    """Idempotently slice the locus and run BOTH tools with parity logging. Cached."""
    wd = f"{CACHE}/{chrom}/{tid}"
    if os.path.exists(f"{wd}/.done"):
        return wd
    ch = ref_chain(chrom, tid)
    if ch is None:
        return None
    lo, hi, _ = ch
    os.makedirs(wd, exist_ok=True)
    sb = f"{wd}/slice.bam"
    reg = f"{chrom}:{max(1, lo - pad)}-{hi + pad}"
    with open(sb, "wb") as fh:
        subprocess.run(["samtools", "view", "-b", f"{PERCHROM}/{chrom}/c.bam", reg],
                       stdout=fh, stderr=subprocess.DEVNULL)
    subprocess.run(["samtools", "index", sb], stderr=subprocess.DEVNULL)
    _run_tool(RUSTLE, "RUSTLE_PARITY_LOG", sb, chrom, f"{wd}/r.gtf", f"{wd}/r.jsonl")
    _run_tool(ST, "STRINGTIE_PARITY_LOG", sb, chrom, f"{wd}/s.gtf", f"{wd}/s.jsonl")
    open(f"{wd}/.done", "w").close()
    return wd
```

Create `bench/flow_recall_phase0/capture.py`:
```python
"""Driver: ensure cached cross-tool logs exist for every st_only miss. Resumable."""
import json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

def main():
    misses = [json.loads(l) for l in open(f"{lib.RG}/st_only_classified.jsonl")]
    done = 0
    for i, m in enumerate(misses):
        wd = lib.ensure_locus_logs(m["chrom"], m["tid"])
        if wd:
            done += 1
        if (i + 1) % 25 == 0:
            print(f"{i+1}/{len(misses)} captured ({done} ok)", flush=True)
    print(f"DONE {done}/{len(misses)}")

if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd bench/flow_recall_phase0 && python3 -m pytest tests/test_capture_smoke.py -q`
Expected: PASS (1 passed).

- [ ] **Step 5: Run a single-locus live smoke check**

Run: `cd bench/flow_recall_phase0 && python3 -c "import lib; wd=lib.ensure_locus_logs('NC_073224.2','rna-XM_004027824.5'); print(wd); import os; print('r.jsonl', os.path.getsize(wd+'/r.jsonl')>0); print('s.jsonl', os.path.getsize(wd+'/s.jsonl')>0)"`
Expected: prints a cache path and `r.jsonl True` / `s.jsonl True` (both logs non-empty). If `s.jsonl` is empty, confirm `tools/stringtie/stringtie` exists and runs.

- [ ] **Step 6: Commit**

```bash
git add bench/flow_recall_phase0/lib.py bench/flow_recall_phase0/capture.py bench/flow_recall_phase0/tests/test_capture_smoke.py
git commit -m "feat(phase0): per-locus cross-tool capture (cached, resumable)"
```

---

## Task 3: Generation census (analysis a)

**Files:**
- Create: `bench/flow_recall_phase0/gen_census.py`
- Create: `bench/flow_recall_phase0/tests/test_gen_census.py`

Classifies each miss as **generated** (ref chain present in rustle `path_extracted`) or **non_generated**, then aggregates by category.

- [ ] **Step 1: Write the failing test**

Create `bench/flow_recall_phase0/tests/test_gen_census.py`:
```python
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import gen_census

def test_classify_generated():
    ref = ((201, 300), (501, 800))
    rustle_chains = {((201, 300), (501, 800))}
    assert gen_census.classify(ref, rustle_chains) == "generated"

def test_classify_non_generated():
    ref = ((201, 300), (501, 800))
    rustle_chains = {((201, 300), (777, 888))}   # different second intron
    assert gen_census.classify(ref, rustle_chains) == "non_generated"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd bench/flow_recall_phase0 && python3 -m pytest tests/test_gen_census.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'gen_census'`.

- [ ] **Step 3: Write minimal implementation**

Create `bench/flow_recall_phase0/gen_census.py`:
```python
"""Analysis (a): split the 668 misses into generated vs non-generated by rustle's flow."""
import json, sys, os
from collections import Counter, defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

def classify(ref_introns, rustle_chains):
    return "generated" if tuple(ref_introns) in rustle_chains else "non_generated"

def rustle_chain_set(wd):
    chains = set()
    for ev in lib.read_parity(f"{wd}/r.jsonl", "path_extracted"):
        s = ev["payload"].get("introns", "")
        if s:
            chains.add(lib.parse_chain_str(s))
    return chains

def main():
    misses = [json.loads(l) for l in open(f"{lib.RG}/st_only_classified.jsonl")]
    by_cat = defaultdict(Counter)
    rows = []
    for m in misses:
        wd = f"{lib.CACHE}/{m['chrom']}/{m['tid']}"
        if not os.path.exists(f"{wd}/.done"):
            by_cat[m["cat"]]["uncaptured"] += 1
            continue
        ch = lib.ref_chain(m["chrom"], m["tid"])
        if ch is None:
            continue
        verdict = classify(ch[2], rustle_chain_set(wd))
        by_cat[m["cat"]][verdict] += 1
        rows.append({"chrom": m["chrom"], "tid": m["tid"], "cat": m["cat"], "verdict": verdict})
    tot = Counter()
    for c, cc in by_cat.items():
        tot.update(cc)
    print("=== generation census (n=%d) ===" % len(rows))
    print(f"  generated:     {tot['generated']}")
    print(f"  non_generated: {tot['non_generated']}")
    print("\nby category:")
    for cat in sorted(by_cat):
        cc = by_cat[cat]
        print(f"  {cat:22s} gen={cc['generated']:4d} nongen={cc['non_generated']:4d}")
    with open(f"{lib.ROOT}/bench/flow_recall_phase0/gen_census.jsonl", "w") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")

if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd bench/flow_recall_phase0 && python3 -m pytest tests/test_gen_census.py -q`
Expected: PASS (2 passed).

- [ ] **Step 5: Commit**

```bash
git add bench/flow_recall_phase0/gen_census.py bench/flow_recall_phase0/tests/test_gen_census.py
git commit -m "feat(phase0): generation census analysis"
```

---

## Task 4: Recall-oracle ceiling (analysis b)

**Files:**
- Create: `bench/flow_recall_phase0/recall_oracle.py`
- Create: `bench/flow_recall_phase0/tests/test_recall_oracle.py`

Computes the genome-wide recall ceiling and precision cost from the per-chrom GTFs + FSM files already on disk (no tool re-runs). Recall numerator if rustle adopts ST's FSM set = `both + rustle_only + st_only` per chrom (the FSM union). Precision cost = number of ST `path_extracted` chains at the miss loci that are NOT annotated (would be added as FPs).

- [ ] **Step 1: Write the failing test**

Create `bench/flow_recall_phase0/tests/test_recall_oracle.py`:
```python
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import recall_oracle

def test_fsm_union_known_chrom():
    # NC_073227.2 results.jsonl: both=970, rustle_only=53, st_only=20.
    agg = recall_oracle.fsm_union_per_chrom()
    rec = agg["NC_073227.2"]
    assert rec["both"] == 970 and rec["rustle_only"] == 53 and rec["st_only"] == 20
    assert rec["union_fsm"] == 970 + 53 + 20      # 1043
    assert rec["rustle_fsm"] == 970 + 53          # 1023 = current rustle recall

def test_delta_recall_positive():
    agg = recall_oracle.fsm_union_per_chrom()
    # adopting ST's misses adds exactly st_only per chrom
    assert agg["NC_073227.2"]["union_fsm"] - agg["NC_073227.2"]["rustle_fsm"] == 20
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd bench/flow_recall_phase0 && python3 -m pytest tests/test_recall_oracle.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'recall_oracle'`.

- [ ] **Step 3: Write minimal implementation**

Create `bench/flow_recall_phase0/recall_oracle.py`:
```python
"""Analysis (b): recall ceiling (FSM union) + precision cost of adopting ST's missing chains."""
import json, sys, os
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

def fsm_union_per_chrom():
    agg = {}
    for line in open(f"{lib.RG}/results.jsonl"):
        r = json.loads(line)
        if not r.get("ok"):
            continue
        both, ro, so = r["both"], r["rustle_only"], r["st_only"]
        agg[r["chrom"]] = {
            "both": both, "rustle_only": ro, "st_only": so,
            "rustle_fsm": both + ro,
            "union_fsm": both + ro + so,
        }
    return agg

def annotated_chain_sets():
    """Per chrom: set of annotated intron chains (for precision-cost test)."""
    out = {}
    for line in open(f"{lib.RG}/results.jsonl"):
        r = json.loads(line)
        if not r.get("ok"):
            continue
        chrom = r["chrom"]
        ref = lib.parse_gtf(f"{lib.PERCHROM}/{chrom}/ref.gtf")
        out[chrom] = {d["introns"] for d in ref.values() if d["introns"]}
    return out

def precision_cost():
    """ST path_extracted chains at miss loci that are NOT annotated -> would-be FPs added."""
    misses = [json.loads(l) for l in open(f"{lib.RG}/st_only_classified.jsonl")]
    ann = annotated_chain_sets()
    extra = 0
    seen = set()
    for m in misses:
        wd = f"{lib.CACHE}/{m['chrom']}/{m['tid']}"
        if not os.path.exists(f"{wd}/.done"):
            continue
        rustle = {lib.parse_chain_str(e["payload"]["introns"])
                  for e in lib.read_parity(f"{wd}/r.jsonl", "path_extracted")
                  if e["payload"].get("introns")}
        for e in lib.read_parity(f"{wd}/s.jsonl", "path_extracted"):
            s = e["payload"].get("introns")
            if not s:
                continue
            ch = lib.parse_chain_str(s)
            key = (m["chrom"], ch)
            if key in seen:
                continue
            seen.add(key)
            if ch not in rustle and ch not in ann.get(m["chrom"], set()):
                extra += 1
    return extra

def main():
    agg = fsm_union_per_chrom()
    tot_rustle = sum(a["rustle_fsm"] for a in agg.values())
    tot_union = sum(a["union_fsm"] for a in agg.values())
    tot_so = sum(a["st_only"] for a in agg.values())
    print("=== recall-oracle ceiling (FSM union) ===")
    print(f"  current rustle FSM (genome-wide): {tot_rustle}")
    print(f"  ceiling if adopt ST misses:       {tot_union}  (+{tot_so})")
    fps = precision_cost()
    print(f"\n=== precision cost ===")
    print(f"  ST-extracted non-annotated chains at miss loci (would-be added FPs): {fps}")
    print(f"  crude recall:FP ratio = {tot_so}:{fps}"
          + (f"  ({tot_so/fps:.2f})" if fps else "  (no FPs)"))

if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd bench/flow_recall_phase0 && python3 -m pytest tests/test_recall_oracle.py -q`
Expected: PASS (2 passed).

- [ ] **Step 5: Commit**

```bash
git add bench/flow_recall_phase0/recall_oracle.py bench/flow_recall_phase0/tests/test_recall_oracle.py
git commit -m "feat(phase0): recall-oracle ceiling + precision-cost analysis"
```

---

## Task 5: Separability test (analysis c) — the gate input

**Files:**
- Create: `bench/flow_recall_phase0/separability.py`
- Create: `bench/flow_recall_phase0/tests/test_separability.py`

Builds feature vectors for **ST-extracted chains rustle does NOT generate**, labelled `real` (chain is annotated) vs `spurious` (not annotated), and computes how well each feature separates them. The label/feature join is what the gate needs: if no feature separates `real` from `spurious`, net-F1-positive is impossible.

- [ ] **Step 1: Write the failing test**

Create `bench/flow_recall_phase0/tests/test_separability.py`:
```python
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import separability as sep

def test_best_threshold_separates():
    # real label tends to higher longcov than spurious -> AUC should be high
    rows = [{"label": "real", "longcov": v} for v in (10, 12, 8, 15)] + \
           [{"label": "spurious", "longcov": v} for v in (1, 2, 1, 3)]
    auc = sep.auc_for_feature(rows, "longcov")
    assert auc > 0.9

def test_no_separation_gives_half_auc():
    rows = [{"label": "real", "longcov": v} for v in (5, 5, 5, 5)] + \
           [{"label": "spurious", "longcov": v} for v in (5, 5, 5, 5)]
    auc = sep.auc_for_feature(rows, "longcov")
    assert abs(auc - 0.5) < 1e-9
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd bench/flow_recall_phase0 && python3 -m pytest tests/test_separability.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'separability'`.

- [ ] **Step 3: Write minimal implementation**

Create `bench/flow_recall_phase0/separability.py`:
```python
"""Analysis (c): can any feature separate real-missing from spurious-missing ST chains?

AUC computed via the Mann-Whitney U statistic (ties counted as 0.5) -> no numpy needed.
"""
import json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

FEATURES = ["longcov", "cov", "entry_abund", "nexons"]

def auc_for_feature(rows, feat):
    pos = [r[feat] for r in rows if r["label"] == "real"]
    neg = [r[feat] for r in rows if r["label"] == "spurious"]
    if not pos or not neg:
        return float("nan")
    wins = 0.0
    for p in pos:
        for n in neg:
            wins += 1.0 if p > n else (0.5 if p == n else 0.0)
    return wins / (len(pos) * len(neg))

def build_rows():
    misses = [json.loads(l) for l in open(f"{lib.RG}/st_only_classified.jsonl")]
    ann_by_chrom = {}
    rows = []
    for m in misses:
        chrom = m["chrom"]
        wd = f"{lib.CACHE}/{chrom}/{m['tid']}"
        if not os.path.exists(f"{wd}/.done"):
            continue
        if chrom not in ann_by_chrom:
            ref = lib.parse_gtf(f"{lib.PERCHROM}/{chrom}/ref.gtf")
            ann_by_chrom[chrom] = {d["introns"] for d in ref.values() if d["introns"]}
        ann = ann_by_chrom[chrom]
        rustle = {lib.parse_chain_str(e["payload"]["introns"])
                  for e in lib.read_parity(f"{wd}/r.jsonl", "path_extracted")
                  if e["payload"].get("introns")}
        seen = set()
        for e in lib.read_parity(f"{wd}/s.jsonl", "path_extracted"):
            s = e["payload"].get("introns")
            if not s:
                continue
            ch = lib.parse_chain_str(s)
            if ch in rustle or ch in seen:
                continue
            seen.add(ch)
            p = e["payload"]
            rows.append({
                "label": "real" if ch in ann else "spurious",
                "longcov": float(p.get("longcov", 0.0)),
                "cov": float(p.get("cov", 0.0)),
                "entry_abund": float(p.get("entry_abund", 0.0)),
                "nexons": float(p.get("nexons", 0)),
            })
    return rows

def main():
    rows = build_rows()
    nr = sum(1 for r in rows if r["label"] == "real")
    ns = sum(1 for r in rows if r["label"] == "spurious")
    print(f"=== separability: ST chains rustle misses (real={nr} spurious={ns}) ===")
    for feat in FEATURES:
        auc = auc_for_feature(rows, feat)
        print(f"  {feat:12s} AUC(real>spurious) = {auc:.3f}")
    print("\nGate input: any AUC >> 0.5 (or << 0.5) means a discriminator exists.")
    with open(f"{lib.ROOT}/bench/flow_recall_phase0/separability_rows.jsonl", "w") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")

if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd bench/flow_recall_phase0 && python3 -m pytest tests/test_separability.py -q`
Expected: PASS (2 passed).

- [ ] **Step 5: Commit**

```bash
git add bench/flow_recall_phase0/separability.py bench/flow_recall_phase0/tests/test_separability.py
git commit -m "feat(phase0): separability test (real vs spurious missing chains)"
```

---

## Task 6: Mechanism attribution (analysis d)

**Files:**
- Create: `bench/flow_recall_phase0/attribute.py`
- Create: `bench/flow_recall_phase0/tests/test_attribute.py`

For each **non-generated** miss, split the cause using junction presence: if any of the ref chain's introns is absent from rustle's `junction_accept` set → `graph_missing` (deeper than flow; out of `parse_trflong` scope). If all introns are present but the chain isn't extracted → `flow_enumeration` (the `parse_trflong` target). This is the cheap, code-free attribution; deeper seed-loop instrumentation is deferred to Phase 1 only if needed.

- [ ] **Step 1: Write the failing test**

Create `bench/flow_recall_phase0/tests/test_attribute.py`:
```python
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import attribute

def test_flow_enumeration_when_all_junctions_present():
    ref = ((201, 300), (501, 800))
    junctions = {(201, 300), (501, 800), (900, 950)}
    assert attribute.attribute_cause(ref, junctions) == "flow_enumeration"

def test_graph_missing_when_a_junction_absent():
    ref = ((201, 300), (501, 800))
    junctions = {(201, 300)}              # second intron not accepted
    assert attribute.attribute_cause(ref, junctions) == "graph_missing"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd bench/flow_recall_phase0 && python3 -m pytest tests/test_attribute.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'attribute'`.

- [ ] **Step 3: Write minimal implementation**

Create `bench/flow_recall_phase0/attribute.py`:
```python
"""Analysis (d): for non-generated misses, split graph_missing vs flow_enumeration."""
import json, sys, os
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

def attribute_cause(ref_introns, rustle_junctions):
    for iv in ref_introns:
        if iv not in rustle_junctions:
            return "graph_missing"
    return "flow_enumeration"

def rustle_junctions(wd):
    return {(j["start"], j["end"]) for j in lib.read_parity(f"{wd}/r.jsonl", "junction_accept")}

def main():
    # Reuse gen_census output if present; else recompute generated/non per miss.
    census = {}
    gc = f"{lib.ROOT}/bench/flow_recall_phase0/gen_census.jsonl"
    if os.path.exists(gc):
        for line in open(gc):
            r = json.loads(line)
            census[(r["chrom"], r["tid"])] = r["verdict"]
    misses = [json.loads(l) for l in open(f"{lib.RG}/st_only_classified.jsonl")]
    cause = Counter()
    rows = []
    for m in misses:
        key = (m["chrom"], m["tid"])
        wd = f"{lib.CACHE}/{m['chrom']}/{m['tid']}"
        if not os.path.exists(f"{wd}/.done"):
            continue
        if census.get(key) == "generated":
            continue   # only attribute the non-generated set
        ch = lib.ref_chain(m["chrom"], m["tid"])
        if ch is None:
            continue
        c = attribute_cause(ch[2], rustle_junctions(wd))
        cause[c] += 1
        rows.append({"chrom": m["chrom"], "tid": m["tid"], "cat": m["cat"], "cause": c})
    print("=== non-generation attribution ===")
    for k, v in cause.most_common():
        print(f"  {k:16s} {v}")
    with open(f"{lib.ROOT}/bench/flow_recall_phase0/attribution.jsonl", "w") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")

if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd bench/flow_recall_phase0 && python3 -m pytest tests/test_attribute.py -q`
Expected: PASS (2 passed).

- [ ] **Step 5: Commit**

```bash
git add bench/flow_recall_phase0/attribute.py bench/flow_recall_phase0/tests/test_attribute.py
git commit -m "feat(phase0): non-generation attribution (graph vs flow)"
```

---

## Task 7: Gate report (aggregate + go/no-go verdict)

**Files:**
- Create: `bench/flow_recall_phase0/gate_report.py`
- Create: `bench/flow_recall_phase0/tests/test_gate_report.py`

Applies the pre-registered gate from the spec: proceed only if (1) recall ceiling ≥ 150 AND (2) a discriminator exists (some feature AUC ≥ 0.70 or ≤ 0.30 on the separability rows). Prints PROCEED or STOP with the numbers.

- [ ] **Step 1: Write the failing test**

Create `bench/flow_recall_phase0/tests/test_gate_report.py`:
```python
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import gate_report as gr

def test_gate_proceed():
    v = gr.decide(recall_ceiling=200, best_auc=0.82)
    assert v["verdict"] == "PROCEED"

def test_gate_stop_no_discriminator():
    v = gr.decide(recall_ceiling=200, best_auc=0.55)
    assert v["verdict"] == "STOP"
    assert "discriminator" in v["reason"]

def test_gate_stop_small_ceiling():
    v = gr.decide(recall_ceiling=40, best_auc=0.9)
    assert v["verdict"] == "STOP"
    assert "ceiling" in v["reason"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd bench/flow_recall_phase0 && python3 -m pytest tests/test_gate_report.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'gate_report'`.

- [ ] **Step 3: Write minimal implementation**

Create `bench/flow_recall_phase0/gate_report.py`:
```python
"""Apply the pre-registered go/no-go gate (spec: recall ceiling + discriminator)."""
import json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

RECALL_MIN = 150
AUC_HI, AUC_LO = 0.70, 0.30

def decide(recall_ceiling, best_auc):
    if recall_ceiling < RECALL_MIN:
        return {"verdict": "STOP", "reason": f"recall ceiling {recall_ceiling} < {RECALL_MIN} (small ceiling)"}
    if not (best_auc >= AUC_HI or best_auc <= AUC_LO):
        return {"verdict": "STOP", "reason": f"no discriminator (best AUC {best_auc:.2f} in [{AUC_LO},{AUC_HI}])"}
    return {"verdict": "PROCEED", "reason": f"ceiling {recall_ceiling} >= {RECALL_MIN}, discriminator AUC {best_auc:.2f}"}

def main():
    import recall_oracle, separability
    agg = recall_oracle.fsm_union_per_chrom()
    ceiling = sum(a["st_only"] for a in agg.values())
    rows = []
    sp = f"{lib.ROOT}/bench/flow_recall_phase0/separability_rows.jsonl"
    if os.path.exists(sp):
        rows = [json.loads(l) for l in open(sp)]
    aucs = {f: separability.auc_for_feature(rows, f) for f in separability.FEATURES} if rows else {}
    best_auc = max((max(a, 1 - a) for a in aucs.values() if a == a), default=0.5)
    v = decide(ceiling, best_auc)
    print("=== PHASE 0 GATE ===")
    print(f"  recall ceiling (st_only genome-wide): {ceiling}")
    print(f"  feature AUCs: " + ", ".join(f"{f}={a:.2f}" for f, a in aucs.items()))
    print(f"  best |discrimination| AUC: {best_auc:.2f}")
    print(f"\n  VERDICT: {v['verdict']} — {v['reason']}")

if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd bench/flow_recall_phase0 && python3 -m pytest tests/test_gate_report.py -q`
Expected: PASS (3 passed).

- [ ] **Step 5: Commit**

```bash
git add bench/flow_recall_phase0/gate_report.py bench/flow_recall_phase0/tests/test_gate_report.py
git commit -m "feat(phase0): go/no-go gate report"
```

---

## Task 8: Run the full diagnostic + record the verdict

**Files:**
- Create: `bench/flow_recall_phase0/PHASE0_RESULTS.md`

- [ ] **Step 1: Capture all 668 loci (the one expensive step)**

Run: `cd bench/flow_recall_phase0 && python3 capture.py 2>&1 | tail -5`
Expected: progress lines ending `DONE <n>/668` with `<n>` close to 668. Resumable — re-run if interrupted. Watch memory: `free -m` should stay well under limits (per-locus slices are tiny). If any tool hangs on a locus, note the tid and skip it (the cache `.done` sentinel makes the rest resumable).

- [ ] **Step 2: Run the four analyses**

Run:
```bash
cd bench/flow_recall_phase0
python3 gen_census.py
python3 recall_oracle.py
python3 separability.py
python3 attribute.py
python3 gate_report.py
```
Expected: each prints its summary; `gate_report.py` ends with `VERDICT: PROCEED` or `STOP`.

- [ ] **Step 3: Record results**

Write `bench/flow_recall_phase0/PHASE0_RESULTS.md` capturing: generation census (gen vs non-gen, by category), recall ceiling + precision cost, per-feature separability AUCs, attribution split (graph vs flow), and the gate verdict with the exact numbers. If PROCEED, name the implicated stub (from attribution: `flow_enumeration` dominant → `parse_trflong` seed/extension; note whether it's `seed_order_st` or `update_abundance_st` from a follow-up trace). If STOP, state which condition failed — this is the deliverable.

- [ ] **Step 4: Commit**

```bash
git add bench/flow_recall_phase0/PHASE0_RESULTS.md bench/flow_recall_phase0/*.jsonl
git commit -m "docs(phase0): diagnostic results + gate verdict"
```

- [ ] **Step 5: Update memory**

Add/update a memory entry (`project_st_only_genomewide_miss.md` or a new `project_flow_recall_phase0.md`) with the verdict, the recall ceiling, the best separability AUC, and — if PROCEED — which stub Phase 1 targets; if STOP, the falsified condition. Add the one-line pointer to `MEMORY.md`.

---

## Phase 1 (NOT in this plan)

Phase 1 (complete the implicated `parse_trflong_st` stub, env-gated, parity-validated) is **contingent on a PROCEED verdict** and gets its own spec+plan once Task 8 names the target stub and confirms a discriminator. Do not start it from this plan.

## Self-review notes

- **Spec coverage:** (a) Task 3, (b) Task 4, (c) Task 5, (d) Task 6, gate Task 7, run+verdict Task 8, OOM-safe locus slicing Task 2, contingent Phase 1 explicitly out-of-plan. All spec sections mapped.
- **Determinism:** all outputs iterate sorted/stable inputs; no timestamps/random.
- **Known-fixture tests:** every analysis has a unit test on hand-computed data (no live tool runs in tests); live runs are smoke-checked in Task 2 Step 5 and executed in Task 8.
