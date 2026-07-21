# Mechanism Transparency Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make Rustle's single method legible end-to-end and disclose every heuristic, so the advisor sees two mechanisms + two consequences (not four approaches) and cannot read incompleteness as withholding.

**Architecture:** A hand-maintained heuristic registry (`bench/mechanism/heuristics.toml`) is verified against source by a generator that opens each `file:line` and asserts the literal value still matches — the anti-drift device. Verification runs (V1–V4) on the Soto RNA families and a tandem removal-recovery positive control produce real numbers and a recovery figure. A single `bench/rustle_mechanism.html`, led by the two-sentence method (§0 of the spec), embeds the registry table and the V4 figure. Three existing artifacts that contradict the code are corrected.

**Tech Stack:** Rust (existing built binaries `target/release/{copy_assign,gw_family_catalog,family_define,asj}`), Python 3.14 (stdlib only — `tomllib` is built in), minimap2/samtools on PATH, self-contained HTML/CSS/JS.

## Global Constraints

- **Disclosure, not obfuscation:** every heuristic in the registry is truthful — the generator asserts each value against its `file:line` and fails loudly on drift. Never hand-edit `heuristics.tsv`; it is generated.
- **No false unification:** Define (exon-sum homology) and Assign (significance certificate) are presented as TWO mechanisms, never one. The "one test, three uses" claim is about the certificate only and must be verified in code before it appears in the artifact.
- **Honest floors shown, not hidden:** the identical-copy K=0 arm (recover number, abstain on sequence) ships as a first-class result.
- **Recovery eval is permutation-invariant set matching:** a family is a set of copies; {copy1,copy3,copy2} recovering {copy1,copy2,copy3} is a perfect recovery. Never penalize reordering.
- **Crash rule (WSL2):** `copy_assign` runs FOREGROUND, serial, small batches, outputs to `/home/juanfra/winloci_scratch` (NOT `/tmp`). No `nohup`, no waiter loops, no `pkill -f`. Heavy genome-wide binary runs may use harness-tracked background only.
- **Alignment presets are fixed project-wide:** simulated/real read alignment uses `-N 50` and the splice preset for RNA; never `--secondary=no` (yields 0 families).
- **Data locations:** reads/assemblies at `/mnt/linuxdisk/home/juanfraitu/winloci_data/` (bind-mounted; if absent, remount per the mount steps). Full gorilla BAM `GGO_mm.bam` + `GGO.fasta` and small `GGO_19_mm.bam` in `/home/juanfra/winloci_scratch/`.
- **Live family edge is `refine_families_exon_sum`** (asm20 id≥0.80 cov-of-shorter≥0.50, ≥2 distinct loci). `detect_families`/`core_recip≥0.13` is DEAD as an edge; `T_CORE=0.13` survives only as thin-locus rescue admission. Artifacts must say this.

---

## File Structure

- `bench/mechanism/heuristics.toml` — the registry: one `[[heuristic]]` per constant (name, file, line, value, stage, tier, kind, rationale).
- `bench/mechanism/gen_heuristics.py` — reads the toml, asserts each value against source, emits `bench/mechanism/heuristics.tsv`; non-zero exit on drift.
- `bench/mechanism/heuristics.tsv` — generated table (committed so the HTML can embed it and reviewers can diff it).
- `tests/mechanism/test_gen_heuristics.py` — unit tests for the generator (drift detection, tsv shape).
- `bench/mechanism/v1_spine_numbers.sh` — re-runs live commands, records the spine's decision numbers.
- `bench/mechanism/v2_v3_guards.py` — measures observed-max per inert-guard on the Soto families; toggles READ_CAP/POA_CAP.
- `bench/mechanism/sim_tandem.py` — V4a: build a 3-copy divergent tandem, simulate reads, degrade, recover, set-match.
- `bench/mechanism/real_tandem.py` — V4b: gorilla A119b vs mGorGor1 same-read-set recovery.
- `bench/mechanism/verification_results.json` — collected numbers from V1–V4 (consumed by the HTML build).
- `bench/rustle_mechanism.html` — the artifact.
- Corrections: `bench/quasiclique.html`, `bench/DEFINITIONS_FORMAL.md`, `figures/01_pipeline_overview.md` (+02/03/04), `analysis/family_graphs/docs/advisor_response.md`.

---

## Task 1: Registry schema + generator with drift detection (TDD)

**Files:**
- Create: `bench/mechanism/gen_heuristics.py`
- Create: `tests/mechanism/test_gen_heuristics.py`
- Create: `bench/mechanism/heuristics.toml` (seed with 3 real entries only; expanded in Task 2)

**Interfaces:**
- Produces: `load_registry(toml_path) -> list[dict]`; `verify_entry(entry, repo_root) -> str|None` (returns an error string if the literal `value` is NOT found on `file:line`, else `None`); `emit_tsv(entries, out_path) -> None`; CLI `python3 gen_heuristics.py --repo-root <root> --toml <path> --out <tsv>` exits non-zero if any entry fails verification.
- TSV columns (tab-separated, this exact order): `stage	tier	kind	name	value	file	line	rationale`.

- [ ] **Step 1: Write the failing test**

```python
# tests/mechanism/test_gen_heuristics.py
import subprocess, sys, tempfile, os, pathlib
GEN = pathlib.Path(__file__).resolve().parents[2] / "bench/mechanism/gen_heuristics.py"
REPO = pathlib.Path(__file__).resolve().parents[2]

def _toml(body):
    f = tempfile.NamedTemporaryFile("w", suffix=".toml", delete=False)
    f.write(body); f.close(); return f.name

def test_verify_passes_on_real_line():
    # denovo_pipeline.rs:2370 really is `min_identity: 0.80,`
    t = _toml('[[heuristic]]\nname="edge_id"\nfile="src/rustle/vg_family/denovo_pipeline.rs"\n'
              'line=2370\nvalue="0.80"\nstage="O1-edge"\ntier="decision"\nkind="arbitrary"\nrationale="edge identity floor"\n')
    out = tempfile.mktemp(suffix=".tsv")
    r = subprocess.run([sys.executable, str(GEN), "--repo-root", str(REPO), "--toml", t, "--out", out],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    rows = open(out).read().splitlines()
    assert rows[0] == "stage\ttier\tkind\tname\tvalue\tfile\tline\trationale"
    assert any("edge_id" in row and "0.80" in row for row in rows[1:])

def test_verify_fails_on_drifted_value():
    # value that is NOT on that line -> non-zero exit, names the entry
    t = _toml('[[heuristic]]\nname="edge_id"\nfile="src/rustle/vg_family/denovo_pipeline.rs"\n'
              'line=2370\nvalue="0.99"\nstage="O1-edge"\ntier="decision"\nkind="arbitrary"\nrationale="x"\n')
    out = tempfile.mktemp(suffix=".tsv")
    r = subprocess.run([sys.executable, str(GEN), "--repo-root", str(REPO), "--toml", t, "--out", out],
                       capture_output=True, text=True)
    assert r.returncode != 0
    assert "edge_id" in (r.stderr + r.stdout)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd /mnt/c/Users/jfris/Desktop/Rustle && python3 -m pytest tests/mechanism/test_gen_heuristics.py -v`
Expected: FAIL (gen_heuristics.py does not exist).

- [ ] **Step 3: Write minimal implementation**

```python
# bench/mechanism/gen_heuristics.py
"""Verify each registry heuristic against its source file:line, then emit a TSV.
Anti-drift: if a constant's literal value is no longer on its recorded line, fail loudly."""
import argparse, sys, tomllib, pathlib

COLS = ["stage", "tier", "kind", "name", "value", "file", "line", "rationale"]

def load_registry(toml_path):
    with open(toml_path, "rb") as f:
        return tomllib.load(f).get("heuristic", [])

def verify_entry(entry, repo_root):
    p = pathlib.Path(repo_root) / entry["file"]
    if not p.exists():
        return f'{entry["name"]}: file not found: {entry["file"]}'
    lines = p.read_text(errors="replace").splitlines()
    ln = entry["line"]
    if ln < 1 or ln > len(lines):
        return f'{entry["name"]}: line {ln} out of range in {entry["file"]}'
    if str(entry["value"]) not in lines[ln - 1]:
        return (f'{entry["name"]}: value {entry["value"]!r} not on '
                f'{entry["file"]}:{ln} — got: {lines[ln-1].strip()!r}')
    return None

def emit_tsv(entries, out_path):
    with open(out_path, "w") as f:
        f.write("\t".join(COLS) + "\n")
        for e in sorted(entries, key=lambda e: (e["stage"], e["tier"], e["name"])):
            f.write("\t".join(str(e.get(c, "")) for c in COLS) + "\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", required=True)
    ap.add_argument("--toml", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    entries = load_registry(a.toml)
    errors = [msg for e in entries if (msg := verify_entry(e, a.repo_root))]
    if errors:
        print("DRIFT DETECTED — registry does not match source:", file=sys.stderr)
        for m in errors:
            print("  " + m, file=sys.stderr)
        sys.exit(1)
    emit_tsv(entries, a.out)
    print(f"verified {len(entries)} heuristics -> {a.out}")

if __name__ == "__main__":
    main()
```

Also create the 3-entry seed so the generator has something real to run in Task 3's smoke:

```toml
# bench/mechanism/heuristics.toml  (seed; expanded in Task 2)
[[heuristic]]
name = "RefineParams::min_identity"
file = "src/rustle/vg_family/denovo_pipeline.rs"
line = 2370
value = "0.80"
stage = "O1-edge"
tier = "decision"
kind = "arbitrary"
rationale = "asm20 exon-sum homology edge identity floor — THE family edge test"

[[heuristic]]
name = "GATE_MIN_READS"
file = "src/rustle/vg_family/denovo_assemble.rs"
line = 27
value = "3"
stage = "pre-O1"
tier = "decision"
kind = "arbitrary"
rationale = "locus admission gate — the pipeline's unit of existence"

[[heuristic]]
name = "copy_assign alpha (Bonferroni thr)"
file = "src/rustle/vg_family/copy_assign.rs"
line = 453
value = "n.saturating_sub(1)"
stage = "O2"
tier = "derived"
kind = "derived"
rationale = "certificate threshold alpha/(n-1), union bound over n-1 competitors"
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /mnt/c/Users/jfris/Desktop/Rustle && python3 -m pytest tests/mechanism/test_gen_heuristics.py -v`
Expected: PASS (2 tests).

- [ ] **Step 5: Commit**

```bash
git add bench/mechanism/gen_heuristics.py bench/mechanism/heuristics.toml tests/mechanism/test_gen_heuristics.py
git commit -m "feat(mechanism): heuristic registry generator with verify-against-source drift detection

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: Populate the full registry (all ~200 heuristics, tiered)

**Files:**
- Modify: `bench/mechanism/heuristics.toml` (expand from 3 seed entries to the full set)
- Create: `bench/mechanism/heuristics.tsv` (generated output, committed)

**Interfaces:**
- Consumes: `gen_heuristics.py` from Task 1.
- Produces: a `heuristics.tsv` covering every heuristic from the design sweep, each with a `tier` in {`decision`, `result-guard`, `inert-guard`, `derived`} and a `stage` in {`pre-O1`, `O1-prefilter`, `O1-edge`, `O1-cluster`, `O1-count`, `O2`, `O3`, `O4`}.

The exhaustive source list is the design-time sweep in the spec's §2.3 and the inventory tables. Enter each as a `[[heuristic]]` block. The tiers are pre-decided:
- **decision:** `RefineParams::min_identity` (denovo_pipeline.rs:2370, 0.80), `RefineParams::min_coverage` (denovo_pipeline.rs:2371, 0.50), `GAMMA` (family_definition.rs:173, 0.20), `GATE_MIN_READS` (denovo_assemble.rs:27, 3), `read_conflict delta` (read_conflict.rs:48), `de_max` (read_conflict.rs:49), `copy_assign alpha` (copy_assign.rs:183 + thr at :453), `junction_weight`, `error_rate` (copy_assign.rs), `min-copies` default 2.
- **result-guard:** `POA_CAP=20000` (copy_assign_pipeline.rs:524), `READ_CAP=6000` (o2_materialize.rs:250), sensitive-tier `0.70` literal (denovo_pipeline.rs:2506) vs `sensitive_identity=0.60` (denovo_pipeline.rs:2377/2645), and each of the five per-base error rates (0.003 copy_assign.rs; 0.01 copy_assign_pipeline.rs:948 QUANT_ERROR; 0.01/0.005/0.05 mosaic.rs; 0.001 read_conflict.rs:80).
- **inert-guard (tier justified by Task 5 measurement):** `MAX_MEMBERS=30` (multi_repeat_bridge.rs:71), `MAX_LOCI=60` (recombinant_split.rs:55), `LEN_CAP` 9000 (family_rescue.rs:50) / 20000 (family_detect.rs:41), `MAX_SPAN`, `PAIR_CAP`, `MAX_PAIRS`, cost caps.
- **derived:** `alpha/(n-1)` Bonferroni (copy_assign.rs:453), exact Poisson-binomial DP (copy_assign.rs:142-165), `tau_from_p` (copy_assign.rs:131-135), BH-FDR q (asj.rs:34), permutation p (linearize.rs), log-gamma numerics.

For any entry whose exact value string is tricky (floats formatted differently, or the value spans the line oddly), set `value` to the shortest literal substring that uniquely appears on that line (e.g. `"0.80"`, `"GAMMA: f64 = 0.20"`, `"6000"`). The generator only checks substring presence.

- [ ] **Step 1: Expand the toml** — add one `[[heuristic]]` per row from the sweep tables (~200). Group in the file by stage for editing sanity (order does not matter; the generator sorts).

- [ ] **Step 2: Generate and verify against source**

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
python3 bench/mechanism/gen_heuristics.py --repo-root . --toml bench/mechanism/heuristics.toml --out bench/mechanism/heuristics.tsv
```
Expected: `verified N heuristics -> bench/mechanism/heuristics.tsv` with exit 0. If any DRIFT line appears, the recorded `line`/`value` is stale — open that `file:line`, correct the registry entry to the current value, re-run. Do NOT change source to match the registry.

- [ ] **Step 3: Sanity-check tier counts**

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
cut -f2 bench/mechanism/heuristics.tsv | tail -n +2 | sort | uniq -c
```
Expected: four tiers present; `decision` ≈ 10–15, `derived` ≈ 15–20, the rest split between `result-guard` and `inert-guard`. If `decision` is much larger than ~15, re-check whether a constant is really a decision or machinery.

- [ ] **Step 4: Commit**

```bash
git add bench/mechanism/heuristics.toml bench/mechanism/heuristics.tsv
git commit -m "feat(mechanism): full heuristic registry (~200 constants, tiered, verified vs source)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: V1 — verify the spine numbers on real data + the "one test, three uses" claim

**Files:**
- Create: `bench/mechanism/v1_spine_numbers.sh`
- Create: `bench/mechanism/verification_results.json` (seed with V1 block)

**Interfaces:**
- Produces: `verification_results.json` with a top-level `"v1"` object: `{"edge": {"id":0.80,"cov":0.50,"loci":">=2"}, "gstm": {"n_copies":<int>}, "asj": {...}, "certificate_call_sites": [3 paths]}`.

- [ ] **Step 1: Confirm the certificate reuse claim in code (must pass before the artifact asserts it)**

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
for f in src/rustle/vg_family/copy_assign.rs src/rustle/vg_family/absent_copy.rs src/rustle/vg_family/collapse_gate.rs; do
  echo "=== $f ==="
  grep -nE "alpha|poisson|saturating_sub\(1\)|1e-3|0.001" "$f" | head -4
done
```
Expected: all three files reference the same significance form (alpha, Bonferroni `saturating_sub(1)` or a shared `min_p`/`min_p_distinct` helper). Record the three paths in the JSON `certificate_call_sites`. If any file does NOT share the certificate, the "three uses" claim is downgraded to the verified subset — do not overstate.

- [ ] **Step 2: Re-run the GSTM assignment (the spine's O2 number)**

Run (foreground, per crash rule):
```bash
cd /home/juanfra/winloci_scratch
/mnt/c/Users/jfris/Desktop/Rustle/target/release/copy_assign \
  --bam GGO_mm.bam --fasta GGO.fasta \
  --region NC_073224.2:129160000-129230000 \
  --homology-primary --min-copies 2 --dump-psv --out mech_gstm
cut -f1-6 mech_gstm.families.tsv | head
```
Expected: a family row with `n_copies=3` (GSTM). Record `n_copies` into the JSON.

- [ ] **Step 3: Confirm the edge values the binary applies**

Run:
```bash
/mnt/c/Users/jfris/Desktop/Rustle/target/release/gw_family_catalog --help 2>&1 | grep -iE "asm20 identity|cov-of-shorter|DISTINCT loci"
```
Expected: the help text states `asm20 identity>=0.80 cov-of-shorter>=0.50` and `>=2 spatially-DISTINCT loci` — the spine's O1-edge numbers, confirmed from the shipped binary. Record into JSON `edge`.

- [ ] **Step 4: Write the driver script that captures the above into JSON**

```bash
# bench/mechanism/v1_spine_numbers.sh
#!/usr/bin/env bash
# V1: verify the spine's decision numbers on real data. Foreground; outputs to winloci_scratch.
set -euo pipefail
RUSTLE=/mnt/c/Users/jfris/Desktop/Rustle
SCRATCH=/home/juanfra/winloci_scratch
BIN=$RUSTLE/target/release
cd "$SCRATCH"
"$BIN/copy_assign" --bam GGO_mm.bam --fasta GGO.fasta \
  --region NC_073224.2:129160000-129230000 \
  --homology-primary --min-copies 2 --out mech_gstm
NCOPIES=$(awk -F'\t' 'NR==2{print $2}' mech_gstm.families.tsv)
echo "GSTM n_copies=$NCOPIES"
# certificate reuse:
for f in copy_assign absent_copy collapse_gate; do
  grep -q "saturating_sub(1)\|min_p\|alpha" "$RUSTLE/src/rustle/vg_family/$f.rs" \
    && echo "certificate present: $f.rs"
done
```

Run: `bash bench/mechanism/v1_spine_numbers.sh` — Expected: `GSTM n_copies=3` and three `certificate present` lines.

- [ ] **Step 5: Commit**

```bash
git add bench/mechanism/v1_spine_numbers.sh bench/mechanism/verification_results.json
git commit -m "test(mechanism): V1 spine-number + certificate-reuse verification on real GGO data

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: V2/V3 — inert-guard proof + result-changer toggles on the Soto families

**Files:**
- Create: `bench/mechanism/v2_v3_guards.py`
- Modify: `bench/mechanism/verification_results.json` (add `"v2"`, `"v3"` blocks)

**Interfaces:**
- Consumes: `bench/soto/a119b_detected_families.tsv` (66 families), `bench/soto/80_fams.chr.bed`.
- Produces: JSON `"v2": {"<guard_name>": {"guard": <num>, "observed_max": <num>, "fires": <bool>}}` for each inert-guard; `"v3": {"READ_CAP": {...}, "POA_CAP": {...}}` deltas.

- [ ] **Step 1: Write the failing test for the guard-classification logic**

```python
# tests/mechanism/test_guards.py
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "bench/mechanism"))
from v2_v3_guards import classify_guard

def test_guard_inert_when_observed_below():
    assert classify_guard(guard=6000, observed_max=812)["fires"] is False

def test_guard_fires_when_observed_meets():
    r = classify_guard(guard=30, observed_max=44)
    assert r["fires"] is True
```

- [ ] **Step 2: Run to verify it fails**

Run: `cd /mnt/c/Users/jfris/Desktop/Rustle && python3 -m pytest tests/mechanism/test_guards.py -v`
Expected: FAIL (module missing).

- [ ] **Step 3: Implement the measurement + classification**

```python
# bench/mechanism/v2_v3_guards.py
"""V2: prove each inert-guard sits outside the data's range on the Soto families
(observed_max < guard => never fires => cannot move the numbers).
V3: toggle the two known result-changers and report the delta honestly."""
import json, subprocess, pathlib, argparse

REPO = pathlib.Path(__file__).resolve().parents[2]
SOTO_FAMS = REPO / "bench/soto/a119b_detected_families.tsv"

def classify_guard(guard, observed_max):
    return {"guard": guard, "observed_max": observed_max, "fires": observed_max >= guard}

def soto_family_sizes():
    """max copies (n_copies col) and max member count across the 66 Soto families."""
    sizes = []
    with open(SOTO_FAMS) as f:
        next(f)
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) > 1 and cols[1].isdigit():
                sizes.append(int(cols[1]))
    return sizes

def measure_v2():
    sizes = soto_family_sizes()
    max_members = max(sizes) if sizes else 0
    # MAX_MEMBERS=30 and MAX_LOCI=60 guard the repeat-bridge / recombinant-split gates.
    return {
        "MAX_MEMBERS(30)": classify_guard(30, max_members),
        "MAX_LOCI(60)": classify_guard(60, max_members),
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(REPO / "bench/mechanism/verification_results.json"))
    a = ap.parse_args()
    out = pathlib.Path(a.out)
    data = json.loads(out.read_text()) if out.exists() else {}
    data["v2"] = measure_v2()
    out.write_text(json.dumps(data, indent=2))
    print("V2 guard classification:")
    for k, v in data["v2"].items():
        print(f"  {k}: observed_max={v['observed_max']} fires={v['fires']}")

if __name__ == "__main__":
    main()
```

Note on scope: `MAX_MEMBERS`/`MAX_LOCI` are measurable directly from the Soto family table (fast, no BAM). READ_CAP/POA_CAP observed-max need per-family read pools — measure those in Step 5 by extracting reads at each Soto region and counting; if the full pass is too heavy, sample the 10 largest Soto families (by `avg_reads`) and label the result as a sample in the JSON.

- [ ] **Step 4: Run tests + the V2 measurement**

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
python3 -m pytest tests/mechanism/test_guards.py -v
python3 bench/mechanism/v2_v3_guards.py
```
Expected: tests PASS; V2 prints `MAX_MEMBERS(30)` and `MAX_LOCI(60)` with `fires=` set from the real max family size. Any guard with `fires=True` MUST be moved from `inert-guard` to `result-guard` in `heuristics.toml` (re-run Task 2's generator afterward).

- [ ] **Step 5: V3 — toggle READ_CAP and POA_CAP on the largest affected Soto locus**

Run (foreground; pick the highest-`avg_reads` Soto family region):
```bash
cd /home/juanfra/winloci_scratch
RUSTLE=/mnt/c/Users/jfris/Desktop/Rustle
# baseline (default caps):
"$RUSTLE/target/release/copy_assign" --bam GGO_mm.bam --fasta GGO.fasta \
  --region <BIG_SOTO_REGION> --homology-primary --min-copies 2 --out mech_v3_base
md5sum mech_v3_base.assignments.tsv
# NOTE: READ_CAP/POA_CAP are compile-time consts (o2_materialize.rs:250, copy_assign_pipeline.rs:524),
# not CLI flags. If the observed pool for this locus is below the cap (from V2-style extraction),
# record fires=false and that the cap could not have altered THIS result. If a genuine toggle is
# needed, that requires a one-line const change + rebuild — defer to deliverable B and record the
# observed pool size + cap here as the honest evidence.
```
Record into JSON `"v3"`: for each cap, the observed pool/length for the biggest Soto locus, the cap value, and whether it could fire. Report the truth (including "not togglable without rebuild; observed below cap so inert on this locus").

- [ ] **Step 6: Commit**

```bash
git add bench/mechanism/v2_v3_guards.py tests/mechanism/test_guards.py bench/mechanism/verification_results.json bench/mechanism/heuristics.toml bench/mechanism/heuristics.tsv
git commit -m "test(mechanism): V2/V3 guard classification on the 66 Soto families

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: V4a — simulation tandem removal-recovery (divergent + identical arms, TDD)

**Files:**
- Create: `bench/mechanism/sim_tandem.py`
- Create: `tests/mechanism/test_sim_tandem.py`
- Modify: `bench/mechanism/verification_results.json` (add `"v4a"` block)

**Interfaces:**
- Produces: `build_tandem(seed, n_copies, divergence, seed_rng) -> (ref_fasta_str, list[copy_seq])`; `set_match(recovered, truth) -> list[(rec_idx, truth_idx, identity)]` (permutation-invariant, greedy best-identity matching); JSON `"v4a": {"divergent": {"n_recovered":int,"copy_number_correct":bool,"min_identity":float}, "identical": {"copy_number_correct":bool,"assign_abstained":bool}}`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/mechanism/test_sim_tandem.py
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "bench/mechanism"))
from sim_tandem import build_tandem, set_match

def test_build_tandem_plants_divergence():
    ref, copies = build_tandem("ACGTACGTACGTACGT"*20, n_copies=3, divergence=0.02, seed_rng=7)
    assert len(copies) == 3
    # divergent copies must differ pairwise (at least one PSV)
    assert copies[0] != copies[1] and copies[1] != copies[2]

def test_build_tandem_identical_when_zero_divergence():
    ref, copies = build_tandem("ACGT"*100, n_copies=3, divergence=0.0, seed_rng=7)
    assert copies[0] == copies[1] == copies[2]

def test_set_match_is_permutation_invariant():
    truth = ["AAAA", "CCCC", "GGGG"]
    recovered = ["GGGG", "AAAA", "CCCC"]   # reordered
    m = set_match(recovered, truth)
    assert len(m) == 3
    assert all(identity == 1.0 for _, _, identity in m)
```

- [ ] **Step 2: Run to verify failure**

Run: `cd /mnt/c/Users/jfris/Desktop/Rustle && python3 -m pytest tests/mechanism/test_sim_tandem.py -v`
Expected: FAIL (module missing).

- [ ] **Step 3: Implement build + set-match (pure functions first)**

```python
# bench/mechanism/sim_tandem.py
"""V4a: simulate a 3-copy tandem array, delete copies from the reference, recover
via the certificate, and set-match recovered vs deleted copies (permutation-invariant).
Divergent arm proves recovery; identical arm proves the honest K=0 abstain floor.
Uses a seeded PRNG (no Math.random / wall-clock) so runs are reproducible."""
import random, pathlib, subprocess, json, argparse

BASES = "ACGT"

def build_tandem(seed, n_copies, divergence, seed_rng):
    """Return (reference_with_n_tandem_copies, [copy_seqs]). Divergence = per-base
    substitution prob planted independently per copy (shared coordinate frame => PSVs)."""
    rng = random.Random(seed_rng)
    base = list(seed)
    copies = []
    for _ in range(n_copies):
        c = base[:]
        if divergence > 0:
            for i in range(len(c)):
                if rng.random() < divergence:
                    c[i] = rng.choice([b for b in BASES if b != c[i]])
        copies.append("".join(c))
    # simple tandem: copies concatenated with a short spacer
    spacer = "N" * 50
    ref = spacer.join(copies)
    return ref, copies

def _identity(a, b):
    n = min(len(a), len(b))
    if n == 0:
        return 0.0
    same = sum(1 for i in range(n) if a[i] == b[i])
    return same / max(len(a), len(b))

def set_match(recovered, truth):
    """Greedy permutation-invariant matching: each recovered copy -> best-identity truth copy."""
    matches = []
    used = set()
    for ri, r in enumerate(recovered):
        best, best_id = None, -1.0
        for ti, t in enumerate(truth):
            if ti in used:
                continue
            idv = _identity(r, t)
            if idv > best_id:
                best, best_id = ti, idv
        if best is not None:
            used.add(best)
            matches.append((ri, best, best_id))
    return matches
```

- [ ] **Step 4: Run pure-function tests**

Run: `cd /mnt/c/Users/jfris/Desktop/Rustle && python3 -m pytest tests/mechanism/test_sim_tandem.py -v`
Expected: PASS (3 tests).

- [ ] **Step 5: Add the end-to-end simulate→align→recover driver (uses real binaries)**

Append a `main()` to `sim_tandem.py` that: (1) builds a divergent 3-copy tandem from a real IsoSeq-derived seed; (2) simulates reads by sampling each copy with a planted per-base HiFi error (~0.003) using the seeded PRNG; (3) writes the degraded reference (delete copy index 2, then copies 1&2); (4) aligns with `minimap2 -ax splice:hq -N 50` (NOT `--secondary=no`) + `samtools sort/index`; (5) runs `copy_assign --absent-copies --linearize` foreground to recover; (6) reads back the recovered copy consensus and `set_match`es against the deleted copies; (7) repeats with `divergence=0.0` for the identical arm and asserts copy-number recovered but assignment abstains (Tied). Write results to JSON `"v4a"`.

```python
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scratch", default="/home/juanfra/winloci_scratch/mech_sim")
    ap.add_argument("--out", default=str(pathlib.Path(__file__).resolve().parents[2] / "bench/mechanism/verification_results.json"))
    a = ap.parse_args()
    work = pathlib.Path(a.scratch); work.mkdir(parents=True, exist_ok=True)
    seed_seq = ("ATGGCACCGTTAGCCTAGGCTAACGGTTAACCGGTACGTAGCTAGCTAGGCTA" * 30)
    result = {}
    for arm, div in (("divergent", 0.02), ("identical", 0.0)):
        ref, copies = build_tandem(seed_seq, 3, div, seed_rng=11)
        # ... write ref.fa, simulate reads.fq (seeded), minimap2 align, copy_assign --absent-copies,
        #     parse <out>.dna_needs.tsv / recovered consensus, set_match vs the deleted copies ...
        # Record: copy_number_correct, min_identity (divergent), assign_abstained (identical).
        result[arm] = {"n_copies_planted": 3, "divergence": div}  # fill with measured fields
    out = pathlib.Path(a.out)
    data = json.loads(out.read_text()) if out.exists() else {}
    data["v4a"] = result
    out.write_text(json.dumps(data, indent=2))
    print("V4a:", json.dumps(result, indent=2))

if __name__ == "__main__":
    main()
```

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle && python3 bench/mechanism/sim_tandem.py
```
Expected: `v4a.divergent.copy_number_correct=true` with `min_identity` high (>0.95 to the deleted copy); `v4a.identical.assign_abstained=true` with copy-number recovered. Record both. If the divergent arm does NOT recover, STOP — that is a real finding about the method, not a script bug; systematic-debugging before proceeding.

- [ ] **Step 6: Commit**

```bash
git add bench/mechanism/sim_tandem.py tests/mechanism/test_sim_tandem.py bench/mechanism/verification_results.json
git commit -m "test(mechanism): V4a simulated tandem removal-recovery (divergent + identical K=0 floor)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 6: V4b — real gorilla tandem removal-recovery (same-read-set ground truth)

**Files:**
- Create: `bench/mechanism/real_tandem.py`
- Modify: `bench/mechanism/verification_results.json` (add `"v4b"` block)

**Interfaces:**
- Consumes: `A119b_ds.bam` + `mGorGor1` index on the mounted disk; a chosen Soto tandem family (co-located, ≥3 copies).
- Produces: JSON `"v4b": {"locus":str,"individual":"A119b","truth_source":"same-read-set vs full mGorGor1","n_copies_supported":int,"recovered":bool,"identity_to_intact":float}`.

- [ ] **Step 1: Pick a co-located Soto tandem family with ≥3 copies**

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
awk -F'\t' 'NR>1 && $3==1 && $2>=3 {print $1"\t"$2"\t"$4}' bench/soto/a119b_detected_families.tsv | head
```
Expected: one or more same-chrom (`n_chroms==1`) families with `n_copies>=3`. Choose one; record its `family_id` and locus. (Same-chrom = a true tandem array, the intended substrate.)

- [ ] **Step 2: Establish same-read-set ground truth (align A119b to full mGorGor1 at that locus)**

Run (foreground; region-restricted so it stays light):
```bash
cd /home/juanfra/winloci_scratch
D=/mnt/linuxdisk/home/juanfraitu/winloci_data
RUSTLE=/mnt/c/Users/jfris/Desktop/Rustle
# extract the family's reads and realign to gorilla, then confirm 3 divergent copies are supported:
"$RUSTLE/target/release/copy_assign" --bam "$D/A119b.t2t.bam" --fasta "$D/mGorGor1.mat.fasta" \
  --region <GORILLA_LOCUS> --homology-primary --min-copies 2 --dump-psv --out mech_real_full
awk -F'\t' 'NR==2{print "supported copies:", $2}' mech_real_full.families.tsv
```
Expected: `supported copies: >=3` for this read set against the intact assembly — that IS the ground truth. (If the gorilla FASTA is only present as `.mmi`, build/obtain the FASTA or pick a locus reachable from the human-coord catalog; record which reference path was used.)

- [ ] **Step 3: Degrade the assembly (delete one copy) and recover**

`real_tandem.py`: extract the intact copies at the locus from `mech_real_full` (the `.copies` / phased output), write a degraded reference FASTA with one copy removed, re-align the SAME extracted reads, run `copy_assign --absent-copies --linearize`, and compare the recovered copy to the removed intact copy via `set_match` (imported from `sim_tandem`).

```python
# bench/mechanism/real_tandem.py
"""V4b: real gorilla A119b tandem removal-recovery. Ground truth = the copies THIS read
set supports against the intact mGorGor1 assembly (same-read-set => no cross-individual
confound, since A119b != Kamilah). Delete a copy, re-align the same reads, recover, compare."""
import json, pathlib, argparse
from sim_tandem import set_match  # permutation-invariant matching, reused

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--locus", required=True)
    ap.add_argument("--intact-prefix", required=True)   # mech_real_full
    ap.add_argument("--scratch", default="/home/juanfra/winloci_scratch/mech_real")
    ap.add_argument("--out", default=str(pathlib.Path(__file__).resolve().parents[2] / "bench/mechanism/verification_results.json"))
    a = ap.parse_args()
    # ... build degraded ref (drop 1 copy), realign same reads, copy_assign --absent-copies,
    #     set_match recovered vs the dropped intact copy ...
    v4b = {"locus": a.locus, "individual": "A119b",
           "truth_source": "same-read-set vs full mGorGor1"}  # fill measured fields
    out = pathlib.Path(a.out)
    data = json.loads(out.read_text()) if out.exists() else {}
    data["v4b"] = v4b
    out.write_text(json.dumps(data, indent=2))
    print("V4b:", json.dumps(v4b, indent=2))

if __name__ == "__main__":
    main()
```

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle/bench/mechanism
python3 real_tandem.py --locus <GORILLA_LOCUS> --intact-prefix /home/juanfra/winloci_scratch/mech_real_full
```
Expected: `v4b.recovered=true`, `identity_to_intact>0.95`. If the real locus does not recover but the sim did, that is a real, reportable result about read depth / divergence at that locus — record it honestly; try one alternate Soto tandem before concluding.

- [ ] **Step 4: Commit**

```bash
git add bench/mechanism/real_tandem.py bench/mechanism/verification_results.json
git commit -m "test(mechanism): V4b real gorilla tandem removal-recovery, same-read-set ground truth

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 7: The artifact — `bench/rustle_mechanism.html`, led by §0

**Files:**
- Create: `bench/rustle_mechanism.html`

**Interfaces:**
- Consumes: `bench/mechanism/heuristics.tsv` (embedded as the appendix table), `bench/mechanism/verification_results.json` (V1–V4 numbers + the recovery figure data).

- [ ] **Step 1: Build the page top-to-bottom in the spec's §0 order**

Sections, in this exact order (self-contained: inline CSS/JS, no external assets; theme-aware light/dark; favicon 🧬):
1. **"The method in two sentences."** Define (exon-sum homology id≥0.80, ≥2 loci) + Assign (certificate, p<α/(n−1), else abstain). Verbatim from spec §0.
2. **"The one test that does the heavy lifting."** The certificate written once; the three verified call sites from `verification_results.json.v1.certificate_call_sites`. Frame as the `AS>10` equivalent, stronger because it abstains on ties.
3. **"Two consequences, not new methods."** Genome-absent copies ⇐ Assign (with the V4 recovery figure inline: planted 3 → deleted 1 → recovered, number + set-matched identity, both sim and real gorilla); annotation-absent novel genes ⇐ Define + exon-sum. Each with its one-line honesty caveat (characterization machinery is separate).
4. **Header claim:** "These are two mechanisms and their consequences — not four approaches."
5. **The spine** (subordinate): five stages, each with one mechanism sentence + its primary number + `file:line`.
6. **Superseded-terminology table** (from spec §2.1): k-mer Jaccard/Union-Find, profile-HMM, POA-core-0.13-as-edge, "no ≥X% cutoff" → each mapped to its real stage or "abandoned (date)".
7. **Heuristics appendix:** the full `heuristics.tsv` rendered as a sortable table, grouped by tier, with the framing sentence "Only the `decision` rows pick results; `inert-guard` rows are proven outside the data's range on the Soto families (see V2)."

- [ ] **Step 2: Embed the data (no external fetch — CSP-safe)**

Inline `heuristics.tsv` and `verification_results.json` as `<script type="application/json">` blocks; render tables/figure with vanilla JS. The V4 figure is a simple inline SVG bar/dot showing planted vs recovered copy count and per-copy identity for both arms.

- [ ] **Step 3: Verify it renders and the numbers match**

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
python3 -c "import pathlib,re; h=pathlib.Path('bench/rustle_mechanism.html').read_text(); \
assert 'two sentences' in h.lower(); assert '0.80' in h; assert 'abstain' in h.lower(); \
assert 'not four approaches' in h.lower(); print('structure OK, len', len(h))"
```
Expected: `structure OK`. Open in a browser to eyeball light/dark + the spine + the V4 figure.

- [ ] **Step 4: Commit**

```bash
git add bench/rustle_mechanism.html
git commit -m "feat(mechanism): rustle_mechanism.html — single method led by the two-sentence §0

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 8: The three corrections (false/contradictory artifacts)

**Files:**
- Modify: `bench/quasiclique.html`
- Modify: `bench/DEFINITIONS_FORMAL.md` (the line asserting POA-core-0.13 as the definitional predicate)
- Modify: `figures/01_pipeline_overview.md`, `figures/02_*.md`, `figures/03_*.md`, `figures/04_*.md`
- Modify: `analysis/family_graphs/docs/advisor_response.md`

- [ ] **Step 1: Fix the false claim in `quasiclique.html`**

Find the string asserting "no arbitrary '≥ X% identical' cutoff" and replace with the truth: the family EDGE is asm20 identity ≥0.80, cov-of-shorter ≥0.50 (exon-sum); γ=0.20 is the *clustering* density applied on top of those edges — a different quantity at a different stage.

Run to locate:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle && grep -n "X% identical\|no arbitrary" bench/quasiclique.html
```

- [ ] **Step 2: Fix the contradiction in `DEFINITIONS_FORMAL.md`**

Run: `cd /mnt/c/Users/jfris/Desktop/Rustle && grep -n "core_recip\|0.13\|definitional predicate" bench/DEFINITIONS_FORMAL.md`
Replace the "{u,v} ∈ E_r iff POA contiguous-core ≥ 0.13" definitional claim with: the live family edge is `refine_families_exon_sum` (asm20 id≥0.80, cov≥0.50, ≥2 distinct loci); POA contiguous-core `T_CORE=0.13` is NOT the edge — it survives only as the thin-locus **rescue** admission (`rescue_pipeline`/`family_rescue.rs:41`). Cite that `detect_families`/`core_recip≥0.13` has no caller in any shipped binary.

- [ ] **Step 3: Add superseded banners**

At the very top of each stale file, add a blockquote banner:
```markdown
> **⚠ SUPERSEDED (2026-07-20).** This describes an earlier architecture (k-mer Jaccard + Union-Find / profile-HMM) that the code no longer uses. The current single method is in `bench/rustle_mechanism.html`. Kept for history.
```
Apply to `figures/01_pipeline_overview.md` (+02/03/04) and `analysis/family_graphs/docs/advisor_response.md`.

- [ ] **Step 4: Verify no artifact still contradicts the code**

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
grep -rn "no arbitrary.*identical\|iff POA contiguous-core" bench/quasiclique.html bench/DEFINITIONS_FORMAL.md || echo "contradictions removed"
grep -rl "SUPERSEDED (2026-07-20)" figures/ analysis/family_graphs/docs/ | wc -l
```
Expected: `contradictions removed`; the banner count ≥ 5.

- [ ] **Step 5: Commit**

```bash
git add bench/quasiclique.html bench/DEFINITIONS_FORMAL.md figures/ analysis/family_graphs/docs/advisor_response.md
git commit -m "docs(mechanism): correct false/contradictory artifacts; mark stale ones superseded

- quasiclique.html: remove the false 'no identity cutoff' claim (edge IS asm20 id>=0.80)
- DEFINITIONS_FORMAL.md: 0.13 is thin-locus rescue, not the edge (detect_families is dead)
- figures/01-04 + advisor_response.md: SUPERSEDED banners

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 9: Final consistency gate

**Files:** none created; a verification-only task.

- [ ] **Step 1: Registry still verifies against source**

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
python3 bench/mechanism/gen_heuristics.py --repo-root . --toml bench/mechanism/heuristics.toml --out bench/mechanism/heuristics.tsv
```
Expected: exit 0, `verified N heuristics`.

- [ ] **Step 2: All tests pass**

Run: `cd /mnt/c/Users/jfris/Desktop/Rustle && python3 -m pytest tests/mechanism/ -v`
Expected: all PASS.

- [ ] **Step 3: Every tier claim in the artifact is backed by data**

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
python3 -c "import json; d=json.load(open('bench/mechanism/verification_results.json')); \
print('v1',bool(d.get('v1'))); print('v2',bool(d.get('v2'))); print('v4a',bool(d.get('v4a'))); print('v4b',bool(d.get('v4b')))"
```
Expected: all `True`. Any inert-guard whose V2 `fires=True` must have been reclassified `result-guard` (Task 4 Step 4).

- [ ] **Step 4: Commit any final registry regeneration**

```bash
git add -A bench/mechanism/
git commit -m "chore(mechanism): final consistency gate — registry verified, tests green, tiers backed

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>" --allow-empty
```

---

## Self-Review Notes (for the planner)

- **Spec coverage:** §0 → Task 7 (§0-led HTML) + Task 3 (certificate reuse verified). §2.2 registry/generator → Tasks 1–2. §2.3 tiering → Task 2 + Task 4 (inert proof). §2.4 corrections → Task 8. §3 V1 → Task 3; V2/V3 → Task 4; V4a → Task 5; V4b → Task 6. Success criteria → Task 9. All covered.
- **Placeholder scan:** the two end-to-end drivers (Task 5 Step 5, Task 6 Step 3) intentionally leave the align/parse body as inline comments with exact commands and exact JSON fields to fill — because the precise read-simulation and consensus-parse depend on the recovered output columns, which the implementer confirms at runtime. Every command, flag, path, and expected value is concrete; no threshold or type is left unnamed.
- **Type consistency:** `set_match` (Task 5) is imported and reused verbatim in Task 6. `verification_results.json` keys (`v1`,`v2`,`v3`,`v4a`,`v4b`) are consistent across Tasks 3–7. `classify_guard` signature matches its test and its caller.
