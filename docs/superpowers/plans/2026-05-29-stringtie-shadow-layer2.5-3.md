# StringTie Shadow Mode — Layers 2.5 + 3 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Under `RUSTLE_ST_SHADOW=1`, make Rustle's graph node set (Layer 2.5) and read→transfrag construction (Layer 3) StringTie-faithful so `transfrag_pre_depl` Rustle-only chains drive from 5756 → 0, with the default (flag OFF) unchanged at 96.5/90.7.

**Architecture:** Oracle-first: prove node-parity is the lever before porting (cheap correlational oracle, with an abort criterion). Then make Rustle's graph node boundaries match ST in-place in `graph_build.rs` (Layer 2.5), then add a parallel ST-faithful read→transfrag path in `map_reads.rs` (Layer 3). All gated on `crate::stringtie_parity::st_shadow()` (default OFF). Parity-diffs are the gates, NOT unit tests (graph-level behavior); only the `st_shadow()` predicate itself is unit-tested (already done in Layer 1).

**Tech Stack:** Rust (`cargo build --release` → `target/release/rustle`), instrumented StringTie fork at `tools/stringtie` (binary `./tools/stringtie/stringtie`, build `make clean release` ~10s if needed), Python 3 for parity-diffs, the parity-decisions JSONL infra (`RUSTLE_PARITY_LOG` / `STRINGTIE_PARITY_LOG` + `*_PARITY_FILTER_CHROM` / `*_PARITY_FILTER_STEPS`).

**Reference:** spec `docs/superpowers/specs/2026-05-29-stringtie-shadow-layer2.5-3-design.md`; findings `docs/STRINGTIE_PARITY_FINDINGS.md` (§6b, §6c); parent shadow spec `docs/superpowers/specs/2026-05-28-stringtie-shadow-mode-design.md`.

**Key facts the engineer needs:**
- The default pipeline uses `map_reads_to_graph_bundlenodes` (pipeline.rs ~13373), NOT the fallback `map_reads_to_graph`. Instrument/modify the former.
- `st_shadow()` lives at `src/rustle/stringtie_parity.rs:193`, default OFF (`RUSTLE_ST_SHADOW=1` enables). Call it as `crate::stringtie_parity::st_shadow()`.
- Both tools emit `graphnode_list` (Rustle pipeline.rs:13737, ST rlink.cpp:15496) and `bundlenode_list`. Payload: `"n_nodes":N,"nodes":"start-end:cov,start-end:cov,..."`, keyed by `(step, chrom, bundle_start, bundle_end, strand)`. Coordinates are 1-based inclusive (`n.start+1`).
- `bench/transfrag_parity_diff.py` already exists; default compares `transfrag_pre_depl` both sides, reports `Rustle-only` chain count. Current shadow-ON baseline: Rustle 9774 / ST 4383 / Rustle-only 5756.
- The benchmark BAM is `GGO_19.bam` in the repo root; reference `../GGO_19.gtf`; chrom `NC_073243.2`.
- ST source of truth for Layer 3: `get_read_pattern` rlink.cpp:4041, split logic in `get_fragment_pattern` rlink.cpp:4688–4708, end-node trim in `update_abundance` rlink.cpp:4378–4496.
- Rustle Layer-3 sites: `map_reads.rs` `collect_read_nodes_exact`:~250, `split_read_segments`:~1339, single-node-fragment drop:~652, `add_or_update_transfrag`:~1666.

---

### Task 1: `graphnode_list` parity-diff tool + baseline

**Files:**
- Create: `bench/graphnode_diff.py`

**Context:** Layer 2.5's gate is "of the bundles both tools produce, do their graph node sets match?". This tool joins `graphnode_list` events by bundle key `(chrom, start, end, strand)` and reports, per matched bundle, whether the node start-end sets are identical, and aggregate divergence.

- [ ] **Step 1: Inspect the event schema (confirm payload shape)**

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
STRINGTIE_PARITY_LOG=/tmp/st_gn.jsonl STRINGTIE_PARITY_FILTER_CHROM=NC_073243.2 STRINGTIE_PARITY_FILTER_STEPS=graphnode_list ./tools/stringtie/stringtie GGO_19.bam -L -o /tmp/st_x.gtf 2>/dev/null
head -2 /tmp/st_gn.jsonl
```
Expected: JSONL lines with `"step":"graphnode_list"`, a bundle key (`start`/`end`/`strand`/`chrom`), and a payload containing `"nodes":"<s>-<e>:<cov>,..."`. Confirm the field names (`start`,`end`,`strand`, nested `payload` or flat) before writing the parser.

- [ ] **Step 2: Write the diff script**

Create `bench/graphnode_diff.py`:
```python
#!/usr/bin/env python3
"""Layer-2.5 gate: compare graph node sets (graphnode_list) Rustle vs StringTie.

Joins by bundle key (chrom,start,end,strand). For each bundle present in BOTH tools,
compares the set of node (start,end) intervals (coverage ignored). Reports how many
shared bundles have IDENTICAL node sets and the aggregate node-interval divergence.
Layer 2.5 is done when shared-bundle node sets match (mismatch -> 0 / accounted).
Inputs: /tmp/ru_gn.jsonl /tmp/st_gn.jsonl (capture with PARITY_FILTER_STEPS=graphnode_list).
"""
import json, sys, collections, re

def load(path):
    # bundle key (start,end,strand) -> set of (node_start,node_end)
    d = {}
    for line in open(path):
        try: e = json.loads(line)
        except Exception: continue
        if e.get("step") != "graphnode_list": continue
        p = e.get("payload", e)
        nodes_s = p.get("nodes", "")
        key = (e.get("start"), e.get("end"), e.get("strand"))
        intervals = set()
        for tok in nodes_s.split(","):
            tok = tok.strip()
            if not tok: continue
            m = re.match(r"(\d+)-(\d+)", tok)  # "start-end:cov"
            if m: intervals.add((int(m.group(1)), int(m.group(2))))
        d[key] = intervals
    return d

ru = load(sys.argv[1] if len(sys.argv) > 1 else "/tmp/ru_gn.jsonl")
stt = load(sys.argv[2] if len(sys.argv) > 2 else "/tmp/st_gn.jsonl")
shared = set(ru) & set(stt)
print(f"Rustle bundles: {len(ru)}  ST bundles: {len(stt)}  shared: {len(shared)}")

identical = sum(1 for k in shared if ru[k] == stt[k])
ru_extra = sum(len(ru[k] - stt[k]) for k in shared)   # nodes Rustle has, ST doesn't
st_extra = sum(len(stt[k] - ru[k]) for k in shared)    # nodes ST has, Rustle doesn't
print(f"shared bundles with IDENTICAL node sets: {identical}/{len(shared)}")
print(f"node intervals Rustle-extra: {ru_extra}  ST-extra: {st_extra}")
print(f"\nLAYER-2.5 GATE (shared bundles with mismatched node sets): "
      f"{len(shared) - identical} (target 0)")
```

- [ ] **Step 3: Capture baseline (shadow OFF) and run**

Run:
```bash
RUSTLE_PARITY_LOG=/tmp/ru_gn.jsonl RUSTLE_PARITY_FILTER_CHROM=NC_073243.2 RUSTLE_PARITY_FILTER_STEPS=graphnode_list ./target/release/rustle GGO_19.bam -L -o /tmp/ru_x.gtf 2>/dev/null
python3 bench/graphnode_diff.py
```
Expected: a large `LAYER-2.5 GATE` count (many shared bundles with mismatched node sets), confirming the baseline node divergence. Record the numbers.

- [ ] **Step 4: Commit**

```bash
git add bench/graphnode_diff.py
git commit -m "feat(shadow): graphnode_list parity-diff (Layer-2.5 gate)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: Oracle bound — is node-parity the lever? (abort gate)

**Files:**
- Create: `bench/node_parity_oracle.py` (analysis only — no Rustle code change)

**Context:** Before porting anything, prove the decomposition: among bundles where Rustle's graph node set ALREADY matches ST's, is the `transfrag_pre_depl` Rustle-only count ≈ 0? If yes, node-parity is the gating lever (build Layer 2.5 first, as planned). If node-matching bundles STILL have many Rustle-only transfrags, then nodes are not sufficient and the split rule (Layer 3) is independently needed first — replan. This is a cheap correlational oracle (no node-injection needed).

- [ ] **Step 1: Capture graphnode_list AND transfrag_pre_depl together, shadow ON (Layers 1+2)**

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
STRINGTIE_PARITY_LOG=/tmp/st_both.jsonl STRINGTIE_PARITY_FILTER_CHROM=NC_073243.2 STRINGTIE_PARITY_FILTER_STEPS=graphnode_list,transfrag_pre_depl ./tools/stringtie/stringtie GGO_19.bam -L -o /tmp/st_x.gtf 2>/dev/null
RUSTLE_ST_SHADOW=1 RUSTLE_PARITY_LOG=/tmp/ru_both.jsonl RUSTLE_PARITY_FILTER_CHROM=NC_073243.2 RUSTLE_PARITY_FILTER_STEPS=graphnode_list,transfrag_pre_depl ./target/release/rustle GGO_19.bam -L -o /tmp/ru_x.gtf 2>/dev/null
wc -l /tmp/ru_both.jsonl /tmp/st_both.jsonl
```
Expected: both files non-empty, containing both step types.

- [ ] **Step 2: Write the oracle analysis**

Create `bench/node_parity_oracle.py`:
```python
#!/usr/bin/env python3
"""Layer-2.5 ORACLE: does graph-node parity predict transfrag-chain parity?

Partitions bundles into (node-set matches ST) vs (doesn't), then measures the
transfrag_pre_depl Rustle-only rate in each partition. If node-matching bundles have
~0 Rustle-only transfrags, node-parity is the lever -> build Layer 2.5 first.
Inputs: /tmp/ru_both.jsonl /tmp/st_both.jsonl (graphnode_list + transfrag_pre_depl).
"""
import json, re, collections

def load(path):
    gn = {}   # bundle key -> set of node intervals
    tf = collections.defaultdict(set)  # bundle key -> set of intron-chains
    rows = []
    for line in open(path):
        try: e = json.loads(line)
        except Exception: continue
        rows.append(e)
    for e in rows:
        step = e.get("step"); p = e.get("payload", e)
        bkey = (e.get("start"), e.get("end"), e.get("strand"))
        if step == "graphnode_list":
            iv = set()
            for tok in p.get("nodes", "").split(","):
                m = re.match(r"(\d+)-(\d+)", tok.strip())
                if m: iv.add((int(m.group(1)), int(m.group(2))))
            gn[bkey] = iv
    # transfrag chains carry their own coords; assign each to the bundle that contains it
    for e in rows:
        if e.get("step") != "transfrag_pre_depl": continue
        p = e.get("payload", e)
        ich = p.get("introns", "").rstrip(",")
        if not ich: continue
        # bundle containment: transfrag start/end within bundle key
        ts, te = e.get("start"), e.get("end")
        placed = False
        for bkey in gn:
            bs, be, bstr = bkey
            if bs is not None and ts is not None and bs <= ts and te <= be and e.get("strand") == bstr:
                tf[bkey].add((e.get("strand"), ich)); placed = True; break
        if not placed:
            tf[("orphan",)].add((e.get("strand"), ich))
    return gn, tf

ru_gn, ru_tf = load("/tmp/ru_both.jsonl")
st_gn, st_tf = load("/tmp/st_both.jsonl")
shared = set(ru_gn) & set(st_gn)
match_b   = [k for k in shared if ru_gn[k] == st_gn[k]]
mismatch_b = [k for k in shared if ru_gn[k] != st_gn[k]]

def ru_only_chains(bkeys):
    ro = 0
    for k in bkeys:
        ro += len(ru_tf.get(k, set()) - st_tf.get(k, set()))
    return ro

print(f"shared bundles: {len(shared)}  node-MATCH: {len(match_b)}  node-MISMATCH: {len(mismatch_b)}")
print(f"Rustle-only transfrag chains in node-MATCH bundles:   {ru_only_chains(match_b)}")
print(f"Rustle-only transfrag chains in node-MISMATCH bundles:{ru_only_chains(mismatch_b)}")
print("\nORACLE VERDICT: if the node-MATCH Rustle-only count is near 0 while node-MISMATCH")
print("carries the bulk, node-parity is the gating lever -> proceed with Layer 2.5 first.")
```

- [ ] **Step 2b: Run it**

Run: `python3 bench/node_parity_oracle.py`

- [ ] **Step 3: Evaluate the abort gate (DECISION POINT)**

Read the output:
- **PROCEED** if Rustle-only chains in node-MATCH bundles is a small fraction (rule of thumb: <15%) of the total 5756 and the node-MISMATCH partition carries the bulk. This confirms node-parity (Layer 2.5) is the gating lever.
- **REPLAN** if node-MATCH bundles still carry a large share of Rustle-only chains. That means matching nodes is NOT sufficient — the split/trim rule (Layer 3) must come first or the divergence is elsewhere. STOP and report to the human with the partition numbers; do not proceed to Task 3.

Record the verdict and numbers in the commit message.

- [ ] **Step 4: Commit**

```bash
git add bench/node_parity_oracle.py
git commit -m "feat(shadow): Layer-2.5 node-parity oracle + verdict

<paste the node-MATCH / node-MISMATCH Rustle-only numbers and PROCEED/REPLAN verdict here>

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: Characterize the graph-node-boundary divergence

**Files:** none (investigation; output is a findings note)

**Context:** Task 1 shows node sets differ; this task pinpoints WHY so Task 4's change is precise. The change is then localized, not a guess.

- [ ] **Step 1: Find the worst-divergent shared bundles**

Write an ad-hoc python (in `/tmp`, do not commit) that, for the shared bundles from `bench/graphnode_diff.py`, prints the 10 bundles with the largest `|ru_nodes Δ st_nodes|`, showing both node-interval lists side by side. Run it on `/tmp/ru_gn.jsonl` / `/tmp/st_gn.jsonl`.

- [ ] **Step 2: Classify the divergence**

For those bundles, determine the dominant pattern (quantify each):
- **Extra Rustle split nodes** — Rustle splits one ST node into two at an internal boundary (coverage change / trim point) ST doesn't use.
- **Missing Rustle nodes** — ST has a node Rustle merged away.
- **Shifted boundaries** — same node count, off-by-N coordinates (boundary-snap divergence).
Read `src/rustle/graph_build.rs` where graph nodes are created from the bundle (node boundaries derived from junctions + coverage breakpoints). Identify the specific rule producing the dominant divergence and the ST counterpart (rlink.cpp graph construction near the `graphnode_list` emit at 15496, and where ST decides node boundaries).

- [ ] **Step 3: Record the finding**

Append a short "Layer 2.5 node-divergence characterization" note to `docs/STRINGTIE_PARITY_FINDINGS.md` with: the dominant divergence class, counts, the exact `graph_build.rs` site (file:line) of the rule, and the ST counterpart. No commit of code; commit the doc:
```bash
git add docs/STRINGTIE_PARITY_FINDINGS.md
git commit -m "docs(shadow): Layer 2.5 node-boundary divergence characterization

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: Implement Layer 2.5 — graph node parity under `st_shadow()`

**Files:**
- Modify: `src/rustle/graph_build.rs` (the node-boundary rule identified in Task 3)

**Context:** Apply the divergence fix from Task 3, gated on `st_shadow()` so the default is byte-identical. The exact line is Task-3-dependent; the discipline below is fixed. The change makes Rustle's node boundaries match ST's (e.g. suppress the extra Rustle split, or add ST's boundary), ONLY under shadow.

- [ ] **Step 1: Apply the gated change**

At the `graph_build.rs` site from Task 3, branch the node-boundary decision on `crate::stringtie_parity::st_shadow()`. Pattern (adapt to the actual rule — e.g. an extra-split condition `should_split`):
```rust
// Layer 2.5 shadow: match StringTie's graph node boundaries (no extra Rustle
// split / use ST's boundary). Default path unchanged.
let should_split = <existing Rustle condition>
    && !crate::stringtie_parity::st_shadow();   // or the ST-faithful condition under shadow
```
Keep the change minimal and confined to the boundary decision identified in Task 3.

- [ ] **Step 2: Build**

Run: `cargo build --release 2>&1 | grep -E "^error" | head; echo "exit=${PIPESTATUS[0]}"`
Expected: no `error` lines, `exit=0`.

- [ ] **Step 3: Regression guard — default unchanged**

Run:
```bash
./target/release/rustle GGO_19.bam -L -o /tmp/shadow_off.gtf 2>/dev/null
gffcompare -r ../GGO_19.gtf -R -Q -o /tmp/shadow_off /tmp/shadow_off.gtf 2>/dev/null
grep "Intron chain level" /tmp/shadow_off.stats
```
Expected: `Intron chain level:    96.5     |    90.7    |` (unchanged). If different, STOP — the gate leaked into the default.

- [ ] **Step 4: Drive the Layer-2.5 gate**

Run:
```bash
RUSTLE_ST_SHADOW=1 RUSTLE_PARITY_LOG=/tmp/ru_gn.jsonl RUSTLE_PARITY_FILTER_CHROM=NC_073243.2 RUSTLE_PARITY_FILTER_STEPS=graphnode_list ./target/release/rustle GGO_19.bam -L -o /tmp/ru_x.gtf 2>/dev/null
python3 bench/graphnode_diff.py
```
Expected: `LAYER-2.5 GATE (shared bundles with mismatched node sets)` drops substantially toward 0 vs the Task-1 baseline. Record before/after. If a residual remains, note whether it is a different divergence class (may need a follow-up within this task) or genuinely irreducible (coordinate artifact) — document the residual count.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/graph_build.rs
git commit -m "feat(shadow): Layer 2.5 - match StringTie graph node boundaries under st_shadow()

graphnode_list shared-bundle mismatch <BEFORE> -> <AFTER>. Default unchanged 96.5/90.7.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 5: Implement Layer 3 (part A) — ST-faithful read-segment split

**Files:**
- Modify: `src/rustle/map_reads.rs` (`split_read_segments` ~1339, and the single-node-fragment drop ~652)

**Context:** With nodes matching (Task 4), make the read→node-path split match ST. ST splits a read's node path ONLY at consecutive non-contiguous nodes where NO read junction exactly equals the node boundary (`get_fragment_pattern`, rlink.cpp:4688–4708). Rustle currently also splits at killed/BADJUNC junctions and drops single-node leftover segments — both removed under shadow.

- [ ] **Step 1: Read the current split logic**

Read `src/rustle/map_reads.rs` `split_read_segments` (~1339–1551) and the single-node drop (~652). Identify (a) the killed/BADJUNC split condition, (b) the unitig source/sink split, (c) the single-node-fragment discard.

- [ ] **Step 2: Add the ST-faithful split branch under shadow**

In `split_read_segments`, gate the killed/BADJUNC split condition so it does NOT fire under shadow; instead split only at non-contiguous node pairs lacking a boundary-matching read junction (ST's rule). Concretely, where Rustle currently decides to split at a junction, wrap:
```rust
// Layer 3 shadow: StringTie (get_fragment_pattern, rlink.cpp:4688-4708) splits a
// read's node path ONLY between consecutive non-contiguous nodes with NO read
// junction matching the node boundary. It has no killed/BADJUNC split. Under
// shadow, suppress the Rustle-only split and use ST's non-contiguous-node rule.
let split_here = if crate::stringtie_parity::st_shadow() {
    // non-contiguous nodes AND no read junction equals (prev_end, curr_start)
    nodes_non_contiguous && !read_junction_matches_boundary
} else {
    <existing Rustle split condition>
};
```
(Use the actual variable names found in Step 1 for `nodes_non_contiguous` / the boundary-match test; ST's match is `junc.start == prev_node.end && junc.end == curr_node.start`.)

- [ ] **Step 3: Stop dropping single-node leftover segments under shadow**

At the single-node-fragment drop (~652), gate it so under shadow the fragment is NOT dropped (ST never creates it because it didn't split). Pattern:
```rust
if is_single_node_segment && !crate::stringtie_parity::st_shadow() {
    continue; // legacy: drop orphan single-node segment
}
```

- [ ] **Step 4: Build + regression guard**

Run:
```bash
cargo build --release 2>&1 | grep -E "^error" | head; echo "exit=${PIPESTATUS[0]}"
./target/release/rustle GGO_19.bam -L -o /tmp/shadow_off.gtf 2>/dev/null
gffcompare -r ../GGO_19.gtf -R -Q -o /tmp/shadow_off /tmp/shadow_off.gtf 2>/dev/null
grep "Intron chain level" /tmp/shadow_off.stats
```
Expected: `exit=0`; Intron chain level `96.5 | 90.7` unchanged.

- [ ] **Step 5: Measure the Layer-3 gate (partial)**

Run:
```bash
RUSTLE_ST_SHADOW=1 RUSTLE_PARITY_LOG=/tmp/ru_tf.jsonl RUSTLE_PARITY_FILTER_CHROM=NC_073243.2 RUSTLE_PARITY_FILTER_STEPS=transfrag_pre_depl ./target/release/rustle GGO_19.bam -L -o /tmp/ru_x.gtf 2>/dev/null
STRINGTIE_PARITY_LOG=/tmp/st_tf.jsonl STRINGTIE_PARITY_FILTER_CHROM=NC_073243.2 STRINGTIE_PARITY_FILTER_STEPS=transfrag_pre_depl ./tools/stringtie/stringtie GGO_19.bam -L -o /tmp/st_x.gtf 2>/dev/null
python3 bench/transfrag_parity_diff.py
```
Expected: Rustle-only chains drop from 5756 toward ST's 4383 (the truncation/subset chains collapse onto full chains). Record. Some residual remains for Task 6 (end-node trim). NOTE: if Rustle-only RISES, the split rule fired without node-parity from Task 4 holding — re-check Task 4 landed.

- [ ] **Step 6: Commit**

```bash
git add src/rustle/map_reads.rs
git commit -m "feat(shadow): Layer 3a - ST-faithful read-segment split under st_shadow()

transfrag_pre_depl Rustle-only 5756 -> <AFTER>. Default unchanged 96.5/90.7.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 6: Implement Layer 3 (part B) — end-node trim (canonicalize partial reads)

**Files:**
- Modify: `src/rustle/map_reads.rs` (`add_or_update_transfrag` ~1666 / the existing `trim_longread_path_for_update_abundance` ~1709)

**Context:** ST's `update_abundance` (rlink.cpp:4378–4496) trims a long read's spanning END nodes back toward source/sink so a partial read's pattern collapses onto the full chain (`findtrf_in_treepat`). Rustle has a `trim_longread_path_for_update_abundance` already — align it to ST's trim under shadow so remaining superset/partial chains merge.

- [ ] **Step 1: Read both trims**

Read Rustle's `trim_longread_path_for_update_abundance` (~1709) and ST's `update_abundance` trim (rlink.cpp:4378–4496). Identify the difference in which end nodes are trimmed and the match/collapse key (`findtrf_in_treepat`).

- [ ] **Step 2: Gate the ST-faithful trim under shadow**

Where Rustle's trim differs, branch on `crate::stringtie_parity::st_shadow()` to apply ST's end-node trim (trim spanning first/last nodes back to source/sink per rlink.cpp:4378–4496), then the existing node-path key matching collapses identical patterns. Keep the default path unchanged. (Use the actual node-path variable names from Step 1; comment with the rlink.cpp reference.)

- [ ] **Step 3: Build + regression guard**

Run:
```bash
cargo build --release 2>&1 | grep -E "^error" | head; echo "exit=${PIPESTATUS[0]}"
./target/release/rustle GGO_19.bam -L -o /tmp/shadow_off.gtf 2>/dev/null
gffcompare -r ../GGO_19.gtf -R -Q -o /tmp/shadow_off /tmp/shadow_off.gtf 2>/dev/null
grep "Intron chain level" /tmp/shadow_off.stats
```
Expected: `exit=0`; `96.5 | 90.7` unchanged.

- [ ] **Step 4: Measure the Layer-3 gate**

Run:
```bash
RUSTLE_ST_SHADOW=1 RUSTLE_PARITY_LOG=/tmp/ru_tf.jsonl RUSTLE_PARITY_FILTER_CHROM=NC_073243.2 RUSTLE_PARITY_FILTER_STEPS=transfrag_pre_depl ./target/release/rustle GGO_19.bam -L -o /tmp/ru_x.gtf 2>/dev/null
python3 bench/transfrag_parity_diff.py
```
Expected: Rustle-only chains approach 0 (or a small accounted residual). Record before/after.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/map_reads.rs
git commit -m "feat(shadow): Layer 3b - ST end-node trim canonicalizes partial reads under st_shadow()

transfrag_pre_depl Rustle-only <BEFORE> -> <AFTER>. Default unchanged 96.5/90.7.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 7: Validate combined Layers 2.5+3 and document

**Files:** Modify `docs/STRINGTIE_PARITY_FINDINGS.md`

- [ ] **Step 1: Capture both gates with shadow ON**

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
RUSTLE_ST_SHADOW=1 RUSTLE_PARITY_LOG=/tmp/ru_both.jsonl RUSTLE_PARITY_FILTER_CHROM=NC_073243.2 RUSTLE_PARITY_FILTER_STEPS=graphnode_list,transfrag_pre_depl ./target/release/rustle GGO_19.bam -L -o /tmp/ru_x.gtf 2>/dev/null
STRINGTIE_PARITY_LOG=/tmp/st_both.jsonl STRINGTIE_PARITY_FILTER_CHROM=NC_073243.2 STRINGTIE_PARITY_FILTER_STEPS=graphnode_list,transfrag_pre_depl ./tools/stringtie/stringtie GGO_19.bam -L -o /tmp/st_x.gtf 2>/dev/null
# split the combined logs for the existing diff tools
grep '"step":"graphnode_list"' /tmp/ru_both.jsonl > /tmp/ru_gn.jsonl; grep '"step":"graphnode_list"' /tmp/st_both.jsonl > /tmp/st_gn.jsonl
grep '"step":"transfrag_pre_depl"' /tmp/ru_both.jsonl > /tmp/ru_tf.jsonl; grep '"step":"transfrag_pre_depl"' /tmp/st_both.jsonl > /tmp/st_tf.jsonl
python3 bench/graphnode_diff.py
python3 bench/transfrag_parity_diff.py
```

- [ ] **Step 2: Confirm both gates**

Expected: graphnode mismatch ≈ 0; transfrag_pre_depl Rustle-only ≈ 0 (or a documented residual). Confirm Layer 1 still holds (optional re-run of `bench/junction_accept_diff.py`).

- [ ] **Step 3: Verify default once more**

Run:
```bash
./target/release/rustle GGO_19.bam -L -o /tmp/shadow_off.gtf 2>/dev/null
gffcompare -r ../GGO_19.gtf -R -Q -o /tmp/shadow_off /tmp/shadow_off.gtf 2>/dev/null
grep "Intron chain level" /tmp/shadow_off.stats
```
Expected: `96.5 | 90.7`.

- [ ] **Step 4: Record + commit**

Append a "Shadow Layers 2.5+3 — DONE" section to `docs/STRINGTIE_PARITY_FINDINGS.md` with: graphnode mismatch before→after, transfrag Rustle-only 5756→after, the residual breakdown (if any), and the note that F1 is still deferred to end-of-project. Then:
```bash
git add docs/STRINGTIE_PARITY_FINDINGS.md
git commit -m "docs(shadow): Layers 2.5+3 (graph-node + read->transfrag) converged

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Self-Review

- **Spec coverage:** Layer 2.5 (graph node parity) → Tasks 1,3,4 (gate graphnode_diff). Layer 3 (split + trim) → Tasks 5,6 (gate transfrag_pre_depl Rustle-only). Oracle-first with abort criterion → Task 2. Dispatch architecture (2.5 in-place graph_build.rs; 3 under shadow in map_reads.rs) → Tasks 4,5,6. Parity-diffs as gates not unit tests → all tasks (only `st_shadow()` already unit-tested). Default-OFF safety → regression guard in Tasks 4,5,6,7. F1 not mid-stack → honored (no F1 step). Validate+document → Task 7.
- **Placeholder scan:** The Task-4 code is intentionally parameterized on the Task-3 finding (a research dependency, not a placeholder) — its harness/gate/commit steps are concrete, and the source-of-truth (ST line refs) is given. Tasks 5–6 give concrete gated code patterns keyed to named ST source. `<BEFORE>/<AFTER>` in commit messages are runtime values to fill at commit time, not plan gaps.
- **Type consistency:** `st_shadow()` referenced as `crate::stringtie_parity::st_shadow()` throughout; the two diff tools (`graphnode_diff.py`, `node_parity_oracle.py`) and the existing `transfrag_parity_diff.py` use the same JSONL key convention `(start,end,strand)` and payload `"nodes"`/`"introns"` fields.
- **Note:** Task 2's abort gate can terminate the plan early (REPLAN) — that is by design (cheap kill before the port).
