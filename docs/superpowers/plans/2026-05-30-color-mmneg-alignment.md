# StringTie-faithful color + mm_negative alignment — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Measure — oracle-first, default-OFF — whether aligning both halves of bad-junction handling to StringTie converges Rustle's junction set, colors, and bundle segmentation to ST's, and decide ship/shelve from the measured prize before building Phase 1.

**Architecture:** Phase 0 is pure analysis — junction-set parity instrumentation (env-gated emits in both tools) + two new bench diff tools (`junction_set_diff.py`, `segmentation_prize.py`) — gated by a hard abort gate (`< 5 net FP/FN`). Phase 1 (only if the gate clears) wires a single default-OFF flag `RUSTLE_ST_BADJUNC` that OR-composes both halves into existing flag paths: Half 1 = the `RUSTLE_SUBBUNDLE_USE_STATS` junction-stats path (ST's `strand==0` color break); Half 2 = the `keep_mm_neg` path (stop the extra mm_negative edge-drop). No strand mutation anywhere.

**Tech Stack:** Rust (Rustle pipeline), C++ (instrumented StringTie at `tools/stringtie/`), Python 3 stdlib (bench diff tools), gffcompare (benchmark).

---

## Design correction from grounding (READ FIRST)

The spec (`2026-05-30-color-mmneg-alignment-design.md` §4 Half 2) said "at `graph_build.rs:845`, set `JunctionStat.strand=0` instead of edge-filtering." **This is infeasible and is corrected here.** A grounding investigation (adversarial verifier, confirmed) established:

- `graph_build.rs` `filter_junctions_for_bundle` receives `junction_stats: Option<&JunctionStats>` — an **immutable** reference. It only *reads* `stat.strand`/`stat.mm`; it cannot set strand. Mutating it breaks the read-only filtering contract and multiple call sites.
- `JunctionStat.strand = Some(0)` is set **upfront in `killed_junctions.rs`** (lines 656, 676, 755, 796, 818, 837, 871) — this is Rustle's good_junc-equivalent (ST `good_junc`, `rlink.cpp:13765-13868`, strand-zeroing at 13771/13782/13797/13827/13847/13862).
- `graph_build.rs:843-849` *separately* edge-drops `mm < 0.0` junctions via `keep_mm_neg` (already gated by `RUSTLE_KEEP_MM_NEG` || `st_shadow()`). This mm_negative edge-drop is Rustle's EXTRA aggression that ST's good_junc does NOT have as a standalone filter.

**Therefore the correct Half 2** is to OR `st_badjunc()` into the existing `keep_mm_neg` predicate at `graph_build.rs:843` — i.e. *stop* the extra mm_negative edge-drop when the flag is on (matching ST, which keeps those edges unless its good_junc witness/long-intron gates fire and set strand=0 in `killed_junctions.rs`). Half 1 and Half 2 then both become one-line OR-compositions with existing flag paths.

**Caveat for Phase 1 (record, do not pre-solve):** keeping mm_negative edges means full ST faithfulness also depends on `killed_junctions.rs` having ST's good_junc witness/long-intron strand-zeroing logic. If Rustle's `killed_junctions.rs` lacks parts of that witness logic, some mm_negative junctions ST would demote (strand=0) will stay ±1 in Rustle — this is the known mm_negative floor (findings §6h). Phase 1 measures the gap; it does not attempt to re-port the witness logic (out of scope).

---

## File structure

**Phase 0 (analysis — no production behavior change):**
- Modify `src/rustle/pipeline.rs` — add two env-gated parity emits: `junction_raw` (after `compute_initial_junction_stats_for_reads`, ~line 1508) and `good_junction` (at good_junctions finalization, ~line 12096, next to the existing `cgroup_feed` dump).
- Modify `tools/stringtie/rlink.cpp` — add matching `junction_raw` emit (after raw junction build, ~line 503, end of `processRead`'s junction loop) and `good_junction` emit (in/after `good_junc`, ~line 13868).
- Create `bench/junction_set_diff.py` — three-layer junction-set classifier.
- Create `bench/segmentation_prize.py` — prize-bound at the divergent loci (per-bundle node-set breakdown + reference FP/FN attribution).
- Create `bench/capture_parity.sh` — canonical capture harness for both tools + gffcompare.

**Phase 1 (production, behind `RUSTLE_ST_BADJUNC`, default OFF — only if gate clears):**
- Modify `src/rustle/stringtie_parity.rs` — add `st_badjunc_from` / `st_badjunc` predicates + unit test.
- Modify `src/rustle/pipeline.rs:14815` — OR `st_badjunc()` into `use_stats_for_subbundle` (Half 1).
- Modify `src/rustle/graph_build.rs:843` — OR `st_badjunc()` into `keep_mm_neg` (Half 2).

---

## PHASE 0 — Analysis & abort gate

### Task 1: Capture harness for both tools

**Files:**
- Create: `bench/capture_parity.sh`

- [ ] **Step 1: Write the capture script**

```bash
#!/bin/bash
# Capture parity events from Rustle + instrumented StringTie on GGO_19 (chr19),
# and gffcompare Rustle's output vs the reference.
# Usage: bash bench/capture_parity.sh
set -euo pipefail
cd /mnt/c/Users/jfris/Desktop/Rustle

CHROM=NC_073243.2
BAM=GGO_19.bam                       # cwd copy (also at ../GGO_19.bam)
REF=/mnt/c/Users/jfris/Desktop/GGO_19.gtf
STEPS=graphnode_list,junction_accept,junction_raw,good_junction

# 1) Rustle de novo + parity log
RUSTLE_PARITY_LOG=/tmp/ru_parity.jsonl \
RUSTLE_PARITY_FILTER_STEPS=$STEPS \
RUSTLE_PARITY_FILTER_CHROM=$CHROM \
  ./target/release/rustle -L "$BAM" -o /tmp/ru_final.gtf 2>/dev/null

# 2) Instrumented StringTie + parity log
STRINGTIE_PARITY_LOG=/tmp/st_parity.jsonl \
STRINGTIE_PARITY_FILTER_STEPS=$STEPS \
STRINGTIE_PARITY_FILTER_CHROM=$CHROM \
  ./tools/stringtie/stringtie -L "$BAM" -o /tmp/st_final.gtf 2>/dev/null

# 3) Split per-step for the diff tools (graphnode_list lines into _gn files)
grep '"step":"graphnode_list"' /tmp/ru_parity.jsonl > /tmp/ru_gn.jsonl || true
grep '"step":"graphnode_list"' /tmp/st_parity.jsonl > /tmp/st_gn.jsonl || true

# 4) gffcompare Rustle vs reference (headline)
gffcompare -r "$REF" -o /tmp/ru_cmp /tmp/ru_final.gtf 2>/dev/null
echo "=== Rustle vs reference ==="
grep -E "Intron chain level:|Transcript level:" /tmp/ru_cmp.stats
echo "(parity logs: /tmp/ru_parity.jsonl /tmp/st_parity.jsonl; finals: /tmp/ru_final.gtf /tmp/st_final.gtf)"
```

- [ ] **Step 2: Make executable and dry-run (Rustle side only first)**

Run: `chmod +x bench/capture_parity.sh && bash bench/capture_parity.sh`
Expected: prints "Intron chain level: 96.2 | 91.7" and "Transcript level: 95.6 | 90.5" (±0.1pp run-to-run nondeterminism). If the StringTie binary at `tools/stringtie/stringtie` does not exist, build it first: `( cd tools/stringtie && make clean release )` (~10s).

- [ ] **Step 3: Verify both parity logs are non-empty and carry the existing steps**

Run: `wc -l /tmp/ru_gn.jsonl /tmp/st_gn.jsonl && head -1 /tmp/ru_parity.jsonl`
Expected: both `_gn` files have hundreds of lines (≈3400 bundles each); first line is a JSON object with `"step"` and `"payload"`. The `junction_raw`/`good_junction` steps will be EMPTY until Tasks 2–3 add them — that is expected at this point.

- [ ] **Step 4: Commit**

```bash
git add bench/capture_parity.sh
git commit -m "bench: capture harness for color+mmneg parity (Rustle+ST+gffcompare)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: Add `junction_raw` + `good_junction` parity emits to Rustle

**Files:**
- Modify: `src/rustle/pipeline.rs` (~1508 raw; ~12096 good_junctions finalization)

The emit API (grounded): `crate::parity::decisions::emit(step, chrom, start, end, strand, payload_json)` with the `pjson!` macro; gated by `crate::parity::decisions::is_enabled()`. Existing `junction_accept` emits `start=j.donor`, `end=j.acceptor+1` (1-based to match ST), payload `{accepted, bundle_strand, jstrand, mm, reason}` (`graph_build.rs:856-882`). Match that coordinate convention.

- [ ] **Step 1: Emit `junction_raw` after raw junction stats are computed**

In `pipeline.rs`, immediately after `compute_initial_junction_stats_for_reads(...)` returns its `(JunctionStats, _)` (around line 1508), add (use the actual returned binding name for the stats map — read the surrounding code to confirm it; assume `junction_stats`):

```rust
// Parity: emit every RAW observed junction (pre-filter) for cross-tool junction-set diff.
if crate::parity::decisions::is_enabled() {
    for (j, st) in junction_stats.iter() {
        let strand_c = match st.strand {
            Some(1) => '+',
            Some(-1) => '-',
            _ => '.',
        };
        crate::parity::decisions::emit(
            "junction_raw",
            None,
            j.donor,
            j.acceptor + 1,
            strand_c,
            &crate::pjson!(
                "jstrand": st.strand.map(|s| s as i64).unwrap_or(-9),
                "mm": format!("{:.4}", st.mm),
                "nreads_good": format!("{:.2}", st.nreads_good),
                "mrcount": format!("{:.2}", st.mrcount)
            ),
        );
    }
}
```

- [ ] **Step 2: Emit `good_junction` at good_junctions finalization**

In `pipeline.rs` near line 12096 (right where the finalized `good_junctions_set` is logged via the existing `cgroup_feed`/`junction_dump` path, BEFORE it is passed to `build_bundlenodes_...` at ~12375), add:

```rust
// Parity: emit the FINAL good_junctions membership set (what has_good_left queries).
if crate::parity::decisions::is_enabled() {
    for j in good_junctions_set.iter() {
        crate::parity::decisions::emit(
            "good_junction",
            None,
            j.donor,
            j.acceptor + 1,
            '.',
            "",
        );
    }
}
```

- [ ] **Step 3: Build**

Run: `cargo build --release 2>&1 | tail -3`
Expected: `Finished release` with no errors. (If the `good_junctions_set` binding name differs at line 12096, read 11980–12099 and use the actual name.)

- [ ] **Step 4: Verify the new steps now emit**

Run: `bash bench/capture_parity.sh >/dev/null 2>&1; grep -c '"step":"junction_raw"' /tmp/ru_parity.jsonl; grep -c '"step":"good_junction"' /tmp/ru_parity.jsonl`
Expected: `junction_raw` count ≈ a few thousand (all observed junctions on chr19); `good_junction` count ≈ 7000–7500 (matches the known 7303 accepted). Both > 0.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/pipeline.rs
git commit -m "parity(rustle): emit junction_raw + good_junction for junction-set diff

Env-gated (RUSTLE_PARITY_LOG); no production behavior change.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: Add matching `junction_raw` + `good_junction` emits to StringTie

**Files:**
- Modify: `tools/stringtie/rlink.cpp` (~503 raw build; ~13868 good_junc)

ST emit API (grounded): `pd_emit(const char* step, const char* chrom, unsigned long start, unsigned long end, char strand, const char* payload_json)`, gated by `pd_enabled()` (`parity_decisions.h:19-26`, `parity_decisions.cc:86-112`). Existing `junction_accept` emit at `rlink.cpp:16896` uses `j->start`, `j->end` (ST is 1-based inclusive — donor=start, acceptor=end; this aligns with Rustle's `donor` / `acceptor+1`). The raw junction vector is built in `processRead` (`AddIfNew` at line 494); the full junction list is `junction` (a `GList<CJunction>`).

- [ ] **Step 1: Emit `junction_raw` for the full raw junction list**

ST builds junctions incrementally per read, so emit once over the complete `junction` list rather than per-read. Find where `count_good_junctions` / `good_junc` is invoked over the bundle (near the existing `junction_accept` loop at `rlink.cpp:16876-16901`); immediately BEFORE the good_junc filtering runs, iterate the raw `junction` list and emit. Read 16850–16905 to locate the loop variable; mirror its `cname`/strand handling:

```cpp
// Parity: emit every RAW junction (pre good_junc) to match Rustle's junction_raw.
if (pd_enabled()) {
  for (int i = 0; i < junction.Count(); i++) {
    CJunction* j = junction[i];
    char strand_c = (j->strand > 0) ? '+' : (j->strand < 0 ? '-' : '.');
    char payload[256];
    snprintf(payload, sizeof(payload),
      "\"jstrand\":%d,\"nreads_good\":%.2f,\"nreads\":%.2f",
      (int)j->strand, (double)j->nreads_good, (double)j->nreads);
    pd_emit("junction_raw", cname, (unsigned long)j->start, (unsigned long)j->end, strand_c, payload);
  }
}
```

- [ ] **Step 2: Emit `good_junction` after good_junc filtering**

After `good_junc` has run (so `j->strand` reflects demotions; the strand-zeroing sites are 13771–13862), in the same post-filter region (near 16876–16901), emit only junctions ST keeps (strand != 0):

```cpp
// Parity: emit the kept (strand!=0) junctions = ST's good-junction membership.
if (pd_enabled()) {
  for (int i = 0; i < junction.Count(); i++) {
    CJunction* j = junction[i];
    if (j->strand == 0) continue;
    char strand_c = (j->strand > 0) ? '+' : '-';
    pd_emit("good_junction", cname, (unsigned long)j->start, (unsigned long)j->end, strand_c, "");
  }
}
```

- [ ] **Step 3: Build StringTie**

Run: `( cd tools/stringtie && make clean release 2>&1 | tail -3 )`
Expected: `stringtie` binary rebuilt, no errors (~10s).

- [ ] **Step 4: Verify ST now emits both steps**

Run: `bash bench/capture_parity.sh >/dev/null 2>&1; grep -c '"step":"junction_raw"' /tmp/st_parity.jsonl; grep -c '"step":"good_junction"' /tmp/st_parity.jsonl`
Expected: both > 0; `good_junction` ≈ 17000+ (ST keeps far more — the known 17477 accepted).

- [ ] **Step 5: Commit**

```bash
git add tools/stringtie/rlink.cpp
git commit -m "parity(stringtie): emit junction_raw + good_junction to match Rustle

Env-gated (STRINGTIE_PARITY_LOG); no behavior change.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: `bench/junction_set_diff.py` — three-layer junction-set classifier

**Files:**
- Create: `bench/junction_set_diff.py`

Template: `bench/junction_accept_diff.py` (reason-bucket pattern) + `bench/flow_oracle.py` `load_accepted` (junction jsonl load). All junctions keyed by `(start, end)` — both tools now emit on the same 1-based convention (Rustle `donor`/`acceptor+1`, ST `start`/`end`).

- [ ] **Step 1: Write the classifier**

```python
#!/usr/bin/env python3
"""Three-layer junction-set parity: raw / accepted / good_junction.
Answers: do we have exactly ST's junctions? only-RU / only-ST / shared-same / shared-divergent.
Usage: python3 bench/junction_set_diff.py /tmp/ru_parity.jsonl /tmp/st_parity.jsonl"""
import json, sys
from collections import defaultdict

def load(path):
    layers = {"junction_raw": {}, "junction_accept": {}, "good_junction": set()}
    with open(path) as fh:
        for line in fh:
            try:
                r = json.loads(line)
            except ValueError:
                continue
            step = r.get("step")
            if step not in layers:
                continue
            key = (r["start"], r["end"])
            if step == "good_junction":
                layers[step].add(key)
            else:
                layers[step][key] = r.get("payload", {})
    return layers

def report(name, ru_set, st_set):
    only_ru = ru_set - st_set
    only_st = st_set - ru_set
    shared = ru_set & st_set
    print(f"  [{name}] shared={len(shared)}  only-RU={len(only_ru)}  only-ST={len(only_st)}")
    return only_ru, only_st, shared

def main():
    ru = load(sys.argv[1] if len(sys.argv) > 1 else "/tmp/ru_parity.jsonl")
    st = load(sys.argv[2] if len(sys.argv) > 2 else "/tmp/st_parity.jsonl")

    print("=== Layer 1: RAW observed junctions ===")
    report("raw", set(ru["junction_raw"]), set(st["junction_raw"]))

    print("=== Layer 2: ACCEPTED (junction_accept, accepted=true) ===")
    ru_acc = {k for k, p in ru["junction_accept"].items() if p.get("accepted")}
    st_acc = {k for k, p in st["junction_accept"].items() if p.get("accepted")}
    report("accepted", ru_acc, st_acc)

    print("=== Layer 3: good_junction membership ===")
    only_ru_g, only_st_g, shared_g = report("good", ru["good_junction"], st["good_junction"])

    print("=== Shared-divergent treatment (in raw of both, kept by ST, dropped by RU good) ===")
    raw_both = set(ru["junction_raw"]) & set(st["junction_raw"])
    divergent = [k for k in raw_both
                 if k in st["good_junction"] and k not in ru["good_junction"]]
    # bucket the divergent ones by RU's reason (from junction_accept) / mm sign
    buckets = defaultdict(int)
    for k in divergent:
        p = ru["junction_accept"].get(k)
        if p is None:
            buckets["not_in_ru_accept"] += 1
        else:
            reason = p.get("reason", "?")
            mmneg = "mm<0" if float(p.get("mm", 0)) < 0 else "mm>=0"
            buckets[f"{reason}|{mmneg}"] += 1
    print(f"  ST-keeps-but-RU-drops: {len(divergent)}")
    for b, n in sorted(buckets.items(), key=lambda x: -x[1]):
        print(f"    {b}: {n}")

if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run it on captured logs**

Run: `bash bench/capture_parity.sh >/dev/null 2>&1; python3 bench/junction_set_diff.py /tmp/ru_parity.jsonl /tmp/st_parity.jsonl`
Expected: Layer-2 accepted shows `only-RU≈0` (we know RU 7303 ⊆ ST 17477); Layer-3 good shows ST keeps thousands RU drops; the divergent bucket breaks those down by reason (`mm_negative|mm<0`, `strand_zero`, etc.). This IS the deliverable — no pass/fail.

- [ ] **Step 3: Commit**

```bash
git add bench/junction_set_diff.py
git commit -m "bench: junction_set_diff.py (3-layer RU vs ST junction parity)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 5: `bench/segmentation_prize.py` — prize-bound at divergent loci

**Files:**
- Create: `bench/segmentation_prize.py`

This is the abort-gate measurement. It needs per-bundle node sets (graphnode_diff only prints aggregates — we re-load `graphnode_list` directly, as graphnode_diff does at lines 12-28) AND reference-based FP/FN attribution (reuse the gffcompare `.tracking`/`.refmap` from Task 1, or recompute chains vs reference like `selection_oracle.py:81-90`). Approach: (1) find bundles where RU's node set ≠ ST's and RU has ≥2 same-strand bundles spanning ST's single envelope (true splits); (2) for each, count RU final transcripts inside that envelope whose intron chain is a reference FP (`class_code` u/p/j/etc., absent from reference), attributable to the split.

- [ ] **Step 1: Write the prize-bound tool**

```python
#!/usr/bin/env python3
"""Prize-bound for color/segmentation convergence: net reachable FP+FN at the
boundary-mismatch bundles. Abort gate: if net < 5, shelve Phase 1.
Usage: python3 bench/segmentation_prize.py /tmp/ru_gn.jsonl /tmp/st_gn.jsonl \
           /tmp/ru_final.gtf /mnt/c/Users/jfris/Desktop/GGO_19.gtf"""
import json, sys, re
from collections import defaultdict

def load_gn(path):
    bundles = []
    with open(path) as fh:
        for line in fh:
            try:
                r = json.loads(line)
            except ValueError:
                continue
            if r.get("step") != "graphnode_list":
                continue
            nodes = tuple(tuple(map(int, m.split("-")))
                          for m in re.findall(r"\d+-\d+", r.get("payload", {}).get("nodes", "")))
            bundles.append((r["start"], r["end"], r["strand"], nodes))
    return bundles

def load_gtf_chains(path):
    """transcript_id -> (strand, (start,end), intron_chain) ; intron_chain=() for single-exon."""
    tx = defaultdict(list)
    strand = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            m = re.search(r'transcript_id "([^"]+)"', f[8])
            if not m:
                continue
            tid = m.group(1)
            tx[tid].append((int(f[3]), int(f[4])))
            strand[tid] = f[6]
    chains = {}
    for tid, ex in tx.items():
        ex.sort()
        introns = tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))
        chains[tid] = (strand[tid], (ex[0][0], ex[-1][1]), introns)
    return chains

def overlaps(a0, a1, b0, b1):
    return a0 <= b1 and b0 <= a1

def main():
    ru_gn = load_gn(sys.argv[1] if len(sys.argv) > 1 else "/tmp/ru_gn.jsonl")
    st_gn = load_gn(sys.argv[2] if len(sys.argv) > 2 else "/tmp/st_gn.jsonl")
    ru_gtf = load_gtf_chains(sys.argv[3] if len(sys.argv) > 3 else "/tmp/ru_final.gtf")
    ref = load_gtf_chains(sys.argv[4] if len(sys.argv) > 4 else "/mnt/c/Users/jfris/Desktop/GGO_19.gtf")
    ref_chains = {(s, ic) for (s, _sp, ic) in ref.values() if ic}

    # Find ST envelopes split by RU (>=2 same-strand RU bundles spanning one ST bundle).
    splits = []
    for (s0, s1, sstrand, snodes) in st_gn:
        ru_inside = [(r0, r1) for (r0, r1, rstrand, _rn) in ru_gn
                     if rstrand == sstrand and r0 >= s0 - 1 and r1 <= s1 + 1]
        # true split: >=2 RU bundles, none spanning the whole ST envelope
        if len(ru_inside) >= 2 and not any(r0 <= s0 + 1 and r1 >= s1 - 1 for (r0, r1) in ru_inside):
            splits.append((s0, s1, sstrand))

    # Attribute FP: RU final chains inside a split envelope that are NOT in the reference.
    fp = 0
    fp_examples = []
    for tid, (strand, (t0, t1), ic) in ru_gtf.items():
        if not ic:
            continue
        if (strand, ic) in ref_chains:
            continue  # a TP, not attributable
        for (s0, s1, sstrand) in splits:
            if sstrand == strand and overlaps(t0, t1, s0, s1):
                fp += 1
                fp_examples.append((tid, strand, t0, t1))
                break

    # FN: reference chains overlapping a split envelope that RU does NOT produce.
    ru_chains = {(s, ic) for (s, _sp, ic) in ru_gtf.values() if ic}
    fn = 0
    for (s, (r0, r1), ic) in ref.values():
        if not ic or (s, ic) in ru_chains:
            continue
        for (s0, s1, sstrand) in splits:
            if sstrand == s and overlaps(r0, r1, s0, s1):
                fn += 1
                break

    print(f"True-split ST envelopes (RU fragments): {len(splits)}")
    print(f"Attributable FP (RU non-ref chains in split envelopes): {fp}")
    print(f"Attributable FN (ref chains RU misses in split envelopes): {fn}")
    print(f"NET reachable prize (FP+FN): {fp + fn}")
    print(f"ABORT GATE: {'ABORT (shelve Phase 1)' if fp + fn < 5 else 'PROCEED to Phase 1'}")
    for tid, strand, t0, t1 in fp_examples[:10]:
        print(f"    FP {tid} {strand} {t0}-{t1}")

if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run the prize-bound**

Run: `python3 bench/segmentation_prize.py /tmp/ru_gn.jsonl /tmp/st_gn.jsonl /tmp/ru_final.gtf /mnt/c/Users/jfris/Desktop/GGO_19.gtf`
Expected: prints true-split count (≈11 from prior analysis), attributable FP (≈2–3: RSTL.103.1 class p, RSTL.331.2 class j), FN (≈0–2), NET prize, and the ABORT-GATE verdict.

- [ ] **Step 3: Cross-check the FP examples against the earlier trace**

Run: confirm the printed FP transcript IDs include the previously-traced split-locus FPs (the class-p / class-j cases near loci `31437809-31559882` and `20532505-20568524`). If `segmentation_prize.py` reports 0 splits, the node-set load/strand keys are wrong — debug against `bench/graphnode_diff.py` output (it should report ≈125 mismatched bundles) before trusting the prize number.

- [ ] **Step 4: Commit**

```bash
git add bench/segmentation_prize.py
git commit -m "bench: segmentation_prize.py (abort-gate prize-bound at split loci)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 6: Run Phase 0, record the abort-gate decision

**Files:**
- Modify: `docs/STRINGTIE_PARITY_FINDINGS.md` (extend §6j)
- Modify: memory `project_color_cgroup_parity.md`

- [ ] **Step 1: Run the full Phase-0 capture + both diff tools**

```bash
bash bench/capture_parity.sh
python3 bench/junction_set_diff.py /tmp/ru_parity.jsonl /tmp/st_parity.jsonl
python3 bench/segmentation_prize.py /tmp/ru_gn.jsonl /tmp/st_gn.jsonl /tmp/ru_final.gtf /mnt/c/Users/jfris/Desktop/GGO_19.gtf
```

- [ ] **Step 2: Record the junction-set picture + the prize number**

Append to `docs/STRINGTIE_PARITY_FINDINGS.md` §6j: the three-layer junction classification (only-RU/only-ST/shared-same/shared-divergent with reason buckets) and the prize-bound (true-split count, attributable FP+FN, net). Update memory `project_color_cgroup_parity.md` with the measured numbers.

- [ ] **Step 3: ABORT-GATE decision**

If NET prize `< 5` (~0.25pp): **STOP — Phase 1 shelved.** Record "color divergence bounded as structurally negligible; junction-set parity answered" as the successful outcome. Do NOT proceed to Phase 1. Report to the user.

If NET prize `>= 5`: proceed to Phase 1 (Tasks 7–10).

- [ ] **Step 4: Commit the findings**

```bash
git add docs/STRINGTIE_PARITY_FINDINGS.md
git commit -m "docs: Phase-0 results — junction-set parity + segmentation prize-bound

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## PHASE 1 — Two-halves alignment (ONLY if abort gate cleared)

### Task 7: Add the `RUSTLE_ST_BADJUNC` flag predicate

**Files:**
- Modify: `src/rustle/stringtie_parity.rs` (predicates + test, follow the `st_flow` template at lines 213/217/240-253)

- [ ] **Step 1: Write the failing test**

Add to the `#[cfg(test)]` module in `stringtie_parity.rs`:

```rust
#[test]
fn st_badjunc_default_off() {
    use super::st_badjunc_from;
    assert!(!st_badjunc_from(None));
    assert!(st_badjunc_from(Some("1")));
    assert!(!st_badjunc_from(Some("0")));
}
```

- [ ] **Step 2: Run it to verify it fails**

Run: `cargo test -p rustle st_badjunc_default_off 2>&1 | tail -5`
Expected: FAIL — `cannot find function st_badjunc_from`.

- [ ] **Step 3: Implement the predicate**

Add next to `st_flow` (around line 213):

```rust
pub fn st_badjunc_from(v: Option<&str>) -> bool {
    matches!(v, Some(s) if s != "0")
}

pub fn st_badjunc() -> bool {
    st_badjunc_from(std::env::var("RUSTLE_ST_BADJUNC").ok().as_deref())
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cargo test -p rustle st_badjunc_default_off 2>&1 | tail -5`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/stringtie_parity.rs
git commit -m "feat(parity): RUSTLE_ST_BADJUNC predicate (default OFF)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 8: Half 1 — OR `st_badjunc()` into the junction-stats path

**Files:**
- Modify: `src/rustle/pipeline.rs:14815`

- [ ] **Step 1: OR the flag into `use_stats_for_subbundle`**

Change line 14815 from:

```rust
let use_stats_for_subbundle = std::env::var_os("RUSTLE_SUBBUNDLE_USE_STATS").is_some();
```

to:

```rust
let use_stats_for_subbundle = std::env::var_os("RUSTLE_SUBBUNDLE_USE_STATS").is_some()
    || crate::stringtie_parity::st_badjunc();
```

- [ ] **Step 2: Build + confirm default unchanged**

Run: `cargo build --release 2>&1 | tail -2 && bash bench/capture_parity.sh`
Expected: builds; default headline still `95.6 / 90.5` (flag OFF → no change).

- [ ] **Step 3: Commit**

```bash
git add src/rustle/pipeline.rs
git commit -m "feat(parity): half-1 — st_badjunc routes color break to ST strand==0

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 9: Half 2 — OR `st_badjunc()` into `keep_mm_neg`

**Files:**
- Modify: `src/rustle/graph_build.rs:843`

- [ ] **Step 1: OR the flag into `keep_mm_neg`**

Change the `keep_mm_neg` definition (line 843) from:

```rust
let keep_mm_neg = std::env::var_os("RUSTLE_KEEP_MM_NEG").is_some()
    || crate::stringtie_parity::st_shadow();
```

to:

```rust
let keep_mm_neg = std::env::var_os("RUSTLE_KEEP_MM_NEG").is_some()
    || crate::stringtie_parity::st_shadow()
    || crate::stringtie_parity::st_badjunc();
```

- [ ] **Step 2: Build + confirm default unchanged**

Run: `cargo build --release 2>&1 | tail -2 && bash bench/capture_parity.sh`
Expected: builds; default headline still `95.6 / 90.5` (flag OFF).

- [ ] **Step 3: Commit**

```bash
git add src/rustle/graph_build.rs
git commit -m "feat(parity): half-2 — st_badjunc keeps mm_negative edges (no extra drop)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 10: Combined parity gate + F1 measurement + decision

**Files:**
- Modify: `docs/STRINGTIE_PARITY_FINDINGS.md`, memory `project_color_cgroup_parity.md`

- [ ] **Step 1: Run with the flag ON and capture parity + F1**

```bash
RUSTLE_ST_BADJUNC=1 RUSTLE_PARITY_LOG=/tmp/ru_bj.jsonl \
RUSTLE_PARITY_FILTER_STEPS=graphnode_list,junction_accept,good_junction \
RUSTLE_PARITY_FILTER_CHROM=NC_073243.2 \
  ./target/release/rustle -L GGO_19.bam -o /tmp/ru_bj.gtf 2>/dev/null
grep '"step":"graphnode_list"' /tmp/ru_bj.jsonl > /tmp/ru_bj_gn.jsonl
python3 bench/graphnode_diff.py /tmp/ru_bj_gn.jsonl /tmp/st_gn.jsonl
gffcompare -r /mnt/c/Users/jfris/Desktop/GGO_19.gtf -o /tmp/ru_bj_cmp /tmp/ru_bj.gtf 2>/dev/null
grep -E "Intron chain level:|Transcript level:" /tmp/ru_bj_cmp.stats
```

- [ ] **Step 2: Evaluate the convergence + F1**

Success criterion (per spec §6): `graphnode_diff` mismatched-bundle count should drop toward 0 (from 125) and the junction divergent-treatment count (re-run `junction_set_diff.py` with `/tmp/ru_bj.jsonl`) should drop. Record the F1 delta vs the `95.6 / 90.5` baseline. Transient regression is accepted — record its magnitude as the cost of exact-ST color parity.

- [ ] **Step 3: Default decision**

If flag-ON F1 holds or improves AND segmentation converged: recommend default-flip to the user (do NOT flip without explicit approval — the default operating point must not change without sign-off). Otherwise: keep `RUSTLE_ST_BADJUNC` as a documented ST-faithful opt-in; record the regression magnitude and the residual mm_negative-floor gap (the good_junc witness logic not present in `killed_junctions.rs`).

- [ ] **Step 4: Commit findings + update memory**

```bash
git add docs/STRINGTIE_PARITY_FINDINGS.md
git commit -m "docs: Phase-1 combined st_badjunc result (convergence + F1 + decision)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Self-review notes

- **Spec coverage:** Phase 0a junction-set parity → Tasks 2–4; Phase 0b prize-bound → Task 5; abort gate → Task 6 Step 3; Phase 1 two-halves → Tasks 7–9; combined gate + F1 + default decision → Task 10. The spec's infeasible Half-2 (graph_build strand mutation) is corrected in the "Design correction" section and Task 9 (keep_mm_neg OR) — recorded, not silently changed.
- **Default-unchanged guard:** Tasks 8 & 9 Step 2 each re-verify `95.6 / 90.5` with the flag OFF.
- **Type consistency:** `st_badjunc_from`/`st_badjunc` match the `st_flow` template; emit steps `junction_raw`/`good_junction` used consistently across Tasks 2/3/4; `keep_mm_neg`/`use_stats_for_subbundle` are the real existing bindings (grounded).
- **No placeholders:** every code/command step shows actual code/commands. Two bindings flagged to confirm against surrounding code at edit time (the stats-map name at `pipeline.rs:1508`, the `good_junctions_set` name at `~12096`, and ST's loop variable at `rlink.cpp:16876`) — these are read-and-confirm instructions, not deferred work.
