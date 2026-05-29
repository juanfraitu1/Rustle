# StringTie Shadow Mode — Scaffold + Layer 1 (Junction Acceptance) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add the `st_shadow()` master predicate (default OFF) and wire Layer 1 — junction acceptance — so that under `RUSTLE_ST_SHADOW=1` Rustle stops rejecting `mm_negative` junctions that StringTie accepts, then validate that the `junction_accept` mm_negative divergence drives to zero.

**Architecture:** `st_shadow()` is a new predicate in `stringtie_parity.rs`, default off, distinct from the always-on `stringtie_exact()`. Layer 1 OR's `st_shadow()` into the existing `mm_negative` reject gate at `graph_build.rs` so the junction-acceptance behavior matches StringTie only in shadow mode. Validation is a parity-diff of the `junction_accept` JSONL events between Rustle and StringTie (not a unit test) — the gate is "Rustle-rejected ∩ ST-accepted, reason=mm_negative → 0."

**Tech Stack:** Rust (cargo, `target/release/rustle`), the instrumented StringTie fork at `tools/stringtie` (build `make clean release`), Python 3 for the parity-diff, `gffcompare` for the final F1 check (only after all layers — not this increment).

**Reference:** spec `docs/superpowers/specs/2026-05-28-stringtie-shadow-mode-design.md`; findings `docs/STRINGTIE_PARITY_FINDINGS.md`.

---

### Task 1: Add the `st_shadow()` predicate (pure helper + env wrapper)

**Files:**
- Modify: `src/rustle/stringtie_parity.rs` (add after `stringtie_exact()`, ~line 179)

- [ ] **Step 1: Write the failing test**

Add to the bottom of `src/rustle/stringtie_parity.rs` (create a `#[cfg(test)] mod tests` block if none exists):

```rust
#[cfg(test)]
mod tests {
    use super::st_shadow_from;

    #[test]
    fn st_shadow_default_off() {
        assert!(!st_shadow_from(None));        // unset -> OFF (opposite of stringtie_exact)
        assert!(st_shadow_from(Some("1")));    // =1 -> ON
        assert!(st_shadow_from(Some("true"))); // any non-"0" -> ON
        assert!(!st_shadow_from(Some("0")));   // =0 -> OFF
    }
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p rustle st_shadow_default_off 2>&1 | tail -20`
Expected: FAIL — `cannot find function st_shadow_from in this scope` (compile error).

- [ ] **Step 3: Write minimal implementation**

Add after `stringtie_exact()` (after line 179) in `src/rustle/stringtie_parity.rs`:

```rust
/// Pure helper for `st_shadow()` — testable without process env.
#[inline]
pub fn st_shadow_from(v: Option<&str>) -> bool {
    matches!(v, Some(s) if s != "0")
}

/// Master predicate for the coherent end-to-end StringTie-faithful SHADOW mode.
/// Default OFF (distinct from the always-on `stringtie_exact()`). Enable with
/// `RUSTLE_ST_SHADOW=1`. When on, every wired layer behaves like StringTie *together*
/// (junction acceptance, graph, read->transfrag, flow, filter, abundance) so behaviors
/// that regress in isolation reinforce. See docs/STRINGTIE_PARITY_FINDINGS.md.
#[inline]
pub fn st_shadow() -> bool {
    st_shadow_from(std::env::var("RUSTLE_ST_SHADOW").ok().as_deref())
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p rustle st_shadow_default_off 2>&1 | tail -20`
Expected: PASS (`test result: ok. 1 passed`).

- [ ] **Step 5: Commit**

```bash
git add src/rustle/stringtie_parity.rs
git commit -m "feat(shadow): add st_shadow() master predicate (default off)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: Layer 1 — gate `mm_negative` junction acceptance on `st_shadow()`

**Files:**
- Modify: `src/rustle/graph_build.rs` (the `mm_negative` reject in `filter_junctions_for_bundle`, ~line 828-840)

**Context:** The gate currently rejects any junction with `stat.mm < 0.0` (a higherr/long-read demotion marker, NOT absence of support — `junctions.rs:20` warns). StringTie accepts these on raw read support. `RUSTLE_KEEP_MM_NEG=1` already keeps them when `nreads_good>0`. Layer 1 makes `st_shadow()` imply that behavior.

- [ ] **Step 1: Find the exact current code**

Run: `grep -n "keep_mm_neg" src/rustle/graph_build.rs`
Expected: shows the `let keep_mm_neg = std::env::var_os("RUSTLE_KEEP_MM_NEG").is_some();` line inside `filter_junctions_for_bundle`.

- [ ] **Step 2: Make the change**

In `src/rustle/graph_build.rs`, replace:

```rust
                    let keep_mm_neg = std::env::var_os("RUSTLE_KEEP_MM_NEG").is_some();
```

with:

```rust
                    // Layer 1 of the StringTie shadow mode: StringTie accepts mm<0
                    // (demotion-marker) junctions on raw read support; rejecting them
                    // splits reads -> over-segmentation. st_shadow() implies keep.
                    let keep_mm_neg = std::env::var_os("RUSTLE_KEEP_MM_NEG").is_some()
                        || crate::stringtie_parity::st_shadow();
```

- [ ] **Step 3: Build to verify it compiles**

Run: `cargo build --release 2>&1 | grep -E "^error" | head; echo "exit=${PIPESTATUS[0]}"`
Expected: no `error` lines, `exit=0`.

- [ ] **Step 4: Verify default behavior is unchanged (regression guard)**

Run:
```bash
./target/release/rustle GGO_19.bam -L -o /tmp/shadow_off.gtf 2>/dev/null
gffcompare -r ../GGO_19.gtf -R -Q -o /tmp/shadow_off /tmp/shadow_off.gtf 2>/dev/null
grep "Intron chain level" /tmp/shadow_off.stats
```
Expected: `Intron chain level:    96.5     |    90.7` (unchanged — `st_shadow()` is off by default).

- [ ] **Step 5: Commit**

```bash
git add src/rustle/graph_build.rs
git commit -m "feat(shadow): Layer 1 - accept mm_negative junctions under st_shadow()

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: Build the `junction_accept` parity-diff (the Layer 1 gate)

**Files:**
- Create: `bench/junction_accept_diff.py`

**Context:** Both tools emit `junction_accept` JSONL events keyed by `(start, end)` (donor, acceptor+1 — same convention, verified). The Layer-1 gate is: of the junctions StringTie ACCEPTS, how many does Rustle REJECT with reason `mm_negative`. Layer 1 done = that count → 0.

- [ ] **Step 1: Write the diff script**

Create `bench/junction_accept_diff.py`:

```python
#!/usr/bin/env python3
"""Layer-1 gate: compare junction_accept between Rustle and StringTie.

Reports, of the junctions StringTie ACCEPTS, how many Rustle REJECTS and why.
The Layer-1 (mm_negative) gate is satisfied when the 'mm_negative' bucket -> 0.
Inputs: /tmp/ru_ja.jsonl /tmp/st_ja.jsonl (capture with PARITY_FILTER_STEPS=junction_accept).
"""
import json, collections

def accepted_and_rejected(path):
    acc = set()
    rej = {}  # (start,end) -> reason
    for line in open(path):
        try: e = json.loads(line)
        except Exception: continue
        if e.get("step") != "junction_accept": continue
        p = e.get("payload", e)
        k = (e.get("start"), e.get("end"))
        if p.get("accepted") in (True, 1, "true"):
            acc.add(k)
        else:
            rej[k] = str(p.get("reason", "?"))
    return acc, rej

ru_acc, ru_rej = accepted_and_rejected("/tmp/ru_ja.jsonl")
st_acc, _ = accepted_and_rejected("/tmp/st_ja.jsonl")
print(f"Rustle accepted: {len(ru_acc)}  ST accepted: {len(st_acc)}")

# of ST-accepted junctions, which does Rustle reject, by reason
buckets = collections.Counter()
for k in st_acc:
    if k in ru_acc:
        continue
    buckets[ru_rej.get(k, "absent_from_rustle")] += 1
print("=== ST-accepted junctions NOT accepted by Rustle, by Rustle reason ===")
for r, n in buckets.most_common():
    print(f"  {n}\t{r}")
print(f"\nLAYER-1 GATE (mm_negative bucket): {buckets.get('mm_negative', 0)} (target 0)")
```

- [ ] **Step 2: Capture baseline (shadow OFF) junction_accept from both tools**

Run:
```bash
RUSTLE_PARITY_LOG=/tmp/ru_ja.jsonl RUSTLE_PARITY_FILTER_CHROM=NC_073243.2 RUSTLE_PARITY_FILTER_STEPS=junction_accept \
  ./target/release/rustle GGO_19.bam -L -o /dev/null 2>/dev/null
STRINGTIE_PARITY_LOG=/tmp/st_ja.jsonl STRINGTIE_PARITY_FILTER_CHROM=NC_073243.2 STRINGTIE_PARITY_FILTER_STEPS=junction_accept \
  ./tools/stringtie/stringtie GGO_19.bam -L -o /tmp/x.gtf 2>/dev/null
python3 bench/junction_accept_diff.py
```
Expected (shadow off): the `mm_negative` bucket is large (~10000+), confirming the baseline divergence.

- [ ] **Step 3: Commit the diff tool**

```bash
git add bench/junction_accept_diff.py
git commit -m "feat(shadow): junction_accept parity-diff (Layer-1 gate)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: Validate Layer 1 — drive `mm_negative` bucket to 0 under shadow

**Files:** none (validation only)

- [ ] **Step 1: Capture junction_accept with shadow ON**

Run:
```bash
RUSTLE_ST_SHADOW=1 RUSTLE_PARITY_LOG=/tmp/ru_ja.jsonl RUSTLE_PARITY_FILTER_CHROM=NC_073243.2 RUSTLE_PARITY_FILTER_STEPS=junction_accept \
  ./target/release/rustle GGO_19.bam -L -o /dev/null 2>/dev/null
python3 bench/junction_accept_diff.py
```

- [ ] **Step 2: Confirm the gate**

Expected: `LAYER-1 GATE (mm_negative bucket): 0 (target 0)`.
If non-zero: some `mm<0` junctions still rejected — investigate whether they have `nreads_good>0` (the keep condition); if `nreads_good==0` they are genuinely unsupported and ST may reject them too (acceptable residual). Note the residual count.
The `strand_mismatch` bucket is EXPECTED to remain non-zero — that is Layer 2 (graph/strand), not Layer 1.

- [ ] **Step 3: Record the Layer-1 result in the findings doc**

Append a "Shadow Layer 1 — DONE" note to `docs/STRINGTIE_PARITY_FINDINGS.md` §5b with the before/after `mm_negative` bucket counts (baseline ~10000+ → 0) and the remaining `strand_mismatch` count (deferred to Layer 2).

- [ ] **Step 4: Commit**

```bash
git add docs/STRINGTIE_PARITY_FINDINGS.md
git commit -m "docs(shadow): Layer 1 (junction mm_negative acceptance) converged

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Self-Review

- **Spec coverage:** Scaffold (`st_shadow()` default-off, distinct from `stringtie_exact()`) → Task 1. Layer 1 (junction acceptance, `mm_negative`) → Task 2. Per-layer parity-diff gate (`junction_accept` → 0) → Tasks 3–4. "F1 only at end / not mid-stack" → honored (this increment does NOT gate on F1; Task 2 Step 4 only checks the *default* is unchanged, which must hold since `st_shadow()` is off there). Safety (default-off) → Task 1 + Task 2 Step 4.
- **Placeholder scan:** none — all code and commands are concrete.
- **Type consistency:** `st_shadow_from(Option<&str>) -> bool` and `st_shadow() -> bool` used consistently in Task 1 and referenced as `crate::stringtie_parity::st_shadow()` in Task 2.
- **Note:** Tasks 2–4 are validated by parity-diff, not unit tests (the behavior change is graph-level; the diff IS the test). Task 1 is the only classically unit-testable piece and uses TDD.
