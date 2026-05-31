# Prediction-Population Parity (2A donor-snap + 2B CSR-suppression) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Measure — oracle-first, default-OFF — whether suppressing RU's chimeric-suffix-rescue (2B) and re-enabling a narrowed junction donor-snap (2A) remove over-enumeration FPs at acceptable TP-cost; decide ship/shelve each from its oracle.

**Architecture:** Two independent, sequenced sub-levers, both acting on the prediction population. 2B first (it has a free probe `RUSTLE_DISABLE_CSR`). Each: oracle Phase 0 → hard abort gate → default-OFF Phase-1 fix. Gffcompare + parity-diffs are the gates. The retained-intron lesson is enforced: an oracle prize is a ceiling; verify the FP disappears AND the read folds into a surviving TP.

**Tech Stack:** Rust (`path_extract.rs`, `junction_correction.rs`, `stringtie_parity.rs`), Python 3 stdlib (bench), gffcompare.

---

## Grounding facts (verified)

**2B — chimeric_suffix_rescue:**
- `is_chimeric_suffix_rescue(tf_nodes, tmatch, kept_paths, tf_longstart, graph) -> bool` at `path_extract.rs:2992-3033`. Fires when `tf_nodes` is a suffix of a longer kept path (`kept_paths[kidx]`, a 4-tuple `(nodes, cov, guide, out_idx)`). Call site `path_extract.rs:9485` gated by `RUSTLE_DISABLE_CSR`.
- The matched containing path is known via `tmatch` → `kept_paths[kidx]`; its kept/emitted status is testable via `out_idx` and the boundary `flow_kept_paths_len = kept_paths.len()` captured at `path_extract.rs:9257` (indices `< flow_kept_paths_len` = full-length FLOW paths; `>=` = checktrf-rescue paths). So "contained in a kept full-length flow path" = the matched path index `< flow_kept_paths_len` AND its `out_idx` valid (`< out.len()`). No new global state.
- RU does NOT emit `rescue_class` in parity JSON, BUT the GTF carries `source "checktrf_rescue"` + `rescue_class "chimeric_suffix_rescue"` (`path_extract.rs:10555/10568`) — `csr_classify.py` reads those from `/tmp/ri_ru.gtf`.

**2A — junction higherr donor-snap:**
- `is_bad` gate `junction_correction.rs:178-189`: `(s.nm>0 && s.mrcount>0 && s.nm>=s.mrcount) && (intron_len > 100_000 && s.nreads_good < 10.0)`. Narrowed (comment 164-175) from ST's blanket `nm>=nreads` which "killed small exons."
- Nearby-stronger search `junction_correction.rs:215-255`: within `window` (25bp), same strand, `is_reliable(cand)`, `cand.leftsupport > cur.leftsupport*tolerance(0.9)`, `ok_to_demote`. Redirect wired `junctions.rs:626` → read-repair `pipeline.rs:11443`.
- Signals in scope (`JunctionStat`, `types.rs:522-560`): `nm, mrcount, nreads_good, leftsupport, rightsupport, strand, mm, guide_match, consleft, consright` (consensus: -1 unknown/0 non-canonical/1 canonical), `donor, acceptor`.
- `junction_accept` parity payload (`graph_build.rs:863-881`): `accepted, bundle_strand, jstrand, mm, reason`; coords `j.donor`, `j.acceptor+1`. (Lacks nreads_good/leftsupport/consleft — `donor_snap_prize.py` derives dominance from RU final-chain coverage + reference; enhance the emit only if needed.)
- Default 95.6/90.5; baseline /tmp/ri_ru.gtf, /tmp/ri_st.jsonl from earlier captures (re-capture via `bench/capture_parity.sh` if stale).

---

## File structure
- Create `bench/csr_classify.py`, `bench/donor_snap_prize.py`.
- Modify `src/rustle/stringtie_parity.rs` — `st_csr_fold` + `st_higherr_snap` predicates.
- Modify `src/rustle/path_extract.rs:9485` (2B suppression guard).
- Modify `src/rustle/junction_correction.rs:178-189` (2A narrowed re-enable).
- Modify `docs/STRINGTIE_PARITY_FINDINGS.md` + memory.

---

## SUB-LEVER 2B — chimeric-suffix-rescue suppression (FIRST)

### Task 1: 2B free probe (whole-CSR net)

**Files:** none (existing toggle).

- [ ] **Step 1: Capture baseline + DISABLE_CSR**
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
bash bench/capture_parity.sh >/dev/null 2>&1   # refreshes /tmp/ri_ru.gtf etc. via STEPS; or reuse /tmp/ri_ru.gtf
./target/release/rustle -L ../GGO_19.bam -o /tmp/csr_base.gtf 2>/dev/null
RUSTLE_DISABLE_CSR=1 ./target/release/rustle -L ../GGO_19.bam -o /tmp/csr_off.gtf 2>/dev/null
gffcompare -r ../GGO_19.gtf /tmp/csr_base.gtf -o /tmp/csr_base_cmp 2>/dev/null
gffcompare -r ../GGO_19.gtf /tmp/csr_off.gtf -o /tmp/csr_off_cmp 2>/dev/null
echo "=== baseline ==="; grep -E "Intron chain level:|Transcript level:" /tmp/csr_base_cmp.stats
echo "=== CSR disabled ==="; grep -E "Intron chain level:|Transcript level:" /tmp/csr_off_cmp.stats
python3 bench/gtf_chain_diff.py /tmp/csr_base.gtf ../GGO_19.gtf 2>/dev/null | tail -4   # baseline RU-only FP
python3 bench/gtf_chain_diff.py /tmp/csr_off.gtf ../GGO_19.gtf 2>/dev/null | tail -4    # CSR-off RU-only FP
```
Expected: CSR-off should have fewer total tx; net = (FP removed) − (TP lost). Record both. This is the WHOLE-CSR net (upper bound on what suppression can do).

- [ ] **Step 2: Record probe result (no commit — informs the gate)**

Note: baseline vs CSR-off FP count, TP count (in-both chains), Sn/Pr. If CSR-off is strongly net-NEGATIVE (loses more real alt-TSS TPs than FPs removed), the rescue earns its keep → only the targeted (flow-parent) subset is suppressible (Task 2 measures it). If net-positive, even the blunt disable helps.

---

### Task 2: `bench/csr_classify.py` — targeted prize/cost

**Files:** Create `bench/csr_classify.py`

- [ ] **Step 1: Write the classifier**

```python
#!/usr/bin/env python3
"""Classify chimeric_suffix_rescue (CSR) predictions: FP vs TP vs reference.
CSR tx are identified by GTF source="checktrf_rescue" or rescue_class="chimeric_suffix_rescue".
Usage: python3 bench/csr_classify.py /tmp/csr_base.gtf /mnt/c/Users/jfris/Desktop/GGO_19.gtf"""
import sys, re
from collections import defaultdict

def load(path, want_csr_only=False):
    tx_ex = defaultdict(list); strand = {}; is_csr = {}
    for line in open(path):
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9: continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m: continue
        tid = m.group(1)
        if f[2] == "transcript" or f[2] == "exon":
            csr = ('checktrf_rescue' in f[1]) or ('chimeric_suffix_rescue' in f[8])
            is_csr[tid] = is_csr.get(tid, False) or csr
        if f[2] == "exon":
            tx_ex[tid].append((int(f[3]), int(f[4]))); strand[tid] = f[6]
    chains = {}
    for t, ex in tx_ex.items():
        ex.sort()
        ic = tuple((ex[i-1][1]+1, ex[i][0]-1) for i in range(1, len(ex)))
        if ic: chains[t] = (strand[t], ic, is_csr.get(t, False))
    return chains

def main():
    ru = load(sys.argv[1] if len(sys.argv) > 1 else "/tmp/csr_base.gtf")
    ref = load(sys.argv[2] if len(sys.argv) > 2 else "/mnt/c/Users/jfris/Desktop/GGO_19.gtf")
    refset = {(s, ic) for (s, ic, _) in ref.values()}
    csr = [(t, s, ic) for t, (s, ic, c) in ru.items() if c]
    fp = [(t, s, ic) for (t, s, ic) in csr if (s, ic) not in refset]
    tp = [(t, s, ic) for (t, s, ic) in csr if (s, ic) in refset]
    print(f"CSR predictions (multi-intron): {len(csr)}")
    print(f"  FP (not in ref) — suppressible prize: {len(fp)}")
    print(f"  TP (in ref) — cost if suppressed:     {len(tp)}")
    net = len(fp) - len(tp)
    print(f"NET (CSR-FP - CSR-TP): {net}")
    print(f"ABORT GATE: {'ABORT (shelve 2B)' if net <= 0 else ('WEAK (<3)' if net < 3 else 'PROCEED')}")
    for t, s, ic in tp[:8]:
        print(f"    CSR-TP {t} {s} {ic[:2]}...")

if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run + interpret**

Run: `python3 bench/csr_classify.py /tmp/csr_base.gtf ../GGO_19.gtf`
Expected: CSR FP (prize) vs CSR TP (cost). NET = FP − TP. The CSR-TP are real alt-TSS isoforms the rescue correctly recovers (the cost). NOTE: this counts ALL CSR; the Phase-1 guard suppresses the FLOW-PARENT subset (which ST folds) — likely a subset of these. If NET here is clearly positive, the targeted guard is at least as good; if NET ≤ 0, the rescue mostly recovers real isoforms → 2B abort.

- [ ] **Step 3: Commit**
```bash
git add bench/csr_classify.py
git commit -m "bench: csr_classify.py (chimeric-suffix-rescue FP/TP prize-bound)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: 2B abort-gate decision + record

**Files:** Modify `docs/STRINGTIE_PARITY_FINDINGS.md`; memory.

- [ ] **Step 1: Read the gate.** From Task 1 (whole-CSR net) + Task 2 (CSR FP/TP). ABORT 2B if CSR-FP − CSR-TP ≤ 0 (the rescue recovers ≥ as many real isoforms as FPs it makes). PROCEED only if net ≥ ~3 AND the free probe (Task 1) shows disabling is not strongly net-negative.
- [ ] **Step 2: Record** the whole-CSR net + CSR FP/TP + decision. If ABORT: record "CSR rescue earns its keep (recovers N real alt-TSS); suppression net ≤ 0"; STOP 2B (proceed to 2A).
- [ ] **Step 3: Commit** the findings update.

---

### Task 4: `RUSTLE_CSR_FOLD` flag predicate (only if 2B gate clears)

**Files:** Modify `src/rustle/stringtie_parity.rs`

- [ ] **Step 1: Failing test** (`#[cfg(test)]`):
```rust
#[test]
fn st_csr_fold_default_off() {
    use super::st_csr_fold_from;
    assert!(!st_csr_fold_from(None));
    assert!(st_csr_fold_from(Some("1")));
    assert!(!st_csr_fold_from(Some("0")));
}
```
- [ ] **Step 2: Run → FAIL** `cargo test -p rustle st_csr_fold_default_off --lib 2>&1 | tail -3`
- [ ] **Step 3: Implement** (next to `st_flow`):
```rust
pub fn st_csr_fold_from(v: Option<&str>) -> bool { matches!(v, Some(s) if s != "0") }
pub fn st_csr_fold() -> bool { st_csr_fold_from(std::env::var("RUSTLE_CSR_FOLD").ok().as_deref()) }
```
- [ ] **Step 4: Run → PASS**; `cargo build --release 2>&1 | tail -2` clean.
- [ ] **Step 5: Commit** `git add src/rustle/stringtie_parity.rs && git commit -m "feat(parity): RUSTLE_CSR_FOLD predicate (default OFF)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"`

---

### Task 5: 2B suppression guard (only if 2B gate clears)

**Files:** Modify `src/rustle/path_extract.rs:9485`

- [ ] **Step 1: Add the flow-parent containment guard at the call site**

At `path_extract.rs:9485`, the call is `if std::env::var_os("RUSTLE_DISABLE_CSR").is_none() && is_chimeric_suffix_rescue(&tf_nodes, &tmatch, &kept_paths, transfrags[t].longstart, graph) { ... }`. Under `RUSTLE_CSR_FOLD`, additionally SUPPRESS (skip the rescue) when the matched containing path is a kept full-length FLOW path. Add, computed before the `if`:

```rust
// RUSTLE_CSR_FOLD: ST folds a suffix into its full-length flow parent rather than re-extracting it.
// Suppress the rescue when the matched containing path is a kept full-length FLOW path
// (index < flow_kept_paths_len) that was actually emitted (out_idx valid).
let csr_fold_suppress = crate::stringtie_parity::st_csr_fold() && tmatch.iter().any(|&kidx| {
    kidx < flow_kept_paths_len
        && kept_paths.get(kidx).map(|kp| kp.3 < out.len()).unwrap_or(false)
});
```
Then change the guard to `if std::env::var_os("RUSTLE_DISABLE_CSR").is_none() && !csr_fold_suppress && is_chimeric_suffix_rescue(...) {`.
(Confirm `flow_kept_paths_len`, `out`, and the `kept_paths` tuple field `.3 == out_idx` names at edit time — grounding cites `path_extract.rs:9257`, `:9214`.)

- [ ] **Step 2: Build + measure flag ON vs default**
```bash
cargo build --release 2>&1 | tail -2
RUSTLE_CSR_FOLD=1 ./target/release/rustle -L ../GGO_19.bam -o /tmp/csr_fold.gtf 2>/dev/null
gffcompare -r ../GGO_19.gtf /tmp/csr_fold.gtf -o /tmp/csr_fold_cmp 2>/dev/null
echo "=== CSR_FOLD ON ==="; grep -E "Intron chain level:|Transcript level:" /tmp/csr_fold_cmp.stats
python3 bench/gtf_chain_diff.py /tmp/csr_fold.gtf ../GGO_19.gtf 2>/dev/null | tail -4
./target/release/rustle -L ../GGO_19.bam -o /tmp/csr_def.gtf 2>/dev/null
gffcompare -r ../GGO_19.gtf /tmp/csr_def.gtf -o /tmp/csr_def_cmp 2>/dev/null
echo "=== default OFF ==="; grep -E "Transcript level:" /tmp/csr_def_cmp.stats
```
Expected: flag-ON RU-only FP drops by the targeted prize, Pr up, Sn drop ≤ the contained-TP count; flag-OFF MUST be 96.2/91.7, 95.6/90.5.

- [ ] **Step 3: Realizability cross-check + decision.** Confirm the removed FPs are CSR predictions (the oracle-predicted set) and Sn loss matches the contained-TP estimate (no surprise collateral). If net F1 ≥ baseline → recommend default-flip (do NOT flip without approval); else keep `RUSTLE_CSR_FOLD` opt-in + record cost. Commit `src/rustle/path_extract.rs` + findings.

---

## SUB-LEVER 2A — junction donor-snap (SECOND)

### Task 6: `bench/donor_snap_prize.py` — near-dup-donor prize/cost

**Files:** Create `bench/donor_snap_prize.py`

- [ ] **Step 1: Write the prize tool** (reference-based; finds RU final chains whose intron's donor is within 25bp of a stronger same-strand donor used by a higher-coverage RU chain; prize = weak intron ∉ ref, cost = weak intron ∈ ref):

```python
#!/usr/bin/env python3
"""2A prize: RU final chains using a weak donor within 25bp of a stronger same-strand donor.
prize = such chains whose weak intron is NOT a reference intron (FP the snap would remove);
cost  = such chains whose weak intron IS a reference intron (real alt-donor the snap would destroy).
Usage: python3 bench/donor_snap_prize.py /tmp/csr_base.gtf /mnt/c/Users/jfris/Desktop/GGO_19.gtf"""
import sys, re
from collections import defaultdict

def load(path):
    tx_ex = defaultdict(list); strand = {}; cov = {}
    for line in open(path):
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9: continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m: continue
        tid = m.group(1)
        if f[2] == "transcript":
            cm = re.search(r'cov "([0-9.]+)"', f[8]); cov[tid] = float(cm.group(1)) if cm else 0.0
        if f[2] == "exon":
            tx_ex[tid].append((int(f[3]), int(f[4]))); strand[tid] = f[6]
    out = {}
    for t, ex in tx_ex.items():
        ex.sort()
        introns = [(ex[i-1][1]+1, ex[i][0]-1) for i in range(1, len(ex))]
        out[t] = (strand[t], introns, cov.get(t, 0.0))
    return out

def main():
    ru = load(sys.argv[1] if len(sys.argv) > 1 else "/tmp/csr_base.gtf")
    ref = load(sys.argv[2] if len(sys.argv) > 2 else "/mnt/c/Users/jfris/Desktop/GGO_19.gtf")
    ref_introns = set()
    for (s, introns, _c) in ref.values():
        for (d, a) in introns:
            ref_introns.add((s, d, a))
    # collect all RU introns with the max coverage of any chain using them (donor dominance proxy)
    donor_cov = defaultdict(float)   # (strand, donor) -> max cov of a chain using that donor
    intron_cov = defaultdict(float)
    for (s, introns, c) in ru.values():
        for (d, a) in introns:
            donor_cov[(s, d)] = max(donor_cov[(s, d)], c)
            intron_cov[(s, d, a)] = max(intron_cov[(s, d, a)], c)
    WIN = 25
    prize_chains = set(); cost_chains = set()
    for t, (s, introns, c) in ru.items():
        for (d, a) in introns:
            # is there a STRONGER same-strand donor within 25bp?
            stronger = any(donor_cov[(s, d2)] > intron_cov[(s, d, a)] * 1.1
                           for d2 in range(d - WIN, d + WIN + 1) if (s, d2) in donor_cov and d2 != d)
            if not stronger:
                continue
            if (s, d, a) in ref_introns:
                cost_chains.add(t)   # weak donor is a REAL reference intron
            else:
                prize_chains.add(t)  # weak donor FP
    refset = {(s, tuple(i)) for (s, i, _c) in ((s, [(d,a) for (d,a) in ii], c) for (s, ii, c) in ref.values())}
    # restrict prize to chains that are RU-only FP overall (not in ref)
    ru_fp = {t for t, (s, ii, c) in ru.items() if (s, tuple(ii)) not in refset and ii}
    prize = prize_chains & ru_fp
    cost = cost_chains
    print(f"chains with a weak donor near a stronger same-strand donor (<=25bp): prize-side {len(prize_chains)}, cost-side {len(cost_chains)}")
    print(f"PRIZE (RU-only-FP w/ weak near-dup donor): {len(prize)}")
    print(f"COST  (chains whose weak intron IS a real reference intron): {len(cost)}")
    net = len(prize) - len(cost)
    print(f"NET (prize - cost): {net}")
    print(f"ABORT GATE: {'ABORT (cost>=prize)' if len(cost) >= len(prize) else ('WEAK (<3)' if net < 3 else 'PROCEED')}")

if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run + interpret**

Run: `python3 bench/donor_snap_prize.py /tmp/csr_base.gtf ../GGO_19.gtf`
Expected: prize (FP removable by snapping) vs cost (real alt-donors within 25bp the snap would destroy = the small-exon regression). If cost ≥ prize → ABORT 2A (the snap destroys as many real alt-donors as FPs). Record. NOTE: this is a coverage-proxy bound; the actual narrowed gate (canonical+guide+ratio) is STRICTER (spares canonical/guided), so the realized cost should be ≤ this — but verify in Phase 1.

- [ ] **Step 3: Commit** `git add bench/donor_snap_prize.py && git commit -m "bench: donor_snap_prize.py (2A near-dup-donor prize/cost)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"`

---

### Task 7: 2A abort-gate decision + record

- [ ] **Step 1: Read the gate.** ABORT 2A if cost ≥ prize or net < ~3. The narrowed gate (Task 9) will be stricter than this proxy, so a borderline-positive proxy can still proceed IF the canonical/guide gates plausibly remove the cost (note which cost chains are canonical/guided — those would be spared).
- [ ] **Step 2: Record** prize/cost/net + decision (findings + memory). If ABORT, STOP 2A.

---

### Task 8: `RUSTLE_HIGHERR_SNAP` flag predicate (only if 2A gate clears)

**Files:** Modify `src/rustle/stringtie_parity.rs` (same TDD pattern as Task 4)

- [ ] **Step 1-5:** Add `st_higherr_snap_from`/`st_higherr_snap` (`RUSTLE_HIGHERR_SNAP`) with the `st_csr_fold_default_off`-style test; fail → implement → pass → build → commit. (Repeat the Task-4 pattern with `higherr_snap`/`RUSTLE_HIGHERR_SNAP`.)

---

### Task 9: 2A narrowed `is_bad` re-enable (only if 2A gate clears)

**Files:** Modify `src/rustle/junction_correction.rs:178-189`

- [ ] **Step 1: Widen `is_bad` under the flag, with canonical+guide+ratio gates**

Replace the `is_bad` closure body's final line so that, under `RUSTLE_HIGHERR_SNAP`, a junction is ALSO bad (snap-eligible) when it's a NOISE donor — non-canonical, unguided, low-coverage — but the existing long-intron rule and default behavior are unchanged when the flag is off:

```rust
let is_bad = |s: &JunctionStat, j: &Junction| {
    if s.nm <= 0.0 || s.mrcount <= 0.0 { return false; }
    if s.nm < s.mrcount { return false; }
    let intron_len = j.acceptor.saturating_sub(j.donor);
    // default (narrowed) rule:
    if intron_len > longintron && s.nreads_good < chi_win_error { return true; }
    // RUSTLE_HIGHERR_SNAP: also snap NOISE donors (non-canonical, unguided, very low support)
    // so they can redirect to a nearby stronger canonical donor (matches ST's higherr snap),
    // while sparing canonical/guided real alt-donors (avoids the documented small-exon regression).
    if crate::stringtie_parity::st_higherr_snap()
        && !s.guide_match
        && s.consleft == 0
        && s.nreads_good < 10.0
    {
        return true;
    }
    false
};
```
(The candidate-side `is_reliable` already requires the stronger donor be `!is_bad` + `leftsupport>0.9*` + `ok_to_demote`; the canonical-strong check is implicit via the stronger donor passing `is_bad`=false. If 2A Phase-0 showed canonical cost remains, ALSO add `&& cand.consleft == 1` at the candidate strength check `:240`. Add only if needed.)

- [ ] **Step 2: Build + measure flag ON vs default**
```bash
cargo build --release 2>&1 | tail -2
RUSTLE_HIGHERR_SNAP=1 ./target/release/rustle -L ../GGO_19.bam -o /tmp/snap_on.gtf 2>/dev/null
gffcompare -r ../GGO_19.gtf /tmp/snap_on.gtf -o /tmp/snap_on_cmp 2>/dev/null
echo "=== HIGHERR_SNAP ON ==="; grep -E "Intron chain level:|Transcript level:" /tmp/snap_on_cmp.stats
python3 bench/gtf_chain_diff.py /tmp/snap_on.gtf ../GGO_19.gtf 2>/dev/null | tail -4
./target/release/rustle -L ../GGO_19.bam -o /tmp/snap_def.gtf 2>/dev/null
gffcompare -r ../GGO_19.gtf /tmp/snap_def.gtf -o /tmp/snap_def_cmp 2>/dev/null
echo "=== default OFF ==="; grep -E "Transcript level:" /tmp/snap_def_cmp.stats
```
Expected: flag-ON RU-only FP drops by the 2A prize, Pr up, Sn drop ≤ cost; flag-OFF MUST be 96.2/91.7, 95.6/90.5.

- [ ] **Step 3: Realizability cross-check + decision.** Confirm removed FPs are the near-dup-donor chains and Sn loss ≤ the (canonical/guided-spared) cost. Ship only on net F1 gain; else opt-in + record. Commit `src/rustle/junction_correction.rs` + findings.

---

## Self-review notes
- **Spec coverage:** 2B free probe → Task 1; 2B targeted prize → Task 2; 2B gate → Task 3; 2B fix → Tasks 4-5. 2A prize → Task 6; 2A gate → Task 7; 2A fix → Tasks 8-9. Realizability cross-check → Tasks 5/9 Step 3.
- **Default-unchanged guard:** Tasks 5 & 9 Step 2 re-verify 95.6/90.5 flag-OFF. All Phase-1 logic behind `st_csr_fold()` / `st_higherr_snap()`.
- **Feasibility confirmed (grounding):** 2B containment computable via `flow_kept_paths_len` + `out_idx` (no new state); 2A signals (`consleft`/`guide_match`/`nreads_good`) in scope. CSR tx findable via GTF source/rescue_class.
- **No placeholders:** complete code/commands per step; the two "add only if needed" clauses (canonical candidate check at `:240`; junction_accept emit enhancement) are explicitly conditional on the Phase-0 measurement, not deferred work.
- **Lesson enforced:** each sub-lever oracle-gated before code; TP-cost is the dominant trigger; realizability cross-check (removed FP = oracle set, read folds into surviving TP) in the ship decision.
