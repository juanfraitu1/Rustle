# Retained-Intron Scope Widening — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Measure — oracle-first, default-OFF — whether widening Rustle's retained-intron check to StringTie's pairwise/local scope removes over-enumeration FPs genome-wide at acceptable Sn cost; diagnose which RU scope under-fires; decide ship/shelve from the oracle.

**Architecture:** Phase 0a bounds the total prize (ST `pred_kill(retained_intron)` attributed onto RU final chains, net of TP-cost). Phase 0b diagnoses flags-vs-pairing via the existing `RUSTLE_LOWINTRON_ORACLE` fed by ST's `pred_intron_low`. A hard abort gate follows. Phase 1 widens the identified scope behind `RUSTLE_RI_LOCAL` only if the gate clears. Both phases require ST to emit the intron CHAIN in `pred_kill` + `pred_intron_low` (current payloads lack it, making cross-tool matching ambiguous/impossible).

**Tech Stack:** C++ (instrumented StringTie, `tools/stringtie/`), Rust (Rustle filter), Python 3 stdlib (bench), gffcompare.

---

## Grounding facts (verified — confirm names at edit time)

- **ST `pred_kill(retained_intron)` emit:** `tools/stringtie/rlink.cpp:18599-18612`; payload `"reason":"retained_intron","cov":%.4f,"nexons":%d,"stage":"pairwise"` — NO chain, NO killer. The killed prediction is `pred[n2]` (its `->exons` are in scope).
- **ST `pred_intron_low` emit:** `tools/stringtie/rlink.cpp:18275-18277` (pred[0]) + `18399-18401` (pred[n]); payload has `intron_low` mask (`"1,1,0"`), `intron_covs`, `exon_covs`, `cov`, `longcov`, `tlen`, `pileup_cov` — but NO chain coords. Exons in scope as `pred->exons`.
- **RU `pred_kill` emit (reference format):** `transcript_filter.rs:3166-3187`, payload includes `reason,cov,nexons,stage,killer_cov,killer_nexons,killer_start,killer_end` (RU already emits killer info; ST should mirror it).
- **`RUSTLE_LOWINTRON_ORACLE` loader:** `transcript_filter.rs:1533-1555`. Reads a file path; each line = `<strand> <mask> <chain>` whitespace-separated. **`mask` is parsed PER-CHARACTER** (`mask.chars().map(|c| c=='1')`, line 1549) → the mask MUST be comma-free (`"110"`, not `"1,1,0"`). `chain` = the key (trailing comma trimmed).
- **RU oracle key format (the chain):** built at `transcript_filter.rs:1648-1653` as `"{exon[j-1].1 + 1}-{exon[j].0}"` per interior intron, comma-joined, NO trailing comma. RU exons are **0-based half-open** `[start,end)`. So for a junction with 1-based donor D / acceptor A: RU token = `"{D+1}-{A-1}"` (because RU `exon[j-1].1 == D`, RU `exon[j].0 == A-1`). **ST native is 1-based inclusive**, so ST must emit `"{exon[z-1].end+1}-{exon[z].start-1}"` to match (note the `-1` on the acceptor). THIS MUST BE VERIFIED byte-for-byte (Task 1 Step 5) — an off-by-one here silently no-ops the oracle.
- **Fix site FLAGS:** `build_lowintron_flags` `transcript_filter.rs:1557-1644` (decides low-cov via the tx's OWN bpcov drop, `:1620-1643`); `ord_neighbors`/`txs` available to widen toward a local-exonic-killer comparison.
- **Fix site PAIRING:** pairwise loop `transcript_filter.rs:3063-3071`; `ord` sorted by `tx_score` desc (`exonic_len*cov`, `:3016-3029`); `overlaps_adj` from `build_significant_overlap_*` (`:3032-3052`); skip `ord_pos[n2] <= oi` (`:3069`).
- **Flag predicate pattern:** `stringtie_parity.rs:213` `matches!(v, Some(s) if s != "0")`.
- **Genome-wide capture:** CHROM filter, no RANGE, all steps → ~4-6 MB log; no size cap. `bench/gtf_chain_diff.py <ru.gtf> <st.gtf>`; gffcompare `-r REF -o cmp out.gtf`; baseline 96.2/91.7 IC, 95.6/90.5 Tx.

---

## File structure
- Modify `tools/stringtie/rlink.cpp` — add intron chain (+ killer info) to `pred_kill(retained_intron)` + intron chain to `pred_intron_low`.
- Create `bench/retained_intron_prize.py` — Phase 0a prize-bound.
- Create `bench/build_lowintron_oracle.py` — Phase 0b oracle-file builder (comma-free mask + chain).
- Modify `src/rustle/stringtie_parity.rs` — `RUSTLE_RI_LOCAL` predicate (Phase 1).
- Modify `src/rustle/transcript_filter.rs` — the widened scope (Phase 1, flags or pairing).
- Modify `docs/STRINGTIE_PARITY_FINDINGS.md` + memory — record.

---

## PHASE 0 — Oracles & abort gate

### Task 1: ST emits the intron chain in `pred_kill` + `pred_intron_low`

**Files:** Modify `tools/stringtie/rlink.cpp` (submodule).

- [ ] **Step 1: Add a chain helper near the pred_kill emit**

In `rlink.cpp`, before the `pred_kill(retained_intron)` emit (~18605), build the killed prediction's chain string from `pred[n2]->exons` in RU's key convention (`{end+1}-{start-1}` per interior intron, comma-joined):

```cpp
// Intron chain in Rustle's RUSTLE_LOWINTRON_ORACLE key convention:
//   token = "{exon[z-1].end+1}-{exon[z].start-1}"  (RU uses 0-based acceptor)
GStr __ri_chain;
for (int z = 1; z < pred[n2]->exons.Count(); z++) {
  if (z > 1) __ri_chain.append(",");
  __ri_chain.appendfmt("%u-%u", pred[n2]->exons[z-1].end + 1, pred[n2]->exons[z].start - 1);
}
```

- [ ] **Step 2: Extend the `pred_kill(retained_intron)` payload**

Change the emit (`rlink.cpp:18607-18609`) to include `chain` + killer info (mirror RU). `pred[n1]` is the killer; use its exon count + span + cov:

```cpp
char __pk[512];
snprintf(__pk, sizeof(__pk),
  "\"reason\":\"retained_intron\",\"cov\":%.4f,\"nexons\":%d,\"stage\":\"pairwise\","
  "\"killer_cov\":%.4f,\"killer_nexons\":%d,\"killer_start\":%u,\"killer_end\":%u,\"chain\":\"%s\"",
  pred[n2]->cov, pred[n2]->exons.Count(),
  pred[n1]->cov, pred[n1]->exons.Count(),
  pred[n1]->exons[0].start, pred[n1]->exons.Last().end,
  __ri_chain.chars());
pd_emit("pred_kill", cname, (unsigned long)pred[n2]->exons[0].start, (unsigned long)pred[n2]->exons.Last().end, sign, __pk);
```
(Confirm the existing emit's variable names — `cname`, `sign`, the pred index `n1`/`n2`, the snprintf buffer — and reuse them; only ADD the killer + chain fields.)

- [ ] **Step 3: Add the chain to `pred_intron_low`**

At both `pred_intron_low` emits (`~18275-18277` and `~18399-18401`), build the same chain string from `pred->exons` (the loop variable; `pred[0]` for the first site, `pred[n]` for the second) and append `,"chain":"<chain>"` to the payload. Use the identical `{end+1}-{start-1}` convention.

- [ ] **Step 4: Build StringTie**

Run: `( cd tools/stringtie && make clean release 2>&1 | tail -3 )`
Expected: clean rebuild (~10s).

- [ ] **Step 5: VERIFY the chain key byte-matches RU's (CRITICAL — guards the oracle)**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
STRINGTIE_PARITY_LOG=/tmp/st_chain.jsonl STRINGTIE_PARITY_FILTER_CHROM=NC_073243.2 ./tools/stringtie/stringtie -L GGO_19.bam -o /tmp/st_chk.gtf 2>/dev/null
# pull one pred_intron_low chain from ST:
grep -m1 '"step":"pred_intron_low"' /tmp/st_chain.jsonl | python3 -c 'import sys,json; print(json.loads(sys.stdin.readline())["payload"]["chain"])'
```
Then capture RU's key for the SAME transcript: run RU with `RUSTLE_PARITY_FILTER_STEPS=pred_intron_low` (RU emits it at pipeline.rs:7346) on the same locus and compare, OR add a one-off `eprintln!` of RU's key at `transcript_filter.rs:1653` for a tx whose span matches the ST one, and confirm the strings are **byte-identical**. If they differ (off-by-one on donor or acceptor), adjust the ST `+1`/`-1` until identical. DO NOT proceed until a shared transcript's ST chain == RU key exactly.

- [ ] **Step 6: Commit (submodule)**

```bash
git -C tools/stringtie add rlink.cpp
git -C tools/stringtie commit -m "parity(stringtie): emit intron chain (+ killer) in pred_kill + pred_intron_low

RU-key convention {end+1}-{start-1}; enables exact cross-tool chain matching for
the retained-intron oracle. Env-gated; no behavior change.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: Genome-wide capture (both tools, all steps)

**Files:** none (uses existing binaries).

- [ ] **Step 1: Capture chr19 genome-wide parity + finals**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
RUSTLE_PARITY_LOG=/tmp/ri_ru.jsonl RUSTLE_PARITY_FILTER_CHROM=NC_073243.2 \
  ./target/release/rustle -L GGO_19.bam -o /tmp/ri_ru.gtf 2>/dev/null
STRINGTIE_PARITY_LOG=/tmp/ri_st.jsonl STRINGTIE_PARITY_FILTER_CHROM=NC_073243.2 \
  ./tools/stringtie/stringtie -L GGO_19.bam -o /tmp/ri_st.gtf 2>/dev/null
wc -l /tmp/ri_ru.jsonl /tmp/ri_st.jsonl
grep -c '"reason":"retained_intron"' /tmp/ri_st.jsonl
grep -c '"step":"pred_intron_low"' /tmp/ri_st.jsonl
```
Expected: logs a few-MB each; ST `retained_intron` pred_kills > 0 (the trace saw them); `pred_intron_low` in the hundreds. Confirm `chain` is present: `grep -m1 '"reason":"retained_intron"' /tmp/ri_st.jsonl | grep -o '"chain":"[^"]*"'`.

- [ ] **Step 2: Baseline RU final vs reference (headline)**

Run: `gffcompare -r ../GGO_19.gtf /tmp/ri_ru.gtf -o /tmp/ri_ru_cmp 2>/dev/null; grep -E "Intron chain level:|Transcript level:" /tmp/ri_ru_cmp.stats`
Expected: 96.2/91.7, 95.6/90.5. (No commit — capture step.)

---

### Task 3: `bench/retained_intron_prize.py` — Phase 0a total-prize bound

**Files:** Create `bench/retained_intron_prize.py`

- [ ] **Step 1: Write the prize tool**

```python
#!/usr/bin/env python3
"""Phase 0a: bound the retained-intron prize. Match ST pred_kill(retained_intron) chains
to RU final chains; net = (RU-only-FP matched) - (RU-TP matched).
Usage: python3 bench/retained_intron_prize.py /tmp/ri_st.jsonl /tmp/ri_ru.gtf /mnt/c/Users/jfris/Desktop/GGO_19.gtf"""
import json, sys, re
from collections import defaultdict

def gtf_chains(path):
    tx = defaultdict(list); strand = {}
    for line in open(path):
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon": continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m: continue
        tx[m.group(1)].append((int(f[3]), int(f[4]))); strand[m.group(1)] = f[6]
    chains = set()
    for t, ex in tx.items():
        ex.sort()
        # chain token in ST/RU key convention: {exon[j-1].end+1}-{exon[j].start-1}
        ic = tuple((ex[i-1][1] + 1, ex[i][0] - 1) for i in range(1, len(ex)))
        if ic: chains.add((strand[t], ic))
    return chains

def st_killed_chains(path):
    out = set()
    for line in open(path):
        if '"step":"pred_kill"' not in line: continue
        try: r = json.loads(line)
        except ValueError: continue
        p = r.get("payload", {})
        if p.get("reason") != "retained_intron": continue
        ch = p.get("chain", "")
        if not ch: continue
        ic = tuple(tuple(map(int, tok.split("-"))) for tok in ch.split(","))
        out.add((r["strand"], ic))
    return out

def main():
    st_killed = st_killed_chains(sys.argv[1] if len(sys.argv) > 1 else "/tmp/ri_st.jsonl")
    ru = gtf_chains(sys.argv[2] if len(sys.argv) > 2 else "/tmp/ri_ru.gtf")
    ref = gtf_chains(sys.argv[3] if len(sys.argv) > 3 else "/mnt/c/Users/jfris/Desktop/GGO_19.gtf")
    print(f"ST retained_intron-killed chains: {len(st_killed)}")
    print(f"RU final chains: {len(ru)}  (RU-only-FP vs ref: {len(ru - ref)})")
    matched = ru & st_killed                       # RU chains ST would kill via retained_intron
    fp_removed = matched - ref                      # RU-only-FP among them (the prize)
    tp_lost = matched & ref                          # RU-TP among them (the cost — ST over-kill)
    print(f"RU chains matched by an ST retained_intron kill: {len(matched)}")
    print(f"  of which FP (prize): {len(fp_removed)}")
    print(f"  of which TP (cost):  {len(tp_lost)}")
    net = len(fp_removed) - len(tp_lost)
    print(f"NET prize (FP_removed - TP_lost): {net}")
    print(f"ABORT GATE: {'ABORT (shelve Phase 1)' if net <= 0 else ('WEAK (<5)' if net < 5 else 'PROCEED')}")
    for c in list(tp_lost)[:8]:
        print(f"    TP-COST chain {c[0]} {c[1][:3]}...")

if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run it**

Run: `python3 bench/retained_intron_prize.py /tmp/ri_st.jsonl /tmp/ri_ru.gtf ../GGO_19.gtf`
Expected: prints matched/FP/TP/net + gate verdict. SANITY: if "ST retained_intron-killed chains" is 0, Task 1's chain emit didn't land — debug before trusting. If matched=0 but ST killed >0 and RU-only-FP >0, suspect the chain-key convention mismatch (Task 1 Step 5) — re-verify byte-match.

- [ ] **Step 3: Commit**

```bash
git add bench/retained_intron_prize.py
git commit -m "bench: retained_intron_prize.py (Phase 0a prize-bound via ST chain match)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: `bench/build_lowintron_oracle.py` + run Phase 0b

**Files:** Create `bench/build_lowintron_oracle.py`

- [ ] **Step 1: Write the oracle-file builder**

```python
#!/usr/bin/env python3
"""Phase 0b: build the RUSTLE_LOWINTRON_ORACLE file from ST pred_intron_low events.
File format (per line, whitespace-separated): "<strand> <comma-FREE mask> <chain>".
The loader (transcript_filter.rs:1549) parses the mask PER-CHARACTER, so commas must be stripped.
Usage: python3 bench/build_lowintron_oracle.py /tmp/ri_st.jsonl /tmp/lowintron.oracle"""
import json, sys

def main():
    src = sys.argv[1] if len(sys.argv) > 1 else "/tmp/ri_st.jsonl"
    out = sys.argv[2] if len(sys.argv) > 2 else "/tmp/lowintron.oracle"
    n = 0
    with open(out, "w") as fh:
        for line in open(src):
            if '"step":"pred_intron_low"' not in line: continue
            try: r = json.loads(line)
            except ValueError: continue
            p = r.get("payload", {})
            chain = p.get("chain", "")
            mask = p.get("intron_low", "")
            if not chain or not mask: continue
            mask_nc = mask.replace(",", "")          # comma-FREE: loader parses per-char
            if len(mask_nc) != chain.count(",") + 1:  # mask length must == #introns
                continue
            fh.write(f"{r['strand']} {mask_nc} {chain}\n")
            n += 1
    print(f"wrote {n} oracle entries to {out}")

if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Build the oracle file + run RU with it**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
python3 bench/build_lowintron_oracle.py /tmp/ri_st.jsonl /tmp/lowintron.oracle
RUSTLE_LOWINTRON_ORACLE=/tmp/lowintron.oracle ./target/release/rustle -L GGO_19.bam -o /tmp/ri_ru_orc.gtf 2>/dev/null
gffcompare -r ../GGO_19.gtf /tmp/ri_ru_orc.gtf -o /tmp/ri_orc_cmp 2>/dev/null
echo "=== oracle ON ==="; grep -E "Intron chain level:|Transcript level:" /tmp/ri_orc_cmp.stats
echo "=== baseline ===";  grep -E "Intron chain level:|Transcript level:" /tmp/ri_ru_cmp.stats
python3 bench/gtf_chain_diff.py /tmp/ri_ru_orc.gtf /tmp/ri_st.gtf | tail -3
```

- [ ] **Step 3: VERIFY the oracle actually fired (guards against silent no-op)**

The oracle ON output MUST differ from baseline (fewer tx and/or higher Pr). If the GTFs are byte-identical, the override silently no-opped — almost certainly the mask/chain key mismatch. Debug: `head -3 /tmp/lowintron.oracle` (mask must be comma-free like `110`); and confirm a RU tx whose chain is a key in the file actually changed. Do NOT trust the 0b number until a difference is confirmed.

- [ ] **Step 4: Record 0b numbers + diagnose**

Compute FP-reduction (baseline RU-only-FP minus oracle RU-only-FP, via `gtf_chain_diff` vs ref), TP-cost (TPs lost), Sn delta. DIAGNOSIS: if 0b FP-reduction ≈ 0a prize → the under-firing is the lowintron FLAGS (Phase-1 fix = `build_lowintron_flags`); if 0b ≪ 0a → the killer-PAIRING scope (Phase-1 fix = the pairwise loop). Record both.

- [ ] **Step 5: Commit**

```bash
git add bench/build_lowintron_oracle.py
git commit -m "bench: build_lowintron_oracle.py (Phase 0b — ST pred_intron_low -> RUSTLE_LOWINTRON_ORACLE)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 5: Abort-gate decision + record

**Files:** Modify `docs/STRINGTIE_PARITY_FINDINGS.md`; memory `project_needy_loci_decision_trace.md`.

- [ ] **Step 1: Read the gate**

Net prize from 0a (FP_removed − TP_lost). **ABORT if ≤ 0 or < ~5.** Cross-check with 0b: if 0a and 0b prize estimates disagree wildly, investigate (chain-key mismatch / oracle no-op) before trusting either. The dominant abort signal is TP-cost (ST over-killing real alt-isoforms RU correctly keeps).

- [ ] **Step 2: Record**

Append to findings §6j area + update memory: 0a net prize + TP-cost, 0b FP-reduction/Sn, the flags-vs-pairing diagnosis. If ABORT: record "retained-intron prize bounded as <gate>; trace's local finding does/doesn't generalize"; STOP — do not start Phase 1.

- [ ] **Step 3: Commit**

```bash
git add docs/STRINGTIE_PARITY_FINDINGS.md
git commit -m "docs: Phase-0 retained-intron oracle result (prize + flags/pairing diagnosis + gate)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## PHASE 1 — Targeted scope fix (ONLY if abort gate cleared)

### Task 6: `RUSTLE_RI_LOCAL` flag predicate

**Files:** Modify `src/rustle/stringtie_parity.rs`

- [ ] **Step 1: Write the failing test** (in the `#[cfg(test)]` module):
```rust
#[test]
fn st_ri_local_default_off() {
    use super::st_ri_local_from;
    assert!(!st_ri_local_from(None));
    assert!(st_ri_local_from(Some("1")));
    assert!(!st_ri_local_from(Some("0")));
}
```
- [ ] **Step 2: Run → FAIL** `cargo test -p rustle st_ri_local_default_off 2>&1 | tail -3`
- [ ] **Step 3: Implement** (next to `st_flow`, ~line 213):
```rust
pub fn st_ri_local_from(v: Option<&str>) -> bool { matches!(v, Some(s) if s != "0") }
pub fn st_ri_local() -> bool { st_ri_local_from(std::env::var("RUSTLE_RI_LOCAL").ok().as_deref()) }
```
- [ ] **Step 4: Run → PASS**; `cargo build --release 2>&1 | tail -2` clean.
- [ ] **Step 5: Commit**
```bash
git add src/rustle/stringtie_parity.rs
git commit -m "feat(parity): RUSTLE_RI_LOCAL predicate (default OFF)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

### Task 7: Widen the identified scope behind `RUSTLE_RI_LOCAL`

**Files:** Modify `src/rustle/transcript_filter.rs` (the site 0b identified).

- [ ] **Step 1: Implement the widening**

**If 0b diagnosed FLAGS:** in `build_lowintron_flags` (`transcript_filter.rs:1620-1643`), under `crate::stringtie_parity::st_ri_local()`, additionally flag a victim intron low when a LOCAL overlapping prediction is exonic over it with much higher coverage. Use `ord_neighbors`/`txs` (in scope per grounding) to find an overlapping prediction whose exon spans the intron interval, and compare `intron_cov < max_local_exon_cov * intronfrac`. Show the exact added block (derive the local-exon-cov from `bpcov` over the overlapping prediction's exon span).

**If 0b diagnosed PAIRING:** in the pairwise loop (`transcript_filter.rs:3063-3071`), under `st_ri_local()`, additionally consider local high-coverage overlapping killers excluded by the `ord_pos[n2] <= oi` skip — i.e. for the retained-intron check specifically, also test victim `n2` against higher-coverage overlapping predictions regardless of `ord` position. Show the exact added pairing (guard it to the retained-intron test only, so it doesn't perturb other pairwise reasons).

(The exact block is conditional on the Task-5 diagnosis; write only the branch 0b selected, with complete code derived from the surrounding logic — do not add both speculatively.)

- [ ] **Step 2: Build + measure flag ON vs default**

```bash
cargo build --release 2>&1 | tail -2
RUSTLE_RI_LOCAL=1 ./target/release/rustle -L ../GGO_19.bam -o /tmp/ri_on.gtf 2>/dev/null
gffcompare -r ../GGO_19.gtf /tmp/ri_on.gtf -o /tmp/ri_on_cmp 2>/dev/null
echo "=== RI_LOCAL ON ==="; grep -E "Intron chain level:|Transcript level:" /tmp/ri_on_cmp.stats
./target/release/rustle -L ../GGO_19.bam -o /tmp/ri_def.gtf 2>/dev/null
gffcompare -r ../GGO_19.gtf /tmp/ri_def.gtf -o /tmp/ri_def_cmp 2>/dev/null
echo "=== default (flag OFF) ==="; grep -E "Transcript level:" /tmp/ri_def_cmp.stats
```
Expected: flag-ON should reduce j-FPs (Pr up) toward the 0a prize; flag-OFF MUST stay 95.6/90.5.

- [ ] **Step 3: Decision + commit**

If flag-ON net F1 ≥ baseline (Pr up, Sn not materially hurt): recommend default-flip to the user (do NOT flip without explicit approval). Else keep `RUSTLE_RI_LOCAL` opt-in; record the cost vs the oracle ceiling. Commit:
```bash
git add src/rustle/transcript_filter.rs docs/STRINGTIE_PARITY_FINDINGS.md
git commit -m "feat(parity): RUSTLE_RI_LOCAL widens retained-intron scope to local killers (default OFF)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Self-review notes

- **Spec coverage:** Phase 0a prize → Tasks 1-3; Phase 0b diagnosis → Tasks 1,4; abort gate → Task 5; Phase 1 fix → Tasks 6-7. The spec's "enhance ST emit if ambiguous" → Task 1 (grounding confirmed it IS ambiguous, so it's required, not optional).
- **Two silent-failure guards baked in:** Task 1 Step 5 (chain key byte-match — guards both phases) and Task 4 Step 3 (oracle-fired check — guards 0b). These are the segmentation-prize lesson applied: validate the oracle moved something before trusting its number.
- **Default-unchanged guard:** Task 7 Step 2 re-verifies 95.6/90.5 flag-OFF. Phase 0 changes are ST-side (parity-gated) + Python + the existing default-OFF `RUSTLE_LOWINTRON_ORACLE` — zero RU default change until Phase 1, which is behind `RUSTLE_RI_LOCAL`.
- **Type/name consistency:** chain convention `{end+1}-{start-1}` used identically in Task 1 (ST emit), Task 3 (RU gtf_chains), and matches the RU key builder; `st_ri_local`/`st_ri_local_from` predicates; comma-free mask in Task 4 matches the loader's per-char parse.
- **No placeholders:** every code/command step is complete; Task 7's branch is explicitly conditional on the Task-5 diagnosis (write only the selected branch), not a deferred TODO.
