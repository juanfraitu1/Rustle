# Predcluster Selection-Parity Port (sub-project 1) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Under `RUSTLE_PREDCLUSTER_ST=1` (default OFF), make Rustle's per-cluster winner-selection reproduce the current StringTie binary's selection decisions, driving the 110 selection errors (81 Type-A FP + 29 Type-B FN) toward 0 on candidate-matching clusters; flip to default if it beats baseline F1 (93.78).

**Architecture:** A parallel ST-faithful selection function `select_predictions_st()` in a new `src/rustle/predcluster_st.rs`, dispatched in place of Rustle's predcluster when the flag is on. It runs StringTie's selection sub-stages in ST's order on the candidate predictions Rustle already produces: prediction-collapse parity → pairwise containment (RI unconditional + included_drop) → per-maxint isofrac. Oracle-first: bound the F1 ceiling before building.

**Tech Stack:** Rust (`cargo build --release` → `target/release/rustle`), instrumented StringTie fork `./tools/stringtie/stringtie`, Python 3 (`bench/gtf_chain_diff.py`, a new `bench/pred_kill_diff.py`), `gffcompare`, the parity-decisions JSONL (`pred_kill`, `path_extracted`, new `pred_merge`).

**Reference:** spec `docs/superpowers/specs/2026-05-30-predcluster-selection-parity-design.md`; map in memory `project_jfp_missr_characterization` (2026-05-29 PREDCLUSTER SELECTION MAP).

**Key facts the engineer needs:**
- Baseline (deterministic): TP=1750 / FN=64 / FP=168 / Sn=96.47 / Pr=91.24 / F1=93.78 vs `../GGO_19.gtf` (1814 multi-intron chains). Current ST scores 1754/60/67.
- chain_counts: multi-intron chains keyed (strand, intron-tuple), intron=(exon_end, next_start); TP=|ru∩ref|, Sn=TP/1814, Pr=TP/|ru|, F1=2TP/(2TP+FN+FP).
- Rustle predcluster entry: `print_predcluster_with_summary_multi` (transcript_filter.rs; called at transcript_filter.rs:7595). Sub-rules: isofrac `isofrac_with_summary` (~2161), RI `retainedintron_like` (~1666), pairwise (~3237-3263/3418).
- ST selection sites (de-novo -L): RI/containment `rlink.cpp:7363-7404` + `has_retained_intron` 7191; `included_drop` ~18438; isofrac/longunder per-maxint `rlink.cpp:18734-18794`. **De-novo prediction-collapse = the keeptrf containment loop `rlink.cpp:5588-5800` (NOT the --merge `merge_transfrags` at 7217)** — Rustle already ports this as `process_transfrags` (transfrag_process.rs:2050-2280).
- `Transcript` fields for the adapter (types.rs): `exons`, `coverage`, `longcov`, `hardstart`, `hardend`; intron-pattern derivable from `exons` + graph nodes. `compatible_long` (transfrag_process.rs:932) already implements ST's compatibility convention.
- Default `-L` de novo (no -G). The flag predicate `st_predcluster()` lives in stringtie_parity.rs (TDD'd in Task 2).
- Selection-divergent clusters & the 110 errors are defined in the map; the validation restricts to **candidate-matching clusters** (both tools' path_extracted sets equal at the locus).

---

### Task 1: Oracle — bound the selection ceiling (abort gate, NO production code)

**Files:** Create `bench/selection_oracle.py`

**Context:** Before building, prove the decomposition: if Rustle kept EXACTLY ST's winners on candidate-matching clusters, would the 110 selection errors close? This filters Rustle's candidate set to ST's final chains, restricted to clusters where the candidate sets already match, and measures the resulting chain-level F1.

- [ ] **Step 1: Capture both tools' candidates + finals**

Run:
```bash
RUSTLE_PARITY_LOG=/tmp/ru_pe.jsonl RUSTLE_PARITY_FILTER_CHROM=NC_073243.2 RUSTLE_PARITY_FILTER_STEPS=path_extracted ./target/release/rustle GGO_19.bam -L -o /tmp/ru.gtf 2>/dev/null
STRINGTIE_PARITY_LOG=/tmp/st_pe.jsonl STRINGTIE_PARITY_FILTER_CHROM=NC_073243.2 STRINGTIE_PARITY_FILTER_STEPS=path_extracted ./tools/stringtie/stringtie GGO_19.bam -L -o /tmp/st.gtf 2>/dev/null
```

- [ ] **Step 2: Write the oracle**

Create `bench/selection_oracle.py`:
```python
#!/usr/bin/env python3
"""Selection-ceiling oracle: if Rustle kept exactly ST's winners on candidate-matching
clusters, what is the chain-level F1? Bounds sub-project-1's payoff before building.
Inputs: /tmp/ru.gtf /tmp/st.gtf ../GGO_19.gtf, /tmp/ru_pe.jsonl /tmp/st_pe.jsonl."""
import json, re, collections

def gtf_chains(p):
    tx=collections.defaultdict(list); strand={}
    for line in open(p):
        if line.startswith('#'): continue
        f=line.rstrip('\n').split('\t')
        if len(f)<9 or f[2]!='exon': continue
        m=re.search(r'transcript_id "([^"]+)"',f[8])
        if not m: continue
        tx[m.group(1)].append((int(f[3]),int(f[4]))); strand[m.group(1)]=f[6]
    out={}
    for t,ex in tx.items():
        ex.sort()
        if len(ex)<2: continue
        out[(strand[t],tuple((ex[i][1],ex[i+1][0]) for i in range(len(ex)-1)))]=(ex[0][0],ex[-1][1])
    return out

def pe_chains(p):
    out=collections.defaultdict(set)  # never used as dict-of-set here; simple set
    s=set()
    for line in open(p):
        try: e=json.loads(line)
        except Exception: continue
        if e.get("step")!="path_extracted": continue
        istr=e["payload"].get("introns","")
        if not istr: continue
        # ST + Rustle path_extracted introns already == GTF (exon_end+1,next_start-1)? verify; map to GTF coords
        ch=tuple((int(d)-1,int(a)+1) for d,a in (tok.split("-") for tok in istr.split(",")))
        s.add((e["strand"],ch))
    return s

ref=gtf_chains('../GGO_19.gtf'); ru=gtf_chains('/tmp/ru.gtf'); st=gtf_chains('/tmp/st.gtf')
ru_pe=pe_chains('/tmp/ru_pe.jsonl'); st_pe=pe_chains('/tmp/st_pe.jsonl')

# cluster by overlapping span+strand over the union of all final+candidate chains
def cluster(chains_with_span):
    items=sorted((s,a,b,ch) for (s,ch),(a,b) in chains_with_span.items())
    cl={}; cur=None; lid=-1
    for s,a,b,ch in items:
        if cur is None or s!=cur[0] or a>cur[2]: lid+=1; cur=[s,a,b,lid]
        else: cur[2]=max(cur[2],b)
        cl[(s,ch)]=lid
    return cl
# build span map for all chains (finals)
allspan={**{(s,ch):v for (s,ch),v in ru.items()}, **{(s,ch):v for (s,ch),v in st.items()}, **{(s,ch):v for (s,ch),v in ref.items()}}
clid=cluster(allspan)

# candidate-matching clusters: ru_pe and st_pe identical within the cluster
pe_by_cl=collections.defaultdict(lambda:[set(),set()])
for k in ru_pe:
    if k in clid: pe_by_cl[clid[k]][0].add(k)
for k in st_pe:
    if k in clid: pe_by_cl[clid[k]][1].add(k)
matching_cl={c for c,(r,s) in pe_by_cl.items() if r==s and r}

# ORACLE output: on matching clusters keep ST's finals; elsewhere keep Rustle's finals
oracle=set()
for (s,ch),(a,b) in ru.items():
    c=clid.get((s,ch))
    if c in matching_cl:
        continue  # drop Rustle's winner; ST's will be added below
    oracle.add((s,ch))
for (s,ch),(a,b) in st.items():
    c=clid.get((s,ch))
    if c in matching_cl:
        oracle.add((s,ch))
refset=set(ref)
tp=len(oracle&refset); fn=len(refset-oracle); fp=len(oracle-refset)
f1=200*tp/(2*tp+fn+fp)
print(f"matching clusters: {len(matching_cl)}")
print(f"BASELINE Rustle: TP={len(set(ru)&refset)} FN={len(refset-set(ru))} FP={len(set(ru)-refset)}")
print(f"ORACLE (ST winners on matching clusters): TP={tp} FN={fn} FP={fp} F1={f1:.2f}")
```

- [ ] **Step 3: Run + evaluate the ABORT GATE**

Run: `python3 bench/selection_oracle.py`
- **PROCEED** if the oracle materially closes the selection errors — i.e. ORACLE FP drops well below baseline 168 (toward ~88 = 168−80 Type-A) and FN drops toward ~37 (64−27 Type-B), with F1 clearly above 93.78. This confirms the 110 errors are genuinely selection-fixable.
- **REPLAN/ABORT** if ORACLE barely moves FP/FN (selection errors are not purely selection — they're coverage/merge artifacts upstream). STOP and report; do not build the port.

- [ ] **Step 4: Commit**

```bash
git add bench/selection_oracle.py
git commit -m "feat(selection): oracle bounding the predcluster-parity ceiling

<paste ORACLE TP/FN/FP/F1 and PROCEED/ABORT verdict>

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: Scaffold `predcluster_st.rs` + `st_predcluster()` flag + pass-through dispatch

**Files:** Create `src/rustle/predcluster_st.rs`; Modify `src/rustle/stringtie_parity.rs`, `src/rustle/transcript_filter.rs`, `src/rustle/lib.rs` (module decl).

- [ ] **Step 1: Failing test for the flag predicate**

Add to `src/rustle/stringtie_parity.rs` tests mod:
```rust
#[test]
fn st_predcluster_default_off() {
    use super::st_predcluster_from;
    assert!(!st_predcluster_from(None));
    assert!(st_predcluster_from(Some("1")));
    assert!(!st_predcluster_from(Some("0")));
}
```

- [ ] **Step 2: Run → fail**

Run: `cargo test -p rustle st_predcluster_default_off 2>&1 | tail -10` → FAIL (function not found).

- [ ] **Step 3: Implement the predicate**

In `src/rustle/stringtie_parity.rs`:
```rust
#[inline]
pub fn st_predcluster_from(v: Option<&str>) -> bool { matches!(v, Some(s) if s != "0") }
/// ST-faithful predcluster selection (default OFF). Enable: RUSTLE_PREDCLUSTER_ST=1.
#[inline]
pub fn st_predcluster() -> bool {
    st_predcluster_from(std::env::var("RUSTLE_PREDCLUSTER_ST").ok().as_deref())
}
```

- [ ] **Step 4: Create the module skeleton (pass-through)**

Create `src/rustle/predcluster_st.rs`:
```rust
//! StringTie-faithful prediction SELECTION (sub-project 1). Default OFF
//! (RUSTLE_PREDCLUSTER_ST=1). Runs ST's selection sub-stages in ST's order on the
//! candidate predictions Rustle's flow already produced. See
//! docs/superpowers/specs/2026-05-30-predcluster-selection-parity-design.md.
use crate::types::Transcript;

/// Pass-through skeleton: returns all candidates kept (no selection yet).
/// Subsequent tasks add containment, isofrac, and collapse parity.
pub fn select_predictions_st(candidates: Vec<Transcript>) -> Vec<Transcript> {
    candidates
}
```
Add `pub mod predcluster_st;` to `src/rustle/lib.rs`.

- [ ] **Step 5: Dispatch (guarded, pass-through)**

At the predcluster call site (transcript_filter.rs:7595 region — inspect to find where the final kept set is produced/returned), add a shadow dispatch that, when `crate::stringtie_parity::st_predcluster()` is true, routes the candidate set through `crate::predcluster_st::select_predictions_st(...)` INSTEAD of the normal predcluster. With the pass-through skeleton this keeps all candidates (so the flag-on output ⊇ baseline) — that's expected at this scaffolding step; the point is the wiring, not parity yet.

- [ ] **Step 6: Build + verify default unchanged**

Run:
```bash
cargo build --release 2>&1 | grep -E "^error" | head; echo "exit=${PIPESTATUS[0]}"
./target/release/rustle GGO_19.bam -L -o /tmp/off.gtf 2>/dev/null
gffcompare -r ../GGO_19.gtf -R -Q -o /tmp/off /tmp/off.gtf 2>/dev/null; grep "Intron chain level" /tmp/off.stats
```
Expected: exit=0; default (flag off) Intron chain `96.5 | 91.2`-ish (unchanged; nondeterministic ±0.1). Flag-on will OVER-produce (no selection yet) — that's fine at this step.

- [ ] **Step 7: Commit**

```bash
git add src/rustle/predcluster_st.rs src/rustle/stringtie_parity.rs src/rustle/transcript_filter.rs src/rustle/lib.rs
git commit -m "feat(selection): scaffold predcluster_st + st_predcluster() flag + pass-through dispatch

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: Port `pairwise_containment_st` (RI unconditional + included_drop)

**Files:** Modify `src/rustle/predcluster_st.rs`; Create `bench/pred_kill_diff.py`.

**Context:** ST rlink.cpp:7363-7404: for prediction pair where i's introns ⊆ j's: if `has_retained_intron(j,i)` (7191: an intron edge in j not in i) AND `cov_i < cov_j` → kill i; ELSE (contained, no RI) → kill i UNCONDITIONALLY. Plus `included_drop` (18438): contained pred with cov < DROP·container → kill. Rustle's current RI (`retainedintron_like`) is bp-coverage-gated and SPARES what ST kills — the `_st` version is unconditional per ST.

- [ ] **Step 1: Build `pred_kill_diff.py`** (the gate) — joins Rustle vs ST `pred_kill` by (reason, intron-chain) restricted to candidate-matching clusters; reports, per reason, how many ST kills Rustle reproduces vs misses. (Mirror `bench/junction_accept_diff.py` structure; ST pred_kill carries no introns — match by coords+nexons as the map did, document the convention.)

- [ ] **Step 2: Implement the containment kill in `select_predictions_st`** — pairwise over candidates sorted by cov DESC; intron-pattern subset test (derive each candidate's intron set from `exons`); apply ST's exact branch (RI+lower-cov OR unconditional-contained); apply included_drop. Emit `pred_kill` reason=`retained_intron`/`included_drop` for each kill (gated on `crate::parity::decisions::is_enabled()`), matching the existing payload shape (introns field).

- [ ] **Step 3: Build + default-unchanged guard** (cargo build; flag-off Intron chain unchanged).

- [ ] **Step 4: Drive the gate**

Run rustle with `RUSTLE_PREDCLUSTER_ST=1` + pred_kill capture and ST pred_kill; `python3 bench/pred_kill_diff.py`. Expected: of ST's `retained_intron`+`included_drop` kills on candidate-matching clusters, Rustle now reproduces the large majority. Record before/after. Also measure flag-on chain-level TP/FP (it will improve on the Type-A RI/included FPs but isofrac FPs remain until Task 4).

- [ ] **Step 5: Commit** (`feat(selection): ST-faithful pairwise containment (RI unconditional + included_drop)`).

---

### Task 4: Port `isofrac_st` (per-maxint longunder)

**Files:** Modify `src/rustle/predcluster_st.rs`.

**Context:** ST rlink.cpp:18734-18794: iterate maximal-coverage intervals (maxint); per interval seed `usedcov[s]` from the single dominant prediction through it; kill `cov < isofraclong·usedcov[s]` (isofraclong relaxed for high-cov via CHI_WIN); accumulate `usedcov[s] += cov` per interval. This is the per-interval reference (NOT Rustle's accumulated cross-interval usedcov).

- [ ] **Step 1: Implement `isofrac_st`** in `select_predictions_st` after the containment stage: build maxint intervals from the kept candidates' exon structure; per interval compute the dominant; kill candidates below `isofraclong·usedcov`; emit `pred_kill` reason=`isofrac`. Use ST's constants (ERROR_PERC=0.1, DROP, CHI_WIN, isofraclong relax) — grep ST for values.

- [ ] **Step 2: Build + default guard.**

- [ ] **Step 3: Drive the gate** — `pred_kill_diff.py`: ST's `isofrac` kills on candidate-matching clusters reproduced by Rustle. Record. Flag-on chain-level TP/FP.

- [ ] **Step 4: Commit** (`feat(selection): ST-faithful per-maxint isofrac`).

---

### Task 5: Prediction-collapse parity (the de-novo keeptrf collapse, NOT --merge)

**Files:** Modify `src/rustle/transfrag_process.rs` (process_transfrags keeptrf), possibly `src/rustle/predcluster_st.rs`; add a `pred_merge` parity event.

**Context:** The 25 "merge-collapse" FP + 27 Type-B FN come from the de-novo containment-collapse. ST does this in the keeptrf loop `rlink.cpp:5588-5800`; Rustle already ports it as `process_transfrags` (transfrag_process.rs:2050-2280). So this task is a DIFF-AND-ALIGN, not a fresh port: instrument both with a `pred_merge`/`transfrag_collapse` event (the shipped `transfrag_collapse` from b0d13ed may already suffice — reuse it), diff on the selection-divergent clusters, and identify whether the 25/27 are (a) a divergence in the existing keeptrf collapse (align it under `st_predcluster()`/`st_shadow()`), or (b) actually candidate-EXTRACTION artifacts → if (b), DROP them from this sub-project's scope and document for sub-project 2.

- [ ] **Step 1: Capture `transfrag_collapse` both tools** (already wired, b0d13ed) on the selection-divergent clusters; quantify how many of the 25 FP / 27 FN are explained by a keeptrf group-membership divergence vs absent upstream.
- [ ] **Step 2:** If keeptrf-collapse divergence: align the divergent rule (transfrag_process.rs) gated on the flag; drive the `transfrag_collapse` diff on these clusters → match. If extraction-artifact: document + defer (no code), commit the finding.
- [ ] **Step 3: Build + default guard; commit** (`feat/docs(selection): keeptrf-collapse parity for de-novo merge-collapse cases`).

---

### Task 6: Validate combined selection parity + decide default

**Files:** Modify `docs/STRINGTIE_PARITY_FINDINGS.md`.

- [ ] **Step 1: Combined gate** — with `RUSTLE_PREDCLUSTER_ST=1`, run `bench/gtf_chain_diff.py(Rustle, ST)` restricted to selection-divergent clusters → the 110 selection errors → 0 (or accounted). Run `bench/pred_kill_diff.py` → all selection reasons match on candidate-matching clusters.
- [ ] **Step 2: Full F1** — `RUSTLE_PREDCLUSTER_ST=1` chain-level TP/FN/FP/Sn/Pr/F1 vs reference, 3 runs (determinism). Compare to baseline 1750/64/168 (F1 93.78). Expected: Pr up (~+3pp from the 80 Type-A FP), Sn ~flat-to-up (Type-B FN recovered minus any new losses).
- [ ] **Step 3: Default decision** — if flag-on F1 strictly > 93.78 with Pr up and Sn not hurt, flip `st_predcluster()` default-on (opt-out `RUSTLE_PREDCLUSTER_ST=0`); else keep opt-in. Verify default-OFF regression unchanged.
- [ ] **Step 4: Record** — append "Predcluster selection-parity (sub-project 1) — DONE" to `docs/STRINGTIE_PARITY_FINDINGS.md`: selection-error before/after, flag-on F1, default decision, residual = candidate-extraction (sub-project 2). Commit.

---

## Self-Review

- **Spec coverage:** §3 architecture (parallel `_st` fn, dispatch, flag) → Task 2. §4 components: containment → Task 3, isofrac → Task 4, collapse → Task 5, orchestrator → Tasks 2-4. §6 validation (pred_kill/pred_merge/gtf_chain_diff gates, candidate-matching restriction) → Tasks 3,4,6. §7 oracle-first + abort → Task 1; default-flip → Task 6. §2 non-goals (candidate-extraction, guided) honored — Task 5 explicitly defers extraction-artifacts.
- **Placeholder note:** Tasks 3-5 specify behavior by ST source-of-truth + the parity gate + the discipline rather than literal final Rust, because the exact port is read-from-ST (a research port); the oracle/scaffold/gates/commands are concrete. The `merge_transfrags`-is-merge-mode correction is captured (Task 5 targets the keeptrf collapse).
- **Type consistency:** `st_predcluster()`/`st_predcluster_from()` and `select_predictions_st(Vec<Transcript>) -> Vec<Transcript>` used consistently; `pred_kill` reasons (`retained_intron`/`included_drop`/`isofrac`) match ST's and the existing payload shape; gates are the established parity-diff pattern.
- **Risk:** Task 1's abort gate can end the plan early (by design — cheap kill). Task 5 may resolve to "defer to sub-project 2" (also fine).
