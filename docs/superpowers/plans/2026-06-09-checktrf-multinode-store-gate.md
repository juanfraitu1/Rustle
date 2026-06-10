# Checktrf Multi-Node Store Gate — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

> **OUTCOME (2026-06-09):** Tasks 1–4 executed. Task 4 validation FALSIFIED the default-on premise
> (gate over-rejects 156 ST-shared real isoforms; ST-only 104→259; −18 RefSeq TP). Narrowing
> investigation: no gate-time feature separates the 43 FPs from the 158 real recoveries (same wall
> as Approach A). Disposition: gate flipped to **opt-in, default-off** (`RUSTLE_CHECKTRF_MULTINODE_DROP=1`),
> default byte-identical. Tasks 5 (genome-wide) and 6 (this is partly done via §6o) are moot for a
> default change. See `docs/STRINGTIE_PARITY_FINDINGS.md` §6o. Next lever: flow-enumeration parity.

**Goal:** Make rustle's `parse_trflong` checktrf path StringTie-faithful — a multi-node long-read transfrag with no kept-path match is dropped (not stored as an independent prediction), eliminating ~48 false-positive intron chains on chr19.

**Architecture:** A single structural gate at the entry of the checktrf independent-rescue block in `src/rustle/path_extract.rs`. A pure predicate function (unit-tested) decides the drop; the call site reads the opt-out env, records a seed outcome, emits a parity trace event, and `continue`s. Mirrors StringTie's `if(!shortread && nodes.Count()>1)` redistribute-only branch (rlink.cpp:10369) whose independent-store `else` (rlink.cpp:10413) is unreachable for multi-node long-read transfrags. Default-on, opt-out via `RUSTLE_CHECKTRF_MULTINODE_RESCUE=1`.

**Tech Stack:** Rust (cargo), Python 3 for validation harness, gffcompare, the instrumented StringTie binary at `tools/stringtie/stringtie`.

**Spec:** `docs/superpowers/specs/2026-06-09-checktrf-multinode-store-gate-design.md`

**Working directory for all commands:** `/mnt/c/Users/jfris/Desktop/Rustle`

---

### Task 1: Pure drop-predicate function + unit tests

The decision logic is extracted into a pure, module-level function so it can be unit-tested in isolation (the call site is buried inside a ~4000-line function). The function takes plain bools/usize — no graph or transfrag types — so tests need no fixtures.

**Files:**
- Modify: `src/rustle/path_extract.rs` (add function after the `SeedOutcome` enum, which ends at line 599; add a test module near the existing `#[cfg(test)] mod capacity_confidence_tests` at line 12074)

- [ ] **Step 1: Write the failing test**

Add a new test module at the end of `src/rustle/path_extract.rs` (after the closing `}` of `mod capacity_confidence_tests`):

```rust
#[cfg(test)]
mod checktrf_gate_tests {
    use super::checktrf_multinode_no_match_drop;

    // multi-node long-read, non-guide, non-csr, opt-out unset -> DROP (the fix)
    #[test]
    fn drops_multinode_longread_no_match() {
        assert!(checktrf_multinode_no_match_drop(false, 2, false, false, false));
        assert!(checktrf_multinode_no_match_drop(false, 9, false, false, false));
    }

    // short-read transfrag -> KEEP (StringTie's else branch stores it)
    #[test]
    fn keeps_shortread() {
        assert!(!checktrf_multinode_no_match_drop(true, 5, false, false, false));
    }

    // single inner node -> KEEP (nodes.Count() <= 1)
    #[test]
    fn keeps_single_node() {
        assert!(!checktrf_multinode_no_match_drop(false, 1, false, false, false));
        assert!(!checktrf_multinode_no_match_drop(false, 0, false, false, false));
    }

    // guide -> KEEP (StringTie always recovers guides)
    #[test]
    fn keeps_guide() {
        assert!(!checktrf_multinode_no_match_drop(false, 4, false, true, false));
    }

    // csr_triggered chimeric-suffix fold -> KEEP (legitimate RUSTLE_CSR_FOLD case)
    #[test]
    fn keeps_csr_triggered() {
        assert!(!checktrf_multinode_no_match_drop(false, 4, true, false, false));
    }

    // opt-out env set -> KEEP (restore old store-it behavior)
    #[test]
    fn keeps_when_opt_out() {
        assert!(!checktrf_multinode_no_match_drop(false, 4, false, false, true));
    }
}
```

- [ ] **Step 2: Run the test to verify it fails (function not defined)**

Run: `cargo test --release checktrf_gate_tests 2>&1 | tail -20`
Expected: compile error `cannot find function 'checktrf_multinode_no_match_drop' in this scope` (or the import fails).

- [ ] **Step 3: Implement the predicate function**

Insert immediately after the `SeedOutcome` enum closing brace (after line 599, `}` of `pub enum SeedOutcome`):

```rust
/// StringTie-faithful checktrf store gate (rlink.cpp:10369/10413).
///
/// In StringTie, a multi-node long-read transfrag with no kept-path match is handled by the
/// `if(!shortread && nodes.Count()>1)` redistribute-only branch; the independent-store `else`
/// branch is structurally unreachable for such transfrags, so they are NEVER stored as a new
/// prediction. Rustle historically fell through and stored them (the checktrf rustle-only FP
/// chains, ~48 on chr19). This predicate returns `true` when rustle should DROP (not store) the
/// transfrag, matching StringTie.
///
/// Returns `false` (keep) for the cases StringTie's `else` branch DOES store/handle: short-read
/// or single-node (`n_inner_nodes <= 1`) transfrags, guides, and `csr_triggered` chimeric-suffix
/// folds. `rescue_opt_out` is the `RUSTLE_CHECKTRF_MULTINODE_RESCUE` env opt-out (`true` restores
/// the old store-it behavior).
fn checktrf_multinode_no_match_drop(
    is_shortread: bool,
    n_inner_nodes: usize,
    csr_triggered: bool,
    is_guide: bool,
    rescue_opt_out: bool,
) -> bool {
    !rescue_opt_out && !is_shortread && n_inner_nodes > 1 && !csr_triggered && !is_guide
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cargo test --release checktrf_gate_tests 2>&1 | tail -20`
Expected: `test result: ok. 6 passed; 0 failed`.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/path_extract.rs
git commit -m "feat(parity): add checktrf multinode drop predicate + unit tests

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 2: Add `SeedOutcome::ChecktrfMultinodeNoMatchDrop` variant and update exhaustive matches

The new outcome variant must be added to the enum AND to the three `match outcome` blocks that are exhaustive (no `_` arm). Two other matches (`pipeline.rs:15179`, `pipeline.rs:2007`) already have wildcard arms; `pipeline.rs:2007` is updated for correct checktrf accounting but would compile without it.

**Files:**
- Modify: `src/rustle/path_extract.rs:596-598` (enum) and `:11227` (match)
- Modify: `src/rustle/pipeline.rs:918`, `:7344`, `:2007` (matches)

- [ ] **Step 1: Add the enum variant**

In `src/rustle/path_extract.rs`, change (lines 595-598):

```rust
    /// Checktrf: rescue exons empty or too short or zero coverage.
    ChecktrfRescueFail,
    /// Successfully stored as transcript at output index.
    Stored(usize),
```

to:

```rust
    /// Checktrf: rescue exons empty or too short or zero coverage.
    ChecktrfRescueFail,
    /// Checktrf: multi-node long-read transfrag with no kept-path match — dropped to match
    /// StringTie (rlink.cpp:10369 redistribute-only branch; never independently stored).
    ChecktrfMultinodeNoMatchDrop,
    /// Successfully stored as transcript at output index.
    Stored(usize),
```

- [ ] **Step 2: Update the exhaustive match in `path_extract.rs:11227`**

Find the arm `SeedOutcome::ChecktrfRescueFail => ("checktrf_rescue_fail", None),` and add the new arm immediately after it:

```rust
                                SeedOutcome::ChecktrfRescueFail => ("checktrf_rescue_fail", None),
                                SeedOutcome::ChecktrfMultinodeNoMatchDrop => ("checktrf_multinode_drop", None),
```

- [ ] **Step 3: Update the exhaustive match in `pipeline.rs:918`**

Find `SeedOutcome::ChecktrfRescueFail => "checktrf_rescue_fail",` and add after it:

```rust
        SeedOutcome::ChecktrfRescueFail => "checktrf_rescue_fail",
        SeedOutcome::ChecktrfMultinodeNoMatchDrop => "checktrf_multinode_drop",
```

- [ ] **Step 4: Update the exhaustive match in `pipeline.rs:7344`**

Find `crate::path_extract::SeedOutcome::ChecktrfRescueFail => "checktrf_rescue_fail".into(),` and add after it:

```rust
                crate::path_extract::SeedOutcome::ChecktrfRescueFail => "checktrf_rescue_fail".into(),
                crate::path_extract::SeedOutcome::ChecktrfMultinodeNoMatchDrop => "checktrf_multinode_drop".into(),
```

- [ ] **Step 5: Update the checktrf-total accounting in `pipeline.rs:2007`**

Change the grouped arm:

```rust
            SeedOutcome::ChecktrfReadthr
            | SeedOutcome::ChecktrfRedistributed
            | SeedOutcome::ChecktrfIncomplete
            | SeedOutcome::ChecktrfRescueFail => {
                out.checktrf_total += 1;
            }
```

to include the new variant:

```rust
            SeedOutcome::ChecktrfReadthr
            | SeedOutcome::ChecktrfRedistributed
            | SeedOutcome::ChecktrfIncomplete
            | SeedOutcome::ChecktrfMultinodeNoMatchDrop
            | SeedOutcome::ChecktrfRescueFail => {
                out.checktrf_total += 1;
            }
```

- [ ] **Step 6: Build to verify all matches compile**

Run: `cargo build --release 2>&1 | grep -E "^error|Finished"`
Expected: `Finished \`release\` profile` with NO `error[E0004]` (non-exhaustive match) lines.

- [ ] **Step 7: Commit**

```bash
git add src/rustle/path_extract.rs src/rustle/pipeline.rs
git commit -m "feat(parity): add ChecktrfMultinodeNoMatchDrop seed outcome + match arms

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 3: Wire the gate at the independent-rescue block

Insert the gate at the top of the independent-rescue block (`path_extract.rs:9779`), before the existing eonly check at line 9786. All inputs are already in scope: `is_shortread_tf` (line 9540), `tf_nodes` (line 9483, inner nodes excluding source/sink), `csr_triggered` (line 9541), `transfrags[t].guide`.

**Files:**
- Modify: `src/rustle/path_extract.rs:9779-9790`

- [ ] **Step 1: Insert the gate**

The current code at lines 9779-9790 is:

```rust
            {
                // Independent rescue: store complete transfrag as its own prediction.
                // Applies to: shortread/<=1-node, AND multi-node longread with no kept-path match.
                // `else if(!eonly || guide)` — in non-eonly mode always rescues.
                // Non-contiguous jumps must have explicit transfrag edge pattern support.

                // in eonly mode, only guide transfrags are independently rescued.
                if config.eonly && !transfrags[t].guide {
                    record_outcome!(t, SeedOutcome::ChecktrfEonlySkip);
                    emit_checktrf_result!(t, "eonly_skip", &tf_nodes);
                    continue;
                }
```

Replace it with (note the corrected comment and the new gate inserted before the eonly check):

```rust
            {
                // Independent rescue: store complete transfrag as its own prediction.
                // Applies to: shortread/<=1-node only (StringTie rlink.cpp:10413 `else` branch).
                // `else if(!eonly || guide)` — in non-eonly mode always rescues.
                // Non-contiguous jumps must have explicit transfrag edge pattern support.

                // ST-faithful checktrf gate (rlink.cpp:10369/10413): a multi-node long-read
                // transfrag with no kept-path match is redistribute-only in StringTie and is
                // NEVER stored as an independent prediction. Drop it here to match ST. Preserves
                // short-read/single-node rescue, guides, and csr_triggered chimeric-suffix folds.
                // Opt-out: RUSTLE_CHECKTRF_MULTINODE_RESCUE=1 restores the old store-it behavior.
                if checktrf_multinode_no_match_drop(
                    is_shortread_tf,
                    tf_nodes.len(),
                    csr_triggered,
                    transfrags[t].guide,
                    std::env::var_os("RUSTLE_CHECKTRF_MULTINODE_RESCUE").is_some(),
                ) {
                    record_outcome!(t, SeedOutcome::ChecktrfMultinodeNoMatchDrop);
                    emit_checktrf_result!(t, "multinode_no_match_drop", &tf_nodes);
                    continue;
                }

                // in eonly mode, only guide transfrags are independently rescued.
                if config.eonly && !transfrags[t].guide {
                    record_outcome!(t, SeedOutcome::ChecktrfEonlySkip);
                    emit_checktrf_result!(t, "eonly_skip", &tf_nodes);
                    continue;
                }
```

- [ ] **Step 2: Build**

Run: `cargo build --release 2>&1 | grep -E "^error|Finished"`
Expected: `Finished \`release\` profile`, no errors.

- [ ] **Step 3: Smoke-check that the gate fires**

Run (regenerates the rustle parity log and counts the new drop events):

```bash
RUSTLE_PARITY_LOG=/tmp/ruP.jsonl RUSTLE_PARITY_FILTER_STEPS=checktrf_result \
  ./target/release/rustle -L GGO_19.bam -o /tmp/ruP.gtf 2>/dev/null
grep -c "multinode_no_match_drop" /tmp/ruP.jsonl
```

Expected: a count `>= 48` (the gate fires on the multi-node long-read no-match transfrags; some loci have more than one).

- [ ] **Step 4: Commit**

```bash
git add src/rustle/path_extract.rs
git commit -m "feat(parity): gate checktrf multinode-longread store to match StringTie

Drops multi-node long-read transfrags with no kept-path match instead of
storing them as independent predictions (StringTie rlink.cpp:10369 is
redistribute-only). Default-on; opt-out RUSTLE_CHECKTRF_MULTINODE_RESCUE=1.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 4: Integration validation (parity gate, TP/FP, opt-out byte-identical)

Commit a reproducible TP/FP classifier and verify the four chr19-level success criteria from the spec.

**Files:**
- Create: `bench/classify_checktrf_tpfp.py`

- [ ] **Step 1: Create the TP/FP classifier**

Create `bench/classify_checktrf_tpfp.py`:

```python
#!/usr/bin/env python3
"""Classify rustle-only intron chains (vs StringTie) as TP/FP against the real RefSeq
annotation, split by extraction source (flow vs checktrf). Verifies the checktrf store-gate
fix impact. Usage: python3 bench/classify_checktrf_tpfp.py <ru.gtf> <st.gtf> <ru.jsonl> <annot.gff> <chrom>
"""
import re, sys, json, collections

def chains_gtf(path):
    tx = collections.defaultdict(list); strand = {}
    for line in open(path):
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon': continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m: continue
        tx[m.group(1)].append((int(f[3]), int(f[4]))); strand[m.group(1)] = f[6]
    out = set()
    for t, ex in tx.items():
        ex.sort()
        if len(ex) < 2: continue
        out.add((strand[t], tuple((ex[i][1], ex[i+1][0]) for i in range(len(ex)-1))))
    return out

def chains_gff(path, chrom):
    tx = collections.defaultdict(list); strand = {}
    for line in open(path):
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != chrom or f[2] != 'exon': continue
        m = re.search(r'Parent=([^;]+)', f[8])
        if not m: continue
        tx[m.group(1)].append((int(f[3]), int(f[4]))); strand[m.group(1)] = f[6]
    out = set()
    for t, ex in tx.items():
        ex.sort()
        if len(ex) < 2: continue
        out.add((strand[t], tuple((ex[i][1], ex[i+1][0]) for i in range(len(ex)-1))))
    return out

def safe(p):
    for l in open(p):
        try: yield json.loads(l)
        except ValueError: continue

ru_gtf, st_gtf, ru_jsonl, annot, chrom = sys.argv[1:6]
ru = chains_gtf(ru_gtf); st = chains_gtf(st_gtf); ann = chains_gff(annot, chrom)
ru_only = ru - st
parse = lambda s: tuple((int(a), int(b)) for a, b in re.findall(r'(\d+)-(\d+)', s))
src = collections.defaultdict(set)
for d in safe(ru_jsonl):
    if d.get('step') != 'path_extracted': continue
    pl = d['payload']; ch = tuple((a-1, b+1) for a, b in parse(pl.get('introns', '')))
    src[(d['strand'], ch)].add(pl.get('source', '?'))

res = collections.Counter()
for k in ru_only:
    s = src.get(k, set())
    has_ck = any('checktrf' in x for x in s); has_flow = any(x == 'flow' for x in s)
    bucket = 'checktrf_only' if (has_ck and not has_flow) else \
             'checktrf+flow' if (has_ck and has_flow) else 'flow_only'
    res[(bucket, 'TP' if k in ann else 'FP')] += 1
print(f"rustle-only chains: {len(ru_only)}  (annotation chains: {len(ann)})")
for b in ('flow_only', 'checktrf_only', 'checktrf+flow'):
    print(f"  {b:14s}: TP={res[(b,'TP')]:3d}  FP={res[(b,'FP')]:3d}")
```

- [ ] **Step 2: Regenerate the StringTie reference trace (deterministic ground truth)**

Run:

```bash
STRINGTIE_PARITY_LOG=/tmp/stP.jsonl tools/stringtie/stringtie -L GGO_19.bam -o /tmp/stP.gtf 2>/dev/null
RUSTLE_PARITY_LOG=/tmp/ruP.jsonl RUSTLE_PARITY_FILTER_STEPS=parse_trflong_seed,path_extracted \
  ./target/release/rustle -L GGO_19.bam -o /tmp/ruP.gtf 2>/dev/null
```

- [ ] **Step 3: Run the primary parity gate**

Run: `python3 bench/gtf_chain_diff.py /tmp/ruP.gtf /tmp/stP.gtf`
Expected: `Rustle-only: 139  ST-only: 104` (rustle-only drops from 187 to 139; **ST-only stays 104** — if ST-only rose, the gate over-rejected, STOP and investigate).

- [ ] **Step 4: Run the TP/FP classifier**

Run: `python3 bench/classify_checktrf_tpfp.py /tmp/ruP.gtf /tmp/stP.gtf /tmp/ruP.jsonl ../GGO_genomic.gff NC_073243.2`
Expected: `checktrf_only: TP=0  FP=0` (the 48 checktrf-only chains are gone), `flow_only` unchanged at `TP=30 FP=107`, `checktrf+flow: TP=0 FP=2` (survive via flow). Net vs the pre-fix run: −46 FP, −2 TP.

- [ ] **Step 5: Verify opt-out is byte-identical to old behavior**

Run:

```bash
RUSTLE_CHECKTRF_MULTINODE_RESCUE=1 ./target/release/rustle -L GGO_19.bam -o /tmp/ru_optout.gtf 2>/dev/null
python3 bench/gtf_chain_diff.py /tmp/ru_optout.gtf /tmp/stP.gtf
```

Expected: `Rustle-only: 187  ST-only: 104` (opt-out restores the pre-fix output).

- [ ] **Step 6: 2-TP safety check (do not erode a shipped recall win)**

Run (lists the 2 real-annotation TP chains the gate drops, so they can be checked against shipped recovery wins):

```bash
python3 - <<'PY'
import re, json, collections
def chains_gtf(p):
    tx=collections.defaultdict(list); s={}
    for l in open(p):
        if l.startswith('#'): continue
        f=l.rstrip().split('\t')
        if len(f)<9 or f[2]!='exon': continue
        m=re.search(r'transcript_id "([^"]+)"',f[8])
        if not m: continue
        tx[m.group(1)].append((int(f[3]),int(f[4]))); s[m.group(1)]=f[6]
    o=set()
    for t,ex in tx.items():
        ex.sort()
        if len(ex)<2: continue
        o.add((s[t],tuple((ex[i][1],ex[i+1][0]) for i in range(len(ex)-1))))
    return o
def chains_gff(p,c):
    tx=collections.defaultdict(list); s={}
    for l in open(p):
        if l.startswith('#'): continue
        f=l.rstrip().split('\t')
        if len(f)<9 or f[0]!=c or f[2]!='exon': continue
        m=re.search(r'Parent=([^;]+)',f[8])
        if not m: continue
        tx[m.group(1)].append((int(f[3]),int(f[4]))); s[m.group(1)]=f[6]
    o=set()
    for t,ex in tx.items():
        ex.sort()
        if len(ex)<2: continue
        o.add((s[t],tuple((ex[i][1],ex[i+1][0]) for i in range(len(ex)-1))))
    return o
def safe(p):
    for l in open(p):
        try: yield json.loads(l)
        except ValueError: continue
ru=chains_gtf('/tmp/ru_optout.gtf'); st=chains_gtf('/tmp/stP.gtf')
ann=chains_gff('../GGO_genomic.gff','NC_073243.2')
parse=lambda s: tuple((int(a),int(b)) for a,b in re.findall(r'(\d+)-(\d+)',s))
src=collections.defaultdict(set)
for d in safe('/tmp/ruP.jsonl'):
    if d.get('step')!='path_extracted': continue
    pl=d['payload']; ch=tuple((a-1,b+1) for a,b in parse(pl.get('introns','')))
    src[(d['strand'],ch)].add(pl.get('source','?'))
ru_only=ru-st
for k in ru_only:
    s=src.get(k,set())
    if any('checktrf' in x for x in s) and not any(x=='flow' for x in s) and k in ann:
        st_,ch=k
        print(f"TP-at-risk {st_} {ch[0][0]}..{ch[-1][1]} ({len(ch)} introns)")
PY
```

Expected: exactly 2 lines. Manually confirm neither chain is a strand-bundling (`RUSTLE_STRAND_PURE_MINORITY`) or read-chain recovery win by checking they are not produced only under those modes (re-run with `RUSTLE_STRAND_PURE_MINORITY_OFF=1` and confirm the same 2 chains were already checktrf-sourced, not strand-recovered). Record the 2 coordinates in the commit message. If either IS a shipped win, STOP — the guard needs a carve-out.

- [ ] **Step 7: Commit**

```bash
git add bench/classify_checktrf_tpfp.py
git commit -m "test(parity): add checktrf TP/FP classifier + record gate validation

chr19: rustle-only 187->139, ST-only 104 unchanged, -46 FP / -2 TP vs RefSeq.
opt-out RUSTLE_CHECKTRF_MULTINODE_RESCUE=1 byte-identical.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 5: Genome-wide validation (gating, since default-on)

The change ships default-on, so it must be net-positive genome-wide, not just on chr19. This task RUNS rustle on the whole-genome BAM (`../GGO.bam`, baseline `-L` — no `--vg`, so no OOM concern) fix-on vs opt-out and compares against the real RefSeq annotation. It changes no source. It is long-running (whole-genome `-L` on a 1.7 GB BAM) and may be handed to the user/controller to run rather than blocking a subagent.

**Files:**
- None (validation only). Records result in the final findings note (Task 6).

- [ ] **Step 1: Run whole-genome `-L`, fix-on vs opt-out**

Run from the repo root:

```bash
# fix ON (default)
./target/release/rustle -L -p 8 ../GGO.bam -o /tmp/gw_fixon.gtf 2>/tmp/gw_fixon.err
# fix OFF (opt-out) for the paired comparison
RUSTLE_CHECKTRF_MULTINODE_RESCUE=1 ./target/release/rustle -L -p 8 ../GGO.bam -o /tmp/gw_fixoff.gtf 2>/tmp/gw_fixoff.err
```

(If memory pressure appears, run per-chromosome instead using the `bench/flow_recall_phase0/perchrom/<chrom>/c.bam` slices, looping over chromosomes; the comparison in Step 2 is per-GTF and works on any subset.)

- [ ] **Step 2: Compute the net FSM / precision delta vs RefSeq**

Run gffcompare for both against the real annotation and compare:

```bash
gffcompare -r ../GGO_genomic.gff -o /tmp/gw_fixon_cmp  /tmp/gw_fixon.gtf  2>/dev/null
gffcompare -r ../GGO_genomic.gff -o /tmp/gw_fixoff_cmp /tmp/gw_fixoff.gtf 2>/dev/null
echo "=== FIX ON ===";  grep -E "Transcript level|Intron chain" /tmp/gw_fixon_cmp.stats
echo "=== FIX OFF ==="; grep -E "Transcript level|Intron chain" /tmp/gw_fixoff_cmp.stats
```

Expected: fix-on shows **higher precision** with intron-chain sensitivity nearly flat (a small FSM/TP cost), F1-positive overall (the chr19 ratio was 23:1 FP:TP removed).
**Abort criterion:** if genome-wide intron-chain sensitivity drops materially (more than the ~chr19-scaled −2 TP per-locus rate) OR precision does not improve, do NOT keep default-on — flip to opt-in (invert the env to default-off, enabling with a new `RUSTLE_CHECKTRF_MULTINODE_DROP=1`) and note it in Task 6.

- [ ] **Step 3: Record the numbers (no commit; feeds Task 6)**

Save the fix-on/fix-off intron-chain Sn/Pr and the transcript-level counts for the findings note.

---

### Task 6: Findings note + memory update

**Files:**
- Modify: `docs/STRINGTIE_PARITY_FINDINGS.md` (append §6o)
- Modify: `/home/juanfra/.claude/projects/-mnt-c-Users-jfris-Desktop/memory/project_post_flow_gate_pin.md` and `MEMORY.md`

- [ ] **Step 1: Append the findings note**

Append a `## §6o` section to `docs/STRINGTIE_PARITY_FINDINGS.md` recording: the divergence (rlink.cpp:10369/10413 vs path_extract.rs:9779), the fix (default-on gate, opt-out `RUSTLE_CHECKTRF_MULTINODE_RESCUE`), the chr19 result (rustle-only 187→139, −46 FP / −2 TP, ST-only unchanged), and the genome-wide result from Task 5. State explicitly that the flow-139 pool remains open and is scoped separately next.

- [ ] **Step 2: Update memory**

Update `project_post_flow_gate_pin.md` with the shipped fix + result, and refresh its `MEMORY.md` index line.

- [ ] **Step 3: Commit**

```bash
git add docs/STRINGTIE_PARITY_FINDINGS.md
git commit -m "docs(parity): record checktrf multinode store-gate fix (§6o)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

## Notes for the implementer

- **Build time:** `cargo build --release` is ~2-3 minutes. Run it in the background if the harness supports it.
- **Trace JSONL gotcha:** the instrumented StringTie log (`/tmp/stP.jsonl`) contains ~0.7% malformed lines; always parse with a `try/except json.JSONDecodeError: continue` guard (the classifier already does).
- **`/tmp` is cleared between sessions** — Tasks 4-5 regenerate the traces from the BAM; do not assume `/tmp/stP.gtf` etc. persist.
- **Do NOT touch the flow loop.** The flow-139 pool (30 TP / 107 FP) is out of scope; four depletion/ordering/coverage hypotheses were already refuted and a downstream coverage fix over-killed 119 TP (§6n). Keep changes strictly to the checktrf store gate.
- **`tf_nodes` is inner nodes** (source/sink filtered, path_extract.rs:9483), so `tf_nodes.len() > 1` correctly means multi-exon, matching StringTie's `nodes.Count()>1`.
- **Do NOT zero `transfrags[t].abundance` in the gate.** The gate only records the outcome, emits the trace, and `continue`s — leaving abundance untouched. StringTie zeroes abundance only inside the *matched* redistribution branch (rlink.cpp:~10410); on a no-match multi-node transfrag it leaves abundance as-is. Zeroing it could starve a later rustle post-pass. This is a deliberate spec decision (Component 2).
