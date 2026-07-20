# Family Variation Graph as the O2 Assignment Substrate — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the per-family variation graph the explicit substrate the O2 copy-assignment decision runs on — an auditable ad-hoc reference (copies = paths) each read is threaded through — with **bit-identical** assignments.

**Architecture:** Add a `BubbleGraph` built once per family from the copy profiles; refactor the per-read likelihood (`read_copy_evidence`) to thread the read as a walk through it (same arithmetic → identical `logl`/`min_p`/`n_decisive`); keep the significance gate and EM untouched; add an opt-in per-read provenance audit emit (certificate tags on the read `W`-lines).

**Tech Stack:** Rust (`src/rustle/vg_family/`), `cargo test`. Byte-identity is the safety property, guarded by a golden characterization test captured before the refactor.

## Global Constraints

- **Bit-identical `.assignments.tsv`** (and EM quant): the decision arithmetic and PSV-column order are preserved exactly; the graph is a re-expression, not a new algorithm.
- The significance certificate / `assign_read_editing` gate and the EM are **not modified**; `read_copy_evidence` is shared by both and must stay bit-identical for both.
- The `BubbleGraph` is built **once per family** and threaded into `read_copy_evidence`/`assign_read_editing`; `assign_read` (the wrapper) builds it internally so existing call sites/tests keep their signatures.
- The audit emit is **opt-in and additive** (read `W`-line tags only when a cert is attached); default `--phase` / `.assignments.tsv` output is unchanged and not slowed.
- Out of scope: any number change, gate/EM logic, O1, mapping/recovery, `assemble_vg` unification.
- Validation runs of `copy_assign` are **foreground + serial + `winloci_scratch`** (WSL2 crash rule); no `copy_assign` background/nohup/`pkill`.

---

### Task 1: `BubbleGraph` type + `from_copies`

**Files:**
- Modify: `src/rustle/vg_family/copy_assign.rs` (add the types near `CopyProfile`, ~line 30)

**Interfaces:**
- Consumes: `CopyProfile { copy_id, alleles: Vec<Option<u8>>, junctions }` (existing).
- Produces:
  - `pub struct Bubble { pub col: usize, pub copy_allele: Vec<Option<u8>>, pub decisive: bool }`
  - `pub struct BubbleGraph { pub bubbles: Vec<Bubble>, pub n_copies: usize }`
  - `BubbleGraph::from_copies(copies: &[CopyProfile]) -> BubbleGraph`

- [ ] **Step 1: Write the failing test**

Add to the `#[cfg(test)] mod tests` in `src/rustle/vg_family/copy_assign.rs`:

```rust
#[test]
fn bubble_graph_from_copies_structure() {
    // 3 copies over 3 columns: col0 all 'A' (not decisive); col1 A/C/A (decisive); col2 copy2 gap
    let mk = |a: Vec<Option<u8>>| CopyProfile { copy_id: 0, alleles: a, junctions: vec![] };
    let copies = vec![
        mk(vec![Some(b'A'), Some(b'A'), Some(b'G')]),
        mk(vec![Some(b'A'), Some(b'C'), Some(b'G')]),
        mk(vec![Some(b'A'), Some(b'A'), None]),
    ];
    let g = BubbleGraph::from_copies(&copies);
    assert_eq!(g.n_copies, 3);
    assert_eq!(g.bubbles.len(), 3);
    assert_eq!(g.bubbles[0].col, 0);
    assert!(!g.bubbles[0].decisive);                 // all 'A'
    assert!(g.bubbles[1].decisive);                  // A vs C
    assert!(!g.bubbles[2].decisive);                 // G, G, gap -> one distinct allele
    assert_eq!(g.bubbles[1].copy_allele, vec![Some(b'A'), Some(b'C'), Some(b'A')]);
    assert_eq!(g.bubbles[2].copy_allele, vec![Some(b'G'), Some(b'G'), None]);
}
```

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test --lib bubble_graph_from_copies_structure`
Expected: FAIL to compile — `cannot find type BubbleGraph`.

- [ ] **Step 3: Implement the types**

Add to `src/rustle/vg_family/copy_assign.rs` (just after the `CopyProfile` struct):

```rust
/// One PSV bubble of the family's ad-hoc-reference variation graph.
#[derive(Clone, Debug)]
pub struct Bubble {
    /// PSV column index (= bubble id). Bubbles are in ascending column order (== the matrix order).
    pub col: usize,
    /// Allele-node each copy PATH visits here (index = copy). `None` = the copy has a gap.
    pub copy_allele: Vec<Option<u8>>,
    /// The copies carry >= 2 distinct non-`None` alleles here (read-independent). Matches the old
    /// per-column `differ` test in `read_copy_evidence`.
    pub decisive: bool,
}

/// The per-family variation graph the O2 decision threads reads through: one bubble per PSV column,
/// each copy a path over the allele-nodes. Built once per family; the ad-hoc, auditable reference.
#[derive(Clone, Debug)]
pub struct BubbleGraph {
    pub bubbles: Vec<Bubble>,
    pub n_copies: usize,
}

impl BubbleGraph {
    /// Build the family's bubble graph from the copy profiles. Deterministic; `decisive` is computed
    /// exactly as `read_copy_evidence`'s inner `differ` loop (>= 2 distinct non-`None` alleles).
    pub fn from_copies(copies: &[CopyProfile]) -> BubbleGraph {
        let n_cols = copies.iter().map(|c| c.alleles.len()).max().unwrap_or(0);
        let mut bubbles = Vec::with_capacity(n_cols);
        for col in 0..n_cols {
            let copy_allele: Vec<Option<u8>> =
                copies.iter().map(|c| c.alleles.get(col).copied().flatten()).collect();
            let mut seen: Option<u8> = None;
            let mut decisive = false;
            for a in copy_allele.iter().flatten() {
                match seen {
                    None => seen = Some(*a),
                    Some(s) => {
                        if s != *a {
                            decisive = true;
                        }
                    }
                }
            }
            bubbles.push(Bubble { col, copy_allele, decisive });
        }
        BubbleGraph { bubbles, n_copies: copies.len() }
    }
}
```

- [ ] **Step 4: Run to verify it passes**

Run: `cargo test --lib bubble_graph_from_copies_structure`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/copy_assign.rs
git commit -m "feat(vg-o2): BubbleGraph ad-hoc-reference type + from_copies"
```

---

### Task 2: Golden `ReadEvidence` characterization test (safety baseline, captured BEFORE the refactor)

**Files:**
- Modify: `src/rustle/vg_family/copy_assign.rs` (add a test in `mod tests`)

**Interfaces:**
- Consumes: the CURRENT `read_copy_evidence(read, copies, p, editing_cols)` (unchanged in this task) — this test freezes its exact `f64` output so Task 3's refactor is provably bit-identical.

- [ ] **Step 1: Write a characterization test that PASSES on the current code**

Add to `mod tests` in `src/rustle/vg_family/copy_assign.rs`:

```rust
// Golden fixture: freezes the EXACT current ReadEvidence so the Task-3 graph refactor is provably
// bit-identical. If this ever changes, the decision changed — which this project forbids.
fn golden_fixture() -> (ReadFeatures, Vec<CopyProfile>, AssignParams) {
    let cp = |a: Vec<Option<u8>>, j: Vec<i64>| CopyProfile { copy_id: 0, alleles: a, junctions: j };
    let copies = vec![
        cp(vec![Some(b'A'), Some(b'C'), Some(b'G'), Some(b'T')], vec![100]),
        cp(vec![Some(b'A'), Some(b'G'), Some(b'G'), None],       vec![]),
        cp(vec![Some(b'A'), Some(b'C'), Some(b'A'), Some(b'T')], vec![100, 250]),
    ];
    let read = ReadFeatures {
        psv_obs: vec![Some(b'A'), Some(b'C'), None, Some(b'T')],
        psv_qual: vec![Some(30), None, None, Some(20)],
        junctions: vec![100],
    };
    (read, copies, AssignParams::default())
}

#[test]
fn read_copy_evidence_golden_is_stable() {
    let (read, copies, p) = golden_fixture();
    let ev = read_copy_evidence(&read, &copies, &p, &[]);
    // Capture the CURRENT bit-exact values: run once, read the assert-failure message, paste them back,
    // then this test freezes them. (Use {:?} on f64 for the exact decimal that round-trips.)
    // EXPECTED (fill from the first run, then it must never change):
    let want_logl: Vec<f64> = EXPECT_LOGL.to_vec();
    let want_min_p: f64 = EXPECT_MIN_P;
    let want_n_decisive: usize = EXPECT_N_DECISIVE;
    assert_eq!(ev.logl, want_logl, "logl drifted -> decision changed");
    assert_eq!(ev.min_p, want_min_p, "min_p drifted");
    assert_eq!(ev.n_decisive, want_n_decisive, "n_decisive drifted");
}
```

Add the frozen constants near the test (values are filled in Step 3):
```rust
const EXPECT_LOGL: [f64; 3] = [/* filled in Step 3 */];
const EXPECT_MIN_P: f64 = /* filled in Step 3 */;
const EXPECT_N_DECISIVE: usize = /* filled in Step 3 */;
```

- [ ] **Step 2: Run to capture the current values**

Run: `cargo test --lib read_copy_evidence_golden_is_stable -- --nocapture`
Expected: FAIL to compile (placeholders) or FAIL the assert. Temporarily add `eprintln!("{:?} {:?} {}", ev.logl, ev.min_p, ev.n_decisive);` before the asserts, run again, and record the exact printed `logl` vector, `min_p`, and `n_decisive`.

- [ ] **Step 3: Freeze the captured values**

Replace the three `const` placeholders with the exact values printed in Step 2 (copy the `f64` debug literals verbatim, e.g. `const EXPECT_LOGL: [f64; 3] = [-1.2039728043259361, -8.940..., -1.203...];`). Remove the temporary `eprintln!`.

- [ ] **Step 4: Run to verify it now PASSES on the unchanged code**

Run: `cargo test --lib read_copy_evidence_golden_is_stable`
Expected: PASS. This is the frozen baseline; Task 3 must keep it green.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/copy_assign.rs
git commit -m "test(vg-o2): freeze golden ReadEvidence baseline before the graph refactor"
```

---

### Task 3: Thread `read_copy_evidence` through the `BubbleGraph` (bit-identical) + update callers

**Files:**
- Modify: `src/rustle/vg_family/copy_assign.rs:267` (`read_copy_evidence`), `:367` (`assign_read_editing`), `:358` (`assign_read`)
- Modify: `src/rustle/vg_family/copy_assign_pipeline.rs:1496,1502` (per-family caller); `src/rustle/vg_family/em_copy_assign.rs:254,525,589` (EM callers)

**Interfaces:**
- Consumes: `BubbleGraph::from_copies` (Task 1); the golden test (Task 2).
- Produces (new signatures):
  - `read_copy_evidence(read: &ReadFeatures, graph: &BubbleGraph, copies: &[CopyProfile], p: &AssignParams, editing_cols: &[bool]) -> ReadEvidence`
  - `assign_read_editing(read: &ReadFeatures, graph: &BubbleGraph, copies: &[CopyProfile], p: &AssignParams, editing_cols: &[bool]) -> Option<Assignment>`
  - `assign_read(read: &ReadFeatures, copies: &[CopyProfile], p: &AssignParams) -> Option<Assignment>` (unchanged signature — builds the graph internally).

- [ ] **Step 1: Update the golden test call + add the graph-vs-inline parity assertion**

In `read_copy_evidence_golden_is_stable`, change the call to build and pass the graph, and assert the graph path equals the golden values (the whole point):
```rust
    let (read, copies, p) = golden_fixture();
    let graph = BubbleGraph::from_copies(&copies);
    let ev = read_copy_evidence(&read, &graph, &copies, &p, &[]);
```
Also update `read_copy_evidence_matches_assignment_internals` (`copy_assign.rs:491`) and the `read_copy_evidence`/`assign_read_editing` test call sites (`:499`, `:503`) to build a `graph` and pass it.

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test --lib read_copy_evidence_golden_is_stable`
Expected: FAIL to compile — `read_copy_evidence` takes 4 args, not 5.

- [ ] **Step 3: Refactor `read_copy_evidence` to thread the graph**

Replace the PSV-term block of `read_copy_evidence` (`copy_assign.rs:267`) so the signature gains `graph: &BubbleGraph` and the column loop walks the bubbles (the junction term, argmax, and `min_p` blocks are UNCHANGED):

```rust
pub(crate) fn read_copy_evidence(
    read: &ReadFeatures,
    graph: &BubbleGraph,
    copies: &[CopyProfile],
    p: &AssignParams,
    editing_cols: &[bool],
) -> ReadEvidence {
    let n = copies.len();
    let mut logl = vec![0.0f64; n];
    let mut n_decisive = 0usize;

    // --- PSV term: thread the read as a WALK through the bubbles it spans (== the old `spanned` loop,
    //     same columns, same order, so the f64 sums are bit-identical). ---
    for b in &graph.bubbles {
        let obs = match read.psv_obs.get(b.col).copied().flatten() {
            Some(o) => o,
            None => continue, // read does not span this bubble
        };
        let e = match read.psv_qual.get(b.col).copied().flatten() {
            Some(q) => super::copy_split::phred_err(q),
            None => p.error_rate,
        };
        let lp_match = (1.0 - e).ln();
        let lp_mis = (e / 3.0).ln();
        if b.decisive {
            n_decisive += 1;
        }
        for (ci, a) in b.copy_allele.iter().enumerate() {
            if let Some(a) = a {
                logl[ci] += if obs == *a { lp_match } else { lp_mis };
            }
        }
    }

    // --- junction (copy-specific splicing) term --- UNCHANGED (still reads `copies`)
    for &jb in &read.junctions {
        let present: Vec<bool> = copies
            .iter()
            .map(|c| boundary_present(jb, &c.junctions, p.boundary_tol))
            .collect();
        let np = present.iter().filter(|&&x| x).count();
        if np > 0 && np < n {
            n_decisive += 1;
        }
        for ci in 0..n {
            logl[ci] += if present[ci] { p.junction_weight } else { -p.junction_weight };
        }
    }

    // --- argmax + min_p (IsoCon identifiability bound) --- UNCHANGED
    let mut best = 0usize;
    for i in 1..n {
        if logl[i] > logl[best] {
            best = i;
        }
    }
    let mut min_p = 0.0f64;
    for c in 0..n {
        if c == best {
            continue;
        }
        let (_pbc, attain) = copy_pair_significance(read, &copies[best], &copies[c], p, editing_cols);
        if attain > min_p {
            min_p = attain;
        }
    }
    ReadEvidence { logl, min_p, n_decisive }
}
```

- [ ] **Step 4: Thread the graph through `assign_read_editing` and `assign_read`**

`assign_read_editing` gains `graph: &BubbleGraph` and forwards it:
```rust
pub fn assign_read_editing(
    read: &ReadFeatures,
    graph: &BubbleGraph,
    copies: &[CopyProfile],
    p: &AssignParams,
    editing_cols: &[bool],
) -> Option<Assignment> {
    let n = copies.len();
    if n == 0 {
        return None;
    }
    let ev = read_copy_evidence(read, graph, copies, p, editing_cols);
    // ... rest of the function UNCHANGED ...
```
`assign_read` builds the graph once and forwards (signature unchanged, so its many callers/tests are untouched):
```rust
pub fn assign_read(read: &ReadFeatures, copies: &[CopyProfile], p: &AssignParams) -> Option<Assignment> {
    let graph = BubbleGraph::from_copies(copies);
    assign_read_editing(read, &graph, copies, p, &[])
}
```

- [ ] **Step 5: Update the two real callers to build the graph once per family**

In `src/rustle/vg_family/copy_assign_pipeline.rs`, build the graph once before the read loop that contains lines 1496/1502 and pass `&graph`:
```rust
    let graph = super::copy_assign::BubbleGraph::from_copies(&fp.profiles);
    // ... inside the per-read loop:
    let Some(combined) = assign_read_editing(&feats, &graph, &fp.profiles, p, &editing_cols) else { ... };
    // ...
    let Some(psv) = assign_read_editing(&psv_feats, &graph, &fp.profiles, p, &editing_cols) else { ... };
```
In `src/rustle/vg_family/em_copy_assign.rs`, build the graph once from `copies` in each of the three scopes (lines 254, 525, 589) before the `.map(|..| read_copy_evidence(...))` and pass `&graph`:
```rust
    let graph = super::copy_assign::BubbleGraph::from_copies(&copies);
    // ... read_copy_evidence(r, &graph, &copies, &params, &editing_cols) ...
```

- [ ] **Step 6: Run the golden test + the full O2 suite (must be bit-identical)**

Run: `cargo test --lib copy_assign`
Then: `cargo test --lib em_copy_assign`
Expected: PASS, including `read_copy_evidence_golden_is_stable` (the frozen values are reproduced through the graph) — proving bit-identity.

- [ ] **Step 7: Commit**

```bash
git add src/rustle/vg_family/copy_assign.rs src/rustle/vg_family/copy_assign_pipeline.rs src/rustle/vg_family/em_copy_assign.rs
git commit -m "refactor(vg-o2): thread read_copy_evidence through the BubbleGraph (bit-identical)"
```

---

### Task 4: Per-read provenance — certificate tags on the read `W`-lines (opt-in, additive)

**Files:**
- Modify: `src/rustle/vg_family/copy_graph.rs:76` (`ReadWalk`), `:209-226` (W-line emit)

**Interfaces:**
- Consumes: `AssignStatus` (from `copy_assign`), the per-read `Assignment` (`p_value`, `min_p_value`, `status`).
- Produces: `pub struct ReadCert { pub p_value: f64, pub min_p_value: f64, pub status: AssignStatus }` and a new `pub cert: Option<ReadCert>` field on `ReadWalk`. When `Some`, the read `W`-line carries `CP:Z:<copy>  PV:f:<p>  MP:f:<min_p>  ST:Z:<status>`.

- [ ] **Step 1: Write the failing test**

Add to `#[cfg(test)] mod tests` in `src/rustle/vg_family/copy_graph.rs`:

```rust
#[test]
fn read_walk_cert_tags_emitted_when_present() {
    use super::super::copy_assign::AssignStatus;
    let mut g = small_two_copy_graph(); // existing test helper that builds a CopyGraph w/ >=1 bubble
    g.reads = vec![
        ReadWalk { name: "r1".into(), obs: vec![Some(b'A'), Some(b'C')], assigned_copy: Some(0),
                   cert: Some(ReadCert { p_value: 0.001, min_p_value: 0.0005, status: AssignStatus::Assigned }) },
        ReadWalk { name: "r2".into(), obs: vec![Some(b'A'), Some(b'C')], assigned_copy: None, cert: None },
    ];
    let gfa = g.to_gfa();
    let r1 = gfa.lines().find(|l| l.starts_with("W\tr1")).expect("r1 walk");
    assert!(r1.contains("CP:Z:copy0") && r1.contains("PV:f:0.001") && r1.contains("MP:f:0.0005") && r1.contains("ST:Z:Assigned"));
    let r2 = gfa.lines().find(|l| l.starts_with("W\tr2")).expect("r2 walk");
    assert!(!r2.contains("CP:Z") && !r2.contains("PV:f"), "no cert -> no tags (backward-compatible)");
}
```
(If no `small_two_copy_graph` helper exists, build the `CopyGraph` inline as the existing W-line test at `copy_graph.rs:607` does — reuse that fixture and add `cert: None` to its `ReadWalk`s.)

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test --lib read_walk_cert_tags_emitted_when_present`
Expected: FAIL to compile — `ReadWalk` has no field `cert` / no `ReadCert`.

- [ ] **Step 3: Add `ReadCert` + the `cert` field + emit the tags**

In `src/rustle/vg_family/copy_graph.rs`, add the struct and field:
```rust
/// Per-read significance certificate carried onto the audit W-line.
#[derive(Clone, Debug)]
pub struct ReadCert {
    pub p_value: f64,
    pub min_p_value: f64,
    pub status: super::copy_assign::AssignStatus,
}
```
and add `pub cert: Option<ReadCert>,` to `ReadWalk` (after `assigned_copy`). Then change the W-line emit (`:225`) to append the tags when a cert is present:
```rust
            let hap = r.assigned_copy.map(|c| c as i64).unwrap_or(-1).max(0);
            let mut line = format!("W\t{}\t{}\t{}\t0\t{}\t{}", r.name, hap, self.family, toks.len(), w);
            if let Some(c) = &r.cert {
                use super::copy_assign::AssignStatus::*;
                let st = match c.status { Assigned => "Assigned", Ambiguous => "Ambiguous", Tied => "Tied" };
                let cp = r.assigned_copy.map(|c| format!("copy{c}")).unwrap_or_else(|| "none".into());
                line.push_str(&format!("\tCP:Z:{cp}\tPV:f:{}\tMP:f:{}\tST:Z:{st}", c.p_value, c.min_p_value));
            }
            out.walks.push(line);
```
Add `cert: None` to EVERY existing `ReadWalk { .. }` construction in the crate (the copy_graph tests at ~`:607`/`:628`, and the `--phase` build path in `src/bin/copy_assign.rs` where read walks are assembled — grep `ReadWalk {` and add the field). In the `src/bin/copy_assign.rs` `--phase` path, populate `cert` from the read's `Assignment` (`Some(ReadCert { p_value: a.p_value, min_p_value: a.min_p_value, status: a.status })`) so a real audit run carries the certificates.

- [ ] **Step 4: Run to verify it passes**

Run: `cargo test --lib copy_graph`
Expected: PASS (the new test + the existing copy_graph tests with `cert: None`).

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/copy_graph.rs src/bin/copy_assign.rs
git commit -m "feat(vg-o2): per-read certificate tags (CP/PV/MP/ST) on audit W-lines"
```

---

### Task 5: Data-gated real-family byte-identity validation

**Files:**
- (No source changes.) A one-off validation run + a recorded result note.

**Interfaces:**
- Consumes: the built `copy_assign` binary and a real family region. This is a belt-and-suspenders end-to-end check on top of the Task-2/3 golden test (which already guarantees bit-identity at the `read_copy_evidence` level).

Note on the invocation: use the smallest known-family regression, **GSTM**. Its exact `copy_assign` command line is recorded in the run log/script beside the existing outputs `/home/juanfra/winloci_scratch/gstm_vg.*` and `b_GSTM_*` (a GSTM BAM + region/families input + `--out <prefix>`); confirm the flags with `./target/release/copy_assign --help`. Substitute that command for `<gstm-invocation>` below. The point is only that Step 1 and Step 2 run the **identical** command against two builds.

- [ ] **Step 1: Capture the PRE-change assignments (build at the branch's merge-base)**

Foreground, serial, output to `winloci_scratch`:
```bash
git stash 2>/dev/null; git checkout $(git merge-base HEAD dna-family-fallback) -- src/ 2>/dev/null || true
cargo build --release --bin copy_assign
PYTHONHASHSEED=0 ./target/release/copy_assign <gstm-invocation> --out /home/juanfra/winloci_scratch/vgo2_base
git checkout HEAD -- src/
md5sum /home/juanfra/winloci_scratch/vgo2_base.assignments.tsv
```
(Simpler alternative if a pre-change GSTM `.assignments.tsv` already exists on disk, e.g. `gstm_vg.assignments.tsv`: use it as the baseline and skip the merge-base rebuild.)

- [ ] **Step 2: Run the POST-change build on the same family**

Run (this branch's HEAD, foreground, serial, `winloci_scratch`, the SAME `<gstm-invocation>`):
```bash
cargo build --release --bin copy_assign
PYTHONHASHSEED=0 ./target/release/copy_assign <gstm-invocation> --out /home/juanfra/winloci_scratch/vgo2_new
md5sum /home/juanfra/winloci_scratch/vgo2_new.assignments.tsv
```

- [ ] **Step 3: Assert byte-identity**

Run: `diff /home/juanfra/winloci_scratch/vgo2_base.assignments.tsv /home/juanfra/winloci_scratch/vgo2_new.assignments.tsv && echo IDENTICAL`
Expected: `IDENTICAL` (zero diff; md5s equal). If not identical, the refactor drifted — STOP and fix (the golden test in Task 2/3 should have caught it; if the diff is real, do not "adjust" thresholds to hide it — report it).

- [ ] **Step 4: Optionally emit + eyeball the audit GFA**

Run the same family with the audit emit (`--phase` or `--vg-audit`), confirm the read `W`-lines carry `CP/PV/MP/ST` tags and every `P`/`W` step is `L`-line-backed (the existing copy_graph invariant test already asserts backing; here just eyeball one family).

- [ ] **Step 5: Commit the validation note**

Record the md5s + IDENTICAL result in the branch's SDD ledger / a short `bench/` note; no source change. (This task's "commit" is the recorded evidence, not code.)

---

## Notes for the implementer

- **Bit-identity is the whole point.** The Task-2 golden test freezes the exact `f64` output; Task 3 must reproduce it through the graph. Same columns, same order, same `lp_match`/`lp_mis` terms → same sum. If the golden test goes red, the refactor is wrong — do not re-capture the golden to make it pass.
- `read_copy_evidence` is shared by the hard gate (`assign_read_editing`) AND the EM (`em_copy_assign.rs`); both must stay bit-identical (the golden test + the EM suite cover this).
- The audit emit is opt-in/additive: a `ReadWalk` with `cert: None` emits exactly the old W-line, so the default `--phase` output and all existing copy_graph tests stay unchanged (once `cert: None` is added to their constructors).
- **WSL2 crash rule:** the Task-5 validation runs `copy_assign` FOREGROUND, serial, small family, output under `winloci_scratch` — never backgrounded/nohup/pkill.
```
