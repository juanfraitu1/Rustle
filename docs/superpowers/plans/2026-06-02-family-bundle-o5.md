# O5 Family-Guided Richer Bundle Implementation Plan

> ## ⛔ BLOCKED — DO NOT EXECUTE (premise falsified 2026-06-02)
> Grounding against the code AFTER this plan was written showed the structure-borrowing it builds (Tasks 3–7) **already exists and is default-on** in `--vg` (`vg::build_bundle_completion_nodes` / `build_bundle_borrow_junctions` / `build_bundle_borrow_coverage`, consumed at pipeline.rs:13089/13122/13194). Re-running `--vg --genome-fasta <chrY>` with that machinery on, RBMY c6 STILL yields 0 transcripts — so the real gap is on the **coverage/flow** side, not structure. See the CORRECTION + `## Re-grounded next step` in `docs/superpowers/specs/2026-06-02-family-bundle-o5-design.md`.
> **Still valid if/when a coverage-side fix lands:** Task 1 (`family_guided` field) + Task 2 (GTF emission) — the honesty labelling. Tasks 3–8 target the wrong layer and must be re-scoped after the diagnostic.

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Let a read-starved/5′-truncated paralog copy assemble its full-length isoform by borrowing the family's exon structure (O5: share evidence across copies via the family graph), while keeping abundance honest (copy's own reads only) and labelling the borrowed structure as inferred.

**Architecture:** Structure-borrow, abundance-own. When a family copy's normal per-bundle assembly yields nothing (or only a fragment) but it has ≥1 of its own reads, build a **completion bundle** = the copy's own reads + the copy's homologous exon scaffold projected from the `FamilyGraph` (`paralog_exon_spans`, already used by O2 rescue). Assemble it via the existing synthesis path (`synthesize_bundles_refined`). Because **no sibling reads are pooled**, abundance stays at the copy's own read level and only the copy's own isoform emits — structurally avoiding the naive-pool over-enumeration. Emitted transcripts are tagged `family_guided "true"` with a low `capacity_confidence`, so inferred structure is never reported as observed.

**Tech Stack:** Rust (rustle), `FamilyGraph` (`src/rustle/vg_hmm/family_graph.rs`), rescue synthesis (`src/rustle/vg_hmm/rescue.rs`), `--vg` pipeline (`src/rustle/pipeline.rs`), cargo test.

**Spec:** `docs/superpowers/specs/2026-06-02-family-bundle-o5-design.md` (draft 2026-06-02).

**Depends on:** the flow-capacity apportionment EM (same branch `vg/flow-capacity-apportionment`). The EM is what assigns each read to its copy; O5 consumes that assignment to know which reads are "the copy's own." O5 refuses to run if `RUSTLE_VG_ANCHOR_PRIOR=0`.

**Task dependency order:** T1 `family_guided` field (foundation) → T2 GTF emission → T3 `family_completion_spans` (pure) → T4 `is_starved_copy` (pure) → T5 flag plumbing → T6 `complete_starved_copies` synthesis wrapper → T7 pipeline wiring → T8 validation. Implement in order; T1 must compile before anything else.

---

### Task 1: Add `family_guided` field to `Transcript`

**Files:**
- Modify: `src/rustle/path_extract.rs:620` (the `Transcript` struct, add field after `abundance_min`)
- Modify: 37 `Transcript` literal sites (compiler-enforced): `src/rustle/path_extract.rs` (10), `src/rustle/pipeline.rs` (22), `src/rustle/transcript_filter.rs` (2), `src/rustle/merge_mode.rs` (1), `src/rustle/single_exon_pileup.rs` (1), `src/rustle/vg.rs` (1)

Notes (verified by reading code):
- `Transcript` has **no** `Default` impl and **no** `..spread` literals — every construction lists all fields explicitly, so the compiler enumerates every missing-field site after the struct changes. This is the same pattern the flow-capacity plan used for `capacity_confidence`.
- Every literal already contains a `capacity_confidence:` line; the new field defaults to `None` at every construction site (it is only set to `Some(true)` in Task 7).

- [ ] **Step 1: Add the field to the struct (after the `abundance_min` field at path_extract.rs:~697)**

```rust
    /// O5 family-guided assembly (spec 2026-06-02): `Some(true)` when this
    /// transcript's structure was completed using exon/junction scaffold
    /// borrowed from sibling paralog copies via the family graph (the copy's
    /// own reads did not span the full structure). Its EXPRESSION is observed
    /// (own-read coverage); its full-length STRUCTURE is inferred from family
    /// homology — `capacity_confidence` is correspondingly low. Emitted as the
    /// `family_guided` GTF attribute when `Some(true)`. `None`/`Some(false)`
    /// for ordinary transcripts. Only ever set in `--vg` mode.
    pub family_guided: Option<bool>,
```

- [ ] **Step 2: Build to get the exhaustive list of literals to fix**

Run: `cargo build 2>&1 | grep -E "missing field|E0063" | head -50`
Expected: ~37 `error[E0063]: missing field family_guided in initializer of ...::Transcript` lines, each with a file:line.

- [ ] **Step 3: At each reported site, insert `family_guided: None,`**

Insert immediately after the existing `abundance_min: ...,` (or `family_verdict: ...,`) entry in each literal. Every site defaults to `None`. Example (path_extract.rs:9190 area):
```rust
                    capacity_confidence: cap_conf_opt, abundance_min: abund_min_opt, family_guided: None, family_verdict: None, intron_low: Vec::new(), synthetic: false, rescue_class: None,
```

- [ ] **Step 4: Build clean**

Run: `cargo build 2>&1 | grep -E "error" | head`
Expected: no output (compiles).

- [ ] **Step 5: Commit**

```bash
git add src/rustle/path_extract.rs src/rustle/pipeline.rs src/rustle/transcript_filter.rs src/rustle/merge_mode.rs src/rustle/single_exon_pileup.rs src/rustle/vg.rs
git commit -m "feat(vg): add Transcript.family_guided field (O5 scaffold)"
```

---

### Task 2: Emit `family_guided` in GTF

**Files:**
- Modify: `src/rustle/gtf.rs:197` (after the `capacity_confidence` emission block, before `abundance_min`)
- Test: `src/rustle/gtf.rs` (test module)

- [ ] **Step 1: Write the failing test (add to the gtf.rs `#[cfg(test)] mod tests`)**

```rust
    #[test]
    fn family_guided_attr_emitted_when_some_true() {
        let mut tx = make_test_transcript(); // existing helper used by capacity_confidence tests
        tx.vg_copy_id = Some(0);
        tx.family_guided = Some(true);
        let line = format_transcript_attrs(&tx, 0.05); // existing attr-formatter used by the cc tests
        assert!(line.contains("family_guided \"true\";"), "got: {line}");
    }

    #[test]
    fn family_guided_attr_absent_when_none_or_false() {
        let mut tx = make_test_transcript();
        tx.vg_copy_id = Some(0);
        tx.family_guided = None;
        assert!(!format_transcript_attrs(&tx, 0.05).contains("family_guided"));
        tx.family_guided = Some(false);
        assert!(!format_transcript_attrs(&tx, 0.05).contains("family_guided"));
    }
```

(If the existing cc tests at gtf.rs:307 use an inline `tx` rather than helpers `make_test_transcript`/`format_transcript_attrs`, mirror their exact construction/call instead — match the surrounding test style.)

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test -p rustle --lib gtf::tests::family_guided -- --nocapture`
Expected: FAIL — attribute not emitted.

- [ ] **Step 3: Add the emission block (immediately after the `capacity_confidence` block, gtf.rs:~202)**

```rust
        // O5 family-guided structure flag (--vg only; spec 2026-06-02). Marks a
        // transcript whose full-length structure was borrowed from sibling
        // copies (expression observed, structure inferred). Only set in VG mode,
        // so non-VG GTF is byte-identical.
        if tx.family_guided == Some(true) {
            tx_attrs.push_str(" family_guided \"true\";");
        }
```

- [ ] **Step 4: Run to verify it passes**

Run: `cargo test -p rustle --lib gtf::tests::family_guided`
Expected: PASS (both tests).

- [ ] **Step 5: Commit**

```bash
git add src/rustle/gtf.rs
git commit -m "feat(vg): emit family_guided GTF attribute"
```

---

### Task 3: Pure helper `family_completion_spans`

Decides whether a starved copy's own exons should be completed to the family scaffold, and returns the completed structure. Pure (no `FamilyGraph`/`Bundle`) so it is trivially unit-testable; the caller feeds it the already-projected scaffold from `FamilyGraph::paralog_exon_spans`.

**Files:**
- Create: `src/rustle/vg_hmm/family_complete.rs`
- Modify: `src/rustle/vg_hmm/mod.rs` (add `pub mod family_complete;`)
- Test: `src/rustle/vg_hmm/family_complete.rs` (inline test module)

- [ ] **Step 1: Create the module with the function and failing tests**

```rust
//! O5 family-guided completion (spec 2026-06-02). Given a starved paralog
//! copy's own observed exons and the family scaffold projected into that copy's
//! coordinates (`FamilyGraph::paralog_exon_spans`), decide whether to complete
//! the copy's structure and return the completed exon list. NEVER fabricates a
//! copy with zero observed exons.

/// Return the family scaffold as the completed structure for this copy, or
/// `None` when completion is not warranted.
///
/// `own_exons`  — exons the copy's OWN reads actually cover (genomic-ascending).
/// `scaffold`   — the family's homologous exon spans in this copy's coordinates.
///
/// Returns `Some(scaffold)` iff:
///   - the copy has ≥1 own exon (no-fabrication guard: empty own → `None`), AND
///   - the scaffold strictly extends the copy's own footprint (more exons OR a
///     wider genomic span than `own_exons`), i.e. the own reads are truncated.
/// Otherwise `None` (own reads already span the structure → no help needed, or
/// nothing to borrow).
pub fn family_completion_spans(
    own_exons: &[(u64, u64)],
    scaffold: &[(u64, u64)],
) -> Option<Vec<(u64, u64)>> {
    if own_exons.is_empty() || scaffold.is_empty() {
        return None; // no observed expression, or nothing to borrow
    }
    let own_span = (own_exons.first().unwrap().0, own_exons.last().unwrap().1);
    let scaf_span = (scaffold.first().unwrap().0, scaffold.last().unwrap().1);
    let scaffold_extends = scaffold.len() > own_exons.len()
        || scaf_span.0 < own_span.0
        || scaf_span.1 > own_span.1;
    if scaffold_extends {
        Some(scaffold.to_vec())
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn completes_truncated_copy_to_scaffold() {
        // own reads cover only the 3' two exons; scaffold has the full 4.
        let own = [(2000, 2100), (2200, 2300)];
        let scaf = [(1000, 1100), (1500, 1600), (2000, 2100), (2200, 2300)];
        assert_eq!(family_completion_spans(&own, &scaf), Some(scaf.to_vec()));
    }

    #[test]
    fn no_completion_when_own_already_full() {
        let own = [(1000, 1100), (1500, 1600), (2000, 2100), (2200, 2300)];
        let scaf = own.clone();
        assert_eq!(family_completion_spans(&own, &scaf), None);
    }

    #[test]
    fn never_fabricates_silent_copy() {
        // zero observed exons -> never invent structure (DAZ3 trap).
        let scaf = [(1000, 1100), (1500, 1600)];
        assert_eq!(family_completion_spans(&[], &scaf), None);
    }

    #[test]
    fn no_completion_without_scaffold() {
        let own = [(2000, 2100)];
        assert_eq!(family_completion_spans(&own, &[]), None);
    }
}
```

- [ ] **Step 2: Wire the module in**

Add to `src/rustle/vg_hmm/mod.rs`:
```rust
pub mod family_complete;
```

- [ ] **Step 3: Run the tests**

Run: `cargo test -p rustle --lib vg_hmm::family_complete`
Expected: PASS (4 tests).

- [ ] **Step 4: Commit**

```bash
git add src/rustle/vg_hmm/family_complete.rs src/rustle/vg_hmm/mod.rs
git commit -m "feat(vg): family_completion_spans helper (O5, no-fabrication guard)"
```

---

### Task 4: Pure helper `is_starved_copy`

**Files:**
- Modify: `src/rustle/vg_hmm/family_complete.rs`
- Test: same file

- [ ] **Step 1: Add the function + failing tests**

```rust
/// Eligibility test for O5 completion. A copy is "starved" (assembly likely
/// truncated/empty despite real expression) when it has at least one own read
/// but its normal assembly covered less than `minfrac` of the family-consensus
/// length. `n_own_reads == 0` is NOT eligible (no observed expression).
///
/// `assembled_len` — total exonic length the copy's own assembly produced (0 if
///                   no transcript). `consensus_len` — family scaffold exonic
///                   length. `minfrac` — env `RUSTLE_VG_FAMILY_BUNDLE_MINFRAC`.
pub fn is_starved_copy(
    n_own_reads: usize,
    assembled_len: u64,
    consensus_len: u64,
    minfrac: f64,
) -> bool {
    if n_own_reads == 0 || consensus_len == 0 {
        return false;
    }
    (assembled_len as f64) < minfrac * (consensus_len as f64)
}
```

```rust
    #[test]
    fn starved_when_zero_assembly_with_reads() {
        assert!(is_starved_copy(2, 0, 14000, 0.6)); // c6: 2 reads, no transcript
    }

    #[test]
    fn starved_when_fragmentary() {
        assert!(is_starved_copy(2, 4000, 14000, 0.6)); // 4000/14000 = 0.29 < 0.6
    }

    #[test]
    fn not_starved_when_full() {
        assert!(!is_starved_copy(40, 13000, 14000, 0.6)); // 0.93 >= 0.6
    }

    #[test]
    fn not_starved_with_zero_reads() {
        assert!(!is_starved_copy(0, 0, 14000, 0.6)); // silent copy: never eligible
    }
```

- [ ] **Step 2: Run**

Run: `cargo test -p rustle --lib vg_hmm::family_complete::tests::starved`
Expected: PASS.

- [ ] **Step 3: Run the whole module**

Run: `cargo test -p rustle --lib vg_hmm::family_complete`
Expected: PASS (8 tests total).

- [ ] **Step 4: Commit**

```bash
git add src/rustle/vg_hmm/family_complete.rs
git commit -m "feat(vg): is_starved_copy eligibility helper (O5)"
```

---

### Task 5: Flag plumbing — `RUSTLE_VG_FAMILY_BUNDLE` + `_MINFRAC`

**Files:**
- Modify: `src/rustle/vg_hmm/family_complete.rs` (config reader)
- Test: same file

- [ ] **Step 1: Add a config reader + failing test**

```rust
/// O5 runtime config, read from env. `enabled` defaults OFF (opt-in prototype);
/// promote to default-ON-in-VG only after the spec's validation gates pass.
#[derive(Debug, Clone, Copy)]
pub struct O5Config {
    pub enabled: bool,
    pub minfrac: f64,
}

impl O5Config {
    pub fn from_env() -> Self {
        let enabled = std::env::var("RUSTLE_VG_FAMILY_BUNDLE")
            .map(|v| v == "1" || v.eq_ignore_ascii_case("true"))
            .unwrap_or(false);
        let minfrac = std::env::var("RUSTLE_VG_FAMILY_BUNDLE_MINFRAC")
            .ok()
            .and_then(|v| v.parse().ok())
            .filter(|f: &f64| *f > 0.0 && *f <= 1.0)
            .unwrap_or(0.6);
        O5Config { enabled, minfrac }
    }
}
```

```rust
    #[test]
    fn o5config_defaults_off_minfrac_point_six() {
        // env unset in this test process by default.
        let c = O5Config::from_env();
        assert!(!c.enabled);
        assert!((c.minfrac - 0.6).abs() < 1e-9);
    }
```

- [ ] **Step 2: Run**

Run: `cargo test -p rustle --lib vg_hmm::family_complete::tests::o5config`
Expected: PASS (assuming the env vars are unset; do not set them in this test).

- [ ] **Step 3: Commit**

```bash
git add src/rustle/vg_hmm/family_complete.rs
git commit -m "feat(vg): O5Config env reader (RUSTLE_VG_FAMILY_BUNDLE, default OFF)"
```

---

### Task 6: Synthesis wrapper `complete_starved_copies`

Bridges the helpers to the existing rescue synthesis. For each starved copy, project its scaffold (`FamilyGraph::paralog_exon_spans`), decide completion (`family_completion_spans`), and emit a completion `Bundle` from the copy's own read node-paths via `synthesize_bundles_refined` (with `min_reads=1`, so a 2-read copy is not dropped). Tag the bundle so its transcripts can be marked `family_guided` downstream.

**Files:**
- Modify: `src/rustle/vg_hmm/family_complete.rs`
- Reference (read, do not modify): `src/rustle/vg_hmm/rescue.rs:1396` (`synthesize_bundles_refined`), `family_graph.rs:118` (`paralog_exon_spans`)
- Test: same file (focused routing test; full bundle synthesis is validated by the RBMY integration test in Task 8)

- [ ] **Step 1: Add the wrapper**

```rust
use crate::types::Bundle;
use crate::vg_hmm::family_graph::{CopyId, FamilyGraph};

/// One starved copy's input: the copy id within `fg`, the number of the copy's
/// OWN reads, the exonic length its normal assembly produced, and the copy's
/// own read node-paths (same tuple shape `synthesize_bundles_refined` consumes:
/// (read_name, node_path, read_spans)).
pub struct StarvedCopy {
    pub cid: CopyId,
    pub n_own_reads: usize,
    pub assembled_len: u64,
    pub read_paths: Vec<(String, Vec<crate::vg_hmm::family_graph::NodeIdx>, Vec<(usize, usize)>)>,
}

/// Build completion bundles for the eligible starved copies of one family.
/// Returns bundles already tagged with `vg_family_id` and
/// `rescue_class = Some(FamilyCompleted)` (see Task 7 for the enum variant).
/// Empty Vec when O5 is disabled or no copy is eligible.
pub fn complete_starved_copies(
    fg: &FamilyGraph,
    family_id: usize,
    starved: &[StarvedCopy],
    cfg: O5Config,
) -> Vec<Bundle> {
    if !cfg.enabled {
        return Vec::new();
    }
    let mut out = Vec::new();
    for sc in starved {
        let scaffold = fg.paralog_exon_spans(sc.cid);
        let consensus_len: u64 = scaffold.iter().map(|(s, e)| e.saturating_sub(*s)).sum();
        if !is_starved_copy(sc.n_own_reads, sc.assembled_len, consensus_len, cfg.minfrac) {
            continue;
        }
        // own_exons proxy: the scaffold nodes the copy's own reads actually hit.
        let own_exons = own_exons_from_read_paths(fg, sc.cid, &sc.read_paths);
        if family_completion_spans(&own_exons, &scaffold).is_none() {
            continue; // own reads already span structure, or nothing to borrow
        }
        // Reuse the proven synthesis primitive with min_reads=1 (1 long read is
        // informative). It clusters the copy's own read node-paths against the
        // family graph and emits a structurally-complete bundle.
        let mut bundles = crate::vg_hmm::rescue::synthesize_bundles_refined(
            fg, &sc.read_paths, /*min_reads=*/ 1, /*kmer_len=*/ 15,
            /*hard_min_reads=*/ 1, /*cluster_window=*/ 20,
        );
        for b in bundles.iter_mut() {
            b.vg_family_id = Some(family_id);
            b.rescue_class = Some(crate::vg_hmm::diagnostic::RescueClass::FamilyCompleted);
        }
        out.append(&mut bundles);
    }
    out
}

/// Exons of the copy's OWN reads, taken as the scaffold spans of the nodes its
/// reads traverse (genomic-ascending, deduped).
fn own_exons_from_read_paths(
    fg: &FamilyGraph,
    cid: CopyId,
    read_paths: &[(String, Vec<crate::vg_hmm::family_graph::NodeIdx>, Vec<(usize, usize)>)],
) -> Vec<(u64, u64)> {
    use std::collections::BTreeSet;
    let mut hit: BTreeSet<(u64, u64)> = BTreeSet::new();
    for (_, path, _) in read_paths {
        for nidx in path {
            if let Some((_, sp)) = fg.nodes[nidx.0]
                .per_copy_spans
                .iter()
                .find(|(c, _)| *c == cid)
            {
                hit.insert(*sp);
            }
        }
    }
    hit.into_iter().collect()
}
```

- [ ] **Step 2: Add a focused test (disabled-config short-circuit; the no-fabrication path)**

```rust
    #[test]
    fn disabled_config_yields_no_bundles() {
        // With O5 disabled, never synthesize regardless of inputs.
        let fg = crate::vg_hmm::family_graph::FamilyGraph { family_id: 0, nodes: vec![], edges: vec![] };
        let cfg = O5Config { enabled: false, minfrac: 0.6 };
        let starved = vec![StarvedCopy { cid: 0, n_own_reads: 2, assembled_len: 0, read_paths: vec![] }];
        assert!(complete_starved_copies(&fg, 0, &starved, cfg).is_empty());
    }

    #[test]
    fn empty_scaffold_yields_no_bundles() {
        // enabled, but the family graph has no nodes -> scaffold empty -> no fabrication.
        let fg = crate::vg_hmm::family_graph::FamilyGraph { family_id: 0, nodes: vec![], edges: vec![] };
        let cfg = O5Config { enabled: true, minfrac: 0.6 };
        let starved = vec![StarvedCopy { cid: 0, n_own_reads: 2, assembled_len: 0, read_paths: vec![] }];
        assert!(complete_starved_copies(&fg, 0, &starved, cfg).is_empty());
    }
```

- [ ] **Step 3: Add the `FamilyCompleted` enum variant (referenced above)**

In `src/rustle/vg_hmm/diagnostic.rs`, add to the `RescueClass` enum (find `pub enum RescueClass`):
```rust
    /// O5: bundle synthesized by completing a starved copy's own reads with the
    /// family scaffold (spec 2026-06-02). Transcripts carry `family_guided`.
    FamilyCompleted,
```
Then `cargo build 2>&1 | grep -E "non-exhaustive|E0004"` — if any `match` on `RescueClass` is now non-exhaustive, add a `FamilyCompleted => ...` arm mirroring the nearest existing rescue variant's arm.

- [ ] **Step 4: Run**

Run: `cargo test -p rustle --lib vg_hmm::family_complete && cargo build 2>&1 | grep -E "error" | head`
Expected: tests PASS; build clean.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_hmm/family_complete.rs src/rustle/vg_hmm/diagnostic.rs
git commit -m "feat(vg): complete_starved_copies synthesis wrapper + FamilyCompleted class"
```

---

### Task 7: Wire O5 into the VG pipeline

Run the completion pass after per-copy VG assembly: detect starved copies, synthesize completion bundles, assemble them through the existing synthetic-bundle path, and mark resulting transcripts `family_guided = Some(true)` with a low `capacity_confidence`.

**Files:**
- Modify: `src/rustle/pipeline.rs` near the family-assembly block (after `family_graphs` is built at :10702 and per-copy bundles have assembled; co-locate with the existing rescue call at :10485 or just after the VG family assembly completes)

- [ ] **Step 1: Read the surrounding code to choose the precise insertion point**

Run: `sed -n '10480,10520p' src/rustle/pipeline.rs && sed -n '10700,10860p' src/rustle/pipeline.rs`
Identify: (a) where `family_graphs` and the per-copy read→node-path mappings are in scope, (b) the `Vec<Bundle>` that feeds the normal assembly, (c) where per-copy assembled transcripts are available to measure `assembled_len`.

- [ ] **Step 2: Add the O5 pass (guarded by O5Config + the EM-dependency check)**

Insert after the VG family assembly produces its transcripts, before they merge into `all_transcripts`:
```rust
        // ── O5: family-guided completion of starved copies (spec 2026-06-02) ──
        // Opt-in (RUSTLE_VG_FAMILY_BUNDLE=1). Requires the anchor-prior EM so
        // read→copy assignment is trustworthy; refuse otherwise.
        let o5 = crate::vg_hmm::family_complete::O5Config::from_env();
        if o5.enabled {
            let anchor_prior_on = std::env::var("RUSTLE_VG_ANCHOR_PRIOR")
                .map(|v| v != "0").unwrap_or(true);
            if !anchor_prior_on {
                eprintln!("[O5] RUSTLE_VG_FAMILY_BUNDLE set but RUSTLE_VG_ANCHOR_PRIOR=0 — refusing (abundance would be untrusted)");
            } else {
                for (pi, fg_opt) in family_graphs.iter().enumerate() {
                    let fg = match fg_opt { Some(f) => f, None => continue };
                    // Build StarvedCopy inputs from this partition's per-copy
                    // read node-paths + measured assembled length (see Step 1
                    // for the exact in-scope variable names).
                    let starved = build_starved_inputs(fg, pi /*, …in-scope maps… */);
                    let completion = crate::vg_hmm::family_complete::complete_starved_copies(
                        fg, fg.family_id, &starved, o5,
                    );
                    for b in &completion {
                        eprintln!("[O5] family {} completion bundle: {} exon spans", fg.family_id, b.reads.len());
                    }
                    o5_bundles.extend(completion);
                }
            }
        }
        // Assemble O5 completion bundles through the normal synthetic path and
        // tag their transcripts.
        for b in &o5_bundles {
            let mut txs = extract_transcripts(/* same args as the synthetic-bundle assembly call */);
            for t in txs.iter_mut() {
                t.family_guided = Some(true);
                // Low capacity_confidence: borrowed structure is unanchored.
                if t.capacity_confidence.is_none() { t.capacity_confidence = Some(0.0); }
                t.vg_family_id = Some(fg_family_id_for(b));
            }
            all_transcripts.append(&mut txs);
        }
```

Notes for the implementer (resolve against the code read in Step 1):
- `build_starved_inputs` is a small local closure/fn assembling `Vec<StarvedCopy>` from the partition's existing read→node-path data and the per-copy assembled-length you measure from the just-assembled transcripts for that copy (sum of `exons` lengths; 0 if none).
- Reuse the **exact** `extract_transcripts(...)` argument list already used for synthetic-bundle assembly elsewhere in this block — do not invent a new call shape.
- `o5_bundles` is a `Vec<Bundle>` declared just before the loop.

- [ ] **Step 3: Build**

Run: `cargo build 2>&1 | grep -E "error" | head`
Expected: no errors (resolve in-scope names until clean).

- [ ] **Step 4: Smoke-run on RBMY (default OFF stays off)**

Run: `cargo build --release 2>&1 | tail -1`
Then with O5 OFF (default): `target/release/rustle -L /tmp/tspy.bam --vg --genome-fasta <genome> -o /tmp/rbmy/o5off.gtf 2>/dev/null; grep -c family_guided /tmp/rbmy/o5off.gtf`
Expected: `0` (no `family_guided` when the flag is off).

- [ ] **Step 5: Commit**

```bash
git add src/rustle/pipeline.rs
git commit -m "feat(vg): wire O5 family-guided completion pass (opt-in)"
```

---

### Task 8: Validation

**Files:**
- No source changes (validation only). Record results in `docs/superpowers/specs/2026-06-02-family-bundle-o5-design.md` under a new `## Validation results` section.

- [ ] **Step 1: Regression guard — default de-novo byte-identical**

Run (no `--vg`):
```bash
target/release/rustle -L bench/<GGO_19 bam> -G bench/GGO_19.gtf -o /tmp/o5_denovo.gtf
# compare to a HEAD-built baseline GTF
diff <(grep -v '^#' /tmp/o5_baseline.gtf) <(grep -v '^#' /tmp/o5_denovo.gtf) | head
```
Expected: empty diff (byte-identical); de-novo metrics unchanged at 95.6/90.5. **Hard gate — must pass.**

- [ ] **Step 2: RBMY c6 primary success (O5 ON)**

Run:
```bash
RUSTLE_VG_FAMILY_BUNDLE=1 RUSTLE_VG_ANCHOR_PRIOR=1 \
  target/release/rustle -L /tmp/tspy.bam --vg --genome-fasta <genome> -o /tmp/rbmy/o5on.gtf
# c6 region 19,717,578-19,730,926
awk '$4>=19717578 && $5<=19731200 && $3=="transcript"' /tmp/rbmy/o5on.gtf
```
Expected: c6 now emits **≥1** transcript spanning the gene, each tagged `family_guided "true"` with a low `capacity_confidence`. Record the transcript count and coverage.

- [ ] **Step 2: Over-enumeration guard (the flaw test)**

From the same `o5on.gtf`, assert c6's transcript count equals its **own-read-supported** isoform count (expected **1**, not the 6 the naive all-reads pool produced):
```bash
awk '$4>=19717578 && $5<=19731200 && $3=="transcript"' /tmp/rbmy/o5on.gtf | wc -l
```
Expected: `1` (or the small own-read count) — **not** 6. If >1, inspect: completion must use the copy's OWN read paths only, never sibling reads.

- [ ] **Step 3: c1–c5 unchanged**

Run: `diff <(grep -v family_guided /tmp/rbmy/o5on.gtf | awk '$4<19717578') <(awk '$4<19717578' /tmp/rbmy/o5off.gtf)`
Expected: empty diff — O5 is a no-op on the already-resolved copies c1–c5.

- [ ] **Step 4: No-fabrication guard (DAZ3 stays silent)**

Run O5 ON on the DAZ region; confirm DAZ3 (the inverted near-identical copy with ~2 decisive reads) does **not** gain a `family_guided` transcript:
```bash
RUSTLE_VG_FAMILY_BUNDLE=1 RUSTLE_VG_ANCHOR_PRIOR=1 \
  target/release/rustle -L /tmp/daz.bam --vg --genome-fasta <genome> -o /tmp/daz_o5.gtf
awk '$4>=42879918 && $5<=42945552' /tmp/daz_o5.gtf | grep -c family_guided
```
Expected: `0` — DAZ3 has no decisive own reads → not eligible (≤ noise), never fabricated.

- [ ] **Step 5: Genome-wide spillover guard**

Re-run `bench/paralog_secondary_scan` with O5 ON; confirm the 89 `pure_spillover` candidates gain **no** `family_guided` transcripts (they lack decisive own reads → ineligible).
Expected: 0 new guided transcripts on spillover loci.

- [ ] **Step 6: Record results + commit**

Append a `## Validation results (2026-06-XX)` section to the design doc with the numbers from Steps 1–5, then:
```bash
git add docs/superpowers/specs/2026-06-02-family-bundle-o5-design.md
git commit -m "docs(vg): record O5 family-bundle validation results"
```

---

## Promotion criteria (flag → permanent)

Flip `RUSTLE_VG_FAMILY_BUNDLE` to default-ON-in-VG only when **all** hold (from the spec):
1. Task 8 Step 1 byte-identical (no de-novo impact) — hard requirement.
2. RBMY c6 emits the correct single full-length isoform at its own-read abundance with honesty labels (Step 2 + over-enum guard).
3. No-fabrication + spillover guards show zero spurious guided transcripts (Steps 4–5).
4. On a real multi-copy benchmark, O5 adds recovered real copies without adding FPs (precision flat or up).

Until then it stays opt-in, documented in `docs/vg_genome_scoping.md`.
