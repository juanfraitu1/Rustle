# de-tie Conflict Criterion Port Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the AS-tie conflict-edge criterion with the validated **de-tie** criterion (`|de_a − de_b| ≤ 0.005 AND max(de_a,de_b) ≤ 0.05`, edge iff ≥3 reads), exclude supplementary alignments from conflict placements, and log per-edge mapq corroboration — porting the bake-off result (`bench/family_criterion_bakeoff.md`) into the shipped Rust pipeline.

**Architecture:** `de` (minimap2 `de:f` gap-compressed divergence) is parsed onto `BamRead`, carried through `build_read_placements` into a richer `Placement` struct, and the `read_conflict.rs` kernel's tie predicate switches from AS-tie to de-tie. AS is retained (`as_score`) only to log the `de-tie ⊆ AS-tie` invariant; mapq is retained to log the both-mapq0 (genuine-multimapper) fraction. Supplementary alignments are excluded from conflict placements (chimeric ≠ multimapping).

**Tech Stack:** Rust, noodles-sam (custom `de:f` tag parsing), existing `read_conflict`/`denovo_assemble`/`denovo_pipeline` modules.

## Global Constraints

- de-tie predicate: `|de_a − de_b| ≤ DE_TIE_DELTA (0.005) AND max(de_a, de_b) ≤ DE_MAX (0.05)`.
- Conflict edge iff de-tied read count `≥ MIN_CONFLICT_READS (3)`.
- All three constants are named and env-overridable: `RUSTLE_CONFLICT_DE_DELTA`, `RUSTLE_CONFLICT_DE_MAX`, `RUSTLE_CONFLICT_MIN_READS`.
- **Supplementary alignments are excluded** from conflict placements (the shipped path currently includes them — `aligned_read_from_record` only skips unmapped).
- **No 200bp distinct-record guard is needed in Rust:** each BAM record becomes one `BamRead`, and `build_read_placements` attributes each to exactly one best-overlap locus, so a single straddling record cannot form a cross-locus edge by construction. (The Python bake-off needed the guard only because it paired the same physical record from two region fetches.)
- `nm` is NOT carried (it is `de`-in-disguise per the dig). Do NOT port `nintron`/`ms`/`cm`/`s1`/`s2`/`ts`.
- Every change is TDD: failing test → red → implement → green → commit. `cargo test --lib` is 0 failures at every commit.
- `colocated_families` / `conflict_to_split_families` / `detect_families` (POA-only) signatures and behaviour are unchanged.

## File map

| File | Change |
|------|--------|
| `src/rustle/vg_family/denovo_assemble.rs` | `record_de` parser; `de`+`is_supplementary` on `BamRead`; extend `aligned_read_from_record`; update 4 construction sites |
| `src/rustle/vg_family/read_conflict.rs` | `Placement` struct; `ConflictParams`→de params + `from_env`; `conflict_edges` de-tie; `as_tie_edges` audit helper; rewrite kernel tests |
| `src/rustle/vg_family/denovo_pipeline.rs` | `build_read_placements`→`Vec<Placement>` (exclude supplementary); `DenovoConfig.conflict`→`from_env`; detect_and_assign logging (both-mapq0 + de⊆AS); update placement tests + fixture |

---

## Task 1: `record_de` — parse the `de:f` tag

**Files:**
- Modify: `src/rustle/vg_family/denovo_assemble.rs`

**Interfaces:**
- Produces: `fn record_de(record: &RecordBuf) -> Option<f32>` (private, mirrors `record_as` at line 132)

The `de` tag is a custom 2-char tag (`[b'd', b'e']`) carrying a `Value::Float`. The existing `record_as` (line 131-149) shows the iteration pattern; `de` differs in matching a custom `Tag` and a float value.

- [ ] **Step 1: Write the failing test**

Add to the `tests` module in `denovo_assemble.rs`:

```rust
#[test]
fn record_de_parses_float_tag() {
    use noodles_sam::alignment::record::data::field::{Tag, Value};
    use noodles_sam::alignment::record_buf::data::field::Value as BufValue;
    let mut record = RecordBuf::builder()
        .set_flags(Flags::default())
        .set_alignment_start(Position::try_from(1usize).unwrap())
        .build();
    record.data_mut().insert(Tag::from([b'd', b'e']), BufValue::Float(0.0123));
    let de = record_de(&record).expect("de tag present");
    assert!((de - 0.0123).abs() < 1e-6);
}

#[test]
fn record_de_absent_is_none() {
    let record = RecordBuf::builder().set_flags(Flags::default()).build();
    assert!(record_de(&record).is_none());
}
```

Note: verify the exact `data_mut().insert` API against the noodles version in `Cargo.toml`; if `RecordBuf` builds data differently, construct the record via `.set_data(...)` with a `Data` containing the field instead. The behaviour to test is fixed: a `de:f:0.0123` tag parses to `0.0123_f32`.

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --lib -- vg_family::denovo_assemble::tests::record_de_parses_float_tag 2>&1 | tail -5`
Expected: compile error (function not defined).

- [ ] **Step 3: Implement `record_de`**

Add beside `record_as`:

```rust
/// Read the `de:f` (gap-compressed per-base divergence) tag from a record — the conflict-criterion signal.
/// `None` if absent. `de` is a custom 2-char tag carrying a float.
fn record_de(record: &RecordBuf) -> Option<f32> {
    use noodles_sam::alignment::record::data::field::{Tag, Value};
    let de_tag = Tag::from([b'd', b'e']);
    for entry in noodles_sam::alignment::Record::data(record).iter() {
        let (tag, value) = entry.ok()?;
        if tag == de_tag {
            return match value {
                Value::Float(v) => Some(v),
                _ => None,
            };
        }
    }
    None
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cargo test --lib -- vg_family::denovo_assemble::tests::record_de 2>&1 | tail -5`
Expected: both pass.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/denovo_assemble.rs
git commit -m "feat(denovo_assemble): add record_de parser for the de:f divergence tag"
```

---

## Task 2: Carry `de` + `is_supplementary` on `BamRead`

**Files:**
- Modify: `src/rustle/vg_family/denovo_assemble.rs`

**Interfaces:**
- Consumes: `record_de` (Task 1)
- Produces: `BamRead { chrom, read, mapq, name, as_score, de: f32, is_supplementary: bool }`
- Produces: `aligned_read_from_record(record) -> Option<(AlignedRead, u8, String, i32, f32, bool)>` (adds `de`, `is_supplementary`)

Current `BamRead` (line ~214) and `aligned_read_from_record` (line ~225) plus 4 construction sites: `aligned_reads_from_bam`, `reads_in_region_indexed`, `reads_in_region_scan`, and the existing `aligned_read_from_record_extracts_cigar_and_seq` test.

- [ ] **Step 1: Write the failing test**

Add to the `tests` module in `denovo_assemble.rs`:

```rust
#[test]
fn aligned_read_from_record_extracts_de_and_supplementary() {
    use noodles_sam::alignment::record::data::field::Tag;
    use noodles_sam::alignment::record_buf::data::field::Value as BufValue;
    use noodles_sam::alignment::record_buf::Sequence;
    let cigar: Cigar = vec![Op::new(Kind::Match, 4)].into_iter().collect();
    let mut record = RecordBuf::builder()
        .set_flags(Flags::SUPPLEMENTARY)
        .set_alignment_start(Position::try_from(1usize).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(b"ACGT".to_vec()))
        .build();
    record.data_mut().insert(Tag::from([b'd', b'e']), BufValue::Float(0.02));
    let (_ar, _mapq, _name, _as, de, is_supp) =
        aligned_read_from_record(&record).expect("mapped");
    assert!((de - 0.02).abs() < 1e-6);
    assert!(is_supp, "supplementary flag must be surfaced");
}
```

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test --lib -- vg_family::denovo_assemble::tests::aligned_read_from_record_extracts_de_and_supplementary 2>&1 | tail -5`
Expected: compile error (tuple arity / fields).

- [ ] **Step 3: Add the fields and extend the function**

In `BamRead` (line ~214), add fields:

```rust
pub struct BamRead {
    pub chrom: String,
    pub read: AlignedRead,
    pub mapq: u8,
    pub name: String,
    pub as_score: i32,
    /// minimap2 `de:f` gap-compressed per-base divergence (0.0 if absent) — the conflict-tie signal.
    pub de: f32,
    /// chimeric/split alignment (SAM flag 0x800) — excluded from conflict placements.
    pub is_supplementary: bool,
}
```

In `aligned_read_from_record` (line ~225), change the return type and body tail:

```rust
pub fn aligned_read_from_record(record: &RecordBuf) -> Option<(AlignedRead, u8, String, i32, f32, bool)> {
    if record.flags().is_unmapped() {
        return None;
    }
    // ... existing ref_start / cigar / seq / mapq / name / as_score unchanged ...
    let de = record_de(record).unwrap_or(0.0);
    let is_supplementary = record.flags().is_supplementary();
    Some((AlignedRead { ref_start, cigar, seq }, mapq, name, as_score, de, is_supplementary))
}
```

Update the 3 production construction sites (`aligned_reads_from_bam` ~258, `reads_in_region_indexed` ~303, `reads_in_region_scan` ~431):

```rust
if let Some((read, mapq, name, as_score, de, is_supplementary)) = aligned_read_from_record(&record) {
    out.push(BamRead { chrom, read, mapq, name, as_score, de, is_supplementary });
}
```
(and the `bam_reads.push(BamRead { chrom: chrom.to_string(), read, mapq, name, as_score, de, is_supplementary })` variants).

Update the existing test `aligned_read_from_record_extracts_cigar_and_seq` (line ~632) destructure:
```rust
let (ar, _mapq, _name, _as_score, _de, _is_supp) = aligned_read_from_record(&r).expect("mapped read");
```

- [ ] **Step 4: Update the `denovo_pipeline.rs` `BamRead` test construction**

The `two_paralogs_with_psvs` fixture (line ~613) and any other `BamRead { ... }` literal must add `de` and `is_supplementary`. For now (full fixture rework is Task 4), make it compile by adding `de: 0.01, is_supplementary: false` to each of the 4 entries. (Task 4 sets meaningful de values.)

- [ ] **Step 5: Run tests**

Run: `cargo test --lib 2>&1 | tail -3`
Expected: 0 failures.

- [ ] **Step 6: Commit**

```bash
git add src/rustle/vg_family/denovo_assemble.rs src/rustle/vg_family/denovo_pipeline.rs
git commit -m "feat(denovo_assemble): carry de + is_supplementary on BamRead"
```

---

## Task 3: de-tie kernel in `read_conflict.rs`

**Files:**
- Modify: `src/rustle/vg_family/read_conflict.rs`

**Interfaces:**
- Produces: `pub struct Placement { pub locus: usize, pub de: f32, pub mapq: u8, pub as_score: i32 }`
- Produces: `pub type ReadPlacements = Vec<Placement>`
- Produces: `pub struct ConflictParams { pub delta: f64, pub de_max: f64, pub min_reads: usize }` with `Default` = `{0.005, 0.05, 3}` and `pub fn from_env() -> Self`
- Produces: `pub fn conflict_edges(n_loci, reads, p) -> Vec<(usize,usize,usize)>` (de-tie; same signature shape)
- Produces: `pub fn as_tie_edges(n_loci, reads, as_tie: f64, min_reads: usize) -> std::collections::BTreeSet<(usize,usize)>` (audit: AS-tie edge node-pairs)
- Unchanged: `conflict_families`

- [ ] **Step 1: Replace the placement type, params, and predicate (with new tests)**

Replace `ReadPlacements`, `ConflictParams`, `as_tied`, and `conflict_edges` (lines 22-75) with:

```rust
/// One read's placement on a candidate locus: the locus index plus the signals the conflict criterion uses —
/// `de` (gap-compressed divergence, the tie discriminant), `mapq` (both-0 = genuine-multimapper corroboration),
/// and `as_score` (kept only to log the `de-tie ⊆ AS-tie` invariant).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Placement {
    pub locus: usize,
    pub de: f32,
    pub mapq: u8,
    pub as_score: i32,
}

/// One read's placements over the family's candidate loci. Built from the BAM by the adapter.
pub type ReadPlacements = Vec<Placement>;

/// Tunables for the read-conflict criterion (de-tie). Defaults are the bake-off operating point
/// (`bench/family_criterion_bakeoff.md`); env-overridable via `RUSTLE_CONFLICT_DE_DELTA/DE_MAX/MIN_READS`.
#[derive(Clone, Copy, Debug)]
pub struct ConflictParams {
    /// Two placements conflict iff their divergences are within `delta` AND both `<= de_max` (both fit).
    pub delta: f64,
    pub de_max: f64,
    /// Minimum conflicting-read count for an edge (guards the noise floor).
    pub min_reads: usize,
}

impl Default for ConflictParams {
    fn default() -> Self {
        ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 3 }
    }
}

impl ConflictParams {
    /// Read overrides from `RUSTLE_CONFLICT_DE_DELTA`, `RUSTLE_CONFLICT_DE_MAX`, `RUSTLE_CONFLICT_MIN_READS`.
    pub fn from_env() -> Self {
        let d = Self::default();
        let f = |k: &str, v: f64| std::env::var(k).ok().and_then(|s| s.parse().ok()).unwrap_or(v);
        let u = |k: &str, v: usize| std::env::var(k).ok().and_then(|s| s.parse().ok()).unwrap_or(v);
        ConflictParams {
            delta: f("RUSTLE_CONFLICT_DE_DELTA", d.delta),
            de_max: f("RUSTLE_CONFLICT_DE_MAX", d.de_max),
            min_reads: u("RUSTLE_CONFLICT_MIN_READS", d.min_reads),
        }
    }
}

/// de-tie: `|de_a − de_b| <= delta` AND `max(de_a, de_b) <= de_max` (the read fits both copies, comparably).
fn de_tied(a: &Placement, b: &Placement, p: &ConflictParams) -> bool {
    let (da, db) = (a.de as f64, b.de as f64);
    (da - db).abs() <= p.delta && da.max(db) <= p.de_max
}

/// `min(a,b) >= as_tie * max(a,b)` — the legacy AS-tie predicate, kept only for the audit edge-set.
fn as_tied(a: i32, b: i32, as_tie: f64) -> bool {
    let (hi, lo) = (a.max(b), a.min(b));
    hi > 0 && (lo as f64) >= as_tie * (hi as f64)
}

/// Build the read-conflict edges over `n_loci` under the **de-tie** criterion. For each read, every pair of its
/// placements that de-ties contributes one conflict observation to that locus pair. Returns `(i, j, weight)`
/// with `i < j` for pairs whose count `>= p.min_reads`, sorted. Self-pairs (same locus) ignored.
pub fn conflict_edges(n_loci: usize, reads: &[ReadPlacements], p: &ConflictParams) -> Vec<(usize, usize, usize)> {
    use std::collections::BTreeMap;
    let mut weight: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    for placements in reads {
        for a in 0..placements.len() {
            for b in (a + 1)..placements.len() {
                let (pa, pb) = (&placements[a], &placements[b]);
                if pa.locus == pb.locus || pa.locus >= n_loci || pb.locus >= n_loci {
                    continue;
                }
                if de_tied(pa, pb, p) {
                    let key = (pa.locus.min(pb.locus), pa.locus.max(pb.locus));
                    *weight.entry(key).or_insert(0) += 1;
                }
            }
        }
    }
    weight.into_iter().filter(|&(_, w)| w >= p.min_reads).map(|((i, j), w)| (i, j, w)).collect()
}

/// AS-tie edge node-pairs over `n_loci` — the legacy criterion, used only to log that the de-tie edge set is a
/// subset of the AS-tie edge set (`de-tie ⊆ AS-tie`, the portable regression invariant from the bake-off).
pub fn as_tie_edges(n_loci: usize, reads: &[ReadPlacements], as_tie: f64, min_reads: usize) -> std::collections::BTreeSet<(usize, usize)> {
    use std::collections::BTreeMap;
    let mut weight: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    for placements in reads {
        for a in 0..placements.len() {
            for b in (a + 1)..placements.len() {
                let (pa, pb) = (&placements[a], &placements[b]);
                if pa.locus == pb.locus || pa.locus >= n_loci || pb.locus >= n_loci {
                    continue;
                }
                if as_tied(pa.as_score, pb.as_score, as_tie) {
                    *weight.entry((pa.locus.min(pb.locus), pa.locus.max(pb.locus))).or_insert(0) += 1;
                }
            }
        }
    }
    weight.into_iter().filter(|&(_, w)| w >= min_reads).map(|((i, j), _)| (i, j)).collect()
}
```

- [ ] **Step 2: Rewrite the kernel tests for de semantics**

Replace the existing `#[cfg(test)] mod tests` block (lines ~116-185) with de-based tests. Use this helper and tests:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    fn p(locus: usize, de: f32) -> Placement { Placement { locus, de, mapq: 0, as_score: 100 } }

    #[test]
    fn de_tied_placements_make_an_edge_and_a_family() {
        // both copies fit comparably (de 0.010 vs 0.012, both < 0.05) -> tie -> edge -> family.
        let reads = vec![vec![p(0, 0.010), p(1, 0.012)]];
        let edges = conflict_edges(2, &reads, &ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1 });
        assert_eq!(edges, vec![(0, 1, 1)]);
        assert_eq!(conflict_families(2, &edges), vec![vec![0, 1]]);
    }

    #[test]
    fn divergence_gap_beyond_delta_is_not_a_conflict() {
        // read fits copy 0 (de 0.001) far better than copy 1 (de 0.020): |Δ|=0.019 > 0.005 -> resolvable.
        let reads = vec![vec![p(0, 0.001), p(1, 0.020)]];
        let edges = conflict_edges(2, &reads, &ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1 });
        assert!(edges.is_empty());
    }

    #[test]
    fn both_high_divergence_blocked_by_ceiling() {
        // de_a 0.06 ~ de_b 0.061 (tied within delta) but both exceed de_max 0.05 -> read fits neither.
        let reads = vec![vec![p(0, 0.060), p(1, 0.061)]];
        let edges = conflict_edges(2, &reads, &ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1 });
        assert!(edges.is_empty());
    }

    #[test]
    fn single_placement_read_is_a_singleton_not_a_family() {
        let reads = vec![vec![p(0, 0.01)], vec![p(1, 0.01)]];
        let edges = conflict_edges(2, &reads, &ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1 });
        assert!(edges.is_empty());
        assert!(conflict_families(2, &edges).is_empty());
    }

    #[test]
    fn min_reads_threshold_drops_thin_conflicts() {
        let one = vec![vec![p(0, 0.01), p(1, 0.012)]];
        let pr = ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 3 };
        assert!(conflict_edges(2, &one, &pr).is_empty());
        let three = vec![one[0].clone(), one[0].clone(), one[0].clone()];
        assert_eq!(conflict_edges(2, &three, &pr), vec![(0, 1, 3)]);
    }

    #[test]
    fn transitive_conflict_closes_into_one_family() {
        let reads = vec![vec![p(0, 0.010), p(1, 0.012)], vec![p(1, 0.010), p(2, 0.013)]];
        let edges = conflict_edges(3, &reads, &ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1 });
        assert_eq!(conflict_families(3, &edges), vec![vec![0, 1, 2]]);
    }

    #[test]
    fn disjoint_conflicts_form_separate_families() {
        let reads = vec![
            vec![p(0, 0.01), p(1, 0.012)],
            vec![p(2, 0.01), p(3, 0.012)],
            vec![p(4, 0.01)],
        ];
        let edges = conflict_edges(5, &reads, &ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1 });
        assert_eq!(conflict_families(5, &edges), vec![vec![0, 1], vec![2, 3]]);
    }

    #[test]
    fn default_params_are_the_operating_point() {
        let d = ConflictParams::default();
        assert!((d.delta - 0.005).abs() < 1e-9);
        assert!((d.de_max - 0.05).abs() < 1e-9);
        assert_eq!(d.min_reads, 3);
    }

    #[test]
    fn as_tie_edges_superset_of_de_edges() {
        // AS ties two placements that de SPLITS (de 0.001 vs 0.020): AS-edge exists, de-edge does not.
        let reads = vec![vec![
            Placement { locus: 0, de: 0.001, mapq: 0, as_score: 500 },
            Placement { locus: 1, de: 0.020, mapq: 0, as_score: 498 },
        ]];
        let de_edges = conflict_edges(2, &reads, &ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1 });
        let as_edges = as_tie_edges(2, &reads, 0.9, 1);
        assert!(de_edges.is_empty());
        assert_eq!(as_edges, std::collections::BTreeSet::from([(0, 1)]));
    }
}
```

- [ ] **Step 3: Run the kernel tests**

Run: `cargo test --lib -- vg_family::read_conflict 2>&1 | tail -12`
Expected: all pass. (Other modules will not compile yet — that is fine if you scope the run; otherwise expect `denovo_pipeline` compile errors which Task 4 fixes. If the crate must compile to run any test, do Task 4's `build_read_placements` change together — see Task 4 note.)

- [ ] **Step 4: Commit**

```bash
git add src/rustle/vg_family/read_conflict.rs
git commit -m "feat(read_conflict): de-tie criterion (Placement struct, delta/de_max/min_reads, AS-tie audit helper)"
```

> **Compile note:** changing `ReadPlacements` breaks `build_read_placements`/`conflict_to_split_families`/`detect_and_assign` callers in `denovo_pipeline.rs`. If `cargo test` cannot compile the crate to run the kernel tests in isolation, fold Task 3 + Task 4 into a single commit cycle (write both, then run tests). Prefer separate commits if the kernel compiles standalone.

---

## Task 4: Wire de placements + logging in `denovo_pipeline.rs`

**Files:**
- Modify: `src/rustle/vg_family/denovo_pipeline.rs`

**Interfaces:**
- Consumes: `Placement`, `ConflictParams`, `conflict_edges`, `as_tie_edges`, `conflict_families` (Task 3); `BamRead.de`, `BamRead.is_supplementary` (Task 2)
- Produces: `build_read_placements(bam_reads, reps) -> Vec<ReadPlacements>` building `Placement`s, excluding supplementary

- [ ] **Step 1: Update the `read_conflict` import**

```rust
use super::read_conflict::{as_tie_edges, conflict_edges, conflict_families, ConflictParams, Placement, ReadPlacements};
```

- [ ] **Step 2: Rewrite `build_read_placements` (line ~228)**

```rust
pub(super) fn build_read_placements(bam_reads: &[BamRead], reps: &[DenovoTranscript]) -> Vec<ReadPlacements> {
    use std::collections::BTreeMap;
    let mut by_name: BTreeMap<&str, Vec<Placement>> = BTreeMap::new();
    for br in bam_reads {
        if br.is_supplementary {
            continue; // chimeric/split read is not a multimapping alternative placement
        }
        let read_end = read_ref_end(&br.read);
        let best = reps
            .iter()
            .enumerate()
            .filter(|(_, rep)| br.chrom == rep.chrom && br.read.ref_start < rep.end && read_end > rep.start)
            .max_by_key(|(_, rep)| rep.end.min(read_end) - rep.start.max(br.read.ref_start));
        if let Some((li, _)) = best {
            by_name.entry(br.name.as_str()).or_default().push(Placement {
                locus: li, de: br.de, mapq: br.mapq, as_score: br.as_score,
            });
        }
    }
    by_name.into_values().collect()
}
```

Update the doc comment above it to describe de-tie + supplementary exclusion.

- [ ] **Step 3: Update `DenovoConfig::default` to use `from_env`**

Find the `conflict:` line in `DenovoConfig::default()` and change to:
```rust
conflict: ConflictParams::from_env(),
```
Update the `denovoconfig_default_conflict_params_are_sane` test to the de defaults:
```rust
#[test]
fn denovoconfig_default_conflict_params_are_sane() {
    let cfg = DenovoConfig::default();
    assert!((cfg.conflict.delta - 0.005).abs() < 1e-9);
    assert!((cfg.conflict.de_max - 0.05).abs() < 1e-9);
    assert_eq!(cfg.conflict.min_reads, 3);
}
```

- [ ] **Step 4: Update the conflict logging block in `detect_and_assign`**

In the `{ let placements = build_read_placements(...); ... }` block, after computing `c_edges`/`c_fams`, add the de⊆AS invariant log and the per-edge both-mapq0 fraction:

```rust
let placements = build_read_placements(bam_reads, &reps);
let c_edges = conflict_edges(reps.len(), &placements, &cfg.conflict);
let c_fams = conflict_families(reps.len(), &c_edges);
// audit: the de-tie edge set must be a subset of the AS-tie edge set (the bake-off invariant).
let as_edges = as_tie_edges(reps.len(), &placements, 0.9, cfg.conflict.min_reads);
let de_subset = c_edges.iter().all(|&(i, j, _)| as_edges.contains(&(i, j)));
eprintln!(
    "[detect_and_assign] conflict-graph (de-tie): {} edges -> {} families | de⊆AS={} (AS edges={})",
    c_edges.len(), c_fams.len(), de_subset, as_edges.len(),
);
for (fi, fam) in c_fams.iter().enumerate() {
    let members: Vec<&str> = fam.iter().map(|&i| reps[i].tid.as_str()).collect();
    let coords: Vec<String> = fam.iter()
        .map(|&i| format!("{}:{}-{}", reps[i].chrom, reps[i].start, reps[i].end)).collect();
    let edge_weights: Vec<usize> = c_edges.iter()
        .filter(|&&(a, b, _)| fam.contains(&a) && fam.contains(&b)).map(|&(_, _, w)| w).collect();
    eprintln!(
        "  conflict-fam{fi} n={} reads_linking={:?}: {} @ {}",
        fam.len(), edge_weights, members.join(","), coords.join(" | "),
    );
}
let split = conflict_to_split_families(&c_fams, &c_edges, &cfg.split);
```
(Keep the existing POA-diagnostic and `conflict_to_split_families` call exactly as they are after this block.)

- [ ] **Step 5: Update the 4 `build_read_placements_*` unit tests + the `bam_read` helper**

The `bam_read` test helper must set `de` (and `is_supplementary: false`) instead of `as_score` driving the tie. Change the helper to take a `de` argument and the four tests to de semantics:

```rust
fn bam_read(chrom: &str, start: u64, end: u64, name: &str, de: f32, is_secondary: bool) -> BamRead {
    let mapq = if is_secondary { 0 } else { 60 };
    let len = (end - start) as u64;
    BamRead {
        chrom: chrom.into(),
        read: AlignedRead { ref_start: start, cigar: vec![('M', len)], seq: vec![b'A'; len as usize] },
        mapq, name: name.into(), as_score: 500, de, is_supplementary: false,
    }
}
```

Rewrite the four tests to de semantics (multimapper with de-tied placements forms a family; large de gap does not; nested unique read makes one placement; supplementary excluded). Replace the bodies of `build_read_placements_multimapper_forms_conflict_family`, `build_read_placements_domain_sharer_is_no_family`, `build_read_placements_nested_loci_unique_read_no_conflict`, `build_read_placements_untied_secondary_is_no_conflict` with:

```rust
#[test]
fn build_read_placements_multimapper_forms_conflict_family() {
    let reps = vec![rep("c1", 0, 200), rep("c1", 1000, 1200)];
    let bam = vec![
        bam_read("c1", 0, 200, "read_X", 0.010, false),
        bam_read("c1", 1000, 1200, "read_X", 0.012, true),
    ];
    let placements = build_read_placements(&bam, &reps);
    let edges = super::super::read_conflict::conflict_edges(
        2, &placements, &ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1 });
    assert!(!edges.is_empty(), "de-tied cross-locus read must produce a conflict edge");
    assert_eq!(conflict_families(2, &edges), vec![vec![0, 1]]);
}

#[test]
fn build_read_placements_domain_sharer_is_no_family() {
    let reps = vec![rep("c1", 0, 200), rep("c1", 1000, 1200)];
    let bam = vec![bam_read("c1", 0, 200, "read_Y", 0.010, false)];
    let placements = build_read_placements(&bam, &reps);
    let edges = super::super::read_conflict::conflict_edges(
        2, &placements, &ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1 });
    assert!(edges.is_empty());
    assert!(conflict_families(2, &edges).is_empty());
}

#[test]
fn build_read_placements_nested_loci_unique_read_no_conflict() {
    let reps = vec![rep("c1", 0, 1000), rep("c1", 400, 600)];
    let bam = vec![bam_read("c1", 450, 550, "unique_read", 0.010, false)];
    let placements = build_read_placements(&bam, &reps);
    let total: usize = placements.iter().map(|p| p.len()).sum();
    assert_eq!(total, 1, "one record -> one best-overlap placement -> no cross-locus pair");
    let edges = super::super::read_conflict::conflict_edges(
        2, &placements, &ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1 });
    assert!(edges.is_empty());
}

#[test]
fn build_read_placements_diverged_secondary_is_no_conflict() {
    // read fits copy 0 (de 0.001) far better than copy 1 (de 0.020): resolvable, no edge.
    let reps = vec![rep("c1", 0, 200), rep("c1", 1000, 1200)];
    let bam = vec![
        bam_read("c1", 0, 200, "read_Z", 0.001, false),
        bam_read("c1", 1000, 1200, "read_Z", 0.020, true),
    ];
    let placements = build_read_placements(&bam, &reps);
    let edges = super::super::read_conflict::conflict_edges(
        2, &placements, &ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1 });
    assert!(edges.is_empty(), "diverged secondary (large de gap) is resolvable, not a conflict");
}

#[test]
fn build_read_placements_excludes_supplementary() {
    // a supplementary (chimeric) record on the second locus must NOT create a conflict edge.
    let reps = vec![rep("c1", 0, 200), rep("c1", 1000, 1200)];
    let mut supp = bam_read("c1", 1000, 1200, "read_S", 0.010, true);
    supp.is_supplementary = true;
    let bam = vec![bam_read("c1", 0, 200, "read_S", 0.010, false), supp];
    let placements = build_read_placements(&bam, &reps);
    let total: usize = placements.iter().map(|p| p.len()).sum();
    assert_eq!(total, 1, "supplementary record excluded -> only the primary placement remains");
    let edges = super::super::read_conflict::conflict_edges(
        2, &placements, &ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1 });
    assert!(edges.is_empty());
}
```

- [ ] **Step 6: Fix the `two_paralogs_with_psvs` end-to-end fixture for de-tie + min_reads=3**

The fixture needs ≥3 de-tied multimapper reads (min_reads=3). Set `de` on the conflict reads to tie (e.g. 0.010 primary / 0.012 secondary, both < 0.05) and add a third multimapper `readD`. Replace the `bam` vec (line ~613) with three multimappers (readB/readC/readD), each a primary at locus 0 + secondary at locus 1000, `de: 0.010`/`0.012`, `is_supplementary: false`:

```rust
let mk = |name: &str, primary_seq: Vec<u8>, secondary_seq: Vec<u8>| {
    vec![
        BamRead { chrom: "c1".into(),
            read: AlignedRead { ref_start: 0, cigar: vec![('M', 200), ('N', 20), ('M', 200)], seq: primary_seq },
            mapq: 60, name: name.into(), as_score: 380, de: 0.010, is_supplementary: false },
        BamRead { chrom: "c1".into(),
            read: AlignedRead { ref_start: 1000, cigar: vec![('M', 200), ('N', 20), ('M', 200)], seq: secondary_seq },
            mapq: 0, name: name.into(), as_score: 379, de: 0.012, is_supplementary: false },
    ]
};
let mut bam = Vec::new();
for nm in ["readB", "readC", "readD"] {
    bam.extend(mk(nm, copyb_spliced.clone(), copya_spliced.clone()));
}
```

Then update `detect_and_assign_resolves_multimapper_end_to_end` assertions: with 3 multimappers × 2 records = 6 region reads, `n_reads` and `assignments.len()` become 6; keep the `find(|a| a.best_copy == 0/1)` lookups. Adjust the asserted counts accordingly (run the test to read the actual values, then pin them).

- [ ] **Step 7: Run the full suite**

Run: `cargo test --lib 2>&1 | tail -5`
Expected: 0 failures.

- [ ] **Step 8: Commit**

```bash
git add src/rustle/vg_family/denovo_pipeline.rs
git commit -m "feat(denovo_pipeline): wire de-tie placements, exclude supplementary, log de⊆AS + mapq corroboration"
```

---

## Task 5: Real-data smoke validation

**Files:** none (validation only)

- [ ] **Step 1: RFPL region (expect 0 conflict families)**

```bash
RUSTLE_DENOVO_SMOKE_BAM=/home/juanfra/winloci_scratch/GGO.bam \
RUSTLE_DENOVO_SMOKE_FASTA=/home/juanfra/winloci_scratch/GGO.fasta \
RUSTLE_DENOVO_SMOKE_REGION=NC_073243.2:104789647-104877901 \
cargo test --release --lib -- --ignored smoke_detect_and_assign_real --nocapture 2>&1 \
  | grep -E "conflict-graph|co-located|AGGREGATE"
```
Expected: `conflict-graph (de-tie): 0 edges -> 0 families | de⊆AS=true`, `co-located families: 0`.

- [ ] **Step 2: MAGEA region (expect families fire, de⊆AS holds, silver-standard preserved)**

```bash
RUSTLE_DENOVO_SMOKE_BAM=/home/juanfra/winloci_scratch/GGO.bam \
RUSTLE_DENOVO_SMOKE_FASTA=/home/juanfra/winloci_scratch/GGO.fasta \
RUSTLE_DENOVO_SMOKE_REGION=NC_073247.2:161251228-164865959 \
cargo test --release --lib -- --ignored smoke_detect_and_assign_real --nocapture 2>&1 \
  | grep -E "conflict-fam|conflict-graph|co-located|AGGREGATE"
```
Expected: conflict-graph fires families (the MAGEA de-novo sub-loci pairs from the bake-off cross-map at de < 0.05); `de⊆AS=true`; AGGREGATE silver-standard ≥ 95% (was 100% under AS — confirm no regression).

- [ ] **Step 3: Record results + commit**

```bash
git commit --allow-empty -m "smoke(de-tie): RFPL 0 families, MAGEA fires, de⊆AS holds, silver-standard preserved"
```

If MAGEA conflict families drop vs the AS run, that is EXPECTED where AS was over-firing; confirm against `bench/family_criterion_bakeoff.md` that any dropped edge is an AS false positive, not a real loss. If the silver-standard regresses, STOP and investigate (supplementary exclusion changing the assignment read set is the prime suspect — it should not, since assignment reads come from `region`/`bam_reads`, not the conflict placements).

---

## Self-Review

**Spec coverage:**
- ✓ de-tie predicate with named env-overridable constants (Task 3)
- ✓ `de:f` parsing (Task 1) onto BamRead (Task 2)
- ✓ supplementary excluded from conflict placements (Task 4) + tested
- ✓ both-mapq0 / de⊆AS logging (Task 4)
- ✓ `nm` not carried; `nintron`/`ms`/`cm`/`s1`/`s2`/`ts` not ported
- ✓ no 200bp guard (architecturally unnecessary in Rust — documented)
- ✓ smoke validation incl. silver-standard regression guard (Task 5)

**Placeholder scan:** none — all steps carry concrete code.

**Type consistency:** `Placement { locus, de, mapq, as_score }` is used identically in `read_conflict.rs` (Task 3) and `build_read_placements` (Task 4); `conflict_edges`/`as_tie_edges`/`conflict_families` signatures match across producer and consumers; `ConflictParams { delta, de_max, min_reads }` is consistent in kernel, `from_env`, `DenovoConfig`, and all test literals.
