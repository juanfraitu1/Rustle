# Reference-Absent (Collapsed) Copy Wiring — Implementation Plan (v1)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Wire read-covariation collapsed-copy discovery into per-read assignment so a paralog absent from the assembly becomes a copy reads can be assigned to — behind a two-stage freeze that keeps the O2 result byte-identical when no collapsed copy is admitted.

**Architecture:** A discovery side (`recover_collapsed_copies` emits candidates instead of counting them; a new `collapsed_copy_to_transcript` overlays the discovered haplotype onto the host's spliced reference to make a synthetic `DenovoTranscript`; an admission gate of RNA-internal discriminators decides copy-vs-DNA-needs). An integration side (a two-stage assignment in the colocated-family loop: Stage 1 = ref-only ≡ O2; Stage 2 re-threads only the Stage-1 abstain pool + single-ref-copy collapsed loci against ref+admitted copies). Everything is gated by a flag and defaults OFF, so the production path is unchanged until explicitly enabled.

**Tech Stack:** Rust, `cargo test --release --lib`; poasta (existing); minimap2 subprocess (existing pattern in `copy_assign_pipeline.rs::minimap2_msa_pair`).

## Global Constraints

- **O2 is byte-frozen.** With the feature flag OFF, every emitted assignment/abundance row MUST be byte-identical to `o2_definitive.assignments.tsv`/`.quant.tsv`. Verified in Task 8.
- **Feature is opt-in:** gated behind `AssignParams`-adjacent config + a CLI flag `--absent-copies` (default OFF). OFF ⇒ no code path change.
- **A copy is admitted ONLY if ALL discriminators pass** (Task 4); any failure ⇒ DNA-needs candidate, never a forced copy.
- **Two-stage freeze (Task 5):** Stage-1-Assigned reads at multi-ref-copy families are frozen; only Tied/Ambiguous reads and single-ref-copy collapsed loci are re-threaded in Stage 2.
- **Locked v1 fork choices** (from the spec §8): consensus-`.seq` = reference+allele OVERLAY; distinctness K = `min_p < alpha`; certificate = FLAG (`discovery_coupled`), not hold-out.
- **Determinism:** all new collections sorted (BTreeMap/sort), no HashMap-into-output (match the `copy_split.rs` convention).
- Tests live in the modified file's `#[cfg(test)] mod tests`. Run with `cargo test --release --lib <name>`.
- Commit after each task. Branch: `vg/flow-capacity-apportionment` (do NOT work on main).

**Pre-req:** commit the verified, uncommitted speedup work first (clean baseline) — see Task 0.

---

## File Structure

| File | Responsibility | Change |
|------|---------------|--------|
| `src/rustle/vg_family/copy_split.rs` | collapsed discovery primitives | NEW `collapsed_copy_to_transcript`, NEW `CollapsedCandidate`, NEW discriminators (`cluster_count`, `strand_symmetric`, `min_p_distinct`) |
| `src/rustle/vg_family/absent_copy.rs` (NEW) | the admission gate + genome-remap | NEW module: `admit_candidate`, `remap_gate` |
| `src/rustle/vg_family/denovo_pipeline.rs` | family assembly + assignment loop | `recover_collapsed_copies` → emit; two-stage freeze in the colocated loop; DNA-needs collection |
| `src/rustle/vg_family/copy_assign.rs` | the Assignment struct | add `discovery_coupled: bool` to `Assignment` |
| `src/bin/copy_assign.rs` | CLI | `--absent-copies` flag; `<out>.dna_needs.tsv` writer |
| `src/rustle/vg_family/mod.rs` | module registration | `pub mod absent_copy;` |

---

## Task 0: Commit the speedup baseline

**Files:** (no source change) — git only.

- [ ] **Step 1: Confirm tests green and working tree is the verified speedup state**

Run: `cargo test --release --lib 2>&1 | tail -2`
Expected: `test result: ok. 659 passed`

- [ ] **Step 2: Commit the speedup work (parallel discover_psvs, POA-skip flag, trio, BAM cache, gated A*, timers)**

```bash
git add src/rustle/vg_family/copy_assign.rs src/rustle/vg_family/copy_assign_pipeline.rs \
        src/rustle/vg_family/denovo_assemble.rs src/rustle/vg_family/denovo_pipeline.rs \
        src/rustle/vg_family/family_graph.rs src/bin/copy_assign.rs \
        bench/PERFORMANCE_OPTIMIZATION.md bench/CLUSTER_FUTURE_WORK.md \
        docs/superpowers/specs/2026-06-27-reference-absent-collapsed-copy-wiring-design.md \
        docs/superpowers/plans/2026-06-27-reference-absent-collapsed-copy-wiring.md
git commit -m "perf: byte-identical copy-assign speedup (POA-skip flag, parallel discover_psvs, spanned-gate, BAM cache) + reference-absent design/plan"
```

Expected: clean commit; `git status` shows only unrelated pre-existing modifications.

---

## Task 1: `collapsed_copy_to_transcript` — the consensus-`.seq` overlay primitive

**Files:**
- Modify: `src/rustle/vg_family/copy_split.rs` (add function + tests near `consensus_haplotype`, ~line 411)

**Interfaces:**
- Consumes: `CopyIsoform{intron_chain, allele_vector}` (copy_split.rs:28); the parallel PSV genome positions `psv_pos: &[u64]` (from `discover_locus_psvs`, copy_split.rs:294); the host `DenovoTranscript` (family_detect.rs:48); `build_spliced_seq` (denovo_assemble.rs:581); `exon_map` (copy_assign_pipeline.rs:62).
- Produces: `fn collapsed_copy_to_transcript(iso: &CopyIsoform, psv_pos: &[u64], host: &DenovoTranscript, genome: &GenomeIndex) -> Option<DenovoTranscript>` — a synthetic transcript whose `.seq` is the host's spliced reference with `iso.allele_vector` overlaid at the PSV positions, satisfying `exon_map(&t).len() == t.seq.len()`.

- [ ] **Step 1: Write the failing test**

```rust
#[test]
fn collapsed_copy_to_transcript_overlays_alleles_at_psv_positions() {
    // host: single-exon transcript on '+', chrom "c1", spliced seq fetched from a tiny genome.
    // Build a GenomeIndex stub via the existing test helper (see how other tests build it);
    // here we exercise the OVERLAY arithmetic with a host whose exon_map is identity.
    let host = DenovoTranscript {
        tid: "H".into(), chrom: "c1".into(), start: 100, end: 110, n_reads: 9,
        strand: '+', introns: vec![], seq: b"AAAAAAAAAA".to_vec(),
    };
    // PSV at genome positions 102 and 107 → spliced offsets 2 and 7 (identity exon_map for a single exon).
    let psv_pos = vec![102u64, 107u64];
    let iso = CopyIsoform {
        intron_chain: vec![],
        allele_vector: vec![Some(b'C'), Some(b'G')],
        read_count: 5,
        identifiable: true,
    };
    // Overlay directly against the host seq (no genome fetch needed for a single-exon identity map):
    let t = collapsed_copy_to_transcript_from_host_seq(&iso, &psv_pos, &host)
        .expect("transcript built");
    assert_eq!(t.seq, b"AACAAAAGAA".to_vec(), "C at offset 2, G at offset 7, rest = host");
    assert_eq!(t.seq.len(), exon_map(&t).len(), "seq/exon_map length invariant");
    assert_eq!(t.chrom, "c1");
    assert_eq!(t.strand, '+');
}
```

> Note: split the function so the overlay arithmetic (`_from_host_seq`, taking the already-spliced host seq + its `exon_map`) is unit-testable without a `GenomeIndex`; the public `collapsed_copy_to_transcript` is a thin wrapper that fetches the host seq via `build_spliced_seq` then calls it.

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --release --lib collapsed_copy_to_transcript_overlays -- --nocapture`
Expected: FAIL — `cannot find function collapsed_copy_to_transcript_from_host_seq`

- [ ] **Step 3: Implement the overlay primitive**

```rust
use crate::vg_family::family_detect::DenovoTranscript;
use crate::vg_family::copy_assign_pipeline::exon_map;

/// Overlay `iso.allele_vector` (parallel to `psv_pos`, genome coords) onto an ALREADY-spliced host sequence
/// `host` (with its `exon_map`), returning a synthetic transcript carrying the collapsed copy's distinguishing
/// bases. Substitution-only (v1): a PSV whose genome position is not in the host's exon map (intron/indel) is
/// SKIPPED and the copy is flagged via `None` return when ANY allele cannot be placed (caller routes to
/// DNA-needs). Forward-genome coords; the host's own strand/RC is already baked into `host.seq`.
fn collapsed_copy_to_transcript_from_host_seq(
    iso: &CopyIsoform,
    psv_pos: &[u64],
    host: &DenovoTranscript,
) -> Option<DenovoTranscript> {
    if iso.allele_vector.len() != psv_pos.len() {
        return None; // parallel-vector invariant violated
    }
    let emap = exon_map(host); // spliced offset -> forward-genome coord
    // invert: forward-genome coord -> spliced offset
    let mut g2o: std::collections::BTreeMap<u64, usize> = std::collections::BTreeMap::new();
    for (off, &g) in emap.iter().enumerate() {
        g2o.insert(g, off);
    }
    let mut seq = host.seq.clone();
    let mut placed = 0usize;
    for (k, &pos) in psv_pos.iter().enumerate() {
        if let Some(base) = iso.allele_vector[k] {
            match g2o.get(&pos) {
                Some(&off) if off < seq.len() => {
                    seq[off] = base.to_ascii_uppercase();
                    placed += 1;
                }
                _ => return None, // PSV not placeable in host exon frame -> DNA-needs (indel/intron)
            }
        }
    }
    if placed == 0 {
        return None; // no distinguishing base placed -> not a usable synthetic copy
    }
    Some(DenovoTranscript {
        tid: format!("AC_{}_{}", host.chrom, host.start),
        chrom: host.chrom.clone(),
        start: host.start,
        end: host.end,
        n_reads: iso.read_count as u32,
        strand: host.strand,
        introns: iso.intron_chain.clone(),
        seq,
    })
}

/// Public wrapper: fetch the host's spliced sequence from the genome (using the COPY's intron chain), then
/// overlay the discovered alleles. Returns None if the host sequence can't be built or any allele can't be placed.
pub fn collapsed_copy_to_transcript(
    iso: &CopyIsoform,
    psv_pos: &[u64],
    host: &DenovoTranscript,
    genome: &crate::genome::GenomeIndex,
) -> Option<DenovoTranscript> {
    use crate::vg_family::denovo_assemble::build_spliced_seq;
    let (seq, strand) = build_spliced_seq(genome, &host.chrom, host.start, host.end, &iso.intron_chain)?;
    let host_spliced = DenovoTranscript { seq, strand, ..host.clone() };
    collapsed_copy_to_transcript_from_host_seq(iso, psv_pos, &host_spliced)
}
```

> If `DenovoTranscript` lacks `#[derive(Clone)]`, add it (it is already `Clone` — used via `.clone()` in `denovo_pipeline.rs`).

- [ ] **Step 4: Run the test to verify it passes**

Run: `cargo test --release --lib collapsed_copy_to_transcript -- --nocapture`
Expected: PASS (2 functions covered).

- [ ] **Step 5: Add the "unplaceable → None" test + commit**

```rust
#[test]
fn collapsed_copy_to_transcript_none_when_allele_unplaceable() {
    let host = DenovoTranscript { tid: "H".into(), chrom: "c1".into(), start: 100, end: 105,
        n_reads: 9, strand: '+', introns: vec![], seq: b"AAAAA".to_vec() };
    let psv_pos = vec![999u64]; // not in host exon frame
    let iso = CopyIsoform { intron_chain: vec![], allele_vector: vec![Some(b'C')], read_count: 5, identifiable: true };
    assert!(collapsed_copy_to_transcript_from_host_seq(&iso, &psv_pos, &host).is_none());
}
```

Run: `cargo test --release --lib collapsed_copy_to_transcript` → PASS.
```bash
git add src/rustle/vg_family/copy_split.rs
git commit -m "feat(absent-copy): collapsed_copy_to_transcript overlay primitive"
```

---

## Task 2: `recover_collapsed_copies` emits candidates (not a count)

**Files:**
- Modify: `src/rustle/vg_family/denovo_pipeline.rs:367-384` (`recover_collapsed_copies`)
- Modify: `src/rustle/vg_family/copy_split.rs` (add `CollapsedCandidate`)

**Interfaces:**
- Produces: `pub struct CollapsedCandidate { pub host_tid: String, pub chrom: String, pub start: u64, pub end: u64, pub iso: CopyIsoform, pub psv_pos: Vec<u64>, pub n_clusters: usize }`
- Produces: `fn recover_collapsed_candidates(reps: &[DenovoTranscript], bam_reads: &[BamRead]) -> Vec<CollapsedCandidate>` (the emit version). KEEP the existing `recover_collapsed_copies(...) -> usize` as a thin wrapper `recover_collapsed_candidates(...).len()`-style for the current `collapsed_copies` metric so no caller breaks.

- [ ] **Step 1: Write the failing test** (in `copy_split.rs` tests, exercising the candidate shape via `split_locus_copies` + `discover_locus_psvs` parallelism)

```rust
#[test]
fn split_and_positions_are_parallel_vectors() {
    // Build N AlignedReads at one locus with two co-varying alleles at 2 positions; assert the
    // discovered PSV positions length == each emitted CopyIsoform.allele_vector length.
    let reads = make_two_copy_locus_reads(); // helper: 6+ reads, A/C split at 2 cols
    let pos = discover_locus_psvs(&reads, 3);
    let copies = split_locus_copies(&reads, 3, 2, 3);
    assert!(copies.len() >= 2, "two identifiable copies");
    for c in &copies {
        assert_eq!(c.allele_vector.len(), pos.len(), "allele_vector parallel to discovered positions");
    }
}
```

- [ ] **Step 2: Run to verify it fails** (helper missing)

Run: `cargo test --release --lib split_and_positions_are_parallel`
Expected: FAIL — `cannot find function make_two_copy_locus_reads`

- [ ] **Step 3: Implement the helper + assert the invariant**, then add `CollapsedCandidate` and `recover_collapsed_candidates`:

```rust
// copy_split.rs
#[derive(Clone, Debug)]
pub struct CollapsedCandidate {
    pub host_tid: String,
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub iso: CopyIsoform,
    pub psv_pos: Vec<u64>,  // parallel to iso.allele_vector (genome coords)
    pub n_clusters: usize,  // # identifiable copies discovered at this locus
}
```

```rust
// denovo_pipeline.rs — replace recover_collapsed_copies body, keep the usize wrapper
fn recover_collapsed_candidates(reps: &[DenovoTranscript], bam_reads: &[BamRead]) -> Vec<CollapsedCandidate> {
    use super::copy_split::{discover_locus_psvs, split_locus_copies, CollapsedCandidate};
    let mut out = Vec::new();
    for rep in reps {
        let reads: Vec<AlignedRead> = bam_reads.iter()
            .filter(|br| br.chrom == rep.chrom && br.read.ref_start < rep.end && read_ref_end(&br.read) > rep.start)
            .map(|br| br.read.clone()).collect();
        if reads.len() < 6 { continue; }
        let psv_pos = discover_locus_psvs(&reads, 3);
        let copies = split_locus_copies(&reads, 3, 2, 3);
        if copies.len() < 2 { continue; }
        let n_clusters = copies.iter().filter(|c| c.identifiable).count();
        // emit every identifiable copy EXCEPT the most-supported one (that is the host/ref copy).
        let mut ids: Vec<&CopyIsoform> = copies.iter().filter(|c| c.identifiable).collect();
        ids.sort_by_key(|c| std::cmp::Reverse(c.read_count));
        for iso in ids.into_iter().skip(1) {
            out.push(CollapsedCandidate {
                host_tid: rep.tid.clone(), chrom: rep.chrom.clone(), start: rep.start, end: rep.end,
                iso: (*iso).clone(), psv_pos: psv_pos.clone(), n_clusters,
            });
        }
    }
    out
}

fn recover_collapsed_copies(reps: &[DenovoTranscript], bam_reads: &[BamRead]) -> usize {
    recover_collapsed_candidates(reps, bam_reads).len()
}
```

- [ ] **Step 4: Run tests** → `cargo test --release --lib split_and_positions_are_parallel` PASS; `cargo build --release --lib` clean.

- [ ] **Step 5: Commit**
```bash
git add src/rustle/vg_family/copy_split.rs src/rustle/vg_family/denovo_pipeline.rs
git commit -m "feat(absent-copy): recover_collapsed_candidates emits CollapsedCandidate (keep usize wrapper)"
```

---

## Task 3: Discovery discriminators — cluster-count, strand-symmetry, min_p-distinctness

**Files:**
- Modify: `src/rustle/vg_family/copy_split.rs` (add `strand_symmetric_spectrum`, `min_p_distinct`)

**Interfaces:**
- Produces: `pub fn strand_symmetric_spectrum(host_alleles: &[Option<u8>], cand_alleles: &[Option<u8>]) -> bool` — true iff the substitution spectrum (host base → candidate base, both in transcription orientation) is NOT consistent with pure A→G (RNA editing). A candidate whose every differing column is A→G (or symmetric T→C) is editing-like ⇒ returns `false` (reject).
- Produces: `pub fn min_p_distinct(cand: &[Option<u8>], reference: &[Option<u8>], error_rate: f64, alpha: f64) -> bool` — true iff, treating the candidate as a read of itself, the identifiability bound `min_p = Π (error_rate/3)` over distinguishing columns is `< alpha` (i.e. ≥ enough distinguishing columns to be certifiably distinct from `reference`).

- [ ] **Step 1: Write failing tests**

```rust
#[test]
fn min_p_distinct_requires_enough_columns() {
    // 1 distinguishing column: min_p = e/3 ≈ 1e-3 -> NOT < alpha=1e-3 -> false
    let cand = vec![Some(b'C'), Some(b'A')];
    let reference = vec![Some(b'A'), Some(b'A')];
    assert!(!min_p_distinct(&cand, &reference, 0.003, 1e-3), "1 column insufficient at alpha=1e-3");
    // 3 distinguishing columns: min_p ≈ (1e-3)^3 << alpha -> true
    let cand3 = vec![Some(b'C'), Some(b'C'), Some(b'C')];
    let ref3 = vec![Some(b'A'), Some(b'A'), Some(b'A')];
    assert!(min_p_distinct(&cand3, &ref3, 0.003, 1e-3), "3 columns sufficient");
}

#[test]
fn strand_symmetric_rejects_pure_a_to_g() {
    // every difference is A->G: editing-like -> reject (false)
    let host = vec![Some(b'A'), Some(b'A')];
    let cand = vec![Some(b'G'), Some(b'G')];
    assert!(!strand_symmetric_spectrum(&host, &cand), "pure A->G is editing-like");
    // mixed spectrum (A->C, C->T): real divergence -> accept (true)
    let host2 = vec![Some(b'A'), Some(b'C')];
    let cand2 = vec![Some(b'C'), Some(b'T')];
    assert!(strand_symmetric_spectrum(&host2, &cand2), "mixed spectrum is divergence-like");
}
```

- [ ] **Step 2: Run to verify they fail.** Run: `cargo test --release --lib min_p_distinct strand_symmetric` → FAIL (functions missing).

- [ ] **Step 3: Implement**

```rust
/// `min_p` identifiability bound (the same construction the assignment gate uses, copy_assign.rs:320):
/// Π over distinguishing columns of (error_rate/3). `< alpha` ⇒ certifiably distinct.
pub fn min_p_distinct(cand: &[Option<u8>], reference: &[Option<u8>], error_rate: f64, alpha: f64) -> bool {
    let eps = (error_rate / 3.0).clamp(0.0, 1.0);
    let mut prod = 1.0f64;
    let mut any = false;
    for (a, b) in cand.iter().zip(reference.iter()) {
        if let (Some(x), Some(y)) = (a, b) {
            if x != y { prod *= eps; any = true; }
        }
    }
    any && prod < alpha
}

/// Reject candidates whose differing columns are ALL A->G (plus-strand) — the RNA-editing signature
/// (Clair3-RNA). PSV alleles are in transcription orientation, so editing shows as A->G uniformly.
/// Returns true (=keep) iff at least one differing column is NOT A->G.
pub fn strand_symmetric_spectrum(host_alleles: &[Option<u8>], cand_alleles: &[Option<u8>]) -> bool {
    let mut diffs = 0usize;
    let mut non_ag = 0usize;
    for (h, c) in host_alleles.iter().zip(cand_alleles.iter()) {
        if let (Some(hb), Some(cb)) = (h, c) {
            if hb != cb {
                diffs += 1;
                if !(*hb == b'A' && *cb == b'G') { non_ag += 1; }
            }
        }
    }
    diffs > 0 && non_ag > 0
}
```

- [ ] **Step 4: Run tests** → PASS.
- [ ] **Step 5: Commit**
```bash
git add src/rustle/vg_family/copy_split.rs
git commit -m "feat(absent-copy): discovery discriminators (min_p-distinct, strand-symmetry)"
```

---

## Task 4: The admission gate module (`absent_copy.rs`) incl. genome-remap

**Files:**
- Create: `src/rustle/vg_family/absent_copy.rs`
- Modify: `src/rustle/vg_family/mod.rs` (add `pub mod absent_copy;`)

**Interfaces:**
- Consumes: `CollapsedCandidate` (Task 2), the discriminators (Task 3), the host `DenovoTranscript`, `GenomeIndex`, `minimap2` (subprocess, mirror `copy_assign_pipeline::minimap2_msa_pair` temp-file pattern).
- Produces: `pub enum Admission { Copy(DenovoTranscript), DnaNeeds(DnaNeedsRecord) }`; `pub struct DnaNeedsRecord { pub chrom: String, pub start: u64, pub end: u64, pub n_clusters: usize, pub reason: String, pub read_count: usize }`; `pub fn admit_candidate(cand: &CollapsedCandidate, host: &DenovoTranscript, genome: &GenomeIndex, p: &AbsentCopyParams) -> Admission`.
- Produces: `pub struct AbsentCopyParams { pub error_rate: f64, pub alpha: f64, pub min_clusters: usize, pub remap_max_identity: f64 }` with `Default` (`error_rate: 0.003, alpha: 1e-3, min_clusters: 3, remap_max_identity: 0.98`).

- [ ] **Step 1: Write the failing test** (gate logic with remap stubbed via an injected closure so it's hermetic)

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::vg_family::copy_split::{CollapsedCandidate, CopyIsoform};
    use crate::vg_family::family_detect::DenovoTranscript;

    fn host() -> DenovoTranscript {
        DenovoTranscript { tid: "H".into(), chrom: "c1".into(), start: 100, end: 110,
            n_reads: 9, strand: '+', introns: vec![], seq: b"AAAAAAAAAA".to_vec() }
    }

    #[test]
    fn admit_rejects_two_cluster_substitution_only_as_dna_needs() {
        let cand = CollapsedCandidate {
            host_tid: "H".into(), chrom: "c1".into(), start: 100, end: 110,
            iso: CopyIsoform { intron_chain: vec![], allele_vector: vec![Some(b'C'), Some(b'C'), Some(b'C')],
                read_count: 5, identifiable: true },
            psv_pos: vec![102, 104, 106], n_clusters: 2, // only 2 clusters
        };
        let p = AbsentCopyParams::default();
        // remap closure returns 0.5 identity (distinct), so ONLY the cluster-count rule fails:
        let got = admit_candidate_with_remap(&cand, &host(), &p, |_seq| Some(0.5));
        match got { Admission::DnaNeeds(r) => assert!(r.reason.contains("clusters")), _ => panic!("expected DnaNeeds") }
    }

    #[test]
    fn admit_accepts_three_cluster_distinct_divergent() {
        let cand = CollapsedCandidate {
            host_tid: "H".into(), chrom: "c1".into(), start: 100, end: 110,
            iso: CopyIsoform { intron_chain: vec![], allele_vector: vec![Some(b'C'), Some(b'T'), Some(b'C')],
                read_count: 5, identifiable: true },
            psv_pos: vec![102, 104, 106], n_clusters: 3,
        };
        let p = AbsentCopyParams::default();
        let got = admit_candidate_with_remap(&cand, &host(), &p, |_seq| Some(0.5));
        match got { Admission::Copy(t) => assert_eq!(t.chrom, "c1"), _ => panic!("expected Copy") }
    }
}
```

- [ ] **Step 2: Run to verify it fails.** Run: `cargo test --release --lib admit_` → FAIL (module missing).

- [ ] **Step 3: Implement `absent_copy.rs`** (the gate, with the host allele vector derived by reading the host's bases at `psv_pos`, the discriminators from Task 3, and a remap that shells minimap2 `-cx asm20 --eqx`; the `_with_remap` variant takes the identity closure for hermetic tests, and `admit_candidate` calls it with the real remap):

```rust
//! Admission gate for reference-ABSENT (collapsed) copies: turns a CollapsedCandidate into either a
//! synthetic copy (admitted) or a flagged DNA-needs candidate. A copy is admitted ONLY if ALL of:
//! (1) >= min_clusters co-varying clusters, (2) min_p-distinct from the host, (3) NOT pure-A->G (editing),
//! (4) the consensus overlay is buildable, (5) genome remap identity < remap_max_identity. Any failure
//! routes to DnaNeeds — never a forced copy. See the spec.
use crate::genome::GenomeIndex;
use crate::vg_family::copy_split::{
    collapsed_copy_to_transcript, min_p_distinct, strand_symmetric_spectrum, CollapsedCandidate,
};
use crate::vg_family::family_detect::DenovoTranscript;

pub struct AbsentCopyParams { pub error_rate: f64, pub alpha: f64, pub min_clusters: usize, pub remap_max_identity: f64 }
impl Default for AbsentCopyParams {
    fn default() -> Self { Self { error_rate: 0.003, alpha: 1e-3, min_clusters: 3, remap_max_identity: 0.98 } }
}

#[derive(Clone, Debug)]
pub struct DnaNeedsRecord { pub chrom: String, pub start: u64, pub end: u64, pub n_clusters: usize, pub reason: String, pub read_count: usize }

pub enum Admission { Copy(DenovoTranscript), DnaNeeds(DnaNeedsRecord) }

fn dna_needs(c: &CollapsedCandidate, reason: &str) -> Admission {
    Admission::DnaNeeds(DnaNeedsRecord { chrom: c.chrom.clone(), start: c.start, end: c.end,
        n_clusters: c.n_clusters, reason: reason.into(), read_count: c.iso.read_count })
}

/// Host allele vector at the candidate's PSV positions = the host's own bases (it is the reference copy).
fn host_alleles_at(host: &DenovoTranscript, psv_pos: &[u64]) -> Vec<Option<u8>> {
    use crate::vg_family::copy_assign_pipeline::exon_map;
    let emap = exon_map(host);
    let mut g2o = std::collections::BTreeMap::new();
    for (off, &g) in emap.iter().enumerate() { g2o.insert(g, off); }
    psv_pos.iter().map(|p| g2o.get(p).and_then(|&o| host.seq.get(o).copied())).collect()
}

pub(crate) fn admit_candidate_with_remap<F: Fn(&[u8]) -> Option<f64>>(
    cand: &CollapsedCandidate, host: &DenovoTranscript, p: &AbsentCopyParams, remap_identity: F,
) -> Admission {
    if cand.n_clusters < p.min_clusters {
        return dna_needs(cand, &format!("<{} clusters (copy-vs-allele needs DNA)", p.min_clusters));
    }
    let host_alleles = host_alleles_at(host, &cand.psv_pos);
    if !min_p_distinct(&cand.iso.allele_vector, &host_alleles, p.error_rate, p.alpha) {
        return dna_needs(cand, "not min_p-distinct from host");
    }
    if !strand_symmetric_spectrum(&host_alleles, &cand.iso.allele_vector) {
        return dna_needs(cand, "pure A->G spectrum (editing-like)");
    }
    let t = match collapsed_copy_to_transcript(&cand.iso, &cand.psv_pos, host, &GENOME_NONE) {
        Some(t) => t, None => return dna_needs(cand, "consensus unplaceable (indel/intron)"),
    };
    match remap_identity(&t.seq) {
        Some(id) if id < p.remap_max_identity => Admission::Copy(t),
        Some(_) => dna_needs(cand, ">=98% remap identity (paralog-leak or het)"),
        None => dna_needs(cand, "no homology on remap"),
    }
}
```

> Implementation notes for the engineer: `collapsed_copy_to_transcript` needs a real `GenomeIndex`; in the
> `_with_remap` test path the candidate's overlay can be built from the host seq directly via the
> `collapsed_copy_to_transcript_from_host_seq` path (expose it `pub(crate)`), so the test does NOT need a
> genome. Refactor `admit_candidate_with_remap` to take the already-built host-spliced transcript, and have the
> public `admit_candidate(cand, host, genome, p)` fetch the host seq then delegate. The real `remap_identity`
> shells minimap2 `-cx asm20 --eqx` of `t.seq` vs the genome (mirror `minimap2_msa_pair`'s temp-file + Cleanup
> pattern, parse `de:f:`/`NM` for identity); return `None` on no alignment.

- [ ] **Step 4: Run tests** → PASS. `cargo build --release --lib` clean.
- [ ] **Step 5: Add the real `admit_candidate` + minimap2 remap, build, commit**
```bash
git add src/rustle/vg_family/absent_copy.rs src/rustle/vg_family/mod.rs
git commit -m "feat(absent-copy): admission gate (clusters/min_p/strand/remap) -> Copy | DnaNeeds"
```

---

## Task 5: Two-stage freeze integration in the colocated-family loop

**Files:**
- Modify: `src/rustle/vg_family/denovo_pipeline.rs` (the colocated loop ~478-601, and `detect_and_assign`'s signature to thread a flag)
- Modify: `src/rustle/vg_family/copy_assign.rs` (add `discovery_coupled: bool` to `Assignment`, default `false`)

**Interfaces:**
- Consumes: `recover_collapsed_candidates` (Task 2), `admit_candidate` (Task 4), `assign_family_detailed` (copy_assign_pipeline.rs).
- Produces: behaviour — when `absent_copies` is OFF, the loop is byte-identical to today; when ON, after Stage-1 `assign_family_detailed(&copies, ...)`, build admitted copies, re-run `assign_family_detailed(&copies_plus_absent, ...)`, and for each read TAKE the Stage-2 result ONLY IF the Stage-1 status is Tied/Ambiguous (or the family had a single ref copy); else keep Stage-1. Mark reads assigned to an absent copy `discovery_coupled = true`.

- [ ] **Step 1: Write the failing test** — a hermetic family with one collapsed locus; assert OFF ⇒ Stage-1 result, ON ⇒ a previously-Tied read resolves to the absent copy with `discovery_coupled=true`, and a Stage-1-Assigned read is unchanged. (Build via the existing `assign_family_detailed` test scaffolding; see `copy_assign_pipeline.rs` tests for how families/reads are constructed.)

```rust
#[test]
fn two_stage_freeze_keeps_stage1_assigned_and_resolves_abstain_pool() {
    // ref family with 2 copies (assigns most reads), plus a collapsed locus producing a 3-cluster absent copy.
    // OFF: assignments == ref-only. ON: a read that was Tied under ref-only becomes Assigned to AC_*,
    //      flagged discovery_coupled; a read Assigned under ref-only is byte-identical.
    // (Construct with the module's existing read/copy builders.)
}
```

- [ ] **Step 2: Run to verify it fails.** Run: `cargo test --release --lib two_stage_freeze` → FAIL.

- [ ] **Step 3: Implement** — add the field, thread `absent_copies: bool` through `detect_and_assign`, and insert the Stage-2 block after the existing `let mut detail = assign_family_detailed(&copies, &region, p);` (line 532). Pseudocode of the insert (engineer fills exact field plumbing):

```rust
// after Stage-1 `detail`:
let mut dna_needs_local: Vec<DnaNeedsRecord> = Vec::new();
if absent_copies {
    let cands = recover_collapsed_candidates(&all_copies, bam_reads); // candidates at this family's loci
    let mut admitted: Vec<DenovoTranscript> = Vec::new();
    for cand in &cands {
        if let Some(host) = all_copies.iter().find(|t| t.tid == cand.host_tid) {
            match admit_candidate(cand, host, genome, &AbsentCopyParams::default()) {
                Admission::Copy(t) => admitted.push(t),
                Admission::DnaNeeds(r) => dna_needs_local.push(r),
            }
        }
    }
    if !admitted.is_empty() {
        let n_ref = all_copies.len();
        let mut copies2_owned = all_copies.clone();
        copies2_owned.extend(admitted.iter().cloned());
        let copies2: Vec<&DenovoTranscript> = copies2_owned.iter().collect();
        let detail2 = assign_family_detailed(&copies2, &region, p);
        // FREEZE: take Stage-2 only for reads Stage-1-classified Tied/Ambiguous (or single ref copy).
        let single_ref = n_ref <= 1;
        for (r1, r2) in detail.results.iter_mut().zip(detail2.results.iter()) {
            let s1 = r1.combined.status;
            let take2 = single_ref || matches!(s1, AssignStatus::Tied | AssignStatus::Ambiguous);
            if take2 && r2.combined.status == AssignStatus::Assigned {
                let mut a = r2.combined.clone();
                if a.best_copy >= n_ref { a.discovery_coupled = true; } // assigned to an absent copy
                r1.combined = a;
            }
        }
    }
}
```

> The engineer threads `absent_copies` from `detect_and_assign`'s signature (add a `bool` param) up to the
> CLI (Task 7). `detail.results`/`detail2.results` are index-aligned because `assign_family_detailed` iterates
> the same `region` reads in order; assert this with a length check.

- [ ] **Step 4: Run** `cargo test --release --lib two_stage_freeze` → PASS; `cargo test --release --lib` → all green (existing tests prove OFF path unchanged).
- [ ] **Step 5: Commit**
```bash
git add src/rustle/vg_family/denovo_pipeline.rs src/rustle/vg_family/copy_assign.rs
git commit -m "feat(absent-copy): two-stage freeze (Stage-2 re-thread of abstain pool only) + discovery_coupled flag"
```

---

## Task 6: Surface DNA-needs candidates + the CLI flag

**Files:**
- Modify: `src/bin/copy_assign.rs` (add `--absent-copies` flag; collect `DnaNeedsRecord`s; write `<out>.dna_needs.tsv`)
- Modify: `src/rustle/vg_family/denovo_pipeline.rs` (return the per-family `DnaNeedsRecord`s from `detect_and_assign` alongside the families)

**Interfaces:**
- Consumes: `DnaNeedsRecord` (Task 4); the `absent_copies` bool (Task 5).
- Produces: CLI `--absent-copies` (default false); a new output file `<out>.dna_needs.tsv` with header `chrom\tstart\tend\tn_clusters\tread_count\treason`.

- [ ] **Step 1: Write the failing test** — a small CLI-level or function-level test asserting `detect_and_assign(..., absent_copies=true)` returns ≥1 `DnaNeedsRecord` for a 2-cluster substitution-only locus, and an empty vec when `absent_copies=false`.
- [ ] **Step 2: Run → FAIL** (return type lacks the records).
- [ ] **Step 3: Implement** — extend `detect_and_assign`'s return to `(Vec<FamilyAssignment>, Vec<FallbackEdge>, Vec<DnaNeedsRecord>)`; thread the flag from the CLI; in the bin, after the region loop, write `<out>.dna_needs.tsv`. Mirror the existing TSV writers in `copy_assign.rs` (e.g. the `quant.tsv` writer).
- [ ] **Step 4: Run** the new test + `cargo test --release --lib` → green; build the bin `cargo build --release --bin copy_assign`.
- [ ] **Step 5: Commit**
```bash
git add src/bin/copy_assign.rs src/rustle/vg_family/denovo_pipeline.rs
git commit -m "feat(absent-copy): --absent-copies CLI flag + <out>.dna_needs.tsv"
```

---

## Task 7: Acceptance — O2 byte-identical with flag OFF + a sim5x admit/reject check

**Files:**
- Create: `bench/absent_copy_accept.sh` (the acceptance runner)

**Interfaces:** consumes the built `copy_assign` binary; the data in `/home/juanfra/winloci_scratch`.

- [ ] **Step 1: Flag-OFF byte-identical guard.** Re-run the O2 substrate WITHOUT `--absent-copies` and diff vs the committed `o2_definitive.assignments.tsv`:

```bash
cd /home/juanfra/winloci_scratch
/mnt/c/Users/jfris/Desktop/Rustle/target/release/copy_assign \
  --bam GGO_mm.bam --fasta GGO.fasta --regions o2_regions.txt --min-copies 2 --skip-poa-diagnostic --out o2_off
diff <(sort o2_definitive.assignments.tsv) <(sort o2_off.assignments.tsv) | wc -l   # MUST be 0
```
Expected: `0`.

- [ ] **Step 2: Flag-ON non-regression.** Run WITH `--absent-copies`; require silver non-decreasing and the assigned set a superset:

```bash
/mnt/c/.../copy_assign --bam GGO_mm.bam --fasta GGO.fasta --regions o2_regions.txt --min-copies 2 \
  --skip-poa-diagnostic --absent-copies --out o2_on
grep "silver-standard" o2_on.log   # >= 24660/24682
# every O2-Assigned read must still be Assigned (freeze guarantee): report violations (expect 0)
python3 bench/absent_copy_check_superset.py o2_definitive.assignments.tsv o2_on.assignments.tsv
```
Expected: silver ≥ old; 0 Assigned→Tied flips at frozen reads.

- [ ] **Step 3: sim5x admit/reject.** Inject a known collapsed copy and a known diploid het into the sim5x harness (`bench/sim_reads.py`); assert the collapsed copy is admitted (`AC_*` appears, its reads assigned) and the het is in `dna_needs.tsv`, NOT admitted.

- [ ] **Step 4: Commit the runner**
```bash
git add bench/absent_copy_accept.sh bench/absent_copy_check_superset.py
git commit -m "test(absent-copy): acceptance — O2 byte-identical OFF, non-regression + sim5x admit/reject ON"
```

---

## Self-Review

**Spec coverage:** §4.1 freeze → Task 5; §4.2 consensus-seq → Task 1; §4.3 cross-frame seed → folded into Task 4 `host_alleles_at` (reads host bases at the candidate PSV positions = the shared frame); §4.4 admission gate → Tasks 3+4; §4.5 certificate flag → Task 5 (`discovery_coupled`); §3 DNA-needs state → Tasks 4+6; §6 acceptance → Task 7. Scope §5 (divergent route, VG-indel) correctly excluded. **Gap noted:** §4.4.5 indel/QV mask before clustering is NOT yet a task — it lives inside `discover_locus_psvs` and is a refinement; ADD as Task 4b if the sim5x reject test (Task 7.3) shows homopolymer FPs, else defer to v2. Flagged here so it isn't silently dropped.

**Placeholder scan:** integration Tasks 5/6 give pseudocode for the wiring (exact field plumbing depends on surrounding code the implementer sees); the NOVEL logic (Tasks 1,3,4) has complete code. This is the intended division — the review loop in subagent-driven-development catches integration mismatches.

**Type consistency:** `CollapsedCandidate` (Task 2) → consumed by `admit_candidate` (Task 4); `Admission::Copy(DenovoTranscript)` → consumed by Task 5; `discovery_coupled` added in Task 5 to `Assignment` (copy_assign.rs:52). `min_p_distinct`/`strand_symmetric_spectrum` (Task 3) → called in Task 4. Consistent.

---

## Execution Handoff

Plan saved to `docs/superpowers/plans/2026-06-27-reference-absent-collapsed-copy-wiring.md`. Two execution options:
1. **Subagent-Driven (recommended)** — fresh subagent per task, spec+quality review between tasks.
2. **Inline Execution** — batch with checkpoints.
