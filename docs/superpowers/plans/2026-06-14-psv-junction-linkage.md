# PSV↔Junction Linkage for Paralog Isoform Recovery — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking. Each task: implementer → spec review → code-quality review before the next.

**Goal:** Recover harder per-copy splice isoforms in paralog families by using within-molecule PSV→copy linkage as positive junction-assignment evidence, gated to identifiable families — behind `--vg-layer2-psv-linkage` (default off), strictly additive, VG⊇baseline floor preserved.

**Architecture:** Per family in `run_layer2`: (E) identifiability gate → targeted second BAM pass that genotypes each read at the family's PSV columns + captures its junctions → per-read copy assignment (≥N agreeing PSVs) → junction→copy linkage (≥K linked reads) → per-copy isoform assembly into `FamilyPath`s tagged `IsoformSource::PsvLinked` with a `PsvCertificate` → existing additive union-by-chain emission. New, opt-in channel alongside Part A.

**Tech Stack:** Rust (the rustle crate). Reuses `enumerate_diagnostic_sites` (vg.rs), family-graph `per_copy_sequences`/`per_copy_spans` (family_graph.rs), `bam.rs` per-position mismatch-vs-FASTA reconstruction, `build_exons_from_chain`/`FamilyPath`/union-by-chain (vg_family/layer2.rs + pipeline.rs). Tests: `cargo test`. Genome eval: bench/ harness + gffcompare. Spec: `docs/superpowers/specs/2026-06-14-psv-junction-linkage-design.md`.

**Standing constraints (every task):** `RAYON_NUM_THREADS=1` for deterministic runs; `DetHashMap`/`DetHashSet` only (no std hash reaching output); NEVER `git add -A`/`-u` (stage explicit paths); NEVER `pkill -f rustle`; do not touch pre-existing dirty files (`bench/copy_recovery_eval/*`, `tools/stringtie`, etc.); commit messages end with `Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>`; do not push. The full suite has 3 PRE-EXISTING unrelated failures in `tests/regression/vg_fingerprint_em_correctness.rs` — ignore them. `--vg` (layer2 off) and `--vg-layer2` without the new flag MUST stay byte-identical to today.

---

## File Structure

- **Create** `src/rustle/vg_family/psv_linkage.rs` — the (C) engine: `PsvColumn`, `family_identifiability` (E), `ReadGenotype`, per-read assignment, `link_junctions`, per-copy isoform assembly, `PsvCertificate`. One module, clear responsibility, pure where possible (genotyping pass takes a BAM path).
- **Modify** `src/rustle/vg_family/mod.rs` — register `pub mod psv_linkage;`.
- **Modify** `src/rustle/vg_family/layer2.rs` — `IsoformSource::PsvLinked { copy_id }` variant; `FamilyPath` carries an optional `PsvCertificate`; `run_layer2` calls the engine (gated) and appends results to `novel_transcripts`.
- **Modify** `src/rustle/vg.rs` — if needed, expose PSV columns from `enumerate_diagnostic_sites` (a sibling fn returning per-column per-copy alleles + offsets) without changing the existing fn's output.
- **Modify** `src/rustle/main.rs` + `src/rustle/types.rs` — `--vg-layer2-psv-linkage` flag + `RunConfig.vg_layer2_psv_linkage` + env `RUSTLE_VG_LAYER2_PSV_LINKAGE`.
- **Modify** `src/rustle/pipeline.rs` — pass the flag into `run_layer2`; (emission path already exists, no change there).
- **Create** `bench/sim/gen_psvlink_fixture.py` + **modify** `bench/layer2_invariant.sh` — synthetic fixture + harness leg.
- **Create** `bench/gw_psvlink.sh` + `bench/gw_psvlink_aggregate.py` (or extend `bench/gw_attribute.py`) — genome-wide PSV-linkage attribution.

---

## Task 1: `PsvColumn` extraction + `IsoformSource::PsvLinked`

**Files:** Create `src/rustle/vg_family/psv_linkage.rs`; modify `src/rustle/vg_family/mod.rs`, `src/rustle/vg_family/layer2.rs`, `src/rustle/vg.rs`.

A `PsvColumn` is one paralog-distinguishing base position: its per-copy allele and the genomic coordinate of that position **in each copy's own frame** (derived from the diagnostic node's `per_copy_spans` start + the within-exon offset).

- [ ] **Step 1: Write the failing test** (in `psv_linkage.rs` `#[cfg(test)] mod tests`)

```rust
// A 2-copy family with one diagnostic exon node: copy0 allele "A...", copy1 "C..."
// at within-exon offset 3. per_copy_spans copy0=(100,110), copy1=(500,510).
// psv_columns_for_family must return one column with offset 3 ->
//   genomic 103 (copy0, allele 'A') and 503 (copy1, allele 'C').
#[test]
fn psv_columns_map_offset_to_each_copy_frame() {
    let fg = tests_support::two_copy_one_psv_fg(); // helper: builds the fg above
    let cols = psv_columns_for_family(&fg);
    assert_eq!(cols.len(), 1);
    let c = &cols[0];
    assert!(c.per_copy.iter().any(|p| p.copy_id == 0 && p.genomic_pos == 103 && p.allele == b'A'));
    assert!(c.per_copy.iter().any(|p| p.copy_id == 1 && p.genomic_pos == 503 && p.allele == b'C'));
}
```

- [ ] **Step 2: Run it, confirm it fails** (`RAYON_NUM_THREADS=1 cargo test --lib vg_family::psv_linkage` → fails: type/fn missing).

- [ ] **Step 3: Implement `PsvColumn` + `psv_columns_for_family`**

```rust
pub struct PsvCopyAllele { pub copy_id: usize, pub genomic_pos: u64, pub allele: u8 }
pub struct PsvColumn { pub node_idx: usize, pub per_copy: Vec<PsvCopyAllele> }

/// Diagnostic columns: positions within diagnostic exon nodes where >=2 copies
/// carry distinct bases. For each such position, record each present copy's allele
/// and the genomic coordinate IN THAT COPY's frame (per_copy_spans start + offset).
/// Skips positions where copies are insertion/deletion-misaligned (only equal-length
/// per_copy_sequences columns in v1 — documented limitation).
pub fn psv_columns_for_family(fg: &FamilyGraph) -> Vec<PsvColumn> { /* iterate nodes;
  for diagnostic nodes (>=2 distinct per_copy_sequences of equal len), for each offset
  where alleles differ, build a PsvColumn using per_copy_spans[copy].0 + offset. */ }
```
Use `DetHashSet`/sorted iteration so output is deterministic (sort columns by (node_idx, offset); sort per_copy by copy_id).

- [ ] **Step 4: Add `IsoformSource::PsvLinked`** in layer2.rs:

```rust
pub enum IsoformSource { Native, Transferred { donor_copy: usize }, PsvLinked { copy_id: usize } }
```
Update the `source_rank` total-order in `emit_family_isoforms` to rank `PsvLinked` deterministically (e.g. `(2, copy_id)`), preserving the both-flags-off regression anchor (all-Native still reduces to decompose order).

- [ ] **Step 5: Register module** — `pub mod psv_linkage;` in `vg_family/mod.rs`.

- [ ] **Step 6: Run tests, confirm pass; full `cargo test --lib` no new failures. Commit** (stage `psv_linkage.rs`, `mod.rs`, `layer2.rs`, `vg.rs` if touched).

---

## Task 2: (E) family identifiability predictor

**Files:** `src/rustle/vg_family/psv_linkage.rs`.

- [ ] **Step 1: Write the failing tests**

```rust
#[test]
fn identifiable_family_passes() {
    // >=3 distinguishable PSV columns on exons, low error -> pass
    let fg = tests_support::rabl2_like_fg(); // many PSVs
    assert!(family_identifiability(&fg, /*error_rate*/0.001, /*min_psv_columns*/3));
}
#[test]
fn unidentifiable_family_fails() {
    // <3 PSV columns (DAZ-like long identical cores) -> fail
    let fg = tests_support::daz_like_fg(); // 1 PSV
    assert!(!family_identifiability(&fg, 0.001, 3));
}
```

> **Naming note (type consistency):** the family gate parameter is `min_psv_columns`
> (the count of distinguishing PSV *columns* in a family; env `RUSTLE_VG_LAYER2_PSV_FAMILY_MIN`).
> It is DISTINCT from the per-read `min_psv` in Task 3 (`assign_read_to_copy`, the count of
> agreeing PSVs a single read needs; env `RUSTLE_VG_LAYER2_PSV_MIN`). Both default to 3 but
> are different parameters — keep the names distinct in code.

- [ ] **Step 2: Run, confirm fail.**

- [ ] **Step 3: Implement**

```rust
/// (E) gate: a family is identifiable iff it has >= min_psv distinguishable PSV columns
/// on exonic sequence AND those columns are distinguishable above the error floor
/// (a PSV is "distinguishable" only if the allele difference is not plausibly an error;
/// v1: count columns from psv_columns_for_family; an error_rate sanity bound is applied
/// by requiring min_psv distinct columns, since min_psv independent error agreements are
/// improbable at typical error_rate). Returns false for families below the floor (DAZ-like).
pub fn family_identifiability(fg: &FamilyGraph, error_rate: f64, min_psv_columns: usize) -> bool {
    let _ = error_rate; // v1: count-based; error_rate reserved for the coverage-weighted refinement
    psv_columns_for_family(fg).len() >= min_psv_columns
}
```
(Doc-comment the v1 count rule + the refinement, per spec §8.)

- [ ] **Step 4: Run, confirm pass. Commit** (`psv_linkage.rs`).

---

## Task 3: PSV genotyping pass + per-read copy assignment

**Files:** `src/rustle/vg_family/psv_linkage.rs` (+ reuse `bam.rs` reconstruction).

`ReadGenotype` = one read's PSV evidence + junctions. The genotyping pass is a focused second BAM read over the family copies' genomic spans.

- [ ] **Step 1: Write failing tests for the PURE assignment logic** (decouple from BAM I/O — test `assign_read_to_copy` on hand-built `ReadGenotype`s):

```rust
#[test]
fn read_with_sibling_alleles_assigned_to_sibling() {
    // read at copy0's locus but 3 PSVs carry copy1's allele -> assigned copy1
    let g = ReadGenotype { read_name_hash: 7, psv_votes: vec![(1,3),(0,0)], junctions: vec![(160,300)] };
    assert_eq!(assign_read_to_copy(&g, /*min_psv*/3, /*margin*/1), Some(1));
}
#[test]
fn ambiguous_read_unassigned() {
    let g = ReadGenotype { read_name_hash: 8, psv_votes: vec![(0,2),(1,2)], junctions: vec![] };
    assert_eq!(assign_read_to_copy(&g, 3, 1), None); // neither reaches min_psv; tie
}
```

- [ ] **Step 2: Run, confirm fail.**

- [ ] **Step 3: Implement the pure assignment**

```rust
pub struct ReadGenotype { pub read_name_hash: u64, pub psv_votes: Vec<(usize,usize)>, pub junctions: Vec<(u64,u64)> }

/// Confidently assign iff one copy has >= min_psv supporting PSVs AND strictly dominates
/// every sibling by >= margin. Else None (ambiguous — never guessed).
pub fn assign_read_to_copy(g: &ReadGenotype, min_psv: usize, margin: usize) -> Option<usize> { /* ... */ }
```

- [ ] **Step 4: Write the genotyping pass** (integration-style; reuse the BAM reader + `bam.rs` mismatch reconstruction). Signature:

```rust
/// Second BAM pass over `family_loci` (the copies' genomic spans). For each read
/// overlapping a copy, for each PsvColumn position the read covers (in that copy's
/// frame), read the read's base (CIGAR-mapped, reuse bam.rs) and vote for whichever
/// copy's allele it matches. Returns per-read votes + junction chain. Deterministic
/// (reads in BAM order; votes in DetHashMap keyed by read_name_hash).
pub fn genotype_family_reads(
    bam_path: &Path, family_loci: &[(String,u64,u64)], psv_columns: &[PsvColumn],
    genome: &GenomeIndex,
) -> Vec<ReadGenotype>;
```
Reuse `collect_secondary_index_from_bam`'s reader scaffolding (bundle.rs) as the model for the focused pass — but read PRIMARIES too (a read's primary at the well copy is the cleanest genotype). Cover the `RUSTLE_VG_INCLUDE_SECONDARY` parsing conventions.

- [ ] **Step 5: Test the pass on a tiny BAM fixture** (a 3-read synthetic SAM→BAM where one read carries sibling alleles) OR defer the I/O test to the synthetic fixture in Task 7 and keep Task 3 tests on the pure logic; if deferring, state it explicitly in the commit.

- [ ] **Step 6: Run tests, confirm pass; full suite no new failures. Commit** (`psv_linkage.rs`).

---

## Task 4: Junction→copy linkage + per-copy isoform assembly + `PsvCertificate`

**Files:** `src/rustle/vg_family/psv_linkage.rs`, `src/rustle/vg_family/layer2.rs` (FamilyPath certificate field).

- [ ] **Step 1: Write failing tests**

```rust
#[test]
fn junction_assigned_when_k_linked_reads_agree() {
    // 2 reads confidently copy0 both span junction (160,300); K=2 -> assigned copy0
    let gens = vec![
        ReadGenotype{read_name_hash:1, psv_votes:vec![(0,3)], junctions:vec![(160,300)]},
        ReadGenotype{read_name_hash:2, psv_votes:vec![(0,3)], junctions:vec![(160,300)]},
    ];
    let m = link_junctions(&gens, /*min_psv*/3, /*margin*/1, /*k*/2);
    assert_eq!(m.get(&0).unwrap(), &vec![(160,300)]);
}
#[test]
fn junction_below_k_not_assigned() {
    let gens = vec![ReadGenotype{read_name_hash:1, psv_votes:vec![(0,3)], junctions:vec![(160,300)]}];
    assert!(link_junctions(&gens, 3, 1, 2).get(&0).is_none());
}
```

- [ ] **Step 2: Run, confirm fail.**

- [ ] **Step 3: Implement linkage + assembly**

```rust
pub struct PsvCertificate { pub copy_id: usize, pub linked_reads: usize, pub n_psvs: usize, pub copy_posterior: f64 }

/// Confidently-assigned reads vote their junctions to their copy; a junction is kept for
/// a copy iff >= k votes. Returns copy_id -> assigned junctions (deterministic order).
pub fn link_junctions(gens: &[ReadGenotype], min_psv: usize, margin: usize, k: usize)
    -> DetHashMap<usize, Vec<(u64,u64)>>;

/// Assemble each copy's assigned junctions + its linked reads' exon coords into exon
/// chains (reuse build_exons_from_chain on the linked-read subset), keep chains with >=k
/// linked support (chains-only). Returns FamilyPaths tagged PsvLinked{copy_id} with a
/// PsvCertificate (copy_posterior = vote fraction in v1).
pub fn assemble_psv_isoforms(...) -> Vec<FamilyPath>;
```
Add `pub certificate: Option<PsvCertificate>` to `FamilyPath` (default `None` everywhere else; only PsvLinked paths set it). Update all existing `FamilyPath { .. }` constructors to `certificate: None`.

- [ ] **Step 4: Run tests, confirm pass; full suite no new failures. Commit** (`psv_linkage.rs`, `layer2.rs`).

---

## Task 5: Wire into `run_layer2` + flag (output-changing, gated default-off)

**Files:** `src/rustle/main.rs`, `src/rustle/types.rs`, `src/rustle/pipeline.rs`, `src/rustle/vg_family/layer2.rs`.

- [ ] **Step 1: Flag plumbing.** Add `--vg-layer2-psv-linkage` (clap bool default false) in main.rs; `RunConfig.vg_layer2_psv_linkage` (types.rs, default false, asserted off in the existing default test); `config.vg_layer2_psv_linkage = args.vg_layer2_psv_linkage || RUSTLE_VG_LAYER2_PSV_LINKAGE`. Mirror the `--vg-layer2-new-copies` plumbing exactly.

- [ ] **Step 2: Write a wiring test** (in layer2.rs tests) `run_layer2_psv_linkage_gated_and_additive`: a family with PSV-linkable reads → flag OFF emits no `PsvLinked` path; flag ON (and (E) passes) emits one; assert gated-off is unchanged.

- [ ] **Step 3: Wire `run_layer2`.** Add `psv_linkage: bool` param (passed from pipeline.rs ~17809, mirroring `new_copies`). When `psv_linkage && genome.is_some()`, per family: `if family_identifiability(fg, err, min_psv_columns) { let cols = psv_columns_for_family(fg); let gens = genotype_family_reads(bam_path, family_loci, &cols, genome); novel_transcripts.extend(assemble_psv_isoforms(...)); }`. Read env thresholds once: per-read `min_psv` = `RUSTLE_VG_LAYER2_PSV_MIN` (N=3); per-junction `min_linked` = `RUSTLE_VG_LAYER2_PSV_MIN_LINKED` (K=2); family `min_psv_columns` = `RUSTLE_VG_LAYER2_PSV_FAMILY_MIN` (3). The BAM path must be threaded into run_layer2 (add a param) — confirm pipeline.rs has it at the call site.

- [ ] **Step 4: Verify default-off byte-identity.** `RAYON_NUM_THREADS=1 cargo test --lib vg_family` passes; `GGO_GENOME=/mnt/c/Users/jfris/Desktop/GGO.fasta bash bench/layer2_invariant.sh` ALL GREEN (the flag is off in the harness → proves no collateral change).

- [ ] **Step 5: Smoke ON on chrY** (manual): `RUSTLE_LAYER2_DEBUG=1 ... ./target/release/rustle -L bench/fixtures/chrY_family.bam --vg --vg-layer2 --vg-layer2-psv-linkage --genome-fasta GGO.fasta -o /tmp/psv_on.gtf` vs without the flag → confirm additive (on ⊇ off via `scripts/coord_signature_superset.py`), no panic, and report how many `PsvLinked` transcripts emitted. Do NOT judge quality (Task 8).

- [ ] **Step 6: Commit** (main.rs, types.rs, pipeline.rs, layer2.rs).

---

## Task 6: `PsvCertificate` in the rescue report

**Files:** the `RUSTLE_VG_RESCUE_REPORT` emitter (locate via `grep -rn RUSTLE_VG_RESCUE_REPORT src/`).

- [ ] **Step 1:** When a transcript carries a `PsvCertificate`, add columns (copy_id, linked_reads, n_psvs, copy_posterior) to the rescue TSV. TDD with a unit test asserting the row format for a PsvLinked transcript.
- [ ] **Step 2:** Default-off byte-identical (no report unless the env var is set). Commit.

---

## Task 7: Synthetic fixture + harness leg

**Files:** Create `bench/sim/gen_psvlink_fixture.py`; modify `bench/layer2_invariant.sh`.

- [ ] **Step 1: Generator** (model on `bench/sim/gen_starved_fixture.py`). A 2-copy family (copy A, copy B) sharing exon structure; **copy A expresses an exon-skip isoform**; the skip-isoform reads carry ≥3 copy-A-distinguishing PSVs spanning both a PSV and the skip junction. Construct so the skip junction is assignable to copy A **only** via the linked PSV (no PSV → ambiguous). Output `bench/fixtures/sim_psvlink.{fa,reads.fa,bam}` (+ .bai/.fai) via minimap2 (`-N` ≥ copy count). Verify: copy A's skip reads carry copy-A alleles at the PSVs.
- [ ] **Step 2: Harness leg** (leg 9 in `bench/layer2_invariant.sh`, model on leg 8): flag-OFF copyA-skip transcript ABSENT (or only via Part A's looser path — set `RUSTLE_VG_LAYER2_NO_MULTI_ISOFORM=1` to isolate (C)); flag-ON the PSV-linked skip isoform RECOVERED and assigned to copy A; additive; floor holds; graceful SKIP if fixture absent.
- [ ] **Step 3: Load-bearing check** — confirm the leg fails if `min_psv` is set absurdly high (linkage suppressed) or if the PSV check is bypassed (mis-assignment). Run full harness ALL GREEN.
- [ ] **Step 4: Commit** (`git add -f` the fixture binaries; stage the generator + harness explicitly).

---

## Task 8: Genome-wide PSV-linkage attribution

**Files:** Create `bench/gw_psvlink.sh` (or extend `bench/gw_threeway.sh`) + `bench/gw_psvlink_aggregate.py` (or extend `bench/gw_attribute.py`).

- [ ] **Step 1:** Per chrom (serial, resumable; reuse the gw harness pattern), produce a `--vg-layer2-psv-linkage`-ON output `pl_$C.gtf` alongside the existing `vg_$C.gtf`; gffcompare `pl_$C` vs the NCBI ref.
- [ ] **Step 2: Aggregate** — report the **PSV-linkage net-new recoveries**: NCBI-corroborated transcripts that the PSV-linkage channel adds beyond the current VG-layer (the 62-transcript baseline from `bench/gw_attribute.py`), at the strict `=/c` and broad `=/c/j` tiers, grouped by gene; plus the precision delta vs the non-PSV VG run. Confirm 100% StringTie parity + 0 regressions preserved.
- [ ] **Step 3:** Update `bench/scorecard.md` with the PSV-linkage results. Commit (run in background + monitor for the ~30-40 min genome pass, per the established protocol).

---

## Final review

After all tasks: dispatch a final whole-implementation code reviewer over the `psv_linkage` commits (composition, determinism, default-off byte-identity, floor/additivity, FP-gating correctness), then `superpowers:finishing-a-development-branch`.
