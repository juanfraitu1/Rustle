# EM Restriction + PSV-FASTA — Implementation Plan

> **For agentic workers:** execute task-by-task; each task is TDD (failing test → impl → build → test).
> Steps use `- [ ]`. Build with `cargo build --release`; tests `cargo test --release --lib <filter>`.

**Goal:** (A) restrict the VG fingerprint-EM to families where it can change the answer (new default,
opt-out `RUSTLE_VG_EM_LEGACY`); (B) add opt-in `RUSTLE_VG_PSV_FASTA` emitting a per-transcript spliced
FASTA + companion `.psv.tsv` of copy-distinguishing variants.

**Spec:** `docs/superpowers/specs/2026-06-12-em-restriction-and-psv-fasta-design.md`

**Standing constraints:** never `git add -A` (stage explicitly); commit only at task-group ends with
`Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>`; `RUSTLE_PRECISE=1` must stay byte-identical
to 4705ab1 (`bench/mini3/check_precise.sh`); RAYON_NUM_THREADS=1 -p1 for deterministic runs;
`*.bam/*.gtf` gitignored.

---

## Part A — EM restriction

### Task A1: gate predicate (pure helper + unit tests)

**Files:** Modify `src/rustle/vg.rs` (add helper near `run_fingerprint_em` ~4838) + `#[cfg(test)]`.

- [ ] **Failing tests** for `em_family_qualifies(n_sites, n_unique, n_strict, n_tied, k) -> bool`:
  - PSV-rich + tie: `(3,0,0,2,2) == true`
  - 0-site + unique anchor + tie: `(0,1,0,2,2) == true`
  - 0-site + strict anchor + tie: `(0,0,1,2,2) == true`
  - 0-site + no anchor + tie: `(0,0,0,2,2) == false`
  - any anchor + no tie: `(5,3,3,0,2) == false`
  - PSV-rich + no tie: `(9,0,0,0,2) == false`
- [ ] **Implement:** `pub(crate) fn em_family_qualifies(n_sites: usize, n_unique: usize, n_strict: usize, n_tied: usize, k: usize) -> bool { (n_sites >= k || n_unique >= 1 || n_strict >= 1) && n_tied >= 1 }`
- [ ] Build + `cargo test --release --lib vg::` (the new tests pass).

### Task A2: wire the gate + iter==0 early-out + opt-out

**Files:** Modify `src/rustle/vg.rs` (`run_fingerprint_em` per-family loop).

- [ ] In the per-family loop, after the fingerprint is built (`fp.n_sites` in scope, ~vg.rs:4909-4944)
  and before `PreEntry` collection (~4958), obtain the ownership tallies. `compute_copy_ownership`
  (vg.rs:1857) returns per-copy `n_unique/n_strict/n_tied`; aggregate to family totals
  (`any_unique`, `any_strict`, `total_tied`). Reuse the pass already feeding `anchored_mass_per_copy`
  (vg.rs:5120) if it exposes these; otherwise call `compute_copy_ownership` once here.
- [ ] Read `k` from `RUSTLE_VG_EM_MIN_DECISIVE_SITES` (existing, vg.rs:4874, default 2).
- [ ] **Legacy opt-out:** `let legacy = std::env::var_os("RUSTLE_VG_EM_LEGACY").is_some();`
- [ ] Gate: `if !legacy && !em_family_qualifies(fp.n_sites, n_unique_total, n_strict_total, total_tied, k) { push EmResult::default(); continue; }` — abstain leaves reads at incoming 1/NH (byte-identical for that family).
- [ ] **iter==0 early-out:** at the very top of the per-family loop add `if max_iter == 0 { push EmResult::default(); continue; }` so `--vg-em-iter 0` truly disables independent of the anchor prior (fixes the silent override at vg.rs:5139).
- [ ] Add a one-line stderr tally: `[VG] EM gate: ran N / skipped M families (legacy=…)`.
- [ ] Build + `cargo test --release --lib` (suite green).

### Task A3: delete the dead `--vg-em-uniform-prior` CLI arg

**Files:** `src/rustle/main.rs` (arg def ~585-586, config map ~843), `src/rustle/types.rs`
(`vg_em_uniform_prior` field ~993), and any other reference.

- [ ] `grep -rn "vg_em_uniform_prior\|vg-em-uniform-prior" src/` — confirm it is parsed but never read
  by `run_fingerprint_em`.
- [ ] Remove the clap arg, the `RunConfig` field, and the assignment. Build.
- [ ] Verify `target/release/rustle --help` no longer lists it and `--vg-em-uniform-prior` now errors
  (expected — it was a no-op lever).
- [ ] Build + full suite green.

### Task A4: measurement gate (decide default vs opt-in)

**Files:** none (measurement only); record outcome in the spec/memory.

- [ ] chr19: run legacy (`RUSTLE_VG_EM_LEGACY=1`) vs new-default on `GGO_19.bam`
  (`--vg --vg-snp --genome-fasta /mnt/c/Users/jfris/Desktop/GGO.fasta`, RAYON_NUM_THREADS=1 -p1).
  Compare transcript count, renumber-invariant coord-signature struct diff, gffcompare Sn/Pr vs
  `GGO_19.gtf`.
- [ ] chrY: run new-default vs legacy on `bench/multi_copy_eval/ggo_Y.bam` (needs its genome — use
  `GGO.fasta`); confirm the EM still fires on RBMY/TSPY/DAZ families (gate passes) and copy
  over-counting is not re-introduced.
- [ ] **Decision:** if new-default is net-neutral-or-better on chr19 AND chrY paralog wins preserved →
  keep default-on. Else flip to opt-in: rename the gate to require `RUSTLE_VG_EM_RESTRICT=1` and leave
  legacy as default. Record which path was taken.
- [ ] Escape hatch: `bash bench/mini3/check_precise.sh` (byte-identical).

### Commit (Part A)
- [ ] `git add src/rustle/vg.rs src/rustle/main.rs src/rustle/types.rs` + commit.

---

## Part B — PSV FASTA (opt-in, default-off)

### Task B1: spliced-sequence + coordinate-map helpers (new module + unit tests)

**Files:** Create `src/rustle/psv_fasta.rs`; register `pub mod psv_fasta;` in `src/rustle/lib.rs`
(or wherever modules are declared).

- [ ] **Failing tests:**
  - `spliced_sequence(&genome, "chrT", '+', &[(0,3),(6,9)])` returns exon1+exon2 concatenated.
  - same with `'-'` returns the reverse-complement of the `+` result.
  - `genomic_to_transcript_pos(&[(0,3),(6,9)], '+', g)` maps: g in exon1 → offset; g in exon2 →
    len(exon1)+offset; g in intron (e.g. 4) → None; for `'-'` the offset is measured from the 3' end so
    it indexes the emitted (RC) sequence.
- [ ] **Implement:**
  - `pub fn spliced_sequence(genome: &GenomeIndex, chrom: &str, strand: char, exons: &[(u64,u64)]) -> Vec<u8>`
    — sort exons, concat `genome.fetch_sequence` (genome.rs:113); if `'-'`, `reverse_complement`
    (reuse `crate::vg::reverse_complement`, vg.rs:5952).
  - `pub fn genomic_to_transcript_pos(exons_sorted: &[(u64,u64)], strand: char, gpos: u64) -> Option<u64>`
    — walk exons accumulating length; return 0-based offset into the emitted (strand-oriented) sequence;
    `None` if gpos is intronic/out of range. (For `'-'`, transcript_pos = total_len-1 - forward_offset.)
- [ ] Build + `cargo test --release --lib psv_fasta::`.

### Task B2: PSV table recompute + per-(family,copy) map (behind flag)

**Files:** `src/rustle/psv_fasta.rs` (a builder fn) — invoked later from pipeline behind the flag.

- [ ] **Implement** `pub fn build_psv_map(vg_families, bundles, genome, k_min_copies) -> HashMap<(usize,usize), Vec<PsvSite>>`
  where `PsvSite { ref_pos: u64, this_base: u8, sibling_bases: Vec<u8> }`. For each family: rebuild the
  graph via `crate::vg_family::family_graph::build_family_graph(fam, bundles, Some(genome), …)` then
  `crate::vg::enumerate_diagnostic_sites(&fg, fam.bundle_indices.len())`; read
  `per_copy_site_refs[copy_id]` → `(ref_pos, this_base)`, and the sibling bases from the same
  `ExonVariantSite`/per-copy refs. Key by `(family_id, copy_id)` (== transcripts' `(vg_family_id, vg_copy_id)`,
  invariant pipeline.rs:11064-11066).
- [ ] **Unit test** with a 2-copy synthetic family graph + genome: assert the PSV map has the expected
  `(ref_pos, this_base, sibling_bases)` for each copy (reuse the family-graph/genome test helpers from
  `vg_family/tandem.rs`).
- [ ] Build + tests.

### Task B3: emission block (FASTA + .psv.tsv) wired into run()

**Files:** `src/rustle/pipeline.rs` (after `write_gtf`, ~19691-19692), reusing the
`RUSTLE_VG_RESCUE_REPORT` idiom (~19761-19799); orchestration fn in `src/rustle/psv_fasta.rs`.

- [ ] **Implement** `pub fn emit_psv_fasta(all_transcripts, genome, psv_map, gene_tx_no, out_base) -> io::Result<()>`:
  - `<out>.transcripts.fa`: one record per transcript, header
    `>{tx_id} [family_id=.. copy_id=..] n_psv=N`; body = `spliced_sequence(...)` wrapped (e.g. 60/line).
  - `<out>.psv.tsv`: header + one row per PSV of each family transcript:
    `tx_id  chrom  genomic_pos(1-based)  transcript_pos(1-based)  this_base  sibling_bases  family_id  copy_id`,
    using `genomic_to_transcript_pos` (skip intronic PSVs). Non-family tx → no rows, `n_psv=0` in FASTA.
  - Use the same `assign_gene_tx_numbers` mapping the GTF writer uses (gtf.rs:38) so ids match the GTF.
- [ ] Wire into `run()`: `if std::env::var_os("RUSTLE_VG_PSV_FASTA").is_some() { let psv_map = build_psv_map(...); emit_psv_fasta(...)?; }` immediately after `write_gtf`. The genome is `vg_snp_genome` (in scope, pipeline.rs:10403). Recompute is behind the flag → default pays nothing.
- [ ] Build.

### Task B4: integration + byte-identical + showcase

**Files:** test under `tests/regression/` (new) or a bench check script.

- [ ] **Default-off byte-identical:** run `GGO_19.bam` VG with and without `RUSTLE_VG_PSV_FASTA`; assert
  no `.transcripts.fa`/`.psv.tsv` produced when unset and the GTF is identical.
- [ ] **Flag-on integration:** with the flag, assert one FASTA record per GTF transcript and that every
  `.psv.tsv` row's `tx_id` exists in the FASTA and is a family transcript.
- [ ] **Showcase:** run on a chrY paralog family (`ggo_Y.bam`); assert `n_psv>0` rows exist (the
  StringTie-is-blind evidence). Record an example PSV.
- [ ] Escape hatch `bench/mini3/check_precise.sh` (byte-identical).

### Commit (Part B)
- [ ] `git add src/rustle/psv_fasta.rs src/rustle/lib.rs src/rustle/pipeline.rs <tests>` + commit.

---

## Final: adversarial review + land

- [ ] Run an adversarial review workflow over the combined diff (correctness, byte-identical/escape-hatch
  safety, EM-abstain-equals-1/NH, FASTA coord-mapping/RC, recompute join correctness). Fix confirmed
  findings.
- [ ] Full suite green; escape hatch byte-identical; A4 + B4 measurements recorded.
- [ ] Update memory (`project_consensus_correction_wired` sibling or new file) with the EM-restriction
  outcome (default vs opt-in) and the PSV-FASTA flag.
