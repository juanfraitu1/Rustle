# VG multimapper copy recovery (family super-graph) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development. Steps use `- [ ]` checkboxes.

**Goal:** Recover paralog copies the `-G` baseline starves, by repositioning the existing family
variation-graph + PSV machinery from an inert additive side-channel into a targeted recovery of
copies ABSENT from the `-G` output, driven by multimapper allele evidence. Behind a new flag
`RUSTLE_VG_RECOVER_COPIES`, default-off, byte-identical when off.

**Architecture:** spec `docs/superpowers/specs/2026-06-14-vg-multimapper-copy-recovery-design.md`.
Family super-graph = `vg_family/family_graph.rs::FamilyGraph` (copies are paths). PSV nodes +
multimapper genotyping (incl. secondaries) + per-copy assignment/assembly =
`vg_family/psv_linkage.rs`. Orchestrator = `vg_family/layer2.rs::run_layer2` (PSV path at lines
~626-634). Validation harness = `bench/copy_recovery_eval/` + `/tmp/cre_guided/headroom.py` (the
23-copy target).

**Tech stack:** Rust (rustle), pysam/SQANTI3 (validation), the GGO chimpanzee long-read dataset.

**Standing constraints:** never `pkill -f rustle`; whole-genome `-L` OOMs → per-chrom serial;
`RUSTLE_PRECISE=1` byte-identical to 4705ab1; commit only when the user asks; stage explicit paths
(no `git add -A`); `RAYON_NUM_THREADS=1` for deterministic runs; test builds `CARGO_PROFILE_TEST_DEBUG=0`.

---

## Phase 0 — Diagnose why the existing PSV path doesn't recover a starved copy

This determines the real fix; do it before writing recovery code.

### Task 0.1: Reproduce a starved flagship copy and trace the PSV path

**Files:** none (investigation). Use `LOC101129569` (4 copies missed; `XM_055375231.2` mm=25/uniq=0).

- [ ] **Step 1:** Get the family's copy loci coords from `/tmp/cre_guided/exons.tsv`
  (`grep -E "XM_055375231.2|XM_055375232.2|XM_055375235.2|XM_055375230.1" /tmp/cre_guided/exons.tsv`).
  Note chrom + the span covering all copies.
- [ ] **Step 2:** Slice that region (+50kb pad) from `/mnt/c/Users/jfris/Desktop/GGO.bam`; confirm it
  contains secondary alignments (`samtools view -c -f 256`).
- [ ] **Step 3:** Run, on the slice, with the PSV super-graph path ON and verbose:
  `RUSTLE_VG_LAYER2_PSV_DEBUG=1 rustle -L slice.bam --vg --vg-layer2 --vg-layer2-psv-linkage -G <stringtie-of-slice>.gtf --genome-fasta GGO.fasta -o out.gtf`
  (generate the `-G` guide with `stringtie -L slice.bam`).
- [ ] **Step 4:** gffcompare `out.gtf` vs the RefSeq for the region; is `XM_055375231.2` recovered (FSM)?
- [ ] **Step 5:** If NOT recovered, pin the failure to ONE of (decision tree → sets the fix task):
  - (a) **family not formed** (the starved copy's locus isn't in the `FamilyGraph`) — because no
    bundle/Layer-1 graph exists there (no primary alignments). → fix in Task 2.x: seed the family with
    the secondary-supported locus.
  - (b) **copy not a path** in the graph (`recover_paralog_path` returns nothing for it). → fix:
    add the copy-path from its reference locus.
  - (c) **genotyping casts no votes for the copy** (`genotype_family_reads` doesn't reach the
    secondaries / PSV columns empty). → fix: ensure the multimapper loci + PSV columns include it.
  - (d) **assigned but not assembled** (`assign_read_to_copy` pins reads but `assemble_psv_isoforms`
    emits nothing, or the `(E)` `family_identifiability` gate rejects the family). → fix: relax/
    correct the gate for this regime; ensure assembly runs on the assigned reads.
  - (e) **assembled but dropped as "already present"** (union-by-chain dedup drops it). → fix: the
    "missed" filter logic (Task 3).
- [ ] **Step 6:** Write the finding (which of a-e, with evidence) to
  `bench/copy_recovery_eval/COPY_RECOVERY_DIAGNOSIS.md`. **This gates the exact code in Phase 2-3.**

---

## Phase 1 — Flag plumbing (mechanical; independent of Phase 0)

### Task 1.1: Add `--vg-layer2-recover-copies` flag

**Files:** Modify `src/rustle/main.rs` (near the other `--vg-layer2-*` args ~594-611), `src/rustle/types.rs`, `src/rustle/pipeline.rs` (run_layer2 call site).

- [ ] **Step 1: Failing test** — `tests/regression/` or a `main.rs` unit test asserting the arg
  parses and sets `RUSTLE_VG_RECOVER_COPIES`. Pattern-match the existing `--vg-layer2-psv-linkage` test.
- [ ] **Step 2:** Run it, confirm it fails (flag absent).
- [ ] **Step 3:** Add `#[arg(long = "vg-layer2-recover-copies", default_value_t = false)] vg_layer2_recover_copies: bool` in `main.rs`; set `std::env::set_var("RUSTLE_VG_RECOVER_COPIES","1")` when true (mirror `--vg-layer2-psv-linkage`); require `--vg`.
- [ ] **Step 4:** Add `recover_copies: bool` param to `run_layer2` (after `psv_filter`); thread it from the pipeline call site (read `RUSTLE_VG_RECOVER_COPIES` or pass the config bool). Update existing `run_layer2` test call sites (pass `false`).
- [ ] **Step 5:** Run tests; build (`cargo build --release`). Commit (when user authorizes).

---

## Phase 2 — Make the super-graph recover the starved copy (Phase-0-directed)

> The exact edits are set by Task 0.1's finding (a-e). The structure below is the additive
> recovery; the implementer fills the specific gap. TDD against the synthetic fixture.

### Task 2.1: Fixture — starved copy recoverable only via secondary allele evidence

**Files:** reuse/extend `bench/sim/` 2-copy IsoSeq fixture (commit 2ec1745) + a `tests/regression/vg_copy_recovery.rs`.

- [ ] **Step 1: Failing test** — load the 2-copy fixture where copy B's reads are primary-at-A but
  carry B's PSV alleles on their secondary alignments. Assert: with `RUSTLE_VG_RECOVER_COPIES=1`,
  `run_layer2` (or the end-to-end pipeline) emits copy B's chain; without it, B absent.
- [ ] **Step 2:** Run, confirm it fails (B not recovered) — reproduces the headroom gap in-test.
- [ ] **Step 3:** Implement the gap fix identified in Task 0.1 (one of a-e) in
  `vg_family/layer2.rs` / `vg_family/psv_linkage.rs` / `vg_family/family_graph.rs`, **gated on
  `recover_copies`** so default behavior is untouched. Keep changes minimal and additive.
- [ ] **Step 4:** Run the fixture test → B recovered with flag, absent without. Run the full
  `vg_family` unit suite for regressions.
- [ ] **Step 5:** Commit (when authorized).

### Task 2.2: "Missed-copy" filter + PSV-decisive phantom guard

**Files:** `src/rustle/vg_family/layer2.rs` (where PSV paths become `novel_transcripts`).

- [ ] **Step 1: Failing tests** — (i) a recovered copy-path whose chain is ALREADY in the baseline
  output is NOT re-emitted (no inflation); (ii) an allele-tied (homology-shadow) copy with no
  PSV-decisive reads is NOT emitted (phantom guard).
- [ ] **Step 2:** Run, confirm failing.
- [ ] **Step 3:** When `recover_copies`, keep a PSV copy-path only if (a) its intron chain is absent
  from the chains already in `novel_transcripts`/baseline AND (b) it has ≥`K` reads phased by
  `assign_read_to_copy` with margin (read `RUSTLE_VG_RECOVER_MIN_DECISIVE`, default 2) AND (c) the
  family passes `family_identifiability` (the `(E)` gate). Tag `source "vg_copy_recovery"`, attach `PsvCertificate`.
- [ ] **Step 4:** Run tests; confirm both guards hold + Task 2.1 fixture still recovers B.
- [ ] **Step 5:** Commit (when authorized).

### Task 2.3: Byte-identical-off + RUSTLE_PRECISE invariants

**Files:** test only.

- [ ] **Step 1:** Test: `run_layer2` with `recover_copies=false` returns identical `novel_transcripts`
  to before this feature (snapshot the fixture output).
- [ ] **Step 2:** chr19 end-to-end: `rustle --vg --vg-layer2 -G st.gtf` (flag off) byte-identical to
  pre-feature; `RUSTLE_PRECISE=1` unaffected. Run, confirm.
- [ ] **Step 3:** Commit (when authorized).

---

## Phase 3 — Genome-wide validation against the 23 (SQANTI3)

### Task 3.1: Recover genome-wide + SQANTI3-validate

**Files:** `bench/` (a small runner, mirror `bench/rc_gen.sh`).

- [ ] **Step 1:** Per-chrom (OOM-safe, serial) `rustle -L $C.bam --vg --vg-layer2 --vg-layer2-recover-copies -G st_$C.gtf --genome-fasta` → `rc2_$C.gtf`; merge (re-ID by chrom).
- [ ] **Step 2:** SQANTI3 the merged output via the copy-recovery harness; run `40_recovery_sets.py` + `50_authenticity_guard.py`.
- [ ] **Step 3:** Diff vs the `gd` baseline recoveries; cross-ref the 23 headroom copies. Report:
  recovered / authentic / phantom counts; how many of the 23; any thin extras.
- [ ] **Step 4:** Write results to `bench/scorecard.md` (new "VG copy recovery" section). Success =
  material fraction of 23 recovered, SQANTI3-FSM + authentic, **0 phantom**.

---

## Final review

- [ ] Dispatch a final code-reviewer over the whole change (spec compliance + quality + the invariants).
- [ ] Confirm: default byte-identical, `RUSTLE_PRECISE` intact, additive (VG ⊇ baseline), 0 phantoms.
- [ ] Use superpowers:finishing-a-development-branch.
