# VG-HMM novel-copy — next-steps brief

Companion to `2026-04-28-vg-novel-copy-hmm.md` (the full implementation
plan). Last updated 2026-04-29 after Blocker 1 (strand split) and Blocker 2
(forward DP banding) fixes; branch `vg-hmm-novel-copy` at `a9a8f02`.

## Where we are

| Phase | Status | What landed |
|---|---|---|
| 0 — Scaffolding | ✅ done | Branch cut, `poasta 0.1.0` pinned (chose pure-Rust over `rust-spoa` after cmake-4 conflict), `src/rustle/vg_hmm/` module tree. |
| 1 — Family graph | ✅ done + smoke validated (1.6.5) | `family_graph.rs`: types, builder, position-overlap + minimizer-Jaccard clustering. Smoke caught a junction-rounding bug; fixed. 9 tests. |
| 2 — Profile HMM | ✅ done (with three deviation fixes) | `profile.rs`: `ProfileHmm`, `from_singleton`, `from_msa`; **real POA via poasta** (`poa_graph_to_fasta` round-trip); family pseudo-counts (α = 1/(1+n_seq) mixing fraction); strict-majority match-column threshold. 5 tests. |
| 3 — Forward + Viterbi | ✅ done + smoke validated (3.2.5) | `scorer.rs`: `profile_forward_with_boundary`, **boundary-threaded `forward_against_family`** (O(N·L·M) instead of O(L³·M)), `viterbi_path` with traceback. Smoke caught a critical `TransRow::default()` bug (zeros = log P=1 every transition); fixed. **Final smoke: +2307 nat separation between positive and random reads on real GOLGA6L7 data.** 9 tests. |
| 4 — CLI/dispatcher | ✅ done | `Bundle.synthetic`, `Transcript.synthetic`, `Bundle.rescue_class`, `Transcript.rescue_class` fields; three CLI flags `--vg-discover-novel-mode {kmer,hmm}`, `--vg-rescue-diagnostic`, `--vg-rescue-min-loglik`; `discover_novel_copies` dispatcher; `run_rescue` initial stub. |
| 5 — Rescue plumbing | ✅ done (+ propagation fix) | `rescue.rs`: prefilter, `run_rescue_in_memory`, `run_rescue_with_bundles`, `synthesize_bundles`, BAM iteration in `run_rescue`. Pipeline wired at `pipeline.rs:8137`. Filter exemptions in `transcript_filter.rs::isofrac_with_summary` and `filter_contained_transcripts`. **Synthetic-flag propagation bug** (transcripts not getting bundle.synthetic) fixed in `pipeline.rs::extract_bundle_transcripts_for_graph` (commit `d1720dc`). 14 tests. |
| 6 — Failure-mode diagnostic | ✅ done | `diagnostic.rs`: `RescueClass` enum (5 buckets + NeedsExternalVerification), `classify_internal` (chain-score + seed-mask), `classify_external` (5-step minimap2 ladder; minimap2 at `/home/juanfra/miniforge3/bin/minimap2`). `rescue_class` GTF attribute via `gtf.rs::write_gtf`. Per-family TSV rescue-class columns in `vg.rs::write_family_report_with_em`. 2 tests. |
| 7.1 — chr19 smoke integration | ✅ done | End-to-end `cargo test --test vg_hmm_integration_golga` against `GGO_19.bam`. Skips if BAM/genome missing. **Runtime: 31 minutes** (perf bottleneck per F1 below). |
| 7.2 — Synthetic ground-truth fixture | ✅ done | `tools/make_vg_hmm_fixture.py` deterministic fixture generator (pure stdlib, no pysam). `test_data/vg_hmm/synthetic_family.fa` (6 KB), `test_data/vg_hmm/synthetic_reads.bam` (5.7 KB). 2 tests. **Note:** synthetic BAM has 0 mapped reads → no family discovery → rescue does not fire. Full rescue-rate assertion deferred to Task 7.4 with a mixed-read fixture. |
| 7.3 — Negative control | ✅ done | `tests/regression/vg_hmm_negative_control.rs`: poly-A reads, high-threshold reads, empty families — all yield 0 candidates. 3 tests. |
| **7.4 — Full benchmark + flip default + docs** | ⏳ **pending — needs (2) and (1) below first** | — |

**81 tests passing.** `cargo build --profile dev-opt` clean. Binaries: `rustle`, `dump_family_graph` (smoke).

## What is empirically true at this checkpoint

- The family-aware HMM, scored over real chr19 BAM data, separates reads
  mapped to GOLGA6L7 from random chr19 reads by **+2307 nats** (likelihood
  ratio ≫ e²⁰⁰⁰). Architectural premise empirically validated.
- HMM-mode end-to-end pipeline runs to completion on chr19 BAM in 31 min,
  produces a non-empty GTF, and exits 0.
- Synthetic-flag propagation is correct; filter exemptions reachable.
- Diagnostic infrastructure (internal classifier + external minimap2)
  is in place; not yet exercised on real data.

## What's blocking Task 7.4

### Blocker 1 — Mixed-strand family rejection ✅ FIXED (commit `9a2cbb6`)
`build_family_graph` rejects mixed-strand families; `is_strand_mirror`'s
filter was too narrow (only same-locus opposite-strand pairs). Fixed by
adding `partition_family_by_strand` in `vg_hmm/rescue.rs`: splits a
mixed-strand `FamilyGroup` into per-strand sub-families before building
the HMM family graph. Sub-families with <2 bundles are dropped.
Result on `GGO_19.bam`: 22 input families → 25 sub-families ready for
HMM scoring; zero `[VG-HMM] skipping family ... mixed strands` messages
(was ~12).

### Blocker 2 — Forward DP perf ✅ MITIGATED (commit `a9a8f02`)
Active-band tracking in `profile_forward_with_boundary_banded` (default
`FORWARD_BANDWIDTH = 100`). Smoke runtime on chr19 dropped from 4m46s
wall to 4m29s wall (1.7× user-time speedup, but wall time dominated by
I/O / BAM iteration). Score separation preserved — actually slightly
improved at +2402 nats (was +2307).

**Known limitation:** downstream nodes in a chain inherit wide
boundary distributions, so the effective band grows through the chain.
Full ~17× speedup would need per-edge band propagation (passing narrow
`r_min/r_max` metadata across edges). Deferred. Acceptable for the
Phase 7.4 benchmark — extrapolated full-BAM runtime ~2.5 hours, not 20.

### Blocker 3 — Synthetic BAM fixture has no mapped reads ⏳ STILL OPEN
The Phase 7.2 fixture exercises the BAM-write path but doesn't fire the
full rescue logic. Phase 7.4 still needs a fixture with realistic input:
several mapped paralog reads (to seed family discovery) + several
unmapped reads matching a perturbed novel copy (to test rescue).

**Workaround:** Phase 7.4 can run directly against `GGO.bam` (full
gorilla IsoSeq, ~10⁵ unmapped reads, real biology) without needing the
synthetic fixture extended. The synthetic fixture extension can stay
deferred to F-future.

## Recommended ordering for next session

Blockers 1 and 2 fixed. Phase 7.4 is ready to run.

1. **Phase 7.4 — Full-BAM benchmark on `GGO.bam`.** ~2.5 hours runtime
   (extrapolated from chr19 wall time × 37). Run from inside a tmux/screen
   session or as a background job:
   ```bash
   ./target/release/rustle -L --vg --vg-discover-novel \
       --vg-discover-novel-mode hmm \
       --vg-report /tmp/fam_hmm.tsv \
       --genome-fasta /mnt/c/Users/jfris/Desktop/GGO.fasta \
       -o /tmp/rustle_hmm_full.gtf \
       /mnt/c/Users/jfris/Desktop/GGO.bam 2>&1 | tee /tmp/rustle_hmm_full.log
   ```
   Build with `cargo build --release` first (the existing `dev-opt`
   profile is fine but `release` is faster).

2. **Compare against `MULTI_COPY_FAMILY_PROOF.md` baseline:**
   ```bash
   python3 tools/classify_family_assembly.py /mnt/c/Users/jfris/Desktop/GGO_genomic.gff \
       "golgin subfamily A member 6||golgin A6 family like" \
       /tmp/rustle_hmm_full.gtf --label Rustle_HMM/GOLGA6
   python3 tools/classify_family_assembly.py /mnt/c/Users/jfris/Desktop/GGO_genomic.gff \
       "golgin subfamily A member 8||golgin A8 family like" \
       /tmp/rustle_hmm_full.gtf --label Rustle_HMM/GOLGA8
   ```
   Baseline: GOLGA6 3 exact, 1 wrong-strand, 9 missing; GOLGA8 5 exact,
   0 wrong-strand, 5 missing. Net 8/33 exact across both.

3. **Decision point:** if HMM ≥ baseline (8/33 exact, 1 wrong-strand)
   AND no regressions on AMY/TBC1D3 families, flip default in
   `src/bin/rustle.rs`:
   ```rust
   #[arg(long = "vg-discover-novel-mode", default_value = "hmm")]  // was "kmer"
   ```
   Otherwise: keep default at "kmer", document the gap, and decide
   whether F-future work (per-edge band propagation, multi-copy POA
   smoke F3, synthetic fixture extension) is worth before re-attempting.

4. **Final docs update.** Update `MULTI_COPY_FAMILY_PROOF.md` (replace
   item #3 "scaffold exists" with shipped status + rescue-class
   breakdown from the TSV). Update `docs/ALGORITHMS.md §10` linking to
   spec + companion. Delete `src/bin/dump_family_graph.rs` smoke binary
   (F5). Squash the `vg-hmm: smoke …` commits into the relevant phase
   commits if a clean PR history is desired.

## Other follow-ups (not blocking 7.4 but worth before genome-wide)

- **F2 — Topological sort instead of insertion order.** Current scorer
  trusts node insertion order as topological. Safe for tandem chains;
  Kahn's algorithm + cycle detection would harden against future input
  shapes.
- **F3 — Multi-copy POA validation on a non-tandem family.** GOLGA6L7
  is all singletons (3 paralogs at non-overlapping positions). The
  `from_msa` POA path has not been exercised on real biological data.
  Pick AMY2A/AMY2B (overlapping paralogs) and run a smoke specifically
  for the multi-copy profile path.
- **F5 — Smoke binary cleanup.** Delete `src/bin/dump_family_graph.rs`
  after Phase 7.4 lands; roll its assertions into a regression test.

## Risk register update (vs end-of-Phase-3 version)

| Risk | Status |
|---|---|
| POA crate maturity (poasta 0.1.0) | ✅ Validated through Phase 6; pinned. |
| Forward DP perf on long reads | ✅ Mitigated (commit `a9a8f02`, F1 banding). Per-edge band propagation deferred for further speedup. |
| Topological-order assumption | Acceptable; F2 deferred. |
| Multi-copy POA path untested on real data | ⚠️ Should land F3 before genome-wide. Not strict 7.4 blocker. |
| Mixed-strand family rejection | ✅ Fixed (commit `9a2cbb6`, partition_family_by_strand). |
| `BundleRead` field drift | ✅ Stable across phases. |
| External minimap2 dependency | ✅ Available; `--vg-rescue-diagnostic` off by default. |
| Synthetic-flag propagation through transcript filtering | ✅ Fixed at end of Phase 5 (commit d1720dc). |

## Resumption checklist

When picking back up:

- [ ] `cd /mnt/c/Users/jfris/Desktop/Rustle && git checkout vg-hmm-novel-copy && git status` — confirm clean tree (HEAD should be `a9a8f02` or descendant).
- [ ] `cargo test --profile quick --lib --test vg_hmm_family_graph --test vg_hmm_profile --test vg_hmm_scorer --test vg_hmm_rescue --test vg_hmm_diagnostic --test vg_hmm_negative_control 2>&1 | grep "test result"` — confirm unit tests pass (75 tests; skip the 30+min `vg_hmm_integration_golga` for fast feedback).
- [ ] `cargo build --release` — release binary for the Phase 7.4 benchmark.
- [ ] Read this file. Pick up at step 1 of "Recommended ordering" above (Phase 7.4 full-BAM benchmark).
- [ ] Open `docs/superpowers/plans/2026-04-28-vg-novel-copy-hmm.md` for the canonical task texts; this file is the index.
