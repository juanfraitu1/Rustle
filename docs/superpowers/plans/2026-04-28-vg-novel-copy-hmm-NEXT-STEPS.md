# VG-HMM novel-copy — next-steps brief

Companion to `2026-04-28-vg-novel-copy-hmm.md` (the full implementation
plan). Last updated 2026-04-29 at the end-of-Phase-7.3 pause; branch
`vg-hmm-novel-copy` at `4f859f3`.

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

Three concerns must be resolved before flipping the `--vg-discover-novel-mode`
default to `hmm`:

### Blocker 1 — All chr19 families rejected as mixed-strand
`build_family_graph` (in `family_graph.rs`) has a strict assertion that
all bundles in a family must share strand. On `GGO_19.bam`, all 22 chr19
families discovered by `discover_family_groups` get rejected. The plan
assumed `discover_family_groups`'s `is_strand_mirror` filter prevents
this — but it doesn't, in practice.

**Investigation needed before Task 7.4:**
- Why does `is_strand_mirror` fail to filter? (Look at `vg.rs:232–258`.)
- Should `build_family_graph` relax the constraint (e.g., split into
  per-strand sub-families)?
- Does GOLGA6L7 specifically pass the check on the full BAM? (We know
  from the GFF that it's all `+` strand at NC_073243.2.)

### Blocker 2 — Forward DP too slow for full-BAM benchmark
Phase 7.1's chr19 integration took 31 minutes. Full `GGO.bam` is 37× the
size: ~20 hours. Need banding before the benchmark is tractable.

**F1 (banding) approach:** add a `bandwidth` parameter to
`profile_forward_with_boundary`; only update DP cells where
`|column − read_position| < bandwidth`. Bandwidth ≈ 2× max-deletion-tolerance
(say 100 bp) is sufficient for biological reads. Expected speedup: 5–10×.

### Blocker 3 — Synthetic BAM fixture has no mapped reads
The Phase 7.2 fixture exercises the BAM-write path but doesn't fire the
full rescue logic. Phase 7.4 needs a fixture with realistic input:
several mapped paralog reads (to seed family discovery) + several
unmapped reads matching a perturbed novel copy (to test rescue).

## Recommended ordering for next session

1. **Investigate mixed-strand rejection (Blocker 1).** ~1–2 hours. Find
   the gap between `is_strand_mirror` and `build_family_graph`'s strand
   check; decide whether to relax or fix the upstream filter. Outcome:
   single commit on the same branch.

2. **Land F1 (banding) (Blocker 2).** ~half-day. Add bandwidth parameter
   to forward DP; preserve unbanded as the default when bandwidth=0
   (backwards compat for existing tests). Outcome: 1–2 commits.

3. **Extend Phase 7.2 fixture for real rescue test (Blocker 3).** ~1
   hour. Generate one or two synthetic "primary" alignments anchoring a
   family, plus several unmapped reads against a perturbed copy. Re-run
   the synthetic ground-truth test with assertions that actually exercise
   rescue (≥80% recovery, rescue_class populated).

4. **Phase 7.4 — Full benchmark + flip default.** ~30–45 min runtime
   after F1, plus ~15 min to compare results vs `MULTI_COPY_FAMILY_PROOF.md`
   baseline. **Decision point:** if HMM ≥ baseline (≥ 8/33 GOLGA-exact
   matches, ≤ 1 wrong-strand) and no regressions on AMY/TBC1D3, flip
   default; otherwise document why and keep `kmer` as default.

5. **Final docs update.** Update `MULTI_COPY_FAMILY_PROOF.md` (replace
   item #3 "scaffold exists" with shipped status + rescue-class
   breakdown). Update `docs/ALGORITHMS.md §10` linking to spec +
   companion. Delete `dump_family_graph` smoke binary (next-steps F5).

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
| Forward DP perf on long reads | ⚠️ **Active blocker** for Task 7.4. F1 mitigates. |
| Topological-order assumption | Acceptable; F2 deferred. |
| Multi-copy POA path untested on real data | ⚠️ Should land F3 before genome-wide. Not strict 7.4 blocker. |
| Mixed-strand family rejection | 🔴 **NEW critical blocker.** Needs investigation in Blocker 1 above. |
| `BundleRead` field drift | ✅ Stable across phases. |
| External minimap2 dependency | ✅ Available; `--vg-rescue-diagnostic` off by default. |
| Synthetic-flag propagation through transcript filtering | ✅ Fixed at end of Phase 5 (commit d1720dc). |

## Resumption checklist

When picking back up:

- [ ] `cd /mnt/c/Users/jfris/Desktop/Rustle && git checkout vg-hmm-novel-copy && git status` — confirm clean tree at `4f859f3` (or whatever HEAD is).
- [ ] `cargo test --profile quick 2>&1 | grep "test result"` — confirm 81 tests still pass.
- [ ] Read this file. Pick up at "Recommended ordering" above.
- [ ] Open `docs/superpowers/plans/2026-04-28-vg-novel-copy-hmm.md` for the canonical task texts; this file is the index.
