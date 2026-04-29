# VG-HMM novel-copy — next-steps brief

Companion to `2026-04-28-vg-novel-copy-hmm.md` (the full implementation
plan). Snapshot at the end of Phase 3 (boundary-threaded forward +
Viterbi); branch `vg-hmm-novel-copy`; HEAD `a1b64a0`.

## Where we are

| Phase | Status | What landed |
|---|---|---|
| 0 — Scaffolding | ✅ done | Branch cut, `poasta 0.1.0` pinned (chose pure-Rust over `rust-spoa` after cmake-4 conflict), `src/rustle/vg_hmm/` module tree. |
| 1 — Family graph | ✅ done + smoke validated | `family_graph.rs`: `ExonClass`, `JunctionEdge`, `FamilyGraph`; builder (position-overlap + minimizer-Jaccard clustering, junction binding). Smoke (Task 1.6.5) caught a junction-rounding bug; fixed. 9 tests. |
| 2 — Profile HMM | ✅ done (with three deviation fixes) | `profile.rs`: `ProfileHmm` + `from_singleton` + `from_msa`; real POA via `poasta::io::fasta::poa_graph_to_fasta`; family pseudo-counts (α = 1/(1+n_seq) mixing fraction); strict-majority match-column threshold. 7 tests. |
| 3 — Forward + Viterbi | ✅ done + smoke validated | `scorer.rs`: `profile_forward_with_boundary`, `forward_against_family` (single-DP boundary threading, O(N·L·M) instead of O(L³·M)), `viterbi_path` with traceback. Smoke (Task 3.2.5) caught a critical `TransRow::default()` bug (zeros in log space = P=1 every transition); fixed. Smoke now shows +2307 nats separation between positive and random reads on real GOLGA6L7 data. |

**What is empirically true at this checkpoint:**
- The family-aware HMM, scored over the real GGO_19.bam, separates reads
  mapped to GOLGA6L7 from reads mapped elsewhere on chr19 by more than
  2000 nats (likelihood ratio ≫ e²⁰⁰⁰). The architectural premise of
  the design is *empirically validated on real data* — independent of
  the python research spike.
- Total runtime on 89 long reads against the 27-node GOLGA6L7 family is
  4 minutes (≈ 2.7 s/read). Tractable but not fast.

## What's left

### Phase 4 — Config + CLI dispatcher
Three small tasks. Mechanical wiring; no algorithmic risk.

- **Task 4.1** — Add `synthetic: bool` field to `Bundle` in `src/rustle/types.rs`, default `false`. Update every `Bundle { … }` construction site (search shows ~20 sites in tests + production code). If the fan-out is >20, consider adding `Default for Bundle` so existing tests use `..Default::default()` instead of churn.
- **Task 4.2** — Add three CLI flags in `src/bin/rustle.rs` and matching fields on `RunConfig` in `types.rs`: `--vg-discover-novel-mode {kmer|hmm}` (default `kmer` until validated), `--vg-rescue-diagnostic` (default off), `--vg-rescue-min-loglik` (default 30.0).
- **Task 4.3** — Make `vg::discover_novel_copies` a dispatcher: rename current body to `discover_novel_copies_kmer`, add a thin top-level dispatcher that branches on `config.vg_discover_novel_mode`. Add a stub `vg_hmm::rescue::run_rescue` returning `Ok(Vec::new())` so the build passes; rescue itself lands in Phase 5.

### Phase 5 — Rescue plumbing
The structural bridge from "read scored against family" to "synthetic bundle in the assembly pipeline." Five tasks.

- **Task 5.1** — `prefilter_read(seq, family_kmers, k, min_hits)` reusing the existing FNV-1a k-mer index pattern from `vg.rs::discover_novel_copies`. Cheap gate before the expensive HMM forward.
- **Task 5.2** — `run_rescue_in_memory(families, bundles, unmapped_reads, min_loglik)`. Build family graph + fit profiles per family; for each unmapped read passing the prefilter, run `forward_against_family`; record candidates above the log-odds threshold.
- **Task 5.3** — `run_rescue(bam_path, …)` BAM-backed driver. Reuse the 4-bit decoding pattern from `vg.rs::discover_novel_copies` (lines 1485–1530). Stream unmapped records, decode sequences, dispatch to in-memory.
- **Task 5.4** — `synthesize_bundles(fg, rescued_with_paths, min_reads)`. Cluster surviving reads by the Viterbi best-path's implied region; emit `Bundle { synthetic: true, source: "vg_rescue", … }` per cluster. Synthetic bundles need a synthesised `junction_stats` derived from the path.
- **Task 5.5** — Wire synthetic bundles into `pipeline.rs` (after the existing `discover_novel_copies` call, around line 14017). Update `transcript_filter.rs` to exempt `synthetic = true` bundles from `transcript_isofrac` and cross-bundle pairwise-contained filters; do NOT exempt them from splice-consensus or min-coverage.

### Phase 6 — Failure-mode diagnostic
Three tasks. The "rescue diagnostic" from `docs/UNMAPPED_FAMILY_RESCUE.md` — classifies each rescued read into one of five buckets explaining why minimap2 missed it.

- **Task 6.1** — `diagnostic.rs::classify_internal(n_kmer_hits, k, read_len, masked_fraction)` returns `RescueClass`. Internal reproduction of minimap2's chain-score and seed-frequency-mask formulas. Buckets 1 (`BelowThresholdChain`), 2 (`SeedMasked`), and `NeedsExternalVerification`.
- **Task 6.2** — `classify_external(ref_fasta, read_seq)` shells out to `minimap2` with stepped-down parameters (default → `-s 20` → `-f 0.0001` → `-N 50` → `-s 10 -f 0`). Buckets 3 (`Divergent`), 4 (`Structural`), 5 (`ReferenceAbsent`).
- **Task 6.3** — Wire diagnostic into rescue: per-read bucket → `rescue_class` GTF attribute on the synthesised transcript; per-family `rescue_class` aggregate added to `--vg-report` TSV.

### Phase 7 — Integration tests + benchmark + flip default + docs
Four tasks. The validation gate before claiming this works in production.

- **Task 7.1** — Integration test: full Rustle run on `GGO_19.bam` with `--vg-discover-novel-mode hmm`. Skip-if-missing fallback for CI machines without the BAM.
- **Task 7.2** — Synthetic ground-truth fixture: take one annotated GOLGA copy, perturb 5% bases, insert a copy-specific exon, simulate 50 reads, assert ≥80% rescued and classed as bucket 3/4 + external minimap2 confirms.
- **Task 7.3** — Negative control: random low-complexity unmapped reads → rescue rate ≤ 1%.
- **Task 7.4** — Full-BAM benchmark on `GGO.bam` (1.67 GB). Compare HMM mode against `MULTI_COPY_FAMILY_PROOF.md` baseline (8/33 exact-match GOLGA paralogs). If HMM ≥ baseline, flip the `--vg-discover-novel-mode` default from `kmer` to `hmm`. Update `MULTI_COPY_FAMILY_PROOF.md` and `docs/ALGORITHMS.md §10`.

## Non-task follow-ups

These are not in the plan but should be considered before flipping the default mode in Phase 7.4.

### F1 — Performance: band the profile forward
Current scorer is O(N·L·M) ≈ 4 min / 89 reads on GOLGA6L7. For genome-wide families × ~10⁵ unmapped reads, this is ~hours. Two cheap mitigations:
- **Banded forward** within `profile_forward_with_boundary`: only update DP cells where |column − read_position| < bandwidth. Bandwidth ≈ 2 × max-deletion-tolerance is sufficient for biological reads. Should give 5–10× speedup with negligible accuracy loss.
- **Read-length cap** at the rescue layer: truncate reads to (e.g.) 500 bp before scoring. PacBio HiFi family-similar reads are highly informative in their first kilobase; full-length scoring buys little additional discrimination once the read clearly belongs to one family.

Pick one before flipping default. Defer the other.

### F2 — Topological sort instead of insertion order
`forward_against_family` and `viterbi_path` walk nodes in insertion order. This is correct for our typical input (single-strand chains from `build_family_graph`'s sequential exon-class construction) but would silently produce wrong answers if ever a family graph contains back-edges or out-of-order insertion. Add Kahn's algorithm with a debug-assert that no cycles exist.

### F3 — Multi-copy profile validation on a non-tandem family
GOLGA6L7 is a tandem paralog triple — every node ends up `copy_specific = true` (singleton profile), so Phase 2's `from_msa` POA path was never exercised on the smoke checkpoint. Pick a family where copies *do* cluster across positions (overlapping paralogs like AMY2A/AMY2B from `MULTI_COPY_FAMILY_PROOF.md` §AMY) and run a smoke. This validates the POA branch, the family-prior smoothing, and the boundary-threading scorer on *non-singleton* profiles before Phase 7's full-BAM benchmark.

### F4 — Decide on `Bundle::synthetic` opt-in
Phase 4.1's plan adds `synthetic: bool` to `Bundle`. The codebase has 20+ `Bundle { … }` construction sites; the migration could be a `Default` implementation (one commit, every callsite uses `..Default::default()`) or a per-site `synthetic: false` insertion (less invasive but more churn). Pick the approach in Task 4.1 — both work, and the choice is reversible.

### F5 — Smoke binary cleanup
`src/bin/dump_family_graph.rs` and the per-task smoke commits are scaffolded as throwaway. Phase 7.4 should delete the binary and roll the assertions into a regression test. Until then, leave it in for future debugging.

## Risk register

| Risk | Mitigation status |
|---|---|
| POA crate maturity (poasta 0.1.0) | Validated in Phase 2; produces correct MSA rows via `poa_graph_to_fasta`. Pin version. |
| Forward DP performance on long reads | Mitigated to O(N·L·M) by boundary threading; F1 (banding / read-cap) recommended before genome-wide rollout. |
| Topological-order assumption | Currently safe for tandem chains; F2 needed for robustness. |
| Multi-copy POA path untested on real data | F3 explicitly addresses; should land before Phase 7. |
| `BundleRead` field drift | Test fixtures already adapted in Phase 1.2; subsequent phases inherit the working pattern. |
| External minimap2 dependency for Phase 6.2 | Optional (`--vg-rescue-diagnostic` off by default); installable via conda/bioconda. Document in README. |

## Recommended ordering for the next session

1. **Phase 4 (CLI/dispatcher)** — mechanical, ~1 day. Done with one batched implementer subagent and one combined review.
2. **Follow-up F1 (banding) OR F2 (Kahn's sort)** before Phase 5 — depends on whether you want to derisk performance or correctness first. F1 is the higher-leverage one.
3. **Phase 5 (rescue plumbing)** — the meatiest phase; ~3 days. Each task gets its own implementer + reviews.
4. **Follow-up F3 (multi-copy POA smoke)** — between Phase 5 and Phase 6.
5. **Phase 6 (diagnostic)** — ~2 days. Phase 6.2 depends on minimap2 being installed.
6. **Phase 7 (integration + benchmark + flip default + docs)** — ~1 day. The Phase 7.4 full-BAM run is the longest-running step (≈30 min on `GGO.bam`).

Total remaining: ~7–8 working days at the current pacing, modulo the flexibility of which follow-ups land where.
