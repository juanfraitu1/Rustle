# M-FLOOR — VG ⊇ baseline floor (pulled-forward union) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development. Pulled forward ahead of M2 (user decision 2026-06-13) so VG⊇baseline holds before discovery/merge/amend are built on top. Behind the existing default-off `--vg-layer2` flag.

**Goal:** Guarantee `--vg --vg-layer2` output ⊇ baseline output (every baseline intron chain present, verbatim), additively and deterministically — recovering the chains VG currently drops (chr19 RSTL.589.1 + 18 chrY chains) — while keeping `--vg` (layer2 OFF) byte-identical.

**Root cause (empirically confirmed):** In `--vg`, the pre-assembly family EM (`run_fingerprint_em` / `run_pre_assembly_em_with_snps`, pipeline.rs:11377/11388) reweights multimapper reads toward the best-fit family copy, zeroing them *in their original bundle*. The parallel assembly loop (12441) filters `weight > 0.0`, so a copy's baseline chain whose reads were reweighted away is never extracted. Proof: `RUSTLE_VG_ANCHOR_PRIOR=0 --vg-em-iter 0` (EM off) → VG is a clean superset of baseline on BOTH chr19 (1927 chains, RSTL.589.1 present) and chrY (231 chains). EM is the sole cause; we KEEP it (it powers copy recovery / the +copy wins) and union the dropped baseline chains back additively.

**Why this is simpler than the original M6 plan:** Phase 1 (664919c) already (a) removed secondaries from bundles and (b) eliminated mega-bundle merging — so family bundles are now the SAME clean per-locus bundles as baseline, just *linked* for EM (not merged). Therefore **no re-split is needed**: cloning a family bundle's reads at pre-EM weights and running them through the existing assembly reproduces the baseline chain exactly. The existing `RUSTLE_VG_UNION_BASELINE` infra is reused (`build_fresh_baseline_subbundle`, the holdout at 17614, the union-back at 19567) but its dead secondary-bearing trigger (12392 — never fires post-Phase-1) is bypassed.

**Architecture:** Behind `config.vg_layer2`: (1) BEFORE the EM reweighting, clone each family-member bundle's reads (still at pre-EM/baseline weights) into a UnionBaseline shadow bundle; (2) inject the shadow bundles into `bundles_vec` so the real assembly produces their baseline chains; (3) the existing holdout pulls UnionBaseline transcripts out of predcluster and the union-back adds any whose intron chain is absent from VG output — with the longcov floor set to 0.0 under layer2 so low-coverage baseline chains are not re-dropped.

**Tech Stack:** Rust. Verified by `bench/layer2_invariant.sh` (leg A: layer2 == default byte-identity; leg B: baseline ⊆ layer2, promoted to HARD under `LAYER2_STRICT=1`).

---

## Conventions
Same as the main plan: `DetHashMap`/`DetHashSet::default()`; `RAYON_NUM_THREADS=1` for deterministic runs; stage files explicitly (never `git add -A`); commit messages end with `Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>`. `RUSTLE_PRECISE=1` must stay byte-identical to 4705ab1. Whole-genome `-L` OOMs → per-chrom serial only.

## Grounding (verified against the live tree)
- `families_for_em: Vec<crate::vg::FamilyGroup>` is populated at pipeline.rs:11133; `FamilyGroup.bundle_indices: Vec<usize>` (vg.rs:35) are global bundle indices. The EM reweights `&mut bundles` at 11377 (`run_fingerprint_em`) / 11388 (`run_pre_assembly_em_with_snps`); the comment at 11395 confirms "the EM only reweighted reads; nothing has pruned them" — so reads survive (at weight 0), but assembly filters `weight > 0.0`.
- `bundles_vec` is built and the assembly runs at 12441 (`bundles_vec.into_par_iter().try_for_each(|(bundle_idx, mut bundle)|...)`). The dead `RUSTLE_VG_UNION_BASELINE` clone-injection block is at 12386-12436 (its 12392 secondary-bearing trigger never fires post-Phase-1).
- `build_fresh_baseline_subbundle(parent: &Bundle, reads: Vec<BundleRead>, config: &RunConfig) -> Bundle` (pipeline.rs:728) rebuilds a clean bundle from `reads` (start/end/strand/junction_stats recomputed) tagged `RescueClass::UnionBaseline`. `stamp_union_baseline_rescue_class` (770) tags the resulting transcripts (called at 17067/17600).
- Holdout: pipeline.rs:17614 (`if config.vg_mode && env RUSTLE_VG_UNION_BASELINE`) drains UnionBaseline transcripts into `union_baseline_holdout`.
- Union-back: pipeline.rs:19567 — for each held-out tx, skip if `exons.len() < 2`, skip if `t.longcov < min_longcov` (default 2.0 via `RUSTLE_VG_UNION_BASELINE_MIN_LONGCOV`), skip if its intron chain is already in `vg_chains`; else clear rescue_class and push (strictly additive). `Bundle.reads: Vec<BundleRead>`, `Bundle.chrom: String`.

---

## Task MF.1 — Capture pre-EM family-bundle clones and inject them (the producer)

**Files:** Modify `src/rustle/pipeline.rs` (capture before the EM call ~11377; inject in the 12386 region; extend the holdout gate at 17614).

This task has no isolated unit test (it is pipeline-integrated, like the M1.5 wiring); its acceptance test is the end-to-end superset check in MF.3. Verify via the diagnostic log + the harness.

- [ ] **Capture (before the EM reweighting).** Find the point AFTER `families_for_em` is fully populated (pipeline.rs:11133+) and BEFORE the EM block at ~11377 where both `families_for_em` (or the equivalent family list in scope) and `bundles` are accessible. Insert, gated by `config.vg_layer2`:

```rust
        // Layer 2 (M-FLOOR): capture a pre-EM baseline clone of every family-member
        // bundle. The family EM (below) zeroes multimapper read weights in their
        // original bundle, so a copy's baseline chain can vanish from VG output. We
        // clone here — BEFORE reweighting — so the clone keeps baseline (pre-EM)
        // weights; assembling it reproduces the dropped baseline chain. Injected into
        // bundles_vec after EM, held out of predcluster, and unioned back by chain
        // (strictly additive ⇒ VG ⊇ baseline). Gated default-off via --vg-layer2.
        let mut layer2_baseline_clones: Vec<crate::types::Bundle> = Vec::new();
        if config.vg_layer2 {
            let mut seen_idx: crate::types::DetHashSet<usize> = crate::types::DetHashSet::default();
            for fam in families_for_em.iter() {
                for &bi in fam.bundle_indices.iter() {
                    if !seen_idx.insert(bi) {
                        continue; // a bundle may appear in only one family, but be safe
                    }
                    if let Some(b) = bundles.get(bi) {
                        if b.reads.is_empty() {
                            continue;
                        }
                        let reads = b.reads.clone(); // pre-EM weights
                        layer2_baseline_clones.push(build_fresh_baseline_subbundle(b, reads, &config));
                    }
                }
            }
            eprintln!(
                "[layer2-floor] captured {} pre-EM baseline clone(s) of family bundles",
                layer2_baseline_clones.len()
            );
        }
```

> If `families_for_em` is out of scope at the chosen point (it is built inside a nested block), capture instead from whichever in-scope structure holds the family→bundle mapping at that point (e.g. the same list passed to the EM), keyed by `FamilyGroup.bundle_indices`. The invariant: clone every family-member bundle's reads BEFORE `run_fingerprint_em`/`run_pre_assembly_em_with_snps` mutates weights. If the only feasible in-scope capture is after EM, instead snapshot each family read's pre-EM `weight` into a `DetHashMap<(usize/*bundle*/, u64/*read_uid*/), f64>` BEFORE EM and restore those weights on the clones — but prefer the clone-before-EM form above.

- [ ] **Inject (in the 12386 region, after the dead `RUSTLE_VG_UNION_BASELINE` block).** Add, gated by `config.vg_layer2`:

```rust
    if config.vg_layer2 && !layer2_baseline_clones.is_empty() {
        let mut next_idx = bundles_vec.len();
        let n = layer2_baseline_clones.len();
        for sub in std::mem::take(&mut layer2_baseline_clones) {
            bundles_vec.push((next_idx, sub));
            next_idx += 1;
        }
        eprintln!("[layer2-floor] injected {n} baseline-floor clone bundle(s) for assembly");
    }
```

> `layer2_baseline_clones` must be in scope at both the capture point and here. If the capture is inside a nested block that ends before 12386, hoist the `let mut layer2_baseline_clones` declaration to the enclosing scope (declare it `let mut ... = Vec::new();` before the family/EM block, fill it inside, consume it at injection). The clones are already tagged `UnionBaseline` by `build_fresh_baseline_subbundle`, so `stamp_union_baseline_rescue_class` (already called in the assembly loop at 17067/17600) tags their transcripts.

- [ ] **Extend the holdout gate** at pipeline.rs:17614 from:
```rust
    if config.vg_mode && std::env::var_os("RUSTLE_VG_UNION_BASELINE").is_some() {
```
to:
```rust
    if config.vg_mode
        && (config.vg_layer2 || std::env::var_os("RUSTLE_VG_UNION_BASELINE").is_some())
    {
```

- [ ] Build: `cargo build --release` → `Finished release`. Confirm no new warnings beyond pre-existing. Run `./target/release/rustle -L GGO_19.bam --vg --vg-layer2 -o /tmp/mf.gtf 2>&1 | grep '\[layer2-floor\]'` → prints non-zero captured + injected counts.
- [ ] `cargo test --lib` → `0 failed`.
- [ ] Commit:
  `git add src/rustle/pipeline.rs`
  `git commit -m "feat(vg): M-FLOOR producer — pre-EM family-bundle baseline clones injected for assembly" -m "Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"`

## Task MF.2 — Union-back longcov floor under layer2

**Files:** Modify `src/rustle/pipeline.rs` (union-back, ~19577).

- [ ] Change the `min_longcov` default so a low-coverage baseline chain is not re-dropped under layer2. Replace (19577-19580):
```rust
        let min_longcov: f64 = std::env::var("RUSTLE_VG_UNION_BASELINE_MIN_LONGCOV")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(2.0);
```
with:
```rust
        // Layer-2 floor must guarantee VG ⊇ baseline → do NOT re-drop a low-longcov
        // baseline chain. Default the floor to 0.0 under --vg-layer2 (legacy env path
        // keeps its 2.0 default).
        let default_min_longcov = if config.vg_layer2 { 0.0 } else { 2.0 };
        let min_longcov: f64 = std::env::var("RUSTLE_VG_UNION_BASELINE_MIN_LONGCOV")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(default_min_longcov);
```

- [ ] Build `cargo build --release`. `cargo test --lib` → `0 failed`.
- [ ] Commit:
  `git add src/rustle/pipeline.rs`
  `git commit -m "feat(vg): M-FLOOR union floor 0.0 under --vg-layer2 (never re-drop a baseline chain)" -m "Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"`

## Task MF.3 — End-to-end superset + byte-identity + harness strict

**Files:** none (verification) + possibly `bench/layer2_invariant.sh` (no change needed; it already measures leg B and supports `LAYER2_STRICT=1`).

- [ ] Build release. Generate the three GTFs per chrom and assert:
  - **Superset (leg B):** `baseline ⊆ --vg --vg-layer2` on chr19 AND chrY → 0 missing.
    ```bash
    export RAYON_NUM_THREADS=1
    for tag in chr19 chrY; do
      bam=$([ $tag = chr19 ] && echo GGO_19.bam || echo bench/fixtures/chrY_family.bam)
      ./target/release/rustle -L "$bam" -o /tmp/mf_${tag}_base.gtf 2>/dev/null
      ./target/release/rustle -L "$bam" --vg --vg-layer2 -o /tmp/mf_${tag}_l2.gtf 2>/dev/null
      python3 scripts/coord_signature_superset.py /tmp/mf_${tag}_l2.gtf /tmp/mf_${tag}_base.gtf
    done
    ```
    Expected: both print `OK: VG superset baseline (... all present)`. (chr19 recovers RSTL.589.1; chrY recovers all 18.)
  - **Additivity / byte-identity (leg A):** `--vg --vg-layer2` (with the producer present but, to prove gating, also) `--vg` alone — confirm layer2-OFF is byte-identical to pre-MF `--vg`. Since all MF code is `config.vg_layer2`-gated, run `./target/release/rustle -L GGO_19.bam --vg -o /tmp/mf_vgonly.gtf` and diff its coord-signature against the M1-era `/tmp/layer2/chr19_vg_default.gtf` (regenerate if absent) — must be mutual-superset (==).
- [ ] Run the standing harness in STRICT mode: `LAYER2_STRICT=1 RAYON_NUM_THREADS=1 bench/layer2_invariant.sh`. Expected: leg (A) OK both chroms, leg (B) now PASSES (no WIP) both chroms, ending `ALL HARD INVARIANTS PASS`. (Under LAYER2_STRICT the previously-WIP leg B is fatal — it must be green now.)
- [ ] `cargo test --lib` → `0 failed`.
- [ ] Update the memory/plan note that M-FLOOR closed the gap; commit any harness tweak if made:
  `git add bench/layer2_invariant.sh` (only if changed)
  `git commit -m "test(vg): M-FLOOR closes VG⊇baseline gap (chr19+chrY, LAYER2_STRICT green)" -m "Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"`

---

## Risks
- **Shadow clone assembles a different chain than baseline.** Detected by MF.3 leg B (must be 0 missing). If a chain still drops, the clone's reads/weights differ from baseline — investigate `build_fresh_baseline_subbundle` strand/junction recomputation vs the real per-locus bundle.
- **Layer2-OFF not byte-identical.** All MF code is `config.vg_layer2`-gated; MF.3 leg A guards this.
- **Double assembly cost** (family bundles assembled twice under layer2). Acceptable for a correctness floor behind a dev flag; per-chrom serial avoids OOM.
- **longcov floor 0.0 re-admits FPs.** Bounded: union is add-only by exact intron chain against the baseline set — it can only add chains baseline already produces (i.e. matches StringTie), never net-new FPs beyond baseline.
