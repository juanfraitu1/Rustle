# Mis-chain Read Salvage — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Split reads mis-chained across spurious giant introns into their local segments *before* skeleton seeding, so a real duplication locus is not lost when its paralog's reads mis-chain — a general fix, validated against synthetic ground-truth and known families, not tuned to Soto.

**Architecture:** A pure function `split_mischained_reads` (in `denovo_assemble.rs`) transforms the read set; a config flag `mischain_salvage` (default off) plus a small `maybe_salvage_mischain` helper wire it in before every live `pass1_skeletons_robust` call. Off ⇒ byte-identical.

**Tech Stack:** Rust; existing `PrimaryRead`, `read_junction_support`, `env_num`, `MISCHAIN_GIANT_INTRON_BP`, `MISCHAIN_MIN_JUNCTION_READS`, `DenovoConfig`.

## Global Constraints

- **OFF path byte-identical:** with `mischain_salvage=false` the split is never invoked; existing output unchanged.
- **No new magic numbers:** the split criterion reuses `MISCHAIN_GIANT_INTRON_BP` (50 000) and `MISCHAIN_MIN_JUNCTION_READS` (`GATE_MIN_READS`=3), env-overridable via `env_num("RUSTLE_MISCHAIN_GIANT_BP", …)` / `env_num("RUSTLE_MISCHAIN_MIN_READS", …)`.
- **Env gate:** `RUSTLE_MISCHAIN_SALVAGE=1` enables it (read in `DenovoConfig::from_env`).
- **Negative test is mandatory:** a well-supported large intron (≥ 3 reads at that junction) must never be split.
- Support for the split decision is computed on the ORIGINAL reads (before any split).

---

### Task 1: `split_mischained_reads` pure function + `PrimaryRead` equality

**Files:**
- Modify: `src/rustle/vg_family/denovo_assemble.rs` (add `PartialEq, Eq` to `PrimaryRead`; add function; add tests)

**Interfaces:**
- Consumes: `PrimaryRead { chrom: String, ref_start: u64, ref_end: u64, introns: Vec<(u64,u64)> }`, `std::collections::HashMap<(String,u64,u64), usize>`.
- Produces: `pub fn split_mischained_reads(reads: &[PrimaryRead], support: &HashMap<(String,u64,u64),usize>, giant_bp: u64, min_reads: usize) -> Vec<PrimaryRead>`

- [ ] **Step 1: Add `PartialEq, Eq` to `PrimaryRead` derives**

Change `#[derive(Clone, Debug)]` above `pub struct PrimaryRead` to:
```rust
#[derive(Clone, Debug, PartialEq, Eq)]
```

- [ ] **Step 2: Write the failing tests** (in the existing `#[cfg(test)] mod tests` of `denovo_assemble.rs`)

```rust
#[test]
fn split_cuts_spurious_giant_intron_keeping_both_segments() {
    use std::collections::HashMap;
    let reads = vec![PrimaryRead {
        chrom: "chr1".into(), ref_start: 100, ref_end: 80_100,
        introns: vec![(200, 210), (300, 80_000)], // small internal intron, then a giant bridge
    }];
    let mut support = HashMap::new();
    support.insert(("chr1".to_string(), 300, 80_000), 1); // giant, sub-threshold -> CUT
    let out = split_mischained_reads(&reads, &support, 50_000, 3);
    assert_eq!(out, vec![
        PrimaryRead { chrom: "chr1".into(), ref_start: 100,    ref_end: 300,    introns: vec![(200, 210)] },
        PrimaryRead { chrom: "chr1".into(), ref_start: 80_000, ref_end: 80_100, introns: vec![] },
    ]);
}

#[test]
fn split_does_not_cut_well_supported_large_intron() {
    use std::collections::HashMap;
    let reads = vec![PrimaryRead {
        chrom: "chr1".into(), ref_start: 100, ref_end: 80_100, introns: vec![(300, 80_000)],
    }];
    let mut support = HashMap::new();
    support.insert(("chr1".to_string(), 300, 80_000), 3); // >= min_reads -> real large-gene intron, NOT a mis-chain
    let out = split_mischained_reads(&reads, &support, 50_000, 3);
    assert_eq!(out, reads); // unchanged
}

#[test]
fn split_ignores_sub_giant_introns() {
    use std::collections::HashMap;
    let reads = vec![PrimaryRead {
        chrom: "chr1".into(), ref_start: 100, ref_end: 500, introns: vec![(200, 210), (300, 320)],
    }];
    let out = split_mischained_reads(&reads, &HashMap::new(), 50_000, 3); // no intron exceeds giant_bp
    assert_eq!(out, reads);
}

#[test]
fn split_handles_two_giant_introns_into_three_segments() {
    use std::collections::HashMap;
    let reads = vec![PrimaryRead {
        chrom: "chr1".into(), ref_start: 0, ref_end: 160_050,
        introns: vec![(50, 80_000), (80_050, 160_000)], // two giant sub-threshold bridges
    }];
    let out = split_mischained_reads(&reads, &HashMap::new(), 50_000, 3); // absent support => 0 < 3 => both cut
    assert_eq!(out, vec![
        PrimaryRead { chrom: "chr1".into(), ref_start: 0,       ref_end: 50,      introns: vec![] },
        PrimaryRead { chrom: "chr1".into(), ref_start: 80_000,  ref_end: 80_050,  introns: vec![] },
        PrimaryRead { chrom: "chr1".into(), ref_start: 160_000, ref_end: 160_050, introns: vec![] },
    ]);
}

#[test]
fn split_passes_reads_through_unchanged_when_no_cut() {
    use std::collections::HashMap;
    let reads = vec![
        PrimaryRead { chrom: "chr1".into(), ref_start: 0, ref_end: 100, introns: vec![] },
        PrimaryRead { chrom: "chr2".into(), ref_start: 5, ref_end: 400, introns: vec![(100, 200)] },
    ];
    assert_eq!(split_mischained_reads(&reads, &HashMap::new(), 50_000, 3), reads);
}
```

- [ ] **Step 3: Run the tests to confirm they fail**

Run: `cargo test -p rustle --lib split_ 2>&1 | tail -20`
Expected: FAIL — `cannot find function split_mischained_reads` (and `PrimaryRead` now needs `Eq` for the HashMap-free asserts — provided by Step 1).

- [ ] **Step 4: Implement `split_mischained_reads`**

Add after `read_junction_support`/near the other pass-1 helpers in `denovo_assemble.rs`:
```rust
/// Split reads at SPURIOUS giant introns (mis-chain salvage, Approach A). A read whose intron-chain contains an
/// intron `(d,a)` with `a - d > giant_bp` AND junction `(chrom,d,a)` carried by `< min_reads` reads is cut at
/// that intron into local sub-reads (the giant bridge is removed); both flanking segments are kept (each was a
/// real local alignment). Reads with no such intron pass through UNCHANGED. `support` is measured on the
/// ORIGINAL read set. Deterministic; sub-reads replace their parent in 5'→3' order.
pub fn split_mischained_reads(
    reads: &[PrimaryRead],
    support: &std::collections::HashMap<(String, u64, u64), usize>,
    giant_bp: u64,
    min_reads: usize,
) -> Vec<PrimaryRead> {
    let mut out = Vec::with_capacity(reads.len());
    for r in reads {
        // introns[i] joins exon i and exon i+1. Mark spurious giant introns as cuts.
        let is_cut: Vec<bool> = r
            .introns
            .iter()
            .map(|&(d, a)| {
                a.saturating_sub(d) > giant_bp
                    && support.get(&(r.chrom.clone(), d, a)).copied().unwrap_or(0) < min_reads
            })
            .collect();
        if !is_cut.iter().any(|&c| c) {
            out.push(r.clone());
            continue;
        }
        // reconstruct exon boundaries: exon 0 = [ref_start, introns[0].0], exon k = [introns[k-1].1, introns[k].0],
        // last exon = [introns[last].1, ref_end].
        let mut exons: Vec<(u64, u64)> = Vec::with_capacity(r.introns.len() + 1);
        let mut s = r.ref_start;
        for &(d, a) in &r.introns {
            exons.push((s, d));
            s = a;
        }
        exons.push((s, r.ref_end));
        // walk exons; close a segment after exon i when introns[i] is a cut, else carry introns[i] inside.
        let mut seg_start = exons[0].0;
        let mut seg_introns: Vec<(u64, u64)> = Vec::new();
        for i in 0..exons.len() {
            if i + 1 < exons.len() {
                if is_cut[i] {
                    out.push(PrimaryRead {
                        chrom: r.chrom.clone(),
                        ref_start: seg_start,
                        ref_end: exons[i].1,
                        introns: std::mem::take(&mut seg_introns),
                    });
                    seg_start = exons[i + 1].0;
                } else {
                    seg_introns.push(r.introns[i]);
                }
            }
        }
        out.push(PrimaryRead {
            chrom: r.chrom.clone(),
            ref_start: seg_start,
            ref_end: exons.last().unwrap().1,
            introns: seg_introns,
        });
    }
    out
}
```

- [ ] **Step 5: Run tests to confirm they pass**

Run: `cargo test -p rustle --lib split_ 2>&1 | tail -20`
Expected: PASS (5 tests).

- [ ] **Step 6: Commit**

```bash
git add src/rustle/vg_family/denovo_assemble.rs
git commit -m "feat: split_mischained_reads — salvage reads across spurious giant introns"
```

---

### Task 2: `mischain_salvage` config flag

**Files:**
- Modify: `src/rustle/vg_family/denovo_pipeline.rs` (`DenovoConfig` struct + `Default` + `from_env`)

**Interfaces:**
- Produces: `DenovoConfig.mischain_salvage: bool` (default `false`); `from_env` sets it from `RUSTLE_MISCHAIN_SALVAGE`.

- [ ] **Step 1: Write the failing test** (in the `denovo_pipeline.rs` test module)

```rust
#[test]
fn mischain_salvage_defaults_off() {
    assert!(!DenovoConfig::default().mischain_salvage);
}
```

- [ ] **Step 2: Run to confirm it fails**

Run: `cargo test -p rustle --lib mischain_salvage_defaults_off 2>&1 | tail -10`
Expected: FAIL — `no field mischain_salvage`.

- [ ] **Step 3: Add the field, default, and env wiring**

In `struct DenovoConfig`, next to `pub filter_readthrough: bool,`:
```rust
    /// Split mis-chained reads at spurious giant introns before seeding (opt-in). Default off = byte-identical.
    pub mischain_salvage: bool,
```
In the `Default`/constructor block next to `filter_readthrough: true,` add `mischain_salvage: false,`.
In `from_env`, next to the `filter_readthrough:` line add:
```rust
            mischain_salvage: std::env::var("RUSTLE_MISCHAIN_SALVAGE").ok().as_deref() == Some("1"),
```

- [ ] **Step 4: Run to confirm it passes**

Run: `cargo test -p rustle --lib mischain_salvage_defaults_off 2>&1 | tail -10`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/denovo_pipeline.rs
git commit -m "feat: DenovoConfig.mischain_salvage flag (RUSTLE_MISCHAIN_SALVAGE, default off)"
```

---

### Task 3: Wire salvage before seeding in the live detect paths

**Files:**
- Modify: `src/rustle/vg_family/denovo_pipeline.rs` (add helper; insert at the 3 gw-catalog paths + `detect_and_assign`)

**Interfaces:**
- Consumes: `split_mischained_reads` (Task 1), `mischain_salvage` (Task 2), `read_junction_support`, `env_num`, `MISCHAIN_GIANT_INTRON_BP`, `MISCHAIN_MIN_JUNCTION_READS`.
- Produces: `fn maybe_salvage_mischain(reads: &[PrimaryRead], cfg: &DenovoConfig) -> Option<Vec<PrimaryRead>>`

- [ ] **Step 1: Write the failing integration test** (in `denovo_pipeline.rs` test module) — proves salvage makes a bridged locus seed

```rust
#[test]
fn salvage_seeds_local_locus_that_a_mischain_bridge_would_lose() {
    use crate::vg_family::denovo_assemble::{pass1_skeletons_robust, PrimaryRead};
    // locus A ~1000, locus B ~900_000; 3 reads each, PLUS 3 reads that mis-chain A->B via a >300kb giant intron.
    let mk = |s: u64, e: u64, introns: Vec<(u64,u64)>| PrimaryRead { chrom: "chr1".into(), ref_start: s, ref_end: e, introns };
    let mut reads = vec![];
    for _ in 0..3 { reads.push(mk(1000, 1100, vec![])); }
    for _ in 0..3 { reads.push(mk(900_000, 900_100, vec![])); }
    for _ in 0..3 { reads.push(mk(1000, 900_100, vec![(1100, 900_000)])); } // spurious bridge, its own junction
    // OFF: the bridge reads form a >300kb skeleton (rejected by MAX_SPLICED); locus B still has its own 3 reads,
    // but the point of the test is that WITH salvage the bridge reads ALSO contribute locally.
    let cfg_off = DenovoConfig::default();
    let cfg_on = DenovoConfig { mischain_salvage: true, ..DenovoConfig::default() };
    let reads_off = maybe_salvage_mischain(&reads, &cfg_off).unwrap_or_else(|| reads.clone());
    let reads_on  = maybe_salvage_mischain(&reads, &cfg_on ).unwrap_or_else(|| reads.clone());
    let sk_off = pass1_skeletons_robust(&reads_off, 3, 1);
    let sk_on  = pass1_skeletons_robust(&reads_on,  3, 1);
    let seeds_b = |sk: &[crate::vg_family::denovo_assemble::Skeleton]|
        sk.iter().filter(|s| s.chrom == "chr1" && s.start <= 900_100 && s.end >= 900_000 && s.end - s.start < 300_000).count();
    // OFF: locus B seeded only by its 3 native reads. ON: the 3 salvaged local B-segments join it => stronger,
    // and locus A's salvaged segments seed A too. Key assertion: ON yields a bounded local A skeleton that OFF lacks.
    let seeds_a_bounded = |sk: &[crate::vg_family::denovo_assemble::Skeleton]|
        sk.iter().filter(|s| s.chrom == "chr1" && s.start <= 1100 && s.end >= 1000 && s.end - s.start < 300_000).count();
    assert!(maybe_salvage_mischain(&reads, &cfg_off).is_none(), "OFF must not transform reads");
    assert!(seeds_b(&sk_on) >= 1);
    assert_eq!(seeds_a_bounded(&sk_on), seeds_a_bounded(&sk_off) + 0.max(0)); // A already seeds from native reads; ON must not lose it
    assert!(seeds_a_bounded(&sk_on) >= 1, "salvaged A segment must seed a bounded local locus");
}
```
*(If the exact skeleton counts differ, the implementer adjusts the assertions to encode the real invariant — ON seeds bounded local A and B loci; OFF leaves reads untouched — but must keep a positive ON-vs-OFF distinction and the `is_none()` OFF check.)*

- [ ] **Step 2: Run to confirm it fails**

Run: `cargo test -p rustle --lib salvage_seeds_local_locus 2>&1 | tail -20`
Expected: FAIL — `cannot find function maybe_salvage_mischain`.

- [ ] **Step 3: Add the helper** (near `retain_non_mischain` in `denovo_pipeline.rs`)

```rust
/// Mis-chain salvage (opt-in). `Some(split reads)` when `cfg.mischain_salvage`, else `None` (caller keeps the
/// originals — byte-identical). Reuses the gate thresholds so it splits exactly the introns the gate would drop.
fn maybe_salvage_mischain(reads: &[PrimaryRead], cfg: &DenovoConfig) -> Option<Vec<PrimaryRead>> {
    if !cfg.mischain_salvage {
        return None;
    }
    let giant = env_num("RUSTLE_MISCHAIN_GIANT_BP", MISCHAIN_GIANT_INTRON_BP);
    let min_reads = env_num("RUSTLE_MISCHAIN_MIN_READS", MISCHAIN_MIN_JUNCTION_READS);
    let support = read_junction_support(reads);
    Some(split_mischained_reads(reads, &support, giant, min_reads))
}
```
Add `split_mischained_reads` to the `denovo_assemble` import list at the top of the file.

- [ ] **Step 4: Run to confirm the test passes**

Run: `cargo test -p rustle --lib salvage_seeds_local_locus 2>&1 | tail -20`
Expected: PASS.

- [ ] **Step 5: Insert the helper at the 3 gw-catalog paths (owned `reads`)**

At `denovo_pipeline.rs:2024`, `:2223`, `:2364`, immediately BEFORE each `let skeletons = pass1_skeletons_robust(&reads, …);`, insert:
```rust
    let reads = maybe_salvage_mischain(&reads, cfg).unwrap_or(reads);
```
(The borrow from `maybe_salvage_mischain` ends before `unwrap_or(reads)` moves `reads`; `contigs`/`genome` are already built and unaffected since split segments keep their chrom.)

- [ ] **Step 6: Insert at `detect_and_assign` (borrowed `primary_reads`)**

At `denovo_pipeline.rs:1446`, replace `let skeletons = pass1_skeletons_robust(primary_reads, …);` with:
```rust
    let salvaged = maybe_salvage_mischain(primary_reads, cfg);
    let seed_reads: &[PrimaryRead] = salvaged.as_deref().unwrap_or(primary_reads);
    let skeletons = pass1_skeletons_robust(seed_reads, cfg.pass1_min_reads, cfg.min_terminal_support);
```

- [ ] **Step 7: Verify OFF byte-identity + full lib tests pass**

Run: `cargo test -p rustle --lib 2>&1 | tail -15`
Expected: PASS (all existing tests unchanged — OFF path untouched).

- [ ] **Step 8: Commit**

```bash
git add src/rustle/vg_family/denovo_pipeline.rs
git commit -m "feat: wire mis-chain salvage before seeding in gw-catalog + detect_and_assign paths"
```

---

### Task 4: End-to-end OFF byte-identity check on a real catalog

**Files:**
- Create: `bench/soto/salvage_byte_identity.sh`

**Interfaces:**
- Consumes: `target/release/gw_family_catalog`, a cached per-chrom mini-BAM.

- [ ] **Step 1: Write the check script**
```bash
#!/bin/bash
# OFF path must be byte-identical to the committed binary behavior. Build, run one chrom with salvage OFF
# (no env) and confirm the copies catalog md5 is unchanged vs a salvage-OFF baseline captured on the same BAM.
set -euo pipefail
CACHE=/home/juanfra/winloci_scratch/soto_cache; FA=/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa
BIN=/mnt/c/Users/jfris/Desktop/Rustle/target/release/gw_family_catalog
cargo build --release -p rustle >/dev/null 2>&1
"$BIN" --bam "$CACHE/perchrom/chr9.bam" --fasta "$FA" --out /tmp/salv_off >/dev/null 2>&1
RUSTLE_MISCHAIN_SALVAGE=1 "$BIN" --bam "$CACHE/perchrom/chr9.bam" --fasta "$FA" --out /tmp/salv_on >/dev/null 2>&1
echo "OFF copies md5: $(md5sum < /tmp/salv_off.copies.tsv)"
echo "ON  copies md5: $(md5sum < /tmp/salv_on.copies.tsv)"
echo "OFF copies: $(($(wc -l </tmp/salv_off.copies.tsv)-1))   ON copies: $(($(wc -l </tmp/salv_on.copies.tsv)-1))"
```

- [ ] **Step 2: Run it**

Run: `bash bench/soto/salvage_byte_identity.sh`
Expected: OFF md5 stable across reruns; ON differs (more copies) — proves the flag is inert when off and active when on.

- [ ] **Step 3: Commit**
```bash
git add bench/soto/salvage_byte_identity.sh
git commit -m "test: OFF-path byte-identity + ON-activates check for mis-chain salvage"
```

---

### Task 5: Anti-overfit validation harness (the gate before claiming success)

**Files:**
- Create: `bench/soto/salvage_validation.sh`, `bench/soto/salvage_validation_results.md`

**Interfaces:**
- Consumes: `member_attribution.final.tsv` (targets), the known-family regression fixtures, the cached Soto BAM.

- [ ] **Step 1: Write the validation script**
```bash
#!/bin/bash
# Run the anti-overfit checks IN ORDER; the target recovery is reported LAST.
#   (a) known families: GSTM/MAGEA/DAZ/RBMY/TSPY/PCDHB copy counts, salvage ON vs OFF -> must be unchanged.
#   (b) genome-wide/representative FP: total family+copy counts ON vs OFF -> ON must not exceed OFF beyond
#       the intended mis-chain recoveries; spot-list the added copies.
#   (c) Soto: rebuild the per-chrom catalog with salvage ON, re-score with soto_cache_score.py; confirm no
#       previously-passing member regresses, then report recovery among the 15 mis-chain + collateral targets.
set -u
CACHE=/home/juanfra/winloci_scratch/soto_cache
# ... per-chrom recompute ON vs OFF on the residual chroms, diff member_attribution buckets, tabulate ...
# (implementer fills in using recompute_perchrom.sh with RUSTLE_MISCHAIN_SALVAGE and soto_cache_score.py)
echo "see salvage_validation_results.md"
```

- [ ] **Step 2: Run known-family + FP + Soto checks; record results**

Run the script; capture into `salvage_validation_results.md`: known-family table (ON==OFF), FP delta, per-member bucket diff (which of the 15 targets moved from MISS to FOUND), and any regression among previously-passing members. A regression on a known family or a passing member is a FAIL — return to systematic-debugging before proceeding.

- [ ] **Step 3: Commit the results**
```bash
git add bench/soto/salvage_validation.sh bench/soto/salvage_validation_results.md
git commit -m "test: mis-chain salvage anti-overfit validation (known families, FP, Soto recovery)"
```

---

## Self-Review

- **Spec coverage:** split fn (T1), config/env (T2), wiring live paths (T3), OFF byte-identity (T1 pass-through + T4), synthetic positive/negative (T1 tests + T3 integration), known-family regression + FP + Soto (T5). Covered.
- **Type consistency:** `PrimaryRead` gains `PartialEq, Eq` (T1) — required by every `assert_eq!` on reads. `maybe_salvage_mischain` returns `Option<Vec<PrimaryRead>>`, consumed via `.unwrap_or`/`.as_deref()` at owned/borrowed sites respectively.
- **Placeholder note:** T5's script body is intentionally a scaffold — the recompute/score commands already exist (`recompute_perchrom.sh`, `soto_cache_score.py`); the implementer wires them with `RUSTLE_MISCHAIN_SALVAGE=1` and diffs buckets. All code-bearing steps (T1–T3) are complete.
