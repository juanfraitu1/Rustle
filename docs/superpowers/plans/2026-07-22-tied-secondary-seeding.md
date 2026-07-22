# Tied-secondary (multi-mapper) seeding — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Seed a candidate locus from AS-tied secondary reads that share an intron chain when no primary anchors it, to recover covered-but-tied K=0 Soto members as detected-but-unassignable loci.

**Architecture:** New pure function `tied_seed_skeletons` groups the already-AS-gated tied secondaries by `(chrom, intron-chain)`, keeps chains with ≥ `min_reads` whose span doesn't overlap a primary skeleton, and returns them as `Skeleton`s. It's called inside `detect_and_assign` — the production detection path, which already receives the tied secondaries as its `rescue_extra: &[PrimaryRead]` param — behind a `cfg.tied_seed` flag (off = byte-identical). A `--tied-seed` CLI flag on `copy_assign` sets it. Phase 1 measures recall gain vs over-seeding on the Soto benchmark by diffing runs with/without the flag.

**Tech Stack:** Rust (cargo), the `rustle` crate; Python3 + samtools for the benchmark eval.

## Global Constraints

- **Integration point is `detect_and_assign`, NOT `detect_families`** (spec said `detect_families`; exploration showed that is a test-only function while `detect_and_assign` is the production path and already carries the tied secondaries in `rescue_extra`). This is a strict improvement over the spec.
- **Off = byte-identical:** all new behavior is behind `cfg.tied_seed` (default `false`). Existing tests must stay green unchanged.
- **Seed only from SPLICED tied reads** (non-empty intron chain) — the shared-intron-chain agreement gate. Extent = simple min `ref_start` / max `ref_end` over the chain group.
- **Phase 1 measures via with/without-`--tied-seed` diff.** Propagating the `tied_seeded` flag all the way to the catalog's per-member status string ("detected (tied), unassignable") is deferred to Phase 2; the flag is set on `Skeleton` now for that future use.
- **CRASH RULE for the benchmark:** `copy_assign` runs FOREGROUND, serial, small batches, outputs to `/home/juanfra/winloci_scratch` (NOT `/tmp`). No nohup/waiter/pkill.
- Determinism: group with `BTreeMap`, mirroring `pass1_skeletons_robust`.

---

### Task 1: `Skeleton.tied_seeded` field

**Files:**
- Modify: `src/rustle/vg_family/denovo_assemble.rs` (struct at :48; literals at :104, :844, :863, :1077)

**Interfaces:**
- Produces: `Skeleton` gains `pub tied_seeded: bool`. Primary/unspliced seeding sets it `false`.

- [ ] **Step 1: Add the field to the struct** (`denovo_assemble.rs:48`)

```rust
pub struct Skeleton {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub n_reads: u32,
    pub introns: Vec<(u64, u64)>,
    pub tied_seeded: bool,
}
```

- [ ] **Step 2: Set `tied_seeded: false` in all four existing literals.**

`denovo_assemble.rs:104` (inside `pass1_skeletons_robust`):
```rust
            Skeleton {
                chrom: chrom.to_string(),
                start: starts[si],
                end: ends[ei],
                n_reads: n,
                introns,
                tied_seeded: false,
            }
```
`denovo_assemble.rs:844` (inside `cluster_unspliced` — the `out.push(Skeleton {` literal): add `tied_seeded: false,` before the closing `}`.
`denovo_assemble.rs:863` (the one-line helper): `Skeleton { chrom: chrom.into(), start, end, n_reads, introns: introns.to_vec(), tied_seeded: false }`.
`denovo_assemble.rs:1077` (test helper): `Skeleton { chrom: chrom.into(), start, end, n_reads: n, introns: introns.to_vec(), tied_seeded: false }`.

- [ ] **Step 3: Build and run existing tests to confirm no regression.**

Run: `cargo test -p rustle --lib vg_family::denovo_assemble 2>&1 | tail -20`
Expected: compiles; all existing `denovo_assemble` tests PASS.

- [ ] **Step 4: Commit.**

```bash
git add src/rustle/vg_family/denovo_assemble.rs
git commit -m "feat(seed): add Skeleton.tied_seeded provenance field (default false)"
```

---

### Task 2: `tied_seed_skeletons` (the core, TDD)

**Files:**
- Modify: `src/rustle/vg_family/denovo_assemble.rs` (add function near `pass1_skeletons_robust`; tests in the `mod tests` at :946, using the existing `pr(chrom,start,end,introns)` helper at :949)

**Interfaces:**
- Consumes: `PrimaryRead { chrom: String, introns: Vec<(u64,u64)>, ref_start: u64, ref_end: u64, .. }`, `Skeleton` (with `tied_seeded`).
- Produces: `pub fn tied_seed_skeletons(tied_reads: &[PrimaryRead], primary_skeletons: &[Skeleton], min_reads: u32) -> Vec<Skeleton>` — every returned `Skeleton` has `tied_seeded == true`.

- [ ] **Step 1: Write the failing tests** (append inside `mod tests`, `denovo_assemble.rs`)

```rust
    #[test]
    fn tied_seed_agreeing_chain_seeds_one() {
        // 3 tied secondaries sharing one intron chain at a locus with no primary skeleton.
        let tied = vec![
            pr("chr1", 100, 500, &[(200, 300)]),
            pr("chr1", 110, 490, &[(200, 300)]),
            pr("chr1", 105, 495, &[(200, 300)]),
        ];
        let primaries: Vec<Skeleton> = vec![]; // nothing seeded here yet
        let out = tied_seed_skeletons(&tied, &primaries, 3);
        assert_eq!(out.len(), 1);
        assert!(out[0].tied_seeded);
        assert_eq!(out[0].introns, vec![(200, 300)]);
        assert_eq!(out[0].start, 100); // min ref_start
        assert_eq!(out[0].end, 500); // max ref_end
        assert_eq!(out[0].n_reads, 3);
    }

    #[test]
    fn tied_seed_scattered_chains_seed_nothing() {
        // 3 reads, 3 different intron chains -> no chain reaches min_reads.
        let tied = vec![
            pr("chr1", 100, 500, &[(200, 300)]),
            pr("chr1", 100, 500, &[(210, 300)]),
            pr("chr1", 100, 500, &[(200, 310)]),
        ];
        assert_eq!(tied_seed_skeletons(&tied, &[], 3).len(), 0);
    }

    #[test]
    fn tied_seed_below_min_reads_seeds_nothing() {
        let tied = vec![
            pr("chr1", 100, 500, &[(200, 300)]),
            pr("chr1", 100, 500, &[(200, 300)]),
        ];
        assert_eq!(tied_seed_skeletons(&tied, &[], 3).len(), 0);
    }

    #[test]
    fn tied_seed_dedups_against_overlapping_primary() {
        // Same locus is already a primary skeleton -> the tied group must NOT re-seed it.
        let tied = vec![
            pr("chr1", 100, 500, &[(200, 300)]),
            pr("chr1", 100, 500, &[(200, 300)]),
            pr("chr1", 100, 500, &[(200, 300)]),
        ];
        let primaries = vec![Skeleton {
            chrom: "chr1".into(), start: 90, end: 480, n_reads: 5,
            introns: vec![(200, 300)], tied_seeded: false,
        }];
        assert_eq!(tied_seed_skeletons(&tied, &primaries, 3).len(), 0);
    }
```

- [ ] **Step 2: Run to verify they fail**

Run: `cargo test -p rustle --lib tied_seed 2>&1 | tail -20`
Expected: FAIL to compile — `tied_seed_skeletons` not found.

- [ ] **Step 3: Implement** (add after `pass1_skeletons_robust`, before `primary_read_from_record`, `denovo_assemble.rs`)

```rust
/// Seed skeletons from AS-tied SECONDARY reads (already gated by [`tied_secondary_reads`]) that AGREE on an
/// intron chain, at loci not already covered by a primary skeleton. Recovers "starved" co-located copies
/// (K=0 members with 0 primaries) as DETECTED-but-unassignable loci: the copy-assignment gate still abstains
/// on them (no copy-specific PSV), so they never falsely acquire read assignments. Groups by
/// `(chrom, intron-chain)` like `pass1_skeletons_robust`, keeps chains with `>= min_reads`, and drops any
/// whose span overlaps a `primary_skeletons` entry on the same chrom (that locus is already seeded). Only
/// spliced reads (non-empty chain) seed — the shared-intron-chain agreement gate. Extent = min start / max
/// end over the group. Deterministic (`BTreeMap` order = sorted by `(chrom, introns)`).
pub fn tied_seed_skeletons(
    tied_reads: &[PrimaryRead],
    primary_skeletons: &[Skeleton],
    min_reads: u32,
) -> Vec<Skeleton> {
    use std::collections::BTreeMap;
    let mut groups: BTreeMap<(&str, Vec<(u64, u64)>), (u32, u64, u64)> = BTreeMap::new();
    for r in tied_reads {
        if r.introns.is_empty() {
            continue; // shared-intron-chain gate: only spliced reads seed
        }
        let e = groups
            .entry((r.chrom.as_str(), r.introns.clone()))
            .or_insert((0, u64::MAX, 0));
        e.0 += 1;
        e.1 = e.1.min(r.ref_start);
        e.2 = e.2.max(r.ref_end);
    }
    groups
        .into_iter()
        .filter(|(_, (n, _, _))| *n >= min_reads)
        .filter_map(|((chrom, introns), (n, start, end))| {
            let overlaps_primary = primary_skeletons
                .iter()
                .any(|s| s.chrom == chrom && s.start < end && start < s.end);
            if overlaps_primary {
                return None;
            }
            Some(Skeleton {
                chrom: chrom.to_string(),
                start,
                end,
                n_reads: n,
                introns,
                tied_seeded: true,
            })
        })
        .collect()
}
```

- [ ] **Step 4: Run to verify they pass**

Run: `cargo test -p rustle --lib tied_seed 2>&1 | tail -20`
Expected: `test result: ok. 4 passed`.

- [ ] **Step 5: Commit.**

```bash
git add src/rustle/vg_family/denovo_assemble.rs
git commit -m "feat(seed): tied_seed_skeletons — seed loci from agreeing AS-tied multimappers"
```

---

### Task 3: `DenovoConfig.tied_seed` flag + wire into `detect_and_assign`

**Files:**
- Modify: `src/rustle/vg_family/denovo_pipeline.rs` (`DenovoConfig` struct + `Default` impl at :120-141; `detect_and_assign` seeding at ~:1396)

**Interfaces:**
- Consumes: `tied_seed_skeletons(&[PrimaryRead], &[Skeleton], u32) -> Vec<Skeleton>` (Task 2); `detect_and_assign`'s existing `rescue_extra: &[PrimaryRead]` param.
- Produces: `DenovoConfig.tied_seed: bool` (default `false`).

- [ ] **Step 1: Write the failing test** (append inside `denovo_pipeline.rs` `mod tests`)

```rust
    #[test]
    fn denovo_config_tied_seed_defaults_off() {
        assert!(!DenovoConfig::default().tied_seed);
    }
```

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test -p rustle --lib denovo_config_tied_seed 2>&1 | tail -10`
Expected: FAIL to compile — no field `tied_seed`.

- [ ] **Step 3: Add the field + default** (`denovo_pipeline.rs`)

In the `DenovoConfig` struct (near the other `pub bool` fields, e.g. after `pub dna_family_fallback: bool,`):
```rust
    pub tied_seed: bool,
```
In `impl Default for DenovoConfig` (after `dna_family_fallback: false,`):
```rust
            tied_seed: false,
```
(If `from_env` also constructs the struct field-by-field, add `tied_seed: false,` there too; if it starts from `Default::default()`, no change needed.)

- [ ] **Step 4: Import + wire the seeding call** into `detect_and_assign`. Ensure `tied_seed_skeletons` is in the `use ...denovo_assemble::{...}` import list at the top of `denovo_pipeline.rs`, then change the seeding line (~:1396):

```rust
    let mut skeletons = pass1_skeletons_robust(primary_reads, cfg.pass1_min_reads, cfg.min_terminal_support);
    if cfg.tied_seed {
        let tied = tied_seed_skeletons(rescue_extra, &skeletons, cfg.pass1_min_reads);
        skeletons.extend(tied);
    }
```
(Note: `skeletons` was previously `let skeletons = ...`; make it `let mut skeletons`. `rescue_extra` is the existing param.)

- [ ] **Step 5: Run the new test + the full detect_and_assign suite (off-path byte-identity proof)**

Run: `cargo test -p rustle --lib denovo_pipeline 2>&1 | tail -20`
Expected: `denovo_config_tied_seed_defaults_off` PASSES; all existing `detect_and_assign_*` tests still PASS (they run with `tied_seed=false` → unchanged skeletons).

- [ ] **Step 6: Commit.**

```bash
git add src/rustle/vg_family/denovo_pipeline.rs
git commit -m "feat(seed): DenovoConfig.tied_seed flag wired into detect_and_assign (off=identical)"
```

---

### Task 4: `--tied-seed` CLI flag on `copy_assign`

**Files:**
- Modify: `src/bin/copy_assign.rs` (`Args` struct ~:88-200; the `cfg` build; the `extra` fetch at ~:968)

**Interfaces:**
- Consumes: `DenovoConfig.tied_seed` (Task 3); the existing `--recover-copies` / `--as-ratio` machinery and `tied_secondary_reads_in_region`.
- Produces: a `--tied-seed` flag that sets `cfg.tied_seed = true` AND ensures the tied secondaries are fetched into `extra`.

- [ ] **Step 1: Add the CLI arg** to `Args` (`copy_assign.rs`, next to the other `#[arg(long, default_value_t = false)]` bools):

```rust
    /// Seed candidate loci from AS-tied secondary reads that share an intron chain, even with no primary
    /// (recovers covered-but-tied K=0 copies as detected-but-unassignable). Implies fetching tied secondaries.
    #[arg(long, default_value_t = false)]
    tied_seed: bool,
```

- [ ] **Step 2: Set the config flag** where `cfg` is built (after `let mut cfg = ...`, near the other `cfg.<flag> = args.<flag>;` lines):

```rust
    cfg.tied_seed = args.tied_seed;
```

- [ ] **Step 3: Ensure tied secondaries are fetched when `--tied-seed`** — change the `extra` guard (~:968) from `if args.recover_copies` to:

```rust
        let extra = if args.recover_copies || args.tied_seed {
            tied_secondary_reads_in_region(&args.bam, contig, lo, hi, args.as_ratio).unwrap_or_default()
        } else {
            Vec::new()
        };
```

- [ ] **Step 4: Build the binary.**

Run: `cargo build --release --bin copy_assign 2>&1 | tail -5`
Expected: compiles; `target/release/copy_assign` updated.

- [ ] **Step 5: Smoke-test the flag is recognized.**

Run: `target/release/copy_assign --help 2>&1 | grep -A1 tied-seed`
Expected: shows `--tied-seed` in the help.

- [ ] **Step 6: Commit.**

```bash
git add src/bin/copy_assign.rs
git commit -m "feat(seed): --tied-seed CLI flag on copy_assign"
```

---

### Task 5: Benchmark evaluation — recall gain vs over-seeding

**Files:**
- Create: `bench/soto/tied_seed_eval.py` (runner + classifier)
- Create (output): `bench/soto/tied_seed_eval.tsv`

**Interfaces:**
- Consumes: `target/release/copy_assign` with/without `--tied-seed`; the Soto member truth (`bench/soto/soto_member_detection.tsv`), the scoped Soto BAM (`/mnt/linuxdisk/home/juanfraitu/winloci_data/soto_reads.bam`), the genome (`chm13v2.0.fa`), and the family region set derived from `bench/soto/80_fams.chr.bed`.
- Produces: `tied_seed_eval.tsv` with the recall-before/after and per-new-locus TP/FP.

- [ ] **Step 1: Write the eval driver** `bench/soto/tied_seed_eval.py`:

```python
#!/usr/bin/env python3
"""Phase-1 benchmark for tied-secondary seeding. Run copy_assign on the Soto family regions WITH and
WITHOUT --tied-seed; the loci detected only WITH it are the tied-seeded ones. Classify each against the
Soto member truth: TP = overlaps a previously-MISSED member; FP = overlaps no member (over-seeding).
CRASH RULE: copy_assign runs FOREGROUND, serial, output to winloci_scratch."""
import csv, subprocess, os
CA = "target/release/copy_assign"
D = "/mnt/linuxdisk/home/juanfraitu/winloci_data"
BAM = f"{D}/soto_reads.bam"; FASTA = f"{D}/chm13v2.0.fa"
SCR = "/home/juanfra/winloci_scratch"
BED = "bench/soto/80_fams.chr.bed"
MEMB = "bench/soto/soto_member_detection.tsv"

def members():
    out = []
    for r in csv.DictReader(open(MEMB), delimiter="\t"):
        out.append((r["chrom"], int(r["start"]), int(r["end"]),
                    r["family_id"], r["gene"], r["detected"].strip().upper() == "Y"))
    return out

def run(tag, tied):
    out = f"{SCR}/tied_{tag}"
    cmd = [CA, "--bam", BAM, "--fasta", FASTA, "--regions", BED,
           "--recover-copies", "--out", out]
    if tied:
        cmd.append("--tied-seed")
    subprocess.run(cmd, check=True)
    # copies TSV: chrom,start,end per detected copy (adapt columns to copy_assign's <out>.copies.tsv header)
    loci = []
    cp = f"{out}.copies.tsv"
    if os.path.exists(cp):
        for r in csv.DictReader(open(cp), delimiter="\t"):
            loci.append((r["chrom"], int(r["start"]), int(r["end"])))
    return set(loci)

def overlaps(a, b):
    return a[0] == b[0] and a[1] < b[2] and b[1] < a[2]

if __name__ == "__main__":
    mem = members()
    missed = [m for m in mem if not m[5]]
    base = run("off", False)
    withts = run("on", True)
    new = withts - base                      # loci detected only with --tied-seed
    rows = []
    tp = fp = 0
    for locus in sorted(new):
        hit = next((m for m in mem if overlaps(locus, m)), None)
        miss_hit = next((m for m in missed if overlaps(locus, m)), None)
        cls = "TP-recovers-missed-member" if miss_hit else ("overlaps-detected-member" if hit else "FP-no-member")
        if miss_hit: tp += 1
        elif not hit: fp += 1
        rows.append((locus[0], locus[1], locus[2], cls,
                     f"{miss_hit[3]}:{miss_hit[4]}" if miss_hit else (f"{hit[3]}:{hit[4]}" if hit else "")))
    with open("bench/soto/tied_seed_eval.tsv", "w") as f:
        f.write("chrom\tstart\tend\tclass\tmember\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    det0 = sum(1 for m in mem if m[5])
    print(f"baseline detected members: {det0}/{len(mem)} = {100*det0/len(mem):.1f}%")
    print(f"tied-seed NEW loci: {len(new)} | TP (recover missed member): {tp} | FP (no member): {fp}")
    print(f"projected recall: {det0}+{tp} = {det0+tp}/{len(mem)} = {100*(det0+tp)/len(mem):.1f}%")
```

- [ ] **Step 2: Confirm `copy_assign`'s region + output flags** match the script. Run:

Run: `target/release/copy_assign --help 2>&1 | grep -iE "regions|--out|copies|--bam|--fasta"`
Expected: shows the real flag names. If `--regions`/`--out` differ, fix the `cmd` list and the `<out>.copies.tsv` filename/columns in Step 1 to match before running.

- [ ] **Step 3: Run the benchmark (FOREGROUND, serial — crash rule).**

Run: `python3 bench/soto/tied_seed_eval.py`
Expected: prints baseline recall (76.2%), the count of new tied-seeded loci, TP/FP split, and projected recall. Writes `bench/soto/tied_seed_eval.tsv`.

- [ ] **Step 4: Read the result and decide.** Inspect `bench/soto/tied_seed_eval.tsv`. Success = recall rises with FP small (report the actual TP/FP). Note the number for the Phase-2 go/no-go.

- [ ] **Step 5: Commit the eval + result.**

```bash
git add bench/soto/tied_seed_eval.py bench/soto/tied_seed_eval.tsv
git commit -m "bench(seed): Phase-1 tied-seed eval — recall gain vs over-seeding on Soto"
```

---

## Self-Review

- **Spec coverage:** `tied_seed_skeletons` (Task 2) ✓; `Skeleton.tied_seeded` (Task 1) ✓; `DenovoConfig.tied_seed` + integration (Task 3, into `detect_and_assign` per the Global Constraints deviation) ✓; `--tied-seed` CLI (Task 4) ✓; three-gate phantom guard — AS-tie (upstream `tied_secondary_reads`) + shared intron-chain (Task 2 grouping) + distinct-locus dedup (Task 2 filter) ✓; benchmark TP/FP (Task 5) ✓; off = byte-identical (Task 3 Step 5) ✓. Deferred per spec "out of scope": Phase 2 genome-wide, assignment changes, full status-string propagation.
- **Placeholder scan:** none — all code and commands are literal. Task 5 Step 2 explicitly reconciles the `copy_assign` flag names before the run (the one spot that depends on the binary's real CLI).
- **Type consistency:** `tied_seed_skeletons(&[PrimaryRead], &[Skeleton], u32) -> Vec<Skeleton>` identical in Tasks 2 and 3; `Skeleton.tied_seeded: bool` and `DenovoConfig.tied_seed: bool` consistent across tasks.
