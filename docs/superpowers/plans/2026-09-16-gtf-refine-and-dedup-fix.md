# Read-Dedup Fix + `--gtf-refine` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix `copy_assign`'s within-window primary-read dedup bug (measuring its O2/de novo impact first), add an opt-in annotation-free `--gtf-refine` bundle for `--gtf` (five components; `tss` added by the 2026-09-16 addendum, Tasks 7–9), reproduce the chr20 simulation numbers exactly, and validate on held-out chr17 under pre-registered rules.

**Architecture:** The dedup fix is a bin-private pure function in `src/bin/copy_assign.rs` (cross-window dedup only, with a `RUSTLE_LEGACY_PLACEMENT_DEDUP` escape hatch). The bundle's logic lives in a new pure lib module `src/rustle/vg_family/gtf_refine.rs`, wired ONLY into the `if args.gtf` block of `copy_assign.rs` behind `--gtf-refine`. Data tasks use generalized bench scripts under `bench/` and write outputs to `/mnt/linuxdisk`.

**Tech Stack:** Rust (existing crate), bash + Python 3 bench scripts, gffcompare, SQANTI3 (conda env `sqanti3`), StringTie, FLAIR (conda env `flair`), samtools.

**Spec:** `docs/superpowers/specs/2026-09-16-gtf-refine-and-dedup-fix-design.md`

## Global Constraints

- Default outputs change ONLY through the dedup fix. `RUSTLE_LEGACY_PLACEMENT_DEDUP=1` must reproduce the pre-fix binary byte-for-byte on every output file.
- `--gtf-refine` unset ⇒ byte-identical to the fixed-dedup binary. `--gtf-refine` requires `--gtf`.
- Never modify `pass1_skeletons_robust`, `pass1_skeletons_robust_with`, `detect_and_assign`, `assemble_gate`, or `gw_family_catalog`. `--gtf-refine` touches only the `if args.gtf { ... }` block.
- No `--gtf-refine` component reads the annotation.
- Thresholds are exactly the spec's: strand margin 0.90; subset overhang ≤ 5 bp only at a non-terminal container exon; mono coverage ≥ 50% of the model length; fragment emission `exact + fragments >= cfg.pass1_min_reads` (default 2); assign-or-abstain fragments.
- O2 (`--families`) byte-identity between fixed and legacy dedup is a HARD STOP condition (Task 2).
- `tss` (addendum): 5' end = densest window of exact-chain read 5' ends, ties to the most upstream, multi-exon stranded models only, applied after `mono`; its window `TSS_WINDOW_BP` is chosen once on chr20 (Task 7) and frozen before the pre-registration.
- No chr17 file may be generated before `docs/PREREG_gtf_refine_chr17_2026-09-16.md` is committed (Task 10). No threshold/rule change after any chr17 number is seen.
- WSL2 rules: build with `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target`, `--release`; redirect cargo/tool output to a file, then tail it; ONE heavy command at a time in the foreground; never background, never `pkill -f`; big outputs under `/mnt/linuxdisk/home/juanfraitu/`.

---

### Task 1: Dedup fix (cross-window only) with escape hatch

**Files:**
- Modify: `src/bin/copy_assign.rs` (the `let (primary, mut bam_reads) = { ... }` block, currently ~lines 2392-2427; the `use rustle::vg_family::denovo_assemble::{...}` list at ~line 28; the `#[cfg(test)] mod tests` block at the end of the file)

**Interfaces:**
- Produces: `fn dedup_primary_across_windows(per_window: Vec<Vec<PrimaryRead>>, legacy_placement_dedup: bool) -> Vec<PrimaryRead>` (bin-private); env var `RUSTLE_LEGACY_PLACEMENT_DEDUP`; and a preserved pre-fix binary at `/mnt/linuxdisk/home/juanfraitu/dedupfix/copy_assign.base` for Task 2.

- [ ] **Step 1: Preserve the pre-fix binary (before ANY source edit)**

```bash
mkdir -p /mnt/linuxdisk/home/juanfraitu/dedupfix
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo build --release --bin copy_assign > /tmp/t1_base_build.log 2>&1; echo "EXIT=$?"; tail -5 /tmp/t1_base_build.log
cp /mnt/linuxdisk/home/juanfraitu/rustle_target/release/copy_assign /mnt/linuxdisk/home/juanfraitu/dedupfix/copy_assign.base
git rev-parse HEAD > /mnt/linuxdisk/home/juanfraitu/dedupfix/copy_assign.base.commit
```
Expected: `EXIT=0`, the copied binary exists, and the commit file names the base commit.

- [ ] **Step 2: Write the failing unit tests** (append inside the existing `#[cfg(test)] mod tests` in `src/bin/copy_assign.rs`, which already has `use super::*;`)

```rust
    fn pr_read(chrom: &str, s: u64, e: u64, introns: Vec<(u64, u64)>) -> PrimaryRead {
        PrimaryRead { chrom: chrom.into(), ref_start: s, ref_end: e, introns, reverse: false }
    }

    #[test]
    fn dedup_keeps_distinct_molecules_with_identical_placement_in_one_window() {
        let w = vec![vec![
            pr_read("c1", 100, 500, vec![(200, 300)]),
            pr_read("c1", 100, 500, vec![(200, 300)]),
            pr_read("c1", 120, 500, vec![(200, 300)]),
        ]];
        assert_eq!(dedup_primary_across_windows(w.clone(), false).len(), 3, "distinct molecules must all count");
        assert_eq!(dedup_primary_across_windows(w, true).len(), 2, "legacy placement dedup collapses identical placements");
    }

    #[test]
    fn dedup_counts_each_boundary_spanning_molecule_once() {
        let a = pr_read("c1", 900, 1100, vec![]);
        let w = vec![
            vec![a.clone(), a.clone(), pr_read("c1", 100, 200, vec![])],
            vec![a.clone(), a.clone(), pr_read("c1", 1200, 1300, vec![])],
        ];
        // two distinct molecules share a placement across the boundary: each window returns both.
        assert_eq!(dedup_primary_across_windows(w, false).len(), 4);
    }

    #[test]
    fn dedup_never_merges_across_chromosomes() {
        let w = vec![vec![pr_read("c1", 100, 200, vec![])], vec![pr_read("c2", 100, 200, vec![])]];
        assert_eq!(dedup_primary_across_windows(w, false).len(), 2);
    }
```

- [ ] **Step 3: Run tests, confirm they fail**

Run: `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo test --release --bin copy_assign dedup_ > /tmp/t1_test.log 2>&1; tail -30 /tmp/t1_test.log`
Expected: compile error (`dedup_primary_across_windows` / `PrimaryRead` not in scope).

- [ ] **Step 4: Implement**

Add `PrimaryRead` to the existing import list:
```rust
use rustle::vg_family::denovo_assemble::{
    assemble_gate, merge_fuzzy_skeletons, pass1_skeletons, reads_in_region, tied_secondary_reads_in_region,
    BamIndexCache, BamRead, PrimaryRead, GATE_MIN_READS,
};
```

Add the function (near the other bin-private helpers, e.g. just above `fn discover_copies_for_family`):
```rust
/// Pool a region's per-window primary reads. A read spanning a window boundary is returned by every window it
/// overlaps and must count once, so dedup runs ACROSS windows only: a read is dropped iff an EARLIER window
/// already contributed the same `(chrom, ref_start, ref_end, introns)`. Identical placement implies the read
/// also overlaps that earlier window, so such a match is exactly a boundary duplicate, while distinct
/// molecules with identical placement are all kept. `legacy_placement_dedup`
/// (`RUSTLE_LEGACY_PLACEMENT_DEDUP=1`) reproduces the pre-2026-09-16 behaviour, which also collapsed distinct
/// molecules inside one window (human chr20, one window: 25,341 primaries -> 11,755).
fn dedup_primary_across_windows(per_window: Vec<Vec<PrimaryRead>>, legacy_placement_dedup: bool) -> Vec<PrimaryRead> {
    let mut out = Vec::new();
    let mut seen: std::collections::HashSet<(String, u64, u64, Vec<(u64, u64)>)> = std::collections::HashSet::new();
    for window in per_window {
        let mut kept_here = Vec::new();
        for x in window {
            let key = (x.chrom.clone(), x.ref_start, x.ref_end, x.introns.clone());
            if legacy_placement_dedup {
                if seen.insert(key) {
                    out.push(x);
                }
            } else if !seen.contains(&key) {
                kept_here.push(key);
                out.push(x);
            }
        }
        seen.extend(kept_here);
    }
    out
}
```

Replace the primary-read part of the gathering block. Current code (confirm it is still exactly this before editing):
```rust
        let (primary, mut bam_reads) = {
            let mut pr: Vec<_> = Vec::new();
            let mut br: Vec<_> = Vec::new();
            let mut seen = std::collections::HashSet::new();
            for (wchrom, wlo, whi) in &wins {
                ...
                // Windows are disjoint after merging, but a read spanning a boundary is returned by
                // both queries; key on (name, start) so one molecule is never two witnesses.
                for x in p {
                    // PrimaryRead has no name; its (chrom, span, intron chain) identifies the placement.
                    if seen.insert((x.chrom.clone(), x.ref_start, x.ref_end, x.introns.clone())) {
                        pr.push(x);
                    }
                }
                for x in b {
                    br.push(x);
                }
            }
```
New code (keep the window-fetch lines and everything after the loop, including the `bseen`/`br.retain(...)` block and its comment, unchanged):
```rust
        let (primary, mut bam_reads) = {
            let mut per_window_primary: Vec<Vec<PrimaryRead>> = Vec::new();
            let mut br: Vec<_> = Vec::new();
            for (wchrom, wlo, whi) in &wins {
                ...
                // A read spanning a window boundary is returned by both queries; the cross-window dedup
                // below counts it once without collapsing distinct molecules (see
                // `dedup_primary_across_windows`).
                per_window_primary.push(p);
                for x in b {
                    br.push(x);
                }
            }
            let legacy_dedup =
                matches!(std::env::var("RUSTLE_LEGACY_PLACEMENT_DEDUP"), Ok(v) if v != "0" && !v.is_empty());
            let pr = dedup_primary_across_windows(per_window_primary, legacy_dedup);
```

- [ ] **Step 5: Run tests and build, confirm pass**

```bash
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo test --release --bin copy_assign dedup_ > /tmp/t1_test.log 2>&1; tail -15 /tmp/t1_test.log
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo build --release --bin copy_assign > /tmp/t1_build.log 2>&1; echo "EXIT=$?"; tail -5 /tmp/t1_build.log
```
Expected: 3/3 pass; clean build (no new warnings).

- [ ] **Step 6: Commit**

```bash
git add src/bin/copy_assign.rs
git commit -m "fix: dedup copy_assign primary reads across windows only (RUSTLE_LEGACY_PLACEMENT_DEDUP restores old)"
```

---

### Task 2: Measure the dedup fix's impact (real data; HARD STOP condition)

**Files:**
- Create: `bench/dedup_fix_impact.sh`
- Create: `bench/DEDUP_FIX_IMPACT.md`

**Interfaces:**
- Consumes: `/mnt/linuxdisk/home/juanfraitu/dedupfix/copy_assign.base` (Task 1 Step 1); the new binary at `/mnt/linuxdisk/home/juanfraitu/rustle_target/release/copy_assign` (build it first).

- [ ] **Step 1: Write `bench/dedup_fix_impact.sh`**

```bash
#!/bin/bash
# Dedup-fix impact measurement (docs/superpowers/specs/2026-09-16-gtf-refine-and-dedup-fix-design.md, Part 1).
# Three arms per substrate: base (pre-fix binary), fixed (new binary, default), legacy (new binary,
# RUSTLE_LEGACY_PLACEMENT_DEDUP=1). Serial, foreground, outputs under /mnt/linuxdisk.
set -euo pipefail
O=/mnt/linuxdisk/home/juanfraitu/dedupfix
BASE=$O/copy_assign.base
NEW=/mnt/linuxdisk/home/juanfraitu/rustle_target/release/copy_assign
BAM=/mnt/linuxdisk/home/juanfraitu/npip_cat/npip3.bam
FA=/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta
CAT=/mnt/linuxdisk/home/juanfraitu/mec/batch.copies.tsv
CFA=/mnt/linuxdisk/home/juanfraitu/npip_cat/arm_f2/cat.copies.fa
REG=/mnt/linuxdisk/home/juanfraitu/mec/regions.txt
C20=/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20

run() { # $1=substrate $2=arm $3=binary $4=legacy(0/1) rest=args
  local sub=$1 arm=$2 bin=$3 leg=$4; shift 4
  mkdir -p "$O/$sub/$arm"
  if [ "$leg" = 1 ]; then
    ( cd "$O/$sub/$arm" && RUSTLE_LEGACY_PLACEMENT_DEDUP=1 "$bin" "$@" --out run > run.stdout 2> run.stderr )
  else
    ( cd "$O/$sub/$arm" && "$bin" "$@" --out run > run.stdout 2> run.stderr )
  fi
  echo "$sub/$arm exit=$?"
}

for arm in base fixed legacy; do
  bin=$NEW; leg=0
  [ $arm = base ] && bin=$BASE
  [ $arm = legacy ] && leg=1
  run o2_families $arm $bin $leg --bam $BAM --fasta $FA --families $CAT --copies-fa $CFA --regions $REG --dump-psv
  run denovo $arm $bin $leg --bam $BAM --fasta $FA --regions $REG
  run chr20_gtf $arm $bin $leg --gtf --bam $C20/chr20.bam --fasta $C20/chr20.fa --region chr20:1-66210255
done
```

- [ ] **Step 2: Build the new binary and run the script** (foreground; expect several minutes)

```bash
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo build --release --bin copy_assign > /tmp/t2_build.log 2>&1; echo "EXIT=$?"
bash bench/dedup_fix_impact.sh > /tmp/t2_run.log 2>&1; echo "EXIT=$?"; cat /tmp/t2_run.log
```

- [ ] **Step 3: Byte-compare**

```bash
O=/mnt/linuxdisk/home/juanfraitu/dedupfix
for sub in o2_families denovo chr20_gtf; do
  echo "=== $sub: base vs legacy (MUST be identical) ==="
  for f in $(cd $O/$sub/base && ls run.* | grep -v -e stdout -e stderr); do cmp -s $O/$sub/base/$f $O/$sub/legacy/$f && echo "same $f" || echo "DIFF $f"; done
done
echo "=== o2_families: fixed vs legacy (pre-registered: MUST be identical for assignments/quant/families/family_join) ==="
for f in run.assignments.tsv run.quant.tsv run.families.tsv run.family_join.tsv; do cmp -s $O/o2_families/fixed/$f $O/o2_families/legacy/$f && echo "same $f" || echo "DIFF $f"; done
```
If `run.params.tsv` differs base vs legacy, show `diff` and state why (it must be explainable, e.g. a recorded flag); any other base-vs-legacy difference is a Task 1 defect. **If any of the four O2 files differ fixed vs legacy: STOP. Report BLOCKED with the diffs; do not continue.**

- [ ] **Step 4: Quantify the de novo and `--gtf` deltas**

```bash
O=/mnt/linuxdisk/home/juanfraitu/dedupfix
for arm in legacy fixed; do
  echo "=== denovo/$arm ==="
  grep -m3 "primary ->" $O/denovo/$arm/run.stderr || true
  echo "families rows: $(tail -n +2 $O/denovo/$arm/run.families.tsv | wc -l)"
  echo "quant rows (copies): $(tail -n +2 $O/denovo/$arm/run.quant.tsv | wc -l)"
  head -1 $O/denovo/$arm/run.assignments.tsv
done
cd /mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20
for arm in legacy fixed; do
  gffcompare -r chr20_ref.gtf -o $O/chr20_gtf/$arm/gffc $O/chr20_gtf/$arm/run.gtf > /dev/null 2>&1
  echo "=== chr20_gtf/$arm ==="; sed -n '/Query mRNAs/p;/Transcript level/p;/Intron chain level/p;/Matching transcripts/p' $O/chr20_gtf/$arm/gffc.stats
done
```
For the de novo arm, also tabulate the assignment status column counts (find the status column in the header printed above) for legacy vs fixed, and list families present in one arm but not the other (by `family_id` and copy coordinates).
Expected for `chr20_gtf/fixed`: 1022 query mRNAs, matching transcripts 350, transcript 7.7 / 34.2, intron chain 8.1 / 42.9 (spec fidelity anchor A1). `chr20_gtf/legacy`: 976 / 347 / 7.6 / 35.6 / 8.0 / 44.6.

- [ ] **Step 5: Write `bench/DEDUP_FIX_IMPACT.md`**

Sections: the bug (code location, chr20 25,341 → 11,755); the fix; the exact commands; a table of the byte-comparison results; the O2 prediction outcome; the de novo deltas (primary counts, families, copies, status counts, differing families) with real numbers; the chr20 `--gtf` anchor outcome; and an explicit line: "Default decision pending user review before merge."

- [ ] **Step 6: Commit**

```bash
git add bench/dedup_fix_impact.sh bench/DEDUP_FIX_IMPACT.md
git commit -m "docs: measure the copy_assign dedup fix on O2, de novo and chr20 --gtf"
```

---

### Task 3: `gtf_refine.rs` — module skeleton, strict chain strand, fragment support

**Files:**
- Create: `src/rustle/vg_family/gtf_refine.rs`
- Modify: `src/rustle/vg_family/mod.rs` (append `pub mod gtf_refine;` at the END of the `pub mod` list, after `pub mod copy_discovery;`)
- Modify: `docs/MODULE_STATUS.md` (add a `gtf_refine.rs` row to the `## OPT-IN` section in the same format as neighbouring rows; bump the OPT-IN heading count and the summary-table count by 1)

**Interfaces:**
- Consumes: `crate::genome::GenomeIndex::is_canonical_junction(&self, chrom: &str, donor: u64, acceptor: u64, strand: char) -> bool`; `crate::vg_family::denovo_assemble::{PrimaryRead, Skeleton}` (Skeleton fields: `chrom, start, end, n_reads, introns, tied_seeded, read_strand, footprint, read_rev, read_tot`).
- Produces: `pub fn strict_chain_strand(genome: &GenomeIndex, chrom: &str, introns: &[(u64, u64)]) -> Option<char>`; `pub fn fragment_supported_spliced(all_chains: &[Skeleton], reads: &[PrimaryRead], min_support: u32, chain_strand: impl Fn(&str, &[(u64, u64)]) -> Option<char>) -> Vec<Skeleton>`.

- [ ] **Step 1: Create the module with its header and the failing tests**

```rust
//! Opt-in refinement of `copy_assign --gtf`'s de novo isoform set (`--gtf-refine`).
//!
//! **STATUS:** OPT-IN
//!
//! Reached only from `src/bin/copy_assign.rs`'s `if args.gtf` block, behind `--gtf-refine` (default empty =
//! byte-identical). Four annotation-free components derived from the human chr20 error-pattern dissection;
//! rules and thresholds are frozen in `docs/superpowers/specs/2026-09-16-gtf-refine-and-dedup-fix-design.md`.
//! Not validated until the pre-registered held-out chr17 run.

use std::collections::{HashMap, HashSet};

use crate::genome::GenomeIndex;
use crate::vg_family::denovo_assemble::{PrimaryRead, Skeleton};
use crate::vg_family::family_detect::DenovoTranscript;

#[cfg(test)]
mod tests {
    use super::*;

    fn read(s: u64, e: u64, introns: Vec<(u64, u64)>) -> PrimaryRead {
        PrimaryRead { chrom: "c1".into(), ref_start: s, ref_end: e, introns, reverse: false }
    }

    fn chain(s: u64, e: u64, introns: Vec<(u64, u64)>, n: u32) -> Skeleton {
        Skeleton {
            chrom: "c1".into(), start: s, end: e, n_reads: n, introns, tied_seeded: false,
            read_strand: Some('+'), footprint: false, read_rev: 0, read_tot: n,
        }
    }

    // full chain: exons [100,200) [300,400) [500,600) [700,800); introns (200,300) (400,500) (600,700)
    fn full() -> Skeleton {
        chain(100, 800, vec![(200, 300), (400, 500), (600, 700)], 1)
    }

    #[test]
    fn plus_strand_suffix_fragment_rescues_a_one_read_chain() {
        let reads = vec![read(100, 800, full().introns.clone()), read(520, 800, vec![(600, 700)])];
        let out = fragment_supported_spliced(&[full()], &reads, 2, |_, _| Some('+'));
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].n_reads, 2, "1 exact + 1 fragment");
        assert_eq!((out[0].start, out[0].end), (100, 800), "boundaries from exact reads only");
    }

    #[test]
    fn plus_strand_prefix_is_not_three_prime_anchored() {
        let reads = vec![read(100, 800, full().introns.clone()), read(100, 350, vec![(200, 300)])];
        assert!(fragment_supported_spliced(&[full()], &reads, 2, |_, _| Some('+')).is_empty());
    }

    #[test]
    fn minus_strand_requires_the_prefix() {
        let reads = vec![read(100, 800, full().introns.clone()), read(100, 350, vec![(200, 300)])];
        assert_eq!(fragment_supported_spliced(&[full()], &reads, 2, |_, _| Some('-')).len(), 1);
    }

    #[test]
    fn a_read_extending_into_the_previous_intron_is_not_compatible() {
        // the read's chain (600,700) is full's suffix, but it starts at 380 < 500 = the acceptor of full's
        // preceding intron (400,500): it reaches into that intron, so it is not intron-compatible
        let reads = vec![read(100, 800, full().introns.clone()), read(380, 800, vec![(600, 700)])];
        assert!(fragment_supported_spliced(&[full()], &reads, 2, |_, _| Some('+')).is_empty());
    }

    #[test]
    fn an_ambiguous_fragment_abstains() {
        let other = chain(550, 900, vec![(600, 700), (750, 850)], 1); // (600,700) is NOT its last intron
        let alt = chain(450, 800, vec![(500, 550), (600, 700)], 1); // (600,700) IS its last intron
        let reads = vec![
            read(100, 800, full().introns.clone()),
            read(450, 800, alt.introns.clone()),
            read(560, 800, vec![(600, 700)]), // suffix of BOTH full and alt -> abstain
        ];
        let out = fragment_supported_spliced(&[full(), alt, other], &reads, 2, |_, _| Some('+'));
        assert!(out.is_empty(), "no chain may take an ambiguous fragment");
    }

    #[test]
    fn unspliced_reads_and_unstranded_chains_give_no_support() {
        let reads = vec![read(100, 800, full().introns.clone()), read(720, 800, vec![])];
        assert!(fragment_supported_spliced(&[full()], &reads, 2, |_, _| Some('+')).is_empty());
        let reads = vec![read(100, 800, full().introns.clone()), read(520, 800, vec![(600, 700)])];
        assert_eq!(fragment_supported_spliced(&[full()], &reads, 2, |_, _| Some('+')).len(), 1, "control");
        assert!(fragment_supported_spliced(&[full()], &reads, 2, |_, _| None).is_empty());
    }

    #[test]
    fn chains_already_at_support_are_kept_and_single_exon_or_footprint_candidates_are_ignored() {
        let two = chain(100, 800, vec![(200, 300), (400, 500), (600, 700)], 2);
        let mono = chain(1000, 1200, vec![], 5);
        let mut fp = chain(2000, 2600, vec![(2100, 2200)], 5);
        fp.footprint = true;
        let out = fragment_supported_spliced(&[two, mono, fp], &[], 2, |_, _| Some('+'));
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].n_reads, 2);
    }
}
```

- [ ] **Step 2: Register the module** (`pub mod gtf_refine;` at the end of `mod.rs`'s `pub mod` list, with a short trailing comment in the style of the neighbouring lines) and add the `docs/MODULE_STATUS.md` OPT-IN row (gate: `--gtf-refine` on `copy_assign`, default empty; evidence: reached only from `copy_assign.rs`'s `if args.gtf` block — the wiring lands in Task 5, so say "wired in Task 5 of docs/superpowers/plans/2026-09-16-gtf-refine-and-dedup-fix.md").

- [ ] **Step 3: Run tests, confirm they fail**

Run: `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo test --release --lib gtf_refine > /tmp/t3_test.log 2>&1; tail -30 /tmp/t3_test.log`
Expected: compile errors (functions not defined).

- [ ] **Step 4: Implement** (above the `#[cfg(test)]` block)

```rust
/// Strand of an intron chain when EVERY junction is canonical on one consistent strand, else `None` —
/// the strict assemble gate's rule (`build_spliced_seq`, strict branch).
pub fn strict_chain_strand(genome: &GenomeIndex, chrom: &str, introns: &[(u64, u64)]) -> Option<char> {
    let mut strand: Option<char> = None;
    for &(d, a) in introns {
        let js = if genome.is_canonical_junction(chrom, d, a, '+') {
            '+'
        } else if genome.is_canonical_junction(chrom, d, a, '-') {
            '-'
        } else {
            return None;
        };
        if strand.is_some_and(|s| s != js) {
            return None;
        }
        strand = Some(js);
    }
    strand
}

/// True iff `read` is a 3'-anchored, intron-compatible fragment of the chain `full` whose strand is
/// `full_strand`: a strictly shorter contiguous sub-chain, not reaching into `full`'s neighbouring introns,
/// ending at `full`'s 3' end ('+': last intron; '-': first intron).
fn is_three_prime_fragment(read: &PrimaryRead, full: &[(u64, u64)], full_strand: Option<char>) -> bool {
    let r = &read.introns;
    let m = r.len();
    let Some(st) = full_strand else { return false };
    if m == 0 || full.len() <= m {
        return false;
    }
    for i in 0..=(full.len() - m) {
        if full[i..i + m] != r[..] {
            continue;
        }
        if i > 0 && read.ref_start < full[i - 1].1 {
            continue;
        }
        if i + m < full.len() && read.ref_end > full[i + m].0 {
            continue;
        }
        if (st == '+' && i + m != full.len()) || (st == '-' && i != 0) {
            continue;
        }
        return true;
    }
    false
}

/// `fragsupport` (spec Part 2). `all_chains` must be `pass1_skeletons(reads, 1)`: every spliced chain with at
/// least one exact read, carrying pass-1's own boundary and strand fields. Returns the spliced, non-footprint
/// chains whose `exact + fragments >= min_support`, with `n_reads` set to that total. A read counts as a
/// fragment only when it is a fragment of EXACTLY ONE candidate (assign-or-abstain, no 1/k splitting);
/// candidates are looked up through the read's first intron. Unspliced reads never count.
pub fn fragment_supported_spliced(
    all_chains: &[Skeleton],
    reads: &[PrimaryRead],
    min_support: u32,
    chain_strand: impl Fn(&str, &[(u64, u64)]) -> Option<char>,
) -> Vec<Skeleton> {
    let cands: Vec<&Skeleton> = all_chains.iter().filter(|s| !s.introns.is_empty() && !s.footprint).collect();
    let mut by_intron: HashMap<(&str, (u64, u64)), Vec<usize>> = HashMap::new();
    for (ci, c) in cands.iter().enumerate() {
        for &it in &c.introns {
            by_intron.entry((c.chrom.as_str(), it)).or_default().push(ci);
        }
    }
    let strands: Vec<Option<char>> = cands.iter().map(|c| chain_strand(&c.chrom, &c.introns)).collect();
    let mut frag = vec![0u32; cands.len()];
    for r in reads {
        let Some(&first) = r.introns.first() else { continue };
        let Some(ids) = by_intron.get(&(r.chrom.as_str(), first)) else { continue };
        let hits: Vec<usize> =
            ids.iter().copied().filter(|&ci| is_three_prime_fragment(r, &cands[ci].introns, strands[ci])).collect();
        if hits.len() == 1 {
            frag[hits[0]] += 1;
        }
    }
    cands
        .iter()
        .enumerate()
        .filter_map(|(ci, c)| {
            let total = c.n_reads + frag[ci];
            (total >= min_support).then(|| Skeleton { n_reads: total, ..(*c).clone() })
        })
        .collect()
}
```
Note: `HashSet` and `DenovoTranscript` are imported for Task 4; if the compiler warns they are unused in this task, add `#[allow(unused_imports)]` on those two `use` lines only for this commit and remove it in Task 4.

- [ ] **Step 5: Run tests, confirm pass; run the module-status tests**

```bash
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo test --release --lib gtf_refine > /tmp/t3_test.log 2>&1; tail -15 /tmp/t3_test.log
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo test --release --lib module_status > /tmp/t3_ms.log 2>&1; tail -15 /tmp/t3_ms.log
```
Expected: 7/7 `gtf_refine` tests pass. `module_status` shows ONLY the 2 pre-existing `shared_definition.rs` failures — no mention of `gtf_refine.rs`.

- [ ] **Step 6: Commit**

```bash
git add src/rustle/vg_family/gtf_refine.rs src/rustle/vg_family/mod.rs docs/MODULE_STATUS.md
git commit -m "feat: gtf_refine module with 3'-anchored assign-or-abstain fragment support"
```

---

### Task 4: `gtf_refine.rs` — `subset` and `mono` post-filters

**Files:**
- Modify: `src/rustle/vg_family/gtf_refine.rs`

**Interfaces:**
- Consumes: `DenovoTranscript` (fields used: `chrom, start, end, strand: char, introns`; implements `Default`), `PrimaryRead`.
- Produces: `pub fn subset_removals(models: &[DenovoTranscript]) -> HashSet<usize>`; `pub fn mono_in_own_exon_removals(models: &[DenovoTranscript], removed: &HashSet<usize>) -> HashSet<usize>`; `pub fn mono_spliced_dominated_removals(models: &[DenovoTranscript], removed: &HashSet<usize>, reads: &[PrimaryRead]) -> HashSet<usize>`; `pub fn apply_post_filters(models: Vec<DenovoTranscript>, reads: &[PrimaryRead], subset: bool, mono: bool) -> Vec<DenovoTranscript>`.

- [ ] **Step 1: Write the failing tests** (append inside `mod tests`)

```rust
    fn tx(s: u64, e: u64, strand: char, introns: Vec<(u64, u64)>) -> DenovoTranscript {
        DenovoTranscript { chrom: "c1".into(), start: s, end: e, strand, introns, ..Default::default() }
    }

    // container S: exons [100,200) [300,400) [500,600) [700,800)
    fn container() -> DenovoTranscript {
        tx(100, 800, '+', vec![(200, 300), (400, 500), (600, 700)])
    }

    #[test]
    fn subset_removes_a_contained_sub_chain_and_keeps_the_container() {
        let t = tx(320, 800, '+', vec![(400, 500), (600, 700)]);
        assert_eq!(subset_removals(&[container(), t]), HashSet::from([1]));
    }

    #[test]
    fn subset_allows_5bp_overhang_only_into_a_container_intron() {
        // p = 1 (not the container's first exon): 5 bp into intron (200,300) is allowed, 6 is not
        assert_eq!(subset_removals(&[container(), tx(295, 800, '+', vec![(400, 500), (600, 700)])]), HashSet::from([1]));
        assert!(subset_removals(&[container(), tx(294, 800, '+', vec![(400, 500), (600, 700)])]).is_empty());
        // p = 0 (the container's first exon): any overhang past the container start is not allowed
        assert_eq!(subset_removals(&[container(), tx(100, 400, '+', vec![(200, 300)])]), HashSet::from([1]), "control");
        assert!(subset_removals(&[container(), tx(99, 400, '+', vec![(200, 300)])]).is_empty());
    }

    #[test]
    fn subset_never_crosses_strand_and_ignores_single_exon() {
        assert!(subset_removals(&[container(), tx(320, 800, '-', vec![(400, 500), (600, 700)])]).is_empty());
        assert!(subset_removals(&[container(), tx(320, 380, '+', vec![])]).is_empty());
    }

    #[test]
    fn mono_inside_own_spliced_exon_is_strand_aware() {
        let models = vec![container(), tx(510, 590, '+', vec![]), tx(510, 590, '-', vec![]), tx(150, 350, '+', vec![])];
        assert_eq!(mono_in_own_exon_removals(&models, &HashSet::new()), HashSet::from([1]));
        // a removed container no longer shelters anything
        assert!(mono_in_own_exon_removals(&models, &HashSet::from([0])).is_empty());
    }

    #[test]
    fn mono_spliced_dominance_counts_intronic_and_same_strand_exonic_reads() {
        let m = tx(420, 480, '+', vec![]); // len 60, inside intron (400,500)
        let spliced = PrimaryRead { chrom: "c1".into(), ref_start: 100, ref_end: 800, introns: container().introns.clone(), reverse: false };
        let unspliced = PrimaryRead { chrom: "c1".into(), ref_start: 410, ref_end: 490, introns: vec![], reverse: false };
        // 1 spliced read spans it with an intron vs 1 unspliced read -> SI (1) >= U (1) -> removed
        assert_eq!(mono_spliced_dominated_removals(&[m.clone()], &HashSet::new(), &[spliced.clone(), unspliced.clone()]), HashSet::from([0]));
        // 2 unspliced reads outvote it -> kept
        assert!(mono_spliced_dominated_removals(&[m.clone()], &HashSet::new(), &[spliced.clone(), unspliced.clone(), unspliced.clone()]).is_empty());
        // exonic evidence only counts on the model's strand
        let m2 = tx(520, 580, '+', vec![]); // inside exon [500,600)
        let rev = PrimaryRead { reverse: true, ..spliced.clone() };
        let u2 = PrimaryRead { chrom: "c1".into(), ref_start: 510, ref_end: 590, introns: vec![], reverse: false };
        assert!(mono_spliced_dominated_removals(&[m2.clone()], &HashSet::new(), &[rev, u2.clone()]).is_empty());
        assert_eq!(mono_spliced_dominated_removals(&[m2], &HashSet::new(), &[spliced, u2]), HashSet::from([0]));
    }

    #[test]
    fn apply_post_filters_runs_subset_then_mono_and_keeps_order() {
        let models = vec![container(), tx(320, 800, '+', vec![(400, 500), (600, 700)]), tx(510, 590, '+', vec![]), tx(5000, 5100, '+', vec![])];
        let u = PrimaryRead { chrom: "c1".into(), ref_start: 5000, ref_end: 5100, introns: vec![], reverse: false };
        let kept = apply_post_filters(models, &[u.clone(), u], true, true);
        let spans: Vec<(u64, u64)> = kept.iter().map(|t| (t.start, t.end)).collect();
        assert_eq!(spans, vec![(100, 800), (5000, 5100)]);
    }
```

- [ ] **Step 2: Run tests, confirm they fail**

Run: `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo test --release --lib gtf_refine > /tmp/t4_test.log 2>&1; tail -30 /tmp/t4_test.log`
Expected: compile errors for the four new functions.

- [ ] **Step 3: Implement** (remove any `#[allow(unused_imports)]` added in Task 3)

```rust
fn exons_of(start: u64, end: u64, introns: &[(u64, u64)]) -> Vec<(u64, u64)> {
    let mut ex = Vec::with_capacity(introns.len() + 1);
    let mut prev = start;
    for &(d, a) in introns {
        ex.push((prev, d));
        prev = a;
    }
    ex.push((prev, end));
    ex
}

/// `subset` (spec Part 2): indices of multi-exon models whose intron chain is a contiguous proper sub-chain of
/// another model on the same chrom and strand, with both terminal exons inside the container's matching exons
/// (up to 5 bp of overhang allowed only into a container INTRON, never past the container's own ends).
/// Computed against the whole input at once; no support condition. (The proper-sub-chain requirement equals
/// the chr20 simulation whenever chains are unique, which the gate guarantees.)
pub fn subset_removals(models: &[DenovoTranscript]) -> HashSet<usize> {
    let mut idx: HashMap<(&str, char, (u64, u64)), Vec<(usize, usize)>> = HashMap::new();
    for (mi, m) in models.iter().enumerate() {
        for (p, &it) in m.introns.iter().enumerate() {
            idx.entry((m.chrom.as_str(), m.strand, it)).or_default().push((mi, p));
        }
    }
    let mut rm = HashSet::new();
    for (ti, t) in models.iter().enumerate() {
        let k = t.introns.len();
        if k == 0 {
            continue;
        }
        let Some(hits) = idx.get(&(t.chrom.as_str(), t.strand, t.introns[0])) else { continue };
        for &(si, p) in hits {
            let s = &models[si];
            if si == ti || s.introns.len() <= k || p + k > s.introns.len() || s.introns[p..p + k] != t.introns[..] {
                continue;
            }
            let sx = exons_of(s.start, s.end, &s.introns);
            let lo = sx[p].0 as i64 - t.start as i64;
            let ro = t.end as i64 - sx[p + k].1 as i64;
            let left_ok = lo <= 0 || (p != 0 && lo <= 5);
            let right_ok = ro <= 0 || (p + k != s.introns.len() && ro <= 5);
            if left_ok && right_ok {
                rm.insert(ti);
                break;
            }
        }
    }
    rm
}

/// `mono` rule 1 (spec Part 2): single-exon models lying fully inside an exon of a not-removed multi-exon
/// model on the same chrom and strand.
pub fn mono_in_own_exon_removals(models: &[DenovoTranscript], removed: &HashSet<usize>) -> HashSet<usize> {
    let mut exons_by: HashMap<(&str, char), Vec<(u64, u64)>> = HashMap::new();
    for (i, m) in models.iter().enumerate() {
        if removed.contains(&i) || m.introns.is_empty() {
            continue;
        }
        exons_by.entry((m.chrom.as_str(), m.strand)).or_default().extend(exons_of(m.start, m.end, &m.introns));
    }
    let mut rm = HashSet::new();
    for (i, m) in models.iter().enumerate() {
        if removed.contains(&i) || !m.introns.is_empty() {
            continue;
        }
        if let Some(ex) = exons_by.get(&(m.chrom.as_str(), m.strand)) {
            if ex.iter().any(|&(s, e)| s <= m.start && e >= m.end) {
                rm.insert(i);
            }
        }
    }
    rm
}

/// `mono` rule 2 (spec Part 2): single-exon model M (not already removed) is removed when spliced reads
/// dominate its locus: `SI >= U || SE >= U`, where over the reads overlapping M, `U` = unspliced reads,
/// `SI` = spliced reads with an intron covering >= 50% of M, `SE` = spliced reads on M's strand (FLAG 0x10)
/// with an aligned exon block covering >= 50% of M.
pub fn mono_spliced_dominated_removals(
    models: &[DenovoTranscript],
    removed: &HashSet<usize>,
    reads: &[PrimaryRead],
) -> HashSet<usize> {
    let mut rm = HashSet::new();
    for (i, m) in models.iter().enumerate() {
        if removed.contains(&i) || !m.introns.is_empty() {
            continue;
        }
        let len = (m.end - m.start) as i64;
        let ov = |a: u64, b: u64| b.min(m.end) as i64 - a.max(m.start) as i64;
        let (mut u, mut si, mut se) = (0usize, 0usize, 0usize);
        for r in reads {
            if r.chrom != m.chrom || r.ref_start >= m.end || r.ref_end <= m.start {
                continue;
            }
            if r.introns.is_empty() {
                u += 1;
                continue;
            }
            let intronic = r.introns.iter().map(|&(d, a)| ov(d, a)).max().unwrap_or(i64::MIN);
            if 2 * intronic >= len {
                si += 1;
            }
            let read_strand = if r.reverse { '-' } else { '+' };
            if read_strand == m.strand {
                let exonic = exons_of(r.ref_start, r.ref_end, &r.introns).iter().map(|&(a, b)| ov(a, b)).max().unwrap_or(i64::MIN);
                if 2 * exonic >= len {
                    se += 1;
                }
            }
        }
        if si >= u || se >= u {
            rm.insert(i);
        }
    }
    rm
}

/// Apply `subset`, then `mono` (rule 1 then rule 2) when enabled; returns the kept models in input order.
pub fn apply_post_filters(
    models: Vec<DenovoTranscript>,
    reads: &[PrimaryRead],
    subset: bool,
    mono: bool,
) -> Vec<DenovoTranscript> {
    let mut removed: HashSet<usize> = if subset { subset_removals(&models) } else { HashSet::new() };
    if mono {
        let r1 = mono_in_own_exon_removals(&models, &removed);
        removed.extend(r1);
        let r2 = mono_spliced_dominated_removals(&models, &removed, reads);
        removed.extend(r2);
    }
    models.into_iter().enumerate().filter(|(i, _)| !removed.contains(i)).map(|(_, m)| m).collect()
}
```

- [ ] **Step 4: Run tests, confirm pass**

Run: `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo test --release --lib gtf_refine > /tmp/t4_test.log 2>&1; tail -20 /tmp/t4_test.log`
Expected: 13/13 pass, no warnings from `gtf_refine.rs`.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/gtf_refine.rs
git commit -m "feat: gtf_refine subset and single-exon post-filters"
```

---

### Task 5: `--gtf-refine` CLI flag and wiring into the `--gtf` block

**Files:**
- Modify: `src/bin/copy_assign.rs` (Args struct near `gtf: bool` ~line 200; the flag validation area near the `--copies-fa is only meaningful with --families` check ~line 1482; the `if args.gtf { ... }` block ~lines 2628-2660; the test module)
- Modify: `tests/copy_assign_families.rs`

**Interfaces:**
- Consumes: Task 3/4 functions; existing `assemble_gate_with(skeletons, genome, p, use_read_strand: bool, strand_margin: f64)` from `denovo_assemble`.
- Produces: `--gtf-refine` (comma list: `strand`, `subset`, `mono`, `fragsupport`, `all`).

- [ ] **Step 1: Write the failing tests**

In `src/bin/copy_assign.rs`'s test module:
```rust
    #[test]
    fn gtf_refine_parses_components_all_and_rejects_unknown() {
        let s = |v: &[&str]| v.iter().map(|x| x.to_string()).collect::<Vec<_>>();
        assert_eq!(parse_gtf_refine(&s(&[])).unwrap(), GtfRefine::default());
        let r = parse_gtf_refine(&s(&["strand", "mono"])).unwrap();
        assert!(r.strand && r.mono && !r.subset && !r.fragsupport);
        assert_eq!(parse_gtf_refine(&s(&["all"])).unwrap(), GtfRefine { strand: true, subset: true, mono: true, fragsupport: true });
        assert!(parse_gtf_refine(&s(&["bogus"])).is_err());
    }
```
In `tests/copy_assign_families.rs` (use the file's existing `scratch`, `run`, `stderr` helpers):
```rust
/// `--gtf-refine` is validated up front: it needs `--gtf`, and an unknown component is an error. Its ON-state
/// behaviour is not exercised here -- this fixture's `--gtf` emits 0 rows -- it is proven by the gtf_refine
/// unit tests and the chr20 fidelity anchors (docs/superpowers/plans/2026-09-16-gtf-refine-and-dedup-fix.md, Task 6).
#[test]
fn gtf_refine_requires_gtf_and_known_components() {
    let d = scratch("gtf_refine_cli");
    let (o, _) = run(&d, &["--no-refine", "--gtf-refine", "all"]);
    assert!(!o.status.success(), "--gtf-refine without --gtf must fail");
    assert!(stderr(&o).contains("--gtf-refine is only meaningful with --gtf"), "{}", stderr(&o));
    let (o, _) = run(&d, &["--no-refine", "--gtf", "--gtf-refine", "bogus"]);
    assert!(!o.status.success(), "an unknown component must fail");
    let (o, _) = run(&d, &["--no-refine", "--gtf", "--gtf-refine", "all"]);
    assert!(o.status.success(), "{}", stderr(&o));
}
```

- [ ] **Step 2: Run tests, confirm they fail**

```bash
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo test --release --bin copy_assign gtf_refine > /tmp/t5_unit.log 2>&1; tail -20 /tmp/t5_unit.log
```
Expected: compile error (`parse_gtf_refine` / `GtfRefine` undefined).

- [ ] **Step 3: Implement**

Args (directly after the `gtf: bool` field):
```rust
    /// Opt-in refinement of `--gtf`'s de novo isoforms (docs/superpowers/specs/2026-09-16-gtf-refine-and-dedup-fix-design.md).
    /// Comma-separated components: `strand` (single-exon strand from read orientation, margin 0.90),
    /// `subset` (drop own truncated sub-chain models), `mono` (drop single-exon models inside an own spliced
    /// exon or at spliced-read-dominated loci), `fragsupport` (count 3'-anchored truncated reads toward a
    /// chain's support, assign-or-abstain), or `all`. Annotation-free. Requires `--gtf`. Default: none
    /// (byte-identical). Thresholds are frozen by the spec; not validated until the held-out chr17 run.
    #[arg(long, value_delimiter = ',')]
    gtf_refine: Vec<String>,
```

Parser (near other bin-private helpers):
```rust
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct GtfRefine {
    strand: bool,
    subset: bool,
    mono: bool,
    fragsupport: bool,
}

fn parse_gtf_refine(items: &[String]) -> Result<GtfRefine> {
    let mut r = GtfRefine::default();
    for it in items.iter().map(|s| s.trim()).filter(|s| !s.is_empty()) {
        match it {
            "strand" => r.strand = true,
            "subset" => r.subset = true,
            "mono" => r.mono = true,
            "fragsupport" => r.fragsupport = true,
            "all" => r = GtfRefine { strand: true, subset: true, mono: true, fragsupport: true },
            other => anyhow::bail!(
                "--gtf-refine: unknown component '{other}' (expected strand, subset, mono, fragsupport, all)"
            ),
        }
    }
    Ok(r)
}
```

Validation, next to the existing `--copies-fa is only meaningful with --families` check in `main`:
```rust
    let gtf_refine = parse_gtf_refine(&args.gtf_refine)?;
    if gtf_refine != GtfRefine::default() && !args.gtf {
        anyhow::bail!("--gtf-refine is only meaningful with --gtf");
    }
```
(If that check lives before `args` is fully available or inside a different scope, place these lines at the first point in `main` where `args` is parsed and before regions are processed; `gtf_refine` is `Copy`, so the per-region `compute` closure can capture it.)

Import (extend the existing `denovo_assemble` use list with `assemble_gate_with`, and add):
```rust
use rustle::vg_family::gtf_refine::{apply_post_filters, fragment_supported_spliced, strict_chain_strand};
```

In the `if args.gtf { ... }` block, replace:
```rust
            let skeletons = pass1_skeletons(&primary, cfg.pass1_min_reads);
```
with:
```rust
            // --gtf-refine fragsupport: unspliced (and footprint) nodes exactly as before; spliced chains from
            // every >=1-read chain plus assign-or-abstain 3'-anchored fragments (gtf_refine.rs).
            let skeletons = if gtf_refine.fragsupport {
                let mut sk: Vec<_> = pass1_skeletons(&primary, cfg.pass1_min_reads)
                    .into_iter()
                    .filter(|s| s.introns.is_empty() || s.footprint)
                    .collect();
                let all_chains = pass1_skeletons(&primary, 1);
                sk.extend(fragment_supported_spliced(&all_chains, &primary, cfg.pass1_min_reads, |c, ch| {
                    strict_chain_strand(&genome, c, ch)
                }));
                sk
            } else {
                pass1_skeletons(&primary, cfg.pass1_min_reads)
            };
```
Keep the existing `RUSTLE_JUNCTION_FUZZ_BP` lines unchanged. Then replace:
```rust
            let iso = assemble_gate(&skeletons, &genome, &cfg.gate);
```
with:
```rust
            let iso = if gtf_refine.strand {
                assemble_gate_with(&skeletons, &genome, &cfg.gate, true, 0.90)
            } else {
                assemble_gate(&skeletons, &genome, &cfg.gate)
            };
            let iso = if gtf_refine.subset || gtf_refine.mono {
                apply_post_filters(iso, &primary, gtf_refine.subset, gtf_refine.mono)
            } else {
                iso
            };
```
(`let groups = collapse_loci_groups(&iso);` and the rest of the block stay unchanged.)

- [ ] **Step 4: Build and run the tests**

```bash
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo build --release --bin copy_assign > /tmp/t5_build.log 2>&1; echo "EXIT=$?"; tail -5 /tmp/t5_build.log
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo test --release --bin copy_assign gtf_refine > /tmp/t5_unit.log 2>&1; tail -10 /tmp/t5_unit.log
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo test --release --test copy_assign_families > /tmp/t5_int.log 2>&1; tail -25 /tmp/t5_int.log
```
Expected: clean build; unit test passes; `copy_assign_families` all pass (previous count + 1).

- [ ] **Step 5: Commit**

```bash
git add src/bin/copy_assign.rs tests/copy_assign_families.rs
git commit -m "feat: wire --gtf-refine into copy_assign's --gtf block"
```

---

### Task 6: chr20 fidelity anchors

**Files:**
- Create: `bench/gtf_refine_chr20_fidelity.sh`
- Modify: `bench/CHR20_ASSEMBLER_COMPARISON.md` (append a dated section)

- [ ] **Step 1: Write `bench/gtf_refine_chr20_fidelity.sh`**

```bash
#!/bin/bash
# chr20 fidelity anchors for --gtf-refine and the dedup fix
# (docs/superpowers/specs/2026-09-16-gtf-refine-and-dedup-fix-design.md, "Fidelity anchors"). Serial, foreground.
set -euo pipefail
W=/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20
F=$W/fidelity
BIN=/mnt/linuxdisk/home/juanfraitu/rustle_target/release/copy_assign
mkdir -p "$F"
arm() { # $1 label, $2 legacy(0/1), $3 refine list ('' = none)
  local label=$1 leg=$2 refine=$3
  mkdir -p "$F/$label"
  local extra=()
  [ -n "$refine" ] && extra=(--gtf-refine "$refine")
  ( cd "$F/$label"
    if [ "$leg" = 1 ]; then export RUSTLE_LEGACY_PLACEMENT_DEDUP=1; fi
    "$BIN" --gtf --bam "$W/chr20.bam" --fasta "$W/chr20.fa" --region chr20:1-66210255 "${extra[@]}" --out ours \
      > ours.stdout 2> ours.stderr )
  ( cd "$F/$label" && gffcompare -r "$W/chr20_ref.gtf" -o gffc ours.gtf > /dev/null 2>&1 )
  echo "== $label"; sed -n '/Query mRNAs/p;/Transcript level/p;/Intron chain level/p;/Locus level/p;/Matching transcripts/p' "$F/$label/gffc.stats"
}
arm legacy_none 1 ''
arm fixed_none 0 ''
arm fixed_fragsupport 0 fragsupport
arm legacy_subset 1 subset
arm legacy_strand 1 strand
arm legacy_strand_subset_mono 1 strand,subset,mono
arm fixed_all 0 all
cmp -s "$F/legacy_none/ours.gtf" "$W/ours/ours.gtf" && echo "legacy_none GTF byte-identical to the 2026-09-15 ours.gtf" || echo "legacy_none GTF DIFFERS from the 2026-09-15 ours.gtf"
```

- [ ] **Step 2: Run it** (foreground; each arm is a whole-chr20 run)

```bash
bash bench/gtf_refine_chr20_fidelity.sh > /tmp/t6_fid.log 2>&1; echo "EXIT=$?"; cat /tmp/t6_fid.log
```

- [ ] **Step 3: Check against the spec's anchors**

EXACT (query mRNAs / matching transcripts / transcript Sn/Pr / intron-chain Sn/Pr):
- `legacy_none`: 976 / 347 / 7.6 / 35.6 / 8.0 / 44.6, and its GTF byte-identical to `human_chr20/ours/ours.gtf`
- `fixed_none`: 1022 / 350 / 7.7 / 34.2 / 8.1 / 42.9
- `fixed_fragsupport`: 1054 / 362 / 7.9 / 34.3 / 8.4 / 42.7
- `legacy_subset`: 913 / 347 / 7.6 / 38.0 / 8.0 / 48.6

APPROXIMATE (explain every deviation with the spec's named differences):
- `legacy_strand`: ≈ 976 / 349, locus Pr ≈ 50.4
- `legacy_strand_subset_mono`: ≈ 820 / 349, transcript Pr ≈ 42.6, locus Pr ≈ 52.2

For any deviation: diff the arm's transcript intron chains and spans against the simulation's GTF under `human_chr20/diagnostics/` — `fixed_none` ↔ `sensloss_sim/A1_no_placement_dedup.gtf`; `fixed_fragsupport` ↔ `sensloss_sim/B3_fragment3p_unique_support2.gtf`; `legacy_subset` ↔ `precdissect_gffc/sub_own_ovh5.gtf`; `legacy_strand` ↔ `precdissect_gffc/mono_strand_from_read_vote_ge90.gtf`; `legacy_strand_subset_mono` ↔ `precdissect_gffc/P5_P1_plus_monoSplicedReadDominated.gtf` — list every differing transcript, and determine whether the Rust code or the simulation is wrong. Report DONE_WITH_CONCERNS with that per-transcript explanation if any EXACT anchor deviates; never adjust a rule to force a match.

- [ ] **Step 4: Append the results section** to `bench/CHR20_ASSEMBLER_COMPARISON.md`: "Follow-up: Rust fidelity of the dedup fix and --gtf-refine (2026-09-16)" with the full table (expected vs observed for every arm), explanations of any deviation, and the `fixed_all` row labelled "development number, never simulated as a combination — not validation".

- [ ] **Step 5: Commit**

```bash
git add bench/gtf_refine_chr20_fidelity.sh bench/CHR20_ASSEMBLER_COMPARISON.md
git commit -m "docs: chr20 fidelity anchors for the dedup fix and --gtf-refine"
```

---

### Task 7: `tss` window selection by simulation on chr20 (ADDENDUM)

**Files:**
- Create: `bench/tss_window_sim.py`
- Modify: `bench/CHR20_ASSEMBLER_COMPARISON.md` (append a dated section)

**Interfaces:**
- Consumes: `/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20/fidelity/fixed_none/ours.gtf` (Task 6, the A1 model set); `/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20/chr20_ref.gtf`; `/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20/diagnostics/sensloss_chr20bam_primary_reads.pkl` (all chr20 primaries, `-F 2308`, no placement dedup; tuples `(name, chrom, start, end, is_reverse, ts_strand, chain)`, 0-based half-open).
- Produces: `/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20/fidelity/tss_sim/{summary.tsv, chosen_W.txt, tss_W10.gtf, tss_W25.gtf, tss_W50.gtf, tss_W100.gtf}`. `chosen_W.txt` holds one integer consumed by Task 8.

- [ ] **Step 1: Write `bench/tss_window_sim.py`**

```python
#!/usr/bin/env python3
"""Choose the `tss` densest-5'-start window W on chr20 (development substrate) and write one simulated GTF per
W for the Rust fidelity check (docs/superpowers/specs/2026-09-16-gtf-refine-and-dedup-fix-design.md, `tss`
addendum).

usage: tss_window_sim.py MODELS_GTF REF_GTF READS_PKL OUT_DIR
       tss_window_sim.py --self-test
MODELS_GTF: fixed-dedup, no-refine `copy_assign --gtf` output (A1). READS_PKL: chr20 primaries (-F 2308, no
placement dedup) as (name, chrom, start, end, is_reverse, ts_strand, chain) tuples, 0-based half-open.
"""
import bisect
import collections
import os
import pickle
import re
import statistics
import sys

GRID = (10, 25, 50, 100)
TOL = 50  # SQANTI3's reference_match tolerance


def densest_five_prime(ends, window, strand):
    """Spec rule. '+': smallest p in ends maximizing #{x: p <= x <= p+W}. '-': largest p maximizing
    #{x: p-W <= x <= p}. None for no ends."""
    if not ends:
        return None
    v = sorted(ends)
    best, best_pos = -1, None
    if strand == '+':
        for lo in v:  # ascending: strict '>' keeps the smallest (most upstream) on ties
            c = bisect.bisect_right(v, lo + window) - bisect.bisect_left(v, lo)
            if c > best:
                best, best_pos = c, lo
    else:
        for hi in reversed(v):  # descending: strict '>' keeps the largest (most upstream on '-') on ties
            c = bisect.bisect_right(v, hi) - bisect.bisect_left(v, hi - window)
            if c > best:
                best, best_pos = c, hi
    return best_pos


def parse_models(path):
    """transcript_id -> dict(chrom, strand, exons, start, end, introns) from exon rows (0-based half-open)."""
    tx = collections.OrderedDict()
    for line in open(path):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon':
            continue
        tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        d = tx.setdefault(tid, dict(chrom=f[0], strand=f[6], exons=[]))
        d['exons'].append((int(f[3]) - 1, int(f[4])))
    for d in tx.values():
        d['exons'].sort()
        d['start'], d['end'] = d['exons'][0][0], d['exons'][-1][1]
        d['introns'] = tuple((a[1], b[0]) for a, b in zip(d['exons'], d['exons'][1:]))
    return tx


def five(m):
    return m['start'] if m['strand'] == '+' else m['end']


def simulate(models, by_chain, window):
    """tid -> new 5' end for every stranded multi-exon model (unchanged when it has no exact reads)."""
    out = {}
    for tid, m in models.items():
        if not m['introns'] or m['strand'] not in ('+', '-'):
            continue
        ex = by_chain.get((m['chrom'], m['introns']), [])
        pos = [s if m['strand'] == '+' else e for s, e in ex]
        out[tid] = densest_five_prime(pos, window, m['strand']) if pos else five(m)
    return out


def metrics(models, new5, ref_by_chain, ref_by_intron):
    within = moved_in = moved_out = guard = changed = n_fsm = 0
    diffs = []
    for tid, p in new5.items():
        m = models[tid]
        cur = five(m)
        changed += p != cur
        refs = ref_by_chain.get((m['chrom'], m['strand'], m['introns']))
        if refs:
            n_fsm += 1
            d_new = min(abs(p - t) for t in refs)
            d_old = min(abs(cur - t) for t in refs)
            diffs.append(d_new)
            within += d_new <= TOL
            moved_in += d_old > TOL and d_new <= TOL
            moved_out += d_old <= TOL and d_new > TOL
        near = set()
        for it in m['introns']:
            near |= ref_by_intron.get((m['chrom'], m['strand'], it), set())
        guard += any(abs(p - t) <= TOL for t in near)
    return dict(n_multi_stranded=len(new5), n_changed=changed, n_fsm_chain=n_fsm, n_within50=within,
                n_moved_in=moved_in, n_moved_out=moved_out,
                median_abs_diff=statistics.median(diffs) if diffs else 'NA', n_guard_within50=guard)


def write_gtf(models_gtf, models, new5, path):
    with open(path, 'w') as out:
        for line in open(models_gtf):
            f = line.rstrip('\n').split('\t')
            if len(f) >= 9 and f[2] in ('transcript', 'exon'):
                tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
                if tid in new5:
                    m = models[tid]
                    if m['strand'] == '+' and int(f[3]) - 1 == m['start']:
                        f[3] = str(new5[tid] + 1)
                    elif m['strand'] == '-' and int(f[4]) == m['end']:
                        f[4] = str(new5[tid])
                    line = '\t'.join(f) + '\n'
            out.write(line)


def main(models_gtf, ref_gtf, reads_pkl, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    models = parse_models(models_gtf)
    ref_by_chain = collections.defaultdict(list)
    ref_by_intron = collections.defaultdict(set)
    for r in parse_models(ref_gtf).values():
        if r['introns']:
            ref_by_chain[(r['chrom'], r['strand'], r['introns'])].append(five(r))
            for it in r['introns']:
                ref_by_intron[(r['chrom'], r['strand'], it)].add(five(r))
    by_chain = collections.defaultdict(list)
    reads = pickle.load(open(reads_pkl, 'rb'))
    for (_name, chrom, s, e, _rev, _ts, chain) in reads:
        if chain:
            by_chain[(chrom, tuple(chain))].append((s, e))
    print(f'models={len(models)} primaries={len(reads)}', file=sys.stderr)
    baseline = {tid: five(m) for tid, m in models.items() if m['introns'] and m['strand'] in ('+', '-')}
    rows = [dict(W='none', **metrics(models, baseline, ref_by_chain, ref_by_intron))]
    for w in GRID:
        new5 = simulate(models, by_chain, w)
        rows.append(dict(W=w, **metrics(models, new5, ref_by_chain, ref_by_intron)))
        write_gtf(models_gtf, models, new5, f'{out_dir}/tss_W{w}.gtf')
    cols = list(rows[0].keys())
    with open(f'{out_dir}/summary.tsv', 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for r in rows:
            fh.write('\t'.join(str(r[c]) for c in cols) + '\n')
    graded = [r for r in rows if r['W'] != 'none']
    chosen = min(graded, key=lambda r: (-r['n_within50'], r['W']))['W']
    open(f'{out_dir}/chosen_W.txt', 'w').write(f'{chosen}\n')
    print(open(f'{out_dir}/summary.tsv').read() + f'chosen_W={chosen}')


def self_test():
    assert densest_five_prime([100, 500, 505, 510], 25, '+') == 500
    assert densest_five_prime([2000, 900, 905, 910], 25, '-') == 910
    assert densest_five_prime([100, 110, 500, 510], 25, '+') == 100
    assert densest_five_prime([100, 110, 500, 510], 25, '-') == 510
    assert densest_five_prime([100], 25, '+') == 100
    assert densest_five_prime([100, 900], 25, '+') == 100
    assert densest_five_prime([100, 900], 25, '-') == 900
    assert densest_five_prime([], 25, '+') is None
    assert densest_five_prime([100, 400, 425], 25, '+') == 400
    assert densest_five_prime([100, 400, 425], 24, '+') == 100
    print('self-test OK')


if __name__ == '__main__':
    if sys.argv[1:] == ['--self-test']:
        self_test()
    else:
        main(*sys.argv[1:5])
```

- [ ] **Step 2: Self-test, then run on chr20**

```bash
python3 bench/tss_window_sim.py --self-test
W=/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20
python3 bench/tss_window_sim.py $W/fidelity/fixed_none/ours.gtf $W/chr20_ref.gtf \
  $W/diagnostics/sensloss_chr20bam_primary_reads.pkl $W/fidelity/tss_sim > /tmp/t7_tss.log 2>&1; echo EXIT=$?; cat /tmp/t7_tss.log
```
Expected: `self-test OK`; `primaries=25341` (if not, report the actual count and stop — the pickle is not the no-dedup chr20 primary set); a summary with rows `none, 10, 25, 50, 100` and `chosen_W=<n>`. Sanity: for W=`none`, `n_changed` = 0; `n_moved_in`/`n_moved_out` = 0.

- [ ] **Step 3: Append the selection section** to `bench/CHR20_ASSEMBLER_COMPARISON.md`: "Follow-up: `tss` window selection (development, chr20, 2026-09-16)" — the summary table verbatim, the chosen W, the selection rule in one line (max `n_within50`, ties → smaller W), and one line stating W is now frozen for chr17.

- [ ] **Step 4: Commit**

```bash
git add bench/tss_window_sim.py bench/CHR20_ASSEMBLER_COMPARISON.md
git commit -m "bench: choose the tss densest-start window on chr20 by simulation"
```

---

### Task 8: `tss` component in Rust (ADDENDUM)

**Files:**
- Modify: `src/rustle/vg_family/gtf_refine.rs` (new constant + two functions + tests)
- Modify: `src/bin/copy_assign.rs` (`GtfRefine.tss`, parser, CLI doc, wiring, the existing parse unit test)
- Modify: `docs/MODULE_STATUS.md` (only if the `gtf_refine.rs` row enumerates components — add `tss`)

**Interfaces:**
- Consumes: Task 7's `chosen_W.txt` value (the controller passes it in the dispatch); Task 3/4 test helpers `read(s, e, introns) -> PrimaryRead` (chrom `c1`, forward) and `tx(s, e, strand, introns) -> DenovoTranscript`; `apply_post_filters` wiring from Task 5.
- Produces: `pub const TSS_WINDOW_BP: u64`; `pub fn densest_five_prime(ends: &[u64], window: u64, strand: char) -> Option<u64>`; `pub fn refine_tss(models: Vec<DenovoTranscript>, reads: &[PrimaryRead], window: u64) -> Vec<DenovoTranscript>`; `--gtf-refine tss` (and `all` = five components).

- [ ] **Step 1: Write the failing tests** (append inside `mod tests` in `gtf_refine.rs`)

```rust
    #[test]
    fn densest_five_prime_skips_a_lone_upstream_outlier_on_both_strands() {
        assert_eq!(densest_five_prime(&[100, 500, 505, 510], 25, '+'), Some(500));
        assert_eq!(densest_five_prime(&[2000, 900, 905, 910], 25, '-'), Some(910));
    }

    #[test]
    fn densest_five_prime_ties_keep_the_most_upstream_window() {
        assert_eq!(densest_five_prime(&[100, 110, 500, 510], 25, '+'), Some(100));
        assert_eq!(densest_five_prime(&[100, 110, 500, 510], 25, '-'), Some(510));
    }

    #[test]
    fn densest_five_prime_never_moves_one_or_two_reads() {
        assert_eq!(densest_five_prime(&[100], 25, '+'), Some(100));
        assert_eq!(densest_five_prime(&[100, 900], 25, '+'), Some(100));
        assert_eq!(densest_five_prime(&[100, 900], 25, '-'), Some(900));
        assert_eq!(densest_five_prime(&[], 25, '+'), None);
    }

    #[test]
    fn densest_five_prime_window_is_inclusive() {
        assert_eq!(densest_five_prime(&[100, 400, 425], 25, '+'), Some(400));
        assert_eq!(densest_five_prime(&[100, 400, 425], 24, '+'), Some(100));
    }

    #[test]
    fn refine_tss_moves_only_stranded_multi_exon_5prime_ends_from_exact_chain_reads() {
        let chain = vec![(200, 300), (400, 500), (600, 700)];
        let minus = vec![(200, 300), (400, 500)];
        let other = vec![(200, 300), (400, 500), (600, 650)];
        let on_c2 = |s: u64| PrimaryRead { chrom: "c2".into(), ref_start: s, ref_end: 800, introns: chain.clone(), reverse: false };
        let mut reads = vec![
            read(10, 800, chain.clone()), // lone upstream outlier
            read(150, 800, chain.clone()),
            read(155, 800, chain.clone()),
            read(160, 800, chain.clone()),
            read(520, 800, vec![(600, 700)]), // fragment: its chain differs, never counts
            read(20, 800, other.clone()),     // other chain: if counted, the window at 10 would win (4 > 3)
            read(21, 800, other.clone()),
            read(22, 800, other.clone()),
        ];
        reads.extend([on_c2(190), on_c2(191), on_c2(192), on_c2(193)]); // other chromosome: would win if counted
        for e in [990, 840, 845, 850] {
            reads.push(read(100, e, minus.clone())); // '-' model: outlier 990, cluster 840..850
        }
        let models = vec![
            tx(10, 800, '+', chain.clone()),
            tx(100, 990, '-', minus.clone()),
            tx(1000, 1200, '+', vec![]),             // single-exon: unchanged
            tx(5000, 6000, '+', vec![(5100, 5200)]), // no exact reads: unchanged
            tx(10, 800, '.', chain.clone()),         // unstranded: unchanged
        ];
        let spans: Vec<(u64, u64)> = refine_tss(models, &reads, 25).iter().map(|t| (t.start, t.end)).collect();
        assert_eq!(spans, vec![(150, 800), (100, 850), (1000, 1200), (5000, 6000), (10, 800)]);
    }
```

In `src/bin/copy_assign.rs`, replace the body of `gtf_refine_parses_components_all_and_rejects_unknown` with:
```rust
        let s = |v: &[&str]| v.iter().map(|x| x.to_string()).collect::<Vec<_>>();
        assert_eq!(parse_gtf_refine(&s(&[])).unwrap(), GtfRefine::default());
        let r = parse_gtf_refine(&s(&["strand", "mono"])).unwrap();
        assert!(r.strand && r.mono && !r.subset && !r.fragsupport && !r.tss);
        assert!(parse_gtf_refine(&s(&["tss"])).unwrap().tss);
        assert_eq!(
            parse_gtf_refine(&s(&["all"])).unwrap(),
            GtfRefine { strand: true, subset: true, mono: true, fragsupport: true, tss: true }
        );
        assert!(parse_gtf_refine(&s(&["bogus"])).is_err());
```

- [ ] **Step 2: Run tests, confirm they fail**

```bash
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo test --release --lib gtf_refine > /tmp/t8_red.log 2>&1; grep -n "error\[" /tmp/t8_red.log | head
```
Expected: compile errors (`densest_five_prime`, `refine_tss` not found).

- [ ] **Step 3: Implement in `gtf_refine.rs`** (above `#[cfg(test)]`; replace `<W>` with the integer from Task 7's `chosen_W.txt`, given in your dispatch)

```rust
/// `tss` window W (bp): chosen once on chr20 by `bench/tss_window_sim.py` (max models within 50 bp of the
/// chain-matched reference TSS, ties → smaller W; grid 10/25/50/100) and frozen before the chr17
/// pre-registration. Do not tune.
pub const TSS_WINDOW_BP: u64 = <W>;

/// Spec `tss` rule: the read 5' end with the most read 5' ends inside the `window`-bp window extending
/// downstream from it (`'+'`: `[p, p + W]`; `'-'`: `[p - W, p]`), ties to the most upstream. `None` when
/// `ends` is empty. One or two ends never move off the most upstream one.
pub fn densest_five_prime(ends: &[u64], window: u64, strand: char) -> Option<u64> {
    let mut v = ends.to_vec();
    v.sort_unstable();
    let mut best: Option<(usize, u64)> = None;
    if strand == '-' {
        for &hi in v.iter().rev() {
            let c = v.partition_point(|&x| x <= hi) - v.partition_point(|&x| x < hi.saturating_sub(window));
            if best.map_or(true, |(b, _)| c > b) {
                best = Some((c, hi));
            }
        }
    } else {
        for &lo in &v {
            let c = v.partition_point(|&x| x <= lo.saturating_add(window)) - v.partition_point(|&x| x < lo);
            if best.map_or(true, |(b, _)| c > b) {
                best = Some((c, lo));
            }
        }
    }
    best.map(|(_, p)| p)
}

/// `tss` (spec addendum): move each stranded multi-exon model's 5' end to `densest_five_prime` over the 5'
/// ends of reads with exactly its intron chain on its chromosome. Single-exon models, unstranded models,
/// models with no exact reads, and every 3' end are unchanged. Only `start`/`end` change — `seq` is not
/// recomputed, because the `--gtf` block consumes coordinates only.
pub fn refine_tss(mut models: Vec<DenovoTranscript>, reads: &[PrimaryRead], window: u64) -> Vec<DenovoTranscript> {
    let mut by_chain: HashMap<(&str, &[(u64, u64)]), Vec<&PrimaryRead>> = HashMap::new();
    for r in reads.iter().filter(|r| !r.introns.is_empty()) {
        by_chain.entry((r.chrom.as_str(), r.introns.as_slice())).or_default().push(r);
    }
    for m in models.iter_mut() {
        if m.introns.is_empty() || (m.strand != '+' && m.strand != '-') {
            continue;
        }
        let Some(exact) = by_chain.get(&(m.chrom.as_str(), m.introns.as_slice())) else { continue };
        let ends: Vec<u64> = exact.iter().map(|r| if m.strand == '+' { r.ref_start } else { r.ref_end }).collect();
        if let Some(p) = densest_five_prime(&ends, window, m.strand) {
            if m.strand == '+' {
                m.start = p;
            } else {
                m.end = p;
            }
        }
    }
    models
}
```

- [ ] **Step 4: Wire it into `src/bin/copy_assign.rs`**

- Import: extend the `gtf_refine` use line to `use rustle::vg_family::gtf_refine::{apply_post_filters, fragment_supported_spliced, refine_tss, strict_chain_strand, TSS_WINDOW_BP};`.
- `GtfRefine`: add field `tss: bool`. `parse_gtf_refine`: add arm `"tss" => r.tss = true,`; make `"all"` set all five fields; the unknown-component message lists `strand, subset, mono, fragsupport, tss, all`.
- CLI doc comment on `gtf_refine`: add "`tss` (5' end = densest window of exact-read 5' ends, W = `TSS_WINDOW_BP`)" to the component list.
- In the `if args.gtf` block, directly after the `apply_post_filters` statement and before `collapse_loci_groups`:
```rust
            let iso = if gtf_refine.tss { refine_tss(iso, &primary, TSS_WINDOW_BP) } else { iso };
```
- Before editing, confirm that after this point `iso` is consumed only by `collapse_loci_groups(&iso)` and the `TranscriptRec` map (coordinates only); if anything reads `seq`, report NEEDS_CONTEXT.

- [ ] **Step 5: Run tests and build**

```bash
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo test --release --lib gtf_refine > /tmp/t8_green.log 2>&1; grep -n "test result\|warning.*gtf_refine\|gtf_refine.rs" /tmp/t8_green.log
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo test --release --bin copy_assign > /tmp/t8_bin.log 2>&1; grep -n "test result" /tmp/t8_bin.log
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo test --release --test copy_assign_families > /tmp/t8_int.log 2>&1; grep -n "test result" /tmp/t8_int.log
```
Expected: gtf_refine 20 passed (15 + 5), no warning lines into `gtf_refine.rs` or `copy_assign.rs`; copy_assign bin 28 passed; copy_assign_families 15 passed.

- [ ] **Step 6: Commit**

```bash
git add src/rustle/vg_family/gtf_refine.rs src/bin/copy_assign.rs docs/MODULE_STATUS.md
git commit -m "feat: tss densest-start 5' end component for --gtf-refine"
```

---

### Task 9: chr20 fidelity + SQANTI3 check for `tss` (ADDENDUM)

**Files:**
- Modify: `bench/gtf_refine_chr20_fidelity.sh` (add arm `fixed_tss` = fixed dedup + `--gtf-refine tss`)
- Modify: `bench/CHR20_ASSEMBLER_COMPARISON.md` (append results to the Task 7 section)

**Interfaces:**
- Consumes: Task 7's `tss_sim/tss_W<chosen>.gtf`; Task 8's binary.
- Produces: `fidelity/fixed_tss/ours.gtf` + gffcompare; re-run `fidelity/fixed_all/` (now five components); SQANTI3 on `fixed_none` and `fixed_tss` under `/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20/sqanti3_tss/<label>/`.

- [ ] **Step 1: Build, add the arm, run `fixed_tss` and `fixed_all`** (each a separate foreground call)

```bash
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo build --release --bin copy_assign > /tmp/t9_build.log 2>&1; echo EXIT=$?
bash bench/gtf_refine_chr20_fidelity.sh fixed_tss > /tmp/t9_fixed_tss.log 2>&1; cat /tmp/t9_fixed_tss.log
bash bench/gtf_refine_chr20_fidelity.sh fixed_all > /tmp/t9_fixed_all.log 2>&1; cat /tmp/t9_fixed_all.log
```
Expected for `fixed_tss`: query mRNAs 1022, matching transcripts 350, transcript 7.7/34.2, intron chain 8.1/42.9 (identical to A1 — multi-exon ends are not scored).

- [ ] **Step 2: EXACT per-transcript check vs the simulation**

```bash
W=/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20; CW=$(cat $W/fidelity/tss_sim/chosen_W.txt)
python3 - "$W/fidelity/fixed_tss/ours.gtf" "$W/fidelity/tss_sim/tss_W${CW}.gtf" <<'EOF'
import collections, re, sys
def keys(path):
    tx = collections.defaultdict(lambda: dict(ex=[]))
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon':
            continue
        t = tx[re.search(r'transcript_id "([^"]+)"', f[8]).group(1)]
        t['chrom'], t['strand'] = f[0], f[6]
        t['ex'].append((int(f[3]) - 1, int(f[4])))
    out = collections.Counter()
    for t in tx.values():
        ex = sorted(t['ex'])
        out[(t['chrom'], t['strand'], tuple((a[1], b[0]) for a, b in zip(ex, ex[1:])), ex[0][0], ex[-1][1])] += 1
    return out
a, b = keys(sys.argv[1]), keys(sys.argv[2])
print('rust', sum(a.values()), 'sim', sum(b.values()), 'only_rust', sum((a - b).values()), 'only_sim', sum((b - a).values()))
for k in list((a - b).elements())[:20]: print('only_rust', k)
for k in list((b - a).elements())[:20]: print('only_sim', k)
EOF
```
Expected: `only_rust 0 only_sim 0`. Any difference: explain every transcript (classify intended spec difference / simulation bug / Rust bug). A Rust bug → DONE_WITH_CONCERNS with evidence, no code change.

- [ ] **Step 3: SQANTI3 on `fixed_none` and `fixed_tss`** (one per call; if a run exceeds the 10-min tool cap, background it and end your turn with `WAITING PID=<pid> <label>`)

```bash
W=/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20
for L in fixed_none fixed_tss; do  # run the loop body once per label, as separate calls
  source /home/juanfra/miniforge3/etc/profile.d/conda.sh; conda activate sqanti3
  mkdir -p $W/sqanti3_tss/$L; cd /mnt/linuxdisk/home/juanfraitu/_from_wsl/tools/SQANTI3
  python sqanti3_qc.py --isoforms $W/fidelity/$L/ours.gtf --refGTF $W/chr20_ref.gtf --refFasta $W/chr20.fa \
    -o $L -d $W/sqanti3_tss/$L --report skip -t 4 > $W/sqanti3_tss/$L.qc.log 2>&1; echo "$L exit=$?"
done
```

- [ ] **Step 4: E5-style metrics on chr20 (descriptive)**

```bash
W=/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20
python3 - $W/sqanti3_tss/fixed_none/fixed_none_classification.txt $W/sqanti3_tss/fixed_tss/fixed_tss_classification.txt <<'EOF'
import csv, sys
for path in sys.argv[1:]:
    n = within = gene = 0
    for r in csv.DictReader(open(path), delimiter='\t'):
        if r['structural_category'] != 'full-splice_match' or r['subcategory'] == 'mono-exon':
            continue
        n += 1
        try: within += abs(float(r['diff_to_TSS'])) <= 50
        except ValueError: pass
        try: gene += abs(float(r['diff_to_gene_TSS'])) <= 50
        except ValueError: pass
    print(path.split('/')[-1], 'n', n, 'p', round(within / n, 4), 'within', within, 'g', gene)
EOF
```
Report both rows (no pass/fail on chr20; this is development confirmation with SQANTI3's own reference choice).

- [ ] **Step 5: Append results and commit**

Append to the Task 7 section of `bench/CHR20_ASSEMBLER_COMPARISON.md`: `fixed_tss` gffcompare line vs A1; the per-transcript check result; the SQANTI3 p/within/g table for `fixed_none` vs `fixed_tss`; the new `fixed_all` dev row (labelled development-only).

```bash
git add bench/gtf_refine_chr20_fidelity.sh bench/CHR20_ASSEMBLER_COMPARISON.md
git commit -m "bench: chr20 fidelity and SQANTI3 5'-end check for tss"
```

---

### Task 10: Pre-register the chr17 evaluation (before any chr17 data exists)

**Files:**
- Create: `docs/PREREG_gtf_refine_chr17_2026-09-16.md`
- Create: `bench/gtf_refine_verdict.py`

- [ ] **Step 1: Confirm no chr17 data exists yet**

Run: `ls /mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr17 2>&1`
Expected: `No such file or directory`. If it exists, STOP and report BLOCKED (the pre-registration would no longer be blind).

- [ ] **Step 2: Write `bench/gtf_refine_verdict.py`**

```python
#!/usr/bin/env python3
"""Pre-registered verdict for the held-out --gtf-refine evaluation (docs/PREREG_gtf_refine_chr17_2026-09-16.md).

usage: gtf_refine_verdict.py BASE_TMAP BASE_STATS BASE_JUNCTIONS BUNDLE_TMAP BUNDLE_STATS BUNDLE_JUNCTIONS \
           BASE_CLASSIFICATION ABL_TSS_CLASSIFICATION
       gtf_refine_verdict.py --self-test
"""
import csv
import json
import re
import sys


def eq_refs(tmap):
    with open(tmap) as fh:
        return {r['ref_id'] for r in csv.DictReader(fh, delimiter='\t') if r['class_code'] == '='}


def precisions(stats):
    txt = open(stats).read()
    tr = float(re.search(r'Transcript level:\s+[\d.]+\s+\|\s+([\d.]+)', txt).group(1))
    ic = float(re.search(r'Intron chain level:\s+[\d.]+\s+\|\s+([\d.]+)', txt).group(1))
    return tr, ic


def novel_canonical_junctions(junctions):
    out = set()
    with open(junctions) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r['junction_category'] == 'novel' and r['canonical'] == 'canonical':
                out.add((r['chrom'], r['strand'], r['genomic_start_coord'], r['genomic_end_coord']))
    return out


def verdict(base_eq, base_pr, base_nj, bun_eq, bun_pr, bun_nj):
    e1 = len(bun_eq) >= len(base_eq)
    e2 = bun_pr[0] > base_pr[0] and bun_pr[1] > base_pr[1]
    lost = base_eq - bun_eq
    e3 = len(lost) <= 0.01 * len(base_eq)
    e4 = len(bun_nj) >= len(base_nj)
    if e1 and e2 and e3 and e4:
        v = 'SUPPORTED'
    elif e3 and e4 and (e1 != e2):
        v = 'PARTIAL'
    else:
        v = 'REFUTED'
    return dict(E1_matches=e1, E2_precision=e2, E3_collateral=e3, E4_novel_junctions=e4, verdict=v,
                n_base_eq=len(base_eq), n_bundle_eq=len(bun_eq), n_lost=len(lost), lost_refs=sorted(lost),
                base_tx_ic_pr=base_pr, bundle_tx_ic_pr=bun_pr,
                n_base_novel_junctions=len(base_nj), n_bundle_novel_junctions=len(bun_nj))


def tss_metrics(classification):
    """E5 inputs over SQANTI3 multi-exon full-splice matches: n, p = fraction |diff_to_TSS| <= 50 (unrounded),
    g = count |diff_to_gene_TSS| <= 50. Rows with a non-numeric diff count in n but not in p/g."""
    n = within = gene = 0
    with open(classification) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r['structural_category'] != 'full-splice_match' or r['subcategory'] == 'mono-exon':
                continue
            n += 1
            try:
                within += abs(float(r['diff_to_TSS'])) <= 50
            except ValueError:
                pass
            try:
                gene += abs(float(r['diff_to_gene_TSS'])) <= 50
            except ValueError:
                pass
    return dict(n=n, within=within, p=within / n if n else 0.0, g=gene)


def tss_verdict(base, abl):
    e5 = abl['p'] > base['p'] and abl['g'] >= base['g']
    return dict(E5_tss=e5, tss_verdict='SUPPORTED' if e5 else 'REFUTED', baseline_tss=base, abl_tss=abl)


def self_test():
    tb = dict(n=100, within=40, p=0.40, g=60)
    assert tss_verdict(tb, dict(n=100, within=45, p=0.45, g=60))['tss_verdict'] == 'SUPPORTED'
    assert tss_verdict(tb, dict(n=100, within=40, p=0.40, g=70))['tss_verdict'] == 'REFUTED'  # p tie fails
    assert tss_verdict(tb, dict(n=100, within=45, p=0.45, g=59))['tss_verdict'] == 'REFUTED'  # guard fails
    b = set(range(100))
    nj = {('c', '+', '1', '2')}
    assert verdict(b, (35.0, 44.0), nj, b | {100}, (40.0, 48.0), nj)['verdict'] == 'SUPPORTED'
    assert verdict(b, (35.0, 44.0), nj, b, (35.0, 48.0), nj)['verdict'] == 'PARTIAL'  # E2 fails on tx Pr tie
    assert verdict(b, (35.0, 44.0), nj, b - {0}, (40.0, 48.0), nj)['verdict'] == 'PARTIAL'  # E1 fails, 1 lost <= 1%
    assert verdict(b, (35.0, 44.0), nj, b - {0, 1}, (40.0, 48.0), nj)['verdict'] == 'REFUTED'  # E3: 2 lost > 1
    assert verdict(b, (35.0, 44.0), nj, b, (40.0, 48.0), set())['verdict'] == 'REFUTED'  # E4 fails
    assert verdict(b, (35.0, 44.0), nj, b - {0}, (34.0, 44.0), nj)['verdict'] == 'REFUTED'  # neither E1 nor E2
    print('self-test OK')


if __name__ == '__main__':
    if sys.argv[1:] == ['--self-test']:
        self_test()
        sys.exit(0)
    bt, bs, bj, nt, ns, nj_, bc, tc = sys.argv[1:9]
    res = verdict(eq_refs(bt), precisions(bs), novel_canonical_junctions(bj),
                  eq_refs(nt), precisions(ns), novel_canonical_junctions(nj_))
    res.update(tss_verdict(tss_metrics(bc), tss_metrics(tc)))
    print(json.dumps(res, indent=1))
```

- [ ] **Step 3: Run the self-test**

Run: `python3 bench/gtf_refine_verdict.py --self-test`
Expected: `self-test OK`.

- [ ] **Step 4: Write the PREREG doc**

Follow the house style of `docs/PREREG_core_definition_2026-09-12.md` (title line stating it was written before any score is computed; `## Question`; method; decision rules). Content must include, verbatim from the spec: the substrate (chr17 of `human_testis.t2t.bam`, CHM13 genome, chr17 subset of `chm13v2.0_RefSeq_full.gff.gz` prepared exactly like `bench/prep_chr20_ref.sh`, outputs under `/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr17/`); the compared arms (BASELINE = fixed dedup, no `--gtf-refine`; BUNDLE = fixed dedup, `--gtf-refine all` = all five components; ABL_TSS = fixed dedup, `--gtf-refine tss`, scored for E5 only); the frozen component rules and thresholds, including the `tss` rule with the literal frozen value of `TSS_WINDOW_BP` (read it from `src/rustle/vg_family/gtf_refine.rs` and cite Task 7's selection section); the descriptive-only runs (legacy, the other four single-component ablations, StringTie `-L`, FLAIR unguided); E1–E4 and their verdict rule, and E5 and the separate TSS verdict, exactly as in the spec; the E5 scoring convention (SQANTI3 `_classification.txt` rows with `structural_category == full-splice_match` and `subcategory != mono-exon`; `p` unrounded; non-numeric diffs count in n only); the scoring conventions: `'='` = distinct `ref_id` with `class_code == '='` in the gffcompare `.tmap`; precision = the one-decimal value in the gffcompare `.stats` line (a tie at one decimal fails E2); novel junction = distinct `(chrom, strand, genomic_start_coord, genomic_end_coord)` with `junction_category == novel` and `canonical == canonical` in SQANTI3's `_junctions.txt`; the verdict is computed by `bench/gtf_refine_verdict.py` (this commit); and the sentence "No threshold, component, or rule may change after any chr17 number is seen."

- [ ] **Step 5: Commit** (this commit must precede any chr17 file)

```bash
git add docs/PREREG_gtf_refine_chr17_2026-09-16.md bench/gtf_refine_verdict.py
git commit -m "docs: pre-register the held-out chr17 --gtf-refine evaluation"
```

---

### Task 11: chr17 substrate and all runs

**Files:**
- Create: `bench/prep_chrom_ref.sh`, `bench/bakeoff_chrom_ours.sh`, `bench/bakeoff_chrom_stringtie.sh`, `bench/bakeoff_chrom_flair.sh`, `bench/chrom_score.sh`

**Interfaces:**
- Consumes: the committed PREREG (Task 10).
- Produces under `/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr17/`: `chr17.bam`, `chr17.fa`, `chr17_ref.gtf`; `<label>/ours.gtf` for labels `baseline`, `bundle`, `legacy`, `abl_strand`, `abl_subset`, `abl_mono`, `abl_fragsupport`, `abl_tss`; `stringtie/st.gtf`; `flair/flair.isoforms.gtf`; `gffcompare/<label>.{stats,tmap-in-query-dir}` for all 10; `sqanti3/<label>/<label>_{classification,junctions}.txt` for `baseline`, `bundle`, `abl_tss`, `stringtie`, `flair`.

- [ ] **Step 1: Confirm the PREREG commit exists**

Run: `git log --oneline -- docs/PREREG_gtf_refine_chr17_2026-09-16.md | head -1`
Expected: one commit. If empty, STOP (BLOCKED).

- [ ] **Step 2: Write the generalized scripts** (do NOT modify the committed chr20 scripts)

`bench/prep_chrom_ref.sh`: copy `bench/prep_chr20_ref.sh` verbatim, then parameterize: `CHROM=${1:?usage: prep_chrom_ref.sh CHROM}`, `W=/mnt/linuxdisk/home/juanfraitu/bakeoff/human_${CHROM}`, and replace every literal `chr20` in paths/filters with `${CHROM}` (keep the embedded Python resort block and the gffread step unchanged).

`bench/bakeoff_chrom_ours.sh`:
```bash
#!/bin/bash
# usage: bakeoff_chrom_ours.sh CHROM LABEL [extra copy_assign args...]
# Runs copy_assign --gtf over the whole chromosome into $W/LABEL/ours.gtf. Env vars (e.g.
# RUSTLE_LEGACY_PLACEMENT_DEDUP) pass through from the caller.
set -euo pipefail
CHROM=${1:?CHROM}; LABEL=${2:?LABEL}; shift 2
W=/mnt/linuxdisk/home/juanfraitu/bakeoff/human_${CHROM}
BIN=/mnt/linuxdisk/home/juanfraitu/rustle_target/release/copy_assign
LEN=$(awk -v c="$CHROM" '$1==c {print $2}' "$W/${CHROM}.fa.fai")
mkdir -p "$W/$LABEL"; cd "$W/$LABEL"
"$BIN" --gtf --bam "$W/${CHROM}.bam" --fasta "$W/${CHROM}.fa" --region "${CHROM}:1-${LEN}" "$@" --out ours \
  > ours.stdout.log 2> ours.stderr.log
echo "$LABEL exit=$? transcripts=$(awk -F'\t' '$3=="transcript"' ours.gtf | wc -l)"
```

`bench/bakeoff_chrom_stringtie.sh` and `bench/bakeoff_chrom_flair.sh`: copies of the chr20 versions with `CHROM=${1:?CHROM}` and `W=/mnt/linuxdisk/home/juanfraitu/bakeoff/human_${CHROM}`, every `chr20` path replaced by `${CHROM}`; keep the FLAIR PATH workaround and the skipped `correct` step with their comments.

`bench/chrom_score.sh`:
```bash
#!/bin/bash
# usage: chrom_score.sh CHROM LABEL GTF [--sqanti]
set -euo pipefail
CHROM=${1:?CHROM}; LABEL=${2:?LABEL}; GTF=${3:?GTF}; SQ_FLAG=${4:-}
W=/mnt/linuxdisk/home/juanfraitu/bakeoff/human_${CHROM}
mkdir -p "$W/gffcompare"; cd "$W/gffcompare"
gffcompare -r "$W/${CHROM}_ref.gtf" -o "$LABEL" "$GTF" > "$LABEL.gffcompare.log" 2>&1
sed -n '/Query mRNAs/p;/Transcript level/p;/Intron chain level/p;/Locus level/p;/Matching transcripts/p' "$LABEL.stats"
if [ "$SQ_FLAG" = "--sqanti" ]; then
  source /home/juanfra/miniforge3/etc/profile.d/conda.sh; conda activate sqanti3
  mkdir -p "$W/sqanti3/$LABEL"; cd /mnt/linuxdisk/home/juanfraitu/_from_wsl/tools/SQANTI3
  python sqanti3_qc.py --isoforms "$GTF" --refGTF "$W/${CHROM}_ref.gtf" --refFasta "$W/${CHROM}.fa" \
    -o "$LABEL" -d "$W/sqanti3/$LABEL" --report skip -t 4 > "$W/sqanti3/$LABEL.qc.log" 2>&1
  echo "SQANTI3 $LABEL exit=$?"
fi
```
(gffcompare writes the `.tmap`/`.refmap` next to the query GTF as `<dir>/<LABEL>.<gtf basename>.tmap`.)

- [ ] **Step 3: Rebuild, then run everything serially in the foreground** (redirect each to a log; expect a long total runtime — FLAIR and the five SQANTI3 runs dominate)

```bash
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo build --release --bin copy_assign > /tmp/t8_build.log 2>&1; echo "EXIT=$?"
bash bench/prep_chrom_ref.sh chr17 > /tmp/t8_prep.log 2>&1; echo "prep EXIT=$?"; tail -8 /tmp/t8_prep.log
bash bench/bakeoff_chrom_ours.sh chr17 baseline > /tmp/t8_baseline.log 2>&1; cat /tmp/t8_baseline.log
bash bench/bakeoff_chrom_ours.sh chr17 bundle --gtf-refine all > /tmp/t8_bundle.log 2>&1; cat /tmp/t8_bundle.log
RUSTLE_LEGACY_PLACEMENT_DEDUP=1 bash bench/bakeoff_chrom_ours.sh chr17 legacy > /tmp/t8_legacy.log 2>&1; cat /tmp/t8_legacy.log
for c in strand subset mono fragsupport tss; do
  bash bench/bakeoff_chrom_ours.sh chr17 abl_$c --gtf-refine $c > /tmp/t8_abl_$c.log 2>&1; cat /tmp/t8_abl_$c.log
done
bash bench/bakeoff_chrom_stringtie.sh chr17 > /tmp/t8_st.log 2>&1; echo "stringtie EXIT=$?"; tail -3 /tmp/t8_st.log
bash bench/bakeoff_chrom_flair.sh chr17 > /tmp/t8_flair.log 2>&1; echo "flair EXIT=$?"; tail -5 /tmp/t8_flair.log
W=/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr17
for l in legacy abl_strand abl_subset abl_mono abl_fragsupport; do bash bench/chrom_score.sh chr17 $l $W/$l/ours.gtf; done
bash bench/chrom_score.sh chr17 baseline $W/baseline/ours.gtf --sqanti
bash bench/chrom_score.sh chr17 bundle $W/bundle/ours.gtf --sqanti
bash bench/chrom_score.sh chr17 abl_tss $W/abl_tss/ours.gtf --sqanti
bash bench/chrom_score.sh chr17 stringtie $W/stringtie/st.gtf --sqanti
bash bench/chrom_score.sh chr17 flair $W/flair/flair.isoforms.gtf --sqanti
```
Do NOT compute or look at the E1–E5 verdicts in this task (Task 12 applies the pre-registered script). If any tool fails, report BLOCKED with the log tail; do not change any rule or threshold.

- [ ] **Step 4: Commit the scripts** (data stays out of git)

```bash
git add bench/prep_chrom_ref.sh bench/bakeoff_chrom_ours.sh bench/bakeoff_chrom_stringtie.sh bench/bakeoff_chrom_flair.sh bench/chrom_score.sh
git commit -m "bench: chromosome-parameterized assembler bakeoff scripts; chr17 held-out runs"
```

---

### Task 12: Apply the pre-registered verdicts and write them up

**Files:**
- Create: `bench/CHR17_GTF_REFINE_VALIDATION.md`
- Modify (only if the E1–E4 verdict is PARTIAL or REFUTED, or the TSS verdict is REFUTED): `docs/NEGATIVE_RESULTS_REGISTER.md` (one new row per failed verdict, in its existing format)

- [ ] **Step 1: Compute the verdicts with the committed script (no edits to it)**

```bash
W=/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr17
BT=$(ls $W/baseline/baseline.*.tmap); NT=$(ls $W/bundle/bundle.*.tmap)
python3 bench/gtf_refine_verdict.py "$BT" $W/gffcompare/baseline.stats $W/sqanti3/baseline/baseline_junctions.txt \
  "$NT" $W/gffcompare/bundle.stats $W/sqanti3/bundle/bundle_junctions.txt \
  $W/sqanti3/baseline/baseline_classification.txt $W/sqanti3/abl_tss/abl_tss_classification.txt > /tmp/t9_verdict.json
cat /tmp/t9_verdict.json
git log --oneline -1 -- bench/gtf_refine_verdict.py docs/PREREG_gtf_refine_chr17_2026-09-16.md
```
(Confirm with `git log` that the script and PREREG are unchanged since Task 10's commit.)

- [ ] **Step 2: Gather descriptive tables**

- gffcompare summary (query mRNAs, matching transcripts, transcript / intron-chain / locus Sn and Pr) for all 10 labels, from `$W/gffcompare/<label>.stats`.
- SQANTI3 structural category counts for `baseline`, `bundle`, `abl_tss`, `stringtie`, `flair` (column 6 of each `_classification.txt`, as in `bench/chr20_score.sh`), plus the E5 metrics (n, within, p, g) for the same five.
- The `lost_refs` list from the verdict JSON, each with what the bundle emitted at that reference (from the bundle `.tmap`/`.refmap`).
- Distinct canonical novel junctions: baseline, bundle, and how many of the bundle's are absent from the baseline.

- [ ] **Step 3: Write `bench/CHR17_GTF_REFINE_VALIDATION.md`**

Sections: purpose and link to the PREREG commit; substrate; exact commands; the verdict JSON verbatim with E1–E4 spelled out and the separate E5/TSS verdict spelled out; the descriptive tables; the ablation table (each component alone vs baseline); StringTie/FLAIR context; an honest one-paragraph reading. State the verdict plainly whatever it is; do not reinterpret a failed endpoint.

- [ ] **Step 4: Register a negative result if applicable**

If the E1–E4 verdict is PARTIAL or REFUTED, or the TSS verdict is REFUTED, add one row per failed verdict to `docs/NEGATIVE_RESULTS_REGISTER.md` in its existing column format (check the last row number), citing the PREREG doc and the failed endpoint(s) with the real numbers.

- [ ] **Step 5: Commit**

```bash
git add bench/CHR17_GTF_REFINE_VALIDATION.md docs/NEGATIVE_RESULTS_REGISTER.md
git commit -m "docs: held-out chr17 verdict for --gtf-refine"
```
