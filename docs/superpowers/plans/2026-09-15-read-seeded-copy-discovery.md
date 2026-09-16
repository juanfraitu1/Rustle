# Read-Seeded Copy Discovery Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give `copy_assign --families` a `--discover-copies` flag that clusters AS-tied reads' out-of-catalog
placements into candidate new copies, reported to a TSV (never auto-merged), closing the §6l5 gap.

**Architecture:** New pure, unit-tested clustering module (`src/rustle/vg_family/copy_discovery.rs`), wired
into `src/bin/copy_assign.rs`'s existing per-region closure using the exact same opt-in-flag-gated pattern
already used for `--flag-missing-copies` (O3), then drained and written the same way `o3_raw_pairs`/
`o3_orphan_loci` already are.

**Tech Stack:** Rust (existing `rustle` crate), no new dependencies.

**Spec:** `docs/superpowers/specs/2026-09-15-read-seeded-copy-discovery-design.md`

## Global Constraints

- `--discover-copies` defaults to `false`; unset, output must be BYTE-IDENTICAL to the current binary
  (verify with a real regression run, Task 6).
- Never mutates the input catalog or the current run's assignments — report file only, two-pass workflow.
- Admission bar: `n_supporting_reads >= 2` (reuses the existing PSV read-support convention, do not
  invent a different constant).
- Cluster merge distance: `500` bp, a named constant (`TIE_PARTNER_MERGE_DISTANCE_BP`), not a magic number
  inline.
- No genome-only discovery: every candidate must trace to a real AS-tied read's own real BAM placement.

---

### Task 1: `DiscoveredCopy` + pure clustering primitive, with unit tests

**Files:**
- Create: `src/rustle/vg_family/copy_discovery.rs`
- Modify: `src/rustle/vg_family/mod.rs` (add `pub mod copy_discovery;` — find the existing `pub mod
  read_conflict;` or similar line and add a new line next to it, alphabetically if the file is already
  sorted that way, otherwise at the end of the `pub mod` block)
- Test: same file, `#[cfg(test)] mod tests` block at the bottom (project convention — every other module
  in `vg_family/` keeps its tests inline; see `read_conflict.rs`'s own `mod tests` for the pattern)

**Interfaces:**
- Produces: `pub struct DiscoveredCopy { family_id: String, chrom: String, start: u64, end: u64,
  n_supporting_reads: usize, read_names: Vec<String>, nearest_copy_tid: String, nearest_copy_distance: u64
  }` and `pub fn cluster_tie_partners(tied_reads: &[(String, Vec<(String, u64, u64)>)], family_id: &str,
  existing_copies: &[(String, u64, u64, String)], merge_distance: u64, min_support: usize) ->
  Vec<DiscoveredCopy>` — both consumed by Task 3. `existing_copies` is `(chrom, start, end, tid)` per
  already-catalogued copy in this family — deliberately a plain tuple slice, not tied to any one upstream
  struct type (`CatalogFamily`'s `CatalogCopy` and `ColocatedFamily`'s `DenovoTranscript` both have these
  same four fields under different names; Task 3 confirmed `ColocatedFamily`/`DenovoTranscript`
  — `src/rustle/vg_family/denovo_pipeline.rs:281-287`, `:65-77` — is what's actually in scope at the call
  site, so this function must not hard-couple to `catalog_input::CatalogFamily`).
- Consumes: nothing from `catalog_input` or `denovo_pipeline` directly — kept decoupled on purpose (see
  above).

- [ ] **Step 1: Write the failing tests**

```rust
#[cfg(test)]
mod tests {
    use super::*;

    fn one_copy(tid: &str, chrom: &str, start: u64, end: u64) -> Vec<(String, u64, u64, String)> {
        vec![(chrom.to_string(), start, end, tid.to_string())]
    }

    #[test]
    fn defensively_excludes_positions_inside_a_catalog_copy() {
        let existing = one_copy("c0", "chr1", 1000, 2000);
        // this "tied" position sits INSIDE c0's span -- must never surface as a discovery
        let tied = vec![("read1".to_string(), vec![("chr1".to_string(), 1200, 1300)])];
        let out = cluster_tie_partners(&tied, "FAM0", &existing, 500, 2);
        assert!(out.is_empty(), "a position already inside a catalog copy must never be reported");
    }

    #[test]
    fn merges_positions_within_merge_distance_and_respects_min_support() {
        let existing = one_copy("c0", "chr1", 1000, 2000);
        let tied = vec![
            ("read1".to_string(), vec![("chr1".to_string(), 5000, 5100)]),
            ("read2".to_string(), vec![("chr1".to_string(), 5050, 5150)]), // within 500bp of read1's site
        ];
        let out = cluster_tie_partners(&tied, "FAM0", &existing, 500, 2);
        assert_eq!(out.len(), 1, "two nearby out-of-catalog positions with 2 supporting reads = 1 cluster");
        assert_eq!(out[0].n_supporting_reads, 2);
        assert_eq!(out[0].nearest_copy_tid, "c0");
        assert_eq!(out[0].nearest_copy_distance, 3000); // 5000 - 2000
    }

    #[test]
    fn keeps_clusters_separate_beyond_merge_distance() {
        let existing = one_copy("c0", "chr1", 1000, 2000);
        let tied = vec![
            ("read1".to_string(), vec![("chr1".to_string(), 5000, 5100)]),
            ("read2".to_string(), vec![("chr1".to_string(), 5100, 5100)]),
            ("read3".to_string(), vec![("chr1".to_string(), 9000, 9100)]),
            ("read4".to_string(), vec![("chr1".to_string(), 9050, 9150)]),
        ];
        let out = cluster_tie_partners(&tied, "FAM0", &existing, 500, 2);
        assert_eq!(out.len(), 2, "two far-apart pairs must stay two separate clusters");
    }

    #[test]
    fn drops_clusters_below_min_support() {
        let existing = one_copy("c0", "chr1", 1000, 2000);
        let tied = vec![("read1".to_string(), vec![("chr1".to_string(), 5000, 5100)])];
        let out = cluster_tie_partners(&tied, "FAM0", &existing, 500, 2);
        assert!(out.is_empty(), "a single supporting read must not clear min_support=2");
    }
}
```

- [ ] **Step 2: Run the tests, confirm they fail to compile** (the module/function don't exist yet)

Run: `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo test --release
copy_discovery 2>&1 | tail -40`
Expected: compile error, `cluster_tie_partners`/`DiscoveredCopy` not found.

- [ ] **Step 3: Implement `cluster_tie_partners`**

```rust
pub const TIE_PARTNER_MERGE_DISTANCE_BP: u64 = 500;
pub const TIE_PARTNER_MIN_SUPPORT: usize = 2;

#[derive(Clone, Debug, PartialEq)]
pub struct DiscoveredCopy {
    pub family_id: String,
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub n_supporting_reads: usize,
    pub read_names: Vec<String>,
    pub nearest_copy_tid: String,
    pub nearest_copy_distance: u64,
}

fn inside_any_copy(existing_copies: &[(String, u64, u64, String)], chrom: &str, start: u64, end: u64) -> bool {
    existing_copies.iter().any(|(c_chrom, c_start, c_end, _)| c_chrom == chrom && start < *c_end && end > *c_start)
}

fn nearest_copy(existing_copies: &[(String, u64, u64, String)], chrom: &str, start: u64, end: u64) -> (String, u64) {
    existing_copies
        .iter()
        .filter(|(c_chrom, ..)| c_chrom == chrom)
        .map(|(_, c_start, c_end, tid)| {
            let d = if end <= *c_start { c_start - end } else if start >= *c_end { start - c_end } else { 0 };
            (tid.clone(), d)
        })
        .min_by_key(|(_, d)| *d)
        .unwrap_or(("NA".to_string(), u64::MAX))
}

/// Cluster out-of-catalog AS-tie-partner positions into candidate copies for one family.
/// `tied_reads`: (read_name, all of that read's max-AS placements as (chrom, start, end)).
/// `existing_copies`: (chrom, start, end, tid) for every copy already catalogued in this family (caller
/// converts from whatever its own copy-list type is -- `ColocatedFamily.copies: Vec<DenovoTranscript>` at
/// the real call site, Task 3).
/// Positions already inside an existing copy span are defensively re-excluded here (never surfaced),
/// even though the caller (Task 3) is expected to have filtered them out already.
pub fn cluster_tie_partners(
    tied_reads: &[(String, Vec<(String, u64, u64)>)],
    family_id: &str,
    existing_copies: &[(String, u64, u64, String)],
    merge_distance: u64,
    min_support: usize,
) -> Vec<DiscoveredCopy> {
    // Flatten to (read_name, chrom, start, end), keeping only positions outside every existing copy.
    let mut sites: Vec<(String, String, u64, u64)> = Vec::new();
    for (name, placements) in tied_reads {
        for (chrom, start, end) in placements {
            if !inside_any_copy(existing_copies, chrom, *start, *end) {
                sites.push((name.clone(), chrom.clone(), *start, *end));
            }
        }
    }
    // Sort by (chrom, start) so overlap/proximity clustering is a single linear pass.
    sites.sort_by(|a, b| (a.1.clone(), a.2).cmp(&(b.1.clone(), b.2)));

    let mut clusters: Vec<(String, u64, u64, Vec<String>)> = Vec::new(); // (chrom, start, end, read_names)
    for (name, chrom, start, end) in sites {
        if let Some(last) = clusters.last_mut() {
            let (lchrom, _lstart, lend, names) = last;
            if *lchrom == chrom && start <= *lend + merge_distance {
                *lend = (*lend).max(end);
                if !names.contains(&name) {
                    names.push(name);
                }
                continue;
            }
        }
        clusters.push((chrom, start, end, vec![name]));
    }

    clusters
        .into_iter()
        .filter(|(_, _, _, names)| names.len() >= min_support)
        .map(|(chrom, start, end, read_names)| {
            let (nearest_copy_tid, nearest_copy_distance) = nearest_copy(existing_copies, &chrom, start, end);
            DiscoveredCopy {
                family_id: family_id.to_string(),
                chrom,
                start,
                end,
                n_supporting_reads: read_names.len(),
                read_names,
                nearest_copy_tid,
                nearest_copy_distance,
            }
        })
        .collect()
}
```

- [ ] **Step 4: Run tests, confirm they pass**

Run: `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo test --release
copy_discovery 2>&1 | tail -40`
Expected: 4 passed, 0 failed.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/copy_discovery.rs src/rustle/vg_family/mod.rs
git commit -m "feat: add copy_discovery::cluster_tie_partners, pure and unit-tested"
```

---

### Task 2: Extract tied reads' full placement sets from `bam_reads`

**Files:**
- Modify: `src/rustle/vg_family/copy_discovery.rs` (append; keep the pure clustering primitive from Task 1
  untouched)
- Test: same file's `mod tests`

**Interfaces:**
- Consumes: `crate::rustle::vg_family::denovo_assemble::BamRead` — confirmed fields `chrom: String, read:
  AlignedRead, mapq: u8, name: String, as_score: i32, de: f32, is_supplementary: bool, is_secondary: bool,
  reverse: bool, ts: Option<char>` (`denovo_assemble.rs:950-969`); `AlignedRead { ref_start: u64, cigar:
  Vec<(char,u64)>, seq: Vec<u8>, qual: Vec<u8> }` (`copy_split.rs:177-182`). Reference-end-from-CIGAR is
  already implemented in this exact shape as `read_ref_end_local` at `src/bin/copy_assign.rs:1333-1336`
  (`ref_start + sum of M/=/X/D/N op lengths`) — it is private to the `bin` target, so this new `lib` module
  cannot import it; the `ref_end` helper below is that same logic, duplicated on purpose (a `bin`-private
  fn is never visible to a `lib` module in this crate layout).
- Produces: `pub fn tie_partner_placements(bam_reads: &[BamRead]) -> Vec<(String, Vec<(String, u64,
  u64)>)]` — one entry per read name that has >=2 placements tied at that read's own maximum `as_score`
  among its non-supplementary records (a read with only one placement can't be tied with anything and is
  excluded). This is the direct input to Task 1's `cluster_tie_partners`, plumbed together in Task 3.

- [ ] **Step 1: Write the failing test**

```rust
#[test]
fn tie_partner_placements_finds_reads_tied_at_their_own_max_as() {
    use crate::rustle::vg_family::denovo_assemble::{AlignedRead, BamRead};
    let mk = |name: &str, chrom: &str, start: u64, as_score: i32| BamRead {
        chrom: chrom.into(),
        read: AlignedRead { ref_start: start, cigar: vec![('M', 100)], seq: vec![], qual: vec![] },
        mapq: 0, name: name.into(), as_score, de: 0.0,
        is_supplementary: false, is_secondary: as_score != 200, reverse: false, ts: None,
    };
    let reads = vec![
        mk("tied_read", "chr1", 1000, 200),   // best
        mk("tied_read", "chr1", 5000, 200),   // tied with the above
        mk("tied_read", "chr1", 9000, 150),   // worse, not part of the tie
        mk("unique_read", "chr1", 2000, 300), // only one placement, never tied
    ];
    let out = tie_partner_placements(&reads);
    assert_eq!(out.len(), 1);
    assert_eq!(out[0].0, "tied_read");
    assert_eq!(out[0].1.len(), 2, "only the 2 max-scoring placements, not the 150-scoring one");
}
```

(`AlignedRead`'s fields are confirmed exactly `{ref_start: u64, cigar: Vec<(char,u64)>, seq: Vec<u8>, qual:
Vec<u8>}` — `src/rustle/vg_family/copy_split.rs:177-182` — the literals above are correct as written.)

- [ ] **Step 2: Run test, confirm it fails to compile / fails**

- [ ] **Step 3: Implement `tie_partner_placements`**

```rust
use std::collections::HashMap;
use crate::rustle::vg_family::denovo_assemble::BamRead;

fn ref_end(br: &BamRead) -> u64 {
    br.read.ref_start
        + br.read.cigar.iter().filter(|(op, _)| matches!(op, 'M' | '=' | 'X' | 'D' | 'N')).map(|(_, n)| *n).sum::<u64>()
}

pub fn tie_partner_placements(bam_reads: &[BamRead]) -> Vec<(String, Vec<(String, u64, u64)>)> {
    let mut by_name: HashMap<&str, Vec<&BamRead>> = HashMap::new();
    for br in bam_reads.iter().filter(|b| !b.is_supplementary) {
        by_name.entry(br.name.as_str()).or_default().push(br);
    }
    let mut out = Vec::new();
    for (name, placements) in by_name {
        if placements.len() < 2 {
            continue;
        }
        let max_as = placements.iter().map(|b| b.as_score).max().unwrap();
        let tied: Vec<&&BamRead> = placements.iter().filter(|b| b.as_score == max_as).collect();
        if tied.len() >= 2 {
            out.push((
                name.to_string(),
                tied.into_iter().map(|b| (b.chrom.clone(), b.read.ref_start, ref_end(b))).collect(),
            ));
        }
    }
    out
}
```

- [ ] **Step 4: Run tests, confirm pass**

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/copy_discovery.rs
git commit -m "feat: extract per-read AS-tied placement sets for copy discovery"
```

---

### Task 3: Wire `--discover-copies` into `copy_assign`'s per-region pass

**Files:**
- Modify: `src/bin/copy_assign.rs`

**Interfaces:**
- Consumes: `copy_discovery::{tie_partner_placements, cluster_tie_partners, DiscoveredCopy,
  TIE_PARTNER_MERGE_DISTANCE_BP, TIE_PARTNER_MIN_SUPPORT}` (Tasks 1-2); the region's existing `bam_reads:
  Vec<BamRead>` and `fams: Vec<FamilyAssignment>` (each with a `.family_id`/copies list reachable the same
  way the O3 block at `copy_assign.rs:2609-2634` already reaches it — read that block in full before
  writing this task, and mirror its exact "gate the whole computation behind the flag, produce an empty
  Vec when off" shape); the region's own `CatalogFamily` list (however `load_supplied_families` /
  `RegionFamilies` exposes it to this closure — grep `fams` and `bound` in the surrounding ~200 lines to
  find the right handle; do not guess a name that isn't actually in scope here).
- Produces: a new `discovered: Vec<DiscoveredCopy>` field on `RegionWork` (mirroring `o3_raw_pairs`/
  `o3_orphan_loci`, which are exactly this shape already — add the field next to them in the struct
  definition at `copy_assign.rs:72` and in both places `RegionWork { .. }` is constructed/destructured,
  lines 2843 and 2887 as of this plan's writing — re-check those line numbers still match before editing,
  this file changes over the course of this same plan).

- [ ] **Step 1: Add the CLI flag**

In the `Args` struct (near the other opt-in bool flags, e.g. right after `dump_psv` at
`copy_assign.rs:337-343`):

```rust
/// Cluster AS-tied reads' out-of-catalog placements into candidate new copies, written to
/// `<out>.discovered_copies.tsv`. Report only -- never mutates the input catalog or this run's own
/// assignments (two-pass: inspect the report, append accepted rows to the catalog by hand, re-run).
/// Default off; unset, output is byte-identical to a run without this flag.
#[arg(long, default_value_t = false)]
discover_copies: bool,
```

- [ ] **Step 2: Add the field to `RegionWork`**

At `copy_assign.rs:72` (the `struct RegionWork` definition), add:
```rust
discovered: Vec<crate::rustle::vg_family::copy_discovery::DiscoveredCopy>,
```

- [ ] **Step 3: Compute it in the per-region closure, gated on the flag**

Traced this session: `fams: Vec<FamilyAssignment>` (the O3 block's own population, `copy_assign.rs:2609`
onward) is `detect_and_assign`'s OUTPUT and only carries `family_id`, not a copy-span list. The copy spans
this task needs are the INPUT to that same call, already in scope a few lines above it:
`supplied: Option<Vec<ColocatedFamily>>` (`copy_assign.rs:2387`), where `ColocatedFamily { family_id,
chrom, start, end, copies: Vec<DenovoTranscript> }` (`denovo_pipeline.rs:281-287`) and `DenovoTranscript`
has `tid: String, chrom: String, start: u64, end: u64` (`family_detect.rs:65-69`) — exactly the four
fields `cluster_tie_partners` needs, joined to each `fams` entry by matching `family_id`.

Immediately after the O3 block that ends around `copy_assign.rs:2634`-ish (find its actual closing brace;
this plan's line numbers are approximate once Task 3's own edits start landing), add a sibling block:

```rust
// Read-seeded copy discovery (opt-in, --discover-copies): cluster AS-tied reads' out-of-catalog
// placements into candidate new copies. Gated the same way as the O3 block above -- empty Vec, no
// allocation, when the flag is unset.
let discovered: Vec<_> = if args.discover_copies {
    let tied = crate::rustle::vg_family::copy_discovery::tie_partner_placements(&bam_reads);
    let empty: Vec<ColocatedFamily> = Vec::new();
    let colocated = supplied.as_ref().unwrap_or(&empty);
    fams.iter()
        .flat_map(|fa| {
            let existing_copies: Vec<(String, u64, u64, String)> = colocated
                .iter()
                .find(|cf| cf.family_id == fa.family_id)
                .map(|cf| cf.copies.iter().map(|c| (c.chrom.clone(), c.start, c.end, c.tid.clone())).collect())
                .unwrap_or_default();
            crate::rustle::vg_family::copy_discovery::cluster_tie_partners(
                &tied,
                &fa.family_id,
                &existing_copies,
                crate::rustle::vg_family::copy_discovery::TIE_PARTNER_MERGE_DISTANCE_BP,
                crate::rustle::vg_family::copy_discovery::TIE_PARTNER_MIN_SUPPORT,
            )
        })
        .collect()
} else {
    Vec::new()
};
```

Confirm `supplied` is still named/typed exactly this way and still in scope at this point in the
function before pasting this in — it was read directly from the live file this session, but Tasks 1-2's
own edits land first and this file may have shifted by the time this task starts.

- [ ] **Step 4: Add `discovered` to both `RegionWork { .. }` sites**

At the construction site (`copy_assign.rs:2843` as of this plan) add `discovered,` to the struct literal;
at the destructuring site (`copy_assign.rs:2887`) add `discovered` to the pattern.

- [ ] **Step 5: Build and fix type errors**

Run: `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo build --release --bin
copy_assign 2>&1 | tee /tmp/build.log | tail -60`
Fix whatever the real field/accessor names turn out to be (this task's Step 3 flagged the one genuinely
uncertain line). Do not proceed until this builds clean.

- [ ] **Step 6: Commit**

```bash
git add src/bin/copy_assign.rs
git commit -m "feat: wire --discover-copies into the per-region assignment pass"
```

---

### Task 4: Write `<out>.discovered_copies.tsv`

**Files:**
- Modify: `src/bin/copy_assign.rs`

**Interfaces:**
- Consumes: the drained `discovered: Vec<DiscoveredCopy>` across all regions -- find where
  `o3_raw_pairs`/`o3_orphan_loci` get collected from each region's `RegionWork` after the parallel/serial
  region loop finishes (grep `o3_raw_pairs` for the accumulation site; it is downstream of line 2887 in
  this same function) and mirror that exact accumulation for `discovered`.
- Produces: `<out>.discovered_copies.tsv`, columns `family_id\tchrom\tstart\tend\tn_supporting_reads\tread_names\tnearest_copy_tid\tnearest_copy_distance`
  (`read_names` comma-joined), written only when `args.discover_copies` is true, following the exact
  writer style already used for e.g. `<out>.psv_copies.tsv` (`copy_assign.rs:4606-4609`ish -- `File::create`
  + `writeln!` header + one `writeln!` per row + an `eprintln!` summary line matching this binary's own
  convention, e.g. `"[copy_assign] --discover-copies: {n} candidate cop{ies/y} -> {out}.discovered_copies.tsv"`).

- [ ] **Step 1: Accumulate `discovered` across all regions**

Wherever the per-region `Vec<RegionWork>` is flattened into whole-run vectors (find the accumulation for
`o3_raw_pairs`/`o3_orphan_loci` and add a sibling `let all_discovered: Vec<_> =
work_results.iter().flat_map(|w| w.discovered.clone()).collect();` -- match whatever the real variable
names are at that site, do not invent new ones that don't match the surrounding code's naming).

- [ ] **Step 2: Write the file**

```rust
if args.discover_copies {
    let mut dh = std::fs::File::create(format!("{}.discovered_copies.tsv", args.out))?;
    writeln!(dh, "family_id\tchrom\tstart\tend\tn_supporting_reads\tread_names\tnearest_copy_tid\tnearest_copy_distance")?;
    for d in &all_discovered {
        writeln!(
            dh, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            d.family_id, d.chrom, d.start, d.end, d.n_supporting_reads,
            d.read_names.join(","), d.nearest_copy_tid, d.nearest_copy_distance
        )?;
    }
    eprintln!(
        "[copy_assign] --discover-copies: {} candidate cop{} -> {}.discovered_copies.tsv",
        all_discovered.len(), if all_discovered.len() == 1 { "y" } else { "ies" }, args.out
    );
}
```

- [ ] **Step 3: Build**

Run: `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo build --release --bin
copy_assign 2>&1 | tail -40`
Expected: clean build.

- [ ] **Step 4: Commit**

```bash
git add src/bin/copy_assign.rs
git commit -m "feat: write <out>.discovered_copies.tsv"
```

---

### Task 5: Regression -- flag off is byte-identical

**Files:**
- Test: `tests/copy_assign_families.rs` (or wherever the existing `--families` regression fixtures live --
  grep for an existing small-fixture test in that file and mirror its setup rather than building a new
  fixture from scratch)

- [ ] **Step 1: Add a test that runs the existing fixture with and without `--discover-copies` and diffs
  every OTHER output file byte-for-byte** (not `.discovered_copies.tsv`, which won't exist in the without
  case)

```rust
#[test]
fn discover_copies_off_by_default_is_byte_identical() {
    // Reuse whatever small real/synthetic --families fixture this file's other tests already use --
    // do not invent a new one. Run copy_assign twice, once with no extra flags, once with
    // --discover-copies added, into two different --out prefixes, then assert
    // std::fs::read(prefix_a + ".assignments.tsv") == std::fs::read(prefix_b + ".assignments.tsv")
    // (and .families.tsv, .quant.tsv) -- discover_copies must never perturb these.
}
```

- [ ] **Step 2: Run the full existing test suite, confirm no regressions**

Run: `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo test --release 2>&1 | tail -60`
Expected: same pass/fail counts as the pre-existing baseline (854 lib passed / 2 failed, unrelated) plus
the new tests from Tasks 1, 2, and this one, all passing.

- [ ] **Step 3: Commit**

```bash
git add tests/
git commit -m "test: --discover-copies is a no-op on existing output when unset"
```

---

### Task 6: Real-data acceptance test (manual, not `cargo test`)

**Files:** none (a verification run, not a code change) -- unless it fails, in which case return to Task 3
to adjust `TIE_PARTNER_MERGE_DISTANCE_BP`/`TIE_PARTNER_MIN_SUPPORT` first, a real code fix second.

- [ ] **Step 1: Run on the real §6l5 substrate**

```bash
BIN=/mnt/linuxdisk/home/juanfraitu/rustle_target/release/copy_assign
$BIN --families /mnt/linuxdisk/home/juanfraitu/mec/psv.copies.tsv \
     --bam /mnt/linuxdisk/home/juanfraitu/npip_cat/npip3.bam \
     --fasta /mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta \
     --out /mnt/linuxdisk/home/juanfraitu/o3_probe_verify/discover_copies_test/run \
     --discover-copies
```
(Confirm `mec/psv.copies.tsv` is the right catalog file name -- this session's own investigations
referenced both `mec/psv.*` and `npip_cat/arm_f2/cat.copies.tsv` for what may or may not be the same
catalog; check before trusting either path, per the spec's own open note.)

- [ ] **Step 2: Check for the known candidate**

```bash
grep GWFAM55 /mnt/linuxdisk/home/juanfraitu/o3_probe_verify/discover_copies_test/run.discovered_copies.tsv
```
Expected: a row near `NC_073242.2:21,674,468` (within ~500bp) with `n_supporting_reads >= 2`.

- [ ] **Step 3: Log the result in the ledger**

Whether it passes or fails, append a new `§6l6` section to `docs/o1_ledger.md` (check the last `## §6l`
section number first, it may have moved past `§6l5` if other work landed in between) stating the real
result: did the known candidate appear, at what distance, with what support count; if not, what was tried
to fix it and whether that worked. Follow this document's own established convention (real numbers, exact
file paths, no rounding).

- [ ] **Step 4: Commit the ledger update**

```bash
git add docs/o1_ledger.md
git commit -m "docs: §6l6 read-seeded copy discovery acceptance test result"
```
