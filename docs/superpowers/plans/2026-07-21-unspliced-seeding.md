# Position-Aware Unspliced-Read Seeding Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stop `pass1_skeletons_robust` from pooling all unspliced reads on a chromosome into one over-length skeleton (rejected by `MAX_SPLICED`), by clustering the empty-intron-chain reads position-aware so each pseudogene/retrocopy copy seeds its own locus — restoring the reads the 28 "distinguishable-but-merged" Soto members lost at seeding.

**Architecture:** A new pure function `cluster_unspliced` does single-linkage overlap clustering of unspliced reads (threshold-free, mirroring `thin_loci`). `pass1_skeletons_robust` splits its input: spliced reads keep the current `(chrom, introns)` keying byte-for-byte; unspliced reads route through `cluster_unspliced`. Downstream χ(H)/assignment is unchanged and still governs which seeded members resolve vs abstain.

**Tech Stack:** Rust (`src/rustle/vg_family/denovo_assemble.rs`), `cargo`, Python3 (Soto eval + floor decompose), real gorilla data (`A119b.t2t.bam` + `chm13v2.0.fa` on `/mnt/linuxdisk/`).

## Global Constraints

- **Result-changing — validated on Soto, not byte-identity.** Success = the 28 merged cases now seed at their loci; Soto member sensitivity **> 76.2%**; recovery precision **= 100%** (every newly-seeded locus overlaps a real annotated Soto member — no noise splits); the 36 K=0 members produce **no false copies**; **spliced-only regions byte-identical**.
- **Threshold-free clustering:** single-linkage on pure span overlap. No new distance/gap constant. The only gate is the existing `min_reads` floor per cluster.
- **Spliced path untouched:** the `(chrom, introns)` keying for non-empty intron chains must be byte-for-byte the current code (the existing tests at `denovo_assemble.rs:1056,1058` use spliced reads and must still pass unchanged).
- **`pass1_skeletons_robust` signature unchanged** (`reads, min_reads, min_terminal_support`) — its 5 callers (`denovo_pipeline.rs:183,1396,1920,2119,2260`) must not change.
- **Crash rule (WSL2):** any `copy_assign` run is FOREGROUND, serial, region-restricted, output to `/home/juanfra/winloci_scratch`. NO nohup/background `&`/waiter loops/pkill.
- Commit on the current branch (`dna-family-fallback`); messages end with `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.

## File Structure

- `src/rustle/vg_family/denovo_assemble.rs` — `cluster_unspliced` (new pure fn + tests) and the split inside `pass1_skeletons_robust`.
- `bench/soto/merge_baseline.txt` — reused from the prior effort (the frozen before-numbers + K=0 family list).

Struct facts (verbatim): `Skeleton { chrom: String, start: u64, end: u64, n_reads: u32, introns: Vec<(u64,u64)> }`; `PrimaryRead` has `chrom: String`, `introns: Vec<(u64,u64)>`, `ref_start: u64`, `ref_end: u64`.

---

## Task 1: `cluster_unspliced` — single-linkage overlap clustering (pure fn, TDD)

**Files:**
- Modify: `src/rustle/vg_family/denovo_assemble.rs` (add the fn + `#[cfg(test)]` cases)

**Interfaces:**
- Produces: `pub fn cluster_unspliced(reads: &[PrimaryRead], min_reads: u32, k: usize) -> Vec<Skeleton>` — clusters ONLY the empty-intron-chain reads by single-linkage span overlap, per chromosome; each cluster with ≥ `min_reads` reads becomes a `Skeleton { introns: vec![], start/end = per-cluster robust boundaries }`. Spliced reads in the input are ignored (they are handled by the caller's spliced path).

- [ ] **Step 1: Write the failing tests**

Construct `PrimaryRead` values the same way the neighbouring tests in `denovo_assemble.rs` do (find the existing read-builder helper via `grep -n "fn read\|PrimaryRead {" src/rustle/vg_family/denovo_assemble.rs`); each test read here has `introns: vec![]` (unspliced).

```rust
// in denovo_assemble.rs #[cfg(test)] mod tests
#[test]
fn cluster_unspliced_separates_distant_loci() {
    // two unspliced read groups 50kb apart -> TWO skeletons, not one chromosome-spanning giant
    let reads = vec![
        unspliced("c1", 1000, 3000), unspliced("c1", 1100, 3100), unspliced("c1", 1200, 3200),
        unspliced("c1", 51000, 53000), unspliced("c1", 51100, 53100), unspliced("c1", 51200, 53200),
    ];
    let sk = cluster_unspliced(&reads, 3, 1);
    assert_eq!(sk.len(), 2, "distant unspliced loci must seed as separate skeletons");
    assert!(sk.iter().all(|s| s.end - s.start < 300_000), "each skeleton spans one locus, not the chromosome");
    assert!(sk.iter().all(|s| s.introns.is_empty()));
}
#[test]
fn cluster_unspliced_merges_overlapping_and_filters_min_reads() {
    // one overlapping pile -> ONE skeleton; a lone read below min_reads -> dropped
    let reads = vec![
        unspliced("c1", 100, 400), unspliced("c1", 150, 450), unspliced("c1", 200, 500),
        unspliced("c1", 90000, 90200), // singleton, below min_reads=2
    ];
    let sk = cluster_unspliced(&reads, 2, 1);
    assert_eq!(sk.len(), 1, "overlapping unspliced reads merge; the singleton is below min_reads");
    assert_eq!(sk[0].n_reads, 3);
}
```
Add a tiny `fn unspliced(chrom: &str, s: u64, e: u64) -> PrimaryRead` test helper next to them if none exists, building a `PrimaryRead` with `introns: vec![]`, `ref_start: s`, `ref_end: e`, and whatever other fields the struct requires at their zero/empty defaults (match the existing read-builder).

- [ ] **Step 2: Run to verify failure**

Run: `cd /mnt/c/Users/jfris/Desktop/Rustle && cargo test --lib cluster_unspliced 2>&1 | tail -5`
Expected: FAIL (function not defined).

- [ ] **Step 3: Implement `cluster_unspliced`**

```rust
/// Position-aware seeding for UNSPLICED (empty-intron-chain) reads: single-linkage span-overlap
/// clustering per chromosome, so each pseudogene/retrocopy copy seeds its OWN skeleton instead of
/// every unspliced read on a chromosome pooling into one `(chrom, [])` group that overruns
/// MAX_SPLICED. Threshold-free (pure overlap), mirroring `thin_loci`'s single-linkage-by-span.
pub fn cluster_unspliced(reads: &[PrimaryRead], min_reads: u32, k: usize) -> Vec<Skeleton> {
    use std::collections::BTreeMap;
    let mut by_chrom: BTreeMap<&str, Vec<&PrimaryRead>> = BTreeMap::new();
    for r in reads {
        if r.introns.is_empty() {
            by_chrom.entry(r.chrom.as_str()).or_default().push(r);
        }
    }
    let mut out = Vec::new();
    for (chrom, mut rs) in by_chrom {
        rs.sort_by_key(|r| (r.ref_start, r.ref_end));
        let mut i = 0;
        while i < rs.len() {
            // single-linkage: extend the cluster while the next read starts before the running max end
            let mut cluster_end = rs[i].ref_end;
            let mut j = i + 1;
            while j < rs.len() && rs[j].ref_start < cluster_end {
                cluster_end = cluster_end.max(rs[j].ref_end);
                j += 1;
            }
            let cluster = &rs[i..j];
            let n = cluster.len() as u32;
            if n >= min_reads {
                // per-cluster robust boundaries: k-th smallest start, k-th largest end (matches
                // pass1_skeletons_robust's boundary rule).
                let mut starts: Vec<u64> = cluster.iter().map(|r| r.ref_start).collect();
                let mut ends: Vec<u64> = cluster.iter().map(|r| r.ref_end).collect();
                starts.sort_unstable();
                ends.sort_unstable_by(|a, b| b.cmp(a));
                let si = k.min(starts.len()).saturating_sub(1);
                let ei = k.min(ends.len()).saturating_sub(1);
                out.push(Skeleton {
                    chrom: chrom.to_string(),
                    start: starts[si],
                    end: ends[ei],
                    n_reads: n,
                    introns: Vec::new(),
                });
            }
            i = j;
        }
    }
    out
}
```

- [ ] **Step 4: Run to verify pass**

Run: `cd /mnt/c/Users/jfris/Desktop/Rustle && cargo test --lib cluster_unspliced 2>&1 | tail -5`
Expected: PASS (2 tests).

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/denovo_assemble.rs
git commit -m "feat(seeding): cluster_unspliced — single-linkage overlap seeding of unspliced reads (TDD)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: Split `pass1_skeletons_robust` — spliced untouched, unspliced routed

**Files:**
- Modify: `src/rustle/vg_family/denovo_assemble.rs` (`pass1_skeletons_robust` body)

**Interfaces:**
- Consumes: `cluster_unspliced` (Task 1).
- Produces: `pass1_skeletons_robust` unchanged signature; spliced reads keep the `(chrom, introns)` keying; unspliced reads (empty chain) are seeded via `cluster_unspliced` and appended.

- [ ] **Step 1: Write the failing test (unspliced reads seed separately)**

```rust
// denovo_assemble.rs tests — the WHOLE-FUNCTION behavior
#[test]
fn pass1_seeds_unspliced_reads_position_aware_not_one_giant() {
    // spliced reads (one intron chain) + two DISTANT unspliced piles on the same chrom
    let mut reads = vec![
        // spliced locus (unchanged behavior)
        read_spliced("c1", 1000, 3000, &[(1500, 2000)]), read_spliced("c1", 1010, 3010, &[(1500, 2000)]),
        // unspliced pile A
        unspliced("c1", 400000, 402000), unspliced("c1", 400100, 402100), unspliced("c1", 400200, 402200),
        // unspliced pile B, 2 Mb away
        unspliced("c1", 2400000, 2402000), unspliced("c1", 2400100, 2402100), unspliced("c1", 2400200, 2402200),
    ];
    let sk = pass1_skeletons_robust(&reads, 2, 1);
    // 1 spliced skeleton + 2 distinct unspliced skeletons = 3; NOT a single giant unspliced one
    let unspliced_sk: Vec<_> = sk.iter().filter(|s| s.introns.is_empty()).collect();
    assert_eq!(unspliced_sk.len(), 2, "the two distant unspliced piles seed separately");
    assert!(unspliced_sk.iter().all(|s| s.end - s.start < 300_000), "no chromosome-spanning unspliced skeleton");
}
```
(Use the existing read-builder for spliced reads — `read_spliced`/whatever the file already has — and the `unspliced` helper from Task 1.)

- [ ] **Step 2: Run to verify failure**

Run: `cd /mnt/c/Users/jfris/Desktop/Rustle && cargo test --lib pass1_seeds_unspliced 2>&1 | tail -6`
Expected: FAIL — the current code pools both unspliced piles into one giant `(c1, [])` skeleton (unspliced_sk.len() == 1, and/or spans > 2 Mb).

- [ ] **Step 3: Split the function**

In `pass1_skeletons_robust`, keep the existing `groups` BTreeMap loop but **skip unspliced reads in it**, then append `cluster_unspliced` for those:

```rust
pub fn pass1_skeletons_robust(reads: &[PrimaryRead], min_reads: u32, min_terminal_support: u32) -> Vec<Skeleton> {
    use std::collections::BTreeMap;
    let k = min_terminal_support.max(1) as usize;
    let mut groups: BTreeMap<(&str, Vec<(u64, u64)>), (u32, Vec<u64>, Vec<u64>)> = BTreeMap::new();
    for r in reads {
        if r.introns.is_empty() {
            continue; // unspliced reads are seeded position-aware below (empty chain would pool chromosome-wide)
        }
        let e = groups
            .entry((r.chrom.as_str(), r.introns.clone()))
            .or_insert((0, Vec::new(), Vec::new()));
        e.0 += 1;
        let pos = e.1.partition_point(|&x| x <= r.ref_start);
        if pos < k { e.1.insert(pos, r.ref_start); e.1.truncate(k); }
        let pos = e.2.partition_point(|&x| x >= r.ref_end);
        if pos < k { e.2.insert(pos, r.ref_end); e.2.truncate(k); }
    }
    let mut skels: Vec<Skeleton> = groups
        .into_iter()
        .filter(|(_, (n, _, _))| *n >= min_reads)
        .map(|((chrom, introns), (n, starts, ends))| {
            let si = k.min(starts.len()).saturating_sub(1);
            let ei = k.min(ends.len()).saturating_sub(1);
            Skeleton { chrom: chrom.to_string(), start: starts[si], end: ends[ei], n_reads: n, introns }
        })
        .collect();
    // position-aware seeding of the unspliced reads (the fix)
    skels.extend(cluster_unspliced(reads, min_reads, k));
    skels
}
```
Note: adding the `if r.introns.is_empty() { continue; }` guard means the spliced branch is byte-for-byte the old code MINUS the unspliced reads it used to (mis-)pool. The existing tests at `:1056,1058` use spliced reads → unaffected.

- [ ] **Step 4: Build + run the whole-function + existing tests**

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
cargo build --release 2>&1 | tail -3
cargo test --lib pass1 2>&1 | tail -8
cargo test --lib 2>&1 | tail -4
```
Expected: builds; the new `pass1_seeds_unspliced...` test PASSES; the existing `:1056,1058` pass1 tests still PASS (they use spliced reads); full suite green except the 1 known pre-existing deleted-fixture failure. If an existing pass1 test that happened to use unspliced reads now changes, inspect it — a pre-existing test asserting the OLD (buggy) giant-skeleton behavior should be updated to the correct per-locus behavior and noted.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/denovo_assemble.rs
git commit -m "fix(seeding): pass1_skeletons_robust seeds unspliced reads position-aware

Unspliced (empty-intron-chain) reads no longer pool into one (chrom,[]) skeleton that
overruns MAX_SPLICED and gets rejected; they route through cluster_unspliced (single-linkage
overlap) so each pseudogene/retrocopy copy seeds its own locus. Spliced path unchanged.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: Validate on Soto (28 seed, precision held, K=0 no false copies, spliced byte-identical)

**Files:** none created; a verification-only task (may append a results file).

**Interfaces:**
- Consumes: `bench/soto/merge_baseline.txt` (frozen before-numbers + the K=0 family list, from the prior effort). OLD binary = commit `fedf8f0` (current HEAD before the Task 1/2 seeding commits); NEW binary = HEAD after Task 2.

- [ ] **Step 1: Old-vs-new on the 28 merged loci (they must now SEED)**

Build BOTH versions and run the same regions (foreground, crash rule). From `bench/soto/soto_floor_decomposition.tsv`, the `MERGED`-cause rows with the highest `unique_own_reads` (col 8): ID_26 SLC9B1P1 chr16:32804124-32821138 (40), ID_261 CNTNAP3P1 chr9:80247804-80257862 (24), ID_260 MEP1AP4 chr9:44369523-44375633 (11).
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
# NEW (current HEAD, seeding fix in): build + run the 3 loci
cargo build --release 2>&1 | tail -2
for reg in chr16:32754124-32871138 chr9:80197804-80307862 chr9:44319523-44425633; do
  cd /home/juanfra/winloci_scratch
  /mnt/c/Users/jfris/Desktop/Rustle/target/release/copy_assign --bam /mnt/linuxdisk/home/juanfraitu/winloci_data/A119b.t2t.bam \
    --fasta /mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa --region "$reg" --homology-primary --min-copies 2 \
    --out "new_$(echo $reg|tr ':' '_')" ; cd /mnt/c/Users/jfris/Desktop/Rustle
done
# OLD (fedf8f0): stash, checkout, build, run same regions
git stash; git checkout fedf8f0 -- src/; cargo build --release 2>&1 | tail -2
for reg in chr16:32754124-32871138 chr9:80197804-80307862 chr9:44319523-44425633; do
  cd /home/juanfra/winloci_scratch
  /mnt/c/Users/jfris/Desktop/Rustle/target/release/copy_assign --bam /mnt/linuxdisk/home/juanfraitu/winloci_data/A119b.t2t.bam \
    --fasta /mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa --region "$reg" --homology-primary --min-copies 2 \
    --out "old_$(echo $reg|tr ':' '_')" ; cd /mnt/c/Users/jfris/Desktop/Rustle
done
git checkout dna-family-fallback -- src/; git stash pop 2>/dev/null; cargo build --release 2>&1 | tail -2
```
Compare per region: does a family/copy now appear at the member's coordinates (chr16:32804124-32821138 etc.) in `new_*.families.tsv` that was absent in `old_*.families.tsv`? Expected: **NEW seeds the member's locus** (a copy at/overlapping the member's coords) where OLD did not.

- [ ] **Step 2: Precision + K=0 check on the same runs**

For each NEW locus that appeared, confirm it overlaps a real annotated Soto member (`bench/soto/80_fams.chr.bed`). And run 3 K=0 (`expressed-K=0`) loci from `merge_baseline.txt` the same old-vs-new way: they must NOT gain a spurious copy (a K=0 member may now be SEEDED as a locus but must not produce a *false* extra copy beyond the annotated members).
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
# for each NEW copy coordinate, check overlap with an annotated member:
awk -F'\t' '{split($4,a,"|"); print $1, $2, $3, a[1]}' bench/soto/80_fams.chr.bed | head  # reference of real members
```
Expected: every new locus overlaps an annotated member (precision held); no K=0 locus emits a copy that does not correspond to a real member.

- [ ] **Step 3: Spliced-region byte-identity**

Pick a spliced Soto family (one whose members are multi-exon — GSTM works: `NC_073224.2` is gorilla; for CHM13 use a spliced family region from the BED). Run NEW vs OLD on it and confirm `families.tsv` md5 is IDENTICAL (the spliced path is unchanged).
```bash
cd /home/juanfra/winloci_scratch
# after building NEW and OLD as above, compare a spliced family region's families.tsv md5
md5sum new_<spliced_region>.families.tsv old_<spliced_region>.families.tsv
```
Expected: identical md5 (spliced path untouched). If different, the split leaked into the spliced path — investigate.

- [ ] **Step 4: Re-run the floor decomposition**

Run:
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
python3 bench/soto/soto_floor_decompose.py 2>&1 | grep -E "MERGED|distinguishable|information floor|unseeded"
```
Expected: the "distinguishable-but-MERGED" bucket shrinks vs `merge_baseline.txt`; the "unseeded" bucket likely also shrinks (same root cause). Record the before/after. NOTE: the floor decompose reads the detection catalog; for the genome-wide number it must run against a catalog built with the new binary — if that full regeneration is out of reach here, report the per-region old-vs-new evidence (Steps 1–2) as the primary result and note the genome-wide number needs a background catalog rebuild.

- [ ] **Step 5: Write the results + commit**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
{ echo "# unspliced-seeding validation (old fedf8f0 vs new HEAD)";
  echo "## 28-case seeding (per-region old-vs-new): <fill: which of the 3 now seed>";
  echo "## precision: <fill: new loci all overlap annotated members? Y/N>";
  echo "## K=0: <fill: no false copies? Y/N>";
  echo "## spliced byte-identity: <fill: md5 match Y/N>";
  echo "## floor decomposition: MERGED <before>-><after>"; } > bench/soto/seeding_validation.txt
git add bench/soto/seeding_validation.txt
git commit -m "test(seeding): Soto old-vs-new validation — unspliced members now seed at their loci

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

**If precision drops (a new locus does NOT map to a real member) or a K=0 case emits a false copy:** STOP and report — that is the over-split failure mode; the `min_reads` floor or the clustering needs adjustment before this is accepted. Do not paper over a precision regression.

---

## Self-Review Notes (for the planner)

- **Spec coverage:** §1 the split → Tasks 1+2. §2 downstream (χ(H) unchanged) → nothing to build (existing machinery); asserted by the K=0 check in Task 3. §3 components → Task 1 (`cluster_unspliced`) + Task 2 (the split). §4 validation → Task 3 (all 6 criteria: 28 seed, precision, K=0, spliced byte-identity, floor decompose, unit tests in Tasks 1–2). §6 success criteria → Task 3. All covered.
- **Placeholder scan:** Task 3's `seeding_validation.txt` uses `<fill: …>` for the measured outcomes the implementer records at runtime — these are result placeholders (the numbers don't exist until the runs happen), not spec placeholders; every command that produces them is concrete. The spliced-region choice in Step 3 is left to the implementer to pick a multi-exon family from the BED (concrete criterion given). No logic/threshold is left unnamed.
- **Type consistency:** `cluster_unspliced(reads: &[PrimaryRead], min_reads: u32, k: usize) -> Vec<Skeleton>` defined in Task 1 and called identically in Task 2. `Skeleton`/`PrimaryRead` fields match the verbatim struct facts. `pass1_skeletons_robust` signature unchanged.
- **No new constant:** clustering is pure overlap; the only gate is the pre-existing `min_reads`.
