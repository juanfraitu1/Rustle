# Fuzzy Junction Skeleton Merging Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Merge de-novo skeletons whose intron chains match in count and are within a pre-registered
per-junction coordinate tolerance, closing part of the gap to StringTie found in
`bench/CHR20_ASSEMBLER_COMPARISON.md`, without touching the exact-match grouping function itself (zero
blast radius on the shared multi-copy family pipeline).

**Architecture:** A new pure function `merge_fuzzy_skeletons` runs AFTER `pass1_skeletons_robust` builds
its skeletons exactly as today, wired in only on the one call site that serves `--gtf`'s pure de novo path
(`detect_and_assign`, `supplied_families.is_none()`), gated behind a new opt-in env var
`RUSTLE_JUNCTION_FUZZ_BP`. The tolerance value is measured from real chr20 alignment jitter and
pre-registered in a committed doc BEFORE its effect on gffcompare is observed.

**Tech Stack:** Rust (existing `denovo_assemble.rs`/`denovo_pipeline.rs`), Python 3 + pysam (existing
`bench/*.py` convention) for the jitter-measurement script.

**Spec:** `docs/superpowers/specs/2026-09-15-fuzzy-junction-merge-design.md`

## Global Constraints

- `RUSTLE_JUNCTION_FUZZ_BP` unset (or `0`) MUST produce byte-identical `--gtf` output to the current binary.
- Never wired into the multi-copy family / O1 detection pipeline in this plan — general/`--gtf` de novo
  path only (the `detect_and_assign` call site where `supplied_families.is_none()`).
- The tolerance value is measured from real chr20 data (Task 1) and fixed BEFORE its effect on the metric
  is observed (Task 5) — never re-derived after seeing a disappointing result without saying so explicitly.
- Merged junction position = most-common EXACT value among pooled reads at that slot (never an average or
  invented coordinate); ties broken by the lowest coordinate. Merged skeleton's `chrom`/`start`/`end`/
  `read_strand` carried from whichever input skeleton has more reads, ties broken by lower `start`.
- Never merge a `footprint: true` skeleton with anything (its `introns` are uncovered read-coverage gaps,
  not real splice junctions — jitter-tolerance semantics do not apply).

---

### Task 1: Measure real chr20 junction jitter and pre-register the tolerance

**Files:**
- Create: `bench/measure_junction_jitter.py`
- Create: `docs/PREREG_junction_fuzz_2026-09-15.md`

**Interfaces:**
- Consumes: `/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20/chr20.bam`,
  `/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20/chr20_ref.gtf` (both already exist from the chr20
  bakeoff work, `bench/CHR20_ASSEMBLER_COMPARISON.md`).
- Produces: a real, printed 90th-percentile jitter value that Task 2/3's `RUSTLE_JUNCTION_FUZZ_BP` default
  documentation and Task 5's real run will use verbatim — copy the ACTUAL number this script prints into
  the PREREG doc; do not estimate it.

- [ ] **Step 1: Write the jitter-measurement script**

```python
#!/usr/bin/env python3
"""PREREG (docs/PREREG_junction_fuzz_2026-09-15.md): measure real per-junction alignment jitter on chr20,
BEFORE any fuzzy-merge code is written or evaluated against a metric. For every real chr20 read's spliced
junction that sits within CAPTURE_BP of a real annotated RefSeq intron boundary (matched by donor-site
distance), record the signed offset at both the donor and acceptor site. The 90th percentile of the
pooled |offset| becomes RUSTLE_JUNCTION_FUZZ_BP's pre-registered value.

usage: python3 bench/measure_junction_jitter.py <chr20.bam> <chr20_ref.gtf>
"""
import sys, re, collections, pysam

bam_p, gtf_p = sys.argv[1], sys.argv[2]
CAPTURE_BP = 500

# annotated introns: for each ref transcript, gaps between consecutive exons
ex = collections.defaultdict(list)
for line in open(gtf_p):
    if line.startswith('#'):
        continue
    f = line.rstrip('\n').split('\t')
    if len(f) < 9 or f[2] != 'exon':
        continue
    m = re.search(r'transcript_id "([^"]+)"', f[8])
    if not m:
        continue
    ex[m.group(1)].append((int(f[3]) - 1, int(f[4]), f[0]))

ann_introns = collections.defaultdict(set)  # chrom -> set of (start, end)
for t, v in ex.items():
    v.sort()
    for a, b in zip(v, v[1:]):
        ann_introns[a[2]].add((a[1], b[0]))
ann_sorted = {c: sorted(v) for c, v in ann_introns.items()}


def introns_of(pos, cig):
    o, p = [], pos
    for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cig):
        n = int(n)
        if op in 'M=XD':
            p += n
        elif op == 'N':
            o.append((p, p + n))
            p += n
    return o


def nearest(chrom, don):
    """closest annotated intron by donor-site distance, within CAPTURE_BP, else None."""
    best, bd = None, CAPTURE_BP + 1
    for a, b in ann_sorted.get(chrom, ()):
        d = abs(a - don)
        if d < bd:
            best, bd = (a, b), d
    return best if bd <= CAPTURE_BP else None


bam = pysam.AlignmentFile(bam_p, 'rb')
offsets = []
n_reads = n_junctions = n_matched = 0
for rec in bam:
    if rec.is_unmapped or rec.is_secondary or rec.is_supplementary or not rec.cigarstring:
        continue
    n_reads += 1
    chrom = rec.reference_name
    for don, acc in introns_of(rec.reference_start, rec.cigarstring):
        n_junctions += 1
        ref = nearest(chrom, don)
        if ref is None:
            continue
        n_matched += 1
        rd, ra = ref
        offsets.append(abs(don - rd))
        offsets.append(abs(acc - ra))

offsets.sort()
n = len(offsets)
p50 = offsets[int(n * 0.50)]
p90 = offsets[int(n * 0.90)]
p95 = offsets[int(n * 0.95)]
print(f'reads scanned: {n_reads}')
print(f'junctions seen: {n_junctions}')
print(f'junctions matched to an annotated intron within {CAPTURE_BP}bp: {n_matched}')
print(f'pooled donor+acceptor |offset| samples: {n}')
print(f'median |offset|: {p50}')
print(f'90th percentile |offset|: {p90}   <-- this is RUSTLE_JUNCTION_FUZZ_BP\'s pre-registered value')
print(f'95th percentile |offset|: {p95}')
```

- [ ] **Step 2: Run it against the real chr20 substrate**

```bash
source /home/juanfra/miniforge3/etc/profile.d/conda.sh; conda activate sqanti3  # has pysam
python3 bench/measure_junction_jitter.py \
  /mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20/chr20.bam \
  /mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20/chr20_ref.gtf \
  > /tmp/jitter_output.txt 2>&1
cat /tmp/jitter_output.txt
```

Expected: real printed numbers (reads scanned should be in the tens of thousands, matching
`bench/CHR20_ASSEMBLER_COMPARISON.md`'s own chr20 read counts). Whatever the 90th percentile prints is the
REAL, final value for Task 2/3/5 — do not adjust it.

- [ ] **Step 3: Write the PREREG doc with the real numbers from Step 2**

Follow the existing `docs/PREREG_*.md` convention (see `docs/PREREG_core_definition_2026-09-12.md` for the
header/section style: `# PREREG — <question> (date, written before any score below is computed)`,
`## Question`, then the method and the real numbers). Content:

```markdown
# PREREG — what junction-jitter tolerance should RUSTLE_JUNCTION_FUZZ_BP use? (2026-09-15, written before
any gffcompare/SQANTI3 score from this feature is computed)

## Question
On real chr20 IsoSeq reads, how far do individual reads' own CIGAR-derived splice junction coordinates
scatter around the nearest real annotated RefSeq intron boundary -- and what tolerance, chosen by ONE
fixed rule stated here, should `merge_fuzzy_skeletons` (docs/superpowers/specs/2026-09-15-fuzzy-junction-
merge-design.md) use?

## Method
`bench/measure_junction_jitter.py` against `/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20/chr20.bam`
and the same substrate's `chr20_ref.gtf`. Every spliced primary read's own junction is matched to the
nearest annotated intron by donor-site distance (capture window 500bp); the signed offset at both donor
and acceptor is pooled into one distribution.

## Result
[PASTE THE REAL OUTPUT OF STEP 2 VERBATIM HERE -- reads scanned, junctions seen, junctions matched, pooled
sample count, median, 90th percentile, 95th percentile]

## Decision (fixed now, before Task 5 runs)
`RUSTLE_JUNCTION_FUZZ_BP` = the 90th percentile value above, an integer number of base pairs. This rule
(90th percentile of real, freshly-measured jitter) was chosen in the design spec BEFORE this script ran and
is not changed based on Task 5's result.
```

- [ ] **Step 4: Commit**

```bash
git add bench/measure_junction_jitter.py docs/PREREG_junction_fuzz_2026-09-15.md
git commit -m "docs: pre-register RUSTLE_JUNCTION_FUZZ_BP from real chr20 jitter, before any score"
```

---

### Task 2: `merge_fuzzy_skeletons` pure function + unit tests

**Files:**
- Modify: `src/rustle/vg_family/denovo_assemble.rs` (append, near `pass1_skeletons_robust`)
- Test: same file's `#[cfg(test)] mod tests` block

**Interfaces:**
- Consumes: `Skeleton { chrom: String, start: u64, end: u64, n_reads: u32, introns: Vec<(u64,u64)>,
  tied_seeded: bool, read_strand: Option<char>, footprint: bool, read_rev: u32, read_tot: u32 }`
  (`denovo_assemble.rs:54-102`, confirmed exact field set this session).
- Produces: `pub fn merge_fuzzy_skeletons(skeletons: Vec<Skeleton>, tolerance_bp: u64) -> Vec<Skeleton>` —
  Task 3 calls this with the value from `RUSTLE_JUNCTION_FUZZ_BP`.

- [ ] **Step 1: Write the failing tests**

```rust
#[test]
fn fuzzy_merge_combines_skeletons_within_tolerance() {
    let mk = |start: u64, end: u64, introns: Vec<(u64, u64)>, n_reads: u32| Skeleton {
        chrom: "chr1".into(), start, end, n_reads, introns, tied_seeded: false,
        read_strand: Some('+'), footprint: false, read_rev: 0, read_tot: n_reads,
    };
    let a = mk(1000, 5000, vec![(1500, 2500), (3000, 4000)], 5);
    let b = mk(1003, 5000, vec![(1502, 2497), (3001, 3998)], 3); // every junction within 5bp of a's
    let out = merge_fuzzy_skeletons(vec![a, b], 5);
    assert_eq!(out.len(), 1, "within-tolerance skeletons must merge into one");
    assert_eq!(out[0].n_reads, 8, "reads sum across the merge");
    assert_eq!(out[0].introns, vec![(1500, 2500), (3000, 4000)], "the higher-n_reads skeleton's exact junctions win each slot (5 reads > 3)");
}

#[test]
fn fuzzy_merge_never_merges_different_intron_counts() {
    let mk = |introns: Vec<(u64, u64)>| Skeleton {
        chrom: "chr1".into(), start: 1000, end: 5000, n_reads: 5, introns, tied_seeded: false,
        read_strand: Some('+'), footprint: false, read_rev: 0, read_tot: 5,
    };
    let a = mk(vec![(1500, 2500)]);
    let b = mk(vec![(1500, 2500), (3000, 4000)]); // one more intron than a
    let out = merge_fuzzy_skeletons(vec![a, b], 1000); // huge tolerance, must still not merge
    assert_eq!(out.len(), 2, "different intron counts must never merge regardless of tolerance");
}

#[test]
fn fuzzy_merge_never_merges_beyond_tolerance() {
    let mk = |don: u64| Skeleton {
        chrom: "chr1".into(), start: 1000, end: 5000, n_reads: 5, introns: vec![(don, don + 1000)],
        tied_seeded: false, read_strand: Some('+'), footprint: false, read_rev: 0, read_tot: 5,
    };
    let a = mk(1500);
    let b = mk(1520); // 20bp away
    let out = merge_fuzzy_skeletons(vec![a, b], 5); // tolerance 5, gap is 20
    assert_eq!(out.len(), 2, "a junction beyond tolerance must never merge");
}

#[test]
fn fuzzy_merge_chains_single_linkage_through_an_intermediate() {
    let mk = |don: u64| Skeleton {
        chrom: "chr1".into(), start: 1000, end: 5000, n_reads: 1, introns: vec![(don, don + 1000)],
        tied_seeded: false, read_strand: Some('+'), footprint: false, read_rev: 0, read_tot: 1,
    };
    let a = mk(1500);
    let b = mk(1504); // within 5 of a
    let c = mk(1508); // within 5 of b, but 8 away from a (beyond tolerance 5 if compared directly)
    let out = merge_fuzzy_skeletons(vec![a, b, c], 5);
    assert_eq!(out.len(), 1, "single-linkage chaining must merge all three through b, matching cluster_tie_partners' own tested behavior");
    assert_eq!(out[0].n_reads, 3);
}

#[test]
fn fuzzy_merge_at_zero_tolerance_is_a_no_op() {
    let mk = |don: u64| Skeleton {
        chrom: "chr1".into(), start: 1000, end: 5000, n_reads: 1, introns: vec![(don, don + 1000)],
        tied_seeded: false, read_strand: Some('+'), footprint: false, read_rev: 0, read_tot: 1,
    };
    let a = mk(1500);
    let b = mk(1500); // exact duplicate, would merge at ANY tolerance >= 0 by the within_tol test
    let out = merge_fuzzy_skeletons(vec![a, b], 0);
    assert_eq!(out.len(), 2, "tolerance_bp=0 must be an explicit no-op, matching every other opt-in flag's off-state contract -- byte-identical to not calling this function at all");
}

#[test]
fn fuzzy_merge_never_merges_a_footprint_skeleton() {
    let mk = |footprint: bool| Skeleton {
        chrom: "chr1".into(), start: 1000, end: 5000, n_reads: 5, introns: vec![(1500, 2500)],
        tied_seeded: false, read_strand: Some('+'), footprint, read_rev: 0, read_tot: 5,
    };
    let a = mk(true);
    let b = mk(false); // identical introns, but a is a footprint (uncovered-gap semantics, not real junctions)
    let out = merge_fuzzy_skeletons(vec![a, b], 1000);
    assert_eq!(out.len(), 2, "a footprint skeleton's 'introns' are read-coverage gaps, not splice junctions -- never eligible for jitter-tolerance merging");
}

#[test]
fn fuzzy_merge_respects_known_strand_disagreement_but_unknown_strand_never_blocks() {
    // Resolves the design spec's open question: Skeleton::read_strand is None for a spliced skeleton
    // whose strand hasn't been determined via junction motifs yet at this stage of the pipeline, so an
    // absent signal must not block an otherwise-valid merge -- but a KNOWN disagreement must.
    let mk = |strand: Option<char>| Skeleton {
        chrom: "chr1".into(), start: 1000, end: 5000, n_reads: 5, introns: vec![(1500, 2500)],
        tied_seeded: false, read_strand: strand, footprint: false, read_rev: 0, read_tot: 5,
    };
    let plus = mk(Some('+'));
    let minus = mk(Some('-')); // identical introns, but KNOWN opposite strand
    let out = merge_fuzzy_skeletons(vec![plus, minus], 1000);
    assert_eq!(out.len(), 2, "two skeletons with KNOWN opposite strand must never merge, even with identical junctions");

    let known = mk(Some('+'));
    let unknown = mk(None); // strand not yet determined on this side
    let out2 = merge_fuzzy_skeletons(vec![known, unknown], 1000);
    assert_eq!(out2.len(), 1, "an absent strand signal on one side must not block a merge on otherwise-matching junctions");
}
```

- [ ] **Step 2: Run tests, confirm they fail to compile / fail**

Run: `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo test --release fuzzy_merge 2>&1 | tail -40`
Expected: compile error (`merge_fuzzy_skeletons` not defined).

- [ ] **Step 3: Implement `merge_fuzzy_skeletons`**

```rust
/// Merge skeletons whose intron chains are structurally identical (same intron count) and differ only by
/// <= `tolerance_bp` at every corresponding donor/acceptor position. Single-linkage: if A merges with B
/// and B merges with C, all three end up in one group even if A and C alone exceed the tolerance (same
/// chaining pattern as `cluster_tie_partners`, `src/rustle/vg_family/copy_discovery.rs`). Pure function
/// over `Skeleton` values -- no BAM/env access, unit-testable with synthetic skeletons.
///
/// `tolerance_bp == 0` is an explicit no-op (returns `skeletons` unchanged) -- every opt-in mechanism in
/// this codebase must be byte-identical to "not called" at its off value.
///
/// A `footprint: true` skeleton's `introns` are uncovered READ-COVERAGE gaps, not real splice junctions
/// (see `Skeleton::footprint`'s own doc comment) -- jitter tolerance is meaningless for it, so it is never
/// merged with anything, including another footprint.
pub fn merge_fuzzy_skeletons(skeletons: Vec<Skeleton>, tolerance_bp: u64) -> Vec<Skeleton> {
    if tolerance_bp == 0 {
        return skeletons;
    }
    fn within_tol(a: &Skeleton, b: &Skeleton, tol: u64) -> bool {
        !a.footprint
            && !b.footprint
            && a.chrom == b.chrom
            && !a.introns.is_empty()
            && a.introns.len() == b.introns.len()
            // A spliced skeleton's strand may not be determined yet at this stage (read_strand is
            // consulted only for unspliced models -- see Skeleton::read_strand's own doc comment), so an
            // absent signal on either side must not block the merge; a KNOWN disagreement must.
            && match (a.read_strand, b.read_strand) {
                (Some(sa), Some(sb)) => sa == sb,
                _ => true,
            }
            && a.introns.iter().zip(&b.introns).all(|(x, y)| {
                x.0.abs_diff(y.0) <= tol && x.1.abs_diff(y.1) <= tol
            })
    }
    let n = skeletons.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(parent: &mut [usize], x: usize) -> usize {
        if parent[x] != x {
            parent[x] = find(parent, parent[x]);
        }
        parent[x]
    }
    for i in 0..n {
        for j in (i + 1)..n {
            if within_tol(&skeletons[i], &skeletons[j], tolerance_bp) {
                let (ri, rj) = (find(&mut parent, i), find(&mut parent, j));
                if ri != rj {
                    parent[ri] = rj;
                }
            }
        }
    }
    let mut groups: std::collections::BTreeMap<usize, Vec<usize>> = std::collections::BTreeMap::new();
    for i in 0..n {
        let r = find(&mut parent, i);
        groups.entry(r).or_default().push(i);
    }
    groups
        .into_values()
        .map(|idxs| {
            if idxs.len() == 1 {
                return skeletons[idxs[0]].clone();
            }
            let members: Vec<&Skeleton> = idxs.iter().map(|&i| &skeletons[i]).collect();
            let n_introns = members[0].introns.len();
            let mut merged_introns = Vec::with_capacity(n_introns);
            for slot in 0..n_introns {
                // most-common EXACT value at this slot, weighted by each skeleton's own n_reads;
                // ties broken by the lowest coordinate (never an average/invented coordinate).
                let mut counts: std::collections::BTreeMap<(u64, u64), u32> = std::collections::BTreeMap::new();
                for m in &members {
                    *counts.entry(m.introns[slot]).or_insert(0) += m.n_reads;
                }
                let mut items: Vec<((u64, u64), u32)> = counts.into_iter().collect();
                items.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
                merged_introns.push(items[0].0);
            }
            // representative for chrom/start/end/read_strand: more reads wins, ties -> lower start.
            let mut reps: Vec<&Skeleton> = members.clone();
            reps.sort_by(|a, b| b.n_reads.cmp(&a.n_reads).then(a.start.cmp(&b.start)));
            let rep = reps[0];
            Skeleton {
                chrom: rep.chrom.clone(),
                start: rep.start,
                end: rep.end,
                n_reads: members.iter().map(|m| m.n_reads).sum(),
                introns: merged_introns,
                tied_seeded: members.iter().any(|m| m.tied_seeded),
                read_strand: rep.read_strand,
                footprint: false, // guaranteed by within_tol excluding footprints from ever reaching here
                read_rev: members.iter().map(|m| m.read_rev).sum(),
                read_tot: members.iter().map(|m| m.read_tot).sum(),
            }
        })
        .collect()
}
```

- [ ] **Step 4: Run tests, confirm pass**

Run: `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo test --release fuzzy_merge 2>&1 | tail -40`
Expected: 7/7 pass.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/denovo_assemble.rs
git commit -m "feat: add merge_fuzzy_skeletons, pure and unit-tested"
```

---

### Task 3: Wire `RUSTLE_JUNCTION_FUZZ_BP` into the pure de novo call site

**Files:**
- Modify: `src/rustle/vg_family/denovo_pipeline.rs`

**Interfaces:**
- Consumes: `merge_fuzzy_skeletons` (Task 2).
- Produces: the `skeletons` binding inside `detect_and_assign` reflects the merge when the flag is set.

- [ ] **Step 1: Locate and confirm the exact call site**

Read `detect_and_assign` (`denovo_pipeline.rs`, starts at line 2165 as of this plan's writing — confirm,
this file changes over the course of this plan) around the block:

```rust
let skeletons = if supplied {
    Vec::new()
} else {
    pass1_skeletons_robust(seed_reads, cfg.pass1_min_reads, cfg.min_terminal_support)
};
```

This is the ONLY call site of the 5 total `pass1_skeletons_robust` sites that serves the pure de novo
`--gtf` path (no `--families`) — confirmed this session by tracing `detect_and_assign`'s own comment
("with a supplied catalog the whole detection front end is DEAD WORK... skeletons is enough to switch it
all off"). Confirm this comment and the `if supplied { Vec::new() } else { ... }` shape are still present
and unchanged before editing; if the surrounding code has shifted, find the real current equivalent rather
than assuming the line number.

- [ ] **Step 2: Add the opt-in merge call**

```rust
let skeletons = if supplied {
    Vec::new()
} else {
    let sk = pass1_skeletons_robust(seed_reads, cfg.pass1_min_reads, cfg.min_terminal_support);
    // Opt-in (RUSTLE_JUNCTION_FUZZ_BP, default off): merge skeletons whose intron chains match in count
    // and differ only by a pre-registered per-junction tolerance -- docs/PREREG_junction_fuzz_2026-09-15.md,
    // docs/superpowers/specs/2026-09-15-fuzzy-junction-merge-design.md. Zero effect when unset (tolerance
    // 0 is `merge_fuzzy_skeletons`'s own explicit no-op), so every existing catalog stays byte-identical.
    let fuzz_bp: u64 = std::env::var("RUSTLE_JUNCTION_FUZZ_BP")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0);
    if fuzz_bp > 0 {
        crate::vg_family::denovo_assemble::merge_fuzzy_skeletons(sk, fuzz_bp)
    } else {
        sk
    }
};
```

- [ ] **Step 3: Build**

Run: `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo build --release --bin copy_assign 2>&1 | tail -60`
Expected: clean build.

- [ ] **Step 4: Commit**

```bash
git add src/rustle/vg_family/denovo_pipeline.rs
git commit -m "feat: wire RUSTLE_JUNCTION_FUZZ_BP into the pure de novo skeleton path"
```

---

### Task 4: Byte-identical-when-unset regression test

**Files:**
- Test: `tests/copy_assign_families.rs` (mirror the existing `discover_copies_off_by_default_is_byte_identical`
  test's structure and helpers in the same file)

- [ ] **Step 1: Write the failing test**

```rust
/// `RUSTLE_JUNCTION_FUZZ_BP` unset and explicitly `0` are the same "off" state -- both must produce
/// byte-identical output. This is the directly-testable form of "unset is a no-op": the fixture this file
/// uses (`same_chrom_supplement`) is small enough that a real merge may or may not trigger at any given
/// tolerance, but unset-vs-zero must ALWAYS agree regardless, since both resolve to `tolerance_bp == 0`
/// (`merge_fuzzy_skeletons`'s own proven no-op, Task 2) before the merge function is ever called.
#[test]
fn junction_fuzz_unset_and_explicit_zero_are_byte_identical() {
    let d_unset = scratch("junction_fuzz_unset");
    let (o_unset, out_unset) = run(&d_unset, &["--no-refine", "--gtf"]);
    assert!(o_unset.status.success(), "unset run failed:\n{}", stderr(&o_unset));

    let d_zero = scratch("junction_fuzz_zero");
    let out_s_zero = d_zero.join("o").to_str().expect("utf-8 path").to_string();
    let o_zero = std::process::Command::new(env!("CARGO_BIN_EXE_copy_assign"))
        .args(["--bam", &format!("{FIX}/reads.bam"), "--fasta", &format!("{FIX}/genome.fa")])
        .args(["--region", "c1:200-600", "--out", &out_s_zero])
        .args(["--no-refine", "--gtf"])
        .env("RUSTLE_JUNCTION_FUZZ_BP", "0")
        .output()
        .expect("copy_assign failed to spawn");
    assert!(o_zero.status.success(), "explicit-zero run failed:\n{}", String::from_utf8_lossy(&o_zero.stderr));

    for ext in ["assignments.tsv", "families.tsv", "quant.tsv", "gtf"] {
        let a = std::fs::read(format!("{out_unset}.{ext}")).unwrap_or_else(|e| panic!("read {out_unset}.{ext}: {e}"));
        let b = std::fs::read(format!("{out_s_zero}.{ext}")).unwrap_or_else(|e| panic!("read {out_s_zero}.{ext}: {e}"));
        assert_eq!(a, b, "RUSTLE_JUNCTION_FUZZ_BP unset vs explicit 0 must be byte-identical for .{ext}");
    }
}
```

- [ ] **Step 2: Run, confirm it fails to compile / fails for the right reason before Task 3's commit; passes after**

Run: `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo test --release --test copy_assign_families junction_fuzz 2>&1 | tail -40`

- [ ] **Step 3: Commit**

```bash
git add tests/copy_assign_families.rs
git commit -m "test: RUSTLE_JUNCTION_FUZZ_BP unset is a no-op"
```

---

### Task 5: Real-data acceptance test using the pre-registered tolerance

**Files:** none (a verification run) — unless a real, non-cosmetic bug surfaces, in which case return to
Task 2/3 for a real code fix, not a tuning pass.

- [ ] **Step 1: Rebuild**

```bash
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo build --release --bin copy_assign > /tmp/task5_build.log 2>&1
tail -20 /tmp/task5_build.log
```

- [ ] **Step 2: Rerun the chr20 "ours" arm with the PRE-REGISTERED tolerance from Task 1**

```bash
W=/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20
BIN=/mnt/linuxdisk/home/juanfraitu/rustle_target/release/copy_assign
mkdir -p "$W/ours_fuzzy"
cd "$W/ours_fuzzy"
RUSTLE_JUNCTION_FUZZ_BP=<VALUE FROM docs/PREREG_junction_fuzz_2026-09-15.md — DO NOT GUESS, READ THE FILE> \
  "$BIN" --gtf --bam "$W/chr20.bam" --fasta "$W/chr20.fa" --region chr20:1-66210255 --out ours \
  > ours.stdout.log 2> ours.stderr.log
echo exit=$?
grep -c $'\ttranscript\t' ours.gtf
```

- [ ] **Step 3: Rescore with gffcompare, compare to the existing baseline**

```bash
cd "$W"
gffcompare -r chr20_ref.gtf -o /tmp/verify_fuzzy ours_fuzzy/ours.gtf
cat /tmp/verify_fuzzy.stats | sed -n '1,20p'
echo "=== baseline for comparison ==="
cat gffcompare/ours.stats | sed -n '1,20p'
```

- [ ] **Step 4: Log the result honestly in `bench/CHR20_ASSEMBLER_COMPARISON.md`**

Append a new dated section (`## Follow-up: fuzzy junction merging at the pre-registered tolerance
(2026-09-15, [POSITIVE/NEGATIVE/MIXED])`) with: the exact tolerance used (from Task 1's PREREG doc), the
number of transcripts before/after merging, the full gffcompare stats table for both baseline and the
fuzzy-merge run side by side (matching the TES-snap follow-up section's own table format, already in this
file), and an honest one-paragraph verdict — whichever way it goes. Do NOT adjust `RUSTLE_JUNCTION_FUZZ_BP`
away from the pre-registered value to chase a better number; if the pre-registered value doesn't help,
that is the answer, reported plainly, matching this project's own standing culture (§6l6's negative
result, the TES-snap follow-up's negative result).

- [ ] **Step 5: Commit**

```bash
git add bench/CHR20_ASSEMBLER_COMPARISON.md
git commit -m "docs: fuzzy junction merge acceptance test at the pre-registered tolerance"
```
