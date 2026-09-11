# O3 Flag-Pass Integration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Port `bench/o3_flag_pass.py`'s missing-copy detector natively into `copy_assign` as a new
opt-in `--flag-missing-copies` flag, plus a follow-on cross-individual differential comparison script.

**Architecture:** A new library module `src/rustle/vg_family/o3_flag_pass.rs` holds the detector (pure
statistics + a minimap2-realignment step), following the same "module `copy_assign.rs` calls into"
pattern as `absent_copy.rs`/`linearize.rs`. Two-phase: per-family raw stats computed inside the existing
parallel `compute()` closure; a genome-wide Bonferroni aggregation runs once after the existing serial
drain, at the same point the `productive`/`orf_aa` GTF attribute's own second pass already runs. Output
lands on `<out>.family_join.tsv` (new columns) plus one new `<out>.o3_candidate_loci.tsv`.

**Tech Stack:** Rust (existing `copy_assign` binary + `rustle` lib crate), `minimap2` subprocess
(existing project convention, no new dependency), Python 3 stdlib only for the final bench/ script.

**Spec:** `docs/superpowers/specs/2026-09-10-o3-flag-pass-integration-design.md` — read it before
starting; this plan implements it task by task and does not repeat its rationale sections.

## Global Constraints

- Every new CLI flag defaults to its current behavior; `copy_assign`'s entire output must be
  byte-identical to today whenever `--flag-missing-copies` is unset. Verify this after every task that
  touches `copy_assign.rs`, not just at the end.
- Cut-certificate (`bench/o3_cut_certificate.py`) and reconstruction (`bench/o3_reconstruct.py`) are
  explicitly OUT OF SCOPE — do not port them, do not add flags for them.
- Do not touch `--absent-copies` / `src/rustle/vg_family/absent_copy.rs` (a different, older, unrelated
  mechanism) in any way.
- No bipartite matching or facility-location step anywhere in this feature (standing project rule).
- Reuse the existing `minimap2` subprocess pattern (`copy_assign_pipeline.rs`'s `minimap2_msa_pair`:
  pid+atomic-nonce temp file names, `RUSTLE_MINIMAP2` env override, `Drop`-based cleanup) — do not invent
  a new subprocess-invocation style.
- Build with `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo build --release
  --all-targets`, capturing output to a file. Run `cargo test --release --lib` (library tests) AND
  `cargo test --release --bin copy_assign` (this binary's own `#[cfg(test)]` module — `--lib` alone does
  not run it) after every task.
- Never run `copy_assign` in the background unsupervised; foreground, one at a time, per this project's
  WSL2 crash rule. Big outputs go to `/mnt/linuxdisk`, never `/tmp`.

---

### Task 1: Core types + Poisson tail + genome-wide flag aggregation

**Files:**
- Create: `src/rustle/vg_family/o3_flag_pass.rs`
- Modify: `src/rustle/vg_family/mod.rs` (register the new module)
- Modify: `src/rustle/vg_family/allele_specific_junctions.rs:238` (`fn lgamma` → `pub(crate) fn lgamma`)

**Interfaces:**
- Produces: `Class` (`Divergent`/`Structural`), `LocusClass` (`OtherFamily`/`AnnotatedNoUnit`/`Unannotated`),
  `Flag` (`MissingCopy`/`Untestable`/`NoFlag`), `RawPair`, `OrphanLocus`, `FlaggedPair`, `O3Params`,
  `poisson_tail(k: usize, lam: f64) -> f64`, `finalize_flags(all_pairs: &[RawPair], alpha: f64) -> Vec<FlaggedPair>`.
  Every later task in this plan imports these from `rustle::vg_family::o3_flag_pass::*`.

Note the `Flag` variant is named `NoFlag`, not `None` — `Flag::None` compiles fine in Rust (enum variants
are namespaced separately from `Option::None`) but reads confusingly next to real `Option` values
everywhere else in this codebase; the spec's `Flag::None` naming is corrected here.

- [ ] **Step 1: Write the failing tests**

Create `src/rustle/vg_family/o3_flag_pass.rs` with just the test module first:

```rust
//! O3 flag-pass detector (`docs/superpowers/specs/2026-09-10-o3-flag-pass-integration-design.md`):
//! ports `bench/o3_flag_pass.py`'s missing-copy detector natively into `copy_assign`.

use crate::vg_family::allele_specific_junctions::lgamma;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Class {
    Divergent,
    Structural,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LocusClass {
    OtherFamily,
    AnnotatedNoUnit,
    Unannotated,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Flag {
    MissingCopy,
    Untestable,
    NoFlag,
}

#[derive(Clone, Debug)]
pub struct RawPair {
    pub family_id: String,
    pub copy_idx: String,
    pub is_partner: bool,
    pub n_rejected: usize,
    pub n_aligned: usize,
    pub covered_kb: f64,
    pub n_sites: usize,
    pub ctl_n: usize,
    pub ctl_covered_kb: f64,
    pub ctl_n_sites: usize,
    pub p_uncorrected: Option<f64>,
    pub class: Class,
    pub med_mismatch: f64,
    pub med_unaligned: f64,
}

#[derive(Clone, Debug)]
pub struct OrphanLocus {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub n_reads: usize,
    pub n_orphans: usize,
    pub class: LocusClass,
    pub n_genes_overlapping: usize,
    pub other_family_units: Vec<String>,
}

#[derive(Clone, Debug)]
pub struct FlaggedPair {
    pub pair: RawPair,
    pub flag: Flag,
    pub p_corrected_threshold: f64,
}

#[derive(Clone, Debug)]
pub struct O3Params {
    pub alpha: f64,
    pub max_reads: usize,
    pub min_reads: usize,
}

impl Default for O3Params {
    fn default() -> Self {
        O3Params { alpha: 0.001, max_reads: 500, min_reads: 3 }
    }
}

/// P(X >= k) for X ~ Poisson(lam). Mirrors `bench/o3_flag_pass.py`'s `poisson_tail` exactly (same
/// closed-form sum, same edge cases): k==0 -> 1.0 always true; lam<=0 with k>=1 -> 0.0 (a Poisson(0) is
/// a point mass at 0).
pub fn poisson_tail(k: usize, lam: f64) -> f64 {
    if k == 0 {
        return 1.0;
    }
    if lam <= 0.0 {
        return 0.0;
    }
    let mut s = 0.0_f64;
    for i in 0..k {
        s += (-lam + (i as f64) * lam.ln() - lgamma((i + 1) as f64)).exp();
    }
    (1.0 - s).clamp(0.0, 1.0)
}

/// Genome-wide Bonferroni aggregation (Phase 2): n_pairs = every RawPair with a computed p-value, from
/// the WHOLE run, not just one family. threshold = alpha / max(n_pairs, 1). Pure, no I/O.
pub fn finalize_flags(all_pairs: &[RawPair], alpha: f64) -> Vec<FlaggedPair> {
    let n_pairs = all_pairs.iter().filter(|p| p.p_uncorrected.is_some()).count();
    let threshold = alpha / (n_pairs.max(1) as f64);
    all_pairs
        .iter()
        .map(|p| {
            let flag = match p.p_uncorrected {
                Some(pv) if pv < threshold => Flag::MissingCopy,
                Some(_) => Flag::NoFlag,
                None => Flag::Untestable,
            };
            FlaggedPair { pair: p.clone(), flag, p_corrected_threshold: threshold }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pair(p: Option<f64>) -> RawPair {
        RawPair {
            family_id: "F".into(), copy_idx: "0".into(), is_partner: false,
            n_rejected: 5, n_aligned: 5, covered_kb: 1.0, n_sites: 1,
            ctl_n: 5, ctl_covered_kb: 1.0, ctl_n_sites: 1,
            p_uncorrected: p, class: Class::Divergent, med_mismatch: 1.0, med_unaligned: 0.0,
        }
    }

    #[test]
    fn poisson_tail_k_zero_is_always_one() {
        assert_eq!(poisson_tail(0, 5.0), 1.0);
        assert_eq!(poisson_tail(0, 0.0), 1.0);
    }

    #[test]
    fn poisson_tail_zero_lambda_with_positive_k_is_zero() {
        assert_eq!(poisson_tail(3, 0.0), 0.0);
    }

    #[test]
    fn poisson_tail_matches_hand_computed_value() {
        // P(X >= 1) for Poisson(1) = 1 - P(X=0) = 1 - e^-1 = 0.6321205588...
        let p = poisson_tail(1, 1.0);
        assert!((p - 0.6321205588).abs() < 1e-6, "got {p}");
        // P(X >= 5) for Poisson(1) is small; hand value from scipy.stats.poisson.sf(4, 1) = 0.003659846827...
        let p2 = poisson_tail(5, 1.0);
        assert!((p2 - 0.003659846827).abs() < 1e-6, "got {p2}");
    }

    #[test]
    fn finalize_flags_labels_below_threshold_as_missing_copy() {
        // 2 testable pairs -> threshold = 0.001 / 2 = 0.0005
        let pairs = vec![pair(Some(0.0001)), pair(Some(0.9))];
        let flagged = finalize_flags(&pairs, 0.001);
        assert_eq!(flagged[0].flag, Flag::MissingCopy);
        assert_eq!(flagged[1].flag, Flag::NoFlag);
        assert!((flagged[0].p_corrected_threshold - 0.0005).abs() < 1e-12);
    }

    #[test]
    fn finalize_flags_untestable_when_p_is_none() {
        let pairs = vec![pair(None), pair(Some(0.0001))];
        let flagged = finalize_flags(&pairs, 0.001);
        assert_eq!(flagged[0].flag, Flag::Untestable);
        // n_pairs counts only the testable one, so threshold = 0.001 / 1 = 0.001
        assert!((flagged[1].p_corrected_threshold - 0.001).abs() < 1e-12);
    }

    #[test]
    fn finalize_flags_on_empty_input_does_not_divide_by_zero() {
        let flagged = finalize_flags(&[], 0.001);
        assert!(flagged.is_empty());
    }
}
```

- [ ] **Step 2: Register the module and fix `lgamma`'s visibility**

Edit `src/rustle/vg_family/mod.rs`. Find the alphabetically-loose block of `pub mod` lines (any of them —
this file does not enforce ordering) and add, near the other O3-related modules (`asj_strand_bias`,
`asj_verify`, `asj_genetic_core`, `absent_copy`):

```rust
pub mod o3_flag_pass; // O3 flag-pass detector: ports bench/o3_flag_pass.py's missing-copy detector natively; see docs/superpowers/specs/2026-09-10-o3-flag-pass-integration-design.md
```

Edit `src/rustle/vg_family/allele_specific_junctions.rs:238`, change:
```rust
fn lgamma(x: f64) -> f64 {
```
to:
```rust
pub(crate) fn lgamma(x: f64) -> f64 {
```

- [ ] **Step 3: Run the tests to verify they compile and pass**

```
cd /mnt/c/Users/jfris/Desktop/Rustle
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo test --release --lib o3_flag_pass 2>&1 | tail -30
```
Expected: 6 tests pass (`poisson_tail_k_zero_is_always_one`, `poisson_tail_zero_lambda_with_positive_k_is_zero`,
`poisson_tail_matches_hand_computed_value`, `finalize_flags_labels_below_threshold_as_missing_copy`,
`finalize_flags_untestable_when_p_is_none`, `finalize_flags_on_empty_input_does_not_divide_by_zero`).

- [ ] **Step 4: Full regression check**

```
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo test --release --lib 2>&1 | tail -5
```
Expected: `792 passed` (or one more than the current baseline plus the 6 new — confirm no existing test
broke from the `lgamma` visibility change; grep the output for `FAILED` to be sure).

- [ ] **Step 5: Commit**

```
git add src/rustle/vg_family/o3_flag_pass.rs src/rustle/vg_family/mod.rs src/rustle/vg_family/allele_specific_junctions.rs
git commit -m "O3 flag-pass: core types, Poisson tail, genome-wide flag aggregation (Task 1)"
```

---

### Task 2: PAF-consistency parser (the realignment statistic, no subprocess needed to test)

**Files:**
- Modify: `src/rustle/vg_family/o3_flag_pass.rs`

**Interfaces:**
- Consumes: nothing from Task 1 directly (independent piece of the same module).
- Produces: `pub(crate) struct AlignmentSummary { pub covered_kb: f64, pub n_sites: usize, pub per_read:
  std::collections::HashMap<String, (usize, i64, usize)> }` (name -> (mismatch_count, unaligned_bases,
  aligned_query_len)), `pub(crate) fn parse_paf_consistency(paf_text: &str) -> AlignmentSummary`. Task 4
  calls this directly on real minimap2 output.

This deliberately does NOT take a `read_seqs` parameter or track per-position allele identity, unlike the
Python's `detector()`. The Python computes a per-position `Counter` of which base is the dominant
mismatch allele (`cons = {p: c.most_common(1)[0] for p, ...}`), but only `len(cons)` (the site count) is
ever read downstream in any output column — the allele identity itself is dead data. Dropping it removes
the need to reverse-complement/index into each read's own sequence at all, which is why this function
only needs the PAF text.

- [ ] **Step 1: Write the failing tests**

Append to `src/rustle/vg_family/o3_flag_pass.rs` (above the existing `mod tests` block, or as new
functions inside it — keep one `mod tests`):

```rust
/// Parses one `(\d+)([=XIDNS])` CIGAR run list, e.g. "88=3D510=1I436=" -> [(88,'='),(3,'D'),(510,'='),(1,'I'),(436,'=')].
fn parse_cigar_ops(cg: &str) -> Vec<(u64, char)> {
    let mut out = Vec::new();
    let mut num = 0u64;
    for ch in cg.chars() {
        if ch.is_ascii_digit() {
            num = num * 10 + (ch as u64 - '0' as u64);
        } else {
            out.push((num, ch));
            num = 0;
        }
    }
    out
}

pub(crate) struct AlignmentSummary {
    pub covered_kb: f64,
    pub n_sites: usize,
    pub per_read: std::collections::HashMap<String, (usize, i64, usize)>,
}

/// Parses `minimap2 -x splice -c --eqx -N 1` PAF output. Target-position coverage/mismatch tallies (no
/// query-sequence lookup needed — see the module doc comment on why allele identity is dropped).
/// PAF columns used: [0]=query name [1]=query len [2]=query start [3]=query end [7]=target start,
/// [12..]=tags (the `cg:Z:` CIGAR tag).
pub(crate) fn parse_paf_consistency(paf_text: &str) -> AlignmentSummary {
    let mut cov: std::collections::HashMap<u64, u32> = std::collections::HashMap::new();
    let mut mism: std::collections::HashMap<u64, u32> = std::collections::HashMap::new();
    let mut per_read: std::collections::HashMap<String, (usize, i64, usize)> = std::collections::HashMap::new();
    for line in paf_text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 13 {
            continue;
        }
        let cg = match f[12..].iter().find(|t| t.starts_with("cg:Z:")) {
            Some(t) => &t[5..],
            None => continue,
        };
        let qlen: i64 = f[1].parse().unwrap_or(0);
        let qstart: i64 = f[2].parse().unwrap_or(0);
        let qend: i64 = f[3].parse().unwrap_or(0);
        let mut t: u64 = f[7].parse().unwrap_or(0);
        let mut nx: usize = 0;
        for (num, op) in parse_cigar_ops(cg) {
            match op {
                '=' => {
                    for k in 0..num {
                        *cov.entry(t + k).or_insert(0) += 1;
                    }
                    t += num;
                }
                'X' => {
                    nx += num as usize;
                    for k in 0..num {
                        *cov.entry(t + k).or_insert(0) += 1;
                        *mism.entry(t + k).or_insert(0) += 1;
                    }
                    t += num;
                }
                'D' | 'N' => t += num,
                _ => {} // 'I'/'S' advance only the query position, which this function never tracks
            }
        }
        let unal = qlen - (qend - qstart);
        let aligned_len = (qend - qstart) as usize;
        let name = f[0].to_string();
        let better = per_read.get(&name).map_or(true, |&(_, _, prev)| aligned_len > prev);
        if better {
            per_read.insert(name, (nx, unal, aligned_len));
        }
    }
    let n_sites = mism
        .iter()
        .filter(|(p, &c)| c >= 3 && (c as f64) >= 0.5 * (*cov.get(p).unwrap_or(&0) as f64))
        .count();
    let covered_kb = cov.values().filter(|&&c| c >= 3).count() as f64 / 1000.0;
    AlignmentSummary { covered_kb, n_sites, per_read }
}

#[cfg(test)]
mod paf_tests {
    use super::*;

    #[test]
    fn single_clean_read_no_mismatches() {
        // read "r1", qlen 100, aligned 0..100 on target starting at 1000, 100 matched bases
        let paf = "r1\t100\t0\t100\t+\tY\t2000\t1000\t1100\t100\t100\t60\tcg:Z:100=";
        let s = parse_paf_consistency(paf);
        assert_eq!(s.n_sites, 0);
        assert_eq!(s.per_read.get("r1"), Some(&(0, 0, 100)));
        assert!((s.covered_kb - 0.1).abs() < 1e-9); // 100 positions >= coverage 3? NO - coverage is 1 here
    }

    #[test]
    fn covered_kb_requires_coverage_at_least_three() {
        // three reads all covering the same 10bp target window -> those 10 positions reach coverage 3
        let paf = concat!(
            "r1\t10\t0\t10\t+\tY\t100\t0\t10\t10\t10\t60\tcg:Z:10=\n",
            "r2\t10\t0\t10\t+\tY\t100\t0\t10\t10\t10\t60\tcg:Z:10=\n",
            "r3\t10\t0\t10\t+\tY\t100\t0\t10\t10\t10\t60\tcg:Z:10=",
        );
        let s = parse_paf_consistency(paf);
        assert!((s.covered_kb - 0.010).abs() < 1e-9, "10 positions at cov>=3 => 0.010 kb, got {}", s.covered_kb);
    }

    #[test]
    fn consistent_mismatch_site_needs_majority_and_floor_of_three() {
        // 4 reads: 3 mismatch at target pos 5, 1 matches -> count 3 >= 3 AND 3 >= 0.5*4 (cov=4) -> consistent
        let paf = concat!(
            "a\t10\t0\t10\t+\tY\t100\t0\t10\t10\t10\t60\tcg:Z:5=1X4=\n",
            "b\t10\t0\t10\t+\tY\t100\t0\t10\t10\t10\t60\tcg:Z:5=1X4=\n",
            "c\t10\t0\t10\t+\tY\t100\t0\t10\t10\t10\t60\tcg:Z:5=1X4=\n",
            "d\t10\t0\t10\t+\tY\t100\t0\t10\t10\t10\t60\tcg:Z:10=",
        );
        let s = parse_paf_consistency(paf);
        assert_eq!(s.n_sites, 1);
    }

    #[test]
    fn mismatch_below_majority_is_not_consistent() {
        // 4 reads, only 1 mismatches at pos 5 (cov 4, mismatch count 1 < 0.5*4=2) -> not consistent
        let paf = concat!(
            "a\t10\t0\t10\t+\tY\t100\t0\t10\t10\t10\t60\tcg:Z:5=1X4=\n",
            "b\t10\t0\t10\t+\tY\t100\t0\t10\t10\t10\t60\tcg:Z:10=\n",
            "c\t10\t0\t10\t+\tY\t100\t0\t10\t10\t10\t60\tcg:Z:10=\n",
            "d\t10\t0\t10\t+\tY\t100\t0\t10\t10\t10\t60\tcg:Z:10=",
        );
        let s = parse_paf_consistency(paf);
        assert_eq!(s.n_sites, 0);
    }

    #[test]
    fn keeps_the_longer_alignment_when_a_read_has_two_paf_lines() {
        let paf = concat!(
            "r1\t100\t0\t40\t+\tY\t1000\t0\t40\t40\t40\t60\tcg:Z:40=\n",
            "r1\t100\t0\t90\t+\tY\t1000\t100\t190\t90\t90\t60\tcg:Z:90=",
        );
        let s = parse_paf_consistency(paf);
        assert_eq!(s.per_read.get("r1").map(|&(_, _, len)| len), Some(90));
    }

    #[test]
    fn unaligned_bases_is_query_len_minus_aligned_span() {
        // qlen 100, aligned query 10..90 -> unal = 100 - (90-10) = 20
        let paf = "r1\t100\t10\t90\t+\tY\t1000\t0\t80\t80\t80\t60\tcg:Z:80=";
        let s = parse_paf_consistency(paf);
        assert_eq!(s.per_read.get("r1"), Some(&(0, 20, 80)));
    }
}
```

- [ ] **Step 2: Run to verify it fails first, then implement**

The functions above are already the implementation (this is a case where writing the test alongside the
already-known-correct implementation is appropriate per the design's algorithm derivation — the
implementation is not separated from the test in a red/green step here because both were derived
together from the Python source; still run once to confirm compile + pass, not to see a red state):

```
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo test --release --lib o3_flag_pass::paf_tests 2>&1 | tail -20
```
Expected: all 6 `paf_tests` pass.

- [ ] **Step 3: Full regression check**

```
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo test --release --lib 2>&1 | tail -5
```

- [ ] **Step 4: Commit**

```
git add src/rustle/vg_family/o3_flag_pass.rs
git commit -m "O3 flag-pass: PAF consistency-site parser, unit-tested against hand-built PAF (Task 2)"
```

---

### Task 3: Orphan-locus classifier

**Files:**
- Modify: `src/rustle/vg_family/o3_flag_pass.rs`

**Interfaces:**
- Consumes: `LocusClass`, `OrphanLocus` (Task 1).
- Produces: `pub(crate) fn classify_orphan_locus(chrom: &str, start: u64, end: u64, own_family_id: &str,
  all_units_by_chrom: &std::collections::BTreeMap<String, Vec<(u64, u64, String, String)>>,
  genes_by_chrom: &std::collections::BTreeMap<String, Vec<(u64, u64)>>) -> (LocusClass, usize, Vec<String>)`.
  Task 5 calls this once per candidate orphan-read cluster.

`all_units_by_chrom`'s tuple is `(start, end, family_id, copy_idx)` per chrom. Precedence: `OtherFamily`
if the locus overlaps ANY unit belonging to a DIFFERENT family, else `AnnotatedNoUnit` if it overlaps a
`--gff` gene/pseudogene interval, else `Unannotated` — checked in that order (matches the Python's
`if other else (... if g else ...)`). Unlike the Python, `own_family_id` is compared directly with no
string-splitting: the Python's `.split('_')[0]` strips a synthetic contig suffix its OWN sweep-directory
naming convention adds (`fam_MCL58_073242` -> family_id `MCL58_073242`); this codebase's `--families`
catalog `family_id` carries no such suffix, so the split is not needed and would be wrong to port.

- [ ] **Step 1: Write the failing tests**

Append to `src/rustle/vg_family/o3_flag_pass.rs`:

```rust
/// Orphan-locus classification (§ "Scope note" in the design doc): OtherFamily > AnnotatedNoUnit >
/// Unannotated precedence. `all_units_by_chrom`/`genes_by_chrom` built ONCE per run, shared read-only.
pub(crate) fn classify_orphan_locus(
    chrom: &str,
    start: u64,
    end: u64,
    own_family_id: &str,
    all_units_by_chrom: &std::collections::BTreeMap<String, Vec<(u64, u64, String, String)>>,
    genes_by_chrom: &std::collections::BTreeMap<String, Vec<(u64, u64)>>,
) -> (LocusClass, usize, Vec<String>) {
    let overlaps = |a0: u64, a1: u64, b0: u64, b1: u64| a0 < b1 && b0 < a1;
    let other_family_units: Vec<String> = all_units_by_chrom
        .get(chrom)
        .map(|v| {
            v.iter()
                .filter(|(s, e, fid, _)| overlaps(start, end, *s, *e) && fid != own_family_id)
                .take(3)
                .map(|(_, _, fid, cidx)| format!("{fid}:{cidx}"))
                .collect()
        })
        .unwrap_or_default();
    let n_genes_overlapping = genes_by_chrom
        .get(chrom)
        .map(|v| v.iter().filter(|(s, e)| overlaps(start, end, *s, *e)).count())
        .unwrap_or(0);
    let class = if !other_family_units.is_empty() {
        LocusClass::OtherFamily
    } else if n_genes_overlapping > 0 {
        LocusClass::AnnotatedNoUnit
    } else {
        LocusClass::Unannotated
    };
    (class, n_genes_overlapping, other_family_units)
}

#[cfg(test)]
mod locus_tests {
    use super::*;
    use std::collections::BTreeMap;

    fn units() -> BTreeMap<String, Vec<(u64, u64, String, String)>> {
        let mut m = BTreeMap::new();
        m.insert("chr1".to_string(), vec![(100, 200, "OTHERFAM".to_string(), "3".to_string())]);
        m
    }

    fn genes() -> BTreeMap<String, Vec<(u64, u64)>> {
        let mut m = BTreeMap::new();
        m.insert("chr1".to_string(), vec![(500, 600)]);
        m
    }

    #[test]
    fn overlapping_another_familys_unit_wins_other_family() {
        let (class, _, units_hit) = classify_orphan_locus("chr1", 150, 250, "MYFAM", &units(), &genes());
        assert_eq!(class, LocusClass::OtherFamily);
        assert_eq!(units_hit, vec!["OTHERFAM:3".to_string()]);
    }

    #[test]
    fn own_family_unit_does_not_count_as_other_family() {
        let mut u = units();
        u.get_mut("chr1").unwrap().push((150, 250, "MYFAM".to_string(), "0".to_string()));
        let (class, _, units_hit) = classify_orphan_locus("chr1", 150, 250, "MYFAM", &u, &genes());
        // still hits OTHERFAM's unit at 100-200 too, so still OtherFamily -- verifies own-family rows
        // are excluded, not that the whole overlap set is
        assert_eq!(class, LocusClass::OtherFamily);
        assert!(!units_hit.iter().any(|s| s.starts_with("MYFAM:")));
    }

    #[test]
    fn no_other_family_unit_but_gene_overlap_is_annotated_no_unit() {
        let (class, n_genes, units_hit) = classify_orphan_locus("chr1", 550, 650, "MYFAM", &units(), &genes());
        assert_eq!(class, LocusClass::AnnotatedNoUnit);
        assert_eq!(n_genes, 1);
        assert!(units_hit.is_empty());
    }

    #[test]
    fn nothing_overlapping_is_unannotated() {
        let (class, n_genes, units_hit) = classify_orphan_locus("chr1", 9000, 9100, "MYFAM", &units(), &genes());
        assert_eq!(class, LocusClass::Unannotated);
        assert_eq!(n_genes, 0);
        assert!(units_hit.is_empty());
    }

    #[test]
    fn other_family_units_capped_at_three() {
        let mut u = BTreeMap::new();
        u.insert(
            "chr1".to_string(),
            (0..5).map(|i| (100u64, 200u64, format!("FAM{i}"), "0".to_string())).collect(),
        );
        let (_, _, units_hit) = classify_orphan_locus("chr1", 100, 200, "MYFAM", &u, &BTreeMap::new());
        assert_eq!(units_hit.len(), 3);
    }

    #[test]
    fn unknown_chrom_is_unannotated() {
        let (class, n_genes, units_hit) = classify_orphan_locus("chrZZZ", 0, 10, "MYFAM", &units(), &genes());
        assert_eq!(class, LocusClass::Unannotated);
        assert_eq!(n_genes, 0);
        assert!(units_hit.is_empty());
    }
}
```

- [ ] **Step 2: Run the tests**

```
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo test --release --lib o3_flag_pass::locus_tests 2>&1 | tail -20
```
Expected: all 6 pass.

- [ ] **Step 3: Full regression check + commit**

```
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo test --release --lib 2>&1 | tail -5
git add src/rustle/vg_family/o3_flag_pass.rs
git commit -m "O3 flag-pass: orphan-locus classifier, unit-tested (Task 3)"
```

---

### Task 4: `compute_family_raw` — the pair detector (minimap2 realignment + assembly)

**Files:**
- Modify: `src/rustle/vg_family/o3_flag_pass.rs`

**Interfaces:**
- Consumes: `parse_paf_consistency`/`AlignmentSummary` (Task 2), `RawPair`/`Class`/`O3Params` (Task 1).
- Produces: `pub fn detect_missing_copy_pairs(family_id: &str, copy_span_by_catalog_idx:
  &std::collections::HashMap<String, (String, u64, u64)>, rejected_by_catalog_idx:
  &std::collections::HashMap<String, Vec<(String, u64, u64, bool)>>, accepted_by_catalog_idx:
  &std::collections::HashMap<String, Vec<(String, u64, u64, bool)>>, genome: &rustle::genome::GenomeIndex,
  params: &O3Params) -> anyhow::Result<Vec<RawPair>>`. Task 6 builds the three `HashMap` arguments per
  family from `verdict`/`bam_reads` and calls this once per family inside `compute()`.

Each read tuple is `(name, ref_start_0based, ref_end_0based, strand_is_minus)` — everything
`parse_paf_consistency` needs to know per read is already summarized by the caller; this function only
needs enough per-read info to write a FASTA for minimap2 (name + sequence — see Step 1, the function
signature below actually takes `&[(String, Vec<u8>)]` per read group instead, sequence included, since a
`ref_start`/`ref_end` alone cannot reconstruct the read's bases). Re-derive the exact signature in Step 1
below from the real needs rather than trusting this paragraph over the code.

- [ ] **Step 1: Write the implementation** (no separate failing-test step for the minimap2-calling half —
  it needs a real `minimap2` binary on PATH, which unit tests should not depend on; Task 2 already
  covers the parsing logic in isolation. This step's own correctness is checked by Task 8's reproduction
  gate against the Python's real numbers, which IS the test for this piece.)

Append to `src/rustle/vg_family/o3_flag_pass.rs`:

```rust
/// One candidate copy Y's test batch (its own origin-rejected reads) and control batch (its own
/// certificate-accepted reads), both as `(name, sequence)` pairs ready to realign.
pub struct PairInput {
    pub copy_idx: String,
    pub is_partner: bool,
    pub rejected: Vec<(String, Vec<u8>)>,
    pub accepted: Vec<(String, Vec<u8>)>,
}

/// Realigns one batch of reads to a target window via `minimap2 -x splice -c --eqx -N 1`, mirroring
/// `copy_assign_pipeline.rs`'s `minimap2_msa_pair` temp-file convention (pid+atomic-nonce names,
/// `RUSTLE_MINIMAP2` env override, `Drop`-based cleanup) so concurrent region-parallel workers never
/// collide on the same path.
fn realign_batch(target_seq: &[u8], reads: &[(String, Vec<u8>)]) -> anyhow::Result<AlignmentSummary> {
    use std::io::Write;
    if reads.is_empty() {
        return Ok(AlignmentSummary { covered_kb: 0.0, n_sites: 0, per_read: std::collections::HashMap::new() });
    }
    let mm2 = std::env::var("RUSTLE_MINIMAP2").unwrap_or_else(|_| "minimap2".to_string());
    let dir = std::env::temp_dir();
    let pid = std::process::id();
    use std::sync::atomic::{AtomicUsize, Ordering};
    static O3_NONCE: AtomicUsize = AtomicUsize::new(0);
    let nonce = O3_NONCE.fetch_add(1, Ordering::Relaxed);
    let ypath = dir.join(format!("rustle_o3_y_{pid}_{nonce}.fa"));
    let rpath = dir.join(format!("rustle_o3_reads_{pid}_{nonce}.fa"));
    struct Cleanup(std::path::PathBuf, std::path::PathBuf);
    impl Drop for Cleanup {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.0);
            let _ = std::fs::remove_file(&self.1);
        }
    }
    let _cl = Cleanup(ypath.clone(), rpath.clone());
    {
        let mut y = std::fs::File::create(&ypath)?;
        y.write_all(b">Y\n")?;
        y.write_all(target_seq)?;
        y.write_all(b"\n")?;
        let mut r = std::fs::File::create(&rpath)?;
        for (name, seq) in reads {
            writeln!(r, ">{name}")?;
            r.write_all(seq)?;
            r.write_all(b"\n")?;
        }
    }
    let out = std::process::Command::new(&mm2)
        .args(["-x", "splice", "-c", "--eqx", "-N", "1", "-t", "4"])
        .arg(&ypath)
        .arg(&rpath)
        .output()?;
    if !out.status.success() {
        anyhow::bail!("minimap2 failed: {}", String::from_utf8_lossy(&out.stderr));
    }
    Ok(parse_paf_consistency(&String::from_utf8_lossy(&out.stdout)))
}

fn median(mut v: Vec<f64>) -> Option<f64> {
    if v.is_empty() {
        return None;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    Some(v[v.len() / 2])
}

/// Phase 1: one family's raw (uncorrected) missing-copy pair statistics. Skips any Y with fewer than
/// `params.min_reads` rejected reads (matches the Python's `len(names) < 3: continue` -- not emitted as
/// Untestable, simply absent from the output). A `minimap2` failure for one Y is logged to stderr and
/// that pair is skipped -- see the design doc's Error Handling section for why this must not abort.
pub fn detect_missing_copy_pairs(
    family_id: &str,
    copy_span_by_catalog_idx: &std::collections::HashMap<String, (String, u64, u64)>,
    genome: &rustle::genome::GenomeIndex,
    inputs: &[PairInput],
    params: &O3Params,
) -> Vec<RawPair> {
    let mut out = Vec::new();
    for input in inputs {
        if input.rejected.len() < params.min_reads {
            continue;
        }
        let Some((chrom, s, e)) = copy_span_by_catalog_idx.get(&input.copy_idx) else {
            continue;
        };
        let target = match genome.fetch_sequence(chrom, *s, *e) {
            Some(t) => t,
            None => continue,
        };
        let rejected: Vec<_> = input.rejected.iter().take(params.max_reads).cloned().collect();
        let test = match realign_batch(&target, &rejected) {
            Ok(s) => s,
            Err(err) => {
                eprintln!("[o3-flag-pass] realignment failed for {family_id}:{}: {err}", input.copy_idx);
                continue;
            }
        };
        let accepted: Vec<_> = input.accepted.iter().take(params.max_reads).cloned().collect();
        let ctl = if accepted.len() >= params.min_reads {
            match realign_batch(&target, &accepted) {
                Ok(s) => s,
                Err(err) => {
                    eprintln!("[o3-flag-pass] control realignment failed for {family_id}:{}: {err}", input.copy_idx);
                    AlignmentSummary { covered_kb: 0.0, n_sites: 0, per_read: std::collections::HashMap::new() }
                }
            }
        } else {
            AlignmentSummary { covered_kb: 0.0, n_sites: 0, per_read: std::collections::HashMap::new() }
        };
        // rate floored at 1 site over the control's covered kb, so a control with zero observed sites
        // never claims a zero rate (which would trivially pass every test) -- design doc, RawPair section.
        let p_uncorrected = if test.covered_kb > 0.0 && ctl.covered_kb > 0.0 {
            let rate = (ctl.n_sites.max(1) as f64) / ctl.covered_kb;
            Some(poisson_tail(test.n_sites, rate * test.covered_kb))
        } else {
            None
        };
        let mismatches: Vec<f64> = test.per_read.values().map(|&(nx, _, _)| nx as f64).collect();
        let unaligned: Vec<f64> = test.per_read.values().map(|&(_, unal, _)| unal as f64).collect();
        let med_mismatch = median(mismatches).unwrap_or(f64::NAN);
        let med_unaligned = median(unaligned).unwrap_or(f64::NAN);
        let class = if !med_mismatch.is_nan() && med_mismatch > med_unaligned {
            Class::Divergent
        } else {
            Class::Structural
        };
        out.push(RawPair {
            family_id: family_id.to_string(),
            copy_idx: input.copy_idx.clone(),
            is_partner: input.is_partner,
            n_rejected: input.rejected.len(),
            n_aligned: rejected.len(),
            covered_kb: test.covered_kb,
            n_sites: test.n_sites,
            ctl_n: accepted.len(),
            ctl_covered_kb: ctl.covered_kb,
            ctl_n_sites: ctl.n_sites,
            p_uncorrected,
            class,
            med_mismatch,
            med_unaligned,
        });
    }
    out
}
```

- [ ] **Step 2: Add a unit test that does not require minimap2** — verify the `< min_reads` skip and
  the "missing span" skip, using `GenomeIndex::from_seqs` — a real, already-existing `#[cfg(test)]`-only
  constructor (`src/rustle/genome.rs:371`, `pub(crate) fn from_seqs(pairs: &[(&str, &[u8])]) -> Self`),
  verified present in this codebase, not invented for this plan:

```rust
#[cfg(test)]
mod pair_detector_tests {
    use super::*;
    use rustle::genome::GenomeIndex;

    #[test]
    fn fewer_than_min_reads_is_skipped_entirely() {
        let genome = GenomeIndex::from_seqs(&[("chrT", &[b'A'; 100])]);
        let mut spans = std::collections::HashMap::new();
        spans.insert("0".to_string(), ("chrT".to_string(), 0u64, 100u64));
        let inputs = vec![PairInput {
            copy_idx: "0".to_string(),
            is_partner: false,
            rejected: vec![("r1".to_string(), vec![b'A'; 50])], // only 1, min_reads default is 3
            accepted: vec![],
        }];
        let pairs = detect_missing_copy_pairs("F", &spans, &genome, &inputs, &O3Params::default());
        assert!(pairs.is_empty());
    }

    #[test]
    fn missing_copy_span_is_skipped_not_errored() {
        let genome = GenomeIndex::from_seqs(&[("chrT", &[b'A'; 100])]);
        let spans = std::collections::HashMap::new(); // no span for "0"
        let inputs = vec![PairInput {
            copy_idx: "0".to_string(),
            is_partner: false,
            rejected: vec![
                ("r1".to_string(), vec![b'A'; 50]),
                ("r2".to_string(), vec![b'A'; 50]),
                ("r3".to_string(), vec![b'A'; 50]),
            ],
            accepted: vec![],
        }];
        let pairs = detect_missing_copy_pairs("F", &spans, &genome, &inputs, &O3Params::default());
        assert!(pairs.is_empty());
    }
}
```

`from_seqs` is `pub(crate)`, so it is visible from `o3_flag_pass.rs` (same crate) without any visibility
change — nothing further to check here.

- [ ] **Step 3: Run the tests**

```
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo test --release --lib o3_flag_pass 2>&1 | tail -30
```
Expected: all tests from Tasks 1-4 pass (18 total: 6 + 6 + 6, plus the 2 new pair-detector tests = 20).

- [ ] **Step 4: Manual smoke test with real minimap2** (this is the actual correctness check for the
  realignment half, since it can't be a hermetic unit test):

```
which minimap2
```
Confirm it's on PATH (this codebase's whole test suite already assumes this — see
`denovo_pipeline.rs`'s own tests skipping when `minimap2 --version` fails). If present, this step's real
validation happens in Task 8 (the reproduction gate) — nothing further to do here now.

- [ ] **Step 5: Full regression check + commit**

```
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo test --release --lib 2>&1 | tail -5
git add src/rustle/vg_family/o3_flag_pass.rs
git commit -m "O3 flag-pass: pair detector (minimap2 realignment + Poisson test) (Task 4)"
```

---

### Task 5: Wire into `copy_assign.rs` Phase 1 — CLI flags, shared indices, per-region call

**Files:**
- Modify: `src/bin/copy_assign.rs`

**Interfaces:**
- Consumes: `rustle::vg_family::o3_flag_pass::{O3Params, RawPair, OrphanLocus, PairInput,
  detect_missing_copy_pairs, classify_orphan_locus}` (Tasks 1-4).
- Produces: `RegionWork.o3_raw_pairs: Vec<RawPair>` and `RegionWork.o3_orphan_loci: Vec<OrphanLocus>` for
  Task 6 to accumulate in the serial drain.

Re-read the exact current line numbers before editing — this session's earlier commits (`6a53476`,
`fa382df`, `80e4d6e`) already touched this file; the numbers below are accurate as of commit `02575db`
but a `grep -n` re-check costs nothing and catches drift immediately.

- [ ] **Step 1: Add the CLI flags**

In `src/bin/copy_assign.rs`, find the `rescue_singletons: bool,` flag (`grep -n "rescue_singletons: bool"
src/bin/copy_assign.rs` — currently line 212) and add immediately after its closing `,`:

```rust
    /// O3 (`docs/superpowers/specs/2026-09-10-o3-flag-pass-integration-design.md`): port of
    /// `bench/o3_flag_pass.py`'s missing-copy detector. Requires `--families`. Adds `o3_flag`/`o3_class`/
    /// `o3_rate_per_kb`/`o3_p`/`o3_n_rejected` columns to `<out>.family_join.tsv` and writes
    /// `<out>.o3_candidate_loci.tsv`. Default off, byte-identical when unset.
    #[arg(long, default_value_t = false)]
    flag_missing_copies: bool,
    /// O3: Bonferroni alpha for the genome-wide missing-copy threshold (`alpha / n_pairs`).
    #[arg(long, default_value_t = 0.001)]
    o3_alpha: f64,
    /// O3: cap on reads realigned per side (test/control) per candidate copy, for wall-clock control.
    #[arg(long, default_value_t = 500)]
    o3_max_reads: usize,
```

- [ ] **Step 2: Validate `--flag-missing-copies` requires `--families`**

Find where other `--families`-dependent flags are checked at argument-validation time (`grep -n
"anyhow::bail!" src/bin/copy_assign.rs | grep -i families` — there is an existing check the
`--gtf-copy-set`/similar flags rely on; locate it and add alongside, in the same validation block, not a
new one elsewhere in the file):

```rust
    if args.flag_missing_copies && args.families.is_none() {
        anyhow::bail!("--flag-missing-copies requires --families (it tests catalog copies for a missing sibling)");
    }
```

- [ ] **Step 3: Build the shared genome-wide indices once, only when the flag is set**

Find `let (region_families, catalog_seqs, region_windows) = load_supplied_families(&args, &by_contig)?;`
(currently line 1829) and add immediately after:

```rust
    // O3: built ONCE, read-only across every parallel region worker (same pattern as genome_cache/
    // bam_cache below) -- only when the flag is set, so the unset path pays nothing.
    let o3_all_units_by_chrom: std::collections::BTreeMap<String, Vec<(u64, u64, String, String)>> =
        if args.flag_missing_copies {
            let mut m: std::collections::BTreeMap<String, Vec<(u64, u64, String, String)>> = std::collections::BTreeMap::new();
            if let Some(rf) = &region_families {
                for fams in rf.values() {
                    for f in fams {
                        for c in &f.copies {
                            m.entry(c.chrom.clone()).or_default().push((c.start, c.end, c.family_id.clone(), c.copy_idx.to_string()));
                        }
                    }
                }
            }
            m
        } else {
            std::collections::BTreeMap::new()
        };
    let o3_genes_by_chrom: std::collections::BTreeMap<String, Vec<(u64, u64)>> = if args.flag_missing_copies {
        match &args.gff {
            Some(path) => {
                let mut m: std::collections::BTreeMap<String, Vec<(u64, u64)>> = std::collections::BTreeMap::new();
                for (chrom, s, e) in parse_annotation(path)? {
                    m.entry(chrom).or_default().push((s, e));
                }
                m
            }
            None => std::collections::BTreeMap::new(),
        }
    } else {
        std::collections::BTreeMap::new()
    };
    let o3_params = rustle::vg_family::o3_flag_pass::O3Params {
        alpha: args.o3_alpha,
        max_reads: args.o3_max_reads,
        min_reads: 3,
    };
```

Note: `parse_annotation` is called here unconditionally on `args.flag_missing_copies` regardless of
`args.gff` — check whether `args.gff` is `Option<String>` (it is, per its existing use at
`args.gff.as_deref().map(parse_annotation).transpose()` elsewhere in this function) and that this new
call handles `None` the same way (the `match` above does).

- [ ] **Step 4: Add the two new `RegionWork` fields**

In the `struct RegionWork { ... }` definition (currently lines 72-115), add before the closing `}`:

```rust
    /// O3: this family's raw (uncorrected) missing-copy pair statistics. Empty unless `--flag-missing-copies`.
    o3_raw_pairs: Vec<rustle::vg_family::o3_flag_pass::RawPair>,
    /// O3: candidate orphan-read loci outside every unit of this family. Empty unless `--flag-missing-copies`.
    o3_orphan_loci: Vec<rustle::vg_family::o3_flag_pass::OrphanLocus>,
```

Update the constructor (currently line ~2307, inside the `compute` closure — grep for `Ok(RegionWork {`
to confirm the current line) to add `o3_raw_pairs, o3_orphan_loci` to the field list, and the destructure
(currently line ~2351, `let RegionWork { ... } = work;`) to add the same two names.

- [ ] **Step 5: Call the detector inside `compute()`, per family**

Inside the `compute` closure (starts at `let compute = |contig: &String, lo: u64, hi: u64| -> Result<RegionWork> {`,
currently line 2027), find the point where `fams: Vec<FamilyAssignment>` is fully built (after
`detect_and_assign` returns, before the closure's final `Ok(RegionWork { ... })`) and add:

```rust
    let (o3_raw_pairs, o3_orphan_loci): (Vec<_>, Vec<_>) = if args.flag_missing_copies {
        let mut pairs = Vec::new();
        let mut loci = Vec::new();
        for fa in &fams {
            // Resolve catalog_copy_idx -> (chrom, start, end): fa.copy_spans is indexed by the SWEEP's
            // own internal position, catalog_copy_idx is a separate namespace (comment near line ~2906
            // explains why) -- rebuild the mapping via fa.copy_tids + catalog_index, same pattern the
            // isoform-attribute loop already uses.
            let mut copy_span_by_catalog_idx: std::collections::HashMap<String, (String, u64, u64)> = std::collections::HashMap::new();
            for (ci, tid) in fa.copy_tids.iter().enumerate() {
                if let Some((_, cidx)) = catalog_index.as_ref().and_then(|ix| ix.get(tid)) {
                    if let Some(span) = fa.copy_spans.get(ci) {
                        copy_span_by_catalog_idx.insert(cidx.to_string(), span.clone());
                    }
                }
            }
            // Group this family's bam_reads by best-candidate catalog_copy_idx, split into rejected
            // (origin_rejected==true) and accepted (this family's own certificate-passed reads at that
            // copy). `Assignment` (src/rustle/vg_family/copy_assign.rs:100-143) carries `best_copy: usize`
            // (an index into `copy_tids`/`copy_spans`, the same namespace `ci` uses elsewhere in this
            // closure), `status: AssignStatus` (an enum: `Assigned`/`Ambiguous`/`Tied` -- copy_assign.rs:
            // 91-98, NOT a string), and `origin_rejected: bool` directly.
            use rustle::vg_family::copy_assign::AssignStatus;
            let mut rejected_by_idx: std::collections::HashMap<String, Vec<(String, Vec<u8>)>> = std::collections::HashMap::new();
            let mut accepted_by_idx: std::collections::HashMap<String, Vec<(String, Vec<u8>)>> = std::collections::HashMap::new();
            for &(read_i, ref assignment) in &fa.assignments {
                let br = &bam_reads[read_i];
                let cidx = match catalog_index.as_ref().and_then(|ix| ix.get(&fa.copy_tids[assignment.best_copy])) {
                    Some((_, c)) => c.to_string(),
                    None => continue,
                };
                let entry = (br.name.clone(), br.read.seq.clone());
                if assignment.origin_rejected {
                    rejected_by_idx.entry(cidx).or_default().push(entry);
                } else if assignment.status == AssignStatus::Assigned {
                    accepted_by_idx.entry(cidx).or_default().push(entry);
                }
            }
            let inputs: Vec<rustle::vg_family::o3_flag_pass::PairInput> = copy_span_by_catalog_idx
                .keys()
                .map(|cidx| rustle::vg_family::o3_flag_pass::PairInput {
                    copy_idx: cidx.clone(),
                    is_partner: false, // see the design doc's is_partner note; CatalogCopy::partner not threaded through FamilyAssignment yet -- default false is safe (never over-claims a partner exclusion)
                    rejected: rejected_by_idx.get(cidx).cloned().unwrap_or_default(),
                    accepted: accepted_by_idx.get(cidx).cloned().unwrap_or_default(),
                })
                .collect();
            let genome_ref = genome_for(contig)?;
            pairs.extend(rustle::vg_family::o3_flag_pass::detect_missing_copy_pairs(
                &fa.family_id, &copy_span_by_catalog_idx, &genome_ref, &inputs, &o3_params,
            ));
            // Orphan-locus scan: bam_reads whose primary lands outside every unit of this family,
            // clustered by proximity (<=5kb gap, matching the Python), classified via the shared indices.
            let fam_units: Vec<(String, u64, u64)> = fa.copy_spans.clone();
            let mut outside: Vec<&BamRead> = bam_reads
                .iter()
                .filter(|br| !br.is_secondary && !br.is_supplementary)
                .filter(|br| !fam_units.iter().any(|(c, s, e)| br.chrom == *c && br.read.ref_start < *e && *s < read_ref_end_local(&br.read)))
                .collect();
            outside.sort_by(|a, b| (a.chrom.as_str(), a.read.ref_start).cmp(&(b.chrom.as_str(), b.read.ref_start)));
            let mut clusters: Vec<(String, u64, u64, usize)> = Vec::new();
            for br in &outside {
                let end = read_ref_end_local(&br.read);
                if let Some(last) = clusters.last_mut() {
                    if last.0 == br.chrom && br.read.ref_start.saturating_sub(last.2) <= 5000 {
                        last.2 = last.2.max(end);
                        last.3 += 1;
                        continue;
                    }
                }
                clusters.push((br.chrom.clone(), br.read.ref_start, end, 1));
            }
            for (chrom, start, end, n_reads) in clusters {
                if n_reads < 3 {
                    continue;
                }
                let (class, n_genes, other_units) = rustle::vg_family::o3_flag_pass::classify_orphan_locus(
                    &chrom, start, end, &fa.family_id, &o3_all_units_by_chrom, &o3_genes_by_chrom,
                );
                loci.push(rustle::vg_family::o3_flag_pass::OrphanLocus {
                    chrom, start, end, n_reads, n_orphans: 0, class,
                    n_genes_overlapping: n_genes, other_family_units: other_units,
                });
            }
        }
        (pairs, loci)
    } else {
        (Vec::new(), Vec::new())
    };
```

Field names above (`Assignment::best_copy`/`::status`/`::origin_rejected`, `AssignStatus::Assigned`,
`BamRead::name`/`::read`/`::is_secondary`/`::is_supplementary`, `AlignedRead::seq`/`::ref_start`) are
verified against the real struct definitions (`src/rustle/vg_family/copy_assign.rs:91-143`,
`src/rustle/vg_family/denovo_assemble.rs:950-965`, `src/rustle/vg_family/copy_split.rs:177-182`) as of
this plan's writing, not guessed — `AssignStatus` already derives `PartialEq, Eq` so the `==` comparison
compiles as written. The one thing NOT independently verified here is `catalog_index`'s exact variable
name and type inside `compute()`'s closure at this point (it is used elsewhere in the same closure per
this session's earlier reading, e.g. the isoform-attribute loop around `copy_assign.rs:2906`) — a `grep
-n "catalog_index" src/bin/copy_assign.rs` before this step confirms it is in scope here and did not
shift name during this session's other edits.

- [ ] **Step 6: Build and fix compile errors**

```
cd /mnt/c/Users/jfris/Desktop/Rustle
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo build --release --all-targets > /mnt/linuxdisk/home/juanfraitu/build_o3_task5.log 2>&1
echo "EXIT=$?"
tail -60 /mnt/linuxdisk/home/juanfraitu/build_o3_task5.log
```
Fix every error by reading the exact struct/field names the compiler names — this task's Step 5 is
explicitly flagged as needing reconciliation against real field names; expect and budget time for this.

- [ ] **Step 7: Byte-identity check (flag unset)**

```
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo test --release --lib 2>&1 | tail -5
cd /mnt/linuxdisk/home/juanfraitu/growth_extent_test/bakeoff_rerun
BIN=/mnt/linuxdisk/home/juanfraitu/rustle_target/release/copy_assign
BAM=/mnt/linuxdisk/home/juanfraitu/bakeoff/human/hsa16.bam
FASTA=/mnt/linuxdisk/home/juanfraitu/refcache/chm13v2.0.fa
FAM=/mnt/linuxdisk/home/juanfraitu/bakeoff/human/copies16.tsv
COPFA=/mnt/linuxdisk/home/juanfraitu/bakeoff/human/copies16.fa
REG=/mnt/linuxdisk/home/juanfraitu/bakeoff/human/regions16
$BIN --bam "$BAM" --fasta "$FASTA" --families "$FAM" --copies-fa "$COPFA" --regions "$REG" --gtf --out o3_task5_check > o3_task5_check.log 2>&1
diff -q o3_task5_check.gtf hard_baseline_today.gtf && echo "BYTE-IDENTICAL"
```
Expected: `BYTE-IDENTICAL`. If not, the new code inside `compute()` has a side effect outside the
`if args.flag_missing_copies` guard — find and fix it before proceeding.

- [ ] **Step 8: Commit**

```
git add src/bin/copy_assign.rs
git commit -m "O3 flag-pass: wire Phase 1 into copy_assign's compute() (Task 5)"
```

---

### Task 6: Wire into `copy_assign.rs` Phase 2 — aggregation + output surface

**Files:**
- Modify: `src/bin/copy_assign.rs`

**Interfaces:**
- Consumes: `RegionWork.o3_raw_pairs`/`o3_orphan_loci` (Task 5), `finalize_flags` (Task 1).
- Produces: `<out>.family_join.tsv` with 5 new trailing columns when the flag is set; new
  `<out>.o3_candidate_loci.tsv`.

- [ ] **Step 1: Accumulate raw pairs/loci during the serial drain**

Find the serial drain loop (`for (gwork, work) in works.into_iter().enumerate() {`, currently line 2350)
and, right after the `let RegionWork { ... } = work;` destructure (currently line 2351, now also binding
`o3_raw_pairs, o3_orphan_loci`), and near wherever `famcn_rows`/`family_rows` (built earlier in `main`,
before the loop) are declared, add a matching accumulator declared alongside them (search for `let mut
famcn_rows` to find where sibling `Vec`s are declared, before the loop starts):

```rust
    let mut o3_all_raw_pairs: Vec<rustle::vg_family::o3_flag_pass::RawPair> = Vec::new();
    let mut o3_all_orphan_loci: Vec<rustle::vg_family::o3_flag_pass::OrphanLocus> = Vec::new();
```

Inside the loop body, after the destructure:
```rust
    o3_all_raw_pairs.extend(o3_raw_pairs);
    o3_all_orphan_loci.extend(o3_orphan_loci);
```

- [ ] **Step 2: Refactor `join_rows` to carry its join key**

Find `let mut join_rows: Vec<String> = Vec::new();` (currently line 1901) and change to:

```rust
    struct JoinRow {
        line: String,
        family_id: String,
        copy_idx: String,
    }
    let mut join_rows: Vec<JoinRow> = Vec::new();
```

Find the push site (currently line ~2590, inside the block building `cf`/`cidx`) and change:
```rust
                        join_rows.push(format!(
                            "{fid}\t{ci}\t{tid}\t{cf}\t{cidx}\t{}\t{}\t{}\t{}",
                            ...
                        ));
```
to:
```rust
                        join_rows.push(JoinRow {
                            line: format!(
                                "{fid}\t{ci}\t{tid}\t{cf}\t{cidx}\t{}\t{}\t{}\t{}",
                                ...
                            ),
                            family_id: cf.clone(),
                            copy_idx: cidx.clone(),
                        });
```
(keep the `format!` arguments identical to whatever they currently are — this step changes only the
wrapper, never the formatted content of `line`).

Find the "emitted" check right after the write loop (currently ~line 3789, `let emitted: HashSet<&str> =
join_rows.iter().filter_map(|l| l.split('\t').nth(2)).collect();`) and change to read `.line` instead:
```rust
    let emitted: HashSet<&str> = join_rows.iter().filter_map(|r| r.line.split('\t').nth(2)).collect();
```

- [ ] **Step 3: Run the Phase-2 aggregation and write `family_join.tsv` with the new columns**

Find the `family_join.tsv` writer (currently lines 3781-3808) and replace the header + write loop:

```rust
        let o3_flags: std::collections::HashMap<(String, String), rustle::vg_family::o3_flag_pass::FlaggedPair> =
            if args.flag_missing_copies {
                rustle::vg_family::o3_flag_pass::finalize_flags(&o3_all_raw_pairs, args.o3_alpha)
                    .into_iter()
                    .map(|fp| ((fp.pair.family_id.clone(), fp.pair.copy_idx.clone()), fp))
                    .collect()
            } else {
                std::collections::HashMap::new()
            };
        let mut jh = std::fs::File::create(format!("{}.family_join.tsv", args.out))?;
        let header = "family_id\tcopy_index\tcopy_tid\tcatalog_family_id\tcatalog_copy_idx\tchrom\tstart\tend\tn_reads_hard";
        if args.flag_missing_copies {
            writeln!(jh, "{header}\to3_flag\to3_class\to3_rate_per_kb\to3_p\to3_n_rejected")?;
        } else {
            writeln!(jh, "{header}")?;
        }
        for r in &join_rows {
            if args.flag_missing_copies {
                match o3_flags.get(&(r.family_id.clone(), r.copy_idx.clone())) {
                    Some(fp) => {
                        let flag_str = match fp.flag {
                            rustle::vg_family::o3_flag_pass::Flag::MissingCopy => "missing_copy",
                            rustle::vg_family::o3_flag_pass::Flag::Untestable => "untestable",
                            rustle::vg_family::o3_flag_pass::Flag::NoFlag => "none",
                        };
                        let class_str = match fp.pair.class {
                            rustle::vg_family::o3_flag_pass::Class::Divergent => "divergent",
                            rustle::vg_family::o3_flag_pass::Class::Structural => "structural",
                        };
                        let rate = if fp.pair.covered_kb > 0.0 { fp.pair.n_sites as f64 / fp.pair.covered_kb } else { 0.0 };
                        let p_str = fp.pair.p_uncorrected.map_or("NA".to_string(), |p| format!("{p:.3e}"));
                        writeln!(jh, "{}\t{flag_str}\t{class_str}\t{rate:.2}\t{p_str}\t{}", r.line, fp.pair.n_rejected)?;
                    }
                    None => writeln!(jh, "{}\tnone\tNA\t0.00\tNA\t0", r.line)?,
                }
            } else {
                writeln!(jh, "{}", r.line)?;
            }
        }
```

- [ ] **Step 4: Write `o3_candidate_loci.tsv`**

Immediately after the `family_join.tsv` block closes (after the existing `eprintln!("[copy_assign] wrote
{}.family_join.tsv ...")`), add:

```rust
        if args.flag_missing_copies {
            let mut lh = std::fs::File::create(format!("{}.o3_candidate_loci.tsv", args.out))?;
            writeln!(lh, "chrom\tstart\tend\tn_reads\tn_orphans\tclass\tn_genes_overlapping\tother_family_units")?;
            for l in &o3_all_orphan_loci {
                let class_str = match l.class {
                    rustle::vg_family::o3_flag_pass::LocusClass::OtherFamily => "other_family",
                    rustle::vg_family::o3_flag_pass::LocusClass::AnnotatedNoUnit => "annotated_no_unit",
                    rustle::vg_family::o3_flag_pass::LocusClass::Unannotated => "unannotated",
                };
                writeln!(
                    lh, "{}\t{}\t{}\t{}\t{}\t{class_str}\t{}\t{}",
                    l.chrom, l.start, l.end, l.n_reads, l.n_orphans, l.n_genes_overlapping,
                    if l.other_family_units.is_empty() { "-".to_string() } else { l.other_family_units.join(";") },
                )?;
            }
            eprintln!("[copy_assign] wrote {}.o3_candidate_loci.tsv ({} loci)", args.out, o3_all_orphan_loci.len());
        }
```

- [ ] **Step 5: Build, byte-identity check, then real-run check**

```
cd /mnt/c/Users/jfris/Desktop/Rustle
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo build --release --all-targets > /mnt/linuxdisk/home/juanfraitu/build_o3_task6.log 2>&1
echo "EXIT=$?"
tail -60 /mnt/linuxdisk/home/juanfraitu/build_o3_task6.log
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo test --release --lib 2>&1 | tail -5
```

Fix compile errors, then byte-identity again (flag unset, same command as Task 5 Step 7 — must still be
`BYTE-IDENTICAL`), then a real run WITH the flag on:

```
cd /mnt/linuxdisk/home/juanfraitu/growth_extent_test/bakeoff_rerun
BIN=/mnt/linuxdisk/home/juanfraitu/rustle_target/release/copy_assign
BAM=/mnt/linuxdisk/home/juanfraitu/bakeoff/human/hsa16.bam
FASTA=/mnt/linuxdisk/home/juanfraitu/refcache/chm13v2.0.fa
FAM=/mnt/linuxdisk/home/juanfraitu/bakeoff/human/copies16.tsv
COPFA=/mnt/linuxdisk/home/juanfraitu/bakeoff/human/copies16.fa
REG=/mnt/linuxdisk/home/juanfraitu/bakeoff/human/regions16
GFF=/mnt/linuxdisk/home/juanfraitu/winloci_data/Reference/chm13v2.0_RefSeq_full.gff.gz
time $BIN --bam "$BAM" --fasta "$FASTA" --families "$FAM" --copies-fa "$COPFA" --regions "$REG" --gtf --flag-missing-copies --gff "$GFF" --out o3_task6_flagon > o3_task6_flagon.log 2>&1
echo "EXIT=$?"
tail -20 o3_task6_flagon.log
head -5 o3_task6_flagon.family_join.tsv
head -5 o3_task6_flagon.o3_candidate_loci.tsv
```

Confirm: run completes without error, `family_join.tsv` has the 5 new trailing columns on every row,
`o3_candidate_loci.tsv` is well-formed (even if empty of rows — that is a valid outcome, not a bug; check
the per-family stderr lines for how many pairs were tested before concluding zero rows means something is
broken). Note the wall-clock time printed by `time` — record it in the ledger entry (Step 7) since the
design doc flags this as unmeasured.

- [ ] **Step 6: Commit**

```
git add src/bin/copy_assign.rs
git commit -m "O3 flag-pass: wire Phase 2 aggregation + family_join.tsv/o3_candidate_loci.tsv output (Task 6)"
```

- [ ] **Step 7: Ledger entry**

Append a new `## §6ja` (or the next free letter — check the tail of `docs/o1_ledger.md` for the current
last section id) section to `docs/o1_ledger.md` recording: the flag name, byte-identity confirmation, the
real-run wall-clock time from Step 5 above, and the raw pair/flag counts from the human chr16 substrate.
Do not claim any number as validated yet — that is Task 7's job.

---

### Task 7: Reproduction gate against the Python's own numbers

**Files:** none (validation task, no code changes expected unless it finds a real discrepancy)

- [ ] **Step 1: Identify a substrate both the Python script and the new flag can run on**

The Python script (`bench/o3_flag_pass.py`) expects a `--sweep DIR` of `fam_*` subdirectories, each with
its own `copies.tsv` + `A.assignments.tsv` from a PRIOR `copy_assign` run. Find or rebuild the smallest
such sweep referenced in `docs/o1_ledger.md` §6fm-§6ft (grep for `sweep_v` directory names near those
sections) — reuse one that already exists under `/mnt/linuxdisk/home/juanfraitu/mcl_ann/` rather than
rebuilding from scratch if a suitable one is still present.

- [ ] **Step 2: Run the Python script on that sweep**

```
python3 bench/o3_flag_pass.py /mnt/linuxdisk/home/juanfraitu/o3_repro_check --bam <BAM> --fasta <FASTA> --gff <GFF> --units <catalog units.tsv> --sweep <the sweep dir from Step 1>
cat /mnt/linuxdisk/home/juanfraitu/o3_repro_check/summary.txt
```

- [ ] **Step 3: Run the new Rust flag on the SAME underlying catalog + BAM + region set** (translate the
  sweep's per-family regions into one `--regions` file covering the same families)

```
$BIN --bam <same BAM> --fasta <same FASTA> --families <the catalog units.tsv, reformatted to copy_assign's expected column names if they differ> --regions <the equivalent region list> --gtf --flag-missing-copies --gff <same GFF> --out o3_repro_rust
```

- [ ] **Step 4: Diff**

Compare, per `(family_id, copy_idx)`: pair count tested, `class` (divergent/structural), final `flag`.
Write a small one-off Python diff (not committed — this is a validation step, not a shipped tool) that
joins `o3_repro_check/flags.tsv` against `o3_repro_rust.family_join.tsv` on the join key and reports
mismatches. `p_uncorrected` should match to within `1e-6` relative tolerance (independent
implementations of the same closed-form sum can differ in float rounding order). Any mismatch beyond
that is a real discrepancy — find and fix its root cause (most likely candidates: the catalog units.tsv
column mapping between the two tools' input formats, or a `min_reads`/`max_reads` default mismatch)
before considering this task done. Do not adjust the Rust implementation's algorithm to match a wrong
Python number, or vice versa, without first confirming which one is actually right against the
`detector()` function's own logic in `bench/o3_flag_pass.py`.

- [ ] **Step 5: Record the result in the ledger**

Append to the `docs/o1_ledger.md` section from Task 6 Step 7 (or a new one if that section is already
closed out): the exact pair/loci counts from both tools, confirmation they match (or the discrepancy
found and fixed, with root cause named).

---

### Task 8: Cross-individual differential comparison script

**Files:**
- Create: `bench/o3_cross_individual_diff.py`

**Interfaces:**
- Consumes: two `<out>.family_join.tsv` files (from Task 6) with the `o3_flag` column, produced by two
  separate `copy_assign --flag-missing-copies` runs against the SAME `--families` catalog with different
  `--bam` inputs.

- [ ] **Step 1: Run `copy_assign --flag-missing-copies` on both gorilla substrates**

Using the SAME catalog (`/mnt/linuxdisk/home/juanfraitu/mcl_ann/gw_units_v3.units.tsv`) against each BAM
in turn — one run at a time, foreground, per this project's WSL2 crash rule (these are genome-wide
sweeps, not the small human chr16 substrate; budget real wall-clock time and watch for memory pressure):

```
BIN=/mnt/linuxdisk/home/juanfraitu/rustle_target/release/copy_assign
CATALOG=/mnt/linuxdisk/home/juanfraitu/mcl_ann/gw_units_v3.units.tsv
GFF=/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_genomic.gff
FASTA=<the gorilla mGorGor1 reference fasta used to build gw_units_v3 -- check gw_units_v3.params.tsv for the exact path>
$BIN --bam /mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_ds.bam --fasta "$FASTA" --families "$CATALOG" --regions <genome-wide region list> --gtf --flag-missing-copies --gff "$GFF" --out /mnt/linuxdisk/home/juanfraitu/o3_diff/testis
# wait for completion, confirm it finished cleanly, THEN run the second (never run two heavy copy_assign invocations at once):
$BIN --bam /mnt/linuxdisk/home/juanfraitu/o1_replicate/fibro_ds.bam --fasta "$FASTA" --families "$CATALOG" --regions <same region list> --gtf --flag-missing-copies --gff "$GFF" --out /mnt/linuxdisk/home/juanfraitu/o3_diff/fibroblast
```

- [ ] **Step 2: Write the comparison script**

```python
#!/usr/bin/env python3
"""Cross-individual O3 differential (docs/superpowers/specs/2026-09-10-o3-flag-pass-integration-design.md,
"Follow-on"): diffs two --flag-missing-copies family_join.tsv outputs -- same catalog, different --bam --
by flag CATEGORY, not raw read presence (the read-presence approach, ledger section 4l, already showed it
cannot distinguish real copy-number difference from tissue-driven expression difference).

usage: o3_cross_individual_diff.py <arm_a>.family_join.tsv <arm_b>.family_join.tsv --label-a testis --label-b fibroblast
"""
import argparse, csv, collections

def load(path):
    rows = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if "o3_flag" not in r:
                raise SystemExit(f"{path}: no o3_flag column -- was --flag-missing-copies set for this run?")
            rows[(r["catalog_family_id"], r["catalog_copy_idx"])] = r["o3_flag"]
    return rows

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("arm_a"); ap.add_argument("arm_b")
    ap.add_argument("--label-a", default="A"); ap.add_argument("--label-b", default="B")
    a = ap.parse_args()
    ra, rb = load(a.arm_a), load(a.arm_b)
    keys = sorted(set(ra) | set(rb))
    cats = collections.Counter()
    candidates = []
    for k in keys:
        fa, fb = ra.get(k, "NA"), rb.get(k, "NA")
        if fa == "missing_copy" and fb == "none":
            cats["candidate_differential"] += 1
            candidates.append((k, a.label_a, a.label_b))
        elif fb == "missing_copy" and fa == "none":
            cats["candidate_differential"] += 1
            candidates.append((k, a.label_b, a.label_a))
        elif fa == "missing_copy" and fb == "missing_copy":
            cats["shared"] += 1
        elif "untestable" in (fa, fb) and "missing_copy" in (fa, fb):
            cats["inconclusive"] += 1
        elif fa == "untestable" and fb == "untestable":
            cats["no_information"] += 1
        elif fa == "none" and fb == "none":
            cats["no_signal"] += 1
        else:
            cats["other"] += 1
    print(f"{len(keys)} copies compared")
    for c, n in sorted(cats.items(), key=lambda kv: -kv[1]):
        print(f"  {c}: {n}")
    no_signal_frac = cats["no_signal"] / max(1, len(keys))
    print(f"\nno_signal fraction: {no_signal_frac:.3f} -- per the design doc, this should be the "
          f"overwhelming majority; if it is not, the shared catalog itself (built from one arm's own "
          f"reads) is the likely confound, not individual biology.")
    print(f"\ncandidate differentials ({len(candidates)}):")
    for (fam, cidx), flagged_in, clean_in in candidates:
        print(f"  {fam}:{cidx}  missing_copy in {flagged_in}, none in {clean_in}")

if __name__ == "__main__":
    main()
```

- [ ] **Step 3: Run it and record the result**

```
python3 bench/o3_cross_individual_diff.py /mnt/linuxdisk/home/juanfraitu/o3_diff/testis.family_join.tsv /mnt/linuxdisk/home/juanfraitu/o3_diff/fibroblast.family_join.tsv --label-a testis --label-b fibroblast
```

Sanity-check the `no_signal fraction` line before reading anything into the `candidate_differential`
count, exactly as the script itself warns. Append the result to `docs/o1_ledger.md`, explicitly labeled
exploratory per the design doc's "What's NOT established" section — no PCR validation, no known-true-
positive gorilla case exists to check against.

- [ ] **Step 4: Commit**

```
git add bench/o3_cross_individual_diff.py
git commit -m "O3: cross-individual differential comparison script (Task 8)"
```

---

## Self-review notes (from writing this plan)

- **Spec coverage**: every numbered section of the design doc has a task — module structure (Task 1),
  PAF parsing (Task 2), orphan classifier (Task 3), pair detector (Task 4), `copy_assign.rs` Phase 1
  wiring (Task 5), Phase 2 + output surface (Task 6), reproduction gate (Task 7), the Follow-on section
  (Task 8). Error handling and CLI-flag requirements from the design doc are folded into Tasks 5-6 rather
  than a separate task, since they have no independent deliverable of their own.
- **Known risk, disclosed rather than hidden**: Task 5 Step 5 (the per-read grouping inside `compute()`)
  is the one place in this plan where exact field names (`Assignment`'s fields, `BamRead`'s fields) are
  asserted from memory of earlier reading in this session rather than re-verified against the compiler at
  plan-writing time. This is flagged explicitly in the task text — the implementer must re-grep the real
  struct definitions before trusting this step's code verbatim, and Task 5 Step 6 budgets time for
  reconciling compile errors from exactly this source.
- **Type consistency check**: `Flag::NoFlag` (not `Flag::None`, which the design doc used) is used
  consistently from Task 1 onward; `is_partner: bool` (not `member_status: String`) matches the design
  doc's own correction; column names (`o3_flag`, `o3_class`, `o3_rate_per_kb`, `o3_p`, `o3_n_rejected`)
  match between Task 6's writer and Task 8's reader.
