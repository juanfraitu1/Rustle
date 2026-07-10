# Collapse gate — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Admit a *collapsed* single-rep locus as a multi-copy family — copy **number** only by default, reads
certified tied — so DAZ and BPY2 stop returning zero families.

**Architecture:** One new decision point. `detect_and_assign` already calls `colocated_families`, which requires
`>= min_copies` assembled **reps**. Reps that belong to no resulting family are passed to a two-leg gate: leg 1
tests whether the locus is *collapsed* (a binomial tail on ambiguously-placed primary reads), and only then leg 2
counts PSV haplotypes (`split_locus_copies` + `chi_h`). Order matters: running leg 2 alone made a single-copy gene
(TSPYL1) report 12 copies. The multi-rep path is untouched.

**Tech Stack:** Rust (`src/rustle/vg_family/`), existing `poisson_binomial_upper_tail`, `split_locus_copies`,
`chi_h`. Python for the planted simulation (`bench/sim_genome.py`).

**Spec:** `docs/superpowers/specs/2026-07-09-collapse-gate-design.md`

## Global Constraints

- **Collapse first, haplotypes second.** Leg 2 must never run unless leg 1 fires. This is the entire lesson from
  `bench/COLLAPSED_COPY_GATE.md`.
- **Primary records only.** Filter `is_supplementary`; a read is *ambiguous* iff `mapq == 0`. Secondary records
  never reach `bam_reads` as separate entries for this purpose — a chimeric segment is adjacency, not ambiguity.
- **`eps_amb` is never zero.** Estimate per region from uniquely-mappable reps with a Jeffreys prior:
  `eps_amb = (k_bg as f64 + 0.5) / (n_bg as f64 + 1.0)`. If there is **no** background rep, the gate **abstains**.
- **Default emits copy NUMBER, not copies.** `n_copies = χ(H)`; every read `AssignStatus::Tied`; no per-copy
  consensus is materialised. `--admit-collapsed-copies` is opt-in and **must not be implemented until Task 6's
  precondition is met**.
- **`min_copies` applies to χ(H)**, not to the rep count. `χ(H) < min_copies` ⇒ no family.
- Gate OFF (`cfg.collapse_gate == false`) ⇒ output byte-identical to today.
- Significance level: reuse `AssignParams::alpha` (default `1e-3`). No new constant.

## File Structure

- `src/rustle/vg_family/collapse_gate.rs` — **CREATE.** Pure functions: ambiguity counting, `eps_amb` estimation,
  the collapse p-value, and the gate verdict. No I/O, no BAM, fully unit-testable.
- `src/rustle/vg_family/mod.rs` — **MODIFY.** Register the module.
- `src/rustle/vg_family/denovo_pipeline.rs` — **MODIFY.** `DenovoConfig.collapse_gate`; call the gate on
  un-familied reps; build the gated `FamilyAssignment`.
- `src/bin/copy_assign.rs` — **MODIFY.** `--no-collapse-gate`.
- `bench/sim_genome.py` — **MODIFY.** Plant het alleles alongside collapsed copies (Task 5).
- `bench/collapse_gate_validation.py`, `bench/COLLAPSE_GATE_VALIDATION.md` — **CREATE** (Task 5).

---

### Task 1: Ambiguity evidence and the collapse p-value

**Files:**
- Create: `src/rustle/vg_family/collapse_gate.rs`
- Modify: `src/rustle/vg_family/mod.rs`
- Test: in-crate `mod tests` in `collapse_gate.rs`

**Interfaces:**
- Consumes: `crate::vg_family::copy_assign::poisson_binomial_upper_tail(k: usize, probs: &[f64]) -> f64`
  (`pub(crate)`, `copy_assign.rs:99`).
- Produces:
  - `pub struct Ambiguity { pub n: usize, pub k: usize }`
  - `pub fn estimate_eps_amb(bg: Ambiguity) -> Option<f64>`
  - `pub fn collapse_pvalue(obs: Ambiguity, eps_amb: f64) -> f64`

- [ ] **Step 1: Write the failing tests**

Create `src/rustle/vg_family/collapse_gate.rs` with only this test module plus `use` lines:

```rust
//! Collapse gate: is a single-rep locus actually several copies the aligner could not separate?
//!
//! SDA (Vollger 2019) detects a collapse by read-depth excess and only THEN defines PSVs, "requiring sequence
//! coverages consistent with a single-copy locus in order to distinguish PSVs from allelic variants". We ran
//! that second stage alone and a single-copy gene (TSPYL1) reported 12 collapsed copies against DAZ's 3.
//!
//! In DNA a collapse shows as excess depth. In RNA depth is copy number x expression (Clair3-RNA: "the coverage
//! is uneven across genomic regions in RNA-seq"), so depth cannot be the instrument. A collapse shows instead as
//! reads the aligner cannot place uniquely. Measured on GGO: 0 MAPQ-0 primaries across 9449 reads at five
//! single-copy loci, against 19/20 at DAZ2 and 30/34 at TSPY -- and expression-invariant (TSPYL1 has 2151 reads
//! and no ambiguity; DAZ2 has 20 reads and 95%).

/// Ambiguously-placed primary reads (`k`) out of primary reads (`n`) at a locus.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Ambiguity {
    pub n: usize,
    pub k: usize,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn eps_amb_is_never_zero_even_when_no_background_read_is_ambiguous() {
        // The five single-copy controls give 0 ambiguous reads in 9449. The MLE is 0, under which ONE stray
        // MAPQ-0 read would be infinitely significant. Jeffreys keeps it positive.
        let eps = estimate_eps_amb(Ambiguity { n: 9449, k: 0 }).unwrap();
        assert!(eps > 0.0, "eps_amb must be strictly positive, got {eps}");
        assert!((eps - 0.5 / 9450.0).abs() < 1e-12, "Jeffreys: (k + 1/2) / (n + 1)");
    }

    #[test]
    fn eps_amb_abstains_without_background_reads() {
        assert_eq!(estimate_eps_amb(Ambiguity { n: 0, k: 0 }), None, "no background => cannot estimate => abstain");
    }

    #[test]
    fn collapse_pvalue_is_significant_for_daz2_and_not_for_a_clean_locus() {
        let eps = estimate_eps_amb(Ambiguity { n: 9449, k: 0 }).unwrap();
        // DAZ2: 19 of 20 primaries ambiguous.
        let p_daz2 = collapse_pvalue(Ambiguity { n: 20, k: 19 }, eps);
        assert!(p_daz2 < 1e-6, "DAZ2 must be overwhelmingly significant, got {p_daz2}");
        // A clean locus with no ambiguous read cannot be significant.
        let p_clean = collapse_pvalue(Ambiguity { n: 2151, k: 0 }, eps);
        assert!((p_clean - 1.0).abs() < 1e-12, "k = 0 => p = 1, got {p_clean}");
    }

    #[test]
    fn collapse_pvalue_of_a_single_stray_read_is_not_significant_at_alpha() {
        // One ambiguous read in a well-covered locus must NOT clear alpha = 1e-3.
        let eps = estimate_eps_amb(Ambiguity { n: 9449, k: 0 }).unwrap();
        let p = collapse_pvalue(Ambiguity { n: 500, k: 1 }, eps);
        assert!(p > 1e-3, "a single stray MAPQ-0 read must not fire the gate, got {p}");
    }
}
```

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test --lib collapse_gate 2>&1 | tail -5`
Expected: FAIL — `cannot find function estimate_eps_amb` / `collapse_pvalue`.

- [ ] **Step 3: Implement**

Add to `src/rustle/vg_family/collapse_gate.rs`, above the test module:

```rust
use crate::vg_family::copy_assign::poisson_binomial_upper_tail;

/// Background per-read ambiguity rate, under a Jeffreys prior so it is never exactly zero.
///
/// `bg` is pooled over the region's uniquely-mappable reps (reps that are not gate candidates). The controls
/// observe `k = 0`, whose MLE is 0; a zero rate would make a single stray MAPQ-0 read infinitely significant.
/// `None` when there is no background to estimate from — the caller must then ABSTAIN, never fire.
pub fn estimate_eps_amb(bg: Ambiguity) -> Option<f64> {
    if bg.n == 0 {
        return None;
    }
    Some((bg.k as f64 + 0.5) / (bg.n as f64 + 1.0))
}

/// Upper-tail probability of seeing `obs.k` or more ambiguous reads among `obs.n`, if the locus were unique.
///
/// `k ~ Binomial(n, eps_amb)` under the null. Reuses the shipped Poisson-binomial tail with a constant
/// probability vector rather than adding a second distribution implementation.
pub fn collapse_pvalue(obs: Ambiguity, eps_amb: f64) -> f64 {
    if obs.k == 0 || obs.n == 0 {
        return 1.0;
    }
    let probs = vec![eps_amb; obs.n];
    poisson_binomial_upper_tail(obs.k, &probs)
}
```

Register the module — in `src/rustle/vg_family/mod.rs`, add alongside the other `pub mod` lines:

```rust
pub mod collapse_gate;
```

`poisson_binomial_upper_tail` is `pub(crate)`, so it is reachable from a sibling module with no visibility change.

- [ ] **Step 4: Run to verify it passes**

Run: `cargo test --lib collapse_gate 2>&1 | tail -5`
Expected: `test result: ok. 4 passed`

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/collapse_gate.rs src/rustle/vg_family/mod.rs
git commit -m "feat(collapse_gate): ambiguity statistic and binomial collapse p-value"
```

---

### Task 2: The gate verdict — collapse first, haplotypes second

**Files:**
- Modify: `src/rustle/vg_family/collapse_gate.rs`
- Test: in-crate

**Interfaces:**
- Consumes: Task 1's `Ambiguity`, `estimate_eps_amb`, `collapse_pvalue`.
- Consumes: `crate::vg_family::readonly_copy_number::chi_h(copy_alleles: &[Vec<Option<u8>>]) -> usize`
  (`readonly_copy_number.rs:27`).
- Produces:
  - `pub enum CollapseVerdict { Fire { chi_h: usize, p_value: f64 }, NotCollapsed { p_value: f64 }, Abstain(&'static str) }`
  - `pub fn collapse_verdict(obs: Ambiguity, bg: Ambiguity, haplotypes: &[Vec<Option<u8>>], alpha: f64, min_copies: usize) -> CollapseVerdict`

`haplotypes` are the `allele_vector`s of the **identifiable** `CopyIsoform`s at the locus. The caller computes
them; the gate stays pure.

- [ ] **Step 1: Write the failing tests**

Append to the `mod tests` block in `collapse_gate.rs`:

```rust
    /// Allele vectors: two haplotypes that differ at a column conflict, so chi_h counts them separately.
    fn haps(rows: &[&[u8]]) -> Vec<Vec<Option<u8>>> {
        rows.iter().map(|r| r.iter().map(|&b| if b == b'.' { None } else { Some(b) }).collect()).collect()
    }

    /// DAZ: 19/20 ambiguous, background clean, two distinguishable haplotypes.
    #[test]
    fn verdict_fires_on_a_collapsed_locus_with_two_haplotypes() {
        let v = collapse_verdict(
            Ambiguity { n: 20, k: 19 },
            Ambiguity { n: 9449, k: 0 },
            &haps(&[b"ACGT", b"ACGA"]),
            1e-3,
            2,
        );
        match v {
            CollapseVerdict::Fire { chi_h, .. } => assert_eq!(chi_h, 2),
            other => panic!("expected Fire, got {other:?}"),
        }
    }

    /// TSPYL1: a single-copy gene whose reads are ALL uniquely placed. Leg 2 would report many haplotypes
    /// (het alleles, editing, isoform noise) -- leg 1 must stop it before leg 2 is ever consulted.
    #[test]
    fn verdict_rejects_a_unique_locus_however_many_haplotypes_it_reports() {
        let twelve: Vec<Vec<Option<u8>>> =
            (0..12u8).map(|i| vec![Some(b'A' + i), Some(b'C')]).collect();
        let v = collapse_verdict(Ambiguity { n: 2151, k: 0 }, Ambiguity { n: 9449, k: 0 }, &twelve, 1e-3, 2);
        assert!(matches!(v, CollapseVerdict::NotCollapsed { .. }), "unique locus must never fire, got {v:?}");
    }

    #[test]
    fn verdict_abstains_without_a_background_estimate() {
        let v = collapse_verdict(Ambiguity { n: 20, k: 19 }, Ambiguity { n: 0, k: 0 }, &haps(&[b"AC", b"AG"]), 1e-3, 2);
        assert!(matches!(v, CollapseVerdict::Abstain(_)), "no background => abstain, got {v:?}");
    }

    /// min_copies applies to chi(H), not to the rep count: one haplotype is not a family.
    #[test]
    fn verdict_rejects_a_collapse_that_resolves_to_one_haplotype() {
        let v = collapse_verdict(Ambiguity { n: 20, k: 19 }, Ambiguity { n: 9449, k: 0 }, &haps(&[b"ACGT"]), 1e-3, 2);
        assert!(matches!(v, CollapseVerdict::NotCollapsed { .. }) || matches!(v, CollapseVerdict::Abstain(_)),
                "chi_h = 1 < min_copies => no family, got {v:?}");
    }
```

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test --lib collapse_gate 2>&1 | tail -5`
Expected: FAIL — `cannot find type CollapseVerdict`.

- [ ] **Step 3: Implement**

Append to `collapse_gate.rs` (above `mod tests`):

```rust
use crate::vg_family::readonly_copy_number::chi_h;

/// What the gate decided about a locus that has only one assembled rep.
#[derive(Clone, Debug, PartialEq)]
pub enum CollapseVerdict {
    /// Collapsed, and its reads resolve into `chi_h >= min_copies` conflicting haplotypes.
    Fire { chi_h: usize, p_value: f64 },
    /// The locus places its reads unambiguously, or its haplotypes do not reach `min_copies`.
    NotCollapsed { p_value: f64 },
    /// Cannot decide — never fire on an unbounded statistic.
    Abstain(&'static str),
}

/// Two legs, in SDA's order. Leg 2 is NOT consulted unless leg 1 fires: a single-copy gene reports plenty of
/// haplotypes (het alleles are haplotypes), and gating on them alone made TSPYL1 report 12 copies against DAZ's 3.
pub fn collapse_verdict(
    obs: Ambiguity,
    bg: Ambiguity,
    haplotypes: &[Vec<Option<u8>>],
    alpha: f64,
    min_copies: usize,
) -> CollapseVerdict {
    // leg 1 -- is this locus collapsed at all?
    let Some(eps_amb) = estimate_eps_amb(bg) else {
        return CollapseVerdict::Abstain("no uniquely-mappable rep to estimate the background ambiguity rate");
    };
    let p_value = collapse_pvalue(obs, eps_amb);
    if p_value >= alpha {
        return CollapseVerdict::NotCollapsed { p_value };
    }
    // leg 2 -- and how many copies collapsed? chi(H) is a LOWER BOUND: two copies x two alleles also yields
    // four haplotypes. Hence the default output is a copy NUMBER with reads tied, never an assignment.
    let chi = chi_h(haplotypes);
    if chi < min_copies {
        return CollapseVerdict::NotCollapsed { p_value };
    }
    CollapseVerdict::Fire { chi_h: chi, p_value }
}
```

- [ ] **Step 4: Run to verify it passes**

Run: `cargo test --lib collapse_gate 2>&1 | tail -5`
Expected: `test result: ok. 8 passed`

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/collapse_gate.rs
git commit -m "feat(collapse_gate): two-leg verdict, collapse before haplotypes"
```

---

### Task 3: Wire the gate into `detect_and_assign`

**Files:**
- Modify: `src/rustle/vg_family/denovo_pipeline.rs`
- Test: in-crate

**Interfaces:**
- Consumes: Task 2's `collapse_verdict`, `Ambiguity`, `CollapseVerdict`.
- Consumes: `crate::vg_family::copy_split::split_locus_copies(reads: &[AlignedRead], min_allele_reads: usize,
  min_psv_k: usize, min_reads_per_copy: usize) -> Vec<CopyIsoform>`; `CopyIsoform { allele_vector: Vec<Option<u8>>,
  identifiable: bool, read_count: usize, .. }` (`copy_split.rs:28`).
- Consumes: `BamRead { chrom, read: AlignedRead, mapq: u8, is_supplementary: bool, .. }`
  (`denovo_assemble.rs:257`); `read_ref_end(&AlignedRead) -> u64` (already used by `recover_collapsed_candidates`).
- Produces: `DenovoConfig.collapse_gate: bool` (default `true`);
  `fn locus_ambiguity(rep: &DenovoTranscript, bam_reads: &[BamRead]) -> Ambiguity`;
  `fn gated_family(rep, reads_idx, chi_h, family_id) -> FamilyAssignment`.

**Where.** In `detect_and_assign`, immediately after the `for cf in colocated_families(...)` loop finishes and
before `out` is returned. Reps already covered by an emitted family are skipped.

- [ ] **Step 1: Write the failing test**

Add to `mod tests` in `denovo_pipeline.rs`:

```rust
    /// A locus whose reads are ALL uniquely placed contributes no ambiguity, whatever its haplotypes.
    #[test]
    fn locus_ambiguity_counts_only_mapq0_primaries_overlapping_the_rep() {
        let rep = rep_s(1_000, 2_000, vec![], 10);
        let mk = |start: u64, mapq: u8, supp: bool| BamRead {
            chrom: "c1".into(),
            read: AlignedRead { ref_start: start, cigar: vec![('M', 100)], seq: vec![b'A'; 100], qual: vec![30; 100] },
            mapq,
            name: format!("r{start}_{mapq}"),
            as_score: 100,
            de: 0.01,
            is_supplementary: supp,
        };
        let reads = vec![
            mk(1_100, 0, false),   // ambiguous, inside  -> counts
            mk(1_200, 60, false),  // unique,    inside  -> denominator only
            mk(1_300, 0, true),    // supplementary      -> ignored entirely
            mk(9_000, 0, false),   // outside the rep    -> ignored entirely
        ];
        let a = locus_ambiguity(&rep, &reads);
        assert_eq!(a, crate::vg_family::collapse_gate::Ambiguity { n: 2, k: 1 });
    }

    #[test]
    fn denovoconfig_default_enables_the_collapse_gate() {
        assert!(DenovoConfig::default().collapse_gate);
    }
```

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test --lib locus_ambiguity 2>&1 | tail -5`
Expected: FAIL — `cannot find function locus_ambiguity`.

- [ ] **Step 3: Implement**

Add `collapse_gate: bool` to `DenovoConfig` (doc comment: *"Admit a collapsed single-rep locus as a multi-copy
family: copy NUMBER only, reads certified tied. Default on."*) and `collapse_gate: true` to its `Default` impl.

Add near `recover_collapsed_candidates` in `denovo_pipeline.rs`:

```rust
use crate::vg_family::collapse_gate::{collapse_verdict, Ambiguity, CollapseVerdict};

/// Ambiguously-placed primary reads over a rep's span. A read is ambiguous iff `mapq == 0`: the aligner found
/// no reason to prefer this placement. Supplementary records are excluded for the same reason they are excluded
/// from conflict edges — a chimeric segment is adjacency, not ambiguity.
fn locus_ambiguity(rep: &DenovoTranscript, bam_reads: &[BamRead]) -> Ambiguity {
    let mut n = 0usize;
    let mut k = 0usize;
    for br in bam_reads {
        if br.is_supplementary || br.chrom != rep.chrom {
            continue;
        }
        if br.read.ref_start < rep.end && read_ref_end(&br.read) > rep.start {
            n += 1;
            if br.mapq == 0 {
                k += 1;
            }
        }
    }
    Ambiguity { n, k }
}
```

Then, after the `colocated_families` loop, add the gate. `familied` is the set of rep indices already emitted;
`bg` pools ambiguity over reps that are NOT gate candidates (the uniquely-mappable background):

```rust
    if cfg.collapse_gate {
        let familied: std::collections::HashSet<&str> =
            out.iter().flat_map(|fa: &FamilyAssignment| fa.copy_tids.iter().map(|s| s.as_str())).collect();
        let candidates: Vec<&DenovoTranscript> = reps.iter().filter(|r| !familied.contains(r.tid.as_str())).collect();
        // background = every rep that is NOT a candidate, pooled
        let bg = reps
            .iter()
            .filter(|r| familied.contains(r.tid.as_str()))
            .map(|r| locus_ambiguity(r, bam_reads))
            .fold(Ambiguity { n: 0, k: 0 }, |a, b| Ambiguity { n: a.n + b.n, k: a.k + b.k });

        for rep in candidates {
            let obs = locus_ambiguity(rep, bam_reads);
            let reads: Vec<AlignedRead> = bam_reads
                .iter()
                .filter(|br| !br.is_supplementary && br.chrom == rep.chrom
                    && br.read.ref_start < rep.end && read_ref_end(&br.read) > rep.start)
                .map(|br| br.read.clone())
                .collect();
            // leg 2's input, computed lazily-enough: split_locus_copies is cheap next to assignment, and the
            // verdict discards it unless leg 1 fired.
            let haplotypes: Vec<Vec<Option<u8>>> = split_locus_copies(&reads, 3, 2, 3)
                .into_iter()
                .filter(|c| c.identifiable)
                .map(|c| c.allele_vector)
                .collect();
            match collapse_verdict(obs, bg, &haplotypes, params.alpha, min_copies) {
                CollapseVerdict::Fire { chi_h, p_value } => {
                    eprintln!(
                        "[detect_and_assign] collapse gate: {}:{}-{} ambiguous {}/{} p={:.2e} -> chi(H)={} copies \
                         (reads certified TIED; chi(H) is a LOWER BOUND)",
                        rep.chrom, rep.start, rep.end, obs.k, obs.n, p_value, chi_h
                    );
                    out.push(gated_family(rep, bam_reads, chi_h, format!("DSFAM{}", out.len())));
                }
                CollapseVerdict::Abstain(why) => eprintln!(
                    "[detect_and_assign] collapse gate ABSTAINS at {}:{}-{}: {why}",
                    rep.chrom, rep.start, rep.end
                ),
                CollapseVerdict::NotCollapsed { .. } => {}
            }
        }
    }
```

And the gated-family constructor. It does **not** go through the assignment pipeline: with no materialised copies
there is nothing to assign to, and that is the point — every read is `Tied` with `min_p_value = 1.0`, exactly the
certificate that says "no distinguishing column was available to this read".

```rust
/// A family whose copies were never assembled: `n_copies` is chi(H), every read is certified TIED, and no
/// per-copy consensus exists. `--admit-collapsed-copies` (Task 6) is what would materialise them.
fn gated_family(rep: &DenovoTranscript, bam_reads: &[BamRead], chi_h: usize, family_id: String) -> FamilyAssignment {
    let assignments: Vec<(usize, Assignment)> = bam_reads
        .iter()
        .enumerate()
        .filter(|(_, br)| !br.is_supplementary && br.chrom == rep.chrom
            && br.read.ref_start < rep.end && read_ref_end(&br.read) > rep.start)
        .map(|(i, _)| (i, Assignment {
            best_copy: 0,
            status: AssignStatus::Tied,
            log_lr_margin: 0.0,
            n_decisive: 0,
            p_value: 1.0,
            min_p_value: 1.0,
            discovery_coupled: false,
            posterior: vec![1.0 / chi_h as f64; chi_h],
        }))
        .collect();
    FamilyAssignment {
        family_id,
        chrom: rep.chrom.clone(),
        n_copies: chi_h,
        n_reads: assignments.len(),
        collapsed_copies: chi_h,          // records HOW the count was obtained
        assignments,
        copy_tids: vec![rep.tid.clone()], // the one assembled rep; the rest have no sequence
        copy_spans: vec![(rep.chrom.clone(), rep.start, rep.end)],
        ..FamilyAssignment::empty()
    }
}
```

If `FamilyAssignment` has no `empty()` constructor, add one in the same commit that returns all-zero/empty fields;
do not spread `Default` onto a struct with meaningful invariants.

- [ ] **Step 4: Run to verify it passes**

Run: `cargo test --lib 2>&1 | tail -3`
Expected: all tests pass (853 + the 2 new + Task 1/2's 8).

- [ ] **Step 5: Verify DAZ and the controls on real data (FOREGROUND, one region at a time)**

⚠ **Crash rule:** `copy_assign` foreground, serial, small batches, outputs under `/home/juanfra/winloci_scratch`.
No `nohup`, no background waiters, no `pkill -f`.

```bash
cd /home/juanfra/winloci_scratch
CA=/mnt/c/Users/jfris/Desktop/Rustle/target/release/copy_assign
timeout 500 $CA --bam GGO_mm.bam --fasta GGO.fasta --region NC_073248.2:42778133-42950552 \
    --min-copies 2 --skip-poa-diagnostic --homology-primary --out cg_DAZ
timeout 500 $CA --bam GGO_mm.bam --fasta GGO.fasta --region NC_073229.2:140471252-140485905 \
    --min-copies 2 --skip-poa-diagnostic --out cg_TSPYL1
```
Expected: DAZ prints `collapse gate: ... ambiguous 22/200 ... -> chi(H)=N copies` and emits **1 family** whose
reads are all `tied`. TSPYL1 emits **0 families** and prints nothing from the gate.

- [ ] **Step 6: Commit**

```bash
git add src/rustle/vg_family/denovo_pipeline.rs
git commit -m "feat(o2): collapse gate admits a collapsed single-rep locus as chi(H) copies, reads tied"
```

---

### Task 4: CLI opt-out and byte-identical OFF path

**Files:**
- Modify: `src/bin/copy_assign.rs`

**Interfaces:**
- Consumes: `DenovoConfig.collapse_gate` from Task 3.

- [ ] **Step 1: Add the flag**

In the `Args` struct, next to `keep_readthrough`:

```rust
    /// Disable the collapse gate. By default a single-rep locus whose reads are ambiguously placed at a rate
    /// incompatible with a unique locus is admitted as a multi-copy family with `n_copies = chi(H)` and its
    /// reads certified tied (no per-copy sequence is materialised, so no assignment can be wrong). Pass this to
    /// reproduce the pre-gate behaviour exactly.
    #[arg(long, default_value_t = false)]
    no_collapse_gate: bool,
```

And where the other config fields are set:

```rust
    cfg.collapse_gate = !args.no_collapse_gate;
```

- [ ] **Step 2: Verify the OFF path is byte-identical**

```bash
cd /home/juanfra/winloci_scratch
CA=/mnt/c/Users/jfris/Desktop/Rustle/target/release/copy_assign
R=NC_073224.2:101578582-101607889
timeout 400 $CA --bam GGO_mm.bam --fasta GGO.fasta --region $R --min-copies 2 --skip-poa-diagnostic \
    --no-collapse-gate --out cg_off
md5sum cg_off.assignments.tsv   # must equal the pre-gate run's md5
```
Expected: identical md5 to a build from commit `d0144e6` on the same region.

- [ ] **Step 3: Commit**

```bash
cargo build --release --bin copy_assign
git add src/bin/copy_assign.rs
git commit -m "feat(copy_assign): --no-collapse-gate"
```

---

### Task 5: The discriminating experiment — 2 copies × 2 alleles vs 4 copies

**Files:**
- Modify: `bench/sim_genome.py`
- Create: `bench/collapse_gate_validation.py`, `bench/COLLAPSE_GATE_VALIDATION.md`

This task decides whether Task 6 is allowed to exist. **A negative result is the deliverable if that is what the
data says** — write it up either way.

- [ ] **Step 1: Plant het alleles**

`bench/sim_genome.py` already has `plant_collapsed(fam, chrom, n_copies, n_psv, nr)` ("the regime the gate is
FOR"). Add a sibling that plants a **diploid** version — each copy carries two alleles differing at `n_het`
columns that are *not* PSV columns:

```python
def plant_collapsed_diploid(fam, chrom, n_copies, n_psv, n_het, nr):
    """n_copies collapsed copies, EACH with two het alleles.

    The gate cannot tell a haplotype from a copy: n_copies x 2 alleles yields 2*n_copies haplotypes. Planting
    both arms with matched read counts and matched column counts is the only way to ask whether chi(H) recovers
    n_copies or 2*n_copies. Returns the truth record for scoring.
    """
    # (implementation mirrors plant_collapsed; after choosing the per-copy PSV allele vector, choose n_het
    # additional columns and give each copy TWO reads-generating haplotypes differing only there)
```

Then plant two loci with **matched total reads and matched column counts**:
- locus `A`: `plant_collapsed(n_copies=4, n_psv=8, nr=200)` → truth `chi_h = 4`
- locus `B`: `plant_collapsed_diploid(n_copies=2, n_psv=8, n_het=8, nr=200)` → truth `chi_h = 2`, haplotypes `4`

- [ ] **Step 2: Score**

`bench/collapse_gate_validation.py` runs `copy_assign` on the planted BAM (foreground, serial) and reports, per
planted locus: did the gate fire, what `chi(H)` did it emit, and what was the truth. Also plant 5 unique loci and
report the gate's false-positive count.

- [ ] **Step 3: Write it up**

`bench/COLLAPSE_GATE_VALIDATION.md`: a table of planted truth vs emitted `chi(H)`; the false-positive rate on
planted unique loci; and an explicit verdict sentence — **"chi(H) distinguishes 2 copies × 2 alleles from 4
copies: YES / NO"**. If NO, state that `--admit-collapsed-copies` must not be implemented, and why: an admitted
"copy" would be a haplotype.

- [ ] **Step 4: Commit**

```bash
git add bench/sim_genome.py bench/collapse_gate_validation.py bench/COLLAPSE_GATE_VALIDATION.md
git commit -m "bench(collapse_gate): planted validation, 2 copies x 2 alleles vs 4 copies"
```

---

### Task 6: `--admit-collapsed-copies` — CONDITIONAL on Task 5

**Do not start this task unless Task 5 reported YES.** If Task 5 reported NO, close it out by recording in
`bench/COLLAPSE_GATE_VALIDATION.md` that copy materialisation is refused, and stop. That is a complete outcome.

If YES: materialise per-copy consensus from the identifiable `CopyIsoform`s (reuse `CollapsedCandidate`, which
`recover_collapsed_candidates` already builds and which carries `iso` + `psv_pos`), add them to the family's copy
set, and let the existing assignment path run. Add `--admit-collapsed-copies`, default off. Re-run Task 3's DAZ
check and report how many of DAZ's 220 reads become assigned.

---

## Self-Review

- **Spec coverage.** Ambiguity instrument → T1. Jeffreys `eps_amb`, abstain-without-background → T1, T2.
  Collapse-before-haplotypes ordering → T2 (enforced by control flow, tested by
  `verdict_rejects_a_unique_locus_however_many_haplotypes_it_reports`). `min_copies` on χ(H) → T2. χ(H)-only
  default with tied reads → T3 (`gated_family`). Byte-identical OFF → T4. Planted validation incl. the
  copy-vs-allele experiment → T5. Opt-in materialisation, licensed by T5 → T6. Non-goal (reference-absent) is
  not implemented anywhere, as intended.
- **Type consistency.** `Ambiguity { n, k }` produced T1, consumed T2/T3. `CollapseVerdict` produced T2, matched
  T3. `chi_h(&[Vec<Option<u8>>])` matches `CopyIsoform::allele_vector: Vec<Option<u8>>`.
- **Known unknown, flagged for the implementer.** `FamilyAssignment` has many fields; T3 uses
  `..FamilyAssignment::empty()` and instructs adding that constructor if absent. Verify the field list against
  `denovo_pipeline.rs:296` before writing `gated_family` — do not guess field names.
- **Risk.** The gate's background `bg` is pooled over *familied* reps. In a region where every rep is a gate
  candidate, `bg.n == 0` and the gate abstains — correct, and tested.
