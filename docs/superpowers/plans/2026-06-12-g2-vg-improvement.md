# `-G2` Annotation Flag + VG-Improvement Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `-G2` annotation flag that defines paralog families (position-agnostic, gated by a tunable shared-exon similarity threshold) to guide VG copy recovery, plus a per-run head-to-head scorecard vs StringTie.

**Architecture:** A new `annotation_families` module loads `-G2`, extracts each transcript's exon sequence from the `--genome-fasta` `GenomeIndex`, and clusters transcripts into paralog families using rustle's existing **sequence-only** similarity primitives (`refine_by_minimizer_jaccard` / `merge_singletons_by_sequence`) — never the same-chrom `cluster_by_position`. When `-G2` is present, VG family discovery is sourced from these families; the existing copy-partition/phasing/EM machinery (already position-agnostic) runs unchanged. Two TSV reports are emitted: `vg_families.tsv` (per-family achieved similarity + per-copy chroms) and `vg_eval.tsv` (per-copy WIN/tie/miss/regression vs `-G`).

**Tech Stack:** Rust; `clap` (CLI), existing modules `reference_gtf`, `genome::GenomeIndex`, `vg_family::family_graph` (similarity primitives), `vg` (family discovery + EM).

**Spec:** `docs/superpowers/specs/2026-06-12-g2-vg-improvement-design.md`

**Invariants (every task):** `-G2` absent ⇒ output byte-identical to current behavior (`RUSTLE_PRECISE=1 bash bench/mini3/check_precise.sh` must print "byte-matches 4705ab1"). Never `git add -A` (the tree carries an unrelated VG-family source track); stage files explicitly. Commit messages end `Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>`.

---

## File Structure

- **Create** `src/rustle/annotation_families.rs` — load `-G2`, extract exon sequences via `GenomeIndex`, cluster into position-agnostic paralog families, define `AnnotationFamily`/`CopyStructure`. One responsibility: annotation → families.
- **Create** `src/rustle/vg_eval.rs` — the two report writers (`write_vg_families_tsv`, `write_vg_eval_tsv`) + chain-match helper. One responsibility: reporting/scoring.
- **Modify** `src/rustle/main.rs` — add `--guide2`/`-G2` and `--family-exon-similarity` CLI flags; `--vg` constraint; thread into `RunConfig`.
- **Modify** `src/rustle/types.rs` — `RunConfig` fields `guide2_path: Option<String>`, `family_exon_similarity: f64`.
- **Modify** `src/rustle/lib.rs` — `pub mod annotation_families; pub mod vg_eval;`.
- **Modify** `src/rustle/vg.rs` — at the family-discovery entry, when `guide2_path.is_some()` build `FamilyGroup`s from `annotation_families` output instead of `discover_family_groups`.
- **Create** `tests/regression/g2_annotation_families.rs` — clustering, trans-chromosomal, far-apart, gate, report, eval, read-based-linker regression.

---

## Task 1: CLI flags + RunConfig wiring

**Files:**
- Modify: `src/rustle/main.rs` (clap struct ~line 44 area; RunConfig build ~line 780-828)
- Modify: `src/rustle/types.rs` (RunConfig struct ~line 818-821; Default ~line 1220)

- [ ] **Step 1: Add the failing test (config defaults + constraint)**

In `tests/regression/g2_annotation_families.rs` (new file), add:

```rust
use rustle::types::RunConfig;

#[test]
fn runconfig_g2_defaults_off() {
    let c = RunConfig::default();
    assert!(c.guide2_path.is_none(), "guide2_path must default to None");
    assert!((c.family_exon_similarity - rustle::vg_family::family_graph::family_merge_jaccard()).abs() < 1e-9,
        "family_exon_similarity default must equal family_merge_jaccard()");
}
```

- [ ] **Step 2: Run it — expect FAIL** (`cargo test --release --test g2_annotation_families runconfig_g2_defaults_off`) with "no field `guide2_path`".

- [ ] **Step 3: Add the RunConfig fields**

In `src/rustle/types.rs`, in the `RunConfig` struct near `pub guide_mode: bool` (~line 821):

```rust
    /// Path to the -G2 reference annotation (drives VG family definition + scorecard). None = off.
    pub guide2_path: Option<String>,
    /// Shared-exon minimizer-Jaccard threshold gating family-graph merging (--family-exon-similarity).
    pub family_exon_similarity: f64,
```

In the `Default` impl (~line 1220, alongside `guide_mode: false`):

```rust
            guide2_path: None,
            family_exon_similarity: crate::vg_family::family_graph::family_merge_jaccard(),
```

- [ ] **Step 4: Add the clap args**

In `src/rustle/main.rs`, after the `-G` arg (`guide: Option<String>` ~line 45):

```rust
    /// Second reference annotation (-G2): defines paralog families that guide VG copy
    /// recovery + the per-run scorecard vs -G. Requires --vg and --genome-fasta. GTF/GFF.
    #[arg(long = "guide2")]
    guide2: Option<String>,

    /// Shared-exon minimizer-Jaccard threshold (0..1) gating whether two copies merge into
    /// one family graph. Default = the built-in family_merge_jaccard bar.
    #[arg(long = "family-exon-similarity")]
    family_exon_similarity: Option<f64>,
```

- [ ] **Step 5: Thread into RunConfig + add the --vg / --genome-fasta constraint**

In `src/rustle/main.rs` where `RunConfig` is built (~line 780-828), add:

```rust
        guide2_path: args.guide2.clone(),
        family_exon_similarity: args.family_exon_similarity
            .unwrap_or_else(rustle::vg_family::family_graph::family_merge_jaccard),
```

And immediately after parsing args (before the run), add the constraint:

```rust
    if args.guide2.is_some() && !args.vg {
        eprintln!("error: --guide2/-G2 requires --vg");
        std::process::exit(2);
    }
    if args.guide2.is_some() && args.genome_fasta.is_none() {
        eprintln!("error: --guide2/-G2 requires --genome-fasta (exon sequences for family clustering)");
        std::process::exit(2);
    }
```

- [ ] **Step 6: Run the test — expect PASS** (`cargo test --release --test g2_annotation_families runconfig_g2_defaults_off`).

- [ ] **Step 7: Escape-hatch invariant**

Run: `RUSTLE_PRECISE=1 bash bench/mini3/check_precise.sh`
Expected: "ESCAPE-HATCH OK: RUSTLE_PRECISE=1 byte-matches 4705ab1".

- [ ] **Step 8: Commit**

```bash
git add src/rustle/types.rs src/rustle/main.rs tests/regression/g2_annotation_families.rs
git commit -m "feat(g2): add -G2/--guide2 + --family-exon-similarity flags + RunConfig"
```

---

## Task 2: `annotation_families` data types + loader

**Files:**
- Create: `src/rustle/annotation_families.rs`
- Modify: `src/rustle/lib.rs` (add `pub mod annotation_families;`)
- Test: `tests/regression/g2_annotation_families.rs`

- [ ] **Step 1: Failing test (load → copies)**

```rust
use rustle::annotation_families::{load_copies, CopyStructure};

#[test]
fn load_copies_parses_two_transcripts() {
    let gtf = "chr1\ttest\texon\t101\t200\t.\t+\t.\ttranscript_id \"A\";\n\
               chr5\ttest\texon\t301\t400\t.\t+\t.\ttranscript_id \"B\";\n";
    let path = std::env::temp_dir().join("g2_load.gtf");
    std::fs::write(&path, gtf).unwrap();
    let copies = load_copies(&path).unwrap();
    assert_eq!(copies.len(), 2);
    assert_eq!(copies[0].copy_id, "A");
    assert_eq!(copies[0].chrom, "chr1");
    assert_eq!(copies[1].chrom, "chr5");
    // 1-based inclusive [101,200] -> 0-based [100,200)
    assert_eq!(copies[0].exons, vec![(100, 200)]);
}
```

- [ ] **Step 2: Run — expect FAIL** ("unresolved import `rustle::annotation_families`").

- [ ] **Step 3: Implement the module skeleton + loader**

Create `src/rustle/annotation_families.rs`:

```rust
//! -G2 annotation -> position-agnostic paralog families.
//! Clusters annotation transcripts into families by SHARED-EXON SEQUENCE similarity
//! (never coordinate proximity), gated by RunConfig.family_exon_similarity.

use anyhow::Result;
use std::path::Path;

/// One annotated transcript = one copy in a (potential) family.
#[derive(Clone, Debug)]
pub struct CopyStructure {
    pub copy_id: String,
    pub chrom: String,
    pub strand: char,
    /// 0-based half-open exons.
    pub exons: Vec<(u64, u64)>,
}

/// A paralog family: >=2 copies merged by shared-exon similarity, with the achieved similarity.
#[derive(Clone, Debug)]
pub struct AnnotationFamily {
    pub family_id: String,
    pub copies: Vec<CopyStructure>,
    pub achieved_min_sim: f64,
    pub achieved_mean_sim: f64,
}

/// Load -G2 GTF into copy structures (one per transcript).
pub fn load_copies<P: AsRef<Path>>(path: P) -> Result<Vec<CopyStructure>> {
    let refs = crate::reference_gtf::parse_reference_gtf(path)?;
    Ok(refs.into_iter()
        .map(|r| CopyStructure { copy_id: r.id, chrom: r.chrom, strand: r.strand, exons: r.exons })
        .collect())
}
```

In `src/rustle/lib.rs`, near the other `pub mod` lines, add:

```rust
pub mod annotation_families; // -G2 annotation -> paralog families
```

- [ ] **Step 4: Run — expect PASS** (`cargo test --release --test g2_annotation_families load_copies_parses_two_transcripts`).

- [ ] **Step 5: Commit**

```bash
git add src/rustle/annotation_families.rs src/rustle/lib.rs tests/regression/g2_annotation_families.rs
git commit -m "feat(g2): annotation_families module + -G2 copy loader"
```

---

## Task 3: Exon-sequence extraction from `GenomeIndex`

**Files:**
- Modify: `src/rustle/annotation_families.rs`
- Test: `tests/regression/g2_annotation_families.rs`

**Verify first:** read `src/rustle/genome.rs` for the `GenomeIndex` accessor that returns a subsequence for `(chrom, start, end)` (used by `discover_sequence_similar_bundles` at `vg.rs:283`). Match its exact signature; the call below assumes `genome.fetch_seq(chrom, start, end) -> Option<Vec<u8>>` — replace with the real method name found in `genome.rs`.

- [ ] **Step 1: Failing test (concatenate exon seqs)**

```rust
use rustle::annotation_families::copy_sequence;
use rustle::genome::GenomeIndex;

#[test]
fn copy_sequence_concatenates_exons() {
    // Build a tiny in-memory genome fixture per genome.rs test helpers, chr1 = "ACGTACGTAC".
    let genome = GenomeIndex::from_pairs(vec![("chr1".to_string(), b"ACGTACGTAC".to_vec())]);
    let copy = rustle::annotation_families::CopyStructure {
        copy_id: "A".into(), chrom: "chr1".into(), strand: '+',
        exons: vec![(0, 4), (6, 10)], // ACGT + CGTAC[6..10]=CTAC
    };
    let seq = copy_sequence(&copy, &genome).unwrap();
    assert_eq!(seq, b"ACGTGTAC");
}
```
> Note: adjust `GenomeIndex::from_pairs` to the actual test constructor in `genome.rs`. If none exists, write the fixture via a temp FASTA + the real loader.

- [ ] **Step 2: Run — expect FAIL** ("no function `copy_sequence`").

- [ ] **Step 3: Implement `copy_sequence`**

In `src/rustle/annotation_families.rs`:

```rust
/// Concatenated exon sequence of a copy, in transcription order (uses genome subsequence fetch).
/// Returns None if any exon coordinate is out of range for its chromosome.
pub fn copy_sequence(copy: &CopyStructure, genome: &crate::genome::GenomeIndex) -> Option<Vec<u8>> {
    let mut seq = Vec::new();
    for &(s, e) in &copy.exons {
        // Replace fetch_seq with the verified GenomeIndex accessor from genome.rs.
        let part = genome.fetch_seq(&copy.chrom, s, e)?;
        seq.extend_from_slice(&part);
    }
    Some(seq)
}
```

- [ ] **Step 4: Run — expect PASS.**

- [ ] **Step 5: Commit**

```bash
git add src/rustle/annotation_families.rs tests/regression/g2_annotation_families.rs
git commit -m "feat(g2): extract copy exon sequences from GenomeIndex"
```

---

## Task 4: Position-agnostic clustering (sequence-only)

**Files:**
- Modify: `src/rustle/annotation_families.rs`
- Test: `tests/regression/g2_annotation_families.rs`

This is the methodological core. Clustering MUST use sequence-only similarity (minimizer-Jaccard), never `cluster_by_position`. Reuse `vg_family::family_graph::minimizers` (`family_graph.rs:226`).

- [ ] **Step 1: Failing test — TRANS-CHROMOSOMAL + far-apart + gate**

```rust
use rustle::annotation_families::{cluster_families, CopyStructure};

fn copy(id: &str, chrom: &str, exons: Vec<(u64,u64)>) -> CopyStructure {
    CopyStructure { copy_id: id.into(), chrom: chrom.into(), strand: '+', exons }
}

#[test]
fn trans_chromosomal_copies_form_one_family() {
    // Two copies on DIFFERENT chromosomes with near-identical sequence.
    let seqs = vec![
        (copy("A", "chr1", vec![(0,30)]), b"ACGTACGTACGTACGTACGTACGTACGTAC".to_vec()),
        (copy("B", "chr5", vec![(0,30)]), b"ACGTACGTACGTACGTACGTACGTACGTAG".to_vec()), // 1 mismatch
    ];
    let fams = cluster_families(seqs, 0.5);
    assert_eq!(fams.len(), 1, "trans-chromosomal homologs must form ONE family");
    assert_eq!(fams[0].copies.len(), 2);
    let chroms: Vec<_> = fams[0].copies.iter().map(|c| c.chrom.as_str()).collect();
    assert!(chroms.contains(&"chr1") && chroms.contains(&"chr5"));
}

#[test]
fn far_apart_same_chrom_form_one_family() {
    let seqs = vec![
        (copy("A", "chr1", vec![(0,30)]), b"ACGTACGTACGTACGTACGTACGTACGTAC".to_vec()),
        (copy("B", "chr1", vec![(20_000_000,20_000_030)]), b"ACGTACGTACGTACGTACGTACGTACGTAC".to_vec()),
    ];
    let fams = cluster_families(seqs, 0.5);
    assert_eq!(fams.len(), 1, ">10Mb-apart homologs must form ONE family");
}

#[test]
fn dissimilar_copies_stay_separate_and_threshold_is_real() {
    let seqs = vec![
        (copy("A", "chr1", vec![(0,30)]), b"ACGTACGTACGTACGTACGTACGTACGTAC".to_vec()),
        (copy("B", "chr5", vec![(0,30)]), b"TTTTGGGGCCCCAAAATTTTGGGGCCCCAA".to_vec()),
    ];
    assert_eq!(cluster_families(seqs.clone(), 0.5).len(), 2, "dissimilar -> 2 families");
    // Raise T above a similar pair's similarity -> they split too (gate is real).
    let sim = vec![
        (copy("A", "chr1", vec![(0,30)]), b"ACGTACGTACGTACGTACGTACGTACGTAC".to_vec()),
        (copy("B", "chr5", vec![(0,30)]), b"ACGTACGTACGTACGTTTTTACGTACGTAC".to_vec()),
    ];
    assert_eq!(cluster_families(sim.clone(), 0.99).len(), 2, "T above achieved sim -> split");
}
```

- [ ] **Step 2: Run — expect FAIL** ("no function `cluster_families`").

- [ ] **Step 3: Implement `cluster_families` (union-find by minimizer-Jaccard, no position key)**

In `src/rustle/annotation_families.rs`:

```rust
use std::collections::HashSet;

fn jaccard(a: &HashSet<u64>, b: &HashSet<u64>) -> f64 {
    if a.is_empty() && b.is_empty() { return 1.0; }
    let inter = a.intersection(b).count() as f64;
    let union = a.union(b).count() as f64;
    if union == 0.0 { 0.0 } else { inter / union }
}

/// Cluster copies into families purely by shared-exon minimizer-Jaccard >= threshold.
/// POSITION-AGNOSTIC: no chrom / coordinate / distance gate. Same-strand is the only
/// structural precondition (enforced upstream where families feed build_family_graph).
pub fn cluster_families(copies_with_seq: Vec<(CopyStructure, Vec<u8>)>, threshold: f64) -> Vec<AnnotationFamily> {
    let n = copies_with_seq.len();
    // Minimizer set per copy (k=15, w=10 mirror refine_by_minimizer_jaccard defaults).
    let mins: Vec<HashSet<u64>> = copies_with_seq.iter()
        .map(|(_, s)| crate::vg_family::family_graph::minimizers(s, 15, 10))
        .collect();
    // Union-find by pairwise similarity (NO position key).
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(p: &mut Vec<usize>, x: usize) -> usize { if p[x]!=x { let r=find(p,p[x]); p[x]=r; } p[x] }
    let mut sim = vec![vec![0.0f64; n]; n];
    for i in 0..n { for j in (i+1)..n {
        let s = jaccard(&mins[i], &mins[j]);
        sim[i][j] = s; sim[j][i] = s;
        if s >= threshold { let (ri,rj)=(find(&mut parent,i),find(&mut parent,j)); if ri!=rj { parent[ri]=rj; } }
    }}
    // Group, compute achieved similarity, drop singletons (a family needs >=2 copies).
    use std::collections::BTreeMap;
    let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for i in 0..n { let r = find(&mut parent, i); groups.entry(r).or_default().push(i); }
    let mut out = Vec::new();
    for (gi, (_, members)) in groups.iter().enumerate() {
        if members.len() < 2 { continue; }
        let mut pair_sims = Vec::new();
        for a in 0..members.len() { for b in (a+1)..members.len() { pair_sims.push(sim[members[a]][members[b]]); } }
        let min_s = pair_sims.iter().cloned().fold(f64::INFINITY, f64::min);
        let mean_s = pair_sims.iter().sum::<f64>() / pair_sims.len() as f64;
        out.push(AnnotationFamily {
            family_id: format!("FAM{}", gi),
            copies: members.iter().map(|&m| copies_with_seq[m].0.clone()).collect(),
            achieved_min_sim: min_s,
            achieved_mean_sim: mean_s,
        });
    }
    out
}
```

- [ ] **Step 4: Run — expect PASS (all four cluster tests).**

- [ ] **Step 5: Commit**

```bash
git add src/rustle/annotation_families.rs tests/regression/g2_annotation_families.rs
git commit -m "feat(g2): position-agnostic family clustering by shared-exon minimizer-Jaccard"
```

---

## Task 5: `build_families` entry (load + seq + cluster, same-strand split)

**Files:**
- Modify: `src/rustle/annotation_families.rs`
- Test: `tests/regression/g2_annotation_families.rs`

- [ ] **Step 1: Failing test (opposite strands never merge)**

```rust
#[test]
fn opposite_strand_copies_never_merge() {
    use rustle::annotation_families::{cluster_families, CopyStructure};
    let a = CopyStructure{ copy_id:"A".into(), chrom:"chr1".into(), strand:'+', exons:vec![(0,30)] };
    let mut b = a.clone(); b.copy_id="B".into(); b.chrom="chr5".into(); b.strand='-';
    let s = b"ACGTACGTACGTACGTACGTACGTACGTAC".to_vec();
    // build_families groups by strand BEFORE clustering, so a + and - copy never share a family.
    let fams = rustle::annotation_families::families_from_grouped(
        vec![(a, s.clone()), (b, s)], 0.5);
    assert_eq!(fams.len(), 0, "opposite-strand copies must not merge (build_family_graph bails on mixed strand)");
}
```

- [ ] **Step 2: Run — expect FAIL.**

- [ ] **Step 3: Implement `families_from_grouped` + `build_families`**

In `src/rustle/annotation_families.rs`:

```rust
/// Cluster within each strand group separately (same-strand is the only structural precondition,
/// per build_family_graph at family_graph.rs:493-497), then concat. Position-agnostic within strand.
pub fn families_from_grouped(copies_with_seq: Vec<(CopyStructure, Vec<u8>)>, threshold: f64) -> Vec<AnnotationFamily> {
    use std::collections::BTreeMap;
    let mut by_strand: BTreeMap<char, Vec<(CopyStructure, Vec<u8>)>> = BTreeMap::new();
    for cs in copies_with_seq { by_strand.entry(cs.0.strand).or_default().push(cs); }
    let mut out = Vec::new();
    for (_, group) in by_strand { out.extend(cluster_families(group, threshold)); }
    out
}

/// Top-level: load -G2, fetch exon sequences from genome, cluster (position-agnostic, same-strand).
pub fn build_families<P: AsRef<Path>>(
    g2_path: P, genome: &crate::genome::GenomeIndex, threshold: f64,
) -> Result<Vec<AnnotationFamily>> {
    let copies = load_copies(g2_path)?;
    let mut with_seq = Vec::new();
    for c in copies {
        if let Some(s) = copy_sequence(&c, genome) { with_seq.push((c, s)); }
    }
    Ok(families_from_grouped(with_seq, threshold))
}
```

- [ ] **Step 4: Run — expect PASS.**

- [ ] **Step 5: Commit**

```bash
git add src/rustle/annotation_families.rs tests/regression/g2_annotation_families.rs
git commit -m "feat(g2): build_families entry (load+seq+cluster, same-strand split)"
```

---

## Task 6: `vg_families.tsv` report (per-copy chroms via per_copy_spans)

**Files:**
- Create: `src/rustle/vg_eval.rs`
- Modify: `src/rustle/lib.rs` (`pub mod vg_eval;`)
- Test: `tests/regression/g2_annotation_families.rs`

- [ ] **Step 1: Failing test (report shows both chroms)**

```rust
#[test]
fn vg_families_report_shows_per_copy_chroms() {
    use rustle::annotation_families::{AnnotationFamily, CopyStructure};
    let fam = AnnotationFamily {
        family_id: "FAM0".into(),
        copies: vec![
            CopyStructure{copy_id:"A".into(),chrom:"chr1".into(),strand:'+',exons:vec![(0,30)]},
            CopyStructure{copy_id:"B".into(),chrom:"chr5".into(),strand:'+',exons:vec![(0,30)]},
        ],
        achieved_min_sim: 0.93, achieved_mean_sim: 0.93,
    };
    let tsv = rustle::vg_eval::render_vg_families(&[fam], 0.5);
    assert!(tsv.contains("chr1;chr5"), "copy_chroms must list both contigs; got:\n{}", tsv);
    assert!(tsv.contains("0.500"), "must record merge_threshold_T");
    assert!(tsv.contains("0.930"), "must record achieved similarity");
}
```

- [ ] **Step 2: Run — expect FAIL.**

- [ ] **Step 3: Implement the report renderer**

Create `src/rustle/vg_eval.rs`:

```rust
//! -G2 reports: vg_families.tsv (family merge audit) + vg_eval.tsv (head-to-head vs -G).
use crate::annotation_families::AnnotationFamily;

/// Render the family-merge audit. copy_chroms is read from each copy's own chrom
/// (NOT a family-graph node label, which collapses multi-chrom families to copies[0]).
pub fn render_vg_families(families: &[AnnotationFamily], threshold: f64) -> String {
    let mut out = String::from("family_id\tn_copies\tcopy_chroms\tmerge_threshold_T\tachieved_min_sim\tachieved_mean_sim\n");
    for f in families {
        let chroms: Vec<&str> = f.copies.iter().map(|c| c.chrom.as_str()).collect();
        out.push_str(&format!("{}\t{}\t{}\t{:.3}\t{:.3}\t{:.3}\n",
            f.family_id, f.copies.len(), chroms.join(";"),
            threshold, f.achieved_min_sim, f.achieved_mean_sim));
    }
    out
}
```

In `src/rustle/lib.rs`: `pub mod vg_eval; // -G2 reports`.

- [ ] **Step 4: Run — expect PASS.**

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_eval.rs src/rustle/lib.rs tests/regression/g2_annotation_families.rs
git commit -m "feat(g2): vg_families.tsv merge-audit report (per-copy chroms)"
```

---

## Task 7: Scorecard (`vg_eval.tsv`) — verdicts vs `-G`

**Files:**
- Modify: `src/rustle/vg_eval.rs`
- Test: `tests/regression/g2_annotation_families.rs`

- [ ] **Step 1: Failing test (verdict logic)**

```rust
#[test]
fn scorecard_verdicts() {
    use rustle::vg_eval::verdict;
    assert_eq!(verdict(true, false), "WIN");        // in_vg, not in stringtie
    assert_eq!(verdict(true, true), "tie");
    assert_eq!(verdict(false, false), "miss");
    assert_eq!(verdict(false, true), "regression"); // stringtie had it, vg lost it
}
```

- [ ] **Step 2: Run — expect FAIL.**

- [ ] **Step 3: Implement `verdict` + chain-match + `render_vg_eval`**

In `src/rustle/vg_eval.rs`:

```rust
/// Per-copy verdict. in_vg = rustle-VG emitted a matching chain; in_st = -G (StringTie) had it.
pub fn verdict(in_vg: bool, in_st: bool) -> &'static str {
    match (in_vg, in_st) {
        (true, false) => "WIN",
        (true, true) => "tie",
        (false, false) => "miss",
        (false, true) => "regression",
    }
}

/// True if `chain` (sorted introns) is in `set` within +/-2 bp per coordinate.
pub fn chain_in_set(chain: &[(u64,u64)], set: &std::collections::HashSet<Vec<(u64,u64)>>) -> bool {
    for cand in set {
        if cand.len() != chain.len() { continue; }
        if chain.iter().zip(cand).all(|(a,b)| a.0.abs_diff(b.0) <= 2 && a.1.abs_diff(b.1) <= 2) {
            return true;
        }
    }
    false
}

/// Render the head-to-head scorecard. Each entry = (copy_id, family_id, in_st, in_vg).
pub fn render_vg_eval(entries: &[(String, String, bool, bool)]) -> String {
    let mut out = String::from("copy_id\tfamily_id\tin_annotation\tin_stringtie\tin_vg\tverdict\n");
    let (mut win, mut tie, mut miss, mut reg) = (0,0,0,0);
    for (cid, fid, in_st, in_vg) in entries {
        let v = verdict(*in_vg, *in_st);
        match v { "WIN"=>win+=1, "tie"=>tie+=1, "miss"=>miss+=1, _=>reg+=1 }
        out.push_str(&format!("{}\t{}\tyes\t{}\t{}\t{}\n", cid, fid, in_st, in_vg, v));
    }
    out.push_str(&format!("SUMMARY\tmode=guided\tWIN={}\ttie={}\tmiss={}\tregression={}\n", win, tie, miss, reg));
    out
}
```

- [ ] **Step 4: Run — expect PASS.**

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_eval.rs tests/regression/g2_annotation_families.rs
git commit -m "feat(g2): scorecard verdicts + chain-match + vg_eval.tsv renderer"
```

---

## Task 8: Wire `-G2` into VG family discovery + emit both reports

**Files:**
- Modify: `src/rustle/vg.rs` (the family-discovery entry — locate the call site of `discover_family_groups`, ~line 269 caller)
- Modify: `src/rustle/lib.rs` or the run driver where VG output GTF is written, to call the report writers.

**Verify first:** read the `FamilyGroup` struct and how `discover_family_groups` results feed `build_family_graph`. The switch must produce the SAME `FamilyGroup` shape from `annotation_families` (each `AnnotationFamily.copies` → the bundles/loci that overlap those copies). Read `vg.rs:269-417` and the `FamilyGroup` definition before implementing; do not invent fields.

- [ ] **Step 1: Failing integration test (RABL2-style or synthetic 2-copy)**

Add `tests/regression/g2_end_to_end.rs` that runs the binary on a small synthetic 2-copy BAM+GTF+FASTA fixture with `--vg --guide stringtie.gtf --guide2 annotation.gtf --genome-fasta g.fa` and asserts `out.vg_families.tsv` and `out.vg_eval.tsv` exist and the family report lists 2 copies. (Reuse an existing VG integration fixture under `tests/regression/` as the template; e.g. the smallest `vg_family_*` fixture.)

- [ ] **Step 2: Run — expect FAIL** (reports not written).

- [ ] **Step 3: Implement the discovery switch + report emission**

At the family-discovery entry in `vg.rs`, gate on `config.guide2_path`:

```rust
let families = if let Some(g2) = config.guide2_path.as_deref() {
    let genome = genome.expect("--guide2 requires --genome-fasta (validated in main)");
    let ann = crate::annotation_families::build_families(g2, genome, config.family_exon_similarity)?;
    // Map each AnnotationFamily to a FamilyGroup over the bundles its copies overlap.
    // (Use the SAME FamilyGroup constructor discover_family_groups uses — verified in Step 0.)
    annotation_families_to_groups(&ann, bundles)
} else {
    discover_family_groups(bundles, /* existing args */)
};
```

After the VG output GTF is written, emit the reports:

```rust
if let Some(g2) = config.guide2_path.as_deref() {
    let genome = genome.expect("validated");
    let ann = crate::annotation_families::build_families(g2, genome, config.family_exon_similarity)?;
    std::fs::write(format!("{}.vg_families.tsv", out_base),
        crate::vg_eval::render_vg_families(&ann, config.family_exon_similarity))?;
    // Build scorecard entries: for each annotation copy, chain-match against -G (stringtie) and out.gtf.
    let entries = crate::vg_eval::build_eval_entries(&ann, &stringtie_chains, &vg_output_chains);
    std::fs::write(format!("{}.vg_eval.tsv", out_base),
        crate::vg_eval::render_vg_eval(&entries))?;
}
```

Add `build_eval_entries` to `vg_eval.rs` (loads `-G` GTF chains + the just-written `out.gtf` chains via `reference_gtf::parse_reference_gtf`, computes `chain_in_set` per copy). Add `annotation_families_to_groups` near the discovery code.

- [ ] **Step 4: Run — expect PASS.**

- [ ] **Step 5: Escape-hatch + suite**

Run: `RUSTLE_PRECISE=1 bash bench/mini3/check_precise.sh` (byte-matches) and `cargo test --release --lib` (green).

- [ ] **Step 6: Commit**

```bash
git add src/rustle/vg.rs src/rustle/vg_eval.rs tests/regression/g2_end_to_end.rs
git commit -m "feat(g2): source VG families from -G2 + emit vg_families/vg_eval reports"
```

---

## Task 9: Regression — read-based trans-chromosomal linker (closes coverage gap)

**Files:**
- Test: `tests/regression/g2_annotation_families.rs`

The audit found `discover_family_groups` (vg.rs:269) — the autonomous trans-chrom path RABL2 relies on — has **zero direct test coverage**. Add one.

- [ ] **Step 1: Failing test**

```rust
#[test]
fn read_based_linker_unions_trans_chromosomal_bundles() {
    // Two bundles on different chroms sharing >= min_shared_reads multimappers by read_name_hash.
    // Build the minimal Bundle fixtures (reuse the helper used by existing vg.rs unit tests at
    // ~vg.rs:7036-7340) with the same read_name_hash present in both bundles, on chr1 and chr5.
    // Assert discover_family_groups returns ONE FamilyGroup spanning both bundles.
    // (Fill in with the verified Bundle/FamilyGroup constructors.)
}
```

- [ ] **Step 2: Run — expect FAIL.**
- [ ] **Step 3: Implement the fixture using the verified constructors (read `vg.rs:7036-7340` for the existing bundle test helper).**
- [ ] **Step 4: Run — expect PASS.**
- [ ] **Step 5: Commit**

```bash
git add tests/regression/g2_annotation_families.rs
git commit -m "test(g2): regression for read-based trans-chromosomal family linker"
```

---

## Task 10: Invariant — `-G2` absent is byte-identical; misuse errors

**Files:**
- Test: `tests/regression/g2_end_to_end.rs`

- [ ] **Step 1: Test — `-G2` absent leaves `--vg` output unchanged**

Run the existing smallest `--vg` fixture WITHOUT `--guide2`, capture the output GTF, and assert it equals a stored golden (or equals a second no-`-G2` run byte-for-byte). Assert no `.vg_eval.tsv` / `.vg_families.tsv` are written.

- [ ] **Step 2: Test — `--guide2` without `--vg` exits non-zero** with the "requires --vg" message; `--guide2` without `--genome-fasta` exits non-zero.

- [ ] **Step 3: Run both — expect PASS.**

- [ ] **Step 4: Final escape-hatch + full suite**

Run: `RUSTLE_PRECISE=1 bash bench/mini3/check_precise.sh` and `cargo test --release` (green).

- [ ] **Step 5: Commit**

```bash
git add tests/regression/g2_end_to_end.rs
git commit -m "test(g2): -G2-absent byte-identical + misuse-error invariants"
```

---

## Notes for the implementer

- **Three "verify first" points** (Tasks 3, 8, 9) require reading an exact existing signature before writing the call: `GenomeIndex` subsequence accessor (`genome.rs`), `FamilyGroup` constructor + `discover_family_groups` caller (`vg.rs:269-417`), and the bundle test helper (`vg.rs:7036-7340`). These are deliberate — do not invent the signatures.
- **Iteration 2 (deferred, not in this plan):** `--g2-eval-only` (load `-G2` for the scorecard but do NOT feed family definition → fair de-novo benchmark, `mode=denovo`); clustering-threshold tuning; genome-wide scorecard aggregation.
- **Test oracle:** use **RABL2** (validated trans-chromosomal showcase). Do NOT use DAZ (retired false positive).
