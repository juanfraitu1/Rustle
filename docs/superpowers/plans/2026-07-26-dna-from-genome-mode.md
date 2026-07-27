# `--from-genome` Mode Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a read-free, annotation-free `--from-genome` front-end to `gw_family_catalog` that discovers segmental-duplication families from the CHM13v2.0 genome alone and feeds the exact same homology-grouping core the RNA `--homology-primary` path uses, so a DNA-vs-RNA comparison on the Soto benchmark isolates the substrate (splicing) as the cause of RNA's lower recall.

**Architecture:** One engine, two front-ends. A new `from_genome::genome_reps` discovers duplicated genomic loci by self-alignment (reusing `genome_projection`'s minimap2 primitive) and emits them as `DenovoTranscript` reps with genomic `seq` and an empty intron chain. A shared, extracted `denovo_pipeline::homology_blocks` (the `homology_edges_all_reps → gamma_quasi_clique_partition` core) groups reps into families for BOTH the RNA and DNA paths. The binary gains a `--from-genome <windows.bed>` flag that dispatches to the DNA front-end, then reuses the identical family/copies emit code.

**Tech Stack:** Rust (cargo, clap, anyhow, rayon), minimap2 (external, invoked via `std::process::Command`), samtools/tabix (external, for the deliverable scripts only). Python 3 for the scoring/deliverable script.

## Global Constraints

- The RNA code path must stay behavior-identical: the three existing tests `tests/gw_family_catalog_homology_primary.rs`, `tests/gw_family_catalog_regression.rs`, `tests/gw_family_catalog_same_chrom_supplement.rs` MUST still pass unchanged after every task.
- `--from-genome` is opt-in; with it absent, the binary behaves exactly as today (`--bam` remains the normal entry).
- DNA reps carry `introns: vec![]` and `seq` = full genomic locus sequence (introns included). This invariant is the scientific claim and is unit-tested.
- Family definition uses only genome sequence; the Soto BED is used only for scoring.
- Reference genome: T2T-CHM13v2.0 at `/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0.fa`. Benchmark BED: `bench/soto/80_fams.chr.bed`.
- minimap2 binary path from `RefineParams.minimap2` (env `RUSTLE_MINIMAP2`, default `minimap2`).
- Build/test on WSL2: run `cargo` in the foreground; do not background the compile.

---

### Task 1: Extract the shared grouping core `homology_blocks`

Extract the homology edges + γ-quasi-clique partition (the "same engine") into one function that both the RNA genome-wide catalog and the new DNA path call. Pure refactor — no behavior change.

**Files:**
- Modify: `src/rustle/vg_family/denovo_pipeline.rs` (add `homology_blocks`; call it from `detect_homology_catalog_genome_wide` at the current lines 2422-2425)
- Test: `src/rustle/vg_family/denovo_pipeline.rs` (`#[cfg(test)]` module) + the existing `tests/gw_family_catalog_homology_primary.rs` as the behavior guard

**Interfaces:**
- Consumes: `homology_edges_all_reps(&[DenovoTranscript], &RefineParams) -> Result<Vec<(usize,usize)>>` (existing, `pub(crate)`), `family_split::gamma_quasi_clique_partition(usize, &[(usize,usize,f64)], f64) -> Vec<Vec<usize>>` (existing, `pub`).
- Produces: `pub(crate) fn homology_blocks(reps: &[DenovoTranscript], refine: &RefineParams, gamma: f64) -> Result<Vec<Vec<usize>>>` — the rep-index blocks (γ-quasi-cliques of the E_r homology graph). Later tasks group reps with it.

- [ ] **Step 1: Write the failing test**

Add to the `#[cfg(test)] mod tests` in `denovo_pipeline.rs`:

```rust
#[test]
fn homology_blocks_groups_identical_reps_and_isolates_unrelated() {
    // three reps: two identical sequences (should share an E_r edge -> one block) + one unrelated.
    if std::process::Command::new("minimap2").arg("--version").output().is_err() { return; }
    let mk = |tid: &str, chrom: &str, start: u64, seq: &[u8]| DenovoTranscript {
        tid: tid.into(), chrom: chrom.into(), start, end: start + seq.len() as u64,
        n_reads: 5, strand: '+', introns: vec![], seq: seq.to_vec(), distinguishing_uniq: 0,
    };
    // a 400 bp "gene" duplicated at two loci, plus an unrelated 400 bp sequence.
    let a: Vec<u8> = b"ACGT".iter().cycle().take(400).copied().collect();
    let b: Vec<u8> = b"TGCA".iter().cycle().take(400).copied().collect();
    let reps = vec![
        mk("d1", "chr1", 1000, &a),
        mk("d2", "chr9", 5000, &a),   // identical to d1 -> same family
        mk("d3", "chr1", 9000, &b),   // unrelated -> its own block
    ];
    let refine = homology_refine_params(None, 2);
    let blocks = homology_blocks(&reps, &refine, 0.20).unwrap();
    // d1 and d2 land in the same block; d3 is alone.
    let block_of = |i: usize| blocks.iter().position(|bl| bl.contains(&i)).unwrap();
    assert_eq!(block_of(0), block_of(1), "identical reps must share a block");
    assert_ne!(block_of(0), block_of(2), "unrelated rep must be a separate block");
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --lib homology_blocks_groups_identical_reps -- --nocapture`
Expected: FAIL to compile with `cannot find function homology_blocks`.

- [ ] **Step 3: Add `homology_blocks` and call it from the RNA path**

Add the function (place it just above `detect_homology_catalog_genome_wide`):

```rust
/// The shared family-grouping ENGINE: E_r homology edges over rep sequences -> γ-quasi-clique blocks
/// (rep-index groups). Both the RNA genome-wide homology catalog and the DNA `--from-genome` path call
/// this, so "same engine, two substrates" is literally the same function. No read/annotation dependency:
/// it consumes only `rep.seq`.
pub(crate) fn homology_blocks(
    reps: &[DenovoTranscript],
    refine: &RefineParams,
    gamma: f64,
) -> Result<Vec<Vec<usize>>> {
    let edges2 = homology_edges_all_reps(reps, refine)?;
    let edges3: Vec<(usize, usize, f64)> = edges2.iter().map(|&(a, b)| (a, b, 1.0)).collect();
    Ok(crate::vg_family::family_split::gamma_quasi_clique_partition(reps.len(), &edges3, gamma))
}
```

Then in `detect_homology_catalog_genome_wide`, replace the three inline lines (currently ~2423-2425):

```rust
    let edges2 = homology_edges_all_reps(&reps, refine)?;
    let edges3: Vec<(usize, usize, f64)> = edges2.iter().map(|&(a, b)| (a, b, 1.0)).collect();
    let blocks = crate::vg_family::family_split::gamma_quasi_clique_partition(reps.len(), &edges3, gamma);
```

with:

```rust
    let blocks = homology_blocks(&reps, refine, gamma)?;
```

- [ ] **Step 4: Run the new test and the RNA behavior guard**

Run: `cargo test --lib homology_blocks_groups_identical_reps -- --nocapture`
Expected: PASS.
Run: `cargo test --test gw_family_catalog_homology_primary --test gw_family_catalog_regression --test gw_family_catalog_same_chrom_supplement`
Expected: PASS (RNA path unchanged — the refactor is byte-identical logic).

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/denovo_pipeline.rs
git commit -m "refactor: extract shared homology_blocks grouping core (RNA path byte-identical)"
```

---

### Task 2: `genome_reps` — the DNA front-end (SD detector + rep construction)

New module that discovers duplicated genomic loci in given windows and returns them as reps. Read-free, annotation-free.

**Files:**
- Create: `src/rustle/vg_family/from_genome.rs`
- Modify: `src/rustle/vg_family/mod.rs` (add `pub mod from_genome;`)
- Test: `src/rustle/vg_family/from_genome.rs` (`#[cfg(test)]` module) + a tiny in-repo FASTA fixture `tests/fixtures/two_copy.fa`

**Interfaces:**
- Consumes: `genome_projection::project_families_batch(consensuses: &[(String, Vec<u8>)], fasta: &str, known: &HashMap<String, Vec<(String,u64,u64)>>, min_identity: f64, min_cov: f64, minimap2: &str, threads: usize) -> Result<HashMap<String, Vec<CopyLocus>>>` (existing — returns, per query, the genome loci where it recurs; `pub struct CopyLocus { pub chrom: String, pub start: u64, pub end: u64, pub identity: f64, pub cov: f64 }` — VERIFIED). `crate::genome::GenomeIndex::from_fasta_contigs(fasta: &str, &HashSet<String>) -> Result<GenomeIndex>` and `GenomeIndex::fetch_sequence(&self, chrom: &str, start: u64, end: u64) -> Option<Vec<u8>>` (existing, VERIFIED — note: `GenomeIndex` upper-cases sequence, which is fine for homology alignment) for sequence extraction.
- Produces: `pub fn genome_reps(fasta_path: &str, windows: &[(String, u64, u64)], p: &GenomeRepParams) -> Result<Vec<DenovoTranscript>>`; `pub struct GenomeRepParams { pub min_identity: f64, pub min_block: u64, pub max_locus_span: u64, pub minimap2: String, pub threads: usize }` with `Default` (`0.90, 1000, 3_000_000, "minimap2", 4`) and `from_env` (`RUSTLE_GENOME_MIN_IDENTITY`, `RUSTLE_GENOME_MIN_BLOCK`, `RUSTLE_MINIMAP2`).

- [ ] **Step 1: Write the failing unit test (rep invariant + discovery)**

Create `tests/fixtures/two_copy.fa` — a reference with the SAME 600 bp sequence at two loci on one contig, plus filler so the two copies are disjoint:

```
>chr_test
<300 bp of ACGT filler><600 bp motif M><300 bp filler><600 bp motif M><300 bp filler>
```

(Generate M as a fixed pseudo-random-but-static 600 bp string; write the exact bases into the fixture so the test is deterministic. Total contig ~2.4 kb. Add `tests/fixtures/two_copy.fa.fai` via `samtools faidx` and commit both.)

Add to `from_genome.rs` tests:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn genome_reps_finds_two_copies_with_genomic_seq_and_no_introns() {
        if std::process::Command::new("minimap2").arg("--version").output().is_err() { return; }
        let fa = "tests/fixtures/two_copy.fa";
        // one window spanning the whole contig; genome_reps self-aligns and finds the duplicated motif twice.
        let windows = vec![("chr_test".to_string(), 0u64, 2400u64)];
        let p = GenomeRepParams { min_identity: 0.95, min_block: 400, ..Default::default() };
        let reps = genome_reps(fa, &windows, &p).unwrap();
        assert!(reps.len() >= 2, "expected >=2 duplicated loci, got {}", reps.len());
        for r in &reps {
            assert!(r.introns.is_empty(), "DNA rep must have an empty intron chain");
            assert_eq!(r.seq.len() as u64, r.end - r.start, "DNA rep seq must be the full genomic span");
        }
    }
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --lib genome_reps_finds_two_copies -- --nocapture`
Expected: FAIL to compile (`from_genome` / `genome_reps` not found).

- [ ] **Step 3: Implement `genome_reps`**

```rust
//! DNA front-end: discover duplicated genomic loci by self-alignment and emit them as reps for the
//! shared homology-grouping core. Read-free and annotation-free — the genome-only counterpart of the
//! RNA read front-end. The reps differ from RNA reps in exactly one way: genomic `seq` (introns included)
//! and an empty intron chain.
use std::collections::{HashMap, HashSet};
use anyhow::Result;
use crate::vg_family::family_detect::DenovoTranscript;
use crate::genome::GenomeIndex;
use crate::vg_family::genome_projection::project_families_batch;

pub struct GenomeRepParams {
    pub min_identity: f64,
    pub min_block: u64,
    pub max_locus_span: u64,
    pub minimap2: String,
    pub threads: usize,
}
impl Default for GenomeRepParams {
    fn default() -> Self {
        Self { min_identity: 0.90, min_block: 1000, max_locus_span: 3_000_000,
                minimap2: "minimap2".into(), threads: 4 }
    }
}
impl GenomeRepParams {
    pub fn from_env() -> Self {
        let mut p = Self::default();
        if let Ok(v) = std::env::var("RUSTLE_GENOME_MIN_IDENTITY") { if let Ok(x) = v.parse() { p.min_identity = x; } }
        if let Ok(v) = std::env::var("RUSTLE_GENOME_MIN_BLOCK") { if let Ok(x) = v.parse() { p.min_block = x; } }
        if let Ok(v) = std::env::var("RUSTLE_MINIMAP2") { p.minimap2 = v; }
        p
    }
}

pub fn genome_reps(
    fasta_path: &str,
    windows: &[(String, u64, u64)],
    p: &GenomeRepParams,
) -> Result<Vec<DenovoTranscript>> {
    let contigs: HashSet<String> = windows.iter().map(|(c, _, _)| c.clone()).collect();
    let genome = GenomeIndex::from_fasta_contigs(fasta_path, &contigs)?;

    // (1) SD detector: each window's sequence is a query; project_families_batch returns every genome
    // locus it recurs at (identity >= min_identity over an aligned block). One batched minimap2 pass.
    let consensuses: Vec<(String, Vec<u8>)> = windows.iter().enumerate().filter_map(|(i, (c, s, e))| {
        genome.fetch_sequence(c, *s, *e).map(|seq| (format!("w{i}"), seq))
    }).collect();
    let known: HashMap<String, Vec<(String, u64, u64)>> = HashMap::new(); // no loci to exclude; keep all hits
    let cov = 0.0_f64; // block length is gated below by min_block, not by fractional coverage of the window
    let hits = project_families_batch(&consensuses, fasta_path, &known, p.min_identity, cov,
                                      &p.minimap2, p.threads)?;

    // (2) rep construction: collect all hit loci (across windows), keep blocks >= min_block and spans
    // <= max_locus_span, merge overlapping loci (single-linkage by genomic overlap), emit one rep each.
    let mut loci: Vec<(String, u64, u64)> = Vec::new();
    for (_q, hs) in hits {
        for h in hs {
            if h.end.saturating_sub(h.start) >= p.min_block
                && h.end.saturating_sub(h.start) <= p.max_locus_span {
                loci.push((h.chrom, h.start, h.end));
            }
        }
    }
    loci.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    let merged = merge_overlapping(&loci);

    let mut reps = Vec::new();
    for (chrom, start, end) in merged {
        if let Some(seq) = genome.fetch_sequence(&chrom, start, end) {
            reps.push(DenovoTranscript {
                tid: format!("DN_{chrom}_{start}_1"),
                chrom, start, end, n_reads: 1, strand: '+',
                introns: vec![], seq, distinguishing_uniq: 0,
            });
        }
    }
    Ok(reps)
}

/// Single-linkage merge of overlapping genomic intervals (sorted input).
fn merge_overlapping(loci: &[(String, u64, u64)]) -> Vec<(String, u64, u64)> {
    let mut out: Vec<(String, u64, u64)> = Vec::new();
    for (c, s, e) in loci.iter().cloned() {
        match out.last_mut() {
            Some((pc, _ps, pe)) if *pc == c && s <= *pe => { *pe = (*pe).max(e); }
            _ => out.push((c, s, e)),
        }
    }
    out
}
```

Add `pub mod from_genome;` to `src/rustle/vg_family/mod.rs`.

- [ ] **Step 4: Run tests**

Run: `cargo test --lib genome_reps_finds_two_copies -- --nocapture`
Expected: PASS (>=2 reps, each with empty introns and genomic-length seq).
Run: `cargo build --bin gw_family_catalog`
Expected: builds clean.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/from_genome.rs src/rustle/vg_family/mod.rs tests/fixtures/two_copy.fa tests/fixtures/two_copy.fa.fai
git commit -m "feat: genome_reps DNA front-end (self-align -> genomic reps, empty intron chain)"
```

---

### Task 3: Wire `--from-genome` into the binary

Add the flag, make `--bam` optional, dispatch to the DNA front-end + shared grouping, reuse the emit path.

**Files:**
- Modify: `src/bin/gw_family_catalog.rs`
- Modify: `src/rustle/vg_family/denovo_pipeline.rs` (add `families_from_reps` convenience wrapper)
- Test: `tests/gw_family_catalog_from_genome.rs` (new integration test)

**Interfaces:**
- Consumes: `homology_blocks` (Task 1), `genome_reps` + `GenomeRepParams` (Task 2), existing `distinct_locus_reps(Vec<DenovoTranscript>, usize) -> Vec<DenovoTranscript>` (private in denovo_pipeline.rs).
- Produces: `pub fn families_from_reps(reps: Vec<DenovoTranscript>, refine: &RefineParams, gamma: f64, min_copies: usize, min_reads: usize) -> Result<Vec<Vec<DenovoTranscript>>>` (denovo_pipeline.rs) — blocks via `homology_blocks`, then `distinct_locus_reps` + `>= min_copies` filter. The DNA binary path calls this and feeds the result to the existing families/copies emit code.

- [ ] **Step 1: Write the failing integration test**

Create `tests/gw_family_catalog_from_genome.rs`:

```rust
// Runs the binary in --from-genome mode on the two-copy fixture and asserts a >=2-copy family is emitted.
use std::process::Command;

#[test]
fn from_genome_emits_a_two_copy_family() {
    if Command::new("minimap2").arg("--version").output().is_err() { return; }
    let dir = std::env::temp_dir().join("rustle_fromgenome_test");
    let _ = std::fs::create_dir_all(&dir);
    let win = dir.join("win.bed");
    std::fs::write(&win, "chr_test\t0\t2400\n").unwrap();
    let out = dir.join("dna");
    let status = Command::new(env!("CARGO_BIN_EXE_gw_family_catalog"))
        .args(["--from-genome", win.to_str().unwrap(),
               "--fasta", "tests/fixtures/two_copy.fa",
               "--min-identity", "0.95",
               "--out", out.to_str().unwrap()])
        .env("RUSTLE_GENOME_MIN_BLOCK", "400")
        .status().unwrap();
    assert!(status.success());
    let fams = std::fs::read_to_string(format!("{}.families.tsv", out.to_str().unwrap())).unwrap();
    // header + >=1 family row; the family has n_copies >= 2 (column 2).
    let rows: Vec<&str> = fams.lines().skip(1).collect();
    assert!(!rows.is_empty(), "expected >=1 family");
    assert!(rows.iter().any(|r| r.split('\t').nth(1).map(|n| n.parse::<usize>().unwrap_or(0) >= 2).unwrap_or(false)),
            "expected a family with >=2 copies");
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --test gw_family_catalog_from_genome`
Expected: FAIL — `--from-genome` is an unknown argument (clap error), non-zero exit.

- [ ] **Step 3: Add `families_from_reps`, then wire the flag**

In `denovo_pipeline.rs` add (below `homology_blocks`):

```rust
/// Group a rep set into families with the shared engine (`homology_blocks`) and keep those spanning
/// >= `min_copies` spatially-distinct loci. Used by the DNA `--from-genome` path; the RNA path keeps its
/// collapse-aware loop but shares `homology_blocks` for the grouping itself.
pub fn families_from_reps(
    reps: Vec<DenovoTranscript>,
    refine: &RefineParams,
    gamma: f64,
    min_copies: usize,
    min_reads: usize,
) -> Result<Vec<Vec<DenovoTranscript>>> {
    let blocks = homology_blocks(&reps, refine, gamma)?;
    let mut out = Vec::new();
    for block in blocks {
        let copies: Vec<DenovoTranscript> = block.iter().map(|&i| reps[i].clone()).collect();
        let loci = distinct_locus_reps(copies, min_reads);
        if block.len() >= min_copies && loci.len() >= min_copies {
            out.push(loci);
        }
    }
    Ok(out)
}
```

In `src/bin/gw_family_catalog.rs`:

1. Change the `bam` field to optional and add the flag:

```rust
    /// Coordinate-sorted BAM (a `.bai` next to it enables the fast per-contig region query). Required
    /// unless --from-genome is given.
    #[arg(long)]
    bam: Option<String>,
    /// GENOME-ONLY mode: discover SD families from the genome alone (no reads, no annotation). Argument is a
    /// BED of search windows (e.g. bench/soto/80_fams.chr.bed). Reps are duplicated genomic loci found by
    /// self-alignment; grouping uses the SAME homology_blocks core as --homology-primary. Mutually exclusive
    /// with --bam.
    #[arg(long)]
    from_genome: Option<String>,
```

2. At the top of `main`, after `Args::parse()`, add the DNA branch (before the single-copy-baseline block). It builds reps, groups with `families_from_reps`, and jumps straight to the existing emit code by binding `fams`:

```rust
    // --from-genome: read-free/annotation-free DNA family catalog.
    if let Some(win_bed) = args.from_genome.as_ref() {
        use rustle::vg_family::from_genome::{genome_reps, GenomeRepParams};
        let windows = read_windows_bed(win_bed)?;           // helper below
        let mut gp = GenomeRepParams::from_env();
        gp.threads = args.threads;
        gp.minimap2 = refine_params.minimap2.clone();
        let reps = genome_reps(&args.fasta, &windows, &gp)?;
        eprintln!("[gw-catalog-genome] {} windows -> {} duplicated-locus reps", windows.len(), reps.len());
        let dna_fams = rustle::vg_family::denovo_pipeline::families_from_reps(
            reps, &refine_params, 0.20, args.min_copies, 0,
        )?;
        emit_catalog(&args.out, dna_fams, &refine_params)?;  // extracted emit (below)
        return Ok(());
    }
```

3. Extract the families/copies/`copies.fa` writing loop (current lines ~235-321) into `fn emit_catalog(out: &str, fams: Vec<Vec<DenovoTranscript>>, refine_params: &RefineParams) -> Result<()>` and call it from both the existing path and the DNA branch. (Pure move — same code, so existing output is unchanged.)

4. Add helpers near the top of the file:

```rust
fn read_windows_bed(path: &str) -> Result<Vec<(String, u64, u64)>> {
    let mut out = Vec::new();
    for line in std::fs::read_to_string(path)?.lines() {
        if line.is_empty() || line.starts_with('#') { continue; }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() >= 3 { out.push((f[0].to_string(), f[1].parse()?, f[2].parse()?)); }
    }
    Ok(out)
}
```

5. Guard the existing `--bam`-consuming paths: replace `&args.bam` with an unwrap that errors clearly when neither `--bam` nor `--from-genome` is set. Add right after the DNA branch:

```rust
    let bam = args.bam.as_deref().ok_or_else(|| anyhow::anyhow!("--bam is required unless --from-genome is given"))?;
```

and use `bam` in place of `&args.bam` in the calls below (single-copy baseline, `detect_homology_catalog_genome_wide`, `detect_conflict_catalog_genome_wide*`, and the `--project-all-families` `reads_in_region`).

- [ ] **Step 4: Run tests**

Run: `cargo test --test gw_family_catalog_from_genome`
Expected: PASS (a >=2-copy family emitted).
Run: `cargo test --test gw_family_catalog_homology_primary --test gw_family_catalog_regression --test gw_family_catalog_same_chrom_supplement`
Expected: PASS (RNA path unchanged — `emit_catalog` is a pure code move; `--bam` still works).

- [ ] **Step 5: Commit**

```bash
git add src/bin/gw_family_catalog.rs src/rustle/vg_family/denovo_pipeline.rs tests/gw_family_catalog_from_genome.rs
git commit -m "feat: --from-genome mode wired into gw_family_catalog (shared grouping engine)"
```

---

### Task 4: Benchmark run + DNA-vs-RNA deliverable

Run DNA mode on the 83-family benchmark, score against Soto with the same rule as the RNA side, and write the comparison + honest framing.

**Files:**
- Create: `bench/soto/score_from_genome.py` (scorer: DNA copies.tsv vs `80_fams.chr.bed`, member-recovery rule identical to the RNA scorer)
- Create: `bench/soto/dna_vs_rna_mode.md` (the deliverable table + framing)
- Create: `bench/soto/run_from_genome.sh` (the end-to-end command, resumable)

**Interfaces:**
- Consumes: the built `target/release/gw_family_catalog` with `--from-genome`; `bench/soto/80_fams.chr.bed`; `bench/soto/member_attribution.final.tsv` (RNA result: 313/362 counting top-up, or 306 strict — cite both).
- Produces: `dna_mode.families.tsv` / `dna_mode.copies.tsv` in `winloci_scratch`, a member-recall number, and `bench/soto/dna_vs_rna_mode.md`.

- [ ] **Step 1: Write `score_from_genome.py`**

```python
#!/usr/bin/env python3
"""Score a --from-genome copies.tsv against the Soto benchmark, same rule as the RNA scorer:
a Soto member is RECOVERED if a >=2-copy DNA family has a copy locus overlapping it."""
import sys, csv
from collections import defaultdict
BED = "bench/soto/80_fams.chr.bed"
copies_tsv = sys.argv[1]  # dna_mode.copies.tsv
# family_id -> list of (chrom,start,end)
fam = defaultdict(list)
for i, ln in enumerate(open(copies_tsv)):
    if i == 0: continue
    f = ln.rstrip("\n").split("\t")  # family_id copy_idx tid chrom start end n_exon strand n_reads
    fam[f[0]].append((f[3], int(f[4]), int(f[5])))
# only >=2-copy families count
copyloci = [(c, s, e) for fid, locs in fam.items() if len(locs) >= 2 for (c, s, e) in locs]
members = []
for ln in open(BED):
    c, s, e, name, *_ = ln.rstrip("\n").split("\t")
    members.append((name.split("|")[1], name.split("|")[0], c, int(s), int(e)))
def hit(c, s, e):
    return any(cc == c and s < ee and e > ss for (cc, ss, ee) in copyloci)
rec = sum(1 for (_f, _g, c, s, e) in members if hit(c, s, e))
print(f"DNA member recovery: {rec}/{len(members)} = {100*rec/len(members):.1f}%")
```

- [ ] **Step 2: Write and run `run_from_genome.sh`**

```bash
#!/usr/bin/env bash
set -euo pipefail
BIN=target/release/gw_family_catalog
FA=/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0.fa
BED=bench/soto/80_fams.chr.bed
OUT=/home/juanfra/winloci_scratch/dna_mode
cargo build --release --bin gw_family_catalog
"$BIN" --from-genome "$BED" --fasta "$FA" --min-identity 0.90 --out "$OUT"
python3 bench/soto/score_from_genome.py "$OUT.copies.tsv"
```

Run: `bash bench/soto/run_from_genome.sh`
Expected: prints a DNA member-recovery percentage (target ~90-98%).

- [ ] **Step 3: Write the deliverable `dna_vs_rna_mode.md`**

Fill in the measured DNA number and the RNA number (306 strict / 313 with top-up = 85-86.5%). Include:
- The comparison table (DNA vs RNA member recovery, same engine).
- One paragraph: "same `homology_blocks` grouping; the DNA front-end feeds genomic loci (copies pre-separated by the assembly), the RNA front-end feeds spliced reads (copies must be deconvolved). The gap is the substrate."
- The caveats verbatim from the spec's "What this proves / does NOT" section.

- [ ] **Step 4: Commit**

```bash
git add bench/soto/score_from_genome.py bench/soto/run_from_genome.sh bench/soto/dna_vs_rna_mode.md
git commit -m "bench: --from-genome Soto benchmark run + DNA-vs-RNA comparison deliverable"
```

---

## Self-Review

**Spec coverage:** goal (Task 3 flag + Task 4 run), one-engine-two-front-ends (Task 1 `homology_blocks` shared; Task 2 DNA front-end), DNA-rep invariant (Task 2 test), scoring vs Soto (Task 4), byte-identity guard (Tasks 1&3 run the three existing tests), two identity thresholds (Task 2 `GenomeRepParams.min_identity` for discovery; `--min-identity` for grouping via `refine_params`), skip read-based steps (Task 3 DNA branch returns before the read paths), honest framing (Task 4 deliverable). All spec sections map to a task.

**Placeholder scan:** every code step has real code; the fixture in Task 2 says to write explicit bases (no "TODO"); the deliverable numbers are filled from the Task 4 run, not left as TBD.

**Type consistency:** `homology_blocks(&[DenovoTranscript], &RefineParams, f64) -> Result<Vec<Vec<usize>>>` used identically in Tasks 1/3; `families_from_reps(Vec<DenovoTranscript>, &RefineParams, f64, usize, usize)` and `genome_reps(&str, &[(String,u64,u64)], &GenomeRepParams)` consistent across Tasks 2/3; `GenomeRepParams` fields match between definition (Task 2) and use (Task 3). `CopyLocus` fields (`chrom,start,end,identity,cov`) match the `project_families_batch` return used in Task 2 — verify the exact field names against `genome_projection.rs` at implementation time and adjust the `h.chrom/h.start/h.end` accessors if they differ.
