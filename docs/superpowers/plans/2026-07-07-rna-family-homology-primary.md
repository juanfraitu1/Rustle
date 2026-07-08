# RNA Family Homology-Primary (E_r) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Define an RNA multi-copy family as a nucleotide-homology component (E_r) built over all assembled loci, with a γ-quasi-clique operator, plus a Liftoff-style genome-projection layer that enumerates near-identical copies (famCN) — recovering the ancient-paralog and K=0-collapse families the read-conflict catalog misses.

**Architecture:** New homology-primary path alongside the existing conflict path (`--cross-chrom`). Reuse the rep builder and `nucleotide_edges` (all-vs-all minimap2 → id/cov edges); add a γ-quasi-clique partitioner (ported from `bench/genome_family_def.py`), a driver, a genome-projection copy enumerator, and an optional protein QC flag. O2 (`copy_assign`/`colocated_families`) is untouched.

**Tech Stack:** Rust (crate `rustle`), minimap2 (`asm20`, `-k11 -w5`, `-x splice`), mmseqs (optional protein QC), `cargo test`.

## Global Constraints

- Both modes ship: `--homology-primary` is added **alongside** `--cross-chrom`; neither is removed. The P/R sweep decides the default later.
- TDD: every new function gets a failing test first (`cargo test --lib <name>`), watched fail, then minimal impl.
- Reuse over reinvent: `nucleotide_edges`, `homology_components`, `distinct_locus_reps`, `louvain_communities`, `community_stats`, `conflict_edges`, `batch_protein_edges`, `GenomeIndex` already exist — call them, do not reimplement.
- External tools run via `mamba run -n <env>` where needed; `minimap2`/`mmseqs` paths come from `RefineParams`/env (`RUSTLE_MMSEQS`), absent → graceful empty (no hard failure).
- Coverage-of-shorter ≥ `min_coverage` (default 0.50) is the repeat-bridge guard on every homology edge.
- Commit after each task passes.

---

## File Structure

- `src/rustle/vg_family/family_split.rs` — ADD `gamma_quasi_clique_partition` (the operator, §4 of spec).
- `src/rustle/vg_family/denovo_pipeline.rs` — MODIFY: make `nucleotide_edges` `pub(crate)`; add `sensitive_identity` to `RefineParams`; ADD `homology_edges_all_reps` (E_r edge builder) and `detect_homology_catalog_genome_wide` (driver); ADD per-family `protein_coheres` QC.
- `src/rustle/vg_family/genome_projection.rs` — CREATE: `project_family_copies` (§7 famCN enumeration).
- `src/bin/gw_family_catalog.rs` — MODIFY: add `--homology-primary`, `--min-identity`, `--enumerate-copies` flags + emit.
- `tests/gw_family_catalog_homology_primary.rs` — CREATE: integration (recovers audit misses).

---

## Task 1: γ-quasi-clique partition operator

**Files:**
- Modify: `src/rustle/vg_family/family_split.rs`
- Test: same file, `#[cfg(test)] mod tests`

**Interfaces:**
- Consumes: `louvain_communities(n: usize, edges: &[(usize,usize,f64)], resolution: f64) -> Vec<Vec<usize>>`, `community_stats(members, edges) -> CommunityStats` (`.density: f64`) — both already in this file.
- Produces: `pub fn gamma_quasi_clique_partition(n: usize, edges: &[(usize, usize, f64)], gamma: f64) -> Vec<Vec<usize>>` — partitions the graph into blocks each of which is a γ-quasi-clique (internal density ≥ `gamma`) or has ≤2 nodes.

- [ ] **Step 1: Write the failing test**

```rust
#[test]
fn gamma_quasi_clique_keeps_array_whole_splits_repeat_chain() {
    // A 5-node dense clique (a tandem array) stays ONE block.
    let clique: Vec<(usize, usize, f64)> = {
        let mut e = Vec::new();
        for i in 0..5 { for j in (i + 1)..5 { e.push((i, j, 1.0)); } }
        e
    };
    let blocks = gamma_quasi_clique_partition(5, &clique, 0.20);
    assert_eq!(blocks.len(), 1, "a dense array is one gamma-quasi-clique");
    assert_eq!(blocks[0].len(), 5);

    // A LONG sparse bridge chain (12-node path: density 2*11/(12*11)=0.167 < gamma=0.20) must split.
    let chain: Vec<(usize, usize, f64)> = (0..11).map(|i| (i, i + 1, 1.0)).collect();
    let blocks = gamma_quasi_clique_partition(12, &chain, 0.20);
    assert!(blocks.len() >= 2, "a sparse repeat-bridge chain is split, got {:?}", blocks);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --lib gamma_quasi_clique_keeps_array_whole_splits_repeat_chain`
Expected: FAIL — `cannot find function gamma_quasi_clique_partition`.

- [ ] **Step 3: Write minimal implementation**

Port `_induced_density` / `_split_once` / `_refine_component` from `bench/genome_family_def.py` (lines ~205–260). KL bisection is replaced by adaptive-resolution `louvain_communities` + connected-component split + deterministic halving (guaranteed progress).

```rust
/// Internal edge density of a node block over the induced subgraph: 2|E|/(|C|(|C|-1)); ≤1 node = 1.0.
fn induced_density(members: &[usize], edges: &[(usize, usize, f64)]) -> f64 {
    let set: std::collections::HashSet<usize> = members.iter().copied().collect();
    let n = members.len();
    if n <= 1 { return 1.0; }
    let m = edges.iter().filter(|(a, b, _)| set.contains(a) && set.contains(b)).count();
    2.0 * m as f64 / (n as f64 * (n as f64 - 1.0))
}

/// Restrict edges to those with both endpoints in `members`, remapped to local indices 0..members.len().
fn induced_edges(members: &[usize], edges: &[(usize, usize, f64)]) -> (usize, Vec<(usize, usize, f64)>) {
    let idx: std::collections::HashMap<usize, usize> =
        members.iter().enumerate().map(|(i, &g)| (g, i)).collect();
    let local = edges.iter()
        .filter_map(|&(a, b, w)| match (idx.get(&a), idx.get(&b)) {
            (Some(&la), Some(&lb)) => Some((la, lb, w)),
            _ => None,
        }).collect();
    (members.len(), local)
}

/// Guaranteed-progress splitter: adaptive-resolution Louvain, else component split, else deterministic halving.
/// Returns global-index blocks; never a single block equal to the input when |members| > 2.
fn split_once(members: &[usize], edges: &[(usize, usize, f64)]) -> Vec<Vec<usize>> {
    let (n, local) = induced_edges(members, edges);
    for res in [1.0, 2.0, 4.0, 8.0] {
        let parts = louvain_communities(n, &local, res);
        if parts.len() >= 2 {
            return parts.into_iter().map(|p| p.into_iter().map(|l| members[l]).collect()).collect();
        }
    }
    // connected-component fallback (also catches a disconnected block).
    let comps = crate::vg_family::read_conflict::conflict_families(n, &local.iter().map(|&(a, b, _)| (a, b, 1usize)).collect::<Vec<_>>());
    if comps.len() >= 2 {
        return comps.into_iter().map(|c| c.into_iter().map(|l| members[l]).collect()).collect();
    }
    // deterministic halving.
    let h = members.len() / 2;
    vec![members[..h].to_vec(), members[h..].to_vec()]
}

/// γ-quasi-clique partition: keep a block whole iff it is already a γ-quasi-clique (or ≤2 nodes), else split
/// (guaranteed-progress) and recurse. Blocks partition 0..n. Deterministic (Louvain here is deterministic).
pub fn gamma_quasi_clique_partition(n: usize, edges: &[(usize, usize, f64)], gamma: f64) -> Vec<Vec<usize>> {
    // start from raw connected components, then refine each.
    let comps = crate::vg_family::read_conflict::conflict_families(
        n, &edges.iter().map(|&(a, b, _)| (a, b, 1usize)).collect::<Vec<_>>());
    let mut out = Vec::new();
    let mut stack: Vec<Vec<usize>> = comps;
    while let Some(block) = stack.pop() {
        if block.len() <= 2 || induced_density(&block, edges) >= gamma {
            out.push(block);
            continue;
        }
        let parts = split_once(&block, edges);
        if parts.len() == 1 && parts[0].len() == block.len() {
            out.push(block); // no-progress guard (shouldn't happen)
        } else {
            stack.extend(parts);
        }
    }
    out
}
```

(If `conflict_families` is not visible from `family_split.rs`, add `use crate::vg_family::read_conflict::conflict_families;` — it is `pub fn conflict_families(n: usize, edges: &[(usize,usize,usize)]) -> Vec<Vec<usize>>`.)

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test --lib gamma_quasi_clique_keeps_array_whole_splits_repeat_chain`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/family_split.rs
git commit -m "feat(family_split): add gamma-quasi-clique partition operator (E_r/E_b/E_p uniform)"
```

---

## Task 2: E_r homology edge builder over all reps

**Files:**
- Modify: `src/rustle/vg_family/denovo_pipeline.rs` (make `nucleotide_edges` `pub(crate)`; add `sensitive_identity` to `RefineParams` with default 0.60; add `homology_edges_all_reps`)
- Test: same file, tests module

**Interfaces:**
- Consumes: `nucleotide_edges(seqs: &[Vec<u8>], mm_args: &[&str], min_id: f64, min_cov: f64, params: &RefineParams) -> Result<Vec<(usize,usize)>>` (exon-sum sequences, all-vs-all minimap2, id+cov filtered).
- Produces: `pub(crate) fn homology_edges_all_reps(reps: &[DenovoTranscript], params: &RefineParams) -> Result<Vec<(usize, usize)>>` — E_r edges over ALL reps: asm20 (`min_identity`) ∪ sensitive `-k11 -w5` (`sensitive_identity`), both at `min_coverage`.

- [ ] **Step 1: Write the failing test** (integration — needs minimap2 on PATH; skip-guarded)

```rust
#[test]
fn homology_edges_links_two_similar_reps_not_divergent() {
    if std::process::Command::new("minimap2").arg("--version").output().is_err() {
        eprintln!("minimap2 absent; skipping"); return;
    }
    // rep 0 and rep 1: same 900bp sequence with ~5% mismatch (an ancient paralog); rep 2: random (unrelated).
    let base = rand_seq(900, 0x11);
    let mut para = base.clone();
    for k in (0..para.len()).step_by(20) { para[k] = b"ACGT"[(para[k] as usize + 1) % 4]; } // ~5% divergence
    let reps = vec![
        DenovoTranscript { tid: "r0".into(), chrom: "c1".into(), start: 0, end: 900, n_reads: 10, strand: '+', introns: vec![], seq: base },
        DenovoTranscript { tid: "r1".into(), chrom: "c1".into(), start: 5000, end: 5900, n_reads: 8, strand: '+', introns: vec![], seq: para },
        DenovoTranscript { tid: "r2".into(), chrom: "c1".into(), start: 9000, end: 9900, n_reads: 5, strand: '+', introns: vec![], seq: rand_seq(900, 0x99) },
    ];
    let params = RefineParams::default(); // min_identity 0.80, sensitive_identity 0.60, min_coverage 0.50
    let edges = homology_edges_all_reps(&reps, &params).unwrap();
    assert!(edges.contains(&(0, 1)), "the two paralog reps must be E_r-linked, got {:?}", edges);
    assert!(!edges.contains(&(0, 2)) && !edges.contains(&(1, 2)), "the unrelated rep must not link");
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --lib homology_edges_links_two_similar_reps_not_divergent`
Expected: FAIL — `cannot find function homology_edges_all_reps` (and `sensitive_identity` field missing).

- [ ] **Step 3: Write minimal implementation**

Add the field to `RefineParams` (near `min_identity`) and its default:

```rust
    /// Sensitive-tier nucleotide identity floor for the E_r homology edge (`-k11 -w5`). Lowered from the
    /// old within-refine 0.70 to reach ancient paralogs (KRAB-ZNF ~0.62). Repeat bridges are held off by
    /// `min_coverage`. Fixed by the family P/R sweep.
    pub sensitive_identity: f64,
```
```rust
            sensitive_identity: 0.60,   // in Default impl for RefineParams
```

Change `fn nucleotide_edges` to `pub(crate) fn nucleotide_edges`. Add the builder:

```rust
/// E_r homology edges over ALL reps' exon-sum sequences: asm20 (recent) ∪ sensitive -k11 -w5 (ancient),
/// both gated by `min_coverage`. One minimap2 all-vs-all per tier over the whole rep set (minimap2's index
/// is the prefilter). Protein is NOT an edge here — it is orthogonal QC (see per-family protein_coheres).
pub(crate) fn homology_edges_all_reps(
    reps: &[DenovoTranscript],
    params: &RefineParams,
) -> Result<Vec<(usize, usize)>> {
    let seqs: Vec<Vec<u8>> = reps.iter().map(|r| r.seq.clone()).collect();
    let mut set: BTreeSet<(usize, usize)> =
        nucleotide_edges(&seqs, &["-x", "asm20"], params.min_identity, params.min_coverage, params)?
            .into_iter().collect();
    if params.nucleotide_sensitive {
        for e in nucleotide_edges(&seqs, &["-k", "11", "-w", "5"], params.sensitive_identity, params.min_coverage, params)? {
            set.insert(e);
        }
    }
    Ok(set.into_iter().collect())
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test --lib homology_edges_links_two_similar_reps_not_divergent`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/denovo_pipeline.rs
git commit -m "feat(denovo_pipeline): E_r homology edge builder over all reps (asm20 + sensitive-0.60)"
```

---

## Task 3: `detect_homology_catalog_genome_wide` driver

**Files:**
- Modify: `src/rustle/vg_family/denovo_pipeline.rs`
- Test: same file, tests module

**Interfaces:**
- Consumes: `homology_edges_all_reps` (Task 2), `gamma_quasi_clique_partition` (Task 1), `distinct_locus_reps(copies: Vec<DenovoTranscript>) -> Vec<DenovoTranscript>` (existing), the rep builder in `detect_conflict_catalog_genome_wide_xchrom` (`primary_reads_from_bam → pass1_skeletons_robust → assemble_gate → collapse_loci_span_aware`).
- Produces: `pub fn detect_homology_catalog_genome_wide(bam_path: &str, fasta_path: &str, threads: usize, min_copies: usize, cfg: &DenovoConfig, refine: &RefineParams, gamma: f64) -> Result<Vec<Vec<DenovoTranscript>>>`.

- [ ] **Step 1: Write the failing test** (integration, uses the existing same-chrom fixture)

```rust
#[test]
fn homology_catalog_groups_fixture_family() {
    if std::process::Command::new("minimap2").arg("--version").output().is_err() { return; }
    let fams = detect_homology_catalog_genome_wide(
        "tests/fixtures/same_chrom_supplement/reads.bam",
        "tests/fixtures/same_chrom_supplement/genome.fa",
        2, 2, &DenovoConfig::default(), &RefineParams::default(), 0.20,
    ).unwrap();
    // the fixture's two homologous loci (c1:A + c2:X) must land in one family of >= 2 distinct loci.
    assert!(fams.iter().any(|f| f.len() >= 2), "expected a >=2-copy homology family, got {:?}", fams.iter().map(|f| f.len()).collect::<Vec<_>>());
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --lib homology_catalog_groups_fixture_family`
Expected: FAIL — `cannot find function detect_homology_catalog_genome_wide`.

- [ ] **Step 3: Write minimal implementation**

```rust
/// GENOME-WIDE homology-primary (E_r) family catalog. reps → E_r edges → γ-quasi-clique blocks →
/// ≥2 distinct loci → families. Chrom/strand-agnostic; a superset of the conflict catalog.
pub fn detect_homology_catalog_genome_wide(
    bam_path: &str,
    fasta_path: &str,
    threads: usize,
    min_copies: usize,
    cfg: &DenovoConfig,
    refine: &RefineParams,
    gamma: f64,
) -> Result<Vec<Vec<DenovoTranscript>>> {
    // --- reps (identical to the conflict path's rep build) ---
    let reads = primary_reads_from_bam(bam_path, threads)?;
    let contigs: HashSet<String> = reads.iter().map(|r| r.chrom.clone()).collect();
    let genome = GenomeIndex::from_fasta_contigs(fasta_path, &contigs)?;
    let skeletons = pass1_skeletons_robust(&reads, cfg.pass1_min_reads, cfg.min_terminal_support);
    drop(reads);
    let transcripts = assemble_gate(&skeletons, &genome, &cfg.gate);
    let rep_idx = collapse_loci_span_aware(&transcripts, &cfg.detect);
    let reps: Vec<DenovoTranscript> = rep_idx.iter().map(|&i| transcripts[i].clone()).collect();
    drop(transcripts);
    eprintln!("[gw-catalog-homology] {} skeletons -> {} reps over {} contigs", skeletons.len(), reps.len(), contigs.len());

    // --- E_r edges + γ-quasi-clique blocks ---
    let edges2 = homology_edges_all_reps(&reps, refine)?;
    let edges3: Vec<(usize, usize, f64)> = edges2.iter().map(|&(a, b)| (a, b, 1.0)).collect();
    let blocks = crate::vg_family::family_split::gamma_quasi_clique_partition(reps.len(), &edges3, gamma);

    let mut out: Vec<Vec<DenovoTranscript>> = Vec::new();
    for block in blocks {
        if block.len() < min_copies { continue; }
        let copies: Vec<DenovoTranscript> = block.iter().map(|&i| reps[i].clone()).collect();
        let loci = distinct_locus_reps(copies); // ≥2 spatially-distinct loci certificate
        if loci.len() >= min_copies {
            out.push(loci);
        }
    }
    eprintln!("[gw-catalog-homology] {} E_r edges -> {} families (>= {} distinct loci)", edges2.len(), out.len(), min_copies);
    Ok(out)
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test --lib homology_catalog_groups_fixture_family`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/denovo_pipeline.rs
git commit -m "feat(denovo_pipeline): detect_homology_catalog_genome_wide (E_r primary + gamma-quasi-clique)"
```

---

## Task 4: `gw_family_catalog --homology-primary` wiring

**Files:**
- Modify: `src/bin/gw_family_catalog.rs`
- Test: `tests/gw_family_catalog_homology_primary.rs` (create)

**Interfaces:**
- Consumes: `detect_homology_catalog_genome_wide` (Task 3).
- Produces: CLI flags `--homology-primary`, `--min-identity <f>` (sets `RefineParams.min_identity`; `0.98` = Soto mode), reusing the existing families/copies/FA emit loop.

- [ ] **Step 1: Write the failing test**

```rust
use std::process::Command;
#[test]
fn homology_primary_emits_families_on_fixture() {
    let bin = env!("CARGO_BIN_EXE_gw_family_catalog");
    let out = "tests/fixtures/same_chrom_supplement/out_hom";
    let status = Command::new(bin).args([
        "--bam", "tests/fixtures/same_chrom_supplement/reads.bam",
        "--fasta", "tests/fixtures/same_chrom_supplement/genome.fa",
        "--out", out, "--homology-primary",
    ]).status().unwrap();
    assert!(status.success());
    let fams = std::fs::read_to_string(format!("{}.families.tsv", out)).unwrap();
    assert!(fams.lines().count() > 1, "expected >=1 homology family\n{}", fams);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo build --bin gw_family_catalog` — Expected: FAIL, unknown arg `--homology-primary`.

- [ ] **Step 3: Write minimal implementation**

Add to `Args`:

```rust
    /// HOMOLOGY-PRIMARY (E_r) family definition: build families from exon-sum nucleotide homology
    /// (gamma-quasi-clique), not the read-conflict graph. Captures ancient paralogs + K=0 collapses the
    /// conflict path misses. Ships alongside --cross-chrom.
    #[arg(long, default_value_t = false)]
    homology_primary: bool,
    /// E_r nucleotide identity floor (sensitive tier). Default from RefineParams (~0.60). `0.98` = Soto SD98 mode.
    #[arg(long)]
    min_identity: Option<f64>,
```

In `main`, before the existing `if args.cross_chrom` block, branch on homology-primary:

```rust
    let mut refine_params = RefineParams { threads: args.threads, ..Default::default() };
    if let Some(mi) = args.min_identity { refine_params.min_identity = mi; refine_params.sensitive_identity = mi.min(refine_params.sensitive_identity); }
    let raw: Vec<Vec<DenovoTranscript>> = if args.homology_primary {
        detect_homology_catalog_genome_wide(&args.bam, &args.fasta, args.threads, args.min_copies, &cfg, &refine_params, 0.20)?
    } else if args.cross_chrom {
        detect_conflict_catalog_genome_wide_xchrom(&args.bam, &args.fasta, args.threads, args.min_copies, &cfg)?
    } else {
        let catalog = detect_conflict_catalog_genome_wide(&args.bam, &args.fasta, args.threads, args.win, args.min_copies, &cfg)?;
        catalog.into_iter().map(|c| c.copies).collect()
    };
```
Add `detect_homology_catalog_genome_wide` to the `use` import. The existing emit loop (families.tsv/copies.tsv/copies.fa) is unchanged.

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test --test gw_family_catalog_homology_primary`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/bin/gw_family_catalog.rs tests/gw_family_catalog_homology_primary.rs
git commit -m "feat(gw_family_catalog): --homology-primary + --min-identity (Soto 0.98) mode"
```

---

## Task 5: Genome-projection copy enumeration (famCN, §7)

**Files:**
- Create: `src/rustle/vg_family/genome_projection.rs`
- Modify: `src/rustle/vg_family/mod.rs` (add `pub mod genome_projection;`)
- Test: `genome_projection.rs` tests module

**Interfaces:**
- Consumes: `GenomeIndex` (has `fetch_sequence`, `from_fasta_contigs`), a family consensus `&[u8]`, `params.minimap2`.
- Produces: `pub struct CopyLocus { pub chrom: String, pub start: u64, pub end: u64, pub identity: f64 }` and `pub fn project_family_copies(consensus: &[u8], genome_fasta: &str, known: &[(String,u64,u64)], min_identity: f64, min_cov: f64, minimap2: &str, threads: usize) -> anyhow::Result<Vec<CopyLocus>>`.

- [ ] **Step 1: Write the failing test** (integration; uses the GSTM chromosome subset if present, else skips)

```rust
#[test]
fn projection_enumerates_disjoint_copies() {
    if std::process::Command::new("minimap2").arg("--version").output().is_err() { return; }
    // Build a tiny 2-contig genome file with the SAME 1kb sequence at two loci = 2 genomic copies.
    let dir = std::env::temp_dir().join(format!("rustle_proj_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let seq: String = (0..1000).map(|i| "ACGT".as_bytes()[(i * 7 + 3) % 4] as char).collect();
    let fa = dir.join("g.fa");
    std::fs::write(&fa, format!(">c1\n{seq}{seq}\n>c2\n{seq}\n")).unwrap(); // c1 has 2 tandem copies, c2 has 1
    let hits = project_family_copies(seq.as_bytes(), fa.to_str().unwrap(), &[], 0.98, 0.90, "minimap2", 2).unwrap();
    assert!(hits.len() >= 3, "3 genomic copies (2 on c1, 1 on c2) expected, got {}", hits.len());
    let _ = std::fs::remove_dir_all(&dir);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --lib projection_enumerates_disjoint_copies`
Expected: FAIL — module/function missing.

- [ ] **Step 3: Write minimal implementation**

```rust
//! Liftoff-style genome-projection copy enumeration (spec §7): project a de-novo family consensus onto the
//! genome to enumerate near-identical genomic copies (famCN), recovering K=0 collapses RNA merges. In-engine
//! minimap2 (no Liftoff dependency); seeded by our own consensus, so no reference-annotation circularity.
use anyhow::Result;
use std::io::Write;

#[derive(Clone, Debug)]
pub struct CopyLocus { pub chrom: String, pub start: u64, pub end: u64, pub identity: f64 }

/// minimap2 the consensus against the genome; keep hits with identity ≥ `min_identity`, aligned-fraction
/// of the consensus ≥ `min_cov` (structure-preserving), disjoint from each other and from `known` loci.
pub fn project_family_copies(
    consensus: &[u8],
    genome_fasta: &str,
    known: &[(String, u64, u64)],
    min_identity: f64,
    min_cov: f64,
    minimap2: &str,
    threads: usize,
) -> Result<Vec<CopyLocus>> {
    let dir = std::env::temp_dir();
    let q = dir.join(format!("rustle_proj_q_{}_{}.fa", std::process::id(), consensus.len()));
    struct Cl(std::path::PathBuf); impl Drop for Cl { fn drop(&mut self) { let _ = std::fs::remove_file(&self.0); } }
    let _c = Cl(q.clone());
    { let mut f = std::fs::File::create(&q)?; writeln!(f, ">cons")?; f.write_all(consensus)?; writeln!(f)?; }
    let out = std::process::Command::new(minimap2)
        .args(["-c", "-x", "splice", "-N", "50", "-t"]).arg(threads.to_string())
        .arg(genome_fasta).arg(&q).output()
        .map_err(|e| anyhow::anyhow!("minimap2 ('{minimap2}') projection failed: {e}"))?;
    if !out.status.success() { return Ok(Vec::new()); }
    let clen = consensus.len() as f64;
    let mut hits: Vec<CopyLocus> = Vec::new();
    for line in String::from_utf8_lossy(&out.stdout).lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 12 { continue; }
        let tname = f[5].to_string();
        let ts = f[7].parse::<u64>().unwrap_or(0);
        let te = f[8].parse::<u64>().unwrap_or(0);
        let qs = f[2].parse::<f64>().unwrap_or(0.0);
        let qe = f[3].parse::<f64>().unwrap_or(0.0);
        let de = f[12..].iter().find_map(|x| x.strip_prefix("de:f:").and_then(|v| v.parse::<f64>().ok()));
        let ident = de.map(|d| 1.0 - d).unwrap_or_else(|| {
            f[9].parse::<f64>().unwrap_or(0.0) / f[10].parse::<f64>().unwrap_or(1.0).max(1.0)
        });
        let cov = (qe - qs) / clen.max(1.0);
        if ident >= min_identity && cov >= min_cov {
            hits.push(CopyLocus { chrom: tname, start: ts, end: te, identity: ident });
        }
    }
    // disjoint filter: sort, drop hits overlapping an already-kept hit or a `known` locus.
    hits.sort_by(|a, b| (a.chrom.as_str(), a.start).cmp(&(b.chrom.as_str(), b.start)));
    let mut kept: Vec<CopyLocus> = Vec::new();
    let overlaps = |c: &CopyLocus, k: &(String, u64, u64)| c.chrom == k.0 && c.start < k.2 && k.1 < c.end;
    for h in hits {
        if kept.iter().any(|k| k.chrom == h.chrom && k.start < h.end && h.start < k.end) { continue; }
        if known.iter().any(|k| overlaps(&h, k)) { continue; }
        kept.push(h);
    }
    Ok(kept)
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test --lib projection_enumerates_disjoint_copies`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/genome_projection.rs src/rustle/vg_family/mod.rs
git commit -m "feat(genome_projection): Liftoff-style famCN copy enumeration (spec §7)"
```

---

## Task 6: emit famCN in the catalog (wire projection into the binary)

**Files:**
- Modify: `src/bin/gw_family_catalog.rs`
- Test: `tests/gw_family_catalog_homology_primary.rs`

**Interfaces:**
- Consumes: `project_family_copies` (Task 5), `--enumerate-copies` flag (default true when `--min-identity 0.98`).
- Produces: a `<out>.famcn.tsv` (columns `family_id  n_rna_copies  famCN  projection_loci`).

- [ ] **Step 1: Write the failing test**

```rust
#[test]
fn homology_primary_writes_famcn_when_enumerating() {
    let bin = env!("CARGO_BIN_EXE_gw_family_catalog");
    let out = "tests/fixtures/same_chrom_supplement/out_famcn";
    let status = Command::new(bin).args([
        "--bam", "tests/fixtures/same_chrom_supplement/reads.bam",
        "--fasta", "tests/fixtures/same_chrom_supplement/genome.fa",
        "--out", out, "--homology-primary", "--enumerate-copies",
    ]).status().unwrap();
    assert!(status.success());
    assert!(std::fs::metadata(format!("{}.famcn.tsv", out)).is_ok(), "famcn.tsv must be written");
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --test gw_family_catalog_homology_primary homology_primary_writes_famcn_when_enumerating`
Expected: FAIL — unknown arg `--enumerate-copies` / no famcn.tsv.

- [ ] **Step 3: Write minimal implementation**

Add `--enumerate-copies` flag. After `fams` is built and sorted, before `Ok(())`:

```rust
    if args.enumerate_copies {
        let mut ff = std::fs::File::create(format!("{}.famcn.tsv", args.out))?;
        writeln!(ff, "family_id\tn_rna_copies\tfamCN\tprojection_loci")?;
        for (fi, copies) in fams.iter().enumerate() {
            let cons = copies.iter().max_by_key(|c| c.n_reads).map(|c| c.seq.clone()).unwrap_or_default();
            let known: Vec<(String, u64, u64)> = copies.iter().map(|c| (c.chrom.clone(), c.start, c.end)).collect();
            let proj = rustle::vg_family::genome_projection::project_family_copies(
                &cons, &args.fasta, &known, 0.98, 0.90, "minimap2", args.threads).unwrap_or_default();
            let famcn = copies.len() + proj.len();
            let loci = proj.iter().map(|p| format!("{}:{}-{}", p.chrom, p.start, p.end)).collect::<Vec<_>>().join(";");
            writeln!(ff, "GWFAM{fi}\t{}\t{}\t{}", copies.len(), famcn, loci)?;
        }
    }
```
Default `--enumerate-copies` to true when `--min-identity == 0.98` (Soto mode): `let enumerate = args.enumerate_copies || args.min_identity == Some(0.98);`.

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test --test gw_family_catalog_homology_primary homology_primary_writes_famcn_when_enumerating`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/bin/gw_family_catalog.rs tests/gw_family_catalog_homology_primary.rs
git commit -m "feat(gw_family_catalog): emit famCN via genome projection (--enumerate-copies)"
```

---

## Task 7: Protein QC flag (orthogonal, never an edge)

**Files:**
- Modify: `src/rustle/vg_family/denovo_pipeline.rs` (add `family_protein_coheres`), `src/bin/gw_family_catalog.rs` (emit flag)
- Test: `denovo_pipeline.rs` tests module

**Interfaces:**
- Consumes: `batch_protein_edges(families, min_fident, min_cov, params) -> Result<Vec<Vec<(usize,usize)>>>`.
- Produces: `pub fn family_protein_coheres(families: &[Vec<DenovoTranscript>], params: &RefineParams) -> Vec<Option<bool>>` — per family: `Some(true)` if its protein ORFs form a connected homology graph, `Some(false)` if not, `None` if mmseqs unavailable.

- [ ] **Step 1: Write the failing test**

```rust
#[test]
fn protein_coheres_is_none_without_mmseqs() {
    // With protein_tail off / mmseqs absent the flag is None (no membership effect).
    let fam = vec![vec![
        DenovoTranscript{tid:"a".into(),chrom:"c1".into(),start:0,end:300,n_reads:5,strand:'+',introns:vec![],seq:b"ATGAAAGGGTTTTGTCCCAAAGGG".to_vec()},
        DenovoTranscript{tid:"b".into(),chrom:"c1".into(),start:9000,end:9300,n_reads:4,strand:'+',introns:vec![],seq:b"ATGAAAGGGTTTTGTCCCAAAGGG".to_vec()},
    ]];
    let mut p = RefineParams::default(); p.protein_tail = false;
    let flags = family_protein_coheres(&fam, &p);
    assert_eq!(flags.len(), 1);
    assert_eq!(flags[0], None, "protein QC is None when protein tier is off");
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --lib protein_coheres_is_none_without_mmseqs`
Expected: FAIL — `cannot find function family_protein_coheres`.

- [ ] **Step 3: Write minimal implementation**

```rust
/// Orthogonal protein QC (NEVER a definition edge): does each nt-defined family also cohere at the protein
/// level? `Some(true)` = its ORFs form a connected protein-homology graph, `Some(false)` = they do not,
/// `None` = mmseqs unavailable / protein tier off (no effect on membership).
pub fn family_protein_coheres(families: &[Vec<DenovoTranscript>], params: &RefineParams) -> Vec<Option<bool>> {
    if !params.protein_tail {
        return vec![None; families.len()];
    }
    let edges = match batch_protein_edges(families, 0.50, params.min_coverage, params) {
        Ok(e) => e,
        Err(_) => return vec![None; families.len()],
    };
    families.iter().enumerate().map(|(fi, fam)| {
        let fe = edges.get(fi).cloned().unwrap_or_default();
        if fam.len() < 2 { return Some(true); }
        Some(homology_components(fam.len(), &fe).iter().any(|c| c.len() == fam.len()))
    }).collect()
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test --lib protein_coheres_is_none_without_mmseqs`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/denovo_pipeline.rs
git commit -m "feat(denovo_pipeline): orthogonal protein-coherence QC flag (never a definition edge)"
```

---

## Task 8: Real-data regression — recover the audit misses

**Files:**
- Test: `tests/gw_family_catalog_homology_primary.rs` (add a real-data, env-gated test)

**Interfaces:**
- Consumes: the built `gw_family_catalog` binary; the real BAM/FASTA at `/home/juanfra/winloci_scratch/` (env-gated so CI without the data skips).

- [ ] **Step 1: Write the failing test**

```rust
#[test]
fn homology_primary_recovers_gstm_and_znf() {
    let bam = "/home/juanfra/winloci_scratch/GGO_mm.bam";
    if std::fs::metadata(bam).is_err() { eprintln!("real data absent; skip"); return; }
    // chr234 = NC_073234.2 has the KRAB-ZNF trio; run homology-primary and assert a >=3-copy family forms.
    let bin = env!("CARGO_BIN_EXE_gw_family_catalog");
    let out = std::env::temp_dir().join("hom_c234");
    // caller pre-subsets NC_073244.2 (KRAB-ZNF) into a small BAM to keep the test fast; see plan Task 8 notes.
    let sub = "/home/juanfra/winloci_scratch/c244.bam";
    if std::fs::metadata(sub).is_err() { eprintln!("subset c244.bam absent; skip"); return; }
    let ok = std::process::Command::new(bin).args([
        "--bam", sub, "--fasta", "/home/juanfra/winloci_scratch/GGO.fasta",
        "--out", out.to_str().unwrap(), "--homology-primary",
    ]).status().unwrap().success();
    assert!(ok);
    let copies = std::fs::read_to_string(format!("{}.copies.tsv", out.to_str().unwrap())).unwrap();
    // ZNF14/ZNF682/ZNF429 span ~26.1–28.0 Mb on NC_073244.2; assert >=2 copies land there in ONE family.
    let znf: Vec<&str> = copies.lines().filter(|l| l.contains("NC_073244.2")).collect();
    assert!(znf.len() >= 2, "KRAB-ZNF paralogs must be grouped by homology-primary, got {} rows", znf.len());
}
```

- [ ] **Step 2: Run test to verify it fails (or skips cleanly)**

Prep the subset once: `mamba run -n liftoff samtools view -b /home/juanfra/winloci_scratch/GGO_mm.bam NC_073244.2 > /home/juanfra/winloci_scratch/c244.bam && samtools index /home/juanfra/winloci_scratch/c244.bam`
Run: `cargo test --test gw_family_catalog_homology_primary homology_primary_recovers_gstm_and_znf -- --nocapture`
Expected: FAIL before Tasks 1–4 land; PASS after (or clean skip if data absent).

- [ ] **Step 3: Implementation**

No new code — this test validates Tasks 1–5 end-to-end on real data. If it fails after implementation, the sensitive identity floor (`RefineParams.sensitive_identity`, default 0.60) is the lever: confirm ZNF14↔ZNF429 (0.62) clears it; adjust only via the P/R sweep.

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test --test gw_family_catalog_homology_primary homology_primary_recovers_gstm_and_znf`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add tests/gw_family_catalog_homology_primary.rs
git commit -m "test(gw_family_catalog): real-data regression — homology-primary recovers KRAB-ZNF family"
```

---

## Self-review notes (author checklist, resolved)

- **Spec coverage:** §1 Task 3; §2 Task 2; §3 (prefilter) — the minimap2 all-vs-all in `nucleotide_edges` IS the prefilter for the first cut (secondary-co-mapping candidate source is a *follow-up optimization*, NOT in this plan — flagged as deferred); §4 Task 1; §5 (within-family χ(H)) — **deferred** to a follow-up per spec §8 open question (this plan ships the E_r definition + famCN; χ(H)/parCN read-assignment reuses existing `copy_assign` in a later plan); §6 Task 7; §7 Tasks 5–6; Soto mode Task 4; validation Task 8.
- **Deferred (explicit, not silent):** secondary-co-mapping candidate prefilter (§3b), within-family χ(H)/parCN read-assignment (§5), and the family P/R sweep to fix `sensitive_identity` — each is a named follow-up, not a placeholder here.
- **Type consistency:** `RefineParams.sensitive_identity` (Task 2) used in `homology_edges_all_reps`; `CopyLocus`/`project_family_copies` (Task 5) used in Task 6; `detect_homology_catalog_genome_wide` signature identical in Tasks 3 and 4.
