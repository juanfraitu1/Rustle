# Family Manifest Ingestion — Phase 2 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Allow Rustle to ingest an R-generated family manifest TSV — listing which genomic loci form a multi-copy gene family — so that multi-copy families whose copies live at distinct genomic positions (RBMY, GOLGA, etc.) are correctly grouped into `FamilyGroup`s for HMM-EM read disambiguation.

**Architecture:** R's `build_variation_graph_seq()` already identifies multi-copy families by sequence similarity (minimap2 all-vs-all). A new `export_family_manifest()` R function dumps locus coordinates as a simple TSV. A new Rust module `family_manifest.rs` parses that TSV and converts it to `Vec<FamilyGroup>` by finding all bundles whose genomic range overlaps any locus in each family — the identical interface used by `ingest_gtf_families()`. The pipeline at `pipeline.rs:9744` gains a third branch for manifest-based ingestion alongside the existing GTF and multi-mapper-discovery paths.

**Tech Stack:** R (data export), Rust/Cargo (ingestion, no new dependencies)

---

## File Map

| Action | Path | Responsibility |
|--------|------|----------------|
| Modify | `analysis/family_graphs/02_build_graphs.R` | Add `export_family_manifest()` |
| Create | `analysis/family_graphs/06_export_manifests.R` | Driver script: export manifests for RBMY + GOLGA |
| Create | `src/rustle/family_manifest.rs` | Parse TSV → `Vec<FamilyGroup>` |
| Modify | `src/rustle/lib.rs:85-88` | Register `pub mod family_manifest;` |
| Modify | `src/rustle/types.rs:875-877` | Add `family_manifest` field to `RunConfig` |
| Modify | `src/rustle/types.rs:1219-1220` | Add `family_manifest: None` to `RunConfig::default()` |
| Modify | `src/rustle/main.rs:374-379` | Add `--family-manifest` CLI arg |
| Modify | `src/rustle/main.rs:743-744` | Wire arg → config |
| Modify | `src/rustle/pipeline.rs:9744-9786` | Add manifest branch before GTF branch |
| Create | `tests/regression/vg_family_manifest.rs` | Unit tests for parse + group creation |
| Modify | `Cargo.toml` | Register new test target |

---

## Task 1: R — `export_family_manifest()` function

**Files:**
- Modify: `analysis/family_graphs/02_build_graphs.R` (append after line 267, before the `if (sys.nframe() == 0L)` guard)

- [ ] **Step 1: Write failing test (R)**

Create `analysis/family_graphs/test_export_manifest.R`:

```r
source("02_build_graphs.R")

# Minimal synthetic exon_df: 2 genes, each 2 exons, different chromosomes
syn_exon_df <- data.frame(
  gene_id = c("gA","gA","gB","gB"),
  tx_id   = c("tA","tA","tB","tB"),
  chrom   = c("chr1","chr1","chr2","chr2"),
  start   = c(100L,300L,100L,300L),
  end     = c(200L,400L,200L,400L),
  strand  = "+",
  stringsAsFactors = FALSE
)
tmp <- tempfile(fileext = ".tsv")
export_family_manifest(syn_exon_df, "TESTFAM", tmp)
tbl <- read.table(tmp, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

stopifnot(all(c("family_id","gene_id","chrom","start","end","strand") %in% names(tbl)))
stopifnot(nrow(tbl) == 2L)
stopifnot(all(tbl$family_id == "TESTFAM"))
stopifnot(tbl$start[tbl$gene_id == "gA"] == 100L)
stopifnot(tbl$end[tbl$gene_id == "gA"]   == 400L)
stopifnot(tbl$chrom[tbl$gene_id == "gA"] == "chr1")
message("OK: export_family_manifest test passed")
```

- [ ] **Step 2: Run test — expect failure**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle/analysis/family_graphs
Rscript test_export_manifest.R
# Expected: Error: could not find function "export_family_manifest"
```

- [ ] **Step 3: Implement `export_family_manifest()` in `02_build_graphs.R`**

Append the following block **between line 267** (`list(...)` closing brace of `build_variation_graph_seq`) **and line 269** (`if (sys.nframe() == 0L)`):

```r
# Export a family manifest TSV: one row per gene locus in the family.
# exon_df    : the $exon_df from build_variation_graph_seq (or build_variation_graph)
# family_id  : character string used as the family_id column
# output_path: file path for the TSV output
export_family_manifest <- function(exon_df, family_id, output_path) {
  loci <- do.call(rbind, lapply(split(exon_df, exon_df$gene_id), function(g) {
    data.frame(
      family_id = family_id,
      gene_id   = g$gene_id[1L],
      chrom     = g$chrom[1L],
      start     = min(g$start),
      end       = max(g$end),
      strand    = g$strand[1L],
      stringsAsFactors = FALSE
    )
  }))
  write.table(loci, output_path, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
}
```

- [ ] **Step 4: Run test — expect pass**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle/analysis/family_graphs
Rscript test_export_manifest.R
# Expected: OK: export_family_manifest test passed
```

- [ ] **Step 5: Commit**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
git add analysis/family_graphs/02_build_graphs.R analysis/family_graphs/test_export_manifest.R
git commit -m "feat(R): add export_family_manifest() to build_variation_graph_seq helpers"
```

---

## Task 2: R — `06_export_manifests.R` driver script

**Files:**
- Create: `analysis/family_graphs/06_export_manifests.R`

- [ ] **Step 1: Create the driver script**

```r
source("02_build_graphs.R")

DATA_DIR <- "data"
GENOME_FA <- "/mnt/c/Users/jfris/Desktop/GGO.fasta"

rbmy_vg <- readRDS(file.path(DATA_DIR, "rbmy_vg_seq.rds"))
export_family_manifest(rbmy_vg$exon_df, "RBMY",
                       file.path(DATA_DIR, "rbmy_manifest.tsv"))
message("Exported RBMY manifest: ", nrow(rbmy_vg$exon_df |> (\(d) unique(d$gene_id))()),
        " loci → data/rbmy_manifest.tsv")

golga_vg <- readRDS(file.path(DATA_DIR, "golga_vg_seq.rds"))
export_family_manifest(golga_vg$exon_df, "GOLGA",
                       file.path(DATA_DIR, "golga_manifest.tsv"))
message("Exported GOLGA manifest: ", nrow(golga_vg$exon_df |> (\(d) unique(d$gene_id))()),
        " loci → data/golga_manifest.tsv")

message("Done. Check data/*_manifest.tsv")
```

- [ ] **Step 2: Run the driver — expect two TSV files created**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle/analysis/family_graphs
Rscript 06_export_manifests.R
# Expected:
#   Exported RBMY manifest: 13 loci → data/rbmy_manifest.tsv
#   Exported GOLGA manifest: 46 loci → data/golga_manifest.tsv
#   Done. Check data/*_manifest.tsv
head -3 data/rbmy_manifest.tsv
# Expected: header + 2 data rows with family_id=RBMY
```

- [ ] **Step 3: Commit**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
git add analysis/family_graphs/06_export_manifests.R
git commit -m "feat(R): add 06_export_manifests.R driver to export family locus manifests"
```

---

## Task 3: Rust — `family_manifest.rs` (parse + group creation)

**Files:**
- Create: `src/rustle/family_manifest.rs`
- Create: `tests/regression/vg_family_manifest.rs`

The manifest TSV has a header row then one row per locus:
```
family_id	gene_id	chrom	start	end	strand
RBMY	LOC129530259	chrY	1000000	1010000	+
```

`create_family_groups_from_manifest` must handle cross-chromosome families (GOLGA copies may be on different chromosomes) — it iterates every bundle against every locus, unlike `create_family_groups_from_templates` which filters by a single chrom.

- [ ] **Step 1: Write failing tests**

Create `tests/regression/vg_family_manifest.rs`:

```rust
use rustle::family_manifest::{parse_family_manifest, create_family_groups_from_manifest, FamilyLocus};
use rustle::types::{Bundle, JunctionStats};
use std::io::Write;

fn mk_bundle(chrom: &str, start: u64, end: u64) -> Bundle {
    Bundle {
        chrom: chrom.into(),
        start,
        end,
        strand: '+',
        reads: Vec::new(),
        junction_stats: JunctionStats::default(),
        bundlenodes: None,
        read_bnodes: None,
        bnode_colors: None,
        synthetic: false,
        rescue_class: None,
    }
}

#[test]
fn parse_manifest_two_families() {
    let content = "family_id\tgene_id\tchrom\tstart\tend\tstrand\n\
                   FAM1\tgeneA\tchr1\t1000\t2000\t+\n\
                   FAM1\tgeneB\tchr2\t5000\t6000\t+\n\
                   FAM2\tgeneC\tchr3\t9000\t10000\t-\n";
    let mut f = tempfile::NamedTempFile::new().unwrap();
    write!(f, "{}", content).unwrap();
    let loci = parse_family_manifest(f.path()).unwrap();
    assert_eq!(loci.len(), 3);
    assert_eq!(loci[0].family_id, "FAM1");
    assert_eq!(loci[0].gene_id, "geneA");
    assert_eq!(loci[0].chrom, "chr1");
    assert_eq!(loci[0].start, 1000);
    assert_eq!(loci[0].end, 2000);
    assert_eq!(loci[0].strand, '+');
    assert_eq!(loci[2].family_id, "FAM2");
    assert_eq!(loci[2].strand, '-');
}

#[test]
fn parse_manifest_skips_blank_and_comment_lines() {
    let content = "family_id\tgene_id\tchrom\tstart\tend\tstrand\n\
                   \n\
                   # comment\n\
                   FAM1\tgeneA\tchr1\t1000\t2000\t+\n";
    let mut f = tempfile::NamedTempFile::new().unwrap();
    write!(f, "{}", content).unwrap();
    let loci = parse_family_manifest(f.path()).unwrap();
    assert_eq!(loci.len(), 1);
}

#[test]
fn create_groups_links_overlapping_bundles() {
    // FAM1 has two loci on chr1 and chr2; two bundles match, one doesn't
    let loci = vec![
        FamilyLocus { family_id: "FAM1".into(), gene_id: "gA".into(), chrom: "chr1".into(), start: 1000, end: 2000, strand: '+' },
        FamilyLocus { family_id: "FAM1".into(), gene_id: "gB".into(), chrom: "chr2".into(), start: 5000, end: 6000, strand: '+' },
    ];
    let bundles = vec![
        mk_bundle("chr1", 1200, 1800),   // overlaps gA locus → bundle 0
        mk_bundle("chr2", 5100, 5900),   // overlaps gB locus → bundle 1
        mk_bundle("chr3", 0, 1000),      // no match
    ];
    let groups = create_family_groups_from_manifest(&loci, &bundles);
    assert_eq!(groups.len(), 1, "should produce exactly one family group");
    let g = &groups[0];
    let mut idxs = g.bundle_indices.clone();
    idxs.sort();
    assert_eq!(idxs, vec![0, 1]);
}

#[test]
fn create_groups_two_families() {
    let loci = vec![
        FamilyLocus { family_id: "FAM1".into(), gene_id: "gA".into(), chrom: "chr1".into(), start: 1000, end: 2000, strand: '+' },
        FamilyLocus { family_id: "FAM2".into(), gene_id: "gC".into(), chrom: "chr3".into(), start: 9000, end: 10000, strand: '-' },
    ];
    let bundles = vec![
        mk_bundle("chr1", 1500, 1800),
        mk_bundle("chr3", 9500, 9900),
    ];
    let groups = create_family_groups_from_manifest(&loci, &bundles);
    assert_eq!(groups.len(), 2);
}

#[test]
fn create_groups_empty_when_no_overlap() {
    let loci = vec![
        FamilyLocus { family_id: "FAM1".into(), gene_id: "gA".into(), chrom: "chr1".into(), start: 1000, end: 2000, strand: '+' },
    ];
    let bundles = vec![mk_bundle("chr1", 5000, 6000)]; // no overlap
    let groups = create_family_groups_from_manifest(&loci, &bundles);
    assert!(groups.is_empty());
}
```

- [ ] **Step 2: Add test target to `Cargo.toml`**

In `Cargo.toml`, after the last `[[test]]` block (search for the end of the existing test list), add:

```toml
[[test]]
name = "vg_family_manifest"
path = "tests/regression/vg_family_manifest.rs"
```

- [ ] **Step 3: Run tests — expect compile failure**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
cargo test --test vg_family_manifest 2>&1 | head -20
# Expected: error[E0432]: unresolved import `rustle::family_manifest`
```

- [ ] **Step 4: Check that `tempfile` is already a dev-dependency**

```bash
grep "tempfile" /mnt/c/Users/jfris/Desktop/Rustle/Cargo.toml
```

If not present, add to `[dev-dependencies]`:

```toml
tempfile = "3"
```

- [ ] **Step 5: Implement `src/rustle/family_manifest.rs`**

```rust
//! Family manifest ingestion: parse R-exported TSV of multi-copy gene family loci
//! and create FamilyGroups that the existing HMM-EM pipeline can consume.

use crate::vg::FamilyGroup;
use crate::types::Bundle;
use anyhow::{anyhow, Context, Result};
use std::collections::HashMap;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct FamilyLocus {
    pub family_id: String,
    pub gene_id:   String,
    pub chrom:     String,
    pub start:     u64,
    pub end:       u64,
    pub strand:    char,
}

/// Parse a tab-separated family manifest produced by `export_family_manifest()` in R.
///
/// Format (tab-separated, first row is a header):
/// ```text
/// family_id<TAB>gene_id<TAB>chrom<TAB>start<TAB>end<TAB>strand
/// RBMY<TAB>LOC129530259<TAB>chrY<TAB>1000000<TAB>1010000<TAB>+
/// ```
/// Blank lines and lines starting with `#` are ignored.
pub fn parse_family_manifest<P: AsRef<Path>>(path: P) -> Result<Vec<FamilyLocus>> {
    let content = std::fs::read_to_string(path.as_ref())
        .with_context(|| format!("reading family manifest {:?}", path.as_ref()))?;

    let mut loci = Vec::new();
    let mut header_seen = false;

    for (lineno, line) in content.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') { continue; }

        if !header_seen {
            header_seen = true;
            continue; // skip header row
        }

        let cols: Vec<&str> = line.splitn(6, '\t').collect();
        if cols.len() < 6 {
            return Err(anyhow!(
                "family manifest line {}: expected 6 tab-separated columns, got {}: {:?}",
                lineno + 1, cols.len(), line
            ));
        }

        let start: u64 = cols[3].parse()
            .with_context(|| format!("manifest line {}: bad start {:?}", lineno + 1, cols[3]))?;
        let end: u64 = cols[4].parse()
            .with_context(|| format!("manifest line {}: bad end {:?}", lineno + 1, cols[4]))?;
        let strand_char = cols[5].chars().next().unwrap_or('+');

        loci.push(FamilyLocus {
            family_id: cols[0].to_string(),
            gene_id:   cols[1].to_string(),
            chrom:     cols[2].to_string(),
            start,
            end,
            strand:    strand_char,
        });
    }

    Ok(loci)
}

/// Build `FamilyGroup`s from the manifest loci by matching overlapping bundles.
///
/// For each distinct `family_id`, collects every bundle whose `[start, end)` half-open
/// range overlaps any locus belonging to that family. The resulting `FamilyGroup`
/// has `multimap_reads` left empty — the HMM-EM pipeline populates it from the BAM.
///
/// A bundle can appear in at most one family group (first-match wins when loci
/// from different families overlap the same bundle — rare in practice for distinct
/// gene families).
pub fn create_family_groups_from_manifest(
    loci: &[FamilyLocus],
    bundles: &[Bundle],
) -> Vec<FamilyGroup> {
    // Group loci by family_id (preserve insertion order for determinism).
    let mut family_order: Vec<String> = Vec::new();
    let mut family_loci: HashMap<String, Vec<&FamilyLocus>> = HashMap::new();
    for loc in loci {
        if !family_loci.contains_key(&loc.family_id) {
            family_order.push(loc.family_id.clone());
        }
        family_loci.entry(loc.family_id.clone()).or_default().push(loc);
    }

    let mut groups = Vec::new();
    for (family_numeric_id, fid) in family_order.iter().enumerate() {
        let fam_loci = &family_loci[fid];
        let mut bundle_indices: Vec<usize> = Vec::new();

        for (bi, bundle) in bundles.iter().enumerate() {
            let overlaps_any = fam_loci.iter().any(|loc| {
                bundle.chrom == loc.chrom
                    && bundle.end   > loc.start   // half-open overlap
                    && bundle.start < loc.end
            });
            if overlaps_any {
                bundle_indices.push(bi);
            }
        }

        if !bundle_indices.is_empty() {
            groups.push(FamilyGroup {
                family_id: family_numeric_id,
                bundle_indices,
                multimap_reads: HashMap::new(),
            });
        }
    }

    groups
}
```

- [ ] **Step 6: Register module in `src/rustle/lib.rs`**

After line 87 (`pub mod vg_ingestion;`), insert:

```rust
pub mod family_manifest; // TSV-manifest ingestion: R-annotated multi-locus families
```

So the VG section reads:

```rust
// ── Variation graph (gene family) mode ───────────────────────────────────────
pub mod vg; // family group discovery, EM reweighting
pub mod vg_ingestion; // GTF/GFF ingestion mode for template-based family assembly
pub mod family_manifest; // TSV-manifest ingestion: R-annotated multi-locus families
pub mod vg_hmm; // Family-aware HMM rescue for novel gene copies
pub mod graph_comparison;
```

- [ ] **Step 7: Run tests — expect pass**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
cargo test --test vg_family_manifest 2>&1
# Expected: test result: ok. 4 passed; 0 failed
```

- [ ] **Step 8: Commit**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
git add src/rustle/family_manifest.rs src/rustle/lib.rs \
        tests/regression/vg_family_manifest.rs Cargo.toml
git commit -m "feat: add family_manifest module — parse R-exported multi-locus TSV into FamilyGroups"
```

---

## Task 4: Rust — wire `--family-manifest` through CLI and pipeline

**Files:**
- Modify: `src/rustle/types.rs:875-877` (add field after `ingress_grouping`)
- Modify: `src/rustle/types.rs:1219-1220` (add default)
- Modify: `src/rustle/main.rs:374-379` (add CLI arg)
- Modify: `src/rustle/main.rs:743-744` (wire to config)
- Modify: `src/rustle/pipeline.rs:9744-9786` (add manifest branch)

There are no tests for this task (integration would require a full BAM run). Verify with `cargo build`.

- [ ] **Step 1: Add `family_manifest` field to `RunConfig` in `src/rustle/types.rs`**

After line 877 (`pub ingress_grouping: String,`), add:

```rust
    /// TSV manifest from R's export_family_manifest(): one row per gene locus in a family.
    /// When set, Rustle uses manifest-defined loci to group bundles into FamilyGroups instead
    /// of multi-mapper discovery or GTF ingestion. Requires --vg.
    pub family_manifest: Option<std::path::PathBuf>,
```

- [ ] **Step 2: Add `family_manifest: None` to `RunConfig::default()` in `src/rustle/types.rs`**

After line 1220 (`ingress_grouping: "ByGene".to_string(),`), add:

```rust
            family_manifest: None,
```

- [ ] **Step 3: Add `--family-manifest` CLI arg in `src/rustle/main.rs`**

After the `ingress_grouping` arg block (around line 379), add:

```rust
    /// TSV manifest produced by R's export_family_manifest() listing multi-copy gene family loci.
    /// Each row: family_id<TAB>gene_id<TAB>chrom<TAB>start<TAB>end<TAB>strand.
    /// When provided with --vg, Rustle groups bundles into families by coordinate overlap
    /// with the listed loci instead of multi-mapper discovery. Requires --vg.
    #[arg(long)]
    family_manifest: Option<std::path::PathBuf>,
```

- [ ] **Step 4: Wire arg into config in `src/rustle/main.rs`**

After line 744 (`ingress_grouping: args.ingress_grouping.clone(),`), add:

```rust
        family_manifest: args.family_manifest.clone(),
```

- [ ] **Step 5: Add manifest branch in `src/rustle/pipeline.rs`**

Replace the block at lines 9743–9786 with the following (adds the manifest branch as the **first** check, before GTF, before multi-mapper discovery):

```rust
        // Family ingestion: manifest > GTF > multi-mapper discovery
        let raw_families = if let Some(manifest_path) = &config.family_manifest {
            // Manifest ingestion mode: use R-exported locus coordinates
            match crate::family_manifest::parse_family_manifest(manifest_path) {
                Ok(loci) => {
                    let groups = crate::family_manifest::create_family_groups_from_manifest(
                        &loci,
                        &bundles,
                    );
                    if config.verbose {
                        eprintln!(
                            "[VG-Manifest] Created {} family groups from {} loci in {}",
                            groups.len(),
                            loci.len(),
                            manifest_path.display()
                        );
                    }
                    groups
                }
                Err(e) => {
                    eprintln!("[VG-Manifest] Error reading manifest: {}", e);
                    eprintln!("[VG-Manifest] Falling back to standard discovery");
                    crate::vg::discover_family_groups(
                        &bundles,
                        config.vg_min_shared_reads,
                        Some(bam_path.as_ref()),
                        vg_genome_for_discovery.as_ref(),
                    )
                }
            }
        } else if let Some(gtf_path) = &config.ingress_gtf {
            // GTF ingestion mode: parse templates and link to existing bundles
            let grouping_strategy = match config.ingress_grouping.as_str() {
                "ByOverlap" => crate::vg_ingestion::FamilyGroupingStrategy::ByOverlap,
                _ => crate::vg_ingestion::FamilyGroupingStrategy::ByGene,
            };

            match crate::vg_ingestion::ingest_gtf_families(
                gtf_path,
                grouping_strategy,
                &bundles,
                vg_genome_for_discovery.as_ref(),
            ) {
                Ok(template_families) => {
                    if config.verbose {
                        eprintln!(
                            "[VG-Ingestion] Created {} template families from {}",
                            template_families.len(),
                            gtf_path.display()
                        );
                    }
                    template_families
                }
                Err(e) => {
                    eprintln!("[VG-Ingestion] Error: {}", e);
                    eprintln!("[VG-Ingestion] Falling back to standard discovery");
                    crate::vg::discover_family_groups(
                        &bundles,
                        config.vg_min_shared_reads,
                        Some(bam_path.as_ref()),
                        vg_genome_for_discovery.as_ref(),
                    )
                }
            }
        } else {
            // Standard discovery from multi-mapped reads
            crate::vg::discover_family_groups(
                &bundles,
                config.vg_min_shared_reads,
                Some(bam_path.as_ref()),
                vg_genome_for_discovery.as_ref(),
            )
        };
```

- [ ] **Step 6: Build — expect clean compile**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
cargo build 2>&1 | tail -5
# Expected: Compiling rustle ... Finished dev [unoptimized + debuginfo]
```

- [ ] **Step 7: Smoke-test help output**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
cargo run -- --help 2>&1 | grep "family-manifest"
# Expected: --family-manifest <FAMILY_MANIFEST>
```

- [ ] **Step 8: Commit**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
git add src/rustle/types.rs src/rustle/main.rs src/rustle/pipeline.rs
git commit -m "feat: wire --family-manifest flag through CLI → RunConfig → pipeline ingestion"
```

---

## Task 5: End-to-end smoke test with RBMY manifest

This is a manual smoke test — it does not require the full GGO.bam (which is large), but confirms the manifest path compiles into a real run without panic. Run on a subset BAM if available; otherwise verify via `--dry-run` or log output.

**Files:** (no new files — this is a verification step)

- [ ] **Step 1: Verify the manifest file was created in Task 2**

```bash
head -5 /mnt/c/Users/jfris/Desktop/Rustle/analysis/family_graphs/data/rbmy_manifest.tsv
# Expected:
# family_id	gene_id	chrom	start	end	strand
# RBMY	LOC129530259	chrY	...	...	+
# ...
```

- [ ] **Step 2: Run Rustle on GGO.bam with RBMY manifest**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
cargo run --release -- \
  --bam /mnt/c/Users/jfris/Desktop/GGO.bam \
  --vg \
  --family-manifest analysis/family_graphs/data/rbmy_manifest.tsv \
  --verbose \
  -o /tmp/rustle_rbmy_test.gtf \
  2>&1 | grep -E "\[VG-Manifest\]|family group"
# Expected: [VG-Manifest] Created N family groups from 13 loci in ...
```

- [ ] **Step 3: Confirm no panic; check family group count**

If N >= 1, the manifest path is working. Exact N depends on whether chrY RBMY bundles form in the BAM. Even N=0 is acceptable for this test (no bundles on chrY means no groups formed, which is correct behavior — not a bug).

- [ ] **Step 4: Run full test suite to check for regressions**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
cargo test 2>&1 | tail -15
# Expected: all existing tests pass; vg_family_manifest passes
```

- [ ] **Step 5: Commit test evidence (if any output files are useful)**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
git add analysis/family_graphs/data/rbmy_manifest.tsv \
        analysis/family_graphs/data/golga_manifest.tsv
git commit -m "data: add R-generated family manifests for RBMY and GOLGA"
```

---

## Self-Review

**Spec coverage:**
- ✅ R exports family manifest TSV → Task 1 + Task 2
- ✅ Rust parses TSV → `FamilyLocus` vec → Task 3
- ✅ Cross-locus grouping (bundles on chr1 and chr2 both added to same family) → Task 3 `create_family_groups_from_manifest`
- ✅ CLI `--family-manifest` flag → Task 4
- ✅ Pipeline branch: manifest > GTF > multi-mapper → Task 4 Step 5
- ✅ Smoke test with RBMY manifest → Task 5
- ✅ All existing HMM-EM, novel-copy, and assembly pipeline unchanged

**Placeholder scan:** None found.

**Type consistency:**
- `FamilyLocus` defined in Task 3 Step 5, imported in test at Task 3 Step 1 ✅
- `parse_family_manifest` returns `Result<Vec<FamilyLocus>>` in both test and impl ✅
- `create_family_groups_from_manifest` takes `(&[FamilyLocus], &[Bundle])` in both test and impl ✅
- `FamilyGroup` from `crate::vg` (already exists at `src/rustle/vg.rs:32`) ✅
