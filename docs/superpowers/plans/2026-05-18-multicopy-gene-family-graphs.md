# Multi-Copy Gene Family Graphs — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build annotation-derived variation graphs for AMY, GOLGA, and NPIP gene families from the GGO gorilla GFF, validate a formal definition of multi-copy gene families, and produce three publication-quality figures.

**Architecture:** Parse GFF with rtracklayer, disjoin exon intervals to build exon-equivalence-class nodes, construct igraph DAGs with shared/copy-specific labels, visualize with ggraph + patchwork. AMY first (validation), then GOLGA and NPIP.

**Tech Stack:** R, rtracklayer, GenomicRanges, igraph, ggraph, ggplot2, patchwork, testthat — all installed via mamba.

---

## File Map

```
Rustle/
  analysis/family_graphs/
    00_install.sh             mamba installs — run once
    01_parse_gff.R            GFF → per-family exon tables (RDS)
    02_build_graphs.R         exon table → igraph DAG + node/path metadata (RDS)
    03_figure_problem.R       Figure 1 — multi-mapping ambiguity (3 panes)
    04_figure_method.R        Figure 2 — variation graph + copy paths (3 panes)
    05_figure_examples.R      Figure 3 — examples & counterexamples (3 panes)
    tests/
      test_build_graphs.R     testthat unit tests for graph construction
    data/                     generated RDS (gitignored)
    figures/                  generated PDF + PNG (gitignored)
```

---

## Task 0: Create directories and install packages

**Files:**
- Create: `analysis/family_graphs/00_install.sh`
- Create: `analysis/family_graphs/data/.gitkeep`
- Create: `analysis/family_graphs/figures/.gitkeep`
- Create: `analysis/family_graphs/tests/.gitkeep`

- [ ] **Step 1: Create directory tree**

```bash
mkdir -p /mnt/c/Users/jfris/Desktop/Rustle/analysis/family_graphs/{data,figures,tests}
touch /mnt/c/Users/jfris/Desktop/Rustle/analysis/family_graphs/data/.gitkeep
touch /mnt/c/Users/jfris/Desktop/Rustle/analysis/family_graphs/figures/.gitkeep
touch /mnt/c/Users/jfris/Desktop/Rustle/analysis/family_graphs/tests/.gitkeep
```

- [ ] **Step 2: Write install script**

Write `analysis/family_graphs/00_install.sh`:

```bash
#!/usr/bin/env bash
set -euo pipefail

mamba install -y \
  -c bioconda -c conda-forge \
  bioconductor-rtracklayer \
  bioconductor-genomicranges \
  r-igraph \
  r-ggraph \
  r-ggplot2 \
  r-dplyr \
  r-tidyr \
  r-stringr \
  r-patchwork \
  r-testthat
```

- [ ] **Step 3: Run install script**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle/analysis/family_graphs
bash 00_install.sh
```

Expected: mamba resolves environment and installs all packages without conflict.

- [ ] **Step 4: Smoke-test R imports**

```bash
Rscript -e "library(rtracklayer); library(igraph); library(ggraph); library(patchwork); cat('OK\n')"
```

Expected output: `OK`

- [ ] **Step 5: Commit**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
git add analysis/family_graphs/
git commit -m "feat: scaffold family_graphs analysis directory"
```

---

## Task 1: Cache GFF as RDS

**Files:**
- Create: `analysis/family_graphs/01_parse_gff.R`

The GFF is 662 MB. Loading it once and caching as RDS saves ~2 minutes on every subsequent run.

- [ ] **Step 1: Write caching block at top of 01_parse_gff.R**

```r
suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(dplyr)
  library(stringr)
})

GFF_PATH  <- "/mnt/c/Users/jfris/Desktop/GGO_genomic.gff"
DATA_DIR  <- "data"

gff_cache <- file.path(DATA_DIR, "gff_cache.rds")

if (!file.exists(gff_cache)) {
  message("Importing GFF (this takes ~2 min)...")
  gff <- import(GFF_PATH)
  saveRDS(gff, gff_cache)
  message("Cached to ", gff_cache)
} else {
  message("Loading cached GFF...")
  gff <- readRDS(gff_cache)
}
message("GFF loaded: ", length(gff), " features")
```

- [ ] **Step 2: Run and verify**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle/analysis/family_graphs
Rscript 01_parse_gff.R
```

Expected output (approximately):
```
Importing GFF (this takes ~2 min)...
Cached to data/gff_cache.rds
GFF loaded: 3500000 features
```

---

## Task 2: Write parse_family() and extract AMY

**Files:**
- Modify: `analysis/family_graphs/01_parse_gff.R`

`parse_family()` takes the cached GFF and a description regex, returns a list with `$genes` (GRanges), `$exon_df` (data.frame). One row per exon, columns: `gene_id`, `tx_id`, `chrom`, `start`, `end`, `strand`.

- [ ] **Step 1: Add helper to unwrap CharacterList Parents**

Append to `01_parse_gff.R`:

```r
# GFF3 Parent attributes come as CharacterList; unwrap to plain character
first_parent <- function(x) {
  sapply(as.list(x), function(p) if (length(p) > 0L) p[[1L]] else NA_character_)
}
```

- [ ] **Step 2: Add parse_family() function**

```r
parse_family <- function(gff, description_pattern) {
  # -- genes --
  genes <- gff[gff$type == "gene"]
  desc  <- mcols(genes)$description
  desc[is.na(desc)] <- ""
  family_genes <- genes[grepl(description_pattern, desc, ignore.case = TRUE)]
  gene_ids     <- family_genes$ID
  message("  genes found: ", length(gene_ids))

  # -- transcripts (mRNA or transcript) --
  txs <- gff[gff$type %in% c("mRNA", "transcript")]
  tx_parent <- first_parent(txs$Parent)
  family_txs <- txs[!is.na(tx_parent) & tx_parent %in% gene_ids]
  tx_ids     <- family_txs$ID

  # tx → gene lookup
  tx_gene <- setNames(
    tx_parent[!is.na(tx_parent) & tx_parent %in% gene_ids],
    tx_ids
  )

  # -- exons --
  exons <- gff[gff$type == "exon"]
  ex_parent <- first_parent(exons$Parent)
  family_exons <- exons[!is.na(ex_parent) & ex_parent %in% tx_ids]

  exon_df <- data.frame(
    gene_id = tx_gene[first_parent(family_exons$Parent)],
    tx_id   = first_parent(family_exons$Parent),
    chrom   = as.character(seqnames(family_exons)),
    start   = start(family_exons),
    end     = end(family_exons),
    strand  = as.character(strand(family_exons)),
    stringsAsFactors = FALSE
  )
  exon_df <- exon_df[!is.na(exon_df$gene_id), ]
  message("  exons: ", nrow(exon_df), " across ", length(unique(exon_df$gene_id)), " genes")
  list(genes = family_genes, exon_df = exon_df)
}
```

- [ ] **Step 3: Call parse_family() for AMY and save**

```r
message("=== AMY family ===")
amy <- parse_family(gff, "amylase")
saveRDS(amy, file.path(DATA_DIR, "amy_family.rds"))
```

- [ ] **Step 4: Run and verify**

```bash
Rscript 01_parse_gff.R
```

Expected output includes:
```
=== AMY family ===
  genes found: 5
  exons: NNN across 5 genes
```

If `genes found` is 0, run this diagnostic to inspect description values:

```bash
Rscript -e "
  gff <- readRDS('data/gff_cache.rds')
  genes <- gff[gff\$type=='gene']
  desc  <- mcols(genes)\$description
  desc[is.na(desc)] <- ''
  print(head(sort(unique(desc[grepl('amy|amyl', desc, ignore.case=TRUE)])), 20))
"
```

- [ ] **Step 5: Commit**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
git add analysis/family_graphs/01_parse_gff.R
git commit -m "feat: parse_family() extracts AMY exon table from GFF"
```

---

## Task 3: Write and test build_variation_graph()

**Files:**
- Create: `analysis/family_graphs/02_build_graphs.R`
- Create: `analysis/family_graphs/tests/test_build_graphs.R`

Core algorithm: `disjoin()` all exons to get maximal non-overlapping segments = nodes. Each node gets `node_type = shared | copy_specific` based on how many distinct gene loci overlap it. Edges come from consecutive nodes within each transcript.

- [ ] **Step 1: Write failing tests first**

Write `tests/test_build_graphs.R`:

```r
library(testthat)
source("../02_build_graphs.R")

# Minimal fixture: 2 genes, one shared exon (100-200), one unique each
make_fixture <- function() {
  data.frame(
    gene_id = c("gA","gA","gB","gB"),
    tx_id   = c("tA","tA","tB","tB"),
    chrom   = rep("chr1", 4),
    start   = c(100L, 300L, 100L, 500L),
    end     = c(200L, 400L, 200L, 600L),
    strand  = rep("+", 4),
    stringsAsFactors = FALSE
  )
}

test_that("shared node has n_copies >= 2", {
  vg <- build_variation_graph(make_fixture())
  shared <- vg$nodes[vg$nodes$node_type == "shared", ]
  expect_true(nrow(shared) >= 1)
  expect_true(all(shared$n_copies >= 2))
})

test_that("copy-specific node belongs to exactly one gene", {
  vg <- build_variation_graph(make_fixture())
  specific <- vg$nodes[vg$nodes$node_type == "copy_specific", ]
  expect_true(nrow(specific) >= 2)   # gA has one, gB has one
  expect_true(all(specific$n_copies == 1))
})

test_that("each copy path traverses at least one copy-specific node", {
  vg <- build_variation_graph(make_fixture())
  for (gene in names(vg$paths)) {
    path_nodes <- vg$nodes[vg$nodes$node_id %in% vg$paths[[gene]], ]
    own_specific <- path_nodes[
      path_nodes$node_type == "copy_specific" &
      sapply(path_nodes$gene_set, function(gs) gene %in% gs), ]
    expect_true(nrow(own_specific) >= 1, info = paste("gene:", gene))
  }
})

test_that("is_multicopy_family() returns FALSE for single gene", {
  single <- make_fixture()[1:2, ]
  single$gene_id <- "gA"
  vg <- build_variation_graph(single)
  expect_false(is_multicopy_family(vg))
})

test_that("is_multicopy_family() returns FALSE for two unrelated genes (no shared exon)", {
  no_share <- data.frame(
    gene_id = c("gA","gB"),
    tx_id   = c("tA","tB"),
    chrom   = rep("chr1", 2),
    start   = c(100L, 1000L),
    end     = c(200L, 1100L),
    strand  = rep("+", 2),
    stringsAsFactors = FALSE
  )
  vg <- build_variation_graph(no_share)
  expect_false(is_multicopy_family(vg))
})
```

- [ ] **Step 2: Run tests — expect FAIL (functions not yet defined)**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle/analysis/family_graphs
Rscript -e "testthat::test_file('tests/test_build_graphs.R')"
```

Expected: errors about missing `build_variation_graph` and `is_multicopy_family`.

- [ ] **Step 3: Implement build_variation_graph()**

Write `02_build_graphs.R`:

```r
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(igraph)
  library(dplyr)
})

# Returns list(graph, nodes, edges, paths, exon_df)
build_variation_graph <- function(exon_df) {
  stopifnot(all(c("gene_id","tx_id","chrom","start","end","strand") %in% names(exon_df)))

  # -- 1. Build GRanges and disjoin into segments (nodes) --
  gr <- makeGRangesFromDataFrame(exon_df, keep.extra.columns = TRUE)
  segs <- disjoin(gr, ignore.strand = FALSE)
  mcols(segs)$node_id <- seq_along(segs)

  # -- 2. Map each segment to which gene_ids overlap it --
  hits <- findOverlaps(segs, gr)
  gene_per_exon <- exon_df$gene_id

  # Explicit loop: guarantees every segment index 1..N is represented
  gene_sets_vec <- vector("list", length(segs))
  for (k in seq_len(length(segs))) {
    exon_idx          <- subjectHits(hits)[queryHits(hits) == k]
    gene_sets_vec[[k]] <- unique(gene_per_exon[exon_idx])
  }

  mcols(segs)$gene_set  <- I(gene_sets_vec)
  mcols(segs)$n_copies  <- sapply(mcols(segs)$gene_set, length)
  mcols(segs)$node_type <- ifelse(mcols(segs)$n_copies >= 2L, "shared", "copy_specific")
  mcols(segs)$exon_len  <- width(segs)

  node_df <- as.data.frame(mcols(segs))
  node_df$chrom  <- as.character(seqnames(segs))
  node_df$start  <- start(segs)
  node_df$end    <- end(segs)
  node_df$strand <- as.character(strand(segs))

  # -- 3. Build edges: consecutive segments within each transcript --
  edge_list <- list()
  for (tx in unique(exon_df$tx_id)) {
    tx_rows  <- exon_df[exon_df$tx_id == tx, ]
    tx_gr    <- makeGRangesFromDataFrame(tx_rows)
    ov       <- findOverlaps(tx_gr, segs)
    seg_ids  <- sort(unique(mcols(segs)$node_id[subjectHits(ov)]))
    if (length(seg_ids) >= 2L) {
      for (i in seq_len(length(seg_ids) - 1L)) {
        edge_list <- c(edge_list, list(c(seg_ids[i], seg_ids[i + 1L])))
      }
    }
  }
  edge_df <- unique(do.call(rbind, edge_list))
  colnames(edge_df) <- c("from", "to")

  # -- 4. Build igraph --
  g <- graph_from_data_frame(
    d        = as.data.frame(edge_df),
    vertices = node_df[, c("node_id", "node_type", "n_copies", "exon_len",
                            "chrom", "start", "end", "strand")],
    directed = TRUE
  )

  # -- 5. Copy paths --
  paths <- list()
  for (gene in unique(exon_df$gene_id)) {
    gene_rows <- exon_df[exon_df$gene_id == gene, ]
    gene_gr   <- makeGRangesFromDataFrame(gene_rows)
    ov        <- findOverlaps(gene_gr, segs)
    paths[[gene]] <- sort(unique(mcols(segs)$node_id[subjectHits(ov)]))
  }

  list(graph = g, nodes = node_df, edges = as.data.frame(edge_df),
       paths = paths, exon_df = exon_df)
}

is_multicopy_family <- function(vg) {
  n_genes    <- length(unique(vg$exon_df$gene_id))
  has_shared <- any(vg$nodes$node_type == "shared")
  has_spec   <- any(vg$nodes$node_type == "copy_specific")
  n_genes >= 2L && has_shared && has_spec
}
```

- [ ] **Step 4: Run tests — expect PASS**

```bash
Rscript -e "testthat::test_file('tests/test_build_graphs.R')"
```

Expected:
```
[ FAIL 0 | WARN 0 | SKIP 0 | PASS 5 ]
```

- [ ] **Step 5: Build AMY graph and save**

Append to `02_build_graphs.R`:

```r
if (!interactive()) {
  DATA_DIR <- "data"
  amy <- readRDS(file.path(DATA_DIR, "amy_family.rds"))
  message("Building AMY variation graph...")
  amy_vg <- build_variation_graph(amy$exon_df)
  message("  nodes: ", nrow(amy_vg$nodes),
          "  shared: ", sum(amy_vg$nodes$node_type == "shared"),
          "  copy-specific: ", sum(amy_vg$nodes$node_type == "copy_specific"))
  message("  is_multicopy_family: ", is_multicopy_family(amy_vg))
  saveRDS(amy_vg, file.path(DATA_DIR, "amy_vg.rds"))
}
```

- [ ] **Step 6: Run and verify**

```bash
Rscript 02_build_graphs.R
```

Expected output:
```
Building AMY variation graph...
  nodes: NN  shared: NN  copy-specific: NN
  is_multicopy_family: TRUE
```

- [ ] **Step 7: Commit**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
git add analysis/family_graphs/02_build_graphs.R \
        analysis/family_graphs/tests/test_build_graphs.R
git commit -m "feat: build_variation_graph() with TDD — AMY graph passes all tests"
```

---

## Task 4: Figure 1 — The Problem

**Files:**
- Create: `analysis/family_graphs/03_figure_problem.R`

Three-pane figure. Pane A: genomic strip of AMY loci (exon blocks). Pane B: same strip + simulated multi-mapping arcs between shared exon positions across copies. Pane C: naive collapsed output (all arcs forced to locus 1). Output: `figures/fig1_problem.pdf` and `.png`.

Simulated arcs: for each pair of AMY copies, draw one arc per shared node between the midpoint of that node in copy i and copy j. Arc alpha = 0.4 to convey ambiguity.

- [ ] **Step 1: Write 03_figure_problem.R**

```r
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

DATA_DIR    <- "data"
FIGURES_DIR <- "figures"

amy    <- readRDS(file.path(DATA_DIR, "amy_family.rds"))
amy_vg <- readRDS(file.path(DATA_DIR, "amy_vg.rds"))

exon_df <- amy_vg$exon_df
nodes   <- amy_vg$nodes
paths   <- amy_vg$paths

# Order genes by genomic start position
gene_order <- exon_df %>%
  group_by(gene_id) %>%
  summarise(locus_start = min(start), .groups = "drop") %>%
  arrange(locus_start)

exon_df$gene_rank <- match(exon_df$gene_id, gene_order$gene_id)

# --- Pane A: genomic strip ---
pA <- ggplot(exon_df, aes(xmin = start, xmax = end,
                           ymin = gene_rank - 0.3, ymax = gene_rank + 0.3,
                           fill = gene_id)) +
  geom_rect(show.legend = FALSE) +
  geom_hline(aes(yintercept = gene_rank), color = "gray40", linewidth = 0.3) +
  scale_y_continuous(breaks = gene_order$gene_rank,
                     labels = gene_order$gene_id) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Genomic position (bp)", y = NULL, title = "A  Annotated loci") +
  theme_classic(base_size = 9) +
  theme(axis.text.y = element_text(size = 7))

# --- Pane B: multi-mapping arcs ---
# For each shared node, draw arcs between every pair of copies that traverse it
shared_nodes <- nodes[nodes$node_type == "shared", ]

arc_df <- do.call(rbind, lapply(seq_len(nrow(shared_nodes)), function(i) {
  nd      <- shared_nodes[i, ]
  mid_pos <- (nd$start + nd$end) / 2
  copies  <- nd$gene_set[[1]]
  if (length(copies) < 2) return(NULL)
  pairs   <- combn(copies, 2, simplify = FALSE)
  do.call(rbind, lapply(pairs, function(p) {
    r1 <- match(p[1], gene_order$gene_id)
    r2 <- match(p[2], gene_order$gene_id)
    data.frame(x = mid_pos, y_from = r1, y_to = r2,
               node_id = nd$node_id, stringsAsFactors = FALSE)
  }))
}))

pB <- pA +
  geom_curve(data = arc_df,
             aes(x = x, xend = x, y = y_from, yend = y_to),
             curvature = 0.5, alpha = 0.35, color = "firebrick",
             linewidth = 0.5, inherit.aes = FALSE) +
  labs(title = "B  Multi-mapping ambiguity")

# --- Pane C: naive collapse (all reads assigned to first locus) ---
# Visually: only locus 1 has exon blocks; others are empty lines
first_gene    <- gene_order$gene_id[1]
collapsed_df  <- exon_df
collapsed_df$gene_rank[collapsed_df$gene_id != first_gene] <- NA

pC <- ggplot(exon_df, aes(y = gene_rank)) +
  geom_hline(aes(yintercept = gene_rank), color = "gray40", linewidth = 0.3) +
  geom_rect(data = exon_df[exon_df$gene_id == first_gene, ],
            aes(xmin = start, xmax = end,
                ymin = gene_rank - 0.3, ymax = gene_rank + 0.3),
            fill = scales::hue_pal()(5)[1], inherit.aes = FALSE) +
  scale_y_continuous(breaks = gene_order$gene_rank,
                     labels = gene_order$gene_id) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Genomic position (bp)", y = NULL,
       title = "C  Naive: reads collapsed to one locus") +
  theme_classic(base_size = 9) +
  theme(axis.text.y = element_text(size = 7))

fig1 <- pA / pB / pC + plot_layout(heights = c(1, 1, 1))

ggsave(file.path(FIGURES_DIR, "fig1_problem.pdf"), fig1,
       width = 7, height = 8, units = "in")
ggsave(file.path(FIGURES_DIR, "fig1_problem.png"), fig1,
       width = 7, height = 8, units = "in", dpi = 300)
message("Figure 1 saved.")
```

- [ ] **Step 2: Run and inspect output**

```bash
Rscript 03_figure_problem.R
```

Expected: `figures/fig1_problem.pdf` and `figures/fig1_problem.png` created. Open PDF and verify three stacked panels with labeled gene tracks.

- [ ] **Step 3: Commit**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
git add analysis/family_graphs/03_figure_problem.R
git commit -m "feat: Figure 1 — multi-mapping ambiguity across AMY loci"
```

---

## Task 5: Figure 2 — The Method (variation graph + copy paths)

**Files:**
- Create: `analysis/family_graphs/04_figure_method.R`

Three-pane figure using ggraph. Pane A: full DAG colored by node_type. Pane B: same DAG with copy-path ribbons overlaid (one color per gene). Pane C: zoomed subgraph of one shared node plus its flanking copy-specific bubbles, with labels anchoring the formal definition.

- [ ] **Step 1: Write 04_figure_method.R**

```r
suppressPackageStartupMessages({
  library(ggraph)
  library(igraph)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

DATA_DIR    <- "data"
FIGURES_DIR <- "figures"

amy_vg <- readRDS(file.path(DATA_DIR, "amy_vg.rds"))
g      <- amy_vg$graph
nodes  <- amy_vg$nodes
paths  <- amy_vg$paths

# Assign copy colors
gene_ids   <- names(paths)
copy_cols  <- setNames(scales::hue_pal()(length(gene_ids)), gene_ids)

# --- Pane A: full DAG ---
pA <- ggraph(g, layout = "sugiyama") +
  geom_edge_link(arrow = arrow(length = unit(2, "mm"), type = "closed"),
                 end_cap = circle(3, "mm"),
                 color = "gray60", linewidth = 0.4) +
  geom_node_point(aes(color = node_type,
                      size  = exon_len)) +
  scale_color_manual(values = c(shared = "gray50", copy_specific = "steelblue"),
                     name = "Node type") +
  scale_size_continuous(range = c(2, 8), guide = "none") +
  labs(title = "A  Variation graph") +
  theme_graph(base_size = 9) +
  theme(legend.position = "bottom")

# --- Pane B: copy paths overlaid ---
# For each copy, add colored edges along its path
path_edges <- do.call(rbind, lapply(gene_ids, function(gid) {
  p      <- paths[[gid]]
  if (length(p) < 2) return(NULL)
  from_v <- p[-length(p)]
  to_v   <- p[-1]
  data.frame(from = from_v, to = to_v, gene_id = gid,
             stringsAsFactors = FALSE)
}))

# Mark graph edges that belong to each path
E(g)$copy_gene <- NA_character_
for (i in seq_len(nrow(path_edges))) {
  # igraph vertex names are character strings (node_id column coerced on graph_from_data_frame)
  eid <- get.edge.ids(g, c(as.character(path_edges$from[i]),
                            as.character(path_edges$to[i])))
  if (eid > 0) E(g)$copy_gene[eid] <- path_edges$gene_id[i]
}

pB <- ggraph(g, layout = "sugiyama") +
  geom_edge_link(aes(color = copy_gene, alpha = !is.na(copy_gene)),
                 arrow = arrow(length = unit(2, "mm"), type = "closed"),
                 end_cap = circle(3, "mm"), linewidth = 0.7) +
  scale_edge_color_manual(values = copy_cols, na.value = "gray85", name = "Copy") +
  scale_edge_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.2), guide = "none") +
  geom_node_point(aes(color = node_type, size = exon_len)) +
  scale_color_manual(values = c(shared = "gray50", copy_specific = "steelblue"),
                     name = "Node type") +
  scale_size_continuous(range = c(2, 8), guide = "none") +
  labs(title = "B  Copy paths") +
  theme_graph(base_size = 9) +
  theme(legend.position = "bottom")

# --- Pane C: zoomed bubble subgraph ---
# Find a shared node that has >=2 copy-specific successors — ideal bubble centre
shared_ids <- nodes$node_id[nodes$node_type == "shared"]
bubble_centre <- NA
for (sid in shared_ids) {
  succs <- neighbors(g, as.character(sid), mode = "out")
  succ_types <- nodes$node_type[nodes$node_id %in% as.integer(names(succs))]
  if (sum(succ_types == "copy_specific") >= 2) {
    bubble_centre <- sid
    break
  }
}

if (!is.na(bubble_centre)) {
  succs  <- as.integer(names(neighbors(g, as.character(bubble_centre), mode = "out")))
  preds  <- as.integer(names(neighbors(g, as.character(bubble_centre), mode = "in")))
  zoom_ids   <- unique(c(preds, bubble_centre, succs))
  zoom_graph <- induced_subgraph(g, vids = as.character(zoom_ids))

  # Determine if any edges are "isoform skip" (from pred directly to succ)
  E(zoom_graph)$edge_type <- "junction"
  skip_pairs <- expand.grid(from = preds, to = succs)
  for (i in seq_len(nrow(skip_pairs))) {
    eid <- get.edge.ids(zoom_graph,
                        c(as.character(skip_pairs$from[[i]]),
                          as.character(skip_pairs$to[[i]])))
    if (eid > 0) E(zoom_graph)$edge_type[eid] <- "isoform_skip"
  }

  pC <- ggraph(zoom_graph, layout = "sugiyama") +
    geom_edge_link(aes(linetype = edge_type),
                   arrow = arrow(length = unit(2, "mm"), type = "closed"),
                   end_cap = circle(4, "mm"),
                   color = "gray40", linewidth = 0.6) +
    scale_edge_linetype_manual(values = c(junction = "solid", isoform_skip = "dashed"),
                               name = NULL) +
    geom_node_point(aes(color = node_type), size = 6) +
    geom_node_label(aes(label = ifelse(node_type == "shared",
                                       "shared\nexon",
                                       "copy-specific\nexon")),
                    size = 2.5, label.padding = unit(0.15, "lines")) +
    scale_color_manual(values = c(shared = "gray50", copy_specific = "steelblue"),
                       guide = "none") +
    labs(title = "C  Bubble detail") +
    theme_graph(base_size = 9) +
    theme(legend.position = "bottom")
} else {
  pC <- ggplot() + annotate("text", x=0.5, y=0.5, label="No bubble found") +
    theme_void() + labs(title = "C  Bubble detail")
}

fig2 <- pA + pB + pC + plot_layout(ncol = 3)

ggsave(file.path(FIGURES_DIR, "fig2_method.pdf"), fig2,
       width = 12, height = 5, units = "in")
ggsave(file.path(FIGURES_DIR, "fig2_method.png"), fig2,
       width = 12, height = 5, units = "in", dpi = 300)
message("Figure 2 saved.")
```

- [ ] **Step 2: Run and inspect**

```bash
Rscript 04_figure_method.R
```

Expected: `figures/fig2_method.pdf` — three panels side by side. Pane C should show a central node with two colored branches departing from it and a dashed skip edge.

- [ ] **Step 3: Commit**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
git add analysis/family_graphs/04_figure_method.R
git commit -m "feat: Figure 2 — variation graph and copy-path visualization"
```

---

## Task 6: Parse GOLGA family and build its graph

**Files:**
- Modify: `analysis/family_graphs/01_parse_gff.R`
- Modify: `analysis/family_graphs/02_build_graphs.R`

GOLGA provides the positive example for Figure 3. Pattern: `golgin`.

- [ ] **Step 1: Add GOLGA parsing to 01_parse_gff.R**

Append after AMY block:

```r
message("=== GOLGA family ===")
golga <- parse_family(gff, "golgin")
saveRDS(golga, file.path(DATA_DIR, "golga_family.rds"))
```

- [ ] **Step 2: Run and verify gene count**

```bash
Rscript 01_parse_gff.R
```

Expected: `genes found: N` where N ≥ 4 (gorilla has multiple GOLGA6L and GOLGA8 paralogs).

- [ ] **Step 3: Add GOLGA graph build to 02_build_graphs.R**

Append before the closing `}` of the `if (!interactive())` block:

```r
  golga <- readRDS(file.path(DATA_DIR, "golga_family.rds"))
  message("Building GOLGA variation graph...")
  golga_vg <- build_variation_graph(golga$exon_df)
  message("  nodes: ", nrow(golga_vg$nodes),
          "  shared: ", sum(golga_vg$nodes$node_type == "shared"),
          "  copy-specific: ", sum(golga_vg$nodes$node_type == "copy_specific"))
  message("  is_multicopy_family: ", is_multicopy_family(golga_vg))
  saveRDS(golga_vg, file.path(DATA_DIR, "golga_vg.rds"))
```

- [ ] **Step 4: Run and verify**

```bash
Rscript 02_build_graphs.R
```

Expected: `is_multicopy_family: TRUE` for GOLGA.

- [ ] **Step 5: Commit**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
git add analysis/family_graphs/01_parse_gff.R \
        analysis/family_graphs/02_build_graphs.R
git commit -m "feat: parse and build GOLGA variation graph"
```

---

## Task 7: Figure 3 — Examples and Counterexamples

**Files:**
- Create: `analysis/family_graphs/05_figure_examples.R`

Pane A: GOLGA graph (positive example). Pane B: single-gene counterexample (build from any single AMY gene). Pane C: two-unrelated-gene counterexample (build from first AMY gene + first GOLGA gene, no shared exons).

- [ ] **Step 1: Write 05_figure_examples.R**

```r
suppressPackageStartupMessages({
  library(ggraph)
  library(igraph)
  library(ggplot2)
  library(patchwork)
})

source("02_build_graphs.R")

DATA_DIR    <- "data"
FIGURES_DIR <- "figures"

golga_vg <- readRDS(file.path(DATA_DIR, "golga_vg.rds"))
amy_vg   <- readRDS(file.path(DATA_DIR, "amy_vg.rds"))

graph_pane <- function(g, title, verdict_label, verdict_color) {
  ggraph(g, layout = "sugiyama") +
    geom_edge_link(arrow = arrow(length = unit(2, "mm"), type = "closed"),
                   end_cap = circle(3, "mm"),
                   color = "gray60", linewidth = 0.4) +
    geom_node_point(aes(color = node_type, size = exon_len)) +
    scale_color_manual(values = c(shared = "gray50", copy_specific = "steelblue"),
                       guide = "none") +
    scale_size_continuous(range = c(2, 7), guide = "none") +
    annotate("label", x = Inf, y = Inf, label = verdict_label,
             hjust = 1.1, vjust = 1.5, size = 3,
             fill = verdict_color, color = "white", fontface = "bold",
             label.padding = unit(0.3, "lines")) +
    labs(title = title) +
    theme_graph(base_size = 9)
}

# --- Pane A: GOLGA (positive example) ---
pA <- graph_pane(golga_vg$graph,
                 "A  GOLGA — multi-copy family",
                 "multi-copy family ✓", "#2166ac")

# --- Pane B: single-gene counterexample ---
first_gene    <- unique(amy_vg$exon_df$gene_id)[1]
single_exons  <- amy_vg$exon_df[amy_vg$exon_df$gene_id == first_gene, ]
single_vg     <- build_variation_graph(single_exons)

pB <- graph_pane(single_vg$graph,
                 "B  Single locus — isoforms only",
                 "isoforms only ✗", "#d6604d")

# --- Pane C: two unrelated genes (no shared exons) ---
# Use first AMY gene + first GOLGA gene — different chromosomes, no overlap
amy_gene1   <- unique(amy_vg$exon_df$gene_id)[1]
golga_gene1 <- unique(golga_vg$exon_df$gene_id)[1]
unrelated_exons <- rbind(
  amy_vg$exon_df[amy_vg$exon_df$gene_id == amy_gene1, ],
  golga_vg$exon_df[golga_vg$exon_df$gene_id == golga_gene1, ]
)
unrelated_vg <- build_variation_graph(unrelated_exons)

pC <- graph_pane(unrelated_vg$graph,
                 "C  Unrelated genes — no shared exons",
                 "no shared exon structure ✗", "#d6604d")

fig3 <- pA + pB + pC + plot_layout(ncol = 3)

ggsave(file.path(FIGURES_DIR, "fig3_examples.pdf"), fig3,
       width = 12, height = 5, units = "in")
ggsave(file.path(FIGURES_DIR, "fig3_examples.png"), fig3,
       width = 12, height = 5, units = "in", dpi = 300)
message("Figure 3 saved.")
```

Note: `build_variation_graph` must be sourced. Add at the top of `05_figure_examples.R`:

```r
source("02_build_graphs.R")
```

But `02_build_graphs.R` currently runs side effects in `if (!interactive())`. Refactor that block so it only runs when called as a script, not when sourced:

Wrap the bottom of `02_build_graphs.R` in:
```r
if (sys.nframe() == 0L) {
  # ... existing data-build block ...
}
```

Replace the existing `if (!interactive())` with `if (sys.nframe() == 0L)`.

- [ ] **Step 2: Refactor 02_build_graphs.R entry guard**

In `02_build_graphs.R`, replace:
```r
if (!interactive()) {
```
with:
```r
if (sys.nframe() == 0L) {
```

- [ ] **Step 3: Run and inspect**

```bash
Rscript 05_figure_examples.R
```

Expected: `figures/fig3_examples.pdf` — three panels. Pane A should have a rich DAG with bubbles. Pane B should be a simple linear or lightly branched graph. Pane C should show two disconnected components.

- [ ] **Step 4: Commit**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
git add analysis/family_graphs/05_figure_examples.R \
        analysis/family_graphs/02_build_graphs.R
git commit -m "feat: Figure 3 — examples and counterexamples of multi-copy gene families"
```

---

## Task 8: Scale to NPIP

**Files:**
- Modify: `analysis/family_graphs/01_parse_gff.R`
- Modify: `analysis/family_graphs/02_build_graphs.R`

- [ ] **Step 1: Add NPIP parsing**

Append to `01_parse_gff.R`:

```r
message("=== NPIP family ===")
npip <- parse_family(gff, "nuclear pore.*interact|NPIP")
saveRDS(npip, file.path(DATA_DIR, "npip_family.rds"))
```

- [ ] **Step 2: Run and check gene count**

```bash
Rscript 01_parse_gff.R
```

If `genes found: 0` for NPIP, run the description diagnostic:

```bash
Rscript -e "
  gff <- readRDS('data/gff_cache.rds')
  genes <- gff[gff\$type=='gene']
  desc  <- mcols(genes)\$description
  desc[is.na(desc)] <- ''
  print(head(sort(unique(desc[grepl('pore|NPIP|nuclear', desc, ignore.case=TRUE)])), 20))
"
```

Adjust regex in `parse_family()` call to match the actual description strings found.

- [ ] **Step 3: Add NPIP graph build**

Append to the `if (sys.nframe() == 0L)` block in `02_build_graphs.R`:

```r
  npip <- readRDS(file.path(DATA_DIR, "npip_family.rds"))
  if (length(unique(npip$exon_df$gene_id)) >= 2) {
    message("Building NPIP variation graph...")
    npip_vg <- build_variation_graph(npip$exon_df)
    message("  is_multicopy_family: ", is_multicopy_family(npip_vg))
    saveRDS(npip_vg, file.path(DATA_DIR, "npip_vg.rds"))
  } else {
    message("NPIP: fewer than 2 genes found, skipping graph build")
  }
```

- [ ] **Step 4: Run and verify**

```bash
Rscript 02_build_graphs.R
```

- [ ] **Step 5: Commit**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
git add analysis/family_graphs/01_parse_gff.R \
        analysis/family_graphs/02_build_graphs.R
git commit -m "feat: parse and build NPIP variation graph"
```

---

## Done

At this point the following exist:
- `data/amy_vg.rds`, `data/golga_vg.rds`, `data/npip_vg.rds` — variation graphs
- `figures/fig1_problem.pdf/png` — multi-mapping ambiguity
- `figures/fig2_method.pdf/png` — variation graph + copy paths
- `figures/fig3_examples.pdf/png` — examples & counterexamples
- `tests/test_build_graphs.R` — 5 passing unit tests covering the formal definition conditions

Phase 2 (BAM-derived read linkage, HMM-EM, novel copy discovery) begins after advisor review.
