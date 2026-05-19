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
    exon_idx           <- subjectHits(hits)[queryHits(hits) == k]
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
  if (length(edge_list) == 0L) {
    edge_df <- matrix(integer(0), nrow = 0L, ncol = 2L)
  } else {
    edge_df <- unique(do.call(rbind, edge_list))
  }
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

if (sys.nframe() == 0L) {
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
