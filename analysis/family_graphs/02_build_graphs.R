suppressPackageStartupMessages({
  library(GenomicRanges)
  library(igraph)
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
  # Pre-collect all edges as two integer vectors, then cbind once
  from_vec <- integer(0)
  to_vec   <- integer(0)

  for (tx in unique(exon_df$tx_id)) {
    tx_rows   <- exon_df[exon_df$tx_id == tx, ]
    tx_gr     <- makeGRangesFromDataFrame(tx_rows)
    ov        <- findOverlaps(tx_gr, segs)
    seg_ids   <- sort(unique(mcols(segs)$node_id[subjectHits(ov)]))
    tx_strand <- unique(tx_rows$strand)
    if (length(tx_strand) == 1L && tx_strand == "-") seg_ids <- rev(seg_ids)
    if (length(seg_ids) >= 2L) {
      from_vec <- c(from_vec, seg_ids[-length(seg_ids)])
      to_vec   <- c(to_vec,   seg_ids[-1])
    }
  }

  if (length(from_vec) == 0L) {
    edge_df <- matrix(integer(0), nrow = 0L, ncol = 2L)
  } else {
    edge_df <- unique(cbind(from = from_vec, to = to_vec))
  }
  colnames(edge_df) <- c("from", "to")

  # -- 4. Build igraph --
  g <- graph_from_data_frame(
    d        = as.data.frame(edge_df),
    vertices = node_df[, c("node_id", "node_type", "n_copies", "exon_len",
                            "chrom", "start", "end", "strand")],
    directed = TRUE
  )

  # -- 5. Copy paths: use canonical (most-exons) transcript per gene --
  paths <- list()
  for (gene in unique(exon_df$gene_id)) {
    gene_txs <- unique(exon_df$tx_id[exon_df$gene_id == gene])
    # pick transcript with most exons
    tx_exon_counts <- sapply(gene_txs, function(tx) sum(exon_df$tx_id == tx))
    canonical_tx   <- gene_txs[which.max(tx_exon_counts)]

    tx_rows  <- exon_df[exon_df$tx_id == canonical_tx, ]
    tx_gr    <- makeGRangesFromDataFrame(tx_rows)
    ov       <- findOverlaps(tx_gr, segs)
    seg_ids  <- sort(unique(mcols(segs)$node_id[subjectHits(ov)]))
    tx_strand <- unique(tx_rows$strand)
    if (length(tx_strand) == 1L && tx_strand == "-") seg_ids <- rev(seg_ids)
    paths[[gene]] <- seg_ids
  }

  list(graph = g, nodes = node_df, edges = as.data.frame(edge_df),
       paths = paths, exon_df = exon_df)
}

# Returns TRUE if the family has both shared exons (structural homology) AND
# copy-specific exons (sequence divergence between paralogs).
# A family where all copies are 100% identical would return FALSE — treated as
# one locus. This is conservative but appropriate for the variation-graph
# formulation (copies must be distinguishable by their paths).
is_multicopy_family <- function(vg) {
  n_genes    <- length(unique(vg$exon_df$gene_id))
  has_shared <- any(vg$nodes$node_type == "shared")
  has_spec   <- any(vg$nodes$node_type == "copy_specific")
  n_genes >= 2L && has_shared && has_spec
}

build_variation_graph_seq <- function(exon_df, genome_fa,
                                       min_identity = 0.85,
                                       min_coverage  = 0.80) {
  suppressPackageStartupMessages({
    library(Rsamtools)
    library(Biostrings)
    library(GenomicRanges)
  })

  # 1. Unique genomic intervals per gene (collapse isoforms)
  intervals <- unique(exon_df[, c("gene_id", "chrom", "start", "end", "strand")])
  intervals$exon_uid <- seq_len(nrow(intervals))

  # 2. Extract sequences (needs indexed FASTA)
  fa_path <- genome_fa
  if (!file.exists(paste0(fa_path, ".fai"))) indexFa(fa_path)
  fa   <- FaFile(fa_path)
  open(fa)
  gr   <- GRanges(seqnames = intervals$chrom,
                  ranges   = IRanges(intervals$start, intervals$end),
                  strand   = intervals$strand)
  seqs <- scanFa(fa, gr)
  close(fa)
  names(seqs) <- paste0("e", intervals$exon_uid, "_g",
                        match(intervals$gene_id, unique(intervals$gene_id)))

  # 3. Write temp FASTA; run minimap2 all-vs-all (-X = skip self, -k8 for short seqs)
  tmp_fa  <- tempfile(fileext = ".fa")
  tmp_paf <- tempfile(fileext = ".paf")
  writeXStringSet(seqs, tmp_fa)
  system2("minimap2",
          args = c("-X", "-c", "--cs", "-k8", "-w1",
                   shQuote(tmp_fa), shQuote(tmp_fa)),
          stdout = tmp_paf, stderr = FALSE)

  # 4. Parse PAF; filter by identity + bidirectional coverage + cross-gene
  paf <- tryCatch(
    read.table(tmp_paf, header = FALSE, sep = "\t", fill = TRUE,
               col.names = c("qname","qlen","qstart","qend","ostrand",
                              "tname","tlen","tstart","tend","nmatch","alen","mapq",
                              paste0("v", 1:20))[seq_len(32)]),
    error = function(e) data.frame()
  )

  node_id_vec <- seq_len(nrow(intervals))   # union-find parent array
  find_root <- function(x) {
    while (node_id_vec[x] != x) {
      node_id_vec[x] <<- node_id_vec[node_id_vec[x]]
      x <- node_id_vec[x]
    }
    x
  }
  unite <- function(x, y) {
    rx <- find_root(x); ry <- find_root(y)
    if (rx != ry) node_id_vec[rx] <<- ry
  }

  if (nrow(paf) > 0 && ncol(paf) >= 12) {
    paf$identity <- paf$nmatch / pmax(paf$alen, 1L)
    paf$q_cov    <- (paf$qend - paf$qstart) / pmax(paf$qlen, 1L)
    paf$t_cov    <- (paf$tend - paf$tstart) / pmax(paf$tlen, 1L)
    paf$q_uid    <- as.integer(sub("e([0-9]+)_.*", "\\1", paf$qname))
    paf$t_uid    <- as.integer(sub("e([0-9]+)_.*", "\\1", paf$tname))
    paf$q_gidx   <- as.integer(sub(".*_g([0-9]+)", "\\1", paf$qname))
    paf$t_gidx   <- as.integer(sub(".*_g([0-9]+)", "\\1", paf$tname))

    good <- paf[!is.na(paf$q_uid) &
                paf$identity >= min_identity &
                pmin(paf$q_cov, paf$t_cov) >= min_coverage &
                paf$q_gidx != paf$t_gidx, ]

    for (i in seq_len(nrow(good))) unite(good$q_uid[i], good$t_uid[i])
  }

  # 5. Assign node IDs (canonical roots, re-mapped to 1..n)
  roots <- sapply(seq_len(nrow(intervals)), find_root)
  node_map <- match(roots, sort(unique(roots)))
  intervals$node_id <- node_map
  n_nodes <- max(node_map)

  # 6. Build node metadata
  node_genes <- lapply(seq_len(n_nodes), function(nid) {
    unique(intervals$gene_id[intervals$node_id == nid])
  })
  exon_len_mean <- sapply(seq_len(n_nodes), function(nid) {
    rows <- intervals[intervals$node_id == nid, ]
    mean(rows$end - rows$start + 1L)
  })
  node_df <- data.frame(
    node_id   = seq_len(n_nodes),
    n_copies  = sapply(node_genes, length),
    exon_len  = exon_len_mean,
    stringsAsFactors = FALSE
  )
  node_df$gene_set  <- I(node_genes)
  node_df$node_type <- ifelse(node_df$n_copies >= 2L, "shared", "copy_specific")

  # 7. Map each exon occurrence in exon_df to its node_id
  exon_to_node <- merge(
    exon_df,
    intervals[, c("gene_id", "chrom", "start", "end", "node_id")],
    by = c("gene_id", "chrom", "start", "end"),
    all.x = TRUE
  )

  # 8. Build edges: consecutive nodes within each transcript (strand-aware)
  from_vec <- integer(0); to_vec <- integer(0)
  for (tx in unique(exon_to_node$tx_id)) {
    rows     <- exon_to_node[exon_to_node$tx_id == tx, ]
    rows     <- rows[order(rows$start), ]
    seg_ids  <- unique(rows$node_id[!is.na(rows$node_id)])
    tx_strand <- unique(rows$strand)
    if (length(tx_strand) == 1L && tx_strand == "-") seg_ids <- rev(seg_ids)
    if (length(seg_ids) < 2L) next
    from_vec <- c(from_vec, seg_ids[-length(seg_ids)])
    to_vec   <- c(to_vec,   seg_ids[-1L])
  }
  edge_df <- unique(data.frame(from = from_vec, to = to_vec, stringsAsFactors = FALSE))

  # 9. Build igraph
  vertex_df <- data.frame(
    name      = as.character(node_df$node_id),
    node_type = node_df$node_type,
    n_copies  = node_df$n_copies,
    exon_len  = node_df$exon_len,
    stringsAsFactors = FALSE
  )
  if (nrow(edge_df) == 0L) {
    g <- igraph::make_empty_graph(n = nrow(vertex_df), directed = TRUE)
    igraph::V(g)$name      <- vertex_df$name
    igraph::V(g)$node_type <- vertex_df$node_type
    igraph::V(g)$n_copies  <- vertex_df$n_copies
    igraph::V(g)$exon_len  <- vertex_df$exon_len
  } else {
    g <- igraph::graph_from_data_frame(
      d        = as.data.frame(lapply(edge_df, as.character)),
      vertices = vertex_df,
      directed = TRUE
    )
  }

  # 10. Copy paths: canonical (most-exon) transcript per gene
  #     transcript_paths: every isoform path (for Phase 2 read assignment)
  paths            <- list()
  transcript_paths <- list()
  for (gene in unique(exon_df$gene_id)) {
    gene_txs       <- unique(exon_df$tx_id[exon_df$gene_id == gene])
    tx_exon_counts <- sapply(gene_txs, function(tx) sum(exon_df$tx_id == tx))
    canonical_tx   <- gene_txs[which.max(tx_exon_counts)]

    for (tx in gene_txs) {
      rows      <- exon_to_node[exon_to_node$tx_id == tx, ]
      rows      <- rows[order(rows$start), ]
      seg_ids   <- unique(rows$node_id[!is.na(rows$node_id)])
      tx_strand <- unique(rows$strand)
      if (length(tx_strand) == 1L && tx_strand == "-") seg_ids <- rev(seg_ids)
      transcript_paths[[tx]] <- seg_ids
      if (tx == canonical_tx) paths[[gene]] <- seg_ids
    }
  }

  list(graph = g, nodes = node_df, edges = edge_df,
       paths = paths, transcript_paths = transcript_paths,
       exon_df = exon_df, intervals = intervals)
}

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

  golga <- readRDS(file.path(DATA_DIR, "golga_family.rds"))
  message("Building GOLGA variation graph...")
  golga_vg <- build_variation_graph(golga$exon_df)
  message("  nodes: ", nrow(golga_vg$nodes),
          "  shared: ", sum(golga_vg$nodes$node_type == "shared"),
          "  copy-specific: ", sum(golga_vg$nodes$node_type == "copy_specific"))
  message("  is_multicopy_family: ", is_multicopy_family(golga_vg))
  saveRDS(golga_vg, file.path(DATA_DIR, "golga_vg.rds"))

  npip <- readRDS(file.path(DATA_DIR, "npip_family.rds"))
  if (length(unique(npip$exon_df$gene_id)) >= 2L) {
    message("Building NPIP variation graph...")
    npip_vg <- build_variation_graph(npip$exon_df)
    message("  nodes: ", nrow(npip_vg$nodes),
            "  shared: ", sum(npip_vg$nodes$node_type == "shared"),
            "  copy-specific: ", sum(npip_vg$nodes$node_type == "copy_specific"))
    message("  is_multicopy_family: ", is_multicopy_family(npip_vg))
    saveRDS(npip_vg, file.path(DATA_DIR, "npip_vg.rds"))
  } else {
    message("NPIP: fewer than 2 genes found, skipping graph build")
  }

  GENOME_FA <- "/mnt/c/Users/jfris/Desktop/GGO.fasta"

  message("Building RBMY variation graph (sequence-based)...")
  rbmy_family <- readRDS(file.path(DATA_DIR, "rbmy_family.rds"))
  rbmy_vg_seq <- build_variation_graph_seq(rbmy_family$exon_df, GENOME_FA)
  cat(sprintf("  nodes: %d  shared: %d  is_multicopy: %s\n",
    nrow(rbmy_vg_seq$nodes),
    sum(rbmy_vg_seq$nodes$node_type == "shared"),
    is_multicopy_family(rbmy_vg_seq)))
  saveRDS(rbmy_vg_seq, file.path(DATA_DIR, "rbmy_vg_seq.rds"))

  message("Building GOLGA variation graph (sequence-based)...")
  golga_vg_seq <- build_variation_graph_seq(golga$exon_df, GENOME_FA)
  cat(sprintf("  nodes: %d  shared: %d  is_multicopy: %s\n",
    nrow(golga_vg_seq$nodes),
    sum(golga_vg_seq$nodes$node_type == "shared"),
    is_multicopy_family(golga_vg_seq)))
  saveRDS(golga_vg_seq, file.path(DATA_DIR, "golga_vg_seq.rds"))
}
