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

# ──────────────────────────────────────────────────────────────────────────────
# Read-equivalence of gene copies (mathematical formulation)
# ──────────────────────────────────────────────────────────────────────────────
# Let G = (V, E) be the variation graph for a gene family F.
# For each gene copy c in F, let T(c) ⊆ Paths(G) be the set of directed paths
# through G that correspond to c's annotated transcripts (one path per isoform).
#
# Two copies c, c' ∈ F are READ-EQUIVALENT (written c ≈ c') iff T(c) = T(c').
#
# Properties:
#   reflexive  : T(c) = T(c)
#   symmetric  : T(c) = T(c') ⇒ T(c') = T(c)
#   transitive : T(c)=T(c') ∧ T(c')=T(c'') ⇒ T(c)=T(c'')
# Hence ≈ is an equivalence relation. The quotient F/≈ is the set of
# RESOLVABLE EQUIVALENCE CLASSES; each class is the minimum unit at which
# the family can be quantified using read evidence on G.
#
# Resolvability theorem. For any pair c ≈ c', no read whose evidence is a
# pattern of nodes V ⊆ V and junctions E ⊆ E can distinguish c from c':
# the read is consistent with c iff it embeds in some path in T(c) = T(c'),
# which is by hypothesis the same set of paths for c'. □
#
# Witness lemma. If c ≁ c', then ∃ a node v or edge e that appears in some
# path of T(c) but no path of T(c') (or vice versa). v (resp. e) is called a
# WITNESS for the distinction; any read that maps to v (resp. spans e) can
# discriminate c from c'.
#
# Implementation notes:
#   • The set T(c) is built from vg$transcript_paths when available (one
#     ordered node-id sequence per annotated isoform). Falls back to the
#     single canonical path vg$paths[[c]] only if transcript_paths is absent.
#   • Two paths are compared by exact sequence equality (DAG paths are
#     uniquely determined by node sequence, so set-of-paths equality coincides
#     with set-of-sequences equality).
#   • Sets are canonicalised by sorting their string representations before
#     comparison, so signature equality is order-independent within T(c).
family_equivalence_classes <- function(vg) {
  copy_ids <- names(vg$paths)

  # Build T(c) for each copy c: set of distinct transcript paths
  tx_paths_per_copy <- if (!is.null(vg$transcript_paths) &&
                            length(vg$transcript_paths) > 0L) {
    # Map tx_id -> gene_id
    tx_gene <- unique(vg$exon_df[, c("gene_id", "tx_id")])
    lapply(copy_ids, function(c) {
      gene_txs <- tx_gene$tx_id[tx_gene$gene_id == c]
      gene_txs <- gene_txs[gene_txs %in% names(vg$transcript_paths)]
      paths <- lapply(gene_txs, function(tx) vg$transcript_paths[[tx]])
      unique(paths)
    })
  } else {
    lapply(copy_ids, function(c) list(vg$paths[[c]]))
  }
  names(tx_paths_per_copy) <- copy_ids

  # Canonical signature: sorted unique transcript paths
  signature <- vapply(tx_paths_per_copy, function(tx_set) {
    s <- vapply(tx_set, function(p) paste(as.character(p), collapse = ","),
                character(1))
    paste(sort(unique(s)), collapse = " | ")
  }, character(1))

  unique_sigs <- unique(unname(signature))
  class_names <- paste0("class_", seq_along(unique_sigs))
  classes <- setNames(
    lapply(unique_sigs, function(s) copy_ids[signature == s]),
    class_names
  )

  # Representative isoform set per class (the first member's T(c))
  class_isoform_sets <- setNames(
    lapply(seq_along(unique_sigs), function(i) {
      rep_copy <- classes[[i]][1]
      tx_paths_per_copy[[rep_copy]]
    }),
    class_names
  )

  # For convenience: the union of nodes/edges used by each class
  class_node_set <- setNames(
    lapply(class_isoform_sets, function(paths) {
      sort(unique(unlist(paths)))
    }),
    class_names
  )
  class_edge_set <- setNames(
    lapply(class_isoform_sets, function(paths) {
      e <- do.call(rbind, lapply(paths, function(p) {
        if (length(p) < 2L) return(matrix(integer(0), 0, 2))
        cbind(p[-length(p)], p[-1L])
      }))
      if (is.null(e) || nrow(e) == 0L) return(matrix(integer(0), 0, 2))
      unique(e)
    }),
    class_names
  )

  list(
    n_copies        = length(copy_ids),
    n_classes       = length(unique_sigs),
    classes         = classes,
    class_isoforms  = class_isoform_sets,
    class_nodes     = class_node_set,
    class_edges     = class_edge_set,
    fully_resolved  = length(unique_sigs) == length(copy_ids)
  )
}

# ──────────────────────────────────────────────────────────────────────────────
# EXTENDED DEFINITION: multi-component multi-copy gene loci
# ──────────────────────────────────────────────────────────────────────────────
# Let L be a set of N copies at clustered genomic positions. The variation
# graph G(L) decomposes into K ≥ 1 weakly-connected components G_1,...,G_K,
# where n_k = # copies whose path lies entirely in G_k (so Σ n_k = N).
#
# Classification:
#   K = 1, n_1 ≥ 2  → COHERENT multi-copy family (passes is_multicopy_family)
#   K > 1, any n_k>1 → MIXED CLUSTER (some sub-families collapse, others don't)
#   K = N            → FULLY FRAGMENTED CLUSTER (every copy is its own
#                     component; sequence-level paralogy too diverged for
#                     ≥85% identity but copies remain genomically clustered)
#
# Key property: reads land on nodes ∈ V(G_k) for exactly one k; therefore
# read-assignment is component-local. Equivalence classes are computed
# independently within each G_k. The total set of resolvable units is
# ∐_{k=1..K} (F_k / ≈_k), the disjoint union of per-component quotients.
#
# Quantification: any read landing in G_k is assigned to a member of F_k / ≈_k
# (the equivalence class quotient for that component). Cross-component
# assignment never occurs by construction.
analyse_locus_cluster <- function(vg) {
  g <- vg$graph
  if (igraph::vcount(g) == 0L) {
    return(list(type = "empty", n_components = 0L, n_copies = 0L,
                n_total_classes = 0L, fully_resolved = TRUE,
                components = data.frame()))
  }

  comps <- igraph::components(g, mode = "weak")
  n_comp   <- comps$no
  n_copies <- length(vg$paths)
  node_to_comp <- setNames(as.integer(comps$membership), igraph::V(g)$name)

  # Each path lies entirely in one component (DAG paths can't cross components)
  copy_to_comp <- vapply(vg$paths, function(p) {
    if (length(p) == 0L) return(NA_integer_)
    as.integer(unname(node_to_comp[as.character(p[1])]))
  }, integer(1))

  # Per-component equivalence-class analysis
  comp_rows <- list()
  total_classes  <- 0L
  fully_resolved <- TRUE

  for (k in seq_len(n_comp)) {
    member_copies <- names(copy_to_comp)[copy_to_comp == k]
    if (length(member_copies) == 0L) next   # component is exon-only orphan

    sub_paths <- vg$paths[member_copies]
    sub_exon  <- vg$exon_df[vg$exon_df$gene_id %in% member_copies, ]
    sub_tx    <- if (!is.null(vg$transcript_paths)) {
      tx_in_comp <- unique(sub_exon$tx_id)
      vg$transcript_paths[tx_in_comp[tx_in_comp %in% names(vg$transcript_paths)]]
    } else NULL

    sub_vg <- list(paths = sub_paths, exon_df = sub_exon,
                   transcript_paths = sub_tx)
    eq <- family_equivalence_classes(sub_vg)

    total_classes <- total_classes + eq$n_classes
    if (!eq$fully_resolved) fully_resolved <- FALSE

    comp_rows[[length(comp_rows) + 1L]] <- data.frame(
      component       = k,
      n_nodes         = comps$csize[k],
      n_copies        = length(member_copies),
      n_classes       = eq$n_classes,
      fully_resolved  = eq$fully_resolved,
      stringsAsFactors = FALSE
    )
  }

  comp_df <- if (length(comp_rows) > 0L) do.call(rbind, comp_rows) else data.frame()

  type <- if (n_comp == 1L && n_copies >= 2L) {
    "coherent_family"
  } else if (n_comp > 1L && any(comp_df$n_copies >= 2L)) {
    "mixed_cluster"
  } else if (n_comp == n_copies && n_copies >= 2L) {
    "fully_fragmented_cluster"
  } else if (n_copies == 1L) {
    "singleton"
  } else {
    "empty"
  }

  list(
    type            = type,
    n_components    = n_comp,
    n_copies        = n_copies,
    n_total_classes = total_classes,
    fully_resolved  = fully_resolved,
    components      = comp_df
  )
}

# Broader predicate that accepts multi-component clusters as valid quantification
# targets. Use this when L is known to be a genomically clustered set of
# paralogs but may be too diverged for a single coherent VG.
is_multicopy_locus <- function(vg) {
  if (is.null(vg) || length(vg$paths) < 2L) return(FALSE)
  cl <- analyse_locus_cluster(vg)
  cl$type %in% c("coherent_family", "mixed_cluster", "fully_fragmented_cluster")
}

# ──────────────────────────────────────────────────────────────────────────────
# Witness identification: for every pair of distinct classes, return the
# discriminating nodes and edges (constructive proof of distinctness).
# A read mapping to a witness node — or spanning a witness edge — can
# discriminate the two classes.
# ──────────────────────────────────────────────────────────────────────────────
class_witnesses <- function(eq) {
  k <- eq$n_classes
  if (k < 2L) return(list())
  pairs <- list()
  for (i in seq_len(k - 1L)) {
    for (j in (i + 1L):k) {
      ni <- eq$class_nodes[[i]]; nj <- eq$class_nodes[[j]]
      ei <- eq$class_edges[[i]]; ej <- eq$class_edges[[j]]
      # node-level witnesses
      nodes_i_only <- setdiff(ni, nj)
      nodes_j_only <- setdiff(nj, ni)
      # edge-level witnesses (require row-wise set diff)
      ei_str <- if (nrow(ei) > 0L) apply(ei, 1, paste, collapse = ">") else character(0)
      ej_str <- if (nrow(ej) > 0L) apply(ej, 1, paste, collapse = ">") else character(0)
      edges_i_only <- setdiff(ei_str, ej_str)
      edges_j_only <- setdiff(ej_str, ei_str)
      pairs[[paste0(names(eq$classes)[i], "_vs_", names(eq$classes)[j])]] <- list(
        nodes_i_only = nodes_i_only,
        nodes_j_only = nodes_j_only,
        edges_i_only = edges_i_only,
        edges_j_only = edges_j_only,
        n_node_witnesses = length(nodes_i_only) + length(nodes_j_only),
        n_edge_witnesses = length(edges_i_only) + length(edges_j_only)
      )
    }
  }
  pairs
}

# Returns TRUE if the variation graph represents a coherent multi-copy family.
# Tightened criteria (post 2026-05-19 methodology review across AMY/RBMY/GOLGA/
# NPIP and 5 GOLGA sub-families):
#   1. At least 2 gene copies
#   2. Has shared exons (structural homology)
#   3. Has copy-specific exons (paths must be distinguishable)
#   4. Single weakly-connected component
#   5. >= 20% of nodes are shared
#   6. Max copies sharing a node >= ceiling(n_copies / 2)
#   7. At least 2 distinct paths through the VG (else copies are all
#      structurally identical and the family collapses to one effective unit)
#
# NOTE: Returning TRUE does NOT imply every copy is individually resolvable.
# Families like RBMY have 5 distinct paths for 13 copies — they're valid
# multi-copy families that must be quantified at the equivalence-class level.
# Use family_equivalence_classes() to check resolvability.
is_multicopy_family <- function(vg) {
  n_genes <- length(vg$paths)
  if (is.null(n_genes) || n_genes < 2L) return(FALSE)

  nd <- vg$nodes
  if (sum(nd$node_type == "shared") == 0L) return(FALSE)
  if (sum(nd$node_type == "copy_specific") == 0L) return(FALSE)

  g <- vg$graph
  if (igraph::vcount(g) > 0L &&
      igraph::components(g, mode = "weak")$no > 1L) return(FALSE)

  pct_shared <- sum(nd$node_type == "shared") / nrow(nd)
  if (pct_shared < 0.20) return(FALSE)

  max_share <- max(nd$n_copies)
  if (max_share < ceiling(n_genes / 2)) return(FALSE)

  eq <- family_equivalence_classes(vg)
  if (eq$n_classes < 2L) return(FALSE)

  TRUE
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
