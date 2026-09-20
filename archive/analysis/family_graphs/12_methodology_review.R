suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
  library(ggplot2)
  library(ggraph)
  library(patchwork)
})

source("02_build_graphs.R")

DATA_DIR    <- "data"
FIGURES_DIR <- "figures"
GENOME_FA   <- "/mnt/c/Users/jfris/Desktop/GGO.fasta"

families <- c("amy", "golga", "npip", "rbmy")

# ── Build seq-VG for any family that doesn't have one yet ─────────────────────
ensure_seq_vg <- function(nm) {
  out <- file.path(DATA_DIR, paste0(nm, "_vg_seq.rds"))
  if (file.exists(out)) return(readRDS(out))
  fam <- readRDS(file.path(DATA_DIR, paste0(nm, "_family.rds")))
  if (length(unique(fam$exon_df$gene_id)) < 2L) {
    message("Skipping ", nm, ": <2 genes")
    return(NULL)
  }
  message("Building seq VG for ", nm, " ...")
  vg <- build_variation_graph_seq(fam$exon_df, GENOME_FA)
  saveRDS(vg, out)
  vg
}

vgs <- setNames(lapply(families, ensure_seq_vg), families)

# ── Structural metrics + methodology checks ───────────────────────────────────
count_bubbles <- function(vg) {
  g  <- vg$graph
  nd <- vg$nodes
  if (vcount(g) == 0L) return(0L)
  # bubble = shared node with >=2 out-neighbours that are copy-specific
  # and those copy-specific nodes share a downstream shared node
  out_list <- ego(g, order = 1, mode = "out")
  n_bubble <- 0L
  for (i in seq_along(out_list)) {
    if (nd$node_type[i] != "shared") next
    succs <- setdiff(as.integer(out_list[[i]]), i)
    spec_succs <- succs[nd$node_type[succs] == "copy_specific"]
    if (length(spec_succs) < 2L) next
    # check if at least 2 of the spec successors converge to a common shared node
    downs <- lapply(spec_succs, function(s) {
      d <- setdiff(as.integer(ego(g, order = 1, mode = "out")[[s]]), s)
      d[nd$node_type[d] == "shared"]
    })
    common <- Reduce(intersect, downs)
    if (length(common) >= 1L) n_bubble <- n_bubble + 1L
  }
  n_bubble
}

metric_row <- function(nm, vg) {
  if (is.null(vg)) {
    return(data.frame(family=toupper(nm), n_copies=NA, status="<2 genes"))
  }
  g  <- vg$graph
  nd <- vg$nodes
  n_genes  <- length(vg$paths)
  n_total  <- nrow(nd)
  n_shared <- sum(nd$node_type == "shared")
  n_spec   <- sum(nd$node_type == "copy_specific")
  pct_spec <- round(100 * n_spec / n_total, 1)
  full_sh  <- sum(nd$n_copies == n_genes)
  max_sh   <- max(nd$n_copies)
  median_sh<- median(nd$n_copies[nd$node_type == "shared"])
  n_edges  <- ecount(g)
  n_comp   <- if (vcount(g) > 0L) igraph::components(g, mode = "weak")$no else NA_integer_
  n_bub    <- count_bubbles(vg)
  ok       <- is_multicopy_family(vg)
  data.frame(
    family            = toupper(nm),
    n_copies          = n_genes,
    n_nodes           = n_total,
    n_shared          = n_shared,
    n_spec            = n_spec,
    pct_spec          = pct_spec,
    full_shared       = full_sh,
    max_copies_shared = max_sh,
    median_n_shared   = median_sh,
    n_edges           = n_edges,
    n_components      = n_comp,
    n_bubbles         = n_bub,
    is_multicopy      = ok,
    stringsAsFactors  = FALSE
  )
}

stats <- do.call(rbind, lapply(families, function(n) metric_row(n, vgs[[n]])))

cat("\n========== METHODOLOGY REVIEW ACROSS FAMILIES ==========\n\n")
print(stats, row.names = FALSE)
cat("\n")

write.table(stats, file.path(FIGURES_DIR, "methodology_review.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Stats written to ", file.path(FIGURES_DIR, "methodology_review.tsv"))

# ── Distribution of n_copies per shared node (does the family really share?) ──
share_dist_df <- do.call(rbind, lapply(families, function(nm) {
  vg <- vgs[[nm]]
  if (is.null(vg)) return(NULL)
  data.frame(family = toupper(nm),
             n_copies_sharing = vg$nodes$n_copies[vg$nodes$node_type == "shared"],
             stringsAsFactors = FALSE)
}))

p_share <- ggplot(share_dist_df, aes(x = n_copies_sharing)) +
  geom_histogram(binwidth = 1, fill = "#1565C0", color = "white") +
  facet_wrap(~ family, scales = "free", ncol = 2) +
  labs(title = "Shared-node depth: how many copies share each shared node?",
       subtitle = "x = number of copies a shared node appears in (>=2 by definition)",
       x = "copies sharing this node", y = "count of shared nodes") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"))

# ── Topology thumbnails: one per family ───────────────────────────────────────
mk_topology <- function(nm, vg) {
  if (is.null(vg) || vcount(vg$graph) == 0L) return(NULL)
  g  <- vg$graph
  nd <- vg$nodes
  # ggraph with sugiyama layout
  ggraph(g, layout = "sugiyama") +
    geom_edge_link(arrow = arrow(length = unit(1.5,"mm"), type="closed"),
                   end_cap = circle(2,"mm"),
                   color = "gray60", linewidth = 0.3) +
    geom_node_point(aes(color = node_type, size = n_copies)) +
    scale_color_manual(values = c(shared = "#e6a817",
                                   copy_specific = "#1565C0"),
                       guide = "none") +
    scale_size_continuous(range = c(1, 4), guide = "none") +
    labs(title = sprintf("%s  (%d copies, %d nodes, %d bubbles)",
                          toupper(nm),
                          length(vg$paths),
                          nrow(nd),
                          count_bubbles(vg))) +
    theme_graph(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 10))
}

topo_panels <- lapply(families, function(n) mk_topology(n, vgs[[n]]))
topo_panels <- topo_panels[!sapply(topo_panels, is.null)]

p_topo <- wrap_plots(topo_panels, ncol = 2)

# ── Save figures ──────────────────────────────────────────────────────────────
ggsave(file.path(FIGURES_DIR, "review_share_dist.png"), p_share,
       width = 10, height = 7, units = "in", dpi = 180)
ggsave(file.path(FIGURES_DIR, "review_topologies.png"), p_topo,
       width = 12, height = 10, units = "in", dpi = 180)
message("Review figures saved.")
