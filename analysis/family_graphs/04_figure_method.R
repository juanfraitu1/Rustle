suppressPackageStartupMessages({
  library(ggraph)
  library(igraph)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
})

DATA_DIR    <- "data"
FIGURES_DIR <- "figures"

amy_vg <- readRDS(file.path(DATA_DIR, "amy_vg.rds"))
g      <- amy_vg$graph
nodes  <- amy_vg$nodes
paths  <- amy_vg$paths

# Assign copy colors (one per gene)
gene_ids  <- names(paths)
copy_cols <- setNames(hue_pal()(length(gene_ids)), gene_ids)

# --- Pane A: full DAG ---
pA <- ggraph(g, layout = "sugiyama") +
  geom_edge_link(arrow = arrow(length = unit(2, "mm"), type = "closed"),
                 end_cap = circle(3, "mm"),
                 color = "gray60", linewidth = 0.4) +
  geom_node_point(aes(color = node_type, size = exon_len)) +
  scale_color_manual(values = c(shared = "gray50", copy_specific = "steelblue"),
                     name = "Node type") +
  scale_size_continuous(range = c(2, 8), guide = "none") +
  labs(title = "A  Variation graph") +
  theme_graph(base_size = 9) +
  theme(legend.position = "bottom")

# --- Pane B: copy paths overlaid ---
path_edges <- do.call(rbind, lapply(gene_ids, function(gid) {
  p <- paths[[gid]]
  if (length(p) < 2) return(NULL)
  data.frame(from = p[-length(p)], to = p[-1], gene_id = gid,
             stringsAsFactors = FALSE)
}))

# Mark each graph edge with which copy's path uses it
E(g)$copy_gene <- NA_character_
if (!is.null(path_edges) && nrow(path_edges) > 0) {
  for (i in seq_len(nrow(path_edges))) {
    eid <- get.edge.ids(g, c(as.character(path_edges$from[i]),
                              as.character(path_edges$to[i])))
    if (eid > 0L) E(g)$copy_gene[eid] <- path_edges$gene_id[i]
  }
}

pB <- ggraph(g, layout = "sugiyama") +
  geom_edge_link(aes(color = copy_gene,
                     alpha = ifelse(is.na(copy_gene), 0.15, 1.0)),
                 arrow = arrow(length = unit(2, "mm"), type = "closed"),
                 end_cap = circle(3, "mm"), linewidth = 0.7) +
  scale_edge_color_manual(values = copy_cols, na.value = "gray85",
                          name = "Copy", na.translate = FALSE) +
  scale_edge_alpha_identity(guide = "none") +
  geom_node_point(aes(color = node_type, size = exon_len)) +
  scale_color_manual(values = c(shared = "gray50", copy_specific = "steelblue"),
                     name = "Node type") +
  scale_size_continuous(range = c(2, 8), guide = "none") +
  labs(title = "B  Copy paths") +
  theme_graph(base_size = 9) +
  theme(legend.position = "bottom")

# --- Pane C: zoomed bubble subgraph ---
# Find shared node with >=2 copy-specific successors
shared_ids    <- nodes$node_id[nodes$node_type == "shared"]
bubble_centre <- NA_integer_

for (sid in shared_ids) {
  succs      <- as.integer(names(neighbors(g, as.character(sid), mode = "out")))
  succ_types <- nodes$node_type[nodes$node_id %in% succs]
  if (sum(succ_types == "copy_specific") >= 2L) {
    bubble_centre <- sid
    break
  }
}

if (!is.na(bubble_centre)) {
  succs    <- as.integer(names(neighbors(g, as.character(bubble_centre), mode = "out")))
  preds    <- as.integer(names(neighbors(g, as.character(bubble_centre), mode = "in")))

  # If too many successors, limit to 4 for readability
  if (length(succs) > 4L) {
    succs <- succs[seq_len(4L)]
  }

  zoom_ids <- unique(c(preds, bubble_centre, succs))
  zoom_g   <- induced_subgraph(g, vids = as.character(zoom_ids))

  E(zoom_g)$edge_type <- "junction"
  skip_pairs <- expand.grid(from = as.character(preds), to = as.character(succs),
                             stringsAsFactors = FALSE)
  for (i in seq_len(nrow(skip_pairs))) {
    eid <- get.edge.ids(zoom_g, c(skip_pairs$from[[i]], skip_pairs$to[[i]]))
    if (eid > 0L) E(zoom_g)$edge_type[eid] <- "isoform_skip"
  }

  pC <- ggraph(zoom_g, layout = "sugiyama") +
    geom_edge_link(aes(linetype = edge_type),
                   arrow = arrow(length = unit(2, "mm"), type = "closed"),
                   end_cap = circle(4, "mm"),
                   color = "gray40", linewidth = 0.6) +
    scale_edge_linetype_manual(
      values = c(junction = "solid", isoform_skip = "dashed"),
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
  # Fallback: no bubble found with 2+ copy-specific successors
  # Show the shared node(s) and their immediate neighbors instead
  if (length(shared_ids) > 0L) {
    sid      <- shared_ids[1]
    nb_ids   <- as.integer(names(neighbors(g, as.character(sid), mode = "all")))
    zoom_ids <- unique(c(sid, nb_ids))
    zoom_g   <- induced_subgraph(g, vids = as.character(zoom_ids))
    E(zoom_g)$edge_type <- "junction"

    pC <- ggraph(zoom_g, layout = "sugiyama") +
      geom_edge_link(aes(linetype = edge_type),
                     arrow = arrow(length = unit(2, "mm"), type = "closed"),
                     end_cap = circle(4, "mm"), color = "gray40", linewidth = 0.6) +
      scale_edge_linetype_manual(values = c(junction = "solid"), name = NULL) +
      geom_node_point(aes(color = node_type), size = 6) +
      geom_node_label(aes(label = ifelse(node_type == "shared",
                                         "shared\nexon", "copy-specific\nexon")),
                      size = 2.5, label.padding = unit(0.15, "lines")) +
      scale_color_manual(values = c(shared = "gray50", copy_specific = "steelblue"),
                         guide = "none") +
      labs(title = "C  Shared node neighborhood") +
      theme_graph(base_size = 9) +
      theme(legend.position = "bottom")
  } else {
    pC <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No shared node found") +
      theme_void() + labs(title = "C  (no shared node)")
  }
}

fig2 <- pA + pB + pC + plot_layout(ncol = 3)

ggsave(file.path(FIGURES_DIR, "fig2_method.pdf"), fig2,
       width = 12, height = 5, units = "in")
ggsave(file.path(FIGURES_DIR, "fig2_method.png"), fig2,
       width = 12, height = 5, units = "in", dpi = 300)
message("Figure 2 saved.")
