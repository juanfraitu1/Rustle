suppressPackageStartupMessages({
  library(ggraph)
  library(igraph)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
})

source("02_build_graphs.R")

DATA_DIR    <- "data"
FIGURES_DIR <- "figures"

rbmy_vg  <- readRDS(file.path(DATA_DIR, "rbmy_vg_seq.rds"))
g        <- rbmy_vg$graph
nodes    <- rbmy_vg$nodes
paths    <- rbmy_vg$paths
gene_ids <- names(paths)

# --- Build representative subgraph ---
# Include genes with copy-specific nodes (show divergence) + 2 backbone-only copies
spec_counts <- sapply(gene_ids, function(gid) {
  p <- paths[[gid]]
  sum(nodes$node_type[match(p, nodes$node_id)] == "copy_specific", na.rm = TRUE)
})
diverging   <- names(spec_counts[spec_counts > 0])
backbone    <- names(spec_counts[spec_counts == 0])
rep_genes   <- unique(c(diverging, head(backbone, max(0L, 5L - length(diverging)))))
rep_node_ids <- unique(unlist(paths[rep_genes]))
rep_g        <- induced_subgraph(g, vids = as.character(rep_node_ids))

copy_cols  <- setNames(hue_pal()(length(rep_genes)), rep_genes)

# --- Pane A: variation graph DAG ---
pA <- ggraph(rep_g, layout = "sugiyama") +
  geom_edge_link(arrow = arrow(length = unit(2, "mm"), type = "closed"),
                 end_cap = circle(3, "mm"),
                 color = "gray60", linewidth = 0.5) +
  geom_node_point(aes(color = node_type, size = exon_len)) +
  scale_color_manual(values = c(shared = "#e6a817", copy_specific = "steelblue"),
                     name = NULL,
                     labels = c(shared = "shared exon", copy_specific = "copy-specific")) +
  scale_size_continuous(range = c(2, 6), guide = "none") +
  labs(title = "A  Variation graph (RBMY)") +
  theme_graph(base_size = 10) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8))

# --- Pane B: copy paths overlaid ---
path_edges <- do.call(rbind, lapply(rep_genes, function(gid) {
  p <- paths[[gid]]
  p <- p[p %in% rep_node_ids]
  if (length(p) < 2L) return(NULL)
  data.frame(from = p[-length(p)], to = p[-1], gene_id = gid,
             stringsAsFactors = FALSE)
}))

E(rep_g)$copy_gene <- NA_character_
if (!is.null(path_edges) && nrow(path_edges) > 0) {
  for (i in seq_len(nrow(path_edges))) {
    eid <- get_edge_ids(rep_g, c(as.character(path_edges$from[i]),
                                  as.character(path_edges$to[i])))
    if (eid > 0L) E(rep_g)$copy_gene[eid] <- path_edges$gene_id[i]
  }
}

gene_labels <- sub("^gene-LOC", "LOC", rep_genes)
named_cols  <- setNames(copy_cols, rep_genes)

pB <- ggraph(rep_g, layout = "sugiyama") +
  geom_edge_link(aes(color = copy_gene,
                     alpha = ifelse(is.na(copy_gene), 0.1, 0.85)),
                 arrow = arrow(length = unit(2, "mm"), type = "closed"),
                 end_cap = circle(3, "mm"), linewidth = 0.8) +
  scale_edge_color_manual(values = named_cols, na.value = "gray90",
                          name = "Copy", na.translate = FALSE,
                          labels = setNames(gene_labels, rep_genes)) +
  scale_edge_alpha_identity(guide = "none") +
  geom_node_point(aes(color = node_type, size = exon_len)) +
  scale_color_manual(values = c(shared = "#e6a817", copy_specific = "steelblue"),
                     guide = "none") +
  scale_size_continuous(range = c(2, 6), guide = "none") +
  labs(title = "B  Copy paths") +
  theme_graph(base_size = 10) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"))

# --- Pane C: bubble — shared node branching to copy-specific successors ---
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
  succs <- as.integer(names(neighbors(g, as.character(bubble_centre), mode = "out")))
  preds <- as.integer(names(neighbors(g, as.character(bubble_centre), mode = "in")))
  if (length(succs) > 4L) succs <- succs[seq_len(4L)]
  if (length(preds) > 2L) preds <- preds[seq_len(2L)]

  zoom_ids <- unique(c(preds, bubble_centre, succs))
  zoom_g   <- induced_subgraph(g, vids = as.character(zoom_ids))
  E(zoom_g)$edge_type <- "junction"

  pC <- ggraph(zoom_g, layout = "sugiyama") +
    geom_edge_link(aes(linetype = edge_type),
                   arrow = arrow(length = unit(2.5, "mm"), type = "closed"),
                   end_cap = circle(5, "mm"),
                   color = "gray35", linewidth = 0.7) +
    scale_edge_linetype_manual(values = c(junction = "solid"), name = NULL) +
    geom_node_point(aes(color = node_type), size = 7) +
    geom_node_label(aes(label = ifelse(node_type == "shared", "shared", "copy-specific")),
                    size = 2.8, label.padding = unit(0.18, "lines"),
                    nudge_y = 0.35, fill = "white", alpha = 0.85) +
    scale_color_manual(values = c(shared = "#e6a817", copy_specific = "steelblue"),
                       guide = "none") +
    labs(title = "C  Shared exon bubble") +
    coord_cartesian(clip = "off") +
    theme_graph(base_size = 10) +
    theme(legend.position = "none",
          plot.margin = margin(t = 20, r = 20, b = 10, l = 10))
} else {
  pC <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "No bubble found", size = 4, color = "gray40") +
    theme_void() + labs(title = "C  Shared exon bubble")
}

fig2 <- pA + pB + pC + plot_layout(ncol = 3)

ggsave(file.path(FIGURES_DIR, "fig2_method.pdf"), fig2,
       width = 15, height = 6, units = "in")
ggsave(file.path(FIGURES_DIR, "fig2_method.png"), fig2,
       width = 15, height = 6, units = "in", dpi = 300)
message("Figure 2 saved.")
