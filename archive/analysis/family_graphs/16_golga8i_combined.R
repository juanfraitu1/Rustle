suppressPackageStartupMessages({
  library(ggplot2)
  library(ggraph)
  library(igraph)
  library(dplyr)
  library(scales)
  library(patchwork)
})

FIGURES_DIR <- "figures"
DATA_DIR    <- "data"

vg        <- readRDS(file.path(DATA_DIR, "golga_GOLGA8I_vg_seq.rds"))
intervals <- vg$intervals
nodes     <- vg$nodes
graph     <- vg$graph

# ── Gene metadata + colors ────────────────────────────────────────────────────
gene_meta <- intervals %>%
  group_by(gene_id) %>%
  summarise(strand = strand[1], chrom = chrom[1],
            min_pos = min(start), max_pos = max(end),
            span_bp = max(end) - min(start),
            .groups = "drop") %>%
  arrange(min_pos)
gene_meta$y_rank     <- seq_len(nrow(gene_meta))
gene_meta$gene_label <- sub("^gene-", "", gene_meta$gene_id)

copy_pal <- setNames(
  c("#1565C0","#B71C1C","#2E7D32","#6A1B9A","#EF6C00","#00838F","#5D4037"),
  gene_meta$gene_id
)
shared_col <- "#e6a817"

# Add node metadata to per-copy intervals
intervals <- merge(intervals, nodes[, c("node_id","node_type","n_copies")],
                    by = "node_id", all.x = TRUE)
intervals <- intervals %>%
  left_join(gene_meta[, c("gene_id","strand","min_pos","max_pos","y_rank")],
            by = c("gene_id","strand")) %>%
  mutate(
    rel_start = ifelse(strand == "+", start - min_pos, max_pos - end),
    rel_end   = ifelse(strand == "+", end - min_pos,   max_pos - start),
    fill_col  = ifelse(node_type == "shared", shared_col, copy_pal[gene_id]),
    exon_w    = rel_end - rel_start
  )

intron_df <- intervals %>%
  arrange(gene_id, rel_start) %>%
  group_by(gene_id) %>%
  mutate(prev_end = lag(rel_end)) %>%
  filter(!is.na(prev_end)) %>%
  ungroup()

# ═══════════════════════════════════════════════════════════════════════════════
# Top panel: gene structure with node-ID labels
# ═══════════════════════════════════════════════════════════════════════════════
exon_h <- 0.32
LABEL_MIN_W <- 250  # only label exons wider than this (bp) inside the block

# For narrow exons, place the node-ID below the gene track
intervals$label_inside <- intervals$exon_w >= LABEL_MIN_W
intervals$label_y      <- ifelse(intervals$label_inside,
                                   intervals$y_rank,
                                   intervals$y_rank + 0.55)

y_labels <- paste0(gene_meta$gene_label,
                   "  (", gene_meta$strand, ", ",
                   gene_meta$span_bp %/% 1000, "kb)")

p_top <- ggplot() +
  geom_segment(data = intron_df,
               aes(x = prev_end, xend = rel_start,
                   y = y_rank,   yend = y_rank),
               color = "gray45", linewidth = 0.5) +
  geom_hline(yintercept = gene_meta$y_rank, color = "gray88", linewidth = 0.3) +
  geom_rect(data = intervals,
            aes(xmin = rel_start, xmax = rel_end,
                ymin = y_rank - exon_h, ymax = y_rank + exon_h,
                fill = fill_col),
            color = "white", linewidth = 0.4) +
  scale_fill_identity() +
  # Node-ID labels inside wide exons
  geom_text(data = filter(intervals, label_inside),
            aes(x = (rel_start + rel_end) / 2, y = label_y,
                label = node_id),
            color = "white", fontface = "bold", size = 2.8) +
  # Node-ID labels above narrow exons (small spec exons especially)
  geom_text(data = filter(intervals, !label_inside),
            aes(x = (rel_start + rel_end) / 2, y = label_y,
                label = node_id, color = fill_col),
            fontface = "bold", size = 2.6) +
  scale_color_identity() +
  scale_y_continuous(breaks = gene_meta$y_rank, labels = y_labels,
                     trans  = "reverse",
                     expand = c(0.07, 0.5)) +
  scale_x_continuous(labels = scales::comma, expand = c(0.01, 200)) +
  labs(title    = "GOLGA8I  (7 copies):  gene structure with VG node IDs",
       subtitle = "Numbers = VG node IDs. Gold = shared exon, colored = copy-specific. Cross-reference with VG below.",
       x = "Position within gene (bp, 5' to 3')", y = NULL) +
  theme_classic(base_size = 10) +
  theme(plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(color = "gray40", size = 9.5),
        axis.text.y   = element_text(face = "italic", size = 8.5),
        axis.text.x   = element_text(size = 8))

# ═══════════════════════════════════════════════════════════════════════════════
# Bottom panel: VG topology with node IDs, sized by exon length (log)
# ═══════════════════════════════════════════════════════════════════════════════
# Assign per-node fill: shared = gold; copy-specific = colored by the owning copy
node_owner <- sapply(nodes$node_id, function(nid) {
  if (nodes$node_type[match(nid, nodes$node_id)] == "shared") return(NA_character_)
  # find which gene owns this node
  ow <- nodes$gene_set[[which(nodes$node_id == nid)]]
  if (length(ow) == 0L) return(NA_character_) else ow[1]
})
V(graph)$type     <- nodes$node_type
V(graph)$exon_len <- nodes$exon_len
V(graph)$n_copies <- nodes$n_copies
V(graph)$owner    <- node_owner
V(graph)$fillc    <- ifelse(nodes$node_type == "shared",
                             shared_col,
                             unname(copy_pal[node_owner]))
V(graph)$label_id <- as.character(nodes$node_id)

# Highlight the 5→{2,3,6,9,22,24}→7 bubble explicitly via an edge attribute
bubble_in  <- as.character(c(5,5,5,5,5,5))
bubble_mid <- as.character(c(2,3,6,9,22,24))
bubble_out <- as.character(c(7,7,7,7,7,7))
bubble_edges <- rbind(cbind(bubble_in, bubble_mid),
                      cbind(bubble_mid, bubble_out))
E(graph)$is_bubble <- FALSE
for (i in seq_len(nrow(bubble_edges))) {
  eid <- get_edge_ids(graph, c(bubble_edges[i,1], bubble_edges[i,2]))
  if (eid > 0L) E(graph)$is_bubble[eid] <- TRUE
}

p_bot <- ggraph(graph, layout = "sugiyama") +
  geom_edge_link(aes(color = is_bubble, linewidth = is_bubble),
                 arrow   = arrow(length = unit(1.7, "mm"), type = "closed"),
                 end_cap = circle(3, "mm")) +
  scale_edge_color_manual(values = c(`FALSE` = "gray60", `TRUE` = "#7B1FA2"),
                          guide = "none") +
  scale_edge_width_manual(values = c(`FALSE` = 0.4, `TRUE` = 1.0),
                          guide = "none") +
  geom_node_point(aes(fill = fillc, size = exon_len),
                  shape = 21, color = "white", stroke = 0.5) +
  scale_fill_identity() +
  scale_size_continuous(trans  = "log10",
                        range  = c(2.5, 11),
                        breaks = c(50, 200, 1000, 3000),
                        name   = "exon len (bp)") +
  geom_node_text(aes(label = label_id),
                 size = 2.6, fontface = "bold", color = "gray10") +
  # Bubble annotation
  annotate("label", x = 0, y = -3.5,
           label = "Purple thick edges = the 5-branch bubble (S5 -> {2,3,6,9,22,24} -> S7)\nEach of the 6 copies takes a different branch; each branch is 39 bp",
           fill = "#F3E5F5", color = "#4A148C", size = 3.0,
           label.padding = unit(0.3, "lines"), hjust = 0) +
  labs(title    = "VG topology:  same node IDs as above",
       subtitle = "Node size = log(exon length). Gold = shared (n_copies >= 2); each color = a copy-specific node owned by that copy.") +
  theme_graph(base_size = 10) +
  theme(plot.title       = element_text(face = "bold", size = 12),
        plot.subtitle    = element_text(color = "gray40", size = 9.5),
        legend.position  = "right",
        legend.title     = element_text(size = 8),
        legend.text      = element_text(size = 7))

# ═══════════════════════════════════════════════════════════════════════════════
# Combine
# ═══════════════════════════════════════════════════════════════════════════════
fig <- p_top / p_bot + plot_layout(heights = c(1.1, 1.2))

ggsave(file.path(FIGURES_DIR, "fig_golga8i_combined.pdf"), fig,
       width = 15, height = 12, units = "in")
ggsave(file.path(FIGURES_DIR, "fig_golga8i_combined.png"), fig,
       width = 15, height = 12, units = "in", dpi = 180)
message("Combined figure saved: figures/fig_golga8i_combined.png")
