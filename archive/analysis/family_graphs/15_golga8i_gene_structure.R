suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
})

FIGURES_DIR <- "figures"
DATA_DIR    <- "data"

vg <- readRDS(file.path(DATA_DIR, "golga_GOLGA8I_vg_seq.rds"))

# Per-gene intervals (collapsed isoforms) with their VG node assignment
intervals <- vg$intervals
nodes     <- vg$nodes

intervals <- merge(intervals,
                    nodes[, c("node_id", "node_type", "n_copies")],
                    by = "node_id", all.x = TRUE)

# Gene-level metadata (strand, span)
gene_meta <- intervals %>%
  group_by(gene_id) %>%
  summarise(strand   = strand[1],
            chrom    = chrom[1],
            min_pos  = min(start),
            max_pos  = max(end),
            span_bp  = max(end) - min(start),
            .groups  = "drop") %>%
  arrange(min_pos)
gene_meta$y_rank    <- seq_len(nrow(gene_meta))
gene_meta$gene_label <- sub("^gene-", "", gene_meta$gene_id)

# Per-copy colors (consistent with the decision-tree figure)
copy_pal <- setNames(
  c("#1565C0","#B71C1C","#2E7D32","#6A1B9A","#EF6C00","#00838F","#5D4037"),
  gene_meta$gene_id
)
shared_col <- "#e6a817"

# Project each exon to a 5'->3' relative coordinate, flipping negative-strand genes
intervals <- intervals %>%
  left_join(gene_meta[, c("gene_id","strand","min_pos","max_pos","y_rank")],
            by = c("gene_id","strand")) %>%
  mutate(
    rel_start = ifelse(strand == "+",
                       start - min_pos,
                       max_pos - end),
    rel_end   = ifelse(strand == "+",
                       end - min_pos,
                       max_pos - start),
    fill_col  = ifelse(node_type == "shared",
                       shared_col,
                       copy_pal[gene_id])
  )

# Introns: between consecutive exons within each gene (in 5'->3' order)
intron_df <- intervals %>%
  arrange(gene_id, rel_start) %>%
  group_by(gene_id) %>%
  mutate(prev_end = lag(rel_end)) %>%
  filter(!is.na(prev_end)) %>%
  ungroup()

exon_h <- 0.34

# Y-axis label: include strand indicator and number of exons
y_labels <- paste0(gene_meta$gene_label,
                   "  (", gene_meta$strand, ", ",
                   gene_meta$span_bp %/% 1000, "kb)")

# Compute summary numbers for the subtitle
n_shared_nodes <- sum(nodes$node_type == "shared")
n_spec_nodes   <- sum(nodes$node_type == "copy_specific")
n_full_shared  <- sum(nodes$n_copies == nrow(gene_meta))

p <- ggplot() +
  # Introns
  geom_segment(data = intron_df,
               aes(x = prev_end, xend = rel_start,
                   y = y_rank, yend = y_rank),
               color = "gray45", linewidth = 0.5) +
  # Horizontal guide line through each gene
  geom_hline(yintercept = gene_meta$y_rank,
             color = "gray85", linewidth = 0.3) +
  # Exon blocks
  geom_rect(data = intervals,
            aes(xmin = rel_start, xmax = rel_end,
                ymin = y_rank - exon_h, ymax = y_rank + exon_h,
                fill = fill_col),
            color = "white", linewidth = 0.4) +
  scale_fill_identity() +
  scale_y_continuous(breaks = gene_meta$y_rank,
                     labels = y_labels,
                     trans = "reverse",  # locus 1 at top
                     expand = c(0.05, 0.5)) +
  scale_x_continuous(labels = scales::comma) +
  labs(title    = paste0("GOLGA8I:  exon-intron structure of all ",
                          nrow(gene_meta), " copies"),
       subtitle = paste0(
         "Gold = shared exon (sequence-similar across >=2 copies).  ",
         "Other colors = copy-specific exons (one color per copy).\n",
         "VG summary: ", n_shared_nodes, " shared / ",
         n_spec_nodes, " copy-specific nodes;  ",
         n_full_shared, " nodes shared by ALL 7 copies (the backbone).  ",
         "Negative-strand genes flipped to 5'->3'."),
       x = "Position within gene (bp, 5' to 3')",
       y = NULL,
       caption = paste0("Genomic locations: all 7 copies on ",
                         gene_meta$chrom[1],
                         " between ", scales::comma(min(gene_meta$min_pos)),
                         "-", scales::comma(max(gene_meta$max_pos)), " bp")) +
  theme_classic(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "gray40", size = 10,
                                     lineheight = 1.3),
        plot.caption  = element_text(color = "gray45", size = 9, hjust = 0),
        axis.text.y   = element_text(face = "italic", size = 9.5),
        axis.text.x   = element_text(size = 9))

ggsave(file.path(FIGURES_DIR, "fig_golga8i_genes.pdf"), p,
       width = 13, height = 5.5, units = "in")
ggsave(file.path(FIGURES_DIR, "fig_golga8i_genes.png"), p,
       width = 13, height = 5.5, units = "in", dpi = 200)
message("Figure saved: figures/fig_golga8i_genes.png")
