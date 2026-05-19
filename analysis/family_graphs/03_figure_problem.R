suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
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

message("Arc pairs: ", if (is.null(arc_df)) 0 else nrow(arc_df))

pB <- pA +
  geom_curve(data = arc_df,
             aes(x = x, xend = x, y = y_from, yend = y_to),
             curvature = 0.5, alpha = 0.35, color = "firebrick",
             linewidth = 0.5, inherit.aes = FALSE) +
  labs(title = "B  Multi-mapping ambiguity")

# --- Pane C: naive collapse (all reads assigned to first locus) ---
first_gene <- gene_order$gene_id[1]

pC <- ggplot(exon_df, aes(y = gene_rank)) +
  geom_hline(aes(yintercept = gene_rank), color = "gray40", linewidth = 0.3) +
  geom_rect(data = exon_df[exon_df$gene_id == first_gene, ],
            aes(xmin = start, xmax = end,
                ymin = gene_rank - 0.3, ymax = gene_rank + 0.3),
            fill = scales::hue_pal()(length(unique(exon_df$gene_id)))[1],
            inherit.aes = FALSE) +
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
