suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
})

source("02_build_graphs.R")

DATA_DIR    <- "data"
FIGURES_DIR <- "figures"

rbmy_vg <- readRDS(file.path(DATA_DIR, "rbmy_vg_seq.rds"))
exon_df  <- rbmy_vg$exon_df
nodes    <- rbmy_vg$nodes
intervals <- rbmy_vg$intervals   # gene_id, chrom, start, end, node_id

# --- Normalize coordinates to gene-relative position ---
gene_offsets <- exon_df %>%
  group_by(gene_id) %>%
  summarise(offset = min(start), .groups = "drop")

exon_df <- exon_df %>%
  left_join(gene_offsets, by = "gene_id") %>%
  mutate(rel_start = start - offset,
         rel_end   = end   - offset)

# Order genes by chromosome then genomic start
gene_order <- exon_df %>%
  group_by(gene_id) %>%
  summarise(chrom = chrom[1], locus_start = min(start), .groups = "drop") %>%
  arrange(chrom, locus_start) %>%
  mutate(gene_rank  = row_number(),
         gene_label = sub("^gene-", "", gene_id))

exon_df <- exon_df %>%
  left_join(gene_order[, c("gene_id", "gene_rank", "gene_label")], by = "gene_id")

# --- Arc data: for each shared node, draw arc between each pair of copies ---
# Use per-gene exon coordinates from intervals (sequence-based shared nodes)
shared_nodes <- nodes[nodes$node_type == "shared", ]

arc_df <- do.call(rbind, lapply(seq_len(nrow(shared_nodes)), function(i) {
  nd     <- shared_nodes[i, ]
  nid    <- nd$node_id
  copies <- nd$gene_set[[1]]
  if (length(copies) < 2L) return(NULL)

  # Get per-gene midpoints for this shared node
  node_intervals <- intervals[intervals$node_id == nid, ]
  rel_mids <- sapply(copies, function(cid) {
    row <- node_intervals[node_intervals$gene_id == cid, ]
    if (nrow(row) == 0L) return(NA_real_)
    off <- gene_offsets$offset[gene_offsets$gene_id == cid]
    ((row$start[1] + row$end[1]) / 2) - off
  })
  copies <- copies[!is.na(rel_mids)]
  rel_mids <- rel_mids[!is.na(rel_mids)]
  if (length(copies) < 2L) return(NULL)

  pairs <- combn(seq_along(copies), 2, simplify = FALSE)
  do.call(rbind, lapply(pairs, function(p) {
    i1 <- p[1]; i2 <- p[2]
    data.frame(
      x    = rel_mids[i1], xend = rel_mids[i2],
      y    = gene_order$gene_rank[gene_order$gene_id == copies[i1]],
      yend = gene_order$gene_rank[gene_order$gene_id == copies[i2]],
      stringsAsFactors = FALSE
    )
  }))
}))

# Limit arcs to a readable subset (sample if too many)
if (!is.null(arc_df) && nrow(arc_df) > 60) {
  set.seed(42)
  arc_df <- arc_df[sample(nrow(arc_df), 60), ]
}

first_gene <- gene_order$gene_id[1]

# --- Pane A: gene models (normalized) ---
pA <- ggplot(exon_df,
             aes(xmin = rel_start, xmax = rel_end,
                 ymin = gene_rank - 0.35, ymax = gene_rank + 0.35,
                 fill = gene_id)) +
  geom_hline(yintercept = unique(exon_df$gene_rank), color = "gray70", linewidth = 0.25) +
  geom_rect(show.legend = FALSE) +
  scale_y_continuous(breaks = gene_order$gene_rank,
                     labels = gene_order$gene_label,
                     expand = c(0.03, 0)) +
  scale_x_continuous(labels = comma, expand = c(0.02, 0)) +
  labs(x = "Position within gene (bp)", y = NULL,
       title = "A  RBMY — 13 copies on chr Y") +
  theme_classic(base_size = 9) +
  theme(axis.text.y = element_text(size = 6, hjust = 1))

# --- Pane B: arcs connecting sequence-similar exons across copies ---
pB <- pA +
  labs(title = "B  Multi-mapping ambiguity") +
  { if (!is.null(arc_df) && nrow(arc_df) > 0)
      geom_curve(data = arc_df,
                 aes(x = x, xend = xend, y = y, yend = yend),
                 curvature = 0.25, color = "firebrick", alpha = 0.5,
                 linewidth = 0.7, inherit.aes = FALSE)
    else
      geom_blank()
  }

# --- Pane C: naive collapse ---
pC <- ggplot(exon_df, aes(y = gene_rank)) +
  geom_hline(yintercept = unique(exon_df$gene_rank), color = "gray70", linewidth = 0.25) +
  geom_rect(data = exon_df[exon_df$gene_id == first_gene, ],
            aes(xmin = rel_start, xmax = rel_end,
                ymin = gene_rank - 0.35, ymax = gene_rank + 0.35),
            fill = hue_pal()(nrow(gene_order))[gene_order$gene_rank[gene_order$gene_id == first_gene]],
            inherit.aes = FALSE) +
  scale_y_continuous(breaks = gene_order$gene_rank,
                     labels = gene_order$gene_label,
                     expand = c(0.03, 0)) +
  scale_x_continuous(labels = comma, expand = c(0.02, 0)) +
  labs(x = "Position within gene (bp)", y = NULL,
       title = "C  Naive: reads collapsed to one locus") +
  theme_classic(base_size = 9) +
  theme(axis.text.y = element_text(size = 6, hjust = 1))

fig1 <- pA / pB / pC + plot_layout(heights = c(1, 1, 1))

ggsave(file.path(FIGURES_DIR, "fig1_problem.pdf"), fig1,
       width = 7, height = 12, units = "in")
ggsave(file.path(FIGURES_DIR, "fig1_problem.png"), fig1,
       width = 7, height = 12, units = "in", dpi = 300)
message("Figure 1 saved.")
