suppressPackageStartupMessages({
  library(ggraph)
  library(igraph)
  library(ggplot2)
  library(patchwork)
})

source("02_build_graphs.R")

DATA_DIR    <- "data"
FIGURES_DIR <- "figures"

golga_vg <- readRDS(file.path(DATA_DIR, "golga_vg.rds"))
amy_vg   <- readRDS(file.path(DATA_DIR, "amy_vg.rds"))

graph_pane <- function(g, title, verdict_label, verdict_color) {
  ggraph(g, layout = "sugiyama") +
    geom_edge_link(arrow = arrow(length = unit(2, "mm"), type = "closed"),
                   end_cap = circle(3, "mm"),
                   color = "gray60", linewidth = 0.4) +
    geom_node_point(aes(color = node_type, size = exon_len)) +
    scale_color_manual(values = c(shared = "gray50", copy_specific = "steelblue"),
                       guide = "none") +
    scale_size_continuous(range = c(2, 7), guide = "none") +
    annotate("label", x = -Inf, y = Inf, label = verdict_label,
             hjust = -0.1, vjust = 1.5, size = 3,
             fill = verdict_color, color = "white", fontface = "bold",
             label.padding = unit(0.3, "lines")) +
    labs(title = title) +
    theme_graph(base_size = 9)
}

# --- Pane A: GOLGA (positive example — large multi-copy tandem family) ---
# Use the largest cluster on NC_073240.2 (33 genes) for a rich graph
golga_big_chrom <- golga_vg$exon_df[golga_vg$exon_df$chrom == "NC_073240.2", ]
golga_cluster_genes <- unique(golga_big_chrom$gene_id)
# Limit to first 8 genes to keep the graph readable
golga_subset_genes <- golga_cluster_genes[seq_len(min(8L, length(golga_cluster_genes)))]
golga_subset_exons <- golga_vg$exon_df[golga_vg$exon_df$gene_id %in% golga_subset_genes, ]
golga_subset_vg    <- build_variation_graph(golga_subset_exons)

pA <- graph_pane(golga_subset_vg$graph,
                 "A  GOLGA — multi-copy tandem family (8 of 46 paralogs)",
                 "multi-copy family", "#2166ac")

# --- Pane B: single-gene counterexample ---
first_gene   <- unique(amy_vg$exon_df$gene_id)[1]
single_exons <- amy_vg$exon_df[amy_vg$exon_df$gene_id == first_gene, ]
single_vg    <- build_variation_graph(single_exons)

pB <- graph_pane(single_vg$graph,
                 "B  Single locus — isoforms only",
                 "isoforms only", "#d6604d")

# --- Pane C: two unrelated genes on different chromosomes (guaranteed disconnected) ---
# gene-LOC101133335 is on NC_073224.2; gene-GOLGA1 is on NC_073237.2 — different chroms
amy_gene1   <- "gene-LOC101133335"
golga_gene1 <- "gene-GOLGA1"
unrelated_exons <- rbind(
  amy_vg$exon_df[amy_vg$exon_df$gene_id == amy_gene1, ],
  golga_vg$exon_df[golga_vg$exon_df$gene_id == golga_gene1, ]
)
unrelated_vg <- build_variation_graph(unrelated_exons)

pC <- graph_pane(unrelated_vg$graph,
                 "C  Unrelated genes — no shared exons",
                 "no shared structure", "#d6604d")

fig3 <- pA + pB + pC + plot_layout(ncol = 3)

ggsave(file.path(FIGURES_DIR, "fig3_examples.pdf"), fig3,
       width = 12, height = 5, units = "in")
ggsave(file.path(FIGURES_DIR, "fig3_examples.png"), fig3,
       width = 12, height = 5, units = "in", dpi = 300)
message("Figure 3 saved.")
