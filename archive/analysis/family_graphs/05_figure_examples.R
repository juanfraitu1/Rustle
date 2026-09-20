suppressPackageStartupMessages({
  library(ggraph)
  library(igraph)
  library(ggplot2)
  library(patchwork)
})

source("02_build_graphs.R")

DATA_DIR    <- "data"
FIGURES_DIR <- "figures"
GENOME_FA   <- "/mnt/c/Users/jfris/Desktop/GGO.fasta"

rbmy_vg  <- readRDS(file.path(DATA_DIR, "rbmy_vg_seq.rds"))
amy_vg   <- readRDS(file.path(DATA_DIR, "amy_vg.rds"))
golga_vg <- readRDS(file.path(DATA_DIR, "golga_vg.rds"))

# Helper: draw one graph pane with a verdict badge
graph_pane <- function(g, title, verdict_label, verdict_color, layout = "sugiyama") {
  ggraph(g, layout = layout) +
    geom_edge_link(arrow = arrow(length = unit(2, "mm"), type = "closed"),
                   end_cap = circle(3, "mm"),
                   color = "gray55", linewidth = 0.45) +
    geom_node_point(aes(color = node_type, size = exon_len)) +
    scale_color_manual(values = c(shared = "#e6a817", copy_specific = "steelblue"),
                       guide = "none") +
    scale_size_continuous(range = c(1.5, 6), guide = "none") +
    annotate("label",
             x = -Inf, y = Inf, label = verdict_label,
             hjust = -0.1, vjust = 1.5, size = 3.2,
             fill = verdict_color, color = "white", fontface = "bold",
             label.padding = unit(0.25, "lines")) +
    labs(title = title) +
    coord_cartesian(clip = "off") +
    theme_graph(base_size = 10) +
    theme(plot.title = element_text(size = 10, margin = margin(b = 4)))
}

# --- Pane A: RBMY — positive example (sequence-based, 13 copies, 11 shared nodes) ---
pA <- graph_pane(rbmy_vg$graph,
                 "A  RBMY — 13 copies (seq-based)",
                 "multi-copy family ✓", "#2166ac",
                 layout = "sugiyama")

# --- Pane B: single AMY gene (isoforms only, no shared node) ---
first_gene   <- unique(amy_vg$exon_df$gene_id)[1]
single_exons <- amy_vg$exon_df[amy_vg$exon_df$gene_id == first_gene, ]
single_vg    <- build_variation_graph(single_exons)

pB <- graph_pane(single_vg$graph,
                 "B  Single locus — isoforms",
                 "isoforms only ✗", "#d6604d")

# --- Pane C: two unrelated genes — no shared structure ---
# AMY gene 1 + smallest GOLGA gene (different chromosomes: guaranteed disconnected)
amy_gene1   <- unique(amy_vg$exon_df$gene_id)[1]
golga_sizes <- sort(table(golga_vg$exon_df$gene_id))
golga_gene1 <- names(golga_sizes)[1]

unrelated_exons <- rbind(
  amy_vg$exon_df[amy_vg$exon_df$gene_id   == amy_gene1,   ],
  golga_vg$exon_df[golga_vg$exon_df$gene_id == golga_gene1, ]
)
unrelated_vg <- build_variation_graph(unrelated_exons)

pC <- graph_pane(unrelated_vg$graph,
                 "C  Unrelated genes",
                 "no shared exons ✗", "#d6604d")

fig3 <- pA + pB + pC +
  plot_layout(ncol = 3, widths = c(2, 1, 1))

ggsave(file.path(FIGURES_DIR, "fig3_examples.pdf"), fig3,
       width = 14, height = 5.5, units = "in")
ggsave(file.path(FIGURES_DIR, "fig3_examples.png"), fig3,
       width = 14, height = 5.5, units = "in", dpi = 300)
message("Figure 3 saved.")
