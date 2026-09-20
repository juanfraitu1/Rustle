suppressPackageStartupMessages({
  library(ggplot2)
  library(ggraph)
  library(igraph)
  library(dplyr)
  library(patchwork)
})

source("02_build_graphs.R")

DATA_DIR    <- "data"
FIGURES_DIR <- "figures"

# Six representative VGs covering the three types
families <- list(
  list(key = "RBMY",     label = "RBMY  (13 copies, 1 component)",        type = "Coherent family"),
  list(key = "GOLGA8I",  label = "GOLGA8I  (7 copies, 1 component)",      type = "Coherent family"),
  list(key = "MAGEA",    label = "MAGEA  (5 copies, 4 components)",       type = "Mixed cluster"),
  list(key = "APOL",     label = "APOL  (6 copies, 3 components)",        type = "Mixed cluster"),
  list(key = "DEFB",     label = "DEFB  (10 copies, 10 components)",      type = "Fully fragmented cluster"),
  list(key = "PRAMEF",   label = "PRAMEF  (4 copies, 4 components)",      type = "Fully fragmented cluster")
)

paths_for_key <- function(k) {
  switch(k,
    RBMY    = file.path(DATA_DIR, "rbmy_vg_seq.rds"),
    GOLGA8I = file.path(DATA_DIR, "golga_GOLGA8I_vg_seq.rds"),
    MAGEA   = file.path(DATA_DIR, "magea_vg_seq.rds"),
    APOL    = file.path(DATA_DIR, "apol_vg_seq.rds"),
    DEFB    = file.path(DATA_DIR, "defb_vg_seq.rds"),
    PRAMEF  = file.path(DATA_DIR, "pramef_vg_seq.rds")
  )
}

type_colors <- c("Coherent family"          = "#2E7D32",
                 "Mixed cluster"            = "#EF6C00",
                 "Fully fragmented cluster" = "#B71C1C")

# Build a per-component palette for visual distinction
component_palette <- c("#1565C0","#B71C1C","#2E7D32","#6A1B9A","#EF6C00",
                       "#00838F","#5D4037","#827717","#4A148C","#01579B")

mk_panel <- function(spec) {
  vg <- readRDS(paths_for_key(spec$key))
  g  <- vg$graph
  nd <- vg$nodes

  comps <- igraph::components(g, mode = "weak")
  V(g)$comp <- as.integer(comps$membership)
  V(g)$node_type <- nd$node_type
  V(g)$n_copies  <- nd$n_copies

  cl <- analyse_locus_cluster(vg)

  # color each node by its component (cycled palette)
  comp_color <- component_palette[((V(g)$comp - 1L) %% length(component_palette)) + 1L]
  V(g)$fillc <- comp_color

  ggraph(g, layout = "sugiyama") +
    geom_edge_link(arrow   = arrow(length = unit(1.3, "mm"), type = "closed"),
                   end_cap = circle(2, "mm"),
                   color   = "gray70", linewidth = 0.25) +
    geom_node_point(aes(fill = fillc,
                          shape = node_type,
                          size  = n_copies),
                    color = "white", stroke = 0.5) +
    scale_fill_identity() +
    scale_shape_manual(values = c(shared = 21, copy_specific = 24), guide = "none") +
    scale_size_continuous(range = c(1.4, 4.0), guide = "none") +
    labs(title    = spec$label,
         subtitle = sprintf("%d total classes  %s",
                              cl$n_total_classes,
                              ifelse(cl$fully_resolved, "(fully resolved)",
                                                         "(some collapse)"))) +
    theme_graph(base_size = 9) +
    theme(plot.title       = element_text(face = "bold", size = 9.5),
          plot.subtitle    = element_text(color = "gray45", size = 8))
}

panels <- lapply(families, mk_panel)

# Add a small banner above each pair indicating the type
banner <- function(text, col) {
  ggplot() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
             fill = col, alpha = 0.85) +
    annotate("text", x = 0.5, y = 0.5, label = text,
             color = "white", fontface = "bold", size = 4.5) +
    coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
    theme_void()
}

# Build a 3-row layout where each row has a banner + 2 family panels
row_layout <- function(banner_text, banner_color, p1, p2) {
  b <- banner(banner_text, banner_color)
  (b | (p1 + p2)) + plot_layout(widths = c(0.18, 0.41, 0.41))
}

row1 <- row_layout("COHERENT\nFAMILY\n(K = 1)",
                    type_colors["Coherent family"],
                    panels[[1]], panels[[2]])
row2 <- row_layout("MIXED\nCLUSTER\n(K > 1, some n_k >= 2)",
                    type_colors["Mixed cluster"],
                    panels[[3]], panels[[4]])
row3 <- row_layout("FULLY\nFRAGMENTED\nCLUSTER\n(K = N)",
                    type_colors["Fully fragmented cluster"],
                    panels[[5]], panels[[6]])

fig <- (row1 / row2 / row3) +
  plot_annotation(
    title    = "Taxonomy of multi-copy gene loci",
    subtitle = "A locus L of N copies is classified by its VG component decomposition. Node fill = component (each component a different color). Triangles = copy-specific exons, circles = shared exons. Size = n_copies sharing that node.",
    theme = theme(
      plot.title    = element_text(size = 15, face = "bold"),
      plot.subtitle = element_text(size = 10.5, color = "gray40",
                                   lineheight = 1.3)
    )
  )

ggsave(file.path(FIGURES_DIR, "fig_family_taxonomy.pdf"), fig,
       width = 14, height = 11, units = "in")
ggsave(file.path(FIGURES_DIR, "fig_family_taxonomy.png"), fig,
       width = 14, height = 11, units = "in", dpi = 180)
message("Taxonomy figure saved.")
