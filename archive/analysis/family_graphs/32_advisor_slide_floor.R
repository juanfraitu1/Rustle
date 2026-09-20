suppressPackageStartupMessages({
  library(ggplot2)
  library(ggraph)
  library(igraph)
  library(dplyr)
  library(patchwork)
  library(tidyr)
})

source("02_build_graphs.R")

FIGURES_DIR <- "figures"
DATA_DIR    <- "data"

# ── Three families to compare for the resolution-floor message ───────────────
fam_files <- list(
  GOLGA8I = "golga_GOLGA8I_vg_seq.rds",
  RBMY    = "rbmy_vg_seq.rds",
  GAGE    = "gage_vg_seq.rds"
)
vgs <- setNames(lapply(fam_files, function(f) readRDS(file.path(DATA_DIR, f))),
                 names(fam_files))

# Also need the strict RBMY for the "dial" panel
rbmy_strict <- readRDS(file.path(DATA_DIR, "rbmy_vg_seq_strict.rds"))

# ── Compute equivalence-class breakdown ───────────────────────────────────────
class_counts <- function(vg) {
  eq <- family_equivalence_classes(vg)
  n_copies  <- eq$n_copies
  n_classes <- eq$n_classes
  sizes <- sapply(eq$classes, length)
  list(n_copies = n_copies, n_classes = n_classes, sizes = sizes, eq = eq)
}
breakdowns <- setNames(lapply(vgs, class_counts), names(vgs))
strict_breakdown <- class_counts(rbmy_strict)

class_pal <- c("#1565C0","#B71C1C","#2E7D32","#6A1B9A","#EF6C00","#00838F",
                "#5D4037","#827717","#4A148C","#01579B","#33691E","#880E4F",
                "#1A237E")

# ── PANEL A: VG topology per family, colored by equivalence class ────────────
mk_class_topo <- function(name, vg, eq) {
  copy_to_class <- setNames(rep(NA_character_, length(vg$paths)), names(vg$paths))
  for (cn in names(eq$classes)) {
    for (mem in eq$classes[[cn]]) copy_to_class[mem] <- cn
  }
  V(vg$graph)$node_type <- vg$nodes$node_type
  V(vg$graph)$n_copies  <- vg$nodes$n_copies

  # Color a node by the class(es) that traverse it; for shared nodes, use gray
  node_class <- vapply(seq_len(igraph::vcount(vg$graph)), function(i) {
    nid <- as.integer(V(vg$graph)$name[i])
    classes_using <- unique(copy_to_class[
      vapply(names(vg$paths), function(g) nid %in% vg$paths[[g]], logical(1))
    ])
    if (length(classes_using) > 1L) return("multi")
    if (length(classes_using) == 1L) return(classes_using)
    "none"
  }, character(1))
  V(vg$graph)$class <- node_class

  class_levels <- names(eq$classes)
  class_pal_named <- setNames(class_pal[seq_along(class_levels)], class_levels)
  fill_map <- c(class_pal_named, multi = "#e6a817", none = "gray80")

  V(vg$graph)$fill <- fill_map[V(vg$graph)$class]

  ggraph(vg$graph, layout = "sugiyama") +
    geom_edge_link(arrow   = arrow(length = unit(1.3, "mm"), type = "closed"),
                   end_cap = circle(2, "mm"),
                   color   = "gray65", linewidth = 0.3) +
    geom_node_point(aes(fill = fill, size = n_copies),
                    shape = 21, color = "white", stroke = 0.4) +
    scale_fill_identity() +
    scale_size_continuous(range = c(1.5, 4.5), guide = "none") +
    labs(title = sprintf("%s  (%d copies -> %d classes)",
                          name, eq$n_copies, eq$n_classes)) +
    theme_graph(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11))
}

p_g8i  <- mk_class_topo("GOLGA8I", vgs$GOLGA8I, breakdowns$GOLGA8I$eq)
p_rbmy <- mk_class_topo("RBMY",    vgs$RBMY,    breakdowns$RBMY$eq)
p_gage <- mk_class_topo("GAGE",    vgs$GAGE,    breakdowns$GAGE$eq)

# ── PANEL B: equivalence-class size distribution per family ──────────────────
class_size_df <- bind_rows(
  data.frame(family = "GOLGA8I",
             class = names(breakdowns$GOLGA8I$sizes),
             size  = unname(breakdowns$GOLGA8I$sizes)),
  data.frame(family = "RBMY (0.85)",
             class = names(breakdowns$RBMY$sizes),
             size  = unname(breakdowns$RBMY$sizes)),
  data.frame(family = "GAGE",
             class = names(breakdowns$GAGE$sizes),
             size  = unname(breakdowns$GAGE$sizes))
)
class_size_df$family <- factor(class_size_df$family,
                                levels = c("GOLGA8I","RBMY (0.85)","GAGE"))
class_size_df$class <- factor(class_size_df$class,
                               levels = paste0("class_", 1:20))

p_sizes <- ggplot(class_size_df, aes(x = class, y = size, fill = family)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.4) +
  geom_text(aes(label = size), vjust = -0.3, size = 3, color = "gray20") +
  scale_fill_manual(values = c("GOLGA8I" = "#1565C0",
                                "RBMY (0.85)" = "#B71C1C",
                                "GAGE" = "#2E7D32"),
                    guide = "none") +
  facet_wrap(~ family, scales = "free", ncol = 3) +
  labs(title    = "Equivalence-class size distribution per family",
       subtitle = "Class size = number of copies that share one VG path. Size 1 = fully resolvable; size > 1 = read-indistinguishable members.",
       x = NULL, y = "Class size (copies)") +
  theme_minimal(base_size = 10) +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(color = "gray40", size = 9.5),
        strip.text    = element_text(face = "bold", size = 10),
        panel.grid.major.x = element_blank(),
        axis.text.x   = element_text(size = 8, angle = 30, hjust = 1))

# ── PANEL C: the "identity threshold dial" — RBMY 0.85 vs 0.95 ───────────────
dial_df <- data.frame(
  threshold = c("0.85 (default)", "0.95 (strict)"),
  n_copies  = c(13, 13),
  n_classes = c(breakdowns$RBMY$n_classes, strict_breakdown$n_classes)
)
dial_df$threshold <- factor(dial_df$threshold,
                             levels = c("0.85 (default)", "0.95 (strict)"))

p_dial <- ggplot(dial_df, aes(x = threshold, y = n_classes,
                                fill = threshold)) +
  geom_col(width = 0.45, color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%d classes\nfor %d copies",
                                  n_classes, n_copies)),
            vjust = -0.3, size = 4, fontface = "bold", color = "gray20") +
  geom_hline(yintercept = 13, linetype = "dashed", color = "#2E7D32") +
  annotate("text", x = 1.5, y = 13.5,
           label = "full resolution (13 copies = 13 classes)",
           color = "#2E7D32", fontface = "italic", size = 3.2) +
  scale_fill_manual(values = c("0.85 (default)" = "#FB8C00",
                                "0.95 (strict)" = "#1B5E20"),
                    guide = "none") +
  scale_y_continuous(limits = c(0, 16), breaks = seq(0, 16, 4)) +
  labs(title    = "RBMY: identity-threshold dial",
       subtitle = "Stricter clustering -> more classes -> finer resolution at the cost of VG complexity.",
       x = NULL, y = "Resolvable classes") +
  theme_classic(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(color = "gray40", size = 9.5),
        axis.text.x   = element_text(face = "bold", size = 10))

# ── Header / theorem callout ─────────────────────────────────────────────────
p_header <- ggplot() +
  annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = 1.2,
           fill = "#FAFAFA", color = "gray50", linewidth = 0.5) +
  annotate("text", x = 0.2, y = 0.85,
           label = "DEFINITION:  Two copies are read-equivalent iff they traverse the same paths through the VG.",
           hjust = 0, vjust = 0.5, fontface = "bold", size = 4.0, color = "gray15") +
  annotate("text", x = 0.2, y = 0.55,
           label = "THEOREM 1:  Members of one equivalence class cannot be distinguished by any read whose evidence is in V(G) U E(G).",
           hjust = 0, vjust = 0.5, fontface = "italic", size = 3.5, color = "gray30") +
  annotate("text", x = 0.2, y = 0.25,
           label = "CONSEQUENCE:  The equivalence class is the minimum unit of resolution — a hard floor set by the VG. Threshold tunes the floor.",
           hjust = 0, vjust = 0.5, size = 3.5, color = "gray30") +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 1.2)) +
  theme_void()

# ── Combine ───────────────────────────────────────────────────────────────────
layout_design <- "
AAAA
BCDE
FFFG
"
fig <- p_header +
       p_g8i + p_rbmy + p_gage + plot_spacer() +
       p_sizes + p_dial +
  plot_layout(design = layout_design,
              heights = c(0.6, 1.2, 1.0)) +
  plot_annotation(
    title    = "Resolution floor:  the VG equivalence class",
    subtitle = "What can and cannot be resolved by ANY read-based method is fixed by the VG topology. Node fills below show class membership; gray nodes are shared by multiple classes.",
    theme = theme(plot.title    = element_text(size = 17, face = "bold"),
                  plot.subtitle = element_text(size = 12, color = "gray40",
                                                lineheight = 1.25))
  )

ggsave(file.path(FIGURES_DIR, "fig_advisor_floor.pdf"), fig,
       width = 17, height = 10, units = "in")
ggsave(file.path(FIGURES_DIR, "fig_advisor_floor.png"), fig,
       width = 17, height = 10, units = "in", dpi = 180)
message("Floor slide saved.")
