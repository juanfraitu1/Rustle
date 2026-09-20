suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggraph)
  library(patchwork)
})

source("02_build_graphs.R")

DATA_DIR    <- "data"
FIGURES_DIR <- "figures"
GENOME_FA   <- "/mnt/c/Users/jfris/Desktop/GGO.fasta"

# ── Combine multiple known-unrelated families and ONE coherent one ───────────
fam_keys <- list(
  RBMY     = "rbmy_family.rds",
  AMY      = "amy_family.rds",
  GOLGA8I  = "golga_family.rds",      # we'll filter to GOLGA8I sub-family
  DEFB     = "defb_family.rds",
  APOL     = "apol_family.rds",
  GAGE     = "gage_family.rds"
)

# For GOLGA8I, restrict to its sub-family members
classify_golga <- function(desc) {
  if (is.na(desc)) return(NA_character_)
  if (grepl("RAB6", desc)) return("GORAB")
  m <- regmatches(desc, regexpr("subfamily A member [0-9]+[A-Z]*", desc, perl=TRUE))
  if (length(m) > 0L) {
    num <- sub("subfamily A member ([0-9]+).*", "\\1", m, perl=TRUE)
    letter <- sub("subfamily A member [0-9]+([A-Z]*).*", "\\1", m, perl=TRUE)
    likes <- ifelse(grepl("-like", desc), "L", "")
    return(paste0("GOLGA", num, letter, likes))
  }
  "OTHER"
}

# Build the combined exon table
combined <- list()
for (key in names(fam_keys)) {
  fam <- readRDS(file.path(DATA_DIR, fam_keys[[key]]))
  ex  <- fam$exon_df
  if (key == "GOLGA8I") {
    # filter to GOLGA8I sub-family only
    gdesc <- as.data.frame(GenomicRanges::mcols(fam$genes))
    gene_subfam <- setNames(sapply(gdesc$description, classify_golga),
                             fam$genes$ID)
    ex <- ex[gene_subfam[ex$gene_id] == "GOLGA8I", ]
  }
  ex$family <- key
  combined[[key]] <- ex
  cat(sprintf("%-8s : %d genes, %d exons\n",
              key, length(unique(ex$gene_id)), nrow(ex)))
}
combined_df <- do.call(rbind, combined)
combined_df <- combined_df[, c("gene_id","tx_id","chrom","start","end","strand","family")]

# ── Save the source-of-truth gene-to-family map ───────────────────────────────
gene_family_map <- combined_df %>%
  group_by(gene_id) %>%
  summarise(family = family[1], .groups = "drop")
cat("\nTotal genes in combined input: ", nrow(gene_family_map), "\n")
cat("Family sizes:\n")
print(table(gene_family_map$family))

# Drop the helper column for VG build
build_input <- combined_df[, c("gene_id","tx_id","chrom","start","end","strand")]

# ── Build combined VG (default 0.85 threshold) ───────────────────────────────
out_path <- file.path(DATA_DIR, "cross_family_vg.rds")
if (!file.exists(out_path)) {
  cat("\nBuilding combined VG (this takes a few minutes)...\n")
  combined_vg <- build_variation_graph_seq(build_input, GENOME_FA,
                                            min_identity = 0.85,
                                            min_coverage = 0.80)
  saveRDS(combined_vg, out_path)
} else {
  combined_vg <- readRDS(out_path)
}

cat("\nCombined VG: ", nrow(combined_vg$nodes), " nodes (",
    sum(combined_vg$nodes$node_type == "shared"), " shared, ",
    sum(combined_vg$nodes$node_type == "copy_specific"), " copy-spec)\n", sep="")

# ── Map every gene to its connected component in the combined VG ─────────────
g  <- combined_vg$graph
nd <- combined_vg$nodes
comps <- igraph::components(g, mode = "weak")
node_to_comp <- setNames(as.integer(comps$membership), V(g)$name)

copy_to_comp <- vapply(combined_vg$paths, function(p) {
  if (length(p) == 0L) return(NA_integer_)
  unname(node_to_comp[as.character(p[1])])
}, integer(1))

gene_comp_df <- data.frame(gene_id = names(copy_to_comp),
                            comp    = unname(copy_to_comp),
                            stringsAsFactors = FALSE) %>%
  left_join(gene_family_map, by = "gene_id")

# ── Per-component: which families appear, and how many genes ─────────────────
comp_summary <- gene_comp_df %>%
  group_by(comp) %>%
  summarise(n_genes     = n(),
             families    = paste(sort(unique(family)), collapse = ","),
             n_families  = length(unique(family)),
             pure        = n_families == 1L,
             .groups = "drop") %>%
  arrange(desc(n_genes))

cat("\n=== COMPONENT BREAKDOWN ===\n")
print(comp_summary, n = Inf)

# Components that contain >1 family = FALSE POSITIVES at the family level
fp_components <- comp_summary %>% filter(!pure)
cat("\n=== FALSE-POSITIVE COMPONENTS (mixing families) ===\n")
if (nrow(fp_components) == 0L) {
  cat("None. Every connected component contains genes from exactly ONE family.\n")
} else {
  print(fp_components, n = Inf)
}

# ── Equivalence classes: check whether any class has cross-family members ────
eq <- family_equivalence_classes(combined_vg)
class_purity <- do.call(rbind, lapply(names(eq$classes), function(cn) {
  members <- eq$classes[[cn]]
  fams    <- gene_family_map$family[match(members, gene_family_map$gene_id)]
  data.frame(class       = cn,
             n_members   = length(members),
             families    = paste(sort(unique(fams)), collapse = ","),
             n_families  = length(unique(fams)),
             pure        = length(unique(fams)) == 1L,
             stringsAsFactors = FALSE)
}))
cat("\n=== EQUIVALENCE-CLASS PURITY ===\n")
cat("Total classes: ", nrow(class_purity), "\n")
cat("Pure classes (single family): ", sum(class_purity$pure), "\n")
cat("Mixed classes (multi-family): ", sum(!class_purity$pure), "\n\n")

cat("Mixed classes (if any):\n")
mixed_classes <- class_purity %>% filter(!pure) %>% arrange(desc(n_members))
if (nrow(mixed_classes) == 0L) {
  cat("None. Every equivalence class contains genes from exactly ONE family.\n")
} else {
  print(mixed_classes, n = Inf)
}

# Write TSVs
write.table(comp_summary, file.path(FIGURES_DIR, "cross_family_components.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(class_purity, file.path(FIGURES_DIR, "cross_family_classes.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# ── Visualisation: combined VG colored by family ─────────────────────────────
fam_pal <- c(RBMY="#B71C1C", AMY="#1565C0", GOLGA8I="#2E7D32",
              DEFB="#EF6C00", APOL="#6A1B9A", GAGE="#00838F")

# For each node, determine which family(ies) use it
node_families <- vapply(seq_len(igraph::vcount(g)), function(i) {
  nid <- as.integer(V(g)$name[i])
  using <- vapply(names(combined_vg$paths),
                  function(gid) nid %in% combined_vg$paths[[gid]],
                  logical(1))
  fams <- unique(gene_family_map$family[match(names(combined_vg$paths)[using],
                                                 gene_family_map$gene_id)])
  if (length(fams) > 1L) return("MIXED")
  if (length(fams) == 1L) return(fams)
  return("none")
}, character(1))
V(g)$family <- node_families
V(g)$fill   <- ifelse(node_families == "MIXED", "#000000",
                      ifelse(node_families == "none", "gray80",
                             fam_pal[node_families]))
V(g)$n_copies <- nd$n_copies

p_topo <- ggraph(g, layout = "sugiyama") +
  geom_edge_link(arrow = arrow(length = unit(1.0, "mm"), type = "closed"),
                 end_cap = circle(1.5, "mm"),
                 color = "gray70", linewidth = 0.2) +
  geom_node_point(aes(fill = fill, size = n_copies),
                  shape = 21, color = "white", stroke = 0.3) +
  scale_fill_identity() +
  scale_size_continuous(range = c(1, 3), guide = "none") +
  labs(title    = "Combined VG of 6 known-unrelated gene families",
       subtitle = sprintf("%d total genes, %d nodes, %d components. Node fill = family (black = node used by >1 family).",
                            nrow(gene_family_map),
                            nrow(nd), comps$no)) +
  theme_graph(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(color = "gray40", size = 10))

# Bar chart: per-family genes and per-family components
summary_per_fam <- gene_comp_df %>%
  group_by(family) %>%
  summarise(n_genes = n(), n_components = length(unique(comp)),
            .groups = "drop") %>%
  pivot_longer(cols = c(n_genes, n_components), names_to = "metric",
                values_to = "count")

p_bars <- ggplot(summary_per_fam,
                  aes(x = family, y = count,
                      fill = fam_pal[family],
                      alpha = ifelse(metric == "n_genes", 1.0, 0.55))) +
  geom_col(width = 0.7, color = "white", position = position_dodge(width = 0.8)) +
  scale_fill_identity() +
  scale_alpha_identity() +
  facet_wrap(~ metric, scales = "free_y",
             labeller = labeller(metric = c(n_genes = "Genes in family",
                                              n_components = "VG components"))) +
  geom_text(aes(label = count),
            position = position_dodge(width = 0.8),
            vjust = -0.3, size = 3.2, color = "gray20") +
  labs(title    = "Per-family decomposition in the combined VG",
       subtitle = "Number of VG connected components a family decomposes into. Should be >= 1.",
       x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(color = "gray40", size = 10),
        axis.text.x   = element_text(face = "bold", size = 10),
        strip.text    = element_text(face = "bold"))

# Verdict box
verdict_lines <- c()
if (nrow(fp_components) == 0L) {
  verdict_lines <- c(verdict_lines,
    "COMPONENT-LEVEL FALSE POSITIVES: 0  -> Every VG component contains only one family.")
} else {
  verdict_lines <- c(verdict_lines,
    sprintf("COMPONENT-LEVEL FALSE POSITIVES: %d  -> %d connected components mix families!",
            nrow(fp_components), nrow(fp_components)))
}
if (nrow(mixed_classes) == 0L) {
  verdict_lines <- c(verdict_lines,
    "EQUIVALENCE-CLASS FALSE POSITIVES: 0  -> Every equivalence class contains only one family.")
} else {
  verdict_lines <- c(verdict_lines,
    sprintf("EQUIVALENCE-CLASS FALSE POSITIVES: %d", nrow(mixed_classes)))
}

p_verdict <- ggplot() +
  annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = 1,
           fill = ifelse(nrow(fp_components) == 0 & nrow(mixed_classes) == 0,
                          "#E8F5E9", "#FFEBEE"),
           color = ifelse(nrow(fp_components) == 0 & nrow(mixed_classes) == 0,
                           "#2E7D32", "#B71C1C"),
           linewidth = 0.7) +
  annotate("text", x = 0.15, y = 0.75,
           label = verdict_lines[1], hjust = 0, vjust = 0.5,
           size = 4.0, fontface = "bold", color = "gray15") +
  annotate("text", x = 0.15, y = 0.30,
           label = verdict_lines[2], hjust = 0, vjust = 0.5,
           size = 4.0, fontface = "bold", color = "gray15") +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 1)) +
  theme_void()

fig <- (p_topo / p_bars / p_verdict) +
  plot_layout(heights = c(2.0, 1.0, 0.4)) +
  plot_annotation(
    title    = "Cross-family specificity: do families stay separate in the combined VG?",
    subtitle = "Combined input = 6 known-unrelated gene families. Test: do they cluster correctly (per-family components, per-family classes)?",
    theme = theme(plot.title    = element_text(size = 14, face = "bold"),
                  plot.subtitle = element_text(size = 11, color = "gray40"))
  )

ggsave(file.path(FIGURES_DIR, "fig_cross_family_specificity.pdf"), fig,
       width = 14, height = 12, units = "in")
ggsave(file.path(FIGURES_DIR, "fig_cross_family_specificity.png"), fig,
       width = 14, height = 12, units = "in", dpi = 180)
message("\nCross-family specificity figure saved.")
