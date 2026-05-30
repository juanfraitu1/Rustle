suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
  library(ggplot2)
  library(ggraph)
  library(patchwork)
  library(GenomicRanges)
  library(S4Vectors)
})

source("02_build_graphs.R")

DATA_DIR    <- "data"
FIGURES_DIR <- "figures"
GENOME_FA   <- "/mnt/c/Users/jfris/Desktop/GGO.fasta"
MIN_COPIES  <- 3L

# ── Sub-family classifiers ────────────────────────────────────────────────────
classify_golga <- function(desc) {
  if (is.na(desc)) return(NA_character_)
  if (grepl("RAB6", desc)) return("GORAB")
  # Old format: "golgin A4", "golgin B1"
  m1 <- regmatches(desc, regexpr("golgin\\s+[AB]\\s*[0-9]+", desc, ignore.case=TRUE, perl=TRUE))
  if (length(m1) > 0L) {
    letter <- toupper(sub("golgin\\s+([AB])\\s*[0-9]+", "\\1", m1, perl=TRUE))
    num    <- sub("golgin\\s+[AB]\\s*([0-9]+)", "\\1", m1, perl=TRUE)
    return(paste0("GOLG", letter, num))
  }
  # New format: "subfamily A member N" + optional letter + optional "-like"
  m2 <- regmatches(desc, regexpr("subfamily A member [0-9]+[A-Z]*(-like)?", desc, perl=TRUE))
  if (length(m2) > 0L) {
    num    <- sub("subfamily A member ([0-9]+).*", "\\1", m2, perl=TRUE)
    letter <- sub("subfamily A member [0-9]+([A-Z]*).*", "\\1", m2, perl=TRUE)
    likes  <- ifelse(grepl("-like", m2), "L", "")
    return(paste0("GOLGA", num, letter, likes))
  }
  "UNCLASSIFIED"
}

classify_npip <- function(desc) {
  if (is.na(desc)) return(NA_character_)
  m <- regmatches(desc, regexpr("NPIP[A-Z][0-9]*[A-Z]*", desc))
  if (length(m) > 0L) return(toupper(m[1]))
  m2 <- regmatches(desc,
    regexpr("nuclear pore complex interacting protein family member [A-Z][0-9]*",
            desc, perl=TRUE))
  if (length(m2) > 0L) {
    return(paste0("NPIP", sub(".*member ", "", m2)))
  }
  "UNCLASSIFIED"
}

# ── Helper to split a family by sub-family + build per-sub-family VGs ─────────
build_subfamily_vgs <- function(family_rds, classifier, parent_name) {
  fam <- readRDS(family_rds)
  gdesc <- as.data.frame(mcols(fam$genes))
  desc_df <- data.frame(gene_id     = fam$genes$ID,
                        description = gdesc$description,
                        stringsAsFactors = FALSE)
  desc_df$subfam <- sapply(desc_df$description, classifier)
  gene_to_subfam <- setNames(desc_df$subfam, desc_df$gene_id)

  exon_df <- fam$exon_df
  exon_df$subfam <- gene_to_subfam[exon_df$gene_id]

  cat("\n", parent_name, " sub-family breakdown:\n", sep="")
  sf_counts <- exon_df %>%
    group_by(subfam) %>%
    summarise(n_copies = length(unique(gene_id)), .groups="drop") %>%
    arrange(desc(n_copies))
  print(sf_counts, row.names=FALSE)

  out <- list()
  for (sf in sf_counts$subfam[sf_counts$n_copies >= MIN_COPIES]) {
    sf_ex <- exon_df[exon_df$subfam == sf, ]
    sf_ex$subfam <- NULL
    n_genes <- length(unique(sf_ex$gene_id))
    cache_path <- file.path(DATA_DIR, sprintf("%s_%s_vg_seq.rds", parent_name, sf))
    if (file.exists(cache_path)) {
      message(sf, " (", n_genes, " copies): cached")
      out[[sf]] <- readRDS(cache_path)
    } else {
      message(sf, " (", n_genes, " copies): building VG ...")
      vg <- build_variation_graph_seq(sf_ex, GENOME_FA)
      saveRDS(vg, cache_path)
      out[[sf]] <- vg
    }
  }
  out
}

# ── Metric helpers ────────────────────────────────────────────────────────────
count_bubbles <- function(vg) {
  g <- vg$graph; nd <- vg$nodes
  if (vcount(g) == 0L) return(0L)
  out_list <- ego(g, order = 1, mode = "out")
  n_bub <- 0L
  for (i in seq_along(out_list)) {
    if (nd$node_type[i] != "shared") next
    succs <- setdiff(as.integer(out_list[[i]]), i)
    spec  <- succs[nd$node_type[succs] == "copy_specific"]
    if (length(spec) < 2L) next
    downs <- lapply(spec, function(s) {
      d <- setdiff(as.integer(ego(g, order = 1, mode = "out")[[s]]), s)
      d[nd$node_type[d] == "shared"]
    })
    if (length(Reduce(intersect, downs)) >= 1L) n_bub <- n_bub + 1L
  }
  n_bub
}

metric_row <- function(parent, sf, vg) {
  g <- vg$graph; nd <- vg$nodes
  data.frame(
    parent       = parent,
    subfamily    = sf,
    n_copies     = length(vg$paths),
    n_nodes      = nrow(nd),
    n_shared     = sum(nd$node_type == "shared"),
    n_spec       = sum(nd$node_type == "copy_specific"),
    pct_spec     = round(100 * sum(nd$node_type == "copy_specific") / nrow(nd), 1),
    full_shared  = sum(nd$n_copies == length(vg$paths)),
    n_components = igraph::components(g, mode = "weak")$no,
    n_bubbles    = count_bubbles(vg),
    is_multicopy = is_multicopy_family(vg),
    stringsAsFactors = FALSE
  )
}

# ── Build per-sub-family VGs ──────────────────────────────────────────────────
golga_subs <- build_subfamily_vgs(file.path(DATA_DIR, "golga_family.rds"),
                                    classify_golga, "golga")
npip_subs  <- build_subfamily_vgs(file.path(DATA_DIR, "npip_family.rds"),
                                    classify_npip, "npip")

# ── Metrics table ─────────────────────────────────────────────────────────────
stats <- bind_rows(
  do.call(rbind, lapply(names(golga_subs),
                          function(sf) metric_row("GOLGA", sf, golga_subs[[sf]]))),
  do.call(rbind, lapply(names(npip_subs),
                          function(sf) metric_row("NPIP",  sf, npip_subs[[sf]])))
)

# Add RBMY + AMY for comparison
ref_vgs <- list(RBMY = readRDS(file.path(DATA_DIR, "rbmy_vg_seq.rds")),
                AMY  = readRDS(file.path(DATA_DIR, "amy_vg_seq.rds")))
ref_stats <- do.call(rbind, lapply(names(ref_vgs),
                                     function(nm) metric_row("REF", nm, ref_vgs[[nm]])))
stats <- bind_rows(ref_stats, stats)

cat("\n=========== PER-SUB-FAMILY VG METRICS ===========\n\n")
print(stats, row.names = FALSE)
write.table(stats, file.path(FIGURES_DIR, "subfamily_review.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ── Topology thumbnails ───────────────────────────────────────────────────────
mk_topo <- function(name, vg) {
  ggraph(vg$graph, layout = "sugiyama") +
    geom_edge_link(arrow   = arrow(length = unit(1.5, "mm"), type = "closed"),
                   end_cap = circle(2, "mm"),
                   color   = "gray60", linewidth = 0.3) +
    geom_node_point(aes(color = node_type, size = n_copies)) +
    scale_color_manual(values = c(shared = "#e6a817", copy_specific = "#1565C0"),
                       guide = "none") +
    scale_size_continuous(range = c(1, 4), guide = "none") +
    labs(title = sprintf("%s  (%d copies | %d nodes | %d bubbles | %d cc)",
                          name,
                          length(vg$paths), nrow(vg$nodes),
                          count_bubbles(vg),
                          igraph::components(vg$graph, mode = "weak")$no)) +
    theme_graph(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 9))
}

all_vgs   <- c(ref_vgs,
                setNames(golga_subs, paste0("GOLGA-", names(golga_subs))),
                setNames(npip_subs,  paste0("NPIP-",  names(npip_subs))))
topos     <- Map(mk_topo, names(all_vgs), all_vgs)

p_topos <- wrap_plots(topos, ncol = 3)
n_rows  <- ceiling(length(topos) / 3)
ggsave(file.path(FIGURES_DIR, "subfamily_topologies.png"), p_topos,
       width = 14, height = n_rows * 3.5, units = "in", dpi = 180)
message("Saved figures/subfamily_topologies.png")
