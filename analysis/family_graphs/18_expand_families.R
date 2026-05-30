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

gff <- readRDS(file.path(DATA_DIR, "gff_cache.rds"))

# ── Helper to parse a family, restricted to the dominant chromosome ───────────
first_parent <- function(x) {
  sapply(as.list(x), function(p) if (length(p) > 0L) p[[1L]] else NA_character_)
}

parse_one_family <- function(gff, name_pattern = NULL, desc_pattern = NULL,
                              restrict_chr = TRUE, label) {
  genes <- gff[gff$type == "gene"]
  mdc <- as.data.frame(mcols(genes))
  nm <- mdc$Name; nm[is.na(nm)] <- ""
  ds <- mdc$description; ds[is.na(ds)] <- ""
  hit <- if (!is.null(name_pattern)) grepl(name_pattern, nm) else logical(length(genes))
  if (!is.null(desc_pattern)) hit <- hit | grepl(desc_pattern, ds, ignore.case = TRUE)
  family_genes <- genes[hit]
  if (length(family_genes) == 0L) {
    message(label, ": 0 hits — skipping")
    return(NULL)
  }
  if (restrict_chr) {
    chrs <- as.character(seqnames(family_genes))
    top_chr <- names(sort(table(chrs), decreasing = TRUE))[1]
    family_genes <- family_genes[chrs == top_chr]
    message(label, ": ", length(family_genes), " genes on ", top_chr)
  } else {
    message(label, ": ", length(family_genes), " genes (multi-chr)")
  }
  gene_ids <- family_genes$ID

  txs <- gff[gff$type %in% c("mRNA","transcript")]
  tx_parent <- first_parent(txs$Parent)
  family_txs <- txs[!is.na(tx_parent) & tx_parent %in% gene_ids]
  tx_gene <- setNames(tx_parent[!is.na(tx_parent) & tx_parent %in% gene_ids],
                      family_txs$ID)

  exons <- gff[gff$type == "exon"]
  ex_parent <- first_parent(exons$Parent)
  family_exons <- exons[!is.na(ex_parent) & ex_parent %in% names(tx_gene)]
  exon_df <- data.frame(
    gene_id = tx_gene[first_parent(family_exons$Parent)],
    tx_id   = first_parent(family_exons$Parent),
    chrom   = as.character(seqnames(family_exons)),
    start   = start(family_exons),
    end     = end(family_exons),
    strand  = as.character(strand(family_exons)),
    stringsAsFactors = FALSE
  )
  exon_df <- exon_df[!is.na(exon_df$gene_id), ]
  list(genes = family_genes, exon_df = exon_df)
}

# ── Candidate families ────────────────────────────────────────────────────────
candidates <- list(
  DEFB    = list(name = "^DEFB[0-9]", label = "DEFB (beta-defensin)"),
  GAGE    = list(desc = "G antigen [0-9]", label = "GAGE family"),
  MAGEA   = list(name = "^MAGEA[0-9]", label = "MAGEA family"),
  APOL    = list(desc = "apolipoprotein L", label = "APOL family"),
  PRAMEF  = list(name = "^PRAMEF", label = "PRAMEF family")
)

# ── Parse + build VGs + validate ─────────────────────────────────────────────
results <- list()
for (key in names(candidates)) {
  spec <- candidates[[key]]
  fam <- parse_one_family(gff,
                          name_pattern = spec$name,
                          desc_pattern = spec$desc,
                          label = spec$label)
  if (is.null(fam) || length(unique(fam$exon_df$gene_id)) < 3L) {
    message(key, ": <3 genes — skipping")
    next
  }
  saveRDS(fam, file.path(DATA_DIR, paste0(tolower(key), "_family.rds")))

  # Loose (0.85) build
  loose_path <- file.path(DATA_DIR, paste0(tolower(key), "_vg_seq.rds"))
  if (!file.exists(loose_path)) {
    message("  Building VG (loose 0.85) ...")
    vg_loose <- build_variation_graph_seq(fam$exon_df, GENOME_FA,
                                            min_identity = 0.85,
                                            min_coverage = 0.80)
    saveRDS(vg_loose, loose_path)
  } else {
    vg_loose <- readRDS(loose_path)
  }
  # Strict (0.95) build
  strict_path <- file.path(DATA_DIR, paste0(tolower(key), "_vg_seq_strict.rds"))
  if (!file.exists(strict_path)) {
    message("  Building VG (strict 0.95) ...")
    vg_strict <- build_variation_graph_seq(fam$exon_df, GENOME_FA,
                                             min_identity = 0.95,
                                             min_coverage = 0.95)
    saveRDS(vg_strict, strict_path)
  } else {
    vg_strict <- readRDS(strict_path)
  }
  results[[key]] <- list(loose = vg_loose, strict = vg_strict, fam = fam)
}

# ── Summarise ─────────────────────────────────────────────────────────────────
summarise_vg <- function(vg) {
  if (is.null(vg) || nrow(vg$nodes) == 0L) return(NULL)
  eq <- family_equivalence_classes(vg)
  data.frame(
    n_copies        = eq$n_copies,
    n_nodes         = nrow(vg$nodes),
    n_shared        = sum(vg$nodes$node_type == "shared"),
    n_spec          = sum(vg$nodes$node_type == "copy_specific"),
    pct_shared      = round(100 * sum(vg$nodes$node_type == "shared") / nrow(vg$nodes), 1),
    n_comp          = igraph::components(vg$graph, mode = "weak")$no,
    n_classes       = eq$n_classes,
    fully_resolved  = eq$fully_resolved,
    is_multicopy    = is_multicopy_family(vg),
    stringsAsFactors = FALSE
  )
}

# Existing reference families
ref_list <- list(
  RBMY     = readRDS(file.path(DATA_DIR, "rbmy_vg_seq.rds")),
  AMY      = readRDS(file.path(DATA_DIR, "amy_vg_seq.rds")),
  GOLGA8I  = readRDS(file.path(DATA_DIR, "golga_GOLGA8I_vg_seq.rds")),
  GOLGA8NL = readRDS(file.path(DATA_DIR, "golga_GOLGA8NL_vg_seq.rds")),
  GOLGA6C  = readRDS(file.path(DATA_DIR, "golga_GOLGA6C_vg_seq.rds"))
)

all_summary <- bind_rows(
  do.call(rbind, lapply(names(ref_list), function(nm) {
    s <- summarise_vg(ref_list[[nm]])
    cbind(family = nm, threshold = "0.85", s)
  })),
  do.call(rbind, lapply(names(results), function(nm) {
    s1 <- summarise_vg(results[[nm]]$loose)
    s2 <- summarise_vg(results[[nm]]$strict)
    rbind(cbind(family = nm, threshold = "0.85", s1),
          cbind(family = nm, threshold = "0.95", s2))
  }))
)

cat("\n========== EXPANDED FAMILY METHODOLOGY VALIDATION ==========\n\n")
print(all_summary, row.names = FALSE)

write.table(all_summary, file.path(FIGURES_DIR, "expanded_families_validation.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ── Topology figure for all new families (at strict threshold) ────────────────
mk_topo <- function(name, vg) {
  ggraph(vg$graph, layout = "sugiyama") +
    geom_edge_link(arrow = arrow(length = unit(1.4, "mm"), type = "closed"),
                   end_cap = circle(2, "mm"),
                   color = "gray60", linewidth = 0.3) +
    geom_node_point(aes(color = node_type, size = n_copies)) +
    scale_color_manual(values = c(shared = "#e6a817",
                                    copy_specific = "#1565C0"),
                       guide = "none") +
    scale_size_continuous(range = c(1, 4), guide = "none") +
    labs(title = sprintf("%s  (%d copies, %d nodes, %d classes%s)",
                          name,
                          length(vg$paths),
                          nrow(vg$nodes),
                          family_equivalence_classes(vg)$n_classes,
                          ifelse(family_equivalence_classes(vg)$fully_resolved,
                                  " ✓", "")) ) +
    theme_graph(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 8.5))
}

new_vgs <- lapply(results, function(r) r$strict)
ref_strict_list <- list()
# also load strict RBMY if it exists
if (file.exists(file.path(DATA_DIR, "rbmy_vg_seq_strict.rds"))) {
  ref_strict_list$RBMY_strict <- readRDS(file.path(DATA_DIR, "rbmy_vg_seq_strict.rds"))
}

all_vgs <- c(ref_list,
              setNames(new_vgs, paste0(names(new_vgs), "_strict")),
              ref_strict_list)
panels <- Map(mk_topo, names(all_vgs), all_vgs)

p_all <- wrap_plots(panels, ncol = 3)
ggsave(file.path(FIGURES_DIR, "expanded_families_topologies.png"),
       p_all, width = 14,
       height = ceiling(length(panels) / 3) * 3.2,
       units = "in", dpi = 180)
message("\nTopologies saved.")
