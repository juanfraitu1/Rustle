suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(patchwork)
})

source("02_build_graphs.R")

DATA_DIR    <- "data"
FIGURES_DIR <- "figures"

# Load all families we have seq-VGs for
fam_files <- list(
  RBMY     = "rbmy_vg_seq.rds",
  AMY      = "amy_vg_seq.rds",
  GOLGA8I  = "golga_GOLGA8I_vg_seq.rds",
  GOLGA8NL = "golga_GOLGA8NL_vg_seq.rds",
  GOLGA6C  = "golga_GOLGA6C_vg_seq.rds"
)
vgs <- lapply(fam_files, function(f) readRDS(file.path(DATA_DIR, f)))

# ── Helper: classify the difference between two paths ────────────────────────
# Returns a character vector of mechanisms present in the pair-wise diff.
classify_pair <- function(pa, pb) {
  # Normalise paths to character
  pa <- as.character(pa)
  pb <- as.character(pb)
  mechs <- c()

  if (identical(pa, pb)) return("IDENTICAL")

  # Alt start: different first node
  if (pa[1] != pb[1]) mechs <- c(mechs, "alt_start")

  # Alt end: different last node
  if (pa[length(pa)] != pb[length(pb)]) mechs <- c(mechs, "alt_end")

  # Truncation: one path is a strict subsequence of the other (much shorter)
  if (abs(length(pa) - length(pb)) >= 3L) mechs <- c(mechs, "truncation")

  # Internal divergence: paths share endpoints but differ in middle
  if (pa[1] == pb[1] && pa[length(pa)] == pb[length(pb)]) {
    # find any internal diff
    common_prefix <- 0L
    while (common_prefix < min(length(pa), length(pb)) &&
           pa[common_prefix + 1L] == pb[common_prefix + 1L]) {
      common_prefix <- common_prefix + 1L
    }
    if (common_prefix < min(length(pa), length(pb))) {
      mechs <- c(mechs, "internal_bubble")
    }
  } else if (length(intersect(pa, pb)) >= 2L &&
             !("alt_start" %in% mechs) && !("alt_end" %in% mechs) &&
             length(setdiff(pa, pb)) > 0L) {
    mechs <- c(mechs, "internal_bubble")
  }

  if (length(mechs) == 0L) mechs <- "other"
  mechs
}

# ── Per-family analysis ───────────────────────────────────────────────────────
analyse_family <- function(name, vg) {
  paths <- vg$paths
  n_copies <- length(paths)

  # Distinctness check
  path_strings <- sapply(paths, function(p) paste(as.character(p), collapse = ","))
  n_distinct   <- length(unique(path_strings))
  dup_pairs    <- if (n_distinct < n_copies) {
    dups <- path_strings[duplicated(path_strings) |
                          duplicated(path_strings, fromLast = TRUE)]
    paste(names(dups), collapse = ", ")
  } else { "" }

  # Pairwise mechanism classification
  copy_ids <- names(paths)
  pair_records <- list()
  for (i in seq_len(n_copies - 1L)) {
    for (j in (i + 1L):n_copies) {
      mechs <- classify_pair(paths[[i]], paths[[j]])
      pair_records[[length(pair_records) + 1L]] <- data.frame(
        family = name,
        copy_a = copy_ids[i],
        copy_b = copy_ids[j],
        mechanism = mechs,
        stringsAsFactors = FALSE
      )
    }
  }
  pair_df <- do.call(rbind, pair_records)

  summary_row <- data.frame(
    family            = name,
    n_copies          = n_copies,
    n_distinct_paths  = n_distinct,
    duplicates        = dup_pairs,
    pairs_with_alt_start = sum(pair_df$mechanism == "alt_start" |
                                 (duplicated(paste(pair_df$copy_a, pair_df$copy_b)) &
                                    pair_df$mechanism == "alt_start")),
    stringsAsFactors  = FALSE
  )

  list(summary = summary_row, pairs = pair_df)
}

results <- lapply(names(vgs), function(nm) analyse_family(nm, vgs[[nm]]))
names(results) <- names(vgs)

# ── Build summary table ──────────────────────────────────────────────────────
mech_levels <- c("alt_start","alt_end","internal_bubble","truncation","other")

mech_counts <- do.call(rbind, lapply(results, function(r) {
  d <- r$pairs
  n_total_pairs <- length(unique(paste(d$copy_a, d$copy_b)))
  out <- data.frame(family = d$family[1],
                    n_copies = NA, n_pairs = n_total_pairs,
                    stringsAsFactors = FALSE)
  for (m in mech_levels) {
    pairs_with_m <- unique(paste(d$copy_a[d$mechanism == m],
                                  d$copy_b[d$mechanism == m]))
    out[[m]] <- length(pairs_with_m)
  }
  out
}))
mech_counts$n_copies <- sapply(results, function(r) r$summary$n_copies)
mech_counts$all_distinct <- sapply(results,
  function(r) r$summary$n_distinct_paths == r$summary$n_copies)

cat("\n=========== PATH-BASED COPY DEFINITION VALIDATION ===========\n\n")
print(mech_counts[, c("family","n_copies","n_pairs","all_distinct",
                      "alt_start","alt_end","internal_bubble","truncation")],
      row.names = FALSE)

# Per-family verbose breakdown
for (nm in names(results)) {
  cat("\n", nm, " ---------\n", sep = "")
  cat("  copies      : ", results[[nm]]$summary$n_copies, "\n")
  cat("  distinct    : ", results[[nm]]$summary$n_distinct_paths, "\n")
  if (nzchar(results[[nm]]$summary$duplicates)) {
    cat("  DUPLICATES  : ", results[[nm]]$summary$duplicates, "\n")
  }
  # which mechanisms dominate?
  mech_freq <- table(results[[nm]]$pairs$mechanism)
  cat("  mechanism distribution (pair-mechanisms, multi-count if pair has >1):\n")
  for (m in names(mech_freq)) cat("    ", m, ": ", mech_freq[[m]], "\n")
}

# ── Heatmap-style figure: stacked bar per family showing mechanism mix ────────
plot_df <- mech_counts %>%
  select(family, n_pairs, alt_start, alt_end, internal_bubble, truncation, other) %>%
  pivot_longer(cols = c(alt_start, alt_end, internal_bubble, truncation, other),
                names_to = "mechanism", values_to = "n_pairs_with") %>%
  mutate(family = factor(family, levels = names(vgs)),
         mechanism = factor(mechanism,
                             levels = c("alt_start","alt_end",
                                         "internal_bubble","truncation","other")))

mech_pal <- c(alt_start = "#1565C0",
              alt_end   = "#B71C1C",
              internal_bubble = "#7B1FA2",
              truncation = "#EF6C00",
              other      = "#9E9E9E")

p_mech <- ggplot(plot_df, aes(x = family, y = n_pairs_with, fill = mechanism)) +
  geom_col(width = 0.6, color = "white") +
  scale_fill_manual(values = mech_pal,
                     labels = c(alt_start = "alt start",
                                 alt_end   = "alt end",
                                 internal_bubble = "internal bubble",
                                 truncation = "truncation",
                                 other = "other")) +
  geom_text(aes(label = ifelse(n_pairs_with > 0, n_pairs_with, "")),
            position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 3.4) +
  labs(title    = "Mechanisms distinguishing copy pairs in each family",
       subtitle = "y = # of copy pairs in which that mechanism contributes (a pair can have multiple)",
       x = NULL, y = "Copy pairs with mechanism",
       fill = NULL) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = "gray40", size = 9.5),
        legend.position = "right")

ggsave(file.path(FIGURES_DIR, "fig_path_definition_mechanisms.png"), p_mech,
       width = 11, height = 5.5, units = "in", dpi = 180)
write.table(mech_counts, file.path(FIGURES_DIR, "path_definition_validation.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Mechanism figure saved.")
