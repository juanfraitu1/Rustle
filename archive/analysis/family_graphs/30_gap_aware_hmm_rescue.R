suppressPackageStartupMessages({
  library(Biostrings)
  library(pwalign)
  library(Rsamtools)
  library(GenomicRanges)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(tidyr)
})

source("02_build_graphs.R")

FIGURES_DIR <- "figures"
DATA_DIR    <- "data"
GENOME_FA   <- "/mnt/c/Users/jfris/Desktop/GGO.fasta"

vg <- readRDS(file.path(DATA_DIR, "golga_GOLGA8I_vg_seq.rds"))
genome <- FaFile(GENOME_FA); open(genome)

get_canonical_tx <- function(gene_id) {
  gene_txs <- unique(vg$exon_df$tx_id[vg$exon_df$gene_id == gene_id])
  cnt <- sapply(gene_txs, function(tx) sum(vg$exon_df$tx_id == tx))
  gene_txs[which.max(cnt)]
}
extract_copy_seq <- function(gene_id) {
  tx <- get_canonical_tx(gene_id)
  ex <- vg$exon_df[vg$exon_df$tx_id == tx, ]; ex <- ex[order(ex$start), ]
  gr <- GRanges(seqnames = ex$chrom, ranges = IRanges(ex$start, ex$end),
                strand = ex$strand)
  seqs <- scanFa(genome, gr)
  if (as.character(ex$strand[1]) == "-") seqs <- rev(reverseComplement(seqs))
  list(seq = unlist(seqs))
}
copy_ids <- names(vg$paths)
copies   <- setNames(lapply(copy_ids, extract_copy_seq), copy_ids)
labels   <- setNames(sub("^gene-LOC", "", copy_ids), copy_ids)
close(genome)

SOURCE   <- "gene-LOC129527171"
src_seq  <- as.character(copies[[SOURCE]]$seq)

set.seed(11)
ERROR_RATE <- 0.02
mutate_seq <- function(s) {
  ch  <- strsplit(s, "")[[1]]
  alt <- c("A","C","G","T")
  mp  <- sample.int(length(ch), round(length(ch) * ERROR_RATE))
  for (p in mp) ch[p] <- sample(setdiff(alt, ch[p]), 1)
  paste(ch, collapse = "")
}
read_seq <- mutate_seq(src_seq)
read_chars <- strsplit(read_seq, "")[[1]]

# ── Pairwise-align source to each non-source copy, get per-position bases ────
SUB_MAT <- pwalign::nucleotideSubstitutionMatrix(match = 1, mismatch = -1,
                                                  baseOnly = TRUE)
get_copy_bases <- function(copy_seq) {
  aln <- pwalign::pairwiseAlignment(
    pattern = DNAString(src_seq), subject = DNAString(copy_seq),
    type = "global", substitutionMatrix = SUB_MAT,
    gapOpening = -2, gapExtension = -1
  )
  ps <- strsplit(as.character(alignedPattern(aln)), "")[[1]]
  ss <- strsplit(as.character(alignedSubject(aln)), "")[[1]]
  # walk through alignment columns; only keep columns where source has a base
  out <- rep("-", nchar(src_seq))
  src_idx <- cumsum(ps != "-")
  mask <- ps != "-"
  out[src_idx[mask]] <- ss[mask]   # may include "-" where copy has gap
  out
}

cat("Computing per-copy bases at source positions (incl. gaps)...\n")
base_mat <- matrix("-", nrow = length(copy_ids), ncol = nchar(src_seq))
rownames(base_mat) <- copy_ids
base_mat[SOURCE, ] <- strsplit(src_seq, "")[[1]]
for (g in setdiff(copy_ids, SOURCE)) {
  base_mat[g, ] <- get_copy_bases(copies[[g]]$seq)
}

# ── Two scoring schemes ──────────────────────────────────────────────────────
P_MATCH   <- 0.95
LOG_MATCH <- log(P_MATCH)
LOG_MISS  <- log((1 - P_MATCH) / 3)
LOG_GAP   <- log(1e-6)   # "copy has no base here" — strong evidence of absence

score_old <- function(copy_bases) {
  # OLD: skip positions where copy has gap (filter them out)
  is_gap <- copy_bases == "-"
  ok <- !is_gap
  L_ok <- sum(ok)
  M    <- sum(read_chars[ok] == copy_bases[ok])
  X    <- L_ok - M
  list(L_total = L_ok, M = M, X = X, G = sum(is_gap),
       logL = M * LOG_MATCH + X * LOG_MISS)
}
score_new <- function(copy_bases) {
  # NEW: count gaps as evidence of absence
  is_gap <- copy_bases == "-"
  is_match <- (!is_gap) & (read_chars == copy_bases)
  is_miss  <- (!is_gap) & (!is_match)
  M <- sum(is_match); X <- sum(is_miss); G <- sum(is_gap)
  list(L_total = M + X + G, M = M, X = X, G = G,
       logL = M * LOG_MATCH + X * LOG_MISS + G * LOG_GAP)
}

results <- data.frame()
for (g in copy_ids) {
  cb <- base_mat[g, ]
  o <- score_old(cb); n <- score_new(cb)
  results <- bind_rows(results, data.frame(
    copy = g, copy_label = labels[g],
    L_total_old = o$L_total, M_old = o$M, X_old = o$X, G_old = o$G,
    logL_old = o$logL,
    L_total_new = n$L_total, M_new = n$M, X_new = n$X, G_new = n$G,
    logL_new = n$logL,
    is_source = g == SOURCE,
    stringsAsFactors = FALSE
  ))
}
results <- results %>% arrange(desc(logL_old))
results$delta_old <- results$logL_old - max(results$logL_old)
results$delta_new <- results$logL_new - max(results$logL_new)

cat("\nComparison of OLD (filter gaps) vs NEW (gap as absence):\n\n")
print(results[, c("copy_label","is_source","M_old","X_old","G_old","delta_old",
                  "M_new","X_new","G_new","delta_new")],
      row.names = FALSE, digits = 4)

# ── Build figure showing the rescue ──────────────────────────────────────────
plot_df <- results %>%
  arrange(desc(logL_new)) %>%
  mutate(copy_label = factor(copy_label, levels = rev(copy_label)),
         bar_fill = ifelse(is_source, "#2E7D32", "#9E9E9E"))

# Panel A: old scoring (filter gaps) — 171 ties with 071
p_old <- ggplot(plot_df, aes(y = copy_label, x = delta_old, fill = bar_fill)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "gray30") +
  geom_text(aes(label = sprintf("%.0f", delta_old)),
            x = -5, color = "white", fontface = "bold", size = 3, hjust = 1) +
  geom_text(data = plot_df %>% filter(is_source),
            aes(y = copy_label, x = max(plot_df$delta_old) + 5,
                label = "SOURCE"),
            color = "#2E7D32", fontface = "bold", size = 3, hjust = 0) +
  scale_fill_identity() +
  labs(title = "OLD scoring  —  filter out gap positions",
       subtitle = "Discards 271 bp of source's alt 3' end (no homolog in 071) -> 171 and 071 tie",
       x = "DlogL (vs best)", y = NULL) +
  theme_classic(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(color = "gray40", size = 10),
        axis.text.y   = element_text(face = "bold", size = 10))

# Panel B: new scoring (gap as absence)
p_new <- ggplot(plot_df, aes(y = copy_label, x = delta_new, fill = bar_fill)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "gray30") +
  geom_text(aes(label = sprintf("%.0f", delta_new)),
            x = ifelse(plot_df$delta_new < -200,
                       plot_df$delta_new + 80,
                       -80),
            color = "white", fontface = "bold", size = 3, hjust = 1) +
  geom_text(data = plot_df %>% filter(is_source),
            aes(y = copy_label, x = 100,
                label = "SOURCE"),
            color = "#2E7D32", fontface = "bold", size = 3, hjust = 0) +
  scale_fill_identity() +
  labs(title = "NEW scoring  —  gap as evidence of absence  (logP_gap = -13.8 per position)",
       subtitle = "Counts the 271 missing bp of node 8 against copies that lack it -> 171 dominates",
       x = "DlogL (vs best)", y = NULL) +
  theme_classic(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(color = "gray40", size = 10),
        axis.text.y   = element_text(face = "bold", size = 10))

# Panel C: gap counts per copy (so we can SEE why)
p_gaps <- ggplot(plot_df, aes(y = copy_label, x = G_new, fill = bar_fill)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%d gaps", G_new)),
            x = pmax(plot_df$G_new - 5, 5), color = "white",
            fontface = "bold", size = 3, hjust = 1) +
  scale_fill_identity() +
  labs(title = "Gap-position count per copy  (alignment 171 -> copy)",
       subtitle = "Higher count = more positions where the copy has no homologous base to 171's transcript.",
       x = "# gap positions (out of 2789 source positions)",
       y = NULL) +
  theme_classic(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(color = "gray40", size = 10),
        axis.text.y   = element_text(face = "bold", size = 10))

fig <- (p_old | p_new) / p_gaps + plot_layout(heights = c(1, 0.9)) +
  plot_annotation(
    title    = "Gap-as-absence rescues GOLGA8I FLNC_129527171",
    subtitle = paste0(
      "Read R = LOC129527171 transcript (2789 bp, truncated copy with alt 3' end exon = node 8) + 2% error.\n",
      "OLD HMM scoring filters out positions where any copy has a gap, throwing away the diagnostic ",
      "271 bp of node 8.\n",
      "NEW scoring treats gaps as evidence the copy lacks a base at that position. 071 (and others) ",
      "have ~270 gaps where node 8 should be -> massive negative score -> source 171 wins decisively."),
    theme = theme(plot.title    = element_text(size = 14, face = "bold"),
                  plot.subtitle = element_text(size = 10.5, color = "gray40",
                                                lineheight = 1.3))
  )

ggsave(file.path(FIGURES_DIR, "fig_gap_aware_rescue.pdf"), fig,
       width = 16, height = 10, units = "in")
ggsave(file.path(FIGURES_DIR, "fig_gap_aware_rescue.png"), fig,
       width = 16, height = 10, units = "in", dpi = 180)
message("Gap-aware HMM rescue figure saved.")
