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

# ── Pairwise-align each copy to every other copy ─────────────────────────────
SUB_MAT <- pwalign::nucleotideSubstitutionMatrix(match = 1, mismatch = -1,
                                                  baseOnly = TRUE)
get_copy_bases <- function(src_seq, copy_seq) {
  aln <- pwalign::pairwiseAlignment(
    pattern = DNAString(src_seq), subject = DNAString(copy_seq),
    type = "global", substitutionMatrix = SUB_MAT,
    gapOpening = -2, gapExtension = -1
  )
  ps <- strsplit(as.character(alignedPattern(aln)), "")[[1]]
  ss <- strsplit(as.character(alignedSubject(aln)), "")[[1]]
  out <- rep("-", nchar(src_seq))
  src_idx <- cumsum(ps != "-"); mask <- ps != "-"
  out[src_idx[mask]] <- ss[mask]
  out
}

cat("Precomputing per-source-per-copy base matrices (gap-aware) ...\n")
base_mat_list <- list()
for (s in copy_ids) {
  src_seq <- as.character(copies[[s]]$seq)
  base_mat_list[[s]] <- matrix("-", nrow = length(copy_ids), ncol = nchar(src_seq))
  rownames(base_mat_list[[s]]) <- copy_ids
  base_mat_list[[s]][s, ] <- strsplit(src_seq, "")[[1]]
  for (c in setdiff(copy_ids, s)) {
    base_mat_list[[s]][c, ] <- get_copy_bases(src_seq, as.character(copies[[c]]$seq))
  }
}

# ── Gap-aware logL scoring ──────────────────────────────────────────────────
P_MATCH   <- 0.95
LOG_MATCH <- log(P_MATCH)
LOG_MISS  <- log((1 - P_MATCH) / 3)
LOG_GAP   <- log(1e-6)

score_gap_aware <- function(read_chars, copy_bases) {
  is_gap   <- copy_bases == "-"
  is_match <- (!is_gap) & (read_chars == copy_bases)
  is_miss  <- (!is_gap) & (!is_match)
  sum(is_match) * LOG_MATCH +
    sum(is_miss)  * LOG_MISS  +
    sum(is_gap)   * LOG_GAP
}

# ── Simulate a library with known abundances ─────────────────────────────────
set.seed(42)
true_theta <- c(0.25, 0.10, 0.18, 0.12, 0.05, 0.20, 0.10)
names(true_theta) <- copy_ids
cat("\nTrue abundances:\n"); print(round(true_theta, 3))

N <- 100
sources <- sample(copy_ids, N, replace = TRUE, prob = true_theta)

ERROR_RATE <- 0.02
mutate_seq <- function(s) {
  ch  <- strsplit(s, "")[[1]]
  alt <- c("A","C","G","T")
  mp  <- sample.int(length(ch), round(length(ch) * ERROR_RATE))
  for (p in mp) ch[p] <- sample(setdiff(alt, ch[p]), 1)
  paste(ch, collapse = "")
}

cat("\nSimulating", N, "FLNC reads...\n")
read_chars_list <- lapply(seq_len(N), function(i) {
  strsplit(mutate_seq(as.character(copies[[sources[i]]]$seq)), "")[[1]]
})

# ── For each read, compute logL vs each copy using gap-aware scoring ─────────
# FLNC-aware scoring: add a length-consistency penalty.
# For an FLNC read R from copy C, |R| should equal |C|'s transcript length.
# Penalty: -|R-len - C-len| / LAMBDA  (LAMBDA controls strictness)
LAMBDA_LEN <- 20   # ~50 bp difference = 2.5 log-units
copy_lens  <- sapply(copies, function(c) length(c$seq))
names(copy_lens) <- copy_ids

cat("Scoring each read against each copy (gap-aware + FLNC length)...\n")
logL_matrix     <- matrix(NA, nrow = N, ncol = length(copy_ids))
logL_no_lenpen  <- matrix(NA, nrow = N, ncol = length(copy_ids))
colnames(logL_matrix)    <- copy_ids
colnames(logL_no_lenpen) <- copy_ids
read_lens <- sapply(read_chars_list, length)
for (i in seq_len(N)) {
  src <- sources[i]
  rc  <- read_chars_list[[i]]
  rl  <- read_lens[i]
  for (j in seq_along(copy_ids)) {
    cb <- base_mat_list[[src]][copy_ids[j], ]
    base_logL <- score_gap_aware(rc, cb)
    len_pen   <- -abs(rl - copy_lens[copy_ids[j]]) / LAMBDA_LEN
    logL_no_lenpen[i, j] <- base_logL
    logL_matrix[i, j]    <- base_logL + len_pen
  }
}

# ── Identify uniquely-assignable reads (large DlogL margin > threshold) ──────
DLOG_UNIQUE <- 50   # margin to declare a read "uniquely assignable"
best_idx <- apply(logL_matrix, 1, which.max)
second_idx <- apply(logL_matrix, 1, function(r) order(-r)[2])
margins <- vapply(seq_len(N), function(i) {
  logL_matrix[i, best_idx[i]] - logL_matrix[i, second_idx[i]]
}, numeric(1))
unique_reads <- which(margins >= DLOG_UNIQUE)
ambig_reads  <- which(margins <  DLOG_UNIQUE)
cat("\nUniquely-assignable reads: ", length(unique_reads),
    "  (margin >= ", DLOG_UNIQUE, ")\n", sep = "")
cat("Ambiguous reads          : ", length(ambig_reads), "\n", sep = "")

# Bootstrap theta from uniquely-assigned reads
unique_counts <- table(factor(copy_ids[best_idx[unique_reads]], levels = copy_ids))
theta_init <- as.numeric(unique_counts) / sum(unique_counts)
theta_init[theta_init == 0] <- 1e-6
theta_init <- theta_init / sum(theta_init)
names(theta_init) <- copy_ids
cat("\nBootstrap theta from unique reads:\n"); print(round(theta_init, 3))

# ── EM iterations ────────────────────────────────────────────────────────────
N_ITER <- 12
theta_trajectory <- matrix(0, nrow = N_ITER + 1, ncol = length(copy_ids))
colnames(theta_trajectory) <- copy_ids
theta_trajectory[1, ] <- theta_init

theta <- theta_init
for (iter in seq_len(N_ITER)) {
  log_theta <- log(theta)
  # E-step: posterior P(copy | read) ∝ θ * P(read | copy)
  log_post <- sweep(logL_matrix, 2, log_theta, "+")
  log_post <- log_post - apply(log_post, 1, max)
  post <- exp(log_post)
  post <- post / rowSums(post)
  # M-step
  theta <- colMeans(post)
  theta[theta < 1e-8] <- 1e-8
  theta <- theta / sum(theta)
  theta_trajectory[iter + 1L, ] <- theta
}

cat("\nFinal theta after EM:\n"); print(round(theta, 3))
cat("True theta:           "); print(round(true_theta, 3))

# ── Assemble figure ──────────────────────────────────────────────────────────
traj_df <- as.data.frame(theta_trajectory)
traj_df$iteration <- 0:N_ITER
traj_long <- traj_df %>%
  pivot_longer(cols = -iteration, names_to = "copy", values_to = "theta") %>%
  mutate(copy_label = labels[copy])

copy_pal <- setNames(
  c("#1565C0","#B71C1C","#2E7D32","#6A1B9A","#EF6C00","#00838F","#5D4037"),
  copy_ids
)
traj_long$copy_label <- factor(traj_long$copy_label, levels = labels[copy_ids])
traj_long$true_val <- true_theta[traj_long$copy]

# Panel A: theta trajectory over iterations
p_traj <- ggplot(traj_long, aes(x = iteration, y = theta, color = copy_label,
                                  group = copy_label)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.0) +
  geom_hline(data = data.frame(copy_label = factor(labels[copy_ids],
                                                     levels = labels[copy_ids]),
                                true_val = true_theta[copy_ids]),
             aes(yintercept = true_val, color = copy_label),
             linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  scale_color_manual(values = setNames(unname(copy_pal), labels[copy_ids]),
                     name = NULL) +
  scale_x_continuous(breaks = 0:N_ITER) +
  labs(title    = "EM trajectory of theta over iterations",
       subtitle = "Solid lines = EM estimate, dashed lines = true abundance (target). Iteration 0 = bootstrap from uniquely-assigned reads only.",
       x = "EM iteration",
       y = "theta (estimated abundance fraction)") +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(color = "gray40", size = 10),
        legend.position = "right",
        panel.grid.minor = element_blank())

# Panel B: bar chart comparing true vs final theta
final_df <- data.frame(
  copy = copy_ids,
  copy_label = labels[copy_ids],
  true_theta = unname(true_theta[copy_ids]),
  final_theta = unname(theta[copy_ids]),
  initial_theta = unname(theta_init[copy_ids])
) %>%
  pivot_longer(cols = c(true_theta, initial_theta, final_theta),
                names_to = "kind", values_to = "value") %>%
  mutate(kind = factor(kind, levels = c("initial_theta","final_theta","true_theta"),
                        labels = c("theta_0 (bootstrap)","theta_final (EM)","theta_true")),
         copy_label = factor(copy_label, levels = labels[copy_ids]))

p_bars <- ggplot(final_df, aes(x = copy_label, y = value, fill = kind)) +
  geom_col(position = "dodge", color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", value)),
            position = position_dodge(width = 0.9), vjust = -0.3,
            size = 2.5, color = "gray20") +
  scale_fill_manual(values = c("theta_0 (bootstrap)" = "#90A4AE",
                                "theta_final (EM)"    = "#1565C0",
                                "theta_true"          = "#2E7D32"),
                    name = NULL) +
  labs(title    = "Final EM estimate vs true abundance vs bootstrap initial",
       subtitle = sprintf("EM with N=%d simulated FLNC reads; %d uniquely-assignable, %d ambiguous reads.",
                            N, length(unique_reads), length(ambig_reads)),
       x = NULL, y = "theta") +
  theme_classic(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(color = "gray40", size = 10),
        axis.text.x   = element_text(face = "bold", size = 10),
        legend.position = "top")

fig <- p_traj / p_bars + plot_layout(heights = c(1.0, 1.0)) +
  plot_annotation(
    title    = "EM bootstrap for GOLGA8I abundance estimation",
    subtitle = paste0(
      "Library of ", N, " simulated FLNC reads from 7 GOLGA8I copies (2% error). ",
      "Gap-aware HMM scoring computes P(read | copy_i). ",
      "Bootstrap from uniquely-assignable reads (DlogL margin >= 50) gives theta_0; ",
      "EM iterates over ambiguous reads using current theta as prior."),
    theme = theme(plot.title    = element_text(size = 14, face = "bold"),
                  plot.subtitle = element_text(size = 10.5, color = "gray40",
                                                lineheight = 1.35))
  )

ggsave(file.path(FIGURES_DIR, "fig_em_bootstrap.pdf"), fig,
       width = 14, height = 10, units = "in")
ggsave(file.path(FIGURES_DIR, "fig_em_bootstrap.png"), fig,
       width = 14, height = 10, units = "in", dpi = 180)
message("EM bootstrap figure saved.")
