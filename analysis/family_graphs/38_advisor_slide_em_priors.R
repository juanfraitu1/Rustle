suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(tidyr)
})

FIGURES_DIR <- "figures"

# ── Toy example representative of FLNC: 3 copies, length 15, 4 SNPs ──────────
# Reads are FULL-LENGTH (cover the entire transcript), as IsoSeq FLNC reads do.
set.seed(7)
L_TX  <- 15L
bases <- c("A","T","C","G")
snps  <- c(3L, 7L, 11L, 14L)

# Build copy sequences sharing a backbone but differing at the SNPs
shared <- sample(bases, L_TX, replace = TRUE)
seq_A <- shared; seq_B <- shared; seq_C <- shared
# Make A, B, C differ at each SNP (all three distinct)
for (p in snps) {
  trio <- sample(bases, 3, replace = FALSE)
  seq_A[p] <- trio[1]; seq_B[p] <- trio[2]; seq_C[p] <- trio[3]
}
copies <- list(A = seq_A, B = seq_B, C = seq_C)
cat("Copy A:", paste(seq_A, collapse=""), "\n")
cat("Copy B:", paste(seq_B, collapse=""), "\n")
cat("Copy C:", paste(seq_C, collapse=""), "\n")
cat("SNP positions:", snps, "\n\n")

# Per-copy HMM emission
p_match <- 0.95; p_other <- (1 - p_match)/3
mk_emit <- function(seq) {
  m <- matrix(p_other, nrow = length(seq), ncol = 4,
              dimnames = list(seq_along(seq), bases))
  for (i in seq_along(seq)) m[i, seq[i]] <- p_match
  m
}
emit <- list(A = mk_emit(seq_A), B = mk_emit(seq_B), C = mk_emit(seq_C))

# ── Simulate 40 full-length FLNC reads with ELEVATED error rate ──────────────
# Real IsoSeq HiFi is ~0.5% error; we use 15% here to expose, in this toy,
# the regime where errors at SNPs create per-read ambiguity. (Same regime
# you would see for real data on copies sharing very high identity at >99%.)
ERROR_RATE <- 0.15
TRUE_THETA <- c(A = 0.50, B = 0.30, C = 0.20)
N_READS <- 40
add_err <- function(seq, rate) {
  ch <- seq
  for (i in seq_along(ch)) {
    if (runif(1) < rate) ch[i] <- sample(setdiff(bases, ch[i]), 1)
  }
  ch
}
reads <- vector("list", N_READS)
for (i in seq_len(N_READS)) {
  src <- sample(names(copies), 1, prob = TRUE_THETA)
  read_seq <- add_err(copies[[src]], ERROR_RATE)
  reads[[i]] <- list(read_id = sprintf("r%02d", i),
                     src = src,
                     bases = read_seq,
                     seq = paste(read_seq, collapse = ""))
}

# ── Score each read against each copy's per-copy HMM (FULL-LENGTH) ──────────
score_full <- function(read_bases, emit_mat) {
  cols <- match(read_bases, bases)
  sum(log(emit_mat[cbind(seq_along(read_bases), cols)]))
}
logL_mat <- t(sapply(reads, function(r) {
  c(A = score_full(r$bases, emit$A),
    B = score_full(r$bases, emit$B),
    C = score_full(r$bases, emit$C))
}))
rownames(logL_mat) <- sapply(reads, function(r) r$read_id)

# Compute hard call + margin
best_idx <- apply(logL_mat, 1, which.max)
sort_logL <- t(apply(logL_mat, 1, function(r) sort(r, decreasing = TRUE)))
margins <- sort_logL[,1] - sort_logL[,2]
hard_call <- colnames(logL_mat)[best_idx]
DLOG_THRESH <- 5
clean <- margins >= DLOG_THRESH

read_summary <- data.frame(
  read = sapply(reads, function(r) r$read_id),
  seq  = sapply(reads, function(r) r$seq),
  src  = sapply(reads, function(r) r$src),
  best = hard_call,
  margin = margins,
  clean = clean,
  stringsAsFactors = FALSE
)
cat("Hard calls (HMM only, margin >=", DLOG_THRESH, "):\n")
cat("  Clean reads: ", sum(clean), " / ", N_READS, "\n")
cat("  Ambiguous reads: ", sum(!clean), "\n")
print(table(predicted = read_summary$best, true_source = read_summary$src))

# Bootstrap theta from clean reads
clean_calls <- read_summary$best[read_summary$clean]
counts <- as.numeric(table(factor(clean_calls, levels = c("A","B","C"))))
theta_init <- counts / max(1, sum(counts)); theta_init[theta_init == 0] <- 1e-6
theta_init <- theta_init / sum(theta_init); names(theta_init) <- c("A","B","C")
cat("\nBootstrap theta from clean reads:\n"); print(round(theta_init, 3))

# EM
N_ITER <- 10
theta_traj <- matrix(0, nrow = N_ITER + 1, ncol = 3,
                      dimnames = list(NULL, c("A","B","C")))
theta_traj[1, ] <- theta_init
theta <- theta_init
for (iter in seq_len(N_ITER)) {
  log_theta <- log(theta)
  log_post <- sweep(logL_mat, 2, log_theta, "+")
  log_post <- log_post - apply(log_post, 1, max)
  post <- exp(log_post); post <- post / rowSums(post)
  theta <- colMeans(post)
  theta[theta < 1e-8] <- 1e-8; theta <- theta / sum(theta)
  names(theta) <- c("A","B","C")
  theta_traj[iter + 1, ] <- theta
}
cat("\nFinal theta after EM:\n"); print(round(theta, 3))
cat("True theta:\n"); print(round(TRUE_THETA, 3))

# Per-read final posterior
log_theta_f <- log(theta)
log_post_f  <- sweep(logL_mat, 2, log_theta_f, "+")
log_post_f  <- log_post_f - apply(log_post_f, 1, max)
post_f      <- exp(log_post_f); post_f <- post_f / rowSums(post_f)
read_summary$P_A <- post_f[,1]
read_summary$P_B <- post_f[,2]
read_summary$P_C <- post_f[,3]

# ── Figure ────────────────────────────────────────────────────────────────────

# Panel A: header explaining the setup (FLNC-realistic)
p_header <- ggplot() +
  annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = 10,
           fill = "#E3F2FD", color = "#0D47A1", linewidth = 0.5) +
  annotate("text", x = 0.2, y = 9.0,
           label = "Per-copy HMMs + EM priors  —  toy example representative of FLNC",
           hjust = 0, fontface = "bold", size = 4.6, color = "#0D47A1") +
  annotate("text", x = 0.2, y = 7.6,
           label = sprintf("3 copies of length %d, %d SNP positions. Reads are FULL-LENGTH (cover the entire transcript), 8%% per-base error.",
                            L_TX, length(snps)),
           hjust = 0, size = 3.4, color = "gray15") +
  annotate("text", x = 0.2, y = 6.4,
           label = sprintf("Copy A: %s", paste(seq_A, collapse = "")),
           hjust = 0, family = "mono", size = 3.4, color = "#1565C0",
           fontface = "bold") +
  annotate("text", x = 0.2, y = 5.6,
           label = sprintf("Copy B: %s", paste(seq_B, collapse = "")),
           hjust = 0, family = "mono", size = 3.4, color = "#B71C1C",
           fontface = "bold") +
  annotate("text", x = 0.2, y = 4.8,
           label = sprintf("Copy C: %s", paste(seq_C, collapse = "")),
           hjust = 0, family = "mono", size = 3.4, color = "#2E7D32",
           fontface = "bold") +
  annotate("text", x = 0.2, y = 3.6,
           label = sprintf("SNPs at positions  %s.   Other positions are shared across copies.",
                            paste(snps, collapse = ", ")),
           hjust = 0, size = 3.2, color = "gray35", fontface = "italic") +
  annotate("text", x = 0.2, y = 2.3,
           label = sprintf("True abundances: theta_A = %.2f,  theta_B = %.2f,  theta_C = %.2f",
                            TRUE_THETA["A"], TRUE_THETA["B"], TRUE_THETA["C"]),
           hjust = 0, size = 3.4, color = "gray15", fontface = "bold") +
  annotate("text", x = 0.2, y = 1.0,
           label = sprintf("Library: %d FLNC reads simulated.", N_READS),
           hjust = 0, size = 3.4, color = "gray15") +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 10)) +
  theme_void()

# Panel B: heatmap of reads (rows) showing bases (cols), errors at SNPs highlighted
read_long <- bind_rows(lapply(seq_along(reads), function(i) {
  r <- reads[[i]]
  src_seq <- copies[[r$src]]
  data.frame(read = factor(r$read_id, levels = sapply(reads, function(rr) rr$read_id)),
             pos  = seq_along(r$bases),
             base = r$bases,
             src  = r$src,
             is_snp = seq_along(r$bases) %in% snps,
             is_error = r$bases != src_seq,
             stringsAsFactors = FALSE)
}))
nuc_pal <- c(A="#43A047", T="#1E88E5", C="#FB8C00", G="#E53935")

p_reads <- ggplot(read_long, aes(x = pos, y = read, fill = base)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = base), color = "white", fontface = "bold", size = 2.6) +
  # mark errors with a black square
  geom_tile(data = filter(read_long, is_error), aes(x = pos, y = read),
            fill = NA, color = "black", linewidth = 0.9) +
  scale_fill_manual(values = nuc_pal, guide = "none") +
  geom_vline(xintercept = snps - 0.5, color = "purple", linetype = "dashed",
             alpha = 0.5) +
  geom_vline(xintercept = snps + 0.5, color = "purple", linetype = "dashed",
             alpha = 0.5) +
  scale_x_continuous(breaks = 1:L_TX, labels = 1:L_TX) +
  scale_y_discrete(limits = rev) +
  labs(title    = sprintf("%d full-length FLNC reads — base-by-base view",
                            N_READS),
       subtitle = "Purple dashed lines bracket the 4 SNP positions. Black border = position where the read differs from its true source (= error).",
       x = "Position", y = NULL) +
  theme_minimal(base_size = 9) +
  theme(panel.grid    = element_blank(),
        plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(color = "gray40", size = 9),
        axis.text.x   = element_text(size = 7),
        axis.text.y   = element_text(size = 7, face = "bold"))

# Panel C: per-read hard call + final EM posterior bar
src_pal <- c(A="#1565C0", B="#B71C1C", C="#2E7D32")
post_long <- read_summary %>%
  select(read, src, best, margin, clean, P_A, P_B, P_C) %>%
  pivot_longer(cols = starts_with("P_"), names_to = "copy", values_to = "post") %>%
  mutate(copy = sub("P_", "", copy),
         read = factor(read, levels = read_summary$read))

p_post <- ggplot(post_long, aes(y = read, x = post, fill = copy)) +
  geom_col(position = "stack", width = 0.7, color = "white", linewidth = 0.2) +
  # annotation: true source on the right
  geom_text(data = read_summary %>%
              mutate(read = factor(read, levels = read_summary$read)),
            aes(y = read, x = 1.05,
                label = sprintf("true=%s  hard=%s  m=%.0f%s",
                                 src, best, margin,
                                 ifelse(clean, "", " (ambig)"))),
            inherit.aes = FALSE,
            hjust = 0, size = 2.3,
            color = ifelse(read_summary$clean, "gray45", "#B71C1C"),
            fontface = ifelse(read_summary$clean, "plain", "bold")) +
  scale_fill_manual(values = src_pal, name = "Copy") +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(limits = c(0, 1.8), breaks = c(0, 0.5, 1)) +
  labs(title    = "Per-read posterior after EM",
       subtitle = "Stacked bar = P(A|r), P(B|r), P(C|r). Right text = true source, HMM hard call, DlogL margin.",
       x = "P(copy | read)", y = NULL) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(color = "gray40", size = 9),
        axis.text.y   = element_text(size = 7, face = "bold"))

# Panel D: theta trajectory + bar comparison
traj_df <- as.data.frame(theta_traj); traj_df$iter <- 0:N_ITER
traj_long <- pivot_longer(traj_df, cols = -iter, names_to = "copy",
                           values_to = "theta")
p_traj <- ggplot(traj_long, aes(x = iter, y = theta, color = copy, group = copy)) +
  geom_line(linewidth = 1.4) +
  geom_point(size = 2.8) +
  geom_hline(yintercept = TRUE_THETA["A"], linetype = "dashed",
             color = "#1565C0", alpha = 0.5) +
  geom_hline(yintercept = TRUE_THETA["B"], linetype = "dashed",
             color = "#B71C1C", alpha = 0.5) +
  geom_hline(yintercept = TRUE_THETA["C"], linetype = "dashed",
             color = "#2E7D32", alpha = 0.5) +
  scale_color_manual(values = src_pal, name = "Copy") +
  scale_x_continuous(breaks = 0:N_ITER) +
  scale_y_continuous(limits = c(0, 0.75)) +
  labs(title    = "EM trajectory of theta",
       subtitle = sprintf("Iter 0 = bootstrap (%d / %d clean reads). Dashed = true theta. Converges quickly.",
                            sum(clean), N_READS),
       x = "EM iteration", y = "theta") +
  theme_classic(base_size = 10) +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(color = "gray40", size = 9))

# Final theta vs true bar chart
final_df <- data.frame(
  copy = c("A","B","C"),
  true = unname(TRUE_THETA),
  estimated = unname(theta)
) %>% pivot_longer(cols = c(true, estimated),
                    names_to = "kind", values_to = "value")

p_bars <- ggplot(final_df, aes(x = copy, y = value, fill = kind)) +
  geom_col(position = "dodge", color = "white", linewidth = 0.3, width = 0.75) +
  geom_text(aes(label = sprintf("%.2f", value)),
            position = position_dodge(width = 0.75),
            vjust = -0.3, size = 3.2, color = "gray20") +
  scale_fill_manual(values = c(true = "#2E7D32", estimated = "#1565C0"),
                    name = NULL,
                    labels = c(true = "true theta",
                                estimated = "theta after EM")) +
  labs(title    = "True vs estimated theta",
       x = NULL, y = NULL) +
  theme_classic(base_size = 10) +
  theme(plot.title  = element_text(face = "bold", size = 11),
        axis.text.x = element_text(face = "bold", size = 11),
        legend.position = "top")

# Layout
top_row <- p_header
mid_row <- p_reads | p_post
bot_row <- p_traj | p_bars
fig <- top_row / mid_row / bot_row +
  plot_layout(heights = c(0.55, 1.6, 0.85)) +
  plot_annotation(
    title    = "Per-copy profile HMMs + EM priors on full-length FLNC-like reads",
    subtitle = "Toy example sized for FLNC (full-length reads, multi-SNP discrimination). HMM alone hard-calls the clear-margin reads; EM uses those calls to soft-assign the rest.",
    theme = theme(plot.title    = element_text(size = 16, face = "bold"),
                  plot.subtitle = element_text(size = 10.5, color = "gray40"))
  )

ggsave(file.path(FIGURES_DIR, "fig_advisor_slide_em_priors.pdf"), fig,
       width = 18, height = 13, units = "in")
ggsave(file.path(FIGURES_DIR, "fig_advisor_slide_em_priors.png"), fig,
       width = 18, height = 13, units = "in", dpi = 170)
message("Slide C (FLNC-representative) saved.")
