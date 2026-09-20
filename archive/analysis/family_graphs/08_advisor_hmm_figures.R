suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(tidyr)
})

FIGURES_DIR <- "figures"
dir.create(FIGURES_DIR, showWarnings = FALSE)

# ── Toy example ───────────────────────────────────────────────────────────────
# 6-position exon; SNPs at positions 3 (C/G) and 5 (A/T)
pos   <- 1:6
seq_A <- c("A","T","C","G","A","T")
seq_B <- c("A","T","G","G","T","T")
snps  <- which(seq_A != seq_B)           # positions 3 and 5

base_pal <- c(A="#43A047", T="#1E88E5", C="#FB8C00", G="#E53935")
SLIDE_BASE <- 13   # base_size for slides

# ── Emission matrix helpers ────────────────────────────────────────────────────
p_match <- 0.95
p_other <- (1 - p_match) / 3
bases   <- c("A","T","C","G")

mk_emit <- function(seq) {
  m <- matrix(p_other, nrow = length(seq), ncol = 4,
              dimnames = list(seq_along(seq), bases))
  for (i in seq_along(seq)) m[i, seq[i]] <- p_match
  m
}
emit_A <- mk_emit(seq_A)
emit_B <- mk_emit(seq_B)

log_lik_read <- function(read_bases, emit, start_pos) {
  rows <- seq_along(read_bases) + start_pos - 1L
  sum(log(emit[cbind(rows, match(read_bases, bases))]))
}

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — The Assignment Problem
# ═══════════════════════════════════════════════════════════════════════════════
#   Top:    Two paralog sequences as colored tiles; SNP positions boxed
#   Bottom: Four reads as bars; annotated "assignable" or "ambiguous"
#           + a small edit-distance table showing why ambiguous reads fail

# ── Sequence tile data ────────────────────────────────────────────────────────
seq_df <- bind_rows(
  data.frame(copy = "Copy A", pos = pos, base = seq_A, is_snp = pos %in% snps),
  data.frame(copy = "Copy B", pos = pos, base = seq_B, is_snp = pos %in% snps)
)

p_seqs <- ggplot(seq_df, aes(x = pos, y = copy)) +
  geom_tile(aes(fill = base), color = "white", linewidth = 2,
            width = 0.88, height = 0.70) +
  geom_tile(data = filter(seq_df, is_snp),
            fill = NA, color = "black", linewidth = 2.2,
            width = 0.88, height = 0.70) +
  geom_text(aes(label = base), color = "white", fontface = "bold", size = 6) +
  scale_fill_manual(values = base_pal) +
  scale_x_continuous(breaks = 1:6, labels = paste0("pos ", 1:6),
                     expand = c(0.08, 0)) +
  labs(title = "Two paralog sequences  (SNPs boxed in black)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = SLIDE_BASE) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 11),
        plot.title = element_text(face = "bold", size = 14))

# ── Read bar data ─────────────────────────────────────────────────────────────
# R1: pos 1-2  (no SNP) → ambiguous
# R2: pos 2-4  (covers SNP at pos 3, base = C → Copy A) → assignable
# R3: pos 4-6  (covers SNP at pos 5, base = A → Copy A) → assignable
# R4: pos 4-6  (covers SNP at pos 5, base = T → Copy B) → assignable
# (R4 has the Copy B base at pos 5)
reads_meta <- data.frame(
  read_id   = c("R1 (ambiguous)", "R2 (assignable)", "R3 (assignable)", "R4 (assignable)"),
  start     = c(1, 2, 4, 4),
  end       = c(2, 4, 6, 6),
  src       = c("ambiguous", "Copy A", "Copy A", "Copy B"),
  # which SNP does this read cover, if any?
  snp_label = c("no diagnostic SNP", "SNP at pos 3: 'C' -> Copy A",
                "SNP at pos 5: 'A' -> Copy A", "SNP at pos 5: 'T' -> Copy B")
)
reads_meta$read_id <- factor(reads_meta$read_id, levels = rev(reads_meta$read_id))

assign_pal <- c("Copy A" = "#1565C0", "Copy B" = "#B71C1C", "ambiguous" = "#757575")

p_reads <- ggplot(reads_meta) +
  # shaded region covered by read
  geom_rect(aes(xmin = start - 0.4, xmax = end + 0.4,
                ymin = as.integer(read_id) - 0.38,
                ymax = as.integer(read_id) + 0.38,
                fill = src), alpha = 0.25) +
  # bold border
  geom_rect(aes(xmin = start - 0.4, xmax = end + 0.4,
                ymin = as.integer(read_id) - 0.38,
                ymax = as.integer(read_id) + 0.38,
                color = src), fill = NA, linewidth = 1.2) +
  # SNP annotation
  geom_text(aes(x = 6.65, y = as.integer(read_id), label = snp_label,
                color = src),
            hjust = 0, size = 3.8, fontface = "italic") +
  # vertical SNP position lines
  geom_vline(xintercept = snps, linetype = "dashed",
             color = "black", alpha = 0.35, linewidth = 0.7) +
  scale_fill_manual(values = assign_pal, name = NULL) +
  scale_color_manual(values = assign_pal, name = NULL) +
  scale_x_continuous(breaks = 1:6, labels = paste0("pos ", 1:6),
                     expand = c(0.08, 0)) +
  scale_y_continuous(breaks = seq_along(levels(reads_meta$read_id)),
                     labels = levels(reads_meta$read_id)) +
  coord_cartesian(xlim = c(0.5, 12), clip = "off") +
  labs(title = "Reads: coverage determines whether they're assignable",
       x = NULL, y = NULL) +
  theme_minimal(base_size = SLIDE_BASE) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 11),
        plot.title = element_text(face = "bold", size = 14),
        plot.margin = margin(r = 180))

fig1 <- p_seqs / p_reads + plot_layout(heights = c(1, 1.4)) +
  plot_annotation(
    title    = "Figure 1 - The Assignment Problem",
    subtitle = "Edit distance fails for R1: 0 mismatches to both copies. Profile HMMs can resolve this.",
    theme = theme(plot.title = element_text(size = 16, face = "bold"),
                  plot.subtitle = element_text(size = 12, color = "gray40"))
  )

ggsave(file.path(FIGURES_DIR, "fig_hmm1_problem.pdf"), fig1,
       width = 11, height = 7, units = "in")
ggsave(file.path(FIGURES_DIR, "fig_hmm1_problem.png"), fig1,
       width = 11, height = 7, units = "in", dpi = 200)
message("Figure 1 saved.")

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — Per-Copy Profile HMM Construction + Scoring
# ═══════════════════════════════════════════════════════════════════════════════
#   Left:   Emission-probability heatmap for HMM_A and HMM_B
#   Right:  For R2 (assignable) and R1 (ambiguous): compute logL_A, logL_B,
#           ΔlogL, and show as a bar (positive = Copy A, negative = Copy B)

# ── Emission heatmap data ─────────────────────────────────────────────────────
emit_long <- bind_rows(
  expand.grid(copy = "HMM - Copy A", pos = 1:6, base = bases,
              stringsAsFactors = FALSE) |>
    mutate(prob = emit_A[cbind(pos, match(base, bases))]),
  expand.grid(copy = "HMM - Copy B", pos = 1:6, base = bases,
              stringsAsFactors = FALSE) |>
    mutate(prob = emit_B[cbind(pos, match(base, bases))])
)
emit_long$base <- factor(emit_long$base, levels = rev(bases))
emit_long$is_snp <- emit_long$pos %in% snps

p_emit <- ggplot(emit_long, aes(x = pos, y = base, fill = prob)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_tile(data = filter(emit_long, is_snp),
            fill = NA, color = "black", linewidth = 1.8) +
  geom_text(aes(label = ifelse(prob > 0.5,
                               sprintf("%.2f", prob),
                               sprintf("%.3f", prob))),
            size = 3.2, color = "white", fontface = "bold") +
  scale_fill_gradient(low = "#EEF2FF", high = "#1A237E",
                      name = "P(base)") +
  scale_x_continuous(breaks = 1:6, labels = paste0("pos ", 1:6)) +
  facet_wrap(~ copy, ncol = 1) +
  labs(title = "Per-copy emission probabilities  (SNP positions boxed)",
       x = NULL, y = "Base") +
  theme_minimal(base_size = SLIDE_BASE) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 13),
        axis.text.x = element_text(size = 10),
        plot.title = element_text(face = "bold", size = 14))

# ── Log-likelihood scoring ─────────────────────────────────────────────────────
# R1: pos 1-2, bases = A T  (same in both copies → ΔlogL = 0)
# R2: pos 2-4, bases = T C G  (Copy A read, covers SNP at pos 3)
# R3: pos 4-6, bases = G A T  (Copy A read, covers SNP at pos 5)
# R4: pos 4-6, bases = G T T  (Copy B read, covers SNP at pos 5)

score_reads <- data.frame(
  read_id = c("R1 (ambiguous,\npos 1-2: A T)",
              "R2 (Copy A,\npos 2-4: T C G)",
              "R3 (Copy A,\npos 4-6: G A T)",
              "R4 (Copy B,\npos 4-6: G T T)"),
  start   = c(1L, 2L, 4L, 4L),
  bases   = I(list(c("A","T"), c("T","C","G"), c("G","A","T"), c("G","T","T")))
)
score_reads$logL_A <- mapply(log_lik_read,
                             score_reads$bases, list(emit_A), score_reads$start)
score_reads$logL_B <- mapply(log_lik_read,
                             score_reads$bases, list(emit_B), score_reads$start)
score_reads$delta  <- score_reads$logL_A - score_reads$logL_B
score_reads$decision <- ifelse(abs(score_reads$delta) < 0.1, "ambiguous",
                                ifelse(score_reads$delta > 0, "Copy A", "Copy B"))
score_reads$read_id <- factor(score_reads$read_id,
                               levels = rev(score_reads$read_id))

p_score <- ggplot(score_reads, aes(y = read_id, x = delta,
                                    fill = decision)) +
  geom_col(width = 0.55) +
  geom_vline(xintercept = 0, linewidth = 0.8, color = "gray30") +
  geom_text(aes(label = sprintf("D = %.2f", delta),
                x = delta + sign(delta) * 0.15),
            hjust = ifelse(score_reads$delta >= 0, 0, 1),
            size = 4, fontface = "bold") +
  scale_fill_manual(values = c("Copy A" = "#1565C0",
                                "Copy B" = "#B71C1C",
                                "ambiguous" = "#757575"),
                    name = NULL) +
  labs(title = "DlogL = log P(read | HMM_A) - log P(read | HMM_B)",
       x = "D log-likelihood  (+ favours Copy A,  - favours Copy B)",
       y = NULL) +
  theme_minimal(base_size = SLIDE_BASE) +
  theme(legend.position = "top",
        axis.text.y = element_text(size = 11),
        plot.title = element_text(face = "bold", size = 14))

fig2 <- p_emit | p_score +
  plot_annotation(
    title    = "Figure 2 - Per-Copy Profile HMM Construction & Scoring",
    subtitle = paste0(
      "Each HMM emits P(base) ~= 0.95 at the copy's own base, ~= 0.017 elsewhere.\n",
      "At SNP positions (boxed) the two HMMs diverge strongly -> large |DlogL|.\n",
      "R1 covers no SNP -> DlogL ~= 0 -> ambiguous. R2-R4 cover a SNP -> clear assignment."
    ),
    theme = theme(plot.title = element_text(size = 16, face = "bold"),
                  plot.subtitle = element_text(size = 11, color = "gray40"))
  )

ggsave(file.path(FIGURES_DIR, "fig_hmm2_scoring.pdf"), fig2,
       width = 14, height = 7.5, units = "in")
ggsave(file.path(FIGURES_DIR, "fig_hmm2_scoring.png"), fig2,
       width = 14, height = 7.5, units = "in", dpi = 200)
message("Figure 2 saved.")

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — EM Resolves Ambiguous Reads
# ═══════════════════════════════════════════════════════════════════════════════
#   Scenario: 10 reads total
#     6 clearly from Copy A (high +ΔlogL)
#     2 clearly from Copy B (high −ΔlogL)
#     2 ambiguous (ΔlogL ≈ 0, cover no SNP)
#   True: Copy A has 6 + some ambiguous fraction, Copy B has 2 + rest
#   Show EM iterates: fractional weight per ambiguous read converges.

# ── EM implementation ─────────────────────────────────────────────────────────
n_reads <- 10L
n_A_clear <- 6L;  n_B_clear <- 2L;  n_ambig <- 2L

# logL differences for each read category
delta_clear_A  <-  3.5   # log-odds strongly favouring Copy A
delta_clear_B  <- -3.5
delta_ambiguous <-  0.0  # truly ambiguous: equal prob from either copy

# EM loop
# State: theta = P(Copy A)   [probability that a random read came from Copy A]
# E-step: for ambiguous read, weight_A = theta * exp(logL_A) / normaliser
#         For clear reads, weight is effectively 0 or 1 (very high |delta|)
# M-step: theta = (sum of Copy A weights) / n_reads

em_iters <- 15L
theta_vec   <- numeric(em_iters + 1L)
w_ambig_vec <- numeric(em_iters + 1L)  # weight assigned to Copy A for each ambig read
theta_vec[1] <- 0.5   # start uniform

for (iter in seq_len(em_iters)) {
  th <- theta_vec[iter]
  # E-step: weights
  # clear_A reads: P(Copy A | read) ≈ 1   (ΔlogL = 3.5 >> 0)
  # clear_B reads: P(Copy A | read) ≈ 0
  # ambig reads:   P(Copy A | read) = th * exp(0) / (th*1 + (1-th)*1) = th
  #   (because exp(delta_ambiguous) = 1 for ambiguous reads)
  w_A_clear  <- n_A_clear  * 1.0
  w_A_ambig  <- n_ambig    * th       # fractional weight per ambiguous read
  w_B_clear  <- n_B_clear  * 0.0
  w_B_ambig  <- n_ambig    * (1 - th)

  # M-step
  total_A <- w_A_clear + w_A_ambig
  total_B <- w_B_clear + w_B_ambig
  theta_new <- total_A / (total_A + total_B)
  theta_vec[iter + 1L]   <- theta_new
  w_ambig_vec[iter + 1L] <- th   # the weight that was USED in this iteration
}

em_df <- data.frame(
  iter    = 0:em_iters,
  theta   = theta_vec,
  w_ambig = w_ambig_vec
)
em_df$count_A <- em_df$theta * n_reads
em_df$count_B <- (1 - em_df$theta) * n_reads

# ── EM convergence line plot ──────────────────────────────────────────────────
p_conv <- ggplot(em_df, aes(x = iter)) +
  geom_hline(yintercept = (n_A_clear + n_ambig * (n_A_clear / (n_A_clear + n_B_clear))),
             linetype = "dashed", color = "#1565C0", linewidth = 0.8, alpha = 0.5) +
  geom_hline(yintercept = (n_B_clear + n_ambig * (n_B_clear / (n_A_clear + n_B_clear))),
             linetype = "dashed", color = "#B71C1C", linewidth = 0.8, alpha = 0.5) +
  geom_line(aes(y = count_A), color = "#1565C0", linewidth = 1.6) +
  geom_point(aes(y = count_A), color = "#1565C0", size = 3) +
  geom_line(aes(y = count_B), color = "#B71C1C", linewidth = 1.6) +
  geom_point(aes(y = count_B), color = "#B71C1C", size = 3) +
  annotate("label", x = 14, y = em_df$count_A[em_iters + 1] + 0.15,
           label = sprintf("Copy A: %.2f reads", em_df$count_A[em_iters + 1]),
           color = "#1565C0", fill = "white", size = 4.2, fontface = "bold",
           hjust = 1) +
  annotate("label", x = 14, y = em_df$count_B[em_iters + 1] - 0.15,
           label = sprintf("Copy B: %.2f reads", em_df$count_B[em_iters + 1]),
           color = "#B71C1C", fill = "white", size = 4.2, fontface = "bold",
           hjust = 1) +
  scale_x_continuous(breaks = 0:em_iters) +
  scale_y_continuous(breaks = 0:10, limits = c(0, 10)) +
  labs(title = "EM convergence: per-copy read counts",
       x = "EM iteration", y = "Effective read count") +
  theme_classic(base_size = SLIDE_BASE) +
  theme(plot.title = element_text(face = "bold", size = 14))

# ── Stacked bar of read allocations ──────────────────────────────────────────
alloc_df <- bind_rows(
  data.frame(iter = 0L,          # initial state
             type = c("clear -> A", "ambig -> A (weight)", "clear -> B", "ambig -> B (weight)"),
             count = c(n_A_clear, n_ambig * 0.5, n_B_clear, n_ambig * 0.5),
             copy  = c("Copy A","Copy A","Copy B","Copy B")),
  data.frame(iter = em_iters,    # converged state
             type = c("clear -> A", "ambig -> A (weight)", "clear -> B", "ambig -> B (weight)"),
             count = c(n_A_clear,
                       n_ambig * theta_vec[em_iters + 1],
                       n_B_clear,
                       n_ambig * (1 - theta_vec[em_iters + 1])),
             copy  = c("Copy A","Copy A","Copy B","Copy B"))
)
alloc_df$iter_label <- factor(
  ifelse(alloc_df$iter == 0, "Initial (iter 0)", "Converged"),
  levels = c("Initial (iter 0)", "Converged")
)
alloc_df$type <- factor(alloc_df$type,
  levels = c("clear -> A","ambig -> A (weight)","ambig -> B (weight)","clear -> B"))

bar_pal <- c("clear -> A"         = "#1565C0",
             "ambig -> A (weight)"= "#90CAF9",
             "ambig -> B (weight)"= "#EF9A9A",
             "clear -> B"         = "#B71C1C")

p_bars <- ggplot(alloc_df, aes(x = iter_label, y = count, fill = type)) +
  geom_col(width = 0.5, color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.1f", count)),
            position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 4.5) +
  scale_fill_manual(values = bar_pal, name = NULL) +
  labs(title = "Effective read allocations",
       x = NULL, y = "Reads (fractional)") +
  theme_classic(base_size = SLIDE_BASE) +
  theme(legend.position = "right",
        legend.text = element_text(size = 10),
        plot.title = element_text(face = "bold", size = 14))

fig3 <- (p_bars | p_conv) +
  plot_annotation(
    title    = "Figure 3 - EM Assigns Fractional Weights to Ambiguous Reads",
    subtitle = paste0(
      "Setup: 6 clearly-Copy-A reads, 2 clearly-Copy-B reads, 2 ambiguous (cover no SNP).\n",
      "EM E-step: weight of ambiguous read to Copy A = current theta_A (current abundance fraction).\n",
      "M-step: update theta_A = (clear_A + sum(ambig weights_A)) / total.  Converges in ~3 iterations."
    ),
    theme = theme(plot.title = element_text(size = 16, face = "bold"),
                  plot.subtitle = element_text(size = 11, color = "gray40"))
  )

ggsave(file.path(FIGURES_DIR, "fig_hmm3_em.pdf"), fig3,
       width = 14, height = 6.5, units = "in")
ggsave(file.path(FIGURES_DIR, "fig_hmm3_em.png"), fig3,
       width = 14, height = 6.5, units = "in", dpi = 200)
message("Figure 3 saved.")
message("All advisor HMM figures written to figures/fig_hmm{1,2,3}_*")
