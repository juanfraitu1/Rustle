suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(tidyr)
})

FIGURES_DIR <- "figures"

# Same toy example as Slide A
pos   <- 1:6
seq_A <- c("A","T","C","G","A","T")
seq_B <- c("A","T","G","G","T","T")
snps  <- which(seq_A != seq_B)
bases <- c("A","T","C","G")
p_match <- 0.95; p_other <- (1 - p_match) / 3
LOG_M <- log(p_match); LOG_X <- log(p_other)

emit_A <- matrix(p_other, nrow = 6, ncol = 4,
                  dimnames = list(seq_along(seq_A), bases))
emit_B <- matrix(p_other, nrow = 6, ncol = 4,
                  dimnames = list(seq_along(seq_B), bases))
for (i in 1:6) emit_A[i, seq_A[i]] <- p_match
for (i in 1:6) emit_B[i, seq_B[i]] <- p_match

# ── Define 3 toy reads ────────────────────────────────────────────────────────
reads <- list(
  R1 = list(label = "R1: read covers pos 2-4 (includes the SNP at pos 3)",
             bases = c("T","C","G"), start = 2L,
             source = "Copy A (read carries A's allele 'C' at pos 3)"),
  R2 = list(label = "R2: read covers pos 4-6 (includes the SNP at pos 5)",
             bases = c("G","T","T"), start = 4L,
             source = "Copy B (read carries B's allele 'T' at pos 5)"),
  R3 = list(label = "R3: read covers pos 1-2 (no SNP coverage)",
             bases = c("A","T"), start = 1L,
             source = "ambiguous (covers no diagnostic position)")
)

# ── Score one read against one copy's HMM ─────────────────────────────────────
score <- function(read_bases, start, emit) {
  rows <- seq_along(read_bases) + start - 1L
  cols <- match(read_bases, bases)
  per_pos <- emit[cbind(rows, cols)]
  list(per_pos_probs = per_pos,
       logL         = sum(log(per_pos)),
       L            = length(read_bases))
}

# Compute and stash
score_records <- list()
for (rn in names(reads)) {
  r <- reads[[rn]]
  sA <- score(r$bases, r$start, emit_A)
  sB <- score(r$bases, r$start, emit_B)
  score_records[[rn]] <- list(
    read = rn,
    seq  = paste(r$bases, collapse = ""),
    start = r$start,
    end   = r$start + length(r$bases) - 1L,
    source = r$source,
    label = r$label,
    pA = sA$per_pos_probs, pB = sB$per_pos_probs,
    logLA = sA$logL, logLB = sB$logL,
    delta = sA$logL - sB$logL
  )
}

# ── Build per-read "calculation card" panels ──────────────────────────────────
mk_card <- function(rec) {
  L <- length(rec$pA)
  pos_idx <- rec$start:(rec$start + L - 1L)
  txt_lines <- c()
  txt_lines <- c(txt_lines,
                  sprintf("read sequence:   %s   (covers positions %d-%d)",
                          rec$seq, rec$start, rec$end))
  txt_lines <- c(txt_lines, "")
  txt_lines <- c(txt_lines,
                  sprintf("vs HMM_A:  P_per_pos = %s",
                          paste(sprintf("%.3f", rec$pA), collapse = " . ")))
  txt_lines <- c(txt_lines,
                  sprintf("           logL_A   = %s", sprintf("%.3f", rec$logLA)))
  txt_lines <- c(txt_lines, "")
  txt_lines <- c(txt_lines,
                  sprintf("vs HMM_B:  P_per_pos = %s",
                          paste(sprintf("%.3f", rec$pB), collapse = " . ")))
  txt_lines <- c(txt_lines,
                  sprintf("           logL_B   = %s", sprintf("%.3f", rec$logLB)))
  txt_lines <- c(txt_lines, "")
  call <- if (abs(rec$delta) < 0.05) "AMBIGUOUS (DlogL ~= 0)"
          else if (rec$delta > 0) "-> assign to Copy A"
          else                    "-> assign to Copy B"
  txt_lines <- c(txt_lines,
                  sprintf("DlogL = logL_A - logL_B = %+.3f       %s",
                          rec$delta, call))
  body <- paste(txt_lines, collapse = "\n")

  ggplot() +
    annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = 10,
             fill = "#FAFAFA", color = "gray50", linewidth = 0.4) +
    annotate("text", x = 0.2, y = 9.4,
             label = rec$label,
             hjust = 0, fontface = "bold", size = 3.5, color = "#0D47A1") +
    annotate("text", x = 0.2, y = 8.6,
             label = paste0("(simulated source: ", rec$source, ")"),
             hjust = 0, size = 3.0, color = "gray40", fontface = "italic") +
    annotate("text", x = 0.3, y = 4.5, label = body, family = "mono",
             hjust = 0, vjust = 0.5, size = 3.1, color = "gray10",
             lineheight = 1.25) +
    coord_cartesian(xlim = c(0, 10), ylim = c(0, 10)) +
    theme_void()
}
cards <- lapply(score_records, mk_card)

# ── Header explanation ────────────────────────────────────────────────────────
p_header <- ggplot() +
  annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = 10,
           fill = "#E3F2FD", color = "#0D47A1", linewidth = 0.5) +
  annotate("text", x = 0.2, y = 9.0,
           label = "How to USE per-copy profile HMMs to assign a read to a copy",
           hjust = 0, fontface = "bold", size = 5.0, color = "#0D47A1") +
  annotate("text", x = 0.3, y = 7.2,
           label = "Two HMMs:  HMM_A (built from Copy A's sequence),  HMM_B (built from Copy B's sequence).",
           hjust = 0, size = 3.6, color = "gray15") +
  annotate("text", x = 0.3, y = 5.7,
           label = "For each match state i:   P(base = S_c[i] | HMM_c) = 0.95   and   P(other base) = 0.017.",
           hjust = 0, size = 3.6, color = "gray15", family = "mono") +
  annotate("text", x = 0.3, y = 4.2,
           label = "Given a read R with bases r_1 ... r_L starting at position k:",
           hjust = 0, size = 3.6, color = "gray15") +
  annotate("text", x = 0.5, y = 2.7,
           label = "logL(R | HMM_c)  =  sum_i log P(r_i | HMM_c at position k+i-1)",
           hjust = 0, size = 4.0, color = "#0D47A1", fontface = "bold",
           family = "mono") +
  annotate("text", x = 0.3, y = 1.0,
           label = "DlogL = logL_A - logL_B  ->  sign of DlogL tells us which copy R came from.",
           hjust = 0, size = 3.6, color = "gray15") +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 10)) +
  theme_void()

# ── Tiny visual of the two sequences and read positions ──────────────────────
base_pal <- c(A="#43A047", T="#1E88E5", C="#FB8C00", G="#E53935")
seq_df <- bind_rows(
  data.frame(copy = "Copy A", pos = 1:6, base = seq_A, is_snp = (1:6) %in% snps),
  data.frame(copy = "Copy B", pos = 1:6, base = seq_B, is_snp = (1:6) %in% snps)
)
read_df <- bind_rows(
  data.frame(read = "R1", y = 3.5, start = 2, end = 4, label_pos = 3),
  data.frame(read = "R2", y = 4.2, start = 4, end = 6, label_pos = 5),
  data.frame(read = "R3", y = 4.9, start = 1, end = 2, label_pos = 1.5)
)

p_visual <- ggplot() +
  geom_tile(data = seq_df, aes(x = pos, y = copy, fill = base),
            color = "white", linewidth = 1.5, width = 0.85, height = 0.65) +
  geom_tile(data = filter(seq_df, is_snp),
            aes(x = pos, y = copy), fill = NA, color = "black",
            linewidth = 1.8, width = 0.85, height = 0.65) +
  geom_text(data = seq_df, aes(x = pos, y = copy, label = base),
            color = "white", fontface = "bold", size = 5.5) +
  scale_fill_manual(values = base_pal, guide = "none") +
  # reads as bars above
  geom_rect(data = read_df,
            aes(xmin = start - 0.4, xmax = end + 0.4,
                ymin = y - 0.15, ymax = y + 0.15),
            fill = "gray20", color = "black", linewidth = 0.5) +
  geom_text(data = read_df, aes(x = label_pos, y = y, label = read),
            color = "white", fontface = "bold", size = 3.5) +
  scale_x_continuous(breaks = 1:6, labels = paste0("pos", 1:6)) +
  scale_y_discrete(limits = c("Copy A", "Copy B")) +
  coord_cartesian(ylim = c(0.5, 5.4)) +
  labs(title = "Toy example: copies A, B and 3 reads (R1, R2, R3)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(fontface = "bold", size = 12))

# ── DlogL bar chart ──────────────────────────────────────────────────────────
bar_df <- data.frame(
  read = names(reads),
  delta = sapply(score_records, function(r) r$delta),
  call  = sapply(score_records, function(r) {
    if (abs(r$delta) < 0.05) "ambiguous"
    else if (r$delta > 0) "Copy A"
    else "Copy B"
  })
)
bar_df$read <- factor(bar_df$read, levels = c("R1","R2","R3"))

p_bars <- ggplot(bar_df, aes(x = read, y = delta,
                              fill = call)) +
  geom_col(width = 0.5, color = "white", linewidth = 0.5) +
  geom_hline(yintercept = 0, color = "gray30") +
  geom_text(aes(label = sprintf("DlogL = %+.2f", delta),
                y = delta + sign(delta) * 0.4),
            size = 3.6, fontface = "bold", color = "gray20") +
  scale_fill_manual(values = c("Copy A" = "#1565C0",
                                "Copy B" = "#B71C1C",
                                "ambiguous" = "#757575"),
                    name = NULL) +
  labs(title = "Resulting assignment per read",
       x = NULL, y = "DlogL  (sign of DlogL -> copy)") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(face = "bold", size = 12),
        legend.position = "top")

# ── Combine ───────────────────────────────────────────────────────────────────
top <- (p_header) | (p_visual)
mid <- cards$R1 | cards$R2 | cards$R3
fig <- top / mid / p_bars +
  plot_layout(heights = c(1.0, 1.5, 0.8)) +
  plot_annotation(
    title    = "Per-copy profile HMMs: worked toy example",
    subtitle = "Addressing your second question: 'how do you USE per-copy profile HMMs to solve the assignment problem?'",
    theme = theme(plot.title    = element_text(size = 17, face = "bold"),
                  plot.subtitle = element_text(size = 12, color = "gray40"))
  )

ggsave(file.path(FIGURES_DIR, "fig_advisor_slide_toy_assignment.pdf"), fig,
       width = 18, height = 12, units = "in")
ggsave(file.path(FIGURES_DIR, "fig_advisor_slide_toy_assignment.png"), fig,
       width = 18, height = 12, units = "in", dpi = 180)
message("Slide B saved.")
