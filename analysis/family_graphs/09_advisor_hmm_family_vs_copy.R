suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

FIGURES_DIR <- "figures"
dir.create(FIGURES_DIR, showWarnings = FALSE)

SLIDE_BASE <- 13
copy_cols  <- c(A = "#1565C0", B = "#B71C1C", C = "#2E7D32")
shared_col <- "#e6a817"
exon_h     <- 0.30   # half-height of exon blocks

copy_y <- c(A = 3.0, B = 2.0, C = 1.0)

# ─── Panel A: Gene family + multimapping read ──────────────────────────────────
# Each copy: [shared1 x=0-2] --- [copy-specific x=3-4.2] --- [shared2 x=5.2-7]
# Read spans part of shared1 (x=0.4-1.6): maps to all 3 with equal MAPQ

exon_df <- bind_rows(
  data.frame(copy = c("A","B","C"), x0 = 0,   x1 = 2,
             etype = "shared",  stringsAsFactors = FALSE),
  data.frame(copy = c("A","B","C"), x0 = 3,   x1 = 4.2,
             etype = c("copy_A","copy_B","copy_C"), stringsAsFactors = FALSE),
  data.frame(copy = c("A","B","C"), x0 = 5.2, x1 = 7,
             etype = "shared",  stringsAsFactors = FALSE)
)
exon_df$cy  <- copy_y[exon_df$copy]
exon_df$fill <- ifelse(exon_df$etype == "shared", shared_col,
                       copy_cols[sub("copy_", "", exon_df$etype)])

intron_df <- bind_rows(
  data.frame(copy = c("A","B","C"), x0 = 2.0, x1 = 3.0),
  data.frame(copy = c("A","B","C"), x0 = 4.2, x1 = 5.2)
)
intron_df$cy <- copy_y[intron_df$copy]

read_x0 <- 0.35
read_x1 <- 1.65
read_y0 <- 3.80
read_y1 <- 3.80 + exon_h * 2.4
read_xm  <- (read_x0 + read_x1) / 2

arrow_df <- data.frame(
  copy = c("A","B","C"),
  xend = c(read_xm, read_xm, read_xm),
  yend = copy_y + exon_h
)

p_cartoon <- ggplot() +
  # Introns
  geom_segment(data = intron_df,
               aes(x = x0, xend = x1, y = cy, yend = cy),
               color = "gray40", linewidth = 0.8) +
  # Exon blocks
  geom_rect(data = exon_df,
            aes(xmin = x0, xmax = x1,
                ymin = cy - exon_h, ymax = cy + exon_h,
                fill = fill),
            color = "white", linewidth = 0.5) +
  scale_fill_identity() +
  # Exon labels inside blocks
  geom_text(data = filter(exon_df, etype == "shared", x0 < 1),
            aes(x = (x0 + x1) / 2, y = cy, label = "shared"),
            size = 3.0, color = "white", fontface = "bold") +
  geom_text(data = filter(exon_df, etype == "shared", x0 > 4),
            aes(x = (x0 + x1) / 2, y = cy, label = "shared"),
            size = 3.0, color = "white", fontface = "bold") +
  geom_text(data = filter(exon_df, etype != "shared"),
            aes(x = (x0 + x1) / 2, y = cy,
                label = paste0("copy ", sub("copy_", "", etype), "\nspecific")),
            size = 2.5, color = "white", fontface = "bold") +
  # Copy labels (left margin)
  geom_text(data = data.frame(copy = c("A","B","C"), x = -0.6,
                               y = unname(copy_y),
                               col = unname(copy_cols)),
            aes(x = x, y = y, label = paste0("Copy ", copy), color = col),
            fontface = "bold", size = 4.5, hjust = 1) +
  scale_color_identity() +
  # Multimapping read
  annotate("rect", xmin = read_x0, xmax = read_x1,
           ymin = read_y0, ymax = read_y1,
           fill = "gray25", color = "black", linewidth = 1.3) +
  annotate("text", x = read_xm, y = (read_y0 + read_y1) / 2,
           label = "read", color = "white", fontface = "bold", size = 4.5) +
  # Equal-MAPQ arrows
  geom_segment(data = arrow_df,
               aes(x = xend, xend = xend, y = read_y0, yend = yend),
               linetype = "dashed", color = "gray45", linewidth = 0.7,
               arrow = arrow(length = unit(2.5, "mm"), type = "closed")) +
  # MAPQ label
  annotate("text", x = read_xm + 0.55, y = (read_y0 + copy_y["B"] + exon_h) / 2,
           label = "MAPQ = 0\n(equal alignment\nto all 3)", color = "gray35",
           size = 3.2, fontface = "italic", hjust = 0) +
  coord_cartesian(xlim = c(-1.4, 8.5), ylim = c(0.45, 4.55), clip = "off") +
  labs(title = "A  Multi-mapping read: no diagnostic position, equal MAPQ to all 3 loci",
       x = NULL, y = NULL) +
  theme_minimal(base_size = SLIDE_BASE) +
  theme(panel.grid  = element_blank(),
        axis.text   = element_blank(),
        axis.ticks  = element_blank(),
        plot.title  = element_text(face = "bold", size = 13,
                                   margin = margin(b = 6)))

# ─── Panel B: Per-copy HMM scoring ────────────────────────────────────────────
# Simulate log-likelihoods for a ~120 bp read aligning to a ~95% shared region.
# Copies are NOT 100% identical in the "shared" region (threshold was 85%),
# so subtle positional differences accumulate across the read.
# DlogL relative to best (Copy A):
scores_df <- data.frame(
  copy  = factor(c("A","B","C"), levels = c("A","B","C")),
  delta = c(0.0, -3.5, -1.9),
  fill  = unname(copy_cols[c("A","B","C")])
)

p_scores <- ggplot(scores_df, aes(x = copy, y = delta, fill = fill)) +
  geom_col(width = 0.55, color = "white", linewidth = 0.8) +
  geom_hline(yintercept = 0, linewidth = 0.8, color = "gray30") +
  geom_text(aes(label = ifelse(delta == 0, "best", sprintf("%.1f", delta)),
                y = delta - 0.12),
            vjust = 1, size = 4.2, fontface = "bold", color = "white") +
  scale_fill_identity() +
  scale_x_discrete(labels = c(A = "HMM\nCopy A", B = "HMM\nCopy B",
                               C = "HMM\nCopy C")) +
  scale_y_continuous(limits = c(-5, 0.8),
                     breaks = c(-4, -2, 0),
                     labels = c("-4", "-2", "0")) +
  annotate("label", x = 1, y = 0.55, label = "assigned",
           fill = copy_cols["A"], color = "white",
           size = 3.8, fontface = "bold",
           label.padding = unit(0.25, "lines")) +
  labs(title = "B  Per-copy profile HMM scores",
       subtitle = "DlogL = log P(read | HMM_i) - log P(read | best HMM)\nNo SNP required - cumulative signal over full read",
       x = NULL,
       y = "DlogL  (0 = best)") +
  theme_classic(base_size = SLIDE_BASE) +
  theme(plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9.5, color = "gray40",
                                     lineheight = 1.3),
        axis.text.x   = element_text(face = "bold", size = 11))

# ─── Combine ──────────────────────────────────────────────────────────────────
fig4 <- p_cartoon + p_scores +
  plot_layout(widths = c(2, 1)) +
  plot_annotation(
    title    = "Per-copy profile HMMs resolve multimapping reads without diagnostic SNPs",
    subtitle = "Each copy's HMM scores the full read alignment. Small per-position differences accumulate across ~100 bp to distinguish copies.",
    theme = theme(
      plot.title    = element_text(size = 15, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "gray40")
    )
  )

ggsave(file.path(FIGURES_DIR, "fig_hmm4_family_vs_copy.pdf"), fig4,
       width = 14, height = 6.5, units = "in")
ggsave(file.path(FIGURES_DIR, "fig_hmm4_family_vs_copy.png"), fig4,
       width = 14, height = 6.5, units = "in", dpi = 200)
message("Figure 4 saved.")
