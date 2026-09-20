suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(tidyr)
})

FIGURES_DIR <- "figures"

# ── Toy example: 2 copies of length 6, SNPs at pos 3 (C/G) and pos 5 (A/T) ──
pos   <- 1:6
seq_A <- c("A","T","C","G","A","T")
seq_B <- c("A","T","G","G","T","T")
snps  <- which(seq_A != seq_B)
bases <- c("A","T","C","G")

p_match <- 0.95
p_other <- (1 - p_match) / 3

# Family HMM emission: at each position, fraction of copies with each base
# For 2 copies:
#   conserved position -> 1.0 for the conserved base
#   SNP position       -> 0.5 each for the two alleles
mk_family_emit <- function(seqA, seqB) {
  n <- length(seqA)
  m <- matrix(0, nrow = n, ncol = 4, dimnames = list(seq_along(seqA), bases))
  for (i in seq_along(seqA)) {
    if (seqA[i] == seqB[i]) {
      m[i, ] <- p_other
      m[i, seqA[i]] <- p_match
    } else {
      m[i, seqA[i]] <- 0.5
      m[i, seqB[i]] <- 0.5
    }
  }
  m
}
mk_per_copy_emit <- function(seq) {
  m <- matrix(p_other, nrow = length(seq), ncol = 4,
              dimnames = list(seq_along(seq), bases))
  for (i in seq_along(seq)) m[i, seq[i]] <- p_match
  m
}
emit_fam    <- mk_family_emit(seq_A, seq_B)
emit_copy_A <- mk_per_copy_emit(seq_A)
emit_copy_B <- mk_per_copy_emit(seq_B)

# Long form for plotting
to_long <- function(m, label) {
  df <- as.data.frame(m)
  df$pos <- 1:nrow(df)
  df <- pivot_longer(df, -pos, names_to = "base", values_to = "prob")
  df$base <- factor(df$base, levels = rev(bases))
  df$model <- label
  df$is_snp <- df$pos %in% snps
  df
}
df_fam <- to_long(emit_fam, "Profile HMM built from MSA of A+B (family HMM)")
df_A   <- to_long(emit_copy_A, "Profile HMM built from Copy A only")
df_B   <- to_long(emit_copy_B, "Profile HMM built from Copy B only")

# ── Helper: emission heatmap with letters ────────────────────────────────────
mk_heatmap <- function(df, title_text, badge_text, badge_color) {
  ggplot(df, aes(x = pos, y = base, fill = prob)) +
    geom_tile(color = "white", linewidth = 1.0) +
    geom_tile(data = filter(df, is_snp),
              fill = NA, color = "black", linewidth = 1.8) +
    geom_text(aes(label = ifelse(prob >= 0.05,
                                   sprintf("%.2f", prob), "")),
              size = 3.8, color = "white", fontface = "bold") +
    scale_fill_gradient(low = "#EEF2FF", high = "#1A237E", limits = c(0, 1),
                        name = "P(base)") +
    scale_x_continuous(breaks = 1:6, labels = paste0("pos", 1:6)) +
    annotate("label", x = 6.7, y = 4, label = badge_text,
             fill = badge_color, color = "white", size = 3.2, fontface = "bold",
             hjust = 0, label.padding = unit(0.3, "lines")) +
    labs(title = title_text, x = NULL, y = "Base") +
    coord_cartesian(xlim = c(0.5, 9), clip = "off") +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(size = 10),
          plot.title  = element_text(face = "bold", size = 11),
          legend.position = "bottom",
          legend.key.width = unit(1.3, "cm"))
}

# ── MSA visualisation panel ──────────────────────────────────────────────────
seq_df <- bind_rows(
  data.frame(copy = "Copy A", pos = 1:6, base = seq_A, is_snp = (1:6) %in% snps),
  data.frame(copy = "Copy B", pos = 1:6, base = seq_B, is_snp = (1:6) %in% snps)
)
base_pal <- c(A="#43A047", T="#1E88E5", C="#FB8C00", G="#E53935")

p_msa <- ggplot(seq_df, aes(x = pos, y = copy)) +
  geom_tile(aes(fill = base), color = "white", linewidth = 1.5,
            width = 0.85, height = 0.7) +
  geom_tile(data = filter(seq_df, is_snp), fill = NA, color = "black",
            linewidth = 2.0, width = 0.85, height = 0.7) +
  geom_text(aes(label = base), color = "white", fontface = "bold", size = 6.5) +
  scale_fill_manual(values = base_pal, guide = "none") +
  scale_x_continuous(breaks = 1:6, labels = paste0("pos", 1:6)) +
  labs(title = "Toy example: two copies of length 6, SNPs at pos 3 and pos 5",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        plot.title  = element_text(face = "bold", size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(fontface = "bold", size = 12))

# ── Top section: Classical (family) profile HMM ──────────────────────────────
p_classical <- mk_heatmap(
  df_fam,
  "1.  Classical profile HMM (built from MSA of A and B)",
  "FAMILY MEMBERSHIP\n(your description, correct)\nCANNOT ASSIGN A vs B",
  "#B71C1C")

# ── Bottom section: two per-copy HMMs ────────────────────────────────────────
p_per_A <- mk_heatmap(
  df_A,
  "2a.  Per-copy profile HMM for Copy A",
  "discriminative\nat SNP positions",
  "#1B5E20")
p_per_B <- mk_heatmap(
  df_B,
  "2b.  Per-copy profile HMM for Copy B",
  "discriminative\nat SNP positions",
  "#1B5E20")

# ── Explanation text boxes ───────────────────────────────────────────────────
mk_text_box <- function(text, fill) {
  ggplot() +
    annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = 1,
             fill = fill, color = "gray60", linewidth = 0.4) +
    annotate("text", x = 0.15, y = 0.5, label = text,
             hjust = 0, vjust = 0.5, size = 3.5, color = "gray15",
             lineheight = 1.35) +
    coord_cartesian(xlim = c(0, 10), ylim = c(0, 1)) +
    theme_void()
}

txt_classical <- mk_text_box(paste(
  "Built from the alignment (MSA) of copies A and B (this is what you described in your email).",
  "At SNP positions, the match state has 0.50 / 0.50 emission for the two alleles (averaging across the MSA).",
  "Viterbi gives P(read fits family). There is only ONE model -> no way to compare P(read | A) vs P(read | B).",
  "-> Useful for asking 'does this sequence belong to the family?'  NOT useful for assigning a read to a specific copy.",
  sep = "\n"),
  fill = "#FFEBEE")

txt_per_copy <- mk_text_box(paste(
  "Build ONE profile HMM per copy, from that copy's sequence alone.",
  "At every position, match state emits 0.95 for that copy's base, 0.017 for any other base.",
  "Now there are TWO models. For a read R:  logL_A = log P(R | HMM_A),  logL_B = log P(R | HMM_B).",
  "DlogL = logL_A - logL_B  ->  sign tells us which copy the read came from.",
  "-> Solves the read-to-paralog assignment problem.",
  sep = "\n"),
  fill = "#E8F5E9")

# ── Layout ────────────────────────────────────────────────────────────────────
fig <- (p_msa) /
       (p_classical) /
       (txt_classical) /
       ((p_per_A | p_per_B)) /
       (txt_per_copy) +
  plot_layout(heights = c(0.55, 1.1, 0.6, 1.1, 0.6)) +
  plot_annotation(
    title    = "Classical (family) profile HMM  vs.  per-copy profile HMMs",
    subtitle = "Addressing your first concern. Toy 2-copy, 6-base example. No mention of any other concept.",
    theme = theme(plot.title    = element_text(size = 17, face = "bold"),
                  plot.subtitle = element_text(size = 12, color = "gray40"))
  )

ggsave(file.path(FIGURES_DIR, "fig_advisor_slide_classical_vs_percopy.pdf"), fig,
       width = 13, height = 14, units = "in")
ggsave(file.path(FIGURES_DIR, "fig_advisor_slide_classical_vs_percopy.png"), fig,
       width = 13, height = 14, units = "in", dpi = 180)
message("Slide A saved.")
