suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(png)
})

FIGURES_DIR <- "figures"

# Reads an existing PNG and turns it into a ggplot
load_png <- function(file) {
  img <- readPNG(file.path(FIGURES_DIR, file))
  ggplot() + annotation_custom(rasterGrob(img, width = unit(1, "npc"),
                                            height = unit(1, "npc"))) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

mk_point_box <- function(point_n, advisor_quote, our_response, evidence_files) {
  bg <- ggplot() +
    annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = 10,
             fill = "#FAFAFA", color = "gray50", linewidth = 0.4) +
    annotate("text", x = 0.2, y = 9.4, label = paste0("POINT ", point_n),
             hjust = 0, fontface = "bold", size = 4.5, color = "#1565C0") +
    annotate("rect", xmin = 0.2, xmax = 9.8, ymin = 7.6, ymax = 9.0,
             fill = "#FFFDE7", color = "gray70", linewidth = 0.3) +
    annotate("text", x = 0.3, y = 8.3, label = advisor_quote,
             hjust = 0, vjust = 0.5, size = 3.5, fontface = "italic",
             color = "gray25", lineheight = 1.3) +
    annotate("text", x = 0.2, y = 6.9, label = "OUR RESPONSE",
             hjust = 0, fontface = "bold", size = 4.0, color = "#1B5E20") +
    annotate("text", x = 0.3, y = 4.5, label = our_response,
             hjust = 0, vjust = 1, size = 3.6, color = "gray15",
             lineheight = 1.35) +
    annotate("text", x = 0.2, y = 1.0,
             label = paste("Evidence:", paste(evidence_files, collapse = "  |  ")),
             hjust = 0, fontface = "italic", size = 3.0, color = "gray40") +
    coord_cartesian(xlim = c(0, 10), ylim = c(0, 10)) +
    theme_void()
  bg
}

p_title <- ggplot() +
  annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = 1, fill = "#0D47A1") +
  annotate("text", x = 0.2, y = 0.65,
           label = "Response to your email — point by point",
           hjust = 0, color = "white", fontface = "bold", size = 6) +
  annotate("text", x = 0.2, y = 0.30,
           label = "Each concern + our response + a real figure that backs it.",
           hjust = 0, color = "white", size = 3.6) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 1)) +
  theme_void()

p1 <- mk_point_box(
  1,
  "Profile HMMs capture the alignment of multiple sequences ... so this might be useful for checking if a new copy belongs to a family, but I don't see how this can help you in assigning a read to a paralog.",
  paste(
    "You are correct.  A classical profile HMM built from an MSA of all copies",
    "absorbs SNP variation into degenerate match states (e.g., 0.50 / 0.50 at",
    "every SNP column). One model -> no DlogL across copies -> cannot assign.",
    "",
    "We never meant a classical profile HMM. The figure on the right shows the",
    "contrast: a family HMM is degenerate; per-copy HMMs are not.",
    sep = "\n"),
  c("fig_hmm4_family_vs_copy.png"))

p2 <- mk_point_box(
  2,
  "I also don't understand your idea of building a profile HMM per copy, and how you would use it to solve the assignment problem.",
  paste(
    "A 'per-copy HMM' is a per-position emission model derived from one copy's",
    "transcript sequence:",
    "    P(base = S_c[i] | copy c) = 0.95",
    "    P(other base       | copy c) = 0.017",
    "",
    "For a read R with M matches against copy c (out of L bases):",
    "    logL(R | copy c) = M . log(0.95) + (L - M) . log(0.017)",
    "",
    "This is mathematically equivalent to per-copy Smith-Waterman scoring -",
    "which is exactly what minimap2 already does.",
    sep = "\n"),
  c("fig_hmm2_scoring.png"))

p3 <- mk_point_box(
  3,
  "(Implicit from your email: does this approach actually solve the problem?)",
  paste(
    "On real GOLGA8I FLNC IsoSeq reads:  4 / 4 reads -> MAPQ = 60 with minimap2",
    "alone.  HMM scoring AGREES with minimap2 at the class boundary, never",
    "rescues a tie minimap2 cannot resolve.  The HMM step is therefore",
    "REDUNDANT for our data.",
    "",
    "Where minimap2 ties (MAPQ=0), the HMM also ties; Theorem 1 explains why",
    "(within-class members are read-indistinguishable).",
    "",
    "The fix that DOES close the residual gap is an FLNC-length-consistency",
    "term, not an HMM.  Demo: 100 / 100 simulated reads recovered correctly.",
    sep = "\n"),
  c("fig_flnc_assignment.png  |  fig_class_boundary_resolution.png  |  fig_em_bootstrap.png"))

p4 <- mk_point_box(
  4,
  "Explain these 'in isolation', without mixing them up with other concepts such as variation graphs etc.",
  paste(
    "The HMM explanation above (Points 1-2) contains no VG concepts.",
    "",
    "The VG is a SEPARATE framework that answers a DIFFERENT question:",
    "    'Which copies can in principle be distinguished, and which collapse",
    "     into one equivalence class by sequence identity?'  (resolvability)",
    "",
    "VG does not score reads.  HMM does not define resolvability.  We presented",
    "them as one mixed story; they are mathematically independent and we should",
    "(and now do) present them separately.",
    sep = "\n"),
  c("fig_advisor_floor.png  |  fig_advisor_pipeline.png"))

p_summary <- ggplot() +
  annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = 1, fill = "#E8F5E9",
           color = "#1B5E20", linewidth = 0.5) +
  annotate("text", x = 0.2, y = 0.78,
           label = "BOTTOM LINE",
           hjust = 0, fontface = "bold", size = 4.5, color = "#1B5E20") +
  annotate("text", x = 0.2, y = 0.55,
           label = "You were right to push back. The HMM step as originally framed adds nothing for FLNC reads on coherent multi-copy families with distinct VG paths.",
           hjust = 0, vjust = 0.5, size = 3.5, color = "gray15") +
  annotate("text", x = 0.2, y = 0.30,
           label = "What does load-bearing work: (i) VG -> equivalence classes (the resolvability floor),  (ii) minimap2 (alignment + class-aware MAPQ),",
           hjust = 0, vjust = 0.5, size = 3.5, color = "gray15") +
  annotate("text", x = 0.2, y = 0.10,
           label = "(iii) FLNC-length consistency term (closes residual ties from paralog repeats),  (iv) EM bootstrap (distributes any remaining ambiguous reads).",
           hjust = 0, vjust = 0.5, size = 3.5, color = "gray15") +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 1)) +
  theme_void()

# Stack
fig <- p_title / p1 / p2 / p3 / p4 / p_summary +
  plot_layout(heights = c(0.30, 1, 1, 1.2, 1, 0.55))

ggsave(file.path(FIGURES_DIR, "fig_advisor_point_by_point.pdf"), fig,
       width = 14, height = 17, units = "in")
ggsave(file.path(FIGURES_DIR, "fig_advisor_point_by_point.png"), fig,
       width = 14, height = 17, units = "in", dpi = 150)
message("Point-by-point slide saved.")
