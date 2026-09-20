suppressPackageStartupMessages({
  library(ggplot2)
  library(ggraph)
  library(igraph)
  library(dplyr)
  library(patchwork)
})

FIGURES_DIR <- "figures"
dir.create(FIGURES_DIR, showWarnings = FALSE)

SLIDE_BASE <- 11
copy_A_col <- "#1565C0"
copy_B_col <- "#B71C1C"
copy_C_col <- "#2E7D32"
tied_A_col <- "#7B1FA2"   # purple
tied_B_col <- "#4A148C"   # dark purple (similar to tied_A)
shared_col <- "#e6a817"

# ═══════════════════════════════════════════════════════════════════════════════
# PANEL A — VG topology with 2 bubble types + a shared region
# ═══════════════════════════════════════════════════════════════════════════════
# s1 -> {bA, bB} -> s2 -> {tA, tB} -> s3 -> s4
#   Bubble 1: divergent branches (bA blue / bB red)
#   Bubble 2: near-identical branches (tA / tB both purple)
#   s3 -> s4: shared region, no bubble

nodes_a <- data.frame(
  name = c("s1","bA","bB","s2","tA","tB","s3","s4"),
  type = c("shared","div_A","div_B","shared","tied_A","tied_B","shared","shared"),
  vx   = c(1, 2, 2, 3, 4, 4, 5, 6),
  vy   = c(0, 0.55,-0.55, 0, 0.42,-0.42, 0, 0),
  stringsAsFactors = FALSE
)
edges_a <- data.frame(
  from = c("s1","s1","bA","bB","s2","s2","tA","tB","s3"),
  to   = c("bA","bB","s2","s2","tA","tB","s3","s3","s4"),
  stringsAsFactors = FALSE
)
g_a   <- graph_from_data_frame(edges_a, directed=TRUE, vertices=nodes_a)
vord  <- V(g_a)$name
ldf_a <- nodes_a[match(vord, nodes_a$name), ]
lay_a <- create_layout(g_a, layout="manual", x=ldf_a$vx, y=ldf_a$vy)

node_pal <- c(shared = shared_col,
              div_A  = copy_A_col,  div_B  = copy_B_col,
              tied_A = tied_A_col,  tied_B = tied_B_col)

p_vg <- ggraph(lay_a) +
  geom_edge_link(arrow=arrow(length=unit(2.5,"mm"), type="closed"),
                 end_cap=circle(6.0,"mm"),
                 color="gray45", linewidth=0.85) +
  geom_node_point(aes(color=type), size=15) +
  scale_color_manual(values=node_pal, guide="none") +
  geom_node_text(aes(label=c("S","A","B","S","A'","B'","S","S")),
                 color="white", fontface="bold", size=4.2) +
  # Bubble 1 label
  annotate("text", x=2.0, y=1.18, label="Bubble 1: divergent branches",
           size=3.3, fontface="bold", color=copy_A_col) +
  annotate("text", x=2.0, y=1.02, label="(~75% identity)",
           size=2.8, color="gray45") +
  # Bubble 2 label
  annotate("text", x=4.0, y=1.18, label="Bubble 2: near-identical branches",
           size=3.3, fontface="bold", color=tied_A_col) +
  annotate("text", x=4.0, y=1.02, label="(~99% identity)",
           size=2.8, color="gray45") +
  # Shared region label
  annotate("text", x=5.5, y=0.55, label="Shared region", size=3.3,
           fontface="bold", color=shared_col) +
  annotate("text", x=5.5, y=0.40, label="(no bubble)", size=2.8, color="gray45") +
  # Read markers below
  annotate("segment", x=2.0, xend=2.0, y=-1.05, yend=-0.78,
           arrow=arrow(length=unit(2.2,"mm")), color=copy_A_col, linewidth=0.8) +
  annotate("text", x=2.0, y=-1.20, label="R1",
           fontface="bold", size=4.2, color=copy_A_col) +
  annotate("segment", x=4.0, xend=4.0, y=-1.05, yend=-0.65,
           arrow=arrow(length=unit(2.2,"mm")), color=tied_A_col, linewidth=0.8) +
  annotate("text", x=4.0, y=-1.20, label="R2",
           fontface="bold", size=4.2, color=tied_A_col) +
  annotate("segment", x=5.5, xend=5.5, y=-1.05, yend=-0.18,
           arrow=arrow(length=unit(2.2,"mm")), color="gray45", linewidth=0.8) +
  annotate("text", x=5.5, y=-1.20, label="R3",
           fontface="bold", size=4.2, color="gray45") +
  coord_cartesian(ylim=c(-1.4, 1.3), xlim=c(0.5, 6.5)) +
  labs(title="A  VG topology: a copy is a branch of a bubble",
       subtitle="Per-branch HMMs are built from each branch's sequence. With IsoSeq, a single long read often spans multiple bubbles.") +
  theme_graph(base_size=SLIDE_BASE) +
  theme(plot.title=element_text(face="bold", size=12),
        plot.subtitle=element_text(size=10, color="gray40"))

# ═══════════════════════════════════════════════════════════════════════════════
# PANEL B — Read R1 at divergent bubble: HMM works
# ═══════════════════════════════════════════════════════════════════════════════
case1_df <- data.frame(
  copy  = factor(c("Copy A","Copy B"), levels=c("Copy A","Copy B")),
  delta = c(0, -4.2),
  fill  = c(copy_A_col, copy_B_col)
)

p_case1 <- ggplot(case1_df, aes(x=copy, y=delta, fill=fill)) +
  geom_col(width=0.5, color="white", linewidth=0.7) +
  scale_fill_identity() +
  geom_hline(yintercept=0, linewidth=0.7, color="gray30") +
  geom_text(aes(label=ifelse(delta==0, "best", sprintf("%.1f", delta)),
                y=delta - 0.2),
            vjust=1, color="white", fontface="bold", size=3.8) +
  annotate("label", x=1, y=0.7, label="ASSIGNED",
           fill=copy_A_col, color="white", size=3.2, fontface="bold",
           label.padding=unit(0.22,"lines")) +
  annotate("text", x=1.5, y=-5.3,
           label="HMM_A and HMM_B differ at many positions\nDlogL is large -> direct, hard assignment",
           size=2.8, color="gray30", fontface="italic", hjust=0.5,
           lineheight=1.2) +
  scale_y_continuous(limits=c(-6, 1.3)) +
  labs(title="B  R1 at divergent bubble: HMM assigns",
       x=NULL, y="DlogL  (0 = best)") +
  theme_classic(base_size=SLIDE_BASE) +
  theme(plot.title=element_text(face="bold", size=11),
        axis.text.x=element_text(size=10, face="bold"))

# ═══════════════════════════════════════════════════════════════════════════════
# PANEL C — Read R2 at tied bubble: three resolution paths
# ═══════════════════════════════════════════════════════════════════════════════
strategy_df <- data.frame(
  x       = c(1, 2, 3),
  title_l = c("1. COLLAPSE",
              "2. INTRA-READ LINKAGE",
              "3. EM + PRIORS"),
  body    = c(paste("Merge the near-identical",
                    "branches into one",
                    "expression unit.",
                    "Report combined",
                    "abundance.",
                    "",
                    "(honest answer when",
                    "data cannot separate)", sep="\n"),
              paste("IsoSeq reads are long.",
                    "If the same read also",
                    "spans a divergent bubble",
                    "in this locus, that",
                    "bubble's HMM signal",
                    "resolves R2 too.",
                    "",
                    "(per-read linkage across",
                    "bubbles)", sep="\n"),
              paste("Estimate theta from",
                    "reads in informative",
                    "bubbles elsewhere",
                    "in the locus.",
                    "Distribute R2 by",
                    "current theta.",
                    "",
                    "(works as long as",
                    "ANY bubble is informative)", sep="\n")),
  fill    = c("#FFCCBC", "#C8E6C9", "#BBDEFB"),
  border  = c("#D84315", "#2E7D32", "#1565C0")
)

p_case2 <- ggplot(strategy_df) +
  # Header banner
  annotate("rect", xmin=0.4, xmax=3.6, ymin=4.4, ymax=5.05,
           fill="#FFF3E0", color="gray35", linewidth=0.7) +
  annotate("text", x=2, y=4.72,
           label="Read R2 in tied bubble:  DlogL ~= 0  -> direct HMM assignment fails",
           size=3.0, fontface="bold", color="#E65100") +
  # Strategy boxes
  geom_rect(aes(xmin=x-0.46, xmax=x+0.46, ymin=1.2, ymax=4.2,
                fill=fill, color=border),
            linewidth=0.9) +
  scale_fill_identity() +
  scale_color_identity() +
  geom_rect(aes(xmin=x-0.46, xmax=x+0.46, ymin=3.85, ymax=4.2, fill=border),
            color=NA) +
  geom_text(aes(x=x, y=4.02, label=title_l),
            fontface="bold", size=3.4, color="white") +
  geom_text(aes(x=x, y=2.55, label=body),
            size=2.85, color="gray15", lineheight=1.15) +
  coord_cartesian(xlim=c(0.3, 3.7), ylim=c(1.0, 5.2)) +
  labs(title="C  R2 at tied bubble: three resolution paths",
       x=NULL, y=NULL) +
  theme_void(base_size=SLIDE_BASE) +
  theme(plot.title=element_text(face="bold", size=11, hjust=0,
                                margin=margin(b=4, l=2)))

# ═══════════════════════════════════════════════════════════════════════════════
# PANEL D — Read R3 in purely shared region: EM with priors only
# ═══════════════════════════════════════════════════════════════════════════════
case3_df <- data.frame(
  copy  = factor(c("Copy A","Copy B","Copy C"), levels=c("Copy A","Copy B","Copy C")),
  delta = c(0, 0, 0),
  fill  = c(copy_A_col, copy_B_col, copy_C_col)
)

p_case3 <- ggplot(case3_df, aes(x=copy, y=delta, fill=fill)) +
  geom_col(width=0.5, color="white", linewidth=0.7) +
  scale_fill_identity() +
  geom_hline(yintercept=0, linewidth=0.7, color="gray30") +
  annotate("label", x=2, y=0.65, label="ALL TIED",
           fill="gray55", color="white", size=3.2, fontface="bold",
           label.padding=unit(0.22,"lines")) +
  annotate("text", x=2, y=-2.6,
           label=paste("HMMs are identical in shared regions",
                       "by construction -> DlogL = 0 for all",
                       "",
                       "Only path: EM with priors",
                       "(set theta from informative bubbles",
                       "elsewhere; distribute R3 by theta)",
                       sep="\n"),
           size=2.8, color="gray25", hjust=0.5, lineheight=1.25, fontface="italic") +
  scale_y_continuous(limits=c(-6, 1.3)) +
  labs(title="D  R3 in shared region: EM is the only option",
       x=NULL, y="DlogL") +
  theme_classic(base_size=SLIDE_BASE) +
  theme(plot.title=element_text(face="bold", size=11),
        axis.text.x=element_text(size=10, face="bold"))

# ═══════════════════════════════════════════════════════════════════════════════
# Combine: A spans top, B/C/D in bottom row
# ═══════════════════════════════════════════════════════════════════════════════
bottom_row <- p_case1 + p_case2 + p_case3 +
  plot_layout(ncol=3, widths=c(1.0, 1.4, 1.0))

fig6 <- p_vg / bottom_row +
  plot_layout(heights=c(1.0, 1.0)) +
  plot_annotation(
    title    = "Read-to-copy assignment: a decision tree",
    subtitle = "A copy = a branch of a VG bubble. HMM works at divergent bubbles; tied bubbles have three resolution paths; shared regions need EM with priors.",
    theme = theme(plot.title=element_text(size=15, face="bold"),
                  plot.subtitle=element_text(size=11, color="gray40",
                                             lineheight=1.4))
  )

ggsave(file.path(FIGURES_DIR, "fig_hmm6_decision_tree.pdf"), fig6,
       width=16, height=10, units="in")
ggsave(file.path(FIGURES_DIR, "fig_hmm6_decision_tree.png"), fig6,
       width=16, height=10, units="in", dpi=180)
message("Figure 6 saved.")
