suppressPackageStartupMessages({
  library(ggplot2)
  library(ggraph)
  library(igraph)
  library(dplyr)
  library(patchwork)
})

FIGURES_DIR <- "figures"
DATA_DIR    <- "data"

SLIDE_BASE <- 11
shared_col <- "#e6a817"

# Per-copy colors (7 copies of GOLGA8I)
copy_ids <- c("LOC101150678","LOC129527057","LOC129527171","LOC129527071",
              "LOC134757232","LOC101151050","LOC129527170")
copy_pal <- setNames(
  c("#1565C0","#B71C1C","#2E7D32","#6A1B9A","#EF6C00","#00838F","#5D4037"),
  copy_ids
)

# ═══════════════════════════════════════════════════════════════════════════════
# PANEL A — Real GOLGA8I VG (simplified): shared backbone + the 5-branch bubble
# ═══════════════════════════════════════════════════════════════════════════════
# Structure derived from data/golga_GOLGA8I_vg_seq.rds:
#   - 10 fully-shared nodes (n_copies = 7) form the backbone (S nodes)
#   - One major bubble: shared S5 (n=6) -> 5 copy-specific branches -> shared S7 (n=6)
#   - Plus a tied/secondary bubble (illustrative): S_alt with 2 near-identical branches
#   - Plus a shared region (no bubble): consecutive S nodes on the right

nodes_a <- data.frame(
  name = c("S1","S2","S3","S5",
           "b1","b2","b3","b4","b5",
           "S7","S8","S9",
           "tA","tB",
           "S10","S11","S12"),
  type = c("shared","shared","shared","shared",
           "spec1","spec2","spec3","spec4","spec5",
           "shared","shared","shared",
           "tied_A","tied_B",
           "shared","shared","shared"),
  vx   = c(0, 1, 2, 3,
           4, 4, 4, 4, 4,
           5, 6, 7,
           8, 8,
           9, 10, 11),
  vy   = c(0, 0, 0, 0,
           1.0, 0.5, 0, -0.5, -1.0,
           0, 0, 0,
           0.45, -0.45,
           0, 0, 0),
  stringsAsFactors = FALSE
)
edges_a <- data.frame(
  from = c("S1","S2","S3",
           "S5","S5","S5","S5","S5",
           "b1","b2","b3","b4","b5",
           "S7","S8",
           "S9","S9",
           "tA","tB",
           "S10","S11"),
  to   = c("S2","S3","S5",
           "b1","b2","b3","b4","b5",
           "S7","S7","S7","S7","S7",
           "S8","S9",
           "tA","tB",
           "S10","S10",
           "S11","S12"),
  stringsAsFactors = FALSE
)

g_a   <- graph_from_data_frame(edges_a, directed=TRUE, vertices=nodes_a)
vord  <- V(g_a)$name
ldf   <- nodes_a[match(vord, nodes_a$name), ]
lay_a <- create_layout(g_a, layout="manual", x=ldf$vx, y=ldf$vy)

node_label <- c(S1="S", S2="S", S3="S", S5="S",
                b1="1", b2="2", b3="3", b4="4", b5="5",
                S7="S", S8="S", S9="S",
                tA="A'", tB="B'",
                S10="S", S11="S", S12="S")

node_color <- c(
  shared = shared_col,
  spec1  = copy_pal[1], spec2 = copy_pal[2], spec3 = copy_pal[3],
  spec4  = copy_pal[4], spec5 = copy_pal[5],
  tied_A = "#7B1FA2",   tied_B = "#4A148C"
)

p_vg <- ggraph(lay_a) +
  geom_edge_link(arrow=arrow(length=unit(2,"mm"), type="closed"),
                 end_cap=circle(5.5,"mm"),
                 color="gray50", linewidth=0.7) +
  geom_node_point(aes(color=type), size=11) +
  scale_color_manual(values=node_color, guide="none") +
  geom_node_text(aes(label=node_label[name]),
                 color="white", fontface="bold", size=3.5) +
  # Bubble 1 (real GOLGA8I): 5 spec branches
  annotate("text", x=4.0, y=1.55, label="Bubble 1 (real GOLGA8I)",
           size=3.4, fontface="bold", color="gray15") +
  annotate("text", x=4.0, y=1.38,
           label="5 copy-specific exons -> 5 of 7 copies diverge here",
           size=2.85, color="gray45") +
  # Bubble 2 (hypothetical near-identical)
  annotate("text", x=8.0, y=1.05, label="Bubble 2 (near-identical, illustrative)",
           size=3.2, fontface="bold", color="#4A148C") +
  annotate("text", x=8.0, y=0.88, label="(~99% identity)",
           size=2.75, color="gray45") +
  # Shared region (right side)
  annotate("text", x=10.0, y=0.55, label="Shared region",
           size=3.2, fontface="bold", color=shared_col) +
  annotate("text", x=10.0, y=0.40, label="(no bubble)", size=2.75, color="gray45") +
  # Read R1 in Bubble 1
  annotate("segment", x=4.0, xend=4.0, y=-1.65, yend=-1.32,
           arrow=arrow(length=unit(2.2,"mm")), color=copy_pal[1], linewidth=0.8) +
  annotate("text", x=4.0, y=-1.82, label="R1",
           fontface="bold", size=4.2, color=copy_pal[1]) +
  # Read R2 in Bubble 2
  annotate("segment", x=8.0, xend=8.0, y=-1.65, yend=-0.7,
           arrow=arrow(length=unit(2.2,"mm")), color="#7B1FA2", linewidth=0.8) +
  annotate("text", x=8.0, y=-1.82, label="R2",
           fontface="bold", size=4.2, color="#7B1FA2") +
  # Read R3 on shared backbone
  annotate("segment", x=10.0, xend=10.0, y=-1.65, yend=-0.2,
           arrow=arrow(length=unit(2.2,"mm")), color="gray45", linewidth=0.8) +
  annotate("text", x=10.0, y=-1.82, label="R3",
           fontface="bold", size=4.2, color="gray45") +
  # Family-taxonomy summary box (inset)
  annotate("rect", xmin=11.55, xmax=14.5, ymin=-2.05, ymax=0.95,
           fill="#F8F8F8", color="gray50", linewidth=0.6) +
  annotate("text", x=13.025, y=0.82,
           label="LOCUS TAXONOMY",
           fontface="bold", size=2.95, color="gray15") +
  annotate("text", x=13.025, y=0.64,
           label="(by VG component decomposition)",
           size=2.3, color="gray45", fontface="italic") +
  annotate("segment", x=11.65, xend=14.4, y=0.45, yend=0.45,
           color="gray60", linewidth=0.4) +
  # Type 1
  annotate("text", x=13.025, y=0.29,
           label="K=1: Coherent family",
           size=2.8, color="#2E7D32", fontface="bold") +
  annotate("text", x=13.025, y=0.13,
           label="RBMY 13->10, GOLGA8I 7->7,",
           size=2.4, color="gray45") +
  annotate("text", x=13.025, y=-0.02,
           label="AMY 3->3, GAGE 12->3",
           size=2.4, color="gray45") +
  # Type 2
  annotate("text", x=13.025, y=-0.30,
           label="K>1 (mixed): sub-families",
           size=2.8, color="#EF6C00", fontface="bold") +
  annotate("text", x=13.025, y=-0.46,
           label="MAGEA 5->5 (4 comps),",
           size=2.4, color="gray45") +
  annotate("text", x=13.025, y=-0.61,
           label="APOL 6->6 (3 comps)",
           size=2.4, color="gray45") +
  # Type 3
  annotate("text", x=13.025, y=-0.90,
           label="K=N: fully fragmented",
           size=2.8, color="#B71C1C", fontface="bold") +
  annotate("text", x=13.025, y=-1.06,
           label="DEFB 10->10, PRAMEF 4->4",
           size=2.4, color="gray45") +
  annotate("text", x=13.025, y=-1.21,
           label="(reads never multi-map)",
           size=2.4, color="gray45", fontface="italic") +
  # Theorem note
  annotate("segment", x=11.65, xend=14.4, y=-1.40, yend=-1.40,
           color="gray60", linewidth=0.4) +
  annotate("text", x=13.025, y=-1.58,
           label="Reads on G_k assignable",
           size=2.5, color="gray25") +
  annotate("text", x=13.025, y=-1.74,
           label="only to copies in G_k.",
           size=2.5, color="gray25") +
  annotate("text", x=13.025, y=-1.92,
           label="Quantify per-component.",
           size=2.5, color="gray25", fontface="italic") +
  coord_cartesian(ylim=c(-2.2, 1.8), xlim=c(-0.6, 14.6)) +
  labs(title="A  GOLGA8I (7 copies, 7 classes):  shared backbone + one real 5-branch bubble",
       subtitle="A copy = a path through the VG. When two copies share a path, they form one equivalence class and must be quantified jointly.\nBubble 1 has 5 copy-specific branches = 5 of the 7 GOLGA8I classes. With IsoSeq, a single long read often spans multiple bubbles.") +
  theme_graph(base_size=SLIDE_BASE) +
  theme(plot.title=element_text(face="bold", size=12),
        plot.subtitle=element_text(size=10, color="gray40"))

# ═══════════════════════════════════════════════════════════════════════════════
# PANEL B — Read R1 at Bubble 1 (5 divergent copy-specific branches): HMM works
# ═══════════════════════════════════════════════════════════════════════════════
# Copy 1 (LOC101150678) is the assigned copy because R1 maps to its specific exon
# DlogL values are illustrative but realistic for a ~150 bp IsoSeq read
case1_df <- data.frame(
  copy  = factor(paste0("Cp", 1:5), levels = paste0("Cp", 1:5)),
  delta = c(0, -8.7, -9.1, -8.5, -9.3),
  fill  = unname(copy_pal[1:5])
)

p_case1 <- ggplot(case1_df, aes(x=copy, y=delta, fill=fill)) +
  geom_col(width=0.55, color="white", linewidth=0.6) +
  scale_fill_identity() +
  geom_hline(yintercept=0, linewidth=0.7, color="gray30") +
  geom_text(aes(label=ifelse(delta==0, "best", sprintf("%.1f", delta)),
                y=delta - 0.3),
            vjust=1, color="white", fontface="bold", size=3.2) +
  annotate("label", x=1, y=1.0, label="ASSIGNED",
           fill=copy_pal[1], color="white", size=3.0, fontface="bold",
           label.padding=unit(0.2,"lines")) +
  annotate("text", x=3, y=-11,
           label=paste("R1 lands on copy-specific exon for Cp1.",
                       "Cp1's HMM scores it perfectly;",
                       "the other 4 copies' HMMs see a wrong exon",
                       "in this position -> large negative DlogL.",
                       sep="\n"),
           size=2.75, color="gray30", fontface="italic", hjust=0.5, lineheight=1.25) +
  scale_y_continuous(limits=c(-13, 2)) +
  labs(title="B  R1 at Bubble 1 (5 divergent branches): HMM assigns",
       x=NULL, y="DlogL  (0 = best)") +
  theme_classic(base_size=SLIDE_BASE) +
  theme(plot.title=element_text(face="bold", size=10.5),
        axis.text.x=element_text(size=9, face="bold"))

# ═══════════════════════════════════════════════════════════════════════════════
# PANEL C — Read R2 at hypothetical tied bubble: three resolution paths
# ═══════════════════════════════════════════════════════════════════════════════
strategy_df <- data.frame(
  x       = c(1, 2, 3),
  title_l = c("1. COLLAPSE",
              "2. INTRA-READ LINKAGE",
              "3. EM + PRIORS"),
  body    = c(paste("Members of one",
                    "equivalence class share",
                    "an isoform set T(c) and",
                    "are read-indistinguishable.",
                    "Report joint abundance.",
                    "",
                    "Real: RBMY 13->10 classes",
                    "(class_8 has 3 copies,",
                    " class_7 has 2)", sep="\n"),
              paste("IsoSeq reads are long.",
                    "If the same read also",
                    "spans Bubble 1 (the",
                    "real 5-branch bubble),",
                    "that bubble's HMM",
                    "signal resolves R2 too.",
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
                    "(works as long as ANY",
                    "bubble is informative)", sep="\n")),
  fill    = c("#FFCCBC", "#C8E6C9", "#BBDEFB"),
  border  = c("#D84315", "#2E7D32", "#1565C0")
)

p_case2 <- ggplot(strategy_df) +
  annotate("rect", xmin=0.4, xmax=3.6, ymin=4.4, ymax=5.05,
           fill="#FFF3E0", color="gray35", linewidth=0.6) +
  annotate("text", x=2, y=4.72,
           label="Read R2 at tied bubble:  DlogL ~= 0  -> direct HMM assignment fails",
           size=2.95, fontface="bold", color="#E65100") +
  geom_rect(aes(xmin=x-0.46, xmax=x+0.46, ymin=1.2, ymax=4.2,
                fill=fill, color=border),
            linewidth=0.85) +
  scale_fill_identity() +
  scale_color_identity() +
  geom_rect(aes(xmin=x-0.46, xmax=x+0.46, ymin=3.85, ymax=4.2, fill=border),
            color=NA) +
  geom_text(aes(x=x, y=4.02, label=title_l),
            fontface="bold", size=3.2, color="white") +
  geom_text(aes(x=x, y=2.55, label=body),
            size=2.7, color="gray15", lineheight=1.15) +
  coord_cartesian(xlim=c(0.3, 3.7), ylim=c(1.0, 5.2)) +
  labs(title="C  R2 at tied bubble: three resolution paths") +
  theme_void(base_size=SLIDE_BASE) +
  theme(plot.title=element_text(face="bold", size=10.5, hjust=0,
                                margin=margin(b=4, l=2)))

# ═══════════════════════════════════════════════════════════════════════════════
# PANEL D — Read R3 on shared backbone: EM with priors only
# ═══════════════════════════════════════════════════════════════════════════════
case3_df <- data.frame(
  copy  = factor(paste0("Cp", 1:7), levels = paste0("Cp", 1:7)),
  delta = rep(0, 7),
  fill  = unname(copy_pal[1:7])
)

p_case3 <- ggplot(case3_df, aes(x=copy, y=delta, fill=fill)) +
  geom_col(width=0.55, color="white", linewidth=0.6) +
  scale_fill_identity() +
  geom_hline(yintercept=0, linewidth=0.7, color="gray30") +
  annotate("label", x=4, y=0.85, label="ALL TIED",
           fill="gray55", color="white", size=3.0, fontface="bold",
           label.padding=unit(0.2,"lines")) +
  annotate("text", x=4, y=-6,
           label=paste("R3 lands on a fully-shared backbone node (n=7)",
                       "All 7 HMMs are identical at that position",
                       "-> DlogL = 0 for every copy",
                       "",
                       "Only path: EM with priors",
                       "(theta from reads in Bubble 1 distributes R3)",
                       sep="\n"),
           size=2.75, color="gray25", hjust=0.5, lineheight=1.25, fontface="italic") +
  scale_y_continuous(limits=c(-13, 2)) +
  labs(title="D  R3 on shared backbone: EM is the only option",
       x=NULL, y="DlogL") +
  theme_classic(base_size=SLIDE_BASE) +
  theme(plot.title=element_text(face="bold", size=10.5),
        axis.text.x=element_text(size=9, face="bold"))

# ═══════════════════════════════════════════════════════════════════════════════
# Combine
# ═══════════════════════════════════════════════════════════════════════════════
bottom_row <- p_case1 + p_case2 + p_case3 +
  plot_layout(ncol=3, widths=c(1.0, 1.4, 1.0))

fig7 <- p_vg / bottom_row +
  plot_layout(heights=c(1.0, 1.0)) +
  plot_annotation(
    title    = "Read-to-copy assignment (GOLGA8I, 7 copies):  decision tree",
    subtitle = "Bubble 1 is real (5 copy-specific exons between shared S5 and S7); Bubble 2 is illustrative. A copy = a path through the VG; per-branch HMMs score reads.",
    theme = theme(plot.title=element_text(size=15, face="bold"),
                  plot.subtitle=element_text(size=11, color="gray40",
                                             lineheight=1.4))
  )

ggsave(file.path(FIGURES_DIR, "fig_hmm7_golga8i_real.pdf"), fig7,
       width=16, height=10, units="in")
ggsave(file.path(FIGURES_DIR, "fig_hmm7_golga8i_real.png"), fig7,
       width=16, height=10, units="in", dpi=180)
message("Figure 7 saved.")
