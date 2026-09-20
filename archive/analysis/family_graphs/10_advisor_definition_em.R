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
copy_cols  <- c(A = "#1565C0", B = "#B71C1C", C = "#2E7D32")
shared_col <- "#e6a817"
eh <- 0.24  # exon half-height

# ═══════════════════════════════════════════════════════════════════════════════
# PANEL A — Multi-copy gene family: copies vs isoforms  (RBMY, 3 of 13 copies)
# ═══════════════════════════════════════════════════════════════════════════════
cy <- c(A = 3.0, B = 2.0, C = 1.0)

# Three exons per copy: [shared | copy-specific | shared]
ex_df <- do.call(rbind, lapply(c("A","B","C"), function(cp) {
  data.frame(x0   = c(0.0, 2.0, 3.2),
             x1   = c(1.5, 3.0, 4.7),
             type = c("shared", paste0("copy_",cp), "shared"),
             copy = cp, cy = cy[cp], stringsAsFactors = FALSE)
}))
ex_df$fill <- ifelse(ex_df$type == "shared", shared_col,
                      copy_cols[sub("copy_","", ex_df$type)])

it_df <- do.call(rbind, lapply(c("A","B","C"), function(cp) {
  data.frame(x0 = c(1.5, 3.0), x1 = c(2.0, 3.2), cy = cy[cp])
}))

# Isoforms of Copy A shown below main rows
iso_y1 <- 0.40; iso_y2 <- -0.22

p_def <- ggplot() +
  # Introns (main copies)
  geom_segment(data = it_df, aes(x=x0, xend=x1, y=cy, yend=cy),
               color="gray40", linewidth=0.7) +
  # Exon blocks (main copies)
  geom_rect(data = ex_df,
            aes(xmin=x0, xmax=x1, ymin=cy-eh, ymax=cy+eh, fill=fill),
            color="white", linewidth=0.4) +
  scale_fill_identity() +
  # Copy labels: bold letter + italic locus ID
  geom_text(data = data.frame(x=-0.42, y=unname(cy)+0.37,
                               lab=paste0("Copy ",c("A","B","C")),
                               col=unname(copy_cols)),
            aes(x=x, y=y, label=lab, color=col),
            fontface="bold", size=3.7, hjust=1) +
  geom_text(data = data.frame(x=-0.42, y=unname(cy),
                               lab=c("LOC389895","LOC729966","LOC101928837"),
                               col=unname(copy_cols)),
            aes(x=x, y=y, label=lab, color=col),
            fontface="italic", size=2.6, hjust=1) +
  scale_color_identity() +
  # Exon text labels (top copy only — same positions on all copies)
  annotate("text", x=0.75, y=cy["A"], label="shared\nexon",
           size=2.2, color="white", fontface="bold") +
  annotate("text", x=3.95, y=cy["A"], label="shared\nexon",
           size=2.2, color="white", fontface="bold") +
  annotate("text", x=2.5, y=cy["A"], label="A-specific", size=2.1, color="white", fontface="bold") +
  annotate("text", x=2.5, y=cy["B"], label="B-specific", size=2.1, color="white", fontface="bold") +
  annotate("text", x=2.5, y=cy["C"], label="C-specific", size=2.1, color="white", fontface="bold") +
  # Shared-across-copies bracket (right side)
  annotate("segment", x=5.02, xend=5.02, y=cy["C"]-eh-0.06, yend=cy["A"]+eh+0.06,
           color=shared_col, linewidth=1.8,
           arrow=arrow(ends="both", length=unit(1.5,"mm"))) +
  annotate("text", x=5.28, y=2.0,
           label="shared\nexons\n(>=85%\nidentity)",
           size=2.4, color=shared_col, fontface="bold", hjust=0) +
  # Family header
  annotate("label", x=2.35, y=3.62,
           label="RBMY family  (13 copies on chrY — 3 shown)",
           fill="gray96", color="gray20", fontface="bold", size=3.0,
           label.padding=unit(0.2,"lines")) +
  # Dotted bracket from isoforms up to Copy A
  annotate("segment", x=-0.06, xend=-0.06, y=iso_y2-eh*0.7, yend=cy["A"]-eh-0.04,
           color=copy_cols["A"], linetype="dotted", linewidth=0.8) +
  annotate("text", x=-0.42, y=(iso_y1+iso_y2)/2,
           label="isoforms of Copy A\n(same locus,\ndiff. transcripts)",
           size=2.5, color=copy_cols["A"], fontface="italic", hjust=1) +
  # Isoform 1: all 3 exons
  geom_segment(data=data.frame(x0=c(1.5,3.0), x1=c(2.0,3.2)),
               aes(x=x0, xend=x1, y=iso_y1, yend=iso_y1),
               color="gray45", linewidth=0.5) +
  geom_rect(data=data.frame(x0=c(0,2.0,3.2), x1=c(1.5,3.0,4.7),
                              fill=c(shared_col, copy_cols["A"], shared_col)),
            aes(xmin=x0, xmax=x1, ymin=iso_y1-eh*.65, ymax=iso_y1+eh*.65, fill=fill),
            color="white", linewidth=0.3) +
  annotate("text", x=4.88, y=iso_y1, label="isoform 1 (all exons)",
           size=2.5, color="gray35", hjust=0) +
  # Isoform 2: skips copy-specific exon (dashed intron across the gap)
  annotate("segment", x=1.5, xend=3.2, y=iso_y2, yend=iso_y2,
           color="gray45", linewidth=0.5, linetype="dashed") +
  geom_rect(data=data.frame(x0=c(0,3.2), x1=c(1.5,4.7),
                              fill=c(shared_col, shared_col)),
            aes(xmin=x0, xmax=x1, ymin=iso_y2-eh*.65, ymax=iso_y2+eh*.65, fill=fill),
            color="white", linewidth=0.3) +
  annotate("text", x=4.88, y=iso_y2, label="isoform 2 (copy-specific exon skipped)",
           size=2.5, color="gray35", hjust=0) +
  coord_cartesian(xlim=c(-2.7, 8.0), ylim=c(-0.65, 4.05), clip="off") +
  labs(title="A  Multi-copy gene family: copies (different loci) vs isoforms (same locus)",
       x=NULL, y=NULL) +
  theme_minimal(base_size=SLIDE_BASE) +
  theme(panel.grid=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),
        plot.title=element_text(face="bold", size=10.5))

# ═══════════════════════════════════════════════════════════════════════════════
# PANEL B — VG bubble types: copy-divergence vs isoform
# ═══════════════════════════════════════════════════════════════════════════════
# Nodes: s1(shared) → {cA, cB} → s2(shared) → {e_iso, skip} → s3(shared)
#   copy-divergence bubble: s1 → cA or cB → s2
#   isoform bubble:          s2 → e_iso or skip → s3

nodes_b <- data.frame(
  name = c("s1","cA","cB","s2","e_iso","s3"),
  type = c("shared","copy_A","copy_B","shared","isoform","shared"),
  vx   = c(1, 2, 2, 3, 4, 5),
  vy   = c(0, 0.55,-0.55, 0, 0.50, 0),
  stringsAsFactors = FALSE
)
edges_b <- data.frame(
  from  = c("s1","s1","cA","cB","s2","s2","e_iso"),
  to    = c("cA","cB","s2","s2","e_iso","s3","s3"),
  etype = c("copy_A","copy_B","copy_A","copy_B","isoform","skip","isoform"),
  stringsAsFactors = FALSE
)

g_b   <- graph_from_data_frame(edges_b, directed=TRUE, vertices=nodes_b)
vord  <- V(g_b)$name
ldf_b <- nodes_b[match(vord, nodes_b$name), ]
lay_b <- create_layout(g_b, layout="manual", x=ldf_b$vx, y=ldf_b$vy)

p_vg <- ggraph(lay_b) +
  geom_edge_link(aes(color=etype, linetype=etype),
                 arrow    = arrow(length=unit(2.5,"mm"), type="closed"),
                 end_cap  = circle(6.5,"mm"),
                 linewidth = 0.95) +
  scale_edge_color_manual(
    values  = c(copy_A=copy_cols["A"], copy_B=copy_cols["B"],
                isoform="gray35", skip="gray35"),
    guide   = "none") +
  scale_edge_linetype_manual(
    values  = c(copy_A="solid", copy_B="solid", isoform="solid", skip="dashed"),
    guide   = "none") +
  geom_node_point(aes(color=type), size=16) +
  scale_color_manual(
    values  = c(shared=shared_col, copy_A=copy_cols["A"],
                copy_B=copy_cols["B"], isoform="#7B1FA2"),
    guide   = "none") +
  geom_node_text(aes(label=c("S","A","B","S","iso","S")),
                 color="white", fontface="bold", size=4.5) +
  # Bubble type annotations
  annotate("text", x=2.0, y=1.12,
           label="copy-divergence bubble", size=3.2, fontface="bold.italic", color="gray15") +
  annotate("text", x=2.0, y=0.96,
           label="Copy A and B have\ndifferent exon sequences", size=2.7, color="gray45") +
  annotate("text", x=4.0, y=1.12,
           label="isoform bubble", size=3.2, fontface="bold.italic", color="gray15") +
  annotate("text", x=4.0, y=0.96,
           label="Copy A: two transcripts\n(one skips this exon)", size=2.7, color="gray45") +
  # Dividing line between bubbles
  annotate("segment", x=3.5, xend=3.5, y=-1.0, yend=0.88,
           linetype="dotted", color="gray55", linewidth=0.8) +
  # Validation warning
  annotate("label", x=3.0, y=-0.88,
           label="Both look structurally identical in the VG\nValidation required per family",
           fill="#FFFDE7", color="gray30", size=2.7,
           label.padding=unit(0.22,"lines"), hjust=0.5) +
  coord_cartesian(ylim=c(-1.08, 1.25)) +
  labs(title="B  VG bubble types: copy divergence vs isoform variation") +
  theme_graph(base_size=SLIDE_BASE) +
  theme(plot.title=element_text(face="bold", size=10.5, margin=margin(b=4)))

# ═══════════════════════════════════════════════════════════════════════════════
# PANEL C — EM: copy-specific reads establish theta, then distribute ambiguous reads
# ═══════════════════════════════════════════════════════════════════════════════
# Concrete numbers:
#   Copy-specific uniquely-mapped reads: A=6, B=3, C=1  -> theta_A=0.60, B=0.30, C=0.10
#   Ambiguous reads in shared exon:      4 reads, equal MAPQ to all 3
#   EM distributes: A gets +2.4, B gets +1.2, C gets +0.4

n_uniq  <- c(A=6, B=3, C=1)
theta   <- n_uniq / sum(n_uniq)
n_ambig <- 4
n_em    <- n_ambig * theta  # A=2.4, B=1.2, C=0.4

em_df <- data.frame(
  copy    = factor(c("A","B","C"), levels=c("A","B","C")),
  unique  = as.numeric(n_uniq),
  em_part = as.numeric(n_em),
  fill    = unname(copy_cols[c("A","B","C")])
)

# Stacked bar: unique (solid) + EM fractional (lighter) per copy
em_long <- rbind(
  data.frame(copy=em_df$copy, reads=em_df$unique,
             src=factor("copy-specific (assigned directly)",
                        levels=c("copy-specific (assigned directly)",
                                 "shared region (EM fractional)")),
             fill=em_df$fill, stringsAsFactors=FALSE),
  data.frame(copy=em_df$copy, reads=em_df$em_part,
             src=factor("shared region (EM fractional)",
                        levels=c("copy-specific (assigned directly)",
                                 "shared region (EM fractional)")),
             fill=em_df$fill, stringsAsFactors=FALSE)
)

# Label positions for text inside bars
em_df$y_uniq_mid <- em_df$unique / 2
em_df$y_em_mid   <- em_df$unique + em_df$em_part / 2
em_df$y_top      <- em_df$unique + em_df$em_part + 0.25
em_df$theta_lab  <- sprintf("theta=%.2f", as.numeric(theta))

p_em <- ggplot(em_long, aes(x=copy, y=reads, fill=fill, alpha=src)) +
  geom_col(width=0.55, color="white", linewidth=0.8, position="stack") +
  scale_fill_identity() +
  scale_alpha_manual(
    values=c("copy-specific (assigned directly)"=1.00,
             "shared region (EM fractional)"=0.40),
    name=NULL,
    guide=guide_legend(override.aes=list(fill="gray50"))) +
  # Unique-read labels
  geom_text(data=em_df, aes(x=copy, y=y_uniq_mid,
                              label=paste0(unique," reads\n(copy-specific)")),
            size=2.9, color="white", fontface="bold", inherit.aes=FALSE) +
  # EM-read labels
  geom_text(data=em_df, aes(x=copy, y=y_em_mid,
                              label=paste0("+",em_part," (EM)")),
            size=2.7, color="gray25", inherit.aes=FALSE) +
  # theta labels above bars
  geom_text(data=em_df, aes(x=copy, y=y_top, label=theta_lab, color=fill),
            size=3.0, fontface="bold", inherit.aes=FALSE) +
  scale_color_identity() +
  # Ambiguous reads annotation
  annotate("label", x=3.62, y=7.0,
           label="4 reads map to shared exon\n(MAPQ=0 to all 3 copies)\n-> EM distributes by theta",
           fill="gray95", color="gray30", size=2.7,
           label.padding=unit(0.22,"lines"), hjust=0.5) +
  annotate("segment", x=3.25, xend=2.8, y=6.1, yend=5.2,
           arrow=arrow(length=unit(2,"mm"), type="closed"),
           color="gray40", linewidth=0.8) +
  # theta derivation note
  annotate("text", x=0.55, y=9.4,
           label="theta derived from 10 copy-specific reads; 4 ambiguous reads then split by theta",
           size=2.6, color="gray45", hjust=0, fontface="italic") +
  scale_x_discrete(labels=c(A="Copy A\n(LOC389895)",
                              B="Copy B\n(LOC729966)",
                              C="Copy C\n(LOC101928837)")) +
  scale_y_continuous(limits=c(0, 10.0), breaks=seq(0,10,2)) +
  labs(title="C  EM assigns ambiguous reads: copy-specific reads anchor theta",
       x=NULL, y="Read count (fractional for EM-assigned)") +
  theme_classic(base_size=SLIDE_BASE) +
  theme(axis.text.x=element_text(size=8.5),
        legend.position="bottom",
        legend.text=element_text(size=9),
        plot.title=element_text(face="bold", size=10.5))

# ═══════════════════════════════════════════════════════════════════════════════
# Combine
# ═══════════════════════════════════════════════════════════════════════════════
fig5 <- p_def + p_vg + p_em +
  plot_layout(ncol=3, widths=c(1.55, 1.05, 1.0)) +
  plot_annotation(
    title="Multi-copy gene family: definition, VG structure, and read assignment",
    theme=theme(plot.title=element_text(size=14, face="bold"))
  )

ggsave(file.path(FIGURES_DIR, "fig_hmm5_definition_em.pdf"), fig5,
       width=18, height=7, units="in")
ggsave(file.path(FIGURES_DIR, "fig_hmm5_definition_em.png"), fig5,
       width=18, height=7, units="in", dpi=180)
message("Figure 5 saved.")
