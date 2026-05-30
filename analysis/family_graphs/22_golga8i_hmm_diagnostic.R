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

# ── Per-copy transcript sequences ────────────────────────────────────────────
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
  list(seq = unlist(seqs), len = sum(width(seqs)))
}
copy_ids <- names(vg$paths)
copies   <- setNames(lapply(copy_ids, extract_copy_seq), copy_ids)
labels   <- setNames(sub("^gene-LOC", "", copy_ids), copy_ids)
close(genome)

# ── Generate a full-length read from LOC101150678 + 2% error ─────────────────
SOURCE <- "gene-LOC101150678"
src_seq <- as.character(copies[[SOURCE]]$seq)
READ_LEN <- nchar(src_seq)

set.seed(7)
ERROR_RATE <- 0.02
mutate_seq <- function(s) {
  ch <- strsplit(s, "")[[1]]
  alt <- c("A","C","G","T")
  mp <- sample.int(length(ch), round(length(ch) * ERROR_RATE))
  for (p in mp) ch[p] <- sample(setdiff(alt, ch[p]), 1)
  paste(ch, collapse = "")
}
read_seq <- mutate_seq(src_seq)
cat("Source:", labels[SOURCE], "  read length:", READ_LEN, "bp\n")

# ── For every non-source copy: pairwise-align to source, find diagnostic cols ─
SUB_MAT <- pwalign::nucleotideSubstitutionMatrix(match = 1, mismatch = -1,
                                                  baseOnly = TRUE)

# Align two strings globally; return per-position information at SOURCE positions
align_to_source <- function(copy_seq) {
  aln <- pwalign::pairwiseAlignment(
    pattern = DNAString(src_seq), subject = DNAString(copy_seq),
    type = "global", substitutionMatrix = SUB_MAT,
    gapOpening = -2, gapExtension = -1
  )
  ps <- strsplit(as.character(alignedPattern(aln)), "")[[1]]
  ss <- strsplit(as.character(alignedSubject(aln)), "")[[1]]
  # Map alignment columns to source positions
  src_idx <- cumsum(ps != "-")
  data.frame(
    src_pos   = src_idx,
    src_base  = ps,
    copy_base = ss,
    is_gap_src  = ps == "-",
    is_gap_copy = ss == "-",
    is_match    = ps == ss & ps != "-",
    is_diag     = ps != ss & ps != "-" & ss != "-",
    is_indel    = ps == "-" | ss == "-",
    stringsAsFactors = FALSE
  )
}

# Pre-compute alignments
cat("Aligning source vs each non-source copy...\n")
non_source_ids <- setdiff(copy_ids, SOURCE)
aln_per_copy   <- lapply(non_source_ids, function(g) align_to_source(copies[[g]]$seq))
names(aln_per_copy) <- non_source_ids

# Read alignment in source coordinates (read covers full source with 2% noise,
# so read position i corresponds to source position i for any non-shifted region)
read_chars <- strsplit(read_seq, "")[[1]]
src_chars  <- strsplit(src_seq, "")[[1]]

# For each non-source copy, classify each diagnostic position:
#   support_src  = read base equals source base (favours source)
#   support_copy = read base equals copy base   (favours that copy)
#   neither      = read base is different from both (read error not on diag)
classify_evidence <- function(aln_df) {
  diag_rows <- aln_df[aln_df$is_diag & aln_df$src_pos <= length(read_chars), ]
  if (nrow(diag_rows) == 0L) return(NULL)
  rb <- read_chars[diag_rows$src_pos]
  diag_rows$read_base <- rb
  diag_rows$evidence <- with(diag_rows,
    ifelse(read_base == src_base,  "supports_source",
    ifelse(read_base == copy_base, "supports_copy", "neither")))
  diag_rows
}

evidence <- lapply(aln_per_copy, classify_evidence)
names(evidence) <- non_source_ids

# Summarise: counts of evidence types per non-source copy
ev_summary <- do.call(rbind, lapply(non_source_ids, function(g) {
  ev <- evidence[[g]]
  total_diag <- nrow(ev)
  data.frame(
    copy = g,
    copy_label = labels[g],
    n_diag_positions = total_diag,
    n_indel = sum(aln_per_copy[[g]]$is_indel),
    supports_source = sum(ev$evidence == "supports_source"),
    supports_copy   = sum(ev$evidence == "supports_copy"),
    neither         = sum(ev$evidence == "neither"),
    stringsAsFactors = FALSE
  )
}))
ev_summary <- ev_summary %>%
  mutate(pct_source = round(100 * supports_source / n_diag_positions, 1),
         pct_copy   = round(100 * supports_copy   / n_diag_positions, 1))

# Also include the source itself for context
self_row <- data.frame(copy = SOURCE, copy_label = labels[SOURCE],
                       n_diag_positions = 0L, n_indel = 0L,
                       supports_source = 0L, supports_copy = 0L, neither = 0L,
                       pct_source = NA_real_, pct_copy = NA_real_)
ev_summary <- bind_rows(self_row, ev_summary)

cat("\nEvidence summary at diagnostic (variant) positions:\n")
print(ev_summary, row.names = FALSE)

# ── PANEL A: source vs reads — context block ──────────────────────────────────
p_A <- ggplot() +
  annotate("rect", xmin = 0, xmax = nchar(src_seq), ymin = 0.45, ymax = 0.85,
           fill = "#e6a817", color = "white", linewidth = 0.4) +
  annotate("text", x = nchar(src_seq)/2, y = 0.65,
           label = sprintf("%s transcript  (%d bp)  —  ~5290 bp mostly conserved",
                            labels[SOURCE], nchar(src_seq)),
           color = "white", fontface = "bold", size = 3.7) +
  annotate("rect", xmin = 0, xmax = nchar(src_seq), ymin = 1.05, ymax = 1.45,
           fill = "#212121", color = "black", linewidth = 0.6) +
  annotate("text", x = nchar(src_seq)/2, y = 1.25,
           label = sprintf("Read R  =  %s transcript  +  2%% substitution noise",
                            labels[SOURCE]),
           color = "white", fontface = "bold", size = 3.7) +
  coord_cartesian(ylim = c(0.25, 1.65)) +
  labs(title    = "A  The read",
       subtitle = "Most of these 5290 bp are identical across all 7 GOLGA8I copies. Skip them — only the diagnostic columns below contribute to assignment.") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(color = "gray40", size = 10))

# ── PANEL B: compressed view of diagnostic positions per non-source copy ─────
# Each row is one non-source copy. Ticks at diagnostic positions, colored by
# whether read supports SOURCE (green) or this COPY (red) or neither (gray).
all_evidence <- bind_rows(lapply(non_source_ids, function(g) {
  ev <- evidence[[g]]
  if (is.null(ev) || nrow(ev) == 0L) return(NULL)
  data.frame(copy = g, copy_label = labels[g],
              src_pos = ev$src_pos, evidence = ev$evidence,
              stringsAsFactors = FALSE)
}))

# Order copies by support_source descending (most-distinguishable last)
order_copies <- ev_summary %>%
  filter(copy != SOURCE) %>%
  arrange(desc(pct_source))
all_evidence$copy_label <- factor(all_evidence$copy_label,
                                    levels = order_copies$copy_label)

evidence_pal <- c(supports_source = "#1B5E20",
                  supports_copy   = "#B71C1C",
                  neither         = "#9E9E9E")

p_B <- ggplot(all_evidence,
               aes(x = src_pos, y = copy_label, color = evidence)) +
  geom_point(shape = "|", size = 2.5, alpha = 0.7) +
  scale_color_manual(values = evidence_pal,
                     name = NULL,
                     labels = c(supports_source = "read base = source base (supports source)",
                                 supports_copy   = "read base = this copy's base (supports copy)",
                                 neither         = "read base = neither (read error at variant)")) +
  coord_cartesian(xlim = c(1, nchar(src_seq))) +
  labs(title    = "B  Diagnostic positions per copy",
       subtitle = "Each tick = one position where this copy's transcript differs from source. Conserved positions are NOT shown.",
       x = "Position in source transcript (bp)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title       = element_text(face = "bold", size = 12),
        plot.subtitle    = element_text(color = "gray40", size = 10),
        panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
        panel.grid.minor   = element_blank(),
        axis.text.y      = element_text(face = "bold", size = 10),
        legend.position  = "bottom")

# ── PANEL C: stacked bar of evidence counts per copy ────────────────────────
bar_df <- ev_summary %>%
  filter(copy != SOURCE) %>%
  select(copy_label, supports_source, supports_copy, neither) %>%
  pivot_longer(cols = -copy_label, names_to = "evidence", values_to = "count") %>%
  mutate(evidence = factor(evidence,
    levels = c("neither","supports_copy","supports_source"))) %>%
  mutate(copy_label = factor(copy_label, levels = order_copies$copy_label))

p_C <- ggplot(bar_df, aes(y = copy_label, x = count, fill = evidence)) +
  geom_col(position = "stack", color = "white", linewidth = 0.4) +
  geom_text(data = ev_summary %>% filter(copy != SOURCE) %>%
              mutate(copy_label = factor(copy_label, levels = order_copies$copy_label)),
            aes(y = copy_label,
                x = supports_source + supports_copy + neither + 5,
                label = sprintf("%d diag,  %d supports source vs %d supports this",
                                  n_diag_positions, supports_source, supports_copy)),
            color = "gray30", size = 3.0, hjust = 0, inherit.aes = FALSE) +
  scale_fill_manual(values = evidence_pal, guide = "none") +
  labs(title    = "C  Evidence count at diagnostic positions",
       subtitle = "All 6 non-source copies have MANY diagnostic positions, and read R supports SOURCE at almost all of them.",
       x = "Number of diagnostic positions", y = NULL) +
  theme_classic(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(color = "gray40", size = 10),
        axis.text.y   = element_text(face = "bold", size = 10))

# ── PANEL D: final DlogL ──────────────────────────────────────────────────────
# logL at diagnostic positions only (more interpretable):
# DlogL approx = (supports_copy - supports_source) * (log(0.0167) - log(0.95))
P_MATCH   <- 0.95
LOG_MATCH <- log(P_MATCH)
LOG_MISS  <- log((1 - P_MATCH) / 3)

ev_summary$logL_diag_src  <- ev_summary$supports_source * LOG_MATCH +
                              (ev_summary$supports_copy + ev_summary$neither) * LOG_MISS
ev_summary$logL_diag_self <- ev_summary$supports_copy * LOG_MATCH +
                              (ev_summary$supports_source + ev_summary$neither) * LOG_MISS
ev_summary$delta_at_diag <- ev_summary$logL_diag_src - ev_summary$logL_diag_self

bar_d <- ev_summary %>%
  filter(copy != SOURCE) %>%
  mutate(copy_label = factor(copy_label, levels = rev(order_copies$copy_label)))

p_D <- ggplot(bar_d, aes(y = copy_label, x = delta_at_diag)) +
  geom_col(fill = "#1565C0", width = 0.55, color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("+%.0f", delta_at_diag), x = delta_at_diag - 30),
            color = "white", fontface = "bold", size = 3.2, hjust = 1) +
  geom_vline(xintercept = 0, color = "gray30", linewidth = 0.5) +
  labs(title    = "D  Net log-odds for SOURCE at diagnostic positions",
       subtitle = "Positive = read supports SOURCE more than this copy. Larger = more confident.",
       x = "log-odds in favour of source (logL_src - logL_this_copy)", y = NULL) +
  theme_classic(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(color = "gray40", size = 10),
        axis.text.y   = element_text(face = "bold", size = 10))

# ── Combine with explicit layout ─────────────────────────────────────────────
layout_design <- "
AA
BB
BB
CD
CD
"
fig <- p_A + p_B + p_C + p_D +
  plot_layout(design = layout_design) +
  plot_annotation(
    title    = "How the HMM actually discriminates: only diagnostic positions matter",
    subtitle = "Most of the 5290 bp transcript is the same in all 7 copies. The discrimination comes from a small subset of positions where the transcripts differ.",
    theme = theme(plot.title    = element_text(size = 15, face = "bold"),
                  plot.subtitle = element_text(size = 11, color = "gray40"))
  )

ggsave(file.path(FIGURES_DIR, "fig_golga8i_diagnostic_view.pdf"), fig,
       width = 16, height = 11, units = "in")
ggsave(file.path(FIGURES_DIR, "fig_golga8i_diagnostic_view.png"), fig,
       width = 16, height = 11, units = "in", dpi = 180)
message("Diagnostic-view figure saved.")
