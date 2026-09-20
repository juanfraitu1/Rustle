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

SOURCE   <- "gene-LOC101150678"
src_seq  <- as.character(copies[[SOURCE]]$seq)
src_chars <- strsplit(src_seq, "")[[1]]

set.seed(7)
ERROR_RATE <- 0.02
mutate_seq <- function(s) {
  ch  <- strsplit(s, "")[[1]]
  alt <- c("A","C","G","T")
  mp  <- sample.int(length(ch), round(length(ch) * ERROR_RATE))
  for (p in mp) ch[p] <- sample(setdiff(alt, ch[p]), 1)
  paste(ch, collapse = "")
}
read_seq  <- mutate_seq(src_seq)
read_chars <- strsplit(read_seq, "")[[1]]

# ── Pairwise-align each copy to source; get the copy's base at every source pos ─
SUB_MAT <- pwalign::nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = TRUE)
get_copy_bases_at_src <- function(copy_seq) {
  aln <- pwalign::pairwiseAlignment(
    pattern = DNAString(src_seq), subject = DNAString(copy_seq),
    type = "global", substitutionMatrix = SUB_MAT,
    gapOpening = -2, gapExtension = -1
  )
  ps <- strsplit(as.character(alignedPattern(aln)), "")[[1]]
  ss <- strsplit(as.character(alignedSubject(aln)), "")[[1]]
  out <- rep("-", nchar(src_seq))
  src_idx <- cumsum(ps != "-")
  mask <- ps != "-"           # only positions where source has a base
  out[src_idx[mask]] <- ss[mask]
  out
}

cat("Computing per-copy bases at all source positions ...\n")
base_mat <- matrix("-", nrow = length(copy_ids), ncol = nchar(src_seq))
rownames(base_mat) <- copy_ids
base_mat[SOURCE, ] <- src_chars
for (g in setdiff(copy_ids, SOURCE)) {
  base_mat[g, ] <- get_copy_bases_at_src(copies[[g]]$seq)
}

# Variant positions: at least one non-gap copy has a different base than source
# AND no gaps among the 7 copies at that position (clean SNVs only)
is_variant <- function(p) {
  bases <- base_mat[, p]
  if (any(bases == "-")) return(FALSE)
  length(unique(bases)) > 1L
}
variant_idx <- which(vapply(seq_len(ncol(base_mat)), is_variant, logical(1)))
cat("Found ", length(variant_idx), " clean variant positions across 7 copies\n")

# Pick N representative positions evenly spread across the transcript,
# preferring positions where MAX copies have non-source bases (more discriminative).
N_DISPLAY <- 24
# Score each variant position by Shannon entropy of the 7-allele distribution
entropy_of <- function(p) {
  bases <- base_mat[, p]
  tab <- prop.table(table(bases))
  -sum(tab * log(tab))
}
# Bin transcript into N_DISPLAY equal bins, take highest-entropy variant per bin
bins <- cut(variant_idx, breaks = N_DISPLAY, labels = FALSE)
display_idx <- sapply(seq_len(N_DISPLAY), function(b) {
  cand <- variant_idx[bins == b]
  if (length(cand) == 0L) return(NA_integer_)
  ent <- vapply(cand, entropy_of, numeric(1))
  cand[which.max(ent)]
})
display_idx <- display_idx[!is.na(display_idx)]
cat("Displaying ", length(display_idx), " variant positions\n")

# Build the fingerprint matrix
ordered_copies <- c(SOURCE, setdiff(copy_ids, SOURCE))
fp_mat <- rbind(
  READ = read_chars[display_idx],
  base_mat[ordered_copies, display_idx]
)
rownames(fp_mat) <- c("READ", labels[ordered_copies])
colnames(fp_mat) <- as.character(display_idx)

# Long form for ggplot
fp_long <- as.data.frame(fp_mat) %>%
  tibble::rownames_to_column("row") %>%
  pivot_longer(cols = -row, names_to = "src_pos", values_to = "base")
fp_long$src_pos  <- as.integer(fp_long$src_pos)
fp_long$row      <- factor(fp_long$row, levels = rev(c("READ", labels[ordered_copies])))

# Mark cells whose base matches the READ's base at that column
read_base_at <- setNames(fp_mat["READ", ], colnames(fp_mat))
fp_long$matches_read <- fp_long$base == read_base_at[as.character(fp_long$src_pos)]
fp_long$is_read_row  <- fp_long$row == "READ"

# Nucleotide palette
nuc_pal <- c(A = "#43A047", C = "#FB8C00", G = "#E53935", T = "#1E88E5", "-" = "#BDBDBD")
fp_long$fill_col <- nuc_pal[fp_long$base]

# ── Score per row: count cells matching READ ─────────────────────────────────
match_count <- fp_long %>%
  filter(!is_read_row) %>%
  group_by(row) %>%
  summarise(matches = sum(matches_read), .groups = "drop") %>%
  arrange(desc(matches))
cat("\nMatches with READ (at displayed variant positions):\n")
print(match_count, n = Inf)

# Reorder rows: READ at top, then copies sorted by match count desc
ordered_labels <- c("READ", as.character(match_count$row))
fp_long$row <- factor(fp_long$row, levels = rev(ordered_labels))

# ── PANEL A: header ───────────────────────────────────────────────────────────
p_A <- ggplot() +
  annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = 1,
           fill = "#F5F5F5", color = "gray60", linewidth = 0.4) +
  annotate("text", x = 5, y = 0.5,
           label = paste0(
             "Read R = ", labels[SOURCE], " transcript (", nchar(src_seq),
             " bp) with 2% noise.    ",
             length(variant_idx),
             " clean variant positions across 7 copies — we show ",
             length(display_idx),
             " of them (one per bin, highest-entropy)."),
           hjust = 0.5, size = 3.7, color = "gray20") +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 1)) +
  theme_void()

# ── PANEL B: the fingerprint matrix ──────────────────────────────────────────
# Position labels along the x-axis (one tick per displayed column)
pos_labels <- as.character(display_idx)

p_B <- ggplot(fp_long,
               aes(x = factor(src_pos, levels = display_idx),
                   y = row, fill = fill_col)) +
  geom_tile(color = "white", linewidth = 1.2) +
  scale_fill_identity() +
  # Highlight cells matching the READ (except in the READ row itself)
  geom_tile(data = filter(fp_long, matches_read & !is_read_row),
            aes(x = factor(src_pos, levels = display_idx), y = row),
            fill = NA, color = "#1B5E20", linewidth = 2.0) +
  # Letter inside cell
  geom_text(aes(label = base,
                color = ifelse(matches_read & !is_read_row, "#1B5E20", "white"),
                fontface = ifelse(matches_read & !is_read_row, "bold", "bold")),
            size = 4.2) +
  scale_color_identity() +
  # Black bar separating READ row from copies
  geom_hline(yintercept = length(ordered_labels) - 0.5, color = "black",
             linewidth = 0.8) +
  scale_x_discrete(labels = pos_labels, position = "top") +
  labs(title    = "B  Fingerprint matrix at variant positions",
       subtitle = "Each column = one variant position (bp coordinate in source transcript). Each row = one copy's nucleotide at that position. Green-bordered cells = same base as the READ.",
       x = "Source-transcript position (bp)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title       = element_text(face = "bold", size = 12),
        plot.subtitle    = element_text(color = "gray40", size = 10),
        panel.grid       = element_blank(),
        axis.text.x      = element_text(size = 8, angle = 60, hjust = 0),
        axis.text.y      = element_text(fontface = "bold", size = 10.5))

# ── PANEL C: per-row total matches with READ ─────────────────────────────────
match_df <- match_count %>%
  mutate(row = factor(row, levels = rev(as.character(match_count$row))),
         pct = matches / length(display_idx) * 100)

p_C <- ggplot(match_df, aes(y = row, x = matches,
                              fill = ifelse(matches == max(matches),
                                              "best","other"))) +
  geom_col(width = 0.55, color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%d / %d  (%.0f%%)",
                                  matches, length(display_idx), pct),
                x = matches - 0.3),
            color = "white", fontface = "bold", size = 3.2, hjust = 1) +
  scale_fill_manual(values = c(best = "#1B5E20", other = "#9E9E9E"),
                    guide = "none") +
  scale_x_continuous(limits = c(0, length(display_idx) + 1)) +
  labs(title    = "C  Matches with READ at displayed variant positions",
       subtitle = "Source copy matches the read at every variant position (within noise). No other copy comes close.",
       x = paste0("# matches  (out of ", length(display_idx), " variant positions)"),
       y = NULL) +
  theme_classic(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(color = "gray40", size = 10),
        axis.text.y   = element_text(fontface = "bold", size = 10.5))

# ── Combine ───────────────────────────────────────────────────────────────────
layout_design <- "
AA
BB
BB
BB
CC
"
fig <- p_A + p_B + p_C +
  plot_layout(design = layout_design) +
  plot_annotation(
    title    = "GOLGA8I read-to-copy assignment: variant-position fingerprint",
    subtitle = "The conserved part of the transcript is ignored. Only positions where the 7 copies disagree carry information. The read's base at each variant position 'votes' for one or more copies — the source wins because it matches at (almost) every position.",
    theme = theme(plot.title    = element_text(size = 15, face = "bold"),
                  plot.subtitle = element_text(size = 11, color = "gray40",
                                                lineheight = 1.35))
  )

ggsave(file.path(FIGURES_DIR, "fig_golga8i_fingerprint.pdf"), fig,
       width = 18, height = 10, units = "in")
ggsave(file.path(FIGURES_DIR, "fig_golga8i_fingerprint.png"), fig,
       width = 18, height = 10, units = "in", dpi = 180)
message("Fingerprint figure saved.")
