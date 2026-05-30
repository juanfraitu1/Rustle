suppressPackageStartupMessages({
  library(Biostrings)
  library(pwalign)
  library(Rsamtools)
  library(GenomicRanges)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(tidyr)
  library(scales)
})

source("02_build_graphs.R")

FIGURES_DIR <- "figures"
DATA_DIR    <- "data"
GENOME_FA   <- "/mnt/c/Users/jfris/Desktop/GGO.fasta"

vg <- readRDS(file.path(DATA_DIR, "golga_GOLGA8I_vg_seq.rds"))
genome <- FaFile(GENOME_FA); open(genome)

# ── Build per-copy transcript sequences ──────────────────────────────────────
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
cat("Source transcript length:", READ_LEN, "bp\n")

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

# ── Align read to all 7 copies (Smith-Waterman local) ────────────────────────
SUB_MAT <- pwalign::nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = TRUE)
align_one <- function(rseq, cseq) {
  aln <- pwalign::pairwiseAlignment(
    pattern = DNAString(rseq), subject = cseq,
    type = "local", substitutionMatrix = SUB_MAT,
    gapOpening = -2, gapExtension = -1
  )
  pat <- as.character(alignedPattern(aln))
  sub <- as.character(alignedSubject(aln))
  pc  <- strsplit(pat, "")[[1]]; sc <- strsplit(sub, "")[[1]]
  matches <- sum(pc == sc & pc != "-")
  list(sub_start = start(aln@subject),
       sub_end   = end(aln@subject),
       aln_len   = nchar(pat),
       matches   = matches,
       mismatches = nchar(rseq) - matches)
}

P_MATCH   <- 0.95
LOG_MATCH <- log(P_MATCH)
LOG_MISS  <- log((1 - P_MATCH) / 3)
score_logL <- function(M, L) M * LOG_MATCH + (L - M) * LOG_MISS

cat("Aligning full read to 7 copies...\n")
alignments <- setNames(lapply(copy_ids, function(g) align_one(read_seq, copies[[g]]$seq)),
                        copy_ids)

score_df <- data.frame(
  copy       = copy_ids,
  copy_label = labels[copy_ids],
  copy_len   = sapply(copies, function(c) c$len),
  sub_start  = sapply(alignments, function(a) a$sub_start),
  sub_end    = sapply(alignments, function(a) a$sub_end),
  matches    = sapply(alignments, function(a) a$matches),
  mismatches = sapply(alignments, function(a) a$mismatches),
  stringsAsFactors = FALSE
)
score_df$logL  <- score_logL(score_df$matches, READ_LEN)
score_df$delta <- score_df$logL - max(score_df$logL)
score_df       <- score_df %>% arrange(desc(logL))
score_df$y_rank <- seq_len(nrow(score_df))

cat("\nPer-copy alignment summary:\n")
print(score_df[, c("copy_label","copy_len","sub_start","sub_end",
                    "matches","mismatches","logL","delta")],
      row.names = FALSE, digits = 4)

# ── PANEL A: the read in its source copy ─────────────────────────────────────
p_A <- ggplot() +
  # source transcript bar
  annotate("rect", xmin = 1, xmax = nchar(src_seq), ymin = 0.4, ymax = 0.8,
           fill = "#e6a817", color = "white", linewidth = 0.5) +
  annotate("text", x = nchar(src_seq)/2, y = 0.6,
           label = sprintf("%s transcript  (%d bp)", labels[SOURCE], nchar(src_seq)),
           color = "white", fontface = "bold", size = 4) +
  # read bar (identical span, slightly offset)
  annotate("rect", xmin = 1, xmax = nchar(src_seq), ymin = 1.05, ymax = 1.45,
           fill = "#212121", color = "black", linewidth = 0.7) +
  annotate("text", x = nchar(src_seq)/2, y = 1.25,
           label = sprintf("Read R  =  %s tx  +  %.0f%% substitution noise",
                            labels[SOURCE], 100 * ERROR_RATE),
           color = "white", fontface = "bold", size = 4) +
  coord_cartesian(ylim = c(0.2, 1.7)) +
  labs(title    = "A  The read",
       subtitle = paste0("Full-length IsoSeq-style read sampled from one copy and ",
                          "perturbed at 2% positions. We will score this read against ",
                          "every copy's HMM.")) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(color = "gray40", size = 10))

# ── PANEL B: contestants — where R aligns in each copy ───────────────────────
# Build a single dataframe with copy tracks
track_df <- do.call(rbind, lapply(score_df$copy, function(g) {
  data.frame(copy = g, y_rank = score_df$y_rank[score_df$copy == g],
             copy_len = score_df$copy_len[score_df$copy == g])
}))

# alignment-window boxes
align_df <- score_df %>%
  mutate(label = sprintf("aln tx[%d-%d],  M=%d / %d,  X=%d",
                          sub_start, sub_end, matches, READ_LEN, mismatches))

p_B <- ggplot() +
  # transcript bar for each copy (length-proportional)
  geom_rect(data = track_df,
            aes(xmin = 1, xmax = copy_len,
                ymin = y_rank - 0.18, ymax = y_rank + 0.18),
            fill = "#e6a817", color = "white", linewidth = 0.4) +
  # best alignment window of R inside that copy
  geom_rect(data = align_df,
            aes(xmin = sub_start, xmax = sub_end,
                ymin = y_rank - 0.3, ymax = y_rank + 0.3),
            fill = NA, color = "#212121", linewidth = 1.0) +
  # left label: copy + alignment stats
  geom_text(data = align_df,
            aes(x = -300, y = y_rank,
                label = sprintf("%s\n%s",
                                 copy_label,
                                 sprintf("aln len %d, M=%d, X=%d",
                                           sub_end - sub_start + 1,
                                           matches, mismatches)),
                color = ifelse(delta == 0, "best", "other")),
            hjust = 1, size = 3.0, fontface = "bold", lineheight = 1.05) +
  # right label: logL + DlogL
  geom_text(data = align_df,
            aes(x = max(track_df$copy_len) + 200, y = y_rank,
                label = sprintf("logL = %.0f\nDlogL = %.0f", logL, delta),
                color = ifelse(delta == 0, "best", "other")),
            hjust = 0, size = 3.0, fontface = "bold", lineheight = 1.05) +
  scale_color_manual(values = c(best = "#1B5E20", other = "gray35"),
                     guide = "none") +
  scale_y_continuous(breaks = NULL, expand = c(0.08, 0.5)) +
  coord_cartesian(xlim = c(-2400, max(track_df$copy_len) + 1900),
                  clip = "off") +
  labs(title    = "B  The 7 contestants — where R best aligns in every copy",
       subtitle = "Each row = one GOLGA8I copy. Yellow bar = full transcript. Black box = read R's best local alignment.\nM = matches, X = mismatches (out of L = read length).",
       x = "Position in copy transcript (bp)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(color = "gray40", size = 9.5, lineheight = 1.2),
        panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
        panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_blank())

# ── PANEL C: the math, step by step ──────────────────────────────────────────
best <- score_df[1, ]; second <- score_df[2, ]; worst <- score_df[nrow(score_df), ]

math_text <- paste0(
  "STEP 1 — Per-position emission model:\n",
  "  P(read base = copy base | match) = p_match = 0.95\n",
  "  P(other base | mismatch)         = (1 - p_match) / 3 = 0.0167\n\n",
  "STEP 2 — Take logs:\n",
  "  log(0.95)   = -0.0513\n",
  "  log(0.0167) = -4.0942\n\n",
  "STEP 3 — Score per copy (L = read length):\n",
  "  logL(R | copy) = M . (-0.0513) + (L - M) . (-4.0942)\n\n",
  "STEP 4 — Apply to GOLGA8I (L = ", READ_LEN, "):\n",
  sprintf("  best  : %s   M=%d  ->  logL = %.1f\n",
          best$copy_label, best$matches, best$logL),
  sprintf("  2nd   : %s   M=%d  ->  logL = %.1f  (DlogL = %.1f)\n",
          second$copy_label, second$matches, second$logL, second$delta),
  sprintf("  worst : %s   M=%d  ->  logL = %.1f  (DlogL = %.1f)\n",
          worst$copy_label, worst$matches, worst$logL, worst$delta),
  "\n",
  "Differences come from many small ones in 'shared' nodes (85-99% identity\n",
  "between copies). They accumulate over the full read length."
)

p_C <- ggplot() +
  annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = 10,
           fill = "#F8F8F8", color = "gray60", linewidth = 0.5) +
  annotate("text", x = 0.3, y = 5, label = math_text,
           hjust = 0, vjust = 0.5, size = 3.3, family = "mono",
           color = "gray15", lineheight = 1.35) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 10)) +
  labs(title = "C  How each number is computed") +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0,
                                   margin = margin(b = 6)))

# ── PANEL D: DlogL ranking bar chart ─────────────────────────────────────────
copy_pal <- setNames(
  c("#1565C0","#B71C1C","#2E7D32","#6A1B9A","#EF6C00","#00838F","#5D4037"),
  copy_ids
)
score_df$bar_color  <- copy_pal[score_df$copy]
score_df$copy_label <- factor(score_df$copy_label, levels = rev(score_df$copy_label))

p_D <- ggplot(score_df, aes(y = copy_label, x = delta, fill = bar_color)) +
  geom_col(width = 0.55, color = "white", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "gray30", linewidth = 0.5) +
  geom_text(aes(label = ifelse(delta == 0, "BEST",
                                  sprintf("%.0f", delta)),
                x = pmin(delta - 15, -15)),
            color = "white", fontface = "bold", size = 3.2, hjust = 1) +
  scale_fill_identity() +
  labs(title    = "D  DlogL ranking",
       subtitle = "0 = best assignment, negatives = worse",
       x = "DlogL (read vs. copy's HMM)", y = NULL) +
  theme_classic(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(color = "gray45", size = 9.5),
        axis.text.y   = element_text(face = "bold", size = 10))

# ── Combine ──────────────────────────────────────────────────────────────────
layout_design <- "
AAAA
BBBB
BBBB
BBBB
CCDD
CCDD
"
fig <- p_A + p_B + p_C + p_D +
  plot_layout(design = layout_design) +
  plot_annotation(
    title    = "How a read gets assigned to a GOLGA8I copy — step by step",
    subtitle = "Real LOC101150678 transcript sampled with 2% noise, scored against all 7 copies' HMMs.",
    theme = theme(plot.title    = element_text(size = 15, face = "bold"),
                  plot.subtitle = element_text(size = 11, color = "gray40"))
  )

ggsave(file.path(FIGURES_DIR, "fig_golga8i_hmm_explainer.pdf"), fig,
       width = 16, height = 11.5, units = "in")
ggsave(file.path(FIGURES_DIR, "fig_golga8i_hmm_explainer.png"), fig,
       width = 16, height = 11.5, units = "in", dpi = 180)
message("Explainer figure saved.")
