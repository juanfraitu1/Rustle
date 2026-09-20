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

DATA_DIR    <- "data"
FIGURES_DIR <- "figures"
GENOME_FA   <- "/mnt/c/Users/jfris/Desktop/GGO.fasta"

vg <- readRDS(file.path(DATA_DIR, "golga_GOLGA8I_vg_seq.rds"))
intervals <- vg$intervals
nodes     <- vg$nodes

genome <- FaFile(GENOME_FA)
open(genome)

# ── 1. Build per-copy canonical transcript sequence ──────────────────────────
get_canonical_tx <- function(gene_id) {
  gene_txs <- unique(vg$exon_df$tx_id[vg$exon_df$gene_id == gene_id])
  exon_counts <- sapply(gene_txs, function(tx) sum(vg$exon_df$tx_id == tx))
  gene_txs[which.max(exon_counts)]
}

extract_copy_seq <- function(gene_id) {
  tx <- get_canonical_tx(gene_id)
  ex <- vg$exon_df[vg$exon_df$tx_id == tx, ]
  ex <- ex[order(ex$start), ]
  gr <- GRanges(seqnames = ex$chrom,
                ranges   = IRanges(ex$start, ex$end),
                strand   = ex$strand)
  seqs <- scanFa(genome, gr)
  # 5' -> 3' orientation: if negative strand, reverse complement and reverse order
  if (as.character(ex$strand[1]) == "-") {
    seqs <- rev(reverseComplement(seqs))
  }
  list(seq = unlist(seqs),
       tx = tx,
       exons = ex,
       gr = gr,
       strand = as.character(ex$strand[1]))
}

copy_ids <- names(vg$paths)
copies <- setNames(lapply(copy_ids, extract_copy_seq), copy_ids)

# Compact display labels (last 6 digits of LOC ID)
short_label <- function(g) sub("^gene-LOC", "", g)
labels <- setNames(short_label(copy_ids), copy_ids)

# ── 2. Resolve node -> per-copy transcript position ───────────────────────────
# For each (copy, node), find the start/end positions within the concatenated
# transcript sequence (5'->3'). Used to know where each "bubble exon" sits in
# each copy's read-able sequence.
copy_node_pos <- list()
for (gid in copy_ids) {
  cp <- copies[[gid]]
  pos <- 1L
  rows <- list()
  # ex is already sorted by genomic start; for negative strand we need 3'->5'
  ex <- cp$exons
  ex <- ex[order(ex$start), ]
  if (cp$strand == "-") ex <- ex[order(-ex$start), ]
  for (i in seq_len(nrow(ex))) {
    w <- ex$end[i] - ex$start[i] + 1L
    # look up node id for this exon (from intervals)
    nid <- intervals$node_id[
      intervals$gene_id == gid & intervals$start == ex$start[i] &
      intervals$end == ex$end[i]
    ][1]
    rows[[i]] <- data.frame(gene_id = gid, exon_idx = i,
                            node_id = nid,
                            tx_start = pos, tx_end = pos + w - 1L,
                            width = w)
    pos <- pos + w
  }
  copy_node_pos[[gid]] <- do.call(rbind, rows)
}

# ── 3. Simulate reads from chosen sources ────────────────────────────────────
# Strategy: pick well-defined source regions in specific copies, sample 150 bp
# centered on the region, add 2% random substitution noise to mimic IsoSeq HiFi.
set.seed(42)
ERROR_RATE <- 0.02
mutate_seq <- function(s, rate = ERROR_RATE) {
  chars <- strsplit(s, "")[[1]]
  alt   <- c("A","C","G","T")
  mut_pos <- sample.int(length(chars), size = round(length(chars) * rate))
  for (p in mut_pos) {
    chars[p] <- sample(setdiff(alt, chars[p]), 1)
  }
  paste0(chars, collapse = "")
}

sample_read_around_node <- function(gene_id, target_node, read_len = 150L) {
  cp_pos <- copy_node_pos[[gene_id]]
  hit <- cp_pos[cp_pos$node_id == target_node, ]
  if (nrow(hit) == 0L) return(NULL)
  # center read on node's middle
  midpoint <- round((hit$tx_start + hit$tx_end) / 2)
  start <- max(1L, midpoint - read_len %/% 2L)
  end   <- min(length(copies[[gene_id]]$seq), start + read_len - 1L)
  if (end - start + 1 < read_len) start <- max(1L, end - read_len + 1L)
  raw <- as.character(subseq(copies[[gene_id]]$seq, start, end))
  list(read_seq = mutate_seq(raw),
       source_gene = gene_id, source_node = target_node,
       tx_start = start, tx_end = end)
}

# Use realistic IsoSeq long reads: full-length transcript + a few short controls.
# IsoSeq HiFi reads are typically full-length cDNAs (1-4 kb in GOLGA8I).
#
# R1-R3 : FULL transcript reads from three copies — the long-read scenario
# R4    : 200 bp read from a fully-shared backbone exon (node 16, n=7)
# R5    : 800 bp read spanning ~3 exons (some shared, one bubble branch)
sample_full_transcript <- function(gene_id) {
  list(read_seq = mutate_seq(as.character(copies[[gene_id]]$seq)),
       source_gene = gene_id, source_node = NA_integer_,
       tx_start = 1L, tx_end = length(copies[[gene_id]]$seq),
       len = length(copies[[gene_id]]$seq))
}

sample_window <- function(gene_id, tx_start, len) {
  L <- length(copies[[gene_id]]$seq)
  start <- max(1L, min(tx_start, L - len + 1L))
  end   <- start + len - 1L
  raw <- as.character(subseq(copies[[gene_id]]$seq, start, end))
  list(read_seq = mutate_seq(raw),
       source_gene = gene_id,
       source_node = NA_integer_,
       tx_start = start, tx_end = end, len = len)
}

# Locate node 16 (shared backbone, n=7) start position in LOC101150678 for R4
n16_rows <- copy_node_pos[["gene-LOC101150678"]]
n16_start_678 <- n16_rows$tx_start[n16_rows$node_id == 16][1]
# Locate node 2 start in LOC101150678 for R5 (covers ~end-of-5 + 2 + start-of-7+12)
n2_start_678  <- n16_rows$tx_start[n16_rows$node_id == 2][1]

read_specs <- list(
  R1 = list(gene = "gene-LOC101150678", kind = "full",
             tag = "R1 (full IsoSeq tx from copy 678)"),
  R2 = list(gene = "gene-LOC129527057", kind = "full",
             tag = "R2 (full IsoSeq tx from copy 057)"),
  R3 = list(gene = "gene-LOC129527170", kind = "full",
             tag = "R3 (full IsoSeq tx from copy 170)"),
  R4 = list(gene = "gene-LOC101150678", kind = "window",
             tx_start = n16_start_678, len = 200L,
             tag = "R4 (200bp on shared node 16)"),
  R5 = list(gene = "gene-LOC101150678", kind = "window",
             tx_start = max(1L, n2_start_678 - 250L), len = 800L,
             tag = "R5 (800bp spanning bubble of 678)")
)

reads <- lapply(read_specs, function(s) {
  if (s$kind == "full") sample_full_transcript(s$gene)
  else                  sample_window(s$gene, s$tx_start, s$len)
})

# ── 4. Score each read against each copy's HMM ────────────────────────────────
# Simple HMM = position-independent emission model with p_match per position.
# logL(read | copy) = best sliding-window log-likelihood on copy's sequence.
P_MATCH    <- 0.95
LOG_MATCH  <- log(P_MATCH)
LOG_MISS   <- log((1 - P_MATCH) / 3)

SUB_MAT <- pwalign::nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = TRUE)

# Score a read against a copy using Smith-Waterman local alignment.
# Best contiguous match gives matches (M) and aligned length (L_aln).
# logL ≈ M * log(p_match) + (read_len - M) * log((1-p_match)/3)
# (unaligned read positions are treated as mismatches.)
score_read_vs_copy <- function(read_seq, copy_seq) {
  aln <- pwalign::pairwiseAlignment(
    pattern = DNAString(read_seq),
    subject = copy_seq,
    type = "local",
    substitutionMatrix = SUB_MAT,
    gapOpening = -2, gapExtension = -1
  )
  aln_len  <- nchar(as.character(alignedPattern(aln)))
  pat_chr  <- strsplit(as.character(alignedPattern(aln)), "")[[1]]
  sub_chr  <- strsplit(as.character(alignedSubject(aln)), "")[[1]]
  matches  <- sum(pat_chr == sub_chr & pat_chr != "-")
  L_read   <- nchar(read_seq)
  mismatches <- L_read - matches
  matches * LOG_MATCH + mismatches * LOG_MISS
}

cat("Scoring", length(reads), "reads x", length(copies), "copies...\n")
score_matrix <- matrix(NA_real_,
                       nrow = length(reads), ncol = length(copies),
                       dimnames = list(names(reads), copy_ids))
for (rn in names(reads)) {
  r <- reads[[rn]]$read_seq
  for (cn in copy_ids) {
    score_matrix[rn, cn] <- score_read_vs_copy(r, copies[[cn]]$seq)
  }
}

# Compute ΔlogL relative to the best copy per read
delta_logL <- score_matrix - apply(score_matrix, 1, max)

# Per-read summary
cat("\nlogL matrix (read x copy):\n")
print(round(score_matrix, 1))
cat("\nΔlogL (= logL - max(logL) per row):\n")
print(round(delta_logL, 1))

# ── 5. Build figure ───────────────────────────────────────────────────────────
plot_df <- as.data.frame(delta_logL) %>%
  mutate(read = rownames(.)) %>%
  pivot_longer(cols = -read, names_to = "copy", values_to = "delta_logL")
plot_df$copy_label <- labels[plot_df$copy]
plot_df$read_tag <- sapply(plot_df$read, function(r) read_specs[[r]]$tag)
plot_df$source_gene <- sapply(plot_df$read, function(r) read_specs[[r]]$gene)
plot_df$source_label <- labels[plot_df$source_gene]
plot_df$is_source <- plot_df$copy == plot_df$source_gene
plot_df$copy_label <- factor(plot_df$copy_label, levels = unname(labels))
plot_df$read_tag <- factor(plot_df$read_tag,
                             levels = sapply(read_specs, function(s) s$tag))

# Heatmap of ΔlogL with copy and source-copy annotations
p_heat <- ggplot(plot_df,
                  aes(x = copy_label, y = read_tag, fill = delta_logL)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = ifelse(delta_logL == 0, "BEST",
                                sprintf("%.0f", delta_logL))),
            color = "white", fontface = "bold", size = 3.2) +
  scale_fill_gradient2(low = "#B71C1C", mid = "#FB8C00", high = "#1565C0",
                        midpoint = -150, name = "DlogL\n(0 = best)") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev) +
  labs(title    = "Per-copy HMM scoring of reads at GOLGA8I",
       subtitle = "Per-copy HMMs built from each copy's transcript sequence (p_match=0.95). 150-bp simulated IsoSeq reads (2% error).",
       x = "GOLGA8I copy (LOC suffix)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(color = "gray40", size = 10),
        axis.text.x   = element_text(face = "bold", size = 10),
        axis.text.y   = element_text(face = "bold", size = 9),
        panel.grid    = element_blank())

# Bar chart per read: assignment confidence (max - second-best ΔlogL)
margin_df <- plot_df %>%
  group_by(read, read_tag, source_label) %>%
  arrange(desc(delta_logL)) %>%
  summarise(best_copy   = copy_label[1],
            second_best = delta_logL[2],
            margin      = -second_best,   # log-odds margin
            confident   = -second_best > 5,
            .groups = "drop")

p_margin <- ggplot(margin_df,
                    aes(x = read_tag, y = margin,
                        fill = ifelse(confident, "confident", "tied"))) +
  geom_col(width = 0.55, color = "white", linewidth = 0.8) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "gray40") +
  geom_text(aes(label = paste0("best: ", best_copy)),
            y = -3, color = "gray20", size = 3, fontface = "italic",
            angle = 0, hjust = 0.5, vjust = 0) +
  scale_fill_manual(values = c(confident = "#2E7D32", tied = "#9E9E9E"),
                     name = NULL,
                     labels = c(confident = "confident (margin > 5)",
                                tied      = "tied (margin <= 5)")) +
  scale_y_continuous(limits = c(-10, NA)) +
  labs(title = "Assignment margin per read",
       subtitle = "Margin = DlogL gap between best and second-best copy. Margin > 5 = confident assignment.",
       x = NULL, y = "DlogL margin (best vs. 2nd best)") +
  theme_classic(base_size = 11) +
  theme(plot.title       = element_text(face = "bold", size = 12),
        plot.subtitle    = element_text(color = "gray40", size = 9.5),
        axis.text.x      = element_text(angle = 25, hjust = 1, size = 9),
        legend.position  = "right")

fig <- p_heat / p_margin + plot_layout(heights = c(1.4, 1))

ggsave(file.path(FIGURES_DIR, "fig_golga8i_hmm_demo.pdf"), fig,
       width = 13, height = 9, units = "in")
ggsave(file.path(FIGURES_DIR, "fig_golga8i_hmm_demo.png"), fig,
       width = 13, height = 9, units = "in", dpi = 180)
message("HMM-demo figure saved.")

close(genome)
