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
  list(seq = unlist(seqs), exons = ex, strand = as.character(ex$strand[1]))
}
copy_ids <- names(vg$paths)
copies   <- setNames(lapply(copy_ids, extract_copy_seq), copy_ids)
labels   <- setNames(sub("^gene-LOC", "", copy_ids), copy_ids)

SOURCE   <- "gene-LOC101150678"
src_seq  <- as.character(copies[[SOURCE]]$seq)

# Per-copy node position (transcript coords) — for sampling reads from specific regions
copy_node_pos <- list()
intervals <- vg$intervals
for (gid in copy_ids) {
  ex <- copies[[gid]]$exons
  ex <- ex[order(ex$start), ]
  if (copies[[gid]]$strand == "-") ex <- ex[order(-ex$start), ]
  pos <- 1L; rows <- list()
  for (i in seq_len(nrow(ex))) {
    w <- ex$end[i] - ex$start[i] + 1L
    nid <- intervals$node_id[intervals$gene_id == gid &
                              intervals$start == ex$start[i] &
                              intervals$end   == ex$end[i]][1]
    rows[[i]] <- data.frame(node_id = nid, tx_start = pos,
                             tx_end = pos + w - 1L, width = w)
    pos <- pos + w
  }
  copy_node_pos[[gid]] <- do.call(rbind, rows)
}
close(genome)

# ── Sample reads at various scales ────────────────────────────────────────────
set.seed(11)
ERROR_RATE <- 0.02
mutate_seq <- function(s) {
  ch  <- strsplit(s, "")[[1]]
  alt <- c("A","C","G","T")
  mp  <- sample.int(length(ch), round(length(ch) * ERROR_RATE))
  for (p in mp) ch[p] <- sample(setdiff(alt, ch[p]), 1)
  paste(ch, collapse = "")
}
sample_window <- function(gene_id, tx_start, len) {
  L <- length(copies[[gene_id]]$seq)
  start <- max(1L, min(tx_start, L - len + 1L))
  raw <- as.character(subseq(copies[[gene_id]]$seq, start, start + len - 1L))
  list(read = mutate_seq(raw), len = len, src = gene_id, tx_start = start,
       tx_end = start + len - 1L)
}

# Locate node 16 (fully shared, n=7) — should be a conserved region
src_pos_df <- copy_node_pos[[SOURCE]]
n16_start  <- src_pos_df$tx_start[src_pos_df$node_id == 16][1]
n2_start   <- src_pos_df$tx_start[src_pos_df$node_id == 2][1]

reads <- list(
  full = list(label = "Full transcript (5288 bp)",
              read = mutate_seq(src_seq),
              len  = nchar(src_seq), src = SOURCE),
  mid  = list(label = "Medium 800 bp (spans bubble of 678)",
              read = sample_window(SOURCE, max(1, n2_start - 250L), 800L)$read,
              len  = 800L, src = SOURCE),
  short_conserved = list(label = "Short 200 bp on shared node 16",
                          read = sample_window(SOURCE, n16_start, 200L)$read,
                          len  = 200L, src = SOURCE),
  short_variant = list(label = "Short 200 bp covering 2-3 variant SNVs",
                        read = sample_window(SOURCE, max(1, n2_start - 80L), 200L)$read,
                        len  = 200L, src = SOURCE)
)

# ── Smith-Waterman alignment for each (read, copy) pair ──────────────────────
# Use match=2, mismatch=-4, gap=-4 — minimap2's default scoring for asm5
SUB_MAT <- pwalign::nucleotideSubstitutionMatrix(match = 2, mismatch = -4,
                                                  baseOnly = TRUE)

align_score <- function(read_seq, copy_seq) {
  aln <- pwalign::pairwiseAlignment(
    pattern = DNAString(read_seq), subject = copy_seq,
    type = "local", substitutionMatrix = SUB_MAT,
    gapOpening = -4, gapExtension = -1
  )
  pat <- as.character(alignedPattern(aln))
  sub <- as.character(alignedSubject(aln))
  pc  <- strsplit(pat, "")[[1]]; sc <- strsplit(sub, "")[[1]]
  matches <- sum(pc == sc & pc != "-")
  list(score    = score(aln),
       sw_score = score(aln),
       matches  = matches,
       mismatches = nchar(read_seq) - matches,
       aln_len  = nchar(pat),
       aln_start = start(aln@subject),
       aln_end   = end(aln@subject))
}

cat("Computing SW alignment scores across reads and copies...\n")
score_records <- list()
for (read_name in names(reads)) {
  r <- reads[[read_name]]
  for (cn in copy_ids) {
    a <- align_score(r$read, copies[[cn]]$seq)
    score_records[[length(score_records) + 1L]] <- data.frame(
      read       = read_name,
      read_label = r$label,
      read_len   = r$len,
      copy       = cn,
      copy_label = labels[cn],
      sw_score   = a$sw_score,
      matches    = a$matches,
      aln_len    = a$aln_len,
      stringsAsFactors = FALSE
    )
  }
}
scores <- do.call(rbind, score_records)

# Compute "tied region": copies whose score is within 5% of the best
scores <- scores %>%
  group_by(read) %>%
  mutate(best_score = max(sw_score),
         score_pct  = sw_score / best_score * 100,
         tied = score_pct >= 95)

cat("\nAlignment-score summary (SW score with minimap2 scoring; match=+2, mismatch=-4):\n\n")
for (read_name in names(reads)) {
  cat("--- ", reads[[read_name]]$label, " ---\n", sep = "")
  d <- scores %>% filter(read == read_name) %>%
    arrange(desc(sw_score)) %>%
    mutate(gap_to_best = sw_score - best_score,
            mapq_proxy = pmin(60, round(-10 * log10(pmax(1e-6,
                              exp((sw_score - best_score) / 5))))))
  print(d[, c("copy_label","sw_score","gap_to_best","matches","aln_len",
              "score_pct","tied")], row.names = FALSE)
  cat("\n")
}

# ── Visualise ─────────────────────────────────────────────────────────────────
scores$read_label <- factor(scores$read_label,
                             levels = sapply(reads, function(r) r$label))
scores$copy_label <- factor(scores$copy_label, levels = unname(labels[copy_ids]))

p_main <- ggplot(scores, aes(x = copy_label, y = sw_score,
                              fill = tied)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.3) +
  geom_text(aes(label = sw_score), vjust = -0.4, size = 2.6, color = "gray20") +
  geom_hline(data = scores %>% group_by(read_label) %>%
                       summarise(best = max(sw_score), .groups = "drop") %>%
                       mutate(thresh = best * 0.95),
             aes(yintercept = thresh),
             linetype = "dashed", color = "#E65100") +
  scale_fill_manual(values = c(`TRUE` = "#E65100", `FALSE` = "#1565C0"),
                     name = NULL,
                     labels = c(`TRUE` = "tied (>=95% of best)",
                                 `FALSE` = "clearly worse")) +
  facet_wrap(~ read_label, scales = "free_y", ncol = 1) +
  labs(title    = "Raw Smith-Waterman alignment scores per copy",
       subtitle = "Scoring: match +2, mismatch -4, gap open -4, gap ext -1 (minimap2 default for asm5).\nDashed orange line = 95% of best score; bars above this line would be 'tied' for a real aligner.",
       x = NULL,
       y = "SW alignment score (raw)") +
  theme_minimal(base_size = 11) +
  theme(plot.title       = element_text(face = "bold", size = 13),
        plot.subtitle    = element_text(color = "gray40", size = 10,
                                         lineheight = 1.3),
        strip.text       = element_text(face = "bold", size = 10),
        axis.text.x      = element_text(size = 9, face = "bold"),
        legend.position  = "top")

ggsave(file.path(FIGURES_DIR, "fig_golga8i_alignment_ties.pdf"), p_main,
       width = 13, height = 11, units = "in")
ggsave(file.path(FIGURES_DIR, "fig_golga8i_alignment_ties.png"), p_main,
       width = 13, height = 11, units = "in", dpi = 180)
message("Alignment-tie figure saved.")
