suppressPackageStartupMessages({
  library(Biostrings)
  library(pwalign)
  library(Rsamtools)
  library(GenomicRanges)
  library(GenomicAlignments)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(tidyr)
})

source("02_build_graphs.R")

FIGURES_DIR <- "figures"
DATA_DIR    <- "data"
BAM_FILE    <- "/mnt/c/Users/jfris/Desktop/GGO.bam"
GENOME_FA   <- "/mnt/c/Users/jfris/Desktop/GGO.fasta"

vg <- readRDS(file.path(DATA_DIR, "golga_GOLGA8I_vg_seq.rds"))

# ── Build per-copy transcript sequences from genome ──────────────────────────
genome <- FaFile(GENOME_FA); open(genome)
get_canonical_tx <- function(gene_id) {
  gene_txs <- unique(vg$exon_df$tx_id[vg$exon_df$gene_id == gene_id])
  cnt <- sapply(gene_txs, function(tx) sum(vg$exon_df$tx_id == tx))
  gene_txs[which.max(cnt)]
}
extract_copy_seq <- function(gene_id) {
  tx <- get_canonical_tx(gene_id)
  ex <- vg$exon_df[vg$exon_df$tx_id == tx, ]
  ex <- ex[order(ex$start), ]
  gr <- GRanges(seqnames = ex$chrom, ranges = IRanges(ex$start, ex$end),
                strand = ex$strand)
  seqs <- scanFa(genome, gr)
  if (as.character(ex$strand[1]) == "-") seqs <- rev(reverseComplement(seqs))
  list(seq = unlist(seqs))
}
copy_ids <- names(vg$paths)
copies   <- setNames(lapply(copy_ids, extract_copy_seq), copy_ids)
labels   <- setNames(sub("^gene-LOC", "", copy_ids), copy_ids)
copy_lens <- sapply(copies, function(c) length(c$seq))

# Approximate genomic START positions per copy (so we can match BAM positions)
copy_starts <- vapply(copy_ids, function(g) {
  ex <- vg$exon_df[vg$exon_df$gene_id == g, ]
  min(ex$start)
}, numeric(1))
names(copy_starts) <- copy_ids

# ── Read the BAM for primary MAPQ=0 reads in GOLGA8I region ──────────────────
golga8i_range <- GRanges("NC_073240.2", IRanges(23000000, 32000000))
param <- ScanBamParam(
  which = golga8i_range,
  what  = c("qname","flag","rname","pos","mapq","cigar","seq"),
  tag   = c("AS","NM"),
  flag  = scanBamFlag(isUnmappedQuery = FALSE)
)
cat("Reading BAM region...\n")
bam <- scanBam(BAM_FILE, param = param)[[1]]

bam_df <- data.frame(
  qname = bam$qname,
  flag  = bam$flag,
  rname = as.character(bam$rname),
  pos   = bam$pos,
  mapq  = bam$mapq,
  cigar = bam$cigar,
  AS    = bam$tag$AS,
  NM    = bam$tag$NM,
  stringsAsFactors = FALSE
)
bam_df$is_primary   <- !bitwAnd(bam_df$flag, 0x900)
bam_df$is_secondary <- bitwAnd(bam_df$flag, 0x100) != 0L

# Read sequences (only on primary records)
read_seqs <- as.character(bam$seq)
names(read_seqs) <- bam_df$qname
read_seqs <- read_seqs[bam_df$is_primary & nchar(read_seqs) > 1]

# Find primary MAPQ=0 reads (= multi-mapping reads as flagged by minimap2)
multimap_qnames <- bam_df %>%
  filter(is_primary & mapq == 0) %>%
  pull(qname)
cat("Multi-mapping (MAPQ=0) primary alignments in GOLGA8I region:", length(multimap_qnames), "\n")

# ── For each multi-mapping read, get the AS for each copy ────────────────────
# Map each alignment to the nearest annotated GOLGA8I copy
assign_copy_by_pos <- function(p) {
  diffs <- abs(p - copy_starts)
  if (min(diffs) > 50000) return(NA_character_)
  copy_ids[which.min(diffs)]
}
bam_df$copy_id <- vapply(bam_df$pos, assign_copy_by_pos, character(1))

mm_alignments <- bam_df %>%
  filter(qname %in% multimap_qnames, !is.na(copy_id)) %>%
  select(qname, copy_id, pos, mapq, AS, NM, is_primary)

cat("Sample of multi-mapping alignments:\n")
print(head(mm_alignments, 20))

# Build a wide table: row = read, columns = AS per copy (NA if not reported)
AS_mat <- mm_alignments %>%
  select(qname, copy_id, AS) %>%
  group_by(qname, copy_id) %>%
  summarise(AS = max(AS, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = copy_id, values_from = AS) %>%
  as.data.frame()

# Keep reads with at least 2 alignments AND have the actual read sequence
AS_mat$has_seq <- AS_mat$qname %in% names(read_seqs)
# Make sure all 7 copy columns exist (add NA columns for any missing)
for (cid in copy_ids) {
  if (!(cid %in% colnames(AS_mat))) AS_mat[[cid]] <- NA_real_
}
AS_mat$n_alignments <- rowSums(!is.na(AS_mat[, copy_ids]))
mm_focal <- AS_mat %>% filter(has_seq, n_alignments >= 3) %>%
  arrange(desc(n_alignments))
cat("\nFocal multi-mapping reads (>=3 GOLGA8I alignments, with sequence):",
    nrow(mm_focal), "\n")
close(genome)

# ── Compute gap-aware HMM logL for each focal read vs each copy ──────────────
SUB_MAT <- pwalign::nucleotideSubstitutionMatrix(match = 1, mismatch = -1,
                                                  baseOnly = TRUE)
P_MATCH <- 0.95; LOG_MATCH <- log(P_MATCH)
LOG_MISS <- log((1 - P_MATCH) / 3); LOG_GAP <- log(1e-6)

cat("\nScoring focal reads against each GOLGA8I transcript ...\n")
hmm_logL <- matrix(NA_real_, nrow = nrow(mm_focal), ncol = length(copy_ids))
rownames(hmm_logL) <- mm_focal$qname
colnames(hmm_logL) <- copy_ids
read_lengths <- integer(nrow(mm_focal))

for (i in seq_len(nrow(mm_focal))) {
  qn <- mm_focal$qname[i]
  read <- DNAString(read_seqs[[qn]])
  read_lengths[i] <- length(read)
  cat("  ", i, "/", nrow(mm_focal), " ", qn, " (", length(read), " bp)\n", sep = "")
  for (g in copy_ids) {
    aln <- pwalign::pairwiseAlignment(
      pattern = read, subject = copies[[g]]$seq,
      type = "local", substitutionMatrix = SUB_MAT,
      gapOpening = -2, gapExtension = -1
    )
    pat <- as.character(alignedPattern(aln))
    sub <- as.character(alignedSubject(aln))
    pc  <- strsplit(pat, "")[[1]]; sc <- strsplit(sub, "")[[1]]
    matches <- sum(pc == sc & pc != "-")
    L_read  <- length(read)
    mismatches <- L_read - matches
    # gap-aware: account for read positions that didn't align
    aln_len <- nchar(pat)
    unaligned <- L_read - aln_len
    hmm_logL[i, g] <- matches * LOG_MATCH +
                       mismatches * LOG_MISS +
                       unaligned * LOG_GAP
  }
}
mm_focal$read_len <- read_lengths

# ── Compute consistency between minimap2 AS and HMM logL ────────────────────
# For each read, get which copy each method picks
mm2_pick <- vapply(seq_len(nrow(mm_focal)), function(i) {
  as_vals <- as.numeric(mm_focal[i, copy_ids])
  if (all(is.na(as_vals))) return(NA_character_)
  copy_ids[which.max(as_vals)]
}, character(1))
hmm_pick <- vapply(seq_len(nrow(mm_focal)), function(i) {
  copy_ids[which.max(hmm_logL[i, ])]
}, character(1))
mm_focal$mm2_pick <- mm2_pick
mm_focal$hmm_pick <- hmm_pick
mm_focal$mm2_max  <- vapply(seq_len(nrow(mm_focal)), function(i)
  max(as.numeric(mm_focal[i, copy_ids]), na.rm = TRUE), numeric(1))
mm_focal$mm2_margin <- vapply(seq_len(nrow(mm_focal)), function(i) {
  v <- sort(as.numeric(mm_focal[i, copy_ids]), decreasing = TRUE, na.last = NA)
  if (length(v) < 2) NA_real_ else v[1] - v[2]
}, numeric(1))
mm_focal$hmm_margin <- apply(hmm_logL, 1, function(r) {
  sort_r <- sort(r, decreasing = TRUE)
  sort_r[1] - sort_r[2]
})
mm_focal$picks_agree <- mm_focal$mm2_pick == mm_focal$hmm_pick

cat("\n=== SUMMARY ===\n")
cat("Multi-mapping reads analyzed:", nrow(mm_focal), "\n")
cat("Reads where mm2 and HMM agree on the primary copy:",
    sum(mm_focal$picks_agree, na.rm = TRUE), "/", nrow(mm_focal), "\n\n")

# Save results
write.table(mm_focal[, c("qname","read_len","mm2_pick","mm2_margin",
                          "hmm_pick","hmm_margin","picks_agree", copy_ids)],
            file.path(FIGURES_DIR, "real_multimap_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(list(mm_focal = mm_focal, hmm_logL = hmm_logL,
              read_seqs = read_seqs[mm_focal$qname]),
        file.path(DATA_DIR, "real_multimap_golga8i.rds"))

# ── Build figure: per-read mm2 AS vs HMM logL ───────────────────────────────
score_long <- bind_rows(
  lapply(seq_len(nrow(mm_focal)), function(i) {
    qn <- mm_focal$qname[i]
    do.call(rbind, lapply(copy_ids, function(g) {
      data.frame(
        qname    = qn,
        copy     = labels[g],
        copy_id  = g,
        AS       = as.numeric(mm_focal[i, g]),
        logL     = hmm_logL[i, g],
        stringsAsFactors = FALSE
      )
    }))
  })
)
score_long <- score_long %>%
  mutate(short_qname = sub("/ccs$", "", sub("^m[0-9_]+/", "", qname)))

# Top-N reads (shortest qnames for legibility)
top_reads <- mm_focal %>% arrange(desc(read_len)) %>% head(8) %>% pull(qname)
score_top <- score_long %>% filter(qname %in% top_reads)
score_top$short_qname <- factor(score_top$short_qname,
                                 levels = sapply(top_reads, function(q)
                                   sub("/ccs$", "", sub("^m[0-9_]+/", "", q))))

# ── Crystal-clear: 4 panels (one per read), AS and HMM side-by-side per copy ──
# Each read has ONE AS per copy and ONE HMM logL per copy. Normalise so the two
# methods are visible at the same scale, but show the actual raw numbers in text.
score_top <- score_top %>%
  group_by(short_qname) %>%
  mutate(AS_max   = max(AS,   na.rm = TRUE),
         logL_max = max(logL, na.rm = TRUE),
         is_mm2_pick = AS   == AS_max,
         is_hmm_pick = logL == logL_max,
         # range-normalise both to [0,1] within each read
         AS_n   = (AS   - min(AS,   na.rm = TRUE)) /
                  pmax(1e-9, max(AS,   na.rm = TRUE) - min(AS,   na.rm = TRUE)),
         logL_n = (logL - min(logL, na.rm = TRUE)) /
                  pmax(1e-9, max(logL, na.rm = TRUE) - min(logL, na.rm = TRUE))) %>%
  ungroup()

# Build the panel-label string with read id + length
score_top <- score_top %>% left_join(
  mm_focal %>% select(qname, read_len),
  by = "qname"
)
score_top$panel_label <- sprintf("Read  %s   (%d bp)",
                                  score_top$short_qname, score_top$read_len)
panel_order <- score_top %>% distinct(panel_label, read_len) %>%
  arrange(desc(read_len)) %>% pull(panel_label)
score_top$panel_label <- factor(score_top$panel_label, levels = panel_order)

score_long <- score_top %>%
  select(panel_label, copy, AS, logL, AS_n, logL_n,
         is_mm2_pick, is_hmm_pick) %>%
  pivot_longer(cols = c(AS_n, logL_n),
                names_to = "method", values_to = "score_n") %>%
  mutate(
    method = ifelse(method == "AS_n",   "minimap2 AS", "gap-aware HMM logL"),
    method = factor(method,
                     levels = c("minimap2 AS", "gap-aware HMM logL")),
    is_pick = (method == "minimap2 AS"      & is_mm2_pick) |
              (method == "gap-aware HMM logL" & is_hmm_pick)
  )

# Per-panel pick label for annotation
pick_df <- score_top %>%
  select(panel_label, copy, is_mm2_pick, is_hmm_pick) %>%
  pivot_longer(cols = c(is_mm2_pick, is_hmm_pick),
                names_to = "method", values_to = "picked") %>%
  filter(picked) %>%
  mutate(method = ifelse(method == "is_mm2_pick",
                         "minimap2 AS", "gap-aware HMM logL"),
         method = factor(method,
                         levels = c("minimap2 AS", "gap-aware HMM logL")))

p_main <- ggplot(score_long,
                  aes(x = copy, y = score_n, fill = method,
                      alpha = ifelse(is_pick, 1.0, 0.55))) +
  geom_col(position = position_dodge(width = 0.8),
           width = 0.75, color = "white", linewidth = 0.4) +
  geom_text(data = pick_df,
            aes(x = copy, y = 1.08,
                label = ifelse(method == "minimap2 AS",
                                 "mm2 pick", "HMM pick"),
                color = method),
            inherit.aes = FALSE,
            position = position_dodge(width = 0.8),
            fontface = "bold", size = 3.0, vjust = 0) +
  scale_fill_manual(values = c("minimap2 AS"        = "#1565C0",
                                "gap-aware HMM logL" = "#1B5E20"),
                    name = NULL) +
  scale_color_manual(values = c("minimap2 AS"        = "#1565C0",
                                 "gap-aware HMM logL" = "#1B5E20"),
                     guide = "none") +
  scale_alpha_identity() +
  scale_y_continuous(limits = c(0, 1.25),
                     breaks = c(0, 0.5, 1.0)) +
  facet_wrap(~ panel_label, ncol = 4) +
  labs(title    = "minimap2 AS  vs  gap-aware HMM logL  on 4 real multi-mapping FLNC reads (GGO.bam, GOLGA8I)",
       subtitle = "Each panel = one full-length read. Each copy has ONE score per method. Scores normalised to [0,1] per read so both methods are visible on one axis.\nLabels above bars mark each method's primary pick (the copy with the highest score).",
       x = NULL, y = "score (normalised per read)") +
  theme_minimal(base_size = 10) +
  theme(plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(color = "gray40", size = 9.5,
                                      lineheight = 1.25),
        strip.text    = element_text(face = "bold", size = 10),
        strip.background = element_rect(fill = "#ECEFF1", color = NA),
        axis.text.x   = element_text(face = "bold", size = 8,
                                       angle = 30, hjust = 1),
        panel.spacing = unit(0.6, "cm"),
        legend.position = "top")

# Stub the unused old plots so the rest of the script still works
p_top <- ggplot() + theme_void()
p_bot <- ggplot() + theme_void()

# Picks table built from explicit row index
picks_df <- mm_focal[mm_focal$qname %in% top_reads, ] %>%
  mutate(short_qname = sub("/ccs$", "", sub("^m[0-9_]+/", "", qname)),
         mm2_short = labels[mm2_pick],
         hmm_short = labels[hmm_pick])
picks_df$row_y <- nrow(picks_df) - seq_len(nrow(picks_df)) + 1L
p_picks <- ggplot(picks_df) +
  annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = nrow(picks_df) + 1,
           fill = "#FAFAFA", color = "gray60", linewidth = 0.3) +
  geom_text(aes(x = 0.2, y = row_y, label = short_qname),
            hjust = 0, size = 3.0, fontface = "bold", color = "gray15") +
  geom_text(aes(x = 3.7, y = row_y, label = paste0("len ", read_len, " bp")),
            hjust = 0, size = 2.8, color = "gray45") +
  geom_text(aes(x = 5.4, y = row_y, label = paste0("mm2 -> ", mm2_short)),
            hjust = 0, size = 3.0, fontface = "bold", color = "#1565C0") +
  geom_text(aes(x = 7.4, y = row_y, label = paste0("HMM -> ", hmm_short)),
            hjust = 0, size = 3.0, fontface = "bold", color = "#1B5E20") +
  geom_text(aes(x = 9.3, y = row_y,
                label = ifelse(picks_agree, "(agree)", "(DIFFER)")),
            hjust = 0, size = 2.7,
            color = ifelse(picks_df$picks_agree, "gray45", "#B71C1C"),
            fontface = "italic") +
  annotate("text", x = 0.2, y = nrow(picks_df) + 0.5,
           label = "read", hjust = 0, fontface = "bold", size = 3.0, color = "gray35") +
  annotate("text", x = 3.7, y = nrow(picks_df) + 0.5,
           label = "length", hjust = 0, fontface = "bold", size = 3.0, color = "gray35") +
  annotate("text", x = 5.4, y = nrow(picks_df) + 0.5,
           label = "minimap2 primary", hjust = 0,
           fontface = "bold", size = 3.0, color = "gray35") +
  annotate("text", x = 7.4, y = nrow(picks_df) + 0.5,
           label = "HMM primary", hjust = 0,
           fontface = "bold", size = 3.0, color = "gray35") +
  coord_cartesian(xlim = c(0, 10), ylim = c(0.5, nrow(picks_df) + 1.2)) +
  theme_void()

# Bottom verdict
n_agree <- sum(mm_focal$picks_agree, na.rm = TRUE)
n_total <- nrow(mm_focal)
verdict_text <- sprintf(
  "On %d of %d multi-mapping FLNC reads in GGO.bam at GOLGA8I, gap-aware HMM picks a different copy than minimap2.  HMM uses the soft-clipped (unaligned) read region as evidence against copies that lack that sequence.",
  n_total - n_agree, n_total)
p_verdict <- ggplot() +
  annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = 1,
           fill = "#FFF3E0", color = "#E65100", linewidth = 0.6) +
  annotate("text", x = 0.15, y = 0.5, label = verdict_text,
           hjust = 0, vjust = 0.5, size = 3.4, color = "gray15",
           fontface = "bold", lineheight = 1.3) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 1)) +
  theme_void()

fig <- p_main / p_picks / p_verdict +
  plot_layout(heights = c(1.8, 0.5, 0.18))

ggsave(file.path(FIGURES_DIR, "fig_real_multimap_mm2_vs_hmm.pdf"), fig,
       width = 16, height = 9, units = "in")
ggsave(file.path(FIGURES_DIR, "fig_real_multimap_mm2_vs_hmm.png"), fig,
       width = 16, height = 9, units = "in", dpi = 180)
message("Real-data figure saved.")
