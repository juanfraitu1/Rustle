suppressPackageStartupMessages({
  library(Biostrings)
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
TMPDIR      <- file.path(DATA_DIR, "minimap2_tie_check")
dir.create(TMPDIR, showWarnings = FALSE)

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

# ── Build the per-copy node-position table to anchor reads ───────────────────
intervals <- vg$intervals
copy_node_pos <- list()
for (gid in copy_ids) {
  ex <- copies[[gid]]$exons
  ex <- ex[order(ex$start), ]
  if (copies[[gid]]$strand == "-") ex <- ex[order(-ex$start), ]
  pos <- 1L; rows <- list()
  for (i in seq_len(nrow(ex))) {
    w   <- ex$end[i] - ex$start[i] + 1L
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

# ── Write the 7 GOLGA8I transcripts as a reference FASTA ─────────────────────
ref_fa <- file.path(TMPDIR, "golga8i_transcripts.fa")
ref_set <- DNAStringSet(setNames(
  lapply(copy_ids, function(g) copies[[g]]$seq),
  labels[copy_ids]
))
writeXStringSet(ref_set, ref_fa)
cat("Reference written: ", ref_fa, " (", length(ref_set), " sequences)\n", sep = "")

# ── Build a small set of simulated reads ─────────────────────────────────────
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
  mutate_seq(raw)
}

SOURCE   <- "gene-LOC101150678"
src_pos  <- copy_node_pos[[SOURCE]]
n16_start <- src_pos$tx_start[src_pos$node_id == 16][1]
n2_start  <- src_pos$tx_start[src_pos$node_id == 2][1]

reads <- list(
  R_full = list(
    label = "R_full  (full transcript, 5288 bp)",
    seq   = mutate_seq(as.character(copies[[SOURCE]]$seq))),
  R_mid = list(
    label = "R_mid   (800 bp, spans bubble of 678)",
    seq   = sample_window(SOURCE, max(1, n2_start - 250L), 800L)),
  R_short_conserved = list(
    label = "R_short_conserved  (200 bp, fully-shared node 16)",
    seq   = sample_window(SOURCE, n16_start, 200L)),
  R_short_variant = list(
    label = "R_short_variant  (200 bp, covers some variant SNVs)",
    seq   = sample_window(SOURCE, max(1, n2_start - 80L), 200L))
)

reads_fa <- file.path(TMPDIR, "reads.fa")
read_vec <- sapply(reads, function(r) r$seq)
names(read_vec) <- names(reads)
reads_set <- DNAStringSet(read_vec)
writeXStringSet(reads_set, reads_fa)
cat("Reads written: ", reads_fa, " (", length(reads_set), " reads)\n", sep = "")

# ── Run minimap2 — two passes so both long & short reads align ───────────────
# Pass 1: map-hifi for full-length / long reads
# Pass 2: sr (short-read) for the 200 bp reads
sam_hifi <- file.path(TMPDIR, "aln_hifi.sam")
sam_sr   <- file.path(TMPDIR, "aln_sr.sam")
common_args <- c("--secondary=yes", "-N", "10", "-p", "0.5",
                 shQuote(ref_fa), shQuote(reads_fa))
ret1 <- system2("minimap2",
                args = c("-ax", "map-hifi", common_args),
                stdout = sam_hifi, stderr = FALSE)
ret2 <- system2("minimap2",
                args = c("-ax", "sr", common_args),
                stdout = sam_sr, stderr = FALSE)
cat("minimap2 map-hifi exit: ", ret1, "\n", sep = "")
cat("minimap2 sr      exit: ", ret2, "\n", sep = "")

# ── Parse SAM (both runs) ────────────────────────────────────────────────────
parse_sam_file <- function(path, preset) {
  lines <- readLines(path)
  lines <- lines[!grepl("^@", lines)]
  if (length(lines) == 0L) return(NULL)
  out <- do.call(rbind, lapply(lines, parse_record))
  out$preset <- preset
  out
}
# Keep the long-read assignments from map-hifi for R_full/R_mid;
# use sr only for R_short_* reads.
sam_lines <- character(0)  # placeholder so old-flow parse stays

parse_record <- function(line) {
  f <- strsplit(line, "\t", fixed = TRUE)[[1]]
  qname <- f[1]; flag <- as.integer(f[2]); rname <- f[3]; pos <- as.integer(f[4])
  mapq  <- as.integer(f[5])
  cigar <- f[6]
  # Optional tags
  tags <- f[12:length(f)]
  get_tag <- function(prefix) {
    tag <- tags[startsWith(tags, prefix)]
    if (length(tag) == 0L) return(NA)
    val <- sub(paste0("^", prefix), "", tag[1])
    as.numeric(val)
  }
  data.frame(
    qname = qname, flag = flag, rname = rname, pos = pos, mapq = mapq,
    cigar = cigar,
    AS = get_tag("AS:i:"),
    NM = get_tag("NM:i:"),
    is_unmapped  = bitwAnd(flag, 0x4) != 0L,
    is_secondary = bitwAnd(flag, 0x100) != 0L,
    is_supplementary = bitwAnd(flag, 0x800) != 0L,
    stringsAsFactors = FALSE
  )
}

rec_hifi <- parse_sam_file(sam_hifi, "map-hifi")
rec_sr   <- parse_sam_file(sam_sr,   "sr")

# Merge: prefer map-hifi for long reads; use sr for the 200 bp reads
long_reads  <- c("R_full","R_mid")
short_reads <- c("R_short_conserved","R_short_variant")
records <- bind_rows(
  if (!is.null(rec_hifi)) rec_hifi %>% filter(qname %in% long_reads) else NULL,
  if (!is.null(rec_sr))   rec_sr   %>% filter(qname %in% short_reads) else NULL
)
records <- records[!records$is_unmapped, ]

cat("\nminimap2 alignment records:\n")
print(records[, c("qname","rname","mapq","AS","NM","is_secondary","is_supplementary")],
      row.names = FALSE)

# Order records by read and AS score
records$read_label <- sapply(records$qname, function(q) reads[[q]]$label)
records$alignment_kind <- ifelse(records$is_secondary, "secondary",
                                  ifelse(records$is_supplementary, "supplementary",
                                         "primary"))

# Best AS per read
best_AS <- records %>% group_by(qname) %>% summarise(best = max(AS, na.rm = TRUE))
records <- records %>% left_join(best_AS, by = "qname") %>%
  mutate(AS_gap = AS - best, score_pct = AS / best * 100,
         tied   = score_pct >= 95)

cat("\nGap to best per record:\n")
print(records[, c("read_label","rname","mapq","AS","AS_gap","score_pct","alignment_kind","tied")],
      row.names = FALSE)

# Save TSV
write.table(records[, c("qname","read_label","rname","mapq","AS","NM","AS_gap","score_pct",
                          "alignment_kind","tied")],
            file.path(FIGURES_DIR, "golga8i_minimap2_alignments.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ── Visualise: AS score per copy per read ────────────────────────────────────
records$read_label <- factor(records$read_label,
                              levels = sapply(reads, function(r) r$label))
records$rname <- factor(records$rname, levels = unname(labels[copy_ids]))

p <- ggplot(records, aes(x = rname, y = AS,
                          fill = tied,
                          alpha = ifelse(alignment_kind == "primary", 1, 0.6))) +
  geom_col(width = 0.55, color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("AS=%d\nMAPQ=%d", AS, mapq)),
            vjust = -0.3, size = 2.5, color = "gray20") +
  geom_hline(data = records %>% group_by(read_label) %>%
                       summarise(best = max(AS), .groups = "drop") %>%
                       mutate(thresh = best * 0.95),
             aes(yintercept = thresh),
             linetype = "dashed", color = "#E65100", linewidth = 0.5) +
  scale_fill_manual(values = c(`TRUE` = "#E65100", `FALSE` = "#1565C0"),
                     name = NULL,
                     labels = c(`TRUE` = "tied (AS >= 95% of best)",
                                 `FALSE` = "clearly worse")) +
  scale_alpha_identity() +
  facet_wrap(~ read_label, scales = "free_y", ncol = 1) +
  labs(title    = "Real minimap2 alignment scores against the 7 GOLGA8I transcripts",
       subtitle = paste0(
         "Reference = the 7 transcripts written as a 7-sequence FASTA. ",
         "minimap2 -ax map-hifi --secondary=yes -N 10 -p 0.5.\n",
         "AS = SAM alignment score; MAPQ = mapping quality (0 = aligner cannot decide). ",
         "Dashed line = 95% of best AS."),
       x = NULL, y = "minimap2 AS score") +
  theme_minimal(base_size = 11) +
  theme(plot.title       = element_text(face = "bold", size = 13),
        plot.subtitle    = element_text(color = "gray40", size = 10,
                                         lineheight = 1.3),
        strip.text       = element_text(face = "bold", size = 10),
        axis.text.x      = element_text(size = 9, face = "bold"),
        legend.position  = "top")

ggsave(file.path(FIGURES_DIR, "fig_golga8i_minimap2_ties.pdf"), p,
       width = 13, height = 11, units = "in")
ggsave(file.path(FIGURES_DIR, "fig_golga8i_minimap2_ties.png"), p,
       width = 13, height = 11, units = "in", dpi = 180)
message("\nFigure saved: figures/fig_golga8i_minimap2_ties.png")
message("Records TSV saved: figures/golga8i_minimap2_alignments.tsv")
