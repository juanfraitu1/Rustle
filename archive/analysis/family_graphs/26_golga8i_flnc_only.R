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
TMPDIR      <- file.path(DATA_DIR, "flnc_check")
dir.create(TMPDIR, showWarnings = FALSE)

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
  list(seq = unlist(seqs), exons = ex, strand = as.character(ex$strand[1]))
}
copy_ids <- names(vg$paths)
copies   <- setNames(lapply(copy_ids, extract_copy_seq), copy_ids)
labels   <- setNames(sub("^gene-LOC", "", copy_ids), copy_ids)
close(genome)

# ── Simulate one FLNC read per copy ──────────────────────────────────────────
# FLNC = full-length non-chimeric cDNA. ~1-5 kb, full transcript with low error.
set.seed(11)
ERROR_RATE <- 0.02
mutate_seq <- function(s) {
  ch  <- strsplit(s, "")[[1]]
  alt <- c("A","C","G","T")
  mp  <- sample.int(length(ch), round(length(ch) * ERROR_RATE))
  for (p in mp) ch[p] <- sample(setdiff(alt, ch[p]), 1)
  paste(ch, collapse = "")
}

# Pick 4 representative copies (4 of 7 to keep the figure compact)
selected_sources <- c("gene-LOC101150678",
                       "gene-LOC129527057",
                       "gene-LOC129527170",
                       "gene-LOC134757232")
reads <- lapply(selected_sources, function(g) {
  list(source = g,
       label  = sprintf("FLNC from %s  (%d bp)",
                         labels[g], length(copies[[g]]$seq)),
       read   = mutate_seq(as.character(copies[[g]]$seq)),
       len    = length(copies[[g]]$seq))
})
names(reads) <- paste0("FLNC_", labels[selected_sources])

# ── Run minimap2 -ax map-hifi ────────────────────────────────────────────────
ref_fa <- file.path(TMPDIR, "ref.fa")
ref_set <- DNAStringSet(setNames(
  sapply(copy_ids, function(g) as.character(copies[[g]]$seq)),
  labels[copy_ids]))
writeXStringSet(ref_set, ref_fa)
reads_fa <- file.path(TMPDIR, "reads.fa")
read_vec <- sapply(reads, function(r) r$read); names(read_vec) <- names(reads)
writeXStringSet(DNAStringSet(read_vec), reads_fa)

sam_out <- file.path(TMPDIR, "aln.sam")
system2("minimap2",
        args = c("-ax", "map-hifi", "--secondary=yes", "-N", "10", "-p", "0.5",
                 shQuote(ref_fa), shQuote(reads_fa)),
        stdout = sam_out, stderr = FALSE)

parse_sam <- function(path) {
  L <- readLines(path); L <- L[!grepl("^@", L)]
  do.call(rbind, lapply(L, function(line) {
    f <- strsplit(line, "\t", fixed = TRUE)[[1]]
    flag <- as.integer(f[2])
    tags <- f[12:length(f)]
    get_tag <- function(p) {
      t <- tags[startsWith(tags, p)]
      if (length(t) == 0L) NA_real_ else as.numeric(sub(p, "", t[1]))
    }
    data.frame(qname = f[1], rname = f[3], mapq = as.integer(f[5]),
               AS = get_tag("AS:i:"), NM = get_tag("NM:i:"),
               is_unmapped = bitwAnd(flag, 0x4) != 0L,
               is_secondary = bitwAnd(flag, 0x100) != 0L,
               stringsAsFactors = FALSE)
  }))
}
mm2 <- parse_sam(sam_out) %>% filter(!is_unmapped)
mm2$copy_label <- factor(mm2$rname, levels = unname(labels[copy_ids]))
mm2 <- mm2 %>% group_by(qname) %>%
  mutate(best_AS = max(AS), AS_gap = AS - best_AS,
         score_pct = AS / best_AS * 100,
         tied = score_pct >= 95) %>% ungroup()

cat("\nminimap2 FLNC alignment summary:\n")
print(mm2 %>% select(qname, rname, mapq, AS, AS_gap, score_pct, tied),
      row.names = FALSE)

# Each read's source — recover from label
mm2$source <- sapply(mm2$qname, function(q) labels[reads[[q]]$source])
mm2$is_correct_primary <- (!mm2$is_secondary) & (as.character(mm2$copy_label) == mm2$source)

# ── Build the figure: one row per read, two panels (AS bar + verdict) ────────
mk_panel <- function(read_name) {
  rspec <- reads[[read_name]]
  src_label <- labels[rspec$source]
  d <- mm2 %>% filter(qname == read_name) %>% arrange(desc(AS))
  # ensure every copy appears
  miss <- setdiff(unname(labels[copy_ids]), as.character(d$copy_label))
  if (length(miss) > 0L) {
    d <- bind_rows(d, data.frame(qname = read_name,
                                   rname = miss,
                                   mapq = NA, AS = 0,
                                   NM = NA, is_unmapped = FALSE,
                                   is_secondary = TRUE, copy_label = factor(miss, levels = unname(labels[copy_ids])),
                                   best_AS = max(d$AS), AS_gap = NA, score_pct = 0, tied = FALSE,
                                   source = src_label, is_correct_primary = FALSE,
                                   stringsAsFactors = FALSE))
  }
  # Order rows by AS desc; rows with AS=0 (missing alignment) go last
  d <- d[order(-d$AS), ]
  d$copy_label <- factor(as.character(d$copy_label), levels = unique(as.character(d$copy_label)))
  primary <- d %>% filter(!is_secondary) %>% slice(1)
  primary_mapq <- primary$mapq
  primary_copy <- as.character(primary$copy_label)
  is_correct <- primary_copy == src_label

  pA <- ggplot(d, aes(x = copy_label, y = AS, fill = tied)) +
    geom_col(width = 0.6, color = "white", linewidth = 0.3) +
    geom_text(aes(label = AS), vjust = -0.3, size = 2.9, color = "gray20") +
    geom_hline(yintercept = max(d$AS) * 0.95,
               linetype = "dashed", color = "#E65100", linewidth = 0.5) +
    scale_fill_manual(values = c(`TRUE` = "#E65100", `FALSE` = "#1565C0"),
                      name = NULL,
                      labels = c(`TRUE` = ">=95% of best (tied)",
                                  `FALSE` = "below 95%")) +
    labs(title = sprintf("minimap2 AS scores  —  primary: %s, MAPQ=%s",
                           primary_copy,
                           ifelse(is.na(primary_mapq), "—", primary_mapq)),
         x = NULL, y = "AS") +
    theme_classic(base_size = 10) +
    theme(plot.title  = element_text(face = "bold", size = 10),
          axis.text.x = element_text(size = 9, face = "bold"),
          legend.position = "none")

  verdict_color <- ifelse(is_correct, "#1B5E20", "#B71C1C")
  verdict <- sprintf("Source = %s   |   minimap2 picked %s   |   MAPQ = %s   |   %s",
                      src_label, primary_copy,
                      ifelse(is.na(primary_mapq), "—", primary_mapq),
                      ifelse(is_correct, "CORRECT ✓", "WRONG ✗"))

  p_verdict <- ggplot() +
    annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = 1,
             fill = ifelse(is_correct, "#E8F5E9", "#FFEBEE"),
             color = verdict_color, linewidth = 0.4) +
    annotate("text", x = 0.15, y = 0.5, label = verdict,
             hjust = 0, vjust = 0.5, size = 3.5,
             color = verdict_color, fontface = "bold") +
    coord_cartesian(xlim = c(0, 10), ylim = c(0, 1)) +
    theme_void()

  # title for the read
  p_title <- ggplot() +
    annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = 1, fill = "#E0E0E0") +
    annotate("text", x = 0.15, y = 0.5, label = rspec$label,
             hjust = 0, vjust = 0.5, size = 4, fontface = "bold") +
    coord_cartesian(xlim = c(0, 10), ylim = c(0, 1)) +
    theme_void()

  list(title = p_title, A = pA, verdict = p_verdict)
}

parts <- lapply(names(reads), mk_panel); names(parts) <- names(reads)

# 4 read blocks; each block = title (0.15) + AS (1.0) + verdict (0.18)
layout_design <- "
AA
BB
CC
DD
EE
FF
GG
HH
II
JJ
KK
LL
"
fig <- parts[[1]]$title + parts[[1]]$A + parts[[1]]$verdict +
       parts[[2]]$title + parts[[2]]$A + parts[[2]]$verdict +
       parts[[3]]$title + parts[[3]]$A + parts[[3]]$verdict +
       parts[[4]]$title + parts[[4]]$A + parts[[4]]$verdict +
  plot_layout(design = layout_design,
              heights = c(0.15, 1.0, 0.20,
                          0.15, 1.0, 0.20,
                          0.15, 1.0, 0.20,
                          0.15, 1.0, 0.20)) +
  plot_annotation(
    title    = "FLNC read-to-copy assignment in GOLGA8I (full-length IsoSeq cDNAs only)",
    subtitle = paste0(
      "Four FLNC reads, one from each of four different GOLGA8I copies. ",
      "minimap2 -ax map-hifi against a 7-transcript FASTA. ",
      "MAPQ=60 means the aligner is confident about the source copy."),
    theme = theme(plot.title    = element_text(size = 14, face = "bold"),
                  plot.subtitle = element_text(size = 11, color = "gray40",
                                                lineheight = 1.3))
  )

ggsave(file.path(FIGURES_DIR, "fig_flnc_assignment.pdf"), fig,
       width = 14, height = 14, units = "in")
ggsave(file.path(FIGURES_DIR, "fig_flnc_assignment.png"), fig,
       width = 14, height = 14, units = "in", dpi = 180)
message("FLNC-only figure saved.")
