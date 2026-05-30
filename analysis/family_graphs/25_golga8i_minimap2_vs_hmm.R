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
TMPDIR      <- file.path(DATA_DIR, "minimap2_vs_hmm")
dir.create(TMPDIR, showWarnings = FALSE)

vg <- readRDS(file.path(DATA_DIR, "golga_GOLGA8I_vg_seq.rds"))
genome <- FaFile(GENOME_FA); open(genome)

# ── Per-copy transcripts + node-position table ────────────────────────────────
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

SOURCE   <- "gene-LOC101150678"
src_seq  <- as.character(copies[[SOURCE]]$seq)
src_chars <- strsplit(src_seq, "")[[1]]

intervals <- vg$intervals
src_pos_df <- (function() {
  ex <- copies[[SOURCE]]$exons
  ex <- ex[order(ex$start), ]
  if (copies[[SOURCE]]$strand == "-") ex <- ex[order(-ex$start), ]
  pos <- 1L; rows <- list()
  for (i in seq_len(nrow(ex))) {
    w <- ex$end[i] - ex$start[i] + 1L
    nid <- intervals$node_id[intervals$gene_id == SOURCE &
                              intervals$start == ex$start[i] &
                              intervals$end   == ex$end[i]][1]
    rows[[i]] <- data.frame(node_id = nid, tx_start = pos,
                             tx_end = pos + w - 1L, width = w)
    pos <- pos + w
  }
  do.call(rbind, rows)
})()

n16_start <- src_pos_df$tx_start[src_pos_df$node_id == 16][1]
n2_start  <- src_pos_df$tx_start[src_pos_df$node_id == 2][1]

# ── Per-copy base matrix at source coords (recompute) ────────────────────────
SUB_MAT <- pwalign::nucleotideSubstitutionMatrix(match = 1, mismatch = -1,
                                                  baseOnly = TRUE)
get_copy_bases_at_src <- function(copy_seq) {
  aln <- pwalign::pairwiseAlignment(
    pattern = DNAString(src_seq), subject = DNAString(copy_seq),
    type = "global", substitutionMatrix = SUB_MAT,
    gapOpening = -2, gapExtension = -1
  )
  ps <- strsplit(as.character(alignedPattern(aln)), "")[[1]]
  ss <- strsplit(as.character(alignedSubject(aln)), "")[[1]]
  out <- rep("-", nchar(src_seq))
  src_idx <- cumsum(ps != "-"); mask <- ps != "-"
  out[src_idx[mask]] <- ss[mask]
  out
}

cat("Building per-copy base matrix...\n")
base_mat <- matrix("-", nrow = length(copy_ids), ncol = nchar(src_seq))
rownames(base_mat) <- copy_ids
base_mat[SOURCE, ] <- src_chars
for (g in setdiff(copy_ids, SOURCE)) {
  base_mat[g, ] <- get_copy_bases_at_src(copies[[g]]$seq)
}

is_variant <- function(p) {
  bases <- base_mat[, p]
  if (any(bases == "-")) return(FALSE)
  length(unique(bases)) > 1L
}
variant_idx_global <- which(vapply(seq_len(ncol(base_mat)), is_variant, logical(1)))
cat("Clean variant positions across 7 copies:", length(variant_idx_global), "\n")

# ── Simulate reads (same as before) ──────────────────────────────────────────
set.seed(11)
ERROR_RATE <- 0.02
mutate_seq <- function(s) {
  ch  <- strsplit(s, "")[[1]]
  alt <- c("A","C","G","T")
  mp  <- sample.int(length(ch), round(length(ch) * ERROR_RATE))
  for (p in mp) ch[p] <- sample(setdiff(alt, ch[p]), 1)
  paste(ch, collapse = "")
}
sample_window <- function(tx_start, len) {
  L <- nchar(src_seq)
  start <- max(1L, min(tx_start, L - len + 1L))
  raw   <- substr(src_seq, start, start + len - 1L)
  list(read = mutate_seq(raw), tx_start = start, tx_end = start + len - 1L,
       raw_unnoised = raw)
}
sample_full <- function() {
  list(read = mutate_seq(src_seq), tx_start = 1L, tx_end = nchar(src_seq),
       raw_unnoised = src_seq)
}

reads <- list(
  R_full = list(label = "R_full  (5288 bp, full transcript)",
                 sample = sample_full()),
  R_mid  = list(label = "R_mid  (800 bp, spans bubble)",
                 sample = sample_window(max(1, n2_start - 250L), 800L)),
  R_short_variant = list(label = "R_short_variant  (200 bp, near bubble)",
                          sample = sample_window(max(1, n2_start - 80L), 200L)),
  R_short_conserved = list(label = "R_short_conserved  (200 bp, fully shared)",
                            sample = sample_window(n16_start, 200L))
)

# ── Run minimap2 (two presets) ───────────────────────────────────────────────
ref_fa <- file.path(TMPDIR, "ref.fa")
ref_set <- DNAStringSet(setNames(sapply(copy_ids, function(g) as.character(copies[[g]]$seq)),
                                  labels[copy_ids]))
writeXStringSet(ref_set, ref_fa)

reads_fa <- file.path(TMPDIR, "reads.fa")
read_vec <- sapply(reads, function(r) r$sample$read); names(read_vec) <- names(reads)
writeXStringSet(DNAStringSet(read_vec), reads_fa)

run_minimap <- function(preset, out_path) {
  args <- c("-ax", preset, "--secondary=yes", "-N", "10", "-p", "0.5",
            shQuote(ref_fa), shQuote(reads_fa))
  system2("minimap2", args = args, stdout = out_path, stderr = FALSE)
}
sam_hifi <- file.path(TMPDIR, "aln_hifi.sam"); run_minimap("map-hifi", sam_hifi)
sam_sr   <- file.path(TMPDIR, "aln_sr.sam");   run_minimap("sr",       sam_sr)

parse_sam <- function(path) {
  L <- readLines(path); L <- L[!grepl("^@", L)]
  if (length(L) == 0L) return(NULL)
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
mm2 <- bind_rows(
  parse_sam(sam_hifi) %>% filter(qname %in% c("R_full","R_mid")),
  parse_sam(sam_sr)   %>% filter(qname %in% c("R_short_variant","R_short_conserved"))
) %>% filter(!is_unmapped)

mm2$copy_label <- factor(mm2$rname, levels = unname(labels[copy_ids]))
mm2 <- mm2 %>% group_by(qname) %>%
  mutate(best_AS = max(AS), AS_gap = AS - best_AS, score_pct = AS / best_AS * 100,
         tied = score_pct >= 95) %>% ungroup()

cat("\nminimap2 alignment summary:\n")
print(mm2 %>% select(qname, rname, mapq, AS, AS_gap, score_pct, tied),
      row.names = FALSE)

# ── For each read: which variant positions fall in its source-coord range? ───
get_variants_in_read <- function(read_spec) {
  rng <- read_spec$sample$tx_start:read_spec$sample$tx_end
  variant_idx_global[variant_idx_global %in% rng]
}
read_variants <- lapply(reads, get_variants_in_read)
read_var_counts <- sapply(read_variants, length)
cat("\n# variant positions in each read's coverage:\n")
print(read_var_counts)

# Build fingerprint matrix per read at THOSE variant positions
build_fp <- function(read_name) {
  vp <- read_variants[[read_name]]
  if (length(vp) == 0L) return(NULL)
  read_chars <- strsplit(reads[[read_name]]$sample$read, "")[[1]]
  rs_start <- reads[[read_name]]$sample$tx_start
  # variant positions in read coordinates
  read_idx <- vp - rs_start + 1L
  read_bases <- read_chars[read_idx]
  m <- rbind(READ = read_bases,
              base_mat[c(SOURCE, setdiff(copy_ids, SOURCE)), vp])
  colnames(m) <- as.character(vp)
  rownames(m) <- c("READ", labels[c(SOURCE, setdiff(copy_ids, SOURCE))])
  m
}
fp_per_read <- lapply(names(reads), build_fp); names(fp_per_read) <- names(reads)

# Per-copy match count at the read's variant positions
match_at_var <- function(read_name) {
  fp <- fp_per_read[[read_name]]
  if (is.null(fp)) return(data.frame(copy_label = labels[c(SOURCE, setdiff(copy_ids, SOURCE))],
                                       matches = NA_integer_, total = 0L))
  read_row <- fp["READ", ]
  copies_only <- fp[setdiff(rownames(fp), "READ"), , drop = FALSE]
  matches <- rowSums(copies_only == matrix(read_row, nrow = nrow(copies_only),
                                              ncol = ncol(copies_only), byrow = TRUE))
  data.frame(copy_label = rownames(copies_only),
             matches    = as.integer(matches),
             total      = ncol(copies_only),
             stringsAsFactors = FALSE)
}
match_summary <- bind_rows(lapply(names(reads), function(r) {
  d <- match_at_var(r)
  d$qname <- r
  d
}))
cat("\nMatch count per copy at the read's variant positions:\n")
print(match_summary)

# ── Build one composite panel per read (minimap2 AS bar + fingerprint) ───────
nuc_pal <- c(A = "#43A047", C = "#FB8C00", G = "#E53935", T = "#1E88E5", "-" = "#BDBDBD")

mk_panel_for_read <- function(read_name) {
  rspec <- reads[[read_name]]
  mm2_r <- mm2 %>% filter(qname == read_name) %>%
    arrange(desc(AS))
  # Make sure all copies appear (sort by AS desc)
  copy_order <- as.character(mm2_r$copy_label)
  # Some copies may not have aligned — pad with AS=0
  missing <- setdiff(unname(labels[copy_ids]), copy_order)
  if (length(missing) > 0L) {
    mm2_r <- bind_rows(mm2_r, data.frame(qname = read_name, rname = missing,
                                           copy_label = factor(missing, levels = unname(labels[copy_ids])),
                                           mapq = NA, AS = 0, AS_gap = NA, NM = NA,
                                           is_unmapped = FALSE, is_secondary = NA,
                                           best_AS = max(mm2_r$AS), score_pct = 0, tied = FALSE))
  }
  mm2_r$copy_label <- factor(mm2_r$copy_label, levels = c(copy_order, missing))

  # Plot A: minimap2 AS bars (sorted by AS)
  primary_mapq <- mm2_r$mapq[!mm2_r$is_secondary & !is.na(mm2_r$mapq)][1]
  pA <- ggplot(mm2_r, aes(x = copy_label, y = AS, fill = tied)) +
    geom_col(width = 0.6, color = "white", linewidth = 0.3) +
    geom_text(aes(label = AS), vjust = -0.3, size = 2.7, color = "gray20") +
    geom_hline(yintercept = max(mm2_r$AS, na.rm = TRUE) * 0.95,
               linetype = "dashed", color = "#E65100", linewidth = 0.5) +
    scale_fill_manual(values = c(`TRUE` = "#E65100", `FALSE` = "#1565C0"),
                      guide = "none") +
    labs(title = sprintf("minimap2 AS scores  |  primary MAPQ = %s",
                          ifelse(is.na(primary_mapq), "—", primary_mapq)),
         x = NULL, y = "AS") +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 9.5),
          axis.text.x = element_text(size = 8, face = "bold", angle = 0))

  # Plot B: fingerprint matrix at the read's variant positions
  fp <- fp_per_read[[read_name]]
  if (is.null(fp)) {
    pB <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = "0 variant positions in this read.\nHMM cannot discriminate —\nuse EM with priors.",
               size = 4, color = "gray30", fontface = "italic") +
      coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
      theme_void() +
      labs(title = "HMM fingerprint at read's variant positions") +
      theme(plot.title = element_text(face = "bold", size = 9.5))
  } else {
    # Order rows: READ on top, then copies by # matches
    ms <- match_at_var(read_name) %>% arrange(desc(matches))
    row_order <- c("READ", ms$copy_label)
    long <- as.data.frame(fp) %>%
      tibble::rownames_to_column("row") %>%
      pivot_longer(cols = -row, names_to = "pos", values_to = "base")
    long$pos <- as.integer(long$pos)
    long$row <- factor(long$row, levels = rev(row_order))
    long$fill_col <- nuc_pal[long$base]
    # Mark cells matching READ
    read_base_at <- setNames(fp["READ", ], colnames(fp))
    long$matches_read <- long$base == read_base_at[as.character(long$pos)]
    long$is_read_row  <- long$row == "READ"

    pB <- ggplot(long, aes(x = factor(pos, levels = sort(unique(pos))),
                            y = row, fill = fill_col)) +
      geom_tile(color = "white", linewidth = 1) +
      scale_fill_identity() +
      geom_tile(data = filter(long, matches_read & !is_read_row),
                fill = NA, color = "#1B5E20", linewidth = 1.6) +
      geom_text(aes(label = base,
                    color = ifelse(matches_read & !is_read_row, "#1B5E20", "white")),
                size = 3.4, fontface = "bold") +
      scale_color_identity() +
      geom_hline(yintercept = length(row_order) - 0.5,
                 color = "black", linewidth = 0.7) +
      labs(title = sprintf("HMM fingerprint  |  %d variant positions covered by read",
                            ncol(fp)),
           x = "source-tx pos (bp)", y = NULL) +
      theme_minimal(base_size = 9) +
      theme(plot.title = element_text(face = "bold", size = 9.5),
            panel.grid = element_blank(),
            axis.text.x = element_text(size = 7, angle = 60, hjust = 0),
            axis.text.y = element_text(face = "bold", size = 8))
  }

  # Verdict text annotation (placed in a row)
  best_copy_mm2 <- as.character(mm2_r$copy_label[which.max(mm2_r$AS)])
  n_tied <- sum(mm2_r$tied)
  if (!is.null(fp)) {
    ms <- match_at_var(read_name) %>% arrange(desc(matches))
    best_copy_hmm <- ms$copy_label[1]
    hmm_total <- ms$total[1]
    hmm_best  <- ms$matches[1]
    hmm_second <- ms$matches[2]
    verdict <- sprintf(
      "minimap2: best=%s, MAPQ=%s, %d copies tied within 5%%   |   HMM at %d variant positions: best=%s (%d/%d), 2nd=%d   |   FINAL: %s",
      best_copy_mm2, ifelse(is.na(primary_mapq), "—", primary_mapq),
      n_tied, hmm_total, best_copy_hmm, hmm_best, hmm_total, hmm_second, best_copy_hmm
    )
  } else {
    verdict <- sprintf(
      "minimap2: best=%s, MAPQ=%s, %d copies tied within 5%%   |   HMM: 0 variant positions — cannot discriminate   |   FINAL: requires EM with priors",
      best_copy_mm2, ifelse(is.na(primary_mapq), "—", primary_mapq), n_tied
    )
  }

  p_verdict <- ggplot() +
    annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = 1,
             fill = "#F5F5F5", color = "gray60", linewidth = 0.3) +
    annotate("text", x = 0.15, y = 0.5, label = verdict,
             hjust = 0, vjust = 0.5, size = 3.0, color = "gray15") +
    coord_cartesian(xlim = c(0, 10), ylim = c(0, 1)) +
    theme_void()

  # Add a title strip for this read at the top
  p_title <- ggplot() +
    annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = 1,
             fill = "#E0E0E0") +
    annotate("text", x = 0.15, y = 0.5,
             label = rspec$label,
             hjust = 0, vjust = 0.5, size = 4.0, fontface = "bold",
             color = "gray10") +
    coord_cartesian(xlim = c(0, 10), ylim = c(0, 1)) +
    theme_void()

  list(title = p_title, A = pA, B = pB, verdict = p_verdict)
}

cat("\nBuilding 4 panels...\n")
parts <- lapply(names(reads), mk_panel_for_read)
names(parts) <- names(reads)

# Lay out with explicit design: each read = title row + (AS | FP) + verdict
# Total = 4 reads * 3 rows = 12 rows; cols = 2 (AS narrow, FP wide)
layout_design <- "
AAAAA
BBCCC
DDDDD
EEEEE
FFGGG
HHHHH
IIIII
JJKKK
LLLLL
MMMMM
NNOOO
PPPPP
"
fig <- parts$R_full$title + parts$R_full$A + parts$R_full$B + parts$R_full$verdict +
       parts$R_mid$title  + parts$R_mid$A  + parts$R_mid$B  + parts$R_mid$verdict  +
       parts$R_short_variant$title + parts$R_short_variant$A + parts$R_short_variant$B +
       parts$R_short_variant$verdict +
       parts$R_short_conserved$title + parts$R_short_conserved$A + parts$R_short_conserved$B +
       parts$R_short_conserved$verdict +
  plot_layout(design = layout_design,
              heights = c(0.18, 1, 0.18,  # R_full
                          0.18, 1, 0.18,  # R_mid
                          0.18, 1, 0.18,  # R_short_variant
                          0.18, 1, 0.18)) +
  plot_annotation(
    title    = "minimap2 alignment vs HMM fingerprint — read by read",
    subtitle = "For each read: (left) minimap2 AS scores per copy + MAPQ; (right) HMM fingerprint at variant positions in the read's coverage. Verdict combines both.",
    theme = theme(plot.title    = element_text(size = 15, face = "bold"),
                  plot.subtitle = element_text(size = 11, color = "gray40"))
  )

ggsave(file.path(FIGURES_DIR, "fig_minimap2_vs_hmm.pdf"), fig,
       width = 18, height = 16, units = "in")
ggsave(file.path(FIGURES_DIR, "fig_minimap2_vs_hmm.png"), fig,
       width = 18, height = 16, units = "in", dpi = 160)
message("Combined minimap2-vs-HMM figure saved.")
