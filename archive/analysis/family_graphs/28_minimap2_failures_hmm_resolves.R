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
TMPDIR      <- file.path(DATA_DIR, "mm2_fail_hmm")
dir.create(TMPDIR, showWarnings = FALSE)

# ── Define the failure cases we found ────────────────────────────────────────
cases <- list(
  GOLGA8I = list(
    vg_path = "golga_GOLGA8I_vg_seq.rds",
    source  = "gene-LOC129527171",
    label   = "GOLGA8I  —  FLNC from 129527171 (truncated copy, alt 3' end)"),
  GAGE = list(
    vg_path = "gage_vg_seq.rds",
    source  = "gene-LOC129529768",
    label   = "GAGE  —  FLNC from 129529768 (class_2)")
)

set.seed(11)
ERROR_RATE <- 0.02
mutate_seq <- function(s) {
  ch  <- strsplit(s, "")[[1]]
  alt <- c("A","C","G","T")
  mp  <- sample.int(length(ch), round(length(ch) * ERROR_RATE))
  for (p in mp) ch[p] <- sample(setdiff(alt, ch[p]), 1)
  paste(ch, collapse = "")
}

get_canonical_tx <- function(vg, gene_id) {
  gene_txs <- unique(vg$exon_df$tx_id[vg$exon_df$gene_id == gene_id])
  cnt <- sapply(gene_txs, function(tx) sum(vg$exon_df$tx_id == tx))
  gene_txs[which.max(cnt)]
}
extract_copy_seq <- function(vg, genome, gene_id) {
  tx <- get_canonical_tx(vg, gene_id)
  ex <- vg$exon_df[vg$exon_df$tx_id == tx, ]; ex <- ex[order(ex$start), ]
  gr <- GRanges(seqnames = ex$chrom, ranges = IRanges(ex$start, ex$end),
                strand = ex$strand)
  seqs <- scanFa(genome, gr)
  if (as.character(ex$strand[1]) == "-") seqs <- rev(reverseComplement(seqs))
  list(seq = unlist(seqs))
}

SUB_MAT <- pwalign::nucleotideSubstitutionMatrix(match = 1, mismatch = -1,
                                                  baseOnly = TRUE)
get_copy_bases_at_src <- function(src_seq, copy_seq) {
  aln <- pwalign::pairwiseAlignment(
    pattern = DNAString(src_seq), subject = DNAString(copy_seq),
    type = "global", substitutionMatrix = SUB_MAT,
    gapOpening = -2, gapExtension = -1
  )
  ps <- strsplit(as.character(alignedPattern(aln)), "")[[1]]
  ss <- strsplit(as.character(alignedSubject(aln)), "")[[1]]
  out <- rep("-", nchar(as.character(src_seq)))
  src_idx <- cumsum(ps != "-"); mask <- ps != "-"
  out[src_idx[mask]] <- ss[mask]
  out
}

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

# ── For each case, run analysis ──────────────────────────────────────────────
nuc_pal <- c(A = "#43A047", C = "#FB8C00", G = "#E53935", T = "#1E88E5", "-" = "#BDBDBD")

genome <- FaFile(GENOME_FA); open(genome)

build_case <- function(case_name, case_spec) {
  cat("\n=== ", case_name, " ===\n", sep = "")
  vg <- readRDS(file.path(DATA_DIR, case_spec$vg_path))
  copy_ids <- names(vg$paths)
  labels   <- setNames(sub("^gene-LOC", "", copy_ids), copy_ids)

  copies <- setNames(lapply(copy_ids, function(g) extract_copy_seq(vg, genome, g)),
                      copy_ids)
  src <- case_spec$source
  src_seq <- as.character(copies[[src]]$seq)
  read_seq <- mutate_seq(src_seq)

  # Per-copy base matrix at source positions
  cat("  Computing per-copy base matrix at source coordinates...\n")
  base_mat <- matrix("-", nrow = length(copy_ids), ncol = nchar(src_seq))
  rownames(base_mat) <- copy_ids
  base_mat[src, ] <- strsplit(src_seq, "")[[1]]
  for (g in setdiff(copy_ids, src)) {
    base_mat[g, ] <- get_copy_bases_at_src(src_seq, copies[[g]]$seq)
  }
  variant_idx <- which(vapply(seq_len(ncol(base_mat)), function(p) {
    bases <- base_mat[, p]; if (any(bases == "-")) return(FALSE)
    length(unique(bases)) > 1L
  }, logical(1)))
  cat("  Variant positions across all copies:", length(variant_idx), "\n")

  # Write reference and read; run minimap2
  cdir <- file.path(TMPDIR, case_name); dir.create(cdir, showWarnings = FALSE)
  ref_fa <- file.path(cdir, "ref.fa")
  reads_fa <- file.path(cdir, "reads.fa")
  sam_out <- file.path(cdir, "aln.sam")
  ref_set <- DNAStringSet(setNames(
    sapply(copy_ids, function(g) as.character(copies[[g]]$seq)),
    labels[copy_ids]))
  writeXStringSet(ref_set, ref_fa)
  read_name <- paste0("FLNC_", labels[src])
  writeXStringSet(DNAStringSet(setNames(read_seq, read_name)), reads_fa)
  system2("minimap2",
          args = c("-ax", "map-hifi", "--secondary=yes",
                   "-N", as.character(length(copy_ids) + 2L),
                   "-p", "0.5",
                   shQuote(ref_fa), shQuote(reads_fa)),
          stdout = sam_out, stderr = FALSE)
  mm2 <- parse_sam(sam_out) %>% filter(!is_unmapped)
  mm2$copy_label <- factor(mm2$rname, levels = unname(labels[copy_ids]))
  mm2 <- mm2 %>% group_by(qname) %>%
    mutate(best_AS = max(AS), AS_gap = AS - best_AS,
           score_pct = AS / best_AS * 100,
           tied = score_pct >= 95) %>% ungroup()
  cat("\nminimap2 result:\n")
  print(mm2 %>% select(rname, mapq, AS, AS_gap, score_pct, tied, is_secondary),
        row.names = FALSE)

  # For each copy, count matches with the READ at variant positions
  read_chars <- strsplit(read_seq, "")[[1]]
  read_at_var  <- read_chars[variant_idx]
  copy_at_var  <- base_mat[, variant_idx, drop = FALSE]
  matches_per_copy <- rowSums(copy_at_var == matrix(read_at_var,
                              nrow = nrow(copy_at_var),
                              ncol = ncol(copy_at_var),
                              byrow = TRUE))
  hmm_df <- data.frame(copy = rownames(copy_at_var),
                        copy_label = labels[rownames(copy_at_var)],
                        matches = as.integer(matches_per_copy),
                        total   = length(variant_idx),
                        is_source = rownames(copy_at_var) == src,
                        stringsAsFactors = FALSE)
  hmm_df <- hmm_df %>% arrange(desc(matches))
  hmm_df$pct <- round(100 * hmm_df$matches / hmm_df$total, 1)
  cat("\nHMM match counts at variant positions:\n")
  print(hmm_df, row.names = FALSE)

  # Verdict comparison
  mm2_primary <- mm2 %>% filter(!is_secondary) %>% pull(copy_label) %>% as.character()
  hmm_primary <- hmm_df$copy_label[1]
  list(case_name = case_name, label = case_spec$label,
       src_label = labels[src], src = src,
       mm2 = mm2, hmm_df = hmm_df, variant_idx = variant_idx,
       read_chars = read_chars, base_mat = base_mat,
       copy_ids = copy_ids, labels = labels,
       mm2_primary = mm2_primary, hmm_primary = hmm_primary)
}

results <- lapply(names(cases), function(n) build_case(n, cases[[n]]))
names(results) <- names(cases)
close(genome)

# ── Build one composite panel per case ───────────────────────────────────────
mk_case_panel <- function(res) {
  # left: minimap2 AS bars
  d <- res$mm2 %>% arrange(desc(AS))
  miss <- setdiff(unname(res$labels[res$copy_ids]),
                  as.character(d$copy_label))
  if (length(miss) > 0L) {
    d_miss <- data.frame(qname = d$qname[1], rname = miss,
                         mapq = NA, AS = 0, NM = NA,
                         is_unmapped = FALSE, is_secondary = TRUE,
                         copy_label = factor(miss, levels = unname(res$labels[res$copy_ids])),
                         best_AS = max(d$AS), AS_gap = NA,
                         score_pct = 0, tied = FALSE,
                         stringsAsFactors = FALSE)
    d <- bind_rows(d, d_miss)
  }
  d <- d[order(-d$AS), ]
  d$copy_label <- factor(as.character(d$copy_label),
                         levels = unique(as.character(d$copy_label)))
  primary_mapq <- d$mapq[!d$is_secondary][1]
  d$is_source <- as.character(d$copy_label) == res$src_label

  pA <- ggplot(d, aes(x = copy_label, y = AS, fill = tied)) +
    geom_col(width = 0.6, color = "white", linewidth = 0.3) +
    geom_text(aes(label = AS), vjust = -0.4, size = 2.6, color = "gray20") +
    geom_hline(yintercept = max(d$AS) * 0.95,
               linetype = "dashed", color = "#E65100", linewidth = 0.5) +
    geom_text(data = d %>% filter(is_source),
              aes(x = copy_label, y = AS + max(d$AS) * 0.08,
                  label = "SOURCE"),
              color = "#1B5E20", fontface = "bold", size = 3.2) +
    scale_fill_manual(values = c(`TRUE` = "#E65100", `FALSE` = "#1565C0"),
                      guide = "none") +
    labs(title = sprintf("minimap2 AS  |  primary: %s  |  MAPQ = %s",
                          res$mm2_primary,
                          ifelse(is.na(primary_mapq), "—", primary_mapq)),
         x = NULL, y = "AS") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 10),
          axis.text.x = element_text(size = 9, face = "bold"))

  # right: HMM match counts per copy at variant positions
  h <- res$hmm_df %>%
    mutate(copy_label = factor(copy_label, levels = rev(copy_label)),
           bar_color = ifelse(matches == max(matches), "#1B5E20", "#9E9E9E"))
  h$is_source <- h$copy == res$src

  pB <- ggplot(h, aes(y = copy_label, x = matches,
                       fill = bar_color)) +
    geom_col(width = 0.6, color = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%d / %d (%.1f%%)", matches, total, pct),
                  x = matches - 2),
              color = "white", fontface = "bold", size = 3, hjust = 1) +
    geom_text(data = h %>% filter(is_source),
              aes(y = copy_label, x = total + 4,
                  label = "<- SOURCE"),
              color = "#1B5E20", fontface = "bold", size = 3.0, hjust = 0) +
    scale_fill_identity() +
    coord_cartesian(xlim = c(0, max(h$total) + 25)) +
    labs(title = sprintf("HMM matches at %d variant positions  |  best: %s",
                          h$total[1], res$hmm_primary),
         x = "# variant positions matching READ", y = NULL) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 10),
          axis.text.y = element_text(face = "bold", size = 9))

  # Verdict
  agreed <- res$mm2_primary == res$hmm_primary
  verdict_color <- ifelse(res$hmm_primary == res$src_label, "#1B5E20", "#B71C1C")
  verdict <- sprintf(
    "Source = %s   |   minimap2 = %s (%s)   |   HMM = %s (%s)",
    res$src_label,
    res$mm2_primary, ifelse(res$mm2_primary == res$src_label, "OK", "WRONG"),
    res$hmm_primary, ifelse(res$hmm_primary == res$src_label, "RECOVERS SOURCE",
                              "still wrong"))
  p_verdict <- ggplot() +
    annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = 1,
             fill = "#FFF8E1", color = "gray60", linewidth = 0.4) +
    annotate("text", x = 0.15, y = 0.5, label = verdict,
             hjust = 0, vjust = 0.5, size = 3.6,
             color = verdict_color, fontface = "bold") +
    coord_cartesian(xlim = c(0, 10), ylim = c(0, 1)) +
    theme_void()

  # title
  p_title <- ggplot() +
    annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = 1, fill = "#E0E0E0") +
    annotate("text", x = 0.15, y = 0.5, label = res$label,
             hjust = 0, vjust = 0.5, size = 4.0, fontface = "bold") +
    coord_cartesian(xlim = c(0, 10), ylim = c(0, 1)) +
    theme_void()

  list(title = p_title, A = pA, B = pB, verdict = p_verdict)
}

parts <- lapply(results, mk_case_panel)

layout_design <- "
AAAA
BBCC
DDDD
EEEE
FFGG
HHHH
"
fig <- parts[[1]]$title + parts[[1]]$A + parts[[1]]$B + parts[[1]]$verdict +
       parts[[2]]$title + parts[[2]]$A + parts[[2]]$B + parts[[2]]$verdict +
  plot_layout(design = layout_design,
              heights = c(0.18, 1.0, 0.22,
                          0.18, 1.0, 0.22)) +
  plot_annotation(
    title    = "FLNC reads where minimap2 fails — and how HMM at variant positions resolves",
    subtitle = paste0(
      "Two real cases found by scanning 5 families. ",
      "Source = the copy the FLNC was simulated from. ",
      "Left = minimap2 picks the wrong copy with low MAPQ; ",
      "right = HMM match counts at variant positions recover the source."),
    theme = theme(plot.title    = element_text(size = 14, face = "bold"),
                  plot.subtitle = element_text(size = 10.5, color = "gray40",
                                                lineheight = 1.3))
  )

ggsave(file.path(FIGURES_DIR, "fig_minimap2_fails_hmm_rescues.pdf"), fig,
       width = 16, height = 11, units = "in")
ggsave(file.path(FIGURES_DIR, "fig_minimap2_fails_hmm_rescues.png"), fig,
       width = 16, height = 11, units = "in", dpi = 180)
message("Figure saved.")
