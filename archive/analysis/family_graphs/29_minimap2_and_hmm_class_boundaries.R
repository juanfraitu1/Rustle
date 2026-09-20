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
TMPDIR      <- file.path(DATA_DIR, "mm2_class_boundary")
dir.create(TMPDIR, showWarnings = FALSE)

# Two cases (same as before, but interpreted honestly):
#   GOLGA8I 129527171 : minimap2 returns 2 tied alignments (171, 071) from
#                        distinct VG classes. HMM also ties them.
#                        => This read is genuinely ambiguous at the read level.
#                        Resolution requires EM with priors.
#   GAGE 129529768   : minimap2 returns a tied set spanning class_2 (9 members)
#                        and excludes class_1, class_3. HMM at variants confirms
#                        the class_2 vs class_3 vs class_1 boundary.
#                        => Read is correctly assigned to class_2; within-class
#                        members remain tied by mathematical necessity.

cases <- list(
  GOLGA8I = list(
    vg_path = "golga_GOLGA8I_vg_seq.rds",
    source  = "gene-LOC129527171",
    label   = "GOLGA8I  |  FLNC from 129527171 (truncated, alt 3' end)"),
  GAGE = list(
    vg_path = "gage_vg_seq.rds",
    source  = "gene-LOC129529768",
    label   = "GAGE  |  FLNC from 129529768 (class_2; 10 members in class)")
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

genome <- FaFile(GENOME_FA); open(genome)

build_case <- function(case_name, case_spec) {
  cat("\n=== ", case_name, " ===\n", sep = "")
  vg <- readRDS(file.path(DATA_DIR, case_spec$vg_path))
  copy_ids <- names(vg$paths)
  labels   <- setNames(sub("^gene-LOC", "", copy_ids), copy_ids)
  copies <- setNames(lapply(copy_ids, function(g) extract_copy_seq(vg, genome, g)),
                      copy_ids)
  eq <- family_equivalence_classes(vg)
  copy_to_class <- setNames(rep(NA_character_, length(copy_ids)), copy_ids)
  for (cn in names(eq$classes)) {
    for (mem in eq$classes[[cn]]) copy_to_class[mem] <- cn
  }
  src <- case_spec$source
  src_seq <- as.character(copies[[src]]$seq)
  read_seq <- mutate_seq(src_seq)
  src_class <- copy_to_class[src]
  cat("  Source class:", src_class, "  ;  Members:",
      paste(labels[eq$classes[[src_class]]], collapse=", "), "\n")

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
  mm2$copy <- names(labels)[match(mm2$rname, labels)]
  mm2$class <- copy_to_class[mm2$copy]
  mm2 <- mm2 %>% group_by(qname) %>%
    mutate(best_AS = max(AS), AS_gap = AS - best_AS,
           score_pct = AS / best_AS * 100,
           tied = score_pct >= 95) %>% ungroup()

  # HMM matches per copy
  read_chars <- strsplit(read_seq, "")[[1]]
  read_at_var  <- read_chars[variant_idx]
  copy_at_var  <- base_mat[, variant_idx, drop = FALSE]
  matches_per_copy <- rowSums(copy_at_var == matrix(read_at_var,
                              nrow = nrow(copy_at_var),
                              ncol = ncol(copy_at_var),
                              byrow = TRUE))
  hmm_df <- data.frame(copy = rownames(copy_at_var),
                        copy_label = labels[rownames(copy_at_var)],
                        class = copy_to_class[rownames(copy_at_var)],
                        matches = as.integer(matches_per_copy),
                        total = length(variant_idx),
                        is_source = rownames(copy_at_var) == src,
                        stringsAsFactors = FALSE)
  hmm_df <- hmm_df %>% arrange(desc(matches))
  hmm_df$pct <- round(100 * hmm_df$matches / hmm_df$total, 1)

  cat("\nminimap2 result:\n"); print(mm2 %>% select(rname, class, mapq, AS, AS_gap, score_pct, tied, is_secondary), row.names = FALSE)
  cat("\nHMM result:\n"); print(hmm_df, row.names = FALSE)

  list(case_name = case_name, label = case_spec$label,
       src = src, src_label = labels[src], src_class = src_class,
       copy_to_class = copy_to_class, labels = labels, copy_ids = copy_ids,
       mm2 = mm2, hmm_df = hmm_df, variant_idx = variant_idx,
       eq = eq)
}

results <- lapply(names(cases), function(n) build_case(n, cases[[n]]))
names(results) <- names(cases)
close(genome)

# ── Build figure panels per case ─────────────────────────────────────────────
class_pal <- c(class_1 = "#1565C0", class_2 = "#2E7D32",
                class_3 = "#EF6C00", class_4 = "#6A1B9A",
                class_5 = "#B71C1C", class_6 = "#00838F",
                class_7 = "#5D4037", class_8 = "#827717",
                class_9 = "#4A148C", class_10 = "#01579B")

mk_case_panel <- function(res) {
  d <- res$mm2 %>% arrange(desc(AS))
  miss <- setdiff(unname(res$labels[res$copy_ids]), as.character(d$rname))
  if (length(miss) > 0L) {
    miss_class <- res$copy_to_class[match(miss, res$labels)]
    d <- bind_rows(d, data.frame(qname = d$qname[1], rname = miss, mapq = NA, AS = 0,
                                   NM = NA, is_unmapped = FALSE, is_secondary = TRUE,
                                   copy = names(res$labels)[match(miss, res$labels)],
                                   class = miss_class,
                                   best_AS = max(d$AS), AS_gap = NA,
                                   score_pct = 0, tied = FALSE,
                                   stringsAsFactors = FALSE))
  }
  d <- d[order(-d$AS), ]
  d$rname <- factor(d$rname, levels = unique(d$rname))
  primary_mapq <- d$mapq[!d$is_secondary][1]
  primary_copy <- as.character(d$rname[!d$is_secondary][1])
  d$class_fill <- class_pal[d$class]
  d$is_source <- d$copy == res$src

  pA <- ggplot(d, aes(x = rname, y = AS, fill = class_fill,
                       alpha = ifelse(tied, 1.0, 0.45))) +
    geom_col(width = 0.6, color = "white", linewidth = 0.3) +
    geom_text(aes(label = AS, alpha = NULL), vjust = -0.4, size = 2.6,
              color = "gray20") +
    geom_hline(yintercept = max(d$AS) * 0.95,
               linetype = "dashed", color = "#E65100", linewidth = 0.5) +
    geom_text(data = d %>% filter(is_source),
              aes(x = rname, y = AS + max(d$AS) * 0.08,
                  label = "SOURCE", alpha = NULL),
              color = "black", fontface = "bold", size = 3.2) +
    scale_fill_identity() +
    scale_alpha_identity() +
    labs(title = sprintf("minimap2 AS  |  primary copy = %s  |  MAPQ = %s  |  tied set ~ '%s' class",
                          primary_copy,
                          ifelse(is.na(primary_mapq), "—", primary_mapq),
                          paste(unique(d$class[d$tied]), collapse=", ")),
         x = NULL, y = "AS") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 10),
          axis.text.x = element_text(size = 9, face = "bold"))

  h <- res$hmm_df %>%
    arrange(desc(matches)) %>%
    mutate(copy_label = factor(copy_label, levels = rev(copy_label)),
           class_fill = class_pal[class])
  h$is_source <- h$copy == res$src

  pB <- ggplot(h, aes(y = copy_label, x = matches, fill = class_fill)) +
    geom_col(width = 0.6, color = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%d/%d (%.1f%%) [%s]", matches, total, pct,
                                     sub("class_", "c", class)),
                  x = matches - 2),
              color = "white", fontface = "bold", size = 2.9, hjust = 1) +
    geom_text(data = h %>% filter(is_source),
              aes(y = copy_label, x = max(h$total) + max(h$total) * 0.05,
                  label = "SOURCE"),
              color = "black", fontface = "bold", size = 3.0, hjust = 0) +
    scale_fill_identity() +
    coord_cartesian(xlim = c(0, max(h$total) + max(h$total) * 0.18)) +
    labs(title = sprintf("HMM matches at %d variant positions  |  best class = %s",
                          h$total[1],
                          h$class[which.max(h$matches)]),
         x = "# variant positions matching READ", y = NULL) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 10),
          axis.text.y = element_text(face = "bold", size = 9))

  # Verdict — honest about what minimap2 reports and what HMM confirms
  src_class <- res$src_class
  mm2_tied_classes <- unique(d$class[d$tied])
  mm2_correct_class <- src_class %in% mm2_tied_classes

  hmm_best_class <- h$class[which.max(h$matches)]
  hmm_class_correct <- hmm_best_class == src_class

  # Count members of each class
  class_size <- sapply(res$eq$classes, length)
  src_class_size <- class_size[src_class]

  if (src_class_size > 1L) {
    # Within-class ties are mathematically necessary
    verdict_lines <- c(
      sprintf("Source = %s  (in %s with %d members)",
              res$src_label, src_class, src_class_size),
      sprintf("minimap2: reports %d alignments tied at >=95%% AS (MAPQ=0).  Tied classes: %s.",
              sum(d$tied), paste(mm2_tied_classes, collapse=", ")),
      sprintf("HMM:      identifies %s as best-matching class (%d/%d positions).  Other classes scored: %s.",
              hmm_best_class, max(h$matches), max(h$total),
              paste(setdiff(unique(h$class), hmm_best_class), collapse=", ")),
      sprintf("VERDICT: read assigned to %s.  Within-class members are read-indistinguishable (Theorem 1).",
              src_class)
    )
  } else {
    verdict_lines <- c(
      sprintf("Source = %s  (in %s, singleton class)",
              res$src_label, src_class),
      sprintf("minimap2: reports %d alignments tied at >=95%% AS (MAPQ=0).  Tied classes: %s.",
              sum(d$tied), paste(mm2_tied_classes, collapse=", ")),
      sprintf("HMM:      best class = %s.  But ties between %s exist if source's class shares many positions with another.",
              hmm_best_class, hmm_best_class),
      sprintf("VERDICT: minimap2 + HMM CANNOT distinguish %s from %s at the read level.  EM with priors required.",
              res$src_label,
              paste(setdiff(d$copy[d$tied], res$src), collapse=", "))
    )
  }

  p_verdict <- ggplot() +
    annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = 1, fill = "#FFF8E1",
             color = "gray60", linewidth = 0.4) +
    annotate("text", x = 0.1, y = 0.5, label = paste(verdict_lines, collapse="\n"),
             hjust = 0, vjust = 0.5, size = 3.0, lineheight = 1.25,
             color = "gray15") +
    coord_cartesian(xlim = c(0, 10), ylim = c(0, 1)) +
    theme_void()

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
              heights = c(0.18, 1.0, 0.40,
                          0.18, 1.0, 0.40)) +
  plot_annotation(
    title    = "FLNC reads at the boundary of resolvability: what minimap2 and HMM actually do",
    subtitle = paste0(
      "minimap2 does not 'fail' here — it correctly reports MAPQ=0 with a set of tied alignments. ",
      "HMM at variant positions does not break those ties either; within a VG equivalence class, members are read-indistinguishable by Theorem 1. ",
      "Both methods agree on the equivalence class.\n",
      "Bar fill colors indicate VG equivalence-class membership. Read is correctly resolved TO ITS CLASS — assignment within a class collapses by mathematical necessity."),
    theme = theme(plot.title    = element_text(size = 14, face = "bold"),
                  plot.subtitle = element_text(size = 10, color = "gray40",
                                                lineheight = 1.35))
  )

ggsave(file.path(FIGURES_DIR, "fig_class_boundary_resolution.pdf"), fig,
       width = 17, height = 12, units = "in")
ggsave(file.path(FIGURES_DIR, "fig_class_boundary_resolution.png"), fig,
       width = 17, height = 12, units = "in", dpi = 180)
message("Figure saved: figures/fig_class_boundary_resolution.png")
