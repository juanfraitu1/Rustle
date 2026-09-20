suppressPackageStartupMessages({
  library(Biostrings)
  library(pwalign)
  library(Rsamtools)
  library(GenomicRanges)
  library(ggplot2)
  library(dplyr)
})

source("02_build_graphs.R")

FIGURES_DIR <- "figures"
DATA_DIR    <- "data"
GENOME_FA   <- "/mnt/c/Users/jfris/Desktop/GGO.fasta"
TMPDIR      <- file.path(DATA_DIR, "scan_ties")
dir.create(TMPDIR, showWarnings = FALSE)

# ── Families to scan: coherent multi-copy families ───────────────────────────
families <- list(
  AMY      = "amy_vg_seq.rds",
  GOLGA8I  = "golga_GOLGA8I_vg_seq.rds",
  GOLGA8NL = "golga_GOLGA8NL_vg_seq.rds",
  GOLGA6C  = "golga_GOLGA6C_vg_seq.rds",
  GAGE     = "gage_vg_seq.rds"
)

genome <- FaFile(GENOME_FA); open(genome)

get_canonical_tx <- function(vg, gene_id) {
  gene_txs <- unique(vg$exon_df$tx_id[vg$exon_df$gene_id == gene_id])
  cnt <- sapply(gene_txs, function(tx) sum(vg$exon_df$tx_id == tx))
  gene_txs[which.max(cnt)]
}
extract_copy_seq <- function(vg, gene_id) {
  tx <- get_canonical_tx(vg, gene_id)
  ex <- vg$exon_df[vg$exon_df$tx_id == tx, ]; ex <- ex[order(ex$start), ]
  gr <- GRanges(seqnames = ex$chrom, ranges = IRanges(ex$start, ex$end),
                strand = ex$strand)
  seqs <- scanFa(genome, gr)
  if (as.character(ex$strand[1]) == "-") seqs <- rev(reverseComplement(seqs))
  list(seq = unlist(seqs))
}

set.seed(11)
ERROR_RATE <- 0.02
mutate_seq <- function(s) {
  ch  <- strsplit(s, "")[[1]]
  alt <- c("A","C","G","T")
  mp  <- sample.int(length(ch), round(length(ch) * ERROR_RATE))
  for (p in mp) ch[p] <- sample(setdiff(alt, ch[p]), 1)
  paste(ch, collapse = "")
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

# ── Per-family: simulate one FLNC from each copy + minimap2 align ────────────
all_results <- list()
for (fam_name in names(families)) {
  cat("\n=== ", fam_name, " ===\n", sep = "")
  vg <- readRDS(file.path(DATA_DIR, families[[fam_name]]))
  copy_ids <- names(vg$paths)
  labels   <- setNames(sub("^gene-LOC", "", copy_ids), copy_ids)

  copies <- setNames(lapply(copy_ids, function(g) extract_copy_seq(vg, g)),
                      copy_ids)
  tx_lens <- sapply(copies, function(c) length(c$seq))
  cat("Transcript lengths: ", paste(tx_lens, collapse = ", "), "\n")

  # Write reference + simulated FLNC reads
  fam_dir <- file.path(TMPDIR, fam_name); dir.create(fam_dir, showWarnings = FALSE)
  ref_fa  <- file.path(fam_dir, "ref.fa")
  reads_fa <- file.path(fam_dir, "reads.fa")
  sam_out <- file.path(fam_dir, "aln.sam")

  ref_set <- DNAStringSet(setNames(
    sapply(copy_ids, function(g) as.character(copies[[g]]$seq)),
    labels[copy_ids]
  ))
  writeXStringSet(ref_set, ref_fa)

  read_seqs <- sapply(copy_ids, function(g) mutate_seq(as.character(copies[[g]]$seq)))
  read_names <- paste0("FLNC_", labels[copy_ids])
  reads_set <- DNAStringSet(setNames(read_seqs, read_names))
  writeXStringSet(reads_set, reads_fa)

  system2("minimap2",
          args = c("-ax", "map-hifi", "--secondary=yes",
                   "-N", as.character(length(copy_ids) + 2L),
                   "-p", "0.5",
                   shQuote(ref_fa), shQuote(reads_fa)),
          stdout = sam_out, stderr = FALSE)

  rec <- parse_sam(sam_out)
  if (is.null(rec) || nrow(rec) == 0L) { next }
  rec <- rec[!rec$is_unmapped, ]
  rec$source <- sapply(rec$qname, function(q) sub("^FLNC_", "", q))
  rec <- rec %>% group_by(qname) %>%
    mutate(best_AS = max(AS, na.rm = TRUE),
           AS_gap = AS - best_AS,
           score_pct = AS / best_AS * 100,
           tied = score_pct >= 95) %>% ungroup()

  # Per-read summary
  per_read <- rec %>% group_by(qname) %>%
    summarise(family = fam_name,
              source = source[1],
              primary_copy = rname[!is_secondary][1],
              primary_AS   = max(AS),
              primary_MAPQ = mapq[!is_secondary][1],
              n_tied       = sum(tied),
              second_best_AS  = sort(AS, decreasing = TRUE)[2],
              gap_to_2nd      = max(AS) - sort(AS, decreasing = TRUE)[2],
              gap_pct_to_2nd  = (max(AS) - sort(AS, decreasing = TRUE)[2]) / max(AS) * 100,
              correct = source == primary_copy,
              .groups = "drop")
  print(per_read, n = Inf)
  all_results[[fam_name]] <- per_read
}
close(genome)

# ── Combine across families and rank by "tied-ness" ──────────────────────────
big <- bind_rows(all_results)
cat("\n\n========== ALL FAMILIES — minimap2 MAPQ per FLNC read ==========\n\n")
print(big, n = Inf)

cat("\nReads with MAPQ < 60  (i.e. minimap2 NOT confident):\n")
ambiguous <- big %>% filter(primary_MAPQ < 60)
print(ambiguous, n = Inf)

cat("\nReads where minimap2 picked the WRONG copy:\n")
wrong <- big %>% filter(!correct)
print(wrong, n = Inf)

write.table(big, file.path(FIGURES_DIR, "minimap2_flnc_ties_scan.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("\nScan results written to figures/minimap2_flnc_ties_scan.tsv")
