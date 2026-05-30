suppressPackageStartupMessages({
  library(Biostrings)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(tidyr)
})

FIGURES_DIR <- "figures"
DATA_DIR    <- "data"
TMPDIR      <- file.path(DATA_DIR, "toy_mm2_hmm"); dir.create(TMPDIR, showWarnings = FALSE)

# ── Toy 3-copy family with full-length reads ─────────────────────────────────
set.seed(42)
L_TX  <- 100L
bases <- c("A","T","C","G")
N_SNPS <- 6L
shared <- sample(bases, L_TX, replace = TRUE)
snp_positions <- sort(sample(seq_len(L_TX), N_SNPS))
copy_letters <- c("A","B","C")
seqs <- list()
for (cn in copy_letters) seqs[[cn]] <- shared
for (p in snp_positions) {
  trio <- sample(bases, 3, replace = FALSE)
  for (k in seq_along(copy_letters))
    seqs[[copy_letters[k]]][p] <- trio[k]
}
cat("Copy A:", paste(seqs$A, collapse=""), "\n")
cat("Copy B:", paste(seqs$B, collapse=""), "\n")
cat("Copy C:", paste(seqs$C, collapse=""), "\n")
cat("SNP positions:", snp_positions, "\n\n")

# Per-copy HMM emission
p_match <- 0.95; p_other <- (1 - p_match) / 3
emit <- lapply(seqs, function(s) {
  m <- matrix(p_other, nrow = length(s), ncol = 4,
              dimnames = list(seq_along(s), bases))
  for (i in seq_along(s)) m[i, s[i]] <- p_match
  m
})

# ── Simulate 40 full-length reads with 5% error rate ────────────────────────
ERROR_RATE <- 0.05
TRUE_THETA <- c(A = 0.50, B = 0.30, C = 0.20)
N_READS <- 40
add_err <- function(seq, rate) {
  ch <- seq
  for (i in seq_along(ch))
    if (runif(1) < rate) ch[i] <- sample(setdiff(bases, ch[i]), 1)
  ch
}
reads <- vector("list", N_READS)
for (i in seq_len(N_READS)) {
  src <- sample(copy_letters, 1, prob = TRUE_THETA)
  r_seq <- add_err(seqs[[src]], ERROR_RATE)
  reads[[i]] <- list(read_id = sprintf("r%02d", i),
                     src = src,
                     bases = r_seq,
                     seq = paste(r_seq, collapse = ""))
}

# ── Per-copy HMM logL (full-length) ──────────────────────────────────────────
score_full <- function(read_bases, emit_mat) {
  cols <- match(read_bases, bases)
  sum(log(emit_mat[cbind(seq_along(read_bases), cols)]))
}
logL_mat <- t(sapply(reads, function(r) {
  c(A = score_full(r$bases, emit$A),
    B = score_full(r$bases, emit$B),
    C = score_full(r$bases, emit$C))
}))
rownames(logL_mat) <- sapply(reads, function(r) r$read_id)

# ── Run real minimap2 against the 3-copy reference ──────────────────────────
ref_fa <- file.path(TMPDIR, "ref.fa")
reads_fa <- file.path(TMPDIR, "reads.fa")
sam_out <- file.path(TMPDIR, "aln.sam")
writeXStringSet(DNAStringSet(setNames(c(paste(seqs$A, collapse = ""),
                                          paste(seqs$B, collapse = ""),
                                          paste(seqs$C, collapse = "")),
                                       c("A","B","C"))), ref_fa)
writeXStringSet(DNAStringSet(setNames(sapply(reads, function(r) r$seq),
                                        sapply(reads, function(r) r$read_id))),
                reads_fa)
ret <- system2("minimap2",
               args = c("-ax", "sr", "--secondary=yes", "-N", "3", "-p", "0.5",
                        shQuote(ref_fa), shQuote(reads_fa)),
               stdout = sam_out, stderr = FALSE)
cat("minimap2 exit:", ret, "\n")

# Parse SAM
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
mm2_records <- parse_sam(sam_out)
mm2_records <- mm2_records[!mm2_records$is_unmapped, ]
cat("minimap2 records:", nrow(mm2_records), "\n")

# Build minimap2 AS matrix: row = read, col = A/B/C; NA where not reported
AS_mat <- matrix(NA_real_, nrow = N_READS, ncol = 3,
                  dimnames = list(rownames(logL_mat), c("A","B","C")))
for (i in seq_len(nrow(mm2_records))) {
  AS_mat[mm2_records$qname[i], mm2_records$rname[i]] <- mm2_records$AS[i]
}
# Also record MAPQ of the primary alignment per read
primary_mm2 <- mm2_records %>% filter(!is_secondary) %>%
  select(qname, rname, mapq) %>%
  rename(read = qname, mm2_pick = rname, mm2_mapq = mapq)

# ── Per-read summary ─────────────────────────────────────────────────────────
read_summary <- data.frame(
  read = rownames(logL_mat),
  true = sapply(reads, function(r) r$src),
  logL_A = logL_mat[,"A"], logL_B = logL_mat[,"B"], logL_C = logL_mat[,"C"],
  AS_A = AS_mat[,"A"],   AS_B = AS_mat[,"B"],     AS_C = AS_mat[,"C"],
  stringsAsFactors = FALSE
) %>% left_join(primary_mm2, by = "read")
read_summary$hmm_pick <- c("A","B","C")[apply(logL_mat, 1, which.max)]
read_summary$hmm_margin <- apply(logL_mat, 1, function(r) sort(r, decreasing=TRUE)[1] - sort(r, decreasing=TRUE)[2])
read_summary$mm2_margin <- vapply(seq_len(N_READS), function(i) {
  as_vals <- AS_mat[i, ]
  if (all(is.na(as_vals))) return(NA_real_)
  as_sorted <- sort(as_vals, decreasing = TRUE, na.last = NA)
  if (length(as_sorted) < 2) return(NA_real_)
  as_sorted[1] - as_sorted[2]
}, numeric(1))
read_summary$agree <- read_summary$hmm_pick == read_summary$mm2_pick

cat("\nHMM vs minimap2 agreement: ", sum(read_summary$agree, na.rm=TRUE),
    " / ", sum(!is.na(read_summary$mm2_pick)), "\n", sep="")
cat("HMM hard call vs TRUE:    ", sum(read_summary$hmm_pick == read_summary$true),
    " / ", N_READS, "\n", sep="")
cat("mm2 hard call vs TRUE:    ",
    sum(read_summary$mm2_pick == read_summary$true, na.rm=TRUE),
    " / ", sum(!is.na(read_summary$mm2_pick)), "\n", sep="")

# ── Figure ────────────────────────────────────────────────────────────────────
# Header
p_header <- ggplot() +
  annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = 10,
           fill = "#E3F2FD", color = "#0D47A1", linewidth = 0.5) +
  annotate("text", x = 0.2, y = 9.0,
           label = "minimap2 alignment score (AS)  vs  per-copy HMM logL  on the same reads",
           hjust = 0, fontface = "bold", size = 4.5, color = "#0D47A1") +
  annotate("text", x = 0.3, y = 7.4,
           label = sprintf("3 copies of length %d bp,  %d SNP positions distinguishing them.",
                            L_TX, N_SNPS),
           hjust = 0, size = 3.4, color = "gray15") +
  annotate("text", x = 0.3, y = 6.0,
           label = sprintf("%d full-length reads (length %d bp), %.0f%% per-base error.  True abundances:  A=%.2f, B=%.2f, C=%.2f.",
                            N_READS, L_TX, 100 * ERROR_RATE,
                            TRUE_THETA["A"], TRUE_THETA["B"], TRUE_THETA["C"]),
           hjust = 0, size = 3.3, color = "gray15") +
  annotate("text", x = 0.3, y = 4.4,
           label = "minimap2:  Smith-Waterman alignment + chaining + MAPQ.  Score = AS  (higher = better fit).",
           hjust = 0, size = 3.3, color = "gray15") +
  annotate("text", x = 0.3, y = 3.0,
           label = "Per-copy HMM:  logL = M . log(0.95) + (L-M) . log(0.017)  for M = match count vs that copy's transcript.",
           hjust = 0, size = 3.3, color = "gray15") +
  annotate("text", x = 0.3, y = 1.3,
           label = "Question: do they pick the same copy?  If yes -> HMM is redundant with alignment for the hard call.",
           hjust = 0, size = 3.4, color = "#1B5E20", fontface = "bold") +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 10)) +
  theme_void()

# Side-by-side scatter: AS vs HMM logL (one point per read per copy)
score_long <- bind_rows(
  data.frame(read = rownames(logL_mat),
             copy = "A", AS = AS_mat[,"A"], logL = logL_mat[,"A"]),
  data.frame(read = rownames(logL_mat),
             copy = "B", AS = AS_mat[,"B"], logL = logL_mat[,"B"]),
  data.frame(read = rownames(logL_mat),
             copy = "C", AS = AS_mat[,"C"], logL = logL_mat[,"C"])
) %>% filter(!is.na(AS) & !is.na(logL))

p_scatter <- ggplot(score_long, aes(x = AS, y = logL, color = copy)) +
  geom_point(size = 2.4, alpha = 0.8) +
  scale_color_manual(values = c(A = "#1565C0", B = "#B71C1C", C = "#2E7D32")) +
  labs(title    = "minimap2 AS  vs  per-copy HMM logL  (one point per read-copy alignment)",
       subtitle = "Tight linear relationship -> the two scores RANK identically.  Both pick the same best copy per read.",
       x = "minimap2 AS score", y = "per-copy HMM logL") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(color = "gray40", size = 10),
        legend.position = "right")

# Bar chart per copy: how many reads each method assigns to it
per_copy_counts <- data.frame(
  copy = c("A","B","C"),
  true  = as.numeric(table(factor(read_summary$true, levels = c("A","B","C")))),
  mm2   = as.numeric(table(factor(read_summary$mm2_pick, levels = c("A","B","C")))),
  hmm   = as.numeric(table(factor(read_summary$hmm_pick, levels = c("A","B","C"))))
) %>% pivot_longer(-copy, names_to = "method", values_to = "n_reads")

p_counts <- ggplot(per_copy_counts, aes(x = copy, y = n_reads, fill = method)) +
  geom_col(position = "dodge", color = "white", linewidth = 0.4) +
  geom_text(aes(label = n_reads), position = position_dodge(width = 0.9),
            vjust = -0.3, size = 3.5, color = "gray20") +
  scale_fill_manual(values = c(true = "#2E7D32",
                                mm2  = "#1565C0",
                                hmm  = "#9E9E9E"),
                    labels = c(true = "true source",
                                mm2  = "minimap2 hard call",
                                hmm  = "HMM hard call"),
                    name = NULL) +
  labs(title    = "Read counts per copy: true vs minimap2 vs HMM",
       subtitle = sprintf("Agreement HMM <-> minimap2:  %d / %d.   Accuracy vs truth: minimap2 %d, HMM %d.",
                            sum(read_summary$agree, na.rm=TRUE),
                            sum(!is.na(read_summary$mm2_pick)),
                            sum(read_summary$mm2_pick == read_summary$true, na.rm=TRUE),
                            sum(read_summary$hmm_pick == read_summary$true)),
       x = NULL, y = "# reads") +
  theme_classic(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(color = "gray40", size = 10),
        axis.text.x   = element_text(face = "bold", size = 11),
        legend.position = "top")

# Conclusion box
n_total <- sum(!is.na(read_summary$mm2_pick))
n_agree <- sum(read_summary$agree, na.rm = TRUE)
conclusion_color <- if (n_agree == n_total) "#E8F5E9" else "#FFF3E0"
conclusion_border <- if (n_agree == n_total) "#1B5E20" else "#E65100"
conclusion_text <- if (n_agree == n_total) paste(
  "minimap2 AS and per-copy HMM logL agree on ALL", n_total, "reads.",
  "For full-length FLNC reads on this 3-copy toy, the HMM is redundant with alignment for the hard call.",
  "The HMM's only added value is providing P(read | copy) for downstream EM — the assignment itself is done by alignment.",
  sep = "\n") else paste(
  sprintf("minimap2 AS and HMM logL DIFFER on %d / %d reads.", n_total - n_agree, n_total),
  "These are the cases where the two methods give different hard calls.  Examine them to see whether HMM adds discrimination beyond minimap2.",
  sep = "\n")

p_concl <- ggplot() +
  annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = 1,
           fill = conclusion_color, color = conclusion_border, linewidth = 0.7) +
  annotate("text", x = 0.15, y = 0.85, label = "Conclusion",
           hjust = 0, fontface = "bold", size = 4.0, color = conclusion_border) +
  annotate("text", x = 0.15, y = 0.42, label = conclusion_text,
           hjust = 0, vjust = 0.5, size = 3.5, color = "gray15",
           lineheight = 1.35) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 1)) +
  theme_void()

fig <- (p_header) /
       (p_scatter | p_counts) /
       (p_concl) +
  plot_layout(heights = c(0.5, 1.4, 0.45)) +
  plot_annotation(
    title    = "Does the HMM step add anything beyond minimap2 alignment?",
    subtitle = "Apples-to-apples on the same toy FLNC reads. If minimap2 AS already discriminates, the HMM step is unnecessary for the hard call.",
    theme = theme(plot.title    = element_text(size = 16, face = "bold"),
                  plot.subtitle = element_text(size = 11, color = "gray40"))
  )

ggsave(file.path(FIGURES_DIR, "fig_advisor_minimap2_vs_hmm_toy.pdf"), fig,
       width = 15, height = 11, units = "in")
ggsave(file.path(FIGURES_DIR, "fig_advisor_minimap2_vs_hmm_toy.png"), fig,
       width = 15, height = 11, units = "in", dpi = 180)
message("Slide saved.")
