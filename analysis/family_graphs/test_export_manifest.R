source("02_build_graphs.R")

# Minimal synthetic exon_df: 2 genes, each 2 exons, different chromosomes
syn_exon_df <- data.frame(
  gene_id = c("gA","gA","gB","gB"),
  tx_id   = c("tA","tA","tB","tB"),
  chrom   = c("chr1","chr1","chr2","chr2"),
  start   = c(100L,300L,100L,300L),
  end     = c(200L,400L,200L,400L),
  strand  = "+",
  stringsAsFactors = FALSE
)
tmp <- tempfile(fileext = ".tsv")
export_family_manifest(syn_exon_df, "TESTFAM", tmp)
tbl <- read.table(tmp, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

stopifnot(all(c("family_id","gene_id","chrom","start","end","strand") %in% names(tbl)))
stopifnot(nrow(tbl) == 2L)
stopifnot(all(tbl$family_id == "TESTFAM"))
stopifnot(tbl$start[tbl$gene_id == "gA"] == 100L)
stopifnot(tbl$end[tbl$gene_id == "gA"]   == 400L)
stopifnot(tbl$chrom[tbl$gene_id == "gA"] == "chr1")
message("OK: export_family_manifest test passed")
