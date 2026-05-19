suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(dplyr)
  library(stringr)
})

GFF_PATH  <- "/mnt/c/Users/jfris/Desktop/GGO_genomic.gff"
DATA_DIR  <- "data"

gff_cache <- file.path(DATA_DIR, "gff_cache.rds")

if (!file.exists(gff_cache)) {
  message("Importing GFF (this takes ~2 min)...")
  gff <- import(GFF_PATH)
  saveRDS(gff, gff_cache)
  message("Cached to ", gff_cache)
} else {
  message("Loading cached GFF...")
  gff <- readRDS(gff_cache)
}
message("GFF loaded: ", length(gff), " features")

# GFF3 Parent attributes come as CharacterList; unwrap to plain character
first_parent <- function(x) {
  sapply(as.list(x), function(p) if (length(p) > 0L) p[[1L]] else NA_character_)
}

parse_family <- function(gff, description_pattern) {
  # -- genes --
  genes <- gff[gff$type == "gene"]
  desc  <- mcols(genes)$description
  desc[is.na(desc)] <- ""
  family_genes <- genes[grepl(description_pattern, desc, ignore.case = TRUE)]
  gene_ids     <- family_genes$ID
  message("  genes found: ", length(gene_ids))

  # -- transcripts (mRNA or transcript) --
  txs <- gff[gff$type %in% c("mRNA", "transcript")]
  tx_parent <- first_parent(txs$Parent)
  family_txs <- txs[!is.na(tx_parent) & tx_parent %in% gene_ids]
  tx_ids     <- family_txs$ID

  # tx -> gene lookup
  tx_gene <- setNames(
    tx_parent[!is.na(tx_parent) & tx_parent %in% gene_ids],
    tx_ids
  )

  # -- exons --
  exons <- gff[gff$type == "exon"]
  ex_parent <- first_parent(exons$Parent)
  family_exons <- exons[!is.na(ex_parent) & ex_parent %in% tx_ids]

  exon_df <- data.frame(
    gene_id = tx_gene[first_parent(family_exons$Parent)],
    tx_id   = first_parent(family_exons$Parent),
    chrom   = as.character(seqnames(family_exons)),
    start   = start(family_exons),
    end     = end(family_exons),
    strand  = as.character(strand(family_exons)),
    stringsAsFactors = FALSE
  )
  exon_df <- exon_df[!is.na(exon_df$gene_id), ]
  message("  exons: ", nrow(exon_df), " across ", length(unique(exon_df$gene_id)), " genes")
  list(genes = family_genes, exon_df = exon_df)
}

message("=== AMY family ===")
amy <- parse_family(gff, "amylase")
saveRDS(amy, file.path(DATA_DIR, "amy_family.rds"))
