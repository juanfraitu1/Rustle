#!/usr/bin/env Rscript
# bambu_run.R BAM GENOME.fa OUT_PREFIX NCORE RC_DIR [GUIDE.gtf]   (register 1641)
# Without GUIDE: de novo (annotations=NULL). With GUIDE: guided. Read classes are built from the BAM once PER MODE in RC_DIR
# (never share RC_DIR between modes: read-class construction uses the annotation). Two discovery calls on those read classes:
# bambu's own default NDR -> OUT.default.gtf, and NDR=1 -> OUT.gtf + OUT.tsv (NDR per transcript, to cut afterwards).
suppressMessages(library(bambu))
a <- commandArgs(TRUE); stopifnot(length(a) >= 5)
fa <- Rsamtools::FaFile(a[2])  # a path string makes bambu query BSgenome::available.genomes() (needs network)
guide <- if (length(a) >= 6) prepareAnnotations(a[6]) else NULL
nc <- as.integer(a[4]); dir.create(a[5], recursive = TRUE, showWarnings = FALSE)
invisible(bambu(reads = a[1], annotations = guide, genome = fa, rcOutDir = a[5], discovery = FALSE, quant = FALSE, ncore = nc, lowMemory = TRUE, verbose = TRUE))
rc <- list.files(a[5], pattern = "\\.rds$", full.names = TRUE); stopifnot(length(rc) == 1)
tsv <- function(x, p) { m <- as.data.frame(mcols(x)); write.table(m[, !sapply(m, is.list), drop = FALSE], p, sep = "\t", quote = FALSE, row.names = FALSE) }
d <- bambu(reads = rc, annotations = guide, genome = fa, quant = FALSE, ncore = nc, verbose = TRUE)
writeToGTF(d, paste0(a[3], ".default.gtf")); tsv(d, paste0(a[3], ".default.tsv"))
x <- bambu(reads = rc, annotations = guide, genome = fa, NDR = 1, quant = FALSE, ncore = nc, verbose = TRUE)
saveRDS(x, paste0(a[3], ".rds")); writeToGTF(x, paste0(a[3], ".gtf")); tsv(x, paste0(a[3], ".tsv"))
cat("bambu", as.character(packageVersion("bambu")), "default transcripts", length(d), "novel", sum(mcols(d)$novelTranscript),
    "| NDR=1 transcripts", length(x), "novel", sum(mcols(x)$novelTranscript), "\n")
