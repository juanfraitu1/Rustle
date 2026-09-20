library(testthat)
source("../02_build_graphs.R")

# Minimal fixture: 2 genes, one shared exon (100-200), one unique each
make_fixture <- function() {
  data.frame(
    gene_id = c("gA","gA","gB","gB"),
    tx_id   = c("tA","tA","tB","tB"),
    chrom   = rep("chr1", 4),
    start   = c(100L, 300L, 100L, 500L),
    end     = c(200L, 400L, 200L, 600L),
    strand  = rep("+", 4),
    stringsAsFactors = FALSE
  )
}

test_that("shared node has n_copies >= 2", {
  vg <- build_variation_graph(make_fixture())
  shared <- vg$nodes[vg$nodes$node_type == "shared", ]
  expect_true(nrow(shared) >= 1)
  expect_true(all(shared$n_copies >= 2))
})

test_that("copy-specific node belongs to exactly one gene", {
  vg <- build_variation_graph(make_fixture())
  specific <- vg$nodes[vg$nodes$node_type == "copy_specific", ]
  expect_true(nrow(specific) >= 2)   # gA has one, gB has one
  expect_true(all(specific$n_copies == 1))
})

test_that("each copy path traverses at least one copy-specific node", {
  vg <- build_variation_graph(make_fixture())
  for (gene in names(vg$paths)) {
    path_nodes <- vg$nodes[vg$nodes$node_id %in% vg$paths[[gene]], ]
    own_specific <- path_nodes[
      path_nodes$node_type == "copy_specific" &
      sapply(path_nodes$gene_set, function(gs) gene %in% gs), ]
    expect_true(nrow(own_specific) >= 1, info = paste("gene:", gene))
  }
})

test_that("is_multicopy_family() returns FALSE for single gene", {
  single <- make_fixture()[1:2, ]
  single$gene_id <- "gA"
  vg <- build_variation_graph(single)
  expect_false(is_multicopy_family(vg))
})

test_that("is_multicopy_family() returns FALSE for two unrelated genes (no shared exon)", {
  no_share <- data.frame(
    gene_id = c("gA","gB"),
    tx_id   = c("tA","tB"),
    chrom   = rep("chr1", 2),
    start   = c(100L, 1000L),
    end     = c(200L, 1100L),
    strand  = rep("+", 2),
    stringsAsFactors = FALSE
  )
  vg <- build_variation_graph(no_share)
  expect_false(is_multicopy_family(vg))
})

test_that("is_multicopy_family() returns TRUE for valid family", {
  vg <- build_variation_graph(make_fixture())
  expect_true(is_multicopy_family(vg))
})

test_that("graph is a DAG", {
  vg <- build_variation_graph(make_fixture())
  expect_true(igraph::is_dag(vg$graph))
})

test_that("minus-strand edges run in transcription order (3prime-to-5prime genomic)", {
  minus_df <- data.frame(
    gene_id = c("gA","gA","gB","gB"),
    tx_id   = c("tA","tA","tB","tB"),
    chrom   = rep("chr1", 4),
    start   = c(100L, 300L, 100L, 500L),
    end     = c(200L, 400L, 200L, 600L),
    strand  = rep("-", 4),
    stringsAsFactors = FALSE
  )
  vg <- build_variation_graph(minus_df)
  # For minus strand, path should go from high node_id to low (genomic 3'->5')
  # meaning edges in paths go from higher to lower node_id
  path_a <- vg$paths[["gA"]]
  expect_true(path_a[1] > path_a[length(path_a)])
})
