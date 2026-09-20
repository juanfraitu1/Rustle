#!/usr/bin/env bash
set -euo pipefail

mamba install -y \
  -c bioconda -c conda-forge \
  bioconductor-rtracklayer \
  bioconductor-genomicranges \
  r-igraph \
  r-ggraph \
  r-ggplot2 \
  r-dplyr \
  r-tidyr \
  r-stringr \
  r-patchwork \
  r-testthat
