#!/usr/bin/env Rscript
# Per-dataset quality control and log-normalisation.
#
# Reads each of the ten public sympathetic ganglia datasets listed in
# data/sample_manifest.csv, applies the same quality control to each, and writes one
# normalised Seurat object per study for the integration step to consume.
#
#   Rscript scripts/01_qc_and_normalise.R --indir data/raw --outdir data/processed
#
# Input layout: data/raw/<study>/ containing either 10x-style matrix files
# (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz) or a single count matrix; see
# scripts/00_fetch_source_data.sh, which arranges the downloads this way.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(optparse)
})

# --- Quality control thresholds -----------------------------------------------------

MIN_CELLS <- 3        # a gene must be seen in at least this many cells
MIN_FEATURES <- 200   # a cell must express at least this many genes
MAX_PERCENT_MT <- 20  # cells above this share of mitochondrial counts are discarded
MIN_COUNTS <- 500

N_VARIABLE_FEATURES <- 2000
SCALE_FACTOR <- 1e4


option_list <- list(
  make_option("--indir", default = "data/raw"),
  make_option("--outdir", default = "data/processed"),
  make_option("--manifest", default = "data/sample_manifest.csv")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
manifest <- read.csv(opt$manifest, stringsAsFactors = FALSE)


#' Read a study directory into a count matrix.
#'
#' Handles the two layouts the download script produces: a 10x triplet, or a dense or
#' sparse matrix serialised as .rds.
read_counts <- function(path) {
  if (file.exists(file.path(path, "matrix.mtx.gz")) ||
      file.exists(file.path(path, "matrix.mtx"))) {
    return(Read10X(path))
  }
  rds <- list.files(path, pattern = "\\.rds$", full.names = TRUE, ignore.case = TRUE)
  if (length(rds) >= 1) {
    return(readRDS(rds[1]))
  }
  stop("no recognised count matrix in ", path)
}


#' Quality control and log-normalisation for one study.
process_study <- function(counts, study, ganglion) {
  obj <- CreateSeuratObject(
    counts = counts,
    project = study,
    min.cells = MIN_CELLS,
    min.features = MIN_FEATURES
  )

  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^Rp[sl]")
  obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = "^Hb[ab]-")

  n_before <- ncol(obj)
  obj <- subset(
    obj,
    subset = nFeature_RNA >= MIN_FEATURES &
      nCount_RNA >= MIN_COUNTS &
      percent.mt <= MAX_PERCENT_MT
  )
  message(sprintf("    %s: %d -> %d cells after QC", study, n_before, ncol(obj)))

  obj$Study <- study
  obj$Ganglion <- ganglion
  obj$Species <- "Mouse"

  obj <- NormalizeData(
    obj,
    normalization.method = "LogNormalize",
    scale.factor = SCALE_FACTOR,
    verbose = FALSE
  )
  obj <- FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = N_VARIABLE_FEATURES,
    verbose = FALSE
  )
  obj
}


message("processing ", nrow(manifest), " studies")

for (i in seq_len(nrow(manifest))) {
  study <- manifest$study[i]
  ganglion <- manifest$ganglion[i]
  path <- file.path(opt$indir, study)

  if (!dir.exists(path)) {
    warning("skipping ", study, ": ", path, " not found")
    next
  }

  message("  ", study)
  counts <- read_counts(path)
  obj <- process_study(counts, study, ganglion)

  out <- file.path(opt$outdir, paste0(study, "_qc.rds"))
  saveRDS(obj, out)
  message("    wrote ", out)
}

message("done.")
