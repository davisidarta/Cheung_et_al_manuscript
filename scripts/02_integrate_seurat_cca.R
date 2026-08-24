#!/usr/bin/env Rscript
# Anchor-based integration of the ten sympathetic ganglia datasets.
#
# Implements the canonical correlation analysis workflow of Stuart et al. (2019): shared
# variable features are selected across datasets, mutual nearest neighbours in CCA space
# are used as anchors, and the datasets are corrected onto a common space.
#
#   Rscript scripts/02_integrate_seurat_cca.R --indir data/processed --output data/integrated.rds
#
# This is the expensive step. With ~239k cells across ten datasets it needs a large-memory
# machine and several hours.

suppressPackageStartupMessages({
  library(Seurat)
  library(future)
  library(optparse)
})

# --- Integration parameters ---------------------------------------------------------

N_INTEGRATION_FEATURES <- 2000  # shared features used to find anchors
CCA_DIMS <- 1:30                # canonical correlation vectors searched for anchors
K_ANCHOR <- 5
K_FILTER <- 200
K_WEIGHT <- 100
N_PCS <- 50                     # principal components computed after correction

MAX_GLOBALS_GB <- 200


option_list <- list(
  make_option("--indir", default = "data/processed"),
  make_option("--output", default = "data/integrated.rds"),
  make_option("--threads", default = 8L, type = "integer")
)
opt <- parse_args(OptionParser(option_list = option_list))

plan("multicore", workers = opt$threads)
options(future.globals.maxSize = MAX_GLOBALS_GB * 1024^3)

files <- list.files(opt$indir, pattern = "_qc\\.rds$", full.names = TRUE)
if (length(files) == 0) {
  stop("no *_qc.rds files in ", opt$indir, " -- run 01_qc_and_normalise.R first")
}

message("loading ", length(files), " datasets")
objects <- lapply(files, readRDS)
names(objects) <- sub("_qc\\.rds$", "", basename(files))
for (n in names(objects)) {
  message(sprintf("  %-26s %6d cells", n, ncol(objects[[n]])))
}
message(sprintf("  total: %d cells", sum(vapply(objects, ncol, integer(1)))))

# --- Anchor-based correction --------------------------------------------------------

message("\nselecting integration features")
features <- SelectIntegrationFeatures(
  object.list = objects,
  nfeatures = N_INTEGRATION_FEATURES
)

message("finding anchors (this is the slow part)")
anchors <- FindIntegrationAnchors(
  object.list = objects,
  anchor.features = features,
  reduction = "cca",
  dims = CCA_DIMS,
  k.anchor = K_ANCHOR,
  k.filter = K_FILTER
)

message("integrating")
integrated <- IntegrateData(
  anchorset = anchors,
  dims = CCA_DIMS,
  k.weight = K_WEIGHT
)

DefaultAssay(integrated) <- "integrated"

# The corrected matrix is centred and scaled before dimensionality reduction; this is
# the matrix the downstream embedding is computed from.
message("scaling and running PCA")
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = N_PCS, verbose = FALSE)

message(sprintf(
  "\nintegrated: %d cells x %d features",
  ncol(integrated), nrow(integrated)
))

saveRDS(integrated, opt$output)
message("wrote ", opt$output)
