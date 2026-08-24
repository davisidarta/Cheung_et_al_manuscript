#!/usr/bin/env Rscript
# Select neurons from the integrated atlas and export for the Python analysis.
#
# Neurons are identified by expression of the pan-neuronal markers Tubb3 (tubulin beta-3)
# and Rbfox3 (which encodes NeuN). Clusters are scored on both markers and those above
# threshold are retained; the remaining cells are satellite glia, Schwann cells, immune
# and vascular populations.
#
#   Rscript scripts/03_select_neurons_export.R --input data/integrated.rds --output data/neurons.h5ad

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(optparse)
})

# --- Neuron selection ---------------------------------------------------------------

NEURONAL_MARKERS <- c("Tubb3", "Rbfox3")

#: Clustering resolution used to define the units that are scored for neuronal identity.
CLUSTER_RESOLUTION <- 0.8
N_PCS <- 50

#: A cluster is called neuronal when this share of its cells detect a marker.
MIN_DETECTION_FRACTION <- 0.30


option_list <- list(
  make_option("--input", default = "data/integrated.rds"),
  make_option("--output", default = "data/neurons.h5ad"),
  make_option("--barcodes", default = "data/neuron_barcodes.txt")
)
opt <- parse_args(OptionParser(option_list = option_list))

message("reading ", opt$input)
integrated <- readRDS(opt$input)
DefaultAssay(integrated) <- "integrated"

message("clustering at resolution ", CLUSTER_RESOLUTION)
integrated <- FindNeighbors(integrated, dims = seq_len(N_PCS), verbose = FALSE)
integrated <- FindClusters(integrated, resolution = CLUSTER_RESOLUTION, verbose = FALSE)

# Marker detection is assessed on the uncorrected RNA assay: the integrated assay holds
# corrected values that are not counts and should not be thresholded at zero.
DefaultAssay(integrated) <- "RNA"

present <- NEURONAL_MARKERS[NEURONAL_MARKERS %in% rownames(integrated)]
if (length(present) == 0) {
  stop("none of the pan-neuronal markers are present: ",
       paste(NEURONAL_MARKERS, collapse = ", "))
}
message("scoring clusters on: ", paste(present, collapse = ", "))

expr <- GetAssayData(integrated, layer = "data")[present, , drop = FALSE]
clusters <- Idents(integrated)

detection <- sapply(levels(clusters), function(cl) {
  cells <- which(clusters == cl)
  rowMeans(expr[, cells, drop = FALSE] > 0)
})
detection <- t(as.matrix(detection))
colnames(detection) <- present

neuronal_clusters <- rownames(detection)[
  apply(detection, 1, function(r) all(r >= MIN_DETECTION_FRACTION))
]

message("\ncluster marker detection:")
for (cl in rownames(detection)) {
  flag <- if (cl %in% neuronal_clusters) " <- neuronal" else ""
  message(sprintf("  cluster %-4s %s%s", cl,
                  paste(sprintf("%s=%.2f", present, detection[cl, ]), collapse = "  "),
                  flag))
}

neurons <- subset(integrated, idents = neuronal_clusters)
message(sprintf(
  "\n%d of %d cells retained as neurons (%.1f%%)",
  ncol(neurons), ncol(integrated), 100 * ncol(neurons) / ncol(integrated)
))

writeLines(colnames(neurons), opt$barcodes)
message("wrote ", opt$barcodes)

# --- Export -------------------------------------------------------------------------
# The Python scripts read the corrected matrix from X and the log-normalised counts from
# raw, so both assays are carried across.

message("converting to h5ad")
tmp <- sub("\\.h5ad$", ".h5Seurat", opt$output)
SaveH5Seurat(neurons, filename = tmp, overwrite = TRUE)
Convert(tmp, dest = "h5ad", overwrite = TRUE)
unlink(tmp)

message("wrote ", opt$output)
