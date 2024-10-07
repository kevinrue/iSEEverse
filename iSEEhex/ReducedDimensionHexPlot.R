library(scRNAseq)

# Example data ----
sce <- ReprocessedAllenData(assays="tophat_counts")
class(sce)

library(scater)
sce <- logNormCounts(sce, exprs_values="tophat_counts")

sce <- runPCA(sce, ncomponents=4)
sce <- runTSNE(sce)
rowData(sce)$ave_count <- rowMeans(assay(sce, "tophat_counts"))
rowData(sce)$n_cells <- rowSums(assay(sce, "tophat_counts") > 0)

# launch the app itself ----

if (interactive()) {
  iSEE(sce, initial=list(
    ReducedDimensionHexPlot(BinResolution=50, PanelWidth = 8L),
    ReducedDimensionPlot(PanelWidth = 4L)
  ))
}
