library(iSEE)
library(iSEEu)
library(scRNAseq)
library(scater)
library(scran)

sce <- ReprocessedAllenData(assays="tophat_counts")
sce <- logNormCounts(sce, exprs_values="tophat_counts")
sce <- runPCA(sce, ncomponents=4)
sce <- runTSNE(sce)

if (interactive()) {
  iSEE(sce, initial=list(
    ReducedDimensionPlot(
      PanelWidth=4L,
      BrushData = list(
        lasso = NULL, closed = TRUE,
        mapping = list(x = "X", y = "Y"),
        coord = structure(c(
          -13.6, -18.8, -33.9, -44.2, -47.1, -13.6,
          -52.3, -27.3, -8.5, -6.3, -31.6, -52.3
          ), dim = c(6L, 2L))
        )
    ),
    DynamicMarkerTable(
      PanelWidth=8L,
      ColumnSelectionSource="ReducedDimensionPlot1"
    )
  ))
}
