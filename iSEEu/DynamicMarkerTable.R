library(iSEE)
library(iSEEu)
library(scRNAseq)
library(scater)
library(scran)

sce <- ReprocessedAllenData(assays="tophat_counts")
sce <- logNormCounts(sce, exprs_values="tophat_counts")
sce <- runPCA(sce, ncomponents=4)

if (interactive()) {
  iSEE(sce, initial=list(
    ReducedDimensionPlot(
      PanelWidth=4L,
      BrushData = list(
        lasso = NULL, closed = TRUE,
        mapping = list(x = "X", y = "Y"),
        coord = structure(c(
          -47.8, -41.9, -14.6, -13.6, -19.1, -27.3, -33.6, -44, -47.8,
          -23.6, -44.1, -56.4, -46.9, -26.4, -17.4, -6.2, -5.4, -23.6),
          dim = c(9L, 2L))
        )
    ),
    DynamicMarkerTable(
      PanelWidth=8L,
      ColumnSelectionSource="ReducedDimensionPlot1"
    )
  ))
}
