library(iSEE)
library(iSEEu)
library(scRNAseq)
library(scater)

sce <- ReprocessedAllenData(assays="tophat_counts")
sce <- logNormCounts(sce, exprs_values="tophat_counts")
sce <- runPCA(sce, ncomponents=4)

if (interactive()) {
  iSEE(sce, initial=list(
    ReducedDimensionPlot(
      PanelWidth = 6L
    ),
    DynamicReducedDimensionPlot(
      PanelWidth = 6L,
      Assay="logcounts",
      ColumnSelectionSource="ReducedDimensionPlot1"
    )
  ))
}
