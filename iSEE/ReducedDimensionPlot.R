library(iSEE)
library(scRNAseq)
library(scater)

# Example data ----

sce <- ReprocessedAllenData(assays="tophat_counts")

sce <- logNormCounts(sce, exprs_values="tophat_counts")

sce <- runPCA(sce, ncomponents=4)
sce <- runUMAP(sce)

# launch the app itself ----

app <- iSEE(sce, initial = list(
  ReducedDimensionPlot(
    PanelWidth = 8L,
    Type = "UMAP",
    ColorBy = "Column data", ColorByColumnData = "driver_1_s")))

if (interactive()) {
  shiny::runApp(app, port=1234)
}
