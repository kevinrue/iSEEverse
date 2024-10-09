library(iSEE)
library(scRNAseq)
library(scater)

# Example data ----
sce <- ReprocessedAllenData(assays="tophat_counts")

sce <- logNormCounts(sce, exprs_values="tophat_counts")

# launch the app itself ----

app <- iSEE(sce, initial = list(
  FeatureAssayPlot(
    PanelWidth = 12L,
    YAxisFeatureName = "Rorb",
    XAxis = "Column data", XAxisColumnData = "driver_1_s",
    ColorBy = "Column data", ColorByColumnData = "driver_1_s",
    FacetColumnBy = "Column data", FacetColumnByColData = "Core.Type"
  )
))

if (interactive()) {
  shiny::runApp(app, port=1234)
}
