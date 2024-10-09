library(iSEE)
library(scRNAseq)
library(scater)

# Example data ----
sce <- ReprocessedAllenData(assays="tophat_counts")

sce <- logNormCounts(sce, exprs_values="tophat_counts")

# launch the app itself ----

app <- iSEE(sce, initial = list(
  ColumnDataTable(
    PanelWidth = 12L
  )
))

if (interactive()) {
  shiny::runApp(app, port=1234)
}
