library(iSEE)
library(scRNAseq)
library(scater)

# Example data ----

sce <- ReprocessedAllenData(assays="tophat_counts")

sce <- logNormCounts(sce, exprs_values="tophat_counts")

rowData(sce)$row_var <- rowVars(assay(sce, "logcounts"))

# launch the app itself ----

app <- iSEE(sce, initial = list(
  SampleAssayPlot(
    PanelWidth = 12L,
    YAxisSampleName = "SRR2140028",
    XAxis = "Sample name", XAxisSampleName = "SRR2140022",
    ColorBy = "Row data", ColorByRowData = "row_var"
  )
))

if (interactive()) {
  shiny::runApp(app, port=1234)
}
