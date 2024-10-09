library(iSEE)
library(scRNAseq)
library(scater)

# Example data ----

sce <- ReprocessedAllenData(assays="tophat_counts")

sce <- logNormCounts(sce, exprs_values="tophat_counts")

rowData(sce)$row_var <- rowVars(assay(sce, "logcounts"))
rowData(sce)$n_cells <- rowSums(assay(sce, "logcounts") > 0)

# launch the app itself ----

app <- iSEE(sce, initial = list(
  RowDataPlot(
    PanelWidth = 12L,
    YAxis = "row_var",
    XAxis = "Row data",
    XAxisRowData = "n_cells"
  )
))

if (interactive()) {
  shiny::runApp(app, port=1234)
}
