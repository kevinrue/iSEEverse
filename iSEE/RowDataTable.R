library(iSEE)
library(scRNAseq)

# Example data ----

sce <- ReprocessedAllenData(assays="tophat_counts")

rowData(sce)$ave_count <- rowMeans(assay(sce, "tophat_counts"))
rowData(sce)$n_cells <- rowSums(assay(sce, "tophat_counts") > 0)

# launch the app itself ----

app <- iSEE(sce, initial = list(
  RowDataTable(
    PanelWidth = 12L
  )
))

if (interactive()) {
  shiny::runApp(app, port=1234)
}
