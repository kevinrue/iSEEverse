library(iSEE)

library(scRNAseq)

# Example data ----
sce <- ReprocessedAllenData(assays="tophat_counts")

library(scater)
sce <- logNormCounts(sce, exprs_values="tophat_counts")

sce <- runPCA(sce, ncomponents=4)
sce <- runTSNE(sce)
sce <- runUMAP(sce)
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
