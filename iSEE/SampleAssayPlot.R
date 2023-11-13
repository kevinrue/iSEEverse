library(iSEE)

library(scRNAseq)

# Example data ----
sce <- ReprocessedAllenData(assays="tophat_counts")
class(sce)

library(scater)
sce <- logNormCounts(sce, exprs_values="tophat_counts")

sce <- runPCA(sce, ncomponents=4)
sce <- runTSNE(sce)
sce <- runUMAP(sce)
rowData(sce)$ave_count <- rowMeans(assay(sce, "tophat_counts"))
rowData(sce)$n_cells <- rowSums(assay(sce, "tophat_counts") > 0)
rowData(sce)$row_var <- rowVars(assay(sce, "logcounts"))
sce

# launch the app itself ----

app <- iSEE(sce, initial = list(
  SampleAssayPlot(
    PanelWidth = 8L,
    YAxisSampleName = "SRR2140028",
    XAxis = "Sample name", XAxisSampleName = "SRR2140022",
    ColorBy = "Row data", ColorByRowData = "row_var"
  )
))

if (interactive()) {
  shiny::runApp(app, port=1234)
}

