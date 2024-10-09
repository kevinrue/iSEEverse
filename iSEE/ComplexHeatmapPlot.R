library(iSEE)
library(scRNAseq)
library(scater)
library(tibble)
library(dplyr)

# Example data ----
sce <- ReprocessedAllenData(assays="tophat_counts")

sce <- logNormCounts(sce, exprs_values="tophat_counts")

rowData(sce)$ave_count <- rowMeans(assay(sce, "tophat_counts"))
rowData(sce)$n_cells <- rowSums(assay(sce, "tophat_counts") > 0)
rowData(sce)$row_var <- rowVars(assay(sce, "logcounts"))

# launch the app itself ----

# top 10 genes with highest variance in logcounts
gene_list <- c("Lamp5", "Fam19a1", "Cnr1", "Rorb", "Sparcl1", "Crym", "Lmo3",  "Serpine2", "Ddah1", "Cux2")

app <- iSEE(sce, initial = list(
  ComplexHeatmapPlot(
    PanelWidth = 12L,
    CustomRows = TRUE,
    CustomRowsText = paste0(paste0(gene_list, collapse = "\n"), "\n"),
    ColumnData = "driver_1_s",
    RowData = "row_var"
  )
))

if (interactive()) {
  shiny::runApp(app, port=1234)
}
