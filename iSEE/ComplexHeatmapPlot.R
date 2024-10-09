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

rowData(sce)$row_var <- rowVars(assay(sce, "logcounts"))

library(tibble)
library(dplyr)
rowData(sce) %>%
  as.data.frame() %>%
  rownames_to_column("rownames") %>%
  as_tibble() %>%
  arrange(desc(row_var)) %>%
  slice_head(n = 10) %>%
  pull(rownames) %>%
  paste0(collapse = "\n")

# launch the app itself ----

gene_list <- c("Lamp5", "Fam19a1", "Cnr1", "Rorb", "Sparcl1", "Crym", "Lmo3",  "Serpine2", "Ddah1", "Cux2")

app <- iSEE(sce, initial = list(
  ComplexHeatmapPlot(
    PanelWidth = 12L,
    CustomRows = TRUE,
    CustomRowsText = paste0(paste0(gene_list, collapse = "\n"), "\n"),
    ColumnData = "driver_1_s"
  )
))

if (interactive()) {
  shiny::runApp(app, port=1234)
}
