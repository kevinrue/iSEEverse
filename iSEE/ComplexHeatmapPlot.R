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
sce

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

app <- iSEE(sce, initial = list(
  ComplexHeatmapPlot(
    PanelWidth = 12L,
    CustomRows = TRUE,
    CustomRowsText = "Lamp5\nFam19a1\nCnr1\nRorb\nSparcl1\nCrym\nLmo3\nSerpine2\nDdah1\nCux2\n"
  )
))

if (interactive()) {
  shiny::runApp(app, port=1234)
}

