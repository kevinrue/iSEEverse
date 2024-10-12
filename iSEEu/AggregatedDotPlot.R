library(iSEE)
library(iSEEu)
library(scRNAseq)
library(scater)

# Example data ----

sce <- ReprocessedAllenData(assays="tophat_counts")

sce <- logNormCounts(sce, exprs_values="tophat_counts")

# launch the app itself ----

if (interactive()) {
  iSEE(
    sce,
    initial = list(
      AggregatedDotPlot(
        ColumnDataLabel="Primary.Type",
        CustomRowsText = "Rorb\nSnap25\nFoxp2",
        # PanelHeight = 500L,
        PanelWidth = 12L
      )
    )
  )
}
