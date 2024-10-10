library(iSEE)
library(iSEEhex)
library(scRNAseq)
library(scater)

# Example data ----

sce <- ReprocessedAllenData(assays="tophat_counts")

sce <- logNormCounts(sce, exprs_values="tophat_counts")

sce <- runPCA(sce, ncomponents=4)

# launch the app itself ----

if (interactive()) {
  iSEE(sce, initial=list(
    ReducedDimensionHexPlot(PanelWidth = 12L, BinResolution=50)
  ))
}
