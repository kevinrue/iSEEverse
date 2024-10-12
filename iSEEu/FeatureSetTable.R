library(iSEE)
library(iSEEu)
library(scRNAseq)
library(scater)
library(scran)
library(org.Mm.eg.db)

sce <- LunSpikeInData(location=FALSE)

sce <- logNormCounts(sce)

rowData(sce) <- cbind(rowData(sce), modelGeneVarWithSpikes(sce, "ERCC"))

cmds <- createGeneSetCommands(collections="GO",
  organism="org.Mm.eg.db", identifier="ENSEMBL")
sce <- registerFeatureSetCommands(sce, cmds)

# Setting up the application.

gst <- FeatureSetTable(
  Selected = "GO:0002020"
)

rdp <- RowDataPlot(
  YAxis="total",
  XAxis="Row data", XAxisRowData="mean",
  ColorBy="Row selection",
  RowSelectionSource="FeatureSetTable1"
)

rdt <- RowDataTable(
  RowSelectionSource="FeatureSetTable1"
)

if (interactive()) {
  iSEE(sce, initial=list(gst, rdp, rdt))
}
