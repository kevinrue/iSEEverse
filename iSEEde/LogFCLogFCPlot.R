library("iSEEde")
library("airway")
library("edgeR")
library("DESeq2")
library("iSEE")

# Example data ----

data("airway")
airway$dex <- relevel(airway$dex, "untrt")

# DESeq2 ----

dds <- DESeqDataSet(airway, ~ 0 + dex + cell)

dds <- DESeq(dds)
res_deseq2 <- results(dds, contrast = list("dextrt", "dexuntrt"))

airway <- embedContrastResults(res_deseq2, airway, name = "DESeq2")

# edgeR ----

design <- model.matrix(~ 0 + dex + cell, data = colData(airway))

fit <- glmFit(airway, design, dispersion = 0.1)
lrt <- glmLRT(fit, contrast = c(-1, 1, 0, 0, 0))
res_edger <- topTags(lrt, n = Inf)

airway <- embedContrastResults(res_edger, airway, name = "edgeR")

# iSEE / iSEEde ---

airway <- registerAppOptions(airway, factor.maxlevels = 30, color.maxlevels = 30)

app <- iSEE(airway, initial = list(
  LogFCLogFCPlot(
    ContrastNameX = "DESeq2", ContrastNameY = "edgeR",
    ColorBy = "Row data",
    ColorByRowData = "gene_biotype",
    BrushData = list(
      xmin = 3.6, xmax = 8.2, ymin = 3.8, ymax = 9.8,
      mapping = list(x = "X", y = "Y"),
      direction = "xy", brushId = "LogFCLogFCPlot1_Brush",
      outputId = "LogFCLogFCPlot1"),
    PanelWidth = 8L)
))

if (interactive()) {
  shiny::runApp(app)
}
