library(iSEE)
library(iSEEde)
library(airway)
library(DESeq2)

# Example data ----

data("airway")
airway$dex <- relevel(airway$dex, "untrt")
rowData(airway)$seq_strand <- factor(rowData(airway)$seq_strand)

dds <- DESeqDataSet(airway, ~ 0 + dex + cell)

dds <- DESeq(dds)
res_deseq2 <- results(dds, contrast = list("dextrt", "dexuntrt"))

# iSEE / iSEEde ---

airway <- embedContrastResults(res_deseq2, airway, name = "dex: trt vs untrt")

app <- iSEE(airway, initial = list(
  VolcanoPlot(
    PanelWidth = 12L,
    ContrastName="dex: trt vs untrt",
    ColorBy = "Row data",
    ColorByRowData = "seq_strand"
  )
))

if (interactive()) {
  shiny::runApp(app)
}
