library(iSEE)
library(iSEEde)
library(airway)
library(DESeq2)

# Example data ----

data("airway")
airway$dex <- relevel(airway$dex, "untrt")

dds <- DESeqDataSet(airway, ~ 0 + dex + cell)

dds <- DESeq(dds)
res_deseq2 <- results(dds, contrast = list("dextrt", "dexuntrt"))

# iSEE / iSEEde ---

airway <- embedContrastResults(res_deseq2, airway, name = "dex: trt vs untrt")

app <- iSEE(airway, initial = list(
  DETable(
    PanelWidth = 12L,
    ContrastName="dex: trt vs untrt",
    RoundDigits = TRUE
  )
))

if (interactive()) {
  shiny::runApp(app)
}
