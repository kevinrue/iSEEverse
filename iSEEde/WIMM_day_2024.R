## ----"start", message=FALSE, warning=FALSE--------------------------------------------------------
library("iSEEde")
library("airway")
library("DESeq2")
library("iSEE")

# Example data ----

data("airway")
airway$dex <- relevel(airway$dex, "untrt")

dds <- DESeqDataSet(airway, ~ 0 + dex + cell)

dds <- DESeq(dds)
res_deseq2 <- results(dds, contrast = list("dextrt", "dexuntrt"))

library(dplyr)
res_deseq2@listData <- res_deseq2@listData %>%
  as.data.frame() %>%
  mutate(across(where(is.numeric), function(x) {
    as.numeric(format(x, digits = 2))
  })) %>%
  as.list()

head(res_deseq2)

# iSEE / iSEEde ---

airway <- embedContrastResults(res_deseq2, airway, name = "dex: trt vs untrt")
contrastResults(airway)

rownames(airway) <- scater::uniquifyFeatureNames(
  ID = rowData(airway)[["gene_id"]],
  names = rowData(airway)[["gene_name"]]
)

panelDefaults(TooltipRowData = c("gene_name", "gene_id", "gene_biotype"))

app <- iSEE(airway, initial = list(
  DETable(ContrastName="dex: trt vs untrt", HiddenColumns = c("baseMean",
    "lfcSE", "stat"), PanelWidth = 4L),
  VolcanoPlot(ContrastName="dex: trt vs untrt", PanelWidth = 4L),
  MAPlot(ContrastName="dex: trt vs untrt", PanelWidth = 4L)
))

if (interactive()) {
  shiny::runApp(app)
}
