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

initial <- list()

################################################################################
# Settings for Differential expression table 1
################################################################################

initial[["DETable1"]] <- new("DETable", ContrastName = "dex: trt vs untrt", RoundDigits = FALSE,
  SignifDigits = 3L, Selected = "SPARCL1", Search = "", SearchColumns = c("",
    "", "", "", "", ""), HiddenColumns = c("baseMean", "lfcSE",
      "stat"), VersionInfo = list(iSEE = structure(list(c(2L, 16L,
        0L)), class = c("package_version", "numeric_version"))),
  PanelId = c(DETable = 1L), PanelHeight = 500L, PanelWidth = 4L,
  SelectionBoxOpen = FALSE, RowSelectionSource = "---", ColumnSelectionSource = "---",
  DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE,
  RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE,
  SelectionHistory = list())

################################################################################
# Settings for MA plot 1
################################################################################

initial[["MAPlot1"]] <- new("MAPlot", ContrastName = "dex: trt vs untrt", FacetRowByRowData = "gene_biotype",
  FacetColumnByRowData = "gene_biotype", ColorByRowData = "gene_id",
  ColorBySampleNameAssay = "counts", ColorByFeatureNameColor = "#FF0000",
  ShapeByRowData = "gene_biotype", SizeByRowData = "entrezid",
  TooltipRowData = c("gene_name", "gene_id", "gene_biotype"
  ), FacetRowBy = "None", FacetColumnBy = "None", ColorBy = "Feature name",
  ColorByDefaultColor = "#000000", ColorByFeatureName = "SPARCL1",
  ColorByFeatureSource = "DETable1", ColorByFeatureDynamicSource = FALSE,
  ColorBySampleName = "SRR1039508", ColorBySampleSource = "---",
  ColorBySampleDynamicSource = FALSE, ShapeBy = "None", SizeBy = "None",
  SelectionAlpha = 0.1, ZoomData = numeric(0), BrushData = list(
    lasso = NULL, closed = TRUE, panelvar1 = NULL, panelvar2 = NULL,
    mapping = list(x = "X", y = "Y"), coord = structure(c(0.192359843265319,
      4.48608873454993, 16.9042262142456, 5.83314093573726,
      0.192359843265319, 4.0415813259785, 1.91060853970459,
      2.20453582056996, 9.95686785339402, 4.0415813259785), dim = c(5L,
        2L))), VisualBoxOpen = FALSE, VisualChoices = "Text",
  ContourAdd = FALSE, ContourColor = "#0000FF", PointSize = 1,
  PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200,
  CustomLabels = FALSE, CustomLabelsText = "TSPAN6", FontSize = 1.5,
  LegendPointSize = 1, LegendPosition = "Bottom", HoverInfo = TRUE,
  LabelCenters = FALSE, LabelCentersBy = "gene_biotype", LabelCentersColor = "#000000",
  VersionInfo = list(iSEE = structure(list(c(2L, 16L, 0L)), class = c("package_version",
    "numeric_version"))), PanelId = c(MAPlot = 1L), PanelHeight = 500L,
  PanelWidth = 4L, SelectionBoxOpen = FALSE, RowSelectionSource = "---",
  ColumnSelectionSource = "---", DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE,
  ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE,
  ColumnSelectionRestrict = FALSE, SelectionHistory = list())

################################################################################
# Settings for Volcano plot 1
################################################################################

initial[["VolcanoPlot1"]] <- new("VolcanoPlot", ContrastName = "dex: trt vs untrt", FacetRowByRowData = "gene_biotype",
  FacetColumnByRowData = "gene_biotype", ColorByRowData = "gene_id",
  ColorBySampleNameAssay = "counts", ColorByFeatureNameColor = "#FF0000",
  ShapeByRowData = "gene_biotype", SizeByRowData = "entrezid",
  TooltipRowData = c("gene_name", "gene_id", "gene_biotype"
  ), FacetRowBy = "None", FacetColumnBy = "None", ColorBy = "Row selection",
  ColorByDefaultColor = "#000000", ColorByFeatureName = "TSPAN6",
  ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE,
  ColorBySampleName = "SRR1039508", ColorBySampleSource = "---",
  ColorBySampleDynamicSource = FALSE, ShapeBy = "None", SizeBy = "None",
  SelectionAlpha = 0.1, ZoomData = numeric(0), BrushData = list(),
  VisualBoxOpen = FALSE, VisualChoices = "Text", ContourAdd = FALSE,
  ContourColor = "#0000FF", PointSize = 1, PointAlpha = 1,
  Downsample = FALSE, DownsampleResolution = 200, CustomLabels = FALSE,
  CustomLabelsText = "TSPAN6", FontSize = 1.5, LegendPointSize = 1,
  LegendPosition = "Bottom", HoverInfo = TRUE, LabelCenters = FALSE,
  LabelCentersBy = "gene_biotype", LabelCentersColor = "#000000",
  VersionInfo = list(iSEE = structure(list(c(2L, 16L, 0L)), class = c("package_version",
    "numeric_version"))), PanelId = c(VolcanoPlot = 1L), PanelHeight = 500L,
  PanelWidth = 4L, SelectionBoxOpen = FALSE, RowSelectionSource = "MAPlot1",
  ColumnSelectionSource = "---", DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE,
  ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE,
  ColumnSelectionRestrict = FALSE, SelectionHistory = list())

app <- iSEE(airway, initial = initial)

if (interactive()) {
  shiny::runApp(app)
}
