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

# launch the app itself ----

app <- iSEE(sce, initial = list(
  ColumnDataPlot(
    PanelWidth = 8L,
    YAxis = "NREADS",
    XAxis = "Column data",
    XAxisColumnData = "driver_1_s",
    ColorBy = "Column data",
    ColorByColumnData = "driver_1_s",
    FacetColumnBy = "Column data",
    FacetColumnByColData = "Core.Type"
    )
))

if (interactive()) {
  shiny::runApp(app, port=1234)
}


initial[["ColumnDataPlot1"]] <- new("ColumnDataPlot", XAxis = "Column data", YAxis = "NREADS",
  XAxisColumnData = "driver_1_s", FacetRowByColData = "driver_1_s",
  FacetColumnByColData = "Core.Type", ColorByColumnData = "driver_1_s",
  ColorByFeatureNameAssay = "logcounts", ColorBySampleNameColor = "#FF0000",
  ShapeByColumnData = "driver_1_s", SizeByColumnData = "NREADS",
  TooltipColumnData = character(0), FacetRowBy = "None", FacetColumnBy = "Column data",
  ColorBy = "Column data", ColorByDefaultColor = "#000000",
  ColorByFeatureName = "0610007P14Rik", ColorByFeatureSource = "---",
  ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "SRR2140028",
  ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE,
  ShapeBy = "None", SizeBy = "None", SelectionAlpha = 0.1,
  ZoomData = numeric(0), BrushData = list(), VisualBoxOpen = FALSE,
  VisualChoices = "Facet", ContourAdd = FALSE, ContourColor = "#0000FF",
  FixAspectRatio = FALSE, ViolinAdd = TRUE, PointSize = 1,
  PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200,
  CustomLabels = FALSE, CustomLabelsText = "SRR2140028", FontSize = 1,
  LegendPointSize = 1, LegendPosition = "Bottom", HoverInfo = TRUE,
  LabelCenters = FALSE, LabelCentersBy = "driver_1_s", LabelCentersColor = "#000000",
  VersionInfo = list(iSEE = structure(list(c(2L, 17L, 4L)), class = c("package_version",
    "numeric_version"))), PanelId = c(ColumnDataPlot = 1L), PanelHeight = 500L,
  PanelWidth = 8L, SelectionBoxOpen = FALSE, RowSelectionSource = "---",
  ColumnSelectionSource = "---", DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE,
  ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE,
  ColumnSelectionRestrict = FALSE, SelectionHistory = list())
