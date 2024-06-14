# Source: <https://bioconductor.org/packages/release/bioc/vignettes/iSEEpathways/inst/doc/integration.R>
# (Adapted!)

library(iSEEpathways)

## ----message=FALSE, warning=FALSE-------------------------------------------------------------------------------------
library("airway")
data("airway")
airway$dex <- relevel(airway$dex, "untrt")

## ----message=FALSE, warning=FALSE-------------------------------------------------------------------------------------
library("org.Hs.eg.db")
library("scater")
rowData(airway)[["ENSEMBL"]] <- rownames(airway)
rowData(airway)[["SYMBOL"]] <- mapIds(org.Hs.eg.db, rownames(airway), "SYMBOL", "ENSEMBL")
rowData(airway)[["uniquifyFeatureNames"]] <- uniquifyFeatureNames(
  ID = rowData(airway)[["ENSEMBL"]],
  names = rowData(airway)[["SYMBOL"]]
)
rownames(airway) <- rowData(airway)[["uniquifyFeatureNames"]]

## ---------------------------------------------------------------------------------------------------------------------
library("scuttle")
airway <- logNormCounts(airway)

## ----message=FALSE, warning=FALSE-------------------------------------------------------------------------------------
library("edgeR")

counts <- assay(airway, "counts")
design <- model.matrix(~ 0 + dex + cell, data = colData(airway))

keep <- filterByExpr(counts, design)
v <- voom(counts[keep,], design, plot=FALSE)
fit <- lmFit(v, design)
contr <- makeContrasts("dextrt - dexuntrt", levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
res_limma <- topTable(tmp, sort.by = "P", n = Inf)
head(res_limma)

## ---------------------------------------------------------------------------------------------------------------------
library("iSEEde")
airway <- iSEEde::embedContrastResults(res_limma, airway, name = "Limma-Voom", class = "limma")
rowData(airway)

## ---------------------------------------------------------------------------------------------------------------------
library("org.Hs.eg.db")
pathways <- AnnotationDbi::select(org.Hs.eg.db, keys(org.Hs.eg.db, "ENSEMBL"), c("GOALL"), keytype = "ENSEMBL")
pathways <- subset(pathways, ONTOLOGYALL == "BP")
pathways <- unique(pathways[, c("ENSEMBL", "GOALL")])
pathways <- merge(pathways, rowData(airway)[, c("ENSEMBL", "uniquifyFeatureNames")])
pathways <- split(pathways$uniquifyFeatureNames, pathways$GOALL)

## ---------------------------------------------------------------------------------------------------------------------
map_GO <- function(pathway_id, se) {
  pathway_ensembl <- mapIds(org.Hs.eg.db, pathway_id, "ENSEMBL", keytype = "GOALL", multiVals = "CharacterList")[[pathway_id]]
  pathway_rownames <- rownames(se)[rowData(se)[["gene_id"]] %in% pathway_ensembl]
  pathway_rownames
}
airway <- registerAppOptions(airway, Pathways.map.functions = list(GO = map_GO))

## ---------------------------------------------------------------------------------------------------------------------
library("fgsea")
set.seed(42)
stats <- na.omit(log2FoldChange(contrastResults(airway, "Limma-Voom")))
fgseaRes <- fgsea(pathways = pathways,
  stats    = stats,
  minSize  = 15,
  maxSize  = 500)
head(fgseaRes[order(pval), ])

library(dplyr)
fgseaRes <- fgseaRes %>%
  mutate(across(where(is.numeric), function(x) {
    as.numeric(format(x, digits = 2))
  }))

## ---------------------------------------------------------------------------------------------------------------------
library("iSEEpathways")
fgseaRes <- fgseaRes[order(pval), ]
airway <- embedPathwaysResults(
  fgseaRes, airway, name = "fgsea (p-value)", class = "fgsea",
  pathwayType = "GO", pathwaysList = pathways, featuresStats = stats)
airway

## ----warning=FALSE----------------------------------------------------------------------------------------------------
stats <- na.omit(
  log2FoldChange(contrastResults(airway, "Limma-Voom")) *
    -log10(pValue(contrastResults(airway, "Limma-Voom")))
)
set.seed(42)
fgseaRes <- fgsea(pathways = pathways,
  stats    = na.omit(stats),
  minSize  = 15,
  maxSize  = 500)
fgseaRes <- fgseaRes[order(pval), ]

library(dplyr)
fgseaRes <- fgseaRes %>%
  mutate(across(where(is.numeric), function(x) {
    as.numeric(format(x, digits = 2))
  }))

airway <- embedPathwaysResults(
  fgseaRes, airway, name = "fgsea (p-value & fold-change)", class = "fgsea",
  pathwayType = "GO", pathwaysList = pathways, featuresStats = stats)
airway

## ----warning=FALSE----------------------------------------------------------------------------------------------------
library("GO.db")
library("shiny")
library("iSEE")
go_details <- function(x) {
  info <- AnnotationDbi::select(GO.db, x, c("TERM", "ONTOLOGY", "DEFINITION"), "GOID")
  html <- list(p(strong(info$GOID), ":", info$TERM, paste0("(", info$ONTOLOGY, ")")))
  if (!is.na(info$DEFINITION)) {
    html <- append(html, list(p(info$DEFINITION)))
  }
  tagList(html)
}
airway <- registerAppOptions(airway, PathwaysTable.select.details = go_details)

initial <- list()

################################################################################
# Settings for Pathways Analysis Table 1
################################################################################

initial[["PathwaysTable1"]] <- new("PathwaysTable", ResultName = "fgsea (p-value)", Selected = "GO:0046324",
  Search = "", SearchColumns = c("", "", "", "", "", "", ""
  ), HiddenColumns = character(0), VersionInfo = list(iSEE = structure(list(
    c(2L, 16L, 0L)), class = c("package_version", "numeric_version"
    ))), PanelId = c(PathwaysTable = 1L), PanelHeight = 500L,
  PanelWidth = 2L, SelectionBoxOpen = FALSE, RowSelectionSource = "---",
  ColumnSelectionSource = "---", DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE,
  ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE,
  ColumnSelectionRestrict = FALSE, SelectionHistory = list())

################################################################################
# Settings for GSEA enrichment plot 1
################################################################################

initial[["FgseaEnrichmentPlot1"]] <- new("FgseaEnrichmentPlot", ResultName = "fgsea (p-value)", PathwayId = "GO:0046323",
  BrushData = list(xmin = -280.88574988977, xmax = 3554.8142539503,
    ymin = -0.075085779852016, ymax = 0.68357088641439, coords_css = list(
      xmin = 57.2244567871094, xmax = 203.224456787109,
      ymin = 22.780647920138, ymax = 467.291584401186),
    coords_img = list(xmin = 62.9469024658203, xmax = 223.54690246582,
      ymin = 25.0587127121518, ymax = 514.020742841304),
    img_css_ratio = list(x = 1.1, y = 1.1), mapping = list(
      x = "rank", y = "ES"), domain = list(left = -843.05,
        right = 17704.05, bottom = -0.0750857798520159, top = 0.683570886414392),
    range = list(left = 39.4091965771124, right = 815.972602739726,
      bottom = 514.020742841304, top = 25.0587127121518),
    log = list(x = NULL, y = NULL), direction = "x", brushId = "FgseaEnrichmentPlot1_Brush",
    outputId = "FgseaEnrichmentPlot1"), VersionInfo = list(
      iSEE = structure(list(c(2L, 16L, 0L)), class = c("package_version",
        "numeric_version"))), PanelId = c(FgseaEnrichmentPlot = 1L),
  PanelHeight = 500L, PanelWidth = 4L, SelectionBoxOpen = FALSE,
  RowSelectionSource = "---", ColumnSelectionSource = "---",
  DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE,
  RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE,
  SelectionHistory = list())

################################################################################
# Settings for Volcano plot 1
################################################################################

initial[["VolcanoPlot1"]] <- new("VolcanoPlot", ContrastName = "Limma-Voom", FacetRowByRowData = "gene_biotype",
  FacetColumnByRowData = "gene_biotype", ColorByRowData = "gene_id",
  ColorBySampleNameAssay = "logcounts", ColorByFeatureNameColor = "#FF0000",
  ShapeByRowData = "gene_biotype", SizeByRowData = "entrezid",
  TooltipRowData = c("gene_name", "gene_id", "gene_biotype"
  ), FacetRowBy = "None", FacetColumnBy = "None", ColorBy = "Row selection",
  ColorByDefaultColor = "#000000", ColorByFeatureName = "TSPAN6",
  ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE,
  ColorBySampleName = "SRR1039508", ColorBySampleSource = "---",
  ColorBySampleDynamicSource = FALSE, ShapeBy = "None", SizeBy = "None",
  SelectionAlpha = 0.1, ZoomData = numeric(0), BrushData = list(),
  VisualBoxOpen = FALSE, VisualChoices = "Color", ContourAdd = FALSE,
  ContourColor = "#0000FF", PointSize = 1, PointAlpha = 1,
  Downsample = FALSE, DownsampleResolution = 200, CustomLabels = FALSE,
  CustomLabelsText = "TSPAN6", FontSize = 1.5, LegendPointSize = 1,
  LegendPosition = "Bottom", HoverInfo = TRUE, LabelCenters = FALSE,
  LabelCentersBy = "gene_biotype", LabelCentersColor = "#000000",
  VersionInfo = list(iSEE = structure(list(c(2L, 16L, 0L)), class = c("package_version",
    "numeric_version"))), PanelId = c(VolcanoPlot = 1L), PanelHeight = 500L,
  PanelWidth = 3L, SelectionBoxOpen = FALSE, RowSelectionSource = "FgseaEnrichmentPlot1",
  ColumnSelectionSource = "---", DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE,
  ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE,
  ColumnSelectionRestrict = FALSE, SelectionHistory = list())

################################################################################
# Settings for Complex heatmap 1
################################################################################

initial[["ComplexHeatmapPlot1"]] <- new("ComplexHeatmapPlot", Assay = "logcounts", CustomRows = FALSE,
  CustomRowsText = "TSPAN6", ClusterRows = TRUE, ClusterRowsDistance = "euclidean",
  ClusterRowsMethod = "ward.D2", DataBoxOpen = FALSE, VisualChoices = "Annotations",
  ColumnData = "dex", RowData = character(0), CustomBounds = FALSE,
  LowerBound = NA_real_, UpperBound = NA_real_, AssayCenterRows = TRUE,
  AssayScaleRows = FALSE, DivergentColormap = "purple < black < yellow",
  ShowDimNames = "Rows", LegendPosition = "Bottom", LegendDirection = "Horizontal",
  VisualBoxOpen = FALSE, NamesRowFontSize = 12, NamesColumnFontSize = 10,
  ShowColumnSelection = FALSE, OrderColumnSelection = TRUE,
  VersionInfo = list(iSEE = structure(list(c(2L, 16L, 0L)), class = c("package_version",
    "numeric_version"))), PanelId = c(ComplexHeatmapPlot = 1L),
  PanelHeight = 700L, PanelWidth = 3L, SelectionBoxOpen = FALSE,
  RowSelectionSource = "FgseaEnrichmentPlot1", ColumnSelectionSource = "---",
  RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE,
  RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE,
  SelectionHistory = list())

## ----"start", message=FALSE, warning=FALSE----------------------------------------------------------------------------
app <- iSEE(airway, initial = initial)

if (interactive()) {
  shiny::runApp(app)
}

