library(iSEE)
library(fgsea)
library(iSEEpathways)

# Example data ----

set.seed(1)
simulated_data <- simulateExampleData()

pathways_list <- simulated_data[["pathwaysList"]]
features_stat <- simulated_data[["featuresStat"]]
se <- simulated_data[["summarizedexperiment"]]

# fgsea ----

set.seed(42)
fgseaRes <- fgsea(pathways = pathways_list,
  stats    = features_stat,
  minSize  = 15,
  maxSize  = 500)
fgseaRes <- fgseaRes[order(pval), ]

# iSEE ---

se <- embedPathwaysResults(fgseaRes, se, name = "fgsea", class = "fgsea", pathwayType = "simulated",
  pathwaysList = pathways_list, featuresStats = features_stat)

app <- iSEE(se, initial = list(
  PathwaysTable(ResultName="fgsea", PanelWidth = 12L)
))

if (interactive()) {
  shiny::runApp(app)
}
