# Making up some results:
se <- SummarizedExperiment(matrix(rnorm(10000), 1000, 10))
rownames(se) <- paste0("GENE_", seq_len(nrow(se)))
rowData(se)$PValue1 <- runif(nrow(se))
rowData(se)$LogFC1 <- rnorm(nrow(se))
rowData(se)$PValue2 <- runif(nrow(se))
rowData(se)$LogFC2 <- rnorm(nrow(se))

if (interactive()) {
  iSEE(se, initial=list(
    LogFCLogFCPlot(
      PanelWidth = 12L,
      XAxisRowData="LogFC1", YAxis="LogFC2",
      XPValueField="PValue1", YPValueField="PValue2"
    )
  ))
}
