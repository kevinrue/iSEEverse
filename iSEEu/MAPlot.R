# Making up some results:
se <- SummarizedExperiment(matrix(rnorm(10000), 1000, 10))
rownames(se) <- paste0("GENE_", seq_len(nrow(se)))
rowData(se)$PValue <- runif(nrow(se))
rowData(se)$LogFC <- rnorm(nrow(se))
rowData(se)$AveExpr <- rnorm(nrow(se))

if (interactive()) {
  iSEE(se, initial=list(
    MAPlot(PanelWidth=12L)
  ))
}
