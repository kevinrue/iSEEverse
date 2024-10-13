if (interactive()) {
  iSEE(SummarizedExperiment(), initial=list(
    MarkdownBoard(
      PanelWidth = 12L,
      Content = paste0(c(
        "# Level 1 header",
        "## Level 2 header",
        "**Bold** and *italic*.",
        "[Link](https://isee.github.io/)",
        "* Bullet point\n"),
        collapse = "\n\n"
      )
    )
  ))
}
