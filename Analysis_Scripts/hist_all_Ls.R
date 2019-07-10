#!/usr/bin/env Rscript

library(Hmisc)
args = commandArgs(trailingOnly = TRUE)
comparisons <- read.csv(args[1], header = TRUE)

columnNames <- colnames( comparisons )
numberOfColumns <- length( columnNames ) + 1

createColumnHist <- function(numColumns, colNames) {
  currentColumn <- 3
  while (currentColumn != numColumns) {
    currentColumnName <- colNames[currentColumn]
    currentPlotTitle = paste("Frequencies of Correlations for Fingerprint Comparisons with", currentColumnName)
    print(paste("Plotting", currentPlotTitle, "..."))
    h = hist(comparisons[[currentColumnName]],
         main = currentPlotTitle,
         xlab = "Correlation",
         breaks = seq(-0.25, 1, by=0.025),
         col = "#999999",
         plot = TRUE)
    currentColumn = currentColumn + 1
  }
}

createColumnHist(numberOfColumns, colnames( comparisons ))
hist.data.frame(comparisons[, 3:15])



