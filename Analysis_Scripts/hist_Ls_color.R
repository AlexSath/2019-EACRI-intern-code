#!/usr/bin/env Rscript
library(ggplot2)
args = commandArgs(trailingOnly = TRUE)

#From MrFlick on StackExchange (https://stackoverflow.com/questions/28777626/how-do-i-combine-aes-and-aes-string-options)
`+.uneval` <- function(a,b) {
  `class<-`(modifyList(a,b), "uneval")
}

all_comps <- read.csv(args[1], header = TRUE)
comp_pairs <- read.csv(args[2], header = TRUE)
bnw <- as.double(args[4])
output_folder <- args[3]
#all_comps <- read.csv("/Users/asathler/Documents/RAnalysis/v050_filtered_50_10_2.csv", header = TRUE)
#comp_pairs <- read.csv("/Users/asathler/Documents/RAnalysis/v050_filtered_50_10_2_tnp.csv", header = TRUE)
#bnw <- 0.00625
#output_folder <- "/Users/asathler/Documents/RAnalysis"
pdf_name <- paste(output_folder, "/Rplots.pdf", sep = "")

col_names <- colnames(all_comps)
curr_col <- 1
bin_col <- 0
for (column in col_names) {
  if (column == 'binary') {
    bin_col <- curr_col
  }
  curr_col <- curr_col + 1
}

final_col <- bin_col + 13
col_names <- col_names[bin_col:final_col]
print(col_names)

plot_list = list()
current_col = 1
for (col in col_names) {
  chart_title <- paste("Correlation for", col)
  plt <- ggplot(all_comps, aes_string(col)) +
    geom_histogram(aes(fill = 'Unpaired'), 
                   alpha = 0.5, 
                   binwidth = bnw) +
    geom_histogram(data = comp_pairs, 
                   aes_string(col) + aes(y=..count..*1000, fill = 'Paired'), 
                   alpha = 0.5, 
                   binwidth = bnw) +
    scale_y_continuous(sec.axis = sec_axis(~./1000, name = 'Paired Count (Red)')) +
    xlim(-0.1, 0.9) +
    scale_colour_manual(values = c('red', 'blue')) +
    labs( x = paste("Correlation for ", col), 
          y = "Unpaired Count (Blue)", 
          fill = "Legend" ) +
    theme(legend.position = c(0.88, 0.9))
  plot_list[[current_col]] <- plt
  current_col <- current_col + 1
}

pdf(pdf_name)
current_plot = 1
while (current_plot <= length(plot_list)) {
  print(plot_list[[current_plot]])
  current_plot <- current_plot + 1
}
dev.off()

