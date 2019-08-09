#!/usr/bin/env Rscript

library(ggplot2)
library(forcats)
library(reshape2)
args = commandArgs(trailingOnly = TRUE)

input_comm_csv = args[1]
comm_data = read.csv(input_comm_csv, sep = '\t', header = TRUE)
colnames(comm_data) <- c('first_smpl', 'sec_smpl', 'first', 'snvs', 'findxs', 'second', 'snvs', 'findxs', 'l_200', 'comparison_type', 'SNVs', '%first', '%sec', 'Fingerprint Indexes', '%first', '%sec')
melted_comm_data <- melt(comm_data, 
                  id.vars = c('first_smpl', 'sec_smpl', 'first', 'second', 'comparison_type', 'l_200'), 
                  measure.vars = c('SNVs', 'Fingerprint Indexes'),
                  variable.name = 'cmntype',
                  value.name = 'numshared')
head(melted_comm_data)
output_folder = args[2]

path_arr <- strsplit(input_comm_csv, "/")
len_path <- length(path_arr[[1]])
input_file_name <- path_arr[[1]][len_path]

pdf_name <- paste(output_folder, "/", 
                  gsub('.{4}$', '', path_arr[[1]][len_path]),
                  "corrVcomm2.pdf", sep = "")

graph_title <- paste("from: ", input_file_name, sep = "")
g <- ggplot(melted_comm_data, aes(x = numshared, y = l_200)) + 
     geom_point(color = 'steelblue') +
     geom_smooth(method=lm) + 
     facet_grid(cmntype ~ comparison_type) + 
     ggtitle(graph_title) + 
     labs( x = 'Number of _ Shared Between Sample Pairs',
           y = 'Correlation at L_200')

pdf(pdf_name)
print(g)
dev.off()