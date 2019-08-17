#!/usr/bin/env Rscript

library(ggplot2)
library(forcats)
library(reshape2)
args = commandArgs(trailingOnly = TRUE)

input_comm_csv = args[1]
comm_data = read.csv(input_comm_csv, sep = '\t', header = TRUE)
colnames(comm_data) <- c('first_smpl', 'sec_smpl', 'first', 'snvs1', 'findxs1', 'second', 'snvs2', 'findxs2', 'l_200', 'comparison_type', 'SNVs', '%firstS', '%secS', 'FingerprintIndexes', '%firstF', '%secF')
#melted_comm_data <- melt(comm_data, 
#                  id.vars = c('first_smpl', 'sec_smpl', 'first', 'second', 'comparison_type', 'l_200'), 
#                  measure.vars = c('SNVs', 'Fingerprint Indexes'),
#                  variable.name = 'cmntype',
#                  value.name = 'numshared')
head(comm_data)
output_folder = args[2]

path_arr <- strsplit(input_comm_csv, "/")
len_path <- length(path_arr[[1]])
input_file_name <- path_arr[[1]][len_path]

pdf_name <- paste(output_folder, "/", 
                  gsub('.{4}$', '', path_arr[[1]][len_path]),
                  ".snvsVSfindxs.pdf", sep = "")

graph_title <- paste("from: ", input_file_name, sep = "")
g <- ggplot(comm_data, aes(x = SNVs, y = FingerprintIndexes)) + 
     geom_point(color = 'steelblue') +
     geom_smooth(method=lm) + 
     #facet_wrap(comparison_type) + 
     ggtitle(graph_title) + 
     labs( x = 'Number of SNVs Shared Between Samples',
           y = 'Number of Fingerprint Indexes Shared Between Samples')

pdf(pdf_name, width = 6, height = 8)
print(g)
dev.off()