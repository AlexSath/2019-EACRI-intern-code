#!/usr/bin/env Rscript

library(ggplot2)
library(forcats)
library(reshape2)
args = commandArgs(trailingOnly = TRUE)

input_comm_csv = args[1]
comm_data = read.csv(input_comm_csv, sep = '\t', header = TRUE)
colnames(comm_data) <- c('first_smpl', 'sec_smpl', 'first', 'snvs', 'findxs', 'second', 'snvs', 'findxs', 'l_200', 'comparison_type', 'Common SNVs', '%first', '%sec', 'Common Fingerprint Indexes', '%first', '%sec')
melted_comm_data <- melt(comm_data, 
                  id.vars = c('first_smpl', 'sec_smpl', 'first', 'second', 'comparison_type', 'l_200'), 
                  measure.vars = c('Common SNVs', 'Common Fingerprint Indexes'),
                  variable.name = 'cmntype',
                  value.name = 'numshared')
head(melted_comm_data)
split_comm_data <- split.data.frame(melted_comm_data, melted_comm_data$comparison_type)
output_folder = args[2]

path_arr <- strsplit(input_comm_csv, "/")
len_path <- length(path_arr[[1]])
input_file_name <- path_arr[[1]][len_path]

pdf_name <- paste(output_folder, "/", 
                  gsub('.{4}$', '', path_arr[[1]][len_path]),
                  "corrVcomm.pdf", sep = "")

graph_arr = list()
curr_plt = 1
for (data_frame in split_comm_data) {
  frame_name <- deparse(substitute(data_frame))
  graph_title <- paste(frame_name, "data from:\n", input_file_name, sep = "")
  g <- ggplot(data_frame, aes(x = numshared, y = l_200)) + 
       geom_point(aes(shape = cmntype, col = cmntype)) +
       labs()
       ggtitle(graph_title)
  graph_arr[[curr_plt]] <- g
  curr_plt <- curr_plt + 1
}

pdf(pdf_name)
for (graph in graph_arr) {
  print(graph)
}
dev.off()