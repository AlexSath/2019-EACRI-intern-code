#!/usr/bin/env Rscript
library(ggplot2)
library(forcats)
library(reshape2)
args = commandArgs(trailingOnly = TRUE)

input_path <- args[1]
all_comps <- read.csv(input_path, header = TRUE)
all_comps <- melt(all_comps, 
                  id.vars = c('First_Sample', 'Second_Sample', 'is_pair'), 
                  measure.vars = c('binary','L_5','L_7','L_11','L_13','L_17','L_19','L_20','L_40','L_50','L_80','L_100','L_120','L_200'),
                  variable.name = 'l_num',
                  value.name = 'correlation')

path_arr <- strsplit(input_path, "/")
len_path <- length(path_arr[[1]])
file_name_arr <- strsplit(path_arr[[1]][len_path], "_")
min_depth <- file_name_arr[[1]][2]
min_feq <- file_name_arr[[1]][3]
min_reads <- gsub('.csv', '', file_name_arr[[1]][4])
max_pop_freq <- 100*as.double(gsub('.{4}$', '', file_name_arr[[1]][5]))

output_folder <- args[2]
pdf_name <- paste(output_folder, "/", 
                  gsub('.{4}$', '', path_arr[[1]][len_path]),
                  ".tnp.vlnplot.pdf", sep = "")

g <- ggplot(all_comps, aes(x = l_num, y = correlation)) + 
     geom_violin(alpha = 0.75, aes(fill = all_comps$is_pair), position=position_dodge(1), width = 1) + 
     geom_boxplot(alpha = 0.0, aes(fill = all_comps$is_pair), color = 'black', outlier.alpha = 0.25, position=position_dodge(1), width = 0.2) + 
     ylim(-0.25, 0.75) +
     labs(x = 'Fingerprint L#', y = 'Correlation', fill = 'Pairs Only:') +
     ggtitle(#paste("Correlations of paired and unpaired samples with min_depth of ", min_depth,
             #      ",\n min_freq of ", min_feq, 
             #      ", min_reads of ", min_reads, 
             #      ", max_pop_freq of ", max_pop_freq, "%",
             #      sep = ""
             path_arr[[1]][len_path]) + 
     theme(legend.position = 'top')

pdf(pdf_name)
print(g)
dev.off()

