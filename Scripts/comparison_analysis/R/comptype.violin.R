#!/usr/bin/env Rscript
library(ggplot2)
library(forcats)
library(reshape2)
args = commandArgs(trailingOnly = TRUE)

input_path <- args[1]
all_comps <- read.csv(input_path, header = TRUE)
all_comps <- melt(all_comps, 
                  id.vars = c('First_Sample', 'Second_Sample', 'comp_type'), 
                  measure.vars = c('binary','L_5','L_13','L_20','L_50','L_100','L_200'),
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
                  ".cmptp.vln.pdf", sep = "")

g <- ggplot(all_comps, aes(x = l_num, y = correlation)) + 
  geom_violin(alpha = 0.75, aes(fill = all_comps$comp_type), position=position_dodge(1), width = 1) + 
  geom_boxplot(alpha = 0.0, aes(fill = all_comps$comp_type), color = 'black', outlier.alpha = 0.25, position=position_dodge(1), width = 0.2) + 
  ylim(-0.25, 0.75) +
  labs(x = 'Fingerprint L#', y = 'Correlation', fill = 'Comparison Type:') +
  ggtitle(paste(path_arr[[1]][len_path])) + 
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),
        legend.position = 'top') + 
  guides(fill=guide_legend(nrow=2))

pdf(pdf_name)
print(g)
dev.off()

