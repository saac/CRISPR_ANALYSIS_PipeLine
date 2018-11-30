
# Rscript ProtoSpacerAcquisitionNetwork.R Protospacers-by-virus_mu1e-7_initialDiffDp1_S10P15_R-12499_Time_1750.txt

library(tidyverse)
library(magrittr)
library(bipartite)
library(vegan)
library(RColorBrewer)
library(igraph)

source('Functions.R')

args <- commandArgs(trailingOnly = TRUE)
argument = args[1]

first = strsplit(argument, '/')
second = first[[1]][length(first[[1]])]
third = strsplit(second, 'Protospacers-by-virus_')
third = third[[1]][length(third[[1]])]
name = strsplit(third, '.txt')
name = name[[1]][length(name[[1]])]

# print(name)


protospacer_by_virus <- load_bipartite_file_1(argument,'\t')
protospacer_by_virus <- protospacer_by_virus[,str_detect(colnames(protospacer_by_virus), 'Ps_')]
protospacer_by_virus <- incidence_matrix_to_list(protospacer_by_virus)
network <- protospacer_by_virus$M
network <- t(network)

out = paste("ProtoSpacerAcquisitionNetwork_",name,".png",sep="");

fortitle = strsplit(name, '_')
# title = paste("Spacer acquisition network from Simlated Data ","\n ",fortitle[[1]][5]," = ",fortitle[[1]][6],sep="");
title = paste("PsAN from Simlated Data ","\n ",fortitle[[1]][5]," = ",fortitle[[1]][6],sep="");
# print(title)

# pdf(out, paper="USr");

# png(out, height = 600, width = 1000, units = "px");
# png(out, height = 768, width = 1024, units = "px",res = 600);
# png(out, height = 800, width = 1380, units = "px",res = 600);
# png(out, height = 800, width = 1380, units = "px");
png(out, height = 1000, width = 2000, units = "px");

# png(out);
plot_matrix(network, layout = 'nested', method = 'ggplot', binary_cols = c('white','darkgreen'), title=title, x_title='Proto Spacers', y_title='Virus strains')+theme(legend.position = 'none')
dev.off();



# # ###########################################################
# 
# protospacer_by_virus <- load_bipartite_file_1('mu3e-7_lE-rdw/Protospacers-by-virus_mu1e-7_initialDiffDp1_S10P15_R-12499_Time_1750.txt','\t')
# protospacer_by_virus <- protospacer_by_virus[,str_detect(colnames(protospacer_by_virus), 'Ps_')]
# protospacer_by_virus <- incidence_matrix_to_list(protospacer_by_virus)
# mu3e7_Dp1_S10P15_1750_virus <- protospacer_by_virus$M
# mu3e7_Dp1_S10P15_1750_virus <- t(mu3e7_Dp1_S10P15_1750_virus)
# plot_matrix(mu3e7_Dp1_S10P15_1750_virus, layout = 'nested', method = 'ggplot', binary_cols = c('white','darkgreen'), title='Virus-protospacer network from Simlated Data  t = 1750 hr', x_title='Proto Spacers', y_title='Virus strains')+theme(legend.position = 'none')
# 
# # ###########################################################
