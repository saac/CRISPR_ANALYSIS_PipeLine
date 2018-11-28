
# Rscript ImmunityNetwork.R Bipartite_MATRIX_mu1e-7_initialDiffDp1_S10P15_R-12499_Time_1750.txt

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
third = strsplit(second, 'Bipartite_MATRIX_')
third = third[[1]][length(third[[1]])]
name = strsplit(third, '.txt')
name = name[[1]][length(name[[1]])]

# print(name)


Bipartite_MATRIX <- load_bipartite_file_1(argument,'\t')
Bipartite_MATRIX <- Bipartite_MATRIX[,-which(colnames(Bipartite_MATRIX)=='X.1')]
network <-  t(Bipartite_MATRIX)



out = paste("ImmunityNetwork_",name,".png",sep="");

fortitle = strsplit(name, '_')
# title = paste("Spacer acquisition network from Simlated Data ","\n ",fortitle[[1]][5]," = ",fortitle[[1]][6],sep="");
title = paste("Immunity Network from Simlated Data ","\n ",fortitle[[1]][5]," = ",fortitle[[1]][6],sep="");
# print(title)

# pdf(out, paper="USr");
png(out, height = 600, width = 1000, units = "px",);
# png(out);
plot_matrix(network, layout = 'nested', method = 'ggplot', title=title, x_title='Bacteria strains', y_title='Virus strains')
dev.off();



# # ###########################################################
# 
# Bipartite_MATRIX <- load_bipartite_file_1('mu3e-7_lE-rdw/Bipartite_MATRIX_mu1e-7_initialDiffDp1_S10P15_R-12499_Time_1500.txt','\t')
# Bipartite_MATRIX <- Bipartite_MATRIX[,-which(colnames(Bipartite_MATRIX)=='X.1')]
# mu3e7_Dp1_S10P15_1500_I <- t(Bipartite_MATRIX)
# plot_matrix(mu3e7_Dp1_S10P15_1500_I, layout = 'nested', method = 'ggplot', title='Immunity network from Simlated Data  t = 1500 hr', x_title='Bacteria strains', y_title='Virus strains')
# 
# # ###########################################################
