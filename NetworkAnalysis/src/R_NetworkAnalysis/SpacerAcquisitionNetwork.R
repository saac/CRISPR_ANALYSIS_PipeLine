
# Rscript SpacerAcquisitionNetwork.R Spacers-by-bacteria_mu1e-7_initialDiffDp1_S10P15_R-12499_Time_1750.txt

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
third = strsplit(second, 'Spacers-by-bacteria_')
third = third[[1]][length(third[[1]])]
name = strsplit(third, '.txt')
name = name[[1]][length(name[[1]])]

spacer_by_bact <- load_bipartite_file_1(argument,'\t')
spacer_by_bact <- spacer_by_bact[,str_detect(colnames(spacer_by_bact), 'Sp_')]
spacer_by_bact <- incidence_matrix_to_list(spacer_by_bact)
network <- spacer_by_bact$M
network <- t(network)

out = paste("SpacerAcquisitionNetwork_",name,".png",sep="");

fortitle = strsplit(name, '_')
# title = paste("Spacer acquisition network from Simlated Data ","\n ",fortitle[[1]][5]," = ",fortitle[[1]][6],sep="");
title = paste("Bacteria Aquisition Network from Simulated Data","\n ",fortitle[[1]][5]," = ",fortitle[[1]][6],sep="");

# pdf(out, paper="USr");
# png(out, height = 1600, width = 1000, units = "px",);
# png(out, height = 1600, width = 1380, units = "px");
png(out, height = 1000, width = 2000, units = "px");

plot_matrix(network, layout = 'diagonal', method = 'ggplot', binary_cols = c('white','orange'), title=title, x_title='Spacers', y_title='Bacteria strains')+theme(legend.position = 'none')
dev.off();

