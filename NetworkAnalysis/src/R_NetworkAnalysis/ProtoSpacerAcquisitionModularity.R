
# Rscript ProtoSpacerAcquisitionNetwork.R Protospacers-by-virus_mu1e-7_initialDiffDp1_S10P15_R-12499_Time_1750.txt

library(tidyverse)
library(magrittr)
library(bipartite)
library(vegan)
library(RColorBrewer)
library(igraph)

source('Functions.R')

# nsim=100
# folder_shuffled <- 'shuffled/'
# shuffle_models <- c('r0','c0','curveball')

args <- commandArgs(trailingOnly = TRUE)
argument = args[1]

first = strsplit(argument, '/')
second = first[[1]][length(first[[1]])]
third = strsplit(second, 'Protospacers-by-virus_')
third = third[[1]][length(third[[1]])]
name = strsplit(third, '.txt')
name = name[[1]][length(name[[1]])]

protospacer_by_virus <- load_bipartite_file_1(argument,'\t')
protospacer_by_virus <- protospacer_by_virus[,str_detect(colnames(protospacer_by_virus), 'Ps_')]
protospacer_by_virus <- incidence_matrix_to_list(protospacer_by_virus)
network <- protospacer_by_virus$M
network <- t(network)


# ###############################################
#
# # Create shuffled networks ------------------------------------------------
# file_prefix <- name
# x <- shuffle_bipartite_matrix(network, shuffle_models, nsim=nsim, burnin=1000, write_files=T, file_prefix=file_prefix, folder=folder_shuffled)
# y <- prob_model_wrapper(network, non_zero_rc_thershold = 0.85, nsim = nsim, write_files = T, file_prefix = file_prefix, folder = folder_shuffled)
# 
# # Read shuffled matrices --------------------------------------------------
# file_prefix <- name
# network_shuffled <- read_shuffled_networks(shuff_methods = shuffle_models, file_prefix = file_prefix, nsim = nsim, folder = folder_shuffled)
# sapply(network_shuffled, dim)


# Infomap -----------------------------------------------------------------

# prt = paste(name,' acquisition:',sep="");
# print(prt)

# x <- Infomap_wrapper(Z = network, shuffled_matrices = network_shuffled, bipartite_groups = c('ProtoSpacer','Virus'), file_prefix = name)
x <- Infomap_wrapper_NoShuffled(Z = network, bipartite_groups = c('ProtoSpacer','Virus'), file_prefix = name)

# x$p_value_table
# nameTitle = paste(name," acquisition",sep="");

out = paste("ProtospacerAcquisitionModularity_",name,".png",sep="");

fortitle = strsplit(name, '_')
title = paste("Virus ProtoSpacer acquisition modularity from Simlated Data ","\n ",fortitle[[1]][5]," = ",fortitle[[1]][6],sep="");
# print(title)

# x$p_value_plot+labs(title=nameTitle)

# pdf(out, paper="USr");

# png(out, height = 600, width = 1000, units = "px");
# png(out, height = 800, width = 1380, units = "px");
png(out, height = 1000, width = 2000, units = "px");

# png(out);
# plot_matrix(network, layout = 'nested', method = 'ggplot', title=title, x_title='Bacteria strains', y_title='Virus strains')
ggplot_bipartite_modules(Z=network, x$node_data_obs, module_numbers = F, color_tick_labels = T, border=F, text_size = 18, title=title, xlab='Protospacer ID', ylab='Virus strain')
dev.off();




