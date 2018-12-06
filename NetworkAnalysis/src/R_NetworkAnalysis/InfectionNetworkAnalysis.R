
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

# nsim=100
# folder_shuffled <- 'shuffled/'
# shuffle_models <- c('r0','c0','curveball')

first = strsplit(argument, '/')
second = first[[1]][length(first[[1]])]
third = strsplit(second, 'BipartieInfection_MATRIX_')
third = third[[1]][length(third[[1]])]
name = strsplit(third, '.txt')

name = name[[1]][length(name[[1]])]
out = paste("InfectionNetwork_",name,".png",sep="");

network = BuildBipartiteNetwork(argument)

fortitle = strsplit(name, '_')
title = paste("Infection Network from Simlated Data ","\n ",fortitle[[1]][5]," = ",fortitle[[1]][6],sep="");

png(out, height = 1000, width = 2000, units = "px");
plot_matrix(network, layout = 'nested', method = 'ggplot', title=title, x_title='Bacteria strains', y_title='Virus strains')
dev.off();

# Infomap (Modularity) -----------------------------------------------------------------

# ###############################################
#
# file_prefix <- name
# x <- shuffle_bipartite_matrix(network, shuffle_models, nsim=nsim, burnin=1000, write_files=T, file_prefix=file_prefix, folder=folder_shuffled)
# y <- prob_model_wrapper(network, non_zero_rc_thershold = 0.85, nsim = nsim, write_files = T, file_prefix = file_prefix, folder = folder_shuffled)
# 
# # Read shuffled matrices --------------------------------------------------
# file_prefix <- name
# network_shuffled <- read_shuffled_networks(shuff_methods = shuffle_models, file_prefix = file_prefix, nsim = nsim, folder = folder_shuffled)
# sapply(network_shuffled, dim)
#
# ###############################################

x <- Infomap_wrapper_NoShuffled(Z = network, bipartite_groups = c('Bacteria', 'Virus'), file_prefix = name)

# nameTitle = paste(name," Infection",sep="");
# x$p_value_plot+labs(title=nameTitle)

out = paste("InfectionNetwork-Modularity_",name,".png",sep="");
fortitle = strsplit(name, '_')
title = paste("Infection Network Modularity from Simlated Data ","\n ",fortitle[[1]][5]," = ",fortitle[[1]][6],sep="");

png(out, height = 1000, width = 2000, units = "px");
ggplot_bipartite_modules(Z=network, x$node_data_obs, module_numbers = F, color_tick_labels = T, border=F, text_size = 18, title=title, xlab='Bacteria strain', ylab='Virus strain', weighted=T)
dev.off();


