
# Rscript SpacerAcquisitionModularity.R Spacers-by-bacteria_mu1e-7_initialDiffDp1_S10P15_R-12499_Time_1750.txt



library(tidyverse)
library(magrittr)
library(bipartite)
library(vegan)
library(RColorBrewer)
library(igraph)

source('Functions.R')


nsim=100
folder_shuffled <- 'shuffled/'
shuffle_models <- c('r0','c0','curveball')

args <- commandArgs(trailingOnly = TRUE)
argument = args[1]

first = strsplit(argument, '/')
second = first[[1]][length(first[[1]])]
third = strsplit(second, 'Spacers-by-bacteria_')
third = third[[1]][length(third[[1]])]
name = strsplit(third, '.txt')
name = name[[1]][length(name[[1]])]

# print(name)

spacer_by_bact <- load_bipartite_file_1(argument,'\t')
spacer_by_bact <- spacer_by_bact[,str_detect(colnames(spacer_by_bact), 'Sp_')]
spacer_by_bact <- incidence_matrix_to_list(spacer_by_bact)
network <- spacer_by_bact$M
network <- t(network)

# Create shuffled networks ------------------------------------------------
file_prefix <- name
x <- shuffle_bipartite_matrix(network, shuffle_models, nsim=nsim, burnin=1000, write_files=T, file_prefix=file_prefix, folder=folder_shuffled)
y <- prob_model_wrapper(network, non_zero_rc_thershold = 0.85, nsim = nsim, write_files = T, file_prefix = file_prefix, folder = folder_shuffled)

# Read shuffled matrices --------------------------------------------------
file_prefix <- name
network_shuffled <- read_shuffled_networks(shuff_methods = shuffle_models, file_prefix = file_prefix, nsim = nsim, folder = folder_shuffled)
sapply(network_shuffled, dim)

# print(network_shuffled)

# Infomap -----------------------------------------------------------------

prt = paste(name,' acquisition:',sep="");
print(prt)
x <- Infomap_wrapper(Z = network, shuffled_matrices = network_shuffled, bipartite_groups = c('Spacer','Bacteria'), file_prefix = name)
x$p_value_table

nameTitle = paste(name," acquisition",sep="");
out = paste("SpacerAcquisitionModularity_",name,".png",sep="");

fortitle = strsplit(name, '_')
title = paste("Bacteria Spacer acquisition modularity from Simlated Data ","\n ",fortitle[[1]][5]," = ",fortitle[[1]][6],sep="");
# print(title)

x$p_value_plot+labs(title=nameTitle)

# pdf(out, paper="USr");
png(out, height = 600, width = 1000, units = "px",);
# png(out);
# plot_matrix(network, layout = 'nested', method = 'ggplot', title=title, x_title='Bacteria strains', y_title='Virus strains')
ggplot_bipartite_modules(Z=network, x$node_data_obs, module_numbers = F, color_tick_labels = T, border=F, text_size = 18, title=title, xlab='Spacer ID', ylab='Bacteria strain')
dev.off();



# # ###########################################################
# # 
# # Create shuffled networks ------------------------------------------------
# ## @knitr shuffle_networks
# file_prefix <- 'mu3e7_Dp1_S10P15_1750'
# x <- shuffle_bipartite_matrix(mu3e7_Dp1_S10P15_1750, shuffle_models, nsim=nsim, burnin=1000, write_files=T, file_prefix=file_prefix, folder=folder_shuffled)
# y <- prob_model_wrapper(mu3e7_Dp1_S10P15_1750, non_zero_rc_thershold = 0.85, nsim = nsim, write_files = T, file_prefix = file_prefix, folder = folder_shuffled)
# 
# # Read shuffled matrices --------------------------------------------------
# ## @knitr get_shuffled_matrices
# file_prefix <- 'mu3e7_Dp1_S10P15_1750'
# mu3e7_Dp1_S10P15_1750_shuffled <- read_shuffled_networks(shuff_methods = shuffle_models, file_prefix = file_prefix, nsim = nsim, folder = folder_shuffled)
# sapply(mu3e7_Dp1_S10P15_1750_shuffled, dim)
# 
# # # Infomap -----------------------------------------------------------------
# #
# ## @knitr Infomap_mu3e7_Dp1_S10P15_1750_A
# print('mu3e7_Dp1_S10P15_1750 acquisition:')
# x <- Infomap_wrapper(Z = mu3e7_Dp1_S10P15_1750, shuffled_matrices = mu3e7_Dp1_S10P15_1750_shuffled, bipartite_groups = c('spacer','archaea'), file_prefix = 'mu3e7_Dp1_S10P15_1750')
# x$p_value_table
# x$p_value_plot+labs(title='mu3e7_Dp1_S10P15_1750 acquisition')
# ggplot_bipartite_modules(Z=mu3e7_Dp1_S10P15_1750, x$node_data_obs, module_numbers = F, color_tick_labels = T, border=F, text_size = 18, title='Spacer acquisition modularity t = 1750 hr', xlab='Spacer ID', ylab='Bacteria strain')
# # 
# # ###########################################################
