library(tidyverse)
library(magrittr)
library(bipartite)
library(vegan)
library(RColorBrewer)
library(igraph)

source('functions.R')

spacer_by_bact <- load_bipartite_file_1('mu3e-7_lE-rdw/Spacers-by-bacteria_mu1e-7_initialDiffDp1_S10P15_R-12499_Time_1500.txt','\t')
spacer_by_bact <- spacer_by_bact[,str_detect(colnames(spacer_by_bact), 'Sp_')]
spacer_by_bact <- incidence_matrix_to_list(spacer_by_bact)
mu3e7_Dp1_S10P15_1500 <- spacer_by_bact$M

# dim(mu3e7_Dp1_S10P15_1500)

mu3e7_Dp1_S10P15_1500 <- t(mu3e7_Dp1_S10P15_1500)

# dim(mu3e7_Dp1_S10P15_1500)

# plot_matrix(mu3e7_Dp1_S10P15_1500, layout = 'diagonal', method = 'ggplot', binary_cols = c('white','orange'), title='Spacer acquisition network from Simlated Data  t = 1500 hr', y_title='Bacteria strains', x_title='Spacers')+theme(legend.position = 'none')
plot_matrix(mu3e7_Dp1_S10P15_1500, layout = 'diagonal', method = 'ggplot', binary_cols = c('white','orange'), title='Spacer acquisition network from Simlated Data  t = 1500 hr', x_title='Spacers', y_title='Bacteria strains')+theme(legend.position = 'none')


###########################################################

spacer_by_bact <- load_bipartite_file_1('mu3e-7_lE-rdw/Spacers-by-bacteria_mu1e-7_initialDiffDp1_S10P15_R-12499_Time_1750.txt','\t')
spacer_by_bact <- spacer_by_bact[,str_detect(colnames(spacer_by_bact), 'Sp_')]
spacer_by_bact <- incidence_matrix_to_list(spacer_by_bact)
mu3e7_Dp1_S10P15_1750 <- spacer_by_bact$M
mu3e7_Dp1_S10P15_1750 <- t(mu3e7_Dp1_S10P15_1750)
plot_matrix(mu3e7_Dp1_S10P15_1750, layout = 'diagonal', method = 'ggplot', binary_cols = c('white','orange'), title='Spacer acquisition network from Simlated Data  t = 1750 hr', x_title='Spacers', y_title='Bacteria strains')+theme(legend.position = 'none')

###########################################################

spacer_by_bact <- load_bipartite_file_1('mu3e-7_lE-rdw/Spacers-by-bacteria_mu1e-7_initialDiffDp1_S10P15_R-12499_Time_2000.txt','\t')
spacer_by_bact <- spacer_by_bact[,str_detect(colnames(spacer_by_bact), 'Sp_')]
spacer_by_bact <- incidence_matrix_to_list(spacer_by_bact)
mu3e7_Dp1_S10P15_2000 <- spacer_by_bact$M
mu3e7_Dp1_S10P15_2000 <- t(mu3e7_Dp1_S10P15_2000)
plot_matrix(mu3e7_Dp1_S10P15_2000, layout = 'diagonal', method = 'ggplot', binary_cols = c('white','orange'), title='Spacer acquisition network from Simlated Data  t = 2000 hr', x_title='Spacers', y_title='Bacteria strains')+theme(legend.position = 'none')

###########################################################
###########################################################

protospacer_by_virus <- load_bipartite_file_1('mu3e-7_lE-rdw/Protospacers-by-virus_mu1e-7_initialDiffDp1_S10P15_R-12499_Time_1500.txt','\t')
protospacer_by_virus <- protospacer_by_virus[,str_detect(colnames(protospacer_by_virus), 'Ps_')]
protospacer_by_virus <- incidence_matrix_to_list(protospacer_by_virus)
mu3e7_Dp1_S10P15_1500_virus <- protospacer_by_virus$M
mu3e7_Dp1_S10P15_1500_virus <- t(mu3e7_Dp1_S10P15_1500_virus)
# plot_matrix(mu3e7_Dp1_S10P15_1500_virus, layout = 'diagonal', method = 'ggplot', binary_cols = c('white','orange'), title='Spacer acquisition network from Simlated Data  t = 2000 hr', x_title='Spacers', y_title='Bacteria strains')+theme(legend.position = 'none')
plot_matrix(mu3e7_Dp1_S10P15_1500_virus, layout = 'nested', method = 'ggplot', binary_cols = c('white','darkgreen'), title='Virus-protospacer network from Simlated Data  t = 1500 hr', x_title='Proto Spacers', y_title='Virus strains')+theme(legend.position = 'none')


###########################################################

protospacer_by_virus <- load_bipartite_file_1('mu3e-7_lE-rdw/Protospacers-by-virus_mu1e-7_initialDiffDp1_S10P15_R-12499_Time_1750.txt','\t')
protospacer_by_virus <- protospacer_by_virus[,str_detect(colnames(protospacer_by_virus), 'Ps_')]
protospacer_by_virus <- incidence_matrix_to_list(protospacer_by_virus)
mu3e7_Dp1_S10P15_1750_virus <- protospacer_by_virus$M
mu3e7_Dp1_S10P15_1750_virus <- t(mu3e7_Dp1_S10P15_1750_virus)
plot_matrix(mu3e7_Dp1_S10P15_1750_virus, layout = 'nested', method = 'ggplot', binary_cols = c('white','darkgreen'), title='Virus-protospacer network from Simlated Data  t = 1750 hr', x_title='Proto Spacers', y_title='Virus strains')+theme(legend.position = 'none')

###########################################################

protospacer_by_virus <- load_bipartite_file_1('mu3e-7_lE-rdw/Protospacers-by-virus_mu1e-7_initialDiffDp1_S10P15_R-12499_Time_2000.txt','\t')
protospacer_by_virus <- protospacer_by_virus[,str_detect(colnames(protospacer_by_virus), 'Ps_')]
protospacer_by_virus <- incidence_matrix_to_list(protospacer_by_virus)
mu3e7_Dp1_S10P15_2000_virus <- protospacer_by_virus$M
mu3e7_Dp1_S10P15_2000_virus <- t(mu3e7_Dp1_S10P15_2000_virus)
plot_matrix(mu3e7_Dp1_S10P15_2000_virus, layout = 'nested', method = 'ggplot', binary_cols = c('white','darkgreen'), title='Virus-protospacer network from Simlated Data  t = 2000 hr', x_title='Proto Spacers', y_title='Virus strains')+theme(legend.position = 'none')



###########################################################
######################################################################################################################
###########################################################


Bipartite_MATRIX <- load_bipartite_file_1('mu3e-7_lE-rdw/Bipartite_MATRIX_mu1e-7_initialDiffDp1_S10P15_R-12499_Time_1500.txt','\t')
Bipartite_MATRIX <- Bipartite_MATRIX[,-which(colnames(Bipartite_MATRIX)=='X.1')]
mu3e7_Dp1_S10P15_1500_I <- t(Bipartite_MATRIX)
plot_matrix(mu3e7_Dp1_S10P15_1500_I, layout = 'nested', method = 'ggplot', title='Immunity network from Simlated Data  t = 1500 hr', x_title='Bacteria strains', y_title='Virus strains')

###########################################################

Bipartite_MATRIX <- load_bipartite_file_1('mu3e-7_lE-rdw/Bipartite_MATRIX_mu1e-7_initialDiffDp1_S10P15_R-12499_Time_1750.txt','\t')
Bipartite_MATRIX <- Bipartite_MATRIX[,-which(colnames(Bipartite_MATRIX)=='X.1')]
mu3e7_Dp1_S10P15_1750_I <- t(Bipartite_MATRIX)
plot_matrix(mu3e7_Dp1_S10P15_1750_I, layout = 'nested', method = 'ggplot', title='Immunity network from Simlated Data  t = 1750 hr', x_title='Bacteria strains', y_title='Virus strains')

###########################################################

Bipartite_MATRIX <- load_bipartite_file_1('mu3e-7_lE-rdw/Bipartite_MATRIX_mu1e-7_initialDiffDp1_S10P15_R-12499_Time_2000.txt','\t')
Bipartite_MATRIX <- Bipartite_MATRIX[,-which(colnames(Bipartite_MATRIX)=='X.1')]
mu3e7_Dp1_S10P15_2000_I <- t(Bipartite_MATRIX)
plot_matrix(mu3e7_Dp1_S10P15_2000_I, layout = 'nested', method = 'ggplot', title='Immunity network from Simlated Data  t = 2000 hr', x_title='Bacteria strains', y_title='Virus strains')

###########################################################


              