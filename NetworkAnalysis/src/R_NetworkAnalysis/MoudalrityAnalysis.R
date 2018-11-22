

# Create shuffled networks ------------------------------------------------

## @knitr shuffle_networks
nsim=100
folder_shuffled <- 'shuffled/'
shuffle_models <- c('r0','c0','curveball')

# file_prefix <- 'Russia2000_A'
# x <- shuffle_bipartite_matrix(Russia2000_A, shuffle_models, nsim=nsim, burnin=1000, write_files=T, file_prefix=file_prefix, folder=folder_shuffled)
# y <- prob_model_wrapper(Russia2000_A, non_zero_rc_thershold = 0.85, nsim = nsim, write_files = T, file_prefix = file_prefix, folder = folder_shuffled)

file_prefix <- 'mu3e7_Dp1_S10P15_1500'
x <- shuffle_bipartite_matrix(mu3e7_Dp1_S10P15_1500, shuffle_models, nsim=nsim, burnin=1000, write_files=T, file_prefix=file_prefix, folder=folder_shuffled)
y <- prob_model_wrapper(mu3e7_Dp1_S10P15_1500, non_zero_rc_thershold = 0.85, nsim = nsim, write_files = T, file_prefix = file_prefix, folder = folder_shuffled)


file_prefix <- 'mu3e7_Dp1_S10P15_1750'
x <- shuffle_bipartite_matrix(mu3e7_Dp1_S10P15_1750, shuffle_models, nsim=nsim, burnin=1000, write_files=T, file_prefix=file_prefix, folder=folder_shuffled)
y <- prob_model_wrapper(mu3e7_Dp1_S10P15_1750, non_zero_rc_thershold = 0.85, nsim = nsim, write_files = T, file_prefix = file_prefix, folder = folder_shuffled)

file_prefix <- 'mu3e7_Dp1_S10P15_2000'
x <- shuffle_bipartite_matrix(mu3e7_Dp1_S10P15_2000, shuffle_models, nsim=nsim, burnin=1000, write_files=T, file_prefix=file_prefix, folder=folder_shuffled)
y <- prob_model_wrapper(mu3e7_Dp1_S10P15_2000, non_zero_rc_thershold = 0.85, nsim = nsim, write_files = T, file_prefix = file_prefix, folder = folder_shuffled)




# file_prefix <- 'Russia2000_V'
# x <- shuffle_bipartite_matrix(Russia2000_V, shuffle_models, nsim=nsim, burnin=1000, write_files=T, file_prefix=file_prefix, folder=folder_shuffled)
# y <- prob_model_wrapper(Russia2000_V, non_zero_rc_thershold = 0.95, nsim = nsim, write_files = T, file_prefix = file_prefix, folder = folder_shuffled)
# 
# file_prefix <- 'Russia2000_I'
# x <- shuffle_bipartite_matrix(Russia2000_I, shuffle_models, nsim=nsim, burnin=1000, write_files=T, file_prefix=file_prefix, folder=folder_shuffled)
# y <- prob_model_wrapper(Russia2000_I, non_zero_rc_thershold = 1, nsim = nsim, write_files = T, file_prefix = file_prefix, folder = folder_shuffled)
# 
# file_prefix <- 'Russia2000_A2'
# x <- shuffle_bipartite_matrix(Russia2000_A2, shuffle_models, nsim=nsim, burnin=1000, write_files=T, file_prefix=file_prefix, folder=folder_shuffled)
# y <- prob_model_wrapper(Russia2000_A2, non_zero_rc_thershold = 0.8, nsim = nsim, write_files = T, file_prefix = file_prefix, folder = folder_shuffled)
# 
# file_prefix <- 'Yellowstone_A'
# x <- shuffle_bipartite_matrix(Yellowstone_A, shuffle_models, nsim=nsim, burnin=1000, write_files=T, file_prefix=file_prefix, folder=folder_shuffled)
# y <- prob_model_wrapper(Yellowstone_A, non_zero_rc_thershold = 0.95, nsim = nsim, write_files = T, file_prefix = file_prefix, folder = folder_shuffled)

# Now also include a model for quantitative networks for the immunity network
# file_prefix <- 'Russia2000_I'
# x <- shuffle_bipartite_matrix(Russia2000_I, shuffle_models, nsim=nsim, burnin=1000, write_files=T, file_prefix=file_prefix, folder=folder_shuffled)


# Read shuffled matrices --------------------------------------------------
## @knitr get_shuffled_matrices

# nsim=100
# folder_shuffled <- 'shuffled/'
# shuffle_models <- c('r0','c0','prob','curveball')
# 
# file_prefix <- 'Russia2000_A'
# Russia2000_A_shuffled <- read_shuffled_networks(shuff_methods = shuffle_models, file_prefix = file_prefix, nsim = nsim, folder = folder_shuffled)
# sapply(Russia2000_A_shuffled, dim)

file_prefix <- 'mu3e7_Dp1_S10P15_1500'
mu3e7_Dp1_S10P15_1500_shuffled <- read_shuffled_networks(shuff_methods = shuffle_models, file_prefix = file_prefix, nsim = nsim, folder = folder_shuffled)
sapply(mu3e7_Dp1_S10P15_1500_shuffled, dim)

file_prefix <- 'mu3e7_Dp1_S10P15_1750'
mu3e7_Dp1_S10P15_1750_shuffled <- read_shuffled_networks(shuff_methods = shuffle_models, file_prefix = file_prefix, nsim = nsim, folder = folder_shuffled)
sapply(mu3e7_Dp1_S10P15_1750_shuffled, dim)

file_prefix <- 'mu3e7_Dp1_S10P15_2000'
mu3e7_Dp1_S10P15_2000_shuffled <- read_shuffled_networks(shuff_methods = shuffle_models, file_prefix = file_prefix, nsim = nsim, folder = folder_shuffled)
sapply(mu3e7_Dp1_S10P15_2000_shuffled, dim)

# file_prefix <- 'Russia2000_V'
# Russia2000_V_shuffled <- read_shuffled_networks(shuff_methods = shuffle_models, file_prefix = file_prefix, nsim = nsim, folder = folder_shuffled)
# sapply(Russia2000_V_shuffled, dim)
# 
# file_prefix <- 'Russia2000_I'
# Russia2000_I_shuffled <- read_shuffled_networks(shuff_methods = shuffle_models, file_prefix = file_prefix, nsim = nsim, folder = folder_shuffled)
# sapply(Russia2000_I_shuffled, dim)
# 
# file_prefix <- 'Russia2000_A2'
# Russia2000_A2_shuffled <- read_shuffled_networks(shuff_methods = shuffle_models, file_prefix = file_prefix, nsim = nsim, folder = folder_shuffled)
# sapply(Russia2000_A2_shuffled, dim)
# 
# file_prefix <- 'Yellowstone_A'
# Yellowstone_A_shuffled <- read_shuffled_networks(shuff_methods = shuffle_models, file_prefix = file_prefix, nsim = nsim, folder = folder_shuffled)
# sapply(Yellowstone_A_shuffled, dim)


# # Infomap -----------------------------------------------------------------
# 
# ## @knitr Infomap_Russia_2000_A
# print('Russia 2000 acquisition:')
# x <- Infomap_wrapper(Z = Russia2000_A, shuffled_matrices = Russia2000_A_shuffled, bipartite_groups = c('spacer','archaea'), file_prefix = 'Russia2000_A')
# x$p_value_table
# x$p_value_plot+labs(title='Russia 2000 acquisition')
# ggplot_bipartite_modules(Z=Russia2000_A, x$node_data_obs, module_numbers = F, color_tick_labels = T, border=F, text_size = 18, title='Spacer acquisition modularity', xlab='Spacer ID', ylab='Archaea strain')
 
## @knitr Infomap_mu3e7_Dp1_S10P15_1500_A
print('mu3e7_Dp1_S10P15_1500 acquisition:')
x <- Infomap_wrapper(Z = mu3e7_Dp1_S10P15_1500, shuffled_matrices = mu3e7_Dp1_S10P15_1500_shuffled, bipartite_groups = c('spacer','archaea'), file_prefix = 'mu3e7_Dp1_S10P15_1500')
x$p_value_table
x$p_value_plot+labs(title='mu3e7_Dp1_S10P15_1500 acquisition')
ggplot_bipartite_modules(Z=mu3e7_Dp1_S10P15_1500, x$node_data_obs, module_numbers = F, color_tick_labels = T, border=F, text_size = 18, title='Spacer acquisition modularity t = 1500 hr', xlab='Spacer ID', ylab='Bacteria strain')

## @knitr Infomap_mu3e7_Dp1_S10P15_1750_A
print('mu3e7_Dp1_S10P15_1750 acquisition:')
x <- Infomap_wrapper(Z = mu3e7_Dp1_S10P15_1750, shuffled_matrices = mu3e7_Dp1_S10P15_1750_shuffled, bipartite_groups = c('spacer','archaea'), file_prefix = 'mu3e7_Dp1_S10P15_1750')
x$p_value_table
x$p_value_plot+labs(title='mu3e7_Dp1_S10P15_1750 acquisition')
ggplot_bipartite_modules(Z=mu3e7_Dp1_S10P15_1750, x$node_data_obs, module_numbers = F, color_tick_labels = T, border=F, text_size = 18, title='Spacer acquisition modularity t = 1750 hr', xlab='Spacer ID', ylab='Bacteria strain')

## @knitr Infomap_mu3e7_Dp1_S10P15_2000_A
print('mu3e7_Dp1_S10P15_2000 acquisition:')
x <- Infomap_wrapper(Z = mu3e7_Dp1_S10P15_2000, shuffled_matrices = mu3e7_Dp1_S10P15_2000_shuffled, bipartite_groups = c('spacer','archaea'), file_prefix = 'mu3e7_Dp1_S10P15_2000')
x$p_value_table
x$p_value_plot+labs(title='mu3e7_Dp1_S10P15_2000 acquisition')
ggplot_bipartite_modules(Z=mu3e7_Dp1_S10P15_2000, x$node_data_obs, module_numbers = F, color_tick_labels = T, border=F, text_size = 18, title='Spacer acquisition modularity t = 2000 hr', xlab='Spacer ID', ylab='Bacteria strain')







# ## @knitr Infomap_Russia_2000_V
# print('Russia 2000 virus-protospacer:')
# x <- Infomap_wrapper(Z = Russia2000_V, shuffled_matrices = Russia2000_V_shuffled, bipartite_groups = c('spacer','virus'), file_prefix = 'Russia2000_V')
# x$p_value_table
# x$p_value_plot+labs(title='Russia 2000 virus-protospacer')
# ggplot_bipartite_modules(Z=Russia2000_V, x$node_data_obs, module_numbers = F, color_tick_labels = T, border=F, text_size = 18, title='Virus-protospacer modularity', xlab='Spacer ID', ylab='Virus strain')
# 
# ## @knitr Infomap_Russia_2000_I
# print('Russia 2000 immunity:')
# x <- Infomap_wrapper(Z = Russia2000_I, shuffled_matrices = Russia2000_I_shuffled, bipartite_groups = c('archaea','virus'), file_prefix = 'Russia2000_I')
# x$p_value_table
# x$p_value_plot+labs(title='Russia 2000 immunity')
# ggplot_bipartite_modules(Z=Russia2000_I, x$node_data_obs, module_numbers = F, color_tick_labels = T, border=F, text_size = 18, title='Immunity network modularity', xlab='Archaea strain', ylab='Virus strain')
# 
# ## @knitr Infomap_Russia_2000_A2
# print('Russia 2000 acquisition (no immunity criteria):')
# x <- Infomap_wrapper(Z = Russia2000_A2, shuffled_matrices = Russia2000_A2_shuffled, bipartite_groups = c('spacer','archaea'), file_prefix = 'Russia2000_A2')
# x$p_value_table
# x$p_value_plot+labs(title='Russia 2000 acquisition (no immunity criteria)')
# # ggplot_bipartite_modules(Z=Russia2000_A2, x$node_data_obs, module_numbers = F, color_tick_labels = T, text_size = 18, title='Spacer acquisition (no immunity criteria) modularity', xlab='Spacer ID', ylab='Archaea strain')
# 
# ## @knitr Infomap_Yellowstone_A
# print('Yellowstone acquisition (no immunity criteria):')
# x <- Infomap_wrapper(Z = Yellowstone_A, shuffled_matrices = Yellowstone_A_shuffled, bipartite_groups = c('spacer','archaea'), file_prefix = 'Yellowstone_A')
# x$p_value_table
# x$p_value_plot+labs(title='Yellowstone acquisition')
# ggplot_bipartite_modules(Z=Yellowstone_A, x$node_data_obs, module_numbers = F, color_tick_labels = T, border=F, text_size = 18, title='Spacer acquisition modularity', xlab='Spacer ID', ylab='Archaea strain')
