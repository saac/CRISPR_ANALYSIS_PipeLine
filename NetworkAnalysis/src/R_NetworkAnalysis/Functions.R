# Plotting ----------------------------------------------------------------

gg_color_hue <- function(n, hue_min = 10, hue_max = 280, l = 62, c = 100) {
  hues = seq(hue_min, hue_max, length=n+1)
  hcl(h=hues, l=l, c=c)[1:n]
}


plot_matrix <- function(M, layout='random', method = 'ggplot', binary_cols=c('gray','red'), title='', x_title='', y_title=''){
  if (layout=='random'){
    M <- M[sample(1:nrow(M), size = nrow(M), replace = F), sample(1:ncol(M), size = ncol(M), replace = F)]
  }
  
  if (layout == "diagonal") {
    ca <- cca(M)
    M <- M[order(summary(ca)$sites[, 1], decreasing = TRUE), 
           order(summary(ca)$species[, 1], decreasing = TRUE)]
  }
  
  if (layout=='nested' & method == 'ggplot'){
    M<- M[order(rowSums(M), decreasing = T), order(colSums(M),decreasing = F)]
  }
  if (layout=='nested' & method == 'heatmap'){
    M<- M[order(rowSums(M), decreasing = T), order(colSums(M),decreasing = T)]
  }
  
  if (max(M)>1){
    colors=c('gray',gg_color_hue(n=max(M)))
  } else {
    colors=binary_cols
  }
  
  if (method=='heatmap'){
    heatmap(M, Rowv = NA, Colv = NA, 
            symm = F, scale = 'none', revC = T, margins=c(7,7), 
            # labRow = F, labCol = F, 
            col=colors, main=title)
    return(invisible())
  }
  
  if (method=='ggplot'){
    M <- reshape2::melt(M)
    p=as_tibble(M) %>% 
      ggplot()+
      geom_tile(aes(Var1,Var2,fill=value))+
      scale_fill_gradientn(colours=colors)+
      theme(
        # axis.text=element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_text(size = 18))+
      labs(title=title, x=x_title,y=y_title)+
      coord_fixed()
    return(p)
  }
}

# Plot with ggplot and color the modules. Must get the matrix Z and the node
# data with columns: nodeName, module, color.

ggplot_bipartite_modules <- function(Z, node_data_obs, module_numbers=F, color_tick_labels=T, border=F, border_width=0.5, text_size=18, title=NULL, xlab=NULL, ylab=NULL, weighted=F){
  
  M_set1 <- M_set2 <- as.tibble(reshape2::melt(Z))
  names(M_set1)[1] <- 'nodeName'
  suppressMessages(suppressWarnings(M_set1 %<>% left_join(node_data_obs[,c('nodeName','color','module')])))
  names(M_set1) <- c('Set1','Set2','value','color1','module1')
  
  names(M_set2)[2] <- 'nodeName'
  suppressMessages(suppressWarnings(M_set2 %<>% left_join(node_data_obs[,c('nodeName','color','module')])))
  names(M_set2) <- c('Set1','Set2','value','color2','module2')
  
  M <- suppressMessages(suppressWarnings(full_join(M_set1, M_set2)))
  Set1_modules <- unique(M[,c('Set1','module1','color1')])
  Set1_modules <- with(Set1_modules, Set1_modules[order(module1,Set1),])
  Set2_modules <- unique(M[,c('Set2','module2','color2')])
  Set2_modules <- with(Set2_modules, Set2_modules[order(module2,Set2),])
#   module_colors <- c('gray50',unique(Set1_modules$color1))
  
  if (weighted){  
    M %<>% mutate(edge_in_out=ifelse(module1==module2,'in','out')) %>% 
        mutate(value_mod=ifelse(edge_in_out=='in',module1,0)) %>% 
        mutate(Set1=factor(Set1, levels=Set1_modules$Set1), Set2=factor(Set2, levels=Set2_modules$Set2)) 
  } else {    
    M %<>% filter(value==1) %>% mutate(edge_in_out=ifelse(module1==module2,'in','out')) %>% 
        mutate(value_mod=ifelse(edge_in_out=='in',module1,0)) %>% 
        mutate(Set1=factor(Set1, levels=Set1_modules$Set1), Set2=factor(Set2, levels=Set2_modules$Set2))     
  }     

  if (all(M$edge_in_out=='in')){
    module_colors <- unique(Set1_modules$color1)
  } else {
    module_colors <- c('gray50',unique(Set1_modules$color1))
  }      
    
    
  if (border){
    p <- M %>% 
      ggplot()+
      geom_tile(aes(Set1,Set2,fill=factor(value_mod)), color='black', size=border_width)
  } else {
    p <- M %>% 
      ggplot()+
      geom_tile(aes(Set1,Set2,fill=factor(value_mod)))
  }
  p <- p +
    scale_fill_manual(values=module_colors)+
    theme(
      panel.background = element_blank(),
      panel.border = element_rect(fill = 'transparent',colour = "black",size = 1),
      axis.text.x=element_text(angle=-90),
      legend.position='none',
      axis.ticks = element_blank(),
      axis.title = element_text(size = text_size))+
    coord_fixed()
  
  if (module_numbers){p <- p+geom_text(aes(Set1,Set2,label = value_mod))}
  if (color_tick_labels){p <- p+theme(axis.text.x=element_text(color=Set1_modules$color1),
                                      axis.text.y=element_text(color=Set2_modules$color2))}
  if (!is.null(title)){p <- p+labs(title=title)}
  if (!is.null(xlab) | !is.null(xlab)){p <- p+labs(x=xlab,y=ylab)}
  
  return(p)
}

# Files -------------------------------------------------------------------

load_bipartite_file_1 <- function(file,s,method='csv'){
  if (method=='csv'){
    M <- read.csv(file,sep = s)
    rownames(M) <- M$X
    M <- data.matrix(M[,-1])
    print(M[1:4,1:4])
    # print(range(M))
  }
  
  if(method=='tab'){
    M <- read.table(file, sep = '\t', header = T, row.names = 1, comment.char = '')
    M <- data.matrix(M)
  }
  return(M)
}

sink.reset <- function(){
  for(i in seq_len(sink.number())){
    sink(NULL)
  }
}

# Data structures ---------------------------------------------------------

incidence_matrix_to_list <- function(M){
  require(bipartite)
  web2edges(M, webName = 'tmp')
  nodes <- read_delim('tmp-names.lut', delim = '\t', col_names = c('nodeID','nodeName'), skip = 1)
  edge_list <- read_delim('tmp.pairs', delim = '\t', col_names = c('i','j','w'))
  web2edges(M, webName = 'tmp', out.files = 'groups')
  groups <- read_delim('tmp.groups', delim = '\t', col_names = c('nodeName','group'))
  nodes <- inner_join(nodes,groups)
  unlink('*.pairs')
  unlink('*.lut')
  unlink('*.groups')
  return(list(M=M, edge_list=edge_list,nodes=nodes))
}


# # A function to calculate the probabilistic model from Bascompte 2003, Olesen
# # 2007, Fortuna 2010, etc. Output is a single shuffled matrix.
# prob_model <- function(M){
#   c <- ncol(M)
#   r <- nrow(M)
#   prob_row <- rowSums(M)/c # calculate intrarow probabilities 
#   prob_column <- colSums(M)/r # calculate intracolumn probabilities 
#   P <- matrix(0,r,c) # matrix of probabilities of interactions
#   for (i in 1:r){
#     for(j in 1:c){
#       P[i,j] <- mean(c(prob_row[i],prob_column[j]))
#     }
#   }
#   
#   R <- matrix(runif(r*c),r,c) # M random matrix to compare to
#   S <- (P >= R)*1 # The shuffled matrix
#   return(S)
# }
# 
# # A wrapper function for the probabilistic model
# prob_model_wrapper <- function(M, non_zero_rc_thershold=1, nsim=100, write_files=F, file_prefix, folder){
#   shuffled_prob <- list()
#   i <- 1
#   while (i <= nsim){
#     print(paste('Creating network with probabilistic model ',i,'/',nsim,sep=''))
#     x <- prob_model(M)
#     # Because the model is probabilistic it can create degenerate matrices with
#     # rows and/or columns what have no interactinos (sum is 0). Only include the
#     # shuffled matrix if the proportion of non-zero rows and columns it has is
#     # beyond the threshold. By default all rows and cols must have at least one
#     # interaction (threshold of 1).
#     prop_non_zero_rows <- sum(rowSums(x)!=0)/nrow(x)
#     prop_non_zero_cols <- sum(colSums(x)!=0)/ncol(x)
#     if(prop_non_zero_rows>=non_zero_rc_thershold && prop_non_zero_cols>=non_zero_rc_thershold){
#       shuffled_prob[[i]] <- x
#       i <- i+1
#     }
#   }
#   
#   if (write_files){
#     for(i in 1:nsim){
#       x <- shuffled_prob[[i]]
#       file <- paste(folder,'/',file_prefix,'_prob_',i,'.csv',sep='')
#       print(file)
#       write.csv(x, file)
#     }
#   }
#   return(shuffled_prob)
# }
# 
# # A function to shuffle a matrix with multiple methods
# shuffle_bipartite_matrix <- function(M, shuff_methods, nsim=10, burnin=1000, write_files=F, file_prefix, folder){
#   shuffled_list <- list()
#   for (sm in shuff_methods){
#     print(sm)
#     null_model <- vegan::nullmodel(M, method = sm)
#     shuffled <- simulate(null_model, nsim = nsim, burnin = burnin)
#     shuffled_list[[which(shuff_methods==sm)]] <- shuffled
#   }
#   names(shuffled_list) <- shuff_methods
#   
#   if (write_files){
#     for (sm in shuff_methods){
#       for(i in 1:nsim){
#         x <- shuffled_list[[which(shuff_methods==sm)]][,,i]
#         file <- paste(folder,'/',file_prefix,'_',sm,'_',i,'.csv',sep='')
#         print(file)
#         write.csv(x, file)
#       }
#     }  
#   }
#   return(shuffled_list)
# }
# 
# # Read shuffled networks to a list
# read_shuffled_networks <- function(shuff_methods, nsim=10, file_prefix, folder){
#   require(data.table)
#   shuffled_list <- vector(mode = 'list', length(shuff_methods))
#   names(shuffled_list) <- shuff_methods
#   for (sm in shuff_methods){
#     for (i in 1:nsim){
#       file <- paste(folder,'/',file_prefix,'_',sm,'_',i,'.csv',sep='')
#       x <- fread(file, sep = ',', header = F, stringsAsFactors = F, skip = 1, drop = 1)
#       # Create an array when reading the first matrix
#       if (is.null(shuffled_list[[which(sm==shuff_methods)]])){
#         shuffled_list[[which(sm==shuff_methods)]] <- array(0, dim=list(nrow(x),ncol(x),nsim))
#       }
#       shuffled_list[[which(sm==shuff_methods)]][,,i] <- data.matrix(x)
#     }
#   }
#   return(shuffled_list)
# }
# 
# 
# # Nestedness --------------------------------------------------------------
# 
# nestedness_analysis <- function(M, nsim = 1000, null_model_method, plotit=F){
#   # Nestedness for binary matrices only
#   M[M>0] <- 1
#   # Get observed value
#   nodf <- nestednodf(M)
#   nodf_obs <- nodf$statistic[3]
#   # Plot
#   if (plotit){
#     print(plot_matrix(M = M, layout = 'nested', method = 'ggplot'))
#   }
#   
#   # If a list of shuffled networks is supplied
#   if (class(null_model_method)=='list'){ 
#     nodf_shuffled <- lapply(null_model_method, nestednodf)
#     nodf_shuffled <- sapply(nodf_shuffled, function(x) x$statistic[3])
#     pvalue <- sum(nodf_shuffled>nodf$statistic[3])/length(nodf_shuffled)
#   }
#   
#   # If an array of shuffled networks is supplied (rows x cols x num simulations)
#   if (class(null_model_method)=='array'){ 
#     nodf_shuffled <- apply(null_model_method, MARGIN = 3, nestednodf)
#     nodf_shuffled <- sapply(nodf_shuffled, function(x) x$statistic[3])
#     pvalue <- sum(nodf_shuffled>nodf$statistic[3])/length(nodf_shuffled)
#   }
#   
#   # If using built-in models of shuffling
#   if (class(null_model_method)=='character'){ #
#     nodf_eval <- oecosimu(M, nestednodf, method = null_model_method, nsimul = nsim, burnin = 1000)
#     pvalue <- nodf_eval$oecosimu$pval[3]
#     nodf_shuffled=nodf_eval$oecosimu$simulated[3,]
#   }
#   
#   # If no shuffling
#   if (class(null_model_method)=='NULL'){ #
#     nodf_eval <- NULL
#     pvalue <- NULL
#     nodf_shuffled <- NULL
#   }
#   
#   return(list(pvalue=pvalue, nodf_obs=nodf_obs, nodf_shuffled=nodf_shuffled))
# }
# 
# # Wrapper functions to automate nestedness calculations
# nestedness_wrapper_prob <- function(M, nsim=100, non_zero_rc_thershold = 0.8){
#   print('Creating shuffled networks...')
#   M_shuffled_prob <- prob_model_wrapper(M, non_zero_rc_thershold = non_zero_rc_thershold, nsim = nsim)
#   print('Calculating nestedness for observed and shuffled networks...')
#   M_nestedness <- nestedness_analysis(M, null_model_method=M_shuffled_prob)
#   pvalue_hist <- as.tibble(M_nestedness$nodf_shuffled) %>% ggplot(aes(value))+geom_histogram()+geom_vline(xintercept = M_nestedness$nodf_obs)
#   return(list(pvalue_hist=pvalue_hist, pvalue=M_nestedness$pvalue))
# }


# # Modularity (Q) ----------------------------------------------------------
# 
# calculate_Q <- function(M){
#   mod <- computeModules(M)
#   slotNames(mod) # see ?moduleWeb for details
#   return(mod@likelihood) # This is the value of the modularity function Q. NOTICE THE @ SIGN (instead of $).
# }
# 
# modularity_analysis <- function(M, num_iterations = 100, null_model_method='curveball'){
#   null_model <- vegan::nullmodel(M, method = null_model_method)
#   shuffled <- simulate(null_model, nsim = num_iterations, burnin=10000)
#   Q_randomized <- apply(shuffled, 3, calculate_Q) # This produces an array with 10 matrices, each of them is a shuffled matrix.
#   Q_obs <- calculate_Q(M)
#   plt_hist <- as.tibble(Q_randomized) %>% 
#     ggplot()+geom_histogram(aes(value))+geom_vline(xintercept = Q_obs, color='red')
#   p_value <- sum(Q_randomized>Q_obs)/num_iterations
#   mod <- computeModules(M)
#   return(list(p_value,plt_hist, mod))
# }


# Infomap -----------------------------------------------------------------

get_node_data_infomap <- function(M, is_bipartite=T, bipartite_groups, metadata=NULL){
  # The function assumes row and column names are present!
  if (is.null(rownames(M)) | is.null(colnames(M))){
    stop('Matrix must have row and column names! Aborting.')
  }
  # Still need to code the unipartite version
  
  if(is_bipartite){
    # Start with the first set (in rows)
    node_data <- tibble(nodeName=rownames(M))
    node_data %<>% arrange(nodeName)
    node_data$group <- bipartite_groups[1]
    node_data$InfomapName <- paste('"',node_data$nodeName,'"',sep = '')
    # Now second set (in columns)
    set2 <- tibble(nodeName=colnames(M))
    set2 %<>% arrange(nodeName)
    set2$group <- bipartite_groups[2]
    set2$InfomapName <- paste('"',set2$nodeName,'"',sep = '')
    # Join sets
    node_data %<>% bind_rows(set2)
    node_data$runningID <- 1:nrow(node_data)
    # Add additional metadata if exists
    if(!is.null(metadata)){
      node_data %<>% bind_cols(metadata)
    }
  }
  return(node_data)
}

write_infomap <- function(M, node_data, is_bipartite=F, directed=F, filename){
  # !!! STILL NEED TO TEST THIS FOR UNIPARTITE NETWORKS AND FOR DIRECTED NETWORKS !!!

  #This function prepares data to use with Infomap's unipartite link list. (see
  #http://www.mapequation.org/code.html#Link-list-format). If this is used for a
  #bipartite network, then the incidence matrix is tranformed into an adjacency
  #matrix via transposing. This is done to find modules which contain nodes from
  #both sets. To find modules in one set based on connections to the other use
  #function write_infomap_bipartite. The book-keeping on the bipartite names is
  #done by using ids 1:m for the firs set and (m+1):k for the second, where m is
  #the number of nodes in the first set and k is the total number of nodes in
  #the network. This hshould be considered when parsing results. node_data must
  #have columns: nodeName (character), and runningID (sequential number from 1
  #to the total number of nodes).
  
  g <- graph.incidence(M, directed = directed, mode = 'all', weighted = T)
  links <- igraph::as_data_frame(g, what = 'edges')
  if (is_bipartite){
    M_transposed <- t(M)
    g <- graph.incidence(M_transposed, directed = directed, mode = 'all', weighted = T)
    links_transposed <- igraph::as_data_frame(g, what = 'edges')
    links <- rbind(links, links_transposed)
  }
  
  links$from <- node_data$runningID[match(links$from, node_data$nodeName)]
  links$to <- node_data$runningID[match(links$to, node_data$nodeName)]
  
  unlink(filename)
  sink(filename, append = T)
  if(is_bipartite){
    cat("# A bipartite network parsed as AA^T matrix");cat('\n')
  } else {
    cat("# A unipartite network");cat('\n')
  }
  cat(paste("*Vertices",max(node_data$runningID)));cat('\n')
  write.table(subset(node_data, select=c(runningID, InfomapName)), filename, sep = ' ', row.names = F, col.names = F, quote = F, append = T)
  # write.table(features[,c(1,3)], 'infomap_input.txt', sep = ' ', row.names = F, col.names = F, quote = F, append = T)
  cat("*Edges ");cat(nrow(links));cat('\n')
  write.table(links, filename, sep = ' ', row.names = F, col.names = F, quote = F, append = T)
  sink.reset()
}

# write_infomap_simple <- function(M, is_bipartite=T, directed=F, filename){
#   # This function is a simple version of write_infomap because it does not
#   # receive any information on the nodes. It is used for the shuffled networks.
#   N <- nrow(M)+ncol(M)
#   rownames(M) <- 1:nrow(M)
#   colnames(M) <- (nrow(M)+1):N
#   g <- graph.incidence(M, directed = directed, mode = 'all', weighted = T)
#   links <- igraph::as_data_frame(g, what = 'edges')
#   if (is_bipartite){
#     M_transposed <- t(M)
#     g <- graph.incidence(M_transposed, directed = directed, mode = 'all', weighted = T)
#     links_transposed <- igraph::as_data_frame(g, what = 'edges')
#     links <- rbind(links, links_transposed)
#   }
#   
#   unlink(filename)
#   sink(filename, append = T)
#   if(is_bipartite){
#     cat("# A bipartite network parsed as AA^T matrix | No node names | shuffled");cat('\n')
#   } else {
#     cat("# A unipartite network");cat('\n')
#   }
#   cat(paste("*Vertices",N));cat('\n')
#   write.table(data.frame(node=1:N,id=1:N), filename, sep = ' ', row.names = F, col.names = F, quote = F, append = T)
#   # write.table(features[,c(1,3)], 'infomap_input.txt', sep = ' ', row.names = F, col.names = F, quote = F, append = T)
#   cat("*Edges ");cat(nrow(links));cat('\n')
#   write.table(links, filename, sep = ' ', row.names = F, col.names = F, quote = F, append = T)
#   sink.reset()
# }

# # This function is a wrapper for write_infomap and write_infomap_simple. It
# # receives an observed matrix and an array of corresponding shuffled matrices
# # (row x col x randomization). It calculates the map equation and returns
# # p-values and histograms.
# # Infomap_wrapper <- function(Z, shuffled_matrices, bipartite_groups, file_prefix, infomap_executable='Infomap_01923'){
# Infomap_wrapper <- function(Z, shuffled_matrices, bipartite_groups, file_prefix, infomap_executable='Infomap'){
#   node_data <- get_node_data_infomap(Z, is_bipartite = T, bipartite_groups = bipartite_groups)
#   file_obs <- paste(file_prefix,'_Infomap.txt',sep='')
#   write_infomap(Z, node_data, is_bipartite = T, directed = F, file_obs)
#   system(paste('./',infomap_executable,' ',file_obs,' . -i pajek -N 20 --tree -2 --silent',sep=''))
#   file_obs_output <- str_replace(file_obs, 'txt','tree')
#   node_data_obs <- parse_modules(file_obs_output, infomap_bipartite_format = F, node_data, T)
#   
#   for (sm in names(shuffled_matrices)){
#     for(i in 1:dim(shuffled_matrices[[1]])[3]){
#       x <- shuffled_matrices[[which(names(shuffled_matrices)==sm)]][,,i]
#       file_shuff <- paste('shuffled/',file_prefix,'_',sm,'_',i,'_infomap.txt',sep='')
#       print(file_shuff)
#       write_infomap_simple(x, is_bipartite = T, directed = F, filename = file_shuff)
#       system(paste('./',infomap_executable,' ',file_shuff,' shuffled/ -i pajek -N 1 --tree -2 --silent',sep=''))
#     }
#   }  
#   
#   # Compare codelength between observed and shuffled networks
#   x <- readLines(file_obs_output)[1]
#   codelength_obs <- parse_number(str_sub(x, str_locate_all(x,'codelength ')[[1]][2,1],str_length(x)))
#   codelength_obs_1level <- parse_number(str_sub(x, str_locate_all(x,'codelength ')[[1]][1,1],str_length(x)))
#   codelength_obs_normalized <- codelength_obs/codelength_obs_1level
#   # tree_files <- list.files(pattern='.tree', path = '~/Documents/CRISPR/shuffled/', full.names = T)
#   tree_files <- list.files(pattern='.tree', path = 'shuffled/', full.names = T)
#   tree_files <- tree_files[str_detect(tree_files, file_prefix)]
#   codelength_table <- map(tree_files,
#                           function(t){
#                             x <- str_split(t,'_')[[1]]
#                             run <- parse_number(x[4])
#                             sm <- x[3]                         
#                             x <- readLines(t)[1]
#                             codelength <- parse_number(str_sub(x, str_locate_all(x,'codelength ')[[1]][2,1],str_length(x)))
#                             codelength_1level <- parse_number(str_sub(x, str_locate_all(x,'codelength ')[[1]][1,1],str_length(x)))
#                             codelength_normalized <- codelength/codelength_1level
#                             return(tibble(shuff_method=sm, run , codelength, codelength_1level, codelength_normalized))
#                           }) %>% bind_rows()
#   
#   # It is crucial that the condition here would be <= because if there is one
#   # module than the normalized value is always 1.
#   p_value_table <- codelength_table %>% 
#     mutate(test=codelength_normalized <= codelength_obs_normalized) %>% 
#     group_by(shuff_method) %>% summarise(cl_obs=codelength_obs_normalized, 
#                                          cln_shuff_mean=mean(codelength_normalized), 
#                                          cln_shuff_sd=sd(codelength_normalized), 
#                                          pvalue=sum(test)/length(test))
#   
#   p_value_plot <- codelength_table %>% ggplot(aes(codelength_normalized,fill=shuff_method))+
#     geom_histogram(binwidth=0.02)+
#     facet_wrap(~shuff_method)+
#     geom_vline(xintercept = codelength_obs_normalized)
#   
#   return(list(node_data_obs=node_data_obs, p_value_table=p_value_table, p_value_plot=p_value_plot))
# }


# 
# write_infomap_bipartite <- function(M, node_data, filename){
#   # This function prepares data to use with Infomap's bipartite structure of
#   # features and nodes (see
#   # http://www.mapequation.org/code.html#Bipartite-format). It assumes nodes are
#   # in the columns and features are in the rows. node_data should be organized
#   # with nodes first and features afterwards. node_data must have columns:
#   # nodeName (character), InfomapCode (nX or fX where n and f are nodes or
#   # features and X is a number (see mapequaion.org)), InfomapGroup (n or f), and
#   # runningID (sequential number from 1 to the total number of nodes).
#   colnames(M) <- node_data$InfomapCode[match(colnames(M), node_data$nodeName)]
#   rownames(M) <- node_data$InfomapCode[match(rownames(M), node_data$nodeName)]
#   g <- graph.incidence(M, directed = F, mode = 'all', weighted = T)
#   links <- igraph::as_data_frame(g, what = 'edges')
#   unlink(filename)
#   sink(filename, append = T)
#   cat("# A bipartite network with node names");cat('\n')
#   cat(paste("*Vertices",ncol(M)));cat('\n')
#   write.table(subset(node_data, InfomapGroup=='n',select=c(runningID, InfomapName)), filename, sep = ' ', row.names = F, col.names = F, quote = F, append = T)
#   # write.table(features[,c(1,3)], 'infomap_input.txt', sep = ' ', row.names = F, col.names = F, quote = F, append = T)
#   cat("*Edges ");cat(nrow(links));cat('\n')
#   write.table(links, filename, sep = ' ', row.names = F, col.names = F, quote = F, append = T)
#   sink.reset()
# }

# write_infomap_bipartite_simple <- function(data, filename){
#   # This function is a simple version of write_infomap_bipartite because it does not
#   # receive any information on the nodes. It is used for the shuffled networks.
#   nodes <- data.frame(id=1:ncol(data),name=colnames(data),name_output=paste('"',colnames(data),'"',sep = ''))
#   features <- data.frame(id=(1:nrow(data)),name=rownames(data), name_output=paste('"',rownames(data),'"',sep = ''))
#   g <- graph.incidence(data, directed = F, mode = 'all', weighted = T)
#   df <- igraph::as_data_frame(g, what = 'edges')
#   df$to <- nodes$id[match(as.character(df$to), as.character(nodes$name))]
#   df$from <- features$id[match(as.character(df$from), as.character(features$name))]
#   df$from <- paste('f',df$from,sep='')
#   df$to <- paste('n',df$to,sep='')
#   unlink(filename)
#   sink(filename, append = T)
#   cat("# A bipartite network with node names");cat('\n')
#   cat(paste("*Vertices",nrow(nodes)));cat('\n')
#   write.table(nodes[,c(1,3)], filename, sep = ' ', row.names = F, col.names = F, quote = F, append = T)
#   # write.table(features[,c(1,3)], 'infomap_input.txt', sep = ' ', row.names = F, col.names = F, quote = F, append = T)
#   cat("*Edges ");cat(nrow(df));cat('\n')
#   write.table(df, filename, sep = ' ', row.names = F, col.names = F, quote = F, append = T)
#   sink.reset()
# }

# Extract module composition and map the modules to nodes
parse_modules <- function(file, infomap_bipartite_format=F, node_data, assign_colors=F){
  # !!! STILL NEED TO PROGRAM THAT FOR UNIPARTITE/DIRECTED NETWORKS !!!
  
  bipartite_groups <- unique(node_data$group)
  
  modules_obs <- read_delim(file, delim = ' ', skip = 2, col_names = c('path', 'flow', 'name', 'node'))
  
  # If parsing modules for runs that were done using Infomap's bipartite format
  # with features and nodes.
  if (infomap_bipartite_format){
    # Make sure this get the features and nodes correctly with the group names and not mixed up.
    modules_obs %<>% select(name, path, InfomapCode=node) %>% rowwise %>% mutate(module=str_split(path,':')[[1]][1]) %>% 
    mutate(group=ifelse(str_detect(InfomapCode,'f'),bipartite_groups[1],bipartite_groups[2])) %>% select(InfomapCode, module)
  } else {
    modules_obs %<>% select(name, path, runningID=node) %>% 
      rowwise %>% 
      mutate(module=str_split(path,':')[[1]][1]) %>% 
      select(runningID, module)
  }
  
  node_data %<>% left_join(modules_obs) %>% mutate(module=as.numeric(module))
  if (assign_colors){
    module_colors <- tibble(module=sort(unique(node_data$module)), color=gg_color_hue(max(node_data$module)))
    node_data %<>% left_join(module_colors)
  }
  return(node_data)
}


# #########################################################################

Infomap_wrapper_NoShuffled <- function(Z, bipartite_groups, file_prefix, infomap_executable='Infomap'){
  node_data <- get_node_data_infomap(Z, is_bipartite = T, bipartite_groups = bipartite_groups)
  file_obs <- paste(file_prefix,'_Infomap.txt',sep='')
  write_infomap(Z, node_data, is_bipartite = T, directed = F, file_obs)
#   system(paste('./',infomap_executable,' ',file_obs,' . -i pajek -N 20 --tree -2 --silent',sep=''))
  system(paste('./',infomap_executable,' ',file_obs,' . -i pajek -N 100 --tree -2 --silent',sep=''))
#   system(paste('./',infomap_executable,' ',file_obs,' . -i pajek -s $(echo $RANDOM) -N 100 --tree -2 --silent',sep=''))
  file_obs_output <- str_replace(file_obs, 'txt','tree')
  node_data_obs <- parse_modules(file_obs_output, infomap_bipartite_format = F, node_data, T)
  
  # for (sm in names(shuffled_matrices)){
  #   for(i in 1:dim(shuffled_matrices[[1]])[3]){
  #     x <- shuffled_matrices[[which(names(shuffled_matrices)==sm)]][,,i]
  #     file_shuff <- paste('shuffled/',file_prefix,'_',sm,'_',i,'_infomap.txt',sep='')
  #     print(file_shuff)
  #     write_infomap_simple(x, is_bipartite = T, directed = F, filename = file_shuff)
  #     system(paste('./',infomap_executable,' ',file_shuff,' shuffled/ -i pajek -N 1 --tree -2 --silent',sep=''))
  #   }
  # }  
  # 
  # # Compare codelength between observed and shuffled networks
  # x <- readLines(file_obs_output)[1]
  # codelength_obs <- parse_number(str_sub(x, str_locate_all(x,'codelength ')[[1]][2,1],str_length(x)))
  # codelength_obs_1level <- parse_number(str_sub(x, str_locate_all(x,'codelength ')[[1]][1,1],str_length(x)))
  # codelength_obs_normalized <- codelength_obs/codelength_obs_1level
  # # tree_files <- list.files(pattern='.tree', path = '~/Documents/CRISPR/shuffled/', full.names = T)
  # tree_files <- list.files(pattern='.tree', path = 'shuffled/', full.names = T)
  # tree_files <- tree_files[str_detect(tree_files, file_prefix)]
  # codelength_table <- map(tree_files,
  #                         function(t){
  #                           x <- str_split(t,'_')[[1]]
  #                           run <- parse_number(x[4])
  #                           sm <- x[3]                         
  #                           x <- readLines(t)[1]
  #                           codelength <- parse_number(str_sub(x, str_locate_all(x,'codelength ')[[1]][2,1],str_length(x)))
  #                           codelength_1level <- parse_number(str_sub(x, str_locate_all(x,'codelength ')[[1]][1,1],str_length(x)))
  #                           codelength_normalized <- codelength/codelength_1level
  #                           return(tibble(shuff_method=sm, run , codelength, codelength_1level, codelength_normalized))
  #                         }) %>% bind_rows()
  # 
  # # It is crucial that the condition here would be <= because if there is one
  # # module than the normalized value is always 1.
  # p_value_table <- codelength_table %>% 
  #   mutate(test=codelength_normalized <= codelength_obs_normalized) %>% 
  #   group_by(shuff_method) %>% summarise(cl_obs=codelength_obs_normalized, 
  #                                        cln_shuff_mean=mean(codelength_normalized), 
  #                                        cln_shuff_sd=sd(codelength_normalized), 
  #                                        pvalue=sum(test)/length(test))
  # 
  # p_value_plot <- codelength_table %>% ggplot(aes(codelength_normalized,fill=shuff_method))+
  #   geom_histogram(binwidth=0.02)+
  #   facet_wrap(~shuff_method)+
  #   geom_vline(xintercept = codelength_obs_normalized)
  # 
  # return(list(node_data_obs=node_data_obs, p_value_table=p_value_table, p_value_plot=p_value_plot))
  return(list(node_data_obs=node_data_obs))
}



