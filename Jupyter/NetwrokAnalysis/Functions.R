
# The following are a set of functions that have been originaly developed by Shai Pilosof (https://github.com/shainova) for the Network Analysis of Bipartite Networks.
# These functions are still under development in Mercedes Pascual Group (https://github.com/pascualgroup).

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
    M %<>% filter(value!=0) %>% mutate(edge_in_out=ifelse(module1==module2,'in','out')) %>% 
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


# ###################################  These functions have been developed by Sergio A. Alcala-Corona (https://github.com/saac) ###################################

Infomap_wrapper_NoShuffled <- function(Z, bipartite_groups, file_prefix, infomap_executable='Infomap'){
  node_data <- get_node_data_infomap(Z, is_bipartite = T, bipartite_groups = bipartite_groups)
  file_obs <- paste(file_prefix,'_Infomap.txt',sep='')
  write_infomap(Z, node_data, is_bipartite = T, directed = F, file_obs)
#   system(paste('./',infomap_executable,' ',file_obs,' . -i pajek -N 20 --tree -2 --silent',sep=''))
  system(paste('./',infomap_executable,' ',file_obs,' . -i pajek -N 100 --tree -2 --silent',sep=''))
#   system(paste('./',infomap_executable,' ',file_obs,' . -i pajek -s $(echo $RANDOM) -N 100 --tree -2 --silent',sep=''))
  file_obs_output <- str_replace(file_obs, 'txt','tree')
  node_data_obs <- parse_modules(file_obs_output, infomap_bipartite_format = F, node_data, T)
  
  return(list(node_data_obs=node_data_obs))
}

BuildNetwork <- function(argument,PSp){
    protospacer_by_virus <- load_bipartite_file_1(argument,'\t')
    protospacer_by_virus <- protospacer_by_virus[,str_detect(colnames(protospacer_by_virus), PSp)]
    protospacer_by_virus <- incidence_matrix_to_list(protospacer_by_virus)
    network <- protospacer_by_virus$M
    network <- t(network)
    return(network)
}

BuildBipartiteNetwork<- function(argument){
    Bipartite_MATRIX <- load_bipartite_file_1(argument,'\t')
    Bipartite_MATRIX <- Bipartite_MATRIX[,-which(colnames(Bipartite_MATRIX)=='X.1')]
    network <-  t(Bipartite_MATRIX)
    return(network)
}

