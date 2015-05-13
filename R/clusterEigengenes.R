clusterPCA <- function(cluster_frame, expression_frame, n) {
  
  # n = vector of time points
  # cluster_frame = dataframe with genes and cluster assignments
  # expression frame = dataframe with genes as rows and timepoints as columns
  
  clusters <- levels(cluster_frame$cluster)
  genes <- rownames(expression_frame)
  input_frame <- data.frame(expression_frame, genes)
  colnames(input_frame) <- c(n, "genes")
  
  merged <- merge(x=input_frame,
                    y=cluster_frame,
                    by="genes")
  
  rownames(merged) <- merged$genes
  merged <- merged[,-1]
  
  # create a dataframe for each cluster with the relevant genes asssigned to it
  eigenList = list()
  counter = 1
  for(color in clusters){
    col_frame = data.frame(merged[merged$cluster == color,])
    
    # centred and scaled PCA to extract PC1, loadings and prop^ variance explained
    col_pca = prcomp(t(col_frame[1:length(n)]), scale=T, center=T)
    prop_var = summary(col_pca)$importance[2,][1]
    pc1 = -1*col_pca$x[,1]
    pc1_load = col_pca$rotation[,1]
    
    # put into a list for each cluster
    cluster_list = list("cluster" = color, 
                        "eigenExpress" = pc1, 
                        "varExp" = prop_var, 
                        "eigenLoading" = pc1_load)
    
    eigenList[[counter]] = cluster_list
    counter = counter + 1
  }
  return(eigenList)
}

eigenExpress <- function(eigenList, n) {
  # eigenList = return value from clusterPCA
  # n = numeric vector of time points
  # return a dataframe of average eigengene expression
  
  eigens_list = list()
  for(i in 1:length(eigenList)) {
    eigens_list[[i]] = c(eigenList[[i]]$cluster, as.vector(eigenList[[i]]$eigenExpress))
  }
  eigen_frame = data.frame(do.call(rbind, eigens_list))
  colnames(eigen_frame) = c("cluster", n)
  clusters = eigen_frame$cluster
  eigen_frame = eigen_frame[,-1]
  # fix conversion of values to numeric
  
  for(i in 1:length(n)) {
    eigen_frame[,i] <- as.numeric(as.character(eigen_frame[,i]))
  }  
  eigen_frame = data.frame(clusters, eigen_frame)
  colnames(eigen_frame) = c("cluster", n)
  return(eigen_frame)
}

eigenLoad <- function(eigenList, image.dir, condition) {
  # eigenList = return value from clusterPCA
  # genes = character vector of gene names
  # image.dir = directory save .png barplot images
  # generates a plot of eigengene loadings

  loadings_list = list()
  for(i in 1:length(eigenList)) {
    load_vec = as.numeric(abs(eigenList[[i]]$eigenLoading))[order(as.numeric(eigenList[[i]]$eigenLoading))]
    save_location = paste0(image.dir,
                           "/",
                           condition,
                           "-",
                           eigenList[[i]]$cluster,
                           "-eigengene_loadings")
        
    png(paste0(save_location, ".png"))
    barplot(load_vec, main = paste0(eigenList[[i]]$cluster , " eigengene loadings"))
    dev.off()
    
    # save loadings as .tsv in images.dir
    
    write.table(eigenList[[i]]$eigenLoading, file=paste0(save_location,".tsv"), sep="\t")
  }
}

eigenPlot <- function(eigenFrame, image.dir, condition) {
  # eigenFrame = dataframe of eigengene expression (output from eigenExpress)
  # image.dir = directory to save images to
  # condition = experimental condition
  # plots and saves ggplots of eigengene expression profiles
  
  suppressPackageStartupMessages(library(reshape2))
  suppressPackageStartupMessages(library(ggplot2))
  
  melted = melt(eigenFrame, id.vars="cluster")
  melted$variable = as.numeric(as.character(melted$variable))
  melted$value = as.numeric(melted$value)
  
  eigen_plot = ggplot(melted, aes(x=variable, y=value, colour=cluster)) + 
    geom_line(aes(group=cluster)) + geom_point(aes(group=cluster)) + 
    labs(x="time (hours)", y="expression", title=paste0(condition, " module eigengene expression")) + 
    scale_x_continuous(limits=c(0, 124), breaks=melted$variable)
  
  png(paste0(image.dir,"/",condition,"-","eigengene_expression.png"))
  print(eigen_plot)
  dev.off()
          
}
  
  
