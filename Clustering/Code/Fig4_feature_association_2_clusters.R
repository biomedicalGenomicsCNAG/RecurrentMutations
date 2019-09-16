library(FactoMineR)
library(ggplot2)
library(reshape2)

#' Correct the p-values for the results of the v tests (features versus clusters) for multiple testing.
#' @param res_pca: PCA object
#' @param num_clusters: number of clusters
#' @return data.frame with all statisticially significant associations (p < 0.05) after multiple testing correction.
correctPvalAssociations4MultTest <- function(res_pca, num_clusters){

    num_features <- nrow(res_pca$var$contrib)
    set.seed(10)
    res_hcpc_th_prob1 <- HCPC(res_pca, nb.clust=num_clusters,max=20,graph = FALSE,consol=TRUE, proba=1)
    
    all_pval <- as.data.frame(matrix(nrow=(num_features * num_clusters), ncol=2))
    colnames(all_pval) <- c("pval", "clust_num")
    
    curRow <- 1
    
    for(i in 1:num_clusters){
      curCluster <- res_hcpc_th_prob1$desc.var$quanti[[i]]
      all_pval[curRow:(curRow+num_features-1),"pval"] <- curCluster[,"p.value"]
      all_pval[curRow:(curRow+num_features-1),"clust_num"] <- i
      curRow <- curRow + num_features
    }
    all_pval$pval_adj <- p.adjust(all_pval$pval, method="BY")
    
    descr_var_hcpc <- res_hcpc_th_prob1$desc.var$quanti
    
    for(i in 1:num_clusters){
      
      descr_var_hcpc[[i]][,"p.value"] <- all_pval[which(all_pval$clust_num == i), "pval_adj"]
      descr_var_hcpc[[i]] <- descr_var_hcpc[[i]][which(descr_var_hcpc[[i]][,"p.value"] < 0.05),]
    }
    
    return(descr_var_hcpc)
}

#' Plot the association between the features and the clusters
#' @param descr_var_hcpc: all statisticially significant associations (p < 0.05) after multiple testing correction.
#' @param feature_names: names of the features
#' @param cluster_labels: labels of the clusters
#' @param resultsDir: directory to store the resulting plot  
plotAssociationOfFeatures2Clusters <- function(descr_var_hcpc, feature_names, cluster_labels, resultsDir, plotname)
{
  num_clusters <- length(descr_var_hcpc)
  
  features2vtestVal <- as.data.frame(matrix(data=0,nrow=length(feature_names), ncol=num_clusters))
  rownames(features2vtestVal) <- feature_names
  
  colnames(features2vtestVal) <- cluster_labels
  
  for(i in 1:num_clusters)
  {
    cur_vtestValues <- descr_var_hcpc[[i]][,"v.test"]
    features2vtestVal[names(cur_vtestValues),cluster_labels[i]] <- cur_vtestValues
  }
  
  features2vtestVal$feature_names <- rownames(features2vtestVal)
  features2vtestVal_m <- melt(features2vtestVal)
  colnames(features2vtestVal_m) <- c("feature_names", "cluster", "vtest")
  
  features2vtestVal_m$feature_names <- factor(as.character(features2vtestVal_m$feature_names), levels=rev(rownames(features2vtestVal)))
  features2vtestVal_m$cluster <- factor(as.character(features2vtestVal_m$cluster))
  
  # force white color for non-significant vtest values
  features2vtestVal_m[which(features2vtestVal_m$vtest == 0),"vtest"]  <- NA
  
  p_heat <- ggplot(features2vtestVal_m, aes(x=cluster,y=feature_names, fill = vtest)) + 
    geom_tile(color="black") +
    geom_text(aes(label = round(vtest, 1)), size=2) +
    scale_fill_gradient2(low="black", high="grey90", na.value = "white", mid="grey45") +
    scale_x_discrete(expand=c(0,0)) +
    theme(axis.title=element_blank())
  ggsave(file=paste(resultsDir, "/", plotname, ".pdf", sep=""))
}