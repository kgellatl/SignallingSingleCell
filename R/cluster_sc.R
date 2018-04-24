#' Cluster Single Cell
#'
#' This will perform clustering on your single cell data.
#'
#' @param input the input ex_sc
#' @param dimension either "Comp" or "2d"
#' @param method can either be "spectral" or "density" which is on 2d
#' @param num_clust the number of clusters
#' @export
#' @details
#' This will perform clustering on either the high dimensional PCA / ICA components if dimension = Comp,
#' or the 2d tsne result if method = density. Typically spectral clustering works much better on higher dimensional data,
#' which density based clustering works better on 2d data.
#' @examples
#' ex_sc_example <- cluster_sc(input = ex_sc_example, dimension = "Comp", method = "spectral", num_clust = 6)

cluster_sc <- function(input, dimension, method, num_clust){
  if(dimension == "Comp"){
    if(method == "spectral"){
      spec <- kknn::specClust(pData(input)[,grep("Comp", colnames(pData(input)))], centers = num_clust, method = 'random-walk')
      cluster <- spec$cluster
      cluster <- paste0("Cluster", cluster)
      pData(input)$Cluster <- cluster
    }
    if(method == "density"){
      dist_tSNE = dist(pData(input)[,grep("Comp", colnames(pData(input)))])
      clust_tSNE = densityClust::densityClust(dist_tSNE, gaussian = T) # calculate density
      comb = as.data.frame(clust_tSNE$rho*clust_tSNE$delta) # combine rho and delta values
      comb = comb[order(comb[,1], decreasing = T), ,drop=F] # order cells by highest rho*delta
      cellcut = rownames(comb)[1:num_clust] # select top num_clust cells to be cluster centers
      cellidx = which(names(clust_tSNE$rho)%in%cellcut) # get index for cluster centers
      clust_tSNE = densityClust::findClusters(clust_tSNE, peaks = cellidx) # define clusters
      pData(input)$Cluster = paste0("Cluster", clust_tSNE$clusters)
    }
  }
  if(dimension == "2d"){
    if(method == "spectral"){
      spec <- kknn::specClust(pData(input)[,c("x", "y")], centers = num_clust, method = 'random-walk')
      cluster <- spec$cluster
      cluster <- paste0("Cluster", cluster)
      pData(input)$Cluster <- cluster
    }
    if(method == "density"){
      dist_tSNE = dist(pData(input)[,c("x","y")]) # select tSNE coordinates
      clust_tSNE = densityClust::densityClust(dist_tSNE, gaussian = T) # calculate density
      comb = as.data.frame(clust_tSNE$rho*clust_tSNE$delta) # combine rho and delta values
      comb = comb[order(comb[,1], decreasing = T), ,drop=F] # order cells by highest rho*delta
      cellcut = rownames(comb)[1:num_clust] # select top num_clust cells to be cluster centers
      cellidx = which(names(clust_tSNE$rho)%in%cellcut) # get index for cluster centers
      clust_tSNE = densityClust::findClusters(clust_tSNE, peaks = cellidx) # define clusters
      pData(input)$Cluster = paste0("Cluster", clust_tSNE$clusters)
    }
  }
  return(input)
}
