#' Cluster_sc
#' This will perform clustering on your single cell data.
#'
#' @param input the input data matrix.
#' @param dimension either "Comp" or "2d"
#' @param method can either be "spectral" or "density"
#' @param num_clust the number of clusters
#' @export
#' @details
#' This will perform clustering on either the high dimensional PCA / ICA components or on the final 2d result from tsne
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
    }
  }
  return(input)
}



