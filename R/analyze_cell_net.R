#' Dimension Reduction
#'
#' This function will do dimensionality reduction.
#'
#' @param input the input ex_sc
#' @param genelist the subset of genes to perform dimensionality reduction on
#' @param pre_reduce the algorithm choice for reduction before tSNE (either "ICA", "PCA", "iPCA").
#' @param nComp the number of components to reduce too before tSNE, 5-20 recommended.
#' @param tSNE_perp number of cells expressed above threshold for a given gene, 10-100 recommended.
#' @param iterations The number of iterations for tSNE to perform.
#' @param print_progress will print progress if TRUE
#' @param nVar cutoff for percent of variance explained from PCs
#' @importFrom fastICA fastICA
#' @importFrom  Rtsne Rtsne
#' @importFrom irlba prcomp_irlba
#' @export
#' @details
#' If the method is ICA, independent component analysis will be performed, and then tSNE will do the final dimension reduction. If PCA is selected, PCA will be performed before on the expression matrix transpose before tSNE. This PCA will use the cells positions on the principal components. If iPCA is selected, PCA will be be performed but without transposing the data. This will create "meta cells" instead of meta genes created in the typical PCA. Then tSNE will be performed on each cells contribution (loading) to the meta cell. We find that iPCA is much more robust and leads to cleaner clusters than traditional PCA.
#' @examples
#' ex_sc_example <- dim_reduce(input = ex_sc_example, genelist = gene_subset, pre_reduce = "iPCA", nComp = 15, tSNE_perp = 30, iterations = 500, print_progress=TRUE)
#'

analyze_cell_net <- function(input_ex_sc, input_cell_net, verbose = T, layout = "fr", cluster_method = "eigen", crossing_edge_weight = 1, internal_edge_weight = 10){

  if(layout == "fr"){
    if(verbose){
      print("Calculating fr layout")
    }
    l <- layout_with_fr(input_cell_net)
    input_ex_sc$fr_x <- l[,1]
    input_ex_sc$fr_y <- l[,2]
  }

  if(layout == "kk"){
  if(verbose){
    print("Calculating kk layout")
  }

  l <- layout_with_kk(input_cell_net)
  input_ex_sc$kk_x <- l[,1]
  input_ex_sc$kk_y <- l[,2]
  }

  if(cluster_method == "eigen"){
    if(verbose){
      print("Calculating leading eigen clusters")
    }
    clusters <- cluster_leading_eigen(as.undirected(input_cell_net))
    weights_clusters <- ifelse(crossing(clusters, input_cell_net), crossing_edge_weight, internal_edge_weight)
    cluster_weighted_layout <- layout_with_fr(input_cell_net, weights = weights_clusters)
    pData(input_ex_sc)$eigen_x <- cluster_weighted_layout[,1]
    pData(input_ex_sc)$eigen_y <- cluster_weighted_layout[,2]
    input_ex_sc$leading_eigen_clusters <- paste0("leading_eigen_", clusters$membership)
  }

  if(cluster_method == "greedy"){
    if(verbose){
      print("Calculating fast greedy clusters")
    }
    clusters <- cluster_fast_greedy(as.undirected(input_cell_net))
    weights_clusters <- ifelse(crossing(clusters, input_cell_net), crossing_edge_weight, internal_edge_weight)
    cluster_weighted_layout <- layout_with_fr(input_cell_net, weights = weights_clusters)
    pData(input_ex_sc)$greedy_x <- cluster_weighted_layout[,1]
    pData(input_ex_sc)$greedy_y <- cluster_weighted_layout[,2]
    input_ex_sc$fast_greedy_clusters <- paste0("fast_greedy_", clusters$membership)
  }


  if(cluster_method == "louvain"){
    if(verbose){
      print("Calculating louvain clusters")
    }

    clusters <- cluster_louvain(as.undirected(input_cell_net))
    weights_clusters <- ifelse(crossing(clusters, input_cell_net), crossing_edge_weight, internal_edge_weight)
    cluster_weighted_layout <- layout_with_fr(input_cell_net, weights = weights_clusters)
    pData(input_ex_sc)$louvain_x <- cluster_weighted_layout[,1]
    pData(input_ex_sc)$louvain_y <- cluster_weighted_layout[,2]
    input_ex_sc$louvain_clusters <- paste0("louvain_", clusters$membership)
  }

  return(input_ex_sc)
}

