#' Initial Dimension Reduction
#'
#' This function will do a preliminary dimension reduction. This is necessary because normalization across cell types requires preliminary clusters.
#'
#' @param input the input data matrix
#' @param threshold UMI threshold for gene detection
#' @param pre_reduce the algorithm choice for reduction before tSNE
#' @param minCells number of cells expressed above threshold for a given gene
#' @param nComp the number of components to reduce too before tSNE, 5-20 recommended
#' @param tSNE_perp number of cells expressed above threshold for a given gene, 10-100 recommended
#' @export
#' @details
#' Before normalization preliminary clustering is needed in order to normalize internally within a cluster before normalizing across clusters.
#' The minCells will be highly dependent on the number of cells in your dataset. In general aim to include genes expressed in no less than 1:100 cells, unless you have extremely rare cell types below that threshold.
#' Possible pre_reduce algorithms are ICA, PCA, iPCA, and sPCA.
#' @examples
#' filtered_data <- dim_reduce(input = filtered_data, threshold = 5, minCells = 50, pre_reduce = "ICA", nComp = 10, tSNE_perp = 30, print_progress=TRUE)

dim_reduce <- function(input, threshold = 5, minCells = 50, pre_reduce = "ICA", nComp = 10, tSNE_perp = 30, print_progress=TRUE){
  gCount <- apply(input,1,function(x) length(which(x>=threshold)))
  gene_subset <- rownames(input[(which(gCount >= minCells)),])
  input_scale <- scale(log2(input[gene_subset,]+2)-1)
  if(pre_reduce == "ICA"){
    if(print_progress == TRUE){
      print("Starting ICA")
    }
    ica <- fastICA::fastICA(t(input_scale), n.comp=(nComp), alg.typ = 'parallel', fun='logcosh', alpha = 1.0, method = 'C', verbose = print_progress)
    colnames(ica$A) <- gene_subset
    rownames(ica$S) <- colnames(input)
    set.seed(100)
    if(print_progress == TRUE){
      print("Starting tSNE")
    }
    tSNE_result <- Rtsne::Rtsne(ica$S, dims = 2, perplexity = tSNE_perp, theta = 0.5, check_duplicates = F, pca = F, max_iter = 1000, verbose = print_progress)
    tSNE_result <- tSNE_result$Y
    row.names(tSNE_result) <- rownames(ica$S)
    colnames(ica$S) <- paste0("IC", seq(1:ncol(ica$S)))
    colnames(tSNE_result) <- c("x", "y")
    tSNE_result[,"x"] <-  abs(min(tSNE_result[,"x"]))+tSNE_result[,"x"]
    tSNE_result[,"x"] <-  tSNE_result[,"x"]/max(tSNE_result[,"x"])
    tSNE_result[,"y"] <-  abs(min(tSNE_result[,"y"]))+tSNE_result[,"y"]
    tSNE_result[,"y"] <-  tSNE_result[,"y"]/max(tSNE_result[,"y"])
    dim_reduce <- cbind(tSNE_result, ica$S)
  }
  if(pre_reduce == "PCA"){
    if(print_progress == TRUE){
      print("Starting PCA")
    }
    PCA <- irlba::prcomp_irlba(t(input_scale), nComp, center = F)
    rownames(PCA$x) <- colnames(input)
    set.seed(100)
    if(print_progress == TRUE){
      print("Starting tSNE")
    }
    tSNE_result <- Rtsne::Rtsne(PCA$x, dims = 2, perplexity = tSNE_perp, theta = 0.5, check_duplicates = F, pca = F, max_iter = 1000, verbose = print_progress)
    tSNE_result <- tSNE_result$Y
    row.names(tSNE_result) <- colnames(input)
    colnames(tSNE_result) <- c("x", "y")
    tSNE_result[,"x"] <-  abs(min(tSNE_result[,"x"]))+tSNE_result[,"x"]
    tSNE_result[,"x"] <-  tSNE_result[,"x"]/max(tSNE_result[,"x"])
    tSNE_result[,"y"] <-  abs(min(tSNE_result[,"y"]))+tSNE_result[,"y"]
    tSNE_result[,"y"] <-  tSNE_result[,"y"]/max(tSNE_result[,"y"])
    dim_reduce <- cbind(tSNE_result, PCA$x)
  }
  if(pre_reduce == "iPCA"){
    if(print_progress == TRUE){
      print("Starting iPCA")
    }
    iPCA <- irlba::prcomp_irlba(input_scale, nComp, center = F)
    rownames(iPCA$rotation) <- colnames(input)
    set.seed(100)
    if(print_progress == TRUE){
      print("Starting tSNE")
    }
    tSNE_result <- Rtsne::Rtsne(iPCA$rotation, dims = 2, perplexity = tSNE_perp, theta = 0.5, check_duplicates = F, pca = F, max_iter = 1000, verbose = print_progress)
    tSNE_result <- tSNE_result$Y
    row.names(tSNE_result) <- colnames(input)
    colnames(tSNE_result) <- c("x", "y")
    tSNE_result[,"x"] <-  abs(min(tSNE_result[,"x"]))+tSNE_result[,"x"]
    tSNE_result[,"x"] <-  tSNE_result[,"x"]/max(tSNE_result[,"x"])
    tSNE_result[,"y"] <-  abs(min(tSNE_result[,"y"]))+tSNE_result[,"y"]
    tSNE_result[,"y"] <-  tSNE_result[,"y"]/max(tSNE_result[,"y"])
    dim_reduce <- cbind(tSNE_result, iPCA$rotation)
  }
  return(dim_reduce)
}

