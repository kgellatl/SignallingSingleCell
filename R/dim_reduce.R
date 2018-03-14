#' Initial Dimension Reduction
#'
#' This function will do a preliminary dimension reduction. This is necessary because normalization across cell types requires preliminary clusters.
#'
#' @param input the input data matrix
#' @param threshold UMI threshold for gene detection
#' @param minCells number of cells expressed above threshold for a given gene
#' @param nComp_ICA the number of ICA components to reduce too, 5-20 recommended
#' @param tSNE_perp number of cells expressed above threshold for a given gene, 10-100 recommended
#' @export
#' @details
#' Before normalization preliminary clustering is needed in order to normalize internally within a cluster before normalizing across clusters. The minCells will be highly dependent on the number of cells in your dataset. In general aim to include genes expressed in no less than 1:100 cells, unless you have extremely rare cell types below that threshold.
#' @examples
#' filtered_data <- initial_dim_reduce(input = filtered_data, threshold = 5, minCells = 50, nComp_ICA = 10, tSNE_perp = 30, print_progress=TRUE)

dim_reduce <- function(input, threshold, minCells, nComp_ICA, tSNE_perp, print_progress=FALSE){
  gCount <- apply(input,1,function(x) length(which(x>=threshold)))
  gene_subset <- rownames(input[(which(gCount >= minCells)),])
  ic <- fastICA::fastICA(t(log2(input[gene_subset,]+2)-1), n.comp=(nComp_ICA), alg.typ = 'parallel', fun='logcosh', alpha = 1.0, method = 'C', verbose = print_progress)
  colnames(ic$A) <- gene_subset
  rownames(ic$S) <- colnames(input)
  set.seed(100)
  tSNE_result <- Rtsne::Rtsne(ic$S, dims = 2, perplexity = tSNE_perp, theta = 0.5, check_duplicates = F, pca = F, max_iter = 1000, verbose = print_progress)
  tSNE_result <- tSNE_result$Y
  row.names(tSNE_result) <- rownames(ic$S)
  colnames(ic$S) <- paste0("ICA Component ", seq(1:ncol(ic$S)))
  colnames(tSNE_result) <- c("x", "y")
  tSNE_result[,"x"] <-  abs(min(tSNE_result[,"x"]))+tSNE_result[,"x"]
  tSNE_result[,"x"] <-  tSNE_result[,"x"]/max(tSNE_result[,"x"])
  tSNE_result[,"y"] <-  abs(min(tSNE_result[,"y"]))+tSNE_result[,"y"]
  tSNE_result[,"y"] <-  tSNE_result[,"y"]/max(tSNE_result[,"y"])
  dim_reduce <- cbind(tSNE_result, ic$S)
  return(dim_reduce)
}
