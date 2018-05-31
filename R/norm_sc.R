#' Normalization
#'
#' This function will perform normalization of your data
#'
#' @param input the input ex_sc
#' @param genelist the subset of genes to use for calculating size factors.
#' @param pool_size A vector of sizes for each pool for the normalization method
#' @param positive enforces positive size factors
#' @export
#' @details
#' If the method is ICA, independent component analysis will be performed, and then tSNE will do the final dimension reduction. If PCA is selected, PCA will be performed before on the expression matrix transpose before tSNE. This PCA will use the cells positions on the principal components. If iPCA is selected, PCA will be be performed but without transposing the data. This will create "meta cells" instead of meta genes created in the typical PCA. Then tSNE will be performed on each cells contribution (loading) to the meta cell. We find that iPCA is much more robust and leads to cleaner clusters than traditional PCA.
#' @examples
#' ex_sc_example <- dim_reduce(input = ex_sc_example, genelist = gene_subset, pre_reduce = "iPCA", nComp = 15, tSNE_perp = 30, iterations = 500, print_progress=TRUE)
#'
norm_sc <- function(input, genelist = gene_subset, pool_sizes = c(20,30,40,50), positive = TRUE){
  clusters <- pData(input)$Cluster
  larger.sce <- SingleCellExperiment::SingleCellExperiment(list(counts = exprs(input)))
  larger.sce <- scran::computeSumFactors(larger.sce, cluster=clusters, subset.row = genelist, sizes=pool_sizes, positive = positive)
  larger.sce <-  scater::normalise(larger.sce)
  norm_counts <- exprs(larger.sce)
  norm_counts[!is.finite(norm_counts)] <- 0
  norm_counts <- norm_counts[apply(norm_counts, 1, function(x) !all(x==0)),]
  norm_counts <- norm_counts[,colSums(norm_counts)>0]
  input_norm <- construct_ex_sc(norm_counts)
  ind <- match(colnames(input_norm), colnames(larger.sce))
  size_factor <- sizeFactors(larger.sce)
  size_factor <- size_factor[ind]
  pData(input_norm)$size_factor <- size_factor
  return(input_norm)
}
