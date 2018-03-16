#' Select Genes
#'
#' This will select genes based on minimum expression and coefficient of variation, or by a preliminary PCA that selects genes based on malhanobis distance from the center.
#'
#' @param input the input data matrix.
#' @param threshold UMI threshold for gene detection
#' @param minCells number of cells expressed above threshold for a given gene
#' @export
#' @details
#' This selects genes.
#' @examples
#' gene_subset <- subset_genes(exprs(ex_sc_example), threshold = 1, minCells = 10)

subset_genes <- function(input, threshold, minCells){
  gCount <- apply(input,1,function(x) length(which(x>=threshold)))
  gene_subset <- rownames(input[(which(gCount >= minCells)),])
}
