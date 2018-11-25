#' Normalization
#'
#' This function will perform normalization of your data
#'
#' @param input the input ex_sc
#' @param gene_set the fraction of cells expressing a given gene to be included in normalization
#' @param gene_names the fraction of cells expressing a given gene to be included in normalization
#' @export
#' @details
#' If the method is ICA, independent component analysis will be performed, and then tSNE will do the final dimension reduction. If PCA is selected, PCA will be performed before on the expression matrix transpose before tSNE. This PCA will use the cells positions on the principal components. If iPCA is selected, PCA will be be performed but without transposing the data. This will create "meta cells" instead of meta genes created in the typical PCA. Then tSNE will be performed on each cells contribution (loading) to the meta cell. We find that iPCA is much more robust and leads to cleaner clusters than traditional PCA.
#' @examples
#' ex_sc_example <- dim_reduce(input = ex_sc_example, genelist = gene_subset, pre_reduce = "iPCA", nComp = 15, tSNE_perp = 30, iterations = 500, print_progress=TRUE)
#'
agg_gene <- function(input, gene_set, gene_name = NA){
  if(is.na(gene_name)){
    gene_name <- paste0(gene_set, collapse = "_")
  }
  if(is.null(fData(input)$agg_gene) == TRUE){
    fData(input)$agg_gene <- FALSE
  }
  tmp <- exprs(input)[gene_set,]
  g_sum <- apply(tmp,2,sum)
  full <- exprs(input)
  full_bound <- rbind(full, g_sum)
  rownames(full_bound) <- c(rownames(full), gene_name)
  input_bound <- construct_ex_sc(full_bound)
  pData(input_bound) <- pData(input)
  fData(input_bound) <- rbind(fData(input), c(rep(NA, ncol(fData(input)))))
  rownames(fData(input_bound)) <- c(rownames(fData(input)), gene_name)
  fData(input_bound)[gene_name,"agg_gene"] <-  paste0(gene_set, collapse = "_")
  return(input_bound)
}
