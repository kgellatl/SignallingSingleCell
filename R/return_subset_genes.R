#' Subset Genes
#'
#' This will select genes based on minimum expression, coefficient of variation,
#' or by a preliminary PCA.
#'
#' @param input the input ex_sc.
#' @param threshold UMI threshold for gene detection
#' @param minCells number of cells expressed above threshold for a given gene
#' @param method can either be "Expression", CV", or "PCA".
#' @param nComp if method = PCA, the number of components to keep
#' @param cutoff the percentile of genes to keep

#' @export
#' @details
#' Genes will be first filtered by minimum expression selecting by subsetting to genes that are
#' expressed above the threshold in more than minCells.
#' If the method is CV, it will first subset the genes based on the expression cutoffs,
#' then find the coefficient of variation across all genes.
#' Next it will select the percentile of genes (cutoff) based on their coefficient of variation. The last method will perform PCA on the cells, and then look at the loadings of each gene. By finding genes that are off center (via malhanoobis distance) we can filter to include only genes that contribute significant variance to the data.
#' @examples
#' gene_subset <- subset_genes(input = exprs(ex_sc_example), method = "PCA", threshold = 3, minCells = 30, nComp = 15, cutoff = 0.75)

return_subset_genes <- function(input, method, cutoff = NULL, num_genes = NULL){
  mat <- fData(input)
  if(method == "CV"){
    if(is.null(cutoff) && is.null(num_genes)){
      gene_set <- rownames(mat)[which(mat$CV_selected == T)]
    } else {
      if(!is.null(num_genes)){
        mat <- mat[rev(order(mat$CV)),]
        gene_set <- rownames(mat)[1:num_genes]
      }
      if(!is.null(cutoff)){
        val <- quantile(mat$CV, cutoff)
        ind <- which(mat$CV > val)
        gene_set <- rownames(mat)[ind]

      }
    }

  }

  if(method == "PCA"){
    if(is.null(cutoff) && is.null(num_genes)){
      gene_set <- rownames(mat)[which(mat$malhanobis_selected == T)]
    } else {
      if(!is.null(num_genes)){
        mat <- mat[rev(order(mat$malhanobis_d)),]
        gene_set <- rownames(mat)[1:num_genes]
      }
      if(!is.null(cutoff)){
        val <- quantile(mat$malhanobis_d, cutoff)
        ind <- which(mat$malhanobis_d > val)
        gene_set <- rownames(mat)[ind]

      }
    }
  }

  return(gene_set)
}



