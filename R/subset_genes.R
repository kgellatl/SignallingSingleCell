#' Select Genes
#'
#' This will select genes based on minimum expression and coefficient of variation, or by a preliminary PCA that selects genes based on malhanobis distance from the center.
#'
#' @param input the input data matrix.
#' @param threshold UMI threshold for gene detection
#' @param minCells number of cells expressed above threshold for a given gene
#' @param method can either be "Expression", CV", or "PCA"
#' @param nComp if method = PCA, the number of components to keep
#' @param cutoff the percentile of genes by coefficient of variation or PCA loadings to keep

#' @export
#' @details
#' This selects genes.
#' @examples
#' gene_subset <- subset_genes(input = exprs(ex_sc_example), method = "PCA", threshold = 3, minCells = 30, nComp = 15, cutoff = 0.75)

subset_genes <- function(input, method, threshold, minCells, nComp, cutoff){
  input <- exprs(input)
  if(method == "Expression"){
    gCount <- apply(input,1,function(x) length(which(x>=threshold)))
    gene_subset <- rownames(input[(which(gCount >= minCells)),])
  }
  if(method == "CV"){
    gCount <- apply(input,1,function(x) length(which(x>=threshold)))
    gene_subset <- rownames(input[(which(gCount >= minCells)),])
    g_exp <- log2(input[gene_subset,]+2)-1
    gmeans <- apply(g_exp,1,mean)
    gsd <- apply(g_exp,1,sd)
    CV <- gsd / gmeans
    CV <- CV[order(CV)]
    gene_subset <- names(CV[round(length(CV)*cutoff):length(CV)])
  }
  if(method == "PCA"){
    gCount <- apply(input,1,function(x) length(which(x>=threshold)))
    gene_subset <- rownames(input[(which(gCount >= minCells)),])
    input_scale <- scale(log2(input[gene_subset,]+2)-1)
    pc <- irlba::prcomp_irlba(t(input_scale), nComp, center = F)
    rownames(pc$rotation) <- gene_subset
    d <- mahalanobis(pc$rotation[,1:nComp], center=rep(0, nComp), cov = cov(pc$rotation[,1:nComp]))
    dThresh <- quantile(d,cutoff)
    gene_subset <- names(which(d>dThresh))
  }
  return(gene_subset)
}



