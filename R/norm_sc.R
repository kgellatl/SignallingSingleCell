#' Normalization
#'
#' This function will perform normalization of your data
#'
#' @param input the input ex_sc
#' @param gene_frac the fraction of cells expressing a given gene to be included in normalization
#' @param gene_var the percentile of least variable genes to keep (ie 0.75 removes the 25 percent most variable genes)
#' @param sf_keep size factors can be greatly skewed in some cells. This will filter (based on z score thresholding of the size factors)
#' cells whose size factor is an outlier.
#' @param genelist the subset of genes to use for calculating size factors. Defaults to null. Provide to overrule internal gene selection.
#' @param norm_by the pdata variable on which to perform internal normalization before normalizing across this variable.
#' @param pool_size A vector of sizes for each pool for the normalization method
#' @param positive enforces positive size factors
#' @export
#' @details
#' If the method is ICA, independent component analysis will be performed, and then tSNE will do the final dimension reduction. If PCA is selected, PCA will be performed before on the expression matrix transpose before tSNE. This PCA will use the cells positions on the principal components. If iPCA is selected, PCA will be be performed but without transposing the data. This will create "meta cells" instead of meta genes created in the typical PCA. Then tSNE will be performed on each cells contribution (loading) to the meta cell. We find that iPCA is much more robust and leads to cleaner clusters than traditional PCA.
#' @examples
#' ex_sc_example <- dim_reduce(input = ex_sc_example, genelist = gene_subset, pre_reduce = "iPCA", nComp = 15, tSNE_perp = 30, iterations = 500, print_progress=TRUE)
#'
norm_sc <- function(input, gene_frac = 0.25, gene_var = 0.75, genelist = NULL, norm_by = "Cluster", pool_sizes = c(20,30,40,50), positive = TRUE){

  clusters <- pData(input)[,norm_by]

  if(is.null(genelist)){
    frac_mat <-  input
    exprs(frac_mat)[which(exprs(frac_mat) > 0)] <- 1

    gene_set <- c()

    for (i in 1:length(unique(clusters))) {
      int <- unique(clusters)[i]
      ind <- grep(int, clusters)
      frac <- apply(exprs(frac_mat)[,ind], 1, mean)
      set <- names(which(frac > gene_frac))
      gene_set <- c(gene_set, set)
    }

    tab <- table(gene_set)
    genelist <- names(tab)[which(tab == length(unique(clusters)))]
    input2 <- input[genelist,]
    gene_subset_var <- subset_genes(input2, method = "PCA", threshold = 0, minCells = 0, nComp = 10, cutoff = gene_var) # filter genes on variability
    stable_genes <- genelist[!genelist %in% gene_subset_var] # remove variable genes
    genelist <- stable_genes
    csums <- apply(exprs(input)[stable_genes,], 2,sum)
    zero_csums <- which(csums == 0)
    if(length(zero_csums) > 0)
      err_result <- vector(mode = "list", length = 2)
    err_result[[1]] <- csums
    err_result[[2]] <- stable_genes
    warning("With provided parameters, some cells have zero expression. Try reducing gene_frac argument. See returned result")
    return(err_result)
  }

  SCE <- SingleCellExperiment::SingleCellExperiment(list(counts = exprs(input)))
  SCE <- scran::computeSumFactors(SCE, cluster=clusters, subset.row = genelist, sizes=pool_sizes, positive = positive, min.mean = 0)
  SCE <-  scater::normalise(SCE)
  norm_counts <- exprs(SCE)
  norm_counts[!is.finite(norm_counts)] <- 0
  norm_counts <- norm_counts[apply(norm_counts, 1, function(x) !all(x==0)),]
  norm_counts <- norm_counts[,colSums(norm_counts)>0]
  input_norm <- construct_ex_sc(norm_counts)

  ind <- match(colnames(input_norm), colnames(SCE))
  pData(input_norm) <- pData(input)[ind,]

  size_factor <- sizeFactors(SCE)
  size_factor <- size_factor[ind]
  pData(input_norm)$size_factor <- size_factor
  ind <- match(stable_genes, rownames(input_norm))

  fData(input_norm)$norm_gene <- 0

  fData(input_norm)$norm_gene[ind] <- "Yes"
  fData(input_norm)$norm_gene[-ind] <- "No"

  return(input_norm)
}
