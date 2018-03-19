#' Dimension Reduction
#'
#' This function will do dimensionality reduction.
#'
#' @param input the input data matrix.
#' @param pre_reduce the algorithm choice for reduction before tSNE (either "ICA", "PCA", "iPCA").
#' @param nComp the number of components to reduce too before tSNE, 5-20 recommended.
#' @param tSNE_perp number of cells expressed above threshold for a given gene, 10-100 recommended.
#' @param iterations The number of iterations for tSNE to perform.
#' @param print_progress Print to the terminal progress information.
#' @export
#' @details
#' If the method is ICA, independent component analysis will be performed, and then tSNE will do the final dimension reduction. If PCA is selected, PCA will be performed before on the expression matrix transpose before tSNE. This PCA will use the cells positions on the principal components. If iPCA is selected, PCA will be be performed but without transposing the data. This will create "meta cells" instead of meta genes created in the typical PCA. Then tSNE will be performed on each cells contribution (loading) to the meta cell. We find that iPCA is much more robust and leads to cleaner clusters than traditional PCA.
#' @examples
#' ex_sc_example <- dim_reduce(input = ex_sc_example, genelist = gene_subset, pre_reduce = "iPCA", nComp = 15, tSNE_perp = 30, iterations = 500, print_progress=TRUE)
#'
dim_reduce <- function(input, genelist = gene_subset, pre_reduce = "iPCA", nComp = 15, tSNE_perp = 30, iterations = 1000, print_progress=TRUE){
  input_exp <- exprs(input)[genelist,]
  input_scale <- scale(log2(input_exp[,]+2)-1)
  check <- grep("Comp", colnames(pData(input)))
  if(length(check) > 0){
    pData(input) <- pData(input)[,-check]
  }
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
    tSNE_result <- Rtsne::Rtsne(ica$S, dims = 2, perplexity = tSNE_perp, theta = 0.5, check_duplicates = F, pca = F, max_iter = iterations, verbose = print_progress)
    tSNE_result <- tSNE_result$Y
    row.names(tSNE_result) <- rownames(ica$S)
    colnames(ica$S) <- paste0("IC-Comp", seq(1:ncol(ica$S)))
    colnames(tSNE_result) <- c("x", "y")
    tSNE_result[,"x"] <-  abs(min(tSNE_result[,"x"]))+tSNE_result[,"x"]
    tSNE_result[,"x"] <-  tSNE_result[,"x"]/max(tSNE_result[,"x"])
    tSNE_result[,"y"] <-  abs(min(tSNE_result[,"y"]))+tSNE_result[,"y"]
    tSNE_result[,"y"] <-  tSNE_result[,"y"]/max(tSNE_result[,"y"])
    prelim_dims <- cbind(tSNE_result, ica$S)
  }
  if(pre_reduce == "PCA"){
    if(print_progress == TRUE){
      print("Starting PCA")
    }
    PCA <- irlba::prcomp_irlba(t(input_scale), nComp, center = F)
    rownames(PCA$x) <- colnames(input)
    colnames(PCA$x) <- paste0("PC-Comp", seq(1:ncol(PCA$x)))
    set.seed(100)
    if(print_progress == TRUE){
      print("Starting tSNE")
    }
    tSNE_result <- Rtsne::Rtsne(PCA$x, dims = 2, perplexity = tSNE_perp, theta = 0.5, check_duplicates = F, pca = F, max_iter = iterations, verbose = print_progress)
    tSNE_result <- tSNE_result$Y
    row.names(tSNE_result) <- colnames(input)
    colnames(tSNE_result) <- c("x", "y")
    tSNE_result[,"x"] <-  abs(min(tSNE_result[,"x"]))+tSNE_result[,"x"]
    tSNE_result[,"x"] <-  tSNE_result[,"x"]/max(tSNE_result[,"x"])
    tSNE_result[,"y"] <-  abs(min(tSNE_result[,"y"]))+tSNE_result[,"y"]
    tSNE_result[,"y"] <-  tSNE_result[,"y"]/max(tSNE_result[,"y"])
    prelim_dims <- cbind(tSNE_result, PCA$x)
  }
  if(pre_reduce == "iPCA"){
    if(print_progress == TRUE){
      print("Starting iPCA")
    }
    iPCA <- irlba::prcomp_irlba(input_scale, nComp, center = F)
    rownames(iPCA$rotation) <- colnames(input)
    colnames(iPCA$rotation) <- paste0("iPC-Comp", seq(1:ncol(iPCA$rotation)))
    set.seed(100)
    if(print_progress == TRUE){
      print("Starting tSNE")
    }
    tSNE_result <- Rtsne::Rtsne(iPCA$rotation, dims = 2, perplexity = tSNE_perp, theta = 0.5, check_duplicates = F, pca = F, max_iter = iterations, verbose = print_progress)
    tSNE_result <- tSNE_result$Y
    row.names(tSNE_result) <- colnames(input)
    colnames(tSNE_result) <- c("x", "y")
    tSNE_result[,"x"] <-  abs(min(tSNE_result[,"x"]))+tSNE_result[,"x"]
    tSNE_result[,"x"] <-  tSNE_result[,"x"]/max(tSNE_result[,"x"])
    tSNE_result[,"y"] <-  abs(min(tSNE_result[,"y"]))+tSNE_result[,"y"]
    tSNE_result[,"y"] <-  tSNE_result[,"y"]/max(tSNE_result[,"y"])
    prelim_dims <- cbind(tSNE_result, iPCA$rotation)
  }
  pData(input)$x <- prelim_dims[,"x"]
  pData(input)$y <- prelim_dims[,"y"]
  index <- grep("Comp", colnames(prelim_dims))
  pData(input) <- cbind(pData(input), prelim_dims[,index])
  return(input)
}

