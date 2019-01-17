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
dim_reduce <- function(input, genelist = gene_subset, pre_reduce = "iPCA", nComp = 15, tSNE_perp = 30, iterations = 1000, print_progress=TRUE, nVar=NA, log = T, scale = T){
  input_exp <- exprs(input)[genelist,]
  if(log){
    input_exp <- log2(input_exp[,]+2)-1
  }
  if(scale){
    input_exp <- scale(input_exp)
  }
  check <- grep("Comp", colnames(pData(input)))
  if(length(check) > 0){
    pData(input) <- pData(input)[,-check]
  }
  if(pre_reduce == "ICA"){
    if(print_progress == TRUE){
      print("Starting ICA")
    }
    ica <- fastICA::fastICA(t(input_exp), n.comp=(nComp), alg.typ = 'parallel', fun='logcosh', alpha = 1.0, method = 'C', verbose = print_progress)
    colnames(ica$A) <- gene_subset
    rownames(ica$S) <- colnames(input)
    set.seed(100)
    if(print_progress == TRUE){
      print("Starting tSNE")
    }
    tsne_input = ica$S
    colnames(ica$S) <- paste0("IC_Comp", seq(1:ncol(ica$S)))
  }
  if(pre_reduce == "PCA"){
    if(print_progress == TRUE){
      print("Starting PCA")
    }
    PCA <- irlba::prcomp_irlba(t(input_exp), nComp, center = F)
    rownames(PCA$x) <- colnames(input)
    colnames(PCA$x) <- paste0("PC_Comp", seq(1:ncol(PCA$x)))
    set.seed(100)
    if(print_progress == TRUE){
      print("Starting tSNE")
    }
    tsne_input = PCA$x
  }
  if(pre_reduce == "iPCA"){
    if(print_progress == TRUE){
      print("Starting iPCA")
    }
    iPCA <- irlba::prcomp_irlba(input_exp, nComp, center = F)
    rownames(iPCA$rotation) <- colnames(input)
    colnames(iPCA$rotation) <- paste0("iPC_Comp", seq(1:ncol(iPCA$rotation)))
    set.seed(100)
    if(print_progress == TRUE){
      print("Starting tSNE")
    }
    tsne_input = iPCA$rotation
  }
  if(pre_reduce == "vPCA"){
    if(print_progress == TRUE){
      print("Starting vPCA")
    }
    vPCA <- irlba::prcomp_irlba(input_exp, n = nComp)
    # sum components until variance is >= x%
    var = vPCA$sdev^2/sum(vPCA$sdev^2)
    totalvar = var[1]
    maxPC = 1
    while (totalvar < nVar) {
      maxPC=maxPC+1
      totalvar = sum(var[1:maxPC])
    }
    vPCA$rotation = vPCA$rotation[,1:maxPC]
    rownames(vPCA$rotation) <- colnames(input)
    colnames(vPCA$rotation) <- paste0("iPC_Comp", seq(1:ncol(vPCA$rotation)))
    set.seed(100)
    if(print_progress == TRUE){
      print("Starting tSNE")
    }
    tsne_input = vPCA$rotation
  }
  tSNE_result <- Rtsne::Rtsne(tsne_input, dims = 2, perplexity = tSNE_perp, theta = 0.5, check_duplicates = F, pca = F, max_iter = iterations, verbose = print_progress)
  tSNE_result <- tSNE_result$Y
  row.names(tSNE_result) <- rownames(tsne_input)
  colnames(tSNE_result) <- c("x", "y")
  tSNE_result[,"x"] <-  abs(min(tSNE_result[,"x"]))+tSNE_result[,"x"]
  tSNE_result[,"x"] <-  tSNE_result[,"x"]/max(tSNE_result[,"x"])
  tSNE_result[,"y"] <-  abs(min(tSNE_result[,"y"]))+tSNE_result[,"y"]
  tSNE_result[,"y"] <-  tSNE_result[,"y"]/max(tSNE_result[,"y"])
  prelim_dims <- cbind(tSNE_result, tsne_input)
  pData(input)$x <- prelim_dims[,"x"]
  pData(input)$y <- prelim_dims[,"y"]
  index <- grep("Comp", colnames(prelim_dims))
  pData(input) <- cbind(pData(input), prelim_dims[,index])
  return(input)
}

