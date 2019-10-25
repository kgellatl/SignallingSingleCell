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

cell_cor <- function(input, genelist, cor = "spearman", n_blocks = 5) {

    bigcor <- function(x, nblocks, verbose = TRUE, ...)
    {
      orig_names <- colnames(x)
      NCOL <- ncol(x)
      NCOL_original <- ncol(x)

      ## test if ncol(x) %% nblocks gives remainder 0
      if (NCOL %% nblocks != 0) {
        block_size <- floor(NCOL/nblocks)
        block_size*(nblocks+1)
        new_cols <- block_size - NCOL %% nblocks
        new_mat <- matrix(0, ncol = new_cols, nrow = nrow(x))
        x <- cbind(x, new_mat)
        nblocks <- nblocks+1
      }

      NCOL <- ncol(x)

      ## preallocate square matrix of dimension
      ## ncol(x) in 'ff' single format
      corMAT <- matrix(ncol = NCOL, nrow = NCOL)

      ## split column numbers into 'nblocks' groups
      SPLIT <- split(1:NCOL, rep(1:nblocks, each = NCOL/nblocks))

      ## create all unique combinations of blocks
      COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
      COMBS <- t(apply(COMBS, 1, sort))
      COMBS <- unique(COMBS)

      ## iterate through each block combination, calculate correlation matrix
      ## between blocks and store them in the preallocated matrix on both
      ## symmetric sides of the diagonal
      for (i in 1:nrow(COMBS)) {
        COMB <- COMBS[i, ]
        G1 <- SPLIT[[COMB[1]]]
        G2 <- SPLIT[[COMB[2]]]
        if (verbose) cat("Block", COMB[1], "with Block", COMB[2], "\n")
        flush.console()
        if(cor == "spearman"){
          COR <- cor(x[, G1], x[, G2], method = "spearman")
        }
        if(cor == "pearsons"){
          COR <- cor(x[, G1], x[, G2])
        }
        corMAT[G1, G2] <- COR
        corMAT[G2, G1] <- t(COR)
        COR <- NULL
      }
      gc()
      corMAT <- corMAT[1:NCOL_original,1:NCOL_original]
      colnames(corMAT) <- orig_names
      rownames(corMAT) <- orig_names
      return(corMAT)
    }

    mat <- exprs(input)[genelist,]
    cor_mat <- bigcor(x = mat, nblocks =  n_blocks)
    return(cor_mat)

}

