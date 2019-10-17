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

cell_net <- function(input, genelist, n_blocks = 5, number_edges = 5, method = "bigcor") {

  #############  #############  #############  #############
  #############  BIGCOR section
  #############  #############  #############  #############
  if (method == "bigcor"){

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
        COR <- cor(x[, G1], x[, G2], method = "spearman")
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

    # cor_mat <- HiClimR::fastCor(x = mat, nSplit =  n_blocks)
    cor_mat <- bigcor(x = mat, nblocks =  n_blocks)

    num_cutoff <- number_edges
    num_edges <- matrix(ncol = 3, nrow = ncol(mat)*num_cutoff)
    positions <- seq(from = 1, to = nrow(num_edges), by = num_cutoff)
    for (i in 1:ncol(cor_mat)) {
      int_cell <- colnames(cor_mat)[i]
      int_vec <- cor_mat[,i]
      int_vec <- rev(sort(int_vec))
      int_vec <- int_vec[2:length(int_vec)]
      num_cells <- int_vec[1:num_cutoff]
      c1 <- rep(int_cell, length(num_cells))
      c2 <- names(num_cells)
      c3 <- as.vector(num_cells)
      start <- positions[i]
      end <- positions[i]+(num_cutoff-1)

      num_edges[start:end,1] <- c1
      num_edges[start:end,2] <- c2
      num_edges[start:end,3] <- c3
    }

    dim(num_edges)
    num_edges <- as.data.frame(num_edges)
    num_edges[,4] <- as.numeric(as.character((num_edges[,3])))
    vals <- as.vector(scale(num_edges[,4]))
    vals <- (abs(min(vals))+vals)+1
    num_graph <- graph_from_data_frame(num_edges)
    E(num_graph)$weight <- vals

  }

  #############  #############  #############  #############
  #############  truncated section
  #############  #############  #############  #############

  if (method == "truncated") {



  }

  return(num_graph)

}

