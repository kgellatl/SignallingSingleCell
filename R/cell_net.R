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

cell_net <- function(input_cor_matrix, number_edges = NULL, correlation_cutoff = NULL, min_edges = 10, max_edges = 100, weight = "standardized") {

  if(!is.null(number_edges)){
    num_cutoff <- number_edges
    num_edges <- matrix(ncol = 3, nrow = ncol(input_cor_matrix)*num_cutoff)
    num_edges <- as.data.frame(num_edges)

    positions <- seq(from = 1, to = nrow(num_edges), by = num_cutoff)

    alerts <- c()
    for (i in 1:20) {
      printi <- floor(nrow(input_cor_matrix)/20)*i
      alerts <- c(alerts, printi)
    }

    for (i in 1:ncol(input_cor_matrix)) {
      if(i %in% alerts){
        ind <- match(i, alerts)
        print(paste0(ind*5, "% Complete"))
      }
      int_cell <- colnames(input_cor_matrix)[i]
      int_vec <- input_cor_matrix[,i]
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
  }

  if(!is.null(correlation_cutoff)){
    num_cutoff <- correlation_cutoff
    num_edges <- matrix(ncol = 3, nrow = max_edges*ncol(input_cor_matrix))
    num_edges <- as.data.frame(num_edges)
    cells_below_cutoff <- 0
    cells_above_cutoff <- 0

    alerts <- c()
    for (i in 1:20) {
      printi <- floor(nrow(input_cor_matrix)/20)*i
      alerts <- c(alerts, printi)
    }


    start <- 1
    for (i in 1:ncol(input_cor_matrix)) {
      if(i %in% alerts){
        ind <- match(i, alerts)
        print(paste0(ind*5, "% Complete"))
      }
      int_cell <- colnames(input_cor_matrix)[i]
      int_vec <- input_cor_matrix[,i]
      int_vec <- sort(int_vec, decreasing = T)
      int_vec <- int_vec[2:length(int_vec)]
      int_vec_ind <- which(int_vec > correlation_cutoff)
      int_vec_trim <- int_vec[int_vec_ind]

      if(length(int_vec_trim) < min_edges){
        cells_below_cutoff <- cells_below_cutoff + 1
        int_cells <- int_vec[1:min_edges]
      } else {
        if(length(int_vec_trim) > max_edges){
          int_cells <- int_vec_trim[1:max_edges]
          cells_above_cutoff <- cells_above_cutoff + 1
        } else {
          int_cells <- int_vec_trim
        }
      }
      c1 <- rep(int_cell, length(int_cells))
      c2 <- names(int_cells)
      c3 <- as.vector(int_cells)

      end <- start + length(int_cells)-1

      num_edges[start:end,1] <- c1
      num_edges[start:end,2] <- c2
      num_edges[start:end,3] <- c3

      start <- end+1

    }

    blank_edges <- which(is.na(num_edges[,3]))
    num_edges <- num_edges[-blank_edges,]
    print(paste0(cells_below_cutoff, " cells were assigned edges below their correlation cutoff"))
    print(paste0(cells_above_cutoff, " cells were trimmed edges below their max edge cutoff"))


  }

  dim(num_edges)

  num_graph <- graph_from_data_frame(num_edges)

  if(weight == "standardized"){
    vals <- as.vector(scale(num_edges[,3]))
    vals <- (abs(min(vals))+vals)+1
    E(num_graph)$weight <- vals
  }

  if(weight == "raw"){
    E(num_graph)$weight <- num_edges[,3]
  }

  return(num_graph)

}

