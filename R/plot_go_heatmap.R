#' Construct Expression Set Class
#'
#' This function will take an input expression matrix and make an Expression Set Class
#'
#' @param input the input Go result from clusterProfiler
#' @param sig_val p.adjust of pvalue
#' @param cutoff the value to be included
#' @param max_categories mox categories included
#' @param subset_group table(input@compareClusterResult$group)
#' @export
#' @details
#' This will take an input matrix and create the expression set class that further analysis
#' will be written to
#' @examples
#' construct_ex_sc(input = sc_dat)

plot_go_heatmap <- function(input, sig_val = "p.adjust", cutoff = 0.1, max_categories = 30, subset_group = F){
  if(class(input)[1] == "enrichResult"){

    genes <- input@result$geneID
    descriptions <- input@result$Description
    sig_vals <- as.vector(unlist(input@result[sig_val]))
    sig <- which(sig_vals < cutoff)
    if(length(sig) > max_categories){
      sig <- sig[1:max_categories]
    }

    top_genes <- unique(unlist(strsplit(genes[sig], split = "/")))
    go_matrix <- matrix(0, nrow = length(top_genes), ncol = length(sig))
    rownames(go_matrix) <- top_genes
    colnames(go_matrix) <- descriptions[sig]
    go_matrix <- as.data.frame(go_matrix)

    for (i in 1:ncol(go_matrix)) {
      desc <- descriptions[i]
      gset <- genes[i]
      gset <- unlist(strsplit(gset, split = "/"))
      ind <- match(gset, rownames(go_matrix))
      go_matrix[ind, i] <- 1
    }
    sums <- apply(go_matrix,1,sum)
    for (i in 1:nrow(go_matrix)) {
      int <- go_matrix[i,]
      ind <- which(int > 0)
      int[ind] <- sums[i]
      go_matrix[i,] <- int
    }

    go_sc <- construct_ex_sc(go_matrix)
    numcol <- dim(table(unlist(go_matrix)))
    colpal <- viridis::viridis(numcol)
    colpal[1] <- "gray"
    plot_heatmap(go_sc, genes = rownames(go_sc), type = "single_cell", cluster_by = "both", scale_by = F, text_angle = 45, color_pal = colpal)
  }


  if(class(input)[1] == "compareClusterResult"){

    genes <- input@compareClusterResult$geneID
    descriptions <- input@compareClusterResult$Description
    sig_vals <- as.vector(unlist(input@compareClusterResult[sig_val]))
    groups <- input@compareClusterResult$group

    if(subset_group != FALSE){
      ind <- which(groups == subset_group)
      genes <- genes[ind]
      descriptions <- descriptions[ind]
      sig_vals <- sig_vals[ind]

    } else {
      stop("For compare cluster results please provide a subset_group")
    }

    sig <- which(sig_vals < cutoff)
    if(length(sig) > max_categories){
      sig <- sig[1:max_categories]
    }

    top_genes <- unique(unlist(strsplit(genes[sig], split = "/")))
    go_matrix <- matrix(0, nrow = length(top_genes), ncol = length(sig))
    rownames(go_matrix) <- top_genes
    colnames(go_matrix) <- descriptions[sig]
    go_matrix <- as.data.frame(go_matrix)

    for (i in 1:ncol(go_matrix)) {
      desc <- descriptions[i]
      gset <- genes[i]
      gset <- unlist(strsplit(gset, split = "/"))
      ind <- match(gset, rownames(go_matrix))
      go_matrix[ind, i] <- 1
    }
    sums <- apply(go_matrix,1,sum)
    for (i in 1:nrow(go_matrix)) {
      int <- go_matrix[i,]
      ind <- which(int > 0)
      int[ind] <- sums[i]
      go_matrix[i,] <- int
    }

    go_sc <- construct_ex_sc(go_matrix)
    numcol <- dim(table(unlist(go_matrix)))
    colpal <- viridis::viridis(numcol)
    colpal[1] <- "gray"
    plot_heatmap(go_sc, genes = rownames(go_sc), type = "single_cell", cluster_by = "both", scale_by = F, text_angle = 45, color_pal = colpal)


  }


}



