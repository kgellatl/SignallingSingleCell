#' Construct Expression Set Class
#'
#' This function will take an input expression matrix and make an Expression Set Class
#'
#' @param input the input Go result from clusterProfiler
#' @param sig_val p.adjust of pvalue
#' @param cutoff the value to be included
#' @param max_categories mox categories included
#' @param color_by what the cells should be shaded by. "count" or "sig"
#' @param subset_group table(input@compareClusterResult$group)
#' @param prune If true will prune the resulting heatmap based on jaccard similarity
#' @param cutoff_jaccard The jaccard similarity index of terms to create a new group
#' @param word_similarity for groups of go terms, the fraction of the original terms that a given word must occur in
#' @export
#' @details
#' This will take the output of a GO enrichment analysis or GSEA object from clusterProfiler
#' and merge terms based on the overlap between the enriched genes
#' @examples
#' plot_go_heatmap(input = gsea_result)

plot_go_heatmap <- function(input,
                            sig_val = "p.adjust",
                            cutoff = 0.1,
                            max_categories = 30,
                            color_by = "sig",
                            subset_group = F,
                            prune = F,
                            cutoff_jaccard = 1,
                            word_similarity = 0.5) {

  if (class(input)[1] == "enrichResult") {
    genes <- input@result$geneID
    descriptions <- input@result$Description
    sig_vals <- as.vector(unlist(input@result[sig_val]))
    sig <- which(sig_vals < cutoff)
  }
  ###### GSEA Result
  if(class(input)[1] == "gseaResult"){
    genes <- input@result$core_enrichment
    descriptions <- input@result$Description
    sig_vals <- as.vector(unlist(input@result[sig_val]))
    sig <- which(sig_vals < cutoff)
  }
  ###### merged table input
  if(class(input)[1] == "data.frame"){
    genes <- input$core_enrichment
    descriptions <- input$Description
    sig_vals <- as.vector(unlist(input[sig_val]))
    sig <- which(sig_vals < cutoff)
  }
  if (class(input)[1] == "compareClusterResult") {
    genes <- input@compareClusterResult$geneID
    descriptions <- input@compareClusterResult$Description
    sig_vals <- as.vector(unlist(input@compareClusterResult[sig_val]))
    groups <- input@compareClusterResult$group
    if (subset_group != FALSE) {
      ind <- which(groups == subset_group)
      genes <- genes[ind]
      descriptions <- descriptions[ind]
      sig_vals <- sig_vals[ind]
    }
    else {
      stop("For compare cluster results please provide a subset_group")
    }
    sig <- which(sig_vals < cutoff)
  }
  if (length(sig) > max_categories) {
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
  if (prune) {
    jaccard_grid <- expand.grid(colnames(go_matrix), colnames(go_matrix),
                                stringsAsFactors = F)
    jaccard_grid <- as.data.frame(jaccard_grid)
    jaccard_grid$similarity <- 0
    jaccard <- function(M) {
      sums = rowSums(M)
      similarity = length(sums[sums == 2])
      total = length(sums[sums == 1]) + similarity
      return(similarity/total)
    }
    for (i in 1:nrow(jaccard_grid)) {
      int_terms <- jaccard_grid[i, ]
      M <- go_matrix[, c(int_terms$Var1, int_terms$Var2)]
      jaccard_val <- jaccard(M)
      jaccard_grid$similarity[i] <- jaccard_val
    }
    jaccard_matrix <- matrix(jaccard_grid$similarity, ncol = ncol(go_matrix),
                             nrow = ncol(go_matrix))
    colnames(jaccard_matrix) <- colnames(go_matrix)
    rownames(jaccard_matrix) <- colnames(go_matrix)
    similar_results <- apply(jaccard_matrix, 1, function(x) which(x >=
                                                                    cutoff_jaccard))
    unique_terms <- c()
    for (i in 1:length(similar_results)) {
      int <- similar_results[[i]]
      int <- as.numeric(int)
      int <- paste0(int, collapse = "_")
      unique_terms <- unique(c(unique_terms, int))
    }
    new_terms <- c()
    for (i in 1:length(unique_terms)) {
      int <- unique_terms[i]
      all_full_terms <- colnames(go_matrix)[as.numeric(unlist(strsplit(int,
                                                                       split = "_")))]
      all_full_terms <- gsub(" ", "-", all_full_terms)
      all_full_terms <- gsub("-", "_", all_full_terms)
      if (length(all_full_terms) > 1) {
        all_words <- unlist(strsplit(all_full_terms,
                                     "_"))
        counts_words <- table(all_words)
        common <- names(which(counts_words >= word_similarity *
                                length(all_full_terms)))
        if (length(common) == 0) {
          common <- names(which(counts_words == max(counts_words)))
        }
        new_term <- unique(all_words[all_words %in%
                                       common])
        new_term <- paste0(new_term, collapse = "_")
        new_terms <- c(new_terms, new_term)
      }
      else {
        new_terms <- c(new_terms, all_full_terms)
      }
    }
    new_go_matrix <- matrix(ncol = length(new_terms), nrow = nrow(go_matrix))
    colnames(new_go_matrix) <- new_terms
    rownames(new_go_matrix) <- rownames(go_matrix)
    for (i in 1:length(unique_terms)) {
      ind <- as.numeric(unlist(strsplit(unique_terms[i],
                                        "_")))
      if (length(ind) == 1) {
        new_go_matrix[, i] <- go_matrix[, ind]
      }
      else {
        sums <- apply(go_matrix[, ind], 1, sum)
        new_go_matrix[, i] <- sums
      }
    }
    colnames(new_go_matrix) <- ave(as.character(colnames(new_go_matrix)),
                                   colnames(new_go_matrix), FUN = function(x) if (length(x) >
                                                                                  1)
                                     paste0(x[1], "(", seq_along(x), ")")
                                   else x[1])
    go_sc <- construct_ex_sc(new_go_matrix)
    numcol <- dim(table(unlist(new_go_matrix)))
    colpal <- viridis::viridis(numcol)
    colpal[1] <- "gray"
    res <- vector(mode = "list", length = 2)
    g <- plot_heatmap(go_sc, genes = rownames(go_sc), type = "single_cell",
                      cluster_by = "both", scale_by = F, text_angle = 60,
                      color_pal = colpal)
    res[[1]] <- g
    res[[2]] <- colnames(new_go_matrix)
    return(res)
  }
  else {
    if (color_by == "count") {
      sums <- apply(go_matrix, 1, sum)
      for (i in 1:nrow(go_matrix)) {
        int <- go_matrix[i, ]
        ind <- which(int > 0)
        int[ind] <- sums[i]
        go_matrix[i, ] <- int
      }
    }
    if (color_by == "sig") {
      for (i in 1:ncol(go_matrix)) {
        int <- go_matrix[, i]
        int[which(int > 0)] <- sig_vals[i]
        int[which(int > 0)] <- log10(int[which(int >
                                                 0)])
        int[which(int == 0)] <- log10(int[which(int ==
                                                  0)] + 10)
        go_matrix[, i] <- int
      }
    }
    all_full_terms <- gsub(" ", "-", colnames(go_matrix))
    all_full_terms <- gsub("-", "_", all_full_terms)
    colnames(go_matrix) <- all_full_terms
    go_sc <- construct_ex_sc(go_matrix)
    numcol <- dim(table(unlist(go_matrix)))
    colpal <- viridis::viridis(numcol)
    if (color_by == "count") {
      colpal[1] <- "gray"
    }
    if (color_by == "sig") {
      colpal <- rev(colpal)
      colpal[length(colpal)] <- "gray"
    }
    g <- plot_heatmap(go_sc, genes = rownames(go_sc), type = "single_cell",
                      cluster_by = "both", scale_by = F, color_pal = colpal)
    return(g)
  }
}
