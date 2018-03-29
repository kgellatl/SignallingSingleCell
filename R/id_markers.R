#' ID markers
#'
#' This will perform marker gene identification
#'
#' @param input the input ex_sc
#' @param print_progress will print progress if TRUE
#' @export
#' @details
#' This will find marker genes for each cluster. First, for each cluster it calculates the fraction of cells expressing a given gene,
#' and the mean expression of that gene (within the cluster).
#' It then compares these values to the values of all other groups looking for genes which maximize the distance in both
#' fraction expressing and the mean expression.
#' @examples
#' ex_sc_example <- id_markers(input = ex_sc_example, print_progress = TRUE)

id_markers <- function(input, print_progress = TRUE){
  marker_input <- input
  fData(marker_input) <- fData(marker_input)[,-grep("marker_score", colnames(fData(marker_input)))]
  fData(marker_input)$tmp <- "tmp"
  if(print_progress == TRUE){
    print("Finding markers based on fraction expressing")
  }
  num_cluster <- length(unique(pData(marker_input)[,"Cluster"]))
  num_marker_genes <- nrow(fData(input))
  non0 <- c()
  percent <- c()
  gene_percent <- c()
  for(n in 1:length(unique(pData(marker_input)[,"Cluster"]))){
    cluster <- (paste0("Cluster", n))
    index <- which(pData(marker_input)[,"Cluster"] == cluster)
    cells <- rownames(pData(marker_input))[index]
    tmp_expr_matrix <- exprs(marker_input)[,cells]
    non0 <- which(tmp_expr_matrix > 0)
    tmp_expr_matrix[non0]<- 1
    for(i in 1:nrow(tmp_expr_matrix)){
      numExpr <- sum(tmp_expr_matrix[i,])
      percent <- numExpr/ncol(tmp_expr_matrix)
      gene_percent <- c(gene_percent, percent)
    }
  }
  gene_markers_fractionExpression <- matrix(data = gene_percent, nrow = nrow(fData(marker_input)), ncol = length(unique(pData(marker_input)[,"Cluster"])))
  colnames(gene_markers_fractionExpression) <- c(seq(1:num_cluster))
  rownames(gene_markers_fractionExpression) <- rownames(fData(marker_input))
  gene_markers_meandiff <- matrix(ncol = ncol(gene_markers_fractionExpression), nrow=nrow(gene_markers_fractionExpression))
  colnames(gene_markers_meandiff) <- c(seq(1:num_cluster))
  rownames(gene_markers_meandiff) <- rownames(fData(marker_input))
  for(i in 1:nrow(gene_markers_fractionExpression)){
    genetest <- gene_markers_fractionExpression[i,]
    for(j in 1:ncol(gene_markers_fractionExpression)){
      refer <- genetest[j]
      test <- genetest[-j]
      meandiff <- mean(refer - test)
      gene_markers_meandiff[i,j] <- meandiff
    }
  }
  cluster_markers <- matrix(ncol = ncol(gene_markers_meandiff), nrow=num_marker_genes)
  colnames(cluster_markers) <- c(seq(1:num_cluster))
  for(i in 1:ncol(gene_markers_meandiff)){
    sorted <- sort(gene_markers_meandiff[,i])
    interested_genes <- names(sorted)
    cluster_markers[,i] <- interested_genes
  }
  if(print_progress == TRUE){
    print("Finding markers based on mean expressing")
  }
  mean_exp <- c()
  mean_data <- c()
  for(n in 1:length(unique(pData(marker_input)[,"Cluster"]))){
    cluster <- (paste0("Cluster", n))
    index <- which(pData(marker_input)[,"Cluster"] == cluster)
    cells <- rownames(pData(marker_input))[index]
    tmp_expr_matrix <- exprs(marker_input)[,cells]
    for(i in 1:nrow(tmp_expr_matrix)){
      mean_exp <- mean(tmp_expr_matrix[i,])
      mean_data <- c(mean_data, mean_exp)
    }
  }
  gene_means <- matrix(data = mean_data, nrow = nrow(fData(marker_input)), ncol = length(unique(pData(marker_input)[,"Cluster"])))
  colnames(gene_means) <- c(seq(1:num_cluster))
  rownames(gene_means) <- rownames(fData(marker_input))
  gene_markers_foldchange <- matrix(ncol = ncol(gene_means), nrow=nrow(gene_means))
  colnames(gene_markers_foldchange) <- c(seq(1:num_cluster))
  rownames(gene_markers_foldchange) <- rownames(fData(marker_input))
  for(i in 1:nrow(gene_means)){
    genetest <- gene_means[i,]+.001
    for(j in 1:ncol(gene_means)){
      refer <- genetest[j]
      test <- genetest[-j]
      meandiff <- mean(log2(refer / test))
      gene_markers_foldchange[i,j] <- meandiff
    }
  }
  cluster_markers_FC <- matrix(ncol = ncol(gene_markers_foldchange), nrow=num_marker_genes)
  colnames(cluster_markers_FC) <- c(seq(1:num_cluster))
  for(i in 1:ncol(gene_markers_foldchange)){
    sorted <- sort(gene_markers_foldchange[,i])
    interested_genes <- names(sorted)
    cluster_markers_FC[,i] <- interested_genes
  }
  if(print_progress == TRUE){
    print("Merging Lists")
  }
  marker_genes <- list()
  for(i in 1:ncol(cluster_markers_FC)){
    cm <- cluster_markers[,i]
    cf <- cluster_markers_FC[,i]
    cm <- as.data.frame(cm)
    cf <- as.data.frame(cf)
    rownames(cm) <- cm$cm
    rownames(cf) <- cf$cf
    cm$rankcm <- seq(1:nrow(cm))
    cf$rankcf <- seq(1:nrow(cf))
    cm <- cm[order(cm$cm),]
    cf <- cf[order(cf$cf),]
    final <- cbind(cm, cf)
    final <- final[,c(2,4)]
    final$score <- apply(final,1,sum)
    final <- final[order(-final$score),]
    cluster_marker_rank <- rownames(final)
    val <- match(rownames(fData(input)), cluster_marker_rank)
    fData(marker_input)$name <- val
    colnames(fData(marker_input))[grep("name", colnames(fData(marker_input)))] <- paste0("Cluster",i,"_marker_score")
  }
  fData(marker_input) <- fData(marker_input)[,-grep("tmp", colnames(fData(marker_input)))]
  return(marker_input)
}



