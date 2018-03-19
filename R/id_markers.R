#' ID markers
#' This will perform marker gene identification
#'
#' @param input the input data matrix.
#' @param num_markers the number of markers to return
#' @export
#' @details
#' This will find marker genes for each cluster.
#' @examples
#' ex_sc_example <- id_markers(input = exprs(ex_sc_example), num_markers = 50, num_reference = 1000)

id_markers <- function(input, num_markers, num_reference){
  print("Finding markers based on fraction expressing")
  marker_input <- input
  num_cluster <- length(unique(pData(marker_input)[,"Cluster"]))
  num_marker_genes <- nrow(fData(input))
  non0 <- c()
  percent <- c()
  gene_perent <- c()
  for(n in 1:length(unique(pData(marker_input)[,"Cluster"]))){
    if (typeof(pData(marker_input)[,"Cluster"]) == "character"){
      cluster <- (paste0("Cluster", n))
    } else {
      cluster <- (n)
    }
    index <- which(pData(marker_input)[,"Cluster"] == cluster)
    cells <- rownames(pData(marker_input))[index]
    tmp_expr_matrix <- exprs(marker_input)[,cells]
    non0 <- which(tmp_expr_matrix > 0)
    tmp_expr_matrix[non0]<- 1
    for(i in 1:nrow(tmp_expr_matrix)){
      numExpr <- sum(tmp_expr_matrix[i,])
      percent <- numExpr/ncol(tmp_expr_matrix)
      gene_perent <- c(gene_perent, percent)
    }
  }

  gene_markers_fractionExpression <- matrix(data = gene_perent, nrow = nrow(fData(marker_input)), ncol = length(unique(pData(marker_input)[,"Cluster"])))
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
    interested_genes <- names((tail(sorted, n=num_marker_genes)))
    cluster_markers[,i] <- interested_genes
  }

  # write.csv(cluster_markers, file = "cluster_markers.csv")

  ###########################################################################################################
  ### Marker Identification Mean Expression
  ###########################################################################################################
  print("Finding markers based on mean expressing")

  mean_exp <- c()
  mean_data <- c()
  for(n in 1:length(unique(pData(marker_input)[,"Cluster"]))){
    if (typeof(pData(marker_input)[,"Cluster"]) == "character"){
      cluster <- (paste0("Cluster", n))
    } else {
      cluster <- (n)
    }
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
    interested_genes <- names((tail(sorted, n=num_marker_genes)))
    cluster_markers_FC[,i] <- interested_genes
  }

  # marker_input_blood_ependymal_oligo_astro_vascular_neuron_foldchange <- cluster_markers_FC
  # write.csv(marker_input_blood_ependymal_oligo_astro_vascular_neuron_foldchange, file = "marker_input_blood_ependymal_oligo_astro_vascular_neuron_foldchange.csv")

  ###########################################################################################################
  ### Merging to find Markers
  ###########################################################################################################

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
    final <- final[order(final$score),]
    marker_genes[[i]] <- rev(rownames(tail(final, num_markers)))
  }

  ###########################################################################################################
  ### Merging to find Common Genes
  ###########################################################################################################
  print("Finding reference genes")

  num_marker_genes <- num_reference
  high_frac_exprs <- matrix(ncol = ncol(gene_markers_meandiff), nrow=num_marker_genes)
  colnames(high_frac_exprs) <- c(seq(1:num_cluster))
  for(i in 1:ncol(gene_markers_fractionExpression)){
    sorted <- sort(gene_markers_fractionExpression[,i])
    interested_genes <- names((tail(sorted, n=num_marker_genes)))
    high_frac_exprs[,i] <- interested_genes
  }

  cluster_no_FC <- matrix(ncol = ncol(gene_markers_foldchange), nrow=num_marker_genes)
  colnames(cluster_no_FC) <- c(seq(1:num_cluster))
  for(i in 1:ncol(gene_markers_foldchange)){
    distance <- sort(abs(gene_markers_foldchange[,i]))
    interested_genes <- names((head(distance, n=num_marker_genes)))
    cluster_no_FC[,i] <- interested_genes
  }

  common_genes <- list()
  for(i in 1:ncol(high_frac_exprs)){
    cluster_common <- intersect(high_frac_exprs[,i], cluster_no_FC[,i])
    common_genes[[i]] <- cluster_common
  }

  common_genes_recurring <- table(unlist(common_genes))
  common_genes_recurring <- rev(sort(common_genes_recurring))
  #  <- common_genes_recurring
  # save(common_genes_recurring, file = "common_genes_recurring.Rdata")

  high_frac_exprs <- matrix(ncol = ncol(gene_markers_meandiff), nrow=num_marker_genes)
  colnames(high_frac_exprs) <- c(seq(1:num_cluster))
  for(i in 1:ncol(gene_markers_fractionExpression)){
    sorted <- sort(gene_markers_fractionExpression[,i])
    interested_genes <- names((tail(sorted, n=num_marker_genes)))
    high_frac_exprs[,i] <- interested_genes
  }


  common_genes_recurring_allexprs <- table(unlist(cluster_no_FC))
  common_genes_recurring_allexprs <- rev(sort(common_genes_recurring_allexprs))
  common_genes_recurring_allexprs <- intersect(names(common_genes_recurring_allexprs), unlist(high_frac_exprs))
  marker_genes

  fData(marker_input)$Markers <- NA
  for(i in 1:length(marker_genes)){
    set <- marker_genes[[i]]
    ind <- match(set, rownames(fData(marker_input)))
    fData(marker_input)$Markers[ind] <- paste0("Cluster", i)
  }
  fData(marker_input)$Reference <- NA
  ind <- match(common_genes_recurring_allexprs, rownames(fData(marker_input)))
  fData(marker_input)$Reference[ind] <- "Reference"
  return(marker_input)
}



