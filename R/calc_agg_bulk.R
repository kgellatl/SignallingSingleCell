#' Calculate mean expression values across pData values
#'
#' This will calculate mean expression values across pData columns. Useful for heatmaps and networking analysis.
#'
#' @param input the input ex_sc.
#' @param aggregate_by The pData variables to break by
#' @export
#' @details
#' Utilize information stored in pData to control the plot display.
#' @examples
#' plot_tsne_metadata(ex_sc_example, color_by = "UMI_sum", title = "UMI_sum across clusters", facet_by = "Cluster", ncol = 3)

calc_agg_bulk <- function(input, aggregate_by){
  check <- grep("bulk", colnames(fData(input)))
  if(length(check) > 0){
    fData(input) <- fData(input)[,-check]
  }
  to_expand <- vector("list", length(aggregate_by))
  for(i in 1:length(aggregate_by)) {
    var <- aggregate_by[i]
    vars <- unique(pData(input)[,var])
    to_expand[[i]] <- sort(vars)
  }
  names(to_expand) <- aggregate_by
  bulks <- expand.grid(to_expand, stringsAsFactors = FALSE)
  colnames(bulks) <- c(aggregate_by)
  mean_vals <- c()
  for (j in 1:nrow(bulks)) {
    int <- bulks[j,]
    full_match <- c()
    for (k in 1:length(int)) {
      ind <- which(pData(input)[,colnames(bulks)[k]] == int[[k]])
      if (k == 1){
        full_match <- c(full_match, ind)
      } else {
        full_match <- intersect(full_match, ind)
      }
    }
    mean <- apply(exprs(input)[,full_match],1,mean)
    mean_vals <- c(mean_vals, mean)
  }
  bulk <- matrix(mean_vals, nrow = nrow(exprs(input)))
  rownames(bulk) <- rownames(exprs(input))
  colnames(bulk) <- seq(1:ncol(bulk))
  for (l in 1:nrow(bulks)) {
    cname <- bulks[l,]
    cname <- paste0(c(cname, "bulk"), collapse = "_")
    colnames(bulk)[l] <- cname
  }
  fData(input) <- cbind(fData(input), bulk)
  return(input)
}

