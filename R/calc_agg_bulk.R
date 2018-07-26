#' Calculate UPM expression values across pData values
#'
#' This will calculate UMIs per million UPM expression values across pData columns. Useful for heatmaps and networking analysis.
#'
#' @param input the input ex_sc.
#' @param aggregate_by The pData variables to break by
#' @param group_by The pData variables that contains distinct groups if included in aggregate_by (eg Timepoint, Condition, etc)
#' This is used to calculate portions internally to this group
#' @param cutoff a  value below which gene expression values will be set to 0 for the mean.
#' This is useful for removing nodes from networks that contain only a couple of cells.
#' @export
#' @details
#' Utilize information stored in pData to control the plot display.
#' @examples
#' plot_tsne_metadata(ex_sc_example, color_by = "UMI_sum", title = "UMI_sum across clusters", facet_by = "Cluster", ncol = 3)

calc_agg_bulk <- function(input, aggregate_by, group_by = FALSE, cutoff = FALSE){
  if(group_by != FALSE){
    ind <- match(group_by, aggregate_by)
    if(ind != 1){
      stop("Please provide the group_by value first in the aggreggate_by argument")
    }
  }
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
  groups <- unique(pData(input)[,group_by])
  upm_vals <- c()
  num_cells_vals <- c()
  num_genes_vals <- c()
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
    num_cells <- length(full_match)
    num_cells_vals <- c(num_cells_vals, num_cells)
    if(length(full_match) > 1){
      tmp <- exprs(input)[,full_match]
      upm <- (apply(tmp, 1, sum)/sum(tmp))*1000000
      if(cutoff != FALSE){
        tmp2 <- tmp
        tmp2[which(tmp2 > 0)] <- 1
        gSums <- apply(tmp2,1,sum)
        frac <- gSums / num_cells
        zero_out <- which(frac < cutoff)
        upm[zero_out] <- 0
      }
    }
    expressed <- length(which(upm > 0))#not right!!!! for the fraction calc is different..
    upm_vals <- c(upm_vals, upm)
    num_genes_vals <- c(num_genes_vals, expressed)
  }
  bulks$numcells <- num_cells_vals
  bulks$numgenes <- num_genes_vals
  bulks$proportion <- 0
  if(group_by != FALSE){
    for (i in 1:length(groups)) {
      ind <- grep(groups[i], bulks[,group_by])
      total <- sum(bulks$numcells[ind])
      bulks$proportion <- round((bulks$numcells/total)*100,2)
    }
  } else {
    total <- sum(bulks$numcells)
    bulks$proportion <- round((bulks$numcells/total)*100,2)
  }
  bulk <- matrix(upm_vals, nrow = nrow(exprs(input)))
  rownames(bulk) <- rownames(exprs(input))
  colnames(bulk) <- seq(1:ncol(bulk))
  for (l in 1:nrow(bulks)) {
    cname <- bulks[l,-c(match(c("numcells", "proportion", "numgenes"), colnames(bulks)))]
    cname2 <- c()
    for (i in 1:length(cname)) {
      cint <- as.character(cname[[i]])
      cname2 <- c(cname2, cint)
    }
    cname <- cname2
    cnum <- bulks[l,"numcells"]
    cpro <- bulks[l,"proportion"]
    cgen <- bulks[l,"numgenes"]
    cname <- paste0(c(cname, "num_genes", cgen, "num_cells", cnum, "percent", cpro, "bulk"), collapse = "_")
    colnames(bulk)[l] <- cname
  }
  fData(input) <- cbind(fData(input), bulk)
  return(input)
}
