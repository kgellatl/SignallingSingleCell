#' Return markers
#'
#' This will return the marker genes.
#'
#' @param input the input ex_sc
#' @param num_markers the number of markers to return for each cell type
#' @export
#' @details
#' This will return the marker and reference genes stored in fData(input) in a convienient list format.
#' @examples
#' marker_list <- return_markers(input = ex_sc, num_markers = 50)

return_markers <- function(input, num_markers){
  dat <- fData(input)
  Clusters <- sort(unique(pData(input)$Cluster))
  markers <- vector(mode = "list", length = length(Clusters))
  for(i in 1:length(Clusters)){
    cint <- Clusters[i]
    ind <- grep(paste0(cint, "_marker"), colnames(dat))
    ind2 <- match(seq(1:num_markers), dat[,ind])
    marker_final <- rownames(dat)[ind2]
    markers[[i]] <- marker_final
  }
  names(markers) <- c(paste0("Cluster", seq(1:length(Clusters))))
  return(markers)
}


