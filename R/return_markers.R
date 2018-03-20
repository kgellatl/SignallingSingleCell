#' Return markers
#'
#' This will return the marker genes.
#'
#' @param input the input data matrix.
#' @export
#' @details
#' This will return the marker and reference genes stored in fData(input) in a convienient list format.
#' @examples
#' marker_list <- return_markers(input)

return_markers <- function(input){
  Clusters <- sort(unique(pData(input)$Cluster))
  markers <- vector(mode = "list", length = length(Clusters)+1)
  for(i in 1:length(Clusters)){
    cint <- Clusters[i]
    ind <- grep(cint, fData(input)$Markers)
    marker <- rownames(fData(input))[ind]
    markers[[i]] <- marker
  }
  reference <- grep("Reference", fData(input)$Reference)
  markers[[length(Clusters)+1]] <- rownames(fData(input))[reference]
  names(markers) <- c(paste0("Cluster", seq(1:length(Clusters))), "Reference")
  return(markers)
}



