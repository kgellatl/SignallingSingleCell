#' This function will filter based on min number of UMIs, max number of UMIs, and genes based on expression across all cells.
#' @param minUMI min number of UMIs per cell
#' @param maxUMI max number of UMIs per cell
#' @param threshold UMI threshold for gene detection
#' @param minCells number of cells expressed above threshold for a given gene
#' @keywords cats
#' @export
#' @examples
#' filter_UMIs()
filter_UMIs <- function(input, minUMI, maxUMI, threshold, minCells){
  cSum <- apply(input,2,sum)
  range(cSum)
  hist(cSum)
  Index <- which(cSum < minUMI)
  input_minUMI <- input[,-Index]
  cSum <- apply(input_minUMI,2,sum)
  Index <- which(cSum > maxUMI)
  input_minUMI_maxUMI <- input_minUMI[,-Index]
  cSum <- apply(input_minUMI_maxUMI,2,sum)
  hist(cSum)
  compThresh <- threshold
  gCount <- apply(input_minUMI_maxUMI,1,function(x) length(which(x>compThresh)))
  input_filtered <- input_minUMI_maxUMI[(which(gCount > minCells)),]
  return(input_filtered)
}
