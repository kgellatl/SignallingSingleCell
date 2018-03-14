#' Filter Cells based on Expression
#'
#' This function will filter based on min number of UMIs, max number of UMIs, and genes based on expression across all cells.
#'
#' @param input the input data matrix
#' @param minUMI min number of UMIs per cell
#' @param maxUMI max number of UMIs per cell
#' @param threshold UMI threshold for gene detection
#' @param minCells number of cells expressed above threshold for a given gene
#' @export
#' @details
#' When processing the data, low quality cells may contain very few UMIs, while some overrepresented cell barcodes may indicate barcode bleedover or celll doublets. Filtering out both low and high UMI count cells is recommended before normalization.
#' @examples
#' filtered_data <- filter_UMIs(input = mDC_0hr_1hr_4hr_CLEAN, minUMI = 1000, maxUMI = 10000, threshold = 1, minCells = 10)

filter_UMIs <- function(input, minUMI, maxUMI=NA, threshold=NA, minCells=NA){
  cSum <- apply(input,2,sum)
  Index <- which(cSum < minUMI)
  input <- input[,-Index]
  if(is.na(maxUMI) != TRUE){
    cSum <- apply(input,2,sum)
    Index <- which(cSum > maxUMI)
    input <- input[,-Index]
  }
  cSum <- apply(input,2,sum)
  if(is.na(threshold) != TRUE){
    if(is.na(threshold) != TRUE){
    }
    compThresh <- threshold
    gCount <- apply(input,1,function(x) length(which(x>compThresh)))
    input <- input[(which(gCount > minCells)),]
  }
  return(input)
}
