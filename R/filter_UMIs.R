#' Filter Cells based on Expression
#'
#' This function will filter based on min number of UMIs, max number of UMIs, and genes based on expression across all cells.
#'
#' @param minUMI min number of UMIs per cell
#' @param maxUMI max number of UMIs per cell
#' @param threshold UMI threshold for gene detection
#' @param minCells number of cells expressed above threshold for a given gene
#' @export
#' @details
#' When processing the data, low quality cells may contain very few UMIs, while some overrepresented cell barcodes may indicate barcode bleedover or celll doublets. Filtering out both low and high UMI count cells is recommended before normalization.
#' @examples
#' filtered_data <- filter_UMIs(input = mDC_0hr_1hr_4hr_CLEAN, minUMI = 1000, maxUMI = 10000, threshold = 1, minCells = 10)

filter_UMIs <- function(input, minUMI, maxUMI, threshold, minCells){
  cSum <- apply(input,2,sum)
  range(cSum)
  hist(cSum, main = "UMI distribution pre Filter")
  Index <- which(cSum < minUMI)
  input_minUMI <- input[,-Index]
  cSum <- apply(input_minUMI,2,sum)
  Index <- which(cSum > maxUMI)
  input_minUMI_maxUMI <- input_minUMI[,-Index]
  cSum <- apply(input_minUMI_maxUMI,2,sum)
  hist(cSum, main = "UMI distribution post Filter")
  compThresh <- threshold
  gCount <- apply(input_minUMI_maxUMI,1,function(x) length(which(x>compThresh)))
  input_filtered <- input_minUMI_maxUMI[(which(gCount > minCells)),]
  return(input_filtered)
}
