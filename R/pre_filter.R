#' Filter Cells based on Expression
#'
#' This function will filter based on min number of UMIs, max number of UMIs, and genes based on expression across all cells.
#'
#' @param input the input data
#' @param minUMI min number of UMIs per cell
#' @param maxUMI max number of UMIs per cell
#' @param threshold UMI threshold for gene detection
#' @param minCells number of cells expressed above threshold for a given gene
#' @export
#' @details
#' When processing the data, low quality cells may contain very few UMIs, while some overrepresented cell barcodes may indicate barcode bleedover or celll doublets. Filtering out both low and high UMI count cells is recommended before normalization.
#' @examples
#' filtered_data <- filter_UMIs(input = ex_sc_example, minUMI = 1000, maxUMI = 10000, threshold = 1, minCells = 10)

pre_filter <- function(input, minUMI, maxUMI=NA, threshold=NA, minCells=NA){
  if(is.na(threshold) != TRUE){
    if(is.na(minCells) != TRUE){
    }
    print("Filtering Genes")
    gCount <- apply(exprs(input),1,function(x) length(which(x>=threshold)))
    input <- input[(which(gCount >= minCells)),]
  }
  print("Filtering Low Cells")
  cSum <- apply(exprs(input),2,sum)
  if(minUMI < min(cSum)){
    stop("minUMI is less than minimum number of UMIs in data")
  }
  low_c <- which(cSum < minUMI)
  if(is.na(maxUMI) != TRUE){
    if(maxUMI > max(cSum)){
      stop("maxUMI is greater than maximum number of UMIs in data")
    }
    print("Filtering High Cells")
    high_c <- which(cSum > maxUMI)
    input <- input[,-c(low_c, high_c)]
  } else (
    input <- input[,-low_c]
  )
  return(input)
}
