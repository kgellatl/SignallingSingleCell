#' Calculate Library Size
#'
#' This function will count the number of reads / UMIs per cell
#'
#' @param input the input data
#' @export
#' @details
#' This will calculate library sizes
#' @examples
#' ex_sc_example <- calc_libsize(ex_sc_example)

calc_libsize <- function(input){
  libsizes <- apply(exprs(input),2,sum)
  pData(input)$UMI_sum <- libsizes
  return(input)
}
