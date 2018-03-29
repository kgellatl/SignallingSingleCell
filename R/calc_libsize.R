#' Calculate Library Size
#'
#' This function will count the number of reads / UMIs per cell
#'
#' @param input the input ex_sc
#' @export
#' @details
#' This will calculate total UMIs on a per cell basis.
#' @examples
#' ex_sc_example <- calc_libsize(input = ex_sc_example)

calc_libsize <- function(input){
  libsizes <- apply(exprs(input),2,sum)
  pData(input)$UMI_sum <- libsizes
  return(input)
}
