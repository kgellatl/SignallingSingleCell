#' Calculate Library Size
#'
#' This function will count the number of reads / UMIs per cell
#'
#' @param input the input ex_sc
#' @export
#' @details
#' This will calculate the total number of UMIs on a per cell basis.
#' @examples
#' ex_sc_example <- calc_libsize(input = ex_sc_example)

calc_libsize <- function(input, suffix = NA){
  libsizes <- apply(exprs(input),2,sum)
  if(!is.na(suffix)){
    column <- paste0("UMI_sum", "_", suffix)
  } else {
    column <- "UMI_sum"
  }
  pData(input)[,column] <- libsizes
  return(input)
}
