#' Calculate Gene Stats
#'
#' This function will calculate basic statistics about gene usage
#'
#' @param input the input ex_sc
#' @export
#' @details
#' This will calculate the total number of UMIs on a per cell basis.
#' @examples
#' ex_sc_example <- calc_genestats(input = ex_sc_example)

calc_genestats <- function(input){

  pData(input)$genes_expressed <- apply(exprs(input), 2, FUN = function(x) length(which(x > 0)))
  fData(input)$cells_expressing <- apply(exprs(input), 1, FUN = function(x) length(which(x > 0)))

  return(input)
}
