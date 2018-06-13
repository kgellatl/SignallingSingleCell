#' Identifies all R / L interactions
#'
#' This function will map all RL interactions
#'
#' @param input the input ex_sc
#' @param break_by the pData columns calc_agg_bulk was calculated on
#' @export
#' @details
#' This will use the calc_agg_bulk results to ID networks
#' @examples
#' ex_sc_example <- id_rl(input = ex_sc_example)

plot_rl_network <- function(input, break_by = FALSE){
  all_dat <- which(!is.na(fData(input)[,"networks_expressed_pairs"]))
  all_dat <- fData(input)[all_dat,]
  all_expr <- grep("_bulk", colnames(all_dat))
  all_expr <- all_dat[,all_expr]
  result <- vector(mode = "list", length = ncol(all_expr))

  list(colnames(all_expr), list(colnames(all_expr)))


}
