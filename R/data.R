#' Receptor Ligand Data
#'
#' @description  Simplified Receptor Ligand database
#' @usage data(Receptor_Ligand_Data)
#'
#' Self note... Create .rda file in /data, then a .R file called data
#' run devtools::load_all()
#' run roxygen2::roxygenise()
#'
#' @references Ramilowski, Jordan A., et al. "A draft network of ligand–receptor-mediated multicellular signalling in human." Nature communications 6 (2015): 7866.
#'
#' A dataset containing a simplified version of the datafile in reference above.
#'
#' @format A data frame with 2557 rows and 5 variables:
#' \itemize{
#'   \item Pair.Name: A specific receptor ligand pair
#'   \item Ligand.ApprovedSymbol: the gene name
#'   \item Ligand.Name: common name of the ligand
#'   \item Receptor.ApprovedSymbol: the gene name
#'   \item Receptor.Name: common name of the receptor
#' }
"Receptor_Ligand_Data"
