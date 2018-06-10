#' identify receptors and ligands
#'
#' This function will id receptors and ligands in your data and mark them in fData
#'
#' @param input the input ex_sc
#' @export
#' @details
#' This will use data(Receptor_Ligand_Data) to write to fData 2 columns, one for
#' receptors and one for ligands. TRUE FALSE statements will be written depending on their presence.
#' Note that this uses HUMAN GENE NAMES!!
#' @examples
#' ex_sc_example <- id_rl(input = ex_sc_example)

id_rl <- function(input){
  data("Receptor_Ligand_Data")
  ligs <- unique(Receptor_Ligand_Data$Ligand.ApprovedSymbol)
  recs <- unique(Receptor_Ligand_Data$Receptor.ApprovedSymbol)
}
