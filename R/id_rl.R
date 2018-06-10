#' identify receptors and ligands
#'
#' This function will id receptors and ligands in your data and mark them in fData
#'
#' @param input the input ex_sc
#' @export
#' @details
#' This will use data(Receptor_Ligand_Data) to write to fData 3 columns, one for
#' receptors, one for ligands, and one for the expressed pairs. TRUE FALSE in the first 2 if the gene is a ligand or receptor respectively.
#' For all pairs where both a receptor and ligand are expressed, a 3rd column named
#' networks_expressed_pairs will be created with all ligands expressed that bind that receptor.
#' Note that this uses HUMAN GENE NAMES!!
#' @examples
#' ex_sc_example <- id_rl(input = ex_sc_example)

id_rl <- function(input){
  data("Receptor_Ligand_Data")
  ligs <- unique(Receptor_Ligand_Data$Ligand.ApprovedSymbol)
  recs <- unique(Receptor_Ligand_Data$Receptor.ApprovedSymbol)
  expressed_genes <- rownames(fData(input))
  fData(input)$networks_Receptors <- FALSE
  fData(input)$networks_ligands <- FALSE
  fData(input)$networks_expressed_pairs <- ""
  fData(input)$networks_ligands[which(is.na(match(expressed_genes, ligs)) == FALSE)] <- TRUE
  fData(input)$networks_Receptors[which(is.na(match(expressed_genes, recs)) == FALSE)] <- TRUE
  for (i in 1:nrow(Receptor_Ligand_Data)) {
    int <- Receptor_Ligand_Data[i,]
    pair <- unlist(strsplit(int$Pair.Name, split = "_"))
    ind_lig <- match(pair[1], rownames(fData(input)))
    ind_rec <- match(pair[2], rownames(fData(input)))
    lig_state <- fData(input)$networks_ligands[ind_lig]
    rec_state <- fData(input)$networks_Receptors[ind_rec]
    if(is.na(lig_state) == TRUE | is.na(rec_state) == TRUE){
    } else
    if(lig_state == TRUE & rec_state == TRUE){
      fData(input)$networks_expressed_pairs[ind_rec] <- paste0(fData(input)$networks_expressed_pairs[ind_rec], pair[1], "_")
    }
  }
  fData(input)$networks_expressed_pairs <- substr(fData(input)$networks_expressed_pairs,1,nchar(fData(input)$networks_expressed_pairs)-1)
  fData(input)$networks_expressed_pairs[which(fData(input)$networks_expressed_pairs == "")] <- NA
  return(input)
}
