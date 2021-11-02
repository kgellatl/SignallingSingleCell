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

id_rl <- function(input, database = "helft"){
  if(database == "helft"){
    Receptor_Ligand_Data <- SignallingSingleCell:::Receptor_Ligand_Data
  }
  else {
    Receptor_Ligand_Data <- database
  }
  ligs <- unique(Receptor_Ligand_Data$Ligand.ApprovedSymbol)
  recs <- unique(Receptor_Ligand_Data$Receptor.ApprovedSymbol)
  ind <- grep("bulk", colnames(fData(input)))
  if(length(ind) ==0){
    stop("Calculate aggregate bulk values first")
  }
  mat <- fData(input)[,ind]
  gsum <- apply(mat,1,sum)
  expressed_genes_bulk <- names(which(gsum >0))
  expressed_genes <- rownames(fData(input))
  fData(input)$networks_Receptors <- FALSE
  fData(input)$networks_ligands <- FALSE
  fData(input)$networks_ligs_to_receptor <- ""
  fData(input)$networks_receptor_to_ligs <- ""
  fData(input)$networks_ligands[which(is.na(match(expressed_genes, ligs)) == FALSE)] <- TRUE
  fData(input)$networks_Receptors[which(is.na(match(expressed_genes, recs)) == FALSE)] <- TRUE
  true_ex <- match(expressed_genes_bulk, rownames(fData(input)))
  fData(input)$networks_ligands[-true_ex] <- FALSE
  fData(input)$networks_Receptors[-true_ex] <- FALSE
  for (i in 1:nrow(Receptor_Ligand_Data)) {
    int <- Receptor_Ligand_Data[i,]
    pair <- unlist(strsplit(as.character(int$Pair.Name), split = "_"))
    ind_lig <- match(pair[1], rownames(fData(input)))
    ind_rec <- match(pair[2], rownames(fData(input)))
    lig_state <- fData(input)$networks_ligands[ind_lig]
    rec_state <- fData(input)$networks_Receptors[ind_rec]
    if(is.na(lig_state) == TRUE | is.na(rec_state) == TRUE){
    } else
    if(lig_state == TRUE & rec_state == TRUE){
      fData(input)$networks_ligs_to_receptor[ind_rec] <- paste0(fData(input)$networks_ligs_to_receptor[ind_rec], pair[1], "_")
      fData(input)$networks_receptor_to_ligs[ind_lig] <- paste0(fData(input)$networks_receptor_to_ligs[ind_lig], pair[2], "_")
    }
  }
  fData(input)$networks_ligs_to_receptor <- substr(fData(input)$networks_ligs_to_receptor,1,nchar(fData(input)$networks_ligs_to_receptor)-1)
  fData(input)$networks_ligs_to_receptor[which(fData(input)$networks_ligs_to_receptor == "")] <- NA
  fData(input)$networks_receptor_to_ligs <- substr(fData(input)$networks_receptor_to_ligs,1,nchar(fData(input)$networks_receptor_to_ligs)-1)
  fData(input)$networks_receptor_to_ligs[which(fData(input)$networks_receptor_to_ligs == "")] <- NA
  return(input)
}
