#' Calculate Library Size
#'
#' This function will count the number of reads / UMIs per cell
#'
#' @param input the input DIRECTORY
#' @export
#' @details
#' This will calculate the total number of UMIs on a per cell basis.
#' @examples
#' ex_sc_example <- calc_libsize(input = ex_sc_example)

construct_ex_sc_10x <- function(bc_path, long_genes = T, feature_assay = F){

  samples <- list.files(bc_path)

  for (i in 1:length(samples)) {
    int_sample <- samples[i]
    matrix_dir <- paste0(bc_path, int_sample, "/")
    barcode.path <- paste0(matrix_dir, "barcodes.tsv")
    features.path <- paste0(matrix_dir, "features.tsv")
    matrix.path <- paste0(matrix_dir, "matrix.mtx")
    mat <- readMM(file = matrix.path)
    feature.names = read.delim(features.path,
                               header = FALSE,
                               stringsAsFactors = FALSE)
    barcode.names = read.delim(barcode.path,
                               header = FALSE,
                               stringsAsFactors = FALSE)
    colnames(mat) = gsub("-1", "", barcode.names$V1)
    colnames(mat) <- paste0(colnames(mat), "-", int_sample)
    rownames(mat) = feature.names$V2
    if( i == 1){
      master_data <- mat
    } else {
      master_data <- cbind(master_data, mat)
    }
  }

  #####
  dups <- names(which(table(rownames(master_data)) > 1))
  for (i in 1:length(dups)) {
    int_dup <- dups[i]
    ind <- grep(paste0("^",int_dup,"$"), rownames(master_data))
    master_data[ind[1],] <- apply(master_data[ind,],2,sum)
    master_data <- master_data[-ind[2],]
  }

  ex_sc <- construct_ex_sc(master_data)

  if(long_genes){
    ind <- match(rownames(ex_sc), feature.names$V2)
    fData(ex_sc)$"gene_long" <- feature.names$V1[ind]
  }

  if(feature_assay){
    ind <- match(rownames(ex_sc), feature.names$V2)
    fData(ex_sc)$"feature_type" <- feature.names$V3[ind]

  }

  ex_sc$Sample <- matrix(unlist(strsplit(colnames(ex_sc), split = "-")),byrow = T, ncol = 2)[,2]

  return(ex_sc)
}
