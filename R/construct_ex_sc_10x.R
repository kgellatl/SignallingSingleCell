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

construct_ex_sc_10x <- function(bc_path, long_genes = T, feature_assay = T){

  samples <- list.files(bc_path)

  for (i in 1:length(samples)) {
    int_sample <- samples[i]
    matrix_dir <- paste0(bc_path, int_sample, "/")
    barcode.path <- paste0(matrix_dir, "barcodes.tsv")
    features.path <- paste0(matrix_dir, "features.tsv")
    matrix.path <- paste0(matrix_dir, "matrix.mtx")
    mat <- Matrix::readMM(file = matrix.path)
    print(dim(mat))
    feature.names = read.delim(features.path,
                               header = FALSE,
                               stringsAsFactors = FALSE)
    feature.names$V2 <- make.unique(feature.names$V2)
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

  master_data <- as.matrix(master_data)

  if(feature_assay){

    ind_fd <- which(feature.names$V3 != "Gene Expression")
    ind_exp <- which(feature.names$V3 == "Gene Expression")

    expr_mat <- master_data[ind_exp,]
    fDat <- master_data[ind_fd,]
    ex_sc <- construct_ex_sc(expr_mat)
    ex_sc$Sample <- matrix(unlist(strsplit(colnames(ex_sc), split = "-")),byrow = T, ncol = 2)[,2]
    pData(ex_sc) <- cbind(as.data.frame(t(fDat)), pData(ex_sc))

  } else {
    ex_sc <- construct_ex_sc(master_data)
    ex_sc$Sample <- matrix(unlist(strsplit(colnames(ex_sc), split = "-")),byrow = T, ncol = 2)[,2]
  }

  if(long_genes){
    ind <- match(rownames(ex_sc), feature.names$V2)
    fData(ex_sc)$"gene_long" <- feature.names$V1[ind]
  }

  return(ex_sc)
}
