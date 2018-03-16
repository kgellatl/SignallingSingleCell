#' Construct Expression Set Class
#'
#' This function will take an input expression matrix and make an ESC
#'
#' @param input the input data matrix.
#' @export
#' @details
#' This will take an input matrix and create the expression set class that further analysis can be written into.
#' @examples
#' construct_ex_sc(filtered_data)

construct_ex_sc <- function(input){
  HSMM_expr_matrix <- input
  num.columns <- ncol(HSMM_expr_matrix)
  name.columns <- colnames(HSMM_expr_matrix)
  num.rows <- nrow(HSMM_expr_matrix)
  name.rows <- rownames(HSMM_expr_matrix)
  HSMM_sample_sheet <- matrix(nrow = num.columns, ncol= 0)
  rownames(HSMM_sample_sheet) <- c(name.columns)
  HSMM_sample_sheet <- data.frame(HSMM_sample_sheet)
  HSMM_gene_annotation <- matrix(nrow = num.rows, ncol= 0)
  rownames(HSMM_gene_annotation) <- c(name.rows)
  HSMM_gene_annotation <- data.frame(HSMM_gene_annotation)
  pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
  fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
  ex_sc <- ExpressionSet(as.matrix(HSMM_expr_matrix), phenoData = pd, featureData = fd)
  return(ex_sc)
}

