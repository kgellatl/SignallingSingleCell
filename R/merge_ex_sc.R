#' This will merge pData and fData
#'
#' This function will take an input ex_sc and produce a subsetted one preserving all metadata
#'
#' @param input the expression set to merge into
#' @param subset_input the expression set to be merged (contains a SUBSET of input1)
#' @param merge_key this key will be appended to the colnames of x,y, and cluster
#' @param include defaults to all
#' @export
#' @details
#' If you have subsetted data and want to merge the tSNE coordinates and cluster information into the master data, this will do so.
#' @examples
#' merged_ex_sc <- merge_ex_sc(input1, input2, merge_key)

merge_ex_sc <- function(input, subset_input, merge_key){

}

