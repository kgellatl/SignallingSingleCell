#' This will setset your expression set by some variable in pData
#'
#' This function will take an input ex_sc and produce a subsetted one preserving all metadata
#'
#' @param input the input ex_sc
#' @param variable the variable of pData to select on
#' @param select the value or values to subset
#' @export
#' @details
#' By providing both a variable, and a selecting argument, only cells of of a certain condition will be selected
#' @examples
#' subset_ex_sc(ex_sc_example, variable = "Cluster", select = c("Cluster1"))

subset_ex_sc <- function(input, variable, select){
 dat <- pData(input)
 ind <- match(variable, colnames(dat))
 val <- dat[,ind]
 selected <- c()
 for(i in 1:length(select)){
   int1 <- select[i]
   ind2 <- grep(int1, val)
   selected <- c(selected, ind2)
 }
 cells <- rownames(dat)[selected]
 ex_sc <- input[,cells]
 return(ex_sc)
}

