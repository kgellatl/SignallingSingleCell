#' This will run edgeR to find differentially expressed genes. Use the wrapping functions findDEgenes or findDEmarkers to use pData columns and get the proper input
#'
#'
#'
#' @param input Input expression set
#' @param groups an integer vector with the group for each cell
#' @param batch a vector with the batch for each cell, if you don't need batches create a vector with the same id for all cells
#' @param sizefactor a pData column containing the scran reported size factor for each cell
#' @param lib_size a pData column containing the library size for each cell
#' @param minCells minimum cell fraction the genes should be expressed in to be tested for DE
#' @param pVal pvalue cutoff for reported results, by default reports all results
#' @param contrast a list of contrasts to be used for the DE analysis, default is a two-class comparison, with the second class as the comparison class
#' @export
#' @details
#' Utilize information stored in pData to control the DE performed
#' The result is a list with the contrast results, access each element for the DE table
#' @examples
#' myDEresults = edgeRDE(input, groups = mygroup, batch = beads, sizefactor = "sizefactor", lib_size = "lib_size")
###
edgeRDE <- function(input,
                    groups,
                    batch,
                    sizefactor,
                    lib_size,
                    minCells = 0.1,
                    pVal = 1,
                    contrast = list(c(-1,1))) {
  ## remove any zero-variance genes
  rvar <- apply(exprs(input),1,var)
  idx.keep <- which(rvar>0)
  input <- input[idx.keep,]
  ## remove genes expressed in less than input% of the cells defined in argument minCells (default = 10%)
  gfrac <- apply(exprs(input),1,function(a) length(a[a>0])/length(a))
  idx.keep <- which(gfrac>minCells)
  input <- input[idx.keep,]
  ## make DGElist
  y <- edgeR::DGEList(counts=exprs(input), group=groups)
  y$samples$norm.factors = pData(input)[,sizefactor]/pData(input)[,lib_size]
  # set batch as coefficient if more than one batch is specified
  if (length(unique(batch))>1) {
    design <- model.matrix(~0+batch+group, data=y$samples)
  } else {
    design <- model.matrix(~0+group, data=y$samples)
  }
  message('Estimating dispersion...');flush.console()
  y <- edgeR::estimateDisp(y, design)
  message('Doing likelihood ratio fit...');flush.console()
  fit <- edgeR::glmFit(y, design)
  tab <- list()
  if (length(unique(batch))>1) {
    lrt <- edgeR::glmLRT(fit)
    tab[['contrast_1']] <- as.data.frame(edgeR::topTags(lrt, p.value=pVal, n=Inf, sort.by='logFC'))
  } else {
    for (j in 1:length(contrast)) {
      cStr <- sprintf('contrast_%d', j)
      message(sprintf('DE analysis for contrast %s...',cStr));flush.console()
      lrt <- edgeR::glmLRT(fit, contrast=contrast[[j]])
      tab[[cStr]] <- as.data.frame(edgeR::topTags(lrt, p.value=pVal, n=Inf, sort.by='logFC'))
    }
  }
  return(tab)
}
