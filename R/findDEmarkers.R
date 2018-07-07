#' This will perform differential expression (DE) using edgeR to find markers for each cluster or group compared to all others
#'
#'
#'
#' @param input Input expression set
#' @param pd pData (optional) by default set to pData(input)
#' @param DEgroup a pData variable
#' @param batchID a pData variable to use for model accounting for batch effects
#' @param sizefactor a pData column containing the scran reported size factor for each cell
#' @param lib_size a pData column containing the library size for each cell
#' @param outdir an output directory (default="DEmarkers/")
#' @param outsuffix an output file suffix (default="DEmarkers.tsv")
#' @param minCells minimum cell fraction the genes should be expressed in to be tested for DE (default = 10%)
#' @param pVal pvalue cutoff for reported results (default=1, reports all results)
#' @param contrast a list of contrasts to be used for the DE analysis (default is a two-class comparison, with the second class as the comparison class)
#' @export
#' @details
#' Utilize information stored in pData to control the DE performed
#' @examples
#' findDEmarkers(input, pd = selected_table, DEgroup = "cluster", batchID = "beads", sizefactor="sizefactor", lib_size="lib_size")
#' ###
findDEmarkers = function(input,
                         pd = pData(input),
                         DEgroup,
                         batchID,
                         sizefactor,
                         lib_size,
                         outdir = "DEmarkers/",
                         outsuffix = "DEmarkers.tsv",
                         minCells = 0.1,
                         pVal = 1,
                         contrast = list(c(-1,1))) {
  for (i in 1:length(unique(pd[,DEgroup]))) {
    name = unique(as.character(pd[,DEgroup]))[i]
    idx = rownames(pd)[which(pd[,DEgroup]==name)]
    z = as.matrix(exprs(input[,rownames(pd)]))
    batch = as.factor(pd[,batchID])
    groupList = rep(0, times=ncol(z))   # all cells are reference
    groupList[which(colnames(z) %in% idx)] = 1 # cells that match id are used as contrast
    group = factor(groupList)
    z = construct_ex_sc(z)
    pData(z) = pd[colnames(z),]
    tab = edgeRDE(z, group, batch, sizefactor, lib_size, minCells)
    DEtableclsall = tab[['contrast_1']]
    outfile = paste(name, outsuffix, sep = "_")
    write.table(DEtablecntr, paste(outdir, outfile, sep = ""), sep = "\t")
  }
}
