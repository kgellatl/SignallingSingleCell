#' This will perform differential expression (DE) using edgeR for pairwise comparisons
#'
#'
#'
#' @param input Input expression set
#' @param pd pData (optional) by default set to pData(input)
#' @param DEgroup a pData variable
#' @param contrastID the pData group to set as reference (example = "WT")
#' @param sizefactor a pData column containing the scran reported size factor for each cell
#' @param lib_size a pData column containing the library size for each cell
#' @param facet_by a pData variable to iterate through and perform DE for a condition within each group
#' @param minCells minimum cell fraction the genes should be expressed in to be tested for DE (default = 10%)
#' @param batchID a pData variable to use for model accounting for batch effects
#' @param outdir an output directory (default="DEresults/")
#' @param outsuffix an output file suffix (default="DEresults.tsv")
#' @param pVal pvalue cutoff for reported results (default=1, reports all results)
#' @param contrast a list of contrasts to be used for the DE analysis (default is a two-class comparison, with the second class as the comparison class)
#' @export
#' @details
#' Utilize information stored in pData to control the DE performed
#' @examples
#' findDEgenes(input, pd = selected_table, DEgroup = "genotype", batchID = "beads", facet_by = "cluster", sizefactor="sizefactor", lib_size="lib_size")
###
findDEgenes = function(input,
                       pd = pData(input),
                       DEgroup,
                       contrastID,
                       sizefactor,
                       lib_size,
                       facet_by,
                       minCells = 0.1,
                       batchID = NULL,
                       outdir = "DEresults/",
                       outsuffix = "DEresults.tsv",
                       pVal = 1,
                       contrast = list(c(-1,1))) {
  # select cells to perform DE
  for (i in 1:length(unique(pd[,facet_by]))) {
    name = unique(as.character(pd[,facet_by]))[i]
    idx = rownames(pd)[which(pd[,facet_by]==name)]
    cntr = rownames(pd)[which(pd[,facet_by]==name & pd[,DEgroup]==contrastID)]
    # select cells
    z = as.matrix(exprs(input))[,idx]
    if (is.null(batchID)) {
      # if batchID is null all cells are in the same batch
      batch = as.factor(rep(1,ncol(input)))
    } else {
      # get batch factors from batchID column in pData
      batch = as.factor(pd[idx,batchID])
    }
    # create reference group
    groupList = rep(0, times=ncol(z))
    # create contrast group
    groupList[which(colnames(z) %in% cntr)] = 1
    group = factor(groupList)
    z = construct_ex_sc(z)
    pData(z) = pd[idx,]
    # perform DE
    tab = edgeRDE(z, group, batch, sizefactor, lib_size)
    # write DE result for pairwise comparison
    DEtablecntr = tab[['contrast_1']]
    outfile = paste(name, contrastID, outsuffix, sep = "_")
    write.table(DEtablecntr, paste(outdir, outfile, sep = ""), sep = "\t")
  }
}
