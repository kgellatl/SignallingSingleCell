#' This will perform differential expression using edgeR for pairwise comparisons
#'
#'
#'
#' @param input Input expression set
#' @param pd pData, optional, by default set to pData of input
#' @param DEgroup a pData variable
#' @param contrastID the condition value from DEgroup to be used as the contrast
#' @param sizefactor a pData column containing the scran reported size factor for each cell
#' @param lib_size a pData column containing the library size for each cell
#' @param facet_by a pData variable to iterate through and perform DE for a condition within each group
#' @param minCells minimum cell fraction the genes should be expressed in to be tested for DE
#' @param batchID a pData variable to use for model accounting for batch effects
#' @param outdir an output directory
#' @param outsuffix an output file suffix
#' @param pVal pvalue cutoff for reported results, by default reports all results
#' @param contrast a list of contrasts to be used for the DE analysis, default is a two-class comparison, with the second class as the comparison class
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
                       sizefactor = NULL,
                       lib_size = NULL,
                       facet_by,
                       minCells = 0.1,
                       batchID = NULL,
                       outdir = "DEresults/",
                       outsuffix = "DEresults.tsv",
                       pVal = 1,
                       contrast = list(c(-1,1))) {
  # select cells to perform DE
  if (is.null(lib_size)) {
    # if lib_size is null sum up UMI counts
    libsize = colSums(exprs(input)[,rownames(pd)])
    pd$lib_size = libsize
    lib_size = "lib_size"
  }
  for (i in 1:length(unique(pd[,facet_by]))) {
    name = sort(unique(as.character(pd[,facet_by])))[i]
    print(paste0("Performing DE for ", name))
    idx = rownames(pd)[which(pd[,facet_by]==name)]
    cntr = rownames(pd)[which(pd[,facet_by]==name & pd[,DEgroup]==contrastID)]
    # select cells
    z = as.matrix(exprs(input))[,idx]
    if(ncol(z) > 2){
      if (is.null(batchID)) {
        # if batchID is null all cells are in the same batch
        batch = as.factor(rep(1,ncol(z)))
      } else {
        # check if there are cells of each DEgroup per batch
        batch_table = as.matrix(table(pd[idx,batchID],pd[idx,DEgroup]))
        if (length(batch_table[batch_table==0])<2) {
          # get batch factors from batchID column in pData
          batch = as.factor(pd[idx,batchID])
        } else {
          message(sprintf('Warning: not enough cells to include batch in model, setting batches to 1'));flush.console()
          batch = as.factor(rep(1,ncol(z)))
        }
      }
      # create reference group
      groupList = rep(0, times=ncol(z))
      # create contrast group
      groupList[which(colnames(z) %in% cntr)] = 1
      if(length(unique(groupList)) > 1){
        group = factor(groupList)
        z = construct_ex_sc(z)
        pData(z) = pd[idx,,drop=F]
        # perform DE
        tab = edgeRDE(input = z, groups = group, batch = batch, sizefactor = sizefactor, lib_size = lib_size,
                      minCells = minCells, pVal = pVal, contrast = contrast)
        # write DE result for pairwise comparison
        DEtablecntr = tab[['contrast_1']]
        outfile = paste(name, contrastID, outsuffix, sep = "_")
        write.table(DEtablecntr, paste(outdir, outfile, sep = ""), sep = "\t")
      }
    }
  }
}
