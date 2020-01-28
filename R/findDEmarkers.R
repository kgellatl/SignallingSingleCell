#' This will perform differential expression using edgeR to find markers for each cluster or group compared to all others
#'
#'
#'
#' @param input Input expression set
#' @param pd pData, optional, by default set to pData of input
#' @param DEgroup a pData variable
#' @param batchID a pData variable to use for model accounting for batch effects
#' @param sizefactor a pData column containing the scran reported size factor for each cell
#' @param lib_size a pData column containing the library size for each cell
#' @param outdir an output directory
#' @param outsuffix an output file suffix
#' @param minCells minimum cell fraction the genes should be expressed in to be tested for DE
#' @param pVal pvalue cutoff for reported results, by default reports all results
#' @param contrast a list of contrasts to be used for the DE analysis default is a two-class comparison, with the second class as the comparison class
#' @export
#' @details
#' Utilize information stored in pData to control the DE performed
#' @examples
#' findDEmarkers(input, pd = selected_table, DEgroup = "cluster", batchID = "beads", sizefactor="sizefactor", lib_size="lib_size")
###
findDEmarkers = function(input,
                         pd = pData(input),
                         DEgroup,
                         sizefactor = NULL,
                         lib_size = NULL,
                         batchID = NULL,
                         outdir = "DEmarkers/",
                         outsuffix = "DEmarkers.tsv",
                         minCells = 0.01,
                         pVal = 1,
                         contrast = list(c(-1,1))) {
	if (is.null(lib_size)) {
		# if lib_size is null sum up UMI counts
		lib_size = colSums(exprs(input))
	} else {
		lib_size = pd[,lib_size]
	}
	if (is.null(batchID)) {
		# if batchID is null all cells are in the same batch
		batch = as.factor(rep(1,ncol(input)))
	} else {
		# get batch factors from batchID column in pData
		batch = as.factor(pd[,batchID])
	}

	for (i in 1:length(unique(pd[,DEgroup]))) {
		name = unique(as.character(pd[,DEgroup]))[i]
		message(paste("Working on group: ", name, sep=""));flush.console()
		idx = rownames(pd)[which(pd[,DEgroup]==name)]
		groupList = rep(0, times=ncol(input))   # all cells are reference
		groupList[which(colnames(input) %in% idx)] = 1 # cells that match id are used as contrast
		group = factor(groupList)
		tab = edgeRDE(input, group, batch, sizefactor, lib_size, minCells, pVal, contrast)
		#Handle results with multiple contrasts
		message('Writing results to output file');flush.console()
		if (length(tab) == 1){
			DEtableclsall = tab[['contrast_1']]
			outfile = paste(name, outsuffix, sep = "_")
			write.table(DEtableclsall, paste(outdir, outfile, sep = ""), sep = "\t")
		} else {
			for (t in length(tab)){
				DEtableclsall = tab[[t]]
				outfile = paste(name, t, outsuffix, sep = "_")
				write.table(DEtableclsall, paste(outdir, outfile, sep = ""), sep = "\t")
			}
		}
	}
}
