#' This will create a violin plot
#'
#' This will plot a given gene via violin plot
#'
#' @param gene Will only plot data for this gene, if the gene paramter is a list the violin plot will display the sum of the expression for all genes in the list
#' @param input Input expression set
#' @param sampleID pData variable to group by on the x axis
#' @param facetID a pData variable
#' @param colourID a pData variable
#' @param cols Personalized colour vector, length should equal number of groups in colourID
#' @param subsetID a pData variable to subset the data by
#' @param subsetName a value from the subsetID variable specified
#' @param facet_scale ggplot2 parameter for plot scales defaults to fixed
#' @export
#' @details
#' Utilize information stored in pData to control the plot display.
#' @examples
#' plotViolin("Actb", input, sampleID = "sample", facetID = "subtype", colourID = "genotype", cols = NULL, subsetID = "celltype", subsetName = "Neurons")
###
plotViolin = function(gene, input, sampleID, facetID, colourID, cols = NULL, subsetID = NA, subsetName = NA, facet_scale = "fixed") { # gene, HSMM object, column with sample ID, column for facet, column for colour, colours, column to subset, name to subset
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  ###
  label.n = function(x) {
    return(c(y=-1, label=length(x)))
  }
  ###
  if (!is.null(subsetID) & !is.na(subsetID)) {
    input = input[,which(pData(input)[,subsetID]==subsetName)]
  }
  inputM = as.data.frame(exprs(input[rownames(input)%in%gene,]))
  geneName = gene;
  if(length(gene)>1) {
    sum = apply(inputM, MARGIN = 2, FUN = sum);
    inputM[1,] = sum;
    inputM = inputM[-(2:length(rownames(inputM))),]
    geneName = gsub(", ","_",toString(gene));
    rownames(inputM) = geneName;
  }
  inputM = log2(inputM+1)
  merged = cbind(pData(input),t(inputM))
  if (is.null(cols) | length(cols)==0) {
    cols = gg_color_hue(length(unique(merged[,colourID])))
  }

  title = geneName;
  if (!is.null(subsetName) & !is.na(subsetName) ) {
    title = paste(subsetName, geneName, sep=" - ");
  }
  ggplot(merged, aes_string(sampleID, geneName, colour = colourID)) +
    facet_wrap(reformulate(facetID), scales = facet_scale) +
    geom_jitter(width = 0.2, alpha=0.5) +
    geom_violin(fill = NA) +
    stat_summary(data = merged, aes_string(x = sampleID, y = geneName), fun.data = label.n, fun.y = "mean", colour = "black", geom = "text", size=3) +
    stat_summary(data = merged, aes_string(x = sampleID, y = geneName), fun.y = mean, colour = "black", geom = "point", size=3, shape = 18) +
    theme_bw() +
    scale_colour_manual(values=cols) +
    ggtitle(title)
}
