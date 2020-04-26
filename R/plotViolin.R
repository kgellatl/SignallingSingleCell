#' This will create a violin plot
#'
#' This will plot a given gene via violin plot
#'
#' @param gene Will only plot data for this gene, if the gene paramter is a list the violin plot will display the sum of the expression for all genes in the list
#' @param input Input expression set
#' @param sampleID pData variable to group by on the x axis
#' @param colourID a pData variable to color by
#' @param facetID an optional pData variable to plot by
#' @param cols Personalized colour vector, length should equal number of groups in colourID
#' @param subsetID a pData variable to subset the original data by
#' @param subsetName a value from the subsetID variable specified, this will select only these cells to plot
#' @param facet_scale ggplot2 parameter for plot scales defaults to fixed
#' @param mean_text ggplot2 parameter for adding the mean expression per group defaults to False
#' @param fraction_text ggplot2 parameter for adding the fraction of cells with expression greater than 0 for that gene defaults to False
#' @param type plot type, default is violin, can also plot density and cdf distribution
#' @param fudge a pseudocount value to avoid dropping cells with zero counts
#' @export
#' @details
#' Utilize information stored in pData to control the plot display.
#' @examples
#' plotViolin("Actb", input, sampleID = "sample", facetID = "subtype", colourID = "genotype", cols = NULL, subsetID = "celltype", subsetName = "Neurons", facet_scale = "fixed")
###
plotViolin = function(gene,
                      input,
                      sampleID,
                      colourID,
                      facetID = NULL,
                      cols = NULL,
                      subsetID = NA,
                      subsetName = NA,
                      facet_scale = "fixed",
                      mean_text = F,
                      fraction_text = F,
                      type = "violin",
                      legend = T,
                      fudge = 0)
{
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  ###
  label.n = function(c) {
    return(c(y=-0.5, label=length(c)))
  }
  meanSC = function(x){
    return(c(y = -1, label = round(10^(mean(x)), 2)))
  }
  fracSC = function(f){
    return(c(y = -1.5, label = round(10^(mean(f)), 2)))
  }
  ###
  if (!is.null(subsetID) & !is.na(subsetID)) {
    input = input[,which(pData(input)[,subsetID]%in%subsetName)]
  }
  inputM = as.data.frame(exprs(input[rownames(input)%in%gene,]))
  inputM = inputM + fudge
  geneName = gene;
  if(length(gene)>1) {
    sum = apply(inputM, MARGIN = 2, FUN = sum);
    inputM[1,] = sum;
    inputM = inputM[-(2:length(rownames(inputM))),]
    geneName = gsub(", ","_",toString(gene));
    rownames(inputM) = geneName;
  }
  merged = cbind(pData(input),t(inputM))
  # calculate fraction of cells per sampleID that express the gene
  fracCells = "frac"
  if (!is.null(facetID)) {
    frac = table(merged[merged[,gene]>fudge,sampleID],merged[merged[,gene]>fudge,facetID])/table(merged[,sampleID],merged[,facetID])
    frac = as.data.frame(frac)
    colnames(frac) = c(sampleID,facetID,"frac")
    merged = merge(merged, frac, by=c(sampleID,facetID))
  } else {
    #frac = table(merged[merged[,gene]>fudge,sampleID])/table(merged[,sampleID])
    #frac = as.data.frame(frac)
    #colnames(frac) = c(sampleID,"frac")
    #merged = merge(merged, frac, by=sampleID)
  }
  if (is.null(cols) | length(cols)==0) {
    cols = gg_color_hue(length(unique(merged[,colourID])))
  }
  title = geneName;
  if (!is.null(subsetName) & !is.na(subsetName) ) {
    title = paste(subsetName, geneName, sep=" - ");
  }
  if (type == "violin") {
    p = ggplot(merged, aes_string(sampleID, geneName, colour = colourID)) +
      geom_jitter(width = 0.2, size=0.5) +
      geom_violin(fill = NA) +
      stat_summary(data = merged, aes_string(x = sampleID, y = geneName), fun.data = label.n, fun.y = "mean", colour = "black", geom = "text", size=3) +
      stat_summary(data = merged, aes_string(x = sampleID, y = geneName), fun.y = mean, colour = "black", geom = "point", size=3, shape = 18) +
      theme_bw() +
      scale_colour_manual(values=cols) +
      scale_y_continuous(trans = scales::pseudo_log_trans()) +
      ggtitle(title)
    if (mean_text == T) {
      p = p + stat_summary(data = merged, aes_string(x = sampleID, y = geneName), fun.data = meanSC, fun.y = "mean", colour = "black", geom = "text", size=3)
    }
    if (fraction_text == T) {
      p = p + stat_summary(data = merged, aes_string(x = sampleID, y = fracCells), fun.data = fracSC, fun.y = "max", colour = "red", geom = "text", size=3)
    }
    if (!is.null(facetID)) {
      p = p + facet_wrap(reformulate(facetID), scales = facet_scale)
    }
    if (legend == F) {
      p = p + theme(legend.position = "none")
    }
  }
  if (type == "density") {
    p = ggplot(merged, aes_string(geneName, colour = colourID)) +
      geom_density(fill = NA) +
      theme_bw() +
      scale_colour_manual(values=cols) +
      scale_x_log10() +
      ggtitle(title)
    if (!is.null(facetID)) {
      p = p + facet_wrap(reformulate(facetID), scales = facet_scale)
    }
  }
  if (type == "cdf") {
    p = ggplot(merged, aes_string(geneName, colour = colourID)) +
      stat_ecdf() +
      theme_bw() +
      scale_colour_manual(values=cols) +
      scale_x_log10() +
      ggtitle(title)
    if (!is.null(facetID)) {
      p = p + facet_wrap(reformulate(facetID), scales = facet_scale)
    }
  }
  return(p)
}
