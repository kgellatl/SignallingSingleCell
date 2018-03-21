#' This will create a violin plot
#'
#' This will plot a given gene via violin plot
#'
#' @param input The input data
#' @param title The title
#' @param gene if provided will color_by the gene
#' @param color_by a pData variable
#' @param facet_by a pData variable
#' @param ncol How many columns if faceting
#' @param size The size of the points
#' @param colors What colors to utilize for categorial data. Be sure it is of the proper length!
#' @export
#' @details
#' Utilize information stored in pData to control the plot display.
#' @examples
#' plot_violin(input, title = "Actb across clusters", gene = "Actb", color_by = "Timepoint", facet_by = "Cluster", size = 1, ncol = 3)

plot_violin <- function(input, title, color_by, gene, facet_by = "NA", ncol = "NA", size = 1, colors = "NA"){
  if(facet_by == "NA"){
    geneColored1 <- pData(input)[,c("x", "y", color_by)]
    geneColored1 <- cbind(geneColored1, log2(exprs(input)[gene,]+2)-1)
    colnames(geneColored1)<-c("x","y", color_by, gene)
  } else {
    geneColored1 <- pData(input)[,c("x", "y", color_by, facet_by)]
    geneColored1 <- cbind(geneColored1, log2(exprs(input)[gene,]+2)-1)
    colnames(geneColored1)<-c("x","y", color_by, facet_by, gene)
  }
  g <- ggplot(geneColored1)
  g <- g + theme_classic()
  g <- g + labs(title= title, x = color_by, y = gene)
  g <- g + theme(plot.title = element_text(size = 20), axis.title = element_text(size = 10), legend.title = element_text(size = 15), legend.text=element_text(size=10))
  g <- g + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
  g <- g + geom_violin(aes_string(x=color_by, y=gene, fill = color_by), trim = T)
  g <- g + geom_jitter(aes_string(x=color_by, y=gene, col = color_by), width = 0.2, alpha=0.5, size = size)
  g <- g + stat_summary(aes_string(x=color_by, y=gene), fun.y=mean, geom="point", size=3, color="black")
  if(facet_by != "NA"){
    if(ncol != "NA"){
      g <- g +  facet_wrap(facets = reformulate(facet_by), ncol = ncol)
    } else {
      g <- g +  facet_wrap(facets = reformulate(facet_by))
    }
  }
  return(g)
}
