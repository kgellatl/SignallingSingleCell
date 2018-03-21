#' Create a Scatter Plot
#'
#' This will plot information onto a 2d scatter (Gene 1 v Gene 2)
#'
#' @param input The input data
#' @param title The plot title
#' @param gene1 The first gene for scatter plot
#' @param gene2 The second gene for scatter plot
#' @param color_by What to color points by, either "UMI_sum" or pData categorial variable
#' @param facet_by What to break the plots by
#' @param ncol How many columns if faceting
#' @param size The size of the points
#' @param colors What colors to utilize for categorial data. Be sure it is of the proper length!
#' @export
#' @details
#' Utilize information stored in pData to control the plot display.
#' @examples
#' plot_scatter(input = ex_sc_example, title = "Plot", gene1 = "Ccl22", gene2 = "Ccl5", color_by = "Cluster", facet_by = "Timepoint")

plot_scatter <- function(input, gene1, title, gene2, color_by, facet_by = "NA", ncol = "NA", size = 2, colors = "NA"){
  dat <- as.data.frame(t(exprs(input)[c(gene1, gene2),]))
  dat[,1] <- log2(dat[,1]+2)-1
  dat[,2] <- log2(dat[,2]+2)-1
  dat$color_by <- pData(input)[,color_by]
  colnames(dat) <- c(gene1, gene2, color_by)
  if (facet_by != "NA"){
    dat$facet_by <- pData(input)[,facet_by]
    colnames(dat) <- c(gene1, gene2, color_by, facet_by)
  }
  g <- ggplot(dat)
  g <- g + ggtitle(paste0(gene1, " vs ", gene2, " scatter plot"))
  g <- g + theme_classic()
  g <- g + labs(x = paste0("log2(", gene1, ")"), y = paste0("log2(", gene2, ")"))
  g <- g + theme(plot.title = element_text(size = 20), axis.title = element_text(size = 10), legend.title = element_text(size = 15), legend.text=element_text(size=10))
  g <- g + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
  g <- g + geom_point(aes_string(x = gene1, y = gene2, col = color_by), size = size, alpha = 0.75)
  if(colors == "NA"){
  } else {
    g <- g + scale_color_manual(values = c(colors))
  }
  if(facet_by != "NA"){
    if(ncol != "NA"){
      g <- g +  facet_wrap(facets = reformulate(facet_by), ncol = ncol)
    } else {
      g <- g +  facet_wrap(facets = reformulate(facet_by))
    }
  }
  g <- g + labs(title= title)
  return(g)
}



