#' tSNE Plot on metadata
#'
#' This will plot information onto a 2d tsne plot
#'
#' @param input the input ex_sc.
#' @param title The title
#' @param color_by What to color points by, either "UMI_sum", or pData categorial variable, ignored if gene is provided
#' @param facet_by What to break the plots by
#' @param ncol How many columns if faceting
#' @param size The size of the points
#' @param colors What colors to utilize for categorial data. Be sure it is of the proper length!
#' @param legend_dot_size Size of dot in legend
#' @param xcol pData column to use for x axis
#' @param ycol pData column to use for y axis
#' @export
#' @details
#' Utilize information stored in pData to control the plot display.
#' @examples
#' plot_tsne_metadata(ex_sc_example, color_by = "UMI_sum", title = "UMI_sum across clusters", facet_by = "Cluster", ncol = 3)

plot_tsne_metadata <- function(input,
                               title = "",
                               color_by,
                               facet_by = NA,
                               ncol = "NA",
                               size = 1.5,
                               colors = NA,
                               theme = "classic",
                               legend_dot_size = 1.5,
                               xcol="x",
                               ycol="y"){

  tmp <- pData(input)
  tmp <- tmp[sample(nrow(tmp)),]
  g <- ggplot(tmp)
  if(theme == "bw") {
    g <- g + theme_bw();
  } else {
    g <- g + theme_classic()
  }

  if(title == ""){
    title <- color_by
    g <- g +  labs(title= title, x = "tSNE[1]", y = "tSNE[2]")
  } else {
    g <- g +  labs(title= title, x = "tSNE[1]", y = "tSNE[2]")
  }
  g <- g + theme(plot.title = element_text(size = 20), axis.title = element_text(size = 10), legend.title = element_text(size = 15), legend.text=element_text(size=10))
  g <- g + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
  if(all(is.na(facet_by)) == FALSE){
    tmp <- pData(input)[c(xcol, ycol)]
    g <- g + geom_point(data = tmp, aes(x=x, y=y), shape = 20, col = "gray", size = size)
  }
  if(class(pData(input)[,color_by]) == "double" || class(pData(input)[,color_by]) == "integer" || class(pData(input)[,color_by]) == "numeric" ){
    g <- g +  geom_point(aes_string(x = xcol, y = ycol, col = color_by), shape = 20, size = size)
    g <- g +  scale_color_gradientn(colours=c('blue', 'red', 'yellow'))
  } else {
    g <- g +  geom_point(aes_string(x = xcol, y = ycol, col = color_by), shape = 20, size = size)
    if(all(is.na(colors)) == FALSE){
      g <- g + scale_color_manual(values = c(colors))
    }
  }
  if(all(is.na(facet_by)) == FALSE){
    if(ncol != "NA"){
      g <- g +  facet_wrap(facets = reformulate(facet_by), ncol = ncol)
    } else {
      g <- g +  facet_wrap(facets = reformulate(facet_by))
    }
  }
  g <- g + guides(colour = guide_legend(override.aes = list(size=legend_dot_size)))
  plot(g)
}

