#' Plot Density
#'
#' This function will draw a distribution for a given value
#'
#' @param input the input ex_sc
#' @param val Either "UMI_sum" or a gene name
#' @param title The title of the plot
#' @param color_by The metadata in pData that will be used to break the distributions into groups. Can be left blank.
#' @importFrom ggridges geom_density_ridges
#' @export
#' @details
#' This will draw a UMI distribution for a given value, broken into groups if color_by != "NA"
#' @examples
#' draw_density(input = ex_sc_example, val = "UMI_sum", color_by = "Cluster", statistic = "mean")

plot_density_ridge <- function(input, val, title = "", color_by){
  dat <- pData(input)
  ind1 <- grep(val, colnames(dat))
  if(length(ind1) == 0){
    dat <- cbind(dat, log2(exprs(input)[val,]+2)-1)
    colnames(dat) <- c(colnames(dat[2:ncol(dat)-1]), val)
  }
  g <- ggplot(dat, aes_string(x = val, y = color_by, col = color_by, fill = color_by))
  g <- g + ggridges::geom_density_ridges(scale = 1, alpha = 0.25)
  g <- g + theme_classic()
  g <- g + theme(plot.title = element_text(size = 20), axis.title = element_text(size = 10), legend.title = element_text(size = 15), legend.text=element_text(size=10))
  g <- g + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
  if(title == ""){
    title <- val
    g <- g + ggtitle(title)
  } else {
    g <- g + ggtitle(title)
  }
  g <- g + xlab(val)
  if(length(ind1) == 1){
    g <- g + scale_x_log10()
  }
  return(g)
}
