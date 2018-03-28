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

plot_density_ridge <- function(input, val, title, color_by){
  dat <- pData(input)
  if(val != "UMI_sum"){
    dat <- cbind(dat, log2(exprs(input)[val,]+2)-1)
    colnames(dat) <- c(colnames(dat[2:ncol(dat)-1]), val)
  }
  g <- ggplot(dat, aes_string(x = val, y = color_by, col = color_by))
  g <- g + ggridges::geom_density_ridges(alpha=0.1, lwd=1)
  g <- g + theme_classic()
  g <- g + theme(plot.title = element_text(size = 20), axis.title = element_text(size = 10), legend.title = element_text(size = 15), legend.text=element_text(size=10))
  g <- g + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
  g <- g + ggtitle(title)
  # g <- g + scale_x_log10()
  g <- g + xlab(val)
  print(g)
}
