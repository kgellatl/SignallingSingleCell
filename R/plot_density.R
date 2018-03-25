#' Plot Density
#'
#' This function will draw a distribution for a given value
#'
#' @param input the input ex_sc
#' @param val Either "UMI_sum" or a gene name
#' @param title The title of the plot
#' @param color_by The metadata in pData that will be used to break the distributions into groups. Can be left blank.
#' @param statistic Either the "mean" or "median" can be calculated
#' @export
#' @details
#' This will draw a UMI distribution for a given value, broken into groups if color_by != "NA"
#' @examples
#' draw_density(input = ex_sc_example, val = "UMI_sum", color_by = "Cluster", statistic = "mean")

plot_density <- function(input, val, title, color_by = "NA", statistic = "mean"){
  gg_color_hue <- function(n) { #ggplot color selection tool. For the mean line
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  dat <- pData(input)
  if(val != "UMI_sum"){
    dat <- cbind(dat, log2(exprs(input)[val,]+2)-1)
    colnames(dat) <- c(colnames(dat[2:ncol(dat)-1]), val)
  }
  if(color_by == "NA"){
    dat$group <- "All"
    clustermean <- aggregate(dat[,val], list(group=dat[,"group"]), statistic)
    clustermean$x = round(clustermean$x)
  } else {
    clustermean <- aggregate(dat[,val], list(group=dat[,color_by]), statistic)
    clustermean$x = round(clustermean$x, 2)
  }
  if(color_by == "NA"){
    g <- ggplot(dat, aes_string(val))
    g <- g + geom_vline(data=clustermean, aes(xintercept=x), linetype="dashed", size=0.5)
  } else {
    g <- ggplot(dat, aes_string(val, colour=color_by, fill=color_by))
    g <- g + geom_vline(data=clustermean, aes(xintercept=x, colour=group), linetype="dashed", size=0.5)
  }
  g <- g + theme_classic()
  g <- g + theme(plot.title = element_text(size = 20), axis.title = element_text(size = 10), legend.title = element_text(size = 15), legend.text=element_text(size=10))
  g <- g + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
  g <- g + ggtitle(title)
  g <- g + geom_density(alpha=0.1, lwd=1)
  g <- g + scale_x_log10()
  g <- g + xlab(val)
  if(color_by == "NA"){
    g <- g + annotate("text",x=clustermean$x, y=0.1,label=clustermean$x, angle = 45)
  } else {
    cols <- gg_color_hue(length(unique(pData(input)[,color_by])))
    g <- g + annotate("text",x=clustermean$x, y=0.1,label=clustermean$x, angle = 45, colour = cols)
  }
  return(g)
}
