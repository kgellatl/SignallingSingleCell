#' Plot Density
#'
#' This function will draw a UMI distribution
#'
#' @param input the input data
#' @param val Either "UMI_sum" or a gene name
#' @param color_by The metadata in pData that will be used to break the distributions into groups. Can be left blank.
#' @param statistic Either the "mean" or "median" can be calculated
#' @export
#' @details
#' This will draw a UMI distribution for each variable in the metadat.
#' @examples
#' draw_density(input = ex_sc_example, val = "UMI_sum", color_by = "Cluster", statistic = "mean")

plot_density <- function(input, val, color_by = "NA", statistic = "mean"){
  if(color_by == "NA"){
    if(val == "UMI_sum"){
      gpd <- pData(input)$UMI_sum
      gpd <- as.data.frame(gpd)
      colnames(gpd) <- "UMI_sum"
      gpd$group <- "All"
      clustermean <- aggregate(gpd[,"UMI_sum"], list(group=gpd[,"group"]), statistic)
      clustermean$x = round(clustermean$x)
      ggplot(gpd, aes(UMI_sum)) +
        geom_density(alpha=0.1, lwd=1) +
        scale_x_log10() +
        xlab("UMI_sum per cell") +
        geom_vline(data=clustermean, aes(xintercept=x),
                   linetype="dashed", size=0.5) +
        annotate("text",x=clustermean$x, y=0.1,label=clustermean$x, angle = 45) +
        ggtitle(paste0("log10(", val, ")")) +
        theme_classic() +
        theme(plot.title = element_text(size = 20), axis.title = element_text(size = 10), legend.title = element_text(size = 15), legend.text=element_text(size=10)) +
        theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
    } else {
      gpd <- log2(exprs(input)[val,]+2)-1
      gpd <- as.data.frame(gpd)
      colnames(gpd) <- "gene"
      gpd$group <- "All"
      clustermean <- aggregate(gpd[,"gene"], list(group=gpd[,"group"]), statistic)
      clustermean$x <- round(clustermean$x, 2)
      ggplot(gpd, aes(gene)) +
        geom_density(alpha=0.1, lwd=1)  +
        geom_vline(data=clustermean, aes(xintercept=x),
                   linetype="dashed", size=0.5) +
        annotate("text",x=clustermean$x, y=0.1,label=clustermean$x, angle = 45) +
        xlab(paste0("log2(", val, ")")) +
        theme_classic() +
        theme(plot.title = element_text(size = 20), axis.title = element_text(size = 10), legend.title = element_text(size = 15), legend.text=element_text(size=10)) +
        theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5)) +
        ggtitle(paste0("log2(", val, ")"))
    }
  } else {
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    cols <- gg_color_hue(length(unique(pData(input)[,color_by])))
    if(val == "UMI_sum"){
      gpd <- pData(input)[,c("UMI_sum", color_by)]
      clustermean <- aggregate(gpd[,"UMI_sum"], list(group=gpd[,color_by]), statistic)
      clustermean$x = round(clustermean$x)
      ggplot(gpd, aes_string("UMI_sum", colour=color_by, fill=color_by)) +
        geom_density(alpha=0.1, lwd=1) +
        scale_x_log10() +
        xlab("log10(UMI_sum)") +
        ggtitle(paste0("log10(UMI_sum) across ", color_by)) +
        geom_vline(data=clustermean, aes(xintercept=x, colour=group),
                   linetype="dashed", size=0.5) +
        annotate("text",x=clustermean$x, y=0.1,label=clustermean$x, angle = 45, colour = cols) +
        theme_classic() +
        theme(plot.title = element_text(size = 20), axis.title = element_text(size = 10), legend.title = element_text(size = 15), legend.text=element_text(size=10)) +
        theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
    } else {
      gg_color_hue <- function(n) {
        hues = seq(15, 375, length = n + 1)
        hcl(h = hues, l = 65, c = 100)[1:n]
      }
      cols <- gg_color_hue(length(unique(pData(input)[,color_by])))
      gpd <- pData(input)[,c("UMI_sum", color_by)]
      gpd$gene <- log2(exprs(input)[val,]+2)-1
      clustermean <- aggregate(gpd[,"gene"], list(group=gpd[,color_by]), statistic)
      clustermean$x <- round(clustermean$x, 2)
      colnames(gpd) <- c("UMI_sum", color_by, val)
      ggplot(gpd, aes_string(val, colour=color_by, fill=color_by)) +
        geom_density(alpha=0.1, lwd=1)  +
        geom_vline(data=clustermean, aes(xintercept=x, colour=group),
                   linetype="dashed", size=0.5) +
        annotate("text",x=clustermean$x, y=0.1,label=clustermean$x, angle = 45, colour = cols) +
        xlab(paste0("log2(", val, ")")) +
        theme_classic() +
        theme(plot.title = element_text(size = 20), axis.title = element_text(size = 10), legend.title = element_text(size = 15), legend.text=element_text(size=10)) +
        theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5)) +
        ggtitle(paste0("log2(", val, ") across ", color_by))
    }
  }
}
