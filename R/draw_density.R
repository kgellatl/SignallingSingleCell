#' Draw Density
#'
#' This function will draw a UMI distribution
#'
#' @param input the input data
#' @param val Either "UMIs" or a gene name
#' @param var1 The metadata in pData that will be used to break the distributions into groups. Can be left blank.
#' @param statistic Either the "mean" or "median" can be calculated
#' @export
#' @details
#' This will draw a UMI distribution for each variable in the metadat.
#' @examples
#' draw_density(input = ex_sc_example, val = "UMIs", var1 = "Cluster", statistic = "mean")

draw_density <- function(input, val, var1 = "NA", statistic = "mean"){
  if(var1 == "NA"){
    if(val == "UMIs"){
      gpd <- pData(input)$UMI_sum
      gpd <- as.data.frame(gpd)
      colnames(gpd) <- "UMI_sum"
      gpd$group <- "All"
      clustermean <- aggregate(gpd[,"UMI_sum"], list(group=gpd[,"group"]), statistic)
      clustermean$x = round(clustermean$x)
      ggplot(gpd, aes_string(UMI_sum)) +
        geom_density(alpha=0.1, lwd=1) +
        scale_x_log10() +
        xlab("UMIs per cell") +
        geom_vline(data=clustermean, aes(xintercept=x),
                   linetype="dashed", size=0.5) +
        annotate("text",x=clustermean$x, y=0.1,label=clustermean$x, angle = 45) +
        ggtitle(paste0("log2(", val, ")")) +
        theme_classic()
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
        ggtitle(paste0("log2(", val, ")"))
    }
  }else {
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    cols <- gg_color_hue(length(unique(pData(input)[,var1])))
    if(val == "UMIs"){
      gpd <- pData(input)[,c("UMI_sum", var1)]
      clustermean <- aggregate(gpd[,"UMI_sum"], list(group=gpd[,var1]), statistic)
      clustermean$x = round(clustermean$x)
      ggplot(gpd, aes_string(UMI_sum, colour=var1, fill=var1)) +
        geom_density(alpha=0.1, lwd=1) +
        scale_x_log10() +
        xlab("log10(UMIs)") +
        ggtitle(paste0("log10(UMIs) across ", var1)) +
        geom_vline(data=clustermean, aes(xintercept=x, colour=group),
                   linetype="dashed", size=0.5) +
        annotate("text",x=clustermean$x, y=0.1,label=clustermean$x, angle = 45, colour = cols) +
        theme_classic()
    } else {
      gg_color_hue <- function(n) {
        hues = seq(15, 375, length = n + 1)
        hcl(h = hues, l = 65, c = 100)[1:n]
      }
      cols <- gg_color_hue(length(unique(pData(input)[,var1])))
      gpd <- pData(input)[,c("UMI_sum", var1)]
      gpd$gene <- log2(exprs(input)[val,]+2)-1
      clustermean <- aggregate(gpd[,"gene"], list(group=gpd[,var1]), statistic)
      clustermean$x <- round(clustermean$x, 2)
      colnames(gpd) <- c("UMI_sum", "Cluster", val)
      ggplot(gpd, aes_string(val, colour=var1, fill=var1)) +
        geom_density(alpha=0.1, lwd=1)  +
        geom_vline(data=clustermean, aes(xintercept=x, colour=group),
                   linetype="dashed", size=0.5) +
        annotate("text",x=clustermean$x, y=0.1,label=clustermean$x, angle = 45, colour = cols) +
        xlab(paste0("log2(", val, ")")) +
        theme_classic() +
        ggtitle(paste0("log2(", val, ") across ", var1))
    }
  }
}
