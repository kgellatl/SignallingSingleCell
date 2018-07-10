#' Plots a heatmap
#'
#' Makes either a single cell or bulk heatmap,
#'
#' @param input the input ex_sc
#' @param genes a vector of genes to go into the heatmap
#' @param type can either be "bulk" or "single_cell. If type is "bulk" it will utilize the bulk info stored in fData
#' @param title A heatmap title
#' @param color_pal The color pallete to be used
#' @param scale_by scale across "row" (genes), "col" (groups), or FALSE
#' @param cluster_by either "row", col, or both
#' @param text_angle The desired angle for text on the group labels
#' @param group_names whether groups should be labelled
#' @param gene_names whether genes should be labelled
#' @param facet_by will create breaks in the heatmap by some pData Variable
#' @export
#' @details
#' Utilize information stored in pData to control the plot display.
#' @examples
#' plot_tsne_metadata(ex_sc_example, color_by = "UMI_sum", title = "UMI_sum across clusters", facet_by = "Cluster", ncol = 3)

plot_heatmap <- function(input, genes, type, title = "Heatmap", scale_by = "row", cluster_by = "row", color_pal = viridis::inferno(256),
                         text_angle = 60, group_names = TRUE, gene_names = TRUE, facet_by = FALSE){
  if(type == "bulk"){
    heat_dat <- fData(input)[,grep("bulk", colnames(fData(input)))]
    heat_dat <- heat_dat[genes,]
  }
  if(type == "single_cell"){
    heat_dat <- exprs(input)[genes,]
  }
  if(scale_by == "row"){
    heat_dat_2 <- t(apply(heat_dat,1,scale))
    colnames(heat_dat_2) <- colnames(heat_dat)
    heat_dat <- heat_dat_2
  }
  if(scale_by == "col"){
    heat_dat_2 <- apply(heat_dat,2,scale)
    rownames(heat_dat_2) <- rownames(heat_dat)
    heat_dat <- heat_dat_2
  }
  if(cluster_by == "row"){
    d <- dist(heat_dat, method = "euclidean")
    hc1 <- hclust(d, method = "complete" )
    heat_dat <- heat_dat[hc1$order,]
  }
  if(cluster_by == "col"){
    d <- dist(t(heat_dat), method = "euclidean")
    hc1 <- hclust(d, method = "complete" )
    heat_dat <- heat_dat[,hc1$order]
  }
  if(cluster_by == "both"){
    d <- dist(heat_dat, method = "euclidean")
    hc1 <- hclust(d, method = "complete" )
    heat_dat <- heat_dat[hc1$order,]
    d <- dist(t(heat_dat), method = "euclidean")
    hc1 <- hclust(d, method = "complete" )
    heat_dat <- heat_dat[,hc1$order]
  }
  #####
  if(type == "single_cell"){
    heat_dat[which(heat_dat > 5)] <- 5
  }
  if(facet_by != FALSE){
   colnames(heat_dat) <- pData(input)[,facet_by]
  }
  heat_dat <- as.data.frame(heat_dat)
  colnames(heat_dat) <- sub("_num_.*", "", colnames(heat_dat) )
  heat_dat_lng <- tidyr::gather(heat_dat, key = "group", "Expression", 1:ncol(heat_dat), factor_key = "TRUE")
  heat_dat_lng$genes <- rep(factor(rownames(heat_dat), levels = rownames(heat_dat)), ncol(heat_dat))
  g <- ggplot(heat_dat_lng, aes(group, genes))
  g <- g + geom_tile(aes(fill = Expression))
  g <- g + theme_classic()
  g <- g + scale_fill_gradientn(colours = color_pal)
  g <- g + theme(plot.title = element_text(size = 20), axis.title = element_text(size = 10), legend.title = element_text(size = 15), legend.text=element_text(size=10))
  g <- g + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
  g <- g +  labs(title= title)
  g <- g + theme(axis.text.x = element_text(angle = text_angle, hjust = 1))
  g <- g +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          line = element_blank())
  if(group_names == FALSE){
    g <- g + theme(axis.text.x=element_blank())
  }
  if(gene_names == FALSE){
    g <- g + theme(axis.text.y=element_blank())
  }
  plot(g)
}

