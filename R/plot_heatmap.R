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
#' @param pdf_format can be "tile" or "raster." tile is generally higher quality while raster is more efficient
#' @param ceiling A value above which to truncate
#' @export
#' @details
#' Utilize information stored in pData to control the plot display.
#' @examples
#' plot_tsne_metadata(ex_sc_example, color_by = "UMI_sum", title = "UMI_sum across clusters", facet_by = "Cluster", ncol = 3)

plot_heatmap <- function(input, genes, type, title = "Heatmap", scale_by = "row", cluster_by = "row", color_pal = viridis::magma(256),
                         text_angle = 60, group_names = TRUE, gene_names = TRUE, facet_by = FALSE, pdf_format = "raster", ceiling = FALSE){
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  if(facet_by != FALSE){
    n = length(unique(pData(input)[,facet_by]))
    dynamic_colors = gg_color_hue(n)
  }
  if(type == "bulk"){
    heat_dat <- fData(input)[,grep("bulk", colnames(fData(input)))]
    heat_dat <- heat_dat[genes,]
  }
  if(type == "single_cell"){
    heat_dat <- exprs(input)[genes,]
    metadata <- pData(input)[,facet_by]
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
    if(type == "single_cell"){
      metadata <- metadata[hc1$order]
    }
  }
  if(cluster_by == "both"){
    d <- dist(heat_dat, method = "euclidean")
    hc1 <- hclust(d, method = "complete" )
    heat_dat <- heat_dat[hc1$order,]
    d <- dist(t(heat_dat), method = "euclidean")
    hc1 <- hclust(d, method = "complete" )
    heat_dat <- heat_dat[,hc1$order]
    if(type == "single_cell"){
      metadata <- metadata[hc1$order]
    }
  }
  #####
  if(ceiling  != FALSE){
    heat_dat[which(heat_dat > ceiling)] <- ceiling
  }

  heat_dat <- as.data.frame(heat_dat)
  if(type == "bulk"){
    colnames(heat_dat) <- sub("_num_.*", "", colnames(heat_dat))
  }
  heat_dat_lng <- tidyr::gather(heat_dat, key = "group", "Expression", 1:ncol(heat_dat), factor_key = "TRUE")
  heat_dat_lng$genes <- rep(factor(rownames(heat_dat), levels = rownames(heat_dat)), ncol(heat_dat))
  if(facet_by != FALSE){
    facs <- unique(pData(input)[,facet_by])
    facs <- sort(facs)
    heat_dat_lng$facet  <- NA
    if(type == "bulk"){
      for (i in 1:length(facs)) {
        int <- facs[i]
        vals <- strsplit(as.character(heat_dat_lng$group), split = "_")
        vals <- matrix(unlist(vals), ncol = length(vals[[1]]), byrow = T)
        for (j in 1:nrow(vals)) {
          int2 <- vals[j,]
          ind <- match(int, int2)
          if(!is.na(ind)){
            heat_dat_lng$facet[j] <- int
          }
        }
      }
      heat_dat_lng$facet <- factor(heat_dat_lng$facet)
      colnames(heat_dat_lng)[ncol(heat_dat_lng)] <- facet_by
    } else {
      heat_dat_lng$facet <- rep(metadata, each = length(genes))
      colnames(heat_dat_lng)[ncol(heat_dat_lng)] <- facet_by
    }
  }
  g <- ggplot(heat_dat_lng, aes(group, genes))
  if(pdf_format == "raster"){
    g <- g + geom_raster(aes(fill = Expression))
  }
  if(pdf_format == "tile"){
    g <- g + geom_tile(aes(fill = Expression, colour=Expression))
  }
  g <- g + theme_classic()
  g <- g + scale_fill_gradientn(colours = color_pal)
  g <- g + scale_color_gradientn(colours = color_pal)
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
  if(facet_by != FALSE){
    g <- g + facet_grid(reformulate(facet_by), scales = "free_x", space = "free_x")
    # g <- g + facet_wrap(reformulate(facet_by), scales = "free_x", ncol = ncol)
   # https://github.com/tidyverse/ggplot2/issues/2096
    #https://stackoverflow.com/questions/19440069/ggplot2-facet-wrap-strip-color-based-on-variable-in-data-set
  }
  plot(g)
}

