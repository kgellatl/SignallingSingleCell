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
#' @param cluster_type "kmeans" or "hierarchical"
#' @param k if cluster type is kmeans must provide k
#' @param text_angle The desired angle for text on the group labels
#' @param text_sizes a vector of title_size, axis_title, axis_text, legend_title, legend_text, facet_text, faults too c(20,10,5,10,5,5)
#' @param group_names whether groups should be labelled
#' @param gene_names whether genes should be labelled
#' @param facet_by will create breaks in the heatmap by some pData Variable
#' @param pdf_format can be "tile" or "raster." tile is generally higher quality while raster is more efficient
#' @param ceiling A value above which to truncate
#' @param color_facets if true will use colors instead of text labels for the facets
#' @param plotly if true will be interactive
#' # note that this option cannot be saved with save_ggplot(), and also is time consuming for single cell heatmaps
#' @export
#' @details
#' Utilize information stored in pData to control the plot display.
#' @examples
#' plot_tsne_metadata(ex_sc_example, color_by = "UMI_sum", title = "UMI_sum across clusters", facet_by = "Cluster", ncol = 3)

plot_heatmap <- function(input,
                         genes,
                         type,
                         facet_by = FALSE,
                         scale_group = F,
                         title = "Heatmap",
                         scale_by = "row",
                         cluster_by = "row",
                         cluster_type = "hierarchical",
                         k = NULL,
                         show_k = F,
                         ceiling = FALSE,
                         color_pal = viridis::magma(256),
                         color_facets = FALSE,
                         group_names = TRUE,
                         gene_names = TRUE,
                         text_angle = 90,
                         pdf_format = "tile",
                         interactive = FALSE,
                         text_sizes = c(20,10,5,10,5,5),
                         gene_labels = NULL,
                         gene_labels_size = 2,
                         gene_labels_nudge  = -0.5,
                         gene_labels_col = 1,
                         gene_labels_force = 1,
                         return_results = F){
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  #####
  if(facet_by != FALSE){
    n = length(unique(pData(input)[,facet_by]))
    dynamic_colors = gg_color_hue(n)
  }
  #####
  if(type == "bulk"){
    heat_dat <- fData(input)[,grep("bulk", colnames(fData(input)))]
    heat_dat <- heat_dat[genes,]
  }
  #####
  if(type == "single_cell"){
    heat_dat <- exprs(input)[genes,]
    if(facet_by != FALSE){
      metadata <- pData(input)[,facet_by]
    }
  }
  #####
  if(scale_by == "row"){
    if(scale_group != F){
      groups <- unique(pData(input)[,scale_group])

      for (i in 1:length(groups)) {
        int_group <- groups[i]
        ind <- grep(int_group, colnames(heat_dat))
        if(length(ind) == 0){
          stop("scale_group was not used to calculate aggregate bulk")
        }
        mat <- heat_dat[,ind]
        mat <- t(apply(mat,1,scale))
        mat[is.na(mat)] <- 0
        heat_dat[,ind] <- mat

      }


    } else {
      heat_dat_2 <- t(apply(heat_dat,1,scale))
      colnames(heat_dat_2) <- colnames(heat_dat)
      heat_dat <- heat_dat_2
    }
  }

  #####
  if(scale_by == "col"){
    heat_dat_2 <- apply(heat_dat,2,scale)
    rownames(heat_dat_2) <- rownames(heat_dat)
    heat_dat <- heat_dat_2
  }
  #####
  if(cluster_by == "row"){
    if(cluster_type == "hierarchical"){
      d <- dist(heat_dat, method = "euclidean")
      hc1 <- hclust(d, method = "complete" )
      heat_dat <- heat_dat[hc1$order,]
    }
    if(cluster_type == "kmeans"){
      if(is.null(k)){
        stop("Must provide a k for kmeans clustering")
      }
      res <- kmeans(heat_dat, centers = k, nstart = 25)
      reord <- order(res$cluster)
      heat_dat <- heat_dat[reord,]
      set.seed(100)
      for (i in 1:length(1:k)) {
        int <- c(1:k)[i]
        s1 <- names(res$cluster[which(res$cluster == int)])
        ind <- match(s1, rownames(heat_dat))
        tmp <- heat_dat[ind,]
        d <- dist(tmp, method = "euclidean")
        hc1 <- hclust(d, method = "complete" )
        tmp <- tmp[hc1$order,]
        if(i == 1){
          new_dat <- tmp
        } else {
          new_dat <- rbind(new_dat, tmp)
        }
      }
      heat_dat <- new_dat
      res$cluster <- sort(res$cluster)
      res$cluster <- res$cluster[match(rownames(heat_dat), names(res$cluster))]
    }
  }
  #####
  if(cluster_by == "col"){
    if(cluster_type == "hierarchical"){
      d <- dist(t(heat_dat), method = "euclidean")
      hc1 <- hclust(d, method = "complete" )
      heat_dat <- heat_dat[,hc1$order]
      if(type == "single_cell"){
        if(facet_by != FALSE){
          metadata <- metadata[hc1$order]
        }
      }
    }
    if(cluster_type == "kmeans"){
      if(is.null(k)){
        stop("Dont do kmeans on cols")
      }
    }
  }
  #####
  if(cluster_by == "both"){
    if(cluster_type == "hierarchical"){
      d <- dist(heat_dat, method = "euclidean")
      hc1 <- hclust(d, method = "complete" )
      heat_dat <- heat_dat[hc1$order,]
      d <- dist(t(heat_dat), method = "euclidean")
      hc1 <- hclust(d, method = "complete" )
      heat_dat <- heat_dat[,hc1$order]
      if(type == "single_cell"){
        if(facet_by != FALSE){
          metadata <- metadata[hc1$order]
        }
      }
    }
    if(cluster_type == "kmeans"){
      stop("Dont do kmeans on both")
    }
  }
  #####
  if(ceiling != FALSE){
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
        heat_dat_lng$facet[grep(paste0(int, "$"), as.character(heat_dat_lng$group))] <- int
        # vals <- strsplit(as.character(heat_dat_lng$group), split = "_")
        # vals <- matrix(unlist(vals), ncol = length(vals[[1]]), byrow = T)
        # for (j in 1:nrow(vals)) {
        #   int2 <- vals[j,]
        #   ind <- match(int, int2)
        #   if(!is.na(ind)){
        #     heat_dat_lng$facet[j] <- int
        #   }
        # }
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
  g <- g + theme(plot.title = element_text(size = text_sizes[1]), axis.title = element_text(size = text_sizes[2]), axis.text = element_text(size = text_sizes[3]), legend.title = element_text(size = text_sizes[4]), legend.text=element_text(size=text_sizes[5]))
  g <- g + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
  g <- g + labs(title= title)
  g <- g + theme(axis.text.x = element_text(angle = text_angle, hjust = 1))
  g <- g + theme(panel.spacing = unit(0.02, "lines"))
  if(cluster_type == "kmeans"){
    if(show_k == T){
      ylabs <- res$cluster
      ulabs <- unique(ylabs)
      ind <- match(ulabs, ylabs)
      ylabs[-ind] <- ""
      g <- g + scale_y_discrete(labels=ylabs)
    }
  }
  g <- g +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          line = element_blank())
  if(group_names == FALSE){
    g <- g + theme(axis.text.x=element_blank())
  }
  if(!is.null(gene_labels)){
    # ligs_reorder_label <- genes
    # ind <- match(gene_labels, genes)
    # ligs_reorder_label[-ind] <- ""
    g <- g + theme(axis.text.y=element_blank())
    heat_dat_lng$label[heat_dat_lng$genes %in% gene_labels] <- as.character(heat_dat_lng$genes[heat_dat_lng$genes %in% gene_labels])
    rows_lab <- seq(from = (gene_labels_col-1)*nrow(heat_dat)+1, to = (gene_labels_col-1)*nrow(heat_dat)+1+nrow(heat_dat)-1)
    g <- g + ggrepel::geom_text_repel(data = heat_dat_lng[rows_lab,], mapping = ggplot2::aes(label = label), na.rm = T, nudge_x = gene_labels_nudge, direction = "y", min.segment.length = 0, size = gene_labels_size, force = gene_labels_force, max.iter = 10000)

    }
  if(gene_names == FALSE){
    g <- g + theme(axis.text.y=element_blank())
  }
  if(facet_by == FALSE){
    plot(g)
  }
  if(facet_by != FALSE){
    g <- g + facet_grid(reformulate(facet_by), scales = "free_x", space = "free_x")
    if(color_facets != TRUE){
      plot(g)
    }
    if(color_facets == TRUE) {
      dummy <- ggplot(data = heat_dat_lng, aes(fill = CellType, colour=CellType), size = 0.5)+ facet_grid(reformulate(facet_by)) +
        geom_rect(aes(colour=CellType), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
        theme_minimal() + scale_colour_manual(values = dynamic_colors)

      g1 <- ggplotGrob(g)
      g2 <- ggplotGrob(dummy)

      gtable_select <- function (x, ...)
      {
        matches <- c(...)
        x$layout <- x$layout[matches, , drop = FALSE]
        x$grobs <- x$grobs[matches]
        x
      }
      panels <- grepl(pattern="panel", g2$layout$name)
      strips <- grepl(pattern="strip_t", g2$layout$name)
      g2$layout$t[panels] <- g2$layout$t[panels] - 1
      g2$layout$b[panels] <- g2$layout$b[panels] - 1

      new_strips <- gtable_select(g2, panels | strips)
      grid::grid.newpage()
      grid::grid.draw(new_strips)

      gtable_stack <- function(g1, g2){
        g1$grobs <- c(g1$grobs, g2$grobs)
        g1$layout <- transform(g1$layout, z= z-max(z), name="g2")
        g1$layout <- rbind(g1$layout, g2$layout)
        g1
      }

      ## ideally you'd remove the old strips, for now they're just covered
      new_plot <- gtable_stack(g1, new_strips)
      grid::grid.newpage()
      grid::grid.draw(new_plot)
      #####
    }
  }
  if(return_results){
    if(cluster_type == "kmeans"){
      kmean_res <- res
      return_result <- vector(mode = "list", length = 2)
      return_result[[1]] <-  g
      return_result[[2]] <-  kmean_res
      return(return_result)
    } else {
      hc_res <- hc1
      return_result <- vector(mode = "list", length = 2)
      return_result[[1]] <-  g
      return_result[[2]] <-  hc_res
      return(return_result)
    }
  }
  if(interactive == TRUE){
    ggplotly(g, source = "master")
  }
}

