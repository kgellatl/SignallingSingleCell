#' tSNE Plot on a gene or gene
#'
#' This will plot gene information onto a 2d tsne plot
#'
#' @param input The input data
#' @param gene the gene or gene to ploy
#' @param log_scale if true will log2(value)
#' @param title The title
#' @param colors_points Colors for the cells
#' @param density If true will color a density of the points
#' @param resolution Control the density drawing resolution. Higher values will increase file size
#' @param cutoff removes points from density below this threshold. Lower values will increase file size
#' @param facet_by If ONE gene is being plotted it can be faceted. You cannot facet multiple gene
#' @param ncol controls the number of columns if faceting
#' @param size The size of the points
#' @param xcol pData column to use for x axis
#' @param ycol pData column to use for y axis
#' @export
#' @details
#' Utilize information stored in pData to control the plot display.
#' @examples
#' plot_tsne_gene(ex_sc_example, gene = "Tnf", title = "Tnf over Time", facet_by = "Timepoint", density = TRUE)
#'
plot_tsne_gene <- function(input,
                           gene,
                           log_scale = F,
                           title = "",
                           density = FALSE,
                           facet_by = "NA",
                           colors_points = c("gray", 'blue', 'red', 'yellow'),
                           size = 1.5,
                           ncol = 2,
                           resolution = 500,
                           theme = "classic",
                           cutoff = 0.2,
                           xcol="x",
                           ycol="y",
                           text_sizes = c(20,10,5,10,5,5)){
  kde2d_weighted <- function (x, y, w, h, n , lims = c(range(x), range(y))) {
    nx <- length(x)
    if (length(y) != nx)
      stop("data vectors must be the same length")
    if (length(w) != nx & length(w) != 1)
      stop("weight vectors must be 1 or length of data")
    gx <- seq(lims[1], lims[2], length = n) # gridpoints x
    gy <- seq(lims[3], lims[4], length = n) # gridpoints y
    if (missing(h))
      h <- c(MASS::bandwidth.nrd(x), MASS::bandwidth.nrd(y));
    if (missing(w))
      w <- numeric(nx)+1;
    h <- h/4
    ax <- outer(gx, x, "-")/h[1] # distance of each point to each grid point in x-direction
    ay <- outer(gy, y, "-")/h[2] # distance of each point to each grid point in y-direction
    z <- (matrix(rep(w,n), nrow=n, ncol=nx, byrow=TRUE)*matrix(dnorm(ax), n, nx)) %*% t(matrix(dnorm(ay), n, nx))/(sum(w) * h[1] * h[2]) # z is the density
    return(list(x = gx, y = gy, z = z))
  }
  if(facet_by == "NA"){ #This will allow plotting of 1 - n gene
    geneColored1 <- pData(input)[,c(xcol, ycol)]
    geneColored1 <- do.call(rbind, replicate(length(gene), geneColored1, simplify=FALSE))
    geneColored3 <- t(exprs(input)[gene,])
    for(i in 1:ncol(geneColored3)){
      vals <- geneColored3[,i]
      if(log_scale == TRUE){
        geneColored3[,i] <- log2(vals+2)-1
      } else {
        geneColored3[,i] <- vals
      }
    }
    vals <- as.vector(geneColored3)
    geneColored1$vals <- vals
    gene_final <- c()
    for(i in 1:length(gene)){
      genevec <- rep(gene[i], ncol(exprs(input)))
      gene_final <- c(gene_final, genevec)
    }
    geneColored1$gene <- gene_final
    final_dfdens_norm <- data.frame()
    if(density == TRUE){
      for(i in 1:length(gene)){
        index <- grep(gene[i], geneColored1[,"gene"])
        subset <- geneColored1[index,]
        dens <- kde2d_weighted(subset[,xcol], subset[,ycol], subset$vals, n = resolution)
        dfdens <- data.frame(expand.grid(x=dens$x, y=dens$y), z=as.vector(dens$z))
        dfdens_norm <- dfdens
        dfdens_norm[,3] <- dfdens_norm[,3]/max(dfdens_norm[,3])
        dfdens_norm$gene <- gene[i]
        final_dfdens_norm <- rbind(dfdens_norm, final_dfdens_norm)
      }
      remove <- which(final_dfdens_norm[,"z"] < cutoff)
      final_dfdens_norm <- final_dfdens_norm[-remove,]
      geneColored1 <- geneColored1[with(geneColored1, order(geneColored1[,"gene"])), ]
      final_dfdens_norm <- final_dfdens_norm[with(final_dfdens_norm, order(final_dfdens_norm[,"gene"])), ]
      geneColored1 <- geneColored1[with(geneColored1, order(geneColored1[,3])), ]
    } else {
      geneColored1 <- geneColored1[with(geneColored1, order(geneColored1[,3])), ]
    }
    g <- ggplot(geneColored1)
    if(theme == "bw") {
      g <- g + theme_bw();
    } else {
      g <- g + theme_classic()
    }
    g <- g + theme(plot.title = element_text(size = text_sizes[1]), axis.title = element_text(size = text_sizes[2]), axis.text = element_text(size = text_sizes[3]), legend.title = element_text(size = text_sizes[4]), legend.text=element_text(size=text_sizes[5]))
    g <- g +  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
    g <- g +  facet_wrap(~gene, ncol = ncol)
    if(density == TRUE){
      g <- g +  geom_point(data=final_dfdens_norm, aes_string(x=xcol, y=ycol, color="z"), shape=20, size=.1, stroke = 1)
    }
    g <- g +  scale_color_gradientn(colours=colors_points)
    g <- g +  geom_point(data= geneColored1, aes_string(x=xcol, y=ycol, col="vals"), shape=20, size = size)
    if(title == ""){
      title <- gene
      g <- g +  labs(title= title, col= "UMIs", x = "tSNE[1]", y = "tSNE[2]")
    } else {
      g <- g +  labs(title= title, col= "UMIs", x = "tSNE[1]", y = "tSNE[2]")
    }
  } else { #This will allow plotting of 1 gene, faceted by variable in pData
    if(length(gene) > 1){
      stop("You cannot facet multiple gene")
    }
    geneColored1 <- pData(input)[,c(xcol, ycol, facet_by)]
    geneColored1 <- as.data.frame(geneColored1)
    values <- log2(t(exprs(input)[gene,])+2)-1
    values <- as.vector(values)
    geneColored1$vals <- values
    final_dfdens_norm <- data.frame()
    facets <- unique(pData(input)[,facet_by])
    facets <- sort(facets)
    if(density == TRUE){
      for(i in 1:length(facets)){
        index <- grep(facets[i], geneColored1[,facet_by])
        subset <- geneColored1[index,]
        dens <- kde2d_weighted(subset[,xcol], subset[,ycol], subset$vals, n = resolution)
        dfdens <- data.frame(expand.grid(x=dens$x, y=dens$y), z=as.vector(dens$z))
        dfdens_norm <- dfdens
        dfdens_norm[,3] <- dfdens_norm[,3]/max(dfdens_norm[,3])
        dfdens_norm$facets <- facets[i]
        final_dfdens_norm <- rbind(dfdens_norm, final_dfdens_norm)
      }
      remove <- which(final_dfdens_norm[,"z"] < cutoff)
      final_dfdens_norm <- final_dfdens_norm[-remove,]
      final_dfdens_norm <- final_dfdens_norm[with(final_dfdens_norm, order(final_dfdens_norm[,3])), ]
      geneColored1 <- geneColored1[with(geneColored1, order(geneColored1[,4])), ]
      colnames(final_dfdens_norm) <- c(xcol, ycol, "z", facet_by)
    } else {
      geneColored1 <- geneColored1[with(geneColored1, order(geneColored1[,4])), ]
    }
    tmp <- pData(input)[c(xcol, ycol)]
    g <- ggplot(geneColored1)
    if(density == TRUE){
      g <- g +  geom_point(data=final_dfdens_norm, aes(x=x, y=y, color=z), shape=20, size=.1, stroke = 1)
    }
    g <- g +  geom_point(data= tmp, aes_string(x=xcol, y=ycol), shape=20, size = size, col = "gray")
    g <- g +  theme_classic()
    g <- g + theme(plot.title = element_text(size = text_sizes[1]), axis.title = element_text(size = text_sizes[2]), axis.text = element_text(size = text_sizes[3]), legend.title = element_text(size = text_sizes[4]), legend.text=element_text(size=text_sizes[5]))
    g <- g +  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
    g <- g +  facet_wrap(facets = reformulate(facet_by), ncol = ncol)
    g <- g +  scale_color_gradientn(colours=colors_points)
    g <- g +  geom_point(data= geneColored1, aes_string(x=xcol, y=ycol, col="vals"), shape=20, size = size)
    if(title == ""){
      title <- gene
      g <- g +  labs(title= title, col= "UMIs", x = "tSNE[1]", y = "tSNE[2]")
    } else {
      g <- g +  labs(title= title, col= "UMIs", x = "tSNE[1]", y = "tSNE[2]")
    }
  }
  return(g)
}
