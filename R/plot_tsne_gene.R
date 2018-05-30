#' tSNE Plot on a gene or genes
#'
#' This will plot gene information onto a 2d tsne plot
#'
#' @param input The input data
#' @param genes the gene or genes to ploy
#' @param title The title
#' @param colors_points Colors for the cells
#' @param density If true will color a density of the points
#' @param resolution Control the density drawing resolution. Higher values will increase file size
#' @param cutoff removes points from density below this threshold. Lower values will increase file size
#' @param facet_by If ONE gene is being plotted it can be faceted. You cannot facet multiple genes
#' @param ncol controls the number of columns if faceting
#' @param size The size of the points
#' @export
#' @details
#' Utilize information stored in pData to control the plot display.
#' @examples
#' plot_tsne_gene(ex_sc_example, gene = "Tnf", title = "Tnf over Time", facet_by = "Timepoint", density = TRUE)
#'
plot_tsne_gene <- function(input, genes, title = "", density = FALSE,  facet_by = "NA",
                           colors_points = c("gray", 'blue', 'red', 'yellow'),
                           size = 1.5, ncol = 2, resolution = 500, cutoff = 0.2){
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
  if(facet_by == "NA"){ #This will allow plotting of 1 - n genes
    geneColored1 <- pData(input)[,c("x", "y")]
    geneColored1 <- do.call(rbind, replicate(length(genes), geneColored1, simplify=FALSE))
    geneColored3 <- t(exprs(input)[genes,])
    for(i in 1:ncol(geneColored3)){
      vals <- geneColored3[,i]
      geneColored3[,i] <- log2(vals+2)-1
    }
    vals <- as.vector(geneColored3)
    geneColored1$vals <- vals
    gene_final <- c()
    for(i in 1:length(genes)){
      genesvec <- rep(genes[i], ncol(exprs(input)))
      gene_final <- c(gene_final, genesvec)
    }
    geneColored1$gene <- gene_final
    final_dfdens_norm <- data.frame()
    if(density == TRUE){
      for(i in 1:length(genes)){
        index <- grep(genes[i], geneColored1[,"gene"])
        subset <- geneColored1[index,]
        dens <- kde2d_weighted(subset$x, subset$y, subset$vals, n = resolution)
        dfdens <- data.frame(expand.grid(x=dens$x, y=dens$y), z=as.vector(dens$z))
        dfdens_norm <- dfdens
        dfdens_norm[,3] <- dfdens_norm[,3]/max(dfdens_norm[,3])
        dfdens_norm$gene <- genes[i]
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
    g <- g +  theme_classic()
    g <- g +  theme(plot.title = element_text(size = 15), axis.title = element_text(size = 10), legend.title = element_text(size = 10), legend.text=element_text(size=5))
    g <- g +  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
    g <- g +  facet_wrap(~gene, ncol = ncol)
    if(density == TRUE){
      g <- g +  geom_point(data=final_dfdens_norm, aes(x=x, y=y, color=z), shape=20, size=.1, stroke = 1)
    }
    g <- g +  scale_color_gradientn(colours=colors_points)
    g <- g +  geom_point(data= geneColored1, aes(x=x, y=y, col=vals), shape=20, size = size)
    if(title == ""){
      title <- genes
      g <- g +  labs(title= title, col= "UMIs", x = "tSNE[1]", y = "tSNE[2]")
    } else {
      g <- g +  labs(title= title, col= "UMIs", x = "tSNE[1]", y = "tSNE[2]")
    }
  } else { #This will allow plotting of 1 gene, faceted by variable in pData
    if(length(genes) > 1){
      stop("You cannot facet multiple genes")
    }
    geneColored1 <- pData(input)[,c("x", "y", facet_by)]
    geneColored1 <- as.data.frame(geneColored1)
    values <- log2(t(exprs(input)[genes,])+2)-1
    values <- as.vector(values)
    geneColored1$vals <- values
    final_dfdens_norm <- data.frame()
    facets <- unique(pData(input)[,facet_by])
    facets <- sort(facets)
    if(density == TRUE){
      for(i in 1:length(facets)){
        index <- grep(facets[i], geneColored1[,facet_by])
        subset <- geneColored1[index,]
        dens <- kde2d_weighted(subset$x, subset$y, subset$vals, n = resolution)
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
      colnames(final_dfdens_norm) <- c("x", "y", "z", facet_by)
    } else {
      geneColored1 <- geneColored1[with(geneColored1, order(geneColored1[,4])), ]
    }
    tmp <- pData(input)[c("x", "y")]
    g <- ggplot(geneColored1)
    if(density == TRUE){
      g <- g +  geom_point(data=final_dfdens_norm, aes(x=x, y=y, color=z), shape=20, size=.1, stroke = 1)
    }
    g <- g +  geom_point(data= tmp, aes(x=x, y=y), shape=20, size = size, col = "gray")
    g <- g +  theme_classic()
    g <- g +  theme(plot.title = element_text(size = 15), axis.title = element_text(size = 10), legend.title = element_text(size = 10), legend.text=element_text(size=5))
    g <- g +  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
    g <- g +  facet_wrap(facets = reformulate(facet_by), ncol = ncol)
    g <- g +  scale_color_gradientn(colours=colors_points)
    g <- g +  geom_point(data= geneColored1, aes(x=x, y=y, col=vals), shape=20, size = size)
    if(title == ""){
      title <- genes
      g <- g +  labs(title= title, col= "UMIs", x = "tSNE[1]", y = "tSNE[2]")
    } else {
      g <- g +  labs(title= title, col= "UMIs", x = "tSNE[1]", y = "tSNE[2]")
    }
  }
  return(g)
}
