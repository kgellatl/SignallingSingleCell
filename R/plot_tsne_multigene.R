#' tSNE Plot with multiple genes
#'
#' This will plot genes onto a 2d tsne plot
#'
#' @param input The input data
#' @param title The title
#' @param genes if provided will color_by the gene
#' @param ncol How many columns if faceting
#' @param size The size of the points
#' @param colors What colors to utilize for categorial data. Be sure it is of the proper length!

#' @export
#' @details
#' Utilize information stored in pData to control the plot display.
#' @examples
#' plot_tsne_multigene(input, title = "UMI_sum across clusters", genes = c("Actb", "Gapdh"), ncol = 4, size = 2)

plot_tsne_multigene <- function(input, title, genes, ncol = 4, size = 1, colors = "NA"){
  geneColored1 <- pData(input)[,c("x", "y")] # Change 4th in facet!
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
  geneColored1 <- geneColored1[with(geneColored1, order(geneColored1[,3])), ]
  ggplot(geneColored1) +
    geom_point(aes(x=x, y=y, col=vals), size = 0.5) + # notice col = dens, column 3 which is expression value
    facet_wrap(~gene, ncol = ncol) +
    scale_color_gradientn(colours=c("gray", 'blue', 'red', 'yellow')) +  # COLORS FOR CONTOUR POINTS, FEEL FREE TO PLAY
    labs(title= title, col= "log2(Normalized\nUMIs)", x = "tSNE[1]", y = "tSNE[2]") +
    theme_classic() +
    theme(plot.title = element_text(size = 25), axis.title = element_text(size = 17.5), legend.title = element_text(size = 17.5), legend.text=element_text(size=10)) +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
}
