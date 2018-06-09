#' Plot of genes by pData variable
#'
#' This will plot information onto a 2d tsne plot
#'
#' @param input the input ex_sc.
#' @param genes What to color points by, either "UMI_sum", or pData categorial variable, ignored if gene is provided
#' @param break_by a pData variable to break the x axis by
#' @param title The title
#' @param scale_by Whether each gene should be scaled by break_by ("rows") or gene ("cols")
#' @export
#' @details
#' Utilize information stored in pData to control the plot display.
#' @examples
#' plot_tsne_metadata(ex_sc_example, color_by = "UMI_sum", title = "UMI_sum across clusters", facet_by = "Cluster", ncol = 3)

plot_gene_dots <- function(input, genes, break_by, title = "", scale_by = FALSE){
  vars <- unique(pData(input)[,break_by])
  vars <- sort(vars)
  plot_dat <- matrix(nrow = length(vars), ncol = length(genes))
  colnames(plot_dat) <- genes
  rownames(plot_dat) <- vars
  for (i in 1:nrow(plot_dat)) {
    int <- vars[i]
    cells <- grep(int, pData(input)[,break_by])
    vals <- apply(exprs(input)[genes,cells],1,mean)
    plot_dat[i,] <- vals
  }
  if(scale_by == "rows"){
    plot_dat <- t(apply(plot_dat,1,scale))
    rownames(plot_dat) <- vars
    colnames(plot_dat) <- genes
  }
  if(scale_by == "cols"){
    plot_dat <- apply(plot_dat,2,scale)
    rownames(plot_dat) <- vars
  }
  plot_dat <- as.data.frame(plot_dat)
  plot_dat_lng <- tidyr::gather(plot_dat, key = "gene", "Expression", 1:ncol(plot_dat), factor_key = "TRUE")
  plot_dat_lng$celltype <- rep(vars, length(genes))
  plot_dat_lng$x <- rep(seq(1:length(vars)), length(genes))
  plot_dat_lng$y <- as.numeric(as.factor(plot_dat_lng$gene))
  g <- ggplot(plot_dat_lng)
  g <- g + geom_point(aes(x = x, y = y, col = Expression, size = Expression))
  g <- g + theme_classic()
  g <- g + scale_color_gradientn(colours=c('blue', 'red', 'yellow'))
  g <- g + scale_x_continuous(breaks = seq(1:nrow(plot_dat)), labels=c(vars))
  g <- g + scale_y_continuous(breaks = seq(1:ncol(plot_dat)), labels=c(genes))
  if(title == ""){
    g <- g + labs(title = paste0("Genes by ", break_by), x = break_by, y = "Genes")
  } else {
    g <- g + labs(title = title, x = break_by, y = "Genes")
  }
  g <- g + theme(plot.title = element_text(size = 20), axis.title = element_text(size = 10), legend.title = element_text(size = 15), legend.text=element_text(size=10))
  g <- g + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
  g <- g + guides(size=FALSE)
  plot(g)
}

